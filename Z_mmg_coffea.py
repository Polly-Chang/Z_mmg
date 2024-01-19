import awkward as ak
from coffea import processor
from coffea.nanoevents.methods import candidate
import hist
import uproot
from coffea.nanoevents import NanoEventsFactory, BaseSchema, NanoAODSchema
import matplotlib.pyplot as plt
from pprint import pprint
import numpy as np

import os
import glob
import sys
import correctionlib
from coffea import lookup_tools
from coffea import analysis_tools
from coffea.lookup_tools import extractor
from coffea.analysis_tools import Weights

import time
import ctypes
import math

from coffea.processor import dict_accumulator, column_accumulator, defaultdict_accumulator
import concurrent.futures

from coffea.lumi_tools import LumiMask

def selec_lumis(events: ak.Array, goldenJson: str)->np.ndarray:
    lumimask = LumiMask(goldenJson)
    return lumimask(events.run, events.luminosityBlock)
#! -------------------Notice--------------------
#  change minBiasXsce = 66000 / 69200 / 72400  
#  and dR region                               
#  sys.argv[1] = minBiasXsec
#  sys.argv[2] = mc_xs
#  sys.argv[3] = try_1121 or try_1121_preEE
#  sys.argv[4] = "mc_file"
#  e.g. "/data5/NanoAOD/Run3Summer22EENanoAODv11/DYGto2LG-1Jets_MLL-50_PTG-10to50_TuneCP5_13p6TeV_amcatnloFXFX-pythia8_postEE/*.root:Events"
#  sys.argv[5] = data output root file
#  e.g. data_EFG.root
#  sys.argv[6] = mc output root file
#  e.g. Zg_postEE_10-50.root
#! ---------------------------------------------
minBiasXsec = str(sys.argv[1])

def Photon_SCEta(EB, EE, photon_eta, photon_phi, PV_x, PV_y, PV_z):
    tg_theta_over_2 = np.exp(-photon_eta)
    tg_theta = 2 * tg_theta_over_2 / (1-tg_theta_over_2*tg_theta_over_2)

    if (EB):
        R = 130
        if (PV_x > 0): angle_x0_y0 = np.arctan(PV_y/PV_x)
        if (PV_x < 0): angle_x0_y0 = math.pi + np.arctan(PV_y/PV_x)
        if (PV_y > 0): angle_x0_y0 = math.pi / 2
        if (PV_y < 0): angle_x0_y0 = -math.pi / 2

        alpha = angle_x0_y0 + (math.pi - photon_phi)
        sin_beta = np.sqrt(PV_x*PV_x + PV_y*PV_y) / R * np.sin(alpha)
        beta = abs(np.arcsin(sin_beta))
        gamma = math.pi/2 - alpha - beta
        l = np.sqrt(R*R + (PV_x*PV_x + PV_y*PV_y) - 2*R*np.sqrt(PV_x*PV_x + PV_y*PV_y)*np.cos(gamma))

        z0_zSC = l / tg_theta
        tg_sctheta = R / (PV_z + z0_zSC)

        sctheta = np.arctan(tg_sctheta)
        if (sctheta < 0): sctheta += math.pi
        tg_sctheta_over_2 = np.tan(sctheta/2)
        SCEta = -np.log(tg_sctheta_over_2)
    if (EE):
        intersection_z = 310 if (photon_eta>0) else -310
        base = intersection_z - PV_z
        r = base * tg_theta

        crystalX = PV_x + r * np.cos(photon_phi)
        crystalY = PV_y + r * np.sin(photon_phi)
        tg_sctheta = np.sqrt(crystalX*crystalX + crystalY*crystalY) /intersection_z
    return photon_eta

def add_pileup_weight(events: ak.Array, pileup_profile: str = None) -> ak.Array:
    """
    Parameter pileup_profile not yet used here. Idea: integrate in metaconditions json which file to read in this function.
    Reading privately produced file here.
    """
    path_pileup_profile = os.path.join(
        os.path.dirname(__file__),
        "Pileup_Prompt_Run2022FG_{}.root".format(minBiasXsec),
    )
    pileup_profile = uproot.open(path_pileup_profile)["pileup"]
    pileup_profile = pileup_profile.to_numpy()[0]
    pileup_profile /= pileup_profile.sum()

    pileup_MC = np.histogram(ak.to_numpy(events.Pileup.nPU), bins=1000, range=(0, 1000))[0].astype("float64")
    # avoid division by zero later
    pileup_MC[pileup_MC == 0.] = 1
    pileup_MC /= pileup_MC.sum()

    # print(len(ak.to_numpy(events.Pileup.nPU)))

    pileup_correction = pileup_profile / pileup_MC
    # remove large MC reweighting factors to prevent artifacts
    pileup_correction[pileup_correction > 10] = 1
    # print(len(pileup_correction))

    # if nPU > 99, weight_pileup = pileup_correction[99]
    npu = ak.to_numpy(events.Pileup.nPU)
    # npu[npu > 99] = 99

    # access value(pileup_correction) from index array(npu)
    # https://stackoverflow.com/questions/74787458/how-to-select-data-by-index-in-array
    weight_pileup = pileup_correction[npu]

    # weight_pileup = ak.to_numpy(pileup_correction)
    events["weight_pileup"] = weight_pileup
    # events["weight_pileup"] = pileup_correction

    return events

year = "2022"
lum_b = 0.09768    #fb^-1
lum_c = 5.0707
lum_d = 3.0063
lum_e = 5.8783
lum_f = 18.0070
lum_g = 3.1219
if str(sys.argv[3]) == "try_NanoAODv12_noBR":
    data = glob.glob("/data1/NanoAOD/22Sep2023-v1/Muon_Run2022E/*.root") + glob.glob("/data1/NanoAOD/22Sep2023-v1/Muon_Run2022F/*.root") + glob.glob("/data1/NanoAOD/22Sep2023-v1/Muon_Run2022G/*.root")
    lum = lum_e + lum_f + lum_g
else:
    data = glob.glob("/data1/NanoAOD/22Sep2023-v1/SingleMuon_Run2022B/*.root") + glob.glob("/data1/NanoAOD/22Sep2023-v1/SingleMuon_Run2022C/*.root") + glob.glob("/data1/NanoAOD/22Sep2023-v1/Muon_Run2022C/*.root") + glob.glob("/data1/NanoAOD/22Sep2023-v1/Muon_Run2022D/*.root")
    lum = lum_c + lum_d

#! times BR when using ttbar
mc_xs = float(sys.argv[2])
# mc_xs = float(sys.argv[2])*0.3254*0.3254   #fb
tstart = time.time()

mc_file = sys.argv[4] + "/*.root:Events"
gen_events = uproot.concatenate([mc_file], ["genWeight"], num_workers=46)
sum_weight = ak.sum(gen_events.genWeight / abs(gen_events.genWeight))
mc_weight = mc_xs * lum / sum_weight
# If calculate mc_weight in processor, we don't know how the executor split events in cpu
# Then the mc_weight will change, so we calculate the mc_weight before processor

here = os.path.dirname(os.path.abspath(__file__))
class MyProcessor (processor.ProcessorABC):
    def __init__(self, year, mc_xs, lum):
        self._year = year
        self._mc_xs = mc_xs
        self._lum = lum

        self._fileset = {
            "2022": {
                "data": data,
                "mc": glob.glob(sys.argv[4] + "/*.root"),
            }
        }

    def process (self, events):
        num_event = len(events)
        dataset = events.metadata["dataset"]
        isRealData = not hasattr(events, "genWeight")

        # Data Lumimask----------------------------------------------
        if isRealData:
            lumimask = selec_lumis(events, "Cert_Collisions2022_355100_362760_Golden.json")
            events = events[lumimask]
            # print(len(events))
            if (len(events) == 0):
                return {
                    dataset: {
                        # "entries": len(events)
                    }
                }

        # Weight------------------------------------------------
        weights = Weights(len(events), storeIndividual=True)
        if isRealData:
            #! weights.add("genweight", np.ones(len(events)))
            weights.add("puweight", np.ones(len(events)))
            weights.add("mcweight", np.ones(len(events)))
        else:
            events_pu = add_pileup_weight(events)
            pu_weight = events_pu.weight_pileup
            weights.add("puweight", pu_weight)
            # mc weight-----------------------
            mc_weight_array = np.ones(len(events))
            mc_weight_array.fill(mc_weight)
            weights.add("mcweight", mc_weight_array)

        weight_array = weights.weight()

        # Cut Setting-------------------------------------------
        weight_array = weight_array[(ak.num(events.Photon) > 0) & (ak.num(events.Muon) > 1)]
        events = events[(ak.num(events.Photon) > 0) & (ak.num(events.Muon) > 1)]

        muons = ak.zip(
            {
                "pt": events.Muon.pt,
                "eta": events.Muon.eta,
                "phi": events.Muon.phi,
                "mass": events.Muon.mass,
                "charge": events.Muon.charge,
                "medium_ID": events.Muon.mediumPromptId,
            },
            with_name = "PtEtaPhiMCandidate",
            behavior = candidate.behavior,
        )
        charge = np.zeros(len(events))
        photons = ak.zip(
            {
                "pt": events.Photon.pt,
                "eta": events.Photon.eta,
                "phi": events.Photon.phi,
                "mass": events.Photon.mass,
                "charge": charge,
                "EB": events.Photon.isScEtaEB,
                "EE": events.Photon.isScEtaEE,
                "electronveto": events.Photon.electronVeto,
                "pixelseed": events.Photon.pixelSeed,
                "r9": events.Photon.r9,
                "cutBased": events.Photon.cutBased,
                "mvaID": events.Photon.mvaID,
                "mvaID_WP80": events.Photon.mvaID_WP80,
                "mvaID_WP90": events.Photon.mvaID_WP90,
                "gain": events.Photon.seedGain,
                "PV_x": events.PV.x,
                "PV_y": events.PV.y,
                "PV_z": events.PV.z,
                "PV_npvs": events.PV.npvs,
                "event_num": events.event,
            },
            with_name = "PtEtaPhiMCandidate",
            behavior = candidate.behavior,
        )
        if isRealData:
            photons["genPartFlav"] = 1  # give it an arbitrary number, cuz only MC have gen particle
            photons["genPartIdx"] = 1
            muons["genPartIdx"] = 1

        else:
            photons["genPartFlav"] = events.Photon.genPartFlav
            photons["genPartIdx"] = events.Photon.genPartIdx
            muons["genPartIdx"] = events.Muon.genPartIdx

        h_mass = (
            hist.Hist.new
            .Reg(10, 80,100, name="mass", label="$m_{\mu\mu\gamma}$ [GeV]")
            .Weight()
        )
        
        # Object Selection--------------------------------------
        # Can also from coffea.analysis_tools import PackedSelection and use PackedSelection() to select events
        muons = muons[(muons.medium_ID) & (muons.pt > 4) & (abs(muons.eta) < 2.4)]
        photons = photons[(photons.pt > 10) & ((photons.EB) | (photons.EE))]

        # Event Selection-------------------------
        event_cut = (
            (ak.num(muons) > 1)
            & (ak.num(photons) > 0)
        )
        events = events[event_cut]
        weight_array = weight_array[event_cut]
        muons = muons[event_cut]
        photons = photons[event_cut]

        sort_idx_mu = ak.argsort(muons.pt, axis = 1, ascending = False)
        sort_mu = muons[sort_idx_mu]
        # Traditionally, 0 will be leading, and 1 will be trailing
        # if (ak.all(muons.pt == sort_mu.pt)):
        #    print("the same")
        sort_idx_pho = ak.argsort(photons.pt, axis = 1, ascending = False)
        sort_pho = photons[sort_idx_pho]
        # if (ak.all(photons.pt == sort_pho.pt)):
        #    print("the same")

        cut_HLT_doublemu = (
            (events.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL)
            | (events.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ)
            | (events.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass8)
        )
        cut_HLT_singlemu = (
            events.HLT.IsoMu24
        )
        cut_doublemu = (
            (sort_mu[:, 0].pt > 20)
            & (sort_mu[:, 1].pt > 10)
        )
        cut_singlemu = (
            (sort_mu[:, 0].pt > 25)
        )

        events = events[(cut_HLT_doublemu & cut_doublemu) | (cut_HLT_singlemu & cut_singlemu)]
        weight_array = weight_array[(cut_HLT_doublemu & cut_doublemu) | (cut_HLT_singlemu & cut_singlemu)]
        sort_mu = sort_mu[(cut_HLT_doublemu & cut_doublemu) | (cut_HLT_singlemu & cut_singlemu)]
        sort_pho = sort_pho[(cut_HLT_doublemu & cut_doublemu) | (cut_HLT_singlemu & cut_singlemu)]

        Z = sort_mu[:, 0] + sort_mu[:, 1] + sort_pho  # same size with sort_pho
        mm = sort_mu[:, 0] + sort_mu[:, 1]            # num=1

        cut_kinematic = (
            ((sort_mu[:, 0].delta_r(sort_pho) < 0.8) | (sort_mu[:, 1].delta_r(sort_pho) < 0.8))
            & ((sort_mu[:, 0].delta_r(sort_pho) > 0.4) & (sort_mu[:, 1].delta_r(sort_pho) > 0.4))
            & ((Z.mass > 80) & (Z.mass < 100))
            & ((mm.mass + Z.mass) < 180)
        )   #same size with sort_pho
        sort_pho = sort_pho[cut_kinematic]

        num_pho = (
            ak.num(sort_pho) > 0
        )

        events = events[num_pho]
        weight_array = weight_array[num_pho]
        sort_mu = sort_mu[num_pho]
        sort_pho = sort_pho[num_pho]

        Z = sort_mu[:, 0] + sort_mu[:, 1] + sort_pho #sort_pho size
        choose_pho = ak.argsort(np.fabs(Z.mass - 91.18), axis = 1, ascending = True) #sort_pho size
        sort_pho = sort_pho[choose_pho]

        Z = sort_mu[:, 0] + sort_mu[:, 1] + sort_pho[:, 0]

        ################### S+S calibration parameter storing ###################
        if (str(sys.argv[3]) == "try_NanoAODv12"):
            SS_json = "SS_EFG.json"
            Smearing_json = "Prompt2022FG_SmearingJSON"
            Scale_json = "Prompt2022FG_ScaleJSON"
        else:
            SS_json = "SS_BCD.json"
            Smearing_json = "Rereco2022BCD_SmearingJSON"
            Scale_json = "Rereco2022BCD_ScaleJSON"

        if (ak.any(ak.num(sort_pho) < 1)):
            print("num of sort_pho < 1")

        run = events.run
        photons_gain = sort_pho[:, 0].gain
        photons_pt = sort_pho[:, 0].pt
        sceta = map(Photon_SCEta, sort_pho[:, 0].EB, sort_pho[:, 0].EE, sort_pho[:, 0].eta, sort_pho[:, 0].phi, sort_pho[:, 0].PV_x, sort_pho[:, 0].PV_y, sort_pho[:, 0].PV_z)
        photons_sceta = np.array(list(sceta))
        photons_r9 = sort_pho[:, 0].r9
        # #------------- Load the correctionlib JSON file -------------#
        evaluator = correctionlib.CorrectionSet.from_file(SS_json)
        if isRealData:
            # Scale #
            evaluator_scale = evaluator[Scale_json]
            scale = evaluator_scale.evaluate("total_correction", photons_gain, run, photons_sceta, photons_r9, photons_pt)
            # Scale is multiplicative correction, unlike smearing, it is deterministic
            photons_pt_corrected = scale * photons_pt

            # Validation: The transverse momenta of the photons actually changed after applying the scale correction
            # print("Nominal/Scale pT for each photon:", photons_pt/photons_pt_corrected)
        else:
            # Smearing #
            evaluator_smearing = evaluator[Smearing_json]
            rho = evaluator_smearing.evaluate("rho", photons_sceta, photons_r9)  # rho is internal scale and smearing lingo
            rng = np.random.default_rng(seed=125)  # The smearing is done statistically, so we need some random numbers
            smearing = rng.normal(loc=1., scale=rho)
            photons_pt_corrected = smearing * photons_pt

            # Validation: The transverse momenta of the photons actually changed
            # print("Nominal/Smeared pT for each photon:", photons_pt/photons_pt_corrected)

        if isRealData:
            pho_genpt = np.zeros(len(events))
            pho_geneta = np.zeros(len(events))
            pho_genphi = np.zeros(len(events))
            pho_genmass = np.zeros(len(events))
            pho_genid = np.zeros(len(events))
            pho_genmotherid = np.zeros(len(events))
            mu1_genpt = np.zeros(len(events))
            mu1_geneta = np.zeros(len(events))
            mu1_genphi = np.zeros(len(events))
            mu1_genmass = np.zeros(len(events))
            mu1_genid = np.zeros(len(events))
            mu1_genmotherid = np.zeros(len(events))
            mu2_genpt = np.zeros(len(events))
            mu2_geneta = np.zeros(len(events))
            mu2_genphi = np.zeros(len(events))
            mu2_genmass = np.zeros(len(events))
            mu2_genid = np.zeros(len(events))
            mu2_genmotherid = np.zeros(len(events))
        else:
            genpart = events.GenPart
            pho_genindices = sort_pho.genPartIdx
            mu_genindices = sort_mu.genPartIdx
            # gen_indices = sort_pho[:, 0].genPartIdx
            # [sort_pho[:, 0].genPartIdx != -1]

            # same concept as pu weight, take the index of gen_indices to get the genpart
            pho_genpart = genpart[pho_genindices]
            mu_genpart = genpart[mu_genindices]
            # genpart = genpart[abs(genpart.genPartIdxMother) == 13]
            # genpart = genpart[abs(genpart[genpart.genPartIdxMother].pdgId) == 13]

            pho_genpt = ak.to_numpy(pho_genpart[:, 0].pt)
            pho_geneta = ak.to_numpy(pho_genpart[:, 0].eta)
            pho_genphi = ak.to_numpy(pho_genpart[:, 0].phi)
            pho_genmass = ak.to_numpy(pho_genpart[:, 0].mass)
            pho_genid = ak.to_numpy(pho_genpart[:, 0].pdgId)
            pho_genmotherid = ak.to_numpy(genpart[pho_genpart.genPartIdxMother][:, 0].pdgId)

            mu1_genpt = ak.to_numpy(mu_genpart[:, 0].pt)
            mu1_geneta = ak.to_numpy(mu_genpart[:, 0].eta)
            mu1_genphi = ak.to_numpy(mu_genpart[:, 0].phi)
            mu1_genmass = ak.to_numpy(mu_genpart[:, 0].mass)
            mu1_genid = ak.to_numpy(mu_genpart[:, 0].pdgId)
            mu1_genmotherid = ak.to_numpy(genpart[mu_genpart.genPartIdxMother][:, 0].pdgId)

            mu2_genpt = ak.to_numpy(mu_genpart[:, 1].pt)
            mu2_geneta = ak.to_numpy(mu_genpart[:, 1].eta)
            mu2_genphi = ak.to_numpy(mu_genpart[:, 1].phi)
            mu2_genmass = ak.to_numpy(mu_genpart[:, 1].mass)
            mu2_genid = ak.to_numpy(mu_genpart[:, 1].pdgId)
            mu2_genmotherid = ak.to_numpy(genpart[mu_genpart.genPartIdxMother][:, 1].pdgId)

        Z = sort_mu[:, 0] + sort_mu[:, 1] + sort_pho[:, 0]

        # ---------------------------------------------------
        h_mass.fill(
            mass = Z.mass,
            weight = weight_array
        )

        sort_m1_pt = ak.to_numpy(sort_mu[:, 0].pt)
        sort_m2_pt = ak.to_numpy(sort_mu[:, 1].pt)
        sort_pho_pt_origin = ak.to_numpy(sort_pho[:, 0].pt)
        sort_pho_pt = ak.to_numpy(photons_pt_corrected)
        sort_m1_eta = ak.to_numpy(sort_mu[:, 0].eta)
        sort_m2_eta = ak.to_numpy(sort_mu[:, 1].eta)
        sort_pho_eta = ak.to_numpy(sort_pho[:, 0].eta)
        sort_m1_phi = ak.to_numpy(sort_mu[:, 0].phi)
        sort_m2_phi = ak.to_numpy(sort_mu[:, 1].phi)
        sort_pho_phi = ak.to_numpy(sort_pho[:, 0].phi)
        sort_m1_m = ak.to_numpy(sort_mu[:, 0].mass)
        sort_m2_m = ak.to_numpy(sort_mu[:, 1].mass)
        sort_pho_m = ak.to_numpy(sort_pho[:, 0].mass)
        pho_electronveto = ak.to_numpy(sort_pho[:, 0].electronveto)
        pho_pixelseed = ak.to_numpy(sort_pho[:, 0].pixelseed)
        pho_r9 = ak.to_numpy(sort_pho[:, 0].r9)
        pho_cutBased = ak.to_numpy(sort_pho[:, 0].cutBased)
        # pho_cutBased_Fall17V2 = ak.to_numpy(sort_pho[:, 0].cutBased_Fall17V2)
        pho_mvaID = ak.to_numpy(sort_pho[:, 0].mvaID)
        # pho_mvaID_Fall17V2 = ak.to_numpy(sort_pho[:, 0].mvaID_Fall17V2)
        # pho_mvaID_Fall17V2_WP80 = ak.to_numpy(sort_pho[:, 0].mvaID_Fall17V2_WP80)
        # pho_mvaID_Fall17V2_WP90 = ak.to_numpy(sort_pho[:, 0].mvaID_Fall17V2_WP90)
        pho_mvaID_WP80 = ak.to_numpy(sort_pho[:, 0].mvaID_WP80)
        pho_mvaID_WP90 = ak.to_numpy(sort_pho[:, 0].mvaID_WP90)
        pho_EB = ak.to_numpy(sort_pho[:, 0].EB)
        pho_EE = ak.to_numpy(sort_pho[:, 0].EE)
        weights_ = ak.to_numpy(weight_array)
        Z_m = ak.to_numpy(Z.mass)

        PV_x = ak.to_numpy(sort_pho[:, 0].PV_x)
        PV_y = ak.to_numpy(sort_pho[:, 0].PV_y)
        PV_z = ak.to_numpy(sort_pho[:, 0].PV_z)
        PV_npvs = ak.to_numpy(sort_pho[:, 0].PV_npvs)

        event_num = ak.to_numpy(sort_pho[:, 0].event_num)

        genPartFlav = ak.to_numpy(sort_pho[:, 0].genPartFlav)
        
        adef = dict_accumulator(
            {
            "m1_pt": column_accumulator(np.zeros(shape=(0,))),
            "m2_pt": column_accumulator(np.zeros(shape=(0,))),
            "pho_pt_origin": column_accumulator(np.zeros(shape=(0,))),
            "pho_pt": column_accumulator(np.zeros(shape=(0,))),
            "m1_eta": column_accumulator(np.zeros(shape=(0,))),
            "m2_eta": column_accumulator(np.zeros(shape=(0,))),
            "pho_eta": column_accumulator(np.zeros(shape=(0,))),
            "m1_phi": column_accumulator(np.zeros(shape=(0,))),
            "m2_phi": column_accumulator(np.zeros(shape=(0,))),
            "pho_phi": column_accumulator(np.zeros(shape=(0,))),
            "m1_m": column_accumulator(np.zeros(shape=(0,))),
            "m2_m": column_accumulator(np.zeros(shape=(0,))),
            "pho_m": column_accumulator(np.zeros(shape=(0,))),
            "pho_electronveto": column_accumulator(np.zeros(shape=(0,))),
            "pho_pixelseed": column_accumulator(np.zeros(shape=(0,))),
            "pho_r9": column_accumulator(np.zeros(shape=(0,))),
            "pho_cutBased": column_accumulator(np.zeros(shape=(0,))),
            "pho_mvaID": column_accumulator(np.zeros(shape=(0,))),
            "pho_mvaID_WP80": column_accumulator(np.zeros(shape=(0,))),
            "pho_mvaID_WP90": column_accumulator(np.zeros(shape=(0,))),
            "pho_EB": column_accumulator(np.zeros(shape=(0,))),
            "pho_EE": column_accumulator(np.zeros(shape=(0,))),
            "weights": column_accumulator(np.zeros(shape=(0,))),
            "Z_m": column_accumulator(np.zeros(shape=(0,))),
            "PV_x": column_accumulator(np.zeros(shape=(0,))),
            "PV_y": column_accumulator(np.zeros(shape=(0,))),
            "PV_z": column_accumulator(np.zeros(shape=(0,))),
            "PV_npvs": column_accumulator(np.zeros(shape=(0,))),
            "event_num": column_accumulator(np.zeros(shape=(0,))),
            "genPartFlav": column_accumulator(np.zeros(shape=(0,))),
            "pho_genpt": column_accumulator(np.zeros(shape=(0,))),
            "pho_geneta": column_accumulator(np.zeros(shape=(0,))),
            "pho_genphi": column_accumulator(np.zeros(shape=(0,))),
            "pho_genmass": column_accumulator(np.zeros(shape=(0,))),
            "pho_genid": column_accumulator(np.zeros(shape=(0,))),
            "pho_genmotherid": column_accumulator(np.zeros(shape=(0,))),
            "mu1_genpt": column_accumulator(np.zeros(shape=(0,))),
            "mu1_geneta": column_accumulator(np.zeros(shape=(0,))),
            "mu1_genphi": column_accumulator(np.zeros(shape=(0,))),
            "mu1_genmass": column_accumulator(np.zeros(shape=(0,))),
            "mu1_genid": column_accumulator(np.zeros(shape=(0,))),
            "mu1_genmotherid": column_accumulator(np.zeros(shape=(0,))),
            "mu2_genpt": column_accumulator(np.zeros(shape=(0,))),
            "mu2_geneta": column_accumulator(np.zeros(shape=(0,))),
            "mu2_genphi": column_accumulator(np.zeros(shape=(0,))),
            "mu2_genmass": column_accumulator(np.zeros(shape=(0,))),
            "mu2_genid": column_accumulator(np.zeros(shape=(0,))),
            "mu2_genmotherid": column_accumulator(np.zeros(shape=(0,))),
            
            }
        )
        acc = adef.identity()
        acc["m1_pt"] += column_accumulator(sort_m1_pt)
        acc["m2_pt"] += column_accumulator(sort_m2_pt)
        acc["pho_pt_origin"] += column_accumulator(sort_pho_pt_origin)
        acc["pho_pt"] += column_accumulator(sort_pho_pt)
        acc["m1_eta"] += column_accumulator(sort_m1_eta)
        acc["m2_eta"] += column_accumulator(sort_m2_eta)
        acc["pho_eta"] += column_accumulator(sort_pho_eta)
        acc["m1_phi"] += column_accumulator(sort_m1_phi)
        acc["m2_phi"] += column_accumulator(sort_m2_phi)
        acc["pho_phi"] += column_accumulator(sort_pho_phi)
        acc["m1_m"] += column_accumulator(sort_m1_m)
        acc["m2_m"] += column_accumulator(sort_m2_m)
        acc["pho_m"] += column_accumulator(sort_pho_m)
        acc["pho_electronveto"] += column_accumulator(pho_electronveto)
        acc["pho_pixelseed"] += column_accumulator(pho_pixelseed)
        acc["pho_r9"] += column_accumulator(pho_r9)
        acc["pho_cutBased"] += column_accumulator(pho_cutBased)
        acc["pho_mvaID"] += column_accumulator(pho_mvaID)
        acc["pho_mvaID_WP80"] += column_accumulator(pho_mvaID_WP80)
        acc["pho_mvaID_WP90"] += column_accumulator(pho_mvaID_WP90)
        acc["pho_EB"] += column_accumulator(pho_EB)
        acc["pho_EE"] += column_accumulator(pho_EE)
        acc["weights"] += column_accumulator(weights_)
        acc["Z_m"] += column_accumulator(Z_m)
        acc["PV_x"] += column_accumulator(PV_x)
        acc["PV_y"] += column_accumulator(PV_y)
        acc["PV_z"] += column_accumulator(PV_z)
        acc["PV_npvs"] += column_accumulator(PV_npvs)
        acc["event_num"] += column_accumulator(event_num)
        acc["genPartFlav"] += column_accumulator(genPartFlav)
        acc["pho_genpt"] += column_accumulator(pho_genpt)
        acc["pho_geneta"] += column_accumulator(pho_geneta)
        acc["pho_genphi"] += column_accumulator(pho_genphi)
        acc["pho_genmass"] += column_accumulator(pho_genmass)
        acc["pho_genid"] += column_accumulator(pho_genid)
        acc["pho_genmotherid"] += column_accumulator(pho_genmotherid)
        acc["mu1_genpt"] += column_accumulator(mu1_genpt)
        acc["mu1_geneta"] += column_accumulator(mu1_geneta)
        acc["mu1_genphi"] += column_accumulator(mu1_genphi)
        acc["mu1_genmass"] += column_accumulator(mu1_genmass)
        acc["mu1_genid"] += column_accumulator(mu1_genid)
        acc["mu1_genmotherid"] += column_accumulator(mu1_genmotherid)
        acc["mu2_genpt"] += column_accumulator(mu2_genpt)
        acc["mu2_geneta"] += column_accumulator(mu2_geneta)
        acc["mu2_genphi"] += column_accumulator(mu2_genphi)
        acc["mu2_genmass"] += column_accumulator(mu2_genmass)
        acc["mu2_genid"] += column_accumulator(mu2_genid)
        acc["mu2_genmotherid"] += column_accumulator(mu2_genmotherid)

        return {
            dataset: {
                "total_num": num_event,
                "entries": len(events),
                "hist": h_mass,
                "acc": acc,
            }
        }
    def postprocess(self, accumulator):
        pass

p = MyProcessor(year, mc_xs, lum)

fileset = p._fileset[p._year]

# iterative_run = processor.Runner(
#     executor = processor.IterativeExecutor(compression=None),
#     schema=BaseSchema,
# )
# out = iterative_run(
#     fileset,
#     treename = "Events",
#     processor_instance = MyProcessor(year),
# )

futures_run = processor.Runner(
    executor = processor.FuturesExecutor(compression=None, workers=46),
    schema = NanoAODSchema,
    # maxchunks = 2,
)
out = futures_run(
    fileset,
    treename = "Events",
    processor_instance = MyProcessor(year, mc_xs, lum)
)

# events = NanoEventsFactory.from_root(
#     file,
#     entry_stop=10000,
#     metadata={"dataset": "Z_mmg"},
#     schemaclass=NanoAODSchema,
# ).events()

# out = p.process(events)
elapsed = time.time() - tstart
print(elapsed)
# print(out)


with uproot.recreate("./{}/{}/{}".format(sys.argv[3], minBiasXsec, sys.argv[5])) as f:
    f["tree"] = {"m1_pt": (out["data"]["acc"]["m1_pt"].value), "m2_pt": (out["data"]["acc"]["m2_pt"].value), "pho_pt_origin": (out["data"]["acc"]["pho_pt_origin"].value), "pho_pt": (out["data"]["acc"]["pho_pt"].value),
                "m1_eta": (out["data"]["acc"]["m1_eta"].value), "m2_eta": (out["data"]["acc"]["m2_eta"].value), "pho_eta": (out["data"]["acc"]["pho_eta"].value),
                "m1_phi": (out["data"]["acc"]["m1_phi"].value), "m2_phi": (out["data"]["acc"]["m2_phi"].value), "pho_phi": (out["data"]["acc"]["pho_phi"].value),
                "m1_m": (out["data"]["acc"]["m1_m"].value), "m2_m": (out["data"]["acc"]["m2_m"].value), "pho_m": (out["data"]["acc"]["pho_m"].value),
                "pho_electronveto": (out["data"]["acc"]["pho_electronveto"].value), "pho_pixelseed": (out["data"]["acc"]["pho_pixelseed"].value), "pho_r9": (out["data"]["acc"]["pho_r9"].value),
                "pho_cutBased": (out["data"]["acc"]["pho_cutBased"].value),
                "pho_mvaID": (out["data"]["acc"]["pho_mvaID"].value), "pho_mvaID_WP80": (out["data"]["acc"]["pho_mvaID_WP80"].value), "pho_mvaID_WP90": (out["data"]["acc"]["pho_mvaID_WP90"].value),
                "pho_EB": (out["data"]["acc"]["pho_EB"].value), "pho_EE": (out["data"]["acc"]["pho_EE"].value),
                "weights": (out["data"]["acc"]["weights"].value), "Z_m": (out["data"]["acc"]["Z_m"].value),
                "PV_x": (out["data"]["acc"]["PV_x"].value), "PV_y": (out["data"]["acc"]["PV_y"].value), "PV_z": (out["data"]["acc"]["PV_z"].value), "PV_npvs": (out["data"]["acc"]["PV_npvs"].value),
                "event_num": (out["data"]["acc"]["event_num"].value)}
    f["tree"].show()
with uproot.recreate("./{}/{}/{}".format(sys.argv[3], minBiasXsec, sys.argv[6])) as f:
    f["tree"] = {"m1_pt": (out["mc"]["acc"]["m1_pt"].value), "m2_pt": (out["mc"]["acc"]["m2_pt"].value), "pho_pt_origin": (out["mc"]["acc"]["pho_pt_origin"].value), "pho_pt": (out["mc"]["acc"]["pho_pt"].value),
                "m1_eta": (out["mc"]["acc"]["m1_eta"].value), "m2_eta": (out["mc"]["acc"]["m2_eta"].value), "pho_eta": (out["mc"]["acc"]["pho_eta"].value),
                "m1_phi": (out["mc"]["acc"]["m1_phi"].value), "m2_phi": (out["mc"]["acc"]["m2_phi"].value), "pho_phi": (out["mc"]["acc"]["pho_phi"].value),
                "m1_m": (out["mc"]["acc"]["m1_m"].value), "m2_m": (out["mc"]["acc"]["m2_m"].value), "pho_m": (out["mc"]["acc"]["pho_m"].value),
                "pho_electronveto": (out["mc"]["acc"]["pho_electronveto"].value), "pho_pixelseed": (out["mc"]["acc"]["pho_pixelseed"].value), "pho_r9": (out["mc"]["acc"]["pho_r9"].value),
                "pho_cutBased": (out["mc"]["acc"]["pho_cutBased"].value),
                "pho_mvaID": (out["mc"]["acc"]["pho_mvaID"].value), "pho_mvaID_WP80": (out["mc"]["acc"]["pho_mvaID_WP80"].value), "pho_mvaID_WP90": (out["mc"]["acc"]["pho_mvaID_WP90"].value),
                "pho_EB": (out["mc"]["acc"]["pho_EB"].value), "pho_EE": (out["mc"]["acc"]["pho_EE"].value),
                "weights": (out["mc"]["acc"]["weights"].value), "Z_m": (out["mc"]["acc"]["Z_m"].value),
                "PV_x": (out["mc"]["acc"]["PV_x"].value), "PV_y": (out["mc"]["acc"]["PV_y"].value), "PV_z": (out["mc"]["acc"]["PV_z"].value), "PV_npvs": (out["mc"]["acc"]["PV_npvs"].value),
                "event_num": (out["mc"]["acc"]["event_num"].value),
                "genPartFlav": (out["mc"]["acc"]["genPartFlav"].value),
                "pho_genpt": (out["mc"]["acc"]["pho_genpt"].value), "pho_geneta": (out["mc"]["acc"]["pho_geneta"].value), "pho_genphi": (out["mc"]["acc"]["pho_genphi"].value), "pho_genmass": (out["mc"]["acc"]["pho_genmass"].value),
                "pho_genid": (out["mc"]["acc"]["pho_genid"].value), "pho_genmotherid": (out["mc"]["acc"]["pho_genmotherid"].value),
                "mu1_genpt": (out["mc"]["acc"]["mu1_genpt"].value), "mu1_geneta": (out["mc"]["acc"]["mu1_geneta"].value), "mu1_genphi": (out["mc"]["acc"]["mu1_genphi"].value), "mu1_genmass": (out["mc"]["acc"]["mu1_genmass"].value),
                "mu1_genid": (out["mc"]["acc"]["mu1_genid"].value), "mu1_genmotherid": (out["mc"]["acc"]["mu1_genmotherid"].value),
                "mu2_genpt": (out["mc"]["acc"]["mu2_genpt"].value), "mu2_geneta": (out["mc"]["acc"]["mu2_geneta"].value), "mu2_genphi": (out["mc"]["acc"]["mu2_genphi"].value), "mu2_genmass": (out["mc"]["acc"]["mu2_genmass"].value),
                "mu2_genid": (out["mc"]["acc"]["mu2_genid"].value), "mu2_genmotherid": (out["mc"]["acc"]["mu2_genmotherid"].value)}
    f["tree"].show()