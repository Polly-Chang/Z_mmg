vector<string> IDs = {"Loose","Medium","Tight", "MVA80", "MVA90"};
vector<string> Cuts = {"HasPix","CSEV"};
vector<string> R9 = {"InclusiveR9","HR9","LR9"};
vector<string> Region = {"EB","EE", ""};
vector<string> category = {"EB Inc.","EB high R9","EB low R9","EE Inc.","EE high R9","EE low R9"};
// vector<string> lumi = {"20.665 fb^{-1}"};//F+G
vector<string> measurement ={"phoPt","phoPhi", "phoEta", "zmass", "mva", "Sceta", "npvs"};
string ColorName[3] = {"#533e2d","#e85d04"};
vector<string> ColorName_sp = {"#303030", "#e21e08", "#228b22", "#4040ff", "#f0768b"};
// vector<string> ColorName_sp = {"#303030", "#e21e08", "#228b22", "#4040ff"};


double Error_Propagation(double f, double sigmaA, double A, double sigmaB, double B){
     double error;
     error = f * sqrt((sigmaA/A)*(sigmaA/A)+(sigmaB/B)*(sigmaB/B));
     return error;
}

vector<float> SummaryPlot(vector<TEfficiency *> eff){
    //---------- eff[0]->Prompt, eff[1]->MC ----------//
    eff[0]->SetStatisticOption(TEfficiency::kBUniform);
    eff[0]->SetConfidenceLevel(0.683);
    eff[0]->SetPosteriorMode(1);
    eff[1]->SetStatisticOption(TEfficiency::kBUniform);
    eff[1]->SetConfidenceLevel(0.683);
    eff[1]->SetPosteriorMode(1);
    //---------- SFs calculation using Zmmg ----------//
    TH1F *heff_prompt = (TH1F*) eff[0]->GetCopyPassedHisto();
    heff_prompt->Divide((TH1F*) eff[0]->GetCopyTotalHisto());
    TH1F *heff_MC = (TH1F*) eff[1]->GetCopyPassedHisto();
    heff_MC->Divide((TH1F*) eff[1]->GetCopyTotalHisto());
    TH1F *hSF_prompt = (TH1F*) heff_prompt->Clone();
    hSF_prompt->Divide(heff_MC);
    vector<float> sf_unc = {};

    //---------- Error Propagation ----------//
    int SFNbins_prompt = heff_prompt->GetNbinsX();
    // cout << "SFNbins_prompt: " << SFNbins_prompt << endl;
    float ey_prompt;
    for (int i = 0; i < SFNbins_prompt; i++){
        float data_error = TMath::Sqrt(0.5*(pow(eff[0]->GetEfficiencyErrorUp(i+1),2) + pow(eff[0]->GetEfficiencyErrorLow(i+1),2)));
        float MC_error = TMath::Sqrt(0.5*(pow(eff[1]->GetEfficiencyErrorUp(i+1),2) + pow(eff[1]->GetEfficiencyErrorLow(i+1),2)));
        //The same as TGraphAssmError->GetErrorY() method, but if there are missing points, GetErrorY does not work.
        
        float f = hSF_prompt->GetBinContent(i+1);
        ey_prompt = Error_Propagation(f, data_error , heff_prompt->GetBinContent(i+1), MC_error, heff_MC->GetBinContent(i+1));
        sf_unc.push_back(f);
        sf_unc.push_back(ey_prompt);
    }

    return sf_unc;
}

void DrawSummaryPlot(vector<vector<float>> sp_value, string pu){

    vector<vector<TH1F*>> vhsp;
    for (int iCuts = 0; iCuts < Cuts.size(); iCuts++){

        // TFile *fsf_out = TFile::Open(Form("./try_NanoAODv12_preEE/%s_SummaryPlot_22.root", Cuts[iCuts].c_str()), "RECREATE");
        // fsf_out->cd();

        vhsp.push_back(vector<TH1F*>());

        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gStyle->SetOptStat(0);
        TCanvas *canvas = new TCanvas("c", "c", 800, 800);
        canvas->SetBottomMargin(0.06);
        canvas->SetRightMargin(0.05);
        // canvas->SetTopMargin(0.095);
        canvas->SetLeftMargin(0.13);

        TLegend *legend = new TLegend(0.15, 0.09, 0.55, 0.22); //(x1,y1,x2,y2)

        //------sp_value.size() = 90; 45 for CSEV and 45 for pixelseed------//
        vector<float> y_loose = {sp_value[0+(iCuts*45)][0], sp_value[1+(iCuts*45)][0], sp_value[2+(iCuts*45)][0], sp_value[3+(iCuts*45)][0], sp_value[4+(iCuts*45)][0], sp_value[5+(iCuts*45)][0]};
        vector<float> y_medium = {sp_value[6+(iCuts*45)][0], sp_value[7+(iCuts*45)][0], sp_value[8+(iCuts*45)][0], sp_value[9+(iCuts*45)][0], sp_value[10+(iCuts*45)][0], sp_value[11+(iCuts*45)][0]};
        vector<float> y_tight = {sp_value[18+(iCuts*45)][0], sp_value[19+(iCuts*45)][0], sp_value[20+(iCuts*45)][0], sp_value[21+(iCuts*45)][0], sp_value[22+(iCuts*45)][0], sp_value[23+(iCuts*45)][0]};
        vector<float> y_mva80 = {sp_value[27+(iCuts*45)][0], sp_value[28+(iCuts*45)][0], sp_value[29+(iCuts*45)][0], sp_value[30+(iCuts*45)][0], sp_value[31+(iCuts*45)][0], sp_value[32+(iCuts*45)][0]};
        vector<float> y_mva90 = {sp_value[36+(iCuts*45)][0], sp_value[37+(iCuts*45)][0], sp_value[38+(iCuts*45)][0], sp_value[39+(iCuts*45)][0], sp_value[40+(iCuts*45)][0], sp_value[41+(iCuts*45)][0]};

        vector<float> y_unc_l = {sp_value[0+(iCuts*45)][1], sp_value[1+(iCuts*45)][1], sp_value[2+(iCuts*45)][1], sp_value[3+(iCuts*45)][1], sp_value[4+(iCuts*45)][1], sp_value[5+(iCuts*45)][1]};
        vector<float> y_unc_m = {sp_value[6+(iCuts*45)][1], sp_value[7+(iCuts*45)][1], sp_value[8+(iCuts*45)][1], sp_value[9+(iCuts*45)][1], sp_value[10+(iCuts*45)][1], sp_value[11+(iCuts*45)][1]};
        vector<float> y_unc_ti = {sp_value[18+(iCuts*45)][1], sp_value[19+(iCuts*45)][1], sp_value[20+(iCuts*45)][1], sp_value[21+(iCuts*45)][1], sp_value[22+(iCuts*45)][1], sp_value[23+(iCuts*45)][1]};
        vector<float> y_unc_80 = {sp_value[27+(iCuts*45)][1], sp_value[28+(iCuts*45)][1], sp_value[29+(iCuts*45)][1], sp_value[30+(iCuts*45)][1], sp_value[31+(iCuts*45)][1], sp_value[32+(iCuts*45)][1]};
        vector<float> y_unc_90 = {sp_value[36+(iCuts*45)][1], sp_value[37+(iCuts*45)][1], sp_value[38+(iCuts*45)][1], sp_value[39+(iCuts*45)][1], sp_value[40+(iCuts*45)][1], sp_value[41+(iCuts*45)][1]};

        //--------Get the value of SF--------//
        cout << Cuts[iCuts] << ", looseID :  " << y_loose[0] << ",  " << y_loose[1] <<  ",  " <<  y_loose[2] <<  ",  " <<  y_loose[3] <<  ",  " <<  y_loose[4] <<  ",  " <<  y_loose[5] << endl;
        cout << Cuts[iCuts] << ", mediumID:  " << y_medium[0] << ",  " << y_medium[1] <<  ",  " <<  y_medium[2] <<  ",  " <<  y_medium[3] <<  ",  " <<  y_medium[4] <<  ",  " <<  y_medium[5] << endl;
        cout << Cuts[iCuts] << ", tightID:  " << y_tight[0] << ",  " << y_tight[1] <<  ",  " <<  y_tight[2] <<  ",  " <<  y_tight[3] <<  ",  " <<  y_tight[4] <<  ",  " <<  y_tight[5] << endl;
        cout << Cuts[iCuts] << ", mva80ID:  " << y_mva80[0] << ",  " << y_mva80[1] <<  ",  " <<  y_mva80[2] <<  ",  " <<  y_mva80[3] <<  ",  " <<  y_mva80[4] <<  ",  " <<  y_mva80[5] << endl;
        cout << Cuts[iCuts] << ", mva90ID:  " << y_mva90[0] << ",  " << y_mva90[1] <<  ",  " <<  y_mva90[2] <<  ",  " <<  y_mva90[3] <<  ",  " <<  y_mva90[4] <<  ",  " <<  y_mva90[5] << endl;

        vector<vector<float>> y_value = {y_loose, y_medium, y_tight, y_mva80, y_mva90};
        vector<vector<float>> unc_value = {y_unc_l, y_unc_m, y_unc_ti, y_unc_80, y_unc_90};

        for (int iIDs = 0; iIDs < IDs.size(); iIDs++){
            vhsp[iCuts].push_back(new TH1F(Form("hsp_%s_%s", Cuts[iCuts].c_str(), IDs[iIDs].c_str()),"",6,0,6));
            vhsp[iCuts][iIDs]->SetStats(0);

            for (int icat = 0; icat < category.size(); icat++){
                vhsp[iCuts][iIDs]->SetBinContent(icat+1, y_value[iIDs][icat]);
                vhsp[iCuts][iIDs]->SetBinError(icat+1, unc_value[iIDs][icat]);
                vhsp[iCuts][iIDs]->GetXaxis()->SetBinLabel(icat+1, category[icat].c_str());
            }
            vhsp[iCuts][iIDs]->GetYaxis()->SetTitle("Scale Factor");
            vhsp[iCuts][iIDs]->GetYaxis()->SetRangeUser(0.79, 1.05);
            vhsp[iCuts][iIDs]->SetMarkerStyle(20);
            // if (iIDs >= 3) vhsp[iCuts][iIDs]->SetMarkerStyle(24);
            vhsp[iCuts][iIDs]->SetMarkerColor(TColor::GetColor(ColorName_sp[iIDs].c_str()));
            vhsp[iCuts][iIDs]->SetMarkerSize(1.5);
            vhsp[iCuts][iIDs]->SetLineColor(TColor::GetColor(ColorName_sp[iIDs].c_str()));
            vhsp[iCuts][iIDs]->SetLineWidth(2);
            cout << "here" << endl;
            if (iIDs == 0) vhsp[iCuts][iIDs]->Draw("le1p");
            else vhsp[iCuts][iIDs]->Draw("le1p same");

            legend->AddEntry(vhsp[iCuts][iIDs], Form("%s_ID", IDs[iIDs].c_str()), "LE1P");

            // vhsp[iCuts][iIDs]->Write();
        }
        legend->SetTextSize(0.03);
        legend->SetBorderSize(0);
        legend->SetNColumns(2);
        legend->Draw();

        TLatex *ltx = new TLatex();
        string information = Form("%s", Cuts[iCuts].c_str());
        ltx->SetNDC();
        ltx->SetTextFont(42);
        ltx->SetTextSize(0.04);
        ltx->DrawLatex(0.15, 0.25, information.c_str());

        TLatex* CMS = new TLatex();
        CMS->SetNDC(true);
        CMS->SetTextFont(42);
        CMS->SetTextSize(0.04);
        CMS->DrawLatex(0.17, 0.91, "#bf{CMS} #it{work-in-progress}");

        TLatex* lumi = new TLatex();
        lumi->SetNDC();
        lumi->SetTextFont(42);
        lumi->SetTextSize(0.03);
        // lumi->DrawLatex(0.69, 0.91, "27.007 fb^{-1} (13.8TeV)");
        lumi->DrawLatex(0.69, 0.91, "8.1747 fb^{-1} (13.8TeV)");

        canvas->Print(Form("./try_NanoAODv12_preEE/69200/sf/SummaryPlot_%s_%s.pdf", Cuts[iCuts].c_str(), pu.c_str()));
        canvas->Close();

        // fsf_out->Close();
    }
}
float err_calc(float nominal_sf, float up_sf, float down_sf, float nominal_unc, float nominal_Zgttb_sf){
    float pu_syserr = (abs(up_sf - nominal_sf) + abs(down_sf - nominal_sf))/2/nominal_sf;
    float Zgttb_syserr = fabs(nominal_sf - nominal_Zgttb_sf) / nominal_sf;
    float toterr = sqrt(pu_syserr*pu_syserr + Zgttb_syserr*Zgttb_syserr + nominal_unc*nominal_unc);
    return toterr;
}
float syserr_calc(float nominal_sf, float up_sf, float down_sf, float nominal_Zgttb_sf){
    float pu_syserr = (abs(up_sf - nominal_sf) + abs(down_sf - nominal_sf))/2/nominal_sf;
    float Zgttb_syserr = fabs(nominal_sf - nominal_Zgttb_sf) / nominal_sf;
    float syserr = sqrt(pu_syserr*pu_syserr + Zgttb_syserr*Zgttb_syserr);
    return syserr;
}
float ttbar_syserr_calc(float nominal_sf, float nominal_Zgttb_sf){
    float Zgttb_syserr = fabs(nominal_sf - nominal_Zgttb_sf) / nominal_sf;
    return Zgttb_syserr;
}

void DrawNominal(vector<vector<float>> nominal_sp_value, vector<vector<float>> up_sp_value, vector<vector<float>> down_sp_value, vector<vector<float>> nominal_sp_Zgttbar_value){

    vector<vector<TH1F*>> vhsp;
    for (int iCuts = 0; iCuts < Cuts.size(); iCuts++){

        TFile *fsf_out = TFile::Open(Form("./try_NanoAODv12_preEE/%s_SummaryPlot_22.root", Cuts[iCuts].c_str()), "RECREATE");
        fsf_out->cd();

        vhsp.push_back(vector<TH1F*>());

        gStyle->SetPadTickX(1);
        gStyle->SetPadTickY(1);
        gStyle->SetOptStat(0);
        TCanvas *canvas = new TCanvas("c", "c", 800, 800);
        canvas->SetBottomMargin(0.06);
        canvas->SetRightMargin(0.05);
        // canvas->SetTopMargin(0.095);
        canvas->SetLeftMargin(0.13);

        TLegend *legend = new TLegend(0.15, 0.09, 0.55, 0.22); //(x1,y1,x2,y2)

        vector<float> nominal_y_loose = {nominal_sp_value[0+(iCuts*45)][0], nominal_sp_value[1+(iCuts*45)][0], nominal_sp_value[2+(iCuts*45)][0], nominal_sp_value[3+(iCuts*45)][0], nominal_sp_value[4+(iCuts*45)][0], nominal_sp_value[5+(iCuts*45)][0]};
        vector<float> nominal_y_medium = {nominal_sp_value[9+(iCuts*45)][0], nominal_sp_value[10+(iCuts*45)][0], nominal_sp_value[11+(iCuts*45)][0], nominal_sp_value[12+(iCuts*45)][0], nominal_sp_value[13+(iCuts*45)][0], nominal_sp_value[14+(iCuts*45)][0]};
        vector<float> nominal_y_tight = {nominal_sp_value[18+(iCuts*45)][0], nominal_sp_value[19+(iCuts*45)][0], nominal_sp_value[20+(iCuts*45)][0], nominal_sp_value[21+(iCuts*45)][0], nominal_sp_value[22+(iCuts*45)][0], nominal_sp_value[23+(iCuts*45)][0]};
        vector<float> nominal_y_mva80 = {nominal_sp_value[27+(iCuts*45)][0], nominal_sp_value[28+(iCuts*45)][0], nominal_sp_value[29+(iCuts*45)][0], nominal_sp_value[30+(iCuts*45)][0], nominal_sp_value[31+(iCuts*45)][0], nominal_sp_value[32+(iCuts*45)][0]};
        vector<float> nominal_y_mva90 = {nominal_sp_value[36+(iCuts*45)][0], nominal_sp_value[37+(iCuts*45)][0], nominal_sp_value[38+(iCuts*45)][0], nominal_sp_value[39+(iCuts*45)][0], nominal_sp_value[40+(iCuts*45)][0], nominal_sp_value[41+(iCuts*45)][0]};

        vector<float> up_y_loose = {up_sp_value[0+(iCuts*45)][0], up_sp_value[1+(iCuts*45)][0], up_sp_value[2+(iCuts*45)][0], up_sp_value[3+(iCuts*45)][0], up_sp_value[4+(iCuts*45)][0], up_sp_value[5+(iCuts*45)][0]};
        vector<float> up_y_medium = {up_sp_value[9+(iCuts*45)][0], up_sp_value[10+(iCuts*45)][0], up_sp_value[11+(iCuts*45)][0], up_sp_value[12+(iCuts*45)][0], up_sp_value[13+(iCuts*45)][0], up_sp_value[14+(iCuts*45)][0]};
        vector<float> up_y_tight = {up_sp_value[18+(iCuts*45)][0], up_sp_value[19+(iCuts*45)][0], up_sp_value[20+(iCuts*45)][0], up_sp_value[21+(iCuts*45)][0], up_sp_value[22+(iCuts*45)][0], up_sp_value[23+(iCuts*45)][0]};
        vector<float> up_y_mva80 = {up_sp_value[27+(iCuts*45)][0], up_sp_value[28+(iCuts*45)][0], up_sp_value[29+(iCuts*45)][0], up_sp_value[30+(iCuts*45)][0], up_sp_value[31+(iCuts*45)][0], up_sp_value[32+(iCuts*45)][0]};
        vector<float> up_y_mva90 = {up_sp_value[36+(iCuts*45)][0], up_sp_value[37+(iCuts*45)][0], up_sp_value[38+(iCuts*45)][0], up_sp_value[39+(iCuts*45)][0], up_sp_value[40+(iCuts*45)][0], up_sp_value[41+(iCuts*45)][0]};

        vector<float> down_y_loose = {down_sp_value[0+(iCuts*45)][0], down_sp_value[1+(iCuts*45)][0], down_sp_value[2+(iCuts*45)][0], down_sp_value[3+(iCuts*45)][0], down_sp_value[4+(iCuts*45)][0], down_sp_value[5+(iCuts*45)][0]};
        vector<float> down_y_medium = {down_sp_value[9+(iCuts*45)][0], down_sp_value[10+(iCuts*45)][0], down_sp_value[11+(iCuts*45)][0], down_sp_value[12+(iCuts*45)][0], down_sp_value[13+(iCuts*45)][0], down_sp_value[14+(iCuts*45)][0]};
        vector<float> down_y_tight = {down_sp_value[18+(iCuts*45)][0], down_sp_value[19+(iCuts*45)][0], down_sp_value[20+(iCuts*45)][0], down_sp_value[21+(iCuts*45)][0], down_sp_value[22+(iCuts*45)][0], down_sp_value[23+(iCuts*45)][0]};
        vector<float> down_y_mva80 = {down_sp_value[27+(iCuts*45)][0], down_sp_value[28+(iCuts*45)][0], down_sp_value[29+(iCuts*45)][0], down_sp_value[30+(iCuts*45)][0], down_sp_value[31+(iCuts*45)][0], down_sp_value[32+(iCuts*45)][0]};
        vector<float> down_y_mva90 = {down_sp_value[36+(iCuts*45)][0], down_sp_value[37+(iCuts*45)][0], down_sp_value[38+(iCuts*45)][0], down_sp_value[39+(iCuts*45)][0], down_sp_value[40+(iCuts*45)][0], down_sp_value[41+(iCuts*45)][0]};

        vector<float> nominal_y_unc_l = {nominal_sp_value[0+(iCuts*45)][1], nominal_sp_value[1+(iCuts*45)][1], nominal_sp_value[2+(iCuts*45)][1], nominal_sp_value[3+(iCuts*45)][1], nominal_sp_value[4+(iCuts*45)][1], nominal_sp_value[5+(iCuts*45)][1]};
        vector<float> nominal_y_unc_m = {nominal_sp_value[9+(iCuts*45)][1], nominal_sp_value[10+(iCuts*45)][1], nominal_sp_value[11+(iCuts*45)][1], nominal_sp_value[12+(iCuts*45)][1], nominal_sp_value[13+(iCuts*45)][1], nominal_sp_value[14+(iCuts*45)][1]};
        vector<float> nominal_y_unc_ti = {nominal_sp_value[18+(iCuts*45)][1], nominal_sp_value[19+(iCuts*45)][1], nominal_sp_value[20+(iCuts*45)][1], nominal_sp_value[21+(iCuts*45)][1], nominal_sp_value[22+(iCuts*45)][1], nominal_sp_value[23+(iCuts*45)][1]};
        vector<float> nominal_y_unc_80 = {nominal_sp_value[27+(iCuts*45)][1], nominal_sp_value[28+(iCuts*45)][1], nominal_sp_value[29+(iCuts*45)][1], nominal_sp_value[30+(iCuts*45)][1], nominal_sp_value[31+(iCuts*45)][1], nominal_sp_value[32+(iCuts*45)][1]};
        vector<float> nominal_y_unc_90 = {nominal_sp_value[36+(iCuts*45)][1], nominal_sp_value[37+(iCuts*45)][1], nominal_sp_value[38+(iCuts*45)][1], nominal_sp_value[39+(iCuts*45)][1], nominal_sp_value[40+(iCuts*45)][1], nominal_sp_value[41+(iCuts*45)][1]};

        vector<float> nominal_Zgttb_y_l = {nominal_sp_Zgttbar_value[0+(iCuts*45)][0], nominal_sp_Zgttbar_value[1+(iCuts*45)][0], nominal_sp_Zgttbar_value[2+(iCuts*45)][0], nominal_sp_Zgttbar_value[3+(iCuts*45)][0], nominal_sp_Zgttbar_value[4+(iCuts*45)][0], nominal_sp_Zgttbar_value[5+(iCuts*45)][0]};
        vector<float> nominal_Zgttb_y_m = {nominal_sp_Zgttbar_value[9+(iCuts*45)][0], nominal_sp_Zgttbar_value[10+(iCuts*45)][0], nominal_sp_Zgttbar_value[11+(iCuts*45)][0], nominal_sp_Zgttbar_value[12+(iCuts*45)][0], nominal_sp_Zgttbar_value[13+(iCuts*45)][0], nominal_sp_Zgttbar_value[14+(iCuts*45)][0]};
        vector<float> nominal_Zgttb_y_ti = {nominal_sp_Zgttbar_value[18+(iCuts*45)][0], nominal_sp_Zgttbar_value[19+(iCuts*45)][0], nominal_sp_Zgttbar_value[20+(iCuts*45)][0], nominal_sp_Zgttbar_value[21+(iCuts*45)][0], nominal_sp_Zgttbar_value[22+(iCuts*45)][0], nominal_sp_Zgttbar_value[23+(iCuts*45)][0]};
        vector<float> nominal_Zgttb_y_80 = {nominal_sp_Zgttbar_value[27+(iCuts*45)][0], nominal_sp_Zgttbar_value[28+(iCuts*45)][0], nominal_sp_Zgttbar_value[29+(iCuts*45)][0], nominal_sp_Zgttbar_value[30+(iCuts*45)][0], nominal_sp_Zgttbar_value[31+(iCuts*45)][0], nominal_sp_Zgttbar_value[32+(iCuts*45)][0]};
        vector<float> nominal_Zgttb_y_90 = {nominal_sp_Zgttbar_value[36+(iCuts*45)][0], nominal_sp_Zgttbar_value[37+(iCuts*45)][0], nominal_sp_Zgttbar_value[38+(iCuts*45)][0], nominal_sp_Zgttbar_value[39+(iCuts*45)][0], nominal_sp_Zgttbar_value[40+(iCuts*45)][0], nominal_sp_Zgttbar_value[41+(iCuts*45)][0]};

        vector<float> y_unc_l = {};
        vector<float> y_unc_m = {};
        vector<float> y_unc_ti = {};
        vector<float> y_unc_80 = {};
        vector<float> y_unc_90 = {};
        for (int i=0; i<nominal_y_loose.size(); i++){
            y_unc_l.push_back(err_calc(nominal_y_loose[i], up_y_loose[i], down_y_loose[i], nominal_y_unc_l[i], nominal_Zgttb_y_l[i]));
            y_unc_m.push_back(err_calc(nominal_y_medium[i], up_y_medium[i], down_y_medium[i], nominal_y_unc_m[i], nominal_Zgttb_y_m[i]));
            y_unc_ti.push_back(err_calc(nominal_y_tight[i], up_y_tight[i], down_y_tight[i], nominal_y_unc_ti[i], nominal_Zgttb_y_ti[i]));
            y_unc_80.push_back(err_calc(nominal_y_mva80[i], up_y_mva80[i], down_y_mva80[i], nominal_y_unc_80[i], nominal_Zgttb_y_80[i]));
            y_unc_90.push_back(err_calc(nominal_y_mva90[i], up_y_mva90[i], down_y_mva90[i], nominal_y_unc_90[i], nominal_Zgttb_y_90[i]));
        }
        cout << Cuts[iCuts] << ", looseID_staterr :  " << nominal_y_unc_l[0] << ",  " << nominal_y_unc_l[1] <<  ",  " <<  nominal_y_unc_l[2] <<  ",  " <<  nominal_y_unc_l[3] <<  ",  " <<  nominal_y_unc_l[4] <<  ",  " <<  nominal_y_unc_l[5] << endl;
        cout << Cuts[iCuts] << ", mediumID_staterr:  " << nominal_y_unc_m[0] << ",  " << nominal_y_unc_m[1] <<  ",  " <<  nominal_y_unc_m[2] <<  ",  " <<  nominal_y_unc_m[3] <<  ",  " <<  nominal_y_unc_m[4] <<  ",  " <<  nominal_y_unc_m[5] << endl;
        cout << Cuts[iCuts] << ", tightID_staterr:  " << nominal_y_unc_ti[0] << ",  " << nominal_y_unc_ti[1] <<  ",  " <<  nominal_y_unc_ti[2] <<  ",  " <<  nominal_y_unc_ti[3] <<  ",  " <<  nominal_y_unc_ti[4] <<  ",  " <<  nominal_y_unc_ti[5] << endl;
        cout << Cuts[iCuts] << ", mva80ID_staterr:  " << nominal_y_unc_80[0] << ",  " << nominal_y_unc_80[1] <<  ",  " <<  nominal_y_unc_80[2] <<  ",  " <<  nominal_y_unc_80[3] <<  ",  " <<  nominal_y_unc_80[4] <<  ",  " <<  nominal_y_unc_80[5] << endl;
        cout << Cuts[iCuts] << ", mva90ID_staterr:  " << nominal_y_unc_90[0] << ",  " << nominal_y_unc_90[1] <<  ",  " <<  nominal_y_unc_90[2] <<  ",  " <<  nominal_y_unc_90[3] <<  ",  " <<  nominal_y_unc_90[4] <<  ",  " <<  nominal_y_unc_90[5] << endl;
        
        //--------Print the value of total syserr--------//
        vector<float> sys_unc_l = {};
        vector<float> sys_unc_m = {};
        vector<float> sys_unc_ti = {};
        vector<float> sys_unc_80 = {};
        vector<float> sys_unc_90 = {};
        for (int i=0; i<nominal_y_loose.size(); i++){
            sys_unc_l.push_back(syserr_calc(nominal_y_loose[i], up_y_loose[i], down_y_loose[i], nominal_Zgttb_y_l[i]));
            sys_unc_m.push_back(syserr_calc(nominal_y_medium[i], up_y_medium[i], down_y_medium[i], nominal_Zgttb_y_m[i]));
            sys_unc_ti.push_back(syserr_calc(nominal_y_tight[i], up_y_tight[i], down_y_tight[i], nominal_Zgttb_y_ti[i]));
            sys_unc_80.push_back(syserr_calc(nominal_y_mva80[i], up_y_mva80[i], down_y_mva80[i], nominal_Zgttb_y_80[i]));
            sys_unc_90.push_back(syserr_calc(nominal_y_mva90[i], up_y_mva90[i], down_y_mva90[i], nominal_Zgttb_y_90[i]));
        }
        cout << Cuts[iCuts] << ", looseID_systerr :  " << sys_unc_l[0] << ",  " << sys_unc_l[1] <<  ",  " <<  sys_unc_l[2] <<  ",  " <<  sys_unc_l[3] <<  ",  " <<  sys_unc_l[4] <<  ",  " <<  sys_unc_l[5] << endl;
        cout << Cuts[iCuts] << ", mediumID_systerr:  " << sys_unc_m[0] << ",  " << sys_unc_m[1] <<  ",  " <<  sys_unc_m[2] <<  ",  " <<  sys_unc_m[3] <<  ",  " <<  sys_unc_m[4] <<  ",  " <<  sys_unc_m[5] << endl;
        cout << Cuts[iCuts] << ", tightID_systerr:  " << sys_unc_ti[0] << ",  " << sys_unc_ti[1] <<  ",  " <<  sys_unc_ti[2] <<  ",  " <<  sys_unc_ti[3] <<  ",  " <<  sys_unc_ti[4] <<  ",  " <<  sys_unc_ti[5] << endl;
        cout << Cuts[iCuts] << ", mva80ID_systerr:  " << sys_unc_80[0] << ",  " << sys_unc_80[1] <<  ",  " <<  sys_unc_80[2] <<  ",  " <<  sys_unc_80[3] <<  ",  " <<  sys_unc_80[4] <<  ",  " <<  sys_unc_80[5] << endl;
        cout << Cuts[iCuts] << ", mva90ID_systerr:  " << sys_unc_90[0] << ",  " << sys_unc_90[1] <<  ",  " <<  sys_unc_90[2] <<  ",  " <<  sys_unc_90[3] <<  ",  " <<  sys_unc_90[4] <<  ",  " <<  sys_unc_90[5] << endl;
        //--------Print the value of ttbar syserr--------//
        vector<float> ttbar_sys_unc_l = {};
        vector<float> ttbar_sys_unc_m = {};
        vector<float> ttbar_sys_unc_ti = {};
        vector<float> ttbar_sys_unc_80 = {};
        vector<float> ttbar_sys_unc_90 = {};
        for (int i=0; i<nominal_y_loose.size(); i++){
            ttbar_sys_unc_l.push_back(ttbar_syserr_calc(nominal_y_loose[i], nominal_Zgttb_y_l[i]));
            ttbar_sys_unc_m.push_back(ttbar_syserr_calc(nominal_y_medium[i], nominal_Zgttb_y_m[i]));
            ttbar_sys_unc_ti.push_back(ttbar_syserr_calc(nominal_y_tight[i], nominal_Zgttb_y_ti[i]));
            ttbar_sys_unc_80.push_back(ttbar_syserr_calc(nominal_y_mva80[i], nominal_Zgttb_y_80[i]));
            ttbar_sys_unc_90.push_back(ttbar_syserr_calc(nominal_y_mva90[i], nominal_Zgttb_y_90[i]));
        }
        cout << Cuts[iCuts] << ", looseID_ttbar_syserr :  " << ttbar_sys_unc_l[0] << ",  " << ttbar_sys_unc_l[1] <<  ",  " <<  ttbar_sys_unc_l[2] <<  ",  " <<  ttbar_sys_unc_l[3] <<  ",  " <<  ttbar_sys_unc_l[4] <<  ",  " <<  ttbar_sys_unc_l[5] << endl;
        cout << Cuts[iCuts] << ", mediumID_ttbar_syserr:  " << ttbar_sys_unc_m[0] << ",  " << ttbar_sys_unc_m[1] <<  ",  " <<  ttbar_sys_unc_m[2] <<  ",  " <<  ttbar_sys_unc_m[3] <<  ",  " <<  ttbar_sys_unc_m[4] <<  ",  " <<  ttbar_sys_unc_m[5] << endl;
        cout << Cuts[iCuts] << ", tightID_ttbar_syserr:  " << ttbar_sys_unc_ti[0] << ",  " << ttbar_sys_unc_ti[1] <<  ",  " <<  ttbar_sys_unc_ti[2] <<  ",  " <<  ttbar_sys_unc_ti[3] <<  ",  " <<  ttbar_sys_unc_ti[4] <<  ",  " <<  ttbar_sys_unc_ti[5] << endl;
        cout << Cuts[iCuts] << ", mva80ID_ttbar_syserr:  " << ttbar_sys_unc_80[0] << ",  " << ttbar_sys_unc_80[1] <<  ",  " <<  ttbar_sys_unc_80[2] <<  ",  " <<  ttbar_sys_unc_80[3] <<  ",  " <<  ttbar_sys_unc_80[4] <<  ",  " <<  ttbar_sys_unc_80[5] << endl;
        cout << Cuts[iCuts] << ", mva90ID_ttbar_syserr:  " << ttbar_sys_unc_90[0] << ",  " << ttbar_sys_unc_90[1] <<  ",  " <<  ttbar_sys_unc_90[2] <<  ",  " <<  ttbar_sys_unc_90[3] <<  ",  " <<  ttbar_sys_unc_90[4] <<  ",  " <<  ttbar_sys_unc_90[5] << endl;
        
        //--------Get the value of SF--------//
        vector<vector<float>> y_value = {nominal_y_loose, nominal_y_medium, nominal_y_tight, nominal_y_mva80, nominal_y_mva90};
        vector<vector<float>> unc_value = {y_unc_l, y_unc_m, y_unc_ti, y_unc_80, y_unc_90};

        for (int iIDs = 0; iIDs < IDs.size(); iIDs++){
            vhsp[iCuts].push_back(new TH1F(Form("hsp_%s_%s", Cuts[iCuts].c_str(), IDs[iIDs].c_str()),"",6,0,6));
            vhsp[iCuts][iIDs]->SetStats(0);

            for (int icat = 0; icat < category.size(); icat++){
                vhsp[iCuts][iIDs]->SetBinContent(icat+1, y_value[iIDs][icat]);
                vhsp[iCuts][iIDs]->SetBinError(icat+1, unc_value[iIDs][icat]);
                vhsp[iCuts][iIDs]->GetXaxis()->SetBinLabel(icat+1, category[icat].c_str());
            }
            vhsp[iCuts][iIDs]->GetYaxis()->SetTitle("Scale Factor");
            vhsp[iCuts][iIDs]->GetYaxis()->SetRangeUser(0.79, 1.05);
            vhsp[iCuts][iIDs]->SetMarkerStyle(20);
            // if (iIDs >= 3) vhsp[iCuts][iIDs]->SetMarkerStyle(24);
            vhsp[iCuts][iIDs]->SetMarkerColor(TColor::GetColor(ColorName_sp[iIDs].c_str()));
            vhsp[iCuts][iIDs]->SetMarkerSize(1.5);
            vhsp[iCuts][iIDs]->SetLineColor(TColor::GetColor(ColorName_sp[iIDs].c_str()));
            vhsp[iCuts][iIDs]->SetLineWidth(2);
            cout << "here" << endl;
            if (iIDs == 0) vhsp[iCuts][iIDs]->Draw("le1p");
            else vhsp[iCuts][iIDs]->Draw("le1p same");

            legend->AddEntry(vhsp[iCuts][iIDs], Form("%s_ID", IDs[iIDs].c_str()), "LE1P");

            vhsp[iCuts][iIDs]->Write();
        }
        legend->SetTextSize(0.03);
        legend->SetBorderSize(0);
        legend->SetNColumns(2);
        legend->Draw();

        TLatex *ltx = new TLatex();
        string information = Form("%s", Cuts[iCuts].c_str());
        ltx->SetNDC();
        ltx->SetTextFont(42);
        ltx->SetTextSize(0.04);
        ltx->DrawLatex(0.15, 0.25, information.c_str());

        TLatex* CMS = new TLatex();
        CMS->SetNDC(true);
        CMS->SetTextFont(42);
        CMS->SetTextSize(0.04);
        CMS->DrawLatex(0.17, 0.91, "#bf{CMS} #it{work-in-progress}");

        TLatex* lumi = new TLatex();
        lumi->SetNDC();
        lumi->SetTextFont(42);
        lumi->SetTextSize(0.03);
        // lumi->DrawLatex(0.69, 0.91, "27.007 fb^{-1} (13.8TeV)");
        lumi->DrawLatex(0.69, 0.91, "8.1747 fb^{-1} (13.8TeV)");

        canvas->Print(Form("./try_NanoAODv12_preEE/69200/sf/SummaryPlot_%s_withsys.pdf", Cuts[iCuts].c_str()));
        canvas->Close();

        fsf_out->Close();
    }
}

vector<TEfficiency*> getHist(string infile_data, string infile_mc, int IDType, int CutType, int R9Type, int region){

    vector<TFile*> vfiles; vfiles.clear();
    vfiles.push_back(new TFile(infile_data.c_str()));
    vfiles.push_back(new TFile(infile_mc.c_str()));

    vector<double> bin_sp = {-TMath::Pi()-1, TMath::Pi()+1};
    int n_bin_sp =  bin_sp.size() - 1;
    vector<TEfficiency*> eff_sp; eff_sp.clear();

    vector<string> fileName = {"Prompt", "mc"};
    for (int ifile = 0; ifile < vfiles.size(); ifile++){

        TTree *inTree = (TTree*) vfiles[ifile]->Get("tree");
        eff_sp.push_back(new TEfficiency(Form("eff_sp_%s", fileName[ifile].c_str()), " ", n_bin_sp, bin_sp.data()));

        Double_t pho_phi, weights;
        Double_t Cutbased, mva, MVA_WP80, MVA_WP90;
        Double_t HasPix, CSEV;
        Double_t pho_r9, pho_EB, pho_EE;
        inTree->SetBranchAddress("pho_phi", &pho_phi);
        inTree->SetBranchAddress("weights", &weights);
        inTree->SetBranchAddress("pho_cutBased", &Cutbased);
        inTree->SetBranchAddress("pho_mvaID", &mva);
        inTree->SetBranchAddress("pho_mvaID_WP80", &MVA_WP80);
        inTree->SetBranchAddress("pho_mvaID_WP90", &MVA_WP90);
        inTree->SetBranchAddress("pho_pixelseed", &HasPix);
        inTree->SetBranchAddress("pho_electronveto", &CSEV);
        inTree->SetBranchAddress("pho_r9", &pho_r9);
        inTree->SetBranchAddress("pho_EB", &pho_EB);
        inTree->SetBranchAddress("pho_EE", &pho_EE);

        for (Long64_t ev=0; ev < inTree->GetEntriesFast(); ev++){
            inTree->GetEntry(ev);

            bool cut;
            if (CutType == 0) cut = (HasPix == 0); //HasPix == 0, pass pixel seed veto
            if (CutType == 1) cut = (CSEV == 1);  //CSEV == 1, pass conversion safe electron veto

            bool passID;
            if (IDType == 0) passID = (Cutbased >= 1); //CutBased = 1: Loose
            if (IDType == 1) passID = (Cutbased >= 2); //CutBased = 2: Medium
            if (IDType == 2) passID = (Cutbased == 3); //CutBased = 3: Tight
            if (IDType == 3) passID = (MVA_WP80 == 1);
            if (IDType == 4) passID = (MVA_WP90 == 1);
            if (!passID) continue;

            bool cat;
            if (R9Type == 0 && region == 0) cat = (pho_EB > 0.5); // > 0.5 means ==1
            if (R9Type == 1 && region == 0) cat = (pho_r9 > 0.96) & (pho_EB > 0.5);
            if (R9Type == 2 && region == 0) cat = (pho_r9 < 0.96) & (pho_EB > 0.5);
            if (R9Type == 0 && region == 1) cat = (pho_EE > 0.5);
            if (R9Type == 1 && region == 1) cat = (pho_r9 > 0.96) & (pho_EE > 0.5);
            if (R9Type == 2 && region == 1) cat = (pho_r9 < 0.96) & (pho_EE > 0.5);
            if (R9Type == 0 && region == 2) cat = true;
            if (R9Type == 1 && region == 2) cat = (pho_r9 > 0.96);
            if (R9Type == 2 && region == 2) cat = (pho_r9 < 0.96);
            if (!cat) continue;

            eff_sp[ifile]->FillWeighted(cut, weights, pho_phi);
        }
    }
    return eff_sp;
}

void pu_syserr(){
    vector<vector<float>> nominal_sp_value = {};
    vector<vector<float>> up_sp_value = {};
    vector<vector<float>> down_sp_value = {};
    vector<vector<float>> nominal_sp_Zgttbar_value = {};
    for (int iCuts = 0; iCuts < Cuts.size(); iCuts++){
        for (int iIDs = 0; iIDs < IDs.size(); iIDs++){
            for (int iRegion = 0; iRegion < Region.size(); iRegion++){
                for (int iR9 = 0; iR9 < R9.size(); iR9++){
                    // getHist("data_try_acc_singlemu.root", "mc_try_acc_singlemu.root", iIDs, iCuts, iR9, iRegion);
                    //----------Create value of summaryplot----------//
                    auto nominal_sf_Zg = SummaryPlot(getHist("./try_NanoAODv12_preEE/69200/data_BCD.root", "./try_NanoAODv12_preEE/69200/DY_preEE_all.root", iIDs, iCuts, iR9, iRegion));
                    //! auto nominal_sf_Zgttb = SummaryPlot(getHist("data.root", "signal.root", iIDs, iCuts, iR9, iRegion));
                    float sf_y =  nominal_sf_Zg[0];
                    float sf_unc = nominal_sf_Zg[1];
                    // float sf_Zgttb_y = nominal_sf_Zgttb[0];
                    //! vector<float> nominal_sf_value = {sf_Zg_y, sf_Zg_unc, sf_Zgttb_y};
                    vector<float> nominal_sf_value = {sf_y, sf_unc};
                    nominal_sp_value.push_back(nominal_sf_value);

                    auto up_sf_all = SummaryPlot(getHist("./try_NanoAODv12_preEE/66000/data_BCD.root", "./try_NanoAODv12_preEE/66000/DY_preEE_all.root", iIDs, iCuts, iR9, iRegion));
                    sf_y =  up_sf_all[0];
                    sf_unc = up_sf_all[1];
                    vector<float> up_sf_value = {sf_y, sf_unc};
                    up_sp_value.push_back(up_sf_value);

                    auto down_sf_all = SummaryPlot(getHist("./try_NanoAODv12_preEE/72400/data_BCD.root", "./try_NanoAODv12_preEE/72400/DY_preEE_all.root", iIDs, iCuts, iR9, iRegion));
                    sf_y =  down_sf_all[0];
                    sf_unc = down_sf_all[1];
                    vector<float> down_sf_value = {sf_y, sf_unc};
                    down_sp_value.push_back(down_sf_value);

                    auto nominal_sf_Zgttb = SummaryPlot(getHist("./try_NanoAODv12_preEE/69200/data_BCD.root", "./try_NanoAODv12_preEE/69200/DY_ttbar_preEE.root", iIDs, iCuts, iR9, iRegion));
                    sf_y =  nominal_sf_Zgttb[0];
                    sf_unc = nominal_sf_Zgttb[1];
                    vector<float> nominal_sf_Zgttbar_value = {sf_y, sf_unc};
                    nominal_sp_Zgttbar_value.push_back(nominal_sf_Zgttbar_value);
                }
            }
        }
    }
    DrawSummaryPlot(nominal_sp_value, "nominal");
    DrawSummaryPlot(up_sp_value, "up");
    DrawSummaryPlot(down_sp_value, "down");

    DrawNominal(nominal_sp_value, up_sp_value, down_sp_value, nominal_sp_Zgttbar_value);
}