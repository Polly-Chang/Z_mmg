vector<string> IDs = {"Loose","Medium","Tight", "MVA80", "MVA90"};
vector<string> Cuts = {"HasPix","CSEV"};
vector<string> R9 = {"InclusiveR9","HR9","LR9"};
vector<string> Region = {"EB","EE",""};
vector<string> category = {"EB Inc.","EB high R9","EB low R9","EE Inc.","EE high R9","EE low R9"};
vector<string> measurement ={"phoPt","phoPhi", "phoEta", "zmass", "mva", "Sceta", "npvs"};
string ColorName[3] = {"#533e2d","#e85d04"};
vector<string> ColorName_sp = {"#303030", "#e21e08", "#228b22", "#4040ff", "#e85d04"};

// TFile *fout = TFile::Open("./SF_eff.root", "RECREATE");

float Photon_SCEta(double pho_eta, double pho_phi, bool pho_EB, bool pho_EE, double PV_x, double PV_y, double PV_z){
  float tg_theta_over_2   = exp(pho_eta);
  float tg_theta          = 2 * tg_theta_over_2 / (1-tg_theta_over_2*tg_theta_over_2);//tan(atan(tg_theta_over_2)*2);//
  float tg_sctheta;

  if (pho_EB){//barrel
    float R              = 130;
    float angle_x0_y0    = 0;
    if      (PV_x>0) angle_x0_y0 = atan(PV_y/PV_x);
    else if (PV_x<0) angle_x0_y0 = M_PI + atan(PV_y/PV_x);
    else if (PV_y>0) angle_x0_y0 = M_PI / 2;
    else             angle_x0_y0 = - M_PI / 2;

    float alpha      = angle_x0_y0 + (M_PI - pho_phi);
    float sin_beta   = sqrt(PV_x*PV_x + PV_y*PV_y) / R * sin(alpha);
    float beta       = abs( asin( sin_beta ) );
    float gamma      = M_PI/2 - alpha - beta;
    float l          = sqrt(R*R + (PV_x*PV_x + PV_y*PV_y) - 2*R*sqrt(PV_x*PV_x + PV_y*PV_y)*cos(gamma));

    float z0_zSC    = l / tg_theta;
    tg_sctheta      = R / (PV_z + z0_zSC);

   
  } else if (pho_EE){//endcap

    float intersection_z = (pho_eta>0)?310:-310;
    float base           = intersection_z - PV_z;//10;
    float r              = base * tg_theta;

    float crystalX       = PV_x + r * cos(pho_phi);
    float crystalY       = PV_y + r * sin(pho_phi);
    tg_sctheta           = sqrt( crystalX*crystalX + crystalY*crystalY ) /intersection_z;
  }
  else return pho_eta;

  float sctheta = atan(tg_sctheta);
  if (sctheta<0) sctheta += M_PI;
  float tg_sctheta_over_2 = tan(sctheta/2);//-sqrt(tg_sctheta*tg_sctheta+1)/tg_sctheta -1/tg_sctheta;
  float SCEta = -log(tg_sctheta_over_2);

  return SCEta;
}

double Error_Propagation(double f, double sigmaA, double A, double sigmaB, double B){
     double error;
     error = f * sqrt((sigmaA/A)*(sigmaA/A)+(sigmaB/B)*(sigmaB/B));
     return error;
}

void DrawEfficiency(vector<TEfficiency*> eff, float xmin, float xmax, int Nmeasure, string XaxisName, int IDType, int CutType, int R9Type, int region, string outName){
    //---------- eff[0]->Prompt, eff[1]->MC ----------//
    eff[0]->SetStatisticOption(TEfficiency::kBUniform);
    eff[0]->SetConfidenceLevel(0.683);
    eff[0]->SetPosteriorMode(1);
    eff[1]->SetStatisticOption(TEfficiency::kBUniform);
    eff[1]->SetConfidenceLevel(0.683);
    eff[1]->SetPosteriorMode(1);
    //---------- Start Drawing ----------//
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("c","c",800,800);
    c->cd();
    //---------Create upper plot (Efficiency)(data & mc)---------//
    TPad* pad1 = new TPad("pad1", " ", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0.05);
    // pad1->SetTopMargin(0.01);
    pad1->SetRightMargin(0.05);
    pad1->SetLeftMargin(0.13);
    pad1->SetBottomMargin(0.03);
    pad1->Draw();             
    pad1->cd();

    eff[0]->SetTitle(";;Efficiency");
    eff[0]->SetMarkerColor(TColor::GetColor(ColorName[0].c_str()));
    eff[0]->SetMarkerSize(1.5);
    eff[0]->SetMarkerStyle(20);
    eff[0]->SetLineColor(TColor::GetColor(ColorName[0].c_str()));
    eff[0]->SetLineWidth(2);
    eff[0]->Draw("AP");
    eff[1]->SetMarkerColor(TColor::GetColor(ColorName[1].c_str()));
    eff[1]->SetMarkerSize(1.5);
    eff[1]->SetMarkerStyle(20);
    eff[1]->SetLineColor(TColor::GetColor(ColorName[1].c_str()));
    eff[1]->SetLineWidth(2);
    eff[1]->Draw("Psame");

    pad1->Update();
    TGraphAsymmErrors *gprompt = eff[0]->GetPaintedGraph();
    gprompt->GetXaxis()->SetLabelOffset(0.05);
    gprompt->GetXaxis()->SetLimits(xmin, xmax);
    // gprompt->GetXaxis()->SetRangeUser(xmin,xmax);
    gprompt->GetXaxis()->SetTickSize(0.03);
    gprompt->GetYaxis()->SetTitleSize(0.07);
    gprompt->GetYaxis()->SetRangeUser(0.2, 1.4);
    gprompt->GetYaxis()->SetTickSize(0.03);
    gprompt->GetYaxis()->SetTitleSize(0.07);
    gprompt->GetYaxis()->SetLabelSize(0.055);
    gprompt->GetYaxis()->SetTitleOffset(0.8);

    TLegend *legend = new TLegend(0.17,0.68,0.88,0.8);//(x1,y1,x2,y2)
    legend->SetNColumns(2);
    legend->SetTextSize(0.05);
    legend->AddEntry(eff[0]->GetPaintedGraph(), "Sep2023 BCD", "LE1P");
    legend->AddEntry(eff[1]->GetPaintedGraph(), "DY_preEE", "LE1P");
    legend->SetFillColor(0);
    legend->SetLineColor(0);
    legend->Draw();

    TLatex *ltx = new TLatex();
    string information = Form("%s, %s, %s ID, %s", Region[region].c_str(), R9[R9Type].c_str(), IDs[IDType].c_str(), Cuts[CutType].c_str());
    if (Nmeasure == 5) information = Form("%s, %s ID, %s", R9[R9Type].c_str(), IDs[IDType].c_str(), Cuts[CutType].c_str());
    ltx->SetNDC();
    ltx->SetTextFont(42);
    ltx->SetTextSize(0.05);
    ltx->DrawLatex(0.15, 0.81, information.c_str());
    
    TLatex* CMS = new TLatex();
    CMS->SetNDC(true);
    CMS->SetTextFont(42);
    CMS->SetTextSize(0.05);
    CMS->DrawLatex(0.13, 0.91, "#bf{CMS} #it{work-in-progress}");
    
    TLatex* lumi = new TLatex();
    lumi->SetNDC();
    lumi->SetTextFont(42);
    lumi->SetTextSize(0.04);
    // lumi->DrawLatex(0.72, 0.91, "27.0072 fb^{-1} (13.6TeV)");
    lumi->DrawLatex(0.72, 0.91, "8.17468 fb^{-1} (13.6TeV)");

    c->cd();

    //---------Create lower plot (SF)---------//
    TPad *pad2 = new TPad("pad2", "", 0, 0, 1, 0.3);
    pad2->SetGridy();
    pad2->SetRightMargin(0.05);
    pad2->SetLeftMargin(0.13);
    pad2->SetTopMargin(0.06);
    pad2->SetBottomMargin(0.4);
    pad2->Draw();
    pad2->cd();

    //---------- SFs calculation using----------//
    TH1F *heff_prompt = (TH1F*) eff[0]->GetCopyPassedHisto();

    for (int i=0; i<heff_prompt->GetNbinsX(); i++){
        cout << heff_prompt->GetBinContent(i+1) << endl;
    }

    heff_prompt->Divide((TH1F*) eff[0]->GetCopyTotalHisto());

    // for (int i=0; i<eff[0]->GetCopyTotalHisto()->GetNbinsX(); i++){
    //     cout << eff[0]->GetCopyTotalHisto()->GetBinContent(i+1) << endl;
    // }
    // for (int i=0; i<eff[1]->GetCopyTotalHisto()->GetNbinsX(); i++){
    //     cout << eff[1]->GetCopyTotalHisto()->GetBinContent(i+1) << endl;
    // }

    // for (int i=0; i<heff_prompt->GetNbinsX(); i++){
    //     cout << heff_prompt->GetBinContent(i+1) << endl;
    // }

    TH1F *heff_MC = (TH1F*) eff[1]->GetCopyPassedHisto();
    heff_MC->Divide((TH1F*) eff[1]->GetCopyTotalHisto());
    TH1F *hSF_prompt = (TH1F*) heff_prompt->Clone();
    hSF_prompt->Divide(heff_MC);

    //---------- For Check ----------//
    // for (int i=0; i<heff_MC->GetNbinsX(); i++){
    //     cout << heff_MC->GetBinContent(i+1) << endl;
    //     // cout << heff_MC->GetBinError(i+1) << endl;
    // }
    
    //---------- Error Propagation ----------//
    int SFNbins_prompt = heff_prompt->GetNbinsX();
    float ey_prompt[SFNbins_prompt];
    for (int i = 0; i < SFNbins_prompt; i++){
        float data_error = TMath::Sqrt(0.5*(pow(eff[0]->GetEfficiencyErrorUp(i+1),2) + pow(eff[0]->GetEfficiencyErrorLow(i+1),2)));
        float MC_error = TMath::Sqrt(0.5*(pow(eff[1]->GetEfficiencyErrorUp(i+1),2) + pow(eff[1]->GetEfficiencyErrorLow(i+1),2)));
        //The same as TGraphAssmError->GetErrorY() method, but if there are missing points, GetErrorY does not work.

        float f = hSF_prompt->GetBinContent(i+1);
        ey_prompt[i] = Error_Propagation(f, data_error , heff_prompt->GetBinContent(i+1), MC_error, heff_MC->GetBinContent(i+1));

        hSF_prompt->SetBinError(i+1, ey_prompt[i]);
    }

    string HistName = Form("SFs_%s_%s_%s_%s_%s", measurement[Nmeasure].c_str(), Region[region].c_str(), R9[R9Type].c_str(), IDs[IDType].c_str(), Cuts[CutType].c_str());

    hSF_prompt->SetTitle("");
    hSF_prompt->GetYaxis()->SetTitle("Scale Factor");
    hSF_prompt->GetXaxis()->SetTitle(XaxisName.c_str());

    hSF_prompt->SetMarkerColor(TColor::GetColor(ColorName[0].c_str()));
    hSF_prompt->SetMarkerSize(1);
    hSF_prompt->SetMarkerStyle(20);
    hSF_prompt->SetLineColor(TColor::GetColor(ColorName[0].c_str()));
    hSF_prompt->SetLineWidth(2);
    hSF_prompt->GetXaxis()->SetRangeUser(xmin,xmax);
    hSF_prompt->GetXaxis()->SetTitleSize(0.16);
    hSF_prompt->GetXaxis()->SetTitleOffset(1);
    hSF_prompt->GetXaxis()->SetLabelSize(0.12);
    hSF_prompt->GetXaxis()->SetLabelOffset(0.03);
    hSF_prompt->GetYaxis()->SetRangeUser(0.6,1.4);
    hSF_prompt->GetYaxis()->SetTitleSize(0.13);
    hSF_prompt->GetYaxis()->SetTitleOffset(0.4);
    hSF_prompt->GetYaxis()->SetLabelSize(0.12);
    hSF_prompt->GetYaxis()->SetNdivisions(502);
    hSF_prompt->Draw("LE1P");

    hSF_prompt->SetName(HistName.c_str());
    c->Print(outName.c_str());
    c->Close();
}
vector<float> SummaryPlot(vector<TEfficiency *> eff){
    //---------- eff[0]->Prompt, eff[1]->MC ----------//
    eff[0]->SetStatisticOption(TEfficiency::kBUniform);
    eff[0]->SetConfidenceLevel(0.683);
    eff[0]->SetPosteriorMode(1);
    eff[1]->SetStatisticOption(TEfficiency::kBUniform);
    eff[1]->SetConfidenceLevel(0.683);
    eff[1]->SetPosteriorMode(1);
    //---------- SFs calculation using Zg ----------//
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

void DrawDistribution(vector<TH1F*> h, float xmin, float xmax, int Nmeasure, bool IDCorrec, string XaxisName, int IDType, string outName){
    //---------- h[0]->Prompt, h[1]->MC ----------//
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);

    TCanvas *canvas = new TCanvas("c", "c", 800, 800);
    canvas->SetBottomMargin(0.12);
    canvas->SetLeftMargin(0.17);
    canvas->SetRightMargin(0.07);
    if (Nmeasure == 0) canvas->SetLogy();

    // h[0]->SetStats(0);
    h[0]->SetLineColor(kBlack);
    h[0]->SetLineWidth(2);
    h[0]->SetMarkerStyle(20);
    h[0]->SetMarkerSize(1.5);
    h[0]->SetMarkerColor(kBlack);
    h[0]->GetXaxis()->SetTitle(XaxisName.c_str());
    h[0]->GetYaxis()->SetTitle("Events");
    h[0]->GetXaxis()->SetTitleSize(0.05);
    h[0]->GetYaxis()->SetTitleSize(0.045);
    float h1_maxbin = h[0]->GetBinContent((h[0]->GetMaximumBin()));
    float h2_maxbin = h[1]->GetBinContent((h[1]->GetMaximumBin()));
    float maxbin = (h1_maxbin > h2_maxbin) ? h1_maxbin : h2_maxbin;
    h[0]->GetYaxis()->SetRangeUser(1, maxbin*1.1);
    h[0]->Draw("LE1P");

    h[1]->Scale(h[0]->Integral(-1, -1)/h[1]->Integral(-1, -1));
    h[1]->SetStats(0);
    h[1]->SetLineColor(kRed);
    h[1]->SetLineWidth(2);
    h[1]->SetMarkerStyle(20);
    h[1]->SetMarkerSize(1.5);
    h[1]->SetMarkerColor(kRed);
    h[1]->SetFillColorAlpha(kRed, 0.1);
    h[1]->Draw("HISTSAME");
    
    string HistName = Form("h_%s_%s", measurement[Nmeasure].c_str(), IDs[IDType].c_str());
    h[0]->SetName(HistName.c_str());

    TLatex *ltx = new TLatex();
    string information = Form("%s ID", IDs[IDType].c_str());
    ltx->SetNDC();
    ltx->SetTextFont(42);
    ltx->SetTextSize(0.04);
    if (IDCorrec) ltx->DrawLatex(0.2, 0.83, information.c_str());

    TLatex* CMS = new TLatex();
    CMS->SetNDC(true);
    CMS->SetTextFont(42);
    CMS->SetTextSize(0.04);
    CMS->DrawLatex(0.17, 0.91, "#bf{CMS} #it{work-in-progress}");
    // CMS->DrawLatex(0.5, 0.915, Form("#bf{#scale[1.2]{%s}}", Xname.c_str()));

    TLegend *legend = new TLegend(0.65, 0.77, 0.85, 0.87);//(x1,y1,x2,y2)
    legend->SetTextSize(0.04);
    legend->AddEntry(h[0], "Sep2023 BCD", "LE1P");
    legend->AddEntry(h[1], "DY_preEE", "LE1P");
    legend->SetBorderSize(0);
    legend->Draw();
    canvas->Print(outName.c_str());
    canvas->Close();

}
void DrawDistribution_wRatio(vector<TH1F *> h, float xmin, float xmax, string XaxisName, int IDType, string outName){
    //---------- h[0]->Prompt, h[1]->MC ----------//
    //---------- Start Drawing ----------//
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetOptStat(0);
    TCanvas* c = new TCanvas("c","c",800,800);
    c->cd();
    //---------Create upper plot---------//
    TPad* pad1 = new TPad("pad1", " ", 0, 0.3, 1, 1.0);
    pad1->SetBottomMargin(0.05);
    // pad1->SetTopMargin(0.01);
    pad1->SetRightMargin(0.05);
    pad1->SetLeftMargin(0.13);
    pad1->SetBottomMargin(0.03);
    pad1->Draw();             
    pad1->cd();

    // h[0]->SetStats(0);
    h[0]->SetLineColor(kBlack);
    h[0]->SetLineWidth(2);
    h[0]->SetMarkerStyle(20);
    h[0]->SetMarkerSize(1.5);
    h[0]->SetMarkerColor(kBlack);
    h[0]->GetXaxis()->SetTitle(XaxisName.c_str());
    h[0]->GetXaxis()->SetLabelOffset(1.5);
    h[0]->GetYaxis()->SetTitle("Events");
    h[0]->GetXaxis()->SetTitleSize(0.05);
    h[0]->GetYaxis()->SetTitleSize(0.045);
    float h1_maxbin = h[0]->GetBinContent((h[0]->GetMaximumBin()));
    float h2_maxbin = h[1]->GetBinContent((h[1]->GetMaximumBin()));
    float maxbin = (h1_maxbin > h2_maxbin) ? h1_maxbin : h2_maxbin;
    h[0]->GetYaxis()->SetRangeUser(1, maxbin*1.1);
    h[0]->Draw("LE1P");

    h[1]->Scale(h[0]->Integral(-1, -1)/h[1]->Integral(-1, -1));
    h[1]->SetStats(0);
    h[1]->SetLineColor(kRed);
    h[1]->SetLineWidth(2);
    h[1]->SetMarkerStyle(20);
    h[1]->SetMarkerSize(1.5);
    h[1]->SetMarkerColor(kRed);
    h[1]->SetFillColorAlpha(kRed, 0.1);
    h[1]->Draw("HISTSAME");

    //------Calculate chi-square of data & MC histogram------
    cout << "chi-square of " << XaxisName.c_str() << " is: " << h[0]->Chi2Test(h[1], "CHI2/NDF") << endl;
    //-------------------------------------------------------

    pad1->Update();

    TLegend *legend = new TLegend(0.67, 0.77, 0.87, 0.87);//(x1,y1,x2,y2)
    legend->SetTextSize(0.05);
    legend->AddEntry(h[0], "Sep2023 BCD", "LE1P");
    legend->AddEntry(h[1], "DY_preEE", "F");
    legend->SetBorderSize(0);
    legend->Draw();
    
    TLatex* CMS = new TLatex();
    CMS->SetNDC(true);
    CMS->SetTextFont(42);
    CMS->SetTextSize(0.05);
    CMS->DrawLatex(0.13, 0.91, "#bf{CMS} #it{work-in-progress}");
    
    TLatex* lumi = new TLatex();
    lumi->SetNDC();
    lumi->SetTextFont(42);
    lumi->SetTextSize(0.04);
    // lumi->DrawLatex(0.72, 0.91, "27.0072 fb^{-1} (13.6TeV)");
    lumi->DrawLatex(0.72, 0.91, "8.17468 fb^{-1} (13.6TeV)");

    c->cd();

    //---------Create lower plot (SF)---------//
    TPad *pad2 = new TPad("pad2", "", 0, 0, 1, 0.3);
    pad2->SetGridy();
    pad2->SetRightMargin(0.05);
    pad2->SetLeftMargin(0.13);
    pad2->SetTopMargin(0.06);
    pad2->SetBottomMargin(0.4);
    pad2->Draw();
    pad2->cd();

    //---------- Ratio calculation----------//
    TH1F *hratio = (TH1F*) h[0]->Clone();
    hratio->Divide(h[1]);

    string HistName = Form("ratio_%s", IDs[IDType].c_str());

    hratio->SetTitle("");
    hratio->GetYaxis()->SetTitle("Ratio");
    hratio->GetXaxis()->SetTitle(XaxisName.c_str());

    hratio->SetMarkerColor(kBlack);
    hratio->SetMarkerSize(1);
    hratio->SetMarkerStyle(20);
    // hratio->SetLineColor(TColor::GetColor(ColorName[0].c_str()));
    hratio->SetLineWidth(2);
    hratio->GetXaxis()->SetRangeUser(xmin,xmax);
    hratio->GetXaxis()->SetTitleSize(0.16);
    hratio->GetXaxis()->SetTitleOffset(1);
    hratio->GetXaxis()->SetLabelSize(0.12);
    hratio->GetXaxis()->SetLabelOffset(0.03);
    hratio->GetYaxis()->SetRangeUser(0.6,1.4);
    hratio->GetYaxis()->SetTitleSize(0.13);
    hratio->GetYaxis()->SetTitleOffset(0.4);
    hratio->GetYaxis()->SetLabelSize(0.12);
    hratio->GetYaxis()->SetNdivisions(502);
    hratio->Draw("LE1P");

    // fout->cd();
    hratio->SetName(HistName.c_str());
    // hratio->Write();
    // fout->WriteObjectAny(&hratio, "TH1F*", "hratio");
    c->Print(outName.c_str());
    c->Close();
}

//! Change path (for different dR region)
string picture_path = "69200";
vector<TEfficiency*> getHist(int IDType, int CutType, int R9Type, int region){

    vector<TFile*> vfiles; vfiles.clear();
    //! Change input file
    vfiles.push_back(new TFile(Form("./try_NanoAODv12_preEE/%s/data_BCD.root", picture_path.c_str())));
    vfiles.push_back(new TFile(Form("./try_NanoAODv12_preEE/%s/DY_preEE.root", picture_path.c_str())));
    //------------Create Efficiency Histogram------------//
    vector<double> bin_phopt = {10, 15, 20, 25, 35, 50, 70, 200};
    vector<double> bin_phosceta = {-2.5, -2, -1.566, -1.4442, -1, -0.5, 0, 0.5, 1, 1.4442, 1.566, 2, 2.5};
    vector<double> bin_npvs = {10, 15, 20, 25, 30, 35, 40, 50, 100};
    int n_bin_phopt =  bin_phopt.size() - 1;
    int n_bin_phosceta =  bin_phosceta.size() - 1;
    int n_bin_npvs =  bin_npvs.size() - 1;
    vector<TEfficiency*> eff_phoPt; eff_phoPt.clear();
    vector<TEfficiency*> eff_phoSceta; eff_phoSceta.clear();
    vector<TEfficiency*> eff_npvs; eff_npvs.clear();
    //------------Create Distribution Histogram------------//
    vector<double> bin_phopt_h = {10, 15, 20, 25, 30, 35, 40, 50, 60, 70, 90, 120, 200};
    vector<double> bin_phophi = {-TMath::Pi(), -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2, 2.5, 3, TMath::Pi()};
    vector<double> bin_phoeta = {-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5};
    vector<double> bin_zmass = {80, 82, 84, 86, 88, 90, 92, 94, 96, 98, 100};
    vector<double> bin_mva = {-1, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1};
    int n_bin_phopt_h =  bin_phopt_h.size() - 1;
    int n_bin_phophi =  bin_phophi.size() - 1;
    int n_bin_phoeta =  bin_phoeta.size() - 1;
    int n_bin_zmass =  bin_zmass.size() - 1;
    int n_bin_mva =  bin_mva.size() - 1;
    vector<TH1F*> h_phoPt; h_phoPt.clear();
    vector<TH1F*> h_phoPhi; h_phoPhi.clear();
    vector<TH1F*> h_phoEta; h_phoEta.clear();
    vector<TH1F*> h_zmass; h_zmass.clear();
    vector<TH1F*> h_mva; h_mva.clear();
    vector<TH1F*> h_mva_EB; h_mva_EB.clear();
    vector<TH1F*> h_mva_EE; h_mva_EE.clear();
    vector<TH1F*> h_mva_90; h_mva_90.clear();
    vector<TH1F*> h_mva_90_EB; h_mva_90_EB.clear();
    vector<TH1F*> h_mva_90_EE; h_mva_90_EE.clear();
    vector<TH1F*> h_dr; h_dr.clear();
    //------------Create SummaryPlot Histogram------------//
    vector<double> bin_sp = {-TMath::Pi()-1, TMath::Pi()+1};
    int n_bin_sp =  bin_sp.size() - 1;
    vector<TEfficiency*> eff_sp; eff_sp.clear();

    vector<string> fileName = {"Prompt", "mc"};
    // cout << vfiles.size() << endl;
    for (int ifile = 0; ifile < vfiles.size(); ifile++){
        
        eff_phoPt.push_back(new TEfficiency(Form("eff_phoPt_%s", fileName[ifile].c_str()), " ", n_bin_phopt, bin_phopt.data()));
        eff_phoSceta.push_back(new TEfficiency(Form("eff_phoSceta_%s", fileName[ifile].c_str()), " ", n_bin_phosceta, bin_phosceta.data()));
        eff_npvs.push_back(new TEfficiency(Form("eff_npvs_%s", fileName[ifile].c_str()), " ", n_bin_npvs, bin_npvs.data()));

        // h_phoPt.push_back(new TH1F(Form("h_phoPt_%s", fileName[ifile].c_str()), " ", 37, 15, 200));
        h_phoPt.push_back(new TH1F(Form("h_phoPt_%s", fileName[ifile].c_str()), " ", n_bin_phopt_h, bin_phopt_h.data()));
        h_phoPhi.push_back(new TH1F(Form("h_phoPhi_%s", fileName[ifile].c_str()), " ", n_bin_phophi, bin_phophi.data()));
        h_phoEta.push_back(new TH1F(Form("h_phoEta_%s", fileName[ifile].c_str()), " ", n_bin_phoeta, bin_phoeta.data()));
        h_zmass.push_back(new TH1F(Form("h_zmass_%s", fileName[ifile].c_str()), " ", 40, 80, 100));
        h_mva.push_back(new TH1F(Form("h_mva_%s", fileName[ifile].c_str()), " ", 40, -1, 1));
        h_mva_EB.push_back(new TH1F(Form("h_mva_EB_%s", fileName[ifile].c_str()), " ", 40, -1, 1));
        h_mva_EE.push_back(new TH1F(Form("h_mva_EE_%s", fileName[ifile].c_str()), " ", 40, -1, 1));
        h_mva_90.push_back(new TH1F(Form("h_mva_90_%s", fileName[ifile].c_str()), " ", 40, -0.3, 1));
        h_mva_90_EB.push_back(new TH1F(Form("h_mva_90_EB_%s", fileName[ifile].c_str()), " ", 40, 0, 1));
        h_mva_90_EE.push_back(new TH1F(Form("h_mva_90_EE_%s", fileName[ifile].c_str()), " ", 40, -0.3, 1));
        h_dr.push_back(new TH1F(Form("h_dr_%s", fileName[ifile].c_str()), " ", 20, 0.1, 1));

        eff_sp.push_back(new TEfficiency(Form("eff_sp_%s", fileName[ifile].c_str()), " ", n_bin_sp, bin_sp.data()));

        TTree *inTree = (TTree*) vfiles[ifile]->Get("tree");
        
        Double_t Cutbased, mva, MVA_WP80, MVA_WP90;
        Double_t HasPix, CSEV;
        Double_t pho_pt, pho_phi, pho_eta, zmass;
        Double_t weights;
        Double_t pho_r9, pho_EB, pho_EE;
        Double_t PV_x, PV_y, PV_z, PV_npvs;
        Double_t mu1_phi, mu1_eta, mu2_phi, mu2_eta;
        inTree->SetBranchAddress("pho_cutBased", &Cutbased);
        inTree->SetBranchAddress("pho_mvaID", &mva);
        inTree->SetBranchAddress("pho_mvaID_WP80", &MVA_WP80);
        inTree->SetBranchAddress("pho_mvaID_WP90", &MVA_WP90);
        inTree->SetBranchAddress("pho_pixelseed", &HasPix);
        inTree->SetBranchAddress("pho_electronveto", &CSEV);
        inTree->SetBranchAddress("pho_pt", &pho_pt);
        inTree->SetBranchAddress("pho_phi", &pho_phi);
        inTree->SetBranchAddress("pho_eta", &pho_eta);
        inTree->SetBranchAddress("m1_phi", &mu1_phi);
        inTree->SetBranchAddress("m1_eta", &mu1_eta);
        inTree->SetBranchAddress("m2_phi", &mu2_phi);
        inTree->SetBranchAddress("m2_eta", &mu2_eta);
        inTree->SetBranchAddress("pho_r9", &pho_r9);
        inTree->SetBranchAddress("pho_EB", &pho_EB);
        inTree->SetBranchAddress("pho_EE", &pho_EE);
        inTree->SetBranchAddress("Z_m", &zmass);
        inTree->SetBranchAddress("weights", &weights);
        inTree->SetBranchAddress("PV_x", &PV_x);
        inTree->SetBranchAddress("PV_y", &PV_y);
        inTree->SetBranchAddress("PV_z", &PV_z);
        inTree->SetBranchAddress("PV_npvs", &PV_npvs);
        for (Long64_t ev=0; ev < inTree->GetEntriesFast(); ev++){              
            inTree->GetEntry(ev);

            bool cut;
            if (CutType == 0) cut = (HasPix == 0); //HasPix == 0, pass pixel seed veto
            if (CutType == 1) cut = (CSEV == 1);  //CSEV == 1, pass conversion safe electron veto
            // CSEV: 
            // This veto requires the absence of charged particle tracks, with a hit in the innermost layer of the 
            // pixel detector not matched to a reconstructed conversion vertex, pointing to the photon cluster in the ECAL
            // Pixel seed: 
            // A more efficient rejection of electrons can be achieved by rejecting any photon for which a pixel detector seed consisting of 
            // at least two hits in the pixel detector points to the ECAL within some window defined around the photon SC position.
            
            h_mva[ifile]->Fill(mva, weights);
            if (pho_EB > 0.5){   // Loose ID ensure include all event
                h_mva_EB[ifile]->Fill(mva, weights);
            }
            if (pho_EE > 0.5){
                h_mva_EE[ifile]->Fill(mva, weights);
            }

            Double_t dr1 = ROOT::VecOps::DeltaR(pho_eta, mu1_eta, pho_phi, mu1_phi);
            Double_t dr2 = ROOT::VecOps::DeltaR(pho_eta, mu2_eta, pho_phi, mu2_phi);
            Double_t dr = (dr1 < dr2) ? dr1 : dr2;
            h_dr[ifile]->Fill(dr, weights);

            bool passID;
            if (IDType == 0) passID = (Cutbased >= 1); //CutBased = 1: Loose
            if (IDType == 1) passID = (Cutbased >= 2); //CutBased = 2: Medium
            if (IDType == 2) passID = (Cutbased == 3); //CutBased = 3: Tight
            if (IDType == 3) passID = (MVA_WP80 == 1);
            if (IDType == 4) passID = (MVA_WP90 == 1);
            if (!passID) continue;

            h_phoPt[ifile]->Fill(pho_pt, weights);
            h_phoPhi[ifile]->Fill(pho_phi, weights);
            h_phoEta[ifile]->Fill(pho_eta, weights);
            h_zmass[ifile]->Fill(zmass, weights);

            if (IDType == 4){
                h_mva_90[ifile]->Fill(mva, weights);
            }

            bool cat;
            if (R9Type == 0 && region == 0) cat = (pho_EB > 0.5); // >0.5 means == 1
            if (R9Type == 1 && region == 0) cat = (pho_r9 > 0.96) & (pho_EB > 0.5);
            if (R9Type == 2 && region == 0) cat = (pho_r9 < 0.96) & (pho_EB > 0.5);
            if (R9Type == 0 && region == 1) cat = (pho_EE > 0.5);
            if (R9Type == 1 && region == 1) cat = (pho_r9 > 0.96) & (pho_EE > 0.5);
            if (R9Type == 2 && region == 1) cat = (pho_r9 < 0.96) & (pho_EE > 0.5);
            if (R9Type == 0 && region == 2) cat = true;
            if (R9Type == 1 && region == 2) cat = (pho_r9 > 0.96);
            if (R9Type == 2 && region == 2) cat = (pho_r9 < 0.96);
            if (!cat) continue;


            if (IDType == 4 && R9Type == 0 && region == 0){
                h_mva_90_EB[ifile]->Fill(mva, weights);       // MVA_WP90 is only bool, mva is the value of mva score
            }
            if (IDType == 4 && R9Type == 0 && region == 1){
                h_mva_90_EE[ifile]->Fill(mva, weights);
            }

            eff_phoPt[ifile]->FillWeighted(cut, weights, pho_pt);
            float pho_sceta = Photon_SCEta(pho_eta, pho_phi, pho_EB, pho_EE, PV_x, PV_y, PV_z);
            eff_phoSceta[ifile]->FillWeighted(cut, weights, pho_sceta);
            eff_npvs[ifile]->FillWeighted(cut, weights, PV_npvs);
            eff_sp[ifile]->FillWeighted(cut, weights, pho_phi);
        }
    }
    // TH1F *heff_prompt = (TH1F*) eff_phoPt[0]->GetCopyPassedHisto();
    // for (int i=0; i<eff_phoPt[0]->GetCopyTotalHisto()->GetNbinsX(); i++){
    //     cout << eff_phoPt[0]->GetCopyTotalHisto()->GetBinContent(i+1) << endl;
    // }
    
    DrawEfficiency(eff_phoPt,10,200,0,"P^{#gamma}_{t} [Gev]",IDType,CutType,R9Type,region, Form("./try_NanoAODv12_preEE/%s/sf/phoPt/%s_%s_%s_%s.pdf", picture_path.c_str(), Region[region].c_str(),R9[R9Type].c_str(), IDs[IDType].c_str(), Cuts[CutType].c_str()));
    DrawEfficiency(eff_phoSceta,-2.5,2.5,5,"#eta^{#gamma}_{SC}",IDType,CutType,R9Type,region, Form("./try_NanoAODv12_preEE/%s/sf/phoSceta/%s_%s_%s_%s.pdf", picture_path.c_str(), Region[region].c_str(),R9[R9Type].c_str(), IDs[IDType].c_str(), Cuts[CutType].c_str()));
    DrawEfficiency(eff_npvs,10,100,6,"Number of primary vertices",IDType,CutType,R9Type,region, Form("./try_NanoAODv12_preEE/%s/sf/npvs/%s_%s_%s_%s.pdf", picture_path.c_str(), Region[region].c_str(),R9[R9Type].c_str(), IDs[IDType].c_str(), Cuts[CutType].c_str()));
    DrawDistribution(h_phoPt,10,200,0,true,"P^{#gamma}_{t} [Gev]",IDType, Form("./try_NanoAODv12_preEE/%s/kinematic/phoPt_%s.pdf", picture_path.c_str(), IDs[IDType].c_str()));
    DrawDistribution(h_phoPhi,-TMath::Pi(),TMath::Pi(),1,true,"#phi^{#gamma}",IDType, Form("./try_NanoAODv12_preEE/%s/kinematic/phoPhi_%s.pdf", picture_path.c_str(), IDs[IDType].c_str()));
    DrawDistribution(h_phoEta,-5,5,2,true,"#eta^{#gamma}",IDType, Form("./try_NanoAODv12_preEE/%s/kinematic/phoEta_%s.pdf", picture_path.c_str(), IDs[IDType].c_str()));
    DrawDistribution(h_zmass,80,100,3,true,"M^{#mu#mu#gamma}",IDType, Form("./try_NanoAODv12_preEE/%s/kinematic/zmass_%s.pdf", picture_path.c_str(), IDs[IDType].c_str()));
    DrawDistribution_wRatio(h_mva,-1,1,"MVA",IDType, Form("./try_NanoAODv12_preEE/%s/kinematic/mva.pdf", picture_path.c_str()));
    DrawDistribution_wRatio(h_mva_EB,-1,1,"MVA_EB",IDType, Form("./try_NanoAODv12_preEE/%s/kinematic/mva_EB.pdf", picture_path.c_str()));
    DrawDistribution_wRatio(h_mva_EE,-1,1,"MVA_EE",IDType, Form("./try_NanoAODv12_preEE/%s/kinematic/mva_EE.pdf", picture_path.c_str()));
    if (IDType == 4){
        DrawDistribution_wRatio(h_mva_90,-0.3,1,"MVA_WP90",IDType, Form("./try_NanoAODv12_preEE/%s/kinematic/mva_WP90.pdf", picture_path.c_str()));
    }
    if (IDType == 4 && R9Type == 0 && region == 0){
        DrawDistribution_wRatio(h_mva_90_EB,0,1,"MVA_WP90_EB",IDType, Form("./try_NanoAODv12_preEE/%s/kinematic/mva_WP90_EB.pdf", picture_path.c_str()));
    }
    if (IDType == 4 && R9Type == 0 && region == 1){
        DrawDistribution_wRatio(h_mva_90_EE,-0.3,1,"MVA_WP90_EE",IDType, Form("./try_NanoAODv12_preEE/%s/kinematic/mva_WP90_EE.pdf", picture_path.c_str()));
    }
    DrawDistribution(h_dr,0.1,1,4,false,"#Delta R (#mu,#gamma)",IDType, Form("./try_NanoAODv12_preEE/%s/kinematic/dr.pdf", picture_path.c_str()));
    return eff_sp;
}
void DrawSummaryPlot(vector<vector<float>> sp_value){

    vector<vector<TH1F*>> vhsp;
    for (int iCuts = 0; iCuts < Cuts.size(); iCuts++){

        TFile *fsf_out = TFile::Open(Form("./%s_SummaryPlot_22_preEE.root", Cuts[iCuts].c_str()), "RECREATE");
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

        TLegend *legend = new TLegend(0.16, 0.7, 0.68, 0.83); //(x1,y1,x2,y2)

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
            vhsp[iCuts][iIDs]->GetYaxis()->SetRangeUser(0.83, 1.10);
            vhsp[iCuts][iIDs]->SetMarkerStyle(20);
            // if (iIDs >= 3) vhsp[iCuts][iIDs]->SetMarkerStyle(24);
            vhsp[iCuts][iIDs]->SetMarkerColor(TColor::GetColor(ColorName_sp[iIDs].c_str()));
            vhsp[iCuts][iIDs]->SetMarkerSize(1.5);
            vhsp[iCuts][iIDs]->SetLineColor(TColor::GetColor(ColorName_sp[iIDs].c_str()));
            vhsp[iCuts][iIDs]->SetLineWidth(2);
            // cout << "here" << endl;
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
        ltx->DrawLatex(0.17, 0.83, information.c_str());

        TLatex* CMS = new TLatex();
        CMS->SetNDC(true);
        CMS->SetTextFont(42);
        CMS->SetTextSize(0.04);
        CMS->DrawLatex(0.17, 0.91, "#bf{CMS} #it{work-in-progress}");

        TLatex* lumi = new TLatex();
        lumi->SetNDC();
        lumi->SetTextFont(42);
        lumi->SetTextSize(0.03);
        // lumi->DrawLatex(0.69, 0.91, "27.0072 fb^{-1} (13.6TeV)");
        lumi->DrawLatex(0.69, 0.91, "8.17468 fb^{-1} (13.6TeV)");

        canvas->Print(Form("./try_NanoAODv12_preEE/%s/sf/SummaryPlot_%s.pdf", picture_path.c_str(), Cuts[iCuts].c_str()));
        canvas->Close();

        fsf_out->Close();
    }
}

void SF_eff(){
    // fout->cd();

    vector<vector<float>> sp_value = {};
    for (int iCuts = 0; iCuts < Cuts.size(); iCuts++){
        for (int iIDs = 0; iIDs < IDs.size(); iIDs++){
            for (int iRegion = 0; iRegion < Region.size(); iRegion++){
                for (int iR9 = 0; iR9 < R9.size(); iR9++){
                    getHist(iIDs, iCuts, iR9, iRegion);
                    //----------Create value of summaryplot----------//
                    SummaryPlot(getHist(iIDs, iCuts, iR9, iRegion));
                    float sf_y =  SummaryPlot(getHist(iIDs, iCuts, iR9, iRegion))[0];
                    float sf_unc = SummaryPlot(getHist(iIDs, iCuts, iR9, iRegion))[1];
                    vector<float> sf_value = {sf_y, sf_unc};
                    sp_value.push_back(sf_value);
                }
            }
        }
    }
    DrawSummaryPlot(sp_value);

    // fout->Close();
    // delete fout;

    // getHist(1, 0, 2, 1);
}
