#include "CCProtonPi0_SideBandTool.h"

using namespace PlotUtils;

CCProtonPi0_SideBandTool::CCProtonPi0_SideBandTool() : CCProtonPi0_NTupleAnalysis()
{
    //current_unv = 0;
    //OpenRootFiles();
    //initSideBands();
}

void CCProtonPi0_SideBandTool::OpenRootFiles()
{
    std::string rootDir;

    rootDir = Folder_List::rootDir_sideBand_Original_mc;
    Original.f_mc_fit = new TFile(rootDir.c_str());
    Original.f_mc_var = new TFile(rootDir.c_str());

    rootDir = Folder_List::rootDir_sideBand_Original_data;
    Original.f_data_fit = new TFile(rootDir.c_str());
    Original.f_data_var = new TFile(rootDir.c_str());
    
    rootDir = Folder_List::rootDir_sideBand_Michel_mc;
    Michel.f_mc_fit = new TFile(rootDir.c_str());
    Michel.f_mc_var = new TFile(rootDir.c_str());

    rootDir = Folder_List::rootDir_sideBand_Michel_data;
    Michel.f_data_fit = new TFile(rootDir.c_str());
    Michel.f_data_var = new TFile(rootDir.c_str());

    rootDir = Folder_List::rootDir_sideBand_pID_mc;
    pID.f_mc_fit = new TFile(rootDir.c_str());
    pID.f_mc_var = new TFile(rootDir.c_str());

    rootDir = Folder_List::rootDir_sideBand_pID_data;
    pID.f_data_fit = new TFile(rootDir.c_str());
    pID.f_data_var = new TFile(rootDir.c_str());
   
    // Low and High Invariant Mass uses Original Side Band for Fit 
    rootDir = Folder_List::rootDir_sideBand_Original_mc;
    LowInvMass.f_mc_fit = new TFile(rootDir.c_str());
    rootDir = Folder_List::rootDir_sideBand_LowInvMass_mc;
    LowInvMass.f_mc_var = new TFile(rootDir.c_str());
    
    rootDir = Folder_List::rootDir_sideBand_Original_data;
    LowInvMass.f_data_fit = new TFile(rootDir.c_str());
    rootDir = Folder_List::rootDir_sideBand_LowInvMass_data;
    LowInvMass.f_data_var = new TFile(rootDir.c_str());

    rootDir = Folder_List::rootDir_sideBand_Original_mc;
    HighInvMass.f_mc_fit = new TFile(rootDir.c_str());
    rootDir = Folder_List::rootDir_sideBand_HighInvMass_mc;
    HighInvMass.f_mc_var = new TFile(rootDir.c_str());

    rootDir = Folder_List::rootDir_sideBand_Original_data;
    HighInvMass.f_data_fit = new TFile(rootDir.c_str());
    rootDir = Folder_List::rootDir_sideBand_HighInvMass_data;
    HighInvMass.f_data_var = new TFile(rootDir.c_str());
}

void CCProtonPi0_SideBandTool::initSideBands()
{
    SetNames(Original, "Original");
    SetNames(Michel, "Michel");
    SetNames(pID, "pID");
    SetNames(LowInvMass, "LowInvMass");
    SetNames(HighInvMass, "HighInvMass");

    initSideBand_FitHistograms(Original);
    initSideBand_FitHistograms(Michel);
    initSideBand_FitHistograms(pID);
    initSideBand_FitHistograms(LowInvMass);
    initSideBand_FitHistograms(HighInvMass);

    initSideBand_AllUniverses(Original);
    initSideBand_AllUniverses(Michel);
    initSideBand_AllUniverses(pID);
    initSideBand_AllUniverses(LowInvMass);
    initSideBand_AllUniverses(HighInvMass);

    initSideBand_XSecHistograms(Original);
    initSideBand_XSecHistograms(Michel);
    initSideBand_XSecHistograms(pID);
    initSideBand_XSecHistograms(LowInvMass);
    initSideBand_XSecHistograms(HighInvMass);

    GetStatistics(Original);
    GetStatistics(Michel);
    GetStatistics(pID);
    GetStatistics(LowInvMass);
    GetStatistics(HighInvMass);

    // Calc Global ChiSq Before Fit
    for (unsigned int i = 0; i < Original.data_all_universes.size(); ++i){
        double Global_ChiSq = calc_Global_ChiSq(i);
        ChiSq_before_fit.push_back(Global_ChiSq);
    }
}

void CCProtonPi0_SideBandTool::initSideBand_FitHistograms(SideBand& sb)
{
    std::string var_name = "hCut_pi0invMass";
    std::string var;

    var = var_name + "_0";
    sb.data = GetMnvH1D(sb.f_data_fit, var);   
    sb.mc_total = GetMnvH1D(sb.f_mc_fit, var);   

    var = var_name + "_1";
    sb.signal[0] = GetMnvH1D(sb.f_mc_fit, var);   
  
    var = var_name + "_3";
    sb.WithPi0[0] = GetMnvH1D(sb.f_mc_fit, var);   

    var = var_name + "_4";
    sb.QELike[0] = GetMnvH1D(sb.f_mc_fit, var);   
   
    var = var_name + "_5";
    sb.SinglePiPlus[0] = GetMnvH1D(sb.f_mc_fit, var);   
   
    var = var_name + "_6";
    sb.Other[0] = GetMnvH1D(sb.f_mc_fit, var);   
}

void CCProtonPi0_SideBandTool::initSideBand_AllUniverses(SideBand& sb)
{
    GetAllUniverses(sb.data, sb.data_all_universes, sb.err_bands_data_all_universes, sb.hist_ind_data_all_universes);
    GetAllUniverses(sb.mc_total, sb.mc_total_all_universes, sb.err_bands_mc_total_all_universes, sb.hist_ind_mc_total_all_universes);
    GetAllUniverses(sb.signal[0], sb.signal_all_universes, sb.err_bands_signal_all_universes, sb.hist_ind_signal_all_universes);
    GetAllUniverses(sb.WithPi0[0], sb.WithPi0_all_universes, sb.err_bands_WithPi0_all_universes, sb.hist_ind_WithPi0_all_universes);
    GetAllUniverses(sb.QELike[0], sb.QELike_all_universes, sb.err_bands_QELike_all_universes, sb.hist_ind_QELike_all_universes);
    GetAllUniverses(sb.SinglePiPlus[0], sb.SinglePiPlus_all_universes, sb.err_bands_SinglePiPlus_all_universes, sb.hist_ind_SinglePiPlus_all_universes);
    GetAllUniverses(sb.Other[0], sb.Other_all_universes, sb.err_bands_Other_all_universes, sb.hist_ind_Other_all_universes);

    N_Universes = sb.data_all_universes.size();
}

void CCProtonPi0_SideBandTool::initSideBand_XSecHistograms(SideBand& sb)
{
    initSideBand_XSecHistograms(sb, sb.muon_P, "muon_P");
    initSideBand_XSecHistograms(sb, sb.muon_theta, "muon_theta");
    initSideBand_XSecHistograms(sb, sb.pi0_P, "pi0_P");
    initSideBand_XSecHistograms(sb, sb.pi0_KE, "pi0_KE");
    initSideBand_XSecHistograms(sb, sb.pi0_theta, "pi0_theta");
    initSideBand_XSecHistograms(sb, sb.neutrino_E, "neutrino_E");
    initSideBand_XSecHistograms(sb, sb.QSq, "QSq");
    initSideBand_XSecHistograms(sb, sb.W, "W");
}

void CCProtonPi0_SideBandTool::initSideBand_XSecHistograms(SideBand& sb, XSec_Var& xsec_var, std::string name )
{
    std::string var_name = "SideBand_" + name;
    std::string var;

    var = var_name + "_0";
    xsec_var.data = GetMnvH1D(sb.f_data_var, var);   
    xsec_var.mc_total = GetMnvH1D(sb.f_mc_var, var);   

    var = var_name + "_1";
    xsec_var.signal[0] = GetMnvH1D(sb.f_mc_var, var);   
  
    var = var_name + "_3";
    xsec_var.WithPi0[0] = GetMnvH1D(sb.f_mc_var, var);   

    var = var_name + "_4";
    xsec_var.QELike[0] = GetMnvH1D(sb.f_mc_var, var);   
   
    var = var_name + "_5";
    xsec_var.SinglePiPlus[0] = GetMnvH1D(sb.f_mc_var, var);   
   
    var = var_name + "_6";
    xsec_var.Other[0] = GetMnvH1D(sb.f_mc_var, var);   
}

void CCProtonPi0_SideBandTool::SetNames(SideBand &sb, std::string name)
{
    sb.name = name;
    sb.model_names[0] = "Signal";
    sb.model_names[1] = "WithPi0";
    sb.model_names[2] = "QELike";
    sb.model_names[3] = "SinglePiPlus";
    sb.model_names[4] = "Other";
    sb.model_names[5] = "Total";
}

CCProtonPi0_SideBandTool::~CCProtonPi0_SideBandTool()
{

}

void CCProtonPi0_SideBandTool::Plot()
{
    // Plot Fit Results on Invariant Mass
    Plot(Original);
    Plot(Michel);
    Plot(pID);
    Plot(LowInvMass);
    Plot(HighInvMass);

    // Plot XSec Variables in Each Side Band
   
    //std::cout<<"Plotting Data MC With Error Band"<<std::endl;
    //DrawDataMCWithErrorBand(Original);
    //DrawDataMCWithErrorBand(Michel);
    //DrawDataMCWithErrorBand(pID);
    //DrawDataMCWithErrorBand(LowInvMass);
    //DrawDataMCWithErrorBand(HighInvMass);

}

void CCProtonPi0_SideBandTool::Plot(SideBand &sb)
{
    std::cout<<"Plotting "<<sb.name<<std::endl;

    // Fit Results
    Plot(sb, 0); // Raw
    Plot(sb, 1); // Modified
//
//    Plot(sb, sb.muon_P, 0, "muon_P");
//    Plot(sb, sb.muon_P, 1, "muon_P");
//    Plot(sb, sb.muon_theta, 0, "muon_theta");
//    Plot(sb, sb.muon_theta, 1, "muon_theta");
//    Plot(sb, sb.pi0_P, 0, "pi0_P");
//    Plot(sb, sb.pi0_P, 1, "pi0_P");
//    Plot(sb, sb.pi0_KE, 0, "pi0_KE");
//    Plot(sb, sb.pi0_KE, 1, "pi0_KE");
//    Plot(sb, sb.pi0_theta, 0, "pi0_theta");
//    Plot(sb, sb.pi0_theta, 1, "pi0_theta");
//    Plot(sb, sb.neutrino_E, 0, "neutrino_E");
//    Plot(sb, sb.neutrino_E, 1, "neutrino_E");
//    Plot(sb, sb.QSq, 0, "QSq");
//    Plot(sb, sb.QSq, 1, "QSq");
//    Plot(sb, sb.W, 0, "W");
//    Plot(sb, sb.W, 1, "W");
}

void CCProtonPi0_SideBandTool::Plot(SideBand &sb, int ind)
{
    Plot(ind, sb.name, "pi0_InvMass",
                sb.data, 
                sb.mc_total, 
                sb.signal[ind], 
                sb.WithPi0[ind], 
                sb.QELike[ind], 
                sb.SinglePiPlus[ind], 
                sb.Other[ind]);
}

void CCProtonPi0_SideBandTool::Plot(SideBand &sb, XSec_Var &xsec_var, int ind, std::string var_name)
{
    Plot(ind, sb.name, var_name,
                xsec_var.data, 
                xsec_var.mc_total, 
                xsec_var.signal[ind], 
                xsec_var.WithPi0[ind], 
                xsec_var.QELike[ind], 
                xsec_var.SinglePiPlus[ind], 
                xsec_var.Other[ind]);
}

void CCProtonPi0_SideBandTool::Plot(int ind, std::string sb_name, std::string var_name, MnvH1D* data, MnvH1D* mc_total, MnvH1D* signal, MnvH1D* WithPi0, MnvH1D* QELike, MnvH1D* SinglePiPlus, MnvH1D* Other)
{
    std::string type;
    if (ind == 0) type = "Raw";
    else type = "Modified";

    std::string norm = "POT";
    std::string plot_title = "Side Band: " + sb_name + " " + type + " " + norm + " Normalized";
 
    // Get Histograms -- Use new Histograms not to change originals
    ColorHists(data, signal, WithPi0, QELike, SinglePiPlus, Other);
    TH1D* h_data = new TH1D(data->GetCVHistoWithStatError());
    TH1D* h_signal = new TH1D(signal->GetCVHistoWithStatError());
    TH1D* h_WithPi0 = new TH1D(WithPi0->GetCVHistoWithStatError());
    TH1D* h_QELike = new TH1D(QELike->GetCVHistoWithStatError());
    TH1D* h_SinglePiPlus = new TH1D(SinglePiPlus->GetCVHistoWithStatError());
    TH1D* h_Other = new TH1D(Other->GetCVHistoWithStatError());
    // MC Total depend on the Modification
    //      If Raws - take the mc_total directly
    //      If Modified - Add all mc models;
    TH1D* h_mc_total;
    if (ind == 0){
        h_mc_total = new TH1D(mc_total->GetCVHistoWithStatError());
    }else{
        h_mc_total = new TH1D(*h_signal);
        h_mc_total->Add(h_WithPi0);
        h_mc_total->Add(h_QELike);
        h_mc_total->Add(h_SinglePiPlus);
        h_mc_total->Add(h_Other);
    }

    // Scale Histograms
    h_data->Scale(1,"width");
    double mc_ratio = POT_ratio; 
    h_mc_total->Scale(mc_ratio,"width");
    h_signal->Scale(mc_ratio,"width");
    h_WithPi0->Scale(mc_ratio,"width");
    h_QELike->Scale(mc_ratio,"width");
    h_SinglePiPlus->Scale(mc_ratio,"width");
    h_Other->Scale(mc_ratio,"width");

    // Create Canvas and Divide it into two
    TCanvas* c = new TCanvas("c","c",1280,1280);
    
    // Upper Pad is the Data vs MC
    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0); 
    pad1->SetBottomMargin(0); // Top and Bottom Plots attached
    pad1->Draw();               
    pad1->cd(); // pad1 is the current pad
   
    // Plot MC Models as THStack
    THStack* hs = new THStack("hs",plot_title.c_str());
    hs->Add(h_WithPi0);  
    hs->Add(h_QELike);  
    hs->Add(h_SinglePiPlus);  
    hs->Add(h_Other);  
    hs->Add(h_signal);  

    // Add MC Stacked 
    hs->Draw("HIST");
    hs->GetYaxis()->SetTitle("N(Events)");
    hs->GetXaxis()->SetTitle("");
    h_data->GetXaxis()->SetTitle("");

    // Add Data Plot
    h_data->Draw("SAME E1 X0");

    // Add Legend
    TLegend *legend = new TLegend(0.7,0.68,0.9,0.9);  
    legend->AddEntry(h_data, "Data");
    legend->AddEntry(h_signal, "Signal", "f");
    legend->AddEntry(h_Other, "Bckg: Other", "f");
    legend->AddEntry(h_SinglePiPlus, "Bckg: 1 #pi^{+}", "f");
    legend->AddEntry(h_QELike, "Bckg: QE Like", "f");
    legend->AddEntry(h_WithPi0, "Bckg: #pi^{0} + X", "f");
    legend->SetTextSize(0.03);
    legend->Draw();

    // Add Pi0 InvMass Lines
    if (var_name.compare("pi0_InvMass") == 0){
        double hist_max = h_data->GetMaximum();
        hs->SetMaximum(hist_max * 1.2);
        TLine pi0Mass;
        pi0Mass.SetLineWidth(2);
        pi0Mass.SetLineColor(kBlue);
        pi0Mass.DrawLine(134.98,0,134.98,hist_max);

        TLine pi0Mass_min;
        pi0Mass_min.SetLineWidth(2);
        pi0Mass_min.SetLineColor(kBlack);
        pi0Mass_min.DrawLine(60.0,0,60.0,hist_max);

        TLine pi0Mass_max;
        pi0Mass_max.SetLineWidth(2);
        pi0Mass_max.SetLineColor(kBlack);
        pi0Mass_max.DrawLine(200.0,0,200.0,hist_max);
    }

    // Add Weights as Text to Modified Plot 
    if (ind != 0){
        int nPars = 3;
        int nPoints = 136;

        TLatex* text = new TLatex;
        text->SetTextSize(0.03);
        text->SetNDC();
        text->DrawLatex(0.6, 0.64, Form("Fit Results with %d points, %d pars", nPoints, nPars));
        text->DrawLatex(0.6, 0.61, Form("Before Fit #chi^{2} = %3.2f", ChiSq_before_fit[0]));
        text->DrawLatex(0.6, 0.58, Form("Before Fit #chi^{2}/dof = %3.2f", ChiSq_before_fit[0]/(nPoints-nPars)));
        text->DrawLatex(0.6, 0.55, Form("After Fit #chi^{2} = %3.2f", ChiSq_after_fit[0]));
        text->DrawLatex(0.6, 0.52, Form("After Fit #chi^{2}/dof = %3.2f", ChiSq_after_fit[0]/(nPoints-nPars)));
        text->DrawLatex(0.6, 0.49, Form("#color[4]{wgt(SinglePiPlus) = %3.2f#pm %3.2f}", wgt_SinglePiPlus[0], err_SinglePiPlus[0]));
        text->DrawLatex(0.6, 0.46, Form("wgt(QELike) = %3.2f#pm %3.2f", wgt_QELike[0], err_QELike[0]));
        text->DrawLatex(0.6, 0.43, Form("#color[2]{wgt(WithPi0) = %3.2f#pm %3.2f}", wgt_WithPi0[0], err_WithPi0[0]));
        delete text;
    }
    
    // Add Plot-Area and Plot-ChiSq
    double area_data = h_data->Integral("width");
    double area_mc = h_mc_total->Integral("width");
    TLatex* areaText = new TLatex;
    areaText->SetNDC();
    areaText->SetTextSize(0.03);
    areaText->SetTextColor(kBlue);
    areaText->DrawLatex(0.15, 0.87,Form("Area(Data)/Area(MC) = %3.2f",area_data/area_mc));
    double plot_chisq = calc_ChiSq(data, signal, WithPi0, QELike, SinglePiPlus, Other);
    double nPoints = h_data->GetNbinsX();
    areaText->DrawLatex(0.15, 0.83, Form("Plot #chi^{2} = %3.2f", plot_chisq));
    areaText->DrawLatex(0.15, 0.79, Form("Plot #chi^{2}/dof = %3.2f", plot_chisq/nPoints));
    delete areaText;

    // Plot Lower Plot: Data vs MC Ratio
    c->cd(); // Go back to default Canvas before creating 2nd Pad
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->Draw();
    pad2->cd();

    // Calculate the Ratio 
    TH1D* h_data_mc_ratio = new TH1D(*(h_data));
    h_data_mc_ratio->Divide(h_mc_total); 

    // Style Ratio Plot
    h_data_mc_ratio->SetTitle("");
    h_data_mc_ratio->GetXaxis()->SetTitle(h_mc_total->GetXaxis()->GetTitle());
    h_data_mc_ratio->GetYaxis()->SetTitle("N(Data)/N(MC)");
    h_data_mc_ratio->SetLineColor(kRed);
    h_data_mc_ratio->SetLineWidth(3);
    h_data_mc_ratio->SetFillColor(kWhite);
    h_data_mc_ratio->SetMinimum(0.5);
    h_data_mc_ratio->SetMaximum(1.5);
    h_data_mc_ratio->SetStats(0);

    // Y axis ratio plot settings
    h_data_mc_ratio->GetYaxis()->SetNdivisions(505);
    h_data_mc_ratio->GetYaxis()->SetTitleSize(30);
    h_data_mc_ratio->GetYaxis()->SetTitleFont(43);
    h_data_mc_ratio->GetYaxis()->SetTitleOffset(1.55);
    h_data_mc_ratio->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h_data_mc_ratio->GetYaxis()->SetLabelSize(30);

    // X axis ratio plot settings
    h_data_mc_ratio->GetXaxis()->SetTitleSize(30);
    h_data_mc_ratio->GetXaxis()->SetTitleFont(43);
    h_data_mc_ratio->GetXaxis()->SetTitleOffset(4.);
    h_data_mc_ratio->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
    h_data_mc_ratio->GetXaxis()->SetLabelSize(30);

    // Add Ratio Plot
    h_data_mc_ratio->Draw("HIST");

    // Add Ratio = 1 Line 
    TLine ratio_1;
    ratio_1.SetLineWidth(2);
    ratio_1.SetLineStyle(7);
    ratio_1.SetLineColor(kBlue);
    double line_min = h_data->GetBinLowEdge(1);
    double line_max = h_data->GetBinLowEdge(h_data->GetNbinsX()+1);
    ratio_1.DrawLine(line_min,1,line_max,1);


    // Plot Output
    gStyle->SetOptStat(0); 
    c->Update();
    std::string plotDir = Folder_List::plotDir_SideBand;
    std::string out_name;
    out_name = plotDir + var_name + "_" + sb_name + "_" + type + "_" + norm + ".png"; 

    c->Print(out_name.c_str(),"png");

    delete legend;
    delete hs;
    delete pad1;
    delete pad2;
    delete c;
}

void CCProtonPi0_SideBandTool::ColorHists(MnvH1D* data, MnvH1D* signal, MnvH1D* WithPi0, MnvH1D* QELike, MnvH1D* SinglePiPlus, MnvH1D* Other)
{
    // MC
    signal->SetFillColor(kGreen);
    signal->SetFillStyle(3001);

    WithPi0->SetFillColor(kRed);
    WithPi0->SetFillStyle(3001);

    QELike->SetFillColor(kOrange);
    QELike->SetFillStyle(3001);

    SinglePiPlus->SetFillColor(kBlue);
    SinglePiPlus->SetFillStyle(3001);

    Other->SetFillColor(kGray);
    Other->SetFillStyle(3001);

    // Data
    data->SetMarkerColor(kBlack);
    data->SetMarkerStyle(20);
    data->SetMarkerSize(1);
    data->SetLineWidth(1);
    data->SetLineColor(kBlack);
    data->SetFillStyle(0);
}

void CCProtonPi0_SideBandTool::SaveFitResults(double chisq, double par_values[3], double par_errors[3])
{
    ChiSq_after_fit.push_back(chisq);
    wgt_WithPi0.push_back(par_values[0]);
    wgt_QELike.push_back(par_values[1]);
    wgt_SinglePiPlus.push_back(par_values[2]);

    err_WithPi0.push_back(par_errors[0]);
    err_QELike.push_back(par_errors[1]);
    err_SinglePiPlus.push_back(par_errors[2]);
}

void CCProtonPi0_SideBandTool::ApplyFitResults()
{
    ApplyFitResults(Original);
    ApplyFitResults(Michel);
    ApplyFitResults(pID);
    ApplyFitResults(LowInvMass);
    ApplyFitResults(HighInvMass);
}

void CCProtonPi0_SideBandTool::ApplyFitResults(SideBand &sb)
{
    std::cout<<"Applying Fit Result to "<<sb.name<<std::endl;
    // Clone Original Histograms
    sb.signal[1] = new MnvH1D (*sb.signal[0]);
    sb.WithPi0[1] = new MnvH1D (*sb.WithPi0[0]);
    sb.QELike[1] = new MnvH1D (*sb.QELike[0]);
    sb.SinglePiPlus[1] = new MnvH1D (*sb.SinglePiPlus[0]);
    sb.Other[1] = new MnvH1D (*sb.Other[0]);

    // Scale 
    sb.WithPi0[1]->Scale(wgt_WithPi0[0]);
    sb.QELike[1]->Scale(wgt_QELike[0]);
    sb.SinglePiPlus[1]->Scale(wgt_SinglePiPlus[0]);

    // Cross Section Variables
    ApplyFitResults(sb.muon_P);
    ApplyFitResults(sb.muon_theta);
    ApplyFitResults(sb.pi0_P);
    ApplyFitResults(sb.pi0_KE);
    ApplyFitResults(sb.pi0_theta);
    ApplyFitResults(sb.neutrino_E);
    ApplyFitResults(sb.QSq);
    ApplyFitResults(sb.W);
}

void CCProtonPi0_SideBandTool::ApplyFitResults(XSec_Var &xsec_var)
{
    // Clone Original Histograms
    xsec_var.signal[1] = new MnvH1D (*xsec_var.signal[0]);
    xsec_var.WithPi0[1] = new MnvH1D (*xsec_var.WithPi0[0]);
    xsec_var.QELike[1] = new MnvH1D (*xsec_var.QELike[0]);
    xsec_var.SinglePiPlus[1] = new MnvH1D (*xsec_var.SinglePiPlus[0]);
    xsec_var.Other[1] = new MnvH1D (*xsec_var.Other[0]);

    // Scale 
    xsec_var.WithPi0[1]->Scale(wgt_WithPi0[0]);
    xsec_var.QELike[1]->Scale(wgt_QELike[0]);
    xsec_var.SinglePiPlus[1]->Scale(wgt_SinglePiPlus[0]);
}

double CCProtonPi0_SideBandTool::calc_ChiSq(MnvH1D* data, MnvH1D* signal, MnvH1D* WithPi0, MnvH1D* QELike, MnvH1D* SinglePiPlus, MnvH1D* Other)
{
    double mc_ratio = POT_ratio;
    double chi_sq = 0.0;

    for (int i = 1; i <= data->GetNbinsX(); ++i){
        // Get N(Events) in Single Bin
        double nData = data->GetBinContent(i);
        if (nData == 0) continue;
        double nSignal = signal->GetBinContent(i);
        double nWithPi0 = WithPi0->GetBinContent(i);
        double nQELike = QELike->GetBinContent(i);
        double nSinglePiPlus = SinglePiPlus->GetBinContent(i);
        double nOther = Other->GetBinContent(i);

        // Add All MC and scale them
        double nMC_total = (nSignal + nWithPi0 + nQELike + nSinglePiPlus + nOther) * mc_ratio;
        
        double delta = pow((nData-nMC_total),2) / nData;
        chi_sq += delta;
    }
    return chi_sq;
}

void CCProtonPi0_SideBandTool::DrawDataMCWithErrorBand(SideBand &sb)
{
    // Get Histograms -- Use new Histograms not to change originals
    MnvH1D* h_data = new MnvH1D(*(sb.data));
    MnvH1D* h_mc_total = new MnvH1D(*(sb.mc_total));

    double mc_ratio = POT_ratio; 

    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas("c","c",1280,1280);

    //plotter->DrawErrorSummary(h_mc_total);
    plotter->DrawDataMCWithErrorBand(h_data, h_mc_total, mc_ratio, "TR", false, NULL, NULL,false,true);

    MnvH1D* h_signal = new MnvH1D(*(sb.signal[1]));
    MnvH1D* h_WithPi0 = new MnvH1D(*(sb.WithPi0[1]));
    MnvH1D* h_QELike = new MnvH1D(*(sb.QELike[1]));
    MnvH1D* h_SinglePiPlus = new MnvH1D(*(sb.SinglePiPlus[1]));
    MnvH1D* h_Other = new MnvH1D(*(sb.Other[1]));
    MnvH1D* h_mc_total2 = new MnvH1D(*h_signal);
    h_mc_total2->Add(h_WithPi0);
    h_mc_total2->Add(h_QELike);
    h_mc_total2->Add(h_SinglePiPlus);
    h_mc_total2->Add(h_Other);
    
    double norm_bin_width = h_mc_total2->GetNormBinWidth();
    h_mc_total2->Scale(mc_ratio);
    h_mc_total2->Scale(norm_bin_width,"width");
    h_mc_total2->SetLineWidth(3);
    h_mc_total2->SetLineColor(kBlue);
    h_mc_total2->SetFillColor(kWhite);

    h_mc_total2->Draw("HIST SAME");

    // Plot Output
    std::string plotDir = Folder_List::plotDir_SideBand;
    std::string out_name;
    out_name = plotDir + sb.name + "_MC_Errors.png"; 

    c->Print(out_name.c_str(),"png");

    delete h_data;
    delete h_mc_total;
    delete c;
    delete plotter;
}

double CCProtonPi0_SideBandTool::calc_ChiSq_SideBand(SideBand &sb, int unv, bool isPartial, int min_bin, int max_bin)
{
    if (!isPartial){
        min_bin = 1;
        max_bin = sb.data->GetNbinsX();
    }
    
    if (min_bin == max_bin){
        std::cout<<"Wrong Range for Fit"<<std::endl;
        exit(EXIT_FAILURE);
    }

    double ChiSq = 0.0;

    for (int i = 1; i <= max_bin; ++i) {
        double nData = sb.data_all_universes[unv]->GetBinContent(i);
        if (nData == 0) continue;

        // Do not use Signal and Other in Fit
        double nSignal = sb.signal_all_universes[unv]->GetBinContent(i) * POT_ratio;
        double nWithPi0 = sb.WithPi0_all_universes[unv]->GetBinContent(i) * POT_ratio;
        double nQELike = sb.QELike_all_universes[unv]->GetBinContent(i) * POT_ratio;
        double nSinglePiPlus = sb.SinglePiPlus_all_universes[unv]->GetBinContent(i) * POT_ratio;
        double nOther = sb.Other_all_universes[unv]->GetBinContent(i) * POT_ratio;
        
        double nTotalMC = nSignal + nWithPi0 + nQELike + nSinglePiPlus + nOther;

        double delta  = std::pow((nData - nTotalMC),2)/nData;
        ChiSq += delta;
    }

    return ChiSq;
}

double CCProtonPi0_SideBandTool::calc_Global_ChiSq(int unv)
{
    double ChiSq = 0;

    // Calculate ChiSq for Michel for ALL Bins
    ChiSq += calc_ChiSq_SideBand(Michel, unv);
   
    // Calculate ChiSq for pID for ALL Bins
    ChiSq += calc_ChiSq_SideBand(pID, unv);

    // Calculate ChiSq for Low Inv Mass 
    //      Inv Mass itself for first 6 bins
    ChiSq += calc_ChiSq_SideBand(LowInvMass, unv, true, 1, 6);
    
    // Calculate ChiSq for High Inv Mass
    //      Inv Mass itself for last 30 bins
    ChiSq += calc_ChiSq_SideBand(HighInvMass, unv, true, 21, 50);

    return ChiSq;
}

void CCProtonPi0_SideBandTool::GetStatistics(SideBand &sb)
{
    for (unsigned int i = 0; i < sb.data_all_universes.size(); ++i){
        sb.nData.push_back(sb.data_all_universes[i]->Integral());
        sb.nMC.push_back(sb.mc_total_all_universes[i]->Integral());
        sb.nSignal.push_back(sb.signal_all_universes[i]->Integral());
        sb.nWithPi0.push_back(sb.WithPi0_all_universes[i]->Integral());
        sb.nQELike.push_back(sb.QELike_all_universes[i]->Integral());
        sb.nSinglePiPlus.push_back(sb.SinglePiPlus_all_universes[i]->Integral());
        sb.nOther.push_back(sb.Other_all_universes[i]->Integral());
    }
}

void CCProtonPi0_SideBandTool::WriteStatistics()
{
    WriteStatistics(Original);
    WriteStatistics(Michel);
    WriteStatistics(pID);
    WriteStatistics(LowInvMass);
    WriteStatistics(HighInvMass);
}

void CCProtonPi0_SideBandTool::WriteStatistics(SideBand &sb)
{
    // Open Text File
    std::string file_name = Folder_List::output + Folder_List::textOut + "SideBand_Statistics_" + sb.name + ".txt";
    ofstream file;
    file.open(file_name.c_str());

    // Write Header
    file<<std::left;
    file.width(12); file<<"Universe"<<" "; 
    file.width(12); file<<"N(Data)"<<" ";    
    file.width(12); file<<"N(MC)"<<" ";    
    file.width(12); file<<"N(Signal)"<<" ";    
    file.width(20); file<<"N(SinglePiPlus)"<<" ";    
    file.width(12); file<<"N(QELike)"<<" ";    
    file.width(12); file<<"N(WithPi0)"<<" ";    
    file.width(12); file<<"N(Other)"<<" ";    
    file<<std::endl;

    for (unsigned int i = 0; i < sb.nData.size(); ++i){
        file.width(12); file<<i<<" "; 
        file.width(12); file<<sb.nData[i]<<" ";    
        file.width(12); file<<sb.nMC[i]<<" ";    
        file.width(12); file<<sb.nSignal[i]<<" ";    
        file.width(20); file<<sb.nSinglePiPlus[i]<<" ";    
        file.width(12); file<<sb.nQELike[i]<<" ";    
        file.width(12); file<<sb.nWithPi0[i]<<" ";    
        file.width(12); file<<sb.nOther[i]<<" ";    
        file<<std::endl;
    }

    std::cout<<"Writing "<<file_name<<std::endl;
    file.close();
}

void CCProtonPi0_SideBandTool::WriteFitResults()
{
    std::cout<<"Writing List of Weights for all Universes"<<std::endl;
    // Open Text File
    std::string file_name = Folder_List::BckgConstraints; 
    ofstream file;
    file.open(file_name.c_str());

    // Write Header
    file<<std::left;
    file.width(32); file<<"Error Band"<<" "; 
    file.width(6); file<<"Hist"<<" "; 
    file.width(20); file<<"ChiSq Before Fit"<<" ";    
    file.width(20); file<<"ChiSq After Fit"<<" ";    
    file.width(20); file<<"wgt(SinglePiPlus)"<<" ";    
    file.width(20); file<<"wgt(QELike)"<<" ";    
    file.width(20); file<<"wgt(WithPi0)"<<" ";    
    file.width(20); file<<"err(SinglePiPlus)"<<" ";    
    file.width(20); file<<"err(QELike)"<<" ";    
    file.width(20); file<<"err(WithPi0)"<<" ";    
    file<<std::endl;

    for (unsigned int i = 0; i < ChiSq_after_fit.size(); ++i){
        file.width(32); file<<Original.err_bands_data_all_universes[i]<<" "; 
        file.width(6); file<<Original.hist_ind_data_all_universes[i]<<" "; 
        file.width(20); file<<ChiSq_before_fit[i]<<" ";    
        file.width(20); file<<ChiSq_after_fit[i]<<" ";    
        file.width(20); file<<wgt_SinglePiPlus[i]<<" ";    
        file.width(20); file<<wgt_QELike[i]<<" ";    
        file.width(20); file<<wgt_WithPi0[i]<<" ";    
        file.width(20); file<<err_SinglePiPlus[i]<<" ";    
        file.width(20); file<<err_QELike[i]<<" ";    
        file.width(20); file<<err_WithPi0[i]<<" ";    
        file<<std::endl;
    }

    std::cout<<"Writing "<<file_name<<std::endl;
    file.close();
}
