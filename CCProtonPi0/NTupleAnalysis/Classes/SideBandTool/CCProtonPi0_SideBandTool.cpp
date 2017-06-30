#include "CCProtonPi0_SideBandTool.h"

using namespace PlotUtils;

CCProtonPi0_SideBandTool::CCProtonPi0_SideBandTool() : CCProtonPi0_NTupleAnalysis()
{
    current_unv = 0;
    OpenRootFiles();
    initSideBands();
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
    //Plot(Michel);
    //Plot(pID);
    //Plot(LowInvMass);
    //Plot(HighInvMass);

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

    //Plot(sb, sb.muon_P, 0, "muon_P");
    //Plot(sb, sb.muon_P, 1, "muon_P");
    //Plot(sb, sb.muon_theta, 0, "muon_theta");
    //Plot(sb, sb.muon_theta, 1, "muon_theta");
    //Plot(sb, sb.pi0_P, 0, "pi0_P");
    //Plot(sb, sb.pi0_P, 1, "pi0_P");
    //Plot(sb, sb.pi0_KE, 0, "pi0_KE");
    //Plot(sb, sb.pi0_KE, 1, "pi0_KE");
    //Plot(sb, sb.pi0_theta, 0, "pi0_theta");
    //Plot(sb, sb.pi0_theta, 1, "pi0_theta");
    //Plot(sb, sb.neutrino_E, 0, "neutrino_E");
    //Plot(sb, sb.neutrino_E, 1, "neutrino_E");
    //Plot(sb, sb.QSq, 0, "QSq");
    //Plot(sb, sb.QSq, 1, "QSq");
//    Plot(sb, sb.W, 0, "W");
//    Plot(sb, sb.W, 1, "W");
}

void CCProtonPi0_SideBandTool::Plot(SideBand &sb, int ind)
{
    //Plot_NoRatio(ind, sb.name, "pi0_InvMass",
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

void CCProtonPi0_SideBandTool::Plot_NoRatio(int ind, std::string sb_name, std::string var_name, MnvH1D* data, MnvH1D* mc_total, MnvH1D* signal, MnvH1D* WithPi0, MnvH1D* QELike, MnvH1D* SinglePiPlus, MnvH1D* Other)
{
    (void) mc_total;
    std::string type;
    if (ind == 0) type = "Nominal";
    else type = "Fitted";

    std::string norm = "POT";
    std::string plot_title = "Side Band: " + sb_name + " " + type + " " + norm + " Normalized";

    // Get Histograms -- Use new Histograms not to change originals
    MnvH1D* h_data = new MnvH1D(*data);
    h_data->GetYaxis()->CenterTitle();
    h_data->GetXaxis()->SetNdivisions(5,5,0);
    h_data->GetYaxis()->SetNdivisions(5,5,0);

    MnvH1D* h_signal = new MnvH1D(*signal);
    h_signal->GetYaxis()->CenterTitle();
    h_signal->SetFillColor(kGreen+1);
    h_signal->SetLineColor(kGreen+1);
    h_signal->SetFillStyle(3002);
    h_signal->SetLineWidth(2);
    h_signal->GetXaxis()->SetNdivisions(5,5,0);
    h_signal->GetYaxis()->SetNdivisions(5,5,0);

    MnvH1D* h_WithPi0 = new MnvH1D(*WithPi0);
    h_WithPi0->GetYaxis()->CenterTitle();
    h_WithPi0->SetFillColor(kRed);
    h_WithPi0->SetLineColor(kRed);
    h_WithPi0->SetFillStyle(3002);
    h_WithPi0->SetLineWidth(2);
    h_WithPi0->GetXaxis()->SetNdivisions(5,5,0);
    h_WithPi0->GetYaxis()->SetNdivisions(5,5,0);
  
    MnvH1D* h_QELike = new MnvH1D(*QELike);
    h_QELike->GetYaxis()->CenterTitle();
    h_QELike->SetFillColor(kOrange-1);
    h_QELike->SetLineColor(kOrange-1);
    h_QELike->SetFillStyle(3002);
    h_QELike->SetLineWidth(2);
    h_QELike->GetXaxis()->SetNdivisions(5,5,0);
    h_QELike->GetYaxis()->SetNdivisions(5,5,0);
 
    MnvH1D* h_SinglePiPlus = new MnvH1D(*SinglePiPlus);
    h_SinglePiPlus->GetYaxis()->CenterTitle();
    h_SinglePiPlus->SetFillColor(kBlue);
    h_SinglePiPlus->SetLineColor(kBlue);
    h_SinglePiPlus->SetFillStyle(3002);
    h_SinglePiPlus->SetLineWidth(2);
    h_SinglePiPlus->GetXaxis()->SetNdivisions(5,5,0);
    h_SinglePiPlus->GetYaxis()->SetNdivisions(5,5,0);

    MnvH1D* h_Other = new MnvH1D(*Other);
    h_Other->GetYaxis()->CenterTitle();
    h_Other->SetFillColor(kGray+2);
    h_Other->SetLineColor(kGray+2);
    h_Other->SetFillStyle(3002);
    h_Other->SetLineWidth(2);
    h_Other->GetXaxis()->SetNdivisions(5,5,0);
    h_Other->GetYaxis()->SetNdivisions(5,5,0);
  
    // Clear Error Bars
    h_signal->ClearAllErrorBands();
    h_WithPi0->ClearAllErrorBands();
    h_QELike->ClearAllErrorBands();
    h_SinglePiPlus->ClearAllErrorBands();
    h_Other->ClearAllErrorBands();

    TObjArray* mc_hists = new TObjArray;
    mc_hists->Add(h_Other);
    mc_hists->Add(h_SinglePiPlus);
    mc_hists->Add(h_QELike);
    mc_hists->Add(h_WithPi0);
    mc_hists->Add(h_signal);
 
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(2);
    gStyle->SetStripDecimals(false);

    plotter->axis_minimum = 0.01;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    plotter->data_marker_size = 0.8;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;
    plotter->headroom = 1.75;

    plotter->mc_line_width = 2;
    
    TCanvas* c = new TCanvas("c");

    /*
       void MnvPlotter::DrawDataStackedMC(
            const MnvH1D *  dataHist,
            const TObjArray *    mcHists,
            const Double_t   mcScale = 1.0,
            const std::string &  legPos = "L",
            const std::string &  dataName = "Data",
            const Int_t  mcBaseColor = 2,
            const Int_t  mcColorOffset = 1,
            const Int_t  mcFillStyle = 3001,
            const char *     xaxislabel = "",
            const char *     yaxislabel = "",
            bool     cov_area_normalize = false   
       )    
    */
    plotter->DrawDataStackedMC(h_data, mc_hists, POT_ratio, "N", "Data (3.33e20 POT)", 0, 0, 3002);

    // Add Legend
    double leg_x_min = 0.55;
    double leg_x_max = 0.90;
    double leg_y_min = 0.55;
    double leg_y_max = 0.85;

    h_data->SetMarkerStyle(plotter->data_marker);
    h_data->SetMarkerSize(plotter->data_marker_size);
    h_data->SetMarkerColor(plotter->data_color);
    h_data->SetLineColor(plotter->data_color);
    h_data->SetLineWidth(plotter->data_line_width);

    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    legend->SetBorderSize(0);
    legend->SetFillColor(-1);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.04);
    legend->AddEntry(h_data, "Data", "lep");
    legend->AddEntry(h_signal, "Signal", "f");
    legend->AddEntry(h_WithPi0, "Bkgrd: #pi^{0} + meson(s)", "f");
    legend->AddEntry(h_QELike, "Bkgrd: zero meson", "f");
    legend->AddEntry(h_SinglePiPlus, "Bkgrd: charged meson(s)", "f");
    legend->AddEntry(h_Other, "Bkgrd: other", "f");
    legend->SetTextSize(0.04);
    legend->SetTextFont(42);
    legend->Draw();

    // Add Alines if Original Side Band
    //if (false) {
    if (sb_name.compare("Original") == 0) {
        double max_bin = h_data->GetMaximumBin();
        double hist_max = h_data->GetBinContent(max_bin)*1.2;

        TLine line;
        TArrow arrow;
        // Cut Line at 60 MeV
        line.SetLineWidth(2);
        line.SetLineStyle(2);
        line.SetLineColor(kBlack);
        arrow.SetLineWidth(2);
        arrow.SetLineColor(kBlack);
        line.DrawLine(60.0, 0.0, 60.0, hist_max);
        arrow.DrawArrow(60.0, hist_max, 50.0, hist_max, 0.01);

        line.DrawLine(200.0, 0.0, 200.0, hist_max);
        arrow.DrawArrow(200.0, hist_max, 210.0, hist_max, 0.01);

        arrow.SetLineColor(kBlue);
        arrow.DrawArrow(134.98,hist_max/8,134.98,0, 0.01); 
    } 

    // Add Normalization Label
    TLatex norm_text;
    norm_text.SetNDC();
    norm_text.SetTextColor(kBlue);
    norm_text.SetTextSize(0.03);
    norm_text.SetTextAlign(22);
    norm_text.DrawLatex(0.35,0.85,"POT Normalized");

    double info_text_x = 0.85;
    double info_text_y = 0.85;
    TLatex info_text;
    info_text.SetNDC();
    info_text.SetTextColor(kBlack);
    info_text.SetTextFont(62);
    info_text.SetTextAlign(12);
    info_text.SetTextSize(0.04);
    if ( ind == 0) info_text.DrawLatex(info_text_x, info_text_y, "(a)");
    else info_text.DrawLatex(info_text_x, info_text_y, "(b)");

    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name;
    out_name = plotDir + var_name + "_" + sb_name + "_" + type + "_" + norm + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete h_data;
    delete h_signal;
    delete h_WithPi0;
    delete h_QELike;
    delete h_SinglePiPlus;
    delete h_Other;
    delete legend;
    delete c;
    delete plotter;
}

void CCProtonPi0_SideBandTool::Plot(int ind, std::string sb_name, std::string var_name, MnvH1D* data, MnvH1D* mc_total, MnvH1D* signal, MnvH1D* WithPi0, MnvH1D* QELike, MnvH1D* SinglePiPlus, MnvH1D* Other)
{
    std::string type;
    if (ind == 0) type = "Nominal";
    else type = "Fitted";

    std::string norm = "POT";
    std::string plot_title = "Side Band: " + sb_name + " " + type + " " + norm + " Normalized";
 
    // Get Histograms -- Use new Histograms not to change originals
    ColorHists(data, signal, WithPi0, QELike, SinglePiPlus, Other);
    TH1D* h_data = new TH1D(data->GetBinNormalizedCopy().GetCVHistoWithError());
    TH1D* h_signal = new TH1D(signal->GetBinNormalizedCopy().GetCVHistoWithError());
    TH1D* h_WithPi0 = new TH1D(WithPi0->GetBinNormalizedCopy().GetCVHistoWithError());
    TH1D* h_QELike = new TH1D(QELike->GetBinNormalizedCopy().GetCVHistoWithError());
    TH1D* h_SinglePiPlus = new TH1D(SinglePiPlus->GetBinNormalizedCopy().GetCVHistoWithError());
    TH1D* h_Other = new TH1D(Other->GetBinNormalizedCopy().GetCVHistoWithError());
  
    // MC Total depend on the Modification
    //      If Raws - take the mc_total directly
    //      If Modified - Add all mc models;
    TH1D* h_mc_total;
    if (ind == 0){
        h_mc_total = new TH1D(mc_total->GetBinNormalizedCopy().GetCVHistoWithError());
    }else{
        h_mc_total = new TH1D(*h_signal);
        h_mc_total->Add(h_WithPi0);
        h_mc_total->Add(h_QELike);
        h_mc_total->Add(h_SinglePiPlus);
        h_mc_total->Add(h_Other);
    }

    // Scale Histograms
    double mc_ratio = POT_ratio; 
    h_mc_total->Scale(mc_ratio);
    h_signal->Scale(mc_ratio);
    h_WithPi0->Scale(mc_ratio);
    h_QELike->Scale(mc_ratio);
    h_SinglePiPlus->Scale(mc_ratio);
    h_Other->Scale(mc_ratio);

    // ------------------------------------------------------------------------
    // Unique Plot for Single Error Band, Single Universe
    //      Comment this section out   
    // ------------------------------------------------------------------------
    //std::string err_band = "GENIE_MaRES";
    //int hist_ind = 1;
    //TH1D* h_data = new TH1D(*(data->GetVertErrorBand(err_band)->GetHist(hist_ind)));
    //TH1D* h_signal = new TH1D(*(signal->GetVertErrorBand(err_band)->GetHist(hist_ind)));
    //TH1D* h_WithPi0 = new TH1D(*(WithPi0->GetVertErrorBand(err_band)->GetHist(hist_ind)));
    //TH1D* h_QELike = new TH1D(*(QELike->GetVertErrorBand(err_band)->GetHist(hist_ind)));
    //TH1D* h_SinglePiPlus = new TH1D(*(SinglePiPlus->GetVertErrorBand(err_band)->GetHist(hist_ind)));
    //TH1D* h_Other = new TH1D(*(Other->GetVertErrorBand(err_band)->GetHist(hist_ind)));
    //TH1D* h_mc_total;
    //if (ind == 0){
    //    h_mc_total = new TH1D(*(mc_total->GetVertErrorBand(err_band)->GetHist(hist_ind)));
    //}else{
    //    h_mc_total = new TH1D(*h_signal);
    //    h_mc_total->Add(h_WithPi0);
    //    h_mc_total->Add(h_QELike);
    //    h_mc_total->Add(h_SinglePiPlus);
    //    h_mc_total->Add(h_Other);
    //}
    //ColorHists(h_data, h_signal, h_WithPi0, h_QELike, h_SinglePiPlus, h_Other);
 
    //// Scale Histograms
    //h_data->Scale(1,"width");
    //double mc_ratio = POT_ratio; 
    //h_mc_total->Scale(mc_ratio,"width");
    //h_signal->Scale(mc_ratio,"width");
    //h_WithPi0->Scale(mc_ratio,"width");
    //h_QELike->Scale(mc_ratio,"width");
    //h_SinglePiPlus->Scale(mc_ratio,"width");
    //h_Other->Scale(mc_ratio,"width");
    // ------------------------------------------------------------------------
   

    // Create Canvas and Divide it into two
    TCanvas* c = new TCanvas("c","c",1280,1280);
    
    // Upper Pad is the Data vs MC
    TPad *pad1 = new TPad("pad1", "pad1", 0.05, 0.3, 1, 1.0); 
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

    hs->Draw("HIST");

    // Styling
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0); 

    hs->GetYaxis()->SetTitle("N(Events)");
    hs->GetYaxis()->SetTitleFont(62);
    hs->GetYaxis()->SetTitleSize(0.06);
    //hs->GetYaxis()->CenterTitle();
    //hs->GetYaxis()->SetTitleOffset(1.2);
    hs->GetYaxis()->SetLabelFont(42);
    hs->GetYaxis()->SetLabelSize(0.05);

    double hist_max = hs->GetMaximum();
    hs->SetMaximum(hist_max * 1.5);

    // Add Data Plot
    h_data->GetXaxis()->SetTitle("");
    h_data->Draw("SAME E1 X0");

    // Add Legend
    TLegend *legend = new TLegend(0.6,0.65,0.9,0.9);  
    legend->AddEntry(h_data, "Data");
    legend->AddEntry(h_signal, "Signal", "f");
    legend->AddEntry(h_Other, "Bckg: Other", "f");
    legend->AddEntry(h_SinglePiPlus, "Bckg: 1 #pi^{+}", "f");
    legend->AddEntry(h_QELike, "Bckg: QE Like", "f");
    legend->AddEntry(h_WithPi0, "Bckg: #pi^{0} + X", "f");
    legend->SetTextSize(0.04);
    legend->SetTextFont(42);
    legend->Draw();

//    // Add Pi0 InvMass Lines
//    TLine line;
//    line.SetLineWidth(3);
//    line.SetLineColor(kBlack);
//
//    TArrow arrow;
//    arrow.SetLineWidth(4);
//    arrow.SetLineStyle(1);
//    arrow.SetLineColor(1); 
//    
//    // Low Inv Mass Region
//    line.DrawLine(60.0,0,60.0,hist_max);
//    arrow.DrawArrow(60.0,hist_max,60-20,hist_max,0.01,">");     
//
//    // High Inv Mass Region
//    line.DrawLine(200.0,0,200.0,hist_max);
//    arrow.DrawArrow(200.0,hist_max,200+20,hist_max,0.01,">");     

    // Add Weights as Text to Modified Plot 
    if (ind != 0){
        int nPars = 3;
        int nPoints = 136;

        TLatex* text = new TLatex;
        text->SetTextSize(0.03);
        text->SetNDC();
        text->DrawLatex(0.55, 0.60, Form("Fit Results with %d points, %d pars", nPoints, nPars));
        text->DrawLatex(0.55, 0.57, Form("Before Fit #chi^{2} = %3.2f", ChiSq_before_fit[0]));
        text->DrawLatex(0.55, 0.54, Form("Before Fit #chi^{2}/dof = %3.2f", ChiSq_before_fit[0]/(nPoints-nPars)));
        text->DrawLatex(0.55, 0.51, Form("After Fit #chi^{2} = %3.2f", ChiSq_after_fit[0]));
        text->DrawLatex(0.55, 0.48, Form("After Fit #chi^{2}/dof = %3.2f", ChiSq_after_fit[0]/(nPoints-nPars)));
        text->DrawLatex(0.55, 0.45, Form("#color[4]{wgt(ChargedPion) = %3.2f#pm %3.2f}", wgt_SinglePiPlus[0], err_SinglePiPlus[0]));
        text->DrawLatex(0.55, 0.42, Form("wgt(QELike) = %3.2f#pm %3.2f", wgt_QELike[0], err_QELike[0]));
        text->DrawLatex(0.55, 0.39, Form("#color[2]{wgt(WithPi0) = %3.2f#pm %3.2f}", wgt_WithPi0[0], err_WithPi0[0]));
        delete text;
    }
    
    // Add Plot-ChiSq
    TLatex* text = new TLatex;
    text->SetNDC();
    text->SetTextSize(0.04);
    text->SetTextColor(kBlue);
    //double plot_chisq = calc_ChiSq(data, signal, WithPi0, QELike, SinglePiPlus, Other);
    double plot_chisq = calc_ChiSq(h_data, h_signal, h_WithPi0, h_QELike, h_SinglePiPlus, h_Other);
    double nPoints = h_data->GetNbinsX()-1;
    //text->DrawLatex(0.15, 0.85, Form("Plot #chi^{2} = %3.2f", plot_chisq));
    text->DrawLatex(0.15, 0.85, Form("Plot #chi^{2}/dof = %3.2f", plot_chisq/nPoints));
    delete text;

    // Plot Lower Plot: Data vs MC Ratio
    c->cd(); // Go back to default Canvas before creating 2nd Pad
    TPad *pad2 = new TPad("pad2", "pad2", 0.05, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->Draw();
    pad2->cd();

    // Calculate the Ratio 
    TH1D* h_data_mc_ratio = new TH1D(*(h_data));
    h_data_mc_ratio->Divide(h_mc_total); 

    // Style Ratio Plot
    h_data_mc_ratio->SetTitle("");
    h_data_mc_ratio->SetLineColor(kRed);
    h_data_mc_ratio->SetLineWidth(3);
    h_data_mc_ratio->SetFillColor(kWhite);
    h_data_mc_ratio->SetMinimum(0.5);
    h_data_mc_ratio->SetMaximum(1.5);

    // X axis ratio plot settings
    //h_data_mc_ratio->GetXaxis()->SetTitle(data->GetXaxis()->GetTitle());
    h_data_mc_ratio->GetXaxis()->SetTitle("#gamma#gamma Invariant Mass [MeV]");
    h_data_mc_ratio->GetXaxis()->SetNdivisions(408);
    h_data_mc_ratio->GetXaxis()->CenterTitle();
    h_data_mc_ratio->GetXaxis()->SetTitleFont(62);
    h_data_mc_ratio->GetXaxis()->SetTitleSize(0.18);
    //h_data_mc_ratio->GetXaxis()->SetTitleOffset(1.2);
    h_data_mc_ratio->GetXaxis()->SetLabelFont(42); // Absolute font size in pixel (precision 3)
    h_data_mc_ratio->GetXaxis()->SetLabelSize(0.12);

    // Y axis ratio plot settings
    h_data_mc_ratio->GetYaxis()->CenterTitle();
    h_data_mc_ratio->GetYaxis()->SetNdivisions(408);
    h_data_mc_ratio->GetYaxis()->SetTitle("Data/MC");
    h_data_mc_ratio->GetYaxis()->SetTitleFont(62);
    h_data_mc_ratio->GetYaxis()->SetTitleSize(0.18);
    h_data_mc_ratio->GetYaxis()->SetTitleOffset(0.35);
    h_data_mc_ratio->GetYaxis()->SetLabelFont(42); // Absolute font size in pixel (precision 3)
    h_data_mc_ratio->GetYaxis()->SetLabelSize(0.12);

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

void CCProtonPi0_SideBandTool::ColorHists(TH1D* data, TH1D* signal, TH1D* WithPi0, TH1D* QELike, TH1D* SinglePiPlus, TH1D* Other)
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

double CCProtonPi0_SideBandTool::calc_ChiSq(TH1D* data, TH1D* signal, TH1D* WithPi0, TH1D* QELike, TH1D* SinglePiPlus, TH1D* Other)
{
    //double mc_ratio = POT_ratio;
    double mc_ratio = 1.0;
    double chi_sq = 0.0;

    for (int i = 1; i <= data->GetNbinsX(); ++i){
        // Get N(Events) in Single Bin
        double err = data->GetBinError(i);
        double nData = data->GetBinContent(i);
        if (nData == 0) continue;
        double nSignal = signal->GetBinContent(i);
        double nWithPi0 = WithPi0->GetBinContent(i);
        double nQELike = QELike->GetBinContent(i);
        double nSinglePiPlus = SinglePiPlus->GetBinContent(i);
        double nOther = Other->GetBinContent(i);

        // Add All MC and scale them
        double nMC_total = (nSignal + nWithPi0 + nQELike + nSinglePiPlus + nOther) * mc_ratio;
        
        double delta = pow((nData-nMC_total),2) / pow(err,2);
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
    int first_bin;
    int last_bin;
    if (sb.name.compare("LowInvMass") == 0){
        first_bin = 1;
        last_bin = 6;
    }else if (sb.name.compare("HighInvMass") == 0){
        first_bin = 21;
        last_bin = 50;
    }else{
        first_bin = 1;
        last_bin = 50;
    }
    for (unsigned int i = 0; i < sb.data_all_universes.size(); ++i){
        sb.nData.push_back(sb.data_all_universes[i]->Integral(first_bin, last_bin));
        sb.nMC.push_back(sb.mc_total_all_universes[i]->Integral(first_bin, last_bin));
        sb.nSignal.push_back(sb.signal_all_universes[i]->Integral(first_bin, last_bin));
        sb.nWithPi0.push_back(sb.WithPi0_all_universes[i]->Integral(first_bin, last_bin));
        sb.nQELike.push_back(sb.QELike_all_universes[i]->Integral(first_bin, last_bin));
        sb.nSinglePiPlus.push_back(sb.SinglePiPlus_all_universes[i]->Integral(first_bin, last_bin));
        sb.nOther.push_back(sb.Other_all_universes[i]->Integral(first_bin, last_bin));
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
