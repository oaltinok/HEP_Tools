#ifndef CCProtonPi0_Plotter_Paper_cpp
#define CCProtonPi0_Plotter_Paper_cpp

#include "CCProtonPi0_Plotter.h"

using namespace PlotUtils;

void CCProtonPi0_Plotter::DrawPaper_Error_pi0_KE()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_Original_Data;
    TFile* f_data = new TFile(root_dir_data.c_str());

    std::string var_name = "pi0_KE_xsec";

    MnvH1D* hist = GetMnvH1D(f_data, var_name);
 
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
    
    gStyle->SetStripDecimals(false);
    gStyle->SetEndErrorSize(4);
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); 
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->data_line_width = 2;

    plotter->axis_minimum = 0.01;
    plotter->axis_maximum = 0.4;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.056;
    plotter->axis_title_offset_y = 0.8;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;
 
    ApplyStyle_Errors(plotter, true);

    TCanvas* c = new TCanvas("c");
    //bool MnvPlotter::DrawErrorSummary   (   MnvH1D *    h,
    //        const std::string &     legPos = "TR",
    //        const bool  includeStat = true,
    //        const bool  solidLinesOnly = true,
    //        const double    ignoreThreshold = 0.00001,
    //        const bool  covAreaNormalize = false,
    //        const std::string &     errorGroupName = "",
    //        const bool  asfrac = true,
    //        const std::string &     Ytitle = "",
    //        bool    ignoreUngrouped = false  
    //        )   
    plotter->DrawErrorSummary(hist,"N", true, true, 0.0);

    // ------------------------------------------------------------------------
    // Legend 
    // ------------------------------------------------------------------------
    TH1D* h_total = GetBinNormalizedTH1D(hist);
    h_total->SetLineWidth(3);
    h_total->SetLineColor(kBlack);
    h_total->SetLineStyle(kSolid);

    TH1D* h_stat = GetBinNormalizedTH1D(hist);
    h_stat->SetLineWidth(3);
    h_stat->SetLineColor(kBlack);
    h_stat->SetLineStyle(kDashed);

    TH1D* h_detector = GetBinNormalizedTH1D(hist);
    h_detector->SetLineWidth(3); 
    //h_detector->SetLineStyle(3);
    //h_detector->SetLineColor(kBlue+2);
    h_detector->SetLineColor(kBlue);
 
    TH1D* h_xsec = GetBinNormalizedTH1D(hist);
    h_xsec->SetLineWidth(3); 
    //h_xsec->SetLineStyle(5);       
    //h_xsec->SetLineColor(kRed+2);
    h_xsec->SetLineColor(kRed);
  
    TH1D* h_fsi = GetBinNormalizedTH1D(hist);
    h_fsi->SetLineWidth(3); 
    h_fsi->SetLineColor(kMagenta+2);
    h_fsi->SetLineStyle(9);       
   
    TH1D* h_flux = GetBinNormalizedTH1D(hist);
    h_flux->SetLineWidth(3); 
    //h_flux->SetLineStyle(6);       
    h_flux->SetLineStyle(7);       
    //h_flux->SetLineColor(TColor::GetColor("#8b4513"));
    h_flux->SetLineColor(kGreen+3);
    
    TH1D* h_other = GetBinNormalizedTH1D(hist);
    h_other->SetLineWidth(3); 
    //h_other->SetLineStyle(7);       
    h_other->SetLineStyle(4);       
    //h_other->SetLineColor(kGreen+3);
    h_other->SetLineColor(kOrange+2);
    
    TLegend *legend = new TLegend(0.25,0.65,0.75,0.85);  
    ApplyStyle_Legend(legend);
    legend->SetNColumns(2);
    legend->SetTextSize(0.05);
    legend->AddEntry(h_total, "Total Error", "l" );
    legend->AddEntry(h_stat, "Statistical", "l" );
    //legend->SetColumnSeparation(0.0);
    legend->SetEntrySeparation(0.12);
    legend->AddEntry(h_xsec, "X-Sec Model", "l");
    legend->AddEntry(h_detector, "Detector", "l");
    legend->AddEntry(h_fsi, "FSI Model", "l");
    legend->AddEntry(h_flux, "Flux", "l");
    legend->AddEntry(h_other, "Other", "l");
    legend->Draw();

    double info_text_x = 0.85;
    double info_text_y = 0.85;
    TLatex info_text;
    info_text.SetNDC();
    info_text.SetTextColor(kBlack);
    info_text.SetTextFont(62);
    info_text.SetTextAlign(12);
    info_text.SetTextSize(0.04);
    info_text.DrawLatex(info_text_x, info_text_y, "(b)");


    // Print Plot
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + "Errors_" + var_name + ".pdf";
    c->Print(out_name.c_str(), "pdf");

    delete h_total;
    delete h_stat;
    delete h_detector;
    delete h_xsec;
    delete h_fsi;
    delete h_flux;
    delete h_other;

    delete c;
    delete plotter;
}

void CCProtonPi0_Plotter::DrawPaper_Error_muon_P()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_Original_Data;
    TFile* f_data = new TFile(root_dir_data.c_str());

    std::string var_name = "muon_P_xsec";

    MnvH1D* hist = GetMnvH1D(f_data, var_name);
 
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
    
    gStyle->SetStripDecimals(false);
    gStyle->SetEndErrorSize(4);
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); 
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->data_line_width = 2;

    plotter->axis_minimum = 0.01;
    plotter->axis_maximum = 0.4;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.056;
    plotter->axis_title_offset_y = 0.8;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;
 
    ApplyStyle_Errors(plotter, true);

    TCanvas* c = new TCanvas("c");
    //bool MnvPlotter::DrawErrorSummary   (   MnvH1D *    h,
    //        const std::string &     legPos = "TR",
    //        const bool  includeStat = true,
    //        const bool  solidLinesOnly = true,
    //        const double    ignoreThreshold = 0.00001,
    //        const bool  covAreaNormalize = false,
    //        const std::string &     errorGroupName = "",
    //        const bool  asfrac = true,
    //        const std::string &     Ytitle = "",
    //        bool    ignoreUngrouped = false  
    //        )   
    plotter->DrawErrorSummary(hist,"N", true, true, 0.0);

    // ------------------------------------------------------------------------
    // Legend 
    // ------------------------------------------------------------------------
    TH1D* h_total = GetBinNormalizedTH1D(hist);
    h_total->SetLineWidth(3);
    h_total->SetLineColor(kBlack);
    h_total->SetLineStyle(kSolid);

    TH1D* h_stat = GetBinNormalizedTH1D(hist);
    h_stat->SetLineWidth(3);
    h_stat->SetLineColor(kBlack);
    h_stat->SetLineStyle(kDashed);

    TH1D* h_detector = GetBinNormalizedTH1D(hist);
    h_detector->SetLineWidth(3); 
    //h_detector->SetLineStyle(3);
    //h_detector->SetLineColor(kBlue+2);
    h_detector->SetLineColor(kBlue);
 
    TH1D* h_xsec = GetBinNormalizedTH1D(hist);
    h_xsec->SetLineWidth(3); 
    //h_xsec->SetLineStyle(5);       
    //h_xsec->SetLineColor(kRed+2);
    h_xsec->SetLineColor(kRed);
  
    TH1D* h_fsi = GetBinNormalizedTH1D(hist);
    h_fsi->SetLineWidth(3); 
    h_fsi->SetLineColor(kMagenta+2);
    h_fsi->SetLineStyle(9);       
   
    TH1D* h_flux = GetBinNormalizedTH1D(hist);
    h_flux->SetLineWidth(3); 
    //h_flux->SetLineStyle(6);       
    h_flux->SetLineStyle(7);       
    //h_flux->SetLineColor(TColor::GetColor("#8b4513"));
    h_flux->SetLineColor(kGreen+3);
    
    TH1D* h_other = GetBinNormalizedTH1D(hist);
    h_other->SetLineWidth(3); 
    //h_other->SetLineStyle(7);       
    h_other->SetLineStyle(4);       
    //h_other->SetLineColor(kGreen+3);
    h_other->SetLineColor(kOrange+2);
    
    TLegend *legend = new TLegend(0.25,0.65,0.75,0.85);  
    ApplyStyle_Legend(legend);
    legend->SetNColumns(2);
    legend->SetTextSize(0.05);
    legend->AddEntry(h_total, "Total Error", "l" );
    legend->AddEntry(h_stat, "Statistical", "l" );
    //legend->SetColumnSeparation(0.0);
    legend->SetEntrySeparation(0.12);
    legend->AddEntry(h_xsec, "X-Sec Model", "l");
    legend->AddEntry(h_detector, "Detector", "l");
    legend->AddEntry(h_fsi, "FSI Model", "l");
    legend->AddEntry(h_flux, "Flux", "l");
    legend->AddEntry(h_other, "Other", "l");
    legend->Draw();

 
    double info_text_x = 0.85;
    double info_text_y = 0.85;
    TLatex info_text;
    info_text.SetNDC();
    info_text.SetTextColor(kBlack);
    info_text.SetTextFont(62);
    info_text.SetTextAlign(12);
    info_text.SetTextSize(0.04);
    info_text.DrawLatex(info_text_x, info_text_y, "(a)");

    // Print Plot
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + "Errors_" + var_name + ".pdf";
    c->Print(out_name.c_str(), "pdf");

    delete h_total;
    delete h_stat;
    delete h_detector;
    delete h_xsec;
    delete h_fsi;
    delete h_flux;
    delete h_other;

    delete c;
    delete plotter;
}

void CCProtonPi0_Plotter::DrawPaper_Enu(bool is1Track)
{
    std::string root_dir_data = Folder_List::rootDir_Interaction_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Interaction_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var = is1Track ? "Enu_1Track_0" : "Enu_2Track_0";

    MnvH1D* tempData = GetMnvH1D(f_data, var);
    MnvH1D* tempMC = GetMnvH1D(f_mc, var);
    
    tempData->GetXaxis()->SetNdivisions(5,4,0);
    tempMC->GetXaxis()->SetNdivisions(5,4,0);
    
    tempData->GetYaxis()->SetNdivisions(5,4,0);
    tempMC->GetYaxis()->SetNdivisions(5,4,0);
 
    tempData->GetYaxis()->CenterTitle();
    tempMC->GetYaxis()->CenterTitle();

    tempData->SetTitle("Data (3.33e20 POT)");
    tempMC->SetTitle("Simulation");
   
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
       
    gStyle->SetStripDecimals(false);
    gStyle->SetEndErrorSize(4);
    gStyle->SetCanvasDefW(480);
    gStyle->SetCanvasDefH(480); 
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.5;
    plotter->legend_text_size = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)

    plotter->data_line_width = 2;

    plotter->axis_minimum = 0.01;
    plotter->axis_maximum = 1700;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.2;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;
     
    TCanvas* c = new TCanvas("c");

    //    void MnvPlotter::DrawDataMCWithErrorBand    (   const MnvH1D *  dataHist,
    //            const MnvH1D *  mcHist,
    //            const Double_t  mcScale = 1.0,
    //            const std::string &     legPos = "L",
    //            const bool  useHistTitles = false,
    //            const MnvH1D *  bkgdHist = NULL,
    //            const MnvH1D *  dataBkgdHist = NULL,
    //            const bool  covAreaNormalize = false,
    //            const bool  statPlusSys = false  
    //            const bool isSmooth = false
    //            )   
    plotter->DrawDataMCWithErrorBand(tempData, tempMC, POT_ratio, "N", true, NULL, NULL, false, true, false);
 
    // ------------------------------------------------------------------------
    // Add Background Histogram
    // ------------------------------------------------------------------------
    std::string var_bckg = is1Track ? "Enu_1Track_2" : "Enu_2Track_2";
    MnvH1D* mnv_bckg = GetMnvH1D(f_mc, var_bckg);
    mnv_bckg->ClearAllErrorBands();
    std::cout<<"MnvBckg Integral = "<<mnv_bckg->Integral("width")<<std::endl;
    TH1D* bckg = GetBinNormalizedTH1D(mnv_bckg, false);
    std::cout<<"bckg Integral = "<<bckg->Integral()<<std::endl;
    bckg->Scale(POT_ratio);
    bckg->SetTitle("Background");
    bckg->SetLineWidth(2);
    bckg->SetLineColor(kGray+2);
    bckg->SetFillColor(kGray+2);
    bckg->SetFillStyle(3010);
    bckg->GetXaxis()->SetNdivisions(5,4,0);
    bckg->GetYaxis()->SetNdivisions(5,4,0);
    bckg->GetYaxis()->CenterTitle();
    bckg->Draw("HIST SAME");

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");
    text.DrawLatex(0.30,0.83,"(3.33e20 POT)");

    // Add Legend
    tempData->SetMarkerStyle(plotter->data_marker);
    tempData->SetMarkerSize(plotter->data_marker_size);
    tempData->SetMarkerColor(plotter->data_color);
    tempData->SetLineColor(plotter->data_color);
    tempData->SetLineWidth(plotter->data_line_width);

    tempMC->SetLineColor(plotter->mc_color);
    tempMC->SetLineWidth(plotter->mc_line_width);
    tempMC->SetLineStyle(plotter->mc_line_style);
    tempMC->SetFillColor(plotter->mc_error_color);
    tempMC->SetFillStyle(plotter->mc_error_style);

    double leg_x_min = 0.50;
    double leg_x_max = 0.85;
    double leg_y_min = 0.55;
    double leg_y_max = 0.80;

    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    ApplyStyle_Legend(legend);
    legend->SetTextSize(0.06);
    legend->AddEntry(tempData, "Data", "lep" );
    legend->AddEntry(tempMC, "Simulation", "fl");
    legend->AddEntry(bckg, "Background", "f");
    legend->Draw();

    double info_text_x = 0.85;
    double info_text_y = 0.85;
    TLatex info_text;
    info_text.SetNDC();
    info_text.SetTextColor(kBlack);
    info_text.SetTextFont(62);
    info_text.SetTextAlign(12);
    info_text.SetTextSize(0.04);
    if (is1Track) info_text.DrawLatex(info_text_x, info_text_y, "(a)");
    else info_text.DrawLatex(info_text_x, info_text_y, "(b)");


    // Plot Output
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var + ".pdf"; 
    c->Print(out_name.c_str(),"pdf");

    delete c;
    delete plotter;
    delete tempData;
    delete tempMC;
}


void CCProtonPi0_Plotter::DrawPaper_extra_total_energy(bool is1Track)
{
    std::string root_dir_data = Folder_List::rootDir_Interaction_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Interaction_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var = is1Track ? "extra_total_energy_1Track_0" : "extra_total_energy_2Track_0";

    MnvH1D* tempData = GetMnvH1D(f_data, var);
    MnvH1D* tempMC = GetMnvH1D(f_mc, var);
    
    tempData->GetXaxis()->SetNdivisions(5,4,0);
    tempMC->GetXaxis()->SetNdivisions(5,4,0);
    
    tempData->GetYaxis()->SetNdivisions(5,4,0);
    tempMC->GetYaxis()->SetNdivisions(5,4,0);
 
    tempData->GetYaxis()->CenterTitle();
    tempMC->GetYaxis()->CenterTitle();

    tempData->SetTitle("Data (3.33e20 POT)");
    tempMC->SetTitle("Simulation");
   
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
       
    gStyle->SetStripDecimals(false);
    gStyle->SetEndErrorSize(4);
    gStyle->SetCanvasDefW(480);
    gStyle->SetCanvasDefH(480); 
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.5;
    plotter->legend_text_size = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)

    plotter->data_line_width = 2;

    plotter->axis_minimum = 0.01;
    plotter->axis_maximum = 800;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.2;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;
     
    TCanvas* c = new TCanvas("c");

    //    void MnvPlotter::DrawDataMCWithErrorBand    (   const MnvH1D *  dataHist,
    //            const MnvH1D *  mcHist,
    //            const Double_t  mcScale = 1.0,
    //            const std::string &     legPos = "L",
    //            const bool  useHistTitles = false,
    //            const MnvH1D *  bkgdHist = NULL,
    //            const MnvH1D *  dataBkgdHist = NULL,
    //            const bool  covAreaNormalize = false,
    //            const bool  statPlusSys = false  
    //            const bool isSmooth = false
    //            )   
    plotter->DrawDataMCWithErrorBand(tempData, tempMC, POT_ratio, "N", true, NULL, NULL, false, true, false);
 
    // ------------------------------------------------------------------------
    // Add Background Histogram
    // ------------------------------------------------------------------------
    std::string var_bckg = is1Track ? "extra_total_energy_1Track_2" : "extra_total_energy_2Track_2";
    MnvH1D* mnv_bckg = GetMnvH1D(f_mc, var_bckg);
    mnv_bckg->ClearAllErrorBands();
    std::cout<<"MnvBckg Integral = "<<mnv_bckg->Integral("width")<<std::endl;
    TH1D* bckg = GetBinNormalizedTH1D(mnv_bckg, false);
    std::cout<<"bckg Integral = "<<bckg->Integral()<<std::endl;
    bckg->Scale(POT_ratio);
    bckg->SetTitle("Background");
    bckg->SetLineWidth(2);
    bckg->SetLineColor(kGray+2);
    bckg->SetFillColor(kGray+2);
    bckg->SetFillStyle(3010);
    bckg->GetXaxis()->SetNdivisions(5,4,0);
    bckg->GetYaxis()->SetNdivisions(5,4,0);
    bckg->GetYaxis()->CenterTitle();
    bckg->Draw("HIST SAME");

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");
    text.DrawLatex(0.30,0.83,"(3.33e20 POT)");

    // Add Legend
    tempData->SetMarkerStyle(plotter->data_marker);
    tempData->SetMarkerSize(plotter->data_marker_size);
    tempData->SetMarkerColor(plotter->data_color);
    tempData->SetLineColor(plotter->data_color);
    tempData->SetLineWidth(plotter->data_line_width);

    tempMC->SetLineColor(plotter->mc_color);
    tempMC->SetLineWidth(plotter->mc_line_width);
    tempMC->SetLineStyle(plotter->mc_line_style);
    tempMC->SetFillColor(plotter->mc_error_color);
    tempMC->SetFillStyle(plotter->mc_error_style);

    double leg_x_min = 0.50;
    double leg_x_max = 0.85;
    double leg_y_min = 0.55;
    double leg_y_max = 0.80;

    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    ApplyStyle_Legend(legend);
    legend->SetTextSize(0.06);
    legend->AddEntry(tempData, "Data", "lep" );
    legend->AddEntry(tempMC, "Simulation", "fl");
    legend->AddEntry(bckg, "Background", "f");
    legend->Draw();

    double info_text_x = 0.85;
    double info_text_y = 0.85;
    TLatex info_text;
    info_text.SetNDC();
    info_text.SetTextColor(kBlack);
    info_text.SetTextFont(62);
    info_text.SetTextAlign(12);
    info_text.SetTextSize(0.04);
    if (is1Track) info_text.DrawLatex(info_text_x, info_text_y, "(a)");
    else info_text.DrawLatex(info_text_x, info_text_y, "(b)");


    // Plot Output
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var + ".pdf"; 
    c->Print(out_name.c_str(),"pdf");

    delete c;
    delete plotter;
    delete tempData;
    delete tempMC;
}

void CCProtonPi0_Plotter::DrawPaper_vertex_energy(bool is1Track)
{
    std::string root_dir_data = Folder_List::rootDir_Interaction_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Interaction_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var = is1Track ? "vertex_energy_1Track_0" : "vertex_energy_2Track_0";

    MnvH1D* tempData = GetMnvH1D(f_data, var);
    MnvH1D* tempMC = GetMnvH1D(f_mc, var);
    
    tempData->GetXaxis()->SetNdivisions(5,4,0);
    tempMC->GetXaxis()->SetNdivisions(5,4,0);
    
    tempData->GetYaxis()->SetNdivisions(5,4,0);
    tempMC->GetYaxis()->SetNdivisions(5,4,0);
 
    tempData->GetYaxis()->CenterTitle();
    tempMC->GetYaxis()->CenterTitle();

    tempData->SetTitle("Data (3.33e20 POT)");
    tempMC->SetTitle("Simulation");
   
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
       
    gStyle->SetStripDecimals(false);
    gStyle->SetEndErrorSize(4);
    gStyle->SetCanvasDefW(480);
    gStyle->SetCanvasDefH(480); 
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.2;
    plotter->legend_text_font = 42; // default 62 (bold)

    plotter->data_line_width = 2;

    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    if(is1Track) plotter->axis_title_offset_y = 1.2;
    else plotter->axis_title_offset_y = 0.8; 
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;
     
    TCanvas* c = new TCanvas("c");

    //    void MnvPlotter::DrawDataMCWithErrorBand    (   const MnvH1D *  dataHist,
    //            const MnvH1D *  mcHist,
    //            const Double_t  mcScale = 1.0,
    //            const std::string &     legPos = "L",
    //            const bool  useHistTitles = false,
    //            const MnvH1D *  bkgdHist = NULL,
    //            const MnvH1D *  dataBkgdHist = NULL,
    //            const bool  covAreaNormalize = false,
    //            const bool  statPlusSys = false  
    //            const bool isSmooth = false
    //            )   
    plotter->DrawDataMCWithErrorBand(tempData, tempMC, POT_ratio, "N", true, NULL, NULL, false, true, false);
 
    // ------------------------------------------------------------------------
    // Add Background Histogram
    // ------------------------------------------------------------------------
    std::string var_bckg = is1Track ? "vertex_energy_1Track_2" : "vertex_energy_2Track_2";
    MnvH1D* mnv_bckg = GetMnvH1D(f_mc, var_bckg);
    mnv_bckg->ClearAllErrorBands();
    std::cout<<"MnvBckg Integral = "<<mnv_bckg->Integral("width")<<std::endl;
    TH1D* bckg = GetBinNormalizedTH1D(mnv_bckg, false);
    std::cout<<"bckg Integral = "<<bckg->Integral()<<std::endl;
    bckg->Scale(POT_ratio);
    bckg->SetTitle("Background");
    bckg->SetLineWidth(2);
    bckg->SetLineColor(kGray+2);
    bckg->SetFillColor(kGray+2);
    bckg->SetFillStyle(3010);
    bckg->GetXaxis()->SetNdivisions(5,4,0);
    bckg->GetYaxis()->SetNdivisions(5,4,0);
    bckg->GetYaxis()->CenterTitle();
    bckg->Draw("HIST SAME");

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");
    text.DrawLatex(0.30,0.83,"(3.33e20 POT)");

    // Add Legend
    tempData->SetMarkerStyle(plotter->data_marker);
    tempData->SetMarkerSize(plotter->data_marker_size);
    tempData->SetMarkerColor(plotter->data_color);
    tempData->SetLineColor(plotter->data_color);
    tempData->SetLineWidth(plotter->data_line_width);

    tempMC->SetLineColor(plotter->mc_color);
    tempMC->SetLineWidth(plotter->mc_line_width);
    tempMC->SetLineStyle(plotter->mc_line_style);
    tempMC->SetFillColor(plotter->mc_error_color);
    tempMC->SetFillStyle(plotter->mc_error_style);

    double leg_x_min = 0.50;
    double leg_x_max = 0.85;
    double leg_y_min = 0.55;
    double leg_y_max = 0.80;

    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    ApplyStyle_Legend(legend);

    legend->SetTextSize(0.06);
    legend->AddEntry(tempData, "Data", "lep" );
    legend->AddEntry(tempMC, "Simulation", "fl");
    legend->AddEntry(bckg, "Background", "f");
    legend->Draw();

    double info_text_x = 0.85;
    double info_text_y = 0.85;
    TLatex info_text;
    info_text.SetNDC();
    info_text.SetTextColor(kBlack);
    info_text.SetTextFont(62);
    info_text.SetTextAlign(12);
    info_text.SetTextSize(0.04);
    if (is1Track) info_text.DrawLatex(info_text_x, info_text_y, "(a)");
    else info_text.DrawLatex(info_text_x, info_text_y, "(b)");


    // Plot Output
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var + ".pdf"; 
    c->Print(out_name.c_str(),"pdf");

    delete c;
    delete plotter;
    delete tempData;
    delete tempMC;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_Enu_IntType()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "Enu";
    std::string data_var = var_name + "_xsec";
    std::string hist_name;

    MnvH1D* data = GetMnvH1D(f_data, data_var);
 
    data->GetXaxis()->SetNdivisions(5,4,0);
    data->GetYaxis()->SetNdivisions(5,4,0);
    data->GetYaxis()->CenterTitle();

    MnvH1D* temp = NULL;
    std::vector<MnvH1D*> mc_IntType;
    for (int i = 0; i < nIntType; ++i){
        hist_name = var_name + "_xsec_IntType_" + std::to_string((long long int)i);
        temp = GetMnvH1D(f_mc, hist_name);
        temp->ClearAllErrorBands();
        mc_IntType.push_back(temp);
    }

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 2.0;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    TCanvas* c = new TCanvas("c");

    MnvH1D* h_data = new MnvH1D (*data);
    TH1D* h_data_stat_only = GetBinNormalizedTH1D(h_data, false);
    h_data_stat_only->SetMarkerStyle(20);
    h_data_stat_only->SetMarkerSize(1);
    h_data_stat_only->SetMarkerColor(kBlack);
    h_data_stat_only->SetLineWidth(1);
    h_data_stat_only->SetLineColor(kBlack);
    h_data->GetXaxis()->SetNdivisions(5,4,0);
    h_data->GetYaxis()->SetNdivisions(5,4,0);
    h_data->GetYaxis()->CenterTitle();

    MnvH1D* h_RES_Delta = new MnvH1D (*mc_IntType[0]);
    h_RES_Delta->SetLineWidth(1);
    h_RES_Delta->SetLineColor(17);
    h_RES_Delta->SetFillColor(17);
    h_RES_Delta->SetFillStyle(1001);
    h_RES_Delta->GetXaxis()->SetNdivisions(5,4,0);
    h_RES_Delta->GetYaxis()->SetNdivisions(5,4,0);
    h_RES_Delta->GetYaxis()->CenterTitle();

    MnvH1D* h_RES_Other = new MnvH1D (*mc_IntType[1]);
    h_RES_Other->SetLineWidth(1);
    h_RES_Other->SetLineColor(kGreen+2);
    h_RES_Other->SetFillColor(kGreen+2);
    h_RES_Other->SetFillStyle(1001);
    h_RES_Other->GetXaxis()->SetNdivisions(5,4,0);
    h_RES_Other->GetYaxis()->SetNdivisions(5,4,0);
    h_RES_Other->GetYaxis()->CenterTitle();

    MnvH1D* h_NonRES = new MnvH1D (*mc_IntType[2]);
    h_NonRES->SetLineWidth(1);
    h_NonRES->SetLineColor(kRed+2);
    h_NonRES->SetFillColor(kRed+2);
    h_NonRES->SetFillStyle(1001);
    h_NonRES->GetXaxis()->SetNdivisions(5,4,0);
    h_NonRES->GetYaxis()->SetNdivisions(5,4,0);
    h_NonRES->GetYaxis()->CenterTitle();
    TObjArray* mc_hists = new TObjArray;
    mc_hists->Add(h_NonRES);
    mc_hists->Add(h_RES_Other);
    mc_hists->Add(h_RES_Delta);

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
    plotter->DrawDataStackedMC(h_data, mc_hists, 1.0, "N", "Data (3.33e20 POT)", 0, 0, 1001);

    h_data_stat_only->Draw("SAME E1 X0");

    // ------------------------------------------------------------------------
    // Legend 
    // ------------------------------------------------------------------------
    
    double leg_x_min = 0.55;
    double leg_x_max = 0.90;
    double leg_y_min = 0.62;
    double leg_y_max = 0.90;
    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);

    ApplyStyle_Legend(legend);
    legend->AddEntry(h_data_stat_only, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(h_RES_Delta, "Delta resonance","f");
    legend->AddEntry(h_RES_Other, "Other resonances","f");
    legend->AddEntry(h_NonRES, "Non-Resonant","f");
    legend->Draw();

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");

    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var_name + "_xsec_IntType" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete h_data;
    delete h_data_stat_only;
    delete h_RES_Delta;
    delete h_RES_Other;
    delete h_NonRES;
    delete mc_hists;
    delete legend;
    delete c;
    delete plotter;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_Delta_pi_phi_2Bin_IntType()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_2Bin_DeltaRich_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_2Bin_DeltaRich_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "Delta_pi_phi";
    std::string data_var = var_name + "_xsec";
    std::string hist_name;

    MnvH1D* data = GetMnvH1D(f_data, data_var);
 
    data->GetXaxis()->SetNdivisions(6,5,0);
    data->GetYaxis()->SetNdivisions(4,4,0);
    data->GetYaxis()->CenterTitle();

    bool isDebug = true;
    std::cout<<"---------------------------------------------------------------"<<std::endl;
    std::cout<<"Phi Asymmetry"<<std::endl;
    TH1D* err_data = new TH1D(data->GetTotalError());
    double n_right = data->GetBinContent(1);
    double err_right = err_data->GetBinContent(1);
    double n_left = data->GetBinContent(2);
    double err_left = err_data->GetBinContent(2);
    if (isDebug){
        std::cout<<"Right = "<<n_right<<" "<<err_right<<std::endl;
        std::cout<<"Left = "<<n_left<<" "<<err_left<<std::endl;
    }
    double diff = n_left - n_right;
    double err_diff = std::sqrt(std::pow(err_left,2) + std::pow(err_right,2)) ;
    double sum = n_left + n_right;
    double err_sum = std::sqrt(std::pow(err_left,2) + std::pow(err_right,2)) ;

    double asymmetry = diff/sum;
    double err_asymmetry = std::sqrt( std::pow(err_diff/diff,2) + std::pow(err_sum/sum,2) ) * asymmetry;

    if (isDebug){
        std::cout<<"Diff = "<<diff<<" "<<err_diff<<std::endl;
        std::cout<<"Sum = "<<sum<<" "<<err_sum<<std::endl;
    }
    std::cout<<"Asymmetry = "<<asymmetry<<" "<<err_asymmetry<<std::endl;
    std::cout<<"---------------------------------------------------------------"<<std::endl;


    MnvH1D* temp = NULL;
    std::vector<MnvH1D*> mc_IntType;
    for (int i = 0; i < nIntType; ++i){
        hist_name = var_name + "_xsec_IntType_" + std::to_string((long long int)i);
        temp = GetMnvH1D(f_mc, hist_name);
        temp->ClearAllErrorBands();
        mc_IntType.push_back(temp);
    }

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.001;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    TCanvas* c = new TCanvas("c");

    MnvH1D* h_data = new MnvH1D (*data);
    TH1D* h_data_stat_only = GetBinNormalizedTH1D(h_data, false);
    h_data_stat_only->SetMarkerStyle(20);
    h_data_stat_only->SetMarkerSize(1);
    h_data_stat_only->SetMarkerColor(kBlack);
    h_data_stat_only->SetLineWidth(1);
    h_data_stat_only->SetLineColor(kBlack);
    h_data->GetXaxis()->SetNdivisions(7,5,0);
    h_data->GetYaxis()->SetNdivisions(4,4,0);
    h_data->GetYaxis()->CenterTitle();

    MnvH1D* h_RES_Delta = new MnvH1D (*mc_IntType[0]);
    h_RES_Delta->SetLineWidth(1);
    h_RES_Delta->SetLineColor(17);
    h_RES_Delta->SetFillColor(17);
    h_RES_Delta->SetFillStyle(1001);
    h_RES_Delta->GetXaxis()->SetNdivisions(7,5,0);
    h_RES_Delta->GetYaxis()->SetNdivisions(4,4,0);
    h_RES_Delta->GetYaxis()->CenterTitle();

    MnvH1D* h_RES_Other = new MnvH1D (*mc_IntType[1]);
    h_RES_Other->SetLineWidth(1);
    h_RES_Other->SetLineColor(kGreen+2);
    h_RES_Other->SetFillColor(kGreen+2);
    h_RES_Other->SetFillStyle(1001);
    h_RES_Other->GetXaxis()->SetNdivisions(7,5,0);
    h_RES_Other->GetYaxis()->SetNdivisions(4,4,0);
    h_RES_Other->GetYaxis()->CenterTitle();

    MnvH1D* h_NonRES = new MnvH1D (*mc_IntType[2]);
    h_NonRES->SetLineWidth(1);
    h_NonRES->SetLineColor(kRed+2);
    h_NonRES->SetFillColor(kRed+2);
    h_NonRES->SetFillStyle(1001);
    h_NonRES->GetXaxis()->SetNdivisions(7,5,0);
    h_NonRES->GetYaxis()->SetNdivisions(4,4,0);
    h_NonRES->GetYaxis()->CenterTitle();
    TObjArray* mc_hists = new TObjArray;
    mc_hists->Add(h_NonRES);
    mc_hists->Add(h_RES_Other);
    mc_hists->Add(h_RES_Delta);

    MnvH1D* mc_total = GetMnvH1D(f_mc, data_var);
    double area_total = mc_total->Integral();
    std::cout<<var_name<<std::endl;
    std::cout<<"h_RES_Delta Area = "<<h_RES_Delta->Integral()/area_total*100<<std::endl;
    std::cout<<"h_RES_Other Area = "<<h_RES_Other->Integral()/area_total*100<<std::endl;
    std::cout<<"h_NonRES Area = "<<h_NonRES->Integral()/area_total*100<<std::endl;
    delete mc_total;

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
    plotter->DrawDataStackedMC(h_data, mc_hists, 1.0, "N", "Data (3.33e20 POT)", 0, 0, 1001);

    h_data_stat_only->Draw("SAME E1 X0");

    // ------------------------------------------------------------------------
    // No Legend for phi
    // ------------------------------------------------------------------------
   
    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(12);
    text.DrawLatex(0.20,0.87,"POT Normalized");
    text.DrawLatex(0.20,0.83,"T_{p} > 0.1 GeV");
    text.DrawLatex(0.20,0.79,"W_{exp} < 1.4 GeV");

    double info_text_x = 0.85;
    double info_text_y = 0.85;
    TLatex info_text;
    info_text.SetNDC();
    info_text.SetTextColor(kBlack);
    info_text.SetTextFont(62);
    info_text.SetTextAlign(12);
    info_text.SetTextSize(0.04);
    info_text.DrawLatex(info_text_x, info_text_y, "(b)");


    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var_name + "_2Bin_xsec_IntType" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete h_data;
    delete h_data_stat_only;
    delete h_RES_Delta;
    delete h_RES_Other;
    delete h_NonRES;
    delete mc_hists;
    delete c;
    delete plotter;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_Delta_pi_phi_IntType()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_DeltaRich_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_DeltaRich_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "Delta_pi_phi";
    std::string data_var = var_name + "_xsec";
    std::string hist_name;

    MnvH1D* data = GetMnvH1D(f_data, data_var);
 
    data->GetXaxis()->SetNdivisions(6,5,0);
    data->GetYaxis()->SetNdivisions(4,4,0);
    data->GetYaxis()->CenterTitle();

    MnvH1D* temp = NULL;
    std::vector<MnvH1D*> mc_IntType;
    for (int i = 0; i < nIntType; ++i){
        hist_name = var_name + "_xsec_IntType_" + std::to_string((long long int)i);
        temp = GetMnvH1D(f_mc, hist_name);
        temp->ClearAllErrorBands();
        mc_IntType.push_back(temp);
    }

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.001;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    TCanvas* c = new TCanvas("c");

    MnvH1D* h_data = new MnvH1D (*data);
    TH1D* h_data_stat_only = GetBinNormalizedTH1D(h_data, false);
    h_data_stat_only->SetMarkerStyle(20);
    h_data_stat_only->SetMarkerSize(1);
    h_data_stat_only->SetMarkerColor(kBlack);
    h_data_stat_only->SetLineWidth(1);
    h_data_stat_only->SetLineColor(kBlack);
    h_data->GetXaxis()->SetNdivisions(7,5,0);
    h_data->GetYaxis()->SetNdivisions(4,4,0);
    h_data->GetYaxis()->CenterTitle();

    MnvH1D* h_RES_Delta = new MnvH1D (*mc_IntType[0]);
    h_RES_Delta->SetLineWidth(2);
    h_RES_Delta->SetLineColor(kGray+1);
    h_RES_Delta->SetFillColor(kGray+1);
    h_RES_Delta->SetFillStyle(3002);
    h_RES_Delta->GetXaxis()->SetNdivisions(7,5,0);
    h_RES_Delta->GetYaxis()->SetNdivisions(4,4,0);
    h_RES_Delta->GetYaxis()->CenterTitle();

    MnvH1D* h_RES_Other = new MnvH1D (*mc_IntType[1]);
    h_RES_Other->SetLineWidth(2);
    h_RES_Other->SetLineColor(kGreen+3);
    h_RES_Other->SetFillColor(kGreen+3);
    h_RES_Other->SetFillStyle(3002);
    h_RES_Other->GetXaxis()->SetNdivisions(7,5,0);
    h_RES_Other->GetYaxis()->SetNdivisions(4,4,0);
    h_RES_Other->GetYaxis()->CenterTitle();

    MnvH1D* h_NonRES = new MnvH1D (*mc_IntType[2]);
    h_NonRES->SetLineWidth(2);
    h_NonRES->SetLineColor(kRed+4);
    h_NonRES->SetFillColor(kRed+4);
    h_NonRES->SetFillStyle(3002);
    h_NonRES->GetXaxis()->SetNdivisions(7,5,0);
    h_NonRES->GetYaxis()->SetNdivisions(4,4,0);
    h_NonRES->GetYaxis()->CenterTitle();
    TObjArray* mc_hists = new TObjArray;
    mc_hists->Add(h_NonRES);
    mc_hists->Add(h_RES_Other);
    mc_hists->Add(h_RES_Delta);

    MnvH1D* mc_total = GetMnvH1D(f_mc, data_var);
    double area_total = mc_total->Integral();
    std::cout<<var_name<<std::endl;
    std::cout<<"h_RES_Delta Area = "<<h_RES_Delta->Integral()/area_total*100<<std::endl;
    std::cout<<"h_RES_Other Area = "<<h_RES_Other->Integral()/area_total*100<<std::endl;
    std::cout<<"h_NonRES Area = "<<h_NonRES->Integral()/area_total*100<<std::endl;
    delete mc_total;

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
    plotter->DrawDataStackedMC(h_data, mc_hists, 1.0, "N", "Data (3.33e20 POT)", 0, 0, 3002);

    h_data_stat_only->Draw("SAME E1 X0");

    // ------------------------------------------------------------------------
    // No Legend for phi
    // ------------------------------------------------------------------------
   
    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(12);
    text.DrawLatex(0.20,0.87,"POT Normalized");
    text.DrawLatex(0.20,0.83,"T_{p} > 0.1 GeV");
    text.DrawLatex(0.20,0.79,"W_{exp} < 1.4 GeV");

    double info_text_x = 0.85;
    double info_text_y = 0.85;
    TLatex info_text;
    info_text.SetNDC();
    info_text.SetTextColor(kBlack);
    info_text.SetTextFont(62);
    info_text.SetTextAlign(12);
    info_text.SetTextSize(0.04);
    info_text.DrawLatex(info_text_x, info_text_y, "(b)");


    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var_name + "_xsec_IntType" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete h_data;
    delete h_data_stat_only;
    delete h_RES_Delta;
    delete h_RES_Other;
    delete h_NonRES;
    delete mc_hists;
    delete c;
    delete plotter;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_Delta_pi_theta_IntType()
{    
    std::string root_dir_data = Folder_List::rootDir_Signal_DeltaRich_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_DeltaRich_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "Delta_pi_theta";
    std::string data_var = var_name + "_xsec";
    std::string hist_name;

    MnvH1D* data = GetMnvH1D(f_data, data_var);
 
    data->GetXaxis()->SetNdivisions(5,4,0);
    data->GetYaxis()->SetNdivisions(5,4,0);
    data->GetYaxis()->CenterTitle();

    MnvH1D* temp = NULL;
    std::vector<MnvH1D*> mc_IntType;
    for (int i = 0; i < nIntType; ++i){
        hist_name = var_name + "_xsec_IntType_" + std::to_string((long long int)i);
        temp = GetMnvH1D(f_mc, hist_name);
        temp->ClearAllErrorBands();
        mc_IntType.push_back(temp);
    }

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    TCanvas* c = new TCanvas("c");

    MnvH1D* h_data = new MnvH1D (*data);
    TH1D* h_data_stat_only = GetBinNormalizedTH1D(h_data, false);
    h_data_stat_only->SetMarkerStyle(20);
    h_data_stat_only->SetMarkerSize(1);
    h_data_stat_only->SetMarkerColor(kBlack);
    h_data_stat_only->SetLineWidth(1);
    h_data_stat_only->SetLineColor(kBlack);
    h_data->GetXaxis()->SetNdivisions(5,4,0);
    h_data->GetYaxis()->SetNdivisions(5,4,0);
    h_data->GetYaxis()->CenterTitle();

    MnvH1D* h_RES_Delta = new MnvH1D (*mc_IntType[0]);
    h_RES_Delta->SetLineWidth(2);
    h_RES_Delta->SetLineColor(kGray+1);
    h_RES_Delta->SetFillColor(kGray+1);
    h_RES_Delta->SetFillStyle(3002);
    h_RES_Delta->GetXaxis()->SetNdivisions(5,4,0);
    h_RES_Delta->GetYaxis()->SetNdivisions(5,4,0);
    h_RES_Delta->GetYaxis()->CenterTitle();

    MnvH1D* h_RES_Other = new MnvH1D (*mc_IntType[1]);
    h_RES_Other->SetLineWidth(2);
    h_RES_Other->SetLineColor(kGreen+3);
    h_RES_Other->SetFillColor(kGreen+3);
    h_RES_Other->SetFillStyle(3002);
    h_RES_Other->GetXaxis()->SetNdivisions(5,4,0);
    h_RES_Other->GetYaxis()->SetNdivisions(5,4,0);
    h_RES_Other->GetYaxis()->CenterTitle();

    MnvH1D* h_NonRES = new MnvH1D (*mc_IntType[2]);
    h_NonRES->SetLineWidth(2);
    h_NonRES->SetLineColor(kRed+4);
    h_NonRES->SetFillColor(kRed+4);
    h_NonRES->SetFillStyle(3002);
    h_NonRES->GetXaxis()->SetNdivisions(5,4,0);
    h_NonRES->GetYaxis()->SetNdivisions(5,4,0);
    h_NonRES->GetYaxis()->CenterTitle();
    TObjArray* mc_hists = new TObjArray;
    mc_hists->Add(h_NonRES);
    mc_hists->Add(h_RES_Other);
    mc_hists->Add(h_RES_Delta);

    MnvH1D* mc_total = GetMnvH1D(f_mc, data_var);
    double area_total = mc_total->Integral();
    std::cout<<var_name<<std::endl;
    std::cout<<"h_RES_Delta Area = "<<h_RES_Delta->Integral()/area_total*100<<std::endl;
    std::cout<<"h_RES_Other Area = "<<h_RES_Other->Integral()/area_total*100<<std::endl;
    std::cout<<"h_NonRES Area = "<<h_NonRES->Integral()/area_total*100<<std::endl;
    delete mc_total;

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
    plotter->DrawDataStackedMC(h_data, mc_hists, 1.0, "N", "Data (3.33e20 POT)", 0, 0, 3002);

    h_data_stat_only->Draw("SAME E1 X0");

    // ------------------------------------------------------------------------
    // Legend 
    // ------------------------------------------------------------------------
    double leg_x_min = 0.45;
    double leg_x_max = 0.80;
    double leg_y_min = 0.59;
    double leg_y_max = 0.87;
    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);

    ApplyStyle_Legend(legend);
    legend->AddEntry(h_data_stat_only, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(h_RES_Delta, "Delta resonance","f");
    legend->AddEntry(h_RES_Other, "Other resonances","f");
    legend->AddEntry(h_NonRES, "Non-Resonant","f");
    legend->Draw();

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(12);
    text.DrawLatex(0.20,0.87,"POT Normalized");
    text.DrawLatex(0.20,0.83,"T_{p} > 0.1 GeV");
    text.DrawLatex(0.20,0.79,"W_{exp} < 1.4 GeV");

    double info_text_x = 0.85;
    double info_text_y = 0.85;
    TLatex info_text;
    info_text.SetNDC();
    info_text.SetTextColor(kBlack);
    info_text.SetTextFont(62);
    info_text.SetTextAlign(12);
    info_text.SetTextSize(0.04);
    info_text.DrawLatex(info_text_x, info_text_y, "(a)");


    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var_name + "_xsec_IntType" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete h_data;
    delete h_data_stat_only;
    delete h_RES_Delta;
    delete h_RES_Other;
    delete h_NonRES;
    delete mc_hists;
    delete legend;
    delete c;
    delete plotter;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_Delta_pi_theta_2Bin_IntType()
{    
    std::string root_dir_data = Folder_List::rootDir_Signal_2Bin_DeltaRich_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_2Bin_DeltaRich_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "Delta_pi_theta";
    std::string data_var = var_name + "_xsec";
    std::string hist_name;

    MnvH1D* data = GetMnvH1D(f_data, data_var);
 
    data->GetXaxis()->SetNdivisions(5,4,0);
    data->GetYaxis()->SetNdivisions(5,4,0);
    data->GetYaxis()->CenterTitle();

    bool isDebug = true;
    std::cout<<"---------------------------------------------------------------"<<std::endl;
    std::cout<<"cos(theta) Asymmetry"<<std::endl;
    TH1D* err_data = new TH1D(data->GetTotalError());
    double n_back = data->GetBinContent(1);
    double err_back = err_data->GetBinContent(1);
    double n_forward = data->GetBinContent(2);
    double err_forward = err_data->GetBinContent(2);
    if (isDebug){
        std::cout<<"Backward = "<<n_back<<" "<<err_back<<std::endl;
        std::cout<<"Forward = "<<n_forward<<" "<<err_forward<<std::endl;
    }
    double diff = n_forward - n_back;
    double err_diff = std::sqrt(std::pow(err_forward,2) + std::pow(err_back,2)) ;
    double sum = n_forward + n_back;
    double err_sum = std::sqrt(std::pow(err_forward,2) + std::pow(err_back,2)) ;

    double asymmetry = diff/sum;
    double err_asymmetry = std::sqrt( std::pow(err_diff/diff,2) + std::pow(err_sum/sum,2) ) * asymmetry;

    if (isDebug){
        std::cout<<"Diff = "<<diff<<" "<<err_diff<<std::endl;
        std::cout<<"Sum = "<<sum<<" "<<err_sum<<std::endl;
    }
    std::cout<<"Asymmetry = "<<asymmetry<<" "<<err_asymmetry<<std::endl;
    std::cout<<"---------------------------------------------------------------"<<std::endl;

    delete err_data; 
    MnvH1D* temp = NULL;
    std::vector<MnvH1D*> mc_IntType;
    for (int i = 0; i < nIntType; ++i){
        hist_name = var_name + "_xsec_IntType_" + std::to_string((long long int)i);
        temp = GetMnvH1D(f_mc, hist_name);
        temp->ClearAllErrorBands();
        mc_IntType.push_back(temp);
    }

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    TCanvas* c = new TCanvas("c");

    MnvH1D* h_data = new MnvH1D (*data);
    TH1D* h_data_stat_only = GetBinNormalizedTH1D(h_data, false);
    h_data_stat_only->SetMarkerStyle(20);
    h_data_stat_only->SetMarkerSize(1);
    h_data_stat_only->SetMarkerColor(kBlack);
    h_data_stat_only->SetLineWidth(1);
    h_data_stat_only->SetLineColor(kBlack);
    h_data->GetXaxis()->SetNdivisions(5,4,0);
    h_data->GetYaxis()->SetNdivisions(5,4,0);
    h_data->GetYaxis()->CenterTitle();

    MnvH1D* h_RES_Delta = new MnvH1D (*mc_IntType[0]);
    h_RES_Delta->SetLineWidth(1);
    h_RES_Delta->SetLineColor(17);
    h_RES_Delta->SetFillColor(17);
    h_RES_Delta->SetFillStyle(1001);
    h_RES_Delta->GetXaxis()->SetNdivisions(5,4,0);
    h_RES_Delta->GetYaxis()->SetNdivisions(5,4,0);
    h_RES_Delta->GetYaxis()->CenterTitle();

    MnvH1D* h_RES_Other = new MnvH1D (*mc_IntType[1]);
    h_RES_Other->SetLineWidth(1);
    h_RES_Other->SetLineColor(kGreen+2);
    h_RES_Other->SetFillColor(kGreen+2);
    h_RES_Other->SetFillStyle(1001);
    h_RES_Other->GetXaxis()->SetNdivisions(5,4,0);
    h_RES_Other->GetYaxis()->SetNdivisions(5,4,0);
    h_RES_Other->GetYaxis()->CenterTitle();

    MnvH1D* h_NonRES = new MnvH1D (*mc_IntType[2]);
    h_NonRES->SetLineWidth(1);
    h_NonRES->SetLineColor(kRed+2);
    h_NonRES->SetFillColor(kRed+2);
    h_NonRES->SetFillStyle(1001);
    h_NonRES->GetXaxis()->SetNdivisions(5,4,0);
    h_NonRES->GetYaxis()->SetNdivisions(5,4,0);
    h_NonRES->GetYaxis()->CenterTitle();
    TObjArray* mc_hists = new TObjArray;
    mc_hists->Add(h_NonRES);
    mc_hists->Add(h_RES_Other);
    mc_hists->Add(h_RES_Delta);

    MnvH1D* mc_total = GetMnvH1D(f_mc, data_var);
    double area_total = mc_total->Integral();
    std::cout<<var_name<<std::endl;
    std::cout<<"h_RES_Delta Area = "<<h_RES_Delta->Integral()/area_total*100<<std::endl;
    std::cout<<"h_RES_Other Area = "<<h_RES_Other->Integral()/area_total*100<<std::endl;
    std::cout<<"h_NonRES Area = "<<h_NonRES->Integral()/area_total*100<<std::endl;
    delete mc_total;

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
    plotter->DrawDataStackedMC(h_data, mc_hists, 1.0, "N", "Data (3.33e20 POT)", 0, 0, 1001);

    h_data_stat_only->Draw("SAME E1 X0");

    // ------------------------------------------------------------------------
    // Legend 
    // ------------------------------------------------------------------------
    double leg_x_min = 0.45;
    double leg_x_max = 0.80;
    double leg_y_min = 0.65;
    double leg_y_max = 0.90;
    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);

    ApplyStyle_Legend(legend);
    legend->AddEntry(h_data_stat_only, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(h_RES_Delta, "Delta resonance","f");
    legend->AddEntry(h_RES_Other, "Other resonances","f");
    legend->AddEntry(h_NonRES, "Non-Resonant","f");
    legend->Draw();

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(12);
    text.DrawLatex(0.20,0.87,"POT Normalized");
    text.DrawLatex(0.20,0.83,"T_{p} > 0.1 GeV");
    text.DrawLatex(0.20,0.79,"W_{exp} < 1.4 GeV");

    double info_text_x = 0.85;
    double info_text_y = 0.85;
    TLatex info_text;
    info_text.SetNDC();
    info_text.SetTextColor(kBlack);
    info_text.SetTextFont(62);
    info_text.SetTextAlign(12);
    info_text.SetTextSize(0.04);
    info_text.DrawLatex(info_text_x, info_text_y, "(a)");


    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var_name + "_2Bin_xsec_IntType" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete h_data;
    delete h_data_stat_only;
    delete h_RES_Delta;
    delete h_RES_Other;
    delete h_NonRES;
    delete mc_hists;
    delete legend;
    delete c;
    delete plotter;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_W_NoLimit_IntType()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_NoWLimit_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_NoWLimit_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "W";
    std::string data_var = var_name + "_xsec";
    std::string hist_name;

    MnvH1D* data = GetMnvH1D(f_data, data_var);
 
    //data->GetXaxis()->SetNdivisions(5,4,0);
    data->GetYaxis()->SetNdivisions(5,4,0);
    data->GetYaxis()->CenterTitle();
    
    data->GetXaxis()->SetRangeUser(0.6,3.0);    
    MnvH1D* temp = NULL;
    std::vector<MnvH1D*> mc_IntType;
    for (int i = 0; i < nIntType; ++i){
        hist_name = var_name + "_xsec_IntType_" + std::to_string((long long int)i);
        temp = GetMnvH1D(f_mc, hist_name);
        temp->ClearAllErrorBands();
        temp->GetXaxis()->SetRangeUser(0.6,3.0);    
        mc_IntType.push_back(temp);
    }

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 2.0;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    TCanvas* c = new TCanvas("c");

    MnvH1D* h_data = new MnvH1D (*data);
    TH1D* h_data_stat_only = GetBinNormalizedTH1D(h_data, false);
    h_data_stat_only->SetMarkerStyle(20);
    h_data_stat_only->SetMarkerSize(1);
    h_data_stat_only->SetMarkerColor(kBlack);
    h_data_stat_only->SetLineWidth(1);
    h_data_stat_only->SetLineColor(kBlack);
    //h_data->GetXaxis()->SetNdivisions(5,4,0);
    h_data->GetYaxis()->SetNdivisions(5,4,0);
    h_data->GetYaxis()->CenterTitle();

    MnvH1D* h_RES_Delta = new MnvH1D (*mc_IntType[0]);
    h_RES_Delta->SetLineWidth(2);
    h_RES_Delta->SetLineColor(kGray+1);
    h_RES_Delta->SetFillColor(kGray+1);
    h_RES_Delta->SetFillStyle(3002);
    //h_RES_Delta->GetXaxis()->SetNdivisions(5,4,0);
    h_RES_Delta->GetYaxis()->SetNdivisions(5,4,0);
    h_RES_Delta->GetYaxis()->CenterTitle();

    MnvH1D* h_RES_Other = new MnvH1D (*mc_IntType[1]);
    h_RES_Other->SetLineWidth(2);
    h_RES_Other->SetLineColor(kGreen+3);
    h_RES_Other->SetFillColor(kGreen+3);
    h_RES_Other->SetFillStyle(3002);
    //h_RES_Other->GetXaxis()->SetNdivisions(5,4,0);
    h_RES_Other->GetYaxis()->SetNdivisions(5,4,0);
    h_RES_Other->GetYaxis()->CenterTitle();

    MnvH1D* h_NonRES = new MnvH1D (*mc_IntType[2]);
    h_NonRES->SetLineWidth(2);
    h_NonRES->SetLineColor(kRed+4);
    h_NonRES->SetFillColor(kRed+4);
    h_NonRES->SetFillStyle(3002);
    //h_NonRES->GetXaxis()->SetNdivisions(5,4,0);
    h_NonRES->GetYaxis()->SetNdivisions(5,4,0);
    h_NonRES->GetYaxis()->CenterTitle();
    TObjArray* mc_hists = new TObjArray;
    mc_hists->Add(h_NonRES);
    mc_hists->Add(h_RES_Other);
    mc_hists->Add(h_RES_Delta);

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
    plotter->DrawDataStackedMC(h_data, mc_hists, 1.0, "N", "Data (3.33e20 POT)", 0, 0, 3002);

    h_data_stat_only->Draw("SAME E1 X0");

    TLine line;
    TArrow arrow;
    line.SetLineWidth(2);
    line.SetLineStyle(2);
    line.SetLineColor(kBlue);
    arrow.SetLineWidth(2);
    arrow.SetLineColor(kBlue);
    line.DrawLine(1.8, 0.0, 1.8, 38);
    arrow.DrawArrow(1.8, 38, 1.7, 38, 0.01);

    // ------------------------------------------------------------------------
    // Legend 
    // ------------------------------------------------------------------------
    double leg_x_min = 0.60;
    double leg_x_max = 0.95;
    double leg_y_min = 0.62;
    double leg_y_max = 0.90;   
    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);

    ApplyStyle_Legend(legend);
    legend->AddEntry(h_data_stat_only, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(h_RES_Delta, "Delta resonance","f");
    legend->AddEntry(h_RES_Other, "Other resonances","f");
    legend->AddEntry(h_NonRES, "Non-Resonant","f");
    legend->Draw();

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");

    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var_name + "_NoLimit_xsec_IntType" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete h_data;
    delete h_data_stat_only;
    delete h_RES_Delta;
    delete h_RES_Other;
    delete h_NonRES;
    delete mc_hists;
    delete legend;
    delete c;
    delete plotter;
}


void CCProtonPi0_Plotter::DrawPaper_xsec_W_IntType()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "W";
    std::string data_var = var_name + "_xsec";
    std::string hist_name;

    MnvH1D* data = GetMnvH1D(f_data, data_var);
 
    //data->GetXaxis()->SetNdivisions(5,4,0);
    data->GetYaxis()->SetNdivisions(5,4,0);
    data->GetYaxis()->CenterTitle();

    MnvH1D* temp = NULL;
    std::vector<MnvH1D*> mc_IntType;
    for (int i = 0; i < nIntType; ++i){
        hist_name = var_name + "_xsec_IntType_" + std::to_string((long long int)i);
        temp = GetMnvH1D(f_mc, hist_name);
        temp->ClearAllErrorBands();
        mc_IntType.push_back(temp);
    }

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 2.0;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    TCanvas* c = new TCanvas("c");

    MnvH1D* h_data = new MnvH1D (*data);
    TH1D* h_data_stat_only = GetBinNormalizedTH1D(h_data, false);
    h_data_stat_only->SetMarkerStyle(20);
    h_data_stat_only->SetMarkerSize(1);
    h_data_stat_only->SetMarkerColor(kBlack);
    h_data_stat_only->SetLineWidth(1);
    h_data_stat_only->SetLineColor(kBlack);
    //h_data->GetXaxis()->SetNdivisions(5,4,0);
    h_data->GetYaxis()->SetNdivisions(5,4,0);
    h_data->GetYaxis()->CenterTitle();

    MnvH1D* h_RES_Delta = new MnvH1D (*mc_IntType[0]);
    h_RES_Delta->SetLineWidth(2);
    h_RES_Delta->SetLineColor(kGray+1);
    h_RES_Delta->SetFillColor(kGray+1);
    h_RES_Delta->SetFillStyle(3002);
    //h_RES_Delta->GetXaxis()->SetNdivisions(5,4,0);
    h_RES_Delta->GetYaxis()->SetNdivisions(5,4,0);
    h_RES_Delta->GetYaxis()->CenterTitle();

    MnvH1D* h_RES_Other = new MnvH1D (*mc_IntType[1]);
    h_RES_Other->SetLineWidth(2);
    h_RES_Other->SetLineColor(kGreen+3);
    h_RES_Other->SetFillColor(kGreen+3);
    h_RES_Other->SetFillStyle(3002);
    //h_RES_Other->GetXaxis()->SetNdivisions(5,4,0);
    h_RES_Other->GetYaxis()->SetNdivisions(5,4,0);
    h_RES_Other->GetYaxis()->CenterTitle();

    MnvH1D* h_NonRES = new MnvH1D (*mc_IntType[2]);
    h_NonRES->SetLineWidth(2);
    h_NonRES->SetLineColor(kRed+4);
    h_NonRES->SetFillColor(kRed+4);
    h_NonRES->SetFillStyle(3002);
    //h_NonRES->GetXaxis()->SetNdivisions(5,4,0);
    h_NonRES->GetYaxis()->SetNdivisions(5,4,0);
    h_NonRES->GetYaxis()->CenterTitle();
    TObjArray* mc_hists = new TObjArray;
    mc_hists->Add(h_NonRES);
    mc_hists->Add(h_RES_Other);
    mc_hists->Add(h_RES_Delta);

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
    plotter->DrawDataStackedMC(h_data, mc_hists, 1.0, "N", "Data (3.33e20 POT)", 0, 0, 3002);

    h_data_stat_only->Draw("SAME E1 X0");

    // ------------------------------------------------------------------------
    // Legend 
    // ------------------------------------------------------------------------
    
    double leg_x_min = 0.2;
    double leg_x_max = 0.57;
    double leg_y_min = 0.57;
    double leg_y_max = 0.82;
    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);

    ApplyStyle_Legend(legend);
    legend->AddEntry(h_data_stat_only, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(h_RES_Delta, "Delta resonance","f");
    legend->AddEntry(h_RES_Other, "Other resonances","f");
    legend->AddEntry(h_NonRES, "Non-Resonant","f");
    legend->Draw();

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");

    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var_name + "_xsec_IntType" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete h_data;
    delete h_data_stat_only;
    delete h_RES_Delta;
    delete h_RES_Other;
    delete h_NonRES;
    delete mc_hists;
    delete legend;
    delete c;
    delete plotter;
}


void CCProtonPi0_Plotter::DrawPaper_xsec_deltaInvMass2_IntType()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_TwoTrack_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_TwoTrack_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "deltaInvMass";
    std::string data_var = var_name + "_xsec";
    std::string hist_name;

    MnvH1D* data = GetMnvH1D(f_data, data_var);
 
    data->GetXaxis()->SetNdivisions(5,4,0);
    data->GetYaxis()->SetNdivisions(5,4,0);
    data->GetYaxis()->CenterTitle();

    MnvH1D* temp = NULL;
    std::vector<MnvH1D*> mc_IntType;
    for (int i = 0; i < nIntType; ++i){
        hist_name = var_name + "_xsec_IntType_" + std::to_string((long long int)i);
        temp = GetMnvH1D(f_mc, hist_name);
        temp->ClearAllErrorBands();
        mc_IntType.push_back(temp);
    }

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    plotter->mc_line_width = 0.0;
    TCanvas* c = new TCanvas("c");

    MnvH1D* h_data = new MnvH1D (*data);
    TH1D* h_data_stat_only = GetBinNormalizedTH1D(h_data, false);
    h_data_stat_only->SetMarkerStyle(20);
    h_data_stat_only->SetMarkerSize(1);
    h_data_stat_only->SetMarkerColor(kBlack);
    h_data_stat_only->SetLineWidth(1);
    h_data_stat_only->SetLineColor(kBlack);
    h_data->GetXaxis()->SetNdivisions(5,4,0);
    h_data->GetYaxis()->SetNdivisions(5,4,0);
    h_data->GetYaxis()->CenterTitle();

    MnvH1D* h_RES_Delta = new MnvH1D (*mc_IntType[0]);
    h_RES_Delta->SetLineWidth(1);
    h_RES_Delta->SetLineColor(17);
    h_RES_Delta->SetFillColor(17);
    h_RES_Delta->SetFillStyle(1001);
    h_RES_Delta->GetXaxis()->SetNdivisions(5,4,0);
    h_RES_Delta->GetYaxis()->SetNdivisions(5,4,0);
    h_RES_Delta->GetYaxis()->CenterTitle();

    MnvH1D* h_RES_Other = new MnvH1D (*mc_IntType[1]);
    h_RES_Other->SetLineWidth(1);
    h_RES_Other->SetLineColor(kGreen+2);
    h_RES_Other->SetFillColor(kGreen+2);
    h_RES_Other->SetFillStyle(1001);
    h_RES_Other->GetXaxis()->SetNdivisions(5,4,0);
    h_RES_Other->GetYaxis()->SetNdivisions(5,4,0);
    h_RES_Other->GetYaxis()->CenterTitle();

    MnvH1D* h_NonRES = new MnvH1D (*mc_IntType[2]);
    h_NonRES->SetLineWidth(1);
    h_NonRES->SetLineColor(kRed+2);
    h_NonRES->SetFillColor(kRed+2);
    h_NonRES->SetFillStyle(1001);
    h_NonRES->GetXaxis()->SetNdivisions(5,4,0);
    h_NonRES->GetYaxis()->SetNdivisions(5,4,0);
    h_NonRES->GetYaxis()->CenterTitle();
    TObjArray* mc_hists = new TObjArray;
    mc_hists->Add(h_NonRES);
    mc_hists->Add(h_RES_Other);
    mc_hists->Add(h_RES_Delta);

    MnvH1D* mc_total = GetMnvH1D(f_mc, data_var);
    double area_total = mc_total->Integral();
    std::cout<<var_name<<std::endl;
    std::cout<<"h_RES_Delta Area = "<<h_RES_Delta->Integral()/area_total*100<<std::endl;
    std::cout<<"h_RES_Other Area = "<<h_RES_Other->Integral()/area_total*100<<std::endl;
    std::cout<<"h_NonRES Area = "<<h_NonRES->Integral()/area_total*100<<std::endl;
    delete mc_total;

    //printBins(h_RES_Delta,"h_RES_Delta", true);
    //printBins(h_RES_Other,"h_RES_Other", true);
    //printBins(h_NonRES,"h_NonRES", true);
    
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
    plotter->DrawDataStackedMC(h_data, mc_hists, 1.0, "N", "Data (3.33e20 POT)", 0, 0, 1001);

    h_data_stat_only->Draw("SAME E1 X0");

    // ------------------------------------------------------------------------
    // Legend 
    // ------------------------------------------------------------------------
    double leg_x_min = 0.55;
    double leg_x_max = 0.85;
    double leg_y_min = 0.55;
    double leg_y_max = 0.80;
    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);

    ApplyStyle_Legend(legend);
    legend->AddEntry(h_data_stat_only, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(h_RES_Delta, "Delta resonance","f");
    legend->AddEntry(h_RES_Other, "Other resonances","f");
    legend->AddEntry(h_NonRES, "Non-Resonant","f");
    legend->Draw();

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(12);
    text.DrawLatex(0.20,0.87,"POT Normalized");
    text.DrawLatex(0.20,0.83,"T_{p} > 0.1 GeV");

    double info_text_x = 0.85;
    double info_text_y = 0.85;
    TLatex info_text;
    info_text.SetNDC();
    info_text.SetTextColor(kBlack);
    info_text.SetTextFont(62);
    info_text.SetTextAlign(12);
    info_text.SetTextSize(0.04);
    info_text.DrawLatex(info_text_x, info_text_y, "(b)");

    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var_name + "2_xsec_IntType" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete h_data;
    delete h_data_stat_only;
    delete h_RES_Delta;
    delete h_RES_Other;
    delete h_NonRES;
    delete mc_hists;
    delete legend;
    delete c;
    delete plotter;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_deltaInvMass_IntType()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_DeltaRich_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_DeltaRich_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "deltaInvMass";
    std::string data_var = var_name + "_xsec";
    std::string hist_name;

    MnvH1D* data = GetMnvH1D(f_data, data_var);
 
    data->GetXaxis()->SetNdivisions(5,4,0);
    data->GetYaxis()->SetNdivisions(5,4,0);
    data->GetYaxis()->CenterTitle();

    MnvH1D* temp = NULL;
    std::vector<MnvH1D*> mc_IntType;
    for (int i = 0; i < nIntType; ++i){
        hist_name = var_name + "_xsec_IntType_" + std::to_string((long long int)i);
        temp = GetMnvH1D(f_mc, hist_name);
        temp->ClearAllErrorBands();
        mc_IntType.push_back(temp);
    }

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    plotter->mc_line_width = 2;
    TCanvas* c = new TCanvas("c");

    MnvH1D* h_data = new MnvH1D (*data);
    TH1D* h_data_stat_only = GetBinNormalizedTH1D(h_data, false);
    h_data_stat_only->SetMarkerStyle(20);
    h_data_stat_only->SetMarkerSize(1);
    h_data_stat_only->SetMarkerColor(kBlack);
    h_data_stat_only->SetLineWidth(1);
    h_data_stat_only->SetLineColor(kBlack);
    h_data->GetXaxis()->SetNdivisions(5,4,0);
    h_data->GetYaxis()->SetNdivisions(5,4,0);
    h_data->GetYaxis()->CenterTitle();

    MnvH1D* h_RES_Delta = new MnvH1D (*mc_IntType[0]);
    h_RES_Delta->SetLineWidth(2);
    h_RES_Delta->SetLineColor(kGray+1);
    h_RES_Delta->SetFillColor(kGray+1);
    h_RES_Delta->SetFillStyle(3002);
    h_RES_Delta->GetXaxis()->SetNdivisions(5,4,0);
    h_RES_Delta->GetYaxis()->SetNdivisions(5,4,0);
    h_RES_Delta->GetYaxis()->CenterTitle();

    MnvH1D* h_RES_Other = new MnvH1D (*mc_IntType[1]);
    h_RES_Other->SetLineWidth(2);
    h_RES_Other->SetLineColor(kGreen+3);
    h_RES_Other->SetFillColor(kGreen+3);
    h_RES_Other->SetFillStyle(3002);
    h_RES_Other->GetXaxis()->SetNdivisions(5,4,0);
    h_RES_Other->GetYaxis()->SetNdivisions(5,4,0);
    h_RES_Other->GetYaxis()->CenterTitle();

    MnvH1D* h_NonRES = new MnvH1D (*mc_IntType[2]);
    h_NonRES->SetLineWidth(2);
    h_NonRES->SetLineColor(kRed+4);
    h_NonRES->SetFillColor(kRed+4);
    h_NonRES->SetFillStyle(3002);
    h_NonRES->GetXaxis()->SetNdivisions(5,4,0);
    h_NonRES->GetYaxis()->SetNdivisions(5,4,0);
    h_NonRES->GetYaxis()->CenterTitle();
    TObjArray* mc_hists = new TObjArray;
    mc_hists->Add(h_NonRES);
    mc_hists->Add(h_RES_Other);
    mc_hists->Add(h_RES_Delta);

    MnvH1D* mc_total = GetMnvH1D(f_mc, data_var);
    double area_total = mc_total->Integral();
    std::cout<<var_name<<std::endl;
    std::cout<<"h_RES_Delta Area = "<<h_RES_Delta->Integral()/area_total*100<<std::endl;
    std::cout<<"h_RES_Other Area = "<<h_RES_Other->Integral()/area_total*100<<std::endl;
    std::cout<<"h_NonRES Area = "<<h_NonRES->Integral()/area_total*100<<std::endl;
    delete mc_total;

    //printBins(h_RES_Delta,"h_RES_Delta", true);
    //printBins(h_RES_Other,"h_RES_Other", true);
    //printBins(h_NonRES,"h_NonRES", true);
    
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
    plotter->DrawDataStackedMC(h_data, mc_hists, 1.0, "N", "Data (3.33e20 POT)", 0, 0, 3002);

    h_data_stat_only->Draw("SAME E1 X0");

    // ------------------------------------------------------------------------
    // Legend 
    // ------------------------------------------------------------------------
    double leg_x_min = 0.55;
    double leg_x_max = 0.85;
    double leg_y_min = 0.55;
    double leg_y_max = 0.80;
    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);

    ApplyStyle_Legend(legend);
    legend->AddEntry(h_data_stat_only, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(h_RES_Delta, "Delta resonance","f");
    legend->AddEntry(h_RES_Other, "Other resonances","f");
    legend->AddEntry(h_NonRES, "Non-Resonant","f");
    legend->Draw();

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(12);
    text.DrawLatex(0.20,0.87,"POT Normalized");
    text.DrawLatex(0.20,0.83,"T_{p} > 0.1 GeV");
    text.DrawLatex(0.20,0.79,"W_{exp} < 1.4 GeV");

    double info_text_x = 0.85;
    double info_text_y = 0.85;
    TLatex info_text;
    info_text.SetNDC();
    info_text.SetTextColor(kBlack);
    info_text.SetTextFont(62);
    info_text.SetTextAlign(12);
    info_text.SetTextSize(0.04);
    info_text.DrawLatex(info_text_x, info_text_y, "(b)");

    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var_name + "_xsec_IntType" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete h_data;
    delete h_data_stat_only;
    delete h_RES_Delta;
    delete h_RES_Other;
    delete h_NonRES;
    delete mc_hists;
    delete legend;
    delete c;
    delete plotter;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_QSq_IntType()
{    
    std::string root_dir_data = Folder_List::rootDir_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "QSq";
    std::string data_var = var_name + "_xsec";
    std::string hist_name;

    MnvH1D* data = GetMnvH1D(f_data, data_var);
 
    data->GetXaxis()->SetNdivisions(5,4,0);
    data->GetYaxis()->SetNdivisions(5,4,0);
    data->GetYaxis()->CenterTitle();

    MnvH1D* temp = NULL;
    std::vector<MnvH1D*> mc_IntType;
    for (int i = 0; i < nIntType; ++i){
        hist_name = var_name + "_xsec_IntType_" + std::to_string((long long int)i);
        temp = GetMnvH1D(f_mc, hist_name);
        temp->ClearAllErrorBands();
        mc_IntType.push_back(temp);
    }

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    TCanvas* c = new TCanvas("c");

    MnvH1D* h_data = new MnvH1D (*data);
    TH1D* h_data_stat_only = GetBinNormalizedTH1D(h_data, false);
    h_data_stat_only->SetMarkerStyle(20);
    h_data_stat_only->SetMarkerSize(1);
    h_data_stat_only->SetMarkerColor(kBlack);
    h_data_stat_only->SetLineWidth(1);
    h_data_stat_only->SetLineColor(kBlack);
    h_data->GetXaxis()->SetNdivisions(5,4,0);
    h_data->GetYaxis()->SetNdivisions(5,4,0);
    h_data->GetYaxis()->CenterTitle();

    MnvH1D* h_RES_Delta = new MnvH1D (*mc_IntType[0]);
    h_RES_Delta->SetLineWidth(2);
    h_RES_Delta->SetLineColor(kGray+1);
    h_RES_Delta->SetFillColor(kGray+1);
    h_RES_Delta->SetFillStyle(3002);
    h_RES_Delta->GetXaxis()->SetNdivisions(5,4,0);
    h_RES_Delta->GetYaxis()->SetNdivisions(5,4,0);
    h_RES_Delta->GetYaxis()->CenterTitle();

    MnvH1D* h_RES_Other = new MnvH1D (*mc_IntType[1]);
    h_RES_Other->SetLineWidth(2);
    h_RES_Other->SetLineColor(kGreen+3);
    h_RES_Other->SetFillColor(kGreen+3);
    h_RES_Other->SetFillStyle(3002);
    h_RES_Other->GetXaxis()->SetNdivisions(5,4,0);
    h_RES_Other->GetYaxis()->SetNdivisions(5,4,0);
    h_RES_Other->GetYaxis()->CenterTitle();

    MnvH1D* h_NonRES = new MnvH1D (*mc_IntType[2]);
    h_NonRES->SetLineWidth(2);
    h_NonRES->SetLineColor(kRed+4);
    h_NonRES->SetFillColor(kRed+4);
    h_NonRES->SetFillStyle(3002);
    h_NonRES->GetXaxis()->SetNdivisions(5,4,0);
    h_NonRES->GetYaxis()->SetNdivisions(5,4,0);
    h_NonRES->GetYaxis()->CenterTitle();
    TObjArray* mc_hists = new TObjArray;
    mc_hists->Add(h_NonRES);
    mc_hists->Add(h_RES_Other);
    mc_hists->Add(h_RES_Delta);

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
    plotter->DrawDataStackedMC(h_data, mc_hists, 1.0, "N", "Data (3.33e20 POT)", 0, 0, 3002);

    h_data_stat_only->Draw("SAME E1 X0");

    // ------------------------------------------------------------------------
    // Legend 
    // ------------------------------------------------------------------------
    
    double leg_x_min = 0.55;
    double leg_x_max = 0.90;
    double leg_y_min = 0.57;
    double leg_y_max = 0.85;
    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);

    ApplyStyle_Legend(legend);
    legend->AddEntry(h_data_stat_only, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(h_RES_Delta, "Delta resonance","f");
    legend->AddEntry(h_RES_Other, "Other resonances","f");
    legend->AddEntry(h_NonRES, "Non-Resonant","f");
    legend->Draw();

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");

    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var_name + "_xsec_IntType" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete h_data;
    delete h_data_stat_only;
    delete h_RES_Delta;
    delete h_RES_Other;
    delete h_NonRES;
    delete mc_hists;
    delete legend;
    delete c;
    delete plotter;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_pi0_theta_IntType()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "pi0_theta";
    std::string data_var = var_name + "_xsec";
    std::string hist_name;

    MnvH1D* data = GetMnvH1D(f_data, data_var);
 
    data->GetXaxis()->SetNdivisions(5,4,0);
    data->GetYaxis()->SetNdivisions(5,4,0);
    data->GetYaxis()->CenterTitle();

    MnvH1D* temp = NULL;
    std::vector<MnvH1D*> mc_IntType;
    for (int i = 0; i < nIntType; ++i){
        hist_name = var_name + "_xsec_IntType_" + std::to_string((long long int)i);
        temp = GetMnvH1D(f_mc, hist_name);
        temp->ClearAllErrorBands();
        mc_IntType.push_back(temp);
    }

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    TCanvas* c = new TCanvas("c");

    MnvH1D* h_data = new MnvH1D (*data);
    TH1D* h_data_stat_only = GetBinNormalizedTH1D(h_data, false);
    h_data_stat_only->SetMarkerStyle(20);
    h_data_stat_only->SetMarkerSize(1);
    h_data_stat_only->SetMarkerColor(kBlack);
    h_data_stat_only->SetLineWidth(1);
    h_data_stat_only->SetLineColor(kBlack);
    h_data->GetXaxis()->SetNdivisions(5,4,0);
    h_data->GetYaxis()->SetNdivisions(5,4,0);
    h_data->GetYaxis()->CenterTitle();

    MnvH1D* h_RES_Delta = new MnvH1D (*mc_IntType[0]);
    h_RES_Delta->SetLineWidth(1);
    h_RES_Delta->SetLineColor(17);
    h_RES_Delta->SetFillColor(17);
    h_RES_Delta->SetFillStyle(1001);
    h_RES_Delta->GetXaxis()->SetNdivisions(5,4,0);
    h_RES_Delta->GetYaxis()->SetNdivisions(5,4,0);
    h_RES_Delta->GetYaxis()->CenterTitle();

    MnvH1D* h_RES_Other = new MnvH1D (*mc_IntType[1]);
    h_RES_Other->SetLineWidth(1);
    h_RES_Other->SetLineColor(kGreen+2);
    h_RES_Other->SetFillColor(kGreen+2);
    h_RES_Other->SetFillStyle(1001);
    h_RES_Other->GetXaxis()->SetNdivisions(5,4,0);
    h_RES_Other->GetYaxis()->SetNdivisions(5,4,0);
    h_RES_Other->GetYaxis()->CenterTitle();

    MnvH1D* h_NonRES = new MnvH1D (*mc_IntType[2]);
    h_NonRES->SetLineWidth(1);
    h_NonRES->SetLineColor(kRed+2);
    h_NonRES->SetFillColor(kRed+2);
    h_NonRES->SetFillStyle(1001);
    h_NonRES->GetXaxis()->SetNdivisions(5,4,0);
    h_NonRES->GetYaxis()->SetNdivisions(5,4,0);
    h_NonRES->GetYaxis()->CenterTitle();
    TObjArray* mc_hists = new TObjArray;
    mc_hists->Add(h_NonRES);
    mc_hists->Add(h_RES_Other);
    mc_hists->Add(h_RES_Delta);

    MnvH1D* mc_total = GetMnvH1D(f_mc, data_var);
    double area_total = mc_total->Integral();
    std::cout<<"h_RES_Delta Area = "<<h_RES_Delta->Integral()/area_total*100<<std::endl;
    std::cout<<"h_RES_Other Area = "<<h_RES_Other->Integral()/area_total*100<<std::endl;
    std::cout<<"h_NonRES Area = "<<h_NonRES->Integral()/area_total*100<<std::endl;
    delete mc_total;

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
    plotter->DrawDataStackedMC(h_data, mc_hists, 1.0, "N", "Data (3.33e20 POT)", 0, 0, 1001);

    h_data_stat_only->Draw("SAME E1 X0");

    // ------------------------------------------------------------------------
    // Legend 
    // ------------------------------------------------------------------------
    
    double leg_x_min = 0.55;
    double leg_x_max = 0.90;
    double leg_y_min = 0.57;
    double leg_y_max = 0.85;
    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);

    ApplyStyle_Legend(legend);
    legend->AddEntry(h_data_stat_only, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(h_RES_Delta, "Delta resonance","f");
    legend->AddEntry(h_RES_Other, "Other resonances","f");
    legend->AddEntry(h_NonRES, "Non-Resonant","f");
    legend->Draw();

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");

    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var_name + "_xsec_IntType" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete h_data;
    delete h_data_stat_only;
    delete h_RES_Delta;
    delete h_RES_Other;
    delete h_NonRES;
    delete mc_hists;
    delete legend;
    delete c;
    delete plotter;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_pi0_KE_IntType()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "pi0_KE";
    std::string data_var = var_name + "_xsec";
    std::string hist_name;

    MnvH1D* data = GetMnvH1D(f_data, data_var);
 
    data->GetXaxis()->SetNdivisions(5,4,0);
    data->GetYaxis()->SetNdivisions(5,4,0);
    data->GetYaxis()->CenterTitle();

    MnvH1D* temp = NULL;
    std::vector<MnvH1D*> mc_IntType;
    for (int i = 0; i < nIntType; ++i){
        hist_name = var_name + "_xsec_IntType_" + std::to_string((long long int)i);
        temp = GetMnvH1D(f_mc, hist_name);
        temp->ClearAllErrorBands();
        mc_IntType.push_back(temp);
    }

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    TCanvas* c = new TCanvas("c");

    MnvH1D* h_data = new MnvH1D (*data);
    TH1D* h_data_stat_only = GetBinNormalizedTH1D(h_data, false);
    h_data_stat_only->SetMarkerStyle(20);
    h_data_stat_only->SetMarkerSize(1);
    h_data_stat_only->SetMarkerColor(kBlack);
    h_data_stat_only->SetLineWidth(1);
    h_data_stat_only->SetLineColor(kBlack);
    h_data->GetXaxis()->SetNdivisions(5,4,0);
    h_data->GetYaxis()->SetNdivisions(5,4,0);
    h_data->GetYaxis()->CenterTitle();

    MnvH1D* h_RES_Delta = new MnvH1D (*mc_IntType[0]);
    h_RES_Delta->SetLineWidth(1);
    h_RES_Delta->SetLineColor(17);
    h_RES_Delta->SetFillColor(17);
    h_RES_Delta->SetFillStyle(1001);
    h_RES_Delta->GetXaxis()->SetNdivisions(5,4,0);
    h_RES_Delta->GetYaxis()->SetNdivisions(5,4,0);
    h_RES_Delta->GetYaxis()->CenterTitle();

    MnvH1D* h_RES_Other = new MnvH1D (*mc_IntType[1]);
    h_RES_Other->SetLineWidth(1);
    h_RES_Other->SetLineColor(kGreen+2);
    h_RES_Other->SetFillColor(kGreen+2);
    h_RES_Other->SetFillStyle(1001);
    h_RES_Other->GetXaxis()->SetNdivisions(5,4,0);
    h_RES_Other->GetYaxis()->SetNdivisions(5,4,0);
    h_RES_Other->GetYaxis()->CenterTitle();

    MnvH1D* h_NonRES = new MnvH1D (*mc_IntType[2]);
    h_NonRES->SetLineWidth(1);
    h_NonRES->SetLineColor(kRed+2);
    h_NonRES->SetFillColor(kRed+2);
    h_NonRES->SetFillStyle(1001);
    h_NonRES->GetXaxis()->SetNdivisions(5,4,0);
    h_NonRES->GetYaxis()->SetNdivisions(5,4,0);
    h_NonRES->GetYaxis()->CenterTitle();
    TObjArray* mc_hists = new TObjArray;
    mc_hists->Add(h_NonRES);
    mc_hists->Add(h_RES_Other);
    mc_hists->Add(h_RES_Delta);

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
    plotter->DrawDataStackedMC(h_data, mc_hists, 1.0, "N", "Data (3.33e20 POT)", 0, 0, 1001);

    h_data_stat_only->Draw("SAME E1 X0");

    // ------------------------------------------------------------------------
    // Legend 
    // ------------------------------------------------------------------------
    
    double leg_x_min = 0.55;
    double leg_x_max = 0.90;
    double leg_y_min = 0.57;
    double leg_y_max = 0.85;
    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);

    ApplyStyle_Legend(legend);
    legend->AddEntry(h_data_stat_only, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(h_RES_Delta, "Delta resonance","f");
    legend->AddEntry(h_RES_Other, "Other resonances","f");
    legend->AddEntry(h_NonRES, "Non-Resonant","f");
    legend->Draw();

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");

    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var_name + "_xsec_IntType" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete h_data;
    delete h_data_stat_only;
    delete h_RES_Delta;
    delete h_RES_Other;
    delete h_NonRES;
    delete mc_hists;
    delete legend;
    delete c;
    delete plotter;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_muon_theta_IntType()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "muon_theta";
    std::string data_var = var_name + "_xsec";
    std::string hist_name;

    MnvH1D* data = GetMnvH1D(f_data, data_var);
 
    data->GetXaxis()->SetNdivisions(5,4,0);
    data->GetYaxis()->SetNdivisions(5,4,0);
    data->GetYaxis()->CenterTitle();

    MnvH1D* temp = NULL;
    std::vector<MnvH1D*> mc_IntType;
    for (int i = 0; i < nIntType; ++i){
        hist_name = var_name + "_xsec_IntType_" + std::to_string((long long int)i);
        temp = GetMnvH1D(f_mc, hist_name);
        temp->ClearAllErrorBands();
        mc_IntType.push_back(temp);
    }

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 2.0;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    TCanvas* c = new TCanvas("c");

    MnvH1D* h_data = new MnvH1D (*data);
    TH1D* h_data_stat_only = GetBinNormalizedTH1D(h_data, false);
    h_data_stat_only->SetMarkerStyle(20);
    h_data_stat_only->SetMarkerSize(1);
    h_data_stat_only->SetMarkerColor(kBlack);
    h_data_stat_only->SetLineWidth(1);
    h_data_stat_only->SetLineColor(kBlack);
    h_data->GetXaxis()->SetNdivisions(5,4,0);
    h_data->GetYaxis()->SetNdivisions(5,4,0);
    h_data->GetYaxis()->CenterTitle();

    MnvH1D* h_RES_Delta = new MnvH1D (*mc_IntType[0]);
    h_RES_Delta->SetLineWidth(2);
    h_RES_Delta->SetLineColor(kGray+1);
    h_RES_Delta->SetFillColor(kGray+1);
    h_RES_Delta->SetFillStyle(3002);
    h_RES_Delta->GetXaxis()->SetNdivisions(5,4,0);
    h_RES_Delta->GetYaxis()->SetNdivisions(5,4,0);
    h_RES_Delta->GetYaxis()->CenterTitle();

    MnvH1D* h_RES_Other = new MnvH1D (*mc_IntType[1]);
    h_RES_Other->SetLineWidth(2);
    h_RES_Other->SetLineColor(kGreen+3);
    h_RES_Other->SetFillColor(kGreen+3);
    h_RES_Other->SetFillStyle(3002);
    h_RES_Other->GetXaxis()->SetNdivisions(5,4,0);
    h_RES_Other->GetYaxis()->SetNdivisions(5,4,0);
    h_RES_Other->GetYaxis()->CenterTitle();

    MnvH1D* h_NonRES = new MnvH1D (*mc_IntType[2]);
    h_NonRES->SetLineWidth(2);
    h_NonRES->SetLineColor(kRed+4);
    h_NonRES->SetFillColor(kRed+4);
    h_NonRES->SetFillStyle(3002);
    h_NonRES->GetXaxis()->SetNdivisions(5,4,0);
    h_NonRES->GetYaxis()->SetNdivisions(5,4,0);
    h_NonRES->GetYaxis()->CenterTitle();
    TObjArray* mc_hists = new TObjArray;
    mc_hists->Add(h_NonRES);
    mc_hists->Add(h_RES_Other);
    mc_hists->Add(h_RES_Delta);

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
    plotter->DrawDataStackedMC(h_data, mc_hists, 1.0, "N", "Data (3.33e20 POT)", 0, 0, 3002);

    h_data_stat_only->Draw("SAME E1 X0");

    // ------------------------------------------------------------------------
    // Legend 
    // ------------------------------------------------------------------------
    
    double leg_x_min = 0.55;
    double leg_x_max = 0.90;
    double leg_y_min = 0.62;
    double leg_y_max = 0.90;
    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);

    ApplyStyle_Legend(legend);
    legend->AddEntry(h_data_stat_only, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(h_RES_Delta, "Delta resonance","f");
    legend->AddEntry(h_RES_Other, "Other resonances","f");
    legend->AddEntry(h_NonRES, "Non-Resonant","f");
    legend->Draw();

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");

    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var_name + "_xsec_IntType" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete h_data;
    delete h_data_stat_only;
    delete h_RES_Delta;
    delete h_RES_Other;
    delete h_NonRES;
    delete mc_hists;
    delete legend;
    delete c;
    delete plotter;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_muon_P_IntType()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "muon_P";
    std::string data_var = var_name + "_xsec";
    std::string hist_name;

    MnvH1D* data = GetMnvH1D(f_data, data_var);
 
    data->GetXaxis()->SetNdivisions(5,4,0);
    data->GetYaxis()->SetNdivisions(5,4,0);
    data->GetYaxis()->CenterTitle();

    MnvH1D* temp = NULL;
    std::vector<MnvH1D*> mc_IntType;
    for (int i = 0; i < nIntType; ++i){
        hist_name = var_name + "_xsec_IntType_" + std::to_string((long long int)i);
        temp = GetMnvH1D(f_mc, hist_name);
        temp->ClearAllErrorBands();
        mc_IntType.push_back(temp);
    }

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 0.4;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    TCanvas* c = new TCanvas("c");

    MnvH1D* h_data = new MnvH1D (*data);
    TH1D* h_data_stat_only = GetBinNormalizedTH1D(h_data, false);
    h_data_stat_only->SetMarkerStyle(20);
    h_data_stat_only->SetMarkerSize(1);
    h_data_stat_only->SetMarkerColor(kBlack);
    h_data_stat_only->SetLineWidth(1);
    h_data_stat_only->SetLineColor(kBlack);
    h_data->GetXaxis()->SetNdivisions(5,4,0);
    h_data->GetYaxis()->SetNdivisions(5,4,0);
    h_data->GetYaxis()->CenterTitle();

    MnvH1D* h_RES_Delta = new MnvH1D (*mc_IntType[0]);
    h_RES_Delta->SetLineWidth(2);
    h_RES_Delta->SetLineColor(kGray+1);
    h_RES_Delta->SetFillColor(kGray+1);
    h_RES_Delta->SetFillStyle(3002);
    h_RES_Delta->GetXaxis()->SetNdivisions(5,4,0);
    h_RES_Delta->GetYaxis()->SetNdivisions(5,4,0);
    h_RES_Delta->GetYaxis()->CenterTitle();

    MnvH1D* h_RES_Other = new MnvH1D (*mc_IntType[1]);
    h_RES_Other->SetLineWidth(2);
    h_RES_Other->SetLineColor(kGreen+3);
    h_RES_Other->SetFillColor(kGreen+3);
    h_RES_Other->SetFillStyle(3002);
    h_RES_Other->GetXaxis()->SetNdivisions(5,4,0);
    h_RES_Other->GetYaxis()->SetNdivisions(5,4,0);
    h_RES_Other->GetYaxis()->CenterTitle();

    MnvH1D* h_NonRES = new MnvH1D (*mc_IntType[2]);
    h_NonRES->SetLineWidth(2);
    h_NonRES->SetLineColor(kRed+4);
    h_NonRES->SetFillColor(kRed+4);
    h_NonRES->SetFillStyle(3002);
    h_NonRES->GetXaxis()->SetNdivisions(5,4,0);
    h_NonRES->GetYaxis()->SetNdivisions(5,4,0);
    h_NonRES->GetYaxis()->CenterTitle();

    TObjArray* mc_hists = new TObjArray;
    mc_hists->Add(h_NonRES);
    mc_hists->Add(h_RES_Other);
    mc_hists->Add(h_RES_Delta);

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
    plotter->DrawDataStackedMC(h_data, mc_hists, 1.0, "N", "Data (3.33e20 POT)", 0, 0, 3002);

    h_data_stat_only->Draw("SAME E1 X0");

    // ------------------------------------------------------------------------
    // Legend 
    // ------------------------------------------------------------------------
    
    double leg_x_min = 0.55;
    double leg_x_max = 0.90;
    double leg_y_min = 0.57;
    double leg_y_max = 0.85;
    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);

    ApplyStyle_Legend(legend);
    legend->AddEntry(h_data_stat_only, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(h_RES_Delta, "Delta resonance","f");
    legend->AddEntry(h_RES_Other, "Other resonances","f");
    legend->AddEntry(h_NonRES, "Non-Resonant","f");
    legend->Draw();

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");

    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var_name + "_xsec_IntType" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete h_data;
    delete h_data_stat_only;
    delete h_RES_Delta;
    delete h_RES_Other;
    delete h_NonRES;
    delete mc_hists;
    delete legend;
    delete c;
    delete plotter;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_Enu_FSIType()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "Enu";
    std::string data_var = var_name + "_xsec";
    std::string mc_var = var_name + "_xsec";
    std::string hist_name;

    MnvH1D* data = GetMnvH1D(f_data, data_var);
    MnvH1D* mc = GetMnvH1D(f_mc, mc_var);

    mc->ClearAllErrorBands();

    MnvH1D* temp = NULL;
    std::vector<MnvH1D*> mc_FSIType;
    for (int i = 0; i < nFSIType; ++i){
        hist_name = var_name + "_xsec_FSIType_" + std::to_string((long long int)i);
        temp = GetMnvH1D(f_mc, hist_name);
        temp->ClearAllErrorBands();
        mc_FSIType.push_back(temp);
    }

    std::string norm_label;
    double mc_ratio = GetMCNormalization(norm_label, false, data, mc);

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 2.50;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    TCanvas* c = new TCanvas("c");

    MnvH1D* h_data = new MnvH1D (*data);
    TH1D* h_data_stat_only = GetBinNormalizedTH1D(h_data, false);
    h_data_stat_only->SetMarkerStyle(20);
    h_data_stat_only->SetMarkerSize(1);
    h_data_stat_only->SetMarkerColor(kBlack);
    h_data_stat_only->SetLineWidth(1);
    h_data_stat_only->SetLineColor(kBlack);
    h_data->GetXaxis()->SetNdivisions(8,4,0);
    h_data->GetYaxis()->SetNdivisions(8,4,0);
    h_data->GetYaxis()->CenterTitle();

    MnvH1D* h_NonInteracting = new MnvH1D (*mc_FSIType[0]);
    h_NonInteracting->SetTitle("#pi^{0} Non-interacting");
    h_NonInteracting->SetLineWidth(1);
    h_NonInteracting->SetLineColor(kRed+1);
    h_NonInteracting->SetFillColor(kRed+1);
    h_NonInteracting->SetFillStyle(1001);
    h_NonInteracting->GetXaxis()->SetNdivisions(8,4,0);
    h_NonInteracting->GetYaxis()->SetNdivisions(8,4,0);
    h_NonInteracting->GetYaxis()->CenterTitle();

    MnvH1D* h_Elastic = new MnvH1D (*mc_FSIType[1]);
    h_Elastic->SetTitle("#pi^{0} Elastic");
    h_Elastic->SetLineWidth(1);
    h_Elastic->SetLineColor(kOrange+6);
    h_Elastic->SetFillColor(kOrange+6);
    h_Elastic->SetFillStyle(1001);
    h_Elastic->GetXaxis()->SetNdivisions(8,4,0);
    h_Elastic->GetYaxis()->SetNdivisions(8,4,0);
    h_Elastic->GetYaxis()->CenterTitle();

    MnvH1D* h_Inelastic = new MnvH1D(*mc_FSIType[2]);
    h_Inelastic->SetTitle("#pi^{0} Inelastic");
    h_Inelastic->SetLineWidth(1);
    h_Inelastic->SetLineColor(kGray+1);
    h_Inelastic->SetFillColor(kGray+1);
    h_Inelastic->SetFillStyle(1001);
    h_Inelastic->GetXaxis()->SetNdivisions(8,4,0);
    h_Inelastic->GetYaxis()->SetNdivisions(8,4,0);
    h_Inelastic->GetYaxis()->CenterTitle();

    MnvH1D* h_Cex = new MnvH1D(*mc_FSIType[3]);
    h_Cex->SetTitle("#pi^{#pm} #rightarrow #pi^{0}");
    h_Cex->SetLineWidth(1);
    h_Cex->SetLineColor(kGreen-2);
    h_Cex->SetFillColor(kGreen-2);
    h_Cex->SetFillStyle(1001);
    h_Cex->GetXaxis()->SetNdivisions(8,4,0);
    h_Cex->GetYaxis()->SetNdivisions(8,4,0);
    h_Cex->GetYaxis()->CenterTitle();

    MnvH1D* h_MultiPion = new MnvH1D(*mc_FSIType[4]);
    h_MultiPion->SetTitle("multi-#pi #rightarrow #pi^{0}");
    h_MultiPion->SetLineWidth(1);
    h_MultiPion->SetLineColor(kCyan-8);
    h_MultiPion->SetFillColor(kCyan-8);
    h_MultiPion->SetFillStyle(1001);
    h_MultiPion->GetXaxis()->SetNdivisions(8,4,0);
    h_MultiPion->GetYaxis()->SetNdivisions(8,4,0);
    h_MultiPion->GetYaxis()->CenterTitle();

    MnvH1D* h_Other = new MnvH1D(*mc_FSIType[5]);
    h_Other->SetTitle("Other #pi^{0} production");
    h_Other->SetLineWidth(1);
    h_Other->SetLineColor(kMagenta+1);
    h_Other->SetFillColor(kMagenta+1);
    h_Other->SetFillStyle(1001);
    h_Other->GetXaxis()->SetNdivisions(8,4,0);
    h_Other->GetYaxis()->SetNdivisions(8,4,0);
    h_Other->GetYaxis()->CenterTitle();


    TObjArray* mc_hists = new TObjArray;
    mc_hists->Add(h_NonInteracting);
    mc_hists->Add(h_Elastic);
    mc_hists->Add(h_Inelastic);
    mc_hists->Add(h_Other);
    mc_hists->Add(h_Cex);
    mc_hists->Add(h_MultiPion);
   

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
    //mc_ratio = 1.0; // For pi0_P and pi0_theta from Trung's Paper
    plotter->DrawDataStackedMC(h_data, mc_hists, mc_ratio , "N", "Data (3.33e20 POT)", 0, 0, 1001);

    h_data_stat_only->Draw("SAME E1 X0");

    // ------------------------------------------------------------------------
    // Legend 
    // ------------------------------------------------------------------------

    double leg_x_min = 0.60;
    double leg_x_max = 0.90;
    double leg_y_min = 0.52;
    double leg_y_max = 0.90;

    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    ApplyStyle_Legend(legend);
    legend->SetTextSize(0.04);
    legend->AddEntry(h_data_stat_only, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(h_MultiPion, "multi-#pi #rightarrow #pi^{0}","f");
    legend->AddEntry(h_Cex, "#pi^{#pm} #rightarrow #pi^{0}","f");
    legend->AddEntry(h_Other, "zero-#pi #rightarrow #pi^{0}","f");
    legend->AddEntry(h_Inelastic, "#pi^{0} Inelastic","f");
    legend->AddEntry(h_Elastic, "#pi^{0} Elastic","f");
    legend->AddEntry(h_NonInteracting, "#pi^{0} Non-interacting","f");
    legend->Draw();

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"Area Normalized");
    //text.DrawLatex(0.30,0.87,"POT Normalized"); // For pi0_P and pi0_theta from Trung's Paper

    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var_name + "_xsec_FSIType" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete h_data;
    delete h_data_stat_only;
    delete h_NonInteracting;
    delete h_Elastic;
    delete h_Inelastic;
    delete h_Cex;
    delete h_MultiPion;
    delete h_Other;
    delete mc_hists;
    delete legend;
    delete c;
    delete plotter;
    delete f_mc;
    delete f_data;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_QSq_FSIType()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "QSq";
    std::string data_var = var_name + "_xsec";
    std::string mc_var = var_name + "_xsec";
    std::string hist_name;

    MnvH1D* data = GetMnvH1D(f_data, data_var);
    MnvH1D* mc = GetMnvH1D(f_mc, mc_var);

    mc->ClearAllErrorBands();

    MnvH1D* temp = NULL;
    std::vector<MnvH1D*> mc_FSIType;
    for (int i = 0; i < nFSIType; ++i){
        hist_name = var_name + "_xsec_FSIType_" + std::to_string((long long int)i);
        temp = GetMnvH1D(f_mc, hist_name);
        temp->ClearAllErrorBands();
        mc_FSIType.push_back(temp);
    }

  
    std::string norm_label;
    double mc_ratio = GetMCNormalization(norm_label, false, data, mc);

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    TCanvas* c = new TCanvas("c");

    MnvH1D* h_data = new MnvH1D (*data);
    TH1D* h_data_stat_only = GetBinNormalizedTH1D(h_data, false);
    h_data_stat_only->SetMarkerStyle(20);
    h_data_stat_only->SetMarkerSize(1);
    h_data_stat_only->SetMarkerColor(kBlack);
    h_data_stat_only->SetLineWidth(1);
    h_data_stat_only->SetLineColor(kBlack);
    h_data->GetXaxis()->SetNdivisions(5,4,0);
    h_data->GetYaxis()->SetNdivisions(5,4,0);
    h_data->GetYaxis()->CenterTitle();

    MnvH1D* h_NonInteracting = new MnvH1D (*mc_FSIType[0]);
    h_NonInteracting->SetTitle("#pi^{0} Non-interacting");
    h_NonInteracting->SetLineWidth(1);
    h_NonInteracting->SetLineColor(kRed+1);
    h_NonInteracting->SetFillColor(kRed+1);
    h_NonInteracting->SetFillStyle(1001);
    h_NonInteracting->GetXaxis()->SetNdivisions(5,4,0);
    h_NonInteracting->GetYaxis()->SetNdivisions(5,4,0);
    h_NonInteracting->GetYaxis()->CenterTitle();

    MnvH1D* h_Elastic = new MnvH1D (*mc_FSIType[1]);
    h_Elastic->SetTitle("#pi^{0} Elastic");
    h_Elastic->SetLineWidth(1);
    h_Elastic->SetLineColor(kOrange+6);
    h_Elastic->SetFillColor(kOrange+6);
    h_Elastic->SetFillStyle(1001);
    h_Elastic->GetXaxis()->SetNdivisions(5,4,0);
    h_Elastic->GetYaxis()->SetNdivisions(5,4,0);
    h_Elastic->GetYaxis()->CenterTitle();

    MnvH1D* h_Inelastic = new MnvH1D(*mc_FSIType[2]);
    h_Inelastic->SetTitle("#pi^{0} Inelastic");
    h_Inelastic->SetLineWidth(1);
    h_Inelastic->SetLineColor(kGray+1);
    h_Inelastic->SetFillColor(kGray+1);
    h_Inelastic->SetFillStyle(1001);
    h_Inelastic->GetXaxis()->SetNdivisions(5,4,0);
    h_Inelastic->GetYaxis()->SetNdivisions(5,4,0);
    h_Inelastic->GetYaxis()->CenterTitle();

    MnvH1D* h_Cex = new MnvH1D(*mc_FSIType[3]);
    h_Cex->SetTitle("#pi^{#pm} #rightarrow #pi^{0}");
    h_Cex->SetLineWidth(1);
    h_Cex->SetLineColor(kGreen-2);
    h_Cex->SetFillColor(kGreen-2);
    h_Cex->SetFillStyle(1001);
    h_Cex->GetXaxis()->SetNdivisions(5,4,0);
    h_Cex->GetYaxis()->SetNdivisions(5,4,0);
    h_Cex->GetYaxis()->CenterTitle();

    MnvH1D* h_MultiPion = new MnvH1D(*mc_FSIType[4]);
    h_MultiPion->SetTitle("multi-#pi #rightarrow #pi^{0}");
    h_MultiPion->SetLineWidth(1);
    h_MultiPion->SetLineColor(kCyan-8);
    h_MultiPion->SetFillColor(kCyan-8);
    h_MultiPion->SetFillStyle(1001);
    h_MultiPion->GetXaxis()->SetNdivisions(5,4,0);
    h_MultiPion->GetYaxis()->SetNdivisions(5,4,0);
    h_MultiPion->GetYaxis()->CenterTitle();

    MnvH1D* h_Other = new MnvH1D(*mc_FSIType[5]);
    h_Other->SetTitle("Other #pi^{0} production");
    h_Other->SetLineWidth(1);
    h_Other->SetLineColor(kMagenta+1);
    h_Other->SetFillColor(kMagenta+1);
    h_Other->SetFillStyle(1001);
    h_Other->GetXaxis()->SetNdivisions(5,4,0);
    h_Other->GetYaxis()->SetNdivisions(5,4,0);
    h_Other->GetYaxis()->CenterTitle();



    TObjArray* mc_hists = new TObjArray;
    mc_hists->Add(h_NonInteracting);
    mc_hists->Add(h_Elastic);
    mc_hists->Add(h_Inelastic);
    mc_hists->Add(h_Other);
    mc_hists->Add(h_Cex);
    mc_hists->Add(h_MultiPion);
   

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
    //mc_ratio = 1.0; // For pi0_P and pi0_theta from Trung's Paper
    plotter->DrawDataStackedMC(h_data, mc_hists, mc_ratio , "N", "Data (3.33e20 POT)", 0, 0, 1001);

    h_data_stat_only->Draw("SAME E1 X0");

    // ------------------------------------------------------------------------
    // Legend 
    // ------------------------------------------------------------------------

    double leg_x_min = 0.55;
    double leg_x_max = 0.90;
    double leg_y_min = 0.36;
    double leg_y_max = 0.85;

    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    ApplyStyle_Legend(legend);
    legend->SetTextSize(0.04);
    legend->AddEntry(h_data_stat_only, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(h_MultiPion, "multi-#pi #rightarrow #pi^{0}","f");
    legend->AddEntry(h_Cex, "#pi^{#pm} #rightarrow #pi^{0}","f");
    legend->AddEntry(h_Other, "zero-#pi #rightarrow #pi^{0}","f");
    legend->AddEntry(h_Inelastic, "#pi^{0} Inelastic","f");
    legend->AddEntry(h_Elastic, "#pi^{0} Elastic","f");
    legend->AddEntry(h_NonInteracting, "#pi^{0} Non-interacting","f");
    legend->Draw();

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"Area Normalized");
    //text.DrawLatex(0.30,0.87,"POT Normalized"); // For pi0_P and pi0_theta from Trung's Paper

    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var_name + "_xsec_FSIType" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete h_data;
    delete h_data_stat_only;
    delete h_NonInteracting;
    delete h_Elastic;
    delete h_Inelastic;
    delete h_Cex;
    delete h_MultiPion;
    delete h_Other;
    delete mc_hists;
    delete legend;
    delete c;
    delete plotter;
    delete f_mc;
    delete f_data;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_pi0_theta_FSIType()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "pi0_theta";
    std::string data_var = var_name + "_xsec";
    std::string mc_var = var_name + "_xsec";
    std::string hist_name;

    MnvH1D* data = GetMnvH1D(f_data, data_var);
    MnvH1D* mc = GetMnvH1D(f_mc, mc_var);
    mc->ClearAllErrorBands();

    MnvH1D* temp = NULL;
    std::vector<MnvH1D*> mc_FSIType;
    for (int i = 0; i < nFSIType; ++i){
        hist_name = var_name + "_xsec_FSIType_" + std::to_string((long long int)i);
        temp = GetMnvH1D(f_mc, hist_name);
        temp->ClearAllErrorBands();
        mc_FSIType.push_back(temp);
    }
  
    std::string norm_label;
    double mc_ratio = GetMCNormalization(norm_label, false, data, mc);

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    TCanvas* c = new TCanvas("c");

    MnvH1D* h_data = new MnvH1D (*data);
    TH1D* h_data_stat_only = GetBinNormalizedTH1D(h_data, false);
    h_data_stat_only->SetMarkerStyle(20);
    h_data_stat_only->SetMarkerSize(1);
    h_data_stat_only->SetMarkerColor(kBlack);
    h_data_stat_only->SetLineWidth(1);
    h_data_stat_only->SetLineColor(kBlack);
    h_data->GetXaxis()->SetNdivisions(5,4,0);
    h_data->GetYaxis()->SetNdivisions(5,4,0);
    h_data->GetYaxis()->CenterTitle();

    MnvH1D* h_NonInteracting = new MnvH1D (*mc_FSIType[0]);
    h_NonInteracting->SetTitle("#pi^{0} Non-interacting");
    h_NonInteracting->SetLineWidth(1);
    h_NonInteracting->SetLineColor(kRed+1);
    h_NonInteracting->SetFillColor(kRed+1);
    h_NonInteracting->SetFillStyle(1001);
    h_NonInteracting->GetXaxis()->SetNdivisions(5,4,0);
    h_NonInteracting->GetYaxis()->SetNdivisions(5,4,0);
    h_NonInteracting->GetYaxis()->CenterTitle();

    MnvH1D* h_Elastic = new MnvH1D (*mc_FSIType[1]);
    h_Elastic->SetTitle("#pi^{0} Elastic");
    h_Elastic->SetLineWidth(1);
    h_Elastic->SetLineColor(kOrange+6);
    h_Elastic->SetFillColor(kOrange+6);
    h_Elastic->SetFillStyle(1001);
    h_Elastic->GetXaxis()->SetNdivisions(5,4,0);
    h_Elastic->GetYaxis()->SetNdivisions(5,4,0);
    h_Elastic->GetYaxis()->CenterTitle();

    MnvH1D* h_Inelastic = new MnvH1D(*mc_FSIType[2]);
    h_Inelastic->SetTitle("#pi^{0} Inelastic");
    h_Inelastic->SetLineWidth(1);
    h_Inelastic->SetLineColor(kGray+1);
    h_Inelastic->SetFillColor(kGray+1);
    h_Inelastic->SetFillStyle(1001);
    h_Inelastic->GetXaxis()->SetNdivisions(5,4,0);
    h_Inelastic->GetYaxis()->SetNdivisions(5,4,0);
    h_Inelastic->GetYaxis()->CenterTitle();

    MnvH1D* h_Cex = new MnvH1D(*mc_FSIType[3]);
    h_Cex->SetTitle("#pi^{#pm} #rightarrow #pi^{0}");
    h_Cex->SetLineWidth(1);
    h_Cex->SetLineColor(kGreen-2);
    h_Cex->SetFillColor(kGreen-2);
    h_Cex->SetFillStyle(1001);
    h_Cex->GetXaxis()->SetNdivisions(5,4,0);
    h_Cex->GetYaxis()->SetNdivisions(5,4,0);
    h_Cex->GetYaxis()->CenterTitle();

    MnvH1D* h_MultiPion = new MnvH1D(*mc_FSIType[4]);
    h_MultiPion->SetTitle("multi-#pi #rightarrow #pi^{0}");
    h_MultiPion->SetLineWidth(1);
    h_MultiPion->SetLineColor(kCyan-8);
    h_MultiPion->SetFillColor(kCyan-8);
    h_MultiPion->SetFillStyle(1001);
    h_MultiPion->GetXaxis()->SetNdivisions(5,4,0);
    h_MultiPion->GetYaxis()->SetNdivisions(5,4,0);
    h_MultiPion->GetYaxis()->CenterTitle();

    MnvH1D* h_Other = new MnvH1D(*mc_FSIType[5]);
    h_Other->SetTitle("Other #pi^{0} production");
    h_Other->SetLineWidth(1);
    h_Other->SetLineColor(kMagenta+1);
    h_Other->SetFillColor(kMagenta+1);
    h_Other->SetFillStyle(1001);
    h_Other->GetXaxis()->SetNdivisions(5,4,0);
    h_Other->GetYaxis()->SetNdivisions(5,4,0);
    h_Other->GetYaxis()->CenterTitle();

    TObjArray* mc_hists = new TObjArray;
    mc_hists->Add(h_NonInteracting);
    mc_hists->Add(h_Elastic);
    mc_hists->Add(h_Inelastic);
    mc_hists->Add(h_Other);
    mc_hists->Add(h_Cex);
    mc_hists->Add(h_MultiPion);
   

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
    //mc_ratio = 1.0; // For pi0_P and pi0_theta from Trung's Paper
    plotter->DrawDataStackedMC(h_data, mc_hists, mc_ratio , "N", "Data (3.33e20 POT)", 0, 0, 1001);

    h_data_stat_only->Draw("SAME E1 X0");

    // ------------------------------------------------------------------------
    // Legend 
    // ------------------------------------------------------------------------

    double leg_x_min = 0.62;
    double leg_x_max = 0.90;
    double leg_y_min = 0.42;
    double leg_y_max = 0.90;

    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    ApplyStyle_Legend(legend);
    legend->SetTextSize(0.04);
    legend->AddEntry(h_data_stat_only, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(h_MultiPion, "multi-#pi #rightarrow #pi^{0}","f");
    legend->AddEntry(h_Cex, "#pi^{#pm} #rightarrow #pi^{0}","f");
    legend->AddEntry(h_Other, "zero-#pi #rightarrow #pi^{0}","f");
    legend->AddEntry(h_Inelastic, "#pi^{0} Inelastic","f");
    legend->AddEntry(h_Elastic, "#pi^{0} Elastic","f");
    legend->AddEntry(h_NonInteracting, "#pi^{0} Non-interacting","f");
    legend->Draw();

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"Area Normalized");
    //text.DrawLatex(0.30,0.87,"POT Normalized"); // For pi0_P and pi0_theta from Trung's Paper

    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var_name + "_xsec_FSIType" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete h_data;
    delete h_data_stat_only;
    delete h_NonInteracting;
    delete h_Elastic;
    delete h_Inelastic;
    delete h_Cex;
    delete h_MultiPion;
    delete h_Other;
    delete mc_hists;
    delete legend;
    delete c;
    delete plotter;
    delete f_mc;
    delete f_data;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_pi0_KE_FSIType()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "pi0_KE";
    std::string data_var = var_name + "_xsec";
    std::string mc_var = var_name + "_xsec";
    std::string hist_name;

    MnvH1D* data = GetMnvH1D(f_data, data_var);
    MnvH1D* mc = GetMnvH1D(f_mc, mc_var);

    mc->ClearAllErrorBands();

    MnvH1D* temp = NULL;
    std::vector<MnvH1D*> mc_FSIType;
    for (int i = 0; i < nFSIType; ++i){
        hist_name = var_name + "_xsec_FSIType_" + std::to_string((long long int)i);
        temp = GetMnvH1D(f_mc, hist_name);
        temp->ClearAllErrorBands();
        mc_FSIType.push_back(temp);
    }

  
    std::string norm_label;
    double mc_ratio = GetMCNormalization(norm_label, false, data, mc);

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    TCanvas* c = new TCanvas("c");

    MnvH1D* h_data = new MnvH1D (*data);
    TH1D* h_data_stat_only = GetBinNormalizedTH1D(h_data, false);
    h_data_stat_only->SetMarkerStyle(20);
    h_data_stat_only->SetMarkerSize(1);
    h_data_stat_only->SetMarkerColor(kBlack);
    h_data_stat_only->SetLineWidth(1);
    h_data_stat_only->SetLineColor(kBlack);
    h_data->GetXaxis()->SetNdivisions(5,4,0);
    h_data->GetYaxis()->SetNdivisions(5,4,0);
    h_data->GetYaxis()->CenterTitle();

    MnvH1D* h_NonInteracting = new MnvH1D (*mc_FSIType[0]);
    h_NonInteracting->SetTitle("#pi^{0} Non-interacting");
    h_NonInteracting->SetLineWidth(2);
    h_NonInteracting->SetLineColor(kRed+1);
    h_NonInteracting->SetFillColor(kRed+1);
    h_NonInteracting->SetFillStyle(3001);
    h_NonInteracting->GetXaxis()->SetNdivisions(5,4,0);
    h_NonInteracting->GetYaxis()->SetNdivisions(5,4,0);
    h_NonInteracting->GetYaxis()->CenterTitle();

    MnvH1D* h_Elastic = new MnvH1D (*mc_FSIType[1]);
    h_Elastic->SetTitle("#pi^{0} Elastic");
    h_Elastic->SetLineWidth(2);
    h_Elastic->SetLineColor(kOrange-9);
    h_Elastic->SetFillColor(kOrange-9);
    h_Elastic->SetFillStyle(3001);
    h_Elastic->GetXaxis()->SetNdivisions(5,4,0);
    h_Elastic->GetYaxis()->SetNdivisions(5,4,0);
    h_Elastic->GetYaxis()->CenterTitle();

    MnvH1D* h_Inelastic = new MnvH1D(*mc_FSIType[2]);
    h_Inelastic->SetTitle("#pi^{0} Inelastic");
    h_Inelastic->SetLineWidth(2);
    h_Inelastic->SetLineColor(kGray+1);
    h_Inelastic->SetFillColor(kGray+1);
    h_Inelastic->SetFillStyle(3001);
    h_Inelastic->GetXaxis()->SetNdivisions(5,4,0);
    h_Inelastic->GetYaxis()->SetNdivisions(5,4,0);
    h_Inelastic->GetYaxis()->CenterTitle();

    MnvH1D* h_Cex = new MnvH1D(*mc_FSIType[3]);
    h_Cex->SetTitle("#pi^{#pm} #rightarrow #pi^{0}");
    h_Cex->SetLineWidth(2);
    h_Cex->SetLineColor(kGreen-2);
    h_Cex->SetFillColor(kGreen-2);
    h_Cex->SetFillStyle(3001);
    h_Cex->GetXaxis()->SetNdivisions(5,4,0);
    h_Cex->GetYaxis()->SetNdivisions(5,4,0);
    h_Cex->GetYaxis()->CenterTitle();

    MnvH1D* h_MultiPion = new MnvH1D(*mc_FSIType[4]);
    h_MultiPion->SetTitle("multi-#pi #rightarrow #pi^{0}");
    h_MultiPion->SetLineWidth(2);
    h_MultiPion->SetLineColor(kCyan-8);
    h_MultiPion->SetFillColor(kCyan-8);
    h_MultiPion->SetFillStyle(3001);
    h_MultiPion->GetXaxis()->SetNdivisions(5,4,0);
    h_MultiPion->GetYaxis()->SetNdivisions(5,4,0);
    h_MultiPion->GetYaxis()->CenterTitle();

    MnvH1D* h_Other = new MnvH1D(*mc_FSIType[5]);
    h_Other->SetTitle("Other #pi^{0} production");
    h_Other->SetLineWidth(2);
    h_Other->SetLineColor(kMagenta+3);
    h_Other->SetFillColor(kMagenta+3);
    h_Other->SetFillStyle(3001);
    h_Other->GetXaxis()->SetNdivisions(5,4,0);
    h_Other->GetYaxis()->SetNdivisions(5,4,0);
    h_Other->GetYaxis()->CenterTitle();


    TObjArray* mc_hists = new TObjArray;
    mc_hists->Add(h_NonInteracting);
    mc_hists->Add(h_Elastic);
    mc_hists->Add(h_Inelastic);
    mc_hists->Add(h_Other);
    mc_hists->Add(h_Cex);
    mc_hists->Add(h_MultiPion);
   

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
    //mc_ratio = 1.0; // For pi0_P and pi0_theta from Trung's Paper
    plotter->DrawDataStackedMC(h_data, mc_hists, mc_ratio , "N", "Data (3.33e20 POT)", 0, 0, 3001);

    h_data_stat_only->Draw("SAME E1 X0");

    // ------------------------------------------------------------------------
    // Legend 
    // ------------------------------------------------------------------------

    double leg_x_min = 0.55;
    double leg_x_max = 0.90;
    double leg_y_min = 0.36;
    double leg_y_max = 0.85;

    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    ApplyStyle_Legend(legend);
    legend->SetTextSize(0.04);
    legend->AddEntry(h_data_stat_only, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(h_MultiPion, "multi-#pi #rightarrow #pi^{0}","f");
    legend->AddEntry(h_Cex, "#pi^{#pm} #rightarrow #pi^{0}","f");
    legend->AddEntry(h_Other, "zero-#pi #rightarrow #pi^{0}","f");
    legend->AddEntry(h_Inelastic, "#pi^{0} Inelastic","f");
    legend->AddEntry(h_Elastic, "#pi^{0} Elastic","f");
    legend->AddEntry(h_NonInteracting, "#pi^{0} Non-interacting","f");
    legend->Draw();

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"Area Normalized");
    //text.DrawLatex(0.30,0.87,"POT Normalized"); // For pi0_P and pi0_theta from Trung's Paper

    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var_name + "_xsec_FSIType" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete h_data;
    delete h_data_stat_only;
    delete h_NonInteracting;
    delete h_Elastic;
    delete h_Inelastic;
    delete h_Cex;
    delete h_MultiPion;
    delete h_Other;
    delete mc_hists;
    delete legend;
    delete c;
    delete plotter;
    delete f_mc;
    delete f_data;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_muon_theta_FSIType()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "muon_theta";
    std::string data_var = var_name + "_xsec";
    std::string mc_var = var_name + "_xsec";
    std::string hist_name;

    MnvH1D* data = GetMnvH1D(f_data, data_var);
    MnvH1D* mc = GetMnvH1D(f_mc, mc_var);

    mc->ClearAllErrorBands();

    MnvH1D* temp = NULL;
    std::vector<MnvH1D*> mc_FSIType;
    for (int i = 0; i < nFSIType; ++i){
        hist_name = var_name + "_xsec_FSIType_" + std::to_string((long long int)i);
        temp = GetMnvH1D(f_mc, hist_name);
        temp->ClearAllErrorBands();
        mc_FSIType.push_back(temp);
    }

    std::string norm_label;
    double mc_ratio = GetMCNormalization(norm_label, false, data, mc);

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 2.50;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    TCanvas* c = new TCanvas("c");

    MnvH1D* h_data = new MnvH1D (*data);
    TH1D* h_data_stat_only = GetBinNormalizedTH1D(h_data, false);
    h_data_stat_only->SetMarkerStyle(20);
    h_data_stat_only->SetMarkerSize(1);
    h_data_stat_only->SetMarkerColor(kBlack);
    h_data_stat_only->SetLineWidth(1);
    h_data_stat_only->SetLineColor(kBlack);
    h_data->GetXaxis()->SetNdivisions(8,4,0);
    h_data->GetYaxis()->SetNdivisions(8,4,0);
    h_data->GetYaxis()->CenterTitle();

    MnvH1D* h_NonInteracting = new MnvH1D (*mc_FSIType[0]);
    h_NonInteracting->SetTitle("#pi^{0} Non-interacting");
    h_NonInteracting->SetLineWidth(1);
    h_NonInteracting->SetLineColor(kRed+1);
    h_NonInteracting->SetFillColor(kRed+1);
    h_NonInteracting->SetFillStyle(1001);
    h_NonInteracting->GetXaxis()->SetNdivisions(8,4,0);
    h_NonInteracting->GetYaxis()->SetNdivisions(8,4,0);
    h_NonInteracting->GetYaxis()->CenterTitle();

    MnvH1D* h_Elastic = new MnvH1D (*mc_FSIType[1]);
    h_Elastic->SetTitle("#pi^{0} Elastic");
    h_Elastic->SetLineWidth(1);
    h_Elastic->SetLineColor(kOrange+6);
    h_Elastic->SetFillColor(kOrange+6);
    h_Elastic->SetFillStyle(1001);
    h_Elastic->GetXaxis()->SetNdivisions(8,4,0);
    h_Elastic->GetYaxis()->SetNdivisions(8,4,0);
    h_Elastic->GetYaxis()->CenterTitle();

    MnvH1D* h_Inelastic = new MnvH1D(*mc_FSIType[2]);
    h_Inelastic->SetTitle("#pi^{0} Inelastic");
    h_Inelastic->SetLineWidth(1);
    h_Inelastic->SetLineColor(kGray+1);
    h_Inelastic->SetFillColor(kGray+1);
    h_Inelastic->SetFillStyle(1001);
    h_Inelastic->GetXaxis()->SetNdivisions(8,4,0);
    h_Inelastic->GetYaxis()->SetNdivisions(8,4,0);
    h_Inelastic->GetYaxis()->CenterTitle();

    MnvH1D* h_Cex = new MnvH1D(*mc_FSIType[3]);
    h_Cex->SetTitle("#pi^{#pm} #rightarrow #pi^{0}");
    h_Cex->SetLineWidth(1);
    h_Cex->SetLineColor(kGreen-2);
    h_Cex->SetFillColor(kGreen-2);
    h_Cex->SetFillStyle(1001);
    h_Cex->GetXaxis()->SetNdivisions(8,4,0);
    h_Cex->GetYaxis()->SetNdivisions(8,4,0);
    h_Cex->GetYaxis()->CenterTitle();

    MnvH1D* h_MultiPion = new MnvH1D(*mc_FSIType[4]);
    h_MultiPion->SetTitle("multi-#pi #rightarrow #pi^{0}");
    h_MultiPion->SetLineWidth(1);
    h_MultiPion->SetLineColor(kCyan-8);
    h_MultiPion->SetFillColor(kCyan-8);
    h_MultiPion->SetFillStyle(1001);
    h_MultiPion->GetXaxis()->SetNdivisions(8,4,0);
    h_MultiPion->GetYaxis()->SetNdivisions(8,4,0);
    h_MultiPion->GetYaxis()->CenterTitle();

    MnvH1D* h_Other = new MnvH1D(*mc_FSIType[5]);
    h_Other->SetTitle("Other #pi^{0} production");
    h_Other->SetLineWidth(1);
    h_Other->SetLineColor(kMagenta+1);
    h_Other->SetFillColor(kMagenta+1);
    h_Other->SetFillStyle(1001);
    h_Other->GetXaxis()->SetNdivisions(8,4,0);
    h_Other->GetYaxis()->SetNdivisions(8,4,0);
    h_Other->GetYaxis()->CenterTitle();


    TObjArray* mc_hists = new TObjArray;
    mc_hists->Add(h_NonInteracting);
    mc_hists->Add(h_Elastic);
    mc_hists->Add(h_Inelastic);
    mc_hists->Add(h_Other);
    mc_hists->Add(h_Cex);
    mc_hists->Add(h_MultiPion);
   

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
    //mc_ratio = 1.0; // For pi0_P and pi0_theta from Trung's Paper
    plotter->DrawDataStackedMC(h_data, mc_hists, mc_ratio , "N", "Data (3.33e20 POT)", 0, 0, 1001);

    h_data_stat_only->Draw("SAME E1 X0");

    // ------------------------------------------------------------------------
    // Legend 
    // ------------------------------------------------------------------------

    double leg_x_min = 0.62;
    double leg_x_max = 0.90;
    double leg_y_min = 0.42;
    double leg_y_max = 0.90;

    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    ApplyStyle_Legend(legend);
    legend->SetTextSize(0.04);
    legend->AddEntry(h_data_stat_only, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(h_MultiPion, "multi-#pi #rightarrow #pi^{0}","f");
    legend->AddEntry(h_Cex, "#pi^{#pm} #rightarrow #pi^{0}","f");
    legend->AddEntry(h_Other, "zero-#pi #rightarrow #pi^{0}","f");
    legend->AddEntry(h_Inelastic, "#pi^{0} Inelastic","f");
    legend->AddEntry(h_Elastic, "#pi^{0} Elastic","f");
    legend->AddEntry(h_NonInteracting, "#pi^{0} Non-interacting","f");
    legend->Draw();

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"Area Normalized");
    //text.DrawLatex(0.30,0.87,"POT Normalized"); // For pi0_P and pi0_theta from Trung's Paper

    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var_name + "_xsec_FSIType" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete h_data;
    delete h_data_stat_only;
    delete h_NonInteracting;
    delete h_Elastic;
    delete h_Inelastic;
    delete h_Cex;
    delete h_MultiPion;
    delete h_Other;
    delete mc_hists;
    delete legend;
    delete c;
    delete plotter;
    delete f_mc;
    delete f_data;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_muon_P_FSIType()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "muon_P";
    std::string data_var = var_name + "_xsec";
    std::string mc_var = var_name + "_xsec";
    std::string hist_name;

    MnvH1D* data = GetMnvH1D(f_data, data_var);
    MnvH1D* mc = GetMnvH1D(f_mc, mc_var);

    mc->ClearAllErrorBands();

    MnvH1D* temp = NULL;
    std::vector<MnvH1D*> mc_FSIType;
    for (int i = 0; i < nFSIType; ++i){
        hist_name = var_name + "_xsec_FSIType_" + std::to_string((long long int)i);
        temp = GetMnvH1D(f_mc, hist_name);
        temp->ClearAllErrorBands();
        mc_FSIType.push_back(temp);
    }

  
    std::string norm_label;
    double mc_ratio = GetMCNormalization(norm_label, false, data, mc);

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 0.8;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    TCanvas* c = new TCanvas("c");

    MnvH1D* h_data = new MnvH1D (*data);
    TH1D* h_data_stat_only = GetBinNormalizedTH1D(h_data, false);
    h_data_stat_only->SetMarkerStyle(20);
    h_data_stat_only->SetMarkerSize(1);
    h_data_stat_only->SetMarkerColor(kBlack);
    h_data_stat_only->SetLineWidth(1);
    h_data_stat_only->SetLineColor(kBlack);
    h_data->GetXaxis()->SetNdivisions(8,4,0);
    h_data->GetYaxis()->SetNdivisions(8,4,0);
    h_data->GetYaxis()->CenterTitle();

    MnvH1D* h_NonInteracting = new MnvH1D (*mc_FSIType[0]);
    h_NonInteracting->SetTitle("#pi^{0} Non-interacting");
    h_NonInteracting->SetLineWidth(1);
    h_NonInteracting->SetLineColor(kRed+1);
    h_NonInteracting->SetFillColor(kRed+1);
    h_NonInteracting->SetFillStyle(1001);
    h_NonInteracting->GetXaxis()->SetNdivisions(8,4,0);
    h_NonInteracting->GetYaxis()->SetNdivisions(8,4,0);
    h_NonInteracting->GetYaxis()->CenterTitle();

    MnvH1D* h_Elastic = new MnvH1D (*mc_FSIType[1]);
    h_Elastic->SetTitle("#pi^{0} Elastic");
    h_Elastic->SetLineWidth(1);
    h_Elastic->SetLineColor(kOrange+6);
    h_Elastic->SetFillColor(kOrange+6);
    h_Elastic->SetFillStyle(1001);
    h_Elastic->GetXaxis()->SetNdivisions(8,4,0);
    h_Elastic->GetYaxis()->SetNdivisions(8,4,0);
    h_Elastic->GetYaxis()->CenterTitle();

    MnvH1D* h_Inelastic = new MnvH1D(*mc_FSIType[2]);
    h_Inelastic->SetTitle("#pi^{0} Inelastic");
    h_Inelastic->SetLineWidth(1);
    h_Inelastic->SetLineColor(kGray+1);
    h_Inelastic->SetFillColor(kGray+1);
    h_Inelastic->SetFillStyle(1001);
    h_Inelastic->GetXaxis()->SetNdivisions(8,4,0);
    h_Inelastic->GetYaxis()->SetNdivisions(8,4,0);
    h_Inelastic->GetYaxis()->CenterTitle();

    MnvH1D* h_Cex = new MnvH1D(*mc_FSIType[3]);
    h_Cex->SetTitle("#pi^{#pm} #rightarrow #pi^{0}");
    h_Cex->SetLineWidth(1);
    h_Cex->SetLineColor(kGreen-2);
    h_Cex->SetFillColor(kGreen-2);
    h_Cex->SetFillStyle(1001);
    h_Cex->GetXaxis()->SetNdivisions(8,4,0);
    h_Cex->GetYaxis()->SetNdivisions(8,4,0);
    h_Cex->GetYaxis()->CenterTitle();

    MnvH1D* h_MultiPion = new MnvH1D(*mc_FSIType[4]);
    h_MultiPion->SetTitle("multi-#pi #rightarrow #pi^{0}");
    h_MultiPion->SetLineWidth(1);
    h_MultiPion->SetLineColor(kCyan-8);
    h_MultiPion->SetFillColor(kCyan-8);
    h_MultiPion->SetFillStyle(1001);
    h_MultiPion->GetXaxis()->SetNdivisions(8,4,0);
    h_MultiPion->GetYaxis()->SetNdivisions(8,4,0);
    h_MultiPion->GetYaxis()->CenterTitle();

    MnvH1D* h_Other = new MnvH1D(*mc_FSIType[5]);
    h_Other->SetTitle("Other #pi^{0} production");
    h_Other->SetLineWidth(1);
    h_Other->SetLineColor(kMagenta+1);
    h_Other->SetFillColor(kMagenta+1);
    h_Other->SetFillStyle(1001);
    h_Other->GetXaxis()->SetNdivisions(8,4,0);
    h_Other->GetYaxis()->SetNdivisions(8,4,0);
    h_Other->GetYaxis()->CenterTitle();

    TObjArray* mc_hists = new TObjArray;
    mc_hists->Add(h_NonInteracting);
    mc_hists->Add(h_Elastic);
    mc_hists->Add(h_Inelastic);
    mc_hists->Add(h_Other);
    mc_hists->Add(h_Cex);
    mc_hists->Add(h_MultiPion);
   

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
    //mc_ratio = 1.0; // For pi0_P and pi0_theta from Trung's Paper
    plotter->DrawDataStackedMC(h_data, mc_hists, mc_ratio , "N", "Data (3.33e20 POT)", 0, 0, 1001);

    h_data_stat_only->Draw("SAME E1 X0");

    // ------------------------------------------------------------------------
    // Legend 
    // ------------------------------------------------------------------------

    double leg_x_min = 0.55;
    double leg_x_max = 0.90;
    double leg_y_min = 0.36;
    double leg_y_max = 0.85;

    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    ApplyStyle_Legend(legend);
    legend->SetTextSize(0.04);
    legend->AddEntry(h_data_stat_only, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(h_MultiPion, "multi-#pi #rightarrow #pi^{0}","f");
    legend->AddEntry(h_Cex, "#pi^{#pm} #rightarrow #pi^{0}","f");
    legend->AddEntry(h_Other, "zero-#pi #rightarrow #pi^{0}","f");
    legend->AddEntry(h_Inelastic, "#pi^{0} Inelastic","f");
    legend->AddEntry(h_Elastic, "#pi^{0} Elastic","f");
    legend->AddEntry(h_NonInteracting, "#pi^{0} Non-interacting","f");
    legend->Draw();

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"Area Normalized");
    //text.DrawLatex(0.30,0.87,"POT Normalized"); // For pi0_P and pi0_theta from Trung's Paper

    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var_name + "_xsec_FSIType" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete h_data;
    delete h_data_stat_only;
    delete h_NonInteracting;
    delete h_Elastic;
    delete h_Inelastic;
    delete h_Cex;
    delete h_MultiPion;
    delete h_Other;
    delete mc_hists;
    delete legend;
    delete c;
    delete plotter;
    delete f_mc;
    delete f_data;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_Enu10_BeforeFSI()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "Enu";
    std::string data_var = var_name + "_xsec";
    std::string mc_var = var_name + "_xsec_AfterFSI";
    std::string mc_var_BeforeFSI = var_name + "_xsec_BeforeFSI";

    MnvH1D* data = GetMnvH1D(f_data, data_var);
    MnvH1D* mc = GetMnvH1D(f_mc, mc_var);
    MnvH1D* mc_BeforeFSI = GetMnvH1D(f_mc, mc_var_BeforeFSI);

    data->GetXaxis()->SetRangeUser(0,10);
    mc->GetXaxis()->SetRangeUser(0,10);
    mc_BeforeFSI->GetXaxis()->SetRangeUser(0,10);

    data->GetXaxis()->SetNdivisions(5,4,0);
    mc->GetXaxis()->SetNdivisions(5,4,0);
    
    data->GetYaxis()->SetNdivisions(5,4,0);
    mc->GetYaxis()->SetNdivisions(5,4,0);
 
    data->GetYaxis()->CenterTitle();
    mc->GetYaxis()->CenterTitle();

    // Remove Error Bands from MC
    mc->ClearAllErrorBands();
    //mc_BeforeFSI->ClearAllErrorBands();
   
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    TCanvas* c = new TCanvas("c");
   
    plotter->DrawDataMCWithErrorBand(data, mc, 1.0, "N", false, NULL, NULL, false, true, true);

    TH1D* h_mc_BeforeFSI = GetBinNormalizedTH1D(mc_BeforeFSI);
    h_mc_BeforeFSI->SetLineColor(kRed);
    h_mc_BeforeFSI->SetLineWidth(2);
    h_mc_BeforeFSI->SetLineStyle(2);
    h_mc_BeforeFSI->SetFillStyle(0);

    h_mc_BeforeFSI->Smooth(2);
    //h_mc_BeforeFSI->Draw("SAME HISTC");

    // ------------------------------------------------------------------------
    // Plot Labels 
    // ------------------------------------------------------------------------
    // Add Legend
    data->SetMarkerStyle(plotter->data_marker);
    data->SetMarkerSize(plotter->data_marker_size);
    data->SetMarkerColor(plotter->data_color);
    data->SetLineColor(plotter->data_color);
    data->SetLineWidth(plotter->data_line_width);

    mc->SetLineColor(plotter->mc_color);
    mc->SetLineWidth(plotter->mc_line_width);

    double leg_x_min = 0.50;
    double leg_x_max = 0.85;
    double leg_y_min = 0.73;
    double leg_y_max = 0.88;

    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    ApplyStyle_Legend(legend);
    legend->SetTextSize(0.05);
    legend->AddEntry(data, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(mc, "GENIE", "l");
    //legend->AddEntry(h_mc_BeforeFSI, "GENIE w/o FSI","l");
    legend->Draw();

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");

    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name;
    out_name = plotDir + var_name + "10_xsec_BeforeFSI" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete legend;
    delete c;
    delete plotter;
    delete data;
    delete mc;
    delete h_mc_BeforeFSI;
    delete f_mc;
    delete f_data;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_Enu_BeforeFSI()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "Enu";
    std::string data_var = var_name + "_xsec";
    std::string mc_var = var_name + "_xsec_AfterFSI";
    std::string mc_var_BeforeFSI = var_name + "_xsec_BeforeFSI";
    std::string mc_var_NuWro = "NuWro_" + var_name;

    MnvH1D* data = GetMnvH1D(f_data, data_var);
    MnvH1D* mc = GetMnvH1D(f_mc, mc_var);
    MnvH1D* mc_BeforeFSI = GetMnvH1D(f_mc, mc_var_BeforeFSI);
    MnvH1D* mc_NuWro = GetMnvH1D(f_mc, mc_var_NuWro);

    // Average Last 10 Bins
    double avg_BeforeFSI = 0;
    double avg_AfterFSI = 0;
    for (int i = 21; i <= 30; ++i){
        avg_BeforeFSI += mc_BeforeFSI->GetBinContent(i);
        avg_AfterFSI += mc->GetBinContent(i);
    }
    avg_BeforeFSI = avg_BeforeFSI/10.0;
    avg_AfterFSI = avg_AfterFSI/10.0;
    for (int i = 21; i <= 30; ++i){
        mc_BeforeFSI->SetBinContent(i, avg_BeforeFSI);
        mc->SetBinContent(i, avg_AfterFSI);
    }

    data->GetXaxis()->SetNdivisions(5,4,0);
    mc->GetXaxis()->SetNdivisions(5,4,0);
    
    data->GetYaxis()->SetNdivisions(5,4,0);
    mc->GetYaxis()->SetNdivisions(5,4,0);
 
    data->GetYaxis()->CenterTitle();
    mc->GetYaxis()->CenterTitle();

    // Remove Error Bands from MC
    mc->ClearAllErrorBands();
    mc_BeforeFSI->ClearAllErrorBands();
   
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    TCanvas* c = new TCanvas("c");
   
    plotter->DrawDataMCWithErrorBand(data, mc, 1.0, "N", false, NULL, NULL, false, true, true);

    TH1D* h_mc_BeforeFSI = GetBinNormalizedTH1D(mc_BeforeFSI);
    h_mc_BeforeFSI->SetLineColor(kRed);
    h_mc_BeforeFSI->SetLineWidth(2);
    h_mc_BeforeFSI->SetLineStyle(2);
    h_mc_BeforeFSI->SetFillStyle(0);

    h_mc_BeforeFSI->Smooth(2);
    //h_mc_BeforeFSI->Draw("SAME HISTC");

    mc_NuWro->SetLineColor(kOrange-6);
    mc_NuWro->SetLineWidth(plotter->mc_line_width);
    mc_NuWro->SetLineStyle(7);
    mc_NuWro->SetFillStyle(0);

    mc_NuWro->Smooth(2);
    mc_NuWro->Draw("SAME HISTC");

    printBins(mc_NuWro, "NuWro_Enu", true);
    // ------------------------------------------------------------------------
    // Plot Labels 
    // ------------------------------------------------------------------------
    // Add Legend
    data->SetMarkerStyle(plotter->data_marker);
    data->SetMarkerSize(plotter->data_marker_size);
    data->SetMarkerColor(plotter->data_color);
    data->SetLineColor(plotter->data_color);
    data->SetLineWidth(plotter->data_line_width);

    mc->SetLineColor(plotter->mc_color);
    mc->SetLineWidth(plotter->mc_line_width);

    double leg_x_min = 0.50;
    double leg_x_max = 0.85;
    double leg_y_min = 0.70;
    double leg_y_max = 0.90;
 
    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    ApplyStyle_Legend(legend);
    legend->SetTextSize(0.05);
    legend->AddEntry(data, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(mc, "Simulation", "l");
    //legend->AddEntry(h_mc_BeforeFSI, "GENIE w/o FSI","l");
    legend->AddEntry(mc_NuWro, "NuWro","l");
    legend->Draw();

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");

    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name;
    out_name = plotDir + var_name + "_xsec_BeforeFSI" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete legend;
    delete c;
    delete plotter;
    delete data;
    delete mc;
    delete h_mc_BeforeFSI;
    delete mc_NuWro;
    delete f_mc;
    delete f_data;
}


void CCProtonPi0_Plotter::DrawPaper_xsec_W_BeforeFSI()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "W";
    std::string data_var = var_name + "_xsec";
    std::string mc_var = var_name + "_xsec_AfterFSI";
    std::string mc_var_BeforeFSI = var_name + "_xsec_BeforeFSI";
    std::string mc_var_NuWro = "NuWro_" + var_name;

    MnvH1D* data = GetMnvH1D(f_data, data_var);
    MnvH1D* mc = GetMnvH1D(f_mc, mc_var);
    MnvH1D* mc_BeforeFSI = GetMnvH1D(f_mc, mc_var_BeforeFSI);
    MnvH1D* mc_NuWro = GetMnvH1D(f_mc, mc_var_NuWro);
 
    data->GetXaxis()->SetNdivisions(5,4,0);
    mc->GetXaxis()->SetNdivisions(5,4,0);
    
    data->GetYaxis()->SetNdivisions(5,4,0);
    mc->GetYaxis()->SetNdivisions(5,4,0);
 
    data->GetYaxis()->CenterTitle();
    mc->GetYaxis()->CenterTitle();

    // Remove Error Bands from MC
    mc->ClearAllErrorBands();
    mc_BeforeFSI->ClearAllErrorBands();
   
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    TCanvas* c = new TCanvas("c");
   
    plotter->DrawDataMCWithErrorBand(data, mc, 1.0, "N", false, NULL, NULL, false, true, true);

    TH1D* h_mc_BeforeFSI = GetBinNormalizedTH1D(mc_BeforeFSI);
    h_mc_BeforeFSI->SetLineColor(kRed);
    h_mc_BeforeFSI->SetLineWidth(2);
    h_mc_BeforeFSI->SetLineStyle(2);
    h_mc_BeforeFSI->SetFillStyle(0);

    h_mc_BeforeFSI->Smooth(2);
    h_mc_BeforeFSI->Draw("SAME HISTC");

    mc_NuWro->SetLineColor(kOrange-6);
    mc_NuWro->SetLineWidth(plotter->mc_line_width);
    mc_NuWro->SetLineStyle(7);
    mc_NuWro->SetFillStyle(0);

    mc_NuWro->Smooth(2);
    mc_NuWro->Draw("SAME HISTC");

    // ------------------------------------------------------------------------
    // Plot Labels 
    // ------------------------------------------------------------------------
    // Add Legend
    data->SetMarkerStyle(plotter->data_marker);
    data->SetMarkerSize(plotter->data_marker_size);
    data->SetMarkerColor(plotter->data_color);
    data->SetLineColor(plotter->data_color);
    data->SetLineWidth(plotter->data_line_width);

    mc->SetLineColor(plotter->mc_color);
    mc->SetLineWidth(plotter->mc_line_width);
 
    double leg_x_min = 0.2;
    double leg_x_max = 0.57;
    double leg_y_min = 0.57;
    double leg_y_max = 0.82;
    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);

    ApplyStyle_Legend(legend);
    legend->AddEntry(data, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(mc, "GENIE w/ FSI", "l");
    legend->AddEntry(h_mc_BeforeFSI, "GENIE w/o FSI","l");
    legend->AddEntry(mc_NuWro, "NuWro","l");
    legend->Draw();


    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");

    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name;
    out_name = plotDir + var_name + "_xsec_BeforeFSI" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete legend;
    delete c;
    delete plotter;
    delete data;
    delete mc;
    delete h_mc_BeforeFSI;
    delete mc_NuWro;
    delete f_mc;
    delete f_data;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_QSq_BeforeFSI()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "QSq";
    std::string data_var = var_name + "_xsec";
    std::string mc_var = var_name + "_xsec_AfterFSI";
    std::string mc_var_BeforeFSI = var_name + "_xsec_BeforeFSI";
    std::string mc_var_NuWro = "NuWro_" + var_name;

    MnvH1D* data = GetMnvH1D(f_data, data_var);
    MnvH1D* mc = GetMnvH1D(f_mc, mc_var);
    MnvH1D* mc_BeforeFSI = GetMnvH1D(f_mc, mc_var_BeforeFSI);
    MnvH1D* mc_NuWro = GetMnvH1D(f_mc, mc_var_NuWro);
 
    data->GetXaxis()->SetNdivisions(5,4,0);
    mc->GetXaxis()->SetNdivisions(5,4,0);
    
    data->GetYaxis()->SetNdivisions(5,4,0);
    mc->GetYaxis()->SetNdivisions(5,4,0);
 
    data->GetYaxis()->CenterTitle();
    mc->GetYaxis()->CenterTitle();

    // Remove Error Bands from MC
    mc->ClearAllErrorBands();
    mc_BeforeFSI->ClearAllErrorBands();
   
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    TCanvas* c = new TCanvas("c");
   
    plotter->DrawDataMCWithErrorBand(data, mc, 1.0, "N", false, NULL, NULL, false, true, true);

    TH1D* h_mc_BeforeFSI = GetBinNormalizedTH1D(mc_BeforeFSI);
    h_mc_BeforeFSI->SetLineColor(kRed);
    h_mc_BeforeFSI->SetLineWidth(2);
    h_mc_BeforeFSI->SetLineStyle(2);
    h_mc_BeforeFSI->SetFillStyle(0);

    h_mc_BeforeFSI->Smooth(2);
    //h_mc_BeforeFSI->Draw("SAME HISTC");

    mc_NuWro->SetLineColor(kOrange-6);
    mc_NuWro->SetLineWidth(plotter->mc_line_width);
    mc_NuWro->SetLineStyle(7);
    mc_NuWro->SetFillStyle(0);

    mc_NuWro->Smooth(2);
    mc_NuWro->Draw("SAME HISTC");

    // ------------------------------------------------------------------------
    // Plot Labels 
    // ------------------------------------------------------------------------
    // Add Legend
    data->SetMarkerStyle(plotter->data_marker);
    data->SetMarkerSize(plotter->data_marker_size);
    data->SetMarkerColor(plotter->data_color);
    data->SetLineColor(plotter->data_color);
    data->SetLineWidth(plotter->data_line_width);

    mc->SetLineColor(plotter->mc_color);
    mc->SetLineWidth(plotter->mc_line_width);

    double leg_x_min = 0.50;
    double leg_x_max = 0.85;
    double leg_y_min = 0.65;
    double leg_y_max = 0.85;

    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    ApplyStyle_Legend(legend);
    legend->SetTextSize(0.05);
    legend->AddEntry(data, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(mc, "GENIE", "l");
    //legend->AddEntry(h_mc_BeforeFSI, "GENIE w/o FSI","l");
    legend->AddEntry(mc_NuWro, "NuWro","l");
    legend->Draw();


    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");

    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name;
    out_name = plotDir + var_name + "_xsec_BeforeFSI" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete legend;
    delete c;
    delete plotter;
    delete data;
    delete mc;
    delete h_mc_BeforeFSI;
    delete mc_NuWro;
    delete f_mc;
    delete f_data;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_Delta_pi_phi_BeforeFSI()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_DeltaRich_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_DeltaRich_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "Delta_pi_phi";
    std::string data_var = var_name + "_xsec";
    std::string mc_var = var_name + "_xsec_AfterFSI";
    std::string mc_var_BeforeFSI = var_name + "_xsec_BeforeFSI";
    std::string mc_var_NuWro = "NuWro_" + var_name;

    MnvH1D* data = GetMnvH1D(f_data, data_var);
    MnvH1D* mc = GetMnvH1D(f_mc, mc_var);
    MnvH1D* mc_BeforeFSI = GetMnvH1D(f_mc, mc_var_BeforeFSI);
    MnvH1D* mc_NuWro = GetMnvH1D(f_mc, mc_var_NuWro);
 
    data->GetXaxis()->SetNdivisions(5,4,0);
    mc->GetXaxis()->SetNdivisions(5,4,0);
    
    data->GetYaxis()->SetNdivisions(5,4,0);
    mc->GetYaxis()->SetNdivisions(5,4,0);
 
    data->GetYaxis()->CenterTitle();
    mc->GetYaxis()->CenterTitle();

    // Remove Error Bands from MC
    mc->ClearAllErrorBands();
    mc_BeforeFSI->ClearAllErrorBands();
   
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.001;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    TCanvas* c = new TCanvas("c");
   
    plotter->DrawDataMCWithErrorBand(data, mc, 1.0, "N", false, NULL, NULL, false, true, true);

    TH1D* h_mc_BeforeFSI = GetBinNormalizedTH1D(mc_BeforeFSI);
    h_mc_BeforeFSI->SetLineColor(kRed);
    h_mc_BeforeFSI->SetLineWidth(2);
    h_mc_BeforeFSI->SetLineStyle(2);
    h_mc_BeforeFSI->SetFillStyle(0);

    h_mc_BeforeFSI->Smooth(2);
    h_mc_BeforeFSI->Draw("SAME HISTC");

    mc_NuWro->SetLineColor(kOrange-6);
    mc_NuWro->SetLineWidth(plotter->mc_line_width);
    mc_NuWro->SetLineStyle(7);
    mc_NuWro->SetFillStyle(0);

    mc_NuWro->Smooth(2);
    mc_NuWro->Draw("SAME HISTC");

    // ------------------------------------------------------------------------
    // Plot Labels 
    // ------------------------------------------------------------------------
    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(12);
    text.DrawLatex(0.20,0.87,"POT Normalized");
    text.DrawLatex(0.20,0.83,"T_{p} > 0.1 GeV");
    text.DrawLatex(0.20,0.79,"W_{exp} < 1.4 GeV");

    double info_text_x = 0.85;
    double info_text_y = 0.85;
    TLatex info_text;
    info_text.SetNDC();
    info_text.SetTextColor(kBlack);
    info_text.SetTextFont(62);
    info_text.SetTextAlign(12);
    info_text.SetTextSize(0.04);
    info_text.DrawLatex(info_text_x, info_text_y, "(b)");

    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name;
    out_name = plotDir + var_name + "_xsec_BeforeFSI" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete c;
    delete plotter;
    delete data;
    delete mc;
    delete h_mc_BeforeFSI;
    delete f_mc;
    delete f_data;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_Delta_pi_theta_BeforeFSI()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_DeltaRich_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_DeltaRich_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "Delta_pi_theta";
    std::string data_var = var_name + "_xsec";
    std::string mc_var = var_name + "_xsec_AfterFSI";
    std::string mc_var_BeforeFSI = var_name + "_xsec_BeforeFSI";
    std::string mc_var_NuWro = "NuWro_" + var_name;

    MnvH1D* data = GetMnvH1D(f_data, data_var);
    MnvH1D* mc = GetMnvH1D(f_mc, mc_var);
    MnvH1D* mc_BeforeFSI = GetMnvH1D(f_mc, mc_var_BeforeFSI);
    MnvH1D* mc_NuWro = GetMnvH1D(f_mc, mc_var_NuWro);
 
    data->GetXaxis()->SetNdivisions(5,4,0);
    mc->GetXaxis()->SetNdivisions(5,4,0);
    
    data->GetYaxis()->SetNdivisions(5,4,0);
    mc->GetYaxis()->SetNdivisions(5,4,0);
 
    data->GetYaxis()->CenterTitle();
    mc->GetYaxis()->CenterTitle();

    // Remove Error Bands from MC
    mc->ClearAllErrorBands();
    mc_BeforeFSI->ClearAllErrorBands();
   
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    TCanvas* c = new TCanvas("c");
   
    plotter->DrawDataMCWithErrorBand(data, mc, 1.0, "N", false, NULL, NULL, false, true, true);

    TH1D* h_mc_BeforeFSI = GetBinNormalizedTH1D(mc_BeforeFSI);
    h_mc_BeforeFSI->SetLineColor(kRed);
    h_mc_BeforeFSI->SetLineWidth(2);
    h_mc_BeforeFSI->SetLineStyle(2);
    h_mc_BeforeFSI->SetFillStyle(0);

    h_mc_BeforeFSI->Smooth(2);
    h_mc_BeforeFSI->Draw("SAME HISTC");

    mc_NuWro->SetLineColor(kOrange-6);
    mc_NuWro->SetLineWidth(plotter->mc_line_width);
    mc_NuWro->SetLineStyle(7);
    mc_NuWro->SetFillStyle(0);

    mc_NuWro->Smooth(2);
    mc_NuWro->Draw("SAME HISTC");

    // ------------------------------------------------------------------------
    // Plot Labels 
    // ------------------------------------------------------------------------
    // Add Legend
    data->SetMarkerStyle(plotter->data_marker);
    data->SetMarkerSize(plotter->data_marker_size);
    data->SetMarkerColor(plotter->data_color);
    data->SetLineColor(plotter->data_color);
    data->SetLineWidth(plotter->data_line_width);

    mc->SetLineColor(plotter->mc_color);
    mc->SetLineWidth(plotter->mc_line_width);

    double leg_x_min = 0.45;
    double leg_x_max = 0.80;
    double leg_y_min = 0.59;
    double leg_y_max = 0.87;
    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);

    ApplyStyle_Legend(legend);    
    legend->AddEntry(data, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(mc, "GENIE w/ FSI", "l");
    legend->AddEntry(h_mc_BeforeFSI, "GENIE w/o FSI","l");
    legend->AddEntry(mc_NuWro, "NuWro","l");
    legend->Draw();

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(12);
    text.DrawLatex(0.20,0.87,"POT Normalized");
    text.DrawLatex(0.20,0.83,"T_{p} > 0.1 GeV");
    text.DrawLatex(0.20,0.79,"W_{exp} < 1.4 GeV");

    double info_text_x = 0.85;
    double info_text_y = 0.85;
    TLatex info_text;
    info_text.SetNDC();
    info_text.SetTextColor(kBlack);
    info_text.SetTextFont(62);
    info_text.SetTextAlign(12);
    info_text.SetTextSize(0.04);
    info_text.DrawLatex(info_text_x, info_text_y, "(a)");

    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name;
    out_name = plotDir + var_name + "_xsec_BeforeFSI" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete legend;
    delete c;
    delete plotter;
    delete data;
    delete mc;
    delete h_mc_BeforeFSI;
    delete f_mc;
    delete f_data;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_pi0_theta_BeforeFSI()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "pi0_theta";
    std::string data_var = var_name + "_xsec";
    std::string mc_var = var_name + "_xsec_AfterFSI";
    std::string mc_var_BeforeFSI = var_name + "_xsec_BeforeFSI";
    std::string mc_var_NuWro = "NuWro_" + var_name;

    MnvH1D* data = GetMnvH1D(f_data, data_var);
    MnvH1D* mc = GetMnvH1D(f_mc, mc_var);
    MnvH1D* mc_BeforeFSI = GetMnvH1D(f_mc, mc_var_BeforeFSI);
    MnvH1D* mc_NuWro = GetMnvH1D(f_mc, mc_var_NuWro);
 
    data->GetXaxis()->SetNdivisions(5,4,0);
    mc->GetXaxis()->SetNdivisions(5,4,0);
    
    data->GetYaxis()->SetNdivisions(5,4,0);
    mc->GetYaxis()->SetNdivisions(5,4,0);
 
    data->GetYaxis()->CenterTitle();
    mc->GetYaxis()->CenterTitle();

    // Remove Error Bands from MC
    mc->ClearAllErrorBands();
    mc_BeforeFSI->ClearAllErrorBands();
   
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    TCanvas* c = new TCanvas("c");
   
    plotter->DrawDataMCWithErrorBand(data, mc, 1.0, "N", false, NULL, NULL, false, true, true);

    TH1D* h_mc_BeforeFSI = GetBinNormalizedTH1D(mc_BeforeFSI);
    h_mc_BeforeFSI->SetLineColor(kRed);
    h_mc_BeforeFSI->SetLineWidth(2);
    h_mc_BeforeFSI->SetLineStyle(2);
    h_mc_BeforeFSI->SetFillStyle(0);

    h_mc_BeforeFSI->Smooth(2);
    h_mc_BeforeFSI->Draw("SAME HISTC");

    mc_NuWro->SetLineColor(kOrange-6);
    mc_NuWro->SetLineWidth(plotter->mc_line_width);
    mc_NuWro->SetLineStyle(7);
    mc_NuWro->SetFillStyle(0);

    mc_NuWro->Smooth(2);
    mc_NuWro->Draw("SAME HISTC");

    // ------------------------------------------------------------------------
    // Plot Labels 
    // ------------------------------------------------------------------------
    // Add Legend
    data->SetMarkerStyle(plotter->data_marker);
    data->SetMarkerSize(plotter->data_marker_size);
    data->SetMarkerColor(plotter->data_color);
    data->SetLineColor(plotter->data_color);
    data->SetLineWidth(plotter->data_line_width);

    mc->SetLineColor(plotter->mc_color);
    mc->SetLineWidth(plotter->mc_line_width);

    TLegend *legend = new TLegend(0.53,0.60,0.80,0.85);  
    ApplyStyle_Legend(legend);    
    legend->SetTextSize(0.05);
    legend->AddEntry(data, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(mc, "GENIE w/ FSI", "l");
    legend->AddEntry(h_mc_BeforeFSI, "GENIE w/o FSI","l");
    legend->AddEntry(mc_NuWro, "NuWro","l");
    legend->Draw();

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");

    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name;
    out_name = plotDir + var_name + "_xsec_BeforeFSI" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete legend;
    delete c;
    delete plotter;
    delete data;
    delete mc;
    delete h_mc_BeforeFSI;
    delete f_mc;
    delete f_data;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_deltaInvMass2_BeforeFSI()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_TwoTrack_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_TwoTrack_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "deltaInvMass";
    std::string data_var = var_name + "_xsec";
    std::string mc_var = var_name + "_xsec_AfterFSI";
    std::string mc_var_BeforeFSI = var_name + "_xsec_BeforeFSI";
    //std::string mc_var_NuWro = "NuWro_" + var_name;

    MnvH1D* data = GetMnvH1D(f_data, data_var);
    MnvH1D* mc = GetMnvH1D(f_mc, mc_var);
    MnvH1D* mc_BeforeFSI = GetMnvH1D(f_mc, mc_var_BeforeFSI);
    //MnvH1D* mc_NuWro = GetMnvH1D(f_mc, mc_var_NuWro);
 
    data->GetXaxis()->SetNdivisions(5,4,0);
    mc->GetXaxis()->SetNdivisions(5,4,0);
    
    data->GetYaxis()->SetNdivisions(5,4,0);
    mc->GetYaxis()->SetNdivisions(5,4,0);
 
    data->GetYaxis()->CenterTitle();
    mc->GetYaxis()->CenterTitle();

    // Remove Error Bands from MC
    mc->ClearAllErrorBands();
    mc_BeforeFSI->ClearAllErrorBands();
   
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    TCanvas* c = new TCanvas("c");
   
    plotter->DrawDataMCWithErrorBand(data, mc, 1.0, "N", false, NULL, NULL, false, true, true);

    TH1D* h_mc_BeforeFSI = GetBinNormalizedTH1D(mc_BeforeFSI);
    h_mc_BeforeFSI->SetLineColor(kRed);
    h_mc_BeforeFSI->SetLineWidth(2);
    h_mc_BeforeFSI->SetLineStyle(2);
    h_mc_BeforeFSI->SetFillStyle(0);

    h_mc_BeforeFSI->Smooth(2);
    h_mc_BeforeFSI->Draw("SAME HISTC");

//    mc_NuWro->SetLineColor(kOrange-6);
//    mc_NuWro->SetLineWidth(plotter->mc_line_width);
//    mc_NuWro->SetLineStyle(7);
//    mc_NuWro->SetFillStyle(0);
//
//    mc_NuWro->Smooth(2);
//    mc_NuWro->Draw("SAME HISTC");

    // ------------------------------------------------------------------------
    // Plot Labels 
    // ------------------------------------------------------------------------
    // Add Legend
    data->SetMarkerStyle(plotter->data_marker);
    data->SetMarkerSize(plotter->data_marker_size);
    data->SetMarkerColor(plotter->data_color);
    data->SetLineColor(plotter->data_color);
    data->SetLineWidth(plotter->data_line_width);

    mc->SetLineColor(plotter->mc_color);
    mc->SetLineWidth(plotter->mc_line_width);

    double leg_x_min = 0.55;
    double leg_x_max = 0.85;
    double leg_y_min = 0.60;
    double leg_y_max = 0.80;
    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    
    ApplyStyle_Legend(legend);    
    legend->AddEntry(data, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(mc, "GENIE w/ FSI", "l");
    legend->AddEntry(h_mc_BeforeFSI, "GENIE w/o FSI","l");
    //legend->AddEntry(mc_NuWro, "NuWro","l");
    legend->Draw();

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(12);
    text.DrawLatex(0.20,0.87,"POT Normalized");
    text.DrawLatex(0.20,0.83,"T_{p} > 0.1 GeV");

    //double info_text_x = 0.85;
    //double info_text_y = 0.85;
    //TLatex info_text;
    //info_text.SetNDC();
    //info_text.SetTextColor(kBlack);
    //info_text.SetTextFont(62);
    //info_text.SetTextAlign(12);
    //info_text.SetTextSize(0.04);
    //info_text.DrawLatex(info_text_x, info_text_y, "(a)");

    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name;
    out_name = plotDir + var_name + "2_xsec_BeforeFSI" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete legend;
    delete c;
    delete plotter;
    delete data;
    delete mc;
    delete h_mc_BeforeFSI;
    //delete mc_NuWro;
    delete f_mc;
    delete f_data;
}


void CCProtonPi0_Plotter::DrawPaper_xsec_deltaInvMass_BeforeFSI()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_DeltaRich_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_DeltaRich_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "deltaInvMass";
    std::string data_var = var_name + "_xsec";
    std::string mc_var = var_name + "_xsec_AfterFSI";
    std::string mc_var_BeforeFSI = var_name + "_xsec_BeforeFSI";
    std::string mc_var_NuWro = "NuWro_" + var_name;

    MnvH1D* data = GetMnvH1D(f_data, data_var);
    MnvH1D* mc = GetMnvH1D(f_mc, mc_var);
    MnvH1D* mc_BeforeFSI = GetMnvH1D(f_mc, mc_var_BeforeFSI);
    MnvH1D* mc_NuWro = GetMnvH1D(f_mc, mc_var_NuWro);
 
    data->GetXaxis()->SetNdivisions(5,4,0);
    mc->GetXaxis()->SetNdivisions(5,4,0);
    
    data->GetYaxis()->SetNdivisions(5,4,0);
    mc->GetYaxis()->SetNdivisions(5,4,0);
 
    data->GetYaxis()->CenterTitle();
    mc->GetYaxis()->CenterTitle();

    // Remove Error Bands from MC
    mc->ClearAllErrorBands();
    mc_BeforeFSI->ClearAllErrorBands();
   
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    TCanvas* c = new TCanvas("c");
   
    plotter->DrawDataMCWithErrorBand(data, mc, 1.0, "N", false, NULL, NULL, false, true, true);

    TH1D* h_mc_BeforeFSI = GetBinNormalizedTH1D(mc_BeforeFSI);
    h_mc_BeforeFSI->SetLineColor(kRed);
    h_mc_BeforeFSI->SetLineWidth(2);
    h_mc_BeforeFSI->SetLineStyle(2);
    h_mc_BeforeFSI->SetFillStyle(0);

    h_mc_BeforeFSI->Smooth(2);
    h_mc_BeforeFSI->Draw("SAME HISTC");

    mc_NuWro->SetLineColor(kOrange-6);
    mc_NuWro->SetLineWidth(plotter->mc_line_width);
    mc_NuWro->SetLineStyle(7);
    mc_NuWro->SetFillStyle(0);

    mc_NuWro->Smooth(2);
    mc_NuWro->Draw("SAME HISTC");

    // ------------------------------------------------------------------------
    // Plot Labels 
    // ------------------------------------------------------------------------
    // Add Legend
    data->SetMarkerStyle(plotter->data_marker);
    data->SetMarkerSize(plotter->data_marker_size);
    data->SetMarkerColor(plotter->data_color);
    data->SetLineColor(plotter->data_color);
    data->SetLineWidth(plotter->data_line_width);

    mc->SetLineColor(plotter->mc_color);
    mc->SetLineWidth(plotter->mc_line_width);

    double leg_x_min = 0.55;
    double leg_x_max = 0.85;
    double leg_y_min = 0.60;
    double leg_y_max = 0.80;
    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    
    ApplyStyle_Legend(legend);    
    legend->AddEntry(data, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(mc, "GENIE w/ FSI", "l");
    legend->AddEntry(h_mc_BeforeFSI, "GENIE w/o FSI","l");
    legend->AddEntry(mc_NuWro, "NuWro","l");
    legend->Draw();

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(12);
    text.DrawLatex(0.20,0.87,"POT Normalized");
    text.DrawLatex(0.20,0.83,"T_{p} > 0.1 GeV");
    text.DrawLatex(0.20,0.79,"W_{exp} < 1.4 GeV");

    double info_text_x = 0.85;
    double info_text_y = 0.85;
    TLatex info_text;
    info_text.SetNDC();
    info_text.SetTextColor(kBlack);
    info_text.SetTextFont(62);
    info_text.SetTextAlign(12);
    info_text.SetTextSize(0.04);
    info_text.DrawLatex(info_text_x, info_text_y, "(a)");

    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name;
    out_name = plotDir + var_name + "_xsec_BeforeFSI" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete legend;
    delete c;
    delete plotter;
    delete data;
    delete mc;
    delete h_mc_BeforeFSI;
    delete mc_NuWro;
    delete f_mc;
    delete f_data;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_pi0_KE_BeforeFSI()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "pi0_KE";
    std::string data_var = var_name + "_xsec";
    std::string mc_var = var_name + "_xsec_AfterFSI";
    std::string mc_var_BeforeFSI = var_name + "_xsec_BeforeFSI";
    std::string mc_var_NuWro = "NuWro_" + var_name;

    MnvH1D* data = GetMnvH1D(f_data, data_var);
    MnvH1D* mc = GetMnvH1D(f_mc, mc_var);
    MnvH1D* mc_BeforeFSI = GetMnvH1D(f_mc, mc_var_BeforeFSI);
    MnvH1D* mc_NuWro = GetMnvH1D(f_mc, mc_var_NuWro);
 
    data->GetXaxis()->SetNdivisions(5,4,0);
    mc->GetXaxis()->SetNdivisions(5,4,0);
    
    data->GetYaxis()->SetNdivisions(5,4,0);
    mc->GetYaxis()->SetNdivisions(5,4,0);
 
    data->GetYaxis()->CenterTitle();
    mc->GetYaxis()->CenterTitle();

    // Remove Error Bands from MC
    mc->ClearAllErrorBands();
    mc_BeforeFSI->ClearAllErrorBands();
   
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    TCanvas* c = new TCanvas("c");
   
    plotter->DrawDataMCWithErrorBand(data, mc, 1.0, "N", false, NULL, NULL, false, true, true);

    TH1D* h_mc_BeforeFSI = GetBinNormalizedTH1D(mc_BeforeFSI);
    h_mc_BeforeFSI->SetLineColor(kRed);
    h_mc_BeforeFSI->SetLineWidth(2);
    h_mc_BeforeFSI->SetLineStyle(2);
    h_mc_BeforeFSI->SetFillStyle(0);

    h_mc_BeforeFSI->Smooth(2);
    h_mc_BeforeFSI->Draw("SAME HISTC");

    mc_NuWro->SetLineColor(kOrange-6);
    mc_NuWro->SetLineWidth(plotter->mc_line_width);
    mc_NuWro->SetLineStyle(7);
    mc_NuWro->SetFillStyle(0);

    mc_NuWro->Smooth(2);
    mc_NuWro->Draw("SAME HISTC");

    // ------------------------------------------------------------------------
    // Plot Labels 
    // ------------------------------------------------------------------------
    // Add Legend
    data->SetMarkerStyle(plotter->data_marker);
    data->SetMarkerSize(plotter->data_marker_size);
    data->SetMarkerColor(plotter->data_color);
    data->SetLineColor(plotter->data_color);
    data->SetLineWidth(plotter->data_line_width);

    mc->SetLineColor(plotter->mc_color);
    mc->SetLineWidth(plotter->mc_line_width);

    TLegend *legend = new TLegend(0.53,0.60,0.80,0.85);  
    ApplyStyle_Legend(legend);    
    legend->SetTextSize(0.05);
    legend->AddEntry(data, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(mc, "GENIE w/ FSI", "l");
    legend->AddEntry(h_mc_BeforeFSI, "GENIE w/o FSI","l");
    legend->AddEntry(mc_NuWro, "NuWro","l");
    legend->Draw();

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");

    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name;
    out_name = plotDir + var_name + "_xsec_BeforeFSI" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete legend;
    delete c;
    delete plotter;
    delete data;
    delete mc;
    delete h_mc_BeforeFSI;
    delete mc_NuWro;
    delete f_mc;
    delete f_data;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_muon_theta_BeforeFSI()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "muon_theta";
    std::string data_var = var_name + "_xsec";
    std::string mc_var = var_name + "_xsec_AfterFSI";
    std::string mc_var_BeforeFSI = var_name + "_xsec_BeforeFSI";
    std::string mc_var_NuWro = "NuWro_" + var_name;

    MnvH1D* data = GetMnvH1D(f_data, data_var);
    MnvH1D* mc = GetMnvH1D(f_mc, mc_var);
    MnvH1D* mc_BeforeFSI = GetMnvH1D(f_mc, mc_var_BeforeFSI);
    MnvH1D* mc_NuWro = GetMnvH1D(f_mc, mc_var_NuWro);
 
    data->GetXaxis()->SetNdivisions(5,4,0);
    mc->GetXaxis()->SetNdivisions(5,4,0);
    
    data->GetYaxis()->SetNdivisions(5,4,0);
    mc->GetYaxis()->SetNdivisions(5,4,0);
 
    data->GetYaxis()->CenterTitle();
    mc->GetYaxis()->CenterTitle();

    // Remove Error Bands from MC
    mc->ClearAllErrorBands();
    mc_BeforeFSI->ClearAllErrorBands();
   
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    TCanvas* c = new TCanvas("c");
   
    plotter->DrawDataMCWithErrorBand(data, mc, 1.0, "N", false, NULL, NULL, false, true, true);

    TH1D* h_mc_BeforeFSI = GetBinNormalizedTH1D(mc_BeforeFSI);
    h_mc_BeforeFSI->SetLineColor(kRed);
    h_mc_BeforeFSI->SetLineWidth(2);
    h_mc_BeforeFSI->SetLineStyle(2);
    h_mc_BeforeFSI->SetFillStyle(0);

    h_mc_BeforeFSI->Smooth(2);
    //h_mc_BeforeFSI->Draw("SAME HISTC");

    mc_NuWro->SetLineColor(kOrange-6);
    mc_NuWro->SetLineWidth(plotter->mc_line_width);
    mc_NuWro->SetLineStyle(7);
    mc_NuWro->SetFillStyle(0);

    mc_NuWro->Smooth(2);
    mc_NuWro->Draw("SAME HISTC");

    // ------------------------------------------------------------------------
    // Plot Labels 
    // ------------------------------------------------------------------------
    // Add Legend
    data->SetMarkerStyle(plotter->data_marker);
    data->SetMarkerSize(plotter->data_marker_size);
    data->SetMarkerColor(plotter->data_color);
    data->SetLineColor(plotter->data_color);
    data->SetLineWidth(plotter->data_line_width);

    mc->SetLineColor(plotter->mc_color);
    mc->SetLineWidth(plotter->mc_line_width);

    double leg_x_min = 0.50;
    double leg_x_max = 0.85;
    double leg_y_min = 0.68;
    double leg_y_max = 0.88;

    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    ApplyStyle_Legend(legend);
    legend->SetTextSize(0.05);
    legend->AddEntry(data, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(mc, "GENIE", "l");
    //legend->AddEntry(h_mc_BeforeFSI, "GENIE w/o FSI","l");
    legend->AddEntry(mc_NuWro, "NuWro","l");
    legend->Draw();

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");

    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name;
    out_name = plotDir + var_name + "_xsec_BeforeFSI" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete legend;
    delete c;
    delete plotter;
    delete data;
    delete mc;
    delete h_mc_BeforeFSI;
    delete mc_NuWro;
    delete f_mc;
    delete f_data;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_muon_P_BeforeFSI()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "muon_P";
    std::string data_var = var_name + "_xsec";
    std::string mc_var = var_name + "_xsec_AfterFSI";
    std::string mc_var_BeforeFSI = var_name + "_xsec_BeforeFSI";
    std::string mc_var_NuWro = "NuWro_" + var_name;

    MnvH1D* data = GetMnvH1D(f_data, data_var);
    MnvH1D* mc = GetMnvH1D(f_mc, mc_var);
    MnvH1D* mc_BeforeFSI = GetMnvH1D(f_mc, mc_var_BeforeFSI);
    MnvH1D* mc_NuWro = GetMnvH1D(f_mc, mc_var_NuWro);
 
    data->GetXaxis()->SetNdivisions(5,4,0);
    mc->GetXaxis()->SetNdivisions(5,4,0);
    
    data->GetYaxis()->SetNdivisions(5,4,0);
    mc->GetYaxis()->SetNdivisions(5,4,0);
 
    data->GetYaxis()->CenterTitle();
    mc->GetYaxis()->CenterTitle();

    // Remove Error Bands from MC
    mc->ClearAllErrorBands();
    mc_BeforeFSI->ClearAllErrorBands();
   
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
 
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);
    gStyle->SetEndErrorSize(4);
    gStyle->SetStripDecimals(false);
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size  = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 0.8;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    TCanvas* c = new TCanvas("c");
   
    plotter->DrawDataMCWithErrorBand(data, mc, 1.0, "N", false, NULL, NULL, false, true, true);

    TH1D* h_mc_BeforeFSI = GetBinNormalizedTH1D(mc_BeforeFSI);
    h_mc_BeforeFSI->SetLineColor(kRed);
    h_mc_BeforeFSI->SetLineWidth(2);
    h_mc_BeforeFSI->SetLineStyle(2);
    h_mc_BeforeFSI->SetFillStyle(0);

    h_mc_BeforeFSI->Smooth(2);
    //h_mc_BeforeFSI->Draw("SAME HISTC");

    mc_NuWro->SetLineColor(kOrange-6);
    mc_NuWro->SetLineWidth(plotter->mc_line_width);
    mc_NuWro->SetLineStyle(7);
    mc_NuWro->SetFillStyle(0);

    mc_NuWro->Smooth(2);
    mc_NuWro->Draw("SAME HISTC");

    // ------------------------------------------------------------------------
    // Plot Labels 
    // ------------------------------------------------------------------------
    // Add Legend
    data->SetMarkerStyle(plotter->data_marker);
    data->SetMarkerSize(plotter->data_marker_size);
    data->SetMarkerColor(plotter->data_color);
    data->SetLineColor(plotter->data_color);
    data->SetLineWidth(plotter->data_line_width);

    mc->SetLineColor(plotter->mc_color);
    mc->SetLineWidth(plotter->mc_line_width);

    double leg_x_min = 0.50;
    double leg_x_max = 0.85;
    double leg_y_min = 0.65;
    double leg_y_max = 0.85;

    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    ApplyStyle_Legend(legend);
    legend->SetTextSize(0.05);
    legend->AddEntry(data, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(mc, "GENIE", "l");
    //legend->AddEntry(h_mc_BeforeFSI, "GENIE w/o FSI","l");
    legend->AddEntry(mc_NuWro, "NuWro","l");
    legend->Draw();

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");

    // Plot Output
    c->Update();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name;
    out_name = plotDir + var_name + "_xsec_BeforeFSI" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete legend;
    delete c;
    delete plotter;
    delete data;
    delete mc;
    delete h_mc_BeforeFSI;
    delete mc_NuWro;
    delete f_mc;
    delete f_data;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_Enu()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var = "Enu_xsec";

    MnvH1D* data = GetMnvH1D(f_data, var);
    MnvH1D* mc = GetMnvH1D(f_mc, var);
   
    data->GetXaxis()->SetNdivisions(5,4,0);
    mc->GetXaxis()->SetNdivisions(5,4,0);
    
    data->GetYaxis()->SetNdivisions(5,4,0);
    mc->GetYaxis()->SetNdivisions(5,4,0);
 
    data->GetYaxis()->CenterTitle();
    mc->GetYaxis()->CenterTitle();

    data->SetTitle("Data (3.33e20 POT)");
    mc->SetTitle("Simulation");

    //data->GetYaxis()->ChangeLabel(0,-1.,0);
    //data->GetYaxis()->ChangeLabel(1,-1.,0);

    mc->ClearAllErrorBands();
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
    
    gStyle->SetStripDecimals(false);
    gStyle->SetEndErrorSize(4);
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); 
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)

    plotter->data_line_width = 2;

    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;
  
    TCanvas* c = new TCanvas("c");

    //    void MnvPlotter::DrawDataMCWithErrorBand    (   const MnvH1D *  dataHist,
    //            const MnvH1D *  mcHist,
    //            const Double_t  mcScale = 1.0,
    //            const std::string &     legPos = "L",
    //            const bool  useHistTitles = false,
    //            const MnvH1D *  bkgdHist = NULL,
    //            const MnvH1D *  dataBkgdHist = NULL,
    //            const bool  covAreaNormalize = false,
    //            const bool  statPlusSys = false  
    //            const bool isSmooth = false
    //            )   
    plotter->DrawDataMCWithErrorBand(data, mc, 1.0, "N", true, NULL, NULL, false, true, false);

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");

    // Add Legend
    data->SetMarkerStyle(plotter->data_marker);
    data->SetMarkerSize(plotter->data_marker_size);
    data->SetMarkerColor(plotter->data_color);
    data->SetLineColor(plotter->data_color);
    data->SetLineWidth(plotter->data_line_width);

    mc->SetLineColor(plotter->mc_color);
    mc->SetLineWidth(plotter->mc_line_width);
    mc->SetLineStyle(plotter->mc_line_style);
    mc->SetFillColor(plotter->mc_error_color);
    mc->SetFillStyle(plotter->mc_error_style);

    double leg_x_min = 0.50;
    double leg_x_max = 0.80;
    double leg_y_min = 0.70;
    double leg_y_max = 0.90;

    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    ApplyStyle_Legend(legend);

    legend->AddEntry(data, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(mc, "Simulation", "fl");
    legend->Draw();

    // Plot Output
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");
 
    delete c;
    delete plotter;
    
    // ------------------------------------------------------------------------
    // Plot Error Summary 
    // ------------------------------------------------------------------------
    //DrawErrorSummary(data, "data_" + var, plotDir);

    delete data;
    delete mc;
    delete f_data;
    delete f_mc;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_W()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var = "W_xsec";

    MnvH1D* data = GetMnvH1D(f_data, var);
    MnvH1D* mc = GetMnvH1D(f_mc, var);
   
    data->GetXaxis()->SetNdivisions(5,4,0);
    mc->GetXaxis()->SetNdivisions(5,4,0);
    
    data->GetYaxis()->SetNdivisions(5,4,0);
    mc->GetYaxis()->SetNdivisions(5,4,0);
 
    data->GetYaxis()->CenterTitle();
    mc->GetYaxis()->CenterTitle();

    data->SetTitle("Data (3.33e20 POT)");
    mc->SetTitle("Simulation");

    //data->GetYaxis()->ChangeLabel(0,-1.,0);
    //data->GetYaxis()->ChangeLabel(1,-1.,0);

    mc->ClearAllErrorBands();
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
    
    gStyle->SetStripDecimals(false);
    gStyle->SetEndErrorSize(4);
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); 
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)

    plotter->data_line_width = 2;

    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;
  
    TCanvas* c = new TCanvas("c");

    //    void MnvPlotter::DrawDataMCWithErrorBand    (   const MnvH1D *  dataHist,
    //            const MnvH1D *  mcHist,
    //            const Double_t  mcScale = 1.0,
    //            const std::string &     legPos = "L",
    //            const bool  useHistTitles = false,
    //            const MnvH1D *  bkgdHist = NULL,
    //            const MnvH1D *  dataBkgdHist = NULL,
    //            const bool  covAreaNormalize = false,
    //            const bool  statPlusSys = false  
    //            const bool isSmooth = false
    //            )   
    plotter->DrawDataMCWithErrorBand(data, mc, 1.0, "N", true, NULL, NULL, false, true, false);

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");

    // Add Legend
    data->SetMarkerStyle(plotter->data_marker);
    data->SetMarkerSize(plotter->data_marker_size);
    data->SetMarkerColor(plotter->data_color);
    data->SetLineColor(plotter->data_color);
    data->SetLineWidth(plotter->data_line_width);

    mc->SetLineColor(plotter->mc_color);
    mc->SetLineWidth(plotter->mc_line_width);
    mc->SetLineStyle(plotter->mc_line_style);
    mc->SetFillColor(plotter->mc_error_color);
    mc->SetFillStyle(plotter->mc_error_style);

    double leg_x_min = 0.50;
    double leg_x_max = 0.80;
    double leg_y_min = 0.70;
    double leg_y_max = 0.90;

    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    ApplyStyle_Legend(legend);

    legend->AddEntry(data, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(mc, "Simulation", "fl");
    legend->Draw();

    // Plot Output
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");
 
    delete c;
    delete plotter;
    
    // ------------------------------------------------------------------------
    // Plot Error Summary 
    // ------------------------------------------------------------------------
    //DrawErrorSummary(data, "data_" + var, plotDir);

    delete data;
    delete mc;
    delete f_data;
    delete f_mc;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_QSq()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var = "QSq_xsec";

    MnvH1D* data = GetMnvH1D(f_data, var);
    MnvH1D* mc = GetMnvH1D(f_mc, var);
   
    data->GetXaxis()->SetNdivisions(5,4,0);
    mc->GetXaxis()->SetNdivisions(5,4,0);
    
    data->GetYaxis()->SetNdivisions(5,4,0);
    mc->GetYaxis()->SetNdivisions(5,4,0);
 
    data->GetYaxis()->CenterTitle();
    mc->GetYaxis()->CenterTitle();

    data->SetTitle("Data (3.33e20 POT)");
    mc->SetTitle("Simulation");

    mc->ClearAllErrorBands();
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
    
    gStyle->SetStripDecimals(false);
    gStyle->SetEndErrorSize(4);
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); 
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->data_line_width = 2;

    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 0.8;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;
  
    TCanvas* c = new TCanvas("c");

    //    void MnvPlotter::DrawDataMCWithErrorBand    (   const MnvH1D *  dataHist,
    //            const MnvH1D *  mcHist,
    //            const Double_t  mcScale = 1.0,
    //            const std::string &     legPos = "L",
    //            const bool  useHistTitles = false,
    //            const MnvH1D *  bkgdHist = NULL,
    //            const MnvH1D *  dataBkgdHist = NULL,
    //            const bool  covAreaNormalize = false,
    //            const bool  statPlusSys = false  
    //            const bool isSmooth = false
    //            )   
    plotter->DrawDataMCWithErrorBand(data, mc, 1.0, "N", true, NULL, NULL, false, true, false);

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");

    // Add Legend
    data->SetMarkerStyle(plotter->data_marker);
    data->SetMarkerSize(plotter->data_marker_size);
    data->SetMarkerColor(plotter->data_color);
    data->SetLineColor(plotter->data_color);
    data->SetLineWidth(plotter->data_line_width);

    mc->SetLineColor(plotter->mc_color);
    mc->SetLineWidth(plotter->mc_line_width);
    mc->SetLineStyle(plotter->mc_line_style);
    mc->SetFillColor(plotter->mc_error_color);
    mc->SetFillStyle(plotter->mc_error_style);

    double leg_x_min = 0.50;
    double leg_x_max = 0.80;
    double leg_y_min = 0.70;
    double leg_y_max = 0.90;

    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    ApplyStyle_Legend(legend);

    legend->AddEntry(data, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(mc, "Simulation", "fl");
    legend->Draw();

    // Plot Output
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");
 
    delete c;
    delete plotter;
    
    // ------------------------------------------------------------------------
    // Plot Error Summary 
    // ------------------------------------------------------------------------
    //DrawErrorSummary(data, "data_" + var, plotDir);

    delete data;
    delete mc;
    delete f_data;
    delete f_mc;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_pi0_theta()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var = "pi0_theta_xsec";

    MnvH1D* data = GetMnvH1D(f_data, var);
    MnvH1D* mc = GetMnvH1D(f_mc, var);
   
    data->GetXaxis()->SetNdivisions(5,4,0);
    mc->GetXaxis()->SetNdivisions(5,4,0);
    
    data->GetYaxis()->SetNdivisions(5,4,0);
    mc->GetYaxis()->SetNdivisions(5,4,0);
 
    data->GetYaxis()->CenterTitle();
    mc->GetYaxis()->CenterTitle();

    data->SetTitle("Data (3.33e20 POT)");
    mc->SetTitle("Simulation");
    
    mc->ClearAllErrorBands();
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
    
    gStyle->SetStripDecimals(false);
    gStyle->SetEndErrorSize(4);
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); 
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)

    plotter->data_line_width = 2;

    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;
  
    TCanvas* c = new TCanvas("c");

    //    void MnvPlotter::DrawDataMCWithErrorBand    (   const MnvH1D *  dataHist,
    //            const MnvH1D *  mcHist,
    //            const Double_t  mcScale = 1.0,
    //            const std::string &     legPos = "L",
    //            const bool  useHistTitles = false,
    //            const MnvH1D *  bkgdHist = NULL,
    //            const MnvH1D *  dataBkgdHist = NULL,
    //            const bool  covAreaNormalize = false,
    //            const bool  statPlusSys = false  
    //            const bool isSmooth = false
    //            )   
    plotter->DrawDataMCWithErrorBand(data, mc, 1.0, "N", true, NULL, NULL, false, true, false);

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");

    // Add Legend
    data->SetMarkerStyle(plotter->data_marker);
    data->SetMarkerSize(plotter->data_marker_size);
    data->SetMarkerColor(plotter->data_color);
    data->SetLineColor(plotter->data_color);
    data->SetLineWidth(plotter->data_line_width);

    mc->SetLineColor(plotter->mc_color);
    mc->SetLineWidth(plotter->mc_line_width);
    mc->SetLineStyle(plotter->mc_line_style);
    mc->SetFillColor(plotter->mc_error_color);
    mc->SetFillStyle(plotter->mc_error_style);

    double leg_x_min = 0.50;
    double leg_x_max = 0.80;
    double leg_y_min = 0.70;
    double leg_y_max = 0.90;

    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    ApplyStyle_Legend(legend);

    legend->AddEntry(data, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(mc, "Simulation", "fl");
    legend->Draw();

 
    // Plot Output
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");
 
    delete c;
    delete plotter;
    
    // ------------------------------------------------------------------------
    // Plot Error Summary 
    // ------------------------------------------------------------------------
    //DrawErrorSummary(data, "data_" + var, plotDir);

    delete data;
    delete mc;
    delete f_data;
    delete f_mc;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_pi0_KE()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var = "pi0_KE_xsec";

    MnvH1D* data = GetMnvH1D(f_data, var);
    MnvH1D* mc = GetMnvH1D(f_mc, var);
   
    data->GetXaxis()->SetNdivisions(5,4,0);
    mc->GetXaxis()->SetNdivisions(5,4,0);
    
    data->GetYaxis()->SetNdivisions(5,4,0);
    mc->GetYaxis()->SetNdivisions(5,4,0);
 
    data->GetYaxis()->CenterTitle();
    mc->GetYaxis()->CenterTitle();

    data->SetTitle("Data (3.33e20 POT)");
    mc->SetTitle("Simulation");
    
    mc->ClearAllErrorBands();
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
    
    gStyle->SetStripDecimals(false);
    gStyle->SetEndErrorSize(4);
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); 
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)

    plotter->data_line_width = 2;

    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;
  
    TCanvas* c = new TCanvas("c");

    //    void MnvPlotter::DrawDataMCWithErrorBand    (   const MnvH1D *  dataHist,
    //            const MnvH1D *  mcHist,
    //            const Double_t  mcScale = 1.0,
    //            const std::string &     legPos = "L",
    //            const bool  useHistTitles = false,
    //            const MnvH1D *  bkgdHist = NULL,
    //            const MnvH1D *  dataBkgdHist = NULL,
    //            const bool  covAreaNormalize = false,
    //            const bool  statPlusSys = false  
    //            const bool isSmooth = false
    //            )   
    plotter->DrawDataMCWithErrorBand(data, mc, 1.0, "N", true, NULL, NULL, false, true, false);

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");

    // Add Legend
    data->SetMarkerStyle(plotter->data_marker);
    data->SetMarkerSize(plotter->data_marker_size);
    data->SetMarkerColor(plotter->data_color);
    data->SetLineColor(plotter->data_color);
    data->SetLineWidth(plotter->data_line_width);

    mc->SetLineColor(plotter->mc_color);
    mc->SetLineWidth(plotter->mc_line_width);
    mc->SetLineStyle(plotter->mc_line_style);
    mc->SetFillColor(plotter->mc_error_color);
    mc->SetFillStyle(plotter->mc_error_style);

    double leg_x_min = 0.50;
    double leg_x_max = 0.80;
    double leg_y_min = 0.70;
    double leg_y_max = 0.90;

    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    ApplyStyle_Legend(legend);

    legend->AddEntry(data, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(mc, "Simulation", "fl");
    legend->Draw();

    // Plot Output
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");
 
    delete c;
    delete plotter;
    
    // ------------------------------------------------------------------------
    // Plot Error Summary 
    // ------------------------------------------------------------------------
    //DrawErrorSummary(data, "data_" + var, plotDir);

    delete data;
    delete mc;
    delete f_data;
    delete f_mc;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_muon_theta()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var = "muon_theta_xsec";

    MnvH1D* data = GetMnvH1D(f_data, var);
    MnvH1D* mc = GetMnvH1D(f_mc, var);
   
    data->GetXaxis()->SetNdivisions(5,4,0);
    mc->GetXaxis()->SetNdivisions(5,4,0);
    
    data->GetYaxis()->SetNdivisions(5,4,0);
    mc->GetYaxis()->SetNdivisions(5,4,0);
 
    data->GetYaxis()->CenterTitle();
    mc->GetYaxis()->CenterTitle();

    data->SetTitle("Data (3.33e20 POT)");
    mc->SetTitle("Simulation");

    //data->GetYaxis()->ChangeLabel(0,-1.,0);
    //data->GetYaxis()->ChangeLabel(1,-1.,0);

    mc->ClearAllErrorBands();
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
    
    gStyle->SetStripDecimals(false);
    gStyle->SetEndErrorSize(4);
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); 
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)

    plotter->data_line_width = 2;

    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;
  
    TCanvas* c = new TCanvas("c");

    //    void MnvPlotter::DrawDataMCWithErrorBand    (   const MnvH1D *  dataHist,
    //            const MnvH1D *  mcHist,
    //            const Double_t  mcScale = 1.0,
    //            const std::string &     legPos = "L",
    //            const bool  useHistTitles = false,
    //            const MnvH1D *  bkgdHist = NULL,
    //            const MnvH1D *  dataBkgdHist = NULL,
    //            const bool  covAreaNormalize = false,
    //            const bool  statPlusSys = false  
    //            const bool isSmooth = false
    //            )   
    plotter->DrawDataMCWithErrorBand(data, mc, 1.0, "N", true, NULL, NULL, false, true, false);

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");

    // Add Legend
    data->SetMarkerStyle(plotter->data_marker);
    data->SetMarkerSize(plotter->data_marker_size);
    data->SetMarkerColor(plotter->data_color);
    data->SetLineColor(plotter->data_color);
    data->SetLineWidth(plotter->data_line_width);

    mc->SetLineColor(plotter->mc_color);
    mc->SetLineWidth(plotter->mc_line_width);
    mc->SetLineStyle(plotter->mc_line_style);
    mc->SetFillColor(plotter->mc_error_color);
    mc->SetFillStyle(plotter->mc_error_style);

    double leg_x_min = 0.50;
    double leg_x_max = 0.80;
    double leg_y_min = 0.70;
    double leg_y_max = 0.90;

    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    ApplyStyle_Legend(legend);

    legend->AddEntry(data, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(mc, "Simulation", "fl");
    legend->Draw();

    // Plot Output
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");
 
    delete c;
    delete plotter;
    
    // ------------------------------------------------------------------------
    // Plot Error Summary 
    // ------------------------------------------------------------------------
    //DrawErrorSummary(data, "data_" + var, plotDir);

    delete data;
    delete mc;
    delete f_data;
    delete f_mc;
}

void CCProtonPi0_Plotter::DrawPaper_xsec_muon_P()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var = "muon_P_xsec";

    MnvH1D* data = GetMnvH1D(f_data, var);
    MnvH1D* mc = GetMnvH1D(f_mc, var);
   
    data->GetXaxis()->SetNdivisions(5,4,0);
    mc->GetXaxis()->SetNdivisions(5,4,0);
    
    data->GetYaxis()->SetNdivisions(5,4,0);
    mc->GetYaxis()->SetNdivisions(5,4,0);
 
    data->GetYaxis()->CenterTitle();
    mc->GetYaxis()->CenterTitle();

    data->SetTitle("Data (3.33e20 POT)");
    mc->SetTitle("Simulation");
   
    mc->ClearAllErrorBands();
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
    
    gStyle->SetStripDecimals(false);
    gStyle->SetEndErrorSize(4);
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); 
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)

    plotter->data_line_width = 2;

    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 0.8;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;
  
    TCanvas* c = new TCanvas("c");

    //    void MnvPlotter::DrawDataMCWithErrorBand    (   const MnvH1D *  dataHist,
    //            const MnvH1D *  mcHist,
    //            const Double_t  mcScale = 1.0,
    //            const std::string &     legPos = "L",
    //            const bool  useHistTitles = false,
    //            const MnvH1D *  bkgdHist = NULL,
    //            const MnvH1D *  dataBkgdHist = NULL,
    //            const bool  covAreaNormalize = false,
    //            const bool  statPlusSys = false  
    //            const bool isSmooth = false
    //            )   
    plotter->DrawDataMCWithErrorBand(data, mc, 1.0, "N", true, NULL, NULL, false, true, false);

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");

    // Add Legend
    data->SetMarkerStyle(plotter->data_marker);
    data->SetMarkerSize(plotter->data_marker_size);
    data->SetMarkerColor(plotter->data_color);
    data->SetLineColor(plotter->data_color);
    data->SetLineWidth(plotter->data_line_width);

    mc->SetLineColor(plotter->mc_color);
    mc->SetLineWidth(plotter->mc_line_width);
    mc->SetLineStyle(plotter->mc_line_style);
    mc->SetFillColor(plotter->mc_error_color);
    mc->SetFillStyle(plotter->mc_error_style);

    double leg_x_min = 0.50;
    double leg_x_max = 0.80;
    double leg_y_min = 0.70;
    double leg_y_max = 0.90;

    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    ApplyStyle_Legend(legend);

    legend->AddEntry(data, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(mc, "Simulation", "fl");
    legend->Draw();

    // Plot Output
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");
 
    delete c;
    delete plotter;
    
    // ------------------------------------------------------------------------
    // Plot Error Summary 
    // ------------------------------------------------------------------------
    //DrawErrorSummary(data, "data_" + var, plotDir);

    delete data;
    delete mc;
    delete f_data;
    delete f_mc;
}

void CCProtonPi0_Plotter::DrawPaper_InvMass_DataMC()
{
    std::string root_dir_data = Folder_List::rootDir_CutHists_Signal_Original_Data;
    std::string root_dir_mc = Folder_List::rootDir_CutHists_Signal_Original_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "hCut_pi0invMass";
    std::string var;

    // ------------------------------------------------------------------------
    // Get Data Histogram
    // ------------------------------------------------------------------------
    var = Form("%s_%d",var_name.c_str(),0);
    MnvH1D* data = (MnvH1D*)f_data->Get(var.c_str()); 
    data->SetMarkerStyle(20);
    data->SetMarkerSize(1);
    data->SetMarkerColor(kBlack);
    data->SetLineWidth(2);
    data->SetLineColor(kBlack);
    data->GetXaxis()->SetNdivisions(5,5,0);
    data->GetYaxis()->SetNdivisions(5,5,0);
    data->GetYaxis()->CenterTitle();
    // ------------------------------------------------------------------------
    // Fill TObjArray - For MC Histograms
    // ------------------------------------------------------------------------
    TObjArray* mc_hists = new TObjArray;

    // Get All Background
    var = Form("%s_%d",var_name.c_str(),2);
    MnvH1D* background = (MnvH1D*)f_mc->Get(var.c_str());
    background->SetTitle("Background");
    background->SetLineWidth(2);
    background->SetLineColor(kGray+2);
    background->SetFillColor(kGray+2);
    background->SetFillStyle(3002);
    background->GetXaxis()->SetNdivisions(5,5,0);
    background->GetYaxis()->SetNdivisions(5,5,0);
    background->GetYaxis()->CenterTitle();
    mc_hists->Add(background);

    // Get Signal
    var = Form("%s_%d",var_name.c_str(),1);
    MnvH1D* signal = (MnvH1D*)f_mc->Get(var.c_str());
    signal->SetTitle("Signal");
    signal->SetLineColor(kGreen+1);
    signal->SetFillColor(kGreen+1);
    signal->SetFillStyle(3002);
    signal->GetXaxis()->SetNdivisions(5,5,0);
    signal->GetYaxis()->SetNdivisions(5,5,0);
    signal->GetYaxis()->CenterTitle();
    mc_hists->Add(signal);

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

    TCanvas* c = new TCanvas("c");
    plotter->DrawDataStackedMC(data, mc_hists, POT_ratio,"N","Data",0,0, 3002);

    // Add Normalization Label
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.35,0.85,"POT Normalized");

    // ------------------------------------------------------------------------
    // Lines 
    // ------------------------------------------------------------------------
    double max_bin = data->GetMaximumBin();
    double hist_max = data->GetBinContent(max_bin)*1.2;

    TLine line;
    TArrow arrow;
    // Cut Line at 60 MeV
    line.SetLineWidth(2);
    line.SetLineStyle(2);
    line.SetLineColor(kBlack);
    arrow.SetLineWidth(2);
    arrow.SetLineColor(kBlack);
    line.DrawLine(60.0, 0.0, 60.0, hist_max);
    arrow.DrawArrow(60.0, hist_max, 70.0, hist_max, 0.01);

    line.DrawLine(200.0, 0.0, 200.0, hist_max);
    arrow.DrawArrow(200.0, hist_max, 190.0, hist_max, 0.01);

    arrow.SetLineColor(kBlue);
    arrow.DrawArrow(134.98,hist_max/8,134.98,0, 0.01); 

    // ------------------------------------------------------------------------
    // Legend 
    // ------------------------------------------------------------------------
    double leg_x_min = 0.60;
    double leg_x_max = 0.90;
    double leg_y_min = 0.58;
    double leg_y_max = 0.78;

    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    ApplyStyle_Legend(legend);
    legend->AddEntry(data, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(signal, "Signal","f");
    legend->AddEntry(background, "Tuned Background","f");
    legend->Draw();
 
    double info_text_x = 0.85;
    double info_text_y = 0.85;
    TLatex info_text;
    info_text.SetNDC();
    info_text.SetTextColor(kBlack);
    info_text.SetTextFont(62);
    info_text.SetTextAlign(12);
    info_text.SetTextSize(0.04);
    info_text.DrawLatex(info_text_x, info_text_y, "(a)");

    // Print Plot
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string plot_out = plotDir + "pi0_invMass_DataMC.pdf";
    c->Print(plot_out.c_str(), "pdf");

    delete c;
    delete plotter;
    delete f_data;
    delete f_mc;
}

void CCProtonPi0_Plotter::DrawPaper_InvMass_SignalMC()
{
    std::string root_dir_mc = Folder_List::rootDir_CutHists_Signal_Original_MC;
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var_name = "hCut_pi0invMass";
    std::string var;

    // ------------------------------------------------------------------------
    // Fill TObjArray - For MC Histograms
    // ------------------------------------------------------------------------
    TObjArray* mc_hists = new TObjArray;

    // Get Signal
    var = Form("%s_%d",var_name.c_str(),1);
    MnvH1D* signal = (MnvH1D*)f_mc->Get(var.c_str());
    signal->SetTitle("Signal");
    signal->GetYaxis()->CenterTitle();
    signal->SetLineColor(kGreen+1);
    signal->SetFillColor(kGreen+1);
    signal->SetFillStyle(3002);
    signal->Scale(POT_ratio);
    signal->GetXaxis()->SetNdivisions(5,5,0);
    signal->GetYaxis()->SetNdivisions(5,5,0);
    mc_hists->Add(signal);

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
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.0;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;
    plotter->headroom = 1.75;

    TCanvas* c = new TCanvas("c");
    plotter->DrawStackedMC(mc_hists, 1.0,"N",0,0,3002);

    // ------------------------------------------------------------------------
    // Lines 
    // ------------------------------------------------------------------------
    double max_bin = signal->GetMaximumBin();
    double hist_max = signal->GetBinContent(max_bin)*1.2;

    TLine line;
    TArrow arrow;
    // Cut Line at 60 MeV
    line.SetLineWidth(2);
    line.SetLineStyle(2);
    line.SetLineColor(kBlack);
    arrow.SetLineWidth(2);
    arrow.SetLineColor(kBlack);
    line.DrawLine(60.0, 0.0, 60.0, hist_max);
    arrow.DrawArrow(60.0, hist_max, 70.0, hist_max, 0.01);

    line.DrawLine(200.0, 0.0, 200.0, hist_max);
    arrow.DrawArrow(200.0, hist_max, 190.0, hist_max, 0.01);

    arrow.SetLineColor(kBlue);
    arrow.DrawArrow(134.98,hist_max/8,134.98,0, 0.01); 

    // Add Normalization Label
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.35,0.85,"POT Normalized");

    text.SetTextFont(42);
    text.SetTextAlign(12);
    text.SetTextColor(kBlack);
    text.SetTextSize(0.05);
    text.DrawLatex(0.60,0.75,"Signal (MC)");
    text.SetTextSize(0.04);
    text.DrawLatex(0.60,0.70,"(after subtraction)");

    double info_text_x = 0.85;
    double info_text_y = 0.85;
    TLatex info_text;
    info_text.SetNDC();
    info_text.SetTextColor(kBlack);
    info_text.SetTextFont(62);
    info_text.SetTextAlign(12);
    info_text.SetTextSize(0.04);
    info_text.DrawLatex(info_text_x, info_text_y, "(b)");

    // Print Plot
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string plot_out = plotDir + "pi0_invMass_SignalMC.pdf";
    c->Print(plot_out.c_str(), "pdf");

    delete c;
    delete plotter;
}

void CCProtonPi0_Plotter::DrawPaper_Delta_pi_phi_2Bin()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_2Bin_DeltaRich_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_2Bin_DeltaRich_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var = "Delta_pi_phi_xsec";

    MnvH1D* data = GetMnvH1D(f_data, var);
    MnvH1D* mc = GetMnvH1D(f_mc, var);
   
    data->GetXaxis()->SetNdivisions(5,4,0);
    mc->GetXaxis()->SetNdivisions(5,4,0);
    
    data->GetYaxis()->SetNdivisions(5,4,0);
    mc->GetYaxis()->SetNdivisions(5,4,0);
 
    data->GetYaxis()->CenterTitle();
    mc->GetYaxis()->CenterTitle();

    data->SetTitle("Data (3.33e20 POT)");
    mc->SetTitle("Simulation");
   
    mc->ClearAllErrorBands();

    std::string norm_label; 
    double mc_ratio = GetMCNormalization(norm_label, false, data, mc);
    (void) norm_label;

    mc->Scale(mc_ratio);

    std::cout<<"---------------------------------------------------------------"<<std::endl;
    std::cout<<"Phi Asymmetry"<<std::endl;
    TH1D* err_data = new TH1D(data->GetTotalError());
    double n_right = data->GetBinContent(1);
    double n_right_mc = mc->GetBinContent(1);
    double err_right = err_data->GetBinContent(1);
    double n_left = data->GetBinContent(2);
    double n_left_mc = mc->GetBinContent(2);
    double err_left = err_data->GetBinContent(2);

    double ratio_data = n_left/n_right;
    double ratio_mc = n_left_mc/n_right_mc;
    double err_ratio = std::sqrt(std::pow(err_left/n_left,2) + std::pow(err_right/n_right,2)) ;

    double asymmetry = ratio_data/ratio_mc;
    double err_asymmetry = err_ratio * asymmetry;

    std::cout<<"left = "<<n_left<<" "<<n_left_mc<<" "<<err_left<<std::endl;
    std::cout<<"right = "<<n_right<<" "<<n_right_mc<<" "<<err_right<<std::endl;
    std::cout<<"Asymmetry = "<<asymmetry<<" "<<err_asymmetry<<std::endl;
    std::cout<<"---------------------------------------------------------------"<<std::endl;

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
    
    gStyle->SetStripDecimals(false);
    gStyle->SetEndErrorSize(4);
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); 
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)

    plotter->data_line_width = 2;

    plotter->axis_minimum = 0.001;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 0.8;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;
  
    //printBins(data,var + "_data",false);
    //printBins(mc,var + "_MC" ,false);

    TCanvas* canvas = new TCanvas("c");

    //    void MnvPlotter::DrawDataMCWithErrorBand    (   const MnvH1D *  dataHist,
    //            const MnvH1D *  mcHist,
    //            const Double_t  mcScale = 1.0,
    //            const std::string &     legPos = "L",
    //            const bool  useHistTitles = false,
    //            const MnvH1D *  bkgdHist = NULL,
    //            const MnvH1D *  dataBkgdHist = NULL,
    //            const bool  covAreaNormalize = false,
    //            const bool  statPlusSys = false  
    //            const bool isSmooth = false
    //            )   
    plotter->DrawDataMCWithErrorBand(data, mc, 1.0, "N", true, NULL, NULL, false, true, false);

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(12);
    text.DrawLatex(0.20,0.87,"Area Normalized");
    text.DrawLatex(0.20,0.83,"T_{p} > 0.1 GeV");
    text.DrawLatex(0.20,0.79,"W_{exp} < 1.4 GeV");

    double info_text_x = 0.85;
    double info_text_y = 0.85;
    TLatex info_text;
    info_text.SetNDC();
    info_text.SetTextColor(kBlack);
    info_text.SetTextFont(62);
    info_text.SetTextAlign(12);
    info_text.SetTextSize(0.04);
    info_text.DrawLatex(info_text_x, info_text_y, "(b)");

    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var + "_2Bin" + ".pdf"; 
    canvas->Print(out_name.c_str(),"pdf");
    
    delete f_data;
    delete f_mc;
    delete data;
    delete mc;
    delete canvas;
    delete plotter;
}

void CCProtonPi0_Plotter::DrawPaper_Delta_pi_theta_2Bin()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_2Bin_DeltaRich_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_2Bin_DeltaRich_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var = "Delta_pi_theta_xsec";

    MnvH1D* data = GetMnvH1D(f_data, var);
    MnvH1D* mc = GetMnvH1D(f_mc, var);
   
    data->GetXaxis()->SetNdivisions(5,4,0);
    mc->GetXaxis()->SetNdivisions(5,4,0);
    
    data->GetYaxis()->SetNdivisions(5,4,0);
    mc->GetYaxis()->SetNdivisions(5,4,0);
 
    data->GetYaxis()->CenterTitle();
    mc->GetYaxis()->CenterTitle();

    data->SetTitle("Data (3.33e20 POT)");
    mc->SetTitle("Simulation");
   
    mc->ClearAllErrorBands();

    std::string norm_label; 
    double mc_ratio = GetMCNormalization(norm_label, false, data, mc);
    (void) norm_label;

    mc->Scale(mc_ratio);

    std::cout<<"---------------------------------------------------------------"<<std::endl;
    std::cout<<"cos(theta) Asymmetry"<<std::endl;
    TH1D* err_data = new TH1D(data->GetTotalError());
    double n_backward = data->GetBinContent(1);
    double n_backward_mc = mc->GetBinContent(1);
    double err_backward = err_data->GetBinContent(1);
    double n_forward = data->GetBinContent(2);
    double n_forward_mc = mc->GetBinContent(2);
    double err_forward = err_data->GetBinContent(2);

    double ratio_data = n_forward/n_backward;
    double ratio_mc = n_forward_mc/n_backward_mc;
    double err_ratio = std::sqrt(std::pow(err_forward/n_forward,2) + std::pow(err_backward/n_backward,2)) ;

    double asymmetry = ratio_data/ratio_mc;
    double err_asymmetry = err_ratio * asymmetry;

    std::cout<<"forward = "<<n_forward<<" "<<n_forward_mc<<" "<<err_forward<<std::endl;
    std::cout<<"backward = "<<n_backward<<" "<<n_backward_mc<<" "<<err_backward<<std::endl;
    std::cout<<"Asymmetry = "<<asymmetry<<" "<<err_asymmetry<<std::endl;
    std::cout<<"---------------------------------------------------------------"<<std::endl;


    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
    
    gStyle->SetStripDecimals(false);
    gStyle->SetEndErrorSize(4);
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); 
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)

    plotter->data_line_width = 2;

    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 0.8;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;
  
    //printBins(data,var + "_data",false);
    //printBins(mc,var + "_MC" ,false);

    TCanvas* canvas = new TCanvas("c");

    //    void MnvPlotter::DrawDataMCWithErrorBand    (   const MnvH1D *  dataHist,
    //            const MnvH1D *  mcHist,
    //            const Double_t  mcScale = 1.0,
    //            const std::string &     legPos = "L",
    //            const bool  useHistTitles = false,
    //            const MnvH1D *  bkgdHist = NULL,
    //            const MnvH1D *  dataBkgdHist = NULL,
    //            const bool  covAreaNormalize = false,
    //            const bool  statPlusSys = false  
    //            const bool isSmooth = false
    //            )   
    plotter->DrawDataMCWithErrorBand(data, mc, 1.0, "N", true, NULL, NULL, false, true, false);

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(12);
    text.DrawLatex(0.20,0.87,"Area Normalized");
    text.DrawLatex(0.20,0.83,"T_{p} > 0.1 GeV");
    text.DrawLatex(0.20,0.79,"W_{exp} < 1.4 GeV");

    double info_text_x = 0.85;
    double info_text_y = 0.85;
    TLatex info_text;
    info_text.SetNDC();
    info_text.SetTextColor(kBlack);
    info_text.SetTextFont(62);
    info_text.SetTextAlign(12);
    info_text.SetTextSize(0.04);
    info_text.DrawLatex(info_text_x, info_text_y, "(a)");

    // Add Legend
    data->SetMarkerStyle(plotter->data_marker);
    data->SetMarkerSize(plotter->data_marker_size);
    data->SetMarkerColor(plotter->data_color);
    data->SetLineColor(plotter->data_color);
    data->SetLineWidth(plotter->data_line_width);

    mc->SetLineColor(plotter->mc_color);
    mc->SetLineWidth(plotter->mc_line_width);
    mc->SetLineStyle(plotter->mc_line_style);
    mc->SetFillColor(plotter->mc_error_color);
    mc->SetFillStyle(plotter->mc_error_style);

    double leg_x_min = 0.45;
    double leg_x_max = 0.75;
    double leg_y_min = 0.68;
    double leg_y_max = 0.85;

    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    ApplyStyle_Legend(legend);
    legend->SetTextSize(0.05);
    legend->AddEntry(data, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(mc, "Simulation", "fl");
    legend->Draw();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var + "_2Bin" + ".pdf"; 

    canvas->Print(out_name.c_str(),"pdf");
    
    delete f_data;
    delete f_mc;
    delete data;
    delete mc;
    delete legend;
    delete canvas;
    delete plotter;
}

void CCProtonPi0_Plotter::Delta_pi_theta_Fit()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_DeltaRich_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_DeltaRich_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var = "Delta_pi_theta_xsec";

    MnvH1D* data = GetMnvH1D(f_data, var);
    MnvH1D* mc = GetMnvH1D(f_mc, var);
   
    data->GetXaxis()->SetNdivisions(5,4,0);
    mc->GetXaxis()->SetNdivisions(5,4,0);
    
    data->GetYaxis()->SetNdivisions(5,4,0);
    mc->GetYaxis()->SetNdivisions(5,4,0);
 
    data->GetYaxis()->CenterTitle();
    mc->GetYaxis()->CenterTitle();

    data->SetTitle("Data (3.33e20 POT)");
    mc->SetTitle("Simulation");
   
    mc->ClearAllErrorBands();
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
    
    gStyle->SetStripDecimals(false);
    gStyle->SetEndErrorSize(4);
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); 
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)

    plotter->data_line_width = 2;

    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 0.8;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;
  
    //printBins(data,var + "_data",false);
    //printBins(mc,var + "_MC" ,false);

    TCanvas* canvas = new TCanvas("c");

    //    void MnvPlotter::DrawDataMCWithErrorBand    (   const MnvH1D *  dataHist,
    //            const MnvH1D *  mcHist,
    //            const Double_t  mcScale = 1.0,
    //            const std::string &     legPos = "L",
    //            const bool  useHistTitles = false,
    //            const MnvH1D *  bkgdHist = NULL,
    //            const MnvH1D *  dataBkgdHist = NULL,
    //            const bool  covAreaNormalize = false,
    //            const bool  statPlusSys = false  
    //            const bool isSmooth = false
    //            )   
    plotter->DrawDataMCWithErrorBand(data, mc, 1.0, "N", true, NULL, NULL, false, true, false);

    // ------------------------------------------------------------------------
    // Fit
    // ------------------------------------------------------------------------
    const double a = 2.931;
    const double a1 = a + 0.18;
    const double a2 = a - 0.18;

    const double b = 0.43;
    const double b1 = b + 0.577;
    const double b2 = b - 0.577;

    std::vector<double> cos_theta_x;
    std::vector<double> cos_theta_y;
    std::vector<double> cos_theta_y1;
    std::vector<double> cos_theta_y2;
    for (double i = -1.0; i <= 1.0; i = i + 0.01){
        cos_theta_x.push_back(i);

        double y = a + b*std::pow(i,2); 
        double y1 = a1 + b1*std::pow(i,2); 
        double y2 = a2 + b2*std::pow(i,2); 

        cos_theta_y.push_back(y);
        cos_theta_y1.push_back(y1);
        cos_theta_y2.push_back(y2);
    }

    int size = cos_theta_x.size();
    double arr_x[size];
    double arr_y[size];
    double arr_y1[size];
    double arr_y2[size];
    for (int i = 0; i < size; i++){
        arr_x[i] = cos_theta_x[i];
        arr_y[i] = cos_theta_y[i];
        arr_y1[i] = cos_theta_y1[i];
        arr_y2[i] = cos_theta_y2[i];
    }

    TGraph* fit_y = new TGraph(size, arr_x, arr_y);
    fit_y->SetLineColor(kBlue);
    fit_y->SetLineWidth(2);
    fit_y->SetLineStyle(1);
    fit_y->Draw("SAMEL");

    TGraph* fit_y1 = new TGraph(size, arr_x, arr_y1);
    fit_y1->SetLineColor(kBlue);
    fit_y1->SetLineWidth(1);
    fit_y1->SetLineStyle(7);
    fit_y1->Draw("SAMEL");

    TGraph* fit_y2 = new TGraph(size, arr_x, arr_y2);
    fit_y2->SetLineColor(kBlue);
    fit_y2->SetLineWidth(1);
    fit_y2->SetLineStyle(7);
    fit_y2->Draw("SAMEL");

    // Calc ChiSq
    int nBins = data->GetNbinsX();
    double ChiSq_Fit = 0.0; 
    for (int i = 1; i <= nBins; ++i){
        double nData = data->GetBinContent(i);
        double bin = data->GetBinCenter(i);
        double err = data->GetBinError(i);
        double nFit = a + b * std::pow(bin,2);
        double ChiSq = std::pow((nData-nFit),2) / std::pow(err,2);
        //std::cout<<"nData = "<<nData<<" nFit = "<<nFit<<" Err = "<<err<<" Chi = "<<ChiSq<<std::endl; 
        ChiSq_Fit += ChiSq;
    }

    double ChiSq = Calc_ChiSq(data, mc);
      
    std::cout<<"Delta_pi_theta ChiSq = "<<ChiSq<<std::endl;
    std::cout<<"Delta_pi_theta ChiSq/dof = "<<ChiSq/9.0<<std::endl;
    std::cout<<"Delta_pi_theta Fit ChiSq = "<<ChiSq_Fit<<std::endl;
    std::cout<<"Delta_pi_theta Fit ChiSq/dof = "<<ChiSq_Fit/8.0<<std::endl;

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(12);
    text.DrawLatex(0.20,0.87,"POT Normalized");
    text.DrawLatex(0.20,0.83,"T_{p} > 0.1 GeV");
    text.DrawLatex(0.20,0.79,"W_{exp} < 1.4 GeV");

    double info_text_x = 0.85;
    double info_text_y = 0.85;
    TLatex info_text;
    info_text.SetNDC();
    info_text.SetTextColor(kBlack);
    info_text.SetTextFont(62);
    info_text.SetTextAlign(12);
    info_text.SetTextSize(0.04);
    info_text.DrawLatex(info_text_x, info_text_y, "(a)");

    // Add Legend
    data->SetMarkerStyle(plotter->data_marker);
    data->SetMarkerSize(plotter->data_marker_size);
    data->SetMarkerColor(plotter->data_color);
    data->SetLineColor(plotter->data_color);
    data->SetLineWidth(plotter->data_line_width);

    mc->SetLineColor(plotter->mc_color);
    mc->SetLineWidth(plotter->mc_line_width);
    mc->SetLineStyle(plotter->mc_line_style);
    mc->SetFillColor(plotter->mc_error_color);
    mc->SetFillStyle(plotter->mc_error_style);

    double leg_x_min = 0.45;
    double leg_x_max = 0.75;
    double leg_y_min = 0.68;
    double leg_y_max = 0.85;

    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    ApplyStyle_Legend(legend);
    legend->SetTextSize(0.05);
    legend->AddEntry(data, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(mc, "Simulation", "fl");
    legend->AddEntry(fit_y, "Fit to Data", "l");
    legend->Draw();
    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var + "_Fit" + ".pdf"; 

    canvas->Print(out_name.c_str(),"pdf");
    
    delete f_data;
    delete f_mc;
    delete data;
    delete mc;
    delete fit_y;
    delete fit_y1;
    delete fit_y2;
    delete legend;
    delete canvas;
    delete plotter;
}

void CCProtonPi0_Plotter::Delta_pi_phi_Fit()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_DeltaRich_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_DeltaRich_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var = "Delta_pi_phi_xsec";

    MnvH1D* data = GetMnvH1D(f_data, var);
    MnvH1D* mc = GetMnvH1D(f_mc, var);
   
    data->GetXaxis()->SetNdivisions(5,4,0);
    mc->GetXaxis()->SetNdivisions(5,4,0);
    
    data->GetYaxis()->SetNdivisions(5,4,0);
    mc->GetYaxis()->SetNdivisions(5,4,0);
 
    data->GetYaxis()->CenterTitle();
    mc->GetYaxis()->CenterTitle();

    data->SetTitle("Data (3.33e20 POT)");
    mc->SetTitle("Simulation");
   
    mc->ClearAllErrorBands();
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
    
    gStyle->SetStripDecimals(false);
    gStyle->SetEndErrorSize(4);
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); 
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)

    plotter->data_line_width = 2;

    plotter->axis_minimum = 0.001;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 0.8;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;
   
    //printBins(data,var + "_data",false);
    //printBins(mc,var + "_MC" ,false);

    TCanvas* canvas = new TCanvas("c1");

    //    void MnvPlotter::DrawDataMCWithErrorBand    (   const MnvH1D *  dataHist,
    //            const MnvH1D *  mcHist,
    //            const Double_t  mcScale = 1.0,
    //            const std::string &     legPos = "L",
    //            const bool  useHistTitles = false,
    //            const MnvH1D *  bkgdHist = NULL,
    //            const MnvH1D *  dataBkgdHist = NULL,
    //            const bool  covAreaNormalize = false,
    //            const bool  statPlusSys = false  
    //            const bool isSmooth = false
    //            )   
    plotter->DrawDataMCWithErrorBand(data, mc, 1.0, "N", true, NULL, NULL, false, true, false);

    // ------------------------------------------------------------------------
    // Fit
    // ------------------------------------------------------------------------
    const double a = 17.64/1000.0;
    const double a1 = a + 0.43/1000.0;
    const double a2 = a - 0.43/1000.0;

    const double b = -4.81/1000.0;
    const double b1 = b + 1.26/1000.0;
    const double b2 = b - 1.26/1000.0;

    std::vector<double> phi_x;
    std::vector<double> phi_y;
    std::vector<double> phi_y1;
    std::vector<double> phi_y2;
    for (double i = 0.0; i <= 360; ++i){
        phi_x.push_back(i);

        double y = a + b*sin(i * TMath::DegToRad());
        double y1 = a1 + b2*sin(i * TMath::DegToRad()); 
        double y2 = a2 + b1*sin(i * TMath::DegToRad()); 

        phi_y.push_back(y);
        phi_y1.push_back(y1);
        phi_y2.push_back(y2);
    }

    int size = phi_x.size();
    double arr_x[size];
    double arr_y[size];
    double arr_y1[size];
    double arr_y2[size];
    for (int i = 0; i < size; i++){
        arr_x[i] = phi_x[i];
        arr_y[i] = phi_y[i];
        arr_y1[i] = phi_y1[i];
        arr_y2[i] = phi_y2[i];
    }

    TGraph* fit_y = new TGraph(size, arr_x, arr_y);
    fit_y->SetLineColor(kBlue);
    fit_y->SetLineWidth(2);
    fit_y->SetLineStyle(1);
    fit_y->Draw("SAMEL");

    TGraph* fit_y1 = new TGraph(size, arr_x, arr_y1);
    fit_y1->SetLineColor(kBlue);
    fit_y1->SetLineWidth(1);
    fit_y1->SetLineStyle(7);
    fit_y1->Draw("SAMEL");

    TGraph* fit_y2 = new TGraph(size, arr_x, arr_y2);
    fit_y2->SetLineColor(kBlue);
    fit_y2->SetLineWidth(1);
    fit_y2->SetLineStyle(7);
    fit_y2->Draw("SAMEL");

    // Calc ChiSq
    int nBins = data->GetNbinsX();
    double ChiSq_Fit = 0.0; 
    for (int i = 1; i <= nBins; ++i){
        double nData = data->GetBinContent(i);
        double bin = data->GetBinCenter(i);
        double err = data->GetBinError(i);
        double nFit = a + b * sin( bin * TMath::DegToRad());
        double ChiSq = std::pow((nData-nFit),2) / std::pow(err,2);
        //std::cout<<"nData = "<<nData<<" nFit = "<<nFit<<" Err = "<<err<<" Chi = "<<ChiSq<<std::endl; 
        ChiSq_Fit += ChiSq;
    }

    double ChiSq = Calc_ChiSq(data, mc);
        
    std::cout<<"Delta_pi_phi ChiSq = "<<ChiSq<<std::endl;
    std::cout<<"Delta_pi_phi ChiSq/dof = "<<ChiSq/7.0<<std::endl;
    std::cout<<"Delta_pi_phi Fit ChiSq = "<<ChiSq_Fit<<std::endl;
    std::cout<<"Delta_pi_phi Fit ChiSq/dof = "<<ChiSq_Fit/6.0<<std::endl;

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(12);
    text.DrawLatex(0.20,0.87,"POT Normalized");
    text.DrawLatex(0.20,0.83,"T_{p} > 0.1 GeV");
    text.DrawLatex(0.20,0.79,"W_{exp} < 1.4 GeV");

    double info_text_x = 0.85;
    double info_text_y = 0.85;
    TLatex info_text;
    info_text.SetNDC();
    info_text.SetTextColor(kBlack);
    info_text.SetTextFont(62);
    info_text.SetTextAlign(12);
    info_text.SetTextSize(0.04);
    info_text.DrawLatex(info_text_x, info_text_y, "(b)");

    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var + "_Fit" + ".pdf"; 

    canvas->Print(out_name.c_str(),"pdf");

    delete f_data;
    delete f_mc;
    delete data;
    delete mc;
    delete fit_y;
    delete fit_y1;
    delete fit_y2;
    delete canvas;
    delete plotter;
}

void CCProtonPi0_Plotter::Delta_pi_theta_diff_Fit()
{
    std::string root_dir_data = Folder_List::rootDir_Signal_DeltaRich_Data;
    std::string root_dir_mc = Folder_List::rootDir_Signal_DeltaRich_MC;

    TFile* f_data = new TFile(root_dir_data.c_str());
    TFile* f_mc = new TFile(root_dir_mc.c_str());

    std::string var = "Delta_pi_theta_xsec";

    MnvH1D* data = GetMnvH1D(f_data, var);
    MnvH1D* mc = GetMnvH1D(f_mc, var);
   
    data->GetXaxis()->SetNdivisions(5,4,0);
    mc->GetXaxis()->SetNdivisions(5,4,0);
    
    data->GetYaxis()->SetNdivisions(5,4,0);
    mc->GetYaxis()->SetNdivisions(5,4,0);
 
    data->GetYaxis()->CenterTitle();
    mc->GetYaxis()->CenterTitle();

    data->SetTitle("Data (3.33e20 POT)");
    mc->SetTitle("Simulation");
   
    mc->ClearAllErrorBands();
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
    
    gStyle->SetStripDecimals(false);
    gStyle->SetEndErrorSize(4);
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); 
    gStyle->SetOptStat(0); 
    gStyle->SetPadRightMargin(0.05);

    plotter->headroom = 1.75;
    plotter->legend_text_size = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)

    plotter->data_line_width = 2;

    plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 0.8;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;
   
    //printBins(data,var + "_data",false);
    //printBins(mc,var + "_MC" ,false);

    TCanvas* canvas = new TCanvas("c1");

    //    void MnvPlotter::DrawDataMCWithErrorBand    (   const MnvH1D *  dataHist,
    //            const MnvH1D *  mcHist,
    //            const Double_t  mcScale = 1.0,
    //            const std::string &     legPos = "L",
    //            const bool  useHistTitles = false,
    //            const MnvH1D *  bkgdHist = NULL,
    //            const MnvH1D *  dataBkgdHist = NULL,
    //            const bool  covAreaNormalize = false,
    //            const bool  statPlusSys = false  
    //            const bool isSmooth = false
    //            )   
    plotter->DrawDataMCWithErrorBand(data, mc, 1.0, "N", true, NULL, NULL, false, true, false);

    // ------------------------------------------------------------------------
    // Fit
    // ------------------------------------------------------------------------
    const double a = 1.7982;
    const double a1 = 1.5671;
    const double a2 = 2.0293;

    const double b = 0.5460;
    const double b1 = 0.4132;
    const double b2 = 0.6788;

    std::vector<double> theta_x;
    std::vector<double> theta_y;
    std::vector<double> theta_y1;
    std::vector<double> theta_y2;
    for (double i = -1.0; i <= 1.0; i = i + 0.01){
        theta_x.push_back(i);

        double y = a*i + b;
        double y1 = a1*i + b1;
        double y2 = a2*i + b2;

        theta_y.push_back(y);
        theta_y1.push_back(y1);
        theta_y2.push_back(y2);
    }

    int size = theta_x.size();
    double arr_x[size];
    double arr_y[size];
    double arr_y1[size];
    double arr_y2[size];
    for (int i = 0; i < size; i++){
        arr_x[i] = theta_x[i];
        arr_y[i] = theta_y[i];
        arr_y1[i] = theta_y1[i];
        arr_y2[i] = theta_y2[i];
    }

    TGraph* fit_y = new TGraph(size, arr_x, arr_y);
    fit_y->SetLineColor(kBlue);
    fit_y->SetLineWidth(2);
    fit_y->SetLineStyle(1);
    fit_y->Draw("SAMEL");

    TGraph* fit_y1 = new TGraph(size, arr_x, arr_y1);
    fit_y1->SetLineColor(kBlue);
    fit_y1->SetLineWidth(1);
    fit_y1->SetLineStyle(7);
    fit_y1->Draw("SAMEL");

    TGraph* fit_y2 = new TGraph(size, arr_x, arr_y2);
    fit_y2->SetLineColor(kBlue);
    fit_y2->SetLineWidth(1);
    fit_y2->SetLineStyle(7);
    fit_y2->Draw("SAMEL");

    // Calc ChiSq
    int nBins = data->GetNbinsX();
    double ChiSq_Fit = 0.0; 
    for (int i = 1; i <= nBins; ++i){
        double nData = data->GetBinContent(i);
        double bin = data->GetBinCenter(i);
        double err = data->GetBinError(i);
        double nFit = a*bin + b;
        double ChiSq = std::pow((nData-nFit),2) / std::pow(err,2);
        //std::cout<<"nData = "<<nData<<" nFit = "<<nFit<<" Err = "<<err<<" Chi = "<<ChiSq<<std::endl; 
        ChiSq_Fit += ChiSq;
    }

    double ChiSq = Calc_ChiSq(data, mc);
        
    std::cout<<"Delta_pi_theta ChiSq = "<<ChiSq<<std::endl;
    std::cout<<"Delta_pi_theta ChiSq/dof = "<<ChiSq/7.0<<std::endl;
    std::cout<<"Delta_pi_theta Fit ChiSq = "<<ChiSq_Fit<<std::endl;
    std::cout<<"Delta_pi_theta Fit ChiSq/dof = "<<ChiSq_Fit/6.0<<std::endl;

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(12);
    text.DrawLatex(0.20,0.87,"POT Normalized");
    text.DrawLatex(0.20,0.83,"T_{p} > 0.1 GeV");
    text.DrawLatex(0.20,0.79,"W_{exp} < 1.4 GeV");

    std::string plotDir = Folder_List::plotDir_Paper;
    std::string out_name = plotDir + var + "_Diff_Fit" + ".pdf"; 

    canvas->Print(out_name.c_str(),"pdf");

    delete f_data;
    delete f_mc;
    delete data;
    delete mc;
    delete fit_y;
    delete fit_y1;
    delete fit_y2;
    delete canvas;
    delete plotter;
}
#endif

