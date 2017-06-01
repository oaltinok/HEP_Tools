/*
   See CCProtonPi0_Plotter.h header for Class Information
*/

#ifndef CCProtonPi0_Plotter_Macros_cpp
#define CCProtonPi0_Plotter_Macros_cpp

#include "CCProtonPi0_Plotter.h"

using namespace PlotUtils;

void CCProtonPi0_Plotter::DrawStackedMC(rootDir &dir, std::string var_name, std::string plotDir, int nCutArrows, CutArrow cutArrow1, CutArrow cutArrow2)
{
    DrawStackedMC_BckgAll(dir, var_name, plotDir, nCutArrows, cutArrow1, cutArrow2);
    //DrawStackedMC_BckgType(dir, var_name, plotDir, nCutArrows, cutArrow1, cutArrow2);
}

void CCProtonPi0_Plotter::DrawDataStackedMC(rootDir &dir, std::string var_name, std::string plotDir, int nCutArrows, CutArrow cutArrow1, CutArrow cutArrow2)
{
    // ------------------------------------------------------------------------
    // POT Normalized Plots 
    // ------------------------------------------------------------------------
    DrawDataStackedMC_BckgAll(dir, var_name,plotDir, true, nCutArrows, cutArrow1, cutArrow2);
    //DrawDataStackedMC_BckgType(dir, var_name,plotDir, true, nCutArrows, cutArrow1, cutArrow2);

    // ------------------------------------------------------------------------
    // Area Normalized Plots 
    // ------------------------------------------------------------------------
    //DrawDataStackedMC_BckgAll(dir, var_name,plotDir, false, nCutArrows, cutArrow1, cutArrow2);
    //DrawDataStackedMC_BckgType(dir, var_name,plotDir, false, nCutArrows, cutArrow1, cutArrow2);
}

void CCProtonPi0_Plotter::DrawDataStackedMC_BckgAll(rootDir &dir, std::string var_name, std::string plotDir, bool isPOTNorm, int nCutArrows, CutArrow cutArrow1, CutArrow cutArrow2)
{
    std::string rootDir_data = dir.data;
    std::string rootDir_mc = dir.mc;

    TFile* f_data = new TFile(rootDir_data.c_str());
    TFile* f_mc = new TFile(rootDir_mc.c_str());

    std::string var;
    double hist_max = 0;
    double max_bin;
    double bin_width;
    double area_mc = 0;

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
    max_bin = data->GetMaximumBin();
    hist_max = (hist_max + data->GetBinContent(max_bin))*1.2;
    bin_width = data->GetBinWidth(1);

    //data->GetXaxis()->SetTitle("Log Likelihood Ratio(LLR) for Proton");

    // ------------------------------------------------------------------------
    // MC Normalization 
    // ------------------------------------------------------------------------
    var = Form("%s_%d",var_name.c_str(),0);
    MnvH1D* mc_all = (MnvH1D*)f_mc->Get(var.c_str()); 
    std::string norm_label; 
    double mc_ratio = GetMCNormalization(norm_label, isPOTNorm, data, mc_all);

    // ------------------------------------------------------------------------
    // Fill TObjArray - For MC Histograms
    // ------------------------------------------------------------------------
    TObjArray* mc_hists = new TObjArray;

    // Get All Background
    var = Form("%s_%d",var_name.c_str(),2);
    MnvH1D* background = (MnvH1D*)f_mc->Get(var.c_str());
    background->SetTitle("Background");
    background->SetLineColor(kRed);
    background->SetFillColor(kRed);
    mc_hists->Add(background);
    area_mc += background->Integral() * mc_ratio;

    std::cout<<"Background = "<<background->Integral("width")<<std::endl;
    // Get Signal
    var = Form("%s_%d",var_name.c_str(),1);
    MnvH1D* signal = (MnvH1D*)f_mc->Get(var.c_str());
    signal->SetTitle("Signal");
    signal->SetLineColor(kGreen);
    signal->SetFillColor(kGreen);
    mc_hists->Add(signal);
    area_mc += signal->Integral() * mc_ratio;

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    ApplyStyle(plotter);
    TCanvas* c = new TCanvas("c","c",1280,800);

    plotter->headroom = 1.75;
    plotter->data_line_width = 2;
    plotter->data_marker_size = 1.5;
    gStyle->SetEndErrorSize(6);

    plotter->DrawDataStackedMC(data, mc_hists,mc_ratio ,"N","Data (3.33e20 POT)",0);

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");
 
    // If Cut Histogram - Add Cut Arrows
    if (nCutArrows == 1){
        AddCutArrow(plotter, cutArrow1, hist_max, bin_width );
    }else if (nCutArrows == 2){
        AddCutArrow(plotter, cutArrow1, hist_max, bin_width);
        AddCutArrow(plotter, cutArrow2, hist_max, bin_width);
    }

    // Add Pi0 InvMass Line
    std::size_t found = var_name.find("invMass");
    if (found != std::string::npos){
        TArrow arrow;
        arrow.SetLineWidth(2);
        arrow.SetLineColor(kBlue);
        arrow.DrawArrow(134.98,hist_max/8,134.98,0, 0.01); 
    }

    // Add Delta+ InvMass Line
    if (    var_name.compare("deltaInvMass") == 0 ||
            var_name.compare("W_Calc") == 0 ){
        TLine deltaMass;
        deltaMass.SetLineColor(kBlue);
        deltaMass.DrawLine(1.232,0,1.232,hist_max);
    }

    // ------------------------------------------------------------------------
    // Legend 
    // ------------------------------------------------------------------------
    double leg_x_min = 0.60;
    double leg_x_max = 0.90;
    double leg_y_min = 0.60;
    double leg_y_max = 0.85;

    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    ApplyStyle_Legend(legend);
    //legend->SetTextSize(0.04);
    legend->AddEntry(data, "Data #scale[0.85]{(3.33e20 POT)}", "lep" );
    legend->AddEntry(signal, "Signal","f");
    legend->AddEntry(background, "Background","f");
    legend->Draw();

 
    // Print Plot
    c->Print(Form("%s%s%s%s%s",plotDir.c_str(),var_name.c_str(),"_bckg_all_",norm_label.c_str(),".pdf"), "pdf");

    delete c;
    delete plotter;
}

void CCProtonPi0_Plotter::DrawDataStackedMC_BckgType(rootDir &dir, std::string var_name, std::string plotDir, bool isPOTNorm, int nCutArrows, CutArrow cutArrow1, CutArrow cutArrow2)
{
    std::string rootDir_data = dir.data;
    std::string rootDir_mc = dir.mc;

    TFile* f_data = new TFile(rootDir_data.c_str());
    TFile* f_mc = new TFile(rootDir_mc.c_str());

    std::string var;
    double hist_max = 0;
    double max_bin;
    double bin_width;

    // ------------------------------------------------------------------------
    // Get Data Histogram
    // ------------------------------------------------------------------------
    var = Form("%s_%d",var_name.c_str(),0);
    MnvH1D* data = (MnvH1D*)f_data->Get(var.c_str()); 

    // ------------------------------------------------------------------------
    // MC Normalization 
    // ------------------------------------------------------------------------
    var = Form("%s_%d",var_name.c_str(),0);
    MnvH1D* mc_all = (MnvH1D*)f_mc->Get(var.c_str()); 
    std::string norm_label; 
    double mc_ratio = GetMCNormalization(norm_label, isPOTNorm, data, mc_all);

    // ------------------------------------------------------------------------
    // Fill TObjArray - For MC Histograms
    // ------------------------------------------------------------------------
    TObjArray* mc_hists = new TObjArray;
    double area_mc = FormTObjArray_BckgType(f_mc, var_name, mc_hists, hist_max, bin_width);

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas("c","c",1280,800);
    ApplyStyle(plotter);
    //plotter->axis_minimum = 0.1;
    plotter->axis_minimum = 0.0;
    plotter->DrawDataStackedMC(data, mc_hists,mc_ratio ,"TR","Data",0);

    // Add Plot Labels
    plotter->AddHistoTitle(data[0].GetTitle());
    AddNormBox(plotter, isPOTNorm, mc_ratio);

    // Get Line Heights from Data
    max_bin = data->GetMaximumBin();
    hist_max = data->GetBinContent(max_bin)*1.2;
    bin_width = data->GetBinWidth(1);

    // If Cut Histogram - Add Cut Arrows
    if (nCutArrows == 1){
        AddCutArrow(plotter, cutArrow1, hist_max, bin_width);
    }else if (nCutArrows == 2){
        AddCutArrow(plotter, cutArrow1, hist_max, bin_width);
        AddCutArrow(plotter, cutArrow2, hist_max, bin_width);
    }

    // Add Pi0 InvMass Line
    std::size_t found = var_name.find("invMass");
    if (found != std::string::npos){
        TLine pi0Mass;
        pi0Mass.SetLineWidth(2);
        pi0Mass.SetLineColor(kBlue);
        pi0Mass.DrawLine(134.98,0,134.98,hist_max);
    }

    // Add Delta+ InvMass Line
    if (    var_name.compare("deltaInvMass") == 0 ||
            var_name.compare("W_Calc") == 0 ){
        TLine deltaMass;
        deltaMass.SetLineWidth(2);
        deltaMass.SetLineColor(kBlue);
        deltaMass.DrawLine(1.232,0,1.232,hist_max);
    }

    // Add Areas
    const double y_pos = 0.88;
    const double text_size = 0.03;
    double area_data = data->Integral();
    area_mc = area_mc * mc_ratio;
    plotter->AddPlotLabel(Form("Area(Data)/Area(MC) = %3.2f",area_data/area_mc),0.3,y_pos,text_size,kBlue); 


    // Print Plot
    c->Print(Form("%s%s%s%s%s",plotDir.c_str(),var_name.c_str(),"_bckg_Type_",norm_label.c_str(),".pdf"), "pdf");

    delete c;
    delete plotter;
}

void CCProtonPi0_Plotter::ApplyStyle(MnvPlotter* plotter)
{
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); // 4x3 aspect ratio
    gStyle->SetPadRightMargin(0.05);

    gStyle->SetStripDecimals(false);
    plotter->legend_text_size  = 0.045;
    plotter->legend_text_font = 42; // default 62 (bold)
    //plotter->data_marker_size = 1.3;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 1.1;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;

    gStyle->SetEndErrorSize(4);

    //gStyle->SetAxisColor(kGray + 1,  "XY");
    //gStyle->SetLabelColor(kGray + 1, "XY");
    //gStyle->SetFrameLineColor(kGray + 1);
    
//    plotter->hist_min_zero = false;
//    plotter->mc_line_width = 0;
//    plotter->axis_title_size_x = 0.04;
//    plotter->axis_title_size_y = 0.04;
//    plotter->axis_label_size = 0.03;
//    plotter->legend_text_size = 0.02;

}

void CCProtonPi0_Plotter::ApplyStyle_Errors(MnvPlotter* plotter, bool groupErrors)
{
    //plotter->hist_min_zero = false;
    //plotter->mc_line_width = 0;
    //plotter->axis_title_size_x = 0.04;
    //plotter->axis_title_size_y = 0.04;
    //plotter->axis_label_size = 0.03;
    plotter->legend_text_size = 0.04;
    plotter->height_nspaces_per_hist = 0.9;
    plotter->width_xspace_per_letter = 0.2;
    plotter->legend_offset_x = -0.2;
    plotter->legend_offset_y = 0.0;

    if (groupErrors){
        //-- define colors of the standard errors
        plotter->error_color_map.clear();
        plotter->error_summary_group_map.clear();

        //plotter->error_color_map["GENIE"] = kGreen+2;

        plotter->error_summary_group_map["Detector"] = detGroup;
        plotter->error_summary_group_map["X-Sec Model"] = genieGroup;   
        plotter->error_summary_group_map["FSI Model"] = fsiGroup;   
        plotter->error_summary_group_map["Flux"] = fluxGroup;
        plotter->error_summary_group_map["Other"] = otherGroup;
    }
}

void CCProtonPi0_Plotter::DrawErrorSummary_PaperStyle(MnvH1D* hist, std::string var_name, std::string plotDir, bool groupErrors)
{
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
    ApplyStyle(plotter);

    TCanvas* c = new TCanvas("c");

    ApplyStyle_Errors(plotter, groupErrors);
    //ApplyStyle_Errors(plotter, false);
 
    plotter->legend_n_columns = 2;
    plotter->axis_maximum = 0.5;

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
    plotter->DrawErrorSummary(hist, "N", true, true, 0.0);

    //// Add Normalization Labels
    TPaveText* text = new TPaveText(0.16, 0.80, 0.51, 0.90, "NDC");
    text->AddText("#scale[1.15]{c)  #nu_{#mu} + CH #rightarrow #mu^{-} + #pi^{0} + X}");
    text->AddText("#scale[0.85]{POT Normalized}");
    text->SetTextColor(kBlue);
    text->SetFillColor(kWhite);
    text->SetTextFont(42);
    text->Draw("SAME");

    // Add Normalization Labels
    //TLatex text;
    //text.SetNDC();
    //text.SetTextColor(kBlue);
    //text.SetTextSize(0.03);
    //text.SetTextAlign(22);
    //text.DrawLatex(0.30,0.87,"POT Normalized");


    // ------------------------------------------------------------------------
    // Legend 
    // ------------------------------------------------------------------------
    TH1D* h_total = GetBinNormalizedTH1D(hist);
    h_total->SetLineWidth(2);
    h_total->SetLineColor(kBlack);
    h_total->SetLineStyle(kSolid);

    TH1D* h_stat = GetBinNormalizedTH1D(hist);
    h_stat->SetLineWidth(2);
    h_stat->SetLineColor(kBlack);
    h_stat->SetLineStyle(kDashed);

    TH1D* h_detector = GetBinNormalizedTH1D(hist);
    h_detector->SetLineWidth(2); 
    h_detector->SetLineColor(kBlue+2);
    h_detector->SetLineStyle(3);
 
    TH1D* h_xsec = GetBinNormalizedTH1D(hist);
    h_xsec->SetLineWidth(2); 
    h_xsec->SetLineColor(kRed+2);
    h_xsec->SetLineStyle(5);       
  
    TH1D* h_fsi = GetBinNormalizedTH1D(hist);
    h_fsi->SetLineWidth(2); 
    h_fsi->SetLineColor(kMagenta+2);
    h_fsi->SetLineStyle(9);       
   
    TH1D* h_flux = GetBinNormalizedTH1D(hist);
    h_flux->SetLineWidth(2); 
    h_flux->SetLineColor(TColor::GetColor("#8b4513"));
    h_flux->SetLineStyle(6);       
    
    TH1D* h_other = GetBinNormalizedTH1D(hist);
    h_other->SetLineWidth(2); 
    h_other->SetLineColor(kGreen+3);
    h_other->SetLineStyle(7);       
   
    TLegend *legend = new TLegend(0.52,0.74,0.94,0.90);  
    ApplyStyle_Legend(legend);
    legend->SetNColumns(2);
    legend->SetTextSize(0.038);
    legend->AddEntry(h_total, "Total Error", "l" );
    legend->AddEntry(h_stat, "Statistical", "l" );
    legend->SetColumnSeparation(0.06);
    legend->SetEntrySeparation(0.12);
    legend->AddEntry(h_detector, "Detector", "l");
    legend->AddEntry(h_xsec, "X-Sec Model", "l");
    legend->AddEntry(h_fsi, "FSI Model", "l");
    legend->AddEntry(h_flux, "Flux", "l");
    legend->AddEntry(h_other, "Other", "l");
    //legend->Draw();

    // Print Plot
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

void CCProtonPi0_Plotter::DrawErrorSummary(MnvH1D* hist, std::string var_name, std::string plotDir, bool groupErrors)
{
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

    plotter->headroom = 1.75;
    plotter->legend_text_size = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)
    
    plotter->data_line_width = 2;

    plotter->axis_minimum = 0.01;
    plotter->axis_maximum = 0.5;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 0.8;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;
 
    ApplyStyle_Errors(plotter, groupErrors);

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
    
    TLegend *legend = new TLegend(0.35,0.74,0.85,0.90);  
    ApplyStyle_Legend(legend);
    legend->SetNColumns(2);
    legend->SetTextSize(0.04);
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

    // Print Plot
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

void CCProtonPi0_Plotter::DrawErrorBand(MnvH1D* hist, std::string error_name, int error_color, std::string var_name, std::string plotDir)
{
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas("c","c",1280,800);

    //ApplyStyle_Error(plotter);

    TH1D* error = hist->GetVertErrorBand(error_name.c_str());
    plotter->DrawErrorBand(error, error_color);

    // Add Plot Labels
    plotter->AddHistoTitle(hist->GetTitle());

    // Print Plot
    c->Print(Form("%s%s%s%s%s",plotDir.c_str(),var_name.c_str(),"_Error_", error_name.c_str() ,".pdf"), "pdf");

    delete c;
    delete plotter;
}

void CCProtonPi0_Plotter::DrawDataMC_Unfolded(rootDir& dir, std::string var_name, std::string plotDir, double nBckg)
{
    std::string rootDir_mc = dir.mc;
    std::string rootDir_data = dir.data;

    TFile* f_mc = new TFile(rootDir_mc.c_str());

    MnvH1D* data = GetBckgSubtractedData(dir, var_name, nBckg);
    MnvH1D* mc = GetMnvH1D(f_mc, Form("%s_%d",var_name.c_str(),1));
    DrawDataMC_Thesis(data, mc, var_name, plotDir, false);

    std::string var = var_name + "_response";
    MnvH2D* response = GetMnvH2D(f_mc, var);

    MnvH1D* data_unfolded = 0;
    MnvH1D* mc_unfolded = 0;

    // Use MnvUnfold to Unfold 
    MinervaUnfold::MnvUnfold::Get().UnfoldHisto(data_unfolded, response, data, RooUnfold::kBayes, 4, true);
    MinervaUnfold::MnvUnfold::Get().UnfoldHisto(mc_unfolded, response, mc, RooUnfold::kBayes, 4, true);

    var_name = var_name + "_unfolded";
    DrawDataMC_Thesis(data_unfolded, mc_unfolded, var_name, plotDir, false);

    delete data;
    delete mc;
    delete response;
    delete f_mc;
}

void CCProtonPi0_Plotter::DrawDataMC_EffCorrected(rootDir& dir, std::string var_name, std::string plotDir, double nBckg)
{
    std::string rootDir_mc = dir.mc;
    std::string rootDir_data = dir.data;

    TFile* f_mc = new TFile(rootDir_mc.c_str());

    // Background Subtraction
    MnvH1D* data = GetBckgSubtractedData(dir, var_name, nBckg);

    // Unfolding
    std::string var = var_name + "_response";
    MnvH2D* response = GetMnvH2D(f_mc, var);

    MnvH1D* data_unfolded = 0;

    MinervaUnfold::MnvUnfold::Get().UnfoldHisto(data_unfolded, response, data, RooUnfold::kBayes, 4, true);

    // Efficiency Correction
    MnvH1D* mc = GetMnvH1D(f_mc, Form("%s_%d",var_name.c_str(),1));

    TFile* f_truth = new TFile(rootDir_Truth.mc.c_str()); 
    var = var_name + "_all_signal";
    MnvH1D* mc_all_signal = GetMnvH1D(f_truth, var);

    MnvH1D* eff_curve = new MnvH1D(*mc);
    eff_curve->Divide(mc, mc_all_signal);

    MnvH1D* data_eff_corrected = new MnvH1D(*data_unfolded);

    data_eff_corrected->Divide(data_unfolded, eff_curve);

    var_name = var_name + "_eff_corrected";
    DrawDataMC_Thesis(data_eff_corrected, mc_all_signal, var_name, plotDir, false);

    printBins(data_eff_corrected, var_name, false);

    delete data;
    delete mc;
    delete response;
    delete data_unfolded;
    delete mc_all_signal;
    delete eff_curve;
    delete data_eff_corrected;
    delete f_truth;
    delete f_mc;
}

void CCProtonPi0_Plotter::DrawDataMC_EffCorrected_Stacked(rootDir& dir, std::string var_name, std::string plotDir, double nBckg)
{
    std::string rootDir_mc = dir.mc;
    std::string rootDir_data = dir.data;

    TFile* f_mc = new TFile(rootDir_mc.c_str());

    // Background Subtraction
    MnvH1D* data = GetBckgSubtractedData(dir, var_name, nBckg);

    // Unfolding
    std::string var = var_name + "_response";
    MnvH2D* response = GetMnvH2D(f_mc, var);

    MnvH1D* data_unfolded = 0;

    MinervaUnfold::MnvUnfold::Get().UnfoldHisto(data_unfolded, response, data, RooUnfold::kBayes, 4, true);

    // ------------------------------------------------------------------------
    // Efficiency Correction
    // ------------------------------------------------------------------------
    MnvH1D* mc = GetMnvH1D(f_mc, Form("%s_%s",var_name.c_str(),"mc_reco_signal"));

    TFile* f_truth = new TFile(rootDir_Truth.mc.c_str()); 
    var = var_name + "_all_signal";
    MnvH1D* mc_all_signal = GetMnvH1D(f_truth, var);

    MnvH1D* eff_curve = new MnvH1D(*mc);
    eff_curve->Divide(mc, mc_all_signal);

    MnvH1D* data_eff_corrected = new MnvH1D(*data_unfolded);

    data_eff_corrected->Divide(data_unfolded, eff_curve);

    std::vector<MnvH1D*> mc_IntType;
    var = var_name + "_delta_res";
    //var = var_name + "_delta_res_BeforeFSI";
    MnvH1D* mc_delta_res = GetMnvH1D(f_truth, var);
    mc_IntType.push_back(mc_delta_res);

    var = var_name + "_other_res";
    //var = var_name + "_other_res_BeforeFSI";
    MnvH1D* mc_other_res = GetMnvH1D(f_truth, var);
    mc_IntType.push_back(mc_other_res);

    //var = var_name + "_non_res";
    var = var_name + "_non_res_BeforeFSI";
    MnvH1D* mc_non_res = GetMnvH1D(f_truth, var);
    mc_IntType.push_back(mc_non_res);

    var_name = var_name + "_eff_corrected_stacked";
    DrawDataMC_IntType(data_eff_corrected, mc_all_signal, mc_IntType, var_name, plotDir);

    delete data;
    delete mc;
    delete response;
    delete data_unfolded;
    delete mc_all_signal;
    delete eff_curve;
    delete data_eff_corrected;
    delete f_truth;
    delete f_mc;
}
void CCProtonPi0_Plotter::DrawDataMC_Thesis(rootDir& dir, std::string var_name, std::string plotDir)
{
    std::string rootDir_mc = dir.mc;
    std::string rootDir_data = dir.data;

    TFile* f_mc = new TFile(rootDir_mc.c_str());
    TFile* f_data = new TFile(rootDir_data.c_str());

    std::string var = Form("%s_%d",var_name.c_str(),0);
    //std::string var = var_name;

    MnvH1D* mc = (MnvH1D*)f_mc->Get(var.c_str());
    MnvH1D* data = (MnvH1D*)f_data->Get(var.c_str()); 
    DrawDataMC_Thesis(data, mc, var_name, plotDir);
}

void CCProtonPi0_Plotter::DrawDataMC_WithBackground(std::string var_name, std::string plotDir)
{
    std::cout<<"Plotting for "<<var_name<<std::endl;

    TFile* f_xsec_data = new TFile(rootDir_CrossSection.data.c_str());
    TFile* f_xsec_mc = new TFile(rootDir_CrossSection.mc.c_str());

    std::string data_var = var_name + "_all";
    std::string mc_var = var_name + "_mc_reco_all";
    std::string bckg_var = var_name + "_mc_reco_bckg";

    MnvH1D* data = GetMnvH1D(f_xsec_data, data_var);
    MnvH1D* mc = GetMnvH1D(f_xsec_mc, mc_var);
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
    // Plot Data vs MC 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
       
    gStyle->SetStripDecimals(false);
    gStyle->SetEndErrorSize(4);
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); 
    gStyle->SetOptStat(0); 

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
     
    if (var_name.compare("pi0_theta_xsec_data_MC") == 0){
        plotter->axis_title_offset_y = 1.2;
    }
 
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
    plotter->DrawDataMCWithErrorBand(data, mc, POT_ratio, "N", true, NULL, NULL, false, true, false);

   
    // ------------------------------------------------------------------------
    // Add Background Histogram
    // ------------------------------------------------------------------------
    MnvH1D* mnv_bckg = GetMnvH1D(f_xsec_mc, bckg_var);
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
    bckg->GetXaxis()->SetNdivisions(5,5,0);
    bckg->GetYaxis()->SetNdivisions(5,5,0);
    bckg->GetYaxis()->CenterTitle();
    bckg->Draw("HIST SAME");

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
    double leg_y_min = 0.65;
    double leg_y_max = 0.90;

    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    ApplyStyle_Legend(legend);

    legend->AddEntry(data, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(mc, "Simulation", "fl");
    legend->AddEntry(bckg, "Background", "f");
    legend->Draw();

    // Plot Output
    std::string out_name = plotDir + var_name + "_POT" + ".pdf"; 
    c->Print(out_name.c_str(),"pdf");

    delete c;
    delete plotter;
    delete data;
    delete mc;
    delete bckg;
    delete mnv_bckg;
}

void CCProtonPi0_Plotter::DrawDataMC_Thesis(MnvH1D* data, MnvH1D* mc, std::string var_name, std::string plotDir, bool isXSec)
{
    std::cout<<"Plotting for "<<var_name<<std::endl;

    double mc_ratio = isXSec ? 1.0 : POT_ratio;

    MnvH1D* tempData = new MnvH1D(*data);
    MnvH1D* tempMC = new MnvH1D(*mc);
    
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
    gStyle->SetCanvasDefW(640);
    gStyle->SetCanvasDefH(480); 
    gStyle->SetOptStat(0); 

    plotter->headroom = 1.75;
    plotter->legend_text_size = 0.04;
    plotter->legend_text_font = 42; // default 62 (bold)

    plotter->data_line_width = 2;

    //plotter->axis_minimum = 0.01;
    plotter->axis_title_font_x   = 42;
    plotter->axis_title_size_x   = 0.06;
    plotter->axis_title_offset_x = 1.1;
    plotter->axis_title_font_y   = 42;
    plotter->axis_title_size_y   = 0.06;
    plotter->axis_title_offset_y = 0.8;
    plotter->axis_label_size = 0.05;
    plotter->axis_label_font = 42;
     
    if (var_name.compare("pi0_theta_xsec_data_MC") == 0){
        plotter->axis_title_offset_y = 1.2;
    }
 
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
    plotter->DrawDataMCWithErrorBand(tempData, tempMC, mc_ratio, "N", true, NULL, NULL, false, true, false);

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");

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
    double leg_x_max = 0.80;
    double leg_y_min = 0.70;
    double leg_y_max = 0.90;

    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    ApplyStyle_Legend(legend);

    legend->AddEntry(tempData, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(tempMC, "Simulation", "fl");
    legend->Draw();

    // Plot Output
    std::string out_name = plotDir + var_name + "_POT" + ".pdf"; 
    c->Print(out_name.c_str(),"pdf");

    delete c;
    delete plotter;
    delete tempData;
    delete tempMC;
}

void CCProtonPi0_Plotter::DrawDataMC(rootDir& dir, std::string var_name, std::string plotDir)
{
    std::string rootDir_mc = dir.mc;
    std::string rootDir_data = dir.data;

    TFile* f_mc = new TFile(rootDir_mc.c_str());
    TFile* f_data = new TFile(rootDir_data.c_str());

    //std::string var = Form("%s_%d",var_name.c_str(),0);
    std::string var = var_name;

    MnvH1D* mc = (MnvH1D*)f_mc->Get(var.c_str());
    MnvH1D* data = (MnvH1D*)f_data->Get(var.c_str()); 
    DrawDataMC(data, mc, var_name, plotDir, false);
}

void CCProtonPi0_Plotter::DrawDataMC_Signal(rootDir& dir, std::string var_name, std::string plotDir, double nBckg)
{
    std::string rootDir_mc = dir.mc;
    std::string rootDir_data = dir.data;

    TFile* f_mc = new TFile(rootDir_mc.c_str());

    std::string var = Form("%s_%d",var_name.c_str(),1);

    MnvH1D* data = GetBckgSubtractedData(dir, var_name, nBckg);
    MnvH1D* mc = GetMnvH1D(f_mc, var);
    DrawDataMC_Thesis(data, mc, var_name, plotDir, false);

    delete data;
    delete mc;
    delete f_mc;
}

void CCProtonPi0_Plotter::DrawDataMC(MnvH1D* data, MnvH1D* mc, std::string var_name, std::string plotDir, bool isXSec)
{
    // POT Normalized
    DrawDataMC_WithRatio(data, mc, var_name, plotDir, true, isXSec);
    // Area Normalized
    DrawDataMC_WithRatio(data, mc, var_name, plotDir, false, isXSec);
}

void CCProtonPi0_Plotter::DrawDataMC_WithRatio(MnvH1D* data, MnvH1D* mc, std::string var_name, std::string plotDir, bool isPOTNorm, bool isXSec)
{
    std::cout<<"Plotting for "<<var_name<<std::endl;

    // ------------------------------------------------------------------------
    // MC Normalization 
    // ------------------------------------------------------------------------
    std::string norm_label; 
    double mc_ratio = GetMCNormalization(norm_label, isPOTNorm, data, mc);

    // Cross Section Calculation Already Includes POT Normalization
    if (isXSec && isPOTNorm) mc_ratio = 1.0;

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas("c","c",1280,1280);

    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0); 
    pad1->SetBottomMargin(0); // Top and Bottom Plots attached
    pad1->Draw();               
    pad1->cd(); // pad1 is the current pad

    MnvH1D* tempData = new MnvH1D(*data);
    MnvH1D* tempMC = new MnvH1D(*mc);

    // ------------------------------------------------------------------------
    // Neutrino Energy Comparison Only
    //tempData->SetMaximum(80);
    //tempMC->SetMaximum(80);
    //tempData->GetXaxis()->SetRangeUser(0,10);
    //tempMC->GetXaxis()->SetRangeUser(0,10);
    // ------------------------------------------------------------------------


    //    void MnvPlotter::DrawDataMCWithErrorBand    (   const MnvH1D *  dataHist,
    //            const MnvH1D *  mcHist,
    //            const Double_t  mcScale = 1.0,
    //            const std::string &     legPos = "L",
    //            const bool  useHistTitles = false,
    //            const MnvH1D *  bkgdHist = NULL,
    //            const MnvH1D *  dataBkgdHist = NULL,
    //            const bool  covAreaNormalize = false,
    //            const bool  statPlusSys = false  
    //            )   
    plotter->headroom = 1.75;
    plotter->data_line_width = 2;
    plotter->data_marker_size = 1.5;
    gStyle->SetEndErrorSize(6);
    plotter->DrawDataMCWithErrorBand(tempData, tempMC, mc_ratio, "TR", false, NULL, NULL, false, true);

    // ------------------------------------------------------------------------
    // Plot Labels 
    // ------------------------------------------------------------------------
    // Add ChiSq Text 
    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.03);
//
//    int ndf;
//    double chiSq = plotter->Chi2DataMC( tempData, tempMC, ndf, mc_ratio);
//    char *chiSq_text = Form("#chi^{2}/ndf = %3.2f/%d = %3.2f", chiSq, ndf, chiSq/(double)ndf); 
//
//    text.DrawLatex(0.20,0.87,chiSq_text);

    tempMC->Scale(mc_ratio);
    double area_data = tempData->Integral("width");
    double area_mc = tempMC->Integral("width");
    char *area_text = Form("Area(Data/MC) = %3.2f", area_data/area_mc); 
    text.DrawLatex(0.20,0.87,area_text);

    // Add Error Explanation Text
    std::string data_err_text = "Data: inner errors statistical";
    std::string mc_err_text = "Simulation: statistical errors only";
    text.DrawLatex(0.55, 0.70, data_err_text.c_str());
    text.DrawLatex(0.55, 0.66, mc_err_text.c_str());

    // Add Normalization Labels
    plotter->AddHistoTitle(data->GetTitle());
    AddNormBox(plotter, isPOTNorm, mc_ratio);

    // Plot Lower Plot: Data vs MC Ratio
    c->cd(); // Go back to default Canvas before creating 2nd Pad
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.2);
    pad2->Draw();
    pad2->cd();

    // Calculate the Ratio 
    TH1D* h_data = new TH1D(tempData->GetCVHistoWithStatError());
    TH1D* h_mc_total = new TH1D(tempMC->GetCVHistoWithStatError());
    // Scale Histograms
    //h_mc_total->Scale(mc_ratio);

    TH1D* h_data_mc_ratio = new TH1D(*h_data);
    h_data_mc_ratio->Divide(h_mc_total); 

    // Style Ratio Plot
    h_data_mc_ratio->SetTitle("");
    h_data_mc_ratio->GetXaxis()->SetTitle(h_mc_total->GetXaxis()->GetTitle());
    h_data_mc_ratio->GetYaxis()->SetTitle("Data/MC");
    h_data_mc_ratio->SetLineColor(kBlue);
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
    ratio_1.SetLineColor(kBlack);
    double line_min = h_data->GetBinLowEdge(1);
    double line_max = h_data->GetBinLowEdge(h_data->GetNbinsX()+1);
    ratio_1.DrawLine(line_min,1,line_max,1);
    //ratio_1.DrawLine(line_min,1,10,1);

    // Plot Output
    gStyle->SetOptStat(0); 
    c->Update();
    std::string out_name;
    out_name = plotDir + var_name + "_" + norm_label + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete pad1;
    delete pad2;
    delete c;
    delete plotter;
    delete tempData;
    delete tempMC;
    delete h_data;
    delete h_mc_total;
    delete h_data_mc_ratio;
}

void CCProtonPi0_Plotter::DrawDataMC_WithOtherData(MnvH1D* data, MnvH1D* mc, TGraph* otherData, std::string var_name, std::string ext_data_name, std::string plotDir)
{
    std::cout<<"Plotting for "<<var_name<<std::endl;

    double mc_ratio = 1.0;

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas("c","c",1280,1280);

    MnvH1D* tempData = new MnvH1D(*data);
    MnvH1D* tempMC = new MnvH1D(*mc);

    plotter->headroom = 1.75;
    plotter->data_line_width = 2;
    plotter->data_marker_size = 1.5;
    gStyle->SetEndErrorSize(6);
    plotter->DrawDataMCWithErrorBand(tempData, tempMC, mc_ratio, "N", false, NULL, NULL, false, true);

    // Draw With Blue Square Markers
    //otherData->SetMarkerColor(4);
    //otherData->SetMarkerSize(1.5);
    //otherData->SetMarkerStyle(21);

    //otherData->Draw("SAMEP");

    // Draw With Blue Line 
    otherData->SetLineColor(4);
    otherData->SetLineWidth(2);

    otherData->Draw("SAMEL");

    // ------------------------------------------------------------------------
    // Plot Labels 
    // ------------------------------------------------------------------------
    // Add Normalization Labels
    plotter->AddHistoTitle(data->GetTitle());
    AddNormBox(plotter, true, mc_ratio);

    // Plot Output
    gStyle->SetOptStat(0); 
    c->Update();
    std::string out_name;
    out_name = plotDir + var_name + "_" + ext_data_name + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete c;
    delete plotter;
    delete tempData;
    delete tempMC;
}

void CCProtonPi0_Plotter::DrawDataMC_FSIType(MnvH1D* data, MnvH1D* mc, std::vector<MnvH1D*> mc_FSIType, std::string var_name,  std::string plotDir)
{
    std::cout<<"Plotting for "<<var_name<<std::endl;

    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
    ApplyStyle(plotter);
    plotter->mc_line_width = 1;
    plotter->headroom = 1.75;
 
    std::string norm_label;
    double mc_ratio = GetMCNormalization(norm_label, false, data, mc);

    if (var_name.compare("Enu") == 0 || var_name.compare("muon_theta") == 0){
        plotter->headroom = 2.50;
    }

    TCanvas* c = new TCanvas("c");

    MnvH1D* h_data = new MnvH1D (*data);
    TH1D* h_data_stat_only = GetBinNormalizedTH1D(h_data, false);
    h_data_stat_only->SetMarkerStyle(20);
    h_data_stat_only->SetMarkerSize(1);
    h_data_stat_only->SetMarkerColor(kBlack);
    h_data_stat_only->SetLineWidth(1);
    h_data_stat_only->SetLineColor(kBlack);

    MnvH1D* h_NonInteracting = new MnvH1D (*mc_FSIType[0]);
    h_NonInteracting->SetTitle("#pi^{0} Non-interacting");
    h_NonInteracting->SetLineWidth(1);
    h_NonInteracting->SetLineColor(kRed+1);
    h_NonInteracting->SetFillColor(kRed+1);
    h_NonInteracting->SetFillStyle(1001);

    MnvH1D* h_Elastic = new MnvH1D (*mc_FSIType[1]);
    h_Elastic->SetTitle("#pi^{0} Elastic");
    h_Elastic->SetLineWidth(1);
    h_Elastic->SetLineColor(kOrange+6);
    h_Elastic->SetFillColor(kOrange+6);
    h_Elastic->SetFillStyle(1001);

    MnvH1D* h_Inelastic = new MnvH1D(*mc_FSIType[2]);
    h_Inelastic->SetTitle("#pi^{0} Inelastic");
    h_Inelastic->SetLineWidth(1);
    h_Inelastic->SetLineColor(kGray+1);
    h_Inelastic->SetFillColor(kGray+1);
    h_Inelastic->SetFillStyle(1001);

    MnvH1D* h_Cex = new MnvH1D(*mc_FSIType[3]);
    h_Cex->SetTitle("#pi^{#pm} #rightarrow #pi^{0}");
    h_Cex->SetLineWidth(1);
    h_Cex->SetLineColor(kGreen-2);
    h_Cex->SetFillColor(kGreen-2);
    h_Cex->SetFillStyle(1001);

    MnvH1D* h_MultiPion = new MnvH1D(*mc_FSIType[4]);
    h_MultiPion->SetTitle("Multi-#pi #rightarrow #pi^{0}");
    h_MultiPion->SetLineWidth(1);
    h_MultiPion->SetLineColor(kCyan-8);
    h_MultiPion->SetFillColor(kCyan-8);
    h_MultiPion->SetFillStyle(1001);

    MnvH1D* h_Other = new MnvH1D(*mc_FSIType[5]);
    h_Other->SetTitle("Other #pi^{0} production");
    h_Other->SetLineWidth(1);
    h_Other->SetLineColor(kMagenta+1);
    h_Other->SetFillColor(kMagenta+1);
    h_Other->SetFillStyle(1001);
 
    if (var_name.compare("Enu") == 0){
        h_data->GetXaxis()->SetRangeUser(0.0, 10.0);
        h_data_stat_only->GetXaxis()->SetRangeUser(0.0, 10.0);
        h_NonInteracting->GetXaxis()->SetRangeUser(0.0, 10.0);
        h_Elastic->GetXaxis()->SetRangeUser(0.0, 10.0);
        h_Inelastic->GetXaxis()->SetRangeUser(0.0, 10.0);
        h_Cex->GetXaxis()->SetRangeUser(0.0, 10.0);
        h_MultiPion->GetXaxis()->SetRangeUser(0.0, 10.0);
        h_Other->GetXaxis()->SetRangeUser(0.0, 10.0);
    }

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

    if (var_name.compare("Enu") == 0){
        leg_x_min = 0.17;
        leg_x_max = 0.60;
        leg_y_min = 0.42;
        leg_y_max = 0.85;
    }

    if (var_name.compare("muon_theta") == 0){
        leg_x_min = 0.62;
        leg_x_max = 0.90;
        leg_y_min = 0.42;
        leg_y_max = 0.90;
    }

    TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    ApplyStyle_Legend(legend);
    legend->SetTextSize(0.04);
    legend->AddEntry(h_data_stat_only, "Data #scale[0.85]{(3.33e20 POT)}", "lep" );
    legend->AddEntry(h_MultiPion, "Multi-#pi #rightarrow #pi^{0}","f");
    legend->AddEntry(h_Cex, "#pi^{#pm} #rightarrow #pi^{0}","f");
    legend->AddEntry(h_Other, "Zero-#pi #rightarrow #pi^{0}","f");
    legend->AddEntry(h_Inelastic, "#pi^{0} Inelastic","f");
    legend->AddEntry(h_Elastic, "#pi^{0} Elastic","f");
    legend->AddEntry(h_NonInteracting, "#pi^{0} Non-interacting","f");
    legend->Draw();

    // Add Normalization Labels
    if (isPaperComparison){
        TPaveText* text = new TPaveText(0.18, 0.80, 0.54, 0.90, "NDC");
        text->AddText("#scale[1.15]{c)  #nu_{#mu} + CH #rightarrow #mu^{-} + #pi^{0} + X}");
        text->AddText("#scale[0.85]{Area Normalized}");
        text->SetTextColor(kBlue);
        text->SetFillColor(kWhite);
        text->SetTextFont(42);
        text->Draw("SAME");
    }else{
        TLatex text;
        text.SetNDC();
        text.SetTextColor(kBlue);
        text.SetTextSize(0.03);
        text.SetTextAlign(22);
        text.DrawLatex(0.30,0.87,"Area Normalized");
        //text.DrawLatex(0.30,0.87,"POT Normalized"); // For pi0_P and pi0_theta from Trung's Paper
    }

    // Plot Output
    c->Update();
    std::string out_name;
    //out_name = plotDir + var_name + "_xsec_FSIType" + ".pdf"; 
    out_name = plotDir + var_name + "_xsec_FSIType_POT" + ".pdf"; 

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
}

void CCProtonPi0_Plotter::DrawDataMC_IntType(MnvH1D* data, MnvH1D* mc, std::vector<MnvH1D*> mc_IntType, std::string var_name,  std::string plotDir)
{
    std::cout<<"Plotting for "<<var_name<<std::endl;

    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
    ApplyStyle(plotter);
    plotter->mc_line_width = 1;
    plotter->headroom = 1.75;

    if (var_name.compare("Enu") == 0 || var_name.compare("muon_theta") == 0){
        plotter->headroom = 2.50;
    }

    TCanvas* c = new TCanvas("c");

    // I may use mc for area normalization later, but for now it is silenced
    (void) mc;

    MnvH1D* h_data = new MnvH1D (*data);
    TH1D* h_data_stat_only = GetBinNormalizedTH1D(h_data, false);
    h_data_stat_only->SetMarkerStyle(20);
    h_data_stat_only->SetMarkerSize(1);
    h_data_stat_only->SetMarkerColor(kBlack);
    h_data_stat_only->SetLineWidth(1);
    h_data_stat_only->SetLineColor(kBlack);

    MnvH1D* h_RES_Delta = new MnvH1D (*mc_IntType[0]);
    h_RES_Delta->SetLineWidth(1);
    h_RES_Delta->SetLineColor(17);
    h_RES_Delta->SetFillColor(17);
    h_RES_Delta->SetFillStyle(1001);

    MnvH1D* h_RES_Other = new MnvH1D (*mc_IntType[1]);
    h_RES_Other->SetLineWidth(1);
    h_RES_Other->SetLineColor(kGreen+2);
    h_RES_Other->SetFillColor(kGreen+2);
    h_RES_Other->SetFillStyle(1001);

    MnvH1D* h_NonRES = new MnvH1D (*mc_IntType[2]);
    h_NonRES->SetLineWidth(1);
    h_NonRES->SetLineColor(kRed+2);
    h_NonRES->SetFillColor(kRed+2);
    h_NonRES->SetFillStyle(1001);
 
    if (var_name.compare("Enu") == 0){
        h_data->GetXaxis()->SetRangeUser(0.0, 10.0);
        h_data_stat_only->GetXaxis()->SetRangeUser(0.0, 10.0);
        h_RES_Delta->GetXaxis()->SetRangeUser(0.0, 10.0);
        h_RES_Other->GetXaxis()->SetRangeUser(0.0, 10.0);
        h_NonRES->GetXaxis()->SetRangeUser(0.0, 10.0);
    }

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
    //TLegend *legend = new TLegend(0.55,0.57,0.90,0.85);  
    TLegend *legend = new TLegend(0.55,0.57,0.90,0.90);  
    ApplyStyle_Legend(legend);
    legend->SetTextSize(0.05);
    legend->AddEntry(h_data_stat_only, "Data #scale[0.85]{(3.33e20 POT)}", "lep" );
    legend->AddEntry(h_RES_Delta, "Delta resonance","f");
    legend->AddEntry(h_RES_Other, "Other resonances","f");
    legend->AddEntry(h_NonRES, "Non-Resonant","f");
    legend->Draw();

    // Add Normalization Labels
    if (isPaperComparison){
        TPaveText* text = new TPaveText(0.18, 0.80, 0.54, 0.90, "NDC");
        text->AddText("#scale[1.15]{c)  #nu_{#mu} + CH #rightarrow #mu^{-} + #pi^{0} + X}");
        text->AddText("#scale[0.85]{POT Normalized}");
        text->SetTextColor(kBlue);
        text->SetFillColor(kWhite);
        text->SetTextFont(42);
        text->Draw("SAME");
    }else{
        TLatex text;
        text.SetNDC();
        text.SetTextColor(kBlue);
        text.SetTextSize(0.03);
        text.SetTextAlign(22);
        text.DrawLatex(0.30,0.87,"POT Normalized");
    }

    // Plot Output
    c->Update();
    std::string out_name;
    out_name = plotDir + var_name + "_xsec_IntType" + ".pdf"; 

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

void CCProtonPi0_Plotter::DrawDataMC_PaperStyle(MnvH1D* data, MnvH1D* mc, MnvH1D* mc_BeforeFSI, std::string var_name,  std::string plotDir)
{
    std::cout<<"Plotting for "<<var_name<<std::endl;

    double mc_ratio = 1.0;

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
    ApplyStyle(plotter);
    plotter->headroom = 1.75;

    TCanvas* c = new TCanvas("c");

    MnvH1D* tempData = new MnvH1D(*data);
    MnvH1D* tempMC = new MnvH1D(*mc);

    if (var_name.compare("Enu") == 0 || var_name.compare("muon_theta") == 0){
        plotter->headroom = 2.50;
    }

    if (var_name.compare("Enu") == 0){
        tempData->GetXaxis()->SetRangeUser(0.0, 10.0);
        tempMC->GetXaxis()->SetRangeUser(0.0, 10.0);
    }

    plotter->DrawDataMCWithErrorBand(tempData, tempMC, mc_ratio, "N", false, NULL, NULL, false, true, true);

    TH1D* h_mc_BeforeFSI = GetBinNormalizedTH1D(mc_BeforeFSI);
    h_mc_BeforeFSI->SetLineColor(kRed);
    h_mc_BeforeFSI->SetLineWidth(2);
    h_mc_BeforeFSI->SetLineStyle(2);
    h_mc_BeforeFSI->SetFillStyle(0);

    h_mc_BeforeFSI->Smooth(2);
    h_mc_BeforeFSI->Draw("SAME HISTC");

    // ------------------------------------------------------------------------
    // Plot Labels 
    // ------------------------------------------------------------------------

    // Add Legend
    tempData->SetMarkerStyle(plotter->data_marker);
    tempData->SetMarkerSize(plotter->data_marker_size);
    tempData->SetMarkerColor(plotter->data_color);
    tempData->SetLineColor(plotter->data_color);
    tempData->SetLineWidth(plotter->data_line_width);

    tempMC->SetLineColor(plotter->mc_color);
    tempMC->SetLineWidth(plotter->mc_line_width);

    TLegend *legend = new TLegend(0.53,0.62,0.80,0.85);  
    ApplyStyle_Legend(legend);    
    legend->SetTextSize(0.06);
    legend->AddEntry(tempData, "Data #scale[0.85]{(3.33e20 POT)}", "lep" );
    legend->AddEntry(tempMC, "GENIE w/ FSI", "l");
    legend->AddEntry(h_mc_BeforeFSI, "GENIE w/o FSI","l");
    legend->Draw();

    // Add Normalization Labels
    if (isPaperComparison){
        TPaveText* text = new TPaveText(0.18, 0.80, 0.54, 0.90, "NDC");
        text->AddText("#scale[1.15]{c)  #nu_{#mu} + CH #rightarrow #mu^{-} + #pi^{0} + X}");
        text->AddText("#scale[0.85]{POT Normalized}");
        text->SetTextColor(kBlue);
        text->SetFillColor(kWhite);
        text->SetTextFont(42);
        text->Draw("SAME");
    }else{
        TLatex text;
        text.SetNDC();
        text.SetTextColor(kBlue);
        text.SetTextSize(0.03);
        text.SetTextAlign(22);
        text.DrawLatex(0.30,0.87,"POT Normalized");
    }

    // Plot Output
    gStyle->SetOptStat(0); 
    c->Update();
    std::string out_name;
    out_name = plotDir + var_name + "_xsec_BeforeFSI" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete legend;
    delete c;
    delete plotter;
    delete tempData;
    delete tempMC;
    delete h_mc_BeforeFSI;
}


void CCProtonPi0_Plotter::DrawDataMC_BeforeFSI(MnvH1D* data, MnvH1D* mc, MnvH1D* mc_BeforeFSI, std::string var_name,  std::string plotDir)
{
    std::cout<<"Plotting for "<<var_name<<std::endl;

    double mc_ratio = 1.0;

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas("c","c",800,800);

    MnvH1D* tempData = new MnvH1D(*data);
    MnvH1D* tempMC = new MnvH1D(*mc);

    plotter->headroom = 1.75;
    plotter->axis_title_offset_y = 1.1;
    plotter->data_line_width = 2;
    plotter->data_marker_size = 1.5;
    gStyle->SetEndErrorSize(6);
    plotter->DrawDataMCWithErrorBand(tempData, tempMC, mc_ratio, "N", false, NULL, NULL, false, true, true);

    TH1D* h_mc_BeforeFSI = GetBinNormalizedTH1D(mc_BeforeFSI);
    h_mc_BeforeFSI->SetLineColor(kBlue);
    h_mc_BeforeFSI->SetLineWidth(3);
    h_mc_BeforeFSI->SetFillStyle(0);

    h_mc_BeforeFSI->Smooth(2);
    h_mc_BeforeFSI->Draw("SAME HISTC");

    // ------------------------------------------------------------------------
    // Plot Labels 
    // ------------------------------------------------------------------------

    // Add Legend
    tempData->SetMarkerStyle(plotter->data_marker);
    tempData->SetMarkerSize(plotter->data_marker_size);
    tempData->SetMarkerColor(plotter->data_color);
    tempData->SetLineColor(plotter->data_color);
    tempData->SetLineWidth(plotter->data_line_width);

    tempMC->SetLineColor(plotter->mc_color);
    tempMC->SetLineWidth(plotter->mc_line_width);

    TLegend *legend = new TLegend(0.40,0.7,0.90,0.9);  
    ApplyStyle_Legend(legend);    
    legend->AddEntry(tempData, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(tempMC, "GENIE w/ FSI", "l");
    legend->AddEntry(h_mc_BeforeFSI, "GENIE w/o FSI","l");
    legend->Draw();

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");

    // Plot Output
    gStyle->SetOptStat(0); 
    c->Update();
    std::string out_name;
    out_name = plotDir + var_name + "_xsec_BeforeFSI" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete legend;
    delete c;
    delete plotter;
    delete tempData;
    delete tempMC;
    delete h_mc_BeforeFSI;
}

void CCProtonPi0_Plotter::DrawDataMCRatio(MnvH1D* data, MnvH1D* mc, std::string var_name, std::string plotDir, bool isPOTNorm, bool isXSec)
{
    // ------------------------------------------------------------------------
    // MC Normalization 
    // ------------------------------------------------------------------------
    std::string norm_label; 
    double mc_ratio = GetMCNormalization(norm_label, isPOTNorm, data, mc);

    // Cross Section Calculation Already Includes POT Normalization
    if (isXSec && isPOTNorm) mc_ratio = 1.0;

    // ------------------------------------------------------------------------
    // Plot
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas("c","c",1280,800);

    plotter->DrawDataMCRatio(data, mc, mc_ratio);

    // Add Plot Labels
    plotter->AddHistoTitle(data[0].GetTitle());
    AddNormBox(plotter, isPOTNorm, mc_ratio);

    // Print Plot
    c->Print(Form("%s%s%s%s%s",plotDir.c_str(),var_name.c_str(),"_ratio_",norm_label.c_str(),".pdf"), "pdf");

    delete c;
    delete plotter;
}

void CCProtonPi0_Plotter::DrawMCWithErrorBand(rootDir& dir, std::string var_name, std::string plotDir)
{
    std::string rootDir_mc = dir.mc;

    TFile* f_mc = new TFile(rootDir_mc.c_str());

    TH1D* mc = (TH1D*)f_mc->Get(var_name.c_str());

    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas("c","c",1280,800);

    // Plot
    plotter->DrawMCWithErrorBand(mc);

    // Print Plot
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".pdf"), "pdf");

    delete c;
    delete plotter;
}

void CCProtonPi0_Plotter::DrawMnvH1D(MnvH1D* hist1D, std::string var_name, std::string plotDir)
{
    TH1D* h = GetBinNormalizedTH1D(hist1D);

    Draw1DHist(h, var_name, plotDir);

}

void CCProtonPi0_Plotter::DrawMnvH2D(MnvH2D* hist2D, std::string var_name, std::string plotDir, bool isMC)
{
    std::string plot_label = isMC ? "MC" : "Data";

    // Canvas
    Double_t w = 800; 
    Double_t h = 800;
    TCanvas* c = new TCanvas("c","c",w,h);
    c->SetWindowSize(w,h);

    if (thesisStyle){
        //hist2D->GetXaxis()->SetTitle("P_{#mu} [GeV]");
        hist2D->GetXaxis()->SetTitleFont(62);
        hist2D->GetXaxis()->SetTitleSize(0.06);
        hist2D->GetXaxis()->CenterTitle();
        hist2D->GetXaxis()->SetTitleOffset(1.15);
        hist2D->GetXaxis()->SetLabelFont(42);
        hist2D->GetXaxis()->SetLabelSize(0.05);
        hist2D->GetXaxis()->SetNdivisions(505);

        hist2D->GetYaxis()->SetTitleFont(62);
        hist2D->GetYaxis()->SetTitleSize(0.06);
        //hist2D->GetYaxis()->CenterTitle();
        hist2D->GetYaxis()->SetTitleOffset(1.2);
        hist2D->GetYaxis()->SetLabelFont(42);
        hist2D->GetYaxis()->SetLabelSize(0.05);
        hist2D->GetYaxis()->SetNdivisions(505);
        TGaxis::SetMaxDigits(3);
    }

    // Pad
    TPad *p = new TPad("p","p",0.05,0.05,0.95,0.95);
    p->Draw();

    p->cd();
    //hist2D->GetYaxis()->SetTitleOffset(1.8);
    hist2D->Draw("colz");
    gPad->Update();

    gStyle->SetOptStat(0); 
    c->Update();

    c->Print(Form("%s%s_%s%s",plotDir.c_str(),var_name.c_str(),plot_label.c_str(), ".pdf"), "pdf");

    delete p;
    delete c;
}

void CCProtonPi0_Plotter::DrawMnvH2D_Signal(rootDir root_dir, std::string var_name, std::string plotDir, double nBckg, bool isMC)
{
    TFile* f;
    MnvH2D* hist2D;

    if (isMC){
        f = new TFile(root_dir.mc.c_str());
        var_name = var_name + "_mc_reco_signal";
        hist2D = GetMnvH2D(f, var_name);
    }else{
        f = new TFile(root_dir.data.c_str());
        hist2D = GetBckgSubtractedData_2D(root_dir, var_name, nBckg);
    }
    DrawMnvH2D(hist2D, var_name, plotDir, isMC);

    delete f;
    delete hist2D;
}

void CCProtonPi0_Plotter::DrawMnvH2D(std::string root_dir, std::string var_name, std::string plotDir)
{
    std::size_t found = root_dir.find("Data");
    
    bool isMC = (found == std::string::npos) ? true : false;

    TFile* f = new TFile(root_dir.c_str());
    MnvH2D* hist2D = GetMnvH2D(f, var_name);

    DrawMnvH2D(hist2D, var_name, plotDir, isMC);

    delete f;
    delete hist2D;
}

void CCProtonPi0_Plotter::DrawMnvH1D(rootDir& dir, std::string var_name, std::string plotDir)
{
    // Get Histogram
    std::string root_dir;
    std::size_t found = var_name.find("data");
    if (found != std::string::npos) root_dir = dir.data;
    else root_dir = dir.mc;

    TFile* f = new TFile(root_dir.c_str());
    MnvH1D* hist1D = (MnvH1D*)f->Get(var_name.c_str());

    DrawMnvH1D(hist1D, var_name, plotDir);
    delete f;
}

void CCProtonPi0_Plotter::Draw1DHist(rootDir& dir, std::string var_name, std::string plotDir, bool isLogScale)
{
    // Get Histogram
    std::string root_dir;
    std::size_t found = var_name.find("data");
    if (found != std::string::npos) root_dir = dir.data;
    else root_dir = dir.mc;

    TFile* f = new TFile(root_dir.c_str());
    TH1D* hist1D = new TH1D( * dynamic_cast<TH1D*>(f->Get(var_name.c_str())) );

    Draw1DHist(hist1D, var_name, plotDir, isLogScale);

    delete f;
}


void CCProtonPi0_Plotter::Draw1DHist(TH1* hist1D, std::string var_name, std::string plotDir, bool isLogScale)
{
    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.03);

    std::cout<<"Plotting "<<var_name<<std::endl;

    // Create Canvas
    TCanvas* c = new TCanvas("c","c",800,800);
    if(isLogScale) c->SetLogy();

    //hist1D->Scale(0.1,"width");

    // Plot Options
    hist1D->SetLineColor(kRed);
    hist1D->SetLineWidth(3);
    hist1D->SetFillColor(kRed);
    hist1D->SetFillStyle(3010);

    if (thesisStyle){
        //hist1D->GetXaxis()->SetTitle("P_{#mu} [GeV]");
        hist1D->GetXaxis()->SetTitleFont(62);
        hist1D->GetXaxis()->SetTitleSize(0.06);
        hist1D->GetXaxis()->CenterTitle();
        hist1D->GetXaxis()->SetTitleOffset(1.15);
        hist1D->GetXaxis()->SetLabelFont(42);
        hist1D->GetXaxis()->SetLabelSize(0.05);
        hist1D->GetXaxis()->SetNdivisions(408);

        hist1D->GetYaxis()->SetTitleFont(62);
        hist1D->GetYaxis()->SetTitleSize(0.06);
        //hist1D->GetYaxis()->CenterTitle();
        hist1D->GetYaxis()->SetTitleOffset(1.2);
        hist1D->GetYaxis()->SetLabelFont(42);
        hist1D->GetYaxis()->SetLabelSize(0.05);
        TGaxis::SetMaxDigits(3);
    }
    
    hist1D->SetMaximum(hist1D->GetMaximum()*1.50);
    hist1D->SetMinimum(0.0);
    hist1D->Draw("HIST");

    // Find Peak
    int max_bin = hist1D->GetMaximumBin();
    double max_value = hist1D->GetBinContent(max_bin);
    double max_bin_location = hist1D->GetBinCenter(max_bin);
    std::cout<<"\tPeak = "<<max_bin_location<<std::endl;

    // Calculate FWHM
    int bin1 = hist1D->FindFirstBinAbove(max_value/2);
    int bin2 = hist1D->FindLastBinAbove(max_value/2);
    double fwhm = (hist1D->GetBinLowEdge(bin2) + hist1D->GetBinWidth(bin2)) - hist1D->GetBinLowEdge(bin1);
    std::cout<<"\tFWHM = "<<fwhm<<std::endl;

    if (!thesisStyle){ 
        text.DrawLatex(0.20,0.85,Form("%s%3.2f", "Peak at ",max_bin_location));
        text.DrawLatex(0.20,0.8,Form("%s%3.2f", "FWHM = ",fwhm));
    }

    // Add Pi0 InvMass Line
    std::size_t found = var_name.find("invMass");
    if (found != std::string::npos){
        TLine pi0Mass;
        pi0Mass.SetLineWidth(2);
        pi0Mass.SetLineColor(kBlue);
        pi0Mass.DrawLine(134.98,0,134.98,max_value);
    }

    gStyle->SetOptStat(0); 

    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".pdf"), "pdf");

    delete c;
    delete hist1D;
}

void CCProtonPi0_Plotter::Draw1DHist_Threshold(rootDir& dir, std::string var_name, std::string plotDir, double threshold, bool isLogScale)
{
    std::string root_dir = dir.mc;

    // Get Histogram
    TFile* f = new TFile(root_dir.c_str());
    TH1D* hist1D = (TH1D*)f->Get(var_name.c_str());

    // Reset Bins below the threshold
    int nBinsX = hist1D->GetNbinsX();
    for (int xBin = 1; xBin <= nBinsX; xBin++ ){
        int nEvents = hist1D->GetBinContent(xBin);
        if (nEvents <= threshold){
            hist1D->SetBinContent(xBin,0);
        }
    }

    // Create Canvas
    TCanvas* c = new TCanvas("c","c",1280,800);
    if(isLogScale) c->SetLogy();

    // Plot Options
    hist1D->SetLineColor(kRed);
    hist1D->SetLineWidth(3);
    hist1D->SetFillColor(kRed);
    hist1D->SetFillStyle(3010);

    hist1D->Draw();
    gPad->Update();
    gStyle->SetOptStat("nemr"); 

    int max_bin = hist1D->GetMaximumBin();
    double max_bin_value = hist1D->GetBinCenter(max_bin);
    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.03);
    text.DrawLatex(0.78,0.7,Form("%s%3.2f", "Peak at ",max_bin_value));

    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".pdf"), "pdf");

    delete c;
    delete f; 
}

void CCProtonPi0_Plotter::Draw2DHist(rootDir& dir, std::string var_name, std::string plotDir, double threshold)
{
    std::string root_dir = dir.mc;

    // Get Histogram
    TFile* f = new TFile(root_dir.c_str());
    TH2D* hist2D = (TH2D*)f->Get(var_name.c_str());

    // Canvas
    Double_t w = 800; 
    Double_t h = 800;
    TCanvas* c = new TCanvas("c","c",w,h);

    // Reset Bins below the threshold
    int nBinsX = hist2D->GetNbinsX();
    int nBinsY = hist2D->GetNbinsY();
    for (int xBin = 1; xBin <= nBinsX; xBin++ ){
        for (int yBin = 1; yBin <=nBinsY; yBin++){
            int nEvents = hist2D->GetBinContent(xBin,yBin);
            if (nEvents <= threshold){
                hist2D->SetBinContent(xBin,yBin,0);
            }
        }
    }

    // Pad
    TPad *p = new TPad("p","p",0.05,0.05,0.95,0.95);
    p->Draw();
    p->cd();

    if (thesisStyle){
        hist2D->GetXaxis()->SetTitleFont(62);
        hist2D->GetXaxis()->SetTitleSize(0.06);
        hist2D->GetXaxis()->CenterTitle();
        hist2D->GetXaxis()->SetTitleOffset(1.15);
        hist2D->GetXaxis()->SetLabelFont(42);
        hist2D->GetXaxis()->SetLabelSize(0.04);
        hist2D->GetXaxis()->SetNdivisions(408);

        //hist2D->GetYaxis()->SetTitle("cos(#theta_{#gamma#gamma})");
        hist2D->GetYaxis()->SetTitleFont(62);
        hist2D->GetYaxis()->SetTitleSize(0.06);
        //hist2D->GetYaxis()->CenterTitle();
        hist2D->GetYaxis()->SetTitleOffset(1.1);
        hist2D->GetYaxis()->SetLabelFont(42);
        hist2D->GetYaxis()->SetLabelSize(0.04);
    }

    gStyle->SetPaintTextFormat("2.0f"); 
    hist2D->GetYaxis()->SetTitleOffset(1.8);
    //hist2D->Draw("colz text");
    hist2D->Draw("colz");

    gPad->Update();

    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".pdf"), "pdf");

    delete p;
    delete c;
    delete f;
}

void CCProtonPi0_Plotter::Draw3DHist(rootDir& dir, std::string var_name, std::string plotDir)
{
    std::string root_dir = dir.mc;

    // Get Histogram
    TFile* f = new TFile(root_dir.c_str());
    TH3D* hist3D = (TH3D*)f->Get(var_name.c_str());

    // Canvas
    Double_t w = 800; 
    Double_t h = 800;
    TCanvas* c = new TCanvas("c","c",w,h);
    c->SetWindowSize(w,h);

    //gStyle->SetCanvasPreferGL(true);
    // Pad
    TPad *p = new TPad("p","p",0.05,0.05,0.95,0.95);
    p->Draw();

    p->cd();
    hist3D->GetXaxis()->SetTitleOffset(1.8);
    hist3D->GetYaxis()->SetTitleOffset(2.4);
    hist3D->GetZaxis()->SetTitleOffset(1.8);
    hist3D->Draw("BOX");
    gPad->Update();

    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".pdf"), "pdf");

    delete p;
    delete c;
    delete f;
}

void CCProtonPi0_Plotter::AddCutArrow(MnvPlotter* plotter, CutArrow &cutArrow, double hist_max, double arrow_length)
{
    double cut_location = cutArrow.cut_location;
    double ymin = 0;
    double ymax = hist_max;
    std::string arrow_direction = cutArrow.arrow_direction;
    plotter->arrow_line_width = 2;
    plotter->arrow_line_style = 1;
    plotter->arrow_type = "";

    plotter->AddCutArrow(cut_location, ymin, ymax, arrow_length, arrow_direction); 
}


void CCProtonPi0_Plotter::DrawStackedMC_BckgAll(rootDir &dir, std::string var_name, std::string plotDir, int nCutArrows, CutArrow cutArrow1, CutArrow cutArrow2)
{
    (void) nCutArrows;
    (void) cutArrow1;
    (void) cutArrow2;

    std::string rootDir_mc = dir.mc;

    TFile* f_mc = new TFile(rootDir_mc.c_str());

    std::string var;

    // ------------------------------------------------------------------------
    // Fill TObjArray - For MC Histograms
    // ------------------------------------------------------------------------
    TObjArray* mc_hists = new TObjArray;
    MnvH1D* temp;

    // Get All Background
    var = Form("%s_%d",var_name.c_str(),2);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Background");
    temp->SetLineColor(kRed);
    temp->SetFillColor(kRed);
    mc_hists->Add(temp);

    // Get Signal
    var = Form("%s_%d",var_name.c_str(),1);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Signal");
    temp->SetLineColor(kGreen);
    temp->SetFillColor(kGreen);
    mc_hists->Add(temp);

    // ------------------------------------------------------------------------
    // Plot  - If you want A Log Plot, axis_minimum = 0.1
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas(var_name.c_str(),var_name.c_str(),1280,800);
    
    plotter->headroom = 1.75;
    plotter->data_line_width = 2;
    plotter->data_marker_size = 1.5;
    gStyle->SetEndErrorSize(6);

    plotter->DrawStackedMC(mc_hists,1,"TR");

    // Print Plot
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),"_mc_bckg_all.pdf"), "pdf");

    delete c;
    delete plotter;
}

void CCProtonPi0_Plotter::DrawStackedMC_BckgType(rootDir &dir, std::string var_name, std::string plotDir, int nCutArrows, CutArrow cutArrow1, CutArrow cutArrow2)
{
    std::string rootDir_mc = dir.mc;

    TFile* f_mc = new TFile(rootDir_mc.c_str());

    // ------------------------------------------------------------------------
    // Fill TObjArray - For MC Histograms
    // ------------------------------------------------------------------------
    double hist_max = 0;
    double bin_width;

    TObjArray* mc_hists = new TObjArray;
    FormTObjArray_BckgType(f_mc, var_name, mc_hists, hist_max, bin_width);

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas("c","c",1280,800);
    ApplyStyle(plotter);
    //plotter->axis_minimum = 0.1;
    plotter->axis_minimum = 0.0;
    plotter->DrawStackedMC(mc_hists,1,"TR",0);

    // Add Pi0 InvMass Line
    std::size_t found = var_name.find("invMass");
    if (found != std::string::npos){
        TLine pi0Mass;
        pi0Mass.SetLineWidth(2);
        pi0Mass.SetLineColor(kBlue);
        pi0Mass.DrawLine(134.98,0,134.98,hist_max);
    }

    // Add Plot Labels
    std::string var = Form("%s_%d",var_name.c_str(),0);
    MnvH1D* temp = (MnvH1D*)f_mc->Get(var.c_str());
    plotter->AddHistoTitle(temp->GetTitle());

    // If Cut Histogram - Add Cut Arrows
    if (nCutArrows == 1){
        AddCutArrow(plotter, cutArrow1, hist_max, bin_width);
    }else if (nCutArrows == 2){
        AddCutArrow(plotter, cutArrow1, hist_max, bin_width);
        AddCutArrow(plotter, cutArrow2, hist_max, bin_width);
    }

    // Print Plot
    c->Print(Form("%s%s%s%s",plotDir.c_str(),var_name.c_str(),"_mc_bckg_Type",".pdf"), "pdf");

    delete c;
    delete plotter;
}

void CCProtonPi0_Plotter::DrawStackedMC_GammaByPDG(rootDir &dir, std::string var_name, int gammaID, std::string plotDir)
{
    // ------------------------------------------------------------------------- 
    // Get Histograms from File 
    // ------------------------------------------------------------------------- 
    std::string rootDir_mc = dir.mc;

    TFile* f_mc = new TFile(rootDir_mc.c_str());

    std::string pi0 = Form("g%d_%s_pi0",gammaID,var_name.c_str());
    std::string pi = Form("g%d_%s_pi",gammaID,var_name.c_str());
    std::string proton = Form("g%d_%s_proton",gammaID,var_name.c_str());
    std::string neutron = Form("g%d_%s_neutron",gammaID,var_name.c_str());
    std::string muon = Form("g%d_%s_muon",gammaID,var_name.c_str());

    TH1D* h_pi0 = new TH1D;   
    TH1D* h_pi = new TH1D;   
    TH1D* h_proton = new TH1D;   
    TH1D* h_neutron = new TH1D;   
    TH1D* h_muon = new TH1D;   

    h_pi0 = (TH1D*)f_mc->Get(pi0.c_str());
    h_pi = (TH1D*)f_mc->Get(pi.c_str());
    h_proton = (TH1D*)f_mc->Get(proton.c_str());
    h_neutron = (TH1D*)f_mc->Get(neutron.c_str());
    h_muon = (TH1D*)f_mc->Get(muon.c_str());

    // Set Histogram Style
    h_pi0->SetFillColor(kGreen);
    h_pi0->SetLineColor(kGreen);
    h_pi0->SetMarkerStyle(21);
    h_pi0->SetMarkerColor(kGreen);

    h_pi->SetFillColor(kRed);
    h_pi->SetLineColor(kRed);
    h_pi->SetMarkerStyle(21);
    h_pi->SetMarkerColor(kRed);

    h_proton->SetFillColor(kBlue);
    h_proton->SetLineColor(kBlue);
    h_proton->SetMarkerStyle(21);
    h_proton->SetMarkerColor(kBlue);

    h_neutron->SetFillColor(kCyan);
    h_neutron->SetLineColor(kCyan);
    h_neutron->SetMarkerStyle(21);
    h_neutron->SetMarkerColor(kCyan);

    h_muon->SetFillColor(kMagenta);
    h_muon->SetLineColor(kMagenta);
    h_muon->SetMarkerStyle(21);
    h_muon->SetMarkerColor(kMagenta);

    // Create Canvas and Form THStack
    TCanvas* c = new TCanvas("c","c",1280,800);
    THStack* hs = new THStack("hs",Form("%s by Particle for Gamma %d",var_name.c_str(),gammaID));
    TLegend* legend = new TLegend(0.7,0.8,0.9,0.9);
    ApplyStyle_Legend(legend);    

    //c->SetLogy();
    // Legend
    legend->AddEntry(h_pi0, "#pi^{0}", "f");
    legend->AddEntry(h_pi, "#pi^{#pm}", "f");
    legend->AddEntry(h_proton, "proton", "f");
    legend->AddEntry(h_neutron, "neutron", "f");
    legend->AddEntry(h_muon, "#mu^{-}", "f");

    hs->Add(h_pi0);
    hs->Add(h_pi);
    hs->Add(h_proton);
    hs->Add(h_neutron);
    hs->Add(h_muon);
    hs->Draw();
    hs->GetXaxis()->SetTitle(Form("%s by Particle for Gamma %d",var_name.c_str(),gammaID));
    hs->GetYaxis()->SetTitle("N(Events)");

    legend->Draw();

    c->Print(Form("%s%s%d%s%s%s",plotDir.c_str(),"gamma_",gammaID,"_mc_",var_name.c_str(),"_by_particle.pdf"), "pdf");

    delete legend;
    delete hs;
    delete c;
}

void CCProtonPi0_Plotter::AddNormBox(MnvPlotter* plotter, bool isPOTNorm, double mc_ratio)
{
    const double y_pos = 0.88;
    const double y_diff = 0.033;
    if (isPOTNorm){
        plotter->AddPOTNormBox(data_POT,mc_POT,0.3,y_pos-y_diff);
    }else{
        plotter->AddAreaNormBox(1.0,mc_ratio,0.3,y_pos-y_diff);
    }
}

double CCProtonPi0_Plotter::GetMCNormalization(std::string &norm_label, bool isPOTNorm, MnvH1D* data, MnvH1D* mc)
{
    double mc_ratio;
    if (isPOTNorm){
        mc_ratio = POT_ratio;
        norm_label = "POT";
    }else{
        double area_data = data->Integral("width");
        double area_mc = mc->Integral("width");
        mc_ratio = area_data/area_mc;
        norm_label = "Area";
    }

    return mc_ratio;
}

void CCProtonPi0_Plotter::SaveRecoRatioPoints(rootDir& dir, std::string var_name, std::string plotDir)     
{   
    // Open Output Text file
    ofstream text;  
    std::string textFile = plotDir + var_name + ".txt";
    text.open(textFile.c_str());    

    // Get Histogram    
    std::string root_dir = dir.mc;
    TFile* f = new TFile(root_dir.c_str());     
    TH2D* hist2D = (TH2D*)f->Get(var_name.c_str());   

    // Get Evis Points
    int nBinsX = hist2D->GetNbinsX();   
    int nBinsY = hist2D->GetNbinsY();   

    for (int xBin = 1; xBin <= nBinsX; xBin++ ){    
        double sum = 0;     
        double n = 0;   
        for (int yBin = 1; yBin <=nBinsY; yBin++){  
            double nEvents = hist2D->GetBinContent(xBin,yBin);  
            double y = hist2D->GetYaxis()->GetBinCenter(yBin);  
            if ( nEvents > 1){  
                sum = sum + y*nEvents;  
                n = n + nEvents;    
            }   
        }   
        text<<std::endl;    

        if (n > 0){     
            double avg = sum / n;   
            double x = hist2D->GetXaxis()->GetBinCenter(xBin);  
            text<<x<<" "<<avg<<" "<<n<<std::endl;   
        }   
    }   

    text.close();   
}

void CCProtonPi0_Plotter::Save2DHistPoints(rootDir& dir, std::string var_name, std::string plotDir)     
{   
    // Open Output Text file
    ofstream text;  
    std::string textFile = plotDir + var_name + ".txt";
    text.open(textFile.c_str());    

    // Get Histogram    
    std::string root_dir = dir.mc;
    TFile* f = new TFile(root_dir.c_str());     
    TH2D* hist2D = (TH2D*)f->Get(var_name.c_str());   

    // Get Evis Points
    int nBinsX = hist2D->GetNbinsX();   
    int nBinsY = hist2D->GetNbinsY();   

    for (int xBin = 1; xBin <= nBinsX; xBin++ ){    
        for (int yBin = 1; yBin <=nBinsY; yBin++){  
            double x = hist2D->GetXaxis()->GetBinLowEdge(xBin);  
            double y = hist2D->GetYaxis()->GetBinLowEdge(yBin);  
            double nEvents = hist2D->GetBinContent(xBin,yBin);  
            if ( nEvents > 5){  
                text<<x<<" "<<y<<" "<<nEvents<<std::endl;
            }   
        }   
    }   

    text.close();   
}

void CCProtonPi0_Plotter::DrawEfficiencyCurve(rootDir& dir, std::string var_name, std::string plotDir)
{
    // Get Histogram
    std::string root_dir;
    std::size_t found = var_name.find("data");
    if (found != std::string::npos) root_dir = dir.data;
    else root_dir = dir.mc;

    TFile* f = new TFile(root_dir.c_str());
    MnvH1D* hist1D = (MnvH1D*)f->Get(var_name.c_str());

    // Create Canvas
    TCanvas* c = new TCanvas("c","c",800,800);

    // Plot Options
    hist1D->GetYaxis()->SetTitle("Reconstruction Efficiency");
    hist1D->SetMinimum(0.0);
    hist1D->SetMaximum(0.15);
    hist1D->SetLineColor(kRed);
    hist1D->SetLineWidth(3);
    hist1D->SetFillColor(kWhite);

    if (thesisStyle){
        //hist1D->GetXaxis()->SetTitle("P_{#mu}");
        hist1D->GetXaxis()->SetTitleFont(62);
        hist1D->GetXaxis()->SetTitleSize(0.06);
        hist1D->GetXaxis()->CenterTitle();
        hist1D->GetXaxis()->SetTitleOffset(1.15);
        hist1D->GetXaxis()->SetLabelFont(42);
        hist1D->GetXaxis()->SetLabelSize(0.05);
        hist1D->GetXaxis()->SetNdivisions(408);

        hist1D->GetYaxis()->SetTitleFont(62);
        hist1D->GetYaxis()->SetTitleSize(0.06);
        //hist1D->GetYaxis()->CenterTitle();
        hist1D->GetYaxis()->SetTitleOffset(1.2);
        hist1D->GetYaxis()->SetLabelFont(42);
        hist1D->GetYaxis()->SetLabelSize(0.05);
        TGaxis::SetMaxDigits(3);
    }

    hist1D->Draw("HIST");
    gPad->Update();
    gStyle->SetOptStat(0); 

    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".pdf"), "pdf");
    delete c;
    delete hist1D;
    delete f;

}


double CCProtonPi0_Plotter::FormTObjArray_BckgType(TFile* f_mc, std::string var_name, TObjArray* mc_hists, double &hist_max, double &bin_width) 
{
    std::string var;
    double area_mc = 0;
    double max_bin;
    MnvH1D* temp;

    // Get Signal
    var = Form("%s_%d",var_name.c_str(),1);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Signal");
    bin_width = temp->GetBinWidth(1);
    max_bin = temp->GetMaximumBin();
    hist_max = hist_max + temp->GetBinContent(max_bin);
    temp->SetLineColor(kGreen);
    temp->SetFillColor(kGreen);
    mc_hists->Add(temp);
    area_mc += temp->Integral();

    // Get Bckg: WithPi0
    var = Form("%s_%d",var_name.c_str(),3);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: #pi^{0} + X");
    max_bin = temp->GetMaximumBin();
    hist_max = hist_max + temp->GetBinContent(max_bin);
    temp->SetLineColor(kRed);
    temp->SetFillColor(kRed);
    mc_hists->Add(temp);
    area_mc += temp->Integral();

    // Get Bckg: QELike
    var = Form("%s_%d",var_name.c_str(),4);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: QELike");
    max_bin = temp->GetMaximumBin();
    hist_max = hist_max + temp->GetBinContent(max_bin);
    temp->SetLineColor(kOrange);
    temp->SetFillColor(kOrange);
    mc_hists->Add(temp);
    area_mc += temp->Integral();

//    // Get Bckg: 2p2h 
//    var = Form("%s_%d",var_name.c_str(),10);
//    temp = (MnvH1D*)f_mc->Get(var.c_str());
//    temp->SetTitle("Bckg: 2p2h");
//    max_bin = temp->GetMaximumBin();
//    hist_max = hist_max + temp->GetBinContent(max_bin);
//    temp->SetLineColor(kBlack);
//    temp->SetFillColor(kBlack);
//    mc_hists->Add(temp);
//    area_mc += temp->Integral();

    // Get Bckg: SinglePiPlus
    var = Form("%s_%d",var_name.c_str(),5);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: #pi^{#pm}");
    max_bin = temp->GetMaximumBin();
    hist_max = hist_max + temp->GetBinContent(max_bin);
    temp->SetLineColor(kBlue);
    temp->SetFillColor(kBlue);
    mc_hists->Add(temp);
    area_mc += temp->Integral();

    // Get Bckg: Other
    var = Form("%s_%d",var_name.c_str(),6);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: Other");
    max_bin = temp->GetMaximumBin();
    hist_max = hist_max + temp->GetBinContent(max_bin);
    temp->SetLineColor(kGray);
    temp->SetFillColor(kGray);
    mc_hists->Add(temp);
    area_mc += temp->Integral();


    return area_mc;
}

void CCProtonPi0_Plotter::DrawSignalMC(rootDir &dir, std::string var_name, std::string plotDir, int nCutArrows, CutArrow cutArrow1, CutArrow cutArrow2)
{
    std::string rootDir_mc = dir.mc;

    TFile* f_mc = new TFile(rootDir_mc.c_str());

    std::string var;
    double hist_max = 0;
    double max_bin;
    double bin_width;

    // ------------------------------------------------------------------------
    // Fill TObjArray - For MC Histograms
    // ------------------------------------------------------------------------
    TObjArray* mc_hists = new TObjArray;
    MnvH1D* signal;

    // ------------------------------------------------------------------------
    // Get Signal
    // ------------------------------------------------------------------------
    var = Form("%s_%d",var_name.c_str(),1);
    signal = (MnvH1D*)f_mc->Get(var.c_str());
    signal->SetTitle("Signal");
    signal->SetLineColor(kGreen);
    signal->SetFillColor(kGreen);

    // Get Stat1s
    max_bin = signal->GetMaximumBin();
    hist_max = hist_max + signal->GetBinContent(max_bin);
    bin_width = signal->GetBinWidth(1);

    // Add to TObjArray
    mc_hists->Add(signal);

    // ------------------------------------------------------------------------
    // Plot  - If you want A Log Plot, axis_minimum = 0.1
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    ApplyStyle(plotter);
    TCanvas* c = new TCanvas(var_name.c_str(),var_name.c_str(),1280,800);
    plotter->DrawStackedMC(mc_hists,1,"N",0);

    // If Cut Histogram - Add Cut Arrows
    if (nCutArrows == 1){
        AddCutArrow(plotter, cutArrow1, hist_max, bin_width);
    }else if (nCutArrows == 2){
        AddCutArrow(plotter, cutArrow1, hist_max, bin_width);
        AddCutArrow(plotter, cutArrow2, hist_max, bin_width);
    }

    // Add Pi0 Mass Line
    std::size_t found = var_name.find("invMass");
    if (found != std::string::npos){
        TArrow arrow;
        arrow.SetLineWidth(2);
        arrow.SetLineColor(kBlue);
        arrow.DrawArrow(134.98,hist_max/8,134.98,0, 0.01); 
    }

    // Add Delta+ InvMass Line
    if (    var_name.compare("deltaInvMass") == 0 ||
            var_name.compare("W_Calc") == 0 ){
        TLine deltaMass;
        deltaMass.SetLineWidth(2);
        deltaMass.SetLineColor(kBlue);
        deltaMass.DrawLine(1.232,0,1.232,hist_max);
    }

    // ------------------------------------------------------------------------
    // Legend 
    // ------------------------------------------------------------------------
    //double leg_x_min = 0.60;
    //double leg_x_max = 0.90;
    //double leg_y_min = 0.80;
    //double leg_y_max = 0.85;

    //TLegend *legend = new TLegend(leg_x_min, leg_y_min, leg_x_max, leg_y_max);
    //ApplyStyle_Legend(legend);
    ////legend->SetTextSize(0.04);
    //legend->AddEntry(signal, "Signal","f");
    //legend->Draw();


    // Print Plot
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),"_mc_signal.pdf"), "pdf");

    delete c;
    delete plotter;
}

void CCProtonPi0_Plotter::DrawNormalizedMigrationHistogram(rootDir &dir, std::string var_name, std::string plotDir)
{
    std::string root_dir = dir.mc;

    // Get Histogram
    TFile* f = new TFile(root_dir.c_str());
    MnvH2D* hist2D = (MnvH2D*)f->Get(var_name.c_str());

    // Canvas
    Double_t w = 800; 
    Double_t h = 800;
    TCanvas* c = new TCanvas("c","c",w,h);

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    //ApplyStyle(plotter);
    plotter->DrawNormalizedMigrationHistogram(hist2D, true, false, false);

    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),"_migration.pdf"), "pdf");

    delete c;
    delete plotter;
    delete f;
}

void CCProtonPi0_Plotter::DrawTGraph(rootDir &dir, std::string var_name, std::string plotDir)
{
    std::string root_dir = dir.mc;

    // Get Histogram
    TFile* f = new TFile(root_dir.c_str());
    TGraph* graph = (TGraph*)f->Get(var_name.c_str());

    // Canvas
    TCanvas* c = new TCanvas("c","c",1280,800);

    // Plot 
    graph->Draw();
    c->Update();
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".pdf"), "pdf");

    delete c;
    delete f;
}

void CCProtonPi0_Plotter::DrawBackgroundSubtraction(bool isMC)
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    TFile* f;
    if (isMC) f = new TFile(rootDir_CrossSection.mc.c_str());
    else f = new TFile(rootDir_CrossSection.data.c_str());

    double hist_max = 0;
    double area_data = 0.0;
    double area_mc = 0.0;
    double max_bin;
    double bin_width;

    // ------------------------------------------------------------------------
    // Get Data Histogram
    // ------------------------------------------------------------------------
    MnvH1D* data;
    if (isMC) data = (MnvH1D*)f->Get("invMass_mc_reco_all"); 
    else data = (MnvH1D*)f->Get("invMass_all"); 
    max_bin = data->GetMaximumBin();
    hist_max = (hist_max + data->GetBinContent(max_bin))*1.2;
    bin_width = data->GetBinWidth(1);
    area_data = data->Integral();

    // ------------------------------------------------------------------------
    // Fill TObjArray - For MC Histograms
    // ------------------------------------------------------------------------
    TObjArray* mc_hists = new TObjArray;
    MnvH1D* temp;

    // Get All Background
    temp = (MnvH1D*)f->Get("invMass_mc_reco_bckg");
    temp->SetTitle("Background");
    area_mc = area_mc + temp->Integral();
    mc_hists->Add(temp);

    // Get Signal
    temp = (MnvH1D*)f->Get("invMass_mc_reco_signal");
    temp->SetTitle("Signal");
    area_mc = area_mc + temp->Integral();
    mc_hists->Add(temp);

    // ------------------------------------------------------------------------
    // MC Normalization 
    // ------------------------------------------------------------------------
    double mc_ratio = area_data/area_mc;

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas("c","c",1280,800);
    ApplyStyle(plotter);
    plotter->axis_minimum = 0.0;
    plotter->DrawDataStackedMC(data,mc_hists,mc_ratio,"TR","Data",2,1);

    TH1D* fit_result = (TH1D*)f->Get("data_fit_result");

    // Plot Options
    fit_result->SetLineColor(kBlue);
    fit_result->SetLineWidth(3);
    fit_result->SetFillStyle(0);

    fit_result->Draw("SAME");

    // Add Plot Labels
    plotter->AddHistoTitle(data->GetTitle());
    const double x_pos = 0.65;
    const double y_pos = 0.7;
    const double text_size = 0.03;
    area_data = data->Integral();
    area_mc = area_mc * mc_ratio;
    double area_fit = fit_result->Integral();
    plotter->AddPlotLabel(Form("Area(Data) = %3.2f",area_data), x_pos, y_pos, text_size, kBlue, 62, 1); 
    plotter->AddPlotLabel(Form("Area(MC) = %3.2f",area_mc), x_pos, y_pos-text_size, text_size, kBlue, 62, 1); 
    plotter->AddPlotLabel(Form("Area(Fit) = %3.2f",area_fit), x_pos, y_pos-(2*text_size), text_size, kBlue, 62, 1); 

    // Add Cut Arrows
    CutArrow pi0invMass_min(60,"R"); 
    CutArrow pi0invMass_max(200,"L"); 
    AddCutArrow(plotter, pi0invMass_min, hist_max, bin_width);
    AddCutArrow(plotter, pi0invMass_max, hist_max, bin_width);

    // Add Pi0 InvMass Line
    TLine pi0Mass;
    pi0Mass.SetLineWidth(2);
    pi0Mass.SetLineColor(kBlue);
    pi0Mass.DrawLine(134.98,0,134.98,hist_max);


    // Print Plot
    std::string tag;
    if(isMC) tag = "_MC";
    else tag = "_Data";
    c->Print(Form("%s%s%s%s",plotDir.c_str(),"Background_Subtraction",tag.c_str(),".pdf"), "pdf");

    delete c;
    delete plotter;
}

void CCProtonPi0_Plotter::NormalizeToNormBinWidth(MnvH1D* hist)
{
    double norm_bin_width = hist->GetNormBinWidth();
    hist->Scale(norm_bin_width,"width");
}


double CCProtonPi0_Plotter::GetSmallestBinWidth(MnvH1D* hist)
{
    double smallest = 99999999;
    int nBins = hist->GetNbinsX();
    for (int i = 0; i <= nBins; ++i){
        double current = hist->GetBinWidth(i);
        if (current < smallest) smallest = current;
    }

    return smallest;
}

TH1D* CCProtonPi0_Plotter::GetBinNormalizedTH1D(MnvH1D* hist, bool WithSystError)
{
    MnvH1D* tempHist = new MnvH1D(*hist);
    TH1D* h = NULL;

    NormalizeToNormBinWidth(tempHist);

    if (WithSystError){
        h = new TH1D(tempHist->GetCVHistoWithError());
    }else{
        h = new TH1D(tempHist->GetCVHistoWithStatError());
    }
    delete tempHist;
    return h;
}

void CCProtonPi0_Plotter::DrawDataStackedMC_WithSignalTypes(rootDir &dir, std::string var_name, std::string plotDir)
{
    std::string var;

    double mc_ratio = POT_ratio;

    // Open ROOT Files
    TFile* f_data = new TFile(dir.data.c_str());
    TFile* f_mc = new TFile(dir.mc.c_str());

    // Get Data Histogram
    var = var_name + "_0";
    TH1D* h_data = GetBinNormalizedTH1D(GetMnvH1D(f_data,var)); 
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(1);
    h_data->SetMarkerColor(kBlack);
    h_data->SetLineWidth(2);
    h_data->SetLineColor(kBlack);

    // Get Signal: Delta RES
    var = var_name + "_7";
    TH1D* h_signal_delta_res = GetBinNormalizedTH1D(GetMnvH1D(f_mc,var));
    h_signal_delta_res->Scale(mc_ratio);
    h_signal_delta_res->SetLineWidth(1);
    h_signal_delta_res->SetLineColor(kGreen);
    h_signal_delta_res->SetFillColor(kGreen);
    h_signal_delta_res->SetFillStyle(3001);

    // Get Signal: Other RES
    var = var_name + "_8";
    TH1D* h_signal_other_res = GetBinNormalizedTH1D(GetMnvH1D(f_mc,var));
    h_signal_other_res->Scale(mc_ratio);
    h_signal_other_res->SetLineWidth(1);
    h_signal_other_res->SetLineColor(kGreen+2);
    h_signal_other_res->SetFillColor(kGreen+2);
    h_signal_other_res->SetFillStyle(3001);

    // Get Signal: Non-RES
    var = var_name + "_9";
    TH1D* h_signal_non_res = GetBinNormalizedTH1D(GetMnvH1D(f_mc,var));
    h_signal_non_res->Scale(mc_ratio);
    h_signal_non_res->SetLineWidth(1);
    h_signal_non_res->SetLineColor(kGreen+4);
    h_signal_non_res->SetFillColor(kGreen+4);
    h_signal_non_res->SetFillStyle(3001);
  
    // Get Background: With Pi0
    var = var_name + "_3";
    TH1D* h_bckg_with_pi0 = GetBinNormalizedTH1D(GetMnvH1D(f_mc,var));
    h_bckg_with_pi0->Scale(mc_ratio);
    h_bckg_with_pi0->SetLineWidth(1);
    h_bckg_with_pi0->SetLineColor(kRed);
    h_bckg_with_pi0->SetFillColor(kRed);
    h_bckg_with_pi0->SetFillStyle(3001);

    // Get Background: QELike 
    var = var_name + "_4";
    TH1D* h_bckg_qelike = GetBinNormalizedTH1D(GetMnvH1D(f_mc,var));
    h_bckg_qelike->Scale(mc_ratio);
    h_bckg_qelike->SetLineWidth(1);
    h_bckg_qelike->SetLineColor(kOrange);
    h_bckg_qelike->SetFillColor(kOrange);
    h_bckg_qelike->SetFillStyle(3001);

    // Get Background: With PiPlus 
    var = var_name + "_5";
    TH1D* h_bckg_piplus = GetBinNormalizedTH1D(GetMnvH1D(f_mc,var));
    h_bckg_piplus->Scale(mc_ratio);
    h_bckg_piplus->SetLineWidth(1);
    h_bckg_piplus->SetLineColor(kBlue);
    h_bckg_piplus->SetFillColor(kBlue);
    h_bckg_piplus->SetFillStyle(3001);

    // Get Background: With PiPlus 
    var = var_name + "_6";
    TH1D* h_bckg_other = GetBinNormalizedTH1D(GetMnvH1D(f_mc,var));
    h_bckg_other->Scale(mc_ratio);
    h_bckg_other->SetLineWidth(1);
    h_bckg_other->SetLineColor(kGray);
    h_bckg_other->SetFillColor(kGray);
    h_bckg_other->SetFillStyle(3001);

    TCanvas* c = new TCanvas("c","c",800,800);

    THStack *hs = new THStack("hs",var_name.c_str());
    hs->Add(h_signal_delta_res);
    hs->Add(h_signal_other_res);
    hs->Add(h_signal_non_res);
    hs->Add(h_bckg_with_pi0);
    hs->Add(h_bckg_qelike);
    hs->Add(h_bckg_piplus);
    hs->Add(h_bckg_other);

    double hs_max = hs->GetMaximum();
    hs->SetMaximum(hs_max * 1.5);
    hs->Draw("HIST");
    hs->SetTitle(h_data->GetTitle());
    hs->GetXaxis()->SetTitle(h_data->GetXaxis()->GetTitle());
    hs->GetYaxis()->SetTitle(h_data->GetYaxis()->GetTitle());

    if (thesisStyle){
        hs->GetXaxis()->SetTitleFont(62);
        hs->GetXaxis()->SetTitleSize(0.06);
        hs->GetXaxis()->CenterTitle();
        hs->GetXaxis()->SetTitleOffset(1.15);
        hs->GetXaxis()->SetLabelFont(42);
        hs->GetXaxis()->SetLabelSize(0.05);
        hs->GetXaxis()->SetNdivisions(408);

        hs->GetYaxis()->SetTitleFont(62);
        hs->GetYaxis()->SetTitleSize(0.06);
        //hs->GetYaxis()->CenterTitle();
        hs->GetYaxis()->SetTitleOffset(1.2);
        hs->GetYaxis()->SetLabelFont(42);
        hs->GetYaxis()->SetLabelSize(0.05);
        TGaxis::SetMaxDigits(3);
    }

    h_data->Draw("SAME E1 X0");

    // ------------------------------------------------------------------------
    // Plot Labels 
    // ------------------------------------------------------------------------
    TLegend *legend = new TLegend(0.6,0.6,0.9,0.9);  
    ApplyStyle_Legend(legend);    
    legend->AddEntry(h_data, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(h_signal_delta_res, "Signal: #Delta Res","f");
    legend->AddEntry(h_signal_other_res, "Signal: Other Res","f");
    legend->AddEntry(h_signal_non_res, "Signal: Non-Res","f");
    legend->AddEntry(h_bckg_with_pi0, "Bckg: #pi^{0} + X","f");
    legend->AddEntry(h_bckg_qelike, "Bckg: QE Like","f");
    legend->AddEntry(h_bckg_piplus, "Bckg: With Charged #pi","f");
    legend->AddEntry(h_bckg_other, "Bckg: Other","f");
    legend->Draw();

    // Add Text
    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.03);
    text.DrawLatex(0.2, 0.85, "#color[4]{POT Normalized}");

    std::size_t found = var_name.find("W");
    if (found != std::string::npos){
        TLine deltaMass;
        deltaMass.SetLineWidth(3);
        deltaMass.SetLineColor(kBlue);
        deltaMass.DrawLine(1.232,0,1.232,hs_max * 1.2);
    }

    // Plot Output
    gStyle->SetEndErrorSize(6);
    gStyle->SetOptStat(0); 
    c->Update();
    std::string out_name;
    out_name = plotDir + var_name + "_xsec_IntType" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete legend;
    delete c;

    delete f_data;
    delete f_mc;
}

// Returns a "new" MnvH1D -- Do not forget to delete
MnvH1D* CCProtonPi0_Plotter::GetBckgSubtractedData(rootDir& dir, std::string var_name, double nBckg)
{
    TFile* f_data = new TFile(dir.data.c_str());
    TFile* f_mc  = new TFile(dir.mc.c_str());

    // Get Data Histogram
    std::string var = var_name + "_all";
    MnvH1D* h_data = GetMnvH1D(f_data,var); 

    // Get MC Background and Normalize
    var = var_name + "_mc_reco_bckg";
    MnvH1D* h_mc_bckg = GetMnvH1D(f_mc,var);
    NormalizeHistogram(h_mc_bckg);
    h_mc_bckg->Scale(nBckg);
    
    // Subtract Background
    h_data->Add(h_mc_bckg, -1);

    delete h_mc_bckg;
    delete f_data;
    delete f_mc;

    return h_data;
}

// Returns a "new" MnvH2D -- Do not forget to delete
MnvH2D* CCProtonPi0_Plotter::GetBckgSubtractedData_2D(rootDir& dir, std::string var_name, double nBckg)
{
    TFile* f_data = new TFile(dir.data.c_str());
    TFile* f_mc  = new TFile(dir.mc.c_str());

    // Get Data Histogram
    std::string var = var_name + "_all";
    MnvH2D* h_data = GetMnvH2D(f_data,var); 

    // Get MC Background and Normalize
    var = var_name + "_mc_reco_bckg";
    MnvH2D* h_mc_bckg = GetMnvH2D(f_mc,var);
    NormalizeHistogram(h_mc_bckg);
    h_mc_bckg->Scale(nBckg);
    
    // Subtract Background
    h_data->Add(h_mc_bckg, -1);

    delete h_mc_bckg;
    delete f_data;
    delete f_mc;

    return h_data;
}

void CCProtonPi0_Plotter::DrawDataStackedMC_Signal(rootDir &dir, std::string var_name, std::string plotDir, double nBckg)
{
    std::string var;

    double mc_ratio = POT_ratio;

    // Open ROOT Files
    TFile* f_data = new TFile(dir.data.c_str());
    TFile* f_mc = new TFile(dir.mc.c_str());

    // Get Bckg Subtracted Data Histogram
    MnvH1D* h_data = GetBckgSubtractedData(dir, var_name, nBckg); 
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(1);
    h_data->SetMarkerColor(kBlack);
    h_data->SetLineWidth(2);
    h_data->SetLineColor(kBlack);

    // Get Signal: Delta RES
    var = var_name + "_7";
    MnvH1D* h_signal_delta_res = GetMnvH1D(f_mc,var);
    h_signal_delta_res->Scale(mc_ratio);
    h_signal_delta_res->SetLineWidth(1);
    h_signal_delta_res->SetLineColor(kGray+1);
    h_signal_delta_res->SetFillColor(kGray+1);
    h_signal_delta_res->SetFillStyle(1001);

    // Get Signal: Other RES
    var = var_name + "_8";
    MnvH1D* h_signal_other_res = GetMnvH1D(f_mc,var);
    h_signal_other_res->Scale(mc_ratio);
    h_signal_other_res->SetLineWidth(1);
    h_signal_other_res->SetLineColor(kGreen+2);
    h_signal_other_res->SetFillColor(kGreen+2);
    h_signal_other_res->SetFillStyle(1001);

    // Get Signal: Non-RES
    var = var_name + "_9";
    MnvH1D* h_signal_non_res = GetMnvH1D(f_mc,var);
    h_signal_non_res->Scale(mc_ratio);
    h_signal_non_res->SetLineWidth(1);
    h_signal_non_res->SetLineColor(kRed+1);
    h_signal_non_res->SetFillColor(kRed+1);
    h_signal_non_res->SetFillStyle(1001);

    TCanvas* c = new TCanvas("c","c",1280,800);

    THStack *hs = new THStack("hs",var_name.c_str());
    hs->Add(h_signal_delta_res);
    hs->Add(h_signal_other_res);
    hs->Add(h_signal_non_res);

    double hs_max = hs->GetMaximum();
    hs->SetMaximum(hs_max * 1.5);
    hs->Draw("HIST");
    hs->SetTitle(h_data->GetTitle());
    hs->GetXaxis()->SetTitle(h_data->GetXaxis()->GetTitle());
    hs->GetYaxis()->SetTitle(h_data->GetYaxis()->GetTitle());

    h_data->Draw("SAME E1 X0");

    // ------------------------------------------------------------------------
    // Plot Labels 
    // ------------------------------------------------------------------------
    TLegend *legend = new TLegend(0.6,0.75,0.9,0.9);  
    ApplyStyle_Legend(legend);    
    legend->AddEntry(h_data, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(h_signal_delta_res, "Signal: #Delta Res","f");
    legend->AddEntry(h_signal_other_res, "Signal: Other Res","f");
    legend->AddEntry(h_signal_non_res, "Signal: Non-Res","f");
    legend->Draw();

    // Add Text
    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.03);

    text.DrawLatex(0.2, 0.85, "#color[4]{POT Normalized}");

    std::size_t found = var_name.find("W");
    if (found != std::string::npos){
        TLine deltaMass;
        deltaMass.SetLineWidth(3);
        deltaMass.SetLineColor(kBlue);
        deltaMass.DrawLine(1.232,0,1.232,hs_max * 1.2);
    }

    // Plot Output
    gStyle->SetEndErrorSize(6);
    gStyle->SetOptStat(0); 
    c->Update();
    std::string out_name;
    out_name = plotDir + var_name + "_Signal" + ".pdf"; 

    c->Print(out_name.c_str(),"pdf");

    delete legend;
    delete c;

    delete f_data;
    delete f_mc;
}

void CCProtonPi0_Plotter::DrawDataMCSignal_Diff(rootDir& dir, std::string var_name, std::string plotDir, double nBckg)
{
    std::string rootDir_mc = dir.mc;
    std::string rootDir_data = dir.data;

    TFile* f_mc = new TFile(rootDir_mc.c_str());

    std::string var = Form("%s_%d",var_name.c_str(),1);

    MnvH1D* data = GetBckgSubtractedData(dir, var_name, nBckg);
    MnvH1D* mc = GetMnvH1D(f_mc, var);
    mc->Scale(POT_ratio);
 
    MnvH1D* diff = new MnvH1D(*data);

    diff->Add(mc,-1);
    MnvH1D* ratio = new MnvH1D(*diff); 
    ratio->Scale(1/mc->Integral());
    ratio->SetLineWidth(3);
    ratio->SetLineColor(kRed);
    ratio->SetFillColor(kWhite);
    //ratio->SetFillStyle(3010);
    ratio->SetFillStyle(0);

    ratio->SetMinimum(0.0);
    ratio->GetYaxis()->SetTitle("(Data-MC)/Area(MC)");

    TCanvas* c = new TCanvas("c","c",1280,800);
    ratio->Draw("HIST");

    // Save Plot 
    //gStyle->SetOptStat(0); 
    c->Update();
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),"_Ratio.pdf"), "pdf");

    delete c;
    delete data;
    delete mc;
    delete ratio;
    delete diff;
    delete f_mc;
}

void CCProtonPi0_Plotter::DrawDataMCSignal_Diff_2D(rootDir& dir, std::string var_name, std::string plotDir, double nBckg)
{
    std::string rootDir_mc = dir.mc;
    std::string rootDir_data = dir.data;

    TFile* f_mc = new TFile(rootDir_mc.c_str());

    std::string var = Form("%s_%s",var_name.c_str(),"mc_reco_signal");

    MnvH2D* data = GetBckgSubtractedData_2D(dir, var_name, nBckg);
    MnvH2D* mc = GetMnvH2D(f_mc, var);
    mc->Scale(POT_ratio);

    data->SetMinimum(0.0);
    mc->SetMinimum(0.0);
    DrawMnvH2D(data, var_name, plotDir, false);
    DrawMnvH2D(mc, var_name, plotDir, true);

    MnvH2D* diff = new MnvH2D (*data);
    diff->Add(mc,-1);
    diff->SetMinimum(0.0);
 
    var = var_name + "_Diff";
    DrawMnvH2D(diff, var, plotDir, false);

    TH2D* errors = new TH2D(data->GetTotalError(true,true));


    std::cout<<"Printing Data"<<std::endl;
    printBins(data,var_name);
    std::cout<<"Printing MC"<<std::endl;
    printBins(mc,var_name);
    std::cout<<"Printing Diff"<<std::endl;
    printBins(diff,var_name);
    std::cout<<"Printing Errors"<<std::endl;
    printBins(errors,var_name);


    delete diff;
    delete data;
    delete errors;
    delete mc;
    delete f_mc;
}

double CCProtonPi0_Plotter::Get2DTotalFlow(MnvH2D* h)
{
    int nBins_X = h->GetNbinsX();
    int nBins_Y = h->GetNbinsY();

    double nTotalFlow = 0.0;
    // First Count Overflow on X
    for (int y = 1; y <= nBins_Y; ++y ){
        nTotalFlow += h->GetBinContent(nBins_X+1, y);
    }
    for (int y = 1; y <= nBins_Y; ++y ){
        nTotalFlow += h->GetBinContent(0, y);
    }

    // Second Count Overflow on Y
    for (int x = 1; x <= nBins_X; ++x ){
        nTotalFlow += h->GetBinContent(x, nBins_Y+1);
    }
    for (int x = 1; x <= nBins_X; ++x){
        nTotalFlow += h->GetBinContent(x, 0);
    }

    // Third Count 4-Corners
    nTotalFlow += h->GetBinContent(0,0);
    nTotalFlow += h->GetBinContent(0,nBins_Y+1);
    nTotalFlow += h->GetBinContent(nBins_X+1,nBins_Y+1);
    nTotalFlow += h->GetBinContent(nBins_X+1,0);

    return nTotalFlow;
}

void CCProtonPi0_Plotter::NormalizeHistogram(MnvH2D* h)
{
    std::cout<<"Normalizing MnvH2D"<<std::endl;
    double area = h->Integral();
    std::cout<<"\tArea Before = "<<area<<std::endl;
    double nTotalFlow = Get2DTotalFlow(h);
    h->Scale(1/(area+nTotalFlow)); // Scale only on CentralValue
    std::cout<<"\tArea After= "<<h->Integral()<<std::endl;
}

void CCProtonPi0_Plotter::NormalizeHistogram(MnvH1D* h)
{
    int NBins = h->GetNbinsX();
    double area = h->Integral();
    double nOverFlow = h->GetBinContent(NBins+1);
    double nUnderFlow = h->GetBinContent(0);
    h->Scale(1/(area+nOverFlow+nUnderFlow),"",false); // Scale only on CentralValue
}

void CCProtonPi0_Plotter::NormalizeHistogram(TH1D* h)
{
    int NBins = h->GetNbinsX();
    double area = h->Integral();
    double nOverFlow = h->GetBinContent(NBins+1);
    double nUnderFlow = h->GetBinContent(0);
    h->Scale(1/(area+nOverFlow+nUnderFlow)); // Scale only on CentralValue
}

void CCProtonPi0_Plotter::ApplyStyle_Thesis()
{
    // Remove Stat Box
    gStyle->SetOptStat(0); 
    
    // Remove Title 
    gStyle->SetOptTitle(0); 

    // Canvas Styles
    gStyle->SetCanvasColor(0);
    gStyle->SetPadBorderMode(0);
    gStyle->SetFrameBorderMode(0);
    gStyle->SetCanvasBorderMode(0);
    gStyle->SetPadTopMargin(0.09);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetFrameLineWidth(2);
    gStyle->SetHistLineWidth(2);

    //gStyle->SetEndErrorSize(2);
}

void CCProtonPi0_Plotter::ApplyStyle_Legend(TLegend* legend)
{
    legend->SetBorderSize(0);
    legend->SetFillColor(-1);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.04);
}


#endif


