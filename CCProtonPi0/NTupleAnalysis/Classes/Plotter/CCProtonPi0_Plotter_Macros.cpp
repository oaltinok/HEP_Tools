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
    DrawStackedMC_BckgType(dir, var_name, plotDir, nCutArrows, cutArrow1, cutArrow2);
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
    max_bin = data->GetMaximumBin();
    hist_max = (hist_max + data->GetBinContent(max_bin))*1.2;
    bin_width = data->GetBinWidth(1);

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
    MnvH1D* temp;

    // Get All Background
    var = Form("%s_%d",var_name.c_str(),2);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Background");
    mc_hists->Add(temp);
    area_mc += temp->Integral() * mc_ratio;
    
    // Get Signal
    var = Form("%s_%d",var_name.c_str(),1);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Signal");
    mc_hists->Add(temp);
    area_mc += temp->Integral() * mc_ratio;

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas("c","c",1280,800);
    ApplyStyle(plotter);
    //plotter->axis_minimum = 0.1;
    plotter->axis_minimum = 0.0;
    plotter->DrawDataStackedMC(data,mc_hists,mc_ratio,"TR","Data",2,1);

    // Add Plot Labels
    plotter->AddHistoTitle(data[0].GetTitle());
    AddNormBox(plotter, isPOTNorm, mc_ratio);

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
        TLine pi0Mass;
        pi0Mass.SetLineWidth(2);
        pi0Mass.SetLineColor(kBlue);
        pi0Mass.DrawLine(134.98,0,134.98,hist_max);
    }

    // Add Delta+ InvMass Line
    if (    var_name.compare("deltaInvMass") == 0 ||
            var_name.compare("W_Calc") == 0 ){
        TLine deltaMass;
        deltaMass.SetLineColor(kBlue);
        deltaMass.DrawLine(1.232,0,1.232,hist_max);
    }

    // Add Areas
    const double y_pos = 0.88;
    const double text_size = 0.03;
    double area_data = data->Integral();
    plotter->AddPlotLabel(Form("Area(Data)/Area(MC) = %3.2f",area_data/area_mc),0.3,y_pos,text_size,kBlue); 

    // Print Plot
    c->Print(Form("%s%s%s%s%s",plotDir.c_str(),var_name.c_str(),"_bckg_all_",norm_label.c_str(),".png"), "png");

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
    c->Print(Form("%s%s%s%s%s",plotDir.c_str(),var_name.c_str(),"_bckg_Type_",norm_label.c_str(),".png"), "png");

    delete c;
    delete plotter;
}

void CCProtonPi0_Plotter::ApplyStyle(MnvPlotter* plotter)
{
    plotter->hist_min_zero = false;
    plotter->mc_line_width = 0;
    plotter->axis_title_size_x = 0.04;
    plotter->axis_title_size_y = 0.04;
    plotter->axis_label_size = 0.03;
    plotter->legend_text_size = 0.02;
}

void CCProtonPi0_Plotter::ApplyStyle_Errors(MnvPlotter* plotter)
{
    //plotter->hist_min_zero = false;
    //plotter->mc_line_width = 0;
    //plotter->axis_title_size_x = 0.04;
    //plotter->axis_title_size_y = 0.04;
    //plotter->axis_label_size = 0.03;
    plotter->legend_text_size = 0.02;
    plotter->legend_n_columns = 2;
    plotter->height_nspaces_per_hist = 0.5;
    plotter->width_xspace_per_letter = 0.1;
    //-- define colors of the standard errors
    plotter->error_color_map.clear();
    plotter->error_summary_group_map.clear();

    plotter->error_color_map["GENIE"] = kGreen+2;

    std::vector<std::string> genieGroup;
    genieGroup.push_back("GENIE_AGKYxF1pi"         );
    genieGroup.push_back("GENIE_AhtBY"             );
    genieGroup.push_back("GENIE_BhtBY"             );
    genieGroup.push_back("GENIE_CCQEPauliSupViaKF" );
    genieGroup.push_back("GENIE_CV1uBY"            );
    genieGroup.push_back("GENIE_CV2uBY"            );
    genieGroup.push_back("GENIE_EtaNCEL"           );
    genieGroup.push_back("GENIE_FrAbs_N"           );
    genieGroup.push_back("GENIE_FrAbs_pi"          );
    genieGroup.push_back("GENIE_FrCEx_N"           );
    genieGroup.push_back("GENIE_FrCEx_pi"          );
    genieGroup.push_back("GENIE_FrElas_N"          );
    genieGroup.push_back("GENIE_FrElas_pi"         );
    genieGroup.push_back("GENIE_FrInel_N"          );
    genieGroup.push_back("GENIE_FrInel_pi"         );
    genieGroup.push_back("GENIE_FrPiProd_N"        );
    genieGroup.push_back("GENIE_FrPiProd_pi"       );
    genieGroup.push_back("GENIE_MFP_N"             );
    genieGroup.push_back("GENIE_MFP_pi"            );
    genieGroup.push_back("GENIE_MaCCQE"            );
    genieGroup.push_back("GENIE_MaCCQEshape"       );
    genieGroup.push_back("GENIE_MaNCEL"            );
    genieGroup.push_back("GENIE_MaRES"             );
    genieGroup.push_back("GENIE_MvRES"             );
    genieGroup.push_back("GENIE_NormCCQE"          );
    genieGroup.push_back("GENIE_NormCCRES"         );
    genieGroup.push_back("GENIE_NormDISCC"         );
    genieGroup.push_back("GENIE_NormNCRES"         );
    genieGroup.push_back("GENIE_RDecBR1gamma"      );
    genieGroup.push_back("GENIE_Rvn1pi"            );
    genieGroup.push_back("GENIE_Rvn2pi"            );
    genieGroup.push_back("GENIE_Rvp1pi"            );
    genieGroup.push_back("GENIE_Rvp2pi"            );
    genieGroup.push_back("GENIE_Theta_Delta2Npi"   );
    genieGroup.push_back("GENIE_VecFFCCQEshape"    );
    plotter->error_summary_group_map["GENIE"] = genieGroup;   
}

void CCProtonPi0_Plotter::DrawErrorSummary(MnvH1D* hist, std::string var_name, std::string plotDir)
{
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas("c","c",1280,800);

    ApplyStyle_Errors(plotter);

    plotter->DrawErrorSummary(hist);

    // Add Plot Labels
    plotter->AddHistoTitle(hist->GetTitle());

    // Print Plot
    std::string out_name = plotDir + "Errors_" + var_name + ".png";
    c->Print(out_name.c_str(), "png");

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
    c->Print(Form("%s%s%s%s%s",plotDir.c_str(),var_name.c_str(),"_Error_", error_name.c_str() ,".png"), "png");

    delete c;
    delete plotter;
}

void CCProtonPi0_Plotter::DrawDataMC(rootDir& dir, std::string var_name, std::string plotDir)
{
    std::string rootDir_mc = dir.mc;
    std::string rootDir_data = dir.data;

    TFile* f_mc = new TFile(rootDir_mc.c_str());
    TFile* f_data = new TFile(rootDir_data.c_str());

    //std::string var = Form("%s_%d",var_name.c_str(),0);
    std::string var = var_name;

    // POT Normalized
    MnvH1D* mc = (MnvH1D*)f_mc->Get(var.c_str());
    MnvH1D* data = (MnvH1D*)f_data->Get(var.c_str()); 
    DrawDataMC(data, mc, var_name, plotDir, true);

    // Area Normalized
    mc = (MnvH1D*)f_mc->Get(var.c_str());
    data = (MnvH1D*)f_data->Get(var.c_str()); 
    DrawDataMC(data, mc, var_name, plotDir, false);
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
    
    // Normalize to Norm Bin Width
    NormalizeToNormBinWidth(tempData);
    NormalizeToNormBinWidth(tempMC);

    //printBins(tempData,"data",false);
    //printBins(tempMC,"MC",true);

    // ------------------------------------------------------------------------
    // Neutrino Energy Comparison Only
    //tempData->SetMaximum(80);
    //tempMC->SetMaximum(80);
    //tempData->GetXaxis()->SetRangeUser(0,10);
    //tempMC->GetXaxis()->SetRangeUser(0,10);
    // ------------------------------------------------------------------------


    plotter->DrawDataMC(tempData, tempMC, mc_ratio, "TR", false);

    // Add Plot Labels
    plotter->AddHistoTitle(data->GetTitle());
    AddNormBox(plotter, isPOTNorm, mc_ratio);

    // Add Areas
    const double y_pos = 0.88;
    const double text_size = 0.03;
    double area_data = data->Integral("width");
    double area_mc = mc->Integral("width") * mc_ratio;
    plotter->AddPlotLabel(Form("Area(Data)/Area(MC) = %3.2f",area_data/area_mc),0.3,y_pos,text_size,kBlue); 


    // Add X = 0 Line 
    TLine line_0;
    line_0.SetLineWidth(2);
    line_0.SetLineStyle(7);
    line_0.SetLineColor(kBlue);
    line_0.DrawLine(0,0,0,250);

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
    h_mc_total->Scale(mc_ratio);
    
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
    out_name = plotDir + var_name + "_" + norm_label + ".png"; 

    c->Print(out_name.c_str(),"png");

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
    c->Print(Form("%s%s%s%s%s",plotDir.c_str(),var_name.c_str(),"_ratio_",norm_label.c_str(),".png"), "png");

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
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".png"), "png");

    delete c;
    delete plotter;
}

void CCProtonPi0_Plotter::DrawMnvH1D(MnvH1D* hist1D, std::string var_name, std::string plotDir)
{
    // Create Canvas
    TCanvas* c = new TCanvas("c","c",1280,800);

    // Plot Options
    hist1D->SetLineColor(kRed);
    hist1D->SetLineWidth(3);
    hist1D->SetFillColor(kRed);
    hist1D->SetFillStyle(3010);

    double norm_bin_width = hist1D->GetNormBinWidth();
    hist1D->Scale(norm_bin_width,"width");

    hist1D->Draw("HIST");
    gPad->Update();
    gStyle->SetOptStat(111111); 

    // Find Peak
    int max_bin = hist1D->GetMaximumBin();
    double max_bin_location = hist1D->GetBinCenter(max_bin);
    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.03);
    text.DrawLatex(0.15,0.85,Form("%s%3.2f", "Peak at ",max_bin_location));

    // Add Pi0 InvMass Line
    std::size_t found = var_name.find("invMass");
    if (found != std::string::npos){
        int max_value = hist1D->GetBinContent(max_bin);
        TLine pi0Mass;
        pi0Mass.SetLineWidth(2);
        pi0Mass.SetLineColor(kBlue);
        pi0Mass.DrawLine(134.98,0,134.98,max_value);
    }

    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".png"), "png");
    delete c;
    delete hist1D;
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
    // Create Canvas
    TCanvas* c = new TCanvas("c","c",1280,800);
    if(isLogScale) c->SetLogy();

    hist1D->Scale(1,"width");

    // Plot Options
    hist1D->SetLineColor(kRed);
    hist1D->SetLineWidth(3);
    hist1D->SetFillColor(kRed);
    hist1D->SetFillStyle(3010);

    hist1D->Draw("HIST");
    gPad->Update();
    gStyle->SetOptStat(111111); 

    // Find Peak
    int max_bin = hist1D->GetMaximumBin();
    int max_value = hist1D->GetBinContent(max_bin);
    double max_bin_location = hist1D->GetBinCenter(max_bin);
    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.03);
    text.DrawLatex(0.15,0.85,Form("%s%3.2f", "Peak at ",max_bin_location));

    // Add Pi0 InvMass Line
    std::size_t found = var_name.find("invMass");
    if (found != std::string::npos){
        TLine pi0Mass;
        pi0Mass.SetLineWidth(2);
        pi0Mass.SetLineColor(kBlue);
        pi0Mass.DrawLine(134.98,0,134.98,max_value);
    }

    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".png"), "png");
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

    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".png"), "png");

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
    c->SetWindowSize(w,h);

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
    hist2D->GetYaxis()->SetTitleOffset(1.8);
    hist2D->Draw("colz");
    gPad->Update();

    //double line_min = hist2D->GetXaxis()->GetBinLowEdge(1);
    //double line_max = hist2D->GetXaxis()->GetBinLowEdge(nBinsX);
    //TLine xy;
    //xy.SetLineWidth(2);
    //xy.SetLineColor(kBlack);
    //xy.DrawLine(line_min,line_min,line_max,line_max);

    //TLine fit;
    //fit.SetLineWidth(2);
    //fit.SetLineColor(kRed);
    //fit.DrawLine(0,80.9,500,585.9);

    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".png"), "png");

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

    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".png"), "png");

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

    plotter->AddCutArrow(cut_location, ymin, ymax, arrow_length, arrow_direction); 
}

void CCProtonPi0_Plotter::DrawStackedMC_BckgAll(rootDir &dir, std::string var_name, std::string plotDir, int nCutArrows, CutArrow cutArrow1, CutArrow cutArrow2)
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
    MnvH1D* temp;

    // ------------------------------------------------------------------------
    // Get All Background
    // ------------------------------------------------------------------------
    var = Form("%s_%d",var_name.c_str(),2);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Background");

    // Get Stats
    double nBckg = temp->GetEntries(); 
    max_bin = temp->GetMaximumBin();
    hist_max = hist_max + temp->GetBinContent(max_bin);
    bin_width = temp->GetBinWidth(1);

    // Add to TObjArray
    mc_hists->Add(temp);

    // ------------------------------------------------------------------------
    // Get Signal
    // ------------------------------------------------------------------------
    var = Form("%s_%d",var_name.c_str(),1);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Signal");

    // Get Stats
    double nSignal = temp->GetEntries(); 
    max_bin = temp->GetMaximumBin();
    hist_max = hist_max + temp->GetBinContent(max_bin);

    // Add to TObjArray
    mc_hists->Add(temp);

    // ------------------------------------------------------------------------
    // Plot  - If you want A Log Plot, axis_minimum = 0.1
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas(var_name.c_str(),var_name.c_str(),1280,800);
    ApplyStyle(plotter);
    //plotter->axis_minimum = 0.1;
    plotter->axis_minimum = 0.0;
    plotter->DrawStackedMC(mc_hists,1,"TR");

    // Add Plot Labels
    var = Form("%s_%d",var_name.c_str(),0);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    plotter->AddHistoTitle(temp->GetTitle());
    const double y_pos = 0.88;
    const double y_diff = 0.033;
    plotter->AddPlotLabel(Form("nSignal = %3.0f, nBckg = %3.0f",nSignal,nBckg),0.3,y_pos,y_diff,kBlue); 

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
        deltaMass.DrawLine(1232,0,1232,hist_max);
    }


    // Print Plot
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),"_mc_bckg_all.png"), "png");

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
    c->Print(Form("%s%s%s%s",plotDir.c_str(),var_name.c_str(),"_mc_bckg_Type_",".png"), "png");

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

    c->Print(Form("%s%s%d%s%s%s",plotDir.c_str(),"gamma_",gammaID,"_mc_",var_name.c_str(),"_by_particle.png"), "png");

    delete legend;
    delete hs;
    delete c;
}

void CCProtonPi0_Plotter::AddNormBox(MnvPlotter* plotter, bool isPOTNorm, double mc_ratio)
{
    const double y_pos = 0.88;
    const double y_diff = 0.033;
    //plotter->AddPlotLabel("Playlist: minerva1",0.3,y_pos,y_diff,kBlue);
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
    TCanvas* c = new TCanvas("c","c",1280,800);

    // Plot Options
    hist1D->GetYaxis()->SetTitle("Reconstruction Efficiency");
    hist1D->SetMinimum(0.0);
    hist1D->SetMaximum(0.15);
    hist1D->SetLineColor(kRed);
    hist1D->SetLineWidth(3);
    hist1D->SetFillColor(kWhite);

    hist1D->Draw("HIST");
    gPad->Update();
    gStyle->SetOptStat(111111); 

    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".png"), "png");
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

    // Get Bckg: SinglePiPlus
    var = Form("%s_%d",var_name.c_str(),5);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: 1#pi^{+}");
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
    MnvH1D* temp;

    // ------------------------------------------------------------------------
    // Get Signal
    // ------------------------------------------------------------------------
    var = Form("%s_%d",var_name.c_str(),1);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Signal");
    temp->SetLineColor(kGreen);
    temp->SetFillColor(kGreen);

    // Get Stat1s
    double nSignal = temp->GetEntries(); 
    max_bin = temp->GetMaximumBin();
    hist_max = hist_max + temp->GetBinContent(max_bin);
    bin_width = temp->GetBinWidth(1);

    // Add to TObjArray
    mc_hists->Add(temp);

    // ------------------------------------------------------------------------
    // Plot  - If you want A Log Plot, axis_minimum = 0.1
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas(var_name.c_str(),var_name.c_str(),1280,800);
    ApplyStyle(plotter);
    //plotter->axis_minimum = 0.1;
    plotter->axis_minimum = 0.0;
    plotter->DrawStackedMC(mc_hists,1,"TR",0);

    // Add Plot Labels
    var = Form("%s_%d",var_name.c_str(),0);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    plotter->AddHistoTitle(temp->GetTitle());
    const double y_pos = 0.88;
    const double y_diff = 0.033;
    plotter->AddPlotLabel(Form("nSignal = %3.0f",nSignal),0.3,y_pos,y_diff,kBlue); 

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



    // Print Plot
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),"_mc_signal.png"), "png");

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
    c->SetWindowSize(w,h);

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    ApplyStyle(plotter);
    plotter->DrawNormalizedMigrationHistogram(hist2D);

    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),"_migration.png"), "png");

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
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".png"), "png");

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
    c->Print(Form("%s%s%s%s",plotDir.c_str(),"Background_Subtraction",tag.c_str(),".png"), "png");

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

#endif


