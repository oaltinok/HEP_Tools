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
    //DrawStackedMC_BckgWithPi0(dir, var_name, plotDir, nCutArrows, cutArrow1, cutArrow2);
    DrawStackedMC_BckgType(dir, var_name, plotDir, nCutArrows, cutArrow1, cutArrow2);
}

void CCProtonPi0_Plotter::DrawDataStackedMC(rootDir &dir, std::string var_name, std::string plotDir, int nCutArrows, CutArrow cutArrow1, CutArrow cutArrow2)
{
    // ------------------------------------------------------------------------
    // POT Normalized Plots 
    // ------------------------------------------------------------------------
    DrawDataStackedMC_BckgAll(dir, var_name,plotDir, true, nCutArrows, cutArrow1, cutArrow2);
    DrawDataStackedMC_BckgWithPi0(dir, var_name,plotDir, false, nCutArrows,  cutArrow1, cutArrow2);
    DrawDataStackedMC_BckgType(dir, var_name,plotDir, false, nCutArrows, cutArrow1, cutArrow2);

    // ------------------------------------------------------------------------
    // Area Normalized Plots 
    // ------------------------------------------------------------------------
    DrawDataStackedMC_BckgAll(dir, var_name,plotDir, false, nCutArrows, cutArrow1, cutArrow2);
    DrawDataStackedMC_BckgWithPi0(dir, var_name,plotDir, false, nCutArrows,  cutArrow1, cutArrow2);
    DrawDataStackedMC_BckgType(dir, var_name,plotDir, false, nCutArrows, cutArrow1, cutArrow2);

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

    // ------------------------------------------------------------------------
    // Get Data Histogram
    // ------------------------------------------------------------------------
    var = Form("%s_%d",var_name.c_str(),0);
    MnvH1D* data = (MnvH1D*)f_data->Get(var.c_str()); 
    max_bin = data->GetMaximumBin();
    hist_max = (hist_max + data->GetBinContent(max_bin))*1.2;
    bin_width = data->GetBinWidth(1);

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

    // Get Signal
    var = Form("%s_%d",var_name.c_str(),1);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Signal");
    mc_hists->Add(temp);

    // ------------------------------------------------------------------------
    // MC Normalization 
    // ------------------------------------------------------------------------
    var = Form("%s_%d",var_name.c_str(),0);
    MnvH1D* mc_all = (MnvH1D*)f_mc->Get(var.c_str()); 
    std::string norm_label; 
    double mc_ratio = GetMCNormalization(norm_label, isPOTNorm, data, mc_all);

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
        deltaMass.SetLineWidth(2);
        deltaMass.SetLineColor(kBlue);
        deltaMass.DrawLine(1.232,0,1.232,hist_max);
    }


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
    double bin_width;

    // ------------------------------------------------------------------------
    // Get Data Histogram
    // ------------------------------------------------------------------------
    var = Form("%s_%d",var_name.c_str(),0);
    MnvH1D* data = (MnvH1D*)f_data->Get(var.c_str()); 

    // ------------------------------------------------------------------------
    // Fill TObjArray - For MC Histograms
    // ------------------------------------------------------------------------
    TObjArray* mc_hists = new TObjArray;
    FormTObjArray_BckgType(f_mc, var_name, mc_hists, hist_max, bin_width);

    // ------------------------------------------------------------------------
    // MC Normalization 
    // ------------------------------------------------------------------------
    var = Form("%s_%d",var_name.c_str(),0);
    MnvH1D* mc_all = (MnvH1D*)f_mc->Get(var.c_str()); 
    std::string norm_label; 
    double mc_ratio = GetMCNormalization(norm_label, isPOTNorm, data, mc_all);

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

    // If Cut Histogram - Add Cut Arrows
    if (nCutArrows == 1){
        AddCutArrow(plotter, cutArrow1, hist_max, bin_width);
    }else if (nCutArrows == 2){
        AddCutArrow(plotter, cutArrow1, hist_max, bin_width);
        AddCutArrow(plotter, cutArrow2, hist_max, bin_width);
    }

    // Print Plot
    c->Print(Form("%s%s%s%s%s",plotDir.c_str(),var_name.c_str(),"_bckg_Type_",norm_label.c_str(),".png"), "png");

    delete c;
    delete plotter;
}

void CCProtonPi0_Plotter::DrawDataStackedMC_BckgWithPi0(rootDir &dir, std::string var_name, std::string plotDir, bool isPOTNorm, int nCutArrows, CutArrow cutArrow1, CutArrow cutArrow2)
{
    std::string rootDir_data = dir.data;
    std::string rootDir_mc = dir.mc;

    TFile* f_data = new TFile(rootDir_data.c_str());
    TFile* f_mc = new TFile(rootDir_mc.c_str());

    std::string var;
    double hist_max = 0;
    double bin_width;

    // ------------------------------------------------------------------------
    // Get Data Histogram
    // ------------------------------------------------------------------------
    var = Form("%s_%d",var_name.c_str(),0);
    MnvH1D* data = (MnvH1D*)f_data->Get(var.c_str()); 

    // ------------------------------------------------------------------------
    // Fill TObjArray - For MC Histograms
    // ------------------------------------------------------------------------
    TObjArray* mc_hists = new TObjArray;
    FormTObjArray_BckgWithPi0(f_mc, var_name, mc_hists, hist_max, bin_width);
 
    // ------------------------------------------------------------------------
    // MC Normalization 
    // ------------------------------------------------------------------------
    var = Form("%s_%d",var_name.c_str(),0);
    MnvH1D* mc_all = (MnvH1D*)f_mc->Get(var.c_str()); 
    std::string norm_label; 
    double mc_ratio = GetMCNormalization(norm_label, isPOTNorm, data, mc_all);

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    TCanvas* c = new TCanvas("c","c",1280,800);
    MnvPlotter* plotter = new MnvPlotter();
    ApplyStyle(plotter);
    //plotter->axis_minimum = 0.1;
    plotter->axis_minimum = 0.0;
    plotter->DrawDataStackedMC(data,mc_hists,mc_ratio ,"TR","Data",0);

    // If Cut Histogram - Add Cut Arrows
    if (nCutArrows == 1){
        AddCutArrow(plotter, cutArrow1, hist_max, bin_width);
    }else if (nCutArrows == 2){
        AddCutArrow(plotter, cutArrow1, hist_max, bin_width);
        AddCutArrow(plotter, cutArrow2, hist_max, bin_width);
    }

    // Add Plot Labels
    plotter->AddHistoTitle(data[0].GetTitle());
    AddNormBox(plotter, isPOTNorm, mc_ratio);

    // Print Plot
    c->Print(Form("%s%s%s%s%s",plotDir.c_str(),var_name.c_str(),"_bckg_Pi0_",norm_label.c_str(),".png"), "png");

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

void CCProtonPi0_Plotter::DrawDataMC(rootDir& dir, std::string var_name, std::string plotDir)
{
    // POT Normalized
    DrawDataMC(dir,var_name,plotDir,true);
    // Area Normalized
    DrawDataMC(dir,var_name,plotDir,false);
}

void CCProtonPi0_Plotter::DrawDataMC(rootDir& dir, std::string var_name, std::string plotDir, bool isPOTNorm)
{
    std::string rootDir_mc = dir.mc;
    std::string rootDir_data = dir.data;

    TFile* f_mc = new TFile(rootDir_mc.c_str());
    TFile* f_data = new TFile(rootDir_data.c_str());

    std::string var = Form("%s_%d",var_name.c_str(),0);
    MnvH1D* mc = (MnvH1D*)f_mc->Get(var.c_str());
    MnvH1D* data = (MnvH1D*)f_data->Get(var.c_str()); 

    // ------------------------------------------------------------------------
    // MC Normalization 
    // ------------------------------------------------------------------------
    std::string norm_label; 
    double mc_ratio = GetMCNormalization(norm_label, isPOTNorm, data, mc);

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas("c","c",1280,800);

    plotter->DrawDataMC(data, mc, mc_ratio, "TR", false);

    // Add Plot Labels
    plotter->AddHistoTitle(data[0].GetTitle());
    AddNormBox(plotter, isPOTNorm, mc_ratio);

    // Add Pi0 InvMass Line
    std::size_t found = var_name.find("invMass");
    if (found != std::string::npos){
        TLine pi0Mass;
        pi0Mass.SetLineWidth(2);
        pi0Mass.SetLineColor(kBlue);
        pi0Mass.DrawLine(134.98,0,134.98,3200);
    }

    // Print Plot
    c->Print(Form("%s%s%s%s%s",plotDir.c_str(),var_name.c_str(),"_",norm_label.c_str(),".png"), "png");

    delete c;
    delete plotter;

    if (isPOTNorm){
        DrawDataMCRatio(dir,var_name,plotDir,isPOTNorm);
    }
}

void CCProtonPi0_Plotter::DrawDataMCRatio(rootDir& dir, std::string var_name, std::string plotDir, bool isPOTNorm)
{
    std::string rootDir_mc = dir.mc;
    std::string rootDir_data = dir.data;

    TFile* f_mc = new TFile(rootDir_mc.c_str());
    TFile* f_data = new TFile(rootDir_data.c_str());

    // ------------------------------------------------------------------------
    // Get data and MC 
    // ------------------------------------------------------------------------
    std::string var = Form("%s_%d",var_name.c_str(),0);
    MnvH1D* mc = (MnvH1D*)f_mc->Get(var.c_str());
    MnvH1D* data = (MnvH1D*)f_data->Get(var.c_str()); 

    // ------------------------------------------------------------------------
    // MC Normalization 
    // ------------------------------------------------------------------------
    std::string norm_label; 
    double mc_ratio = GetMCNormalization(norm_label, isPOTNorm, data, mc);

    // ------------------------------------------------------------------------
    // Plot
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas("c","c",1280,800);

    plotter->DrawDataMCRatio(data, mc, mc_ratio);

    // Add Plot Labels
    plotter->AddHistoTitle(data[0].GetTitle());
    AddNormBox(plotter, isPOTNorm, mc_ratio);

    const double y_pos = 0.88;
    const double y_diff = 0.033;
    plotter->AddPlotLabel("Playlist: minerva1",0.3,y_pos,y_diff,kBlue);
    plotter->AddPOTNormBox(data_POT,mc_POT,0.3,y_pos-y_diff);

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

void CCProtonPi0_Plotter::Draw1DHist(rootDir& dir, std::string var_name, std::string plotDir, bool isLogScale)
{
    std::string root_dir = dir.mc;

    // Get Histogram
    TFile* f = new TFile(root_dir.c_str());
    TH1D* hist1D = (TH1D*)f->Get(var_name.c_str());

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
    gStyle->SetOptStat(111111); 

    // Find Peak
    int max_bin = hist1D->GetMaximumBin();
    double max_bin_location = hist1D->GetBinCenter(max_bin);
    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.03);
    text.DrawLatex(0.15,0.85,Form("%s%3.2f", "Peak at ",max_bin_location));

    // Error Ranges
    //double max_bin_value = hist1D->GetBinContent(max_bin);
    //double error = 0.33;
    //TLine err;
    //err.SetLineWidth(2);
    //err.SetLineColor(kBlack);
    //err.DrawLine(error,0,error,max_bin_value);
    //err.DrawLine(-error,0,-error,max_bin_value);


    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".png"), "png");
    delete c;
    delete hist1D;
    delete f;
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

    double line_min = hist2D->GetXaxis()->GetBinLowEdge(1);
    double line_max = hist2D->GetXaxis()->GetBinLowEdge(nBinsX);
    TLine xy;
    xy.SetLineWidth(2);
    xy.SetLineColor(kBlack);
    xy.DrawLine(line_min,line_min,line_max,line_max);

    //TLine fit;
    //fit.SetLineWidth(2);
    //fit.SetLineColor(kRed);
    //fit.DrawLine(0,80.9,500,585.9);

    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".png"), "png");

    delete p;
    delete c;
    delete f;
}

void CCProtonPi0_Plotter::AddCutArrow(MnvPlotter* plotter, CutArrow &cutArrow, double hist_max, double arrow_length)
{
    double cut_location = cutArrow.cut_location;
    double ymin = 0;
    double ymax = hist_max + hist_max*0.1;
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

void CCProtonPi0_Plotter::DrawStackedMC_BckgWithPi0(rootDir &dir, std::string var_name, std::string plotDir, int nCutArrows, CutArrow cutArrow1, CutArrow cutArrow2)
{
    std::string rootDir_mc = dir.mc;

    TFile* f_mc = new TFile(rootDir_mc.c_str());

    double hist_max = 0;
    double bin_width;

    // ------------------------------------------------------------------------
    // Fill TObjArray - For MC Histograms
    // ------------------------------------------------------------------------
    TObjArray* mc_hists = new TObjArray;
    FormTObjArray_BckgWithPi0(f_mc, var_name, mc_hists, hist_max, bin_width);
    
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    TCanvas* c = new TCanvas("c","c",1280,800);
    MnvPlotter* plotter = new MnvPlotter();
    ApplyStyle(plotter);
    //plotter->axis_minimum = 0.1;
    plotter->axis_minimum = 0.0;
    plotter->DrawStackedMC(mc_hists,1,"TR",0);

    // Add Plot Labels
    std::string var = Form("%s_%d",var_name.c_str(),0);
    MnvH1D* temp = (MnvH1D*)f_mc->Get(var.c_str());
    plotter->AddHistoTitle(temp->GetTitle());

    // Add Pi0 InvMass Line
    std::size_t found = var_name.find("invMass");
    if (found != std::string::npos){
        TLine pi0Mass;
        pi0Mass.SetLineWidth(2);
        pi0Mass.SetLineColor(kBlue);
        pi0Mass.DrawLine(134.98,0,134.98,hist_max);
    }

    // If Cut Histogram - Add Cut Arrows
    if (nCutArrows == 1){
        AddCutArrow(plotter, cutArrow1, hist_max, bin_width);
    }else if (nCutArrows == 2){
        AddCutArrow(plotter, cutArrow1, hist_max, bin_width);
        AddCutArrow(plotter, cutArrow2, hist_max, bin_width);
    }

    // Print Plot
    c->Print(Form("%s%s%s%s",plotDir.c_str(),var_name.c_str(),"_mc_bckg_Pi0_",".png"), "png");

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
        mc_ratio = POT_Ratio_data_mc;
        norm_label = "POT";
    }else{
        mc_ratio = mc->GetAreaNormFactor(data);
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

void CCProtonPi0_Plotter::DrawEfficiencyCurve(std::string var_name, std::string plotDir, TH1D* all_signal, TH1D* signal)
{
    TH1D* h_eff = new TH1D();
    h_eff = signal;
    h_eff->Divide(all_signal);
    h_eff->GetYaxis()->SetTitle("Efficiency");
    h_eff->SetMaximum(0.15);

    // Create Canvas
    TCanvas* c = new TCanvas("c","c",1280,800);

    // Plot Options
    h_eff->SetLineColor(kRed);
    h_eff->SetLineWidth(3);
    h_eff->SetFillColor(kWhite);

    h_eff->Draw();
    gPad->Update();
    gStyle->SetOptStat(111111); 

    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".png"), "png");
    delete c;
    delete h_eff;
}


void CCProtonPi0_Plotter::FormTObjArray_BckgType(TFile* f_mc, std::string var_name, TObjArray* mc_hists, double &hist_max, double &bin_width) 
{
    std::string var;
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

    // Get Bckg: QELike
    var = Form("%s_%d",var_name.c_str(),8);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: QELike");
    max_bin = temp->GetMaximumBin();
    hist_max = hist_max + temp->GetBinContent(max_bin);
    temp->SetLineColor(kOrange);
    temp->SetFillColor(kOrange);
    mc_hists->Add(temp);

    // Get Bckg: SingleChargedPion
    var = Form("%s_%d",var_name.c_str(),9);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: SingleChargedPion");
    max_bin = temp->GetMaximumBin();
    hist_max = hist_max + temp->GetBinContent(max_bin);
    temp->SetLineColor(kRed);
    temp->SetFillColor(kRed);
    mc_hists->Add(temp);

    // Get Bckg: SingleChargedPion Charge Exchanged
    var = Form("%s_%d",var_name.c_str(),10);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: SingleChargedPion_ChargeExc");
    max_bin = temp->GetMaximumBin();
    hist_max = hist_max + temp->GetBinContent(max_bin);
    temp->SetLineColor(kTeal);
    temp->SetFillColor(kTeal);
    mc_hists->Add(temp);

    // Get Bckg: DoublePionWithPi0
    var = Form("%s_%d",var_name.c_str(),11);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: DoublePionWithPi0");
    max_bin = temp->GetMaximumBin();
    hist_max = hist_max + temp->GetBinContent(max_bin);
    temp->SetLineColor(kMagenta);
    temp->SetFillColor(kMagenta);
    mc_hists->Add(temp);

    // Get Bckg: DoublePionWithoutPi0
    var = Form("%s_%d",var_name.c_str(),12);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: DoublePionWithoutPi0");
    max_bin = temp->GetMaximumBin();
    hist_max = hist_max + temp->GetBinContent(max_bin);
    temp->SetLineColor(kYellow);
    temp->SetFillColor(kYellow);
    mc_hists->Add(temp);

    // Get Bckg: MultiPionWithPi0
    var = Form("%s_%d",var_name.c_str(),13);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: MultiPionWithPi0");
    max_bin = temp->GetMaximumBin();
    hist_max = hist_max + temp->GetBinContent(max_bin);
    temp->SetLineColor(kBlue);
    temp->SetFillColor(kBlue);
    mc_hists->Add(temp);

    // Get Bckg: MultiPionWithoutPi0
    var = Form("%s_%d",var_name.c_str(),14);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: MultiPionWithoutPi0");
    max_bin = temp->GetMaximumBin();
    hist_max = hist_max + temp->GetBinContent(max_bin);
    temp->SetLineColor(kBlack);
    temp->SetFillColor(kBlack);
    mc_hists->Add(temp);

    // Get Bckg: NC
    var = Form("%s_%d",var_name.c_str(),6);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: NC");
    max_bin = temp->GetMaximumBin();
    hist_max = hist_max + temp->GetBinContent(max_bin);
    temp->SetLineColor(kGray);
    temp->SetFillColor(kGray);
    mc_hists->Add(temp);

    // Get Bckg: AntiNeutrino
    var = Form("%s_%d",var_name.c_str(),7);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: AntiNeutrino");
    max_bin = temp->GetMaximumBin();
    hist_max = hist_max + temp->GetBinContent(max_bin);
    temp->SetLineColor(kGray);
    temp->SetFillColor(kGray);
    mc_hists->Add(temp);

    // Get Bckg: Other
    var = Form("%s_%d",var_name.c_str(),15);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: Other");
    max_bin = temp->GetMaximumBin();
    hist_max = hist_max + temp->GetBinContent(max_bin);
    temp->SetLineColor(kGray);
    temp->SetFillColor(kGray);
    mc_hists->Add(temp);

}

void CCProtonPi0_Plotter::FormTObjArray_BckgWithPi0(TFile* f_mc, std::string var_name, TObjArray* mc_hists, double &hist_max, double &bin_width) 
{
    std::string var;
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

    // Get Bckg: NoPi0
    var = Form("%s_%d",var_name.c_str(),3);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: NoPi0");
    max_bin = temp->GetMaximumBin();
    hist_max = hist_max + temp->GetBinContent(max_bin);
    temp->SetLineColor(kYellow);
    temp->SetFillColor(kYellow);
    mc_hists->Add(temp);

    // Get Bckg: SinglePi0
    var = Form("%s_%d",var_name.c_str(),4);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: SinglePi0");
    max_bin = temp->GetMaximumBin();
    hist_max = hist_max + temp->GetBinContent(max_bin);
    temp->SetLineColor(kRed);
    temp->SetFillColor(kRed);
    mc_hists->Add(temp);

    // Get Bckg: MultiPi0
    var = Form("%s_%d",var_name.c_str(),5);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: MultiPi0");
    max_bin = temp->GetMaximumBin();
    hist_max = hist_max + temp->GetBinContent(max_bin);
    temp->SetLineColor(kBlue);
    temp->SetFillColor(kBlue);
    mc_hists->Add(temp);
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

    // Get Stats
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
        deltaMass.DrawLine(1232,0,1232,hist_max);
    }



    // Print Plot
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),"_mc_signal.png"), "png");

    delete c;
    delete plotter;
}


#endif


