#ifndef CCProtonPi0_Plotter_W_Fit_cpp
#define CCProtonPi0_Plotter_W_Fit_cpp

#include "CCProtonPi0_Plotter.h"

using namespace PlotUtils;

void CCProtonPi0_Plotter::init_W_FitResults()
{
    pars_MC_deltaRES[0] = 31.47;
    pars_MC_deltaRES[1] = 1.215;
    pars_MC_deltaRES[2] = 0.246;

    pars_MC_otherRES[0] = 68.48;
    pars_MC_otherRES[1] = 1.473;
    pars_MC_otherRES[2] = 0.2864;

    pars_MC_nonRES_G1[0] = 28.66;
    pars_MC_nonRES_G1[1] = 1.197;
    pars_MC_nonRES_G1[2] = 0.1498;

    pars_MC_nonRES_G2[0] = 55.57;
    pars_MC_nonRES_G2[1] = 1.439;
    pars_MC_nonRES_G2[2] = 0.3686;
}

void CCProtonPi0_Plotter::plot_W_FitMinuit(double wgt_DeltaRES, double wgt_OtherRES, double wgt_NonRES)
{
    // Get Background Subtracted Data
    TFile* f_data = new TFile(rootDir_Interaction.data.c_str());
    TFile* f_mc = new TFile(rootDir_Interaction.mc.c_str());
    MnvH1D* data = GetMnvH1D(f_data, "W_All_0");
    MnvH1D* mc_bckg = GetMnvH1D(f_mc, "W_All_2");

    NormalizeHistogram(mc_bckg);
    
    const double nBckg = 2997.1; // Got number from Cross Section Calculation Log File
    mc_bckg->Scale(nBckg);

    data->Add(mc_bckg, -1); 

    // Get MC Signal Types
    MnvH1D* mc_DeltaRES = GetMnvH1D(f_mc, "W_All_7");
    mc_DeltaRES->Scale(POT_ratio);
    mc_DeltaRES->Scale(wgt_DeltaRES);
    mc_DeltaRES->SetFillColor(kGray+1);
    mc_DeltaRES->SetLineColor(kGray+1);
    mc_DeltaRES->SetFillStyle(1001);

    MnvH1D* mc_OtherRES = GetMnvH1D(f_mc, "W_All_8");
    mc_OtherRES->Scale(POT_ratio);
    mc_OtherRES->Scale(wgt_OtherRES);
    mc_OtherRES->SetFillColor(kGreen+2);
    mc_OtherRES->SetLineColor(kGreen+2);
    mc_OtherRES->SetFillStyle(1001);

    MnvH1D* mc_NonRES = GetMnvH1D(f_mc, "W_All_9");
    mc_NonRES->Scale(POT_ratio);
    mc_NonRES->Scale(wgt_NonRES);
    mc_NonRES->SetFillColor(kRed+1);
    mc_NonRES->SetLineColor(kRed+1);
    mc_NonRES->SetFillStyle(1001);

    TCanvas* c = new TCanvas("c","c",1280,800);
    THStack *hs = new THStack("hs","W Minuit Fit");
    hs->Add(mc_DeltaRES);
    hs->Add(mc_OtherRES);
    hs->Add(mc_NonRES);

    double hs_max = hs->GetMaximum();
    hs->SetMaximum(hs_max * 1.75);
    hs->Draw("HIST");
    hs->GetXaxis()->SetTitle(data->GetXaxis()->GetTitle());
    hs->GetYaxis()->SetTitle(data->GetYaxis()->GetTitle());

    data->SetMarkerStyle(20);
    data->SetMarkerSize(1);
    data->SetMarkerColor(kBlack);
    data->SetLineWidth(2);
    data->SetLineColor(kBlack);
    data->Draw("SAME E1 X0");

    // ------------------------------------------------------------------------
    // Plot Labels 
    // ------------------------------------------------------------------------
    TLegend *legend = new TLegend(0.6,0.7,0.8,0.9);  
    legend->AddEntry(data, "Data (3.33e20 POT)", "lep" );
    legend->AddEntry(mc_DeltaRES, "#Delta(1232) resonance","f");
    legend->AddEntry(mc_OtherRES, "Other resonances","f");
    legend->AddEntry(mc_NonRES, "Non-Resonant","f");
    legend->Draw();

    // Add Text
    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.03);

    text.DrawLatex(0.7, 0.64, Form("wgt(DeltaRES) = %3.2f", wgt_DeltaRES));
    text.DrawLatex(0.7, 0.60, Form("wgt(OtherRES) = %3.2f", wgt_OtherRES));
    text.DrawLatex(0.7, 0.56, Form("wgt(NonRES) = %3.2f", wgt_NonRES));

    // Plot Output
    gStyle->SetEndErrorSize(6);
    gStyle->SetOptStat(0); 
    c->Update();
    std::string plotDir = Folder_List::plotDir_OtherStudies;
    std::string out_name;
    std::string type;
    if (wgt_DeltaRES == 1.0 && wgt_OtherRES == 1.0) type = "Nominal";
    else type = "Fitted";
    out_name = plotDir + "W_Fit_Minuit_" + type + ".png"; 

    c->Print(out_name.c_str(),"png");

//    delete mc_DeltaRES;
//    delete mc_OtherRES;
//    delete mc_NonRES;
//    delete mc_bckg;
//    delete data;
//    delete f_data;
//    delete f_mc;
//    delete legend;
}

void CCProtonPi0_Plotter::plot_W_FitResults()
{    
    // Fit Results from MATLAB
    // Fixed Non RES
    double pars_fixed_non_deltaRES[3] = {39.23, 1.167, 0.3365};
    double pars_fixed_non_otherRES[3] = {50, 1.311, 0.248};
    double pars_fixed_non_nonRES_G1[3] = {pars_MC_nonRES_G1[0], pars_MC_nonRES_G1[1], pars_MC_nonRES_G1[2]};
    double pars_fixed_non_nonRES_G2[3] = {pars_MC_nonRES_G2[0], pars_MC_nonRES_G2[1], pars_MC_nonRES_G2[2]};
    W_Fit_Data("FixedNonRES", pars_fixed_non_deltaRES, pars_fixed_non_otherRES, pars_fixed_non_nonRES_G1, pars_fixed_non_nonRES_G2, 6);

    W_Fit_MC();
    W_Fit_MC_DeltaRES();
    W_Fit_MC_OtherRES();
    W_Fit_MC_NonRES();
    //W_Fit_MC_Errors();
}

void CCProtonPi0_Plotter::printBins_W()
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    // Print Data Bins for Fit
    TFile* f_data = new TFile(rootDir_Interaction.data.c_str());
    TFile* f_mc = new TFile(rootDir_Interaction.mc.c_str());
    MnvH1D* h_data_W_All = GetMnvH1D(f_data, "W_All_0");
    MnvH1D* h_mc_bckg_W_All = GetMnvH1D(f_mc, "W_All_2");

    NormalizeHistogram(h_mc_bckg_W_All);
    
    const double nBckg = 2997.1; // Got number from Cross Section Calculation Log File
    h_mc_bckg_W_All->Scale(nBckg);

    h_data_W_All->Add(h_mc_bckg_W_All, -1); 
    printBins(h_data_W_All,"Data: W All");
    
    delete h_mc_bckg_W_All;
    delete h_data_W_All;
    delete f_data;

    // Print MC Bins for Fit
    MnvH1D* h_mc_W_All = GetMnvH1D(f_mc, "W_All_1");
    h_mc_W_All->Scale(POT_ratio);
    printBins(h_mc_W_All,"MC: W All");
  
    MnvH1D* h_mc_W_All_7 = GetMnvH1D(f_mc, "W_All_7");
    h_mc_W_All_7->Scale(POT_ratio);
    printBins(h_mc_W_All_7,"MC: W Delta RES");

    MnvH1D* h_mc_W_All_8 = GetMnvH1D(f_mc, "W_All_8");
    h_mc_W_All_8->Scale(POT_ratio);
    printBins(h_mc_W_All_8,"MC: W Other RES");
   
    MnvH1D* h_mc_W_All_9 = GetMnvH1D(f_mc, "W_All_9");
    h_mc_W_All_9->Scale(POT_ratio);
    printBins(h_mc_W_All_9,"MC: W Non RES");

    delete h_mc_W_All;
    delete h_mc_W_All_8;
    delete h_mc_W_All_9;
    
    //// Delta Res
    //MnvH1D* h_mc_W_All_delta_res = GetMnvH1D(f_mc, "W_All_7");
    //MnvH1D* h_mc_W_1_delta_res = GetMnvH1D(f_mc, "W_1_7");
    //MnvH1D* h_mc_W_2_delta_res = GetMnvH1D(f_mc, "W_2_7");
    //h_mc_W_All_delta_res->Scale(POT_ratio);
    //h_mc_W_1_delta_res->Scale(POT_ratio);
    //h_mc_W_2_delta_res->Scale(POT_ratio);
    //printBins(h_mc_W_All_delta_res,"Delta RES: W All");
    //printBins(h_mc_W_1_delta_res,"Delta RES: W 1 Track");
    //printBins(h_mc_W_2_delta_res,"Delta RES: W 2 Track");
    //delete h_mc_W_All_delta_res;
    //delete h_mc_W_1_delta_res;
    //delete h_mc_W_2_delta_res;
    delete f_mc;
}

void CCProtonPi0_Plotter::NormalizeHistogram(MnvH1D* h)
{
    std::cout<<"\tNormalizing Background Shape on MnvH1D"<<std::endl;
    int NBins = h->GetNbinsX();
    double area = h->Integral();
    double nOverFlow = h->GetBinContent(NBins+1);
    double nUnderFlow = h->GetBinContent(0);
    std::cout<<"\t\tBefore Norm = "<<area<<std::endl;
    h->Scale(1/(area+nOverFlow+nUnderFlow),"",false); // Scale only on CentralValue
    std::cout<<"\t\tAfter Norm = "<<h->Integral()<<std::endl;
    std::cout<<"\tDone!"<<std::endl;
}

void CCProtonPi0_Plotter::W_Fit_Data(std::string fit_name, double* pars_deltaRES, double* pars_otherRES, double* pars_nonRES_G1, double* pars_nonRES_G2, int nPars)
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    // ------------------------------------------------------------------------
    // First Plot Data vs MC using MnvPlotter
    // ------------------------------------------------------------------------
    TFile* f_data = new TFile(rootDir_Interaction.data.c_str());
    TFile* f_mc = new TFile(rootDir_Interaction.mc.c_str());

    // Get Data Background Subtracted
    MnvH1D* h_data = GetMnvH1D(f_data, "W_All_0");
    MnvH1D* h_mc_bckg = GetMnvH1D(f_mc, "W_All_2");

    NormalizeHistogram(h_mc_bckg);
    
    const double nBckg = 2997.1; // Got number from Cross Section Calculation Log File
    h_mc_bckg->Scale(nBckg);
    h_data->Add(h_mc_bckg, -1); 
    
    // Get MC Signal
    MnvH1D* h_mc_signal = GetMnvH1D(f_mc, "W_All_1");
   
    // Plot MC Signal and Bckg Subtracted Data
    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas("c","c",1280,800);

    plotter->headroom = 1.75;
    plotter->data_line_width = 2;
    plotter->data_marker_size = 1.5;
    gStyle->SetEndErrorSize(6);
    plotter->DrawDataMCWithErrorBand(h_data, h_mc_signal, POT_ratio, "N", false, NULL, NULL, false, true);

    // ------------------------------------------------------------------------
    // Second Add Fit Results
    // ------------------------------------------------------------------------
    // Arrays for Fit Results 
    int nPoints = 150;
    double fit_x[nPoints];
    double fit_deltaRES_y[nPoints];
    double fit_otherRES_y[nPoints];
    double fit_nonRES_y[nPoints];
    double fit_nonRES_G1_y[nPoints];
    double fit_nonRES_G2_y[nPoints];
    double fit_total_y[nPoints];

    // X Axis (W)
    double xStart = 0.5;
    double inc = 0.01;
    for (int i = 0; i < nPoints; ++i){
        fit_x[i] = xStart + i*inc; 
    }
  
    Get_BreitWigner(pars_deltaRES, fit_x, fit_deltaRES_y, nPoints);
    Get_Gaussian(pars_otherRES, fit_x, fit_otherRES_y, nPoints);
    Get_Gaussian(pars_nonRES_G1, fit_x, fit_nonRES_G1_y, nPoints);
    Get_Gaussian(pars_nonRES_G2, fit_x, fit_nonRES_G2_y, nPoints);
   
    for (int i = 0; i < nPoints; ++i){
        fit_nonRES_y[i] = fit_nonRES_G1_y[i] + fit_nonRES_G2_y[i];
        fit_total_y[i] = fit_deltaRES_y[i] + fit_otherRES_y[i] + fit_nonRES_y[i];
    }

    TGraph* fit_deltaRES = new TGraph(nPoints, fit_x, fit_deltaRES_y);
    TGraph* fit_otherRES = new TGraph(nPoints, fit_x, fit_otherRES_y);
    TGraph* fit_nonRES = new TGraph(nPoints, fit_x, fit_nonRES_y);
    TGraph* fit_total = new TGraph(nPoints, fit_x, fit_total_y);
    
    fit_deltaRES->SetLineColor(6); // Magenta
    fit_deltaRES->SetLineWidth(4);
    fit_deltaRES->Draw("SAMEL");

    fit_otherRES->SetLineColor(4); // Blue 
    fit_otherRES->SetLineWidth(4);
    fit_otherRES->Draw("SAMEL");

    fit_nonRES->SetLineColor(8); // Green
    fit_nonRES->SetLineWidth(4);
    fit_nonRES->Draw("SAMEL");

    fit_total->SetLineColor(1); // Black
    fit_total->SetLineWidth(4);
    fit_total->Draw("SAMEL");

    // Add Legend
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(1);
    h_data->SetMarkerColor(kBlack);
    h_data->SetLineWidth(2);
    h_data->SetLineColor(kBlack);

    h_mc_signal->SetLineWidth(3);
    h_mc_signal->SetLineColor(kRed);
    h_mc_signal->SetFillColor(kWhite);

    TLegend *legend = new TLegend(0.65,0.65,0.9,0.9);  
    legend->AddEntry(h_data, "Data (3.33e20 POT)", "lep");
    legend->AddEntry(h_mc_signal, "Simulation", "l");
    legend->AddEntry(fit_total, "Fit to Data", "l");
    legend->AddEntry(fit_deltaRES, "Fit: #Delta(1232) RES", "l");
    legend->AddEntry(fit_otherRES, "Fit: Other RES", "l");
    legend->AddEntry(fit_nonRES, "Fit: Non-RES", "l");
    legend->Draw();
 
    // ------------------------------------------------------------------------
    // Calc ChiSq
    // ------------------------------------------------------------------------
    // Arrays for Fit Results 
    nPoints = 20; // Last 4 bins are zero ( W > 1.8 GeV)
    double ChiSq_fit_x[nPoints];
    double ChiSq_data_y[nPoints];
    double ChiSq_fit_deltaRES_y[nPoints];
    double ChiSq_fit_otherRES_y[nPoints];
    double ChiSq_fit_nonRES_y[nPoints];
    double ChiSq_fit_nonRES_G1_y[nPoints];
    double ChiSq_fit_nonRES_G2_y[nPoints];
    double ChiSq_fit_total_y[nPoints];

    // X Axis (W)
    for (int i = 7; i <= nPoints+6; ++i){
        ChiSq_fit_x[i-7] = h_data->GetBinCenter(i); 
        ChiSq_data_y[i-7] = h_data->GetBinContent(i); 
    }

    Get_BreitWigner(pars_deltaRES, ChiSq_fit_x, ChiSq_fit_deltaRES_y, nPoints);
    Get_Gaussian(pars_otherRES, ChiSq_fit_x, ChiSq_fit_otherRES_y, nPoints);
    Get_Gaussian(pars_nonRES_G1, ChiSq_fit_x, ChiSq_fit_nonRES_G1_y, nPoints);
    Get_Gaussian(pars_nonRES_G2, ChiSq_fit_x, ChiSq_fit_nonRES_G2_y, nPoints);
   
    for (int i = 0; i < nPoints; ++i){
        ChiSq_fit_nonRES_y[i] = ChiSq_fit_nonRES_G1_y[i] + ChiSq_fit_nonRES_G2_y[i];
        ChiSq_fit_total_y[i] = ChiSq_fit_deltaRES_y[i] + ChiSq_fit_otherRES_y[i] + ChiSq_fit_nonRES_y[i];
    }
   
    double ChiSq_dof = Calc_ChiSq_dof(ChiSq_data_y, ChiSq_fit_total_y, nPoints, nPars);
 
    // Add Text
    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.03);
    text.DrawLatex(0.20,0.85,Form("#chi^{2}/dof = %3.2f", ChiSq_dof));
    text.DrawLatex(0.20,0.82,Form("dof = %d-%d", nPoints, nPars));


    // Plot 
    c->Update();
    c->Print(Form("%s%s%s%s",plotDir.c_str(),"W_Fit_Data_", fit_name.c_str(), ".png"), "png");

    delete h_data;
    delete h_mc_signal;
    delete h_mc_bckg;
    delete f_data;
    delete f_mc;
    delete fit_deltaRES;
    delete fit_otherRES;
    delete fit_nonRES;
    delete fit_total;
    delete plotter;
    delete legend;
    delete c;

}

void CCProtonPi0_Plotter::W_Fit_MC_DeltaRES()
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    // ------------------------------------------------------------------------
    // First Plot MC using MnvPlotter
    // ------------------------------------------------------------------------
    TFile* f_mc = new TFile(rootDir_Interaction.mc.c_str());

    // Get MC Signal deltaRES
    MnvH1D* h_mc_signal = GetMnvH1D(f_mc, "W_All_7");
    h_mc_signal->Scale(POT_ratio);
   
    // Plot MC Signal deltaRES 
    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas("c","c",1280,800);

    plotter->headroom = 1.75;
    plotter->DrawMCWithErrorBand(h_mc_signal);

    // ------------------------------------------------------------------------
    // Second Add Fit Results
    // ------------------------------------------------------------------------
    // Arrays for Fit Results 
    int nPoints = 150;
    double fit_x[nPoints];
    double fit_deltaRES_y[nPoints];

    // X Axis (W)
    double xStart = 0.5;
    double inc = 0.01;
    for (int i = 0; i < nPoints; ++i){
        fit_x[i] = xStart + i*inc; 
    }

    Get_BreitWigner(pars_MC_deltaRES, fit_x, fit_deltaRES_y, nPoints);
   
    TGraph* fit_deltaRES = new TGraph(nPoints, fit_x, fit_deltaRES_y);
    
    fit_deltaRES->SetLineColor(6); // Magenta
    fit_deltaRES->SetLineWidth(4);
    fit_deltaRES->Draw("SAMEL");

    h_mc_signal->SetLineWidth(3);
    h_mc_signal->SetLineColor(kRed);
    h_mc_signal->SetFillColor(kWhite);

    TLegend *legend = new TLegend(0.65,0.80,0.9,0.9);  
    legend->AddEntry(h_mc_signal, "Signal: #Delta(1232) RES", "le");
    legend->AddEntry(fit_deltaRES, "Fit: #Delta(1232) RES", "l");
    legend->Draw();

    // ------------------------------------------------------------------------
    // Calc ChiSq
    // ------------------------------------------------------------------------
    // Arrays for Fit Results 
    nPoints = 20; // Last 4 bins are zero ( W > 1.8 GeV)
    double ChiSq_fit_x[nPoints];
    double ChiSq_data_y[nPoints];
    double ChiSq_fit_deltaRES_y[nPoints];

    // X Axis (W)
    for (int i = 7; i <= nPoints+6; ++i){
        ChiSq_fit_x[i-7] =h_mc_signal->GetBinCenter(i); 
        ChiSq_data_y[i-7] = h_mc_signal->GetBinContent(i); 
    }

    Get_BreitWigner(pars_MC_deltaRES, ChiSq_fit_x, ChiSq_fit_deltaRES_y, nPoints);
   
    double ChiSq_dof = Calc_ChiSq_dof(ChiSq_data_y, ChiSq_fit_deltaRES_y, nPoints, 3);
 
    // Add Text
    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.03);
    text.DrawLatex(0.20,0.85,Form("#chi^{2}/dof = %3.2f", ChiSq_dof));
    text.DrawLatex(0.20,0.82,Form("dof = %d-%d", nPoints, 3));

    // Plot 
    c->Update();
    c->Print(Form("%s%s",plotDir.c_str(),"W_Fit_MC_deltaRES.png"), "png");

    delete h_mc_signal;
    delete f_mc;
    delete fit_deltaRES;
    delete plotter;
    delete legend;
    delete c;
}

void CCProtonPi0_Plotter::W_Fit_MC_OtherRES()
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    // ------------------------------------------------------------------------
    // First Plot MC using MnvPlotter
    // ------------------------------------------------------------------------
    TFile* f_mc = new TFile(rootDir_Interaction.mc.c_str());

    // Get MC Signal deltaRES
    MnvH1D* h_mc_signal = GetMnvH1D(f_mc, "W_All_8");
    h_mc_signal->Scale(POT_ratio);
   
    // Plot MC Signal deltaRES 
    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas("c","c",1280,800);

    plotter->headroom = 1.75;
    plotter->DrawMCWithErrorBand(h_mc_signal);

    // ------------------------------------------------------------------------
    // Second Add Fit Results
    // ------------------------------------------------------------------------
    // Arrays for Fit Results 
    int nPoints = 150;
    double fit_x[nPoints];
    double fit_otherRES_y[nPoints];

    // X Axis (W)
    double xStart = 0.5;
    double inc = 0.01;
    for (int i = 0; i < nPoints; ++i){
        fit_x[i] = xStart + i*inc; 
    }

    Get_Gaussian(pars_MC_otherRES, fit_x, fit_otherRES_y, nPoints);
   
    TGraph* fit_otherRES = new TGraph(nPoints, fit_x, fit_otherRES_y);
   
    fit_otherRES->SetLineColor(4); // Blue 
    fit_otherRES->SetLineWidth(4);
    fit_otherRES->Draw("SAMEL");

    h_mc_signal->SetLineWidth(3);
    h_mc_signal->SetLineColor(kRed);
    h_mc_signal->SetFillColor(kWhite);

    TLegend *legend = new TLegend(0.65,0.80,0.9,0.9);  
    legend->AddEntry(h_mc_signal, "Signal: Other RES", "le");
    legend->AddEntry(fit_otherRES, "Fit: Other RES", "l");
    legend->Draw();
    // ------------------------------------------------------------------------
    // Calc ChiSq
    // ------------------------------------------------------------------------
    // Arrays for Fit Results 
    nPoints = 20; // Last 4 bins are zero ( W > 1.8 GeV)
    double ChiSq_fit_x[nPoints];
    double ChiSq_data_y[nPoints];
    double ChiSq_fit_otherRES_y[nPoints];

    // X Axis (W)
    for (int i = 7; i <= nPoints+6; ++i){
        ChiSq_fit_x[i-7] = h_mc_signal->GetBinCenter(i); 
        ChiSq_data_y[i-7] = h_mc_signal->GetBinContent(i); 
    }

    Get_Gaussian(pars_MC_otherRES, ChiSq_fit_x, ChiSq_fit_otherRES_y, nPoints);
   
    double ChiSq_dof = Calc_ChiSq_dof(ChiSq_data_y, ChiSq_fit_otherRES_y, nPoints, 3);
 
    // Add Text
    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.03);
    text.DrawLatex(0.20,0.85,Form("#chi^{2}/dof = %3.2f", ChiSq_dof));
    text.DrawLatex(0.20,0.82,Form("dof = %d-%d", nPoints, 3));


    // Plot 
    c->Update();
    c->Print(Form("%s%s",plotDir.c_str(),"W_Fit_MC_otherRES.png"), "png");

    delete h_mc_signal;
    delete f_mc;
    delete fit_otherRES;
    delete plotter;
    delete legend;
    delete c;
}

void CCProtonPi0_Plotter::W_Fit_MC_NonRES()
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    // ------------------------------------------------------------------------
    // First Plot MC using MnvPlotter
    // ------------------------------------------------------------------------
    TFile* f_mc = new TFile(rootDir_Interaction.mc.c_str());

    // Get MC Signal deltaRES
    MnvH1D* h_mc_signal = GetMnvH1D(f_mc, "W_All_9");
    h_mc_signal->Scale(POT_ratio);
   
    // Plot MC Signal deltaRES 
    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas("c","c",1280,800);

    plotter->headroom = 1.75;
    plotter->DrawMCWithErrorBand(h_mc_signal);

    // ------------------------------------------------------------------------
    // Second Add Fit Results
    // ------------------------------------------------------------------------
    // Arrays for Fit Results 
    int nPoints = 150;
    double fit_x[nPoints];
    double fit_nonRES_y[nPoints];
    double fit_nonRES_G1_y[nPoints];
    double fit_nonRES_G2_y[nPoints];

    // X Axis (W)
    double xStart = 0.5;
    double inc = 0.01;
    for (int i = 0; i < nPoints; ++i){
        fit_x[i] = xStart + i*inc; 
    }

    // Non-RES Fit is a Double Gaussian
    Get_Gaussian(pars_MC_nonRES_G1, fit_x, fit_nonRES_G1_y, nPoints);
    Get_Gaussian(pars_MC_nonRES_G2, fit_x, fit_nonRES_G2_y, nPoints);

    for ( int i = 0; i < nPoints; ++i){
        fit_nonRES_y[i] = fit_nonRES_G1_y[i] + fit_nonRES_G2_y[i];
    }
   
    TGraph* fit_nonRES = new TGraph(nPoints, fit_x, fit_nonRES_y);
   
    fit_nonRES->SetLineColor(8); // Green
    fit_nonRES->SetLineWidth(4);
    fit_nonRES->Draw("SAMEL");

    h_mc_signal->SetLineWidth(3);
    h_mc_signal->SetLineColor(kRed);
    h_mc_signal->SetFillColor(kWhite);

    TLegend *legend = new TLegend(0.65,0.80,0.9,0.9);  
    legend->AddEntry(h_mc_signal, "Signal: Non-RES", "le");
    legend->AddEntry(fit_nonRES, "Fit: Non-RES", "l");
    legend->Draw();
 
    // ------------------------------------------------------------------------
    // Calc ChiSq
    // ------------------------------------------------------------------------
    // Arrays for Fit Results 
    nPoints = 20; // Last 4 bins are zero ( W > 1.8 GeV)
    double ChiSq_fit_x[nPoints];
    double ChiSq_data_y[nPoints];
    double ChiSq_fit_nonRES_y[nPoints];
    double ChiSq_fit_nonRES_G1_y[nPoints];
    double ChiSq_fit_nonRES_G2_y[nPoints];

    // X Axis (W)
    for (int i = 7; i <= nPoints+6; ++i){
        ChiSq_fit_x[i-7] = h_mc_signal->GetBinCenter(i); 
        ChiSq_data_y[i-7] = h_mc_signal->GetBinContent(i); 
    }

    // Non-RES Fit is a Double Gaussian
    Get_Gaussian(pars_MC_nonRES_G1, ChiSq_fit_x, ChiSq_fit_nonRES_G1_y, nPoints);
    Get_Gaussian(pars_MC_nonRES_G2, ChiSq_fit_x, ChiSq_fit_nonRES_G2_y, nPoints);
 
    for (int i = 7; i <= nPoints+6; ++i){
        ChiSq_fit_nonRES_y[i-7] = ChiSq_fit_nonRES_G1_y[i-7] + ChiSq_fit_nonRES_G2_y[i-7];    
    }
  
    double ChiSq_dof = Calc_ChiSq_dof(ChiSq_data_y, ChiSq_fit_nonRES_y, nPoints, 3);
 
    // Add Text
    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.03);
    text.DrawLatex(0.20,0.85,Form("#chi^{2}/dof = %3.2f", ChiSq_dof));
    text.DrawLatex(0.20,0.82,Form("dof = %d-%d", nPoints, 3));

    // Plot 
    c->Update();
    c->Print(Form("%s%s",plotDir.c_str(),"W_Fit_MC_nonRES.png"), "png");

    delete h_mc_signal;
    delete f_mc;
    delete fit_nonRES;
    delete plotter;
    delete legend;
    delete c;
}
void CCProtonPi0_Plotter::W_Fit_MC()
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    // ------------------------------------------------------------------------
    // First Plot MC using MnvPlotter
    // ------------------------------------------------------------------------
    TFile* f_mc = new TFile(rootDir_Interaction.mc.c_str());

    // Get MC Signal
    MnvH1D* h_mc_signal = GetMnvH1D(f_mc, "W_All_1");
    h_mc_signal->Scale(POT_ratio);
   
    // Plot MC Signal 
    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas("c","c",1280,800);

    plotter->headroom = 1.75;
    plotter->DrawMCWithErrorBand(h_mc_signal);

    // ------------------------------------------------------------------------
    // Second Add Fit Results
    // ------------------------------------------------------------------------
    // Arrays for Fit Results 
    int nPoints = 150;
    double fit_x[nPoints];
    double fit_deltaRES_y[nPoints];
    double fit_otherRES_y[nPoints];
    double fit_nonRES_y[nPoints];
    double fit_nonRES_G1_y[nPoints];
    double fit_nonRES_G2_y[nPoints];
    double fit_total_y[nPoints];

    // X Axis (W)
    double xStart = 0.5;
    double inc = 0.01;
    for (int i = 0; i < nPoints; ++i){
        fit_x[i] = xStart + i*inc; 
    }

    Get_BreitWigner(pars_MC_deltaRES, fit_x, fit_deltaRES_y, nPoints);
    Get_Gaussian(pars_MC_otherRES, fit_x, fit_otherRES_y, nPoints);
    Get_Gaussian(pars_MC_nonRES_G1, fit_x, fit_nonRES_G1_y, nPoints);
    Get_Gaussian(pars_MC_nonRES_G2, fit_x, fit_nonRES_G2_y, nPoints);
   
    for (int i = 0; i < nPoints; ++i){
        fit_nonRES_y[i] = fit_nonRES_G1_y[i] + fit_nonRES_G2_y[i];
        fit_total_y[i] = fit_deltaRES_y[i] + fit_otherRES_y[i] + fit_nonRES_y[i];
    }

    TGraph* fit_deltaRES = new TGraph(nPoints, fit_x, fit_deltaRES_y);
    TGraph* fit_otherRES = new TGraph(nPoints, fit_x, fit_otherRES_y);
    TGraph* fit_nonRES = new TGraph(nPoints, fit_x, fit_nonRES_y);
    TGraph* fit_total = new TGraph(nPoints, fit_x, fit_total_y);
    
    fit_deltaRES->SetLineColor(6); // Magenta
    fit_deltaRES->SetLineWidth(4);
    fit_deltaRES->Draw("SAMEL");

    fit_otherRES->SetLineColor(4); // Blue 
    fit_otherRES->SetLineWidth(4);
    fit_otherRES->Draw("SAMEL");

    fit_nonRES->SetLineColor(8); // Green
    fit_nonRES->SetLineWidth(4);
    fit_nonRES->Draw("SAMEL");

    fit_total->SetLineColor(1); // Black
    fit_total->SetLineWidth(4);
    fit_total->Draw("SAMEL");

    h_mc_signal->SetLineWidth(3);
    h_mc_signal->SetLineColor(kRed);
    h_mc_signal->SetFillColor(kWhite);

    TLegend *legend = new TLegend(0.65,0.65,0.9,0.9);  
    legend->AddEntry(h_mc_signal, "Signal", "le");
    legend->AddEntry(fit_total, "Fit to Simulation", "l");
    legend->AddEntry(fit_deltaRES, "Fit: #Delta(1232) RES", "l");
    legend->AddEntry(fit_otherRES, "Fit: Other RES", "l");
    legend->AddEntry(fit_nonRES, "Fit: Non-RES", "l");
    legend->Draw();
    
    // ------------------------------------------------------------------------
    // Calc ChiSq
    // ------------------------------------------------------------------------
    // Arrays for Fit Results 
    nPoints = 20; // Last 4 bins are zero ( W > 1.8 GeV)
    double ChiSq_fit_x[nPoints];
    double ChiSq_data_y[nPoints];
    double ChiSq_fit_deltaRES_y[nPoints];
    double ChiSq_fit_otherRES_y[nPoints];
    double ChiSq_fit_nonRES_y[nPoints];
    double ChiSq_fit_nonRES_G1_y[nPoints];
    double ChiSq_fit_nonRES_G2_y[nPoints];
    double ChiSq_fit_total_y[nPoints];

    // X Axis (W)
    for (int i = 7; i <= nPoints+6; ++i){
        ChiSq_fit_x[i-7] = h_mc_signal->GetBinCenter(i); 
        ChiSq_data_y[i-7] = h_mc_signal->GetBinContent(i); 
    }

    Get_BreitWigner(pars_MC_deltaRES, ChiSq_fit_x, ChiSq_fit_deltaRES_y, nPoints);
    Get_Gaussian(pars_MC_otherRES, ChiSq_fit_x, ChiSq_fit_otherRES_y, nPoints);
    Get_Gaussian(pars_MC_nonRES_G1, ChiSq_fit_x, ChiSq_fit_nonRES_G1_y, nPoints);
    Get_Gaussian(pars_MC_nonRES_G2, ChiSq_fit_x, ChiSq_fit_nonRES_G2_y, nPoints);
   
    for (int i = 0; i < nPoints; ++i){
        ChiSq_fit_nonRES_y[i] = ChiSq_fit_nonRES_G1_y[i] + ChiSq_fit_nonRES_G2_y[i];
        ChiSq_fit_total_y[i] = ChiSq_fit_deltaRES_y[i] + ChiSq_fit_otherRES_y[i] + ChiSq_fit_nonRES_y[i];
    }
   
    double ChiSq_dof = Calc_ChiSq_dof(ChiSq_data_y, ChiSq_fit_total_y, nPoints, 9);
 
    // Add Text
    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.03);
    text.DrawLatex(0.20,0.85,Form("#chi^{2}/dof = %3.2f", ChiSq_dof));
    text.DrawLatex(0.20,0.82,Form("dof = %d-%d", nPoints, 9));

    // Plot 
    c->Update();
    c->Print(Form("%s%s",plotDir.c_str(),"W_Fit_MC.png"), "png");

    delete h_mc_signal;
    delete f_mc;
    delete fit_deltaRES;
    delete fit_otherRES;
    delete fit_nonRES;
    delete fit_total;
    delete plotter;
    delete legend;
    delete c;

}

void CCProtonPi0_Plotter::Get_BreitWigner(double* pars, double* x, double* y, int nPoints)
{
    for (int i = 0; i < nPoints; ++i){  
       y[i] = pars[0] * std::pow((x[i]/0.938), 3) * ( (pars[1]*pars[2]) / (std::pow((x[i]*x[i] - pars[1]*pars[1]), 2) + std::pow((pars[1]*pars[2]),2)));
    }
}

void CCProtonPi0_Plotter::Get_Gaussian(double* pars, double* x, double* y, int nPoints)
{
    for (int i = 0; i < nPoints; ++i){  
        y[i] = pars[0] * std::exp(-(std::pow(((x[i]-pars[1])/pars[2]),2)));  
    }
}

#endif

