#include "CCProtonPi0_Plotter.h"

using namespace PlotUtils;


void CCProtonPi0_Plotter::Draw_QSq_EnuLimit()
{
    // Draw Data vs MC Stacked
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    // Draw Fit Results -- Fit Results from MATLAB
    double pars_LowEnu[6];
    pars_LowEnu[0] = 53.49;
    pars_LowEnu[1] = -2.638; 

    // 1 Sigma Down
    pars_LowEnu[2] = 47.0723;
    pars_LowEnu[3] = -2.9084;

    // 1 Sigma Up
    pars_LowEnu[4] = 59.9148;
    pars_LowEnu[5] = -2.3684;

    std::string data_dir = "/minerva/data/users/oaltinok/NTupleAnalysis_QSq_LowEnu/Data/Analyzed/CrossSection.root";
    std::string mc_dir = "/minerva/data/users/oaltinok/NTupleAnalysis_QSq_LowEnu/MC/Analyzed/CrossSection.root";
    Draw_QSq_EnuFit(data_dir, mc_dir, "LowEnu", pars_LowEnu);

    double pars_HighEnu[6];
    pars_HighEnu[0] = 18.27;
    pars_HighEnu[1] = -0.8254; 

    // 1 Sigma Down
    pars_HighEnu[2] = 15.9141;
    pars_HighEnu[3] = -1.0026;

    // 1 Sigma Up
    pars_HighEnu[4] = 20.6191;
    pars_HighEnu[5] = -0.6482;
 
    data_dir = "/minerva/data/users/oaltinok/NTupleAnalysis_QSq_HighEnu/Data/Analyzed/CrossSection.root";
    mc_dir = "/minerva/data/users/oaltinok/NTupleAnalysis_QSq_HighEnu/MC/Analyzed/CrossSection.root";
    Draw_QSq_EnuFit(data_dir, mc_dir, "HighEnu", pars_HighEnu);
}

void CCProtonPi0_Plotter::Draw_QSq_EnuFit(std::string data_dir, std::string mc_dir, std::string var_name, double* pars)
{
    // ------------------------------------------------------------------------
    // First Plot Data vs MC using MnvPlotter
    // ------------------------------------------------------------------------
    TFile* f_data = new TFile(data_dir.c_str());
    TFile* f_mc = new TFile(mc_dir.c_str());

    MnvH1D* h_data = GetMnvH1D(f_data,"QSq_xsec" );
    MnvH1D* h_mc = GetMnvH1D(f_mc, "QSq_xsec" );

    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas("c","c",1280,800);

    double mc_ratio = 1.0; // Drawing Cross Sections

    //TH1D* hist_data = GetBinNormalizedTH1D(h_data);
    //printBins(hist_data, "BinNormalized", false);

    plotter->headroom = 1.75;
    plotter->data_line_width = 2;
    plotter->data_marker_size = 1.5;
    gStyle->SetEndErrorSize(6);
    plotter->DrawDataMCWithErrorBand(h_data, h_mc, mc_ratio, "N", false, NULL, NULL, false, true);

    // ------------------------------------------------------------------------
    // Second Add Fit Results
    // ------------------------------------------------------------------------
    int nPoints = 201;
    double fit_x[nPoints];
    double fit_y[nPoints];
    double fit_y_1Sigma_dn[nPoints];
    double fit_y_1Sigma_up[nPoints];

    // X Axis (W)
    double xStart = 0.0;
    double inc = 0.01;
    for (int i = 0; i < nPoints; ++i){
        fit_x[i] = xStart + i*inc; 
    }

    Get_Exponential(fit_x, fit_y, nPoints, pars[0], pars[1]);
    Get_Exponential(fit_x, fit_y_1Sigma_dn, nPoints, pars[2], pars[3]);
    Get_Exponential(fit_x, fit_y_1Sigma_up, nPoints, pars[4], pars[5]);

    TGraph* fit_CV = new TGraph(nPoints, fit_x, fit_y);
    TGraph* fit_1Sigma_dn = new TGraph(nPoints, fit_x, fit_y_1Sigma_dn);
    TGraph* fit_1Sigma_up = new TGraph(nPoints, fit_x, fit_y_1Sigma_up);

    fit_CV->SetLineColor(kBlue);
    fit_CV->SetLineWidth(3);
    fit_CV->Draw("SAMEL");

    fit_1Sigma_dn->SetLineColor(kBlue);
    fit_1Sigma_dn->SetLineWidth(2);
    fit_1Sigma_dn->SetLineStyle(7);
    fit_1Sigma_dn->Draw("SAMEL");

    fit_1Sigma_up->SetLineColor(kBlue);
    fit_1Sigma_up->SetLineWidth(2);
    fit_1Sigma_up->SetLineStyle(7);
    fit_1Sigma_up->Draw("SAMEL");

    // Add Legend
    h_data->SetMarkerStyle(20);
    h_data->SetMarkerSize(1);
    h_data->SetMarkerColor(kBlack);
    h_data->SetLineWidth(2);
    h_data->SetLineColor(kBlack);

    h_mc->SetLineWidth(3);
    h_mc->SetLineColor(kRed);
    h_mc->SetFillColor(kWhite);

    TLegend *legend = new TLegend(0.65,0.75,0.9,0.9);  
    legend->AddEntry(h_data, "Data (3.33e20 POT)", "lep");
    legend->AddEntry(h_mc, "Simulation", "l");
    legend->AddEntry(fit_CV, "Fit to Data", "l");
    legend->AddEntry(fit_1Sigma_dn, "Fit #pm1#sigma", "l");
    legend->Draw();
 
    // Plot 
    std::string plotDir = Folder_List::plotDir_OtherStudies;
    c->Update();
    c->Print(Form("%s%s%s%s",plotDir.c_str(),"QSq_Enu_Fit_", var_name.c_str(), ".png"), "png");

    delete h_data;
    delete h_mc;
    delete f_data;
    delete f_mc;
    delete fit_CV;
    delete fit_1Sigma_dn;
    delete fit_1Sigma_up;
    delete plotter;
    delete legend;
    delete c;
}


void CCProtonPi0_Plotter::Get_Exponential(double* x, double* y, int nPoints, double a, double b)
{
    for (int i = 0; i < nPoints; ++i){  
        y[i] = a * std::exp(b * x[i]);  
    }
}

void CCProtonPi0_Plotter::Draw_QSq_MaRES_Plots()
{
    // Read Files
    std::string data_dir_SB_Michel = "/minerva/data/users/oaltinok/NTupleAnalysis/Data/Analyzed/CutHistograms_Michel.root";
    std::string data_dir_SB_pID = "/minerva/data/users/oaltinok/NTupleAnalysis/Data/Analyzed/CutHistograms_pID.root";
    std::string data_dir_SB_LowInvMass = "/minerva/data/users/oaltinok/NTupleAnalysis/Data/Analyzed/CutHistograms_LowInvMass.root";
    std::string data_dir_SB_HighInvMass = "/minerva/data/users/oaltinok/NTupleAnalysis/Data/Analyzed/CutHistograms_HighInvMass.root";
    std::string data_dir = "/minerva/data/users/oaltinok/NTupleAnalysis/Data/Analyzed/Interaction.root";
 
    std::string mc_dir_SB_Michel = "/minerva/data/users/oaltinok/NTupleAnalysis/MC/Analyzed/CutHistograms_Michel.root";
    std::string mc_dir_SB_pID = "/minerva/data/users/oaltinok/NTupleAnalysis/MC/Analyzed/CutHistograms_pID.root";
    std::string mc_dir_SB_LowInvMass = "/minerva/data/users/oaltinok/NTupleAnalysis/MC/Analyzed/CutHistograms_LowInvMass.root";
    std::string mc_dir_SB_HighInvMass = "/minerva/data/users/oaltinok/NTupleAnalysis/MC/Analyzed/CutHistograms_HighInvMass.root";
    std::string mc_dir = "/minerva/data/users/oaltinok/NTupleAnalysis/MC/Analyzed/Interaction.root";
   
    TFile* f_data_SB_Michel = new TFile(data_dir_SB_Michel.c_str());
    TFile* f_data_SB_pID = new TFile(data_dir_SB_pID.c_str());
    TFile* f_data_SB_LowInvMass = new TFile(data_dir_SB_LowInvMass.c_str());
    TFile* f_data_SB_HighInvMass = new TFile(data_dir_SB_HighInvMass.c_str());
    TFile* f_data = new TFile(data_dir.c_str());
 
    TFile* f_mc_SB_Michel = new TFile(mc_dir_SB_Michel.c_str());
    TFile* f_mc_SB_pID = new TFile(mc_dir_SB_pID.c_str());
    TFile* f_mc_SB_LowInvMass = new TFile(mc_dir_SB_LowInvMass.c_str());
    TFile* f_mc_SB_HighInvMass = new TFile(mc_dir_SB_HighInvMass.c_str());
    TFile* f_mc = new TFile(mc_dir.c_str());


    // Get Histograms
    MnvH1D* h_data_SB_Michel = GetMnvH1D(f_data_SB_Michel, "SideBand_QSq_0");
    MnvH1D* h_data_SB_pID = GetMnvH1D(f_data_SB_pID, "SideBand_QSq_0");
    MnvH1D* h_data_SB_LowInvMass = GetMnvH1D(f_data_SB_LowInvMass, "SideBand_QSq_0");
    MnvH1D* h_data_SB_HighInvMass = GetMnvH1D(f_data_SB_HighInvMass, "SideBand_QSq_0");
    MnvH1D* h_data = GetMnvH1D(f_data, "QSq_MaRES_0");

    MnvH1D* h_mc_SB_Michel = GetMnvH1D(f_mc_SB_Michel, "SideBand_QSq_0");
    MnvH1D* h_mc_SB_pID = GetMnvH1D(f_mc_SB_pID, "SideBand_QSq_0");
    MnvH1D* h_mc_SB_LowInvMass = GetMnvH1D(f_mc_SB_LowInvMass, "SideBand_QSq_0");
    MnvH1D* h_mc_SB_HighInvMass = GetMnvH1D(f_mc_SB_HighInvMass, "SideBand_QSq_0");
    MnvH1D* h_mc = GetMnvH1D(f_mc, "QSq_MaRES_0");

    // Plot
    std::string plotDir = Folder_List::plotDir_OtherStudies;
    DrawDataMC(h_data_SB_Michel, h_mc_SB_Michel, "QSq_Michel", plotDir, false);
    DrawDataMC(h_data_SB_pID, h_mc_SB_pID, "QSq_pID", plotDir, false);
    DrawDataMC(h_data_SB_LowInvMass, h_mc_SB_LowInvMass, "QSq_LowInvMass", plotDir, false);
    DrawDataMC(h_data_SB_HighInvMass, h_mc_SB_HighInvMass, "QSq_HighInvMass", plotDir, false);
    DrawDataMC(h_data, h_mc, "QSq_MaRES", plotDir, false);

    // Clean Memory
    delete h_data_SB_Michel;
    delete h_data_SB_pID;
    delete h_data_SB_LowInvMass;
    delete h_data_SB_HighInvMass;
    delete h_data;

    delete h_mc_SB_Michel;
    delete h_mc_SB_pID;
    delete h_mc_SB_LowInvMass;
    delete h_mc_SB_HighInvMass;
    delete h_mc;

    delete f_data_SB_Michel;
    delete f_data_SB_pID;
    delete f_data_SB_LowInvMass;
    delete f_data_SB_HighInvMass;
    delete f_data;

    delete f_mc_SB_Michel;
    delete f_mc_SB_pID;
    delete f_mc_SB_LowInvMass;
    delete f_mc_SB_HighInvMass;
    delete f_mc;
}



void CCProtonPi0_Plotter::Draw_QSq_MaRES_AreaNorm()
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    TFile* f_mc = new TFile(rootDir_Interaction.mc.c_str());

    std::string var_name = "QSq_All";
    std::string var = Form("%s_%d",var_name.c_str(),1);

    MnvH1D* data = GetBckgSubtractedData(rootDir_Interaction, var_name, 2997.1);
    MnvH1D* mc = GetMnvH1D(f_mc, var);

    double area_data = data->Integral(3,20);
    double area_mc = mc->Integral(3,20);

    mc->Scale(area_data/area_mc);

    data->SetMarkerStyle(20);
    data->SetMarkerSize(1);
    data->SetMarkerColor(kBlack);
    data->SetLineWidth(2);
    data->SetLineColor(kBlack);

    mc->SetLineWidth(3);
    mc->SetLineColor(kRed);
    mc->SetFillColor(kWhite);

    TCanvas* c = new TCanvas("c","c",1280,800);

    data->SetMaximum(data->GetMaximum()*1.5);
    data->Draw("E1 X0");
    mc->Draw("HIST SAME");

    // TLegend
    TLegend *legend = new TLegend(0.65,0.80,0.9,0.9);  
    legend->AddEntry(data, "Bckg Subt. Data CV ", "lep");
    legend->AddEntry(mc, "GENIE Signal CV", "l");
    legend->Draw();
 
    // Add Text
    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.03);
    text.DrawLatex(0.65,0.70,Form("%s", "Bins used [3, 20]"));
    text.DrawLatex(0.65,0.66,Form("%s%3.2f", "Area(Data/MC) = ", area_data/area_mc));
 
    // Save Plot 
    gStyle->SetOptStat(0); 
    c->Update();
    c->Print(Form("%s%s",plotDir.c_str(),"QSq_MaRES_AreaNorm.png"), "png");

    delete legend;
    delete c;
    delete data;
    delete mc;
    delete f_mc;
}

void CCProtonPi0_Plotter::Draw_QSq_MaRES_Fit(bool isAreaNorm)
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;
 
    std::string data_dir = "/minerva/data/users/oaltinok/NTupleAnalysis_MaRES_Fit_XSecs/Data/Analyzed/CrossSection.root";
    std::string mc_dir = "/minerva/data/users/oaltinok/NTupleAnalysis_MaRES_Fit_XSecs/MC/Analyzed/CrossSection.root";

    TFile* f_data = new TFile(data_dir.c_str());
    TFile* f_mc = new TFile(mc_dir.c_str());

    CCProtonPi0_QSqFitter QSqFitter;
    int ind = QSqFitter.GetMinChiSq(isAreaNorm);

    // --------------------------------------------------------------------
    // Get Histograms 
    // --------------------------------------------------------------------
    // Central Value
    MnvH1D* data_cv = GetMnvH1D(f_data, "QSq_xsec");
    MnvH1D* mc_cv = GetMnvH1D(f_mc, "QSq_xsec");

    std::string norm_label;
    if (isAreaNorm) norm_label = "Area"; 
    else norm_label = "POT"; 

    TH1D* h_data_cv = GetBinNormalizedTH1D(data_cv, true);
    TH1D* h_mc_cv = GetBinNormalizedTH1D(mc_cv);

    // Best MaRES -- Lowest Global ChiSq
    MnvVertErrorBand* err_data = data_cv->GetVertErrorBand("HighMaRES");
    MnvVertErrorBand* err_mc = mc_cv->GetVertErrorBand("HighMaRES");

    std::vector<TH1D*> unv_data = err_data->GetHists();
    std::vector<TH1D*> unv_mc = err_mc->GetHists();

    TH1D* data_best = (TH1D*)unv_data[ind]->Clone();
    TH1D* mc_best = (TH1D*)unv_mc[ind]->Clone();
   
    // Bin Normalization 
    double norm_bin_width = GetSmallestBinWidth(data_cv);
    data_best->Scale(norm_bin_width, "width");
    mc_best->Scale(norm_bin_width, "width");

    // Area Normalization
    int nBins = h_data_cv->GetNbinsX();
    if (isAreaNorm){
        double area_data = h_data_cv->Integral(3, nBins, "width");
        double area_mc = h_mc_cv->Integral(3, nBins, "width");
        double mc_ratio = area_data/area_mc;

        h_mc_cv->Scale(mc_ratio);
    }
   
    if (isAreaNorm){
        double area_data = data_best->Integral(3, nBins, "width");
        double area_mc = mc_best->Integral(3, nBins, "width");
        double mc_ratio = area_data/area_mc;

        mc_best->Scale(mc_ratio);
    }

    h_data_cv->SetMarkerStyle(20);
    h_data_cv->SetMarkerSize(1);
    h_data_cv->SetMarkerColor(kBlack);
    h_data_cv->SetLineWidth(1);
    h_data_cv->SetLineColor(kBlack);

    h_mc_cv->SetLineWidth(3);
    h_mc_cv->SetLineColor(kRed);
    h_mc_cv->SetFillColor(kWhite);

    mc_best->SetLineWidth(3);
    mc_best->SetLineColor(kBlue);
    mc_best->SetLineStyle(1);
    mc_best->SetFillColor(kWhite);

    data_best->SetMarkerStyle(20);
    data_best->SetMarkerSize(1);
    data_best->SetMarkerColor(kBlue);
    data_best->SetLineWidth(1);
    data_best->SetLineStyle(1);
    data_best->SetLineColor(kBlue);

    TCanvas* c = new TCanvas("c","c",640,480);

    if (thesisStyle){
        h_data_cv->GetXaxis()->SetTitleFont(62);
        h_data_cv->GetXaxis()->SetTitleSize(0.06);
        h_data_cv->GetXaxis()->CenterTitle();
        h_data_cv->GetXaxis()->SetTitleOffset(1.15);
        h_data_cv->GetXaxis()->SetLabelFont(42);
        h_data_cv->GetXaxis()->SetLabelSize(0.05);
        h_data_cv->GetXaxis()->SetNdivisions(408);

        h_data_cv->GetYaxis()->SetTitleFont(62);
        h_data_cv->GetYaxis()->SetTitleSize(0.06);
        //h_data_cv->GetYaxis()->CenterTitle();
        h_data_cv->GetYaxis()->SetTitleOffset(1.1);
        h_data_cv->GetYaxis()->SetLabelFont(42);
        h_data_cv->GetYaxis()->SetLabelSize(0.05);
        TGaxis::SetMaxDigits(3);
        
        gStyle->SetEndErrorSize(4);
        gStyle->SetStripDecimals(false);
        gStyle->SetPadRightMargin(0.05);
    }

    h_data_cv->SetMaximum(data_cv->GetMaximum()*1.5);
    h_data_cv->Draw("E1 X0");
    h_mc_cv->Draw("HIST SAME");
    data_best->Draw("E1 X0 SAME");
    mc_best->Draw("HIST SAME");

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");

    // TLegend
    TLegend *legend = new TLegend(0.45,0.60,0.9,0.9);  
    ApplyStyle_Legend(legend);
    legend->AddEntry(h_data_cv, "Data (M_{A}^{RES} = 1.12 GeV)", "lep");
    legend->AddEntry(h_mc_cv, "GENIE (M_{A}^{RES} = 1.12 GeV)", "l");
    legend->AddEntry(data_best, "Data (M_{A}^{RES} = 1.50 GeV)", "lep");
    legend->AddEntry(mc_best, "GENIE (M_{A}^{RES} = 1.50 GeV)", "l");
    legend->Draw();
 
    // Add Text
    //text.DrawLatex(0.65,0.66,Form("%s%3.2f%s", "GENIE MaRES = ", 1.12," GeV"));
    //text.DrawLatex(0.65,0.62,Form("%s%3.2f", "GENIE MaRES #chi^{2} = ", QSqFitter.ChiSqVector_up[0]));
    //text.DrawLatex(0.65,0.58,Form("%s%3.2f%s", "Best MaRES = ", QSqFitter.MaRESVector_up[ind]," GeV"));
    //text.DrawLatex(0.65,0.54,Form("%s%3.2f", "Best MaRES #chi^{2} = ", QSqFitter.ChiSqVector_up[ind]));

    // Save Plot 
    gStyle->SetOptStat(0); 
    c->Update();
    c->Print(Form("%s%s_%s%s",plotDir.c_str(),"QSq_MaRES_Fit",norm_label.c_str(), ".pdf"), "pdf");

    delete data_cv;
    delete data_best;
    delete mc_cv;
    delete mc_best;
    delete c;
    delete f_data;
    delete f_mc;
}

void CCProtonPi0_Plotter::Draw_QSq_MaRES_Fit_SB()
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    TFile* f_data = new TFile(Folder_List::rootDir_Interaction_data.c_str());
    TFile* f_mc = new TFile(Folder_List::rootDir_Interaction_mc.c_str());

    CCProtonPi0_QSqFitter QSqFitter;
    int ind = QSqFitter.GetMinChiSq(false);

    // --------------------------------------------------------------------
    // Get Histograms 
    // --------------------------------------------------------------------
    // Central Value
    MnvH1D* data_cv = GetMnvH1D(f_data, "QSq_MaRES_0");
    MnvH1D* mc_cv = GetMnvH1D(f_mc, "QSq_MaRES_0");
    
    TH1D* h_data_cv = GetBinNormalizedTH1D(data_cv);
    TH1D* h_mc_cv = GetBinNormalizedTH1D(mc_cv);

    // Best MaRES -- Lowest Global ChiSq
    MnvVertErrorBand* err_data = data_cv->GetVertErrorBand("HighMaRES");
    MnvVertErrorBand* err_mc = mc_cv->GetVertErrorBand("HighMaRES");

    std::vector<TH1D*> unv_data = err_data->GetHists();
    std::vector<TH1D*> unv_mc = err_mc->GetHists();

    TH1D* data_best = (TH1D*)unv_data[ind]->Clone();
    TH1D* mc_best = (TH1D*)unv_mc[ind]->Clone();
   
    // Bin Normalization 
    double norm_bin_width = GetSmallestBinWidth(data_cv);
    data_best->Scale(norm_bin_width, "width");
    mc_best->Scale(norm_bin_width, "width");

    h_data_cv->SetMarkerStyle(20);
    h_data_cv->SetMarkerSize(1);
    h_data_cv->SetMarkerColor(kBlack);
    h_data_cv->SetLineWidth(2);
    h_data_cv->SetLineColor(kBlack);

    h_mc_cv->Scale(POT_ratio);
    h_mc_cv->SetLineWidth(3);
    h_mc_cv->SetLineColor(kRed);
    h_mc_cv->SetFillColor(kWhite);

    mc_best->Scale(POT_ratio);
    mc_best->SetLineWidth(3);
    mc_best->SetLineStyle(1);
    mc_best->SetLineColor(kGreen+2);
    mc_best->SetFillColor(kWhite);

    data_best->SetMarkerStyle(20);
    data_best->SetMarkerSize(1);
    data_best->SetLineStyle(1);
    data_best->SetMarkerColor(kBlue);
    data_best->SetLineWidth(2);
    data_best->SetLineColor(kBlue);

    TCanvas* c = new TCanvas("c","c",1280,800);

    h_data_cv->SetMaximum(data_cv->GetMaximum()*1.5);
    h_data_cv->Draw("E1 X0");
    h_mc_cv->Draw("HIST SAME");
    data_best->Draw("E1 X0 SAME");
    mc_best->Draw("HIST SAME");

    // TLegend
    TLegend *legend = new TLegend(0.65,0.75,0.9,0.9);  
    legend->AddEntry(h_data_cv, "Bckg Subt. Data CV ", "lep");
    legend->AddEntry(h_mc_cv, "GENIE Signal CV", "l");
    legend->AddEntry(data_best, "Bckg Subt. Data Best Fit", "lep");
    legend->AddEntry(mc_best, "GENIE Signal Best Fit", "l");
    legend->Draw();
 
    // Add Text
    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.03);
    text.DrawLatex(0.65,0.66,Form("%s%3.2f%s", "GENIE MaRES = ", 1.12," GeV"));
    text.DrawLatex(0.65,0.62,Form("%s%3.2f", "GENIE MaRES #chi^{2} = ", QSqFitter.ChiSqVector_up[0]));
    text.DrawLatex(0.65,0.58,Form("%s%3.2f%s", "Best MaRES = ", QSqFitter.MaRESVector_up[ind]," GeV"));
    text.DrawLatex(0.65,0.54,Form("%s%3.2f", "Best MaRES #chi^{2} = ", QSqFitter.ChiSqVector_up[ind]));

    // Save Plot 
    gStyle->SetOptStat(0); 
    c->Update();
    c->Print(Form("%s%s",plotDir.c_str(),"QSq_MaRES_SB_Fit.png"), "png");

    delete data_cv;
    delete data_best;
    delete mc_cv;
    delete mc_best;
    delete c;
    delete f_data;
    delete f_mc;
}

void CCProtonPi0_Plotter::Draw_QSq_DeltaSuppression_v2(std::string var_name)
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    std::string data_dir = "/minerva/data/users/oaltinok/NTupleAnalysis_Best/Data/Analyzed/CrossSection.root";
    std::string mc_dir = "/minerva/data/users/oaltinok/NTupleAnalysis_Best/MC/Analyzed/CrossSection.root";

    TFile* f_data_best = new TFile(rootDir_CrossSection.data.c_str());
    TFile* f_mc_best = new TFile(rootDir_CrossSection.mc.c_str());

    TFile* f_data_cv = new TFile(data_dir.c_str());
    TFile* f_mc_cv = new TFile(mc_dir.c_str());

    // --------------------------------------------------------------------
    // Get Histograms 
    // --------------------------------------------------------------------
    MnvH1D* data_cv = GetMnvH1D(f_data_cv, var_name);
    MnvH1D* mc_cv = GetMnvH1D(f_mc_cv, var_name);
  
    TH1D* h_data_cv = GetBinNormalizedTH1D(data_cv);
    TH1D* h_mc_cv = GetBinNormalizedTH1D(mc_cv);

    MnvH1D* data_best = GetMnvH1D(f_data_best, var_name);
    MnvH1D* mc_best = GetMnvH1D(f_mc_best, var_name);
  
    TH1D* h_data_best = GetBinNormalizedTH1D(data_best);
    TH1D* h_mc_best = GetBinNormalizedTH1D(mc_best);

    h_data_cv->SetMarkerStyle(20);
    h_data_cv->SetMarkerSize(1);
    h_data_cv->SetMarkerColor(kBlack);
    h_data_cv->SetLineWidth(2);
    h_data_cv->SetLineColor(kBlack);

    h_mc_cv->SetLineWidth(3);
    h_mc_cv->SetLineColor(kRed);
    h_mc_cv->SetFillColor(kWhite);

    h_data_best->SetMarkerStyle(20);
    h_data_best->SetLineStyle(1);
    h_data_best->SetMarkerSize(1);
    h_data_best->SetMarkerColor(kBlue);
    h_data_best->SetLineWidth(2);
    h_data_best->SetLineColor(kBlue);

    h_mc_best->SetLineWidth(3);
    h_mc_best->SetLineStyle(1);
    h_mc_best->SetLineColor(kBlue);
    h_mc_best->SetFillColor(kWhite);
    
    TCanvas* c = new TCanvas("c","c",800,800);

    if (thesisStyle){
        h_data_cv->GetXaxis()->SetTitleFont(62);
        h_data_cv->GetXaxis()->SetTitleSize(0.06);
        h_data_cv->GetXaxis()->CenterTitle();
        h_data_cv->GetXaxis()->SetTitleOffset(1.15);
        h_data_cv->GetXaxis()->SetLabelFont(42);
        h_data_cv->GetXaxis()->SetLabelSize(0.05);
        h_data_cv->GetXaxis()->SetNdivisions(408);

        h_data_cv->GetYaxis()->SetTitleFont(62);
        h_data_cv->GetYaxis()->SetTitleSize(0.06);
        //h_data_cv->GetYaxis()->CenterTitle();
        h_data_cv->GetYaxis()->SetTitleOffset(1.1);
        h_data_cv->GetYaxis()->SetLabelFont(42);
        h_data_cv->GetYaxis()->SetLabelSize(0.05);
        TGaxis::SetMaxDigits(3);
    }

    h_data_cv->SetMaximum(h_data_cv->GetMaximum()*1.75);
    h_data_cv->SetMinimum(0.0);
    h_data_cv->Draw("E1 X0");
    h_mc_cv->Draw("HIST SAME");
    h_data_best->Draw("E1 X0 SAME");
    h_mc_best->Draw("HIST SAME");

    // TLegend
    TLegend *legend = new TLegend(0.50,0.70,0.9,0.9);  
    legend->SetTextSize(0.03);
    legend->AddEntry(h_data_cv, "CV Data d#sigma/dQ^{2}", "lep");
    legend->AddEntry(h_mc_cv, "CV GENIE d#sigma/dQ^{2}", "l");
    legend->AddEntry(h_data_best, "#Delta Supp. Data d#sigma/dQ^{2}", "lep");
    legend->AddEntry(h_mc_best, "#Delta Supp. GENIE d#sigma/dQ^{2}", "l");
    legend->Draw();

    // Add Text
    int nBins = h_data_cv->GetNbinsX();
    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.03);
    text.DrawLatex(0.17,0.85,Form("%s%3.2f", "GENIE CV #chi^{2} = ", Calc_ChiSq(h_data_cv, h_mc_cv,1,nBins)));
    text.DrawLatex(0.17,0.81,Form("%s%3.2f", "#Delta Supp. #chi^{2} = ", Calc_ChiSq(h_data_best, h_mc_best,1,nBins)));

    // Save Plot 
    gStyle->SetOptStat(0); 
    c->Update();
    std::string plot_name = var_name + "_DeltaSuppressed.png";
    c->Print(Form("%s%s",plotDir.c_str(),plot_name.c_str()), "png");

    delete h_data_cv;
    delete data_cv;
    delete h_data_best;
    delete h_mc_cv;
    delete mc_cv;
    delete h_mc_best;
    delete c;
    delete f_data_cv;
    delete f_data_best;
    delete f_mc_cv;
    delete f_mc_best;
}

void CCProtonPi0_Plotter::Draw_QSq_DeltaSuppression()
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    TFile* f_data = new TFile(rootDir_CrossSection.data.c_str());
    TFile* f_mc = new TFile(rootDir_CrossSection.mc.c_str());

    // --------------------------------------------------------------------
    // Get Histograms 
    // --------------------------------------------------------------------
    MnvH1D* data_cv = GetMnvH1D(f_data, "QSq_xsec");
    MnvH1D* mc_cv = GetMnvH1D(f_mc, "QSq_xsec");
  
    TH1D* h_data_cv = GetBinNormalizedTH1D(data_cv);
    TH1D* h_mc_cv = GetBinNormalizedTH1D(mc_cv);

    // MINOS Factor -- First Universe 
    MnvVertErrorBand* err_data = data_cv->GetVertErrorBand("DeltaFactor");
    MnvVertErrorBand* err_mc = mc_cv->GetVertErrorBand("DeltaFactor");

    std::vector<TH1D*> unv_data = err_data->GetHists();
    std::vector<TH1D*> unv_mc = err_mc->GetHists();

    TH1D* h_data_MINOS= (TH1D*)unv_data[0]->Clone();
    TH1D* h_mc_MINOS = (TH1D*)unv_mc[0]->Clone();
 
    // Best Numbers -- Lowest ChiSq 
    CCProtonPi0_QSqFitter QSqFitter;
    int ind = QSqFitter.GetMinChiSq_DeltaFactor();
    TH1D* h_data_Best = (TH1D*)unv_data[ind]->Clone();
    TH1D* h_mc_Best = (TH1D*)unv_mc[ind]->Clone();
   
    // Bin Normalization 
    double norm_bin_width = GetSmallestBinWidth(data_cv);
    h_data_MINOS->Scale(norm_bin_width, "width");
    h_data_Best->Scale(norm_bin_width, "width");
    h_mc_MINOS->Scale(norm_bin_width, "width");
    h_mc_Best->Scale(norm_bin_width, "width");

    h_data_cv->SetMarkerStyle(20);
    h_data_cv->SetMarkerSize(1);
    h_data_cv->SetMarkerColor(kBlack);
    h_data_cv->SetLineWidth(2);
    h_data_cv->SetLineColor(kBlack);

    h_mc_cv->SetLineWidth(3);
    h_mc_cv->SetLineColor(kRed);
    h_mc_cv->SetFillColor(kWhite);

    h_mc_MINOS->SetLineWidth(3);
    h_mc_MINOS->SetLineColor(kBlue);
    h_mc_MINOS->SetFillColor(kWhite);

    h_data_MINOS->SetMarkerStyle(20);
    h_data_MINOS->SetMarkerSize(1);
    h_data_MINOS->SetMarkerColor(kBlue);
    h_data_MINOS->SetLineWidth(2);
    h_data_MINOS->SetLineColor(kBlue);

    h_mc_Best->SetLineWidth(3);
    h_mc_Best->SetLineStyle(1);
    h_mc_Best->SetLineColor(kGreen+2);
    h_mc_Best->SetFillColor(kWhite);

    h_data_Best->SetMarkerStyle(20);
    h_data_Best->SetLineStyle(1);
    h_data_Best->SetMarkerSize(1);
    h_data_Best->SetMarkerColor(kGreen+2);
    h_data_Best->SetLineWidth(2);
    h_data_Best->SetLineColor(kGreen+2);

    TCanvas* c = new TCanvas("c","c",1280,800);

    h_data_cv->SetMaximum(h_data_cv->GetMaximum()*1.5);
    h_data_cv->SetMinimum(0.0);
    h_data_cv->Draw("E1 X0");
    h_mc_cv->Draw("HIST SAME");
    h_data_MINOS->Draw("E1 X0 SAME");
    h_mc_MINOS->Draw("HIST SAME");
    h_data_Best->Draw("E1 X0 SAME");
    h_mc_Best->Draw("HIST SAME");


    // TLegend
    TLegend *legend = new TLegend(0.50,0.70,0.9,0.9);  
    legend->AddEntry(h_data_cv, "CV Data d#sigma/dQ^{2}", "lep");
    legend->AddEntry(h_mc_cv, "CV GENIE d#sigma/dQ^{2}", "l");
    legend->AddEntry(h_data_MINOS, "#Delta Suppressed (MINOS) Data d#sigma/dQ^{2}", "lep");
    legend->AddEntry(h_mc_MINOS, "#Delta Suppressed (MINOS) GENIE d#sigma/dQ^{2}", "l");
    legend->AddEntry(h_data_Best, "#Delta Suppressed (Best) Data d#sigma/dQ^{2}", "lep");
    legend->AddEntry(h_mc_Best, "#Delta Suppressed (Best) GENIE d#sigma/dQ^{2}", "l");
    legend->Draw();

     // Add Text
    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.03);
    text.DrawLatex(0.55,0.66,Form("%s%3.2f", "GENIE CV #chi^{2} = ", Calc_ChiSq(h_data_cv, h_mc_cv,1,3)));
    text.DrawLatex(0.55,0.62,Form("%s%3.2f", "#Delta Suppressed (MINOS) #chi^{2} = ", Calc_ChiSq(h_data_MINOS, h_mc_MINOS,1,3)));
    text.DrawLatex(0.55,0.58,Form("%s%3.2f", "#Delta Suppressed (Best) #chi^{2} = ", Calc_ChiSq(h_data_Best, h_mc_Best,1,3)));

    // Save Plot 
    gStyle->SetOptStat(0); 
    c->Update();
    c->Print(Form("%s%s",plotDir.c_str(),"QSq_DeltaSuppressed.png"), "png");

    delete h_data_cv;
    delete data_cv;
    delete h_data_MINOS;
    delete h_data_Best;
    delete h_mc_cv;
    delete mc_cv;
    delete h_mc_MINOS;
    delete h_mc_Best;
    delete c;
    delete f_data;
    delete f_mc;
}

void CCProtonPi0_Plotter::Draw_QSq_DeltaSuppression_AllPlots()
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    TFile* f_data = new TFile(rootDir_CrossSection.data.c_str());
    TFile* f_mc = new TFile(rootDir_CrossSection.mc.c_str());

    // --------------------------------------------------------------------
    // Get Histograms 
    // --------------------------------------------------------------------
    MnvH1D* data_cv = GetMnvH1D(f_data, "QSq_xsec");
    MnvH1D* mc_cv = GetMnvH1D(f_mc, "QSq_xsec");
  
    // MINOS Factor -- First Universe 
    MnvVertErrorBand* err_data = data_cv->GetVertErrorBand("DeltaFactor");
    MnvVertErrorBand* err_mc = mc_cv->GetVertErrorBand("DeltaFactor");

    std::vector<TH1D*> unv_data = err_data->GetHists();
    std::vector<TH1D*> unv_mc = err_mc->GetHists();

    // Bin Normalization 
    double norm_bin_width = GetSmallestBinWidth(data_cv);
    for (unsigned int i = 0; i < unv_mc.size(); ++i){
        TH1D* h_data = (TH1D*)unv_data[i]->Clone();
        TH1D* h_mc = (TH1D*)unv_mc[i]->Clone();

        h_data->Scale(norm_bin_width, "width");
        h_mc->Scale(norm_bin_width, "width");
        
        h_data->SetMarkerStyle(20);
        h_data->SetMarkerSize(1);
        h_data->SetMarkerColor(kBlack);
        h_data->SetLineWidth(2);
        h_data->SetLineColor(kBlack);

        h_mc->SetLineWidth(3);
        h_mc->SetLineStyle(1);
        h_mc->SetLineColor(kRed);
        h_mc->SetFillColor(kWhite);

        TCanvas* c = new TCanvas("c","c",1280,800);

        h_data->SetMaximum(h_data->GetMaximum()*1.5);
        h_data->SetMinimum(0.0);
        h_data->Draw("E1 X0");
        h_mc->Draw("HIST SAME");

        // TLegend
        TLegend *legend = new TLegend(0.65,0.80,0.9,0.9);  
        legend->AddEntry(h_data, "Data d#sigma/dQ^{2}", "lep");
        legend->AddEntry(h_mc, "GENIE d#sigma/dQ^{2}", "l");
        legend->Draw();

        // Add Text
        TLatex text;
        text.SetNDC();
        text.SetTextSize(0.03);
        text.DrawLatex(0.65,0.70,Form("%s%3.2f", "All Bins #chi^{2} = ", Calc_ChiSq(h_data, h_mc, 1, h_data->GetNbinsX())));
        text.DrawLatex(0.65,0.66,Form("%s%3.2f", "Low Bins #chi^{2} = ", Calc_ChiSq(h_data, h_mc, 1, 2)));

        // Save Plot 
        gStyle->SetOptStat(0); 
        c->Update();
        std::string plot_number;
        if (i < 10) plot_number = "00" + std::to_string((long long int) i);
        else if (i < 100) plot_number = "0" + std::to_string((long long int)i);
        else plot_number = std::to_string((long long int)i);
        c->Print(Form("%s%s_%s%s",plotDir.c_str(),"QSq_DeltaSuppressed",plot_number.c_str(),".png"), "png");

        delete h_data;
        delete h_mc;
        delete legend;
        delete c;
    }
    
    delete mc_cv;
    delete data_cv;
    delete f_data;
    delete f_mc;
}

double CCProtonPi0_Plotter::user_expo(double* x, double* par)
{
    return par[0]*exp(-par[1]*x[0]);
}

void CCProtonPi0_Plotter::expo_fit(const std::string& fileName, const std::string& histName)
{
    gStyle->SetErrorX(0.5);

    TFile* f = new TFile(fileName.c_str(), "read");
    MnvH1D* mnvh1d = static_cast<MnvH1D*>(f->Get(histName.c_str()));

    //assert(mnvh1d);

    std::cout << "NormBinWidth: " << mnvh1d->GetNormBinWidth() << std::endl;
    
    TCanvas* c1 = new TCanvas("c1");
    TH1D* h1d = (TH1D*) mnvh1d->GetBinNormalizedCopy().GetCVHistoWithError(true,false).Clone("h1d");

    for (int i = 1; i <= h1d->GetNbinsX(); ++i) {
        printf("\t g->SetPoint(%d, %10.2e, %10.2e);\n", i, h1d->GetBinCenter(i), h1d->GetBinContent(i));
    }

    for (int i = 1; i <= h1d->GetNbinsX(); ++i) {
        printf("\t g->SetPointError(%d, %10.2e);\n", i, h1d->GetBinError(i));
    }

    h1d->SetStats(false);
    h1d->SetMarkerStyle(20);
    h1d->SetMarkerSize(1.5);
    h1d->SetLineColor(kBlack);
    h1d->SetMarkerColor(kBlack);
    h1d->Draw("E1X0");
    h1d->GetYaxis()->SetRangeUser(0.001,30.0);
    h1d->GetYaxis()->SetTitleSize(0.06);
    h1d->GetYaxis()->SetLabelSize(0.06);
    h1d->GetYaxis()->CenterTitle();

    h1d->GetXaxis()->SetTitleSize(0.06);
    h1d->GetXaxis()->SetLabelSize(0.06);
    
    //TF1* f1_expo = new TF1("user_expo", user_expo,0.0,2.0,2);
    //h1d->Fit(f1_expo, "0", "", 0.25, 2.0);
    //
    //f1_expo->SetLineColor(kOrange-3);
    //f1_expo->Draw("same");
    //
    //TPaveText* expo_label = new TPaveText(0.6, 0.6, 0.8, 0.8, "NDC");
    //expo_label->AddText("f(x) =  a exp(-bx)");
    //expo_label->AddText(Form("a = %1.2f #pm %1.2f",  f1_expo->GetParameter(0), f1_expo->GetParError(0)));
    //expo_label->AddText(Form("b = %1.2f #pm %1.2f",  f1_expo->GetParameter(1), f1_expo->GetParError(1)));
    //expo_label->SetFillColor(0);
    //expo_label->SetTextColor(kBlue);
    //expo_label->Draw("same");

    //TH1D* expo_error_band = new TH1D("q2_fit", "Error band", 100, 0.25, 2.00); // same as fitting range
    //TVirtualFitter::GetFitter()->GetConfidenceIntervals(expo_error_band, 0.6827);
    //expo_error_band->SetLineWidth(3);
    //expo_error_band->SetLineColor(kOrange-3);
    //expo_error_band->SetFillColor(kOrange-3);
    //expo_error_band->SetFillStyle(3001);
    //expo_error_band->SetMarkerStyle(kDot);
    //expo_error_band->SetMarkerColor(kOrange-3);
    //expo_error_band->Draw("e2 same");
    
    c1->Print("q2-data-mc-xs-multiflux-pot-fit.png");

    
}


