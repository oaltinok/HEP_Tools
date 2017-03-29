// ================================================================
// ROOT Fitter -- Thanks to Trung!
// ================================================================

#include <iostream>
#include <cassert>
#include <string>
#include <vector>
#include <numeric>

#include <TFile.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>
#include <TF1.h>
#include <TPaveText.h>
#include <TVirtualFitter.h>
#include <TStyle.h>

#include <PlotUtils/MnvH1D.h>
#include <PlotUtils/MnvPlotter.h>
#include "Cintex/Cintex.h"

using namespace PlotUtils;

double user_expo(double* x, double* par)
{
    return par[0]*exp(-par[1]*x[0]);
}

void ApplyStyle(MnvPlotter* plotter)
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
}

void Fit(std::string rootDir_data, std::string rootDir_mc, std::string plotName)
{
    std::string histName = "QSq_xsec";
    //std::string histName = "q2-xsec-data";

    TFile* f_data = new TFile(rootDir_data.c_str(), "read");
    TFile* f_mc = new TFile(rootDir_mc.c_str(), "read");

    MnvH1D* data = new MnvH1D( * dynamic_cast<MnvH1D*>(f_data->Get(histName.c_str())) );
    MnvH1D* mc = new MnvH1D( * dynamic_cast<MnvH1D*>(f_mc->Get(histName.c_str())) );
   
    // ------------------------------------------------------------------------
    //  Plot Data MC using MnvPlotter
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
    ApplyStyle(plotter);

    TCanvas* c1 = new TCanvas("c1");
 
    data->SetTitle("Data (3.33e20 POT)");
    data->GetXaxis()->SetTitle("Q^{2} (GeV^{2})");
    mc->SetTitle("Simulation");

    std::cout<<"Data Integral = "<<data->Integral()<<std::endl;
    // Remove Error Bands from MC
    mc->ClearAllErrorBands();

    plotter->headroom = 1.75;
    plotter->DrawDataMCWithErrorBand(data, mc, 1.0, "TR", true, NULL, NULL, false, true, false);

    // ------------------------------------------------------------------------
    // Fit Data Points
    // ------------------------------------------------------------------------
    TH1D* h1d = (TH1D*) data->GetBinNormalizedCopy().GetCVHistoWithError(true,false).Clone("h1d");

    for (int i = 1; i <= h1d->GetNbinsX(); ++i) {
        printf("\t g->SetPoint(%d, %10.2e, %10.2e);\n", i, h1d->GetBinCenter(i), h1d->GetBinContent(i));
    }

    for (int i = 1; i <= h1d->GetNbinsX(); ++i) {
        printf("\t g->SetPointError(%d, %10.2e);\n", i, h1d->GetBinError(i));
    }
    
    TF1* f1_expo = new TF1("user_expo", user_expo,0.0,2.0,2);
    h1d->Fit(f1_expo, "0", "", 0.25, 2.0);
    
    f1_expo->SetLineColor(kOrange-3);
    f1_expo->Draw("same");
    
    TPaveText* expo_label = new TPaveText(0.6, 0.5, 0.8, 0.7, "NDC");
    expo_label->AddText("f(x) =  a exp(-bx)");
    expo_label->AddText(Form("a = %1.2f #pm %1.2f",  f1_expo->GetParameter(0), f1_expo->GetParError(0)));
    expo_label->AddText(Form("b = %1.2f #pm %1.2f",  f1_expo->GetParameter(1), f1_expo->GetParError(1)));
    expo_label->SetFillColor(0);
    expo_label->SetTextColor(kBlue);
    expo_label->Draw("same");

    TH1D* expo_error_band = new TH1D("q2_fit", "Error band", 100, 0.25, 2.00); // same as fitting range
    TVirtualFitter::GetFitter()->GetConfidenceIntervals(expo_error_band, 0.6827);
    expo_error_band->SetLineWidth(3);
    expo_error_band->SetLineColor(kOrange-3);
    expo_error_band->SetFillColor(kOrange-3);
    expo_error_band->SetFillStyle(3001);
    expo_error_band->SetMarkerStyle(kDot);
    expo_error_band->SetMarkerColor(kOrange-3);
    expo_error_band->Draw("e2 same");
   
    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");

    std::string output_folder = "/minerva/app/users/oaltinok/cmtuser/Minerva_v10r8p9/Ana/CCProtonPi0/NTupleAnalysis/Output/Plots/OtherStudies/";
    std::string output = output_folder + plotName;    
 
    c1->Print(output.c_str(),"pdf");
}

void Fit_QSq()
{
    ROOT::Cintex::Cintex::Enable();

    std::string data_Trung_HighEnu = "/minerva/app/users/ltrung/cmtuser/Minerva_v10r6p13/Ana/CCPi0/ana/make_hists/q2-xsec-data-high-enu.root"; 
    std::string data_Trung_LowEnu = "/minerva/app/users/ltrung/cmtuser/Minerva_v10r6p13/Ana/CCPi0/ana/make_hists/q2-xsec-data-low-enu.root"; 

    //std::string data_rootDir_HighEnu = "/minerva/data/users/oaltinok/NTupleAnalysis_QSq_HighEnu_4/Data/Analyzed/CrossSection.root";
    //std::string data_rootDir_LowEnu  = "/minerva/data/users/oaltinok/NTupleAnalysis_QSq_LowEnu_4/Data/Analyzed/CrossSection.root";
    //std::string mc_rootDir_HighEnu = "/minerva/data/users/oaltinok/NTupleAnalysis_QSq_HighEnu_4/MC/Analyzed/CrossSection.root";
    //std::string mc_rootDir_LowEnu  = "/minerva/data/users/oaltinok/NTupleAnalysis_QSq_LowEnu_4/MC/Analyzed/CrossSection.root";

    std::string plot_HighEnu = "QSq_HighEnu_ROOT_Fit.pdf"; 
    std::string plot_LowEnu = "QSq_LowEnu_ROOT_Fit.pdf"; 

    //Fit(data_rootDir_HighEnu, mc_rootDir_HighEnu, plot_HighEnu);
    //Fit(data_rootDir_LowEnu, mc_rootDir_LowEnu, plot_LowEnu);
    
    Fit(data_Trung_HighEnu, data_Trung_HighEnu, plot_HighEnu);
    Fit(data_Trung_LowEnu, data_Trung_LowEnu, plot_LowEnu);
}
