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
#include <TF1.h>
#include <TPaveText.h>
#include <TVirtualFitter.h>
#include <TStyle.h>

#include <PlotUtils/MnvH1D.h>

using PlotUtils::MnvH1D;

double user_expo(double* x, double* par)
{
    return par[0]*exp(-par[1]*x[0]);
}


void Fit(std::string fileName, std::string plotName)
{
    gStyle->SetErrorX(0.5);

    std::string histName = "QSq_xsec";

    TFile* f = new TFile(fileName.c_str(), "read");
    MnvH1D* mnvh1d = static_cast<MnvH1D*>(f->Get(histName.c_str()));

    assert(mnvh1d);

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
    
    TF1* f1_expo = new TF1("user_expo", user_expo,0.0,2.0,2);
    h1d->Fit(f1_expo, "0", "", 0.25, 2.0);
    
    f1_expo->SetLineColor(kOrange-3);
    f1_expo->Draw("same");
    
    TPaveText* expo_label = new TPaveText(0.6, 0.6, 0.8, 0.8, "NDC");
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
   
    std::string output_folder = "/minerva/app/users/oaltinok/cmtuser/Minerva_v10r8p9/Ana/CCProtonPi0/NTupleAnalysis/Output/Plots/OtherStudies/";
    std::string output = output_folder + plotName;    
 
    // Remove Stat Box
    gStyle->SetOptStat(0); 
    
    // Remove Title 
    gStyle->SetOptTitle(0); 


    c1->Print(output.c_str());
}


void Fit_QSq()
{
    std::string rootDir_HighEnu = "/minerva/data/users/oaltinok/NTupleAnalysis_QSq_HighEnu_4/Data/Analyzed/CrossSection.root";
    std::string rootDir_LowEnu  = "/minerva/data/users/oaltinok/NTupleAnalysis_QSq_LowEnu_4/Data/Analyzed/CrossSection.root";

    std::string plot_HighEnu = "QSq_HighEnu_ROOT_Fit.png"; 
    std::string plot_LowEnu = "QSq_LowEnu_ROOT_Fit.png"; 

    Fit(rootDir_HighEnu, plot_HighEnu);
    Fit(rootDir_LowEnu, plot_LowEnu);
}
