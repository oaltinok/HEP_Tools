#include "Plotter.h"

using namespace std;

Plotter::Plotter()
{

}

void Plotter::plotHistograms(char* mcFile, string plotDir)
{

    TFile* f_mc = new TFile(mcFile);

    //------------------------------------------------------------------------
    //  2D Histograms
    //------------------------------------------------------------------------

     // Sample
    TH2F* h_mc_pion_score_proton_mom = f_mc->Get("pion_score_proton_mom");
    plot2DHist(h_mc_pion_score_proton_mom,"pion_score_proton_mom.png",plotDir);

    //------------------------------------------------------------------------
    // 1D Histograms
    //------------------------------------------------------------------------

    // Measured dEdX
//     TH1F* h_mc_measured_pion_dEdX= f_mc->Get("measured_pion_dEdX");
//     plotSingleHist_1D(h_mc_measured_pion_dEdX,"measured_pion_dEdX.png",plotDir);
//     
//     TH1F* h_mc_dEdX_difference= f_mc->Get("dEdX_difference");
//     plotSingleHist_1D(h_mc_dEdX_difference,"dEdX_difference.png",plotDir);

       // pID Scores
    TH1F* h_mc_pion_score = f_mc->Get("pion_score");
    plotSingleHist_1D(h_mc_pion_score,"pion_score.png",plotDir);
    
    TH1F* h_mc_proton_score = f_mc->Get("proton_score");
    plotSingleHist_1D(h_mc_proton_score,"proton_score.png",plotDir);
    
//     TH1F* h_mc_proton_score = f_mc->Get("proton_score");
//     plotSingleHist_1D(h_mc_proton_score,"proton_score.png",plotDir);
}

void Plotter::plotSingleHist_1D(TH1F* singleHist, string fileName, string plotDir)
{
    TCanvas* c1 = new TCanvas();
    singleHist->SetLineColor(kRed);
    singleHist->SetLineWidth(3);
    singleHist->SetFillColor(kRed);
    singleHist->SetFillStyle(3010);
    singleHist->Draw();
    c1->Print(Form("%s/%s",plotDir.c_str(),fileName.c_str()),"png");
}

void Plotter::plot2DHist(TH2F* hist2D, string fileName, string plotDir)
{
    TCanvas* c1 = new TCanvas();
    hist2D->Draw("colz");
    c1->Print(Form("%s/%s",plotDir.c_str(),fileName.c_str()),"png");
}





