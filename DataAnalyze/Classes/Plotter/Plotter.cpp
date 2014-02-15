/*
    See Plotter.h header for Class Information
*/
#include "Plotter.h"

using namespace std;

Plotter::Plotter()
{

}

void Plotter::plotHistograms(string mcFile, string plotDir)
{

    TFile* f_mc = new TFile(mcFile.c_str());

    //------------------------------------------------------------------------
    //  2D Histograms
    //------------------------------------------------------------------------

/*    
    TH2F* h_reco_proton_energy_true_proton_energy = f_mc->Get("reco_proton_energy_true_proton_energy");
    plot2DHist(h_reco_proton_energy_true_proton_energy,"reco_proton_energy_true_proton_energy.png",plotDir);    
    
*/

    //------------------------------------------------------------------------
    // 1D Histograms
    //------------------------------------------------------------------------

    TH1F* h_beamEnergy= f_mc->Get("beamEnergy");
    plotSingleHist_1D(h_beamEnergy,"beamEnergy.png",plotDir);
     
    TH1F* h_P_muon= f_mc->Get("P_muon");
    plotSingleHist_1D(h_P_muon,"P_muon.png",plotDir);
    
    TH1F* h_P_proton= f_mc->Get("P_proton");
    plotSingleHist_1D(h_P_proton,"P_proton.png",plotDir);
    
    TH1F* h_P_pion= f_mc->Get("P_pion");
    plotSingleHist_1D(h_P_pion,"P_pion.png",plotDir);
    
    // -------------------------------------------------------------------------
    //     Angles
    //-------------------------------------------------------------------------- 
    TH1F* h_Angle_muon= f_mc->Get("Angle_muon");
    plotSingleHist_1D(h_Angle_muon,"Angle_muon.png",plotDir);
    
    TH1F* h_Angle_proton= f_mc->Get("Angle_proton");
    plotSingleHist_1D(h_Angle_proton,"Angle_proton.png",plotDir);
    
    TH1F* h_Angle_pion= f_mc->Get("Angle_pion");
    plotSingleHist_1D(h_Angle_pion,"Angle_pion.png",plotDir);
    
    TH1F* h_AngleMuon_muon= f_mc->Get("AngleMuon_muon");
    plotSingleHist_1D(h_AngleMuon_muon,"AngleMuon_muon.png",plotDir);
    
    TH1F* h_AngleMuon_proton= f_mc->Get("AngleMuon_proton");
    plotSingleHist_1D(h_AngleMuon_proton,"AngleMuon_proton.png",plotDir);
    
    TH1F* h_AngleMuon_pion= f_mc->Get("AngleMuon_pion");
    plotSingleHist_1D(h_AngleMuon_pion,"AngleMuon_pion.png",plotDir);


}

void Plotter::plotSingleHist_1D(TH1F* singleHist, string fileName, string plotDir)
{
    TCanvas* c1 = new TCanvas();
    singleHist->SetLineColor(kRed);
    singleHist->SetLineWidth(3);
    singleHist->SetFillColor(kRed);
    singleHist->SetFillStyle(3010);
    singleHist->Draw();
    c1->Print(Form("%s%s",plotDir.c_str(),fileName.c_str()),"png");
    
}


void Plotter::plot2DHist(TH2F* hist2D, string fileName, string plotDir)
{
    TCanvas* c1 = new TCanvas();
    hist2D->Draw("colz");
    c1->Print(Form("%s%s",plotDir.c_str(),fileName.c_str()),"png");
}




