/*
    See Plotter.h header for Class Information
*/

#ifndef PLOTTER_CPP
#define PLOTTER_CPP
#include "Plotter.h"

using namespace std;


void Plotter::plotHistograms(bool isMC, bool isReco, bool is2D)
{
    
//     plotCCProtonPi0(isMC,isReco, is2D);

    plotMuon(isMC,isReco, is2D);
    
//     plotProton(isMC,isReco, is2D);
    
//     plotPion(isMC,isReco, is2D);
    
}


void Plotter::plotCCProtonPi0(bool isMC, bool isReco, bool is2D)
{

    string rootDir = "Output/RootFiles/CCProtonPi0.root";
    string plotDir = "Output/Plots/CCProtonPi0/";
    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());
    
    // Plot Only MC Values
    if( isMC ){
        TH1F* h_beamEnergy_mc= (TH1F*)f_Root->Get("beamEnergy_mc");
        plot1D_Hist(h_beamEnergy_mc,"beamEnergy_mc.png",plotDir);
        
        TH1F* h_q2_mc= (TH1F*)f_Root->Get("q2_mc");
        plot1D_Hist(h_q2_mc,"q2_mc.png",plotDir);
        
        TH1F* h_int_channel= (TH1F*)f_Root->Get("int_channel");
        plot1D_Hist(h_int_channel,"int_channel.png",plotDir);
        
        TH1F* h_vertex_z_true= (TH1F*)f_Root->Get("vertex_z_true");
        plot1D_Hist(h_vertex_z_true,"vertex_z_true.png",plotDir);
        
        TH2F* h_vertex_x_y_true= (TH2F*)f_Root->Get("vertex_x_y_true");
        plot2D_Hist(h_vertex_x_y_true,"vertex_x_y_true.png",plotDir);
        
        TH1F* h_n_FSParticles= (TH1F*)f_Root->Get("n_FSParticles");
        plot1D_Hist(h_n_FSParticles,"n_FSParticles.png",plotDir);
        
        TH1F* h_n_gammas= (TH1F*)f_Root->Get("n_gammas");
        plot1D_Hist(h_n_gammas,"n_gammas.png",plotDir);        
    }
    
    // Plot Only Reco Values
    if( isReco ){
        TH1F* h_vertex_z_reco= (TH1F*)f_Root->Get("vertex_z_reco");
        plot1D_Hist(h_vertex_z_reco,"vertex_z_reco.png",plotDir);
        
        TH2F* h_vertex_x_y_reco= (TH2F*)f_Root->Get("vertex_x_y_reco");
        plot2D_Hist(h_vertex_x_y_reco,"vertex_x_y_reco.png",plotDir);
        
        TH1F* h_beamEnergy_reco= (TH1F*)f_Root->Get("beamEnergy_reco");
        plot1D_Hist(h_beamEnergy_reco,"beamEnergy_reco.png",plotDir);
        
        TH1F* h_q2_reco= (TH1F*)f_Root->Get("q2_reco");
        plot1D_Hist(h_q2_reco,"q2_reco.png",plotDir);    

    }

    //  Plot 2D Comparison Plots and Error Plots
    if ( is2D ){
        TH2F* h_beamEnergy_reco_mc= (TH2F*)f_Root->Get("beamEnergy_reco_mc");
        plot2D_Hist(h_beamEnergy_reco_mc,"beamEnergy_reco_mc.png",plotDir);

        TH1F* h_beamEnergy_error= (TH1F*)f_Root->Get("beamEnergy_error");
        plot1D_Hist(h_beamEnergy_error,"beamEnergy_error.png",plotDir);
        
        TH2F* h_q2_reco_mc= (TH2F*)f_Root->Get("q2_reco_mc");
        plot2D_Hist(h_q2_reco_mc,"q2_reco_mc.png",plotDir);
        
        TH1F* h_q2_error= (TH1F*)f_Root->Get("q2_error");
        plot1D_Hist(h_q2_error,"q2_error.png",plotDir);

        TH2F* h_vertex_z_reco_mc= (TH2F*)f_Root->Get("vertex_z_reco_mc");
        plot2D_Hist(h_vertex_z_reco_mc,"vertex_z_reco_mc.png",plotDir);
        
        TH1F* h_vertex_z_error= (TH1F*)f_Root->Get("vertex_z_error");
        plot1D_Hist(h_vertex_z_error,"vertex_z_error.png",plotDir);
    }
    
    
}


void Plotter::plotParticleInfo(  string rootDir, string plotDir, 
                                bool isMC, bool isReco, bool is2D)
{

    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());

    // Plot Only MC Values
    if( isMC ){
        TH1F* h_P_mc= (TH1F*)f_Root->Get("P_mc");
        plot1D_Hist(h_P_mc,"P_mc.png",plotDir);
        
        TH1F* h_KE_mc= (TH1F*)f_Root->Get("KE_mc");
        plot1D_Hist(h_KE_mc,"KE_mc.png",plotDir);
        
        TH1F* h_angleBeam_mc= (TH1F*)f_Root->Get("angleBeam_mc");
        plot1D_Hist(h_angleBeam_mc,"angleBeam_mc.png",plotDir);
        
        TH1F* h_angleMuon_mc= (TH1F*)f_Root->Get("angleMuon_mc");
        plot1D_Hist(h_angleMuon_mc,"angleMuon_mc.png",plotDir);
    }
    
    // Plot Only Reco Values
    if( isReco ){
        TH1F* h_partScore= (TH1F*)f_Root->Get("partScore");
        plot1D_Hist(h_partScore,"partScore.png",plotDir);
        
        TH1F* h_P_reco= (TH1F*)f_Root->Get("P_reco");
        plot1D_Hist(h_P_reco,"P_reco.png",plotDir);
        
        TH1F* h_KE_reco= (TH1F*)f_Root->Get("KE_reco");
        plot1D_Hist(h_KE_reco,"KE_reco.png",plotDir);    
        
        TH1F* h_angleBeam_reco= (TH1F*)f_Root->Get("angleBeam_reco");
        plot1D_Hist(h_angleBeam_reco,"angleBeam_reco.png",plotDir);
        
        TH1F* h_angleMuon_reco= (TH1F*)f_Root->Get("angleMuon_reco");
        plot1D_Hist(h_angleMuon_reco,"angleMuon_reco.png",plotDir);

    }

    //  Plot 2D Comparison Plots and Error Plots
    if ( is2D ){
        TH2F* h_P_reco_mc= (TH2F*)f_Root->Get("P_reco_mc");
        plot2D_Hist(h_P_reco_mc,"P_reco_mc.png",plotDir);

        TH1F* h_P_error= (TH1F*)f_Root->Get("P_error");
        plot1D_Hist(h_P_error,"P_error.png",plotDir);
        
        TH2F* h_KE_reco_mc= (TH2F*)f_Root->Get("KE_reco_mc");
        plot2D_Hist(h_KE_reco_mc,"KE_reco_mc.png",plotDir);
        
        TH1F* h_KE_error= (TH1F*)f_Root->Get("KE_error");
        plot1D_Hist(h_KE_error,"KE_error.png",plotDir);
        
        TH2F* h_angleBeam_reco_mc= (TH2F*)f_Root->Get("angleBeam_reco_mc");
        plot2D_Hist(h_angleBeam_reco_mc,"angleBeam_reco_mc.png",plotDir);
        
        TH1F* h_angleBeam_error= (TH1F*)f_Root->Get("angleBeam_error");
        plot1D_Hist(h_angleBeam_error,"angleBeam_error.png",plotDir);
        
        TH2F* h_angleMuon_reco_mc= (TH2F*)f_Root->Get("angleMuon_reco_mc");
        plot2D_Hist(h_angleMuon_reco_mc,"angleMuon_reco_mc.png",plotDir);
        
        TH1F* h_angleMuon_error= (TH1F*)f_Root->Get("angleMuon_error");
        plot1D_Hist(h_angleMuon_error,"angleMuon_error.png",plotDir);
    }

}

void Plotter::plotMuon(bool isMC, bool isReco, bool is2D)
{
    string rootDir = "Output/RootFiles/Muon.root";
    string plotDir = "Output/Plots/Muon/";
    
    plotParticleInfo(rootDir, plotDir, isMC, isReco, is2D);
}

void Plotter::plotProton(bool isMC, bool isReco, bool is2D)
{
    string rootDir = "Output/RootFiles/Proton.root";
    string plotDir = "Output/Plots/Proton/";
    
    plotParticleInfo(rootDir, plotDir, isMC, isReco, is2D);
}

void Plotter::plotPion(bool isMC, bool isReco, bool is2D)
{
    string rootDir = "Output/RootFiles/Pion.root";
    string plotDir = "Output/Plots/Pion/";
    
    plotParticleInfo(rootDir, plotDir, isMC, isReco, is2D);
}

void Plotter::inform(string rootDir, string plotDir)
{
    cout<<"------------ Plotting ------------"<<endl;
    cout<<"Input File: "<<rootDir<<endl;
    cout<<"Output Folder: "<<plotDir<<endl;

}

void Plotter::plot1D_Hist(TH1F* hist1D, string fileName, string plotDir)
{
    TCanvas* c1 = new TCanvas();
    hist1D->SetLineColor(kRed);
    hist1D->SetLineWidth(3);
    hist1D->SetFillColor(kRed);
    hist1D->SetFillStyle(3010);
    hist1D->Draw();
    c1->Print(Form("%s%s",plotDir.c_str(),fileName.c_str()),"png");
    
}

void Plotter::plot2D_Hist(TH2F* hist2D, string fileName, string plotDir)
{
    Double_t w = 600;
    Double_t h = 600;
    TCanvas* c1 = new TCanvas();
    c1->SetWindowSize(w,h);
    hist2D->Draw("colz");
    c1->Print(Form("%s%s",plotDir.c_str(),fileName.c_str()),"png");
}

#endif





