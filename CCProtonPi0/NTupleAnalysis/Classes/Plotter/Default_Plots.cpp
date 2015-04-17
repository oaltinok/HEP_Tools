#ifndef Single_Plots_cpp
#define Single_Plots_cpp

#include "Plotter.h"

using namespace std;

void Plotter::plotInteraction(bool isMC, bool isReco, bool is2D)
{

    string rootDir = rootDir_Interaction[branchInd];
    string plotDir = plotDir_Interaction[branchInd];
    
    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());

    // Signal only Plots
    if( branchInd == 0){
        TH1D* h_status_Pi0 = (TH1D*)f_Root->Get("status_Pi0");
        plot1D_Hist(h_status_Pi0,"status_Pi0.png",plotDir);
        
        TH1D* h_status_Pi0_Mother = (TH1D*)f_Root->Get("status_Pi0_Mother");
        plot1D_Hist(h_status_Pi0_Mother,"status_Pi0_Mother.png",plotDir);
        
        TH1D* h_status_Pi0_GrandMother = (TH1D*)f_Root->Get("status_Pi0_GrandMother");
        plot1D_Hist(h_status_Pi0_GrandMother,"status_Pi0_GrandMother.png",plotDir);
        
    }
    
    // Plot Only MC Values
    if( isMC ){
            
        TH1D* h_beamEnergy_mc= (TH1D*)f_Root->Get("beamEnergy_mc");
        plot1D_Hist(h_beamEnergy_mc,"beamEnergy_mc.png",plotDir);
        
        TH1D* h_q2_mc= (TH1D*)f_Root->Get("q2_mc");
        plot1D_Hist(h_q2_mc,"q2_mc.png",plotDir);
        
        TH1D* h_w_mc= (TH1D*)f_Root->Get("w_mc");
        plot1D_Hist(h_w_mc,"w_mc.png",plotDir);
        
        TH1D* h_int_channel= (TH1D*)f_Root->Get("int_channel");
        plot1D_Hist(h_int_channel,"int_channel.png",plotDir);
        
        TH1D* h_vertex_z_true= (TH1D*)f_Root->Get("vertex_z_true");
        plot1D_Hist(h_vertex_z_true,"vertex_z_true.png",plotDir);
        
        TH2D* h_vertex_x_y_true= (TH2D*)f_Root->Get("vertex_x_y_true");
        plot2D_Hist(h_vertex_x_y_true,"vertex_x_y_true.png",plotDir);
    }
    
    // Plot Only Reco Values
    if( isReco ){
        TH1D* h_total_E = (TH1D*)f_Root->Get("total_E");
        plot1D_Hist(h_total_E,"total_E.png",plotDir);
        
        TH1D* h_nProngs_hist = (TH1D*)f_Root->Get("nProngs_hist");
        plot1D_Hist(h_nProngs_hist,"nProngs_hist.png",plotDir);
        
        TH1D* h_deltaInvMass_reco = (TH1D*)f_Root->Get("deltaInvMass_reco");
        plot1D_Hist(h_deltaInvMass_reco,"deltaInvMass_reco.png",plotDir);
        
        TH1D* h_vertex_z_reco= (TH1D*)f_Root->Get("vertex_z_reco");
        plot1D_Hist(h_vertex_z_reco,"vertex_z_reco.png",plotDir);
        
        TH2D* h_vertex_x_y_reco= (TH2D*)f_Root->Get("vertex_x_y_reco");
        plot2D_Hist(h_vertex_x_y_reco,"vertex_x_y_reco.png",plotDir);
        
        TH1D* h_pFilter_Status = (TH1D*)f_Root->Get("pFilter_Status");
        plot1D_Hist(h_pFilter_Status,"pFilter_Status.png",plotDir);
        
        TH1D* h_pFilter_RejectedEnergy = (TH1D*)f_Root->Get("pFilter_RejectedEnergy");
        plot1D_Hist(h_pFilter_RejectedEnergy,"pFilter_RejectedEnergy.png",plotDir);
        
        TH1D* h_beamEnergy_reco= (TH1D*)f_Root->Get("beamEnergy_reco");
        plot1D_Hist(h_beamEnergy_reco,"beamEnergy_reco.png",plotDir);
        
        TH1D* h_beamEnergyCal_reco= (TH1D*)f_Root->Get("beamEnergyCal_reco");
        plot1D_Hist(h_beamEnergyCal_reco,"beamEnergyCal_reco.png",plotDir);
        
        TH1D* h_q2_reco= (TH1D*)f_Root->Get("q2_reco");
        plot1D_Hist(h_q2_reco,"q2_reco.png",plotDir);    
        
        TH1D* h_w_reco= (TH1D*)f_Root->Get("w_reco");
        plot1D_Hist(h_w_reco,"w_reco.png",plotDir); 
        
        TH1D* h_wSq_reco= (TH1D*)f_Root->Get("wSq_reco");
        plot1D_Hist(h_wSq_reco,"wSq_reco.png",plotDir); 

    }

    //  Plot 2D Comparison Plots and Error Plots
    if ( is2D ){
        TH2D* h_total_E_neutrinoE = (TH2D*)f_Root->Get("total_E_neutrinoE");
        plot2D_Hist(h_total_E_neutrinoE,"total_E_neutrinoE.png",plotDir);
        
        TH2D* h_deltaInvMass_reco_mc = (TH2D*)f_Root->Get("deltaInvMass_reco_mc");
        plot2D_Hist(h_deltaInvMass_reco_mc,"deltaInvMass_reco_mc.png",plotDir);
        
        TH1D* h_deltaInvMass_error = (TH1D*)f_Root->Get("deltaInvMass_error");
        plot1D_Hist(h_deltaInvMass_error,"deltaInvMass_error.png",plotDir);
         
        TH2D* h_beamEnergy_reco_mc= (TH2D*)f_Root->Get("beamEnergy_reco_mc");
        plot2D_Hist(h_beamEnergy_reco_mc,"beamEnergy_reco_mc.png",plotDir);

        TH1D* h_beamEnergy_error= (TH1D*)f_Root->Get("beamEnergy_error");
        plot1D_Hist(h_beamEnergy_error,"beamEnergy_error.png",plotDir);
        
        TH2D* h_beamEnergyCal_reco_mc= (TH2D*)f_Root->Get("beamEnergyCal_reco_mc");
        plot2D_Hist(h_beamEnergyCal_reco_mc,"beamEnergyCal_reco_mc.png",plotDir);
        
        TH1D* h_beamEnergyCal_error= (TH1D*)f_Root->Get("beamEnergyCal_error");
        plot1D_Hist(h_beamEnergyCal_error,"beamEnergyCal_error.png",plotDir);
        
        TH2D* h_beamEnergy_beamEnergyCal= (TH2D*)f_Root->Get("beamEnergy_beamEnergyCal");
        plot2D_Hist(h_beamEnergy_beamEnergyCal,"beamEnergy_beamEnergyCal.png",plotDir);
        
        TH2D* h_q2_reco_mc= (TH2D*)f_Root->Get("q2_reco_mc");
        plot2D_Hist(h_q2_reco_mc,"q2_reco_mc.png",plotDir);
        
        TH1D* h_q2_error= (TH1D*)f_Root->Get("q2_error");
        plot1D_Hist(h_q2_error,"q2_error.png",plotDir);
        
        TH2D* h_w_reco_mc= (TH2D*)f_Root->Get("w_reco_mc");
        plot2D_Hist(h_w_reco_mc,"w_reco_mc.png",plotDir);
        
        TH1D* h_w_error= (TH1D*)f_Root->Get("w_error");
        plot1D_Hist(h_w_error,"w_error.png",plotDir);

        TH2D* h_vertex_z_reco_mc= (TH2D*)f_Root->Get("vertex_z_reco_mc");
        plot2D_Hist(h_vertex_z_reco_mc,"vertex_z_reco_mc.png",plotDir);
        
        TH1D* h_vertex_z_error= (TH1D*)f_Root->Get("vertex_z_error");
        plot1D_Hist(h_vertex_z_error,"vertex_z_error.png",plotDir);
    }
    
    delete f_Root;
    
}


void Plotter::plotParticleInfo(  string rootDir, string plotDir, 
                                bool isMC, bool isReco, bool is2D)
{

    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());

    // Plot Only MC Values
    if( isMC ){
        TH1D* h_E_mc = (TH1D*)f_Root->Get("E_mc");
        plot1D_Hist(h_E_mc,"P_mc.png",plotDir);
        
        TH1D* h_P_mc = (TH1D*)f_Root->Get("P_mc");
        plot1D_Hist(h_P_mc,"P_mc.png",plotDir);
        
        TH1D* h_KE_mc = (TH1D*)f_Root->Get("KE_mc");
        plot1D_Hist(h_KE_mc,"KE_mc.png",plotDir);
        
        TH1D* h_angleBeam_mc = (TH1D*)f_Root->Get("angleBeam_mc");
        plot1D_Hist(h_angleBeam_mc,"angleBeam_mc.png",plotDir);
        
        TH1D* h_angleMuon_mc = (TH1D*)f_Root->Get("angleMuon_mc");
        plot1D_Hist(h_angleMuon_mc,"angleMuon_mc.png",plotDir);
    }
    
    // Plot Only Reco Values
    if( isReco ){
        TH1D* h_partScore = (TH1D*)f_Root->Get("partScore");
        plot1D_Hist(h_partScore,"partScore.png",plotDir);
        
        TH1D* h_E_reco = (TH1D*)f_Root->Get("E_reco");
        plot1D_Hist(h_E_reco,"E_reco.png",plotDir);
        plot1D_HistLogScale(h_E_reco,"E_reco_log.png",plotDir);
        
        TH1D* h_P_reco = (TH1D*)f_Root->Get("P_reco");
        plot1D_Hist(h_P_reco,"P_reco.png",plotDir);
        plot1D_HistLogScale(h_P_reco,"P_reco_log.png",plotDir);
        
        TH1D* h_KE_reco = (TH1D*)f_Root->Get("KE_reco");
        plot1D_Hist(h_KE_reco,"KE_reco.png",plotDir);    
        
        TH1D* h_angleBeam_reco = (TH1D*)f_Root->Get("angleBeam_reco");
        plot1D_Hist(h_angleBeam_reco,"angleBeam_reco.png",plotDir);
        
        TH1D* h_angleMuon_reco = (TH1D*)f_Root->Get("angleMuon_reco");
        plot1D_Hist(h_angleMuon_reco,"angleMuon_reco.png",plotDir);
        
    }

    //  Plot 2D Comparison Plots and Error Plots
    if ( is2D ){
        TH2D* h_E_reco_mc = (TH2D*)f_Root->Get("E_reco_mc");
        plot2D_Hist(h_E_reco_mc,"E_reco_mc.png",plotDir);
        
        TH2D* h_P_reco_mc = (TH2D*)f_Root->Get("P_reco_mc");
        plot2D_Hist(h_P_reco_mc,"P_reco_mc.png",plotDir);
        
        TH1D* h_E_error = (TH1D*)f_Root->Get("E_error");
        plot1D_Hist(h_E_error,"E_error.png",plotDir);

        TH1D* h_P_error = (TH1D*)f_Root->Get("P_error");
        plot1D_Hist(h_P_error,"P_error.png",plotDir);
        
        TH2D* h_KE_reco_mc = (TH2D*)f_Root->Get("KE_reco_mc");
        plot2D_Hist(h_KE_reco_mc,"KE_reco_mc.png",plotDir);
        
        TH1D* h_KE_error = (TH1D*)f_Root->Get("KE_error");
        plot1D_Hist(h_KE_error,"KE_error.png",plotDir);
        
        TH2D* h_angleBeam_reco_mc = (TH2D*)f_Root->Get("angleBeam_reco_mc");
        plot2D_Hist(h_angleBeam_reco_mc,"angleBeam_reco_mc.png",plotDir);
        
        TH1D* h_angleBeam_error = (TH1D*)f_Root->Get("angleBeam_error");
        plot1D_Hist(h_angleBeam_error,"angleBeam_error.png",plotDir);
        
        TH2D* h_angleMuon_reco_mc = (TH2D*)f_Root->Get("angleMuon_reco_mc");
        plot2D_Hist(h_angleMuon_reco_mc,"angleMuon_reco_mc.png",plotDir);
        
        TH1D* h_angleMuon_error = (TH1D*)f_Root->Get("angleMuon_error");
        plot1D_Hist(h_angleMuon_error,"angleMuon_error.png",plotDir);
    }

    delete f_Root;
}

void Plotter::plotMuon(bool isMC, bool isReco, bool is2D)
{
    plotParticleInfo(rootDir_Muon[branchInd], plotDir_Muon[branchInd], isMC, isReco, is2D);
}

void Plotter::plotProton(bool isMC, bool isReco, bool is2D)
{    
    plotParticleInfo(rootDir_Proton[branchInd], plotDir_Proton[branchInd], isMC, isReco, is2D);
}

void Plotter::plotPion(bool isMC, bool isReco, bool is2D)
{
    string rootDir = rootDir_Pion[branchInd];
    string plotDir = plotDir_Pion[branchInd];
    
    // Standard Plots
    plotParticleInfo(rootDir, plotDir, isMC, isReco, is2D);
    
    // Unique Plots
    TFile* f_Root = new TFile(rootDir.c_str());
    
    TH1D* h_invMass = (TH1D*)f_Root->Get("invMass");
    plot1D_Hist(h_invMass,"invMass.png",plotDir);
    
    // Photon Conversion Length
    TH1D* h_gamma1_ConvLength = (TH1D*)f_Root->Get("gamma1_ConvLength");
    plot1D_Hist(h_gamma1_ConvLength,"gamma1_ConvLength.png",plotDir);
    
    TH1D* h_gamma2_ConvLength = (TH1D*)f_Root->Get("gamma2_ConvLength");
    plot1D_Hist(h_gamma2_ConvLength,"gamma2_ConvLength.png",plotDir);
    
    TH2D* h_ConvLength_gamma2_gamma1 = (TH2D*)f_Root->Get("ConvLength_gamma2_gamma1");
    plot2D_Hist(h_ConvLength_gamma2_gamma1,"ConvLength_gamma2_gamma1.png",plotDir);
    
    // nClusters
    TH2D* h_nClusters_All_gamma2_gamma1 = (TH2D*)f_Root->Get("nClusters_All_gamma2_gamma1");
    plot2D_Hist(h_nClusters_All_gamma2_gamma1,"nClusters_All_gamma2_gamma1.png",plotDir);
    
    // Gamma Energy
    TH2D* h_Energy_gamma2_gamma1 = (TH2D*)f_Root->Get("Energy_gamma2_gamma1");
    plot2D_Hist(h_Energy_gamma2_gamma1,"Energy_gamma2_gamma1.png",plotDir);
    
    
    
    delete f_Root;
}


#endif

