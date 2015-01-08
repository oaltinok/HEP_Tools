#ifndef PionFunctions_cpp
#define PionFunctions_cpp

#include "CCProtonPi0.h"

using namespace std;

void CCProtonPi0::fillPionReco()
{
    double photon_E_asym;
    
    // Fill 4-Momentum
    pion.set_p4(    pi0_px,
                    pi0_py,
                    pi0_pz,
                    pi0_E,
                    false);
    
    
    // set Angle wrt Beam
    pion.set_angleBeam(beam_p3, false);
    
    // set Angle wrt Muon
    pion.set_angleMuon(muon, false);
    
    // Set Invariant Mass
    pion.invMass->Fill(pi0_invMass);
    
    // Set Photon Conversion Length in [cm]
    pion.gamma1_ConvLength->Fill(gamma1_dist_vtx);
    pion.gamma2_ConvLength->Fill(gamma2_dist_vtx);
    pion.ConvLength_gamma2_gamma1->Fill(gamma2_dist_vtx, gamma1_dist_vtx);
    
    // Set Photon N(Clusters)
    pion.gamma1_nClusters_All->Fill(gamma1_blob_nclusters);
    pion.gamma2_nClusters_All->Fill(gamma2_blob_nclusters);
    pion.nClusters_All_gamma2_gamma1->Fill(gamma2_blob_nclusters,gamma1_blob_nclusters);
    
    int gamma1_xCLusters;
    int gamma2_xCLusters;
    if(final_blob_nc[0] >= final_blob_nc[1]){
        gamma1_xCLusters = final_blob_nc[0];
        gamma2_xCLusters = final_blob_nc[1];
    }else{ 
        gamma1_xCLusters = final_blob_nc[1];
        gamma2_xCLusters = final_blob_nc[0];
    }
    pion.gamma1_nClusters_X->Fill(gamma1_xCLusters);
    pion.gamma2_nClusters_X->Fill(gamma2_xCLusters);
    pion.nClusters_X_gamma2_gamma1->Fill(gamma2_xCLusters,gamma1_xCLusters);
    
    // Set Photon Energy [GeV]
    pion.gamma1_Energy->Fill(gamma1_E * HEP_Functions::MeV_to_GeV);
    pion.gamma2_Energy->Fill(gamma2_E * HEP_Functions::MeV_to_GeV);
    pion.Energy_gamma2_gamma1->Fill(gamma2_E * HEP_Functions::MeV_to_GeV,gamma1_E * HEP_Functions::MeV_to_GeV );
    
    // Set Photon Energy Assymmetry
    photon_E_asym = abs((gamma1_E - gamma2_E) / (gamma1_E + gamma2_E));  
    pion.photonEnergy_Asymmetry->Fill(photon_E_asym);
}

void CCProtonPi0::fillPionTrue()
{
    double P_true;
    double P_reco;
    
    // Fill 4-Momentum
    if (truth_pi0_E != -1){
        pion.set_p4(truth_pi0_px,
                    truth_pi0_py,
                    truth_pi0_pz,
                    truth_pi0_E,
                    true);
    }else{
        pion.set_p4(truth_gamma_px[0]+truth_gamma_px[1] ,
                    truth_gamma_py[0]+truth_gamma_py[1],
                    truth_gamma_pz[0]+truth_gamma_pz[1],
                    truth_gamma_E[0]+truth_gamma_E[1],
                    true);
    }
    
    P_reco = HEP_Functions::calcMomentum(pi0_px,pi0_py,pi0_pz);
   
    // Momentum and Invariant Mass Information with True Pi0 Count
    if(truth_N_pi0 == 0){
        pion.P_reco_0Pi0->Fill(P_reco);
        pion.invMass_0Pi0->Fill(pi0_invMass);
    }else if(truth_N_pi0 == 1){
        pion.P_reco_1Pi0->Fill(P_reco);
        pion.invMass_1Pi0->Fill(pi0_invMass);
        P_true = getBestPi0Momentum();
        pion.P_reco_mc_1Pi0->Fill(P_reco,P_true);
    }else if(truth_N_pi0 > 1){
        pion.P_reco_MultPi0->Fill(P_reco);
        pion.invMass_MultPi0->Fill(pi0_invMass);
        P_true = getBestPi0Momentum();
        pion.P_reco_mc_MultPi0->Fill(P_reco,P_true);
    }
}

// Loops over all FS Particles and returns the most energetic pi0
double CCProtonPi0::getBestPi0Momentum()
{
    TVector3 p3;
    double tempP = 0;
    
    for(int i = 0; i < mc_nFSPart && i < max_nFSPart; i++ ){
        if( mc_FSPartPDG[i] == 111){
            p3.SetXYZ(mc_FSPartPx[i],mc_FSPartPy[i],mc_FSPartPz[i]);
            if(p3.Mag() > tempP) tempP = p3.Mag();
        }
    }
    return tempP;
}

int CCProtonPi0::getBestPi0()
{
    TVector3 p3;
    double tempP = 0;
    double ind = -1;
    
    for(int i = 0; i < mc_nFSPart && i < max_nFSPart; i++ ){
        if( mc_FSPartPDG[i] == 111){
            p3.SetXYZ(mc_FSPartPx[i],mc_FSPartPy[i],mc_FSPartPz[i]);
            if(p3.Mag() > tempP){
                tempP = p3.Mag();
                ind = i;
            }
        }
    }
    return ind;
}

bool CCProtonPi0::isPhotonDistanceLow()
{
    
    if (gamma1_dist_vtx < minPhotonDistance && gamma2_dist_vtx < minPhotonDistance){
        return true;
    }else{
        return false;
    }
}

#endif
