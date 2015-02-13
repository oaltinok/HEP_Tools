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
    pion.gamma1_ConvLength->Fill(gamma1_dist_vtx * 0.1);
    pion.gamma2_ConvLength->Fill(gamma2_dist_vtx * 0.1);
    pion.ConvLength_gamma2_gamma1->Fill(gamma2_dist_vtx, gamma1_dist_vtx);
        
    // Set Photon N(Clusters)
    pion.gamma1_nClusters_All->Fill(gamma1_blob_nclusters);
    pion.gamma2_nClusters_All->Fill(gamma2_blob_nclusters);
    pion.nClusters_All_gamma2_gamma1->Fill(gamma2_blob_nclusters,gamma1_blob_nclusters);
    
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
    
    // set Angle wrt Beam
    pion.set_angleBeam(beam_p3, true);
    
    // set Angle wrt Muon
    pion.set_angleMuon(muon, true);
}

#endif
