#ifndef PionFunctions_cpp
#define PionFunctions_cpp

#include "CCProtonPi0.h"

using namespace std;

void CCProtonPi0::fillPionReco()
{
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
}

void CCProtonPi0::fillPionTrue()
{
    double P_true;
    double P_reco;
    
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
