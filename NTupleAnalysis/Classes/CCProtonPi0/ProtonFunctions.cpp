#ifndef ProtonFunctions_cpp
#define ProtonFunctions_cpp

#include "CCProtonPi0.h"

using namespace std;

/*
--------------------------------------------------------------------------------
 countParticles:
    Returns the number of particles in the Final State
        Input 
            int targetPDG
            bool applyPCut - Variable for selecting particles with momentum 
                                (no particle at rest)
--------------------------------------------------------------------------------
*/
int CCProtonPi0::countParticles(int targetPDG, bool applyPCut)
{
    int count = 0;
    TVector3 p3;
    
    for(int i = 0; i < mc_nFSPart && i < max_nFSPart; i++ ){
        if( mc_FSPartPDG[i] == targetPDG){
            if(applyPCut){
                p3.SetXYZ(mc_FSPartPx[i],mc_FSPartPy[i],mc_FSPartPz[i]);
                if(p3.Mag() > 0){
                    count++;
                }
            }
            else{
                count++;
            }
        }
    }
    
    return count;

}

void CCProtonPi0::fillProtonTrue(int ind)
{    
    // Fill 4-Momentum
    proton.set_p4(  CCProtonPi0_trajProtonProngPx[ind],
                    CCProtonPi0_trajProtonProngPy[ind],
                    CCProtonPi0_trajProtonProngPz[ind],
                    -1.0, 
                    true);
       
    // set Angle wrt Beam
    proton.set_angleBeam(beam_p3, true);
    
    // set Angle wrt Muon
    proton.set_angleMuon(muon, true);
    
}

void CCProtonPi0::fillProtonReco(int ind)
{
    // Set Particle Score
    proton.particleScore = CCProtonPi0_proton_score[ind];
    
    // Fill 4-Momentum
    proton.set_p4(  CCProtonPi0_proton_px[ind],
                    CCProtonPi0_proton_py[ind],
                    CCProtonPi0_proton_pz[ind],
                    CCProtonPi0_proton_E[ind],
                    false);
                    
    // set Angle wrt Beam
    proton.set_angleBeam(beam_p3, false);
    
    // set Angle wrt Muon
    proton.set_angleMuon(muon, false);
}

int CCProtonPi0::findBestProton()
{
    double tempScore = CCProtonPi0_proton_score[0];
    int tempInd = 0;
    
    for( int i = 0; i < 10; i++){
        if( CCProtonPi0_proton_score[i] == -1 ) break;
        if( CCProtonPi0_proton_score[i] > tempScore){
            tempScore = CCProtonPi0_proton_score[i];
            tempInd = i;
        }
    }
    
    return tempInd;

}

#endif