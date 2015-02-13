#ifndef ProtonFunctions_cpp
#define ProtonFunctions_cpp

#include "CCProtonPi0.h"

using namespace std;


void CCProtonPi0::fillProtonTrue()
{    
    // Fill 4-Momentum
    proton.set_p4(  truth_proton_px[indTrueProton],
                    truth_proton_py[indTrueProton],
                    truth_proton_pz[indTrueProton],
                    truth_proton_E[indTrueProton], 
                    true);
       
    // set Angle wrt Beam
    proton.set_angleBeam(beam_p3, true);
    
    // set Angle wrt Muon
    proton.set_angleMuon(muon, true);
    
}

void CCProtonPi0::fillProtonReco()
{   
    // Set Particle Score
    proton.particleScore = CCProtonPi0_protonScore[indRecoProton];

    // Fill 4-Momentum
    proton.set_p4(  CCProtonPi0_proton_px[indRecoProton],
                    CCProtonPi0_proton_py[indRecoProton],
                    CCProtonPi0_proton_pz[indRecoProton],
                    CCProtonPi0_proton_E[indRecoProton],
                    false);
                    
    // set Angle wrt Beam
    proton.set_angleBeam(beam_p3, false);
    
    // set Angle wrt Muon
    proton.set_angleMuon(muon, false);
    
    // Set Particle Score - proton and pion
    proton.protonScore->Fill(CCProtonPi0_protonScore[indRecoProton]);
    proton.pionScore->Fill(CCProtonPi0_pionScore[indRecoProton]);
}

//--------------------------------------------------------------------------
//  findTrueProton()
//      finds Truth Proton using Proton Energy - Highest Energy Proton
//--------------------------------------------------------------------------
void CCProtonPi0::findTrueProton()
{
    double tempInd = 0;
    double tempE = truth_proton_E[0];
    
    for ( int i = 1; i < 20; i++){
        if (truth_proton_E[i] == SENTINEL) break;
        if (truth_proton_E[i] > tempE){
            tempE = truth_proton_E[i];
            tempInd = 0;
        }
    }
    
    indTrueProton = tempInd;
}

//--------------------------------------------------------------------------
//  findRecoProton()
//      finds Reco Proton using Proton Score
//--------------------------------------------------------------------------
void CCProtonPi0::findRecoProton()
{
    double tempScore = CCProtonPi0_protonScore[0];
    int tempInd = 0;
    
    for( int i = 0; i < 10; i++){
        if( CCProtonPi0_protonScore[i] == -1 ) break;
        if( CCProtonPi0_protonScore[i] > tempScore){
            tempScore = CCProtonPi0_protonScore[i];
            tempInd = i;
        }
    }
    
    indRecoProton = tempInd;

}

#endif