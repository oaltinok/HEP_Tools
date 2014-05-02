/*
================================================================================
Function List: Proton_Functions.cpp
    Contains Proton Specific Functions
    
    Functions Defined under CCDeltaPlus Namespace

    Usage:
        > #include "CCDeltaPlus.cpp"
        > #ifdef CCDeltaPlus_cxx 
    
    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision:  2014_05_02
================================================================================
*/

#include "CCDeltaPlus.cpp"
#ifdef CCDeltaPlus_cxx

void CCDeltaPlus::fillProtonTrue()
{
    int ind = proton.ind;
    
    // Fill 4-Momentum
    proton.set_p4(  CCDeltaPlusAna_trajProtonProngPx[ind],
                    CCDeltaPlusAna_trajProtonProngPy[ind],
                    CCDeltaPlusAna_trajProtonProngPz[ind],
                    -1.0, 
                    true);
       
    // set Angle wrt Beam
    proton.set_angleBeam(beam_p3, true);
    
    // set Angle wrt Muon
    proton.set_angleMuon(muon, true);
    
}

void CCDeltaPlus::fillProtonReco(int ind)
{
    // Set Particle Score
    proton.particleScore = CCDeltaPlusAna_proton_score[ind];
    
    // Fill 4-Momentum
    proton.set_p4(  CCDeltaPlusAna_proton_px[ind],
                    CCDeltaPlusAna_proton_py[ind],
                    CCDeltaPlusAna_proton_pz[ind],
                    CCDeltaPlusAna_proton_E[ind],
                    false);
                    
    // set Angle wrt Beam
    proton.set_angleBeam(beam_p3, false);
    
    // set Angle wrt Muon
    proton.set_angleMuon(muon, false);
}

int CCDeltaPlus::findBestProton()
{
    double tempScore = CCDeltaPlusAna_proton_score[0];
    int tempInd = 0;
    
    for( int i = 1; i < 10; i++){
        if ( CCDeltaPlusAna_proton_score[i] == -1 ) break;
        if( CCDeltaPlusAna_proton_score[i] > tempScore){
            tempScore = CCDeltaPlusAna_proton_score[i];
            tempInd = i;
        }
    }
    
    return tempInd;

}

bool CCDeltaPlus::isProtonShort(int ind)
{
    const double protonMass = 938; 
    const double minProtonKE = 120;
    double protonKE;
    double protonE;
    
    protonE = mc_FSPartE[ind];
    protonKE =  protonE - protonMass;
    if( protonKE < minProtonKE){
        return true;
    }else{
        return false;
    }
    
}

#endif // #ifdef CCDeltaPlus_cxx