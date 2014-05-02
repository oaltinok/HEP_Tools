/*
================================================================================
Function List: Muon_Functions.cpp
    Contains Muon Specific Functions
    
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


void CCDeltaPlus::fillMuonTrue()
{

    // Fill 4-Momentum
    muon.set_p4(    CCDeltaPlusAna_trajMuonProngPx,
                    CCDeltaPlusAna_trajMuonProngPy,
                    CCDeltaPlusAna_trajMuonProngPz,
                    -1.0, 
                    true);
       
    // set Angle wrt Beam
    muon.set_angleBeam(beam_p3, true);
    
    // set Angle wrt Muon
    muon.set_angleMuon(muon, true);
    
}

void CCDeltaPlus::fillMuonReco()
{
    // Set Particle Score
    muon.particleScore = CCDeltaPlusAna_muon_muScore;
    
    // Fill 4-Momentum
    muon.set_p4(    CCDeltaPlusAna_muon_px,
                    CCDeltaPlusAna_muon_py,
                    CCDeltaPlusAna_muon_pz,
                    CCDeltaPlusAna_muon_E,
                    false);
    
    // set Angle wrt Beam
    muon.set_angleBeam(beam_p3, false);
    
    // set Angle wrt Muon
    muon.set_angleMuon(muon, false);
}

#endif // #ifdef CCDeltaPlus_cxx