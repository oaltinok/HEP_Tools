#ifndef MuonFunctions_cpp
#define MuonFunctions_cpp

#include "CCProtonPi0.h"

using namespace std;

void CCProtonPi0::fillMuonTrue()
{

    // Fill 4-Momentum
    muon.set_p4(    truth_muon_px * HEP_Functions::MeV_to_GeV,
                    truth_muon_py * HEP_Functions::MeV_to_GeV,
                    truth_muon_pz * HEP_Functions::MeV_to_GeV,
                    truth_muon_E, 
                    true);
       
    // set Angle wrt Beam
    muon.set_angleBeam(beam_p3, true);
    
    // set Angle wrt Muon
    muon.set_angleMuon(muon, true);
    
}

void CCProtonPi0::fillMuonReco()
{
    // Set Particle Score
    muon.particleScore = CCProtonPi0_muon_muScore;
    
    // Fill 4-Momentum
    muon.set_p4(    CCProtonPi0_muon_px * HEP_Functions::MeV_to_GeV,
                    CCProtonPi0_muon_py * HEP_Functions::MeV_to_GeV,
                    CCProtonPi0_muon_pz * HEP_Functions::MeV_to_GeV,
                    CCProtonPi0_muon_E * HEP_Functions::MeV_to_GeV,
                    false);
    
    // set Angle wrt Beam
    muon.set_angleBeam(beam_p3, false);
    
    // set Angle wrt Muon
    muon.set_angleMuon(muon, false);
}


#endif
