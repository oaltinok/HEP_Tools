/*
================================================================================
Class: CCProtonPi0_Muon -> Derived Class from Particle Base Clas
    CCProtonPi0_Muon Class  inherits Particle Behaviours and 
                extends base class with Muon specific parameters
                
    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
================================================================================
*/

#ifndef CCProtonPi0_Muon_h
#define CCProtonPi0_Muon_h

#include "../Particle/CCProtonPi0_Particle.h"

class CCProtonPi0_Muon : public CCProtonPi0_Particle
{
    public:
        CCProtonPi0_Muon(bool isModeReduce, bool isMC);
        void initHistograms();
        void writeHistograms();

        // Muon Momentum 
        MnvH1D* data_all_muon_P;
        MnvH1D* mc_reco_bckg_muon_P;
        MnvH1D* mc_reco_signal_muon_P;
        MnvH1D* mc_truth_signal_muon_P;
        MnvH2D* response_P;
        
        // Muon Theta
        MnvH2D* response_theta;

    private:
        CCProtonPi0_SingleBin bin_muonTheta;
};

#endif
