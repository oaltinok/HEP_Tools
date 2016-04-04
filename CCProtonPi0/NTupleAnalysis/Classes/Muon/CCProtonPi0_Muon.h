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

        MnvH2D* theta_theta_test;
        MnvH2D* thetaX_thetaX_test;
        MnvH2D* thetaY_thetaY_test;

        TH1D* thetaX_Diff;
        TH1D* thetaY_Diff;

        // Muon Momentum 
        MnvH1D* muon_P_all;
        MnvH1D* muon_P_mc_reco_all;
        MnvH1D* muon_P_mc_reco_signal;
        MnvH1D* muon_P_mc_reco_bckg;
        MnvH1D* muon_P_mc_truth_signal;
        MnvH2D* muon_P_response;
 
        // Muon Theta 
        MnvH1D* muon_theta_all;
        MnvH1D* muon_theta_mc_reco_all;
        MnvH1D* muon_theta_mc_reco_signal;
        MnvH1D* muon_theta_mc_reco_bckg;
        MnvH1D* muon_theta_mc_truth_signal;
        MnvH2D* muon_theta_response;

        // Muon Theta 
        MnvH1D* muon_cos_theta_all;
        MnvH1D* muon_cos_theta_mc_reco_all;
        MnvH1D* muon_cos_theta_mc_reco_signal;
        MnvH1D* muon_cos_theta_mc_reco_bckg;
        MnvH1D* muon_cos_theta_mc_truth_signal;
        MnvH2D* muon_cos_theta_response;

    private:
        CCProtonPi0_SingleBin bin_thetaX_Diff;
        CCProtonPi0_SingleBin bin_thetaY_Diff;
};

#endif
