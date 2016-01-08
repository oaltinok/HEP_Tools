/*
================================================================================
Class: CCProtonPi0_Pion -> Derived Class from Particle Base Clas
    CCProtonPi0_Pion Class  inherits Particle Behaviours and 
                extends base class with pion specific parameters

    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
================================================================================
*/

#ifndef CCProtonPi0_Pion_h
#define CCProtonPi0_Pion_h

#include "../Particle/CCProtonPi0_Particle.h"

using namespace PlotUtils;

class CCProtonPi0_Pion : public CCProtonPi0_Particle
{
    public:
        // Leading Photon - Energetic Photon
        vector<MnvH1D*> gamma1_E;
        vector<MnvH1D*> gamma1_theta;
        vector<MnvH1D*> gamma1_ConvLength;
        TH1D* gamma1_true_E;
        TH1D* gamma1_reco_error_E;
        TH2D* gamma1_reco_E_true_E;
        TH2D* gamma1_true_E_reco_E_error;

        // Secondary Photon
        vector<MnvH1D*> gamma2_E;
        vector<MnvH1D*> gamma2_theta;
        vector<MnvH1D*> gamma2_ConvLength;
        TH1D* gamma2_true_E;
        TH1D* gamma2_reco_error_E;
        TH2D* gamma2_reco_E_true_E;
        TH2D* gamma2_true_E_reco_E_error;
      
        vector<MnvH1D*> invMass;
        vector<MnvH1D*> photonEnergy_Asymmetry;
        
        TH1D* signal_P;
        TH1D* signal_theta;
        TH2D* gamma1_E_gamma2_E;
        TH2D* gamma1_convLength_gamma2_convLength;
                
        CCProtonPi0_Pion(bool isModeReduce, bool isMC, std::string ana_folder);
        void initHistograms();
        void writeHistograms();
        
    private:
        CCProtonPi0_SingleBin bin_invMass;
        CCProtonPi0_SingleBin bin_photonEnergy_Asymmetry;
        CCProtonPi0_SingleBin bin_photonConvLength;
        CCProtonPi0_SingleBin bin_photonP;
};

#endif

