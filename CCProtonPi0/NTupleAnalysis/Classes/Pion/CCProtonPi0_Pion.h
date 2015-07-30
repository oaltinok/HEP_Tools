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
        vector<MnvH1D*> invMass;
        vector<MnvH1D*> photonEnergy_Asymmetry;

        // Leading Photon - Energetic Photon
        vector<MnvH1D*> gamma1_P;
        vector<MnvH1D*> gamma1_theta;
        vector<MnvH1D*> gamma1_ConvLength;
        // Secondary Photon
        vector<MnvH1D*> gamma2_P;
        vector<MnvH1D*> gamma2_theta;
        vector<MnvH1D*> gamma2_ConvLength;

        // Photon Comparison
        TH2D* gamma1_P_gamma2_P;
        TH2D* gamma1_convLength_gamma2_convLength;
        TH2D* gamma1_reco_P_true_P;
        TH2D* gamma2_reco_P_true_P;
        TH1D* gamma2_P_error;
        TH1D* gamma1_P_error;

        CCProtonPi0_Pion(bool isModeReduce, bool isMC);
        void initHistograms();
        void writeHistograms();

    private:
        CCProtonPi0_SingleBin bin_invMass;
        CCProtonPi0_SingleBin bin_photonEnergy_Asymmetry;
        CCProtonPi0_SingleBin bin_photonConvLength;
        CCProtonPi0_SingleBin bin_photonP;

};

#endif

