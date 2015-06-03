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

class CCProtonPi0_Pion : public CCProtonPi0_Particle
{
    public:
        CCProtonPi0_Pion(int nMode);
        ~CCProtonPi0_Pion();
        void set_kineticEnergy(bool isMC);

        TH1D* gamma1_ConvLength;
        TH1D* gamma2_ConvLength;
        TH2D* ConvLength_gamma2_gamma1;
           
        TH1D* photonEnergy_Asymmetry;
        TH1D* photonEnergy_Asymmetry_true;

        TH1D* invMass;
 
    private:
        void initHistograms();

        static const double restMass = 134.98;
        
        CCProtonPi0_SingleBin bin_invMass;
        CCProtonPi0_SingleBin bin_blob_energy;
        CCProtonPi0_SingleBin bin_photonConvLength;

};

#endif

