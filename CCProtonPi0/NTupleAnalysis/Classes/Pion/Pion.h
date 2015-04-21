/*
================================================================================
Class: Pion -> Derived Class from Particle Base Clas
    Pion Class  inherits Particle Behaviours and 
                extends base class with pion specific parameters

    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision: 2015_04_21
================================================================================
*/

#ifndef Pion_h
#define Pion_h

#include "../Particle/Particle.h"

class Pion : public Particle
{
    public:
        Pion(int nMode);
        ~Pion();
        void set_kineticEnergy(bool isMC);
        
        SingleBin bin_photonConvLength;
        SingleBin bin_gammaClusters;
        SingleBin bin_gammaEnergy;
        
        TH1D* gamma1_nClusters_All;
        TH1D* gamma2_nClusters_All;
        TH2D* nClusters_All_gamma2_gamma1;
        
        TH1D* gamma1_ConvLength;
        TH1D* gamma2_ConvLength;
        TH2D* ConvLength_gamma2_gamma1;

        TH1D* gamma1_Energy;
        TH1D* gamma2_Energy;
        TH2D* Energy_gamma2_gamma1;
           
        TH1D* photonEnergy_Asymmetry;

        TH1D* invMass;
    private:
        static const double restMass = 134.98;
        
        SingleBin bin_invMass;

};



#endif
