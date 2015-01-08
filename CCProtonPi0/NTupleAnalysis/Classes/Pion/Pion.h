/*
================================================================================
Class: Pion -> Derived Class from Particle Base Clas
    Pion Class  inherits Particle Behaviours and 
                extends base class with pion specific parameters

    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision: 2014_12_07
================================================================================
*/

#ifndef Pion_h
#define Pion_h

#include "../Particle/Particle.h"

class Pion : public Particle
{
    public:
        Pion();
        ~Pion();
        void set_kineticEnergy(bool isMC);
        void initialize(int nMode);
        
        SingleBin bin_photonConvLength;
        SingleBin bin_gammaClusters;
        SingleBin bin_gammaEnergy;
        
        TH1D* gamma1_nClusters_All;
        TH1D* gamma2_nClusters_All;
        TH2D* nClusters_All_gamma2_gamma1;
        
        TH1D* gamma1_nClusters_X;
        TH1D* gamma2_nClusters_X;
        TH2D* nClusters_X_gamma2_gamma1;
        
        TH1D* gamma1_ConvLength;
        TH1D* gamma2_ConvLength;
        TH2D* ConvLength_gamma2_gamma1;
        
        TH1D* gamma1_Energy;
        TH1D* gamma2_Energy;
        TH2D* Energy_gamma2_gamma1;
           
        TH1D* photonEnergy_Asymmetry;

        TH1D* invMass;
        TH1D* invMass_0Pi0;
        TH1D* invMass_1Pi0;
        TH1D* invMass_MultPi0;
        
        TH1D* P_reco_0Pi0;
        TH1D* P_reco_1Pi0;
        TH1D* P_reco_MultPi0;
        
        TH2D* P_reco_mc_1Pi0;
        TH2D* P_reco_mc_MultPi0;
        
        TH1D* P_error_1Pi0;
        TH1D* P_error_MultPi0;
        

    private:
        static const double restMass = 134.98;
        
        SingleBin bin_invMass;

};



#endif
