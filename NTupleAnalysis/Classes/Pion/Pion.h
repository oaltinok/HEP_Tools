/*
================================================================================
Class: Pion -> Derived Class from Particle Base Clas
    Pion Class  inherits Particle Behaviours and 
                extends base class with pion specific parameters

    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision: 2014_12_01
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
        
        TH1D* photonConvLength;
        
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
