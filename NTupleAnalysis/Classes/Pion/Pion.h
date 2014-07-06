/*
================================================================================
Class: Pion -> Derived Class from Particle Base Clas
    Pion Class  inherits Particle Behaviours and 
                extends base class with pion specific parameters

    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision: 2014_07_06
================================================================================
*/

#ifndef Pion_h
#define Pion_h

#include "../Particle/Particle.h"

class Pion : public Particle
{
    public:
        Pion();
        void set_kineticEnergy(bool isMC);
        
        TH1F* invMass;
        
        TH1F* P_reco_0Pi0;
        TH1F* P_reco_1Pi0;
        TH1F* P_reco_MultPi0;
        
        TH2F* P_reco_mc_1Pi0;
        TH2F* P_reco_mc_MultPi0;
        
        TH1F* P_error_1Pi0;
        TH1F* P_error_MultPi0;
        

    private:
        static const double restMass = 134.98;
        
        SingleBin bin_invMass;

};



#endif
