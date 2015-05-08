/*
================================================================================
Class: CCProtonPi0_Muon -> Derived Class from Particle Base Clas
    CCProtonPi0_Muon Class  inherits Particle Behaviours and 
                extends base class with Muon specific parameters
                
    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision: 2015_05_07
================================================================================
*/

#ifndef CCProtonPi0_Muon_h
#define CCProtonPi0_Muon_h

#include "../Particle/CCProtonPi0_Particle.h"

class CCProtonPi0_Muon : public CCProtonPi0_Particle
{
    public:
        CCProtonPi0_Muon(int nMode);
        bool get_isMinosMatched();
        void set_isMinosMatched(bool input);
        void set_angleMuon(CCProtonPi0_Particle &mu, bool isMC);
        void set_kineticEnergy(bool isMC);
        
    private:
        static const double restMass = 105.66;
        bool isMinosMatched;
        CCProtonPi0_SingleBin bin_AngleBeam;
        
};

#endif
