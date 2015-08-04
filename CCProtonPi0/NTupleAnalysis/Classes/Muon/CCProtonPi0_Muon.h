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
        CCProtonPi0_Muon(bool isModeReduce, bool isMC, std::string ana_folder);
        void initHistograms();
        void writeHistograms();

    private:
        CCProtonPi0_SingleBin bin_muonTheta;
};

#endif
