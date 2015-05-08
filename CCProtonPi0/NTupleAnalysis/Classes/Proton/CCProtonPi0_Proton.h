/*
================================================================================
Class: CCProtonPi0_Proton -> Derived Class from Particle Base Clas
    CCProtonPi0_Proton Class  inherits Particle Behaviours and 
                extends base class with proton specific parameters

    Author:        Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision: 2015_05_07
================================================================================
*/
#ifndef CCProtonPi0_Proton_h
#define CCProtonPi0_Proton_h

#include "../Particle/CCProtonPi0_Particle.h"

class CCProtonPi0_Proton : public CCProtonPi0_Particle
{
    public:
        CCProtonPi0_Proton(int nMode);
        void set_kineticEnergy(bool isMC);
        
        TH1D* trackLength;
        TH1D* trackKinked;
                
    private:
        static const double restMass = 938.27;
        CCProtonPi0_SingleBin bin_trackLength;
        CCProtonPi0_SingleBin bin_trackKinked;
        
        
};


#endif
