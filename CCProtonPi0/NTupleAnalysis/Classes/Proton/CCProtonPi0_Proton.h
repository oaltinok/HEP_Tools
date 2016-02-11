/*
================================================================================
Class: CCProtonPi0_Proton -> Derived Class from Particle Base Clas
    CCProtonPi0_Proton Class  inherits Particle Behaviours and 
                extends base class with proton specific parameters

    Author:        Ozgur Altinok  - ozgur.altinok@tufts.edu
================================================================================
*/
#ifndef CCProtonPi0_Proton_h
#define CCProtonPi0_Proton_h

#include "../Particle/CCProtonPi0_Particle.h"

using namespace PlotUtils;

class CCProtonPi0_Proton : public CCProtonPi0_Particle
{
    public:
        std::vector<MnvH1D*> trackLength;
        std::vector<MnvH1D*> trackKinked;
        std::vector<MnvH1D*> partScore;

        CCProtonPi0_Proton(bool isModeReduce, bool isMC);
        void initHistograms();
        void writeHistograms();
    
    private:
        CCProtonPi0_SingleBin bin_trackLength;
        CCProtonPi0_SingleBin bin_trackKinked;
};


#endif
