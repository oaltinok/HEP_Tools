#ifndef CCProtonPi0_Pi0BlobTool_h
#define CCProtonPi0_Pi0BlobTool_h

#include <iostream>
#include <vector>

//-- Forward Declarations
#include "Event/MinervaEventFwd.h"

class Pi0BlobTool
{
    public:
        Pi0BlobTool();
        bool isBlob1_Good(Minerva::IDBlob* blob);
        bool isBlob2_Good(Minerva::IDBlob* blob);

    private:
        static const double minNDigits_blob1 = 8;
        static const double minNDigits_blob2 = 6;
        static const double minEnergy_blob1 = 60;     // MeV
        
        bool isNDigitsHigh(Minerva::IDBlob* blob, const double minNDigits);
        bool isEnergyHigh(Minerva::IDBlob* blob, const double minEnergy);
};

#endif

