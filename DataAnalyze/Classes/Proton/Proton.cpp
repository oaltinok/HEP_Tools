/*
    See Proton.h header for Class Information
*/
#include "Proton.h"

void Proton::set_kineticEnergy(bool isMC)
{
    int type = getDataType(isMC);
    
    kineticEnergy[type] = p4[type].Energy() - restMass;
}


