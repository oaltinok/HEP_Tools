/*
    See Pion.h header for Class Information
*/

#include "Pion.h"

void Pion::set_kineticEnergy(bool isMC)
{
    int type = getDataType(isMC);
    
    kineticEnergy[type] = p4[type].Energy() - restMass;
}

