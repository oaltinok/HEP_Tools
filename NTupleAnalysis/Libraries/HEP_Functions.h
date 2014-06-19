/*
================================================================================
Namespace: HEP_Functions
    Common High Energy Functions 
    Miscellaneous functions which are not specific to any analysis package 
    
    Main Directory:
        Libraries/
    
    Usage:
        > #include "Libraries/HEP_Functions.h" 
        > HEP_Functions::calcEnergy(momentum,restMass);
    
    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision:  2014_06_18
================================================================================
*/

#ifndef HEP_Functions_h
#define HEP_Functions_h

namespace HEP_Functions
{

inline double calcEnergy(double P, double restMass)
{
    double energy;
    energy = sqrt(P*P + restMass*restMass);
    return energy;
}

inline double calcMomentum(double px, double py, double pz)
{
    double P;
    P = sqrt(px*px + py*py + pz*pz);
    return P;
}

    
}

#endif

