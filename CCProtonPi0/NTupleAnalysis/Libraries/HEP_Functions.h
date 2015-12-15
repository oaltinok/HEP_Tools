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
    Last Revision:  2014_06_19
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

inline double calcDistance( double x1, double y1, double z1,
                            double x2, double y2, double z2)
{
    double d;
    d = sqrt(   (x1-x2)*(x1-x2) + 
                (y1-y2)*(y1-y2) + 
                (z1-z2)*(z1-z2));
    return d;
}
    
}

#endif

