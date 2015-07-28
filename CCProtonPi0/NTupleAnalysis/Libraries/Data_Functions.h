/*
================================================================================
Namespace: Data_Functions
    Common Data Analysis Functions
    Miscellaneous functions which are not specific to any analysis package 
    
    Main Directory:
        Libraries/
    
    Usage:
        > #include "Libraries/Data_Functions.h" 
        > Data_Functions::getError(trueValue,recoValue)
    
    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision:  2014_05_10
================================================================================
*/

#ifndef Data_Functions_h
#define Data_Functions_h

namespace Data_Functions
{

inline double getError(double trueValue, double recoValue)
{
    double error;
    error = (recoValue - trueValue) / trueValue;
    return error;
}

    
}

#endif

