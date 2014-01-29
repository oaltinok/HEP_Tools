/*
================================================================================
Class: HEP_Misc
    High Energy Physics Miscellaneous Functions
    Miscellaneous functions which are not specific to any analysis package 
    
    Main Directory:
        Classes/HEP_Misc/
    
    Usage:
        > #include "Classes/HEP_Misc/HEP_Misc.cpp" 
        > HEP_Misc misc;
        > misc.getError(trueValue,recoValue)
    
    Last Revision: 2014_01_29
================================================================================
*/

#ifndef HEP_Misc_h
#define HEP_Misc_h

class HEP_Misc
{
    public:
        HEP_Misc();
        
        double getError(double trueValue, double recoValue);

    private:


};

#endif
