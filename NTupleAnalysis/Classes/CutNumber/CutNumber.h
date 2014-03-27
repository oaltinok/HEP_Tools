/*
================================================================================
Class: CutNumber
    Class for Cut Number
        Standard Class implementation for a member variable in CutNumberList
        Each Cut Number have private variables: value and label
        Value can be updated using increment() function
    
    Main Directory:
        Classes/CutNumber
    
    Usage:
        > #include "Classes/CutNumber/CutNumber.cpp"
        > CutNumber nAll;
        > value = nAll.getValue();
    
    Last Revision: 2014_02_14
================================================================================
*/

#ifndef CutNumber_h
#define CutNumber_h

#include<string>

const int initCount = 0.0;

class CutNumber
{
    public:
        CutNumber();
        CutNumber(string inputLabel);
        
        void increment();
        
        void setValue(double input);
        void setLabel(string input);
        
        double getValue();
        string getLabel();

    private:
        double value;
        string label;

};

#endif
