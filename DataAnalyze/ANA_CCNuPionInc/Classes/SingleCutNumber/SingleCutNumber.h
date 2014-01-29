/*
================================================================================
Class: SingleCutNumber
    Class for Cut Number
        Standard Class implementation for a member variable in CutNumberList
        Each Cut Number have private variables: value and label
        Value can be updated using increment() function
    
    Main Directory:
        Classes/SingleCutNumber
    
    Usage:
        > #include "Classes/SingleCutNumber/SingleCutNumber.cpp"
        > SingleCutNumber nAll;
        > value = nAll.getValue();
    
    Last Revision: 2014_01_29
================================================================================
*/

#ifndef SingleCutNumber_h
#define SingleCutNumber_h

#include<string>

class SingleCutNumber
{
    public:
        SingleCutNumber();
        
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
