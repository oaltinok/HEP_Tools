/*
================================================================================
Class: Muon -> Derived Class from Particle Base Clas
    Muon Class  inherits Particle Behaviours and 
                extends base class with muon specific parameters
    
    Last Revision: 2014_02_15
================================================================================
*/

#ifndef Muon_h
#define Muon_h

#include "Classes/Particle/Particle.cpp"

class Muon : public Particle
{
    public:
        Muon();
        bool get_isMinosMatched();
        void set_isMinosMatched(bool input);
        
    private:
        bool isMinosMatched;
        
};

#endif
