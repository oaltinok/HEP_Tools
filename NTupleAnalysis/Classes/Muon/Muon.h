/*
================================================================================
Class: Muon -> Derived Class from Particle Base Clas
    Muon Class  inherits Particle Behaviours and 
                extends base class with muon specific parameters
    
    Last Revision: 2014_05_08
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
        static const double restMass = 105.66;
        bool isMinosMatched;
        SingleBin bin_AngleBeam;
        

         
        
        
};

#endif
