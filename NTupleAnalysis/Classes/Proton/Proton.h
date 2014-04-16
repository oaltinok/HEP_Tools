/*
================================================================================
Class: Proton -> Derived Class from Particle Base Clas
    Proton Class  inherits Particle Behaviours and 
                extends base class with proton specific parameters
    
    Last Revision: 2014_04_16
================================================================================
*/
#ifndef Proton_h
#define Proton_h

#include "Classes/Particle/Particle.cpp"

class Proton : public Particle
{
    public:
        Proton() : Particle() {};
        
    private:
        static const double restMass = 938.27;
        

};


#endif
