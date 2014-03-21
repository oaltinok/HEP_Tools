/*
================================================================================
Class: Pion -> Derived Class from Particle Base Clas
    Pion Class  inherits Particle Behaviours and 
                extends base class with pion specific parameters
    
    Last Revision: 2014_03_20
================================================================================
*/

#ifndef Pion_h
#define Pion_h

#include "Classes/Particle/Particle.cpp"

class Pion : public Particle
{
    public:
        Pion() : Particle() {};


    private:
        static const double restMass = 134.98;


};



#endif
