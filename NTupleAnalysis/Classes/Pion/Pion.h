/*
================================================================================
Class: Pion -> Derived Class from Particle Base Clas
    Pion Class  inherits Particle Behaviours and 
                extends base class with pion specific parameters

    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision: 2014_05_10
================================================================================
*/

#ifndef Pion_h
#define Pion_h

#include "../Particle/Particle.h"

class Pion : public Particle
{
    public:
        Pion();
        void set_kineticEnergy(bool isMC);
        

    private:
        static const double restMass = 134.98;

};



#endif
