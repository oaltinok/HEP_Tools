/*
================================================================================
Class: Proton -> Derived Class from Particle Base Clas
    Proton Class  inherits Particle Behaviours and 
                extends base class with proton specific parameters

    Author:        Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision: 2014_05_08
================================================================================
*/
#ifndef Proton_h
#define Proton_h

#include "Classes/Particle/Particle.cpp"

class Proton : public Particle
{
    public:
        Proton();
        
    private:
        static const double restMass = 938.27;
        

};


#endif
