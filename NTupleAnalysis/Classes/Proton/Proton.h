/*
================================================================================
Class: Proton -> Derived Class from Particle Base Clas
    Proton Class  inherits Particle Behaviours and 
                extends base class with proton specific parameters

    Author:        Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision: 2014_07_16
================================================================================
*/
#ifndef Proton_h
#define Proton_h

#include "../Particle/Particle.h"

class Proton : public Particle
{
    public:
        Proton();
        void set_kineticEnergy(bool isMC);
        
    private:
        static const double restMass = 938.27;
        

};


#endif
