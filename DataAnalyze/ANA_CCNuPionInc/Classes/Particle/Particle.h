/*
================================================================================
Class: Particle
    Particle Class defines a particle which will be used in the analysis
    Contains analysis specific information such as 4-Momentum and Angle wrt Beam
   
    Uses ROOT Specific classes
    
    Last Revision: 2014_01_28
================================================================================
*/

#ifndef PARTICLE_H
#define PARTICLE_H

#include<TMath.h>

// Two Channels in each variable
// Channel = 0 for Reco -- Default is Reconstructed Information
// Channel = 1 for True
const int N_CHANNELS = 2;

class Particle
{
    public:
        Particle();
        
        // Variables
        TLorentzVector* p4;     // 4-Momentum of the Particle (Px,Py,Pz,E)
        double* angleBeam;      // Angle wrt Beam in rads
        double* angleMuon;      // Angle wrt Muon in rads
        int ind;                // indice for MC truth information
    
    
    
    private:

};


#endif 
