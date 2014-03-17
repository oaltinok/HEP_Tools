/*
================================================================================
Class: Particle
    Particle Class defines a particle which will be used in the analysis
    Contains analysis specific information such as 4-Momentum and Angle wrt Beam
   
    Uses ROOT Specific classes
    
    Last Revision: 2014_03_06
================================================================================
*/

#ifndef PARTICLE_H
#define PARTICLE_H

#include<TMath.h>


class Particle
{
    public:
        // Two Data Data Types
        // type = 0 -- Reco -- Reconstructed 
        // type = 1 -- True -- Monte Carlo
        // type = 2 -- Error -- (True - Reco) / True
        static const int N_DATA_TYPE = 3;
        
        // Monte Carlo(MC) and Reconstructed(Reco)
        TLorentzVector p4[N_DATA_TYPE];     // 4-Momentum of the Particle (Px,Py,Pz,E)
        double momentum[N_DATA_TYPE];       // 3 Momentum of the Particle (P)
        double kineticEnergy[N_DATA_TYPE];  // Kinetic energy of the Particle
        double angleBeam[N_DATA_TYPE];      // Angle wrt Beam in rads
        double angleMuon[N_DATA_TYPE];      // Angle wrt Muon in rads
        
        // MC Only
        int ind;                // indice for MC truth information (mc_FSPart)
        
        // Reco Only
        double pID;             // Particle Score from Reconstructed Values
        double trackLength;     // Track Length in [mm]
        
        
        // Histograms
        
        // Functions
        Particle();
        ~Particle();
        int getDataType(bool isMC); // Returns the indice of arrays reco = 0, mc = 1
        virtual void set_angleMuon(Particle &mu, bool isMC);
        virtual void set_kineticEnergy(bool isMC) = 0;
        void set_angleBeam(TVector3 beamp3, bool isMC);
        void set_p4(double px, double py, double pz, double E, bool isMC);
        void set_errors();
        
        
        // Other
        double calc_error(double trueValue, double recoValue);
        
       
        
    
    protected:
        double restMass;
        
        

};


#endif 
