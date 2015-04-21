/*
================================================================================
Class: Particle
    Particle Class defines a particle which will be used in the analysis
    Contains analysis specific information such as 4-Momentum and Angle wrt Beam
   
    Uses ROOT Specific classes
    
    Author:        Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision: 2015_04_21
================================================================================
*/

#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>
#include <string>
#include <cstdlib>

#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TFile.h>

// Libraries
#include "../../Libraries/Data_Functions.h"

//Classes
#include "../NTupleAnalysis/NTupleAnalysis.h"
#include "../SingleBin/SingleBin.h"

using namespace std;

class Particle : public NTupleAnalysis
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
        double energy[N_DATA_TYPE];         // Total Energy
        double kineticEnergy[N_DATA_TYPE];  // Kinetic energy of the Particle
        double angleBeam[N_DATA_TYPE];      // Angle wrt Beam in rads
        double angleMuon[N_DATA_TYPE];      // Angle wrt Muon in rads
        
        // MC Only
        int ind;                // indice for MC truth information (mc_FSPart)
        
        // Reco Only
        double particleScore;       // Particle Score from Reconstructed Values
        double trackLength;         // Track Length in [mm]
        
        TFile* f;
        
        // Histograms
        TH1D* partScore;
        
        TH1D* E_reco;
        TH1D* E_mc;
        TH1D* E_error;
        TH2D* E_reco_mc;
        
        TH1D* P_mc;
        TH1D* P_reco;
        TH1D* P_error;
        TH2D* P_reco_mc;
        
        TH1D* KE_mc;
        TH1D* KE_reco;
        TH1D* KE_error;
        TH2D* KE_reco_mc;
        
        TH1D* angleBeam_mc;
        TH1D* angleBeam_reco;
        TH1D* angleBeam_error;
        TH2D* angleBeam_reco_mc;
        
        TH1D* angleMuon_mc;
        TH1D* angleMuon_reco;
        TH1D* angleMuon_error;
        TH2D* angleMuon_reco_mc;
        
        // Bins for Histograms
        SingleBin bin_error;
        SingleBin bin_E;
        SingleBin bin_P;
        SingleBin bin_KE;
        SingleBin bin_angle;
        SingleBin bin_partScore;
          
        // File Locations
        std::string rootDir;
          
        // Functions
        Particle(int nMode);
        ~Particle();
        int getDataType(bool isMC); // Returns the indice of arrays reco = 0, mc = 1
        virtual void set_angleMuon(Particle &mu, bool isMC);
        virtual void set_kineticEnergy(bool isMC);
        void set_angleBeam(TVector3 beamp3, bool isMC);
        void set_p4(double px, double py, double pz, double E, bool isMC);
        void set_errors();
        void fill_Histograms();
        void write_RootFile();
        void plot_data();
        
    
    protected:
        double restMass;
        
        

};


#endif 
