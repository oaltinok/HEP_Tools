/*
================================================================================
Class: CCProtonPi0_Particle
    CCProtonPi0_Particle Class defines a particle which will be used in the analysis
    Contains analysis specific information such as 4-Momentum and Angle wrt Beam
   
    Uses ROOT Specific classes
    
    Author:        Ozgur Altinok  - ozgur.altinok@tufts.edu
================================================================================
*/

#ifndef CCProtonPi0_Particle_h
#define CCProtonPi0_Particle_h

#include <iostream>
#include <string>
#include <cstdlib>

#include <TFile.h>
#include <PlotUtils/MnvH1D.h>
#include <PlotUtils/MnvH2D.h>

//Classes
#include "../NTupleAnalysis/CCProtonPi0_NTupleAnalysis.h"
#include "../BinList/CCProtonPi0_BinList.h"

using namespace PlotUtils;

class CCProtonPi0_Particle : public CCProtonPi0_NTupleAnalysis
{
    public:
        TFile* f;
        
        // Standard Histograms
        MnvH1D* E;
        MnvH1D* P;
        MnvH1D* KE;
        MnvH1D* theta;
        MnvH1D* phi;
        
        // Bins for Histograms
        CCProtonPi0_BinList binList;
        CCProtonPi0_SingleBin bin_E;
        CCProtonPi0_SingleBin bin_P;
        CCProtonPi0_SingleBin bin_KE;
          
        // File Locations
        std::string rootDir;
          
        // Functions
        CCProtonPi0_Particle(int nMode);
        ~CCProtonPi0_Particle();
        virtual void writeHistograms() = 0;
        virtual void initHistograms() = 0; 
};


#endif 
