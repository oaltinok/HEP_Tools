/*
================================================================================
Class: CCProtonPi0_Particle
    CCProtonPi0_Particle Class defines a particle which will be used in the analysis
   
    Uses ROOT and MINERvA Specific classes
    
    Author:        Ozgur Altinok  - ozgur.altinok@tufts.edu
================================================================================
*/

#ifndef CCProtonPi0_Particle_h
#define CCProtonPi0_Particle_h

//Classes
#include "../NTupleAnalysis/CCProtonPi0_NTupleAnalysis.h"
#include "../BinList/CCProtonPi0_BinList.h"

using namespace PlotUtils;

class CCProtonPi0_Particle : public CCProtonPi0_NTupleAnalysis
{
    public:
        TFile* f;
        
        // Standard Histograms
        vector<MnvH1D*> E;
        vector<MnvH1D*> P;
        vector<MnvH1D*> KE;
        vector<MnvH1D*> theta;
        vector<MnvH1D*> phi;

        TH2D* reco_P_true_P;
        TH1D* P_error;

        // Bins for Histograms
        CCProtonPi0_BinList binList;
        CCProtonPi0_SingleBin bin_E;
        CCProtonPi0_SingleBin bin_P;
        CCProtonPi0_SingleBin bin_KE;
          
        // File Locations
        std::string rootDir;
          
        // Functions
        CCProtonPi0_Particle();
        ~CCProtonPi0_Particle();
        virtual void writeHistograms() = 0;
        virtual void initHistograms() = 0; 
};


#endif 
