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
        std::vector<MnvH1D*> E;
        std::vector<MnvH1D*> P;
        std::vector<MnvH1D*> P_1Track;
        std::vector<MnvH1D*> P_2Track;
        std::vector<MnvH1D*> KE;
        std::vector<MnvH1D*> theta;
        std::vector<MnvH1D*> theta_1Track;
        std::vector<MnvH1D*> theta_2Track;
        std::vector<MnvH1D*> cos_theta;
        std::vector<MnvH1D*> phi;

        TH1D* theta_error;
        TH1D* cos_theta_error;
        TH1D* KE_error;
        TH1D* P_error;
        TH1D* E_error;

        TH1D* theta_Diff;
        TH1D* P_Diff;
        TH1D* E_Diff;
    
        TH2D* reco_P_true_P;
        TH2D* reco_E_true_E;
        TH1D* E_true;
        TH1D* E_reco;

        // Bins for Histograms
        CCProtonPi0_BinList binList;
        CCProtonPi0_SingleBin bin_E;
        CCProtonPi0_SingleBin bin_P;
        CCProtonPi0_SingleBin bin_KE;
        CCProtonPi0_SingleBin bin_E_Diff;
        CCProtonPi0_SingleBin bin_P_Diff;
        CCProtonPi0_SingleBin bin_theta_Diff;
          
        // File Locations
        std::string rootDir;
          
        // Functions
        CCProtonPi0_Particle();
        ~CCProtonPi0_Particle();
        virtual void writeHistograms() = 0;
        virtual void initHistograms() = 0; 
};


#endif 
