/*
================================================================================
Class: CCProtonPi0_Pi0Blob
    CCProtonPi0_Pi0Blob Class Contains Blob Histograms
    All Histograms declared public and can be accessed by Analyzer
                
    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
================================================================================
*/
#ifndef CCProtonPi0_Pi0Blob_h
#define CCProtonPi0_Pi0Blob_h

#include "../NTupleAnalysis/CCProtonPi0_NTupleAnalysis.h"
#include "../BinList/CCProtonPi0_BinList.h"

using namespace PlotUtils;

class CCProtonPi0_Pi0Blob : public CCProtonPi0_NTupleAnalysis
{
    public:
        CCProtonPi0_Pi0Blob(bool isModeReduce, bool isMC, std::string ana_folder);

        void writeHistograms();

        // Histograms
       
        // Truth Match
        std::vector<MnvH1D*> g1_evis_most_pdg;
        std::vector<MnvH1D*> g1_evis_total_truth;
        std::vector<MnvH1D*> g1_evis_frac_pizero;
        std::vector<MnvH1D*> g1_evis_frac_piplus;
        std::vector<MnvH1D*> g1_evis_frac_piminus;
        std::vector<MnvH1D*> g1_evis_frac_proton;
        std::vector<MnvH1D*> g1_evis_frac_neutron;
        std::vector<MnvH1D*> g1_evis_frac_muon;
       
        std::vector<MnvH1D*> g2_evis_most_pdg;
        std::vector<MnvH1D*> g2_evis_total_truth;
        std::vector<MnvH1D*> g2_evis_frac_pizero;
        std::vector<MnvH1D*> g2_evis_frac_piplus;
        std::vector<MnvH1D*> g2_evis_frac_piminus;
        std::vector<MnvH1D*> g2_evis_frac_proton;
        std::vector<MnvH1D*> g2_evis_frac_neutron;
        std::vector<MnvH1D*> g2_evis_frac_muon;
        
        std::vector<MnvH1D*> captured_evis_frac_all; 
        std::vector<MnvH1D*> captured_evis_frac_signal; 

        // Reco Values
        std::vector<MnvH1D*> g1_nPlanes; 
        std::vector<MnvH1D*> g1_cluster_Z; 
        std::vector<MnvH1D*> g1_max_cluster_Z; 
        std::vector<MnvH1D*> g1_min_cluster_Z; 
        std::vector<MnvH1D*> g1_strips; 
        std::vector<MnvH1D*> g1_max_strip; 
        std::vector<MnvH1D*> g1_min_strip; 
        
        std::vector<MnvH1D*> g2_nPlanes; 
        std::vector<MnvH1D*> g2_cluster_Z; 
        std::vector<MnvH1D*> g2_max_cluster_Z; 
        std::vector<MnvH1D*> g2_min_cluster_Z; 
        std::vector<MnvH1D*> g2_strips; 
        std::vector<MnvH1D*> g2_max_strip; 
        std::vector<MnvH1D*> g2_min_strip; 
    
    private:
        void initBins();
        void initHistograms();

        TFile* f;
        string rootDir;
        
        CCProtonPi0_BinList binList;
        
        // Histogram Binning
        CCProtonPi0_SingleBin bin_blob_energy;

};

#endif

