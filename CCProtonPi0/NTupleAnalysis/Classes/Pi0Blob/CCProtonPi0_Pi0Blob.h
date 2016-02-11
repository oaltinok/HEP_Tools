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
        CCProtonPi0_Pi0Blob(bool isModeReduce, bool isMC);

        void writeHistograms();

        // Histograms

        // Truth Match
        // Pi0 Capture Performance
        std::vector<MnvH1D*> evis_frac_true_pi0_reco_all;
        std::vector<MnvH1D*> evis_frac_reco_pi0_true_pi0;
        std::vector<MnvH1D*> evis_frac_reco_pi0_reco_all;
        std::vector<MnvH1D*> evis_frac_reco_nonpi0_reco_all;
         
        // Evis Fractions and PDG of Particle with most Evis
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

        // Evis from each Particle for Stacked Plot
        TH1D* g1_evis_proton;
        TH1D* g1_evis_neutron;
        TH1D* g1_evis_pi;
        TH1D* g1_evis_pi0;
        TH1D* g1_evis_muon;

        TH1D* g2_evis_proton;
        TH1D* g2_evis_neutron;
        TH1D* g2_evis_pi;
        TH1D* g2_evis_pi0;
        TH1D* g2_evis_muon;

        TH1D* g3_evis_proton;
        TH1D* g3_evis_neutron;
        TH1D* g3_evis_pi;
        TH1D* g3_evis_pi0;
        TH1D* g3_evis_muon;
       
        // Reco Values
        std::vector<MnvH1D*> g1_nPlanes; 
        std::vector<MnvH1D*> g2_nPlanes; 
            
    private:
        void initBins();
        void initHistograms();

        TFile* f;
        std::string rootDir;
        
        CCProtonPi0_BinList binList;
        
        // Histogram Binning
        CCProtonPi0_SingleBin bin_blob_energy;

};

#endif

