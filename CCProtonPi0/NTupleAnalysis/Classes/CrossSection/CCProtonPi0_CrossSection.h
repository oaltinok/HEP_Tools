/*
 * =============================================================================
 * Class: CCProtonPi0_CrossSection
 *   CCProtonPi0_CrossSection Class contains ALL required methods to calculate
 *       Cross Section for CCProtonPi0 Analysis
 *  
 *   Uses ROOT and MINERvA Specific classes
 *   
 *   Author: Ozgur Altinok  - ozgur.altinok@tufts.edu
 * =============================================================================
 */
#ifndef CCProtonPi0_CrossSection_h
#define CCProtonPi0_CrossSection_h

// Classes
#include "../NTupleAnalysis/CCProtonPi0_NTupleAnalysis.h"
#include "TObjArray.h"
#include "TFractionFitter.h"
#include <MinervaUnfold/MnvUnfold.h>

using namespace PlotUtils;

class CCProtonPi0_CrossSection : public CCProtonPi0_NTupleAnalysis
{
    public:
        CCProtonPi0_CrossSection();
        void Calc_Crossections();
    
    private:
        int iteration;
        double min_invMass;
        double max_invMass;
        double N_Background_Data;
        double Uncertainity_Bckg;

        // Muon Momentum
        MnvH1D* data_all_muon_P;
        MnvH1D* data_bckg_subtracted_muon_P;
        MnvH1D* data_unfolded_muon_P;
        MnvH1D* data_efficiency_corrected_muon_P;

        MnvH1D* mc_truth_all_signal_muon_P;
        MnvH1D* mc_truth_signal_muon_P;
        MnvH1D* mc_reco_all_muon_P;
        MnvH1D* mc_reco_signal_muon_P;
        MnvH1D* mc_reco_bckg_muon_P;
        
        MnvH2D* response_muon_P;
        MnvH1D* eff_muon_P;
        
        // Data Histograms
        MnvH1D* data_all_pi0_P;
        MnvH1D* data_bckg_subtracted_pi0_P;
        MnvH1D* data_unfolded_pi0_P;
        MnvH1D* data_efficiency_corrected_pi0_P;

        MnvH1D* mc_truth_all_signal_pi0_P;
        MnvH1D* mc_truth_signal_pi0_P;
        MnvH1D* mc_reco_all_pi0_P;
        MnvH1D* mc_reco_signal_pi0_P;
        MnvH1D* mc_reco_bckg_pi0_P;

        MnvH2D* response_pi0_P;
        MnvH1D* eff_pi0_P;

        // ROOT Files    
        TFile* f_out;
        TFile* f_truth;
        TFile* f_data_muon;
        TFile* f_data_pi0;
        TFile* f_mc_muon;
        TFile* f_mc_pi0;
        
        std::string rootDir_out;

        // Functions
        void Calc_CrossSection_Muon_P();
        void Calc_CrossSection_Pi0_P();
        void Subtract_Background(MnvH1D* &bckg_subtracted, MnvH1D* data, MnvH1D* mc_bckg, std::string var_name);
        void Unfold_Data(MnvH1D* &unfolded, MnvH1D* bckg_subtracted, MnvH2D* response, std::string var_name);
        void Efficiency_Divide(MnvH1D* &efficiency_corrected, MnvH1D* unfolded, MnvH1D* eff, std::string var_name);
        void Calc_Normalized_NBackground();
        void NormalizeHistogram(MnvH1D* h);
        double Integrate_SignalRegion(TH1D* h);
        void writeHistograms();
        void OpenRootFiles();
        void initHistograms();
};


#endif


