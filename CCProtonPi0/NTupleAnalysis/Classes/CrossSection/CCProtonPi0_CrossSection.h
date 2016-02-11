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
#include <PlotUtils/TargetUtils.h>

using namespace PlotUtils;

class CCProtonPi0_CrossSection : public CCProtonPi0_NTupleAnalysis
{
    public:
        CCProtonPi0_CrossSection();
        void Calc_CrossSections();
    
    private:
        int iteration;
        double min_invMass;
        double max_invMass;
        double N_Background_Data;
        double Uncertainity_Bckg;
        double data_POT;
        double mc_POT;

        // Muon Momentum
        MnvH1D* data_xsec_muon_P;
        MnvH1D* data_all_muon_P;
        MnvH1D* data_bckg_subtracted_muon_P;
        MnvH1D* data_bckg_estimated_muon_P;
        MnvH1D* data_unfolded_muon_P;
        MnvH1D* data_efficiency_corrected_muon_P;
        MnvH1D* data_integrated_flux_muon_P;

        MnvH1D* mc_truth_xsec_muon_P;
        MnvH1D* mc_truth_all_signal_muon_P;
        MnvH1D* mc_truth_signal_muon_P;
        MnvH1D* mc_truth_integrated_flux_muon_P;
        MnvH1D* mc_reco_all_muon_P;
        MnvH1D* mc_reco_signal_muon_P;
        MnvH1D* mc_reco_bckg_muon_P;
        
        MnvH2D* response_muon_P;
        MnvH1D* eff_muon_P;

        // Data Histograms
        MnvH1D* data_xsec_pi0_P;
        MnvH1D* data_all_pi0_P;
        MnvH1D* data_bckg_subtracted_pi0_P;
        MnvH1D* data_bckg_estimated_pi0_P;
        MnvH1D* data_unfolded_pi0_P;
        MnvH1D* data_efficiency_corrected_pi0_P;
        MnvH1D* data_integrated_flux_pi0_P;

        MnvH1D* mc_truth_xsec_pi0_P;
        MnvH1D* mc_truth_all_signal_pi0_P;
        MnvH1D* mc_truth_signal_pi0_P;
        MnvH1D* mc_truth_integrated_flux_pi0_P;
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
        std::string rootDir_flux;

        // Functions
        void Label_CrossSection_Hists();
        void Calc_CrossSection_Muon_P();
        void Calc_CrossSection_Pi0_P();
        void Calc_Normalized_NBackground();
        void NormalizeHistogram(MnvH1D* h);
        double Integrate_SignalRegion(TH1D* h);
        void writeHistograms();
        void OpenRootFiles();
        void initHistograms();
        MnvH1D* Subtract_Background(MnvH1D* data, MnvH1D* mc_bckg, MnvH1D* &bckg_estimated, std::string var_name);
        MnvH1D* Unfold_Data(MnvH1D* bckg_subtracted, MnvH2D* response, std::string var_name);
        MnvH1D* Efficiency_Divide(MnvH1D* unfolded, MnvH1D* eff, std::string var_name);
        MnvH1D* Integrate_Flux(MnvH1D* data_efficiency_corrected, std::string var_name, bool isEv);
        MnvH1D* Calc_CrossSection(MnvH1D* data_efficiency_corrected, MnvH1D* integrated_flux, double pot, bool isMC);
        MnvH1D* calc_flux( MnvH1D* mnvh1d_template,          // Template histogram to copy the binning for the flux histogram 
                const std::string& flux_filename, // The flux file name
                bool enu_cut,                     // whether to apply the enu cut
                bool isEv,                        // if this variable is Enu since it is treated differently
                bool __reweight_flux = false,     // whether to do the flux reweight study
                double __reweight_emin = 0.0,     // lower bound of the reweighted region 
                double __reweight_emax = 0.0,     // upper bound of the reweighted region
                double __reweight_amount = 0.0,   // amount to reweight
                bool flux_smooth_curve = false,   // reweight the flux over full neutrino energy range
                bool flux_piecewise_curve = false);


};


#endif


