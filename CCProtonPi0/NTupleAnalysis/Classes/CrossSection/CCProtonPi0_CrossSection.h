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
        CCProtonPi0_CrossSection(bool isMC);
        void Calc_CrossSections();
    
    private:
        bool m_isMC;
        int iteration;
        double min_invMass;
        double max_invMass;
        double N_Background_Data;
        double Uncertainity_Bckg;
        double data_POT;
        double mc_POT;

        // Pi0 Invariant Mass
        TH1F* fit_result;
        MnvH1D* invMass_all;
        MnvH1D* invMass_mc_reco_signal;
        MnvH1D* invMass_mc_reco_bckg;

        // --------------------------------------------------------------------
        // Muon Momentum
        // --------------------------------------------------------------------
        // Data
        MnvH1D* muon_P_all;
        MnvH1D* muon_P_bckg_subtracted;
        MnvH1D* muon_P_bckg_estimated;
        MnvH1D* muon_P_unfolded;
        MnvH1D* muon_P_efficiency_corrected;
        MnvH1D* muon_P_integrated_flux;
        MnvH1D* muon_P_xsec;
        // MC Truth 
        MnvH1D* muon_P_mc_truth_all_signal;
        MnvH1D* muon_P_mc_truth_signal;
        MnvH1D* muon_P_mc_reco_signal;
        MnvH1D* muon_P_mc_reco_bckg;
        MnvH2D* muon_P_response;
        MnvH1D* muon_P_eff;

        // --------------------------------------------------------------------
        // Muon Theta 
        // --------------------------------------------------------------------
        // Data
        MnvH1D* muon_theta_all;
        MnvH1D* muon_theta_bckg_subtracted;
        MnvH1D* muon_theta_bckg_estimated;
        MnvH1D* muon_theta_unfolded;
        MnvH1D* muon_theta_efficiency_corrected;
        MnvH1D* muon_theta_integrated_flux;
        MnvH1D* muon_theta_xsec;
        // MC Truth 
        MnvH1D* muon_theta_mc_truth_all_signal;
        MnvH1D* muon_theta_mc_truth_signal;
        MnvH1D* muon_theta_mc_reco_signal;
        MnvH1D* muon_theta_mc_reco_bckg;
        MnvH2D* muon_theta_response;
        MnvH1D* muon_theta_eff;

        // --------------------------------------------------------------------
        // Pi0 Momentum 
        // --------------------------------------------------------------------
        // Data
        MnvH1D* pi0_P_all;
        MnvH1D* pi0_P_bckg_subtracted;
        MnvH1D* pi0_P_bckg_estimated;
        MnvH1D* pi0_P_unfolded;
        MnvH1D* pi0_P_efficiency_corrected;
        MnvH1D* pi0_P_integrated_flux;
        MnvH1D* pi0_P_xsec;
        // MC Truth
        MnvH1D* pi0_P_mc_truth_all_signal;
        MnvH1D* pi0_P_mc_truth_signal;
        MnvH1D* pi0_P_mc_reco_signal;
        MnvH1D* pi0_P_mc_reco_bckg;
        MnvH2D* pi0_P_response;
        MnvH1D* pi0_P_eff;
        
        // ROOT Files    
        TFile* f_out;
        TFile* f_truth;
        TFile* f_data_cutHists;
        TFile* f_data_muon;
        TFile* f_data_pi0;
        TFile* f_mc_cutHists;
        TFile* f_mc_muon;
        TFile* f_mc_pi0;
        
        std::string rootDir_out;
        std::string rootDir_flux;

        // Functions
        void Calc_CrossSection_muon_P();
        void Calc_CrossSection_muon_theta();
        void Calc_CrossSection_pi0_P();
        void Style_muon_P();
        void Style_muon_theta();
        void Style_pi0_P();
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
        MnvH1D* Calc_FinalCrossSection(MnvH1D* data_efficiency_corrected, MnvH1D* integrated_flux, std::string var_name);
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


