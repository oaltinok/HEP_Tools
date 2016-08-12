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
#include "../BinList/CCProtonPi0_BinList.h"
#include "TObjArray.h"
#include "TFractionFitter.h"
#include <MinervaUnfold/MnvUnfold.h>
#include <PlotUtils/TargetUtils.h>

using namespace PlotUtils;

struct XSec
{
    std::string name;
    bool isEv;
    
    // Style
    std::string plot_title;
    std::string plot_xlabel;
    std::string plot_ylabel;

    double smallest_bin_width;

    // ROOT Files
    TFile* f_data;
    TFile* f_mc;

    // Data Histograms
    MnvH1D* all;
    MnvH1D* bckg_subtracted;
    MnvH1D* bckg_estimated;
    MnvH1D* unfolded;
    MnvH1D* efficiency_corrected;
    MnvH1D* flux_integrated;
    MnvH1D* xsec;

    // MC Truth Histograms 
    MnvH1D* mc_truth_all_signal;
    MnvH1D* mc_truth_signal;
    MnvH1D* mc_reco_signal;
    MnvH1D* mc_reco_bckg;
    MnvH2D* response;
    MnvH1D* eff;
};

class CCProtonPi0_CrossSection : public CCProtonPi0_NTupleAnalysis
{
    public:
        CCProtonPi0_CrossSection(bool isMC);
        void Calc_CrossSections();
    
    private:
        CCProtonPi0_BinList binList;

        bool RemoveErrorBands;
        bool m_isMC;
        int iteration;
        double min_invMass;
        double max_invMass;
        std::vector<double> N_Background_Data;

        // Flux File
        MnvH1D* h_flux_minervaLE_FHC;
        MnvH1D* h_flux_rebinned;
        double cv_flux_integral;
        std::vector<double> unv_flux_integral;

        // Pi0 Invariant Mass
        MnvH1D* invMass_all;
        MnvH1D* invMass_mc_reco_signal;
        MnvH1D* invMass_mc_reco_bckg;

        XSec muon_P;
        XSec muon_theta;
        XSec pi0_P;
        XSec pi0_KE;
        XSec pi0_theta;
        XSec QSq;
        XSec W;
        XSec Enu;
        
        // Output Text File
        std::string text_out_name;
        ofstream text_out;

        // ROOT Files    
        TFile* f_out;
        TFile* f_truth;
        TFile* f_data_cutHists;
        TFile* f_mc_cutHists;
        
        std::string rootDir_out;
        
        // Functions
        void IntegrateAllFluxUniverses();
        void AddErrorBands_FluxHistogram();
        void RebinFluxHistogram();
        void RebinFluxHistogram(TH1* rebinned, TH1* reference);
        void Calc_CrossSection(XSec &var);
        void Style_XSec(XSec &var);
        void Calc_Normalized_NBackground();
        void NormalizeHistogram(TH1D* h);
        void NormalizeHistogram(MnvH1D* h);
        double Integrate_SignalRegion(TH1D* h);
        double GetFluxHistContent(TH1* hist, double low1, double low2);
        double GetSmallestBinWidth(MnvH1D* hist);
        void OpenRootFiles();
        void writeHistograms();
        void writeHistograms(XSec &var);
        void initHistograms();
        void initHistograms(XSec &var);
        void initFluxHistograms();
        void initXSecs();
        void init_muon_P();
        void init_muon_theta();
        void init_pi0_P();
        void init_pi0_KE();
        void init_pi0_theta();
        void init_W();
        void init_QSq();
        void init_Enu();
        MnvH1D* Subtract_Background(MnvH1D* data, MnvH1D* mc_bckg, MnvH1D* &bckg_estimated, std::string var_name);
        MnvH1D* Unfold_Data(MnvH1D* bckg_subtracted, MnvH2D* response, std::string var_name);
        MnvH1D* Efficiency_Divide(MnvH1D* unfolded, MnvH1D* eff, std::string var_name);
        MnvH1D* Integrate_Flux(MnvH1D* data_efficiency_corrected, std::string var_name);
        MnvH1D* Calc_FinalCrossSection(MnvH1D* flux_integrated, std::string var_name);
};


#endif


