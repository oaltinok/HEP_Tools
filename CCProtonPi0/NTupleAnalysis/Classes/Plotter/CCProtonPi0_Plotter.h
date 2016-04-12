/*
   ================================================================================
Class: CCProtonPi0_Plotter
CCProtonPi0_Plotter class includes specific functions for plotting 1D or 2D histograms

Main Directory:
Classes/Plotter/

Usage:
> #include "Classes/Plotter/CCProtonPi0_Plotter.cpp" 
> CCProtonPi0_Plotter p;
> p.plotHistograms(std::string mcFile, std::string plotDir)


Author:        Ozgur Altinok  - ozgur.altinok@tufts.edu
================================================================================
*/

#ifndef CCProtonPi0_Plotter_H
#define CCProtonPi0_Plotter_h

#include <TVectorD.h>
#include <TStyle.h>
#include <TPaveStats.h>
#include <TPad.h>
#include <TLegend.h>
#include <THStack.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <PlotUtils/MnvPlotter.h>
#include <PlotUtils/MnvFluxConstraint.h>
#include <PlotUtils/POTCounter.h>
#include <MinervaUnfold/MnvUnfold.h>
#include "Cintex/Cintex.h"

#include "../NTupleAnalysis/CCProtonPi0_NTupleAnalysis.h"
#include "../BinList/CCProtonPi0_BinList.h"

using namespace PlotUtils;

struct rootDir
{
    std::string data;
    std::string mc;
};

class CutArrow
{
    public:
        CutArrow(){ /* Do Nothing */ }
        CutArrow(double cut, std::string arrow_d )
        {
            cut_location = cut;
            arrow_direction = arrow_d;
        }

        double cut_location;
        std::string arrow_direction;
};

class CCProtonPi0_Plotter : public CCProtonPi0_NTupleAnalysis
{
    public:
        CCProtonPi0_Plotter();
        void plotHistograms();

    private:
        rootDir rootDir_GENIEXSec;
        rootDir rootDir_Truth;
        rootDir rootDir_CrossSection;
        rootDir rootDir_OtherStudies;
        rootDir rootDir_CutHists;
        rootDir rootDir_Interaction;
        rootDir rootDir_Muon;
        rootDir rootDir_Proton;
        rootDir rootDir_Pion;

        // POT Stats
        double POT_Ratio_data_mc;

        void setRootDirs();
        void getPOT_MC();
        void getPOT_Data();

        // Cross Section Plots
        void plotErrorSummary();
        void plotOriginalData();
        void plotBackgroundSubtracted();
        void plotBackgroundEstimated();
        void plotEfficiencyCorrected();
        void plotUnfolded();
        void plotIntegratedFlux();
        void plotCrossSection();
        void plotCrossSection_Check();

        // Data vs MC
        void plotInteraction_DataMC();
        void plotMuon_DataMC();
        void plotProton_DataMC();
        void plotPion_DataMC();
        void plotCutHistograms_DataMC();

        // MC Only
        void plotInteraction_MCOnly();
        void plotMuon_MCOnly();
        void plotProton_MCOnly();
        void plotPion_MCOnly();
        void plotCutHistograms_MCOnly();

        // True Signal Events
        void plotPion_True();

        // Other Plots 
        void plotGENIEXSec();
        void plotOtherStudies();
        void plot_InvMass_TruthMatch_Stacked(bool isSignal, bool isStacked);
        void plot_Michel_TruthMatch(std::string var);
        void plot_SignalKinematics_Stacked(std::string var, bool isStacked);
        void plot_final_mc_w_Stacked();
        void plot_stacked_pi0_P();
        void plot_stacked_pi0_theta();
        void plotStandardHistograms(rootDir &dir, std::string plotDir);

        // Helper Functions
        void ApplyStyle(MnvPlotter* plotter);
        void ApplyStyle_Errors(MnvPlotter* plotter);
        void AddCutArrow(MnvPlotter* plotter, CutArrow &cutArrow, double hist_max, double arrow_length);
        void SavePi0InvMassPoints();

        // --------------------------------------------------------------------
        // Plottting Macros - Implemented in CCProtonPi0_Plotter_Macros.cpp
        // --------------------------------------------------------------------
        // Default
        void DrawMnvH1D(rootDir& dir, std::string var_name, std::string plotDir);
        void Draw1DHist(rootDir &dir, std::string var_name, std::string plotDir, bool isLogScale = false);
        void Draw1DHist_Threshold(rootDir &dir, std::string var_name, std::string plotDir, double threshold = 0, bool isLogScale = false);
        void Draw2DHist(rootDir& dir, std::string var_name, std::string plotDir, double threshold = 0);
        void Draw3DHist(rootDir& dir, std::string var_name, std::string plotDir);

        // MC Only
        void DrawNormalizedMigrationHistogram(rootDir &dir, std::string var_name, std::string plotDir);
        void DrawMCWithErrorBand(rootDir& dir, std::string var_name, std::string plotDir);
        void DrawSignalMC(rootDir &dir, std::string var_name, std::string plotDir, int nCutArrows = 0, CutArrow cutArrow1 = CutArrow(), CutArrow cutArrow2 = CutArrow());
        void DrawErrorBand(MnvH1D* hist, std::string error_name, int error_color, std::string var_name, std::string plotDir);
        void DrawStackedMC(rootDir &dir, std::string var_name, std::string plotDir, int nCutArrows = 0, CutArrow cutArrow1 = CutArrow(), CutArrow cutArrow2 = CutArrow());
        void DrawStackedMC_BckgAll(rootDir &dir, std::string var_name, std::string plotDir, int nCutArrows = 0, CutArrow cutArrow1 = CutArrow(), CutArrow cutArrow2 = CutArrow());
        void DrawStackedMC_BckgType(rootDir &dir, std::string var_name, std::string plotDir, int nCutArrows = 0, CutArrow cutArrow1 = CutArrow(), CutArrow cutArrow2 = CutArrow());

        // Data vs MC
        void DrawDataMC(rootDir& dir, std::string var_name, std::string plotDir);
        void DrawDataMC(MnvH1D* data, MnvH1D* mc, std::string var_name, std::string plotDir);
        void DrawDataMC(MnvH1D* data, MnvH1D* mc, std::string var_name, std::string plotDir, bool isPOTNorm);
        void DrawDataMC_CrossSection(MnvH1D* data, MnvH1D* mc, std::string var_name, std::string plotDir);
        void DrawDataMC_CrossSection(MnvH1D* data, MnvH1D* mc, std::string var_name, std::string plotDir, bool isPOTNorm);
        void DrawDataMCRatio(MnvH1D* data, MnvH1D* mc, std::string var_name, std::string plotDir, bool isPOTNorm);
        void DrawDataMCRatio_CrossSection(MnvH1D* data, MnvH1D* mc, std::string var_name, std::string plotDir, bool isPOTNorm);
        void DrawDataStackedMC(rootDir &dir, std::string var_name, std::string plotDir, int nCutArrows = 0, CutArrow cutArrow1 = CutArrow(), CutArrow cutArrow2 = CutArrow());
        void DrawDataStackedMC_BckgAll(rootDir &dir, std::string var_name, std::string plotDir, bool isPOTNorm, int nCutArrows = 0, CutArrow cutArrow1 = CutArrow(), CutArrow cutArrow2 = CutArrow());
        void DrawDataStackedMC_BckgType(rootDir &dir, std::string var_name, std::string plotDir, bool isPOTNorm, int nCutArrows = 0, CutArrow cutArrow1 = CutArrow(), CutArrow cutArrow2 = CutArrow());

        // Other
        void DrawBackgroundSubtraction(bool isMC);
        void DrawErrorSummary(MnvH1D* hist, std::string var_name, std::string plotDir);
        void DrawTGraph(rootDir &dir, std::string var_name, std::string plotDir);
        void DrawEfficiencyCurve(rootDir& dir, std::string var_name, std::string plotDir);
        void DrawStackedMC_GammaEvis(rootDir &dir, int gammaID, std::string plotDir);
        void DrawStackedMC_GammaByPDG(rootDir &dir, std::string var_name, int gammaID, std::string plotDir);

        // Helper
        double FormTObjArray_BckgType(TFile* f_mc, std::string var_name, TObjArray* mc_hists, double &hist_max, double &bin_width); 
        void AddNormBox(MnvPlotter* plotter, bool isPOTNorm, double mc_ratio);
        void SaveRecoRatioPoints(rootDir& dir, std::string var_name, std::string plotDir);
        void Save2DHistPoints(rootDir& dir, std::string var_name, std::string plotDir);     
        double GetMCNormalization(std::string &norm_label, bool isPOTNorm, MnvH1D* data, MnvH1D* mc);
        void printBins(const MnvH1D* hist, const std::string var_name);

        // Unfolding Study
        void UnfoldingStudy();
        void UnfoldingStudy_muon_P();
        void UnfoldingStudy_muon_theta();
        void UnfoldingStudy_muon_cos_theta();
        void UnfoldingStudy_pi0_P();
        void UnfoldingStudy_pi0_KE();
        void UnfoldingStudy_pi0_theta();
        void init_UnfoldingHistograms(std::vector<MnvH1D*> &unfolded, std::vector<MnvH1D*> &error, std::vector<MnvH1D*> &diff);
        MnvH1D* CalcUnfoldingError(MnvH1D* diff, MnvH1D* truth);
        MnvH1D*  CalcUnfoldingDiff(MnvH1D* unfolded, MnvH1D* truth);
        void FillUnfoldingHistograms(MnvH1D* &unfolded, MnvH1D* &error, MnvH1D* &diff, MnvH2D* response, MnvH1D* mc_reco, MnvH1D* mc_true, int niter);
        void StyleUnfoldingHistograms(std::vector<MnvH1D*> &hists);
        void PlotUnfolding_StatErrors(const std::vector<MnvH1D*> &hists, const MnvH1D* truth, std::string var_name);
        void PlotUnfolding_Unfolded(const std::vector<MnvH1D*> &hists, const MnvH1D* truth, std::string var_name);
        void PlotUnfolding_Error(const std::vector<MnvH1D*> &hists, std::string var_name);
        void PlotUnfolding_Diff(const std::vector<MnvH1D*> &hists, std::string var_name);
        void PlotUnfolding_TruthComparison();
        void PlotUnfolding_Migration();
};

#endif

