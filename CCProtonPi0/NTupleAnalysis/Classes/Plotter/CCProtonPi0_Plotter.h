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
#include <PlotUtils/TargetUtils.h>
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
        bool plot_muon_P;
        bool plot_muon_theta;
        bool plot_pi0_P;
        bool plot_pi0_KE;
        bool plot_pi0_theta;
        bool plot_QSq;
        bool plot_W;
        bool plot_Enu;

        void PlotXSecVar(std::string var_name, std::string data_var, std::string mc_var, std::string plotDir, std::string plotName);
        void plotCrossSection();
        void plotOriginalData();
        void plotBackgroundSubtracted();
        void plotBackgroundEstimated();
        void plotEfficiencyCorrected();
        void plotUnfolded();
        void plotFluxIntegrated();
        void plotXSec();
        void plotCrossSection_Check();
        void plotCrossSection_Check(std::string var_name, std::string plotDir);

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
        void plotTruth_Pion();
        void plotTruth_Enu();
        void plotTruth_QSq();
        void plotTruth_W();
        void plotTruth_ShortProton();

        // Other Plots 
        void PlotTotalEnuXSec();
        void plotGENIEXSec();
        void plotOtherStudies();
        void plot_InvMass_TruthMatch_Stacked(bool isSignal, bool isStacked);
        void plot_Michel_TruthMatch(std::string var);
        void plot_SignalKinematics(std::string var, std::string type, bool isStacked);
        void plot_SignalKinematics();
        void plot_stacked_pi0_P();
        void plot_stacked_pi0_theta();
        void plotStandardHistograms(rootDir &dir, std::string plotDir);

        // Helper Functions
        void ApplyStyle(MnvPlotter* plotter);
        void ApplyStyle_Errors(MnvPlotter* plotter, bool groupGENIE);
        void AddCutArrow(MnvPlotter* plotter, CutArrow &cutArrow, double hist_max, double arrow_length);
        void SavePi0InvMassPoints();
        void NormalizeToNormBinWidth(MnvH1D* hist);
        double GetSmallestBinWidth(MnvH1D* hist);
        
        // Flux Study
        double GetFluxHistContent(MnvH1D* hist, double low1, double low2);
        double Integrate(MnvH1D* hist, int start, int end);
        void GetFlux();
        void PlotFluxHistograms();
        void PlotFluxComparison(std::string plotDir);
        void PlotFluxRatio(std::string plotDir);
        void PlotFluxRebinned(std::string plotDir);
        void PlotFluxRebinned(MnvH1D* original, MnvH1D* rebinned, std::string plotDir);

        // --------------------------------------------------------------------
        // Plottting Macros - Implemented in CCProtonPi0_Plotter_Macros.cpp
        // --------------------------------------------------------------------
        // Default
        void DrawMnvH1D(rootDir& dir, std::string var_name, std::string plotDir);
        void DrawMnvH1D(MnvH1D* hist1D, std::string var_name, std::string plotDir);
        void Draw1DHist(TH1* hist1D, std::string var_name, std::string plotDir, bool isLogScale = false);

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
        void DrawDataMC(MnvH1D* data, MnvH1D* mc, std::string var_name, std::string plotDir, bool isXSec = false);
        void DrawDataMC_WithRatio(MnvH1D* data, MnvH1D* mc, std::string var_name, std::string plotDir, bool isPOTNorm, bool isXSec = false);
        void DrawDataMCRatio(MnvH1D* data, MnvH1D* mc, std::string var_name, std::string plotDir, bool isPOTNorm, bool isXSec = false);
        void DrawDataStackedMC(rootDir &dir, std::string var_name, std::string plotDir, int nCutArrows = 0, CutArrow cutArrow1 = CutArrow(), CutArrow cutArrow2 = CutArrow());
        void DrawDataStackedMC_BckgAll(rootDir &dir, std::string var_name, std::string plotDir, bool isPOTNorm, int nCutArrows = 0, CutArrow cutArrow1 = CutArrow(), CutArrow cutArrow2 = CutArrow());
        void DrawDataStackedMC_BckgType(rootDir &dir, std::string var_name, std::string plotDir, bool isPOTNorm, int nCutArrows = 0, CutArrow cutArrow1 = CutArrow(), CutArrow cutArrow2 = CutArrow());

        // Other
        void PlotDelta();
        void DrawBackgroundSubtraction(bool isMC);
        void DrawErrorSummary(MnvH1D* hist, std::string var_name, std::string plotDir, bool groupGENIE = true);
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

        // Systematics
        void Systematics();
        void Systematics_Practice();
        void Systematics_Practice(std::string root_dir, std::string var_name, std::string err_name);
        void Systematics_RawData();
        void Systematics_DrawErrorSummary(std::string data_var, std::string mc_var);
        void Systematics_DrawErrorBand_GENIE(std::string mc_var);
        void Systematics_DrawErrorSummary_GENIE(MnvH1D* hist, std::string var_name, std::string error_name, double err_genie_total);
        void Systematics_WriteTable(MnvH1D* hist, std::string var_name);

};

#endif

