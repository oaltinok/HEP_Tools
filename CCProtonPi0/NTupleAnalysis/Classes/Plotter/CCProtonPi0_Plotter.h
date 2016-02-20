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
#include "Cintex/Cintex.h"

#include "../NTupleAnalysis/CCProtonPi0_NTupleAnalysis.h"

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
        rootDir rootDir_Truth;
        rootDir rootDir_CrossSection;
        rootDir rootDir_OtherStudies;
        rootDir rootDir_CutHists;
        rootDir rootDir_Interaction;
        rootDir rootDir_Muon;
        rootDir rootDir_Proton;
        rootDir rootDir_Pion;
        rootDir rootDir_Pi0Blob;

        // POT Stats
        double data_POT;
        double mc_POT;
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
        void plotCrossSection();
        void plotCrossSection_Check();

        // Data vs MC
        void plotInteraction_DataMC();
        void plotMuon_DataMC();
        void plotProton_DataMC();
        void plotPion_DataMC();
        void plotPi0Blob_DataMC();
        void plotCutHistograms_DataMC();

        // MC Only
        void plotInteraction_MCOnly();
        void plotMuon_MCOnly();
        void plotProton_MCOnly();
        void plotPion_MCOnly();
        void plotPi0Blob_MCOnly();
        void plotCutHistograms_MCOnly();

        // True Signal Events
        void plotPion_True();

        // Other Plots 
        void plotOtherStudies();
        void plot_mc_w_Stacked();
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
        void DrawStackedMC_BckgWithPi0(rootDir &dir, std::string var_name, std::string plotDir, int nCutArrows = 0, CutArrow cutArrow1 = CutArrow(), CutArrow cutArrow2 = CutArrow());

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
        void DrawDataStackedMC_BckgWithPi0(rootDir &dir, std::string var_name, std::string plotDir, bool isPOTNorm, int nCutArrows = 0, CutArrow cutArrow1 = CutArrow(), CutArrow cutArrow2 = CutArrow());

        // Other
        void DrawErrorSummary(MnvH1D* hist, std::string var_name, std::string plotDir);
        void DrawTGraph(rootDir &dir, std::string var_name, std::string plotDir);
        void DrawEfficiencyCurve(rootDir& dir, std::string var_name, std::string plotDir);
        void DrawStackedMC_GammaEvis(rootDir &dir, int gammaID, std::string plotDir);
        void DrawStackedMC_GammaByPDG(rootDir &dir, std::string var_name, int gammaID, std::string plotDir);

        // Helper
        void FormTObjArray_BckgType(TFile* f_mc, std::string var_name, TObjArray* mc_hists, double &hist_max, double &bin_width); 
        void FormTObjArray_BckgWithPi0(TFile* f_mc, std::string var_name, TObjArray* mc_hists, double &hist_max, double &bin_width); 
        void AddNormBox(MnvPlotter* plotter, bool isPOTNorm, double mc_ratio);
        void SaveRecoRatioPoints(rootDir& dir, std::string var_name, std::string plotDir);
        void Save2DHistPoints(rootDir& dir, std::string var_name, std::string plotDir);     
        double GetMCNormalization(std::string &norm_label, bool isPOTNorm, MnvH1D* data, MnvH1D* mc);
};

#endif

