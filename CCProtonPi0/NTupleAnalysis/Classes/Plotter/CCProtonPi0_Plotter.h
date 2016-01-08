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

#include <fstream>
#include <string>
#include <iostream>
#include <TVectorD.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TPaveStats.h>
#include <TPad.h>
#include <TMath.h>
#include <TLegend.h>
#include <THStack.h>
#include <TH2.h>
#include <TH1.h>
#include <TGraph.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <PlotUtils/MnvPlotter.h>
#include <PlotUtils/MnvH2D.h>
#include <PlotUtils/MnvH1D.h>
#include <PlotUtils/MnvFluxConstraint.h>
#include <PlotUtils/POTCounter.h>
#include "Cintex/Cintex.h"

#include "../../Libraries/Folder_List.h"

using namespace PlotUtils;

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

struct rootDir
{
    std::string data;
    std::string mc;
};

class CCProtonPi0_Plotter
{
    public:
        CCProtonPi0_Plotter(std::string ana_folder);
        void plotHistograms();
        
    private:
        // POT Stats
        double data_POT;
        double mc_POT;
        double POT_Ratio_data_mc;
        
        rootDir rootDir_OtherStudies;
        rootDir rootDir_CutHists;
        rootDir rootDir_Interaction;
        rootDir rootDir_Muon;
        rootDir rootDir_Proton;
        rootDir rootDir_Pion;
        rootDir rootDir_Pi0Blob;
        
        std::string plotDir_OtherStudies;
        std::string plotDir_CutHists;
        std::string plotDir_Interaction;
        std::string plotDir_Muon;
        std::string plotDir_Proton;
        std::string plotDir_Pion;
        std::string plotDir_Pi0Blob;
        std::string plotDir_Other;
        
        void setRootDirs(std::string ana_folder);
        void setPlotDirs(std::string ana_folder);
        void getPOT_MC();
        void getPOT_Data();

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
        void plotEfficiencyCurves();
      
        // True Signal Events
        void plotPion_True();

        // Other Plots 
        void plotOtherStudies();
        void plot_mc_w_Stacked();
        void plot_final_mc_w_Stacked();
        void plotStandardHistograms(rootDir &dir, std::string plotDir);

        // Helper Functions
        void ApplyStyle(MnvPlotter* plotter);
        void AddCutArrow(MnvPlotter* plotter, CutArrow &cutArrow, double hist_max, double arrow_length);
        void SavePi0InvMassPoints();
        
        // Plottting Macros
        void Draw1DHist(rootDir &dir, std::string var_name, std::string plotDir, bool isLogScale = false);
        void Draw1DHist_Threshold(rootDir &dir, std::string var_name, std::string plotDir, double threshold = 0, bool isLogScale = false);
        void Draw2DHist(rootDir& dir, std::string var_name, std::string plotDir, double threshold = 0);
        void DrawMCWithErrorBand(rootDir& dir, std::string var_name, std::string plotDir);
        void DrawStackedMC_BckgAll(rootDir &dir, std::string var_name, std::string plotDir, int nCutArrows = 0, CutArrow cutArrow1 = CutArrow(), CutArrow cutArrow2 = CutArrow());
        void DrawDataStackedMC(rootDir &dir, std::string var_name, std::string plotDir, int nCutArrows = 0, CutArrow cutArrow1 = CutArrow(), CutArrow cutArrow2 = CutArrow());
        void DrawDataStackedMC_BckgAll(rootDir &dir, std::string var_name, std::string plotDir, bool isPOTNorm, int nCutArrows = 0, CutArrow cutArrow1 = CutArrow(), CutArrow cutArrow2 = CutArrow());
        void DrawDataStackedMC_BckgWithPi0(rootDir &dir, std::string var_name, std::string plotDir, bool isPOTNorm, int nCutArrows = 0, CutArrow cutArrow1 = CutArrow(), CutArrow cutArrow2 = CutArrow());
        void DrawDataStackedMC_BckgType(rootDir &dir, std::string var_name, std::string plotDir, bool isPOTNorm, int nCutArrows = 0, CutArrow cutArrow1 = CutArrow(), CutArrow cutArrow2 = CutArrow());
        void DrawDataMC(rootDir& dir, std::string var_name, std::string plotDir);
        void DrawDataMC(rootDir& dir, std::string var_name, std::string plotDir, bool isPOTNorm);
        void DrawDataMCRatio(rootDir& dir, std::string var_name, std::string plotDir, bool isPOTNorm);
        void DrawStackedMC_GammaEvis(rootDir &dir, int gammaID, std::string plotDir);
        void DrawStackedMC_GammaByPDG(rootDir &dir, std::string var_name, int gammaID, std::string plotDir);
        void AddNormBox(MnvPlotter* plotter, bool isPOTNorm, double mc_ratio);
        void SaveRecoRatioPoints(rootDir& dir, std::string var_name, std::string plotDir);
        void Save2DHistPoints(rootDir& dir, std::string var_name, std::string plotDir);     
        double GetMCNormalization(std::string &norm_label, bool isPOTNorm, MnvH1D* data, MnvH1D* mc);
        void DrawEfficiencyCurve(std::string var_name, std::string plotDir, TH1D* all_signal, TH1D* signal);


};




#endif
