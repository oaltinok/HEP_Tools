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

#include <iostream>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TFile.h>
#include <THStack.h>
#include <TPaveStats.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TVectorD.h>
#include <PlotUtils/MnvPlotter.h>
#include <PlotUtils/MnvH1D.h>
#include <PlotUtils/MnvH2D.h>
#include <PlotUtils/MnvFluxConstraint.h>
#include "Cintex/Cintex.h"

#include "../../Libraries/Folder_List.h"

using namespace PlotUtils;

struct rootDir
{
    std::string mc_signal;
    std::string mc_background;
    std::string mc_all;
    std::string data;
};

const int nBranches = 3;

class CCProtonPi0_Plotter
{
    public:
   
        // -------------------------------------------------------------------------
        //     void plotHistograms: Generates plots for the output of run() function
        //         mcFile-> name of the .root file
        //         plotDir -> folder name for the plots will be created         
        //--------------------------------------------------------------------------
        CCProtonPi0_Plotter(int nMode, bool isMC);
        void plotHistograms();
        
    private:
        bool m_isMC;
        bool isSignalvsBackground;
        int branchInd;
        std::string branchDir;
        std::string otherDir;

        rootDir rootDir_Muon;
        rootDir rootDir_Proton;
        rootDir rootDir_Pion;

        std::string rootDir_Interaction[nBranches];
        std::string plotDir_Interaction[nBranches];
        
        std::string rootDir_PID[nBranches];
        std::string plotDir_PID[nBranches];
        
        std::string plotDir_Muon[nBranches];
        std::string plotDir_Proton[nBranches];
        std::string plotDir_Pion[nBranches];
        
        std::string rootDir_Pi0Blob[nBranches];
        std::string plotDir_Pi0Blob[nBranches];
  
        void inform(std::string rootDir_mc, std::string rootDir_data, std::string plotDir);
        void inform(std::string rootDir, std::string plotDir);
        void setRootDirs(rootDir& dirs, std::string fileName );
        void setPlotDirs();
       
        // Default Plots - File: Default_Plots.cpp
        void plotInteraction();
        void plotMuon();
        void plotProton();
        void plotPion();
        
        // pID Plots - File: pID_Plots.cpp
        void plotPID();
        void pID_proton();
        void pID_proton_LLR();
        void plot_2D_pID();
        void pIDDiff();
        void pIDStats();
        void KE();
        
        // Other Plots - File: Other_Plots.cpp
        void plot_mc_w_Stacked();
        void plot_final_mc_w_Stacked();
        void plotSignalBackground();
        void plotSignalBackground_Pi0Blob();
        void plotCutHistograms();
        void plotMichel();
        void MichelTool(TH1D* vertex, TH1D* track, TH1D* track2, TH1D* missed, std::string plotName, std::string fileName, std::string plotDir);
        void plotStandardHistograms(rootDir &dir, std::string plotDir);

         // Plottting Macros
        void plot1D_Hist(TH1D* hist1D, std::string fileName, std::string plotDir);
        void plot2D_Hist(TH2D* hist2D, std::string fileName, std::string plotDir);
        void plot1D_HistLogScale(TH1D* hist1D, std::string fileName, std::string plotDir);
        void DrawDataMC(rootDir& dir, std::string var_name, std::string plotDir);
        void plotStacked(TH1D* h_signal, TH1D* h_background, 
                            std::string plotName, std::string fileName, std::string plotDir, 
                            std::string signal_label = "Signal", std::string background_label = "Background",
                            bool isRatioReversed = false);
        void plotStackedLogScale(TH1D* h_signal, TH1D* h_background, std::string plotName, std::string fileName, std::string plotDir);
        void plotSignalRatio(TH1D* h_signal, TH1D* h_background, std::string fileName, std::string plotDir, bool isReversed = false);
        void plot_purity_efficiency(TH1D* h_signal, TH1D* h_background, std::string fileName, std::string plotDir, bool keepEventstoRight = true);
        

};




#endif

