/*
================================================================================
Class: CCProtonPi0_Plotter
    CCProtonPi0_Plotter class includes specific functions for plotting 1D or 2D histograms
    
    Main Directory:
        Classes/Plotter/
    
    Usage:
        > #include "Classes/Plotter/CCProtonPi0_Plotter.cpp" 
        > CCProtonPi0_Plotter p;
        > p.plotHistograms(string mcFile, string plotDir)
            
    
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

#include "../../Libraries/Folder_List.h"

using namespace std;

const int nBranches = 3;

class CCProtonPi0_Plotter
{
    public:
   
        // -------------------------------------------------------------------------
        //     void plotHistograms: Generates plots for the output of run() function
        //         mcFile-> name of the .root file
        //         plotDir -> folder name for the plots will be created         
        //--------------------------------------------------------------------------
        CCProtonPi0_Plotter(int nMode);
        void setFolders();
        void plotHistograms();
        void inform(string rootDir, string plotDir);
        
        // Plottting Macros
        void plot1D_Hist(TH1D* hist1D, string fileName, string plotDir);
        void plot2D_Hist(TH2D* hist2D, string fileName, string plotDir);
        void plot1D_HistLogScale(TH1D* hist1D, string fileName, string plotDir);
        void plotStacked(TH1D* h_signal, TH1D* h_background, 
                            string plotName, string fileName, string plotDir, 
                            string signal_label = "Signal", string background_label = "Background",
                            bool isRatioReversed = false);
        void plotStackedLogScale(TH1D* h_signal, TH1D* h_background, string plotName, string fileName, string plotDir);
        void plotSignalRatio(TH1D* h_signal, TH1D* h_background, string fileName, string plotDir, bool isReversed = false);
        void plot_purity_efficiency(TH1D* h_signal, TH1D* h_background, string fileName, string plotDir, bool keepEventstoRight = true);
        
         // Default Plots - File: Default_Plots.cpp
        void plotInteraction();
        void plotMuon();
        void plotProton();
        void plotPion();
        void plotParticleInfo(  string rootDir, string plotDir);
                                
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
        void plotCutHistograms();
        void plotMichel();
        void MichelTool(TH1D* vertex, TH1D* track, TH1D* track2, TH1D* missed,
                         string plotName, string fileName, string plotDir);
        
    
    private:
        bool isSignalvsBackground;
        int branchInd;
        string branchDir;
        string otherDir;
        
        string rootDir_Interaction[nBranches];
        string plotDir_Interaction[nBranches];
        
        string rootDir_PID[nBranches];
        string plotDir_PID[nBranches];
        
        string rootDir_Muon[nBranches];
        string plotDir_Muon[nBranches];
        
        string rootDir_Proton[nBranches];
        string plotDir_Proton[nBranches];
        
        string rootDir_Pion[nBranches];
        string plotDir_Pion[nBranches];
    
};




#endif

