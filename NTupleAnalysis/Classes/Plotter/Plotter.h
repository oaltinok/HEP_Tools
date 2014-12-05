/*
================================================================================
Class: Plotter
    Plotter class includes specific functions for plotting 1D or 2D histograms
    
    Main Directory:
        Classes/Plotter/
    
    Usage:
        > #include "Classes/Plotter/Plotter.cpp" 
        > Plotter p;
        > p.plotHistograms(string mcFile, string plotDir)
            
    
    Author:        Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision: 2014_12_03
================================================================================
*/

#ifndef PLOTTER_H
#define PLOTTER_H

#include <iostream>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TFile.h>
#include <THStack.h>
#include <TPaveStats.h>
#include <TLegend.h>

#include "../../Libraries/Folder_List.h"

using namespace std;

const int nBranches = 3;

class Plotter
{
    public:
   
        // -------------------------------------------------------------------------
        //     void plotHistograms: Generates plots for the output of run() function
        //         mcFile-> name of the .root file
        //         plotDir -> folder name for the plots will be created         
        //--------------------------------------------------------------------------
        Plotter(int nMode);
        void setFolders();
        void plotHistograms(bool isMC, bool isReco, bool is2D);
        void inform(string rootDir, string plotDir);
        
        // Plottting Macros
        void plot1D_Hist(TH1D* hist1D, string fileName, string plotDir);
        void plot2D_Hist(TH2D* hist2D, string fileName, string plotDir);
        void plotStacked(TH1D* h_signal, TH1D* h_background, string plotName, string fileName, string plotDir);
        
        // Default Plots - File: Default_Plots.cpp
        void plotInteraction(bool isMC, bool isReco, bool is2D);
        void plotMuon(bool isMC, bool isReco, bool is2D);
        void plotProton(bool isMC, bool isReco, bool is2D);
        void plotPion(bool isMC, bool isReco, bool is2D);
        void plotParticleInfo(  string rootDir, string plotDir, 
                                bool isMC, bool isReco, bool is2D);
                            
               
        // Other Plots - File: Other_Plots.cpp
        void plotPID();
        void pID_proton();
        void pID_pion();
        void plot_2D_pID();
        void plot_mc_w_Stacked();
        void plotSignalBackground();
        void plotDebuggingPlots();

        
    
    private:
        int branchInd;
        string branchDir;
        string otherDir;
        
        string rootDir_Interaction[nBranches];
        string plotDir_Interaction[nBranches];
        
        string rootDir_Muon[nBranches];
        string plotDir_Muon[nBranches];
        
        string rootDir_Proton[nBranches];
        string plotDir_Proton[nBranches];
        
        string rootDir_Pion[nBranches];
        string plotDir_Pion[nBranches];
    
};




#endif

