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
    Last Revision: 2014_05_07
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

using namespace std;

class Plotter
{
    public:
   
        // -------------------------------------------------------------------------
        //     void plotHistograms: Generates plots for the output of run() function
        //         mcFile-> name of the .root file
        //         plotDir -> folder name for the plots will be created         
        //--------------------------------------------------------------------------
        void plotHistograms(bool isMC, bool isReco, bool is2D);
        void plotParticleInfo(  string rootDir, string plotDir, 
                                bool isMC, bool isReco, bool is2D);
                            
        void plotCCDeltaPlus(bool isMC, bool isReco, bool is2D);
        void plotMuon(bool isMC, bool isReco, bool is2D);
        void plotProton(bool isMC, bool isReco, bool is2D);
        void plotPion(bool isMC, bool isReco, bool is2D);
                            
        void inform(string rootDir, string plotDir);
        
        // Plottting Macros
        void plot1D_Hist(TH1F* hist1D, string fileName, string plotDir);
        void plot2D_Hist(TH2F* hist2D, string fileName, string plotDir);
    
    private:
    

};




#endif

