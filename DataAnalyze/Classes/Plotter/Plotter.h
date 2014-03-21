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
            
    
    Last Revision: 2014_03_20
================================================================================
*/

#ifndef PLOTTER_H
#define PLOTTER_H

using namespace std;

#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "Libraries/Folder_List.h" // Contains File and Folder Paths


class Plotter
{
    public:
    
        Plotter();
        
        // -------------------------------------------------------------------------
        //     void plotHistograms: Generates plots for the output of run() function
        //         mcFile-> name of the .root file
        //         plotDir -> folder name for the plots will be created         
        //--------------------------------------------------------------------------
        void plotHistograms(string mcFile, string plotDir);
        void plotDefault(  TFile* f_mc, string plotDir, 
                            bool isMC, bool isReco, bool is2D);
                            
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

