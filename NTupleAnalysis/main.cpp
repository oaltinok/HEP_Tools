/*
================================================================================
main.cpp
    Main Function that controls the Analysis Package
    
    Classes Used:
        ANA_CCNuPionInc core class for the package
        Plotter class includes functions specific to generate Plots
    
    Input and Output folders assigned in this function
    Edit isTest variable to run a Test job or all playlist
    
    Terminal Usage:
        > root -l
        > .L main.cpp
            > main(false) for Running the Analysis
            > main(true) for Generating Plots only
    
    Last Revision: 2014_04_16
================================================================================
*/

// Include Required Classes & Libraries
#include "Classes/CCDeltaPlus/CCDeltaPlus.cpp"
#include "Classes/Plotter/Plotter.cpp"
#include "Libraries/Folder_List.h"

#include <string>

using namespace std;


void main(bool onlyPlot)
{
    // Edit isTest Variable for running Test Samples or complete playlist
    bool isTest = true;
    
    // Plot Selection
    bool plotsMC    = true;
    bool plotsReco  = false;
    bool plots2D    = false;
    
    if(onlyPlot){
        Plotter p;
        p.plotHistograms(plotsMC,plotsReco,plots2D); 
    }else{
        CCDeltaPlus t;
        if(isTest){
            t.run(f_PL_Test);
        }else{
            t.run(f_PL_Full);
        }
    }
    
}



