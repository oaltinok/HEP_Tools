/*
================================================================================
main.cpp
    Main Function that controls the Analysis Package
    
    Classes Used:
        CDeltaPlus core class for the package
        Plotter class includes functions specific to generate Plots
    
    Input and Output folders assigned in this function
    Edit isTest variable to run a Test job or all playlist
    
    Terminal Usage:
        > root -l
        > .L main.cpp
            > main(false) for Running the Analysis
            > main(true) for Generating Plots only
    
    Author:        Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision: 2014_05_13
================================================================================
*/

// Include Required Classes
#include "Classes/CCProtonPi0/CCProtonPi0.h"
#include "Classes/Plotter/Plotter.h"

#include <string>

const    string option1 = "run";
const    string option2 = "plot";

void showInputError(char *argv[]);

using namespace std;

int main(int argc, char *argv[] )
{
    bool onlyPlot;
    string input;
    
    if ( argc != 2 || argc > 2) /* argc should be 2 for correct execution */
    {
        /* We print argv[0] assuming it is the program name */
        showInputError(argv);
        return 0;
    }else{
        input = argv[1];
    }
    
    if(input.compare(option1) == 0){
        cout<<"Run Option Selected!"<<endl;
        onlyPlot = false;
    }
    else if(input.compare(option2) == 0){
        cout<<"Plot Option Selected!"<<endl;
        onlyPlot = true;
    }
    else{ 
        showInputError(argv);
        return 0;
    }
    

    // Edit isTest Variable for running Test Samples or complete playlist
    bool isTest = true;
    
    // Plot Selection
    bool plotsMC    = false;
    bool plotsReco  = false;
    bool plots2D    = true;
    
    if(onlyPlot){
        Plotter p;
        p.plotHistograms(plotsMC,plotsReco,plots2D); 
    }else{
        CCProtonPi0 t;
        if(isTest){
            t.run("Input/Playlists/pl_MC_Test.dat");
        }else{
            t.run("Input/Playlists/pl_MC_minerva1.dat");
        }
    }
    
    return 0;
}

void showInputError(char *argv[])
{
    cout<<"Not a valid option!"<<endl;
    cout<<"Run: "<<argv[0]<<" "<<option1<<endl;
    cout<<"Plot: "<<argv[0]<<" "<<option1<<endl;
}



