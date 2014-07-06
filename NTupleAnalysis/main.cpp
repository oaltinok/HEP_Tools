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
    Last Revision: 2014_07_06
================================================================================
*/

// Include Required Classes
#include "Classes/CCProtonPi0/CCProtonPi0.h"
#include "Classes/Plotter/Plotter.h"

#include <string>
#include <ctime>

const    string option1 = "run";
const    string option2 = "plot";

void showInputError(char *argv[]);

using namespace std;

int main(int argc, char *argv[] )
{
    clock_t timeStart;
    timeStart = clock();
    clock_t timeEnd;
    double timeDiff;
    int timeDiff_m;
    int timeDiff_s;
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
    bool isTest     = true;
    bool isComplete = true;
    
    // Plot Selection
    bool plotsMC    = true;
    bool plotsReco  = true;
    bool plots2D    = true;
    
    if(onlyPlot){
        Plotter p;
        p.plotHistograms(plotsMC,plotsReco,plots2D); 
    }else{
        CCProtonPi0 t;
        if(isTest){
//             t.run("Input/Playlists/pl_MC_Test_Sample.dat");
            t.run("Input/Playlists/pl_MC_Test.dat");
        }else if(isComplete){
            t.run("Input/Playlists/MC_minerva1_CCProtonPi0.dat");
        }else{
            t.run("Input/Playlists/pl_MC_Test.dat");
        }
    }
    
    
    timeEnd = clock();
    timeDiff = ( timeEnd - timeStart ) / (double)CLOCKS_PER_SEC;
    

    timeDiff_m = (int)timeDiff / 60;
    timeDiff_s = (int)timeDiff % 60;
    
    cout<<"Time to Complete = "<<timeDiff_m<<"\' "<<timeDiff_s<<"\" or "<<timeDiff<<" seconds"<<endl;
    return 0;
}

void showInputError(char *argv[])
{
    cout<<"Not a valid option!"<<endl;
    cout<<"Run: "<<argv[0]<<" "<<option1<<endl;
    cout<<"Plot: "<<argv[0]<<" "<<option2<<endl;
}



