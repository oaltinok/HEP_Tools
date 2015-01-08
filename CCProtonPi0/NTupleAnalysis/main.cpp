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
//     Last Revision: 2015_01_08
================================================================================
*/

// Include Required Classes
#include "Classes/CCProtonPi0/CCProtonPi0.h"
#include "Classes/Plotter/Plotter.h"

#include <string>
#include <ctime>

const string runOption1 = "run";
const string runOption2 = "plot";
const string modeOption1 = "1";
const string modeOption2 = "2";
const string modeOption3 = "3";
const string modeOption4 = "4";

void showInputError(char *argv[]);
bool get_runSelection(string runSelect, bool &onlyPlot);
bool get_modeSelection(string modeSelect, int &nMode);

using namespace std;

int main(int argc, char *argv[] )
{

    clock_t timeStart;
    timeStart = clock();
    clock_t timeEnd;
    double timeDiff;
    int timeDiff_m;
    int timeDiff_s;
    int nMode;
    bool onlyPlot;
    string runSelect;
    string modeSelect;
    bool isRunSelected;
    bool isModeSelected;
    
    // Edit isTest Variable for running Test Samples or complete playlist
    bool isTest     = false;
    bool isPlaylist = false;
    
    // Plot Selection
    bool plotsMC    = true;
    bool plotsReco  = true;
    bool plots2D    = true;

    // Check User Input
    if ( argc != 3 || argc > 3) /* argc should be 3 for correct execution */
    {
        /* We print argv[0] assuming it is the program name */
        showInputError(argv);
        return 0;
    }else{
        runSelect = argv[1];
        modeSelect = argv[2];
    }
    
    
    // Get Run Selection   
    isRunSelected = get_runSelection(runSelect, onlyPlot);
    if( !isRunSelected ){
        showInputError(argv);
        return 0;
    }
    
    // Get Mode Selection   
    isModeSelected = get_modeSelection(modeSelect, nMode);
    if( !isModeSelected ){
        showInputError(argv);
        return 0;
    }
    
    if(onlyPlot){
        Plotter p(nMode);
        p.plotHistograms(plotsMC,plotsReco,plots2D); 
    }else{
        CCProtonPi0 t;
        t.setAnalysisMode(nMode);
        if(isTest){
            t.setChannelTag("test");
            t.run("Input/Playlists/pl_MC_Test_Sample.dat");
//             t.run("Input/Playlists/pl_MC_Test_Sample_nogrid.dat");
        }else if(isPlaylist){
            t.setChannelTag("minerva13C");
            t.run("Input/Playlists/pl_MC_minerva13C.dat");
        }else{
            t.setChannelTag("runset");
            t.run("Input/Playlists/pl_MC_Run.dat");
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
    cout<<std::left;
    cout<<"----------------------------------------------------------------------"<<endl;
    cout<<"Not a valid syntax!"<<endl;
    cout<<"----------------------------------------------------------------------"<<endl;
    cout<<"Correct Syntax: "<<endl;
    cout<<"\t"<<argv[0]<<" runOption modeOption"<<"\n"<<endl;

    cout<<"runOptions:"<<endl;
    cout.width(16);cout<<"\tRun:"<<runOption1<<endl;
    cout.width(16);cout<<"\tPlot:"<<runOption2<<"\n"<<endl;
    
    cout<<"modeOptions:"<<endl;
    cout.width(16); cout<<"\tSignal:"<<modeOption1<<endl;
    cout.width(16); cout<<"\tBackground:"<<modeOption2<<endl;
    cout.width(16); cout<<"\tEverything:"<<modeOption3<<endl;
    cout<<"----------------------------------------------------------------------"<<endl;
}

bool get_runSelection(string runSelect, bool &onlyPlot)
{
    // Get Run Selection
    if(runSelect.compare(runOption1) == 0){
        cout<<"Run Option Selected!"<<endl;
        onlyPlot = false;
    }
    else if(runSelect.compare(runOption2) == 0){
        cout<<"Plot Option Selected!"<<endl;
        onlyPlot = true;
    }
    else{ 
        return false;
    }
    
    return true;
}

bool get_modeSelection(string modeSelect, int &nMode)
{
    // Get Mode Selection
    if(modeSelect.compare(modeOption1) == 0){
        cout<<"Signal Mode Selected!"<<endl;
        nMode = 1;
    }
    else if(modeSelect.compare(modeOption2) == 0){
        cout<<"Background Mode Selected!"<<endl;
        nMode = 2;
    }
    else if(modeSelect.compare(modeOption3) == 0){
        cout<<"All Events Mode Selected!"<<endl;
        nMode = 3;
    }else if (modeSelect.compare(modeOption4) == 0){
        cout<<"Plot Only! Signal vs Background Mode"<<endl;
        nMode = 4;  
    }else{
        return false;   
    }
    
    return true;
}


