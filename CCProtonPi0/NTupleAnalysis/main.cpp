/*
================================================================================
main.cpp
    Main Function that controls the Analysis Package
    
    Classes Used:
        CCProtonPi0_Analyzer core class for the package
        CCProtonPi0_Plotter class includes functions specific to generate Plots
    
    Input and Output folders assigned in this function
    Edit isTest variable to run a Test job or all playlist
    
    Build Package using "make"
    
    Usage:
    Reduce NTuple to a Single File
            > ./main.exe reduce

    Analyze NTuple Files or Reduced File
        Signal Events
            > ./main.exe run 1
            > ./main.exe plot 1
        Background Events
            > ./main.exe run 2
            > ./main.exe plot 2
        All Events
            > ./main.exe run 3
            > ./main.exe plot 3
        Plot Only: Signal vs Background - Requires both signal and background
            > ./main.exe plot 4
        
    Author:        Ozgur Altinok  - ozgur.altinok@tufts.edu
================================================================================
*/

// Include Required Classes
#include "Classes/Analyzer/CCProtonPi0_Analyzer.h"
#include "Classes/Plotter/CCProtonPi0_Plotter.h"
#include "Cintex/Cintex.h"

#include <string>
#include <ctime>

using namespace std;

bool isMC = false;

const string runOption_Run = "run";
const string runOption_Plot = "plot";
const string runOption_Reduce = "reduce";
const string modeOption_Signal = "1";
const string modeOption_Background = "2";
const string modeOption_All = "3";
const string modeOption_Other = "4";

void showInputError(char *argv[]);
bool get_runSelection(string runSelect, bool &onlyPlot);
bool get_modeSelection(string modeSelect, int &nMode);
void GetAnalysisMode( char* argv[], int &nMode, bool &onlyPlot );
int GetMode(string modeSelect);
bool isModeReduce(int argc, char* argv[]);
bool isModePlot(string runSelect);
bool isValidCommand(int argc, char* argv[]);
bool isValidCommand_Reduce(string runSelect);
bool isValidCommand_Analyze(string runSelect, string modeSelect);
void Plot(int nMode);
void Reduce(string playlist);
void Analyze(string playlist, int nMode);

int main(int argc, char *argv[] )
{
    time_t timeStart; time(&timeStart);
    time_t timeEnd;
    double timeDiff;
    int timeDiff_m;
    int timeDiff_s;
    string pl_reduce;
    string pl_analyze;

    // Check User Command
    if ( !isValidCommand(argc, argv)){
        showInputError(argv);
        return 0;
    }
    
    if (isMC){
        cout<<"MC Playlists Selected!"<<endl;
        pl_reduce = "Input/Playlists/pl_MC_Merged.dat"; 
        pl_analyze = "Input/Playlists/pl_MC_Reduced.dat"; 
    }else{
        cout<<"Data Playlists Selected!"<<endl;
        pl_reduce = "Input/Playlists/pl_Data_Merged.dat"; 
        pl_analyze = "Input/Playlists/pl_Data_Reduced.dat"; 
    }
   
    ROOT::Cintex::Cintex::Enable();
  
    if (isModeReduce(argc,argv)){
        Reduce(pl_reduce);
    }else{
        int nMode;
        bool onlyPlot;
        GetAnalysisMode(argv, nMode, onlyPlot);
   
        if(onlyPlot){
            Plot(nMode);
        }else{
            Analyze(pl_analyze, nMode);
        }
    }

    time(&timeEnd);
    timeDiff = ( timeEnd - timeStart );
    
    timeDiff_m = (int)timeDiff / 60;
    timeDiff_s = (int)timeDiff % 60;
    
    cout<<"Time to Complete = "<<timeDiff_m<<"\' "<<timeDiff_s<<"\" or "<<timeDiff<<" seconds"<<endl;
    return 0;
}


void Reduce(string playlist)
{
    // nMode == 0 for Reduce Mode
    int nMode = 0;
    CCProtonPi0_Analyzer t(nMode,isMC);
    t.reduce(playlist);
}

void Plot(int nMode)
{
    CCProtonPi0_Plotter p(nMode, isMC);
    p.plotHistograms();
}

void Analyze(string playlist, int nMode)
{
    CCProtonPi0_Analyzer t(nMode, isMC);
    t.analyze(playlist);
}

void GetAnalysisMode( char* argv[], int &nMode, bool &onlyPlot )
{
    string runSelect = argv[1];
    string modeSelect = argv[2];

    onlyPlot = isModePlot(runSelect);
    nMode = GetMode(modeSelect);
}

bool isModeReduce(int argc, char* argv[])
{
    if (argc == 3) return false;

    string runSelect = argv[1];
    if (runSelect.compare(runOption_Reduce) == 0){ 
        cout<<"Reduce Mode Selected!"<<endl;
        return true;
    }else{ 
        return false;

    }
}

bool isValidCommand(int argc, char* argv[])
{
    string runSelect;
    string modeSelect;

    // argc can only be 2 or 3
    if (argc != 2 && argc != 3) return false;

    // Check for Reduce Command
    if (argc == 2){
        runSelect = argv[1];
        return isValidCommand_Reduce(runSelect);
    }

    // Check for Analyze Command
    if (argc == 3) {
        runSelect = argv[1];
        modeSelect = argv[2];
        return isValidCommand_Analyze(runSelect,modeSelect);
    }

    return false;
}

bool isValidCommand_Analyze(string runSelect, string modeSelect)
{
    if (runSelect.compare(runOption_Run) == 0){
        if( modeSelect.compare(modeOption_All) == 0 || modeSelect.compare(modeOption_Signal) == 0 || modeSelect.compare(modeOption_Background) == 0 ){
                return true;
        }else{
                return false;
        }
    }else if (runSelect.compare(runOption_Plot) == 0){
        if( modeSelect.compare(modeOption_All) == 0 || modeSelect.compare(modeOption_Signal) == 0 || modeSelect.compare(modeOption_Background) == 0 || modeSelect.compare(modeOption_Other) == 0){
                return true;
        }else{
                return false;
        }
    }else{
        return false;
    }
}

bool isValidCommand_Reduce(string runSelect)
{
    if (runSelect.compare(runOption_Reduce) == 0) return true;
    else return false;
}

void showInputError(char *argv[])
{
    cout<<std::left;
    cout<<"----------------------------------------------------------------------"<<endl;
    cout<<"Not a valid syntax!"<<endl;
    cout<<"----------------------------------------------------------------------"<<endl;
    cout<<"Correct Syntax for NTuple Reduce"<<endl;
    cout<<"\t"<<argv[0]<<" "<<runOption_Reduce<<"\n"<<endl;
    
    cout<<"Correct Syntax for NTuple Analysis"<<endl;
    cout<<"\t"<<argv[0]<<" runOption modeOption\n"<<endl;
    cout<<"runOptions:"<<endl;
    cout.width(16);cout<<"\tRun:"<<runOption_Run<<endl;
    cout.width(16);cout<<"\tPlot:"<<runOption_Plot<<"\n"<<endl;
    
    cout<<"modeOptions:"<<endl;
    cout.width(16); cout<<"\tSignal:"<<modeOption_Signal<<endl;
    cout.width(16); cout<<"\tBackground:"<<modeOption_Background<<endl;
    cout.width(16); cout<<"\tEverything:"<<modeOption_All<<endl;
    cout.width(16); cout<<"\tPlot Only: Signal vs Background:"<<modeOption_Other<<endl;
    cout<<"----------------------------------------------------------------------"<<endl;
}

bool isModePlot(string runSelect)
{
    if(runSelect.compare(runOption_Run) == 0){
        cout<<"Run Option Selected!"<<endl;
        return false;
    }
    else{
        cout<<"Plot Option Selected!"<<endl;
        return true;
    }
}

int GetMode(string modeSelect)
{
    int nMode;

    if(modeSelect.compare(modeOption_Signal) == 0){
        cout<<"Signal Mode Selected!"<<endl;
        nMode = 1;
    }
    else if(modeSelect.compare(modeOption_Background) == 0){
        cout<<"Background Mode Selected!"<<endl;
        nMode = 2;
    }
    else if(modeSelect.compare(modeOption_All) == 0){
        cout<<"All Events Mode Selected!"<<endl;
        nMode = 3;
    }else if (modeSelect.compare(modeOption_Other) == 0){
        cout<<"Plot Only! Signal vs Background Mode"<<endl;
        nMode = 4;  
    }

    return nMode;
}



