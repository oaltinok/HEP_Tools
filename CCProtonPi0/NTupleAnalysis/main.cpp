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
    Analyze Reduced File
            > ./main.exe run
    Plot Analyzed File
            > ./main.exe plot

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

bool isMC = true;

const string runOption_Run = "run";
const string runOption_Plot = "plot";
const string runOption_Reduce = "reduce";

int GetMode(int argc, char* argv[]);
void showInputError(char *argv[]);
void Plot();
void Reduce(string playlist);
void Analyze(string playlist);

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
    int nMode = GetMode(argc,argv);
    if (nMode == -1){
        showInputError(argv);
        return 0;
    }
    
    if (isMC && nMode != 2){
        cout<<"MC Playlists Selected!"<<endl;
        pl_reduce = "Input/Playlists/pl_MC_Merged.dat"; 
        pl_analyze = "Input/Playlists/pl_MC_Reduced.dat"; 
    }else{
        cout<<"Data Playlists Selected!"<<endl;
        pl_reduce = "Input/Playlists/pl_Data_Merged.dat"; 
        pl_analyze = "Input/Playlists/pl_Data_Reduced.dat"; 
    }
   
    ROOT::Cintex::Cintex::Enable();


    if (nMode == 0) Reduce(pl_reduce);
    else if (nMode == 1) Analyze(pl_analyze);
    else Plot();
    
    time(&timeEnd);
    timeDiff = ( timeEnd - timeStart );
    
    timeDiff_m = (int)timeDiff / 60;
    timeDiff_s = (int)timeDiff % 60;
    
    cout<<"Time to Complete = "<<timeDiff_m<<"\' "<<timeDiff_s<<"\" or "<<timeDiff<<" seconds"<<endl;
    return 0;
}

void Reduce(string playlist)
{
    bool isModeReduce = true;
    CCProtonPi0_Analyzer t(isModeReduce, isMC);
    t.reduce(playlist);
}

void Analyze(string playlist)
{
    bool isModeReduce = false;
    // First Analyze 1Track Events
    cout<<"======================================================================"<<endl;
    cout<<"Analyzing 1 Track Events..."<<endl;
    cout<<"======================================================================"<<endl;
    CCProtonPi0_Analyzer analyzer_1Track(isModeReduce, isMC, "1Track/");
    analyzer_1Track.analyze(playlist);

    // Second Analyze 2Track Events 
    cout<<"\n"<<endl;
    cout<<"======================================================================"<<endl;
    cout<<"Analyzing 2+ Track Events..."<<endl;
    cout<<"======================================================================"<<endl;
    CCProtonPi0_Analyzer analyzer_2Track(isModeReduce, isMC, "2Track/");
    analyzer_2Track.analyze(playlist);

    // Finally Analyze All Events 
    cout<<"\n"<<endl;
    cout<<"======================================================================"<<endl;
    cout<<"Analyzing All Events..."<<endl;
    cout<<"======================================================================"<<endl;
    CCProtonPi0_Analyzer analyzer_All(isModeReduce, isMC, "All/");
    analyzer_All.analyze(playlist);
}

void Plot()
{
    // First Plot 1Track Events
    cout<<"======================================================================"<<endl;
    cout<<"Plotting 1 Track Events..."<<endl;
    cout<<"======================================================================"<<endl;
    CCProtonPi0_Plotter plotter_1Track("1Track/");
    plotter_1Track.plotHistograms();

    // Second Plot 2+Track Events
    cout<<"======================================================================"<<endl;
    cout<<"Plotting 2+ Track Events..."<<endl;
    cout<<"======================================================================"<<endl;
    CCProtonPi0_Plotter plotter_2Track("2Track/");
    plotter_2Track.plotHistograms();

    // Finally Plot All Events
    cout<<"======================================================================"<<endl;
    cout<<"Plotting All Events..."<<endl;
    cout<<"======================================================================"<<endl;
    CCProtonPi0_Plotter plotter_all("All/");
    plotter_all.plotHistograms();
}

int GetMode(int argc, char* argv[])
{
    // argc can only be 2
    if (argc != 2) return -1;

    std::string runSelect = argv[1];
    
    if (runSelect.compare(runOption_Reduce) == 0) return 0;
    else if (runSelect.compare(runOption_Run) == 0) return 1;
    else if (runSelect.compare(runOption_Plot) == 0) return 2;
    else return -1;
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
    cout<<"\t"<<argv[0]<<" "<<runOption_Run<<"\n"<<endl;
    cout<<"Correct Syntax for Plotting"<<endl;
    cout<<"\t"<<argv[0]<<" "<<runOption_Plot<<"\n"<<endl;
    cout<<"----------------------------------------------------------------------"<<endl;
}



