/*
================================================================================
main.cpp
    Main Function that controls the Analysis Package
    
    Build Package using "make"
    
Author:        Ozgur Altinok  - ozgur.altinok@tufts.edu
================================================================================
*/
#include <string>
#include <ctime>
#include <iostream>
#include <string>
#include <TROOT.h>

#include "Cintex/Cintex.h"
#include "Classes/SideBandTool/CCProtonPi0_SideBandTool.h"
#include <PlotUtils/MnvPlotter.h>

using namespace std;

// See src/Main_functions.cpp for definitions
int GetMode(int argc, char* argv[]);
void showInputError(char *argv[]);
void Plot();
void Reduce(string playlist, bool isMC);
void Analyze(string playlist, bool isMC);
void Calculate_CrossSection(bool isMC);

// See src/Minuit_Functions.cpp for definitions 
void FitSideBands();
void FitMinuit();
void FitMinuit_W();
void ApplyStyle(PlotUtils::MnvPlotter* plotter);
void Fit(std::string fileName, std::string plotName);
void Fit(std::string rootDir_data, std::string rootDir_mc, std::string plotName);
void Fit_QSq();
void calc_ChiSq(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void calc_ChiSq_W(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
double calc_ChiSq_SideBand(SideBand &sb, Double_t *par, int unv = 0, bool isPartial = false, int min_bin = 1, int max_bin = 1);

int main(int argc, char *argv[] )
{
    time_t timeStart; time(&timeStart);
    time_t timeEnd;
    double timeDiff;
    int timeDiff_m;
    int timeDiff_s;
    string pl_reduce;
    string pl_analyze;
    bool isMC;

    TH1::AddDirectory(false);
    
    // Check User Command
    int nMode = GetMode(argc,argv);
    if (nMode == 0){
        showInputError(argv);
        return 0;
    }

    if ( nMode < 0) isMC = true;
    else isMC = false;

    if (nMode != 10 || nMode != 20 || nMode != 30 || nMode != 40){
        if (isMC){
            cout<<"MC Playlists Selected!\n"<<endl;
            pl_reduce = "Input/Playlists/pl_MC_Merged.dat"; 
            pl_analyze = "Input/Playlists/pl_MC_Reduced.dat"; 
        }else{
            cout<<"Data Playlists Selected!\n"<<endl;
            pl_reduce = "Input/Playlists/pl_Data_Merged.dat"; 
            pl_analyze = "Input/Playlists/pl_Data_Reduced.dat"; 
        }
    }

    ROOT::Cintex::Cintex::Enable();

    if ( abs(nMode) == 1) Reduce(pl_reduce, isMC);
    else if ( abs(nMode) == 2) Analyze(pl_analyze, isMC);
    else if ( abs(nMode) == 3) Calculate_CrossSection(isMC);
    else if ( nMode == 10) Plot();
    else if ( nMode == 20) FitSideBands();
    else if ( nMode == 30) FitMinuit_W();
    else if ( nMode == 40) Fit_QSq();
    else{
        cout<<"Problem on Mode!, Returning"<<endl;
        return 0;
    }

    time(&timeEnd);
    timeDiff = ( timeEnd - timeStart );
    
    timeDiff_m = (int)timeDiff / 60;
    timeDiff_s = (int)timeDiff % 60;
    
    cout<<"Time to Complete = "<<timeDiff_m<<"\' "<<timeDiff_s<<"\" or "<<timeDiff<<" seconds"<<endl;
    return 0;
}

