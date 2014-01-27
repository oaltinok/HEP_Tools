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
    
    Last Revision: 2014_01_27
================================================================================
*/

// Include Required Classes & Libraries
#include "Classes/ANA_CCNuPionInc/ANA_CCNuPionInc.cpp"
#include "Classes/Plotter/Plotter.cpp"
#include "Libraries/Folder_List.h"

#include <string>

using namespace std;

string getFileLocation(string s1, string s2);
string getFileLocation(string s1, string s2, string s3);

void main(bool onlyPlot)
{
    ANA_CCNuPionInc t;
    Plotter p;
    
    // Edit isTest Variable for running Test Samples or complete playlist
    bool isTest;
    isTest = true;
    
    if(onlyPlot){
            p.plotHistograms(getFileLocation(Folder_List::OUT, Folder_List::ROOT,Folder_List::F_ROOT_MC),
            getFileLocation(Folder_List::OUT, Folder_List::PLOT));
        
    }else{
        if(isTest){
            t.run(getFileLocation(Folder_List::IN, Folder_List::PLAYLIST, Folder_List::F_PL_TEST),
            getFileLocation(Folder_List::OUT, Folder_List::ROOT, Folder_List::F_ROOT_MC),
            getFileLocation(Folder_List::OUT, Folder_List::TEXT, Folder_List::F_TEXT_CUT),
            getFileLocation(Folder_List::OUT, Folder_List::TEXT, Folder_List::F_TEXT_README) );
            
        }else{
            t.run(getFileLocation(Folder_List::IN, Folder_List::PLAYLIST, Folder_List::F_PL_FULL),
            getFileLocation(Folder_List::OUT, Folder_List::ROOT, Folder_List::F_ROOT_MC),
            getFileLocation(Folder_List::OUT, Folder_List::TEXT, Folder_List::F_TEXT_CUT),
            getFileLocation(Folder_List::OUT, Folder_List::TEXT, Folder_List::F_TEXT_README) );
        }
    }
    
}

string getFileLocation(string s1, string s2)
{
    return s1+s2;
}

string getFileLocation(string s1, string s2, string s3)
{
    return s1+s2+s3;
}

