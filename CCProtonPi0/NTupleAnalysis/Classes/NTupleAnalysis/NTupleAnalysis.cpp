#ifndef NTupleAnalysis_cpp
#define NTupleanalysis_cpp

#include "NTupleAnalysis.h"

using namespace std;

NTupleAnalysis::NTupleAnalysis(int nMode)
{
    SetAnalysisMode(nMode);
}

void NTupleAnalysis::SetAnalysisMode(int nMode)
{
    // Save Analysis Mode
    anaMode = nMode;
    isAnalysisModeSelected = true;
    
    // Set branchDir - Output Folder
    if ( nMode == 1) branchDir = Folder_List::signal;
    else if ( nMode == 2) branchDir = Folder_List::background;
    else branchDir = Folder_List::allEvents;
    
}


#endif
