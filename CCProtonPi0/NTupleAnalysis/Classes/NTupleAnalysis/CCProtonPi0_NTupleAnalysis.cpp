#ifndef CCProtonPi0_NTupleAnalysis_cpp
#define CCProtonPi0_NTupleanalysis_cpp

#include "CCProtonPi0_NTupleAnalysis.h"

using namespace std;

CCProtonPi0_NTupleAnalysis::CCProtonPi0_NTupleAnalysis(int nMode)
{
    ROOT::Cintex::Cintex::Enable();
    
    SetAnalysisMode(nMode);
}

void CCProtonPi0_NTupleAnalysis::SetAnalysisMode(int nMode)
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
