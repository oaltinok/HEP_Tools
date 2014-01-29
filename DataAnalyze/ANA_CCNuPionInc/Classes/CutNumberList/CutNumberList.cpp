/*
    See CutNumberList.h header for Class Information
*/
#include "CutNumberList.h"

using namespace std;

CutNumberList::CutNumberList()
{
    // Set Initial Values to Zero
    nAll.setValue(initCount);
    nVolume.setValue(initCount);
    nMuon.setValue(initCount);
    nProton.setValue(initCount);
    nPion.setValue(initCount);
    
    // Set Labels
    nAll.setLabel("nAll:\t\t");
    nVolume.setLabel("nVolume:\t");
    nMuon.setLabel("nMuon:\t");
    nProton.setLabel("nProton:\t");
    nPion.setLabel("nPion:\t");
    

}


double CutNumberList::getPercent(double nCurrent)
{
    double percent;

    percent = ( nCurrent / nAll.getValue() ) * 100;
    
    return percent;

}

