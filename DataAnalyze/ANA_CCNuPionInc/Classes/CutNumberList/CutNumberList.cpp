/*
    See CutNumberList.h header for Class Information
*/
#include "CutNumberList.h"

using namespace std;

CutNumberList::CutNumberList()
{

    // -------------------------------------------------------------------------
    //     File for Cut Table
    //--------------------------------------------------------------------------
    // Open Cut File
    cutFile = getFileLocation(Folder_List::OUT, Folder_List::TEXT, Folder_List::F_TEXT_CUT);
    cutText.open( cutFile.c_str() );
    if( !cutText.is_open() ){
        cerr<<"Cannot open output text file: "<<cutFile<<endl;
        exit(1);
    }
    
    // -------------------------------------------------------------------------
    //     Memory Allocation
    //--------------------------------------------------------------------------
    
    nAll = new SingleCutNumber;
    nVolume = new SingleCutNumber;
    nMuon = new SingleCutNumber;
    nProton = new SingleCutNumber;
    nPion = new SingleCutNumber;
    
    // -------------------------------------------------------------------------
    //     Initialization
    //--------------------------------------------------------------------------
    // Set Initial Values to Zero
    nAll->setValue(initCount);
    nVolume->setValue(initCount);
    nMuon->setValue(initCount);
    nProton->setValue(initCount);
    nPion->setValue(initCount);
    
    // Set Labels
    nAll->setLabel("nAll:\t\t");
    nVolume->setLabel("nVolume:\t");
    nMuon->setLabel("nMuon:\t\t");
    nProton->setLabel("nProton:\t");
    nPion->setLabel("nPion:\t\t");
    

}


double CutNumberList::getPercent(double nCurrent)
{
    double percent;

    percent = ( nCurrent / nAll->getValue() ) * 100;
    
    return percent;

}

void CutNumberList::writeCutTable()
{
    cout<<">> Writing "<<cutFile<<endl;
    
    writeCut(nAll);
    writeCut(nVolume);
    writeCut(nMuon);
    writeCut(nProton);
    writeCut(nPion);
    
    cutText.close();
}

void CutNumberList::writeCut(SingleCutNumber* cut)
{
    cutText<<cut->getLabel()<<cut->getValue();
    cutText<<"\tPercent : "<<"\t"<<getPercent(cut->getValue())<<endl;
}

