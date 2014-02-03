/*
    See CutNumberList.h header for Class Information
*/
#include "CutNumberList.h"

using namespace std;

CutNumberList::CutNumberList()
{
    // Open Cut File
    cutFile = getFileLocation(Folder_List::OUT, Folder_List::TEXT, Folder_List::F_TEXT_CUT);
    cutText.open( cutFile.c_str() );
    if( !cutText.is_open() ){
        cerr<<"Cannot open output text file: "<<cutFile<<endl;
        exit(1);
    }
    
    // Set Initial Values to Zero
    nAll.setValue(initCount);
    nVolume.setValue(initCount);
    nMuon.setValue(initCount);
    nProton.setValue(initCount);
    nPion.setValue(initCount);
    
    // Set Labels
    nAll.setLabel("nAll:\t\t");
    nVolume.setLabel("nVolume:\t");
    nMuon.setLabel("nMuon:\t\t");
    nProton.setLabel("nProton:\t");
    nPion.setLabel("nPion:\t\t");
    

}


double CutNumberList::getPercent(double nCurrent)
{
    double percent;

    percent = ( nCurrent / nAll.getValue() ) * 100;
    
    return percent;

}

void CutNumberList::writeCutTable()
{
    cout<<">> Writing "<<cutFile<<endl;

    cutText<<nAll.getLabel()<<nAll.getValue();
    cutText<<"\tPercent : "<<"\t"<<getPercent(nAll.getValue())<<endl;
    
    cutText<<nVolume.getLabel()<<nVolume.getValue();
    cutText<<"\tPercent : "<<"\t"<<getPercent(nVolume.getValue())<<endl;
    
    cutText<<nMuon.getLabel()<<nMuon.getValue();
    cutText<<"\tPercent : "<<"\t"<<getPercent(nMuon.getValue())<<endl;
    
    cutText<<nProton.getLabel()<<nProton.getValue();
    cutText<<"\tPercent : "<<"\t"<<getPercent(nProton.getValue())<<endl;
    
    cutText<<nPion.getLabel()<<nPion.getValue();
    cutText<<"\tPercent : "<<"\t"<<getPercent(nPion.getValue())<<endl;
    
    cutText.close();
}

