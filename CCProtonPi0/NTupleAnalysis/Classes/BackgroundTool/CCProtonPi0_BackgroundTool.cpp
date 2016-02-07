#ifndef CCProtonPi0_BackgroundTool_cpp
#define CCProtonPi0_BackgroundTool_cpp

#include "CCProtonPi0_BackgroundTool.h"

using namespace std;

CCProtonPi0_BackgroundTool::CCProtonPi0_BackgroundTool(bool isModeReduce) : CCProtonPi0_NTupleAnalysis()
{
    cout<<"Initializing CCProtonPi0_BackgroundTool"<<endl;
    if(isModeReduce){
        cout<<"\tNTuple Reduce Mode -- Will not create Text Files"<<endl;
    }else{
        OpenTextFile();
    
        initBackgrounds();
    }
    cout<<"Done!"<<endl;
}

void CCProtonPi0_BackgroundTool::initBackgrounds()
{   
    // Background With Pi0
    initSingleBackground(bckg_NoPi0, "NoPi0");
    initSingleBackground(bckg_SinglePi0, "SinglePi0");
    initSingleBackground(bckg_MultiPi0, "MultiPi0");
    initSingleBackground(bckg_Total_WithPi0, "TotalWithPi0");

    // Background Type
    initSingleBackground(bckg_NC, "NC");
    initSingleBackground(bckg_AntiNeutrino, "AntiNeutrino");
    initSingleBackground(bckg_QELike, "QELike");
    initSingleBackground(bckg_SingleChargedPion, "SingleChargedPion");
    initSingleBackground(bckg_SingleChargedPion_ChargeExchanged, "SingleChargedPion_ChargeExc");
    initSingleBackground(bckg_DoublePionWithPi0, "DoublePionWithPi0");
    initSingleBackground(bckg_DoublePionWithoutPi0, "DoublePionWithoutPi0");
    initSingleBackground(bckg_MultiPionWithPi0, "MultiPionWithPi0");
    initSingleBackground(bckg_MultiPionWithoutPi0, "MultiPionWithoutPi0");
    initSingleBackground(bckg_Other, "Other");
    initSingleBackground(bckg_Total, "Total");
}

void CCProtonPi0_BackgroundTool::initSingleBackground(Background &b, string input_name)
{
    b.name = input_name;
    b.nAll = 0.0;
    b.nWithMichel = 0.0;
}

void CCProtonPi0_BackgroundTool::updateBackground(Background &b, bool withMichel)
{
    // Increment Events with that Background
    b.nAll++;

    // Increment Events with Michel
    if(withMichel) b.nWithMichel++;
}

void CCProtonPi0_BackgroundTool::fillBackgroundWithPi0(bool NoPi0, bool SinglePi0, bool MultiPi0, bool withMichel)
{
    // Update Total Background
    updateBackground(bckg_Total_WithPi0,withMichel);

    // Update Each Background
    if (NoPi0) updateBackground(bckg_NoPi0,withMichel);
    else if (SinglePi0) updateBackground(bckg_SinglePi0,withMichel);
    else if (MultiPi0) updateBackground(bckg_MultiPi0,withMichel);
    else cout<<"WARNING! No BackgroundWithPi0 Found"<<endl;
}

void CCProtonPi0_BackgroundTool::fillBackground(bool NC, bool AntiNeutrino, bool QELike, bool SingleChargedPion, bool SingleChargedPion_ChargeExchanged, bool DoublePionWithPi0, bool DoublePionWithoutPi0,  bool MultiPionWithPi0, bool MultiPionWithoutPi0, bool Other, bool withMichel)
{
    // Update Total Background
    updateBackground(bckg_Total,withMichel);
    
    // Update Each Background
    if (NC) updateBackground(bckg_NC,withMichel);
    else if (AntiNeutrino) updateBackground(bckg_AntiNeutrino,withMichel);
    else if (QELike) updateBackground(bckg_QELike,withMichel);
    else if (SingleChargedPion) updateBackground(bckg_SingleChargedPion,withMichel);
    else if (SingleChargedPion_ChargeExchanged) updateBackground(bckg_SingleChargedPion_ChargeExchanged,withMichel);
    else if (DoublePionWithPi0) updateBackground(bckg_DoublePionWithPi0,withMichel);
    else if (DoublePionWithoutPi0) updateBackground(bckg_DoublePionWithoutPi0,withMichel);
    else if (MultiPionWithPi0) updateBackground(bckg_MultiPionWithPi0,withMichel);
    else if (MultiPionWithoutPi0) updateBackground(bckg_MultiPionWithoutPi0,withMichel);
    else if (Other) updateBackground(bckg_Other,withMichel);
    else cout<<"WARNING! No Background Type Found!"<<endl;
}

void CCProtonPi0_BackgroundTool::writeBackgroundTableHeader()
{
    textFile<<std::left;
    
    // Table Header
    textFile.width(40); textFile<< "Background Type";
    textFile.width(16); textFile<< "N(Events)"; 
    textFile.width(16); textFile<< "Percent"; 
    //textFile.width(16); textFile<< "N(WithMichel)"; 
    //textFile.width(16); textFile<< "Percent"; 
    textFile<<endl;
}


void CCProtonPi0_BackgroundTool::formBackgroundVectors()
{
    // Background With Pi0
    BackgroundWithPi0Vector.push_back(bckg_NoPi0);
    BackgroundWithPi0Vector.push_back(bckg_SinglePi0);
    BackgroundWithPi0Vector.push_back(bckg_MultiPi0);
    BackgroundWithPi0Vector.push_back(bckg_Total_WithPi0);
    
    // Background Type
    BackgroundTypeVector.push_back(bckg_NC);
    BackgroundTypeVector.push_back(bckg_AntiNeutrino);
    BackgroundTypeVector.push_back(bckg_QELike);
    BackgroundTypeVector.push_back(bckg_SingleChargedPion);
    BackgroundTypeVector.push_back(bckg_SingleChargedPion_ChargeExchanged);
    BackgroundTypeVector.push_back(bckg_DoublePionWithPi0);
    BackgroundTypeVector.push_back(bckg_DoublePionWithoutPi0);
    BackgroundTypeVector.push_back(bckg_MultiPionWithPi0);
    BackgroundTypeVector.push_back(bckg_MultiPionWithoutPi0);
    BackgroundTypeVector.push_back(bckg_Other);
    BackgroundTypeVector.push_back(bckg_Total);
}

double CCProtonPi0_BackgroundTool::calcPercent(double nEvents, double nBase)
{
    return (nEvents / nBase * 100.0);
}

void CCProtonPi0_BackgroundTool::writeBackgroundTableRows(vector< Background > &bckgVector)
{
    // Write Background With Pi0
    for (unsigned int i = 0; i < bckgVector.size(); i++){
        Background temp = bckgVector[i];

        textFile.width(40); textFile<<temp.name;  
        textFile.width(16); textFile<<temp.nAll;    
        textFile.width(16); textFile<<calcPercent(temp.nAll,bckg_Total.nAll);    
        //textFile.width(16); textFile<<temp.nWithMichel;   
        //textFile.width(16); textFile<<calcPercent(temp.nWithMichel,bckg_Total.nWithMichel);    
        textFile<<endl;
    }

    textFile<<endl;
}

void CCProtonPi0_BackgroundTool::writeBackgroundTable()
{
    formBackgroundVectors();

    cout<<">> Writing "<<fileName<<endl;

    writeBackgroundTableHeader();
    writeBackgroundTableRows(BackgroundWithPi0Vector);
    writeBackgroundTableRows(BackgroundTypeVector);
}

void CCProtonPi0_BackgroundTool::OpenTextFile()
{
    // Open Background Files
    fileName = Folder_List::output + Folder_List::textOut + "BackgroundTable.txt";

    textFile.open( fileName.c_str() );
    if( !textFile.is_open() ){
        cerr<<"Cannot open output text file: "<<fileName<<endl;
        exit(1);
    }else{
        cout<<"\t"<<fileName<<endl;
    }
}

CCProtonPi0_BackgroundTool::~CCProtonPi0_BackgroundTool()
{
    textFile.close(); 
}



#endif
