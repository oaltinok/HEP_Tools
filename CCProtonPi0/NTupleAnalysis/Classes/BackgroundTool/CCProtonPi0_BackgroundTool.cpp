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
        OpenTextFiles();

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

    // Background Type Compact
    initSingleBackground(bckg_compact_WithPi0,"With Pi0");
    initSingleBackground(bckg_compact_QELike,"QE Like");
    initSingleBackground(bckg_compact_SinglePiPlus,"Single Pi Plus");
    initSingleBackground(bckg_compact_Other,"Other");
    initSingleBackground(bckg_compact_Total,"Total");

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

void CCProtonPi0_BackgroundTool::initSingleBackground(vector<Background> &b, string input_name)
{
    Background temp;
    temp.name = input_name;
    temp.nAll = 0.0;
    temp.nWithMichel = 0.0;

    for(int i = 0; i < nTables; ++i){
        b.push_back(temp);
    }
}

void CCProtonPi0_BackgroundTool::updateBackground(Background &b, bool withMichel, double wgt)
{
    // Increment Events with that Background
    b.nAll = b.nAll + wgt;

    // Increment Events with Michel
    if(withMichel) b.nWithMichel = b.nWithMichel + wgt;
}


void CCProtonPi0_BackgroundTool::fillBackgroundCompact(bool WithPi0, bool QELike, bool SinglePiPlus, bool Other, double wgt)
{
    // ----------------------------------------------------------------------------------
    // Fill Table 0 for Inclusive
    // ----------------------------------------------------------------------------------
    int ind = 0;
    // Update Total Background
    updateBackground(bckg_compact_Total[ind], false, wgt);

    // Update Each Background
    if (WithPi0) updateBackground(bckg_compact_WithPi0[ind], false, wgt);
    else if (QELike) updateBackground(bckg_compact_QELike[ind], false, wgt);
    else if (SinglePiPlus) updateBackground(bckg_compact_SinglePiPlus[ind], false, wgt);
    else if (Other) updateBackground(bckg_compact_Other[ind], false, wgt);
    else cout<<"WARNING! No BackgroundCompact Found"<<endl;

    // ----------------------------------------------------------------------------------
    // Fill Table nTracks for Topology Dependency
    // ----------------------------------------------------------------------------------
    ind = nTracks;
    // Update Total Background
    updateBackground(bckg_compact_Total[ind], false, wgt);

    // Update Each Background
    if (WithPi0) updateBackground(bckg_compact_WithPi0[ind], false, wgt);
    else if (QELike) updateBackground(bckg_compact_QELike[ind], false, wgt);
    else if (SinglePiPlus) updateBackground(bckg_compact_SinglePiPlus[ind], false, wgt);
    else if (Other) updateBackground(bckg_compact_Other[ind], false, wgt);
    else cout<<"WARNING! No BackgroundCompact Found"<<endl;
}

void CCProtonPi0_BackgroundTool::fillBackgroundWithPi0(bool NoPi0, bool SinglePi0, bool MultiPi0, bool withMichel, double wgt)
{
    // ----------------------------------------------------------------------------------
    // Fill Table 0 for Inclusive
    // ----------------------------------------------------------------------------------
    int ind = 0;
    // Update Total Background
    updateBackground(bckg_Total_WithPi0[ind],withMichel, wgt);

    // Update Each Background
    if (NoPi0) updateBackground(bckg_NoPi0[ind],withMichel, wgt);
    else if (SinglePi0) updateBackground(bckg_SinglePi0[ind],withMichel, wgt);
    else if (MultiPi0) updateBackground(bckg_MultiPi0[ind],withMichel, wgt);
    else cout<<"WARNING! No BackgroundWithPi0 Found"<<endl;

    // ----------------------------------------------------------------------------------
    // Fill Table nTracks for Topology Dependency
    // ----------------------------------------------------------------------------------
    ind = nTracks;
    // Update Total Background
    updateBackground(bckg_Total_WithPi0[ind],withMichel, wgt);

    // Update Each Background
    if (NoPi0) updateBackground(bckg_NoPi0[ind],withMichel, wgt);
    else if (SinglePi0) updateBackground(bckg_SinglePi0[ind],withMichel, wgt);
    else if (MultiPi0) updateBackground(bckg_MultiPi0[ind],withMichel, wgt);
    else cout<<"WARNING! No BackgroundWithPi0 Found"<<endl;
}

void CCProtonPi0_BackgroundTool::fillBackground(bool NC, bool AntiNeutrino, bool QELike, bool SingleChargedPion, bool SingleChargedPion_ChargeExchanged, bool DoublePionWithPi0, bool DoublePionWithoutPi0,  bool MultiPionWithPi0, bool MultiPionWithoutPi0, bool Other, bool withMichel, double wgt)
{
    // ----------------------------------------------------------------------------------
    // Fill Table 0 for Inclusive
    // ----------------------------------------------------------------------------------
    int ind = 0;
    // Update Total Background
    updateBackground(bckg_Total[ind],withMichel, wgt);

    // Update Each Background
    if (NC) updateBackground(bckg_NC[ind],withMichel, wgt);
    else if (AntiNeutrino) updateBackground(bckg_AntiNeutrino[ind],withMichel, wgt);
    else if (QELike) updateBackground(bckg_QELike[ind],withMichel, wgt);
    else if (SingleChargedPion) updateBackground(bckg_SingleChargedPion[ind],withMichel, wgt);
    else if (SingleChargedPion_ChargeExchanged) updateBackground(bckg_SingleChargedPion_ChargeExchanged[ind],withMichel, wgt);
    else if (DoublePionWithPi0) updateBackground(bckg_DoublePionWithPi0[ind],withMichel, wgt);
    else if (DoublePionWithoutPi0) updateBackground(bckg_DoublePionWithoutPi0[ind],withMichel, wgt);
    else if (MultiPionWithPi0) updateBackground(bckg_MultiPionWithPi0[ind],withMichel, wgt);
    else if (MultiPionWithoutPi0) updateBackground(bckg_MultiPionWithoutPi0[ind],withMichel, wgt);
    else if (Other) updateBackground(bckg_Other[ind],withMichel, wgt);
    else cout<<"WARNING! No Background Type Found!"<<endl;

    // ----------------------------------------------------------------------------------
    // Fill Table nTracks for Topology Dependency
    // ----------------------------------------------------------------------------------
    ind = nTracks;
    // Update Total Background
    updateBackground(bckg_Total[ind],withMichel, wgt);

    // Update Each Background
    if (NC) updateBackground(bckg_NC[ind],withMichel, wgt);
    else if (AntiNeutrino) updateBackground(bckg_AntiNeutrino[ind],withMichel, wgt);
    else if (QELike) updateBackground(bckg_QELike[ind],withMichel, wgt);
    else if (SingleChargedPion) updateBackground(bckg_SingleChargedPion[ind],withMichel, wgt);
    else if (SingleChargedPion_ChargeExchanged) updateBackground(bckg_SingleChargedPion_ChargeExchanged[ind],withMichel, wgt);
    else if (DoublePionWithPi0) updateBackground(bckg_DoublePionWithPi0[ind],withMichel, wgt);
    else if (DoublePionWithoutPi0) updateBackground(bckg_DoublePionWithoutPi0[ind],withMichel, wgt);
    else if (MultiPionWithPi0) updateBackground(bckg_MultiPionWithPi0[ind],withMichel, wgt);
    else if (MultiPionWithoutPi0) updateBackground(bckg_MultiPionWithoutPi0[ind],withMichel, wgt);
    else if (Other) updateBackground(bckg_Other[ind],withMichel, wgt);
    else cout<<"WARNING! No Background Type Found!"<<endl;
}

void CCProtonPi0_BackgroundTool::writeBackgroundTableHeader()
{
    for (int i = 0; i < nTables; ++i){
        textFile[i]<<std::left;

        // Table Header
        textFile[i].width(40); textFile[i]<<"Background Type";
        textFile[i].width(16); textFile[i]<<"N(Events)"; 
        textFile[i].width(16); textFile[i]<<"Percent"; 
        textFile[i]<<endl;
    }
}


void CCProtonPi0_BackgroundTool::formBackgroundVectors()
{
    // Background With Pi0
    BackgroundWithPi0Vector.push_back(bckg_NoPi0);
    BackgroundWithPi0Vector.push_back(bckg_SinglePi0);
    BackgroundWithPi0Vector.push_back(bckg_MultiPi0);
    BackgroundWithPi0Vector.push_back(bckg_Total_WithPi0);

    // Background Type Compact
    BackgroundCompactVector.push_back(bckg_compact_WithPi0);
    BackgroundCompactVector.push_back(bckg_compact_QELike);
    BackgroundCompactVector.push_back(bckg_compact_SinglePiPlus);
    BackgroundCompactVector.push_back(bckg_compact_Other);
    BackgroundCompactVector.push_back(bckg_compact_Total);

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

void CCProtonPi0_BackgroundTool::writeBackgroundTableRows(vector< vector<Background> > &bckgVector)
{
    for (unsigned int i = 0; i < bckgVector.size(); ++i){
        vector<Background> temp = bckgVector[i];

        for (int j = 0; j < nTables; ++j){
            textFile[j].width(40); textFile[j]<<temp[j].name;  
            textFile[j].width(16); textFile[j]<<temp[j].nAll;    
            textFile[j].width(16); textFile[j]<<calcPercent(temp[j].nAll, bckg_Total[j].nAll);    
            textFile[j]<<endl;
        }
    }

    // Add an empty line after each background type
    for (int j = 0; j < nTables; ++j){
        textFile[j]<<endl;
    }
}

void CCProtonPi0_BackgroundTool::writeBackgroundTable()
{
    formBackgroundVectors();

    cout<<">> Writing Background Tables"<<endl;

    writeBackgroundTableHeader();
    writeBackgroundTableRows(BackgroundCompactVector);
    writeBackgroundTableRows(BackgroundWithPi0Vector);
    writeBackgroundTableRows(BackgroundTypeVector);
}

void CCProtonPi0_BackgroundTool::OpenTextFiles()
{
    // Open Background Files
    fileName[0] = Folder_List::output + Folder_List::textOut + "BackgroundTable_All.txt";
    fileName[1] = Folder_List::output + Folder_List::textOut + "BackgroundTable_1Track.txt";
    fileName[2] = Folder_List::output + Folder_List::textOut + "BackgroundTable_2Track.txt";

    for (int i = 0; i < nTables; ++i){
        OpenTextFile(fileName[i], textFile[i]);
    }
}

CCProtonPi0_BackgroundTool::~CCProtonPi0_BackgroundTool()
{
    for (int i = 0; i < nTables; ++i){
        textFile[i].close(); 
    }
}

void CCProtonPi0_BackgroundTool::set_nTracks(int input)
{
    nTracks = input;
}

#endif
