#ifndef CCProtonPi0_BackgroundTool_cpp
#define CCProtonPi0_BackgroundTool_cpp

#include "CCProtonPi0_BackgroundTool.h"

using namespace std;

CCProtonPi0_BackgroundTool::CCProtonPi0_BackgroundTool(int nMode) : CCProtonPi0_NTupleAnalysis(nMode)
{
    cout<<"Initializing CCProtonPi0_BackgroundTool"<<endl;
    if(nMode == 0){
        cout<<"\tNTuple Reduce Mode -- Will not create Text Files"<<endl;
    }else{
        OpenTextFiles();
    
        initBackgroundVectors();
    }
    cout<<"Done!"<<endl;
    
}

void CCProtonPi0_BackgroundTool::initBackgroundVectors()
{   
    // Background With Pi0
    initSingleBackgroundVector(bckg_NoPi0, "NoPi0");
    initSingleBackgroundVector(bckg_SinglePi0, "SinglePi0");
    initSingleBackgroundVector(bckg_MultiPi0, "MultiPi0");
    initSingleBackgroundVector(bckg_Total_WithPi0, "TotalWithPi0");

    // Background Type
    initSingleBackgroundVector(bckg_NC, "NC");
    initSingleBackgroundVector(bckg_AntiNeutrino, "AntiNeutrino");
    initSingleBackgroundVector(bckg_QELike, "QELike");
    initSingleBackgroundVector(bckg_SinglePion, "SinglePion");
    initSingleBackgroundVector(bckg_DoublePion, "DoublePion");
    initSingleBackgroundVector(bckg_MultiPion, "MultiPion");
    initSingleBackgroundVector(bckg_Other, "Other");
    initSingleBackgroundVector(bckg_Total, "Total");
}

void CCProtonPi0_BackgroundTool::initSingleBackgroundVector(vector<Background> &b, string input_name)
{
    // Declare and Initialize a Temporary Background Object
    Background temp;
    temp.name = input_name;
    temp.nAll = 0.0;
    temp.nWithMichel = 0.0;

    // Push Temporary Background Object to BackgroundVector
    for (int i = 0; i < nTopologies; i++){
        b.push_back(temp);
    }

}

int CCProtonPi0_BackgroundTool::setArrayInd(int nProngs)
{
    if (nProngs == 1) return 0;
    else return 1;
}

void CCProtonPi0_BackgroundTool::updateBackground(Background &b, bool withMichel)
{
    // Increment Events with that Background
    b.nAll++;

    // Increment Events with Michel
    if(withMichel) b.nWithMichel++;
}

void CCProtonPi0_BackgroundTool::fillBackgroundWithPi0( int nProngs, bool NoPi0, bool SinglePi0, bool MultiPi0, bool withMichel)
{
    int ind = setArrayInd(nProngs);

    // Update Total Background
    updateBackground(bckg_Total_WithPi0[ind],withMichel);

    // Update Each Background
    if (NoPi0) updateBackground(bckg_NoPi0[ind],withMichel);
    else if (SinglePi0) updateBackground(bckg_SinglePi0[ind],withMichel);
    else if (MultiPi0) updateBackground(bckg_MultiPi0[ind],withMichel);
    else cout<<"WARNING! No BackgroundWithPi0 Found"<<endl;
}

void CCProtonPi0_BackgroundTool::fillBackground(int nProngs, bool NC, bool AntiNeutrino, bool QELike, bool SinglePion, bool DoublePion, bool MultiPion, bool Other, bool withMichel)
{
    int ind = setArrayInd(nProngs);
    
    // Update Total Background
    updateBackground(bckg_Total[ind],withMichel);
    
    // Update Each Background
    if (NC) updateBackground(bckg_NC[ind],withMichel);
    else if (AntiNeutrino) updateBackground(bckg_AntiNeutrino[ind],withMichel);
    else if (QELike) updateBackground(bckg_QELike[ind],withMichel);
    else if (SinglePion) updateBackground(bckg_SinglePion[ind],withMichel);
    else if (DoublePion) updateBackground(bckg_DoublePion[ind],withMichel);
    else if (MultiPion) updateBackground(bckg_MultiPion[ind],withMichel);
    else if (Other) updateBackground(bckg_Other[ind],withMichel);
    else cout<<"WARNING! No Background Type Found!"<<endl;
}

void CCProtonPi0_BackgroundTool::writeBackgroundTableHeader(int nProngs)
{
    int t = setArrayInd(nProngs);

    textFile[t]<<std::left;
    // Table Header
    textFile[t].width(20); textFile[t]<< "Background Type";
    textFile[t].width(16); textFile[t]<< "N(Events)"; 
    textFile[t].width(16); textFile[t]<< "Percent"; 
    textFile[t].width(16); textFile[t]<< "N(WithMichel)"; 
    textFile[t].width(16); textFile[t]<< "Percent"; 
    textFile[t]<<endl;
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
    BackgroundTypeVector.push_back(bckg_SinglePion);
    BackgroundTypeVector.push_back(bckg_DoublePion);
    BackgroundTypeVector.push_back(bckg_MultiPion);
    BackgroundTypeVector.push_back(bckg_Other);
    BackgroundTypeVector.push_back(bckg_Total);
}

double CCProtonPi0_BackgroundTool::calcPercent(double nEvents, double nBase)
{
    return (nEvents / nBase * 100.0);
}

void CCProtonPi0_BackgroundTool::writeBackgroundTableRows(vector< vector<Background> > &bckgVector, int nProngs)
{
    int t = setArrayInd(nProngs);
       
    // Write Background With Pi0
    for (unsigned int i = 0; i < bckgVector.size(); i++){
        Background temp = bckgVector[i][t];

        textFile[t].width(20); textFile[t]<<temp.name;  
        textFile[t].width(16); textFile[t]<<temp.nAll;    
        textFile[t].width(16); textFile[t]<<calcPercent(temp.nAll,bckg_Total[t].nAll);    
        textFile[t].width(16); textFile[t]<<temp.nWithMichel;   
        textFile[t].width(16); textFile[t]<<calcPercent(temp.nWithMichel,bckg_Total[t].nWithMichel);    
        textFile[t]<<endl;
    }

    textFile[t]<<endl;
}

void CCProtonPi0_BackgroundTool::writeBackgroundTable()
{
    formBackgroundVectors();

    for (int t = 0; t < nTopologies; t++){
        cout<<">> Writing "<<fileName[t]<<endl;

        writeBackgroundTableHeader(t);
        writeBackgroundTableRows(BackgroundWithPi0Vector, t);
        writeBackgroundTableRows(BackgroundTypeVector, t);
    }
}

void CCProtonPi0_BackgroundTool::OpenTextFiles()
{
    // Open Background Files
    fileName[0] = Folder_List::output + Folder_List::textOut + branchDir + "BackgroundTable_1Prong.txt";
    fileName[1] = Folder_List::output + Folder_List::textOut + branchDir + "BackgroundTable_2Prong.txt";
    
    for (int i = 0; i < nTopologies; i++){
        textFile[i].open( fileName[i].c_str() );
        if( !textFile[i].is_open() ){
            cerr<<"Cannot open output text file: "<<fileName[i]<<endl;
            exit(1);
        }else{
            cout<<"\t"<<fileName[i]<<endl;
        }
    }    
}

CCProtonPi0_BackgroundTool::~CCProtonPi0_BackgroundTool()
{
    for (int i = 0; i < nTopologies; i++){
        textFile[i].close(); 
    }
}



#endif
