#ifndef CCProtonPi0_BackgroundTool_cpp
#define CCProtonPi0_BackgroundTool_cpp

#include "CCProtonPi0_BackgroundTool.h"

using namespace std;

CCProtonPi0_BackgroundTool::CCProtonPi0_BackgroundTool(int nMode) : CCProtonPi0_NTupleAnalysis(nMode)
{
    cout<<"Initializing CCProtonPi0_BackgroundTool"<<endl;
    
    OpenTextFiles();
    
    // Init Background Vectors
    for(int i = 0; i < nBckgBranch; i++){
            bckg_1Prong_QELike.push_back(0.0);
            bckg_1Prong_SinglePiPlus.push_back(0.0);
            bckg_1Prong_SinglePiMinus.push_back(0.0);
            bckg_1Prong_MultiPion.push_back(0.0);
            bckg_1Prong_MultiPiZero.push_back(0.0);
            bckg_1Prong_Other.push_back(0.0);
            bckg_2Prong_QELike.push_back(0.0);
            bckg_2Prong_SinglePiPlus.push_back(0.0);
            bckg_2Prong_SinglePiMinus.push_back(0.0);
            bckg_2Prong_MultiPion.push_back(0.0);
            bckg_2Prong_MultiPiZero.push_back(0.0);
            bckg_2Prong_Other.push_back(0.0);
    }
    
    cout<<"Done!"<<endl;
    
}

void CCProtonPi0_BackgroundTool::fillBackgroundBranches(    int nProngs,
                                                bool truth_isBckg_QELike, 
                                                bool truth_isBckg_SinglePiPlus, 
                                                bool truth_isBckg_SinglePiMinus, 
                                                bool truth_isBckg_MultiPion, 
                                                bool truth_isBckg_MultiPiZero, 
                                                bool truth_isBckg_Other,
                                                bool truth_isBckg_withAntiMuon,
                                                bool truth_isBckg_withMichel,
                                                bool truth_isBckg_withPrimaryPi0,
                                                bool truth_isBckg_withSecondaryPi0)
{
    if (nProngs == 1){
        if(truth_isBckg_QELike){ 
            setBackgroundBranch(    bckg_1Prong_QELike,
                                    truth_isBckg_withAntiMuon,
                                    truth_isBckg_withMichel,
                                    truth_isBckg_withPrimaryPi0,
                                    truth_isBckg_withSecondaryPi0);
        }
        
        if(truth_isBckg_SinglePiPlus){ 
            setBackgroundBranch(    bckg_1Prong_SinglePiPlus,
                                    truth_isBckg_withAntiMuon,
                                    truth_isBckg_withMichel,
                                    truth_isBckg_withPrimaryPi0,
                                    truth_isBckg_withSecondaryPi0);
        }
        
        if(truth_isBckg_SinglePiMinus){ 
            setBackgroundBranch(    bckg_1Prong_SinglePiMinus,
                                    truth_isBckg_withAntiMuon,
                                    truth_isBckg_withMichel,
                                    truth_isBckg_withPrimaryPi0,
                                    truth_isBckg_withSecondaryPi0);
        }
        
        if(truth_isBckg_MultiPion){ 
            setBackgroundBranch(    bckg_1Prong_MultiPion,
                                    truth_isBckg_withAntiMuon,
                                    truth_isBckg_withMichel,
                                    truth_isBckg_withPrimaryPi0,
                                    truth_isBckg_withSecondaryPi0);
        }
        
        if(truth_isBckg_MultiPiZero){ 
            setBackgroundBranch(    bckg_1Prong_MultiPiZero,
                                    truth_isBckg_withAntiMuon,
                                    truth_isBckg_withMichel,
                                    truth_isBckg_withPrimaryPi0,
                                    truth_isBckg_withSecondaryPi0);
        }
        
        if(truth_isBckg_Other){ 
            setBackgroundBranch(    bckg_1Prong_Other,
                                    truth_isBckg_withAntiMuon,
                                    truth_isBckg_withMichel,
                                    truth_isBckg_withPrimaryPi0,
                                    truth_isBckg_withSecondaryPi0);
        }
    }else if (nProngs == 2){
            if(truth_isBckg_QELike){ 
            setBackgroundBranch(    bckg_2Prong_QELike,
                                    truth_isBckg_withAntiMuon,
                                    truth_isBckg_withMichel,
                                    truth_isBckg_withPrimaryPi0,
                                    truth_isBckg_withSecondaryPi0);
        }
        
        if(truth_isBckg_SinglePiPlus){ 
            setBackgroundBranch(    bckg_2Prong_SinglePiPlus,
                                    truth_isBckg_withAntiMuon,
                                    truth_isBckg_withMichel,
                                    truth_isBckg_withPrimaryPi0,
                                    truth_isBckg_withSecondaryPi0);
        }
        
        if(truth_isBckg_SinglePiMinus){ 
            setBackgroundBranch(    bckg_2Prong_SinglePiMinus,
                                    truth_isBckg_withAntiMuon,
                                    truth_isBckg_withMichel,
                                    truth_isBckg_withPrimaryPi0,
                                    truth_isBckg_withSecondaryPi0);
        }
        
        if(truth_isBckg_MultiPion){ 
            setBackgroundBranch(    bckg_2Prong_MultiPion,
                                    truth_isBckg_withAntiMuon,
                                    truth_isBckg_withMichel,
                                    truth_isBckg_withPrimaryPi0,
                                    truth_isBckg_withSecondaryPi0);
        }
        
        if(truth_isBckg_MultiPiZero){ 
            setBackgroundBranch(    bckg_2Prong_MultiPiZero,
                                    truth_isBckg_withAntiMuon,
                                    truth_isBckg_withMichel,
                                    truth_isBckg_withPrimaryPi0,
                                    truth_isBckg_withSecondaryPi0);
        }
        
        if(truth_isBckg_Other){ 
            setBackgroundBranch(    bckg_2Prong_Other,
                                    truth_isBckg_withAntiMuon,
                                    truth_isBckg_withMichel,
                                    truth_isBckg_withPrimaryPi0,
                                    truth_isBckg_withSecondaryPi0);
        }
    }
}



void CCProtonPi0_BackgroundTool::setBackgroundBranch(vector<double>& background, bool hasAntiMuon, bool hasMichel, bool hasPrimaryPi0, bool hasSecondaryPi0 )
{
    bool hasPi0Candidate = hasPrimaryPi0 || hasSecondaryPi0;
    background[0]++;                        // Total
    if(hasAntiMuon) background[1]++;        // AntiMuon
    if(hasMichel) background[2]++;          // Michel
    if(hasPrimaryPi0) background[3]++;      // Primary Pi0
    if(hasSecondaryPi0) background[4]++;    // Secondary Pi0
    if(hasPi0Candidate) background[5]++;    // Pi0 Candidate
    
}

void CCProtonPi0_BackgroundTool::formBackgroundVectors()
{
    bckgVector_1Prong.push_back(bckg_1Prong_QELike);
    bckgVector_1Prong.push_back(bckg_1Prong_SinglePiPlus); 
    bckgVector_1Prong.push_back(bckg_1Prong_SinglePiMinus); 
    bckgVector_1Prong.push_back(bckg_1Prong_MultiPion); 
    bckgVector_1Prong.push_back(bckg_1Prong_MultiPiZero);
    bckgVector_1Prong.push_back(bckg_1Prong_Other);
    
    bckgVector_2Prong.push_back(bckg_2Prong_QELike);
    bckgVector_2Prong.push_back(bckg_2Prong_SinglePiPlus); 
    bckgVector_2Prong.push_back(bckg_2Prong_SinglePiMinus); 
    bckgVector_2Prong.push_back(bckg_2Prong_MultiPion); 
    bckgVector_2Prong.push_back(bckg_2Prong_MultiPiZero);
    bckgVector_2Prong.push_back(bckg_2Prong_Other);
}

void CCProtonPi0_BackgroundTool::writeBackgroundTableHeader()
{
    
        for (int t = 0; t < nTopologies; t++){ 
            backgroundText[t]<<std::left;
            // Table Header
            backgroundText[t].width(20); backgroundText[t]<< "Background Type"; backgroundText[t]<<"|";
            backgroundText[t].width(16); backgroundText[t]<< "Total"; backgroundText[t]<<"|";
            backgroundText[t].width(16); backgroundText[t]<< "AntiMuon"; backgroundText[t]<<"|";
            backgroundText[t].width(16); backgroundText[t]<< "Michel"; backgroundText[t]<<"|";
            backgroundText[t].width(16); backgroundText[t]<< "Primary Pi0"; backgroundText[t]<<"|";
            backgroundText[t].width(16); backgroundText[t]<< "Secondary Pi0"; backgroundText[t]<<"|";
            backgroundText[t].width(16); backgroundText[t]<< "Pi0 Candidate"; backgroundText[t]<<"|";
            backgroundText[t]<<endl;
        }
}


void CCProtonPi0_BackgroundTool::writeBackgroundTableRows(vector< vector<double> > bckgVector, int nProngs)
{
    int t;
    if( nProngs == 1 ) t = 0;
    if( nProngs == 2 ) t = 1;
    
    vector<double> bckg_QELike = bckgVector[0];
    vector<double> bckg_SinglePiPlus = bckgVector[1];
    vector<double> bckg_SinglePiMinus = bckgVector[2];
    vector<double> bckg_MultiPion = bckgVector[3];
    vector<double> bckg_MultiPiZero = bckgVector[4];
    vector<double> bckg_Other = bckgVector[5];
    
    backgroundText[t].width(20); backgroundText[t]<<"QE Like"; backgroundText[t]<<"|"; 
    for(int i = 0; i < nBckgBranch; i++){
        backgroundText[t].width(16); backgroundText[t]<<bckg_QELike[i]; backgroundText[t]<<"|";   
    }
    backgroundText[t]<<endl;
    
    backgroundText[t].width(20); backgroundText[t]<<"SinglePiPlus"; backgroundText[t]<<"|"; 
    for(int i = 0; i < nBckgBranch; i++){
        backgroundText[t].width(16); backgroundText[t]<<bckg_SinglePiPlus[i]; backgroundText[t]<<"|";    
    }
    backgroundText[t]<<endl;
    
    backgroundText[t].width(20); backgroundText[t]<<"SinglePiMinus"; backgroundText[t]<<"|"; 
    for(int i = 0; i < nBckgBranch; i++){
        backgroundText[t].width(16); backgroundText[t]<<bckg_SinglePiMinus[i]; backgroundText[t]<<"|";    
    }
    backgroundText[t]<<endl;

    backgroundText[t].width(20); backgroundText[t]<<"MultiChargedPion"; backgroundText[t]<<"|"; 
    for(int i = 0; i < nBckgBranch; i++){
        backgroundText[t].width(16); backgroundText[t]<<bckg_MultiPion[i]; backgroundText[t]<<"|";    
    }
    backgroundText[t]<<endl;
    
    backgroundText[t].width(20); backgroundText[t]<<"MultiPiZero"; backgroundText[t]<<"|"; 
    for(int i = 0; i < nBckgBranch; i++){
        backgroundText[t].width(16); backgroundText[t]<<bckg_MultiPiZero[i]; backgroundText[t]<<"|";    
    }
    backgroundText[t]<<endl;
    
    backgroundText[t].width(20); backgroundText[t]<<"Other"; backgroundText[t]<<"|"; 
    for(int i = 0; i < nBckgBranch; i++){
        backgroundText[t].width(16); backgroundText[t]<<bckg_Other[i]; backgroundText[t]<<"|";  
    }
    backgroundText[t]<<endl;
    
}

void CCProtonPi0_BackgroundTool::writeBackgroundTable()
{
    formBackgroundVectors();
    for (int t = 0; t < nTopologies; t++){
        cout<<">> Writing "<<backgroundFile[t]<<endl;
    }
    writeBackgroundTableHeader();
    writeBackgroundTableRows(bckgVector_1Prong, 1);
    writeBackgroundTableRows(bckgVector_2Prong, 2);
}

void CCProtonPi0_BackgroundTool::OpenTextFiles()
{
    // Open Background Files
    backgroundFile[0] = Folder_List::output + Folder_List::textOut + branchDir + "BackgroundTable_1Prong.txt";
    backgroundFile[1] = Folder_List::output + Folder_List::textOut + branchDir + "BackgroundTable_2Prong.txt";
    
    for (int i = 0; i < nTopologies; i++){
        backgroundText[i].open( backgroundFile[i].c_str() );
        if( !backgroundText[i].is_open() ){
            cerr<<"Cannot open output text file: "<<backgroundFile[i]<<endl;
            exit(1);
        }else{
            cout<<"\t"<<backgroundFile[i]<<endl;
        }
    }    
}

CCProtonPi0_BackgroundTool::~CCProtonPi0_BackgroundTool()
{
    for (int i = 0; i < nTopologies; i++){
        backgroundText[i].close(); 
    }
}



#endif
