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
    cutFile = cutFile + "_" + channelTag + ".txt";
    cutText.open( cutFile.c_str() );
    if( !cutText.is_open() ){
        cerr<<"Cannot open output text file: "<<cutFile<<endl;
        exit(1);
    }
    
    // Initialize Linked List
    head = NULL;
    tail = NULL;
   
    // Add Default CutNumbers to the Linked List
    addCutNumber("nAll");
    addCutNumber("nFSPart");
    addCutNumber("nVolume");
    addCutNumber("nBeamEnergy");
    addCutNumber("nMuon");
    addCutNumber("nMinos");
    addCutNumber("nProton");
    addCutNumber("nPion");
    addCutNumber("nBeamEnergyFail");
    
    
    // Initialize Default CutNumber Pointers
    nAll = getCutNumber("nAll");
    nFSPart = getCutNumber("nFSPart");
    nVolume = getCutNumber("nVolume");
    nBeamEnergy = getCutNumber("nBeamEnergy");
    nMuon = getCutNumber("nMuon");
    nMinos = getCutNumber("nMinos");
    nProton = getCutNumber("nProton");
    nPion = getCutNumber("nPion");
    nBeamEnergyFail = getCutNumber("nBeamEnergyFail");
    
}

bool CutNumberList::addCutNumber(string label)
{
    Node* temp = new Node;
    temp->nCut  = new CutNumber(label);
    
    // Empty List
    if (head == NULL){
        head = temp;
        head->next = NULL;
        tail = head;
        return true;
    }
    // 1 Item in the List
    else if (head->next == NULL){
        head->next = temp;
        temp->next = NULL;
        tail = temp;
        return true;
    }
    // Full List
    else{
        temp->next = NULL;
        tail->next = temp;
        tail = temp;
        return true;
    }
    
    return false;

}

CutNumber* CutNumberList::getCutNumber(string targetLabel)
{
    Node* temp;
    temp = head;
    
    while(  temp != NULL ){
        if( temp->nCut->getLabel() == targetLabel){
            return temp->nCut;
        }else{
            temp = temp->next;
        }
    }
        
    return NULL;
}

void CutNumberList::printList()
{
    Node* temp;
    temp = head;
    
    cout<<"\t------ Cut Number List ------"<<endl;
    while( temp != NULL){
        cout<<"\t"<<temp->nCut->getLabel()<<endl;
        temp = temp->next;
    }
    cout<<"\t-----------------------------"<<endl;
}


double CutNumberList::getPercent(double nCurrent)
{
    double percent;
    
    if(head->nCut->getValue() == 0 ){
        return 0;
    }else{
        percent = ( nCurrent / head->nCut->getValue() ) * 100;
        return percent;
    }
}

void CutNumberList::writeCutTable()
{
    Node* temp = head;
    string label;
    double value;

    cout<<">> Writing "<<cutFile<<endl;
    
    while( temp != NULL ){
        label = temp->nCut->getLabel();
        value = temp->nCut->getValue();
        cutText<<label<<"\t|\t"<<value<<"\t|\t"<<getPercent(value)<<"%"<<endl;
        temp = temp->next;
    }
    
    cutText.close();
}


