/*
================================================================================
Class: CutNumberList
    Class which holds CutNumbers and functions related to cutNumbers 
    Uses CutNumber Class to define Cut Numbers
    
    
    Main Directory:
        Classes/CutNumberList
    
    Usage:
        > #include "Classes/CutNumberList/CutNumberList.cpp" 
        > CutNumberList nCutList;
        > nCutList.nAll.increment();
        > percent = nCutList.getPercent(nCutList.nMuon.getValue());
    
    Last Revision: 2014_03_19
================================================================================
*/
#ifndef CutNumberList_h
#define CutNumberList_h

#include <iostream>
#include "Classes/CutNumber/CutNumber.cpp"
#include "Libraries/Folder_List.h"

struct Node
{
    CutNumber* nCut;
    Node* next;
};

class CutNumberList
{
    public:
        CutNumberList();
        
        // Returns the Percent of input number wrt All Events (No Cut)
        double getPercent(double nCurrent);
        void writeCutTable();
        void addCutNumber(string label);
        void printList();
        CutNumber* getCutNumber(string targetLabel);
        
        // Linked List Implementation
        Node* head;
        Node* tail;
        
        // Predefine pointers for Default Cut Numbers
        // Increases Performance - No Search Needed for Default cut Numbers
        CutNumber* nAll;
        CutNumber* nFSPart;
        CutNumber* nVolume;
        CutNumber* nBeamEnergy;
        CutNumber* nMuon;        
        CutNumber* nMinos;
        CutNumber* nProton;                
        CutNumber* nPion;                
        CutNumber* nBeamEnergyFail;
        
    private:
        string cutFile;
        ofstream cutText;

};


#endif

