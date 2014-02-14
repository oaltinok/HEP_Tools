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
    
    Last Revision: 2014_02_14
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
        
        
        Node* head;
        Node* tail;
        
    private:
        string cutFile;
        ofstream cutText;

};


#endif

