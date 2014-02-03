/*
================================================================================
Class: CutNumberList
    Class which holds CutNumbers and functions related to cutNumbers 
    Uses SingleCutNumber Class to define Cut Numbers
    
    
    Main Directory:
        Classes/CutNumberList
    
    Usage:
        > #include "Classes/CutNumberList/CutNumberList.cpp" 
        > CutNumberList nCutList;
        > nCutList.nAll.increment();
        > percent = nCutList.getPercent(nCutList.nMuon.getValue());
    
    Last Revision: 2014_02_03
================================================================================
*/
#ifndef CutNumberList_h
#define CutNumberList_h

#include <iostream>
#include "Classes/SingleCutNumber/SingleCutNumber.cpp"
#include "Libraries/Folder_List.h"

const int initCount = 0.0;

class CutNumberList
{
    public:
        CutNumberList();
        
        // Returns the Percent of input number wrt All Events (No Cut)
        double getPercent(double nCurrent);
        void writeCutTable();
        void writeCut(SingleCutNumber* cut);
        
        SingleCutNumber* nAll;       // Number of All Events
        SingleCutNumber* nVolume;    // Number of Events After Volume Cut
        SingleCutNumber* nMuon;      // Number of Events After Muon Cut
        SingleCutNumber* nProton;    // Number of Events After Proton Cut
        SingleCutNumber* nPion;      // Number of Events After Pion Cut
        
    private:
        string cutFile;
        ofstream cutText;

};


#endif

