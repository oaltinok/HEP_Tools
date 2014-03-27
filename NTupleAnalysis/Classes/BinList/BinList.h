/*
================================================================================
Class: BinList
    List of the Bins used in Analysis
    Uses SingleBin Class as Member Variable
    
    Main Directory:
        Classes/BinList
    
    Usage:
        > #include "Classes/BinList/BinList.h" 
        > BinList* binList;
        > binList->error.get_nBins();
        > binList->error.get_min();
        > binList->error.get_max();
            
    
    Last Revision: 2014_03_17
================================================================================
*/
#ifndef BinList_h
#define BinList_h

#include "Classes/SingleBin/SingleBin.cpp"

class BinList
{
    public:
        BinList();
        ~BinList();
        
        // Standard Bins
        SingleBin error;
        SingleBin angle;
        SingleBin trackLength;
        SingleBin particleScore;
        SingleBin multiplicity;
        SingleBin int_channel;
        SingleBin vertex_z;
        SingleBin vertex_x_y;
        
        // Analysis
        SingleBin beamE;
        SingleBin q2;
        
    
    private:
        

};

#endif
