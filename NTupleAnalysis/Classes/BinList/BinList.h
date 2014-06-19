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
            
    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision: 2014_06_18
================================================================================
*/
#ifndef BinList_h
#define BinList_h

#include "../SingleBin/SingleBin.h"

class BinList
{
    public:
        BinList();
        
        // Standard Bins
        SingleBin error;
        SingleBin angle;
        SingleBin trackLength;
        SingleBin particleScore;
        SingleBin multiplicity;
        SingleBin int_channel;
        SingleBin vertex_z;
        SingleBin vertex_x_y;
        SingleBin mgg_reco;
        
        // Analysis
        SingleBin beamE;
        SingleBin q2;
        
    
    private:
        

};

#endif
