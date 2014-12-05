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
    Last Revision: 2014_11_25
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
        SingleBin particleStatus;
        SingleBin multiplicity;
        SingleBin int_channel;
        SingleBin vertex_z;
        SingleBin vertex_x_y;
        SingleBin deltaInvMass;
        SingleBin preFilter_Status;
        SingleBin preFilter_RejectedEnergy;
        
        SingleBin muonE;
        SingleBin pionE;
        SingleBin protonKE;
        
        
        // Analysis
        SingleBin beamE;
        SingleBin q2;
        SingleBin w2;
        SingleBin w2fail;
        
    
    private:
        

};

#endif
