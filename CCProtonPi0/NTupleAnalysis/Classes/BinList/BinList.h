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
    Last Revision: 2014_12_10
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
        SingleBin objectCount;
        SingleBin error;
        SingleBin angle;
        SingleBin trackLength;
        SingleBin particleScore;
        SingleBin particleScoreSum;
        SingleBin particleScoreDiff;
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
        SingleBin wSq;
        SingleBin w;
        SingleBin wfail;
        
        SingleBin UsedE;
        SingleBin UnusedE;
        SingleBin time;
        
        // Cut Histograms
        SingleBin eVis_nuclearTarget;
        SingleBin eVis_other;
        SingleBin michelID;
        SingleBin pi0_invMass;
        SingleBin bin_photonConvLength;
        
    
    private:
        

};

#endif
