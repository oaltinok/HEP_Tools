/*
================================================================================
Class: BinList
    List of the Bins used in Analysis
    Uses SingleBin Class as Member Variable
    
    Main Directory:
        Classes/BinList
    
    Usage:
        > #include "Classes/BinList/BinList.h" 
        > BinList binList;
        > binList.error.get_nBins();
        > binList.error.get_min();
        > binList.error.get_max();
            
    
    Last Revision: 2014_01_30
================================================================================
*/
#ifndef BinList_h
#define BinList_h

#include "Classes/SingleBin/SingleBin.cpp"

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
        
        // Neutrino Specific Bins
        SingleBin neutrinoP;
        
        // Muon Specific Bins
        SingleBin muonP;
        
        // Proton Specific Bins
        SingleBin protonP;
        
        // Pion Specific Bins
        SingleBin pionP;
    
    private:
        
        
    


};

#endif
