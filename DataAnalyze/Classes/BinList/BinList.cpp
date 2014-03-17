/*
    See BinList.h header for Class Information
*/
#include "BinList.h"

using namespace std;

BinList::BinList()
{
    // -------------------------------------------------------------------------
    //     Memory Allocation
    //--------------------------------------------------------------------------
    // Standard Bins
    error = new SingleBin;
    angle = new SingleBin;
    trackLength = new SingleBin;
    particleScore = new SingleBin;
    multiplicity = new SingleBin;
    int_channel = new SingleBin;
    
    // Neutrino Specific Bins
    neutrinoP = new SingleBin;
    
    // Muon Specific Bins
    muonP = new SingleBin;
    
    // Proton Specific Bins
    protonP = new SingleBin;
    
    // Pion Specific Bins
    pionP = new SingleBin;
    
    // -------------------------------------------------------------------------
    //     Initialization
    //--------------------------------------------------------------------------
    // Standard Bins
    error->setBin(400, -2.0, 2.0);
    angle->setBin(180, 0.0, 180.0);
    trackLength->setBin(200, 0.0, 2000.0 );
    particleScore->setBin(40, 0.0, 1.0);
    multiplicity->setBin(11,0.0,11.0);
    int_channel->setBin(9,0.0,9.0);
    
    // Neutrino Specific Bins
    neutrinoP->setBin(120,0.0,20000.0);
    
    // Muon Specific Bins
    muonP->setBin(100,0.0,10000.0);
    
    // Proton Specific Bins
    protonP->setBin(40, 0.0, 4000.0);
    
    // Pion Specific Bins
    pionP->setBin(30, 0.0, 3000.0);
}

BinList::~BinList()
{
        delete error;
        delete angle;
        delete trackLength;
        delete particleScore;
        delete multiplicity;
        delete int_channel;
        
        // Neutrino Specific Bins
        delete neutrinoP;
        
        // Muon Specific Bins
        delete muonP;
        
        // Proton Specific Bins
        delete protonP;
        
        // Pion Specific Bins
        delete pionP;
}




