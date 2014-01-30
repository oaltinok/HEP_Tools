/*
    See BinList.h header for Class Information
*/
#include "BinList.h"

using namespace std;

BinList::BinList()
{
    // Standard Bins
    error.setBin(400, -2.0, 2.0);
    angle.setBin(180, 0.0, 180.0);
    trackLength.setBin(200, 0.0, 2000.0 );
    particleScore.setBin(40, 0.0, 1.0);
    multiplicity.setBin(11,0.0,11.0);
    
    // Neutrino Specific Bins
    neutrinoP.setBin(50,0.0,5000.0);
    
    // Muon Specific Bins
    muonP.setBin(100,0.0,10000.0);
    
    // Proton Specific Bins
    protonP.setBin(40, 0.0, 4000.0);
    
    // Pion Specific Bins
    pionP.setBin(30, 0.0, 3000.0);
}




