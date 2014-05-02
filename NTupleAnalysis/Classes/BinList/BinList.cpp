/*
    See BinList.h header for Class Information
*/
#include "BinList.h"

using namespace std;

BinList::BinList()
{
    // -------------------------------------------------------------------------
    //     Initialization
    //--------------------------------------------------------------------------
    // Standard Bins
    error.setBin(400, -2.0, 2.0);
    angle.setBin(180, 0.0, 180.0);
    trackLength.setBin(200, 0.0, 2000.0 );
    particleScore.setBin(20,0.0,1.0);
    multiplicity.setBin(10,0.0,10.0);
    int_channel.setBin(9,0.0,9.0);
    vertex_z.setBin(470,4300.0,9000.0);
    vertex_x_y.setBin(200,-1000.0,1000.0);
    
    // Neutrino Specific Bins
    beamE.setBin(20,0,20000.0);
    q2.setBin(60,0.0,3.0);
}

BinList::~BinList()
{

}




