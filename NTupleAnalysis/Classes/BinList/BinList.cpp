/*
    See BinList.h header for Class Information
*/

#ifndef BinList_cpp
#define BinList_cpp

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
    particleScore.setBin(100,0.0,1.0);
    particleStatus.setBin(17,-1,16);
    multiplicity.setBin(10,0.0,10.0);
    int_channel.setBin(9,0.0,9.0);
    vertex_z.setBin(470,4300.0,9000.0);
    vertex_x_y.setBin(200,-1000.0,1000.0);
    deltaInvMass.setBin(20,800.0,2800.0);
    preFilter_Status.setBin(4,0,4);
    preFilter_RejectedEnergy.setBin(50,0.0,5000.0);
    
    muonE.setBin(100,0,10.0);
    pionE.setBin(30, 0.0, 3.0);
    protonKE.setBin(40, -1.0, 3.0);
    
    // Neutrino Specific Bins
    beamE.setBin(200,0,20.0);
    q2.setBin(40,0.0,4.0);
    w2.setBin(150,0.0,3.0);
    w2fail.setBin(300,-3.0,3.0);
}

#endif




