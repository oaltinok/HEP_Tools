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
    objectCount.setBin(10,0.0,10.0);
    error.setBin(400, -2.0, 2.0);
    angle.setBin(180, 0.0, 180.0);
    trackLength.setBin(200, 0.0, 2000.0 );
    particleScore.setBin(100,0.0,1.0);
    particleScore_LLR.setBin(100,-100.0,100.0);
    particleScoreSum.setBin(200,0.0,2.0);
    particleScoreDiff.setBin(100,-1.0,1.0);
    particleStatus.setBin(17,-1,16);
    multiplicity.setBin(10,0.0,10.0);
    int_channel.setBin(9,0.0,9.0);
    vertex_z.setBin(94,4300.0,9000.0);
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
    w.setBin(150,0.0,3.0);
    wSq.setBin(300,-5.0,5.0);
    wfail.setBin(300,-3.0,3.0);
    
    UsedE.setBin(40,0.0,2000.0);
    UnusedE.setBin(20,0.0,1000.0);
    time.setBin(1000,0.0,10000.0);
    
    // Cut Histograms
    eVis_nuclearTarget.setBin(50,0.0,25.0);
    eVis_other.setBin(60,0.0,3000.0);
    michelID.setBin(2,0.0,2.0);
    pi0_invMass.setBin(60,0.0,600.0);
    bin_photonConvLength.setBin(50,0.0,100.0);
    
    // Michel Tool
    michelMuon_P.setBin(50,0.0,250.0);
    michelMuon_end_dist_vtx.setBin(40,0.0,2000.0);
    michelMuon_length.setBin(25,0.0,50.0);
    michelMuon_Z_vtx.setBin(80,-2000.0,2000.0);
    michelPion_P.setBin(100,0.0,2000.0);
    michelPion_begin_dist_vtx.setBin(40,0.0,2000.0);
    michelPion_length.setBin(40,0.0,2000.0);

}

#endif




