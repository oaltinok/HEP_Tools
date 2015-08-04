/*
    See CCProtonPi0_BinList.h header for Class Information
*/

#ifndef CCProtonPi0_BinList_cpp
#define CCProtonPi0_BinList_cpp

#include "CCProtonPi0_BinList.h"

using namespace std;

CCProtonPi0_BinList::CCProtonPi0_BinList()
{
    // -------------------------------------------------------------------------
    //     Initialization
    //--------------------------------------------------------------------------
    // Standard Bins
    error.setBin(40, -2.0, 2.0);
    angle.setBin(18, 0.0, 180.0);
    particleScore.setBin(100,0.0,1.0);
    particleScore_LLR.setBin(40,-100.0,100.0);
    particleScoreSum.setBin(200,0.0,2.0);
    particleScoreDiff.setBin(100,-1.0,1.0);
    particleStatus.setBin(17,-1,16);
    multiplicity.setBin(10,0.0,10.0);
    vertex_z.setBin(94,4300.0,9000.0);
    vertex_x_y.setBin(200,-1000.0,1000.0);
    deltaInvMass.setBin(20,800.0,2800.0);
    preFilter_Status.setBin(4,0,4);
    preFilter_RejectedEnergy.setBin(50,0.0,5000.0);

    // Neutrino Specific Bins
    beamE.setBin(20,0,20.0);
    q2.setBin(40,0.0,4.0);
    w.setBin(30,0.0,3.0);
    wSq.setBin(50,0.0,5.0);
    
    UsedE.setBin(40,0.0,2000.0);
    UnusedE.setBin(20,0.0,1000.0);
    
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
    michel_time_diff.setBin(100,0.0,10000.0);
}

#endif




