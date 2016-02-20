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
    true_false.setBin(2,0.0,2.0);
    error.setBin(80, -1.0, 1.0);
    ratio.setBin(50, 0.0, 4.0);
    angle.setBin(18, 0.0, 180.0);
    particleScore.setBin(100,0.0,1.0);
    particleScore_LLR.setBin(40,-100.0,100.0);
    particleScoreSum.setBin(200,0.0,2.0);
    particleScoreDiff.setBin(100,-1.0,1.0);
    particleStatus.setBin(17,-1,16);
    multiplicity.setBin(15,0.0,15.0);
    fraction.setBin(11,0,1.1);
    fraction2.setBin(200,-0.5,1.5);
    vertex_z.setBin(114,4300.0,10000.0);
    vertex_x_y.setBin(200,-1000.0,1000.0);
    deltaInvMass.setBin(30,0.0,3.0);
    preFilter_Status.setBin(4,0,4);
    preFilter_RejectedEnergy.setBin(50,0.0,5000.0);
    strip_numbers.setBin(127,1,127);
    shower_length.setBin(100,1,100);
    gamma_evis.setBin(50,0,300);
    gamma_evis_pdg.setBin(30,0,300);
    kE.setBin(50,1,5);
    pi0_evis_pdg.setBin(80,0,300);
    digit_E.setBin(50,0,50);

    pi0_P.setBin(17,0.0,1.7);
    // Event Kinematics Bins
    beamE.setBin(20,0,20.0);
    beamE_Diff_True.setBin(100,1.5,7.5);
    beamE_Diff.setBin(100,-3,3);
    q2.setBin(40,0.0,4.0);
    w.setBin(30,0.0,3.0);
    wSq.setBin(50,0.0,5.0);
    vertex_energy.setBin(10,0.0,500.0);
    vertex_evis.setBin(10,0.0,500.0);
    extra_energy.setBin(50,0.0,20.0);
    short_proton_P.setBin(50,0.0,1000.0);
    short_proton_KE.setBin(50,0.0,500.0);
    
    UsedE.setBin(40,0.0,2000.0);
    UnusedE.setBin(20,0.0,1000.0);
    
    // Cut Histograms
    mc_w.setBin(100,0.8,3.0);
    eVis_nuclearTarget.setBin(50,0.0,25.0);
    eVis_other.setBin(60,0.0,3000.0);
    pi0_invMass.setBin(50,0.0,500.0);
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




