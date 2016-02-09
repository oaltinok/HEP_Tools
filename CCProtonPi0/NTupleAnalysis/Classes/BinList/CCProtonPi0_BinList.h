/*
================================================================================
Class: CCProtonPi0_BinList
    List of the Bins used in Analysis
    Uses CCProtonPi0_SingleBin Class as Member Variable
    
    Main Directory:
        Classes/BinList
    
            
    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
================================================================================
*/
#ifndef CCProtonPi0_BinList_h
#define CCProtonPi0_BinList_h

#include "../SingleBin/CCProtonPi0_SingleBin.h"

class CCProtonPi0_BinList
{
    public:
        CCProtonPi0_BinList();
        
        // Standard Bins
        CCProtonPi0_SingleBin true_false;
        CCProtonPi0_SingleBin error;
        CCProtonPi0_SingleBin ratio;
        CCProtonPi0_SingleBin angle;
        CCProtonPi0_SingleBin particleScore;
        CCProtonPi0_SingleBin particleScore_LLR;
        CCProtonPi0_SingleBin particleScoreSum;
        CCProtonPi0_SingleBin particleScoreDiff;
        CCProtonPi0_SingleBin particleStatus;
        CCProtonPi0_SingleBin multiplicity;
        CCProtonPi0_SingleBin fraction;
        CCProtonPi0_SingleBin fraction2;
        CCProtonPi0_SingleBin vertex_z;
        CCProtonPi0_SingleBin vertex_x_y;
        CCProtonPi0_SingleBin deltaInvMass;
        CCProtonPi0_SingleBin preFilter_Status;
        CCProtonPi0_SingleBin preFilter_RejectedEnergy;
        CCProtonPi0_SingleBin strip_numbers;
        CCProtonPi0_SingleBin shower_length; // nPlanes in Shower
        CCProtonPi0_SingleBin gamma_evis;
        CCProtonPi0_SingleBin gamma_evis_pdg;
        CCProtonPi0_SingleBin pi0_evis_pdg;
        CCProtonPi0_SingleBin digit_E;
        CCProtonPi0_SingleBin kE;
        
        CCProtonPi0_SingleBin pi0_P;
        // Analysis
        CCProtonPi0_SingleBin beamE;
        CCProtonPi0_SingleBin beamE_Diff_True;
        CCProtonPi0_SingleBin beamE_Diff;
        CCProtonPi0_SingleBin q2;
        CCProtonPi0_SingleBin wSq;
        CCProtonPi0_SingleBin w;
        CCProtonPi0_SingleBin vertex_energy;
        CCProtonPi0_SingleBin vertex_evis;
        CCProtonPi0_SingleBin extra_energy;
        CCProtonPi0_SingleBin short_proton_P;
        CCProtonPi0_SingleBin short_proton_KE;
        

        CCProtonPi0_SingleBin UsedE;
        CCProtonPi0_SingleBin UnusedE;
        
        // Cut Histograms
        CCProtonPi0_SingleBin mc_w;
        CCProtonPi0_SingleBin eVis_nuclearTarget;
        CCProtonPi0_SingleBin eVis_other;
        CCProtonPi0_SingleBin pi0_invMass;
        CCProtonPi0_SingleBin bin_photonConvLength;
        
        // Michel Tool
        CCProtonPi0_SingleBin michelMuon_P;
        CCProtonPi0_SingleBin michelMuon_end_dist_vtx;
        CCProtonPi0_SingleBin michelMuon_length;
        CCProtonPi0_SingleBin michelMuon_Z_vtx;
        CCProtonPi0_SingleBin michelPion_P;
        CCProtonPi0_SingleBin michelPion_begin_dist_vtx;
        CCProtonPi0_SingleBin michelPion_length;
        CCProtonPi0_SingleBin michel_time_diff;
        
    
    private:
        

};

#endif
