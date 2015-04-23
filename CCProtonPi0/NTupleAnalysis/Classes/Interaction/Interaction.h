/*
================================================================================
Class: Interaction
    Interaction Class Contains interaction related Histograms
    All Histograms declared public and can be accessed by Analyzer
        Histograms:
            Event Kinematics
            Cut Histograms - May move these histograms to CutTool
            Study Histograms - Temporary (after the study will be removed)
                
    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision: 2015_04_22
================================================================================
*/
#ifndef Interaction_h
#define Interaction_h

#include <iostream>
#include <string>
#include <cstdlib>

#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TFile.h>

#include "../NTupleAnalysis/NTupleAnalysis.h"
#include "../BinList/BinList.h"

using namespace std;

class Interaction : public NTupleAnalysis
{
    public:
        Interaction(int nMode);
        
        void write_RootFile();
        
        //--------------------------------------------------------------------------
        //     Histograms
        //--------------------------------------------------------------------------
        // Event Kinematics
        TH1D* beamEnergy_mc;
        TH1D* beamEnergy_reco;
        TH1D* beamEnergy_error;
        TH2D* beamEnergy_reco_mc;
        
        TH1D* beamEnergyCal_mc;
        TH1D* beamEnergyCal_reco;
        TH1D* beamEnergyCal_error;
        TH2D* beamEnergyCal_reco_mc;
        
        TH2D* beamEnergy_beamEnergyCal;
        
        TH1D* deltaInvMass_mc;
        TH1D* deltaInvMass_reco;
        TH1D* deltaInvMass_error;
        TH2D* deltaInvMass_reco_mc;
        
        TH1D* nProngs_hist;
        
        TH1D* pFilter_Status;
        TH1D* pFilter_RejectedEnergy;
        
        TH1D* q2_mc;
        TH1D* q2_reco;
        TH1D* q2_error;
        TH2D* q2_reco_mc;
        
        TH1D* w_mc;
        TH1D* w_reco;
        TH1D* w_error;
        TH2D* w_reco_mc;
        TH1D* wSq_reco;
        
        TH1D* vertex_count;
        
        TH1D* vertex_z_true;
        TH1D* vertex_z_reco;
        TH1D* vertex_z_error;
        TH2D* vertex_z_reco_mc;
    
        TH1D* int_channel;
        TH2D* vertex_x_y_true;
        TH2D* vertex_x_y_reco;
        
        TH1D* mc_w_DIS;
        TH1D* mc_w_RES;
        TH1D* mc_w_CCQE;
        
        TH1D* final_mc_w_DIS;
        TH1D* final_mc_w_RES;
        TH1D* final_mc_w_CCQE;
        
        TH1D* status_Pi0;
        TH1D* status_Pi0_Mother;
        TH1D* status_Pi0_GrandMother;
        
        TH1D* E_Unused_afterReco;
        TH1D* E_Used_afterReco;
        TH1D* time_AllClusters;
        
        TH1D* total_E;
        TH2D* total_E_neutrinoE;
        
        // Michel Study
        TH1D* N_michelElectrons;
        TH1D* michelMuon_P[4];
        TH1D* michelMuon_end_dist_vtx[4];
        TH1D* michelMuon_length[4];
        TH1D* michelMuon_Z[4];
        TH2D* michelMuon_X_Y[4];
        TH1D* michelMuon_Z_vtx[4];
        TH1D* michelPion_P[4];
        TH1D* michelPion_begin_dist_vtx[4];
        TH1D* michelPion_length[4];
        TH1D* michelElectron_E[5];
        TH2D* michelPion_length_dist_vtx[4];
        TH2D* michelMuon_dist_michelPion_length[4];
        TH2D* michel_dist_reco_true;
        
        TH1D* trueMichel_dist_reco;
        TH1D* fakeMichel_dist_reco;
        TH1D* trueMichel_end_Z;
        TH1D* fakeMichel_end_Z;
        TH1D* trueMichel_end_Z_vtx_Z;
        TH1D* fakeMichel_end_Z_vtx_Z;
        TH1D* trueMichel_energy;
        TH1D* fakeMichel_energy;
        TH1D* trueMichel_time_diff;
        TH1D* fakeMichel_time_diff;

        
        // -------------------------------------------------------------------------
        // Cut Histograms
        // -------------------------------------------------------------------------
        // Common
        TH1D* hCut_vertexCount;
        TH1D* hCut_nProngs;
        
        // Topology Dependent
        TH1D* hCut_1Prong_Michel;
        TH1D* hCut_2Prong_Michel;
        TH1D* hCut_1Prong_eVis_nuclearTarget;
        TH1D* hCut_2Prong_eVis_nuclearTarget;
        TH1D* hCut_1Prong_eVis_other;
        TH1D* hCut_2Prong_eVis_other;
        TH1D* hCut_1Prong_pi0invMass;
        TH1D* hCut_2Prong_pi0invMass;
        TH1D* hCut_1Prong_gamma1ConvDist;
        TH1D* hCut_2Prong_gamma1ConvDist;
        TH1D* hCut_1Prong_gamma2ConvDist;
        TH1D* hCut_2Prong_gamma2ConvDist;
        TH1D* hCut_1Prong_neutrinoE;
        TH1D* hCut_2Prong_neutrinoE;
        TH1D* hCut_1Prong_UnusedE;
        TH1D* hCut_2Prong_UnusedE;

        
        // 2 Prong Specific
        TH1D* hCut_pIDDiff;
        TH1D* hCut_protonScore_LLR;
        TH1D* hCut_deltaInvMass;
    

    private:
        void initHistograms();
        
        TFile* f;
        string rootDir;
        
        BinList binList;
        
    


};


#endif
