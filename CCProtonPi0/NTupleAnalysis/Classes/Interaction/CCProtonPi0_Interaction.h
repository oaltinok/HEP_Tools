/*
================================================================================
Class: CCProtonPi0_Interaction
    CCProtonPi0_Interaction Class Contains interaction related Histograms
    All Histograms declared public and can be accessed by Analyzer
        Histograms:
            Event Kinematics
            Cut Histograms - May move these histograms to CutTool
            Study Histograms - Temporary (after the study will be removed)
                
    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
================================================================================
*/
#ifndef CCProtonPi0_Interaction_h
#define CCProtonPi0_Interaction_h

#include <iostream>
#include <string>
#include <cstdlib>

#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TFile.h>

#include "../NTupleAnalysis/CCProtonPi0_NTupleAnalysis.h"
#include "../BinList/CCProtonPi0_BinList.h"

using namespace std;

class CCProtonPi0_Interaction : public CCProtonPi0_NTupleAnalysis
{
    public:
        CCProtonPi0_Interaction(int nMode);
        
        void write_RootFile();
        
        //--------------------------------------------------------------------------
        //     Histograms
        //--------------------------------------------------------------------------
        // Event Kinematics
        TH1D* Enu_1Track_mc;
        TH1D* Enu_1Track_reco;
        TH1D* Enu_1Track_error;
        TH2D* Enu_1Track_reco_mc;
         
        TH1D* Enu_2Track_mc;
        TH1D* Enu_2Track_reco;
        TH1D* Enu_2Track_error;
        TH2D* Enu_2Track_reco_mc;
       
        TH1D* Enu_Cal_mc;
        TH1D* Enu_Cal_reco;
        TH1D* Enu_Cal_error;
        TH2D* Enu_Cal_reco_mc;
        
        TH2D* Enu_1Track_Enu_Cal;
        
        TH1D* deltaInvMass_mc;
        TH1D* deltaInvMass_reco;
        TH1D* deltaInvMass_error;
        TH2D* deltaInvMass_reco_mc;
        
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
       
        TH1D* nProngs_hist;
        TH1D* mc_w_DIS;
        TH1D* mc_w_RES;
        TH1D* mc_w_CCQE;
        
        TH1D* proton_p;
        TH1D* proton_p_shifted;
        
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
        

        // -------------------------------------------------------------------------
        // Cut Histograms
        // -------------------------------------------------------------------------
        // Common
        TH1D* hCut_vertexCount;
        TH1D* hCut_Michel;
        TH1D* hCut_eVis_nuclearTarget;
        TH1D* hCut_eVis_other;
        TH1D* hCut_pi0invMass;
        TH1D* hCut_gamma1ConvDist;
        TH1D* hCut_gamma2ConvDist;
        
        TH1D* hCut_1Prong_neutrinoE;
        TH1D* hCut_2Prong_neutrinoE;
        TH1D* hCut_1Prong_UnusedE;
        TH1D* hCut_2Prong_UnusedE;

        
        // 2 Prong Specific
        TH1D* hCut_protonScore_pIDDiff;
        TH1D* hCut_protonScore_LLR;
        TH1D* hCut_deltaInvMass;
    

    private:
        void initHistograms();
        
        TFile* f;
        string rootDir;
        
        CCProtonPi0_BinList binList;
        
    


};


#endif
