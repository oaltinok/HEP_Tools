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

#include "../NTupleAnalysis/CCProtonPi0_NTupleAnalysis.h"
#include "../BinList/CCProtonPi0_BinList.h"

using namespace PlotUtils;

class CCProtonPi0_Interaction : public CCProtonPi0_NTupleAnalysis
{
    public:
        CCProtonPi0_Interaction(bool isModeReduce, bool isMC);
        void writeHistograms();
       
        bool isErrHistFilled_NeutronResponse;
        bool isErrHistFilled_PionResponse;
        bool isErrHistFilled_MuonTracking;

        //--------------------------------------------------------------------------
        //     Histograms
        //--------------------------------------------------------------------------
        std::vector<MnvH1D*> CV_weight;
        std::vector<MnvH1D*> CV_weight_Flux;
        std::vector<MnvH1D*> CV_weight_2p2h;
        std::vector<MnvH1D*> CV_weight_Delta;
        std::vector<MnvH1D*> CV_weight_CCRES;
        std::vector<MnvH1D*> CV_weight_NonRes1pi;
        std::vector<MnvH1D*> err_2p2h;
        std::vector<MnvH1D*> genie_wgt_VecFFCCQEshape;
        std::vector<MnvH1D*> genie_wgt_NormDISCC;
        std::vector<MnvH1D*> genie_wgt_Theta_Delta2Npi;
        std::vector<MnvH1D*> updated_wgt_Theta_Delta2Npi;
        std::vector<MnvH1D*> genie_wgt_MaRES;
        std::vector<MnvH1D*> updated_wgt_MaRES;
        std::vector<MnvH1D*> genie_wgt_MvRES;
        std::vector<MnvH1D*> updated_wgt_MvRES;
        std::vector<MnvH1D*> genie_wgt_Rvn1pi;
        std::vector<MnvH1D*> updated_wgt_Rvn1pi;

        // Event Kinematics
        std::vector<MnvH1D*> Enu_1Track;
        std::vector<MnvH1D*> Enu_2Track;
        std::vector<MnvH1D*> Enu;
        std::vector<MnvH1D*> QSq;
        std::vector<MnvH1D*> QSq_1Track;
        std::vector<MnvH1D*> QSq_2Track;
        std::vector<MnvH1D*> WSq;
        std::vector<MnvH1D*> WSq_1Track;
        std::vector<MnvH1D*> WSq_2Track;
        std::vector<MnvH1D*> W;
        std::vector<MnvH1D*> W_1Track;
        std::vector<MnvH1D*> W_2Track;
  
        // Extra Energy
        std::vector<MnvH1D*> extra_leftover_energy_1Track;
        std::vector<MnvH1D*> extra_muon_energy_1Track;
        std::vector<MnvH1D*> extra_rejected_energy_1Track;
        std::vector<MnvH1D*> extra_total_energy_1Track;
        
        std::vector<MnvH1D*> extra_leftover_energy_2Track;
        std::vector<MnvH1D*> extra_muon_energy_2Track;
        std::vector<MnvH1D*> extra_rejected_energy_2Track;
        std::vector<MnvH1D*> extra_total_energy_2Track;

        // Background Subtraction for Studies
        std::vector<MnvH1D*> pi0_invMass_All; 
        std::vector<MnvH1D*> pi0_invMass_1Track;
        std::vector<MnvH1D*> pi0_invMass_2Track;
        std::vector<MnvH1D*> pi0_invMass_DeltaRES;

        // W Study
        std::vector<MnvH1D*> W_p_pi0;
        std::vector<MnvH1D*> W_All;
        std::vector<MnvH1D*> W_1;
        std::vector<MnvH1D*> W_2;
        std::vector<MnvH1D*> W_Shift;
        std::vector<MnvH1D*> W_Shift_Bckg;
        std::vector<MnvH1D*> W_Shift_Signal;

        // QSq Study
        std::vector<MnvH1D*> QSq_CV;
        std::vector<MnvH1D*> QSq_MaRES; // With Vert Error Bands

        // 2p2h Study
        std::vector<MnvH1D*> vertex_energy_All;
        std::vector<MnvH1D*> vertex_energy_1Track;
        std::vector<MnvH1D*> vertex_energy_2Track;

        std::vector<MnvH1D*> vertex_evis_All;
        std::vector<MnvH1D*> vertex_evis_1Track;
        std::vector<MnvH1D*> vertex_evis_2Track;

        std::vector<MnvH2D*> muon_theta_muon_KE;
        std::vector<MnvH2D*> q3_q0;
        std::vector<MnvH2D*> W_QSq;

        // Flux Study
        TH2D* Enu_flux_wgt;
        TH2D* Enu_cvweight;

        // QSq
        MnvH1D* QSq_all;
        MnvH1D* QSq_mc_reco_all;
        MnvH1D* QSq_mc_reco_signal;
        MnvH1D* QSq_mc_reco_bckg;
        MnvH1D* QSq_mc_truth_signal;
        MnvH2D* QSq_response;

        // Enu
        MnvH1D* Enu_all;
        MnvH1D* Enu_mc_reco_all;
        MnvH1D* Enu_mc_reco_signal;
        MnvH1D* Enu_mc_reco_bckg;
        MnvH1D* Enu_mc_truth_signal;
        MnvH2D* Enu_response;

        // W
        MnvH1D* W_all;
        MnvH1D* W_mc_reco_all;
        MnvH1D* W_mc_reco_signal;
        MnvH1D* W_mc_reco_bckg;
        MnvH1D* W_mc_truth_signal;
        MnvH2D* W_response;

        // Vertical Error Band -- Errors 
        TH1D* Err_NeutronResponse;  
        TH1D* Err_PionResponse;  
        TH1D* Err_MuonTracking;  

        // Vertex
        std::vector<MnvH1D*> vertex_z;
        
        // Other Event Parameters 
        TH1D* normal_rand_numbers;
        TH1D* em_shift_rand_numbers;
        TH1D* muonP_shift_rand_numbers;
        TH1D* muon_theta_shift_rand_numbers;
        TH1D* Birks_shift_rand_numbers;

        // Delta Resonance
        std::vector<MnvH1D*> deltaInvMass;
        std::vector<MnvH1D*> Polarization;
        std::vector<MnvH1D*> Delta_pi_theta;
        std::vector<MnvH1D*> Delta_pi_phi;
        MnvH2D* deltaInvMass_response;
        MnvH2D* Delta_pi_theta_response;
        MnvH2D* Delta_pi_phi_response;
        TH1D* resID;
        TH1D* resID_theta;

        MnvH1D* DeltaTransverse_data;
        MnvH1D* DeltaTransverse_mc;
        MnvH2D* DeltaTransverse_mc_res;
        
        TH1D* h_extra_muon_energy;
        TH1D* h_extra_leftover_energy;
        TH1D* h_extra_rejected_energy;
       
        // Short Proton
        MnvH1D* nProtons;
        TH1D* proton_true_P_1Track;
        TH1D* proton_true_KE_1Track;
        TH1D* proton_true_theta_1Track;
      
        // Ejected Nucleons
        TH1D* n_ejected_nucleons_1Track;
        TH1D* n_ejected_nucleons_2Track;
       
        // QSq Error, Difference
        TH2D* WSq_QSq_Diff;
        MnvH2D* QSq_All_response;
        MnvH2D* QSq_1Track_response;
        MnvH2D* QSq_2Track_response;

        TH1D* QSq_Error;
        TH1D* QSq_1Track_Error;
        TH1D* QSq_2Track_Error;

        TH1D* QSq_Diff;
        TH1D* QSq_1Track_Diff;
        TH1D* QSq_2Track_Diff;

        // Neutrino Energy: Error, Difference
        MnvH2D* Enu_All_response;
        MnvH2D* Enu_1Track_response;
        MnvH2D* Enu_2Track_response;

        TH1D* Enu_Error;
        TH1D* Enu_1Track_Error;
        TH1D* Enu_2Track_Error;
 
        TH1D* Enu_Diff;
        TH1D* Enu_1Track_Diff;
        TH1D* Enu_2Track_Diff;

        // W: Error, Difference
        TH1D* W_Error;
        TH1D* W_Diff;

        // Selected Signal 
        // Signal Q2
        TH1D* mc_Q2_QE;

        TH1D* mc_Q2_RES_1232;
        TH1D* mc_Q2_RES_1535;
        TH1D* mc_Q2_RES_1520;
        TH1D* mc_Q2_RES_Other;

        TH1D* mc_Q2_DIS;
        TH1D* mc_Q2_2p2h;
        TH1D* mc_Q2_Non_RES;

        // Signal Truth Q2
        TH1D* truth_QSq_QE;

        TH1D* truth_QSq_RES_1232;
        TH1D* truth_QSq_RES_1535;
        TH1D* truth_QSq_RES_1520;
        TH1D* truth_QSq_RES_Other;

        TH1D* truth_QSq_DIS;
        TH1D* truth_QSq_2p2h;
        TH1D* truth_QSq_Non_RES;

        // Background Truth Q2
        TH1D* reco_bckg_QSq_QE;

        TH1D* reco_bckg_QSq_RES_1232;
        TH1D* reco_bckg_QSq_RES_1535;
        TH1D* reco_bckg_QSq_RES_1520;
        TH1D* reco_bckg_QSq_RES_Other;

        TH1D* reco_bckg_QSq_DIS;
        TH1D* reco_bckg_QSq_2p2h;
        TH1D* reco_bckg_QSq_Non_RES;
        TH1D* reco_bckg_QSq_Coh;

        // Signal incomingE
        TH1D* mc_incomingE_QE;

        TH1D* mc_incomingE_RES_1232;
        TH1D* mc_incomingE_RES_1535;
        TH1D* mc_incomingE_RES_1520;
        TH1D* mc_incomingE_RES_Other;

        TH1D* mc_incomingE_DIS;
        TH1D* mc_incomingE_2p2h;
        TH1D* mc_incomingE_Non_RES;

        // Signal Truth incomingE
        TH1D* truth_Enu_QE;

        TH1D* truth_Enu_RES_1232;
        TH1D* truth_Enu_RES_1535;
        TH1D* truth_Enu_RES_1520;
        TH1D* truth_Enu_RES_Other;

        TH1D* truth_Enu_DIS;
        TH1D* truth_Enu_2p2h;
        TH1D* truth_Enu_Non_RES;

        // Background Truth incomingE
        TH1D* reco_bckg_Enu_QE;

        TH1D* reco_bckg_Enu_RES_1232;
        TH1D* reco_bckg_Enu_RES_1535;
        TH1D* reco_bckg_Enu_RES_1520;
        TH1D* reco_bckg_Enu_RES_Other;

        TH1D* reco_bckg_Enu_DIS;
        TH1D* reco_bckg_Enu_2p2h;
        TH1D* reco_bckg_Enu_Non_RES;
        TH1D* reco_bckg_Enu_Coh;

        // Signal w
        TH1D* mc_w_QE;

        TH1D* mc_w_RES_1232;
        TH1D* mc_w_RES_1535;
        TH1D* mc_w_RES_1520;
        TH1D* mc_w_RES_Other;

        TH1D* mc_w_DIS;
        TH1D* mc_w_2p2h;
        TH1D* mc_w_Non_RES;

        // Signal w
        TH1D* truth_w_QE;

        TH1D* truth_w_RES_1232;
        TH1D* truth_w_RES_1535;
        TH1D* truth_w_RES_1520;
        TH1D* truth_w_RES_Other;

        TH1D* truth_w_DIS;
        TH1D* truth_w_2p2h;
        TH1D* truth_w_Non_RES;

        // Background w
        TH1D* reco_bckg_w_QE;

        TH1D* reco_bckg_w_RES_1232;
        TH1D* reco_bckg_w_RES_1535;
        TH1D* reco_bckg_w_RES_1520;
        TH1D* reco_bckg_w_RES_Other;

        TH1D* reco_bckg_w_DIS;
        TH1D* reco_bckg_w_2p2h;
        TH1D* reco_bckg_w_Non_RES;
        TH1D* reco_bckg_w_Coh;

        // Signal reco w
        TH1D* reco_w_QE;

        TH1D* reco_w_RES_1232;
        TH1D* reco_w_RES_1535;
        TH1D* reco_w_RES_1520;
        TH1D* reco_w_RES_Other;

        TH1D* reco_w_DIS;
        TH1D* reco_w_2p2h;
        TH1D* reco_w_Non_RES;

    private:
        bool m_isMC;
        void initHistograms();
        
        TFile* f;
        std::string rootDir;
        
        CCProtonPi0_BinList binList;
};


#endif
