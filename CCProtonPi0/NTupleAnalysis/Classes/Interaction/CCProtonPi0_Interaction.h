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
        
        //--------------------------------------------------------------------------
        //     Histograms
        //--------------------------------------------------------------------------
        std::vector<MnvH1D*> CV_weight;
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
  
        // Vertex & Extra Energy
        std::vector<MnvH1D*> vertex_energy_1Track;
        std::vector<MnvH1D*> vertex_evis_1Track;
        std::vector<MnvH1D*> extra_leftover_energy_1Track;
        std::vector<MnvH1D*> extra_muon_energy_1Track;
        std::vector<MnvH1D*> extra_rejected_energy_1Track;
        std::vector<MnvH1D*> extra_total_energy_1Track;
        
        std::vector<MnvH1D*> vertex_energy_2Track;
        std::vector<MnvH1D*> vertex_evis_2Track;
        std::vector<MnvH1D*> extra_leftover_energy_2Track;
        std::vector<MnvH1D*> extra_muon_energy_2Track;
        std::vector<MnvH1D*> extra_rejected_energy_2Track;
        std::vector<MnvH1D*> extra_total_energy_2Track;

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

        // Vertex
        std::vector<MnvH1D*> vertex_z;
        
        // Other Event Parameters 
        std::vector<MnvH1D*> deltaInvMass;
        TH1D* normal_rand_numbers;
        TH1D* em_shift_rand_numbers;
        TH1D* muonP_shift_rand_numbers;

        // Delta Resonance
        MnvH1D* Polarization_data;
        MnvH1D* Polarization_mc;
        
        MnvH1D* DeltaTransverse_data;
        MnvH1D* DeltaTransverse_mc;
        MnvH2D* DeltaTransverse_mc_res;
        
        TH1D* h_extra_muon_energy;
        TH1D* h_extra_leftover_energy;
        TH1D* h_extra_rejected_energy;
       
        // Short Proton
        TH1D* proton_true_P_1Track;
        TH1D* proton_true_KE_1Track;
        TH1D* proton_true_theta_1Track;
      
        // Ejected Nucleons
        TH1D* n_ejected_nucleons_1Track;
        TH1D* n_ejected_nucleons_2Track;
       
        // QSq Error, Difference
        TH2D* WSq_QSq_Diff;
        MnvH2D* QSq_1Track_response;
        MnvH2D* QSq_2Track_response;

        TH1D* QSq_Error;
        TH1D* QSq_1Track_Error;
        TH1D* QSq_2Track_Error;

        TH1D* QSq_Diff;
        TH1D* QSq_1Track_Diff;
        TH1D* QSq_2Track_Diff;

        // Neutrino Energy: Error, Difference
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

        TH1D* mc_Q2_DIS_1_pi;
        TH1D* mc_Q2_DIS_2_pi;
        TH1D* mc_Q2_DIS_Multi_pi;
        TH1D* mc_Q2_Non_RES;

        // Signal Truth Q2
        TH1D* truth_QSq_QE;

        TH1D* truth_QSq_RES_1232;
        TH1D* truth_QSq_RES_1535;
        TH1D* truth_QSq_RES_1520;
        TH1D* truth_QSq_RES_Other;

        TH1D* truth_QSq_DIS_1_pi;
        TH1D* truth_QSq_DIS_2_pi;
        TH1D* truth_QSq_DIS_Multi_pi;
        TH1D* truth_QSq_Non_RES;

        // Signal incomingE
        TH1D* mc_incomingE_QE;

        TH1D* mc_incomingE_RES_1232;
        TH1D* mc_incomingE_RES_1535;
        TH1D* mc_incomingE_RES_1520;
        TH1D* mc_incomingE_RES_Other;

        TH1D* mc_incomingE_DIS_1_pi;
        TH1D* mc_incomingE_DIS_2_pi;
        TH1D* mc_incomingE_DIS_Multi_pi;
        TH1D* mc_incomingE_Non_RES;

        // Signal Truth incomingE
        TH1D* truth_Enu_QE;

        TH1D* truth_Enu_RES_1232;
        TH1D* truth_Enu_RES_1535;
        TH1D* truth_Enu_RES_1520;
        TH1D* truth_Enu_RES_Other;

        TH1D* truth_Enu_DIS_1_pi;
        TH1D* truth_Enu_DIS_2_pi;
        TH1D* truth_Enu_DIS_Multi_pi;
        TH1D* truth_Enu_Non_RES;

        // Signal w
        TH1D* mc_w_QE;

        TH1D* mc_w_RES_1232;
        TH1D* mc_w_RES_1535;
        TH1D* mc_w_RES_1520;
        TH1D* mc_w_RES_Other;

        TH1D* mc_w_DIS_1_pi;
        TH1D* mc_w_DIS_2_pi;
        TH1D* mc_w_DIS_Multi_pi;
        TH1D* mc_w_Non_RES;

        // Signal w
        TH1D* truth_w_QE;

        TH1D* truth_w_RES_1232;
        TH1D* truth_w_RES_1535;
        TH1D* truth_w_RES_1520;
        TH1D* truth_w_RES_Other;

        TH1D* truth_w_DIS_1_pi;
        TH1D* truth_w_DIS_2_pi;
        TH1D* truth_w_DIS_Multi_pi;
        TH1D* truth_w_Non_RES;

        // Signal reco w
        TH1D* reco_w_QE;

        TH1D* reco_w_RES_1232;
        TH1D* reco_w_RES_1535;
        TH1D* reco_w_RES_1520;
        TH1D* reco_w_RES_Other;

        TH1D* reco_w_DIS_1_pi;
        TH1D* reco_w_DIS_2_pi;
        TH1D* reco_w_DIS_Multi_pi;
        TH1D* reco_w_Non_RES;

    private:
        void initHistograms();
        
        TFile* f;
        std::string rootDir;
        
        CCProtonPi0_BinList binList;
};


#endif
