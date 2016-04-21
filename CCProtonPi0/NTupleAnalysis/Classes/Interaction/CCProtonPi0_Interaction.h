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
        // Event Kinematics
        std::vector<MnvH1D*> Enu_1Track;
        std::vector<MnvH1D*> Enu_1Track_Alt;
        std::vector<MnvH1D*> Enu_2Track;
        std::vector<MnvH1D*> Enu;
        std::vector<MnvH1D*> QSq;
        std::vector<MnvH1D*> WSq;
        std::vector<MnvH1D*> W;
        std::vector<MnvH1D*> W_Calc;
  
        // Vertex & Extra Energy
        std::vector<MnvH1D*> vertex_energy_1Track;
        std::vector<MnvH1D*> vertex_evis_1Track;
        std::vector<MnvH1D*> extra_evis_1Track;
        std::vector<MnvH1D*> extra_muon_energy_1Track;
        std::vector<MnvH1D*> extra_dispersed_energy_1Track;
        std::vector<MnvH1D*> extra_rejected_energy_1Track;
        std::vector<MnvH1D*> extra_total_energy_1Track;
        
        std::vector<MnvH1D*> vertex_energy_2Track;
        std::vector<MnvH1D*> vertex_evis_2Track;
        std::vector<MnvH1D*> extra_evis_2Track;
        std::vector<MnvH1D*> extra_muon_energy_2Track;
        std::vector<MnvH1D*> extra_dispersed_energy_2Track;
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

        // Vertex
        std::vector<MnvH1D*> vertex_z;
        
        // Other Event Parameters 
        std::vector<MnvH1D*> deltaInvMass;

        TH1D* h_extra_muon_energy;
        TH1D* h_extra_dispersed_energy;
        TH1D* h_extra_rejected_energy;
       
        // MC Only Histograms
        TH1D* final_mc_w_DIS;
        TH1D* final_mc_w_RES;

        TH1D* final_mc_Q2_DIS;
        TH1D* final_mc_Q2_RES;

        // Short Proton
        TH1D* proton_true_P_1Track;
        TH1D* proton_true_KE_1Track;
      
        // Ejected Nucleons
        TH1D* n_ejected_nucleons_1Track;
        TH1D* n_ejected_nucleons_2Track;
       
        // QSq Truth, Error, Difference
        TH1D* QSq_True;
        TH1D* QSq_Error;
        TH1D* QSq_Diff;

        // Neutrino Energy: Truth, Error, Difference
        TH1D* Enu_True_1Track;
        TH1D* Enu_True_2Track;
        
        TH1D* Enu_1Track_Error;
        TH1D* Enu_1Track_Alt_Error;
        TH1D* Enu_2Track_Error;
 
        TH1D* Enu_1Track_Diff;
        TH1D* Enu_2Track_Diff;

    private:
        void initHistograms();
        
        TFile* f;
        std::string rootDir;
        
        CCProtonPi0_BinList binList;
};


#endif
