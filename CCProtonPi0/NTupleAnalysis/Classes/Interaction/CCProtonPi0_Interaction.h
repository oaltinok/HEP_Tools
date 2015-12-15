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
        CCProtonPi0_Interaction(bool isModeReduce, bool isMC, std::string ana_folder);
        void writeHistograms();
        
        //--------------------------------------------------------------------------
        //     Histograms
        //--------------------------------------------------------------------------
        // Event Kinematics
        std::vector<MnvH1D*> Enu_True;
        std::vector<MnvH1D*> Enu_1Track;
        std::vector<MnvH1D*> Enu_1Track_Alt;
        std::vector<MnvH1D*> Enu_2Track;
        std::vector<MnvH1D*> QSq;
        std::vector<MnvH1D*> WSq;
        std::vector<MnvH1D*> W;
  
        // Vertex & Extra Energy
        std::vector<MnvH1D*> vertex_energy_1Track;
        std::vector<MnvH1D*> vertex_energy_2Track;
        std::vector<MnvH1D*> vertex_energy_Corrected_1Track;
        std::vector<MnvH1D*> vertex_energy_Corrected_2Track;
        std::vector<MnvH1D*> vertex_evis_1Track;
        std::vector<MnvH1D*> vertex_evis_2Track;
        std::vector<MnvH1D*> extra_energy_1Track;
        std::vector<MnvH1D*> extra_energy_2Track;


        // Truth Information
        TH1D* Enu_1Track_Error;
        TH1D* Enu_1Track_Alt_Error;
        TH1D* Enu_2Track_Error;
        TH1D* Enu_1Track_Corrected_Error;
        TH1D* Enu_2Track_Corrected_Error;
       
        TH2D* Enu_1Track_True;
        TH2D* Enu_1Track_Alt_True;
        TH2D* Enu_1Track_1Track_Alt;
        TH2D* Enu_2Track_True;

        // Short Protons
        TH1D* proton_true_P_1Track;
        TH1D* proton_true_KE_1Track;

        // Vertex
        std::vector<MnvH1D*> vertex_z;

        // Reconstruction 
        std::vector<MnvH1D*> E_Unused_afterReco;
        std::vector<MnvH1D*> E_Used_afterReco;
        
        // Other Event Parameters 
        std::vector<MnvH1D*> deltaInvMass;
       
        // MC Only Histograms
        TH1D* final_mc_w_DIS;
        TH1D* final_mc_w_RES;
        TH1D* final_mc_w_CCQE;

        // --------------------------------------------------------------------
        // Neutrino Energy Study
        // --------------------------------------------------------------------
        TH1D* extra_energy_50_1Track;
        TH1D* extra_energy_100_1Track;
        TH1D* extra_energy_150_1Track;
        TH1D* extra_energy_200_1Track;
        TH1D* extra_energy_300_1Track;
        TH1D* extra_energy_500_1Track;

        TH1D* extra_energy_50_2Track;
        TH1D* extra_energy_100_2Track;
        TH1D* extra_energy_150_2Track;
        TH1D* extra_energy_200_2Track;
        TH1D* extra_energy_300_2Track;
        TH1D* extra_energy_500_2Track;

        TH2D* vertex_evis_true_proton_KE;
        TH2D* vertex_evis_vertex_evis_ratio;
        

        // 1 Track
        std::vector<MnvH1D*> evis_total_1Track;
        std::vector<MnvH1D*> evis_muon_1Track;
        std::vector<MnvH1D*> evis_pi0_1Track;
      
        std::vector<MnvH1D*> evis_hadron_1Track;
        std::vector<MnvH1D*> energy_hadron_true_1Track;
        std::vector<MnvH1D*> energy_hadron_reco_1Track;
        TH2D* energy_hadron_reco_true_1Track;
        
        std::vector<MnvH1D*> evis_extra_1Track;
        std::vector<MnvH1D*> energy_extra_true_1Track;
        std::vector<MnvH1D*> energy_extra_reco_1Track;

        TH2D* evis_extra_energy_extra_true_1Track;
        
        // 2 Track
        std::vector<MnvH1D*> evis_total_2Track;
        std::vector<MnvH1D*> evis_muon_2Track;
        std::vector<MnvH1D*> evis_pi0_2Track;
        std::vector<MnvH1D*> evis_proton_2Track;

        std::vector<MnvH1D*> evis_hadron_2Track;
        std::vector<MnvH1D*> energy_hadron_true_2Track;
        std::vector<MnvH1D*> energy_hadron_reco_2Track;
        TH2D* energy_hadron_reco_true_2Track;

        std::vector<MnvH1D*> evis_hadron_nopi0_2Track;
        std::vector<MnvH1D*> energy_hadron_nopi0_true_2Track;
        std::vector<MnvH1D*> energy_hadron_nopi0_reco_2Track;
        
        std::vector<MnvH1D*> evis_extra_2Track;
        std::vector<MnvH1D*> energy_extra_true_2Track;
        TH2D* evis_extra_energy_extra_true_2Track;
        // --------------------------------------------------------------------
 
    private:
        void initHistograms();
        
        TFile* f;
        string rootDir;
        
        CCProtonPi0_BinList binList;
};


#endif
