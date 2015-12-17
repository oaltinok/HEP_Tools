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
        std::vector<MnvH1D*> vertex_evis_1Track;
        std::vector<MnvH1D*> vertex_evis_2Track;
        std::vector<MnvH1D*> extra_energy_1Track;
        std::vector<MnvH1D*> extra_energy_2Track;
       
        // Vertex
        std::vector<MnvH1D*> vertex_z;
        
        // Other Event Parameters 
        std::vector<MnvH1D*> deltaInvMass;
       
        // MC Only Histograms
        TH1D* final_mc_w_DIS;
        TH1D* final_mc_w_RES;
        TH1D* final_mc_w_CCQE;
        
        TH1D* n_ejected_nucleons_1Track;
        TH1D* n_ejected_nucleons_2Track;
       
        // Truth Information
        TH1D* Enu_1Track_Error;
        TH1D* Enu_1Track_Alt_Error;
        TH1D* Enu_2Track_Error;
       
        TH2D* Enu_1Track_True;
        TH2D* Enu_1Track_Alt_True;
        TH2D* Enu_1Track_1Track_Alt;
        TH2D* Enu_2Track_True;
      
        TH1D* proton_true_P_1Track;
        TH1D* proton_true_KE_1Track;

    private:
        void initHistograms();
        
        TFile* f;
        string rootDir;
        
        CCProtonPi0_BinList binList;
};


#endif
