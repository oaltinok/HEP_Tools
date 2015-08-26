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
        std::vector<MnvH1D*> Enu_1Track;
        std::vector<MnvH1D*> Enu_2Track;
        std::vector<MnvH1D*> Enu_Cal;
        std::vector<MnvH1D*> q2;
        std::vector<MnvH1D*> w;
        std::vector<MnvH1D*> wSq;

        // Vertex
        std::vector<MnvH1D*> vertex_z;

        // Reconstruction 
        std::vector<MnvH1D*> E_Unused_afterReco;
        std::vector<MnvH1D*> E_Used_afterReco;
        
        // Other Event Parameters 
        std::vector<MnvH1D*> deltaInvMass;
        std::vector<MnvH1D*> nProngs_hist;
       
        // MC Only Histograms
        TH1D* final_mc_w_DIS;
        TH1D* final_mc_w_RES;
        TH1D* final_mc_w_CCQE;
       
    private:
        void initHistograms();
        
        TFile* f;
        string rootDir;
        
        CCProtonPi0_BinList binList;
};


#endif
