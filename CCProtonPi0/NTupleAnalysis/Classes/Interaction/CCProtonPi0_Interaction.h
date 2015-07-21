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
#include <TFile.h>
#include <PlotUtils/MnvH1D.h>
#include <PlotUtils/MnvH2D.h>

#include "../NTupleAnalysis/CCProtonPi0_NTupleAnalysis.h"
#include "../BinList/CCProtonPi0_BinList.h"

using namespace PlotUtils;

class CCProtonPi0_Interaction : public CCProtonPi0_NTupleAnalysis
{
    public:
        CCProtonPi0_Interaction(int nMode, bool isMC);
        void writeHistograms();
        
        //--------------------------------------------------------------------------
        //     Histograms
        //--------------------------------------------------------------------------
        // Event Kinematics
        MnvH1D* Enu_1Track;
        MnvH1D* Enu_2Track;
        MnvH1D* Enu_Cal;
        MnvH1D* q2;
        MnvH1D* w;
        MnvH1D* wSq;
       
        // Reconstruction 
        MnvH1D* E_Unused_afterReco;
        MnvH1D* E_Used_afterReco;
        
        // Other Event Parameters 
        MnvH1D* deltaInvMass;
        MnvH1D* nProngs_hist;
       
        // MC Only Histograms
        TH1D* mc_w_DIS;
        TH1D* mc_w_RES;
        TH1D* mc_w_CCQE;
        
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
