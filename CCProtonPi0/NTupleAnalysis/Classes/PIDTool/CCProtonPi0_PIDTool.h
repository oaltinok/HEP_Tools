/*
================================================================================
Class: CCProtonPi0_PIDTool
    CCProtonPi0_PIDTool class responsible for Particle Identification Optimization
    Generates pID specific Histograms and Statistics
    
    Author:        Ozgur Altinok  - ozgur.altinok@tufts.edu
================================================================================
*/
#ifndef CCProtonPi0_PIDTool_h
#define CCProtonPi0_PIDTool_h

#include <iostream>
#include <cstdlib>

// ROOT Libraries
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

// Classes
#include "../NTupleAnalysis/CCProtonPi0_NTupleAnalysis.h"
#include "../BinList/CCProtonPi0_BinList.h"

// Libraries
#include "../../Libraries/PDG_List.h"

class CCProtonPi0_PIDTool : public CCProtonPi0_NTupleAnalysis
{
    public:
        // Functions
        CCProtonPi0_PIDTool(int nMode, bool isMC);
        void FillHistograms(double protonScore_LLR, double protonScore, double pionScore,
                            int truthPDG, double prongE);
        void get_pID_Stats();
        void get_pID_Stats_LLR();
        void get_pID_Stats_pIDDiff();
        void write_RootFile();
        
        // Histograms
        TH1D* purity_LLR;
        TH1D* efficiency_LLR;
        TH1D* purityXefficiency_LLR;

        TH1D* purity_pIDDiff;
        TH1D* efficiency_pIDDiff;
        TH1D* purityXefficiency_pIDDiff;
        
        TH1D* proton_protonScore_LLR;
        TH1D* proton_protonScore;
        TH2D* proton_pionScore_protonScore;
        TH1D* proton_pIDDiff;
        TH2D* proton_protonScore_protonScore_LLR;
        
        TH1D* piplus_protonScore_LLR;
        TH1D* piplus_protonScore;
        TH2D* piplus_pionScore_protonScore;
        TH1D* piplus_pIDDiff;
        TH2D* piplus_protonScore_protonScore_LLR;
        
        TH1D* piminus_protonScore_LLR;
        TH1D* piminus_protonScore;
        TH2D* piminus_pionScore_protonScore;
        TH1D* piminus_pIDDiff;
        TH2D* piminus_protonScore_protonScore_LLR;
        
        TH1D* other_protonScore_LLR;
        TH1D* other_pIDDiff;
        
        TH1D* KE_proton_pIDDiff;
        TH1D* KE_other_pIDDiff;
        TH1D* KE_proton_LLR;
        TH1D* KE_other_LLR;
 
    private:
        TFile* f;
        std::string rootDir;
        
        CCProtonPi0_BinList binList;

        void initHistograms();
};


#endif

