/*
================================================================================
Class: MichelTool
    Micheltool class responsible for Michel Tool Optimization
    Generates Michel Tool specific Histograms
    
    Author:        Ozgur Altinok  - ozgur.altinok@tufts.edu
================================================================================
*/

#ifndef CCProtonPi0_MichelTool_h
#define CCProtonPi0_MichelTool_h

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


class CCProtonPi0_MichelTool : public CCProtonPi0_NTupleAnalysis
{
    public: 
        CCProtonPi0_MichelTool(int nMode, bool isMC);
        void write_RootFile();

        // Histograms 
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
        
        TH1D* trueMichel_dist_vtx;
        TH1D* fakeMichel_dist_vtx;
        TH1D* trueMichel_dist_end_point;
        TH1D* fakeMichel_dist_end_point;
        TH1D* trueMichel_end_Z;
        TH1D* fakeMichel_end_Z;
        TH1D* trueMichel_end_Z_vtx_Z;
        TH1D* fakeMichel_end_Z_vtx_Z;
        TH1D* trueMichel_energy;
        TH1D* fakeMichel_energy;
        TH1D* trueMichel_time_diff;
        TH1D* fakeMichel_time_diff;
        

    private:
        void initHistograms();
        
        TFile* f;
        std::string rootDir;
        CCProtonPi0_BinList binList;

        double N_trueMichel_before;
        double N_trueMichel_after;
        double N_trueMichel_afterAll;
        double N_noMichel_before;
        double N_noMichel_after;
        double N_detectedMichel_true;
        double N_detectedMichel_fake;
        double N_missedMichel_true;
        double N_missedMichel_fake;
        double N_detectedMichel_true_signal;
        double N_detectedMichel_fake_signal;
        double N_missedMichel_true_signal;
        double N_missedMichel_fake_signal;
        double N_selected_detectedMichel_true;
        double N_selected_detectedMichel_fake;
        double N_selected_missedMichel_true;
        double N_selected_missedMichel_fake;


};



#endif
