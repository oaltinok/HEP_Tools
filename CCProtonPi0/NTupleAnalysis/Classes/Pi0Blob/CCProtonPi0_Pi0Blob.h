/*
================================================================================
Class: CCProtonPi0_Pi0Blob
    CCProtonPi0_Pi0Blob Class Contains Blob Histograms
    All Histograms declared public and can be accessed by Analyzer
                
    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
================================================================================
*/
#ifndef CCProtonPi0_Pi0Blob_h
#define CCProtonPi0_Pi0Blob_h

#include <iostream>
#include <string>
#include <cstdlib>

#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>

#include "../NTupleAnalysis/CCProtonPi0_NTupleAnalysis.h"
#include "../SingleBin/CCProtonPi0_SingleBin.h"

using namespace std;

class CCProtonPi0_Pi0Blob : public CCProtonPi0_NTupleAnalysis
{
    public:
        CCProtonPi0_Pi0Blob(int nMode, bool isMC);

        void write_RootFile();

        // Histograms
        
        // Standard Blob Variables
        TH1D* g1_blob_ndigits;
        TH1D* g1_blob_nclusters;
        TH1D* g1_blob_energy;
        TH1D* g1_blob_minsep;
        TH1D* g2_blob_ndigits;
        TH1D* g2_blob_nclusters;
        TH1D* g2_blob_energy;
        TH1D* g2_blob_minsep;


        // Fit Results
        TH1D* g1_blob_ndof;
        TH1D* g1_blob_fval;

        // dEdX 
        TH1D* g1_blob_dEdx_doublet;
        TH1D* g1_blob_dEdx_empty_plane;
        TH1D* g1_blob_dEdx;
        TH1D* g1_blob_dEdx1;
        TH1D* g1_blob_dEdx_nplane;
        TH1D* g1_blob_dEdx_cluster_energy;

        private:
        void initHistograms();

        TFile* f;
        string rootDir;
        
        // Histogram Binning
        CCProtonPi0_SingleBin bin_blob_nclusters;
        CCProtonPi0_SingleBin bin_blob_energy;
        CCProtonPi0_SingleBin bin_blob_ndigits;
        CCProtonPi0_SingleBin bin_blob_minsep;
        CCProtonPi0_SingleBin bin_blob_ndof;
        CCProtonPi0_SingleBin bin_blob_fval;
        CCProtonPi0_SingleBin bin_blob_dEdx_doublet;
        CCProtonPi0_SingleBin bin_blob_dEdx;
        CCProtonPi0_SingleBin bin_blob_dEdx_nplane;
        CCProtonPi0_SingleBin bin_blob_dEdx_cluster_energy;

};

#endif

