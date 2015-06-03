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
        CCProtonPi0_Pi0Blob(int nMode);

        void write_RootFile();

        // Histograms
        TH1D* g1_blob_ndigits;
        TH1D* g1_blob_nclusters;
        TH1D* g1_blob_energy;
        TH1D* g1_blob_minsep;

        TH1D* g2_blob_ndigits;
        TH1D* g2_blob_nclusters;
        TH1D* g2_blob_energy;
        TH1D* g2_blob_minsep;

    private:
        void initHistograms();

        TFile* f;
        string rootDir;
        
        // Histogram Binning
        CCProtonPi0_SingleBin bin_blob_nclusters;
        CCProtonPi0_SingleBin bin_blob_energy;
        CCProtonPi0_SingleBin bin_blob_ndigits;
        CCProtonPi0_SingleBin bin_blob_minsep;

};

#endif

