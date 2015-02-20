/*
================================================================================
Class: CAMAC_Data
    Reads Camac Output Data File
    Generates and Fill Histograms for the variables listed in Data File
    Updates the NearlineCurrentHistos.root
    NearlineCurrentHistos.root is used by GMBrowser to generate plots
    
    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision:  2015_02_20
================================================================================
*/
#ifndef CAMAC_Data_h
#define CAMAC_Data_h

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>

// ROOT Libraries
#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TString.h>

using namespace std;

class CAMAC_Data
{
    public:
        CAMAC_Data();
        ~CAMAC_Data();
        void ProcessData();
        void SetRootFileDir(string input);
        void SetDataFileDir(string input);
        void SetBranchName(string input);
        
    private:
        // Functions
        void OpenFiles();
            void ReadDataFile();
        void WriteRootFile();
        void ReadHeader(string header);
        void ReadScalar(string line, int varInd);
        void OpenDataFileError();
        void PrintVariables();
        void initHistograms();
        
        // I/O Files
        string rootDir;
        string dataDir;
        string branchName;
        
        ifstream dataFile;
        TFile* f_Root;
        
        // Variable Name Vector
        vector< string > vars;
        
        // Histogram Vector
        vector < TH1D* > hists;
        
        // Bins for the Histograms
        int TOF_nBins;
        double TOF_min;
        double TOF_max;

};

#endif
