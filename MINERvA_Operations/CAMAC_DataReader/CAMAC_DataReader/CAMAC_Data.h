/*
================================================================================
Class: CAMAC_Data
    Reads Camac Output Data File
    Generates and Fill Histograms for the variables listed in Data File
    Updates the NearlineCurrentHistos.root
    NearlineCurrentHistos.root is used by GMBrowser to generate plots
    
    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision:  2015_02_23
================================================================================
*/
#ifndef CAMAC_Data_h
#define CAMAC_Data_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>
#include <vector>

// ROOT Libraries
#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TString.h>
#include <TVectorD.h>
#include <TGraph.h>

// Local Libraries
#include "PlotFormats.h"

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
        void ReadInfoLine(string line, vector<string> &v);
        void ReadHistogramFormat(string line);
        void ReadGraphFormat(string line);
        void ReadFrequencyFormat(string line);
        void ReadDataLine(string line);
        void ProcessDataLine(string var_name, double value);
        void ProcessConvFactors();
        void OpenDataFileError();
        void PrintVector(vector<string> &v);
        void initHistograms();
        void initGraphVectors();
        void FillHistograms(string var_name, double value);
        void FillGraphVectors(string var_name, double value);
        int GetVariableInd(string var_name);
        
        
        bool isDebugging; 
       
        // I/O Files
        string rootDir;
        string dataDir;
        string branchName;
        
        ifstream dataFile;
        TFile* f_Root;
        
        // Vectors for Information Section
        vector<string> version;
        vector<string> vars;
        vector<string> convFactors_str;
        vector<double> convFactors;
        vector<string> units;
        
        // Vectors for Plot Formats
        vector<format_hist> f_hists;
        vector<format_graph> f_graphs;
        vector<format_frequency> f_frequencies;
        
        
        // Histogram Vector
        vector < TH1D* > hists;
        
        // Vectors for Graphs
        vector < TGraph* > graphs;
        vector < vector<double> > x_axis_vectors;
        vector < vector<double> > y_axis_vectors;
        
};

#endif
