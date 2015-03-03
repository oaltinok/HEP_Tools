/*
================================================================================
Class: CAMAC_Data
    Reads Camac Output Data File
    Generates and Fill Histograms for the variables listed in Data File
    Updates the NearlineCurrentHistos.root
    NearlineCurrentHistos.root is used by GMBrowser to generate plots
    
    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision:  2015_03_03
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
#include <TCanvas.h>

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
        bool ChangeDataFileDir(string input);
        string GetRunNumber();
        string GetSubrunNumber();
        string GetTimeStamp();
        string GetFileName();
        
    private:
        // ---------------------------------------------------------------------
        // Functions
        // ---------------------------------------------------------------------
        // File Operations
        void OpenFiles();
        void ReadDataFile();
        void WriteRootFile();
        void OpenDataFileError();
        bool IsNewFile();
        bool IsNewSubrun();
        
        // Reading File
        void ReadInfoLine(string line, vector<string> &v);
        void ReadHistogramFormat(string line);
        void ReadGraphFormat(string line);
        void ReadFrequencyFormat(string line);
        void ReadDataLine(string line);
        
        // Initialize Data Structures
        void initHistograms();
        void initGraphVectors();
        void initFreqVectors();
        void ProcessFileName(string input);
        
        // Process Read Data
        int GetVariableInd(string var_name);
        void ProcessDataLine(string var_name, double value);
        void ProcessConvFactors();
        void ClearVectors();
        void ClearSubrunVectors();
        
        // Fill Plots using Processed Data
        void FillHistograms(string var_name, double value);
        void FillGraphVectors(string var_name, double value);
        void FillFrequencyPlots(int ind);
        void FillSubrunFrequencyPlotVectors();
        void FillSubrunGraphVectors();
     
        // Debugging
        void PrintVector(vector<string> &v);
        
        // ---------------------------------------------------------------------
        // Variables
        // ---------------------------------------------------------------------
        bool isDebugging; // Make true for debugging
        bool isNewSubrun;
       
        string runNumber;
        string subrunNumber;
        string timeStamp;
        
        string latest_runNumber;
        string latest_subrunNumber;
        string latest_timeStamp;
        
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
        vector<double> convFactors; // convFactor read as string then converted to double  
        vector<string> units;
        
        // Vectors for Plot Formats
        vector<format_hist> f_hists;
        vector<format_graph> f_graphs;
        vector<format_frequency> f_frequencies;
        
        // Histogram Vector
        vector < TH1D* > hists;
        
        // Vectors for Graphs
        vector < TGraph* > graphs;
        vector < vector<double> > graph_x_axis_vectors;
        vector < vector<double> > graph_y_axis_vectors;
        vector < vector<double> > subrun_graph_x_axis_vectors;
        vector < vector<double> > subrun_graph_y_axis_vectors;
        
        // Vectors for Frequency Plots
        vector < TGraph* > freq_plots;
        vector < vector<double> > freq_x_axis_vectors;
        vector < vector<double> > freq_y_axis_vectors;
        vector < vector<double> > subrun_freq_x_axis_vectors;
        vector < vector<double> > subrun_freq_y_axis_vectors;
};

#endif
