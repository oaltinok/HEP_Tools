/*
================================================================================
Class: CAMAC_DataReader
    Main Class that controls CAMAC_DataReader
    
    Reader Objects:
        Config  - Reads Config File and initializes Plotters
        Data    - Reads Data File and Fills Plotters
    
    Plotter Objects:
        Hist    - 1D histograms
        Graph   - x vs y Graphs
        Freq    - Frequency Plots
    
    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision:  2015_03_10
================================================================================
*/
#ifndef CAMAC_DataReader_h
#define CAMAC_DataReader_h

#include <iostream>
#include <string>

// Include All Classes
#include "../Plotter/Plotter.h"
#include "../Hist/Hist.h"
#include "../Graph/Graph.h"
#include "../Freq/Freq.h"
#include "../Reader/Reader.h"
#include "../Data/Data.h"
#include "../Config/Config.h"

using namespace std;

/*
    Constants for the Package
        branchName      = Folder name inside ROOT File
        configDir       = Configuration file for CAMAC Readout
        dataDir_Auto    = Input Data File for Auto Mode
        rootDir_Auto    = Output Root File for Auto Mode
*/
const string branchName = "CAMAC_Data";
const string configDir = "camac_config.dat";
const string dataDir_Auto = "/home/nfs/minerva/daq/daqdata/lastspill_camac.dat";
const string rootDir_Auto = "/minerva/data/testbeam2/nearonline/CAMACDataHistos.root";

class CAMAC_DataReader{
    public:
        CAMAC_DataReader();
        CAMAC_DataReader(string mode);
        CAMAC_DataReader(string mode, string run_number, string subrun_number);
        ~CAMAC_DataReader();
        void RunManual();
        void RunAuto();

        // Plotters
        Hist HistPlotter;
        Graph GraphPlotter;
        Freq FreqPlotter;
        
        TFile* f_Root;

        // Configuration Vectors
        vector<string> version;
        vector<string> vectrun;
        vector<string> vectsubrun;
        vector<string> vars;
        vector<string> convFactors_str;
        vector<double> convFactors; // convFactor read as string then converted to double  
        vector<string> units;
        
        
    private:
        // Init Functions
        void InitParamaters();
        void InitReaders();
        void InitPlotters();
        void ResetPlots();
        void InitPlots();
        void SetConfigDir();
        void SetConfigDir_Manual();
        void SetDataDir();
        void SetRootDir();
        void OpenRootFile();
        string Get_dataDir_Manual();
        string Get_rootDir_Manual();
        
        // Flow Functions
        void CheckDataFile();
        void CheckRunSubrun();
        
        // Run Functions
	void ReadFirstTime();
	void ReadRunSubrun();
        void ReadConfigFile();
        void ReadDataFile();
        void WriteRootFile();
        
        // Other 
        void ErrorMode();
	void ResetVectors();
        
        // Readers
        Config ConfigReader;
        Data DataReader;
        
        bool isModeAuto;
        bool newFile;
        bool resetPlots;

        string mode;
        
        string rootDir;
        string dataDir;
        ifstream dataFile;
        
        string run;
        string subrun;
        string time;

        string latest_run;
        string latest_subrun;
        string latest_time;
        
};

#endif

