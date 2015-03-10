#ifndef Reader_h
#define Reader_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <vector>
#include <TROOT.h>
#include <TString.h>

// Include Plot Objects
#include "../Plotter/Plotter.h"
#include "../Hist/Hist.h"
#include "../Graph/Graph.h"
#include "../Freq/Freq.h"

using namespace std;

const string SENTINEL_EOF = "CAMAC_Readout_End";

class Reader{
    
    public:
        bool isDebugging;
        string dataDir;
        ifstream dataFile;

        ~Reader();
        bool IsFileComplete();
        bool OpenFile(string input);
        void OpenDataFileError(); 
        
        void SetVectors(    vector<string>* p_version,
                            vector<string>* p_vars,
                            vector<string>* p_convFactors_str,
                            vector<double>* p_convFactors, // convFactor read as string then converted to double  
                            vector<string>* p_units);
        
        void SetPlots( Hist* p_hist, Graph* p_graph, Freq* p_freq );
        
        
        // Pointers to Configuration Vectors
        vector<string>* version;
        vector<string>* vars;
        vector<string>* convFactors_str;
        vector<double>* convFactors; // convFactor read as string then converted to double  
        vector<string>* units;
        
                
        // Pointers to Plots
        Hist* hist;
        Graph* graph;
        Freq* freq;
        
    private:


        


};


#endif
