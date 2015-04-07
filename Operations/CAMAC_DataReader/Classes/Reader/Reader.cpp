#ifndef Reader_cpp
#define Reader_cpp

#include "Reader.h"

using namespace std;

bool Reader::IsFileComplete()
{
    string line;
    
    while (!dataFile.eof() ){
        getline (dataFile,line);
        if(line.compare(SENTINEL_EOF) == 0){
            return true;   
        }
    }
     
    return false;
}

bool Reader::OpenFile(string input)
{
    dataDir = input;
    // Open Data File
    dataFile.open(dataDir.c_str());
    if (!dataFile.is_open()) OpenDataFileError();
    
    // Check File is Completely Written by CAMAC Readout
    if(!IsFileComplete()){
        cout<<"File is NOT ready for CAMAC_DataReader"<<endl;
        dataFile.close();
        return false;
    }else{
        // Go to the beginning of the FILE
        dataFile.clear();
        dataFile.seekg (0, ios::beg);
        return true;
    }
}

void Reader::SetVectors(    vector<string>* p_version,
                            vector<string>* p_run,
                            vector<string>* p_subrun,
                            vector<string>* p_vars,
                            vector<string>* p_convFactors_str,
                            vector<double>* p_convFactors, // convFactor read as string then converted to double  
                            vector<string>* p_units)
{
        version = p_version;
	run = p_run;
	subrun = p_subrun;
        vars = p_vars;
        convFactors_str = p_convFactors_str;
        convFactors = p_convFactors; // convFactor read as string then converted to double  
        units = p_units;
}


void Reader::SetPlots( Hist* p_hist, Graph* p_graph, Freq* p_freq )
{
    hist = p_hist;
    graph = p_graph;
    freq = p_freq;
}

Reader::~Reader()
{
    // Close Data File
    dataFile.close();
}
/*
    OpenDataFileError()
        Called if the Input Data File can not be opened
*/
void Reader::OpenDataFileError()
{
    cerr<<"Cannot Open File = "<<dataDir<<endl;
    exit(EXIT_FAILURE);
}




#endif
