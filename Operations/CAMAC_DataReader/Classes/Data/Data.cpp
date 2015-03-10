#ifndef Data_cpp
#define Data_cpp

#include "Data.h"

using namespace std;


Data::Data()
{
    isDebugging = false;
    
    newFile = true;
    resetPlots = true;
    

}

void Data::ReadFile()
{
    bool isLine_RunSubrun = true;
    
    string line;

    while (!dataFile.eof())
    {
        getline (dataFile,line);
        
        // Skip Empty Lines
        if (line.size() == 0) continue;
        
        // Break if you read SENTINEL for End of File
        if (line.compare(SENTINEL_EOF) == 0) break;
        
        // First Line is Run Subrun TimeStamp
        if (isLine_RunSubrun){
            ReadRunSubrun(line);
            isLine_RunSubrun = false;
            cout<<"Processing "<<run<<"/"<<subrun<<" "<<time<<endl;
            continue;
        }
        
        ReadDataLine(line);
    }
    
    // Close File after finishing
    dataFile.close();

}



void Data::ReadRunSubrun(string line)
{
    stringstream line_stream(line);
    
    line_stream >> run >> subrun >> time;
}
 
/*
    Reads a single Data line and processes data  
*/
void Data::ReadDataLine(string line)
{
    if(isDebugging){
        cout<<"----"<<endl;
        cout<<"Reading Data:"<<endl; 
        cout<<line<<"\n\n";
    }
    
    stringstream line_stream(line);
    string var_name;
    double value;
    
    line_stream >> var_name >> value;
    
    // Skip Summary Line
    if(var_name.compare("Summary") == 0) return;

    ProcessDataLine(var_name,value);    
}

/*
    Process given data
        Inputs:
            variable name
            value
*/
void Data::ProcessDataLine(string var_name, double value)
{
    // Get Variable Indice
    int var_ind = GetVariableInd(var_name);
    
    // Get Variable Conversion Factor
    double var_convFactor =  (*convFactors)[var_ind]; 
    
    // Convert Scalar Value to REAL Value using Conversion Factor
    value = value * var_convFactor;
    
    // Fill Histograms - Search Histograms for the Variable and fill if necessary
    hist->Fill(var_name,value);
    
    // Fill Graphs - Search Graphs for the Variable and fill if necessary
    graph->Fill(var_name,value);
    
    // Fill Frequency Plots - Check Frequency Arrays to count the Variable
    freq->Fill(var_ind);
}

int Data::GetVariableInd(string var_name)
{
    if ( vars->size() == 0 ){
        cout<<"ERROR: No Variables, Did you read config file correctly!"<<endl;
        exit(EXIT_FAILURE);
    }
    
    int ind = -1;
    
    for(unsigned int i = 0; i < vars->size(); i++){
        if ( var_name.compare((*vars)[i]) == 0){ 
            ind = i;
            break;
        }
    }
    
    if (ind != -1){ 
        return ind;
    }else{
        cout<<"ERROR: No Variable = "<<var_name<<endl;
        exit(EXIT_FAILURE);
    }
}







#endif
