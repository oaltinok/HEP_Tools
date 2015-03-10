#ifndef Config_cpp
#define Config_cpp

#include "Config.h"

using namespace std;

Config::Config()
{
    isDebugging = false;   
}

void Config::ReadFile()
{
    bool isLine_Version = true;
    bool isLine_VariableNames = true;
    bool isLine_ConvFactors = true;
    bool isLine_Units = true;
    
    string line;

    while (!dataFile.eof())
    {
        getline (dataFile,line);
        
        // Skip Empty Lines
        if (line.size() == 0) continue;
        
        // Break if you read SENTINEL for End of File
        if (line.compare(SENTINEL_EOF) == 0) break;
        
        if (isLine_Version){
            ReadConfigLine(line,version);
            isLine_Version = false;
            continue;
        }
        else if (isLine_VariableNames){
            ReadConfigLine(line,vars);
            isLine_VariableNames = false;
            continue;
        }
        else if (isLine_ConvFactors){
            ReadConfigLine(line,convFactors_str);
            ProcessConvFactors();
            isLine_ConvFactors = false;
            continue;
        }
        else if (isLine_Units){
            ReadConfigLine(line,units);
            isLine_Units = false;
            continue;
        }
            
        // Read Plot Formats
        if(line[0] == 'h'){ 
            ReadHistogramFormat(line);
            continue;
        }
        else if( line[0] == 'g'){ 
            ReadGraphFormat(line);
            continue;
        }
        else if( line[0] == 'f'){
            ReadFrequencyFormat(line);
            continue;
        }
    }

}
        

/*
    Macro Function for reading and processing single information line 
*/
void Config::ReadConfigLine(string line, vector<string>* v)
{
    if(isDebugging){
        cout<<"----"<<endl;
        cout<<"Reading Info Line:"<<endl; 
        cout<<line<<"\n\n";
    }
    
    stringstream line_stream(line);
    string word;
    
    while(line_stream >> word){
        if(isDebugging) cout<<word<<" | ";
        v->push_back(word);
    }
    
    if(isDebugging){
        cout<<"Saved Vector: "<<endl;
        PrintVector(v);
    }
}

/*
    Conversion Factors read as strings
        This function converts strings to double
*/
void Config::ProcessConvFactors()
{
    if ( convFactors_str->size() == 0 ){
        cout<<"ERROR: No Conversion Factors to Process! Did you read data file correctly!"<<endl;
        exit(EXIT_FAILURE);
    }
    
    for (unsigned int i = 0; i < convFactors_str->size(); i++){
        TString s((*convFactors_str)[i]);
        convFactors->push_back(s.Atof());
    }
    
    if ( convFactors_str->size() != convFactors->size() ){
        cout<<"ERROR: Conversion Factor Process!"<<endl;
        exit(EXIT_FAILURE);
    }
}

/*
    Reads and Adds Format to Hist Object
*/
void Config::ReadHistogramFormat(string line)
{
    if(isDebugging){
        cout<<"----"<<endl;
        cout<<"Reading Histogram Format Line:"<<endl; 
        cout<<line<<"\n\n";
    }
    stringstream line_stream(line);
    string word;
    string key; // First word is the Key which is not saved
    string title;
    string var_name;
    int nbins;
    double low;
    double high;
        
    
    line_stream >> key 
                >> title >> var_name
                >> nbins >> low >> high;
    
    hist->AddFormat(title, var_name, nbins, low, high);
   
}

/*
    See CAMAC_DataReader/PlotsFormats.h for the structure of format_graph
*/
void Config::ReadGraphFormat(string line)
{
    if(isDebugging){
        cout<<"----"<<endl;
        cout<<"Reading Graph Format Line:"<<endl; 
        cout<<line<<"\n\n";
    }
    stringstream line_stream(line);
    string word;
    string key; // First word is the Key which is not saved
    string title;
    string x_axis;
    string y_axis;
        
    line_stream >> key 
                >> title >> x_axis >>y_axis;
      
    graph->AddFormat(title, x_axis, y_axis);
   
}

/*
    See CAMAC_DataReader/PlotsFormats.h for the structure of format_frequency
*/
void Config::ReadFrequencyFormat(string line)
{
    if(isDebugging){
        cout<<"----"<<endl;
        cout<<"Reading Frequency Format Line:"<<endl; 
        cout<<line<<"\n\n";
    }
    stringstream line_stream(line);
    string word;
    string key; // First word is the Key which is not saved
    string title;
    vector<int> vars;
    int value;
     
    
    line_stream >> key 
                >> title;
                
    while(line_stream >> value){
        vars.push_back(value);
    }
    
    freq->AddFormat(title, vars);   
}



void Config::PrintVector(vector<string>* v)
{
    for(unsigned int i = 0; i < v->size(); i++){
        cout<<i<<" "<<(*v)[i]<<endl;   
    }
}


#endif
