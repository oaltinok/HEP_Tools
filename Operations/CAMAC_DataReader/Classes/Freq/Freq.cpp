#ifndef Freq_cpp
#define Freq_cpp

#include "Freq.h"


void Freq::WriteRootFile()
{
    // -------------------------------------------------------------------------
    // Create Frequency Plots
    // -------------------------------------------------------------------------
    TGraph* tempFreq;
    
    for (unsigned int i = 0; i < freq_x_axis_vectors.size(); i++){
        // Construt TGraph(size, &x_array_first_element, &y_array_first_element)
        tempFreq = new TGraph(freq_x_axis_vectors[i].size(),&freq_x_axis_vectors[i][0],&freq_y_axis_vectors[i][0]);
        tempFreq->SetName(formats[i].title.c_str());
        tempFreq->SetTitle(formats[i].title.c_str());
        tempFreq->GetXaxis()->SetTitle("Variable Codes");
        tempFreq->GetYaxis()->SetTitle("N(Events)");
        freq_plots.push_back(tempFreq);
    }
    
    // Write Frequency Plots
    for (unsigned int i = 0; i < freq_plots.size(); i++){
        freq_plots[i]->Write("",TObject::kOverwrite); 
    }
    
       
}

void Freq::Fill(int ind)
{
    if ( formats.size() == 0 ){
        cout<<"ERROR: No Frequency Plots Initialized... Did you initilize correctly!"<<endl;
        exit(EXIT_FAILURE);
    }
       
    // Check Frequency Format's vars array - Which Holds the information
    for(unsigned int i = 0; i< formats.size(); i++){
        if (formats[i].vars[ind] == 1 ){
            freq_y_axis_vectors[i][ind] = freq_y_axis_vectors[i][ind]+1;   
        }
    }

}


/*
    initFreqVectors()
        Initializes Frequency Vectors used to generate TGraphs before writing file
        Uses format vector to get Frequency Plots Information
*/
void Freq::Init()
{
    // Sanity Check
    if ( formats.size() == 0 ){
        cout<<"ERROR: No Frequency Plot Format. Did You Read File Correctly?"<<endl;
        exit(EXIT_FAILURE);
    }
    
    vector<double> temp_vector;

    // Fill freq_x_axis_vectors and freq_y_axis_vectors with empty vectors
    for (unsigned int i = 0; i < formats.size(); i++){
        freq_x_axis_vectors.push_back(temp_vector);
        freq_y_axis_vectors.push_back(temp_vector);
    }
    
    // Initialize Axis Vectors
    //      x_axis = Variable Indice
    //      y_axis = All Zeros (will be incremented inside ProcessDataLine)
    for(unsigned int i = 0; i < freq_x_axis_vectors.size(); i++){
        for(unsigned int j = 0; j < formats[i].vars.size(); j++){
            freq_x_axis_vectors[i].push_back(j);
            freq_y_axis_vectors[i].push_back(0.0);
        }
    }
    
    if (freq_x_axis_vectors.size() == 0 || freq_y_axis_vectors.size() == 0){
        cout<<"ERROR: Frequency Vectors Initialization"<<endl;
        exit(EXIT_FAILURE);
    }
    
}

void Freq::AddFormat(string in_title, vector<int> in_vars)
{
    Format tempFormat;
    
    tempFormat.title = in_title;
    tempFormat.vars = in_vars;

    formats.push_back(tempFormat);
}



#endif
