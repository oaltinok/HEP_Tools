#ifndef Hist_cpp
#define Hist_cpp

#include "Hist.h"

using namespace std;


/*
    initHistograms()
        Initializes Histograms inside hists vector
        Uses f_hists vector to get Histogram Information
*/
void Hist::Init()
{
    // Sanity Check
    if ( formats.size() == 0 ){
        cout<<"ERROR: No Histogram Format. Did You Read File Correctly?"<<endl;
        exit(EXIT_FAILURE);
    }
    
    TH1D* tempHist;
    string title;
    string var_name;
    int nbins;
    double low;
    double high;
   
    // -------------------------------------------------------------------------
    // Get Existing Histogram or Create a NEW one
    // -------------------------------------------------------------------------
    for (unsigned int i = 0; i < formats.size(); i++){
        
        // Get Histogram Format Info
        title = formats[i].title;
        var_name = formats[i].var_name;
        nbins = formats[i].nbins;
        low = formats[i].low;
        high = formats[i].high;
      
        // Get Existing Histogram
        tempHist = (TH1D*)gDirectory->Get(var_name.c_str());
        
        // If NO Existing Histogram
        //    Create a New One
        if(tempHist == NULL){
            tempHist = new TH1D( var_name.c_str(),title.c_str(),nbins, low, high );
            tempHist->SetBins(nbins, low, high);
            tempHist->GetXaxis()->SetTitle(var_name.c_str());
            tempHist->GetYaxis()->SetTitle("N(Events)");
        }
        
        // Push Histogram to hists Vector
        hists.push_back(tempHist);
    }
    
    if (hists.size() == 0 ){
        cout<<"ERROR: Histogram Initialization"<<endl;
        exit(EXIT_FAILURE);
    }
}


void Hist::AddFormat(   string in_title, string in_var_name,
                        int in_nbins, double in_low, double in_high)
{
    Format tempFormat;
    
    tempFormat.title = in_title;
    tempFormat.var_name = in_var_name;
    tempFormat.nbins = in_nbins;
    tempFormat.low = in_low;
    tempFormat.high = in_high;

    formats.push_back(tempFormat);
}

void Hist::Fill(string var_name, double value)
{
    if ( formats.size() == 0 ){
        cout<<"ERROR: No Histograms Initialized... Did you initilize correctly!"<<endl;
        exit(EXIT_FAILURE);
    }
    
    // Search Histogram Format Vector for the given var_name
    //      Fill Corresponding Histogram if that variable has a histogram
    for ( unsigned int i = 0; i < formats.size(); i++){

        if (var_name.compare(formats[i].var_name) == 0) {
            hists[i]->Fill(value);
            break;
        }
    }
}

void Hist::WriteRootFile()
{
    for(unsigned int i = 0; i < hists.size(); i++){
        hists[i]->Write("",TObject::kOverwrite);   
    }
}




#endif
