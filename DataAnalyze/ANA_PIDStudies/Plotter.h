#ifndef PLOTTER_H
#define PLOTTER_H

using namespace std;

#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h>


class Plotter
{
    public:
    
        Plotter();
        
        // -------------------------------------------------------------------------
        //     void plotHistograms: Generates plots for the output of run() function
        //         mcFile-> name of the .root file
        //         plotDir -> folder name for the plots will be created         
        //--------------------------------------------------------------------------
        void plotHistograms(char* mcFile, string plotDir);
        void plotSingleHist_1D(TH1F* singleHist, string fileName, string plotDir);
        void plot2DHist(TH2F* hist2D, string fileName, string plotDir);
    
    private:
        
    


};




#endif

