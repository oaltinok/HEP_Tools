#ifndef Plotter_h
#define Plotter_h

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>

// ROOT Libraries
#include <TROOT.h>
#include <TFile.h>
#include <TH1.h>
#include <TString.h>
#include <TVectorD.h>
#include <TGraph.h>
#include <TCanvas.h>

using namespace std;



class Plotter{

    public:
        Plotter();
        void SetRootFile(TFile* input);
        
        bool isDebugging;
        bool isNewSubrun;
        bool isNewTimeStamp;
        
        TFile* f_Root;

    private:
        

};



#endif
