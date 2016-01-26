/*
================================================================================
Class: CCProtonPi0_NTupleAnalysis
    Base Class for the CCProtonPi0_NTupleAnalysis Package
    Includes Common Member Variables for all Classes
        
    Author:        Ozgur Altinok  - ozgur.altinok@tufts.edu
================================================================================
*/

#ifndef CCProtonPi0_NTupleAnalysis_h
#define CCProtonPi0_NTupleAnalysis_h

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>

#include "Cintex/Cintex.h"
#include <TROOT.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <PlotUtils/MnvH1D.h>
#include <PlotUtils/MnvH2D.h>

// Libraries
#include "../../Libraries/PDG_List.h"
#include "../../Libraries/Folder_List.h"
#include "../../Libraries/Data_Functions.h"
#include "../../Libraries/HEP_Functions.h"

using namespace std;

class CCProtonPi0_NTupleAnalysis
{
    public:
        CCProtonPi0_NTupleAnalysis();
       
        string ana_folder;
        // Constants 
        static const int nHistograms = 16;
        static const int nTopologies = 2;
        
        static const string version;
        static const double SENTINEL;
        static const double MeV_to_GeV; 
        static const double MeVSq_to_GeVSq;
        static const double mm_to_cm;
        static const double rad_to_deg;

        void OpenTextFile(string file_name, std::ofstream &file);

    private:
};


#endif
