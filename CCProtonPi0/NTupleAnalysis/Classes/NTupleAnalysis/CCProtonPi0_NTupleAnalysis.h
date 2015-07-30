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
        
        // Constants 
        static const int nHistograms = 13;
        static const int nTopologies = 2;
        static const double SENTINEL = -9.9;
        
    private:
};


#endif
