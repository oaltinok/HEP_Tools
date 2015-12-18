/*
================================================================================
Class: PlotUtils::POTCounter 
    A Class designed to Count POT_Used for a given 
        1) TChain
        2) Playlist - A list of root files

Example Usage:
    std::string playlist = "Input/Playlists/pl_MC_All.dat";
    
    POTCounter pot_counter;
    double totalPOT = pot_counter.getPOTfromPlaylist(playlist);
    
    std::cout<<"Total POT = "<<totalPOT<<std::endl;

Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
================================================================================
*/

#ifndef POTCounter_h
#define POTCounter_h 1

// C++ Libraries
#include <string>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <cstdlib>

// ROOT Libraries
#include <TROOT.h>
#include <TFile.h>
#include <TChain.h>
#include <TChainElement.h>
#include <TCollection.h>
#include <TLeaf.h>
#include <TMath.h>
#include <TObjArray.h>

namespace PlotUtils{

    class POTCounter {

        public:
            // Default Constructor
            POTCounter();
            
            double getPOTfromTChain(TChain* ch);
            double getPOTfromPlaylist(std::string playlist);

        private:
            void Init(std::string playlist, TChain* fChain);

            // List of Variables in Meta Tree 
            Double_t        POT_Total;
            Double_t        POT_Used;
            Int_t           nEntries_Header;
            Int_t           nEntries_Truth;

            // List of branches
            TBranch        *b_POT_Total;   //!
            TBranch        *b_POT_Used;   //!
            TBranch        *b_nEntries_Header;   //!
            TBranch        *b_nEntries_Truth;   //!

    }; // end of class POTCounter
} // end of namespace PlotUtils

#endif

