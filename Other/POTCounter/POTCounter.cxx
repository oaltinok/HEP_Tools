#ifndef POTCounter_cxx
#define POTCounter_cxx 1

#include "POTCounter.h"

using namespace PlotUtils;

double POTCounter::getPOTfromPlaylist(std::string playlist)
{
    // Create TChain form Playlist and Initialize
    TChain* fChain = new TChain("Meta");
    Init(playlist, fChain);

    // Check TChain
    if (!fChain || fChain == 0){
        std::cout<<"PlotUtils::POTCounter -- Can NOT find TChain from input playlist! -- Returning -1"<<std::endl;
        return -1;
    }else{
        std::cout<<"PlotUtils::POTCounter -- Initialized fChain using the Input playlist!"<<std::endl;
    }

    // Get POT_Used and return
    double sumPOTUsed = getPOTfromTChain(fChain);

    return sumPOTUsed;
}

double POTCounter::getPOTfromTChain(TChain* ch)
{
    double sumPOTUsed = 0;

    TObjArray* fileElements = ch->GetListOfFiles();
    TIter next(fileElements);
    TChainElement* chEl=0;

    std::cout<<"PlotUtils::POTCounter -- Counting POT"<<std::endl;
    std::cout<<"PlotUtils::POTCounter -- It can take some time depending on the size of the TChain"<<std::endl;

    while (( chEl=(TChainElement*)next() )) {
        TFile f(chEl->GetTitle());
        TTree* t = (TTree*)f.Get("Meta");
        if (!t){
            std::cout<<"PlotUtils::POTCounter -- No Meta tree in file "<<chEl->GetTitle()<<std::endl;
            continue;
        }

        // Loop Over all Entries 
        int n_entries = t->GetEntries();
        for (int i=0;i<n_entries;i++){
            t->GetEntry(i);
            TLeaf* POT_Used = t->GetLeaf("POT_Used");

            if (POT_Used){
                sumPOTUsed = sumPOTUsed + POT_Used->GetValue();
            }
        }
    }

    return sumPOTUsed;
}


/*
 *  Initializes the fChain for a given playlist
 *      Reads the playlist and adds all ROOT files to the fChain
 */
void POTCounter::Init(std::string playlist, TChain* fChain)
{
    ifstream input_pl(playlist.c_str());
    std::string filename;

    if( !input_pl.is_open() ){
        std::cerr<<"PlotUtils::POTCounter -- Cannot open Playlist File!"<<std::endl;
        exit(1);
    }else{
        std::cout<<"PlotUtils::POTCounter -- Reading Playlist: "<<playlist.c_str()<<std::endl;
    }

    /* 
     * Loop input playlist and get file names 
     *     Assumption: file names start with '/' character 
     */
    while (!input_pl.eof()) {
        getline(input_pl,filename);

        if (filename[0] != '/') continue;

        fChain->Add( filename.c_str() );
        //std::cout<<"PlotUtils::POTCounter -- Added "<<filename.c_str()<<std::endl;
    }

    fChain->SetMakeClass(1);
    fChain->SetBranchAddress("POT_Total", &POT_Total, &b_POT_Total);
    fChain->SetBranchAddress("POT_Used", &POT_Used, &b_POT_Used);
    fChain->SetBranchAddress("nEntries_Header", &nEntries_Header, &b_nEntries_Header);
    fChain->SetBranchAddress("nEntries_Truth", &nEntries_Truth, &b_nEntries_Truth);

    input_pl.close();
}

POTCounter::POTCounter()
{
    // Do Nothing!
}

#endif

