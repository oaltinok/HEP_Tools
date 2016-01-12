#ifndef CCProtonPi0_TruthAnalyzer_cpp
#define CCProtonPi0_TruthAnalyzer_cpp

#include "CCProtonPi0_TruthAnalyzer.h"

using namespace std;

void CCProtonPi0_TruthAnalyzer::Loop(string playlist)
{
    // Control Flow
    bool applyMaxEvents = false;
    double nMaxEvents = 100000;
    
    //------------------------------------------------------------------------
    // Create chain
    //------------------------------------------------------------------------
    TChain* fChain = new TChain("Truth");
    Init(playlist, fChain);

    if (!fChain) return;
    if (fChain == 0) return;

    //------------------------------------------------------------------------
    // Loop over Chain
    //------------------------------------------------------------------------
    cout<<"Looping over all entries"<<endl;

    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;

    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        
        nb = fChain->GetEntry(jentry);   nbytes += nb;    
        Long64_t ientry = fChain->GetEntry(jentry);

        if (ientry == 0) {
            cout<<"\tGetEntry failure "<<jentry<<endl;
            break;
        }

        // Progress Message on Terminal
        if (jentry % 500000 == 0 ) cout<<"Entry = "<<jentry<<endl;

        if (applyMaxEvents && jentry == nMaxEvents){
            cout<<"\tReached Event Limit!"<<endl;
            break;
        }

        nAll++;
        if (truth_isFidVol) nFidVol++;
        else nNoFidVol++;

        // Count Signal and Background
        if (truth_isSignal) nSignal++;

        // Background With Pi0
        if (truth_isBckg_NoPi0) nBckg_NoPi0++;
        else if (truth_isBckg_SinglePi0) nBckg_SinglePi0++;
        else if (truth_isBckg_MultiPi0) nBckg_MultiPi0++;

        // Background Types
        if (truth_isBckg_NC) nBckg_NC++;
        else if (truth_isBckg_AntiNeutrino) nBckg_AntiNeutrino++;
        else if (truth_isBckg_QELike) nBckg_QELike++;
        else if (truth_isBckg_SinglePion) nBckg_SinglePion++;
        else if (truth_isBckg_DoublePion) nBckg_DoublePion++;
        else if (truth_isBckg_MultiPion) nBckg_MultiPion++;
        else if (truth_isBckg_Other) nBckg_Other++;

    }

    double totalBackgroundWithPi0 = nBckg_NoPi0 + nBckg_SinglePi0 + nBckg_MultiPi0;
    double totalBackground = nBckg_NC + nBckg_AntiNeutrino + nBckg_QELike + nBckg_SinglePion + nBckg_DoublePion + nBckg_MultiPion + nBckg_Other;
    double totalEvents1 = nSignal + totalBackgroundWithPi0;
    double totalEvents2 = nSignal + totalBackground;

    // Formatting for Text Output
    textFile<<fixed;
    textFile<<setprecision(2);

    // Write to Text File
    textFile<<"nAll = "<<nAll<<endl;
    textFile<<"nFidVol = "<<nFidVol<<" "<<GetPercent(nAll,nFidVol)<<"%"<<endl;
    textFile<<"nNoFidVol = "<<nNoFidVol<<" "<<GetPercent(nAll,nNoFidVol)<<"%"<<endl;
    textFile<<"================================================================"<<endl;
    textFile<<"nSignal = "<<nSignal<<" "<<GetPercent(nFidVol,nSignal)<<"%"<<endl;
    textFile<<"----------------------------------------------------------------"<<endl;
    textFile<<"nBckg_NoPi0 = "<<nBckg_NoPi0<<" "<<GetPercent(nFidVol,nBckg_NoPi0)<<"%"<<endl;
    textFile<<"nBckg_SinglePi0 = "<<nBckg_SinglePi0<<" "<<GetPercent(nFidVol,nBckg_SinglePi0)<<"%"<<endl;
    textFile<<"nBckg_MultiPi0 = "<<nBckg_MultiPi0<<" "<<GetPercent(nFidVol,nBckg_MultiPi0)<<"%"<<endl;
    textFile<<"----------------------------------------------------------------"<<endl;
    textFile<<"nBckg_NC = "<<nBckg_NC<<" "<<GetPercent(nFidVol,nBckg_NC)<<"%"<<endl;
    textFile<<"nBckg_AntiNeutrino = "<<nBckg_AntiNeutrino<<" "<<GetPercent(nFidVol,nBckg_AntiNeutrino)<<"%"<<endl;
    textFile<<"nBckg_QELike = "<<nBckg_QELike<<" "<<GetPercent(nFidVol,nBckg_QELike)<<"%"<<endl;
    textFile<<"nBckg_SinglePion = "<<nBckg_SinglePion<<" "<<GetPercent(nFidVol,nBckg_SinglePion)<<"%"<<endl;
    textFile<<"nBckg_DoublePion = "<<nBckg_DoublePion<<" "<<GetPercent(nFidVol,nBckg_DoublePion)<<"%"<<endl;
    textFile<<"nBckg_MultiPion = "<<nBckg_MultiPion<<" "<<GetPercent(nFidVol,nBckg_MultiPion)<<"%"<<endl;
    textFile<<"nBckg_Other = "<<nBckg_Other<<" "<<GetPercent(nFidVol,nBckg_Other)<<"%"<<endl;
    textFile<<"----------------------------------------------------------------"<<endl;

    textFile<<"Total Signal & Background With Pi0 = "<<totalEvents1<<" "<<GetPercent(nFidVol,totalEvents1)<<"%"<<endl;
    textFile<<"Total Signal & Background = "<<totalEvents2<<" "<<GetPercent(nFidVol,totalEvents2)<<"%"<<endl;

    textFile.close();
}

double CCProtonPi0_TruthAnalyzer::GetPercent(double nAll, double nOther)
{
    double percent = (nOther/nAll) * 100;
    return percent;
}

CCProtonPi0_TruthAnalyzer::CCProtonPi0_TruthAnalyzer()
{
    cout<<"Initializing TruthAnalyzer!"<<endl;
    
    // Open TextFiles
    file_name = "TruthInfo.dat";
    textFile.open(file_name.c_str());
    if (!textFile.is_open()){
        cerr<<"Cannot open output text file: "<<file_name<<endl;
        exit(1);
    }else{
        cout<<"\t"<<file_name<<endl;
    }

    // Reset Counters
    nAll = 0.0;
    nFidVol = 0.0;
    nNoFidVol = 0.0;
    nSignal = 0.0;
    // Background With Pi0
    nBckg_NoPi0 = 0.0;
    nBckg_SinglePi0 = 0.0;
    nBckg_MultiPi0 = 0.0;
    // Background Types
    nBckg_NC = 0.0;
    nBckg_AntiNeutrino = 0.0;
    nBckg_QELike = 0.0;
    nBckg_SinglePion = 0.0;
    nBckg_DoublePion = 0.0;
    nBckg_MultiPion = 0.0;
    nBckg_Other = 0.0;

    cout<<"Finished"<<endl;
}

#endif


