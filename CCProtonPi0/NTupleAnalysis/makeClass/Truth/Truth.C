#define Truth_Branch_cxx
#include "Truth_Branch.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <iostream>
#include <iomanip>
#include <fstream>

using namespace std;

void Truth_Branch::Loop()
{
    if (fChain == 0) return;
   
    // Open Text File
    ofstream textFile;
    textFile.open("TruthInfo.dat");
    
    // Reset Counters
    double nAll = 0.0;
    double nFidVol = 0.0;
    double nNoFidVol = 0.0;
    double nSignal = 0.0;
    double nBckg_QELike = 0.0;
    double nBckg_SinglePiPlus = 0.0;
    double nBckg_SinglePiMinus = 0.0;
    double nBckg_MultiPion = 0.0;
    double nBckg_MultiPiZero = 0.0;
    double nBckg_Other = 0.0;

    // Control Flow
    bool applyMaxEvents = false;
    double maxEvents = 1000000;

    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
    
    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        
        if (jentry % 250000 == 0 ) cout<<"Entry = "<<jentry<<endl;
        
        if (applyMaxEvents && jentry >= maxEvents) break;
        
        nAll++;
        if (truth_isFidVol) nFidVol++;
        else nNoFidVol++;

        // Count Signal and Background
        if (truth_isSignal) nSignal++;
        if (truth_isBckg_QELike) nBckg_QELike++;
        if (truth_isBckg_SinglePiPlus) nBckg_SinglePiPlus++;
        if (truth_isBckg_SinglePiMinus) nBckg_SinglePiMinus++;
        if (truth_isBckg_MultiPion) nBckg_MultiPion++;
        if (truth_isBckg_MultiPiZero) nBckg_MultiPiZero++;
        if (truth_isBckg_Other) nBckg_Other++;

    }
    
    double totalAnalyzedEvent = nSignal + nBckg_QELike + nBckg_SinglePiPlus + nBckg_SinglePiMinus + nBckg_MultiPion + nBckg_MultiPiZero + nBckg_Other;

    // Formatting for Text Output
    textFile<<fixed;
    textFile<<setprecision(2);

    // Write to Text File
    textFile<<"nAll = "<<nAll<<endl;
    textFile<<"nFidVol = "<<nFidVol<<" "<<GetPercent(nAll,nFidVol)<<"%"<<endl;
    textFile<<"nNoFidVol = "<<nNoFidVol<<" "<<GetPercent(nAll,nNoFidVol)<<"%"<<endl;
    textFile<<"----------------------------------------------------------------"<<endl;
    textFile<<"nSignal = "<<nSignal<<" "<<GetPercent(nFidVol,nSignal)<<"%"<<endl;
    textFile<<"nBckg_QELike = "<<nBckg_QELike<<" "<<GetPercent(nFidVol,nBckg_QELike)<<"%"<<endl;
    textFile<<"nBckg_SinglePiPlus = "<<nBckg_SinglePiPlus<<" "<<GetPercent(nFidVol,nBckg_SinglePiPlus)<<"%"<<endl;
    textFile<<"nBckg_SinglePiMinus = "<<nBckg_SinglePiMinus<<" "<<GetPercent(nFidVol,nBckg_SinglePiMinus)<<"%"<<endl;
    textFile<<"nBckg_MultiPion = "<<nBckg_MultiPion<<" "<<GetPercent(nFidVol,nBckg_MultiPion)<<"%"<<endl;
    textFile<<"nBckg_MultiPiZero = "<<nBckg_MultiPiZero<<" "<<GetPercent(nFidVol,nBckg_MultiPiZero)<<"%"<<endl;
    textFile<<"nBckg_Other = "<<nBckg_Other<<" "<<GetPercent(nFidVol,nBckg_Other)<<"%"<<endl;
    textFile<<"Total Signal & Background = "<<totalAnalyzedEvent<<" "<<GetPercent(nFidVol,totalAnalyzedEvent)<<"%"<<endl;

    textFile.close();
}

double Truth_Branch::GetPercent(double nAll, double nOther)
{
     double percent = (nOther/nAll) * 100;
     return percent;
}





