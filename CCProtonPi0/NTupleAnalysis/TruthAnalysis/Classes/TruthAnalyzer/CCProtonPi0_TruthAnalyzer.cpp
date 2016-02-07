#ifndef CCProtonPi0_TruthAnalyzer_cpp
#define CCProtonPi0_TruthAnalyzer_cpp

#include "CCProtonPi0_TruthAnalyzer.h"

using namespace std;

// Initialize Constants
const double CCProtonPi0_TruthAnalyzer::MeV_to_GeV = pow(10,-3);
const double CCProtonPi0_TruthAnalyzer::MeVSq_to_GeVSq = pow(10,-6);
const double CCProtonPi0_TruthAnalyzer::mm_to_cm = pow(10,-1);
const double CCProtonPi0_TruthAnalyzer::rad_to_deg = 180.0/M_PI;

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

    // Disable Branches for Performance
    fChain->SetBranchStatus("*", false);
    fChain->SetBranchStatus("truth_is*", true);
    fChain->SetBranchStatus("truth_pi0_*", true);
    fChain->SetBranchStatus("truth_muon_*", true);
    fChain->SetBranchStatus("mc_Q2", true);
    fChain->SetBranchStatus("mc_incomingE", true);
    fChain->SetBranchStatus("wgt", true);

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
        if (jentry % 1000000 == 0 ) cout<<"Entry = "<<jentry<<endl;

        if (applyMaxEvents && jentry == nMaxEvents){
            cout<<"\tReached Event Limit!"<<endl;
            break;
        }

        nAll++;
        if (truth_isFidVol) nFidVol++;
        else{
            nNoFidVol++;
            continue;
        }


        // Count Signal and Background
        if (truth_isSignal){
            FillSignalHistograms();
            nSignal++;
        }

        // Background With Pi0
        if (truth_isBckg_NoPi0) nBckg_NoPi0++;
        else if (truth_isBckg_SinglePi0) nBckg_SinglePi0++;
        else if (truth_isBckg_MultiPi0) nBckg_MultiPi0++;

        // Background Types
        if (truth_isBckg_NC) nBckg_NC++;
        else if (truth_isBckg_AntiNeutrino) nBckg_AntiNeutrino++;
        else if (truth_isBckg_QELike) nBckg_QELike++;
        else if (truth_isBckg_SingleChargedPion) nBckg_SingleChargedPion++;
        else if (truth_isBckg_SingleChargedPion_ChargeExchanged) nBckg_SingleChargedPion_ChargeExchanged++;
        else if (truth_isBckg_DoublePionWithPi0) nBckg_DoublePion_WithPi0++;
        else if (truth_isBckg_DoublePionWithoutPi0) nBckg_DoublePion_WithoutPi0++;
        else if (truth_isBckg_MultiPionWithPi0) nBckg_MultiPion_WithPi0++;
        else if (truth_isBckg_MultiPionWithoutPi0) nBckg_MultiPion_WithoutPi0++;
        else if (truth_isBckg_Other) nBckg_Other++;
        else if (truth_isFidVol && !truth_isSignal) cout<<"WARNING! No Background Type"<<endl;
    }

    writeTextFile();
    writeHistograms();
}

void CCProtonPi0_TruthAnalyzer::writeTextFile() 
{
    double totalBackgroundWithPi0 = nBckg_NoPi0 + nBckg_SinglePi0 + nBckg_MultiPi0;
    double totalBackground = nBckg_NC + nBckg_AntiNeutrino + nBckg_QELike + nBckg_SingleChargedPion + nBckg_SingleChargedPion_ChargeExchanged + nBckg_DoublePion_WithPi0 + nBckg_DoublePion_WithoutPi0 + nBckg_MultiPion_WithPi0 + nBckg_MultiPion_WithoutPi0 + nBckg_Other;
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
    textFile<<"nBckg_SingleChargedPion = "<<nBckg_SingleChargedPion<<" "<<GetPercent(nFidVol,nBckg_SingleChargedPion)<<"%"<<endl;
    textFile<<"nBckg_SingleChargedPion_ChargeExchanged = "<<nBckg_SingleChargedPion_ChargeExchanged<<" "<<GetPercent(nFidVol,nBckg_SingleChargedPion_ChargeExchanged)<<"%"<<endl;
    textFile<<"nBckg_DoublePion_WithPi0 = "<<nBckg_DoublePion_WithPi0<<" "<<GetPercent(nFidVol,nBckg_DoublePion_WithPi0)<<"%"<<endl;
    textFile<<"nBckg_DoublePion_WithoutPi0 = "<<nBckg_DoublePion_WithoutPi0<<" "<<GetPercent(nFidVol,nBckg_DoublePion_WithoutPi0)<<"%"<<endl;
    textFile<<"nBckg_MultiPion_WithPi0 = "<<nBckg_MultiPion_WithPi0<<" "<<GetPercent(nFidVol,nBckg_MultiPion_WithPi0)<<"%"<<endl;
    textFile<<"nBckg_MultiPion_WithoutPi0 = "<<nBckg_MultiPion_WithoutPi0<<" "<<GetPercent(nFidVol,nBckg_MultiPion_WithoutPi0)<<"%"<<endl;
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

    // Open ROOT File
    rootDir = Folder_List::rootDir_Truth_mc;

    cout<<"\tRoot File: "<<rootDir<<endl;
    f = new TFile(rootDir.c_str(),"RECREATE");

    initHistograms();

    openTextFiles();

    resetCounters();

    cout<<"Finished"<<endl;
}

void CCProtonPi0_TruthAnalyzer::resetCounters() 
{
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
    nBckg_SingleChargedPion = 0.0;
    nBckg_SingleChargedPion_ChargeExchanged = 0.0;
    nBckg_DoublePion_WithPi0 = 0.0;
    nBckg_DoublePion_WithoutPi0 = 0.0;
    nBckg_MultiPion_WithPi0 = 0.0;
    nBckg_MultiPion_WithoutPi0 = 0.0;
    nBckg_Other = 0.0;
}


void CCProtonPi0_TruthAnalyzer::openTextFiles() 
{
    // Open TextFiles
    file_name = "TruthInfo.dat";
    textFile.open(file_name.c_str());
    if (!textFile.is_open()){
        cerr<<"Cannot open output text file: "<<file_name<<endl;
        exit(1);
    }else{
        cout<<"\t"<<file_name<<endl;
    }
}

void CCProtonPi0_TruthAnalyzer::initHistograms()
{
    // ------------------------------------------------------------------------
    // Muon Variables
    // ------------------------------------------------------------------------
    int nBins_muon_P = 10;
    double min_muon_P = 0.0;
    double max_muon_P = 10.0;
    all_signal_muon_P = new TH1D( "all_signal_muon_P","Muon Momentum for Signal Events",nBins_muon_P, min_muon_P, max_muon_P);
    all_signal_muon_P->GetXaxis()->SetTitle("Momentum [GeV]");
    all_signal_muon_P->GetYaxis()->SetTitle("N(Events)");

    int nBins_muon_theta = 12;
    double min_muon_theta = 0.0;
    double max_muon_theta = 25.0;
    all_signal_muon_theta = new TH1D( "all_signal_muon_theta","Pi0 Muon Theta for Signal Events",nBins_muon_theta,min_muon_theta,max_muon_theta);
    all_signal_muon_theta->GetXaxis()->SetTitle("Theta");
    all_signal_muon_theta->GetYaxis()->SetTitle("N(Events)");

    // ------------------------------------------------------------------------
    // Pi0 Variables
    // ------------------------------------------------------------------------
    double binsP[11] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.5};
    all_signal_pi0_P = new TH1D( "all_signal_pi0_P","Pi0 Momentum for Signal Events",10,binsP);
    all_signal_pi0_P->GetXaxis()->SetTitle("Momentum [GeV]");
    all_signal_pi0_P->GetYaxis()->SetTitle("N(Events)");

    int nBins_pi0_theta = 18;
    double min_pi0_theta = 0.0;
    double max_pi0_theta = 180.0;
    all_signal_pi0_theta = new TH1D( "all_signal_pi0_theta","Theta for Signal Events",nBins_pi0_theta, min_pi0_theta, max_pi0_theta);
    all_signal_pi0_theta->GetXaxis()->SetTitle("Theta");
    all_signal_pi0_theta->GetYaxis()->SetTitle("N(Events)");

    // ------------------------------------------------------------------------
    // Neutrino Energy & Q2
    // ------------------------------------------------------------------------
    int nBins_neutrino_E = 20;
    double min_neutrino_E = 0.0;
    double max_neutrino_E = 20.0;
    all_signal_neutrino_E = new TH1D( "all_signal_neutrino_E","Neutrino Energy for Signal Events",nBins_neutrino_E, min_neutrino_E, max_neutrino_E);
    all_signal_neutrino_E->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
    all_signal_neutrino_E->GetYaxis()->SetTitle("N(Events)");

    int nBins_QSq = 40;
    double min_QSq = 0.0;
    double max_QSq = 4.0;
    all_signal_QSq = new TH1D( "all_signal_QSq","Q^{2} for Signal Events",nBins_QSq,min_QSq,max_QSq);
    all_signal_QSq->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    all_signal_QSq->GetYaxis()->SetTitle("N(Events)");
}

void CCProtonPi0_TruthAnalyzer::FillHistogram(TH1D* hist, double var)
{
    hist->Fill(var, wgt);
}

void CCProtonPi0_TruthAnalyzer::FillSignalHistograms()
{
    FillHistogram(all_signal_muon_P, truth_muon_P * MeV_to_GeV);
    FillHistogram(all_signal_muon_theta, truth_muon_theta * rad_to_deg);
    FillHistogram(all_signal_pi0_P, truth_pi0_P * MeV_to_GeV);
    FillHistogram(all_signal_pi0_theta, truth_pi0_theta * rad_to_deg);
    FillHistogram(all_signal_neutrino_E, mc_incomingE * MeV_to_GeV);
    FillHistogram(all_signal_QSq, mc_Q2 * MeVSq_to_GeVSq);
}

void CCProtonPi0_TruthAnalyzer::writeHistograms()
{
    all_signal_muon_P->Write();
    all_signal_muon_theta->Write();
    all_signal_pi0_P->Write();
    all_signal_pi0_theta->Write();
    all_signal_neutrino_E->Write();
    all_signal_QSq->Write();

    f->Close();
}

#endif


