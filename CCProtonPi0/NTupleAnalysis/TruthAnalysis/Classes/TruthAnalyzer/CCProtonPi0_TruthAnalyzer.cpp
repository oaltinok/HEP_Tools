#ifndef CCProtonPi0_TruthAnalyzer_cpp
#define CCProtonPi0_TruthAnalyzer_cpp

#include "CCProtonPi0_TruthAnalyzer.h"

using namespace PlotUtils;

void CCProtonPi0_TruthAnalyzer::Loop(std::string playlist)
{
    // Control Flow
    bool applyMaxEvents = false;
    double nMaxEvents = 100000;

    applyBckgConstraints_CV = true;
    applyBckgConstraints_Unv = true;
    fillErrors_ByHand = true;

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
    std::cout<<"Looping over all entries"<<std::endl;

    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;

    for (Long64_t jentry=0; jentry<nentries;jentry++) {

        nb = fChain->GetEntry(jentry);   nbytes += nb;    
        Long64_t ientry = fChain->GetEntry(jentry);

        if (ientry == 0) {
            std::cout<<"\tGetEntry failure "<<jentry<<std::endl;
            break;
        }

        // Progress Message on Terminal
        if (jentry % 1000000 == 0 ) std::cout<<"Entry = "<<jentry<<std::endl;

        if (applyMaxEvents && jentry == nMaxEvents){
            std::cout<<"\tReached Event Limit!"<<std::endl;
            break;
        }

        nAll.increment();

        if (truth_isFidVol) nFidVol.increment();
        else nFidVol_Out.increment();

        // Count Events Only inside the Fidicual Volume
        if (truth_isFidVol){
            FillMichelHistograms();

            if (truth_isSignalOut_Acceptance) nSignalOut_Acceptance.increment();
            if (truth_isSignalOut_Kinematics) nSignalOut_Kinematics.increment();

            // Count Signal and Background inside Fiducial Volume
            if (truth_isSignal) nSignal.increment();
            else nBckg.increment();

            // Count Signal Type & Fill Histograms
            if (truth_isSignal){
                FillSignalHistograms();
            }
        }
    }

    // Add Other Error Bands and Fill With CV
    AddOtherErrorBands_FillWithCV();

    writeTextFile();
    writeHistograms();
}

void CCProtonPi0_TruthAnalyzer::writeTextFile() 
{
    // Formatting for Text Output
    textFile<<std::fixed;
    textFile<<std::setprecision(1);

    // All Events
    WriteCounter(nAll, nAll);
    textFile<<std::endl;
    
    WriteCounter(nFidVol, nAll);
    WriteCounter(nFidVol_Out, nAll);
    textFile<<std::endl;

    // Events inside Fiducial Volume
    WriteCounter(nSignal, nFidVol);
    WriteCounter(nSignalOut_Acceptance, nFidVol);
    WriteCounter(nSignalOut_Kinematics, nFidVol);
    WriteCounter(nBckg, nFidVol);
    textFile<<std::endl;
 
    // Signal Types
    WriteCounter(nQE, nSignal);
    textFile<<std::endl;

    WriteCounter(nRES_1232, nSignal);
    WriteCounter(nRES_1535, nSignal);
    WriteCounter(nRES_1520, nSignal);
    WriteCounter(nRES_Other, nSignal);
    textFile<<std::endl;
   
    WriteCounter(nDIS_1_pi, nSignal);
    WriteCounter(nDIS_2_pi, nSignal);
    WriteCounter(nDIS_Multi_pi, nSignal);
    WriteCounter(nNon_RES, nSignal);
    textFile<<std::endl;
   
    textFile.close();
}

void CCProtonPi0_TruthAnalyzer::WriteCounter(CCProtonPi0_Counter Counter, CCProtonPi0_Counter PercentBase)
{
    textFile<<Counter.getName()<<"\t"<<Counter.getCount()<<"\t"<<GetPercent(PercentBase,Counter)<<std::endl;
}

double CCProtonPi0_TruthAnalyzer::GetPercent(CCProtonPi0_Counter nAll, CCProtonPi0_Counter nOther)
{
    double percent = (nOther.getCount()/nAll.getCount()) * 100;
    return percent;
}

CCProtonPi0_TruthAnalyzer::CCProtonPi0_TruthAnalyzer() : 
    CCProtonPi0_NTupleAnalysis(),
    BckgConstrainer(Folder_List::BckgConstraints_TruthAnalysis)
{
    std::cout<<"Initializing TruthAnalyzer!"<<std::endl;

    // Required for MINERvA Framework Classes
    ROOT::Cintex::Cintex::Enable();

    // Open ROOT File
    rootDir = Folder_List::rootDir_Truth_mc;

    std::cout<<"\tRoot File: "<<rootDir<<std::endl;
    f = new TFile(rootDir.c_str(),"RECREATE");
 
    initHistograms();
    
    openTextFiles();

    resetCounters();
    
    initCVWeights();

    std::cout<<"Finished"<<std::endl;
}

void CCProtonPi0_TruthAnalyzer::resetCounters() 
{
    nAll.setName("nAll");
   
    nFidVol.setName("nFidVol");
    nFidVol_Out.setName("nFidVol_Out");

    nSignal.setName("nSignal");
    nSignalOut_Acceptance.setName("nSignalOut_Acceptance");
    nSignalOut_Kinematics.setName("nSignalOut_Kinematics");
    nBckg.setName("nBckg");

    // Signal Type
    nQE.setName("nQuasi_Elastic");
    
    nRES_1232.setName("nSignal_RES_Delta");
    nRES_1535.setName("nSignal_RES_1535");
    nRES_1520.setName("nSignal_RES_1520");
    nRES_Other.setName("nSignal_RES_Other");
    
    nDIS_1_pi.setName("nSignal_DIS_1pi");
    nDIS_2_pi.setName("nSignal_DIS_2pi");
    nDIS_Multi_pi.setName("nSignal_DIS_Multi_pi");
    nNon_RES.setName("nSignal_Non_RES");
}


void CCProtonPi0_TruthAnalyzer::openTextFiles() 
{
    // Open TextFiles
    file_name = "TruthInfo.dat";
    textFile.open(file_name.c_str());
    if (!textFile.is_open()){
        std::cerr<<"Cannot open output text file: "<<file_name<<std::endl;
        exit(1);
    }else{
        std::cout<<"\t"<<file_name<<std::endl;
    }
}

void CCProtonPi0_TruthAnalyzer::initHistograms()
{
    // Michel Electron
    michel_ZX = new TH2D("michel_ZX","Michel Distance to Vertex",50,-300,300,50,-300,300);
    michel_ZX->GetXaxis()->SetTitle("Michel_Z - Vertex_Z Position [mm]");
    michel_ZX->GetZaxis()->SetTitle("Michel_X - Vertex_X Position [mm]");

    michel_ZY = new TH2D("michel_ZY","Michel Distance to Vertex",50,-300,300,50,-300,300);
    michel_ZY->GetXaxis()->SetTitle("Michel_Z - Vertex_Z Position [mm]");
    michel_ZY->GetYaxis()->SetTitle("Michel_Y - Vertex_Y Position [mm]");

    michel_dist_X = new TH1D("michel_dist_X","Michel Distance to Vertex",50,-300,300);
    michel_dist_X->GetXaxis()->SetTitle("Michel_X - Vertex_X Position [mm]");
    michel_dist_X->GetYaxis()->SetTitle("N(Events)");

    michel_dist_Y = new TH1D("michel_dist_Y","Michel Distance to Vertex",50,-300,300);
    michel_dist_Y->GetXaxis()->SetTitle("Michel_Y - Vertex_Y Position [mm]");
    michel_dist_Y->GetYaxis()->SetTitle("N(Events)");

    michel_dist_Z = new TH1D("michel_dist_Z","Michel Distance to Vertex",50,-300,300);
    michel_dist_Z->GetXaxis()->SetTitle("Michel_Z - Vertex_Z Position [mm]");
    michel_dist_Z->GetYaxis()->SetTitle("N(Events)");

    michel_dist_total = new TH1D("michel_dist_total","Michel Distance to Vertex",50,0,1000);
    michel_dist_total->GetXaxis()->SetTitle("Michel Begin Position - Vertex Position [mm]");
    michel_dist_total->GetYaxis()->SetTitle("N(Events)");

    michel_pionP = new TH1D("michel_pionP","Michel Pion Momentum",50,0,1000);
    michel_pionP->GetXaxis()->SetTitle("Michel Pion Momentum [MeV]");
    michel_pionP->GetYaxis()->SetTitle("N(Events)");

    michel_pion_dist = new TH1D("michel_pion_dist","Michel Pion Distance to Vertex",50,0,1000);
    michel_pion_dist->GetXaxis()->SetTitle("Michel Pion Distance to Vertex [mm]");
    michel_pion_dist->GetYaxis()->SetTitle("N(Events)");

    // ------------------------------------------------------------------------
    // Muon Variables
    // ------------------------------------------------------------------------
    muon_P_mc_truth_all_signal = new MnvH1D( "muon_P_mc_truth_all_signal","Muon Momentum for Signal Events",binList.size_muon_P, binList.a_muon_P);
    muon_P_mc_truth_all_signal->GetXaxis()->SetTitle("Momentum [GeV]");
    muon_P_mc_truth_all_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBand_Flux(muon_P_mc_truth_all_signal);
    AddVertErrorBand_Genie(muon_P_mc_truth_all_signal);

    muon_theta_mc_truth_all_signal = new MnvH1D( "muon_theta_mc_truth_all_signal","Pi0 Muon Theta for Signal Events",binList.size_muon_theta, binList.a_muon_theta);
    muon_theta_mc_truth_all_signal->GetXaxis()->SetTitle("Theta");
    muon_theta_mc_truth_all_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBand_Flux(muon_theta_mc_truth_all_signal);
    AddVertErrorBand_Genie(muon_theta_mc_truth_all_signal);

    // ------------------------------------------------------------------------
    // Pi0 Variables
    // ------------------------------------------------------------------------
    pi0_P_mc_truth_all_signal = new MnvH1D( "pi0_P_mc_truth_all_signal","Pi0 Momentum for Signal Events",binList.size_pi0_P, binList.a_pi0_P);
    pi0_P_mc_truth_all_signal->GetXaxis()->SetTitle("Momentum [GeV]");
    pi0_P_mc_truth_all_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBand_Flux(pi0_P_mc_truth_all_signal);
    AddVertErrorBand_Genie(pi0_P_mc_truth_all_signal);

    pi0_KE_mc_truth_all_signal = new MnvH1D( "pi0_KE_mc_truth_all_signal","Pi0 Kinetic Energy for Signal Events",binList.size_pi0_KE, binList.a_pi0_KE);
    pi0_KE_mc_truth_all_signal->GetXaxis()->SetTitle("Kinetic Energy [GeV]");
    pi0_KE_mc_truth_all_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBand_Flux(pi0_KE_mc_truth_all_signal);
    AddVertErrorBand_Genie(pi0_KE_mc_truth_all_signal);

    pi0_theta_mc_truth_all_signal = new MnvH1D( "pi0_theta_mc_truth_all_signal","Theta for Signal Events",binList.size_pi0_theta, binList.a_pi0_theta);
    pi0_theta_mc_truth_all_signal->GetXaxis()->SetTitle("Theta");
    pi0_theta_mc_truth_all_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBand_Flux(pi0_theta_mc_truth_all_signal);
    AddVertErrorBand_Genie(pi0_theta_mc_truth_all_signal);

    // ------------------------------------------------------------------------
    // Neutrino Energy, Q2 & W
    // ------------------------------------------------------------------------
    Enu_mc_truth_all_signal = new MnvH1D( "Enu_mc_truth_all_signal","Neutrino Energy for Signal Events",binList.size_Enu, binList.a_Enu);
    Enu_mc_truth_all_signal->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
    Enu_mc_truth_all_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBand_Flux(Enu_mc_truth_all_signal);
    AddVertErrorBand_Genie(Enu_mc_truth_all_signal);

    QSq_mc_truth_all_signal = new MnvH1D( "QSq_mc_truth_all_signal","Q^{2} for Signal Events",binList.size_QSq, binList.a_QSq);
    QSq_mc_truth_all_signal->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    QSq_mc_truth_all_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBand_Flux(QSq_mc_truth_all_signal);
    AddVertErrorBand_Genie(QSq_mc_truth_all_signal);

    W_mc_truth_all_signal = new MnvH1D( "W_mc_truth_all_signal","W for Signal Events",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max());
    W_mc_truth_all_signal->GetXaxis()->SetTitle("W [GeV]");
    W_mc_truth_all_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBand_Flux(W_mc_truth_all_signal);
    AddVertErrorBand_Genie(W_mc_truth_all_signal);

    // ------------------------------------------------------------------------
    // Signal Q2
    // ------------------------------------------------------------------------
    mc_Q2_QE = new TH1D("mc_Q2_QE","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    mc_Q2_QE->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    mc_Q2_QE->GetYaxis()->SetTitle("N(Events)");

    mc_Q2_RES_1232 = new TH1D("mc_Q2_RES_1232","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    mc_Q2_RES_1232->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    mc_Q2_RES_1232->GetYaxis()->SetTitle("N(Events)");

    mc_Q2_RES_1535 = new TH1D("mc_Q2_RES_1535","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    mc_Q2_RES_1535->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    mc_Q2_RES_1535->GetYaxis()->SetTitle("N(Events)");

    mc_Q2_RES_1520 = new TH1D("mc_Q2_RES_1520","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    mc_Q2_RES_1520->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    mc_Q2_RES_1520->GetYaxis()->SetTitle("N(Events)");

    mc_Q2_RES_Other = new TH1D("mc_Q2_RES_Other","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    mc_Q2_RES_Other->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    mc_Q2_RES_Other->GetYaxis()->SetTitle("N(Events)");

    mc_Q2_DIS_1_pi = new TH1D("mc_Q2_DIS_1_pi","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    mc_Q2_DIS_1_pi->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    mc_Q2_DIS_1_pi->GetYaxis()->SetTitle("N(Events)");

    mc_Q2_DIS_2_pi = new TH1D("mc_Q2_DIS_2_pi","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    mc_Q2_DIS_2_pi->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    mc_Q2_DIS_2_pi->GetYaxis()->SetTitle("N(Events)");

    mc_Q2_DIS_Multi_pi = new TH1D("mc_Q2_DIS_Multi_pi","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    mc_Q2_DIS_Multi_pi->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    mc_Q2_DIS_Multi_pi->GetYaxis()->SetTitle("N(Events)");

    mc_Q2_Non_RES = new TH1D("mc_Q2_Non_RES","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    mc_Q2_Non_RES->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    mc_Q2_Non_RES->GetYaxis()->SetTitle("N(Events)");

    // ------------------------------------------------------------------------
    // Signal incomingE
    // ------------------------------------------------------------------------
    mc_incomingE_QE = new TH1D("mc_incomingE_QE","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    mc_incomingE_QE->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    mc_incomingE_QE->GetYaxis()->SetTitle("N(Events)");

    mc_incomingE_RES_1232 = new TH1D("mc_incomingE_RES_1232","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    mc_incomingE_RES_1232->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    mc_incomingE_RES_1232->GetYaxis()->SetTitle("N(Events)");

    mc_incomingE_RES_1535 = new TH1D("mc_incomingE_RES_1535","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    mc_incomingE_RES_1535->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    mc_incomingE_RES_1535->GetYaxis()->SetTitle("N(Events)");

    mc_incomingE_RES_1520 = new TH1D("mc_incomingE_RES_1520","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    mc_incomingE_RES_1520->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    mc_incomingE_RES_1520->GetYaxis()->SetTitle("N(Events)");

    mc_incomingE_RES_Other = new TH1D("mc_incomingE_RES_Other","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    mc_incomingE_RES_Other->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    mc_incomingE_RES_Other->GetYaxis()->SetTitle("N(Events)");

    mc_incomingE_DIS_1_pi = new TH1D("mc_incomingE_DIS_1_pi","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    mc_incomingE_DIS_1_pi->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    mc_incomingE_DIS_1_pi->GetYaxis()->SetTitle("N(Events)");

    mc_incomingE_DIS_2_pi = new TH1D("mc_incomingE_DIS_2_pi","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    mc_incomingE_DIS_2_pi->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    mc_incomingE_DIS_2_pi->GetYaxis()->SetTitle("N(Events)");

    mc_incomingE_DIS_Multi_pi = new TH1D("mc_incomingE_DIS_Multi_pi","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    mc_incomingE_DIS_Multi_pi->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    mc_incomingE_DIS_Multi_pi->GetYaxis()->SetTitle("N(Events)");

    mc_incomingE_Non_RES = new TH1D("mc_incomingE_Non_RES","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    mc_incomingE_Non_RES->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    mc_incomingE_Non_RES->GetYaxis()->SetTitle("N(Events)");

    // ------------------------------------------------------------------------
    // Signal w
    // ------------------------------------------------------------------------
    mc_w_QE = new TH1D("mc_w_QE","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    mc_w_QE->GetXaxis()->SetTitle("W [GeV]");
    mc_w_QE->GetYaxis()->SetTitle("N(Events)");
 
    mc_w_RES_1232 = new TH1D("mc_w_RES_1232","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    mc_w_RES_1232->GetXaxis()->SetTitle("W [GeV]");
    mc_w_RES_1232->GetYaxis()->SetTitle("N(Events)");

    mc_w_RES_1535 = new TH1D("mc_w_RES_1535","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    mc_w_RES_1535->GetXaxis()->SetTitle("W [GeV]");
    mc_w_RES_1535->GetYaxis()->SetTitle("N(Events)");

    mc_w_RES_1520 = new TH1D("mc_w_RES_1520","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    mc_w_RES_1520->GetXaxis()->SetTitle("W [GeV]");
    mc_w_RES_1520->GetYaxis()->SetTitle("N(Events)");

    mc_w_RES_Other = new TH1D("mc_w_RES_Other","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    mc_w_RES_Other->GetXaxis()->SetTitle("W [GeV]");
    mc_w_RES_Other->GetYaxis()->SetTitle("N(Events)");

    mc_w_DIS_1_pi = new TH1D("mc_w_DIS_1_pi","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    mc_w_DIS_1_pi->GetXaxis()->SetTitle("W [GeV]");
    mc_w_DIS_1_pi->GetYaxis()->SetTitle("N(Events)");

    mc_w_DIS_2_pi = new TH1D("mc_w_DIS_2_pi","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    mc_w_DIS_2_pi->GetXaxis()->SetTitle("W [GeV]");
    mc_w_DIS_2_pi->GetYaxis()->SetTitle("N(Events)");

    mc_w_DIS_Multi_pi = new TH1D("mc_w_DIS_Multi_pi","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    mc_w_DIS_Multi_pi->GetXaxis()->SetTitle("W [GeV]");
    mc_w_DIS_Multi_pi->GetYaxis()->SetTitle("N(Events)");

    mc_w_Non_RES = new TH1D("mc_w_Non_RES","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    mc_w_Non_RES->GetXaxis()->SetTitle("W [GeV]");
    mc_w_Non_RES->GetYaxis()->SetTitle("N(Events)");

}

void CCProtonPi0_TruthAnalyzer::FillHistogram(MnvH1D* hist, double var)
{
    hist->Fill(var, cvweight);

    if (fillErrors_ByHand){
        FillVertErrorBand_Flux_ByHand(hist, var);
        FillVertErrorBand_Genie_ByHand(hist, var);
    }else{
        FillVertErrorBand_Flux(hist, var);
        FillVertErrorBand_Genie(hist, var);
    }
}

void CCProtonPi0_TruthAnalyzer::FillHistogram(TH1D* hist, double var)
{
    hist->Fill(var, cvweight);
}

void CCProtonPi0_TruthAnalyzer::FillVertErrorBand_Flux(MnvH1D* h, double var)
{
    std::vector<double> flux_errors = GetFluxError(mc_incomingE * MeV_to_GeV, mc_incoming);
    h->FillVertErrorBand("Flux",  var, &flux_errors[0],  cvweight);
}

void CCProtonPi0_TruthAnalyzer::FillVertErrorBand_Flux_ByHand(MnvH1D* h, double var)
{
    std::vector<double> flux_errors = GetFluxError(mc_incomingE * MeV_to_GeV, mc_incoming);
    FillVertErrorBand_ByHand(h, var, "Flux", flux_errors);
}

void CCProtonPi0_TruthAnalyzer::FillVertErrorBand_Genie(MnvH1D* h, double var)
{
    h->FillVertErrorBand("GENIE_AGKYxF1pi"         ,var, truth_genie_wgt_AGKYxF1pi[2]        , truth_genie_wgt_AGKYxF1pi[4]        , cvweight);
    h->FillVertErrorBand("GENIE_AhtBY"             ,var, truth_genie_wgt_AhtBY[2]            , truth_genie_wgt_AhtBY[4]            , cvweight);
    h->FillVertErrorBand("GENIE_BhtBY"             ,var, truth_genie_wgt_BhtBY[2]            , truth_genie_wgt_BhtBY[4]            , cvweight);
    h->FillVertErrorBand("GENIE_CCQEPauliSupViaKF" ,var, truth_genie_wgt_CCQEPauliSupViaKF[2], truth_genie_wgt_CCQEPauliSupViaKF[4], cvweight);
    h->FillVertErrorBand("GENIE_CV1uBY"            ,var, truth_genie_wgt_CV1uBY[2]           , truth_genie_wgt_CV1uBY[4]           , cvweight);
    h->FillVertErrorBand("GENIE_CV2uBY"            ,var, truth_genie_wgt_CV2uBY[2]           , truth_genie_wgt_CV2uBY[4]           , cvweight);
    h->FillVertErrorBand("GENIE_EtaNCEL"           ,var, truth_genie_wgt_EtaNCEL[2]          , truth_genie_wgt_EtaNCEL[4]          , cvweight);
    h->FillVertErrorBand("GENIE_FrAbs_N"           ,var, truth_genie_wgt_FrAbs_N[2]          , truth_genie_wgt_FrAbs_N[4]          , cvweight);
    h->FillVertErrorBand("GENIE_FrAbs_pi"          ,var, truth_genie_wgt_FrAbs_pi[2]         , truth_genie_wgt_FrAbs_pi[4]         , cvweight);
    h->FillVertErrorBand("GENIE_FrCEx_N"           ,var, truth_genie_wgt_FrCEx_N[2]          , truth_genie_wgt_FrCEx_N[4]          , cvweight);
    h->FillVertErrorBand("GENIE_FrCEx_pi"          ,var, truth_genie_wgt_FrCEx_pi[2]         , truth_genie_wgt_FrCEx_pi[4]         , cvweight);
    h->FillVertErrorBand("GENIE_FrElas_N"          ,var, truth_genie_wgt_FrElas_N[2]         , truth_genie_wgt_FrElas_N[4]         , cvweight);
    h->FillVertErrorBand("GENIE_FrElas_pi"         ,var, truth_genie_wgt_FrElas_pi[2]        , truth_genie_wgt_FrElas_pi[4]        , cvweight);
    h->FillVertErrorBand("GENIE_FrInel_N"          ,var, truth_genie_wgt_FrInel_N[2]         , truth_genie_wgt_FrInel_N[4]         , cvweight);
    h->FillVertErrorBand("GENIE_FrInel_pi"         ,var, truth_genie_wgt_FrInel_pi[2]        , truth_genie_wgt_FrInel_pi[4]        , cvweight);
    h->FillVertErrorBand("GENIE_FrPiProd_N"        ,var, truth_genie_wgt_FrPiProd_N[2]       , truth_genie_wgt_FrPiProd_N[4]       , cvweight);
    h->FillVertErrorBand("GENIE_FrPiProd_pi"       ,var, truth_genie_wgt_FrPiProd_pi[2]      , truth_genie_wgt_FrPiProd_pi[4]      , cvweight);
    h->FillVertErrorBand("GENIE_MFP_N"             ,var, truth_genie_wgt_MFP_N[2]            , truth_genie_wgt_MFP_N[4]            , cvweight);
    h->FillVertErrorBand("GENIE_MFP_pi"            ,var, truth_genie_wgt_MFP_pi[2]           , truth_genie_wgt_MFP_pi[4]           , cvweight);
    h->FillVertErrorBand("GENIE_MaCCQE"            ,var, truth_genie_wgt_MaCCQE[2]           , truth_genie_wgt_MaCCQE[4]           , cvweight);
    //h->FillVertErrorBand("GENIE_MaCCQEshape"       ,var, truth_genie_wgt_MaCCQEshape[2]      , truth_genie_wgt_MaCCQEshape[4]      , cvweight);
    h->FillVertErrorBand("GENIE_MaNCEL"            ,var, truth_genie_wgt_MaNCEL[2]           , truth_genie_wgt_MaNCEL[4]           , cvweight);
    h->FillVertErrorBand("GENIE_MaRES"             ,var, truth_genie_wgt_MaRES[2]            , truth_genie_wgt_MaRES[4]            , cvweight);
    h->FillVertErrorBand("GENIE_MvRES"             ,var, truth_genie_wgt_MvRES[2]            , truth_genie_wgt_MvRES[4]            , cvweight);
    //h->FillVertErrorBand("GENIE_NormCCQE"          ,var, truth_genie_wgt_NormCCQE[2]         , truth_genie_wgt_NormCCQE[4]         , cvweight);
    //h->FillVertErrorBand("GENIE_NormCCRES"         ,var, truth_genie_wgt_NormCCRES[2]        , truth_genie_wgt_NormCCRES[4]        , cvweight);
    h->FillVertErrorBand("GENIE_NormDISCC"         ,var, truth_genie_wgt_NormDISCC[2]        , truth_genie_wgt_NormDISCC[4]        , cvweight);
    h->FillVertErrorBand("GENIE_NormNCRES"         ,var, truth_genie_wgt_NormNCRES[2]        , truth_genie_wgt_NormNCRES[4]        , cvweight);
    h->FillVertErrorBand("GENIE_RDecBR1gamma"      ,var, truth_genie_wgt_RDecBR1gamma[2]     , truth_genie_wgt_RDecBR1gamma[4]     , cvweight);
    h->FillVertErrorBand("GENIE_Rvn1pi"            ,var, truth_genie_wgt_Rvn1pi[2]           , truth_genie_wgt_Rvn1pi[4]           , cvweight);
    h->FillVertErrorBand("GENIE_Rvn2pi"            ,var, truth_genie_wgt_Rvn2pi[2]           , truth_genie_wgt_Rvn2pi[4]           , cvweight);
    h->FillVertErrorBand("GENIE_Rvp1pi"            ,var, truth_genie_wgt_Rvp1pi[2]           , truth_genie_wgt_Rvp1pi[4]           , cvweight);
    h->FillVertErrorBand("GENIE_Rvp2pi"            ,var, truth_genie_wgt_Rvp2pi[2]           , truth_genie_wgt_Rvp2pi[4]           , cvweight);
    h->FillVertErrorBand("GENIE_Theta_Delta2Npi"   ,var, truth_genie_wgt_Theta_Delta2Npi[2]  , truth_genie_wgt_Theta_Delta2Npi[4]  , cvweight);
    h->FillVertErrorBand("GENIE_VecFFCCQEshape"    ,var, truth_genie_wgt_VecFFCCQEshape[2]   , truth_genie_wgt_VecFFCCQEshape[4]   , cvweight);
}

void CCProtonPi0_TruthAnalyzer::FillVertErrorBand_Genie_ByHand(MnvH1D* h, double var)
{
    FillVertErrorBand_ByHand(h, var, "GENIE_AGKYxF1pi"         , truth_genie_wgt_AGKYxF1pi[2]        , truth_genie_wgt_AGKYxF1pi[4]        );
    FillVertErrorBand_ByHand(h, var, "GENIE_AhtBY"             , truth_genie_wgt_AhtBY[2]            , truth_genie_wgt_AhtBY[4]            );
    FillVertErrorBand_ByHand(h, var, "GENIE_BhtBY"             , truth_genie_wgt_BhtBY[2]            , truth_genie_wgt_BhtBY[4]            );
    FillVertErrorBand_ByHand(h, var, "GENIE_CCQEPauliSupViaKF" , truth_genie_wgt_CCQEPauliSupViaKF[2], truth_genie_wgt_CCQEPauliSupViaKF[4]);
    FillVertErrorBand_ByHand(h, var, "GENIE_CV1uBY"            , truth_genie_wgt_CV1uBY[2]           , truth_genie_wgt_CV1uBY[4]           );
    FillVertErrorBand_ByHand(h, var, "GENIE_CV2uBY"            , truth_genie_wgt_CV2uBY[2]           , truth_genie_wgt_CV2uBY[4]           );
    FillVertErrorBand_ByHand(h, var, "GENIE_EtaNCEL"           , truth_genie_wgt_EtaNCEL[2]          , truth_genie_wgt_EtaNCEL[4]          );
    FillVertErrorBand_ByHand(h, var, "GENIE_FrAbs_N"           , truth_genie_wgt_FrAbs_N[2]          , truth_genie_wgt_FrAbs_N[4]          );
    FillVertErrorBand_ByHand(h, var, "GENIE_FrAbs_pi"          , truth_genie_wgt_FrAbs_pi[2]         , truth_genie_wgt_FrAbs_pi[4]         );
    FillVertErrorBand_ByHand(h, var, "GENIE_FrCEx_N"           , truth_genie_wgt_FrCEx_N[2]          , truth_genie_wgt_FrCEx_N[4]          );
    FillVertErrorBand_ByHand(h, var, "GENIE_FrCEx_pi"          , truth_genie_wgt_FrCEx_pi[2]         , truth_genie_wgt_FrCEx_pi[4]         );
    FillVertErrorBand_ByHand(h, var, "GENIE_FrElas_N"          , truth_genie_wgt_FrElas_N[2]         , truth_genie_wgt_FrElas_N[4]         );
    FillVertErrorBand_ByHand(h, var, "GENIE_FrElas_pi"         , truth_genie_wgt_FrElas_pi[2]        , truth_genie_wgt_FrElas_pi[4]        );
    FillVertErrorBand_ByHand(h, var, "GENIE_FrInel_N"          , truth_genie_wgt_FrInel_N[2]         , truth_genie_wgt_FrInel_N[4]         );
    FillVertErrorBand_ByHand(h, var, "GENIE_FrInel_pi"         , truth_genie_wgt_FrInel_pi[2]        , truth_genie_wgt_FrInel_pi[4]        );
    FillVertErrorBand_ByHand(h, var, "GENIE_FrPiProd_N"        , truth_genie_wgt_FrPiProd_N[2]       , truth_genie_wgt_FrPiProd_N[4]       );
    FillVertErrorBand_ByHand(h, var, "GENIE_FrPiProd_pi"       , truth_genie_wgt_FrPiProd_pi[2]      , truth_genie_wgt_FrPiProd_pi[4]      );
    FillVertErrorBand_ByHand(h, var, "GENIE_MFP_N"             , truth_genie_wgt_MFP_N[2]            , truth_genie_wgt_MFP_N[4]            );
    FillVertErrorBand_ByHand(h, var, "GENIE_MFP_pi"            , truth_genie_wgt_MFP_pi[2]           , truth_genie_wgt_MFP_pi[4]           );
    FillVertErrorBand_ByHand(h, var, "GENIE_MaCCQE"            , truth_genie_wgt_MaCCQE[2]           , truth_genie_wgt_MaCCQE[4]           );
    //FillVertErrorBand_ByHand(h, var, "GENIE_MaCCQEshape"       , truth_genie_wgt_MaCCQEshape[2]      , truth_genie_wgt_MaCCQEshape[4]      );
    FillVertErrorBand_ByHand(h, var, "GENIE_MaNCEL"            , truth_genie_wgt_MaNCEL[2]           , truth_genie_wgt_MaNCEL[4]           );
    FillVertErrorBand_ByHand(h, var, "GENIE_MaRES"             , truth_genie_wgt_MaRES[2]            , truth_genie_wgt_MaRES[4]            );
    FillVertErrorBand_ByHand(h, var, "GENIE_MvRES"             , truth_genie_wgt_MvRES[2]            , truth_genie_wgt_MvRES[4]            );
    //FillVertErrorBand_ByHand(h, var, "GENIE_NormCCQE"          , truth_genie_wgt_NormCCQE[2]         , truth_genie_wgt_NormCCQE[4]         );
    //FillVertErrorBand_ByHand(h, var, "GENIE_NormCCRES"         , truth_genie_wgt_NormCCRES[2]        , truth_genie_wgt_NormCCRES[4]        );
    FillVertErrorBand_ByHand(h, var, "GENIE_NormDISCC"         , truth_genie_wgt_NormDISCC[2]        , truth_genie_wgt_NormDISCC[4]        );
    FillVertErrorBand_ByHand(h, var, "GENIE_NormNCRES"         , truth_genie_wgt_NormNCRES[2]        , truth_genie_wgt_NormNCRES[4]        );
    FillVertErrorBand_ByHand(h, var, "GENIE_RDecBR1gamma"      , truth_genie_wgt_RDecBR1gamma[2]     , truth_genie_wgt_RDecBR1gamma[4]     );
    FillVertErrorBand_ByHand(h, var, "GENIE_Rvn1pi"            , truth_genie_wgt_Rvn1pi[2]           , truth_genie_wgt_Rvn1pi[4]           );
    FillVertErrorBand_ByHand(h, var, "GENIE_Rvn2pi"            , truth_genie_wgt_Rvn2pi[2]           , truth_genie_wgt_Rvn2pi[4]           );
    FillVertErrorBand_ByHand(h, var, "GENIE_Rvp1pi"            , truth_genie_wgt_Rvp1pi[2]           , truth_genie_wgt_Rvp1pi[4]           );
    FillVertErrorBand_ByHand(h, var, "GENIE_Rvp2pi"            , truth_genie_wgt_Rvp2pi[2]           , truth_genie_wgt_Rvp2pi[4]           );
    FillVertErrorBand_ByHand(h, var, "GENIE_Theta_Delta2Npi"   , truth_genie_wgt_Theta_Delta2Npi[2]  , truth_genie_wgt_Theta_Delta2Npi[4]  );
    FillVertErrorBand_ByHand(h, var, "GENIE_VecFFCCQEshape"    , truth_genie_wgt_VecFFCCQEshape[2]   , truth_genie_wgt_VecFFCCQEshape[4]   );
}

void CCProtonPi0_TruthAnalyzer::FillVertErrorBand_ByHand(MnvH1D* h, double var, std::string error_name, double err_down, double err_up)
{
    std::vector<double> errors;
    errors.push_back(err_down);
    errors.push_back(err_up);

    FillVertErrorBand_ByHand(h, var, error_name, errors);
}

void CCProtonPi0_TruthAnalyzer::FillVertErrorBand_ByHand(MnvH1D* h, double var, std::string error_name, std::vector<double> errors)
{
    // Get a Pointer to Error Band
    MnvVertErrorBand* err_band =  h->GetVertErrorBand(error_name);

    // Get a Pointer to Histograms 
    std::vector<TH1D*> err_hists = err_band->GetHists();
   
    // Sanity Check
    if (err_hists.size() != errors.size()) {
        std::cout<<"WARNING! Can not Fill Vertical Error Band: "<<error_name<<std::endl;
        exit(1);
    }

    // Fill Error Band Base Histogram with Default cvweight
    int cvbin = err_band->TH1D::Fill(var, cvweight);

    // Update Errors
    if (cvbin == -1){
        cvbin = err_band->FindBin(var);
    }

    for (unsigned int i = 0; i < err_hists.size(); ++i ) {
        
        // wgt_bckg is universe_wgt / cv_wgt 
        double wgt_bckg = applyBckgConstraints_Unv ? GetBckgConstraint(error_name, i) : 1.0;
      
        const double applyWeight = cvweight * wgt_bckg;
        const double wgtU = errors[i]*applyWeight;
        err_hists[i]->AddBinContent( cvbin, wgtU );

        const double err = err_hists[i]->GetBinError(cvbin);
        const double newerr2 = err*err + wgtU*wgtU; 
        const double newerr = (0.<newerr2) ? sqrt(newerr2) : 0.;
        err_hists[i]->SetBinError( cvbin, newerr );
    }
}

double CCProtonPi0_TruthAnalyzer::GetBckgConstraint(std::string error_name, int hist_ind)
{
    // Find the Bckg Constraint if it is one of the constrained events
    if (truth_isBckg_Compact_SinglePiPlus){
        double bckg_wgt_SinglePiPlus = BckgConstrainer.GetBckgConstraint(error_name, hist_ind, "SinglePiPlus");
        return bckg_wgt_SinglePiPlus / cv_wgt_SinglePiPlus;
    }else if (truth_isBckg_Compact_QELike){
        double bckg_wgt_QELike = BckgConstrainer.GetBckgConstraint(error_name, hist_ind, "QELike");
        return bckg_wgt_QELike / cv_wgt_QELike;
    }else if (truth_isBckg_Compact_WithPi0){
        double bckg_wgt_WithPi0 = BckgConstrainer.GetBckgConstraint(error_name, hist_ind, "WithPi0");
        return bckg_wgt_WithPi0 / cv_wgt_WithPi0;
    }else{
        return 1.0;
    }
}

void CCProtonPi0_TruthAnalyzer::initCVWeights()
{
    cv_wgt_SinglePiPlus = BckgConstrainer.GetBckgConstraint("CentralValue", 0, "SinglePiPlus");
    cv_wgt_QELike = BckgConstrainer.GetBckgConstraint("CentralValue", 0, "QELike");
    cv_wgt_WithPi0 = BckgConstrainer.GetBckgConstraint("CentralValue", 0, "WithPi0");

    cv_err_SinglePiPlus = BckgConstrainer.GetBckgConstraintErr("CentralValue", 0, "SinglePiPlus");
    cv_err_QELike = BckgConstrainer.GetBckgConstraintErr("CentralValue", 0, "QELike");
    cv_err_WithPi0 = BckgConstrainer.GetBckgConstraintErr("CentralValue", 0, "WithPi0");

    std::cout<<"CV Background Weights = "<<cv_wgt_SinglePiPlus<<" "<<cv_wgt_QELike<<" "<<cv_wgt_WithPi0<<std::endl;
    std::cout<<"CV Background Weight Errors = "<<cv_err_SinglePiPlus<<" "<<cv_err_QELike<<" "<<cv_err_WithPi0<<std::endl;
}

void CCProtonPi0_TruthAnalyzer::Calc_WeightFromSystematics()
{
    UpdateFluxReweighter(mc_run); 
        
    // Replace cvweight with Flux Weight
    cvweight = GetFluxWeight(mc_incomingE * MeV_to_GeV, mc_incoming);

    // Apply Background Constraints
    if (applyBckgConstraints_CV){
        if (truth_isBckg_Compact_SinglePiPlus) cvweight *= cv_wgt_SinglePiPlus;
        else if (truth_isBckg_Compact_QELike) cvweight *= cv_wgt_QELike;
        else if (truth_isBckg_Compact_WithPi0) cvweight *= cv_wgt_WithPi0;
    }

}
void CCProtonPi0_TruthAnalyzer::AddOtherErrorBands_FillWithCV()
{
    // Add Vertical Error Bands
    AddVertErrorBands_TruthTree(muon_P_mc_truth_all_signal);
    AddVertErrorBands_TruthTree(muon_theta_mc_truth_all_signal);
    AddVertErrorBands_TruthTree(pi0_P_mc_truth_all_signal);
    AddVertErrorBands_TruthTree(pi0_KE_mc_truth_all_signal);
    AddVertErrorBands_TruthTree(pi0_theta_mc_truth_all_signal);
    AddVertErrorBands_TruthTree(Enu_mc_truth_all_signal);
    AddVertErrorBands_TruthTree(QSq_mc_truth_all_signal);
    AddVertErrorBands_TruthTree(W_mc_truth_all_signal);

    // Add Lateral Error Bands
    AddLatErrorBands_TruthTree(muon_P_mc_truth_all_signal);
    AddLatErrorBands_TruthTree(muon_theta_mc_truth_all_signal);
    AddLatErrorBands_TruthTree(pi0_P_mc_truth_all_signal);
    AddLatErrorBands_TruthTree(pi0_KE_mc_truth_all_signal);
    AddLatErrorBands_TruthTree(pi0_theta_mc_truth_all_signal);
    AddLatErrorBands_TruthTree(Enu_mc_truth_all_signal);
    AddLatErrorBands_TruthTree(QSq_mc_truth_all_signal);
    AddLatErrorBands_TruthTree(W_mc_truth_all_signal);
}

void CCProtonPi0_TruthAnalyzer::FillSignalHistograms()
{
    Calc_WeightFromSystematics();

    // Cross Section Variables
    FillHistogram(muon_P_mc_truth_all_signal, truth_muon_P * MeV_to_GeV);
    FillHistogram(muon_theta_mc_truth_all_signal, truth_muon_theta * rad_to_deg);
    FillHistogram(pi0_P_mc_truth_all_signal, truth_pi0_P * MeV_to_GeV);
    FillHistogram(pi0_KE_mc_truth_all_signal, truth_pi0_KE * MeV_to_GeV);
    FillHistogram(pi0_theta_mc_truth_all_signal, truth_pi0_theta * rad_to_deg);
    FillHistogram(Enu_mc_truth_all_signal, mc_incomingE * MeV_to_GeV);
    FillHistogram(QSq_mc_truth_all_signal, truth_QSq_exp * MeVSq_to_GeVSq);
    FillHistogram(W_mc_truth_all_signal, truth_W_exp * MeV_to_GeV);
    
    // Signal Characteristics
    if (mc_intType == 1){
        nQE.increment();
        FillHistogram(mc_Q2_QE, truth_QSq_exp * MeVSq_to_GeVSq);
        FillHistogram(mc_incomingE_QE, mc_incomingE * MeV_to_GeV);
        FillHistogram(mc_w_QE, truth_W_exp * MeV_to_GeV);
    }else if (mc_intType == 2){
        if (mc_resID == 0){ 
            nRES_1232.increment(); 
            FillHistogram(mc_Q2_RES_1232, truth_QSq_exp * MeVSq_to_GeVSq);
            FillHistogram(mc_incomingE_RES_1232, mc_incomingE * MeV_to_GeV);
            FillHistogram(mc_w_RES_1232, truth_W_exp * MeV_to_GeV);
        }else if (mc_resID == 1){
            nRES_1535.increment(); 
            FillHistogram(mc_Q2_RES_1535, truth_QSq_exp * MeVSq_to_GeVSq);
            FillHistogram(mc_incomingE_RES_1535, mc_incomingE * MeV_to_GeV);
            FillHistogram(mc_w_RES_1535, truth_W_exp * MeV_to_GeV);
        }else if (mc_resID == 2){
            nRES_1520.increment(); 
            FillHistogram(mc_Q2_RES_1520, truth_QSq_exp * MeVSq_to_GeVSq);
            FillHistogram(mc_incomingE_RES_1520, mc_incomingE * MeV_to_GeV);
            FillHistogram(mc_w_RES_1520, truth_W_exp * MeV_to_GeV);
        }else{
            nRES_Other.increment(); 
            FillHistogram(mc_Q2_RES_Other, truth_QSq_exp * MeVSq_to_GeVSq);
            FillHistogram(mc_incomingE_RES_Other, mc_incomingE * MeV_to_GeV);
            FillHistogram(mc_w_RES_Other, truth_W_exp * MeV_to_GeV);
        }
    }else if (mc_intType == 3){
        int nFS_pions = Get_nFS_pions();
        if (mc_w*MeV_to_GeV < 1.7 ){
            nNon_RES.increment();
            FillHistogram(mc_Q2_Non_RES, truth_QSq_exp * MeVSq_to_GeVSq);
            FillHistogram(mc_incomingE_Non_RES, mc_incomingE * MeV_to_GeV);
            FillHistogram(mc_w_Non_RES, truth_W_exp * MeV_to_GeV);
        }else if (nFS_pions == 1){
            nDIS_1_pi.increment();
            FillHistogram(mc_Q2_DIS_1_pi, truth_QSq_exp * MeVSq_to_GeVSq);
            FillHistogram(mc_incomingE_DIS_1_pi, mc_incomingE * MeV_to_GeV);
            FillHistogram(mc_w_DIS_1_pi, truth_W_exp * MeV_to_GeV);
        }else if (nFS_pions == 2){
            nDIS_2_pi.increment();
            FillHistogram(mc_Q2_DIS_2_pi, truth_QSq_exp * MeVSq_to_GeVSq);
            FillHistogram(mc_incomingE_DIS_2_pi, mc_incomingE * MeV_to_GeV);
            FillHistogram(mc_w_DIS_2_pi, truth_W_exp * MeV_to_GeV);
        }else if (nFS_pions > 2){
            nDIS_Multi_pi.increment();
            FillHistogram(mc_Q2_DIS_Multi_pi, truth_QSq_exp * MeVSq_to_GeVSq);
            FillHistogram(mc_incomingE_DIS_Multi_pi, mc_incomingE * MeV_to_GeV);
            FillHistogram(mc_w_DIS_Multi_pi, truth_W_exp * MeV_to_GeV);
        }
    }else{
        std::cout<<"WARNING! Signal Event with different interaction Type!"<<std::endl;
    }
}

int CCProtonPi0_TruthAnalyzer::Get_nFS_pions()
{
    int nFS_pions = 0;

    for (int i = 0; i < mc_er_nPart; ++i){
        if (std::abs(mc_er_ID[i]) == 211 || mc_er_ID[i] == 111){
            if (isMother_DIS_Fragment(i)) nFS_pions++;
        }
    }

    return nFS_pions;
}

bool CCProtonPi0_TruthAnalyzer::isMother_DIS_Fragment(int ind)
{
    int mother_ind = mc_er_mother[ind];

    if (mc_er_status[mother_ind] == 12) return true;
    else return false;
}

void CCProtonPi0_TruthAnalyzer::PrintEventRecord()
{
    std::cout<<"-------------"<<std::endl;
    for (int i = 0; i < mc_er_nPart; ++i){
        std::cout<<mc_er_ID[i]<<" "<<mc_er_status[i]<<" "<<mc_er_mother[i]<<std::endl;
    }
}

void CCProtonPi0_TruthAnalyzer::FillMichelHistograms()
{
    if (truth_isBckg_withMichel && truth_michelPion_begin_dist_vtx < 25.0){
        double dist_x = mc_vtx[0]-truth_michelMuon_endPoint[0];
        double dist_y = mc_vtx[1]-truth_michelMuon_endPoint[1];
        double dist_z = mc_vtx[2]-truth_michelMuon_endPoint[2];
        michel_ZX->Fill(dist_z, dist_x);
        michel_ZY->Fill(dist_z, dist_y);
        michel_dist_X->Fill(dist_x);
        michel_dist_Y->Fill(dist_y);
        michel_dist_Z->Fill(dist_z);
        michel_dist_total->Fill(truth_michelMuon_end_dist_vtx);
        michel_pionP->Fill(truth_michelPion_P);
        michel_pion_dist->Fill(truth_michelPion_begin_dist_vtx);
    }
}

void CCProtonPi0_TruthAnalyzer::writeHistograms()
{
    f->cd();

    // Cross Section Variables
    muon_P_mc_truth_all_signal->Write();
    muon_theta_mc_truth_all_signal->Write();
    pi0_P_mc_truth_all_signal->Write();
    pi0_KE_mc_truth_all_signal->Write();
    pi0_theta_mc_truth_all_signal->Write();
    Enu_mc_truth_all_signal->Write();
    QSq_mc_truth_all_signal->Write();
    W_mc_truth_all_signal->Write();

    // Michel Electron
    michel_ZX->Write();
    michel_ZY->Write();
    michel_dist_X->Write();
    michel_dist_Y->Write();
    michel_dist_Z->Write();
    michel_dist_total->Write();
    michel_pionP->Write();
    michel_pion_dist->Write();

    // Signal Q2
    mc_Q2_QE->Write();
    mc_Q2_RES_1232->Write();
    mc_Q2_RES_1535->Write();
    mc_Q2_RES_1520->Write();
    mc_Q2_RES_Other->Write();
   
    mc_Q2_DIS_1_pi->Write();
    mc_Q2_DIS_2_pi->Write();
    mc_Q2_DIS_Multi_pi->Write();
    mc_Q2_Non_RES->Write();
 
    // Signal incomingE
    mc_incomingE_QE->Write();
    mc_incomingE_RES_1232->Write();
    mc_incomingE_RES_1535->Write();
    mc_incomingE_RES_1520->Write();
    mc_incomingE_RES_Other->Write();
   
    mc_incomingE_DIS_1_pi->Write();
    mc_incomingE_DIS_2_pi->Write();
    mc_incomingE_DIS_Multi_pi->Write();
    mc_incomingE_Non_RES->Write();
 
    // Signal w
    mc_w_QE->Write();
    mc_w_RES_1232->Write();
    mc_w_RES_1535->Write();
    mc_w_RES_1520->Write();
    mc_w_RES_Other->Write();
   
    mc_w_DIS_1_pi->Write();
    mc_w_DIS_2_pi->Write();
    mc_w_DIS_Multi_pi->Write();
    mc_w_Non_RES->Write();
  
    
    f->Close();
}

#endif


