#ifndef CCProtonPi0_TruthAnalyzer_cpp
#define CCProtonPi0_TruthAnalyzer_cpp

#include "CCProtonPi0_TruthAnalyzer.h"

using namespace PlotUtils;

void CCProtonPi0_TruthAnalyzer::Loop(std::string playlist)
{
    // Control Flow
    bool applyMaxEvents = false;
    double nMaxEvents = 1000000;

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
            
            Calc_WeightFromSystematics();
            
            // Count Signal Out due to Acceptance or Kinematics
            if (truth_isSignalOut_Acceptance) nSignalOut_Acceptance.increment();
            if (truth_isSignalOut_Kinematics) nSignalOut_Kinematics.increment();

            // Count Signal and Background inside Fiducial Volume
            if (truth_isSignal) nSignal.increment();
            else nBckg.increment();

            if (truth_isSignal){
                // Detailed Interaction breakout
                FillSignal_InteractionType(); // Also counts the N(Signal) in each type
               
                // ------------------------------------------------------------
                // Fill XSec Variables for different kinds
                // ------------------------------------------------------------
                // Default Signal Histograms
                FillSignal_XSec_Variables();
                
                // FSI Type
                int FSIType = GetFSIType();
                CountFSIType(FSIType);
                FillSignal_XSec_Variables_FSIType(FSIType);
               
                // Int Type
                int IntType = GetIntType();
                FillSignal_XSec_Variables_IntType(IntType);

                //FillSignal_Test();
            }

            if (truth_isSignal_BeforeFSI){
                nSignal_BeforeFSI.increment();
                FillSignal_XSec_Variables_BeforeFSI();
            }
        }
    }

    // Add Other Error Bands and Fill With CV
    AddOtherErrorBands_FillWithCV();

    writeTextFile();
    writeHistograms();
}

void CCProtonPi0_TruthAnalyzer::CountFSIType(int type)
{
    switch (type) {
        case 0: nFSI_NonInteracting.increment();
                break;
        case 1: nFSI_Elastic.increment();
                break;
        case 2: nFSI_Inelastic.increment();
                break;
        case 3: nFSI_ChargeExchange.increment();
                break;
        case 4: nFSI_MultiPi.increment();
                break;
        case 5: nFSI_NucleonToPi.increment();
                break;
        default: RunTimeError("Wrong FSI Type");
    }
}

int CCProtonPi0_TruthAnalyzer::GetIntType()
{
    int IntType = -1;
    
    if (mc_intType == 2){
        if (mc_resID == 0){ 
            // Delta Resonance
            IntType = 0;
        }else{
            // Other Resonance
            IntType = 1;
        }
    }else{
        // Non-Resonance 
        IntType = 2;
    }

    return IntType;
}

int CCProtonPi0_TruthAnalyzer::GetFSIType()
{
    const double eps = 0.99;
    //-------------------------------------------------------------------------
    // FSI Based on NPions
    //-------------------------------------------------------------------------
    int nPi0 = 0;
    int nChargedPi = 0;

    // Count Number of Pions
    for (int i = 0; i < mc_er_nPart; ++i){
        if (mc_er_status[i] != 14) continue;
        if (mc_er_ID[i] == 111) nPi0++;    
        else if (std::abs(mc_er_ID[i]) == 211 ) nChargedPi++;
    }

    if ((nChargedPi + nPi0) > 1){
        return 4;
    }else if (nPi0 == 0 && nChargedPi == 1){
        return 3;
    }else if (nPi0 == 0){
        return 5;
    }

    //-------------------------------------------------------------------------
    // FSI Based on Pion Scattering 
    //-------------------------------------------------------------------------
    if (nPi0 == 1){
        double E[3];
        double KE[3];
        double Px[3];
        double Py[3];
        double Pz[3];
        // Loop again to get Pi0 Initial and Final Energy
        for (int i = 0; i < mc_er_nPart; ++i){
            if (mc_er_ID[i] != 111) continue;
            
            if (mc_er_status[i] == 14){
                E[0] = mc_er_E[i];
                KE[0] = E[0]-pi0_mass;
                Px[0] = mc_er_Px[i];
                Py[0] = mc_er_Py[i];
                Pz[0] = mc_er_Pz[i];
            }else if (mc_er_status[i] == 1){
                E[1] = mc_er_E[i];
                KE[1] = E[1]-pi0_mass;
                Px[1] = mc_er_Px[i];
                Py[1] = mc_er_Py[i];
                Pz[1] = mc_er_Pz[i];
            }
        }

        // Calculate the Initial and Final Difference
        E[2] = std::abs(E[0]-E[1]);
        KE[2] = std::abs(KE[0]-KE[1]);
        Px[2] = std::abs(Px[0]-Px[1]);
        Py[2] = std::abs(Py[0]-Py[1]);
        Pz[2] = std::abs(Pz[0]-Pz[1]);

        bool isKE_Changed = KE[2] > eps;
        bool isDirection_Changed = (Px[2] > eps || Py[2] > eps || Pz[2] > eps);

        if (isKE_Changed){
            // Kinetic Energy is Changed
            return 2;
        }else if (isDirection_Changed){
            // KE Conserved, direction changed
            return 1;
        }else{
            // Nothing Changed
            return 0;
        }
    } 
    
    std::cout<<"WARNING! In this stage there must be single Pi0"<<std::endl;
    return -1;
    
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
    WriteCounter(nSignal_BeforeFSI, nFidVol);
    WriteCounter(nSignalOut_Acceptance, nFidVol);
    WriteCounter(nSignalOut_Kinematics, nFidVol);
    WriteCounter(nBckg, nFidVol);
    textFile<<std::endl;

    // FSI Types
    WriteCounter(nFSI_NonInteracting, nSignal);
    WriteCounter(nFSI_Elastic, nSignal);
    WriteCounter(nFSI_Inelastic, nSignal);
    WriteCounter(nFSI_ChargeExchange, nSignal);
    WriteCounter(nFSI_MultiPi, nSignal);
    WriteCounter(nFSI_NucleonToPi, nSignal);
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
    rootDir = Folder_List::rootDir_Truth_mc_Test;

    std::cout<<"\tRoot File: "<<rootDir<<std::endl;
    f = new TFile(rootDir.c_str(),"RECREATE");
 
    initHistograms();
    
    openTextFiles();

    resetCounters();
    
    std::cout<<"Finished"<<std::endl;
}

void CCProtonPi0_TruthAnalyzer::resetCounters() 
{
    nAll.setName("nAll");
   
    nFidVol.setName("nFidVol");
    nFidVol_Out.setName("nFidVol_Out");

    nSignal.setName("nSignal");
    nSignal_BeforeFSI.setName("nSignal_BeforeFSI");
    nSignalOut_Acceptance.setName("nSignalOut_Acceptance");
    nSignalOut_Kinematics.setName("nSignalOut_Kinematics");
    nBckg.setName("nBckg");

    // FSI 
    nFSI_NonInteracting.setName("nFSI_NonInteracting");
    nFSI_Elastic.setName("nFSI_Elastic");
    nFSI_Inelastic.setName("nFSI_Inelastic");
    nFSI_ChargeExchange.setName("nFSI_ChargeExchange");
    nFSI_MultiPi.setName("nFSI_MultiPi");
    nFSI_NucleonToPi.setName("nFSI_NucleonToPi");

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
    OpenTextFile("TruthInfo.txt", textFile);
    OpenTextFile("EventRecord.txt", logFile);

    logFile<<std::left;
    logFile.width(8); logFile<<"PDG"<<" "; 
    logFile.width(8); logFile<<"Status"<<" "; 
    logFile.width(8); logFile<<"Mother"<<" ";    
    logFile.width(12); logFile<<"ER E"<<" ";      
    logFile.width(12); logFile<<"Traj E"<<" ";      
    logFile<<endl;

}

void CCProtonPi0_TruthAnalyzer::initHistograms()
{
    // ------------------------------------------------------------------------
    // Test Variables
    // ------------------------------------------------------------------------
    Test_pi0_P = new TH1D("Test_pi0_P","Pi0 Momentum Check",40,-100,100);
    Test_pi0_P->GetXaxis()->SetTitle("P_{Trajectory} - P_{EventRecord} [Mev]");
    Test_pi0_P->GetYaxis()->SetTitle("NEvents()");

    Test_pi0_theta = new TH1D("Test_pi0_theta","Pi0 Angle Check",40,-10,10);
    Test_pi0_theta->GetXaxis()->SetTitle("#theta_{Trajectory} - #theta_{EventRecord} [degree]");
    Test_pi0_theta->GetYaxis()->SetTitle("NEvents()");

    // ------------------------------------------------------------------------
    // Muon Variables
    // ------------------------------------------------------------------------
    muon_P_mc_truth_all_signal = new MnvH1D( "muon_P_mc_truth_all_signal","Muon Momentum for Signal Events",binList.size_muon_P, binList.a_muon_P);
    muon_P_mc_truth_all_signal->GetXaxis()->SetTitle("Momentum [GeV]");
    muon_P_mc_truth_all_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBand_Flux(muon_P_mc_truth_all_signal);
    AddVertErrorBand_Genie(muon_P_mc_truth_all_signal);

    muon_theta_mc_truth_all_signal = new MnvH1D( "muon_theta_mc_truth_all_signal","Muon Theta for Signal Events",binList.size_muon_theta, binList.a_muon_theta);
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

    pi0_theta_mc_truth_all_signal = new MnvH1D( "pi0_theta_mc_truth_all_signal","Pi0 Theta for Signal Events",binList.size_pi0_theta, binList.a_pi0_theta);
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
    // Muon Variables -- Before FSI
    // ------------------------------------------------------------------------
    muon_P_mc_truth_all_signal_BeforeFSI = new MnvH1D( "muon_P_mc_truth_all_signal_BeforeFSI","Muon Momentum for Signal Events",binList.size_muon_P, binList.a_muon_P);
    muon_P_mc_truth_all_signal_BeforeFSI->GetXaxis()->SetTitle("Momentum [GeV]");
    muon_P_mc_truth_all_signal_BeforeFSI->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBand_Flux(muon_P_mc_truth_all_signal_BeforeFSI);
    AddVertErrorBand_Genie(muon_P_mc_truth_all_signal_BeforeFSI);

    muon_theta_mc_truth_all_signal_BeforeFSI = new MnvH1D( "muon_theta_mc_truth_all_signal_BeforeFSI","Muon Theta for Signal Events",binList.size_muon_theta, binList.a_muon_theta);
    muon_theta_mc_truth_all_signal_BeforeFSI->GetXaxis()->SetTitle("Theta");
    muon_theta_mc_truth_all_signal_BeforeFSI->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBand_Flux(muon_theta_mc_truth_all_signal_BeforeFSI);
    AddVertErrorBand_Genie(muon_theta_mc_truth_all_signal_BeforeFSI);

    // ------------------------------------------------------------------------
    // Pi0 Variables -- Before FSI
    // ------------------------------------------------------------------------
    pi0_P_mc_truth_all_signal_BeforeFSI = new MnvH1D( "pi0_P_mc_truth_all_signal_BeforeFSI","Pi0 Momentum for Signal Events",binList.size_pi0_P, binList.a_pi0_P);
    pi0_P_mc_truth_all_signal_BeforeFSI->GetXaxis()->SetTitle("Momentum [GeV]");
    pi0_P_mc_truth_all_signal_BeforeFSI->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBand_Flux(pi0_P_mc_truth_all_signal_BeforeFSI);
    AddVertErrorBand_Genie(pi0_P_mc_truth_all_signal_BeforeFSI);

    pi0_KE_mc_truth_all_signal_BeforeFSI = new MnvH1D( "pi0_KE_mc_truth_all_signal_BeforeFSI","Pi0 Kinetic Energy for Signal Events",binList.size_pi0_KE, binList.a_pi0_KE);
    pi0_KE_mc_truth_all_signal_BeforeFSI->GetXaxis()->SetTitle("Kinetic Energy [GeV]");
    pi0_KE_mc_truth_all_signal_BeforeFSI->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBand_Flux(pi0_KE_mc_truth_all_signal_BeforeFSI);
    AddVertErrorBand_Genie(pi0_KE_mc_truth_all_signal_BeforeFSI);

    pi0_theta_mc_truth_all_signal_BeforeFSI = new MnvH1D( "pi0_theta_mc_truth_all_signal_BeforeFSI","Pi0 Theta for Signal Events",binList.size_pi0_theta, binList.a_pi0_theta);
    pi0_theta_mc_truth_all_signal_BeforeFSI->GetXaxis()->SetTitle("Theta");
    pi0_theta_mc_truth_all_signal_BeforeFSI->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBand_Flux(pi0_theta_mc_truth_all_signal_BeforeFSI);
    AddVertErrorBand_Genie(pi0_theta_mc_truth_all_signal_BeforeFSI);

    // ------------------------------------------------------------------------
    // Neutrino Energy, Q2 & W -- Before FSI
    // ------------------------------------------------------------------------
    Enu_mc_truth_all_signal_BeforeFSI = new MnvH1D( "Enu_mc_truth_all_signal_BeforeFSI","Neutrino Energy for Signal Events",binList.size_Enu, binList.a_Enu);
    Enu_mc_truth_all_signal_BeforeFSI->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
    Enu_mc_truth_all_signal_BeforeFSI->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBand_Flux(Enu_mc_truth_all_signal_BeforeFSI);
    AddVertErrorBand_Genie(Enu_mc_truth_all_signal_BeforeFSI);

    QSq_mc_truth_all_signal_BeforeFSI = new MnvH1D( "QSq_mc_truth_all_signal_BeforeFSI","Q^{2} for Signal Events",binList.size_QSq, binList.a_QSq);
    QSq_mc_truth_all_signal_BeforeFSI->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    QSq_mc_truth_all_signal_BeforeFSI->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBand_Flux(QSq_mc_truth_all_signal_BeforeFSI);
    AddVertErrorBand_Genie(QSq_mc_truth_all_signal_BeforeFSI);

    W_mc_truth_all_signal_BeforeFSI = new MnvH1D( "W_mc_truth_all_signal_BeforeFSI","W for Signal Events",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max());
    W_mc_truth_all_signal_BeforeFSI->GetXaxis()->SetTitle("W [GeV]");
    W_mc_truth_all_signal_BeforeFSI->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBand_Flux(W_mc_truth_all_signal_BeforeFSI);
    AddVertErrorBand_Genie(W_mc_truth_all_signal_BeforeFSI);

    // ------------------------------------------------------------------------
    // FSI Type 
    // ------------------------------------------------------------------------
    MnvH1D* temp;
    for (int i = 0; i < nFSIType; ++i){
        temp = new MnvH1D( Form("%s_%d","muon_P_mc_truth_all_signal_FSIType",i),"Muon Momentum for Signal Events",binList.size_muon_P, binList.a_muon_P);
        temp->GetXaxis()->SetTitle("Momentum [GeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBand_Flux(temp);
        AddVertErrorBand_Genie(temp);
        muon_P_mc_truth_all_signal_FSIType.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","muon_theta_mc_truth_all_signal_FSIType",i),"Muon Theta for Signal Events",binList.size_muon_theta, binList.a_muon_theta);
        temp->GetXaxis()->SetTitle("Theta");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBand_Flux(temp);
        AddVertErrorBand_Genie(temp);
        muon_theta_mc_truth_all_signal_FSIType.push_back(temp);
        
        temp = new MnvH1D( Form("%s_%d","pi0_P_mc_truth_all_signal_FSIType",i),"Pi0 Momentum for Signal Events",binList.size_pi0_P, binList.a_pi0_P);
        temp->GetXaxis()->SetTitle("Momentum [GeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBand_Flux(temp);
        AddVertErrorBand_Genie(temp);
        pi0_P_mc_truth_all_signal_FSIType.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","pi0_KE_mc_truth_all_signal_FSIType",i),"Pi0 Kinetic Energy for Signal Events",binList.size_pi0_KE, binList.a_pi0_KE);
        temp->GetXaxis()->SetTitle("Kinetic Energy [GeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBand_Flux(temp);
        AddVertErrorBand_Genie(temp);
        pi0_KE_mc_truth_all_signal_FSIType.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","pi0_theta_mc_truth_all_signal_FSIType",i), "Pi0 Theta for Signal Events",binList.size_pi0_theta, binList.a_pi0_theta);
        temp->GetXaxis()->SetTitle("Theta");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBand_Flux(temp);
        AddVertErrorBand_Genie(temp);
        pi0_theta_mc_truth_all_signal_FSIType.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","Enu_mc_truth_all_signal_FSIType",i),"Neutrino Energy for Signal Events",binList.size_Enu, binList.a_Enu);
        temp->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBand_Flux(temp);
        AddVertErrorBand_Genie(temp);
        Enu_mc_truth_all_signal_FSIType.push_back(temp);

        temp = new MnvH1D( Form("%s_%d", "QSq_mc_truth_all_signal_FSIType",i),"Q^{2} for Signal Events",binList.size_QSq, binList.a_QSq);
        temp->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBand_Flux(temp);
        AddVertErrorBand_Genie(temp);
        QSq_mc_truth_all_signal_FSIType.push_back(temp);

        temp = new MnvH1D( Form("%s_%d", "W_mc_truth_all_signal_FSIType",i), "W for Signal Events",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max());
        temp->GetXaxis()->SetTitle("W [GeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBand_Flux(temp);
        AddVertErrorBand_Genie(temp);
        W_mc_truth_all_signal_FSIType.push_back(temp);
    }

    // ------------------------------------------------------------------------
    // Interaction Type 
    // ------------------------------------------------------------------------
    for (int i = 0; i < nIntType; i++){
        temp = new MnvH1D( Form("%s_%d","muon_P_mc_truth_all_signal_IntType",i),"Muon Momentum for Signal Events",binList.size_muon_P, binList.a_muon_P);
        temp->GetXaxis()->SetTitle("Momentum [GeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBand_Flux(temp);
        AddVertErrorBand_Genie(temp);
        muon_P_mc_truth_all_signal_IntType.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","muon_theta_mc_truth_all_signal_IntType",i),"Muon Theta for Signal Events",binList.size_muon_theta, binList.a_muon_theta);
        temp->GetXaxis()->SetTitle("Theta");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBand_Flux(temp);
        AddVertErrorBand_Genie(temp);
        muon_theta_mc_truth_all_signal_IntType.push_back(temp);
        
        temp = new MnvH1D( Form("%s_%d","pi0_P_mc_truth_all_signal_IntType",i),"Pi0 Momentum for Signal Events",binList.size_pi0_P, binList.a_pi0_P);
        temp->GetXaxis()->SetTitle("Momentum [GeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBand_Flux(temp);
        AddVertErrorBand_Genie(temp);
        pi0_P_mc_truth_all_signal_IntType.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","pi0_KE_mc_truth_all_signal_IntType",i),"Pi0 Kinetic Energy for Signal Events",binList.size_pi0_KE, binList.a_pi0_KE);
        temp->GetXaxis()->SetTitle("Kinetic Energy [GeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBand_Flux(temp);
        AddVertErrorBand_Genie(temp);
        pi0_KE_mc_truth_all_signal_IntType.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","pi0_theta_mc_truth_all_signal_IntType",i), "Pi0 Theta for Signal Events",binList.size_pi0_theta, binList.a_pi0_theta);
        temp->GetXaxis()->SetTitle("Theta");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBand_Flux(temp);
        AddVertErrorBand_Genie(temp);
        pi0_theta_mc_truth_all_signal_IntType.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","Enu_mc_truth_all_signal_IntType",i),"Neutrino Energy for Signal Events",binList.size_Enu, binList.a_Enu);
        temp->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBand_Flux(temp);
        AddVertErrorBand_Genie(temp);
        Enu_mc_truth_all_signal_IntType.push_back(temp);

        temp = new MnvH1D( Form("%s_%d", "QSq_mc_truth_all_signal_IntType",i),"Q^{2} for Signal Events",binList.size_QSq, binList.a_QSq);
        temp->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBand_Flux(temp);
        AddVertErrorBand_Genie(temp);
        QSq_mc_truth_all_signal_IntType.push_back(temp);

        temp = new MnvH1D( Form("%s_%d", "W_mc_truth_all_signal_IntType",i), "W for Signal Events",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max());
        temp->GetXaxis()->SetTitle("W [GeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBand_Flux(temp);
        AddVertErrorBand_Genie(temp);
        W_mc_truth_all_signal_IntType.push_back(temp);
    }

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
        
        // wgt_bckg is  1.0 -- there is no Background
        double wgt_bckg = 1.0;
      
        const double applyWeight = cvweight * wgt_bckg;
        const double wgtU = errors[i]*applyWeight;
        err_hists[i]->AddBinContent( cvbin, wgtU );

        const double err = err_hists[i]->GetBinError(cvbin);
        const double newerr2 = err*err + wgtU*wgtU; 
        const double newerr = (0.<newerr2) ? sqrt(newerr2) : 0.;
        err_hists[i]->SetBinError( cvbin, newerr );
    }
}

void CCProtonPi0_TruthAnalyzer::Calc_WeightFromSystematics()
{
    UpdateFluxReweighter(mc_run); 
        
    // Replace cvweight with Flux Weight
    cvweight = GetFluxWeight(mc_incomingE * MeV_to_GeV, mc_incoming);
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

    // Add Vertical Error Bands
    AddVertErrorBands_TruthTree(muon_P_mc_truth_all_signal_BeforeFSI);
    AddVertErrorBands_TruthTree(muon_theta_mc_truth_all_signal_BeforeFSI);
    AddVertErrorBands_TruthTree(pi0_P_mc_truth_all_signal_BeforeFSI);
    AddVertErrorBands_TruthTree(pi0_KE_mc_truth_all_signal_BeforeFSI);
    AddVertErrorBands_TruthTree(pi0_theta_mc_truth_all_signal_BeforeFSI);
    AddVertErrorBands_TruthTree(Enu_mc_truth_all_signal_BeforeFSI);
    AddVertErrorBands_TruthTree(QSq_mc_truth_all_signal_BeforeFSI);
    AddVertErrorBands_TruthTree(W_mc_truth_all_signal_BeforeFSI);

    // Add Lateral Error Bands
    AddLatErrorBands_TruthTree(muon_P_mc_truth_all_signal_BeforeFSI);
    AddLatErrorBands_TruthTree(muon_theta_mc_truth_all_signal_BeforeFSI);
    AddLatErrorBands_TruthTree(pi0_P_mc_truth_all_signal_BeforeFSI);
    AddLatErrorBands_TruthTree(pi0_KE_mc_truth_all_signal_BeforeFSI);
    AddLatErrorBands_TruthTree(pi0_theta_mc_truth_all_signal_BeforeFSI);
    AddLatErrorBands_TruthTree(Enu_mc_truth_all_signal_BeforeFSI);
    AddLatErrorBands_TruthTree(QSq_mc_truth_all_signal_BeforeFSI);
    AddLatErrorBands_TruthTree(W_mc_truth_all_signal_BeforeFSI);

    for (int i = 0; i < nFSIType; ++i){
        // Add Vertical Error Bands
        AddVertErrorBands_TruthTree(muon_P_mc_truth_all_signal_FSIType[i]);
        AddVertErrorBands_TruthTree(muon_theta_mc_truth_all_signal_FSIType[i]);
        AddVertErrorBands_TruthTree(pi0_P_mc_truth_all_signal_FSIType[i]);
        AddVertErrorBands_TruthTree(pi0_KE_mc_truth_all_signal_FSIType[i]);
        AddVertErrorBands_TruthTree(pi0_theta_mc_truth_all_signal_FSIType[i]);
        AddVertErrorBands_TruthTree(Enu_mc_truth_all_signal_FSIType[i]);
        AddVertErrorBands_TruthTree(QSq_mc_truth_all_signal_FSIType[i]);
        AddVertErrorBands_TruthTree(W_mc_truth_all_signal_FSIType[i]);

        // Add Lateral Error Bands
        AddLatErrorBands_TruthTree(muon_P_mc_truth_all_signal_FSIType[i]);
        AddLatErrorBands_TruthTree(muon_theta_mc_truth_all_signal_FSIType[i]);
        AddLatErrorBands_TruthTree(pi0_P_mc_truth_all_signal_FSIType[i]);
        AddLatErrorBands_TruthTree(pi0_KE_mc_truth_all_signal_FSIType[i]);
        AddLatErrorBands_TruthTree(pi0_theta_mc_truth_all_signal_FSIType[i]);
        AddLatErrorBands_TruthTree(Enu_mc_truth_all_signal_FSIType[i]);
        AddLatErrorBands_TruthTree(QSq_mc_truth_all_signal_FSIType[i]);
        AddLatErrorBands_TruthTree(W_mc_truth_all_signal_FSIType[i]);
    }

    for (int i = 0; i < nIntType; ++i){
        // Add Vertical Error Bands
        AddVertErrorBands_TruthTree(muon_P_mc_truth_all_signal_IntType[i]);
        AddVertErrorBands_TruthTree(muon_theta_mc_truth_all_signal_IntType[i]);
        AddVertErrorBands_TruthTree(pi0_P_mc_truth_all_signal_IntType[i]);
        AddVertErrorBands_TruthTree(pi0_KE_mc_truth_all_signal_IntType[i]);
        AddVertErrorBands_TruthTree(pi0_theta_mc_truth_all_signal_IntType[i]);
        AddVertErrorBands_TruthTree(Enu_mc_truth_all_signal_IntType[i]);
        AddVertErrorBands_TruthTree(QSq_mc_truth_all_signal_IntType[i]);
        AddVertErrorBands_TruthTree(W_mc_truth_all_signal_IntType[i]);

        // Add Lateral Error Bands
        AddLatErrorBands_TruthTree(muon_P_mc_truth_all_signal_IntType[i]);
        AddLatErrorBands_TruthTree(muon_theta_mc_truth_all_signal_IntType[i]);
        AddLatErrorBands_TruthTree(pi0_P_mc_truth_all_signal_IntType[i]);
        AddLatErrorBands_TruthTree(pi0_KE_mc_truth_all_signal_IntType[i]);
        AddLatErrorBands_TruthTree(pi0_theta_mc_truth_all_signal_IntType[i]);
        AddLatErrorBands_TruthTree(Enu_mc_truth_all_signal_IntType[i]);
        AddLatErrorBands_TruthTree(QSq_mc_truth_all_signal_IntType[i]);
        AddLatErrorBands_TruthTree(W_mc_truth_all_signal_IntType[i]);
    }

}

void CCProtonPi0_TruthAnalyzer::FillSignal_Test()
{
    for (int i = 0; i < mc_er_nPart; ++i){
        if (mc_er_ID[i] != 111) continue;
        logFile<<std::left;
        logFile.width(8); logFile<<mc_er_ID[i]<<" "; 
        logFile.width(8); logFile<<mc_er_status[i]<<" "; 
        logFile.width(8); logFile<<mc_er_mother[i]<<" ";    
        logFile.width(12); logFile<<mc_er_E[i]<<" ";      
        logFile.width(12); logFile<<truth_pi0_4P[3]<<" ";      
        logFile<<std::endl;
    }
    logFile<<std::endl;
    std::cout<<std::endl;
}

void CCProtonPi0_TruthAnalyzer::FillSignal_XSec_Variables_IntType(int type)
{
    FillHistogram(muon_P_mc_truth_all_signal_IntType[type], truth_muon_P * MeV_to_GeV);
    FillHistogram(muon_theta_mc_truth_all_signal_IntType[type], truth_muon_theta_beam * rad_to_deg);
    FillHistogram(pi0_P_mc_truth_all_signal_IntType[type], truth_pi0_P * MeV_to_GeV);
    FillHistogram(pi0_KE_mc_truth_all_signal_IntType[type], truth_pi0_KE * MeV_to_GeV);
    FillHistogram(pi0_theta_mc_truth_all_signal_IntType[type], truth_pi0_theta_beam * rad_to_deg);
    FillHistogram(Enu_mc_truth_all_signal_IntType[type], mc_incomingE * MeV_to_GeV);
    FillHistogram(QSq_mc_truth_all_signal_IntType[type], truth_QSq_exp * MeVSq_to_GeVSq);
    FillHistogram(W_mc_truth_all_signal_IntType[type], truth_W_exp * MeV_to_GeV);
}

void CCProtonPi0_TruthAnalyzer::FillSignal_XSec_Variables_FSIType(int type)
{
    FillHistogram(muon_P_mc_truth_all_signal_FSIType[type], truth_muon_P * MeV_to_GeV);
    FillHistogram(muon_theta_mc_truth_all_signal_FSIType[type], truth_muon_theta_beam * rad_to_deg);
    FillHistogram(pi0_P_mc_truth_all_signal_FSIType[type], truth_pi0_P * MeV_to_GeV);
    FillHistogram(pi0_KE_mc_truth_all_signal_FSIType[type], truth_pi0_KE * MeV_to_GeV);
    FillHistogram(pi0_theta_mc_truth_all_signal_FSIType[type], truth_pi0_theta_beam * rad_to_deg);
    FillHistogram(Enu_mc_truth_all_signal_FSIType[type], mc_incomingE * MeV_to_GeV);
    FillHistogram(QSq_mc_truth_all_signal_FSIType[type], truth_QSq_exp * MeVSq_to_GeVSq);
    FillHistogram(W_mc_truth_all_signal_FSIType[type], truth_W_exp * MeV_to_GeV);
}

void CCProtonPi0_TruthAnalyzer::FillSignal_XSec_Variables_BeforeFSI()
{
    FillHistogram(muon_P_mc_truth_all_signal_BeforeFSI, truth_muon_P_BeforeFSI * MeV_to_GeV);
    FillHistogram(muon_theta_mc_truth_all_signal_BeforeFSI, truth_muon_theta_beam_BeforeFSI * rad_to_deg);
    FillHistogram(pi0_P_mc_truth_all_signal_BeforeFSI, truth_pi0_P_BeforeFSI * MeV_to_GeV);
    FillHistogram(pi0_KE_mc_truth_all_signal_BeforeFSI, truth_pi0_KE_BeforeFSI * MeV_to_GeV);
    FillHistogram(pi0_theta_mc_truth_all_signal_BeforeFSI, truth_pi0_theta_beam_BeforeFSI * rad_to_deg);
    FillHistogram(Enu_mc_truth_all_signal_BeforeFSI, mc_incomingE * MeV_to_GeV);
    FillHistogram(QSq_mc_truth_all_signal_BeforeFSI, truth_QSq_exp_BeforeFSI * MeVSq_to_GeVSq);
    FillHistogram(W_mc_truth_all_signal_BeforeFSI, truth_W_exp_BeforeFSI * MeV_to_GeV);
}

void CCProtonPi0_TruthAnalyzer::FillSignal_XSec_Variables()
{
    FillHistogram(muon_P_mc_truth_all_signal, truth_muon_P * MeV_to_GeV);
    FillHistogram(muon_theta_mc_truth_all_signal, truth_muon_theta_beam * rad_to_deg);
    FillHistogram(pi0_P_mc_truth_all_signal, truth_pi0_P * MeV_to_GeV);
    FillHistogram(pi0_KE_mc_truth_all_signal, truth_pi0_KE * MeV_to_GeV);
    FillHistogram(pi0_theta_mc_truth_all_signal, truth_pi0_theta_beam * rad_to_deg);
    FillHistogram(Enu_mc_truth_all_signal, mc_incomingE * MeV_to_GeV);
    FillHistogram(QSq_mc_truth_all_signal, truth_QSq_exp * MeVSq_to_GeVSq);
    FillHistogram(W_mc_truth_all_signal, truth_W_exp * MeV_to_GeV);
}

void CCProtonPi0_TruthAnalyzer::FillSignal_InteractionType()
{
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
    std::cout<<std::endl;
    std::cout<<std::left;
    std::cout.width(4); std::cout<<"ID"<<" "; 
    std::cout.width(12); std::cout<<"PDG"<<" "; 
    std::cout.width(8); std::cout<<"Status"<<" ";    
    std::cout.width(8); std::cout<<"Mother"<<" ";      
    std::cout.width(12); std::cout<<"Energy"<<" ";
    std::cout<<std::endl;
    for (int i = 0; i < mc_er_nPart; ++i){
        std::cout.width(4); std::cout<<i<<" "; 
        std::cout.width(12); std::cout<<mc_er_ID[i]<<" "; 
        std::cout.width(8); std::cout<<mc_er_status[i]<<" ";    
        std::cout.width(8); std::cout<<mc_er_mother[i]<<" ";      
        std::cout.width(12); std::cout<<mc_er_E[i]<<" ";
        std::cout<<std::endl;
    }
}

void CCProtonPi0_TruthAnalyzer::writeHistograms()
{
    f->cd();

    // Test Variables
    Test_pi0_P->Write();
    Test_pi0_theta->Write();

    // Cross Section Variables
    muon_P_mc_truth_all_signal->Write();
    muon_theta_mc_truth_all_signal->Write();
    pi0_P_mc_truth_all_signal->Write();
    pi0_KE_mc_truth_all_signal->Write();
    pi0_theta_mc_truth_all_signal->Write();
    Enu_mc_truth_all_signal->Write();
    QSq_mc_truth_all_signal->Write();
    W_mc_truth_all_signal->Write();

    // Cross Section Variables -- Before FSI
    muon_P_mc_truth_all_signal_BeforeFSI->Write();
    muon_theta_mc_truth_all_signal_BeforeFSI->Write();
    pi0_P_mc_truth_all_signal_BeforeFSI->Write();
    pi0_KE_mc_truth_all_signal_BeforeFSI->Write();
    pi0_theta_mc_truth_all_signal_BeforeFSI->Write();
    Enu_mc_truth_all_signal_BeforeFSI->Write();
    QSq_mc_truth_all_signal_BeforeFSI->Write();
    W_mc_truth_all_signal_BeforeFSI->Write();

    // FSI Type
    for (int i = 0; i < nFSIType; ++i){
        muon_P_mc_truth_all_signal_FSIType[i]->Write();
        muon_theta_mc_truth_all_signal_FSIType[i]->Write();
        pi0_P_mc_truth_all_signal_FSIType[i]->Write();
        pi0_KE_mc_truth_all_signal_FSIType[i]->Write();
        pi0_theta_mc_truth_all_signal_FSIType[i]->Write();
        QSq_mc_truth_all_signal_FSIType[i]->Write();
        W_mc_truth_all_signal_FSIType[i]->Write();
        Enu_mc_truth_all_signal_FSIType[i]->Write();
    }

    // Interaction Type
    for (int i = 0; i < nIntType; ++i){
        muon_P_mc_truth_all_signal_IntType[i]->Write();
        muon_theta_mc_truth_all_signal_IntType[i]->Write();
        pi0_P_mc_truth_all_signal_IntType[i]->Write();
        pi0_KE_mc_truth_all_signal_IntType[i]->Write();
        pi0_theta_mc_truth_all_signal_IntType[i]->Write();
        QSq_mc_truth_all_signal_IntType[i]->Write();
        W_mc_truth_all_signal_IntType[i]->Write();
        Enu_mc_truth_all_signal_IntType[i]->Write();
    }

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


