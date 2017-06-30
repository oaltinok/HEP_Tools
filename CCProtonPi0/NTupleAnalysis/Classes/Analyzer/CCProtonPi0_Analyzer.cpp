#ifndef CCProtonPi0_Analyzer_cpp
#define CCProtonPi0_Analyzer_cpp

#include "CCProtonPi0_Analyzer.h"

using namespace std;

void CCProtonPi0_Analyzer::specifyRunTime()
{
    applyMaxEvents = false; 
    nMaxEvents = 10000;
    if(!m_isMC) nMaxEvents = nMaxEvents * POT_ratio;

    // Control Flow
    isDataAnalysis = true;
    isScanRun = false;
    fillErrors_ByHand = true; // Affects only Vertical Error Bands - Lateral Bands always filled ByHand

    applyGENIETuning_Complete = true;
    
    applyGENIETuning_DeltaSuppression = false;

    applyBckgConstraints_CV = true;
    applyBckgConstraints_Unv = true;

    // Side Band Control
    NoSideBand = true;
    sideBand_Michel = false;
    sideBand_PID = false;
    sideBand_LowInvMass = false;
    sideBand_HighInvMass = false;

    // EM Calibration Correction
    // Found the correction using a Double Gaussian Fit in MATLAB
    EM_MC_peak = 130.28;
    EM_Data_peak = 137.76;
    if (m_isMC) EM_correction = pi0_mass/EM_MC_peak;
    else EM_correction = pi0_mass/EM_Data_peak;

    // Event Selections
    applyProtonScore = true;
    minProtonScore_LLR = -5.0;

    applyPhotonDistance = true;
    minPhotonDistance_1 = 14; //cm
    minPhotonDistance_2 = 0; //cm

    min_Pi0_invMass = 60.0;
    max_Pi0_invMass = 200.0;

    applyDeltaInvMass = false;
    min_Delta_invMass = 40.0;
    max_Delta_invMass = 200.0;

    latest_ScanID = 0.0;

    avg_Enu = 0;
    avg_Enu_10 = 0;
    n_data = 0;
    n_data_10 = 0;
    // Counter Names
    setCounterNames();
}

void CCProtonPi0_Analyzer::reduce(string playlist)
{
    string rootDir;

    if(m_isMC) rootDir = Folder_List::rootOut + Folder_List::MC + Folder_List::reduced + "ReducedNTuple_" + version + ".root";
    if(m_isMC) rootDir = Folder_List::rootOut + Folder_List::MC + Folder_List::reduced + "ReducedNTuple_" + version + ".root";
    else rootDir = Folder_List::rootOut + Folder_List::Data + Folder_List::reduced + "ReducedNTuple_" + version + ".root";

    cout<<"Reducing NTuple Files to a single file"<<endl;
    cout<<"\tRoot File: "<<rootDir<<endl;
    TFile* f = new TFile(rootDir.c_str(),"RECREATE");
    if (!f->IsOpen()){
        cout<<"File already exists! Exiting!..."<<endl;
        exit(1);
    }

    // Create Chain and Initialize
    TChain* fChain = new TChain("CCProtonPi0");
    Init(playlist, fChain);
    if (!fChain) return;
    if (fChain == 0) return;

    // Clone Tree from Chain
    TTree* tree = fChain->CloneTree(0);

    // Get First Line for the first File
    getline(DSTFileList,scanFileName);

    //------------------------------------------------------------------------
    // Loop over Chain
    //------------------------------------------------------------------------
    cout<<"Looping over all entries"<<endl;

    Long64_t nbytes = 0, nb = 0;
    Long64_t nentries = fChain->GetEntriesFast();

    for (Long64_t jentry=0; jentry < nentries; jentry++) {

        Long64_t ientry = fChain->GetEntry(jentry);
        if (ientry < 0) break;

        nb = fChain->GetEntry(jentry);   nbytes += nb;    

        if (ientry == 0) {
            cout<<"\tGetEntry failure "<<jentry<<endl;
            break;
        }

        // Reset Counters isCounted
        counter1.isCounted = false;
        counter2.isCounted = false;
        counter3.isCounted = false;
        counter4.isCounted = false;

        // Update scanFileName if running for scan
        //if(isScanRun) UpdateScanFileName();

        // Progress Message on Terminal
        int msg_entry;
        if (m_isMC) msg_entry = 50000;
        else msg_entry = 10000;
        if (jentry%msg_entry == 0) cout<<"\tEntry "<<jentry<<endl;

        if (applyMaxEvents && jentry >= nMaxEvents){
            cout<<"\tReached Event Limit!"<<endl;
            break;
        }

        CorrectEMShowerCalibration();

        Calc_EventKinematics();
        if (truth_isSignal) Calc_EventKinematics_Truth();

        // Revise NTuple Variables -- These functions save updated versions to reduced.root
        if(m_isMC){
            if (truth_isSignal) ReviseSignal();
            if (!truth_isSignal) ReviseBackground();
        }

        CalcEventWeight();

        // Get Cut Statistics
        isPassedAllCuts = getCutStatistics();
        if( !isPassedAllCuts ) continue;

        tree->Fill();
    }

    if (!m_isMC){
        AddVertErrorBands_Data(cutList.invMass_all);
        AddVertErrorBands_Data(cutList.hCut_pi0invMass[0]);
        AddVertErrorBands_Data(cutList.SideBand_muon_P[0]);
        AddVertErrorBands_Data(cutList.SideBand_muon_theta[0]);
        AddVertErrorBands_Data(cutList.SideBand_pi0_P[0]);
        AddVertErrorBands_Data(cutList.SideBand_pi0_KE[0]);
        AddVertErrorBands_Data(cutList.SideBand_pi0_theta[0]);
        AddVertErrorBands_Data(cutList.SideBand_QSq[0]);
        AddVertErrorBands_Data(cutList.SideBand_W[0]);
        AddVertErrorBands_Data(cutList.SideBand_neutrino_E[0]);

        AddLatErrorBands_Data(cutList.invMass_all);
        AddLatErrorBands_Data(cutList.hCut_pi0invMass[0]);
        //        AddLatErrorBands_Data(cutList.SideBand_muon_P[0]);
        //        AddLatErrorBands_Data(cutList.SideBand_muon_theta[0]);
        //        AddLatErrorBands_Data(cutList.SideBand_pi0_P[0]);
        //        AddLatErrorBands_Data(cutList.SideBand_pi0_KE[0]);
        //        AddLatErrorBands_Data(cutList.SideBand_pi0_theta[0]);
        //        AddLatErrorBands_Data(cutList.SideBand_QSq[0]);
        //        AddLatErrorBands_Data(cutList.SideBand_W[0]);
        //        AddLatErrorBands_Data(cutList.SideBand_neutrino_E[0]);
    }
    cutList.writeCutTable();
    cutList.writeHistograms();

    cout<<">> Writing "<<rootDir<<endl;
    tree->AutoSave();    
    f->Write();

    //--------------------------------------------------------------------------
    // Counters
    //--------------------------------------------------------------------------
    n2p2h.print();
    nMichel_Truth.print();
    nMichel_Total_Found.print();
    nMichel_Total_Found_Improved.print();
    nMichel_Truth_Found.print();
    nMichel_Truth_Found_Improved.print();
    counter1.print();
    counter2.print();
    counter3.print();
    counter4.print();
}


void CCProtonPi0_Analyzer::analyze(string playlist)
{

    //------------------------------------------------------------------------
    // Create chain
    //------------------------------------------------------------------------
    TChain* fChain = new TChain("CCProtonPi0");
    Init(playlist, fChain);

    if (!fChain) return;
    if (fChain == 0) return;

    //------------------------------------------------------------------------
    // Loop over Chain
    //------------------------------------------------------------------------
    cout<<"Looping over all entries"<<endl;

    // Get First Line for the first File
    getline(DSTFileList,scanFileName);

    Long64_t nbytes = 0, nb = 0;
    Long64_t nentries = fChain->GetEntriesFast();
    for (Long64_t jentry=0; jentry < nentries; jentry++) {

        nb = fChain->GetEntry(jentry);   nbytes += nb;    
        Long64_t ientry = fChain->GetEntry(jentry);

        if (ientry == 0) {
            cout<<"\tGetEntry failure "<<jentry<<endl;
            break;
        }

        // Progress Message on Terminal
        int msg_entry;
        //if (m_isMC) msg_entry = 25000;
        if (m_isMC) msg_entry = 2500;
        else msg_entry = 5000;
        if (jentry%msg_entry == 0) cout<<"\tEntry "<<jentry<<endl;

        if (applyMaxEvents && jentry >= nMaxEvents){
            cout<<"\tReached Event Limit!"<<endl;
            break;
        }

        // Reset Counters isCounted
        counter1.isCounted = false;
        counter2.isCounted = false;
        counter3.isCounted = false;
        counter4.isCounted = false;

        // Reset Interaction Hist Filled
        interaction.isErrHistFilled_NeutronResponse = false;
        interaction.isErrHistFilled_PionResponse = false;
        interaction.isErrHistFilled_MuonTracking = false;

        // Update scanFileName if running for scan
        //if(isScanRun) UpdateScanFileName();

        CalcEventWeight();
        Calc_EventKinematics();
        if (truth_isSignal || truth_isSignalOut_Kinematics) Calc_EventKinematics_Truth();

        //----------------------------------------------------------------------
        // Fill Background Branches for Background Events
        //----------------------------------------------------------------------
        if(m_isMC && !truth_isSignal) {
            if (nProtonCandidates == 0) bckgTool.set_nTracks(1);
            else bckgTool.set_nTracks(2);
            bckgTool.fillBackgroundCompact(truth_isBckg_Compact_WithPi0, truth_isBckg_Compact_QELike, truth_isBckg_Compact_SinglePiPlus, truth_isBckg_Compact_Other, cvweight);                                    
            //bckgTool.fillBackgroundWithPi0(truth_isBckg_NoPi0, truth_isBckg_SinglePi0, truth_isBckg_MultiPi0, truth_isBckg_withMichel, cvweight);                                    
            //bckgTool.fillBackground(truth_isBckg_NC, truth_isBckg_AntiNeutrino, truth_isBckg_QELike, truth_isBckg_SingleChargedPion,truth_isBckg_SingleChargedPion_ChargeExchanged, truth_isBckg_DoublePionWithPi0, truth_isBckg_DoublePionWithoutPi0, truth_isBckg_MultiPionWithPi0, truth_isBckg_MultiPionWithoutPi0, truth_isBckg_Other, truth_isBckg_withMichel, cvweight);                                    
        }

        // Update Counters
        if (m_isMC) updateCounters(); 

        // Fill Reconstructed Information
        fillInteractionReco();
        fillMuonReco();
        fillPi0Reco();
        fillPi0BlobReco();
        if(nProtonCandidates > 0){
            fillProtonReco();
            deltaInvMass_reco = calcDeltaInvariantMass(false) * MeV_to_GeV;
            GetDelta_pi_angles(false); // Reco Values
        }
        
        // Fill Truth Information if Exist and Set Errors
        if( m_isMC ){
            fillInteractionMC();
            fillMuonMC();
            fillPi0MC();
            fillPi0BlobMC();
            if(nProtonCandidates > 0){
                fillProtonMC();
                deltaInvMass_true = calcDeltaInvariantMass(true) * MeV_to_GeV;
                GetDelta_pi_angles(true);  // True Values
            }
        }

        //----------------------------------------------------------------------
        // Data Analysis -- Fills Histograms with Error Bands (Slow)
        //----------------------------------------------------------------------
        if (isDataAnalysis){
            fill_XSecVars();
            if (m_isMC){
                FillLatErrorBands_ByHand();
            }
        }

        //--------------------------------------------------------------------------
        // Studies during for loop
        //--------------------------------------------------------------------------
        if (writeFSParticleMomentum) writeFSParticle4P(jentry);
        if ( isScanRun && nProtonCandidates > 0 && pi0_KE > 400) {
            writeScanFile(); 
        }

        if (truth_isSignal){
            n_data++;
            avg_Enu += m_Enu * MeV_to_GeV;
            if((m_Enu * MeV_to_GeV) < 10.0){
                n_data_10++;
                avg_Enu_10 += m_Enu * MeV_to_GeV;
            }
        }

        if (truth_isBckg_Compact_SinglePiPlus || truth_isBckg_Compact_QELike || truth_isBckg_Compact_Other){
            bool print = false;
            counter1.increment(cvweight);
            for (int i = 0; i < detmc_traj_id_sz; ++i){
                int pdg = detmc_traj_pdg[i];
                //int mother = detmc_traj_mother[i];

                if (pdg == 111) counter2.increment(cvweight);

                //if (pdg == 111 && isMotherPiPlus(mother)) {
                //    counter2.increment(cvweight);
                //    //print = true;
                //}
            }

            if (print){
                for (int i = 0; i < detmc_traj_id_sz; ++i){
                    std::cout<<i<<" "<<detmc_traj_id[i]<<" "<<detmc_traj_pdg[i]<<" "<<detmc_traj_mother[i]<<std::endl;
                }
                std::cout<<"------------------------"<<std::endl;
                std::cout<<"------------------------"<<std::endl;
                std::cout<<"------------------------"<<std::endl;
            }
        }
        //Study_Paper();
        //Study_Flux();
        //Study_BckgSubtraction();
        Study_W();
        //Study_QSq();
        //Study_ProtonSystematics();
        //Study_GENIE_Weights();
        //Study_DeltaResonance();
        //Study_NonRES();
    } // end for-loop

    if (!m_isMC) AddErrorBands_Data();

    //--------------------------------------------------------------------------
    // Studies Outside Loop 
    //--------------------------------------------------------------------------
    avg_Enu = avg_Enu / n_data;
    avg_Enu_10 = avg_Enu_10 / n_data_10;
    std::cout<<"Average Neutrino Energy = "<<avg_Enu<<std::endl;
    std::cout<<"Average Neutrino Energy less than 10 GeV = "<<avg_Enu_10<<std::endl;

    //--------------------------------------------------------------------------
    // Write Text Files
    //--------------------------------------------------------------------------
    bckgTool.writeBackgroundTable();
    //getPi0Family(); // Pi0 Family Information written inside FailFile

    //--------------------------------------------------------------------------
    // Write Root Files
    //--------------------------------------------------------------------------
    interaction.writeHistograms();
    muon.writeHistograms();
    proton.writeHistograms();
    pi0.writeHistograms();
    pi0Blob.writeHistograms();

    //--------------------------------------------------------------------------
    // Counters
    //--------------------------------------------------------------------------
    printCounters();
    writeEventTypeTable();
}

//------------------------------------------------------------------------------
//  Constructor
//------------------------------------------------------------------------------
CCProtonPi0_Analyzer::CCProtonPi0_Analyzer(bool isModeReduce, bool isMC) : 
    CCProtonPi0_NTupleAnalysis(),
    interaction(isModeReduce, isMC),
    muon(isModeReduce, isMC),
    proton(isModeReduce, isMC),
    pi0(isModeReduce, isMC),
    pi0Blob(isModeReduce, isMC),
    bckgTool(isModeReduce),
    cutList(isModeReduce, isMC),
    BckgConstrainer(Folder_List::BckgConstraints),
    QSqFitter()
{   
    cout<<"Initializing CCProtonPi0_Analyzer"<<endl;

    m_isMC = isMC;
    m_isModeReduce = isModeReduce;

    specifyRunTime();

    openTextFiles();

    initLateralErrorBandShifts(isModeReduce);

    initCVWeights();

    InformSignalType();

    cout<<"Initialization Finished!\n"<<endl;
}


//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
//
//     Specific Functions
//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

bool CCProtonPi0_Analyzer::isMotherPiPlus(int mother)
{
    for(int i = 0; i < detmc_traj_id_sz; ++i){
        int id = detmc_traj_id[i];
        int pdg = detmc_traj_pdg[i];
        if (mother == id){
            if ( std::abs(pdg) == 211) return true;
            else return false;
        }
    }

    return false;
}

void CCProtonPi0_Analyzer::UpdateScanFileName()
{
    /*
       EventID's are increasing in a file.
       If event ID is less than the latest ID then the file is changed
       */
    if(latest_ScanID >= truth_eventID){
        getline(DSTFileList,scanFileName);
    }
    latest_ScanID = truth_eventID;   
}

void CCProtonPi0_Analyzer::getPi0Family()
{         
    double nMother_DIS = 0; // dis hadronic system before hadronization

    double nMother_pi_plus = 0;
    double nMother_pi_minus = 0;

    double nMother_Delta_p_1232 = 0;
    double nMother_Delta_p_1620 = 0;
    double nMother_Delta_p_1700 = 0;

    double nMother_N_p_1440 = 0;
    double nMother_N_p_1520 = 0;
    double nMother_N_p_1535 = 0;
    double nMother_N_p_1650 = 0;
    double nMother_N_p_1675 = 0;
    double nMother_N_p_1680 = 0;

    double nMother_NoPDG = 0;
    double nMother_Other = 0;

    double nGrandMother_neutron = 0;
    double nGrandMother_proton = 0;
    double nGrandMother_Delta_pp_1232 = 0;
    double nGrandMother_Delta_pp_1620 = 0;
    double nGrandMother_Delta_pp_1700 = 0;
    double nGrandMother_NoPDG = 0;
    double nGrandMother_Other = 0;

    for(unsigned int i = 0; i < PDG_pi0_Mother.size(); i++){
        // Count Mothers
        if(PDG_pi0_Mother[i] == -5) nMother_DIS++;
        else if (PDG_pi0_Mother[i] == PDG_List::Delta_p_1232) nMother_Delta_p_1232++;
        else if(PDG_pi0_Mother[i] == PDG_List::Delta_p_1620) nMother_Delta_p_1620++;
        else if(PDG_pi0_Mother[i] == PDG_List::Delta_p_1700) nMother_Delta_p_1700++;
        else if(PDG_pi0_Mother[i] == PDG_List::N_p_1440) nMother_N_p_1440++;
        else if(PDG_pi0_Mother[i] == PDG_List::N_p_1520) nMother_N_p_1520++;
        else if(PDG_pi0_Mother[i] == PDG_List::N_p_1535) nMother_N_p_1535++;
        else if(PDG_pi0_Mother[i] == PDG_List::N_p_1650) nMother_N_p_1650++;
        else if(PDG_pi0_Mother[i] == PDG_List::N_p_1675) nMother_N_p_1675++;
        else if(PDG_pi0_Mother[i] == PDG_List::N_p_1680) nMother_N_p_1680++;
        else if(PDG_pi0_Mother[i] == PDG_List::pi_plus) nMother_pi_plus++;
        else if(PDG_pi0_Mother[i] == PDG_List::pi_minus) nMother_pi_minus++;
        else if(PDG_pi0_Mother[i] == -9) nMother_NoPDG++; 
        else{ 
            nMother_Other ++;
        }

        // Count GrandMothers
        if (PDG_pi0_GrandMother[i] == PDG_List::neutron) nGrandMother_neutron++;
        else if (PDG_pi0_GrandMother[i] == PDG_List::proton) nGrandMother_proton++;
        else if (PDG_pi0_GrandMother[i] == PDG_List::Delta_pp_1232) nGrandMother_Delta_pp_1232++;
        else if (PDG_pi0_GrandMother[i] == PDG_List::Delta_pp_1620) nGrandMother_Delta_pp_1620++;
        else if (PDG_pi0_GrandMother[i] == PDG_List::Delta_pp_1700) nGrandMother_Delta_pp_1700++;
        else if (PDG_pi0_GrandMother[i] == -9) nGrandMother_NoPDG++; 
        else{ 
            nGrandMother_Other++;
        }
    }

    cout<<">> Writing "<<failFile<<endl;

    failText<<std::left;
    failText<<"-----------------------------------------------------------------"<<endl;
    failText.width(20); failText<<"Mother"<<endl;
    failText.width(20); failText<<"DIS"<<" = "<<nMother_DIS<<endl;
    failText.width(20); failText<<"Delta_p_1232"<<" = "<<nMother_Delta_p_1232<<endl;
    failText.width(20); failText<<"Delta_p_1620"<<" = "<<nMother_Delta_p_1620<<endl;
    failText.width(20); failText<<"Delta_p_1700"<<" = "<<nMother_Delta_p_1700<<endl;
    failText.width(20); failText<<"N_p_1440"<<" = "<<nMother_N_p_1440<<endl;
    failText.width(20); failText<<"N_p_1520"<<" = "<<nMother_N_p_1520<<endl;
    failText.width(20); failText<<"N_p_1535"<<" = "<<nMother_N_p_1535<<endl;
    failText.width(20); failText<<"N_p_1650"<<" = "<<nMother_N_p_1650<<endl;
    failText.width(20); failText<<"N_p_1675"<<" = "<<nMother_N_p_1675<<endl;
    failText.width(20); failText<<"N_p_1680"<<" = "<<nMother_N_p_1680<<endl;
    failText.width(20); failText<<"pi_plus"<<" = "<<nMother_pi_plus<<endl;
    failText.width(20); failText<<"pi_minus"<<" = "<<nMother_pi_minus<<endl;
    failText.width(20); failText<<"No PDG"<<" = "<<nMother_NoPDG<<endl;
    failText.width(20); failText<<"Other"<<" = "<<nMother_Other<<endl;

    failText<<endl;
    failText.width(20); failText<<"GrandMother"<<endl;
    failText.width(20); failText<<"neutron"<<" = "<<nGrandMother_neutron<<endl;
    failText.width(20); failText<<"proton"<<" = "<<nGrandMother_proton<<endl;
    failText.width(20); failText<<"Delta_pp_1232"<<" = "<<nGrandMother_Delta_pp_1232<<endl;
    failText.width(20); failText<<"Delta_pp_1620"<<" = "<<nGrandMother_Delta_pp_1620<<endl;
    failText.width(20); failText<<"Delta_pp_1700"<<" = "<<nGrandMother_Delta_pp_1700<<endl;
    failText.width(20); failText<<"No PDG"<<" = "<<nGrandMother_NoPDG<<endl;
    failText.width(20); failText<<"Other"<<" = "<<nGrandMother_Other<<endl;
    failText<<"-----------------------------------------------------------------"<<endl;
}

void CCProtonPi0_Analyzer::fill_XSecVars()
{          
    // Fill Cross Section Variables
    fill_pi0_P();
    fill_pi0_KE();
    fill_pi0_theta();
    fill_muon_P();
    fill_muon_theta();
    fill_QSq();
    fill_Enu();
    fill_W();
    fill_deltaInvMass();
    fill_Delta_pi_theta();
    fill_Delta_pi_phi();
   
    // 2p2h Search
    fill_muon_theta_muon_KE();
    fill_q3_q0();
}

void CCProtonPi0_Analyzer::fillPi0BlobReco()
{
    // Do Nothing
}

void CCProtonPi0_Analyzer::fillPi0BlobMC()
{
    // Fill Truth Match Results
    fillPi0Blob_Pi0EvisRatio();
    fillPi0Blob_EvisStacked();
    fillPi0Blob_Evis_Fractions();
    fillPi0Blob_Evis_MostPDG();
    fillPi0Blob_Evis_Total(); 
}

void CCProtonPi0_Analyzer::fillPi0Blob_Evis_Total()
{
    double captured_fraction;
    if (truth_allClusters_evis_pizero == 0) captured_fraction = -0.5;
    else captured_fraction = truth_total_captured_evis_pizero / truth_allClusters_evis_pizero;

    FillHistogram(pi0Blob.captured_evis_frac_all, captured_fraction);
    if (truth_isSignal) FillHistogram(pi0Blob.captured_evis_frac_signal, captured_fraction);
}

void CCProtonPi0_Analyzer::fillPi0Blob_Evis_MostPDG()
{
    // Fill Most PDG -- PDG Codes modified to fit in a Histogram
    double pdg = -1;

    // Gamma 1
    if (truth_blob1_evis_most_pdg == PDG_List::pi_zero) pdg = 0; 
    else if (truth_blob1_evis_most_pdg == PDG_List::pi_plus) pdg = 1; 
    else if (truth_blob1_evis_most_pdg == PDG_List::pi_minus) pdg = 2; 
    else if (truth_blob1_evis_most_pdg == PDG_List::neutron) pdg = 3; 
    else if (truth_blob1_evis_most_pdg == PDG_List::proton) pdg = 4; 
    else if (truth_blob1_evis_most_pdg == PDG_List::mu_minus) pdg = 5;
    else pdg = 6;

    FillHistogram(pi0Blob.g1_evis_most_pdg, pdg); 

    // Gamma 2
    if (truth_blob2_evis_most_pdg == PDG_List::pi_zero) pdg = 0; 
    else if (truth_blob2_evis_most_pdg == PDG_List::pi_plus) pdg = 1; 
    else if (truth_blob2_evis_most_pdg == PDG_List::pi_minus) pdg = 2; 
    else if (truth_blob2_evis_most_pdg == PDG_List::neutron) pdg = 3; 
    else if (truth_blob2_evis_most_pdg == PDG_List::proton) pdg = 4; 
    else if (truth_blob2_evis_most_pdg == PDG_List::mu_minus) pdg = 5;
    else pdg = 6;

    FillHistogram(pi0Blob.g2_evis_most_pdg, pdg);

}

void CCProtonPi0_Analyzer::fillPi0Blob_Evis_Fractions()
{
    // Save variables in arrays (easier to Fill Histograms)
    double total[2];
    double pizero[2];
    double piplus[2];
    double piminus[2];
    double proton[2];
    double neutron[2];
    double muon[2];

    total[0] = truth_blob1_evis_total_truth;
    pizero[0] = truth_blob1_evis_pizero;
    piplus[0] = truth_blob1_evis_piplus;
    piminus[0] = truth_blob1_evis_piminus;
    proton[0] = truth_blob1_evis_proton;
    neutron[0] = truth_blob1_evis_neutron;
    muon[0] = truth_blob1_evis_muon;

    total[1] = truth_blob2_evis_total_truth;
    pizero[1] = truth_blob2_evis_pizero;
    piplus[1] = truth_blob2_evis_piplus;
    piminus[1] = truth_blob2_evis_piminus;
    proton[1] = truth_blob2_evis_proton;
    neutron[1] = truth_blob2_evis_neutron;
    muon[1] = truth_blob2_evis_muon;

    // Gamma 1
    int i = 0;
    FillHistogram(pi0Blob.g1_evis_total_truth, total[i] * MeV_to_GeV);
    FillHistogram(pi0Blob.g1_evis_frac_pizero, pizero[i]/total[i]);
    FillHistogram(pi0Blob.g1_evis_frac_piplus, piplus[i]/total[i]);
    FillHistogram(pi0Blob.g1_evis_frac_piminus, piminus[i]/total[i]);
    FillHistogram(pi0Blob.g1_evis_frac_proton, proton[i]/total[i]);
    FillHistogram(pi0Blob.g1_evis_frac_neutron, neutron[i]/total[i]);
    FillHistogram(pi0Blob.g1_evis_frac_muon, muon[i]/total[i]);

    // Gamma 2 
    i = 1;
    FillHistogram(pi0Blob.g2_evis_total_truth, total[i] * MeV_to_GeV);
    FillHistogram(pi0Blob.g2_evis_frac_pizero, pizero[i]/total[i]);
    FillHistogram(pi0Blob.g2_evis_frac_piplus, piplus[i]/total[i]);
    FillHistogram(pi0Blob.g2_evis_frac_piminus, piminus[i]/total[i]);
    FillHistogram(pi0Blob.g2_evis_frac_proton, proton[i]/total[i]);
    FillHistogram(pi0Blob.g2_evis_frac_neutron, neutron[i]/total[i]);
    FillHistogram(pi0Blob.g2_evis_frac_muon, muon[i]/total[i]);
}

void CCProtonPi0_Analyzer::writeFSParticle4P(Long64_t nEntry)
{
    for (int t = 0; t < nTopologies; t++){
        // Particle NTuple Info after All Cuts
        failText<<"----------------------------------------------------------------------"<<endl;
        failText<<nEntry<<endl;
        failText<<"Muon 4-P = ( "
            <<muon_px<<", "
            <<muon_py<<", "
            <<muon_pz<<", "
            <<muon_E<<" )"
            <<endl;
        failText<<"Proton 4-P = ( "
            <<proton_px<<", "
            <<proton_py<<", "
            <<proton_pz<<", "
            <<proton_E<<" )"
            <<" Score = "<<proton_LLRScore
            <<endl;
        failText<<"Pi0 4-P = ( "
            <<pi0_px<<", "
            <<pi0_py<<", "
            <<pi0_pz<<", "
            <<pi0_E<<" )"
            <<endl;   
    }
}

bool CCProtonPi0_Analyzer::getCutStatistics()
{
    /*
       Selection Studies 
       Assign the selection to parameters study1 and study2
       Cut Objects will count them and print them in the Cut Table
       */

    // Study 1 - Detected Michels
    // Study 2 - Missed Michels

    // Initial Study Values
    bool study1 = false;
    bool study2 = truth_isBckg_withMichel;

    //==========================================================================
    //
    // Reconstruction Cuts - Basic Selections
    //
    //==========================================================================
    cutList.nCut_All.increment(truth_isSignal, study1, study2, cvweight);

    // ------------------------------------------------------------------------
    // Vertex Cut -- If Cut_Vertex_None == 1 --> No Event Vertex
    // ------------------------------------------------------------------------
    if( Cut_Vertex_None == 1) return false;
    cutList.nCut_Vertex_None.increment(truth_isSignal, study1, study2, cvweight);

    // ------------------------------------------------------------------------
    // Vertex Reconstructable Cut
    // ------------------------------------------------------------------------
    if( Cut_Vertex_Not_Reconstructable == 1) return false;
    cutList.nCut_Vertex_Not_Reconstructable.increment(truth_isSignal, study1, study2, cvweight);

    // ------------------------------------------------------------------------
    // Vertex Fiducial Cut
    // ------------------------------------------------------------------------
    if( Cut_Vertex_Not_Fiducial == 1) return false;
    cutList.nCut_Vertex_Not_Fiducial.increment(truth_isSignal, study1, study2, cvweight);

    // ------------------------------------------------------------------------
    // Muon Cut -- If Cut_Muon_None == 1 --> No MINOS Matched Muon
    // ------------------------------------------------------------------------
    if( Cut_Muon_None == 1) return false;
    cutList.nCut_Muon_None.increment(truth_isSignal, study1, study2, cvweight);

    // Fill Truth W & Q2 for MINOS Matched Signal Events
    if (m_isMC && truth_isSignal){
        // Fill Signal Characteristics
        // true means MINOS match
        FillSignalCharacteristics(true); 
    }

    // ------------------------------------------------------------------------
    // Anti-Muon Cut
    // ------------------------------------------------------------------------
    if( Cut_Muon_Charge == 1) return false;
    cutList.nCut_Muon_Charge.increment(truth_isSignal, study1, study2, cvweight);

    double reco_muon_theta = GetCorrectedMuonTheta();
    if( reco_muon_theta*TMath::RadToDeg() > max_muon_theta) return false;
    cutList.nCut_Muon_Angle.increment(truth_isSignal, study1, study2, cvweight);

    GetMichelStatistics();

    // ------------------------------------------------------------------------
    // Michel Study 
    // ------------------------------------------------------------------------
    if ( Cut_Vertex_Michel_Exist == 1 && truth_vtx_michel_large_evis_most_pdg != -1){
        if (truth_vtx_michel_evis_most_pdg == 211){ 
            FillHistogram(cutList.michel_piplus_time_diff, vtx_michelProng_Large_time_diff);
            FillHistogram(cutList.michel_piplus_energy, vtx_michelProng_Large_energy);
            FillHistogram(cutList.michel_piplus_distance, vtx_michelProng_Large_distance);
            FillHistogram(cutList.michel_piplus_distance_z, vtx_michelProng_Large_begin_Z - vtx_z);
        }else if (truth_vtx_michel_evis_most_pdg == -211){ 
            FillHistogram(cutList.michel_piminus_time_diff, vtx_michelProng_Large_time_diff);
            FillHistogram(cutList.michel_piminus_energy, vtx_michelProng_Large_energy);
            FillHistogram(cutList.michel_piminus_distance, vtx_michelProng_Large_distance);
            FillHistogram(cutList.michel_piminus_distance_z, vtx_michelProng_Large_begin_Z - vtx_z);
        }else if (truth_vtx_michel_evis_most_pdg == 2112){ 
            FillHistogram(cutList.michel_neutron_time_diff, vtx_michelProng_Large_time_diff);
            FillHistogram(cutList.michel_neutron_energy, vtx_michelProng_Large_energy);
            FillHistogram(cutList.michel_neutron_distance, vtx_michelProng_Large_distance);
            FillHistogram(cutList.michel_neutron_distance_z, vtx_michelProng_Large_begin_Z - vtx_z);
        }else if (truth_vtx_michel_evis_most_pdg == 2212){ 
            FillHistogram(cutList.michel_proton_time_diff, vtx_michelProng_Large_time_diff);
            FillHistogram(cutList.michel_proton_energy, vtx_michelProng_Large_energy);
            FillHistogram(cutList.michel_proton_distance, vtx_michelProng_Large_distance);
            FillHistogram(cutList.michel_proton_distance_z, vtx_michelProng_Large_begin_Z - vtx_z);
        }else{ 
            FillHistogram(cutList.michel_other_time_diff, vtx_michelProng_Large_time_diff);
            FillHistogram(cutList.michel_other_energy, vtx_michelProng_Large_energy);
            FillHistogram(cutList.michel_other_distance, vtx_michelProng_Large_distance);
            FillHistogram(cutList.michel_other_distance_z, vtx_michelProng_Large_begin_Z - vtx_z);
        }
    }

    // ------------------------------------------------------------------------
    // Michel Cuts 
    // ------------------------------------------------------------------------
    isMichelEvent = (Cut_Vertex_Michel_Exist == 1) || (Cut_EndPoint_Michel_Exist == 1) || (Cut_secEndPoint_Michel_Exist == 1);
    if( isMichelEvent){
        FillHistogram(cutList.hCut_Michel,1);
    }else{
        FillHistogram(cutList.hCut_Michel,0);
    } 
    if( Cut_Vertex_Michel_Exist == 1 && !sideBand_Michel ) return false;
    cutList.nCut_Vertex_Michel_Exist.increment(truth_isSignal, study1, study2, cvweight);

    if( Cut_EndPoint_Michel_Exist == 1 && !sideBand_Michel) return false;
    cutList.nCut_EndPoint_Michel_Exist.increment(truth_isSignal, study1, study2, cvweight);

    if( Cut_secEndPoint_Michel_Exist == 1 && !sideBand_Michel) return false;
    cutList.nCut_secEndPoint_Michel_Exist.increment(truth_isSignal, study1, study2, cvweight);

    // After Michel I try to save far tracks which changes the Number of Vertices
    // Check nVertices
    FillHistogram(cutList.hCut_nVertices, vtx_total_count);

    // ------------------------------------------------------------------------
    // Tracked Particle Reconstruction Fails
    // ------------------------------------------------------------------------
    if( Cut_Particle_None == 1) return false;
    cutList.nCut_Particle_None.increment(truth_isSignal, study1, study2, cvweight);

    // ------------------------------------------------------------------------
    // Cannot Find Proton in Tracked Particles
    // ------------------------------------------------------------------------
    if( Cut_Proton_None == 1) return false;
    cutList.nCut_Proton_None.increment(truth_isSignal, study1, study2, cvweight);

    // ------------------------------------------------------------------------
    // Proton Momentum NaN
    // ------------------------------------------------------------------------
    if ( Cut_Proton_Bad == 1) return false;
    cutList.nCut_Proton_Bad.increment(truth_isSignal, study1, study2, cvweight);

    // Check nProtonCandidates
    FillHistogram(cutList.hCut_nProtonCandidates, nProtonCandidates);

    // After this stage we can analyze different topologies
    //      1Track = No Proton Events (Only Muon Track)
    //      2Track = With Proton (Muon + Proton)
    if (nProtonCandidates == 0){
        cutList.nCut_1Track_All.increment(truth_isSignal, study1, study2, cvweight);
    }else{
        cutList.nCut_2Track_All.increment(truth_isSignal, study1, study2, cvweight);
    }

    if ((isSignalDeltaRich || isSignalTwoTrack) && nProtonCandidates == 0) return false;

    // ------------------------------------------------------------------------
    // Apply Proton Score to All Proton Candidates
    // ------------------------------------------------------------------------
    isPionTrack = false;
    if (nProtonCandidates > 0){
        for( int i = 0; i < nProtonCandidates; i++){
            if ( applyProtonScore ){
                FillHistogram(cutList.hCut_2Track_protonScore_LLR,all_protons_LLRScore[i]);
                if ( all_protons_LLRScore[i] < minProtonScore_LLR ){
                    if(!sideBand_PID){
                        return false;
                    }else{ 
                        isPionTrack = true;
                        break;
                    }
                }
            }
        }
        cutList.nCut_2Track_ProtonScore.increment(truth_isSignal, study1, study2, cvweight);
    }
    cutList.nCut_ProtonScore.increment(truth_isSignal, study1, study2, cvweight);

    // Apply Proton KE Threshold
    if ((isSignalDeltaRich || isSignalTwoTrack) && !IsProtonLong(proton_KE)) return false;

    // Check nTracks 
    FillHistogram(cutList.hCut_nTracks, nTracks);
    FillHistogram(cutList.hCut_nTracks2, nTracks_Close + nTracks_Far);
    FillHistogram(cutList.hCut_nTracks_Close, nTracks_Close);
    FillHistogram(cutList.hCut_nTracks_Far, nTracks_Far);
    FillHistogram(cutList.hCut_nTracks_Discarded, nTracks_Discarded);

    // ------------------------------------------------------------------------
    // PreFilter Cut
    // ------------------------------------------------------------------------
    if(nProtonCandidates == 0){
        FillHistogram(cutList.hCut_1Track_eVis_nuclearTarget,preFilter_evis_NuclearTarget);
        FillHistogram(cutList.hCut_1Track_eVis_other,preFilter_evis_TotalExceptNuclearTarget);
    }else{
        FillHistogram(cutList.hCut_2Track_eVis_nuclearTarget,preFilter_evis_NuclearTarget);
        FillHistogram(cutList.hCut_2Track_eVis_other,preFilter_evis_TotalExceptNuclearTarget);
    }

    if( Cut_PreFilter_Pi0 == 1) return false;
    cutList.nCut_PreFilter_Pi0.increment(truth_isSignal, study1, study2, cvweight);
    if (nProtonCandidates == 0) cutList.nCut_1Track_PreFilter_Pi0.increment(truth_isSignal, study1, study2, cvweight);
    else cutList.nCut_2Track_PreFilter_Pi0.increment(truth_isSignal, study1, study2, cvweight);

    // ------------------------------------------------------------------------
    // ConeBlobs Cut -- If Cut_ConeBlobs == 1 --> Failed Pi0 Reconstruction
    // ------------------------------------------------------------------------
    FillHistogram(cutList.hCut_nShowerCandidates,anglescan_ncand); 
    if (nProtonCandidates == 0) FillHistogram(cutList.hCut_1Track_nShowerCandidates,anglescan_ncand); 
    else FillHistogram(cutList.hCut_2Track_nShowerCandidates,anglescan_ncand); 

    if( Cut_ConeBlobs == 1 ) return false;
    cutList.nCut_ConeBlobs.increment(truth_isSignal, study1, study2, cvweight);
    if (nProtonCandidates == 0) cutList.nCut_1Track_ConeBlobs.increment(truth_isSignal, study1, study2, cvweight);
    else cutList.nCut_2Track_ConeBlobs.increment(truth_isSignal, study1, study2, cvweight);

    // ------------------------------------------------------------------------
    // Blob Direction Bad Cut
    // ------------------------------------------------------------------------
    if ( Cut_BlobDirectionBad == 1 ) return false;
    cutList.nCut_BlobDirectionBad.increment(truth_isSignal, study1, study2, cvweight);
    if (nProtonCandidates == 0) cutList.nCut_1Track_BlobDirectionBad.increment(truth_isSignal, study1, study2, cvweight);
    else cutList.nCut_2Track_BlobDirectionBad.increment(truth_isSignal, study1, study2, cvweight);

    // ------------------------------------------------------------------------
    // Pi0 Momentum NaN
    // ------------------------------------------------------------------------
    if ( Cut_Pi0_Bad == 1) return false;
    cutList.nCut_Pi0_Bad.increment(truth_isSignal, study1, study2, cvweight);
    if (nProtonCandidates == 0) cutList.nCut_1Track_Pi0_Bad.increment(truth_isSignal, study1, study2, cvweight);
    else cutList.nCut_2Track_Pi0_Bad.increment(truth_isSignal, study1, study2, cvweight);

    // ------------------------------------------------------------------------
    // Check for Michel Electrons at Begin & End Points of the Showers
    // ------------------------------------------------------------------------
    bool isGamma1_Michel = gamma1_isMichel_begin || gamma1_isMichel_end;
    bool isGamma2_Michel = gamma2_isMichel_begin || gamma2_isMichel_end;
    isShower_Michel_Exist = isGamma1_Michel || isGamma2_Michel;

    if (isShower_Michel_Exist && !sideBand_Michel) return false;

    // No End point for vertical showers - so use the following
    double g1_long_dist = abs(gamma1_vertex[2] - vtx_z);
    double g2_long_dist = abs(gamma2_vertex[2] - vtx_z);
    if (Cut_Vertex_Large_Michel_Exist == 1 && g1_long_dist <= 125 ) return false;
    if (Cut_Vertex_Large_Michel_Exist == 1 && g2_long_dist <= 125 ) return false;
    cutList.nCut_Shower_Michel_Exist.increment(truth_isSignal, study1, study2, cvweight);

    // ------------------------------------------------------------------------
    // Gamma1 Conv Length Cut Hist
    // ------------------------------------------------------------------------
    if (nProtonCandidates == 0) FillHistogram(cutList.hCut_1Track_gamma1ConvDist,gamma1_dist_vtx * 0.1);
    else FillHistogram(cutList.hCut_2Track_gamma1ConvDist,gamma1_dist_vtx * 0.1);

    if (applyPhotonDistance && gamma1_dist_vtx * 0.1 < minPhotonDistance_1) return false;
    cutList.nCut_Photon1DistanceLow.increment(truth_isSignal, study1, study2, cvweight);
    if (nProtonCandidates == 0) cutList.nCut_1Track_Photon1DistanceLow.increment(truth_isSignal, study1, study2, cvweight);
    else cutList.nCut_2Track_Photon1DistanceLow.increment(truth_isSignal, study1, study2, cvweight);

    // ------------------------------------------------------------------------
    // Gamma2 Conv Length Cut
    // ------------------------------------------------------------------------
    if (nProtonCandidates == 0) FillHistogram(cutList.hCut_1Track_gamma2ConvDist,gamma2_dist_vtx * 0.1);
    else FillHistogram(cutList.hCut_2Track_gamma2ConvDist,gamma2_dist_vtx * 0.1);

    if (applyPhotonDistance && gamma2_dist_vtx * 0.1 < minPhotonDistance_2) return false;
    cutList.nCut_Photon2DistanceLow.increment(truth_isSignal, study1, study2, cvweight);
    if (nProtonCandidates == 0) cutList.nCut_1Track_Photon2DistanceLow.increment(truth_isSignal, study1, study2, cvweight);
    else cutList.nCut_2Track_Photon2DistanceLow.increment(truth_isSignal, study1, study2, cvweight);

    // Gamma Comparison
    if (truth_isSignal){
        FillHistogram(cutList.signal_gamma_E_cos_openingAngle, (gamma1_E+gamma2_E)*MeV_to_GeV, pi0_cos_openingAngle);
        FillHistogram(cutList.signal_E_cosTheta_convLength, (gamma1_E+gamma2_E)*MeV_to_GeV, pi0_cos_openingAngle, (gamma1_dist_vtx+gamma2_dist_vtx)*0.1);
    }else{
        FillHistogram(cutList.bckg_gamma_E_cos_openingAngle, (gamma1_E+gamma2_E)*MeV_to_GeV, pi0_cos_openingAngle);
        FillHistogram(cutList.bckg_E_cosTheta_convLength, (gamma1_E+gamma2_E)*MeV_to_GeV, pi0_cos_openingAngle, (gamma1_dist_vtx+gamma2_dist_vtx)*0.1);
    }
    // ------------------------------------------------------------------------
    // Low Gamma Energies AND Small Opening Angle Cut 
    // ------------------------------------------------------------------------
    if (IsOpeningAngleSmallAndEnergyLow(gamma1_E, gamma2_E)) return false; 
    cutList.nCut_LowE_SmallAngle.increment(truth_isSignal, study1, study2, cvweight);

    // ------------------------------------------------------------------------
    // Neutrino Energy Cut
    // ------------------------------------------------------------------------
    if (nProtonCandidates == 0) FillHistogram(cutList.hCut_1Track_neutrinoE,m_Enu * MeV_to_GeV);
    else FillHistogram(cutList.hCut_2Track_neutrinoE,m_Enu * MeV_to_GeV); 

    if ( !IsEnuInRange(m_Enu) ) return false;
    cutList.nCut_beamEnergy.increment(truth_isSignal, study1, study2, cvweight);
    if (nProtonCandidates == 0) cutList.nCut_1Track_beamEnergy.increment(truth_isSignal, study1, study2, cvweight);
    else cutList.nCut_2Track_beamEnergy.increment(truth_isSignal, study1, study2, cvweight);

    // ------------------------------------------------------------------------
    // W Cut
    // ------------------------------------------------------------------------
    if (!IsWInRange(m_W)) return false;

    cutList.nCut_W.increment(truth_isSignal, study1, study2, cvweight);
    if (nProtonCandidates == 0) cutList.nCut_1Track_W.increment(truth_isSignal, study1, study2, cvweight);
    else cutList.nCut_2Track_W.increment(truth_isSignal, study1, study2, cvweight);

    // ------------------------------------------------------------------------
    // Pi0 Invariant Mass Cut
    // ------------------------------------------------------------------------

    // Fill Invariant Mass Histograms
    //      If there is no Side Band, Fill for Every Event
    //      Else fill according to Side Band
    if (NoSideBand){
        FillInvMass_TruthMatch();
        FillHistogramWithVertErrors(cutList.hCut_pi0invMass, pi0_invMass);
        // Fill Lateral Error Bands on hCut_pi0invMass
        if (m_isMC){
            FillLatErrorBand_EM_EnergyScale_SideBand_invMass();
            FillLatErrorBand_MuonMomentum_SideBand_invMass();
            FillLatErrorBand_MuonTheta_SideBand_invMass();
            FillLatErrorBand_ProtonEnergy_Birks_SideBand_invMass();
            FillLatErrorBand_ProtonEnergy_SideBand_invMass("ProtonEnergy_MassModel");
            FillLatErrorBand_ProtonEnergy_SideBand_invMass("ProtonEnergy_MEU");
            FillLatErrorBand_ProtonEnergy_SideBand_invMass("ProtonEnergy_BetheBloch");
        }
        if (nProtonCandidates == 0){
            FillHistogram(cutList.pi0_invMass_1Track, pi0_invMass);
            FillHistogram(cutList.hCut_1Track_pi0invMass,pi0_invMass);
        }else{ 
            FillHistogram(cutList.pi0_invMass_2Track, pi0_invMass);
            FillHistogram(cutList.hCut_2Track_pi0invMass, pi0_invMass);
        }
    }else{
        fill_SideBand_InvMass();
    }

    if (pi0_invMass < min_Pi0_invMass) isLowInvMassEvent = true;
    else isLowInvMassEvent = false;

    if (pi0_invMass > max_Pi0_invMass) isHighInvMassEvent = true;
    else isHighInvMassEvent = false;

    if( isLowInvMassEvent && !sideBand_LowInvMass) return false;
    if( isHighInvMassEvent && !sideBand_HighInvMass) return false;

    fill_BackgroundSubtractionHists();

    cutList.nCut_Pi0_invMass.increment(truth_isSignal, study1, study2, cvweight);
    if (nProtonCandidates == 0) cutList.nCut_1Track_Pi0_invMass.increment(truth_isSignal, study1, study2, cvweight);
    else cutList.nCut_2Track_Pi0_invMass.increment(truth_isSignal, study1, study2, cvweight);

    // Fill Other Side Bands 
    //      If there is no Side Band, Fill for Every Event
    //      Else fill according to Side Band
    if (NoSideBand){
        double reco_theta = GetCorrectedMuonTheta();
        FillHistogramWithVertErrors(cutList.SideBand_muon_P, muon_P*MeV_to_GeV);
        FillHistogramWithVertErrors(cutList.SideBand_muon_theta, reco_theta*TMath::RadToDeg());
        FillHistogramWithVertErrors(cutList.SideBand_pi0_P, pi0_P*MeV_to_GeV);
        FillHistogramWithVertErrors(cutList.SideBand_pi0_KE, pi0_KE*MeV_to_GeV);
        FillHistogramWithVertErrors(cutList.SideBand_pi0_theta, pi0_theta_beam*TMath::RadToDeg());
        FillHistogramWithVertErrors(cutList.SideBand_neutrino_E, m_Enu*MeV_to_GeV);
        FillHistogramWithVertErrors(cutList.SideBand_QSq, m_QSq*MeVSq_to_GeVSq);
        FillHistogramWithVertErrors(cutList.SideBand_W, m_W*MeV_to_GeV);
    }else{
        fill_SideBand_Other();
    }

    //-------------------------------------------------------------------------
    // Satisfied All Cuts
    //-------------------------------------------------------------------------
    return true;
}

void CCProtonPi0_Analyzer::fillInteractionMC()
{
    if(truth_isSignal){
        // Fill for selected sample
        //  false means -- Not MINOS match only
        FillSignalCharacteristics(false);

        // Short Proton True Information
        if (nProtonCandidates == 0){
            double proton_mass = 938.27; // MeV
            double proton_true_P = HEP_Functions::calcMomentum(truth_proton_4P[0],truth_proton_4P[1],truth_proton_4P[2]);
            double proton_true_KE = truth_proton_4P[3] - proton_mass;
            FillHistogram(interaction.proton_true_P_1Track, proton_true_P);
            FillHistogram(interaction.proton_true_KE_1Track, proton_true_KE);
            FillHistogram(interaction.proton_true_theta_1Track, truth_proton_theta * rad_to_deg);
        }

        // Ejected Nucleon Count
        double n_nucleons = GetEjectedNucleonCount();
        if (nProtonCandidates == 0){
            FillHistogram(interaction.n_ejected_nucleons_1Track, n_nucleons);
        }else{
            FillHistogram(interaction.n_ejected_nucleons_2Track, n_nucleons);
        }

        // W: Error, Difference
        double W_true = m_W_Truth * MeV_to_GeV;
        double W_reco = m_W * MeV_to_GeV;
        double W_Error = Data_Functions::getError(W_true, W_reco);

        FillHistogram(interaction.W_Error, W_Error);
        FillHistogram(interaction.W_Diff, W_reco-W_true);

        // Neutrino Energy: Truth, Error, Difference
        double E_true = mc_incomingE * MeV_to_GeV;
        double E_reco = m_Enu * MeV_to_GeV;
        double E_Error = Data_Functions::getError(E_true, E_reco);

        FillHistogram(interaction.Enu_All_response, E_reco, E_true);
        FillHistogram(interaction.Enu_Error, E_Error);
        FillHistogram(interaction.Enu_Diff, E_reco-E_true);

        // Fill 1Track and 2 Track Enu 
        if(nProtonCandidates == 0){
            FillHistogram(interaction.Enu_1Track_response, E_reco, E_true);
            FillHistogram(interaction.Enu_1Track_Error, E_Error);
            FillHistogram(interaction.Enu_1Track_Diff, E_reco-E_true);
        }else{ 
            FillHistogram(interaction.Enu_2Track_response, E_reco, E_true);
            FillHistogram(interaction.Enu_2Track_Error, E_Error);
            FillHistogram(interaction.Enu_2Track_Diff, E_reco-E_true);
        }  

        // QSq True, Error and Difference
        double QSq_true = m_QSq_Truth * MeVSq_to_GeVSq;
        double QSq_reco = m_QSq * MeVSq_to_GeVSq;
        double QSq_error = Data_Functions::getError(QSq_true,QSq_reco);

        FillHistogram(interaction.QSq_All_response, QSq_reco, QSq_true);
        FillHistogram(interaction.QSq_Error, QSq_error);
        FillHistogram(interaction.QSq_Diff, QSq_reco-QSq_true);

        if(nProtonCandidates == 0){
            FillHistogram(interaction.QSq_1Track_response, QSq_reco, QSq_true);
            FillHistogram(interaction.QSq_1Track_Error, QSq_error);
            FillHistogram(interaction.QSq_1Track_Diff, QSq_reco-QSq_true);
        }else{ 
            FillHistogram(interaction.QSq_2Track_response, QSq_reco, QSq_true);
            FillHistogram(interaction.QSq_2Track_Error, QSq_error);
            FillHistogram(interaction.QSq_2Track_Diff, QSq_reco-QSq_true);
        }  
    }else{
        FillBackgroundCharacteristics();
    }
}

void CCProtonPi0_Analyzer::fillInteractionReco()
{
    FillHistogram(interaction.CV_weight, cvweight);
    FillHistogram(interaction.CV_weight_Flux, cvweight_Flux);
    FillHistogram(interaction.CV_weight_2p2h, cvweight_2p2h);
    FillHistogram(interaction.CV_weight_Delta, cvweight_Delta);
    FillHistogram(interaction.CV_weight_CCRES, cvweight_CCRES);
    FillHistogram(interaction.CV_weight_NonRes1pi, cvweight_NonRes1pi);

    FillHistogram(interaction.err_2p2h, 1-err_2p2h);
    FillHistogram(interaction.err_2p2h, 1+err_2p2h);

    FillHistogram(interaction.genie_wgt_VecFFCCQEshape, truth_genie_wgt_VecFFCCQEshape[2]);
    FillHistogram(interaction.genie_wgt_VecFFCCQEshape, truth_genie_wgt_VecFFCCQEshape[4]);

    FillHistogram(interaction.genie_wgt_NormDISCC, truth_genie_wgt_NormDISCC[2]);
    FillHistogram(interaction.genie_wgt_NormDISCC, truth_genie_wgt_NormDISCC[4]);

    FillHistogram(interaction.genie_wgt_Theta_Delta2Npi, truth_genie_wgt_Theta_Delta2Npi[2]);
    FillHistogram(interaction.genie_wgt_Theta_Delta2Npi, truth_genie_wgt_Theta_Delta2Npi[4]);

    FillHistogram(interaction.updated_wgt_Theta_Delta2Npi, updated_genie_wgt_Theta_Delta2Npi[2]);
    FillHistogram(interaction.updated_wgt_Theta_Delta2Npi, updated_genie_wgt_Theta_Delta2Npi[4]);

    FillHistogram(interaction.genie_wgt_MaRES, truth_genie_wgt_MaRES[2]);
    FillHistogram(interaction.genie_wgt_MaRES, truth_genie_wgt_MaRES[4]);

    FillHistogram(interaction.updated_wgt_MaRES, updated_genie_wgt_MaRES[2]);
    FillHistogram(interaction.updated_wgt_MaRES, updated_genie_wgt_MaRES[4]);

    FillHistogram(interaction.genie_wgt_MvRES, truth_genie_wgt_MvRES[2]);
    FillHistogram(interaction.genie_wgt_MvRES, truth_genie_wgt_MvRES[4]);

    FillHistogram(interaction.updated_wgt_MvRES, updated_genie_wgt_MvRES[2]);
    FillHistogram(interaction.updated_wgt_MvRES, updated_genie_wgt_MvRES[4]);

    FillHistogram(interaction.genie_wgt_Rvn1pi, truth_genie_wgt_Rvn1pi[2]);
    FillHistogram(interaction.genie_wgt_Rvn1pi, truth_genie_wgt_Rvn1pi[4]);

    FillHistogram(interaction.updated_wgt_Rvn1pi, updated_genie_wgt_Rvn1pi[2]);
    FillHistogram(interaction.updated_wgt_Rvn1pi, updated_genie_wgt_Rvn1pi[4]);

    // Inclusive -- All Events
    FillHistogram(interaction.Enu, m_Enu * MeV_to_GeV);
    FillHistogram(interaction.QSq, m_QSq * MeVSq_to_GeVSq);
    FillHistogram(interaction.WSq, m_WSq * MeVSq_to_GeVSq);
    FillHistogram(interaction.W, m_W * MeV_to_GeV);
    FillHistogram(interaction.WSq_QSq_Diff, m_WSq*MeVSq_to_GeVSq, (m_QSq-mc_Q2)*MeVSq_to_GeVSq);

    // Different Topologies
    if (nProtonCandidates == 0){ 
        FillHistogram(interaction.Enu_1Track, m_Enu * MeV_to_GeV);
        FillHistogram(interaction.QSq_1Track, m_QSq * MeVSq_to_GeVSq);
        FillHistogram(interaction.WSq_1Track, m_WSq * MeVSq_to_GeVSq);
        FillHistogram(interaction.W_1Track, m_W * MeV_to_GeV);
    }else{ 
        FillHistogram(interaction.Enu_2Track, m_Enu * MeV_to_GeV);
        FillHistogram(interaction.QSq_2Track, m_QSq * MeVSq_to_GeVSq);
        FillHistogram(interaction.WSq_2Track, m_WSq * MeVSq_to_GeVSq);
        FillHistogram(interaction.W_2Track, m_W * MeV_to_GeV);
    }

    // Vertex Energy & Evis 
    FillHistogram(interaction.vertex_energy_All, vertex_blob_energy);
    FillHistogram(interaction.vertex_evis_All, vertex_blob_evis );    

    if (nProtonCandidates == 0){
        FillHistogram(interaction.vertex_energy_1Track, vertex_blob_energy);
        FillHistogram(interaction.vertex_evis_1Track, vertex_blob_evis );
    }else{
        FillHistogram(interaction.vertex_energy_2Track, vertex_blob_energy);
        FillHistogram(interaction.vertex_evis_2Track, vertex_blob_evis );
    }

    // Extra Energy
    if (nProtonCandidates == 0){
        FillHistogram(interaction.extra_leftover_energy_1Track, Extra_Energy_Leftover); 
        FillHistogram(interaction.extra_muon_energy_1Track, Extra_Energy_Muon); 
        FillHistogram(interaction.extra_rejected_energy_1Track, Extra_Energy_Rejected); 
        FillHistogram(interaction.extra_total_energy_1Track, Extra_Energy_Total);
    }else{
        FillHistogram(interaction.extra_leftover_energy_2Track, Extra_Energy_Leftover); 
        FillHistogram(interaction.extra_muon_energy_2Track, Extra_Energy_Muon); 
        FillHistogram(interaction.extra_rejected_energy_2Track, Extra_Energy_Rejected); 
        FillHistogram(interaction.extra_total_energy_2Track, Extra_Energy_Total);
    }

    // Extra Energy
    FillHistogram(interaction.h_extra_leftover_energy, Extra_Energy_Leftover); 
    FillHistogram(interaction.h_extra_muon_energy, Extra_Energy_Muon); 
    FillHistogram(interaction.h_extra_rejected_energy, Extra_Energy_Rejected); 

    FillHistogram(interaction.vertex_z, vtx_z);

    if (m_isMC && truth_isSignal) FillSignalCharacteristics_Reco(); 
}

double CCProtonPi0_Analyzer::calcDeltaInvariantMass(TLorentzVector proton, TLorentzVector pi0)
{
    TLorentzVector delta_4P = proton + pi0;
    double M_delta = delta_4P.M(); 

    return M_delta;
}

double CCProtonPi0_Analyzer::calcDeltaInvariantMass(bool isTruth)
{
    TLorentzVector delta_4P = GetDelta4P(isTruth);
    double M_delta = delta_4P.M(); 

    return M_delta;
}

void CCProtonPi0_Analyzer::writeScanFile()
{
    if(isScanRun){
        //// Constants for Roundup List
        //const string arachne_html = "http://minerva05.fnal.gov/Arachne/arachne.html?filename=";
        //const string entryString  = "&entry=";
        //const string other        = "&slice=-1&filetype=dst";
        ////http://minerva05.fnal.gov/Arachne/arachne.html?det=MV&recoVer=v10r6p13&run=3596&subrun=6&gate=597&slice=7
        //roundupText<<arachne_html<<scanFileName<<entryString<<truth_eventID<<other<<" ";
        //roundupText<<mc_run<<" ^ "<<mc_subrun<<" ^ "<<ev_gate<<" ^ "<<mc_incomingE<<endl;


        //http://minerva05.fnal.gov/Arachne/arachne.html?det=MV&recoVer=v10r8p4&run=2000&subrun=2&gate=179&slice=-1
        // Constants for Roundup List
        const string arachne_html = "http://minerva05.fnal.gov/Arachne/arachne.html?det=MV&recoVer=v10r8p4";
        roundupText<<arachne_html<<"&run="<<ev_run<<"&subrun="<<ev_subrun<<"&gate="<<ev_gate<<"&slice="<<slice_numbers[0]<<" ";
        roundupText<<ev_run<<" | "<<ev_subrun<<" | "<<ev_gate<<" "<<gamma1_E<<" "<<gamma2_E<<endl;

    }else{
        cout<<"WARNING! ScanRun is NOT Activated! Are you sure what you are doing?"<<endl;    
    }
}

void CCProtonPi0_Analyzer::closeTextFiles()
{
    if(isScanRun){
        roundupText.close();
        DSTFileList.close();
    }

    failText.close();
    eventTypeTable.close();
}

void CCProtonPi0_Analyzer::openTextFiles()
{
    cout<<"Opening Text Files:"<<endl;

    std::string temp_fName;

    // Open Fail-Check File
    temp_fName = Folder_List::output + Folder_List::textOut + "FailChecks.txt";
    OpenTextFile(temp_fName, failText);

    temp_fName = Folder_List::output + Folder_List::textOut + "EventTypeTable.txt";
    OpenTextFile(temp_fName, eventTypeTable);

    if(isScanRun){
        // Open Roundup Text for Arachne Scanning
        temp_fName = Folder_List::output + Folder_List::textOut + "ArachneRoundup.txt";
        OpenTextFile(temp_fName, roundupText);

        //std::string playlistDST = "Input/Playlists/pl_Scan.dat";
        //DSTFileList.open( playlistDST.c_str() );
        //if( !DSTFileList.is_open() ){
            //cerr<<"Cannot open input text file: "<<playlistDST<<endl;
            //exit(1);
        //}else{
            //cout<<"\t"<<playlistDST<<endl;
        //}
    }

    cout<<"Done!"<<endl;
}

void CCProtonPi0_Analyzer::fillProtonMC()
{  
    if (truth_isSignal){
        // Momentum
        double reco_P = proton_P * MeV_to_GeV;
        double true_P = truth_proton_P * MeV_to_GeV;
        double error_P = Data_Functions::getError(true_P, reco_P);
        FillHistogram(proton.P_error, error_P);
        FillHistogram(proton.P_Diff, reco_P-true_P);
        FillHistogram(proton.proton_P_response, reco_P, true_P);

        // Energy
        double reco_E = proton_E * MeV_to_GeV;
        double true_E = truth_proton_4P[3] * MeV_to_GeV; 
        double error_E = Data_Functions::getError(true_E, reco_E);
        FillHistogram(proton.reco_E_true_E, reco_E,true_E);
        FillHistogram(proton.E_error, error_E);
        FillHistogram(proton.E_Diff, reco_E-true_E);

        // Proton Theta
        double reco_theta = proton_theta_beam * TMath::RadToDeg();
        double true_theta = truth_proton_theta_beam * TMath::RadToDeg();
        double error_theta = Data_Functions::getError(true_theta, reco_theta);

        FillHistogram(proton.theta_error, error_theta);
        FillHistogram(proton.theta_diff, reco_theta-true_theta);
        FillHistogram(proton.proton_theta_response, proton_theta_beam * TMath::RadToDeg(), truth_proton_theta_beam * TMath::RadToDeg());
    
        // Proton Kinetic Energy
        double proton_KE = (truth_proton_4P[3] - proton_mass) * MeV_to_GeV;
        FillHistogram(proton.proton_KE_mc_truth_signal, proton_KE);
    }

}

void CCProtonPi0_Analyzer::fillProtonReco()
{  
    // Unique Histograms
    FillHistogram(proton.partScore, proton_LLRScore);
    FillHistogram(proton.trackLength, proton_length * mm_to_cm);
    FillHistogram(proton.trackKinked, proton_kinked);

    // Standard Histograms
    FillHistogram(proton.E, proton_E * MeV_to_GeV);
    FillHistogram(proton.P, proton_P * MeV_to_GeV);
    FillHistogram(proton.KE, proton_KE * MeV_to_GeV);
    FillHistogram(proton.theta, proton_theta * TMath::RadToDeg());
    FillHistogram(proton.phi, proton_phi * TMath::RadToDeg());
}

void CCProtonPi0_Analyzer::fillPi0MC()
{
    // EM Shower Energy Variables
    if ( truth_isSignal ){
        // Gamma 1
        double g1_reco_E = gamma1_E * MeV_to_GeV;
        double g1_true_E = truth_gamma1_4P[3] * MeV_to_GeV;
        double g1_E_error = Data_Functions::getError(g1_true_E, g1_reco_E);

        FillHistogram(pi0.gamma1_true_E, g1_true_E);
        FillHistogram(pi0.gamma1_reco_error_E, g1_E_error);
        FillHistogram(pi0.gamma1_reco_E_true_E, g1_reco_E, g1_true_E);
        FillHistogram(pi0.gamma1_true_E_reco_E_error, g1_true_E,g1_E_error);

        // Gamma 2 
        double g2_reco_E = gamma2_E * MeV_to_GeV;
        double g2_true_E = truth_gamma2_4P[3] * MeV_to_GeV;
        double g2_E_error = Data_Functions::getError(g2_true_E, g2_reco_E);

        FillHistogram(pi0.gamma2_true_E, g2_true_E);
        FillHistogram(pi0.gamma2_reco_error_E, g2_E_error);
        FillHistogram(pi0.gamma2_reco_E_true_E, g2_reco_E, g2_true_E);
        FillHistogram(pi0.gamma2_true_E_reco_E_error, g2_true_E, g2_E_error);

        // Pi0 Momentum
        double pi0_reco_P = pi0_P * MeV_to_GeV;
        double pi0_true_P = truth_pi0_P * MeV_to_GeV;
        double pi0_P_error = Data_Functions::getError(pi0_true_P, pi0_reco_P);
        double pi0_P_Diff = pi0_reco_P - pi0_true_P;

        FillHistogram(pi0.P_error, pi0_P_error); 
        FillHistogram(pi0.P_Diff, pi0_P_Diff); 
        FillHistogram(pi0.reco_P_true_P, pi0_reco_P, pi0_true_P); 

        // Pi0 Energy
        double reco_E = pi0_E * MeV_to_GeV;
        //double reco_E = pi0_E_Cal * MeV_to_GeV;
        double true_E = truth_pi0_4P[3] * MeV_to_GeV; 
        double error_E = Data_Functions::getError(true_E, reco_E);

        FillHistogram(pi0.E_true, true_E);
        FillHistogram(pi0.E_reco, reco_E);
        FillHistogram(pi0.reco_E_true_E, reco_E,true_E);
        FillHistogram(pi0.E_error, error_E);
        FillHistogram(pi0.E_Diff, reco_E-true_E);

        double reco_KE = pi0_KE * MeV_to_GeV;
        double true_KE = truth_pi0_KE * MeV_to_GeV; 
        double error_KE = Data_Functions::getError(true_KE, reco_KE);

        FillHistogram(pi0.KE_true, true_KE);
        FillHistogram(pi0.KE_reco, reco_KE);
        FillHistogram(pi0.reco_KE_true_KE, reco_KE,true_KE);
        FillHistogram(pi0.KE_error, error_KE);
        FillHistogram(pi0.KE_Diff, reco_KE-true_KE);

        // Pi0 Theta
        double reco_theta = pi0_theta_beam * TMath::RadToDeg();
        double true_theta = truth_pi0_theta_beam * TMath::RadToDeg();
        double error_theta = Data_Functions::getError(true_theta, reco_theta);

        FillHistogram(pi0.theta_error, error_theta);
        FillHistogram(pi0.theta_diff, reco_theta-true_theta);
    }

    // Gamma Comparison
    if (truth_isSignal){
        FillHistogram(pi0.signal_gamma1_E_gamma2_E, gamma1_E, gamma2_E);
        FillHistogram(pi0.signal_gamma1_convLength_gamma2_convLength, gamma1_dist_vtx * 0.1, gamma2_dist_vtx * 0.1);
    }else{
        FillHistogram(pi0.bckg_gamma1_E_gamma2_E, gamma1_E, gamma2_E);
        FillHistogram(pi0.bckg_gamma1_convLength_gamma2_convLength, gamma1_dist_vtx * 0.1, gamma2_dist_vtx * 0.1);
    }
}

void CCProtonPi0_Analyzer::fillPi0Blob_EvisStacked()
{   
    const double min_evis = 10;

    // Gamma 1
    double g1_pi0 = truth_blob1_evis_pizero;
    double g1_pi = truth_blob1_evis_piplus + truth_blob1_evis_piminus;
    double g1_proton = truth_blob1_evis_proton;
    double g1_neutron = truth_blob1_evis_neutron;
    double g1_muon = truth_blob1_evis_muon;

    if (g1_pi0 > min_evis) FillHistogram(pi0Blob.g1_evis_pi0, g1_pi0);
    if (g1_pi > min_evis) FillHistogram(pi0Blob.g1_evis_pi, g1_pi);
    if (g1_proton > min_evis) FillHistogram(pi0Blob.g1_evis_proton, g1_proton);
    if (g1_neutron > min_evis) FillHistogram(pi0Blob.g1_evis_neutron, g1_neutron);
    if (g1_muon > min_evis) FillHistogram(pi0Blob.g1_evis_muon, g1_muon);

    // Gamma 2
    double g2_pi0 = truth_blob2_evis_pizero;
    double g2_pi = truth_blob2_evis_piplus + truth_blob2_evis_piminus;
    double g2_proton = truth_blob2_evis_proton;
    double g2_neutron = truth_blob2_evis_neutron;
    double g2_muon = truth_blob2_evis_muon;

    if (g2_pi0 > min_evis) FillHistogram(pi0Blob.g2_evis_pi0, g2_pi0);
    if (g2_pi > min_evis) FillHistogram(pi0Blob.g2_evis_pi, g2_pi);
    if (g2_proton > min_evis) FillHistogram(pi0Blob.g2_evis_proton, g2_proton);
    if (g2_neutron > min_evis) FillHistogram(pi0Blob.g2_evis_neutron, g2_neutron);
    if (g2_muon > min_evis) FillHistogram(pi0Blob.g2_evis_muon, g2_muon);

    // Gamma 1 + Gamma 2 (Pi0) written as Gamma 3
    double g3_pi0 = g1_pi0 + g2_pi0;
    double g3_pi = g1_pi + g2_pi;
    double g3_proton = g1_proton + g2_proton;
    double g3_neutron = g1_neutron + g2_neutron;
    double g3_muon = g1_muon + g2_muon;

    if (g3_pi0 > min_evis) FillHistogram(pi0Blob.g3_evis_pi0, g3_pi0);
    if (g3_pi > min_evis) FillHistogram(pi0Blob.g3_evis_pi, g3_pi);
    if (g3_proton > min_evis) FillHistogram(pi0Blob.g3_evis_proton, g3_proton);
    if (g3_neutron > min_evis) FillHistogram(pi0Blob.g3_evis_neutron, g3_neutron);
    if (g3_muon > min_evis) FillHistogram(pi0Blob.g3_evis_muon, g3_muon);
}

void CCProtonPi0_Analyzer::fillPi0Blob_Pi0EvisRatio()
{
    double reco_pi0_true_pi0;
    double true_pi0_reco_all;
    double reco_pi0_reco_all;

    /*
     *  true_pi0 = True Pi0 Visible Energy Inside Detector
     *  reco_pi0 = Reco Pi0 Visible Energy Coming from Pi0
     *  reco_all = Reco Total Visible Energy 
     */
    if (truth_allClusters_evis_pizero == 0){
        reco_pi0_true_pi0 = 0;
    }else{
        reco_pi0_true_pi0 = truth_total_captured_evis_pizero / truth_allClusters_evis_pizero;   
    }

    if ( truth_total_captured_evis_total_truth == 0){
        reco_pi0_reco_all = 0;
        true_pi0_reco_all = 0;
    }else{
        reco_pi0_reco_all = truth_total_captured_evis_pizero / truth_total_captured_evis_total_truth;
        true_pi0_reco_all = truth_allClusters_evis_pizero / truth_total_captured_evis_total_truth;
    }

    FillHistogram(pi0Blob.evis_frac_reco_pi0_true_pi0, reco_pi0_true_pi0);
    FillHistogram(pi0Blob.evis_frac_true_pi0_reco_all, true_pi0_reco_all);
    FillHistogram(pi0Blob.evis_frac_reco_pi0_reco_all, reco_pi0_reco_all);
    FillHistogram(pi0Blob.evis_frac_reco_nonpi0_reco_all, 1-reco_pi0_reco_all);

}



void CCProtonPi0_Analyzer::fillPi0Reco()
{
    // Unique Histograms
    FillHistogram(pi0.invMass, pi0_invMass);
    FillHistogram(pi0.cos_openingAngle, pi0_cos_openingAngle);

    // Leading Photon - Energetic Photon
    FillHistogram(pi0.gamma1_ConvLength, gamma1_dist_vtx * 0.1);
    FillHistogram(pi0.gamma1_E, gamma1_E * MeV_to_GeV);
    FillHistogram(pi0.gamma1_theta, gamma1_theta * TMath::RadToDeg());

    // Secondary Photon
    FillHistogram(pi0.gamma2_ConvLength, gamma2_dist_vtx * 0.1);
    FillHistogram(pi0.gamma2_E, gamma2_E * MeV_to_GeV);
    FillHistogram(pi0.gamma2_theta, gamma2_theta * TMath::RadToDeg());

    double photon_E_asym = abs((gamma1_E - gamma2_E) / (gamma1_E + gamma2_E));  
    FillHistogram(pi0.photonEnergy_Asymmetry, photon_E_asym);

    // Standard Histograms
    FillHistogram(pi0.E, pi0_E * MeV_to_GeV);
    FillHistogram(pi0.P, pi0_P * MeV_to_GeV);
    FillHistogram(pi0.KE, pi0_KE * MeV_to_GeV);
    FillHistogram(pi0.theta, pi0_theta * TMath::RadToDeg());
    FillHistogram(pi0.phi, pi0_phi * TMath::RadToDeg());

}

void CCProtonPi0_Analyzer::fillMuonMC()
{
    muon.muon_P_shift->Fill(muon_E_shift);

    if(truth_isSignal){ 
        // Momentum
        double reco_P = muon_P * MeV_to_GeV;
        double true_P = truth_muon_P * MeV_to_GeV;
        double error_P = Data_Functions::getError(true_P, reco_P);

        FillHistogram(muon.reco_P_true_P, reco_P,true_P);
        FillHistogram(muon.P_error, error_P);

        if (true_P < 1.0){
            //std::cout<<"true_P = "<<true_P<<" mc_incomingE = "<<mc_incomingE*MeV_to_GeV<<" recoP = "<<reco_P<<std::endl;
        }
        // Energy
        double reco_E = muon_E * MeV_to_GeV;
        double true_E = truth_muon_4P[3] * MeV_to_GeV; 
        double error_E = Data_Functions::getError(true_E, reco_E);

        FillHistogram(muon.reco_E_true_E, reco_E,true_E);
        FillHistogram(muon.E_error, error_E);
        FillHistogram(muon.E_Diff, reco_E-true_E);

        // Theta
        double reco_theta = GetCorrectedMuonTheta() * TMath::RadToDeg();
        double true_theta = truth_muon_theta_beam * TMath::RadToDeg();
        double error_theta = Data_Functions::getError(true_theta, reco_theta);

        FillHistogram(muon.theta_error, error_theta);
        FillHistogram(muon.theta_diff, reco_theta-true_theta);

        // Cosine Theta
        double reco_cos_theta = cos(reco_theta);
        double true_cos_theta = cos(truth_muon_theta_beam);
        double error_cos_theta = Data_Functions::getError(true_cos_theta, reco_cos_theta);

        FillHistogram(muon.cos_theta_error, error_cos_theta);
    }
}

void CCProtonPi0_Analyzer::fill_BackgroundSubtractionHists() 
{
    if (m_isMC){
        // Fill CV and Vertical Error Bands
        FillHistogramWithVertErrors(cutList.invMass_mc_reco_all, pi0_invMass);
        if (truth_isSignal)  FillHistogramWithVertErrors(cutList.invMass_mc_reco_signal, pi0_invMass);
        else FillHistogramWithVertErrors(cutList.invMass_mc_reco_bckg, pi0_invMass);

        // Fill Lateral Error Bands
        FillLatErrorBand_EM_EnergyScale_invMass();
        FillLatErrorBand_MuonMomentum_invMass();
        FillLatErrorBand_MuonTheta_invMass();
        FillLatErrorBand_ProtonEnergy_Birks_invMass();
        FillLatErrorBand_ProtonEnergy_invMass("ProtonEnergy_MassModel");
        FillLatErrorBand_ProtonEnergy_invMass("ProtonEnergy_MEU");
        FillLatErrorBand_ProtonEnergy_invMass("ProtonEnergy_BetheBloch");
    }else{
        FillHistogram(cutList.invMass_all, pi0_invMass);
    }
}

void CCProtonPi0_Analyzer::fill_W() 
{
    if (m_isMC){
        // MC Reco All
        FillHistogramWithVertErrors(interaction.W_mc_reco_all, m_W * MeV_to_GeV);
        if (truth_isSignal){
            // MC Truth Signal
            FillHistogramWithVertErrors(interaction.W_mc_truth_signal, m_W_Truth * MeV_to_GeV);
            // MC Reco Signal
            FillHistogramWithVertErrors(interaction.W_mc_reco_signal, m_W * MeV_to_GeV);
            // MC Reco vs True -- Response
            FillHistogramWithVertErrors(interaction.W_response, m_W * MeV_to_GeV, m_W_Truth * MeV_to_GeV);
        }else{
            // MC Reco Background
            FillHistogramWithVertErrors(interaction.W_mc_reco_bckg, m_W * MeV_to_GeV);
        }
    }else{
        // Data
        FillHistogram(interaction.W_all, m_W * MeV_to_GeV);
    }
}

void CCProtonPi0_Analyzer::fill_Enu() 
{
    if (m_isMC){
        // MC Reco All
        FillHistogramWithVertErrors(interaction.Enu_mc_reco_all, m_Enu * MeV_to_GeV);
        if (truth_isSignal){
            // MC Truth Signal
            FillHistogramWithVertErrors(interaction.Enu_mc_truth_signal, m_Enu_Truth * MeV_to_GeV);
            // MC Reco Signal
            FillHistogramWithVertErrors(interaction.Enu_mc_reco_signal, m_Enu * MeV_to_GeV);
            // MC Reco vs True -- Response
            FillHistogramWithVertErrors(interaction.Enu_response, m_Enu * MeV_to_GeV, m_Enu_Truth * MeV_to_GeV);
        }else{
            // MC Reco Background
            FillHistogramWithVertErrors(interaction.Enu_mc_reco_bckg, m_Enu * MeV_to_GeV);
        }
    }else{
        // Data
        FillHistogram(interaction.Enu_all, m_Enu * MeV_to_GeV);
    }
}

void CCProtonPi0_Analyzer::fill_QSq() 
{
    if (m_isMC){
        // MC Reco All
        FillHistogramWithVertErrors(interaction.QSq_mc_reco_all, m_QSq * MeVSq_to_GeVSq);
        if (truth_isSignal){
            // MC Truth Signal
            FillHistogramWithVertErrors(interaction.QSq_mc_truth_signal, m_QSq_Truth * MeVSq_to_GeVSq);
            // MC Reco Signal
            FillHistogramWithVertErrors(interaction.QSq_mc_reco_signal, m_QSq * MeVSq_to_GeVSq);
            // MC Reco vs True -- Response
            FillHistogramWithVertErrors(interaction.QSq_response, m_QSq * MeVSq_to_GeVSq, m_QSq_Truth * MeVSq_to_GeVSq);
        }else{
            // MC Reco Background
            FillHistogramWithVertErrors(interaction.QSq_mc_reco_bckg, m_QSq * MeVSq_to_GeVSq);
        }
    }else{
        // Data
        FillHistogram(interaction.QSq_all, m_QSq * MeVSq_to_GeVSq);
    }
}

void CCProtonPi0_Analyzer::fill_muon_P() 
{
    if (m_isMC){
        // MC Reco All
        FillHistogramWithVertErrors(muon.muon_P_mc_reco_all, muon_P * MeV_to_GeV);
        if (truth_isSignal){
            // MC Truth Signal
            FillHistogramWithVertErrors(muon.muon_P_mc_truth_signal, truth_muon_P * MeV_to_GeV);
            // MC Reco Signal
            FillHistogramWithVertErrors(muon.muon_P_mc_reco_signal, muon_P * MeV_to_GeV);
            // MC Reco vs True -- Response
            FillHistogramWithVertErrors(muon.muon_P_response, muon_P * MeV_to_GeV, truth_muon_P * MeV_to_GeV);
        }else{
            // MC Reco Background
            FillHistogramWithVertErrors(muon.muon_P_mc_reco_bckg, muon_P * MeV_to_GeV);
        }
    }else{
        // Data
        FillHistogram(muon.muon_P_all, muon_P * MeV_to_GeV);
    }
}

void CCProtonPi0_Analyzer::fill_muon_theta() 
{
    // Get Muon Theta From Nodes
    double reco_muon_theta = GetCorrectedMuonTheta();

    if (m_isMC){
        // MC Reco All
        FillHistogramWithVertErrors(muon.muon_theta_mc_reco_all, reco_muon_theta * TMath::RadToDeg());
        if (truth_isSignal){
            // MC Truth Signal
            FillHistogramWithVertErrors(muon.muon_theta_mc_truth_signal, truth_muon_theta_beam * TMath::RadToDeg());
            // MC Reco Signal
            FillHistogramWithVertErrors(muon.muon_theta_mc_reco_signal, reco_muon_theta * TMath::RadToDeg());
            // MC Reco vs True -- Response
            FillHistogramWithVertErrors(muon.muon_theta_response, reco_muon_theta * TMath::RadToDeg(), truth_muon_theta_beam * TMath::RadToDeg());
        }else{
            // MC Reco Background
            FillHistogramWithVertErrors(muon.muon_theta_mc_reco_bckg, reco_muon_theta * TMath::RadToDeg());
        }
    }else{
        // Data
        FillHistogram(muon.muon_theta_all, reco_muon_theta * TMath::RadToDeg());
    }
}

void CCProtonPi0_Analyzer::fill_pi0_P() 
{
    if (m_isMC){
        // MC Reco All
        FillHistogramWithVertErrors(pi0.pi0_P_mc_reco_all, pi0_P * MeV_to_GeV);
        if (truth_isSignal){
            // MC Truth Signal
            FillHistogramWithVertErrors(pi0.pi0_P_mc_truth_signal, truth_pi0_P * MeV_to_GeV);
            // MC Reco Signal
            FillHistogramWithVertErrors(pi0.pi0_P_mc_reco_signal, pi0_P * MeV_to_GeV);
            // MC Reco vs True -- Response
            FillHistogramWithVertErrors(pi0.pi0_P_response, pi0_P * MeV_to_GeV, truth_pi0_P * MeV_to_GeV);
        }else{
            // MC Reco Background
            FillHistogramWithVertErrors(pi0.pi0_P_mc_reco_bckg, pi0_P * MeV_to_GeV);
        }
    }else{
        // Data
        FillHistogram(pi0.pi0_P_all, pi0_P * MeV_to_GeV);
    }
}

void CCProtonPi0_Analyzer::fill_pi0_KE() 
{
    if (m_isMC){
        // MC Reco All
        FillHistogramWithVertErrors(pi0.pi0_KE_mc_reco_all, pi0_KE * MeV_to_GeV);
        if (truth_isSignal){
            // MC Truth Signal
            FillHistogramWithVertErrors(pi0.pi0_KE_mc_truth_signal, truth_pi0_KE  * MeV_to_GeV);
            // MC Reco Signal
            FillHistogramWithVertErrors(pi0.pi0_KE_mc_reco_signal, pi0_KE * MeV_to_GeV);
            // MC Reco vs True -- Response
            FillHistogramWithVertErrors(pi0.pi0_KE_response, pi0_KE * MeV_to_GeV, truth_pi0_KE * MeV_to_GeV);
        }else{
            // MC Reco Background
            FillHistogramWithVertErrors(pi0.pi0_KE_mc_reco_bckg, pi0_KE * MeV_to_GeV);
        }
    }else{
        // Data
        FillHistogram(pi0.pi0_KE_all, pi0_KE * MeV_to_GeV);
    }
}

void CCProtonPi0_Analyzer::fill_pi0_theta() 
{
    if (m_isMC){
        // MC Reco All
        FillHistogramWithVertErrors(pi0.pi0_theta_mc_reco_all, pi0_theta_beam * TMath::RadToDeg());
        if (truth_isSignal){
            // MC Truth Signal
            FillHistogramWithVertErrors(pi0.pi0_theta_mc_truth_signal, truth_pi0_theta_beam  * TMath::RadToDeg());
            // MC Reco Signal
            FillHistogramWithVertErrors(pi0.pi0_theta_mc_reco_signal, pi0_theta_beam * TMath::RadToDeg());
            // MC Reco vs True -- Response
            FillHistogramWithVertErrors(pi0.pi0_theta_response, pi0_theta_beam * TMath::RadToDeg(), truth_pi0_theta_beam * TMath::RadToDeg());
        }else{
            // MC Reco Background
            FillHistogramWithVertErrors(pi0.pi0_theta_mc_reco_bckg, pi0_theta_beam * TMath::RadToDeg());
        }
    }else{
        // Data
        FillHistogram(pi0.pi0_theta_all, pi0_theta_beam * TMath::RadToDeg());
    }
}

void CCProtonPi0_Analyzer::fill_deltaInvMass()
{
    if (m_isMC){
        // MC Reco All
        FillHistogramWithVertErrors(interaction.deltaInvMass_mc_reco_all, deltaInvMass_reco);
        if (truth_isSignal){
            // MC Truth Signal
            FillHistogramWithVertErrors(interaction.deltaInvMass_mc_truth_signal, deltaInvMass_true);
            // MC Reco Signal
            FillHistogramWithVertErrors(interaction.deltaInvMass_mc_reco_signal, deltaInvMass_reco);
            // MC Reco vs True -- Response
            FillHistogramWithVertErrors(interaction.deltaInvMass_response, deltaInvMass_reco, deltaInvMass_true);
        }else{
            // MC Reco Background
            FillHistogramWithVertErrors(interaction.deltaInvMass_mc_reco_bckg, deltaInvMass_reco);
        }
    }else{
        // Data
        FillHistogram(interaction.deltaInvMass_all, deltaInvMass_reco);
    }
}

void CCProtonPi0_Analyzer::fill_Delta_pi_theta()
{
    if (m_isMC){
        // MC Reco All
        FillHistogramWithVertErrors(interaction.Delta_pi_theta_mc_reco_all, Delta_pi_theta_reco);
        if (truth_isSignal){
            // MC Truth Signal
            FillHistogramWithVertErrors(interaction.Delta_pi_theta_mc_truth_signal, Delta_pi_theta_true);
            // MC Reco Signal
            FillHistogramWithVertErrors(interaction.Delta_pi_theta_mc_reco_signal, Delta_pi_theta_reco);
            // MC Reco vs True -- Response
            FillHistogramWithVertErrors(interaction.Delta_pi_theta_response, Delta_pi_theta_reco, Delta_pi_theta_true);
        }else{
            // MC Reco Background
            FillHistogramWithVertErrors(interaction.Delta_pi_theta_mc_reco_bckg, Delta_pi_theta_reco);
        }
    }else{
        // Data
        FillHistogram(interaction.Delta_pi_theta_all, Delta_pi_theta_reco);
    }
}

void CCProtonPi0_Analyzer::fill_Delta_pi_phi()
{
    if (m_isMC){
        // MC Reco All
        FillHistogramWithVertErrors(interaction.Delta_pi_phi_mc_reco_all, Delta_pi_phi_reco);
        if (truth_isSignal){
            // MC Truth Signal
            FillHistogramWithVertErrors(interaction.Delta_pi_phi_mc_truth_signal, Delta_pi_phi_true);
            // MC Reco Signal
            FillHistogramWithVertErrors(interaction.Delta_pi_phi_mc_reco_signal, Delta_pi_phi_reco);
            // MC Reco vs True -- Response
            FillHistogramWithVertErrors(interaction.Delta_pi_phi_response, Delta_pi_phi_reco, Delta_pi_phi_true);
        }else{
            // MC Reco Background 
            FillHistogramWithVertErrors(interaction.Delta_pi_phi_mc_reco_bckg, Delta_pi_phi_reco);
        }
    }else{
        // Data
        FillHistogram(interaction.Delta_pi_phi_all, Delta_pi_phi_reco);
    }
}

void CCProtonPi0_Analyzer::fill_q3_q0()
{
    double q0_reco = (m_Enu - muon_E) * MeV_to_GeV;
    double q3_reco = Calc_q3(false) * MeV_to_GeV; // Reco Value

    if (m_isMC){
        // MC Reco All
        FillHistogramWithLeadingErrors(interaction.q3_q0_mc_reco_all, q3_reco, q0_reco);
        if (truth_isSignal){
            double q0_true = (m_Enu_Truth - truth_muon_4P[3]) * MeV_to_GeV;
            double q3_true = Calc_q3(true) * MeV_to_GeV; // True Value
            // MC Truth Signal
            FillHistogramWithLeadingErrors(interaction.q3_q0_mc_truth_signal, q3_true, q0_true);
            // MC Reco Signal
            FillHistogramWithLeadingErrors(interaction.q3_q0_mc_reco_signal, q3_reco, q0_reco);
        }else{
            // MC Reco Background
            FillHistogramWithLeadingErrors(interaction.q3_q0_mc_reco_bckg, q3_reco, q0_reco);
        }
    }else{
        // Data
        FillHistogram(interaction.q3_q0_all, q3_reco, q0_reco);
    }
}

void CCProtonPi0_Analyzer::fill_muon_theta_muon_KE()
{
    double muon_theta_reco = cos( GetCorrectedMuonTheta() );
    double muon_KE_reco = muon_KE * MeV_to_GeV;

    if (m_isMC){
        // MC Reco All
        FillHistogramWithLeadingErrors(interaction.muon_theta_muon_KE_mc_reco_all, muon_theta_reco, muon_KE_reco);
        if (truth_isSignal){
            double muon_theta_true = cos(truth_muon_theta_beam);
            double muon_KE_true = (truth_muon_4P[3] - muon_mass) * MeV_to_GeV;
            // MC Truth Signal
            FillHistogramWithLeadingErrors(interaction.muon_theta_muon_KE_mc_truth_signal, muon_theta_true, muon_KE_true);
            // MC Reco Signal
            FillHistogramWithLeadingErrors(interaction.muon_theta_muon_KE_mc_reco_signal, muon_theta_reco, muon_KE_reco);
        }else{
            // MC Reco Background
            FillHistogramWithLeadingErrors(interaction.muon_theta_muon_KE_mc_reco_bckg, muon_theta_reco, muon_KE_reco);
        }
    }else{
        // Data
        FillHistogram(interaction.muon_theta_muon_KE_all, muon_theta_reco, muon_KE_reco);
    }
}

void CCProtonPi0_Analyzer::fillMuonReco()
{
    // Get Muon Theta From Nodes
    double reco_muon_theta = GetCorrectedMuonTheta();

    FillHistogram(muon.E, muon_E * MeV_to_GeV);
    FillHistogram(muon.P, muon_P * MeV_to_GeV);
    FillHistogram(muon.KE, muon_KE * MeV_to_GeV);
    FillHistogram(muon.theta, reco_muon_theta * TMath::RadToDeg());
    FillHistogram(muon.cos_theta, cos(reco_muon_theta));
    FillHistogram(muon.phi, muon_phi * TMath::RadToDeg());
}

void CCProtonPi0_Analyzer::FillHistogram(TH1D* hist, double var)
{
    hist->Fill(var, cvweight);
}

void CCProtonPi0_Analyzer::FillHistogram(TH2D* hist, double var1, double var2)
{
    hist->Fill(var1,var2, cvweight);
}

void CCProtonPi0_Analyzer::FillHistogram(TH3D* hist, double var1, double var2, double var3)
{
    hist->Fill(var1, var2, var3, cvweight);
}

void CCProtonPi0_Analyzer::FillHistogram(MnvH1D* hist, double var)
{
    hist->Fill(var, cvweight);
}

void CCProtonPi0_Analyzer::FillHistogram(MnvH2D* hist, double var1, double var2)
{
    hist->Fill(var1,var2, cvweight);
}

void CCProtonPi0_Analyzer::FillHistogram(vector<MnvH1D*> &hist, double var)
{
    // Always Fill hist[0]
    hist[0]->Fill(var, cvweight);

    // Fill others only if Analyzing MC
    if (m_isMC){
        // Fill Signal
        if (truth_isSignal){
            hist[1]->Fill(var, cvweight);
 
            int ind = GetSignalTypeInd();
            hist[ind]->Fill(var,cvweight);
        }else{
            // Fill Background
            hist[2]->Fill(var, cvweight); // Always Fill ind == 2 -- All Background

            // Fill Background Type
            int ind = GetBackgroundTypeInd();
            hist[ind]->Fill(var, cvweight);
        }
    }
}

void CCProtonPi0_Analyzer::FillHistogram(vector<MnvH1D*> &hist, double var, double wgt)
{
    // Always Fill hist[0]
    hist[0]->Fill(var, wgt);

    // Fill others only if Analyzing MC
    if (m_isMC){
        // Fill Signal
        if (truth_isSignal){
            hist[1]->Fill(var, wgt);

            int ind = GetSignalTypeInd();
            hist[ind]->Fill(var, wgt);
        }else{
            // Fill Background
            hist[2]->Fill(var, wgt); // Always Fill ind == 2 -- All Background

            // Fill Background Type
            int ind = GetBackgroundTypeInd();
            hist[ind]->Fill(var, wgt);
        }
    }
}

void CCProtonPi0_Analyzer::FillHistogram(vector<MnvH2D*> &hist, double var1, double var2)
{
    // Always Fill hist[0]
    hist[0]->Fill(var1, var2, cvweight);

    // Fill others only if Analyzing MC
    if (m_isMC){
        // Fill Signal
        if (truth_isSignal){
            hist[1]->Fill(var1, var2, cvweight);
 
            int ind = GetSignalTypeInd();
            hist[ind]->Fill(var1, var2, cvweight);
        }else{
            // Fill Background
            hist[2]->Fill(var1, var2,  cvweight); // Always Fill ind == 2 -- All Background

            // Fill Background Type
            int ind = GetBackgroundTypeInd();
            hist[ind]->Fill(var1, var2, cvweight);
        }
    }
}

int CCProtonPi0_Analyzer::GetBackgroundTypeInd()
{
    // Check For Signal
    if (truth_isSignal){
        cout<<"WARNING! Signal Event requested Background Ind! - Returning -1"<<endl;
        return -1;
    }

    if (truth_isBckg_Compact_WithPi0) return 3;
    else if (truth_isBckg_Compact_QELike) return 4;
    else if (truth_isBckg_Compact_SinglePiPlus) return 5;
    else if (truth_isBckg_Compact_Other) return 6;
    else{
        cout<<"WARNING! No Background Type Found - Returning -1"<<endl;
        return -1;
    }
}

int CCProtonPi0_Analyzer::GetSignalTypeInd()
{
    // Check for Background
    if (!truth_isSignal){
        RunTimeError("WARNING! Background Event requested Signal Ind! -- Exiting!");
    }

    int ind = -1;
    if (mc_intType == 2){
        if (mc_resID == 0){
            //Delta Resonance
            ind = 7;
        }else{
            // Other Resonance
            ind = 8;
        }
    }else{
        // Non-Resonance
        ind = 9;
    }

   return ind; 
}

int CCProtonPi0_Analyzer::GetEjectedNucleonCount()
{
    int n_nucleons = 0;
    for(int i = 0; i < mc_nFSPart; i++ ){
        if (mc_FSPartPDG[i] == PDG_List::proton || mc_FSPartPDG[i] == PDG_List::neutron){
            n_nucleons++;
        }
    }

    return n_nucleons;
}

// Function Reserved for Correcting a NTupleVariables
void CCProtonPi0_Analyzer::CorrectNTupleVariables()
{
    // Do Nothing
}

void CCProtonPi0_Analyzer::CorrectEMShowerCalibration()
{
    // Pi0 Variables
    pi0_invMass *= EM_correction;
    pi0_E_Cal *= EM_correction;
    pi0_P *= EM_correction;
    pi0_px *= EM_correction;
    pi0_py *= EM_correction;
    pi0_pz *= EM_correction;
    pi0_E = sqrt(pi0_P*pi0_P + pi0_mass*pi0_mass);
    pi0_KE = pi0_E - pi0_mass;

    // gamma1 
    gamma1_E *= EM_correction;
    gamma1_E_Old *= EM_correction;
    gamma1_P *= EM_correction;
    gamma1_px *= EM_correction;
    gamma1_py *= EM_correction;
    gamma1_pz *= EM_correction;

    // gamma2 
    gamma2_E *= EM_correction;
    gamma2_E_Old *= EM_correction;
    gamma2_P *= EM_correction;
    gamma2_px *= EM_correction;
    gamma2_py *= EM_correction;
    gamma2_pz *= EM_correction;
}

void CCProtonPi0_Analyzer::CalcEventWeight()
{
    if (m_isMC){

        // Reset cvweight
        cvweight = 1.0;
        cvweight_Flux = 1.0;
        cvweight_Delta = 1.0;
        cvweight_CCRES = 1.0;
        cvweight_NonRes1pi = 1.0;
        cvweight_2p2h = 1.0;

        UpdateFluxReweighter(mc_run, mc_intType); 

        // Replace cvweight with Flux Weight
        cvweight_Flux = GetFluxWeight(mc_incomingE * MeV_to_GeV, mc_incoming);
        cvweight *= cvweight_Flux;

        // Update cvweight with MINOS Efficiency Correction
        const double minos_eff_correction = GetMINOSCorrection();
        cvweight *= minos_eff_correction;

        // 2p2h Event Weights
        if ( IsEvent2p2h(mc_intType) ){
            cvweight_2p2h = Get_2p2h_wgt(mc_incomingPartVec, mc_primFSLepton, fit_2p2h_CV) * POT_ratio_2p2h;
            cvweight *= cvweight_2p2h;
            // Calc 2p2h Uncertainty 
            Get2p2hErr();
        }else{
            // 2p2h Uncertainty is 0.0 for other events 
            err_2p2h = 0.0;
        }

        // Apply Background Constraints
        if (applyBckgConstraints_CV){
            if (truth_isBckg_Compact_SinglePiPlus) cvweight *= cv_wgt_SinglePiPlus;
            else if (truth_isBckg_Compact_QELike) cvweight *= cv_wgt_QELike;
            else if (truth_isBckg_Compact_WithPi0) cvweight *= cv_wgt_WithPi0;
        }
            
        if (applyGENIETuning_Complete){

            // Delta decay anisotropy weighting per DocDB 9850.  Weight is 1.0 for non-delta resonance interactions.
            cvweight_Delta *= ( 1.0 + truth_genie_wgt_Theta_Delta2Npi[4] ) / 2.0;
            cvweight *= cvweight_Delta;

            if ( IsGenieRvn1pi() || IsGenieRvp1pi() ){
                cvweight_NonRes1pi *= deuteriumNonResNorm;
                cvweight *= cvweight_NonRes1pi;
            }
        
            if ( IsGenieCCRes() ){
                cvweight_CCRES *= deuteriumResNorm; 
                cvweight_CCRES *= GetMaResWeight(deuteriumResNorm);
                cvweight *= cvweight_CCRES;
            }
        }

        if ( applyGENIETuning_DeltaSuppression && IsGenieCCRes() ){
            // Delta Suppression Factor
            double QSq = truth_QSq_exp * MeVSq_to_GeVSq;
            cvweight_CCRES *= GetDeltaFactor(QSq, DeltaFactor_A, DeltaFactor_Q0);
            cvweight *= cvweight_CCRES;
        }

        UpdateGENIESystematics();
    }else{
        cvweight = 1.0; 
    }
}

void CCProtonPi0_Analyzer::FillInvMass_TruthMatch()
{
    int pdg;
    if (truth_blob1_evis_total_truth > truth_blob2_evis_total_truth){
        pdg = truth_blob1_evis_most_pdg;
    }else{
        pdg = truth_blob2_evis_most_pdg;
    }

    if (truth_isSignal){
        if (pdg == 111) FillHistogram(cutList.signal_invMass_pizero, pi0_invMass);
        else if (pdg == 211) FillHistogram(cutList.signal_invMass_piplus, pi0_invMass);
        else if (pdg == 2212) FillHistogram(cutList.signal_invMass_proton, pi0_invMass);
        else if (pdg == 2112) FillHistogram(cutList.signal_invMass_neutron, pi0_invMass);
        else FillHistogram(cutList.signal_invMass_other, pi0_invMass);
    }else{
        if (pdg == 111) FillHistogram(cutList.background_invMass_pizero, pi0_invMass);
        else if (pdg == 211) FillHistogram(cutList.background_invMass_piplus, pi0_invMass);
        else if (pdg == 2212) FillHistogram(cutList.background_invMass_proton, pi0_invMass);
        else if (pdg == 2112) FillHistogram(cutList.background_invMass_neutron, pi0_invMass);
        else FillHistogram(cutList.background_invMass_other, pi0_invMass);
    }
}


double CCProtonPi0_Analyzer::GetMINOSCorrectionErr()
{
    std::string playlist = GetPlaylist(mc_run, mc_intType);
    if (playlist.compare("minerva_2p2h") == 0) playlist = "minerva13C";

    MnvNormalizer normalizer("Eroica", playlist);
    double correctionErr = normalizer.GetCorrectionErr(CCProtonPi0_minos_trk_p);
    return correctionErr;
}

double CCProtonPi0_Analyzer::GetMINOSCorrection()
{
    std::string playlist = GetPlaylist(mc_run, mc_intType);
    if (playlist.compare("minerva_2p2h") == 0) playlist = "minerva13C";

    MnvNormalizer normalizer("Eroica", playlist);
    double correction = normalizer.GetCorrection(CCProtonPi0_minos_trk_p);
    return correction;
}

double CCProtonPi0_Analyzer::Calc_MuonCosTheta()
{
    double P_muon = truth_muon_P;
    double P_beam = HEP_Functions::calcMomentum(mc_incomingPartVec[0],mc_incomingPartVec[1],mc_incomingPartVec[2]);
    double dot_product = 0.0;

    for (int i = 0; i < 3; ++i){
        dot_product += truth_muon_4P[i]*mc_incomingPartVec[i];
    }


    double cos_theta = dot_product / (P_muon * P_beam);
    return cos_theta;
}

double CCProtonPi0_Analyzer::Calc_Enu()
{
    double Enu;
    if (nProtonCandidates == 0){
        Enu = muon_E + pi0_E + vertex_blob_energy + Extra_Energy_Total;
    }else{
        Enu = muon_E + m_total_proton_KE + pi0_E + vertex_blob_energy + Extra_Energy_Total;
    }

    return Enu;
}

double CCProtonPi0_Analyzer::Calc_Enu_Cal()
{
    double Enu = muon_E + CCProtonPi0_hadron_recoil_CCInc;

    return Enu;
}

double CCProtonPi0_Analyzer::Calc_Enu_shifted(double muon_E_shifted, double pi0_E_shifted, double total_proton_KE_shifted)
{
    double Enu_shifted;
    if (nProtonCandidates == 0){
        Enu_shifted = muon_E_shifted + pi0_E_shifted + vertex_blob_energy + Extra_Energy_Total;
    }else{
        Enu_shifted = muon_E_shifted + total_proton_KE_shifted + pi0_E_shifted + vertex_blob_energy + Extra_Energy_Total;
    }

    return Enu_shifted;
}

TLorentzVector CCProtonPi0_Analyzer::Get_Neutrino_4P(const double Enu) const
{
    const double theta = -0.05887;

    double Px = 0.0;
    double Py = Enu*sin(theta);
    double Pz = Enu*cos(theta); 

    TLorentzVector beam_4P(Px,Py,Pz,Enu);

    return beam_4P;
}

void CCProtonPi0_Analyzer::Calc_EventKinematics()
{
    // Get Corrected Muon Angle
    double reco_muon_theta = GetCorrectedMuonTheta();

    // Get Total Proton KE
    m_total_proton_KE = 0.0;
    for (int i = 0; i < nProtonCandidates; ++i ){
        m_total_proton_KE += all_protons_KE[i];
    }

    m_Enu = Calc_Enu();
    //m_Enu = Calc_Enu_Cal(); // Test results are not good! -- Use Calc_Enu()
    m_QSq = Calc_QSq(m_Enu, muon_E, muon_P, reco_muon_theta);
    m_WSq = Calc_WSq(m_Enu, m_QSq, muon_E);
    m_W = m_WSq > 0 ? sqrt(m_WSq) : -1.0;
}

void CCProtonPi0_Analyzer::Calc_EventKinematics_Truth()
{
    m_Enu_Truth = mc_incomingE;
    m_QSq_Truth = truth_QSq_exp;
    m_WSq_Truth = truth_WSq_exp;
    m_W_Truth = truth_W_exp > 0.0 ? truth_W_exp : -1.0;
}

void CCProtonPi0_Analyzer::FillBackgroundCharacteristics()
{
    if (mc_intType == 1){
        FillHistogram(interaction.reco_bckg_QSq_QE, m_QSq * MeVSq_to_GeVSq);
        FillHistogram(interaction.reco_bckg_Enu_QE, m_Enu * MeV_to_GeV);
        FillHistogram(interaction.reco_bckg_w_QE, m_W * MeV_to_GeV);
    }else if (mc_intType == 2){
        if (mc_resID == 0){ 
            FillHistogram(interaction.reco_bckg_QSq_RES_1232, m_QSq * MeVSq_to_GeVSq);
            FillHistogram(interaction.reco_bckg_Enu_RES_1232, m_Enu * MeV_to_GeV);
            FillHistogram(interaction.reco_bckg_w_RES_1232, m_W * MeV_to_GeV);
        }else if (mc_resID == 1){
            FillHistogram(interaction.reco_bckg_QSq_RES_1535, m_QSq * MeVSq_to_GeVSq);
            FillHistogram(interaction.reco_bckg_Enu_RES_1535, m_Enu * MeV_to_GeV);
            FillHistogram(interaction.reco_bckg_w_RES_1535, m_W * MeV_to_GeV);
        }else if (mc_resID == 2){
            FillHistogram(interaction.reco_bckg_QSq_RES_1520, m_QSq * MeVSq_to_GeVSq);
            FillHistogram(interaction.reco_bckg_Enu_RES_1520, m_Enu * MeV_to_GeV);
            FillHistogram(interaction.reco_bckg_w_RES_1520, m_W * MeV_to_GeV);
        }else{
            FillHistogram(interaction.reco_bckg_QSq_RES_Other, m_QSq * MeVSq_to_GeVSq);
            FillHistogram(interaction.reco_bckg_Enu_RES_Other, m_Enu * MeV_to_GeV);
            FillHistogram(interaction.reco_bckg_w_RES_Other, m_W * MeV_to_GeV);
        }
    }else if ( IsGenieNonRES() ){
        FillHistogram(interaction.reco_bckg_QSq_Non_RES, m_QSq * MeVSq_to_GeVSq);
        FillHistogram(interaction.reco_bckg_Enu_Non_RES, m_Enu * MeV_to_GeV);
        FillHistogram(interaction.reco_bckg_w_Non_RES, m_W * MeV_to_GeV);
    }else if (mc_intType == 3){
        FillHistogram(interaction.reco_bckg_QSq_DIS, m_QSq * MeVSq_to_GeVSq);
        FillHistogram(interaction.reco_bckg_Enu_DIS, m_Enu * MeV_to_GeV);
        FillHistogram(interaction.reco_bckg_w_DIS, m_W * MeV_to_GeV);
    }else if (mc_intType == 4){
        FillHistogram(interaction.reco_bckg_QSq_Coh, m_QSq * MeVSq_to_GeVSq);
        FillHistogram(interaction.reco_bckg_Enu_Coh, m_Enu * MeV_to_GeV);
        FillHistogram(interaction.reco_bckg_w_Coh, m_W * MeV_to_GeV);
    }else if ( IsEvent2p2h(mc_intType) ){
        FillHistogram(interaction.reco_bckg_QSq_2p2h, m_QSq * MeVSq_to_GeVSq);
        FillHistogram(interaction.reco_bckg_Enu_2p2h, m_Enu * MeV_to_GeV);
        FillHistogram(interaction.reco_bckg_w_2p2h, m_W * MeV_to_GeV);
    }else{
        std::cout<<"WARNING! Background Event with different interaction Type!"<<std::endl;
        std::cout<<"mc_intType = "<<mc_intType<<std::endl;
    }
}
void CCProtonPi0_Analyzer::FillSignalCharacteristics_Reco()
{
    if (mc_intType == 1){
        FillHistogram(interaction.reco_w_QE, m_W * MeV_to_GeV);
    }else if (mc_intType == 2){
        if (mc_resID == 0){ 
            FillHistogram(interaction.reco_w_RES_1232, m_W * MeV_to_GeV);
        }else if (mc_resID == 1){
            FillHistogram(interaction.reco_w_RES_1535, m_W * MeV_to_GeV);
        }else if (mc_resID == 2){
            FillHistogram(interaction.reco_w_RES_1520, m_W * MeV_to_GeV);
        }else{
            FillHistogram(interaction.reco_w_RES_Other, m_W * MeV_to_GeV);
        }
    }else if ( IsGenieNonRES() ){
        FillHistogram(interaction.reco_w_Non_RES, m_W * MeV_to_GeV);
    }else if (mc_intType == 3){
        FillHistogram(interaction.reco_w_DIS, m_W * MeV_to_GeV);
    }else if ( IsEvent2p2h(mc_intType) ){
        FillHistogram(interaction.reco_w_2p2h, m_W * MeV_to_GeV);
    }else{
        std::cout<<"WARNING! Signal Event with different interaction Type!"<<std::endl;
    }
}

void CCProtonPi0_Analyzer::FillSignalCharacteristics(bool isMinosMatched)
{
    if (isMinosMatched) {
        // Signal Characteristics
        if (mc_intType == 1){
            FillHistogram(cutList.truth_QSq_QE, m_QSq_Truth * MeVSq_to_GeVSq);
            FillHistogram(cutList.truth_Enu_QE, m_Enu_Truth * MeV_to_GeV);
            FillHistogram(cutList.truth_w_QE, m_W_Truth * MeV_to_GeV);
        }else if (mc_intType == 2){
            if (mc_resID == 0){ 
                FillHistogram(cutList.truth_QSq_RES_1232, m_QSq_Truth * MeVSq_to_GeVSq);
                FillHistogram(cutList.truth_Enu_RES_1232, m_Enu_Truth * MeV_to_GeV);
                FillHistogram(cutList.truth_w_RES_1232, m_W_Truth * MeV_to_GeV);
            }else if (mc_resID == 1){
                FillHistogram(cutList.truth_QSq_RES_1535, m_QSq_Truth * MeVSq_to_GeVSq);
                FillHistogram(cutList.truth_Enu_RES_1535, m_Enu_Truth * MeV_to_GeV);
                FillHistogram(cutList.truth_w_RES_1535, m_W_Truth * MeV_to_GeV);
            }else if (mc_resID == 2){
                FillHistogram(cutList.truth_QSq_RES_1520, m_QSq_Truth * MeVSq_to_GeVSq);
                FillHistogram(cutList.truth_Enu_RES_1520, m_Enu_Truth * MeV_to_GeV);
                FillHistogram(cutList.truth_w_RES_1520, m_W_Truth * MeV_to_GeV);
            }else{
                FillHistogram(cutList.truth_QSq_RES_Other, m_QSq_Truth * MeVSq_to_GeVSq);
                FillHistogram(cutList.truth_Enu_RES_Other, m_Enu_Truth * MeV_to_GeV);
                FillHistogram(cutList.truth_w_RES_Other, m_W_Truth * MeV_to_GeV);
            }
        }else if ( IsGenieNonRES() ){
            FillHistogram(cutList.truth_QSq_Non_RES, m_QSq_Truth * MeVSq_to_GeVSq);
            FillHistogram(cutList.truth_Enu_Non_RES, m_Enu_Truth * MeV_to_GeV);
            FillHistogram(cutList.truth_w_Non_RES, m_W_Truth * MeV_to_GeV);
        }else if (mc_intType == 3){
            FillHistogram(cutList.truth_QSq_DIS, m_QSq_Truth * MeVSq_to_GeVSq);
            FillHistogram(cutList.truth_Enu_DIS, m_Enu_Truth * MeV_to_GeV);
            FillHistogram(cutList.truth_w_DIS, m_W_Truth * MeV_to_GeV);
        }else if ( IsEvent2p2h(mc_intType) ){
            FillHistogram(cutList.truth_QSq_2p2h, m_QSq_Truth * MeVSq_to_GeVSq);
            FillHistogram(cutList.truth_Enu_2p2h, m_Enu_Truth * MeV_to_GeV);
            FillHistogram(cutList.truth_w_2p2h, m_W_Truth * MeV_to_GeV);
        }else{
            std::cout<<"WARNING! Signal Event with different interaction Type!"<<std::endl;
        }
    }else{
        // Signal Characteristics
        if (mc_intType == 1){
            FillHistogram(interaction.mc_Q2_QE, mc_Q2 * MeVSq_to_GeVSq);
            FillHistogram(interaction.mc_incomingE_QE, mc_incomingE * MeV_to_GeV);
            FillHistogram(interaction.mc_w_QE, mc_w * MeV_to_GeV);

            FillHistogram(interaction.truth_QSq_QE, m_QSq_Truth * MeVSq_to_GeVSq);
            FillHistogram(interaction.truth_Enu_QE, m_Enu_Truth * MeV_to_GeV);
            FillHistogram(interaction.truth_w_QE, m_W_Truth * MeV_to_GeV);
        }else if (mc_intType == 2){
            if (mc_resID == 0){ 
                FillHistogram(interaction.mc_Q2_RES_1232, mc_Q2 * MeVSq_to_GeVSq);
                FillHistogram(interaction.mc_incomingE_RES_1232, mc_incomingE * MeV_to_GeV);
                FillHistogram(interaction.mc_w_RES_1232, mc_w * MeV_to_GeV);

                FillHistogram(interaction.truth_QSq_RES_1232, m_QSq_Truth * MeVSq_to_GeVSq);
                FillHistogram(interaction.truth_Enu_RES_1232, m_Enu_Truth * MeV_to_GeV);
                FillHistogram(interaction.truth_w_RES_1232, m_W_Truth * MeV_to_GeV);
            }else if (mc_resID == 1){
                FillHistogram(interaction.mc_Q2_RES_1535, mc_Q2 * MeVSq_to_GeVSq);
                FillHistogram(interaction.mc_incomingE_RES_1535, mc_incomingE * MeV_to_GeV);
                FillHistogram(interaction.mc_w_RES_1535, mc_w * MeV_to_GeV);

                FillHistogram(interaction.truth_QSq_RES_1535, m_QSq_Truth * MeVSq_to_GeVSq);
                FillHistogram(interaction.truth_Enu_RES_1535, m_Enu_Truth * MeV_to_GeV);
                FillHistogram(interaction.truth_w_RES_1535, m_W_Truth * MeV_to_GeV);
            }else if (mc_resID == 2){
                FillHistogram(interaction.mc_Q2_RES_1520, mc_Q2 * MeVSq_to_GeVSq);
                FillHistogram(interaction.mc_incomingE_RES_1520, mc_incomingE * MeV_to_GeV);
                FillHistogram(interaction.mc_w_RES_1520, mc_w * MeV_to_GeV);

                FillHistogram(interaction.truth_QSq_RES_1520, m_QSq_Truth * MeVSq_to_GeVSq);
                FillHistogram(interaction.truth_Enu_RES_1520, m_Enu_Truth * MeV_to_GeV);
                FillHistogram(interaction.truth_w_RES_1520, m_W_Truth * MeV_to_GeV);
            }else{
                FillHistogram(interaction.mc_Q2_RES_Other, mc_Q2 * MeVSq_to_GeVSq);
                FillHistogram(interaction.mc_incomingE_RES_Other, mc_incomingE * MeV_to_GeV);
                FillHistogram(interaction.mc_w_RES_Other, mc_w * MeV_to_GeV);

                FillHistogram(interaction.truth_QSq_RES_Other, m_QSq_Truth * MeVSq_to_GeVSq);
                FillHistogram(interaction.truth_Enu_RES_Other, m_Enu_Truth * MeV_to_GeV);
                FillHistogram(interaction.truth_w_RES_Other, m_W_Truth * MeV_to_GeV);
            }
        }else if ( IsGenieNonRES() ){
            FillHistogram(interaction.mc_Q2_Non_RES, mc_Q2 * MeVSq_to_GeVSq);
            FillHistogram(interaction.mc_incomingE_Non_RES, mc_incomingE * MeV_to_GeV);
            FillHistogram(interaction.mc_w_Non_RES, mc_w * MeV_to_GeV);

            FillHistogram(interaction.truth_QSq_Non_RES, m_QSq_Truth * MeVSq_to_GeVSq);
            FillHistogram(interaction.truth_Enu_Non_RES, m_Enu_Truth * MeV_to_GeV);
            FillHistogram(interaction.truth_w_Non_RES, m_W_Truth * MeV_to_GeV);
        }else if (mc_intType == 3){
            FillHistogram(interaction.mc_Q2_DIS, mc_Q2 * MeVSq_to_GeVSq);
            FillHistogram(interaction.mc_incomingE_DIS, mc_incomingE * MeV_to_GeV);
            FillHistogram(interaction.mc_w_DIS, mc_w * MeV_to_GeV);

            FillHistogram(interaction.truth_QSq_DIS, m_QSq_Truth * MeVSq_to_GeVSq);
            FillHistogram(interaction.truth_Enu_DIS, m_Enu_Truth * MeV_to_GeV);
            FillHistogram(interaction.truth_w_DIS, m_W_Truth * MeV_to_GeV);
        }else if ( IsEvent2p2h(mc_intType) ){
            FillHistogram(interaction.mc_Q2_2p2h, mc_Q2 * MeVSq_to_GeVSq);
            FillHistogram(interaction.mc_incomingE_2p2h, mc_incomingE * MeV_to_GeV);
            FillHistogram(interaction.mc_w_2p2h, mc_w * MeV_to_GeV);

            FillHistogram(interaction.truth_QSq_2p2h, m_QSq_Truth * MeVSq_to_GeVSq);
            FillHistogram(interaction.truth_Enu_2p2h, m_Enu_Truth * MeV_to_GeV);
            FillHistogram(interaction.truth_w_2p2h, m_W_Truth * MeV_to_GeV);
        }else{
            std::cout<<"WARNING! Signal Event with different interaction Type!"<<std::endl;
        }
    }
}

int CCProtonPi0_Analyzer::Get_nFS_pions()
{
    int nFS_pions = 0;

    for (int i = 0; i < mc_er_nPart; ++i){
        if (std::abs(mc_er_ID[i]) == 211 || mc_er_ID[i] == 111){
            if (isMother_DIS_Fragment(i)) nFS_pions++;
        }
    }

    return nFS_pions;
}

bool CCProtonPi0_Analyzer::isMother_DIS_Fragment(int ind)
{
    int mother_ind = mc_er_mother[ind];

    if (mc_er_status[mother_ind] == 12) return true;
    else return false;
}

void CCProtonPi0_Analyzer::fill_SideBand_InvMass()
{
    bool fill_Michel = sideBand_Michel &&  (isMichelEvent || isShower_Michel_Exist );
    bool fill_pID = sideBand_PID && isPionTrack;

    if (fill_Michel || fill_pID){
        FillInvMass_TruthMatch();
        FillHistogramWithVertErrors(cutList.hCut_pi0invMass, pi0_invMass);
        // Fill Lateral Error Bands
        if (m_isMC){
            FillLatErrorBand_EM_EnergyScale_SideBand_invMass();
            FillLatErrorBand_MuonMomentum_SideBand_invMass();
            FillLatErrorBand_MuonTheta_SideBand_invMass();
            FillLatErrorBand_ProtonEnergy_Birks_SideBand_invMass();
            FillLatErrorBand_ProtonEnergy_SideBand_invMass("ProtonEnergy_MassModel");
            FillLatErrorBand_ProtonEnergy_SideBand_invMass("ProtonEnergy_MEU");
            FillLatErrorBand_ProtonEnergy_SideBand_invMass("ProtonEnergy_BetheBloch");
        }
        if (nProtonCandidates == 0){
            FillHistogram(cutList.pi0_invMass_1Track, pi0_invMass);
            FillHistogram(cutList.hCut_1Track_pi0invMass,pi0_invMass);
        }else{ 
            FillHistogram(cutList.pi0_invMass_2Track, pi0_invMass);
            FillHistogram(cutList.hCut_2Track_pi0invMass, pi0_invMass);
        }
    }
}

void CCProtonPi0_Analyzer::fill_SideBand_Other()
{
    bool fill_Michel = sideBand_Michel &&  (isMichelEvent || isShower_Michel_Exist );
    bool fill_pID = sideBand_PID && isPionTrack;
    bool fill_LowInvMass = sideBand_LowInvMass && isLowInvMassEvent;
    bool fill_HighInvMass = sideBand_HighInvMass && isHighInvMassEvent;

    if (fill_Michel || fill_pID || fill_LowInvMass || fill_HighInvMass){
        double reco_theta = GetCorrectedMuonTheta();
        FillHistogramWithVertErrors(cutList.SideBand_muon_P, muon_P*MeV_to_GeV);
        FillHistogramWithVertErrors(cutList.SideBand_muon_theta, reco_theta*TMath::RadToDeg());
        FillHistogramWithVertErrors(cutList.SideBand_pi0_P, pi0_P*MeV_to_GeV);
        FillHistogramWithVertErrors(cutList.SideBand_pi0_KE, pi0_KE*MeV_to_GeV);
        FillHistogramWithVertErrors(cutList.SideBand_pi0_theta, pi0_theta_beam*TMath::RadToDeg());
        FillHistogramWithVertErrors(cutList.SideBand_neutrino_E, m_Enu*MeV_to_GeV);
        FillHistogramWithVertErrors(cutList.SideBand_QSq, m_QSq*MeVSq_to_GeVSq);
        FillHistogramWithVertErrors(cutList.SideBand_W, m_W*MeV_to_GeV);
        // I am not filling Lateral Error Bands on SideBand Other Variables
    }
}

double CCProtonPi0_Analyzer::GetCorrectedMuonTheta()
{
    double corrected_theta;

    if (muon_theta_allNodes_sz >= 28) corrected_theta = muon_theta_allNodes[28];
    else if (muon_theta_allNodes_sz >= 19) corrected_theta = muon_theta_allNodes[19];
    else if (muon_theta_allNodes_sz >= 9) corrected_theta = muon_theta_allNodes[9];
    else corrected_theta = muon_theta_beam;

    return corrected_theta;
}

std::vector<double> CCProtonPi0_Analyzer::Calc_Delta_pi_angles(TLorentzVector proton_4P, TLorentzVector pi0_4P, TLorentzVector muon_4P, TLorentzVector beam_4P)
{
    bool isDebug = false;

    using namespace std;

    // --------------------------------------------------------------------
    // Find Boost from LAB Frame to Delta Rest Frame
    // --------------------------------------------------------------------
    TLorentzVector delta_4P = proton_4P + pi0_4P;
    double M_delta = delta_4P.M(); 
    if (isDebug) cout<<"Delta Mass = "<<M_delta * MeV_to_GeV<<endl;

    double gamma = delta_4P.Gamma();
    double boost_x = -(delta_4P.Px()/ (gamma*M_delta));
    double boost_y = -(delta_4P.Py() / (gamma*M_delta));
    double boost_z = -(delta_4P.Pz() / (gamma*M_delta));

    if (isDebug){
        cout<<"Delta P (Before) = "<<delta_4P.Vect().Mag()<<endl;
        cout<<"Delta P (Before) = "<<delta_4P.Vect().x();
        cout<<" "<<delta_4P.Vect().y();
        cout<<" "<<delta_4P.Vect().z()<<endl;
    }

    delta_4P.Boost(boost_x,boost_y,boost_z);

    if (isDebug){
        cout<<"Delta P (After) = "<<delta_4P.Vect().Mag()<<std::endl;
        cout<<"Delta P (After) = "<<delta_4P.Vect().x();
        cout<<" "<<delta_4P.Vect().y();
        cout<<" "<<delta_4P.Vect().z()<<endl;
        cout<<"Delta Mass (After) = "<<delta_4P.M() * MeV_to_GeV<<endl;
    }

    // --------------------------------------------------------------------
    // Boost All 4-Momentums to Delta Rest Frame
    // --------------------------------------------------------------------
    if (isDebug){
        cout<<"Before Boost"<<endl;
        cout<<"beam 4P = ("<<beam_4P.Px()<<", "<<beam_4P.Py()<<", "<<beam_4P.Pz()<<", "<<beam_4P.E()<<")"<<endl;
        cout<<"muon 4P = ("<<muon_4P.Px()<<", "<<muon_4P.Py()<<", "<<muon_4P.Pz()<<", "<<muon_4P.E()<<")"<<endl;
        cout<<"pi0 4P = ("<<pi0_4P.Px()<<", "<<pi0_4P.Py()<<", "<<pi0_4P.Pz()<<", "<<pi0_4P.E()<<")"<<endl;
    }

    beam_4P.Boost(boost_x,boost_y,boost_z);
    muon_4P.Boost(boost_x,boost_y,boost_z);
    pi0_4P.Boost(boost_x,boost_y,boost_z);

    if (isDebug){
        cout<<"After Boost"<<endl;
        cout<<"beam 4P = ("<<beam_4P.Px()<<", "<<beam_4P.Py()<<", "<<beam_4P.Pz()<<", "<<beam_4P.E()<<")"<<endl;
        cout<<"muon 4P = ("<<muon_4P.Px()<<", "<<muon_4P.Py()<<", "<<muon_4P.Pz()<<", "<<muon_4P.E()<<")"<<endl;
        cout<<"pi0 4P = ("<<pi0_4P.Px()<<", "<<pi0_4P.Py()<<", "<<pi0_4P.Pz()<<", "<<pi0_4P.E()<<")"<<endl;
    }

    // --------------------------------------------------------------------
    // Form Axes
    // --------------------------------------------------------------------
    TVector3 beam_3P = beam_4P.Vect();
    TVector3 muon_3P = muon_4P.Vect();

    // Z Axis -- Momentum Transfer Axis
    TVector3 axis_z(beam_3P.Px() - muon_3P.Px(), beam_3P.Py() - muon_3P.Py(), beam_3P.Pz() - muon_3P.Pz());
    TVector3 unit_axis_z = axis_z.Unit();

    // Y Axis -- Beam x Muon
    TVector3 axis_y = beam_3P.Cross(muon_3P);
    TVector3 unit_axis_y = axis_y.Unit();

    // X Axis -- Right Handed Coordinate System y_axis x z_axis = x_axis
    TVector3 axis_x = axis_y.Cross(axis_z);
    TVector3 unit_axis_x = axis_x.Unit();

    // --------------------------------------------------------------------
    // Calculate Angles
    // --------------------------------------------------------------------
    TVector3 pi0_3P = pi0_4P.Vect();
    TVector3 unit_pi0_3P = pi0_3P.Unit();
    double pion_P = pi0_3P.Mag();

    double cos_theta_z = unit_pi0_3P.Dot(unit_axis_z);

    double theta_x = acos(unit_pi0_3P.Dot(unit_axis_x)); 
    double theta_y = acos(unit_pi0_3P.Dot(unit_axis_y)); 
    double theta_z = acos(cos_theta_z); 

    double Px = pion_P*cos(theta_x);
    double Py = pion_P*cos(theta_y);

    double phi_x = acos(Px/(pion_P*sin(theta_z)));
    double phi = phi_x * TMath::RadToDeg();    

    // Convert to 360
    if (Py < 0 ){
        phi = 360 - phi;
    }

    if (isDebug){
        cout<<"pion_theta = "<<cos_theta_z<<endl;
        cout<<"pion_phi = "<<phi<<endl;
        cout<<endl;
    }

    std::vector<double> angles;
    angles.push_back(cos_theta_z);
    angles.push_back(phi);

    return angles;
}

void CCProtonPi0_Analyzer::GetDelta_pi_angles(bool isTruth)
{
    bool isDebug = false;

    using namespace std;

    // --------------------------------------------------------------------
    // Find Boost from LAB Frame to Delta Rest Frame
    // --------------------------------------------------------------------
    TLorentzVector delta_4P = GetDelta4P(isTruth);
    double M_delta = delta_4P.M(); 
    if (isDebug) cout<<"Delta Mass = "<<M_delta * MeV_to_GeV<<endl;

    double gamma = delta_4P.Gamma();
    double boost_x = -(delta_4P.Px()/ (gamma*M_delta));
    double boost_y = -(delta_4P.Py() / (gamma*M_delta));
    double boost_z = -(delta_4P.Pz() / (gamma*M_delta));

    if (isDebug){
        if (isTruth) cout<<"-------Truth-----------"<<endl;
        cout<<"Delta P (Before) = "<<delta_4P.Vect().Mag()<<endl;
        cout<<"Delta P (Before) = "<<delta_4P.Vect().x();
        cout<<" "<<delta_4P.Vect().y();
        cout<<" "<<delta_4P.Vect().z()<<endl;
    }

    delta_4P.Boost(boost_x,boost_y,boost_z);

    if (isDebug){
        cout<<"Delta P (After) = "<<delta_4P.Vect().Mag()<<std::endl;
        cout<<"Delta P (After) = "<<delta_4P.Vect().x();
        cout<<" "<<delta_4P.Vect().y();
        cout<<" "<<delta_4P.Vect().z()<<endl;
        cout<<"Delta Mass (After) = "<<delta_4P.M() * MeV_to_GeV<<endl;
    }

    // --------------------------------------------------------------------
    // Boost All 4-Momentums to Delta Rest Frame
    // --------------------------------------------------------------------
    TLorentzVector beam_4P = isTruth ? Get_Neutrino_4P(mc_incomingE) : Get_Neutrino_4P(m_Enu);

    double mu_px = isTruth ? truth_muon_4P[0] : muon_px;
    double mu_py = isTruth ? truth_muon_4P[1] : muon_py;
    double mu_pz = isTruth ? truth_muon_4P[2] : muon_pz;
    double mu_E = isTruth ? truth_muon_4P[3] : muon_E;
    TLorentzVector muon_4P(mu_px, mu_py, mu_pz, mu_E); 

    double pi_px = isTruth ? truth_pi0_4P[0] : pi0_px;
    double pi_py = isTruth ? truth_pi0_4P[1] : pi0_py;
    double pi_pz = isTruth ? truth_pi0_4P[2] : pi0_pz;
    double pi_E = isTruth ? truth_pi0_4P[3] : pi0_E;
    TLorentzVector pi0_4P(pi_px, pi_py, pi_pz, pi_E); 

    if (isDebug){
        cout<<"Before Boost"<<endl;
        cout<<"beam 4P = ("<<beam_4P.Px()<<", "<<beam_4P.Py()<<", "<<beam_4P.Pz()<<", "<<beam_4P.E()<<")"<<endl;
        cout<<"muon 4P = ("<<muon_4P.Px()<<", "<<muon_4P.Py()<<", "<<muon_4P.Pz()<<", "<<muon_4P.E()<<")"<<endl;
        cout<<"pi0 4P = ("<<pi0_4P.Px()<<", "<<pi0_4P.Py()<<", "<<pi0_4P.Pz()<<", "<<pi0_4P.E()<<")"<<endl;
    }

    beam_4P.Boost(boost_x,boost_y,boost_z);
    muon_4P.Boost(boost_x,boost_y,boost_z);
    pi0_4P.Boost(boost_x,boost_y,boost_z);

    if (isDebug){
        cout<<"After Boost"<<endl;
        cout<<"beam 4P = ("<<beam_4P.Px()<<", "<<beam_4P.Py()<<", "<<beam_4P.Pz()<<", "<<beam_4P.E()<<")"<<endl;
        cout<<"muon 4P = ("<<muon_4P.Px()<<", "<<muon_4P.Py()<<", "<<muon_4P.Pz()<<", "<<muon_4P.E()<<")"<<endl;
        cout<<"pi0 4P = ("<<pi0_4P.Px()<<", "<<pi0_4P.Py()<<", "<<pi0_4P.Pz()<<", "<<pi0_4P.E()<<")"<<endl;
    }

    // --------------------------------------------------------------------
    // Form Axes
    // --------------------------------------------------------------------
    TVector3 beam_3P = beam_4P.Vect();
    TVector3 muon_3P = muon_4P.Vect();

    // Z Axis -- Momentum Transfer Axis
    TVector3 axis_z(beam_3P.Px() - muon_3P.Px(), beam_3P.Py() - muon_3P.Py(), beam_3P.Pz() - muon_3P.Pz());
    TVector3 unit_axis_z = axis_z.Unit();

    // Y Axis -- Beam x Muon
    TVector3 axis_y = beam_3P.Cross(muon_3P);
    TVector3 unit_axis_y = axis_y.Unit();

    // X Axis -- Right Handed Coordinate System y_axis x z_axis = x_axis
    TVector3 axis_x = axis_y.Cross(axis_z);
    TVector3 unit_axis_x = axis_x.Unit();

    // --------------------------------------------------------------------
    // Calculate Angles
    // --------------------------------------------------------------------
    TVector3 pi0_3P = pi0_4P.Vect();
    TVector3 unit_pi0_3P = pi0_3P.Unit();
    double pion_P = pi0_3P.Mag();

    double cos_theta_z = unit_pi0_3P.Dot(unit_axis_z);

    double theta_x = acos(unit_pi0_3P.Dot(unit_axis_x)); 
    double theta_y = acos(unit_pi0_3P.Dot(unit_axis_y)); 
    double theta_z = acos(cos_theta_z); 

    double Px = pion_P*cos(theta_x);
    double Py = pion_P*cos(theta_y);

    double phi_x = acos(Px/(pion_P*sin(theta_z)));
    double phi = phi_x * TMath::RadToDeg();    

    // Convert to 360
    if (Py < 0 ){
        phi = 360 - phi;
    }

    if (isTruth){
        Delta_pi_theta_true = cos_theta_z;
        Delta_pi_phi_true = phi;
    }else{
        Delta_pi_theta_reco = cos_theta_z;
        Delta_pi_phi_reco = phi;
    }

    if (isDebug){
        cout<<"pion_theta = "<<cos_theta_z<<endl;
        cout<<"pion_phi = "<<phi<<endl;
        cout<<endl;
    }
}

void CCProtonPi0_Analyzer::GetDeltaTransverse()
{
    if (nProtonCandidates > 0){
        if (m_W > 1000 && m_W < 1400 ){
            double delta_transverse_reco = 0;
            if (m_isMC && truth_isSignal && mc_intType == 2 && mc_resID == 0){
                TVector3 unit_beam_muon_reco = GetNeutrinoCrossMuon(false);

                TLorentzVector delta_4P_reco = GetDelta4P(false);
                TVector3 delta_3P_reco = delta_4P_reco.Vect();

                delta_transverse_reco = unit_beam_muon_reco.Dot(delta_3P_reco);


                TVector3 unit_beam_muon_true = GetNeutrinoCrossMuon(true);

                TLorentzVector delta_4P_true = GetDelta4P(true);
                TVector3 delta_3P_true = delta_4P_true.Vect();

                double delta_transverse_true = unit_beam_muon_true.Dot(delta_3P_true);

                FillHistogram(interaction.DeltaTransverse_mc, delta_transverse_true);
                FillHistogram(interaction.DeltaTransverse_mc_res, delta_transverse_true, delta_transverse_true-delta_transverse_reco);
            }else{
                FillHistogram(interaction.DeltaTransverse_data, delta_transverse_reco);
            }
        }
    }
}

TVector3 CCProtonPi0_Analyzer::GetNeutrinoCrossMuon(bool isTruth)
{
    if (isTruth){
        TLorentzVector beam_4P = Get_Neutrino_4P(m_Enu_Truth);
        TLorentzVector muon_4P(truth_muon_4P);

        TVector3 beam_P = beam_4P.Vect();
        TVector3 muon_P = muon_4P.Vect();

        TVector3 beam_muon = beam_P.Cross(muon_P);
        TVector3 unit_beam_muon = beam_muon.Unit();

        return unit_beam_muon;
    }else{
        TLorentzVector beam_4P = Get_Neutrino_4P(m_Enu);
        TVector3 beam_P = beam_4P.Vect();
        TVector3 muon_P(muon_px, muon_py, muon_pz);
        TVector3 beam_muon = beam_P.Cross(muon_P);
        TVector3 unit_beam_muon = beam_muon.Unit();

        return unit_beam_muon;

    }
}


TLorentzVector CCProtonPi0_Analyzer::GetDelta4P(bool isTruth)
{
    double delta_4P[4];

    if (isTruth){
        for(int i = 0; i < 4; i++){
            delta_4P[i] = truth_proton_4P[i] + truth_pi0_4P[i];
        }
    }else{
        delta_4P[0] = proton_px + pi0_px;
        delta_4P[1] = proton_py + pi0_py;
        delta_4P[2] = proton_pz + pi0_pz;
        delta_4P[3] = proton_E + pi0_E;
    }

    TLorentzVector delta_lorentz(delta_4P); 

    return delta_lorentz;
}

bool CCProtonPi0_Analyzer::IsInvMassInRange(double invMass)
{
    return ((invMass >= min_Pi0_invMass) && (invMass <= max_Pi0_invMass)); 
}

bool CCProtonPi0_Analyzer::IsOpeningAngleSmallAndEnergyLow(double E_g1, double E_g2)
{
    bool isAngleSmall = pi0_cos_openingAngle > 0.95;
    bool isEnergyLow = (E_g1 + E_g2) < 400.0;
    if (isEnergyLow && isAngleSmall) return true;
    else return false;
}

int CCProtonPi0_Analyzer::CountEventRecordParticles(int pdg)
{
    int count = 0;
    for (int i = 0; i < mc_er_nPart; ++i){
        // Count only Primary Particles
        if (mc_er_status[i] != 14 ) continue;

        if (mc_er_ID[i] == pdg) count++;
    }
    return count;
}

int CCProtonPi0_Analyzer::CountFSParticles(int pdg, double P_limit)
{
    int count = 0;
    for (int i = 0; i < mc_nFSPart; ++i){
        if (mc_FSPartPDG[i] == pdg ){ 
            double P = HEP_Functions::calcMomentum(mc_FSPartPx[i], mc_FSPartPy[i], mc_FSPartPz[i]);
            if (P > P_limit){
                count++;
            }
        }
    }
    return count;
}

void CCProtonPi0_Analyzer::PrintFSParticles()
{
    for (int i = 0; i < mc_er_nPart; ++i){
        if (mc_er_status[i] == 14){ 
            std::cout<<"\t"<<i<<" "<<mc_er_ID[i]<<std::endl;
        }
    }
}

void CCProtonPi0_Analyzer::GetMichelStatistics()
{
    bool isGamma1_Michel = gamma1_isMichel_begin || gamma1_isMichel_end;
    bool isGamma2_Michel = gamma2_isMichel_begin || gamma2_isMichel_end;

    bool isShowerMichelFound = isGamma1_Michel || isGamma2_Michel; 
    //bool isShowerMichelFound_Improved = ImprovedMichel_Gamma1HasMichel || ImprovedMichel_Gamma2HasMichel;

    bool isMichelFound = (Cut_Vertex_Michel_Exist == 1) || (Cut_EndPoint_Michel_Exist == 1) || (Cut_secEndPoint_Michel_Exist == 1) || isShowerMichelFound;
    //bool isMichelFound_Improved = ImprovedMichel_VertexHasMichel || ImprovedMichel_EndPointHasMichel || ImprovedMichel_secEndPointHasMichel || isShowerMichelFound_Improved; 
    bool isMichelFound_Improved = ImprovedMichel_EventHasMatchedMichel;

    // Total Number of Truth Michels
    if (truth_isBckg_withMichel) nMichel_Truth.increment(cvweight);

    // Total Number of Found Michels
    if (isMichelFound) nMichel_Total_Found.increment(cvweight);
    if (isMichelFound_Improved) nMichel_Total_Found_Improved.increment(cvweight);

    // Truth Found Michels
    if (truth_isBckg_withMichel && isMichelFound) nMichel_Truth_Found.increment(cvweight);
    if (truth_isBckg_withMichel && isMichelFound_Improved) nMichel_Truth_Found_Improved.increment(cvweight);

}

void CCProtonPi0_Analyzer::Study_ProtonSystematics()
{
    FillHistogram(proton.trackLength, proton_length * mm_to_cm);

    // Fill nProtons
    FillHistogram(interaction.nProtons, nProtonCandidates);

    // Fill Energy Shift
    if (nProtonCandidates > 0){
        proton.energy_shift_BetheBloch->Fill(proton_energy_shift_BetheBloch_Down);
        proton.energy_shift_BetheBloch->Fill(proton_energy_shift_BetheBloch_Up);
        proton.energy_shift_Birks->Fill(proton_energy_shift_Birks);
        proton.energy_shift_MEU->Fill(proton_energy_shift_MEU_Down);
        proton.energy_shift_MEU->Fill(proton_energy_shift_MEU_Up);
        proton.energy_shift_Mass->Fill(proton_energy_shift_Mass_Down);
        proton.energy_shift_Mass->Fill(proton_energy_shift_Mass_Up);
        proton.energy_shift_Nominal->Fill(proton_energy_shift_Nominal);
    }

}

void CCProtonPi0_Analyzer::setCounterNames()
{
    // Single Use Counters 
    counter1.setName("N(All)");
    counter2.setName("N(CEX in Detector)");
    counter3.setName("N/A");
    counter4.setName("N/A");

    n2p2h.setName("N(2p2h)");
    for (int i = 0; i < 3; ++i){
        nAll[i].setName("All");
        nAll_Signal[i].setName("TotalSignal");
        nAll_Bckg[i].setName("TotalBckg");
        nSignal_Delta_RES[i].setName("Signal:DeltaRES");
        nSignal_Other_RES[i].setName("Signal:OtherRES");
        nSignal_Non_RES[i].setName("Signal:NonRES");
        nBckg_WithPi0[i].setName("Bckg:WithPi0");
        nBckg_QELike[i].setName("Bckg:QELike");
        nBckg_PiPlus[i].setName("Bckg:ChargedPi");
        nBckg_Other[i].setName("Bckg:Other");
    }
    // Signal Counters
    nSignalOut_Acceptance.setName("SignalOut_Acceptance");
    nSignalOut_Kinematics.setName("SignalOut_Kinematics");

    // Michel Counters
    nMichel_Truth.setName("nMichel_Truth");
    nMichel_Total_Found.setName("nMichel_Total_Found");
    nMichel_Total_Found_Improved.setName("nMichel_Total_Found_Improved");
    nMichel_Truth_Found.setName("nMichel_Truth_Found");
    nMichel_Truth_Found_Improved.setName("nMichel_Truth_Found_Improved");
}

void CCProtonPi0_Analyzer::initCVWeights()
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

void CCProtonPi0_Analyzer::Study_GENIE_Weights()
{
    std::cout<<"MaRES, MvRES and Rvn1pi Shifts"<<std::endl;
    for (int i = 0; i < genie_wgt_n_shifts; ++i){
        std::cout<<i<<"\t"<<truth_genie_wgt_MaRES[i]<<" "<<truth_genie_wgt_MvRES[i]<<" "<<truth_genie_wgt_Rvn1pi[i]<<std::endl; 
    }
}

double CCProtonPi0_Analyzer::GetMaResWeight( double newMaRes )
{
    // Indices 2 and 4 correspond to -1 and +1 sigma shifts, respectively
    int weightIndex = newMaRes < genieMaRes ? 2 : 4;
    return 1.0 + fabs( newMaRes - genieMaRes ) * ( truth_genie_wgt_MaRES[weightIndex] - 1.0 ) / genieMaRes1sig;
}

double CCProtonPi0_Analyzer::GetMvResWeight(double newMvRes )
{
  // Indices 2 and 4 correspond to -1 and +1 sigma shifts, respectively
  int weightIndex = newMvRes < genieMvRes ? 2 : 4;
  return 1.0 + fabs( newMvRes - genieMvRes ) * ( truth_genie_wgt_MvRES[weightIndex] - 1.0 ) / genieMvRes1sig;
}

bool CCProtonPi0_Analyzer::IsGenieCCRes()
{
    return mc_current == 1 && mc_intType == 2;
}

bool CCProtonPi0_Analyzer::IsGenieRvn1pi()
{
    // Weight is NOT equal to 1 for Rvn1pi events
    return truth_genie_wgt_Rvn1pi[2] < 1;
}

bool CCProtonPi0_Analyzer::IsGenieRvp1pi()
{
    // Weight is NOT equal to 1 for Rvp1pi events
    return truth_genie_wgt_Rvp1pi[2] < 1;
}

bool CCProtonPi0_Analyzer::IsGenieNonRES()
{
    // NonRES 1pi 
    if ( IsGenieRvn1pi() || IsGenieRvp1pi() ) return true;

    // NonRES 2pi
    if ( truth_genie_wgt_Rvn2pi[2] < 1 || truth_genie_wgt_Rvp2pi[2] < 1) return true;

    return false;
}

bool CCProtonPi0_Analyzer::IsGenieDeltaRES()
{
    if (truth_genie_wgt_Theta_Delta2Npi[4] == 1.0) return false;
    else return true;
}

bool CCProtonPi0_Analyzer::IsGenie_NonRES_n_piplus()
{
    if ( IsGenieRvn1pi() || IsGenieRvp1pi() ){
        std::vector<int> prim_part = GetPrimaryParticles();

        if (prim_part.empty()) return false;

        if (prim_part[0] == 2112 && prim_part[1] == 211) return true;
        else return false;
    }else{
        return false;
    }
}

void CCProtonPi0_Analyzer::Test_GENIE_DIS()
{
    if ( mc_intType != 3 ) return;

    if ( mc_w > 2000 ) return;

    if ( IsGenieNonRES() ) return; 

    std::cout<<"isSignal = "<<truth_isSignal<<" mc_intType = "<<mc_intType<<" mc_w = "<<mc_w<<std::endl;
    std::cout<<"\ttruth_genie_wgt_Rvn1pi[2] = "<<truth_genie_wgt_Rvn1pi[2]<<std::endl;
    std::cout<<"\ttruth_genie_wgt_Rvp1pi[2] = "<<truth_genie_wgt_Rvp1pi[2]<<std::endl;
    std::cout<<"\ttruth_genie_wgt_Rvn2pi[2] = "<<truth_genie_wgt_Rvn2pi[2]<<std::endl;
    std::cout<<"\ttruth_genie_wgt_Rvp2pi[2] = "<<truth_genie_wgt_Rvp2pi[2]<<std::endl;
}

std::vector<int> CCProtonPi0_Analyzer::GetPrimaryParticles()
{
    std::vector<int> prim_part;

    for (int i = 0; i < mc_er_nPart; ++i){
        if (mc_er_status[i] == 14){ 
            prim_part.push_back(mc_er_ID[i]);
        }
    }

    return prim_part;
}

void CCProtonPi0_Analyzer::Study_BckgSubtraction()
{
    // Fill Pi0 Inv Mass used in Background Subtraction during Studies
    FillHistogram(interaction.pi0_invMass_All, pi0_invMass);
    if ( nProtonCandidates == 0){
        FillHistogram(interaction.pi0_invMass_1Track, pi0_invMass);
    }else{
        FillHistogram(interaction.pi0_invMass_2Track, pi0_invMass);
    }
}

void CCProtonPi0_Analyzer::Study_W()
{
    double W_reco = m_W * MeV_to_GeV;
    double deltaInvMass = calcDeltaInvariantMass(false) * MeV_to_GeV;

    FillHistogram(interaction.W_p_pi0, deltaInvMass);
    FillHistogram(interaction.W_All, W_reco);

    if (nProtonCandidates == 0){
        FillHistogram(interaction.W_1, W_reco);
    }else{
        FillHistogram(interaction.W_2, W_reco);
    }
    
    // Shifted W
    for (int i = 0; i < 151; ++i){
        double W_shifted = m_isMC ? (m_W - i) * MeV_to_GeV : W_reco;
        FillHistogram(interaction.W_Shift[i], W_shifted);
        if (m_isMC){
            if (truth_isSignal) {
                FillHistogram(interaction.W_Shift_Signal[i], W_shifted);
            }else{
                FillHistogram(interaction.W_Shift_Bckg[i], W_shifted);
            }
        }
    }
}

void CCProtonPi0_Analyzer::Study_NonRES()
{
    if ( IsGenieRvn1pi() || IsGenieRvp1pi() ){
        cout<<"GENIE Tag!";
        cout<<"\tmc_intType = "<<mc_intType<<" mc_w = "<<mc_w<<endl;
    }else if (mc_intType == 3){
        cout<<"mc_intType = "<<mc_intType<<" mc_w = "<<mc_w<<endl;
    }
}

void CCProtonPi0_Analyzer::Study_DeltaResonance()
{
    //if (IsGenieCCRes()){
    //    double wgt = truth_genie_wgt_Theta_Delta2Npi[4];
    //    if (wgt == 1.0){
    //        FillHistogram(interaction.resID, mc_resID);
    //    }else{
    //        FillHistogram(interaction.resID_theta, mc_resID);
    //        std::vector<int> prim_part = GetPrimaryParticles();
    //        for (unsigned int i = 0; i < prim_part.size(); ++i){
    //            std::cout<<prim_part[i]<<std::endl;
    //        }
    //        std::cout<<" "<<std::endl;
    //    }
    //}
}

void CCProtonPi0_Analyzer::Study_QSq()
{
    double QSq_reco = m_QSq * MeVSq_to_GeVSq;

    // Central Value Fill
    FillHistogram(interaction.QSq_CV, QSq_reco);

    // ------------------------------------------------------------------------
    //  MaRES Fit -- MaRES Change handled in Vertical Error Bands
    // ------------------------------------------------------------------------
    //FillHistogramWithVertErrors(interaction.QSq_MaRES, QSq_reco);

        
    // Vary Q0
    std::vector<double> Q0_Vector;
    Q0_Vector.push_back(0.156); // MINOS
    for (double q0 = 0.050; q0 <= 0.155; q0 += 0.001){
        Q0_Vector.push_back(q0);
    }

    for (unsigned int i = 0; i <= Q0_Vector.size(); ++i){
        std::cout<<i<<" "<<Q0_Vector[i]<<std::endl;
    }


}

double CCProtonPi0_Analyzer::GetDeltaFactor(double QSq, double A, double Q0)
{
    double factor;

    if (IsGenieCCRes()){
        factor = A / (1 + exp(1-(sqrt(QSq)/Q0)));
    }else{
        factor = 1.0;
    }

    return factor;
}

void CCProtonPi0_Analyzer::updateCounters()
{
    if(truth_isSignalOut_Acceptance) nSignalOut_Acceptance.increment(cvweight);
    if(truth_isSignalOut_Kinematics) nSignalOut_Kinematics.increment(cvweight);

    if (IsEvent2p2h(mc_intType)) n2p2h.increment(cvweight);

    // Counters for ALL Events
    update_EventTypeCounters(0);
  
    // Counters for Enu < 4 GeV and Enu > 4 GeV
    if (m_Enu < 4000){
        update_EventTypeCounters(1);
    }else{ 
        update_EventTypeCounters(2);
    }
}

void CCProtonPi0_Analyzer::update_EventTypeCounters(int i)
{
    nAll[i].increment(cvweight);

    if (truth_isSignal) nAll_Signal[i].increment(cvweight);
    else nAll_Bckg[i].increment(cvweight);

    if (truth_isSignal){
        if (mc_intType == 2){
            if (mc_resID == 0){
                // Delta Resonance
                nSignal_Delta_RES[i].increment(cvweight);
            }else{
                // Other Resonance
                nSignal_Other_RES[i].increment(cvweight);
            }
        }else{
            // Non-Resonance
            nSignal_Non_RES[i].increment(cvweight);
        }
    }else{
        if (truth_isBckg_Compact_WithPi0) nBckg_WithPi0[i].increment(cvweight);
        else if (truth_isBckg_Compact_QELike) nBckg_QELike[i].increment(cvweight);
        else if (truth_isBckg_Compact_SinglePiPlus) nBckg_PiPlus[i].increment(cvweight);
        else if (truth_isBckg_Compact_Other) nBckg_Other[i].increment(cvweight);
        else{
            RunTimeError("No Background Type found in updateCounters()");
        }
    }
}

void CCProtonPi0_Analyzer::printCounters() 
{
    n2p2h.print();
    nSignalOut_Acceptance.print();
    nSignalOut_Kinematics.print();
    counter1.print(); 
    counter2.print(); 
    counter3.print(); 
    counter4.print(); 
}

void CCProtonPi0_Analyzer::writeEventTypeTable()
{
    eventTypeTable<<std::fixed;
    eventTypeTable<<std::setprecision(1);
    eventTypeTable<<std::left;

    for (int i = 0; i < 3; ++i){
        writeEventTypeTableLine(nAll[i], nAll[i]);        
        writeEventTypeTableLine(nAll_Signal[i], nAll[i]);        
        writeEventTypeTableLine(nAll_Bckg[i], nAll[i]);        
        eventTypeTable<<std::endl; 
        
        writeEventTypeTableLine(nSignal_Delta_RES[i], nAll_Signal[i]);        
        writeEventTypeTableLine(nSignal_Other_RES[i], nAll_Signal[i]);        
        writeEventTypeTableLine(nSignal_Non_RES[i], nAll_Signal[i]);        
        writeEventTypeTableLine(nBckg_WithPi0[i], nAll_Bckg[i]);        
        writeEventTypeTableLine(nBckg_PiPlus[i], nAll_Bckg[i]);        
        writeEventTypeTableLine(nBckg_QELike[i], nAll_Bckg[i]);        
        writeEventTypeTableLine(nBckg_Other[i], nAll_Bckg[i]);        
       
        eventTypeTable<<std::endl; 
        eventTypeTable<<std::endl; 
    }
}

void CCProtonPi0_Analyzer::writeEventTypeTableLine(CCProtonPi0_Counter &counter, CCProtonPi0_Counter &base)
{
    eventTypeTable.width(20); eventTypeTable<<counter.getName(); 
    eventTypeTable.width(8); eventTypeTable<<counter.getCount();
    eventTypeTable.width(8); eventTypeTable<<counter.calcPercent(base.getCount());
    eventTypeTable<<std::endl; 
}

double CCProtonPi0_Analyzer::Calc_q3(bool isTruth)
{
    double t_Enu = isTruth ? m_Enu_Truth : m_Enu;
    double t_muon_Px = isTruth ? truth_muon_4P[0] : muon_px;
    double t_muon_Py = isTruth ? truth_muon_4P[1] : muon_py;
    double t_muon_Pz = isTruth ? truth_muon_4P[2] : muon_pz;
    
    double nu_Px = 0.0;
    double nu_Py = t_Enu * sin(beam_theta);
    double nu_Pz = t_Enu * cos(beam_theta);

    double dPx = nu_Px - t_muon_Px;
    double dPy = nu_Py - t_muon_Py;
    double dPz = nu_Pz - t_muon_Pz;

    double dP = HEP_Functions::calcMomentum(dPx, dPy, dPz);

    return dP;
}

void CCProtonPi0_Analyzer::Study_Flux()
{
    double Enu = mc_incomingE * MeV_to_GeV;
    FillHistogram(interaction.Enu_flux_wgt, Enu, cvweight_Flux);
    FillHistogram(interaction.Enu_cvweight, Enu, cvweight);
}

void CCProtonPi0_Analyzer::Study_Paper()
{

}

bool CCProtonPi0_Analyzer::IsProtonPi0(std::vector<int> hadrons)
{
    int nPi0 = 0;
    int nProton = 0;
    int nOther = 0;
    for (unsigned int i = 0; i < hadrons.size(); ++i){
        //std::cout<<hadrons[i]<<std::endl;
        if (hadrons[i] == 111) nPi0++;
        else if (hadrons[i] == 2212) nProton++;
        else nOther++;
    }

    if (nPi0 == 1 && nProton == 1) return true;
    else return false;
}

void CCProtonPi0_Analyzer::ReviseSignal()
{
    if (isSignalOriginal) return;

    if (isSignalOriginal_NoWLimit){
        truth_isSignal = truth_isSignal || (truth_isSignalOut_Kinematics && IsEnuInRange(m_Enu_Truth));
     }

    if (isSignalTwoTrack){
        double Tp = truth_proton_4P[3] - proton_mass;
        truth_isSignal = IsProtonLong(Tp);
    }
    
    if (isSignalDeltaRich){
        truth_isSignal = isDeltaRichSample_Truth();
    }

    if (isSignalLowEnu){
        truth_isSignal = IsEnuInRange_LowEnu(m_Enu_Truth);
    }

    if (isSignalHighEnu){
        truth_isSignal = IsEnuInRange_HighEnu(m_Enu_Truth);
    }
}

void CCProtonPi0_Analyzer::ReviseBackground()
{

    bool isDebug = false;

    // Make All False
    bool WithPi0 = false;
    bool QELike = false;
    bool ChargedMeson = false;
    bool Other = false;
    
    // -------------------------------------------------------------------------
    // Count Final State Particles 
    // -------------------------------------------------------------------------
    int nAntiMuon = 0;
    int nPiZero = 0;
    int nNucleon = 0;
    int nChargedMeson = 0;
    int nOther = 0;

    if(isDebug) std::cout<<ev_run<<" "<<ev_subrun<<" "<<ev_gate<<std::endl;

    for (int i = 0; i < mc_er_nPart; ++i){

        // Check only Stable Final State Particles
        if (mc_er_status[i] != 1) continue;

        int PDG = mc_er_ID[i];

        if (isDebug){
            std::cout<<PDG<<" ";
        }

        if (std::abs(PDG) == 13 || PDG > 10000) continue;

        if (PDG == 2212 || PDG == 2112) nNucleon++;
        else if (PDG == 111) nPiZero++;
        else if (std::abs(PDG) == 211 || std::abs(PDG) == 311) nChargedMeson++;
        else nOther++;
    }

    if (isDebug) std::cout<<std::endl;

    // -------------------------------------------------------------------------
    // Background Types 
    // -------------------------------------------------------------------------
    if (truth_isSignalOut_Acceptance) Other = true;
    else if (nAntiMuon > 0) Other = true;
    else if (nPiZero > 0 || truth_isSignalOut_Kinematics) WithPi0 = true;
    else if (nNucleon > 0 && nChargedMeson == 0 && nOther == 0) QELike = true;
    else if (nChargedMeson > 0) ChargedMeson = true;
    else Other = true;

    if (isDebug){
        if (truth_isBckg_Compact_WithPi0 && !WithPi0){
            std::cout<<"Printing Event Record"<<std::endl;
            for (int i = 0; i < mc_er_nPart; ++i){
                std::cout<<"\t"<<mc_er_status[i]<<"\t"<<mc_er_ID[i]<<std::endl;
            }
            std::cout<<"Printing Detector Trajectory"<<std::endl;
            std::cout<<"ID\tStatus\tPDG\tMother"<<std::endl;
            for (int i = 0; i < detmc_traj_id_sz; ++i){
                std::cout<<detmc_traj_id[i]<<"\t";
                std::cout<<detmc_traj_status[i]<<"\t";
                std::cout<<detmc_traj_pdg[i]<<"\t";
                std::cout<<detmc_traj_mother[i]<<std::endl;
            }
        }
        std::cout<<"\tWithPi0 = "<<truth_isBckg_Compact_WithPi0<<" "<<WithPi0<<std::endl;
        std::cout<<"\tQELike = "<<truth_isBckg_Compact_QELike<<" "<<QELike<<std::endl;
        std::cout<<"\tChargedMeson = "<<truth_isBckg_Compact_SinglePiPlus<<" "<<ChargedMeson<<std::endl;
        std::cout<<"\tOther = "<<truth_isBckg_Compact_Other<<" "<<Other<<std::endl;
    }

    truth_isBckg_Compact_WithPi0 = WithPi0;
    truth_isBckg_Compact_QELike = QELike;
    truth_isBckg_Compact_SinglePiPlus = ChargedMeson;
    truth_isBckg_Compact_Other = Other;

}

bool CCProtonPi0_Analyzer::isDeltaRichSample()
{
    return (nProtonCandidates > 0) && (m_W > 0 && m_W < 1400);
}

bool CCProtonPi0_Analyzer::isDeltaRichSample_Truth()
{
    double Tp = truth_proton_4P[3] - proton_mass;
    return (IsWInRange_DeltaRich(m_W_Truth) && IsProtonLong(Tp));
}


#endif //CCProtonPi0_Analyzer_cpp


