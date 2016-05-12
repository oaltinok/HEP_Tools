/*
 * See CCProtonPi0_Analyzer.h header or Class Information
 */
#ifndef CCProtonPi0_Analyzer_cpp
#define CCProtonPi0_Analyzer_cpp

#include "CCProtonPi0_Analyzer.h"

using namespace std;

void CCProtonPi0_Analyzer::specifyRunTime()
{
    applyMaxEvents = false;
    nMaxEvents = 10000;

    // Control Flow
    isDataAnalysis  = true;
    isScanRun = false;
    applyBckgConstraints = true;
    writeFSParticleMomentum = false;

    // Side Band Control
    NoSideBand = true;
    sideBand_Michel = false;
    sideBand_PID = false;
    sideBand_LowInvMass = false;
    sideBand_HighInvMass = false;

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

    // Counter Names
    counter1.name = "N(+WSq 1 Track) = ";
    counter2.name = "N(+WSq 2 Track) = ";
    counter3.name = "N(-WSq 1 Track) = ";
    counter4.name = "N(-WSq 2 Track) = ";
}

void CCProtonPi0_Analyzer::reduce(string playlist)
{
    string rootDir;
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

        // Update scanFileName if running for scan
        if(isScanRun) UpdateScanFileName();

        // Progress Message on Terminal
        if (jentry%50000 == 0) cout<<"\tEntry "<<jentry<<endl;


        if (applyMaxEvents && jentry == nMaxEvents){
            cout<<"\tReached Event Limit!"<<endl;
            break;
        }

        UpdateSignalDef();
        
        Calc_WeightFromSystematics();

        CorrectEMShowerCalibration();
        Calc_EventKinematics();

        // Get Cut Statistics
        isPassedAllCuts = getCutStatistics();
        if( !isPassedAllCuts ) continue;

        tree->Fill();
    }
    
    if (!m_isMC) AddVertErrorBands_Data(cutList.invMass_all);
    
    cutList.writeCutTable();
    cutList.writeHistograms();

    cout<<">> Writing "<<rootDir<<endl;
    tree->AutoSave();    
    f->Write();

    //--------------------------------------------------------------------------
    // Counters
    //--------------------------------------------------------------------------
    cout<<counter1.name<<counter1.count<<endl;
    cout<<counter2.name<<counter2.count<<endl;
    cout<<counter3.name<<counter3.count<<endl;
    cout<<counter4.name<<counter4.count<<endl;
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
        if (jentry%25000 == 0) cout<<"\tEntry "<<jentry<<endl;

        if (applyMaxEvents && jentry == nMaxEvents){
            cout<<"\tReached Event Limit!"<<endl;
            break;
        }
       
        UpdateSignalDef();

        // Update scanFileName if running for scan
        if(isScanRun) UpdateScanFileName();

        //----------------------------------------------------------------------
        // Fill Background Branches for Background Events
        //----------------------------------------------------------------------
        if(m_isMC && !truth_isSignal) {
            if (nProtonCandidates == 0) bckgTool.set_nTracks(1);
            else bckgTool.set_nTracks(2);
            bckgTool.fillBackgroundCompact(truth_isBckg_Compact_WithPi0, truth_isBckg_Compact_QELike, truth_isBckg_Compact_SinglePiPlus, truth_isBckg_Compact_Other);                                    
            bckgTool.fillBackgroundWithPi0(truth_isBckg_NoPi0, truth_isBckg_SinglePi0, truth_isBckg_MultiPi0, truth_isBckg_withMichel);                                    
            bckgTool.fillBackground(truth_isBckg_NC, truth_isBckg_AntiNeutrino, truth_isBckg_QELike, truth_isBckg_SingleChargedPion,truth_isBckg_SingleChargedPion_ChargeExchanged, truth_isBckg_DoublePionWithPi0, truth_isBckg_DoublePionWithoutPi0, truth_isBckg_MultiPionWithPi0, truth_isBckg_MultiPionWithoutPi0, truth_isBckg_Other, truth_isBckg_withMichel);                                    
        }

        Calc_WeightFromSystematics();
        Calc_EventKinematics();
        
        //----------------------------------------------------------------------
        // Data Analysis and Other Studies
        //----------------------------------------------------------------------
        if (isDataAnalysis) fillData();
        if (writeFSParticleMomentum) writeFSParticle4P(jentry);

    } // end for-loop

    if (!m_isMC) AddErrorBands_Data();
    //--------------------------------------------------------------------------
    // Studies
    //--------------------------------------------------------------------------

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
    cout<<counter1.name<<counter1.count<<endl;
    cout<<counter2.name<<counter2.count<<endl;
    cout<<counter3.name<<counter3.count<<endl;
    cout<<counter4.name<<counter4.count<<endl;
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
    cutList(isModeReduce, isMC)
{   
    cout<<"Initializing CCProtonPi0_Analyzer"<<endl;

    m_isMC = isMC;

    cvweight = 1.0;
   
    ResetCounters();

    specifyRunTime();

    openTextFiles();

    cout<<"Initialization Finished!\n"<<endl;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
//
//     Specific Functions
//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

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

void CCProtonPi0_Analyzer::fillData()
{          
    // Fill Cross Section Variables
    fill_pi0_P();
    fill_pi0_KE();
    fill_pi0_theta();
    fill_muon_P();
    fill_muon_theta();
    fill_muon_cos_theta();
    fill_QSq();
    fill_Enu();

    // Fill Reconstructed Information
    fillInteractionReco();
    fillMuonReco();
    fillPi0Reco();
    fillPi0BlobReco();
    if(nProtonCandidates > 0) fillProtonReco();

    // Fill Truth Information if Exist and Set Errors
    if( m_isMC ){
        fillInteractionMC();
        fillMuonMC();
        if(nProtonCandidates > 0) fillProtonMC();
        fillPi0MC();
        fillPi0BlobMC();
    }
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
    cutList.nCut_All.increment(truth_isSignal, study1, study2);

    // ------------------------------------------------------------------------
    // Vertex Cut -- If Cut_Vertex_None == 1 --> No Event Vertex
    // ------------------------------------------------------------------------
    if( Cut_Vertex_None == 1) return false;
    cutList.nCut_Vertex_None.increment(truth_isSignal, study1, study2);

    // ------------------------------------------------------------------------
    // Vertex Reconstructable Cut
    // ------------------------------------------------------------------------
    if( Cut_Vertex_Not_Reconstructable == 1) return false;
    cutList.nCut_Vertex_Not_Reconstructable.increment(truth_isSignal, study1, study2);

    // ------------------------------------------------------------------------
    // Vertex Fiducial Cut
    // ------------------------------------------------------------------------
    if( Cut_Vertex_Not_Fiducial == 1) return false;
    cutList.nCut_Vertex_Not_Fiducial.increment(truth_isSignal, study1, study2);

    // ------------------------------------------------------------------------
    // Muon Cut -- If Cut_Muon_None == 1 --> No MINOS Matched Muon
    // ------------------------------------------------------------------------
    if( Cut_Muon_None == 1) return false;
    cutList.nCut_Muon_None.increment(truth_isSignal, study1, study2);

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
    cutList.nCut_Muon_Charge.increment(truth_isSignal, study1, study2);

    // ------------------------------------------------------------------------
    // Michel Cuts
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
  
    isMichelEvent = (Cut_Vertex_Michel_Exist == 1) || (Cut_EndPoint_Michel_Exist == 1) || (Cut_secEndPoint_Michel_Exist == 1);
    if( isMichelEvent){
        FillHistogram(cutList.hCut_Michel,1);
    }else{
        FillHistogram(cutList.hCut_Michel,0);
    } 
    if( Cut_Vertex_Michel_Exist == 1 && !sideBand_Michel ) return false;
    cutList.nCut_Vertex_Michel_Exist.increment(truth_isSignal, study1, study2);

    if( Cut_EndPoint_Michel_Exist == 1 && !sideBand_Michel) return false;
    cutList.nCut_EndPoint_Michel_Exist.increment(truth_isSignal, study1, study2);

    if( Cut_secEndPoint_Michel_Exist == 1 && !sideBand_Michel) return false;
    cutList.nCut_secEndPoint_Michel_Exist.increment(truth_isSignal, study1, study2);
    
    // After Michel I try to save far tracks which changes the Number of Vertices
    // Check nVertices
    FillHistogram(cutList.hCut_nVertices, vtx_total_count);

    // ------------------------------------------------------------------------
    // Tracked Particle Reconstruction Fails
    // ------------------------------------------------------------------------
    if( Cut_Particle_None == 1) return false;
    cutList.nCut_Particle_None.increment(truth_isSignal, study1, study2);

    // ------------------------------------------------------------------------
    // Cannot Find Proton in Tracked Particles
    // ------------------------------------------------------------------------
    if( Cut_Proton_None == 1) return false;
    cutList.nCut_Proton_None.increment(truth_isSignal, study1, study2);

    // ------------------------------------------------------------------------
    // Proton Momentum NaN
    // ------------------------------------------------------------------------
    if ( Cut_Proton_Bad == 1) return false;
    cutList.nCut_Proton_Bad.increment(truth_isSignal, study1, study2);

    // Check nProtonCandidates
    FillHistogram(cutList.hCut_nProtonCandidates, nProtonCandidates);

    // After this stage we can analyze different topologies
    //      1Track = No Proton Events (Only Muon Track)
    //      2Track = With Proton (Muon + Proton)
    if (nProtonCandidates == 0){
        cutList.nCut_1Track_All.increment(truth_isSignal, study1, study2);
    }else{
        cutList.nCut_2Track_All.increment(truth_isSignal, study1, study2);
    }

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
        cutList.nCut_2Track_ProtonScore.increment(truth_isSignal, study1, study2);
    }
    cutList.nCut_ProtonScore.increment(truth_isSignal, study1, study2);

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
    cutList.nCut_PreFilter_Pi0.increment(truth_isSignal, study1, study2);
    if (nProtonCandidates == 0) cutList.nCut_1Track_PreFilter_Pi0.increment(truth_isSignal, study1, study2);
    else cutList.nCut_2Track_PreFilter_Pi0.increment(truth_isSignal, study1, study2);

    // ------------------------------------------------------------------------
    // ConeBlobs Cut -- If Cut_ConeBlobs == 1 --> Failed Pi0 Reconstruction
    // ------------------------------------------------------------------------
    FillHistogram(cutList.hCut_nShowerCandidates,anglescan_ncand); 
    if (nProtonCandidates == 0) FillHistogram(cutList.hCut_1Track_nShowerCandidates,anglescan_ncand); 
    else FillHistogram(cutList.hCut_2Track_nShowerCandidates,anglescan_ncand); 

    if( Cut_ConeBlobs == 1 ) return false;
    cutList.nCut_ConeBlobs.increment(truth_isSignal, study1, study2);
    if (nProtonCandidates == 0) cutList.nCut_1Track_ConeBlobs.increment(truth_isSignal, study1, study2);
    else cutList.nCut_2Track_ConeBlobs.increment(truth_isSignal, study1, study2);

    // ------------------------------------------------------------------------
    // Blob Direction Bad Cut
    // ------------------------------------------------------------------------
    if ( Cut_BlobDirectionBad == 1 ) return false;
    cutList.nCut_BlobDirectionBad.increment(truth_isSignal, study1, study2);
    if (nProtonCandidates == 0) cutList.nCut_1Track_BlobDirectionBad.increment(truth_isSignal, study1, study2);
    else cutList.nCut_2Track_BlobDirectionBad.increment(truth_isSignal, study1, study2);

    // ------------------------------------------------------------------------
    // Pi0 Momentum NaN
    // ------------------------------------------------------------------------
    if ( Cut_Pi0_Bad == 1) return false;
    cutList.nCut_Pi0_Bad.increment(truth_isSignal, study1, study2);
    if (nProtonCandidates == 0) cutList.nCut_1Track_Pi0_Bad.increment(truth_isSignal, study1, study2);
    else cutList.nCut_2Track_Pi0_Bad.increment(truth_isSignal, study1, study2);

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
    cutList.nCut_Shower_Michel_Exist.increment(truth_isSignal, study1, study2);
    
    // ------------------------------------------------------------------------
    // Gamma1 Conv Length Cut Hist
    // ------------------------------------------------------------------------
    if (nProtonCandidates == 0) FillHistogram(cutList.hCut_1Track_gamma1ConvDist,gamma1_dist_vtx * 0.1);
    else FillHistogram(cutList.hCut_2Track_gamma1ConvDist,gamma1_dist_vtx * 0.1);

    if (applyPhotonDistance && gamma1_dist_vtx * 0.1 < minPhotonDistance_1) return false;
    cutList.nCut_Photon1DistanceLow.increment(truth_isSignal, study1, study2);
    if (nProtonCandidates == 0) cutList.nCut_1Track_Photon1DistanceLow.increment(truth_isSignal, study1, study2);
    else cutList.nCut_2Track_Photon1DistanceLow.increment(truth_isSignal, study1, study2);

    // ------------------------------------------------------------------------
    // Gamma2 Conv Length Cut
    // ------------------------------------------------------------------------
    if (nProtonCandidates == 0) FillHistogram(cutList.hCut_1Track_gamma2ConvDist,gamma2_dist_vtx * 0.1);
    else FillHistogram(cutList.hCut_2Track_gamma2ConvDist,gamma2_dist_vtx * 0.1);

    if (applyPhotonDistance && gamma2_dist_vtx * 0.1 < minPhotonDistance_2) return false;
    cutList.nCut_Photon2DistanceLow.increment(truth_isSignal, study1, study2);
    if (nProtonCandidates == 0) cutList.nCut_1Track_Photon2DistanceLow.increment(truth_isSignal, study1, study2);
    else cutList.nCut_2Track_Photon2DistanceLow.increment(truth_isSignal, study1, study2);

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
    bool isAngleSmall = pi0_cos_openingAngle > 0.95;
    bool isEnergyLow = (gamma1_E+gamma2_E) < 400.0;
    if (isEnergyLow && isAngleSmall) return false;
    cutList.nCut_LowE_SmallAngle.increment(truth_isSignal, study1, study2);

    // ------------------------------------------------------------------------
    // Neutrino Energy Cut
    // ------------------------------------------------------------------------
    if (nProtonCandidates == 0) FillHistogram(cutList.hCut_1Track_neutrinoE,m_Enu * MeV_to_GeV);
    else FillHistogram(cutList.hCut_2Track_neutrinoE,m_Enu * MeV_to_GeV); 

    bool isEnu_inRange = m_Enu >= min_Enu && m_Enu <= max_Enu;
    if ( !isEnu_inRange ) return false;

    cutList.nCut_beamEnergy.increment(truth_isSignal, study1, study2);
    if (nProtonCandidates == 0) cutList.nCut_1Track_beamEnergy.increment(truth_isSignal, study1, study2);
    else cutList.nCut_2Track_beamEnergy.increment(truth_isSignal, study1, study2);

    // ------------------------------------------------------------------------
    // Pi0 Invariant Mass Cut
    // ------------------------------------------------------------------------
    fill_BackgroundSubtractionHists();

    // Fill Invariant Mass Histograms
    //      If there is no Side Band, Fill for Every Event
    //      Else fill according to Side Band
    if (NoSideBand){
        FillInvMass_TruthMatch();
        FillHistogram(cutList.hCut_pi0invMass, pi0_invMass);

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
    cutList.nCut_Pi0_invMass.increment(truth_isSignal, study1, study2);
    if (nProtonCandidates == 0) cutList.nCut_1Track_Pi0_invMass.increment(truth_isSignal, study1, study2);
    else cutList.nCut_2Track_Pi0_invMass.increment(truth_isSignal, study1, study2);

    // Fill Other Side Bands 
    //      If there is no Side Band, Fill for Every Event
    //      Else fill according to Side Band
    if (NoSideBand){
        FillHistogram(cutList.SideBand_muon_P, muon_P*MeV_to_GeV);
        FillHistogram(cutList.SideBand_muon_theta, muon_theta_beam*TMath::RadToDeg());
        FillHistogram(cutList.SideBand_pi0_P, pi0_P*MeV_to_GeV);
        FillHistogram(cutList.SideBand_pi0_KE, pi0_KE*MeV_to_GeV);
        FillHistogram(cutList.SideBand_pi0_theta, pi0_theta_beam*TMath::RadToDeg());
        FillHistogram(cutList.SideBand_neutrino_E, m_Enu*MeV_to_GeV);
        FillHistogram(cutList.SideBand_QSq, m_QSq*MeVSq_to_GeVSq);
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

        // Neutrino Energy: Truth, Error, Difference
        double E_true = mc_incomingE * MeV_to_GeV;
        double E_reco = m_Enu * MeV_to_GeV;
        double E_Error = Data_Functions::getError(E_true, E_reco);

        FillHistogram(interaction.Enu_Error, E_Error);
        FillHistogram(interaction.Enu_Diff, E_reco-E_true);

        // Fill 1Track and 2 Track Enu 
        if(nProtonCandidates == 0){
            FillHistogramWithDefaultErrors(interaction.Enu_1Track_response, E_reco, E_true);
            FillHistogram(interaction.Enu_1Track_Error, E_Error);
            FillHistogram(interaction.Enu_1Track_Diff, E_reco-E_true);
        }else{ 
            FillHistogramWithDefaultErrors(interaction.Enu_2Track_response, E_reco, E_true);
            FillHistogram(interaction.Enu_2Track_Error, E_Error);
            FillHistogram(interaction.Enu_2Track_Diff, E_reco-E_true);
        }  
    
        // QSq True, Error and Difference
        double QSq_true = mc_Q2 * MeVSq_to_GeVSq;
        double QSq_reco = m_QSq * MeVSq_to_GeVSq;
        double QSq_error = Data_Functions::getError(QSq_true,QSq_reco);

        FillHistogram(interaction.QSq_Error, QSq_error);
        FillHistogram(interaction.QSq_Diff, QSq_reco-QSq_true);
 
        if(nProtonCandidates == 0){
            FillHistogramWithDefaultErrors(interaction.QSq_1Track_response, QSq_reco, QSq_true);
            FillHistogram(interaction.QSq_1Track_Error, QSq_error);
            FillHistogram(interaction.QSq_1Track_Diff, QSq_reco-QSq_true);
        }else{ 
            FillHistogramWithDefaultErrors(interaction.QSq_2Track_response, QSq_reco, QSq_true);
            FillHistogram(interaction.QSq_2Track_Error, QSq_error);
            FillHistogram(interaction.QSq_2Track_Diff, QSq_reco-QSq_true);
        }  
    }
}

void CCProtonPi0_Analyzer::fillInteractionReco()
{
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

    // Vertex & Extra Energy
    if (nProtonCandidates == 0){
        FillHistogram(interaction.vertex_energy_1Track, vertex_blob_energy);
        FillHistogram(interaction.vertex_evis_1Track, vertex_blob_evis );
        FillHistogram(interaction.extra_leftover_energy_1Track, Extra_Energy_Leftover); 
        FillHistogram(interaction.extra_muon_energy_1Track, Extra_Energy_Muon); 
        FillHistogram(interaction.extra_rejected_energy_1Track, Extra_Energy_Rejected); 
        FillHistogram(interaction.extra_total_energy_1Track, Extra_Energy_Total);
    }else{
        FillHistogram(interaction.vertex_evis_2Track, vertex_blob_evis );
        FillHistogram(interaction.vertex_energy_2Track, vertex_blob_energy);
        FillHistogram(interaction.extra_leftover_energy_2Track, Extra_Energy_Leftover); 
        FillHistogram(interaction.extra_muon_energy_2Track, Extra_Energy_Muon); 
        FillHistogram(interaction.extra_rejected_energy_2Track, Extra_Energy_Rejected); 
        FillHistogram(interaction.extra_total_energy_2Track, Extra_Energy_Total);
    }

    // Other Event Parameters
    if (nProtonCandidates > 0 && mc_intType == 2){
        FillHistogram(interaction.deltaInvMass, calcDeltaInvariantMass() * MeV_to_GeV);
    }

    // Extra Energy
    FillHistogram(interaction.h_extra_leftover_energy, Extra_Energy_Leftover); 
    FillHistogram(interaction.h_extra_muon_energy, Extra_Energy_Muon); 
    FillHistogram(interaction.h_extra_rejected_energy, Extra_Energy_Rejected); 

    
    if (m_isMC && truth_isSignal) FillSignalCharacteristics_Reco(); 
}

double CCProtonPi0_Analyzer::calcDeltaInvariantMass()
{
    double invMassSq;

    invMassSq = 
        (pi0_E + proton_E) * 
        (pi0_E + proton_E) -
        (   (pi0_px + proton_px) * 
            (pi0_px + proton_px) + 
            (pi0_py + proton_py) * 
            (pi0_py + proton_py) +
            (pi0_pz + proton_pz) * 
            (pi0_pz + proton_pz));

    return sqrt(invMassSq);
}

void CCProtonPi0_Analyzer::writeScanFile()
{
    if(isScanRun){

        // Constants for Roundup List
        const string arachne_html = "http://minerva05.fnal.gov/Arachne/arachne.html?filename=";
        const string entryString  = "&entry=";
        const string other        = "&slice=-1&filetype=dst";
        //http://minerva05.fnal.gov/Arachne/arachne.html?det=MV&recoVer=v10r6p13&run=3596&subrun=6&gate=597&slice=7
        roundupText<<arachne_html<<scanFileName<<entryString<<truth_eventID<<other<<" ";
        roundupText<<mc_run<<" ^ "<<mc_subrun<<" ^ "<<ev_gate<<" ^ "<<mc_incomingE<<endl;
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

    for (int i = 0; i < nTopologies; i++){
        failText.close();
    }
}

void CCProtonPi0_Analyzer::openTextFiles()
{
    cout<<"Opening Text Files:"<<endl;

    // Open Fail-Check File
    failFile = Folder_List::output + Folder_List::textOut + "FailChecks.txt";

    failText.open( failFile.c_str() );
    if( !failText.is_open() ){
        cerr<<"Cannot open output text file: "<<failFile<<endl;
        exit(1);
    }else{
        cout<<"\t"<<failFile<<endl;
    }

    if(isScanRun){
        // Open Roundup Text for Arachne Scanning
        string roundupFile = Folder_List::output + Folder_List::textOut + "ArachneRoundup.txt";
        roundupText.open(roundupFile.c_str() );
        if( !roundupText.is_open() ){
            cerr<<"Cannot open output text file: "<<roundupFile<<endl;
            exit(1);
        }else{
            cout<<"\t"<<roundupFile<<endl; 
        }

        string playlistDST = "Input/Playlists/pl_Scan.dat";
        DSTFileList.open( playlistDST.c_str() );
        if( !DSTFileList.is_open() ){
            cerr<<"Cannot open input text file: "<<playlistDST<<endl;
            exit(1);
        }else{
            cout<<"\t"<<playlistDST<<endl;
        }
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
        FillHistogramWithDefaultErrors(proton.proton_P_response, reco_P, true_P);

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
        FillHistogram(proton.theta_Diff, reco_theta-true_theta);
        FillHistogramWithDefaultErrors(proton.proton_theta_response, proton_theta_beam * TMath::RadToDeg(), truth_proton_theta_beam * TMath::RadToDeg());
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
        FillHistogram(pi0.theta_Diff, reco_theta-true_theta);
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
    if (nProtonCandidates == 0){
        FillHistogram(pi0.P_1Track, pi0_P * MeV_to_GeV);
        FillHistogram(pi0.theta_1Track, pi0_theta * TMath::RadToDeg());
    }else{
        FillHistogram(pi0.P_2Track, pi0_P * MeV_to_GeV);
        FillHistogram(pi0.theta_2Track, pi0_theta * TMath::RadToDeg());
    }
}

void CCProtonPi0_Analyzer::fillMuonMC()
{
    if(truth_isSignal){ 
        // Momentum
        double reco_P = muon_P * MeV_to_GeV;
        double true_P = truth_muon_P * MeV_to_GeV;
        double error_P = Data_Functions::getError(true_P, reco_P);

        FillHistogram(muon.reco_P_true_P, reco_P,true_P);
        FillHistogram(muon.P_error, error_P);

        // Energy
        double reco_E = muon_E * MeV_to_GeV;
        double true_E = truth_muon_4P[3] * MeV_to_GeV; 
        double error_E = Data_Functions::getError(true_E, reco_E);

        FillHistogram(muon.reco_E_true_E, reco_E,true_E);
        FillHistogram(muon.E_error, error_E);
        FillHistogram(muon.E_Diff, reco_E-true_E);

        // Theta
        double reco_theta = muon_theta_beam * TMath::RadToDeg();
        double true_theta = truth_muon_theta_beam * TMath::RadToDeg();
        double error_theta = Data_Functions::getError(true_theta, reco_theta);

        FillHistogram(muon.theta_error, error_theta);
        FillHistogram(muon.theta_Diff, reco_theta-true_theta);

        // thetaX
        double reco_thetaX = muon_thetaX_beam;
        double true_thetaX = truth_muon_thetaX_beam;
        FillHistogram(muon.thetaX_Diff, reco_thetaX-true_thetaX);
        FillHistogram(muon.thetaX_thetaX_test, reco_thetaX, true_thetaX);
 
        // thetaY
        double reco_thetaY = muon_thetaY_beam;
        double true_thetaY = truth_muon_thetaY_beam;
        FillHistogram(muon.thetaY_Diff, reco_thetaY-true_thetaY);
        FillHistogram(muon.thetaY_thetaY_test, reco_thetaY, true_thetaY);
        
        // Cosine Theta
        double reco_cos_theta = cos(muon_theta_beam);
        double true_cos_theta = cos(truth_muon_theta_beam);
        double error_cos_theta = Data_Functions::getError(true_cos_theta, reco_cos_theta);

        FillHistogram(muon.cos_theta_error, error_cos_theta);
 
        FillHistogramWithDefaultErrors(muon.theta_theta_test, muon_theta_beam * TMath::RadToDeg(), truth_muon_theta_beam * TMath::RadToDeg());
        
    }
}

void CCProtonPi0_Analyzer::fill_BackgroundSubtractionHists() 
{
    if (m_isMC){
        FillHistogramWithDefaultErrors(cutList.invMass_mc_reco_all, pi0_invMass);
        if (truth_isSignal)  FillHistogramWithDefaultErrors(cutList.invMass_mc_reco_signal, pi0_invMass);
        else FillHistogramWithDefaultErrors(cutList.invMass_mc_reco_bckg, pi0_invMass);
    }else{
        FillHistogram(cutList.invMass_all, pi0_invMass);
    }
}

void CCProtonPi0_Analyzer::fill_Enu() 
{
    if (m_isMC){
        // MC Reco All
        FillHistogramWithDefaultErrors(interaction.Enu_mc_reco_all, m_Enu * MeV_to_GeV);
        if (truth_isSignal){
            // MC Truth Signal
            FillHistogramWithDefaultErrors(interaction.Enu_mc_truth_signal, mc_incomingE * MeV_to_GeV);
            // MC Reco Signal
            FillHistogramWithDefaultErrors(interaction.Enu_mc_reco_signal, m_Enu * MeV_to_GeV);
            // MC Reco vs True -- Response
            FillHistogramWithDefaultErrors(interaction.Enu_response, m_Enu * MeV_to_GeV, mc_incomingE * MeV_to_GeV);
        }else{
            // MC Reco Background
            FillHistogramWithDefaultErrors(interaction.Enu_mc_reco_bckg, m_Enu * MeV_to_GeV);
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
        FillHistogramWithDefaultErrors(interaction.QSq_mc_reco_all, m_QSq * MeVSq_to_GeVSq);
        if (truth_isSignal){
            // MC Truth Signal
            FillHistogramWithDefaultErrors(interaction.QSq_mc_truth_signal, mc_Q2 * MeVSq_to_GeVSq);
            // MC Reco Signal
            FillHistogramWithDefaultErrors(interaction.QSq_mc_reco_signal, m_QSq * MeVSq_to_GeVSq);
            // MC Reco vs True -- Response
            FillHistogramWithDefaultErrors(interaction.QSq_response, m_QSq * MeVSq_to_GeVSq, mc_Q2 * MeVSq_to_GeVSq);
        }else{
            // MC Reco Background
            FillHistogramWithDefaultErrors(interaction.QSq_mc_reco_bckg, m_QSq * MeVSq_to_GeVSq);
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
        FillHistogramWithDefaultErrors(muon.muon_P_mc_reco_all, muon_P * MeV_to_GeV);
        if (truth_isSignal){
            // MC Truth Signal
            FillHistogramWithDefaultErrors(muon.muon_P_mc_truth_signal, truth_muon_P * MeV_to_GeV);
            // MC Reco Signal
            FillHistogramWithDefaultErrors(muon.muon_P_mc_reco_signal, muon_P * MeV_to_GeV);
            // MC Reco vs True -- Response
            FillHistogramWithDefaultErrors(muon.muon_P_response, muon_P * MeV_to_GeV, truth_muon_P * MeV_to_GeV);
        }else{
            // MC Reco Background
            FillHistogramWithDefaultErrors(muon.muon_P_mc_reco_bckg, muon_P * MeV_to_GeV);
        }
    }else{
        // Data
        FillHistogram(muon.muon_P_all, muon_P * MeV_to_GeV);
    }
}

void CCProtonPi0_Analyzer::fill_muon_theta() 
{
    if (m_isMC){
        // MC Reco All
        FillHistogramWithDefaultErrors(muon.muon_theta_mc_reco_all, muon_theta_beam * TMath::RadToDeg());
        if (truth_isSignal){
            // MC Truth Signal
            FillHistogramWithDefaultErrors(muon.muon_theta_mc_truth_signal, truth_muon_theta_beam * TMath::RadToDeg());
            // MC Reco Signal
            FillHistogramWithDefaultErrors(muon.muon_theta_mc_reco_signal, muon_theta_beam * TMath::RadToDeg());
            // MC Reco vs True -- Response
            FillHistogramWithDefaultErrors(muon.muon_theta_response, muon_theta_beam * TMath::RadToDeg(), truth_muon_theta_beam * TMath::RadToDeg());
        }else{
            // MC Reco Background
            FillHistogramWithDefaultErrors(muon.muon_theta_mc_reco_bckg, muon_theta_beam * TMath::RadToDeg());
        }
    }else{
        // Data
        FillHistogram(muon.muon_theta_all, muon_theta_beam * TMath::RadToDeg());
    }
}

void CCProtonPi0_Analyzer::fill_muon_cos_theta() 
{
    if (m_isMC){
        // MC Reco All
        FillHistogramWithDefaultErrors(muon.muon_cos_theta_mc_reco_all, cos(muon_theta_beam));
        if (truth_isSignal){
            // MC Truth Signal
            FillHistogramWithDefaultErrors(muon.muon_cos_theta_mc_truth_signal, cos(truth_muon_theta_beam));
            // MC Reco Signal
            FillHistogramWithDefaultErrors(muon.muon_cos_theta_mc_reco_signal, cos(muon_theta_beam));
            // MC Reco vs True -- Response
            FillHistogramWithDefaultErrors(muon.muon_cos_theta_response, cos(muon_theta_beam), cos(truth_muon_theta_beam));
        }else{
            // MC Reco Background
            FillHistogramWithDefaultErrors(muon.muon_cos_theta_mc_reco_bckg, cos(muon_theta_beam));
        }
    }else{
        // Data
        FillHistogram(muon.muon_cos_theta_all, cos(muon_theta_beam));
    }
}

void CCProtonPi0_Analyzer::fill_pi0_P() 
{
    if (m_isMC){
        // MC Reco All
        FillHistogramWithDefaultErrors(pi0.pi0_P_mc_reco_all, pi0_P * MeV_to_GeV);
        if (truth_isSignal){
            // MC Truth Signal
            FillHistogramWithDefaultErrors(pi0.pi0_P_mc_truth_signal, truth_pi0_P * MeV_to_GeV);
            // MC Reco Signal
            FillHistogramWithDefaultErrors(pi0.pi0_P_mc_reco_signal, pi0_P * MeV_to_GeV);
            // MC Reco vs True -- Response
            FillHistogramWithDefaultErrors(pi0.pi0_P_response, pi0_P * MeV_to_GeV, truth_pi0_P * MeV_to_GeV);
        }else{
            // MC Reco Background
            FillHistogramWithDefaultErrors(pi0.pi0_P_mc_reco_bckg, pi0_P * MeV_to_GeV);
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
        FillHistogramWithDefaultErrors(pi0.pi0_KE_mc_reco_all, pi0_KE * MeV_to_GeV);
        if (truth_isSignal){
            // MC Truth Signal
            FillHistogramWithDefaultErrors(pi0.pi0_KE_mc_truth_signal, truth_pi0_KE  * MeV_to_GeV);
            // MC Reco Signal
            FillHistogramWithDefaultErrors(pi0.pi0_KE_mc_reco_signal, pi0_KE * MeV_to_GeV);
            // MC Reco vs True -- Response
            FillHistogramWithDefaultErrors(pi0.pi0_KE_response, pi0_KE * MeV_to_GeV, truth_pi0_KE * MeV_to_GeV);
        }else{
            // MC Reco Background
            FillHistogramWithDefaultErrors(pi0.pi0_KE_mc_reco_bckg, pi0_KE * MeV_to_GeV);
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
        FillHistogramWithDefaultErrors(pi0.pi0_theta_mc_reco_all, pi0_theta_beam * TMath::RadToDeg());
        if (truth_isSignal){
            // MC Truth Signal
            FillHistogramWithDefaultErrors(pi0.pi0_theta_mc_truth_signal, truth_pi0_theta_beam  * TMath::RadToDeg());
            // MC Reco Signal
            FillHistogramWithDefaultErrors(pi0.pi0_theta_mc_reco_signal, pi0_theta_beam * TMath::RadToDeg());
            // MC Reco vs True -- Response
            FillHistogramWithDefaultErrors(pi0.pi0_theta_response, pi0_theta_beam * TMath::RadToDeg(), truth_pi0_theta_beam * TMath::RadToDeg());
        }else{
            // MC Reco Background
            FillHistogramWithDefaultErrors(pi0.pi0_theta_mc_reco_bckg, pi0_theta_beam * TMath::RadToDeg());
        }
    }else{
        // Data
        FillHistogram(pi0.pi0_theta_all, pi0_theta_beam * TMath::RadToDeg());
    }
}

void CCProtonPi0_Analyzer::fillMuonReco()
{
    FillHistogram(muon.E, muon_E * MeV_to_GeV);
    FillHistogram(muon.P, muon_P * MeV_to_GeV);
    FillHistogram(muon.KE, muon_KE * MeV_to_GeV);
    FillHistogram(muon.theta, muon_theta_beam * TMath::RadToDeg());
    FillHistogram(muon.cos_theta, cos(muon_theta_beam));
    FillHistogram(muon.phi, muon_phi * TMath::RadToDeg());
    if (nProtonCandidates == 0){
        FillHistogram(muon.P_1Track, muon_P * MeV_to_GeV);
        FillHistogram(muon.theta_1Track, muon_P * MeV_to_GeV);
    }else{
        FillHistogram(muon.P_2Track, muon_P * MeV_to_GeV);
        FillHistogram(muon.theta_2Track, muon_P * MeV_to_GeV);
    }
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

void CCProtonPi0_Analyzer::FillHistogramWithDefaultErrors(MnvH1D* hist, double var)
{
    hist->Fill(var, cvweight);
    FillVertErrorBand_Flux(hist, var);
    FillVertErrorBand_Genie(hist, var);
    FillVertErrorBand_MuonTracking(hist, var);
}

void CCProtonPi0_Analyzer::FillHistogram(MnvH1D* hist, double var)
{
    hist->Fill(var, cvweight);
}

void CCProtonPi0_Analyzer::FillHistogramWithDefaultErrors(MnvH2D* hist, double var1, double var2)
{
    hist->Fill(var1,var2, cvweight);
    FillVertErrorBand_Flux(hist, var1, var2);
    FillVertErrorBand_Genie(hist, var1, var2);
    FillVertErrorBand_MuonTracking(hist, var1, var2);
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
        }else{
            // Fill Background
            hist[2]->Fill(var, cvweight); // Always Fill ind == 2 -- All Background

            // Fill Background Type
            int ind = GetBackgroundTypeInd();
            hist[ind]->Fill(var, cvweight);
        }
    }
}

void CCProtonPi0_Analyzer::FillVertErrorBand_MuonTracking(MnvH1D* h, double var)
{
    double correctionErr = GetMINOSCorrectionErr();
    h->FillVertErrorBand("MuonTracking", var, 1-correctionErr, 1+correctionErr, cvweight);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_MuonTracking(MnvH2D* h, double var1, double var2)
{
    double correctionErr = GetMINOSCorrectionErr();
    h->FillVertErrorBand("MuonTracking", var1, var2, 1-correctionErr, 1+correctionErr, cvweight);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_Flux(MnvH1D* h, double var)
{
    std::vector<double> flux_errors = GetFluxError();
    h->FillVertErrorBand("Flux",  var, &flux_errors[0],  cvweight, 1.0);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_Flux(MnvH2D* h, double var1, double var2)
{
    std::vector<double> flux_errors = GetFluxError();
    h->FillVertErrorBand("Flux",  var1, var2,  &flux_errors[0],  cvweight, 1.0);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_Genie(MnvH1D* h, double var)
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
    h->FillVertErrorBand("GENIE_MaCCQEshape"       ,var, truth_genie_wgt_MaCCQEshape[2]      , truth_genie_wgt_MaCCQEshape[4]      , cvweight);
    h->FillVertErrorBand("GENIE_MaNCEL"            ,var, truth_genie_wgt_MaNCEL[2]           , truth_genie_wgt_MaNCEL[4]           , cvweight);
    h->FillVertErrorBand("GENIE_MaRES"             ,var, truth_genie_wgt_MaRES[2]            , truth_genie_wgt_MaRES[4]            , cvweight);
    h->FillVertErrorBand("GENIE_MvRES"             ,var, truth_genie_wgt_MvRES[2]            , truth_genie_wgt_MvRES[4]            , cvweight);
    h->FillVertErrorBand("GENIE_NormCCQE"          ,var, truth_genie_wgt_NormCCQE[2]         , truth_genie_wgt_NormCCQE[4]         , cvweight);
    h->FillVertErrorBand("GENIE_NormCCRES"         ,var, truth_genie_wgt_NormCCRES[2]        , truth_genie_wgt_NormCCRES[4]        , cvweight);
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

void CCProtonPi0_Analyzer::FillVertErrorBand_Genie(MnvH2D* h, double var1, double var2)
{
    h->FillVertErrorBand("GENIE_AGKYxF1pi"         ,var1, var2, truth_genie_wgt_AGKYxF1pi[2]        , truth_genie_wgt_AGKYxF1pi[4]        , cvweight);
    h->FillVertErrorBand("GENIE_AhtBY"             ,var1, var2, truth_genie_wgt_AhtBY[2]            , truth_genie_wgt_AhtBY[4]            , cvweight);
    h->FillVertErrorBand("GENIE_BhtBY"             ,var1, var2, truth_genie_wgt_BhtBY[2]            , truth_genie_wgt_BhtBY[4]            , cvweight);
    h->FillVertErrorBand("GENIE_CCQEPauliSupViaKF" ,var1, var2, truth_genie_wgt_CCQEPauliSupViaKF[2], truth_genie_wgt_CCQEPauliSupViaKF[4], cvweight);
    h->FillVertErrorBand("GENIE_CV1uBY"            ,var1, var2, truth_genie_wgt_CV1uBY[2]           , truth_genie_wgt_CV1uBY[4]           , cvweight);
    h->FillVertErrorBand("GENIE_CV2uBY"            ,var1, var2, truth_genie_wgt_CV2uBY[2]           , truth_genie_wgt_CV2uBY[4]           , cvweight);
    h->FillVertErrorBand("GENIE_EtaNCEL"           ,var1, var2, truth_genie_wgt_EtaNCEL[2]          , truth_genie_wgt_EtaNCEL[4]          , cvweight);
    h->FillVertErrorBand("GENIE_FrAbs_N"           ,var1, var2, truth_genie_wgt_FrAbs_N[2]          , truth_genie_wgt_FrAbs_N[4]          , cvweight);
    h->FillVertErrorBand("GENIE_FrAbs_pi"          ,var1, var2, truth_genie_wgt_FrAbs_pi[2]         , truth_genie_wgt_FrAbs_pi[4]         , cvweight);
    h->FillVertErrorBand("GENIE_FrCEx_N"           ,var1, var2, truth_genie_wgt_FrCEx_N[2]          , truth_genie_wgt_FrCEx_N[4]          , cvweight);
    h->FillVertErrorBand("GENIE_FrCEx_pi"          ,var1, var2, truth_genie_wgt_FrCEx_pi[2]         , truth_genie_wgt_FrCEx_pi[4]         , cvweight);
    h->FillVertErrorBand("GENIE_FrElas_N"          ,var1, var2, truth_genie_wgt_FrElas_N[2]         , truth_genie_wgt_FrElas_N[4]         , cvweight);
    h->FillVertErrorBand("GENIE_FrElas_pi"         ,var1, var2, truth_genie_wgt_FrElas_pi[2]        , truth_genie_wgt_FrElas_pi[4]        , cvweight);
    h->FillVertErrorBand("GENIE_FrInel_N"          ,var1, var2, truth_genie_wgt_FrInel_N[2]         , truth_genie_wgt_FrInel_N[4]         , cvweight);
    h->FillVertErrorBand("GENIE_FrInel_pi"         ,var1, var2, truth_genie_wgt_FrInel_pi[2]        , truth_genie_wgt_FrInel_pi[4]        , cvweight);
    h->FillVertErrorBand("GENIE_FrPiProd_N"        ,var1, var2, truth_genie_wgt_FrPiProd_N[2]       , truth_genie_wgt_FrPiProd_N[4]       , cvweight);
    h->FillVertErrorBand("GENIE_FrPiProd_pi"       ,var1, var2, truth_genie_wgt_FrPiProd_pi[2]      , truth_genie_wgt_FrPiProd_pi[4]      , cvweight);
    h->FillVertErrorBand("GENIE_MFP_N"             ,var1, var2, truth_genie_wgt_MFP_N[2]            , truth_genie_wgt_MFP_N[4]            , cvweight);
    h->FillVertErrorBand("GENIE_MFP_pi"            ,var1, var2, truth_genie_wgt_MFP_pi[2]           , truth_genie_wgt_MFP_pi[4]           , cvweight);
    h->FillVertErrorBand("GENIE_MaCCQE"            ,var1, var2, truth_genie_wgt_MaCCQE[2]           , truth_genie_wgt_MaCCQE[4]           , cvweight);
    h->FillVertErrorBand("GENIE_MaCCQEshape"       ,var1, var2, truth_genie_wgt_MaCCQEshape[2]      , truth_genie_wgt_MaCCQEshape[4]      , cvweight);
    h->FillVertErrorBand("GENIE_MaNCEL"            ,var1, var2, truth_genie_wgt_MaNCEL[2]           , truth_genie_wgt_MaNCEL[4]           , cvweight);
    h->FillVertErrorBand("GENIE_MaRES"             ,var1, var2, truth_genie_wgt_MaRES[2]            , truth_genie_wgt_MaRES[4]            , cvweight);
    h->FillVertErrorBand("GENIE_MvRES"             ,var1, var2, truth_genie_wgt_MvRES[2]            , truth_genie_wgt_MvRES[4]            , cvweight);
    h->FillVertErrorBand("GENIE_NormCCQE"          ,var1, var2, truth_genie_wgt_NormCCQE[2]         , truth_genie_wgt_NormCCQE[4]         , cvweight);
    h->FillVertErrorBand("GENIE_NormCCRES"         ,var1, var2, truth_genie_wgt_NormCCRES[2]        , truth_genie_wgt_NormCCRES[4]        , cvweight);
    h->FillVertErrorBand("GENIE_NormDISCC"         ,var1, var2, truth_genie_wgt_NormDISCC[2]        , truth_genie_wgt_NormDISCC[4]        , cvweight);
    h->FillVertErrorBand("GENIE_NormNCRES"         ,var1, var2, truth_genie_wgt_NormNCRES[2]        , truth_genie_wgt_NormNCRES[4]        , cvweight);
    h->FillVertErrorBand("GENIE_RDecBR1gamma"      ,var1, var2, truth_genie_wgt_RDecBR1gamma[2]     , truth_genie_wgt_RDecBR1gamma[4]     , cvweight);
    h->FillVertErrorBand("GENIE_Rvn1pi"            ,var1, var2, truth_genie_wgt_Rvn1pi[2]           , truth_genie_wgt_Rvn1pi[4]           , cvweight);
    h->FillVertErrorBand("GENIE_Rvn2pi"            ,var1, var2, truth_genie_wgt_Rvn2pi[2]           , truth_genie_wgt_Rvn2pi[4]           , cvweight);
    h->FillVertErrorBand("GENIE_Rvp1pi"            ,var1, var2, truth_genie_wgt_Rvp1pi[2]           , truth_genie_wgt_Rvp1pi[4]           , cvweight);
    h->FillVertErrorBand("GENIE_Rvp2pi"            ,var1, var2, truth_genie_wgt_Rvp2pi[2]           , truth_genie_wgt_Rvp2pi[4]           , cvweight);
    h->FillVertErrorBand("GENIE_Theta_Delta2Npi"   ,var1, var2, truth_genie_wgt_Theta_Delta2Npi[2]  , truth_genie_wgt_Theta_Delta2Npi[4]  , cvweight);
    h->FillVertErrorBand("GENIE_VecFFCCQEshape"    ,var1, var2, truth_genie_wgt_VecFFCCQEshape[2]   , truth_genie_wgt_VecFFCCQEshape[4]   , cvweight);
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

double CCProtonPi0_Analyzer::GetVertexEnergy()
{
    double energy;
    if (nProtonCandidates == 0){
        energy = vertex_blob_evis*1.01 + 80.9;
    }else{
        energy = vertex_blob_evis + 29.78;
    }

    return energy;
}

double CCProtonPi0_Analyzer::Calc_TruePi0OpeningAngle()
{
    double g1_P = HEP_Functions::calcMomentum(truth_gamma1_4P[0],truth_gamma1_4P[1],truth_gamma1_4P[2]);
    double g2_P = HEP_Functions::calcMomentum(truth_gamma2_4P[0],truth_gamma2_4P[1],truth_gamma2_4P[2]);
    double g1_g2 = truth_gamma1_4P[0]*truth_gamma2_4P[0] + truth_gamma1_4P[1]*truth_gamma2_4P[1] + truth_gamma1_4P[2]*truth_gamma2_4P[2];
    double cos_theta = g1_g2 / (g1_P*g2_P);
    double theta = acos(cos_theta) * rad_to_deg;

    return theta;
}

// Function Reserved for Correcting a NTupleVariables
void CCProtonPi0_Analyzer::CorrectNTupleVariables()
{
    // Do Nothing
}

void CCProtonPi0_Analyzer::CorrectEMShowerCalibration()
{
    // Found the correction using a Double Gaussian Fit in MATLAB
    const double pi0_mass = 134.98;
    double correction;
    if (m_isMC) correction = pi0_mass/130.28;
    else correction = pi0_mass/137.76;

    // Pi0 Variables
    pi0_invMass *= correction;
    pi0_E *= correction;
    pi0_E_Cal *= correction;
    pi0_P *= correction;
    pi0_KE = pi0_E - pi0_mass;
    pi0_px *= correction;
    pi0_py *= correction;
    pi0_pz *= correction;

    // gamma1 
    gamma1_E *= correction;
    gamma1_E_Old *= correction;
    gamma1_P *= correction;
    gamma1_px *= correction;
    gamma1_py *= correction;
    gamma1_pz *= correction;

    // gamma2 
    gamma2_E *= correction;
    gamma2_E_Old *= correction;
    gamma2_P *= correction;
    gamma2_px *= correction;
    gamma2_py *= correction;
    gamma2_pz *= correction;

}

void CCProtonPi0_Analyzer::Calc_WeightFromSystematics()
{
    if (m_isMC){
        // No need to Update FluxReweighter for minerva1 and after minerva13
        UpdateFluxReweighter(mc_run); 
        
        // Replace cvweight with Flux Weight
        cvweight = GetFluxWeight();

        // Update cvweight with MINOS Efficiency Correction
        const double minos_eff_correction = GetMINOSCorrection();
        cvweight *= minos_eff_correction;
        
        // Apply Background Constraints
        if (applyBckgConstraints){
            if (truth_isBckg_Compact_WithPi0) cvweight *= 0.90;
            else if (truth_isBckg_Compact_QELike) cvweight *= 0.89;
            else if (truth_isBckg_Compact_SinglePiPlus) cvweight *= 0.97;
        }
    }else{
        cvweight = 1.0; 
    }
}

void CCProtonPi0_Analyzer::AddErrorBands_Data()
{
    AddVertErrorBands_Data(pi0.pi0_P_all);
    AddVertErrorBands_Data(pi0.pi0_KE_all);
    AddVertErrorBands_Data(pi0.pi0_theta_all);
    AddVertErrorBands_Data(muon.muon_P_all);
    AddVertErrorBands_Data(muon.muon_theta_all);
    AddVertErrorBands_Data(interaction.QSq_all);
}


void CCProtonPi0_Analyzer::UpdateSignalDef()
{
    // Signal Definition with Neutrino Energy
    if (m_isMC && truth_isSignal){
        bool isEnu_inRange = mc_incomingE >= min_Enu && mc_incomingE <= max_Enu; 
        truth_isSignal = isEnu_inRange;
        // If event no longer a signal due to Enu Range -- it is background
        if (!truth_isSignal){
            truth_isBckg_Compact_WithPi0 = true;
            truth_isBckg_SinglePi0 = true;
            truth_isBckg_Other = true;
        }
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
    std::string playlist = GetPlaylist(mc_run);
    MnvNormalizer normalizer("Eroica", playlist);
    double correctionErr = normalizer.GetCorrectionErr(CCProtonPi0_minos_trk_p);
    return correctionErr;
}

double CCProtonPi0_Analyzer::GetMINOSCorrection()
{
    std::string playlist = GetPlaylist(mc_run);
    MnvNormalizer normalizer("Eroica", playlist);
    double correction = normalizer.GetCorrection(CCProtonPi0_minos_trk_p);
    return correction;
}

double CCProtonPi0_Analyzer::GetFluxWeight()
{
    double Enu = mc_incomingE * MeV_to_GeV; // true neutrino energy (GeV)
    int nuPDG = mc_incoming; //neutrino PDG code

    double flux_weight = frw->GetFluxCVWeight(Enu, nuPDG);

    return flux_weight;
}

std::vector<double> CCProtonPi0_Analyzer::GetFluxError()
{
    double Enu = mc_incomingE * MeV_to_GeV; // true neutrino energy (GeV)
    int nuPDG = mc_incoming; //neutrino PDG code

    std::vector<double> flux_error = frw->GetFluxErrorWeights(Enu, nuPDG);

    return flux_error;
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

void CCProtonPi0_Analyzer::ResetCounters()
{
    counter1.count = 0.0;
    counter2.count = 0.0;
    counter3.count = 0.0;
    counter4.count = 0.0;
}


double CCProtonPi0_Analyzer::Calc_Enu() const
{
    double Enu;
    if (nProtonCandidates == 0){
        Enu = muon_E + pi0_E + vertex_blob_energy + Extra_Energy_Total;
    }else{
        Enu = muon_E + proton_KE + pi0_E + vertex_blob_energy + Extra_Energy_Total;
    }

    return Enu;
}

double CCProtonPi0_Analyzer::Calc_QSq(const double Enu) const
{
    const double Mmu = 105.66;  // Muon Rest Mass [MeV]
    TLorentzVector beam_4P = Get_Neutrino_4P(Enu);
    TLorentzVector muon_4P(muon_px, muon_py, muon_pz, muon_E);

    double qSq = (Mmu * Mmu) - 2 * beam_4P.Dot(muon_4P);
    double QSq = -qSq;

    return QSq;
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

double CCProtonPi0_Analyzer::Calc_WSq(const double Enu, const double QSq) const
{
    const double Mn = 939.57;    // Neutron Rest Mass [MeV]
    double Emu = muon_E;

    // Calculate WSq - Use eq. in Research Logbook page 31
    double WSq = Mn*Mn + 2*Mn*(Enu - Emu) - QSq; 

    return WSq;
}

void CCProtonPi0_Analyzer::Calc_EventKinematics()
{
    m_Enu = Calc_Enu();
    m_QSq = Calc_QSq(m_Enu);
    m_WSq = Calc_WSq(m_Enu, m_QSq);
    if (m_WSq > 0) m_W = sqrt(m_WSq);
    else m_W = 0;

    failText<<m_WSq*MeVSq_to_GeVSq<<" "<<m_Enu*MeV_to_GeV<<" ";
    failText<<m_QSq*MeVSq_to_GeVSq<<" "<<mc_Q2*MeVSq_to_GeVSq<<" ";
    failText<<muon_E*MeV_to_GeV<<" "<<muon_P*MeV_to_GeV<<" ";
    failText<<muon_theta_beam*TMath::RadToDeg()<<" "<<truth_muon_theta_beam*TMath::RadToDeg()<<endl;

    if (m_WSq > 0){
        if (nProtonCandidates == 0) counter1.count++;
        else counter2.count++; 
    }else{
        if (nProtonCandidates == 0) counter3.count++;
        else counter4.count++; 
    }
}

void CCProtonPi0_Analyzer::FillSignalCharacteristics_Reco()
{
// Signal Characteristics
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
        }else if (mc_intType == 3){
            int nFS_pions = Get_nFS_pions();
            if (m_W * MeV_to_GeV < 1.7){
                FillHistogram(interaction.reco_w_Non_RES, m_W * MeV_to_GeV);
            }else if (nFS_pions == 1){
                FillHistogram(interaction.reco_w_DIS_1_pi, m_W * MeV_to_GeV);
            }else if (nFS_pions == 2){
                FillHistogram(interaction.reco_w_DIS_2_pi, m_W * MeV_to_GeV);
            }else if (nFS_pions > 2){
                FillHistogram(interaction.reco_w_DIS_Multi_pi, m_W * MeV_to_GeV);
            }
        }else{
            std::cout<<"WARNING! Signal Event with different interaction Type!"<<std::endl;
        }
}

void CCProtonPi0_Analyzer::FillSignalCharacteristics(bool isMinosMatched)
{
    if (isMinosMatched)
    {
        // Signal Characteristics
        if (mc_intType == 1){
            FillHistogram(cutList.mc_Q2_QE, mc_Q2 * MeVSq_to_GeVSq);
            FillHistogram(cutList.mc_incomingE_QE, mc_incomingE * MeV_to_GeV);
            FillHistogram(cutList.mc_w_QE, mc_w * MeV_to_GeV);
        }else if (mc_intType == 2){
            if (mc_resID == 0){ 
                FillHistogram(cutList.mc_Q2_RES_1232, mc_Q2 * MeVSq_to_GeVSq);
                FillHistogram(cutList.mc_incomingE_RES_1232, mc_incomingE * MeV_to_GeV);
                FillHistogram(cutList.mc_w_RES_1232, mc_w * MeV_to_GeV);
            }else if (mc_resID == 1){
                FillHistogram(cutList.mc_Q2_RES_1535, mc_Q2 * MeVSq_to_GeVSq);
                FillHistogram(cutList.mc_incomingE_RES_1535, mc_incomingE * MeV_to_GeV);
                FillHistogram(cutList.mc_w_RES_1535, mc_w * MeV_to_GeV);
            }else if (mc_resID == 2){
                FillHistogram(cutList.mc_Q2_RES_1520, mc_Q2 * MeVSq_to_GeVSq);
                FillHistogram(cutList.mc_incomingE_RES_1520, mc_incomingE * MeV_to_GeV);
                FillHistogram(cutList.mc_w_RES_1520, mc_w * MeV_to_GeV);
            }else{
                FillHistogram(cutList.mc_Q2_RES_Other, mc_Q2 * MeVSq_to_GeVSq);
                FillHistogram(cutList.mc_incomingE_RES_Other, mc_incomingE * MeV_to_GeV);
                FillHistogram(cutList.mc_w_RES_Other, mc_w * MeV_to_GeV);
            }
        }else if (mc_intType == 3){
            int nFS_pions = Get_nFS_pions();
            if (mc_w * MeV_to_GeV < 1.7){
                FillHistogram(cutList.mc_Q2_Non_RES, mc_Q2 * MeVSq_to_GeVSq);
                FillHistogram(cutList.mc_incomingE_Non_RES, mc_incomingE * MeV_to_GeV);
                FillHistogram(cutList.mc_w_Non_RES, mc_w * MeV_to_GeV);
            }else if (nFS_pions == 1){
                FillHistogram(cutList.mc_Q2_DIS_1_pi, mc_Q2 * MeVSq_to_GeVSq);
                FillHistogram(cutList.mc_incomingE_DIS_1_pi, mc_incomingE * MeV_to_GeV);
                FillHistogram(cutList.mc_w_DIS_1_pi, mc_w * MeV_to_GeV);
            }else if (nFS_pions == 2){
                FillHistogram(cutList.mc_Q2_DIS_2_pi, mc_Q2 * MeVSq_to_GeVSq);
                FillHistogram(cutList.mc_incomingE_DIS_2_pi, mc_incomingE * MeV_to_GeV);
                FillHistogram(cutList.mc_w_DIS_2_pi, mc_w * MeV_to_GeV);
            }else if (nFS_pions > 2){
                FillHistogram(cutList.mc_Q2_DIS_Multi_pi, mc_Q2 * MeVSq_to_GeVSq);
                FillHistogram(cutList.mc_incomingE_DIS_Multi_pi, mc_incomingE * MeV_to_GeV);
                FillHistogram(cutList.mc_w_DIS_Multi_pi, mc_w * MeV_to_GeV);
            }
        }else{
            std::cout<<"WARNING! Signal Event with different interaction Type!"<<std::endl;
        }
    }else{
        // Signal Characteristics
        if (mc_intType == 1){
            FillHistogram(interaction.mc_Q2_QE, mc_Q2 * MeVSq_to_GeVSq);
            FillHistogram(interaction.mc_incomingE_QE, mc_incomingE * MeV_to_GeV);
            FillHistogram(interaction.mc_w_QE, mc_w * MeV_to_GeV);
        }else if (mc_intType == 2){
            if (mc_resID == 0){ 
                FillHistogram(interaction.mc_Q2_RES_1232, mc_Q2 * MeVSq_to_GeVSq);
                FillHistogram(interaction.mc_incomingE_RES_1232, mc_incomingE * MeV_to_GeV);
                FillHistogram(interaction.mc_w_RES_1232, mc_w * MeV_to_GeV);
            }else if (mc_resID == 1){
                FillHistogram(interaction.mc_Q2_RES_1535, mc_Q2 * MeVSq_to_GeVSq);
                FillHistogram(interaction.mc_incomingE_RES_1535, mc_incomingE * MeV_to_GeV);
                FillHistogram(interaction.mc_w_RES_1535, mc_w * MeV_to_GeV);
            }else if (mc_resID == 2){
                FillHistogram(interaction.mc_Q2_RES_1520, mc_Q2 * MeVSq_to_GeVSq);
                FillHistogram(interaction.mc_incomingE_RES_1520, mc_incomingE * MeV_to_GeV);
                FillHistogram(interaction.mc_w_RES_1520, mc_w * MeV_to_GeV);
            }else{
                FillHistogram(interaction.mc_Q2_RES_Other, mc_Q2 * MeVSq_to_GeVSq);
                FillHistogram(interaction.mc_incomingE_RES_Other, mc_incomingE * MeV_to_GeV);
                FillHistogram(interaction.mc_w_RES_Other, mc_w * MeV_to_GeV);
            }
        }else if (mc_intType == 3){
            int nFS_pions = Get_nFS_pions();
            if (mc_w * MeV_to_GeV < 1.7){
                FillHistogram(interaction.mc_Q2_Non_RES, mc_Q2 * MeVSq_to_GeVSq);
                FillHistogram(interaction.mc_incomingE_Non_RES, mc_incomingE * MeV_to_GeV);
                FillHistogram(interaction.mc_w_Non_RES, mc_w * MeV_to_GeV);
            }else if (nFS_pions == 1){
                FillHistogram(interaction.mc_Q2_DIS_1_pi, mc_Q2 * MeVSq_to_GeVSq);
                FillHistogram(interaction.mc_incomingE_DIS_1_pi, mc_incomingE * MeV_to_GeV);
                FillHistogram(interaction.mc_w_DIS_1_pi, mc_w * MeV_to_GeV);
            }else if (nFS_pions == 2){
                FillHistogram(interaction.mc_Q2_DIS_2_pi, mc_Q2 * MeVSq_to_GeVSq);
                FillHistogram(interaction.mc_incomingE_DIS_2_pi, mc_incomingE * MeV_to_GeV);
                FillHistogram(interaction.mc_w_DIS_2_pi, mc_w * MeV_to_GeV);
            }else if (nFS_pions > 2){
                FillHistogram(interaction.mc_Q2_DIS_Multi_pi, mc_Q2 * MeVSq_to_GeVSq);
                FillHistogram(interaction.mc_incomingE_DIS_Multi_pi, mc_incomingE * MeV_to_GeV);
                FillHistogram(interaction.mc_w_DIS_Multi_pi, mc_w * MeV_to_GeV);
            }
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
        FillHistogram(cutList.hCut_pi0invMass, pi0_invMass);

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
        FillHistogram(cutList.SideBand_muon_P, muon_P*MeV_to_GeV);
        FillHistogram(cutList.SideBand_muon_theta, muon_theta_beam*TMath::RadToDeg());
        FillHistogram(cutList.SideBand_pi0_P, pi0_P*MeV_to_GeV);
        FillHistogram(cutList.SideBand_pi0_KE, pi0_KE*MeV_to_GeV);
        FillHistogram(cutList.SideBand_pi0_theta, pi0_theta_beam*TMath::RadToDeg());
        FillHistogram(cutList.SideBand_neutrino_E, m_Enu*MeV_to_GeV);
        FillHistogram(cutList.SideBand_QSq, m_QSq*MeVSq_to_GeVSq);
    }
}



#endif //CCProtonPi0_Analyzer_cpp


