/*
 * See CCProtonPi0_Analyzer.h header or Class Information
 */
#ifndef CCProtonPi0_Analyzer_cpp
#define CCProtonPi0_Analyzer_cpp

// Trung's Implementations
#include "new_flux.h"

#include "CCProtonPi0_Analyzer.h"

using namespace std;

void CCProtonPi0_Analyzer::specifyRunTime()
{
    applyMaxEvents = false;
    nMaxEvents = 100000;

    // Initialize Flux File if it is MC 
    if (m_isMC){
        new_flux::get().read_oldflux_histogram(Folder_List::rootDir_Flux_old);
        new_flux::get().read_newflux_histogram(Folder_List::rootDir_Flux_new);
        new_flux::get().calc_weights();
    }
 
    // Control Flow
    isDataAnalysis  = true;
    isScanRun = true;
    writeFSParticleMomentum = false;
   
    // Event Selections
    applyProtonScore = true;
    minProtonScore_LLR = -10.0;

    applyPhotonDistance = true;
    minPhotonDistance_1 = 14; //cm
    minPhotonDistance_2 = 0; //cm

    applyBeamEnergy = true;
    max_beamEnergy = 20.0; // GeV

    min_Pi0_invMass = 60.0;
    max_Pi0_invMass = 200.0;

    applyDeltaInvMass = false;
    min_Delta_invMass = 40.0;
    max_Delta_invMass = 200.0;

    latest_ScanID = 0.0;

    counter1 = 0;
    counter2 = 0;

    cvweight = 1.0;
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

        Calc_WeightFromSystematics();

        CorrectEMShowerCalibration();

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
    cout<<"counter1 = "<<counter1<<endl;
    cout<<"counter2 = "<<counter2<<endl;
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
       
        CorrectNTupleVariables();

        // Update scanFileName if running for scan
        if(isScanRun) UpdateScanFileName();

        //----------------------------------------------------------------------
        // Fill Background Branches for Background Events
        //----------------------------------------------------------------------
        if(m_isMC && !truth_isSignal) {
            bckgTool.fillBackgroundWithPi0(truth_isBckg_NoPi0, truth_isBckg_SinglePi0, truth_isBckg_MultiPi0, truth_isBckg_withMichel);                                    
            bckgTool.fillBackground(truth_isBckg_NC, truth_isBckg_AntiNeutrino, truth_isBckg_QELike, truth_isBckg_SingleChargedPion,truth_isBckg_SingleChargedPion_ChargeExchanged, truth_isBckg_DoublePionWithPi0, truth_isBckg_DoublePionWithoutPi0, truth_isBckg_MultiPionWithPi0, truth_isBckg_MultiPionWithoutPi0, truth_isBckg_Other, truth_isBckg_withMichel);                                    
        }

        Calc_WeightFromSystematics();
        
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
    cout<<"counter1 = "<<counter1<<endl;
    cout<<"counter2 = "<<counter2<<endl;
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
    normalizer("Eroica","minerva1")
{   
    cout<<"Initializing CCProtonPi0_Analyzer"<<endl;

    m_isMC = isMC;

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
    fill_muon_P();
    fill_muon_theta();

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

    // Fill Truth_W for MINOS Matched Signal Events
    if (m_isMC && truth_isSignal){ 
        fill_mc_w(); 
    }

    // ------------------------------------------------------------------------
    // Anti-Muon Cut
    // ------------------------------------------------------------------------
    if( Cut_Muon_Charge == 1) return false;
    cutList.nCut_Muon_Charge.increment(truth_isSignal, study1, study2);


    // ------------------------------------------------------------------------
    // Michel Cuts
    // ------------------------------------------------------------------------
    if( Cut_Vertex_Michel_Exist == 1 || Cut_EndPoint_Michel_Exist == 1 || Cut_secEndPoint_Michel_Exist == 1 ){
        FillHistogram(cutList.hCut_Michel,1);
    }else{
        FillHistogram(cutList.hCut_Michel,0);
    } 
    if( Cut_Vertex_Michel_Exist == 1) return false;
    cutList.nCut_Vertex_Michel_Exist.increment(truth_isSignal, study1, study2);

    if( Cut_EndPoint_Michel_Exist == 1) return false;
    cutList.nCut_EndPoint_Michel_Exist.increment(truth_isSignal, study1, study2);

    if( Cut_secEndPoint_Michel_Exist == 1) return false;
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
    if (nProtonCandidates > 0){
        for( int i = 0; i < nProtonCandidates; i++){
            if ( applyProtonScore ){
                FillHistogram(cutList.hCut_2Track_protonScore_LLR,all_protons_LLRScore[i]);
                if ( all_protons_LLRScore[i] < minProtonScore_LLR ) return false;
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

    if (preFilter_evis_TotalExceptNuclearTarget > 1500 && truth_isSignal )  counter1++;
    if (preFilter_evis_TotalExceptNuclearTarget > 2000 && truth_isSignal )  counter2++;
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
    if (nProtonCandidates == 0) FillHistogram(cutList.hCut_1Track_neutrinoE,CCProtonPi0_neutrino_E * MeV_to_GeV);
    else FillHistogram(cutList.hCut_2Track_neutrinoE,CCProtonPi0_neutrino_E * MeV_to_GeV); 

    if ( applyBeamEnergy && ((CCProtonPi0_neutrino_E * MeV_to_GeV) > max_beamEnergy)) return false;

    cutList.nCut_beamEnergy.increment(truth_isSignal, study1, study2);
    if (nProtonCandidates == 0) cutList.nCut_1Track_beamEnergy.increment(truth_isSignal, study1, study2);
    else cutList.nCut_2Track_beamEnergy.increment(truth_isSignal, study1, study2);

    // ------------------------------------------------------------------------
    // Pi0 Invariant Mass Cut
    // ------------------------------------------------------------------------
    fill_BackgroundSubtractionHists();

    FillHistogram(cutList.hCut_pi0invMass, pi0_invMass);
    FillHistogram(cutList.hCut_pi0invMass_Old, pi0_invMass_Old);
  
    if (nProtonCandidates == 0){
        FillHistogram(cutList.pi0_invMass_1Track, pi0_invMass);
        FillHistogram(cutList.hCut_1Track_pi0invMass,pi0_invMass);
        FillHistogram(cutList.hCut_1Track_pi0invMass_Old,pi0_invMass_Old);
    }else{ 
        FillHistogram(cutList.pi0_invMass_2Track, pi0_invMass);
        FillHistogram(cutList.hCut_2Track_pi0invMass, pi0_invMass);
        FillHistogram(cutList.hCut_2Track_pi0invMass_Old, pi0_invMass_Old); 
    }
    if( pi0_invMass < min_Pi0_invMass || pi0_invMass > max_Pi0_invMass ) return false;
    cutList.nCut_Pi0_invMass.increment(truth_isSignal, study1, study2);
    if (nProtonCandidates == 0) cutList.nCut_1Track_Pi0_invMass.increment(truth_isSignal, study1, study2);
    else cutList.nCut_2Track_Pi0_invMass.increment(truth_isSignal, study1, study2);

    //-------------------------------------------------------------------------
    // Satisfied All Cuts
    //-------------------------------------------------------------------------
    return true;
}

void CCProtonPi0_Analyzer::fill_mc_w()
{
    if(mc_intType == 1) FillHistogram(cutList.mc_w_CCQE, mc_w * MeV_to_GeV);
    if(mc_intType == 2) FillHistogram(cutList.mc_w_RES, mc_w * MeV_to_GeV);
    if(mc_intType == 3) FillHistogram(cutList.mc_w_DIS, mc_w * MeV_to_GeV);
}

void CCProtonPi0_Analyzer::fillInteractionMC()
{
    if(truth_isSignal){
        if(mc_intType == 1) FillHistogram(interaction.final_mc_w_CCQE, mc_w * MeV_to_GeV);
        if(mc_intType == 2) FillHistogram(interaction.final_mc_w_RES, mc_w * MeV_to_GeV);
        if(mc_intType == 3) FillHistogram(interaction.final_mc_w_DIS, mc_w * MeV_to_GeV);

        // Short Proton True Information
        if (nProtonCandidates == 0){
            double proton_mass = 938.27; // MeV
            double proton_true_P = HEP_Functions::calcMomentum(truth_proton_4P[0],truth_proton_4P[1],truth_proton_4P[2]);
            double proton_true_KE = truth_proton_4P[3] - proton_mass;
            FillHistogram(interaction.proton_true_P_1Track, proton_true_P);
            FillHistogram(interaction.proton_true_KE_1Track, proton_true_KE);
        }

        // Ejected Nucleon Count
        double n_nucleons = GetEjectedNucleonCount();
        if (nProtonCandidates == 0){
            FillHistogram(interaction.n_ejected_nucleons_1Track, n_nucleons);
        }else{
            FillHistogram(interaction.n_ejected_nucleons_2Track, n_nucleons);
        }

        // Neutrino Energy: Truth, Error, Difference
        double E_true = mc_incomingPartVec[3] * MeV_to_GeV;
        double E_reco = CCProtonPi0_neutrino_E * MeV_to_GeV;
        double E_reco_1Track_Alt = CCProtonPi0_neutrino_E_1Track_Alt * MeV_to_GeV;
        double E_Error = Data_Functions::getError(E_true, E_reco);
        double E_1Track_Alt_Error = Data_Functions::getError(E_true, E_reco_1Track_Alt);

        // Fill 1Track and 2 Track Enu 
        if(nProtonCandidates == 0){
            FillHistogram(interaction.Enu_True_1Track, E_true);

            FillHistogram(interaction.Enu_1Track_Error, E_Error);
            FillHistogram(interaction.Enu_1Track_Diff, E_reco-E_true);
            FillHistogram(interaction.Enu_1Track_Alt_Error, E_1Track_Alt_Error);
        }else{ 
            FillHistogram(interaction.Enu_True_2Track, E_true);

            FillHistogram(interaction.Enu_2Track_Error, E_Error);
            FillHistogram(interaction.Enu_2Track_Diff, E_reco-E_true);
        }  

        // Response
        FillHistogram(interaction.response_neutrino_E, CCProtonPi0_neutrino_E * MeV_to_GeV, mc_incomingE * MeV_to_GeV);
        FillHistogram(interaction.response_QSq,CCProtonPi0_QSq * MeVSq_to_GeVSq, mc_Q2 * MeVSq_to_GeVSq);
    }
}

void CCProtonPi0_Analyzer::fillInteractionReco()
{
    // Event Kinematics
    if (nProtonCandidates == 0){ 
        FillHistogram(interaction.Enu_1Track, CCProtonPi0_neutrino_E * MeV_to_GeV);
        FillHistogram(interaction.Enu_1Track_Alt, CCProtonPi0_neutrino_E_1Track_Alt * MeV_to_GeV);
    }else{ 
        FillHistogram(interaction.Enu_2Track, CCProtonPi0_neutrino_E * MeV_to_GeV);
    }
    const double Mn = 939.57; //MeV

    double calc_WSq = -CCProtonPi0_QSq + Mn*Mn + 2*Mn*(CCProtonPi0_neutrino_E - muon_E);
    double calc_W = sqrt(calc_WSq);

    FillHistogram(interaction.Enu, CCProtonPi0_neutrino_E * MeV_to_GeV);
    FillHistogram(interaction.QSq, CCProtonPi0_QSq * MeVSq_to_GeVSq);
    FillHistogram(interaction.WSq, CCProtonPi0_WSq * MeVSq_to_GeVSq);
    FillHistogram(interaction.W, std::sqrt(CCProtonPi0_WSq) * MeV_to_GeV);
    FillHistogram(interaction.W_Calc, calc_W * MeV_to_GeV);

    // Vertex & Extra Energy
    if (nProtonCandidates == 0){
        FillHistogram(interaction.vertex_energy_1Track, CCProtonPi0_vertex_energy);
        FillHistogram(interaction.vertex_evis_1Track, vertex_blob_evis );
        FillHistogram(interaction.extra_evis_1Track,Extra_Evis_Leftover);
        FillHistogram(interaction.extra_muon_energy_1Track, Extra_Energy_Muon); 
        FillHistogram(interaction.extra_dispersed_energy_1Track, Extra_Energy_Dispersed); 
        FillHistogram(interaction.extra_rejected_energy_1Track, Extra_Energy_Rejected); 
        FillHistogram(interaction.extra_total_energy_1Track, Extra_Energy_Dispersed + Extra_Energy_Muon + Extra_Energy_Rejected);
    }else{
        FillHistogram(interaction.vertex_energy_2Track, CCProtonPi0_vertex_energy);
        FillHistogram(interaction.vertex_evis_2Track, vertex_blob_evis );
        FillHistogram(interaction.extra_evis_2Track,Extra_Evis_Leftover);
        FillHistogram(interaction.extra_muon_energy_2Track, Extra_Energy_Muon); 
        FillHistogram(interaction.extra_dispersed_energy_2Track, Extra_Energy_Dispersed); 
        FillHistogram(interaction.extra_rejected_energy_2Track, Extra_Energy_Rejected); 
        FillHistogram(interaction.extra_total_energy_2Track, Extra_Energy_Dispersed + Extra_Energy_Muon + Extra_Energy_Rejected);    }

        // Other Event Parameters
        if (nProtonCandidates > 0){
            FillHistogram(interaction.deltaInvMass, calcDeltaInvariantMass() * MeV_to_GeV);
        }

        // Recovered Pi0
        if (is_blobs_recovered){
            FillHistogram(interaction.recovered_Pi0_P, pi0_P * MeV_to_GeV);
            FillHistogram(interaction.recovered_Pi0_theta, pi0_theta * TMath::RadToDeg());
            FillHistogram(interaction.h_recovered_Pi0_P, pi0_P * MeV_to_GeV);
            FillHistogram(interaction.h_recovered_Pi0_theta, pi0_theta * TMath::RadToDeg());
        }else{
            FillHistogram(interaction.h_original_Pi0_P, pi0_P * MeV_to_GeV);
            FillHistogram(interaction.h_original_Pi0_theta, pi0_theta * TMath::RadToDeg());
        }

        // Extra Energy
        FillHistogram(interaction.h_extra_muon_energy, Extra_Energy_Muon); 
        FillHistogram(interaction.h_extra_dispersed_energy, Extra_Energy_Dispersed); 
        FillHistogram(interaction.h_extra_rejected_energy, Extra_Energy_Rejected); 

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
        roundupText<<ev_subrun<<" ^ "<<ev_gate<<" ^ "<<truth_isGamma1_conv_inside<<" ^ "<<truth_isGamma2_conv_inside<<endl;
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
        double true_P = HEP_Functions::calcMomentum(truth_proton_4P[0],truth_proton_4P[1],truth_proton_4P[2]);
        true_P = true_P * MeV_to_GeV;
        double error_P = (reco_P - true_P) / true_P;

        FillHistogram(proton.reco_P_true_P, reco_P,true_P);
        FillHistogram(proton.P_error, error_P);

        // Energy
        double reco_E = proton_E * MeV_to_GeV;
        double true_E = truth_proton_4P[3] * MeV_to_GeV; 
        double error_E = Data_Functions::getError(true_E, reco_E);

        FillHistogram(proton.reco_E_true_E, reco_E,true_E);
        FillHistogram(proton.E_error, error_E);
        FillHistogram(proton.E_Diff, reco_E-true_E);
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
        double pi0_true_P = HEP_Functions::calcMomentum(truth_pi0_4P[0],truth_pi0_4P[1],truth_pi0_4P[2]);
        pi0_true_P = pi0_true_P * MeV_to_GeV; 
        double pi0_P_error = Data_Functions::getError(pi0_true_P, pi0_reco_P);

        FillHistogram(pi0.P_error, pi0_P_error); 
        FillHistogram(pi0.reco_P_true_P, pi0_reco_P, pi0_true_P); 

        // Pi0 Energy
        double reco_E = pi0_E * MeV_to_GeV;
        double true_E = truth_pi0_4P[3] * MeV_to_GeV; 
        double error_E = Data_Functions::getError(true_E, reco_E);

        FillHistogram(pi0.reco_E_true_E, reco_E,true_E);
        FillHistogram(pi0.E_error, error_E);
        FillHistogram(pi0.E_Diff, reco_E-true_E);

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
    if(truth_isSignal){ 
        // Momentum
        double reco_P = muon_P * MeV_to_GeV;
        double true_P = HEP_Functions::calcMomentum(truth_muon_4P[0],truth_muon_4P[1],truth_muon_4P[2]);
        true_P = true_P * MeV_to_GeV;
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

        // Response

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
        FillHistogramWithDefaultErrors(muon.muon_theta_mc_reco_all, muon_theta * TMath::RadToDeg());
        if (truth_isSignal){
            // MC Truth Signal
            FillHistogramWithDefaultErrors(muon.muon_theta_mc_truth_signal, truth_muon_theta * TMath::RadToDeg());
            // MC Reco Signal
            FillHistogramWithDefaultErrors(muon.muon_theta_mc_reco_signal, muon_theta * TMath::RadToDeg());
            // MC Reco vs True -- Response
            FillHistogramWithDefaultErrors(muon.muon_theta_response, muon_theta * TMath::RadToDeg(), truth_muon_theta * TMath::RadToDeg());
        }else{
            // MC Reco Background
            FillHistogramWithDefaultErrors(muon.muon_theta_mc_reco_bckg, muon_theta * TMath::RadToDeg());
        }
    }else{
        // Data
        FillHistogram(muon.muon_theta_all, muon_theta * TMath::RadToDeg());
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

void CCProtonPi0_Analyzer::fillMuonReco()
{
    FillHistogram(muon.E, muon_E * MeV_to_GeV);
    FillHistogram(muon.P, muon_P * MeV_to_GeV);
    FillHistogram(muon.KE, muon_KE * MeV_to_GeV);
    FillHistogram(muon.theta, muon_theta * TMath::RadToDeg());
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

            // Fill Background with Pi0
            int ind = GetBackgroundWithPi0Ind();
            hist[ind]->Fill(var, cvweight);

            // Fill Background Type
            ind = GetBackgroundTypeInd();
            hist[ind]->Fill(var, cvweight);
        }
    }
}

void CCProtonPi0_Analyzer::FillVertErrorBand_MuonTracking(MnvH1D* h, double var)
{
    double correctionErr = normalizer.GetCorrectionErr(CCProtonPi0_minos_trk_p);
    assert(correctionErr < 1.0);
    h->FillVertErrorBand("MuonTracking", var, 1-correctionErr, 1+correctionErr, cvweight);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_MuonTracking(MnvH2D* h, double var1, double var2)
{
    double correctionErr = normalizer.GetCorrectionErr(CCProtonPi0_minos_trk_p);
    assert(correctionErr < 1.0);
    h->FillVertErrorBand("MuonTracking", var1, var2, 1-correctionErr, 1+correctionErr, cvweight);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_Flux(MnvH1D* h, double var)
{
    double enu0 = mc_incomingE/1.e3;
    std::vector<double> random_weights = new_flux::get().get_random_weights(enu0);
    h->FillVertErrorBand("Flux",  var, &random_weights[0],  cvweight, 1.0);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_Flux(MnvH2D* h, double var1, double var2)
{
    double enu0 = mc_incomingE/1.e3;
    std::vector<double> random_weights = new_flux::get().get_random_weights(enu0);
    h->FillVertErrorBand("Flux",  var1, var2,  &random_weights[0],  cvweight, 1.0);
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

int CCProtonPi0_Analyzer::GetBackgroundWithPi0Ind()
{
    // Check For Signal
    if (truth_isSignal){
        cout<<"WARNING! Signal Event requested Background Ind! - Returning -1"<<endl;
        return -1;
    }

    // Background With Pi0
    if (truth_isBckg_NoPi0) return 3;
    else if (truth_isBckg_SinglePi0) return 4;
    else if (truth_isBckg_MultiPi0) return 5;
    else{
        cout<<"WARNING! No Background Type Found - Returning -1"<<endl;
        return -1;
    }
}

int CCProtonPi0_Analyzer::GetBackgroundTypeInd()
{
    // Check For Signal
    if (truth_isSignal){
        cout<<"WARNING! Signal Event requested Background Ind! - Returning -1"<<endl;
        return -1;
    }

    // Background With Pi0
    if (truth_isBckg_NC) return 6;
    else if (truth_isBckg_AntiNeutrino) return 7;
    else if (truth_isBckg_QELike) return 8;
    else if (truth_isBckg_SingleChargedPion) return 9;
    else if (truth_isBckg_SingleChargedPion_ChargeExchanged) return 10;
    else if (truth_isBckg_DoublePionWithPi0) return 11;
    else if (truth_isBckg_DoublePionWithoutPi0) return 12;
    else if (truth_isBckg_MultiPionWithPi0) return 13;
    else if (truth_isBckg_MultiPionWithoutPi0) return 14;
    else if (truth_isBckg_Other) return 15;
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
    // Nothing for now
}

void CCProtonPi0_Analyzer::CorrectEMShowerCalibration()
{
    if (m_isMC){
        double c = 134.98/125;
        pi0_invMass = pi0_invMass*c;
    }
}

void CCProtonPi0_Analyzer::Calc_WeightFromSystematics()
{
    if (m_isMC){
        // Update cvweight with Flux
        double enu0 = mc_incomingE * MeV_to_GeV;
        cvweight = new_flux::get().get_cvweight(enu0);
        
        // Update cvweight with MINOS Efficiency Correction
        const double minos_eff_correction = normalizer.GetCorrection(CCProtonPi0_minos_trk_p);
        cvweight *= minos_eff_correction;
    }else{
        cvweight = 1.0; 
    }
}

void CCProtonPi0_Analyzer::AddErrorBands_Data()
{
    AddVertErrorBands_Data(pi0.pi0_P_all);
    AddVertErrorBands_Data(muon.muon_P_all);
    AddVertErrorBands_Data(muon.muon_theta_all);
}

#endif //CCProtonPi0_Analyzer_cpp


