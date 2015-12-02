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
    nMaxEvents = 100000;

    // Means for EM Energy Correction
    mean_MC_1Track = 118.93;
    mean_MC_2Track = 116.29;
    mean_Data_1Track = 124.27;
    mean_Data_2Track = 125.64;

    // Control Flow
    isDataAnalysis  = true;
    isScanRun = false;
    writeFSParticleMomentum = false;

    applyProtonScore = true;
    pID_KE_Limit = 300.0;
    minProtonScore_LLR = 10.0;
    minPIDDiff = 0.45;

    applyPhotonDistance = true;
    minPhotonDistance = 15; //cm

    applyBeamEnergy = true;
    max_beamEnergy = 20.0; // GeV

    applyUnusedE = true;
    maxUnusedE = 300;

    min_Pi0_invMass = 75.0;
    max_Pi0_invMass = 195.0;

    applyDeltaInvMass = false;
    min_Delta_invMass = 40.0;
    max_Delta_invMass = 200.0;

    latest_ScanID = 0.0;

    counter1 = 0;
    counter2 = 0;

    nTrkrDominate = 0;
    nECALDominate = 0;
    nSCALDominate = 0;
    nTrkrSCAL = 0;
    nTrkrECAL = 0;
    nECALSCAL = 0;
    nOther = 0;
}

void CCProtonPi0_Analyzer::reduce(string playlist)
{
    string rootDir;
    if (m_isMC) rootDir = Folder_List::rootOut + Folder_List::MC + Folder_List::reduced + "ReducedNTuple_v2_40d.root";
    else rootDir = Folder_List::rootOut + Folder_List::Data + Folder_List::reduced + "ReducedNTuple_v2_40d.root";

    cout<<"Reducing NTuple Files to a single file"<<endl;
    cout<<"\tRoot File: "<<rootDir<<endl;
    TFile* f = new TFile(rootDir.c_str(),"CREATE");
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


        // Progress Message on Terminal
        if (jentry%1000000 == 0) cout<<"\tEntry "<<jentry<<endl;

        if (applyMaxEvents && jentry == nMaxEvents){
            cout<<"\tReached Event Limit!"<<endl;
            break;
        }

        // Topology Check - Don't look events without a track
        if (n_prongs < 1) continue;
        
        // Get Cut Statistics
        isPassedAllCuts = getCutStatistics();
        if( !isPassedAllCuts ) continue;
        
        tree->Fill();
    }
    cutList.writeCutTable();
    cutList.writeHistograms();

    cout<<">> Writing "<<rootDir<<endl;
    tree->AutoSave();    
    f->Write();
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
    //getline(DSTFileList,scanFileName);

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

        // Analyze Event or NOT -- Depend on the 1Track or 2 Track Analysis
        if (!AnalyzeEvent() ) continue;

        // Update scanFileName if running for scan
        //if(isScanRun) UpdateScanFileName();
        
        //----------------------------------------------------------------------
        // Fill Background Branches for Background Events
        //----------------------------------------------------------------------
        if(m_isMC && !truth_isSignal) {
            bckgTool.fillBackgroundWithPi0(truth_isBckg_NoPi0, truth_isBckg_SinglePi0, truth_isBckg_MultiPi0, truth_isBckg_withMichel);                                    
            bckgTool.fillBackground(truth_isBckg_NC, truth_isBckg_AntiNeutrino, truth_isBckg_QELike, truth_isBckg_SinglePion, truth_isBckg_DoublePion, truth_isBckg_MultiPion, truth_isBckg_Other, truth_isBckg_withMichel);                                    
        }

        //----------------------------------------------------------------------
        // Data Analysis and Other Studies
        //----------------------------------------------------------------------
        if (isDataAnalysis) fillData();
        if (writeFSParticleMomentum) writeFSParticle4P(jentry);

        // EM Shower Energy Study
        if (truth_isSignal){
            double pi0_match = truth_total_captured_evis_pizero / truth_allClusters_evis_pizero;
            if (pi0_match >= 0.9){ 
                //fillPi0True();
                //fill_EMEnergy_StudyPlots(); 
                //Study_EMShowerTopology();
                Fill_SCAL_nHits();
            }
        }

        

    } // end for-loop

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

    cout<<"nTrkrDominate = "<<nTrkrDominate<<endl;
    cout<<"nECALDominate = "<<nECALDominate<<endl;
    cout<<"nSCALDominate = "<<nSCALDominate<<endl;
    cout<<"nTrkrSCAL = "<<nTrkrSCAL<<endl;
    cout<<"nTrkrECAL = "<<nTrkrECAL<<endl;
    cout<<"nECALSCAL = "<<nECALSCAL<<endl;
    cout<<"nOther = "<<nOther<<endl;


    model1.close();
    model2.close();
    model3.close();


}


//------------------------------------------------------------------------------
//  Constructor
//------------------------------------------------------------------------------
CCProtonPi0_Analyzer::CCProtonPi0_Analyzer(bool isModeReduce, bool isMC, std::string ana_folder) : 
    CCProtonPi0_NTupleAnalysis(),
    interaction(isModeReduce, isMC, ana_folder),
    muon(isModeReduce, isMC, ana_folder),
    proton(isModeReduce, isMC, ana_folder),
    pi0(isModeReduce, isMC, ana_folder),
    pi0Blob(isModeReduce, isMC, ana_folder),
    bckgTool(isModeReduce, ana_folder),
    cutList(isModeReduce, isMC)
{   
    cout<<"Initializing CCProtonPi0_Analyzer"<<endl;

    m_isMC = isMC;
    m_ana_folder = ana_folder;

    specifyRunTime();

    if (isModeReduce){
        cout<<"\tNTuple Reduce Mode -- Will not create Text Files"<<endl;
    }else{
        openTextFiles();
    }

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

bool CCProtonPi0_Analyzer::AnalyzeEvent()
{
    const string ana_folder_all = "All/";
    const string ana_folder_1 = "1Track/";
    const string ana_folder_2 = "2Track/";

    if (m_ana_folder.compare(ana_folder_all) == 0) return true;
    else if (m_ana_folder.compare(ana_folder_1) == 0 ){
        if (nProngs == 1) return true;
        else return false;
    }else if (m_ana_folder.compare(ana_folder_2) == 0 ){
        if (nProngs > 1 ) return true;
        else return false;
    }else{
        cout<<"WARNING! None of the analysis modes matched!"<<endl;
        return false;
    }
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
    // Fill Reconstructed Information
    fillInteractionReco();
    fillMuonReco();
    fillPi0Reco();
    fillPi0BlobReco();
    if(nProngs > 1) fillProtonReco();

    // Fill Truth Information if Exist and Set Errors
    if( m_isMC ){
        fillInteractionTrue();
        fillMuonTrue();
        if(nProngs > 1) fillProtonTrue();
        fillPi0True();
        fillPi0BlobTrue();
    }
}

void CCProtonPi0_Analyzer::fillPi0BlobReco()
{
    // Do Nothing
}

void CCProtonPi0_Analyzer::fillPi0BlobTrue()
{
    // Fill Truth Match Results
    FillEvis_Fractions();
    FillEvis_MostPDG();
    FillEvis_Total(); 
}

void CCProtonPi0_Analyzer::FillEvis_Total()
{
    double captured_fraction;
    if (truth_allClusters_evis_pizero == 0) captured_fraction = -0.5;
    else captured_fraction = truth_total_captured_evis_pizero / truth_allClusters_evis_pizero;

    FillHistogram(pi0Blob.captured_evis_frac_all, captured_fraction);
    if (truth_isSignal) FillHistogram(pi0Blob.captured_evis_frac_signal, captured_fraction);
}

void CCProtonPi0_Analyzer::FillEvis_MostPDG()
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

void CCProtonPi0_Analyzer::FillEvis_Fractions()
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
    FillHistogram(pi0Blob.g1_evis_total_truth, total[i] * HEP_Functions::MeV_to_GeV);
    FillHistogram(pi0Blob.g1_evis_frac_pizero, pizero[i]/total[i]);
    FillHistogram(pi0Blob.g1_evis_frac_piplus, piplus[i]/total[i]);
    FillHistogram(pi0Blob.g1_evis_frac_piminus, piminus[i]/total[i]);
    FillHistogram(pi0Blob.g1_evis_frac_proton, proton[i]/total[i]);
    FillHistogram(pi0Blob.g1_evis_frac_neutron, neutron[i]/total[i]);
    FillHistogram(pi0Blob.g1_evis_frac_muon, muon[i]/total[i]);

    // Gamma 2 
    i = 1;
    FillHistogram(pi0Blob.g2_evis_total_truth, total[i] * HEP_Functions::MeV_to_GeV);
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
            <<CCProtonPi0_muon_px<<", "
            <<CCProtonPi0_muon_py<<", "
            <<CCProtonPi0_muon_pz<<", "
            <<CCProtonPi0_muon_E<<" )"
            <<endl;
        failText<<"Proton 4-P = ( "
            <<CCProtonPi0_proton_px<<", "
            <<CCProtonPi0_proton_py<<", "
            <<CCProtonPi0_proton_pz<<", "
            <<CCProtonPi0_proton_E<<" )"
            <<" Score = "<<CCProtonPi0_proton_LLRScore
            <<endl;
        failText<<"Pi0 4-P = ( "
            <<CCProtonPi0_pi0_px<<", "
            <<CCProtonPi0_pi0_py<<", "
            <<CCProtonPi0_pi0_pz<<", "
            <<CCProtonPi0_pi0_E<<" )"
            <<endl;   
    }
}

void CCProtonPi0_Analyzer::Increment_nCut(vector<CCProtonPi0_Cut> &nCut, bool study1, bool study2)
{
    int ind;

    if (n_prongs == 1) ind = 0;
    else if (n_prongs > 1) ind = 1;
    else cout<<"WARNING No Topology while Incrementing nCut"<<endl;

    nCut[ind].increment(truth_isSignal, study1, study2);
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

    FillHistogram(cutList.hCut_nProngs, n_prongs);
    cutList.nCut_All[0].increment(truth_isSignal, study1, study2);
    cutList.nCut_All[1].increment(truth_isSignal, study1, study2);

    // Vertex Cut -- If Cut_Vertex_None == 1 --> No Event Vertex
    if( Cut_Vertex_None == 1) return false;
    cutList.nCut_Vertex_None[0].increment(truth_isSignal, study1, study2);
    cutList.nCut_Vertex_None[1].increment(truth_isSignal, study1, study2);

    // Vertex Reconstructable Cut
    if( Cut_Vertex_Not_Reconstructable == 1) return false;
    cutList.nCut_Vertex_Not_Reconstructable[0].increment(truth_isSignal, study1, study2);
    cutList.nCut_Vertex_Not_Reconstructable[1].increment(truth_isSignal, study1, study2);

    // Vertex Fiducial Cut
    if( Cut_Vertex_Not_Fiducial == 1) return false;
    cutList.nCut_Vertex_Not_Fiducial[0].increment(truth_isSignal, study1, study2);
    cutList.nCut_Vertex_Not_Fiducial[1].increment(truth_isSignal, study1, study2);

    // Check nVertices
    FillHistogram(cutList.hCut_nVertices, vtx_total_count);

    // Muon Cut -- If Cut_Muon_None == 1 --> No MINOS Matched Muon
    if( Cut_Muon_None == 1) return false;
    Increment_nCut(cutList.nCut_Muon_None, study1, study2);

    // Fill Truth_W for MINOS Matched Signal Events
    if (m_isMC && truth_isSignal) fill_mc_w(); 

    // Anti-Muon Cut
    if( Cut_Muon_Charge == 1) return false;
    Increment_nCut(cutList.nCut_Muon_Charge, study1, study2);

    // Michel Cuts
    if( Cut_Vertex_Michel_Exist == 1 || Cut_EndPoint_Michel_Exist == 1 || Cut_secEndPoint_Michel_Exist == 1 ){
        if (n_prongs == 1) FillHistogram(cutList.hCut_1Track_Michel,1);
        else FillHistogram(cutList.hCut_2Track_Michel,1);
    }else{
        if (n_prongs == 1) FillHistogram(cutList.hCut_1Track_Michel,0);
        else FillHistogram(cutList.hCut_2Track_Michel,0);
    } 
    if( Cut_Vertex_Michel_Exist == 1) return false;
    Increment_nCut(cutList.nCut_Vertex_Michel_Exist, study1, study2);

    if( Cut_EndPoint_Michel_Exist == 1) return false;
    Increment_nCut(cutList.nCut_EndPoint_Michel_Exist, study1, study2);

    if( Cut_secEndPoint_Michel_Exist == 1) return false;
    Increment_nCut(cutList.nCut_secEndPoint_Michel_Exist, study1, study2);

    // Check nTracks After Michel
    FillHistogram(cutList.hCut_nProngs2, n_prongs);
    FillHistogram(cutList.hCut_nTracks, nTracks);
    FillHistogram(cutList.hCut_nTracks2, nTracks_Close + nTracks_Far);
    FillHistogram(cutList.hCut_nTracks_Close, nTracks_Close);
    FillHistogram(cutList.hCut_nTracks_Far, nTracks_Far);
    FillHistogram(cutList.hCut_nTracks_Discarded, nTracks_Discarded);

    // PreFilter Cut
    if(n_prongs == 1){
        FillHistogram(cutList.hCut_1Track_eVis_nuclearTarget,evis_NuclearTarget);
        FillHistogram(cutList.hCut_1Track_eVis_other,evis_TotalExceptNuclearTarget);
    }else{
        FillHistogram(cutList.hCut_2Track_eVis_nuclearTarget,evis_NuclearTarget);
        FillHistogram(cutList.hCut_2Track_eVis_other,evis_TotalExceptNuclearTarget);
    }
    if( Cut_PreFilter_Pi0 == 1) return false;
    Increment_nCut(cutList.nCut_PreFilter_Pi0, study1, study2);

    // ConeBlobs Cut -- If Cut_ConeBlobs == 1 --> Failed Pi0 Reconstruction
    if( Cut_ConeBlobs == 1 || is_houghtransform_applied ) return false;
    Increment_nCut(cutList.nCut_ConeBlobs, study1, study2);

    // Blob Direction Bad Cut
    if ( Cut_BlobDirectionBad == 1 ) return false;
    Increment_nCut(cutList.nCut_BlobDirectionBad, study1, study2);

    // Blobs Bad Cut
    //if ( Cut_BlobsBad == 1 ) return false;
    Increment_nCut(cutList.nCut_BlobsBad, study1, study2);

    // Gamma1 Conv Length Cut Hist
    if (n_prongs == 1) FillHistogram(cutList.hCut_1Track_gamma1ConvDist,CCProtonPi0_gamma1_dist_vtx * 0.1);
    else FillHistogram(cutList.hCut_2Track_gamma1ConvDist,CCProtonPi0_gamma1_dist_vtx * 0.1);

    if (applyPhotonDistance && CCProtonPi0_gamma1_dist_vtx * 0.1 < minPhotonDistance) return false;
    Increment_nCut(cutList.nCut_Photon1DistanceLow, study1, study2);

    // Gamma2 Conv Length Cut
    if (n_prongs == 1) FillHistogram(cutList.hCut_1Track_gamma2ConvDist,CCProtonPi0_gamma2_dist_vtx * 0.1);
    else FillHistogram(cutList.hCut_2Track_gamma2ConvDist,CCProtonPi0_gamma2_dist_vtx * 0.1);

    if (applyPhotonDistance && CCProtonPi0_gamma2_dist_vtx * 0.1 < minPhotonDistance) return false;
    Increment_nCut(cutList.nCut_Photon2DistanceLow, study1, study2);

    // Pi0 Invariant Mass Cut
    if (n_prongs == 1){
        cutList.pi0_invMass_1Track->Fill( ApplyEMEnergyCorrection(CCProtonPi0_pi0_invMass));
        FillHistogram(cutList.hCut_1Track_pi0invMass,ApplyEMEnergyCorrection(CCProtonPi0_pi0_invMass) );
        FillHistogram(cutList.hCut_1Track_pi0invMass_1,ApplyEMEnergyCorrection(CCProtonPi0_pi0_invMass_Old) );
    }else{ 
        cutList.pi0_invMass_2Track->Fill( ApplyEMEnergyCorrection(CCProtonPi0_pi0_invMass) );
        FillHistogram(cutList.hCut_2Track_pi0invMass, ApplyEMEnergyCorrection(CCProtonPi0_pi0_invMass) );
        FillHistogram(cutList.hCut_2Track_pi0invMass_1, ApplyEMEnergyCorrection(CCProtonPi0_pi0_invMass_Old) );
    }

    if( ApplyEMEnergyCorrection(CCProtonPi0_pi0_invMass) < min_Pi0_invMass ||  
            ApplyEMEnergyCorrection(CCProtonPi0_pi0_invMass) > max_Pi0_invMass ) return false;
    Increment_nCut(cutList.nCut_Pi0_invMass, study1, study2);

    // ------------------------------------------------------------------------
    // Proton Candidate Related Cuts 
    // ------------------------------------------------------------------------
    // 1 Track Events Always Satisfies Proton Cuts - No Proton Candidate to apply Cut!
    if (n_prongs == 1){
        Increment_nCut(cutList.nCut_Particle_None, study1, study2);
        Increment_nCut(cutList.nCut_Proton_None, study1, study2);
        Increment_nCut(cutList.nCut_ProtonScore, study1, study2);
        Increment_nCut(cutList.nCut_DeltaInvMass, study1, study2);
    }

    // 2+ Prong Events Must pass the following Cuts
    if ( nProngs > 1 ){
        if( Cut_Particle_None == 1) return false;
        Increment_nCut(cutList.nCut_Particle_None, study1, study2);

        if( Cut_Proton_None == 1) return false;
        Increment_nCut(cutList.nCut_Proton_None, study1, study2);

        // Apply Proton Score to All Proton Candidates
        for( unsigned int i = 0; i < 10 && CCProtonPi0_all_protons_LLRScore[i] != -9.9; i++){
            if ( applyProtonScore ){
                // Use pID Difference for KE < pID_KE_Limit 
                // Use LLR for KE > pID_KE_Limit
                if (CCProtonPi0_all_protons_KE[i] < pID_KE_Limit ){
                    double pIDDiff = CCProtonPi0_all_protons_protonScore[i] - CCProtonPi0_all_protons_pionScore[i];
                    FillHistogram(cutList.hCut_2Track_protonScore_pIDDiff,pIDDiff);
                    if ( pIDDiff < minPIDDiff ) return false;
                }else{
                    FillHistogram(cutList.hCut_2Track_protonScore_LLR,CCProtonPi0_all_protons_LLRScore[i]);
                    if ( CCProtonPi0_all_protons_LLRScore[i] < minProtonScore_LLR ) return false;
                }
            }
        }
        Increment_nCut(cutList.nCut_ProtonScore, study1, study2);

        double delta_invMass = calcDeltaInvariantMass();
        FillHistogram(cutList.hCut_2Track_deltaInvMass,delta_invMass);
        if (applyDeltaInvMass){
            if ( delta_invMass < min_Delta_invMass || delta_invMass > max_Delta_invMass) return false;  
        }
        Increment_nCut(cutList.nCut_DeltaInvMass, study1, study2);
    }

    // Neutrino Energy Cut
    if (n_prongs == 1) FillHistogram(cutList.hCut_1Track_neutrinoE,CCProtonPi0_neutrino_E_Cal * HEP_Functions::MeV_to_GeV);
    else FillHistogram(cutList.hCut_2Track_neutrinoE,CCProtonPi0_neutrino_E_Cal * HEP_Functions::MeV_to_GeV); 

    if ( applyBeamEnergy && ((CCProtonPi0_neutrino_E_Cal * HEP_Functions::MeV_to_GeV) > max_beamEnergy)) return false;
    Increment_nCut(cutList.nCut_beamEnergy, study1, study2);

    // Unused Energy Cut  
    if (n_prongs == 1) FillHistogram(cutList.hCut_1Track_UnusedE,energyUnused_afterReco);
    else FillHistogram(cutList.hCut_2Track_UnusedE,energyUnused_afterReco);

    if ( applyUnusedE && (energyUnused_afterReco > maxUnusedE)) return false;
    Increment_nCut(cutList.nCut_UnusedE, study1, study2);

    //-------------------------------------------------------------------------
    // Satisfied All Cuts
    //-------------------------------------------------------------------------
    return true;
}

void CCProtonPi0_Analyzer::fill_mc_w()
{
    if(mc_intType == 1) cutList.mc_w_CCQE->Fill(mc_w * HEP_Functions::MeV_to_GeV);
    if(mc_intType == 2) cutList.mc_w_RES->Fill(mc_w * HEP_Functions::MeV_to_GeV);
    if(mc_intType == 3) cutList.mc_w_DIS->Fill(mc_w * HEP_Functions::MeV_to_GeV);
}

void CCProtonPi0_Analyzer::fillInteractionTrue()
{
    if(truth_isSignal){
        if(mc_intType == 1) interaction.final_mc_w_CCQE->Fill(mc_w * HEP_Functions::MeV_to_GeV);
        if(mc_intType == 2) interaction.final_mc_w_RES->Fill(mc_w * HEP_Functions::MeV_to_GeV);
        if(mc_intType == 3) interaction.final_mc_w_DIS->Fill(mc_w * HEP_Functions::MeV_to_GeV);
    }

    // True Neutrino Energy
    double E_true = mc_incomingPartVec[3] * HEP_Functions::MeV_to_GeV;
    double E_1Track = CCProtonPi0_neutrino_E_1Track * HEP_Functions::MeV_to_GeV;
    double E_2Track = CCProtonPi0_neutrino_E_2Track * HEP_Functions::MeV_to_GeV;
    double E_Cal = CCProtonPi0_neutrino_E_Cal * HEP_Functions::MeV_to_GeV;
    double E_1Track_Error = Data_Functions::getError(E_true, E_1Track);
    double E_2Track_Error = Data_Functions::getError(E_true, E_2Track);
    double E_Cal_Error = Data_Functions::getError(E_true, E_Cal);

    FillHistogram(interaction.Enu_True, E_true);
    if(nProngs == 1)interaction.Enu_1Track_Error->Fill(E_1Track_Error);
    else interaction.Enu_2Track_Error->Fill(E_2Track_Error);
    interaction.Enu_Cal_Error->Fill(E_Cal_Error);
}

void CCProtonPi0_Analyzer::fillInteractionReco()
{
    // Event Kinematics
    if (nProngs == 1) FillHistogram(interaction.Enu_1Track, CCProtonPi0_neutrino_E_1Track * HEP_Functions::MeV_to_GeV);
    else FillHistogram(interaction.Enu_2Track, CCProtonPi0_neutrino_E_2Track * HEP_Functions::MeV_to_GeV);
    FillHistogram(interaction.Enu_Cal, CCProtonPi0_neutrino_E_Cal * HEP_Functions::MeV_to_GeV);
    FillHistogram(interaction.q2, CCProtonPi0_QSq_Cal * HEP_Functions::MeVSq_to_GeVSq);
    FillHistogram(interaction.w, CCProtonPi0_W_Cal* HEP_Functions::MeV_to_GeV);
    FillHistogram(interaction.wSq, CCProtonPi0_WSq_Cal * HEP_Functions::MeVSq_to_GeVSq);

    // Reconstruction 
    FillHistogram(interaction.E_Unused_afterReco, energyUnused_afterReco);
    FillHistogram(interaction.E_Used_afterReco, energyUsed_afterReco);

    // Other Event Parameters
    if (nProngs >= 2) FillHistogram(interaction.deltaInvMass, calcDeltaInvariantMass() * HEP_Functions::MeV_to_GeV);
    FillHistogram(interaction.nProngs_hist, nProngs);
}

double CCProtonPi0_Analyzer::calcDeltaInvariantMass()
{
    double invMassSq;

    invMassSq = 
        (ApplyEMEnergyCorrection(CCProtonPi0_pi0_E) + CCProtonPi0_proton_E) * 
        (ApplyEMEnergyCorrection(CCProtonPi0_pi0_E) + CCProtonPi0_proton_E) -
        (   (ApplyEMEnergyCorrection(CCProtonPi0_pi0_px) + CCProtonPi0_proton_px) * 
            (ApplyEMEnergyCorrection(CCProtonPi0_pi0_px) + CCProtonPi0_proton_px) + 
            (ApplyEMEnergyCorrection(CCProtonPi0_pi0_py) + CCProtonPi0_proton_py) * 
            (ApplyEMEnergyCorrection(CCProtonPi0_pi0_py) + CCProtonPi0_proton_py) +
            (ApplyEMEnergyCorrection(CCProtonPi0_pi0_pz) + CCProtonPi0_proton_pz) * 
            (ApplyEMEnergyCorrection(CCProtonPi0_pi0_pz) + CCProtonPi0_proton_pz));

    return sqrt(invMassSq);
}

void CCProtonPi0_Analyzer::writeScanFile()
{
    if(isScanRun){
        // Constants for Roundup List
        const string arachne_html = "http://minerva05.fnal.gov/Arachne/arachne.html?det=MV&recoVer=v10r6p13&run=";
        const string entryString  = "&entry=";
        const string other        = "&slice=-1&filetype=dst";
//http://minerva05.fnal.gov/Arachne/arachne.html?det=MV&recoVer=v10r6p13&run=3596&subrun=6&gate=597&slice=7
        roundupText<<arachne_html<<scanFileName<<entryString<<truth_eventID<<other<<" ";
        roundupText<<CCProtonPi0_gamma1_dist_vtx<<"^"<<CCProtonPi0_gamma2_dist_vtx<<"^"<<mc_incomingE<<endl;
    }else{
        cout<<"WARNING! ScanRun is NOT Activated! Are you sure what you are doing?"<<endl;    
    }
}

void CCProtonPi0_Analyzer::closeTextFiles()
{
    logFile.close();

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

    string model1_text = "/minerva/app/users/oaltinok/cmtuser/Minerva_v10r8p7/Ana/CCProtonPi0/NTupleAnalysis/Output/TextFiles/model1.dat"; 
    string model2_text = "/minerva/app/users/oaltinok/cmtuser/Minerva_v10r8p7/Ana/CCProtonPi0/NTupleAnalysis/Output/TextFiles/model2.dat"; 
    string model3_text = "/minerva/app/users/oaltinok/cmtuser/Minerva_v10r8p7/Ana/CCProtonPi0/NTupleAnalysis/Output/TextFiles/model3.dat"; 

    model1.open(model1_text.c_str());
    model2.open(model2_text.c_str());
    model3.open(model3_text.c_str());

    logFileName = Folder_List::output + Folder_List::textOut + m_ana_folder + "LogFile.txt";
    logFile.open(logFileName.c_str());
    if( !logFile.is_open() ){
        cerr<<"Cannot open output text file: "<<logFileName<<endl;
        exit(EXIT_FAILURE);    
    }else{
        cout<<"\t"<<logFileName<<endl;
    }

    // Open Fail-Check File
    failFile = Folder_List::output + Folder_List::textOut + m_ana_folder + "FailChecks.txt";

    failText.open( failFile.c_str() );
    if( !failText.is_open() ){
        cerr<<"Cannot open output text file: "<<failFile<<endl;
        exit(1);
    }else{
        cout<<"\t"<<failFile<<endl;
    }

    if(isScanRun){
        // Open Roundup Text for Arachne Scanning
        string roundupFile = Folder_List::output + Folder_List::textOut + m_ana_folder + "ArachneRoundup.txt";
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

void CCProtonPi0_Analyzer::fillProtonTrue()
{   
    double reco_P = CCProtonPi0_proton_P * HEP_Functions::MeV_to_GeV;
    double true_P = CCProtonPi0_trajProtonProngMomentum[CCProtonPi0_proton_leadingIndice] * HEP_Functions::MeV_to_GeV;
    double error_P = (reco_P - true_P) / true_P;

    proton.reco_P_true_P->Fill(reco_P,true_P);
    proton.P_error->Fill(error_P);
}

void CCProtonPi0_Analyzer::fillProtonReco()
{  
    // Unique Histograms
    FillHistogram(proton.partScore, CCProtonPi0_proton_LLRScore);
    FillHistogram(proton.trackLength, CCProtonPi0_proton_length * HEP_Functions::mm_to_cm);
    FillHistogram(proton.trackKinked, CCProtonPi0_proton_kinked);

    // Standard Histograms
    FillHistogram(proton.E, CCProtonPi0_proton_E * HEP_Functions::MeV_to_GeV);
    FillHistogram(proton.P, CCProtonPi0_proton_P * HEP_Functions::MeV_to_GeV);
    FillHistogram(proton.KE, CCProtonPi0_proton_KE * HEP_Functions::MeV_to_GeV);
    FillHistogram(proton.theta, CCProtonPi0_proton_theta * TMath::RadToDeg());
    FillHistogram(proton.phi, CCProtonPi0_proton_phi * TMath::RadToDeg());
}

void CCProtonPi0_Analyzer::fillPi0TruthMatch()
{
    EvisRatio();
    EvisStacked();
}

void CCProtonPi0_Analyzer::EvisStacked()
{   
    const double min_evis = 10;

    // Gamma 1
    double g1_pi0 = truth_blob1_evis_pizero;
    double g1_pi = truth_blob1_evis_piplus + truth_blob1_evis_piminus;
    double g1_proton = truth_blob1_evis_proton;
    double g1_neutron = truth_blob1_evis_neutron;
    double g1_muon = truth_blob1_evis_muon;

    if (g1_pi0 > min_evis) pi0.g1_evis_pi0->Fill(g1_pi0);
    if (g1_pi > min_evis) pi0.g1_evis_pi->Fill(g1_pi);
    if (g1_proton > min_evis) pi0.g1_evis_proton->Fill(g1_proton);
    if (g1_neutron > min_evis) pi0.g1_evis_neutron->Fill(g1_neutron);
    if (g1_muon > min_evis) pi0.g1_evis_muon->Fill(g1_muon);

    // Gamma 2
    double g2_pi0 = truth_blob2_evis_pizero;
    double g2_pi = truth_blob2_evis_piplus + truth_blob2_evis_piminus;
    double g2_proton = truth_blob2_evis_proton;
    double g2_neutron = truth_blob2_evis_neutron;
    double g2_muon = truth_blob2_evis_muon;

    if (g2_pi0 > min_evis) pi0.g2_evis_pi0->Fill(g2_pi0);
    if (g2_pi > min_evis) pi0.g2_evis_pi->Fill(g2_pi);
    if (g2_proton > min_evis) pi0.g2_evis_proton->Fill(g2_proton);
    if (g2_neutron > min_evis) pi0.g2_evis_neutron->Fill(g2_neutron);
    if (g2_muon > min_evis) pi0.g2_evis_muon->Fill(g2_muon);

    // Gamma 1 + Gamma 2 (Pi0) written as Gamma 3
    double g3_pi0 = g1_pi0 + g2_pi0;
    double g3_pi = g1_pi + g2_pi;
    double g3_proton = g1_proton + g2_proton;
    double g3_neutron = g1_neutron + g2_neutron;
    double g3_muon = g1_muon + g2_muon;

    if (g3_pi0 > min_evis) pi0.g3_evis_pi0->Fill(g3_pi0);
    if (g3_pi > min_evis) pi0.g3_evis_pi->Fill(g3_pi);
    if (g3_proton > min_evis) pi0.g3_evis_proton->Fill(g3_proton);
    if (g3_neutron > min_evis) pi0.g3_evis_neutron->Fill(g3_neutron);
    if (g3_muon > min_evis) pi0.g3_evis_muon->Fill(g3_muon);
}

void CCProtonPi0_Analyzer::EvisRatio()
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

    FillHistogram(pi0.evis_frac_reco_pi0_true_pi0, reco_pi0_true_pi0);
    FillHistogram(pi0.evis_frac_true_pi0_reco_all, true_pi0_reco_all);
    FillHistogram(pi0.evis_frac_reco_pi0_reco_all, reco_pi0_reco_all);
    FillHistogram(pi0.evis_frac_reco_nonpi0_reco_all, 1-reco_pi0_reco_all);

}


void CCProtonPi0_Analyzer::fillPi0True()
{

    //fillPi0TruthMatch();

    // EM Shower Energy Variables
    if ( truth_isSignal ){
        // Gamma 1
        double g1_reco_E = ApplyEMEnergyCorrection(CCProtonPi0_gamma1_E) * HEP_Functions::MeV_to_GeV;
        double g1_reco_E_Old = CCProtonPi0_gamma1_E_Old * HEP_Functions::MeV_to_GeV;
        double g1_true_E = truth_gamma1_4P[3] * HEP_Functions::MeV_to_GeV;
        double g1_E_error = Data_Functions::getError(g1_true_E, g1_reco_E);
        double g1_E_error_Old = Data_Functions::getError(g1_true_E, g1_reco_E_Old);

        pi0.gamma1_true_E->Fill(g1_true_E);
        pi0.gamma1_reco_error_E->Fill(g1_E_error);
        pi0.gamma1_reco_Old_error_E->Fill(g1_E_error_Old);
        pi0.gamma1_reco_E_true_E->Fill(g1_reco_E, g1_true_E);
        pi0.gamma1_true_E_reco_E_error->Fill(g1_true_E,g1_E_error);

        // Gamma 2 
        double g2_reco_E = ApplyEMEnergyCorrection(CCProtonPi0_gamma2_E) * HEP_Functions::MeV_to_GeV;
        double g2_reco_E_Old = CCProtonPi0_gamma2_E_Old * HEP_Functions::MeV_to_GeV;
        double g2_true_E = truth_gamma2_4P[3] * HEP_Functions::MeV_to_GeV;
        double g2_E_error = Data_Functions::getError(g2_true_E, g2_reco_E);
        double g2_E_error_Old = Data_Functions::getError(g2_true_E, g2_reco_E_Old);

        pi0.gamma2_true_E->Fill(g2_true_E);
        pi0.gamma2_reco_error_E->Fill(g2_E_error);
        pi0.gamma2_reco_Old_error_E->Fill(g2_E_error_Old);
        pi0.gamma2_reco_E_true_E->Fill(g2_reco_E, g2_true_E);
        pi0.gamma2_true_E_reco_E_error->Fill(g2_true_E, g2_E_error);
    
        pi0.reco_E_true_E->Fill(g1_reco_E,g1_true_E);
        pi0.reco_E_true_E->Fill(g2_reco_E,g2_true_E);
    }

}

void CCProtonPi0_Analyzer::fillPi0Reco()
{
    // Unique Histograms
    FillHistogram(pi0.invMass, ApplyEMEnergyCorrection(CCProtonPi0_pi0_invMass));
    FillHistogram(pi0.invMass_Old, ApplyEMEnergyCorrection(CCProtonPi0_pi0_invMass_Old));
    
    // Leading Photon - Energetic Photon
    FillHistogram(pi0.gamma1_ConvLength, CCProtonPi0_gamma1_dist_vtx * 0.1);
    FillHistogram(pi0.gamma1_E, ApplyEMEnergyCorrection(CCProtonPi0_gamma1_E) * HEP_Functions::MeV_to_GeV);
    FillHistogram(pi0.gamma1_E_Old, ApplyEMEnergyCorrection(CCProtonPi0_gamma1_E_Old) * HEP_Functions::MeV_to_GeV);
    FillHistogram(pi0.gamma1_theta, CCProtonPi0_gamma1_theta * TMath::RadToDeg());

    // Secondary Photon
    FillHistogram(pi0.gamma2_ConvLength, CCProtonPi0_gamma2_dist_vtx * 0.1);
    FillHistogram(pi0.gamma2_E, ApplyEMEnergyCorrection(CCProtonPi0_gamma2_E) * HEP_Functions::MeV_to_GeV);
    FillHistogram(pi0.gamma2_E_Old, ApplyEMEnergyCorrection(CCProtonPi0_gamma2_E_Old) * HEP_Functions::MeV_to_GeV);
    FillHistogram(pi0.gamma2_theta, CCProtonPi0_gamma2_theta * TMath::RadToDeg());

    double photon_E_asym = abs((CCProtonPi0_gamma1_E - CCProtonPi0_gamma2_E) / (CCProtonPi0_gamma1_E + CCProtonPi0_gamma2_E));  
    FillHistogram(pi0.photonEnergy_Asymmetry, photon_E_asym);

    // Standard Histograms
    FillHistogram(pi0.E, ApplyEMEnergyCorrection(CCProtonPi0_pi0_E) * HEP_Functions::MeV_to_GeV);
    FillHistogram(pi0.P, ApplyEMEnergyCorrection(CCProtonPi0_pi0_P) * HEP_Functions::MeV_to_GeV);
    FillHistogram(pi0.KE, CCProtonPi0_pi0_KE * HEP_Functions::MeV_to_GeV);
    FillHistogram(pi0.theta, CCProtonPi0_pi0_theta * TMath::RadToDeg());
    FillHistogram(pi0.phi, CCProtonPi0_pi0_phi * TMath::RadToDeg());

    // Photon Comparison
    pi0.gamma1_E_gamma2_E->Fill(ApplyEMEnergyCorrection(CCProtonPi0_gamma1_E) * HEP_Functions::MeV_to_GeV, ApplyEMEnergyCorrection(CCProtonPi0_gamma2_E) * HEP_Functions::MeV_to_GeV);
    pi0.gamma1_convLength_gamma2_convLength->Fill(CCProtonPi0_gamma1_dist_vtx * 0.1, CCProtonPi0_gamma2_dist_vtx * 0.1);
}

void CCProtonPi0_Analyzer::fillMuonTrue()
{
    // Do Nothing
}

void CCProtonPi0_Analyzer::fillMuonReco()
{
    FillHistogram(muon.E, CCProtonPi0_muon_E * HEP_Functions::MeV_to_GeV);
    FillHistogram(muon.P, CCProtonPi0_muon_P * HEP_Functions::MeV_to_GeV);
    FillHistogram(muon.KE, CCProtonPi0_muon_KE * HEP_Functions::MeV_to_GeV);
    FillHistogram(muon.theta, CCProtonPi0_muon_theta * TMath::RadToDeg());
    FillHistogram(muon.phi, CCProtonPi0_muon_phi * TMath::RadToDeg());
}

void CCProtonPi0_Analyzer::FillHistogram(vector<MnvH1D*> &hist, double par)
{
    // Always Fill hist[0]
    hist[0]->Fill(par);

    // Fill others only if Analyzing MC
    if (m_isMC){
        // Fill Signal
        if (truth_isSignal){
            hist[1]->Fill(par);
        }else{
            // Fill Background
            hist[2]->Fill(par); // Always Fill ind == 2 -- All Background

            // Fill Background with Pi0
            int ind = GetBackgroundWithPi0Ind();
            hist[ind]->Fill(par);

            // Fill Background Type
            ind = GetBackgroundTypeInd();
            hist[ind]->Fill(par);
        }
    }
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
    else if (truth_isBckg_SinglePion) return 9;
    else if (truth_isBckg_DoublePion) return 10;
    else if (truth_isBckg_MultiPion) return 11;
    else if (truth_isBckg_Other) return 12;
    else{
        cout<<"WARNING! No Background Type Found - Returning -1"<<endl;
        return -1;
    }
}

double CCProtonPi0_Analyzer::ApplyEMEnergyCorrection(double var)
{
    const double Pi0InvMass = 134.98; // MeV
    double correction;

    if (n_prongs == 1){
        if (m_isMC) correction = Pi0InvMass/mean_MC_1Track;
        else correction = Pi0InvMass/mean_Data_1Track;
    }else{ 
        if (m_isMC) correction = Pi0InvMass/mean_MC_2Track;
        else correction = Pi0InvMass/mean_Data_2Track;
    }

    correction = 1;
    return (var * correction);
}


void CCProtonPi0_Analyzer::fill_EMEnergy_StudyPlots()
{
    // Fill Gamma 1
    if (CCProtonPi0_gamma1_evis_hcal < 0.1 && CCProtonPi0_gamma1_evis_scal_X < 0.1){
        // Contained in Tracker
        if (CCProtonPi0_gamma1_evis_ecal < 0.1){
            double evis = CCProtonPi0_gamma1_evis_trkr;
            double true_E = truth_gamma1_4P[3];
            double reco_E = CCProtonPi0_gamma1_E;
            double reco_error = Data_Functions::getError(true_E, reco_E);

            double evis_ratio = true_E/CCProtonPi0_gamma1_evis_trkr;

            // Calculate Energy with New Method
            double calc_E = get_kT(CCProtonPi0_gamma1_evis_trkr) * CCProtonPi0_gamma1_evis_trkr;
            double calc_error = Data_Functions::getError(true_E, calc_E);

            // Fill Histograms
            //pi0.reco_error_trkr->Fill(reco_error);
            //pi0.calc_error_trkr->Fill(calc_error);
            pi0.true_E_evis_trkr_ratio->Fill(true_E,evis_ratio);
            pi0.evis_evis_trkr_ratio->Fill(CCProtonPi0_gamma1_evis_trkr,evis_ratio);
          
            if (evis >= 50 && evis < 100){ 
                pi0.reco_error_trkr->Fill(reco_error);
                pi0.calc_error_trkr->Fill(calc_error);
            }

            if ( evis < 50) pi0.evis_trkr_ratio_1->Fill(evis_ratio);
            else if ( evis < 100) pi0.evis_trkr_ratio_2->Fill(evis_ratio);
            else if ( evis < 150) pi0.evis_trkr_ratio_3->Fill(evis_ratio);
            else if ( evis < 200) pi0.evis_trkr_ratio_4->Fill(evis_ratio);
            else if ( evis < 300) pi0.evis_trkr_ratio_5->Fill(evis_ratio);
            else if ( evis < 700) pi0.evis_trkr_ratio_6->Fill(evis_ratio);
        }

        // Contained in ECAL
        if (CCProtonPi0_gamma1_evis_trkr < 0.1){
            double evis = CCProtonPi0_gamma1_evis_ecal;
            double true_E = truth_gamma1_4P[3];
            double reco_E = CCProtonPi0_gamma1_E;
            double reco_error = Data_Functions::getError(true_E,reco_E);

            double evis_ratio = true_E/CCProtonPi0_gamma1_evis_ecal;

            // Calculate Energy with New Method
            double calc_E = CCProtonPi0_gamma1_evis_ecal * 3.2367;
            double calc_error = Data_Functions::getError(true_E, calc_E);

            // Fill Histograms
            pi0.reco_error_ecal->Fill(reco_error);
            pi0.calc_error_ecal->Fill(calc_error);
            pi0.evis_evis_ecal_ratio->Fill(CCProtonPi0_gamma1_evis_ecal,evis_ratio);
            pi0.gamma_nPlanes->Fill(g1dedx_nplane);
         
            
            if ( evis < 50) pi0.evis_ecal_ratio_1->Fill(evis_ratio);
            else if ( evis < 100) pi0.evis_ecal_ratio_2->Fill(evis_ratio);
            else if ( evis < 150) pi0.evis_ecal_ratio_3->Fill(evis_ratio);
            else if ( evis < 200) pi0.evis_ecal_ratio_4->Fill(evis_ratio);
            else if ( evis < 300) pi0.evis_ecal_ratio_5->Fill(evis_ratio);
        
        }   

        // Partially Contained in Tracker and ECAL
        double true_E = truth_gamma1_4P[3];
        double reco_E = CCProtonPi0_gamma1_E;
        double reco_error = Data_Functions::getError(true_E,reco_E);

        // Calculate Energy with New Method
        double trkr_evis = CCProtonPi0_gamma1_evis_trkr;
        double ecal_evis = CCProtonPi0_gamma1_evis_ecal;
        double calc_E = trkr_evis * get_kT(trkr_evis) + ecal_evis * 3.2367;
        double calc_error = Data_Functions::getError(true_E, calc_E);

        // Fill Histograms
        pi0.reco_error_trkr_ecal->Fill(reco_error);
        pi0.calc_error_trkr_ecal->Fill(calc_error);

    }

    // Fill Gamma 2
    if (CCProtonPi0_gamma2_evis_hcal < 0.1 && CCProtonPi0_gamma2_evis_scal_X < 0.1){
        // Contained in Tracker
        if (CCProtonPi0_gamma2_evis_ecal < 0.1){
            double evis = CCProtonPi0_gamma2_evis_trkr;
            double true_E = truth_gamma2_4P[3];
            double reco_E = CCProtonPi0_gamma2_E;
            double reco_error = Data_Functions::getError(true_E, reco_E);

            double evis_ratio = true_E/CCProtonPi0_gamma2_evis_trkr;

            // Calculate Energy with New Method
            double calc_E = get_kT(CCProtonPi0_gamma2_evis_trkr) * CCProtonPi0_gamma2_evis_trkr;
            double calc_error = Data_Functions::getError(true_E, calc_E);

            // Fill Histograms
            //pi0.reco_error_trkr->Fill(reco_error);
            //pi0.calc_error_trkr->Fill(calc_error);
            pi0.true_E_evis_trkr_ratio->Fill(true_E,evis_ratio);
            pi0.evis_evis_trkr_ratio->Fill(CCProtonPi0_gamma2_evis_trkr,evis_ratio);

            if (evis >= 50 && evis < 100){ 
                pi0.reco_error_trkr->Fill(reco_error);
                pi0.calc_error_trkr->Fill(calc_error);
            }
            if ( evis < 50) pi0.evis_trkr_ratio_1->Fill(evis_ratio);
            else if ( evis < 100) pi0.evis_trkr_ratio_2->Fill(evis_ratio);
            else if ( evis < 150) pi0.evis_trkr_ratio_3->Fill(evis_ratio);
            else if ( evis < 200) pi0.evis_trkr_ratio_4->Fill(evis_ratio);
            else if ( evis < 300) pi0.evis_trkr_ratio_5->Fill(evis_ratio);
            else if ( evis < 700) pi0.evis_trkr_ratio_6->Fill(evis_ratio);
        }

        // Contained in ECAL
        if (CCProtonPi0_gamma2_evis_trkr < 0.1){
            double evis = CCProtonPi0_gamma2_evis_ecal;
            double true_E = truth_gamma2_4P[3];
            double reco_E = CCProtonPi0_gamma2_E;
            double reco_error = Data_Functions::getError(true_E,reco_E);

            double evis_ratio = true_E/CCProtonPi0_gamma2_evis_ecal;

            // Calculate Energy with New Method
            double calc_E = CCProtonPi0_gamma2_evis_ecal * 3.2367;
            double calc_error = Data_Functions::getError(true_E, calc_E);

            // Fill Histograms
            pi0.reco_error_ecal->Fill(reco_error);
            pi0.calc_error_ecal->Fill(calc_error);
            pi0.evis_evis_ecal_ratio->Fill(CCProtonPi0_gamma2_evis_ecal,evis_ratio);
            pi0.gamma_nPlanes->Fill(g2dedx_nplane);
        
            if ( evis < 50) pi0.evis_ecal_ratio_1->Fill(evis_ratio);
            else if ( evis < 100) pi0.evis_ecal_ratio_2->Fill(evis_ratio);
            else if ( evis < 150) pi0.evis_ecal_ratio_3->Fill(evis_ratio);
            else if ( evis < 200) pi0.evis_ecal_ratio_4->Fill(evis_ratio);
            else if ( evis < 300) pi0.evis_ecal_ratio_5->Fill(evis_ratio);
        }       

        // Partially Contained in Tracker and ECAL
        double true_E = truth_gamma2_4P[3];
        double reco_E = CCProtonPi0_gamma2_E;
        double reco_error = Data_Functions::getError(true_E,reco_E);

        // Calculate Energy with New Method
        double trkr_evis = CCProtonPi0_gamma2_evis_trkr;
        double ecal_evis = CCProtonPi0_gamma2_evis_ecal;
        double calc_E = trkr_evis * get_kT(trkr_evis) + ecal_evis * 3.2367;
        double calc_error = Data_Functions::getError(true_E, calc_E);

        // Fill Histograms
        pi0.reco_error_trkr_ecal->Fill(reco_error);
        pi0.calc_error_trkr_ecal->Fill(calc_error);
    
    }

    // Total Energy
    double gamma1_true_E = truth_gamma1_4P[3];
    double gamma1_calc_E = Calc_Gamma_Energy(CCProtonPi0_gamma1_evis_trkr, CCProtonPi0_gamma1_evis_ecal, CCProtonPi0_gamma1_evis_hcal, CCProtonPi0_gamma1_evis_scal_X);
    double g1_error = Data_Functions::getError(gamma1_true_E, gamma1_calc_E);

    double gamma2_true_E = truth_gamma2_4P[3];
    double gamma2_calc_E = Calc_Gamma_Energy(CCProtonPi0_gamma2_evis_trkr, CCProtonPi0_gamma2_evis_ecal, CCProtonPi0_gamma2_evis_hcal, CCProtonPi0_gamma2_evis_scal_X);
    double g2_error = Data_Functions::getError(gamma2_true_E, gamma2_calc_E);

    double mgg_calc = sqrt(2 * gamma1_calc_E * gamma2_calc_E * (1 - CCProtonPi0_pi0_cos_openingAngle) ); 
    // Fill Histograms
    pi0.gamma1_calc_error->Fill(g1_error);
    pi0.gamma2_calc_error->Fill(g2_error);
    pi0.mgg_calc->Fill(mgg_calc);

}

double CCProtonPi0_Analyzer::get_kT(double evis)
{
    const double evis_limit = 697.77; // MeV
    const double kT_default = 1.326;
    double kT;

    if ( evis < evis_limit){
        const double m = -0.0002336;
        const double c = 1.489;
        kT = m*evis + c;
    }else{ 
        kT = kT_default; 
    }

    return kT;
}


double CCProtonPi0_Analyzer::Calc_Gamma_Energy(double evis_trkr, double evis_ecal, double evis_hcal, double evis_scal)
{
    double energy = evis_trkr * get_kT(evis_trkr) + evis_ecal * 3.2367 + evis_hcal * 12.65 + evis_scal * 3.26 * 1.75;

    return energy;
}

void CCProtonPi0_Analyzer::Study_EMShowerTopology()
{
    // ------------------------------------------------------------------------
    // Visible Energy - reco vs true
    // ------------------------------------------------------------------------
    // Fill Gamma1
    pi0.evis_trkr_reco_true->Fill(CCProtonPi0_gamma1_evis_trkr, gamma1_trkr_true_evis);
    pi0.evis_ecal_reco_true->Fill(CCProtonPi0_gamma1_evis_ecal, gamma1_ecal_true_evis);
    pi0.evis_hcal_reco_true->Fill(CCProtonPi0_gamma1_evis_hcal, gamma1_hcal_true_evis);
    pi0.evis_scal_reco_true->Fill(CCProtonPi0_gamma1_evis_scal_X + CCProtonPi0_gamma1_evis_scal_UV, gamma1_scal_true_evis);
   
    // Fill Gamma2
    pi0.evis_trkr_reco_true->Fill(CCProtonPi0_gamma2_evis_trkr, gamma2_trkr_true_evis);
    pi0.evis_ecal_reco_true->Fill(CCProtonPi0_gamma2_evis_ecal, gamma2_ecal_true_evis);
    pi0.evis_hcal_reco_true->Fill(CCProtonPi0_gamma2_evis_hcal, gamma2_hcal_true_evis);
    pi0.evis_scal_reco_true->Fill(CCProtonPi0_gamma2_evis_scal_X + CCProtonPi0_gamma2_evis_scal_UV, gamma2_scal_true_evis);

    // ------------------------------------------------------------------------
    // Energy - reco 1D 
    // ------------------------------------------------------------------------
    // Fill Gamma1
    pi0.energy_trkr->Fill(CCProtonPi0_gamma1_energy_trkr);
    pi0.energy_ecal->Fill(CCProtonPi0_gamma1_energy_ecal);
    pi0.energy_hcal->Fill(CCProtonPi0_gamma1_energy_hcal);
    pi0.energy_scal->Fill(CCProtonPi0_gamma1_energy_scal_X + CCProtonPi0_gamma1_energy_scal_UV);
 
    // Fill Gamma2
    pi0.energy_trkr->Fill(CCProtonPi0_gamma2_energy_trkr);
    pi0.energy_ecal->Fill(CCProtonPi0_gamma2_energy_ecal);
    pi0.energy_hcal->Fill(CCProtonPi0_gamma2_energy_hcal);
    pi0.energy_scal->Fill(CCProtonPi0_gamma2_energy_scal_X + CCProtonPi0_gamma2_energy_scal_UV);
    
    // ------------------------------------------------------------------------
    // Energy - reco vs true 
    // ------------------------------------------------------------------------
    // Fill Gamma1
    double gamma1_reco_trkr = CCProtonPi0_gamma1_energy_trkr; 
    double gamma1_reco_ecal = CCProtonPi0_gamma1_energy_ecal; 
    double gamma1_reco_hcal = CCProtonPi0_gamma1_energy_hcal; 
    double gamma1_reco_scal = CCProtonPi0_gamma1_energy_scal_X + CCProtonPi0_gamma1_energy_scal_UV; 
   
    double gamma1_true_trkr = gamma1_trkr_true_evis * get_kT(gamma1_trkr_true_evis);
    double gamma1_true_ecal = gamma1_ecal_true_evis * 3.2367;
    double gamma1_true_hcal = gamma1_hcal_true_evis * 12.65;
    double gamma1_true_scal = gamma1_scal_true_evis * 3.2367;

    pi0.energy_trkr_reco_true->Fill(gamma1_reco_trkr, gamma1_true_trkr);
    pi0.energy_ecal_reco_true->Fill(gamma1_reco_ecal, gamma1_true_ecal);
    pi0.energy_hcal_reco_true->Fill(gamma1_reco_hcal, gamma1_true_hcal);
    pi0.energy_scal_reco_true->Fill(gamma1_reco_scal, gamma1_true_scal);

    // Fill Gamma2
    double gamma2_reco_trkr = CCProtonPi0_gamma2_energy_trkr; 
    double gamma2_reco_ecal = CCProtonPi0_gamma2_energy_ecal; 
    double gamma2_reco_hcal = CCProtonPi0_gamma2_energy_hcal; 
    double gamma2_reco_scal = CCProtonPi0_gamma2_energy_scal_X + CCProtonPi0_gamma2_energy_scal_UV; 
    
    double gamma2_true_trkr = gamma2_trkr_true_evis * get_kT(gamma2_trkr_true_evis);
    double gamma2_true_ecal = gamma2_ecal_true_evis * 3.2367;
    double gamma2_true_hcal = gamma2_hcal_true_evis * 12.65;
    double gamma2_true_scal = gamma2_scal_true_evis * 3.2367;

    pi0.energy_trkr_reco_true->Fill(gamma2_reco_trkr, gamma2_true_trkr);
    pi0.energy_ecal_reco_true->Fill(gamma2_reco_ecal, gamma2_true_ecal);
    pi0.energy_hcal_reco_true->Fill(gamma2_reco_hcal, gamma2_true_hcal);
    pi0.energy_scal_reco_true->Fill(gamma2_reco_scal, gamma2_true_scal);

    // ------------------------------------------------------------------------
    // Energy Fractions
    // ------------------------------------------------------------------------
    // Fill Gamma1
    double gamma1_trkr_frac = CCProtonPi0_gamma1_energy_trkr/CCProtonPi0_gamma1_E; 
    double gamma1_ecal_frac = CCProtonPi0_gamma1_energy_ecal/CCProtonPi0_gamma1_E; 
    double gamma1_hcal_frac = CCProtonPi0_gamma1_energy_hcal/CCProtonPi0_gamma1_E; 
    double gamma1_scal_frac = gamma1_reco_scal/CCProtonPi0_gamma1_E; 
    
    pi0.energy_frac_trkr->Fill(gamma1_trkr_frac);
    pi0.energy_frac_ecal->Fill(gamma1_ecal_frac);
    pi0.energy_frac_hcal->Fill(gamma1_hcal_frac);
    pi0.energy_frac_scal->Fill(gamma1_scal_frac);

    // Fill Gamma2
    double gamma2_trkr_frac = CCProtonPi0_gamma2_energy_trkr/CCProtonPi0_gamma2_E; 
    double gamma2_ecal_frac = CCProtonPi0_gamma2_energy_ecal/CCProtonPi0_gamma2_E; 
    double gamma2_hcal_frac = CCProtonPi0_gamma2_energy_hcal/CCProtonPi0_gamma2_E; 
    double gamma2_scal_frac = gamma2_reco_scal/CCProtonPi0_gamma2_E; 
    
    pi0.energy_frac_trkr->Fill(gamma2_trkr_frac);
    pi0.energy_frac_ecal->Fill(gamma2_ecal_frac);
    pi0.energy_frac_hcal->Fill(gamma2_hcal_frac);
    pi0.energy_frac_scal->Fill(gamma2_scal_frac);

    // Shower Topology based on Energy Fractions
    Fill_ShowerTopology(gamma1_trkr_frac, gamma1_ecal_frac, gamma1_hcal_frac, gamma1_scal_frac);
    Fill_ShowerTopology(gamma2_trkr_frac, gamma2_ecal_frac, gamma2_hcal_frac, gamma2_scal_frac);

}


void CCProtonPi0_Analyzer::Fill_ShowerTopology(double trkr, double ecal, double hcal, double scal)
{
    double topology = 6;
    double dominate_threshold = 0.80;
    double min_frac = 0.20;

    if (trkr > dominate_threshold){ 
        topology = 0;
        nTrkrDominate++;
    }else if (ecal > dominate_threshold){
        topology = 1;
        nECALDominate++;
    }else if (scal > dominate_threshold){ 
        topology = 2;
        nSCALDominate++;
    }else if (trkr > min_frac && ecal < min_frac && scal > min_frac && hcal < min_frac){
        topology = 3;
        nTrkrSCAL++;
    }else if (trkr > min_frac && ecal > min_frac && scal < min_frac && hcal < min_frac){
        topology = 4;
        nTrkrECAL++;
    }else if (trkr < min_frac && ecal > min_frac && scal > min_frac && hcal < min_frac){ 
        topology = 5;
        nECALSCAL++;
    }else{ 
        topology = 6;
        nOther++;
    }
    pi0.Shower_Topology->Fill(topology);
}

void CCProtonPi0_Analyzer::Fill_SCAL_nHits()
{
    // Fill gamma1
    double gamma1_trkr_nHits_true = gamma1_center_nHits_trkr + gamma1_side_nHits_trkr;
    double gamma1_scal_nHits_true = gamma1_center_nHits_scal + gamma1_side_nHits_scal;

    pi0.trkr_nHits_reco_correct->Fill(gamma1_center_nHits_all,gamma1_center_nHits_trkr); 
    pi0.trkr_nHits_reco_true->Fill(gamma1_center_nHits_all,gamma1_trkr_nHits_true); 
    pi0.scal_nHits_reco_correct->Fill(gamma1_side_nHits_all,gamma1_side_nHits_scal); 
    pi0.scal_nHits_reco_true->Fill(gamma1_side_nHits_all,gamma1_scal_nHits_true); 

    pi0.trkr_nHits_2_reco_true->Fill(gamma1_trkr_nHits_reco,gamma1_trkr_nHits_true);
    pi0.scal_nHits_2_reco_true->Fill(gamma1_scal_nHits_reco,gamma1_scal_nHits_true);
    
    // Fill gamma2
    double gamma2_trkr_nHits_true = gamma2_center_nHits_trkr + gamma2_side_nHits_trkr;
    double gamma2_scal_nHits_true = gamma2_center_nHits_scal + gamma2_side_nHits_scal;

    pi0.trkr_nHits_reco_correct->Fill(gamma2_center_nHits_all,gamma2_center_nHits_trkr); 
    pi0.trkr_nHits_reco_true->Fill(gamma2_center_nHits_all,gamma2_trkr_nHits_true); 
    pi0.scal_nHits_reco_correct->Fill(gamma2_side_nHits_all,gamma2_side_nHits_scal); 
    pi0.scal_nHits_reco_true->Fill(gamma2_side_nHits_all,gamma2_scal_nHits_true); 

    pi0.trkr_nHits_2_reco_true->Fill(gamma2_trkr_nHits_reco,gamma2_trkr_nHits_true);
    pi0.scal_nHits_2_reco_true->Fill(gamma2_scal_nHits_reco,gamma2_scal_nHits_true);
   
    // model1 = Old Model    
    model1<<gamma1_trkr_nHits_true<<" "<<gamma1_center_nHits_all<<" "<<gamma1_scal_nHits_true<<" "<<gamma1_side_nHits_all<<endl;
    model1<<gamma2_trkr_nHits_true<<" "<<gamma2_center_nHits_all<<" "<<gamma2_scal_nHits_true<<" "<<gamma2_side_nHits_all<<endl;

    // model2 = Geometrical Model
    model2<<gamma1_trkr_nHits_true<<" "<<gamma1_trkr_nHits_reco<<" "<<gamma1_scal_nHits_true<<" "<<gamma1_scal_nHits_reco<<endl;
    model2<<gamma2_trkr_nHits_true<<" "<<gamma2_trkr_nHits_reco<<" "<<gamma2_scal_nHits_true<<" "<<gamma2_scal_nHits_reco<<endl;

    // model3 = Improved Geometrical Model
    model3<<gamma1_improved_trkr_nHits_true<<" "<<gamma1_improved_trkr_nHits_reco<<" "<<gamma1_improved_scal_nHits_true<<" "<<gamma1_improved_scal_nHits_reco<<endl;
    model3<<gamma2_improved_trkr_nHits_true<<" "<<gamma2_improved_trkr_nHits_reco<<" "<<gamma2_improved_scal_nHits_true<<" "<<gamma2_improved_scal_nHits_reco<<endl;


    // Fill Error Histograms
    double trkr_nHits_true;
    double scal_nHits_true;

    if (gamma1_trkr_nHits_true == 0) trkr_nHits_true = 0.01;
    else trkr_nHits_true = gamma1_trkr_nHits_true;

    double gamma1_old_model_trkr_error = Data_Functions::getError(trkr_nHits_true,gamma1_center_nHits_all);
    double gamma1_new_model_trkr_error = Data_Functions::getError(trkr_nHits_true,gamma1_trkr_nHits_reco);
    double gamma1_improved_model_trkr_error = Data_Functions::getError(trkr_nHits_true,gamma1_improved_trkr_nHits_reco);
    pi0.old_model_trkr_nHits_error->Fill(gamma1_old_model_trkr_error);
    pi0.new_model_trkr_nHits_error->Fill(gamma1_new_model_trkr_error);
    pi0.improved_model_trkr_nHits_error->Fill(gamma1_improved_model_trkr_error);
 
    if (gamma1_scal_nHits_true == 0) scal_nHits_true = 0.01;
    else scal_nHits_true = gamma1_scal_nHits_true;

    double gamma1_old_model_scal_error = Data_Functions::getError(scal_nHits_true,gamma1_center_nHits_all);
    double gamma1_new_model_scal_error = Data_Functions::getError(scal_nHits_true,gamma1_scal_nHits_reco);
    double gamma1_improved_model_scal_error = Data_Functions::getError(scal_nHits_true,gamma1_improved_scal_nHits_reco);
    pi0.old_model_scal_nHits_error->Fill(gamma1_old_model_scal_error);
    pi0.new_model_scal_nHits_error->Fill(gamma1_new_model_scal_error);
    pi0.improved_model_scal_nHits_error->Fill(gamma1_improved_model_scal_error);


    // gamma2
    if (gamma2_trkr_nHits_true == 0) trkr_nHits_true = 0.01;
    else trkr_nHits_true = gamma2_trkr_nHits_true;

    double gamma2_old_model_trkr_error = Data_Functions::getError(trkr_nHits_true,gamma2_center_nHits_all);
    double gamma2_new_model_trkr_error = Data_Functions::getError(trkr_nHits_true,gamma2_trkr_nHits_reco);
    double gamma2_improved_model_trkr_error = Data_Functions::getError(trkr_nHits_true,gamma2_improved_trkr_nHits_reco);
    pi0.old_model_trkr_nHits_error->Fill(gamma2_old_model_trkr_error);
    pi0.new_model_trkr_nHits_error->Fill(gamma2_new_model_trkr_error);
    pi0.improved_model_trkr_nHits_error->Fill(gamma2_improved_model_trkr_error);
 
    if (gamma2_scal_nHits_true == 0) scal_nHits_true = 0.01;
    else scal_nHits_true = gamma2_scal_nHits_true;

    double gamma2_old_model_scal_error = Data_Functions::getError(scal_nHits_true,gamma2_center_nHits_all);
    double gamma2_new_model_scal_error = Data_Functions::getError(scal_nHits_true,gamma2_scal_nHits_reco);
    double gamma2_improved_model_scal_error = Data_Functions::getError(scal_nHits_true,gamma2_improved_scal_nHits_reco);
    pi0.old_model_scal_nHits_error->Fill(gamma2_old_model_scal_error);
    pi0.new_model_scal_nHits_error->Fill(gamma2_new_model_scal_error);
    pi0.improved_model_scal_nHits_error->Fill(gamma2_improved_model_scal_error);


    // Fill minZ Information
    
    if ( std::abs(gamma1_new_model_trkr_error) < 0.33 ){
        pi0.scal_minZ_evis->Fill(gamma1_scal_minZ_evis);
        pi0.scal_minZ_nDigits->Fill(gamma1_scal_minZ_nDigits);
    }
    if ( std::abs(gamma2_new_model_trkr_error) < 0.33){
        pi0.scal_minZ_evis->Fill(gamma2_scal_minZ_evis);
        pi0.scal_minZ_nDigits->Fill(gamma2_scal_minZ_nDigits);
    }
}

#endif //CCProtonPi0_Analyzer_cpp

