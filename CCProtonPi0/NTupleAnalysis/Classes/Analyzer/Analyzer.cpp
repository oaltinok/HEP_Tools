/*
    See Analyzer.h header or Class Information
*/

#ifndef Analyzer_cpp
#define Analyzer_cpp

#include "Analyzer.h"

using namespace std;

void Analyzer::specifyRunTime()
{
    applyMaxEvents = false;
    nMaxEvents = 100000;
    
    // Control Flow
    analyze_NoProtonEvents = true;
    isDataAnalysis  = true;
    is_pID_Studies  = true;
    isMichelStudy = false;
    writeFSParticleMomentum = false;
    
    applyProtonScore = true;
    pID_KE_Limit = 300.0;
    minProtonScore_LLR = 10.0;
    minPIDDiff = 0.45;
    
    applyPhotonDistance = true;
    minPhotonDistance = 15; //cm
    
    applyBeamEnergy = true;
    max_beamEnergy = 10.0; // GeV
    
    applyQSq = false;
    max_QSq = 1.5;
    
    applyUnusedE = true;
    maxUnusedE = 300;
    
    min_Pi0_invMass = 40.0;
    max_Pi0_invMass = 200.0;
    
    applyDeltaInvMass = false;
    min_Delta_invMass = 40.0;
    max_Delta_invMass = 200.0;
    
    SENTINEL = -9.9;
    
    latest_ScanID = 0.0;

    //Select Branches to Activate
    m_ActivateMC            = true;
    m_ActivateInteraction   = true;
    m_ActivatePi0           = true;
 
    debug_Counter = 0.0;
    
    N_trueMichel_before = 0.0;
    N_trueMichel_after = 0.0;
    N_trueMichel_afterAll = 0.0;
    N_noMichel_before = 0.0;
    N_noMichel_after = 0.0;
    
    max_nFSPart = 10;

    hasParticleTruthInfo = false;
    
    // Default Beam Configuration
    beam_p3.SetXYZ(1.0,1.0,1.0);
    beam_p3.SetPhi(-1.554);
    beam_p3.SetTheta(0.059);
}

void Analyzer::run(string playlist)
{
    if (!isAnalysisModeSelected){
        cout<<"Warning No Analysis Mode Selected!"<<endl;
        cout<<"Returning!..."<<endl;
        return;
    }
    
    //------------------------------------------------------------------------
    // Create chain
    //------------------------------------------------------------------------
    TChain* fChain = new TChain("CCProtonPi0");
    Init(playlist, fChain);

    if (!fChain) return;
    if (fChain == 0) return;
      
    //------------------------------------------------------------------------
    // Branch Selection for Performance
    //------------------------------------------------------------------------
    fChain->SetBranchStatus("*",0);  // disable all branches
    
    fChain->SetBranchStatus("reco_eventID",1);
    fChain->SetBranchStatus("Cut_*",1);  // Cut List Activated by Default
    fChain->SetBranchStatus("ev_*",1);
    
    if(m_ActivateMC){
        fChain->SetBranchStatus("truth_*",1);
        fChain->SetBranchStatus("mc_*",1);
        fChain->SetBranchStatus("CCProtonPi0*",1);
    }
    
    if(m_ActivateInteraction){
        fChain->SetBranchStatus("preFilter_*",1);
        fChain->SetBranchStatus("nProngs",1);
        fChain->SetBranchStatus("vtx*",1);
        fChain->SetBranchStatus("evis_*",1);
        fChain->SetBranchStatus("energyUnused*",1);
        fChain->SetBranchStatus("energyUsed*",1);
        fChain->SetBranchStatus("CCProtonPi0_total_*",1);
    }
    
    if(m_ActivatePi0){
        fChain->SetBranchStatus("pi0_*",1);
        fChain->SetBranchStatus("gamma*",1);
        fChain->SetBranchStatus("g1blob*",1);
        fChain->SetBranchStatus("g2blob*",1);
    }
    
        
    //------------------------------------------------------------------------
    // Loop over Chain
    //------------------------------------------------------------------------
    Long64_t nbytes = 0, nb = 0;
    
    cout<<"Looping over all entries"<<endl;
    
    Long64_t nentries = fChain->GetEntriesFast();
    
    std::string testOutput;
    
    // Get First Line for the first File
    getline(DSTFileList,scanFileName);

    for (Long64_t jentry=0; jentry<nentries;jentry++) {

        nb = fChain->GetEntry(jentry);   nbytes += nb;    
        Long64_t ientry = fChain->GetEntry(jentry);
        if (ientry == 0) {
            cout<<"\tGetEntry failure "<<jentry<<endl;
            break;
        }
    
        // Progress Message on Terminal
        if (jentry%25000 == 0){
            cout<<"\tEntry "<<jentry<<endl;
        }
        //----------------------------------------------------------------------
        // Seek & Destroy Problematic Events
        
//         if (jentry == 318) continue;
//         if (jentry == 1721) continue;
        
        //----------------------------------------------------------------------
        
//         cout<<jentry<<endl;
        
        /*
            EventID's are increasing in a file.
            If event ID is less than the latest ID then the file is changed
        */
        if(latest_ScanID >= truth_eventID){
            getline(DSTFileList,scanFileName);;
        }
        
        latest_ScanID = truth_eventID;
       
        if (applyMaxEvents && jentry == nMaxEvents){
            cout<<"\tReached Event Limit!"<<endl;
            break;
        }
        
        // Decide to Analyze Event or NOT
        if( !analyzeEvent() ) continue;
        
        
        //----------------------------------------------------------------------
        // Variable Test Site -- Before Cuts
        //----------------------------------------------------------------------

        //----------------------------------------------------------------------
        // Get Cut Statistics
        //----------------------------------------------------------------------
        isPassedAllCuts = getCutStatistics();
        if( !isPassedAllCuts ) continue;
        
        //----------------------------------------------------------------------
        // Variable Test Site -- After Cuts
        //----------------------------------------------------------------------
/*        
        if ( nProngs == 2){
                
            cout<<CCProtonPi0_proton_ekin[indRecoProton]<<" "<<CCProtonPi0_proton_E[indRecoProton]-938.27<<endl;       
        }
        
        
        */

        //----------------------------------------------------------------------
        // Fill Background Branches
        //----------------------------------------------------------------------
        fillBackgroundBranches();

        //----------------------------------------------------------------------
        // Data Analysis and Other Studies
        //----------------------------------------------------------------------
        if (isDataAnalysis) fillData();
        if (is_pID_Studies && nProngs == 2){ 
            pIDTool.FillHistograms(CCProtonPi0_protonScore_LLR[indRecoProton],CCProtonPi0_protonScore[indRecoProton],CCProtonPi0_pionScore[indRecoProton],
                                    CCProtonPi0_trajProtonProngPDG[indRecoProton],CCProtonPi0_proton_E[indRecoProton] );
        }
        if (writeFSParticleMomentum) writeFSParticle4P(jentry);
        
        //----------------------------------------------------------------------
        // Fail Checks
        //----------------------------------------------------------------------

    } // end for-loop
    
    cout<<"debug_Counter = "<<debug_Counter<<endl;
    cout<<"N_trueMichel_before = "<<N_trueMichel_before<<endl;
    cout<<"N_trueMichel_after  = "<<N_trueMichel_after <<endl;
    cout<<"N_trueMichel_afterAll = "<<N_trueMichel_afterAll<<endl;
    cout<<"N_noMichel_before = "<<N_noMichel_before<<endl;
    cout<<"N_noMichel_after  = "<<N_noMichel_after <<endl;
    
    if(is_pID_Studies) pIDTool.get_pID_Stats();
    
    //--------------------------------------------------------------------------
    // Write Text Files
    //--------------------------------------------------------------------------
    cutList.writeCutTable();
    writeBackgroundTable();
    getPi0Family(); // Pi0 Family Information written inside FailFile
    
    //--------------------------------------------------------------------------
    // Write Root Files
    //--------------------------------------------------------------------------
    write_RootFile();           //CCProtonPi0
    muon.write_RootFile();
    proton.write_RootFile();
    pion.write_RootFile();
    if(is_pID_Studies) pIDTool.write_RootFile();
    
    
    closeTextFiles();
    
}

//------------------------------------------------------------------------------
//  Constructor
//------------------------------------------------------------------------------
Analyzer::Analyzer(int nMode) : 
    muon(nMode),
    proton(nMode),
    pion(nMode),
    pIDTool(nMode),
    cutList(nMode)
{
    setAnalysisMode(nMode);
    
    initInteraction();
    
    openTextFiles();
    
    initBackgroundStudy();

    specifyRunTime();
    
    cout<<"Initialization Finished!\n"<<endl;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
//
//     Specific Functions
//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

void Analyzer::fillBackgroundBranches()
{
    if (!truth_isSignal){
        if (nProngs == 1){
            if(truth_isBckg_QELike){ 
                setBackgroundBranch(    bckg_1Prong_QELike,
                                        truth_isBckg_withAntiMuon,
                                        truth_isBckg_withMichel,
                                        truth_isBckg_withPrimaryPi0,
                                        truth_isBckg_withSecondaryPi0);
            }
            
            if(truth_isBckg_SinglePiPlus){ 
                setBackgroundBranch(    bckg_1Prong_SinglePiPlus,
                                        truth_isBckg_withAntiMuon,
                                        truth_isBckg_withMichel,
                                        truth_isBckg_withPrimaryPi0,
                                        truth_isBckg_withSecondaryPi0);
            }
            
            if(truth_isBckg_SinglePiMinus){ 
                setBackgroundBranch(    bckg_1Prong_SinglePiMinus,
                                        truth_isBckg_withAntiMuon,
                                        truth_isBckg_withMichel,
                                        truth_isBckg_withPrimaryPi0,
                                        truth_isBckg_withSecondaryPi0);
            }
            
            if(truth_isBckg_MultiPion){ 
                setBackgroundBranch(    bckg_1Prong_MultiPion,
                                        truth_isBckg_withAntiMuon,
                                        truth_isBckg_withMichel,
                                        truth_isBckg_withPrimaryPi0,
                                        truth_isBckg_withSecondaryPi0);
            }
            
            if(truth_isBckg_MultiPiZero){ 
                setBackgroundBranch(    bckg_1Prong_MultiPiZero,
                                        truth_isBckg_withAntiMuon,
                                        truth_isBckg_withMichel,
                                        truth_isBckg_withPrimaryPi0,
                                        truth_isBckg_withSecondaryPi0);
            }
            
            if(truth_isBckg_Other){ 
                setBackgroundBranch(    bckg_1Prong_Other,
                                        truth_isBckg_withAntiMuon,
                                        truth_isBckg_withMichel,
                                        truth_isBckg_withPrimaryPi0,
                                        truth_isBckg_withSecondaryPi0);
        }
        }else if (nProngs == 2){
                if(truth_isBckg_QELike){ 
                setBackgroundBranch(    bckg_2Prong_QELike,
                                        truth_isBckg_withAntiMuon,
                                        truth_isBckg_withMichel,
                                        truth_isBckg_withPrimaryPi0,
                                        truth_isBckg_withSecondaryPi0);
            }
            
            if(truth_isBckg_SinglePiPlus){ 
                setBackgroundBranch(    bckg_2Prong_SinglePiPlus,
                                        truth_isBckg_withAntiMuon,
                                        truth_isBckg_withMichel,
                                        truth_isBckg_withPrimaryPi0,
                                        truth_isBckg_withSecondaryPi0);
            }
            
            if(truth_isBckg_SinglePiMinus){ 
                setBackgroundBranch(    bckg_2Prong_SinglePiMinus,
                                        truth_isBckg_withAntiMuon,
                                        truth_isBckg_withMichel,
                                        truth_isBckg_withPrimaryPi0,
                                        truth_isBckg_withSecondaryPi0);
            }
            
            if(truth_isBckg_MultiPion){ 
                setBackgroundBranch(    bckg_2Prong_MultiPion,
                                        truth_isBckg_withAntiMuon,
                                        truth_isBckg_withMichel,
                                        truth_isBckg_withPrimaryPi0,
                                        truth_isBckg_withSecondaryPi0);
            }
            
            if(truth_isBckg_MultiPiZero){ 
                setBackgroundBranch(    bckg_2Prong_MultiPiZero,
                                        truth_isBckg_withAntiMuon,
                                        truth_isBckg_withMichel,
                                        truth_isBckg_withPrimaryPi0,
                                        truth_isBckg_withSecondaryPi0);
            }
            
            if(truth_isBckg_Other){ 
                setBackgroundBranch(    bckg_2Prong_Other,
                                        truth_isBckg_withAntiMuon,
                                        truth_isBckg_withMichel,
                                        truth_isBckg_withPrimaryPi0,
                                        truth_isBckg_withSecondaryPi0);
            }
        }
    } 
}

void Analyzer::setBackgroundBranch(vector<double>& background, bool hasAntiMuon, bool hasMichel, bool hasPrimaryPi0, bool hasSecondaryPi0 )
{
    bool hasPi0Candidate = hasPrimaryPi0 || hasSecondaryPi0;
    background[0]++;                        // Total
    if(hasAntiMuon) background[1]++;        // AntiMuon
    if(hasMichel) background[2]++;          // Michel
    if(hasPrimaryPi0) background[3]++;      // Primary Pi0
    if(hasSecondaryPi0) background[4]++;    // Secondary Pi0
    if(hasPi0Candidate) background[5]++;    // Pi0 Candidate
    
}

void Analyzer::getPi0Family()
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
//             failText[t]<<"Mother = "<<PDG_pi0_Mother[i]<<endl;
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
//             failText[t]<<"GrandMother = "<<PDG_pi0_GrandMother[i]<<endl;
        }
    }
    
    for(int t = 0; t < nTopologies; t++){
        cout<<">> Writing "<<failFile[t]<<endl;
        
        failText[t]<<std::left;
        failText[t]<<"-----------------------------------------------------------------"<<endl;
        failText[t].width(20); failText[t]<<"Mother"<<endl;
        failText[t].width(20); failText[t]<<"DIS"<<" = "<<nMother_DIS<<endl;
        failText[t].width(20); failText[t]<<"Delta_p_1232"<<" = "<<nMother_Delta_p_1232<<endl;
        failText[t].width(20); failText[t]<<"Delta_p_1620"<<" = "<<nMother_Delta_p_1620<<endl;
        failText[t].width(20); failText[t]<<"Delta_p_1700"<<" = "<<nMother_Delta_p_1700<<endl;
        failText[t].width(20); failText[t]<<"N_p_1440"<<" = "<<nMother_N_p_1440<<endl;
        failText[t].width(20); failText[t]<<"N_p_1520"<<" = "<<nMother_N_p_1520<<endl;
        failText[t].width(20); failText[t]<<"N_p_1535"<<" = "<<nMother_N_p_1535<<endl;
        failText[t].width(20); failText[t]<<"N_p_1650"<<" = "<<nMother_N_p_1650<<endl;
        failText[t].width(20); failText[t]<<"N_p_1675"<<" = "<<nMother_N_p_1675<<endl;
        failText[t].width(20); failText[t]<<"N_p_1680"<<" = "<<nMother_N_p_1680<<endl;
        failText[t].width(20); failText[t]<<"pi_plus"<<" = "<<nMother_pi_plus<<endl;
        failText[t].width(20); failText[t]<<"pi_minus"<<" = "<<nMother_pi_minus<<endl;
        failText[t].width(20); failText[t]<<"No PDG"<<" = "<<nMother_NoPDG<<endl;
        failText[t].width(20); failText[t]<<"Other"<<" = "<<nMother_Other<<endl;
        
        failText[t]<<endl;
        failText[t].width(20); failText[t]<<"GrandMother"<<endl;
        failText[t].width(20); failText[t]<<"neutron"<<" = "<<nGrandMother_neutron<<endl;
        failText[t].width(20); failText[t]<<"proton"<<" = "<<nGrandMother_proton<<endl;
        failText[t].width(20); failText[t]<<"Delta_pp_1232"<<" = "<<nGrandMother_Delta_pp_1232<<endl;
        failText[t].width(20); failText[t]<<"Delta_pp_1620"<<" = "<<nGrandMother_Delta_pp_1620<<endl;
        failText[t].width(20); failText[t]<<"Delta_pp_1700"<<" = "<<nGrandMother_Delta_pp_1700<<endl;
        failText[t].width(20); failText[t]<<"No PDG"<<" = "<<nGrandMother_NoPDG<<endl;
        failText[t].width(20); failText[t]<<"Other"<<" = "<<nGrandMother_Other<<endl;
        failText[t]<<"-----------------------------------------------------------------"<<endl;
    }
}


void Analyzer::fillData()
{          
    // Fill Reconstructed Information
    fillInteractionReco();
    fillMuonReco();
    fillPionReco();
    if(nProngs == 2) fillProtonReco();
    
    fillInteractionTrue();
    // Fill Truth Information if Exist and Set Errors
    if( hasParticleTruthInfo ){
        findTrueProton();

        fillMuonTrue();
        fillProtonTrue();
        fillPionTrue();
        
        muon.set_errors();
        proton.set_errors();
        pion.set_errors();
    }
    
    // Fill Histograms
    fillHistograms();            
   
}

//------------------------------------------------------------------------------
// analyzeEvent()
//      Decides to Analyze or Skip the Event using Analysis Mode
//------------------------------------------------------------------------------
bool Analyzer::analyzeEvent()
{
    bool isAnalyzable;
    
    if ( anaMode == 1 ){
        if( truth_isSignal) isAnalyzable = true;
        else isAnalyzable = false;
        hasParticleTruthInfo = true;
    }else if ( anaMode == 2 ){
        if( truth_isSignal) isAnalyzable = false;
        else isAnalyzable = true;
    }else{
        isAnalyzable = true;
    }
    
    return isAnalyzable;
    
}

//------------------------------------------------------------------------------
// Analysis Mode
//      1) Signal Events
//      2) Background Events
//      Other) All Events
//------------------------------------------------------------------------------
void Analyzer::setAnalysisMode(int nMode)
{  
    if ( nMode == 1) {
        cout<<"----------------------------------------------------------------------"<<endl;
        cout<<"Analysis Mode: Signal - Only Signal Events will be Analyzed"<<endl;
        branchDir = Folder_List::signal;
    }else if ( nMode == 2){
        cout<<"----------------------------------------------------------------------"<<endl;
        cout<<"Analysis Mode: Background - Only Background Events will be Analyzed"<<endl;
        branchDir = Folder_List::background;
    }else{
        cout<<"----------------------------------------------------------------------"<<endl;
        cout<<"Analysis Mode: All - All Events will be Analyzed"<<endl;
        branchDir = Folder_List::allEvents;
    }
    
    cout<<"----------------------------------------------------------------------"<<endl;
    
    anaMode = nMode;
    isAnalysisModeSelected = true;
}

void Analyzer::writeFSParticle4P(Long64_t nEntry)
{
    for (int t = 0; t < nTopologies; t++){
        // Particle NTuple Info after All Cuts
        failText[t]<<"----------------------------------------------------------------------"<<endl;
        failText[t]<<nEntry<<endl;
        failText[t]<<"Muon 4-P = ( "
                <<CCProtonPi0_muon_px<<", "
                <<CCProtonPi0_muon_py<<", "
                <<CCProtonPi0_muon_pz<<", "
                <<CCProtonPi0_muon_E<<" )"
                <<endl;
        failText[t]<<"Proton 4-P = ( "
                <<CCProtonPi0_proton_px[indRecoProton]<<", "
                <<CCProtonPi0_proton_py[indRecoProton]<<", "
                <<CCProtonPi0_proton_pz[indRecoProton]<<", "
                <<CCProtonPi0_proton_E[indRecoProton]<<" )"
                <<" Score = "<<CCProtonPi0_protonScore_LLR[indRecoProton]
                <<endl;
        failText[t]<<"Pi0 4-P = ( "
                <<pi0_px<<", "
                <<pi0_py<<", "
                <<pi0_pz<<", "
                <<pi0_E<<" )"
                <<endl;   
    }
}

bool Analyzer::getCutStatistics()
{
    /*
     Selection Studies 
        Assign the selection to parameters study1 and study2
        Cut Objects will count them and print them in the Cut Table
    */
    bool michelFound_Vertex = false;
    bool michelFound_Track = false;
    bool michelFound_Track2 = false;
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

    // Increment Cut Number for Both Topologies at the Same Time Until nProngs Cut   
    cutList.nCut_All.increment(truth_isSignal, study1, study2); // Count All Events before Cuts

    if( Cut_Vertex_None == 1) return false;
    cutList.nCut_Vertex_None.increment(truth_isSignal, study1, study2);
        
    if( Cut_Vertex_Not_Reconstructable == 1) return false;
    cutList.nCut_Vertex_Not_Reconstructable.increment(truth_isSignal, study1, study2);
    
    if( Cut_Vertex_Not_Fiducial == 1) return false;
    cutList.nCut_Vertex_Not_Fiducial.increment(truth_isSignal, study1, study2);
    
    hCut_vertexCount->Fill(vtx_total_count);
    if( vtx_total_count > 1) return false;
    cutList.nCut_Vertex_Count.increment(truth_isSignal, study1, study2);
    
    hCut_nProngs->Fill(nProngs);
    if( Cut_nProngs == 1 || nProngs > 2) return false;
    cutList.nCut_nProngs.increment(truth_isSignal, study1, study2);
    
    // -------------------------------------------------------------------------
    // Separate Topologies
    //      Different Statistics for 1 Prong and 2 Prong Events
    // -------------------------------------------------------------------------
    if( Cut_Muon_None == 1) return false;
    if (nProngs == 1) cutList.nCut_1Prong_Muon_None.increment(truth_isSignal, study1, study2);
    if (nProngs == 2) cutList.nCut_2Prong_Muon_None.increment(truth_isSignal, study1, study2);
    
    if (truth_isSignal) fill_mc_w(); // Fill W for MINOS Matched Signal Events
        
    if( Cut_Muon_Not_Plausible == 1) return false;
    if (nProngs == 1) cutList.nCut_1Prong_Muon_Not_Plausible.increment(truth_isSignal, study1, study2);
    if (nProngs == 2) cutList.nCut_2Prong_Muon_Not_Plausible.increment(truth_isSignal, study1, study2);
        
    if( Cut_Muon_Score_Low == 1) return false;
    if (nProngs == 1) cutList.nCut_1Prong_Muon_Score_Low.increment(truth_isSignal, study1, study2);
    if (nProngs == 2) cutList.nCut_2Prong_Muon_Score_Low.increment(truth_isSignal, study1, study2);
    
    if( Cut_Muon_Charge == 1) return false;
    if (nProngs == 1) cutList.nCut_1Prong_Muon_Charge.increment(truth_isSignal, study1, study2);
    if (nProngs == 2) cutList.nCut_2Prong_Muon_Charge.increment(truth_isSignal, study1, study2);
    
    // Cout Events before Michel Selections
    if(truth_isBckg_withMichel) N_trueMichel_before++;
    else N_noMichel_before++;
    
    // Fill  hCut_Michel
    if (nProngs == 1){
        if( Cut_Vertex_Michel_Exist == 1 || Cut_EndPoint_Michel_Exist == 1 || Cut_secEndPoint_Michel_Exist == 1 ){
            hCut_1Prong_Michel->Fill(1);
        }else{
            hCut_1Prong_Michel->Fill(0);
        }
    }else if (nProngs == 2){
        if( Cut_Vertex_Michel_Exist == 1 || Cut_EndPoint_Michel_Exist == 1 || Cut_secEndPoint_Michel_Exist == 1 ){
            hCut_2Prong_Michel->Fill(1);
        }else{
            hCut_2Prong_Michel->Fill(0);
        }
    }
    
    if( Cut_Vertex_Michel_Exist == 1){
        if(isMichelStudy) michelFound_Vertex = true;
        else return false;
    }
    if (nProngs == 1) cutList.nCut_1Prong_Vertex_Michel_Exist.increment(truth_isSignal, study1, study2);
    if (nProngs == 2) cutList.nCut_2Prong_Vertex_Michel_Exist.increment(truth_isSignal, study1, study2);
    
    if( Cut_EndPoint_Michel_Exist == 1){
        if(isMichelStudy) michelFound_Track = true;
        else return false;
    }
    if (nProngs == 1) cutList.nCut_1Prong_EndPoint_Michel_Exist.increment(truth_isSignal, study1, study2);
    if (nProngs == 2) cutList.nCut_2Prong_EndPoint_Michel_Exist.increment(truth_isSignal, study1, study2);
    
    if( Cut_secEndPoint_Michel_Exist == 1){
        if(isMichelStudy) michelFound_Track2 = true;
        else return false;
    }  
    if (nProngs == 1) cutList.nCut_1Prong_secEndPoint_Michel_Exist.increment(truth_isSignal, study1, study2);
    if (nProngs == 2) cutList.nCut_2Prong_secEndPoint_Michel_Exist.increment(truth_isSignal, study1, study2);
    
    if(truth_isBckg_withMichel) N_trueMichel_after++;
    else N_noMichel_after++;
       
    // Fill PreFilter Plots
    if (nProngs == 1){
        hCut_1Prong_eVis_nuclearTarget->Fill(evis_NuclearTarget);
        hCut_1Prong_eVis_other->Fill(evis_TotalExceptNuclearTarget);
    }else if(nProngs == 2){
        hCut_2Prong_eVis_nuclearTarget->Fill(evis_NuclearTarget);
        hCut_2Prong_eVis_other->Fill(evis_TotalExceptNuclearTarget);
    }
    if( Cut_PreFilter_Pi0 == 1) return false;
    if (nProngs == 1) cutList.nCut_1Prong_PreFilter_Pi0.increment(truth_isSignal, study1, study2);
    if (nProngs == 2) cutList.nCut_2Prong_PreFilter_Pi0.increment(truth_isSignal, study1, study2); 
    
    if( Cut_VtxBlob == 1) return false;
    if (nProngs == 1) cutList.nCut_1Prong_VtxBlob.increment(truth_isSignal, study1, study2);
    if (nProngs == 2) cutList.nCut_2Prong_VtxBlob.increment(truth_isSignal, study1, study2);
    
    if( Cut_ConeBlobs == 1) return false;
    if (nProngs == 1) cutList.nCut_1Prong_ConeBlobs.increment(truth_isSignal, study1, study2);
    if (nProngs == 2) cutList.nCut_2Prong_ConeBlobs.increment(truth_isSignal, study1, study2);

    //  Photon Conversion Length
    if (nProngs == 1) hCut_1Prong_gamma1ConvDist->Fill(gamma1_dist_vtx * 0.1);
    if (nProngs == 2) hCut_2Prong_gamma1ConvDist->Fill(gamma1_dist_vtx * 0.1);
    if (applyPhotonDistance && gamma1_dist_vtx * 0.1 < minPhotonDistance) return false;
    if (nProngs == 1) cutList.nCut_1Prong_Photon1DistanceLow.increment(truth_isSignal, study1, study2);
    if (nProngs == 2) cutList.nCut_2Prong_Photon1DistanceLow.increment(truth_isSignal, study1, study2);

    if (nProngs == 1) hCut_1Prong_gamma2ConvDist->Fill(gamma2_dist_vtx * 0.1);
    if (nProngs == 2) hCut_2Prong_gamma2ConvDist->Fill(gamma2_dist_vtx * 0.1);
    if (applyPhotonDistance && gamma2_dist_vtx * 0.1 < minPhotonDistance) return false;
    if (nProngs == 1) cutList.nCut_1Prong_Photon2DistanceLow.increment(truth_isSignal, study1, study2);
    if (nProngs == 2) cutList.nCut_2Prong_Photon2DistanceLow.increment(truth_isSignal, study1, study2);


    // Pi0 Invariant Mass
    if(nProngs == 1) hCut_1Prong_pi0invMass->Fill(pi0_invMass);
    if(nProngs == 2) hCut_2Prong_pi0invMass->Fill(pi0_invMass);
    if( pi0_invMass < min_Pi0_invMass || pi0_invMass > max_Pi0_invMass ) return false;
    if(nProngs == 1) cutList.nCut_1Prong_Pi0_invMass.increment(truth_isSignal, study1, study2);
    if(nProngs == 2) cutList.nCut_2Prong_Pi0_invMass.increment(truth_isSignal, study1, study2);
    
    // Proton Related Cuts only for nProngs == 2
    if ( nProngs == 2 ){
        if( Cut_Particle_None == 1) return false;
        cutList.nCut_2Prong_Particle_None.increment(truth_isSignal, study1, study2);
        
        if( Cut_Proton_None == 1) return false;
        cutList.nCut_2Prong_Proton_None.increment(truth_isSignal, study1, study2);
        
        // Find Best Proton in Reco
        findRecoProton(); 
        
        // Fill Cut Histograms
        if(CCProtonPi0_proton_ekin[indRecoProton] < pID_KE_Limit){ 
            hCut_pIDDiff->Fill(CCProtonPi0_protonScore[indRecoProton] - CCProtonPi0_pionScore[indRecoProton]);
        }else{
            hCut_protonScore_LLR->Fill(CCProtonPi0_protonScore_LLR[indRecoProton]);
        }
        if ( applyProtonScore ){
            // Use pID Difference for KE < pID_KE_Limit 
            // Use LLR for KE > pID_KE_Limit
            if (CCProtonPi0_proton_ekin[indRecoProton] < pID_KE_Limit ){
                if ( (CCProtonPi0_protonScore[indRecoProton] - CCProtonPi0_pionScore[indRecoProton]) < minPIDDiff ) return false;
            }else{
                if ( CCProtonPi0_protonScore_LLR[indRecoProton] < minProtonScore_LLR ) return false;
            }
        }
        cutList.nCut_2Prong_ProtonScore.increment(truth_isSignal, study1, study2);
        
        double delta_invMass = calcDeltaInvariantMass();
        hCut_deltaInvMass->Fill(delta_invMass);
        if (applyDeltaInvMass){
            if ( delta_invMass < min_Delta_invMass || delta_invMass > max_Delta_invMass) return false;  
        }
        cutList.nCut_2Prong_DeltaInvMass.increment(truth_isSignal, study1, study2);
    }
    
    //==========================================================================
    //
    // Analysis Cuts - Improves Purity Dramatically
    //
    //==========================================================================
    
    if(nProngs == 1) hCut_1Prong_neutrinoE->Fill(CCProtonPi0_neutrino_E * HEP_Functions::MeV_to_GeV);
    if(nProngs == 2) hCut_2Prong_neutrinoE->Fill(CCProtonPi0_neutrino_E * HEP_Functions::MeV_to_GeV);
    if( applyBeamEnergy && ((CCProtonPi0_neutrino_E * HEP_Functions::MeV_to_GeV) > max_beamEnergy)) return false;
    if(nProngs == 1) cutList.nCut_1Prong_beamEnergy.increment(truth_isSignal, study1, study2);
    if(nProngs == 2) cutList.nCut_2Prong_beamEnergy.increment(truth_isSignal, study1, study2);
    
    if(nProngs == 1) hCut_1Prong_UnusedE->Fill(energyUnused_afterReco);
    if(nProngs == 2) hCut_2Prong_UnusedE->Fill(energyUnused_afterReco);
    if( applyUnusedE && (energyUnused_afterReco > maxUnusedE)) return false;
    if(nProngs == 1) cutList.nCut_1Prong_UnusedE.increment(truth_isSignal, study1, study2);
    if(nProngs == 2) cutList.nCut_2Prong_UnusedE.increment(truth_isSignal, study1, study2);
        

    // -------------------------------------------------------------------------
    // Michel Study
    //      Save information for events have true michel which are failed 
    //          to be found with Michel Tool
    // -------------------------------------------------------------------------
    if(truth_isBckg_withMichel){ 
        if(truth_michelMuon_end_dist_vtx > 200){
            writeScanFile();
        }
        
        int ind;
        
        if(michelFound_Vertex) ind = 0;
        else if(michelFound_Track) ind = 1;
        else if(michelFound_Track2) ind = 2;
        else ind = 3;

        N_michelElectrons->Fill(truth_N_trueMichelElectrons);
        // All Found Events
        if(michelFound_Vertex || michelFound_Track || michelFound_Track2 ){
            michelElectron_E[4]->Fill(truth_michelElectron_E);  
        }
        michelElectron_E[ind]->Fill(truth_michelElectron_E);
        michelMuon_X_Y[ind]->Fill(truth_michelMuon_endPoint[0],truth_michelMuon_endPoint[1]);
        michelMuon_Z[ind]->Fill(truth_michelMuon_endPoint[2]);
        michelMuon_Z_vtx[ind]->Fill(truth_michelMuon_endPoint[2]-mc_vtx[2]);
        michelMuon_P[ind]->Fill(truth_michelMuon_P);
        michelMuon_end_dist_vtx[ind]->Fill(truth_michelMuon_end_dist_vtx);
        michelMuon_length[ind]->Fill(truth_michelMuon_length);
        michelPion_P[ind]->Fill(truth_michelPion_P);
        michelPion_begin_dist_vtx[ind]->Fill(truth_michelPion_begin_dist_vtx);
        michelPion_length[ind]->Fill(truth_michelPion_length); 
        michelMuon_dist_michelPion_length[ind]->Fill(truth_michelMuon_end_dist_vtx,truth_michelPion_length);

    } 
    
    if(truth_isBckg_withMichel) N_trueMichel_afterAll++;
    return true;
    
}

void Analyzer::fill_mc_w()
{
    if(mc_intType == 1) mc_w_CCQE->Fill(mc_w * HEP_Functions::MeV_to_GeV);
    if(mc_intType == 2) mc_w_RES->Fill(mc_w * HEP_Functions::MeV_to_GeV);
    if(mc_intType == 3) mc_w_DIS->Fill(mc_w * HEP_Functions::MeV_to_GeV);
}

void Analyzer::formBackgroundVectors()
{
    bckgVector_1Prong.push_back(bckg_1Prong_QELike);
    bckgVector_1Prong.push_back(bckg_1Prong_SinglePiPlus); 
    bckgVector_1Prong.push_back(bckg_1Prong_SinglePiMinus); 
    bckgVector_1Prong.push_back(bckg_1Prong_MultiPion); 
    bckgVector_1Prong.push_back(bckg_1Prong_MultiPiZero);
    bckgVector_1Prong.push_back(bckg_1Prong_Other);
    
    bckgVector_2Prong.push_back(bckg_2Prong_QELike);
    bckgVector_2Prong.push_back(bckg_2Prong_SinglePiPlus); 
    bckgVector_2Prong.push_back(bckg_2Prong_SinglePiMinus); 
    bckgVector_2Prong.push_back(bckg_2Prong_MultiPion); 
    bckgVector_2Prong.push_back(bckg_2Prong_MultiPiZero);
    bckgVector_2Prong.push_back(bckg_2Prong_Other);
}

void Analyzer::writeBackgroundTableHeader()
{
        for (int t = 0; t < nTopologies; t++){ 
            backgroundText[t]<<std::left;
            // Table Header
            backgroundText[t].width(20); backgroundText[t]<< "Background Type"; backgroundText[t]<<"|";
            backgroundText[t].width(16); backgroundText[t]<< "Total"; backgroundText[t]<<"|";
            backgroundText[t].width(16); backgroundText[t]<< "AntiMuon"; backgroundText[t]<<"|";
            backgroundText[t].width(16); backgroundText[t]<< "Michel"; backgroundText[t]<<"|";
            backgroundText[t].width(16); backgroundText[t]<< "Primary Pi0"; backgroundText[t]<<"|";
            backgroundText[t].width(16); backgroundText[t]<< "Secondary Pi0"; backgroundText[t]<<"|";
            backgroundText[t].width(16); backgroundText[t]<< "Pi0 Candidate"; backgroundText[t]<<"|";
            backgroundText[t]<<endl;
        }
}


void Analyzer::writeBackgroundTableRows(vector< vector<double> > bckgVector, int nProngs)
{
    int t;
    if( nProngs == 1 ) t = 0;
    if( nProngs == 2 ) t = 1;
    
    vector<double> bckg_QELike = bckgVector[0];
    vector<double> bckg_SinglePiPlus = bckgVector[1];
    vector<double> bckg_SinglePiMinus = bckgVector[2];
    vector<double> bckg_MultiPion = bckgVector[3];
    vector<double> bckg_MultiPiZero = bckgVector[4];
    vector<double> bckg_Other = bckgVector[5];
    
    backgroundText[t].width(20); backgroundText[t]<<"QE Like"; backgroundText[t]<<"|"; 
    for(int i = 0; i < nBckgBranch; i++){
        backgroundText[t].width(16); backgroundText[t]<<bckg_QELike[i]; backgroundText[t]<<"|";   
    }
    backgroundText[t]<<endl;
    
    backgroundText[t].width(20); backgroundText[t]<<"SinglePiPlus"; backgroundText[t]<<"|"; 
    for(int i = 0; i < nBckgBranch; i++){
        backgroundText[t].width(16); backgroundText[t]<<bckg_SinglePiPlus[i]; backgroundText[t]<<"|";    
    }
    backgroundText[t]<<endl;
    
    backgroundText[t].width(20); backgroundText[t]<<"SinglePiMinus"; backgroundText[t]<<"|"; 
    for(int i = 0; i < nBckgBranch; i++){
        backgroundText[t].width(16); backgroundText[t]<<bckg_SinglePiMinus[i]; backgroundText[t]<<"|";    
    }
    backgroundText[t]<<endl;

    backgroundText[t].width(20); backgroundText[t]<<"MultiChargedPion"; backgroundText[t]<<"|"; 
    for(int i = 0; i < nBckgBranch; i++){
        backgroundText[t].width(16); backgroundText[t]<<bckg_MultiPion[i]; backgroundText[t]<<"|";    
    }
    backgroundText[t]<<endl;
    
    backgroundText[t].width(20); backgroundText[t]<<"MultiPiZero"; backgroundText[t]<<"|"; 
    for(int i = 0; i < nBckgBranch; i++){
        backgroundText[t].width(16); backgroundText[t]<<bckg_MultiPiZero[i]; backgroundText[t]<<"|";    
    }
    backgroundText[t]<<endl;
    
    backgroundText[t].width(20); backgroundText[t]<<"Other"; backgroundText[t]<<"|"; 
    for(int i = 0; i < nBckgBranch; i++){
        backgroundText[t].width(16); backgroundText[t]<<bckg_Other[i]; backgroundText[t]<<"|";  
    }
    backgroundText[t]<<endl;
    
}

void Analyzer::writeBackgroundTable()
{
    formBackgroundVectors();
    writeBackgroundTableHeader();
    writeBackgroundTableRows(bckgVector_1Prong, 1);
    writeBackgroundTableRows(bckgVector_2Prong, 2);
}

void Analyzer::fillInteractionTrue()
{
    
    int_channel->Fill(mc_intType);

    vertex_z_true->Fill(mc_vtx[2]);
    vertex_z_error->Fill( Data_Functions::getError(mc_vtx[2],CCProtonPi0_vtx[2]) );
    vertex_z_reco_mc->Fill(CCProtonPi0_vtx[2],mc_vtx[2]);
    
    vertex_x_y_true->Fill(mc_vtx[0],mc_vtx[1]);
    
}

void Analyzer::fillInteractionReco()
{
    if (nProngs == 2) deltaInvMass_reco->Fill(calcDeltaInvariantMass());

    beamEnergy_mc->Fill(mc_incomingE * HEP_Functions::MeV_to_GeV);    
    beamEnergy_reco->Fill(CCProtonPi0_neutrino_E * HEP_Functions::MeV_to_GeV);
    beamEnergy_error->Fill( Data_Functions::getError(mc_incomingE,CCProtonPi0_neutrino_E) );
    beamEnergy_reco_mc->Fill(CCProtonPi0_neutrino_E * HEP_Functions::MeV_to_GeV,mc_incomingE * HEP_Functions::MeV_to_GeV);
    
    beamEnergyCal_mc->Fill(mc_incomingE * HEP_Functions::MeV_to_GeV);    
    beamEnergyCal_reco->Fill(CCProtonPi0_neutrino_E_Cal * HEP_Functions::MeV_to_GeV);
    beamEnergyCal_error->Fill( Data_Functions::getError(mc_incomingE,CCProtonPi0_neutrino_E_Cal) );
    beamEnergyCal_reco_mc->Fill(CCProtonPi0_neutrino_E_Cal * HEP_Functions::MeV_to_GeV,mc_incomingE * HEP_Functions::MeV_to_GeV);
    
    beamEnergy_beamEnergyCal->Fill(CCProtonPi0_neutrino_E * HEP_Functions::MeV_to_GeV, CCProtonPi0_neutrino_E_Cal * HEP_Functions::MeV_to_GeV);

    q2_mc->Fill(mc_Q2 * HEP_Functions::MeVSq_to_GeVSq);
    q2_reco->Fill(CCProtonPi0_QSq * HEP_Functions::MeVSq_to_GeVSq);
    q2_error->Fill( Data_Functions::getError(mc_Q2,CCProtonPi0_QSq) );
    q2_reco_mc->Fill(CCProtonPi0_QSq * HEP_Functions::MeVSq_to_GeVSq,mc_Q2* HEP_Functions::MeVSq_to_GeVSq);
    
    w_mc->Fill(mc_w * HEP_Functions::MeV_to_GeV);
    w_reco->Fill(CCProtonPi0_W* HEP_Functions::MeV_to_GeV);
    w_error->Fill( Data_Functions::getError(mc_w,CCProtonPi0_W) );
    w_reco_mc->Fill(CCProtonPi0_W* HEP_Functions::MeV_to_GeV,mc_w * HEP_Functions::MeV_to_GeV);
    
    wSq_reco->Fill(CCProtonPi0_WSq * HEP_Functions::MeVSq_to_GeVSq);
    
    vertex_count->Fill(vtx_total_count);
    
    vertex_x_y_reco->Fill(CCProtonPi0_vtx[0],CCProtonPi0_vtx[1]);
    vertex_z_reco->Fill(CCProtonPi0_vtx[2]);
    
    nProngs_hist->Fill(nProngs);
    
    total_E->Fill(CCProtonPi0_total_E * HEP_Functions::MeV_to_GeV);
    total_E_neutrinoE->Fill(CCProtonPi0_total_E * HEP_Functions::MeV_to_GeV, CCProtonPi0_neutrino_E * HEP_Functions::MeV_to_GeV);
    
//     time_AllClusters->Fill(AllClustersTime);
    
    E_Unused_afterReco->Fill(energyUnused_afterReco);
    E_Used_afterReco->Fill(energyUsed_afterReco);
    
    if(truth_isSignal){
        if(mc_intType == 1) final_mc_w_CCQE->Fill(mc_w * HEP_Functions::MeV_to_GeV);
        if(mc_intType == 2) final_mc_w_RES->Fill(mc_w * HEP_Functions::MeV_to_GeV);
        if(mc_intType == 3) final_mc_w_DIS->Fill(mc_w * HEP_Functions::MeV_to_GeV);
    }
    
//     if (nProngs == 2){ 
//         cout<<prong_showerScore<<endl;
//     }
}

double Analyzer::calcDeltaInvariantMass()
{
    double invMassSq;
        
    invMassSq = (pi0_E + CCProtonPi0_proton_E[indRecoProton]) * (pi0_E + CCProtonPi0_proton_E[indRecoProton]) -
                ((pi0_px + CCProtonPi0_proton_px[indRecoProton])*(pi0_px + CCProtonPi0_proton_px[indRecoProton]) + 
                 (pi0_py + CCProtonPi0_proton_py[indRecoProton])*(pi0_py + CCProtonPi0_proton_py[indRecoProton]) +
                 (pi0_pz + CCProtonPi0_proton_pz[indRecoProton])*(pi0_pz + CCProtonPi0_proton_pz[indRecoProton]));
        
    return sqrt(invMassSq);
}



void Analyzer::writeScanFile()
{
    // Constants for Roundup List
    const string arachne_html = "http://minerva05.fnal.gov/Arachne/arachne.html?filename=";
    const string entryString  = "&entry=";
    const string other        = "&slice=-1&filetype=dst";
    
    roundupText<<arachne_html<<scanFileName<<entryString<<truth_eventID<<other<<" ";
    roundupText<<truth_eventID<<"^"<<truth_michelMuon_end_dist_vtx<<"^s"<<mc_incomingE<<endl;
}


void Analyzer::fillHistograms()
{
    muon.fill_Histograms();
    pion.fill_Histograms();
    proton.fill_Histograms();
}

void Analyzer::write_RootFile()
{
    cout<<">> Writing "<<rootDir<<endl;
    f->Write();
}

void Analyzer::closeTextFiles()
{
    readme.close(); 
    
    for (int i = 0; i < nTopologies; i++){
        backgroundText[i].close();
        failText[i].close();
    }
}

void Analyzer::openTextFiles()
{
    
    cout<<"\tText Files for Output:"<<endl;
    // Open Readme File
    readmeFile = Folder_List::output + Folder_List::textOut + branchDir + "readme.txt";
    readme.open( readmeFile.c_str() );
    if( !readme.is_open() ){
        cerr<<"Cannot open output text file: "<<readmeFile<<endl;
        exit(1);
    }else{
        cout<<"\t"<<readmeFile<<endl;
    }
    writeReadme();
    

    
    // Open Background Files
    backgroundFile[0] = Folder_List::output + Folder_List::textOut + branchDir + "BackgroundTable_1Prong.txt";
    backgroundFile[1] = Folder_List::output + Folder_List::textOut + branchDir + "BackgroundTable_2Prong.txt";
    
    for (int i = 0; i < nTopologies; i++){
        backgroundText[i].open( backgroundFile[i].c_str() );
        if( !backgroundText[i].is_open() ){
            cerr<<"Cannot open output text file: "<<backgroundFile[i]<<endl;
            exit(1);
        }else{
            cout<<"\t"<<backgroundFile[i]<<endl;
        }
    }
    
    // Open Fail-Check Files
    failFile[0] = Folder_List::output + Folder_List::textOut + branchDir + "FailChecks_1Prong.txt";
    failFile[1] = Folder_List::output + Folder_List::textOut + branchDir + "FailChecks_2Prong.txt";
    
    for (int i = 0; i < nTopologies; i++){
        failText[i].open( failFile[i].c_str() );
        if( !failText[i].is_open() ){
            cerr<<"Cannot open output text file: "<<failFile[i]<<endl;
            exit(1);
        }else{
            cout<<"\t"<<failFile[i]<<endl;
        }
    }
    
    // Open Roundup Text for Arachne Scanning
    string roundupFile = Folder_List::output + Folder_List::textOut + branchDir + "ArachneRoundup.txt";
    roundupText.open(roundupFile.c_str() );
    if( !roundupText.is_open() ){
        cerr<<"Cannot open output text file: "<<roundupFile<<endl;
        exit(1);
    }
    
    string playlistDST = "Input/Playlists/pl_Scan.dat";
    DSTFileList.open( playlistDST.c_str() );
    if( !DSTFileList.is_open() ){
        cerr<<"Cannot open input text file: "<<playlistDST<<endl;
        exit(1);
    }
}

void Analyzer::writeReadme()
{
    readme<<"Test"<<endl;
}



#endif //CCProtonPi0_cpp
