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
    isDataAnalysis  = false;
    isScanRun = false;
    is_pID_Studies  = false;
    isMichelStudy = false;
    applyOtherMichelCuts = false;
    writeFSParticleMomentum = false;
    
    applyProtonScore = true;
    pID_KE_Limit = 300.0;
    minProtonScore_LLR = 10.0;
    minPIDDiff = 0.45;
    
    applyPhotonDistance = true;
    minPhotonDistance = 15; //cm
    
    applyBeamEnergy = true;
    max_beamEnergy = 10.0; // GeV
    
    applyUnusedE = true;
    maxUnusedE = 300;
    
    min_Pi0_invMass = 40.0;
    max_Pi0_invMass = 200.0;
    
    applyDeltaInvMass = false;
    min_Delta_invMass = 40.0;
    max_Delta_invMass = 200.0;
        
    latest_ScanID = 0.0;
    
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
    fChain->SetBranchStatus("truth_*",1);
    fChain->SetBranchStatus("mc_*",1);
    fChain->SetBranchStatus("CCProtonPi0*",1);
    fChain->SetBranchStatus("preFilter_*",1);
    fChain->SetBranchStatus("nProngs",1);
    fChain->SetBranchStatus("vtx*",1);
    fChain->SetBranchStatus("evis_*",1);
    fChain->SetBranchStatus("energyUnused*",1);
    fChain->SetBranchStatus("energyUsed*",1);
    fChain->SetBranchStatus("CCProtonPi0_total_*",1);
    fChain->SetBranchStatus("pi0_*",1);
    fChain->SetBranchStatus("gamma*",1);
    fChain->SetBranchStatus("g1blob*",1);
    fChain->SetBranchStatus("g2blob*",1);
    fChain->SetBranchStatus("michel*",1);
    fChain->SetBranchStatus("endpoint_michel_distance",1);
    
    
    
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
        
        // Update scanFileName if running for scan
        if(isScanRun) UpdateScanFileName();
       
        // Decide to Analyze Event or NOT
        if( !analyzeEvent() ) continue;
        
        //----------------------------------------------------------------------
        // Get Cut Statistics
        //----------------------------------------------------------------------
        isPassedAllCuts = getCutStatistics();
        if( !isPassedAllCuts ) continue;
        
        //----------------------------------------------------------------------
        // Fill Background Branches for Background Events
        //----------------------------------------------------------------------
        if(!truth_isSignal) {
            bckgTool.fillBackgroundBranches(    nProngs,
                                                truth_isBckg_QELike, 
                                                truth_isBckg_SinglePiPlus, 
                                                truth_isBckg_SinglePiMinus, 
                                                truth_isBckg_MultiPion, 
                                                truth_isBckg_MultiPiZero, 
                                                truth_isBckg_Other,
                                                truth_isBckg_withAntiMuon,
                                                truth_isBckg_withMichel,
                                                truth_isBckg_withPrimaryPi0,
                                                truth_isBckg_withSecondaryPi0);
        }

        //----------------------------------------------------------------------
        // Data Analysis and Other Studies
        //----------------------------------------------------------------------
        if (isDataAnalysis) fillData();
        if (is_pID_Studies && nProngs == 2){ 
            pIDTool.FillHistograms(CCProtonPi0_protonScore_LLR[indRecoProton],CCProtonPi0_protonScore[indRecoProton],CCProtonPi0_pionScore[indRecoProton],
                                    CCProtonPi0_trajProtonProngPDG[indRecoProton],CCProtonPi0_proton_E[indRecoProton] );
        }
        if (writeFSParticleMomentum) writeFSParticle4P(jentry);
        

    } // end for-loop
    
    
    //--------------------------------------------------------------------------
    // Studies
    //--------------------------------------------------------------------------
    if(is_pID_Studies) pIDTool.get_pID_Stats();
    
    //--------------------------------------------------------------------------
    // Write Text Files
    //--------------------------------------------------------------------------
    cutList.writeCutTable();
    bckgTool.writeBackgroundTable();
    getPi0Family(); // Pi0 Family Information written inside FailFile
    
    //--------------------------------------------------------------------------
    // Write Root Files
    //--------------------------------------------------------------------------
    interaction.write_RootFile();
    muon.write_RootFile();
    proton.write_RootFile();
    pion.write_RootFile();
    if(is_pID_Studies) pIDTool.write_RootFile();
    if(isMichelStudy) michelTool.write_RootFile();
    
    
}

//------------------------------------------------------------------------------
//  Constructor
//------------------------------------------------------------------------------
Analyzer::Analyzer(int nMode) : 
    NTupleAnalysis(nMode),
    interaction(nMode),
    muon(nMode),
    proton(nMode),
    pion(nMode),
    pIDTool(nMode),
    bckgTool(nMode),
    michelTool(nMode),
    cutList(nMode)
{    
    openTextFiles();
    
    specifyRunTime();
    
    cout<<"Initialization Finished!\n"<<endl;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
//
//     Specific Functions
//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

void Analyzer::UpdateScanFileName()
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
    
    interaction.hCut_vertexCount->Fill(vtx_total_count);
    if( vtx_total_count > 1) return false;
    cutList.nCut_Vertex_Count.increment(truth_isSignal, study1, study2);
    
    interaction.hCut_nProngs->Fill(nProngs);
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
    
    // Fill  hCut_Michel
    if (nProngs == 1){
        if( Cut_Vertex_Michel_Exist == 1 || Cut_EndPoint_Michel_Exist == 1 || Cut_secEndPoint_Michel_Exist == 1 ){
            interaction.hCut_1Prong_Michel->Fill(1);
        }else{
            interaction.hCut_1Prong_Michel->Fill(0);
        }
    }else if (nProngs == 2){
        if( Cut_Vertex_Michel_Exist == 1 || Cut_EndPoint_Michel_Exist == 1 || Cut_secEndPoint_Michel_Exist == 1 ){
            interaction.hCut_2Prong_Michel->Fill(1);
        }else{
            interaction.hCut_2Prong_Michel->Fill(0);
        }
    }
    
    if( Cut_Vertex_Michel_Exist == 1){
        if( !isMichelStudy ) return false;
    }
    if (nProngs == 1) cutList.nCut_1Prong_Vertex_Michel_Exist.increment(truth_isSignal, study1, study2);
    if (nProngs == 2) cutList.nCut_2Prong_Vertex_Michel_Exist.increment(truth_isSignal, study1, study2);
    
    if( Cut_EndPoint_Michel_Exist == 1){
        if( !isMichelStudy ) return false;
    }
    if (nProngs == 1) cutList.nCut_1Prong_EndPoint_Michel_Exist.increment(truth_isSignal, study1, study2);
    if (nProngs == 2) cutList.nCut_2Prong_EndPoint_Michel_Exist.increment(truth_isSignal, study1, study2);
    
    if( Cut_secEndPoint_Michel_Exist == 1){
        if( !isMichelStudy ) return false;
    }  
    if (nProngs == 1) cutList.nCut_1Prong_secEndPoint_Michel_Exist.increment(truth_isSignal, study1, study2);
    if (nProngs == 2) cutList.nCut_2Prong_secEndPoint_Michel_Exist.increment(truth_isSignal, study1, study2);
    
    
    bool isMichelFound =    (Cut_Vertex_Michel_Exist == 1) || 
                            (Cut_EndPoint_Michel_Exist == 1) ||
                            (Cut_secEndPoint_Michel_Exist == 1);
    
    // Apply Other Michel Cuts if michelFound
    // If KeepMichelEvent() Returns False --> Cut the Event
    if(isMichelStudy && (applyOtherMichelCuts && isMichelFound) ){ 
        if( !KeepMichelEvent() ) return false;
    }
    
    // Save Found Michel Information
    if(isMichelStudy && isMichelFound){
        if (truth_isBckg_withMichel){ 
            if (Cut_Vertex_Michel_Exist == 1 ) michelTool.trueMichel_dist_vtx->Fill(vtx_michel_distance);
            if (Cut_Vertex_Michel_Exist == 1 ) michelTool.trueMichel_end_Z_vtx_Z->Fill(michelProng_end_Z-CCProtonPi0_vtx_z);
            if (Cut_EndPoint_Michel_Exist == 1) michelTool.trueMichel_dist_end_point->Fill(endpoint_michel_distance);
            
            michelTool.trueMichel_end_Z->Fill(michelProng_end_Z);
            michelTool.trueMichel_energy->Fill(michelProng_energy);
            michelTool.trueMichel_time_diff->Fill(michelProng_time_diff);
        }
        else{
            if (Cut_Vertex_Michel_Exist == 1 ) michelTool.fakeMichel_dist_vtx->Fill(vtx_michel_distance);
            if (Cut_Vertex_Michel_Exist == 1 ) michelTool.fakeMichel_end_Z_vtx_Z->Fill(michelProng_end_Z-CCProtonPi0_vtx_z);
            if (Cut_EndPoint_Michel_Exist == 1) michelTool.fakeMichel_dist_end_point->Fill(endpoint_michel_distance);
            
            michelTool.fakeMichel_end_Z->Fill(michelProng_end_Z);
            michelTool.fakeMichel_energy->Fill(michelProng_energy);
            michelTool.fakeMichel_time_diff->Fill(michelProng_time_diff);
        }
    }
    
    // Cut Events with Michel after the Study
    if (isMichelStudy && isMichelFound) return false;
    
       
    // Fill PreFilter Plots
    if (nProngs == 1){
        interaction.hCut_1Prong_eVis_nuclearTarget->Fill(evis_NuclearTarget);
        interaction.hCut_1Prong_eVis_other->Fill(evis_TotalExceptNuclearTarget);
    }else if(nProngs == 2){
        interaction.hCut_2Prong_eVis_nuclearTarget->Fill(evis_NuclearTarget);
        interaction.hCut_2Prong_eVis_other->Fill(evis_TotalExceptNuclearTarget);
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
    if (nProngs == 1) interaction.hCut_1Prong_gamma1ConvDist->Fill(gamma1_dist_vtx * 0.1);
    if (nProngs == 2) interaction.hCut_2Prong_gamma1ConvDist->Fill(gamma1_dist_vtx * 0.1);
    if (applyPhotonDistance && gamma1_dist_vtx * 0.1 < minPhotonDistance) return false;
    if (nProngs == 1) cutList.nCut_1Prong_Photon1DistanceLow.increment(truth_isSignal, study1, study2);
    if (nProngs == 2) cutList.nCut_2Prong_Photon1DistanceLow.increment(truth_isSignal, study1, study2);

    if (nProngs == 1) interaction.hCut_1Prong_gamma2ConvDist->Fill(gamma2_dist_vtx * 0.1);
    if (nProngs == 2) interaction.hCut_2Prong_gamma2ConvDist->Fill(gamma2_dist_vtx * 0.1);
    if (applyPhotonDistance && gamma2_dist_vtx * 0.1 < minPhotonDistance) return false;
    if (nProngs == 1) cutList.nCut_1Prong_Photon2DistanceLow.increment(truth_isSignal, study1, study2);
    if (nProngs == 2) cutList.nCut_2Prong_Photon2DistanceLow.increment(truth_isSignal, study1, study2);


    // Pi0 Invariant Mass
    if(nProngs == 1) interaction.hCut_1Prong_pi0invMass->Fill(pi0_invMass);
    if(nProngs == 2) interaction.hCut_2Prong_pi0invMass->Fill(pi0_invMass);
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
            interaction.hCut_pIDDiff->Fill(CCProtonPi0_protonScore[indRecoProton] - CCProtonPi0_pionScore[indRecoProton]);
        }else{
            interaction.hCut_protonScore_LLR->Fill(CCProtonPi0_protonScore_LLR[indRecoProton]);
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
        interaction.hCut_deltaInvMass->Fill(delta_invMass);
        if (applyDeltaInvMass){
            if ( delta_invMass < min_Delta_invMass || delta_invMass > max_Delta_invMass) return false;  
        }
        cutList.nCut_2Prong_DeltaInvMass.increment(truth_isSignal, study1, study2);
    }
    
    if(nProngs == 1) interaction.hCut_1Prong_neutrinoE->Fill(CCProtonPi0_neutrino_E * HEP_Functions::MeV_to_GeV);
    if(nProngs == 2) interaction.hCut_2Prong_neutrinoE->Fill(CCProtonPi0_neutrino_E * HEP_Functions::MeV_to_GeV);
    if( applyBeamEnergy && ((CCProtonPi0_neutrino_E * HEP_Functions::MeV_to_GeV) > max_beamEnergy)) return false;
    if(nProngs == 1) cutList.nCut_1Prong_beamEnergy.increment(truth_isSignal, study1, study2);
    if(nProngs == 2) cutList.nCut_2Prong_beamEnergy.increment(truth_isSignal, study1, study2);
    
    if(nProngs == 1) interaction.hCut_1Prong_UnusedE->Fill(energyUnused_afterReco);
    if(nProngs == 2) interaction.hCut_2Prong_UnusedE->Fill(energyUnused_afterReco);
    if( applyUnusedE && (energyUnused_afterReco > maxUnusedE)) return false;
    if(nProngs == 1) cutList.nCut_1Prong_UnusedE.increment(truth_isSignal, study1, study2);
    if(nProngs == 2) cutList.nCut_2Prong_UnusedE.increment(truth_isSignal, study1, study2);
        

    // -------------------------------------------------------------------------
    // Michel Study
    //      Save information for events have true michel which are failed 
    //          to be found with Michel Tool
    // -------------------------------------------------------------------------
    if(isMichelStudy && truth_isBckg_withMichel){         
        int ind;
        
        if(Cut_Vertex_Michel_Exist == 1) ind = 0;
        else if(Cut_EndPoint_Michel_Exist == 1) ind = 1;
        else if(Cut_secEndPoint_Michel_Exist == 1) ind = 2;
        else ind = 3;

        michelTool.N_michelElectrons->Fill(truth_N_trueMichelElectrons);
        // All Found Events
        if(Cut_Vertex_Michel_Exist == 1 || Cut_EndPoint_Michel_Exist == 1 || Cut_secEndPoint_Michel_Exist == 1 ){
            michelTool.michelElectron_E[4]->Fill(truth_michelElectron_E);  
        }
        michelTool.michelElectron_E[ind]->Fill(truth_michelElectron_E);
        michelTool.michelMuon_X_Y[ind]->Fill(truth_michelMuon_endPoint[0],truth_michelMuon_endPoint[1]);
        michelTool.michelMuon_Z[ind]->Fill(truth_michelMuon_endPoint[2]);
        michelTool.michelMuon_Z_vtx[ind]->Fill(truth_michelMuon_endPoint[2]-mc_vtx[2]);
        michelTool.michelMuon_P[ind]->Fill(truth_michelMuon_P);
        michelTool.michelMuon_end_dist_vtx[ind]->Fill(truth_michelMuon_end_dist_vtx);
        michelTool.michelMuon_length[ind]->Fill(truth_michelMuon_length);
        michelTool.michelPion_P[ind]->Fill(truth_michelPion_P);
        michelTool.michelPion_begin_dist_vtx[ind]->Fill(truth_michelPion_begin_dist_vtx);
        michelTool.michelPion_length[ind]->Fill(truth_michelPion_length); 
        michelTool.michelMuon_dist_michelPion_length[ind]->Fill(truth_michelMuon_end_dist_vtx,truth_michelPion_length);

    } 
    
    return true;
    
}

void Analyzer::fill_mc_w()
{
    if(mc_intType == 1) interaction.mc_w_CCQE->Fill(mc_w * HEP_Functions::MeV_to_GeV);
    if(mc_intType == 2) interaction.mc_w_RES->Fill(mc_w * HEP_Functions::MeV_to_GeV);
    if(mc_intType == 3) interaction.mc_w_DIS->Fill(mc_w * HEP_Functions::MeV_to_GeV);
}

/*
    pre:
        Event with detected Michel by the Tool
    
    Applies series of other selections to save some Events
*/
bool Analyzer::KeepMichelEvent()
{
    // Apply Time Difference Cut - Save Events with Short time Difference
     double maxTimeDiff = 300.0;
     if(michelProng_time_diff < maxTimeDiff){ 
        // Apply Energy Cut
        if(michelProng_energy > 100) return true;  
     }
    
    // Event did not satisfied any of the conditions - will not Keep the Event
    return false;
    
}

void Analyzer::fillInteractionTrue()
{
    
    interaction.int_channel->Fill(mc_intType);

    interaction.vertex_z_true->Fill(mc_vtx[2]);
    interaction.vertex_z_error->Fill( Data_Functions::getError(mc_vtx[2],CCProtonPi0_vtx[2]) );
    interaction.vertex_z_reco_mc->Fill(CCProtonPi0_vtx[2],mc_vtx[2]);
    
    interaction.vertex_x_y_true->Fill(mc_vtx[0],mc_vtx[1]);
    
}

void Analyzer::fillInteractionReco()
{
    if (nProngs == 2) interaction.deltaInvMass_reco->Fill(calcDeltaInvariantMass());

    interaction.beamEnergy_mc->Fill(mc_incomingE * HEP_Functions::MeV_to_GeV);    
    interaction.beamEnergy_reco->Fill(CCProtonPi0_neutrino_E * HEP_Functions::MeV_to_GeV);
    interaction.beamEnergy_error->Fill( Data_Functions::getError(mc_incomingE,CCProtonPi0_neutrino_E) );
    interaction.beamEnergy_reco_mc->Fill(CCProtonPi0_neutrino_E * HEP_Functions::MeV_to_GeV,mc_incomingE * HEP_Functions::MeV_to_GeV);
    
    interaction.beamEnergyCal_mc->Fill(mc_incomingE * HEP_Functions::MeV_to_GeV);    
    interaction.beamEnergyCal_reco->Fill(CCProtonPi0_neutrino_E_Cal * HEP_Functions::MeV_to_GeV);
    interaction.beamEnergyCal_error->Fill( Data_Functions::getError(mc_incomingE,CCProtonPi0_neutrino_E_Cal) );
    interaction.beamEnergyCal_reco_mc->Fill(CCProtonPi0_neutrino_E_Cal * HEP_Functions::MeV_to_GeV,mc_incomingE * HEP_Functions::MeV_to_GeV);
    
    interaction.beamEnergy_beamEnergyCal->Fill(CCProtonPi0_neutrino_E * HEP_Functions::MeV_to_GeV, CCProtonPi0_neutrino_E_Cal * HEP_Functions::MeV_to_GeV);

    interaction.q2_mc->Fill(mc_Q2 * HEP_Functions::MeVSq_to_GeVSq);
    interaction.q2_reco->Fill(CCProtonPi0_QSq * HEP_Functions::MeVSq_to_GeVSq);
    interaction.q2_error->Fill( Data_Functions::getError(mc_Q2,CCProtonPi0_QSq) );
    interaction.q2_reco_mc->Fill(CCProtonPi0_QSq * HEP_Functions::MeVSq_to_GeVSq,mc_Q2* HEP_Functions::MeVSq_to_GeVSq);
    
    interaction.w_mc->Fill(mc_w * HEP_Functions::MeV_to_GeV);
    interaction.w_reco->Fill(CCProtonPi0_W* HEP_Functions::MeV_to_GeV);
    interaction.w_error->Fill( Data_Functions::getError(mc_w,CCProtonPi0_W) );
    interaction.w_reco_mc->Fill(CCProtonPi0_W* HEP_Functions::MeV_to_GeV,mc_w * HEP_Functions::MeV_to_GeV);
    
    interaction.wSq_reco->Fill(CCProtonPi0_WSq * HEP_Functions::MeVSq_to_GeVSq);
    
    interaction.vertex_count->Fill(vtx_total_count);
    
    interaction.vertex_x_y_reco->Fill(CCProtonPi0_vtx[0],CCProtonPi0_vtx[1]);
    interaction.vertex_z_reco->Fill(CCProtonPi0_vtx[2]);
    
    interaction.nProngs_hist->Fill(nProngs);
    
    interaction.total_E->Fill(CCProtonPi0_total_E * HEP_Functions::MeV_to_GeV);
    interaction.total_E_neutrinoE->Fill(CCProtonPi0_total_E * HEP_Functions::MeV_to_GeV, CCProtonPi0_neutrino_E * HEP_Functions::MeV_to_GeV);
    
//     time_AllClusters->Fill(AllClustersTime);
    
    interaction.E_Unused_afterReco->Fill(energyUnused_afterReco);
    interaction.E_Used_afterReco->Fill(energyUsed_afterReco);
    
    if(truth_isSignal){
        if(mc_intType == 1) interaction.final_mc_w_CCQE->Fill(mc_w * HEP_Functions::MeV_to_GeV);
        if(mc_intType == 2) interaction.final_mc_w_RES->Fill(mc_w * HEP_Functions::MeV_to_GeV);
        if(mc_intType == 3) interaction.final_mc_w_DIS->Fill(mc_w * HEP_Functions::MeV_to_GeV);
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
    if(isScanRun){
        // Constants for Roundup List
        const string arachne_html = "http://minerva05.fnal.gov/Arachne/arachne.html?filename=";
        const string entryString  = "&entry=";
        const string other        = "&slice=-1&filetype=dst";
        
        roundupText<<arachne_html<<scanFileName<<entryString<<truth_eventID<<other<<" ";
        roundupText<<truth_eventID<<"^"<<truth_michelMuon_end_dist_vtx<<"^"<<mc_incomingE<<endl;
    }else{
        cout<<"WARNING! ScanRun is NOT Activated! Are you sure what you are doing?"<<endl;    
    }
}


void Analyzer::fillHistograms()
{
    muon.fill_Histograms();
    pion.fill_Histograms();
    proton.fill_Histograms();
}

void Analyzer::closeTextFiles()
{
    logFile.close();
    
    if(isScanRun){
        roundupText.close();
        DSTFileList.close();
    }
    
    for (int i = 0; i < nTopologies; i++){
        failText[i].close();
    }
}

void Analyzer::openTextFiles()
{
    cout<<"Opening Text Files:"<<endl;
    
    logFileName = Folder_List::output + Folder_List::textOut + branchDir + "LogFile.txt";
    logFile.open(logFileName.c_str());
    if( !logFile.is_open() ){
        cerr<<"Cannot open output text file: "<<logFileName<<endl;
        exit(EXIT_FAILURE);    
    }else{
        cout<<"\t"<<logFileName<<endl;
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
    
    if(isScanRun){
        // Open Roundup Text for Arachne Scanning
        string roundupFile = Folder_List::output + Folder_List::textOut + branchDir + "ArachneRoundup.txt";
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

void Analyzer::fillProtonTrue()
{    
    // Fill 4-Momentum
    proton.set_p4(  truth_proton_px[indTrueProton],
                    truth_proton_py[indTrueProton],
                    truth_proton_pz[indTrueProton],
                    truth_proton_E[indTrueProton], 
                    true);
       
    // set Angle wrt Beam
    proton.set_angleBeam(beam_p3, true);
    
    // set Angle wrt Muon
    proton.set_angleMuon(muon, true);
    
}

void Analyzer::fillProtonReco()
{   
    // Set Particle Score
    proton.particleScore = CCProtonPi0_protonScore_LLR[indRecoProton];

    // Fill 4-Momentum
    proton.set_p4(  CCProtonPi0_proton_px[indRecoProton],
                    CCProtonPi0_proton_py[indRecoProton],
                    CCProtonPi0_proton_pz[indRecoProton],
                    CCProtonPi0_proton_E[indRecoProton],
                    false);
                    
    // set Angle wrt Beam
    proton.set_angleBeam(beam_p3, false);
    
    // set Angle wrt Muon
    proton.set_angleMuon(muon, false);
    
}

//--------------------------------------------------------------------------
//  findTrueProton()
//      finds Truth Proton using Proton Energy - Highest Energy Proton
//--------------------------------------------------------------------------
void Analyzer::findTrueProton()
{
    double tempInd = 0;
    double tempE = truth_proton_E[0];
    
    for ( int i = 1; i < 20; i++){
        if (truth_proton_E[i] == SENTINEL) break;
        if (truth_proton_E[i] > tempE){
            tempE = truth_proton_E[i];
            tempInd = 0;
        }
    }
    
    indTrueProton = tempInd;
}

//--------------------------------------------------------------------------
//  findRecoProton()
//      finds Reco Proton using Proton Score
//--------------------------------------------------------------------------
void Analyzer::findRecoProton()
{
    double tempScore = CCProtonPi0_protonScore_LLR[0];
    double currentScore;
    int tempInd = 0;
   
    for( int i = 0; i < 10; i++){
        currentScore = CCProtonPi0_protonScore_LLR[i];
        if( currentScore == SENTINEL ) break;
        if( currentScore > tempScore){
            tempScore = currentScore;
            tempInd = i;
        }
    }
    
    indRecoProton = tempInd;
}

void Analyzer::fillPionReco()
{
    double photon_E_asym;
    
    // Fill 4-Momentum
    pion.set_p4(    pi0_px,
                    pi0_py,
                    pi0_pz,
                    pi0_E,
                    false);
    
    
    // set Angle wrt Beam
    pion.set_angleBeam(beam_p3, false);
    
    // set Angle wrt Muon
    pion.set_angleMuon(muon, false);
    
    // Set Invariant Mass
    pion.invMass->Fill(pi0_invMass);
    
    // Set Photon Conversion Length in [cm]
    pion.gamma1_ConvLength->Fill(gamma1_dist_vtx * 0.1);
    pion.gamma2_ConvLength->Fill(gamma2_dist_vtx * 0.1);
    pion.ConvLength_gamma2_gamma1->Fill(gamma2_dist_vtx, gamma1_dist_vtx);
        
    // Set Photon N(Clusters)
    pion.gamma1_nClusters_All->Fill(gamma1_blob_nclusters);
    pion.gamma2_nClusters_All->Fill(gamma2_blob_nclusters);
    pion.nClusters_All_gamma2_gamma1->Fill(gamma2_blob_nclusters,gamma1_blob_nclusters);
    
    // Set Photon Energy [GeV]
    pion.gamma1_Energy->Fill(gamma1_E * HEP_Functions::MeV_to_GeV);
    pion.gamma2_Energy->Fill(gamma2_E * HEP_Functions::MeV_to_GeV);
    pion.Energy_gamma2_gamma1->Fill(gamma2_E * HEP_Functions::MeV_to_GeV,gamma1_E * HEP_Functions::MeV_to_GeV );
    
    // Set Photon Energy Assymmetry
    photon_E_asym = abs((gamma1_E - gamma2_E) / (gamma1_E + gamma2_E));  
    pion.photonEnergy_Asymmetry->Fill(photon_E_asym);
}

void Analyzer::fillPionTrue()
{
    // Fill 4-Momentum
    if (truth_pi0_E != -1){
        pion.set_p4(truth_pi0_px,
                    truth_pi0_py,
                    truth_pi0_pz,
                    truth_pi0_E,
                    true);
    }else{
        pion.set_p4(truth_gamma_px[0]+truth_gamma_px[1] ,
                    truth_gamma_py[0]+truth_gamma_py[1],
                    truth_gamma_pz[0]+truth_gamma_pz[1],
                    truth_gamma_E[0]+truth_gamma_E[1],
                    true);
    }
    
    // set Angle wrt Beam
    pion.set_angleBeam(beam_p3, true);
    
    // set Angle wrt Muon
    pion.set_angleMuon(muon, true);
}

void Analyzer::fillMuonTrue()
{

    // Fill 4-Momentum
    muon.set_p4(    truth_muon_px * HEP_Functions::MeV_to_GeV,
                    truth_muon_py * HEP_Functions::MeV_to_GeV,
                    truth_muon_pz * HEP_Functions::MeV_to_GeV,
                    truth_muon_E, 
                    true);
       
    // set Angle wrt Beam
    muon.set_angleBeam(beam_p3, true);
    
    // set Angle wrt Muon
    muon.set_angleMuon(muon, true);
    
}

void Analyzer::fillMuonReco()
{
    // Set Particle Score
    muon.particleScore = CCProtonPi0_muon_muScore;
    
    // Fill 4-Momentum
    muon.set_p4(    CCProtonPi0_muon_px * HEP_Functions::MeV_to_GeV,
                    CCProtonPi0_muon_py * HEP_Functions::MeV_to_GeV,
                    CCProtonPi0_muon_pz * HEP_Functions::MeV_to_GeV,
                    CCProtonPi0_muon_E * HEP_Functions::MeV_to_GeV,
                    false);
    
    // set Angle wrt Beam
    muon.set_angleBeam(beam_p3, false);
    
    // set Angle wrt Muon
    muon.set_angleMuon(muon, false);
}


void Analyzer::Init(string playlist, TChain* fChain)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

    ifstream input_pl(playlist.c_str());
    string filename;
    
    if( !input_pl.is_open() ){
        cerr<<"Cannot open Playlist File!"<<endl;
        exit(1);
    }else{
        cout<<"Reading Playlist: "<<playlist.c_str()<<endl;
    }
    
    
   while (input_pl) {
     input_pl>>filename;
     
     if (!input_pl) break;
    
     if (filename[0] != '/') break;
    
     fChain->Add( filename.c_str() );
//      cout<<" Added "<<filename.c_str()<<endl;
   }

   // Set branch addresses and branch pointers
   fChain->SetMakeClass(1);
   
   fChain->SetBranchAddress("eventID", &eventID, &b_eventID);
   fChain->SetBranchAddress("physEvtNum", &physEvtNum, &b_physEvtNum);
   fChain->SetBranchAddress("n_hyps", &n_hyps, &b_n_hyps);
   fChain->SetBranchAddress("processType", &processType, &b_processType);
   fChain->SetBranchAddress("primaryPart", &primaryPart, &b_primaryPart);
   fChain->SetBranchAddress("n_slices", &n_slices, &b_n_slices);
   fChain->SetBranchAddress("slice_numbers", slice_numbers, &b_slice_numbers);
   fChain->SetBranchAddress("shared_slice", &shared_slice, &b_shared_slice);
   fChain->SetBranchAddress("vtx", vtx, &b_vtx);
   fChain->SetBranchAddress("vtxErr", vtxErr, &b_vtxErr);
   fChain->SetBranchAddress("E", E, &b_E);
   fChain->SetBranchAddress("found_truth", &found_truth, &b_found_truth);
   fChain->SetBranchAddress("phys_front_activity", &phys_front_activity, &b_phys_front_activity);
   fChain->SetBranchAddress("phys_energy_in_road_upstream_is_rockmuon_consistent", &phys_energy_in_road_upstream_is_rockmuon_consistent, &b_phys_energy_in_road_upstream_is_rockmuon_consistent);
   fChain->SetBranchAddress("rock_muons_removed", &rock_muons_removed, &b_rock_muons_removed);
   fChain->SetBranchAddress("minos_track_match", &minos_track_match, &b_minos_track_match);
   fChain->SetBranchAddress("minos_stub_match", &minos_stub_match, &b_minos_stub_match);
   fChain->SetBranchAddress("unknown_helicity", &unknown_helicity, &b_unknown_helicity);
   fChain->SetBranchAddress("minos_track_inside_partial_plane", &minos_track_inside_partial_plane, &b_minos_track_inside_partial_plane);
   fChain->SetBranchAddress("prim_vtx_has_misassigned_track_direction", &prim_vtx_has_misassigned_track_direction, &b_prim_vtx_has_misassigned_track_direction);
   fChain->SetBranchAddress("prim_vtx_has_broken_track", &prim_vtx_has_broken_track, &b_prim_vtx_has_broken_track);
   fChain->SetBranchAddress("gamma1_isGoodDirection", &gamma1_isGoodDirection, &b_gamma1_isGoodDirection);
   fChain->SetBranchAddress("gamma1_isGoodPosition", &gamma1_isGoodPosition, &b_gamma1_isGoodPosition);
   fChain->SetBranchAddress("gamma1_isGoodBlob", &gamma1_isGoodBlob, &b_gamma1_isGoodBlob);
   fChain->SetBranchAddress("gamma2_isGoodDirection", &gamma2_isGoodDirection, &b_gamma2_isGoodDirection);
   fChain->SetBranchAddress("gamma2_isGoodPosition", &gamma2_isGoodPosition, &b_gamma2_isGoodPosition);
   fChain->SetBranchAddress("gamma2_isGoodBlob", &gamma2_isGoodBlob, &b_gamma2_isGoodBlob);
   fChain->SetBranchAddress("Cut_ConeBlobs", &Cut_ConeBlobs, &b_Cut_ConeBlobs);
   fChain->SetBranchAddress("Cut_EndPoint_Michel_Exist", &Cut_EndPoint_Michel_Exist, &b_Cut_EndPoint_Michel_Exist);
   fChain->SetBranchAddress("Cut_Muon_Charge", &Cut_Muon_Charge, &b_Cut_Muon_Charge);
   fChain->SetBranchAddress("Cut_Muon_None", &Cut_Muon_None, &b_Cut_Muon_None);
   fChain->SetBranchAddress("Cut_Muon_Not_Plausible", &Cut_Muon_Not_Plausible, &b_Cut_Muon_Not_Plausible);
   fChain->SetBranchAddress("Cut_Muon_Score_Low", &Cut_Muon_Score_Low, &b_Cut_Muon_Score_Low);
   fChain->SetBranchAddress("Cut_Particle_None", &Cut_Particle_None, &b_Cut_Particle_None);
   fChain->SetBranchAddress("Cut_PreFilter_Pi0", &Cut_PreFilter_Pi0, &b_Cut_PreFilter_Pi0);
   fChain->SetBranchAddress("Cut_Proton_None", &Cut_Proton_None, &b_Cut_Proton_None);
   fChain->SetBranchAddress("Cut_Vertex_Michel_Exist", &Cut_Vertex_Michel_Exist, &b_Cut_Vertex_Michel_Exist);
   fChain->SetBranchAddress("Cut_Vertex_None", &Cut_Vertex_None, &b_Cut_Vertex_None);
   fChain->SetBranchAddress("Cut_Vertex_Not_Fiducial", &Cut_Vertex_Not_Fiducial, &b_Cut_Vertex_Not_Fiducial);
   fChain->SetBranchAddress("Cut_Vertex_Not_Reconstructable", &Cut_Vertex_Not_Reconstructable, &b_Cut_Vertex_Not_Reconstructable);
   fChain->SetBranchAddress("Cut_VtxBlob", &Cut_VtxBlob, &b_Cut_VtxBlob);
   fChain->SetBranchAddress("Cut_nProngs", &Cut_nProngs, &b_Cut_nProngs);
   fChain->SetBranchAddress("Cut_secEndPoint_Michel_Exist", &Cut_secEndPoint_Michel_Exist, &b_Cut_secEndPoint_Michel_Exist);
   fChain->SetBranchAddress("anglescan_ncand", &anglescan_ncand, &b_anglescan_ncand);
   fChain->SetBranchAddress("anglescan_ncandx", &anglescan_ncandx, &b_anglescan_ncandx);
   fChain->SetBranchAddress("blob_ndof_1", &blob_ndof_1, &b_blob_ndof_1);
   fChain->SetBranchAddress("blob_ndof_2", &blob_ndof_2, &b_blob_ndof_2);
   fChain->SetBranchAddress("broken_track_most_us_plane", &broken_track_most_us_plane, &b_broken_track_most_us_plane);
   fChain->SetBranchAddress("g1dedx_doublet", &g1dedx_doublet, &b_g1dedx_doublet);
   fChain->SetBranchAddress("g1dedx_empty_plane", &g1dedx_empty_plane, &b_g1dedx_empty_plane);
   fChain->SetBranchAddress("g1dedx_nplane", &g1dedx_nplane, &b_g1dedx_nplane);
   fChain->SetBranchAddress("g2dedx_doublet", &g2dedx_doublet, &b_g2dedx_doublet);
   fChain->SetBranchAddress("g2dedx_empty_plane", &g2dedx_empty_plane, &b_g2dedx_empty_plane);
   fChain->SetBranchAddress("g2dedx_nplane", &g2dedx_nplane, &b_g2dedx_nplane);
   fChain->SetBranchAddress("gamma1_blob_nclusters", &gamma1_blob_nclusters, &b_gamma1_blob_nclusters);
   fChain->SetBranchAddress("gamma1_blob_ndigits", &gamma1_blob_ndigits, &b_gamma1_blob_ndigits);
   fChain->SetBranchAddress("gamma2_blob_nclusters", &gamma2_blob_nclusters, &b_gamma2_blob_nclusters);
   fChain->SetBranchAddress("gamma2_blob_ndigits", &gamma2_blob_ndigits, &b_gamma2_blob_ndigits);
   fChain->SetBranchAddress("nProngs", &nProngs, &b_nProngs);
   fChain->SetBranchAddress("n_anchored_long_trk_prongs", &n_anchored_long_trk_prongs, &b_n_anchored_long_trk_prongs);
   fChain->SetBranchAddress("n_anchored_short_trk_prongs", &n_anchored_short_trk_prongs, &b_n_anchored_short_trk_prongs);
   fChain->SetBranchAddress("n_iso_trk_prongs", &n_iso_trk_prongs, &b_n_iso_trk_prongs);
   fChain->SetBranchAddress("n_vtx_michel_views", &n_vtx_michel_views, &b_n_vtx_michel_views);
   fChain->SetBranchAddress("od_energeticTower", &od_energeticTower, &b_od_energeticTower);
   fChain->SetBranchAddress("phys_energy_in_road_downstream_nplanes", &phys_energy_in_road_downstream_nplanes, &b_phys_energy_in_road_downstream_nplanes);
   fChain->SetBranchAddress("phys_energy_in_road_upstream_nplanes", &phys_energy_in_road_upstream_nplanes, &b_phys_energy_in_road_upstream_nplanes);
   fChain->SetBranchAddress("phys_n_dead_discr_pair", &phys_n_dead_discr_pair, &b_phys_n_dead_discr_pair);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_in_prim_track_region", &phys_n_dead_discr_pair_in_prim_track_region, &b_phys_n_dead_discr_pair_in_prim_track_region);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_two_mod_downstream_prim_track", &phys_n_dead_discr_pair_two_mod_downstream_prim_track, &b_phys_n_dead_discr_pair_two_mod_downstream_prim_track);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_two_mod_upstream_prim_vtx", &phys_n_dead_discr_pair_two_mod_upstream_prim_vtx, &b_phys_n_dead_discr_pair_two_mod_upstream_prim_vtx);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_upstream_prim_track_proj", &phys_n_dead_discr_pair_upstream_prim_track_proj, &b_phys_n_dead_discr_pair_upstream_prim_track_proj);
   fChain->SetBranchAddress("phys_vertex_is_fiducial", &phys_vertex_is_fiducial, &b_phys_vertex_is_fiducial);
   fChain->SetBranchAddress("preFilter_Result", &preFilter_Result, &b_preFilter_Result);
   fChain->SetBranchAddress("vtx_primary_index", &vtx_primary_index, &b_vtx_primary_index);
   fChain->SetBranchAddress("vtx_primary_multiplicity", &vtx_primary_multiplicity, &b_vtx_primary_multiplicity);
   fChain->SetBranchAddress("vtx_secondary_count", &vtx_secondary_count, &b_vtx_secondary_count);
   fChain->SetBranchAddress("vtx_total_count", &vtx_total_count, &b_vtx_total_count);
   fChain->SetBranchAddress("Filament_Vertex_energy", &Filament_Vertex_energy, &b_Filament_Vertex_energy);
   fChain->SetBranchAddress("RE_energy_ECAL", &RE_energy_ECAL, &b_RE_energy_ECAL);
   fChain->SetBranchAddress("RE_energy_HCAL", &RE_energy_HCAL, &b_RE_energy_HCAL);
   fChain->SetBranchAddress("RE_energy_Tracker", &RE_energy_Tracker, &b_RE_energy_Tracker);
   fChain->SetBranchAddress("Sphere_Vertex_energy", &Sphere_Vertex_energy, &b_Sphere_Vertex_energy);
   fChain->SetBranchAddress("Vertex_blob_energy", &Vertex_blob_energy, &b_Vertex_blob_energy);
   fChain->SetBranchAddress("blob_fval_1", &blob_fval_1, &b_blob_fval_1);
   fChain->SetBranchAddress("blob_fval_2", &blob_fval_2, &b_blob_fval_2);
   fChain->SetBranchAddress("endpoint_michel_distance", &endpoint_michel_distance, &b_endpoint_michel_distance);
   fChain->SetBranchAddress("energyUnused_afterReco", &energyUnused_afterReco, &b_energyUnused_afterReco);
   fChain->SetBranchAddress("energyUsed_afterReco", &energyUsed_afterReco, &b_energyUsed_afterReco);
   fChain->SetBranchAddress("energy_from_mc", &energy_from_mc, &b_energy_from_mc);
   fChain->SetBranchAddress("energy_from_mc_fraction", &energy_from_mc_fraction, &b_energy_from_mc_fraction);
   fChain->SetBranchAddress("energy_from_mc_fraction_of_highest", &energy_from_mc_fraction_of_highest, &b_energy_from_mc_fraction_of_highest);
   fChain->SetBranchAddress("evis_ECAL", &evis_ECAL, &b_evis_ECAL);
   fChain->SetBranchAddress("evis_HCAL", &evis_HCAL, &b_evis_HCAL);
   fChain->SetBranchAddress("evis_NuclearTarget", &evis_NuclearTarget, &b_evis_NuclearTarget);
   fChain->SetBranchAddress("evis_TotalExceptNuclearTarget", &evis_TotalExceptNuclearTarget, &b_evis_TotalExceptNuclearTarget);
   fChain->SetBranchAddress("evis_Tracker", &evis_Tracker, &b_evis_Tracker);
   fChain->SetBranchAddress("evis_nearvtx", &evis_nearvtx, &b_evis_nearvtx);
   fChain->SetBranchAddress("evis_total", &evis_total, &b_evis_total);
   fChain->SetBranchAddress("g1blob_minsep", &g1blob_minsep, &b_g1blob_minsep);
   fChain->SetBranchAddress("g1dedx", &g1dedx, &b_g1dedx);
   fChain->SetBranchAddress("g1dedx1", &g1dedx1, &b_g1dedx1);
   fChain->SetBranchAddress("g1dedx_total", &g1dedx_total, &b_g1dedx_total);
   fChain->SetBranchAddress("g1dedx_total1", &g1dedx_total1, &b_g1dedx_total1);
   fChain->SetBranchAddress("g2blob_minsep", &g2blob_minsep, &b_g2blob_minsep);
   fChain->SetBranchAddress("g2dedx", &g2dedx, &b_g2dedx);
   fChain->SetBranchAddress("g2dedx1", &g2dedx1, &b_g2dedx1);
   fChain->SetBranchAddress("g2dedx_total", &g2dedx_total, &b_g2dedx_total);
   fChain->SetBranchAddress("g2dedx_total1", &g2dedx_total1, &b_g2dedx_total1);
   fChain->SetBranchAddress("gamma1_E", &gamma1_E, &b_gamma1_E);
   fChain->SetBranchAddress("gamma1_dEdx", &gamma1_dEdx, &b_gamma1_dEdx);
   fChain->SetBranchAddress("gamma1_dist_vtx", &gamma1_dist_vtx, &b_gamma1_dist_vtx);
   fChain->SetBranchAddress("gamma1_evis_ecal", &gamma1_evis_ecal, &b_gamma1_evis_ecal);
   fChain->SetBranchAddress("gamma1_evis_hcal", &gamma1_evis_hcal, &b_gamma1_evis_hcal);
   fChain->SetBranchAddress("gamma1_evis_scal", &gamma1_evis_scal, &b_gamma1_evis_scal);
   fChain->SetBranchAddress("gamma1_evis_trkr", &gamma1_evis_trkr, &b_gamma1_evis_trkr);
   fChain->SetBranchAddress("gamma1_phi", &gamma1_phi, &b_gamma1_phi);
   fChain->SetBranchAddress("gamma1_px", &gamma1_px, &b_gamma1_px);
   fChain->SetBranchAddress("gamma1_py", &gamma1_py, &b_gamma1_py);
   fChain->SetBranchAddress("gamma1_pz", &gamma1_pz, &b_gamma1_pz);
   fChain->SetBranchAddress("gamma1_score", &gamma1_score, &b_gamma1_score);
   fChain->SetBranchAddress("gamma1_theta", &gamma1_theta, &b_gamma1_theta);
   fChain->SetBranchAddress("gamma1_time", &gamma1_time, &b_gamma1_time);
   fChain->SetBranchAddress("gamma2_E", &gamma2_E, &b_gamma2_E);
   fChain->SetBranchAddress("gamma2_dEdx", &gamma2_dEdx, &b_gamma2_dEdx);
   fChain->SetBranchAddress("gamma2_dist_vtx", &gamma2_dist_vtx, &b_gamma2_dist_vtx);
   fChain->SetBranchAddress("gamma2_evis_ecal", &gamma2_evis_ecal, &b_gamma2_evis_ecal);
   fChain->SetBranchAddress("gamma2_evis_hcal", &gamma2_evis_hcal, &b_gamma2_evis_hcal);
   fChain->SetBranchAddress("gamma2_evis_scal", &gamma2_evis_scal, &b_gamma2_evis_scal);
   fChain->SetBranchAddress("gamma2_evis_trkr", &gamma2_evis_trkr, &b_gamma2_evis_trkr);
   fChain->SetBranchAddress("gamma2_phi", &gamma2_phi, &b_gamma2_phi);
   fChain->SetBranchAddress("gamma2_px", &gamma2_px, &b_gamma2_px);
   fChain->SetBranchAddress("gamma2_py", &gamma2_py, &b_gamma2_py);
   fChain->SetBranchAddress("gamma2_pz", &gamma2_pz, &b_gamma2_pz);
   fChain->SetBranchAddress("gamma2_score", &gamma2_score, &b_gamma2_score);
   fChain->SetBranchAddress("gamma2_theta", &gamma2_theta, &b_gamma2_theta);
   fChain->SetBranchAddress("gamma2_time", &gamma2_time, &b_gamma2_time);
   fChain->SetBranchAddress("hadronVisibleE", &hadronVisibleE, &b_hadronVisibleE);
   fChain->SetBranchAddress("michelProng_begin_Z", &michelProng_begin_Z, &b_michelProng_begin_Z);
   fChain->SetBranchAddress("michelProng_end_Z", &michelProng_end_Z, &b_michelProng_end_Z);
   fChain->SetBranchAddress("michelProng_energy", &michelProng_energy, &b_michelProng_energy);
   fChain->SetBranchAddress("michelProng_time_diff", &michelProng_time_diff, &b_michelProng_time_diff);
   fChain->SetBranchAddress("muonVisibleE", &muonVisibleE, &b_muonVisibleE);
   fChain->SetBranchAddress("muon_phi", &muon_phi, &b_muon_phi);
   fChain->SetBranchAddress("muon_theta", &muon_theta, &b_muon_theta);
   fChain->SetBranchAddress("muon_thetaX", &muon_thetaX, &b_muon_thetaX);
   fChain->SetBranchAddress("muon_thetaY", &muon_thetaY, &b_muon_thetaY);
   fChain->SetBranchAddress("od_downstreamFrame", &od_downstreamFrame, &b_od_downstreamFrame);
   fChain->SetBranchAddress("od_downstreamFrame_z", &od_downstreamFrame_z, &b_od_downstreamFrame_z);
   fChain->SetBranchAddress("od_highStory", &od_highStory, &b_od_highStory);
   fChain->SetBranchAddress("od_highStory_t", &od_highStory_t, &b_od_highStory_t);
   fChain->SetBranchAddress("od_lowStory", &od_lowStory, &b_od_lowStory);
   fChain->SetBranchAddress("od_lowStory_t", &od_lowStory_t, &b_od_lowStory_t);
   fChain->SetBranchAddress("od_maxEnergy", &od_maxEnergy, &b_od_maxEnergy);
   fChain->SetBranchAddress("od_upstreamFrame", &od_upstreamFrame, &b_od_upstreamFrame);
   fChain->SetBranchAddress("od_upstreamFrame_z", &od_upstreamFrame_z, &b_od_upstreamFrame_z);
   fChain->SetBranchAddress("phys_energy_dispersed", &phys_energy_dispersed, &b_phys_energy_dispersed);
   fChain->SetBranchAddress("phys_energy_in_road_downstream", &phys_energy_in_road_downstream, &b_phys_energy_in_road_downstream);
   fChain->SetBranchAddress("phys_energy_in_road_upstream", &phys_energy_in_road_upstream, &b_phys_energy_in_road_upstream);
   fChain->SetBranchAddress("phys_energy_unattached", &phys_energy_unattached, &b_phys_energy_unattached);
   fChain->SetBranchAddress("pi0_E", &pi0_E, &b_pi0_E);
   fChain->SetBranchAddress("pi0_cos_openingAngle", &pi0_cos_openingAngle, &b_pi0_cos_openingAngle);
   fChain->SetBranchAddress("pi0_invMass", &pi0_invMass, &b_pi0_invMass);
   fChain->SetBranchAddress("pi0_openingAngle", &pi0_openingAngle, &b_pi0_openingAngle);
   fChain->SetBranchAddress("pi0_phi", &pi0_phi, &b_pi0_phi);
   fChain->SetBranchAddress("pi0_px", &pi0_px, &b_pi0_px);
   fChain->SetBranchAddress("pi0_py", &pi0_py, &b_pi0_py);
   fChain->SetBranchAddress("pi0_pz", &pi0_pz, &b_pi0_pz);
   fChain->SetBranchAddress("pi0_theta", &pi0_theta, &b_pi0_theta);
   fChain->SetBranchAddress("pi0_thetaX", &pi0_thetaX, &b_pi0_thetaX);
   fChain->SetBranchAddress("pi0_thetaY", &pi0_thetaY, &b_pi0_thetaY);
   fChain->SetBranchAddress("preFilter_rejectedEnergy", &preFilter_rejectedEnergy, &b_preFilter_rejectedEnergy);
   fChain->SetBranchAddress("prim_vtx_smallest_opening_angle", &prim_vtx_smallest_opening_angle, &b_prim_vtx_smallest_opening_angle);
   fChain->SetBranchAddress("prong_showerScore", &prong_showerScore, &b_prong_showerScore);
   fChain->SetBranchAddress("reco_eventID", &reco_eventID, &b_reco_eventID);
   fChain->SetBranchAddress("time", &time, &b_time);
   fChain->SetBranchAddress("totalIDVisibleE", &totalIDVisibleE, &b_totalIDVisibleE);
   fChain->SetBranchAddress("totalODVisibleE", &totalODVisibleE, &b_totalODVisibleE);
   fChain->SetBranchAddress("totalVisibleE", &totalVisibleE, &b_totalVisibleE);
   fChain->SetBranchAddress("vtx_michel_distance", &vtx_michel_distance, &b_vtx_michel_distance);
   fChain->SetBranchAddress("g1dedx_cluster_occupancy_sz", &g1dedx_cluster_occupancy_sz, &b_g1dedx_cluster_occupancy_sz);
   fChain->SetBranchAddress("g1dedx_cluster_occupancy", g1dedx_cluster_occupancy, &b_g1dedx_cluster_occupancy);
   fChain->SetBranchAddress("g2dedx_cluster_occupancy_sz", &g2dedx_cluster_occupancy_sz, &b_g2dedx_cluster_occupancy_sz);
   fChain->SetBranchAddress("g2dedx_cluster_occupancy", g2dedx_cluster_occupancy, &b_g2dedx_cluster_occupancy);
   fChain->SetBranchAddress("Vertex_energy_radii_sz", &Vertex_energy_radii_sz, &b_Vertex_energy_radii_sz);
   fChain->SetBranchAddress("Vertex_energy_radii", Vertex_energy_radii, &b_Vertex_energy_radii);
   fChain->SetBranchAddress("blob_cluster_energy1_sz", &blob_cluster_energy1_sz, &b_blob_cluster_energy1_sz);
   fChain->SetBranchAddress("blob_cluster_energy1", blob_cluster_energy1, &b_blob_cluster_energy1);
   fChain->SetBranchAddress("blob_cluster_energy2_sz", &blob_cluster_energy2_sz, &b_blob_cluster_energy2_sz);
   fChain->SetBranchAddress("blob_cluster_energy2", blob_cluster_energy2, &b_blob_cluster_energy2);
   fChain->SetBranchAddress("g1dedx_cluster_energy_sz", &g1dedx_cluster_energy_sz, &b_g1dedx_cluster_energy_sz);
   fChain->SetBranchAddress("g1dedx_cluster_energy", g1dedx_cluster_energy, &b_g1dedx_cluster_energy);
   fChain->SetBranchAddress("g1dedx_rev_cluster_energy_sz", &g1dedx_rev_cluster_energy_sz, &b_g1dedx_rev_cluster_energy_sz);
   fChain->SetBranchAddress("g1dedx_rev_cluster_energy", g1dedx_rev_cluster_energy, &b_g1dedx_rev_cluster_energy);
   fChain->SetBranchAddress("g2dedx_cluster_energy_sz", &g2dedx_cluster_energy_sz, &b_g2dedx_cluster_energy_sz);
   fChain->SetBranchAddress("g2dedx_cluster_energy", g2dedx_cluster_energy, &b_g2dedx_cluster_energy);
   fChain->SetBranchAddress("g2dedx_rev_cluster_energy_sz", &g2dedx_rev_cluster_energy_sz, &b_g2dedx_rev_cluster_energy_sz);
   fChain->SetBranchAddress("g2dedx_rev_cluster_energy", g2dedx_rev_cluster_energy, &b_g2dedx_rev_cluster_energy);
   fChain->SetBranchAddress("gamma1_direction", gamma1_direction, &b_gamma1_direction);
   fChain->SetBranchAddress("gamma1_vertex", gamma1_vertex, &b_gamma1_vertex);
   fChain->SetBranchAddress("gamma2_direction", gamma2_direction, &b_gamma2_direction);
   fChain->SetBranchAddress("gamma2_vertex", gamma2_vertex, &b_gamma2_vertex);
   fChain->SetBranchAddress("od_distanceBlobTower_sz", &od_distanceBlobTower_sz, &b_od_distanceBlobTower_sz);
   fChain->SetBranchAddress("od_distanceBlobTower", od_distanceBlobTower, &b_od_distanceBlobTower);
   fChain->SetBranchAddress("od_idBlobTime_sz", &od_idBlobTime_sz, &b_od_idBlobTime_sz);
   fChain->SetBranchAddress("od_idBlobTime", od_idBlobTime, &b_od_idBlobTime);
   fChain->SetBranchAddress("od_towerEnergy_sz", &od_towerEnergy_sz, &b_od_towerEnergy_sz);
   fChain->SetBranchAddress("od_towerEnergy", od_towerEnergy, &b_od_towerEnergy);
   fChain->SetBranchAddress("od_towerNClusters_sz", &od_towerNClusters_sz, &b_od_towerNClusters_sz);
   fChain->SetBranchAddress("od_towerNClusters", od_towerNClusters, &b_od_towerNClusters);
   fChain->SetBranchAddress("od_towerTime_sz", &od_towerTime_sz, &b_od_towerTime_sz);
   fChain->SetBranchAddress("od_towerTime", od_towerTime, &b_od_towerTime);
   fChain->SetBranchAddress("od_towerTimeBlobMuon_sz", &od_towerTimeBlobMuon_sz, &b_od_towerTimeBlobMuon_sz);
   fChain->SetBranchAddress("od_towerTimeBlobMuon", od_towerTimeBlobMuon, &b_od_towerTimeBlobMuon);
   fChain->SetBranchAddress("od_towerTimeBlobOD_sz", &od_towerTimeBlobOD_sz, &b_od_towerTimeBlobOD_sz);
   fChain->SetBranchAddress("od_towerTimeBlobOD", od_towerTimeBlobOD, &b_od_towerTimeBlobOD);
   fChain->SetBranchAddress("truth_has_physics_event", &truth_has_physics_event, &b_truth_has_physics_event);
   fChain->SetBranchAddress("truth_isSignal", &truth_isSignal, &b_truth_isSignal);
   fChain->SetBranchAddress("truth_isSignal_1Pi0", &truth_isSignal_1Pi0, &b_truth_isSignal_1Pi0);
   fChain->SetBranchAddress("truth_isSignal_2Gamma", &truth_isSignal_2Gamma, &b_truth_isSignal_2Gamma);
   fChain->SetBranchAddress("truth_isFidVol", &truth_isFidVol, &b_truth_isFidVol);
   fChain->SetBranchAddress("truth_AnalyzeEvent", &truth_AnalyzeEvent, &b_truth_AnalyzeEvent);
   fChain->SetBranchAddress("truth_isBckg_QELike", &truth_isBckg_QELike, &b_truth_isBckg_QELike);
   fChain->SetBranchAddress("truth_isBckg_SinglePiPlus", &truth_isBckg_SinglePiPlus, &b_truth_isBckg_SinglePiPlus);
   fChain->SetBranchAddress("truth_isBckg_SinglePiMinus", &truth_isBckg_SinglePiMinus, &b_truth_isBckg_SinglePiMinus);
   fChain->SetBranchAddress("truth_isBckg_MultiPion", &truth_isBckg_MultiPion, &b_truth_isBckg_MultiPion);
   fChain->SetBranchAddress("truth_isBckg_MultiPiZero", &truth_isBckg_MultiPiZero, &b_truth_isBckg_MultiPiZero);
   fChain->SetBranchAddress("truth_isBckg_Other", &truth_isBckg_Other, &b_truth_isBckg_Other);
   fChain->SetBranchAddress("truth_isBckg_withAntiMuon", &truth_isBckg_withAntiMuon, &b_truth_isBckg_withAntiMuon);
   fChain->SetBranchAddress("truth_isBckg_withMichel", &truth_isBckg_withMichel, &b_truth_isBckg_withMichel);
   fChain->SetBranchAddress("truth_isBckg_withPrimaryPi0", &truth_isBckg_withPrimaryPi0, &b_truth_isBckg_withPrimaryPi0);
   fChain->SetBranchAddress("truth_isBckg_withSecondaryPi0", &truth_isBckg_withSecondaryPi0, &b_truth_isBckg_withSecondaryPi0);
   fChain->SetBranchAddress("truth_N_FSParticles", &truth_N_FSParticles, &b_truth_N_FSParticles);
   fChain->SetBranchAddress("truth_N_gamma", &truth_N_gamma, &b_truth_N_gamma);
   fChain->SetBranchAddress("truth_N_pi0", &truth_N_pi0, &b_truth_N_pi0);
   fChain->SetBranchAddress("truth_N_proton", &truth_N_proton, &b_truth_N_proton);
   fChain->SetBranchAddress("truth_N_trueMichelElectrons", &truth_N_trueMichelElectrons, &b_truth_N_trueMichelElectrons);
   fChain->SetBranchAddress("truth_muon_charge", &truth_muon_charge, &b_truth_muon_charge);
   fChain->SetBranchAddress("truth_pi0_GrandMother", &truth_pi0_GrandMother, &b_truth_pi0_GrandMother);
   fChain->SetBranchAddress("truth_pi0_GrandMotherStatus", &truth_pi0_GrandMotherStatus, &b_truth_pi0_GrandMotherStatus);
   fChain->SetBranchAddress("truth_pi0_Mother", &truth_pi0_Mother, &b_truth_pi0_Mother);
   fChain->SetBranchAddress("truth_pi0_MotherStatus", &truth_pi0_MotherStatus, &b_truth_pi0_MotherStatus);
   fChain->SetBranchAddress("truth_pi0_status", &truth_pi0_status, &b_truth_pi0_status);
   fChain->SetBranchAddress("truth_target_material", &truth_target_material, &b_truth_target_material);
   fChain->SetBranchAddress("truth_vertex_module", &truth_vertex_module, &b_truth_vertex_module);
   fChain->SetBranchAddress("truth_vertex_plane", &truth_vertex_plane, &b_truth_vertex_plane);
   fChain->SetBranchAddress("truth_eventID", &truth_eventID, &b_truth_eventID);
   fChain->SetBranchAddress("truth_michelElectron_E", &truth_michelElectron_E, &b_truth_michelElectron_E);
   fChain->SetBranchAddress("truth_michelElectron_P", &truth_michelElectron_P, &b_truth_michelElectron_P);
   fChain->SetBranchAddress("truth_michelMuon_P", &truth_michelMuon_P, &b_truth_michelMuon_P);
   fChain->SetBranchAddress("truth_michelMuon_end_dist_vtx", &truth_michelMuon_end_dist_vtx, &b_truth_michelMuon_end_dist_vtx);
   fChain->SetBranchAddress("truth_michelMuon_length", &truth_michelMuon_length, &b_truth_michelMuon_length);
   fChain->SetBranchAddress("truth_michelPion_P", &truth_michelPion_P, &b_truth_michelPion_P);
   fChain->SetBranchAddress("truth_michelPion_begin_dist_vtx", &truth_michelPion_begin_dist_vtx, &b_truth_michelPion_begin_dist_vtx);
   fChain->SetBranchAddress("truth_michelPion_length", &truth_michelPion_length, &b_truth_michelPion_length);
   fChain->SetBranchAddress("truth_muon_E", &truth_muon_E, &b_truth_muon_E);
   fChain->SetBranchAddress("truth_muon_px", &truth_muon_px, &b_truth_muon_px);
   fChain->SetBranchAddress("truth_muon_py", &truth_muon_py, &b_truth_muon_py);
   fChain->SetBranchAddress("truth_muon_pz", &truth_muon_pz, &b_truth_muon_pz);
   fChain->SetBranchAddress("truth_muon_theta_wrtbeam", &truth_muon_theta_wrtbeam, &b_truth_muon_theta_wrtbeam);
   fChain->SetBranchAddress("truth_pi0_E", &truth_pi0_E, &b_truth_pi0_E);
   fChain->SetBranchAddress("truth_pi0_px", &truth_pi0_px, &b_truth_pi0_px);
   fChain->SetBranchAddress("truth_pi0_py", &truth_pi0_py, &b_truth_pi0_py);
   fChain->SetBranchAddress("truth_pi0_pz", &truth_pi0_pz, &b_truth_pi0_pz);
   fChain->SetBranchAddress("truth_pi0_theta_wrtbeam", &truth_pi0_theta_wrtbeam, &b_truth_pi0_theta_wrtbeam);
   fChain->SetBranchAddress("truth_gamma_E", truth_gamma_E, &b_truth_gamma_E);
   fChain->SetBranchAddress("truth_gamma_px", truth_gamma_px, &b_truth_gamma_px);
   fChain->SetBranchAddress("truth_gamma_py", truth_gamma_py, &b_truth_gamma_py);
   fChain->SetBranchAddress("truth_gamma_pz", truth_gamma_pz, &b_truth_gamma_pz);
   fChain->SetBranchAddress("truth_gamma_theta_wrtbeam", truth_gamma_theta_wrtbeam, &b_truth_gamma_theta_wrtbeam);
   fChain->SetBranchAddress("genie_wgt_n_shifts", &genie_wgt_n_shifts, &b_genie_wgt_n_shifts);
   fChain->SetBranchAddress("truth_genie_wgt_AGKYxF1pi", truth_genie_wgt_AGKYxF1pi, &b_truth_genie_wgt_AGKYxF1pi);
   fChain->SetBranchAddress("truth_genie_wgt_AhtBY", truth_genie_wgt_AhtBY, &b_truth_genie_wgt_AhtBY);
   fChain->SetBranchAddress("truth_genie_wgt_BhtBY", truth_genie_wgt_BhtBY, &b_truth_genie_wgt_BhtBY);
   fChain->SetBranchAddress("truth_genie_wgt_CCQEPauliSupViaKF", truth_genie_wgt_CCQEPauliSupViaKF, &b_truth_genie_wgt_CCQEPauliSupViaKF);
   fChain->SetBranchAddress("truth_genie_wgt_CV1uBY", truth_genie_wgt_CV1uBY, &b_truth_genie_wgt_CV1uBY);
   fChain->SetBranchAddress("truth_genie_wgt_CV2uBY", truth_genie_wgt_CV2uBY, &b_truth_genie_wgt_CV2uBY);
   fChain->SetBranchAddress("truth_genie_wgt_EtaNCEL", truth_genie_wgt_EtaNCEL, &b_truth_genie_wgt_EtaNCEL);
   fChain->SetBranchAddress("truth_genie_wgt_FrAbs_N", truth_genie_wgt_FrAbs_N, &b_truth_genie_wgt_FrAbs_N);
   fChain->SetBranchAddress("truth_genie_wgt_FrAbs_pi", truth_genie_wgt_FrAbs_pi, &b_truth_genie_wgt_FrAbs_pi);
   fChain->SetBranchAddress("truth_genie_wgt_FrCEx_N", truth_genie_wgt_FrCEx_N, &b_truth_genie_wgt_FrCEx_N);
   fChain->SetBranchAddress("truth_genie_wgt_FrCEx_pi", truth_genie_wgt_FrCEx_pi, &b_truth_genie_wgt_FrCEx_pi);
   fChain->SetBranchAddress("truth_genie_wgt_FrElas_N", truth_genie_wgt_FrElas_N, &b_truth_genie_wgt_FrElas_N);
   fChain->SetBranchAddress("truth_genie_wgt_FrElas_pi", truth_genie_wgt_FrElas_pi, &b_truth_genie_wgt_FrElas_pi);
   fChain->SetBranchAddress("truth_genie_wgt_FrInel_N", truth_genie_wgt_FrInel_N, &b_truth_genie_wgt_FrInel_N);
   fChain->SetBranchAddress("truth_genie_wgt_FrInel_pi", truth_genie_wgt_FrInel_pi, &b_truth_genie_wgt_FrInel_pi);
   fChain->SetBranchAddress("truth_genie_wgt_FrPiProd_N", truth_genie_wgt_FrPiProd_N, &b_truth_genie_wgt_FrPiProd_N);
   fChain->SetBranchAddress("truth_genie_wgt_FrPiProd_pi", truth_genie_wgt_FrPiProd_pi, &b_truth_genie_wgt_FrPiProd_pi);
   fChain->SetBranchAddress("truth_genie_wgt_MFP_N", truth_genie_wgt_MFP_N, &b_truth_genie_wgt_MFP_N);
   fChain->SetBranchAddress("truth_genie_wgt_MFP_pi", truth_genie_wgt_MFP_pi, &b_truth_genie_wgt_MFP_pi);
   fChain->SetBranchAddress("truth_genie_wgt_MaCCQE", truth_genie_wgt_MaCCQE, &b_truth_genie_wgt_MaCCQE);
   fChain->SetBranchAddress("truth_genie_wgt_MaCCQEshape", truth_genie_wgt_MaCCQEshape, &b_truth_genie_wgt_MaCCQEshape);
   fChain->SetBranchAddress("truth_genie_wgt_MaNCEL", truth_genie_wgt_MaNCEL, &b_truth_genie_wgt_MaNCEL);
   fChain->SetBranchAddress("truth_genie_wgt_MaRES", truth_genie_wgt_MaRES, &b_truth_genie_wgt_MaRES);
   fChain->SetBranchAddress("truth_genie_wgt_MvRES", truth_genie_wgt_MvRES, &b_truth_genie_wgt_MvRES);
   fChain->SetBranchAddress("truth_genie_wgt_NormCCQE", truth_genie_wgt_NormCCQE, &b_truth_genie_wgt_NormCCQE);
   fChain->SetBranchAddress("truth_genie_wgt_NormCCRES", truth_genie_wgt_NormCCRES, &b_truth_genie_wgt_NormCCRES);
   fChain->SetBranchAddress("truth_genie_wgt_NormDISCC", truth_genie_wgt_NormDISCC, &b_truth_genie_wgt_NormDISCC);
   fChain->SetBranchAddress("truth_genie_wgt_NormNCRES", truth_genie_wgt_NormNCRES, &b_truth_genie_wgt_NormNCRES);
   fChain->SetBranchAddress("truth_genie_wgt_RDecBR1gamma", truth_genie_wgt_RDecBR1gamma, &b_truth_genie_wgt_RDecBR1gamma);
   fChain->SetBranchAddress("truth_genie_wgt_Rvn1pi", truth_genie_wgt_Rvn1pi, &b_truth_genie_wgt_Rvn1pi);
   fChain->SetBranchAddress("truth_genie_wgt_Rvn2pi", truth_genie_wgt_Rvn2pi, &b_truth_genie_wgt_Rvn2pi);
   fChain->SetBranchAddress("truth_genie_wgt_Rvp1pi", truth_genie_wgt_Rvp1pi, &b_truth_genie_wgt_Rvp1pi);
   fChain->SetBranchAddress("truth_genie_wgt_Rvp2pi", truth_genie_wgt_Rvp2pi, &b_truth_genie_wgt_Rvp2pi);
   fChain->SetBranchAddress("truth_genie_wgt_Theta_Delta2Npi", truth_genie_wgt_Theta_Delta2Npi, &b_truth_genie_wgt_Theta_Delta2Npi);
   fChain->SetBranchAddress("truth_genie_wgt_VecFFCCQEshape", truth_genie_wgt_VecFFCCQEshape, &b_truth_genie_wgt_VecFFCCQEshape);
   fChain->SetBranchAddress("truth_genie_wgt_shifts", truth_genie_wgt_shifts, &b_truth_genie_wgt_shifts);
   fChain->SetBranchAddress("truth_michelMuon_endPoint", truth_michelMuon_endPoint, &b_truth_michelMuon_endPoint);
   fChain->SetBranchAddress("truth_proton_E", truth_proton_E, &b_truth_proton_E);
   fChain->SetBranchAddress("truth_proton_px", truth_proton_px, &b_truth_proton_px);
   fChain->SetBranchAddress("truth_proton_py", truth_proton_py, &b_truth_proton_py);
   fChain->SetBranchAddress("truth_proton_pz", truth_proton_pz, &b_truth_proton_pz);
   fChain->SetBranchAddress("truth_proton_theta_wrtbeam", truth_proton_theta_wrtbeam, &b_truth_proton_theta_wrtbeam);
   fChain->SetBranchAddress("CCProtonPi0_nuFlavor", &CCProtonPi0_nuFlavor, &b_CCProtonPi0_nuFlavor);
   fChain->SetBranchAddress("CCProtonPi0_nuHelicity", &CCProtonPi0_nuHelicity, &b_CCProtonPi0_nuHelicity);
   fChain->SetBranchAddress("CCProtonPi0_intCurrent", &CCProtonPi0_intCurrent, &b_CCProtonPi0_intCurrent);
   fChain->SetBranchAddress("CCProtonPi0_intType", &CCProtonPi0_intType, &b_CCProtonPi0_intType);
   fChain->SetBranchAddress("CCProtonPi0_E", &CCProtonPi0_E, &b_CCProtonPi0_E);
   fChain->SetBranchAddress("CCProtonPi0_Q2", &CCProtonPi0_Q2, &b_CCProtonPi0_Q2);
   fChain->SetBranchAddress("CCProtonPi0_x", &CCProtonPi0_x, &b_CCProtonPi0_x);
   fChain->SetBranchAddress("CCProtonPi0_y", &CCProtonPi0_y, &b_CCProtonPi0_y);
   fChain->SetBranchAddress("CCProtonPi0_W", &CCProtonPi0_W, &b_CCProtonPi0_W);
   fChain->SetBranchAddress("CCProtonPi0_score", &CCProtonPi0_score, &b_CCProtonPi0_score);
   fChain->SetBranchAddress("CCProtonPi0_leptonE", CCProtonPi0_leptonE, &b_CCProtonPi0_leptonE);
   fChain->SetBranchAddress("CCProtonPi0_vtx", CCProtonPi0_vtx, &b_CCProtonPi0_vtx);
   fChain->SetBranchAddress("CCProtonPi0_minos_trk_is_contained", &CCProtonPi0_minos_trk_is_contained, &b_CCProtonPi0_minos_trk_is_contained);
   fChain->SetBranchAddress("CCProtonPi0_minos_trk_is_ok", &CCProtonPi0_minos_trk_is_ok, &b_CCProtonPi0_minos_trk_is_ok);
   fChain->SetBranchAddress("CCProtonPi0_minos_used_range", &CCProtonPi0_minos_used_range, &b_CCProtonPi0_minos_used_range);
   fChain->SetBranchAddress("CCProtonPi0_minos_used_curvature", &CCProtonPi0_minos_used_curvature, &b_CCProtonPi0_minos_used_curvature);
   fChain->SetBranchAddress("CCProtonPi0_isMuonInsideOD", &CCProtonPi0_isMuonInsideOD, &b_CCProtonPi0_isMuonInsideOD);
   fChain->SetBranchAddress("CCProtonPi0_minos_trk_end_plane", &CCProtonPi0_minos_trk_end_plane, &b_CCProtonPi0_minos_trk_end_plane);
   fChain->SetBranchAddress("CCProtonPi0_minos_trk_quality", &CCProtonPi0_minos_trk_quality, &b_CCProtonPi0_minos_trk_quality);
   fChain->SetBranchAddress("CCProtonPi0_muon_N_minosTracks", &CCProtonPi0_muon_N_minosTracks, &b_CCProtonPi0_muon_N_minosTracks);
   fChain->SetBranchAddress("CCProtonPi0_muon_charge", &CCProtonPi0_muon_charge, &b_CCProtonPi0_muon_charge);
   fChain->SetBranchAddress("CCProtonPi0_muon_hasMinosMatchStub", &CCProtonPi0_muon_hasMinosMatchStub, &b_CCProtonPi0_muon_hasMinosMatchStub);
   fChain->SetBranchAddress("CCProtonPi0_muon_hasMinosMatchTrack", &CCProtonPi0_muon_hasMinosMatchTrack, &b_CCProtonPi0_muon_hasMinosMatchTrack);
   fChain->SetBranchAddress("CCProtonPi0_muon_minervaTrack_types", &CCProtonPi0_muon_minervaTrack_types, &b_CCProtonPi0_muon_minervaTrack_types);
   fChain->SetBranchAddress("CCProtonPi0_muon_minosTrackQuality", &CCProtonPi0_muon_minosTrackQuality, &b_CCProtonPi0_muon_minosTrackQuality);
   fChain->SetBranchAddress("CCProtonPi0_muon_roadUpstreamPlanes", &CCProtonPi0_muon_roadUpstreamPlanes, &b_CCProtonPi0_muon_roadUpstreamPlanes);
   fChain->SetBranchAddress("CCProtonPi0_ntrajMuonProng", &CCProtonPi0_ntrajMuonProng, &b_CCProtonPi0_ntrajMuonProng);
   fChain->SetBranchAddress("CCProtonPi0_r_minos_trk_vtx_plane", &CCProtonPi0_r_minos_trk_vtx_plane, &b_CCProtonPi0_r_minos_trk_vtx_plane);
   fChain->SetBranchAddress("CCProtonPi0_t_minos_trk_numFSMuons", &CCProtonPi0_t_minos_trk_numFSMuons, &b_CCProtonPi0_t_minos_trk_numFSMuons);
   fChain->SetBranchAddress("CCProtonPi0_t_minos_trk_primFSLeptonPDG", &CCProtonPi0_t_minos_trk_primFSLeptonPDG, &b_CCProtonPi0_t_minos_trk_primFSLeptonPDG);
   fChain->SetBranchAddress("CCProtonPi0_trajMuonProngPDG", &CCProtonPi0_trajMuonProngPDG, &b_CCProtonPi0_trajMuonProngPDG);
   fChain->SetBranchAddress("CCProtonPi0_trajMuonProngPrimary", &CCProtonPi0_trajMuonProngPrimary, &b_CCProtonPi0_trajMuonProngPrimary);
   fChain->SetBranchAddress("CCProtonPi0_vtx_module", &CCProtonPi0_vtx_module, &b_CCProtonPi0_vtx_module);
   fChain->SetBranchAddress("CCProtonPi0_vtx_plane", &CCProtonPi0_vtx_plane, &b_CCProtonPi0_vtx_plane);
   fChain->SetBranchAddress("CCProtonPi0_QSq", &CCProtonPi0_QSq, &b_CCProtonPi0_QSq);
   fChain->SetBranchAddress("CCProtonPi0_WSq", &CCProtonPi0_WSq, &b_CCProtonPi0_WSq);
   fChain->SetBranchAddress("CCProtonPi0_endMuonTrajMomentum", &CCProtonPi0_endMuonTrajMomentum, &b_CCProtonPi0_endMuonTrajMomentum);
   fChain->SetBranchAddress("CCProtonPi0_endMuonTrajXPosition", &CCProtonPi0_endMuonTrajXPosition, &b_CCProtonPi0_endMuonTrajXPosition);
   fChain->SetBranchAddress("CCProtonPi0_endMuonTrajYPosition", &CCProtonPi0_endMuonTrajYPosition, &b_CCProtonPi0_endMuonTrajYPosition);
   fChain->SetBranchAddress("CCProtonPi0_endMuonTrajZPosition", &CCProtonPi0_endMuonTrajZPosition, &b_CCProtonPi0_endMuonTrajZPosition);
   fChain->SetBranchAddress("CCProtonPi0_minos_trk_bave", &CCProtonPi0_minos_trk_bave, &b_CCProtonPi0_minos_trk_bave);
   fChain->SetBranchAddress("CCProtonPi0_minos_trk_chi2", &CCProtonPi0_minos_trk_chi2, &b_CCProtonPi0_minos_trk_chi2);
   fChain->SetBranchAddress("CCProtonPi0_minos_trk_end_u", &CCProtonPi0_minos_trk_end_u, &b_CCProtonPi0_minos_trk_end_u);
   fChain->SetBranchAddress("CCProtonPi0_minos_trk_end_v", &CCProtonPi0_minos_trk_end_v, &b_CCProtonPi0_minos_trk_end_v);
   fChain->SetBranchAddress("CCProtonPi0_minos_trk_end_x", &CCProtonPi0_minos_trk_end_x, &b_CCProtonPi0_minos_trk_end_x);
   fChain->SetBranchAddress("CCProtonPi0_minos_trk_end_y", &CCProtonPi0_minos_trk_end_y, &b_CCProtonPi0_minos_trk_end_y);
   fChain->SetBranchAddress("CCProtonPi0_minos_trk_end_z", &CCProtonPi0_minos_trk_end_z, &b_CCProtonPi0_minos_trk_end_z);
   fChain->SetBranchAddress("CCProtonPi0_minos_trk_eqp", &CCProtonPi0_minos_trk_eqp, &b_CCProtonPi0_minos_trk_eqp);
   fChain->SetBranchAddress("CCProtonPi0_minos_trk_eqp_qp", &CCProtonPi0_minos_trk_eqp_qp, &b_CCProtonPi0_minos_trk_eqp_qp);
   fChain->SetBranchAddress("CCProtonPi0_minos_trk_fit_pass", &CCProtonPi0_minos_trk_fit_pass, &b_CCProtonPi0_minos_trk_fit_pass);
   fChain->SetBranchAddress("CCProtonPi0_minos_trk_ndf", &CCProtonPi0_minos_trk_ndf, &b_CCProtonPi0_minos_trk_ndf);
   fChain->SetBranchAddress("CCProtonPi0_minos_trk_p", &CCProtonPi0_minos_trk_p, &b_CCProtonPi0_minos_trk_p);
   fChain->SetBranchAddress("CCProtonPi0_minos_trk_p_curvature", &CCProtonPi0_minos_trk_p_curvature, &b_CCProtonPi0_minos_trk_p_curvature);
   fChain->SetBranchAddress("CCProtonPi0_minos_trk_p_range", &CCProtonPi0_minos_trk_p_range, &b_CCProtonPi0_minos_trk_p_range);
   fChain->SetBranchAddress("CCProtonPi0_minos_trk_qp", &CCProtonPi0_minos_trk_qp, &b_CCProtonPi0_minos_trk_qp);
   fChain->SetBranchAddress("CCProtonPi0_minos_trk_vtx_x", &CCProtonPi0_minos_trk_vtx_x, &b_CCProtonPi0_minos_trk_vtx_x);
   fChain->SetBranchAddress("CCProtonPi0_minos_trk_vtx_y", &CCProtonPi0_minos_trk_vtx_y, &b_CCProtonPi0_minos_trk_vtx_y);
   fChain->SetBranchAddress("CCProtonPi0_minos_trk_vtx_z", &CCProtonPi0_minos_trk_vtx_z, &b_CCProtonPi0_minos_trk_vtx_z);
   fChain->SetBranchAddress("CCProtonPi0_muon_E", &CCProtonPi0_muon_E, &b_CCProtonPi0_muon_E);
   fChain->SetBranchAddress("CCProtonPi0_muon_E_shift", &CCProtonPi0_muon_E_shift, &b_CCProtonPi0_muon_E_shift);
   fChain->SetBranchAddress("CCProtonPi0_muon_muScore", &CCProtonPi0_muon_muScore, &b_CCProtonPi0_muon_muScore);
   fChain->SetBranchAddress("CCProtonPi0_muon_p", &CCProtonPi0_muon_p, &b_CCProtonPi0_muon_p);
   fChain->SetBranchAddress("CCProtonPi0_muon_px", &CCProtonPi0_muon_px, &b_CCProtonPi0_muon_px);
   fChain->SetBranchAddress("CCProtonPi0_muon_py", &CCProtonPi0_muon_py, &b_CCProtonPi0_muon_py);
   fChain->SetBranchAddress("CCProtonPi0_muon_pz", &CCProtonPi0_muon_pz, &b_CCProtonPi0_muon_pz);
   fChain->SetBranchAddress("CCProtonPi0_muon_qp", &CCProtonPi0_muon_qp, &b_CCProtonPi0_muon_qp);
   fChain->SetBranchAddress("CCProtonPi0_muon_qpqpe", &CCProtonPi0_muon_qpqpe, &b_CCProtonPi0_muon_qpqpe);
   fChain->SetBranchAddress("CCProtonPi0_muon_roadUpstreamEnergy", &CCProtonPi0_muon_roadUpstreamEnergy, &b_CCProtonPi0_muon_roadUpstreamEnergy);
   fChain->SetBranchAddress("CCProtonPi0_muon_theta", &CCProtonPi0_muon_theta, &b_CCProtonPi0_muon_theta);
   fChain->SetBranchAddress("CCProtonPi0_muon_theta_biasDown", &CCProtonPi0_muon_theta_biasDown, &b_CCProtonPi0_muon_theta_biasDown);
   fChain->SetBranchAddress("CCProtonPi0_muon_theta_biasUp", &CCProtonPi0_muon_theta_biasUp, &b_CCProtonPi0_muon_theta_biasUp);
   fChain->SetBranchAddress("CCProtonPi0_neutrino_E", &CCProtonPi0_neutrino_E, &b_CCProtonPi0_neutrino_E);
   fChain->SetBranchAddress("CCProtonPi0_neutrino_E_Cal", &CCProtonPi0_neutrino_E_Cal, &b_CCProtonPi0_neutrino_E_Cal);
   fChain->SetBranchAddress("CCProtonPi0_r_minos_trk_bdL", &CCProtonPi0_r_minos_trk_bdL, &b_CCProtonPi0_r_minos_trk_bdL);
   fChain->SetBranchAddress("CCProtonPi0_r_minos_trk_end_dcosx", &CCProtonPi0_r_minos_trk_end_dcosx, &b_CCProtonPi0_r_minos_trk_end_dcosx);
   fChain->SetBranchAddress("CCProtonPi0_r_minos_trk_end_dcosy", &CCProtonPi0_r_minos_trk_end_dcosy, &b_CCProtonPi0_r_minos_trk_end_dcosy);
   fChain->SetBranchAddress("CCProtonPi0_r_minos_trk_end_dcosz", &CCProtonPi0_r_minos_trk_end_dcosz, &b_CCProtonPi0_r_minos_trk_end_dcosz);
   fChain->SetBranchAddress("CCProtonPi0_r_minos_trk_vtx_dcosx", &CCProtonPi0_r_minos_trk_vtx_dcosx, &b_CCProtonPi0_r_minos_trk_vtx_dcosx);
   fChain->SetBranchAddress("CCProtonPi0_r_minos_trk_vtx_dcosy", &CCProtonPi0_r_minos_trk_vtx_dcosy, &b_CCProtonPi0_r_minos_trk_vtx_dcosy);
   fChain->SetBranchAddress("CCProtonPi0_r_minos_trk_vtx_dcosz", &CCProtonPi0_r_minos_trk_vtx_dcosz, &b_CCProtonPi0_r_minos_trk_vtx_dcosz);
   fChain->SetBranchAddress("CCProtonPi0_t_minos_trk_primFSLepMinosInitProjPx", &CCProtonPi0_t_minos_trk_primFSLepMinosInitProjPx, &b_CCProtonPi0_t_minos_trk_primFSLepMinosInitProjPx);
   fChain->SetBranchAddress("CCProtonPi0_t_minos_trk_primFSLepMinosInitProjPy", &CCProtonPi0_t_minos_trk_primFSLepMinosInitProjPy, &b_CCProtonPi0_t_minos_trk_primFSLepMinosInitProjPy);
   fChain->SetBranchAddress("CCProtonPi0_t_minos_trk_primFSLepMinosInitProjPz", &CCProtonPi0_t_minos_trk_primFSLepMinosInitProjPz, &b_CCProtonPi0_t_minos_trk_primFSLepMinosInitProjPz);
   fChain->SetBranchAddress("CCProtonPi0_t_minos_trk_primFSLepMinosInitProjX", &CCProtonPi0_t_minos_trk_primFSLepMinosInitProjX, &b_CCProtonPi0_t_minos_trk_primFSLepMinosInitProjX);
   fChain->SetBranchAddress("CCProtonPi0_t_minos_trk_primFSLepMinosInitProjY", &CCProtonPi0_t_minos_trk_primFSLepMinosInitProjY, &b_CCProtonPi0_t_minos_trk_primFSLepMinosInitProjY);
   fChain->SetBranchAddress("CCProtonPi0_t_minos_trk_primFSLepMinosInitProjZ", &CCProtonPi0_t_minos_trk_primFSLepMinosInitProjZ, &b_CCProtonPi0_t_minos_trk_primFSLepMinosInitProjZ);
   fChain->SetBranchAddress("CCProtonPi0_t_minos_trk_primFSLepMnvFinalPx", &CCProtonPi0_t_minos_trk_primFSLepMnvFinalPx, &b_CCProtonPi0_t_minos_trk_primFSLepMnvFinalPx);
   fChain->SetBranchAddress("CCProtonPi0_t_minos_trk_primFSLepMnvFinalPy", &CCProtonPi0_t_minos_trk_primFSLepMnvFinalPy, &b_CCProtonPi0_t_minos_trk_primFSLepMnvFinalPy);
   fChain->SetBranchAddress("CCProtonPi0_t_minos_trk_primFSLepMnvFinalPz", &CCProtonPi0_t_minos_trk_primFSLepMnvFinalPz, &b_CCProtonPi0_t_minos_trk_primFSLepMnvFinalPz);
   fChain->SetBranchAddress("CCProtonPi0_t_minos_trk_primFSLepMnvFinalX", &CCProtonPi0_t_minos_trk_primFSLepMnvFinalX, &b_CCProtonPi0_t_minos_trk_primFSLepMnvFinalX);
   fChain->SetBranchAddress("CCProtonPi0_t_minos_trk_primFSLepMnvFinalY", &CCProtonPi0_t_minos_trk_primFSLepMnvFinalY, &b_CCProtonPi0_t_minos_trk_primFSLepMnvFinalY);
   fChain->SetBranchAddress("CCProtonPi0_t_minos_trk_primFSLepMnvFinalZ", &CCProtonPi0_t_minos_trk_primFSLepMnvFinalZ, &b_CCProtonPi0_t_minos_trk_primFSLepMnvFinalZ);
   fChain->SetBranchAddress("CCProtonPi0_t_minos_trk_primFSLepMnvInitPx", &CCProtonPi0_t_minos_trk_primFSLepMnvInitPx, &b_CCProtonPi0_t_minos_trk_primFSLepMnvInitPx);
   fChain->SetBranchAddress("CCProtonPi0_t_minos_trk_primFSLepMnvInitPy", &CCProtonPi0_t_minos_trk_primFSLepMnvInitPy, &b_CCProtonPi0_t_minos_trk_primFSLepMnvInitPy);
   fChain->SetBranchAddress("CCProtonPi0_t_minos_trk_primFSLepMnvInitPz", &CCProtonPi0_t_minos_trk_primFSLepMnvInitPz, &b_CCProtonPi0_t_minos_trk_primFSLepMnvInitPz);
   fChain->SetBranchAddress("CCProtonPi0_t_minos_trk_primFSLepMnvInitX", &CCProtonPi0_t_minos_trk_primFSLepMnvInitX, &b_CCProtonPi0_t_minos_trk_primFSLepMnvInitX);
   fChain->SetBranchAddress("CCProtonPi0_t_minos_trk_primFSLepMnvInitY", &CCProtonPi0_t_minos_trk_primFSLepMnvInitY, &b_CCProtonPi0_t_minos_trk_primFSLepMnvInitY);
   fChain->SetBranchAddress("CCProtonPi0_t_minos_trk_primFSLepMnvInitZ", &CCProtonPi0_t_minos_trk_primFSLepMnvInitZ, &b_CCProtonPi0_t_minos_trk_primFSLepMnvInitZ);
   fChain->SetBranchAddress("CCProtonPi0_total_E", &CCProtonPi0_total_E, &b_CCProtonPi0_total_E);
   fChain->SetBranchAddress("CCProtonPi0_total_px", &CCProtonPi0_total_px, &b_CCProtonPi0_total_px);
   fChain->SetBranchAddress("CCProtonPi0_total_py", &CCProtonPi0_total_py, &b_CCProtonPi0_total_py);
   fChain->SetBranchAddress("CCProtonPi0_total_pz", &CCProtonPi0_total_pz, &b_CCProtonPi0_total_pz);
   fChain->SetBranchAddress("CCProtonPi0_trajMuonPhi", &CCProtonPi0_trajMuonPhi, &b_CCProtonPi0_trajMuonPhi);
   fChain->SetBranchAddress("CCProtonPi0_trajMuonProngEnergy", &CCProtonPi0_trajMuonProngEnergy, &b_CCProtonPi0_trajMuonProngEnergy);
   fChain->SetBranchAddress("CCProtonPi0_trajMuonProngMomentum", &CCProtonPi0_trajMuonProngMomentum, &b_CCProtonPi0_trajMuonProngMomentum);
   fChain->SetBranchAddress("CCProtonPi0_trajMuonProngPSelf", &CCProtonPi0_trajMuonProngPSelf, &b_CCProtonPi0_trajMuonProngPSelf);
   fChain->SetBranchAddress("CCProtonPi0_trajMuonProngPx", &CCProtonPi0_trajMuonProngPx, &b_CCProtonPi0_trajMuonProngPx);
   fChain->SetBranchAddress("CCProtonPi0_trajMuonProngPy", &CCProtonPi0_trajMuonProngPy, &b_CCProtonPi0_trajMuonProngPy);
   fChain->SetBranchAddress("CCProtonPi0_trajMuonProngPz", &CCProtonPi0_trajMuonProngPz, &b_CCProtonPi0_trajMuonProngPz);
   fChain->SetBranchAddress("CCProtonPi0_trajMuonTheta", &CCProtonPi0_trajMuonTheta, &b_CCProtonPi0_trajMuonTheta);
   fChain->SetBranchAddress("CCProtonPi0_vtx_x", &CCProtonPi0_vtx_x, &b_CCProtonPi0_vtx_x);
   fChain->SetBranchAddress("CCProtonPi0_vtx_y", &CCProtonPi0_vtx_y, &b_CCProtonPi0_vtx_y);
   fChain->SetBranchAddress("CCProtonPi0_vtx_z", &CCProtonPi0_vtx_z, &b_CCProtonPi0_vtx_z);
   fChain->SetBranchAddress("CCProtonPi0_isProtonInsideOD", CCProtonPi0_isProtonInsideOD, &b_CCProtonPi0_isProtonInsideOD);
   fChain->SetBranchAddress("CCProtonPi0_ntrajProtonProng", CCProtonPi0_ntrajProtonProng, &b_CCProtonPi0_ntrajProtonProng);
   fChain->SetBranchAddress("CCProtonPi0_proton_isRecoGood", CCProtonPi0_proton_isRecoGood, &b_CCProtonPi0_proton_isRecoGood);
   fChain->SetBranchAddress("CCProtonPi0_proton_kinked", CCProtonPi0_proton_kinked, &b_CCProtonPi0_proton_kinked);
   fChain->SetBranchAddress("CCProtonPi0_proton_odMatch", CCProtonPi0_proton_odMatch, &b_CCProtonPi0_proton_odMatch);
   fChain->SetBranchAddress("CCProtonPi0_trajProtonProngPDG", CCProtonPi0_trajProtonProngPDG, &b_CCProtonPi0_trajProtonProngPDG);
   fChain->SetBranchAddress("CCProtonPi0_trajProtonProngPrimary", CCProtonPi0_trajProtonProngPrimary, &b_CCProtonPi0_trajProtonProngPrimary);
   fChain->SetBranchAddress("CCProtonPi0_endProtonTrajMomentum", CCProtonPi0_endProtonTrajMomentum, &b_CCProtonPi0_endProtonTrajMomentum);
   fChain->SetBranchAddress("CCProtonPi0_endProtonTrajXPosition", CCProtonPi0_endProtonTrajXPosition, &b_CCProtonPi0_endProtonTrajXPosition);
   fChain->SetBranchAddress("CCProtonPi0_endProtonTrajYPosition", CCProtonPi0_endProtonTrajYPosition, &b_CCProtonPi0_endProtonTrajYPosition);
   fChain->SetBranchAddress("CCProtonPi0_endProtonTrajZPosition", CCProtonPi0_endProtonTrajZPosition, &b_CCProtonPi0_endProtonTrajZPosition);
   fChain->SetBranchAddress("CCProtonPi0_pionScore", CCProtonPi0_pionScore, &b_CCProtonPi0_pionScore);
   fChain->SetBranchAddress("CCProtonPi0_protonScore", CCProtonPi0_protonScore, &b_CCProtonPi0_protonScore);
   fChain->SetBranchAddress("CCProtonPi0_protonScore_LLR", CCProtonPi0_protonScore_LLR, &b_CCProtonPi0_protonScore_LLR);
   fChain->SetBranchAddress("CCProtonPi0_proton_E", CCProtonPi0_proton_E, &b_CCProtonPi0_proton_E);
   fChain->SetBranchAddress("CCProtonPi0_proton_chi2_ndf", CCProtonPi0_proton_chi2_ndf, &b_CCProtonPi0_proton_chi2_ndf);
   fChain->SetBranchAddress("CCProtonPi0_proton_ekin", CCProtonPi0_proton_ekin, &b_CCProtonPi0_proton_ekin);
   fChain->SetBranchAddress("CCProtonPi0_proton_endPointX", CCProtonPi0_proton_endPointX, &b_CCProtonPi0_proton_endPointX);
   fChain->SetBranchAddress("CCProtonPi0_proton_endPointY", CCProtonPi0_proton_endPointY, &b_CCProtonPi0_proton_endPointY);
   fChain->SetBranchAddress("CCProtonPi0_proton_endPointZ", CCProtonPi0_proton_endPointZ, &b_CCProtonPi0_proton_endPointZ);
   fChain->SetBranchAddress("CCProtonPi0_proton_p", CCProtonPi0_proton_p, &b_CCProtonPi0_proton_p);
   fChain->SetBranchAddress("CCProtonPi0_proton_p_calCorrection", CCProtonPi0_proton_p_calCorrection, &b_CCProtonPi0_proton_p_calCorrection);
   fChain->SetBranchAddress("CCProtonPi0_proton_p_dEdXTool", CCProtonPi0_proton_p_dEdXTool, &b_CCProtonPi0_proton_p_dEdXTool);
   fChain->SetBranchAddress("CCProtonPi0_proton_p_visEnergy", CCProtonPi0_proton_p_visEnergy, &b_CCProtonPi0_proton_p_visEnergy);
   fChain->SetBranchAddress("CCProtonPi0_proton_phi", CCProtonPi0_proton_phi, &b_CCProtonPi0_proton_phi);
   fChain->SetBranchAddress("CCProtonPi0_proton_px", CCProtonPi0_proton_px, &b_CCProtonPi0_proton_px);
   fChain->SetBranchAddress("CCProtonPi0_proton_py", CCProtonPi0_proton_py, &b_CCProtonPi0_proton_py);
   fChain->SetBranchAddress("CCProtonPi0_proton_pz", CCProtonPi0_proton_pz, &b_CCProtonPi0_proton_pz);
   fChain->SetBranchAddress("CCProtonPi0_proton_startPointX", CCProtonPi0_proton_startPointX, &b_CCProtonPi0_proton_startPointX);
   fChain->SetBranchAddress("CCProtonPi0_proton_startPointY", CCProtonPi0_proton_startPointY, &b_CCProtonPi0_proton_startPointY);
   fChain->SetBranchAddress("CCProtonPi0_proton_startPointZ", CCProtonPi0_proton_startPointZ, &b_CCProtonPi0_proton_startPointZ);
   fChain->SetBranchAddress("CCProtonPi0_proton_theta", CCProtonPi0_proton_theta, &b_CCProtonPi0_proton_theta);
   fChain->SetBranchAddress("CCProtonPi0_proton_thetaX", CCProtonPi0_proton_thetaX, &b_CCProtonPi0_proton_thetaX);
   fChain->SetBranchAddress("CCProtonPi0_proton_thetaY", CCProtonPi0_proton_thetaY, &b_CCProtonPi0_proton_thetaY);
   fChain->SetBranchAddress("CCProtonPi0_trajProtonPhi", CCProtonPi0_trajProtonPhi, &b_CCProtonPi0_trajProtonPhi);
   fChain->SetBranchAddress("CCProtonPi0_trajProtonProngEnergy", CCProtonPi0_trajProtonProngEnergy, &b_CCProtonPi0_trajProtonProngEnergy);
   fChain->SetBranchAddress("CCProtonPi0_trajProtonProngMomentum", CCProtonPi0_trajProtonProngMomentum, &b_CCProtonPi0_trajProtonProngMomentum);
   fChain->SetBranchAddress("CCProtonPi0_trajProtonProngPSelf", CCProtonPi0_trajProtonProngPSelf, &b_CCProtonPi0_trajProtonProngPSelf);
   fChain->SetBranchAddress("CCProtonPi0_trajProtonProngPx", CCProtonPi0_trajProtonProngPx, &b_CCProtonPi0_trajProtonProngPx);
   fChain->SetBranchAddress("CCProtonPi0_trajProtonProngPy", CCProtonPi0_trajProtonProngPy, &b_CCProtonPi0_trajProtonProngPy);
   fChain->SetBranchAddress("CCProtonPi0_trajProtonProngPz", CCProtonPi0_trajProtonProngPz, &b_CCProtonPi0_trajProtonProngPz);
   fChain->SetBranchAddress("CCProtonPi0_trajProtonTheta", CCProtonPi0_trajProtonTheta, &b_CCProtonPi0_trajProtonTheta);
   fChain->SetBranchAddress("ev_run", &ev_run, &b_ev_run);
   fChain->SetBranchAddress("ev_subrun", &ev_subrun, &b_ev_subrun);
   fChain->SetBranchAddress("ev_detector", &ev_detector, &b_ev_detector);
   fChain->SetBranchAddress("ev_triggerType", &ev_triggerType, &b_ev_triggerType);
   fChain->SetBranchAddress("ev_gate", &ev_gate, &b_ev_gate);
   fChain->SetBranchAddress("ev_global_gate", &ev_global_gate, &b_ev_global_gate);
   fChain->SetBranchAddress("ev_gps_time_sec", &ev_gps_time_sec, &b_ev_gps_time_sec);
   fChain->SetBranchAddress("ev_gps_time_usec", &ev_gps_time_usec, &b_ev_gps_time_usec);
   fChain->SetBranchAddress("mc_run", &mc_run, &b_mc_run);
   fChain->SetBranchAddress("mc_subrun", &mc_subrun, &b_mc_subrun);
   fChain->SetBranchAddress("mc_nInteractions", &mc_nInteractions, &b_mc_nInteractions);
   fChain->SetBranchAddress("mc_MIState", &mc_MIState, &b_mc_MIState);
   fChain->SetBranchAddress("mc_pot", &mc_pot, &b_mc_pot);
   fChain->SetBranchAddress("mc_beamConfig", &mc_beamConfig, &b_mc_beamConfig);
   fChain->SetBranchAddress("mc_processType", &mc_processType, &b_mc_processType);
   fChain->SetBranchAddress("mc_nthEvtInSpill", &mc_nthEvtInSpill, &b_mc_nthEvtInSpill);
   fChain->SetBranchAddress("mc_nthEvtInFile", &mc_nthEvtInFile, &b_mc_nthEvtInFile);
   fChain->SetBranchAddress("mc_intType", &mc_intType, &b_mc_intType);
   fChain->SetBranchAddress("mc_current", &mc_current, &b_mc_current);
   fChain->SetBranchAddress("mc_charm", &mc_charm, &b_mc_charm);
   fChain->SetBranchAddress("mc_weight", &mc_weight, &b_mc_weight);
   fChain->SetBranchAddress("mc_XSec", &mc_XSec, &b_mc_XSec);
   fChain->SetBranchAddress("mc_diffXSec", &mc_diffXSec, &b_mc_diffXSec);
   fChain->SetBranchAddress("mc_incoming", &mc_incoming, &b_mc_incoming);
   fChain->SetBranchAddress("mc_fluxDriverProb", &mc_fluxDriverProb, &b_mc_fluxDriverProb);
   fChain->SetBranchAddress("mc_targetNucleus", &mc_targetNucleus, &b_mc_targetNucleus);
   fChain->SetBranchAddress("mc_targetZ", &mc_targetZ, &b_mc_targetZ);
   fChain->SetBranchAddress("mc_targetA", &mc_targetA, &b_mc_targetA);
   fChain->SetBranchAddress("mc_targetNucleon", &mc_targetNucleon, &b_mc_targetNucleon);
   fChain->SetBranchAddress("mc_struckQuark", &mc_struckQuark, &b_mc_struckQuark);
   fChain->SetBranchAddress("mc_seaQuark", &mc_seaQuark, &b_mc_seaQuark);
   fChain->SetBranchAddress("mc_resID", &mc_resID, &b_mc_resID);
   fChain->SetBranchAddress("mc_primaryLepton", &mc_primaryLepton, &b_mc_primaryLepton);
   fChain->SetBranchAddress("mc_incomingE", &mc_incomingE, &b_mc_incomingE);
   fChain->SetBranchAddress("mc_Bjorkenx", &mc_Bjorkenx, &b_mc_Bjorkenx);
   fChain->SetBranchAddress("mc_Bjorkeny", &mc_Bjorkeny, &b_mc_Bjorkeny);
   fChain->SetBranchAddress("mc_Q2", &mc_Q2, &b_mc_Q2);
   fChain->SetBranchAddress("mc_nuT", &mc_nuT, &b_mc_nuT);
   fChain->SetBranchAddress("mc_w", &mc_w, &b_mc_w);
   fChain->SetBranchAddress("mc_vtx", mc_vtx, &b_mc_vtx);
   fChain->SetBranchAddress("mc_incomingPartVec", mc_incomingPartVec, &b_mc_incomingPartVec);
   fChain->SetBranchAddress("mc_initNucVec", mc_initNucVec, &b_mc_initNucVec);
   fChain->SetBranchAddress("mc_primFSLepton", mc_primFSLepton, &b_mc_primFSLepton);
   fChain->SetBranchAddress("mc_nFSPart", &mc_nFSPart, &b_mc_nFSPart);
   fChain->SetBranchAddress("mc_FSPartPx", mc_FSPartPx, &b_mc_FSPartPx);
   fChain->SetBranchAddress("mc_FSPartPy", mc_FSPartPy, &b_mc_FSPartPy);
   fChain->SetBranchAddress("mc_FSPartPz", mc_FSPartPz, &b_mc_FSPartPz);
   fChain->SetBranchAddress("mc_FSPartE", mc_FSPartE, &b_mc_FSPartE);
   fChain->SetBranchAddress("mc_FSPartPDG", mc_FSPartPDG, &b_mc_FSPartPDG);
   fChain->SetBranchAddress("mc_er_nPart", &mc_er_nPart, &b_mc_er_nPart);
   fChain->SetBranchAddress("mc_er_ID", mc_er_ID, &b_mc_er_ID);
   fChain->SetBranchAddress("mc_er_status", mc_er_status, &b_mc_er_status);
   fChain->SetBranchAddress("mc_er_posInNucX", mc_er_posInNucX, &b_mc_er_posInNucX);
   fChain->SetBranchAddress("mc_er_posInNucY", mc_er_posInNucY, &b_mc_er_posInNucY);
   fChain->SetBranchAddress("mc_er_posInNucZ", mc_er_posInNucZ, &b_mc_er_posInNucZ);
   fChain->SetBranchAddress("mc_er_Px", mc_er_Px, &b_mc_er_Px);
   fChain->SetBranchAddress("mc_er_Py", mc_er_Py, &b_mc_er_Py);
   fChain->SetBranchAddress("mc_er_Pz", mc_er_Pz, &b_mc_er_Pz);
   fChain->SetBranchAddress("mc_er_E", mc_er_E, &b_mc_er_E);
   fChain->SetBranchAddress("mc_er_FD", mc_er_FD, &b_mc_er_FD);
   fChain->SetBranchAddress("mc_er_LD", mc_er_LD, &b_mc_er_LD);
   fChain->SetBranchAddress("mc_er_mother", mc_er_mother, &b_mc_er_mother);
   fChain->SetBranchAddress("mc_fr_nNuAncestorIDs", &mc_fr_nNuAncestorIDs, &b_mc_fr_nNuAncestorIDs);
   fChain->SetBranchAddress("mc_fr_nuAncestorIDs", mc_fr_nuAncestorIDs, &b_mc_fr_nuAncestorIDs);
   fChain->SetBranchAddress("mc_fr_nuParentID", &mc_fr_nuParentID, &b_mc_fr_nuParentID);
   fChain->SetBranchAddress("mc_fr_decMode", &mc_fr_decMode, &b_mc_fr_decMode);
   fChain->SetBranchAddress("mc_fr_primProtonVtx", mc_fr_primProtonVtx, &b_mc_fr_primProtonVtx);
   fChain->SetBranchAddress("mc_fr_primProtonP", mc_fr_primProtonP, &b_mc_fr_primProtonP);
   fChain->SetBranchAddress("mc_fr_nuParentDecVtx", mc_fr_nuParentDecVtx, &b_mc_fr_nuParentDecVtx);
   fChain->SetBranchAddress("mc_fr_nuParentProdVtx", mc_fr_nuParentProdVtx, &b_mc_fr_nuParentProdVtx);
   fChain->SetBranchAddress("mc_fr_nuParentProdP", mc_fr_nuParentProdP, &b_mc_fr_nuParentProdP);
   fChain->SetBranchAddress("mc_cvweight_total", &mc_cvweight_total, &b_mc_cvweight_total);
   fChain->SetBranchAddress("wgt", &wgt, &b_wgt);
   fChain->SetBranchAddress("mc_cvweight_totalFlux", &mc_cvweight_totalFlux, &b_mc_cvweight_totalFlux);
   fChain->SetBranchAddress("mc_cvweight_totalXsec", &mc_cvweight_totalXsec, &b_mc_cvweight_totalXsec);
   fChain->SetBranchAddress("mc_cvweight_NA49", &mc_cvweight_NA49, &b_mc_cvweight_NA49);
   fChain->SetBranchAddress("mc_wgt_GENIE_sz", &mc_wgt_GENIE_sz, &b_mc_wgt_GENIE_sz);
   fChain->SetBranchAddress("mc_wgt_GENIE", mc_wgt_GENIE, &b_mc_wgt_GENIE);
   fChain->SetBranchAddress("mc_wgt_Flux_Tertiary_sz", &mc_wgt_Flux_Tertiary_sz, &b_mc_wgt_Flux_Tertiary_sz);
   fChain->SetBranchAddress("mc_wgt_Flux_Tertiary", mc_wgt_Flux_Tertiary, &b_mc_wgt_Flux_Tertiary);
   fChain->SetBranchAddress("mc_wgt_Flux_BeamFocus_sz", &mc_wgt_Flux_BeamFocus_sz, &b_mc_wgt_Flux_BeamFocus_sz);
   fChain->SetBranchAddress("mc_wgt_Flux_BeamFocus", mc_wgt_Flux_BeamFocus, &b_mc_wgt_Flux_BeamFocus);
   fChain->SetBranchAddress("mc_wgt_Flux_NA49_sz", &mc_wgt_Flux_NA49_sz, &b_mc_wgt_Flux_NA49_sz);
   fChain->SetBranchAddress("mc_wgt_Flux_NA49", mc_wgt_Flux_NA49, &b_mc_wgt_Flux_NA49);
   fChain->SetBranchAddress("n_prongs", &n_prongs, &b_n_prongs);
   fChain->SetBranchAddress("prong_nParticles", prong_nParticles, &b_prong_nParticles);
   fChain->SetBranchAddress("prong_part_score", prong_part_score, &b_prong_part_score);
   fChain->SetBranchAddress("prong_part_mass", prong_part_mass, &b_prong_part_mass);
   fChain->SetBranchAddress("prong_part_charge", prong_part_charge, &b_prong_part_charge);
   fChain->SetBranchAddress("prong_part_pid", prong_part_pid, &b_prong_part_pid);
   fChain->SetBranchAddress("prong_part_E", &prong_part_E, &b_prong_part_E);
   fChain->SetBranchAddress("prong_part_pos", &prong_part_pos, &b_prong_part_pos);
}

Analyzer::~Analyzer()
{
    closeTextFiles();
}

Int_t Analyzer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t Analyzer::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
   }
   return centry;
}



#endif //Analyzer_cpp
