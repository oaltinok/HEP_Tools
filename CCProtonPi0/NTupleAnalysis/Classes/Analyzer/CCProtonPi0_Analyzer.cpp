/*
    See CCProtonPi0_Analyzer.h header or Class Information
*/
#ifndef CCProtonPi0_Analyzer_cpp
#define CCProtonPi0_Analyzer_cpp

#include "CCProtonPi0_Analyzer.h"

using namespace std;

void CCProtonPi0_Analyzer::specifyRunTime()
{
    applyMaxEvents = false;

    nMaxEvents = 100000;
    
    // Control Flow
    isDataAnalysis  = true;
    isScanRun = false;
    is_pID_Studies  = false;
    isMichelStudy = false;
    applyOtherMichelCuts = false;
    writeFSParticleMomentum = false;
    
    applyProtonScore = true;
    pID_KE_Limit = 300.0;
    minProtonScore_LLR = 10.0;
    minPIDDiff = 0.45;

    applyPi0BlobCuts = true;

    applyPhotonDistance = true;
    minPhotonDistance = 14; //cm
    
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
    
    hasParticleTruthInfo = true;
}

void CCProtonPi0_Analyzer::reduce(string playlist)
{

    string rootDir;
    if (m_isMC) rootDir = Folder_List::rootOut + Folder_List::MC + Folder_List::reduced + "ReducedNTuple2.root";
    else rootDir = Folder_List::rootOut + Folder_List::Data + Folder_List::reduced + "ReducedNTuple.root";

    cout<<"Reducing NTuple Files to a single file"<<endl;
    cout<<"\tRoot File: "<<rootDir<<endl;
    TFile* f = new TFile(rootDir.c_str(),"RECREATE");
    
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
        if (jentry%25000 == 0) cout<<"\tEntry "<<jentry<<endl;
        
        // Get Cut Statistics
        isPassedAllCuts = getCutStatistics();
        if( !isPassedAllCuts ) continue;

        tree->Fill();
    }
    
    cutList.writeCutTable();
    cutList.writeRootFile();
   
    cout<<">> Writing "<<rootDir<<endl;
    tree->AutoSave();    
    f->Write();
}


void CCProtonPi0_Analyzer::analyze(string playlist)
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
        
        // Update scanFileName if running for scan
        //if(isScanRun) UpdateScanFileName();
       
        // Decide to Analyze Event or NOT
        if( !analyzeEvent() ) continue;
       
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
            pIDTool.FillHistograms(CCProtonPi0_proton_LLRScore,CCProtonPi0_proton_protonScore,CCProtonPi0_proton_pionScore,
                                    CCProtonPi0_trajProtonProngPDG[indRecoProton],CCProtonPi0_proton_E );
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
    bckgTool.writeBackgroundTable();
    getPi0Family(); // Pi0 Family Information written inside FailFile
    
    //--------------------------------------------------------------------------
    // Write Root Files
    //--------------------------------------------------------------------------
    interaction.writeHistograms();
    muon.writeHistograms();
    proton.writeHistograms();
    pi0.writeHistograms();
    pi0Blob.write_RootFile();
    if(is_pID_Studies) pIDTool.write_RootFile();
    if(isMichelStudy) michelTool.write_RootFile();
    
}

    
//------------------------------------------------------------------------------
//  Constructor
//------------------------------------------------------------------------------
CCProtonPi0_Analyzer::CCProtonPi0_Analyzer(int nMode, bool isMC) : 
    CCProtonPi0_NTupleAnalysis(nMode),
    interaction(nMode, isMC),
    muon(nMode, isMC),
    proton(nMode, isMC),
    pi0(nMode, isMC),
    pi0Blob(nMode, isMC),
    pIDTool(nMode, isMC),
    bckgTool(nMode),
    michelTool(nMode, isMC),
    cutList(nMode)
{   
    cout<<"Initializing CCProtonPi0_Analyzer"<<endl;
    
    specifyRunTime();
    
    m_isMC = isMC;
    
    if (nMode == 0){
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


void CCProtonPi0_Analyzer::fillData()
{          
    // Fill Reconstructed Information
    fillInteractionReco();
    fillMuonReco();
    fillPi0Reco();
    if(nProngs > 1) fillProtonReco();
    
    // Fill Truth Information if Exist and Set Errors
    if( hasParticleTruthInfo ){
        //fillInteractionTrue();
        //fillMuonTrue();
        if(nProngs > 1) fillProtonTrue();
        //fillPi0True();
    }
}

//------------------------------------------------------------------------------
// analyzeEvent()
//      Decides to Analyze or Skip the Event using Analysis Mode
//------------------------------------------------------------------------------
bool CCProtonPi0_Analyzer::analyzeEvent()
{
    bool isAnalyzable;
    
    if ( anaMode == 1 ){
        if( truth_isSignal) isAnalyzable = true;
        else isAnalyzable = false;
    }else if ( anaMode == 2 ){
        if( truth_isSignal) isAnalyzable = false;
        else isAnalyzable = true;
    }else{
        isAnalyzable = true;
    }
    
    return isAnalyzable;
    
}

void CCProtonPi0_Analyzer::writeFSParticle4P(Long64_t nEntry)
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
                <<CCProtonPi0_proton_px<<", "
                <<CCProtonPi0_proton_py<<", "
                <<CCProtonPi0_proton_pz<<", "
                <<CCProtonPi0_proton_E<<" )"
                <<" Score = "<<CCProtonPi0_proton_LLRScore
                <<endl;
        failText[t]<<"Pi0 4-P = ( "
                <<CCProtonPi0_pi0_px<<", "
                <<CCProtonPi0_pi0_py<<", "
                <<CCProtonPi0_pi0_pz<<", "
                <<CCProtonPi0_pi0_E<<" )"
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

    // Increment Cut Number for Both Topologies at the Same Time Until nProngs Cut   
    cutList.nCut_All.increment(truth_isSignal, study1, study2); // Count All Events before Cuts

    if( Cut_Vertex_None == 1) return false;
    cutList.nCut_Vertex_None.increment(truth_isSignal, study1, study2);
        
    if( Cut_Vertex_Not_Reconstructable == 1) return false;
    cutList.nCut_Vertex_Not_Reconstructable.increment(truth_isSignal, study1, study2);
    
    if( Cut_Vertex_Not_Fiducial == 1) return false;
    cutList.nCut_Vertex_Not_Fiducial.increment(truth_isSignal, study1, study2);
   
    // Fill Vertex Count Cut Hist 
    cutList.hCut_vertexCount->Fill(vtx_total_count);
    if( vtx_total_count > 1) return false;
    cutList.nCut_Vertex_Count.increment(truth_isSignal, study1, study2);
    
    if( Cut_Muon_None == 1) return false;
    cutList.nCut_Muon_None.increment(truth_isSignal, study1, study2);

    // Fill Truth_W for MINOS Matched Signal Events
    //if (truth_isSignal) fill_mc_w(); 
        
    if( Cut_Muon_Not_Plausible == 1) return false;
    cutList.nCut_Muon_Not_Plausible.increment(truth_isSignal, study1, study2);
    
    if( Cut_Muon_Charge == 1) return false;
    cutList.nCut_Muon_Charge.increment(truth_isSignal, study1, study2);
   
    // Fill Michel Cut Hist
    if( Cut_Vertex_Michel_Exist == 1 || Cut_EndPoint_Michel_Exist == 1 || Cut_secEndPoint_Michel_Exist == 1 ){
        cutList.hCut_Michel->Fill(1);
    }else{
        cutList.hCut_Michel->Fill(0);
    } 
    if( Cut_Vertex_Michel_Exist == 1) return false;
    cutList.nCut_Vertex_Michel_Exist.increment(truth_isSignal, study1, study2);
    
    if( Cut_EndPoint_Michel_Exist == 1) return false;
    cutList.nCut_EndPoint_Michel_Exist.increment(truth_isSignal, study1, study2);
    
    if( Cut_secEndPoint_Michel_Exist == 1) return false;
    cutList.nCut_secEndPoint_Michel_Exist.increment(truth_isSignal, study1, study2);
    
    // Fill PreFilter Cut Hists
    cutList.hCut_eVis_nuclearTarget->Fill(evis_NuclearTarget);
    cutList.hCut_eVis_other->Fill(evis_TotalExceptNuclearTarget);
    if( Cut_PreFilter_Pi0 == 1) return false;
    cutList.nCut_PreFilter_Pi0.increment(truth_isSignal, study1, study2);
    
    if( Cut_ConeBlobs == 1) return false;
    cutList.nCut_ConeBlobs.increment(truth_isSignal, study1, study2);

    if ( Cut_BlobsBad == 1 ) return false;
    cutList.nCut_BlobsBad.increment(truth_isSignal, study1, study2);

    // Fill Blob Data
    //fillBlobData();

    if ( applyPi0BlobCuts && isPi0BlobBad() ) return false;
    cutList.nCut_Pi0BlobCuts.increment(truth_isSignal, study1, study2);

    // Fill Gamma1 Conv Length Cut Hist
    cutList.hCut_gamma1ConvDist->Fill(CCProtonPi0_gamma1_dist_vtx * 0.1);
    if (applyPhotonDistance && CCProtonPi0_gamma1_dist_vtx * 0.1 < minPhotonDistance) return false;
    cutList.nCut_Photon1DistanceLow.increment(truth_isSignal, study1, study2);

    // Fill Gamma2 Conv Length Cut Hist
    cutList.hCut_gamma2ConvDist->Fill(CCProtonPi0_gamma2_dist_vtx * 0.1);
    if (applyPhotonDistance && CCProtonPi0_gamma2_dist_vtx * 0.1 < minPhotonDistance) return false;
    cutList.nCut_Photon2DistanceLow.increment(truth_isSignal, study1, study2);

    // Fill Pi0 Invariant Mass Cut Hist
    cutList.hCut_pi0invMass->Fill(CCProtonPi0_pi0_invMass);
    if( CCProtonPi0_pi0_invMass < min_Pi0_invMass || CCProtonPi0_pi0_invMass > max_Pi0_invMass ) return false;
    cutList.nCut_Pi0_invMass.increment(truth_isSignal, study1, study2);
   
    // ------------------------------------------------------------------------
    // Proton Candidate Related Cuts 
    // ------------------------------------------------------------------------
    // 1 Prong Events Always Satisfies Proton Cuts - No Proton Candidate to apply Cut!
    if (nProngs == 1) cutList.nCut_1Prong_Particle_None.increment(truth_isSignal, study1, study2);
    if (nProngs == 1) cutList.nCut_1Prong_Proton_None.increment(truth_isSignal, study1, study2);
    if (nProngs == 1) cutList.nCut_1Prong_ProtonScore.increment(truth_isSignal, study1, study2);
    if (nProngs == 1) cutList.nCut_1Prong_DeltaInvMass.increment(truth_isSignal, study1, study2);
    
    // 2+ Prong Events Must pass the following Cuts
    if ( nProngs > 1 ){
        if( Cut_Particle_None == 1) return false;
        cutList.nCut_2Prong_Particle_None.increment(truth_isSignal, study1, study2);
        
        if( Cut_Proton_None == 1) return false;
        cutList.nCut_2Prong_Proton_None.increment(truth_isSignal, study1, study2);
        
        // Find Best Proton in Reco
        findRecoProton(); 
       
        // Apply Proton Score to All Proton Candidates
        for( unsigned int i = 0; i < 10 && CCProtonPi0_all_protons_LLRScore[i] != -9.9; i++){
            if ( applyProtonScore ){
                // Use pID Difference for KE < pID_KE_Limit 
                // Use LLR for KE > pID_KE_Limit
                if (CCProtonPi0_all_protons_ekin[i] < pID_KE_Limit ){
                    double pIDDiff = CCProtonPi0_all_protons_protonScore[i] - CCProtonPi0_all_protons_pionScore[i];
                    cutList.hCut_protonScore_pIDDiff->Fill(pIDDiff);
                    if ( pIDDiff < minPIDDiff ) return false;
                }else{
                    cutList.hCut_protonScore_LLR->Fill(CCProtonPi0_all_protons_LLRScore[i]);
                    if ( CCProtonPi0_all_protons_LLRScore[i] < minProtonScore_LLR ) return false;
                }
            }
        }
        cutList.nCut_2Prong_ProtonScore.increment(truth_isSignal, study1, study2);
        
        double delta_invMass = calcDeltaInvariantMass();
        cutList.hCut_deltaInvMass->Fill(delta_invMass);
        if (applyDeltaInvMass){
            if ( delta_invMass < min_Delta_invMass || delta_invMass > max_Delta_invMass) return false;  
        }
        cutList.nCut_2Prong_DeltaInvMass.increment(truth_isSignal, study1, study2);
    }

    //-------------------------------------------------------------------------
    // Different Statistics for Different Topologies
    //      Muon + Pi0
    //      Muon + Pi0 + X (No Meson)
    //-------------------------------------------------------------------------
    // Fill Neutrino Energy Cut Hists
    if(nProngs == 1) cutList.hCut_1Prong_neutrinoE->Fill(CCProtonPi0_neutrino_E_Cal * HEP_Functions::MeV_to_GeV);
    if(nProngs >= 2) cutList.hCut_2Prong_neutrinoE->Fill(CCProtonPi0_neutrino_E_Cal * HEP_Functions::MeV_to_GeV); 
    if( applyBeamEnergy && ((CCProtonPi0_neutrino_E_Cal * HEP_Functions::MeV_to_GeV) > max_beamEnergy)) return false;
    if(nProngs == 1) cutList.nCut_1Prong_beamEnergy.increment(truth_isSignal, study1, study2);
    if(nProngs >= 2) cutList.nCut_2Prong_beamEnergy.increment(truth_isSignal, study1, study2);
   
    // Fill Unused Energy Cut Hists 
    if(nProngs == 1) cutList.hCut_1Prong_UnusedE->Fill(energyUnused_afterReco);
    if(nProngs >= 2) cutList.hCut_2Prong_UnusedE->Fill(energyUnused_afterReco);
    if( applyUnusedE && (energyUnused_afterReco > maxUnusedE)) return false;
    if(nProngs == 1) cutList.nCut_1Prong_UnusedE.increment(truth_isSignal, study1, study2);
    if(nProngs >= 2) cutList.nCut_2Prong_UnusedE.increment(truth_isSignal, study1, study2);
        
    return true;
    
}

bool CCProtonPi0_Analyzer::isPi0BlobBad()
{
    // NDigits
    if (gamma1_blob_ndigits < 7 || gamma2_blob_ndigits < 5 ) return true;
    
    // Blob Energy
    //if (gamma1_blob_energy < 60) return true;

    // Fit Results
    if ( blob_ndof_1 < 5) return true;

    // dEdx
    if ( g1dedx_nplane < 3 ) return true;
    if ( g1dedx >= 35) return true;

    return false;
}

void CCProtonPi0_Analyzer::fill_mc_w()
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
bool CCProtonPi0_Analyzer::KeepMichelEvent()
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

void CCProtonPi0_Analyzer::fillInteractionTrue()
{
    // Do Nothing!
}

void CCProtonPi0_Analyzer::fillInteractionReco()
{
    // Event Kinematics
    if (nProngs == 1) interaction.Enu_1Track->Fill(CCProtonPi0_neutrino_E_1Track * HEP_Functions::MeV_to_GeV);
    else interaction.Enu_2Track->Fill(CCProtonPi0_neutrino_E_2Track * HEP_Functions::MeV_to_GeV);
    interaction.Enu_Cal->Fill(CCProtonPi0_neutrino_E_Cal * HEP_Functions::MeV_to_GeV);
    interaction.q2->Fill(CCProtonPi0_QSq_Cal * HEP_Functions::MeVSq_to_GeVSq);
    interaction.w->Fill(CCProtonPi0_W_Cal* HEP_Functions::MeV_to_GeV);
    interaction.wSq->Fill(CCProtonPi0_WSq_Cal * HEP_Functions::MeVSq_to_GeVSq);
   
    // Reconstruction 
    interaction.E_Unused_afterReco->Fill(energyUnused_afterReco);
    interaction.E_Used_afterReco->Fill(energyUsed_afterReco);

    // Other Event Parameters
    if (nProngs >= 2) interaction.deltaInvMass->Fill(calcDeltaInvariantMass());
    interaction.nProngs_hist->Fill(nProngs);
    
    // MC Only
    if(truth_isSignal){
        if(mc_intType == 1) interaction.final_mc_w_CCQE->Fill(mc_w * HEP_Functions::MeV_to_GeV);
        if(mc_intType == 2) interaction.final_mc_w_RES->Fill(mc_w * HEP_Functions::MeV_to_GeV);
        if(mc_intType == 3) interaction.final_mc_w_DIS->Fill(mc_w * HEP_Functions::MeV_to_GeV);
    }
}

double CCProtonPi0_Analyzer::calcDeltaInvariantMass()
{
    double invMassSq;
        
    invMassSq = (CCProtonPi0_pi0_E + CCProtonPi0_proton_E) * (CCProtonPi0_pi0_E + CCProtonPi0_proton_E) -
                ((CCProtonPi0_pi0_px + CCProtonPi0_proton_px)*(CCProtonPi0_pi0_px + CCProtonPi0_proton_px) + 
                 (CCProtonPi0_pi0_py + CCProtonPi0_proton_py)*(CCProtonPi0_pi0_py + CCProtonPi0_proton_py) +
                 (CCProtonPi0_pi0_pz + CCProtonPi0_proton_pz)*(CCProtonPi0_pi0_pz + CCProtonPi0_proton_pz));
        
    return sqrt(invMassSq);
}

void CCProtonPi0_Analyzer::writeScanFile()
{
    if(isScanRun){
        // Constants for Roundup List
        const string arachne_html = "http://minerva05.fnal.gov/Arachne/arachne.html?filename=";
        const string entryString  = "&entry=";
        const string other        = "&slice=-1&filetype=dst";
        
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
        failText[i].close();
    }
}

void CCProtonPi0_Analyzer::openTextFiles()
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

void CCProtonPi0_Analyzer::fillProtonTrue()
{   
    double reco_P = CCProtonPi0_proton_P * HEP_Functions::MeV_to_GeV;
    double true_P = CCProtonPi0_trajProtonProngMomentum[indRecoProton] * HEP_Functions::MeV_to_GeV;
    double error_P = (reco_P - true_P) / true_P;

    proton.reco_P_true_P->Fill(reco_P,true_P);
    proton.P_error->Fill(error_P);
}

void CCProtonPi0_Analyzer::fillProtonReco()
{  
    // Unique Histograms
    proton.partScore->Fill(CCProtonPi0_proton_LLRScore);
    proton.trackLength->Fill(CCProtonPi0_proton_length * HEP_Functions::mm_to_cm);
    proton.trackKinked->Fill(CCProtonPi0_proton_kinked);
 
    // Standard Histograms
    proton.E->Fill(CCProtonPi0_proton_E * HEP_Functions::MeV_to_GeV);
    proton.P->Fill(CCProtonPi0_proton_P * HEP_Functions::MeV_to_GeV);
    proton.KE->Fill(CCProtonPi0_proton_ekin * HEP_Functions::MeV_to_GeV);
    proton.theta->Fill(CCProtonPi0_proton_theta * TMath::RadToDeg());
    proton.phi->Fill(CCProtonPi0_proton_phi * TMath::RadToDeg());
}


//--------------------------------------------------------------------------
//  findRecoProton()
//      finds Reco Proton using Proton Score
//--------------------------------------------------------------------------
void CCProtonPi0_Analyzer::findRecoProton()
{
    double temp_p = CCProtonPi0_all_protons_P[0];
    double current_p;
    int tempInd = 0;
   
    for( int i = 0; i < 10; i++){
        current_p = CCProtonPi0_all_protons_P[i];
        //if( current_p == SENTINEL ) break;
        if( current_p > temp_p){
            temp_p = current_p;
            tempInd = i;
        }
    }
    
    indRecoProton = tempInd;
}



void CCProtonPi0_Analyzer::fillBlobData()
{
    // Set Gamma Blob Information
    pi0Blob.g1_blob_energy->Fill(gamma1_blob_energy);
    pi0Blob.g1_blob_minsep->Fill(gamma1_blob_minsep);
    pi0Blob.g1_blob_nclusters->Fill(gamma1_blob_nclusters);
    pi0Blob.g1_blob_ndigits->Fill(gamma1_blob_ndigits);
    
    pi0Blob.g2_blob_energy->Fill(gamma2_blob_energy);
    pi0Blob.g2_blob_minsep->Fill(gamma2_blob_minsep);
    pi0Blob.g2_blob_nclusters->Fill(gamma2_blob_nclusters);
    pi0Blob.g2_blob_ndigits->Fill(gamma2_blob_ndigits);

    pi0Blob.g1_blob_ndof->Fill(blob_ndof_1);
    pi0Blob.g1_blob_fval->Fill(blob_fval_1);
    pi0Blob.g1_blob_dEdx_doublet->Fill(g1dedx_doublet);
    pi0Blob.g1_blob_dEdx_empty_plane->Fill(g1dedx_empty_plane);
    pi0Blob.g1_blob_dEdx->Fill(g1dedx);
    pi0Blob.g1_blob_dEdx1->Fill(g1dedx1);
    pi0Blob.g1_blob_dEdx_nplane->Fill(g1dedx_nplane);
    //pi0Blob.g1_blob_dEdx_cluster_energy->Fill(g1dedx_cluster_energy);
  
    //cout<<"blob_ndof_1 = "<<blob_ndof_1<<endl;
    //cout<<"blob_fval_1 = "<<blob_fval_1<<endl;
    //cout<<"g1dedx_doublet = "<<g1dedx_doublet<<endl;
    //cout<<"g1dedx_empty_plane = "<<g1dedx_empty_plane<<endl;
    //cout<<"g1dedx = "<<g1dedx<<endl;
    //cout<<"g1dedx1 = "<<g1dedx1<<endl;
    //cout<<"g1dedx_nplane = "<<g1dedx_nplane<<endl;
    //cout<<"-----"<<endl;

}

void CCProtonPi0_Analyzer::fillPi0True()
{
  // Do Nothing 
}

void CCProtonPi0_Analyzer::fillPi0Reco()
{
    // Unique Histograms
    pi0.invMass->Fill(CCProtonPi0_pi0_invMass);
    pi0.gamma1_ConvLength->Fill(CCProtonPi0_gamma1_dist_vtx * 0.1);
    pi0.gamma2_ConvLength->Fill(CCProtonPi0_gamma2_dist_vtx * 0.1);
    pi0.ConvLength_gamma2_gamma1->Fill(CCProtonPi0_gamma2_dist_vtx, CCProtonPi0_gamma1_dist_vtx);
    
    double photon_E_asym = abs((CCProtonPi0_gamma1_E - CCProtonPi0_gamma2_E) / (CCProtonPi0_gamma1_E + CCProtonPi0_gamma2_E));  
    pi0.photonEnergy_Asymmetry->Fill(photon_E_asym);
    
    // Standard Histograms
    pi0.E->Fill(CCProtonPi0_pi0_E * HEP_Functions::MeV_to_GeV);
    pi0.P->Fill(CCProtonPi0_pi0_P * HEP_Functions::MeV_to_GeV);
    pi0.KE->Fill(0 * HEP_Functions::MeV_to_GeV);
    pi0.theta->Fill(CCProtonPi0_pi0_theta);
    pi0.phi->Fill(CCProtonPi0_pi0_phi);
        
}
void CCProtonPi0_Analyzer::fillMuonTrue()
{
    // Do Nothing
}

void CCProtonPi0_Analyzer::fillMuonReco()
{
    muon.E->Fill(CCProtonPi0_muon_E * 0.001);
    muon.P->Fill(CCProtonPi0_muon_P * 0.001);
    muon.KE->Fill(0);
    muon.theta->Fill(CCProtonPi0_muon_theta * TMath::RadToDeg());
    muon.phi->Fill(0);
}


void CCProtonPi0_Analyzer::Init(string playlist, TChain* fChain)
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
        //cout<<" Added "<<filename.c_str()<<endl;
    }

    // Set object pointer
    prong_part_E = 0;
    prong_part_pos = 0;

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
   fChain->SetBranchAddress("Cut_BlobsBad", &Cut_BlobsBad, &b_Cut_BlobsBad);
   fChain->SetBranchAddress("Cut_ConeBlobs", &Cut_ConeBlobs, &b_Cut_ConeBlobs);
   fChain->SetBranchAddress("Cut_EndPoint_Michel_Exist", &Cut_EndPoint_Michel_Exist, &b_Cut_EndPoint_Michel_Exist);
   fChain->SetBranchAddress("Cut_Muon_Charge", &Cut_Muon_Charge, &b_Cut_Muon_Charge);
   fChain->SetBranchAddress("Cut_Muon_None", &Cut_Muon_None, &b_Cut_Muon_None);
   fChain->SetBranchAddress("Cut_Muon_Not_Plausible", &Cut_Muon_Not_Plausible, &b_Cut_Muon_Not_Plausible);
   fChain->SetBranchAddress("Cut_Particle_None", &Cut_Particle_None, &b_Cut_Particle_None);
   fChain->SetBranchAddress("Cut_PreFilter_Pi0", &Cut_PreFilter_Pi0, &b_Cut_PreFilter_Pi0);
   fChain->SetBranchAddress("Cut_Proton_None", &Cut_Proton_None, &b_Cut_Proton_None);
   fChain->SetBranchAddress("Cut_Vertex_Michel_Exist", &Cut_Vertex_Michel_Exist, &b_Cut_Vertex_Michel_Exist);
   fChain->SetBranchAddress("Cut_Vertex_None", &Cut_Vertex_None, &b_Cut_Vertex_None);
   fChain->SetBranchAddress("Cut_Vertex_Not_Fiducial", &Cut_Vertex_Not_Fiducial, &b_Cut_Vertex_Not_Fiducial);
   fChain->SetBranchAddress("Cut_Vertex_Not_Reconstructable", &Cut_Vertex_Not_Reconstructable, &b_Cut_Vertex_Not_Reconstructable);
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
   fChain->SetBranchAddress("nTracks", &nTracks, &b_nTracks);
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
   fChain->SetBranchAddress("vtx_fit_converged", &vtx_fit_converged, &b_vtx_fit_converged);
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
   fChain->SetBranchAddress("g1dedx", &g1dedx, &b_g1dedx);
   fChain->SetBranchAddress("g1dedx1", &g1dedx1, &b_g1dedx1);
   fChain->SetBranchAddress("g1dedx_total", &g1dedx_total, &b_g1dedx_total);
   fChain->SetBranchAddress("g1dedx_total1", &g1dedx_total1, &b_g1dedx_total1);
   fChain->SetBranchAddress("g2dedx", &g2dedx, &b_g2dedx);
   fChain->SetBranchAddress("g2dedx1", &g2dedx1, &b_g2dedx1);
   fChain->SetBranchAddress("g2dedx_total", &g2dedx_total, &b_g2dedx_total);
   fChain->SetBranchAddress("g2dedx_total1", &g2dedx_total1, &b_g2dedx_total1);
   fChain->SetBranchAddress("gamma1_blob_energy", &gamma1_blob_energy, &b_gamma1_blob_energy);
   fChain->SetBranchAddress("gamma1_blob_minsep", &gamma1_blob_minsep, &b_gamma1_blob_minsep);
   fChain->SetBranchAddress("gamma2_blob_energy", &gamma2_blob_energy, &b_gamma2_blob_energy);
   fChain->SetBranchAddress("gamma2_blob_minsep", &gamma2_blob_minsep, &b_gamma2_blob_minsep);
   fChain->SetBranchAddress("hadronVisibleE", &hadronVisibleE, &b_hadronVisibleE);
   fChain->SetBranchAddress("michelProng_begin_Z", &michelProng_begin_Z, &b_michelProng_begin_Z);
   fChain->SetBranchAddress("michelProng_distance", &michelProng_distance, &b_michelProng_distance);
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
   fChain->SetBranchAddress("preFilter_rejectedEnergy", &preFilter_rejectedEnergy, &b_preFilter_rejectedEnergy);
   fChain->SetBranchAddress("prim_vtx_smallest_opening_angle", &prim_vtx_smallest_opening_angle, &b_prim_vtx_smallest_opening_angle);
   fChain->SetBranchAddress("reco_eventID", &reco_eventID, &b_reco_eventID);
   fChain->SetBranchAddress("time", &time, &b_time);
   fChain->SetBranchAddress("totalIDVisibleE", &totalIDVisibleE, &b_totalIDVisibleE);
   fChain->SetBranchAddress("totalODVisibleE", &totalODVisibleE, &b_totalODVisibleE);
   fChain->SetBranchAddress("totalVisibleE", &totalVisibleE, &b_totalVisibleE);
   fChain->SetBranchAddress("vtx_fit_chi2", &vtx_fit_chi2, &b_vtx_fit_chi2);
   fChain->SetBranchAddress("g1dedx_cluster_occupancy_sz", &g1dedx_cluster_occupancy_sz, &b_g1dedx_cluster_occupancy_sz);
   fChain->SetBranchAddress("g1dedx_cluster_occupancy", g1dedx_cluster_occupancy, &b_g1dedx_cluster_occupancy);
   fChain->SetBranchAddress("g2dedx_cluster_occupancy_sz", &g2dedx_cluster_occupancy_sz, &b_g2dedx_cluster_occupancy_sz);
   fChain->SetBranchAddress("g2dedx_cluster_occupancy", g2dedx_cluster_occupancy, &b_g2dedx_cluster_occupancy);
   fChain->SetBranchAddress("Vertex_energy_radii_sz", &Vertex_energy_radii_sz, &b_Vertex_energy_radii_sz);
   fChain->SetBranchAddress("Vertex_energy_radii", Vertex_energy_radii, &b_Vertex_energy_radii);
   fChain->SetBranchAddress("fit_vtx", fit_vtx, &b_fit_vtx);
   fChain->SetBranchAddress("g1dedx_cluster_energy_sz", &g1dedx_cluster_energy_sz, &b_g1dedx_cluster_energy_sz);
   fChain->SetBranchAddress("g1dedx_cluster_energy", g1dedx_cluster_energy, &b_g1dedx_cluster_energy);
   fChain->SetBranchAddress("g1dedx_rev_cluster_energy_sz", &g1dedx_rev_cluster_energy_sz, &b_g1dedx_rev_cluster_energy_sz);
   fChain->SetBranchAddress("g1dedx_rev_cluster_energy", g1dedx_rev_cluster_energy, &b_g1dedx_rev_cluster_energy);
   fChain->SetBranchAddress("g2dedx_cluster_energy_sz", &g2dedx_cluster_energy_sz, &b_g2dedx_cluster_energy_sz);
   fChain->SetBranchAddress("g2dedx_cluster_energy", g2dedx_cluster_energy, &b_g2dedx_cluster_energy);
   fChain->SetBranchAddress("g2dedx_rev_cluster_energy_sz", &g2dedx_rev_cluster_energy_sz, &b_g2dedx_rev_cluster_energy_sz);
   fChain->SetBranchAddress("g2dedx_rev_cluster_energy", g2dedx_rev_cluster_energy, &b_g2dedx_rev_cluster_energy);
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
   fChain->SetBranchAddress("truth_gamma1_4P", truth_gamma1_4P, &b_truth_gamma1_4P);
   fChain->SetBranchAddress("truth_gamma2_4P", truth_gamma2_4P, &b_truth_gamma2_4P);
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
   fChain->SetBranchAddress("truth_muon_4P", truth_muon_4P, &b_truth_muon_4P);
   fChain->SetBranchAddress("truth_pi0_4P", truth_pi0_4P, &b_truth_pi0_4P);
   fChain->SetBranchAddress("truth_proton_4P", truth_proton_4P, &b_truth_proton_4P);
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
   fChain->SetBranchAddress("CCProtonPi0_proton_kinked", &CCProtonPi0_proton_kinked, &b_CCProtonPi0_proton_kinked);
   fChain->SetBranchAddress("CCProtonPi0_r_minos_trk_vtx_plane", &CCProtonPi0_r_minos_trk_vtx_plane, &b_CCProtonPi0_r_minos_trk_vtx_plane);
   fChain->SetBranchAddress("CCProtonPi0_t_minos_trk_numFSMuons", &CCProtonPi0_t_minos_trk_numFSMuons, &b_CCProtonPi0_t_minos_trk_numFSMuons);
   fChain->SetBranchAddress("CCProtonPi0_t_minos_trk_primFSLeptonPDG", &CCProtonPi0_t_minos_trk_primFSLeptonPDG, &b_CCProtonPi0_t_minos_trk_primFSLeptonPDG);
   fChain->SetBranchAddress("CCProtonPi0_trajMuonProngPDG", &CCProtonPi0_trajMuonProngPDG, &b_CCProtonPi0_trajMuonProngPDG);
   fChain->SetBranchAddress("CCProtonPi0_trajMuonProngPrimary", &CCProtonPi0_trajMuonProngPrimary, &b_CCProtonPi0_trajMuonProngPrimary);
   fChain->SetBranchAddress("CCProtonPi0_vtx_module", &CCProtonPi0_vtx_module, &b_CCProtonPi0_vtx_module);
   fChain->SetBranchAddress("CCProtonPi0_vtx_plane", &CCProtonPi0_vtx_plane, &b_CCProtonPi0_vtx_plane);
   fChain->SetBranchAddress("CCProtonPi0_QSq_1Track", &CCProtonPi0_QSq_1Track, &b_CCProtonPi0_QSq_1Track);
   fChain->SetBranchAddress("CCProtonPi0_QSq_Cal", &CCProtonPi0_QSq_Cal, &b_CCProtonPi0_QSq_Cal);
   fChain->SetBranchAddress("CCProtonPi0_WSq_1Track", &CCProtonPi0_WSq_1Track, &b_CCProtonPi0_WSq_1Track);
   fChain->SetBranchAddress("CCProtonPi0_WSq_Cal", &CCProtonPi0_WSq_Cal, &b_CCProtonPi0_WSq_Cal);
   fChain->SetBranchAddress("CCProtonPi0_W_1Track", &CCProtonPi0_W_1Track, &b_CCProtonPi0_W_1Track);
   fChain->SetBranchAddress("CCProtonPi0_W_Cal", &CCProtonPi0_W_Cal, &b_CCProtonPi0_W_Cal);
   fChain->SetBranchAddress("CCProtonPi0_endMuonTrajMomentum", &CCProtonPi0_endMuonTrajMomentum, &b_CCProtonPi0_endMuonTrajMomentum);
   fChain->SetBranchAddress("CCProtonPi0_endMuonTrajXPosition", &CCProtonPi0_endMuonTrajXPosition, &b_CCProtonPi0_endMuonTrajXPosition);
   fChain->SetBranchAddress("CCProtonPi0_endMuonTrajYPosition", &CCProtonPi0_endMuonTrajYPosition, &b_CCProtonPi0_endMuonTrajYPosition);
   fChain->SetBranchAddress("CCProtonPi0_endMuonTrajZPosition", &CCProtonPi0_endMuonTrajZPosition, &b_CCProtonPi0_endMuonTrajZPosition);
   fChain->SetBranchAddress("CCProtonPi0_gamma1_E", &CCProtonPi0_gamma1_E, &b_CCProtonPi0_gamma1_E);
   fChain->SetBranchAddress("CCProtonPi0_gamma1_P", &CCProtonPi0_gamma1_P, &b_CCProtonPi0_gamma1_P);
   fChain->SetBranchAddress("CCProtonPi0_gamma1_dEdx", &CCProtonPi0_gamma1_dEdx, &b_CCProtonPi0_gamma1_dEdx);
   fChain->SetBranchAddress("CCProtonPi0_gamma1_dist_vtx", &CCProtonPi0_gamma1_dist_vtx, &b_CCProtonPi0_gamma1_dist_vtx);
   fChain->SetBranchAddress("CCProtonPi0_gamma1_evis_ecal", &CCProtonPi0_gamma1_evis_ecal, &b_CCProtonPi0_gamma1_evis_ecal);
   fChain->SetBranchAddress("CCProtonPi0_gamma1_evis_hcal", &CCProtonPi0_gamma1_evis_hcal, &b_CCProtonPi0_gamma1_evis_hcal);
   fChain->SetBranchAddress("CCProtonPi0_gamma1_evis_scal", &CCProtonPi0_gamma1_evis_scal, &b_CCProtonPi0_gamma1_evis_scal);
   fChain->SetBranchAddress("CCProtonPi0_gamma1_evis_trkr", &CCProtonPi0_gamma1_evis_trkr, &b_CCProtonPi0_gamma1_evis_trkr);
   fChain->SetBranchAddress("CCProtonPi0_gamma1_phi", &CCProtonPi0_gamma1_phi, &b_CCProtonPi0_gamma1_phi);
   fChain->SetBranchAddress("CCProtonPi0_gamma1_px", &CCProtonPi0_gamma1_px, &b_CCProtonPi0_gamma1_px);
   fChain->SetBranchAddress("CCProtonPi0_gamma1_py", &CCProtonPi0_gamma1_py, &b_CCProtonPi0_gamma1_py);
   fChain->SetBranchAddress("CCProtonPi0_gamma1_pz", &CCProtonPi0_gamma1_pz, &b_CCProtonPi0_gamma1_pz);
   fChain->SetBranchAddress("CCProtonPi0_gamma1_theta", &CCProtonPi0_gamma1_theta, &b_CCProtonPi0_gamma1_theta);
   fChain->SetBranchAddress("CCProtonPi0_gamma1_time", &CCProtonPi0_gamma1_time, &b_CCProtonPi0_gamma1_time);
   fChain->SetBranchAddress("CCProtonPi0_gamma2_E", &CCProtonPi0_gamma2_E, &b_CCProtonPi0_gamma2_E);
   fChain->SetBranchAddress("CCProtonPi0_gamma2_P", &CCProtonPi0_gamma2_P, &b_CCProtonPi0_gamma2_P);
   fChain->SetBranchAddress("CCProtonPi0_gamma2_dEdx", &CCProtonPi0_gamma2_dEdx, &b_CCProtonPi0_gamma2_dEdx);
   fChain->SetBranchAddress("CCProtonPi0_gamma2_dist_vtx", &CCProtonPi0_gamma2_dist_vtx, &b_CCProtonPi0_gamma2_dist_vtx);
   fChain->SetBranchAddress("CCProtonPi0_gamma2_evis_ecal", &CCProtonPi0_gamma2_evis_ecal, &b_CCProtonPi0_gamma2_evis_ecal);
   fChain->SetBranchAddress("CCProtonPi0_gamma2_evis_hcal", &CCProtonPi0_gamma2_evis_hcal, &b_CCProtonPi0_gamma2_evis_hcal);
   fChain->SetBranchAddress("CCProtonPi0_gamma2_evis_scal", &CCProtonPi0_gamma2_evis_scal, &b_CCProtonPi0_gamma2_evis_scal);
   fChain->SetBranchAddress("CCProtonPi0_gamma2_evis_trkr", &CCProtonPi0_gamma2_evis_trkr, &b_CCProtonPi0_gamma2_evis_trkr);
   fChain->SetBranchAddress("CCProtonPi0_gamma2_phi", &CCProtonPi0_gamma2_phi, &b_CCProtonPi0_gamma2_phi);
   fChain->SetBranchAddress("CCProtonPi0_gamma2_px", &CCProtonPi0_gamma2_px, &b_CCProtonPi0_gamma2_px);
   fChain->SetBranchAddress("CCProtonPi0_gamma2_py", &CCProtonPi0_gamma2_py, &b_CCProtonPi0_gamma2_py);
   fChain->SetBranchAddress("CCProtonPi0_gamma2_pz", &CCProtonPi0_gamma2_pz, &b_CCProtonPi0_gamma2_pz);
   fChain->SetBranchAddress("CCProtonPi0_gamma2_theta", &CCProtonPi0_gamma2_theta, &b_CCProtonPi0_gamma2_theta);
   fChain->SetBranchAddress("CCProtonPi0_gamma2_time", &CCProtonPi0_gamma2_time, &b_CCProtonPi0_gamma2_time);
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
   fChain->SetBranchAddress("CCProtonPi0_muon_P", &CCProtonPi0_muon_P, &b_CCProtonPi0_muon_P);
   fChain->SetBranchAddress("CCProtonPi0_muon_muScore", &CCProtonPi0_muon_muScore, &b_CCProtonPi0_muon_muScore);
   fChain->SetBranchAddress("CCProtonPi0_muon_px", &CCProtonPi0_muon_px, &b_CCProtonPi0_muon_px);
   fChain->SetBranchAddress("CCProtonPi0_muon_py", &CCProtonPi0_muon_py, &b_CCProtonPi0_muon_py);
   fChain->SetBranchAddress("CCProtonPi0_muon_pz", &CCProtonPi0_muon_pz, &b_CCProtonPi0_muon_pz);
   fChain->SetBranchAddress("CCProtonPi0_muon_qp", &CCProtonPi0_muon_qp, &b_CCProtonPi0_muon_qp);
   fChain->SetBranchAddress("CCProtonPi0_muon_qpqpe", &CCProtonPi0_muon_qpqpe, &b_CCProtonPi0_muon_qpqpe);
   fChain->SetBranchAddress("CCProtonPi0_muon_roadUpstreamEnergy", &CCProtonPi0_muon_roadUpstreamEnergy, &b_CCProtonPi0_muon_roadUpstreamEnergy);
   fChain->SetBranchAddress("CCProtonPi0_muon_theta", &CCProtonPi0_muon_theta, &b_CCProtonPi0_muon_theta);
   fChain->SetBranchAddress("CCProtonPi0_muon_theta_biasDown", &CCProtonPi0_muon_theta_biasDown, &b_CCProtonPi0_muon_theta_biasDown);
   fChain->SetBranchAddress("CCProtonPi0_muon_theta_biasUp", &CCProtonPi0_muon_theta_biasUp, &b_CCProtonPi0_muon_theta_biasUp);
   fChain->SetBranchAddress("CCProtonPi0_neutrino_E_1Track", &CCProtonPi0_neutrino_E_1Track, &b_CCProtonPi0_neutrino_E_1Track);
   fChain->SetBranchAddress("CCProtonPi0_neutrino_E_2Track", &CCProtonPi0_neutrino_E_2Track, &b_CCProtonPi0_neutrino_E_2Track);
   fChain->SetBranchAddress("CCProtonPi0_neutrino_E_Cal", &CCProtonPi0_neutrino_E_Cal, &b_CCProtonPi0_neutrino_E_Cal);
   fChain->SetBranchAddress("CCProtonPi0_pi0_E", &CCProtonPi0_pi0_E, &b_CCProtonPi0_pi0_E);
   fChain->SetBranchAddress("CCProtonPi0_pi0_P", &CCProtonPi0_pi0_P, &b_CCProtonPi0_pi0_P);
   fChain->SetBranchAddress("CCProtonPi0_pi0_cos_openingAngle", &CCProtonPi0_pi0_cos_openingAngle, &b_CCProtonPi0_pi0_cos_openingAngle);
   fChain->SetBranchAddress("CCProtonPi0_pi0_invMass", &CCProtonPi0_pi0_invMass, &b_CCProtonPi0_pi0_invMass);
   fChain->SetBranchAddress("CCProtonPi0_pi0_openingAngle", &CCProtonPi0_pi0_openingAngle, &b_CCProtonPi0_pi0_openingAngle);
   fChain->SetBranchAddress("CCProtonPi0_pi0_phi", &CCProtonPi0_pi0_phi, &b_CCProtonPi0_pi0_phi);
   fChain->SetBranchAddress("CCProtonPi0_pi0_px", &CCProtonPi0_pi0_px, &b_CCProtonPi0_pi0_px);
   fChain->SetBranchAddress("CCProtonPi0_pi0_py", &CCProtonPi0_pi0_py, &b_CCProtonPi0_pi0_py);
   fChain->SetBranchAddress("CCProtonPi0_pi0_pz", &CCProtonPi0_pi0_pz, &b_CCProtonPi0_pi0_pz);
   fChain->SetBranchAddress("CCProtonPi0_pi0_theta", &CCProtonPi0_pi0_theta, &b_CCProtonPi0_pi0_theta);
   fChain->SetBranchAddress("CCProtonPi0_pi0_thetaX", &CCProtonPi0_pi0_thetaX, &b_CCProtonPi0_pi0_thetaX);
   fChain->SetBranchAddress("CCProtonPi0_pi0_thetaY", &CCProtonPi0_pi0_thetaY, &b_CCProtonPi0_pi0_thetaY);
   fChain->SetBranchAddress("CCProtonPi0_proton_E", &CCProtonPi0_proton_E, &b_CCProtonPi0_proton_E);
   fChain->SetBranchAddress("CCProtonPi0_proton_LLRScore", &CCProtonPi0_proton_LLRScore, &b_CCProtonPi0_proton_LLRScore);
   fChain->SetBranchAddress("CCProtonPi0_proton_P", &CCProtonPi0_proton_P, &b_CCProtonPi0_proton_P);
   fChain->SetBranchAddress("CCProtonPi0_proton_ekin", &CCProtonPi0_proton_ekin, &b_CCProtonPi0_proton_ekin);
   fChain->SetBranchAddress("CCProtonPi0_proton_length", &CCProtonPi0_proton_length, &b_CCProtonPi0_proton_length);
   fChain->SetBranchAddress("CCProtonPi0_proton_phi", &CCProtonPi0_proton_phi, &b_CCProtonPi0_proton_phi);
   fChain->SetBranchAddress("CCProtonPi0_proton_pionScore", &CCProtonPi0_proton_pionScore, &b_CCProtonPi0_proton_pionScore);
   fChain->SetBranchAddress("CCProtonPi0_proton_protonScore", &CCProtonPi0_proton_protonScore, &b_CCProtonPi0_proton_protonScore);
   fChain->SetBranchAddress("CCProtonPi0_proton_px", &CCProtonPi0_proton_px, &b_CCProtonPi0_proton_px);
   fChain->SetBranchAddress("CCProtonPi0_proton_py", &CCProtonPi0_proton_py, &b_CCProtonPi0_proton_py);
   fChain->SetBranchAddress("CCProtonPi0_proton_pz", &CCProtonPi0_proton_pz, &b_CCProtonPi0_proton_pz);
   fChain->SetBranchAddress("CCProtonPi0_proton_theta", &CCProtonPi0_proton_theta, &b_CCProtonPi0_proton_theta);
   fChain->SetBranchAddress("CCProtonPi0_proton_thetaX", &CCProtonPi0_proton_thetaX, &b_CCProtonPi0_proton_thetaX);
   fChain->SetBranchAddress("CCProtonPi0_proton_thetaY", &CCProtonPi0_proton_thetaY, &b_CCProtonPi0_proton_thetaY);
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
   fChain->SetBranchAddress("CCProtonPi0_all_protons_kinked", CCProtonPi0_all_protons_kinked, &b_CCProtonPi0_all_protons_kinked);
   fChain->SetBranchAddress("CCProtonPi0_all_protons_odMatch", CCProtonPi0_all_protons_odMatch, &b_CCProtonPi0_all_protons_odMatch);
   fChain->SetBranchAddress("CCProtonPi0_isProtonInsideOD", CCProtonPi0_isProtonInsideOD, &b_CCProtonPi0_isProtonInsideOD);
   fChain->SetBranchAddress("CCProtonPi0_ntrajProtonProng", CCProtonPi0_ntrajProtonProng, &b_CCProtonPi0_ntrajProtonProng);
   fChain->SetBranchAddress("CCProtonPi0_trajProtonProngPDG", CCProtonPi0_trajProtonProngPDG, &b_CCProtonPi0_trajProtonProngPDG);
   fChain->SetBranchAddress("CCProtonPi0_trajProtonProngPrimary", CCProtonPi0_trajProtonProngPrimary, &b_CCProtonPi0_trajProtonProngPrimary);
   fChain->SetBranchAddress("CCProtonPi0_all_protons_E", CCProtonPi0_all_protons_E, &b_CCProtonPi0_all_protons_E);
   fChain->SetBranchAddress("CCProtonPi0_all_protons_LLRScore", CCProtonPi0_all_protons_LLRScore, &b_CCProtonPi0_all_protons_LLRScore);
   fChain->SetBranchAddress("CCProtonPi0_all_protons_P", CCProtonPi0_all_protons_P, &b_CCProtonPi0_all_protons_P);
   fChain->SetBranchAddress("CCProtonPi0_all_protons_chi2_ndf", CCProtonPi0_all_protons_chi2_ndf, &b_CCProtonPi0_all_protons_chi2_ndf);
   fChain->SetBranchAddress("CCProtonPi0_all_protons_ekin", CCProtonPi0_all_protons_ekin, &b_CCProtonPi0_all_protons_ekin);
   fChain->SetBranchAddress("CCProtonPi0_all_protons_endPointX", CCProtonPi0_all_protons_endPointX, &b_CCProtonPi0_all_protons_endPointX);
   fChain->SetBranchAddress("CCProtonPi0_all_protons_endPointY", CCProtonPi0_all_protons_endPointY, &b_CCProtonPi0_all_protons_endPointY);
   fChain->SetBranchAddress("CCProtonPi0_all_protons_endPointZ", CCProtonPi0_all_protons_endPointZ, &b_CCProtonPi0_all_protons_endPointZ);
   fChain->SetBranchAddress("CCProtonPi0_all_protons_length", CCProtonPi0_all_protons_length, &b_CCProtonPi0_all_protons_length);
   fChain->SetBranchAddress("CCProtonPi0_all_protons_p_calCorrection", CCProtonPi0_all_protons_p_calCorrection, &b_CCProtonPi0_all_protons_p_calCorrection);
   fChain->SetBranchAddress("CCProtonPi0_all_protons_p_dEdXTool", CCProtonPi0_all_protons_p_dEdXTool, &b_CCProtonPi0_all_protons_p_dEdXTool);
   fChain->SetBranchAddress("CCProtonPi0_all_protons_p_visEnergy", CCProtonPi0_all_protons_p_visEnergy, &b_CCProtonPi0_all_protons_p_visEnergy);
   fChain->SetBranchAddress("CCProtonPi0_all_protons_phi", CCProtonPi0_all_protons_phi, &b_CCProtonPi0_all_protons_phi);
   fChain->SetBranchAddress("CCProtonPi0_all_protons_pionScore", CCProtonPi0_all_protons_pionScore, &b_CCProtonPi0_all_protons_pionScore);
   fChain->SetBranchAddress("CCProtonPi0_all_protons_protonScore", CCProtonPi0_all_protons_protonScore, &b_CCProtonPi0_all_protons_protonScore);
   fChain->SetBranchAddress("CCProtonPi0_all_protons_px", CCProtonPi0_all_protons_px, &b_CCProtonPi0_all_protons_px);
   fChain->SetBranchAddress("CCProtonPi0_all_protons_py", CCProtonPi0_all_protons_py, &b_CCProtonPi0_all_protons_py);
   fChain->SetBranchAddress("CCProtonPi0_all_protons_pz", CCProtonPi0_all_protons_pz, &b_CCProtonPi0_all_protons_pz);
   fChain->SetBranchAddress("CCProtonPi0_all_protons_startPointX", CCProtonPi0_all_protons_startPointX, &b_CCProtonPi0_all_protons_startPointX);
   fChain->SetBranchAddress("CCProtonPi0_all_protons_startPointY", CCProtonPi0_all_protons_startPointY, &b_CCProtonPi0_all_protons_startPointY);
   fChain->SetBranchAddress("CCProtonPi0_all_protons_startPointZ", CCProtonPi0_all_protons_startPointZ, &b_CCProtonPi0_all_protons_startPointZ);
   fChain->SetBranchAddress("CCProtonPi0_all_protons_theta", CCProtonPi0_all_protons_theta, &b_CCProtonPi0_all_protons_theta);
   fChain->SetBranchAddress("CCProtonPi0_all_protons_thetaX", CCProtonPi0_all_protons_thetaX, &b_CCProtonPi0_all_protons_thetaX);
   fChain->SetBranchAddress("CCProtonPi0_all_protons_thetaY", CCProtonPi0_all_protons_thetaY, &b_CCProtonPi0_all_protons_thetaY);
   fChain->SetBranchAddress("CCProtonPi0_endProtonTrajMomentum", CCProtonPi0_endProtonTrajMomentum, &b_CCProtonPi0_endProtonTrajMomentum);
   fChain->SetBranchAddress("CCProtonPi0_endProtonTrajXPosition", CCProtonPi0_endProtonTrajXPosition, &b_CCProtonPi0_endProtonTrajXPosition);
   fChain->SetBranchAddress("CCProtonPi0_endProtonTrajYPosition", CCProtonPi0_endProtonTrajYPosition, &b_CCProtonPi0_endProtonTrajYPosition);
   fChain->SetBranchAddress("CCProtonPi0_endProtonTrajZPosition", CCProtonPi0_endProtonTrajZPosition, &b_CCProtonPi0_endProtonTrajZPosition);
   fChain->SetBranchAddress("CCProtonPi0_gamma1_direction", CCProtonPi0_gamma1_direction, &b_CCProtonPi0_gamma1_direction);
   fChain->SetBranchAddress("CCProtonPi0_gamma1_vertex", CCProtonPi0_gamma1_vertex, &b_CCProtonPi0_gamma1_vertex);
   fChain->SetBranchAddress("CCProtonPi0_gamma2_direction", CCProtonPi0_gamma2_direction, &b_CCProtonPi0_gamma2_direction);
   fChain->SetBranchAddress("CCProtonPi0_gamma2_vertex", CCProtonPi0_gamma2_vertex, &b_CCProtonPi0_gamma2_vertex);
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
   fChain->SetBranchAddress("ppfx_cv_wgt", &ppfx_cv_wgt, &b_ppfx_cv_wgt);
   fChain->SetBranchAddress("mc_wgt_GENIE_sz", &mc_wgt_GENIE_sz, &b_mc_wgt_GENIE_sz);
   fChain->SetBranchAddress("mc_wgt_GENIE", mc_wgt_GENIE, &b_mc_wgt_GENIE);
   fChain->SetBranchAddress("mc_wgt_Flux_Tertiary_sz", &mc_wgt_Flux_Tertiary_sz, &b_mc_wgt_Flux_Tertiary_sz);
   fChain->SetBranchAddress("mc_wgt_Flux_Tertiary", mc_wgt_Flux_Tertiary, &b_mc_wgt_Flux_Tertiary);
   fChain->SetBranchAddress("mc_wgt_Flux_BeamFocus_sz", &mc_wgt_Flux_BeamFocus_sz, &b_mc_wgt_Flux_BeamFocus_sz);
   fChain->SetBranchAddress("mc_wgt_Flux_BeamFocus", mc_wgt_Flux_BeamFocus, &b_mc_wgt_Flux_BeamFocus);
   fChain->SetBranchAddress("mc_wgt_Flux_NA49_sz", &mc_wgt_Flux_NA49_sz, &b_mc_wgt_Flux_NA49_sz);
   fChain->SetBranchAddress("mc_wgt_Flux_NA49", mc_wgt_Flux_NA49, &b_mc_wgt_Flux_NA49);
   fChain->SetBranchAddress("mc_wgt_Norm_sz", &mc_wgt_Norm_sz, &b_mc_wgt_Norm_sz);
   fChain->SetBranchAddress("mc_wgt_Norm", &mc_wgt_Norm, &b_mc_wgt_Norm);
   fChain->SetBranchAddress("mc_wgt_ppfx_MIPPNumiYields_sz", &mc_wgt_ppfx_MIPPNumiYields_sz, &b_mc_wgt_ppfx_MIPPNumiYields_sz);
   fChain->SetBranchAddress("mc_wgt_ppfx_MIPPNumiYields", &mc_wgt_ppfx_MIPPNumiYields, &b_mc_wgt_ppfx_MIPPNumiYields);
   fChain->SetBranchAddress("mc_wgt_ppfx_TargetAttenuation_sz", &mc_wgt_ppfx_TargetAttenuation_sz, &b_mc_wgt_ppfx_TargetAttenuation_sz);
   fChain->SetBranchAddress("mc_wgt_ppfx_TargetAttenuation", &mc_wgt_ppfx_TargetAttenuation, &b_mc_wgt_ppfx_TargetAttenuation);
   fChain->SetBranchAddress("mc_wgt_ppfx_NA49_sz", &mc_wgt_ppfx_NA49_sz, &b_mc_wgt_ppfx_NA49_sz);
   fChain->SetBranchAddress("mc_wgt_ppfx_NA49", &mc_wgt_ppfx_NA49, &b_mc_wgt_ppfx_NA49);
   fChain->SetBranchAddress("mc_wgt_ppfx_MIPPKaonsYields_sz", &mc_wgt_ppfx_MIPPKaonsYields_sz, &b_mc_wgt_ppfx_MIPPKaonsYields_sz);
   fChain->SetBranchAddress("mc_wgt_ppfx_MIPPKaonsYields", &mc_wgt_ppfx_MIPPKaonsYields, &b_mc_wgt_ppfx_MIPPKaonsYields);
   fChain->SetBranchAddress("mc_wgt_ppfx_MIPPThinTarget_sz", &mc_wgt_ppfx_MIPPThinTarget_sz, &b_mc_wgt_ppfx_MIPPThinTarget_sz);
   fChain->SetBranchAddress("mc_wgt_ppfx_MIPPThinTarget", &mc_wgt_ppfx_MIPPThinTarget, &b_mc_wgt_ppfx_MIPPThinTarget);
   fChain->SetBranchAddress("mc_wgt_ppfx_Absorption_sz", &mc_wgt_ppfx_Absorption_sz, &b_mc_wgt_ppfx_Absorption_sz);
   fChain->SetBranchAddress("mc_wgt_ppfx_Absorption", &mc_wgt_ppfx_Absorption, &b_mc_wgt_ppfx_Absorption);
   fChain->SetBranchAddress("mc_wgt_ppfx_Others_sz", &mc_wgt_ppfx_Others_sz, &b_mc_wgt_ppfx_Others_sz);
   fChain->SetBranchAddress("mc_wgt_ppfx_Others", &mc_wgt_ppfx_Others, &b_mc_wgt_ppfx_Others);
   fChain->SetBranchAddress("n_prongs", &n_prongs, &b_n_prongs);
   fChain->SetBranchAddress("prong_nParticles", prong_nParticles, &b_prong_nParticles);
   fChain->SetBranchAddress("prong_part_score", prong_part_score, &b_prong_part_score);
   fChain->SetBranchAddress("prong_part_mass", prong_part_mass, &b_prong_part_mass);
   fChain->SetBranchAddress("prong_part_charge", prong_part_charge, &b_prong_part_charge);
   fChain->SetBranchAddress("prong_part_pid", prong_part_pid, &b_prong_part_pid);
   fChain->SetBranchAddress("prong_part_E", &prong_part_E, &b_prong_part_E);
   fChain->SetBranchAddress("prong_part_pos", &prong_part_pos, &b_prong_part_pos);


}

CCProtonPi0_Analyzer::~CCProtonPi0_Analyzer()
{
    closeTextFiles();
}

Int_t CCProtonPi0_Analyzer::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t CCProtonPi0_Analyzer::LoadTree(Long64_t entry)
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



#endif //CCProtonPi0_Analyzer_cpp
