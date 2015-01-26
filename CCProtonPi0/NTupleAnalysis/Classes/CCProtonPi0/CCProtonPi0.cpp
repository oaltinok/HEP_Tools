/*
    See CCProtonPi0.h header or Class Information
*/

#ifndef CCProtonPi0_cpp
#define CCProtonPi0_cpp

#include "CCProtonPi0.h"

using namespace std;

void CCProtonPi0::specifyRunTime()
{
    applyMaxEvents = false;
    nMaxEvents = 1000000;
    
    // Control Flow
    analyze_NoProtonEvents = true;
    isDataAnalysis  = true;
    is_pID_Studies  = true;
    writeFSParticleMomentum = false;
    
    applyProtonScore = false;
    minProtonScore = 0.3;
    
    applyProtonScore_LLR = true;
    minProtonScore_LLR = 0.0;
    
    applyPionScore = false;
    maxPionScore = 0.3;
    
    applyPIDDiff = true;
    minPIDDiff = 0.45;
    
    applyPhotonDistance = false;
    minPhotonDistance = 150; //mm
    
    applyBeamEnergy = true;
    max_beamEnergy = 10.0; // GeV
    
    applyQSq = false;
    max_QSq = 1.5;
    
    applyUnusedE = true;
    maxUnusedE = 300;
    
    min_Pi0_invMass = 40.0;
    max_Pi0_invMass = 200.0;
    
    SENTINEL = -9.9;

    //Select Branches to Activate
    m_ActivateMC            = true;
    m_ActivateInteraction   = true;
    m_ActivatePi0           = true;
   
}

void CCProtonPi0::run(string playlist)
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
    // Initialize the Analysis Variables and Histograms
    //------------------------------------------------------------------------
    specifyRunTime();
    initVariables();
      
    //------------------------------------------------------------------------
    // Branch Selection for Performance
    //------------------------------------------------------------------------
    fChain->SetBranchStatus("*",0);  // disable all branches
    
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
        fChain->SetBranchStatus("evis_ntgt",1);
        fChain->SetBranchStatus("evis_other",1);
        fChain->SetBranchStatus("energyUnused*",1);
        fChain->SetBranchStatus("energyUsed*",1);
        fChain->SetBranchStatus("UsedClustersTime",1);
        fChain->SetBranchStatus("UnusedClustersTime",1);
        fChain->SetBranchStatus("AllClustersTime",1);
        fChain->SetBranchStatus("CCProtonPi0_total_*",1);
    }
    
    if(m_ActivatePi0){
        fChain->SetBranchStatus("pi0_*",1);
        fChain->SetBranchStatus("gamma*",1);
        fChain->SetBranchStatus("final_blob_nc",1);
    }
    
        
    //------------------------------------------------------------------------
    // Loop over Chain
    //------------------------------------------------------------------------
    Long64_t nbytes = 0, nb = 0;
    
    cout<<"Looping over all entries"<<endl;
    
    Long64_t nentries = fChain->GetEntriesFast();
    
    std::string testOutput;

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
        // Remove Problematic Event
        if (jentry == 1516) continue;
       
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
        // Fill Background Branches
        //----------------------------------------------------------------------
        fillBackgroundBranches();

        //----------------------------------------------------------------------
        // Data Analysis and Other Studies
        //----------------------------------------------------------------------
        if (isDataAnalysis) fillData();
        if (is_pID_Studies) fill_pID_Data();
        if (writeFSParticleMomentum) writeFSParticle4P(jentry);
        
        //----------------------------------------------------------------------
        // Fail Checks
        //----------------------------------------------------------------------

    } // end for-loop
    
    
    //--------------------------------------------------------------------------
    // Write Text Files
    //--------------------------------------------------------------------------
    writeCutTable();
    writeBackgroundTable();
    getPi0Family(); // Pi0 Family Information written inside FailFile
    
    //--------------------------------------------------------------------------
    // Write Root Files
    //--------------------------------------------------------------------------
    write_RootFile();           //CCProtonPi0
    muon.write_RootFile();
    proton.write_RootFile();
    pion.write_RootFile();
    
    
    closeTextFiles();
    
//     if(is_pID_Studies) get_pID_Stats();
    
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
//
//     Specific Functions
//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

void CCProtonPi0::fillBackgroundBranches()
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

void CCProtonPi0::setBackgroundBranch(vector<double>& background, bool hasAntiMuon, bool hasMichel, bool hasPrimaryPi0, bool hasSecondaryPi0 )
{
    bool hasPi0Candidate = hasPrimaryPi0 || hasSecondaryPi0;
    background[0]++;                        // Total
    if(hasAntiMuon) background[1]++;        // AntiMuon
    if(hasMichel) background[2]++;          // Michel
    if(hasPrimaryPi0) background[3]++;      // Primary Pi0
    if(hasSecondaryPi0) background[4]++;    // Secondary Pi0
    if(hasPi0Candidate) background[5]++;    // Pi0 Candidate
    
}

void CCProtonPi0::getPi0Family()
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

void CCProtonPi0::fill_pID_Data()
{
    double protonScore_LLR = CCProtonPi0_protonScore_LLR[indRecoProton];
    double protonScore = CCProtonPi0_protonScore[indRecoProton];
    double pionScore_LLR = CCProtonPi0_pionScore_LLR[indRecoProton];
    double pionScore = CCProtonPi0_pionScore[indRecoProton];
    
    if(CCProtonPi0_trajProtonProngPDG[indRecoProton] == 2212){
        pID_proton_protonScore_LLR->Fill(protonScore_LLR);
        pID_proton_protonScore->Fill(protonScore);
        pID_proton_pionScore_LLR->Fill(pionScore_LLR);
        pID_proton_pionScore->Fill(pionScore);
        pID_proton_pionScore_protonScore->Fill(pionScore,protonScore);
        pID_proton_pIDSum->Fill(protonScore + pionScore);
        pID_proton_pIDDiff->Fill(protonScore - pionScore);
        pID_proton_protonScore_protonScore_LLR->Fill(protonScore,protonScore_LLR);
        pID_proton_pionScore_pionScore_LLR->Fill(pionScore,pionScore_LLR);
    }else if(CCProtonPi0_trajProtonProngPDG[indRecoProton] == 211){
        pID_piplus_protonScore_LLR->Fill(protonScore_LLR);
        pID_piplus_protonScore->Fill(protonScore);
        pID_piplus_pionScore_LLR->Fill(pionScore_LLR);
        pID_piplus_pionScore->Fill(pionScore);
        pID_piplus_pionScore_protonScore->Fill(CCProtonPi0_pionScore[indRecoProton],protonScore);
        pID_piplus_pIDSum->Fill(protonScore + pionScore);
        pID_piplus_pIDDiff->Fill(protonScore - pionScore);
        pID_piplus_protonScore_protonScore_LLR->Fill(protonScore,protonScore_LLR);
        pID_piplus_pionScore_pionScore_LLR->Fill(pionScore,pionScore_LLR);
    }else if(CCProtonPi0_trajProtonProngPDG[indRecoProton] == -211){
        pID_piminus_protonScore_LLR->Fill(protonScore_LLR);
        pID_piminus_protonScore->Fill(protonScore);
        pID_piminus_pionScore_LLR->Fill(pionScore_LLR);
        pID_piminus_pionScore->Fill(pionScore);
        pID_piminus_pionScore_protonScore->Fill(pionScore,protonScore);
        pID_piminus_pIDSum->Fill(protonScore + pionScore);
        pID_piminus_pIDDiff->Fill(protonScore - pionScore);
        pID_piminus_protonScore_protonScore_LLR->Fill(protonScore,protonScore_LLR);
        pID_piminus_pionScore_pionScore_LLR->Fill(pionScore,pionScore_LLR);
    }
}

void CCProtonPi0::fillData()
{          
    // Fill Reconstructed Information
    fillInteractionReco();
    fillMuonReco();
    fillProtonReco();
    fillPionReco();
    
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
bool CCProtonPi0::analyzeEvent()
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
void CCProtonPi0::setAnalysisMode(int nMode)
{  
    if ( nMode == 1) {
        cout<<"----------------------------------------------------------------------"<<endl;
        cout<<"Analysis Mode: Signal - Only Signal Events will be Analysed"<<endl;
        branchDir = Folder_List::signal;
    }else if ( nMode == 2){
        cout<<"----------------------------------------------------------------------"<<endl;
        cout<<"Analysis Mode: Background - Only Background Events will be Analysed"<<endl;
        branchDir = Folder_List::background;
    }else{
        cout<<"----------------------------------------------------------------------"<<endl;
        cout<<"Analysis Mode: All - All Events will be Analysed"<<endl;
        branchDir = Folder_List::allEvents;
    }
    
    cout<<"----------------------------------------------------------------------"<<endl;
    
    anaMode = nMode;
    isAnalysisModeSelected = true;
}

void CCProtonPi0::writeFSParticle4P(Long64_t nEntry)
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
                <<" Score = "<<CCProtonPi0_protonScore[indRecoProton]
                <<endl;
        failText[t]<<"Pi0 4-P = ( "
                <<pi0_px<<", "
                <<pi0_py<<", "
                <<pi0_pz<<", "
                <<pi0_E<<" )"
                <<endl;   
    }
}

bool CCProtonPi0::getCutStatistics()
{
    //==========================================================================
    //
    // Reconstruction Cuts - Basic Selections
    //
    //==========================================================================
    
    // Increment Cut Number for Both Topologies at the Same Time Until nProngs Cut
   
    nCut_All.increment(truth_isSignal); // Count All Events before Cuts

    if( Cut_Vertex_None == 1) return false;
    nCut_Vertex_None.increment(truth_isSignal);
    
    if( Cut_Vertex_Null == 1) return false;
    nCut_Vertex_Null.increment(truth_isSignal);
    
    if( Cut_Vertex_Not_Reconstructable == 1) return false;
    nCut_Vertex_Not_Reconstructable.increment(truth_isSignal);
    
    if( Cut_Vertex_Not_Fiducial == 1) return false;
    nCut_Vertex_Not_Fiducial.increment(truth_isSignal);
    
    if( vtx_total_count > 1) return false;
    nCut_Vertex_Count.increment(truth_isSignal);
    
    if( Cut_nProngs == 1 || nProngs > 2) return false;
    nCut_nProngs.increment(truth_isSignal);
    
    // Separate Topologies
    if( Cut_Muon_None == 1) return false;
    if (nProngs == 1) nCut_1Prong_Muon_None.increment(truth_isSignal);
    if (nProngs == 2) nCut_2Prong_Muon_None.increment(truth_isSignal);
        
    if( Cut_Muon_Not_Plausible == 1) return false;
    if (nProngs == 1) nCut_1Prong_Muon_Not_Plausible.increment(truth_isSignal);
    if (nProngs == 2) nCut_2Prong_Muon_Not_Plausible.increment(truth_isSignal);
        
    if( Cut_Muon_Score_Low == 1) return false;
    if (nProngs == 1) nCut_1Prong_Muon_Score_Low.increment(truth_isSignal);
    if (nProngs == 2) nCut_2Prong_Muon_Score_Low.increment(truth_isSignal);
    
    if( Cut_Muon_Charge == 1) return false;
    if (nProngs == 1) nCut_1Prong_Muon_Charge.increment(truth_isSignal);
    if (nProngs == 2) nCut_2Prong_Muon_Charge.increment(truth_isSignal);

    if( Cut_Vertex_Michel_Exist == 1) return false;
    if (nProngs == 1) nCut_1Prong_Vertex_Michel_Exist.increment(truth_isSignal);
    if (nProngs == 2) nCut_2Prong_Vertex_Michel_Exist.increment(truth_isSignal);
    
    if( Cut_EndPoint_Michel_Exist == 1) return false;
    if (nProngs == 1) nCut_1Prong_EndPoint_Michel_Exist.increment(truth_isSignal);
    if (nProngs == 2) nCut_2Prong_EndPoint_Michel_Exist.increment(truth_isSignal);
    
    if( Cut_secEndPoint_Michel_Exist == 1) return false;
    if (nProngs == 1) nCut_1Prong_secEndPoint_Michel_Exist.increment(truth_isSignal);
    if (nProngs == 2) nCut_2Prong_secEndPoint_Michel_Exist.increment(truth_isSignal);

    if( Cut_PreFilter_Pi0 == 1) return false;
    if (nProngs == 1) nCut_1Prong_PreFilter_Pi0.increment(truth_isSignal);
    if (nProngs == 2) nCut_2Prong_PreFilter_Pi0.increment(truth_isSignal);
    
    if( Cut_VtxBlob == 1) return false;
    if (nProngs == 1) nCut_1Prong_VtxBlob.increment(truth_isSignal);
    if (nProngs == 2) nCut_2Prong_VtxBlob.increment(truth_isSignal);
    
    if( Cut_ConeBlobs == 1) return false;
    if (nProngs == 1) nCut_1Prong_ConeBlobs.increment(truth_isSignal);
    if (nProngs == 2) nCut_2Prong_ConeBlobs.increment(truth_isSignal);
    
    // Proton Related Cuts only for nProngs == 2
    if ( nProngs == 2 ){
        
        if( Cut_Particle_None == 1) return false;
        nCut_2Prong_Particle_None.increment(truth_isSignal);
        
        if( Cut_Proton_None == 1) return false;
        nCut_2Prong_Proton_None.increment(truth_isSignal);
        
        // Find Best Proton in Reco
        findRecoProton(); 

        if ( applyProtonScore && (CCProtonPi0_protonScore[indRecoProton] < minProtonScore) ) return false;
        nCut_2Prong_Proton_Score.increment(truth_isSignal);
        
        if ( applyPionScore && (CCProtonPi0_pionScore[indRecoProton] > maxPionScore) ) return false;
        nCut_2Prong_Pion_Score.increment(truth_isSignal);
        
        hCut_pIDDiff->Fill(CCProtonPi0_protonScore[indRecoProton] - CCProtonPi0_pionScore[indRecoProton]);
        if ( applyPIDDiff && (CCProtonPi0_protonScore[indRecoProton] - CCProtonPi0_pionScore[indRecoProton]) < minPIDDiff ) return false;
        nCut_2Prong_pIDDiff.increment(truth_isSignal);
        
        hCut_protonScore_LLR->Fill(CCProtonPi0_protonScore_LLR[indRecoProton]);
        if ( applyProtonScore_LLR && (CCProtonPi0_protonScore_LLR[indRecoProton] < minProtonScore_LLR) ) return false;
        nCut_2Prong_Proton_Score_LLR.increment(truth_isSignal);

    }
    
    
    //==========================================================================
    //
    // Analysis Cuts - Improves Purity Dramatically
    //
    //==========================================================================
    if( pi0_invMass < min_Pi0_invMass || pi0_invMass > max_Pi0_invMass ) return false;
    if(nProngs == 1) nCut_1Prong_Pi0_invMass.increment(truth_isSignal);
    if(nProngs == 2) nCut_2Prong_Pi0_invMass.increment(truth_isSignal);
        
    if(applyPhotonDistance && isPhotonDistanceLow()) return false;
    if(nProngs == 1) nCut_1Prong_PhotonDistanceLow.increment(truth_isSignal);
    if(nProngs == 2) nCut_2Prong_PhotonDistanceLow.increment(truth_isSignal);
    
    if( applyBeamEnergy && ((CCProtonPi0_neutrino_E * HEP_Functions::MeV_to_GeV) > max_beamEnergy)) return false;
    if(nProngs == 1) nCut_1Prong_beamEnergy.increment(truth_isSignal);
    if(nProngs == 2) nCut_2Prong_beamEnergy.increment(truth_isSignal);
        
    if( applyUnusedE && (energyUnused_postProton > maxUnusedE)) return false;
    if(nProngs == 1) nCut_1Prong_UnusedE.increment(truth_isSignal);
    if(nProngs == 2) nCut_2Prong_UnusedE.increment(truth_isSignal);
    
    return true;
    
}

// bool CCProtonPi0::getCutStatistics()
// {
//     //==========================================================================
//     //
//     // Reconstruction Cuts - Basic Selections
//     //
//     //==========================================================================
//     
//     
//     nCut_All[t].increment(truth_isSignal); // Count All Events before Cuts
//   
//     if( Cut_Vertex_None == 1) return false;
//     nCut_Vertex_None.increment(truth_isSignal);
//     
//     if( Cut_Vertex_Null == 1) return false;
//     nCut_Vertex_Null.increment(truth_isSignal);
//     
//     if( Cut_Vertex_Not_Reconstructable == 1) return false;
//     nCut_Vertex_Not_Reconstructable.increment(truth_isSignal);
//     
//     if( Cut_Vertex_Not_Fiducial == 1) return false;
//     nCut_Vertex_Not_Fiducial.increment(truth_isSignal);
//     
//     hCut_vertexCount->Fill(vtx_total_count);
//     if( vtx_total_count > 1) return false;
//     nCut_Vertex_Count.increment(truth_isSignal);
//         
//     if( Cut_Muon_None == 1) return false;
//     nCut_Muon_None.increment(truth_isSignal);
//          
//     if( Cut_Muon_Not_Plausible == 1) return false;
//     nCut_Muon_Not_Plausible.increment(truth_isSignal);
//     
//     //--------------------------------------------------------------------------
//     //  Fill MC Truth
//     //--------------------------------------------------------------------------
//         // Fill mc_intType
//         if(truth_isSignal){
//             if(mc_intType == 1) mc_w_CCQE->Fill(mc_w * HEP_Functions::MeV_to_GeV);
//             if(mc_intType == 2) mc_w_RES->Fill(mc_w * HEP_Functions::MeV_to_GeV);
//             if(mc_intType == 3) mc_w_DIS->Fill(mc_w * HEP_Functions::MeV_to_GeV);
//             status_Pi0->Fill(truth_pi0_status);
//             status_Pi0_Mother->Fill(truth_pi0_MotherStatus);
//             status_Pi0_GrandMother->Fill(truth_pi0_GrandMotherStatus);
//             PDG_pi0_Mother.push_back(truth_pi0_Mother);
//             PDG_pi0_GrandMother.push_back(truth_pi0_GrandMother);
//         }
//         
//     //--------------------------------------------------------------------------
//     
//     if( Cut_Muon_Score_Low == 1) return false;
//     nCut_Muon_Score_Low.increment(truth_isSignal);
//     
//     if( Cut_Muon_Charge == 1) return false;
//     nCut_Muon_Charge.increment(truth_isSignal);
//     
//     // Fill hCut_nProngs
//     hCut_nProngs->Fill(nProngs);
//     if( Cut_nProngs == 1 || nProngs > 2) return false;
//     nCut_nProngs.increment(truth_isSignal);
//     
//     // Fill  hCut_Michel
//     if( Cut_Vertex_Michel_Exist == 1 || Cut_EndPoint_Michel_Exist == 1 || Cut_secEndPoint_Michel_Exist == 1 ){
//         hCut_Michel->Fill(1);
//     }else{
//         hCut_Michel->Fill(0);
//     }
//     if( Cut_Vertex_Michel_Exist == 1) return false;
//     nCut_Vertex_Michel_Exist.increment(truth_isSignal);
//     
//     if( Cut_EndPoint_Michel_Exist == 1) return false;
//     nCut_EndPoint_Michel_Exist.increment(truth_isSignal);
//   
//     if( Cut_secEndPoint_Michel_Exist == 1) return false;
//     nCut_secEndPoint_Michel_Exist.increment(truth_isSignal);
//     
//     //--------------------------------------------------------------------------
//     // fill Pre Filter   
//         pFilter_Status->Fill(preFilter_Result);
//         pFilter_RejectedEnergy->Fill(preFilter_rejectedEnergy);
//         hCut_eVis_nuclearTarget->Fill(evis_ntgt);
//         hCut_eVis_other->Fill(evis_other);
//     //--------------------------------------------------------------------------
//     if( Cut_PreFilter_Pi0 == 1) return false;
//     nCut_PreFilter_Pi0.increment(truth_isSignal);
//     
//     if( Cut_VtxBlob == 1) return false;
//     nCut_VtxBlob.increment(truth_isSignal);
//     
//     if( Cut_ConeBlobs == 1) return false;
//     nCut_ConeBlobs.increment(truth_isSignal);
//     
//     
//     if( Cut_Particle_None == 1 && !analyze_NoProtonEvents) return false;
//     nCut_Particle_None.increment(truth_isSignal);
//     
//     if( Cut_Proton_None == 1 && !analyze_NoProtonEvents) return false;
//     nCut_Proton_None.increment(truth_isSignal);
//     
//     if ( Cut_Particle_None == 1 || Cut_Proton_None == 1 || CCProtonPi0_protonScore[0] == -9.9  ){
//         nCut_Proton_Score.increment(truth_isSignal);
//         nCut_Pion_Score.increment(truth_isSignal);
//         nCut_pIDDiff.increment(truth_isSignal);
//         
//     }else{
//         
//         // Find Best Proton in Reco
//         findRecoProton(); 
//         
//         // Fill  hCut pID Histograms
//         hCut_protonScore->Fill(CCProtonPi0_protonScore[indRecoProton]);
//         hCut_pionScore->Fill(CCProtonPi0_pionScore[indRecoProton]);
//         hCut_pIDDiff->Fill(CCProtonPi0_protonScore[indRecoProton] - CCProtonPi0_pionScore[indRecoProton]);
//         
//         if ( applyProtonScore && (CCProtonPi0_protonScore[indRecoProton] < minProtonScore) ) return false;
//         nCut_Proton_Score.increment(truth_isSignal);
//         
//         if ( applyPionScore && (CCProtonPi0_pionScore[indRecoProton] > maxPionScore) ) return false;
//         nCut_Pion_Score.increment(truth_isSignal);
//         
//         if ( applyPIDDiff && (CCProtonPi0_protonScore[indRecoProton] - CCProtonPi0_pionScore[indRecoProton]) < minPIDDiff ) return false;
//         nCut_pIDDiff.increment(truth_isSignal);
//         
//         
//     }
//     
//     //==========================================================================
//     //
//     // Analysis Cuts - Improves Purity Dramatically
//     //
//     //==========================================================================
//     hCut_pi0invMass->Fill(pi0_invMass);
//     if( pi0_invMass < min_Pi0_invMass || pi0_invMass > max_Pi0_invMass ) return false;
//     nCut_Pi0_invMass.increment(truth_isSignal);
//     
//     hCut_gamma1ConvDist->Fill(gamma1_dist_vtx);
//     hCut_gamma2ConvDist->Fill(gamma2_dist_vtx);
//     if( applyPhotonDistance && isPhotonDistanceLow()) return false;
//     nCut_PhotonDistanceLow.increment(truth_isSignal);
//     
//     hCut_neutrinoE->Fill(CCProtonPi0_neutrino_E * HEP_Functions::MeV_to_GeV);
//     if( applyBeamEnergy && ((CCProtonPi0_neutrino_E * HEP_Functions::MeV_to_GeV) > max_beamEnergy)) return false;
//     nCut_beamEnergy.increment(truth_isSignal);
//     
//     hCut_QSq->Fill(CCProtonPi0_QSq * HEP_Functions::MeVSq_to_GeVSq);
//     if( applyQSq && ((CCProtonPi0_QSq * HEP_Functions::MeVSq_to_GeVSq) > max_QSq)) return false;
//     nCut_QSq.increment(truth_isSignal);
//     
//     hCut_UnusedE->Fill(energyUnused_postProton);
//     if( applyUnusedE && (energyUnused_postProton > maxUnusedE)) return false;
//     nCut_UnusedE.increment(truth_isSignal);
//     
//     return true;
//     
// }

void CCProtonPi0::formBackgroundVectors()
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

void CCProtonPi0::writeBackgroundTableHeader()
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


void CCProtonPi0::writeBackgroundTableRows(vector< vector<double> > bckgVector, int nProngs)
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

void CCProtonPi0::writeBackgroundTable()
{
    formBackgroundVectors();
    writeBackgroundTableHeader();
    writeBackgroundTableRows(bckgVector_1Prong, 1);
    writeBackgroundTableRows(bckgVector_2Prong, 2);
}

void CCProtonPi0::writeCutTableHeader()
{
    for (int i = 0; i < nTopologies; i++){
        cutText[i]<<std::left;
        
        cutText[i].width(35); cutText[i]<<"Cut"<<" | "; 
        
        cutText[i].width(12); cutText[i]<<"N(Events)"<<" | ";    
        cutText[i].width(12); cutText[i]<<"N(Signal)"<<" | ";      
        cutText[i].width(12); cutText[i]<<"Efficiency"<<" | ";      
        cutText[i].width(12); cutText[i]<<"Purity"<<" | ";
        cutText[i]<<endl;
    }

}


double CCProtonPi0::getCutEfficiency(double nSig, double effBase)
{
    double eff;   
    eff = (nSig / effBase) * 100.0;
    return eff;   
}

double CCProtonPi0::getCutPurity(double nSig, double nEvents)
{
    double purity;   
    purity = (nSig / nEvents) * 100.0;
    return purity;   
}

void CCProtonPi0::formCutVectors()
{   
    //--------------------------------------------------------------------------
    // Form nCutVector for nProngs == 1
    //--------------------------------------------------------------------------
    // Complete List
    nCutVector_1Prong.push_back(nCut_All);
    nCutVector_1Prong.push_back(nCut_Vertex_None);
    nCutVector_1Prong.push_back(nCut_Vertex_Null);
    nCutVector_1Prong.push_back(nCut_Vertex_Not_Reconstructable); 
    nCutVector_1Prong.push_back(nCut_Vertex_Not_Fiducial);
    nCutVector_1Prong.push_back(nCut_Vertex_Count); 
    nCutVector_1Prong.push_back(nCut_nProngs);
    nCutVector_1Prong.push_back(nCut_1Prong_Muon_None);              
    nCutVector_1Prong.push_back(nCut_1Prong_Muon_Not_Plausible);
    nCutVector_1Prong.push_back(nCut_1Prong_Muon_Score_Low);
    nCutVector_1Prong.push_back(nCut_1Prong_Muon_Charge);
    nCutVector_1Prong.push_back(nCut_1Prong_Vertex_Michel_Exist); 
    nCutVector_1Prong.push_back(nCut_1Prong_EndPoint_Michel_Exist);
    nCutVector_1Prong.push_back(nCut_1Prong_secEndPoint_Michel_Exist);
    nCutVector_1Prong.push_back(nCut_1Prong_PreFilter_Pi0);
    nCutVector_1Prong.push_back(nCut_1Prong_VtxBlob);
    nCutVector_1Prong.push_back(nCut_1Prong_ConeBlobs);
    nCutVector_1Prong.push_back(nCut_1Prong_Pi0_invMass);
    nCutVector_1Prong.push_back(nCut_1Prong_PhotonDistanceLow);
    nCutVector_1Prong.push_back(nCut_1Prong_beamEnergy);
    nCutVector_1Prong.push_back(nCut_1Prong_UnusedE);
    
    // Short List
    nCutVector_1Prong_ShortList.push_back(nCut_All);
    nCutVector_1Prong_ShortList.push_back(nCut_Vertex_Not_Fiducial);
    nCutVector_1Prong_ShortList.push_back(nCut_nProngs);
    nCutVector_1Prong_ShortList.push_back(nCut_1Prong_Muon_None);
    nCutVector_1Prong_ShortList.push_back(nCut_1Prong_secEndPoint_Michel_Exist);
    nCutVector_1Prong_ShortList.push_back(nCut_1Prong_PreFilter_Pi0);
    nCutVector_1Prong_ShortList.push_back(nCut_1Prong_ConeBlobs);
    nCutVector_1Prong_ShortList.push_back(nCut_1Prong_Pi0_invMass);
    nCutVector_1Prong_ShortList.push_back(nCut_1Prong_beamEnergy);
    nCutVector_1Prong_ShortList.push_back(nCut_1Prong_UnusedE);
    
    //--------------------------------------------------------------------------
    // Form nCutVector for nProngs == 2
    //--------------------------------------------------------------------------
    // Complete List
    nCutVector_2Prong.push_back(nCut_All);
    nCutVector_2Prong.push_back(nCut_Vertex_None);
    nCutVector_2Prong.push_back(nCut_Vertex_Null);
    nCutVector_2Prong.push_back(nCut_Vertex_Not_Reconstructable); 
    nCutVector_2Prong.push_back(nCut_Vertex_Not_Fiducial);
    nCutVector_2Prong.push_back(nCut_Vertex_Count); 
    nCutVector_2Prong.push_back(nCut_nProngs);
    nCutVector_2Prong.push_back(nCut_2Prong_Muon_None);              
    nCutVector_2Prong.push_back(nCut_2Prong_Muon_Not_Plausible);
    nCutVector_2Prong.push_back(nCut_2Prong_Muon_Score_Low);
    nCutVector_2Prong.push_back(nCut_2Prong_Muon_Charge);
    nCutVector_2Prong.push_back(nCut_2Prong_Vertex_Michel_Exist); 
    nCutVector_2Prong.push_back(nCut_2Prong_EndPoint_Michel_Exist);
    nCutVector_2Prong.push_back(nCut_2Prong_secEndPoint_Michel_Exist);
    nCutVector_2Prong.push_back(nCut_2Prong_PreFilter_Pi0);
    nCutVector_2Prong.push_back(nCut_2Prong_VtxBlob);
    nCutVector_2Prong.push_back(nCut_2Prong_ConeBlobs);
    nCutVector_2Prong.push_back(nCut_2Prong_Particle_None);
    nCutVector_2Prong.push_back(nCut_2Prong_Proton_None);
    nCutVector_2Prong.push_back(nCut_2Prong_Proton_Score);
    nCutVector_2Prong.push_back(nCut_2Prong_Pion_Score);
    nCutVector_2Prong.push_back(nCut_2Prong_pIDDiff);
    nCutVector_2Prong.push_back(nCut_2Prong_Proton_Score_LLR);
    nCutVector_2Prong.push_back(nCut_2Prong_Pi0_invMass);
    nCutVector_2Prong.push_back(nCut_2Prong_PhotonDistanceLow);
    nCutVector_2Prong.push_back(nCut_2Prong_beamEnergy);
    nCutVector_2Prong.push_back(nCut_2Prong_UnusedE);
    
    // Short List
    nCutVector_2Prong_ShortList.push_back(nCut_All);
    nCutVector_2Prong_ShortList.push_back(nCut_Vertex_Not_Fiducial);
    nCutVector_2Prong_ShortList.push_back(nCut_nProngs);
    nCutVector_2Prong_ShortList.push_back(nCut_2Prong_Muon_None);
    nCutVector_2Prong_ShortList.push_back(nCut_2Prong_secEndPoint_Michel_Exist);
    nCutVector_2Prong_ShortList.push_back(nCut_2Prong_PreFilter_Pi0);
    nCutVector_2Prong_ShortList.push_back(nCut_2Prong_ConeBlobs);
    nCutVector_2Prong_ShortList.push_back(nCut_2Prong_Proton_None);
    nCutVector_2Prong_ShortList.push_back(nCut_2Prong_pIDDiff);
    nCutVector_2Prong_ShortList.push_back(nCut_2Prong_Proton_Score_LLR);
    nCutVector_2Prong_ShortList.push_back(nCut_2Prong_Pi0_invMass);
    nCutVector_2Prong_ShortList.push_back(nCut_2Prong_beamEnergy);
    nCutVector_2Prong_ShortList.push_back(nCut_2Prong_UnusedE);

}


void CCProtonPi0::writeCutTable()
{
    formCutVectors();
    
    for (int t = 0; t < nTopologies; t++){
        cout<<">> Writing "<<cutFile[t]<<endl;
    }
    
    // Write Short List
    writeCutTableHeader();
    writeCutTableRows(nCutVector_1Prong_ShortList, 1, true);
    writeCutTableRows(nCutVector_2Prong_ShortList, 2, true);
    
    for (int t = 0; t < nTopologies; t++){
        cutText[t]<<endl; cutText[t]<<endl;
    }
    
    // Write Complete List
    writeCutTableHeader();
    writeCutTableRows(nCutVector_1Prong, 1, false);
    writeCutTableRows(nCutVector_2Prong, 2, false);
}

void CCProtonPi0::writeCutTableRows(vector<Cut> nCutVector, int nProngs, bool isShortList)
{
    int t; // Topology
    double efficiency;
    double purity;    
    double efficiencyBase;
    int efficiencyInd;
    
    if ( isShortList ) efficiencyInd = 3;
    else efficiencyInd = 7;
    
    if (nProngs == 1) t = 0;
    if (nProngs == 2) t = 1;

    // Efficiency Base is MINOS Matched Muons
    efficiencyBase = nCutVector[efficiencyInd].nSignal.getCount();
      
    
    // Loop over Each Cut    
    for( unsigned int i = 0; i < nCutVector.size(); i++){

            efficiency = getCutEfficiency(nCutVector[i].nSignal.getCount(),efficiencyBase);
            purity = getCutPurity(nCutVector[i].nSignal.getCount(),nCutVector[i].nEvent.getCount());
            
            cutText[t].unsetf( std::ios::floatfield ); 
            cutText[t].width(35); cutText[t]<<nCutVector[i].get_Name()<<" | ";
            cutText[t].width(12); cutText[t]<<nCutVector[i].nEvent.getCount()<<" | ";
            
            // Total Signal
            cutText[t].width(12); cutText[t]<<nCutVector[i].nSignal.getCount()<<" | ";

            cutText[t].precision(4); 
        
            if ( efficiency <= 100){
                cutText[t].width(12); cutText[t]<<efficiency<<" | ";
            }else{
                cutText[t].width(12); cutText[t]<<"N/A"<<" | ";    
            }

            if ( efficiency <= 100){
                cutText[t].width(12); cutText[t]<<purity<<" | ";
            }else{
                cutText[t].width(12);  cutText[t]<<"N/A"<<" | ";    
            }
            
            cutText[t]<<endl;
    }
     
}

void CCProtonPi0::get_pID_Stats()
{
    cout<<"=== Calculationg pID Statistics ==="<<endl;
    
    string rootDir = "Output/RootFiles/Interaction.root";
    TFile* f_Root = new TFile(rootDir.c_str());
    
    TH1D* pID_proton  = (TH1D*)f_Root->Get("pID_proton_protonScore");
    TH1D* pID_piplus  = (TH1D*)f_Root->Get("pID_piplus_protonScore");
    TH1D* pID_piminus = (TH1D*)f_Root->Get("pID_piminus_protonScore");
    TH1D* pID_other   = (TH1D*)f_Root->Get("pID_other_protonScore");
    
    double nProton = 0;
    double nTotalProton = 0;
    double nCapturedEvents = 0;
    double nEvents = 0;
    double purity;
    double efficiency;
    int nBins = 20;
    

    //Get Total Proton
    for(int i = nBins; i >= 1; i--){
        nTotalProton = nTotalProton + pID_proton->GetBinContent(i);
    }
    cout<<"nTotalProton = "<<nTotalProton<<endl;
    
    
    for(int i = nBins; i >= 1; i--){
        nProton = nProton + pID_proton->GetBinContent(i);
        nCapturedEvents =   nCapturedEvents+
                            pID_proton->GetBinContent(i) +
                            pID_piplus->GetBinContent(i) +
                            pID_piminus->GetBinContent(i) +
                            pID_other->GetBinContent(i);
        nEvents =           pID_proton->GetBinContent(i) +
                            pID_piplus->GetBinContent(i) +
                            pID_piminus->GetBinContent(i) +
                            pID_other->GetBinContent(i);
                            
        purity = nProton / nCapturedEvents;
        efficiency = nProton / nTotalProton;
//         cout<<"pID = "<<pID_proton->GetBinLowEdge(i)<<" Purity = "<<purity<<" Efficiency = "<<efficiency<<endl;
        cout<<pID_proton->GetBinLowEdge(i)<<" "<<purity<<" "<<efficiency<<endl;
    }
    
}


void CCProtonPi0::fillInteractionTrue()
{
    
    int_channel->Fill(mc_intType);

    vertex_z_true->Fill(mc_vtx[2]);
    vertex_z_error->Fill( Data_Functions::getError(mc_vtx[2],CCProtonPi0_vtx[2]) );
    vertex_z_reco_mc->Fill(CCProtonPi0_vtx[2],mc_vtx[2]);
    
    vertex_x_y_true->Fill(mc_vtx[0],mc_vtx[1]);
    
}

void CCProtonPi0::fillInteractionReco()
{
//     double invMass_reco;
    
//     invMass_reco = calcDeltaInvariantMass(  pi0_px, 
//                                             pi0_py, 
//                                             pi0_pz,
//                                             pi0_E,
//                                             CCProtonPi0_proton_px[indProton],
//                                             CCProtonPi0_proton_py[indProton],
//                                             CCProtonPi0_proton_pz[indProton],
//                                             CCProtonPi0_proton_E[indProton]);
//                                             
//                                                 
//     deltaInvMass_reco->Fill(invMass_reco);

    // Debugging
//     if(CCProtonPi0_WSq <= 0){
//         double Tp;
//         const double Mn = 939.57;
//         double term1;
//         
//         if (proton.kineticEnergy[0] > 0) Tp = proton.kineticEnergy[0];
//         else Tp = 0.0;
//         
//         term1 = Mn*Mn + 2*Mn*(pion.p4[0].E()+Tp);
//         
//         cout.width(16); cout<<"P(Muon)  = "<<"( "<<muon.p4[0].Px()<<", "<<muon.p4[0].Py()<<", "<<muon.p4[0].Pz()<<", "<<muon.p4[0].E()<<" )"<<endl;        
//         cout.width(16); cout<<"P(Proton) = "<<"( "<<proton.p4[0].Px()<<", "<<proton.p4[0].Py()<<", "<<proton.p4[0].Pz()<<", "<<proton.p4[0].E()<<" )"<<endl;
//         cout.width(16); cout<<"P(Pi0) = "<<"( "<<pion.p4[0].Px()<<", "<<pion.p4[0].Py()<<", "<<pion.p4[0].Pz()<<", "<<pion.p4[0].E()<<" )"<<endl;
//         cout.width(16); cout<<"True Enu = "<<mc_incomingE* HEP_Functions::MeV_to_GeV<<" Reco Enu = "<<CCProtonPi0_neutrino_E* HEP_Functions::MeV_to_GeV<<endl;
//         cout.width(16); cout<<"True Enu_Cal = "<<mc_incomingE* HEP_Functions::MeV_to_GeV<<" Reco Enu_Cal = "<<CCProtonPi0_neutrino_E_Cal* HEP_Functions::MeV_to_GeV<<endl;
//         cout.width(16); cout<<"term1 = "<<term1<<endl;
//         cout.width(16); cout<<"True QSq = "<<mc_Q2<<" Reco QSq = "<<CCProtonPi0_QSq<<endl;
//         cout.width(16); cout<<"True WSq = "<<mc_w*mc_w<<" Reco WSq = "<<CCProtonPi0_WSq<<endl;
//         cout<<endl;
//         
//         w2fail_q2->Fill(CCProtonPi0_QSq * HEP_Functions::MeVSq_to_GeVSq);
//         w2fail_Enu->Fill(CCProtonPi0_neutrino_E * HEP_Functions::MeV_to_GeV);
//         w2fail_muon_E->Fill(muon.p4[0].E() );
//         w2fail_muon_Pz->Fill(muon.p4[0].Pz());
//         w2fail_pion_E->Fill(pion.p4[0].E() * HEP_Functions::MeV_to_GeV);
//         w2fail_proton_KE->Fill(proton.kineticEnergy[0] * HEP_Functions::MeV_to_GeV);
//         w2fail_term1->Fill(term1 * HEP_Functions::MeVSq_to_GeVSq);
//         w2fail_term1_q2->Fill(term1 * HEP_Functions::MeVSq_to_GeVSq,CCProtonPi0_QSq * HEP_Functions::MeVSq_to_GeVSq );
//         w2fail_w2->Fill((term1-CCProtonPi0_QSq)* HEP_Functions::MeVSq_to_GeVSq );
//     }

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
    
    time_AllClusters->Fill(AllClustersTime);
    
    E_Unused_postProton->Fill(energyUnused_postProton);
    E_Used_postProton->Fill(energyUsed_postProton);
    
    if(truth_isSignal){
        if(mc_intType == 1) final_mc_w_CCQE->Fill(mc_w * HEP_Functions::MeV_to_GeV);
        if(mc_intType == 2) final_mc_w_RES->Fill(mc_w * HEP_Functions::MeV_to_GeV);
        if(mc_intType == 3) final_mc_w_DIS->Fill(mc_w * HEP_Functions::MeV_to_GeV);
    }
}

double CCProtonPi0::calcDeltaInvariantMass(double px1, double py1, double pz1, double E1,
                                           double px2, double py2, double pz2, double E2)
{
    double invMassSq;
    
    invMassSq = (E1 + E2) * (E1 + E2) -
                ((px1 + px2)*(px1 + px2) + 
                 (py1 + py2)*(py1 + py2) +
                 (pz1 + pz2)*(pz1 + pz2));
        
    return sqrt(invMassSq);
}

void CCProtonPi0::initVariables()
{
    cout<<"Initializing Interaction"<<endl;
    
    rootDir = Folder_List::output + Folder_List::rootOut + branchDir + "Interaction.root";
    plotDir = Folder_List::output + Folder_List::plotOut + branchDir + "Interaction/";
    
    cout<<"\tRoot File: "<<rootDir<<endl;
    cout<<"\tPlot Output Folder: "<<plotDir<<endl;
 
    // Create Root File 
    f = new TFile(rootDir.c_str(),"RECREATE");
    
    initHistograms();
    
    // Open Text Files
    openTextFiles();
    
    max_nFSPart = 10;
    nBckgBranch = 6;
    
    // -------------------------------------------------------------------------
    //     Initialization
    //--------------------------------------------------------------------------
    initCutNumbers(); 
    
    hasParticleTruthInfo = false;
    
    // Default Beam Configuration
    beam_p3.SetXYZ(1.0,1.0,1.0);
    beam_p3.SetPhi(-1.554);
    beam_p3.SetTheta(0.059);
    
    // Background Types
    for(int i = 0; i < nBckgBranch; i++){
            bckg_1Prong_QELike.push_back(0.0);
            bckg_1Prong_SinglePiPlus.push_back(0.0);
            bckg_1Prong_SinglePiMinus.push_back(0.0);
            bckg_1Prong_MultiPion.push_back(0.0);
            bckg_1Prong_MultiPiZero.push_back(0.0);
            bckg_1Prong_Other.push_back(0.0);
            bckg_2Prong_QELike.push_back(0.0);
            bckg_2Prong_SinglePiPlus.push_back(0.0);
            bckg_2Prong_SinglePiMinus.push_back(0.0);
            bckg_2Prong_MultiPion.push_back(0.0);
            bckg_2Prong_MultiPiZero.push_back(0.0);
            bckg_2Prong_Other.push_back(0.0);
    }
    
       

    
    
    cout<<"Done!"<<endl;
    
    muon.initialize(anaMode);
    proton.initialize(anaMode);
    pion.initialize(anaMode);
}

void CCProtonPi0::initCutNumbers()
{
    // Common Cut Numbers
    nCut_All.set_Name("All");
    nCut_Vertex_None.set_Name("Vertex_None");
    nCut_Vertex_Null.set_Name("Vertex_Null");
    nCut_Vertex_Not_Reconstructable.set_Name("Vertex_Not_Reconstructable"); 
    nCut_Vertex_Not_Fiducial.set_Name("Vertex_Not_Fiducial");
    nCut_Vertex_Count.set_Name("Vertex_Count");  
    nCut_nProngs.set_Name("nProngs");
    
    // nProngs == 1 Cut Numbers
    nCut_1Prong_Muon_None.set_Name("Muon_None");              
    nCut_1Prong_Muon_Not_Plausible.set_Name("Muon_Not_Plausible");
    nCut_1Prong_Muon_Score_Low.set_Name("Muon_Score_Low");
    nCut_1Prong_Muon_Charge.set_Name("Muon_Charge");
    nCut_1Prong_Vertex_Michel_Exist.set_Name("Vertex_Michel_Exist"); 
    nCut_1Prong_EndPoint_Michel_Exist.set_Name("EndPoint_Michel_Exist");
    nCut_1Prong_secEndPoint_Michel_Exist.set_Name("secEndPoint_Michel_Exist");
    nCut_1Prong_PreFilter_Pi0.set_Name("PreFilter_Pi0");
    nCut_1Prong_VtxBlob.set_Name("VtxBlob");
    nCut_1Prong_ConeBlobs.set_Name("ConeBlobs");
    nCut_1Prong_PhotonDistanceLow.set_Name("PhotonDistanceLow");
    nCut_1Prong_Pi0_invMass.set_Name("Pi0_invMass");
    nCut_1Prong_beamEnergy.set_Name("beamEnergy");
    nCut_1Prong_UnusedE.set_Name("UnusedE");
    
    // nProngs == 1 Cut Numbers
    nCut_2Prong_Muon_None.set_Name("Muon_None");              
    nCut_2Prong_Muon_Not_Plausible.set_Name("Muon_Not_Plausible");
    nCut_2Prong_Muon_Score_Low.set_Name("Muon_Score_Low");
    nCut_2Prong_Muon_Charge.set_Name("Muon_Charge");
    nCut_2Prong_Vertex_Michel_Exist.set_Name("Vertex_Michel_Exist"); 
    nCut_2Prong_EndPoint_Michel_Exist.set_Name("EndPoint_Michel_Exist");
    nCut_2Prong_secEndPoint_Michel_Exist.set_Name("secEndPoint_Michel_Exist");
    nCut_2Prong_Particle_None.set_Name("Particle_None");
    nCut_2Prong_Proton_None.set_Name("Proton_None");            
    nCut_2Prong_Proton_Score.set_Name("Proton_Score");
    nCut_2Prong_Pion_Score.set_Name("Pion_Score");
    nCut_2Prong_pIDDiff.set_Name("nCut_pIDDiff");
    nCut_2Prong_Proton_Score_LLR.set_Name("Proton_Score_LLR");
    nCut_2Prong_PreFilter_Pi0.set_Name("PreFilter_Pi0");
    nCut_2Prong_VtxBlob.set_Name("VtxBlob");
    nCut_2Prong_ConeBlobs.set_Name("ConeBlobs");
    nCut_2Prong_PhotonDistanceLow.set_Name("PhotonDistanceLow");
    nCut_2Prong_Pi0_invMass.set_Name("Pi0_invMass");
    nCut_2Prong_beamEnergy.set_Name("beamEnergy");
    nCut_2Prong_UnusedE.set_Name("UnusedE");
    
}


void CCProtonPi0::fillHistograms()
{
    muon.fill_Histograms();
    pion.fill_Histograms();
    proton.fill_Histograms();
}

void CCProtonPi0::write_RootFile()
{
    cout<<">> Writing "<<rootDir<<endl;
    f->Write();
}

void CCProtonPi0::closeTextFiles()
{
    readme.close(); 
    
    for (int i = 0; i < nTopologies; i++){
        cutText[i].close();
        backgroundText[i].close();
        failText[i].close();
    }
}

void CCProtonPi0::openTextFiles()
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
    
    // Open Cut Files
    cutFile[0] = Folder_List::output + Folder_List::textOut + branchDir + "CutTable_1Prong.txt";
    cutFile[1] = Folder_List::output + Folder_List::textOut + branchDir + "CutTable_2Prongs.txt";
    
    for (int i = 0; i < nTopologies; i++){
        cutText[i].open( cutFile[i].c_str() );
        if( !cutText[i].is_open() ){
            cerr<<"Cannot open output text file: "<<cutFile[i]<<endl;
            exit(1);
        }else{
            cout<<"\t"<<cutFile[i]<<endl;
        }
    } 
    
    // Open Background Files
    backgroundFile[0] = Folder_List::output + Folder_List::textOut + branchDir + "BackgroundTable_1Prong.txt";
    backgroundFile[1] = Folder_List::output + Folder_List::textOut + branchDir + "BackgroundTable_2Prongs.txt";
    
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
    failFile[1] = Folder_List::output + Folder_List::textOut + branchDir + "FailChecks_2Prongs.txt";
    
    for (int i = 0; i < nTopologies; i++){
        failText[i].open( failFile[i].c_str() );
        if( !failText[i].is_open() ){
            cerr<<"Cannot open output text file: "<<failFile[i]<<endl;
            exit(1);
        }else{
            cout<<"\t"<<failFile[i]<<endl;
        }
    }
    
    // Open Arachne RoundupFile
    roundupFile = Folder_List::output + Folder_List::textOut + branchDir + "ArachneRoundup.txt";
    roundupText.open(roundupFile.c_str() );
    if( !roundupText.is_open() ){
        cerr<<"Cannot open Output File!"<<endl;
        exit(1);
    }

    
}

void CCProtonPi0::writeReadme()
{
    readme<<"Test"<<endl;
}

void CCProtonPi0::writeScanList(Long64_t entryNo)
{

string arachne_html = "http://minerva05.fnal.gov/Arachne/arachne.html?filename=";
string entryString  = "&entry=";;
string other        = "&slice=-1&filetype=dst";
string filename     = "/minerva/data/users/oaltinok/CCProtonPi0/MC/v1_05/test/grid/central_value/minerva/dst/v10r6p13/00/01/32/00/SIM_minerva_00013200_Subruns_0034-0035-0036-0037-0038_CCProtonPi0_Ana_DST_v10r6p13.root";

roundupText<<arachne_html<<filename<<entryString<<entryNo<<other<<" ";
roundupText<<ev_run<<"|"<<ev_subrun<<"|"<<ev_gate<<endl;
   


// http://minerva05.fnal.gov/Arachne/arachne.html?filename=/minerva/data/users/oaltinok/CCProtonPi0/MC/v1_05/test/grid/central_value/minerva/dst/v10r6p13/00/01/32/01/SIM_minerva_00013201_Subruns_0001-0002-0003-0004-0005_CCProtonPi0_Ana_DST_v10r6p13.root&entry=2262&slice=-1&filetype=dst	SIM_minerva|3597|57|208|All Slices

// http://minerva05.fnal.gov/Arachne/arachne.html?filename=/minerva/data/users/oaltinok/CCProtonPi0/MC/v1_05/test/grid/central_value/minerva/ana/v10r6p13/00/01/32/01/SIM_minerva_00013201_Subruns_0001-0002-0003-0004-0005_CCProtonPi0_Ana_Tuple_v10r6p13.root&entry=550&slice=-1&filetype=dst 3687|31|201
}


/*
--------------------------------------------------------------------------------
countParticles:
Returns the number of particles in the Final State
Input 
int targetPDG
bool applyPCut - Variable for selecting particles with momentum 
(no particle at rest)
--------------------------------------------------------------------------------
*/
int CCProtonPi0::countParticles(int targetPDG, bool applyPCut)
{
    int count = 0;
    TVector3 p3;
    
    for(int i = 0; i < mc_nFSPart && i < max_nFSPart; i++ ){
        if( mc_FSPartPDG[i] == targetPDG){
            if(applyPCut){
                p3.SetXYZ(mc_FSPartPx[i],mc_FSPartPy[i],mc_FSPartPz[i]);
                if(p3.Mag() > 0){
                    count++;
                }
            }
            else{
                count++;
            }
        }
    }
    
    return count;
    
}

void CCProtonPi0::setChannelTag(string input)
{
    channelTag = input;
    cout<<"Channel Tag = "<<input<<endl;
}


#endif //CCProtonPi0_cpp
