/*
    See CCProtonPi0.h header for Class Information
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
    
    applyProtonScore = true;
    minProtonScore = 0.3;
    
    applyPionScore = true;
    maxPionScore = 0.3;
    
    applyPhotonDistance = false;
    minPhotonDistance = 150; //mm
    
    applyBeamEnergy = true;
    max_beamEnergy = 8.0; // GeV
    
    applyQSq = true;
    max_QSq = 1.5;
    
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
    }
    
    if(m_ActivatePi0){
        fChain->SetBranchStatus("pi0_*",1);
        fChain->SetBranchStatus("gamma*",1);
    }
    
        
    // Fail Checks
        double nAntiMuon = 0;
        double nMuonChargeDiff = 0;
        // True Number of Pi0
        double n0Pi0 = 0;
        double n1Pi0 = 0;
        double nMultPi0 = 0;
        
        int nPi0Count;
        double n0Pi0_Test = 0;
        double n1Pi0_Test = 0;
        double nMultPi0_Test = 0;
        
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
        
        if (applyMaxEvents && jentry == nMaxEvents){
            cout<<"\tReached Event Limit!"<<endl;
            break;
        }
        
        // Decide to Analyze Event or NOT
        if( !analyzeEvent() ) continue;
        
        //----------------------------------------------------------------------
        // Variable Test Site
        //----------------------------------------------------------------------
//         if (truth_isSignal){
//             cout<<"------------------------------"<<endl;
//             cout<<"Pi0 Status = "<<truth_pi0_status<<endl;
//             cout<<"Pi0 Mother = "<<truth_pi0_Mother<<endl;
//             cout<<"Pi0 Mother Status = "<<truth_pi0_MotherStatus<<endl;
//             cout<<"Pi0 GrandMother = "<<truth_pi0_GrandMother<<endl;
//             cout<<"Pi0 GrandMother Status = "<<truth_pi0_GrandMotherStatus<<endl;
//         }
       
       //----------------------------------------------------------------------
       // Get Cut Statistics
       //----------------------------------------------------------------------
        isPassedAllCuts = getCutStatistics();
        if( !isPassedAllCuts ) continue;
        
        
        //----------------------------------------------------------------------
        // Data Analysis and Other Studies
        //----------------------------------------------------------------------
        if (isDataAnalysis) fillData();
        if (is_pID_Studies) fill_pID_Data();
        if (writeFSParticleMomentum) writeFSParticle4P(jentry);
        
        //----------------------------------------------------------------------
        // Fail Checks
        //----------------------------------------------------------------------
        // Anti-Muon Contamination
        if((truth_muon_charge - CCProtonPi0_muon_charge) != 0){
            nMuonChargeDiff++;
        }
        if(truth_muon_charge == 1) nAntiMuon++;
        
        // True Pi0 Count
        if(truth_N_pi0 == 0) n0Pi0++;
        if(truth_N_pi0 == 1) n1Pi0++;
        if(truth_N_pi0 > 1) nMultPi0++;
        
        nPi0Count = countParticles(111,false);
        
        if(nPi0Count == 0) n0Pi0_Test++;
        if(nPi0Count == 1) n1Pi0_Test++;
        if(nPi0Count > 1) nMultPi0_Test++;

    } // end for-loop
    
    getPi0Family();
     
    failText<<"N(Muon Charge Diff) = "<<nMuonChargeDiff<<endl;
    failText<<"N(AntiMuon) = "<<nAntiMuon<<endl;
    failText<<"N(n0Pi0) = "<<n0Pi0<<" | N(n0Pi0_Test) = "<<n0Pi0_Test<<endl;
    failText<<"N(n1Pi0) = "<<n1Pi0<<" | N(n1Pi0_Test) = "<<n1Pi0_Test<<endl;
    failText<<"N(nMultPi0) = "<<nMultPi0<<" | N(nMultPi0_Test) = "<<nMultPi0_Test<<endl;

    //--------------------------------------------------------------------------
    // Write Files
    //--------------------------------------------------------------------------
    writeCutTable();
    
    // Write the Root Files
    write_RootFile();           //CCProtonPi0
    muon.write_RootFile();
    proton.write_RootFile();
    pion.write_RootFile();
    
    
    closeFiles();
    
//     if(is_pID_Studies) get_pID_Stats();
    
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
//
//     Specific Functions
//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

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
            failText<<"Mother = "<<PDG_pi0_Mother[i]<<endl;
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
            failText<<"GrandMother = "<<PDG_pi0_GrandMother[i]<<endl;
        }
    }
    
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

void CCProtonPi0::fill_pID_Data()
{
    if(CCProtonPi0_trajProtonProngPDG[indRecoProton] == 2212){
        pID_proton_protonScore->Fill(CCProtonPi0_protonScore[indRecoProton]);
        pID_proton_pionScore->Fill(CCProtonPi0_pionScore[indRecoProton]);
        pID_proton_pionScore_protonScore->Fill(CCProtonPi0_pionScore[indRecoProton],CCProtonPi0_protonScore[indRecoProton]);
    }else if(CCProtonPi0_trajProtonProngPDG[indRecoProton] == 211){
        pID_piplus_protonScore->Fill(CCProtonPi0_protonScore[indRecoProton]);
        pID_piplus_pionScore->Fill(CCProtonPi0_pionScore[indRecoProton]);
        pID_piplus_pionScore_protonScore->Fill(CCProtonPi0_pionScore[indRecoProton],CCProtonPi0_protonScore[indRecoProton]);
    }else if(CCProtonPi0_trajProtonProngPDG[indRecoProton] == -211){
        pID_piminus_protonScore->Fill(CCProtonPi0_protonScore[indRecoProton]);
        pID_piminus_pionScore->Fill(CCProtonPi0_pionScore[indRecoProton]);
        pID_piminus_pionScore_protonScore->Fill(CCProtonPi0_pionScore[indRecoProton],CCProtonPi0_protonScore[indRecoProton]);
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
            <<CCProtonPi0_proton_px[indRecoProton]<<", "
            <<CCProtonPi0_proton_py[indRecoProton]<<", "
            <<CCProtonPi0_proton_pz[indRecoProton]<<", "
            <<CCProtonPi0_proton_E[indRecoProton]<<" )"
            <<" Score = "<<CCProtonPi0_proton_score[indRecoProton]
            <<endl;
    failText<<"Pi0 4-P = ( "
            <<pi0_px<<", "
            <<pi0_py<<", "
            <<pi0_pz<<", "
            <<pi0_E<<" )"
            <<endl;   
}

bool CCProtonPi0::getCutStatistics()
{
    bool isDIS = (mc_intType == 3);
    //----------------------------------------------------------------------
    // Count Events after each Reconstruction Cut
    //----------------------------------------------------------------------
    // Count All Events before Cuts
    nCut_All.increment(truth_isSignal, truth_isSignalGold, truth_isSignalSilver1, isDIS );
  
    if ( Cut_Event_Has_BadObject == 1) return false;
    nCut_Event_Has_BadObject.increment(truth_isSignal, truth_isSignalGold, truth_isSignalSilver1, isDIS );
  
    if( Cut_Vertex_None == 1) return false;
    nCut_Vertex_None.increment(truth_isSignal, truth_isSignalGold, truth_isSignalSilver1, isDIS );
    
    if( Cut_Vertex_Null == 1) return false;
    nCut_Vertex_Null.increment(truth_isSignal, truth_isSignalGold, truth_isSignalSilver1, isDIS );
    
    if( Cut_Vertex_Not_Reconstructable == 1) return false;
    nCut_Vertex_Not_Reconstructable.increment(truth_isSignal, truth_isSignalGold, truth_isSignalSilver1, isDIS );
    
    if( Cut_Vertex_Not_Fiducial == 1) return false;
    nCut_Vertex_Not_Fiducial.increment(truth_isSignal, truth_isSignalGold, truth_isSignalSilver1, isDIS );
        
    if( Cut_Muon_None == 1) return false;
    nCut_Muon_None.increment(truth_isSignal, truth_isSignalGold, truth_isSignalSilver1, isDIS );
         
    if( Cut_Muon_Not_Plausible == 1) return false;
    nCut_Muon_Not_Plausible.increment(truth_isSignal, truth_isSignalGold, truth_isSignalSilver1, isDIS );
    
    //--------------------------------------------------------------------------
    //  Fill Cut Parameters
    //--------------------------------------------------------------------------
        // Fill nProngs
        nProngs_hist_initial->Fill(nProngs);

        // Fill mc_intType
        if(truth_isSignal){
            if(mc_intType == 1) mc_w_CCQE->Fill(mc_w * HEP_Functions::MeV_to_GeV);
            if(mc_intType == 2) mc_w_RES->Fill(mc_w * HEP_Functions::MeV_to_GeV);
            if(mc_intType == 3) mc_w_DIS->Fill(mc_w * HEP_Functions::MeV_to_GeV);
            status_Pi0->Fill(truth_pi0_status);
            status_Pi0_Mother->Fill(truth_pi0_MotherStatus);
            status_Pi0_GrandMother->Fill(truth_pi0_GrandMotherStatus);
            PDG_pi0_Mother.push_back(truth_pi0_Mother);
            PDG_pi0_GrandMother.push_back(truth_pi0_GrandMother);
        }
    //--------------------------------------------------------------------------
    
    if( Cut_Muon_Score_Low == 1) return false;
    nCut_Muon_Score_Low.increment(truth_isSignal, truth_isSignalGold, truth_isSignalSilver1, isDIS );
    
    if( Cut_Muon_Charge == 1) return false;
    nCut_Muon_Charge.increment(truth_isSignal, truth_isSignalGold, truth_isSignalSilver1, isDIS );
    
    if( Cut_nProngs == 1 || nProngs > 3) return false;
    nCut_nProngs.increment(truth_isSignal, truth_isSignalGold, truth_isSignalSilver1, isDIS );
    
    if( Cut_Vertex_Michel_Exist == 1) return false;
    nCut_Vertex_Michel_Exist.increment(truth_isSignal, truth_isSignalGold, truth_isSignalSilver1, isDIS );
    
    if( Cut_EndPoint_Michel_Exist == 1) return false;
    nCut_EndPoint_Michel_Exist.increment(truth_isSignal, truth_isSignalGold, truth_isSignalSilver1, isDIS );
  
    if( Cut_secEndPoint_Michel_Exist == 1) return false;
    nCut_secEndPoint_Michel_Exist.increment(truth_isSignal, truth_isSignalGold, truth_isSignalSilver1, isDIS );
    
    if( Cut_Particle_None == 1 && !analyze_NoProtonEvents) return false;
    nCut_Particle_None.increment(truth_isSignal, truth_isSignalGold, truth_isSignalSilver1, isDIS );
    
    if( Cut_Proton_None == 1 && !analyze_NoProtonEvents) return false;
    nCut_Proton_None.increment(truth_isSignal, truth_isSignalGold, truth_isSignalSilver1, isDIS );
    
    if ( Cut_Particle_None == 1 || Cut_Proton_None == 1 ){
        nCut_Proton_Score.increment(truth_isSignal, truth_isSignalGold, truth_isSignalSilver1, isDIS );
        nCut_Pion_Score.increment(truth_isSignal, truth_isSignalGold, truth_isSignalSilver1, isDIS );
        //         cout<<"WARNING: No Proton!"<<endl;
    }else{
        // Find Best Proton in Reco
        findRecoProton(); 
        
        if ( applyProtonScore && (CCProtonPi0_protonScore[indRecoProton] < minProtonScore) ) return false;
        nCut_Proton_Score.increment(truth_isSignal, truth_isSignalGold, truth_isSignalSilver1, isDIS );
        
        if ( applyPionScore && (CCProtonPi0_pionScore[indRecoProton] > maxPionScore) ) return false;
        nCut_Pion_Score.increment(truth_isSignal, truth_isSignalGold, truth_isSignalSilver1, isDIS );
        
        proton.protonScore->Fill(CCProtonPi0_protonScore[indRecoProton]);
        proton.pionScore->Fill(CCProtonPi0_pionScore[indRecoProton]);
        
    }
    

    //--------------------------------------------------------------------------
    // fill Pre Filter   
        pFilter_Status->Fill(preFilter_Result);
        pFilter_RejectedEnergy->Fill(preFilter_rejectedEnergy);
    //--------------------------------------------------------------------------
    
    if( Cut_PreFilter_Pi0 == 1) return false;
    nCut_PreFilter_Pi0.increment(truth_isSignal, truth_isSignalGold, truth_isSignalSilver1, isDIS );
    
    if( Cut_VtxBlob == 1) return false;
    nCut_VtxBlob.increment(truth_isSignal, truth_isSignalGold, truth_isSignalSilver1, isDIS );
    
    if( Cut_ConeBlobs == 1) return false;
    nCut_ConeBlobs.increment(truth_isSignal, truth_isSignalGold, truth_isSignalSilver1, isDIS );
    
    if( pi0_E == SENTINEL) return false;
    nCut_Other.increment(truth_isSignal, truth_isSignalGold, truth_isSignalSilver1, isDIS );
    
    if( pi0_invMass < min_Pi0_invMass || pi0_invMass > max_Pi0_invMass ) return false;
    nCut_Pi0_invMass.increment(truth_isSignal, truth_isSignalGold, truth_isSignalSilver1, isDIS );
    
    if( applyPhotonDistance && isPhotonDistanceLow()) return false;
    nCut_PhotonDistanceLow.increment(truth_isSignal, truth_isSignalGold, truth_isSignalSilver1, isDIS );
    
    if( applyBeamEnergy && (CCProtonPi0_neutrino_E * HEP_Functions::MeV_to_GeV > max_beamEnergy)) return false;
    nCut_beamEnergy.increment(truth_isSignal, truth_isSignalGold, truth_isSignalSilver1, isDIS );
    
    if( applyQSq && (CCProtonPi0_QSq * HEP_Functions::MeVSq_to_GeVSq > max_QSq)) return false;
    nCut_QSq.increment(truth_isSignal, truth_isSignalGold, truth_isSignalSilver1, isDIS );
    
    
    
    return true;
    
}

void CCProtonPi0::writeCutTableHeader()
{
    cutText<<std::left;
    
    cutText.width(35); cutText<<"Cut"<<" | "; 
    
    cutText.width(12); cutText<<"N(Events)"<<" | ";    
    cutText.width(12); cutText<<"N(Signal)"<<" | ";      
    cutText.width(12); cutText<<"Efficiency"<<" | ";      
    cutText.width(12); cutText<<"Purity"<<" | ";
    
    cutText.width(12); cutText<<"N(Gold)"<<" | ";    
    cutText.width(12); cutText<<"Efficiency"<<" | ";   
    cutText.width(12); cutText<<"Purity"<<" | ";
    
    cutText.width(12); cutText<<"N(Silver1)"<<" | ";   
    cutText.width(12); cutText<<"Efficiency"<<" | ";     
    cutText.width(12); cutText<<"Purity"<<" | "; 
    
    cutText.width(12); cutText<<"N(DIS)"<<" | ";   
    cutText.width(12); cutText<<"Efficiency"<<" | ";     
    cutText.width(12); cutText<<"Purity"<<" | "<<endl;
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


void CCProtonPi0::writeCutTable()
{
    double efficiency_Total;
    double purity_Total;
    double efficiency_Gold;
    double purity_Gold;
    double efficiency_Silver1;
    double purity_Silver1;
    
    double efficiencyBase_Total;
    double efficiencyBase_Gold;
    double efficiencyBase_Silver1;
    
    efficiencyBase_Total = nCut_Muon_None.nSignal.getCount();
    efficiencyBase_Gold = nCut_Muon_None.nSignal_Gold.getCount();
    efficiencyBase_Silver1 = nCut_Muon_None.nSignal_Silver1.getCount();
    
    // Form nCutVector
    nCutVector.push_back(nCut_All);
    nCutVector.push_back(nCut_Event_Has_BadObject);
    nCutVector.push_back(nCut_Vertex_None);
    nCutVector.push_back(nCut_Vertex_Null);
    nCutVector.push_back(nCut_Vertex_Not_Reconstructable); 
    nCutVector.push_back(nCut_Vertex_Not_Fiducial);    
    nCutVector.push_back(nCut_Muon_None);              
    nCutVector.push_back(nCut_Muon_Not_Plausible);
    nCutVector.push_back(nCut_Muon_Score_Low);
    nCutVector.push_back(nCut_Muon_Charge);
    nCutVector.push_back(nCut_nProngs);
    nCutVector.push_back(nCut_Vertex_Michel_Exist); 
    nCutVector.push_back(nCut_EndPoint_Michel_Exist);
    nCutVector.push_back(nCut_secEndPoint_Michel_Exist);
    nCutVector.push_back(nCut_Particle_None);
    nCutVector.push_back(nCut_Proton_None);
    nCutVector.push_back(nCut_Proton_Score);
    nCutVector.push_back(nCut_Pion_Score);
    nCutVector.push_back(nCut_PreFilter_Pi0);
    nCutVector.push_back(nCut_VtxBlob);
    nCutVector.push_back(nCut_ConeBlobs);
    nCutVector.push_back(nCut_Other);
    nCutVector.push_back(nCut_Pi0_invMass);
    nCutVector.push_back(nCut_PhotonDistanceLow);
    nCutVector.push_back(nCut_beamEnergy);
    nCutVector.push_back(nCut_QSq);

    
    cout<<">> Writing "<<cutFile<<endl;
    
    writeCutTableHeader();
 
    for( unsigned int i = 0; i < nCutVector.size(); i++){
        
        efficiency_Total = getCutEfficiency(nCutVector[i].nSignal.getCount(),efficiencyBase_Total);
        efficiency_Gold = getCutEfficiency(nCutVector[i].nSignal_Gold.getCount(),efficiencyBase_Gold);
        efficiency_Silver1 = getCutEfficiency(nCutVector[i].nSignal_Silver1.getCount(),efficiencyBase_Silver1);
        
        purity_Total = getCutPurity(nCutVector[i].nSignal.getCount(),nCutVector[i].nEvent.getCount());
        purity_Gold = getCutPurity(nCutVector[i].nSignal_Gold.getCount(),nCutVector[i].nEvent.getCount());
        purity_Silver1 = getCutPurity(nCutVector[i].nSignal_Silver1.getCount(),nCutVector[i].nEvent.getCount());
        
        cutText.unsetf( std::ios::floatfield ); 
        cutText.width(35); cutText<<nCutVector[i].get_Name()<<" | ";
        cutText.width(12); cutText<<nCutVector[i].nEvent.getCount()<<" | ";
        
        // Total Signal
        cutText.width(12); cutText<<nCutVector[i].nSignal.getCount()<<" | ";

        cutText.precision(4); 
       
        if ( efficiency_Total <= 100){
            cutText.width(12); cutText<<efficiency_Total<<" | ";
        }else{
            cutText.width(12); cutText<<"N/A"<<" | ";    
        }

        if ( efficiency_Total <= 100){
            cutText.width(12); cutText<<purity_Total<<" | ";
        }else{
            cutText.width(12);  cutText<<"N/A"<<" | ";    
        }
        
        // Signal: Gold
        cutText.unsetf( std::ios::floatfield );  // Important
        cutText.width(12); cutText<<nCutVector[i].nSignal_Gold.getCount()<<" | ";
        
        cutText.precision(4); 
        
        if ( efficiency_Gold <= 100){
            cutText.width(12); cutText<<efficiency_Gold<<" | ";
        }else{
            cutText.width(12); cutText<<"N/A"<<" | ";    
        }
        
        if ( efficiency_Gold <= 100){
            cutText.width(12); cutText<<purity_Gold<<" | ";
        }else{
            cutText.width(12);  cutText<<"N/A"<<" | ";    
        }
        
        // Signal: Silver1
        cutText.unsetf( std::ios::floatfield );  // Important
        cutText.width(12); cutText<<nCutVector[i].nSignal_Silver1.getCount()<<" | ";
        
        cutText.precision(4); 
        
        if ( efficiency_Silver1 <= 100){
            cutText.width(12); cutText<<efficiency_Silver1<<" | ";
        }else{
            cutText.width(12); cutText<<"N/A"<<" | ";    
        }
        
        if ( efficiency_Silver1 <= 100){
            cutText.width(12); cutText<<purity_Silver1<<" | ";
        }else{
            cutText.width(12);  cutText<<"N/A"<<" | ";    
        }
        
        // DIS
        cutText.unsetf( std::ios::floatfield );  // Important
        cutText.width(12); cutText<<nCutVector[i].nDIS.getCount()<<" | ";
        
        cutText<<endl;
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

    double WSq;
    if( CCProtonPi0_WSq < 0 ){ 
        WSq = 0;
    }else{
        WSq = CCProtonPi0_WSq;
    }
    
    w2_mc->Fill(mc_w * HEP_Functions::MeV_to_GeV);
    w2_reco->Fill(sqrt(WSq)* HEP_Functions::MeV_to_GeV);
    w2_error->Fill( Data_Functions::getError(mc_w,sqrt(WSq)) );
    w2_reco_mc->Fill(sqrt(WSq)* HEP_Functions::MeV_to_GeV,mc_w * HEP_Functions::MeV_to_GeV);
    
    vertex_x_y_reco->Fill(CCProtonPi0_vtx[0],CCProtonPi0_vtx[1]);
    vertex_z_reco->Fill(CCProtonPi0_vtx[2]);
    
    nProngs_hist->Fill(nProngs);
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
    openFiles();
    
    max_nFSPart = 10;
    
    // -------------------------------------------------------------------------
    //     Initialization
    //--------------------------------------------------------------------------
    hasParticleTruthInfo = false;
    
    // Default Beam Configuration
    beam_p3.SetXYZ(1.0,1.0,1.0);
    beam_p3.SetPhi(-1.554);
    beam_p3.SetTheta(0.059);
    
    // Cut Numbers
    nCut_All.set_Name("All");
    nCut_Event_Has_BadObject.set_Name("Event_Has_BadObjec");
    nCut_Vertex_None.set_Name("Vertex_None");
    nCut_Vertex_Null.set_Name("Vertex_Null");
    nCut_Vertex_Not_Reconstructable.set_Name("Vertex_Not_Reconstructable"); 
    nCut_Vertex_Not_Fiducial.set_Name("Vertex_Not_Fiducial");    
    nCut_Muon_None.set_Name("Muon_None");              
    nCut_Muon_Not_Plausible.set_Name("Muon_Not_Plausible");
    nCut_Muon_Score_Low.set_Name("Muon_Score_Low");
    nCut_Muon_Charge.set_Name("Muon_Charge");
    nCut_nProngs.set_Name("nProngs");
    nCut_Vertex_Michel_Exist.set_Name("Vertex_Michel_Exist"); 
    nCut_EndPoint_Michel_Exist.set_Name("EndPoint_Michel_Exist");
    nCut_secEndPoint_Michel_Exist.set_Name("secEndPoint_Michel_Exist");
    nCut_Particle_None.set_Name("Particle_None");
    nCut_Proton_None.set_Name("Proton_None");            
    nCut_Proton_Score.set_Name("Proton_Score");
    nCut_Pion_Score.set_Name("Pion_Score");
    nCut_PreFilter_Pi0.set_Name("PreFilter_Pi0");
    nCut_VtxBlob.set_Name("VtxBlob");
    nCut_ConeBlobs.set_Name("ConeBlobs");
    nCut_PhotonDistanceLow.set_Name("PhotonDistanceLow");
    nCut_Other.set_Name("Other");
    nCut_Pi0_invMass.set_Name("Pi0_invMass");
    nCut_beamEnergy.set_Name("beamEnergy");
    nCut_QSq.set_Name("QSq");
    
    cout<<"Done!"<<endl;
    
    muon.initialize(anaMode);
    proton.initialize(anaMode);
    pion.initialize(anaMode);
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

void CCProtonPi0::closeFiles()
{
    readme.close();
    cutText.close();
    failText.close();
}

void CCProtonPi0::openFiles()
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
    
    // Open Cut File
    cutFile = Folder_List::output + Folder_List::textOut + branchDir + "CutTable";
    cutFile = cutFile + "_" + channelTag + ".txt";
    cutText.open( cutFile.c_str() );
    if( !cutText.is_open() ){
        cerr<<"Cannot open output text file: "<<cutFile<<endl;
        exit(1);
    }else{
        cout<<"\t"<<cutFile<<endl;
    }
    
    // Open Fail File
    failFile = Folder_List::output + Folder_List::textOut + branchDir + "FailChecks.txt";
    failText.open( failFile.c_str() );
    if( !failText.is_open() ){
        cerr<<"Cannot open output text file: "<<failFile<<endl;
        exit(1);
    }else{
        cout<<"\t"<<failFile<<endl;
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
