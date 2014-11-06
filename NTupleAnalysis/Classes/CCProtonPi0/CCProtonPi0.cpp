/*
    See CCProtonPi0.h header for Class Information
*/

#ifndef CCProtonPi0_cpp
#define CCProtonPi0_cpp

#include "CCProtonPi0.h"

using namespace std;

void CCProtonPi0::specifyRunTime()
{
    channelTag = "test2";
    
    // Control Flow
    isMC            = true;
    isDataAnalysis  = true;
    is_pID_Studies  = false;
    writeFSParticleMomentum = false;
    
    applyProtonScore = false;
    minProtonScore = 0.3;
    
    applyPhotonDistance = false;
    minPhotonDistance = 150; //mm
    
    SENTINEL = -9.9;

    //Select Branches to Activate
    m_ActivateMC            = true;
    m_ActivateInteraction   = true;
    m_ActivateMuon          = true;
    m_ActivateProton        = true;
    m_ActivatePi0           = true;
   
}

void CCProtonPi0::run(string playlist)
{
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
    initHistograms();
    
    //------------------------------------------------------------------------
    // Open files for writing
    //------------------------------------------------------------------------
    openFiles();
    
    //------------------------------------------------------------------------
    // Branch Selection for Performance
    //------------------------------------------------------------------------
    fChain->SetBranchStatus("*",0);  // disable all branches
    
    fChain->SetBranchStatus("Cut_*",1);  // Cut List Activated by Default
    fChain->SetBranchStatus("ev_*",1);
    
    if(m_ActivateMC){
        fChain->SetBranchStatus("truth_*",1);
        fChain->SetBranchStatus("mc_*",1);
        fChain->SetBranchStatus("CCProtonPi0_traj*",1);
    }
    
    if(m_ActivateInteraction){
        fChain->SetBranchStatus("CCProtonPi0_vtx*",1);
        fChain->SetBranchStatus("preFilter_*",1);
        fChain->SetBranchStatus("nProngs",1);
    }
    
    if(m_ActivateMuon){
        fChain->SetBranchStatus("CCProtonPi0_muon_*",1);
    }
    
    if(m_ActivateProton){
        fChain->SetBranchStatus("CCProtonPi0_proton_*",1);
    }
    
    if(m_ActivatePi0){
        fChain->SetBranchStatus("pi0_*",1);
        fChain->SetBranchStatus("gamma*",1);
    }
    
    
//     int indTruePion;
    
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
        // Variable Test Site
        //----------------------------------------------------------------------
//         cout<<"nProngs = "<<nProngs<<endl;
       
       //----------------------------------------------------------------------
       // Get Cut Statistics
       //----------------------------------------------------------------------
        isPassedAllCuts = getCutStatistics();
        
        if( !isPassedAllCuts ) continue;
        if(writeFSParticleMomentum) writeFSParticle4P(jentry);
        
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
        
        
        //----------------------------------------------------------------------
        // pID Studies
        //----------------------------------------------------------------------
        if( is_pID_Studies){
            if(CCProtonPi0_trajProtonProngPDG[indRecoProton] == 2212){
                pID_proton->Fill(CCProtonPi0_proton_score[indRecoProton]);
            }else if(CCProtonPi0_trajProtonProngPDG[indRecoProton] == 211){
                pID_piplus->Fill(CCProtonPi0_proton_score[indRecoProton]);
            }else if(CCProtonPi0_trajProtonProngPDG[indRecoProton] == -211){
                pID_piminus->Fill(CCProtonPi0_proton_score[indRecoProton]);
            }else{
                pID_other->Fill(CCProtonPi0_proton_score[indRecoProton]);
            }
        }
        
        if ( isDataAnalysis){
            //------------------------------------------------------------------
            // Fill Particles
            //------------------------------------------------------------------
            if( isMC ){
                fillInteractionTrue(indRecoProton);
                fillMuonTrue();
                fillProtonTrue(indRecoProton);
                fillPionTrue();
            }
            
            // Fill Reconstructed Information
            fillInteractionReco(indRecoProton);
            fillMuonReco();
            fillProtonReco(indRecoProton);
            fillPionReco();
            
            muon.set_errors();
            proton.set_errors();
//             pion.set_errors();

            //------------------------------------------------------------------
            // Fill Histograms
            //------------------------------------------------------------------
            fillHistograms();            
        }

    } // end for-loop
    
    // Form nCutArray and write CutTable
    nCutVector.push_back(nCut_All);
    nCutVector.push_back(nCut_Event_Not_Plausible);
    nCutVector.push_back(nCut_Event_Has_BadObject);
    nCutVector.push_back(nCut_Vertex_None);
    nCutVector.push_back(nCut_Vertex_Null);
    nCutVector.push_back(nCut_Vertex_Not_Reconstructable); 
    nCutVector.push_back(nCut_Vertex_Not_Fiducial);    
    nCutVector.push_back(nCut_nProngs);
    nCutVector.push_back(nCut_Muon_None);              
    nCutVector.push_back(nCut_Muon_Not_Plausible);
    nCutVector.push_back(nCut_Muon_Score_Low);
    nCutVector.push_back(nCut_Muon_Charge);
    nCutVector.push_back(nCut_Vertex_Michel_Exist); 
    nCutVector.push_back(nCut_EndPoint_Michel_Exist);
    nCutVector.push_back(nCut_secEndPoint_Michel_Exist);
    nCutVector.push_back(nCut_Particle_None);
    nCutVector.push_back(nCut_Proton_None);            
    nCutVector.push_back(nCut_Proton_Score);
    nCutVector.push_back(nCut_PreFilter_Pi0);
    nCutVector.push_back(nCut_VtxBlob);
    nCutVector.push_back(nCut_ConeBlobs);
    nCutVector.push_back(nCut_Other);
    nCutVector.push_back(nCut_PhotonDistanceLow);
    
    writeCutTable();
    
    failText<<"N(Muon Charge Diff) = "<<nMuonChargeDiff<<endl;
    failText<<"N(AntiMuon) = "<<nAntiMuon<<endl;
    failText<<"N(n0Pi0) = "<<n0Pi0<<" | N(n0Pi0_Test) = "<<n0Pi0_Test<<endl;
    failText<<"N(n1Pi0) = "<<n1Pi0<<" | N(n1Pi0_Test) = "<<n1Pi0_Test<<endl;
    failText<<"N(nMultPi0) = "<<nMultPi0<<" | N(nMultPi0_Test) = "<<nMultPi0_Test<<endl;

    
    // Write the Root Files
    write_RootFile();           //CCProtonPi0
    muon.write_RootFile();
    proton.write_RootFile();
    pion.write_RootFile();
    
    
    closeFiles();
    
    if(is_pID_Studies) get_pID_Stats();
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
//
//     Specific Functions
//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

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
    //----------------------------------------------------------------------
    // Count Events after each Reconstruction Cut
    //----------------------------------------------------------------------
    // Count All Events before Cuts
    nCut_All.inc_nEvent();
    if (truth_isSignal) nCut_All.inc_nSignal();
    if (truth_isSignalGold) nCut_All.inc_nSignal_Gold();
    if (truth_isSignalSilver1) nCut_All.inc_nSignal_Silver1();
    
    if( Cut_Event_Not_Plausible == 1) return false;
    nCut_Event_Not_Plausible.inc_nEvent();
    if (truth_isSignal) nCut_Event_Not_Plausible.inc_nSignal();
    if (truth_isSignalGold) nCut_Event_Not_Plausible.inc_nSignal_Gold();
    if (truth_isSignalSilver1) nCut_Event_Not_Plausible.inc_nSignal_Silver1();
    
    if ( Cut_Event_Has_BadObject == 1) return false;
    nCut_Event_Has_BadObject.inc_nEvent();
    if (truth_isSignal) nCut_Event_Has_BadObject.inc_nSignal();
    if (truth_isSignalGold) nCut_Event_Has_BadObject.inc_nSignal_Gold();
    if (truth_isSignalSilver1) nCut_Event_Has_BadObject.inc_nSignal_Silver1();
    
    if( Cut_Vertex_None == 1) return false;
    nCut_Vertex_None.inc_nEvent();
    if (truth_isSignal) nCut_Vertex_None.inc_nSignal();
    if (truth_isSignalGold) nCut_Vertex_None.inc_nSignal_Gold();
    if (truth_isSignalSilver1) nCut_Vertex_None.inc_nSignal_Silver1();
    
    if( Cut_Vertex_Null == 1) return false;
    nCut_Vertex_Null.inc_nEvent();
    if (truth_isSignal) nCut_Vertex_Null.inc_nSignal();
    if (truth_isSignalGold) nCut_Vertex_Null.inc_nSignal_Gold();
    if (truth_isSignalSilver1) nCut_Vertex_Null.inc_nSignal_Silver1();
    
    if( Cut_Vertex_Not_Reconstructable == 1) return false;
    nCut_Vertex_Not_Reconstructable.inc_nEvent();
    if (truth_isSignal) nCut_Vertex_Not_Reconstructable.inc_nSignal();
    if (truth_isSignalGold) nCut_Vertex_Not_Reconstructable.inc_nSignal_Gold();
    if (truth_isSignalSilver1) nCut_Vertex_Not_Reconstructable.inc_nSignal_Silver1();
    
    if( Cut_Vertex_Not_Fiducial == 1) return false;
    nCut_Vertex_Not_Fiducial.inc_nEvent();
    if (truth_isSignal) nCut_Vertex_Not_Fiducial.inc_nSignal();
    if (truth_isSignalGold) nCut_Vertex_Not_Fiducial.inc_nSignal_Gold();
    if (truth_isSignalSilver1) nCut_Vertex_Not_Fiducial.inc_nSignal_Silver1();
    
    if( Cut_nProngs == 1) return false;
    nCut_nProngs.inc_nEvent();
    if (truth_isSignal) nCut_nProngs.inc_nSignal();
    if (truth_isSignalGold) nCut_nProngs.inc_nSignal_Gold();
    if (truth_isSignalSilver1) nCut_nProngs.inc_nSignal_Silver1();
    
    if( Cut_Muon_None == 1) return false;
    nCut_Muon_None.inc_nEvent();
    if (truth_isSignal) nCut_Muon_None.inc_nSignal();
    if (truth_isSignalGold) nCut_Muon_None.inc_nSignal_Gold();
    if (truth_isSignalSilver1) nCut_Muon_None.inc_nSignal_Silver1();
    
    if( Cut_Muon_Not_Plausible == 1) return false;
    nCut_Muon_Not_Plausible.inc_nEvent();
    if (truth_isSignal) nCut_Muon_Not_Plausible.inc_nSignal();
    if (truth_isSignalGold) nCut_Muon_Not_Plausible.inc_nSignal_Gold();
    if (truth_isSignalSilver1) nCut_Muon_Not_Plausible.inc_nSignal_Silver1();
    
    if( Cut_Muon_Score_Low == 1) return false;
    nCut_Muon_Score_Low.inc_nEvent();
    if (truth_isSignal) nCut_Muon_Score_Low.inc_nSignal();
    if (truth_isSignalGold) nCut_Muon_Score_Low.inc_nSignal_Gold();
    if (truth_isSignalSilver1) nCut_Muon_Score_Low.inc_nSignal_Silver1();
    
    if( Cut_Muon_Charge == 1) return false;
    nCut_Muon_Charge.inc_nEvent();
    if (truth_isSignal) nCut_Muon_Charge.inc_nSignal();
    if (truth_isSignalGold) nCut_Muon_Charge.inc_nSignal_Gold();
    if (truth_isSignalSilver1) nCut_Muon_Charge.inc_nSignal_Silver1();
    
    if( Cut_Vertex_Michel_Exist == 1) return false;
    nCut_Vertex_Michel_Exist.inc_nEvent();
    if (truth_isSignal) nCut_Vertex_Michel_Exist.inc_nSignal();
    if (truth_isSignalGold) nCut_Vertex_Michel_Exist.inc_nSignal_Gold();
    if (truth_isSignalSilver1) nCut_Vertex_Michel_Exist.inc_nSignal_Silver1();
    
    if( Cut_EndPoint_Michel_Exist == 1) return false;
    nCut_EndPoint_Michel_Exist.inc_nEvent();
    if (truth_isSignal) nCut_EndPoint_Michel_Exist.inc_nSignal();
    if (truth_isSignalGold) nCut_EndPoint_Michel_Exist.inc_nSignal_Gold();
    if (truth_isSignalSilver1) nCut_EndPoint_Michel_Exist.inc_nSignal_Silver1();
    
    if( Cut_secEndPoint_Michel_Exist == 1) return false;
    nCut_secEndPoint_Michel_Exist.inc_nEvent();
    if (truth_isSignal) nCut_secEndPoint_Michel_Exist.inc_nSignal();
    if (truth_isSignalGold) nCut_secEndPoint_Michel_Exist.inc_nSignal_Gold();
    if (truth_isSignalSilver1) nCut_secEndPoint_Michel_Exist.inc_nSignal_Silver1();
    
    if( Cut_Particle_None == 1) return false;
    nCut_Particle_None.inc_nEvent();
    if (truth_isSignal) nCut_Particle_None.inc_nSignal();
    if (truth_isSignalGold) nCut_Particle_None.inc_nSignal_Gold();
    if (truth_isSignalSilver1) nCut_Particle_None.inc_nSignal_Silver1();
    
    if( Cut_Proton_None == 1) return false;
    nCut_Proton_None.inc_nEvent();
    if (truth_isSignal)  nCut_Proton_None.inc_nSignal();
    if (truth_isSignalGold) nCut_Proton_None.inc_nSignal_Gold();
    if (truth_isSignalSilver1) nCut_Proton_None.inc_nSignal_Silver1();
    
    // Find Best Proton in Reco
    indRecoProton = findBestProton();   
    if ( applyProtonScore && (CCProtonPi0_proton_score[indRecoProton] < minProtonScore) ) return false;
    nCut_Proton_Score.inc_nEvent();
    if (truth_isSignal) nCut_Proton_Score.inc_nSignal();
    if (truth_isSignalGold) nCut_Proton_Score.inc_nSignal_Gold();
    if (truth_isSignalSilver1) nCut_Proton_Score.inc_nSignal_Silver1();
    
    pFilter_Status->Fill(preFilter_Result);
    pFilter_RejectedEnergy->Fill(preFilter_rejectedEnergy);
    
    if( Cut_PreFilter_Pi0 == 1) return false;
    nCut_PreFilter_Pi0.inc_nEvent();
    if (truth_isSignal) nCut_PreFilter_Pi0.inc_nSignal();
    if (truth_isSignalGold) nCut_PreFilter_Pi0.inc_nSignal_Gold();
    if (truth_isSignalSilver1) nCut_PreFilter_Pi0.inc_nSignal_Silver1();
    
    if( Cut_VtxBlob == 1) return false;
    nCut_VtxBlob.inc_nEvent();
    if (truth_isSignal) nCut_VtxBlob.inc_nSignal();
    if (truth_isSignalGold)  nCut_VtxBlob.inc_nSignal_Gold();
    if (truth_isSignalSilver1)  nCut_VtxBlob.inc_nSignal_Silver1();
    
    if( Cut_ConeBlobs == 1) return false;
    nCut_ConeBlobs.inc_nEvent();
    if (truth_isSignal) nCut_ConeBlobs.inc_nSignal();
    if (truth_isSignalGold) nCut_ConeBlobs.inc_nSignal_Gold();
    if (truth_isSignalSilver1)  nCut_ConeBlobs.inc_nSignal_Silver1();
    
    if( pi0_E == SENTINEL) return false;
    nCut_Other.inc_nEvent();
    if (truth_isSignal) nCut_Other.inc_nSignal();
    if (truth_isSignalGold) nCut_Other.inc_nSignal_Gold();
    if (truth_isSignalSilver1) nCut_Other.inc_nSignal_Silver1();
    
    if( applyPhotonDistance && isPhotonDistanceLow()) return false;
    nCut_PhotonDistanceLow.inc_nEvent();
    if (truth_isSignal) nCut_PhotonDistanceLow.inc_nSignal();
    if (truth_isSignalGold) nCut_PhotonDistanceLow.inc_nSignal_Gold();
    if (truth_isSignalSilver1) nCut_PhotonDistanceLow.inc_nSignal_Silver1();
    
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
    
    efficiencyBase_Total = nCut_Vertex_Not_Fiducial.get_nSignal();
    efficiencyBase_Gold = nCut_Vertex_Not_Fiducial.get_nSignal_Gold();
    efficiencyBase_Silver1 = nCut_Vertex_Not_Fiducial.get_nSignal_Silver1();
    
    cout<<">> Writing "<<cutFile<<endl;
    
    writeCutTableHeader();
 
    for( unsigned int i = 0; i < nCutVector.size(); i++){
        
        efficiency_Total = getCutEfficiency(nCutVector[i].get_nSignal(),efficiencyBase_Total);
        efficiency_Gold = getCutEfficiency(nCutVector[i].get_nSignal_Gold(),efficiencyBase_Gold);
        efficiency_Silver1 = getCutEfficiency(nCutVector[i].get_nSignal_Silver1(),efficiencyBase_Silver1);
        
        purity_Total = getCutPurity(nCutVector[i].get_nSignal(),nCutVector[i].get_nEvent());
        purity_Gold = getCutPurity(nCutVector[i].get_nSignal_Gold(),nCutVector[i].get_nEvent());
        purity_Silver1 = getCutPurity(nCutVector[i].get_nSignal_Silver1(),nCutVector[i].get_nEvent());
        
        cutText.unsetf( std::ios::floatfield ); 
        cutText.width(35); cutText<<nCutVector[i].get_Name()<<" | ";
        cutText.width(12); cutText<<nCutVector[i].get_nEvent()<<" | ";
        
        // Total Signal
        cutText.width(12); cutText<<nCutVector[i].get_nSignal()<<" | ";

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
        cutText.width(12); cutText<<nCutVector[i].get_nSignal_Gold()<<" | ";
        
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
        cutText.width(12); cutText<<nCutVector[i].get_nSignal_Silver1()<<" | ";
        
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
        
        cutText<<endl;
    }  
}

void CCProtonPi0::get_pID_Stats()
{
    cout<<"=== Calculationg pID Statistics ==="<<endl;
    
    string rootDir = "Output/RootFiles/Interaction.root";
    TFile* f_Root = new TFile(rootDir.c_str());
    
    TH1D* pID_proton  = (TH1D*)f_Root->Get("pID_proton");
    TH1D* pID_piplus  = (TH1D*)f_Root->Get("pID_piplus");
    TH1D* pID_piminus = (TH1D*)f_Root->Get("pID_piminus");
    TH1D* pID_other   = (TH1D*)f_Root->Get("pID_other");
    
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


void CCProtonPi0::fillInteractionTrue(int indProton)
{
    int indPion;
    double invMass_true;
    double invMass_reco;
    
    int_channel->Fill(mc_intType);
    
    if( truth_N_pi0 == 1 && CCProtonPi0_trajProtonProngPDG[indProton] == 2212){
        indPion = getBestPi0();
        invMass_true = calcDeltaInvariantMass(  mc_FSPartPx[indPion], 
                                                mc_FSPartPy[indPion], 
                                                mc_FSPartPz[indPion],
                                                mc_FSPartE[indPion],
                                                CCProtonPi0_trajProtonProngPx[indProton],
                                                CCProtonPi0_trajProtonProngPy[indProton],
                                                CCProtonPi0_trajProtonProngPz[indProton],
                                                CCProtonPi0_trajProtonProngEnergy[indProton]);
        
        invMass_reco = calcDeltaInvariantMass(  pi0_px, 
                                                pi0_py, 
                                                pi0_pz,
                                                pi0_E,
                                                CCProtonPi0_proton_px[indProton],
                                                CCProtonPi0_proton_py[indProton],
                                                CCProtonPi0_proton_pz[indProton],
                                                CCProtonPi0_proton_E[indProton]);
        
        deltaInvMass_mc->Fill(invMass_true);
        deltaInvMass_reco_mc->Fill(invMass_reco,invMass_true);
        deltaInvMass_error->Fill(Data_Functions::getError(invMass_true,invMass_reco));
    }
    
    beamEnergy_mc->Fill(mc_incomingE);
    
//     beamEnergy_reco->Fill(Erec);
//     beamEnergy_error->Fill( Data_Functions::getError(mc_incomingE,Erec) );
//     beamEnergy_reco_mc->Fill(Erec,mc_incomingE);
    
    q2_mc->Fill(mc_Q2 / mevSq_to_gevSq);
//     q2_reco->Fill(Q2/ mevSq_to_gevSq);
//     q2_error->Fill( Data_Functions::getError(mc_Q2,Q2) );
//     q2_reco_mc->Fill(Q2/mevSq_to_gevSq,mc_Q2 /mevSq_to_gevSq);

    vertex_z_true->Fill(mc_vtx[2]);
    vertex_z_error->Fill( Data_Functions::getError(mc_vtx[2],CCProtonPi0_vtx[2]) );
    vertex_z_reco_mc->Fill(CCProtonPi0_vtx[2],mc_vtx[2]);
    
    vertex_x_y_true->Fill(mc_vtx[0],mc_vtx[1]);
    

}

void CCProtonPi0::fillInteractionReco(int indProton)
{
    double invMass_reco;
    
    invMass_reco = calcDeltaInvariantMass(  pi0_px, 
                                            pi0_py, 
                                            pi0_pz,
                                            pi0_E,
                                            CCProtonPi0_proton_px[indProton],
                                            CCProtonPi0_proton_py[indProton],
                                            CCProtonPi0_proton_pz[indProton],
                                            CCProtonPi0_proton_E[indProton]);
                                            
                                                
    deltaInvMass_reco->Fill(invMass_reco);
    
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
    
    // File Locations
    rootDir =   "Output/RootFiles/Interaction.root";
    plotDir =   "Output/Plots/Interaction/";
    
    max_nFSPart = 10;
    cout<<"\tRoot File: "<<rootDir<<endl;
    cout<<"\tPlot Output Folder: "<<plotDir<<endl;
    
    // Create Root File 
    f = new TFile(rootDir.c_str(),"RECREATE");
    
    // -------------------------------------------------------------------------
    //     Initialization
    //--------------------------------------------------------------------------
    // Default Beam Configuration
    beam_p3.SetXYZ(1.0,1.0,1.0);
    beam_p3.SetPhi(-1.554);
    beam_p3.SetTheta(0.059);
    
    // Cut Numbers
    nCut_All.set_Name("All");
    nCut_Event_Not_Plausible.set_Name("Event_Not_Plausible");
    nCut_Event_Has_BadObject.set_Name("Event_Has_BadObjec");
    nCut_Vertex_None.set_Name("Vertex_None");
    nCut_Vertex_Null.set_Name("Vertex_Null");
    nCut_Vertex_Not_Reconstructable.set_Name("Vertex_Not_Reconstructable"); 
    nCut_Vertex_Not_Fiducial.set_Name("Vertex_Not_Fiducial");    
    nCut_Vertex_Michel_Exist.set_Name("Vertex_Michel_Exist");           
    nCut_nProngs.set_Name("nProngs");
    nCut_Muon_None.set_Name("Muon_None");              
    nCut_Muon_Not_Plausible.set_Name("Muon_Not_Plausible");
    nCut_Muon_Score_Low.set_Name("Muon_Score_Low");
    nCut_Muon_Charge.set_Name("Muon_Charge");
    nCut_EndPoint_Michel_Exist.set_Name("EndPoint_Michel_Exist");
    nCut_secEndPoint_Michel_Exist.set_Name("secEndPoint_Michel_Exist");
    nCut_Particle_None.set_Name("Particle_None");
    nCut_Proton_None.set_Name("Proton_None");            
    nCut_Proton_Score.set_Name("Proton_Score");
    nCut_PreFilter_Pi0.set_Name("PreFilter_Pi0");
    nCut_VtxBlob.set_Name("VtxBlob");
    nCut_ConeBlobs.set_Name("ConeBlobs");
    nCut_PhotonDistanceLow.set_Name("PhotonDistanceLow");
    nCut_Other.set_Name("Other");
    
    cout<<"Done!"<<endl;
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
    cout<<"Text Files for Output:"<<endl;
    // Open Readme File
    readmeFile = "Output/TextFiles/readme.txt";
    readme.open( readmeFile.c_str() );
    if( !readme.is_open() ){
        cerr<<"Cannot open output text file: "<<readmeFile<<endl;
        exit(1);
    }else{
        cout<<"\t"<<readmeFile<<endl;
    }
    writeReadme();
    
    // Open Cut File
    cutFile = "Output/TextFiles/CutTable";
    cutFile = cutFile + "_" + channelTag + ".txt";
    cutText.open( cutFile.c_str() );
    if( !cutText.is_open() ){
        cerr<<"Cannot open output text file: "<<cutFile<<endl;
        exit(1);
    }else{
        cout<<"\t"<<cutFile<<endl;
    }
    
    // Open Fail File
    failFile = "Output/TextFiles/FailChecks.txt";
    failText.open( failFile.c_str() );
    if( !failText.is_open() ){
        cerr<<"Cannot open output text file: "<<failFile<<endl;
        exit(1);
    }else{
        cout<<"\t"<<failFile<<endl;
    }
    
    // Open Arachne RoundupFile
    roundupFile = "Output/TextFiles/ArachneRoundup.txt";
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



#endif //CCProtonPi0_cpp
