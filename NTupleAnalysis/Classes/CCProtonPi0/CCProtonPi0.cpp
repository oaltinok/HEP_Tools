/*
    See CCProtonPi0.h header for Class Information
*/

#ifndef CCProtonPi0_cpp
#define CCProtonPi0_cpp

#include "CCProtonPi0.h"

using namespace std;

void CCProtonPi0::specifyRunTime()
{
    // Control Flow
    isDataAnalysis = true;
    isMC = true;
    applyProtonScore = true;
    minProtonScore = 0.0;
    is_pID_Studies = false;

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
    
    if(m_ActivateMC){
        fChain->SetBranchStatus("truth_*",1);
        fChain->SetBranchStatus("mc_*",1);
        fChain->SetBranchStatus("CCProtonPi0_traj*",1);
    }
    
    fChain->SetBranchStatus("mgg",1);
    
    if(m_ActivateInteraction){
        fChain->SetBranchStatus("CCProtonPi0_vtx*",1);
    }
    
    if(m_ActivateMuon){
        fChain->SetBranchStatus("CCProtonPi0_muon_*",1);
    }
    
    if(m_ActivateProton){
        fChain->SetBranchStatus("CCProtonPi0_proton_*",1);
    }
    
    if(m_ActivatePi0){
        fChain->SetBranchStatus("pi0_*",1);
    }
   

    // Cut Statistics
    double nAll = 0;

    double nCut_Vertex_None = 0;
    double nCut_Vertex_Null = 0;
    double nCut_Vertex_Not_Reconstructable = 0; 
    double nCut_Vertex_Not_Fiducial = 0;    
    double nCut_Vertex_Michel_Exist = 0;           
    double nCut_Muon_None = 0;              
    double nCut_Muon_Score_Low = 0;
    double nCut_EndPoint_Michel_Exist = 0;
    double nCut_secEndPoint_Michel_Exist = 0;
    double nCut_Proton_None = 0;            
    double nCut_Proton_Score = 0;
    double nCut_PreFilter_Pi0 = 0;
    double nCut_Reco_Muon_NoProblem = 0;
    double nCut_Reco_Pi0_NoProblem = 0;
    
    double nSignal = 0;
    double nSignal_Vertex_None = 0;
    double nSignal_Vertex_Null = 0;
    double nSignal_Vertex_Not_Reconstructable = 0; 
    double nSignal_Vertex_Not_Fiducial = 0;    
    double nSignal_Vertex_Michel_Exist = 0;           
    double nSignal_Muon_None = 0;              
    double nSignal_Muon_Score_Low = 0;
    double nSignal_EndPoint_Michel_Exist = 0;
    double nSignal_secEndPoint_Michel_Exist = 0;
    double nSignal_Proton_None = 0;            
    double nSignal_Proton_Score = 0;
    double nSignal_PreFilter_Pi0 = 0;
    double nSignal_Reco_Muon_NoProblem = 0;
    double nSignal_Reco_Pi0_NoProblem = 0;
    
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
    

    
    int indRecoProton;

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
        // Plausibility Cuts
        //----------------------------------------------------------------------
        if (mc_nFSPart > max_nFSPart) continue;

        //----------------------------------------------------------------------
        // Count Events after each Reconstruction Cut
        //----------------------------------------------------------------------
        
        // Count All Events before Cuts
        nAll++;
        if (truth_isSignal) nSignal++;
        
        if( Cut_Vertex_None == 1) continue;
        nCut_Vertex_None++;
        if(truth_isSignal) nSignal_Vertex_None ++;
        
        if( Cut_Vertex_Null == 1) continue;
        nCut_Vertex_Null++;
        if(truth_isSignal) nSignal_Vertex_Null++;
        
        if( Cut_Vertex_Not_Reconstructable == 1) continue;
        nCut_Vertex_Not_Reconstructable++;
        if(truth_isSignal) nSignal_Vertex_Not_Reconstructable++;
        
        if( Cut_Vertex_Not_Fiducial == 1) continue;
        nCut_Vertex_Not_Fiducial++;
        if(truth_isSignal) nSignal_Vertex_Not_Fiducial++;
        
        if( Cut_Muon_None == 1) continue;
        nCut_Muon_None++;
        if(truth_isSignal) nSignal_Muon_None++;
        
        if( Cut_Muon_Score_Low == 1) continue;
        nCut_Muon_Score_Low++;
        if(truth_isSignal) nSignal_Muon_Score_Low++;
        
        if( Cut_Vertex_Michel_Exist == 1) continue;
        nCut_Vertex_Michel_Exist++;
        if(truth_isSignal) nSignal_Vertex_Michel_Exist++;

        if( Cut_EndPoint_Michel_Exist == 1) continue;
        nCut_EndPoint_Michel_Exist++;
        if(truth_isSignal) nSignal_EndPoint_Michel_Exist++;
 
        if( Cut_secEndPoint_Michel_Exist == 1) continue;
        nCut_secEndPoint_Michel_Exist++;
        if(truth_isSignal) nSignal_secEndPoint_Michel_Exist++;
        
        if( Cut_Proton_None == 1 || CCProtonPi0_proton_E[0] == -1) continue;
        nCut_Proton_None++;
        if(truth_isSignal) nSignal_Proton_None++;
        
        // Find Best Proton in Reco
        indRecoProton = findBestProton();
        
        if ( applyProtonScore && (CCProtonPi0_proton_score[indRecoProton] < minProtonScore) ) continue;
        nCut_Proton_Score++;
        if(truth_isSignal) nSignal_Proton_Score++;
        
        if( Cut_PreFilter_Pi0 == 1) continue;
        nCut_PreFilter_Pi0++;
        if(truth_isSignal) nSignal_PreFilter_Pi0++;
        
        if(CCProtonPi0_muon_pz == 0){
            cout<<"Failed Reco Muon E = "<<CCProtonPi0_muon_E<<endl;
            continue;
        };
        nCut_Reco_Muon_NoProblem++;
        if(truth_isSignal) nSignal_Reco_Muon_NoProblem++;
        
        if(pi0_E == -1) continue;
        nCut_Reco_Pi0_NoProblem++;
        if(truth_isSignal) nSignal_Reco_Pi0_NoProblem++;
        
        

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
            for(int i = 0; i < 10; i++){
                if(CCProtonPi0_proton_score[i] == -1) break;
                if(CCProtonPi0_trajProtonProngPDG[i] == 2212){
                    pID_proton->Fill(CCProtonPi0_proton_score[i]);
                }else if(CCProtonPi0_trajProtonProngPDG[i] == 211){
                    pID_piplus->Fill(CCProtonPi0_proton_score[i]);
                }else if(CCProtonPi0_trajProtonProngPDG[i] == -211){
                    pID_piminus->Fill(CCProtonPi0_proton_score[i]);
                }else{
                    pID_other->Fill(CCProtonPi0_proton_score[i]);
                }
            }
        }
        
        
        if ( isDataAnalysis){
            //------------------------------------------------------------------
            // Fill Particles
            //------------------------------------------------------------------
            if( isMC ){
                fillInteractionTrue();
                fillMuonTrue();
                fillProtonTrue();
                fillPionTrue();
            }
            
            // Fill Reconstructed Information
//             fillInteractionReco();
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
    
    if(is_pID_Studies) get_pID_Stats();
    
    failText<<"N(Muon Charge Diff) = "<<nMuonChargeDiff<<endl;
    failText<<"N(AntiMuon) = "<<nAntiMuon<<endl;
    failText<<"N(n0Pi0) = "<<n0Pi0<<" | N(n0Pi0_Test) = "<<n0Pi0_Test<<endl;
    failText<<"N(n1Pi0) = "<<n1Pi0<<" | N(n1Pi0_Test) = "<<n1Pi0_Test<<endl;
    failText<<"N(nMultPi0) = "<<nMultPi0<<" | N(nMultPi0_Test) = "<<nMultPi0_Test<<endl;

    
    cout<<">> Writing "<<cutFile<<endl;
    cutText<<"nAll                          "<<nAll<<endl;
    cutText<<"Cut_Vertex_None               "<<nCut_Vertex_None<<endl;
    cutText<<"Cut_Vertex_Null               "<<nCut_Vertex_Null<<endl;
    cutText<<"Cut_Vertex_Not_Reconstructable     "<<nCut_Vertex_Not_Reconstructable<<endl;
    cutText<<"Cut_Vertex_Not_Fiducial       "<<nCut_Vertex_Not_Fiducial<<endl;
    cutText<<"Cut_Muon_None                 "<<nCut_Muon_None<<endl;
    cutText<<"Cut_Muon_Score_Low            "<<nCut_Muon_Score_Low<<endl; 
    cutText<<"Cut_Vertex_Michel_Exist       "<<nCut_Vertex_Michel_Exist<<endl;
    cutText<<"Cut_EndPoint_Michel_Exist     "<<nCut_EndPoint_Michel_Exist<<endl;
    cutText<<"Cut_secEndPoint_Michel_Exist  "<<nCut_secEndPoint_Michel_Exist<<endl;
    cutText<<"Cut_Proton_None               "<<nCut_Proton_None<<endl;
    cutText<<"Cut_Proton_Score              "<<nCut_Proton_Score<<endl;
    cutText<<"Cut_PreFilter_Pi0             "<<nCut_PreFilter_Pi0<<endl;
    cutText<<"Cut_Reco_Muon_NoProblem       "<<nCut_Reco_Muon_NoProblem<<endl;
    cutText<<"Cut_Reco_Pi0_NoProblem       "<<nCut_Reco_Pi0_NoProblem<<endl;
    cutText<<endl;
    cutText<<"nSignal                          "<<nSignal<<endl;
    cutText<<"Signal_Vertex_None               "<<nSignal_Vertex_None<<endl;
    cutText<<"Signal_Vertex_Null               "<<nSignal_Vertex_Null<<endl;
    cutText<<"Signal_Vertex_Not_Reconstructable     "<<nSignal_Vertex_Not_Reconstructable<<endl;
    cutText<<"Signal_Vertex_Not_Fiducial       "<<nSignal_Vertex_Not_Fiducial<<endl;
    cutText<<"Signal_Muon_None                 "<<nSignal_Muon_None<<endl;
    cutText<<"Signal_Muon_Score_Low            "<<nSignal_Muon_Score_Low<<endl; 
    cutText<<"Signal_Vertex_Michel_Exist       "<<nSignal_Vertex_Michel_Exist<<endl;
    cutText<<"Signal_EndPoint_Michel_Exist     "<<nSignal_EndPoint_Michel_Exist<<endl;
    cutText<<"Signal_secEndPoint_Michel_Exist  "<<nSignal_secEndPoint_Michel_Exist<<endl;
    cutText<<"Signal_Proton_None               "<<nSignal_Proton_None<<endl;
    cutText<<"Signal_Proton_Score              "<<nSignal_Proton_Score<<endl;
    cutText<<"Signal_PreFilter_Pi0             "<<nSignal_PreFilter_Pi0<<endl;
    cutText<<"Signal_Reco_Muon_NoProblem       "<<nSignal_Reco_Muon_NoProblem<<endl;
    cutText<<"Signal_Reco_Pi0_NoProblem       "<<nSignal_Reco_Pi0_NoProblem<<endl;
    
    
    // Write the Root Files
    write_RootFile();           //CCProtonPi0
    muon.write_RootFile();
    proton.write_RootFile();
    pion.write_RootFile();
    
    
    closeFiles();
    
    
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
//
//     Specific Functions
//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

void CCProtonPi0::get_pID_Stats()
{
    cout<<"Calculating pID Statististics"<<endl;
    double nProton = 0;
    double nTotalProton = 0;
    double nCapturedEvents = 0;
    double purity;
    double efficiency;
    //Get Total Proton
    for(int i = 39; i >= 0; i--){
        nTotalProton = nTotalProton + pID_proton->GetBinContent(i);
    }
    cout<<"nTotalProton = "<<nTotalProton<<endl;
    
    for(int i = 39; i >= 0; i--){
        nProton = nProton + pID_proton->GetBinContent(i);
        nCapturedEvents =   nCapturedEvents +
                            pID_proton->GetBinContent(i) +
                            pID_piplus->GetBinContent(i) +
                            pID_piminus->GetBinContent(i) +
                            pID_other->GetBinContent(i);
                            
        purity = nProton / nCapturedEvents;
        efficiency = nProton / nTotalProton;
        cout<<"pID = "<<0.025*i<<" Purity = "<<purity<<" Efficiency "<<efficiency<<endl;
    }
}

void CCProtonPi0::fillInteractionTrue()
{
    beamEnergy_mc->Fill(mc_incomingE);
//     beamEnergy_reco->Fill(Erec);
//     beamEnergy_error->Fill( Data_Functions::getError(mc_incomingE,Erec) );
//     beamEnergy_reco_mc->Fill(Erec,mc_incomingE);
    
    q2_mc->Fill(mc_Q2 / mevSq_to_gevSq);
//     q2_reco->Fill(Q2/ mevSq_to_gevSq);
//     q2_error->Fill( Data_Functions::getError(mc_Q2,Q2) );
//     q2_reco_mc->Fill(Q2/mevSq_to_gevSq,mc_Q2 /mevSq_to_gevSq);

    vertex_z_true->Fill(mc_vtx[2]);
    vertex_z_reco->Fill(CCProtonPi0_vtx[2]);
    vertex_z_error->Fill( Data_Functions::getError(mc_vtx[2],CCProtonPi0_vtx[2]) );
    vertex_z_reco_mc->Fill(CCProtonPi0_vtx[2],mc_vtx[2]);
    
    vertex_x_y_true->Fill(mc_vtx[0],mc_vtx[1]);
    vertex_x_y_reco->Fill(CCProtonPi0_vtx[0],CCProtonPi0_vtx[1]);
    
    int_channel->Fill(mc_intType);
    
    mgg_reco->Fill(mgg);

}

void CCProtonPi0::fillInteractionReco()
{
// Reserved for Future Version
}


void CCProtonPi0::initVariables()
{
    cout<<"Initializing Interaction"<<endl;
    
    channelTag = "Test";
    
    // File Locations
    rootDir =   "Output/RootFiles/Interaction.root";
    plotDir =   "Output/Plots/Interaction/";
    
    
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
    
    max_nFSPart = 15;
    maxBeamEnergy = 20000; //MeV
    
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

    
}

void CCProtonPi0::writeReadme()
{
    readme<<"Test"<<endl;
}

/*
--------------------------------------------------------------------------------
 Beam Energy Cu: isBeamEnergyLow(double maxEnergy)
    Incoming Beam Energy must be lower than a maximum Energy
--------------------------------------------------------------------------------
*/
bool CCProtonPi0::isBeamEnergyLow(double maxEnergy)
{
    if(mc_incomingE > maxEnergy){
        return false;
    }else{
        return true;
    }

}

/*
--------------------------------------------------------------------------------
 findParticle:
    Returns the array indice of the Final State Particle for given PDG
--------------------------------------------------------------------------------
*/
int CCProtonPi0::findTrueParticle(int targetPDG)
{
    int current_ind = -1;
    double current_P = 0;
    
    TVector3 p3;
    
    for(int i = 0; i < mc_nFSPart && i < max_nFSPart; i++ ){
    
        // Momentum of Particle (No Particle at rest)
        // Update Particle Momentum until you find the fastest particle
        if( mc_FSPartPDG[i] == targetPDG){
            p3.SetXYZ(mc_FSPartPx[i],mc_FSPartPy[i],mc_FSPartPz[i]);
            if (current_P < p3.Mag()){
                current_P = p3.Mag();
                current_ind = i;
            }
        }
    }
    
    if (current_P > 0){
        return current_ind;
    }else{
        return -1;
    }
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

void CCProtonPi0::fillProtonTrue()
{
    int ind = proton.ind;
    
    // Fill 4-Momentum
    proton.set_p4(  CCProtonPi0_trajProtonProngPx[ind],
                    CCProtonPi0_trajProtonProngPy[ind],
                    CCProtonPi0_trajProtonProngPz[ind],
                    -1.0, 
                    true);
       
    // set Angle wrt Beam
    proton.set_angleBeam(beam_p3, true);
    
    // set Angle wrt Muon
    proton.set_angleMuon(muon, true);
    
}

void CCProtonPi0::fillProtonReco(int ind)
{
    // Set Particle Score
    proton.particleScore = CCProtonPi0_proton_score[ind];
    
    // Fill 4-Momentum
    proton.set_p4(  CCProtonPi0_proton_px[ind],
                    CCProtonPi0_proton_py[ind],
                    CCProtonPi0_proton_pz[ind],
                    CCProtonPi0_proton_E[ind],
                    false);
                    
    // set Angle wrt Beam
    proton.set_angleBeam(beam_p3, false);
    
    // set Angle wrt Muon
    proton.set_angleMuon(muon, false);
}

int CCProtonPi0::findBestProton()
{
    double tempScore = CCProtonPi0_proton_score[0];
    int tempInd = 0;
    
    for( int i = 1; i < 10; i++){
        if( CCProtonPi0_proton_score[i] == -1 ) break;
        if( CCProtonPi0_proton_score[i] > tempScore){
            tempScore = CCProtonPi0_proton_score[i];
            tempInd = i;
        }
    }
    
    return tempInd;

}

bool CCProtonPi0::isProtonShort(int ind)
{
    const double protonMass = 938; 
    const double minProtonKE = 120;
    double protonKE;
    double protonE;
    
    protonE = mc_FSPartE[ind];
    protonKE =  protonE - protonMass;
    if( protonKE < minProtonKE){
        return true;
    }else{
        return false;
    }
    
}


void CCProtonPi0::fillPionReco()
{
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
}

void CCProtonPi0::fillPionTrue()
{
    double P_true;
    double P_reco;
    
    P_reco = HEP_Functions::calcMomentum(pi0_px,pi0_py,pi0_pz);
   
    // Momentum Information with True Pi0 Count
    if(truth_N_pi0 == 0){
        pion.P_reco_0Pi0->Fill(P_reco);
    }else if(truth_N_pi0 == 1){
        pion.P_reco_1Pi0->Fill(P_reco);
        P_true = getBestPi0Momentum();
        pion.P_reco_mc_1Pi0->Fill(P_reco,P_true);
    }else if(truth_N_pi0 > 1){
        pion.P_reco_MultPi0->Fill(P_reco);
        P_true = getBestPi0Momentum();
        pion.P_reco_mc_MultPi0->Fill(P_reco,P_true);
    }
}

// Loops over all FS Particles and returns the most energetic pi0
double CCProtonPi0::getBestPi0Momentum()
{
    TVector3 p3;
    double tempP = 0;
    
    for(int i = 0; i < mc_nFSPart && i < max_nFSPart; i++ ){
        if( mc_FSPartPDG[i] == 111){
            p3.SetXYZ(mc_FSPartPx[i],mc_FSPartPy[i],mc_FSPartPz[i]);
            if(p3.Mag() > tempP) tempP = p3.Mag();
        }
    }
    return tempP;
}


void CCProtonPi0::fillMuonTrue()
{

    // Fill 4-Momentum
    muon.set_p4(    CCProtonPi0_trajMuonProngPx,
                    CCProtonPi0_trajMuonProngPy,
                    CCProtonPi0_trajMuonProngPz,
                    -1.0, 
                    true);
       
    // set Angle wrt Beam
    muon.set_angleBeam(beam_p3, true);
    
    // set Angle wrt Muon
    muon.set_angleMuon(muon, true);
    
}

void CCProtonPi0::fillMuonReco()
{
    // Set Particle Score
    muon.particleScore = CCProtonPi0_muon_muScore;
    
    // Fill 4-Momentum
    muon.set_p4(    CCProtonPi0_muon_px,
                    CCProtonPi0_muon_py,
                    CCProtonPi0_muon_pz,
                    CCProtonPi0_muon_E,
                    false);
    
    // set Angle wrt Beam
    muon.set_angleBeam(beam_p3, false);
    
    // set Angle wrt Muon
    muon.set_angleMuon(muon, false);
}

void CCProtonPi0::Init(string playlist, TChain* fChain)
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
    cout<<"Initializing Playlist"<<endl;
    
    if( !input_pl.is_open() ){
        cerr<<"Cannot open Playlist File!"<<endl;
        exit(1);
    }else{
        cout<<"\tPlaylist: "<<playlist.c_str()<<endl;
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
   fChain->SetBranchAddress("isMinosMatchTrack", &isMinosMatchTrack, &b_isMinosMatchTrack);
   fChain->SetBranchAddress("isMinosMatchStub", &isMinosMatchStub, &b_isMinosMatchStub);
   fChain->SetBranchAddress("well_fit_vertex", &well_fit_vertex, &b_well_fit_vertex);
   fChain->SetBranchAddress("isBrokenTrack", &isBrokenTrack, &b_isBrokenTrack);
   fChain->SetBranchAddress("is_GoodDirection1", &is_GoodDirection1, &b_is_GoodDirection1);
   fChain->SetBranchAddress("is_GoodPosition1", &is_GoodPosition1, &b_is_GoodPosition1);
   fChain->SetBranchAddress("is_GoodDirection2", &is_GoodDirection2, &b_is_GoodDirection2);
   fChain->SetBranchAddress("is_GoodPosition2", &is_GoodPosition2, &b_is_GoodPosition2);
   fChain->SetBranchAddress("is_GoodBlob1", &is_GoodBlob1, &b_is_GoodBlob1);
   fChain->SetBranchAddress("is_GoodBlob2", &is_GoodBlob2, &b_is_GoodBlob2);
   fChain->SetBranchAddress("Cut_EndPoint_Michel_Exist", &Cut_EndPoint_Michel_Exist, &b_Cut_EndPoint_Michel_Exist);
   fChain->SetBranchAddress("Cut_Muon_Charge", &Cut_Muon_Charge, &b_Cut_Muon_Charge);
   fChain->SetBranchAddress("Cut_Muon_None", &Cut_Muon_None, &b_Cut_Muon_None);
   fChain->SetBranchAddress("Cut_Muon_Score_Low", &Cut_Muon_Score_Low, &b_Cut_Muon_Score_Low);
   fChain->SetBranchAddress("Cut_PreFilter_Pi0", &Cut_PreFilter_Pi0, &b_Cut_PreFilter_Pi0);
   fChain->SetBranchAddress("Cut_Proton_None", &Cut_Proton_None, &b_Cut_Proton_None);
   fChain->SetBranchAddress("Cut_Vertex_Michel_Exist", &Cut_Vertex_Michel_Exist, &b_Cut_Vertex_Michel_Exist);
   fChain->SetBranchAddress("Cut_Vertex_None", &Cut_Vertex_None, &b_Cut_Vertex_None);
   fChain->SetBranchAddress("Cut_Vertex_Not_Fiducial", &Cut_Vertex_Not_Fiducial, &b_Cut_Vertex_Not_Fiducial);
   fChain->SetBranchAddress("Cut_Vertex_Not_Reconstructable", &Cut_Vertex_Not_Reconstructable, &b_Cut_Vertex_Not_Reconstructable);
   fChain->SetBranchAddress("Cut_Vertex_Null", &Cut_Vertex_Null, &b_Cut_Vertex_Null);
   fChain->SetBranchAddress("Cut_secEndPoint_Michel_Exist", &Cut_secEndPoint_Michel_Exist, &b_Cut_secEndPoint_Michel_Exist);
   fChain->SetBranchAddress("anglescan_ncand", &anglescan_ncand, &b_anglescan_ncand);
   fChain->SetBranchAddress("anglescan_ncandx", &anglescan_ncandx, &b_anglescan_ncandx);
   fChain->SetBranchAddress("broken_track_most_us_plane", &broken_track_most_us_plane, &b_broken_track_most_us_plane);
   fChain->SetBranchAddress("g1blob_ncluster", &g1blob_ncluster, &b_g1blob_ncluster);
   fChain->SetBranchAddress("g1blob_ndigit", &g1blob_ndigit, &b_g1blob_ndigit);
   fChain->SetBranchAddress("g1dedx_doublet", &g1dedx_doublet, &b_g1dedx_doublet);
   fChain->SetBranchAddress("g1dedx_empty_plane", &g1dedx_empty_plane, &b_g1dedx_empty_plane);
   fChain->SetBranchAddress("g1dedx_nplane", &g1dedx_nplane, &b_g1dedx_nplane);
   fChain->SetBranchAddress("g1mostevispdg", &g1mostevispdg, &b_g1mostevispdg);
   fChain->SetBranchAddress("g2blob_ncluster", &g2blob_ncluster, &b_g2blob_ncluster);
   fChain->SetBranchAddress("g2blob_ndigit", &g2blob_ndigit, &b_g2blob_ndigit);
   fChain->SetBranchAddress("g2dedx_doublet", &g2dedx_doublet, &b_g2dedx_doublet);
   fChain->SetBranchAddress("g2dedx_empty_plane", &g2dedx_empty_plane, &b_g2dedx_empty_plane);
   fChain->SetBranchAddress("g2dedx_nplane", &g2dedx_nplane, &b_g2dedx_nplane);
   fChain->SetBranchAddress("g2mostevispdg", &g2mostevispdg, &b_g2mostevispdg);
   fChain->SetBranchAddress("n_anchored_long_trk_prongs", &n_anchored_long_trk_prongs, &b_n_anchored_long_trk_prongs);
   fChain->SetBranchAddress("n_anchored_short_trk_prongs", &n_anchored_short_trk_prongs, &b_n_anchored_short_trk_prongs);
   fChain->SetBranchAddress("n_dsp_blob_prongs", &n_dsp_blob_prongs, &b_n_dsp_blob_prongs);
   fChain->SetBranchAddress("n_iso_blob_prongs", &n_iso_blob_prongs, &b_n_iso_blob_prongs);
   fChain->SetBranchAddress("n_iso_trk_prongs", &n_iso_trk_prongs, &b_n_iso_trk_prongs);
   fChain->SetBranchAddress("n_long_tracks", &n_long_tracks, &b_n_long_tracks);
   fChain->SetBranchAddress("n_short_tracks", &n_short_tracks, &b_n_short_tracks);
   fChain->SetBranchAddress("n_startpoint_vertices", &n_startpoint_vertices, &b_n_startpoint_vertices);
   fChain->SetBranchAddress("n_us_muon_clusters", &n_us_muon_clusters, &b_n_us_muon_clusters);
   fChain->SetBranchAddress("n_vtx_michel_views", &n_vtx_michel_views, &b_n_vtx_michel_views);
   fChain->SetBranchAddress("n_vtx_prongs", &n_vtx_prongs, &b_n_vtx_prongs);
   fChain->SetBranchAddress("nblob_anglescan", &nblob_anglescan, &b_nblob_anglescan);
   fChain->SetBranchAddress("nblob_hough", &nblob_hough, &b_nblob_hough);
   fChain->SetBranchAddress("od_energeticTower", &od_energeticTower, &b_od_energeticTower);
   fChain->SetBranchAddress("phys_energy_in_road_downstream_nplanes", &phys_energy_in_road_downstream_nplanes, &b_phys_energy_in_road_downstream_nplanes);
   fChain->SetBranchAddress("phys_energy_in_road_upstream_nplanes", &phys_energy_in_road_upstream_nplanes, &b_phys_energy_in_road_upstream_nplanes);
   fChain->SetBranchAddress("phys_n_dead_discr_pair", &phys_n_dead_discr_pair, &b_phys_n_dead_discr_pair);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_in_prim_track_region", &phys_n_dead_discr_pair_in_prim_track_region, &b_phys_n_dead_discr_pair_in_prim_track_region);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_two_mod_downstream_prim_track", &phys_n_dead_discr_pair_two_mod_downstream_prim_track, &b_phys_n_dead_discr_pair_two_mod_downstream_prim_track);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_two_mod_upstream_prim_vtx", &phys_n_dead_discr_pair_two_mod_upstream_prim_vtx, &b_phys_n_dead_discr_pair_two_mod_upstream_prim_vtx);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_upstream_prim_track_proj", &phys_n_dead_discr_pair_upstream_prim_track_proj, &b_phys_n_dead_discr_pair_upstream_prim_track_proj);
   fChain->SetBranchAddress("phys_vertex_is_fiducial", &phys_vertex_is_fiducial, &b_phys_vertex_is_fiducial);
   fChain->SetBranchAddress("vtx_primary_index", &vtx_primary_index, &b_vtx_primary_index);
   fChain->SetBranchAddress("vtx_primary_multiplicity", &vtx_primary_multiplicity, &b_vtx_primary_multiplicity);
   fChain->SetBranchAddress("vtx_secondary_count", &vtx_secondary_count, &b_vtx_secondary_count);
   fChain->SetBranchAddress("vtx_total_count", &vtx_total_count, &b_vtx_total_count);
   fChain->SetBranchAddress("Filament_Vertex_energy", &Filament_Vertex_energy, &b_Filament_Vertex_energy);
   fChain->SetBranchAddress("RE_energy_ECAL", &RE_energy_ECAL, &b_RE_energy_ECAL);
   fChain->SetBranchAddress("RE_energy_HCAL", &RE_energy_HCAL, &b_RE_energy_HCAL);
   fChain->SetBranchAddress("RE_energy_Tracker", &RE_energy_Tracker, &b_RE_energy_Tracker);
   fChain->SetBranchAddress("RE_photon_dEdx_1", &RE_photon_dEdx_1, &b_RE_photon_dEdx_1);
   fChain->SetBranchAddress("RE_photon_dEdx_2", &RE_photon_dEdx_2, &b_RE_photon_dEdx_2);
   fChain->SetBranchAddress("RE_photon_energy_1", &RE_photon_energy_1, &b_RE_photon_energy_1);
   fChain->SetBranchAddress("RE_photon_energy_2", &RE_photon_energy_2, &b_RE_photon_energy_2);
   fChain->SetBranchAddress("RE_photon_time_1", &RE_photon_time_1, &b_RE_photon_time_1);
   fChain->SetBranchAddress("RE_photon_time_2", &RE_photon_time_2, &b_RE_photon_time_2);
   fChain->SetBranchAddress("RE_scalar", &RE_scalar, &b_RE_scalar);
   fChain->SetBranchAddress("Sphere_Vertex_energy", &Sphere_Vertex_energy, &b_Sphere_Vertex_energy);
   fChain->SetBranchAddress("Vertex_blob_energy", &Vertex_blob_energy, &b_Vertex_blob_energy);
   fChain->SetBranchAddress("dispersedExtraE", &dispersedExtraE, &b_dispersedExtraE);
   fChain->SetBranchAddress("energy_from_mc", &energy_from_mc, &b_energy_from_mc);
   fChain->SetBranchAddress("energy_from_mc_fraction", &energy_from_mc_fraction, &b_energy_from_mc_fraction);
   fChain->SetBranchAddress("energy_from_mc_fraction_of_highest", &energy_from_mc_fraction_of_highest, &b_energy_from_mc_fraction_of_highest);
   fChain->SetBranchAddress("evis_ecal", &evis_ecal, &b_evis_ecal);
   fChain->SetBranchAddress("evis_ecal_u", &evis_ecal_u, &b_evis_ecal_u);
   fChain->SetBranchAddress("evis_ecal_v", &evis_ecal_v, &b_evis_ecal_v);
   fChain->SetBranchAddress("evis_ecal_x", &evis_ecal_x, &b_evis_ecal_x);
   fChain->SetBranchAddress("evis_hcal", &evis_hcal, &b_evis_hcal);
   fChain->SetBranchAddress("evis_hcal_u", &evis_hcal_u, &b_evis_hcal_u);
   fChain->SetBranchAddress("evis_hcal_v", &evis_hcal_v, &b_evis_hcal_v);
   fChain->SetBranchAddress("evis_hcal_x", &evis_hcal_x, &b_evis_hcal_x);
   fChain->SetBranchAddress("evis_nearvtx_total", &evis_nearvtx_total, &b_evis_nearvtx_total);
   fChain->SetBranchAddress("evis_nearvtx_u", &evis_nearvtx_u, &b_evis_nearvtx_u);
   fChain->SetBranchAddress("evis_nearvtx_v", &evis_nearvtx_v, &b_evis_nearvtx_v);
   fChain->SetBranchAddress("evis_nearvtx_x", &evis_nearvtx_x, &b_evis_nearvtx_x);
   fChain->SetBranchAddress("evis_ntgt", &evis_ntgt, &b_evis_ntgt);
   fChain->SetBranchAddress("evis_ntgt_u", &evis_ntgt_u, &b_evis_ntgt_u);
   fChain->SetBranchAddress("evis_ntgt_v", &evis_ntgt_v, &b_evis_ntgt_v);
   fChain->SetBranchAddress("evis_ntgt_x", &evis_ntgt_x, &b_evis_ntgt_x);
   fChain->SetBranchAddress("evis_other", &evis_other, &b_evis_other);
   fChain->SetBranchAddress("evis_total", &evis_total, &b_evis_total);
   fChain->SetBranchAddress("evis_total_u", &evis_total_u, &b_evis_total_u);
   fChain->SetBranchAddress("evis_total_v", &evis_total_v, &b_evis_total_v);
   fChain->SetBranchAddress("evis_total_x", &evis_total_x, &b_evis_total_x);
   fChain->SetBranchAddress("evis_trkr", &evis_trkr, &b_evis_trkr);
   fChain->SetBranchAddress("evis_trkr_u", &evis_trkr_u, &b_evis_trkr_u);
   fChain->SetBranchAddress("evis_trkr_v", &evis_trkr_v, &b_evis_trkr_v);
   fChain->SetBranchAddress("evis_trkr_x", &evis_trkr_x, &b_evis_trkr_x);
   fChain->SetBranchAddress("g1blob_edge_distance", &g1blob_edge_distance, &b_g1blob_edge_distance);
   fChain->SetBranchAddress("g1blob_minsep", &g1blob_minsep, &b_g1blob_minsep);
   fChain->SetBranchAddress("g1blob_vtx_distance", &g1blob_vtx_distance, &b_g1blob_vtx_distance);
   fChain->SetBranchAddress("g1dedx", &g1dedx, &b_g1dedx);
   fChain->SetBranchAddress("g1dedx1", &g1dedx1, &b_g1dedx1);
   fChain->SetBranchAddress("g1dedx_total", &g1dedx_total, &b_g1dedx_total);
   fChain->SetBranchAddress("g1dedx_total1", &g1dedx_total1, &b_g1dedx_total1);
   fChain->SetBranchAddress("g1e", &g1e, &b_g1e);
   fChain->SetBranchAddress("g1g1evis", &g1g1evis, &b_g1g1evis);
   fChain->SetBranchAddress("g1g2evis", &g1g2evis, &b_g1g2evis);
   fChain->SetBranchAddress("g1gmevis", &g1gmevis, &b_g1gmevis);
   fChain->SetBranchAddress("g1mostevisfrac", &g1mostevisfrac, &b_g1mostevisfrac);
   fChain->SetBranchAddress("g1muevis", &g1muevis, &b_g1muevis);
   fChain->SetBranchAddress("g1neutronevis", &g1neutronevis, &b_g1neutronevis);
   fChain->SetBranchAddress("g1otherevis", &g1otherevis, &b_g1otherevis);
   fChain->SetBranchAddress("g1phi", &g1phi, &b_g1phi);
   fChain->SetBranchAddress("g1pi0evis", &g1pi0evis, &b_g1pi0evis);
   fChain->SetBranchAddress("g1pimevis", &g1pimevis, &b_g1pimevis);
   fChain->SetBranchAddress("g1pipevis", &g1pipevis, &b_g1pipevis);
   fChain->SetBranchAddress("g1protonevis", &g1protonevis, &b_g1protonevis);
   fChain->SetBranchAddress("g1recoecalevis", &g1recoecalevis, &b_g1recoecalevis);
   fChain->SetBranchAddress("g1recohcalevis", &g1recohcalevis, &b_g1recohcalevis);
   fChain->SetBranchAddress("g1recoscalevis", &g1recoscalevis, &b_g1recoscalevis);
   fChain->SetBranchAddress("g1recotrkrevis", &g1recotrkrevis, &b_g1recotrkrevis);
   fChain->SetBranchAddress("g1sharedevis", &g1sharedevis, &b_g1sharedevis);
   fChain->SetBranchAddress("g1theta", &g1theta, &b_g1theta);
   fChain->SetBranchAddress("g1totalevis", &g1totalevis, &b_g1totalevis);
   fChain->SetBranchAddress("g2blob_edge_distance", &g2blob_edge_distance, &b_g2blob_edge_distance);
   fChain->SetBranchAddress("g2blob_minsep", &g2blob_minsep, &b_g2blob_minsep);
   fChain->SetBranchAddress("g2blob_vtx_distance", &g2blob_vtx_distance, &b_g2blob_vtx_distance);
   fChain->SetBranchAddress("g2dedx", &g2dedx, &b_g2dedx);
   fChain->SetBranchAddress("g2dedx1", &g2dedx1, &b_g2dedx1);
   fChain->SetBranchAddress("g2dedx_total", &g2dedx_total, &b_g2dedx_total);
   fChain->SetBranchAddress("g2dedx_total1", &g2dedx_total1, &b_g2dedx_total1);
   fChain->SetBranchAddress("g2e", &g2e, &b_g2e);
   fChain->SetBranchAddress("g2g1evis", &g2g1evis, &b_g2g1evis);
   fChain->SetBranchAddress("g2g2evis", &g2g2evis, &b_g2g2evis);
   fChain->SetBranchAddress("g2gmevis", &g2gmevis, &b_g2gmevis);
   fChain->SetBranchAddress("g2mostevisfrac", &g2mostevisfrac, &b_g2mostevisfrac);
   fChain->SetBranchAddress("g2muevis", &g2muevis, &b_g2muevis);
   fChain->SetBranchAddress("g2neutronevis", &g2neutronevis, &b_g2neutronevis);
   fChain->SetBranchAddress("g2otherevis", &g2otherevis, &b_g2otherevis);
   fChain->SetBranchAddress("g2phi", &g2phi, &b_g2phi);
   fChain->SetBranchAddress("g2pi0evis", &g2pi0evis, &b_g2pi0evis);
   fChain->SetBranchAddress("g2pimevis", &g2pimevis, &b_g2pimevis);
   fChain->SetBranchAddress("g2pipevis", &g2pipevis, &b_g2pipevis);
   fChain->SetBranchAddress("g2protonevis", &g2protonevis, &b_g2protonevis);
   fChain->SetBranchAddress("g2recoecalevis", &g2recoecalevis, &b_g2recoecalevis);
   fChain->SetBranchAddress("g2recohcalevis", &g2recohcalevis, &b_g2recohcalevis);
   fChain->SetBranchAddress("g2recoscalevis", &g2recoscalevis, &b_g2recoscalevis);
   fChain->SetBranchAddress("g2recotrkrevis", &g2recotrkrevis, &b_g2recotrkrevis);
   fChain->SetBranchAddress("g2sharedevis", &g2sharedevis, &b_g2sharedevis);
   fChain->SetBranchAddress("g2theta", &g2theta, &b_g2theta);
   fChain->SetBranchAddress("g2totalevis", &g2totalevis, &b_g2totalevis);
   fChain->SetBranchAddress("gamma1_E", &gamma1_E, &b_gamma1_E);
   fChain->SetBranchAddress("gamma1_px", &gamma1_px, &b_gamma1_px);
   fChain->SetBranchAddress("gamma1_py", &gamma1_py, &b_gamma1_py);
   fChain->SetBranchAddress("gamma1_pz", &gamma1_pz, &b_gamma1_pz);
   fChain->SetBranchAddress("gamma2_E", &gamma2_E, &b_gamma2_E);
   fChain->SetBranchAddress("gamma2_px", &gamma2_px, &b_gamma2_px);
   fChain->SetBranchAddress("gamma2_py", &gamma2_py, &b_gamma2_py);
   fChain->SetBranchAddress("gamma2_pz", &gamma2_pz, &b_gamma2_pz);
   fChain->SetBranchAddress("hadronVisibleE", &hadronVisibleE, &b_hadronVisibleE);
   fChain->SetBranchAddress("mgg", &mgg, &b_mgg);
   fChain->SetBranchAddress("muonVisibleE", &muonVisibleE, &b_muonVisibleE);
   fChain->SetBranchAddress("muon_phi", &muon_phi, &b_muon_phi);
   fChain->SetBranchAddress("muon_theta", &muon_theta, &b_muon_theta);
   fChain->SetBranchAddress("muon_thetaX", &muon_thetaX, &b_muon_thetaX);
   fChain->SetBranchAddress("muon_thetaY", &muon_thetaY, &b_muon_thetaY);
   fChain->SetBranchAddress("oangle", &oangle, &b_oangle);
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
   fChain->SetBranchAddress("pi0_px", &pi0_px, &b_pi0_px);
   fChain->SetBranchAddress("pi0_py", &pi0_py, &b_pi0_py);
   fChain->SetBranchAddress("pi0_pz", &pi0_pz, &b_pi0_pz);
   fChain->SetBranchAddress("pienergy", &pienergy, &b_pienergy);
   fChain->SetBranchAddress("piphi", &piphi, &b_piphi);
   fChain->SetBranchAddress("piphib", &piphib, &b_piphib);
   fChain->SetBranchAddress("pitheta", &pitheta, &b_pitheta);
   fChain->SetBranchAddress("pithetab", &pithetab, &b_pithetab);
   fChain->SetBranchAddress("pithetax", &pithetax, &b_pithetax);
   fChain->SetBranchAddress("pithetaxb", &pithetaxb, &b_pithetaxb);
   fChain->SetBranchAddress("pithetay", &pithetay, &b_pithetay);
   fChain->SetBranchAddress("pithetayb", &pithetayb, &b_pithetayb);
   fChain->SetBranchAddress("prim_vtx_smallest_opening_angle", &prim_vtx_smallest_opening_angle, &b_prim_vtx_smallest_opening_angle);
   fChain->SetBranchAddress("time", &time, &b_time);
   fChain->SetBranchAddress("totalIDVisibleE", &totalIDVisibleE, &b_totalIDVisibleE);
   fChain->SetBranchAddress("totalODVisibleE", &totalODVisibleE, &b_totalODVisibleE);
   fChain->SetBranchAddress("totalVisibleE", &totalVisibleE, &b_totalVisibleE);
   fChain->SetBranchAddress("unattachedExtraE", &unattachedExtraE, &b_unattachedExtraE);
   fChain->SetBranchAddress("vtxBlobExtraE", &vtxBlobExtraE, &b_vtxBlobExtraE);
   fChain->SetBranchAddress("vtx_michel_distance", &vtx_michel_distance, &b_vtx_michel_distance);
   fChain->SetBranchAddress("well_fit_vertex_angle", &well_fit_vertex_angle, &b_well_fit_vertex_angle);
   fChain->SetBranchAddress("anglescan_blob_nc_sz", &anglescan_blob_nc_sz, &b_anglescan_blob_nc_sz);
   fChain->SetBranchAddress("anglescan_blob_nc", anglescan_blob_nc, &b_anglescan_blob_nc);
   fChain->SetBranchAddress("anglescan_blob_ncu_sz", &anglescan_blob_ncu_sz, &b_anglescan_blob_ncu_sz);
   fChain->SetBranchAddress("anglescan_blob_ncu", anglescan_blob_ncu, &b_anglescan_blob_ncu);
   fChain->SetBranchAddress("anglescan_blob_ncv_sz", &anglescan_blob_ncv_sz, &b_anglescan_blob_ncv_sz);
   fChain->SetBranchAddress("anglescan_blob_ncv", anglescan_blob_ncv, &b_anglescan_blob_ncv);
   fChain->SetBranchAddress("anglescan_blob_ncx_sz", &anglescan_blob_ncx_sz, &b_anglescan_blob_ncx_sz);
   fChain->SetBranchAddress("anglescan_blob_ncx", anglescan_blob_ncx, &b_anglescan_blob_ncx);
   fChain->SetBranchAddress("anglescan_blob_nd_sz", &anglescan_blob_nd_sz, &b_anglescan_blob_nd_sz);
   fChain->SetBranchAddress("anglescan_blob_nd", anglescan_blob_nd, &b_anglescan_blob_nd);
   fChain->SetBranchAddress("anglescan_blob_ndu_sz", &anglescan_blob_ndu_sz, &b_anglescan_blob_ndu_sz);
   fChain->SetBranchAddress("anglescan_blob_ndu", anglescan_blob_ndu, &b_anglescan_blob_ndu);
   fChain->SetBranchAddress("anglescan_blob_ndv_sz", &anglescan_blob_ndv_sz, &b_anglescan_blob_ndv_sz);
   fChain->SetBranchAddress("anglescan_blob_ndv", anglescan_blob_ndv, &b_anglescan_blob_ndv);
   fChain->SetBranchAddress("anglescan_blob_ndx_sz", &anglescan_blob_ndx_sz, &b_anglescan_blob_ndx_sz);
   fChain->SetBranchAddress("anglescan_blob_ndx", anglescan_blob_ndx, &b_anglescan_blob_ndx);
   fChain->SetBranchAddress("anglescan_cand_nc_sz", &anglescan_cand_nc_sz, &b_anglescan_cand_nc_sz);
   fChain->SetBranchAddress("anglescan_cand_nc", anglescan_cand_nc, &b_anglescan_cand_nc);
   fChain->SetBranchAddress("anglescan_cand_ncu_sz", &anglescan_cand_ncu_sz, &b_anglescan_cand_ncu_sz);
   fChain->SetBranchAddress("anglescan_cand_ncu", anglescan_cand_ncu, &b_anglescan_cand_ncu);
   fChain->SetBranchAddress("anglescan_cand_ncv_sz", &anglescan_cand_ncv_sz, &b_anglescan_cand_ncv_sz);
   fChain->SetBranchAddress("anglescan_cand_ncv", anglescan_cand_ncv, &b_anglescan_cand_ncv);
   fChain->SetBranchAddress("anglescan_cand_ncx_sz", &anglescan_cand_ncx_sz, &b_anglescan_cand_ncx_sz);
   fChain->SetBranchAddress("anglescan_cand_ncx", anglescan_cand_ncx, &b_anglescan_cand_ncx);
   fChain->SetBranchAddress("anglescan_cand_nd_sz", &anglescan_cand_nd_sz, &b_anglescan_cand_nd_sz);
   fChain->SetBranchAddress("anglescan_cand_nd", anglescan_cand_nd, &b_anglescan_cand_nd);
   fChain->SetBranchAddress("anglescan_cand_ndu_sz", &anglescan_cand_ndu_sz, &b_anglescan_cand_ndu_sz);
   fChain->SetBranchAddress("anglescan_cand_ndu", anglescan_cand_ndu, &b_anglescan_cand_ndu);
   fChain->SetBranchAddress("anglescan_cand_ndv_sz", &anglescan_cand_ndv_sz, &b_anglescan_cand_ndv_sz);
   fChain->SetBranchAddress("anglescan_cand_ndv", anglescan_cand_ndv, &b_anglescan_cand_ndv);
   fChain->SetBranchAddress("anglescan_cand_ndx_sz", &anglescan_cand_ndx_sz, &b_anglescan_cand_ndx_sz);
   fChain->SetBranchAddress("anglescan_cand_ndx", anglescan_cand_ndx, &b_anglescan_cand_ndx);
   fChain->SetBranchAddress("anglescan_candx_nc_sz", &anglescan_candx_nc_sz, &b_anglescan_candx_nc_sz);
   fChain->SetBranchAddress("anglescan_candx_nc", anglescan_candx_nc, &b_anglescan_candx_nc);
   fChain->SetBranchAddress("anglescan_candx_nd_sz", &anglescan_candx_nd_sz, &b_anglescan_candx_nd_sz);
   fChain->SetBranchAddress("anglescan_candx_nd", anglescan_candx_nd, &b_anglescan_candx_nd);
   fChain->SetBranchAddress("final_blob_nc_sz", &final_blob_nc_sz, &b_final_blob_nc_sz);
   fChain->SetBranchAddress("final_blob_nc", final_blob_nc, &b_final_blob_nc);
   fChain->SetBranchAddress("final_blob_ncu_sz", &final_blob_ncu_sz, &b_final_blob_ncu_sz);
   fChain->SetBranchAddress("final_blob_ncu", final_blob_ncu, &b_final_blob_ncu);
   fChain->SetBranchAddress("final_blob_ncv_sz", &final_blob_ncv_sz, &b_final_blob_ncv_sz);
   fChain->SetBranchAddress("final_blob_ncv", final_blob_ncv, &b_final_blob_ncv);
   fChain->SetBranchAddress("final_blob_ncx_sz", &final_blob_ncx_sz, &b_final_blob_ncx_sz);
   fChain->SetBranchAddress("final_blob_ncx", final_blob_ncx, &b_final_blob_ncx);
   fChain->SetBranchAddress("final_blob_nd_sz", &final_blob_nd_sz, &b_final_blob_nd_sz);
   fChain->SetBranchAddress("final_blob_nd", final_blob_nd, &b_final_blob_nd);
   fChain->SetBranchAddress("final_blob_ndu_sz", &final_blob_ndu_sz, &b_final_blob_ndu_sz);
   fChain->SetBranchAddress("final_blob_ndu", final_blob_ndu, &b_final_blob_ndu);
   fChain->SetBranchAddress("final_blob_ndv_sz", &final_blob_ndv_sz, &b_final_blob_ndv_sz);
   fChain->SetBranchAddress("final_blob_ndv", final_blob_ndv, &b_final_blob_ndv);
   fChain->SetBranchAddress("final_blob_ndx_sz", &final_blob_ndx_sz, &b_final_blob_ndx_sz);
   fChain->SetBranchAddress("final_blob_ndx", final_blob_ndx, &b_final_blob_ndx);
   fChain->SetBranchAddress("g1dedx_cluster_occupancy_sz", &g1dedx_cluster_occupancy_sz, &b_g1dedx_cluster_occupancy_sz);
   fChain->SetBranchAddress("g1dedx_cluster_occupancy", g1dedx_cluster_occupancy, &b_g1dedx_cluster_occupancy);
   fChain->SetBranchAddress("g2dedx_cluster_occupancy_sz", &g2dedx_cluster_occupancy_sz, &b_g2dedx_cluster_occupancy_sz);
   fChain->SetBranchAddress("g2dedx_cluster_occupancy", g2dedx_cluster_occupancy, &b_g2dedx_cluster_occupancy);
   fChain->SetBranchAddress("hough_blob_nc_sz", &hough_blob_nc_sz, &b_hough_blob_nc_sz);
   fChain->SetBranchAddress("hough_blob_nc", hough_blob_nc, &b_hough_blob_nc);
   fChain->SetBranchAddress("hough_blob_ncu_sz", &hough_blob_ncu_sz, &b_hough_blob_ncu_sz);
   fChain->SetBranchAddress("hough_blob_ncu", hough_blob_ncu, &b_hough_blob_ncu);
   fChain->SetBranchAddress("hough_blob_ncv_sz", &hough_blob_ncv_sz, &b_hough_blob_ncv_sz);
   fChain->SetBranchAddress("hough_blob_ncv", hough_blob_ncv, &b_hough_blob_ncv);
   fChain->SetBranchAddress("hough_blob_ncx_sz", &hough_blob_ncx_sz, &b_hough_blob_ncx_sz);
   fChain->SetBranchAddress("hough_blob_ncx", hough_blob_ncx, &b_hough_blob_ncx);
   fChain->SetBranchAddress("hough_blob_nd_sz", &hough_blob_nd_sz, &b_hough_blob_nd_sz);
   fChain->SetBranchAddress("hough_blob_nd", hough_blob_nd, &b_hough_blob_nd);
   fChain->SetBranchAddress("hough_blob_ndu_sz", &hough_blob_ndu_sz, &b_hough_blob_ndu_sz);
   fChain->SetBranchAddress("hough_blob_ndu", hough_blob_ndu, &b_hough_blob_ndu);
   fChain->SetBranchAddress("hough_blob_ndv_sz", &hough_blob_ndv_sz, &b_hough_blob_ndv_sz);
   fChain->SetBranchAddress("hough_blob_ndv", hough_blob_ndv, &b_hough_blob_ndv);
   fChain->SetBranchAddress("hough_blob_ndx_sz", &hough_blob_ndx_sz, &b_hough_blob_ndx_sz);
   fChain->SetBranchAddress("hough_blob_ndx", hough_blob_ndx, &b_hough_blob_ndx);
   fChain->SetBranchAddress("RE_photon_direction_1_sz", &RE_photon_direction_1_sz, &b_RE_photon_direction_1_sz);
   fChain->SetBranchAddress("RE_photon_direction_1", RE_photon_direction_1, &b_RE_photon_direction_1);
   fChain->SetBranchAddress("RE_photon_direction_2_sz", &RE_photon_direction_2_sz, &b_RE_photon_direction_2_sz);
   fChain->SetBranchAddress("RE_photon_direction_2", RE_photon_direction_2, &b_RE_photon_direction_2);
   fChain->SetBranchAddress("RE_photon_vertex_1_sz", &RE_photon_vertex_1_sz, &b_RE_photon_vertex_1_sz);
   fChain->SetBranchAddress("RE_photon_vertex_1", RE_photon_vertex_1, &b_RE_photon_vertex_1);
   fChain->SetBranchAddress("RE_photon_vertex_2_sz", &RE_photon_vertex_2_sz, &b_RE_photon_vertex_2_sz);
   fChain->SetBranchAddress("RE_photon_vertex_2", RE_photon_vertex_2, &b_RE_photon_vertex_2);
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
   fChain->SetBranchAddress("g1mom", g1mom, &b_g1mom);
   fChain->SetBranchAddress("g2dedx_cluster_energy_sz", &g2dedx_cluster_energy_sz, &b_g2dedx_cluster_energy_sz);
   fChain->SetBranchAddress("g2dedx_cluster_energy", g2dedx_cluster_energy, &b_g2dedx_cluster_energy);
   fChain->SetBranchAddress("g2dedx_rev_cluster_energy_sz", &g2dedx_rev_cluster_energy_sz, &b_g2dedx_rev_cluster_energy_sz);
   fChain->SetBranchAddress("g2dedx_rev_cluster_energy", g2dedx_rev_cluster_energy, &b_g2dedx_rev_cluster_energy);
   fChain->SetBranchAddress("g2mom", g2mom, &b_g2mom);
   fChain->SetBranchAddress("good_mgg_vector_sz", &good_mgg_vector_sz, &b_good_mgg_vector_sz);
   fChain->SetBranchAddress("good_mgg_vector", &good_mgg_vector, &b_good_mgg_vector);
   fChain->SetBranchAddress("mgg_vector_sz", &mgg_vector_sz, &b_mgg_vector_sz);
   fChain->SetBranchAddress("mgg_vector", &mgg_vector, &b_mgg_vector);
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
   fChain->SetBranchAddress("pimom", pimom, &b_pimom);
   fChain->SetBranchAddress("truth_has_physics_event", &truth_has_physics_event, &b_truth_has_physics_event);
   fChain->SetBranchAddress("truth_reco_hasGoodObjects", &truth_reco_hasGoodObjects, &b_truth_reco_hasGoodObjects);
   fChain->SetBranchAddress("truth_reco_isGoodVertex", &truth_reco_isGoodVertex, &b_truth_reco_isGoodVertex);
   fChain->SetBranchAddress("truth_reco_isWellFitVertex", &truth_reco_isWellFitVertex, &b_truth_reco_isWellFitVertex);
   fChain->SetBranchAddress("truth_reco_isFidVol", &truth_reco_isFidVol, &b_truth_reco_isFidVol);
   fChain->SetBranchAddress("truth_reco_isFidVol_smeared", &truth_reco_isFidVol_smeared, &b_truth_reco_isFidVol_smeared);
   fChain->SetBranchAddress("truth_reco_isMinosMatch", &truth_reco_isMinosMatch, &b_truth_reco_isMinosMatch);
   fChain->SetBranchAddress("truth_reco_isBrokenTrack", &truth_reco_isBrokenTrack, &b_truth_reco_isBrokenTrack);
   fChain->SetBranchAddress("truth_isSignal", &truth_isSignal, &b_truth_isSignal);
   fChain->SetBranchAddress("truth_isFidVol", &truth_isFidVol, &b_truth_isFidVol);
   fChain->SetBranchAddress("truth_isPlausible", &truth_isPlausible, &b_truth_isPlausible);
   fChain->SetBranchAddress("truth_N_deltaplus", &truth_N_deltaplus, &b_truth_N_deltaplus);
   fChain->SetBranchAddress("truth_N_gamma", &truth_N_gamma, &b_truth_N_gamma);
   fChain->SetBranchAddress("truth_N_muminus", &truth_N_muminus, &b_truth_N_muminus);
   fChain->SetBranchAddress("truth_N_muplus", &truth_N_muplus, &b_truth_N_muplus);
   fChain->SetBranchAddress("truth_N_neutron", &truth_N_neutron, &b_truth_N_neutron);
   fChain->SetBranchAddress("truth_N_other", &truth_N_other, &b_truth_N_other);
   fChain->SetBranchAddress("truth_N_pi0", &truth_N_pi0, &b_truth_N_pi0);
   fChain->SetBranchAddress("truth_N_piminus", &truth_N_piminus, &b_truth_N_piminus);
   fChain->SetBranchAddress("truth_N_piplus", &truth_N_piplus, &b_truth_N_piplus);
   fChain->SetBranchAddress("truth_N_proton", &truth_N_proton, &b_truth_N_proton);
   fChain->SetBranchAddress("truth_muon_charge", &truth_muon_charge, &b_truth_muon_charge);
   fChain->SetBranchAddress("truth_reco_muonCharge", &truth_reco_muonCharge, &b_truth_reco_muonCharge);
   fChain->SetBranchAddress("truth_target_material", &truth_target_material, &b_truth_target_material);
   fChain->SetBranchAddress("truth_vertex_module", &truth_vertex_module, &b_truth_vertex_module);
   fChain->SetBranchAddress("truth_vertex_plane", &truth_vertex_plane, &b_truth_vertex_plane);
   fChain->SetBranchAddress("truth_muon_E", &truth_muon_E, &b_truth_muon_E);
   fChain->SetBranchAddress("truth_muon_px", &truth_muon_px, &b_truth_muon_px);
   fChain->SetBranchAddress("truth_muon_py", &truth_muon_py, &b_truth_muon_py);
   fChain->SetBranchAddress("truth_muon_pz", &truth_muon_pz, &b_truth_muon_pz);
   fChain->SetBranchAddress("truth_muon_theta_wrtbeam", &truth_muon_theta_wrtbeam, &b_truth_muon_theta_wrtbeam);
   fChain->SetBranchAddress("truth_gamma_parentID", truth_gamma_parentID, &b_truth_gamma_parentID);
   fChain->SetBranchAddress("truth_gamma_trackID", truth_gamma_trackID, &b_truth_gamma_trackID);
   fChain->SetBranchAddress("truth_pi0_parentID", truth_pi0_parentID, &b_truth_pi0_parentID);
   fChain->SetBranchAddress("truth_pi0_trackID", truth_pi0_trackID, &b_truth_pi0_trackID);
   fChain->SetBranchAddress("truth_proton_parentID", truth_proton_parentID, &b_truth_proton_parentID);
   fChain->SetBranchAddress("truth_proton_trackID", truth_proton_trackID, &b_truth_proton_trackID);
   fChain->SetBranchAddress("truth_gamma_E", truth_gamma_E, &b_truth_gamma_E);
   fChain->SetBranchAddress("truth_gamma_px", truth_gamma_px, &b_truth_gamma_px);
   fChain->SetBranchAddress("truth_gamma_py", truth_gamma_py, &b_truth_gamma_py);
   fChain->SetBranchAddress("truth_gamma_pz", truth_gamma_pz, &b_truth_gamma_pz);
   fChain->SetBranchAddress("truth_gamma_theta_wrtbeam", truth_gamma_theta_wrtbeam, &b_truth_gamma_theta_wrtbeam);
   fChain->SetBranchAddress("truth_gamma_vtx_x", truth_gamma_vtx_x, &b_truth_gamma_vtx_x);
   fChain->SetBranchAddress("truth_gamma_vtx_y", truth_gamma_vtx_y, &b_truth_gamma_vtx_y);
   fChain->SetBranchAddress("truth_gamma_vtx_z", truth_gamma_vtx_z, &b_truth_gamma_vtx_z);
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
   fChain->SetBranchAddress("truth_pi0_E", truth_pi0_E, &b_truth_pi0_E);
   fChain->SetBranchAddress("truth_pi0_px", truth_pi0_px, &b_truth_pi0_px);
   fChain->SetBranchAddress("truth_pi0_py", truth_pi0_py, &b_truth_pi0_py);
   fChain->SetBranchAddress("truth_pi0_pz", truth_pi0_pz, &b_truth_pi0_pz);
   fChain->SetBranchAddress("truth_pi0_theta_wrtbeam", truth_pi0_theta_wrtbeam, &b_truth_pi0_theta_wrtbeam);
   fChain->SetBranchAddress("truth_pi0_vtx_x", truth_pi0_vtx_x, &b_truth_pi0_vtx_x);
   fChain->SetBranchAddress("truth_pi0_vtx_y", truth_pi0_vtx_y, &b_truth_pi0_vtx_y);
   fChain->SetBranchAddress("truth_pi0_vtx_z", truth_pi0_vtx_z, &b_truth_pi0_vtx_z);
   fChain->SetBranchAddress("truth_proton_E", truth_proton_E, &b_truth_proton_E);
   fChain->SetBranchAddress("truth_proton_px", truth_proton_px, &b_truth_proton_px);
   fChain->SetBranchAddress("truth_proton_py", truth_proton_py, &b_truth_proton_py);
   fChain->SetBranchAddress("truth_proton_pz", truth_proton_pz, &b_truth_proton_pz);
   fChain->SetBranchAddress("truth_proton_theta_wrtbeam", truth_proton_theta_wrtbeam, &b_truth_proton_theta_wrtbeam);
   fChain->SetBranchAddress("truth_proton_vtx_x", truth_proton_vtx_x, &b_truth_proton_vtx_x);
   fChain->SetBranchAddress("truth_proton_vtx_y", truth_proton_vtx_y, &b_truth_proton_vtx_y);
   fChain->SetBranchAddress("truth_proton_vtx_z", truth_proton_vtx_z, &b_truth_proton_vtx_z);
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
   fChain->SetBranchAddress("CCProtonPi0_isProtonInsideOD", CCProtonPi0_isProtonInsideOD, &b_CCProtonPi0_isProtonInsideOD);
   fChain->SetBranchAddress("CCProtonPi0_ntrajProtonProng", CCProtonPi0_ntrajProtonProng, &b_CCProtonPi0_ntrajProtonProng);
   fChain->SetBranchAddress("CCProtonPi0_proton_kinked", CCProtonPi0_proton_kinked, &b_CCProtonPi0_proton_kinked);
   fChain->SetBranchAddress("CCProtonPi0_proton_odMatch", CCProtonPi0_proton_odMatch, &b_CCProtonPi0_proton_odMatch);
   fChain->SetBranchAddress("CCProtonPi0_proton_trk_pat_history", CCProtonPi0_proton_trk_pat_history, &b_CCProtonPi0_proton_trk_pat_history);
   fChain->SetBranchAddress("CCProtonPi0_trajProtonProngPDG", CCProtonPi0_trajProtonProngPDG, &b_CCProtonPi0_trajProtonProngPDG);
   fChain->SetBranchAddress("CCProtonPi0_trajProtonProngPrimary", CCProtonPi0_trajProtonProngPrimary, &b_CCProtonPi0_trajProtonProngPrimary);
   fChain->SetBranchAddress("CCProtonPi0_endProtonTrajMomentum", CCProtonPi0_endProtonTrajMomentum, &b_CCProtonPi0_endProtonTrajMomentum);
   fChain->SetBranchAddress("CCProtonPi0_endProtonTrajXPosition", CCProtonPi0_endProtonTrajXPosition, &b_CCProtonPi0_endProtonTrajXPosition);
   fChain->SetBranchAddress("CCProtonPi0_endProtonTrajYPosition", CCProtonPi0_endProtonTrajYPosition, &b_CCProtonPi0_endProtonTrajYPosition);
   fChain->SetBranchAddress("CCProtonPi0_endProtonTrajZPosition", CCProtonPi0_endProtonTrajZPosition, &b_CCProtonPi0_endProtonTrajZPosition);
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
   fChain->SetBranchAddress("CCProtonPi0_proton_score", CCProtonPi0_proton_score, &b_CCProtonPi0_proton_score);
   fChain->SetBranchAddress("CCProtonPi0_proton_score1", CCProtonPi0_proton_score1, &b_CCProtonPi0_proton_score1);
   fChain->SetBranchAddress("CCProtonPi0_proton_score2", CCProtonPi0_proton_score2, &b_CCProtonPi0_proton_score2);
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
   Notify();
}


CCProtonPi0::~CCProtonPi0()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t CCProtonPi0::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t CCProtonPi0::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

Bool_t CCProtonPi0::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

    cout<<"Done!"<<endl;
   return kTRUE;
}

void CCProtonPi0::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

Int_t CCProtonPi0::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void CCProtonPi0::initHistograms()
{
    
    pID_purity = new TH1F( "pID_purity","Proton Purity",binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max() );
    pID_purity->GetXaxis()->SetTitle("Proton Purity = Captured Proton / Captured Total Events");
    pID_purity->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList.particleScore.get_width()));
    
    pID_efficiency = new TH1F( "pID_efficiency","Proton Efficiency",binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max() );
    pID_efficiency->GetXaxis()->SetTitle("Proton Efficiency = Captured Proton / Total Protons");
    pID_efficiency->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList.particleScore.get_width()));
    
    pID_piplus = new TH1F( "pID_piplus","Pi Plus",binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max() );
    pID_piplus->GetXaxis()->SetTitle("Pi Plus");
    pID_piplus->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList.particleScore.get_width()));
    
    pID_piminus = new TH1F( "pID_piminus","Pi Minus",binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max() );
    pID_piminus->GetXaxis()->SetTitle("Pi Minus");
    pID_piminus->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList.particleScore.get_width()));
    
    pID_proton = new TH1F( "pID_proton","Proton",binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max() );
    pID_proton->GetXaxis()->SetTitle("Proton");
    pID_proton->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList.particleScore.get_width()));
    
    pID_other = new TH1F( "pID_other","Other",binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max() );
    pID_other->GetXaxis()->SetTitle("Other");
    pID_other->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList.particleScore.get_width()));
    
    beamEnergy_mc = new TH1F( "beamEnergy_mc","True Beam Energy",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    beamEnergy_mc->GetXaxis()->SetTitle("True Beam Energy MeV");
    beamEnergy_mc->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList.beamE.get_width()));
    
    beamEnergy_reco = new TH1F( "beamEnergy_reco","Reconstructed Beam Energy",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    beamEnergy_reco->GetXaxis()->SetTitle("Reconstructed Beam Energy MeV");
    beamEnergy_reco->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList.beamE.get_width()));
    
    beamEnergy_error = new TH1F( "beamEnergy_error","Error on Beam Energy",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    beamEnergy_error->GetXaxis()->SetTitle("(Reco- True) / True");
    beamEnergy_error->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList.error.get_width()));
    
    beamEnergy_reco_mc = new TH2F( "beamEnergy_reco_mc","True vs Reconstructed Beam Energy",
                                binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max(),
                                binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max());
    beamEnergy_reco_mc->GetXaxis()->SetTitle("Reconstructed Beam Energy MeV");
    beamEnergy_reco_mc->GetYaxis()->SetTitle("True Beam Energy MeV");
    
    q2_mc = new TH1F( "q2_mc","True Q^{2}",binList.q2.get_nBins(), binList.q2.get_min(), binList.q2.get_max() );
    q2_mc->GetXaxis()->SetTitle("True Q^{2} MeV");
    q2_mc->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList.q2.get_width()));
    
    q2_reco = new TH1F( "q2_reco","Reconstructed Q^{2}",binList.q2.get_nBins(), binList.q2.get_min(), binList.q2.get_max() );
    q2_reco->GetXaxis()->SetTitle("Reconstructed Q^{2} MeV");
    q2_reco->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList.q2.get_width()));
    
    q2_error = new TH1F( "q2_error","Error on Q^{2}",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    q2_error->GetXaxis()->SetTitle("(Reco- True) / True");
    q2_error->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList.error.get_width()));
    
    q2_reco_mc = new TH2F( "q2_reco_mc","True vs Reconstructed Q^{2}",
                                binList.q2.get_nBins(), binList.q2.get_min(), binList.q2.get_max(),
                                binList.q2.get_nBins(), binList.q2.get_min(), binList.q2.get_max());
    q2_reco_mc->GetXaxis()->SetTitle("Reconstructed Q^{2} MeV");
    q2_reco_mc->GetYaxis()->SetTitle("True Q^{2} MeV");
    
    int_channel = new TH1F( "int_channel","Interaction Channel",binList.int_channel.get_nBins(), binList.int_channel.get_min(), binList.int_channel.get_max() );
    int_channel->GetXaxis()->SetTitle("1 = QE, 2 = Resonant, 3 = DIS, 4 = Coh pi");
    int_channel->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList.int_channel.get_width()));
    
    vertex_x_y_true = new TH2F( "vertex_x_y_true","True Vertex X vs Y",   binList.vertex_x_y.get_nBins(), binList.vertex_x_y.get_min(), binList.vertex_x_y.get_max(),
                                                                binList.vertex_x_y.get_nBins(), binList.vertex_x_y.get_min(), binList.vertex_x_y.get_max());
    vertex_x_y_true->GetXaxis()->SetTitle("True Vertex X [mm]");
    vertex_x_y_true->GetYaxis()->SetTitle("True Vertex Y [mm]");

    vertex_x_y_reco = new TH2F( "vertex_x_y_reco","Reconstructed Vertex X vs Y",   binList.vertex_x_y.get_nBins(), binList.vertex_x_y.get_min(), binList.vertex_x_y.get_max(),
                                                                binList.vertex_x_y.get_nBins(), binList.vertex_x_y.get_min(), binList.vertex_x_y.get_max());
    vertex_x_y_reco->GetXaxis()->SetTitle("Reconstructed Vertex X [mm]");
    vertex_x_y_reco->GetYaxis()->SetTitle("Reconstructed Vertex Y [mm]");
    
    vertex_z_true = new TH1F( "vertex_z_true","True Vertex Z",binList.vertex_z.get_nBins(), binList.vertex_z.get_min(), binList.vertex_z.get_max() );
    vertex_z_true->GetXaxis()->SetTitle("z = 4293 Target, #bf{z = 5810 Interaction Region}, z = 8614 ECAL, z = 9088 HCAL");
    vertex_z_true->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList.vertex_z.get_width()));
    
    vertex_z_reco = new TH1F( "vertex_z_reco","Reconstructed Vertex Z",binList.vertex_z.get_nBins(), binList.vertex_z.get_min(), binList.vertex_z.get_max() );
    vertex_z_reco->GetXaxis()->SetTitle("z = 4293 Target, #bf{z = 5810 Interaction Region}, z = 8614 ECAL, z = 9088 HCAL");
    vertex_z_reco->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList.vertex_z.get_width()));
    
    vertex_z_error = new TH1F( "vertex_z_error","Error on Vertex Z",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    vertex_z_error->GetXaxis()->SetTitle("(Reco- True) / True");
    vertex_z_error->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList.error.get_width()));
    
    mgg_reco = new TH1F( "mgg_reco","Reconstructed Invariant Mass",binList.mgg_reco.get_nBins(), binList.mgg_reco.get_min(), binList.mgg_reco.get_max() );
    mgg_reco ->GetXaxis()->SetTitle("Reconstructed Invariant Mass");
    mgg_reco ->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList.mgg_reco.get_width()));
    
    vertex_z_reco_mc = new TH2F( "vertex_z_reco_mc","True vs Reconstructed Vertex Z",
                                binList.vertex_z.get_nBins(), binList.vertex_z.get_min(), binList.vertex_z.get_max(),
                                binList.vertex_z.get_nBins(), binList.vertex_z.get_min(), binList.vertex_z.get_max());
    vertex_z_reco_mc->GetXaxis()->SetTitle("Reconstructed Vertex Z");
    vertex_z_reco_mc->GetYaxis()->SetTitle("True Vertex Z");
    
}

#endif //CCProtonPi0_cpp
