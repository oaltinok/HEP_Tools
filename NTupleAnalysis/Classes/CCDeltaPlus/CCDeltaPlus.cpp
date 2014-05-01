/*
    See CCDeltaPlus.h header for Class Information
*/

#define CCDeltaPlus_cxx

#include "CCDeltaPlus.h"
#include "initHistograms.cpp"
#include "Default_Functions.cpp"
#include "Muon_Functions.cpp"
#include "Proton_Functions.cpp"
#include "Pion_Functions.cpp"

void CCDeltaPlus::specifyRunTime()
{
    isDataAnalysis = true;
    isMC = true;
    applyProtonScore = false;
    minProtonScore = 0.35;
    is_pID_Studies = false;
}

void CCDeltaPlus::run(string playlist)
{
    //------------------------------------------------------------------------
    // Open files for writing
    //------------------------------------------------------------------------
    openFiles();
   
    //------------------------------------------------------------------------
    // Create chain
    //------------------------------------------------------------------------

    TChain* fChain = new TChain("CCDeltaPlusAna");
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
    // Branch Selection for Performance
    //------------------------------------------------------------------------
//     fChain->SetBranchStatus("*",0);  // disable all branches

//     fChain->SetBranchStatus("ev_run",1);  // activate

    // Cut Statistics
    double nAll = 0;
    double nCut_Vertex_None = 0;
    double nCut_Vertex_Null = 0;
    double nCut_Vertex_Not_Analyzable = 0; 
    double nCut_Vertex_Not_Fiducial = 0;    
    double nCut_Vertex_Michel_Exist = 0;           
    double nCut_Muon_None = 0;              
    double nCut_Muon_Score_Low = 0;
    double nCut_EndPoint_Michel_Exist = 0;
    double nCut_secEndPoint_Michel_Exist = 0;
    double nCut_Proton_None = 0;            
    double nCut_Proton_Score = 0;
    double nCut_True_Muon = 0;
    double nCut_True_Proton = 0;
    double nCut_Reco_Muon_NoProblem = 0;
    
    double nTrue_Muon_None = 0;             double nReco_Muon_None = 0;
    double nTrue_Michel_Exist = 0;          double nReco_Michel_Exist = 0;
    double nTrue_Proton_None = 0;           double nReco_Proton_None = 0;
    
    int tempInd;
    int indRecoProton;
    int nRecoProtons = 0;
    int nTrueProtons = 0;
    

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
        
        if( Cut_Vertex_None == 1) continue;
        nCut_Vertex_None++;
        
        if( Cut_Vertex_Null == 1) continue;
        nCut_Vertex_Null++;
        
        if( Cut_Vertex_Not_Analyzable == 1) continue;
        nCut_Vertex_Not_Analyzable++;
        
        if( Cut_Vertex_Not_Fiducial == 1) continue;
        nCut_Vertex_Not_Fiducial++;
        
        if( Cut_Muon_None == 1) continue;
        nCut_Muon_None++;
        
        if( Cut_Muon_Score_Low == 1) continue;
        nCut_Muon_Score_Low++;
           
        if( Cut_Vertex_Michel_Exist == 1) continue;
        nCut_Vertex_Michel_Exist++;

        if( Cut_EndPoint_Michel_Exist == 1) continue;
        nCut_EndPoint_Michel_Exist++;
 
        if( Cut_secEndPoint_Michel_Exist == 1) continue;
        nCut_secEndPoint_Michel_Exist++;
        
        if( Cut_Proton_None == 1) continue;
        nCut_Proton_None++;
        
        // Find Best Proton in Reco
        indRecoProton = findBestProton();
        
        if ( applyProtonScore && CCDeltaPlusAna_proton_score[indRecoProton] < minProtonScore ) continue;
        nCut_Proton_Score++;
        
        //------------------------------------------------------------------
        // Get True Particle Indices and Some Sanity Checks
        //------------------------------------------------------------------
        muon.ind = findTrueParticle(PDG_List::mu_minus);
        if(muon.ind == -1) continue;
        nCut_True_Muon++;
        
        proton.ind = findTrueParticle(PDG_List::proton);
        if(proton.ind == -1) continue;
        nCut_True_Proton++;
        
        if(CCDeltaPlusAna_muon_pz == 0) continue;
        nCut_Reco_Muon_NoProblem++;
      
        //------------------------------------------------------------------
        // pID Studies
        //------------------------------------------------------------------
        if( is_pID_Studies){
            cout<<"Collecting Data for pID Studies"<<endl;
            for(int i = 0; i < 10; i++){
                if(CCDeltaPlusAna_proton_score[i] == -1) break;
                if(CCDeltaPlusAna_trajProtonProngPDG[i] == 2212){
                    pID_proton->Fill(CCDeltaPlusAna_proton_score[i]);
                }else if(CCDeltaPlusAna_trajProtonProngPDG[i] == 211){
                    pID_piplus->Fill(CCDeltaPlusAna_proton_score[i]);
                }else if(CCDeltaPlusAna_trajProtonProngPDG[i] == -211){
                    pID_piminus->Fill(CCDeltaPlusAna_proton_score[i]);
                }else{
                    pID_other->Fill(CCDeltaPlusAna_proton_score[i]);
                }
            }
        }
        
        
        if ( isDataAnalysis){
            //------------------------------------------------------------------
            // Fill Particles
            //------------------------------------------------------------------
            if( isMC ){
                fillMuonTrue();
                fillProtonTrue();
//                 fillPionTrue();
            }
        
            double truthMatchMomentum = CCDeltaPlusAna_trajProtonProngMomentum[indRecoProton];
            proton.momentum[1] = truthMatchMomentum;
            
            // Fill Reconstructed Information
            fillMuonReco();
            fillProtonReco(indRecoProton);
//             fillPionReco();
            
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
    
    cout<<"Done!"<<endl;
    
    cout<<">> Writing "<<cutFile<<endl;
    cutText<<"nAll                          "<<nAll<<endl;
    cutText<<"Cut_Vertex_None               "<<nCut_Vertex_None<<endl;
    cutText<<"Cut_Vertex_Null               "<<nCut_Vertex_Null<<endl;
    cutText<<"Cut_Vertex_Not_Analyzable     "<<nCut_Vertex_Not_Analyzable<<endl;
    cutText<<"Cut_Vertex_Not_Fiducial       "<<nCut_Vertex_Not_Fiducial<<endl;
    cutText<<"Cut_Muon_None                 "<<nCut_Muon_None<<endl;
    cutText<<"Cut_Muon_Score_Low            "<<nCut_Muon_Score_Low<<endl; 
    cutText<<"Cut_Vertex_Michel_Exist       "<<nCut_Vertex_Michel_Exist<<endl;
    cutText<<"Cut_EndPoint_Michel_Exist     "<<nCut_EndPoint_Michel_Exist<<endl;
    cutText<<"Cut_secEndPoint_Michel_Exist  "<<nCut_secEndPoint_Michel_Exist<<endl;
    cutText<<"Cut_Proton_None               "<<nCut_Proton_None<<endl;
    cutText<<"Cut_Proton_Score              "<<nCut_Proton_Score<<endl;
    cutText<<"Cut_True_Muon                 "<<nCut_True_Muon<<endl;
    cutText<<"Cut_True_Proton               "<<nCut_True_Proton<<endl;
    cutText<<"Cut_Reco_Muon_NoProblem       "<<nCut_Reco_Muon_NoProblem<<endl;
    
    
    // Write the Root Files
    write_RootFile();           //CCDeltaPlus
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

int CCDeltaPlus::get_pID_Stats()
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
        nCapturedEvents =   nCapturedEvents+
                            pID_proton->GetBinContent(i) +
                            pID_piplus->GetBinContent(i) +
                            pID_piminus->GetBinContent(i) +
                            pID_other->GetBinContent(i);
                            
        purity = nProton / nCapturedEvents;
        efficiency = nProton / nTotalProton;
        cout<<"pID = "<<0.025*i<<" Purity = "<<purity<<" Efficiency "<<efficiency<<endl;
    }
}

void CCDeltaPlus::fillCCDeltaPlus()
{
    beamEnergy_mc->Fill(mc_incomingE);
//     beamEnergy_reco->Fill(Erec);
//     beamEnergy_error->Fill( (mc_incomingE - Erec) / mc_incomingE );
//     beamEnergy_reco_mc->Fill(Erec,mc_incomingE);
    
    q2_mc->Fill(mc_Q2 / mevSq_to_gevSq);
//     q2_reco->Fill(Q2/ mevSq_to_gevSq);
//     q2_error->Fill( (mc_Q2 - Q2) / mc_Q2 );
//     q2_reco_mc->Fill(Q2/mevSq_to_gevSq,mc_Q2 /mevSq_to_gevSq);

    vertex_z_true->Fill(mc_vtx[2]);
    vertex_z_reco->Fill(CCDeltaPlusAna_vtx[2]);
    vertex_z_error->Fill((mc_vtx[2] - CCDeltaPlusAna_vtx[2])/mc_vtx[2]);
    vertex_z_reco_mc->Fill(CCDeltaPlusAna_vtx[2],mc_vtx[2]);
    
    vertex_x_y_true->Fill(mc_vtx[0],mc_vtx[1]);
    vertex_x_y_reco->Fill(CCDeltaPlusAna_vtx[0],CCDeltaPlusAna_vtx[1]);
    
    int_channel->Fill(mc_intType);

    n_FSParticles->Fill(mc_nFSPart);
//     n_gammas->Fill();

}

void CCDeltaPlus::initVariables()
{
    cout<<"Initializing CCDeltaPlus Class"<<endl;
    
    // File Locations
    rootDir =   Folder_List::f_Root_CCDeltaPlus;
    plotDir =   Folder_List::f_Plot_CCDeltaPlus;
    
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


void CCDeltaPlus::fillHistograms()
{
    fillCCDeltaPlus();
    muon.fill_Histograms();
    pion.fill_Histograms();
    proton.fill_Histograms();

}

void CCDeltaPlus::write_RootFile()
{
    cout<<">> Writing "<<rootDir<<endl;
    f->Write();
}

void CCDeltaPlus::closeFiles()
{
    readme.close();
}

void CCDeltaPlus::openFiles()
{

    // Open Readme File
    readmeFile = getFileLocation(Folder_List::OUT, Folder_List::TEXT, Folder_List::F_TEXT_README);
    readme.open( readmeFile.c_str() );
    
    if( !readme.is_open() ){
        cerr<<"Cannot open output text file: "<<readmeFile<<endl;
        exit(1);
    }
    
    // Open Cut File
    cutFile = getFileLocation(Folder_List::OUT, Folder_List::TEXT, Folder_List::F_TEXT_CUT);
    cutFile = cutFile + "_" + channelTag + ".txt";
    cutText.open( cutFile.c_str() );
    
    if( !cutText.is_open() ){
        cerr<<"Cannot open output text file: "<<cutFile<<endl;
        exit(1);
    }

    writeReadme();
}

void CCDeltaPlus::writeReadme()
{
    readme<<"Test"<<endl;
}

/*
--------------------------------------------------------------------------------
 Beam Energy Cu: isBeamEnergyLow(double maxEnergy)
    Incoming Beam Energy must be lower than a maximum Energy
--------------------------------------------------------------------------------
*/
bool CCDeltaPlus::isBeamEnergyLow(double maxEnergy)
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
int CCDeltaPlus::findTrueParticle(int targetPDG)
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
int CCDeltaPlus::countParticles(int targetPDG, bool applyPCut)
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

