/*
    See CCDeltaPlus.h header for Class Information
*/

#define CCDeltaPlus_cxx

#include "CCDeltaPlus.h"
#include "Cuts.cpp"



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

    initVariables();
    initHistograms();
    

    //------------------------------------------------------------------------
    // Branch Selection for Performance
    //------------------------------------------------------------------------
//     fChain->SetBranchStatus("*",0);  // disable all branches
// 
//     // Analysis Variables
//     fChain->SetBranchStatus("ev_run",1);  // activate
//     fChain->SetBranchStatus("ev_subrun",1);  // activate
//     fChain->SetBranchStatus("ev_gate",1);  // activate
//     fChain->SetBranchStatus("truth_reco_minos_match",1);  // activate
//     fChain->SetBranchStatus("mc_intType",1);  // activate
//     
//     // Q2
//     fChain->SetBranchStatus("mc_Q2",1);  // activate
//     fChain->SetBranchStatus("Q2",1);  // activate
//     
//     // Cut Variables
//     fChain->SetBranchStatus("mc_vtx",1);
//     
//     // Incoming Particle
//     fChain->SetBranchStatus("mc_incomingPartVec",1);  // activate
//     fChain->SetBranchStatus("mc_incomingE",1);  // activate
//    
//    
//     // Final State Particles
//     fChain->SetBranchStatus("mc_FSPartPx",1);  // activate
//     fChain->SetBranchStatus("mc_FSPartPy",1);  // activate
//     fChain->SetBranchStatus("mc_FSPartPz",1);  // activate
//     fChain->SetBranchStatus("mc_FSPartE",1);  // activate
//     fChain->SetBranchStatus("mc_nFSPart",1);  // activate
//     fChain->SetBranchStatus("mc_FSPartPDG",1);
//     
//     fChain->SetBranchStatus("pimom",1);
//     fChain->SetBranchStatus("pienergy",1);  // activate
//     
//     fChain->SetBranchStatus("mumom",1);  // activate
//     
//     fChain->SetBranchStatus("tfiducial",1);  // activate
//     fChain->SetBranchStatus("survive_minos_match",1);  // activate\
//     fChain->SetBranchStatus("survive_fiducial",1);  // activate
//     fChain->SetBranchStatus("Erec",1);  // activate
//     fChain->SetBranchStatus("mc_intType",1);  // activate

    

    // Cut Statistics
    double nAll = 0;
    double nVolume = 0;
    double nBeamEnergy = 0;
    double nMuon = 0;
    double nMinos = 0;
    double nProton = 0;
    double nPion = 0;
    double nFSPart = 0;
    double nBeamEnergyFail = 0;
    
    

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
    
//         cout<<ev_run<<" | "<<ev_subrun<<" | "<<ev_gate<<endl;
//         cout<<CCDeltaPlusAna_tammy_proton_score<<" | "<<CCDeltaPlusAna_tammy_proton_4p[3]<<endl;


        

        
//         if (jentry >= 20000){
//             break;
//         }
    
        
//         ----------------------------------------------------------------------
//         Apply Cuts
//         ----------------------------------------------------------------------
        
        // Count All Events before Cuts
        nAll++;
        
//         if( mc_nFSPart > max_nFSPart){
//             continue;
//         }
//         nFSPart++;
        
        // Volume Cut
        if( !isVertexContained()){
            continue;
        }
        nVolume++;

        // Incoming Beam Energy Cut
        if ( !isBeamEnergyLow(maxBeamEnergy) ){
            continue;
        }
        nBeamEnergy++;
        
        // Muon Cut
        muon.ind = findParticle(PDG_List::mu_minus);
        if(muon.ind == -1){
            continue;
        }
        nMuon++;
        

        // Minos Match
        if ( CCDeltaPlusAna_tammy_muon_minos_track == 0 ){
            continue;
        }
        nMinos++;
        
                    cout<<CCDeltaPlusAna_tammy_muon_minos_track<<" | "<<CCDeltaPlusAna_leptonE[3]<<" | "<<mc_FSPartE[0]<<endl;
        
        // Proton Cut
        proton.ind = findProton();
        if(proton.ind == -1){
            continue;
        }
        if(isProtonShort(proton.ind)){
            continue;
        }
        nProton++;
        
//         Pion Cut
        pion.ind = findPion();
        if(pion.ind == -1){
            continue;
        }
        nPion++;
        
//         if(Erec == -1){
//             continue;
//         }
//         nBeamEnergyFail++;
        
        

        if ( isDataAnalysis){
//         ----------------------------------------------------------------------
//         Fill Particles
//         ----------------------------------------------------------------------
        
            if( isMC ){
                fillParticleTrue(muon);
                fillParticleTrue(proton);
//                 fillParticleTrue(pion);
            }
        
            // Fill Reconstructed Information
            fillMuon();
            fillProton();
//             fillPion();
            
            muon.set_errors();
            proton.set_errors();
//             pion.set_errors();
            
//         ----------------------------------------------------------------------
//         Fill Histograms
//         ----------------------------------------------------------------------
       
            fillHistograms();            
        }
        
        


    }    
    cout<<"Done!"<<endl;
    
  
    nCutList->nAll->setValue(nAll);
    nCutList->nFSPart->setValue(nFSPart);
    nCutList->nVolume->setValue(nVolume);
    nCutList->nBeamEnergy->setValue(nBeamEnergy);
    nCutList->nMuon->setValue(nMuon);
    nCutList->nMinos->setValue(nMinos);
    nCutList->nProton->setValue(nProton);
    nCutList->nPion->setValue(nPion);
    nCutList->nBeamEnergyFail->setValue(nBeamEnergyFail);
    
    // Write the Root Files
    write_RootFile();           //CCDeltaPlus
    muon.write_RootFile();
    proton.write_RootFile();
    pion.write_RootFile();
    
    

    
    nCutList->writeCutTable();
    
    closeFiles();
    
    
    // Delete Dynamic Variables - I will modify the destructor for
    // better Memory Management
    delete nCutList;

}


void CCDeltaPlus::fillCCDeltaPlus()
{


    beamEnergy_mc->Fill(mc_incomingE);
    beamEnergy_reco->Fill(Erec);
    beamEnergy_error->Fill( (mc_incomingE - Erec) / mc_incomingE );
    beamEnergy_reco_mc->Fill(Erec,mc_incomingE);
    
    q2_mc->Fill(mc_Q2 / mevSq_to_gevSq);
    q2_reco->Fill(Q2/ mevSq_to_gevSq);
    q2_error->Fill( (mc_Q2 - Q2) / mc_Q2 );
    q2_reco_mc->Fill(Q2/mevSq_to_gevSq,mc_Q2 /mevSq_to_gevSq);
    
    int_channel->Fill(mc_intType);
    vertex_z->Fill(mc_vtx[2]);
    vertex_x_y->Fill(mc_vtx[0],mc_vtx[1]);
    n_FSParticles->Fill(mc_nFSPart);
//     n_gammas->Fill();


}

// -------------------------------------------------------------------------
//     Specific Functions
//--------------------------------------------------------------------------


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

    isDataAnalysis = true;
    isMC = true;

    // -------------------------------------------------------------------------
    //     Memory Allocation
    //--------------------------------------------------------------------------
    // Allocate Memory
    
    nCutList = new CutNumberList;
    
    nCutList->printList();
    
    // -------------------------------------------------------------------------
    //     Initialization
    //--------------------------------------------------------------------------
    // Default Beam Configuration
    beam_p3.SetXYZ(1.0,1.0,1.0);
    beam_p3.SetPhi(-1.554);
    beam_p3.SetTheta(0.059);
    
    max_nFSPart = 3;
    maxBeamEnergy = 20000; //MeV
    

    cout<<"Done!"<<endl;

}


void CCDeltaPlus::fillMuon()
{

    // Fill 4-Momentum
    muon.set_p4(    CCDeltaPlusAna_leptonE[0],
                    CCDeltaPlusAna_leptonE[1],
                    CCDeltaPlusAna_leptonE[2],
                    CCDeltaPlusAna_leptonE[3],
                    false);
    
    // set Angle wrt Beam
    muon.set_angleBeam(beam_p3, false);
    
    // set Angle wrt Muon
    muon.set_angleMuon(muon, false);
    

}

void CCDeltaPlus::fillProton()
{
    // Fill 4-Momentum
    proton.set_p4(  CCDeltaPlusAna_tammy_proton_4p[0],
                    CCDeltaPlusAna_tammy_proton_4p[1],
                    CCDeltaPlusAna_tammy_proton_4p[2],
                    CCDeltaPlusAna_tammy_proton_4p[3],
                    false);
    
    // set Angle wrt Beam
    proton.set_angleBeam(beam_p3, false);
    
    // set Angle wrt Muon
    proton.set_angleMuon(muon, false);


}

void CCDeltaPlus::fillPion()
{
    // Fill 4-Momentum
    pion.set_p4(    pimom[0],
                    pimom[1],
                    pimom[2],
                    pienergy,
                    false);
    
    // set Angle wrt Beam
    pion.set_angleBeam(beam_p3, false);
    
    // set Angle wrt Muon
    pion.set_angleMuon(muon, false);
}

void CCDeltaPlus::fillHistograms()
{
//     fillCCDeltaPlus();
    muon.fill_Histograms();
    pion.fill_Histograms();
    proton.fill_Histograms();

}

void CCDeltaPlus::write_RootFile()
{
    cout<<">> Writing "<<rootDir<<endl;
    f->Write();
}

void CCDeltaPlus::fillParticleTrue(Particle& part)
{
    int ind = part.ind;
    
    // Fill 4-Momentum
    part.set_p4(    mc_FSPartPx[ind],
                    mc_FSPartPy[ind],
                    mc_FSPartPz[ind],
                    mc_FSPartE[ind], 
                    true);
       
    // set Angle wrt Beam
    part.set_angleBeam(beam_p3, true);
    
    // set Angle wrt Muon
    part.set_angleMuon(muon, true);
    
    
}


void CCDeltaPlus::initHistograms()
{
    cout<<"Initializing Histograms"<<endl;
    
    beamEnergy_mc = new TH1F( "beamEnergy_mc","True Beam Energy",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    beamEnergy_mc->GetXaxis()->SetTitle("True Beam Energy MeV");
    beamEnergy_mc->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList.beamE.get_width()));
    
    beamEnergy_reco = new TH1F( "beamEnergy_reco","Reconstructed Beam Energy",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    beamEnergy_reco->GetXaxis()->SetTitle("Reconstructed Beam Energy MeV");
    beamEnergy_reco->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList.beamE.get_width()));
    
    beamEnergy_error = new TH1F( "beamEnergy_error","Error on Beam Energy",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    beamEnergy_error->GetXaxis()->SetTitle("(True - Reco) / True");
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
    q2_error->GetXaxis()->SetTitle("(True - Reco) / True");
    q2_error->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList.error.get_width()));
    
    q2_reco_mc = new TH2F( "q2_reco_mc","True vs Reconstructed Q^{2}",
                                binList.q2.get_nBins(), binList.q2.get_min(), binList.q2.get_max(),
                                binList.q2.get_nBins(), binList.q2.get_min(), binList.q2.get_max());
    q2_reco_mc->GetXaxis()->SetTitle("Reconstructed Q^{2} MeV");
    q2_reco_mc->GetYaxis()->SetTitle("True Q^{2} MeV");
    
    int_channel = new TH1F( "int_channel","Interaction Channel",binList.int_channel.get_nBins(), binList.int_channel.get_min(), binList.int_channel.get_max() );
    int_channel->GetXaxis()->SetTitle("1 = QE, 2 = Resonant, 3 = DIS, 4 = Coh pi");
    int_channel->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList.int_channel.get_width()));
    
    vertex_z = new TH1F( "vertex_z","True Vertex Z",binList.vertex_z.get_nBins(), binList.vertex_z.get_min(), binList.vertex_z.get_max() );
    vertex_z->GetXaxis()->SetTitle("z = 4293 Target, #bf{z = 5810 Interaction Region}, z = 8614 ECAL, z = 9088 HCAL");
    vertex_z->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList.vertex_z.get_width()));
    
    vertex_x_y = new TH2F( "vertex_x_y","True Vertex X vs Y",   binList.vertex_x_y.get_nBins(), binList.vertex_x_y.get_min(), binList.vertex_x_y.get_max(),
                                                                binList.vertex_x_y.get_nBins(), binList.vertex_x_y.get_min(), binList.vertex_x_y.get_max());
    vertex_x_y->GetXaxis()->SetTitle("True Vertex X [mm]");
    vertex_x_y->GetYaxis()->SetTitle("True Vertex Y [mm]");

    n_FSParticles = new TH1F( "n_FSParticles","Number of Final State Particles",binList.multiplicity.get_nBins(), binList.multiplicity.get_min(), binList.multiplicity.get_max() );
    n_FSParticles->GetXaxis()->SetTitle("Number of Final State Particles");
    n_FSParticles->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList.multiplicity.get_width()));
    
    n_gammas = new TH1F( "n_gammas","Number of Gammas",binList.multiplicity.get_nBins(), binList.multiplicity.get_min(), binList.multiplicity.get_max() );
    n_gammas->GetXaxis()->SetTitle("Number of Gammas");
    n_gammas->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList.multiplicity.get_width()));    
    

    
    cout<<"Done!"<<endl;
}




void CCDeltaPlus::closeFiles()
{
    readme.close();
}

void CCDeltaPlus::openFiles()
{
    readmeFile = getFileLocation(Folder_List::OUT, Folder_List::TEXT, Folder_List::F_TEXT_README);

    readme.open( readmeFile.c_str() );
    if( !readme.is_open() ){
        cerr<<"Cannot open output text file: "<<readmeFile<<endl;
        exit(1);
    }

    writeReadme();
}

void CCDeltaPlus::writeReadme()
{

readme<<"Test"<<endl;

}

// -------------------------------------------------------------------------
//     Default Functions
//--------------------------------------------------------------------------

#ifdef CCDeltaPlus_cxx
CCDeltaPlus::CCDeltaPlus()
{
    
}

CCDeltaPlus::~CCDeltaPlus()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t CCDeltaPlus::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t CCDeltaPlus::LoadTree(Long64_t entry)
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

void CCDeltaPlus::Init(string playlist, TChain* fChain)
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
   fChain->SetBranchAddress("is_GoodDirection1", &is_GoodDirection1, &b_is_GoodDirection1);
   fChain->SetBranchAddress("is_GoodPosition1", &is_GoodPosition1, &b_is_GoodPosition1);
   fChain->SetBranchAddress("is_GoodDirection2", &is_GoodDirection2, &b_is_GoodDirection2);
   fChain->SetBranchAddress("is_GoodPosition2", &is_GoodPosition2, &b_is_GoodPosition2);
   fChain->SetBranchAddress("is_GoodBlob1", &is_GoodBlob1, &b_is_GoodBlob1);
   fChain->SetBranchAddress("is_GoodBlob2", &is_GoodBlob2, &b_is_GoodBlob2);
   fChain->SetBranchAddress("is_anglescan", &is_anglescan, &b_is_anglescan);
   fChain->SetBranchAddress("is_anglescan_applied", &is_anglescan_applied, &b_is_anglescan_applied);
   fChain->SetBranchAddress("is_houghtransform", &is_houghtransform, &b_is_houghtransform);
   fChain->SetBranchAddress("is_houghtransform_applied", &is_houghtransform_applied, &b_is_houghtransform_applied);
   fChain->SetBranchAddress("is_twoDBlob", &is_twoDBlob, &b_is_twoDBlob);
   fChain->SetBranchAddress("anglescan_ncand", &anglescan_ncand, &b_anglescan_ncand);
   fChain->SetBranchAddress("anglescan_ncandx", &anglescan_ncandx, &b_anglescan_ncandx);
   fChain->SetBranchAddress("blob_ndof_1", &blob_ndof_1, &b_blob_ndof_1);
   fChain->SetBranchAddress("blob_ndof_2", &blob_ndof_2, &b_blob_ndof_2);
   fChain->SetBranchAddress("broken_track_most_us_plane", &broken_track_most_us_plane, &b_broken_track_most_us_plane);
   fChain->SetBranchAddress("discard_track_count", &discard_track_count, &b_discard_track_count);
   fChain->SetBranchAddress("dmode", &dmode, &b_dmode);
   fChain->SetBranchAddress("g1blob_ncluster", &g1blob_ncluster, &b_g1blob_ncluster);
   fChain->SetBranchAddress("g1blob_ndigit", &g1blob_ndigit, &b_g1blob_ndigit);
   fChain->SetBranchAddress("g1convidet", &g1convidet, &b_g1convidet);
   fChain->SetBranchAddress("g1mostevispdg", &g1mostevispdg, &b_g1mostevispdg);
   fChain->SetBranchAddress("g2blob_ncluster", &g2blob_ncluster, &b_g2blob_ncluster);
   fChain->SetBranchAddress("g2blob_ndigit", &g2blob_ndigit, &b_g2blob_ndigit);
   fChain->SetBranchAddress("g2convidet", &g2convidet, &b_g2convidet);
   fChain->SetBranchAddress("g2mostevispdg", &g2mostevispdg, &b_g2mostevispdg);
   fChain->SetBranchAddress("minos_trk_end_plane", &minos_trk_end_plane, &b_minos_trk_end_plane);
   fChain->SetBranchAddress("minos_trk_is_contained", &minos_trk_is_contained, &b_minos_trk_is_contained);
   fChain->SetBranchAddress("minos_trk_is_ok", &minos_trk_is_ok, &b_minos_trk_is_ok);
   fChain->SetBranchAddress("minos_trk_quality", &minos_trk_quality, &b_minos_trk_quality);
   fChain->SetBranchAddress("minos_trk_used_curvature", &minos_trk_used_curvature, &b_minos_trk_used_curvature);
   fChain->SetBranchAddress("minos_trk_used_range", &minos_trk_used_range, &b_minos_trk_used_range);
   fChain->SetBranchAddress("nblob_anglescan", &nblob_anglescan, &b_nblob_anglescan);
   fChain->SetBranchAddress("nblob_hough", &nblob_hough, &b_nblob_hough);
   fChain->SetBranchAddress("nmeson", &nmeson, &b_nmeson);
   fChain->SetBranchAddress("nmumcapture", &nmumcapture, &b_nmumcapture);
   fChain->SetBranchAddress("nmumdecay", &nmumdecay, &b_nmumdecay);
   fChain->SetBranchAddress("nmupdecay", &nmupdecay, &b_nmupdecay);
   fChain->SetBranchAddress("nn", &nn, &b_nn);
   fChain->SetBranchAddress("np", &np, &b_np);
   fChain->SetBranchAddress("npi0", &npi0, &b_npi0);
   fChain->SetBranchAddress("npi02", &npi02, &b_npi02);
   fChain->SetBranchAddress("npim", &npim, &b_npim);
   fChain->SetBranchAddress("npim2", &npim2, &b_npim2);
   fChain->SetBranchAddress("npimcapture", &npimcapture, &b_npimcapture);
   fChain->SetBranchAddress("npimdecay", &npimdecay, &b_npimdecay);
   fChain->SetBranchAddress("npiminelastic", &npiminelastic, &b_npiminelastic);
   fChain->SetBranchAddress("npip", &npip, &b_npip);
   fChain->SetBranchAddress("npip2", &npip2, &b_npip2);
   fChain->SetBranchAddress("npipcapture", &npipcapture, &b_npipcapture);
   fChain->SetBranchAddress("npipdecay", &npipdecay, &b_npipdecay);
   fChain->SetBranchAddress("npipinelastic", &npipinelastic, &b_npipinelastic);
   fChain->SetBranchAddress("npipm", &npipm, &b_npipm);
   fChain->SetBranchAddress("od_energeticTower", &od_energeticTower, &b_od_energeticTower);
   fChain->SetBranchAddress("phys_energy_in_road_downstream_nplanes", &phys_energy_in_road_downstream_nplanes, &b_phys_energy_in_road_downstream_nplanes);
   fChain->SetBranchAddress("phys_energy_in_road_upstream_nplanes", &phys_energy_in_road_upstream_nplanes, &b_phys_energy_in_road_upstream_nplanes);
   fChain->SetBranchAddress("phys_n_dead_discr_pair", &phys_n_dead_discr_pair, &b_phys_n_dead_discr_pair);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_in_prim_track_region", &phys_n_dead_discr_pair_in_prim_track_region, &b_phys_n_dead_discr_pair_in_prim_track_region);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_two_mod_downstream_prim_track", &phys_n_dead_discr_pair_two_mod_downstream_prim_track, &b_phys_n_dead_discr_pair_two_mod_downstream_prim_track);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_two_mod_upstream_prim_vtx", &phys_n_dead_discr_pair_two_mod_upstream_prim_vtx, &b_phys_n_dead_discr_pair_two_mod_upstream_prim_vtx);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_upstream_prim_track_proj", &phys_n_dead_discr_pair_upstream_prim_track_proj, &b_phys_n_dead_discr_pair_upstream_prim_track_proj);
   fChain->SetBranchAddress("phys_vertex_is_fiducial", &phys_vertex_is_fiducial, &b_phys_vertex_is_fiducial);
   fChain->SetBranchAddress("primary_index", &primary_index, &b_primary_index);
   fChain->SetBranchAddress("primary_multiplicity", &primary_multiplicity, &b_primary_multiplicity);
   fChain->SetBranchAddress("primary_multiplicity2", &primary_multiplicity2, &b_primary_multiplicity2);
   fChain->SetBranchAddress("survive_all", &survive_all, &b_survive_all);
   fChain->SetBranchAddress("survive_do_muon", &survive_do_muon, &b_survive_do_muon);
   fChain->SetBranchAddress("survive_do_vertex", &survive_do_vertex, &b_survive_do_vertex);
   fChain->SetBranchAddress("survive_fiducial", &survive_fiducial, &b_survive_fiducial);
   fChain->SetBranchAddress("survive_gammatrack", &survive_gammatrack, &b_survive_gammatrack);
   fChain->SetBranchAddress("survive_has_vertex", &survive_has_vertex, &b_survive_has_vertex);
   fChain->SetBranchAddress("survive_minos_match", &survive_minos_match, &b_survive_minos_match);
   fChain->SetBranchAddress("survive_onetrackpervtx", &survive_onetrackpervtx, &b_survive_onetrackpervtx);
   fChain->SetBranchAddress("survive_plausible", &survive_plausible, &b_survive_plausible);
   fChain->SetBranchAddress("survive_prefilter", &survive_prefilter, &b_survive_prefilter);
   fChain->SetBranchAddress("survive_three_vertex", &survive_three_vertex, &b_survive_three_vertex);
   fChain->SetBranchAddress("survive_vtx_blob", &survive_vtx_blob, &b_survive_vtx_blob);
   fChain->SetBranchAddress("tammy_CreatedShortTracks", &tammy_CreatedShortTracks, &b_tammy_CreatedShortTracks);
   fChain->SetBranchAddress("tammy_FailContainedProng", &tammy_FailContainedProng, &b_tammy_FailContainedProng);
   fChain->SetBranchAddress("tammy_FailExitingProng", &tammy_FailExitingProng, &b_tammy_FailExitingProng);
   fChain->SetBranchAddress("tammy_FailFidVolume", &tammy_FailFidVolume, &b_tammy_FailFidVolume);
   fChain->SetBranchAddress("tammy_FailOutTracks", &tammy_FailOutTracks, &b_tammy_FailOutTracks);
   fChain->SetBranchAddress("tammy_FailRefitFidVolume", &tammy_FailRefitFidVolume, &b_tammy_FailRefitFidVolume);
   fChain->SetBranchAddress("tammy_FailShortOutTrack", &tammy_FailShortOutTrack, &b_tammy_FailShortOutTrack);
   fChain->SetBranchAddress("tammy_NoInteractionVertex", &tammy_NoInteractionVertex, &b_tammy_NoInteractionVertex);
   fChain->SetBranchAddress("tammy_NullVertex", &tammy_NullVertex, &b_tammy_NullVertex);
   fChain->SetBranchAddress("tammy_UnattachedProngsWithTracks", &tammy_UnattachedProngsWithTracks, &b_tammy_UnattachedProngsWithTracks);
   fChain->SetBranchAddress("tammy_classification", &tammy_classification, &b_tammy_classification);
   fChain->SetBranchAddress("tammy_genie_n_charms", &tammy_genie_n_charms, &b_tammy_genie_n_charms);
   fChain->SetBranchAddress("tammy_genie_n_heavy_baryons", &tammy_genie_n_heavy_baryons, &b_tammy_genie_n_heavy_baryons);
   fChain->SetBranchAddress("tammy_genie_n_kaons", &tammy_genie_n_kaons, &b_tammy_genie_n_kaons);
   fChain->SetBranchAddress("tammy_genie_n_mesons", &tammy_genie_n_mesons, &b_tammy_genie_n_mesons);
   fChain->SetBranchAddress("tammy_genie_n_muons", &tammy_genie_n_muons, &b_tammy_genie_n_muons);
   fChain->SetBranchAddress("tammy_genie_n_neutrinos", &tammy_genie_n_neutrinos, &b_tammy_genie_n_neutrinos);
   fChain->SetBranchAddress("tammy_genie_n_neutrons", &tammy_genie_n_neutrons, &b_tammy_genie_n_neutrons);
   fChain->SetBranchAddress("tammy_genie_n_others", &tammy_genie_n_others, &b_tammy_genie_n_others);
   fChain->SetBranchAddress("tammy_genie_n_particles", &tammy_genie_n_particles, &b_tammy_genie_n_particles);
   fChain->SetBranchAddress("tammy_genie_n_photons", &tammy_genie_n_photons, &b_tammy_genie_n_photons);
   fChain->SetBranchAddress("tammy_genie_n_pi_zeros", &tammy_genie_n_pi_zeros, &b_tammy_genie_n_pi_zeros);
   fChain->SetBranchAddress("tammy_genie_n_pions", &tammy_genie_n_pions, &b_tammy_genie_n_pions);
   fChain->SetBranchAddress("tammy_genie_n_protons", &tammy_genie_n_protons, &b_tammy_genie_n_protons);
   fChain->SetBranchAddress("tammy_intraNukeDeltaPlusPlusDecay", &tammy_intraNukeDeltaPlusPlusDecay, &b_tammy_intraNukeDeltaPlusPlusDecay);
   fChain->SetBranchAddress("tammy_intraNukeNParticles", &tammy_intraNukeNParticles, &b_tammy_intraNukeNParticles);
   fChain->SetBranchAddress("tammy_intraNukeNeutronQuasiElasticScatter", &tammy_intraNukeNeutronQuasiElasticScatter, &b_tammy_intraNukeNeutronQuasiElasticScatter);
   fChain->SetBranchAddress("tammy_intraNukeOtherProcess", &tammy_intraNukeOtherProcess, &b_tammy_intraNukeOtherProcess);
   fChain->SetBranchAddress("tammy_muon_enters_front", &tammy_muon_enters_front, &b_tammy_muon_enters_front);
   fChain->SetBranchAddress("tammy_n_odClusters", &tammy_n_odClusters, &b_tammy_n_odClusters);
   fChain->SetBranchAddress("tammy_n_odClustersWithTimeCut", &tammy_n_odClustersWithTimeCut, &b_tammy_n_odClustersWithTimeCut);
   fChain->SetBranchAddress("tammy_passVertexZCut", &tammy_passVertexZCut, &b_tammy_passVertexZCut);
   fChain->SetBranchAddress("tammy_proton_enters_front", &tammy_proton_enters_front, &b_tammy_proton_enters_front);
   fChain->SetBranchAddress("tammy_timeSlice", &tammy_timeSlice, &b_tammy_timeSlice);
   fChain->SetBranchAddress("tammy_vtx_fit_converged", &tammy_vtx_fit_converged, &b_tammy_vtx_fit_converged);
   fChain->SetBranchAddress("tfiducial", &tfiducial, &b_tfiducial);
   fChain->SetBranchAddress("vertex_count", &vertex_count, &b_vertex_count);
   fChain->SetBranchAddress("vertex_count2", &vertex_count2, &b_vertex_count2);
   fChain->SetBranchAddress("Dispersed_blob_energy", &Dispersed_blob_energy, &b_Dispersed_blob_energy);
   fChain->SetBranchAddress("Erec", &Erec, &b_Erec);
   fChain->SetBranchAddress("Erec2", &Erec2, &b_Erec2);
   fChain->SetBranchAddress("Filament_Vertex_energy", &Filament_Vertex_energy, &b_Filament_Vertex_energy);
   fChain->SetBranchAddress("Muon_blob_energy", &Muon_blob_energy, &b_Muon_blob_energy);
   fChain->SetBranchAddress("Q2", &Q2, &b_Q2);
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
   fChain->SetBranchAddress("Rejected_blob_vis_energy", &Rejected_blob_vis_energy, &b_Rejected_blob_vis_energy);
   fChain->SetBranchAddress("Sphere_Vertex_energy", &Sphere_Vertex_energy, &b_Sphere_Vertex_energy);
   fChain->SetBranchAddress("Tn", &Tn, &b_Tn);
   fChain->SetBranchAddress("Tn2", &Tn2, &b_Tn2);
   fChain->SetBranchAddress("Vertex_blob_energy", &Vertex_blob_energy, &b_Vertex_blob_energy);
   fChain->SetBranchAddress("W", &W, &b_W);
   fChain->SetBranchAddress("W2", &W2, &b_W2);
   fChain->SetBranchAddress("blob_fval_1", &blob_fval_1, &b_blob_fval_1);
   fChain->SetBranchAddress("blob_fval_2", &blob_fval_2, &b_blob_fval_2);
   fChain->SetBranchAddress("ecalevis", &ecalevis, &b_ecalevis);
   fChain->SetBranchAddress("energy_from_mc", &energy_from_mc, &b_energy_from_mc);
   fChain->SetBranchAddress("energy_from_mc_fraction", &energy_from_mc_fraction, &b_energy_from_mc_fraction);
   fChain->SetBranchAddress("energy_from_mc_fraction_of_highest", &energy_from_mc_fraction_of_highest, &b_energy_from_mc_fraction_of_highest);
   fChain->SetBranchAddress("g1blob_edge_distance", &g1blob_edge_distance, &b_g1blob_edge_distance);
   fChain->SetBranchAddress("g1blob_minsep", &g1blob_minsep, &b_g1blob_minsep);
   fChain->SetBranchAddress("g1blob_vtx_distance", &g1blob_vtx_distance, &b_g1blob_vtx_distance);
   fChain->SetBranchAddress("g1convdist", &g1convdist, &b_g1convdist);
   fChain->SetBranchAddress("g1e", &g1e, &b_g1e);
   fChain->SetBranchAddress("g1e0", &g1e0, &b_g1e0);
   fChain->SetBranchAddress("g1ecaledep", &g1ecaledep, &b_g1ecaledep);
   fChain->SetBranchAddress("g1ecalo", &g1ecalo, &b_g1ecalo);
   fChain->SetBranchAddress("g1g1evis", &g1g1evis, &b_g1g1evis);
   fChain->SetBranchAddress("g1g2evis", &g1g2evis, &b_g1g2evis);
   fChain->SetBranchAddress("g1gmevis", &g1gmevis, &b_g1gmevis);
   fChain->SetBranchAddress("g1hcaledep", &g1hcaledep, &b_g1hcaledep);
   fChain->SetBranchAddress("g1idetedep", &g1idetedep, &b_g1idetedep);
   fChain->SetBranchAddress("g1mostevisfrac", &g1mostevisfrac, &b_g1mostevisfrac);
   fChain->SetBranchAddress("g1muevis", &g1muevis, &b_g1muevis);
   fChain->SetBranchAddress("g1neutronevis", &g1neutronevis, &b_g1neutronevis);
   fChain->SetBranchAddress("g1nukeedep", &g1nukeedep, &b_g1nukeedep);
   fChain->SetBranchAddress("g1odetedep", &g1odetedep, &b_g1odetedep);
   fChain->SetBranchAddress("g1otherevis", &g1otherevis, &b_g1otherevis);
   fChain->SetBranchAddress("g1othersubdetedep", &g1othersubdetedep, &b_g1othersubdetedep);
   fChain->SetBranchAddress("g1phi", &g1phi, &b_g1phi);
   fChain->SetBranchAddress("g1phi0", &g1phi0, &b_g1phi0);
   fChain->SetBranchAddress("g1pi0evis", &g1pi0evis, &b_g1pi0evis);
   fChain->SetBranchAddress("g1pimevis", &g1pimevis, &b_g1pimevis);
   fChain->SetBranchAddress("g1pipevis", &g1pipevis, &b_g1pipevis);
   fChain->SetBranchAddress("g1protonevis", &g1protonevis, &b_g1protonevis);
   fChain->SetBranchAddress("g1sharedevis", &g1sharedevis, &b_g1sharedevis);
   fChain->SetBranchAddress("g1sideedep", &g1sideedep, &b_g1sideedep);
   fChain->SetBranchAddress("g1theta", &g1theta, &b_g1theta);
   fChain->SetBranchAddress("g1theta0", &g1theta0, &b_g1theta0);
   fChain->SetBranchAddress("g1totalevis", &g1totalevis, &b_g1totalevis);
   fChain->SetBranchAddress("g1trkredep", &g1trkredep, &b_g1trkredep);
   fChain->SetBranchAddress("g2blob_edge_distance", &g2blob_edge_distance, &b_g2blob_edge_distance);
   fChain->SetBranchAddress("g2blob_minsep", &g2blob_minsep, &b_g2blob_minsep);
   fChain->SetBranchAddress("g2blob_vtx_distance", &g2blob_vtx_distance, &b_g2blob_vtx_distance);
   fChain->SetBranchAddress("g2convdist", &g2convdist, &b_g2convdist);
   fChain->SetBranchAddress("g2e", &g2e, &b_g2e);
   fChain->SetBranchAddress("g2e0", &g2e0, &b_g2e0);
   fChain->SetBranchAddress("g2ecaledep", &g2ecaledep, &b_g2ecaledep);
   fChain->SetBranchAddress("g2ecalo", &g2ecalo, &b_g2ecalo);
   fChain->SetBranchAddress("g2g1evis", &g2g1evis, &b_g2g1evis);
   fChain->SetBranchAddress("g2g2evis", &g2g2evis, &b_g2g2evis);
   fChain->SetBranchAddress("g2gmevis", &g2gmevis, &b_g2gmevis);
   fChain->SetBranchAddress("g2hcaledep", &g2hcaledep, &b_g2hcaledep);
   fChain->SetBranchAddress("g2idetedep", &g2idetedep, &b_g2idetedep);
   fChain->SetBranchAddress("g2mostevisfrac", &g2mostevisfrac, &b_g2mostevisfrac);
   fChain->SetBranchAddress("g2muevis", &g2muevis, &b_g2muevis);
   fChain->SetBranchAddress("g2neutronevis", &g2neutronevis, &b_g2neutronevis);
   fChain->SetBranchAddress("g2nukeedep", &g2nukeedep, &b_g2nukeedep);
   fChain->SetBranchAddress("g2odetedep", &g2odetedep, &b_g2odetedep);
   fChain->SetBranchAddress("g2otherevis", &g2otherevis, &b_g2otherevis);
   fChain->SetBranchAddress("g2othersubdetedep", &g2othersubdetedep, &b_g2othersubdetedep);
   fChain->SetBranchAddress("g2phi", &g2phi, &b_g2phi);
   fChain->SetBranchAddress("g2phi0", &g2phi0, &b_g2phi0);
   fChain->SetBranchAddress("g2pi0evis", &g2pi0evis, &b_g2pi0evis);
   fChain->SetBranchAddress("g2pimevis", &g2pimevis, &b_g2pimevis);
   fChain->SetBranchAddress("g2pipevis", &g2pipevis, &b_g2pipevis);
   fChain->SetBranchAddress("g2protonevis", &g2protonevis, &b_g2protonevis);
   fChain->SetBranchAddress("g2sharedevis", &g2sharedevis, &b_g2sharedevis);
   fChain->SetBranchAddress("g2sideedep", &g2sideedep, &b_g2sideedep);
   fChain->SetBranchAddress("g2theta", &g2theta, &b_g2theta);
   fChain->SetBranchAddress("g2theta0", &g2theta0, &b_g2theta0);
   fChain->SetBranchAddress("g2totalevis", &g2totalevis, &b_g2totalevis);
   fChain->SetBranchAddress("g2trkredep", &g2trkredep, &b_g2trkredep);
   fChain->SetBranchAddress("hcalevis", &hcalevis, &b_hcalevis);
   fChain->SetBranchAddress("mgg", &mgg, &b_mgg);
   fChain->SetBranchAddress("minos_trk_eqp", &minos_trk_eqp, &b_minos_trk_eqp);
   fChain->SetBranchAddress("minos_trk_fit_pass", &minos_trk_fit_pass, &b_minos_trk_fit_pass);
   fChain->SetBranchAddress("minos_trk_p", &minos_trk_p, &b_minos_trk_p);
   fChain->SetBranchAddress("minos_trk_p_curvature", &minos_trk_p_curvature, &b_minos_trk_p_curvature);
   fChain->SetBranchAddress("minos_trk_p_range", &minos_trk_p_range, &b_minos_trk_p_range);
   fChain->SetBranchAddress("minos_trk_qp", &minos_trk_qp, &b_minos_trk_qp);
   fChain->SetBranchAddress("muon_phi", &muon_phi, &b_muon_phi);
   fChain->SetBranchAddress("muon_theta", &muon_theta, &b_muon_theta);
   fChain->SetBranchAddress("muon_thetaX", &muon_thetaX, &b_muon_thetaX);
   fChain->SetBranchAddress("muon_thetaY", &muon_thetaY, &b_muon_thetaY);
   fChain->SetBranchAddress("neutronecaledep", &neutronecaledep, &b_neutronecaledep);
   fChain->SetBranchAddress("neutronecalo", &neutronecalo, &b_neutronecalo);
   fChain->SetBranchAddress("neutronhcaledep", &neutronhcaledep, &b_neutronhcaledep);
   fChain->SetBranchAddress("neutronidetedep", &neutronidetedep, &b_neutronidetedep);
   fChain->SetBranchAddress("neutronnukeedep", &neutronnukeedep, &b_neutronnukeedep);
   fChain->SetBranchAddress("neutronodetedep", &neutronodetedep, &b_neutronodetedep);
   fChain->SetBranchAddress("neutronothersubdetedep", &neutronothersubdetedep, &b_neutronothersubdetedep);
   fChain->SetBranchAddress("neutronsideedep", &neutronsideedep, &b_neutronsideedep);
   fChain->SetBranchAddress("neutrontrkredep", &neutrontrkredep, &b_neutrontrkredep);
   fChain->SetBranchAddress("nke", &nke, &b_nke);
   fChain->SetBranchAddress("ntgtevis", &ntgtevis, &b_ntgtevis);
   fChain->SetBranchAddress("oangle", &oangle, &b_oangle);
   fChain->SetBranchAddress("oangle0", &oangle0, &b_oangle0);
   fChain->SetBranchAddress("oangle0x", &oangle0x, &b_oangle0x);
   fChain->SetBranchAddress("od_downstreamFrame", &od_downstreamFrame, &b_od_downstreamFrame);
   fChain->SetBranchAddress("od_downstreamFrame_z", &od_downstreamFrame_z, &b_od_downstreamFrame_z);
   fChain->SetBranchAddress("od_highStory", &od_highStory, &b_od_highStory);
   fChain->SetBranchAddress("od_highStory_t", &od_highStory_t, &b_od_highStory_t);
   fChain->SetBranchAddress("od_lowStory", &od_lowStory, &b_od_lowStory);
   fChain->SetBranchAddress("od_lowStory_t", &od_lowStory_t, &b_od_lowStory_t);
   fChain->SetBranchAddress("od_maxEnergy", &od_maxEnergy, &b_od_maxEnergy);
   fChain->SetBranchAddress("od_upstreamFrame", &od_upstreamFrame, &b_od_upstreamFrame);
   fChain->SetBranchAddress("od_upstreamFrame_z", &od_upstreamFrame_z, &b_od_upstreamFrame_z);
   fChain->SetBranchAddress("otherevis", &otherevis, &b_otherevis);
   fChain->SetBranchAddress("phys_energy_dispersed", &phys_energy_dispersed, &b_phys_energy_dispersed);
   fChain->SetBranchAddress("phys_energy_in_road_downstream", &phys_energy_in_road_downstream, &b_phys_energy_in_road_downstream);
   fChain->SetBranchAddress("phys_energy_in_road_upstream", &phys_energy_in_road_upstream, &b_phys_energy_in_road_upstream);
   fChain->SetBranchAddress("phys_energy_unattached", &phys_energy_unattached, &b_phys_energy_unattached);
   fChain->SetBranchAddress("pi0_evis_dispersed_blob", &pi0_evis_dispersed_blob, &b_pi0_evis_dispersed_blob);
   fChain->SetBranchAddress("pi0_evis_muon_blob", &pi0_evis_muon_blob, &b_pi0_evis_muon_blob);
   fChain->SetBranchAddress("pi0_evis_outtime_blob", &pi0_evis_outtime_blob, &b_pi0_evis_outtime_blob);
   fChain->SetBranchAddress("pi0_evis_vtx_blob", &pi0_evis_vtx_blob, &b_pi0_evis_vtx_blob);
   fChain->SetBranchAddress("pi0_evisfrac_dispersed_blob", &pi0_evisfrac_dispersed_blob, &b_pi0_evisfrac_dispersed_blob);
   fChain->SetBranchAddress("pi0_evisfrac_muon_blob", &pi0_evisfrac_muon_blob, &b_pi0_evisfrac_muon_blob);
   fChain->SetBranchAddress("pi0_evisfrac_outtime_blob", &pi0_evisfrac_outtime_blob, &b_pi0_evisfrac_outtime_blob);
   fChain->SetBranchAddress("pi0_evisfrac_vtx_blob", &pi0_evisfrac_vtx_blob, &b_pi0_evisfrac_vtx_blob);
   fChain->SetBranchAddress("pi0ecaledep", &pi0ecaledep, &b_pi0ecaledep);
   fChain->SetBranchAddress("pi0ecalo", &pi0ecalo, &b_pi0ecalo);
   fChain->SetBranchAddress("pi0hcaledep", &pi0hcaledep, &b_pi0hcaledep);
   fChain->SetBranchAddress("pi0idetedep", &pi0idetedep, &b_pi0idetedep);
   fChain->SetBranchAddress("pi0nukeedep", &pi0nukeedep, &b_pi0nukeedep);
   fChain->SetBranchAddress("pi0odetedep", &pi0odetedep, &b_pi0odetedep);
   fChain->SetBranchAddress("pi0othersubdetedep", &pi0othersubdetedep, &b_pi0othersubdetedep);
   fChain->SetBranchAddress("pi0sideedep", &pi0sideedep, &b_pi0sideedep);
   fChain->SetBranchAddress("pi0trkredep", &pi0trkredep, &b_pi0trkredep);
   fChain->SetBranchAddress("pienergy", &pienergy, &b_pienergy);
   fChain->SetBranchAddress("pienergy0", &pienergy0, &b_pienergy0);
   fChain->SetBranchAddress("pimlength", &pimlength, &b_pimlength);
   fChain->SetBranchAddress("piphi", &piphi, &b_piphi);
   fChain->SetBranchAddress("piphi0", &piphi0, &b_piphi0);
   fChain->SetBranchAddress("piplength", &piplength, &b_piplength);
   fChain->SetBranchAddress("pitheta", &pitheta, &b_pitheta);
   fChain->SetBranchAddress("pitheta0", &pitheta0, &b_pitheta0);
   fChain->SetBranchAddress("pke", &pke, &b_pke);
   fChain->SetBranchAddress("prim_vtx_smallest_opening_angle", &prim_vtx_smallest_opening_angle, &b_prim_vtx_smallest_opening_angle);
   fChain->SetBranchAddress("tammy_endPointEnergy", &tammy_endPointEnergy, &b_tammy_endPointEnergy);
   fChain->SetBranchAddress("tammy_hadronic_energy", &tammy_hadronic_energy, &b_tammy_hadronic_energy);
   fChain->SetBranchAddress("tammy_intraNukeProtonMomentum", &tammy_intraNukeProtonMomentum, &b_tammy_intraNukeProtonMomentum);
   fChain->SetBranchAddress("tammy_isolatedEnergy", &tammy_isolatedEnergy, &b_tammy_isolatedEnergy);
   fChain->SetBranchAddress("tammy_isolatedEnergy_ecal", &tammy_isolatedEnergy_ecal, &b_tammy_isolatedEnergy_ecal);
   fChain->SetBranchAddress("tammy_isolatedEnergy_hcal", &tammy_isolatedEnergy_hcal, &b_tammy_isolatedEnergy_hcal);
   fChain->SetBranchAddress("tammy_isolatedEnergy_targets", &tammy_isolatedEnergy_targets, &b_tammy_isolatedEnergy_targets);
   fChain->SetBranchAddress("tammy_isolatedEnergy_tracker", &tammy_isolatedEnergy_tracker, &b_tammy_isolatedEnergy_tracker);
   fChain->SetBranchAddress("tammy_muonFuzzEnergy", &tammy_muonFuzzEnergy, &b_tammy_muonFuzzEnergy);
   fChain->SetBranchAddress("tammy_odEnergy", &tammy_odEnergy, &b_tammy_odEnergy);
   fChain->SetBranchAddress("tammy_odEnergyWithTimeCut", &tammy_odEnergyWithTimeCut, &b_tammy_odEnergyWithTimeCut);
   fChain->SetBranchAddress("tammy_primaryVertexEnergy", &tammy_primaryVertexEnergy, &b_tammy_primaryVertexEnergy);
   fChain->SetBranchAddress("tammy_protonFuzzEnergy", &tammy_protonFuzzEnergy, &b_tammy_protonFuzzEnergy);
   fChain->SetBranchAddress("tammy_secondaryVertexEnergy", &tammy_secondaryVertexEnergy, &b_tammy_secondaryVertexEnergy);
   fChain->SetBranchAddress("tammy_vtx_fit_chi2", &tammy_vtx_fit_chi2, &b_tammy_vtx_fit_chi2);
   fChain->SetBranchAddress("totalevis", &totalevis, &b_totalevis);
   fChain->SetBranchAddress("trkrevis", &trkrevis, &b_trkrevis);
   fChain->SetBranchAddress("anglescan_blob_nc_sz", &anglescan_blob_nc_sz, &b_anglescan_blob_nc_sz);
   fChain->SetBranchAddress("anglescan_blob_nc", &anglescan_blob_nc, &b_anglescan_blob_nc);
   fChain->SetBranchAddress("anglescan_blob_ncu_sz", &anglescan_blob_ncu_sz, &b_anglescan_blob_ncu_sz);
   fChain->SetBranchAddress("anglescan_blob_ncu", &anglescan_blob_ncu, &b_anglescan_blob_ncu);
   fChain->SetBranchAddress("anglescan_blob_ncv_sz", &anglescan_blob_ncv_sz, &b_anglescan_blob_ncv_sz);
   fChain->SetBranchAddress("anglescan_blob_ncv", &anglescan_blob_ncv, &b_anglescan_blob_ncv);
   fChain->SetBranchAddress("anglescan_blob_ncx_sz", &anglescan_blob_ncx_sz, &b_anglescan_blob_ncx_sz);
   fChain->SetBranchAddress("anglescan_blob_ncx", &anglescan_blob_ncx, &b_anglescan_blob_ncx);
   fChain->SetBranchAddress("anglescan_blob_nd_sz", &anglescan_blob_nd_sz, &b_anglescan_blob_nd_sz);
   fChain->SetBranchAddress("anglescan_blob_nd", &anglescan_blob_nd, &b_anglescan_blob_nd);
   fChain->SetBranchAddress("anglescan_blob_ndu_sz", &anglescan_blob_ndu_sz, &b_anglescan_blob_ndu_sz);
   fChain->SetBranchAddress("anglescan_blob_ndu", &anglescan_blob_ndu, &b_anglescan_blob_ndu);
   fChain->SetBranchAddress("anglescan_blob_ndv_sz", &anglescan_blob_ndv_sz, &b_anglescan_blob_ndv_sz);
   fChain->SetBranchAddress("anglescan_blob_ndv", &anglescan_blob_ndv, &b_anglescan_blob_ndv);
   fChain->SetBranchAddress("anglescan_blob_ndx_sz", &anglescan_blob_ndx_sz, &b_anglescan_blob_ndx_sz);
   fChain->SetBranchAddress("anglescan_blob_ndx", &anglescan_blob_ndx, &b_anglescan_blob_ndx);
   fChain->SetBranchAddress("anglescan_cand_nc_sz", &anglescan_cand_nc_sz, &b_anglescan_cand_nc_sz);
   fChain->SetBranchAddress("anglescan_cand_nc", &anglescan_cand_nc, &b_anglescan_cand_nc);
   fChain->SetBranchAddress("anglescan_cand_ncu_sz", &anglescan_cand_ncu_sz, &b_anglescan_cand_ncu_sz);
   fChain->SetBranchAddress("anglescan_cand_ncu", &anglescan_cand_ncu, &b_anglescan_cand_ncu);
   fChain->SetBranchAddress("anglescan_cand_ncv_sz", &anglescan_cand_ncv_sz, &b_anglescan_cand_ncv_sz);
   fChain->SetBranchAddress("anglescan_cand_ncv", &anglescan_cand_ncv, &b_anglescan_cand_ncv);
   fChain->SetBranchAddress("anglescan_cand_ncx_sz", &anglescan_cand_ncx_sz, &b_anglescan_cand_ncx_sz);
   fChain->SetBranchAddress("anglescan_cand_ncx", &anglescan_cand_ncx, &b_anglescan_cand_ncx);
   fChain->SetBranchAddress("anglescan_cand_nd_sz", &anglescan_cand_nd_sz, &b_anglescan_cand_nd_sz);
   fChain->SetBranchAddress("anglescan_cand_nd", &anglescan_cand_nd, &b_anglescan_cand_nd);
   fChain->SetBranchAddress("anglescan_cand_ndu_sz", &anglescan_cand_ndu_sz, &b_anglescan_cand_ndu_sz);
   fChain->SetBranchAddress("anglescan_cand_ndu", &anglescan_cand_ndu, &b_anglescan_cand_ndu);
   fChain->SetBranchAddress("anglescan_cand_ndv_sz", &anglescan_cand_ndv_sz, &b_anglescan_cand_ndv_sz);
   fChain->SetBranchAddress("anglescan_cand_ndv", &anglescan_cand_ndv, &b_anglescan_cand_ndv);
   fChain->SetBranchAddress("anglescan_cand_ndx_sz", &anglescan_cand_ndx_sz, &b_anglescan_cand_ndx_sz);
   fChain->SetBranchAddress("anglescan_cand_ndx", &anglescan_cand_ndx, &b_anglescan_cand_ndx);
   fChain->SetBranchAddress("anglescan_candx_nc_sz", &anglescan_candx_nc_sz, &b_anglescan_candx_nc_sz);
   fChain->SetBranchAddress("anglescan_candx_nc", &anglescan_candx_nc, &b_anglescan_candx_nc);
   fChain->SetBranchAddress("anglescan_candx_nd_sz", &anglescan_candx_nd_sz, &b_anglescan_candx_nd_sz);
   fChain->SetBranchAddress("anglescan_candx_nd", &anglescan_candx_nd, &b_anglescan_candx_nd);
   fChain->SetBranchAddress("final_blob_nc_sz", &final_blob_nc_sz, &b_final_blob_nc_sz);
   fChain->SetBranchAddress("final_blob_nc", &final_blob_nc, &b_final_blob_nc);
   fChain->SetBranchAddress("final_blob_ncu_sz", &final_blob_ncu_sz, &b_final_blob_ncu_sz);
   fChain->SetBranchAddress("final_blob_ncu", &final_blob_ncu, &b_final_blob_ncu);
   fChain->SetBranchAddress("final_blob_ncv_sz", &final_blob_ncv_sz, &b_final_blob_ncv_sz);
   fChain->SetBranchAddress("final_blob_ncv", &final_blob_ncv, &b_final_blob_ncv);
   fChain->SetBranchAddress("final_blob_ncx_sz", &final_blob_ncx_sz, &b_final_blob_ncx_sz);
   fChain->SetBranchAddress("final_blob_ncx", &final_blob_ncx, &b_final_blob_ncx);
   fChain->SetBranchAddress("final_blob_nd_sz", &final_blob_nd_sz, &b_final_blob_nd_sz);
   fChain->SetBranchAddress("final_blob_nd", &final_blob_nd, &b_final_blob_nd);
   fChain->SetBranchAddress("final_blob_ndu_sz", &final_blob_ndu_sz, &b_final_blob_ndu_sz);
   fChain->SetBranchAddress("final_blob_ndu", &final_blob_ndu, &b_final_blob_ndu);
   fChain->SetBranchAddress("final_blob_ndv_sz", &final_blob_ndv_sz, &b_final_blob_ndv_sz);
   fChain->SetBranchAddress("final_blob_ndv", &final_blob_ndv, &b_final_blob_ndv);
   fChain->SetBranchAddress("final_blob_ndx_sz", &final_blob_ndx_sz, &b_final_blob_ndx_sz);
   fChain->SetBranchAddress("final_blob_ndx", &final_blob_ndx, &b_final_blob_ndx);
   fChain->SetBranchAddress("hough_blob_nc_sz", &hough_blob_nc_sz, &b_hough_blob_nc_sz);
   fChain->SetBranchAddress("hough_blob_nc", &hough_blob_nc, &b_hough_blob_nc);
   fChain->SetBranchAddress("hough_blob_ncu_sz", &hough_blob_ncu_sz, &b_hough_blob_ncu_sz);
   fChain->SetBranchAddress("hough_blob_ncu", &hough_blob_ncu, &b_hough_blob_ncu);
   fChain->SetBranchAddress("hough_blob_ncv_sz", &hough_blob_ncv_sz, &b_hough_blob_ncv_sz);
   fChain->SetBranchAddress("hough_blob_ncv", &hough_blob_ncv, &b_hough_blob_ncv);
   fChain->SetBranchAddress("hough_blob_ncx_sz", &hough_blob_ncx_sz, &b_hough_blob_ncx_sz);
   fChain->SetBranchAddress("hough_blob_ncx", &hough_blob_ncx, &b_hough_blob_ncx);
   fChain->SetBranchAddress("hough_blob_nd_sz", &hough_blob_nd_sz, &b_hough_blob_nd_sz);
   fChain->SetBranchAddress("hough_blob_nd", &hough_blob_nd, &b_hough_blob_nd);
   fChain->SetBranchAddress("hough_blob_ndu_sz", &hough_blob_ndu_sz, &b_hough_blob_ndu_sz);
   fChain->SetBranchAddress("hough_blob_ndu", &hough_blob_ndu, &b_hough_blob_ndu);
   fChain->SetBranchAddress("hough_blob_ndv_sz", &hough_blob_ndv_sz, &b_hough_blob_ndv_sz);
   fChain->SetBranchAddress("hough_blob_ndv", &hough_blob_ndv, &b_hough_blob_ndv);
   fChain->SetBranchAddress("hough_blob_ndx_sz", &hough_blob_ndx_sz, &b_hough_blob_ndx_sz);
   fChain->SetBranchAddress("hough_blob_ndx", &hough_blob_ndx, &b_hough_blob_ndx);
   fChain->SetBranchAddress("multiplicities_sz", &multiplicities_sz, &b_multiplicities_sz);
   fChain->SetBranchAddress("multiplicities", multiplicities, &b_multiplicities);
   fChain->SetBranchAddress("primary_truth_counts_sz", &primary_truth_counts_sz, &b_primary_truth_counts_sz);
   fChain->SetBranchAddress("primary_truth_counts", primary_truth_counts, &b_primary_truth_counts);
   fChain->SetBranchAddress("primary_truth_pdgs1_sz", &primary_truth_pdgs1_sz, &b_primary_truth_pdgs1_sz);
   fChain->SetBranchAddress("primary_truth_pdgs1", primary_truth_pdgs1, &b_primary_truth_pdgs1);
   fChain->SetBranchAddress("primary_truth_pdgs2_sz", &primary_truth_pdgs2_sz, &b_primary_truth_pdgs2_sz);
   fChain->SetBranchAddress("primary_truth_pdgs2", primary_truth_pdgs2, &b_primary_truth_pdgs2);
   fChain->SetBranchAddress("primary_truth_pdgs3_sz", &primary_truth_pdgs3_sz, &b_primary_truth_pdgs3_sz);
   fChain->SetBranchAddress("primary_truth_pdgs3", primary_truth_pdgs3, &b_primary_truth_pdgs3);
   fChain->SetBranchAddress("tammy_has_michel_category_sz", &tammy_has_michel_category_sz, &b_tammy_has_michel_category_sz);
   fChain->SetBranchAddress("tammy_has_michel_category", tammy_has_michel_category, &b_tammy_has_michel_category);
   fChain->SetBranchAddress("tammy_has_michel_in_vertex_point_sz", &tammy_has_michel_in_vertex_point_sz, &b_tammy_has_michel_in_vertex_point_sz);
   fChain->SetBranchAddress("tammy_has_michel_in_vertex_point", &tammy_has_michel_in_vertex_point, &b_tammy_has_michel_in_vertex_point);
   fChain->SetBranchAddress("tammy_has_michel_ndigits_sz", &tammy_has_michel_ndigits_sz, &b_tammy_has_michel_ndigits_sz);
   fChain->SetBranchAddress("tammy_has_michel_ndigits", &tammy_has_michel_ndigits, &b_tammy_has_michel_ndigits);
   fChain->SetBranchAddress("tammy_has_michel_vertex_type_sz", &tammy_has_michel_vertex_type_sz, &b_tammy_has_michel_vertex_type_sz);
   fChain->SetBranchAddress("tammy_has_michel_vertex_type", &tammy_has_michel_vertex_type, &b_tammy_has_michel_vertex_type);
   fChain->SetBranchAddress("RE_photon_direction_1_sz", &RE_photon_direction_1_sz, &b_RE_photon_direction_1_sz);
   fChain->SetBranchAddress("RE_photon_direction_1", &RE_photon_direction_1, &b_RE_photon_direction_1);
   fChain->SetBranchAddress("RE_photon_direction_2_sz", &RE_photon_direction_2_sz, &b_RE_photon_direction_2_sz);
   fChain->SetBranchAddress("RE_photon_direction_2", &RE_photon_direction_2, &b_RE_photon_direction_2);
   fChain->SetBranchAddress("RE_photon_vertex_1_sz", &RE_photon_vertex_1_sz, &b_RE_photon_vertex_1_sz);
   fChain->SetBranchAddress("RE_photon_vertex_1", &RE_photon_vertex_1, &b_RE_photon_vertex_1);
   fChain->SetBranchAddress("RE_photon_vertex_2_sz", &RE_photon_vertex_2_sz, &b_RE_photon_vertex_2_sz);
   fChain->SetBranchAddress("RE_photon_vertex_2", &RE_photon_vertex_2, &b_RE_photon_vertex_2);
   fChain->SetBranchAddress("deviations_sz", &deviations_sz, &b_deviations_sz);
   fChain->SetBranchAddress("deviations", deviations, &b_deviations);
   fChain->SetBranchAddress("g1convpos_sz", &g1convpos_sz, &b_g1convpos_sz);
   fChain->SetBranchAddress("g1convpos", g1convpos, &b_g1convpos);
   fChain->SetBranchAddress("g1mom_sz", &g1mom_sz, &b_g1mom_sz);
   fChain->SetBranchAddress("g1mom", &g1mom, &b_g1mom);
   fChain->SetBranchAddress("g1mom0_sz", &g1mom0_sz, &b_g1mom0_sz);
   fChain->SetBranchAddress("g1mom0", g1mom0, &b_g1mom0);
   fChain->SetBranchAddress("g2convpos_sz", &g2convpos_sz, &b_g2convpos_sz);
   fChain->SetBranchAddress("g2convpos", g2convpos, &b_g2convpos);
   fChain->SetBranchAddress("g2mom_sz", &g2mom_sz, &b_g2mom_sz);
   fChain->SetBranchAddress("g2mom", &g2mom, &b_g2mom);
   fChain->SetBranchAddress("g2mom0_sz", &g2mom0_sz, &b_g2mom0_sz);
   fChain->SetBranchAddress("g2mom0", g2mom0, &b_g2mom0);
   fChain->SetBranchAddress("michel_mom_sz", &michel_mom_sz, &b_michel_mom_sz);
   fChain->SetBranchAddress("michel_mom", michel_mom, &b_michel_mom);
   fChain->SetBranchAddress("michel_pos_sz", &michel_pos_sz, &b_michel_pos_sz);
   fChain->SetBranchAddress("michel_pos", michel_pos, &b_michel_pos);
   fChain->SetBranchAddress("mumom_sz", &mumom_sz, &b_mumom_sz);
   fChain->SetBranchAddress("mumom", &mumom, &b_mumom);
   fChain->SetBranchAddress("od_distanceBlobTower_sz", &od_distanceBlobTower_sz, &b_od_distanceBlobTower_sz);
   fChain->SetBranchAddress("od_distanceBlobTower", &od_distanceBlobTower, &b_od_distanceBlobTower);
   fChain->SetBranchAddress("od_idBlobTime_sz", &od_idBlobTime_sz, &b_od_idBlobTime_sz);
   fChain->SetBranchAddress("od_idBlobTime", &od_idBlobTime, &b_od_idBlobTime);
   fChain->SetBranchAddress("od_towerEnergy_sz", &od_towerEnergy_sz, &b_od_towerEnergy_sz);
   fChain->SetBranchAddress("od_towerEnergy", &od_towerEnergy, &b_od_towerEnergy);
   fChain->SetBranchAddress("od_towerNClusters_sz", &od_towerNClusters_sz, &b_od_towerNClusters_sz);
   fChain->SetBranchAddress("od_towerNClusters", &od_towerNClusters, &b_od_towerNClusters);
   fChain->SetBranchAddress("od_towerTime_sz", &od_towerTime_sz, &b_od_towerTime_sz);
   fChain->SetBranchAddress("od_towerTime", &od_towerTime, &b_od_towerTime);
   fChain->SetBranchAddress("od_towerTimeBlobMuon_sz", &od_towerTimeBlobMuon_sz, &b_od_towerTimeBlobMuon_sz);
   fChain->SetBranchAddress("od_towerTimeBlobMuon", &od_towerTimeBlobMuon, &b_od_towerTimeBlobMuon);
   fChain->SetBranchAddress("od_towerTimeBlobOD_sz", &od_towerTimeBlobOD_sz, &b_od_towerTimeBlobOD_sz);
   fChain->SetBranchAddress("od_towerTimeBlobOD", &od_towerTimeBlobOD, &b_od_towerTimeBlobOD);
   fChain->SetBranchAddress("pimom_sz", &pimom_sz, &b_pimom_sz);
   fChain->SetBranchAddress("pimom", &pimom, &b_pimom);
   fChain->SetBranchAddress("pimom0_sz", &pimom0_sz, &b_pimom0_sz);
   fChain->SetBranchAddress("pimom0", pimom0, &b_pimom0);
   fChain->SetBranchAddress("primary_separations_sz", &primary_separations_sz, &b_primary_separations_sz);
   fChain->SetBranchAddress("primary_separations", primary_separations, &b_primary_separations);
   fChain->SetBranchAddress("primary_trklengths_sz", &primary_trklengths_sz, &b_primary_trklengths_sz);
   fChain->SetBranchAddress("primary_trklengths", primary_trklengths, &b_primary_trklengths);
   fChain->SetBranchAddress("primary_truth_fractions1_sz", &primary_truth_fractions1_sz, &b_primary_truth_fractions1_sz);
   fChain->SetBranchAddress("primary_truth_fractions1", primary_truth_fractions1, &b_primary_truth_fractions1);
   fChain->SetBranchAddress("primary_truth_fractions2_sz", &primary_truth_fractions2_sz, &b_primary_truth_fractions2_sz);
   fChain->SetBranchAddress("primary_truth_fractions2", primary_truth_fractions2, &b_primary_truth_fractions2);
   fChain->SetBranchAddress("primary_truth_fractions3_sz", &primary_truth_fractions3_sz, &b_primary_truth_fractions3_sz);
   fChain->SetBranchAddress("primary_truth_fractions3", primary_truth_fractions3, &b_primary_truth_fractions3);
   fChain->SetBranchAddress("primary_truth_shareds_sz", &primary_truth_shareds_sz, &b_primary_truth_shareds_sz);
   fChain->SetBranchAddress("primary_truth_shareds", primary_truth_shareds, &b_primary_truth_shareds);
   fChain->SetBranchAddress("tammy_fit_vtx", tammy_fit_vtx, &b_tammy_fit_vtx);
   fChain->SetBranchAddress("tammy_has_michel_distance_sz", &tammy_has_michel_distance_sz, &b_tammy_has_michel_distance_sz);
   fChain->SetBranchAddress("tammy_has_michel_distance", &tammy_has_michel_distance, &b_tammy_has_michel_distance);
   fChain->SetBranchAddress("tammy_has_michel_energy_sz", &tammy_has_michel_energy_sz, &b_tammy_has_michel_energy_sz);
   fChain->SetBranchAddress("tammy_has_michel_energy", &tammy_has_michel_energy, &b_tammy_has_michel_energy);
   fChain->SetBranchAddress("tammy_has_michel_time_diff_sz", &tammy_has_michel_time_diff_sz, &b_tammy_has_michel_time_diff_sz);
   fChain->SetBranchAddress("tammy_has_michel_time_diff", &tammy_has_michel_time_diff, &b_tammy_has_michel_time_diff);
   fChain->SetBranchAddress("tammy_intraNukeProtonMomentumVec", tammy_intraNukeProtonMomentumVec, &b_tammy_intraNukeProtonMomentumVec);
   fChain->SetBranchAddress("truth_has_physics_event", &truth_has_physics_event, &b_truth_has_physics_event);
   fChain->SetBranchAddress("truth_reco_minos_match", &truth_reco_minos_match, &b_truth_reco_minos_match);
   fChain->SetBranchAddress("truth_is_fiducial", &truth_is_fiducial, &b_truth_is_fiducial);
   fChain->SetBranchAddress("truth_pass_plausible", &truth_pass_plausible, &b_truth_pass_plausible);
   fChain->SetBranchAddress("truth_is_ccpi0", &truth_is_ccpi0, &b_truth_is_ccpi0);
   fChain->SetBranchAddress("truth_is_cc1pi0", &truth_is_cc1pi0, &b_truth_is_cc1pi0);
   fChain->SetBranchAddress("truth_is_ccpi0secondary", &truth_is_ccpi0secondary, &b_truth_is_ccpi0secondary);
   fChain->SetBranchAddress("truth_is_by_pim", &truth_is_by_pim, &b_truth_is_by_pim);
   fChain->SetBranchAddress("truth_is_ccpi0x", &truth_is_ccpi0x, &b_truth_is_ccpi0x);
   fChain->SetBranchAddress("truth_is_other", &truth_is_other, &b_truth_is_other);
   fChain->SetBranchAddress("truth_tammy_has_michel_electron", &truth_tammy_has_michel_electron, &b_truth_tammy_has_michel_electron);
   fChain->SetBranchAddress("truth_MC_photon_energy_1", &truth_MC_photon_energy_1, &b_truth_MC_photon_energy_1);
   fChain->SetBranchAddress("truth_MC_photon_energy_2", &truth_MC_photon_energy_2, &b_truth_MC_photon_energy_2);
   fChain->SetBranchAddress("truth_MC_pi0_energy", &truth_MC_pi0_energy, &b_truth_MC_pi0_energy);
   fChain->SetBranchAddress("truth_MC_scalar", &truth_MC_scalar, &b_truth_MC_scalar);
   fChain->SetBranchAddress("truth_fslepton_E", &truth_fslepton_E, &b_truth_fslepton_E);
   fChain->SetBranchAddress("truth_fslepton_P", &truth_fslepton_P, &b_truth_fslepton_P);
   fChain->SetBranchAddress("truth_fslepton_T", &truth_fslepton_T, &b_truth_fslepton_T);
   fChain->SetBranchAddress("truth_fslepton_phi", &truth_fslepton_phi, &b_truth_fslepton_phi);
   fChain->SetBranchAddress("truth_fslepton_theta", &truth_fslepton_theta, &b_truth_fslepton_theta);
   fChain->SetBranchAddress("truth_fslepton_theta_x", &truth_fslepton_theta_x, &b_truth_fslepton_theta_x);
   fChain->SetBranchAddress("truth_fslepton_theta_y", &truth_fslepton_theta_y, &b_truth_fslepton_theta_y);
   fChain->SetBranchAddress("truth_MC_photon_direction_1_sz", &truth_MC_photon_direction_1_sz, &b_truth_MC_photon_direction_1_sz);
   fChain->SetBranchAddress("truth_MC_photon_direction_1", &truth_MC_photon_direction_1, &b_truth_MC_photon_direction_1);
   fChain->SetBranchAddress("truth_MC_photon_direction_2_sz", &truth_MC_photon_direction_2_sz, &b_truth_MC_photon_direction_2_sz);
   fChain->SetBranchAddress("truth_MC_photon_direction_2", &truth_MC_photon_direction_2, &b_truth_MC_photon_direction_2);
   fChain->SetBranchAddress("truth_MC_photon_vertex_1_sz", &truth_MC_photon_vertex_1_sz, &b_truth_MC_photon_vertex_1_sz);
   fChain->SetBranchAddress("truth_MC_photon_vertex_1", &truth_MC_photon_vertex_1, &b_truth_MC_photon_vertex_1);
   fChain->SetBranchAddress("truth_MC_photon_vertex_2_sz", &truth_MC_photon_vertex_2_sz, &b_truth_MC_photon_vertex_2_sz);
   fChain->SetBranchAddress("truth_MC_photon_vertex_2", &truth_MC_photon_vertex_2, &b_truth_MC_photon_vertex_2);
   fChain->SetBranchAddress("truth_MC_pi0_momentum_sz", &truth_MC_pi0_momentum_sz, &b_truth_MC_pi0_momentum_sz);
   fChain->SetBranchAddress("truth_MC_pi0_momentum", &truth_MC_pi0_momentum, &b_truth_MC_pi0_momentum);
   fChain->SetBranchAddress("genie_wgt_n_shifts", &genie_wgt_n_shifts, &b_genie_wgt_n_shifts);
   fChain->SetBranchAddress("truth_genie_wgt_AGKYxF1pi", &truth_genie_wgt_AGKYxF1pi, &b_truth_genie_wgt_AGKYxF1pi);
   fChain->SetBranchAddress("truth_genie_wgt_AhtBY", &truth_genie_wgt_AhtBY, &b_truth_genie_wgt_AhtBY);
   fChain->SetBranchAddress("truth_genie_wgt_BhtBY", &truth_genie_wgt_BhtBY, &b_truth_genie_wgt_BhtBY);
   fChain->SetBranchAddress("truth_genie_wgt_CCQEPauliSupViaFK", &truth_genie_wgt_CCQEPauliSupViaFK, &b_truth_genie_wgt_CCQEPauliSupViaFK);
   fChain->SetBranchAddress("truth_genie_wgt_CV1uBY", &truth_genie_wgt_CV1uBY, &b_truth_genie_wgt_CV1uBY);
   fChain->SetBranchAddress("truth_genie_wgt_CV2uBY", &truth_genie_wgt_CV2uBY, &b_truth_genie_wgt_CV2uBY);
   fChain->SetBranchAddress("truth_genie_wgt_EtaNCEL", &truth_genie_wgt_EtaNCEL, &b_truth_genie_wgt_EtaNCEL);
   fChain->SetBranchAddress("truth_genie_wgt_FrAbs_N", &truth_genie_wgt_FrAbs_N, &b_truth_genie_wgt_FrAbs_N);
   fChain->SetBranchAddress("truth_genie_wgt_FrAbs_pi", &truth_genie_wgt_FrAbs_pi, &b_truth_genie_wgt_FrAbs_pi);
   fChain->SetBranchAddress("truth_genie_wgt_FrCEx_N", &truth_genie_wgt_FrCEx_N, &b_truth_genie_wgt_FrCEx_N);
   fChain->SetBranchAddress("truth_genie_wgt_FrCEx_pi", &truth_genie_wgt_FrCEx_pi, &b_truth_genie_wgt_FrCEx_pi);
   fChain->SetBranchAddress("truth_genie_wgt_FrElas_N", &truth_genie_wgt_FrElas_N, &b_truth_genie_wgt_FrElas_N);
   fChain->SetBranchAddress("truth_genie_wgt_FrElas_pi", &truth_genie_wgt_FrElas_pi, &b_truth_genie_wgt_FrElas_pi);
   fChain->SetBranchAddress("truth_genie_wgt_FrInel_N", &truth_genie_wgt_FrInel_N, &b_truth_genie_wgt_FrInel_N);
   fChain->SetBranchAddress("truth_genie_wgt_FrInel_pi", &truth_genie_wgt_FrInel_pi, &b_truth_genie_wgt_FrInel_pi);
   fChain->SetBranchAddress("truth_genie_wgt_FrPiProd_N", &truth_genie_wgt_FrPiProd_N, &b_truth_genie_wgt_FrPiProd_N);
   fChain->SetBranchAddress("truth_genie_wgt_FrPiProd_pi", &truth_genie_wgt_FrPiProd_pi, &b_truth_genie_wgt_FrPiProd_pi);
   fChain->SetBranchAddress("truth_genie_wgt_MFP_N", &truth_genie_wgt_MFP_N, &b_truth_genie_wgt_MFP_N);
   fChain->SetBranchAddress("truth_genie_wgt_MFP_pi", &truth_genie_wgt_MFP_pi, &b_truth_genie_wgt_MFP_pi);
   fChain->SetBranchAddress("truth_genie_wgt_MaCCQE", &truth_genie_wgt_MaCCQE, &b_truth_genie_wgt_MaCCQE);
   fChain->SetBranchAddress("truth_genie_wgt_MaCCQEshape", &truth_genie_wgt_MaCCQEshape, &b_truth_genie_wgt_MaCCQEshape);
   fChain->SetBranchAddress("truth_genie_wgt_MaNCEL", &truth_genie_wgt_MaNCEL, &b_truth_genie_wgt_MaNCEL);
   fChain->SetBranchAddress("truth_genie_wgt_MaRES", &truth_genie_wgt_MaRES, &b_truth_genie_wgt_MaRES);
   fChain->SetBranchAddress("truth_genie_wgt_MvRES", &truth_genie_wgt_MvRES, &b_truth_genie_wgt_MvRES);
   fChain->SetBranchAddress("truth_genie_wgt_NormCCQE", &truth_genie_wgt_NormCCQE, &b_truth_genie_wgt_NormCCQE);
   fChain->SetBranchAddress("truth_genie_wgt_NormCCRES", &truth_genie_wgt_NormCCRES, &b_truth_genie_wgt_NormCCRES);
   fChain->SetBranchAddress("truth_genie_wgt_NormDISCC", &truth_genie_wgt_NormDISCC, &b_truth_genie_wgt_NormDISCC);
   fChain->SetBranchAddress("truth_genie_wgt_NormNCRES", &truth_genie_wgt_NormNCRES, &b_truth_genie_wgt_NormNCRES);
   fChain->SetBranchAddress("truth_genie_wgt_RDecBR1gamma", &truth_genie_wgt_RDecBR1gamma, &b_truth_genie_wgt_RDecBR1gamma);
   fChain->SetBranchAddress("truth_genie_wgt_Rvn1pi", &truth_genie_wgt_Rvn1pi, &b_truth_genie_wgt_Rvn1pi);
   fChain->SetBranchAddress("truth_genie_wgt_Rvn2pi", &truth_genie_wgt_Rvn2pi, &b_truth_genie_wgt_Rvn2pi);
   fChain->SetBranchAddress("truth_genie_wgt_Rvp1pi", &truth_genie_wgt_Rvp1pi, &b_truth_genie_wgt_Rvp1pi);
   fChain->SetBranchAddress("truth_genie_wgt_Rvp2pi", &truth_genie_wgt_Rvp2pi, &b_truth_genie_wgt_Rvp2pi);
   fChain->SetBranchAddress("truth_genie_wgt_Theta_Delta2Npi", &truth_genie_wgt_Theta_Delta2Npi, &b_truth_genie_wgt_Theta_Delta2Npi);
   fChain->SetBranchAddress("truth_genie_wgt_VecFFCCQEshape", &truth_genie_wgt_VecFFCCQEshape, &b_truth_genie_wgt_VecFFCCQEshape);
   fChain->SetBranchAddress("truth_genie_wgt_shifts", &truth_genie_wgt_shifts, &b_truth_genie_wgt_shifts);
   fChain->SetBranchAddress("truth_tammy_has_michel_from_tammy_pion_minus_momentum_sz", &truth_tammy_has_michel_from_tammy_pion_minus_momentum_sz, &b_truth_tammy_has_michel_from_tammy_pion_minus_momentum_sz);
   fChain->SetBranchAddress("truth_tammy_has_michel_from_tammy_pion_minus_momentum", &truth_tammy_has_michel_from_tammy_pion_minus_momentum, &b_truth_tammy_has_michel_from_tammy_pion_minus_momentum);
   fChain->SetBranchAddress("truth_tammy_has_michel_from_tammy_pion_plus_momentum_sz", &truth_tammy_has_michel_from_tammy_pion_plus_momentum_sz, &b_truth_tammy_has_michel_from_tammy_pion_plus_momentum_sz);
   fChain->SetBranchAddress("truth_tammy_has_michel_from_tammy_pion_plus_momentum", &truth_tammy_has_michel_from_tammy_pion_plus_momentum, &b_truth_tammy_has_michel_from_tammy_pion_plus_momentum);
   fChain->SetBranchAddress("CCDeltaPlusAna_nuFlavor", &CCDeltaPlusAna_nuFlavor, &b_CCDeltaPlusAna_nuFlavor);
   fChain->SetBranchAddress("CCDeltaPlusAna_nuHelicity", &CCDeltaPlusAna_nuHelicity, &b_CCDeltaPlusAna_nuHelicity);
   fChain->SetBranchAddress("CCDeltaPlusAna_intCurrent", &CCDeltaPlusAna_intCurrent, &b_CCDeltaPlusAna_intCurrent);
   fChain->SetBranchAddress("CCDeltaPlusAna_intType", &CCDeltaPlusAna_intType, &b_CCDeltaPlusAna_intType);
   fChain->SetBranchAddress("CCDeltaPlusAna_E", &CCDeltaPlusAna_E, &b_CCDeltaPlusAna_E);
   fChain->SetBranchAddress("CCDeltaPlusAna_Q2", &CCDeltaPlusAna_Q2, &b_CCDeltaPlusAna_Q2);
   fChain->SetBranchAddress("CCDeltaPlusAna_x", &CCDeltaPlusAna_x, &b_CCDeltaPlusAna_x);
   fChain->SetBranchAddress("CCDeltaPlusAna_y", &CCDeltaPlusAna_y, &b_CCDeltaPlusAna_y);
   fChain->SetBranchAddress("CCDeltaPlusAna_W", &CCDeltaPlusAna_W, &b_CCDeltaPlusAna_W);
   fChain->SetBranchAddress("CCDeltaPlusAna_score", &CCDeltaPlusAna_score, &b_CCDeltaPlusAna_score);
   fChain->SetBranchAddress("CCDeltaPlusAna_leptonE", CCDeltaPlusAna_leptonE, &b_CCDeltaPlusAna_leptonE);
   fChain->SetBranchAddress("CCDeltaPlusAna_vtx", CCDeltaPlusAna_vtx, &b_CCDeltaPlusAna_vtx);
   fChain->SetBranchAddress("CCDeltaPlusAna_minos_trk_is_contained", &CCDeltaPlusAna_minos_trk_is_contained, &b_CCDeltaPlusAna_minos_trk_is_contained);
   fChain->SetBranchAddress("CCDeltaPlusAna_minos_trk_is_ok", &CCDeltaPlusAna_minos_trk_is_ok, &b_CCDeltaPlusAna_minos_trk_is_ok);
   fChain->SetBranchAddress("CCDeltaPlusAna_minos_used_range", &CCDeltaPlusAna_minos_used_range, &b_CCDeltaPlusAna_minos_used_range);
   fChain->SetBranchAddress("CCDeltaPlusAna_minos_used_curvature", &CCDeltaPlusAna_minos_used_curvature, &b_CCDeltaPlusAna_minos_used_curvature);
   fChain->SetBranchAddress("CCDeltaPlusAna_minos_trk_end_plane", &CCDeltaPlusAna_minos_trk_end_plane, &b_CCDeltaPlusAna_minos_trk_end_plane);
   fChain->SetBranchAddress("CCDeltaPlusAna_minos_trk_quality", &CCDeltaPlusAna_minos_trk_quality, &b_CCDeltaPlusAna_minos_trk_quality);
   fChain->SetBranchAddress("CCDeltaPlusAna_r_minos_trk_vtx_plane", &CCDeltaPlusAna_r_minos_trk_vtx_plane, &b_CCDeltaPlusAna_r_minos_trk_vtx_plane);
   fChain->SetBranchAddress("CCDeltaPlusAna_t_minos_trk_numFSMuons", &CCDeltaPlusAna_t_minos_trk_numFSMuons, &b_CCDeltaPlusAna_t_minos_trk_numFSMuons);
   fChain->SetBranchAddress("CCDeltaPlusAna_t_minos_trk_primFSLeptonPDG", &CCDeltaPlusAna_t_minos_trk_primFSLeptonPDG, &b_CCDeltaPlusAna_t_minos_trk_primFSLeptonPDG);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_inside_minos_partial_plane", &CCDeltaPlusAna_tammy_inside_minos_partial_plane, &b_CCDeltaPlusAna_tammy_inside_minos_partial_plane);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_isMuonInsideOD", &CCDeltaPlusAna_tammy_isMuonInsideOD, &b_CCDeltaPlusAna_tammy_isMuonInsideOD);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_isProtonInsideOD", &CCDeltaPlusAna_tammy_isProtonInsideOD, &b_CCDeltaPlusAna_tammy_isProtonInsideOD);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muon_down_hcal", &CCDeltaPlusAna_tammy_muon_down_hcal, &b_CCDeltaPlusAna_tammy_muon_down_hcal);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muon_minos_stub", &CCDeltaPlusAna_tammy_muon_minos_stub, &b_CCDeltaPlusAna_tammy_muon_minos_stub);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muon_minos_track", &CCDeltaPlusAna_tammy_muon_minos_track, &b_CCDeltaPlusAna_tammy_muon_minos_track);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muon_n_values", &CCDeltaPlusAna_tammy_muon_n_values, &b_CCDeltaPlusAna_tammy_muon_n_values);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muon_odLastFrame", &CCDeltaPlusAna_tammy_muon_odLastFrame, &b_CCDeltaPlusAna_tammy_muon_odLastFrame);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muon_odLastStory", &CCDeltaPlusAna_tammy_muon_odLastStory, &b_CCDeltaPlusAna_tammy_muon_odLastStory);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muon_od_track", &CCDeltaPlusAna_tammy_muon_od_track, &b_CCDeltaPlusAna_tammy_muon_od_track);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muon_side_ecal", &CCDeltaPlusAna_tammy_muon_side_ecal, &b_CCDeltaPlusAna_tammy_muon_side_ecal);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muon_trk_pat_history", &CCDeltaPlusAna_tammy_muon_trk_pat_history, &b_CCDeltaPlusAna_tammy_muon_trk_pat_history);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_ntrajMuonProng", &CCDeltaPlusAna_tammy_ntrajMuonProng, &b_CCDeltaPlusAna_tammy_ntrajMuonProng);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_ntrajProngProng", &CCDeltaPlusAna_tammy_ntrajProngProng, &b_CCDeltaPlusAna_tammy_ntrajProngProng);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_pOK", &CCDeltaPlusAna_tammy_pOK, &b_CCDeltaPlusAna_tammy_pOK);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_proton_kinked", &CCDeltaPlusAna_tammy_proton_kinked, &b_CCDeltaPlusAna_tammy_proton_kinked);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_proton_n_values", &CCDeltaPlusAna_tammy_proton_n_values, &b_CCDeltaPlusAna_tammy_proton_n_values);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_proton_odMatch", &CCDeltaPlusAna_tammy_proton_odMatch, &b_CCDeltaPlusAna_tammy_proton_odMatch);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_proton_trk_pat_history", &CCDeltaPlusAna_tammy_proton_trk_pat_history, &b_CCDeltaPlusAna_tammy_proton_trk_pat_history);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_targetID", &CCDeltaPlusAna_tammy_targetID, &b_CCDeltaPlusAna_tammy_targetID);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_targetZ", &CCDeltaPlusAna_tammy_targetZ, &b_CCDeltaPlusAna_tammy_targetZ);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_trajMuonProngPDG", &CCDeltaPlusAna_tammy_trajMuonProngPDG, &b_CCDeltaPlusAna_tammy_trajMuonProngPDG);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_trajMuonProngPrimary", &CCDeltaPlusAna_tammy_trajMuonProngPrimary, &b_CCDeltaPlusAna_tammy_trajMuonProngPrimary);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_trajProtonProngPDG", &CCDeltaPlusAna_tammy_trajProtonProngPDG, &b_CCDeltaPlusAna_tammy_trajProtonProngPDG);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_trajProtonProngPrimary", &CCDeltaPlusAna_tammy_trajProtonProngPrimary, &b_CCDeltaPlusAna_tammy_trajProtonProngPrimary);
   fChain->SetBranchAddress("CCDeltaPlusAna_minos_trk_bave", &CCDeltaPlusAna_minos_trk_bave, &b_CCDeltaPlusAna_minos_trk_bave);
   fChain->SetBranchAddress("CCDeltaPlusAna_minos_trk_chi2", &CCDeltaPlusAna_minos_trk_chi2, &b_CCDeltaPlusAna_minos_trk_chi2);
   fChain->SetBranchAddress("CCDeltaPlusAna_minos_trk_end_u", &CCDeltaPlusAna_minos_trk_end_u, &b_CCDeltaPlusAna_minos_trk_end_u);
   fChain->SetBranchAddress("CCDeltaPlusAna_minos_trk_end_v", &CCDeltaPlusAna_minos_trk_end_v, &b_CCDeltaPlusAna_minos_trk_end_v);
   fChain->SetBranchAddress("CCDeltaPlusAna_minos_trk_end_x", &CCDeltaPlusAna_minos_trk_end_x, &b_CCDeltaPlusAna_minos_trk_end_x);
   fChain->SetBranchAddress("CCDeltaPlusAna_minos_trk_end_y", &CCDeltaPlusAna_minos_trk_end_y, &b_CCDeltaPlusAna_minos_trk_end_y);
   fChain->SetBranchAddress("CCDeltaPlusAna_minos_trk_end_z", &CCDeltaPlusAna_minos_trk_end_z, &b_CCDeltaPlusAna_minos_trk_end_z);
   fChain->SetBranchAddress("CCDeltaPlusAna_minos_trk_eqp", &CCDeltaPlusAna_minos_trk_eqp, &b_CCDeltaPlusAna_minos_trk_eqp);
   fChain->SetBranchAddress("CCDeltaPlusAna_minos_trk_eqp_qp", &CCDeltaPlusAna_minos_trk_eqp_qp, &b_CCDeltaPlusAna_minos_trk_eqp_qp);
   fChain->SetBranchAddress("CCDeltaPlusAna_minos_trk_fit_pass", &CCDeltaPlusAna_minos_trk_fit_pass, &b_CCDeltaPlusAna_minos_trk_fit_pass);
   fChain->SetBranchAddress("CCDeltaPlusAna_minos_trk_ndf", &CCDeltaPlusAna_minos_trk_ndf, &b_CCDeltaPlusAna_minos_trk_ndf);
   fChain->SetBranchAddress("CCDeltaPlusAna_minos_trk_p", &CCDeltaPlusAna_minos_trk_p, &b_CCDeltaPlusAna_minos_trk_p);
   fChain->SetBranchAddress("CCDeltaPlusAna_minos_trk_p_curvature", &CCDeltaPlusAna_minos_trk_p_curvature, &b_CCDeltaPlusAna_minos_trk_p_curvature);
   fChain->SetBranchAddress("CCDeltaPlusAna_minos_trk_p_range", &CCDeltaPlusAna_minos_trk_p_range, &b_CCDeltaPlusAna_minos_trk_p_range);
   fChain->SetBranchAddress("CCDeltaPlusAna_minos_trk_qp", &CCDeltaPlusAna_minos_trk_qp, &b_CCDeltaPlusAna_minos_trk_qp);
   fChain->SetBranchAddress("CCDeltaPlusAna_minos_trk_vtx_x", &CCDeltaPlusAna_minos_trk_vtx_x, &b_CCDeltaPlusAna_minos_trk_vtx_x);
   fChain->SetBranchAddress("CCDeltaPlusAna_minos_trk_vtx_y", &CCDeltaPlusAna_minos_trk_vtx_y, &b_CCDeltaPlusAna_minos_trk_vtx_y);
   fChain->SetBranchAddress("CCDeltaPlusAna_minos_trk_vtx_z", &CCDeltaPlusAna_minos_trk_vtx_z, &b_CCDeltaPlusAna_minos_trk_vtx_z);
   fChain->SetBranchAddress("CCDeltaPlusAna_r_minos_trk_bdL", &CCDeltaPlusAna_r_minos_trk_bdL, &b_CCDeltaPlusAna_r_minos_trk_bdL);
   fChain->SetBranchAddress("CCDeltaPlusAna_r_minos_trk_end_dcosx", &CCDeltaPlusAna_r_minos_trk_end_dcosx, &b_CCDeltaPlusAna_r_minos_trk_end_dcosx);
   fChain->SetBranchAddress("CCDeltaPlusAna_r_minos_trk_end_dcosy", &CCDeltaPlusAna_r_minos_trk_end_dcosy, &b_CCDeltaPlusAna_r_minos_trk_end_dcosy);
   fChain->SetBranchAddress("CCDeltaPlusAna_r_minos_trk_end_dcosz", &CCDeltaPlusAna_r_minos_trk_end_dcosz, &b_CCDeltaPlusAna_r_minos_trk_end_dcosz);
   fChain->SetBranchAddress("CCDeltaPlusAna_r_minos_trk_vtx_dcosx", &CCDeltaPlusAna_r_minos_trk_vtx_dcosx, &b_CCDeltaPlusAna_r_minos_trk_vtx_dcosx);
   fChain->SetBranchAddress("CCDeltaPlusAna_r_minos_trk_vtx_dcosy", &CCDeltaPlusAna_r_minos_trk_vtx_dcosy, &b_CCDeltaPlusAna_r_minos_trk_vtx_dcosy);
   fChain->SetBranchAddress("CCDeltaPlusAna_r_minos_trk_vtx_dcosz", &CCDeltaPlusAna_r_minos_trk_vtx_dcosz, &b_CCDeltaPlusAna_r_minos_trk_vtx_dcosz);
   fChain->SetBranchAddress("CCDeltaPlusAna_t_minos_trk_primFSLepMinosInitProjPx", &CCDeltaPlusAna_t_minos_trk_primFSLepMinosInitProjPx, &b_CCDeltaPlusAna_t_minos_trk_primFSLepMinosInitProjPx);
   fChain->SetBranchAddress("CCDeltaPlusAna_t_minos_trk_primFSLepMinosInitProjPy", &CCDeltaPlusAna_t_minos_trk_primFSLepMinosInitProjPy, &b_CCDeltaPlusAna_t_minos_trk_primFSLepMinosInitProjPy);
   fChain->SetBranchAddress("CCDeltaPlusAna_t_minos_trk_primFSLepMinosInitProjPz", &CCDeltaPlusAna_t_minos_trk_primFSLepMinosInitProjPz, &b_CCDeltaPlusAna_t_minos_trk_primFSLepMinosInitProjPz);
   fChain->SetBranchAddress("CCDeltaPlusAna_t_minos_trk_primFSLepMinosInitProjX", &CCDeltaPlusAna_t_minos_trk_primFSLepMinosInitProjX, &b_CCDeltaPlusAna_t_minos_trk_primFSLepMinosInitProjX);
   fChain->SetBranchAddress("CCDeltaPlusAna_t_minos_trk_primFSLepMinosInitProjY", &CCDeltaPlusAna_t_minos_trk_primFSLepMinosInitProjY, &b_CCDeltaPlusAna_t_minos_trk_primFSLepMinosInitProjY);
   fChain->SetBranchAddress("CCDeltaPlusAna_t_minos_trk_primFSLepMinosInitProjZ", &CCDeltaPlusAna_t_minos_trk_primFSLepMinosInitProjZ, &b_CCDeltaPlusAna_t_minos_trk_primFSLepMinosInitProjZ);
   fChain->SetBranchAddress("CCDeltaPlusAna_t_minos_trk_primFSLepMnvFinalPx", &CCDeltaPlusAna_t_minos_trk_primFSLepMnvFinalPx, &b_CCDeltaPlusAna_t_minos_trk_primFSLepMnvFinalPx);
   fChain->SetBranchAddress("CCDeltaPlusAna_t_minos_trk_primFSLepMnvFinalPy", &CCDeltaPlusAna_t_minos_trk_primFSLepMnvFinalPy, &b_CCDeltaPlusAna_t_minos_trk_primFSLepMnvFinalPy);
   fChain->SetBranchAddress("CCDeltaPlusAna_t_minos_trk_primFSLepMnvFinalPz", &CCDeltaPlusAna_t_minos_trk_primFSLepMnvFinalPz, &b_CCDeltaPlusAna_t_minos_trk_primFSLepMnvFinalPz);
   fChain->SetBranchAddress("CCDeltaPlusAna_t_minos_trk_primFSLepMnvFinalX", &CCDeltaPlusAna_t_minos_trk_primFSLepMnvFinalX, &b_CCDeltaPlusAna_t_minos_trk_primFSLepMnvFinalX);
   fChain->SetBranchAddress("CCDeltaPlusAna_t_minos_trk_primFSLepMnvFinalY", &CCDeltaPlusAna_t_minos_trk_primFSLepMnvFinalY, &b_CCDeltaPlusAna_t_minos_trk_primFSLepMnvFinalY);
   fChain->SetBranchAddress("CCDeltaPlusAna_t_minos_trk_primFSLepMnvFinalZ", &CCDeltaPlusAna_t_minos_trk_primFSLepMnvFinalZ, &b_CCDeltaPlusAna_t_minos_trk_primFSLepMnvFinalZ);
   fChain->SetBranchAddress("CCDeltaPlusAna_t_minos_trk_primFSLepMnvInitPx", &CCDeltaPlusAna_t_minos_trk_primFSLepMnvInitPx, &b_CCDeltaPlusAna_t_minos_trk_primFSLepMnvInitPx);
   fChain->SetBranchAddress("CCDeltaPlusAna_t_minos_trk_primFSLepMnvInitPy", &CCDeltaPlusAna_t_minos_trk_primFSLepMnvInitPy, &b_CCDeltaPlusAna_t_minos_trk_primFSLepMnvInitPy);
   fChain->SetBranchAddress("CCDeltaPlusAna_t_minos_trk_primFSLepMnvInitPz", &CCDeltaPlusAna_t_minos_trk_primFSLepMnvInitPz, &b_CCDeltaPlusAna_t_minos_trk_primFSLepMnvInitPz);
   fChain->SetBranchAddress("CCDeltaPlusAna_t_minos_trk_primFSLepMnvInitX", &CCDeltaPlusAna_t_minos_trk_primFSLepMnvInitX, &b_CCDeltaPlusAna_t_minos_trk_primFSLepMnvInitX);
   fChain->SetBranchAddress("CCDeltaPlusAna_t_minos_trk_primFSLepMnvInitY", &CCDeltaPlusAna_t_minos_trk_primFSLepMnvInitY, &b_CCDeltaPlusAna_t_minos_trk_primFSLepMnvInitY);
   fChain->SetBranchAddress("CCDeltaPlusAna_t_minos_trk_primFSLepMnvInitZ", &CCDeltaPlusAna_t_minos_trk_primFSLepMnvInitZ, &b_CCDeltaPlusAna_t_minos_trk_primFSLepMnvInitZ);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_calc_muon_p", &CCDeltaPlusAna_tammy_calc_muon_p, &b_CCDeltaPlusAna_tammy_calc_muon_p);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_coplanarAngle", &CCDeltaPlusAna_tammy_coplanarAngle, &b_CCDeltaPlusAna_tammy_coplanarAngle);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_endMuonTrajMomentum", &CCDeltaPlusAna_tammy_endMuonTrajMomentum, &b_CCDeltaPlusAna_tammy_endMuonTrajMomentum);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_endMuonTrajXPosition", &CCDeltaPlusAna_tammy_endMuonTrajXPosition, &b_CCDeltaPlusAna_tammy_endMuonTrajXPosition);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_endMuonTrajYPosition", &CCDeltaPlusAna_tammy_endMuonTrajYPosition, &b_CCDeltaPlusAna_tammy_endMuonTrajYPosition);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_endMuonTrajZPosition", &CCDeltaPlusAna_tammy_endMuonTrajZPosition, &b_CCDeltaPlusAna_tammy_endMuonTrajZPosition);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_endProtonTrajMomentum", &CCDeltaPlusAna_tammy_endProtonTrajMomentum, &b_CCDeltaPlusAna_tammy_endProtonTrajMomentum);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_endProtonTrajXPosition", &CCDeltaPlusAna_tammy_endProtonTrajXPosition, &b_CCDeltaPlusAna_tammy_endProtonTrajXPosition);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_endProtonTrajYPosition", &CCDeltaPlusAna_tammy_endProtonTrajYPosition, &b_CCDeltaPlusAna_tammy_endProtonTrajYPosition);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_endProtonTrajZPosition", &CCDeltaPlusAna_tammy_endProtonTrajZPosition, &b_CCDeltaPlusAna_tammy_endProtonTrajZPosition);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_exit_muon_p", &CCDeltaPlusAna_tammy_exit_muon_p, &b_CCDeltaPlusAna_tammy_exit_muon_p);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muon_enu", &CCDeltaPlusAna_tammy_muon_enu, &b_CCDeltaPlusAna_tammy_muon_enu);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muon_odElossMomentum", &CCDeltaPlusAna_tammy_muon_odElossMomentum, &b_CCDeltaPlusAna_tammy_muon_odElossMomentum);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muon_odEndX", &CCDeltaPlusAna_tammy_muon_odEndX, &b_CCDeltaPlusAna_tammy_muon_odEndX);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muon_odEndY", &CCDeltaPlusAna_tammy_muon_odEndY, &b_CCDeltaPlusAna_tammy_muon_odEndY);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muon_odEndZ", &CCDeltaPlusAna_tammy_muon_odEndZ, &b_CCDeltaPlusAna_tammy_muon_odEndZ);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muon_odFaceX", &CCDeltaPlusAna_tammy_muon_odFaceX, &b_CCDeltaPlusAna_tammy_muon_odFaceX);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muon_odFaceY", &CCDeltaPlusAna_tammy_muon_odFaceY, &b_CCDeltaPlusAna_tammy_muon_odFaceY);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muon_odFaceZ", &CCDeltaPlusAna_tammy_muon_odFaceZ, &b_CCDeltaPlusAna_tammy_muon_odFaceZ);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muon_odLastClusZ", &CCDeltaPlusAna_tammy_muon_odLastClusZ, &b_CCDeltaPlusAna_tammy_muon_odLastClusZ);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muon_odStopDistMomentum", &CCDeltaPlusAna_tammy_muon_odStopDistMomentum, &b_CCDeltaPlusAna_tammy_muon_odStopDistMomentum);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muon_odTrackAvgTime", &CCDeltaPlusAna_tammy_muon_odTrackAvgTime, &b_CCDeltaPlusAna_tammy_muon_odTrackAvgTime);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muon_phi", &CCDeltaPlusAna_tammy_muon_phi, &b_CCDeltaPlusAna_tammy_muon_phi);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muon_q2", &CCDeltaPlusAna_tammy_muon_q2, &b_CCDeltaPlusAna_tammy_muon_q2);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muon_score", &CCDeltaPlusAna_tammy_muon_score, &b_CCDeltaPlusAna_tammy_muon_score);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muon_theta", &CCDeltaPlusAna_tammy_muon_theta, &b_CCDeltaPlusAna_tammy_muon_theta);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muon_thetaX", &CCDeltaPlusAna_tammy_muon_thetaX, &b_CCDeltaPlusAna_tammy_muon_thetaX);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muon_thetaY", &CCDeltaPlusAna_tammy_muon_thetaY, &b_CCDeltaPlusAna_tammy_muon_thetaY);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muotammy_n_odClustersAvgTime", &CCDeltaPlusAna_tammy_muotammy_n_odClustersAvgTime, &b_CCDeltaPlusAna_tammy_muotammy_n_odClustersAvgTime);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_open_angle", &CCDeltaPlusAna_tammy_open_angle, &b_CCDeltaPlusAna_tammy_open_angle);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_pion_chi2_ndf", &CCDeltaPlusAna_tammy_pion_chi2_ndf, &b_CCDeltaPlusAna_tammy_pion_chi2_ndf);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_pion_score", &CCDeltaPlusAna_tammy_pion_score, &b_CCDeltaPlusAna_tammy_pion_score);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_pion_score1", &CCDeltaPlusAna_tammy_pion_score1, &b_CCDeltaPlusAna_tammy_pion_score1);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_pion_score2", &CCDeltaPlusAna_tammy_pion_score2, &b_CCDeltaPlusAna_tammy_pion_score2);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_proton_chi2_ndf", &CCDeltaPlusAna_tammy_proton_chi2_ndf, &b_CCDeltaPlusAna_tammy_proton_chi2_ndf);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_proton_enu", &CCDeltaPlusAna_tammy_proton_enu, &b_CCDeltaPlusAna_tammy_proton_enu);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_proton_p_calCorrection", &CCDeltaPlusAna_tammy_proton_p_calCorrection, &b_CCDeltaPlusAna_tammy_proton_p_calCorrection);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_proton_p_visEnergy", &CCDeltaPlusAna_tammy_proton_p_visEnergy, &b_CCDeltaPlusAna_tammy_proton_p_visEnergy);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_proton_phi", &CCDeltaPlusAna_tammy_proton_phi, &b_CCDeltaPlusAna_tammy_proton_phi);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_proton_q2", &CCDeltaPlusAna_tammy_proton_q2, &b_CCDeltaPlusAna_tammy_proton_q2);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_proton_score", &CCDeltaPlusAna_tammy_proton_score, &b_CCDeltaPlusAna_tammy_proton_score);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_proton_score1", &CCDeltaPlusAna_tammy_proton_score1, &b_CCDeltaPlusAna_tammy_proton_score1);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_proton_score2", &CCDeltaPlusAna_tammy_proton_score2, &b_CCDeltaPlusAna_tammy_proton_score2);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_proton_theta", &CCDeltaPlusAna_tammy_proton_theta, &b_CCDeltaPlusAna_tammy_proton_theta);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_proton_thetaX", &CCDeltaPlusAna_tammy_proton_thetaX, &b_CCDeltaPlusAna_tammy_proton_thetaX);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_proton_thetaY", &CCDeltaPlusAna_tammy_proton_thetaY, &b_CCDeltaPlusAna_tammy_proton_thetaY);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_targetZPos", &CCDeltaPlusAna_tammy_targetZPos, &b_CCDeltaPlusAna_tammy_targetZPos);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_trajMuonPhi", &CCDeltaPlusAna_tammy_trajMuonPhi, &b_CCDeltaPlusAna_tammy_trajMuonPhi);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_trajMuonProngMomentum", &CCDeltaPlusAna_tammy_trajMuonProngMomentum, &b_CCDeltaPlusAna_tammy_trajMuonProngMomentum);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_trajMuonTheta", &CCDeltaPlusAna_tammy_trajMuonTheta, &b_CCDeltaPlusAna_tammy_trajMuonTheta);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_trajProtonPhi", &CCDeltaPlusAna_tammy_trajProtonPhi, &b_CCDeltaPlusAna_tammy_trajProtonPhi);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_trajProtonProngMomentum", &CCDeltaPlusAna_tammy_trajProtonProngMomentum, &b_CCDeltaPlusAna_tammy_trajProtonProngMomentum);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_trajProtonTheta", &CCDeltaPlusAna_tammy_trajProtonTheta, &b_CCDeltaPlusAna_tammy_trajProtonTheta);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_endPointVtxPDG", CCDeltaPlusAna_tammy_endPointVtxPDG, &b_CCDeltaPlusAna_tammy_endPointVtxPDG);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_endPointVtxParentId", CCDeltaPlusAna_tammy_endPointVtxParentId, &b_CCDeltaPlusAna_tammy_endPointVtxParentId);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_isolatedPDG", CCDeltaPlusAna_tammy_isolatedPDG, &b_CCDeltaPlusAna_tammy_isolatedPDG);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_isolatedParentId", CCDeltaPlusAna_tammy_isolatedParentId, &b_CCDeltaPlusAna_tammy_isolatedParentId);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muonFuzzPDG", CCDeltaPlusAna_tammy_muonFuzzPDG, &b_CCDeltaPlusAna_tammy_muonFuzzPDG);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muonFuzzParentId", CCDeltaPlusAna_tammy_muonFuzzParentId, &b_CCDeltaPlusAna_tammy_muonFuzzParentId);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_primaryVtxPDG", CCDeltaPlusAna_tammy_primaryVtxPDG, &b_CCDeltaPlusAna_tammy_primaryVtxPDG);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_primaryVtxParentId", CCDeltaPlusAna_tammy_primaryVtxParentId, &b_CCDeltaPlusAna_tammy_primaryVtxParentId);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_protonFuzzPDG", CCDeltaPlusAna_tammy_protonFuzzPDG, &b_CCDeltaPlusAna_tammy_protonFuzzPDG);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_protonFuzzParentId", CCDeltaPlusAna_tammy_protonFuzzParentId, &b_CCDeltaPlusAna_tammy_protonFuzzParentId);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_secondaryVtxPDG", CCDeltaPlusAna_tammy_secondaryVtxPDG, &b_CCDeltaPlusAna_tammy_secondaryVtxPDG);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_secondaryVtxParentId", CCDeltaPlusAna_tammy_secondaryVtxParentId, &b_CCDeltaPlusAna_tammy_secondaryVtxParentId);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_endPointVtxEnergy", CCDeltaPlusAna_tammy_endPointVtxEnergy, &b_CCDeltaPlusAna_tammy_endPointVtxEnergy);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_endPointVtxTrueEnergy", CCDeltaPlusAna_tammy_endPointVtxTrueEnergy, &b_CCDeltaPlusAna_tammy_endPointVtxTrueEnergy);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_isolatedEnergy", CCDeltaPlusAna_tammy_isolatedEnergy, &b_CCDeltaPlusAna_tammy_isolatedEnergy);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_isolatedTrueEnergy", CCDeltaPlusAna_tammy_isolatedTrueEnergy, &b_CCDeltaPlusAna_tammy_isolatedTrueEnergy);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_minos_uv", CCDeltaPlusAna_tammy_minos_uv, &b_CCDeltaPlusAna_tammy_minos_uv);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muonFuzzEnergy", CCDeltaPlusAna_tammy_muonFuzzEnergy, &b_CCDeltaPlusAna_tammy_muonFuzzEnergy);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muonFuzzTrueEnergy", CCDeltaPlusAna_tammy_muonFuzzTrueEnergy, &b_CCDeltaPlusAna_tammy_muonFuzzTrueEnergy);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muon_chi2_values", CCDeltaPlusAna_tammy_muon_chi2_values, &b_CCDeltaPlusAna_tammy_muon_chi2_values);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muon_endPoint", CCDeltaPlusAna_tammy_muon_endPoint, &b_CCDeltaPlusAna_tammy_muon_endPoint);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muon_p_values", CCDeltaPlusAna_tammy_muon_p_values, &b_CCDeltaPlusAna_tammy_muon_p_values);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_muon_startPoint", CCDeltaPlusAna_tammy_muon_startPoint, &b_CCDeltaPlusAna_tammy_muon_startPoint);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_primaryVtxEnergy", CCDeltaPlusAna_tammy_primaryVtxEnergy, &b_CCDeltaPlusAna_tammy_primaryVtxEnergy);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_primaryVtxTrueEnergy", CCDeltaPlusAna_tammy_primaryVtxTrueEnergy, &b_CCDeltaPlusAna_tammy_primaryVtxTrueEnergy);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_protonFuzzEnergy", CCDeltaPlusAna_tammy_protonFuzzEnergy, &b_CCDeltaPlusAna_tammy_protonFuzzEnergy);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_protonFuzzTrueEnergy", CCDeltaPlusAna_tammy_protonFuzzTrueEnergy, &b_CCDeltaPlusAna_tammy_protonFuzzTrueEnergy);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_proton_4p", CCDeltaPlusAna_tammy_proton_4p, &b_CCDeltaPlusAna_tammy_proton_4p);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_proton_chi2_values", CCDeltaPlusAna_tammy_proton_chi2_values, &b_CCDeltaPlusAna_tammy_proton_chi2_values);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_proton_endPoint", CCDeltaPlusAna_tammy_proton_endPoint, &b_CCDeltaPlusAna_tammy_proton_endPoint);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_proton_p_values", CCDeltaPlusAna_tammy_proton_p_values, &b_CCDeltaPlusAna_tammy_proton_p_values);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_proton_startPoint", CCDeltaPlusAna_tammy_proton_startPoint, &b_CCDeltaPlusAna_tammy_proton_startPoint);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_secondaryVtxEnergy", CCDeltaPlusAna_tammy_secondaryVtxEnergy, &b_CCDeltaPlusAna_tammy_secondaryVtxEnergy);
   fChain->SetBranchAddress("CCDeltaPlusAna_tammy_secondaryVtxTrueEnergy", CCDeltaPlusAna_tammy_secondaryVtxTrueEnergy, &b_CCDeltaPlusAna_tammy_secondaryVtxTrueEnergy);
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

Bool_t CCDeltaPlus::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

    cout<<"Done!"<<endl;
   return kTRUE;
}

void CCDeltaPlus::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

Int_t CCDeltaPlus::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef CCDeltaPlus_cxx

