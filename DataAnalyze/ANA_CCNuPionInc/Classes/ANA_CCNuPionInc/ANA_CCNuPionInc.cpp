/*
    See ANA_CCNuPionInc.h header for Class Information
*/

#define ANA_CCNuPionInc_cxx

#include "ANA_CCNuPionInc.h"
#include "Cuts.cpp"


ANA_CCNuPionInc::ANA_CCNuPionInc()
{



}

void ANA_CCNuPionInc::run(string playlist, string rootFileName)
{

   
 
    //------------------------------------------------------------------------
    // Open files for writing
    //------------------------------------------------------------------------

    openFiles();

    TFile* f = new TFile(rootFileName.c_str(),"RECREATE");

    //------------------------------------------------------------------------
    // Create chain
    //------------------------------------------------------------------------

    TChain* fChain = new TChain("CCNuPionInc");
    Init(playlist, fChain);

    if (!fChain) return;
    if (fChain == 0) return;
    
    //------------------------------------------------------------------------
    // Branch Selection for Performance
    //------------------------------------------------------------------------
//     fChain->SetBranchStatus("*",0);  // disable all branches
//     
//     // Analysis Variables
//     fChain->SetBranchStatus("ev_run",1);  // activate
//     fChain->SetBranchStatus("ev_subrun",1);  // activate
//     fChain->SetBranchStatus("ev_gate",1);  // activate
//     
//     // Cut Variables
//     fChain->SetBranchStatus("mc_vtx",1);
//     
//     // Incoming Particle
//     fChain->SetBranchStatus("mc_incomingPartVec",1);  // activate
//     
//     // Final State Particles
//     fChain->SetBranchStatus("mc_FSPartPx",1);  // activate
//     fChain->SetBranchStatus("mc_FSPartPy",1);  // activate
//     fChain->SetBranchStatus("mc_FSPartPz",1);  // activate
//     fChain->SetBranchStatus("mc_FSPartE",1);  // activate
//     fChain->SetBranchStatus("mc_FSPartPDG",1);
//     
//     
//     
//     // Reconstruction Variables
//     fChain->SetBranchStatus("CCNuPionInc_muon_px",1);  // activate
//     fChain->SetBranchStatus("CCNuPionInc_muon_py",1);  // activate
//     fChain->SetBranchStatus("CCNuPionInc_muon_pz",1);  // activate
//     fChain->SetBranchStatus("CCNuPionInc_muon_E",1);  // activate


    //------------------------------------------------------------------------
    // Initialize the Analysis Variables and Histograms
    //------------------------------------------------------------------------

    initVariables();
    initHistograms();
    
    //------------------------------------------------------------------------
    // Declare Other Variables
    //------------------------------------------------------------------------



    //------------------------------------------------------------------------
    // Loop over Chain
    //------------------------------------------------------------------------
    Long64_t nbytes = 0, nb = 0;
    
    Long64_t nentries = fChain->GetEntriesFast();
    cout<<"There are "<<nentries<<" entries!"<<endl;

    for (Long64_t jentry=0; jentry<nentries;jentry++) {
    
        nb = fChain->GetEntry(jentry);   nbytes += nb;    
        Long64_t ientry = fChain->GetEntry(jentry);
        if (ientry == 0) {
            cout<<"GetEntry failure "<<jentry<<endl;
            break;
        }
    
        // Progress Message on Terminal
        if (jentry%5000 == 0){
            cout<<" Entry "<<jentry<<endl;
        }
    
    
        //----------------------------------------------------------------------
        // Apply Cuts
        //----------------------------------------------------------------------
        
        // Count All Events before Cuts
        nCutList->nAll->increment();
        
        // Volume Cut
        if( !isVertexContained()){
            continue;
        }
        nCutList->nVolume->increment();
        
        // Muon Cut
        muon->ind = findParticle(PDG_List::mu_minus);
        if(muon->ind == -1){
            continue;
        }
        nCutList->nMuon->increment();
        
        // Proton Cut
        proton->ind = findProton();
        if(proton->ind == -1){
            continue;
        }
        nCutList->nProton->increment();
        
        // Pion Cut
        pion->ind = findPion();
        if(pion->ind == -1){
            continue;
        }
        nCutList->nPion->increment();
        

        //----------------------------------------------------------------------
        // Fill Particles
        //----------------------------------------------------------------------
    
        fillParticleTrue(muon);
        fillParticleTrue(proton);
        fillParticleTrue(pion);
        
        //----------------------------------------------------------------------
        // Fill Histograms
        //----------------------------------------------------------------------
       
        fillHistograms();
        
      

    }    
    
    // Write the Root File
    cout<<">> Writing "<<rootFileName<<endl;
    f->Write();
    
    nCutList->writeCutTable();
    
    
    closeFiles();

}

// -------------------------------------------------------------------------
//     Specific Functions
//--------------------------------------------------------------------------

void ANA_CCNuPionInc::fillHistograms()
{
    P_muon->Fill(muon->p4[1].P());
    P_proton->Fill(proton->p4[1].P());
    P_pion->Fill(pion->p4[1].P());
    
    Angle_muon->Fill(muon->angleBeam[1] * TMath::RadToDeg());
    Angle_proton->Fill(proton->angleBeam[1] * TMath::RadToDeg());
    Angle_pion->Fill(pion->angleBeam[1] * TMath::RadToDeg());

}

void ANA_CCNuPionInc::fillParticleTrue(Particle* part)
{
    int ind = part->ind;
    double angleBeam;
    
    // Fill 4-Momentum
    part->p4[1].SetPxPyPzE(mc_FSPartPx[ind],mc_FSPartPy[ind],mc_FSPartPz[ind],mc_FSPartE[ind]);
       
    // Calculate and Fill Angle wrt Beam
    angleBeam = part->p4[1].Angle(*beam_p3);
    part->angleBeam[1] = angleBeam;
    
}


void ANA_CCNuPionInc::initHistograms()
{
    P_muon = new TH1F( "P_muon","True Muon Momentum",binList->muonP->get_nBins(), binList->muonP->get_min(), binList->muonP->get_max() );
    P_muon->GetXaxis()->SetTitle("True Muon Momentum MeV");
    P_muon->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList->muonP->get_width()));
    
    P_proton = new TH1F( "P_proton","True Proton Momentum",binList->protonP->get_nBins(), binList->protonP->get_min(), binList->protonP->get_max() );
    P_proton->GetXaxis()->SetTitle("True Proton Momentum MeV");
    P_proton->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList->protonP->get_width()));
    
    P_pion = new TH1F( "P_pion","True Pion Momentum",binList->pionP->get_nBins(), binList->pionP->get_min(), binList->pionP->get_max() );
    P_pion->GetXaxis()->SetTitle("True Pion Momentum MeV");
    P_pion->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList->pionP->get_width()));
    
    Angle_muon = new TH1F( "Angle_muon","Angle: Beam vs Muon",binList->angle->get_nBins(), binList->angle->get_min(), binList->angle->get_max() );
    Angle_muon->GetXaxis()->SetTitle("Angle");
    Angle_muon->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList->angle->get_width()));
    
    Angle_proton = new TH1F( "Angle_proton","Angle: Beam vs Proton",binList->angle->get_nBins(), binList->angle->get_min(), binList->angle->get_max() );
    Angle_proton->GetXaxis()->SetTitle("Angle");
    Angle_proton->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList->angle->get_width()));
    
    Angle_pion = new TH1F( "Angle_pion","Angle: Beam vs Pion",binList->angle->get_nBins(), binList->angle->get_min(), binList->angle->get_max() );
    Angle_pion->GetXaxis()->SetTitle("Angle");
    Angle_pion->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList->angle->get_width()));
    
    
    cout<<"Histograms are Initialized!"<<endl;
}

void ANA_CCNuPionInc::initSingleHistogram(TH1F* hist, string histName, string title, string xLabel, string yLabel, SingleBin* bin)
{
    // Reserved for a future version - Still needs testing
    
    cout<<bin->get_nBins()<<endl;
    hist = new TH1F(histName.c_str(),title.c_str(),bin->get_nBins(), bin->get_min(), bin->get_max());
    hist->GetXaxis()->SetTitle(xLabel.c_str());
    hist->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin->get_width()));
}




void ANA_CCNuPionInc::initVariables()
{

    // -------------------------------------------------------------------------
    //     Memory Allocation
    //--------------------------------------------------------------------------
    // Allocate Memory
    beam_p3 = new TVector3;
    muon = new Particle;
    proton = new Particle;
    pion = new Particle;
    
    nCutList = new CutNumberList;
    binList = new BinList;
    
    // -------------------------------------------------------------------------
    //     Initialization
    //--------------------------------------------------------------------------
    // Default Beam Configuration
    beam_p3->SetXYZ(1.0,1.0,1.0);
    beam_p3->SetPhi(-1.554);
    beam_p3->SetTheta(0.059);
    
    

    cout<<"Variables are Initialized!"<<endl;

}

void ANA_CCNuPionInc::closeFiles()
{
    readme.close();
}

void ANA_CCNuPionInc::openFiles()
{
    readmeFile = getFileLocation(Folder_List::OUT, Folder_List::TEXT, Folder_List::F_TEXT_README);

    readme.open( readmeFile.c_str() );
    if( !readme.is_open() ){
        cerr<<"Cannot open output text file: "<<readmeFile<<endl;
        exit(1);
    }

    writeReadme();
}

void ANA_CCNuPionInc::writeReadme()
{

readme<<"Test"<<endl;

}

// -------------------------------------------------------------------------
//     Default Functions
//--------------------------------------------------------------------------

#ifdef ANA_CCNuPionInc_cxx
ANA_CCNuPionInc::ANA_CCNuPionInc()
{
    
}

ANA_CCNuPionInc::~ANA_CCNuPionInc()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ANA_CCNuPionInc::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t ANA_CCNuPionInc::LoadTree(Long64_t entry)
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

void ANA_CCNuPionInc::Init(string playlist, TChain* fChain)
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
        cout<<"Playlist: "<<playlist.c_str()<<endl;
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
   fChain->SetBranchAddress("broken_track_most_us_plane", &broken_track_most_us_plane, &b_broken_track_most_us_plane);
   fChain->SetBranchAddress("ddead", &ddead, &b_ddead);
   fChain->SetBranchAddress("dead", &dead, &b_dead);
   fChain->SetBranchAddress("has_vtx_michel", &has_vtx_michel, &b_has_vtx_michel);
   fChain->SetBranchAddress("n_anchored_long_trk_prongs", &n_anchored_long_trk_prongs, &b_n_anchored_long_trk_prongs);
   fChain->SetBranchAddress("n_anchored_short_trk_prongs", &n_anchored_short_trk_prongs, &b_n_anchored_short_trk_prongs);
   fChain->SetBranchAddress("n_dsp_blob_prongs", &n_dsp_blob_prongs, &b_n_dsp_blob_prongs);
   fChain->SetBranchAddress("n_iso_blob_prongs", &n_iso_blob_prongs, &b_n_iso_blob_prongs);
   fChain->SetBranchAddress("n_iso_trk_prongs", &n_iso_trk_prongs, &b_n_iso_trk_prongs);
   fChain->SetBranchAddress("n_long_tracks", &n_long_tracks, &b_n_long_tracks);
   fChain->SetBranchAddress("n_sepTrk_permutations", &n_sepTrk_permutations, &b_n_sepTrk_permutations);
   fChain->SetBranchAddress("n_short_tracks", &n_short_tracks, &b_n_short_tracks);
   fChain->SetBranchAddress("n_startpoint_vertices", &n_startpoint_vertices, &b_n_startpoint_vertices);
   fChain->SetBranchAddress("n_twoTrk_permutations", &n_twoTrk_permutations, &b_n_twoTrk_permutations);
   fChain->SetBranchAddress("n_us_muon_clusters", &n_us_muon_clusters, &b_n_us_muon_clusters);
   fChain->SetBranchAddress("n_vtx_michel_views", &n_vtx_michel_views, &b_n_vtx_michel_views);
   fChain->SetBranchAddress("n_vtx_prongs", &n_vtx_prongs, &b_n_vtx_prongs);
   fChain->SetBranchAddress("phys_energy_in_road_downstream_nplanes", &phys_energy_in_road_downstream_nplanes, &b_phys_energy_in_road_downstream_nplanes);
   fChain->SetBranchAddress("phys_energy_in_road_upstream_nplanes", &phys_energy_in_road_upstream_nplanes, &b_phys_energy_in_road_upstream_nplanes);
   fChain->SetBranchAddress("phys_n_dead_discr_pair", &phys_n_dead_discr_pair, &b_phys_n_dead_discr_pair);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_in_prim_track_region", &phys_n_dead_discr_pair_in_prim_track_region, &b_phys_n_dead_discr_pair_in_prim_track_region);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_two_mod_downstream_prim_track", &phys_n_dead_discr_pair_two_mod_downstream_prim_track, &b_phys_n_dead_discr_pair_two_mod_downstream_prim_track);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_two_mod_upstream_prim_vtx", &phys_n_dead_discr_pair_two_mod_upstream_prim_vtx, &b_phys_n_dead_discr_pair_two_mod_upstream_prim_vtx);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_upstream_prim_track_proj", &phys_n_dead_discr_pair_upstream_prim_track_proj, &b_phys_n_dead_discr_pair_upstream_prim_track_proj);
   fChain->SetBranchAddress("phys_vertex_is_fiducial", &phys_vertex_is_fiducial, &b_phys_vertex_is_fiducial);
   fChain->SetBranchAddress("tdead", &tdead, &b_tdead);
   fChain->SetBranchAddress("twoTrk_road_planes", &twoTrk_road_planes, &b_twoTrk_road_planes);
   fChain->SetBranchAddress("udead", &udead, &b_udead);
   fChain->SetBranchAddress("upstream_plane_num", &upstream_plane_num, &b_upstream_plane_num);
   fChain->SetBranchAddress("dispersedExtraE", &dispersedExtraE, &b_dispersedExtraE);
   fChain->SetBranchAddress("energy_from_mc", &energy_from_mc, &b_energy_from_mc);
   fChain->SetBranchAddress("energy_from_mc_fraction", &energy_from_mc_fraction, &b_energy_from_mc_fraction);
   fChain->SetBranchAddress("energy_from_mc_fraction_of_highest", &energy_from_mc_fraction_of_highest, &b_energy_from_mc_fraction_of_highest);
   fChain->SetBranchAddress("hadronVisibleE", &hadronVisibleE, &b_hadronVisibleE);
   fChain->SetBranchAddress("muonVisibleE", &muonVisibleE, &b_muonVisibleE);
   fChain->SetBranchAddress("muon_phi", &muon_phi, &b_muon_phi);
   fChain->SetBranchAddress("muon_theta", &muon_theta, &b_muon_theta);
   fChain->SetBranchAddress("muon_thetaX", &muon_thetaX, &b_muon_thetaX);
   fChain->SetBranchAddress("muon_thetaY", &muon_thetaY, &b_muon_thetaY);
   fChain->SetBranchAddress("phys_energy_dispersed", &phys_energy_dispersed, &b_phys_energy_dispersed);
   fChain->SetBranchAddress("phys_energy_in_road_downstream", &phys_energy_in_road_downstream, &b_phys_energy_in_road_downstream);
   fChain->SetBranchAddress("phys_energy_in_road_upstream", &phys_energy_in_road_upstream, &b_phys_energy_in_road_upstream);
   fChain->SetBranchAddress("phys_energy_unattached", &phys_energy_unattached, &b_phys_energy_unattached);
   fChain->SetBranchAddress("prim_vtx_smallest_opening_angle", &prim_vtx_smallest_opening_angle, &b_prim_vtx_smallest_opening_angle);
   fChain->SetBranchAddress("time", &time, &b_time);
   fChain->SetBranchAddress("totalIDVisibleE", &totalIDVisibleE, &b_totalIDVisibleE);
   fChain->SetBranchAddress("totalODVisibleE", &totalODVisibleE, &b_totalODVisibleE);
   fChain->SetBranchAddress("totalVisibleE", &totalVisibleE, &b_totalVisibleE);
   fChain->SetBranchAddress("twoTrk_energy_in_road", &twoTrk_energy_in_road, &b_twoTrk_energy_in_road);
   fChain->SetBranchAddress("unattachedExtraE", &unattachedExtraE, &b_unattachedExtraE);
   fChain->SetBranchAddress("vtxBlobExtraE", &vtxBlobExtraE, &b_vtxBlobExtraE);
   fChain->SetBranchAddress("vtx_michel_distance", &vtx_michel_distance, &b_vtx_michel_distance);
   fChain->SetBranchAddress("well_fit_vertex_angle", &well_fit_vertex_angle, &b_well_fit_vertex_angle);
   fChain->SetBranchAddress("sepTrk_plane", sepTrk_plane, &b_sepTrk_plane);
   fChain->SetBranchAddress("twoTrk_plane", twoTrk_plane, &b_twoTrk_plane);
   fChain->SetBranchAddress("sepTrk_midPtDiff", sepTrk_midPtDiff, &b_sepTrk_midPtDiff);
   fChain->SetBranchAddress("sepTrk_opening_angle", sepTrk_opening_angle, &b_sepTrk_opening_angle);
   fChain->SetBranchAddress("twoTrk_opening_angle", twoTrk_opening_angle, &b_twoTrk_opening_angle);
   fChain->SetBranchAddress("truth_has_physics_event", &truth_has_physics_event, &b_truth_has_physics_event);
   fChain->SetBranchAddress("truth_isNuSignal", &truth_isNuSignal, &b_truth_isNuSignal);
   fChain->SetBranchAddress("truth_isNuBarSignal", &truth_isNuBarSignal, &b_truth_isNuBarSignal);
   fChain->SetBranchAddress("truth_isTrackablePion", &truth_isTrackablePion, &b_truth_isTrackablePion);
   fChain->SetBranchAddress("truth_isFidVol", &truth_isFidVol, &b_truth_isFidVol);
   fChain->SetBranchAddress("truth_isPlausible", &truth_isPlausible, &b_truth_isPlausible);
   fChain->SetBranchAddress("truth_reco_hasGoodObjects", &truth_reco_hasGoodObjects, &b_truth_reco_hasGoodObjects);
   fChain->SetBranchAddress("truth_reco_isGoodVertex", &truth_reco_isGoodVertex, &b_truth_reco_isGoodVertex);
   fChain->SetBranchAddress("truth_reco_isWellFitVertex", &truth_reco_isWellFitVertex, &b_truth_reco_isWellFitVertex);
   fChain->SetBranchAddress("truth_reco_isFidVol", &truth_reco_isFidVol, &b_truth_reco_isFidVol);
   fChain->SetBranchAddress("truth_reco_isFidVol_smeared", &truth_reco_isFidVol_smeared, &b_truth_reco_isFidVol_smeared);
   fChain->SetBranchAddress("truth_reco_isMinosMatch", &truth_reco_isMinosMatch, &b_truth_reco_isMinosMatch);
   fChain->SetBranchAddress("truth_reco_isBrokenTrack", &truth_reco_isBrokenTrack, &b_truth_reco_isBrokenTrack);
   fChain->SetBranchAddress("truth_N_gamma", &truth_N_gamma, &b_truth_N_gamma);
   fChain->SetBranchAddress("truth_N_mum", &truth_N_mum, &b_truth_N_mum);
   fChain->SetBranchAddress("truth_N_mup", &truth_N_mup, &b_truth_N_mup);
   fChain->SetBranchAddress("truth_N_neu", &truth_N_neu, &b_truth_N_neu);
   fChain->SetBranchAddress("truth_N_other", &truth_N_other, &b_truth_N_other);
   fChain->SetBranchAddress("truth_N_pbar", &truth_N_pbar, &b_truth_N_pbar);
   fChain->SetBranchAddress("truth_N_pi0", &truth_N_pi0, &b_truth_N_pi0);
   fChain->SetBranchAddress("truth_N_pim", &truth_N_pim, &b_truth_N_pim);
   fChain->SetBranchAddress("truth_N_pip", &truth_N_pip, &b_truth_N_pip);
   fChain->SetBranchAddress("truth_N_pro", &truth_N_pro, &b_truth_N_pro);
   fChain->SetBranchAddress("truth_Nd_gamma", &truth_Nd_gamma, &b_truth_Nd_gamma);
   fChain->SetBranchAddress("truth_Nd_mum", &truth_Nd_mum, &b_truth_Nd_mum);
   fChain->SetBranchAddress("truth_Nd_mup", &truth_Nd_mup, &b_truth_Nd_mup);
   fChain->SetBranchAddress("truth_Nd_neu", &truth_Nd_neu, &b_truth_Nd_neu);
   fChain->SetBranchAddress("truth_Nd_pbar", &truth_Nd_pbar, &b_truth_Nd_pbar);
   fChain->SetBranchAddress("truth_Nd_pi0", &truth_Nd_pi0, &b_truth_Nd_pi0);
   fChain->SetBranchAddress("truth_Nd_pim", &truth_Nd_pim, &b_truth_Nd_pim);
   fChain->SetBranchAddress("truth_Nd_pip", &truth_Nd_pip, &b_truth_Nd_pip);
   fChain->SetBranchAddress("truth_Nd_primaries", &truth_Nd_primaries, &b_truth_Nd_primaries);
   fChain->SetBranchAddress("truth_Nd_pro", &truth_Nd_pro, &b_truth_Nd_pro);
   fChain->SetBranchAddress("truth_mu_charge", &truth_mu_charge, &b_truth_mu_charge);
   fChain->SetBranchAddress("truth_mu_primaryPlanes", &truth_mu_primaryPlanes, &b_truth_mu_primaryPlanes);
   fChain->SetBranchAddress("truth_mu_totPlanes", &truth_mu_totPlanes, &b_truth_mu_totPlanes);
   fChain->SetBranchAddress("truth_mu_trackID", &truth_mu_trackID, &b_truth_mu_trackID);
   fChain->SetBranchAddress("truth_reco_muonCharge", &truth_reco_muonCharge, &b_truth_reco_muonCharge);
   fChain->SetBranchAddress("truth_target_material", &truth_target_material, &b_truth_target_material);
   fChain->SetBranchAddress("truth_vertex_module", &truth_vertex_module, &b_truth_vertex_module);
   fChain->SetBranchAddress("truth_vertex_plane", &truth_vertex_plane, &b_truth_vertex_plane);
   fChain->SetBranchAddress("truth_mu_E", &truth_mu_E, &b_truth_mu_E);
   fChain->SetBranchAddress("truth_mu_px", &truth_mu_px, &b_truth_mu_px);
   fChain->SetBranchAddress("truth_mu_py", &truth_mu_py, &b_truth_mu_py);
   fChain->SetBranchAddress("truth_mu_pz", &truth_mu_pz, &b_truth_mu_pz);
   fChain->SetBranchAddress("truth_mu_theta_wrtbeam", &truth_mu_theta_wrtbeam, &b_truth_mu_theta_wrtbeam);
   fChain->SetBranchAddress("truth_pi_charge", truth_pi_charge, &b_truth_pi_charge);
   fChain->SetBranchAddress("truth_pi_primaryPlanes", truth_pi_primaryPlanes, &b_truth_pi_primaryPlanes);
   fChain->SetBranchAddress("truth_pi_totPlanes", truth_pi_totPlanes, &b_truth_pi_totPlanes);
   fChain->SetBranchAddress("truth_pi_trackID", truth_pi_trackID, &b_truth_pi_trackID);
   fChain->SetBranchAddress("truth_pro_charge", truth_pro_charge, &b_truth_pro_charge);
   fChain->SetBranchAddress("truth_pro_primaryPlanes", truth_pro_primaryPlanes, &b_truth_pro_primaryPlanes);
   fChain->SetBranchAddress("truth_pro_totPlanes", truth_pro_totPlanes, &b_truth_pro_totPlanes);
   fChain->SetBranchAddress("truth_pro_trackID", truth_pro_trackID, &b_truth_pro_trackID);
   fChain->SetBranchAddress("genie_wgt_n_shifts", &genie_wgt_n_shifts, &b_genie_wgt_n_shifts);
   fChain->SetBranchAddress("truth_genie_wgt_AGKYxF1pi", truth_genie_wgt_AGKYxF1pi, &b_truth_genie_wgt_AGKYxF1pi);
   fChain->SetBranchAddress("truth_genie_wgt_AhtBY", truth_genie_wgt_AhtBY, &b_truth_genie_wgt_AhtBY);
   fChain->SetBranchAddress("truth_genie_wgt_BhtBY", truth_genie_wgt_BhtBY, &b_truth_genie_wgt_BhtBY);
   fChain->SetBranchAddress("truth_genie_wgt_CCQEPauliSupViaFK", truth_genie_wgt_CCQEPauliSupViaFK, &b_truth_genie_wgt_CCQEPauliSupViaFK);
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
   fChain->SetBranchAddress("truth_pi_E", truth_pi_E, &b_truth_pi_E);
   fChain->SetBranchAddress("truth_pi_px", truth_pi_px, &b_truth_pi_px);
   fChain->SetBranchAddress("truth_pi_py", truth_pi_py, &b_truth_pi_py);
   fChain->SetBranchAddress("truth_pi_pz", truth_pi_pz, &b_truth_pi_pz);
   fChain->SetBranchAddress("truth_pi_theta_wrtbeam", truth_pi_theta_wrtbeam, &b_truth_pi_theta_wrtbeam);
   fChain->SetBranchAddress("truth_pro_E", truth_pro_E, &b_truth_pro_E);
   fChain->SetBranchAddress("truth_pro_px", truth_pro_px, &b_truth_pro_px);
   fChain->SetBranchAddress("truth_pro_py", truth_pro_py, &b_truth_pro_py);
   fChain->SetBranchAddress("truth_pro_pz", truth_pro_pz, &b_truth_pro_pz);
   fChain->SetBranchAddress("truth_pro_theta_wrtbeam", truth_pro_theta_wrtbeam, &b_truth_pro_theta_wrtbeam);
   fChain->SetBranchAddress("CCNuPionInc_nuFlavor", &CCNuPionInc_nuFlavor, &b_CCNuPionInc_nuFlavor);
   fChain->SetBranchAddress("CCNuPionInc_nuHelicity", &CCNuPionInc_nuHelicity, &b_CCNuPionInc_nuHelicity);
   fChain->SetBranchAddress("CCNuPionInc_intCurrent", &CCNuPionInc_intCurrent, &b_CCNuPionInc_intCurrent);
   fChain->SetBranchAddress("CCNuPionInc_intType", &CCNuPionInc_intType, &b_CCNuPionInc_intType);
   fChain->SetBranchAddress("CCNuPionInc_E", &CCNuPionInc_E, &b_CCNuPionInc_E);
   fChain->SetBranchAddress("CCNuPionInc_Q2", &CCNuPionInc_Q2, &b_CCNuPionInc_Q2);
   fChain->SetBranchAddress("CCNuPionInc_x", &CCNuPionInc_x, &b_CCNuPionInc_x);
   fChain->SetBranchAddress("CCNuPionInc_y", &CCNuPionInc_y, &b_CCNuPionInc_y);
   fChain->SetBranchAddress("CCNuPionInc_W", &CCNuPionInc_W, &b_CCNuPionInc_W);
   fChain->SetBranchAddress("CCNuPionInc_score", &CCNuPionInc_score, &b_CCNuPionInc_score);
   fChain->SetBranchAddress("CCNuPionInc_leptonE", CCNuPionInc_leptonE, &b_CCNuPionInc_leptonE);
   fChain->SetBranchAddress("CCNuPionInc_vtx", CCNuPionInc_vtx, &b_CCNuPionInc_vtx);
   fChain->SetBranchAddress("CCNuPionInc_minos_trk_is_contained", &CCNuPionInc_minos_trk_is_contained, &b_CCNuPionInc_minos_trk_is_contained);
   fChain->SetBranchAddress("CCNuPionInc_minos_trk_is_ok", &CCNuPionInc_minos_trk_is_ok, &b_CCNuPionInc_minos_trk_is_ok);
   fChain->SetBranchAddress("CCNuPionInc_minos_used_range", &CCNuPionInc_minos_used_range, &b_CCNuPionInc_minos_used_range);
   fChain->SetBranchAddress("CCNuPionInc_minos_used_curvature", &CCNuPionInc_minos_used_curvature, &b_CCNuPionInc_minos_used_curvature);
   fChain->SetBranchAddress("CCNuPionInc_hadron_number", &CCNuPionInc_hadron_number, &b_CCNuPionInc_hadron_number);
   fChain->SetBranchAddress("CCNuPionInc_minos_trk_end_plane", &CCNuPionInc_minos_trk_end_plane, &b_CCNuPionInc_minos_trk_end_plane);
   fChain->SetBranchAddress("CCNuPionInc_minos_trk_quality", &CCNuPionInc_minos_trk_quality, &b_CCNuPionInc_minos_trk_quality);
   fChain->SetBranchAddress("CCNuPionInc_muon_N_minosTracks", &CCNuPionInc_muon_N_minosTracks, &b_CCNuPionInc_muon_N_minosTracks);
   fChain->SetBranchAddress("CCNuPionInc_muon_minervaTrack_types", &CCNuPionInc_muon_minervaTrack_types, &b_CCNuPionInc_muon_minervaTrack_types);
   fChain->SetBranchAddress("CCNuPionInc_muon_minosTrackQuality", &CCNuPionInc_muon_minosTrackQuality, &b_CCNuPionInc_muon_minosTrackQuality);
   fChain->SetBranchAddress("CCNuPionInc_muon_roadUpstreamPlanes", &CCNuPionInc_muon_roadUpstreamPlanes, &b_CCNuPionInc_muon_roadUpstreamPlanes);
   fChain->SetBranchAddress("CCNuPionInc_r_minos_trk_vtx_plane", &CCNuPionInc_r_minos_trk_vtx_plane, &b_CCNuPionInc_r_minos_trk_vtx_plane);
   fChain->SetBranchAddress("CCNuPionInc_t_minos_trk_numFSMuons", &CCNuPionInc_t_minos_trk_numFSMuons, &b_CCNuPionInc_t_minos_trk_numFSMuons);
   fChain->SetBranchAddress("CCNuPionInc_t_minos_trk_primFSLeptonPDG", &CCNuPionInc_t_minos_trk_primFSLeptonPDG, &b_CCNuPionInc_t_minos_trk_primFSLeptonPDG);
   fChain->SetBranchAddress("CCNuPionInc_vtx_module", &CCNuPionInc_vtx_module, &b_CCNuPionInc_vtx_module);
   fChain->SetBranchAddress("CCNuPionInc_vtx_plane", &CCNuPionInc_vtx_plane, &b_CCNuPionInc_vtx_plane);
   fChain->SetBranchAddress("CCNuPionInc_hadron_recoil", &CCNuPionInc_hadron_recoil, &b_CCNuPionInc_hadron_recoil);
   fChain->SetBranchAddress("CCNuPionInc_hadron_recoil_CCInc", &CCNuPionInc_hadron_recoil_CCInc, &b_CCNuPionInc_hadron_recoil_CCInc);
   fChain->SetBranchAddress("CCNuPionInc_hadron_recoil_CCInc_em", &CCNuPionInc_hadron_recoil_CCInc_em, &b_CCNuPionInc_hadron_recoil_CCInc_em);
   fChain->SetBranchAddress("CCNuPionInc_hadron_recoil_CCInc_kaon", &CCNuPionInc_hadron_recoil_CCInc_kaon, &b_CCNuPionInc_hadron_recoil_CCInc_kaon);
   fChain->SetBranchAddress("CCNuPionInc_hadron_recoil_CCInc_neutron", &CCNuPionInc_hadron_recoil_CCInc_neutron, &b_CCNuPionInc_hadron_recoil_CCInc_neutron);
   fChain->SetBranchAddress("CCNuPionInc_hadron_recoil_CCInc_other", &CCNuPionInc_hadron_recoil_CCInc_other, &b_CCNuPionInc_hadron_recoil_CCInc_other);
   fChain->SetBranchAddress("CCNuPionInc_hadron_recoil_CCInc_pion", &CCNuPionInc_hadron_recoil_CCInc_pion, &b_CCNuPionInc_hadron_recoil_CCInc_pion);
   fChain->SetBranchAddress("CCNuPionInc_hadron_recoil_CCInc_proton", &CCNuPionInc_hadron_recoil_CCInc_proton, &b_CCNuPionInc_hadron_recoil_CCInc_proton);
   fChain->SetBranchAddress("CCNuPionInc_hadron_recoil_CCInc_xtalk", &CCNuPionInc_hadron_recoil_CCInc_xtalk, &b_CCNuPionInc_hadron_recoil_CCInc_xtalk);
   fChain->SetBranchAddress("CCNuPionInc_minos_trk_bave", &CCNuPionInc_minos_trk_bave, &b_CCNuPionInc_minos_trk_bave);
   fChain->SetBranchAddress("CCNuPionInc_minos_trk_chi2", &CCNuPionInc_minos_trk_chi2, &b_CCNuPionInc_minos_trk_chi2);
   fChain->SetBranchAddress("CCNuPionInc_minos_trk_end_u", &CCNuPionInc_minos_trk_end_u, &b_CCNuPionInc_minos_trk_end_u);
   fChain->SetBranchAddress("CCNuPionInc_minos_trk_end_v", &CCNuPionInc_minos_trk_end_v, &b_CCNuPionInc_minos_trk_end_v);
   fChain->SetBranchAddress("CCNuPionInc_minos_trk_end_x", &CCNuPionInc_minos_trk_end_x, &b_CCNuPionInc_minos_trk_end_x);
   fChain->SetBranchAddress("CCNuPionInc_minos_trk_end_y", &CCNuPionInc_minos_trk_end_y, &b_CCNuPionInc_minos_trk_end_y);
   fChain->SetBranchAddress("CCNuPionInc_minos_trk_end_z", &CCNuPionInc_minos_trk_end_z, &b_CCNuPionInc_minos_trk_end_z);
   fChain->SetBranchAddress("CCNuPionInc_minos_trk_eqp", &CCNuPionInc_minos_trk_eqp, &b_CCNuPionInc_minos_trk_eqp);
   fChain->SetBranchAddress("CCNuPionInc_minos_trk_eqp_qp", &CCNuPionInc_minos_trk_eqp_qp, &b_CCNuPionInc_minos_trk_eqp_qp);
   fChain->SetBranchAddress("CCNuPionInc_minos_trk_fit_pass", &CCNuPionInc_minos_trk_fit_pass, &b_CCNuPionInc_minos_trk_fit_pass);
   fChain->SetBranchAddress("CCNuPionInc_minos_trk_ndf", &CCNuPionInc_minos_trk_ndf, &b_CCNuPionInc_minos_trk_ndf);
   fChain->SetBranchAddress("CCNuPionInc_minos_trk_p", &CCNuPionInc_minos_trk_p, &b_CCNuPionInc_minos_trk_p);
   fChain->SetBranchAddress("CCNuPionInc_minos_trk_p_curvature", &CCNuPionInc_minos_trk_p_curvature, &b_CCNuPionInc_minos_trk_p_curvature);
   fChain->SetBranchAddress("CCNuPionInc_minos_trk_p_range", &CCNuPionInc_minos_trk_p_range, &b_CCNuPionInc_minos_trk_p_range);
   fChain->SetBranchAddress("CCNuPionInc_minos_trk_qp", &CCNuPionInc_minos_trk_qp, &b_CCNuPionInc_minos_trk_qp);
   fChain->SetBranchAddress("CCNuPionInc_minos_trk_vtx_x", &CCNuPionInc_minos_trk_vtx_x, &b_CCNuPionInc_minos_trk_vtx_x);
   fChain->SetBranchAddress("CCNuPionInc_minos_trk_vtx_y", &CCNuPionInc_minos_trk_vtx_y, &b_CCNuPionInc_minos_trk_vtx_y);
   fChain->SetBranchAddress("CCNuPionInc_minos_trk_vtx_z", &CCNuPionInc_minos_trk_vtx_z, &b_CCNuPionInc_minos_trk_vtx_z);
   fChain->SetBranchAddress("CCNuPionInc_muon_E", &CCNuPionInc_muon_E, &b_CCNuPionInc_muon_E);
   fChain->SetBranchAddress("CCNuPionInc_muon_E_shift", &CCNuPionInc_muon_E_shift, &b_CCNuPionInc_muon_E_shift);
   fChain->SetBranchAddress("CCNuPionInc_muon_muScore", &CCNuPionInc_muon_muScore, &b_CCNuPionInc_muon_muScore);
   fChain->SetBranchAddress("CCNuPionInc_muon_p", &CCNuPionInc_muon_p, &b_CCNuPionInc_muon_p);
   fChain->SetBranchAddress("CCNuPionInc_muon_px", &CCNuPionInc_muon_px, &b_CCNuPionInc_muon_px);
   fChain->SetBranchAddress("CCNuPionInc_muon_py", &CCNuPionInc_muon_py, &b_CCNuPionInc_muon_py);
   fChain->SetBranchAddress("CCNuPionInc_muon_pz", &CCNuPionInc_muon_pz, &b_CCNuPionInc_muon_pz);
   fChain->SetBranchAddress("CCNuPionInc_muon_qp", &CCNuPionInc_muon_qp, &b_CCNuPionInc_muon_qp);
   fChain->SetBranchAddress("CCNuPionInc_muon_qpqpe", &CCNuPionInc_muon_qpqpe, &b_CCNuPionInc_muon_qpqpe);
   fChain->SetBranchAddress("CCNuPionInc_muon_roadUpstreamEnergy", &CCNuPionInc_muon_roadUpstreamEnergy, &b_CCNuPionInc_muon_roadUpstreamEnergy);
   fChain->SetBranchAddress("CCNuPionInc_muon_theta", &CCNuPionInc_muon_theta, &b_CCNuPionInc_muon_theta);
   fChain->SetBranchAddress("CCNuPionInc_muon_theta_biasDown", &CCNuPionInc_muon_theta_biasDown, &b_CCNuPionInc_muon_theta_biasDown);
   fChain->SetBranchAddress("CCNuPionInc_muon_theta_biasUp", &CCNuPionInc_muon_theta_biasUp, &b_CCNuPionInc_muon_theta_biasUp);
   fChain->SetBranchAddress("CCNuPionInc_r_minos_trk_bdL", &CCNuPionInc_r_minos_trk_bdL, &b_CCNuPionInc_r_minos_trk_bdL);
   fChain->SetBranchAddress("CCNuPionInc_r_minos_trk_end_dcosx", &CCNuPionInc_r_minos_trk_end_dcosx, &b_CCNuPionInc_r_minos_trk_end_dcosx);
   fChain->SetBranchAddress("CCNuPionInc_r_minos_trk_end_dcosy", &CCNuPionInc_r_minos_trk_end_dcosy, &b_CCNuPionInc_r_minos_trk_end_dcosy);
   fChain->SetBranchAddress("CCNuPionInc_r_minos_trk_end_dcosz", &CCNuPionInc_r_minos_trk_end_dcosz, &b_CCNuPionInc_r_minos_trk_end_dcosz);
   fChain->SetBranchAddress("CCNuPionInc_r_minos_trk_vtx_dcosx", &CCNuPionInc_r_minos_trk_vtx_dcosx, &b_CCNuPionInc_r_minos_trk_vtx_dcosx);
   fChain->SetBranchAddress("CCNuPionInc_r_minos_trk_vtx_dcosy", &CCNuPionInc_r_minos_trk_vtx_dcosy, &b_CCNuPionInc_r_minos_trk_vtx_dcosy);
   fChain->SetBranchAddress("CCNuPionInc_r_minos_trk_vtx_dcosz", &CCNuPionInc_r_minos_trk_vtx_dcosz, &b_CCNuPionInc_r_minos_trk_vtx_dcosz);
   fChain->SetBranchAddress("CCNuPionInc_t_minos_trk_primFSLepMinosInitProjPx", &CCNuPionInc_t_minos_trk_primFSLepMinosInitProjPx, &b_CCNuPionInc_t_minos_trk_primFSLepMinosInitProjPx);
   fChain->SetBranchAddress("CCNuPionInc_t_minos_trk_primFSLepMinosInitProjPy", &CCNuPionInc_t_minos_trk_primFSLepMinosInitProjPy, &b_CCNuPionInc_t_minos_trk_primFSLepMinosInitProjPy);
   fChain->SetBranchAddress("CCNuPionInc_t_minos_trk_primFSLepMinosInitProjPz", &CCNuPionInc_t_minos_trk_primFSLepMinosInitProjPz, &b_CCNuPionInc_t_minos_trk_primFSLepMinosInitProjPz);
   fChain->SetBranchAddress("CCNuPionInc_t_minos_trk_primFSLepMinosInitProjX", &CCNuPionInc_t_minos_trk_primFSLepMinosInitProjX, &b_CCNuPionInc_t_minos_trk_primFSLepMinosInitProjX);
   fChain->SetBranchAddress("CCNuPionInc_t_minos_trk_primFSLepMinosInitProjY", &CCNuPionInc_t_minos_trk_primFSLepMinosInitProjY, &b_CCNuPionInc_t_minos_trk_primFSLepMinosInitProjY);
   fChain->SetBranchAddress("CCNuPionInc_t_minos_trk_primFSLepMinosInitProjZ", &CCNuPionInc_t_minos_trk_primFSLepMinosInitProjZ, &b_CCNuPionInc_t_minos_trk_primFSLepMinosInitProjZ);
   fChain->SetBranchAddress("CCNuPionInc_t_minos_trk_primFSLepMnvFinalPx", &CCNuPionInc_t_minos_trk_primFSLepMnvFinalPx, &b_CCNuPionInc_t_minos_trk_primFSLepMnvFinalPx);
   fChain->SetBranchAddress("CCNuPionInc_t_minos_trk_primFSLepMnvFinalPy", &CCNuPionInc_t_minos_trk_primFSLepMnvFinalPy, &b_CCNuPionInc_t_minos_trk_primFSLepMnvFinalPy);
   fChain->SetBranchAddress("CCNuPionInc_t_minos_trk_primFSLepMnvFinalPz", &CCNuPionInc_t_minos_trk_primFSLepMnvFinalPz, &b_CCNuPionInc_t_minos_trk_primFSLepMnvFinalPz);
   fChain->SetBranchAddress("CCNuPionInc_t_minos_trk_primFSLepMnvFinalX", &CCNuPionInc_t_minos_trk_primFSLepMnvFinalX, &b_CCNuPionInc_t_minos_trk_primFSLepMnvFinalX);
   fChain->SetBranchAddress("CCNuPionInc_t_minos_trk_primFSLepMnvFinalY", &CCNuPionInc_t_minos_trk_primFSLepMnvFinalY, &b_CCNuPionInc_t_minos_trk_primFSLepMnvFinalY);
   fChain->SetBranchAddress("CCNuPionInc_t_minos_trk_primFSLepMnvFinalZ", &CCNuPionInc_t_minos_trk_primFSLepMnvFinalZ, &b_CCNuPionInc_t_minos_trk_primFSLepMnvFinalZ);
   fChain->SetBranchAddress("CCNuPionInc_t_minos_trk_primFSLepMnvInitPx", &CCNuPionInc_t_minos_trk_primFSLepMnvInitPx, &b_CCNuPionInc_t_minos_trk_primFSLepMnvInitPx);
   fChain->SetBranchAddress("CCNuPionInc_t_minos_trk_primFSLepMnvInitPy", &CCNuPionInc_t_minos_trk_primFSLepMnvInitPy, &b_CCNuPionInc_t_minos_trk_primFSLepMnvInitPy);
   fChain->SetBranchAddress("CCNuPionInc_t_minos_trk_primFSLepMnvInitPz", &CCNuPionInc_t_minos_trk_primFSLepMnvInitPz, &b_CCNuPionInc_t_minos_trk_primFSLepMnvInitPz);
   fChain->SetBranchAddress("CCNuPionInc_t_minos_trk_primFSLepMnvInitX", &CCNuPionInc_t_minos_trk_primFSLepMnvInitX, &b_CCNuPionInc_t_minos_trk_primFSLepMnvInitX);
   fChain->SetBranchAddress("CCNuPionInc_t_minos_trk_primFSLepMnvInitY", &CCNuPionInc_t_minos_trk_primFSLepMnvInitY, &b_CCNuPionInc_t_minos_trk_primFSLepMnvInitY);
   fChain->SetBranchAddress("CCNuPionInc_t_minos_trk_primFSLepMnvInitZ", &CCNuPionInc_t_minos_trk_primFSLepMnvInitZ, &b_CCNuPionInc_t_minos_trk_primFSLepMnvInitZ);
   fChain->SetBranchAddress("CCNuPionInc_vtx_x", &CCNuPionInc_vtx_x, &b_CCNuPionInc_vtx_x);
   fChain->SetBranchAddress("CCNuPionInc_vtx_y", &CCNuPionInc_vtx_y, &b_CCNuPionInc_vtx_y);
   fChain->SetBranchAddress("CCNuPionInc_vtx_z", &CCNuPionInc_vtx_z, &b_CCNuPionInc_vtx_z);
   fChain->SetBranchAddress("CCNuPionInc_hadron_1stTrackPatRec", CCNuPionInc_hadron_1stTrackPatRec, &b_CCNuPionInc_hadron_1stTrackPatRec);
   fChain->SetBranchAddress("CCNuPionInc_hadron_endMichel_bgmodule", CCNuPionInc_hadron_endMichel_bgmodule, &b_CCNuPionInc_hadron_endMichel_bgmodule);
   fChain->SetBranchAddress("CCNuPionInc_hadron_endMichel_category", CCNuPionInc_hadron_endMichel_category, &b_CCNuPionInc_hadron_endMichel_category);
   fChain->SetBranchAddress("CCNuPionInc_hadron_endMichel_edmodule", CCNuPionInc_hadron_endMichel_edmodule, &b_CCNuPionInc_hadron_endMichel_edmodule);
   fChain->SetBranchAddress("CCNuPionInc_hadron_endMichel_ndigits", CCNuPionInc_hadron_endMichel_ndigits, &b_CCNuPionInc_hadron_endMichel_ndigits);
   fChain->SetBranchAddress("CCNuPionInc_hadron_endMichel_nmodules", CCNuPionInc_hadron_endMichel_nmodules, &b_CCNuPionInc_hadron_endMichel_nmodules);
   fChain->SetBranchAddress("CCNuPionInc_hadron_endMichel_nplanes", CCNuPionInc_hadron_endMichel_nplanes, &b_CCNuPionInc_hadron_endMichel_nplanes);
   fChain->SetBranchAddress("CCNuPionInc_hadron_endMichel_tm_parentpdg", CCNuPionInc_hadron_endMichel_tm_parentpdg, &b_CCNuPionInc_hadron_endMichel_tm_parentpdg);
   fChain->SetBranchAddress("CCNuPionInc_hadron_endMichel_tm_pdg", CCNuPionInc_hadron_endMichel_tm_pdg, &b_CCNuPionInc_hadron_endMichel_tm_pdg);
   fChain->SetBranchAddress("CCNuPionInc_hadron_endMichel_tm_primarypdg", CCNuPionInc_hadron_endMichel_tm_primarypdg, &b_CCNuPionInc_hadron_endMichel_tm_primarypdg);
   fChain->SetBranchAddress("CCNuPionInc_hadron_hasEndpointMichel", CCNuPionInc_hadron_hasEndpointMichel, &b_CCNuPionInc_hadron_hasEndpointMichel);
   fChain->SetBranchAddress("CCNuPionInc_hadron_hasSecondaryMichel", CCNuPionInc_hadron_hasSecondaryMichel, &b_CCNuPionInc_hadron_hasSecondaryMichel);
   fChain->SetBranchAddress("CCNuPionInc_hadron_isDsECAL", CCNuPionInc_hadron_isDsECAL, &b_CCNuPionInc_hadron_isDsECAL);
   fChain->SetBranchAddress("CCNuPionInc_hadron_isExiting", CCNuPionInc_hadron_isExiting, &b_CCNuPionInc_hadron_isExiting);
   fChain->SetBranchAddress("CCNuPionInc_hadron_isForked", CCNuPionInc_hadron_isForked, &b_CCNuPionInc_hadron_isForked);
   fChain->SetBranchAddress("CCNuPionInc_hadron_isHCAL", CCNuPionInc_hadron_isHCAL, &b_CCNuPionInc_hadron_isHCAL);
   fChain->SetBranchAddress("CCNuPionInc_hadron_isKinked", CCNuPionInc_hadron_isKinked, &b_CCNuPionInc_hadron_isKinked);
   fChain->SetBranchAddress("CCNuPionInc_hadron_isNuclTargs", CCNuPionInc_hadron_isNuclTargs, &b_CCNuPionInc_hadron_isNuclTargs);
   fChain->SetBranchAddress("CCNuPionInc_hadron_isODMatch", CCNuPionInc_hadron_isODMatch, &b_CCNuPionInc_hadron_isODMatch);
   fChain->SetBranchAddress("CCNuPionInc_hadron_isSideECAL", CCNuPionInc_hadron_isSideECAL, &b_CCNuPionInc_hadron_isSideECAL);
   fChain->SetBranchAddress("CCNuPionInc_hadron_isTracker", CCNuPionInc_hadron_isTracker, &b_CCNuPionInc_hadron_isTracker);
   fChain->SetBranchAddress("CCNuPionInc_hadron_piFit_fails", CCNuPionInc_hadron_piFit_fails, &b_CCNuPionInc_hadron_piFit_fails);
   fChain->SetBranchAddress("CCNuPionInc_hadron_piFit_isMM", CCNuPionInc_hadron_piFit_isMM, &b_CCNuPionInc_hadron_piFit_isMM);
   fChain->SetBranchAddress("CCNuPionInc_hadron_piFit_ok", CCNuPionInc_hadron_piFit_ok, &b_CCNuPionInc_hadron_piFit_ok);
   fChain->SetBranchAddress("CCNuPionInc_hadron_piFit_range", CCNuPionInc_hadron_piFit_range, &b_CCNuPionInc_hadron_piFit_range);
   fChain->SetBranchAddress("CCNuPionInc_hadron_proFit_fails", CCNuPionInc_hadron_proFit_fails, &b_CCNuPionInc_hadron_proFit_fails);
   fChain->SetBranchAddress("CCNuPionInc_hadron_proFit_isMM", CCNuPionInc_hadron_proFit_isMM, &b_CCNuPionInc_hadron_proFit_isMM);
   fChain->SetBranchAddress("CCNuPionInc_hadron_proFit_ok", CCNuPionInc_hadron_proFit_ok, &b_CCNuPionInc_hadron_proFit_ok);
   fChain->SetBranchAddress("CCNuPionInc_hadron_proFit_range", CCNuPionInc_hadron_proFit_range, &b_CCNuPionInc_hadron_proFit_range);
   fChain->SetBranchAddress("CCNuPionInc_hadron_secMichel_bgmodule", CCNuPionInc_hadron_secMichel_bgmodule, &b_CCNuPionInc_hadron_secMichel_bgmodule);
   fChain->SetBranchAddress("CCNuPionInc_hadron_secMichel_category", CCNuPionInc_hadron_secMichel_category, &b_CCNuPionInc_hadron_secMichel_category);
   fChain->SetBranchAddress("CCNuPionInc_hadron_secMichel_edmodule", CCNuPionInc_hadron_secMichel_edmodule, &b_CCNuPionInc_hadron_secMichel_edmodule);
   fChain->SetBranchAddress("CCNuPionInc_hadron_secMichel_ndigits", CCNuPionInc_hadron_secMichel_ndigits, &b_CCNuPionInc_hadron_secMichel_ndigits);
   fChain->SetBranchAddress("CCNuPionInc_hadron_secMichel_nmodules", CCNuPionInc_hadron_secMichel_nmodules, &b_CCNuPionInc_hadron_secMichel_nmodules);
   fChain->SetBranchAddress("CCNuPionInc_hadron_secMichel_nplanes", CCNuPionInc_hadron_secMichel_nplanes, &b_CCNuPionInc_hadron_secMichel_nplanes);
   fChain->SetBranchAddress("CCNuPionInc_hadron_tm_PDGCode", CCNuPionInc_hadron_tm_PDGCode, &b_CCNuPionInc_hadron_tm_PDGCode);
   fChain->SetBranchAddress("CCNuPionInc_hadron_tm_PDGCodeDaughter", CCNuPionInc_hadron_tm_PDGCodeDaughter, &b_CCNuPionInc_hadron_tm_PDGCodeDaughter);
   fChain->SetBranchAddress("CCNuPionInc_hadron_tm_destructCode", CCNuPionInc_hadron_tm_destructCode, &b_CCNuPionInc_hadron_tm_destructCode);
   fChain->SetBranchAddress("CCNuPionInc_hadron_tm_endpt_pdg", CCNuPionInc_hadron_tm_endpt_pdg, &b_CCNuPionInc_hadron_tm_endpt_pdg);
   fChain->SetBranchAddress("CCNuPionInc_hadron_tm_firstInelasticCode", CCNuPionInc_hadron_tm_firstInelasticCode, &b_CCNuPionInc_hadron_tm_firstInelasticCode);
   fChain->SetBranchAddress("CCNuPionInc_hadron_tm_isTruthMatched", CCNuPionInc_hadron_tm_isTruthMatched, &b_CCNuPionInc_hadron_tm_isTruthMatched);
   fChain->SetBranchAddress("CCNuPionInc_hadron_tm_trackID", CCNuPionInc_hadron_tm_trackID, &b_CCNuPionInc_hadron_tm_trackID);
   fChain->SetBranchAddress("CCNuPionInc_hadron_trackNodes", CCNuPionInc_hadron_trackNodes, &b_CCNuPionInc_hadron_trackNodes);
   fChain->SetBranchAddress("CCNuPionInc_hadron_calE_CCInc", CCNuPionInc_hadron_calE_CCInc, &b_CCNuPionInc_hadron_calE_CCInc);
   fChain->SetBranchAddress("CCNuPionInc_hadron_coneE", CCNuPionInc_hadron_coneE, &b_CCNuPionInc_hadron_coneE);
   fChain->SetBranchAddress("CCNuPionInc_hadron_endAvdEdX", CCNuPionInc_hadron_endAvdEdX, &b_CCNuPionInc_hadron_endAvdEdX);
   fChain->SetBranchAddress("CCNuPionInc_hadron_endMichel_bgz", CCNuPionInc_hadron_endMichel_bgz, &b_CCNuPionInc_hadron_endMichel_bgz);
   fChain->SetBranchAddress("CCNuPionInc_hadron_endMichel_bgzl", CCNuPionInc_hadron_endMichel_bgzl, &b_CCNuPionInc_hadron_endMichel_bgzl);
   fChain->SetBranchAddress("CCNuPionInc_hadron_endMichel_distance", CCNuPionInc_hadron_endMichel_distance, &b_CCNuPionInc_hadron_endMichel_distance);
   fChain->SetBranchAddress("CCNuPionInc_hadron_endMichel_edz", CCNuPionInc_hadron_endMichel_edz, &b_CCNuPionInc_hadron_endMichel_edz);
   fChain->SetBranchAddress("CCNuPionInc_hadron_endMichel_edzl", CCNuPionInc_hadron_endMichel_edzl, &b_CCNuPionInc_hadron_endMichel_edzl);
   fChain->SetBranchAddress("CCNuPionInc_hadron_endMichel_energy", CCNuPionInc_hadron_endMichel_energy, &b_CCNuPionInc_hadron_endMichel_energy);
   fChain->SetBranchAddress("CCNuPionInc_hadron_endMichel_energy_uncorrected", CCNuPionInc_hadron_endMichel_energy_uncorrected, &b_CCNuPionInc_hadron_endMichel_energy_uncorrected);
   fChain->SetBranchAddress("CCNuPionInc_hadron_endMichel_firedFraction", CCNuPionInc_hadron_endMichel_firedFraction, &b_CCNuPionInc_hadron_endMichel_firedFraction);
   fChain->SetBranchAddress("CCNuPionInc_hadron_endMichel_mcfrac", CCNuPionInc_hadron_endMichel_mcfrac, &b_CCNuPionInc_hadron_endMichel_mcfrac);
   fChain->SetBranchAddress("CCNuPionInc_hadron_endMichel_slice_energy", CCNuPionInc_hadron_endMichel_slice_energy, &b_CCNuPionInc_hadron_endMichel_slice_energy);
   fChain->SetBranchAddress("CCNuPionInc_hadron_endMichel_time_diff", CCNuPionInc_hadron_endMichel_time_diff, &b_CCNuPionInc_hadron_endMichel_time_diff);
   fChain->SetBranchAddress("CCNuPionInc_hadron_endMichel_tm_otherE", CCNuPionInc_hadron_endMichel_tm_otherE, &b_CCNuPionInc_hadron_endMichel_tm_otherE);
   fChain->SetBranchAddress("CCNuPionInc_hadron_endpointE", CCNuPionInc_hadron_endpointE, &b_CCNuPionInc_hadron_endpointE);
   fChain->SetBranchAddress("CCNuPionInc_hadron_matRange", CCNuPionInc_hadron_matRange, &b_CCNuPionInc_hadron_matRange);
   fChain->SetBranchAddress("CCNuPionInc_hadron_piFit_chi2", CCNuPionInc_hadron_piFit_chi2, &b_CCNuPionInc_hadron_piFit_chi2);
   fChain->SetBranchAddress("CCNuPionInc_hadron_piFit_lastAllScore1", CCNuPionInc_hadron_piFit_lastAllScore1, &b_CCNuPionInc_hadron_piFit_lastAllScore1);
   fChain->SetBranchAddress("CCNuPionInc_hadron_piFit_lastHalfScore1", CCNuPionInc_hadron_piFit_lastHalfScore1, &b_CCNuPionInc_hadron_piFit_lastHalfScore1);
   fChain->SetBranchAddress("CCNuPionInc_hadron_piFit_lastSixScore1", CCNuPionInc_hadron_piFit_lastSixScore1, &b_CCNuPionInc_hadron_piFit_lastSixScore1);
   fChain->SetBranchAddress("CCNuPionInc_hadron_piFit_score1", CCNuPionInc_hadron_piFit_score1, &b_CCNuPionInc_hadron_piFit_score1);
   fChain->SetBranchAddress("CCNuPionInc_hadron_piFit_score1_BetheBloch_biasDown", CCNuPionInc_hadron_piFit_score1_BetheBloch_biasDown, &b_CCNuPionInc_hadron_piFit_score1_BetheBloch_biasDown);
   fChain->SetBranchAddress("CCNuPionInc_hadron_piFit_score1_BetheBloch_biasUp", CCNuPionInc_hadron_piFit_score1_BetheBloch_biasUp, &b_CCNuPionInc_hadron_piFit_score1_BetheBloch_biasUp);
   fChain->SetBranchAddress("CCNuPionInc_hadron_piFit_score1_Birks", CCNuPionInc_hadron_piFit_score1_Birks, &b_CCNuPionInc_hadron_piFit_score1_Birks);
   fChain->SetBranchAddress("CCNuPionInc_hadron_piFit_score1_MEU_biasDown", CCNuPionInc_hadron_piFit_score1_MEU_biasDown, &b_CCNuPionInc_hadron_piFit_score1_MEU_biasDown);
   fChain->SetBranchAddress("CCNuPionInc_hadron_piFit_score1_MEU_biasUp", CCNuPionInc_hadron_piFit_score1_MEU_biasUp, &b_CCNuPionInc_hadron_piFit_score1_MEU_biasUp);
   fChain->SetBranchAddress("CCNuPionInc_hadron_piFit_score1_Mass_biasDown", CCNuPionInc_hadron_piFit_score1_Mass_biasDown, &b_CCNuPionInc_hadron_piFit_score1_Mass_biasDown);
   fChain->SetBranchAddress("CCNuPionInc_hadron_piFit_score1_Mass_biasUp", CCNuPionInc_hadron_piFit_score1_Mass_biasUp, &b_CCNuPionInc_hadron_piFit_score1_Mass_biasUp);
   fChain->SetBranchAddress("CCNuPionInc_hadron_piFit_score2", CCNuPionInc_hadron_piFit_score2, &b_CCNuPionInc_hadron_piFit_score2);
   fChain->SetBranchAddress("CCNuPionInc_hadron_piFit_zDiff", CCNuPionInc_hadron_piFit_zDiff, &b_CCNuPionInc_hadron_piFit_zDiff);
   fChain->SetBranchAddress("CCNuPionInc_hadron_pimFit_SPID", CCNuPionInc_hadron_pimFit_SPID, &b_CCNuPionInc_hadron_pimFit_SPID);
   fChain->SetBranchAddress("CCNuPionInc_hadron_pimFit_newChi2", CCNuPionInc_hadron_pimFit_newChi2, &b_CCNuPionInc_hadron_pimFit_newChi2);
   fChain->SetBranchAddress("CCNuPionInc_hadron_pion_E", CCNuPionInc_hadron_pion_E, &b_CCNuPionInc_hadron_pion_E);
   fChain->SetBranchAddress("CCNuPionInc_hadron_pion_E_BetheBloch_biasDown", CCNuPionInc_hadron_pion_E_BetheBloch_biasDown, &b_CCNuPionInc_hadron_pion_E_BetheBloch_biasDown);
   fChain->SetBranchAddress("CCNuPionInc_hadron_pion_E_BetheBloch_biasUp", CCNuPionInc_hadron_pion_E_BetheBloch_biasUp, &b_CCNuPionInc_hadron_pion_E_BetheBloch_biasUp);
   fChain->SetBranchAddress("CCNuPionInc_hadron_pion_E_Birks", CCNuPionInc_hadron_pion_E_Birks, &b_CCNuPionInc_hadron_pion_E_Birks);
   fChain->SetBranchAddress("CCNuPionInc_hadron_pion_E_MEU_biasDown", CCNuPionInc_hadron_pion_E_MEU_biasDown, &b_CCNuPionInc_hadron_pion_E_MEU_biasDown);
   fChain->SetBranchAddress("CCNuPionInc_hadron_pion_E_MEU_biasUp", CCNuPionInc_hadron_pion_E_MEU_biasUp, &b_CCNuPionInc_hadron_pion_E_MEU_biasUp);
   fChain->SetBranchAddress("CCNuPionInc_hadron_pion_E_Mass_biasDown", CCNuPionInc_hadron_pion_E_Mass_biasDown, &b_CCNuPionInc_hadron_pion_E_Mass_biasDown);
   fChain->SetBranchAddress("CCNuPionInc_hadron_pion_E_Mass_biasUp", CCNuPionInc_hadron_pion_E_Mass_biasUp, &b_CCNuPionInc_hadron_pion_E_Mass_biasUp);
   fChain->SetBranchAddress("CCNuPionInc_hadron_pion_openingAngle", CCNuPionInc_hadron_pion_openingAngle, &b_CCNuPionInc_hadron_pion_openingAngle);
   fChain->SetBranchAddress("CCNuPionInc_hadron_pion_px", CCNuPionInc_hadron_pion_px, &b_CCNuPionInc_hadron_pion_px);
   fChain->SetBranchAddress("CCNuPionInc_hadron_pion_py", CCNuPionInc_hadron_pion_py, &b_CCNuPionInc_hadron_pion_py);
   fChain->SetBranchAddress("CCNuPionInc_hadron_pion_pz", CCNuPionInc_hadron_pion_pz, &b_CCNuPionInc_hadron_pion_pz);
   fChain->SetBranchAddress("CCNuPionInc_hadron_pion_theta", CCNuPionInc_hadron_pion_theta, &b_CCNuPionInc_hadron_pion_theta);
   fChain->SetBranchAddress("CCNuPionInc_hadron_pion_theta_biasDown", CCNuPionInc_hadron_pion_theta_biasDown, &b_CCNuPionInc_hadron_pion_theta_biasDown);
   fChain->SetBranchAddress("CCNuPionInc_hadron_pion_theta_biasUp", CCNuPionInc_hadron_pion_theta_biasUp, &b_CCNuPionInc_hadron_pion_theta_biasUp);
   fChain->SetBranchAddress("CCNuPionInc_hadron_pipFit_SPID", CCNuPionInc_hadron_pipFit_SPID, &b_CCNuPionInc_hadron_pipFit_SPID);
   fChain->SetBranchAddress("CCNuPionInc_hadron_pipFit_newChi2", CCNuPionInc_hadron_pipFit_newChi2, &b_CCNuPionInc_hadron_pipFit_newChi2);
   fChain->SetBranchAddress("CCNuPionInc_hadron_proFit_SPID", CCNuPionInc_hadron_proFit_SPID, &b_CCNuPionInc_hadron_proFit_SPID);
   fChain->SetBranchAddress("CCNuPionInc_hadron_proFit_chi2", CCNuPionInc_hadron_proFit_chi2, &b_CCNuPionInc_hadron_proFit_chi2);
   fChain->SetBranchAddress("CCNuPionInc_hadron_proFit_newChi2", CCNuPionInc_hadron_proFit_newChi2, &b_CCNuPionInc_hadron_proFit_newChi2);
   fChain->SetBranchAddress("CCNuPionInc_hadron_proFit_score1", CCNuPionInc_hadron_proFit_score1, &b_CCNuPionInc_hadron_proFit_score1);
   fChain->SetBranchAddress("CCNuPionInc_hadron_proFit_score2", CCNuPionInc_hadron_proFit_score2, &b_CCNuPionInc_hadron_proFit_score2);
   fChain->SetBranchAddress("CCNuPionInc_hadron_proFit_zDiff", CCNuPionInc_hadron_proFit_zDiff, &b_CCNuPionInc_hadron_proFit_zDiff);
   fChain->SetBranchAddress("CCNuPionInc_hadron_rawRange", CCNuPionInc_hadron_rawRange, &b_CCNuPionInc_hadron_rawRange);
   fChain->SetBranchAddress("CCNuPionInc_hadron_secMichel_bgz", CCNuPionInc_hadron_secMichel_bgz, &b_CCNuPionInc_hadron_secMichel_bgz);
   fChain->SetBranchAddress("CCNuPionInc_hadron_secMichel_bgzl", CCNuPionInc_hadron_secMichel_bgzl, &b_CCNuPionInc_hadron_secMichel_bgzl);
   fChain->SetBranchAddress("CCNuPionInc_hadron_secMichel_distance", CCNuPionInc_hadron_secMichel_distance, &b_CCNuPionInc_hadron_secMichel_distance);
   fChain->SetBranchAddress("CCNuPionInc_hadron_secMichel_edz", CCNuPionInc_hadron_secMichel_edz, &b_CCNuPionInc_hadron_secMichel_edz);
   fChain->SetBranchAddress("CCNuPionInc_hadron_secMichel_edzl", CCNuPionInc_hadron_secMichel_edzl, &b_CCNuPionInc_hadron_secMichel_edzl);
   fChain->SetBranchAddress("CCNuPionInc_hadron_secMichel_energy", CCNuPionInc_hadron_secMichel_energy, &b_CCNuPionInc_hadron_secMichel_energy);
   fChain->SetBranchAddress("CCNuPionInc_hadron_secMichel_energy_uncorrected", CCNuPionInc_hadron_secMichel_energy_uncorrected, &b_CCNuPionInc_hadron_secMichel_energy_uncorrected);
   fChain->SetBranchAddress("CCNuPionInc_hadron_secMichel_firedFraction", CCNuPionInc_hadron_secMichel_firedFraction, &b_CCNuPionInc_hadron_secMichel_firedFraction);
   fChain->SetBranchAddress("CCNuPionInc_hadron_secMichel_slice_energy", CCNuPionInc_hadron_secMichel_slice_energy, &b_CCNuPionInc_hadron_secMichel_slice_energy);
   fChain->SetBranchAddress("CCNuPionInc_hadron_secMichel_time_diff", CCNuPionInc_hadron_secMichel_time_diff, &b_CCNuPionInc_hadron_secMichel_time_diff);
   fChain->SetBranchAddress("CCNuPionInc_hadron_tm_beginMomentum", CCNuPionInc_hadron_tm_beginMomentum, &b_CCNuPionInc_hadron_tm_beginMomentum);
   fChain->SetBranchAddress("CCNuPionInc_hadron_tm_endMomentum", CCNuPionInc_hadron_tm_endMomentum, &b_CCNuPionInc_hadron_tm_endMomentum);
   fChain->SetBranchAddress("CCNuPionInc_hadron_tm_endpt_fraction", CCNuPionInc_hadron_tm_endpt_fraction, &b_CCNuPionInc_hadron_tm_endpt_fraction);
   fChain->SetBranchAddress("CCNuPionInc_hadron_tm_endpt_otherE", CCNuPionInc_hadron_tm_endpt_otherE, &b_CCNuPionInc_hadron_tm_endpt_otherE);
   fChain->SetBranchAddress("CCNuPionInc_hadron_tm_firstInelasticDistance", CCNuPionInc_hadron_tm_firstInelasticDistance, &b_CCNuPionInc_hadron_tm_firstInelasticDistance);
   fChain->SetBranchAddress("CCNuPionInc_hadron_tm_firstInelasticEndP", CCNuPionInc_hadron_tm_firstInelasticEndP, &b_CCNuPionInc_hadron_tm_firstInelasticEndP);
   fChain->SetBranchAddress("CCNuPionInc_hadron_tm_fraction", CCNuPionInc_hadron_tm_fraction, &b_CCNuPionInc_hadron_tm_fraction);
   fChain->SetBranchAddress("CCNuPionInc_hadron_tm_otherE", CCNuPionInc_hadron_tm_otherE, &b_CCNuPionInc_hadron_tm_otherE);
   fChain->SetBranchAddress("CCNuPionInc_hadron_tm_pathlength", CCNuPionInc_hadron_tm_pathlength, &b_CCNuPionInc_hadron_tm_pathlength);
   fChain->SetBranchAddress("CCNuPionInc_hadron_visibleE", CCNuPionInc_hadron_visibleE, &b_CCNuPionInc_hadron_visibleE);
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

Bool_t ANA_CCNuPionInc::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

    cout<<"Initialization is OK!"<<endl;
   return kTRUE;
}

void ANA_CCNuPionInc::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

Int_t ANA_CCNuPionInc::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ANA_CCNuPionInc_cxx

