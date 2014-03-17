/*
    See CCInclusive.h header for Class Information
*/

#define CCInclusive_cxx

#include "CCInclusive.h"
#include "Cuts.cpp"


CCInclusive::CCInclusive()
{



}


void CCInclusive::run(string playlist, string rootFileName)
{

   
 
    //------------------------------------------------------------------------
    // Open files for writing
    //------------------------------------------------------------------------

    openFiles();

    TFile* f = new TFile(rootFileName.c_str(),"RECREATE");

    //------------------------------------------------------------------------
    // Create chain
    //------------------------------------------------------------------------

    TChain* fChain = new TChain("CCInclusiveReco");
    Init(playlist, fChain);

    if (!fChain) return;
    if (fChain == 0) return;
    

    //------------------------------------------------------------------------
    // Initialize the Analysis Variables and Histograms
    //------------------------------------------------------------------------

    initVariables();
    initHistograms();
    
    //------------------------------------------------------------------------
    // Declare Other Variables
    //------------------------------------------------------------------------

    //------------------------------------------------------------------------
    // Branch Selection for Performance
    //------------------------------------------------------------------------
    fChain->SetBranchStatus("*",0);  // disable all branches
    
    // Analysis Variables
    fChain->SetBranchStatus("ev_run",1);  // activate
    fChain->SetBranchStatus("ev_subrun",1);  // activate
    fChain->SetBranchStatus("ev_gate",1);  // activate
    fChain->SetBranchStatus("minos_track_match",1);  // activate
    fChain->SetBranchStatus("mc_Q2",1);  // activate
    
    
    // Cut Variables
    fChain->SetBranchStatus("mc_vtx",1);
    
    // Incoming Particle
    fChain->SetBranchStatus("mc_incomingPartVec",1);  // activate
    fChain->SetBranchStatus("mc_incomingE",1);  // activate
   
   
    // Final State Particles
    fChain->SetBranchStatus("mc_FSPartPx",1);  // activate
    fChain->SetBranchStatus("mc_FSPartPy",1);  // activate
    fChain->SetBranchStatus("mc_FSPartPz",1);  // activate
    fChain->SetBranchStatus("mc_FSPartE",1);  // activate
    fChain->SetBranchStatus("mc_nFSPart",1);  // activate
    fChain->SetBranchStatus("mc_FSPartPDG",1);


    double nMinosMatch = 0;
    double nMinosMatch_LongProton = 0;
    double nMinosMatch_ShortProton = 0;
    
    double nMinosMiss = 0;
    double nMinosMiss_LongProton = 0;
    double nMinosMiss_ShortProton = 0;
    
    double nAll = 0;
    double nVolume = 0;
    double nBeamEnergy = 0;
    double nMuon = 0;
    double nProton = 0;
    double nPion = 0;
    double nMinos = 0;
    
    int indMuon;
    int indProton;
    int indPion;

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
        if (jentry%5000 == 0){
            cout<<"\tEntry "<<jentry<<endl;
        }
    
    
        //----------------------------------------------------------------------
        // Apply Cuts
        //----------------------------------------------------------------------
        
        // Count All Events before Cuts
        nAll++;
        
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
        indMuon = findParticle(PDG_List::mu_minus);
        if(indMuon == -1){
            continue;
        }
        nMuon++;
        
        // Proton Cut
        indProton = findProton();
        if(indProton == -1){
            continue;
        }
        nProton++;
        
        
//         // Pion Cut
//         indPion = findPion();
//         if(indPion == -1){
//             continue;
//         }
//         nPion++;
        
        if(!isNoMeson()){
            continue;
        }
        nPion++;

        
        // Minos Match Cut
        if( isMinosMatched() ){
            continue;
        }
        nMinos++;


        if ( isDataAnalysis){
        //----------------------------------------------------------------------
        // Fill Particles
        //----------------------------------------------------------------------
        
            if( isMC ){
                fillParticleTrue(muon);
                fillParticleTrue(proton);
                fillParticleTrue(pion);
            }
        
            // Fill Reconstructed Information
            fillMuon();
            fillProton();
            fillPion();
            
            
        //----------------------------------------------------------------------
        // Fill Histograms
        //----------------------------------------------------------------------
       
            fillHistograms();            
        }
        
        
        
        cout<<mc_Q2<<endl;
        
        

        
    }    
    cout<<"Done!"<<endl;
    
    nCutList->nAll->setValue(nAll);
    nCutList->nVolume->setValue(nVolume);
    nCutList->nBeamEnergy->setValue(nBeamEnergy);
    nCutList->nMuon->setValue(nMuon);
    nCutList->nProton->setValue(nProton);
    nCutList->nPion->setValue(nPion);
    nCutList->nMinos->setValue(nMinos);
    
    
    // Write the Root File
    cout<<">> Writing "<<rootFileName<<endl;
//     f->Write();
    
    nCutList->writeCutTable();
    
    
    closeFiles();
   
    
    
    
    
    // Delete Dynamic Variables - I will modify the destructor for
    // better Memory Management
    delete muon;
    delete proton;
    delete pion;
    delete binList;
    delete nCutList;

}

// -------------------------------------------------------------------------
//     Specific Functions
//--------------------------------------------------------------------------


void CCInclusive::fillMuon()
{
    

}

void CCInclusive::fillProton()
{

}

void CCInclusive::fillPion()
{

}

void CCInclusive::fillHistograms()
{
    beamEnergy->Fill(mc_incomingE);
    int_channel->Fill(mc_intType);
    
    P_muon->Fill(muon->p4[1].P());
    P_proton->Fill(proton->p4[1].P());
    P_pion->Fill(pion->p4[1].P());
    
    KE_muon->Fill(muon->get_KineticEnergy(isMC));
    KE_proton->Fill(proton->get_KineticEnergy(isMC));
    KE_pion->Fill(pion->get_KineticEnergy(isMC));
    
    Angle_muon->Fill(muon->angleBeam[1] * TMath::RadToDeg());
    Angle_proton->Fill(proton->angleBeam[1] * TMath::RadToDeg());
    Angle_pion->Fill(pion->angleBeam[1] * TMath::RadToDeg());
    
    AngleMuon_muon->Fill(muon->angleMuon[1] * TMath::RadToDeg());
    AngleMuon_proton->Fill(proton->angleMuon[1] * TMath::RadToDeg());
    AngleMuon_pion->Fill(pion->angleMuon[1] * TMath::RadToDeg());

}

void CCInclusive::fillParticleTrue(Particle* part)
{
    int ind = part->ind;
    double angleBeam;
    double angleMuon;
    
    // Fill 4-Momentum
    part->p4[1].SetPxPyPzE(mc_FSPartPx[ind],mc_FSPartPy[ind],mc_FSPartPz[ind],mc_FSPartE[ind]);
       
    // set Angle wrt Beam
    part->set_angleBeam(*beam_p3);
    
    // set Angle wrt Muon
    part->set_angleMuon(muon->p4[1].Vect());
    
    
}


void CCInclusive::initHistograms()
{
    cout<<"Initializing Histograms"<<endl;
    // -------------------------------------------------------------------------
    //     True Analysis
    //--------------------------------------------------------------------------
    beamEnergy = new TH1F( "beamEnergy","True Neutrino Energy",binList->neutrinoP->get_nBins(), binList->neutrinoP->get_min(), binList->neutrinoP->get_max() );
    beamEnergy->GetXaxis()->SetTitle("True Neutrino Energy MeV");
    beamEnergy->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList->neutrinoP->get_width()));
    
    int_channel = new TH1F( "int_channel","Interaction Channel",binList->int_channel->get_nBins(), binList->int_channel->get_min(), binList->int_channel->get_max() );
    int_channel->GetXaxis()->SetTitle("1 = QE, 2 = Resonant, 3 = DIS, 4 = Coh pi");
    int_channel->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList->int_channel->get_width()));

    // -------------------------------------------------------------------------
    //     Momentum
    //--------------------------------------------------------------------------
    P_muon = new TH1F( "P_muon","True Muon Momentum",binList->muonP->get_nBins(), binList->muonP->get_min(), binList->muonP->get_max() );
    P_muon->GetXaxis()->SetTitle("True Muon Momentum MeV");
    P_muon->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList->muonP->get_width()));
    
    P_proton = new TH1F( "P_proton","True Proton Momentum",binList->protonP->get_nBins(), binList->protonP->get_min(), binList->protonP->get_max() );
    P_proton->GetXaxis()->SetTitle("True Proton Momentum MeV");
    P_proton->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList->protonP->get_width()));
    
    P_pion = new TH1F( "P_pion","True Pion Momentum",binList->pionP->get_nBins(), binList->pionP->get_min(), binList->pionP->get_max() );
    P_pion->GetXaxis()->SetTitle("True Pion Momentum MeV");
    P_pion->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList->pionP->get_width()));
    
    
    // -------------------------------------------------------------------------
    //     Kinetic Energy
    //--------------------------------------------------------------------------
    KE_muon = new TH1F( "KE_muon","True Muon Kinetic Energy",binList->muonP->get_nBins(), binList->muonP->get_min(), binList->muonP->get_max() );
    KE_muon->GetXaxis()->SetTitle("True Muon Kinetic Energy MeV");
    KE_muon->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList->muonP->get_width()));
    
    KE_proton = new TH1F( "KE_proton","True Proton Kinetic Energy",binList->protonP->get_nBins(), binList->protonP->get_min(), binList->protonP->get_max() );
    KE_proton->GetXaxis()->SetTitle("True Proton Kinetic Energy MeV");
    KE_proton->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList->protonP->get_width()));
    
    KE_pion = new TH1F( "KE_pion","True Pion Kinetic Energy",binList->pionP->get_nBins(), binList->pionP->get_min(), binList->pionP->get_max() );
    KE_pion->GetXaxis()->SetTitle("True Pion Kinetic Energy MeV");
    KE_pion->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList->pionP->get_width()));
    
    
    // -------------------------------------------------------------------------
    //     Angles
    //--------------------------------------------------------------------------
    Angle_muon = new TH1F( "Angle_muon","Angle: Beam vs Muon",binList->angle->get_nBins(), binList->angle->get_min(), binList->angle->get_max() );
    Angle_muon->GetXaxis()->SetTitle("Angle");
    Angle_muon->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList->angle->get_width()));
    
    Angle_proton = new TH1F( "Angle_proton","Angle: Beam vs Proton",binList->angle->get_nBins(), binList->angle->get_min(), binList->angle->get_max() );
    Angle_proton->GetXaxis()->SetTitle("Angle");
    Angle_proton->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList->angle->get_width()));
    
    Angle_pion = new TH1F( "Angle_pion","Angle: Beam vs Pion",binList->angle->get_nBins(), binList->angle->get_min(), binList->angle->get_max() );
    Angle_pion->GetXaxis()->SetTitle("Angle");
    Angle_pion->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList->angle->get_width()));
    
    AngleMuon_muon = new TH1F( "AngleMuon_muon","Angle: Muon vs Muon",binList->angle->get_nBins(), binList->angle->get_min(), binList->angle->get_max() );
    AngleMuon_muon->GetXaxis()->SetTitle("Angle");
    AngleMuon_muon->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList->angle->get_width()));
    
    AngleMuon_proton = new TH1F( "AngleMuon_proton","Angle: Muon vs Proton",binList->angle->get_nBins(), binList->angle->get_min(), binList->angle->get_max() );
    AngleMuon_proton->GetXaxis()->SetTitle("Angle");
    AngleMuon_proton->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList->angle->get_width()));
    
    AngleMuon_pion = new TH1F( "AngleMuon_pion","Angle: Muon vs Pion",binList->angle->get_nBins(), binList->angle->get_min(), binList->angle->get_max() );
    AngleMuon_pion->GetXaxis()->SetTitle("Angle");
    AngleMuon_pion->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",binList->angle->get_width()));
    
    
    cout<<"Done!"<<endl;
}

void CCInclusive::initSingleHistogram(TH1F* hist, string histName, string title, string xLabel, string yLabel, SingleBin* bin)
{
    // Reserved for a future version - Still needs testing
    
    cout<<bin->get_nBins()<<endl;
    hist = new TH1F(histName.c_str(),title.c_str(),bin->get_nBins(), bin->get_min(), bin->get_max());
    hist->GetXaxis()->SetTitle(xLabel.c_str());
    hist->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin->get_width()));
}




void CCInclusive::initVariables()
{
    cout<<"Initializing Variables"<<endl;
    
    isDataAnalysis = true;
    isMC = true;

    // -------------------------------------------------------------------------
    //     Memory Allocation
    //--------------------------------------------------------------------------
    // Allocate Memory
    beam_p3 = new TVector3;
    muon = new Muon;
    proton = new Proton;
    pion = new Pion;
    
    nCutList = new CutNumberList;
    binList = new BinList;
    
    nCutList->printList();
    
    // -------------------------------------------------------------------------
    //     Initialization
    //--------------------------------------------------------------------------
    // Default Beam Configuration
    beam_p3->SetXYZ(1.0,1.0,1.0);
    beam_p3->SetPhi(-1.554);
    beam_p3->SetTheta(0.059);
    
    max_nFSPart = 15;
    maxBeamEnergy = 20000; //MeV
    

    cout<<"Done!"<<endl;

}

void CCInclusive::closeFiles()
{
    readme.close();
}

void CCInclusive::openFiles()
{
    readmeFile = getFileLocation(Folder_List::OUT, Folder_List::TEXT, Folder_List::F_TEXT_README);

    readme.open( readmeFile.c_str() );
    if( !readme.is_open() ){
        cerr<<"Cannot open output text file: "<<readmeFile<<endl;
        exit(1);
    }

    writeReadme();
}

void CCInclusive::writeReadme()
{

readme<<"Test"<<endl;

}

// -------------------------------------------------------------------------
//     Default Functions
//--------------------------------------------------------------------------

#ifdef CCInclusive_cxx
CCInclusive::CCInclusive()
{
    
}

CCInclusive::~CCInclusive()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t CCInclusive::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t CCInclusive::LoadTree(Long64_t entry)
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

void CCInclusive::Init(string playlist, TChain* fChain)
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
   fChain->SetBranchAddress("broken_track_most_us_plane", &broken_track_most_us_plane, &b_broken_track_most_us_plane);
   fChain->SetBranchAddress("n_associated_recoil_subprongs", &n_associated_recoil_subprongs, &b_n_associated_recoil_subprongs);
   fChain->SetBranchAddress("n_muon_blobs", &n_muon_blobs, &b_n_muon_blobs);
   fChain->SetBranchAddress("n_recoil_iso2Dblobs", &n_recoil_iso2Dblobs, &b_n_recoil_iso2Dblobs);
   fChain->SetBranchAddress("n_recoil_iso3Dblobs", &n_recoil_iso3Dblobs, &b_n_recoil_iso3Dblobs);
   fChain->SetBranchAddress("n_recoil_subprongs", &n_recoil_subprongs, &b_n_recoil_subprongs);
   fChain->SetBranchAddress("n_recoil_vtxblobs", &n_recoil_vtxblobs, &b_n_recoil_vtxblobs);
   fChain->SetBranchAddress("n_tracks", &n_tracks, &b_n_tracks);
   fChain->SetBranchAddress("n_tracks_non_prim", &n_tracks_non_prim, &b_n_tracks_non_prim);
   fChain->SetBranchAddress("n_tracks_prim", &n_tracks_prim, &b_n_tracks_prim);
   fChain->SetBranchAddress("n_tracks_prim_forked", &n_tracks_prim_forked, &b_n_tracks_prim_forked);
   fChain->SetBranchAddress("n_tracks_prim_kinked", &n_tracks_prim_kinked, &b_n_tracks_prim_kinked);
   fChain->SetBranchAddress("n_vertices_startpoint", &n_vertices_startpoint, &b_n_vertices_startpoint);
   fChain->SetBranchAddress("pass_canonical_cut", &pass_canonical_cut, &b_pass_canonical_cut);
   fChain->SetBranchAddress("pass_minosmatch_cut", &pass_minosmatch_cut, &b_pass_minosmatch_cut);
   fChain->SetBranchAddress("pass_nearccinclusive_cut", &pass_nearccinclusive_cut, &b_pass_nearccinclusive_cut);
   fChain->SetBranchAddress("phys_energy_in_road_downstream_nplanes", &phys_energy_in_road_downstream_nplanes, &b_phys_energy_in_road_downstream_nplanes);
   fChain->SetBranchAddress("phys_energy_in_road_upstream_nplanes", &phys_energy_in_road_upstream_nplanes, &b_phys_energy_in_road_upstream_nplanes);
   fChain->SetBranchAddress("phys_n_dead_discr_pair", &phys_n_dead_discr_pair, &b_phys_n_dead_discr_pair);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_in_prim_track_region", &phys_n_dead_discr_pair_in_prim_track_region, &b_phys_n_dead_discr_pair_in_prim_track_region);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_two_mod_downstream_prim_track", &phys_n_dead_discr_pair_two_mod_downstream_prim_track, &b_phys_n_dead_discr_pair_two_mod_downstream_prim_track);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_two_mod_upstream_prim_vtx", &phys_n_dead_discr_pair_two_mod_upstream_prim_vtx, &b_phys_n_dead_discr_pair_two_mod_upstream_prim_vtx);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_upstream_prim_track_proj", &phys_n_dead_discr_pair_upstream_prim_track_proj, &b_phys_n_dead_discr_pair_upstream_prim_track_proj);
   fChain->SetBranchAddress("phys_prim_vertex_planeid", &phys_prim_vertex_planeid, &b_phys_prim_vertex_planeid);
   fChain->SetBranchAddress("phys_vertex_is_fiducial", &phys_vertex_is_fiducial, &b_phys_vertex_is_fiducial);
   fChain->SetBranchAddress("energy_from_mc", &energy_from_mc, &b_energy_from_mc);
   fChain->SetBranchAddress("energy_from_mc_fraction", &energy_from_mc_fraction, &b_energy_from_mc_fraction);
   fChain->SetBranchAddress("energy_from_mc_fraction_of_highest", &energy_from_mc_fraction_of_highest, &b_energy_from_mc_fraction_of_highest);
   fChain->SetBranchAddress("muon_blob_energy", &muon_blob_energy, &b_muon_blob_energy);
   fChain->SetBranchAddress("muon_phi", &muon_phi, &b_muon_phi);
   fChain->SetBranchAddress("muon_theta", &muon_theta, &b_muon_theta);
   fChain->SetBranchAddress("muon_thetaX", &muon_thetaX, &b_muon_thetaX);
   fChain->SetBranchAddress("muon_thetaY", &muon_thetaY, &b_muon_thetaY);
   fChain->SetBranchAddress("nu_topological_energy_recoil", &nu_topological_energy_recoil, &b_nu_topological_energy_recoil);
   fChain->SetBranchAddress("numi_horn_curr", &numi_horn_curr, &b_numi_horn_curr);
   fChain->SetBranchAddress("numi_pot", &numi_pot, &b_numi_pot);
   fChain->SetBranchAddress("numi_x", &numi_x, &b_numi_x);
   fChain->SetBranchAddress("numi_x_width", &numi_x_width, &b_numi_x_width);
   fChain->SetBranchAddress("numi_y", &numi_y, &b_numi_y);
   fChain->SetBranchAddress("numi_y_width", &numi_y_width, &b_numi_y_width);
   fChain->SetBranchAddress("phys_energy_dispersed", &phys_energy_dispersed, &b_phys_energy_dispersed);
   fChain->SetBranchAddress("phys_energy_in_road_downstream", &phys_energy_in_road_downstream, &b_phys_energy_in_road_downstream);
   fChain->SetBranchAddress("phys_energy_in_road_upstream", &phys_energy_in_road_upstream, &b_phys_energy_in_road_upstream);
   fChain->SetBranchAddress("phys_energy_unattached", &phys_energy_unattached, &b_phys_energy_unattached);
   fChain->SetBranchAddress("prim_vtx_smallest_opening_angle", &prim_vtx_smallest_opening_angle, &b_prim_vtx_smallest_opening_angle);
   fChain->SetBranchAddress("primary_track_minerva_energy", &primary_track_minerva_energy, &b_primary_track_minerva_energy);
   fChain->SetBranchAddress("primary_track_minerva_phi", &primary_track_minerva_phi, &b_primary_track_minerva_phi);
   fChain->SetBranchAddress("primary_track_minerva_theta", &primary_track_minerva_theta, &b_primary_track_minerva_theta);
   fChain->SetBranchAddress("vtxprong_energy_cal", &vtxprong_energy_cal, &b_vtxprong_energy_cal);
   fChain->SetBranchAddress("vtxprong_energy_visible", &vtxprong_energy_visible, &b_vtxprong_energy_visible);
   fChain->SetBranchAddress("vtxprong_energy_visible_isoblobs", &vtxprong_energy_visible_isoblobs, &b_vtxprong_energy_visible_isoblobs);
   fChain->SetBranchAddress("vtxprong_energy_visible_vtxblobs", &vtxprong_energy_visible_vtxblobs, &b_vtxprong_energy_visible_vtxblobs);
   fChain->SetBranchAddress("plane_id_sz", &plane_id_sz, &b_plane_id_sz);
   fChain->SetBranchAddress("plane_id", plane_id, &b_plane_id);
   fChain->SetBranchAddress("n_vtxprong_isoblobs", &n_vtxprong_isoblobs, &b_n_vtxprong_isoblobs);
   fChain->SetBranchAddress("vtxprong_isoblob_nclusters", vtxprong_isoblob_nclusters, &b_vtxprong_isoblob_nclusters);
   fChain->SetBranchAddress("n_vtxprong_vtxblobs", &n_vtxprong_vtxblobs, &b_n_vtxprong_vtxblobs);
   fChain->SetBranchAddress("vtxprong_vtxblob_nclusters", vtxprong_vtxblob_nclusters, &b_vtxprong_vtxblob_nclusters);
   fChain->SetBranchAddress("clusterU_Angle_sz", &clusterU_Angle_sz, &b_clusterU_Angle_sz);
   fChain->SetBranchAddress("clusterU_Angle", clusterU_Angle, &b_clusterU_Angle);
   fChain->SetBranchAddress("clusterU_Radius_sz", &clusterU_Radius_sz, &b_clusterU_Radius_sz);
   fChain->SetBranchAddress("clusterU_Radius", clusterU_Radius, &b_clusterU_Radius);
   fChain->SetBranchAddress("clusterU_timeDiff_sz", &clusterU_timeDiff_sz, &b_clusterU_timeDiff_sz);
   fChain->SetBranchAddress("clusterU_timeDiff", clusterU_timeDiff, &b_clusterU_timeDiff);
   fChain->SetBranchAddress("clusterU_viewDist_sz", &clusterU_viewDist_sz, &b_clusterU_viewDist_sz);
   fChain->SetBranchAddress("clusterU_viewDist", clusterU_viewDist, &b_clusterU_viewDist);
   fChain->SetBranchAddress("clusterU_visE_binned_sz", &clusterU_visE_binned_sz, &b_clusterU_visE_binned_sz);
   fChain->SetBranchAddress("clusterU_visE_binned", clusterU_visE_binned, &b_clusterU_visE_binned);
   fChain->SetBranchAddress("clusterU_visEnergy_sz", &clusterU_visEnergy_sz, &b_clusterU_visEnergy_sz);
   fChain->SetBranchAddress("clusterU_visEnergy", clusterU_visEnergy, &b_clusterU_visEnergy);
   fChain->SetBranchAddress("clusterU_zDist_sz", &clusterU_zDist_sz, &b_clusterU_zDist_sz);
   fChain->SetBranchAddress("clusterU_zDist", clusterU_zDist, &b_clusterU_zDist);
   fChain->SetBranchAddress("clusterV_Angle_sz", &clusterV_Angle_sz, &b_clusterV_Angle_sz);
   fChain->SetBranchAddress("clusterV_Angle", clusterV_Angle, &b_clusterV_Angle);
   fChain->SetBranchAddress("clusterV_Radius_sz", &clusterV_Radius_sz, &b_clusterV_Radius_sz);
   fChain->SetBranchAddress("clusterV_Radius", clusterV_Radius, &b_clusterV_Radius);
   fChain->SetBranchAddress("clusterV_timeDiff_sz", &clusterV_timeDiff_sz, &b_clusterV_timeDiff_sz);
   fChain->SetBranchAddress("clusterV_timeDiff", clusterV_timeDiff, &b_clusterV_timeDiff);
   fChain->SetBranchAddress("clusterV_viewDist_sz", &clusterV_viewDist_sz, &b_clusterV_viewDist_sz);
   fChain->SetBranchAddress("clusterV_viewDist", clusterV_viewDist, &b_clusterV_viewDist);
   fChain->SetBranchAddress("clusterV_visE_binned_sz", &clusterV_visE_binned_sz, &b_clusterV_visE_binned_sz);
   fChain->SetBranchAddress("clusterV_visE_binned", clusterV_visE_binned, &b_clusterV_visE_binned);
   fChain->SetBranchAddress("clusterV_visEnergy_sz", &clusterV_visEnergy_sz, &b_clusterV_visEnergy_sz);
   fChain->SetBranchAddress("clusterV_visEnergy", clusterV_visEnergy, &b_clusterV_visEnergy);
   fChain->SetBranchAddress("clusterV_zDist_sz", &clusterV_zDist_sz, &b_clusterV_zDist_sz);
   fChain->SetBranchAddress("clusterV_zDist", clusterV_zDist, &b_clusterV_zDist);
   fChain->SetBranchAddress("clusterX_Angle_sz", &clusterX_Angle_sz, &b_clusterX_Angle_sz);
   fChain->SetBranchAddress("clusterX_Angle", clusterX_Angle, &b_clusterX_Angle);
   fChain->SetBranchAddress("clusterX_Radius_sz", &clusterX_Radius_sz, &b_clusterX_Radius_sz);
   fChain->SetBranchAddress("clusterX_Radius", clusterX_Radius, &b_clusterX_Radius);
   fChain->SetBranchAddress("clusterX_timeDiff_sz", &clusterX_timeDiff_sz, &b_clusterX_timeDiff_sz);
   fChain->SetBranchAddress("clusterX_timeDiff", clusterX_timeDiff, &b_clusterX_timeDiff);
   fChain->SetBranchAddress("clusterX_viewDist_sz", &clusterX_viewDist_sz, &b_clusterX_viewDist_sz);
   fChain->SetBranchAddress("clusterX_viewDist", clusterX_viewDist, &b_clusterX_viewDist);
   fChain->SetBranchAddress("clusterX_visE_binned_sz", &clusterX_visE_binned_sz, &b_clusterX_visE_binned_sz);
   fChain->SetBranchAddress("clusterX_visE_binned", clusterX_visE_binned, &b_clusterX_visE_binned);
   fChain->SetBranchAddress("clusterX_visEnergy_sz", &clusterX_visEnergy_sz, &b_clusterX_visEnergy_sz);
   fChain->SetBranchAddress("clusterX_visEnergy", clusterX_visEnergy, &b_clusterX_visEnergy);
   fChain->SetBranchAddress("clusterX_zDist_sz", &clusterX_zDist_sz, &b_clusterX_zDist_sz);
   fChain->SetBranchAddress("clusterX_zDist", clusterX_zDist, &b_clusterX_zDist);
   fChain->SetBranchAddress("plane_visible_energy_sz", &plane_visible_energy_sz, &b_plane_visible_energy_sz);
   fChain->SetBranchAddress("plane_visible_energy", plane_visible_energy, &b_plane_visible_energy);
   fChain->SetBranchAddress("primary_track_minerva_end_position", primary_track_minerva_end_position, &b_primary_track_minerva_end_position);
   fChain->SetBranchAddress("primary_track_minerva_start_position", primary_track_minerva_start_position, &b_primary_track_minerva_start_position);
   fChain->SetBranchAddress("vtxprong_isoblob_visenergy", vtxprong_isoblob_visenergy, &b_vtxprong_isoblob_visenergy);
   fChain->SetBranchAddress("vtxprong_vtxblob_visenergy", vtxprong_vtxblob_visenergy, &b_vtxprong_vtxblob_visenergy);
   fChain->SetBranchAddress("truth_has_physics_event", &truth_has_physics_event, &b_truth_has_physics_event);
   fChain->SetBranchAddress("truth_pass_CCInclusiveReco", &truth_pass_CCInclusiveReco, &b_truth_pass_CCInclusiveReco);
   fChain->SetBranchAddress("truth_pass_plausible", &truth_pass_plausible, &b_truth_pass_plausible);
   fChain->SetBranchAddress("truth_pass_fiducial", &truth_pass_fiducial, &b_truth_pass_fiducial);
   fChain->SetBranchAddress("truth_pass_analyzable", &truth_pass_analyzable, &b_truth_pass_analyzable);
   fChain->SetBranchAddress("truth_pass_fiducial_apothem", &truth_pass_fiducial_apothem, &b_truth_pass_fiducial_apothem);
   fChain->SetBranchAddress("truth_pass_analyzable_apothem", &truth_pass_analyzable_apothem, &b_truth_pass_analyzable_apothem);
   fChain->SetBranchAddress("truth_pass_ClassifyVertex", &truth_pass_ClassifyVertex, &b_truth_pass_ClassifyVertex);
   fChain->SetBranchAddress("truth_pass_ClassifyCC", &truth_pass_ClassifyCC, &b_truth_pass_ClassifyCC);
   fChain->SetBranchAddress("truth_pass_MuonEReconstructed", &truth_pass_MuonEReconstructed, &b_truth_pass_MuonEReconstructed);
   fChain->SetBranchAddress("truth_pass_MuMinusSignSelection", &truth_pass_MuMinusSignSelection, &b_truth_pass_MuMinusSignSelection);
   fChain->SetBranchAddress("truth_pass_MuPlusSignSelection", &truth_pass_MuPlusSignSelection, &b_truth_pass_MuPlusSignSelection);
   fChain->SetBranchAddress("truth_pass_ClassifyTracks", &truth_pass_ClassifyTracks, &b_truth_pass_ClassifyTracks);
   fChain->SetBranchAddress("truth_pass_ClassifyRecoil", &truth_pass_ClassifyRecoil, &b_truth_pass_ClassifyRecoil);
   fChain->SetBranchAddress("truth_pass_Canonical", &truth_pass_Canonical, &b_truth_pass_Canonical);
   fChain->SetBranchAddress("truth_reco_vertex_fiducial", &truth_reco_vertex_fiducial, &b_truth_reco_vertex_fiducial);
   fChain->SetBranchAddress("truth_smeared_reco_vertex_fiducial", &truth_smeared_reco_vertex_fiducial, &b_truth_smeared_reco_vertex_fiducial);
   fChain->SetBranchAddress("truth_pass_nearInclusive", &truth_pass_nearInclusive, &b_truth_pass_nearInclusive);
   fChain->SetBranchAddress("truth_muon_fraction_accepted", &truth_muon_fraction_accepted, &b_truth_muon_fraction_accepted);
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
   fChain->SetBranchAddress("CCInclusiveReco_nuFlavor", &CCInclusiveReco_nuFlavor, &b_CCInclusiveReco_nuFlavor);
   fChain->SetBranchAddress("CCInclusiveReco_nuHelicity", &CCInclusiveReco_nuHelicity, &b_CCInclusiveReco_nuHelicity);
   fChain->SetBranchAddress("CCInclusiveReco_intCurrent", &CCInclusiveReco_intCurrent, &b_CCInclusiveReco_intCurrent);
   fChain->SetBranchAddress("CCInclusiveReco_intType", &CCInclusiveReco_intType, &b_CCInclusiveReco_intType);
   fChain->SetBranchAddress("CCInclusiveReco_E", &CCInclusiveReco_E, &b_CCInclusiveReco_E);
   fChain->SetBranchAddress("CCInclusiveReco_Q2", &CCInclusiveReco_Q2, &b_CCInclusiveReco_Q2);
   fChain->SetBranchAddress("CCInclusiveReco_x", &CCInclusiveReco_x, &b_CCInclusiveReco_x);
   fChain->SetBranchAddress("CCInclusiveReco_y", &CCInclusiveReco_y, &b_CCInclusiveReco_y);
   fChain->SetBranchAddress("CCInclusiveReco_W", &CCInclusiveReco_W, &b_CCInclusiveReco_W);
   fChain->SetBranchAddress("CCInclusiveReco_score", &CCInclusiveReco_score, &b_CCInclusiveReco_score);
   fChain->SetBranchAddress("CCInclusiveReco_leptonE", CCInclusiveReco_leptonE, &b_CCInclusiveReco_leptonE);
   fChain->SetBranchAddress("CCInclusiveReco_vtx", CCInclusiveReco_vtx, &b_CCInclusiveReco_vtx);
   fChain->SetBranchAddress("CCInclusiveReco_minos_trk_is_contained", &CCInclusiveReco_minos_trk_is_contained, &b_CCInclusiveReco_minos_trk_is_contained);
   fChain->SetBranchAddress("CCInclusiveReco_minos_trk_is_ok", &CCInclusiveReco_minos_trk_is_ok, &b_CCInclusiveReco_minos_trk_is_ok);
   fChain->SetBranchAddress("CCInclusiveReco_minos_used_range", &CCInclusiveReco_minos_used_range, &b_CCInclusiveReco_minos_used_range);
   fChain->SetBranchAddress("CCInclusiveReco_minos_used_curvature", &CCInclusiveReco_minos_used_curvature, &b_CCInclusiveReco_minos_used_curvature);
   fChain->SetBranchAddress("CCInclusiveReco_minos_trk_end_plane", &CCInclusiveReco_minos_trk_end_plane, &b_CCInclusiveReco_minos_trk_end_plane);
   fChain->SetBranchAddress("CCInclusiveReco_minos_trk_quality", &CCInclusiveReco_minos_trk_quality, &b_CCInclusiveReco_minos_trk_quality);
   fChain->SetBranchAddress("CCInclusiveReco_r_minos_trk_vtx_plane", &CCInclusiveReco_r_minos_trk_vtx_plane, &b_CCInclusiveReco_r_minos_trk_vtx_plane);
   fChain->SetBranchAddress("CCInclusiveReco_t_minos_trk_numFSMuons", &CCInclusiveReco_t_minos_trk_numFSMuons, &b_CCInclusiveReco_t_minos_trk_numFSMuons);
   fChain->SetBranchAddress("CCInclusiveReco_t_minos_trk_primFSLeptonPDG", &CCInclusiveReco_t_minos_trk_primFSLeptonPDG, &b_CCInclusiveReco_t_minos_trk_primFSLeptonPDG);
   fChain->SetBranchAddress("CCInclusiveReco_Q2_ccqe", &CCInclusiveReco_Q2_ccqe, &b_CCInclusiveReco_Q2_ccqe);
   fChain->SetBranchAddress("CCInclusiveReco_minos_trk_bave", &CCInclusiveReco_minos_trk_bave, &b_CCInclusiveReco_minos_trk_bave);
   fChain->SetBranchAddress("CCInclusiveReco_minos_trk_chi2", &CCInclusiveReco_minos_trk_chi2, &b_CCInclusiveReco_minos_trk_chi2);
   fChain->SetBranchAddress("CCInclusiveReco_minos_trk_end_u", &CCInclusiveReco_minos_trk_end_u, &b_CCInclusiveReco_minos_trk_end_u);
   fChain->SetBranchAddress("CCInclusiveReco_minos_trk_end_v", &CCInclusiveReco_minos_trk_end_v, &b_CCInclusiveReco_minos_trk_end_v);
   fChain->SetBranchAddress("CCInclusiveReco_minos_trk_end_x", &CCInclusiveReco_minos_trk_end_x, &b_CCInclusiveReco_minos_trk_end_x);
   fChain->SetBranchAddress("CCInclusiveReco_minos_trk_end_y", &CCInclusiveReco_minos_trk_end_y, &b_CCInclusiveReco_minos_trk_end_y);
   fChain->SetBranchAddress("CCInclusiveReco_minos_trk_end_z", &CCInclusiveReco_minos_trk_end_z, &b_CCInclusiveReco_minos_trk_end_z);
   fChain->SetBranchAddress("CCInclusiveReco_minos_trk_eqp", &CCInclusiveReco_minos_trk_eqp, &b_CCInclusiveReco_minos_trk_eqp);
   fChain->SetBranchAddress("CCInclusiveReco_minos_trk_eqp_qp", &CCInclusiveReco_minos_trk_eqp_qp, &b_CCInclusiveReco_minos_trk_eqp_qp);
   fChain->SetBranchAddress("CCInclusiveReco_minos_trk_fit_pass", &CCInclusiveReco_minos_trk_fit_pass, &b_CCInclusiveReco_minos_trk_fit_pass);
   fChain->SetBranchAddress("CCInclusiveReco_minos_trk_ndf", &CCInclusiveReco_minos_trk_ndf, &b_CCInclusiveReco_minos_trk_ndf);
   fChain->SetBranchAddress("CCInclusiveReco_minos_trk_p", &CCInclusiveReco_minos_trk_p, &b_CCInclusiveReco_minos_trk_p);
   fChain->SetBranchAddress("CCInclusiveReco_minos_trk_p_curvature", &CCInclusiveReco_minos_trk_p_curvature, &b_CCInclusiveReco_minos_trk_p_curvature);
   fChain->SetBranchAddress("CCInclusiveReco_minos_trk_p_range", &CCInclusiveReco_minos_trk_p_range, &b_CCInclusiveReco_minos_trk_p_range);
   fChain->SetBranchAddress("CCInclusiveReco_minos_trk_qp", &CCInclusiveReco_minos_trk_qp, &b_CCInclusiveReco_minos_trk_qp);
   fChain->SetBranchAddress("CCInclusiveReco_minos_trk_vtx_x", &CCInclusiveReco_minos_trk_vtx_x, &b_CCInclusiveReco_minos_trk_vtx_x);
   fChain->SetBranchAddress("CCInclusiveReco_minos_trk_vtx_y", &CCInclusiveReco_minos_trk_vtx_y, &b_CCInclusiveReco_minos_trk_vtx_y);
   fChain->SetBranchAddress("CCInclusiveReco_minos_trk_vtx_z", &CCInclusiveReco_minos_trk_vtx_z, &b_CCInclusiveReco_minos_trk_vtx_z);
   fChain->SetBranchAddress("CCInclusiveReco_nu_energy_recoil", &CCInclusiveReco_nu_energy_recoil, &b_CCInclusiveReco_nu_energy_recoil);
   fChain->SetBranchAddress("CCInclusiveReco_r_minos_trk_bdL", &CCInclusiveReco_r_minos_trk_bdL, &b_CCInclusiveReco_r_minos_trk_bdL);
   fChain->SetBranchAddress("CCInclusiveReco_r_minos_trk_end_dcosx", &CCInclusiveReco_r_minos_trk_end_dcosx, &b_CCInclusiveReco_r_minos_trk_end_dcosx);
   fChain->SetBranchAddress("CCInclusiveReco_r_minos_trk_end_dcosy", &CCInclusiveReco_r_minos_trk_end_dcosy, &b_CCInclusiveReco_r_minos_trk_end_dcosy);
   fChain->SetBranchAddress("CCInclusiveReco_r_minos_trk_end_dcosz", &CCInclusiveReco_r_minos_trk_end_dcosz, &b_CCInclusiveReco_r_minos_trk_end_dcosz);
   fChain->SetBranchAddress("CCInclusiveReco_r_minos_trk_vtx_dcosx", &CCInclusiveReco_r_minos_trk_vtx_dcosx, &b_CCInclusiveReco_r_minos_trk_vtx_dcosx);
   fChain->SetBranchAddress("CCInclusiveReco_r_minos_trk_vtx_dcosy", &CCInclusiveReco_r_minos_trk_vtx_dcosy, &b_CCInclusiveReco_r_minos_trk_vtx_dcosy);
   fChain->SetBranchAddress("CCInclusiveReco_r_minos_trk_vtx_dcosz", &CCInclusiveReco_r_minos_trk_vtx_dcosz, &b_CCInclusiveReco_r_minos_trk_vtx_dcosz);
   fChain->SetBranchAddress("CCInclusiveReco_t_minos_trk_primFSLepMinosInitProjPx", &CCInclusiveReco_t_minos_trk_primFSLepMinosInitProjPx, &b_CCInclusiveReco_t_minos_trk_primFSLepMinosInitProjPx);
   fChain->SetBranchAddress("CCInclusiveReco_t_minos_trk_primFSLepMinosInitProjPy", &CCInclusiveReco_t_minos_trk_primFSLepMinosInitProjPy, &b_CCInclusiveReco_t_minos_trk_primFSLepMinosInitProjPy);
   fChain->SetBranchAddress("CCInclusiveReco_t_minos_trk_primFSLepMinosInitProjPz", &CCInclusiveReco_t_minos_trk_primFSLepMinosInitProjPz, &b_CCInclusiveReco_t_minos_trk_primFSLepMinosInitProjPz);
   fChain->SetBranchAddress("CCInclusiveReco_t_minos_trk_primFSLepMinosInitProjX", &CCInclusiveReco_t_minos_trk_primFSLepMinosInitProjX, &b_CCInclusiveReco_t_minos_trk_primFSLepMinosInitProjX);
   fChain->SetBranchAddress("CCInclusiveReco_t_minos_trk_primFSLepMinosInitProjY", &CCInclusiveReco_t_minos_trk_primFSLepMinosInitProjY, &b_CCInclusiveReco_t_minos_trk_primFSLepMinosInitProjY);
   fChain->SetBranchAddress("CCInclusiveReco_t_minos_trk_primFSLepMinosInitProjZ", &CCInclusiveReco_t_minos_trk_primFSLepMinosInitProjZ, &b_CCInclusiveReco_t_minos_trk_primFSLepMinosInitProjZ);
   fChain->SetBranchAddress("CCInclusiveReco_t_minos_trk_primFSLepMnvFinalPx", &CCInclusiveReco_t_minos_trk_primFSLepMnvFinalPx, &b_CCInclusiveReco_t_minos_trk_primFSLepMnvFinalPx);
   fChain->SetBranchAddress("CCInclusiveReco_t_minos_trk_primFSLepMnvFinalPy", &CCInclusiveReco_t_minos_trk_primFSLepMnvFinalPy, &b_CCInclusiveReco_t_minos_trk_primFSLepMnvFinalPy);
   fChain->SetBranchAddress("CCInclusiveReco_t_minos_trk_primFSLepMnvFinalPz", &CCInclusiveReco_t_minos_trk_primFSLepMnvFinalPz, &b_CCInclusiveReco_t_minos_trk_primFSLepMnvFinalPz);
   fChain->SetBranchAddress("CCInclusiveReco_t_minos_trk_primFSLepMnvFinalX", &CCInclusiveReco_t_minos_trk_primFSLepMnvFinalX, &b_CCInclusiveReco_t_minos_trk_primFSLepMnvFinalX);
   fChain->SetBranchAddress("CCInclusiveReco_t_minos_trk_primFSLepMnvFinalY", &CCInclusiveReco_t_minos_trk_primFSLepMnvFinalY, &b_CCInclusiveReco_t_minos_trk_primFSLepMnvFinalY);
   fChain->SetBranchAddress("CCInclusiveReco_t_minos_trk_primFSLepMnvFinalZ", &CCInclusiveReco_t_minos_trk_primFSLepMnvFinalZ, &b_CCInclusiveReco_t_minos_trk_primFSLepMnvFinalZ);
   fChain->SetBranchAddress("CCInclusiveReco_t_minos_trk_primFSLepMnvInitPx", &CCInclusiveReco_t_minos_trk_primFSLepMnvInitPx, &b_CCInclusiveReco_t_minos_trk_primFSLepMnvInitPx);
   fChain->SetBranchAddress("CCInclusiveReco_t_minos_trk_primFSLepMnvInitPy", &CCInclusiveReco_t_minos_trk_primFSLepMnvInitPy, &b_CCInclusiveReco_t_minos_trk_primFSLepMnvInitPy);
   fChain->SetBranchAddress("CCInclusiveReco_t_minos_trk_primFSLepMnvInitPz", &CCInclusiveReco_t_minos_trk_primFSLepMnvInitPz, &b_CCInclusiveReco_t_minos_trk_primFSLepMnvInitPz);
   fChain->SetBranchAddress("CCInclusiveReco_t_minos_trk_primFSLepMnvInitX", &CCInclusiveReco_t_minos_trk_primFSLepMnvInitX, &b_CCInclusiveReco_t_minos_trk_primFSLepMnvInitX);
   fChain->SetBranchAddress("CCInclusiveReco_t_minos_trk_primFSLepMnvInitY", &CCInclusiveReco_t_minos_trk_primFSLepMnvInitY, &b_CCInclusiveReco_t_minos_trk_primFSLepMnvInitY);
   fChain->SetBranchAddress("CCInclusiveReco_t_minos_trk_primFSLepMnvInitZ", &CCInclusiveReco_t_minos_trk_primFSLepMnvInitZ, &b_CCInclusiveReco_t_minos_trk_primFSLepMnvInitZ);
   fChain->SetBranchAddress("CCInclusiveReco_sys_muon_energy_shift", CCInclusiveReco_sys_muon_energy_shift, &b_CCInclusiveReco_sys_muon_energy_shift);
   fChain->SetBranchAddress("CCInclusiveReco_sys_muon_qSquared_shift", CCInclusiveReco_sys_muon_qSquared_shift, &b_CCInclusiveReco_sys_muon_qSquared_shift);
   fChain->SetBranchAddress("CCInclusiveReco_sys_muon_wSquared_shift", CCInclusiveReco_sys_muon_wSquared_shift, &b_CCInclusiveReco_sys_muon_wSquared_shift);
   fChain->SetBranchAddress("CCInclusiveReco_sys_muon_xbj_shift", CCInclusiveReco_sys_muon_xbj_shift, &b_CCInclusiveReco_sys_muon_xbj_shift);
   fChain->SetBranchAddress("CCInclusiveReco_sys_muon_y_shift", CCInclusiveReco_sys_muon_y_shift, &b_CCInclusiveReco_sys_muon_y_shift);
   fChain->SetBranchAddress("CCInclusiveReco_sys_nu_energy_shift", CCInclusiveReco_sys_nu_energy_shift, &b_CCInclusiveReco_sys_nu_energy_shift);
   fChain->SetBranchAddress("CCInclusiveReco_sys_recoil_energy_shift", CCInclusiveReco_sys_recoil_energy_shift, &b_CCInclusiveReco_sys_recoil_energy_shift);
   fChain->SetBranchAddress("CCInclusiveReco_sys_recoil_qSquared_shift", CCInclusiveReco_sys_recoil_qSquared_shift, &b_CCInclusiveReco_sys_recoil_qSquared_shift);
   fChain->SetBranchAddress("CCInclusiveReco_sys_recoil_wSquared_shift", CCInclusiveReco_sys_recoil_wSquared_shift, &b_CCInclusiveReco_sys_recoil_wSquared_shift);
   fChain->SetBranchAddress("CCInclusiveReco_sys_recoil_xbj_shift", CCInclusiveReco_sys_recoil_xbj_shift, &b_CCInclusiveReco_sys_recoil_xbj_shift);
   fChain->SetBranchAddress("CCInclusiveReco_sys_recoil_y_shift", CCInclusiveReco_sys_recoil_y_shift, &b_CCInclusiveReco_sys_recoil_y_shift);
   fChain->SetBranchAddress("CCInclusiveReco_sys_total_qSquared_shift", CCInclusiveReco_sys_total_qSquared_shift, &b_CCInclusiveReco_sys_total_qSquared_shift);
   fChain->SetBranchAddress("CCInclusiveReco_sys_total_wSquared_shift", CCInclusiveReco_sys_total_wSquared_shift, &b_CCInclusiveReco_sys_total_wSquared_shift);
   fChain->SetBranchAddress("CCInclusiveReco_sys_total_xbj_shift", CCInclusiveReco_sys_total_xbj_shift, &b_CCInclusiveReco_sys_total_xbj_shift);
   fChain->SetBranchAddress("CCInclusiveReco_sys_total_y_shift", CCInclusiveReco_sys_total_y_shift, &b_CCInclusiveReco_sys_total_y_shift);
   fChain->SetBranchAddress("n_prongs", &n_prongs, &b_n_prongs);
   fChain->SetBranchAddress("prong_nParticles", prong_nParticles, &b_prong_nParticles);
   fChain->SetBranchAddress("prong_PDGCode", prong_PDGCode, &b_prong_PDGCode);
   fChain->SetBranchAddress("prong_EnergyWeightedMeanClusterTime", prong_EnergyWeightedMeanClusterTime, &b_prong_EnergyWeightedMeanClusterTime);
   fChain->SetBranchAddress("prong_VisibleEnergy", prong_VisibleEnergy, &b_prong_VisibleEnergy);
   fChain->SetBranchAddress("prong_nu_energy_dispersed_cal", prong_nu_energy_dispersed_cal, &b_prong_nu_energy_dispersed_cal);
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
   fChain->SetBranchAddress("prong_part_score", prong_part_score, &b_prong_part_score);
   fChain->SetBranchAddress("prong_part_mass", prong_part_mass, &b_prong_part_mass);
   fChain->SetBranchAddress("prong_part_charge", prong_part_charge, &b_prong_part_charge);
   fChain->SetBranchAddress("prong_part_pid", prong_part_pid, &b_prong_part_pid);
   fChain->SetBranchAddress("prong_part_E", &prong_part_E, &b_prong_part_E);
   fChain->SetBranchAddress("prong_part_pos", &prong_part_pos, &b_prong_part_pos);
   Notify();
}

Bool_t CCInclusive::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

    cout<<"Done!"<<endl;
   return kTRUE;
}

void CCInclusive::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

Int_t CCInclusive::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef CCInclusive_cxx

