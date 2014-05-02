/*
================================================================================
Function List: Default_Functions.cpp
    Contains auto-generated ROOT Specific Functions for
    CCDeltaPlus Class
    
    Functions Defined under CCDeltaPlus Namespace

    Usage:
        > #include "CCDeltaPlus.cpp"
        > #ifdef CCDeltaPlus_cxx 
    
    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision:  2014_05_01
================================================================================
*/

#include "CCDeltaPlus.cpp"
#ifdef CCDeltaPlus_cxx

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
   fChain->SetBranchAddress("isMinosMatchTrack", &isMinosMatchTrack, &b_isMinosMatchTrack);
   fChain->SetBranchAddress("isMinosMatchStub", &isMinosMatchStub, &b_isMinosMatchStub);
   fChain->SetBranchAddress("well_fit_vertex", &well_fit_vertex, &b_well_fit_vertex);
   fChain->SetBranchAddress("isBrokenTrack", &isBrokenTrack, &b_isBrokenTrack);
   fChain->SetBranchAddress("Cut_EndPoint_Michel_Exist", &Cut_EndPoint_Michel_Exist, &b_Cut_EndPoint_Michel_Exist);
   fChain->SetBranchAddress("Cut_Muon_None", &Cut_Muon_None, &b_Cut_Muon_None);
   fChain->SetBranchAddress("Cut_Muon_Score_Low", &Cut_Muon_Score_Low, &b_Cut_Muon_Score_Low);
   fChain->SetBranchAddress("Cut_Proton_None", &Cut_Proton_None, &b_Cut_Proton_None);
   fChain->SetBranchAddress("Cut_Vertex_Michel_Exist", &Cut_Vertex_Michel_Exist, &b_Cut_Vertex_Michel_Exist);
   fChain->SetBranchAddress("Cut_Vertex_None", &Cut_Vertex_None, &b_Cut_Vertex_None);
   fChain->SetBranchAddress("Cut_Vertex_Not_Analyzable", &Cut_Vertex_Not_Analyzable, &b_Cut_Vertex_Not_Analyzable);
   fChain->SetBranchAddress("Cut_Vertex_Not_Fiducial", &Cut_Vertex_Not_Fiducial, &b_Cut_Vertex_Not_Fiducial);
   fChain->SetBranchAddress("Cut_Vertex_Null", &Cut_Vertex_Null, &b_Cut_Vertex_Null);
   fChain->SetBranchAddress("Cut_secEndPoint_Michel_Exist", &Cut_secEndPoint_Michel_Exist, &b_Cut_secEndPoint_Michel_Exist);
   fChain->SetBranchAddress("broken_track_most_us_plane", &broken_track_most_us_plane, &b_broken_track_most_us_plane);
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
   fChain->SetBranchAddress("phys_energy_in_road_downstream_nplanes", &phys_energy_in_road_downstream_nplanes, &b_phys_energy_in_road_downstream_nplanes);
   fChain->SetBranchAddress("phys_energy_in_road_upstream_nplanes", &phys_energy_in_road_upstream_nplanes, &b_phys_energy_in_road_upstream_nplanes);
   fChain->SetBranchAddress("phys_n_dead_discr_pair", &phys_n_dead_discr_pair, &b_phys_n_dead_discr_pair);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_in_prim_track_region", &phys_n_dead_discr_pair_in_prim_track_region, &b_phys_n_dead_discr_pair_in_prim_track_region);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_two_mod_downstream_prim_track", &phys_n_dead_discr_pair_two_mod_downstream_prim_track, &b_phys_n_dead_discr_pair_two_mod_downstream_prim_track);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_two_mod_upstream_prim_vtx", &phys_n_dead_discr_pair_two_mod_upstream_prim_vtx, &b_phys_n_dead_discr_pair_two_mod_upstream_prim_vtx);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_upstream_prim_track_proj", &phys_n_dead_discr_pair_upstream_prim_track_proj, &b_phys_n_dead_discr_pair_upstream_prim_track_proj);
   fChain->SetBranchAddress("phys_vertex_is_fiducial", &phys_vertex_is_fiducial, &b_phys_vertex_is_fiducial);
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
   fChain->SetBranchAddress("unattachedExtraE", &unattachedExtraE, &b_unattachedExtraE);
   fChain->SetBranchAddress("vtxBlobExtraE", &vtxBlobExtraE, &b_vtxBlobExtraE);
   fChain->SetBranchAddress("vtx_michel_distance", &vtx_michel_distance, &b_vtx_michel_distance);
   fChain->SetBranchAddress("well_fit_vertex_angle", &well_fit_vertex_angle, &b_well_fit_vertex_angle);
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
   fChain->SetBranchAddress("truth_pi0_trackID", truth_pi0_trackID, &b_truth_pi0_trackID);
   fChain->SetBranchAddress("truth_proton_trackID", truth_proton_trackID, &b_truth_proton_trackID);
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
   fChain->SetBranchAddress("truth_pi0_E", truth_pi0_E, &b_truth_pi0_E);
   fChain->SetBranchAddress("truth_pi0_px", truth_pi0_px, &b_truth_pi0_px);
   fChain->SetBranchAddress("truth_pi0_py", truth_pi0_py, &b_truth_pi0_py);
   fChain->SetBranchAddress("truth_pi0_pz", truth_pi0_pz, &b_truth_pi0_pz);
   fChain->SetBranchAddress("truth_pi0_theta_wrtbeam", truth_pi0_theta_wrtbeam, &b_truth_pi0_theta_wrtbeam);
   fChain->SetBranchAddress("truth_proton_E", truth_proton_E, &b_truth_proton_E);
   fChain->SetBranchAddress("truth_proton_px", truth_proton_px, &b_truth_proton_px);
   fChain->SetBranchAddress("truth_proton_py", truth_proton_py, &b_truth_proton_py);
   fChain->SetBranchAddress("truth_proton_pz", truth_proton_pz, &b_truth_proton_pz);
   fChain->SetBranchAddress("truth_proton_theta_wrtbeam", truth_proton_theta_wrtbeam, &b_truth_proton_theta_wrtbeam);
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
   fChain->SetBranchAddress("CCDeltaPlusAna_isMuonInsideOD", &CCDeltaPlusAna_isMuonInsideOD, &b_CCDeltaPlusAna_isMuonInsideOD);
   fChain->SetBranchAddress("CCDeltaPlusAna_minos_trk_end_plane", &CCDeltaPlusAna_minos_trk_end_plane, &b_CCDeltaPlusAna_minos_trk_end_plane);
   fChain->SetBranchAddress("CCDeltaPlusAna_minos_trk_quality", &CCDeltaPlusAna_minos_trk_quality, &b_CCDeltaPlusAna_minos_trk_quality);
   fChain->SetBranchAddress("CCDeltaPlusAna_muon_N_minosTracks", &CCDeltaPlusAna_muon_N_minosTracks, &b_CCDeltaPlusAna_muon_N_minosTracks);
   fChain->SetBranchAddress("CCDeltaPlusAna_muon_minervaTrack_types", &CCDeltaPlusAna_muon_minervaTrack_types, &b_CCDeltaPlusAna_muon_minervaTrack_types);
   fChain->SetBranchAddress("CCDeltaPlusAna_muon_minosTrackQuality", &CCDeltaPlusAna_muon_minosTrackQuality, &b_CCDeltaPlusAna_muon_minosTrackQuality);
   fChain->SetBranchAddress("CCDeltaPlusAna_muon_roadUpstreamPlanes", &CCDeltaPlusAna_muon_roadUpstreamPlanes, &b_CCDeltaPlusAna_muon_roadUpstreamPlanes);
   fChain->SetBranchAddress("CCDeltaPlusAna_ntrajMuonProng", &CCDeltaPlusAna_ntrajMuonProng, &b_CCDeltaPlusAna_ntrajMuonProng);
   fChain->SetBranchAddress("CCDeltaPlusAna_r_minos_trk_vtx_plane", &CCDeltaPlusAna_r_minos_trk_vtx_plane, &b_CCDeltaPlusAna_r_minos_trk_vtx_plane);
   fChain->SetBranchAddress("CCDeltaPlusAna_t_minos_trk_numFSMuons", &CCDeltaPlusAna_t_minos_trk_numFSMuons, &b_CCDeltaPlusAna_t_minos_trk_numFSMuons);
   fChain->SetBranchAddress("CCDeltaPlusAna_t_minos_trk_primFSLeptonPDG", &CCDeltaPlusAna_t_minos_trk_primFSLeptonPDG, &b_CCDeltaPlusAna_t_minos_trk_primFSLeptonPDG);
   fChain->SetBranchAddress("CCDeltaPlusAna_trajMuonProngPDG", &CCDeltaPlusAna_trajMuonProngPDG, &b_CCDeltaPlusAna_trajMuonProngPDG);
   fChain->SetBranchAddress("CCDeltaPlusAna_trajMuonProngPrimary", &CCDeltaPlusAna_trajMuonProngPrimary, &b_CCDeltaPlusAna_trajMuonProngPrimary);
   fChain->SetBranchAddress("CCDeltaPlusAna_vtx_module", &CCDeltaPlusAna_vtx_module, &b_CCDeltaPlusAna_vtx_module);
   fChain->SetBranchAddress("CCDeltaPlusAna_vtx_plane", &CCDeltaPlusAna_vtx_plane, &b_CCDeltaPlusAna_vtx_plane);
   fChain->SetBranchAddress("CCDeltaPlusAna_endMuonTrajMomentum", &CCDeltaPlusAna_endMuonTrajMomentum, &b_CCDeltaPlusAna_endMuonTrajMomentum);
   fChain->SetBranchAddress("CCDeltaPlusAna_endMuonTrajXPosition", &CCDeltaPlusAna_endMuonTrajXPosition, &b_CCDeltaPlusAna_endMuonTrajXPosition);
   fChain->SetBranchAddress("CCDeltaPlusAna_endMuonTrajYPosition", &CCDeltaPlusAna_endMuonTrajYPosition, &b_CCDeltaPlusAna_endMuonTrajYPosition);
   fChain->SetBranchAddress("CCDeltaPlusAna_endMuonTrajZPosition", &CCDeltaPlusAna_endMuonTrajZPosition, &b_CCDeltaPlusAna_endMuonTrajZPosition);
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
   fChain->SetBranchAddress("CCDeltaPlusAna_muon_E", &CCDeltaPlusAna_muon_E, &b_CCDeltaPlusAna_muon_E);
   fChain->SetBranchAddress("CCDeltaPlusAna_muon_E_shift", &CCDeltaPlusAna_muon_E_shift, &b_CCDeltaPlusAna_muon_E_shift);
   fChain->SetBranchAddress("CCDeltaPlusAna_muon_muScore", &CCDeltaPlusAna_muon_muScore, &b_CCDeltaPlusAna_muon_muScore);
   fChain->SetBranchAddress("CCDeltaPlusAna_muon_p", &CCDeltaPlusAna_muon_p, &b_CCDeltaPlusAna_muon_p);
   fChain->SetBranchAddress("CCDeltaPlusAna_muon_px", &CCDeltaPlusAna_muon_px, &b_CCDeltaPlusAna_muon_px);
   fChain->SetBranchAddress("CCDeltaPlusAna_muon_py", &CCDeltaPlusAna_muon_py, &b_CCDeltaPlusAna_muon_py);
   fChain->SetBranchAddress("CCDeltaPlusAna_muon_pz", &CCDeltaPlusAna_muon_pz, &b_CCDeltaPlusAna_muon_pz);
   fChain->SetBranchAddress("CCDeltaPlusAna_muon_qp", &CCDeltaPlusAna_muon_qp, &b_CCDeltaPlusAna_muon_qp);
   fChain->SetBranchAddress("CCDeltaPlusAna_muon_qpqpe", &CCDeltaPlusAna_muon_qpqpe, &b_CCDeltaPlusAna_muon_qpqpe);
   fChain->SetBranchAddress("CCDeltaPlusAna_muon_roadUpstreamEnergy", &CCDeltaPlusAna_muon_roadUpstreamEnergy, &b_CCDeltaPlusAna_muon_roadUpstreamEnergy);
   fChain->SetBranchAddress("CCDeltaPlusAna_muon_theta", &CCDeltaPlusAna_muon_theta, &b_CCDeltaPlusAna_muon_theta);
   fChain->SetBranchAddress("CCDeltaPlusAna_muon_theta_biasDown", &CCDeltaPlusAna_muon_theta_biasDown, &b_CCDeltaPlusAna_muon_theta_biasDown);
   fChain->SetBranchAddress("CCDeltaPlusAna_muon_theta_biasUp", &CCDeltaPlusAna_muon_theta_biasUp, &b_CCDeltaPlusAna_muon_theta_biasUp);
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
   fChain->SetBranchAddress("CCDeltaPlusAna_trajMuonPhi", &CCDeltaPlusAna_trajMuonPhi, &b_CCDeltaPlusAna_trajMuonPhi);
   fChain->SetBranchAddress("CCDeltaPlusAna_trajMuonProngMomentum", &CCDeltaPlusAna_trajMuonProngMomentum, &b_CCDeltaPlusAna_trajMuonProngMomentum);
   fChain->SetBranchAddress("CCDeltaPlusAna_trajMuonProngPx", &CCDeltaPlusAna_trajMuonProngPx, &b_CCDeltaPlusAna_trajMuonProngPx);
   fChain->SetBranchAddress("CCDeltaPlusAna_trajMuonProngPy", &CCDeltaPlusAna_trajMuonProngPy, &b_CCDeltaPlusAna_trajMuonProngPy);
   fChain->SetBranchAddress("CCDeltaPlusAna_trajMuonProngPz", &CCDeltaPlusAna_trajMuonProngPz, &b_CCDeltaPlusAna_trajMuonProngPz);
   fChain->SetBranchAddress("CCDeltaPlusAna_trajMuonTheta", &CCDeltaPlusAna_trajMuonTheta, &b_CCDeltaPlusAna_trajMuonTheta);
   fChain->SetBranchAddress("CCDeltaPlusAna_vtx_x", &CCDeltaPlusAna_vtx_x, &b_CCDeltaPlusAna_vtx_x);
   fChain->SetBranchAddress("CCDeltaPlusAna_vtx_y", &CCDeltaPlusAna_vtx_y, &b_CCDeltaPlusAna_vtx_y);
   fChain->SetBranchAddress("CCDeltaPlusAna_vtx_z", &CCDeltaPlusAna_vtx_z, &b_CCDeltaPlusAna_vtx_z);
   fChain->SetBranchAddress("CCDeltaPlusAna_isProtonInsideOD", CCDeltaPlusAna_isProtonInsideOD, &b_CCDeltaPlusAna_isProtonInsideOD);
   fChain->SetBranchAddress("CCDeltaPlusAna_ntrajProtonProng", CCDeltaPlusAna_ntrajProtonProng, &b_CCDeltaPlusAna_ntrajProtonProng);
   fChain->SetBranchAddress("CCDeltaPlusAna_proton_kinked", CCDeltaPlusAna_proton_kinked, &b_CCDeltaPlusAna_proton_kinked);
   fChain->SetBranchAddress("CCDeltaPlusAna_proton_odMatch", CCDeltaPlusAna_proton_odMatch, &b_CCDeltaPlusAna_proton_odMatch);
   fChain->SetBranchAddress("CCDeltaPlusAna_proton_trk_pat_history", CCDeltaPlusAna_proton_trk_pat_history, &b_CCDeltaPlusAna_proton_trk_pat_history);
   fChain->SetBranchAddress("CCDeltaPlusAna_trajProtonProngPDG", CCDeltaPlusAna_trajProtonProngPDG, &b_CCDeltaPlusAna_trajProtonProngPDG);
   fChain->SetBranchAddress("CCDeltaPlusAna_trajProtonProngPrimary", CCDeltaPlusAna_trajProtonProngPrimary, &b_CCDeltaPlusAna_trajProtonProngPrimary);
   fChain->SetBranchAddress("CCDeltaPlusAna_endProtonTrajMomentum", CCDeltaPlusAna_endProtonTrajMomentum, &b_CCDeltaPlusAna_endProtonTrajMomentum);
   fChain->SetBranchAddress("CCDeltaPlusAna_endProtonTrajXPosition", CCDeltaPlusAna_endProtonTrajXPosition, &b_CCDeltaPlusAna_endProtonTrajXPosition);
   fChain->SetBranchAddress("CCDeltaPlusAna_endProtonTrajYPosition", CCDeltaPlusAna_endProtonTrajYPosition, &b_CCDeltaPlusAna_endProtonTrajYPosition);
   fChain->SetBranchAddress("CCDeltaPlusAna_endProtonTrajZPosition", CCDeltaPlusAna_endProtonTrajZPosition, &b_CCDeltaPlusAna_endProtonTrajZPosition);
   fChain->SetBranchAddress("CCDeltaPlusAna_proton_E", CCDeltaPlusAna_proton_E, &b_CCDeltaPlusAna_proton_E);
   fChain->SetBranchAddress("CCDeltaPlusAna_proton_chi2_ndf", CCDeltaPlusAna_proton_chi2_ndf, &b_CCDeltaPlusAna_proton_chi2_ndf);
   fChain->SetBranchAddress("CCDeltaPlusAna_proton_ekin", CCDeltaPlusAna_proton_ekin, &b_CCDeltaPlusAna_proton_ekin);
   fChain->SetBranchAddress("CCDeltaPlusAna_proton_endPointX", CCDeltaPlusAna_proton_endPointX, &b_CCDeltaPlusAna_proton_endPointX);
   fChain->SetBranchAddress("CCDeltaPlusAna_proton_endPointY", CCDeltaPlusAna_proton_endPointY, &b_CCDeltaPlusAna_proton_endPointY);
   fChain->SetBranchAddress("CCDeltaPlusAna_proton_endPointZ", CCDeltaPlusAna_proton_endPointZ, &b_CCDeltaPlusAna_proton_endPointZ);
   fChain->SetBranchAddress("CCDeltaPlusAna_proton_p", CCDeltaPlusAna_proton_p, &b_CCDeltaPlusAna_proton_p);
   fChain->SetBranchAddress("CCDeltaPlusAna_proton_p_calCorrection", CCDeltaPlusAna_proton_p_calCorrection, &b_CCDeltaPlusAna_proton_p_calCorrection);
   fChain->SetBranchAddress("CCDeltaPlusAna_proton_p_dEdXTool", CCDeltaPlusAna_proton_p_dEdXTool, &b_CCDeltaPlusAna_proton_p_dEdXTool);
   fChain->SetBranchAddress("CCDeltaPlusAna_proton_p_visEnergy", CCDeltaPlusAna_proton_p_visEnergy, &b_CCDeltaPlusAna_proton_p_visEnergy);
   fChain->SetBranchAddress("CCDeltaPlusAna_proton_phi", CCDeltaPlusAna_proton_phi, &b_CCDeltaPlusAna_proton_phi);
   fChain->SetBranchAddress("CCDeltaPlusAna_proton_px", CCDeltaPlusAna_proton_px, &b_CCDeltaPlusAna_proton_px);
   fChain->SetBranchAddress("CCDeltaPlusAna_proton_py", CCDeltaPlusAna_proton_py, &b_CCDeltaPlusAna_proton_py);
   fChain->SetBranchAddress("CCDeltaPlusAna_proton_pz", CCDeltaPlusAna_proton_pz, &b_CCDeltaPlusAna_proton_pz);
   fChain->SetBranchAddress("CCDeltaPlusAna_proton_score", CCDeltaPlusAna_proton_score, &b_CCDeltaPlusAna_proton_score);
   fChain->SetBranchAddress("CCDeltaPlusAna_proton_score1", CCDeltaPlusAna_proton_score1, &b_CCDeltaPlusAna_proton_score1);
   fChain->SetBranchAddress("CCDeltaPlusAna_proton_score2", CCDeltaPlusAna_proton_score2, &b_CCDeltaPlusAna_proton_score2);
   fChain->SetBranchAddress("CCDeltaPlusAna_proton_startPointX", CCDeltaPlusAna_proton_startPointX, &b_CCDeltaPlusAna_proton_startPointX);
   fChain->SetBranchAddress("CCDeltaPlusAna_proton_startPointY", CCDeltaPlusAna_proton_startPointY, &b_CCDeltaPlusAna_proton_startPointY);
   fChain->SetBranchAddress("CCDeltaPlusAna_proton_startPointZ", CCDeltaPlusAna_proton_startPointZ, &b_CCDeltaPlusAna_proton_startPointZ);
   fChain->SetBranchAddress("CCDeltaPlusAna_proton_theta", CCDeltaPlusAna_proton_theta, &b_CCDeltaPlusAna_proton_theta);
   fChain->SetBranchAddress("CCDeltaPlusAna_proton_thetaX", CCDeltaPlusAna_proton_thetaX, &b_CCDeltaPlusAna_proton_thetaX);
   fChain->SetBranchAddress("CCDeltaPlusAna_proton_thetaY", CCDeltaPlusAna_proton_thetaY, &b_CCDeltaPlusAna_proton_thetaY);
   fChain->SetBranchAddress("CCDeltaPlusAna_trajProtonPhi", CCDeltaPlusAna_trajProtonPhi, &b_CCDeltaPlusAna_trajProtonPhi);
   fChain->SetBranchAddress("CCDeltaPlusAna_trajProtonProngMomentum", CCDeltaPlusAna_trajProtonProngMomentum, &b_CCDeltaPlusAna_trajProtonProngMomentum);
   fChain->SetBranchAddress("CCDeltaPlusAna_trajProtonProngPx", CCDeltaPlusAna_trajProtonProngPx, &b_CCDeltaPlusAna_trajProtonProngPx);
   fChain->SetBranchAddress("CCDeltaPlusAna_trajProtonProngPy", CCDeltaPlusAna_trajProtonProngPy, &b_CCDeltaPlusAna_trajProtonProngPy);
   fChain->SetBranchAddress("CCDeltaPlusAna_trajProtonProngPz", CCDeltaPlusAna_trajProtonProngPz, &b_CCDeltaPlusAna_trajProtonProngPz);
   fChain->SetBranchAddress("CCDeltaPlusAna_trajProtonTheta", CCDeltaPlusAna_trajProtonTheta, &b_CCDeltaPlusAna_trajProtonTheta);
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