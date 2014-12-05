//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Nov 25 09:30:20 2014 by ROOT version 5.34/05
// from TChain CCProtonPi0/
//////////////////////////////////////////////////////////

#ifndef MC_Sample_h
#define MC_Sample_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "/grid/fermiapp/minerva/software_releases/lcgcmake/build/lcg_61/projects/ROOT-5.34.05/src/ROOT/5.34.05/cint/cint/lib/prec_stl/vector"

// Fixed size dimensions of array or collections stored in the TTree if any.

class MC_Sample {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        eventID;
   Int_t           physEvtNum;
   Int_t           n_hyps;
   Int_t           processType;
   Int_t           primaryPart;
   Int_t           n_slices;
   Int_t           slice_numbers[1];   //[n_slices]
   Int_t           shared_slice;
   Double_t        vtx[4];
   Double_t        vtxErr[4];
   Double_t        E[4];
   Bool_t          found_truth;
   Bool_t          phys_front_activity;
   Bool_t          phys_energy_in_road_upstream_is_rockmuon_consistent;
   Bool_t          rock_muons_removed;
   Bool_t          minos_track_match;
   Bool_t          minos_stub_match;
   Bool_t          unknown_helicity;
   Bool_t          minos_track_inside_partial_plane;
   Bool_t          prim_vtx_has_misassigned_track_direction;
   Bool_t          prim_vtx_has_broken_track;
   Bool_t          gamma1_isGoodDirection;
   Bool_t          gamma1_isGoodPosition;
   Bool_t          gamma1_isGoodBlob;
   Bool_t          gamma2_isGoodDirection;
   Bool_t          gamma2_isGoodPosition;
   Bool_t          gamma2_isGoodBlob;
   Int_t           Cut_ConeBlobs;
   Int_t           Cut_EndPoint_Michel_Exist;
   Int_t           Cut_Event_Has_BadObject;
   Int_t           Cut_Muon_Charge;
   Int_t           Cut_Muon_None;
   Int_t           Cut_Muon_Not_Plausible;
   Int_t           Cut_Muon_Score_Low;
   Int_t           Cut_Particle_None;
   Int_t           Cut_PreFilter_Pi0;
   Int_t           Cut_Proton_None;
   Int_t           Cut_Vertex_Michel_Exist;
   Int_t           Cut_Vertex_None;
   Int_t           Cut_Vertex_Not_Fiducial;
   Int_t           Cut_Vertex_Not_Reconstructable;
   Int_t           Cut_Vertex_Null;
   Int_t           Cut_VtxBlob;
   Int_t           Cut_nProngs;
   Int_t           Cut_secEndPoint_Michel_Exist;
   Int_t           anglescan_ncand;
   Int_t           anglescan_ncandx;
   Int_t           broken_track_most_us_plane;
   Int_t           g1dedx_doublet;
   Int_t           g1dedx_empty_plane;
   Int_t           g1dedx_nplane;
   Int_t           g2dedx_doublet;
   Int_t           g2dedx_empty_plane;
   Int_t           g2dedx_nplane;
   Int_t           gamma1_blob_nclusters;
   Int_t           gamma1_blob_ndigits;
   Int_t           gamma2_blob_nclusters;
   Int_t           gamma2_blob_ndigits;
   Int_t           nProngs;
   Int_t           n_anchored_long_trk_prongs;
   Int_t           n_anchored_short_trk_prongs;
   Int_t           n_iso_trk_prongs;
   Int_t           n_long_tracks;
   Int_t           n_short_tracks;
   Int_t           n_vtx_michel_views;
   Int_t           nblob_anglescan;
   Int_t           nblob_hough;
   Int_t           od_energeticTower;
   Int_t           phys_energy_in_road_downstream_nplanes;
   Int_t           phys_energy_in_road_upstream_nplanes;
   Int_t           phys_n_dead_discr_pair;
   Int_t           phys_n_dead_discr_pair_in_prim_track_region;
   Int_t           phys_n_dead_discr_pair_two_mod_downstream_prim_track;
   Int_t           phys_n_dead_discr_pair_two_mod_upstream_prim_vtx;
   Int_t           phys_n_dead_discr_pair_upstream_prim_track_proj;
   Int_t           phys_vertex_is_fiducial;
   Int_t           preFilter_Result;
   Int_t           vtx_primary_index;
   Int_t           vtx_primary_multiplicity;
   Int_t           vtx_secondary_count;
   Int_t           vtx_total_count;
   Double_t        Filament_Vertex_energy;
   Double_t        RE_energy_ECAL;
   Double_t        RE_energy_HCAL;
   Double_t        RE_energy_Tracker;
   Double_t        Sphere_Vertex_energy;
   Double_t        Vertex_blob_energy;
   Double_t        dispersedExtraE;
   Double_t        energy_from_mc;
   Double_t        energy_from_mc_fraction;
   Double_t        energy_from_mc_fraction_of_highest;
   Double_t        evis_ecal;
   Double_t        evis_ecal_u;
   Double_t        evis_ecal_v;
   Double_t        evis_ecal_x;
   Double_t        evis_hcal;
   Double_t        evis_hcal_u;
   Double_t        evis_hcal_v;
   Double_t        evis_hcal_x;
   Double_t        evis_nearvtx_total;
   Double_t        evis_nearvtx_u;
   Double_t        evis_nearvtx_v;
   Double_t        evis_nearvtx_x;
   Double_t        evis_ntgt;
   Double_t        evis_ntgt_u;
   Double_t        evis_ntgt_v;
   Double_t        evis_ntgt_x;
   Double_t        evis_other;
   Double_t        evis_total;
   Double_t        evis_total_u;
   Double_t        evis_total_v;
   Double_t        evis_total_x;
   Double_t        evis_trkr;
   Double_t        evis_trkr_u;
   Double_t        evis_trkr_v;
   Double_t        evis_trkr_x;
   Double_t        g1blob_minsep;
   Double_t        g1dedx;
   Double_t        g1dedx1;
   Double_t        g1dedx_total;
   Double_t        g1dedx_total1;
   Double_t        g2blob_minsep;
   Double_t        g2dedx;
   Double_t        g2dedx1;
   Double_t        g2dedx_total;
   Double_t        g2dedx_total1;
   Double_t        gamma1_E;
   Double_t        gamma1_dEdx;
   Double_t        gamma1_dist_exit;
   Double_t        gamma1_dist_vtx;
   Double_t        gamma1_evis_ecal;
   Double_t        gamma1_evis_hcal;
   Double_t        gamma1_evis_scal;
   Double_t        gamma1_evis_trkr;
   Double_t        gamma1_phi;
   Double_t        gamma1_px;
   Double_t        gamma1_py;
   Double_t        gamma1_pz;
   Double_t        gamma1_theta;
   Double_t        gamma1_time;
   Double_t        gamma2_E;
   Double_t        gamma2_dEdx;
   Double_t        gamma2_dist_exit;
   Double_t        gamma2_dist_vtx;
   Double_t        gamma2_evis_ecal;
   Double_t        gamma2_evis_hcal;
   Double_t        gamma2_evis_scal;
   Double_t        gamma2_evis_trkr;
   Double_t        gamma2_phi;
   Double_t        gamma2_px;
   Double_t        gamma2_py;
   Double_t        gamma2_pz;
   Double_t        gamma2_theta;
   Double_t        gamma2_time;
   Double_t        hadronVisibleE;
   Double_t        muonVisibleE;
   Double_t        muon_phi;
   Double_t        muon_theta;
   Double_t        muon_thetaX;
   Double_t        muon_thetaY;
   Double_t        od_downstreamFrame;
   Double_t        od_downstreamFrame_z;
   Double_t        od_highStory;
   Double_t        od_highStory_t;
   Double_t        od_lowStory;
   Double_t        od_lowStory_t;
   Double_t        od_maxEnergy;
   Double_t        od_upstreamFrame;
   Double_t        od_upstreamFrame_z;
   Double_t        phys_energy_dispersed;
   Double_t        phys_energy_in_road_downstream;
   Double_t        phys_energy_in_road_upstream;
   Double_t        phys_energy_unattached;
   Double_t        pi0_E;
   Double_t        pi0_cos_openingAngle;
   Double_t        pi0_invMass;
   Double_t        pi0_openingAngle;
   Double_t        pi0_phi;
   Double_t        pi0_px;
   Double_t        pi0_py;
   Double_t        pi0_pz;
   Double_t        pi0_theta;
   Double_t        pi0_thetaX;
   Double_t        pi0_thetaY;
   Double_t        preFilter_rejectedEnergy;
   Double_t        prim_vtx_smallest_opening_angle;
   Double_t        time;
   Double_t        totalIDVisibleE;
   Double_t        totalODVisibleE;
   Double_t        totalVisibleE;
   Double_t        unattachedExtraE;
   Double_t        vtxBlobExtraE;
   Double_t        vtx_michel_distance;
   Int_t           anglescan_blob_nc_sz;
   Int_t           anglescan_blob_nc[11];   //[anglescan_blob_nc_sz]
   Int_t           anglescan_blob_ncu_sz;
   Int_t           anglescan_blob_ncu[11];   //[anglescan_blob_ncu_sz]
   Int_t           anglescan_blob_ncv_sz;
   Int_t           anglescan_blob_ncv[11];   //[anglescan_blob_ncv_sz]
   Int_t           anglescan_blob_ncx_sz;
   Int_t           anglescan_blob_ncx[11];   //[anglescan_blob_ncx_sz]
   Int_t           anglescan_blob_nd_sz;
   Int_t           anglescan_blob_nd[11];   //[anglescan_blob_nd_sz]
   Int_t           anglescan_blob_ndu_sz;
   Int_t           anglescan_blob_ndu[11];   //[anglescan_blob_ndu_sz]
   Int_t           anglescan_blob_ndv_sz;
   Int_t           anglescan_blob_ndv[11];   //[anglescan_blob_ndv_sz]
   Int_t           anglescan_blob_ndx_sz;
   Int_t           anglescan_blob_ndx[11];   //[anglescan_blob_ndx_sz]
   Int_t           anglescan_cand_nc_sz;
   Int_t           anglescan_cand_nc[11];   //[anglescan_cand_nc_sz]
   Int_t           anglescan_cand_ncu_sz;
   Int_t           anglescan_cand_ncu[11];   //[anglescan_cand_ncu_sz]
   Int_t           anglescan_cand_ncv_sz;
   Int_t           anglescan_cand_ncv[11];   //[anglescan_cand_ncv_sz]
   Int_t           anglescan_cand_ncx_sz;
   Int_t           anglescan_cand_ncx[11];   //[anglescan_cand_ncx_sz]
   Int_t           anglescan_cand_nd_sz;
   Int_t           anglescan_cand_nd[11];   //[anglescan_cand_nd_sz]
   Int_t           anglescan_cand_ndu_sz;
   Int_t           anglescan_cand_ndu[11];   //[anglescan_cand_ndu_sz]
   Int_t           anglescan_cand_ndv_sz;
   Int_t           anglescan_cand_ndv[11];   //[anglescan_cand_ndv_sz]
   Int_t           anglescan_cand_ndx_sz;
   Int_t           anglescan_cand_ndx[11];   //[anglescan_cand_ndx_sz]
   Int_t           anglescan_candx_nc_sz;
   Int_t           anglescan_candx_nc[11];   //[anglescan_candx_nc_sz]
   Int_t           anglescan_candx_nd_sz;
   Int_t           anglescan_candx_nd[11];   //[anglescan_candx_nd_sz]
   Int_t           final_blob_nc_sz;
   Int_t           final_blob_nc[2];   //[final_blob_nc_sz]
   Int_t           final_blob_ncu_sz;
   Int_t           final_blob_ncu[2];   //[final_blob_ncu_sz]
   Int_t           final_blob_ncv_sz;
   Int_t           final_blob_ncv[2];   //[final_blob_ncv_sz]
   Int_t           final_blob_ncx_sz;
   Int_t           final_blob_ncx[2];   //[final_blob_ncx_sz]
   Int_t           final_blob_nd_sz;
   Int_t           final_blob_nd[2];   //[final_blob_nd_sz]
   Int_t           final_blob_ndu_sz;
   Int_t           final_blob_ndu[2];   //[final_blob_ndu_sz]
   Int_t           final_blob_ndv_sz;
   Int_t           final_blob_ndv[2];   //[final_blob_ndv_sz]
   Int_t           final_blob_ndx_sz;
   Int_t           final_blob_ndx[2];   //[final_blob_ndx_sz]
   Int_t           g1dedx_cluster_occupancy_sz;
   Int_t           g1dedx_cluster_occupancy[6];   //[g1dedx_cluster_occupancy_sz]
   Int_t           g2dedx_cluster_occupancy_sz;
   Int_t           g2dedx_cluster_occupancy[6];   //[g2dedx_cluster_occupancy_sz]
   Int_t           hough_blob_nc_sz;
   Int_t           hough_blob_nc[3];   //[hough_blob_nc_sz]
   Int_t           hough_blob_ncu_sz;
   Int_t           hough_blob_ncu[3];   //[hough_blob_ncu_sz]
   Int_t           hough_blob_ncv_sz;
   Int_t           hough_blob_ncv[3];   //[hough_blob_ncv_sz]
   Int_t           hough_blob_ncx_sz;
   Int_t           hough_blob_ncx[3];   //[hough_blob_ncx_sz]
   Int_t           hough_blob_nd_sz;
   Int_t           hough_blob_nd[3];   //[hough_blob_nd_sz]
   Int_t           hough_blob_ndu_sz;
   Int_t           hough_blob_ndu[3];   //[hough_blob_ndu_sz]
   Int_t           hough_blob_ndv_sz;
   Int_t           hough_blob_ndv[3];   //[hough_blob_ndv_sz]
   Int_t           hough_blob_ndx_sz;
   Int_t           hough_blob_ndx[3];   //[hough_blob_ndx_sz]
   Int_t           Vertex_energy_radii_sz;
   Double_t        Vertex_energy_radii[7];   //[Vertex_energy_radii_sz]
   Int_t           blob_cluster_energy1_sz;
   Double_t        blob_cluster_energy1[6];   //[blob_cluster_energy1_sz]
   Int_t           blob_cluster_energy2_sz;
   Double_t        blob_cluster_energy2[6];   //[blob_cluster_energy2_sz]
   Int_t           g1dedx_cluster_energy_sz;
   Double_t        g1dedx_cluster_energy[6];   //[g1dedx_cluster_energy_sz]
   Int_t           g1dedx_rev_cluster_energy_sz;
   Double_t        g1dedx_rev_cluster_energy[44];   //[g1dedx_rev_cluster_energy_sz]
   Int_t           g2dedx_cluster_energy_sz;
   Double_t        g2dedx_cluster_energy[6];   //[g2dedx_cluster_energy_sz]
   Int_t           g2dedx_rev_cluster_energy_sz;
   Double_t        g2dedx_rev_cluster_energy[46];   //[g2dedx_rev_cluster_energy_sz]
   Double_t        gamma1_direction[3];
   Double_t        gamma1_vertex[3];
   Double_t        gamma2_direction[3];
   Double_t        gamma2_vertex[3];
   Int_t           good_mgg_vector_sz;
   Double_t        good_mgg_vector[1];   //[good_mgg_vector_sz]
   Int_t           mgg_vector_sz;
   Double_t        mgg_vector[1];   //[mgg_vector_sz]
   Int_t           od_distanceBlobTower_sz;
   Double_t        od_distanceBlobTower[2];   //[od_distanceBlobTower_sz]
   Int_t           od_idBlobTime_sz;
   Double_t        od_idBlobTime[2];   //[od_idBlobTime_sz]
   Int_t           od_towerEnergy_sz;
   Double_t        od_towerEnergy[6];   //[od_towerEnergy_sz]
   Int_t           od_towerNClusters_sz;
   Double_t        od_towerNClusters[6];   //[od_towerNClusters_sz]
   Int_t           od_towerTime_sz;
   Double_t        od_towerTime[6];   //[od_towerTime_sz]
   Int_t           od_towerTimeBlobMuon_sz;
   Double_t        od_towerTimeBlobMuon[2];   //[od_towerTimeBlobMuon_sz]
   Int_t           od_towerTimeBlobOD_sz;
   Double_t        od_towerTimeBlobOD[2];   //[od_towerTimeBlobOD_sz]
   Bool_t          truth_has_physics_event;
   Bool_t          truth_isSignal;
   Bool_t          truth_isSignalGold;
   Bool_t          truth_isSignalSilver1;
   Bool_t          truth_isSignal_1Pi0;
   Bool_t          truth_isSignal_2Gamma;
   Bool_t          truth_isFidVol;
   Int_t           truth_N_FSParticles;
   Int_t           truth_N_gamma;
   Int_t           truth_N_pi0;
   Int_t           truth_N_proton;
   Int_t           truth_N_proton_Long;
   Int_t           truth_N_proton_Short;
   Int_t           truth_muon_charge;
   Int_t           truth_pi0_GrandMother;
   Int_t           truth_pi0_GrandMotherStatus;
   Int_t           truth_pi0_Mother;
   Int_t           truth_pi0_MotherStatus;
   Int_t           truth_pi0_status;
   Int_t           truth_target_material;
   Int_t           truth_vertex_module;
   Int_t           truth_vertex_plane;
   Double_t        truth_muon_E;
   Double_t        truth_muon_px;
   Double_t        truth_muon_py;
   Double_t        truth_muon_pz;
   Double_t        truth_muon_theta_wrtbeam;
   Double_t        truth_pi0_E;
   Double_t        truth_pi0_px;
   Double_t        truth_pi0_py;
   Double_t        truth_pi0_pz;
   Double_t        truth_pi0_theta_wrtbeam;
   Double_t        truth_gamma_E[2];
   Double_t        truth_gamma_px[2];
   Double_t        truth_gamma_py[2];
   Double_t        truth_gamma_pz[2];
   Double_t        truth_gamma_theta_wrtbeam[2];
   Int_t           genie_wgt_n_shifts;
   Double_t        truth_genie_wgt_AGKYxF1pi[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_AhtBY[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_BhtBY[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_CCQEPauliSupViaKF[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_CV1uBY[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_CV2uBY[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_EtaNCEL[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_FrAbs_N[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_FrAbs_pi[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_FrCEx_N[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_FrCEx_pi[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_FrElas_N[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_FrElas_pi[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_FrInel_N[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_FrInel_pi[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_FrPiProd_N[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_FrPiProd_pi[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_MFP_N[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_MFP_pi[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_MaCCQE[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_MaCCQEshape[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_MaNCEL[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_MaRES[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_MvRES[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_NormCCQE[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_NormCCRES[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_NormDISCC[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_NormNCRES[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_RDecBR1gamma[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_Rvn1pi[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_Rvn2pi[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_Rvp1pi[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_Rvp2pi[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_Theta_Delta2Npi[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_VecFFCCQEshape[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_shifts[7];   //[genie_wgt_n_shifts]
   Double_t        truth_proton_E[10];
   Double_t        truth_proton_px[10];
   Double_t        truth_proton_py[10];
   Double_t        truth_proton_pz[10];
   Double_t        truth_proton_theta_wrtbeam[10];
   Int_t           CCProtonPi0_nuFlavor;
   Int_t           CCProtonPi0_nuHelicity;
   Int_t           CCProtonPi0_intCurrent;
   Int_t           CCProtonPi0_intType;
   Double_t        CCProtonPi0_E;
   Double_t        CCProtonPi0_Q2;
   Double_t        CCProtonPi0_x;
   Double_t        CCProtonPi0_y;
   Double_t        CCProtonPi0_W;
   Double_t        CCProtonPi0_score;
   Double_t        CCProtonPi0_leptonE[4];
   Double_t        CCProtonPi0_vtx[4];
   Bool_t          CCProtonPi0_minos_trk_is_contained;
   Bool_t          CCProtonPi0_minos_trk_is_ok;
   Bool_t          CCProtonPi0_minos_used_range;
   Bool_t          CCProtonPi0_minos_used_curvature;
   Int_t           CCProtonPi0_isMuonInsideOD;
   Int_t           CCProtonPi0_minos_trk_end_plane;
   Int_t           CCProtonPi0_minos_trk_quality;
   Int_t           CCProtonPi0_muon_N_minosTracks;
   Int_t           CCProtonPi0_muon_charge;
   Int_t           CCProtonPi0_muon_hasMinosMatchStub;
   Int_t           CCProtonPi0_muon_hasMinosMatchTrack;
   Int_t           CCProtonPi0_muon_minervaTrack_types;
   Int_t           CCProtonPi0_muon_minosTrackQuality;
   Int_t           CCProtonPi0_muon_roadUpstreamPlanes;
   Int_t           CCProtonPi0_ntrajMuonProng;
   Int_t           CCProtonPi0_r_minos_trk_vtx_plane;
   Int_t           CCProtonPi0_t_minos_trk_numFSMuons;
   Int_t           CCProtonPi0_t_minos_trk_primFSLeptonPDG;
   Int_t           CCProtonPi0_trajMuonProngPDG;
   Int_t           CCProtonPi0_trajMuonProngPrimary;
   Int_t           CCProtonPi0_vtx_module;
   Int_t           CCProtonPi0_vtx_plane;
   Double_t        CCProtonPi0_QSq;
   Double_t        CCProtonPi0_WSq;
   Double_t        CCProtonPi0_endMuonTrajMomentum;
   Double_t        CCProtonPi0_endMuonTrajXPosition;
   Double_t        CCProtonPi0_endMuonTrajYPosition;
   Double_t        CCProtonPi0_endMuonTrajZPosition;
   Double_t        CCProtonPi0_minos_trk_bave;
   Double_t        CCProtonPi0_minos_trk_chi2;
   Double_t        CCProtonPi0_minos_trk_end_u;
   Double_t        CCProtonPi0_minos_trk_end_v;
   Double_t        CCProtonPi0_minos_trk_end_x;
   Double_t        CCProtonPi0_minos_trk_end_y;
   Double_t        CCProtonPi0_minos_trk_end_z;
   Double_t        CCProtonPi0_minos_trk_eqp;
   Double_t        CCProtonPi0_minos_trk_eqp_qp;
   Double_t        CCProtonPi0_minos_trk_fit_pass;
   Double_t        CCProtonPi0_minos_trk_ndf;
   Double_t        CCProtonPi0_minos_trk_p;
   Double_t        CCProtonPi0_minos_trk_p_curvature;
   Double_t        CCProtonPi0_minos_trk_p_range;
   Double_t        CCProtonPi0_minos_trk_qp;
   Double_t        CCProtonPi0_minos_trk_vtx_x;
   Double_t        CCProtonPi0_minos_trk_vtx_y;
   Double_t        CCProtonPi0_minos_trk_vtx_z;
   Double_t        CCProtonPi0_muon_E;
   Double_t        CCProtonPi0_muon_E_shift;
   Double_t        CCProtonPi0_muon_muScore;
   Double_t        CCProtonPi0_muon_p;
   Double_t        CCProtonPi0_muon_px;
   Double_t        CCProtonPi0_muon_py;
   Double_t        CCProtonPi0_muon_pz;
   Double_t        CCProtonPi0_muon_qp;
   Double_t        CCProtonPi0_muon_qpqpe;
   Double_t        CCProtonPi0_muon_roadUpstreamEnergy;
   Double_t        CCProtonPi0_muon_theta;
   Double_t        CCProtonPi0_muon_theta_biasDown;
   Double_t        CCProtonPi0_muon_theta_biasUp;
   Double_t        CCProtonPi0_neutrino_E;
   Double_t        CCProtonPi0_neutrino_E_Cal;
   Double_t        CCProtonPi0_r_minos_trk_bdL;
   Double_t        CCProtonPi0_r_minos_trk_end_dcosx;
   Double_t        CCProtonPi0_r_minos_trk_end_dcosy;
   Double_t        CCProtonPi0_r_minos_trk_end_dcosz;
   Double_t        CCProtonPi0_r_minos_trk_vtx_dcosx;
   Double_t        CCProtonPi0_r_minos_trk_vtx_dcosy;
   Double_t        CCProtonPi0_r_minos_trk_vtx_dcosz;
   Double_t        CCProtonPi0_t_minos_trk_primFSLepMinosInitProjPx;
   Double_t        CCProtonPi0_t_minos_trk_primFSLepMinosInitProjPy;
   Double_t        CCProtonPi0_t_minos_trk_primFSLepMinosInitProjPz;
   Double_t        CCProtonPi0_t_minos_trk_primFSLepMinosInitProjX;
   Double_t        CCProtonPi0_t_minos_trk_primFSLepMinosInitProjY;
   Double_t        CCProtonPi0_t_minos_trk_primFSLepMinosInitProjZ;
   Double_t        CCProtonPi0_t_minos_trk_primFSLepMnvFinalPx;
   Double_t        CCProtonPi0_t_minos_trk_primFSLepMnvFinalPy;
   Double_t        CCProtonPi0_t_minos_trk_primFSLepMnvFinalPz;
   Double_t        CCProtonPi0_t_minos_trk_primFSLepMnvFinalX;
   Double_t        CCProtonPi0_t_minos_trk_primFSLepMnvFinalY;
   Double_t        CCProtonPi0_t_minos_trk_primFSLepMnvFinalZ;
   Double_t        CCProtonPi0_t_minos_trk_primFSLepMnvInitPx;
   Double_t        CCProtonPi0_t_minos_trk_primFSLepMnvInitPy;
   Double_t        CCProtonPi0_t_minos_trk_primFSLepMnvInitPz;
   Double_t        CCProtonPi0_t_minos_trk_primFSLepMnvInitX;
   Double_t        CCProtonPi0_t_minos_trk_primFSLepMnvInitY;
   Double_t        CCProtonPi0_t_minos_trk_primFSLepMnvInitZ;
   Double_t        CCProtonPi0_trajMuonPhi;
   Double_t        CCProtonPi0_trajMuonProngEnergy;
   Double_t        CCProtonPi0_trajMuonProngMomentum;
   Double_t        CCProtonPi0_trajMuonProngPSelf;
   Double_t        CCProtonPi0_trajMuonProngPx;
   Double_t        CCProtonPi0_trajMuonProngPy;
   Double_t        CCProtonPi0_trajMuonProngPz;
   Double_t        CCProtonPi0_trajMuonTheta;
   Double_t        CCProtonPi0_vtx_x;
   Double_t        CCProtonPi0_vtx_y;
   Double_t        CCProtonPi0_vtx_z;
   Int_t           CCProtonPi0_isProtonInsideOD[10];
   Int_t           CCProtonPi0_ntrajProtonProng[10];
   Int_t           CCProtonPi0_proton_isRecoGood[10];
   Int_t           CCProtonPi0_proton_kinked[10];
   Int_t           CCProtonPi0_proton_odMatch[10];
   Int_t           CCProtonPi0_proton_trk_pat_history[10];
   Int_t           CCProtonPi0_trajProtonProngPDG[10];
   Int_t           CCProtonPi0_trajProtonProngPrimary[10];
   Double_t        CCProtonPi0_endProtonTrajMomentum[10];
   Double_t        CCProtonPi0_endProtonTrajXPosition[10];
   Double_t        CCProtonPi0_endProtonTrajYPosition[10];
   Double_t        CCProtonPi0_endProtonTrajZPosition[10];
   Double_t        CCProtonPi0_pionScore[10];
   Double_t        CCProtonPi0_protonScore[10];
   Double_t        CCProtonPi0_proton_E[10];
   Double_t        CCProtonPi0_proton_chi2_ndf[10];
   Double_t        CCProtonPi0_proton_ekin[10];
   Double_t        CCProtonPi0_proton_endPointX[10];
   Double_t        CCProtonPi0_proton_endPointY[10];
   Double_t        CCProtonPi0_proton_endPointZ[10];
   Double_t        CCProtonPi0_proton_p[10];
   Double_t        CCProtonPi0_proton_p_calCorrection[10];
   Double_t        CCProtonPi0_proton_p_dEdXTool[10];
   Double_t        CCProtonPi0_proton_p_visEnergy[10];
   Double_t        CCProtonPi0_proton_phi[10];
   Double_t        CCProtonPi0_proton_px[10];
   Double_t        CCProtonPi0_proton_py[10];
   Double_t        CCProtonPi0_proton_pz[10];
   Double_t        CCProtonPi0_proton_score[10];
   Double_t        CCProtonPi0_proton_score1[10];
   Double_t        CCProtonPi0_proton_score2[10];
   Double_t        CCProtonPi0_proton_startPointX[10];
   Double_t        CCProtonPi0_proton_startPointY[10];
   Double_t        CCProtonPi0_proton_startPointZ[10];
   Double_t        CCProtonPi0_proton_theta[10];
   Double_t        CCProtonPi0_proton_thetaX[10];
   Double_t        CCProtonPi0_proton_thetaY[10];
   Double_t        CCProtonPi0_trajProtonPhi[10];
   Double_t        CCProtonPi0_trajProtonProngEnergy[10];
   Double_t        CCProtonPi0_trajProtonProngMomentum[10];
   Double_t        CCProtonPi0_trajProtonProngPSelf[10];
   Double_t        CCProtonPi0_trajProtonProngPx[10];
   Double_t        CCProtonPi0_trajProtonProngPy[10];
   Double_t        CCProtonPi0_trajProtonProngPz[10];
   Double_t        CCProtonPi0_trajProtonTheta[10];
   Int_t           ev_run;
   Int_t           ev_subrun;
   Int_t           ev_detector;
   Int_t           ev_triggerType;
   Int_t           ev_gate;
   Int_t           ev_global_gate;
   Int_t           ev_gps_time_sec;
   Int_t           ev_gps_time_usec;
   Int_t           mc_run;
   Int_t           mc_subrun;
   Int_t           mc_nInteractions;
   Int_t           mc_MIState;
   Double_t        mc_pot;
   Int_t           mc_beamConfig;
   Int_t           mc_processType;
   Int_t           mc_nthEvtInSpill;
   Int_t           mc_nthEvtInFile;
   Int_t           mc_intType;
   Int_t           mc_current;
   Int_t           mc_charm;
   Double_t        mc_weight;
   Double_t        mc_XSec;
   Double_t        mc_diffXSec;
   Int_t           mc_incoming;
   Double_t        mc_fluxDriverProb;
   Int_t           mc_targetNucleus;
   Int_t           mc_targetZ;
   Int_t           mc_targetA;
   Int_t           mc_targetNucleon;
   Int_t           mc_struckQuark;
   Int_t           mc_seaQuark;
   Int_t           mc_resID;
   Int_t           mc_primaryLepton;
   Double_t        mc_incomingE;
   Double_t        mc_Bjorkenx;
   Double_t        mc_Bjorkeny;
   Double_t        mc_Q2;
   Double_t        mc_nuT;
   Double_t        mc_w;
   Double_t        mc_vtx[4];
   Double_t        mc_incomingPartVec[4];
   Double_t        mc_initNucVec[4];
   Double_t        mc_primFSLepton[4];
   Int_t           mc_nFSPart;
   Double_t        mc_FSPartPx[24];   //[mc_nFSPart]
   Double_t        mc_FSPartPy[24];   //[mc_nFSPart]
   Double_t        mc_FSPartPz[24];   //[mc_nFSPart]
   Double_t        mc_FSPartE[24];   //[mc_nFSPart]
   Int_t           mc_FSPartPDG[24];   //[mc_nFSPart]
   Int_t           mc_er_nPart;
   Int_t           mc_er_ID[51];   //[mc_er_nPart]
   Int_t           mc_er_status[51];   //[mc_er_nPart]
   Double_t        mc_er_posInNucX[51];   //[mc_er_nPart]
   Double_t        mc_er_posInNucY[51];   //[mc_er_nPart]
   Double_t        mc_er_posInNucZ[51];   //[mc_er_nPart]
   Double_t        mc_er_Px[51];   //[mc_er_nPart]
   Double_t        mc_er_Py[51];   //[mc_er_nPart]
   Double_t        mc_er_Pz[51];   //[mc_er_nPart]
   Double_t        mc_er_E[51];   //[mc_er_nPart]
   Int_t           mc_er_FD[51];   //[mc_er_nPart]
   Int_t           mc_er_LD[51];   //[mc_er_nPart]
   Int_t           mc_er_mother[51];   //[mc_er_nPart]
   Int_t           mc_fr_nNuAncestorIDs;
   Int_t           mc_fr_nuAncestorIDs[7];   //[mc_fr_nNuAncestorIDs]
   Int_t           mc_fr_nuParentID;
   Int_t           mc_fr_decMode;
   Double_t        mc_fr_primProtonVtx[3];
   Double_t        mc_fr_primProtonP[4];
   Double_t        mc_fr_nuParentDecVtx[3];
   Double_t        mc_fr_nuParentProdVtx[3];
   Double_t        mc_fr_nuParentProdP[4];
   Double_t        mc_cvweight_total;
   Double_t        wgt;
   Double_t        mc_cvweight_totalFlux;
   Double_t        mc_cvweight_totalXsec;
   Double_t        mc_cvweight_NA49;
   Int_t           mc_wgt_GENIE_sz;
   Double_t        mc_wgt_GENIE[100];   //[mc_wgt_GENIE_sz]
   Int_t           mc_wgt_Flux_Tertiary_sz;
   Double_t        mc_wgt_Flux_Tertiary[100];   //[mc_wgt_Flux_Tertiary_sz]
   Int_t           mc_wgt_Flux_BeamFocus_sz;
   Double_t        mc_wgt_Flux_BeamFocus[100];   //[mc_wgt_Flux_BeamFocus_sz]
   Int_t           mc_wgt_Flux_NA49_sz;
   Double_t        mc_wgt_Flux_NA49[100];   //[mc_wgt_Flux_NA49_sz]
   Int_t           n_prongs;
   Int_t           prong_nParticles[5];   //[n_prongs]
   Double_t        prong_part_score[5];   //[n_prongs]
   Double_t        prong_part_mass[5];   //[n_prongs]
   Int_t           prong_part_charge[5];   //[n_prongs]
   Int_t           prong_part_pid[5];   //[n_prongs]
   vector<vector<double> > *prong_part_E;
   vector<vector<double> > *prong_part_pos;

   // List of branches
   TBranch        *b_eventID;   //!
   TBranch        *b_physEvtNum;   //!
   TBranch        *b_n_hyps;   //!
   TBranch        *b_processType;   //!
   TBranch        *b_primaryPart;   //!
   TBranch        *b_n_slices;   //!
   TBranch        *b_slice_numbers;   //!
   TBranch        *b_shared_slice;   //!
   TBranch        *b_vtx;   //!
   TBranch        *b_vtxErr;   //!
   TBranch        *b_E;   //!
   TBranch        *b_found_truth;   //!
   TBranch        *b_phys_front_activity;   //!
   TBranch        *b_phys_energy_in_road_upstream_is_rockmuon_consistent;   //!
   TBranch        *b_rock_muons_removed;   //!
   TBranch        *b_minos_track_match;   //!
   TBranch        *b_minos_stub_match;   //!
   TBranch        *b_unknown_helicity;   //!
   TBranch        *b_minos_track_inside_partial_plane;   //!
   TBranch        *b_prim_vtx_has_misassigned_track_direction;   //!
   TBranch        *b_prim_vtx_has_broken_track;   //!
   TBranch        *b_gamma1_isGoodDirection;   //!
   TBranch        *b_gamma1_isGoodPosition;   //!
   TBranch        *b_gamma1_isGoodBlob;   //!
   TBranch        *b_gamma2_isGoodDirection;   //!
   TBranch        *b_gamma2_isGoodPosition;   //!
   TBranch        *b_gamma2_isGoodBlob;   //!
   TBranch        *b_Cut_ConeBlobs;   //!
   TBranch        *b_Cut_EndPoint_Michel_Exist;   //!
   TBranch        *b_Cut_Event_Has_BadObject;   //!
   TBranch        *b_Cut_Muon_Charge;   //!
   TBranch        *b_Cut_Muon_None;   //!
   TBranch        *b_Cut_Muon_Not_Plausible;   //!
   TBranch        *b_Cut_Muon_Score_Low;   //!
   TBranch        *b_Cut_Particle_None;   //!
   TBranch        *b_Cut_PreFilter_Pi0;   //!
   TBranch        *b_Cut_Proton_None;   //!
   TBranch        *b_Cut_Vertex_Michel_Exist;   //!
   TBranch        *b_Cut_Vertex_None;   //!
   TBranch        *b_Cut_Vertex_Not_Fiducial;   //!
   TBranch        *b_Cut_Vertex_Not_Reconstructable;   //!
   TBranch        *b_Cut_Vertex_Null;   //!
   TBranch        *b_Cut_VtxBlob;   //!
   TBranch        *b_Cut_nProngs;   //!
   TBranch        *b_Cut_secEndPoint_Michel_Exist;   //!
   TBranch        *b_anglescan_ncand;   //!
   TBranch        *b_anglescan_ncandx;   //!
   TBranch        *b_broken_track_most_us_plane;   //!
   TBranch        *b_g1dedx_doublet;   //!
   TBranch        *b_g1dedx_empty_plane;   //!
   TBranch        *b_g1dedx_nplane;   //!
   TBranch        *b_g2dedx_doublet;   //!
   TBranch        *b_g2dedx_empty_plane;   //!
   TBranch        *b_g2dedx_nplane;   //!
   TBranch        *b_gamma1_blob_nclusters;   //!
   TBranch        *b_gamma1_blob_ndigits;   //!
   TBranch        *b_gamma2_blob_nclusters;   //!
   TBranch        *b_gamma2_blob_ndigits;   //!
   TBranch        *b_nProngs;   //!
   TBranch        *b_n_anchored_long_trk_prongs;   //!
   TBranch        *b_n_anchored_short_trk_prongs;   //!
   TBranch        *b_n_iso_trk_prongs;   //!
   TBranch        *b_n_long_tracks;   //!
   TBranch        *b_n_short_tracks;   //!
   TBranch        *b_n_vtx_michel_views;   //!
   TBranch        *b_nblob_anglescan;   //!
   TBranch        *b_nblob_hough;   //!
   TBranch        *b_od_energeticTower;   //!
   TBranch        *b_phys_energy_in_road_downstream_nplanes;   //!
   TBranch        *b_phys_energy_in_road_upstream_nplanes;   //!
   TBranch        *b_phys_n_dead_discr_pair;   //!
   TBranch        *b_phys_n_dead_discr_pair_in_prim_track_region;   //!
   TBranch        *b_phys_n_dead_discr_pair_two_mod_downstream_prim_track;   //!
   TBranch        *b_phys_n_dead_discr_pair_two_mod_upstream_prim_vtx;   //!
   TBranch        *b_phys_n_dead_discr_pair_upstream_prim_track_proj;   //!
   TBranch        *b_phys_vertex_is_fiducial;   //!
   TBranch        *b_preFilter_Result;   //!
   TBranch        *b_vtx_primary_index;   //!
   TBranch        *b_vtx_primary_multiplicity;   //!
   TBranch        *b_vtx_secondary_count;   //!
   TBranch        *b_vtx_total_count;   //!
   TBranch        *b_Filament_Vertex_energy;   //!
   TBranch        *b_RE_energy_ECAL;   //!
   TBranch        *b_RE_energy_HCAL;   //!
   TBranch        *b_RE_energy_Tracker;   //!
   TBranch        *b_Sphere_Vertex_energy;   //!
   TBranch        *b_Vertex_blob_energy;   //!
   TBranch        *b_dispersedExtraE;   //!
   TBranch        *b_energy_from_mc;   //!
   TBranch        *b_energy_from_mc_fraction;   //!
   TBranch        *b_energy_from_mc_fraction_of_highest;   //!
   TBranch        *b_evis_ecal;   //!
   TBranch        *b_evis_ecal_u;   //!
   TBranch        *b_evis_ecal_v;   //!
   TBranch        *b_evis_ecal_x;   //!
   TBranch        *b_evis_hcal;   //!
   TBranch        *b_evis_hcal_u;   //!
   TBranch        *b_evis_hcal_v;   //!
   TBranch        *b_evis_hcal_x;   //!
   TBranch        *b_evis_nearvtx_total;   //!
   TBranch        *b_evis_nearvtx_u;   //!
   TBranch        *b_evis_nearvtx_v;   //!
   TBranch        *b_evis_nearvtx_x;   //!
   TBranch        *b_evis_ntgt;   //!
   TBranch        *b_evis_ntgt_u;   //!
   TBranch        *b_evis_ntgt_v;   //!
   TBranch        *b_evis_ntgt_x;   //!
   TBranch        *b_evis_other;   //!
   TBranch        *b_evis_total;   //!
   TBranch        *b_evis_total_u;   //!
   TBranch        *b_evis_total_v;   //!
   TBranch        *b_evis_total_x;   //!
   TBranch        *b_evis_trkr;   //!
   TBranch        *b_evis_trkr_u;   //!
   TBranch        *b_evis_trkr_v;   //!
   TBranch        *b_evis_trkr_x;   //!
   TBranch        *b_g1blob_minsep;   //!
   TBranch        *b_g1dedx;   //!
   TBranch        *b_g1dedx1;   //!
   TBranch        *b_g1dedx_total;   //!
   TBranch        *b_g1dedx_total1;   //!
   TBranch        *b_g2blob_minsep;   //!
   TBranch        *b_g2dedx;   //!
   TBranch        *b_g2dedx1;   //!
   TBranch        *b_g2dedx_total;   //!
   TBranch        *b_g2dedx_total1;   //!
   TBranch        *b_gamma1_E;   //!
   TBranch        *b_gamma1_dEdx;   //!
   TBranch        *b_gamma1_dist_exit;   //!
   TBranch        *b_gamma1_dist_vtx;   //!
   TBranch        *b_gamma1_evis_ecal;   //!
   TBranch        *b_gamma1_evis_hcal;   //!
   TBranch        *b_gamma1_evis_scal;   //!
   TBranch        *b_gamma1_evis_trkr;   //!
   TBranch        *b_gamma1_phi;   //!
   TBranch        *b_gamma1_px;   //!
   TBranch        *b_gamma1_py;   //!
   TBranch        *b_gamma1_pz;   //!
   TBranch        *b_gamma1_theta;   //!
   TBranch        *b_gamma1_time;   //!
   TBranch        *b_gamma2_E;   //!
   TBranch        *b_gamma2_dEdx;   //!
   TBranch        *b_gamma2_dist_exit;   //!
   TBranch        *b_gamma2_dist_vtx;   //!
   TBranch        *b_gamma2_evis_ecal;   //!
   TBranch        *b_gamma2_evis_hcal;   //!
   TBranch        *b_gamma2_evis_scal;   //!
   TBranch        *b_gamma2_evis_trkr;   //!
   TBranch        *b_gamma2_phi;   //!
   TBranch        *b_gamma2_px;   //!
   TBranch        *b_gamma2_py;   //!
   TBranch        *b_gamma2_pz;   //!
   TBranch        *b_gamma2_theta;   //!
   TBranch        *b_gamma2_time;   //!
   TBranch        *b_hadronVisibleE;   //!
   TBranch        *b_muonVisibleE;   //!
   TBranch        *b_muon_phi;   //!
   TBranch        *b_muon_theta;   //!
   TBranch        *b_muon_thetaX;   //!
   TBranch        *b_muon_thetaY;   //!
   TBranch        *b_od_downstreamFrame;   //!
   TBranch        *b_od_downstreamFrame_z;   //!
   TBranch        *b_od_highStory;   //!
   TBranch        *b_od_highStory_t;   //!
   TBranch        *b_od_lowStory;   //!
   TBranch        *b_od_lowStory_t;   //!
   TBranch        *b_od_maxEnergy;   //!
   TBranch        *b_od_upstreamFrame;   //!
   TBranch        *b_od_upstreamFrame_z;   //!
   TBranch        *b_phys_energy_dispersed;   //!
   TBranch        *b_phys_energy_in_road_downstream;   //!
   TBranch        *b_phys_energy_in_road_upstream;   //!
   TBranch        *b_phys_energy_unattached;   //!
   TBranch        *b_pi0_E;   //!
   TBranch        *b_pi0_cos_openingAngle;   //!
   TBranch        *b_pi0_invMass;   //!
   TBranch        *b_pi0_openingAngle;   //!
   TBranch        *b_pi0_phi;   //!
   TBranch        *b_pi0_px;   //!
   TBranch        *b_pi0_py;   //!
   TBranch        *b_pi0_pz;   //!
   TBranch        *b_pi0_theta;   //!
   TBranch        *b_pi0_thetaX;   //!
   TBranch        *b_pi0_thetaY;   //!
   TBranch        *b_preFilter_rejectedEnergy;   //!
   TBranch        *b_prim_vtx_smallest_opening_angle;   //!
   TBranch        *b_time;   //!
   TBranch        *b_totalIDVisibleE;   //!
   TBranch        *b_totalODVisibleE;   //!
   TBranch        *b_totalVisibleE;   //!
   TBranch        *b_unattachedExtraE;   //!
   TBranch        *b_vtxBlobExtraE;   //!
   TBranch        *b_vtx_michel_distance;   //!
   TBranch        *b_anglescan_blob_nc_sz;   //!
   TBranch        *b_anglescan_blob_nc;   //!
   TBranch        *b_anglescan_blob_ncu_sz;   //!
   TBranch        *b_anglescan_blob_ncu;   //!
   TBranch        *b_anglescan_blob_ncv_sz;   //!
   TBranch        *b_anglescan_blob_ncv;   //!
   TBranch        *b_anglescan_blob_ncx_sz;   //!
   TBranch        *b_anglescan_blob_ncx;   //!
   TBranch        *b_anglescan_blob_nd_sz;   //!
   TBranch        *b_anglescan_blob_nd;   //!
   TBranch        *b_anglescan_blob_ndu_sz;   //!
   TBranch        *b_anglescan_blob_ndu;   //!
   TBranch        *b_anglescan_blob_ndv_sz;   //!
   TBranch        *b_anglescan_blob_ndv;   //!
   TBranch        *b_anglescan_blob_ndx_sz;   //!
   TBranch        *b_anglescan_blob_ndx;   //!
   TBranch        *b_anglescan_cand_nc_sz;   //!
   TBranch        *b_anglescan_cand_nc;   //!
   TBranch        *b_anglescan_cand_ncu_sz;   //!
   TBranch        *b_anglescan_cand_ncu;   //!
   TBranch        *b_anglescan_cand_ncv_sz;   //!
   TBranch        *b_anglescan_cand_ncv;   //!
   TBranch        *b_anglescan_cand_ncx_sz;   //!
   TBranch        *b_anglescan_cand_ncx;   //!
   TBranch        *b_anglescan_cand_nd_sz;   //!
   TBranch        *b_anglescan_cand_nd;   //!
   TBranch        *b_anglescan_cand_ndu_sz;   //!
   TBranch        *b_anglescan_cand_ndu;   //!
   TBranch        *b_anglescan_cand_ndv_sz;   //!
   TBranch        *b_anglescan_cand_ndv;   //!
   TBranch        *b_anglescan_cand_ndx_sz;   //!
   TBranch        *b_anglescan_cand_ndx;   //!
   TBranch        *b_anglescan_candx_nc_sz;   //!
   TBranch        *b_anglescan_candx_nc;   //!
   TBranch        *b_anglescan_candx_nd_sz;   //!
   TBranch        *b_anglescan_candx_nd;   //!
   TBranch        *b_final_blob_nc_sz;   //!
   TBranch        *b_final_blob_nc;   //!
   TBranch        *b_final_blob_ncu_sz;   //!
   TBranch        *b_final_blob_ncu;   //!
   TBranch        *b_final_blob_ncv_sz;   //!
   TBranch        *b_final_blob_ncv;   //!
   TBranch        *b_final_blob_ncx_sz;   //!
   TBranch        *b_final_blob_ncx;   //!
   TBranch        *b_final_blob_nd_sz;   //!
   TBranch        *b_final_blob_nd;   //!
   TBranch        *b_final_blob_ndu_sz;   //!
   TBranch        *b_final_blob_ndu;   //!
   TBranch        *b_final_blob_ndv_sz;   //!
   TBranch        *b_final_blob_ndv;   //!
   TBranch        *b_final_blob_ndx_sz;   //!
   TBranch        *b_final_blob_ndx;   //!
   TBranch        *b_g1dedx_cluster_occupancy_sz;   //!
   TBranch        *b_g1dedx_cluster_occupancy;   //!
   TBranch        *b_g2dedx_cluster_occupancy_sz;   //!
   TBranch        *b_g2dedx_cluster_occupancy;   //!
   TBranch        *b_hough_blob_nc_sz;   //!
   TBranch        *b_hough_blob_nc;   //!
   TBranch        *b_hough_blob_ncu_sz;   //!
   TBranch        *b_hough_blob_ncu;   //!
   TBranch        *b_hough_blob_ncv_sz;   //!
   TBranch        *b_hough_blob_ncv;   //!
   TBranch        *b_hough_blob_ncx_sz;   //!
   TBranch        *b_hough_blob_ncx;   //!
   TBranch        *b_hough_blob_nd_sz;   //!
   TBranch        *b_hough_blob_nd;   //!
   TBranch        *b_hough_blob_ndu_sz;   //!
   TBranch        *b_hough_blob_ndu;   //!
   TBranch        *b_hough_blob_ndv_sz;   //!
   TBranch        *b_hough_blob_ndv;   //!
   TBranch        *b_hough_blob_ndx_sz;   //!
   TBranch        *b_hough_blob_ndx;   //!
   TBranch        *b_Vertex_energy_radii_sz;   //!
   TBranch        *b_Vertex_energy_radii;   //!
   TBranch        *b_blob_cluster_energy1_sz;   //!
   TBranch        *b_blob_cluster_energy1;   //!
   TBranch        *b_blob_cluster_energy2_sz;   //!
   TBranch        *b_blob_cluster_energy2;   //!
   TBranch        *b_g1dedx_cluster_energy_sz;   //!
   TBranch        *b_g1dedx_cluster_energy;   //!
   TBranch        *b_g1dedx_rev_cluster_energy_sz;   //!
   TBranch        *b_g1dedx_rev_cluster_energy;   //!
   TBranch        *b_g2dedx_cluster_energy_sz;   //!
   TBranch        *b_g2dedx_cluster_energy;   //!
   TBranch        *b_g2dedx_rev_cluster_energy_sz;   //!
   TBranch        *b_g2dedx_rev_cluster_energy;   //!
   TBranch        *b_gamma1_direction;   //!
   TBranch        *b_gamma1_vertex;   //!
   TBranch        *b_gamma2_direction;   //!
   TBranch        *b_gamma2_vertex;   //!
   TBranch        *b_good_mgg_vector_sz;   //!
   TBranch        *b_good_mgg_vector;   //!
   TBranch        *b_mgg_vector_sz;   //!
   TBranch        *b_mgg_vector;   //!
   TBranch        *b_od_distanceBlobTower_sz;   //!
   TBranch        *b_od_distanceBlobTower;   //!
   TBranch        *b_od_idBlobTime_sz;   //!
   TBranch        *b_od_idBlobTime;   //!
   TBranch        *b_od_towerEnergy_sz;   //!
   TBranch        *b_od_towerEnergy;   //!
   TBranch        *b_od_towerNClusters_sz;   //!
   TBranch        *b_od_towerNClusters;   //!
   TBranch        *b_od_towerTime_sz;   //!
   TBranch        *b_od_towerTime;   //!
   TBranch        *b_od_towerTimeBlobMuon_sz;   //!
   TBranch        *b_od_towerTimeBlobMuon;   //!
   TBranch        *b_od_towerTimeBlobOD_sz;   //!
   TBranch        *b_od_towerTimeBlobOD;   //!
   TBranch        *b_truth_has_physics_event;   //!
   TBranch        *b_truth_isSignal;   //!
   TBranch        *b_truth_isSignalGold;   //!
   TBranch        *b_truth_isSignalSilver1;   //!
   TBranch        *b_truth_isSignal_1Pi0;   //!
   TBranch        *b_truth_isSignal_2Gamma;   //!
   TBranch        *b_truth_isFidVol;   //!
   TBranch        *b_truth_N_FSParticles;   //!
   TBranch        *b_truth_N_gamma;   //!
   TBranch        *b_truth_N_pi0;   //!
   TBranch        *b_truth_N_proton;   //!
   TBranch        *b_truth_N_proton_Long;   //!
   TBranch        *b_truth_N_proton_Short;   //!
   TBranch        *b_truth_muon_charge;   //!
   TBranch        *b_truth_pi0_GrandMother;   //!
   TBranch        *b_truth_pi0_GrandMotherStatus;   //!
   TBranch        *b_truth_pi0_Mother;   //!
   TBranch        *b_truth_pi0_MotherStatus;   //!
   TBranch        *b_truth_pi0_status;   //!
   TBranch        *b_truth_target_material;   //!
   TBranch        *b_truth_vertex_module;   //!
   TBranch        *b_truth_vertex_plane;   //!
   TBranch        *b_truth_muon_E;   //!
   TBranch        *b_truth_muon_px;   //!
   TBranch        *b_truth_muon_py;   //!
   TBranch        *b_truth_muon_pz;   //!
   TBranch        *b_truth_muon_theta_wrtbeam;   //!
   TBranch        *b_truth_pi0_E;   //!
   TBranch        *b_truth_pi0_px;   //!
   TBranch        *b_truth_pi0_py;   //!
   TBranch        *b_truth_pi0_pz;   //!
   TBranch        *b_truth_pi0_theta_wrtbeam;   //!
   TBranch        *b_truth_gamma_E;   //!
   TBranch        *b_truth_gamma_px;   //!
   TBranch        *b_truth_gamma_py;   //!
   TBranch        *b_truth_gamma_pz;   //!
   TBranch        *b_truth_gamma_theta_wrtbeam;   //!
   TBranch        *b_genie_wgt_n_shifts;   //!
   TBranch        *b_truth_genie_wgt_AGKYxF1pi;   //!
   TBranch        *b_truth_genie_wgt_AhtBY;   //!
   TBranch        *b_truth_genie_wgt_BhtBY;   //!
   TBranch        *b_truth_genie_wgt_CCQEPauliSupViaKF;   //!
   TBranch        *b_truth_genie_wgt_CV1uBY;   //!
   TBranch        *b_truth_genie_wgt_CV2uBY;   //!
   TBranch        *b_truth_genie_wgt_EtaNCEL;   //!
   TBranch        *b_truth_genie_wgt_FrAbs_N;   //!
   TBranch        *b_truth_genie_wgt_FrAbs_pi;   //!
   TBranch        *b_truth_genie_wgt_FrCEx_N;   //!
   TBranch        *b_truth_genie_wgt_FrCEx_pi;   //!
   TBranch        *b_truth_genie_wgt_FrElas_N;   //!
   TBranch        *b_truth_genie_wgt_FrElas_pi;   //!
   TBranch        *b_truth_genie_wgt_FrInel_N;   //!
   TBranch        *b_truth_genie_wgt_FrInel_pi;   //!
   TBranch        *b_truth_genie_wgt_FrPiProd_N;   //!
   TBranch        *b_truth_genie_wgt_FrPiProd_pi;   //!
   TBranch        *b_truth_genie_wgt_MFP_N;   //!
   TBranch        *b_truth_genie_wgt_MFP_pi;   //!
   TBranch        *b_truth_genie_wgt_MaCCQE;   //!
   TBranch        *b_truth_genie_wgt_MaCCQEshape;   //!
   TBranch        *b_truth_genie_wgt_MaNCEL;   //!
   TBranch        *b_truth_genie_wgt_MaRES;   //!
   TBranch        *b_truth_genie_wgt_MvRES;   //!
   TBranch        *b_truth_genie_wgt_NormCCQE;   //!
   TBranch        *b_truth_genie_wgt_NormCCRES;   //!
   TBranch        *b_truth_genie_wgt_NormDISCC;   //!
   TBranch        *b_truth_genie_wgt_NormNCRES;   //!
   TBranch        *b_truth_genie_wgt_RDecBR1gamma;   //!
   TBranch        *b_truth_genie_wgt_Rvn1pi;   //!
   TBranch        *b_truth_genie_wgt_Rvn2pi;   //!
   TBranch        *b_truth_genie_wgt_Rvp1pi;   //!
   TBranch        *b_truth_genie_wgt_Rvp2pi;   //!
   TBranch        *b_truth_genie_wgt_Theta_Delta2Npi;   //!
   TBranch        *b_truth_genie_wgt_VecFFCCQEshape;   //!
   TBranch        *b_truth_genie_wgt_shifts;   //!
   TBranch        *b_truth_proton_E;   //!
   TBranch        *b_truth_proton_px;   //!
   TBranch        *b_truth_proton_py;   //!
   TBranch        *b_truth_proton_pz;   //!
   TBranch        *b_truth_proton_theta_wrtbeam;   //!
   TBranch        *b_CCProtonPi0_nuFlavor;   //!
   TBranch        *b_CCProtonPi0_nuHelicity;   //!
   TBranch        *b_CCProtonPi0_intCurrent;   //!
   TBranch        *b_CCProtonPi0_intType;   //!
   TBranch        *b_CCProtonPi0_E;   //!
   TBranch        *b_CCProtonPi0_Q2;   //!
   TBranch        *b_CCProtonPi0_x;   //!
   TBranch        *b_CCProtonPi0_y;   //!
   TBranch        *b_CCProtonPi0_W;   //!
   TBranch        *b_CCProtonPi0_score;   //!
   TBranch        *b_CCProtonPi0_leptonE;   //!
   TBranch        *b_CCProtonPi0_vtx;   //!
   TBranch        *b_CCProtonPi0_minos_trk_is_contained;   //!
   TBranch        *b_CCProtonPi0_minos_trk_is_ok;   //!
   TBranch        *b_CCProtonPi0_minos_used_range;   //!
   TBranch        *b_CCProtonPi0_minos_used_curvature;   //!
   TBranch        *b_CCProtonPi0_isMuonInsideOD;   //!
   TBranch        *b_CCProtonPi0_minos_trk_end_plane;   //!
   TBranch        *b_CCProtonPi0_minos_trk_quality;   //!
   TBranch        *b_CCProtonPi0_muon_N_minosTracks;   //!
   TBranch        *b_CCProtonPi0_muon_charge;   //!
   TBranch        *b_CCProtonPi0_muon_hasMinosMatchStub;   //!
   TBranch        *b_CCProtonPi0_muon_hasMinosMatchTrack;   //!
   TBranch        *b_CCProtonPi0_muon_minervaTrack_types;   //!
   TBranch        *b_CCProtonPi0_muon_minosTrackQuality;   //!
   TBranch        *b_CCProtonPi0_muon_roadUpstreamPlanes;   //!
   TBranch        *b_CCProtonPi0_ntrajMuonProng;   //!
   TBranch        *b_CCProtonPi0_r_minos_trk_vtx_plane;   //!
   TBranch        *b_CCProtonPi0_t_minos_trk_numFSMuons;   //!
   TBranch        *b_CCProtonPi0_t_minos_trk_primFSLeptonPDG;   //!
   TBranch        *b_CCProtonPi0_trajMuonProngPDG;   //!
   TBranch        *b_CCProtonPi0_trajMuonProngPrimary;   //!
   TBranch        *b_CCProtonPi0_vtx_module;   //!
   TBranch        *b_CCProtonPi0_vtx_plane;   //!
   TBranch        *b_CCProtonPi0_QSq;   //!
   TBranch        *b_CCProtonPi0_WSq;   //!
   TBranch        *b_CCProtonPi0_endMuonTrajMomentum;   //!
   TBranch        *b_CCProtonPi0_endMuonTrajXPosition;   //!
   TBranch        *b_CCProtonPi0_endMuonTrajYPosition;   //!
   TBranch        *b_CCProtonPi0_endMuonTrajZPosition;   //!
   TBranch        *b_CCProtonPi0_minos_trk_bave;   //!
   TBranch        *b_CCProtonPi0_minos_trk_chi2;   //!
   TBranch        *b_CCProtonPi0_minos_trk_end_u;   //!
   TBranch        *b_CCProtonPi0_minos_trk_end_v;   //!
   TBranch        *b_CCProtonPi0_minos_trk_end_x;   //!
   TBranch        *b_CCProtonPi0_minos_trk_end_y;   //!
   TBranch        *b_CCProtonPi0_minos_trk_end_z;   //!
   TBranch        *b_CCProtonPi0_minos_trk_eqp;   //!
   TBranch        *b_CCProtonPi0_minos_trk_eqp_qp;   //!
   TBranch        *b_CCProtonPi0_minos_trk_fit_pass;   //!
   TBranch        *b_CCProtonPi0_minos_trk_ndf;   //!
   TBranch        *b_CCProtonPi0_minos_trk_p;   //!
   TBranch        *b_CCProtonPi0_minos_trk_p_curvature;   //!
   TBranch        *b_CCProtonPi0_minos_trk_p_range;   //!
   TBranch        *b_CCProtonPi0_minos_trk_qp;   //!
   TBranch        *b_CCProtonPi0_minos_trk_vtx_x;   //!
   TBranch        *b_CCProtonPi0_minos_trk_vtx_y;   //!
   TBranch        *b_CCProtonPi0_minos_trk_vtx_z;   //!
   TBranch        *b_CCProtonPi0_muon_E;   //!
   TBranch        *b_CCProtonPi0_muon_E_shift;   //!
   TBranch        *b_CCProtonPi0_muon_muScore;   //!
   TBranch        *b_CCProtonPi0_muon_p;   //!
   TBranch        *b_CCProtonPi0_muon_px;   //!
   TBranch        *b_CCProtonPi0_muon_py;   //!
   TBranch        *b_CCProtonPi0_muon_pz;   //!
   TBranch        *b_CCProtonPi0_muon_qp;   //!
   TBranch        *b_CCProtonPi0_muon_qpqpe;   //!
   TBranch        *b_CCProtonPi0_muon_roadUpstreamEnergy;   //!
   TBranch        *b_CCProtonPi0_muon_theta;   //!
   TBranch        *b_CCProtonPi0_muon_theta_biasDown;   //!
   TBranch        *b_CCProtonPi0_muon_theta_biasUp;   //!
   TBranch        *b_CCProtonPi0_neutrino_E;   //!
   TBranch        *b_CCProtonPi0_neutrino_E_Cal;   //!
   TBranch        *b_CCProtonPi0_r_minos_trk_bdL;   //!
   TBranch        *b_CCProtonPi0_r_minos_trk_end_dcosx;   //!
   TBranch        *b_CCProtonPi0_r_minos_trk_end_dcosy;   //!
   TBranch        *b_CCProtonPi0_r_minos_trk_end_dcosz;   //!
   TBranch        *b_CCProtonPi0_r_minos_trk_vtx_dcosx;   //!
   TBranch        *b_CCProtonPi0_r_minos_trk_vtx_dcosy;   //!
   TBranch        *b_CCProtonPi0_r_minos_trk_vtx_dcosz;   //!
   TBranch        *b_CCProtonPi0_t_minos_trk_primFSLepMinosInitProjPx;   //!
   TBranch        *b_CCProtonPi0_t_minos_trk_primFSLepMinosInitProjPy;   //!
   TBranch        *b_CCProtonPi0_t_minos_trk_primFSLepMinosInitProjPz;   //!
   TBranch        *b_CCProtonPi0_t_minos_trk_primFSLepMinosInitProjX;   //!
   TBranch        *b_CCProtonPi0_t_minos_trk_primFSLepMinosInitProjY;   //!
   TBranch        *b_CCProtonPi0_t_minos_trk_primFSLepMinosInitProjZ;   //!
   TBranch        *b_CCProtonPi0_t_minos_trk_primFSLepMnvFinalPx;   //!
   TBranch        *b_CCProtonPi0_t_minos_trk_primFSLepMnvFinalPy;   //!
   TBranch        *b_CCProtonPi0_t_minos_trk_primFSLepMnvFinalPz;   //!
   TBranch        *b_CCProtonPi0_t_minos_trk_primFSLepMnvFinalX;   //!
   TBranch        *b_CCProtonPi0_t_minos_trk_primFSLepMnvFinalY;   //!
   TBranch        *b_CCProtonPi0_t_minos_trk_primFSLepMnvFinalZ;   //!
   TBranch        *b_CCProtonPi0_t_minos_trk_primFSLepMnvInitPx;   //!
   TBranch        *b_CCProtonPi0_t_minos_trk_primFSLepMnvInitPy;   //!
   TBranch        *b_CCProtonPi0_t_minos_trk_primFSLepMnvInitPz;   //!
   TBranch        *b_CCProtonPi0_t_minos_trk_primFSLepMnvInitX;   //!
   TBranch        *b_CCProtonPi0_t_minos_trk_primFSLepMnvInitY;   //!
   TBranch        *b_CCProtonPi0_t_minos_trk_primFSLepMnvInitZ;   //!
   TBranch        *b_CCProtonPi0_trajMuonPhi;   //!
   TBranch        *b_CCProtonPi0_trajMuonProngEnergy;   //!
   TBranch        *b_CCProtonPi0_trajMuonProngMomentum;   //!
   TBranch        *b_CCProtonPi0_trajMuonProngPSelf;   //!
   TBranch        *b_CCProtonPi0_trajMuonProngPx;   //!
   TBranch        *b_CCProtonPi0_trajMuonProngPy;   //!
   TBranch        *b_CCProtonPi0_trajMuonProngPz;   //!
   TBranch        *b_CCProtonPi0_trajMuonTheta;   //!
   TBranch        *b_CCProtonPi0_vtx_x;   //!
   TBranch        *b_CCProtonPi0_vtx_y;   //!
   TBranch        *b_CCProtonPi0_vtx_z;   //!
   TBranch        *b_CCProtonPi0_isProtonInsideOD;   //!
   TBranch        *b_CCProtonPi0_ntrajProtonProng;   //!
   TBranch        *b_CCProtonPi0_proton_isRecoGood;   //!
   TBranch        *b_CCProtonPi0_proton_kinked;   //!
   TBranch        *b_CCProtonPi0_proton_odMatch;   //!
   TBranch        *b_CCProtonPi0_proton_trk_pat_history;   //!
   TBranch        *b_CCProtonPi0_trajProtonProngPDG;   //!
   TBranch        *b_CCProtonPi0_trajProtonProngPrimary;   //!
   TBranch        *b_CCProtonPi0_endProtonTrajMomentum;   //!
   TBranch        *b_CCProtonPi0_endProtonTrajXPosition;   //!
   TBranch        *b_CCProtonPi0_endProtonTrajYPosition;   //!
   TBranch        *b_CCProtonPi0_endProtonTrajZPosition;   //!
   TBranch        *b_CCProtonPi0_pionScore;   //!
   TBranch        *b_CCProtonPi0_protonScore;   //!
   TBranch        *b_CCProtonPi0_proton_E;   //!
   TBranch        *b_CCProtonPi0_proton_chi2_ndf;   //!
   TBranch        *b_CCProtonPi0_proton_ekin;   //!
   TBranch        *b_CCProtonPi0_proton_endPointX;   //!
   TBranch        *b_CCProtonPi0_proton_endPointY;   //!
   TBranch        *b_CCProtonPi0_proton_endPointZ;   //!
   TBranch        *b_CCProtonPi0_proton_p;   //!
   TBranch        *b_CCProtonPi0_proton_p_calCorrection;   //!
   TBranch        *b_CCProtonPi0_proton_p_dEdXTool;   //!
   TBranch        *b_CCProtonPi0_proton_p_visEnergy;   //!
   TBranch        *b_CCProtonPi0_proton_phi;   //!
   TBranch        *b_CCProtonPi0_proton_px;   //!
   TBranch        *b_CCProtonPi0_proton_py;   //!
   TBranch        *b_CCProtonPi0_proton_pz;   //!
   TBranch        *b_CCProtonPi0_proton_score;   //!
   TBranch        *b_CCProtonPi0_proton_score1;   //!
   TBranch        *b_CCProtonPi0_proton_score2;   //!
   TBranch        *b_CCProtonPi0_proton_startPointX;   //!
   TBranch        *b_CCProtonPi0_proton_startPointY;   //!
   TBranch        *b_CCProtonPi0_proton_startPointZ;   //!
   TBranch        *b_CCProtonPi0_proton_theta;   //!
   TBranch        *b_CCProtonPi0_proton_thetaX;   //!
   TBranch        *b_CCProtonPi0_proton_thetaY;   //!
   TBranch        *b_CCProtonPi0_trajProtonPhi;   //!
   TBranch        *b_CCProtonPi0_trajProtonProngEnergy;   //!
   TBranch        *b_CCProtonPi0_trajProtonProngMomentum;   //!
   TBranch        *b_CCProtonPi0_trajProtonProngPSelf;   //!
   TBranch        *b_CCProtonPi0_trajProtonProngPx;   //!
   TBranch        *b_CCProtonPi0_trajProtonProngPy;   //!
   TBranch        *b_CCProtonPi0_trajProtonProngPz;   //!
   TBranch        *b_CCProtonPi0_trajProtonTheta;   //!
   TBranch        *b_ev_run;   //!
   TBranch        *b_ev_subrun;   //!
   TBranch        *b_ev_detector;   //!
   TBranch        *b_ev_triggerType;   //!
   TBranch        *b_ev_gate;   //!
   TBranch        *b_ev_global_gate;   //!
   TBranch        *b_ev_gps_time_sec;   //!
   TBranch        *b_ev_gps_time_usec;   //!
   TBranch        *b_mc_run;   //!
   TBranch        *b_mc_subrun;   //!
   TBranch        *b_mc_nInteractions;   //!
   TBranch        *b_mc_MIState;   //!
   TBranch        *b_mc_pot;   //!
   TBranch        *b_mc_beamConfig;   //!
   TBranch        *b_mc_processType;   //!
   TBranch        *b_mc_nthEvtInSpill;   //!
   TBranch        *b_mc_nthEvtInFile;   //!
   TBranch        *b_mc_intType;   //!
   TBranch        *b_mc_current;   //!
   TBranch        *b_mc_charm;   //!
   TBranch        *b_mc_weight;   //!
   TBranch        *b_mc_XSec;   //!
   TBranch        *b_mc_diffXSec;   //!
   TBranch        *b_mc_incoming;   //!
   TBranch        *b_mc_fluxDriverProb;   //!
   TBranch        *b_mc_targetNucleus;   //!
   TBranch        *b_mc_targetZ;   //!
   TBranch        *b_mc_targetA;   //!
   TBranch        *b_mc_targetNucleon;   //!
   TBranch        *b_mc_struckQuark;   //!
   TBranch        *b_mc_seaQuark;   //!
   TBranch        *b_mc_resID;   //!
   TBranch        *b_mc_primaryLepton;   //!
   TBranch        *b_mc_incomingE;   //!
   TBranch        *b_mc_Bjorkenx;   //!
   TBranch        *b_mc_Bjorkeny;   //!
   TBranch        *b_mc_Q2;   //!
   TBranch        *b_mc_nuT;   //!
   TBranch        *b_mc_w;   //!
   TBranch        *b_mc_vtx;   //!
   TBranch        *b_mc_incomingPartVec;   //!
   TBranch        *b_mc_initNucVec;   //!
   TBranch        *b_mc_primFSLepton;   //!
   TBranch        *b_mc_nFSPart;   //!
   TBranch        *b_mc_FSPartPx;   //!
   TBranch        *b_mc_FSPartPy;   //!
   TBranch        *b_mc_FSPartPz;   //!
   TBranch        *b_mc_FSPartE;   //!
   TBranch        *b_mc_FSPartPDG;   //!
   TBranch        *b_mc_er_nPart;   //!
   TBranch        *b_mc_er_ID;   //!
   TBranch        *b_mc_er_status;   //!
   TBranch        *b_mc_er_posInNucX;   //!
   TBranch        *b_mc_er_posInNucY;   //!
   TBranch        *b_mc_er_posInNucZ;   //!
   TBranch        *b_mc_er_Px;   //!
   TBranch        *b_mc_er_Py;   //!
   TBranch        *b_mc_er_Pz;   //!
   TBranch        *b_mc_er_E;   //!
   TBranch        *b_mc_er_FD;   //!
   TBranch        *b_mc_er_LD;   //!
   TBranch        *b_mc_er_mother;   //!
   TBranch        *b_mc_fr_nNuAncestorIDs;   //!
   TBranch        *b_mc_fr_nuAncestorIDs;   //!
   TBranch        *b_mc_fr_nuParentID;   //!
   TBranch        *b_mc_fr_decMode;   //!
   TBranch        *b_mc_fr_primProtonVtx;   //!
   TBranch        *b_mc_fr_primProtonP;   //!
   TBranch        *b_mc_fr_nuParentDecVtx;   //!
   TBranch        *b_mc_fr_nuParentProdVtx;   //!
   TBranch        *b_mc_fr_nuParentProdP;   //!
   TBranch        *b_mc_cvweight_total;   //!
   TBranch        *b_wgt;   //!
   TBranch        *b_mc_cvweight_totalFlux;   //!
   TBranch        *b_mc_cvweight_totalXsec;   //!
   TBranch        *b_mc_cvweight_NA49;   //!
   TBranch        *b_mc_wgt_GENIE_sz;   //!
   TBranch        *b_mc_wgt_GENIE;   //!
   TBranch        *b_mc_wgt_Flux_Tertiary_sz;   //!
   TBranch        *b_mc_wgt_Flux_Tertiary;   //!
   TBranch        *b_mc_wgt_Flux_BeamFocus_sz;   //!
   TBranch        *b_mc_wgt_Flux_BeamFocus;   //!
   TBranch        *b_mc_wgt_Flux_NA49_sz;   //!
   TBranch        *b_mc_wgt_Flux_NA49;   //!
   TBranch        *b_n_prongs;   //!
   TBranch        *b_prong_nParticles;   //!
   TBranch        *b_prong_part_score;   //!
   TBranch        *b_prong_part_mass;   //!
   TBranch        *b_prong_part_charge;   //!
   TBranch        *b_prong_part_pid;   //!
   TBranch        *b_prong_part_E;   //!
   TBranch        *b_prong_part_pos;   //!

   MC_Sample(TTree *tree=0);
   virtual ~MC_Sample();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MC_Sample_cxx
MC_Sample::MC_Sample(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {

#ifdef SINGLE_TREE
      // The following code should be used if you want this class to access
      // a single tree instead of a chain
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("Memory Directory");
      if (!f || !f->IsOpen()) {
         f = new TFile("Memory Directory");
      }
      f->GetObject("CCProtonPi0",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("CCProtonPi0","");
      chain->Add("/minerva/data/users/oaltinok/CCProtonPi0/MC/v1_07/test/nogrid/central_value/minerva/ana/v10r8/00/01/32/00/SIM_minerva_00013200_Subruns_0005_CCProtonPi0_Ana_Tuple_v10r8.root/CCProtonPi0");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

MC_Sample::~MC_Sample()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MC_Sample::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MC_Sample::LoadTree(Long64_t entry)
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

void MC_Sample::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   prong_part_E = 0;
   prong_part_pos = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
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
   fChain->SetBranchAddress("gamma1_isGoodDirection", &gamma1_isGoodDirection, &b_gamma1_isGoodDirection);
   fChain->SetBranchAddress("gamma1_isGoodPosition", &gamma1_isGoodPosition, &b_gamma1_isGoodPosition);
   fChain->SetBranchAddress("gamma1_isGoodBlob", &gamma1_isGoodBlob, &b_gamma1_isGoodBlob);
   fChain->SetBranchAddress("gamma2_isGoodDirection", &gamma2_isGoodDirection, &b_gamma2_isGoodDirection);
   fChain->SetBranchAddress("gamma2_isGoodPosition", &gamma2_isGoodPosition, &b_gamma2_isGoodPosition);
   fChain->SetBranchAddress("gamma2_isGoodBlob", &gamma2_isGoodBlob, &b_gamma2_isGoodBlob);
   fChain->SetBranchAddress("Cut_ConeBlobs", &Cut_ConeBlobs, &b_Cut_ConeBlobs);
   fChain->SetBranchAddress("Cut_EndPoint_Michel_Exist", &Cut_EndPoint_Michel_Exist, &b_Cut_EndPoint_Michel_Exist);
   fChain->SetBranchAddress("Cut_Event_Has_BadObject", &Cut_Event_Has_BadObject, &b_Cut_Event_Has_BadObject);
   fChain->SetBranchAddress("Cut_Muon_Charge", &Cut_Muon_Charge, &b_Cut_Muon_Charge);
   fChain->SetBranchAddress("Cut_Muon_None", &Cut_Muon_None, &b_Cut_Muon_None);
   fChain->SetBranchAddress("Cut_Muon_Not_Plausible", &Cut_Muon_Not_Plausible, &b_Cut_Muon_Not_Plausible);
   fChain->SetBranchAddress("Cut_Muon_Score_Low", &Cut_Muon_Score_Low, &b_Cut_Muon_Score_Low);
   fChain->SetBranchAddress("Cut_Particle_None", &Cut_Particle_None, &b_Cut_Particle_None);
   fChain->SetBranchAddress("Cut_PreFilter_Pi0", &Cut_PreFilter_Pi0, &b_Cut_PreFilter_Pi0);
   fChain->SetBranchAddress("Cut_Proton_None", &Cut_Proton_None, &b_Cut_Proton_None);
   fChain->SetBranchAddress("Cut_Vertex_Michel_Exist", &Cut_Vertex_Michel_Exist, &b_Cut_Vertex_Michel_Exist);
   fChain->SetBranchAddress("Cut_Vertex_None", &Cut_Vertex_None, &b_Cut_Vertex_None);
   fChain->SetBranchAddress("Cut_Vertex_Not_Fiducial", &Cut_Vertex_Not_Fiducial, &b_Cut_Vertex_Not_Fiducial);
   fChain->SetBranchAddress("Cut_Vertex_Not_Reconstructable", &Cut_Vertex_Not_Reconstructable, &b_Cut_Vertex_Not_Reconstructable);
   fChain->SetBranchAddress("Cut_Vertex_Null", &Cut_Vertex_Null, &b_Cut_Vertex_Null);
   fChain->SetBranchAddress("Cut_VtxBlob", &Cut_VtxBlob, &b_Cut_VtxBlob);
   fChain->SetBranchAddress("Cut_nProngs", &Cut_nProngs, &b_Cut_nProngs);
   fChain->SetBranchAddress("Cut_secEndPoint_Michel_Exist", &Cut_secEndPoint_Michel_Exist, &b_Cut_secEndPoint_Michel_Exist);
   fChain->SetBranchAddress("anglescan_ncand", &anglescan_ncand, &b_anglescan_ncand);
   fChain->SetBranchAddress("anglescan_ncandx", &anglescan_ncandx, &b_anglescan_ncandx);
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
   fChain->SetBranchAddress("n_anchored_long_trk_prongs", &n_anchored_long_trk_prongs, &b_n_anchored_long_trk_prongs);
   fChain->SetBranchAddress("n_anchored_short_trk_prongs", &n_anchored_short_trk_prongs, &b_n_anchored_short_trk_prongs);
   fChain->SetBranchAddress("n_iso_trk_prongs", &n_iso_trk_prongs, &b_n_iso_trk_prongs);
   fChain->SetBranchAddress("n_long_tracks", &n_long_tracks, &b_n_long_tracks);
   fChain->SetBranchAddress("n_short_tracks", &n_short_tracks, &b_n_short_tracks);
   fChain->SetBranchAddress("n_vtx_michel_views", &n_vtx_michel_views, &b_n_vtx_michel_views);
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
   fChain->SetBranchAddress("preFilter_Result", &preFilter_Result, &b_preFilter_Result);
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
   fChain->SetBranchAddress("g1blob_minsep", &g1blob_minsep, &b_g1blob_minsep);
   fChain->SetBranchAddress("g1dedx", &g1dedx, &b_g1dedx);
   fChain->SetBranchAddress("g1dedx1", &g1dedx1, &b_g1dedx1);
   fChain->SetBranchAddress("g1dedx_total", &g1dedx_total, &b_g1dedx_total);
   fChain->SetBranchAddress("g1dedx_total1", &g1dedx_total1, &b_g1dedx_total1);
   fChain->SetBranchAddress("g2blob_minsep", &g2blob_minsep, &b_g2blob_minsep);
   fChain->SetBranchAddress("g2dedx", &g2dedx, &b_g2dedx);
   fChain->SetBranchAddress("g2dedx1", &g2dedx1, &b_g2dedx1);
   fChain->SetBranchAddress("g2dedx_total", &g2dedx_total, &b_g2dedx_total);
   fChain->SetBranchAddress("g2dedx_total1", &g2dedx_total1, &b_g2dedx_total1);
   fChain->SetBranchAddress("gamma1_E", &gamma1_E, &b_gamma1_E);
   fChain->SetBranchAddress("gamma1_dEdx", &gamma1_dEdx, &b_gamma1_dEdx);
   fChain->SetBranchAddress("gamma1_dist_exit", &gamma1_dist_exit, &b_gamma1_dist_exit);
   fChain->SetBranchAddress("gamma1_dist_vtx", &gamma1_dist_vtx, &b_gamma1_dist_vtx);
   fChain->SetBranchAddress("gamma1_evis_ecal", &gamma1_evis_ecal, &b_gamma1_evis_ecal);
   fChain->SetBranchAddress("gamma1_evis_hcal", &gamma1_evis_hcal, &b_gamma1_evis_hcal);
   fChain->SetBranchAddress("gamma1_evis_scal", &gamma1_evis_scal, &b_gamma1_evis_scal);
   fChain->SetBranchAddress("gamma1_evis_trkr", &gamma1_evis_trkr, &b_gamma1_evis_trkr);
   fChain->SetBranchAddress("gamma1_phi", &gamma1_phi, &b_gamma1_phi);
   fChain->SetBranchAddress("gamma1_px", &gamma1_px, &b_gamma1_px);
   fChain->SetBranchAddress("gamma1_py", &gamma1_py, &b_gamma1_py);
   fChain->SetBranchAddress("gamma1_pz", &gamma1_pz, &b_gamma1_pz);
   fChain->SetBranchAddress("gamma1_theta", &gamma1_theta, &b_gamma1_theta);
   fChain->SetBranchAddress("gamma1_time", &gamma1_time, &b_gamma1_time);
   fChain->SetBranchAddress("gamma2_E", &gamma2_E, &b_gamma2_E);
   fChain->SetBranchAddress("gamma2_dEdx", &gamma2_dEdx, &b_gamma2_dEdx);
   fChain->SetBranchAddress("gamma2_dist_exit", &gamma2_dist_exit, &b_gamma2_dist_exit);
   fChain->SetBranchAddress("gamma2_dist_vtx", &gamma2_dist_vtx, &b_gamma2_dist_vtx);
   fChain->SetBranchAddress("gamma2_evis_ecal", &gamma2_evis_ecal, &b_gamma2_evis_ecal);
   fChain->SetBranchAddress("gamma2_evis_hcal", &gamma2_evis_hcal, &b_gamma2_evis_hcal);
   fChain->SetBranchAddress("gamma2_evis_scal", &gamma2_evis_scal, &b_gamma2_evis_scal);
   fChain->SetBranchAddress("gamma2_evis_trkr", &gamma2_evis_trkr, &b_gamma2_evis_trkr);
   fChain->SetBranchAddress("gamma2_phi", &gamma2_phi, &b_gamma2_phi);
   fChain->SetBranchAddress("gamma2_px", &gamma2_px, &b_gamma2_px);
   fChain->SetBranchAddress("gamma2_py", &gamma2_py, &b_gamma2_py);
   fChain->SetBranchAddress("gamma2_pz", &gamma2_pz, &b_gamma2_pz);
   fChain->SetBranchAddress("gamma2_theta", &gamma2_theta, &b_gamma2_theta);
   fChain->SetBranchAddress("gamma2_time", &gamma2_time, &b_gamma2_time);
   fChain->SetBranchAddress("hadronVisibleE", &hadronVisibleE, &b_hadronVisibleE);
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
   fChain->SetBranchAddress("pi0_E", &pi0_E, &b_pi0_E);
   fChain->SetBranchAddress("pi0_cos_openingAngle", &pi0_cos_openingAngle, &b_pi0_cos_openingAngle);
   fChain->SetBranchAddress("pi0_invMass", &pi0_invMass, &b_pi0_invMass);
   fChain->SetBranchAddress("pi0_openingAngle", &pi0_openingAngle, &b_pi0_openingAngle);
   fChain->SetBranchAddress("pi0_phi", &pi0_phi, &b_pi0_phi);
   fChain->SetBranchAddress("pi0_px", &pi0_px, &b_pi0_px);
   fChain->SetBranchAddress("pi0_py", &pi0_py, &b_pi0_py);
   fChain->SetBranchAddress("pi0_pz", &pi0_pz, &b_pi0_pz);
   fChain->SetBranchAddress("pi0_theta", &pi0_theta, &b_pi0_theta);
   fChain->SetBranchAddress("pi0_thetaX", &pi0_thetaX, &b_pi0_thetaX);
   fChain->SetBranchAddress("pi0_thetaY", &pi0_thetaY, &b_pi0_thetaY);
   fChain->SetBranchAddress("preFilter_rejectedEnergy", &preFilter_rejectedEnergy, &b_preFilter_rejectedEnergy);
   fChain->SetBranchAddress("prim_vtx_smallest_opening_angle", &prim_vtx_smallest_opening_angle, &b_prim_vtx_smallest_opening_angle);
   fChain->SetBranchAddress("time", &time, &b_time);
   fChain->SetBranchAddress("totalIDVisibleE", &totalIDVisibleE, &b_totalIDVisibleE);
   fChain->SetBranchAddress("totalODVisibleE", &totalODVisibleE, &b_totalODVisibleE);
   fChain->SetBranchAddress("totalVisibleE", &totalVisibleE, &b_totalVisibleE);
   fChain->SetBranchAddress("unattachedExtraE", &unattachedExtraE, &b_unattachedExtraE);
   fChain->SetBranchAddress("vtxBlobExtraE", &vtxBlobExtraE, &b_vtxBlobExtraE);
   fChain->SetBranchAddress("vtx_michel_distance", &vtx_michel_distance, &b_vtx_michel_distance);
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
   fChain->SetBranchAddress("g2dedx_cluster_energy_sz", &g2dedx_cluster_energy_sz, &b_g2dedx_cluster_energy_sz);
   fChain->SetBranchAddress("g2dedx_cluster_energy", g2dedx_cluster_energy, &b_g2dedx_cluster_energy);
   fChain->SetBranchAddress("g2dedx_rev_cluster_energy_sz", &g2dedx_rev_cluster_energy_sz, &b_g2dedx_rev_cluster_energy_sz);
   fChain->SetBranchAddress("g2dedx_rev_cluster_energy", g2dedx_rev_cluster_energy, &b_g2dedx_rev_cluster_energy);
   fChain->SetBranchAddress("gamma1_direction", gamma1_direction, &b_gamma1_direction);
   fChain->SetBranchAddress("gamma1_vertex", gamma1_vertex, &b_gamma1_vertex);
   fChain->SetBranchAddress("gamma2_direction", gamma2_direction, &b_gamma2_direction);
   fChain->SetBranchAddress("gamma2_vertex", gamma2_vertex, &b_gamma2_vertex);
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
   fChain->SetBranchAddress("truth_has_physics_event", &truth_has_physics_event, &b_truth_has_physics_event);
   fChain->SetBranchAddress("truth_isSignal", &truth_isSignal, &b_truth_isSignal);
   fChain->SetBranchAddress("truth_isSignalGold", &truth_isSignalGold, &b_truth_isSignalGold);
   fChain->SetBranchAddress("truth_isSignalSilver1", &truth_isSignalSilver1, &b_truth_isSignalSilver1);
   fChain->SetBranchAddress("truth_isSignal_1Pi0", &truth_isSignal_1Pi0, &b_truth_isSignal_1Pi0);
   fChain->SetBranchAddress("truth_isSignal_2Gamma", &truth_isSignal_2Gamma, &b_truth_isSignal_2Gamma);
   fChain->SetBranchAddress("truth_isFidVol", &truth_isFidVol, &b_truth_isFidVol);
   fChain->SetBranchAddress("truth_N_FSParticles", &truth_N_FSParticles, &b_truth_N_FSParticles);
   fChain->SetBranchAddress("truth_N_gamma", &truth_N_gamma, &b_truth_N_gamma);
   fChain->SetBranchAddress("truth_N_pi0", &truth_N_pi0, &b_truth_N_pi0);
   fChain->SetBranchAddress("truth_N_proton", &truth_N_proton, &b_truth_N_proton);
   fChain->SetBranchAddress("truth_N_proton_Long", &truth_N_proton_Long, &b_truth_N_proton_Long);
   fChain->SetBranchAddress("truth_N_proton_Short", &truth_N_proton_Short, &b_truth_N_proton_Short);
   fChain->SetBranchAddress("truth_muon_charge", &truth_muon_charge, &b_truth_muon_charge);
   fChain->SetBranchAddress("truth_pi0_GrandMother", &truth_pi0_GrandMother, &b_truth_pi0_GrandMother);
   fChain->SetBranchAddress("truth_pi0_GrandMotherStatus", &truth_pi0_GrandMotherStatus, &b_truth_pi0_GrandMotherStatus);
   fChain->SetBranchAddress("truth_pi0_Mother", &truth_pi0_Mother, &b_truth_pi0_Mother);
   fChain->SetBranchAddress("truth_pi0_MotherStatus", &truth_pi0_MotherStatus, &b_truth_pi0_MotherStatus);
   fChain->SetBranchAddress("truth_pi0_status", &truth_pi0_status, &b_truth_pi0_status);
   fChain->SetBranchAddress("truth_target_material", &truth_target_material, &b_truth_target_material);
   fChain->SetBranchAddress("truth_vertex_module", &truth_vertex_module, &b_truth_vertex_module);
   fChain->SetBranchAddress("truth_vertex_plane", &truth_vertex_plane, &b_truth_vertex_plane);
   fChain->SetBranchAddress("truth_muon_E", &truth_muon_E, &b_truth_muon_E);
   fChain->SetBranchAddress("truth_muon_px", &truth_muon_px, &b_truth_muon_px);
   fChain->SetBranchAddress("truth_muon_py", &truth_muon_py, &b_truth_muon_py);
   fChain->SetBranchAddress("truth_muon_pz", &truth_muon_pz, &b_truth_muon_pz);
   fChain->SetBranchAddress("truth_muon_theta_wrtbeam", &truth_muon_theta_wrtbeam, &b_truth_muon_theta_wrtbeam);
   fChain->SetBranchAddress("truth_pi0_E", &truth_pi0_E, &b_truth_pi0_E);
   fChain->SetBranchAddress("truth_pi0_px", &truth_pi0_px, &b_truth_pi0_px);
   fChain->SetBranchAddress("truth_pi0_py", &truth_pi0_py, &b_truth_pi0_py);
   fChain->SetBranchAddress("truth_pi0_pz", &truth_pi0_pz, &b_truth_pi0_pz);
   fChain->SetBranchAddress("truth_pi0_theta_wrtbeam", &truth_pi0_theta_wrtbeam, &b_truth_pi0_theta_wrtbeam);
   fChain->SetBranchAddress("truth_gamma_E", truth_gamma_E, &b_truth_gamma_E);
   fChain->SetBranchAddress("truth_gamma_px", truth_gamma_px, &b_truth_gamma_px);
   fChain->SetBranchAddress("truth_gamma_py", truth_gamma_py, &b_truth_gamma_py);
   fChain->SetBranchAddress("truth_gamma_pz", truth_gamma_pz, &b_truth_gamma_pz);
   fChain->SetBranchAddress("truth_gamma_theta_wrtbeam", truth_gamma_theta_wrtbeam, &b_truth_gamma_theta_wrtbeam);
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
   fChain->SetBranchAddress("truth_proton_E", truth_proton_E, &b_truth_proton_E);
   fChain->SetBranchAddress("truth_proton_px", truth_proton_px, &b_truth_proton_px);
   fChain->SetBranchAddress("truth_proton_py", truth_proton_py, &b_truth_proton_py);
   fChain->SetBranchAddress("truth_proton_pz", truth_proton_pz, &b_truth_proton_pz);
   fChain->SetBranchAddress("truth_proton_theta_wrtbeam", truth_proton_theta_wrtbeam, &b_truth_proton_theta_wrtbeam);
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
   fChain->SetBranchAddress("CCProtonPi0_r_minos_trk_vtx_plane", &CCProtonPi0_r_minos_trk_vtx_plane, &b_CCProtonPi0_r_minos_trk_vtx_plane);
   fChain->SetBranchAddress("CCProtonPi0_t_minos_trk_numFSMuons", &CCProtonPi0_t_minos_trk_numFSMuons, &b_CCProtonPi0_t_minos_trk_numFSMuons);
   fChain->SetBranchAddress("CCProtonPi0_t_minos_trk_primFSLeptonPDG", &CCProtonPi0_t_minos_trk_primFSLeptonPDG, &b_CCProtonPi0_t_minos_trk_primFSLeptonPDG);
   fChain->SetBranchAddress("CCProtonPi0_trajMuonProngPDG", &CCProtonPi0_trajMuonProngPDG, &b_CCProtonPi0_trajMuonProngPDG);
   fChain->SetBranchAddress("CCProtonPi0_trajMuonProngPrimary", &CCProtonPi0_trajMuonProngPrimary, &b_CCProtonPi0_trajMuonProngPrimary);
   fChain->SetBranchAddress("CCProtonPi0_vtx_module", &CCProtonPi0_vtx_module, &b_CCProtonPi0_vtx_module);
   fChain->SetBranchAddress("CCProtonPi0_vtx_plane", &CCProtonPi0_vtx_plane, &b_CCProtonPi0_vtx_plane);
   fChain->SetBranchAddress("CCProtonPi0_QSq", &CCProtonPi0_QSq, &b_CCProtonPi0_QSq);
   fChain->SetBranchAddress("CCProtonPi0_WSq", &CCProtonPi0_WSq, &b_CCProtonPi0_WSq);
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
   fChain->SetBranchAddress("CCProtonPi0_neutrino_E", &CCProtonPi0_neutrino_E, &b_CCProtonPi0_neutrino_E);
   fChain->SetBranchAddress("CCProtonPi0_neutrino_E_Cal", &CCProtonPi0_neutrino_E_Cal, &b_CCProtonPi0_neutrino_E_Cal);
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
   fChain->SetBranchAddress("CCProtonPi0_proton_isRecoGood", CCProtonPi0_proton_isRecoGood, &b_CCProtonPi0_proton_isRecoGood);
   fChain->SetBranchAddress("CCProtonPi0_proton_kinked", CCProtonPi0_proton_kinked, &b_CCProtonPi0_proton_kinked);
   fChain->SetBranchAddress("CCProtonPi0_proton_odMatch", CCProtonPi0_proton_odMatch, &b_CCProtonPi0_proton_odMatch);
   fChain->SetBranchAddress("CCProtonPi0_proton_trk_pat_history", CCProtonPi0_proton_trk_pat_history, &b_CCProtonPi0_proton_trk_pat_history);
   fChain->SetBranchAddress("CCProtonPi0_trajProtonProngPDG", CCProtonPi0_trajProtonProngPDG, &b_CCProtonPi0_trajProtonProngPDG);
   fChain->SetBranchAddress("CCProtonPi0_trajProtonProngPrimary", CCProtonPi0_trajProtonProngPrimary, &b_CCProtonPi0_trajProtonProngPrimary);
   fChain->SetBranchAddress("CCProtonPi0_endProtonTrajMomentum", CCProtonPi0_endProtonTrajMomentum, &b_CCProtonPi0_endProtonTrajMomentum);
   fChain->SetBranchAddress("CCProtonPi0_endProtonTrajXPosition", CCProtonPi0_endProtonTrajXPosition, &b_CCProtonPi0_endProtonTrajXPosition);
   fChain->SetBranchAddress("CCProtonPi0_endProtonTrajYPosition", CCProtonPi0_endProtonTrajYPosition, &b_CCProtonPi0_endProtonTrajYPosition);
   fChain->SetBranchAddress("CCProtonPi0_endProtonTrajZPosition", CCProtonPi0_endProtonTrajZPosition, &b_CCProtonPi0_endProtonTrajZPosition);
   fChain->SetBranchAddress("CCProtonPi0_pionScore", CCProtonPi0_pionScore, &b_CCProtonPi0_pionScore);
   fChain->SetBranchAddress("CCProtonPi0_protonScore", CCProtonPi0_protonScore, &b_CCProtonPi0_protonScore);
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

Bool_t MC_Sample::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MC_Sample::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MC_Sample::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MC_Sample_cxx
