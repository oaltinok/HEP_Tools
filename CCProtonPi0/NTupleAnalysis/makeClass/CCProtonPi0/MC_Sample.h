//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sat Feb  6 11:37:49 2016 by ROOT version 5.34/05
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
   Bool_t          is_blobs_recovered;
   Bool_t          is_blobs_recovered_direction;
   Bool_t          is_blobs_recovered_invMass;
   Bool_t          is_blobs_recovered_small_angle;
   Bool_t          is_blobs_recovered_search_view_U;
   Bool_t          is_blobs_recovered_search_view_V;
   Int_t           Cut_BlobDirectionBad;
   Int_t           Cut_ConeBlobs;
   Int_t           Cut_EndPoint_Michel_Exist;
   Int_t           Cut_Muon_Charge;
   Int_t           Cut_Muon_None;
   Int_t           Cut_Particle_None;
   Int_t           Cut_Pi0_Bad;
   Int_t           Cut_PreFilter_Pi0;
   Int_t           Cut_Proton_Bad;
   Int_t           Cut_Proton_None;
   Int_t           Cut_Vertex_Michel_Exist;
   Int_t           Cut_Vertex_None;
   Int_t           Cut_Vertex_Not_Fiducial;
   Int_t           Cut_Vertex_Not_Reconstructable;
   Int_t           Cut_secEndPoint_Michel_Exist;
   Int_t           OneShower_nClusters;
   Int_t           ThreeShower_s1_nClusters;
   Int_t           ThreeShower_s2_nClusters;
   Int_t           ThreeShower_s3_nClusters;
   Int_t           anglescan_ncand;
   Int_t           anglescan_ncandx;
   Int_t           anglescan_nfoundBlobs;
   Int_t           broken_track_most_us_plane;
   Int_t           g1blob_1ParFit_ndof;
   Int_t           g1dedx_doublet;
   Int_t           g1dedx_empty_plane;
   Int_t           g1dedx_nplane;
   Int_t           g2blob_1ParFit_ndof;
   Int_t           g2dedx_doublet;
   Int_t           g2dedx_empty_plane;
   Int_t           g2dedx_nplane;
   Int_t           gamma1_blob_nclusters;
   Int_t           gamma1_blob_ndigits;
   Int_t           gamma2_blob_nclusters;
   Int_t           gamma2_blob_ndigits;
   Int_t           muon_N_minosTracks;
   Int_t           muon_charge;
   Int_t           muon_hasMinosMatchStub;
   Int_t           muon_hasMinosMatchTrack;
   Int_t           muon_minervaTrack_types;
   Int_t           muon_minosTrackQuality;
   Int_t           muon_roadUpstreamPlanes;
   Int_t           nProtonCandidates;
   Int_t           nTracks;
   Int_t           nTracks_Close;
   Int_t           nTracks_Discarded;
   Int_t           nTracks_Far;
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
   Int_t           proton_kinked;
   Int_t           proton_leadingIndice;
   Int_t           vtx_fit_converged;
   Int_t           vtx_module;
   Int_t           vtx_plane;
   Int_t           vtx_primary_index;
   Int_t           vtx_primary_multiplicity;
   Int_t           vtx_secondary_count;
   Int_t           vtx_total_count;
   Double_t        ConeBlobs_usable_evis_Tracker;
   Double_t        Coneblobs_usable_evis_ECAL;
   Double_t        Coneblobs_usable_evis_HCAL;
   Double_t        Extra_Energy_Dispersed;
   Double_t        Extra_Energy_Muon;
   Double_t        Extra_Energy_Rejected;
   Double_t        Extra_Evis_Leftover;
   Double_t        OneShower_dist_vtx;
   Double_t        OneShower_energy;
   Double_t        OneShower_theta;
   Double_t        ThreeShower_s1_dist_vtx;
   Double_t        ThreeShower_s1_energy;
   Double_t        ThreeShower_s1_theta;
   Double_t        ThreeShower_s2_dist_vtx;
   Double_t        ThreeShower_s2_energy;
   Double_t        ThreeShower_s2_theta;
   Double_t        ThreeShower_s3_dist_vtx;
   Double_t        ThreeShower_s3_energy;
   Double_t        ThreeShower_s3_theta;
   Double_t        energy_from_mc;
   Double_t        energy_from_mc_fraction;
   Double_t        energy_from_mc_fraction_of_highest;
   Double_t        g1blob_1ParFit_fval;
   Double_t        g1blob_2ParFit_vtx_distance;
   Double_t        g1dedx;
   Double_t        g1dedx1;
   Double_t        g1dedx_total;
   Double_t        g1dedx_total1;
   Double_t        g2blob_1ParFit_fval;
   Double_t        g2blob_2ParFit_vtx_distance;
   Double_t        g2dedx;
   Double_t        g2dedx1;
   Double_t        g2dedx_total;
   Double_t        g2dedx_total1;
   Double_t        gamma1_E;
   Double_t        gamma1_E_Old;
   Double_t        gamma1_P;
   Double_t        gamma1_blob_energy;
   Double_t        gamma1_blob_minsep;
   Double_t        gamma1_dEdx;
   Double_t        gamma1_dist_vtx;
   Double_t        gamma1_energy_ecal;
   Double_t        gamma1_energy_hcal;
   Double_t        gamma1_energy_scal_UV;
   Double_t        gamma1_energy_scal_X;
   Double_t        gamma1_energy_trkr;
   Double_t        gamma1_evis_ecal;
   Double_t        gamma1_evis_hcal;
   Double_t        gamma1_evis_scal_UV;
   Double_t        gamma1_evis_scal_X;
   Double_t        gamma1_evis_trkr;
   Double_t        gamma1_phi;
   Double_t        gamma1_px;
   Double_t        gamma1_py;
   Double_t        gamma1_pz;
   Double_t        gamma1_theta;
   Double_t        gamma1_time;
   Double_t        gamma2_E;
   Double_t        gamma2_E_Old;
   Double_t        gamma2_P;
   Double_t        gamma2_blob_energy;
   Double_t        gamma2_blob_minsep;
   Double_t        gamma2_dEdx;
   Double_t        gamma2_dist_vtx;
   Double_t        gamma2_energy_ecal;
   Double_t        gamma2_energy_hcal;
   Double_t        gamma2_energy_scal_UV;
   Double_t        gamma2_energy_scal_X;
   Double_t        gamma2_energy_trkr;
   Double_t        gamma2_evis_ecal;
   Double_t        gamma2_evis_hcal;
   Double_t        gamma2_evis_scal_UV;
   Double_t        gamma2_evis_scal_X;
   Double_t        gamma2_evis_trkr;
   Double_t        gamma2_phi;
   Double_t        gamma2_px;
   Double_t        gamma2_py;
   Double_t        gamma2_pz;
   Double_t        gamma2_theta;
   Double_t        gamma2_time;
   Double_t        michelProng_begin_Z;
   Double_t        michelProng_distance;
   Double_t        michelProng_end_Z;
   Double_t        michelProng_energy;
   Double_t        michelProng_time_diff;
   Double_t        muon_E;
   Double_t        muon_E_shift;
   Double_t        muon_KE;
   Double_t        muon_P;
   Double_t        muon_muScore;
   Double_t        muon_phi;
   Double_t        muon_px;
   Double_t        muon_py;
   Double_t        muon_pz;
   Double_t        muon_qp;
   Double_t        muon_qpqpe;
   Double_t        muon_roadUpstreamEnergy;
   Double_t        muon_theta;
   Double_t        muon_thetaX;
   Double_t        muon_thetaY;
   Double_t        muon_theta_biasDown;
   Double_t        muon_theta_biasUp;
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
   Double_t        pi0_KE;
   Double_t        pi0_P;
   Double_t        pi0_cos_openingAngle;
   Double_t        pi0_invMass;
   Double_t        pi0_invMass_Old;
   Double_t        pi0_openingAngle;
   Double_t        pi0_phi;
   Double_t        pi0_px;
   Double_t        pi0_py;
   Double_t        pi0_pz;
   Double_t        pi0_theta;
   Double_t        pi0_thetaX;
   Double_t        pi0_thetaY;
   Double_t        preFilter_evis_ECAL;
   Double_t        preFilter_evis_HCAL;
   Double_t        preFilter_evis_NuclearTarget;
   Double_t        preFilter_evis_TotalExceptNuclearTarget;
   Double_t        preFilter_evis_Tracker;
   Double_t        preFilter_evis_nearvtx;
   Double_t        preFilter_evis_total;
   Double_t        preFilter_rejectedEnergy;
   Double_t        prim_vtx_smallest_opening_angle;
   Double_t        proton_E;
   Double_t        proton_KE;
   Double_t        proton_LLRScore;
   Double_t        proton_P;
   Double_t        proton_length;
   Double_t        proton_phi;
   Double_t        proton_pionScore;
   Double_t        proton_protonScore;
   Double_t        proton_px;
   Double_t        proton_py;
   Double_t        proton_pz;
   Double_t        proton_theta;
   Double_t        proton_thetaX;
   Double_t        proton_thetaY;
   Double_t        reco_eventID;
   Double_t        time;
   Double_t        vertex_blob_evis;
   Double_t        vtx_fit_chi2;
   Double_t        vtx_x;
   Double_t        vtx_y;
   Double_t        vtx_z;
   Int_t           all_protons_kinked[10];
   Int_t           all_protons_odMatch[10];
   Int_t           g1dedx_cluster_occupancy_sz;
   Int_t           g1dedx_cluster_occupancy[6];   //[g1dedx_cluster_occupancy_sz]
   Int_t           g2dedx_cluster_occupancy_sz;
   Int_t           g2dedx_cluster_occupancy[6];   //[g2dedx_cluster_occupancy_sz]
   Int_t           nTracks_Secondary_Vtx_sz;
   Int_t           nTracks_Secondary_Vtx[1];   //[nTracks_Secondary_Vtx_sz]
   Double_t        all_protons_E[10];
   Double_t        all_protons_KE[10];
   Double_t        all_protons_LLRScore[10];
   Double_t        all_protons_P[10];
   Double_t        all_protons_chi2_ndf[10];
   Double_t        all_protons_endPointX[10];
   Double_t        all_protons_endPointY[10];
   Double_t        all_protons_endPointZ[10];
   Double_t        all_protons_length[10];
   Double_t        all_protons_p_calCorrection[10];
   Double_t        all_protons_p_dEdXTool[10];
   Double_t        all_protons_p_visEnergy[10];
   Double_t        all_protons_phi[10];
   Double_t        all_protons_pionScore[10];
   Double_t        all_protons_protonScore[10];
   Double_t        all_protons_px[10];
   Double_t        all_protons_py[10];
   Double_t        all_protons_pz[10];
   Double_t        all_protons_startPointX[10];
   Double_t        all_protons_startPointY[10];
   Double_t        all_protons_startPointZ[10];
   Double_t        all_protons_theta[10];
   Double_t        all_protons_thetaX[10];
   Double_t        all_protons_thetaY[10];
   Double_t        fit_vtx[3];
   Int_t           g1dedx_cluster_energy_sz;
   Double_t        g1dedx_cluster_energy[6];   //[g1dedx_cluster_energy_sz]
   Int_t           g1dedx_rev_cluster_energy_sz;
   Double_t        g1dedx_rev_cluster_energy[92];   //[g1dedx_rev_cluster_energy_sz]
   Int_t           g2dedx_cluster_energy_sz;
   Double_t        g2dedx_cluster_energy[6];   //[g2dedx_cluster_energy_sz]
   Int_t           g2dedx_rev_cluster_energy_sz;
   Double_t        g2dedx_rev_cluster_energy[36];   //[g2dedx_rev_cluster_energy_sz]
   Double_t        gamma1_direction[3];
   Double_t        gamma1_vertex[3];
   Double_t        gamma2_direction[3];
   Double_t        gamma2_vertex[3];
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
   Bool_t          truth_isGamma1_conv_inside;
   Bool_t          truth_isGamma2_conv_inside;
   Bool_t          truth_isSignal;
   Bool_t          truth_isFidVol;
   Bool_t          truth_isNC;
   Bool_t          truth_ReconstructEvent;
   Bool_t          truth_isBckg_NoPi0;
   Bool_t          truth_isBckg_SinglePi0;
   Bool_t          truth_isBckg_MultiPi0;
   Bool_t          truth_isBckg_NC;
   Bool_t          truth_isBckg_AntiNeutrino;
   Bool_t          truth_isBckg_QELike;
   Bool_t          truth_isBckg_SingleChargedPion;
   Bool_t          truth_isBckg_SingleChargedPion_ChargeExchanged;
   Bool_t          truth_isBckg_DoublePionWithPi0;
   Bool_t          truth_isBckg_DoublePionWithoutPi0;
   Bool_t          truth_isBckg_MultiPionWithPi0;
   Bool_t          truth_isBckg_MultiPionWithoutPi0;
   Bool_t          truth_isBckg_Other;
   Bool_t          truth_isBckg_withMichel;
   Int_t           truth_Bckg_nOther;
   Int_t           truth_Bckg_nPi0_Primary;
   Int_t           truth_Bckg_nPi0_Secondary;
   Int_t           truth_Bckg_nPi0_Total;
   Int_t           truth_Bckg_nPiCharged;
   Int_t           truth_Bckg_nPiCharged_ChargeExchanged;
   Int_t           truth_N_FSParticles;
   Int_t           truth_N_other;
   Int_t           truth_N_pi0;
   Int_t           truth_N_proton;
   Int_t           truth_N_trueMichelElectrons;
   Int_t           truth_OneShower_evis_most_pdg;
   Int_t           truth_Rejected_unused_evis_most_pdg;
   Int_t           truth_ThreeShower_s1_evis_most_pdg;
   Int_t           truth_ThreeShower_s2_evis_most_pdg;
   Int_t           truth_ThreeShower_s3_evis_most_pdg;
   Int_t           truth_blob1_evis_most_pdg;
   Int_t           truth_blob2_evis_most_pdg;
   Int_t           truth_dispersed_unused_evis_most_pdg;
   Int_t           truth_pi0_GrandMother;
   Int_t           truth_pi0_GrandMotherStatus;
   Int_t           truth_pi0_Mother;
   Int_t           truth_pi0_MotherStatus;
   Int_t           truth_pi0_status;
   Int_t           truth_target_material;
   Int_t           truth_vertex_module;
   Int_t           truth_vertex_plane;
   Int_t           truth_vertex_unused_evis_most_pdg;
   Double_t        truth_OneShower_evis_muon;
   Double_t        truth_OneShower_evis_neutron;
   Double_t        truth_OneShower_evis_piminus;
   Double_t        truth_OneShower_evis_piplus;
   Double_t        truth_OneShower_evis_pizero;
   Double_t        truth_OneShower_evis_proton;
   Double_t        truth_OneShower_evis_total_norm;
   Double_t        truth_OneShower_evis_total_truth;
   Double_t        truth_Rejected_unused_evis_total_norm;
   Double_t        truth_Rejected_unused_evis_total_truth;
   Double_t        truth_ThreeShower_s1_evis_muon;
   Double_t        truth_ThreeShower_s1_evis_neutron;
   Double_t        truth_ThreeShower_s1_evis_piminus;
   Double_t        truth_ThreeShower_s1_evis_piplus;
   Double_t        truth_ThreeShower_s1_evis_pizero;
   Double_t        truth_ThreeShower_s1_evis_proton;
   Double_t        truth_ThreeShower_s1_evis_total_norm;
   Double_t        truth_ThreeShower_s1_evis_total_truth;
   Double_t        truth_ThreeShower_s2_evis_muon;
   Double_t        truth_ThreeShower_s2_evis_neutron;
   Double_t        truth_ThreeShower_s2_evis_piminus;
   Double_t        truth_ThreeShower_s2_evis_piplus;
   Double_t        truth_ThreeShower_s2_evis_pizero;
   Double_t        truth_ThreeShower_s2_evis_proton;
   Double_t        truth_ThreeShower_s2_evis_total_norm;
   Double_t        truth_ThreeShower_s2_evis_total_truth;
   Double_t        truth_ThreeShower_s3_evis_muon;
   Double_t        truth_ThreeShower_s3_evis_neutron;
   Double_t        truth_ThreeShower_s3_evis_piminus;
   Double_t        truth_ThreeShower_s3_evis_piplus;
   Double_t        truth_ThreeShower_s3_evis_pizero;
   Double_t        truth_ThreeShower_s3_evis_proton;
   Double_t        truth_ThreeShower_s3_evis_total_norm;
   Double_t        truth_ThreeShower_s3_evis_total_truth;
   Double_t        truth_allClusters_evis_pizero;
   Double_t        truth_blob1_evis_muon;
   Double_t        truth_blob1_evis_neutron;
   Double_t        truth_blob1_evis_piminus;
   Double_t        truth_blob1_evis_piplus;
   Double_t        truth_blob1_evis_pizero;
   Double_t        truth_blob1_evis_proton;
   Double_t        truth_blob1_evis_total_norm;
   Double_t        truth_blob1_evis_total_truth;
   Double_t        truth_blob2_evis_muon;
   Double_t        truth_blob2_evis_neutron;
   Double_t        truth_blob2_evis_piminus;
   Double_t        truth_blob2_evis_piplus;
   Double_t        truth_blob2_evis_pizero;
   Double_t        truth_blob2_evis_proton;
   Double_t        truth_blob2_evis_total_norm;
   Double_t        truth_blob2_evis_total_truth;
   Double_t        truth_dispersed_unused_evis_gamma;
   Double_t        truth_dispersed_unused_evis_muon;
   Double_t        truth_dispersed_unused_evis_neutron;
   Double_t        truth_dispersed_unused_evis_piminus;
   Double_t        truth_dispersed_unused_evis_piplus;
   Double_t        truth_dispersed_unused_evis_pizero;
   Double_t        truth_dispersed_unused_evis_proton;
   Double_t        truth_dispersed_unused_evis_total_norm;
   Double_t        truth_dispersed_unused_evis_total_truth;
   Double_t        truth_ecal_unused_evis_muon;
   Double_t        truth_ecal_unused_evis_neutron;
   Double_t        truth_ecal_unused_evis_piminus;
   Double_t        truth_ecal_unused_evis_piplus;
   Double_t        truth_ecal_unused_evis_pizero;
   Double_t        truth_ecal_unused_evis_proton;
   Double_t        truth_ecal_unused_evis_total_norm;
   Double_t        truth_ecal_unused_evis_total_truth;
   Double_t        truth_eventID;
   Double_t        truth_hcal_unused_evis_muon;
   Double_t        truth_hcal_unused_evis_neutron;
   Double_t        truth_hcal_unused_evis_piminus;
   Double_t        truth_hcal_unused_evis_piplus;
   Double_t        truth_hcal_unused_evis_pizero;
   Double_t        truth_hcal_unused_evis_proton;
   Double_t        truth_hcal_unused_evis_total_norm;
   Double_t        truth_hcal_unused_evis_total_truth;
   Double_t        truth_michelElectron_E;
   Double_t        truth_michelElectron_P;
   Double_t        truth_michelMuon_P;
   Double_t        truth_michelMuon_end_dist_vtx;
   Double_t        truth_michelMuon_length;
   Double_t        truth_michelPion_P;
   Double_t        truth_michelPion_begin_dist_vtx;
   Double_t        truth_michelPion_length;
   Double_t        truth_muon_P;
   Double_t        truth_muon_theta;
   Double_t        truth_other_unused_evis_muon;
   Double_t        truth_other_unused_evis_neutron;
   Double_t        truth_other_unused_evis_piminus;
   Double_t        truth_other_unused_evis_piplus;
   Double_t        truth_other_unused_evis_pizero;
   Double_t        truth_other_unused_evis_proton;
   Double_t        truth_other_unused_evis_total_norm;
   Double_t        truth_other_unused_evis_total_truth;
   Double_t        truth_pi0_P;
   Double_t        truth_pi0_theta;
   Double_t        truth_proton_P;
   Double_t        truth_proton_theta;
   Double_t        truth_total_captured_evis_pizero;
   Double_t        truth_total_captured_evis_total_norm;
   Double_t        truth_total_captured_evis_total_truth;
   Double_t        truth_vertex_unused_evis_gamma;
   Double_t        truth_vertex_unused_evis_muon;
   Double_t        truth_vertex_unused_evis_neutron;
   Double_t        truth_vertex_unused_evis_piminus;
   Double_t        truth_vertex_unused_evis_piplus;
   Double_t        truth_vertex_unused_evis_pizero;
   Double_t        truth_vertex_unused_evis_proton;
   Double_t        truth_vertex_unused_evis_total_norm;
   Double_t        truth_vertex_unused_evis_total_truth;
   Double_t        truth_gamma1_4P[4];
   Double_t        truth_gamma1_final_pos[3];
   Double_t        truth_gamma1_init_pos[3];
   Double_t        truth_gamma2_4P[4];
   Double_t        truth_gamma2_final_pos[3];
   Double_t        truth_gamma2_init_pos[3];
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
   Double_t        truth_michelMuon_endPoint[3];
   Double_t        truth_muon_4P[4];
   Double_t        truth_pi0_4P[4];
   Double_t        truth_proton_4P[4];
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
   Int_t           CCProtonPi0_ntrajMuonProng;
   Int_t           CCProtonPi0_r_minos_trk_vtx_plane;
   Int_t           CCProtonPi0_t_minos_trk_numFSMuons;
   Int_t           CCProtonPi0_t_minos_trk_primFSLeptonPDG;
   Int_t           CCProtonPi0_trajMuonProngPDG;
   Int_t           CCProtonPi0_trajMuonProngPrimary;
   Double_t        CCProtonPi0_Extra_Energy_Total;
   Double_t        CCProtonPi0_QSq;
   Double_t        CCProtonPi0_QSq_1Track_Alt;
   Double_t        CCProtonPi0_WSq;
   Double_t        CCProtonPi0_WSq_1Track_Alt;
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
   Double_t        CCProtonPi0_neutrino_E;
   Double_t        CCProtonPi0_neutrino_E_1Track_Alt;
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
   Double_t        CCProtonPi0_vertex_energy;
   Int_t           CCProtonPi0_isProtonInsideOD[10];
   Int_t           CCProtonPi0_ntrajProtonProng[10];
   Int_t           CCProtonPi0_trajProtonProngPDG[10];
   Int_t           CCProtonPi0_trajProtonProngPrimary[10];
   Double_t        CCProtonPi0_endProtonTrajMomentum[10];
   Double_t        CCProtonPi0_endProtonTrajXPosition[10];
   Double_t        CCProtonPi0_endProtonTrajYPosition[10];
   Double_t        CCProtonPi0_endProtonTrajZPosition[10];
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
   Double_t        mc_FSPartPx[10];   //[mc_nFSPart]
   Double_t        mc_FSPartPy[10];   //[mc_nFSPart]
   Double_t        mc_FSPartPz[10];   //[mc_nFSPart]
   Double_t        mc_FSPartE[10];   //[mc_nFSPart]
   Int_t           mc_FSPartPDG[10];   //[mc_nFSPart]
   Int_t           mc_er_nPart;
   Int_t           mc_er_ID[21];   //[mc_er_nPart]
   Int_t           mc_er_status[21];   //[mc_er_nPart]
   Double_t        mc_er_posInNucX[21];   //[mc_er_nPart]
   Double_t        mc_er_posInNucY[21];   //[mc_er_nPart]
   Double_t        mc_er_posInNucZ[21];   //[mc_er_nPart]
   Double_t        mc_er_Px[21];   //[mc_er_nPart]
   Double_t        mc_er_Py[21];   //[mc_er_nPart]
   Double_t        mc_er_Pz[21];   //[mc_er_nPart]
   Double_t        mc_er_E[21];   //[mc_er_nPart]
   Int_t           mc_er_FD[21];   //[mc_er_nPart]
   Int_t           mc_er_LD[21];   //[mc_er_nPart]
   Int_t           mc_er_mother[21];   //[mc_er_nPart]
   Int_t           mc_fr_nNuAncestorIDs;
   Int_t           mc_fr_nuAncestorIDs[6];   //[mc_fr_nNuAncestorIDs]
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
   Double_t        mc_ppfx1_cvweight;
   Double_t        mc_gen1_cvweight_total;
   Double_t        gen1_wgt;
   Double_t        mc_gen1_cvweight_totalFlux;
   Double_t        mc_gen1_cvweight_NA49;
   Int_t           mc_wgt_Flux_BeamFocus_sz;
   Double_t        mc_wgt_Flux_BeamFocus[100];   //[mc_wgt_Flux_BeamFocus_sz]
   Int_t           mc_wgt_gen1_Flux_Tertiary_sz;
   Double_t        mc_wgt_gen1_Flux_Tertiary[100];   //[mc_wgt_gen1_Flux_Tertiary_sz]
   Int_t           mc_wgt_gen1_Flux_NA49_sz;
   Double_t        mc_wgt_gen1_Flux_NA49[100];   //[mc_wgt_gen1_Flux_NA49_sz]
   Int_t           mc_wgt_Norm_sz;
   Double_t        mc_wgt_Norm[1];   //[mc_wgt_Norm_sz]
   Int_t           mc_wgt_ppfx1_Total_sz;
   Double_t        mc_wgt_ppfx1_Total[100];   //[mc_wgt_ppfx1_Total_sz]
   Int_t           n_prongs;
   Int_t           prong_nParticles[4];   //[n_prongs]
   Double_t        prong_part_score[4];   //[n_prongs]
   Double_t        prong_part_mass[4];   //[n_prongs]
   Int_t           prong_part_charge[4];   //[n_prongs]
   Int_t           prong_part_pid[4];   //[n_prongs]
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
   TBranch        *b_is_blobs_recovered;   //!
   TBranch        *b_is_blobs_recovered_direction;   //!
   TBranch        *b_is_blobs_recovered_invMass;   //!
   TBranch        *b_is_blobs_recovered_small_angle;   //!
   TBranch        *b_is_blobs_recovered_search_view_U;   //!
   TBranch        *b_is_blobs_recovered_search_view_V;   //!
   TBranch        *b_Cut_BlobDirectionBad;   //!
   TBranch        *b_Cut_ConeBlobs;   //!
   TBranch        *b_Cut_EndPoint_Michel_Exist;   //!
   TBranch        *b_Cut_Muon_Charge;   //!
   TBranch        *b_Cut_Muon_None;   //!
   TBranch        *b_Cut_Particle_None;   //!
   TBranch        *b_Cut_Pi0_Bad;   //!
   TBranch        *b_Cut_PreFilter_Pi0;   //!
   TBranch        *b_Cut_Proton_Bad;   //!
   TBranch        *b_Cut_Proton_None;   //!
   TBranch        *b_Cut_Vertex_Michel_Exist;   //!
   TBranch        *b_Cut_Vertex_None;   //!
   TBranch        *b_Cut_Vertex_Not_Fiducial;   //!
   TBranch        *b_Cut_Vertex_Not_Reconstructable;   //!
   TBranch        *b_Cut_secEndPoint_Michel_Exist;   //!
   TBranch        *b_OneShower_nClusters;   //!
   TBranch        *b_ThreeShower_s1_nClusters;   //!
   TBranch        *b_ThreeShower_s2_nClusters;   //!
   TBranch        *b_ThreeShower_s3_nClusters;   //!
   TBranch        *b_anglescan_ncand;   //!
   TBranch        *b_anglescan_ncandx;   //!
   TBranch        *b_anglescan_nfoundBlobs;   //!
   TBranch        *b_broken_track_most_us_plane;   //!
   TBranch        *b_g1blob_1ParFit_ndof;   //!
   TBranch        *b_g1dedx_doublet;   //!
   TBranch        *b_g1dedx_empty_plane;   //!
   TBranch        *b_g1dedx_nplane;   //!
   TBranch        *b_g2blob_1ParFit_ndof;   //!
   TBranch        *b_g2dedx_doublet;   //!
   TBranch        *b_g2dedx_empty_plane;   //!
   TBranch        *b_g2dedx_nplane;   //!
   TBranch        *b_gamma1_blob_nclusters;   //!
   TBranch        *b_gamma1_blob_ndigits;   //!
   TBranch        *b_gamma2_blob_nclusters;   //!
   TBranch        *b_gamma2_blob_ndigits;   //!
   TBranch        *b_muon_N_minosTracks;   //!
   TBranch        *b_muon_charge;   //!
   TBranch        *b_muon_hasMinosMatchStub;   //!
   TBranch        *b_muon_hasMinosMatchTrack;   //!
   TBranch        *b_muon_minervaTrack_types;   //!
   TBranch        *b_muon_minosTrackQuality;   //!
   TBranch        *b_muon_roadUpstreamPlanes;   //!
   TBranch        *b_nProtonCandidates;   //!
   TBranch        *b_nTracks;   //!
   TBranch        *b_nTracks_Close;   //!
   TBranch        *b_nTracks_Discarded;   //!
   TBranch        *b_nTracks_Far;   //!
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
   TBranch        *b_proton_kinked;   //!
   TBranch        *b_proton_leadingIndice;   //!
   TBranch        *b_vtx_fit_converged;   //!
   TBranch        *b_vtx_module;   //!
   TBranch        *b_vtx_plane;   //!
   TBranch        *b_vtx_primary_index;   //!
   TBranch        *b_vtx_primary_multiplicity;   //!
   TBranch        *b_vtx_secondary_count;   //!
   TBranch        *b_vtx_total_count;   //!
   TBranch        *b_ConeBlobs_usable_evis_Tracker;   //!
   TBranch        *b_Coneblobs_usable_evis_ECAL;   //!
   TBranch        *b_Coneblobs_usable_evis_HCAL;   //!
   TBranch        *b_Extra_Energy_Dispersed;   //!
   TBranch        *b_Extra_Energy_Muon;   //!
   TBranch        *b_Extra_Energy_Rejected;   //!
   TBranch        *b_Extra_Evis_Leftover;   //!
   TBranch        *b_OneShower_dist_vtx;   //!
   TBranch        *b_OneShower_energy;   //!
   TBranch        *b_OneShower_theta;   //!
   TBranch        *b_ThreeShower_s1_dist_vtx;   //!
   TBranch        *b_ThreeShower_s1_energy;   //!
   TBranch        *b_ThreeShower_s1_theta;   //!
   TBranch        *b_ThreeShower_s2_dist_vtx;   //!
   TBranch        *b_ThreeShower_s2_energy;   //!
   TBranch        *b_ThreeShower_s2_theta;   //!
   TBranch        *b_ThreeShower_s3_dist_vtx;   //!
   TBranch        *b_ThreeShower_s3_energy;   //!
   TBranch        *b_ThreeShower_s3_theta;   //!
   TBranch        *b_energy_from_mc;   //!
   TBranch        *b_energy_from_mc_fraction;   //!
   TBranch        *b_energy_from_mc_fraction_of_highest;   //!
   TBranch        *b_g1blob_1ParFit_fval;   //!
   TBranch        *b_g1blob_2ParFit_vtx_distance;   //!
   TBranch        *b_g1dedx;   //!
   TBranch        *b_g1dedx1;   //!
   TBranch        *b_g1dedx_total;   //!
   TBranch        *b_g1dedx_total1;   //!
   TBranch        *b_g2blob_1ParFit_fval;   //!
   TBranch        *b_g2blob_2ParFit_vtx_distance;   //!
   TBranch        *b_g2dedx;   //!
   TBranch        *b_g2dedx1;   //!
   TBranch        *b_g2dedx_total;   //!
   TBranch        *b_g2dedx_total1;   //!
   TBranch        *b_gamma1_E;   //!
   TBranch        *b_gamma1_E_Old;   //!
   TBranch        *b_gamma1_P;   //!
   TBranch        *b_gamma1_blob_energy;   //!
   TBranch        *b_gamma1_blob_minsep;   //!
   TBranch        *b_gamma1_dEdx;   //!
   TBranch        *b_gamma1_dist_vtx;   //!
   TBranch        *b_gamma1_energy_ecal;   //!
   TBranch        *b_gamma1_energy_hcal;   //!
   TBranch        *b_gamma1_energy_scal_UV;   //!
   TBranch        *b_gamma1_energy_scal_X;   //!
   TBranch        *b_gamma1_energy_trkr;   //!
   TBranch        *b_gamma1_evis_ecal;   //!
   TBranch        *b_gamma1_evis_hcal;   //!
   TBranch        *b_gamma1_evis_scal_UV;   //!
   TBranch        *b_gamma1_evis_scal_X;   //!
   TBranch        *b_gamma1_evis_trkr;   //!
   TBranch        *b_gamma1_phi;   //!
   TBranch        *b_gamma1_px;   //!
   TBranch        *b_gamma1_py;   //!
   TBranch        *b_gamma1_pz;   //!
   TBranch        *b_gamma1_theta;   //!
   TBranch        *b_gamma1_time;   //!
   TBranch        *b_gamma2_E;   //!
   TBranch        *b_gamma2_E_Old;   //!
   TBranch        *b_gamma2_P;   //!
   TBranch        *b_gamma2_blob_energy;   //!
   TBranch        *b_gamma2_blob_minsep;   //!
   TBranch        *b_gamma2_dEdx;   //!
   TBranch        *b_gamma2_dist_vtx;   //!
   TBranch        *b_gamma2_energy_ecal;   //!
   TBranch        *b_gamma2_energy_hcal;   //!
   TBranch        *b_gamma2_energy_scal_UV;   //!
   TBranch        *b_gamma2_energy_scal_X;   //!
   TBranch        *b_gamma2_energy_trkr;   //!
   TBranch        *b_gamma2_evis_ecal;   //!
   TBranch        *b_gamma2_evis_hcal;   //!
   TBranch        *b_gamma2_evis_scal_UV;   //!
   TBranch        *b_gamma2_evis_scal_X;   //!
   TBranch        *b_gamma2_evis_trkr;   //!
   TBranch        *b_gamma2_phi;   //!
   TBranch        *b_gamma2_px;   //!
   TBranch        *b_gamma2_py;   //!
   TBranch        *b_gamma2_pz;   //!
   TBranch        *b_gamma2_theta;   //!
   TBranch        *b_gamma2_time;   //!
   TBranch        *b_michelProng_begin_Z;   //!
   TBranch        *b_michelProng_distance;   //!
   TBranch        *b_michelProng_end_Z;   //!
   TBranch        *b_michelProng_energy;   //!
   TBranch        *b_michelProng_time_diff;   //!
   TBranch        *b_muon_E;   //!
   TBranch        *b_muon_E_shift;   //!
   TBranch        *b_muon_KE;   //!
   TBranch        *b_muon_P;   //!
   TBranch        *b_muon_muScore;   //!
   TBranch        *b_muon_phi;   //!
   TBranch        *b_muon_px;   //!
   TBranch        *b_muon_py;   //!
   TBranch        *b_muon_pz;   //!
   TBranch        *b_muon_qp;   //!
   TBranch        *b_muon_qpqpe;   //!
   TBranch        *b_muon_roadUpstreamEnergy;   //!
   TBranch        *b_muon_theta;   //!
   TBranch        *b_muon_thetaX;   //!
   TBranch        *b_muon_thetaY;   //!
   TBranch        *b_muon_theta_biasDown;   //!
   TBranch        *b_muon_theta_biasUp;   //!
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
   TBranch        *b_pi0_KE;   //!
   TBranch        *b_pi0_P;   //!
   TBranch        *b_pi0_cos_openingAngle;   //!
   TBranch        *b_pi0_invMass;   //!
   TBranch        *b_pi0_invMass_Old;   //!
   TBranch        *b_pi0_openingAngle;   //!
   TBranch        *b_pi0_phi;   //!
   TBranch        *b_pi0_px;   //!
   TBranch        *b_pi0_py;   //!
   TBranch        *b_pi0_pz;   //!
   TBranch        *b_pi0_theta;   //!
   TBranch        *b_pi0_thetaX;   //!
   TBranch        *b_pi0_thetaY;   //!
   TBranch        *b_preFilter_evis_ECAL;   //!
   TBranch        *b_preFilter_evis_HCAL;   //!
   TBranch        *b_preFilter_evis_NuclearTarget;   //!
   TBranch        *b_preFilter_evis_TotalExceptNuclearTarget;   //!
   TBranch        *b_preFilter_evis_Tracker;   //!
   TBranch        *b_preFilter_evis_nearvtx;   //!
   TBranch        *b_preFilter_evis_total;   //!
   TBranch        *b_preFilter_rejectedEnergy;   //!
   TBranch        *b_prim_vtx_smallest_opening_angle;   //!
   TBranch        *b_proton_E;   //!
   TBranch        *b_proton_KE;   //!
   TBranch        *b_proton_LLRScore;   //!
   TBranch        *b_proton_P;   //!
   TBranch        *b_proton_length;   //!
   TBranch        *b_proton_phi;   //!
   TBranch        *b_proton_pionScore;   //!
   TBranch        *b_proton_protonScore;   //!
   TBranch        *b_proton_px;   //!
   TBranch        *b_proton_py;   //!
   TBranch        *b_proton_pz;   //!
   TBranch        *b_proton_theta;   //!
   TBranch        *b_proton_thetaX;   //!
   TBranch        *b_proton_thetaY;   //!
   TBranch        *b_reco_eventID;   //!
   TBranch        *b_time;   //!
   TBranch        *b_vertex_blob_evis;   //!
   TBranch        *b_vtx_fit_chi2;   //!
   TBranch        *b_vtx_x;   //!
   TBranch        *b_vtx_y;   //!
   TBranch        *b_vtx_z;   //!
   TBranch        *b_all_protons_kinked;   //!
   TBranch        *b_all_protons_odMatch;   //!
   TBranch        *b_g1dedx_cluster_occupancy_sz;   //!
   TBranch        *b_g1dedx_cluster_occupancy;   //!
   TBranch        *b_g2dedx_cluster_occupancy_sz;   //!
   TBranch        *b_g2dedx_cluster_occupancy;   //!
   TBranch        *b_nTracks_Secondary_Vtx_sz;   //!
   TBranch        *b_nTracks_Secondary_Vtx;   //!
   TBranch        *b_all_protons_E;   //!
   TBranch        *b_all_protons_KE;   //!
   TBranch        *b_all_protons_LLRScore;   //!
   TBranch        *b_all_protons_P;   //!
   TBranch        *b_all_protons_chi2_ndf;   //!
   TBranch        *b_all_protons_endPointX;   //!
   TBranch        *b_all_protons_endPointY;   //!
   TBranch        *b_all_protons_endPointZ;   //!
   TBranch        *b_all_protons_length;   //!
   TBranch        *b_all_protons_p_calCorrection;   //!
   TBranch        *b_all_protons_p_dEdXTool;   //!
   TBranch        *b_all_protons_p_visEnergy;   //!
   TBranch        *b_all_protons_phi;   //!
   TBranch        *b_all_protons_pionScore;   //!
   TBranch        *b_all_protons_protonScore;   //!
   TBranch        *b_all_protons_px;   //!
   TBranch        *b_all_protons_py;   //!
   TBranch        *b_all_protons_pz;   //!
   TBranch        *b_all_protons_startPointX;   //!
   TBranch        *b_all_protons_startPointY;   //!
   TBranch        *b_all_protons_startPointZ;   //!
   TBranch        *b_all_protons_theta;   //!
   TBranch        *b_all_protons_thetaX;   //!
   TBranch        *b_all_protons_thetaY;   //!
   TBranch        *b_fit_vtx;   //!
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
   TBranch        *b_truth_isGamma1_conv_inside;   //!
   TBranch        *b_truth_isGamma2_conv_inside;   //!
   TBranch        *b_truth_isSignal;   //!
   TBranch        *b_truth_isFidVol;   //!
   TBranch        *b_truth_isNC;   //!
   TBranch        *b_truth_ReconstructEvent;   //!
   TBranch        *b_truth_isBckg_NoPi0;   //!
   TBranch        *b_truth_isBckg_SinglePi0;   //!
   TBranch        *b_truth_isBckg_MultiPi0;   //!
   TBranch        *b_truth_isBckg_NC;   //!
   TBranch        *b_truth_isBckg_AntiNeutrino;   //!
   TBranch        *b_truth_isBckg_QELike;   //!
   TBranch        *b_truth_isBckg_SingleChargedPion;   //!
   TBranch        *b_truth_isBckg_SingleChargedPion_ChargeExchanged;   //!
   TBranch        *b_truth_isBckg_DoublePionWithPi0;   //!
   TBranch        *b_truth_isBckg_DoublePionWithoutPi0;   //!
   TBranch        *b_truth_isBckg_MultiPionWithPi0;   //!
   TBranch        *b_truth_isBckg_MultiPionWithoutPi0;   //!
   TBranch        *b_truth_isBckg_Other;   //!
   TBranch        *b_truth_isBckg_withMichel;   //!
   TBranch        *b_truth_Bckg_nOther;   //!
   TBranch        *b_truth_Bckg_nPi0_Primary;   //!
   TBranch        *b_truth_Bckg_nPi0_Secondary;   //!
   TBranch        *b_truth_Bckg_nPi0_Total;   //!
   TBranch        *b_truth_Bckg_nPiCharged;   //!
   TBranch        *b_truth_Bckg_nPiCharged_ChargeExchanged;   //!
   TBranch        *b_truth_N_FSParticles;   //!
   TBranch        *b_truth_N_other;   //!
   TBranch        *b_truth_N_pi0;   //!
   TBranch        *b_truth_N_proton;   //!
   TBranch        *b_truth_N_trueMichelElectrons;   //!
   TBranch        *b_truth_OneShower_evis_most_pdg;   //!
   TBranch        *b_truth_Rejected_unused_evis_most_pdg;   //!
   TBranch        *b_truth_ThreeShower_s1_evis_most_pdg;   //!
   TBranch        *b_truth_ThreeShower_s2_evis_most_pdg;   //!
   TBranch        *b_truth_ThreeShower_s3_evis_most_pdg;   //!
   TBranch        *b_truth_blob1_evis_most_pdg;   //!
   TBranch        *b_truth_blob2_evis_most_pdg;   //!
   TBranch        *b_truth_dispersed_unused_evis_most_pdg;   //!
   TBranch        *b_truth_pi0_GrandMother;   //!
   TBranch        *b_truth_pi0_GrandMotherStatus;   //!
   TBranch        *b_truth_pi0_Mother;   //!
   TBranch        *b_truth_pi0_MotherStatus;   //!
   TBranch        *b_truth_pi0_status;   //!
   TBranch        *b_truth_target_material;   //!
   TBranch        *b_truth_vertex_module;   //!
   TBranch        *b_truth_vertex_plane;   //!
   TBranch        *b_truth_vertex_unused_evis_most_pdg;   //!
   TBranch        *b_truth_OneShower_evis_muon;   //!
   TBranch        *b_truth_OneShower_evis_neutron;   //!
   TBranch        *b_truth_OneShower_evis_piminus;   //!
   TBranch        *b_truth_OneShower_evis_piplus;   //!
   TBranch        *b_truth_OneShower_evis_pizero;   //!
   TBranch        *b_truth_OneShower_evis_proton;   //!
   TBranch        *b_truth_OneShower_evis_total_norm;   //!
   TBranch        *b_truth_OneShower_evis_total_truth;   //!
   TBranch        *b_truth_Rejected_unused_evis_total_norm;   //!
   TBranch        *b_truth_Rejected_unused_evis_total_truth;   //!
   TBranch        *b_truth_ThreeShower_s1_evis_muon;   //!
   TBranch        *b_truth_ThreeShower_s1_evis_neutron;   //!
   TBranch        *b_truth_ThreeShower_s1_evis_piminus;   //!
   TBranch        *b_truth_ThreeShower_s1_evis_piplus;   //!
   TBranch        *b_truth_ThreeShower_s1_evis_pizero;   //!
   TBranch        *b_truth_ThreeShower_s1_evis_proton;   //!
   TBranch        *b_truth_ThreeShower_s1_evis_total_norm;   //!
   TBranch        *b_truth_ThreeShower_s1_evis_total_truth;   //!
   TBranch        *b_truth_ThreeShower_s2_evis_muon;   //!
   TBranch        *b_truth_ThreeShower_s2_evis_neutron;   //!
   TBranch        *b_truth_ThreeShower_s2_evis_piminus;   //!
   TBranch        *b_truth_ThreeShower_s2_evis_piplus;   //!
   TBranch        *b_truth_ThreeShower_s2_evis_pizero;   //!
   TBranch        *b_truth_ThreeShower_s2_evis_proton;   //!
   TBranch        *b_truth_ThreeShower_s2_evis_total_norm;   //!
   TBranch        *b_truth_ThreeShower_s2_evis_total_truth;   //!
   TBranch        *b_truth_ThreeShower_s3_evis_muon;   //!
   TBranch        *b_truth_ThreeShower_s3_evis_neutron;   //!
   TBranch        *b_truth_ThreeShower_s3_evis_piminus;   //!
   TBranch        *b_truth_ThreeShower_s3_evis_piplus;   //!
   TBranch        *b_truth_ThreeShower_s3_evis_pizero;   //!
   TBranch        *b_truth_ThreeShower_s3_evis_proton;   //!
   TBranch        *b_truth_ThreeShower_s3_evis_total_norm;   //!
   TBranch        *b_truth_ThreeShower_s3_evis_total_truth;   //!
   TBranch        *b_truth_allClusters_evis_pizero;   //!
   TBranch        *b_truth_blob1_evis_muon;   //!
   TBranch        *b_truth_blob1_evis_neutron;   //!
   TBranch        *b_truth_blob1_evis_piminus;   //!
   TBranch        *b_truth_blob1_evis_piplus;   //!
   TBranch        *b_truth_blob1_evis_pizero;   //!
   TBranch        *b_truth_blob1_evis_proton;   //!
   TBranch        *b_truth_blob1_evis_total_norm;   //!
   TBranch        *b_truth_blob1_evis_total_truth;   //!
   TBranch        *b_truth_blob2_evis_muon;   //!
   TBranch        *b_truth_blob2_evis_neutron;   //!
   TBranch        *b_truth_blob2_evis_piminus;   //!
   TBranch        *b_truth_blob2_evis_piplus;   //!
   TBranch        *b_truth_blob2_evis_pizero;   //!
   TBranch        *b_truth_blob2_evis_proton;   //!
   TBranch        *b_truth_blob2_evis_total_norm;   //!
   TBranch        *b_truth_blob2_evis_total_truth;   //!
   TBranch        *b_truth_dispersed_unused_evis_gamma;   //!
   TBranch        *b_truth_dispersed_unused_evis_muon;   //!
   TBranch        *b_truth_dispersed_unused_evis_neutron;   //!
   TBranch        *b_truth_dispersed_unused_evis_piminus;   //!
   TBranch        *b_truth_dispersed_unused_evis_piplus;   //!
   TBranch        *b_truth_dispersed_unused_evis_pizero;   //!
   TBranch        *b_truth_dispersed_unused_evis_proton;   //!
   TBranch        *b_truth_dispersed_unused_evis_total_norm;   //!
   TBranch        *b_truth_dispersed_unused_evis_total_truth;   //!
   TBranch        *b_truth_ecal_unused_evis_muon;   //!
   TBranch        *b_truth_ecal_unused_evis_neutron;   //!
   TBranch        *b_truth_ecal_unused_evis_piminus;   //!
   TBranch        *b_truth_ecal_unused_evis_piplus;   //!
   TBranch        *b_truth_ecal_unused_evis_pizero;   //!
   TBranch        *b_truth_ecal_unused_evis_proton;   //!
   TBranch        *b_truth_ecal_unused_evis_total_norm;   //!
   TBranch        *b_truth_ecal_unused_evis_total_truth;   //!
   TBranch        *b_truth_eventID;   //!
   TBranch        *b_truth_hcal_unused_evis_muon;   //!
   TBranch        *b_truth_hcal_unused_evis_neutron;   //!
   TBranch        *b_truth_hcal_unused_evis_piminus;   //!
   TBranch        *b_truth_hcal_unused_evis_piplus;   //!
   TBranch        *b_truth_hcal_unused_evis_pizero;   //!
   TBranch        *b_truth_hcal_unused_evis_proton;   //!
   TBranch        *b_truth_hcal_unused_evis_total_norm;   //!
   TBranch        *b_truth_hcal_unused_evis_total_truth;   //!
   TBranch        *b_truth_michelElectron_E;   //!
   TBranch        *b_truth_michelElectron_P;   //!
   TBranch        *b_truth_michelMuon_P;   //!
   TBranch        *b_truth_michelMuon_end_dist_vtx;   //!
   TBranch        *b_truth_michelMuon_length;   //!
   TBranch        *b_truth_michelPion_P;   //!
   TBranch        *b_truth_michelPion_begin_dist_vtx;   //!
   TBranch        *b_truth_michelPion_length;   //!
   TBranch        *b_truth_muon_P;   //!
   TBranch        *b_truth_muon_theta;   //!
   TBranch        *b_truth_other_unused_evis_muon;   //!
   TBranch        *b_truth_other_unused_evis_neutron;   //!
   TBranch        *b_truth_other_unused_evis_piminus;   //!
   TBranch        *b_truth_other_unused_evis_piplus;   //!
   TBranch        *b_truth_other_unused_evis_pizero;   //!
   TBranch        *b_truth_other_unused_evis_proton;   //!
   TBranch        *b_truth_other_unused_evis_total_norm;   //!
   TBranch        *b_truth_other_unused_evis_total_truth;   //!
   TBranch        *b_truth_pi0_P;   //!
   TBranch        *b_truth_pi0_theta;   //!
   TBranch        *b_truth_proton_P;   //!
   TBranch        *b_truth_proton_theta;   //!
   TBranch        *b_truth_total_captured_evis_pizero;   //!
   TBranch        *b_truth_total_captured_evis_total_norm;   //!
   TBranch        *b_truth_total_captured_evis_total_truth;   //!
   TBranch        *b_truth_vertex_unused_evis_gamma;   //!
   TBranch        *b_truth_vertex_unused_evis_muon;   //!
   TBranch        *b_truth_vertex_unused_evis_neutron;   //!
   TBranch        *b_truth_vertex_unused_evis_piminus;   //!
   TBranch        *b_truth_vertex_unused_evis_piplus;   //!
   TBranch        *b_truth_vertex_unused_evis_pizero;   //!
   TBranch        *b_truth_vertex_unused_evis_proton;   //!
   TBranch        *b_truth_vertex_unused_evis_total_norm;   //!
   TBranch        *b_truth_vertex_unused_evis_total_truth;   //!
   TBranch        *b_truth_gamma1_4P;   //!
   TBranch        *b_truth_gamma1_final_pos;   //!
   TBranch        *b_truth_gamma1_init_pos;   //!
   TBranch        *b_truth_gamma2_4P;   //!
   TBranch        *b_truth_gamma2_final_pos;   //!
   TBranch        *b_truth_gamma2_init_pos;   //!
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
   TBranch        *b_truth_michelMuon_endPoint;   //!
   TBranch        *b_truth_muon_4P;   //!
   TBranch        *b_truth_pi0_4P;   //!
   TBranch        *b_truth_proton_4P;   //!
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
   TBranch        *b_CCProtonPi0_ntrajMuonProng;   //!
   TBranch        *b_CCProtonPi0_r_minos_trk_vtx_plane;   //!
   TBranch        *b_CCProtonPi0_t_minos_trk_numFSMuons;   //!
   TBranch        *b_CCProtonPi0_t_minos_trk_primFSLeptonPDG;   //!
   TBranch        *b_CCProtonPi0_trajMuonProngPDG;   //!
   TBranch        *b_CCProtonPi0_trajMuonProngPrimary;   //!
   TBranch        *b_CCProtonPi0_Extra_Energy_Total;   //!
   TBranch        *b_CCProtonPi0_QSq;   //!
   TBranch        *b_CCProtonPi0_QSq_1Track_Alt;   //!
   TBranch        *b_CCProtonPi0_WSq;   //!
   TBranch        *b_CCProtonPi0_WSq_1Track_Alt;   //!
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
   TBranch        *b_CCProtonPi0_neutrino_E;   //!
   TBranch        *b_CCProtonPi0_neutrino_E_1Track_Alt;   //!
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
   TBranch        *b_CCProtonPi0_vertex_energy;   //!
   TBranch        *b_CCProtonPi0_isProtonInsideOD;   //!
   TBranch        *b_CCProtonPi0_ntrajProtonProng;   //!
   TBranch        *b_CCProtonPi0_trajProtonProngPDG;   //!
   TBranch        *b_CCProtonPi0_trajProtonProngPrimary;   //!
   TBranch        *b_CCProtonPi0_endProtonTrajMomentum;   //!
   TBranch        *b_CCProtonPi0_endProtonTrajXPosition;   //!
   TBranch        *b_CCProtonPi0_endProtonTrajYPosition;   //!
   TBranch        *b_CCProtonPi0_endProtonTrajZPosition;   //!
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
   TBranch        *b_mc_ppfx1_cvweight;   //!
   TBranch        *b_mc_gen1_cvweight_total;   //!
   TBranch        *b_gen1_wgt;   //!
   TBranch        *b_mc_gen1_cvweight_totalFlux;   //!
   TBranch        *b_mc_gen1_cvweight_NA49;   //!
   TBranch        *b_mc_wgt_Flux_BeamFocus_sz;   //!
   TBranch        *b_mc_wgt_Flux_BeamFocus;   //!
   TBranch        *b_mc_wgt_gen1_Flux_Tertiary_sz;   //!
   TBranch        *b_mc_wgt_gen1_Flux_Tertiary;   //!
   TBranch        *b_mc_wgt_gen1_Flux_NA49_sz;   //!
   TBranch        *b_mc_wgt_gen1_Flux_NA49;   //!
   TBranch        *b_mc_wgt_Norm_sz;   //!
   TBranch        *b_mc_wgt_Norm;   //!
   TBranch        *b_mc_wgt_ppfx1_Total_sz;   //!
   TBranch        *b_mc_wgt_ppfx1_Total;   //!
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
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0001-0002-0003-0004_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0005-0006-0007-0008_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0009-0010-0011-0012-0013_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0014-0015-0016-0017-0018_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0019-0020-0021-0023-0024_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0025-0026-0027-0028-0029_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0030-0031-0032-0033-0034_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0035-0036-0037-0038-0039_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0040-0041-0042-0043-0044_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0045-0046-0047-0048-0049_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0050-0051-0052-0053-0054_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0055-0056-0057-0058-0059_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0060-0061-0062-0063-0064_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0065-0066-0067-0068-0069_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0070-0071-0072-0073-0074_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0075-0076-0077-0078-0079_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0080-0081-0082-0083-0084_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0085-0086-0087-0088-0089_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0090-0091-0092-0093-0094_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0095-0096-0097-0098-0099_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0100-0101-0102-0103-0104_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0105-0106-0107-0108-0109_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0110-0111-0112-0113-0114_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0115-0116-0117-0118-0119_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0120-0121-0122-0123-0124_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0125-0126-0127-0128-0129_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0130-0131-0132-0133-0134_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0135-0136-0137-0138-0139_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0140-0141-0142-0143-0144_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0145-0147-0148-0149-0150_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0151-0152-0153-0154-0155_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0156-0157-0158-0159-0160_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0161-0162-0163-0164-0165_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0166-0167-0168-0169-0170_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0171-0172-0173-0174-0175_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0176-0177-0178-0179-0180_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0181-0182-0183-0184-0185_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0186-0187-0188-0189-0190_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0191-0192-0193-0194-0195_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/55/SIM_minerva_00010255_Subruns_0196-0197-0198-0199-0200_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0001-0002-0003-0004_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0005-0006-0007-0008_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0009-0010-0011-0012_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0013-0014-0015-0016-0017_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0018-0019-0020-0021-0022_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0023-0024-0025-0026-0027_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0028-0029-0030-0031-0032_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0033-0034-0035-0036-0037_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0038-0039-0040-0041-0042_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0043-0044-0045-0046-0047_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0048-0049-0050-0051-0052_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0053-0054-0055-0056-0057_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0058-0059-0060-0061-0062_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0063-0064-0065-0066-0067_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0068-0069-0070-0071-0072_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0073-0074-0075-0076-0077_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0078-0079-0080-0081-0082_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0083-0084-0085-0086-0087_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0088-0089-0090-0091-0092_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0093-0094-0095-0096-0097_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0098-0099-0100-0101-0102_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0103-0104-0105-0106-0107_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0108-0109-0110-0111-0112_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0113-0114-0115-0116-0117_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0118-0119-0120-0121-0122_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0123-0124-0125-0126-0127_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0128-0129-0130-0131-0132_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0133-0134-0135-0136-0137_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0138-0139-0140-0141-0142_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0143-0145-0146-0147-0148_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0149-0150-0151-0152-0153_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0154-0155-0156-0157-0158_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0159-0160-0161-0162-0163_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0164-0165-0166-0167-0168_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0169-0170-0172-0173-0174_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0175-0176-0178-0179-0180_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0181-0182-0183-0184-0185_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0186-0187-0188-0189-0190_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0191-0192-0193-0194-0195_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/56/SIM_minerva_00010256_Subruns_0196-0197-0198-0199-0200_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0001-0002-0003-0004_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0005-0006-0007-0008-0009_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0010-0011-0012-0014-0016_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0017-0018-0019-0020-0021_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0022-0023-0024-0025-0026_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0027-0028-0029-0030-0031_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0032-0033-0034-0035-0036_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0042-0043-0044-0045-0046_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0047-0048-0049-0050-0051_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0052-0053-0054-0055-0057_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0058-0059-0060-0061-0064_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0065-0066-0067-0068-0070_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0071-0072-0073-0074-0075_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0076-0077-0078-0079-0080_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0081-0082-0083-0084-0085_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0086-0087-0088-0089-0090_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0091-0092-0093-0094-0095_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0096-0097-0098-0099-0100_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0101-0102-0103-0104-0105_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0106-0107-0108-0109-0110_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0111-0112-0113-0114-0115_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0121-0122-0123-0124-0125_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0126-0127-0128-0129-0130_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0131-0132-0133-0134-0135_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0136-0137-0138-0139-0140_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0141-0142-0143-0144-0145_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0146-0147-0148-0149-0150_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0151-0152-0153-0154-0155_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0156-0157-0158-0159-0160_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0161-0162-0163-0164-0165_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0166-0167-0168-0169-0170_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0171-0172-0173-0174-0175_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0176-0177-0178-0179-0180_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0181-0182-0183-0184-0185_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0186-0187-0188-0189-0190_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0191-0192-0193-0194-0195_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/57/SIM_minerva_00010257_Subruns_0196-0197-0198-0199-0200_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0001-0002-0003-0004_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0009-0010-0011-0012_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0013-0014-0015-0016_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0017-0018-0019-0020-0021_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0022-0023-0024-0025-0026_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0027-0028-0029-0030-0031_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0032-0033-0034-0035-0036_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0037-0038-0039-0040-0041_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0042-0043-0044-0045-0046_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0047-0048-0049-0050-0051_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0052-0053-0054-0055-0056_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0057-0058-0059-0060-0061_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0062-0063-0064-0065-0066_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0067-0068-0069-0070-0071_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0072-0073-0074-0075-0076_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0077-0078-0079-0080-0081_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0082-0084-0085-0086-0087_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0088-0089-0090-0091-0092_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0093-0094-0095-0096-0097_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0098-0099-0100-0101-0102_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0103-0104-0105-0106-0107_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0114-0115-0116-0117-0118_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0119-0120-0121-0122-0123_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0124-0125-0126-0127-0128_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0129-0130-0131-0132-0134_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0135-0136-0137-0138-0139_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0140-0141-0142-0143-0144_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0150-0151-0152-0153-0154_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0155-0156-0157-0158-0159_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0160-0161-0162-0163-0165_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0166-0167-0168-0169-0170_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0171-0172-0173-0174-0175_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0176-0177-0178-0179-0180_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0181-0182-0183-0184-0185_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0186-0187-0188-0189-0190_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0191-0192-0193-0194-0195_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/58/SIM_minerva_00010258_Subruns_0196-0197-0198-0199-0200_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0001-0002-0003-0004_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0005-0006-0007-0008_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0009-0010-0011-0012-0013_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0014-0015-0016-0017-0018_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0019-0020-0021-0022-0023_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0024-0025-0026-0027-0028_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0034-0035-0036-0037-0038_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0039-0040-0041-0042-0043_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0049-0050-0051-0052-0053_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0054-0055-0056-0057-0058_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0059-0060-0061-0062-0063_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0064-0065-0066-0068-0069_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0070-0071-0072-0073-0074_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0075-0076-0077-0078-0079_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0080-0081-0082-0083-0084_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0085-0086-0087-0088-0089_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0090-0091-0092-0093-0094_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0095-0096-0097-0098-0099_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0100-0101-0103-0104-0105_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0106-0107-0108-0109-0110_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0111-0112-0113-0114-0115_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0116-0117-0118-0119-0120_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0121-0122-0123-0124-0125_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0126-0127-0128-0129-0130_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0131-0132-0133-0134-0135_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0136-0137-0138-0139-0140_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0141-0142-0143-0144-0145_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0146-0147-0148-0149-0150_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0151-0152-0153-0154-0155_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0156-0157-0158-0159-0160_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0161-0162-0163-0164-0165_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0166-0167-0168-0169-0170_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0171-0172-0173-0174-0175_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0176-0177-0178-0179-0180_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0181-0182-0183-0184-0185_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0186-0187-0188-0189-0190_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0191-0192-0193-0194-0195_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/minerva9/grid/central_value/minerva/ana/v10r8p7/00/01/02/59/SIM_minerva_00010259_Subruns_0196-0197-0198-0199-0200_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
      chain->Add("/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/v2_55/run/grid/central_value/minerva/ana/v10r8p7/00/01/02/00/SIM_minerva_00010200_Subruns_0001-0002-0003-0004_CCProtonPi0_Ana_Tuple_v10r8p7-oaltinok.root/CCProtonPi0");
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
   fChain->SetBranchAddress("is_blobs_recovered", &is_blobs_recovered, &b_is_blobs_recovered);
   fChain->SetBranchAddress("is_blobs_recovered_direction", &is_blobs_recovered_direction, &b_is_blobs_recovered_direction);
   fChain->SetBranchAddress("is_blobs_recovered_invMass", &is_blobs_recovered_invMass, &b_is_blobs_recovered_invMass);
   fChain->SetBranchAddress("is_blobs_recovered_small_angle", &is_blobs_recovered_small_angle, &b_is_blobs_recovered_small_angle);
   fChain->SetBranchAddress("is_blobs_recovered_search_view_U", &is_blobs_recovered_search_view_U, &b_is_blobs_recovered_search_view_U);
   fChain->SetBranchAddress("is_blobs_recovered_search_view_V", &is_blobs_recovered_search_view_V, &b_is_blobs_recovered_search_view_V);
   fChain->SetBranchAddress("Cut_BlobDirectionBad", &Cut_BlobDirectionBad, &b_Cut_BlobDirectionBad);
   fChain->SetBranchAddress("Cut_ConeBlobs", &Cut_ConeBlobs, &b_Cut_ConeBlobs);
   fChain->SetBranchAddress("Cut_EndPoint_Michel_Exist", &Cut_EndPoint_Michel_Exist, &b_Cut_EndPoint_Michel_Exist);
   fChain->SetBranchAddress("Cut_Muon_Charge", &Cut_Muon_Charge, &b_Cut_Muon_Charge);
   fChain->SetBranchAddress("Cut_Muon_None", &Cut_Muon_None, &b_Cut_Muon_None);
   fChain->SetBranchAddress("Cut_Particle_None", &Cut_Particle_None, &b_Cut_Particle_None);
   fChain->SetBranchAddress("Cut_Pi0_Bad", &Cut_Pi0_Bad, &b_Cut_Pi0_Bad);
   fChain->SetBranchAddress("Cut_PreFilter_Pi0", &Cut_PreFilter_Pi0, &b_Cut_PreFilter_Pi0);
   fChain->SetBranchAddress("Cut_Proton_Bad", &Cut_Proton_Bad, &b_Cut_Proton_Bad);
   fChain->SetBranchAddress("Cut_Proton_None", &Cut_Proton_None, &b_Cut_Proton_None);
   fChain->SetBranchAddress("Cut_Vertex_Michel_Exist", &Cut_Vertex_Michel_Exist, &b_Cut_Vertex_Michel_Exist);
   fChain->SetBranchAddress("Cut_Vertex_None", &Cut_Vertex_None, &b_Cut_Vertex_None);
   fChain->SetBranchAddress("Cut_Vertex_Not_Fiducial", &Cut_Vertex_Not_Fiducial, &b_Cut_Vertex_Not_Fiducial);
   fChain->SetBranchAddress("Cut_Vertex_Not_Reconstructable", &Cut_Vertex_Not_Reconstructable, &b_Cut_Vertex_Not_Reconstructable);
   fChain->SetBranchAddress("Cut_secEndPoint_Michel_Exist", &Cut_secEndPoint_Michel_Exist, &b_Cut_secEndPoint_Michel_Exist);
   fChain->SetBranchAddress("OneShower_nClusters", &OneShower_nClusters, &b_OneShower_nClusters);
   fChain->SetBranchAddress("ThreeShower_s1_nClusters", &ThreeShower_s1_nClusters, &b_ThreeShower_s1_nClusters);
   fChain->SetBranchAddress("ThreeShower_s2_nClusters", &ThreeShower_s2_nClusters, &b_ThreeShower_s2_nClusters);
   fChain->SetBranchAddress("ThreeShower_s3_nClusters", &ThreeShower_s3_nClusters, &b_ThreeShower_s3_nClusters);
   fChain->SetBranchAddress("anglescan_ncand", &anglescan_ncand, &b_anglescan_ncand);
   fChain->SetBranchAddress("anglescan_ncandx", &anglescan_ncandx, &b_anglescan_ncandx);
   fChain->SetBranchAddress("anglescan_nfoundBlobs", &anglescan_nfoundBlobs, &b_anglescan_nfoundBlobs);
   fChain->SetBranchAddress("broken_track_most_us_plane", &broken_track_most_us_plane, &b_broken_track_most_us_plane);
   fChain->SetBranchAddress("g1blob_1ParFit_ndof", &g1blob_1ParFit_ndof, &b_g1blob_1ParFit_ndof);
   fChain->SetBranchAddress("g1dedx_doublet", &g1dedx_doublet, &b_g1dedx_doublet);
   fChain->SetBranchAddress("g1dedx_empty_plane", &g1dedx_empty_plane, &b_g1dedx_empty_plane);
   fChain->SetBranchAddress("g1dedx_nplane", &g1dedx_nplane, &b_g1dedx_nplane);
   fChain->SetBranchAddress("g2blob_1ParFit_ndof", &g2blob_1ParFit_ndof, &b_g2blob_1ParFit_ndof);
   fChain->SetBranchAddress("g2dedx_doublet", &g2dedx_doublet, &b_g2dedx_doublet);
   fChain->SetBranchAddress("g2dedx_empty_plane", &g2dedx_empty_plane, &b_g2dedx_empty_plane);
   fChain->SetBranchAddress("g2dedx_nplane", &g2dedx_nplane, &b_g2dedx_nplane);
   fChain->SetBranchAddress("gamma1_blob_nclusters", &gamma1_blob_nclusters, &b_gamma1_blob_nclusters);
   fChain->SetBranchAddress("gamma1_blob_ndigits", &gamma1_blob_ndigits, &b_gamma1_blob_ndigits);
   fChain->SetBranchAddress("gamma2_blob_nclusters", &gamma2_blob_nclusters, &b_gamma2_blob_nclusters);
   fChain->SetBranchAddress("gamma2_blob_ndigits", &gamma2_blob_ndigits, &b_gamma2_blob_ndigits);
   fChain->SetBranchAddress("muon_N_minosTracks", &muon_N_minosTracks, &b_muon_N_minosTracks);
   fChain->SetBranchAddress("muon_charge", &muon_charge, &b_muon_charge);
   fChain->SetBranchAddress("muon_hasMinosMatchStub", &muon_hasMinosMatchStub, &b_muon_hasMinosMatchStub);
   fChain->SetBranchAddress("muon_hasMinosMatchTrack", &muon_hasMinosMatchTrack, &b_muon_hasMinosMatchTrack);
   fChain->SetBranchAddress("muon_minervaTrack_types", &muon_minervaTrack_types, &b_muon_minervaTrack_types);
   fChain->SetBranchAddress("muon_minosTrackQuality", &muon_minosTrackQuality, &b_muon_minosTrackQuality);
   fChain->SetBranchAddress("muon_roadUpstreamPlanes", &muon_roadUpstreamPlanes, &b_muon_roadUpstreamPlanes);
   fChain->SetBranchAddress("nProtonCandidates", &nProtonCandidates, &b_nProtonCandidates);
   fChain->SetBranchAddress("nTracks", &nTracks, &b_nTracks);
   fChain->SetBranchAddress("nTracks_Close", &nTracks_Close, &b_nTracks_Close);
   fChain->SetBranchAddress("nTracks_Discarded", &nTracks_Discarded, &b_nTracks_Discarded);
   fChain->SetBranchAddress("nTracks_Far", &nTracks_Far, &b_nTracks_Far);
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
   fChain->SetBranchAddress("proton_kinked", &proton_kinked, &b_proton_kinked);
   fChain->SetBranchAddress("proton_leadingIndice", &proton_leadingIndice, &b_proton_leadingIndice);
   fChain->SetBranchAddress("vtx_fit_converged", &vtx_fit_converged, &b_vtx_fit_converged);
   fChain->SetBranchAddress("vtx_module", &vtx_module, &b_vtx_module);
   fChain->SetBranchAddress("vtx_plane", &vtx_plane, &b_vtx_plane);
   fChain->SetBranchAddress("vtx_primary_index", &vtx_primary_index, &b_vtx_primary_index);
   fChain->SetBranchAddress("vtx_primary_multiplicity", &vtx_primary_multiplicity, &b_vtx_primary_multiplicity);
   fChain->SetBranchAddress("vtx_secondary_count", &vtx_secondary_count, &b_vtx_secondary_count);
   fChain->SetBranchAddress("vtx_total_count", &vtx_total_count, &b_vtx_total_count);
   fChain->SetBranchAddress("ConeBlobs_usable_evis_Tracker", &ConeBlobs_usable_evis_Tracker, &b_ConeBlobs_usable_evis_Tracker);
   fChain->SetBranchAddress("Coneblobs_usable_evis_ECAL", &Coneblobs_usable_evis_ECAL, &b_Coneblobs_usable_evis_ECAL);
   fChain->SetBranchAddress("Coneblobs_usable_evis_HCAL", &Coneblobs_usable_evis_HCAL, &b_Coneblobs_usable_evis_HCAL);
   fChain->SetBranchAddress("Extra_Energy_Dispersed", &Extra_Energy_Dispersed, &b_Extra_Energy_Dispersed);
   fChain->SetBranchAddress("Extra_Energy_Muon", &Extra_Energy_Muon, &b_Extra_Energy_Muon);
   fChain->SetBranchAddress("Extra_Energy_Rejected", &Extra_Energy_Rejected, &b_Extra_Energy_Rejected);
   fChain->SetBranchAddress("Extra_Evis_Leftover", &Extra_Evis_Leftover, &b_Extra_Evis_Leftover);
   fChain->SetBranchAddress("OneShower_dist_vtx", &OneShower_dist_vtx, &b_OneShower_dist_vtx);
   fChain->SetBranchAddress("OneShower_energy", &OneShower_energy, &b_OneShower_energy);
   fChain->SetBranchAddress("OneShower_theta", &OneShower_theta, &b_OneShower_theta);
   fChain->SetBranchAddress("ThreeShower_s1_dist_vtx", &ThreeShower_s1_dist_vtx, &b_ThreeShower_s1_dist_vtx);
   fChain->SetBranchAddress("ThreeShower_s1_energy", &ThreeShower_s1_energy, &b_ThreeShower_s1_energy);
   fChain->SetBranchAddress("ThreeShower_s1_theta", &ThreeShower_s1_theta, &b_ThreeShower_s1_theta);
   fChain->SetBranchAddress("ThreeShower_s2_dist_vtx", &ThreeShower_s2_dist_vtx, &b_ThreeShower_s2_dist_vtx);
   fChain->SetBranchAddress("ThreeShower_s2_energy", &ThreeShower_s2_energy, &b_ThreeShower_s2_energy);
   fChain->SetBranchAddress("ThreeShower_s2_theta", &ThreeShower_s2_theta, &b_ThreeShower_s2_theta);
   fChain->SetBranchAddress("ThreeShower_s3_dist_vtx", &ThreeShower_s3_dist_vtx, &b_ThreeShower_s3_dist_vtx);
   fChain->SetBranchAddress("ThreeShower_s3_energy", &ThreeShower_s3_energy, &b_ThreeShower_s3_energy);
   fChain->SetBranchAddress("ThreeShower_s3_theta", &ThreeShower_s3_theta, &b_ThreeShower_s3_theta);
   fChain->SetBranchAddress("energy_from_mc", &energy_from_mc, &b_energy_from_mc);
   fChain->SetBranchAddress("energy_from_mc_fraction", &energy_from_mc_fraction, &b_energy_from_mc_fraction);
   fChain->SetBranchAddress("energy_from_mc_fraction_of_highest", &energy_from_mc_fraction_of_highest, &b_energy_from_mc_fraction_of_highest);
   fChain->SetBranchAddress("g1blob_1ParFit_fval", &g1blob_1ParFit_fval, &b_g1blob_1ParFit_fval);
   fChain->SetBranchAddress("g1blob_2ParFit_vtx_distance", &g1blob_2ParFit_vtx_distance, &b_g1blob_2ParFit_vtx_distance);
   fChain->SetBranchAddress("g1dedx", &g1dedx, &b_g1dedx);
   fChain->SetBranchAddress("g1dedx1", &g1dedx1, &b_g1dedx1);
   fChain->SetBranchAddress("g1dedx_total", &g1dedx_total, &b_g1dedx_total);
   fChain->SetBranchAddress("g1dedx_total1", &g1dedx_total1, &b_g1dedx_total1);
   fChain->SetBranchAddress("g2blob_1ParFit_fval", &g2blob_1ParFit_fval, &b_g2blob_1ParFit_fval);
   fChain->SetBranchAddress("g2blob_2ParFit_vtx_distance", &g2blob_2ParFit_vtx_distance, &b_g2blob_2ParFit_vtx_distance);
   fChain->SetBranchAddress("g2dedx", &g2dedx, &b_g2dedx);
   fChain->SetBranchAddress("g2dedx1", &g2dedx1, &b_g2dedx1);
   fChain->SetBranchAddress("g2dedx_total", &g2dedx_total, &b_g2dedx_total);
   fChain->SetBranchAddress("g2dedx_total1", &g2dedx_total1, &b_g2dedx_total1);
   fChain->SetBranchAddress("gamma1_E", &gamma1_E, &b_gamma1_E);
   fChain->SetBranchAddress("gamma1_E_Old", &gamma1_E_Old, &b_gamma1_E_Old);
   fChain->SetBranchAddress("gamma1_P", &gamma1_P, &b_gamma1_P);
   fChain->SetBranchAddress("gamma1_blob_energy", &gamma1_blob_energy, &b_gamma1_blob_energy);
   fChain->SetBranchAddress("gamma1_blob_minsep", &gamma1_blob_minsep, &b_gamma1_blob_minsep);
   fChain->SetBranchAddress("gamma1_dEdx", &gamma1_dEdx, &b_gamma1_dEdx);
   fChain->SetBranchAddress("gamma1_dist_vtx", &gamma1_dist_vtx, &b_gamma1_dist_vtx);
   fChain->SetBranchAddress("gamma1_energy_ecal", &gamma1_energy_ecal, &b_gamma1_energy_ecal);
   fChain->SetBranchAddress("gamma1_energy_hcal", &gamma1_energy_hcal, &b_gamma1_energy_hcal);
   fChain->SetBranchAddress("gamma1_energy_scal_UV", &gamma1_energy_scal_UV, &b_gamma1_energy_scal_UV);
   fChain->SetBranchAddress("gamma1_energy_scal_X", &gamma1_energy_scal_X, &b_gamma1_energy_scal_X);
   fChain->SetBranchAddress("gamma1_energy_trkr", &gamma1_energy_trkr, &b_gamma1_energy_trkr);
   fChain->SetBranchAddress("gamma1_evis_ecal", &gamma1_evis_ecal, &b_gamma1_evis_ecal);
   fChain->SetBranchAddress("gamma1_evis_hcal", &gamma1_evis_hcal, &b_gamma1_evis_hcal);
   fChain->SetBranchAddress("gamma1_evis_scal_UV", &gamma1_evis_scal_UV, &b_gamma1_evis_scal_UV);
   fChain->SetBranchAddress("gamma1_evis_scal_X", &gamma1_evis_scal_X, &b_gamma1_evis_scal_X);
   fChain->SetBranchAddress("gamma1_evis_trkr", &gamma1_evis_trkr, &b_gamma1_evis_trkr);
   fChain->SetBranchAddress("gamma1_phi", &gamma1_phi, &b_gamma1_phi);
   fChain->SetBranchAddress("gamma1_px", &gamma1_px, &b_gamma1_px);
   fChain->SetBranchAddress("gamma1_py", &gamma1_py, &b_gamma1_py);
   fChain->SetBranchAddress("gamma1_pz", &gamma1_pz, &b_gamma1_pz);
   fChain->SetBranchAddress("gamma1_theta", &gamma1_theta, &b_gamma1_theta);
   fChain->SetBranchAddress("gamma1_time", &gamma1_time, &b_gamma1_time);
   fChain->SetBranchAddress("gamma2_E", &gamma2_E, &b_gamma2_E);
   fChain->SetBranchAddress("gamma2_E_Old", &gamma2_E_Old, &b_gamma2_E_Old);
   fChain->SetBranchAddress("gamma2_P", &gamma2_P, &b_gamma2_P);
   fChain->SetBranchAddress("gamma2_blob_energy", &gamma2_blob_energy, &b_gamma2_blob_energy);
   fChain->SetBranchAddress("gamma2_blob_minsep", &gamma2_blob_minsep, &b_gamma2_blob_minsep);
   fChain->SetBranchAddress("gamma2_dEdx", &gamma2_dEdx, &b_gamma2_dEdx);
   fChain->SetBranchAddress("gamma2_dist_vtx", &gamma2_dist_vtx, &b_gamma2_dist_vtx);
   fChain->SetBranchAddress("gamma2_energy_ecal", &gamma2_energy_ecal, &b_gamma2_energy_ecal);
   fChain->SetBranchAddress("gamma2_energy_hcal", &gamma2_energy_hcal, &b_gamma2_energy_hcal);
   fChain->SetBranchAddress("gamma2_energy_scal_UV", &gamma2_energy_scal_UV, &b_gamma2_energy_scal_UV);
   fChain->SetBranchAddress("gamma2_energy_scal_X", &gamma2_energy_scal_X, &b_gamma2_energy_scal_X);
   fChain->SetBranchAddress("gamma2_energy_trkr", &gamma2_energy_trkr, &b_gamma2_energy_trkr);
   fChain->SetBranchAddress("gamma2_evis_ecal", &gamma2_evis_ecal, &b_gamma2_evis_ecal);
   fChain->SetBranchAddress("gamma2_evis_hcal", &gamma2_evis_hcal, &b_gamma2_evis_hcal);
   fChain->SetBranchAddress("gamma2_evis_scal_UV", &gamma2_evis_scal_UV, &b_gamma2_evis_scal_UV);
   fChain->SetBranchAddress("gamma2_evis_scal_X", &gamma2_evis_scal_X, &b_gamma2_evis_scal_X);
   fChain->SetBranchAddress("gamma2_evis_trkr", &gamma2_evis_trkr, &b_gamma2_evis_trkr);
   fChain->SetBranchAddress("gamma2_phi", &gamma2_phi, &b_gamma2_phi);
   fChain->SetBranchAddress("gamma2_px", &gamma2_px, &b_gamma2_px);
   fChain->SetBranchAddress("gamma2_py", &gamma2_py, &b_gamma2_py);
   fChain->SetBranchAddress("gamma2_pz", &gamma2_pz, &b_gamma2_pz);
   fChain->SetBranchAddress("gamma2_theta", &gamma2_theta, &b_gamma2_theta);
   fChain->SetBranchAddress("gamma2_time", &gamma2_time, &b_gamma2_time);
   fChain->SetBranchAddress("michelProng_begin_Z", &michelProng_begin_Z, &b_michelProng_begin_Z);
   fChain->SetBranchAddress("michelProng_distance", &michelProng_distance, &b_michelProng_distance);
   fChain->SetBranchAddress("michelProng_end_Z", &michelProng_end_Z, &b_michelProng_end_Z);
   fChain->SetBranchAddress("michelProng_energy", &michelProng_energy, &b_michelProng_energy);
   fChain->SetBranchAddress("michelProng_time_diff", &michelProng_time_diff, &b_michelProng_time_diff);
   fChain->SetBranchAddress("muon_E", &muon_E, &b_muon_E);
   fChain->SetBranchAddress("muon_E_shift", &muon_E_shift, &b_muon_E_shift);
   fChain->SetBranchAddress("muon_KE", &muon_KE, &b_muon_KE);
   fChain->SetBranchAddress("muon_P", &muon_P, &b_muon_P);
   fChain->SetBranchAddress("muon_muScore", &muon_muScore, &b_muon_muScore);
   fChain->SetBranchAddress("muon_phi", &muon_phi, &b_muon_phi);
   fChain->SetBranchAddress("muon_px", &muon_px, &b_muon_px);
   fChain->SetBranchAddress("muon_py", &muon_py, &b_muon_py);
   fChain->SetBranchAddress("muon_pz", &muon_pz, &b_muon_pz);
   fChain->SetBranchAddress("muon_qp", &muon_qp, &b_muon_qp);
   fChain->SetBranchAddress("muon_qpqpe", &muon_qpqpe, &b_muon_qpqpe);
   fChain->SetBranchAddress("muon_roadUpstreamEnergy", &muon_roadUpstreamEnergy, &b_muon_roadUpstreamEnergy);
   fChain->SetBranchAddress("muon_theta", &muon_theta, &b_muon_theta);
   fChain->SetBranchAddress("muon_thetaX", &muon_thetaX, &b_muon_thetaX);
   fChain->SetBranchAddress("muon_thetaY", &muon_thetaY, &b_muon_thetaY);
   fChain->SetBranchAddress("muon_theta_biasDown", &muon_theta_biasDown, &b_muon_theta_biasDown);
   fChain->SetBranchAddress("muon_theta_biasUp", &muon_theta_biasUp, &b_muon_theta_biasUp);
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
   fChain->SetBranchAddress("pi0_KE", &pi0_KE, &b_pi0_KE);
   fChain->SetBranchAddress("pi0_P", &pi0_P, &b_pi0_P);
   fChain->SetBranchAddress("pi0_cos_openingAngle", &pi0_cos_openingAngle, &b_pi0_cos_openingAngle);
   fChain->SetBranchAddress("pi0_invMass", &pi0_invMass, &b_pi0_invMass);
   fChain->SetBranchAddress("pi0_invMass_Old", &pi0_invMass_Old, &b_pi0_invMass_Old);
   fChain->SetBranchAddress("pi0_openingAngle", &pi0_openingAngle, &b_pi0_openingAngle);
   fChain->SetBranchAddress("pi0_phi", &pi0_phi, &b_pi0_phi);
   fChain->SetBranchAddress("pi0_px", &pi0_px, &b_pi0_px);
   fChain->SetBranchAddress("pi0_py", &pi0_py, &b_pi0_py);
   fChain->SetBranchAddress("pi0_pz", &pi0_pz, &b_pi0_pz);
   fChain->SetBranchAddress("pi0_theta", &pi0_theta, &b_pi0_theta);
   fChain->SetBranchAddress("pi0_thetaX", &pi0_thetaX, &b_pi0_thetaX);
   fChain->SetBranchAddress("pi0_thetaY", &pi0_thetaY, &b_pi0_thetaY);
   fChain->SetBranchAddress("preFilter_evis_ECAL", &preFilter_evis_ECAL, &b_preFilter_evis_ECAL);
   fChain->SetBranchAddress("preFilter_evis_HCAL", &preFilter_evis_HCAL, &b_preFilter_evis_HCAL);
   fChain->SetBranchAddress("preFilter_evis_NuclearTarget", &preFilter_evis_NuclearTarget, &b_preFilter_evis_NuclearTarget);
   fChain->SetBranchAddress("preFilter_evis_TotalExceptNuclearTarget", &preFilter_evis_TotalExceptNuclearTarget, &b_preFilter_evis_TotalExceptNuclearTarget);
   fChain->SetBranchAddress("preFilter_evis_Tracker", &preFilter_evis_Tracker, &b_preFilter_evis_Tracker);
   fChain->SetBranchAddress("preFilter_evis_nearvtx", &preFilter_evis_nearvtx, &b_preFilter_evis_nearvtx);
   fChain->SetBranchAddress("preFilter_evis_total", &preFilter_evis_total, &b_preFilter_evis_total);
   fChain->SetBranchAddress("preFilter_rejectedEnergy", &preFilter_rejectedEnergy, &b_preFilter_rejectedEnergy);
   fChain->SetBranchAddress("prim_vtx_smallest_opening_angle", &prim_vtx_smallest_opening_angle, &b_prim_vtx_smallest_opening_angle);
   fChain->SetBranchAddress("proton_E", &proton_E, &b_proton_E);
   fChain->SetBranchAddress("proton_KE", &proton_KE, &b_proton_KE);
   fChain->SetBranchAddress("proton_LLRScore", &proton_LLRScore, &b_proton_LLRScore);
   fChain->SetBranchAddress("proton_P", &proton_P, &b_proton_P);
   fChain->SetBranchAddress("proton_length", &proton_length, &b_proton_length);
   fChain->SetBranchAddress("proton_phi", &proton_phi, &b_proton_phi);
   fChain->SetBranchAddress("proton_pionScore", &proton_pionScore, &b_proton_pionScore);
   fChain->SetBranchAddress("proton_protonScore", &proton_protonScore, &b_proton_protonScore);
   fChain->SetBranchAddress("proton_px", &proton_px, &b_proton_px);
   fChain->SetBranchAddress("proton_py", &proton_py, &b_proton_py);
   fChain->SetBranchAddress("proton_pz", &proton_pz, &b_proton_pz);
   fChain->SetBranchAddress("proton_theta", &proton_theta, &b_proton_theta);
   fChain->SetBranchAddress("proton_thetaX", &proton_thetaX, &b_proton_thetaX);
   fChain->SetBranchAddress("proton_thetaY", &proton_thetaY, &b_proton_thetaY);
   fChain->SetBranchAddress("reco_eventID", &reco_eventID, &b_reco_eventID);
   fChain->SetBranchAddress("time", &time, &b_time);
   fChain->SetBranchAddress("vertex_blob_evis", &vertex_blob_evis, &b_vertex_blob_evis);
   fChain->SetBranchAddress("vtx_fit_chi2", &vtx_fit_chi2, &b_vtx_fit_chi2);
   fChain->SetBranchAddress("vtx_x", &vtx_x, &b_vtx_x);
   fChain->SetBranchAddress("vtx_y", &vtx_y, &b_vtx_y);
   fChain->SetBranchAddress("vtx_z", &vtx_z, &b_vtx_z);
   fChain->SetBranchAddress("all_protons_kinked", all_protons_kinked, &b_all_protons_kinked);
   fChain->SetBranchAddress("all_protons_odMatch", all_protons_odMatch, &b_all_protons_odMatch);
   fChain->SetBranchAddress("g1dedx_cluster_occupancy_sz", &g1dedx_cluster_occupancy_sz, &b_g1dedx_cluster_occupancy_sz);
   fChain->SetBranchAddress("g1dedx_cluster_occupancy", g1dedx_cluster_occupancy, &b_g1dedx_cluster_occupancy);
   fChain->SetBranchAddress("g2dedx_cluster_occupancy_sz", &g2dedx_cluster_occupancy_sz, &b_g2dedx_cluster_occupancy_sz);
   fChain->SetBranchAddress("g2dedx_cluster_occupancy", g2dedx_cluster_occupancy, &b_g2dedx_cluster_occupancy);
   fChain->SetBranchAddress("nTracks_Secondary_Vtx_sz", &nTracks_Secondary_Vtx_sz, &b_nTracks_Secondary_Vtx_sz);
   fChain->SetBranchAddress("nTracks_Secondary_Vtx", nTracks_Secondary_Vtx, &b_nTracks_Secondary_Vtx);
   fChain->SetBranchAddress("all_protons_E", all_protons_E, &b_all_protons_E);
   fChain->SetBranchAddress("all_protons_KE", all_protons_KE, &b_all_protons_KE);
   fChain->SetBranchAddress("all_protons_LLRScore", all_protons_LLRScore, &b_all_protons_LLRScore);
   fChain->SetBranchAddress("all_protons_P", all_protons_P, &b_all_protons_P);
   fChain->SetBranchAddress("all_protons_chi2_ndf", all_protons_chi2_ndf, &b_all_protons_chi2_ndf);
   fChain->SetBranchAddress("all_protons_endPointX", all_protons_endPointX, &b_all_protons_endPointX);
   fChain->SetBranchAddress("all_protons_endPointY", all_protons_endPointY, &b_all_protons_endPointY);
   fChain->SetBranchAddress("all_protons_endPointZ", all_protons_endPointZ, &b_all_protons_endPointZ);
   fChain->SetBranchAddress("all_protons_length", all_protons_length, &b_all_protons_length);
   fChain->SetBranchAddress("all_protons_p_calCorrection", all_protons_p_calCorrection, &b_all_protons_p_calCorrection);
   fChain->SetBranchAddress("all_protons_p_dEdXTool", all_protons_p_dEdXTool, &b_all_protons_p_dEdXTool);
   fChain->SetBranchAddress("all_protons_p_visEnergy", all_protons_p_visEnergy, &b_all_protons_p_visEnergy);
   fChain->SetBranchAddress("all_protons_phi", all_protons_phi, &b_all_protons_phi);
   fChain->SetBranchAddress("all_protons_pionScore", all_protons_pionScore, &b_all_protons_pionScore);
   fChain->SetBranchAddress("all_protons_protonScore", all_protons_protonScore, &b_all_protons_protonScore);
   fChain->SetBranchAddress("all_protons_px", all_protons_px, &b_all_protons_px);
   fChain->SetBranchAddress("all_protons_py", all_protons_py, &b_all_protons_py);
   fChain->SetBranchAddress("all_protons_pz", all_protons_pz, &b_all_protons_pz);
   fChain->SetBranchAddress("all_protons_startPointX", all_protons_startPointX, &b_all_protons_startPointX);
   fChain->SetBranchAddress("all_protons_startPointY", all_protons_startPointY, &b_all_protons_startPointY);
   fChain->SetBranchAddress("all_protons_startPointZ", all_protons_startPointZ, &b_all_protons_startPointZ);
   fChain->SetBranchAddress("all_protons_theta", all_protons_theta, &b_all_protons_theta);
   fChain->SetBranchAddress("all_protons_thetaX", all_protons_thetaX, &b_all_protons_thetaX);
   fChain->SetBranchAddress("all_protons_thetaY", all_protons_thetaY, &b_all_protons_thetaY);
   fChain->SetBranchAddress("fit_vtx", fit_vtx, &b_fit_vtx);
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
   fChain->SetBranchAddress("truth_isGamma1_conv_inside", &truth_isGamma1_conv_inside, &b_truth_isGamma1_conv_inside);
   fChain->SetBranchAddress("truth_isGamma2_conv_inside", &truth_isGamma2_conv_inside, &b_truth_isGamma2_conv_inside);
   fChain->SetBranchAddress("truth_isSignal", &truth_isSignal, &b_truth_isSignal);
   fChain->SetBranchAddress("truth_isFidVol", &truth_isFidVol, &b_truth_isFidVol);
   fChain->SetBranchAddress("truth_isNC", &truth_isNC, &b_truth_isNC);
   fChain->SetBranchAddress("truth_ReconstructEvent", &truth_ReconstructEvent, &b_truth_ReconstructEvent);
   fChain->SetBranchAddress("truth_isBckg_NoPi0", &truth_isBckg_NoPi0, &b_truth_isBckg_NoPi0);
   fChain->SetBranchAddress("truth_isBckg_SinglePi0", &truth_isBckg_SinglePi0, &b_truth_isBckg_SinglePi0);
   fChain->SetBranchAddress("truth_isBckg_MultiPi0", &truth_isBckg_MultiPi0, &b_truth_isBckg_MultiPi0);
   fChain->SetBranchAddress("truth_isBckg_NC", &truth_isBckg_NC, &b_truth_isBckg_NC);
   fChain->SetBranchAddress("truth_isBckg_AntiNeutrino", &truth_isBckg_AntiNeutrino, &b_truth_isBckg_AntiNeutrino);
   fChain->SetBranchAddress("truth_isBckg_QELike", &truth_isBckg_QELike, &b_truth_isBckg_QELike);
   fChain->SetBranchAddress("truth_isBckg_SingleChargedPion", &truth_isBckg_SingleChargedPion, &b_truth_isBckg_SingleChargedPion);
   fChain->SetBranchAddress("truth_isBckg_SingleChargedPion_ChargeExchanged", &truth_isBckg_SingleChargedPion_ChargeExchanged, &b_truth_isBckg_SingleChargedPion_ChargeExchanged);
   fChain->SetBranchAddress("truth_isBckg_DoublePionWithPi0", &truth_isBckg_DoublePionWithPi0, &b_truth_isBckg_DoublePionWithPi0);
   fChain->SetBranchAddress("truth_isBckg_DoublePionWithoutPi0", &truth_isBckg_DoublePionWithoutPi0, &b_truth_isBckg_DoublePionWithoutPi0);
   fChain->SetBranchAddress("truth_isBckg_MultiPionWithPi0", &truth_isBckg_MultiPionWithPi0, &b_truth_isBckg_MultiPionWithPi0);
   fChain->SetBranchAddress("truth_isBckg_MultiPionWithoutPi0", &truth_isBckg_MultiPionWithoutPi0, &b_truth_isBckg_MultiPionWithoutPi0);
   fChain->SetBranchAddress("truth_isBckg_Other", &truth_isBckg_Other, &b_truth_isBckg_Other);
   fChain->SetBranchAddress("truth_isBckg_withMichel", &truth_isBckg_withMichel, &b_truth_isBckg_withMichel);
   fChain->SetBranchAddress("truth_Bckg_nOther", &truth_Bckg_nOther, &b_truth_Bckg_nOther);
   fChain->SetBranchAddress("truth_Bckg_nPi0_Primary", &truth_Bckg_nPi0_Primary, &b_truth_Bckg_nPi0_Primary);
   fChain->SetBranchAddress("truth_Bckg_nPi0_Secondary", &truth_Bckg_nPi0_Secondary, &b_truth_Bckg_nPi0_Secondary);
   fChain->SetBranchAddress("truth_Bckg_nPi0_Total", &truth_Bckg_nPi0_Total, &b_truth_Bckg_nPi0_Total);
   fChain->SetBranchAddress("truth_Bckg_nPiCharged", &truth_Bckg_nPiCharged, &b_truth_Bckg_nPiCharged);
   fChain->SetBranchAddress("truth_Bckg_nPiCharged_ChargeExchanged", &truth_Bckg_nPiCharged_ChargeExchanged, &b_truth_Bckg_nPiCharged_ChargeExchanged);
   fChain->SetBranchAddress("truth_N_FSParticles", &truth_N_FSParticles, &b_truth_N_FSParticles);
   fChain->SetBranchAddress("truth_N_other", &truth_N_other, &b_truth_N_other);
   fChain->SetBranchAddress("truth_N_pi0", &truth_N_pi0, &b_truth_N_pi0);
   fChain->SetBranchAddress("truth_N_proton", &truth_N_proton, &b_truth_N_proton);
   fChain->SetBranchAddress("truth_N_trueMichelElectrons", &truth_N_trueMichelElectrons, &b_truth_N_trueMichelElectrons);
   fChain->SetBranchAddress("truth_OneShower_evis_most_pdg", &truth_OneShower_evis_most_pdg, &b_truth_OneShower_evis_most_pdg);
   fChain->SetBranchAddress("truth_Rejected_unused_evis_most_pdg", &truth_Rejected_unused_evis_most_pdg, &b_truth_Rejected_unused_evis_most_pdg);
   fChain->SetBranchAddress("truth_ThreeShower_s1_evis_most_pdg", &truth_ThreeShower_s1_evis_most_pdg, &b_truth_ThreeShower_s1_evis_most_pdg);
   fChain->SetBranchAddress("truth_ThreeShower_s2_evis_most_pdg", &truth_ThreeShower_s2_evis_most_pdg, &b_truth_ThreeShower_s2_evis_most_pdg);
   fChain->SetBranchAddress("truth_ThreeShower_s3_evis_most_pdg", &truth_ThreeShower_s3_evis_most_pdg, &b_truth_ThreeShower_s3_evis_most_pdg);
   fChain->SetBranchAddress("truth_blob1_evis_most_pdg", &truth_blob1_evis_most_pdg, &b_truth_blob1_evis_most_pdg);
   fChain->SetBranchAddress("truth_blob2_evis_most_pdg", &truth_blob2_evis_most_pdg, &b_truth_blob2_evis_most_pdg);
   fChain->SetBranchAddress("truth_dispersed_unused_evis_most_pdg", &truth_dispersed_unused_evis_most_pdg, &b_truth_dispersed_unused_evis_most_pdg);
   fChain->SetBranchAddress("truth_pi0_GrandMother", &truth_pi0_GrandMother, &b_truth_pi0_GrandMother);
   fChain->SetBranchAddress("truth_pi0_GrandMotherStatus", &truth_pi0_GrandMotherStatus, &b_truth_pi0_GrandMotherStatus);
   fChain->SetBranchAddress("truth_pi0_Mother", &truth_pi0_Mother, &b_truth_pi0_Mother);
   fChain->SetBranchAddress("truth_pi0_MotherStatus", &truth_pi0_MotherStatus, &b_truth_pi0_MotherStatus);
   fChain->SetBranchAddress("truth_pi0_status", &truth_pi0_status, &b_truth_pi0_status);
   fChain->SetBranchAddress("truth_target_material", &truth_target_material, &b_truth_target_material);
   fChain->SetBranchAddress("truth_vertex_module", &truth_vertex_module, &b_truth_vertex_module);
   fChain->SetBranchAddress("truth_vertex_plane", &truth_vertex_plane, &b_truth_vertex_plane);
   fChain->SetBranchAddress("truth_vertex_unused_evis_most_pdg", &truth_vertex_unused_evis_most_pdg, &b_truth_vertex_unused_evis_most_pdg);
   fChain->SetBranchAddress("truth_OneShower_evis_muon", &truth_OneShower_evis_muon, &b_truth_OneShower_evis_muon);
   fChain->SetBranchAddress("truth_OneShower_evis_neutron", &truth_OneShower_evis_neutron, &b_truth_OneShower_evis_neutron);
   fChain->SetBranchAddress("truth_OneShower_evis_piminus", &truth_OneShower_evis_piminus, &b_truth_OneShower_evis_piminus);
   fChain->SetBranchAddress("truth_OneShower_evis_piplus", &truth_OneShower_evis_piplus, &b_truth_OneShower_evis_piplus);
   fChain->SetBranchAddress("truth_OneShower_evis_pizero", &truth_OneShower_evis_pizero, &b_truth_OneShower_evis_pizero);
   fChain->SetBranchAddress("truth_OneShower_evis_proton", &truth_OneShower_evis_proton, &b_truth_OneShower_evis_proton);
   fChain->SetBranchAddress("truth_OneShower_evis_total_norm", &truth_OneShower_evis_total_norm, &b_truth_OneShower_evis_total_norm);
   fChain->SetBranchAddress("truth_OneShower_evis_total_truth", &truth_OneShower_evis_total_truth, &b_truth_OneShower_evis_total_truth);
   fChain->SetBranchAddress("truth_Rejected_unused_evis_total_norm", &truth_Rejected_unused_evis_total_norm, &b_truth_Rejected_unused_evis_total_norm);
   fChain->SetBranchAddress("truth_Rejected_unused_evis_total_truth", &truth_Rejected_unused_evis_total_truth, &b_truth_Rejected_unused_evis_total_truth);
   fChain->SetBranchAddress("truth_ThreeShower_s1_evis_muon", &truth_ThreeShower_s1_evis_muon, &b_truth_ThreeShower_s1_evis_muon);
   fChain->SetBranchAddress("truth_ThreeShower_s1_evis_neutron", &truth_ThreeShower_s1_evis_neutron, &b_truth_ThreeShower_s1_evis_neutron);
   fChain->SetBranchAddress("truth_ThreeShower_s1_evis_piminus", &truth_ThreeShower_s1_evis_piminus, &b_truth_ThreeShower_s1_evis_piminus);
   fChain->SetBranchAddress("truth_ThreeShower_s1_evis_piplus", &truth_ThreeShower_s1_evis_piplus, &b_truth_ThreeShower_s1_evis_piplus);
   fChain->SetBranchAddress("truth_ThreeShower_s1_evis_pizero", &truth_ThreeShower_s1_evis_pizero, &b_truth_ThreeShower_s1_evis_pizero);
   fChain->SetBranchAddress("truth_ThreeShower_s1_evis_proton", &truth_ThreeShower_s1_evis_proton, &b_truth_ThreeShower_s1_evis_proton);
   fChain->SetBranchAddress("truth_ThreeShower_s1_evis_total_norm", &truth_ThreeShower_s1_evis_total_norm, &b_truth_ThreeShower_s1_evis_total_norm);
   fChain->SetBranchAddress("truth_ThreeShower_s1_evis_total_truth", &truth_ThreeShower_s1_evis_total_truth, &b_truth_ThreeShower_s1_evis_total_truth);
   fChain->SetBranchAddress("truth_ThreeShower_s2_evis_muon", &truth_ThreeShower_s2_evis_muon, &b_truth_ThreeShower_s2_evis_muon);
   fChain->SetBranchAddress("truth_ThreeShower_s2_evis_neutron", &truth_ThreeShower_s2_evis_neutron, &b_truth_ThreeShower_s2_evis_neutron);
   fChain->SetBranchAddress("truth_ThreeShower_s2_evis_piminus", &truth_ThreeShower_s2_evis_piminus, &b_truth_ThreeShower_s2_evis_piminus);
   fChain->SetBranchAddress("truth_ThreeShower_s2_evis_piplus", &truth_ThreeShower_s2_evis_piplus, &b_truth_ThreeShower_s2_evis_piplus);
   fChain->SetBranchAddress("truth_ThreeShower_s2_evis_pizero", &truth_ThreeShower_s2_evis_pizero, &b_truth_ThreeShower_s2_evis_pizero);
   fChain->SetBranchAddress("truth_ThreeShower_s2_evis_proton", &truth_ThreeShower_s2_evis_proton, &b_truth_ThreeShower_s2_evis_proton);
   fChain->SetBranchAddress("truth_ThreeShower_s2_evis_total_norm", &truth_ThreeShower_s2_evis_total_norm, &b_truth_ThreeShower_s2_evis_total_norm);
   fChain->SetBranchAddress("truth_ThreeShower_s2_evis_total_truth", &truth_ThreeShower_s2_evis_total_truth, &b_truth_ThreeShower_s2_evis_total_truth);
   fChain->SetBranchAddress("truth_ThreeShower_s3_evis_muon", &truth_ThreeShower_s3_evis_muon, &b_truth_ThreeShower_s3_evis_muon);
   fChain->SetBranchAddress("truth_ThreeShower_s3_evis_neutron", &truth_ThreeShower_s3_evis_neutron, &b_truth_ThreeShower_s3_evis_neutron);
   fChain->SetBranchAddress("truth_ThreeShower_s3_evis_piminus", &truth_ThreeShower_s3_evis_piminus, &b_truth_ThreeShower_s3_evis_piminus);
   fChain->SetBranchAddress("truth_ThreeShower_s3_evis_piplus", &truth_ThreeShower_s3_evis_piplus, &b_truth_ThreeShower_s3_evis_piplus);
   fChain->SetBranchAddress("truth_ThreeShower_s3_evis_pizero", &truth_ThreeShower_s3_evis_pizero, &b_truth_ThreeShower_s3_evis_pizero);
   fChain->SetBranchAddress("truth_ThreeShower_s3_evis_proton", &truth_ThreeShower_s3_evis_proton, &b_truth_ThreeShower_s3_evis_proton);
   fChain->SetBranchAddress("truth_ThreeShower_s3_evis_total_norm", &truth_ThreeShower_s3_evis_total_norm, &b_truth_ThreeShower_s3_evis_total_norm);
   fChain->SetBranchAddress("truth_ThreeShower_s3_evis_total_truth", &truth_ThreeShower_s3_evis_total_truth, &b_truth_ThreeShower_s3_evis_total_truth);
   fChain->SetBranchAddress("truth_allClusters_evis_pizero", &truth_allClusters_evis_pizero, &b_truth_allClusters_evis_pizero);
   fChain->SetBranchAddress("truth_blob1_evis_muon", &truth_blob1_evis_muon, &b_truth_blob1_evis_muon);
   fChain->SetBranchAddress("truth_blob1_evis_neutron", &truth_blob1_evis_neutron, &b_truth_blob1_evis_neutron);
   fChain->SetBranchAddress("truth_blob1_evis_piminus", &truth_blob1_evis_piminus, &b_truth_blob1_evis_piminus);
   fChain->SetBranchAddress("truth_blob1_evis_piplus", &truth_blob1_evis_piplus, &b_truth_blob1_evis_piplus);
   fChain->SetBranchAddress("truth_blob1_evis_pizero", &truth_blob1_evis_pizero, &b_truth_blob1_evis_pizero);
   fChain->SetBranchAddress("truth_blob1_evis_proton", &truth_blob1_evis_proton, &b_truth_blob1_evis_proton);
   fChain->SetBranchAddress("truth_blob1_evis_total_norm", &truth_blob1_evis_total_norm, &b_truth_blob1_evis_total_norm);
   fChain->SetBranchAddress("truth_blob1_evis_total_truth", &truth_blob1_evis_total_truth, &b_truth_blob1_evis_total_truth);
   fChain->SetBranchAddress("truth_blob2_evis_muon", &truth_blob2_evis_muon, &b_truth_blob2_evis_muon);
   fChain->SetBranchAddress("truth_blob2_evis_neutron", &truth_blob2_evis_neutron, &b_truth_blob2_evis_neutron);
   fChain->SetBranchAddress("truth_blob2_evis_piminus", &truth_blob2_evis_piminus, &b_truth_blob2_evis_piminus);
   fChain->SetBranchAddress("truth_blob2_evis_piplus", &truth_blob2_evis_piplus, &b_truth_blob2_evis_piplus);
   fChain->SetBranchAddress("truth_blob2_evis_pizero", &truth_blob2_evis_pizero, &b_truth_blob2_evis_pizero);
   fChain->SetBranchAddress("truth_blob2_evis_proton", &truth_blob2_evis_proton, &b_truth_blob2_evis_proton);
   fChain->SetBranchAddress("truth_blob2_evis_total_norm", &truth_blob2_evis_total_norm, &b_truth_blob2_evis_total_norm);
   fChain->SetBranchAddress("truth_blob2_evis_total_truth", &truth_blob2_evis_total_truth, &b_truth_blob2_evis_total_truth);
   fChain->SetBranchAddress("truth_dispersed_unused_evis_gamma", &truth_dispersed_unused_evis_gamma, &b_truth_dispersed_unused_evis_gamma);
   fChain->SetBranchAddress("truth_dispersed_unused_evis_muon", &truth_dispersed_unused_evis_muon, &b_truth_dispersed_unused_evis_muon);
   fChain->SetBranchAddress("truth_dispersed_unused_evis_neutron", &truth_dispersed_unused_evis_neutron, &b_truth_dispersed_unused_evis_neutron);
   fChain->SetBranchAddress("truth_dispersed_unused_evis_piminus", &truth_dispersed_unused_evis_piminus, &b_truth_dispersed_unused_evis_piminus);
   fChain->SetBranchAddress("truth_dispersed_unused_evis_piplus", &truth_dispersed_unused_evis_piplus, &b_truth_dispersed_unused_evis_piplus);
   fChain->SetBranchAddress("truth_dispersed_unused_evis_pizero", &truth_dispersed_unused_evis_pizero, &b_truth_dispersed_unused_evis_pizero);
   fChain->SetBranchAddress("truth_dispersed_unused_evis_proton", &truth_dispersed_unused_evis_proton, &b_truth_dispersed_unused_evis_proton);
   fChain->SetBranchAddress("truth_dispersed_unused_evis_total_norm", &truth_dispersed_unused_evis_total_norm, &b_truth_dispersed_unused_evis_total_norm);
   fChain->SetBranchAddress("truth_dispersed_unused_evis_total_truth", &truth_dispersed_unused_evis_total_truth, &b_truth_dispersed_unused_evis_total_truth);
   fChain->SetBranchAddress("truth_ecal_unused_evis_muon", &truth_ecal_unused_evis_muon, &b_truth_ecal_unused_evis_muon);
   fChain->SetBranchAddress("truth_ecal_unused_evis_neutron", &truth_ecal_unused_evis_neutron, &b_truth_ecal_unused_evis_neutron);
   fChain->SetBranchAddress("truth_ecal_unused_evis_piminus", &truth_ecal_unused_evis_piminus, &b_truth_ecal_unused_evis_piminus);
   fChain->SetBranchAddress("truth_ecal_unused_evis_piplus", &truth_ecal_unused_evis_piplus, &b_truth_ecal_unused_evis_piplus);
   fChain->SetBranchAddress("truth_ecal_unused_evis_pizero", &truth_ecal_unused_evis_pizero, &b_truth_ecal_unused_evis_pizero);
   fChain->SetBranchAddress("truth_ecal_unused_evis_proton", &truth_ecal_unused_evis_proton, &b_truth_ecal_unused_evis_proton);
   fChain->SetBranchAddress("truth_ecal_unused_evis_total_norm", &truth_ecal_unused_evis_total_norm, &b_truth_ecal_unused_evis_total_norm);
   fChain->SetBranchAddress("truth_ecal_unused_evis_total_truth", &truth_ecal_unused_evis_total_truth, &b_truth_ecal_unused_evis_total_truth);
   fChain->SetBranchAddress("truth_eventID", &truth_eventID, &b_truth_eventID);
   fChain->SetBranchAddress("truth_hcal_unused_evis_muon", &truth_hcal_unused_evis_muon, &b_truth_hcal_unused_evis_muon);
   fChain->SetBranchAddress("truth_hcal_unused_evis_neutron", &truth_hcal_unused_evis_neutron, &b_truth_hcal_unused_evis_neutron);
   fChain->SetBranchAddress("truth_hcal_unused_evis_piminus", &truth_hcal_unused_evis_piminus, &b_truth_hcal_unused_evis_piminus);
   fChain->SetBranchAddress("truth_hcal_unused_evis_piplus", &truth_hcal_unused_evis_piplus, &b_truth_hcal_unused_evis_piplus);
   fChain->SetBranchAddress("truth_hcal_unused_evis_pizero", &truth_hcal_unused_evis_pizero, &b_truth_hcal_unused_evis_pizero);
   fChain->SetBranchAddress("truth_hcal_unused_evis_proton", &truth_hcal_unused_evis_proton, &b_truth_hcal_unused_evis_proton);
   fChain->SetBranchAddress("truth_hcal_unused_evis_total_norm", &truth_hcal_unused_evis_total_norm, &b_truth_hcal_unused_evis_total_norm);
   fChain->SetBranchAddress("truth_hcal_unused_evis_total_truth", &truth_hcal_unused_evis_total_truth, &b_truth_hcal_unused_evis_total_truth);
   fChain->SetBranchAddress("truth_michelElectron_E", &truth_michelElectron_E, &b_truth_michelElectron_E);
   fChain->SetBranchAddress("truth_michelElectron_P", &truth_michelElectron_P, &b_truth_michelElectron_P);
   fChain->SetBranchAddress("truth_michelMuon_P", &truth_michelMuon_P, &b_truth_michelMuon_P);
   fChain->SetBranchAddress("truth_michelMuon_end_dist_vtx", &truth_michelMuon_end_dist_vtx, &b_truth_michelMuon_end_dist_vtx);
   fChain->SetBranchAddress("truth_michelMuon_length", &truth_michelMuon_length, &b_truth_michelMuon_length);
   fChain->SetBranchAddress("truth_michelPion_P", &truth_michelPion_P, &b_truth_michelPion_P);
   fChain->SetBranchAddress("truth_michelPion_begin_dist_vtx", &truth_michelPion_begin_dist_vtx, &b_truth_michelPion_begin_dist_vtx);
   fChain->SetBranchAddress("truth_michelPion_length", &truth_michelPion_length, &b_truth_michelPion_length);
   fChain->SetBranchAddress("truth_muon_P", &truth_muon_P, &b_truth_muon_P);
   fChain->SetBranchAddress("truth_muon_theta", &truth_muon_theta, &b_truth_muon_theta);
   fChain->SetBranchAddress("truth_other_unused_evis_muon", &truth_other_unused_evis_muon, &b_truth_other_unused_evis_muon);
   fChain->SetBranchAddress("truth_other_unused_evis_neutron", &truth_other_unused_evis_neutron, &b_truth_other_unused_evis_neutron);
   fChain->SetBranchAddress("truth_other_unused_evis_piminus", &truth_other_unused_evis_piminus, &b_truth_other_unused_evis_piminus);
   fChain->SetBranchAddress("truth_other_unused_evis_piplus", &truth_other_unused_evis_piplus, &b_truth_other_unused_evis_piplus);
   fChain->SetBranchAddress("truth_other_unused_evis_pizero", &truth_other_unused_evis_pizero, &b_truth_other_unused_evis_pizero);
   fChain->SetBranchAddress("truth_other_unused_evis_proton", &truth_other_unused_evis_proton, &b_truth_other_unused_evis_proton);
   fChain->SetBranchAddress("truth_other_unused_evis_total_norm", &truth_other_unused_evis_total_norm, &b_truth_other_unused_evis_total_norm);
   fChain->SetBranchAddress("truth_other_unused_evis_total_truth", &truth_other_unused_evis_total_truth, &b_truth_other_unused_evis_total_truth);
   fChain->SetBranchAddress("truth_pi0_P", &truth_pi0_P, &b_truth_pi0_P);
   fChain->SetBranchAddress("truth_pi0_theta", &truth_pi0_theta, &b_truth_pi0_theta);
   fChain->SetBranchAddress("truth_proton_P", &truth_proton_P, &b_truth_proton_P);
   fChain->SetBranchAddress("truth_proton_theta", &truth_proton_theta, &b_truth_proton_theta);
   fChain->SetBranchAddress("truth_total_captured_evis_pizero", &truth_total_captured_evis_pizero, &b_truth_total_captured_evis_pizero);
   fChain->SetBranchAddress("truth_total_captured_evis_total_norm", &truth_total_captured_evis_total_norm, &b_truth_total_captured_evis_total_norm);
   fChain->SetBranchAddress("truth_total_captured_evis_total_truth", &truth_total_captured_evis_total_truth, &b_truth_total_captured_evis_total_truth);
   fChain->SetBranchAddress("truth_vertex_unused_evis_gamma", &truth_vertex_unused_evis_gamma, &b_truth_vertex_unused_evis_gamma);
   fChain->SetBranchAddress("truth_vertex_unused_evis_muon", &truth_vertex_unused_evis_muon, &b_truth_vertex_unused_evis_muon);
   fChain->SetBranchAddress("truth_vertex_unused_evis_neutron", &truth_vertex_unused_evis_neutron, &b_truth_vertex_unused_evis_neutron);
   fChain->SetBranchAddress("truth_vertex_unused_evis_piminus", &truth_vertex_unused_evis_piminus, &b_truth_vertex_unused_evis_piminus);
   fChain->SetBranchAddress("truth_vertex_unused_evis_piplus", &truth_vertex_unused_evis_piplus, &b_truth_vertex_unused_evis_piplus);
   fChain->SetBranchAddress("truth_vertex_unused_evis_pizero", &truth_vertex_unused_evis_pizero, &b_truth_vertex_unused_evis_pizero);
   fChain->SetBranchAddress("truth_vertex_unused_evis_proton", &truth_vertex_unused_evis_proton, &b_truth_vertex_unused_evis_proton);
   fChain->SetBranchAddress("truth_vertex_unused_evis_total_norm", &truth_vertex_unused_evis_total_norm, &b_truth_vertex_unused_evis_total_norm);
   fChain->SetBranchAddress("truth_vertex_unused_evis_total_truth", &truth_vertex_unused_evis_total_truth, &b_truth_vertex_unused_evis_total_truth);
   fChain->SetBranchAddress("truth_gamma1_4P", truth_gamma1_4P, &b_truth_gamma1_4P);
   fChain->SetBranchAddress("truth_gamma1_final_pos", truth_gamma1_final_pos, &b_truth_gamma1_final_pos);
   fChain->SetBranchAddress("truth_gamma1_init_pos", truth_gamma1_init_pos, &b_truth_gamma1_init_pos);
   fChain->SetBranchAddress("truth_gamma2_4P", truth_gamma2_4P, &b_truth_gamma2_4P);
   fChain->SetBranchAddress("truth_gamma2_final_pos", truth_gamma2_final_pos, &b_truth_gamma2_final_pos);
   fChain->SetBranchAddress("truth_gamma2_init_pos", truth_gamma2_init_pos, &b_truth_gamma2_init_pos);
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
   fChain->SetBranchAddress("truth_michelMuon_endPoint", truth_michelMuon_endPoint, &b_truth_michelMuon_endPoint);
   fChain->SetBranchAddress("truth_muon_4P", truth_muon_4P, &b_truth_muon_4P);
   fChain->SetBranchAddress("truth_pi0_4P", truth_pi0_4P, &b_truth_pi0_4P);
   fChain->SetBranchAddress("truth_proton_4P", truth_proton_4P, &b_truth_proton_4P);
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
   fChain->SetBranchAddress("CCProtonPi0_ntrajMuonProng", &CCProtonPi0_ntrajMuonProng, &b_CCProtonPi0_ntrajMuonProng);
   fChain->SetBranchAddress("CCProtonPi0_r_minos_trk_vtx_plane", &CCProtonPi0_r_minos_trk_vtx_plane, &b_CCProtonPi0_r_minos_trk_vtx_plane);
   fChain->SetBranchAddress("CCProtonPi0_t_minos_trk_numFSMuons", &CCProtonPi0_t_minos_trk_numFSMuons, &b_CCProtonPi0_t_minos_trk_numFSMuons);
   fChain->SetBranchAddress("CCProtonPi0_t_minos_trk_primFSLeptonPDG", &CCProtonPi0_t_minos_trk_primFSLeptonPDG, &b_CCProtonPi0_t_minos_trk_primFSLeptonPDG);
   fChain->SetBranchAddress("CCProtonPi0_trajMuonProngPDG", &CCProtonPi0_trajMuonProngPDG, &b_CCProtonPi0_trajMuonProngPDG);
   fChain->SetBranchAddress("CCProtonPi0_trajMuonProngPrimary", &CCProtonPi0_trajMuonProngPrimary, &b_CCProtonPi0_trajMuonProngPrimary);
   fChain->SetBranchAddress("CCProtonPi0_Extra_Energy_Total", &CCProtonPi0_Extra_Energy_Total, &b_CCProtonPi0_Extra_Energy_Total);
   fChain->SetBranchAddress("CCProtonPi0_QSq", &CCProtonPi0_QSq, &b_CCProtonPi0_QSq);
   fChain->SetBranchAddress("CCProtonPi0_QSq_1Track_Alt", &CCProtonPi0_QSq_1Track_Alt, &b_CCProtonPi0_QSq_1Track_Alt);
   fChain->SetBranchAddress("CCProtonPi0_WSq", &CCProtonPi0_WSq, &b_CCProtonPi0_WSq);
   fChain->SetBranchAddress("CCProtonPi0_WSq_1Track_Alt", &CCProtonPi0_WSq_1Track_Alt, &b_CCProtonPi0_WSq_1Track_Alt);
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
   fChain->SetBranchAddress("CCProtonPi0_neutrino_E", &CCProtonPi0_neutrino_E, &b_CCProtonPi0_neutrino_E);
   fChain->SetBranchAddress("CCProtonPi0_neutrino_E_1Track_Alt", &CCProtonPi0_neutrino_E_1Track_Alt, &b_CCProtonPi0_neutrino_E_1Track_Alt);
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
   fChain->SetBranchAddress("CCProtonPi0_vertex_energy", &CCProtonPi0_vertex_energy, &b_CCProtonPi0_vertex_energy);
   fChain->SetBranchAddress("CCProtonPi0_isProtonInsideOD", CCProtonPi0_isProtonInsideOD, &b_CCProtonPi0_isProtonInsideOD);
   fChain->SetBranchAddress("CCProtonPi0_ntrajProtonProng", CCProtonPi0_ntrajProtonProng, &b_CCProtonPi0_ntrajProtonProng);
   fChain->SetBranchAddress("CCProtonPi0_trajProtonProngPDG", CCProtonPi0_trajProtonProngPDG, &b_CCProtonPi0_trajProtonProngPDG);
   fChain->SetBranchAddress("CCProtonPi0_trajProtonProngPrimary", CCProtonPi0_trajProtonProngPrimary, &b_CCProtonPi0_trajProtonProngPrimary);
   fChain->SetBranchAddress("CCProtonPi0_endProtonTrajMomentum", CCProtonPi0_endProtonTrajMomentum, &b_CCProtonPi0_endProtonTrajMomentum);
   fChain->SetBranchAddress("CCProtonPi0_endProtonTrajXPosition", CCProtonPi0_endProtonTrajXPosition, &b_CCProtonPi0_endProtonTrajXPosition);
   fChain->SetBranchAddress("CCProtonPi0_endProtonTrajYPosition", CCProtonPi0_endProtonTrajYPosition, &b_CCProtonPi0_endProtonTrajYPosition);
   fChain->SetBranchAddress("CCProtonPi0_endProtonTrajZPosition", CCProtonPi0_endProtonTrajZPosition, &b_CCProtonPi0_endProtonTrajZPosition);
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
   fChain->SetBranchAddress("mc_ppfx1_cvweight", &mc_ppfx1_cvweight, &b_mc_ppfx1_cvweight);
   fChain->SetBranchAddress("mc_gen1_cvweight_total", &mc_gen1_cvweight_total, &b_mc_gen1_cvweight_total);
   fChain->SetBranchAddress("gen1_wgt", &gen1_wgt, &b_gen1_wgt);
   fChain->SetBranchAddress("mc_gen1_cvweight_totalFlux", &mc_gen1_cvweight_totalFlux, &b_mc_gen1_cvweight_totalFlux);
   fChain->SetBranchAddress("mc_gen1_cvweight_NA49", &mc_gen1_cvweight_NA49, &b_mc_gen1_cvweight_NA49);
   fChain->SetBranchAddress("mc_wgt_Flux_BeamFocus_sz", &mc_wgt_Flux_BeamFocus_sz, &b_mc_wgt_Flux_BeamFocus_sz);
   fChain->SetBranchAddress("mc_wgt_Flux_BeamFocus", mc_wgt_Flux_BeamFocus, &b_mc_wgt_Flux_BeamFocus);
   fChain->SetBranchAddress("mc_wgt_gen1_Flux_Tertiary_sz", &mc_wgt_gen1_Flux_Tertiary_sz, &b_mc_wgt_gen1_Flux_Tertiary_sz);
   fChain->SetBranchAddress("mc_wgt_gen1_Flux_Tertiary", mc_wgt_gen1_Flux_Tertiary, &b_mc_wgt_gen1_Flux_Tertiary);
   fChain->SetBranchAddress("mc_wgt_gen1_Flux_NA49_sz", &mc_wgt_gen1_Flux_NA49_sz, &b_mc_wgt_gen1_Flux_NA49_sz);
   fChain->SetBranchAddress("mc_wgt_gen1_Flux_NA49", mc_wgt_gen1_Flux_NA49, &b_mc_wgt_gen1_Flux_NA49);
   fChain->SetBranchAddress("mc_wgt_Norm_sz", &mc_wgt_Norm_sz, &b_mc_wgt_Norm_sz);
   fChain->SetBranchAddress("mc_wgt_Norm", &mc_wgt_Norm, &b_mc_wgt_Norm);
   fChain->SetBranchAddress("mc_wgt_ppfx1_Total_sz", &mc_wgt_ppfx1_Total_sz, &b_mc_wgt_ppfx1_Total_sz);
   fChain->SetBranchAddress("mc_wgt_ppfx1_Total", mc_wgt_ppfx1_Total, &b_mc_wgt_ppfx1_Total);
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
