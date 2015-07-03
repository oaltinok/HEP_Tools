//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Jul  1 11:29:45 2015 by ROOT version 5.34/05
// from TChain NukeCCQETwoTrack/
//////////////////////////////////////////////////////////

#ifndef NukeCCQETwoTrack_default_h
#define NukeCCQETwoTrack_default_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "/grid/fermiapp/minerva/software_releases/lcgcmake/build/lcg_61/projects/ROOT-5.34.05/src/ROOT/5.34.05/cint/cint/lib/prec_stl/vector"

// Fixed size dimensions of array or collections stored in the TTree if any.

class NukeCCQETwoTrack_default {
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
   Int_t           CreatedShortTracks;
   Int_t           FailContainedProng;
   Int_t           FailExitingProng;
   Int_t           FailFidVolume;
   Int_t           FailOutTracks;
   Int_t           FailRefitFidVolume;
   Int_t           FailShortOutTrack;
   Int_t           NoInteractionVertex;
   Int_t           NullVertex;
   Int_t           UnattachedProngsWithTracks;
   Int_t           broken_track_most_us_plane;
   Int_t           genie_n_charms;
   Int_t           genie_n_heavy_baryons;
   Int_t           genie_n_kaons;
   Int_t           genie_n_mesons;
   Int_t           genie_n_muons;
   Int_t           genie_n_neutrinos;
   Int_t           genie_n_neutrons;
   Int_t           genie_n_others;
   Int_t           genie_n_particles;
   Int_t           genie_n_photons;
   Int_t           genie_n_pi_zeros;
   Int_t           genie_n_pions;
   Int_t           genie_n_protons;
   Int_t           intraNukeDeltaPlusPlusDecay;
   Int_t           intraNukeNParticles;
   Int_t           intraNukeNeutronQuasiElasticScatter;
   Int_t           intraNukeOtherProcess;
   Int_t           is_true_muon;
   Int_t           muon_enters_front;
   Int_t           n_odClusters;
   Int_t           n_odClustersWithTimeCut;
   Int_t           ntracks;
   Int_t           part_response_recoil_nClus;
   Int_t           passVertexZCut;
   Int_t           phys_energy_in_road_downstream_nplanes;
   Int_t           phys_energy_in_road_upstream_nplanes;
   Int_t           phys_n_dead_discr_pair;
   Int_t           phys_n_dead_discr_pair_in_prim_track_region;
   Int_t           phys_n_dead_discr_pair_two_mod_downstream_prim_track;
   Int_t           phys_n_dead_discr_pair_two_mod_upstream_prim_vtx;
   Int_t           phys_n_dead_discr_pair_upstream_prim_track_proj;
   Int_t           phys_vertex_is_fiducial;
   Int_t           reco_muon_topology;
   Int_t           timeSlice;
   Int_t           vtx_fit_converged;
   Double_t        endPointEnergy;
   Double_t        energy_from_mc;
   Double_t        energy_from_mc_fraction;
   Double_t        energy_from_mc_fraction_of_highest;
   Double_t        hadronic_energy;
   Double_t        intraNukeProtonMomentum;
   Double_t        isolatedEnergy;
   Double_t        isolatedEnergy_ecal;
   Double_t        isolatedEnergy_hcal;
   Double_t        isolatedEnergy_targets;
   Double_t        isolatedEnergy_tracker;
   Double_t        muonFuzzEnergy;
   Double_t        muon_phi;
   Double_t        muon_theta;
   Double_t        muon_thetaX;
   Double_t        muon_thetaY;
   Double_t        numi_horn_curr;
   Double_t        numi_pot;
   Double_t        numi_x;
   Double_t        numi_x_width;
   Double_t        numi_y;
   Double_t        numi_y_width;
   Double_t        odEnergy;
   Double_t        odEnergyWithTimeCut;
   Double_t        part_response_recoil_em_id;
   Double_t        part_response_recoil_em_id_err;
   Double_t        part_response_recoil_em_od;
   Double_t        part_response_recoil_em_od_err;
   Double_t        part_response_recoil_high_neutron_id;
   Double_t        part_response_recoil_high_neutron_id_err;
   Double_t        part_response_recoil_high_neutron_od;
   Double_t        part_response_recoil_high_neutron_od_err;
   Double_t        part_response_recoil_low_neutron_id;
   Double_t        part_response_recoil_low_neutron_id_err;
   Double_t        part_response_recoil_low_neutron_od;
   Double_t        part_response_recoil_low_neutron_od_err;
   Double_t        part_response_recoil_meson_id;
   Double_t        part_response_recoil_meson_id_err;
   Double_t        part_response_recoil_meson_od;
   Double_t        part_response_recoil_meson_od_err;
   Double_t        part_response_recoil_mid_neutron_id;
   Double_t        part_response_recoil_mid_neutron_id_err;
   Double_t        part_response_recoil_mid_neutron_od;
   Double_t        part_response_recoil_mid_neutron_od_err;
   Double_t        part_response_recoil_muon_id;
   Double_t        part_response_recoil_muon_id_err;
   Double_t        part_response_recoil_muon_od;
   Double_t        part_response_recoil_muon_od_err;
   Double_t        part_response_recoil_other_id;
   Double_t        part_response_recoil_other_id_err;
   Double_t        part_response_recoil_other_od;
   Double_t        part_response_recoil_other_od_err;
   Double_t        part_response_recoil_proton_id;
   Double_t        part_response_recoil_proton_id_err;
   Double_t        part_response_recoil_proton_od;
   Double_t        part_response_recoil_proton_od_err;
   Double_t        part_response_recoil_xtalk_id;
   Double_t        part_response_recoil_xtalk_id_err;
   Double_t        part_response_recoil_xtalk_od;
   Double_t        part_response_recoil_xtalk_od_err;
   Double_t        phys_energy_dispersed;
   Double_t        phys_energy_in_road_downstream;
   Double_t        phys_energy_in_road_upstream;
   Double_t        phys_energy_unattached;
   Double_t        prim_vtx_smallest_opening_angle;
   Double_t        primaryVertexEnergy;
   Double_t        protonFuzzEnergy;
   Double_t        secondaryVertexEnergy;
   Double_t        vtx_fit_chi2;
   Double_t        vtxprong_energy_cal;
   Double_t        vtxprong_energy_visible;
   Double_t        vtxprong_energy_visible_isoblobs;
   Double_t        vtxprong_energy_visible_vtxblobs;
   Int_t           has_michel_category_sz;
   Int_t           has_michel_category[5];   //[has_michel_category_sz]
   Int_t           has_michel_in_vertex_point_sz;
   Int_t           has_michel_in_vertex_point[5];   //[has_michel_in_vertex_point_sz]
   Int_t           has_michel_ndigits_sz;
   Int_t           has_michel_ndigits[5];   //[has_michel_ndigits_sz]
   Int_t           has_michel_vertex_type_sz;
   Int_t           has_michel_vertex_type[5];   //[has_michel_vertex_type_sz]
   Int_t           proton_enters_front[10];
   Int_t           n_vtxprong_isoblobs;
   Int_t           vtxprong_isoblob_nclusters[2];   //[n_vtxprong_isoblobs]
   Int_t           n_vtxprong_vtxblobs;
   Int_t           vtxprong_vtxblob_nclusters[1];   //[n_vtxprong_vtxblobs]
   Int_t           clusterU_Angle_sz;
   Double_t        clusterU_Angle[80];   //[clusterU_Angle_sz]
   Int_t           clusterU_Radius_sz;
   Double_t        clusterU_Radius[80];   //[clusterU_Radius_sz]
   Int_t           clusterU_timeDiff_sz;
   Double_t        clusterU_timeDiff[80];   //[clusterU_timeDiff_sz]
   Int_t           clusterU_viewDist_sz;
   Double_t        clusterU_viewDist[80];   //[clusterU_viewDist_sz]
   Int_t           clusterU_visE_binned_sz;
   Double_t        clusterU_visE_binned[10];   //[clusterU_visE_binned_sz]
   Int_t           clusterU_visEnergy_sz;
   Double_t        clusterU_visEnergy[80];   //[clusterU_visEnergy_sz]
   Int_t           clusterU_zDist_sz;
   Double_t        clusterU_zDist[80];   //[clusterU_zDist_sz]
   Int_t           clusterV_Angle_sz;
   Double_t        clusterV_Angle[86];   //[clusterV_Angle_sz]
   Int_t           clusterV_Radius_sz;
   Double_t        clusterV_Radius[86];   //[clusterV_Radius_sz]
   Int_t           clusterV_timeDiff_sz;
   Double_t        clusterV_timeDiff[86];   //[clusterV_timeDiff_sz]
   Int_t           clusterV_viewDist_sz;
   Double_t        clusterV_viewDist[86];   //[clusterV_viewDist_sz]
   Int_t           clusterV_visE_binned_sz;
   Double_t        clusterV_visE_binned[10];   //[clusterV_visE_binned_sz]
   Int_t           clusterV_visEnergy_sz;
   Double_t        clusterV_visEnergy[86];   //[clusterV_visEnergy_sz]
   Int_t           clusterV_zDist_sz;
   Double_t        clusterV_zDist[86];   //[clusterV_zDist_sz]
   Int_t           clusterX_Angle_sz;
   Double_t        clusterX_Angle[185];   //[clusterX_Angle_sz]
   Int_t           clusterX_Radius_sz;
   Double_t        clusterX_Radius[185];   //[clusterX_Radius_sz]
   Int_t           clusterX_timeDiff_sz;
   Double_t        clusterX_timeDiff[185];   //[clusterX_timeDiff_sz]
   Int_t           clusterX_viewDist_sz;
   Double_t        clusterX_viewDist[185];   //[clusterX_viewDist_sz]
   Int_t           clusterX_visE_binned_sz;
   Double_t        clusterX_visE_binned[10];   //[clusterX_visE_binned_sz]
   Int_t           clusterX_visEnergy_sz;
   Double_t        clusterX_visEnergy[185];   //[clusterX_visEnergy_sz]
   Int_t           clusterX_zDist_sz;
   Double_t        clusterX_zDist[185];   //[clusterX_zDist_sz]
   Double_t        fit_vtx[3];
   Int_t           has_michel_distance_sz;
   Double_t        has_michel_distance[5];   //[has_michel_distance_sz]
   Int_t           has_michel_energy_sz;
   Double_t        has_michel_energy[5];   //[has_michel_energy_sz]
   Int_t           has_michel_time_diff_sz;
   Double_t        has_michel_time_diff[5];   //[has_michel_time_diff_sz]
   Double_t        intraNukeProtonMomentumVec[4];
   Double_t        reco_muon_vertex[3];
   Double_t        vtxprong_isoblob_visenergy[2];   //[n_vtxprong_isoblobs]
   Double_t        vtxprong_vtxblob_visenergy[1];   //[n_vtxprong_vtxblobs]
   Bool_t          truth_has_physics_event;
   Int_t           truth_has_michel_electron;
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
   Int_t           truth_has_michel_from_pion_minus_momentum_sz;
   Double_t        truth_has_michel_from_pion_minus_momentum[1];   //[truth_has_michel_from_pion_minus_momentum_sz]
   Int_t           truth_has_michel_from_pion_plus_momentum_sz;
   Double_t        truth_has_michel_from_pion_plus_momentum[13];   //[truth_has_michel_from_pion_plus_momentum_sz]
   Int_t           NukeCCQETwoTrack_nuFlavor;
   Int_t           NukeCCQETwoTrack_nuHelicity;
   Int_t           NukeCCQETwoTrack_intCurrent;
   Int_t           NukeCCQETwoTrack_intType;
   Double_t        NukeCCQETwoTrack_E;
   Double_t        NukeCCQETwoTrack_Q2;
   Double_t        NukeCCQETwoTrack_x;
   Double_t        NukeCCQETwoTrack_y;
   Double_t        NukeCCQETwoTrack_W;
   Double_t        NukeCCQETwoTrack_score;
   Double_t        NukeCCQETwoTrack_leptonE[4];
   Double_t        NukeCCQETwoTrack_vtx[4];
   Bool_t          NukeCCQETwoTrack_minos_trk_is_contained;
   Bool_t          NukeCCQETwoTrack_minos_trk_is_ok;
   Bool_t          NukeCCQETwoTrack_minos_used_range;
   Bool_t          NukeCCQETwoTrack_minos_used_curvature;
   Int_t           NukeCCQETwoTrack_isMuonInsideOD;
   Int_t           NukeCCQETwoTrack_minos_trk_end_plane;
   Int_t           NukeCCQETwoTrack_minos_trk_quality;
   Int_t           NukeCCQETwoTrack_muon_down_hcal;
   Int_t           NukeCCQETwoTrack_muon_minos_stub;
   Int_t           NukeCCQETwoTrack_muon_minos_track;
   Int_t           NukeCCQETwoTrack_muon_odLastFrame;
   Int_t           NukeCCQETwoTrack_muon_odLastStory;
   Int_t           NukeCCQETwoTrack_muon_od_track;
   Int_t           NukeCCQETwoTrack_muon_side_ecal;
   Int_t           NukeCCQETwoTrack_muon_trk_pat_history;
   Int_t           NukeCCQETwoTrack_ntrajMuonProng;
   Int_t           NukeCCQETwoTrack_r_minos_trk_vtx_plane;
   Int_t           NukeCCQETwoTrack_t_minos_trk_numFSMuons;
   Int_t           NukeCCQETwoTrack_t_minos_trk_primFSLeptonPDG;
   Int_t           NukeCCQETwoTrack_targetID;
   Int_t           NukeCCQETwoTrack_targetZ;
   Int_t           NukeCCQETwoTrack_trajMuonProngPDG;
   Int_t           NukeCCQETwoTrack_trajMuonProngPrimary;
   Double_t        NukeCCQETwoTrack_endMuonTrajMomentum;
   Double_t        NukeCCQETwoTrack_endMuonTrajXPosition;
   Double_t        NukeCCQETwoTrack_endMuonTrajYPosition;
   Double_t        NukeCCQETwoTrack_endMuonTrajZPosition;
   Double_t        NukeCCQETwoTrack_minos_trk_bave;
   Double_t        NukeCCQETwoTrack_minos_trk_chi2;
   Double_t        NukeCCQETwoTrack_minos_trk_end_u;
   Double_t        NukeCCQETwoTrack_minos_trk_end_v;
   Double_t        NukeCCQETwoTrack_minos_trk_end_x;
   Double_t        NukeCCQETwoTrack_minos_trk_end_y;
   Double_t        NukeCCQETwoTrack_minos_trk_end_z;
   Double_t        NukeCCQETwoTrack_minos_trk_eqp;
   Double_t        NukeCCQETwoTrack_minos_trk_eqp_qp;
   Double_t        NukeCCQETwoTrack_minos_trk_fit_pass;
   Double_t        NukeCCQETwoTrack_minos_trk_ndf;
   Double_t        NukeCCQETwoTrack_minos_trk_p;
   Double_t        NukeCCQETwoTrack_minos_trk_p_curvature;
   Double_t        NukeCCQETwoTrack_minos_trk_p_range;
   Double_t        NukeCCQETwoTrack_minos_trk_qp;
   Double_t        NukeCCQETwoTrack_minos_trk_vtx_x;
   Double_t        NukeCCQETwoTrack_minos_trk_vtx_y;
   Double_t        NukeCCQETwoTrack_minos_trk_vtx_z;
   Double_t        NukeCCQETwoTrack_muon_enu;
   Double_t        NukeCCQETwoTrack_muon_odClustersAvgTime;
   Double_t        NukeCCQETwoTrack_muon_odElossMomentum;
   Double_t        NukeCCQETwoTrack_muon_odEndX;
   Double_t        NukeCCQETwoTrack_muon_odEndY;
   Double_t        NukeCCQETwoTrack_muon_odEndZ;
   Double_t        NukeCCQETwoTrack_muon_odFaceX;
   Double_t        NukeCCQETwoTrack_muon_odFaceY;
   Double_t        NukeCCQETwoTrack_muon_odFaceZ;
   Double_t        NukeCCQETwoTrack_muon_odLastClusZ;
   Double_t        NukeCCQETwoTrack_muon_odStopDistMomentum;
   Double_t        NukeCCQETwoTrack_muon_odTrackAvgTime;
   Double_t        NukeCCQETwoTrack_muon_phi;
   Double_t        NukeCCQETwoTrack_muon_q2;
   Double_t        NukeCCQETwoTrack_muon_score;
   Double_t        NukeCCQETwoTrack_muon_theta;
   Double_t        NukeCCQETwoTrack_muon_thetaX;
   Double_t        NukeCCQETwoTrack_muon_thetaY;
   Double_t        NukeCCQETwoTrack_r_minos_trk_bdL;
   Double_t        NukeCCQETwoTrack_r_minos_trk_end_dcosx;
   Double_t        NukeCCQETwoTrack_r_minos_trk_end_dcosy;
   Double_t        NukeCCQETwoTrack_r_minos_trk_end_dcosz;
   Double_t        NukeCCQETwoTrack_r_minos_trk_vtx_dcosx;
   Double_t        NukeCCQETwoTrack_r_minos_trk_vtx_dcosy;
   Double_t        NukeCCQETwoTrack_r_minos_trk_vtx_dcosz;
   Double_t        NukeCCQETwoTrack_t_minos_trk_primFSLepMinosInitProjPx;
   Double_t        NukeCCQETwoTrack_t_minos_trk_primFSLepMinosInitProjPy;
   Double_t        NukeCCQETwoTrack_t_minos_trk_primFSLepMinosInitProjPz;
   Double_t        NukeCCQETwoTrack_t_minos_trk_primFSLepMinosInitProjX;
   Double_t        NukeCCQETwoTrack_t_minos_trk_primFSLepMinosInitProjY;
   Double_t        NukeCCQETwoTrack_t_minos_trk_primFSLepMinosInitProjZ;
   Double_t        NukeCCQETwoTrack_t_minos_trk_primFSLepMnvFinalPx;
   Double_t        NukeCCQETwoTrack_t_minos_trk_primFSLepMnvFinalPy;
   Double_t        NukeCCQETwoTrack_t_minos_trk_primFSLepMnvFinalPz;
   Double_t        NukeCCQETwoTrack_t_minos_trk_primFSLepMnvFinalX;
   Double_t        NukeCCQETwoTrack_t_minos_trk_primFSLepMnvFinalY;
   Double_t        NukeCCQETwoTrack_t_minos_trk_primFSLepMnvFinalZ;
   Double_t        NukeCCQETwoTrack_t_minos_trk_primFSLepMnvInitPx;
   Double_t        NukeCCQETwoTrack_t_minos_trk_primFSLepMnvInitPy;
   Double_t        NukeCCQETwoTrack_t_minos_trk_primFSLepMnvInitPz;
   Double_t        NukeCCQETwoTrack_t_minos_trk_primFSLepMnvInitX;
   Double_t        NukeCCQETwoTrack_t_minos_trk_primFSLepMnvInitY;
   Double_t        NukeCCQETwoTrack_t_minos_trk_primFSLepMnvInitZ;
   Double_t        NukeCCQETwoTrack_targetZPos;
   Double_t        NukeCCQETwoTrack_trajMuonPhi;
   Double_t        NukeCCQETwoTrack_trajMuonProngMomentum;
   Double_t        NukeCCQETwoTrack_trajMuonTheta;
   Int_t           NukeCCQETwoTrack_isProtonInsideOD[10];
   Int_t           NukeCCQETwoTrack_ntrajProngProng[10];
   Int_t           NukeCCQETwoTrack_proton_kinked[10];
   Int_t           NukeCCQETwoTrack_proton_odMatch[10];
   Int_t           NukeCCQETwoTrack_proton_trk_pat_history[10];
   Int_t           NukeCCQETwoTrack_trajProtonProngPDG[10];
   Int_t           NukeCCQETwoTrack_trajProtonProngPrimary[10];
   Double_t        NukeCCQETwoTrack_coplanarAngle[10];
   Double_t        NukeCCQETwoTrack_endProtonTrajMomentum[10];
   Double_t        NukeCCQETwoTrack_endProtonTrajXPosition[10];
   Double_t        NukeCCQETwoTrack_endProtonTrajYPosition[10];
   Double_t        NukeCCQETwoTrack_endProtonTrajZPosition[10];
   Double_t        NukeCCQETwoTrack_muon_endPoint[3];
   Double_t        NukeCCQETwoTrack_muon_startPoint[3];
   Double_t        NukeCCQETwoTrack_open_angle[10];
   Double_t        NukeCCQETwoTrack_pion_chi2_ndf[10];
   Double_t        NukeCCQETwoTrack_pion_score[10];
   Double_t        NukeCCQETwoTrack_pion_score1[10];
   Double_t        NukeCCQETwoTrack_pion_score2[10];
   Double_t        NukeCCQETwoTrack_proton_E[10];
   Double_t        NukeCCQETwoTrack_proton_chi2_ndf[10];
   Double_t        NukeCCQETwoTrack_proton_ekin[10];
   Double_t        NukeCCQETwoTrack_proton_endPointX[10];
   Double_t        NukeCCQETwoTrack_proton_endPointY[10];
   Double_t        NukeCCQETwoTrack_proton_endPointZ[10];
   Double_t        NukeCCQETwoTrack_proton_enu[10];
   Double_t        NukeCCQETwoTrack_proton_momentum_shift_BetheBloch_Down[10];
   Double_t        NukeCCQETwoTrack_proton_momentum_shift_BetheBloch_Up[10];
   Double_t        NukeCCQETwoTrack_proton_momentum_shift_Birks[10];
   Double_t        NukeCCQETwoTrack_proton_momentum_shift_MEU_Down[10];
   Double_t        NukeCCQETwoTrack_proton_momentum_shift_MEU_Up[10];
   Double_t        NukeCCQETwoTrack_proton_momentum_shift_Mass_Down[10];
   Double_t        NukeCCQETwoTrack_proton_momentum_shift_Mass_Up[10];
   Double_t        NukeCCQETwoTrack_proton_p[10];
   Double_t        NukeCCQETwoTrack_proton_p_calCorrection[10];
   Double_t        NukeCCQETwoTrack_proton_p_dEdXTool[10];
   Double_t        NukeCCQETwoTrack_proton_p_visEnergy[10];
   Double_t        NukeCCQETwoTrack_proton_phi[10];
   Double_t        NukeCCQETwoTrack_proton_px[10];
   Double_t        NukeCCQETwoTrack_proton_py[10];
   Double_t        NukeCCQETwoTrack_proton_pz[10];
   Double_t        NukeCCQETwoTrack_proton_q2[10];
   Double_t        NukeCCQETwoTrack_proton_score[10];
   Double_t        NukeCCQETwoTrack_proton_score1[10];
   Double_t        NukeCCQETwoTrack_proton_score1_shift_BetheBloch_Down[10];
   Double_t        NukeCCQETwoTrack_proton_score1_shift_BetheBloch_Up[10];
   Double_t        NukeCCQETwoTrack_proton_score1_shift_Birks[10];
   Double_t        NukeCCQETwoTrack_proton_score1_shift_MEU_Down[10];
   Double_t        NukeCCQETwoTrack_proton_score1_shift_MEU_Up[10];
   Double_t        NukeCCQETwoTrack_proton_score1_shift_Mass_Down[10];
   Double_t        NukeCCQETwoTrack_proton_score1_shift_Mass_Up[10];
   Double_t        NukeCCQETwoTrack_proton_score2[10];
   Double_t        NukeCCQETwoTrack_proton_startPointX[10];
   Double_t        NukeCCQETwoTrack_proton_startPointY[10];
   Double_t        NukeCCQETwoTrack_proton_startPointZ[10];
   Double_t        NukeCCQETwoTrack_proton_theta[10];
   Double_t        NukeCCQETwoTrack_proton_thetaX[10];
   Double_t        NukeCCQETwoTrack_proton_thetaY[10];
   Double_t        NukeCCQETwoTrack_sys_muon_curve_energy_shift[2];
   Double_t        NukeCCQETwoTrack_sys_muon_energy_shift[2];
   Double_t        NukeCCQETwoTrack_sys_muon_minerva_energy_shift[2];
   Double_t        NukeCCQETwoTrack_sys_muon_qSquared_shift[2];
   Double_t        NukeCCQETwoTrack_sys_muon_range_energy_shift[2];
   Double_t        NukeCCQETwoTrack_sys_muon_wSquared_shift[2];
   Double_t        NukeCCQETwoTrack_sys_muon_xbj_shift[2];
   Double_t        NukeCCQETwoTrack_sys_muon_y_shift[2];
   Double_t        NukeCCQETwoTrack_sys_nu_energy_shift[2];
   Double_t        NukeCCQETwoTrack_sys_recoil_energy_shift[2];
   Double_t        NukeCCQETwoTrack_sys_recoil_qSquared_shift[2];
   Double_t        NukeCCQETwoTrack_sys_recoil_wSquared_shift[2];
   Double_t        NukeCCQETwoTrack_sys_recoil_xbj_shift[2];
   Double_t        NukeCCQETwoTrack_sys_recoil_y_shift[2];
   Double_t        NukeCCQETwoTrack_sys_total_qSquared_shift[2];
   Double_t        NukeCCQETwoTrack_sys_total_wSquared_shift[2];
   Double_t        NukeCCQETwoTrack_sys_total_xbj_shift[2];
   Double_t        NukeCCQETwoTrack_sys_total_y_shift[2];
   Double_t        NukeCCQETwoTrack_trajProtonPhi[10];
   Double_t        NukeCCQETwoTrack_trajProtonProngMomentum[10];
   Double_t        NukeCCQETwoTrack_trajProtonTheta[10];
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
   Double_t        mc_FSPartPx[152];   //[mc_nFSPart]
   Double_t        mc_FSPartPy[152];   //[mc_nFSPart]
   Double_t        mc_FSPartPz[152];   //[mc_nFSPart]
   Double_t        mc_FSPartE[152];   //[mc_nFSPart]
   Int_t           mc_FSPartPDG[152];   //[mc_nFSPart]
   Int_t           mc_er_nPart;
   Int_t           mc_er_ID[181];   //[mc_er_nPart]
   Int_t           mc_er_status[181];   //[mc_er_nPart]
   Double_t        mc_er_posInNucX[181];   //[mc_er_nPart]
   Double_t        mc_er_posInNucY[181];   //[mc_er_nPart]
   Double_t        mc_er_posInNucZ[181];   //[mc_er_nPart]
   Double_t        mc_er_Px[181];   //[mc_er_nPart]
   Double_t        mc_er_Py[181];   //[mc_er_nPart]
   Double_t        mc_er_Pz[181];   //[mc_er_nPart]
   Double_t        mc_er_E[181];   //[mc_er_nPart]
   Int_t           mc_er_FD[181];   //[mc_er_nPart]
   Int_t           mc_er_LD[181];   //[mc_er_nPart]
   Int_t           mc_er_mother[181];   //[mc_er_nPart]
   Int_t           mc_fr_nNuAncestorIDs;
   Int_t           mc_fr_nuAncestorIDs[8];   //[mc_fr_nNuAncestorIDs]
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
   Int_t           mc_wgt_Norm_sz;
   Double_t        mc_wgt_Norm[1];   //[mc_wgt_Norm_sz]
   Int_t           mc_wgt_ppfx_MIPPNumiYields_sz;
   Double_t        mc_wgt_ppfx_MIPPNumiYields[1];   //[mc_wgt_ppfx_MIPPNumiYields_sz]
   Int_t           mc_wgt_ppfx_TargetAttenuation_sz;
   Double_t        mc_wgt_ppfx_TargetAttenuation[1];   //[mc_wgt_ppfx_TargetAttenuation_sz]
   Int_t           mc_wgt_ppfx_NA49_sz;
   Double_t        mc_wgt_ppfx_NA49[1];   //[mc_wgt_ppfx_NA49_sz]
   Int_t           mc_wgt_ppfx_MIPPKaonsYields_sz;
   Double_t        mc_wgt_ppfx_MIPPKaonsYields[1];   //[mc_wgt_ppfx_MIPPKaonsYields_sz]
   Int_t           mc_wgt_ppfx_MIPPThinTarget_sz;
   Double_t        mc_wgt_ppfx_MIPPThinTarget[1];   //[mc_wgt_ppfx_MIPPThinTarget_sz]
   Int_t           mc_wgt_ppfx_Absorption_sz;
   Double_t        mc_wgt_ppfx_Absorption[1];   //[mc_wgt_ppfx_Absorption_sz]
   Int_t           mc_wgt_ppfx_Others_sz;
   Double_t        mc_wgt_ppfx_Others[1];   //[mc_wgt_ppfx_Others_sz]
   Int_t           n_prongs;
   Int_t           prong_nParticles[6];   //[n_prongs]
   Double_t        prong_part_score[6];   //[n_prongs]
   Double_t        prong_part_mass[6];   //[n_prongs]
   Int_t           prong_part_charge[6];   //[n_prongs]
   Int_t           prong_part_pid[6];   //[n_prongs]
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
   TBranch        *b_CreatedShortTracks;   //!
   TBranch        *b_FailContainedProng;   //!
   TBranch        *b_FailExitingProng;   //!
   TBranch        *b_FailFidVolume;   //!
   TBranch        *b_FailOutTracks;   //!
   TBranch        *b_FailRefitFidVolume;   //!
   TBranch        *b_FailShortOutTrack;   //!
   TBranch        *b_NoInteractionVertex;   //!
   TBranch        *b_NullVertex;   //!
   TBranch        *b_UnattachedProngsWithTracks;   //!
   TBranch        *b_broken_track_most_us_plane;   //!
   TBranch        *b_genie_n_charms;   //!
   TBranch        *b_genie_n_heavy_baryons;   //!
   TBranch        *b_genie_n_kaons;   //!
   TBranch        *b_genie_n_mesons;   //!
   TBranch        *b_genie_n_muons;   //!
   TBranch        *b_genie_n_neutrinos;   //!
   TBranch        *b_genie_n_neutrons;   //!
   TBranch        *b_genie_n_others;   //!
   TBranch        *b_genie_n_particles;   //!
   TBranch        *b_genie_n_photons;   //!
   TBranch        *b_genie_n_pi_zeros;   //!
   TBranch        *b_genie_n_pions;   //!
   TBranch        *b_genie_n_protons;   //!
   TBranch        *b_intraNukeDeltaPlusPlusDecay;   //!
   TBranch        *b_intraNukeNParticles;   //!
   TBranch        *b_intraNukeNeutronQuasiElasticScatter;   //!
   TBranch        *b_intraNukeOtherProcess;   //!
   TBranch        *b_is_true_muon;   //!
   TBranch        *b_muon_enters_front;   //!
   TBranch        *b_n_odClusters;   //!
   TBranch        *b_n_odClustersWithTimeCut;   //!
   TBranch        *b_ntracks;   //!
   TBranch        *b_part_response_recoil_nClus;   //!
   TBranch        *b_passVertexZCut;   //!
   TBranch        *b_phys_energy_in_road_downstream_nplanes;   //!
   TBranch        *b_phys_energy_in_road_upstream_nplanes;   //!
   TBranch        *b_phys_n_dead_discr_pair;   //!
   TBranch        *b_phys_n_dead_discr_pair_in_prim_track_region;   //!
   TBranch        *b_phys_n_dead_discr_pair_two_mod_downstream_prim_track;   //!
   TBranch        *b_phys_n_dead_discr_pair_two_mod_upstream_prim_vtx;   //!
   TBranch        *b_phys_n_dead_discr_pair_upstream_prim_track_proj;   //!
   TBranch        *b_phys_vertex_is_fiducial;   //!
   TBranch        *b_reco_muon_topology;   //!
   TBranch        *b_timeSlice;   //!
   TBranch        *b_vtx_fit_converged;   //!
   TBranch        *b_endPointEnergy;   //!
   TBranch        *b_energy_from_mc;   //!
   TBranch        *b_energy_from_mc_fraction;   //!
   TBranch        *b_energy_from_mc_fraction_of_highest;   //!
   TBranch        *b_hadronic_energy;   //!
   TBranch        *b_intraNukeProtonMomentum;   //!
   TBranch        *b_isolatedEnergy;   //!
   TBranch        *b_isolatedEnergy_ecal;   //!
   TBranch        *b_isolatedEnergy_hcal;   //!
   TBranch        *b_isolatedEnergy_targets;   //!
   TBranch        *b_isolatedEnergy_tracker;   //!
   TBranch        *b_muonFuzzEnergy;   //!
   TBranch        *b_muon_phi;   //!
   TBranch        *b_muon_theta;   //!
   TBranch        *b_muon_thetaX;   //!
   TBranch        *b_muon_thetaY;   //!
   TBranch        *b_numi_horn_curr;   //!
   TBranch        *b_numi_pot;   //!
   TBranch        *b_numi_x;   //!
   TBranch        *b_numi_x_width;   //!
   TBranch        *b_numi_y;   //!
   TBranch        *b_numi_y_width;   //!
   TBranch        *b_odEnergy;   //!
   TBranch        *b_odEnergyWithTimeCut;   //!
   TBranch        *b_part_response_recoil_em_id;   //!
   TBranch        *b_part_response_recoil_em_id_err;   //!
   TBranch        *b_part_response_recoil_em_od;   //!
   TBranch        *b_part_response_recoil_em_od_err;   //!
   TBranch        *b_part_response_recoil_high_neutron_id;   //!
   TBranch        *b_part_response_recoil_high_neutron_id_err;   //!
   TBranch        *b_part_response_recoil_high_neutron_od;   //!
   TBranch        *b_part_response_recoil_high_neutron_od_err;   //!
   TBranch        *b_part_response_recoil_low_neutron_id;   //!
   TBranch        *b_part_response_recoil_low_neutron_id_err;   //!
   TBranch        *b_part_response_recoil_low_neutron_od;   //!
   TBranch        *b_part_response_recoil_low_neutron_od_err;   //!
   TBranch        *b_part_response_recoil_meson_id;   //!
   TBranch        *b_part_response_recoil_meson_id_err;   //!
   TBranch        *b_part_response_recoil_meson_od;   //!
   TBranch        *b_part_response_recoil_meson_od_err;   //!
   TBranch        *b_part_response_recoil_mid_neutron_id;   //!
   TBranch        *b_part_response_recoil_mid_neutron_id_err;   //!
   TBranch        *b_part_response_recoil_mid_neutron_od;   //!
   TBranch        *b_part_response_recoil_mid_neutron_od_err;   //!
   TBranch        *b_part_response_recoil_muon_id;   //!
   TBranch        *b_part_response_recoil_muon_id_err;   //!
   TBranch        *b_part_response_recoil_muon_od;   //!
   TBranch        *b_part_response_recoil_muon_od_err;   //!
   TBranch        *b_part_response_recoil_other_id;   //!
   TBranch        *b_part_response_recoil_other_id_err;   //!
   TBranch        *b_part_response_recoil_other_od;   //!
   TBranch        *b_part_response_recoil_other_od_err;   //!
   TBranch        *b_part_response_recoil_proton_id;   //!
   TBranch        *b_part_response_recoil_proton_id_err;   //!
   TBranch        *b_part_response_recoil_proton_od;   //!
   TBranch        *b_part_response_recoil_proton_od_err;   //!
   TBranch        *b_part_response_recoil_xtalk_id;   //!
   TBranch        *b_part_response_recoil_xtalk_id_err;   //!
   TBranch        *b_part_response_recoil_xtalk_od;   //!
   TBranch        *b_part_response_recoil_xtalk_od_err;   //!
   TBranch        *b_phys_energy_dispersed;   //!
   TBranch        *b_phys_energy_in_road_downstream;   //!
   TBranch        *b_phys_energy_in_road_upstream;   //!
   TBranch        *b_phys_energy_unattached;   //!
   TBranch        *b_prim_vtx_smallest_opening_angle;   //!
   TBranch        *b_primaryVertexEnergy;   //!
   TBranch        *b_protonFuzzEnergy;   //!
   TBranch        *b_secondaryVertexEnergy;   //!
   TBranch        *b_vtx_fit_chi2;   //!
   TBranch        *b_vtxprong_energy_cal;   //!
   TBranch        *b_vtxprong_energy_visible;   //!
   TBranch        *b_vtxprong_energy_visible_isoblobs;   //!
   TBranch        *b_vtxprong_energy_visible_vtxblobs;   //!
   TBranch        *b_has_michel_category_sz;   //!
   TBranch        *b_has_michel_category;   //!
   TBranch        *b_has_michel_in_vertex_point_sz;   //!
   TBranch        *b_has_michel_in_vertex_point;   //!
   TBranch        *b_has_michel_ndigits_sz;   //!
   TBranch        *b_has_michel_ndigits;   //!
   TBranch        *b_has_michel_vertex_type_sz;   //!
   TBranch        *b_has_michel_vertex_type;   //!
   TBranch        *b_proton_enters_front;   //!
   TBranch        *b_n_vtxprong_isoblobs;   //!
   TBranch        *b_vtxprong_isoblob_nclusters;   //!
   TBranch        *b_n_vtxprong_vtxblobs;   //!
   TBranch        *b_vtxprong_vtxblob_nclusters;   //!
   TBranch        *b_clusterU_Angle_sz;   //!
   TBranch        *b_clusterU_Angle;   //!
   TBranch        *b_clusterU_Radius_sz;   //!
   TBranch        *b_clusterU_Radius;   //!
   TBranch        *b_clusterU_timeDiff_sz;   //!
   TBranch        *b_clusterU_timeDiff;   //!
   TBranch        *b_clusterU_viewDist_sz;   //!
   TBranch        *b_clusterU_viewDist;   //!
   TBranch        *b_clusterU_visE_binned_sz;   //!
   TBranch        *b_clusterU_visE_binned;   //!
   TBranch        *b_clusterU_visEnergy_sz;   //!
   TBranch        *b_clusterU_visEnergy;   //!
   TBranch        *b_clusterU_zDist_sz;   //!
   TBranch        *b_clusterU_zDist;   //!
   TBranch        *b_clusterV_Angle_sz;   //!
   TBranch        *b_clusterV_Angle;   //!
   TBranch        *b_clusterV_Radius_sz;   //!
   TBranch        *b_clusterV_Radius;   //!
   TBranch        *b_clusterV_timeDiff_sz;   //!
   TBranch        *b_clusterV_timeDiff;   //!
   TBranch        *b_clusterV_viewDist_sz;   //!
   TBranch        *b_clusterV_viewDist;   //!
   TBranch        *b_clusterV_visE_binned_sz;   //!
   TBranch        *b_clusterV_visE_binned;   //!
   TBranch        *b_clusterV_visEnergy_sz;   //!
   TBranch        *b_clusterV_visEnergy;   //!
   TBranch        *b_clusterV_zDist_sz;   //!
   TBranch        *b_clusterV_zDist;   //!
   TBranch        *b_clusterX_Angle_sz;   //!
   TBranch        *b_clusterX_Angle;   //!
   TBranch        *b_clusterX_Radius_sz;   //!
   TBranch        *b_clusterX_Radius;   //!
   TBranch        *b_clusterX_timeDiff_sz;   //!
   TBranch        *b_clusterX_timeDiff;   //!
   TBranch        *b_clusterX_viewDist_sz;   //!
   TBranch        *b_clusterX_viewDist;   //!
   TBranch        *b_clusterX_visE_binned_sz;   //!
   TBranch        *b_clusterX_visE_binned;   //!
   TBranch        *b_clusterX_visEnergy_sz;   //!
   TBranch        *b_clusterX_visEnergy;   //!
   TBranch        *b_clusterX_zDist_sz;   //!
   TBranch        *b_clusterX_zDist;   //!
   TBranch        *b_fit_vtx;   //!
   TBranch        *b_has_michel_distance_sz;   //!
   TBranch        *b_has_michel_distance;   //!
   TBranch        *b_has_michel_energy_sz;   //!
   TBranch        *b_has_michel_energy;   //!
   TBranch        *b_has_michel_time_diff_sz;   //!
   TBranch        *b_has_michel_time_diff;   //!
   TBranch        *b_intraNukeProtonMomentumVec;   //!
   TBranch        *b_reco_muon_vertex;   //!
   TBranch        *b_vtxprong_isoblob_visenergy;   //!
   TBranch        *b_vtxprong_vtxblob_visenergy;   //!
   TBranch        *b_truth_has_physics_event;   //!
   TBranch        *b_truth_has_michel_electron;   //!
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
   TBranch        *b_truth_has_michel_from_pion_minus_momentum_sz;   //!
   TBranch        *b_truth_has_michel_from_pion_minus_momentum;   //!
   TBranch        *b_truth_has_michel_from_pion_plus_momentum_sz;   //!
   TBranch        *b_truth_has_michel_from_pion_plus_momentum;   //!
   TBranch        *b_NukeCCQETwoTrack_nuFlavor;   //!
   TBranch        *b_NukeCCQETwoTrack_nuHelicity;   //!
   TBranch        *b_NukeCCQETwoTrack_intCurrent;   //!
   TBranch        *b_NukeCCQETwoTrack_intType;   //!
   TBranch        *b_NukeCCQETwoTrack_E;   //!
   TBranch        *b_NukeCCQETwoTrack_Q2;   //!
   TBranch        *b_NukeCCQETwoTrack_x;   //!
   TBranch        *b_NukeCCQETwoTrack_y;   //!
   TBranch        *b_NukeCCQETwoTrack_W;   //!
   TBranch        *b_NukeCCQETwoTrack_score;   //!
   TBranch        *b_NukeCCQETwoTrack_leptonE;   //!
   TBranch        *b_NukeCCQETwoTrack_vtx;   //!
   TBranch        *b_NukeCCQETwoTrack_minos_trk_is_contained;   //!
   TBranch        *b_NukeCCQETwoTrack_minos_trk_is_ok;   //!
   TBranch        *b_NukeCCQETwoTrack_minos_used_range;   //!
   TBranch        *b_NukeCCQETwoTrack_minos_used_curvature;   //!
   TBranch        *b_NukeCCQETwoTrack_isMuonInsideOD;   //!
   TBranch        *b_NukeCCQETwoTrack_minos_trk_end_plane;   //!
   TBranch        *b_NukeCCQETwoTrack_minos_trk_quality;   //!
   TBranch        *b_NukeCCQETwoTrack_muon_down_hcal;   //!
   TBranch        *b_NukeCCQETwoTrack_muon_minos_stub;   //!
   TBranch        *b_NukeCCQETwoTrack_muon_minos_track;   //!
   TBranch        *b_NukeCCQETwoTrack_muon_odLastFrame;   //!
   TBranch        *b_NukeCCQETwoTrack_muon_odLastStory;   //!
   TBranch        *b_NukeCCQETwoTrack_muon_od_track;   //!
   TBranch        *b_NukeCCQETwoTrack_muon_side_ecal;   //!
   TBranch        *b_NukeCCQETwoTrack_muon_trk_pat_history;   //!
   TBranch        *b_NukeCCQETwoTrack_ntrajMuonProng;   //!
   TBranch        *b_NukeCCQETwoTrack_r_minos_trk_vtx_plane;   //!
   TBranch        *b_NukeCCQETwoTrack_t_minos_trk_numFSMuons;   //!
   TBranch        *b_NukeCCQETwoTrack_t_minos_trk_primFSLeptonPDG;   //!
   TBranch        *b_NukeCCQETwoTrack_targetID;   //!
   TBranch        *b_NukeCCQETwoTrack_targetZ;   //!
   TBranch        *b_NukeCCQETwoTrack_trajMuonProngPDG;   //!
   TBranch        *b_NukeCCQETwoTrack_trajMuonProngPrimary;   //!
   TBranch        *b_NukeCCQETwoTrack_endMuonTrajMomentum;   //!
   TBranch        *b_NukeCCQETwoTrack_endMuonTrajXPosition;   //!
   TBranch        *b_NukeCCQETwoTrack_endMuonTrajYPosition;   //!
   TBranch        *b_NukeCCQETwoTrack_endMuonTrajZPosition;   //!
   TBranch        *b_NukeCCQETwoTrack_minos_trk_bave;   //!
   TBranch        *b_NukeCCQETwoTrack_minos_trk_chi2;   //!
   TBranch        *b_NukeCCQETwoTrack_minos_trk_end_u;   //!
   TBranch        *b_NukeCCQETwoTrack_minos_trk_end_v;   //!
   TBranch        *b_NukeCCQETwoTrack_minos_trk_end_x;   //!
   TBranch        *b_NukeCCQETwoTrack_minos_trk_end_y;   //!
   TBranch        *b_NukeCCQETwoTrack_minos_trk_end_z;   //!
   TBranch        *b_NukeCCQETwoTrack_minos_trk_eqp;   //!
   TBranch        *b_NukeCCQETwoTrack_minos_trk_eqp_qp;   //!
   TBranch        *b_NukeCCQETwoTrack_minos_trk_fit_pass;   //!
   TBranch        *b_NukeCCQETwoTrack_minos_trk_ndf;   //!
   TBranch        *b_NukeCCQETwoTrack_minos_trk_p;   //!
   TBranch        *b_NukeCCQETwoTrack_minos_trk_p_curvature;   //!
   TBranch        *b_NukeCCQETwoTrack_minos_trk_p_range;   //!
   TBranch        *b_NukeCCQETwoTrack_minos_trk_qp;   //!
   TBranch        *b_NukeCCQETwoTrack_minos_trk_vtx_x;   //!
   TBranch        *b_NukeCCQETwoTrack_minos_trk_vtx_y;   //!
   TBranch        *b_NukeCCQETwoTrack_minos_trk_vtx_z;   //!
   TBranch        *b_NukeCCQETwoTrack_muon_enu;   //!
   TBranch        *b_NukeCCQETwoTrack_muon_odClustersAvgTime;   //!
   TBranch        *b_NukeCCQETwoTrack_muon_odElossMomentum;   //!
   TBranch        *b_NukeCCQETwoTrack_muon_odEndX;   //!
   TBranch        *b_NukeCCQETwoTrack_muon_odEndY;   //!
   TBranch        *b_NukeCCQETwoTrack_muon_odEndZ;   //!
   TBranch        *b_NukeCCQETwoTrack_muon_odFaceX;   //!
   TBranch        *b_NukeCCQETwoTrack_muon_odFaceY;   //!
   TBranch        *b_NukeCCQETwoTrack_muon_odFaceZ;   //!
   TBranch        *b_NukeCCQETwoTrack_muon_odLastClusZ;   //!
   TBranch        *b_NukeCCQETwoTrack_muon_odStopDistMomentum;   //!
   TBranch        *b_NukeCCQETwoTrack_muon_odTrackAvgTime;   //!
   TBranch        *b_NukeCCQETwoTrack_muon_phi;   //!
   TBranch        *b_NukeCCQETwoTrack_muon_q2;   //!
   TBranch        *b_NukeCCQETwoTrack_muon_score;   //!
   TBranch        *b_NukeCCQETwoTrack_muon_theta;   //!
   TBranch        *b_NukeCCQETwoTrack_muon_thetaX;   //!
   TBranch        *b_NukeCCQETwoTrack_muon_thetaY;   //!
   TBranch        *b_NukeCCQETwoTrack_r_minos_trk_bdL;   //!
   TBranch        *b_NukeCCQETwoTrack_r_minos_trk_end_dcosx;   //!
   TBranch        *b_NukeCCQETwoTrack_r_minos_trk_end_dcosy;   //!
   TBranch        *b_NukeCCQETwoTrack_r_minos_trk_end_dcosz;   //!
   TBranch        *b_NukeCCQETwoTrack_r_minos_trk_vtx_dcosx;   //!
   TBranch        *b_NukeCCQETwoTrack_r_minos_trk_vtx_dcosy;   //!
   TBranch        *b_NukeCCQETwoTrack_r_minos_trk_vtx_dcosz;   //!
   TBranch        *b_NukeCCQETwoTrack_t_minos_trk_primFSLepMinosInitProjPx;   //!
   TBranch        *b_NukeCCQETwoTrack_t_minos_trk_primFSLepMinosInitProjPy;   //!
   TBranch        *b_NukeCCQETwoTrack_t_minos_trk_primFSLepMinosInitProjPz;   //!
   TBranch        *b_NukeCCQETwoTrack_t_minos_trk_primFSLepMinosInitProjX;   //!
   TBranch        *b_NukeCCQETwoTrack_t_minos_trk_primFSLepMinosInitProjY;   //!
   TBranch        *b_NukeCCQETwoTrack_t_minos_trk_primFSLepMinosInitProjZ;   //!
   TBranch        *b_NukeCCQETwoTrack_t_minos_trk_primFSLepMnvFinalPx;   //!
   TBranch        *b_NukeCCQETwoTrack_t_minos_trk_primFSLepMnvFinalPy;   //!
   TBranch        *b_NukeCCQETwoTrack_t_minos_trk_primFSLepMnvFinalPz;   //!
   TBranch        *b_NukeCCQETwoTrack_t_minos_trk_primFSLepMnvFinalX;   //!
   TBranch        *b_NukeCCQETwoTrack_t_minos_trk_primFSLepMnvFinalY;   //!
   TBranch        *b_NukeCCQETwoTrack_t_minos_trk_primFSLepMnvFinalZ;   //!
   TBranch        *b_NukeCCQETwoTrack_t_minos_trk_primFSLepMnvInitPx;   //!
   TBranch        *b_NukeCCQETwoTrack_t_minos_trk_primFSLepMnvInitPy;   //!
   TBranch        *b_NukeCCQETwoTrack_t_minos_trk_primFSLepMnvInitPz;   //!
   TBranch        *b_NukeCCQETwoTrack_t_minos_trk_primFSLepMnvInitX;   //!
   TBranch        *b_NukeCCQETwoTrack_t_minos_trk_primFSLepMnvInitY;   //!
   TBranch        *b_NukeCCQETwoTrack_t_minos_trk_primFSLepMnvInitZ;   //!
   TBranch        *b_NukeCCQETwoTrack_targetZPos;   //!
   TBranch        *b_NukeCCQETwoTrack_trajMuonPhi;   //!
   TBranch        *b_NukeCCQETwoTrack_trajMuonProngMomentum;   //!
   TBranch        *b_NukeCCQETwoTrack_trajMuonTheta;   //!
   TBranch        *b_NukeCCQETwoTrack_isProtonInsideOD;   //!
   TBranch        *b_NukeCCQETwoTrack_ntrajProngProng;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_kinked;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_odMatch;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_trk_pat_history;   //!
   TBranch        *b_NukeCCQETwoTrack_trajProtonProngPDG;   //!
   TBranch        *b_NukeCCQETwoTrack_trajProtonProngPrimary;   //!
   TBranch        *b_NukeCCQETwoTrack_coplanarAngle;   //!
   TBranch        *b_NukeCCQETwoTrack_endProtonTrajMomentum;   //!
   TBranch        *b_NukeCCQETwoTrack_endProtonTrajXPosition;   //!
   TBranch        *b_NukeCCQETwoTrack_endProtonTrajYPosition;   //!
   TBranch        *b_NukeCCQETwoTrack_endProtonTrajZPosition;   //!
   TBranch        *b_NukeCCQETwoTrack_muon_endPoint;   //!
   TBranch        *b_NukeCCQETwoTrack_muon_startPoint;   //!
   TBranch        *b_NukeCCQETwoTrack_open_angle;   //!
   TBranch        *b_NukeCCQETwoTrack_pion_chi2_ndf;   //!
   TBranch        *b_NukeCCQETwoTrack_pion_score;   //!
   TBranch        *b_NukeCCQETwoTrack_pion_score1;   //!
   TBranch        *b_NukeCCQETwoTrack_pion_score2;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_E;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_chi2_ndf;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_ekin;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_endPointX;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_endPointY;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_endPointZ;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_enu;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_momentum_shift_BetheBloch_Down;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_momentum_shift_BetheBloch_Up;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_momentum_shift_Birks;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_momentum_shift_MEU_Down;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_momentum_shift_MEU_Up;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_momentum_shift_Mass_Down;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_momentum_shift_Mass_Up;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_p;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_p_calCorrection;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_p_dEdXTool;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_p_visEnergy;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_phi;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_px;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_py;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_pz;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_q2;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_score;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_score1;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_score1_shift_BetheBloch_Down;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_score1_shift_BetheBloch_Up;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_score1_shift_Birks;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_score1_shift_MEU_Down;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_score1_shift_MEU_Up;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_score1_shift_Mass_Down;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_score1_shift_Mass_Up;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_score2;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_startPointX;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_startPointY;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_startPointZ;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_theta;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_thetaX;   //!
   TBranch        *b_NukeCCQETwoTrack_proton_thetaY;   //!
   TBranch        *b_NukeCCQETwoTrack_sys_muon_curve_energy_shift;   //!
   TBranch        *b_NukeCCQETwoTrack_sys_muon_energy_shift;   //!
   TBranch        *b_NukeCCQETwoTrack_sys_muon_minerva_energy_shift;   //!
   TBranch        *b_NukeCCQETwoTrack_sys_muon_qSquared_shift;   //!
   TBranch        *b_NukeCCQETwoTrack_sys_muon_range_energy_shift;   //!
   TBranch        *b_NukeCCQETwoTrack_sys_muon_wSquared_shift;   //!
   TBranch        *b_NukeCCQETwoTrack_sys_muon_xbj_shift;   //!
   TBranch        *b_NukeCCQETwoTrack_sys_muon_y_shift;   //!
   TBranch        *b_NukeCCQETwoTrack_sys_nu_energy_shift;   //!
   TBranch        *b_NukeCCQETwoTrack_sys_recoil_energy_shift;   //!
   TBranch        *b_NukeCCQETwoTrack_sys_recoil_qSquared_shift;   //!
   TBranch        *b_NukeCCQETwoTrack_sys_recoil_wSquared_shift;   //!
   TBranch        *b_NukeCCQETwoTrack_sys_recoil_xbj_shift;   //!
   TBranch        *b_NukeCCQETwoTrack_sys_recoil_y_shift;   //!
   TBranch        *b_NukeCCQETwoTrack_sys_total_qSquared_shift;   //!
   TBranch        *b_NukeCCQETwoTrack_sys_total_wSquared_shift;   //!
   TBranch        *b_NukeCCQETwoTrack_sys_total_xbj_shift;   //!
   TBranch        *b_NukeCCQETwoTrack_sys_total_y_shift;   //!
   TBranch        *b_NukeCCQETwoTrack_trajProtonPhi;   //!
   TBranch        *b_NukeCCQETwoTrack_trajProtonProngMomentum;   //!
   TBranch        *b_NukeCCQETwoTrack_trajProtonTheta;   //!
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
   TBranch        *b_mc_wgt_Norm_sz;   //!
   TBranch        *b_mc_wgt_Norm;   //!
   TBranch        *b_mc_wgt_ppfx_MIPPNumiYields_sz;   //!
   TBranch        *b_mc_wgt_ppfx_MIPPNumiYields;   //!
   TBranch        *b_mc_wgt_ppfx_TargetAttenuation_sz;   //!
   TBranch        *b_mc_wgt_ppfx_TargetAttenuation;   //!
   TBranch        *b_mc_wgt_ppfx_NA49_sz;   //!
   TBranch        *b_mc_wgt_ppfx_NA49;   //!
   TBranch        *b_mc_wgt_ppfx_MIPPKaonsYields_sz;   //!
   TBranch        *b_mc_wgt_ppfx_MIPPKaonsYields;   //!
   TBranch        *b_mc_wgt_ppfx_MIPPThinTarget_sz;   //!
   TBranch        *b_mc_wgt_ppfx_MIPPThinTarget;   //!
   TBranch        *b_mc_wgt_ppfx_Absorption_sz;   //!
   TBranch        *b_mc_wgt_ppfx_Absorption;   //!
   TBranch        *b_mc_wgt_ppfx_Others_sz;   //!
   TBranch        *b_mc_wgt_ppfx_Others;   //!
   TBranch        *b_n_prongs;   //!
   TBranch        *b_prong_nParticles;   //!
   TBranch        *b_prong_part_score;   //!
   TBranch        *b_prong_part_mass;   //!
   TBranch        *b_prong_part_charge;   //!
   TBranch        *b_prong_part_pid;   //!
   TBranch        *b_prong_part_E;   //!
   TBranch        *b_prong_part_pos;   //!

   NukeCCQETwoTrack_default(TTree *tree=0);
   virtual ~NukeCCQETwoTrack_default();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef NukeCCQETwoTrack_default_cxx
NukeCCQETwoTrack_default::NukeCCQETwoTrack_default(TTree *tree) : fChain(0) 
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
      f->GetObject("NukeCCQETwoTrack",tree);

#else // SINGLE_TREE

      // The following code should be used if you want this class to access a chain
      // of trees.
      TChain * chain = new TChain("NukeCCQETwoTrack","");
      chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_modified/run/grid/central_value/minerva/ana/v10r8/00/01/02/all_10200.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0005-0006-0007-0008_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0009-0010-0011-0012-0013_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0014-0015-0016-0017-0018_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0019-0020-0021-0022-0023_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0024-0025-0026-0027-0028_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0029-0030-0031-0032-0034_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0035-0036-0037-0038-0039_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0040-0041-0042-0043-0044_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0045-0046-0047-0048-0049_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0050-0051-0052-0053-0054_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0055-0056-0057-0058-0059_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0060-0061-0062-0063-0064_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0065-0066-0067-0068-0069_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0070-0071-0072-0073-0074_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0075-0076-0077-0078-0079_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0080-0081-0082-0083-0084_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0085-0086-0087-0088-0089_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0090-0091-0092-0093-0094_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0095-0096-0097-0098-0099_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0100-0101-0102-0103-0104_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0105-0106-0107-0108-0109_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0110-0111-0112-0113-0114_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0115-0116-0117-0118-0119_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0120-0121-0122-0123-0124_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0125-0126-0127-0128-0129_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0130-0131-0132-0133-0134_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0135-0136-0137-0138-0139_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0140-0141-0142-0143-0144_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0145-0146-0147-0148-0149_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0150-0151-0152-0153-0154_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0155-0156-0157-0158-0159_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0160-0161-0163-0164-0165_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0166-0167-0168-0169-0170_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0171-0172-0173-0174-0175_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0176-0177-0178-0179-0180_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0181-0182-0183-0184-0185_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0186-0187-0188-0189-0190_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0191-0192-0193-0194-0195_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/00/SIM_minerva_00010200_Subruns_0196-0197-0198-0199-0200_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0001-0002-0003-0004_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0005-0006-0007-0008_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0009-0010-0011-0012-0013_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0015-0016-0017-0018-0019_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0020-0021-0022-0023-0024_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0025-0026-0027-0028-0029_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0030-0031-0032-0033-0034_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0035-0036-0037-0038-0039_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0040-0041-0042-0043-0044_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0045-0046-0047-0048-0049_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0050-0051-0052-0053-0054_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0055-0056-0057-0058-0059_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0060-0061-0062-0063-0064_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0065-0066-0067-0068-0069_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0070-0071-0072-0073-0074_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0075-0076-0077-0078-0079_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0080-0081-0082-0083-0084_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0085-0086-0087-0088-0089_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0090-0091-0092-0093-0094_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0095-0096-0097-0098-0099_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0100-0101-0102-0103-0104_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0105-0106-0107-0108-0109_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0110-0111-0112-0113-0114_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0116-0117-0118-0119-0120_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0121-0122-0123-0124-0125_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0126-0127-0128-0129-0130_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0131-0132-0133-0134-0135_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0136-0137-0138-0139-0140_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0141-0142-0143-0144-0145_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0146-0147-0148-0149-0150_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0151-0152-0153-0154-0155_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0156-0157-0158-0159-0160_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0161-0162-0163-0164-0165_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0166-0167-0168-0169-0170_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0171-0172-0173-0174-0175_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0176-0177-0178-0179-0180_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0181-0182-0183-0184-0185_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0186-0187-0188-0189-0190_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0191-0192-0193-0194-0195_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/02/SIM_minerva_00010202_Subruns_0196-0197-0198-0199-0200_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0001-0002-0003-0004_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0005-0006-0007-0008_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0009-0010-0011-0012-0013_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0014-0015-0016-0017-0018_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0019-0020-0021-0022-0023_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0024-0025-0026-0027-0028_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0029-0030-0031-0032-0033_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0034-0035-0036-0037-0038_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0039-0040-0041-0042-0043_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0044-0045-0046-0047-0048_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0049-0050-0051-0052-0053_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0054-0055-0056-0057-0058_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0059-0060-0061-0062-0063_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0064-0065-0066-0067-0068_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0069-0070-0071-0072-0073_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0075-0076-0077-0078-0079_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0080-0081-0082-0083-0084_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0085-0086-0087-0088-0089_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0090-0091-0092-0093-0094_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0095-0096-0097-0098-0099_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0100-0101-0102-0103-0104_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0105-0106-0107-0108-0109_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0110-0111-0112-0113-0114_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0115-0116-0117-0118-0119_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0120-0121-0122-0123-0124_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0125-0126-0127-0128-0129_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0130-0131-0132-0133-0134_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0135-0136-0137-0138-0139_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0140-0141-0142-0143-0144_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0145-0146-0147-0148-0149_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0150-0151-0152-0153-0154_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0155-0156-0157-0158-0159_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0160-0161-0162-0163-0164_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0165-0166-0167-0168-0169_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0170-0171-0172-0173-0174_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0175-0176-0177-0178-0179_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0180-0181-0182-0183-0184_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0185-0186-0187-0188-0189_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0190-0191-0192-0193-0194_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/04/SIM_minerva_00010204_Subruns_0195-0197-0198-0199-0200_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0001-0002-0003-0004_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0005-0006-0007-0008-0009_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0010-0011-0012-0013-0014_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0015-0016-0017-0018-0019_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0020-0021-0022-0023-0024_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0025-0026-0027-0028-0029_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0030-0031-0032-0033-0034_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0035-0036-0037-0038-0039_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0040-0041-0042-0043-0044_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0045-0046-0047-0048-0049_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0050-0051-0052-0053-0054_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0055-0056-0057-0058-0059_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0060-0061-0062-0063-0064_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0065-0066-0067-0068-0069_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0070-0071-0072-0073-0074_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0075-0076-0077-0078-0079_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0080-0081-0082-0083-0084_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0085-0086-0087-0088-0089_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0090-0091-0092-0093-0094_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0095-0096-0097-0098-0099_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0100-0101-0102-0103-0104_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0105-0106-0107-0108-0109_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0110-0111-0112-0113-0114_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0115-0116-0117-0118-0119_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0120-0121-0122-0123-0124_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0125-0126-0127-0128-0129_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0130-0131-0132-0133-0134_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0135-0136-0137-0138-0139_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0140-0141-0142-0143-0144_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0145-0146-0147-0148-0149_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0150-0151-0152-0154-0155_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0156-0157-0158-0159-0160_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0161-0162-0163-0164-0165_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0166-0167-0168-0169-0170_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0171-0172-0173-0174-0175_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0176-0177-0178-0179-0180_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0181-0182-0183-0184-0185_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0186-0187-0188-0189-0190_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0191-0192-0193-0194-0195_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/06/SIM_minerva_00010206_Subruns_0196-0197-0198-0199-0200_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0001-0002-0003-0004-0005_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0006-0007-0008-0009-0010_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0011-0012-0013-0014-0015_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0016-0017-0018-0019-0020_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0021-0022-0023-0024-0025_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0026-0027-0028-0029-0030_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0031-0032-0033-0034-0035_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0036-0037-0038-0039-0040_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0041-0042-0043-0044-0045_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0046-0047-0048-0049-0050_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0051-0052-0053-0054-0055_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0056-0057-0058-0059-0060_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0061-0062-0063-0064-0065_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0066-0067-0068-0069-0070_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0071-0072-0073-0074-0075_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0076-0077-0078-0079-0080_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0081-0082-0083-0084-0085_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0086-0087-0088-0089-0090_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0091-0092-0093-0094-0095_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0096-0097-0098-0099-0100_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0101-0102-0103-0104-0105_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0106-0107-0108-0109-0110_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0111-0112-0113-0114-0115_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0116-0117-0118-0119-0120_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0121-0122-0123-0124-0125_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0126-0127-0128-0129-0130_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0131-0132-0133-0134-0135_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0136-0137-0138-0139-0140_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0141-0142-0143-0144-0145_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0146-0147-0148-0149-0150_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0151-0152-0153-0154-0155_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0156-0157-0158-0159-0160_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0161-0162-0163-0164-0165_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0166-0167-0168-0169-0170_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0171-0172-0173-0174-0175_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0176-0177-0178-0179-0180_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0181-0182-0183-0184-0185_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0186-0187-0188-0189-0190_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0191-0192-0193-0194-0195_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/08/SIM_minerva_00010208_Subruns_0196-0197-0198-0199-0200_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0001-0002-0003-0004_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0005-0006-0007-0008_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0009-0010-0011-0012-0013_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0014-0015-0016-0017-0018_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0019-0020-0021-0022-0023_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0025-0026-0027-0028-0030_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0032-0033-0034-0035-0036_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0037-0038-0039-0040-0041_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0042-0043-0044-0045-0046_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0047-0048-0049-0050-0051_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0052-0053-0054-0055-0056_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0057-0058-0059-0060-0061_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0062-0063-0064-0065-0066_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0067-0068-0069-0070-0071_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0072-0074-0075-0076-0077_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0078-0079-0080-0081-0082_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0083-0084-0085-0086-0087_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0088-0089-0090-0091-0092_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0093-0094-0095-0096-0097_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0098-0099-0100-0101-0102_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0103-0104-0105-0106-0107_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0108-0109-0110-0111-0112_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0113-0114-0115-0116-0117_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0118-0119-0120-0122-0123_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0124-0125-0126-0127-0128_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0129-0130-0131-0132-0133_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0134-0135-0136-0137-0138_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0139-0140-0141-0142-0143_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0144-0145-0147-0148-0149_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0150-0151-0152-0153-0154_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0155-0156-0157-0158-0159_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0160-0161-0162-0163-0164_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0165-0166-0167-0168-0169_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0170-0171-0172-0173-0174_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0175-0176-0177-0178-0179_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0180-0181-0182-0183-0184_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0185-0186-0187-0188-0189_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0191-0192-0193-0194-0195_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/10/SIM_minerva_00010210_Subruns_0196-0197-0198-0199-0200_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0001-0002-0003-0004_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0005-0006-0007-0008_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0009-0010-0011-0012_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0013-0014-0015-0016_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0017-0018-0019-0020-0021_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0022-0023-0024-0025-0026_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0027-0028-0029-0030-0031_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0032-0033-0034-0035-0036_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0037-0038-0039-0040-0041_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0042-0043-0044-0045-0046_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0047-0048-0049-0050-0051_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0052-0053-0054-0055-0056_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0057-0058-0059-0060-0061_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0062-0063-0064-0065-0066_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0067-0068-0069-0070-0071_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0072-0073-0074-0075-0076_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0077-0078-0079-0080-0081_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0083-0084-0085-0086-0087_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0088-0089-0090-0091-0092_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0093-0094-0095-0097-0098_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0099-0100-0101-0102-0103_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0104-0105-0106-0107-0108_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0109-0110-0111-0112-0113_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0114-0115-0116-0117-0118_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0119-0120-0121-0122-0123_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0124-0125-0126-0127-0128_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0129-0130-0131-0132-0133_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0134-0135-0136-0138-0139_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0140-0141-0142-0143-0144_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0145-0146-0147-0148-0149_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0150-0151-0152-0153-0154_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0155-0156-0157-0158-0159_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0160-0161-0162-0163-0164_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0165-0166-0167-0168-0169_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0170-0171-0172-0173-0174_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0175-0176-0177-0178-0179_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0180-0181-0182-0183-0184_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0185-0186-0187-0188-0189_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0190-0191-0192-0193-0194_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/12/SIM_minerva_00010212_Subruns_0195-0196-0197-0198-0200_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0001-0002-0003-0004-0005_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0006-0007-0008-0009-0010_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0011-0012-0013-0014-0015_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0016-0017-0018-0019-0020_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0021-0022-0023-0024-0025_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0026-0027-0028-0029-0030_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0031-0032-0033-0034-0035_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0036-0037-0038-0039-0040_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0041-0042-0043-0044-0045_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0046-0047-0048-0049-0050_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0051-0052-0053-0054-0055_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0056-0057-0058-0059-0060_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0061-0062-0063-0064-0065_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0067-0068-0069-0070-0071_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0072-0073-0074-0075-0076_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0077-0078-0080-0081-0083_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0084-0086-0087-0088-0089_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0090-0091-0092-0093-0094_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0095-0096-0097-0098-0099_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0100-0101-0102-0103-0104_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0105-0106-0107-0108-0109_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0110-0111-0112-0113-0114_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0115-0116-0117-0118-0119_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0120-0121-0122-0123-0124_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0125-0126-0127-0128-0129_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0130-0131-0132-0133-0134_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0135-0136-0137-0138-0139_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0140-0141-0142-0143-0144_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0145-0146-0147-0148-0149_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0150-0151-0152-0153-0154_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0155-0156-0157-0158-0159_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0160-0161-0162-0164-0165_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0166-0167-0168-0169-0170_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0171-0172-0173-0174-0175_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0176-0177-0178-0179-0180_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0181-0182-0183-0184-0185_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0186-0187-0188-0189-0190_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0191-0192-0193-0194-0195_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      //chain->Add("/minerva/data/users/oaltinok/NukeCCQETwoTrack/dEdX_default/run/grid/central_value/minerva/ana/v10r8/00/01/02/14/SIM_minerva_00010214_Subruns_0196-0197-0198-0199-0200_NukeCCQETwoTrack_Ana_Tuple_v10r8.root/NukeCCQETwoTrack");
      tree = chain;
#endif // SINGLE_TREE

   }
   Init(tree);
}

NukeCCQETwoTrack_default::~NukeCCQETwoTrack_default()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t NukeCCQETwoTrack_default::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t NukeCCQETwoTrack_default::LoadTree(Long64_t entry)
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

void NukeCCQETwoTrack_default::Init(TTree *tree)
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
   fChain->SetBranchAddress("CreatedShortTracks", &CreatedShortTracks, &b_CreatedShortTracks);
   fChain->SetBranchAddress("FailContainedProng", &FailContainedProng, &b_FailContainedProng);
   fChain->SetBranchAddress("FailExitingProng", &FailExitingProng, &b_FailExitingProng);
   fChain->SetBranchAddress("FailFidVolume", &FailFidVolume, &b_FailFidVolume);
   fChain->SetBranchAddress("FailOutTracks", &FailOutTracks, &b_FailOutTracks);
   fChain->SetBranchAddress("FailRefitFidVolume", &FailRefitFidVolume, &b_FailRefitFidVolume);
   fChain->SetBranchAddress("FailShortOutTrack", &FailShortOutTrack, &b_FailShortOutTrack);
   fChain->SetBranchAddress("NoInteractionVertex", &NoInteractionVertex, &b_NoInteractionVertex);
   fChain->SetBranchAddress("NullVertex", &NullVertex, &b_NullVertex);
   fChain->SetBranchAddress("UnattachedProngsWithTracks", &UnattachedProngsWithTracks, &b_UnattachedProngsWithTracks);
   fChain->SetBranchAddress("broken_track_most_us_plane", &broken_track_most_us_plane, &b_broken_track_most_us_plane);
   fChain->SetBranchAddress("genie_n_charms", &genie_n_charms, &b_genie_n_charms);
   fChain->SetBranchAddress("genie_n_heavy_baryons", &genie_n_heavy_baryons, &b_genie_n_heavy_baryons);
   fChain->SetBranchAddress("genie_n_kaons", &genie_n_kaons, &b_genie_n_kaons);
   fChain->SetBranchAddress("genie_n_mesons", &genie_n_mesons, &b_genie_n_mesons);
   fChain->SetBranchAddress("genie_n_muons", &genie_n_muons, &b_genie_n_muons);
   fChain->SetBranchAddress("genie_n_neutrinos", &genie_n_neutrinos, &b_genie_n_neutrinos);
   fChain->SetBranchAddress("genie_n_neutrons", &genie_n_neutrons, &b_genie_n_neutrons);
   fChain->SetBranchAddress("genie_n_others", &genie_n_others, &b_genie_n_others);
   fChain->SetBranchAddress("genie_n_particles", &genie_n_particles, &b_genie_n_particles);
   fChain->SetBranchAddress("genie_n_photons", &genie_n_photons, &b_genie_n_photons);
   fChain->SetBranchAddress("genie_n_pi_zeros", &genie_n_pi_zeros, &b_genie_n_pi_zeros);
   fChain->SetBranchAddress("genie_n_pions", &genie_n_pions, &b_genie_n_pions);
   fChain->SetBranchAddress("genie_n_protons", &genie_n_protons, &b_genie_n_protons);
   fChain->SetBranchAddress("intraNukeDeltaPlusPlusDecay", &intraNukeDeltaPlusPlusDecay, &b_intraNukeDeltaPlusPlusDecay);
   fChain->SetBranchAddress("intraNukeNParticles", &intraNukeNParticles, &b_intraNukeNParticles);
   fChain->SetBranchAddress("intraNukeNeutronQuasiElasticScatter", &intraNukeNeutronQuasiElasticScatter, &b_intraNukeNeutronQuasiElasticScatter);
   fChain->SetBranchAddress("intraNukeOtherProcess", &intraNukeOtherProcess, &b_intraNukeOtherProcess);
   fChain->SetBranchAddress("is_true_muon", &is_true_muon, &b_is_true_muon);
   fChain->SetBranchAddress("muon_enters_front", &muon_enters_front, &b_muon_enters_front);
   fChain->SetBranchAddress("n_odClusters", &n_odClusters, &b_n_odClusters);
   fChain->SetBranchAddress("n_odClustersWithTimeCut", &n_odClustersWithTimeCut, &b_n_odClustersWithTimeCut);
   fChain->SetBranchAddress("ntracks", &ntracks, &b_ntracks);
   fChain->SetBranchAddress("part_response_recoil_nClus", &part_response_recoil_nClus, &b_part_response_recoil_nClus);
   fChain->SetBranchAddress("passVertexZCut", &passVertexZCut, &b_passVertexZCut);
   fChain->SetBranchAddress("phys_energy_in_road_downstream_nplanes", &phys_energy_in_road_downstream_nplanes, &b_phys_energy_in_road_downstream_nplanes);
   fChain->SetBranchAddress("phys_energy_in_road_upstream_nplanes", &phys_energy_in_road_upstream_nplanes, &b_phys_energy_in_road_upstream_nplanes);
   fChain->SetBranchAddress("phys_n_dead_discr_pair", &phys_n_dead_discr_pair, &b_phys_n_dead_discr_pair);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_in_prim_track_region", &phys_n_dead_discr_pair_in_prim_track_region, &b_phys_n_dead_discr_pair_in_prim_track_region);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_two_mod_downstream_prim_track", &phys_n_dead_discr_pair_two_mod_downstream_prim_track, &b_phys_n_dead_discr_pair_two_mod_downstream_prim_track);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_two_mod_upstream_prim_vtx", &phys_n_dead_discr_pair_two_mod_upstream_prim_vtx, &b_phys_n_dead_discr_pair_two_mod_upstream_prim_vtx);
   fChain->SetBranchAddress("phys_n_dead_discr_pair_upstream_prim_track_proj", &phys_n_dead_discr_pair_upstream_prim_track_proj, &b_phys_n_dead_discr_pair_upstream_prim_track_proj);
   fChain->SetBranchAddress("phys_vertex_is_fiducial", &phys_vertex_is_fiducial, &b_phys_vertex_is_fiducial);
   fChain->SetBranchAddress("reco_muon_topology", &reco_muon_topology, &b_reco_muon_topology);
   fChain->SetBranchAddress("timeSlice", &timeSlice, &b_timeSlice);
   fChain->SetBranchAddress("vtx_fit_converged", &vtx_fit_converged, &b_vtx_fit_converged);
   fChain->SetBranchAddress("endPointEnergy", &endPointEnergy, &b_endPointEnergy);
   fChain->SetBranchAddress("energy_from_mc", &energy_from_mc, &b_energy_from_mc);
   fChain->SetBranchAddress("energy_from_mc_fraction", &energy_from_mc_fraction, &b_energy_from_mc_fraction);
   fChain->SetBranchAddress("energy_from_mc_fraction_of_highest", &energy_from_mc_fraction_of_highest, &b_energy_from_mc_fraction_of_highest);
   fChain->SetBranchAddress("hadronic_energy", &hadronic_energy, &b_hadronic_energy);
   fChain->SetBranchAddress("intraNukeProtonMomentum", &intraNukeProtonMomentum, &b_intraNukeProtonMomentum);
   fChain->SetBranchAddress("isolatedEnergy", &isolatedEnergy, &b_isolatedEnergy);
   fChain->SetBranchAddress("isolatedEnergy_ecal", &isolatedEnergy_ecal, &b_isolatedEnergy_ecal);
   fChain->SetBranchAddress("isolatedEnergy_hcal", &isolatedEnergy_hcal, &b_isolatedEnergy_hcal);
   fChain->SetBranchAddress("isolatedEnergy_targets", &isolatedEnergy_targets, &b_isolatedEnergy_targets);
   fChain->SetBranchAddress("isolatedEnergy_tracker", &isolatedEnergy_tracker, &b_isolatedEnergy_tracker);
   fChain->SetBranchAddress("muonFuzzEnergy", &muonFuzzEnergy, &b_muonFuzzEnergy);
   fChain->SetBranchAddress("muon_phi", &muon_phi, &b_muon_phi);
   fChain->SetBranchAddress("muon_theta", &muon_theta, &b_muon_theta);
   fChain->SetBranchAddress("muon_thetaX", &muon_thetaX, &b_muon_thetaX);
   fChain->SetBranchAddress("muon_thetaY", &muon_thetaY, &b_muon_thetaY);
   fChain->SetBranchAddress("numi_horn_curr", &numi_horn_curr, &b_numi_horn_curr);
   fChain->SetBranchAddress("numi_pot", &numi_pot, &b_numi_pot);
   fChain->SetBranchAddress("numi_x", &numi_x, &b_numi_x);
   fChain->SetBranchAddress("numi_x_width", &numi_x_width, &b_numi_x_width);
   fChain->SetBranchAddress("numi_y", &numi_y, &b_numi_y);
   fChain->SetBranchAddress("numi_y_width", &numi_y_width, &b_numi_y_width);
   fChain->SetBranchAddress("odEnergy", &odEnergy, &b_odEnergy);
   fChain->SetBranchAddress("odEnergyWithTimeCut", &odEnergyWithTimeCut, &b_odEnergyWithTimeCut);
   fChain->SetBranchAddress("part_response_recoil_em_id", &part_response_recoil_em_id, &b_part_response_recoil_em_id);
   fChain->SetBranchAddress("part_response_recoil_em_id_err", &part_response_recoil_em_id_err, &b_part_response_recoil_em_id_err);
   fChain->SetBranchAddress("part_response_recoil_em_od", &part_response_recoil_em_od, &b_part_response_recoil_em_od);
   fChain->SetBranchAddress("part_response_recoil_em_od_err", &part_response_recoil_em_od_err, &b_part_response_recoil_em_od_err);
   fChain->SetBranchAddress("part_response_recoil_high_neutron_id", &part_response_recoil_high_neutron_id, &b_part_response_recoil_high_neutron_id);
   fChain->SetBranchAddress("part_response_recoil_high_neutron_id_err", &part_response_recoil_high_neutron_id_err, &b_part_response_recoil_high_neutron_id_err);
   fChain->SetBranchAddress("part_response_recoil_high_neutron_od", &part_response_recoil_high_neutron_od, &b_part_response_recoil_high_neutron_od);
   fChain->SetBranchAddress("part_response_recoil_high_neutron_od_err", &part_response_recoil_high_neutron_od_err, &b_part_response_recoil_high_neutron_od_err);
   fChain->SetBranchAddress("part_response_recoil_low_neutron_id", &part_response_recoil_low_neutron_id, &b_part_response_recoil_low_neutron_id);
   fChain->SetBranchAddress("part_response_recoil_low_neutron_id_err", &part_response_recoil_low_neutron_id_err, &b_part_response_recoil_low_neutron_id_err);
   fChain->SetBranchAddress("part_response_recoil_low_neutron_od", &part_response_recoil_low_neutron_od, &b_part_response_recoil_low_neutron_od);
   fChain->SetBranchAddress("part_response_recoil_low_neutron_od_err", &part_response_recoil_low_neutron_od_err, &b_part_response_recoil_low_neutron_od_err);
   fChain->SetBranchAddress("part_response_recoil_meson_id", &part_response_recoil_meson_id, &b_part_response_recoil_meson_id);
   fChain->SetBranchAddress("part_response_recoil_meson_id_err", &part_response_recoil_meson_id_err, &b_part_response_recoil_meson_id_err);
   fChain->SetBranchAddress("part_response_recoil_meson_od", &part_response_recoil_meson_od, &b_part_response_recoil_meson_od);
   fChain->SetBranchAddress("part_response_recoil_meson_od_err", &part_response_recoil_meson_od_err, &b_part_response_recoil_meson_od_err);
   fChain->SetBranchAddress("part_response_recoil_mid_neutron_id", &part_response_recoil_mid_neutron_id, &b_part_response_recoil_mid_neutron_id);
   fChain->SetBranchAddress("part_response_recoil_mid_neutron_id_err", &part_response_recoil_mid_neutron_id_err, &b_part_response_recoil_mid_neutron_id_err);
   fChain->SetBranchAddress("part_response_recoil_mid_neutron_od", &part_response_recoil_mid_neutron_od, &b_part_response_recoil_mid_neutron_od);
   fChain->SetBranchAddress("part_response_recoil_mid_neutron_od_err", &part_response_recoil_mid_neutron_od_err, &b_part_response_recoil_mid_neutron_od_err);
   fChain->SetBranchAddress("part_response_recoil_muon_id", &part_response_recoil_muon_id, &b_part_response_recoil_muon_id);
   fChain->SetBranchAddress("part_response_recoil_muon_id_err", &part_response_recoil_muon_id_err, &b_part_response_recoil_muon_id_err);
   fChain->SetBranchAddress("part_response_recoil_muon_od", &part_response_recoil_muon_od, &b_part_response_recoil_muon_od);
   fChain->SetBranchAddress("part_response_recoil_muon_od_err", &part_response_recoil_muon_od_err, &b_part_response_recoil_muon_od_err);
   fChain->SetBranchAddress("part_response_recoil_other_id", &part_response_recoil_other_id, &b_part_response_recoil_other_id);
   fChain->SetBranchAddress("part_response_recoil_other_id_err", &part_response_recoil_other_id_err, &b_part_response_recoil_other_id_err);
   fChain->SetBranchAddress("part_response_recoil_other_od", &part_response_recoil_other_od, &b_part_response_recoil_other_od);
   fChain->SetBranchAddress("part_response_recoil_other_od_err", &part_response_recoil_other_od_err, &b_part_response_recoil_other_od_err);
   fChain->SetBranchAddress("part_response_recoil_proton_id", &part_response_recoil_proton_id, &b_part_response_recoil_proton_id);
   fChain->SetBranchAddress("part_response_recoil_proton_id_err", &part_response_recoil_proton_id_err, &b_part_response_recoil_proton_id_err);
   fChain->SetBranchAddress("part_response_recoil_proton_od", &part_response_recoil_proton_od, &b_part_response_recoil_proton_od);
   fChain->SetBranchAddress("part_response_recoil_proton_od_err", &part_response_recoil_proton_od_err, &b_part_response_recoil_proton_od_err);
   fChain->SetBranchAddress("part_response_recoil_xtalk_id", &part_response_recoil_xtalk_id, &b_part_response_recoil_xtalk_id);
   fChain->SetBranchAddress("part_response_recoil_xtalk_id_err", &part_response_recoil_xtalk_id_err, &b_part_response_recoil_xtalk_id_err);
   fChain->SetBranchAddress("part_response_recoil_xtalk_od", &part_response_recoil_xtalk_od, &b_part_response_recoil_xtalk_od);
   fChain->SetBranchAddress("part_response_recoil_xtalk_od_err", &part_response_recoil_xtalk_od_err, &b_part_response_recoil_xtalk_od_err);
   fChain->SetBranchAddress("phys_energy_dispersed", &phys_energy_dispersed, &b_phys_energy_dispersed);
   fChain->SetBranchAddress("phys_energy_in_road_downstream", &phys_energy_in_road_downstream, &b_phys_energy_in_road_downstream);
   fChain->SetBranchAddress("phys_energy_in_road_upstream", &phys_energy_in_road_upstream, &b_phys_energy_in_road_upstream);
   fChain->SetBranchAddress("phys_energy_unattached", &phys_energy_unattached, &b_phys_energy_unattached);
   fChain->SetBranchAddress("prim_vtx_smallest_opening_angle", &prim_vtx_smallest_opening_angle, &b_prim_vtx_smallest_opening_angle);
   fChain->SetBranchAddress("primaryVertexEnergy", &primaryVertexEnergy, &b_primaryVertexEnergy);
   fChain->SetBranchAddress("protonFuzzEnergy", &protonFuzzEnergy, &b_protonFuzzEnergy);
   fChain->SetBranchAddress("secondaryVertexEnergy", &secondaryVertexEnergy, &b_secondaryVertexEnergy);
   fChain->SetBranchAddress("vtx_fit_chi2", &vtx_fit_chi2, &b_vtx_fit_chi2);
   fChain->SetBranchAddress("vtxprong_energy_cal", &vtxprong_energy_cal, &b_vtxprong_energy_cal);
   fChain->SetBranchAddress("vtxprong_energy_visible", &vtxprong_energy_visible, &b_vtxprong_energy_visible);
   fChain->SetBranchAddress("vtxprong_energy_visible_isoblobs", &vtxprong_energy_visible_isoblobs, &b_vtxprong_energy_visible_isoblobs);
   fChain->SetBranchAddress("vtxprong_energy_visible_vtxblobs", &vtxprong_energy_visible_vtxblobs, &b_vtxprong_energy_visible_vtxblobs);
   fChain->SetBranchAddress("has_michel_category_sz", &has_michel_category_sz, &b_has_michel_category_sz);
   fChain->SetBranchAddress("has_michel_category", has_michel_category, &b_has_michel_category);
   fChain->SetBranchAddress("has_michel_in_vertex_point_sz", &has_michel_in_vertex_point_sz, &b_has_michel_in_vertex_point_sz);
   fChain->SetBranchAddress("has_michel_in_vertex_point", has_michel_in_vertex_point, &b_has_michel_in_vertex_point);
   fChain->SetBranchAddress("has_michel_ndigits_sz", &has_michel_ndigits_sz, &b_has_michel_ndigits_sz);
   fChain->SetBranchAddress("has_michel_ndigits", has_michel_ndigits, &b_has_michel_ndigits);
   fChain->SetBranchAddress("has_michel_vertex_type_sz", &has_michel_vertex_type_sz, &b_has_michel_vertex_type_sz);
   fChain->SetBranchAddress("has_michel_vertex_type", has_michel_vertex_type, &b_has_michel_vertex_type);
   fChain->SetBranchAddress("proton_enters_front", proton_enters_front, &b_proton_enters_front);
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
   fChain->SetBranchAddress("fit_vtx", fit_vtx, &b_fit_vtx);
   fChain->SetBranchAddress("has_michel_distance_sz", &has_michel_distance_sz, &b_has_michel_distance_sz);
   fChain->SetBranchAddress("has_michel_distance", has_michel_distance, &b_has_michel_distance);
   fChain->SetBranchAddress("has_michel_energy_sz", &has_michel_energy_sz, &b_has_michel_energy_sz);
   fChain->SetBranchAddress("has_michel_energy", has_michel_energy, &b_has_michel_energy);
   fChain->SetBranchAddress("has_michel_time_diff_sz", &has_michel_time_diff_sz, &b_has_michel_time_diff_sz);
   fChain->SetBranchAddress("has_michel_time_diff", has_michel_time_diff, &b_has_michel_time_diff);
   fChain->SetBranchAddress("intraNukeProtonMomentumVec", intraNukeProtonMomentumVec, &b_intraNukeProtonMomentumVec);
   fChain->SetBranchAddress("reco_muon_vertex", reco_muon_vertex, &b_reco_muon_vertex);
   fChain->SetBranchAddress("vtxprong_isoblob_visenergy", vtxprong_isoblob_visenergy, &b_vtxprong_isoblob_visenergy);
   fChain->SetBranchAddress("vtxprong_vtxblob_visenergy", vtxprong_vtxblob_visenergy, &b_vtxprong_vtxblob_visenergy);
   fChain->SetBranchAddress("truth_has_physics_event", &truth_has_physics_event, &b_truth_has_physics_event);
   fChain->SetBranchAddress("truth_has_michel_electron", &truth_has_michel_electron, &b_truth_has_michel_electron);
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
   fChain->SetBranchAddress("truth_has_michel_from_pion_minus_momentum_sz", &truth_has_michel_from_pion_minus_momentum_sz, &b_truth_has_michel_from_pion_minus_momentum_sz);
   fChain->SetBranchAddress("truth_has_michel_from_pion_minus_momentum", truth_has_michel_from_pion_minus_momentum, &b_truth_has_michel_from_pion_minus_momentum);
   fChain->SetBranchAddress("truth_has_michel_from_pion_plus_momentum_sz", &truth_has_michel_from_pion_plus_momentum_sz, &b_truth_has_michel_from_pion_plus_momentum_sz);
   fChain->SetBranchAddress("truth_has_michel_from_pion_plus_momentum", truth_has_michel_from_pion_plus_momentum, &b_truth_has_michel_from_pion_plus_momentum);
   fChain->SetBranchAddress("NukeCCQETwoTrack_nuFlavor", &NukeCCQETwoTrack_nuFlavor, &b_NukeCCQETwoTrack_nuFlavor);
   fChain->SetBranchAddress("NukeCCQETwoTrack_nuHelicity", &NukeCCQETwoTrack_nuHelicity, &b_NukeCCQETwoTrack_nuHelicity);
   fChain->SetBranchAddress("NukeCCQETwoTrack_intCurrent", &NukeCCQETwoTrack_intCurrent, &b_NukeCCQETwoTrack_intCurrent);
   fChain->SetBranchAddress("NukeCCQETwoTrack_intType", &NukeCCQETwoTrack_intType, &b_NukeCCQETwoTrack_intType);
   fChain->SetBranchAddress("NukeCCQETwoTrack_E", &NukeCCQETwoTrack_E, &b_NukeCCQETwoTrack_E);
   fChain->SetBranchAddress("NukeCCQETwoTrack_Q2", &NukeCCQETwoTrack_Q2, &b_NukeCCQETwoTrack_Q2);
   fChain->SetBranchAddress("NukeCCQETwoTrack_x", &NukeCCQETwoTrack_x, &b_NukeCCQETwoTrack_x);
   fChain->SetBranchAddress("NukeCCQETwoTrack_y", &NukeCCQETwoTrack_y, &b_NukeCCQETwoTrack_y);
   fChain->SetBranchAddress("NukeCCQETwoTrack_W", &NukeCCQETwoTrack_W, &b_NukeCCQETwoTrack_W);
   fChain->SetBranchAddress("NukeCCQETwoTrack_score", &NukeCCQETwoTrack_score, &b_NukeCCQETwoTrack_score);
   fChain->SetBranchAddress("NukeCCQETwoTrack_leptonE", NukeCCQETwoTrack_leptonE, &b_NukeCCQETwoTrack_leptonE);
   fChain->SetBranchAddress("NukeCCQETwoTrack_vtx", NukeCCQETwoTrack_vtx, &b_NukeCCQETwoTrack_vtx);
   fChain->SetBranchAddress("NukeCCQETwoTrack_minos_trk_is_contained", &NukeCCQETwoTrack_minos_trk_is_contained, &b_NukeCCQETwoTrack_minos_trk_is_contained);
   fChain->SetBranchAddress("NukeCCQETwoTrack_minos_trk_is_ok", &NukeCCQETwoTrack_minos_trk_is_ok, &b_NukeCCQETwoTrack_minos_trk_is_ok);
   fChain->SetBranchAddress("NukeCCQETwoTrack_minos_used_range", &NukeCCQETwoTrack_minos_used_range, &b_NukeCCQETwoTrack_minos_used_range);
   fChain->SetBranchAddress("NukeCCQETwoTrack_minos_used_curvature", &NukeCCQETwoTrack_minos_used_curvature, &b_NukeCCQETwoTrack_minos_used_curvature);
   fChain->SetBranchAddress("NukeCCQETwoTrack_isMuonInsideOD", &NukeCCQETwoTrack_isMuonInsideOD, &b_NukeCCQETwoTrack_isMuonInsideOD);
   fChain->SetBranchAddress("NukeCCQETwoTrack_minos_trk_end_plane", &NukeCCQETwoTrack_minos_trk_end_plane, &b_NukeCCQETwoTrack_minos_trk_end_plane);
   fChain->SetBranchAddress("NukeCCQETwoTrack_minos_trk_quality", &NukeCCQETwoTrack_minos_trk_quality, &b_NukeCCQETwoTrack_minos_trk_quality);
   fChain->SetBranchAddress("NukeCCQETwoTrack_muon_down_hcal", &NukeCCQETwoTrack_muon_down_hcal, &b_NukeCCQETwoTrack_muon_down_hcal);
   fChain->SetBranchAddress("NukeCCQETwoTrack_muon_minos_stub", &NukeCCQETwoTrack_muon_minos_stub, &b_NukeCCQETwoTrack_muon_minos_stub);
   fChain->SetBranchAddress("NukeCCQETwoTrack_muon_minos_track", &NukeCCQETwoTrack_muon_minos_track, &b_NukeCCQETwoTrack_muon_minos_track);
   fChain->SetBranchAddress("NukeCCQETwoTrack_muon_odLastFrame", &NukeCCQETwoTrack_muon_odLastFrame, &b_NukeCCQETwoTrack_muon_odLastFrame);
   fChain->SetBranchAddress("NukeCCQETwoTrack_muon_odLastStory", &NukeCCQETwoTrack_muon_odLastStory, &b_NukeCCQETwoTrack_muon_odLastStory);
   fChain->SetBranchAddress("NukeCCQETwoTrack_muon_od_track", &NukeCCQETwoTrack_muon_od_track, &b_NukeCCQETwoTrack_muon_od_track);
   fChain->SetBranchAddress("NukeCCQETwoTrack_muon_side_ecal", &NukeCCQETwoTrack_muon_side_ecal, &b_NukeCCQETwoTrack_muon_side_ecal);
   fChain->SetBranchAddress("NukeCCQETwoTrack_muon_trk_pat_history", &NukeCCQETwoTrack_muon_trk_pat_history, &b_NukeCCQETwoTrack_muon_trk_pat_history);
   fChain->SetBranchAddress("NukeCCQETwoTrack_ntrajMuonProng", &NukeCCQETwoTrack_ntrajMuonProng, &b_NukeCCQETwoTrack_ntrajMuonProng);
   fChain->SetBranchAddress("NukeCCQETwoTrack_r_minos_trk_vtx_plane", &NukeCCQETwoTrack_r_minos_trk_vtx_plane, &b_NukeCCQETwoTrack_r_minos_trk_vtx_plane);
   fChain->SetBranchAddress("NukeCCQETwoTrack_t_minos_trk_numFSMuons", &NukeCCQETwoTrack_t_minos_trk_numFSMuons, &b_NukeCCQETwoTrack_t_minos_trk_numFSMuons);
   fChain->SetBranchAddress("NukeCCQETwoTrack_t_minos_trk_primFSLeptonPDG", &NukeCCQETwoTrack_t_minos_trk_primFSLeptonPDG, &b_NukeCCQETwoTrack_t_minos_trk_primFSLeptonPDG);
   fChain->SetBranchAddress("NukeCCQETwoTrack_targetID", &NukeCCQETwoTrack_targetID, &b_NukeCCQETwoTrack_targetID);
   fChain->SetBranchAddress("NukeCCQETwoTrack_targetZ", &NukeCCQETwoTrack_targetZ, &b_NukeCCQETwoTrack_targetZ);
   fChain->SetBranchAddress("NukeCCQETwoTrack_trajMuonProngPDG", &NukeCCQETwoTrack_trajMuonProngPDG, &b_NukeCCQETwoTrack_trajMuonProngPDG);
   fChain->SetBranchAddress("NukeCCQETwoTrack_trajMuonProngPrimary", &NukeCCQETwoTrack_trajMuonProngPrimary, &b_NukeCCQETwoTrack_trajMuonProngPrimary);
   fChain->SetBranchAddress("NukeCCQETwoTrack_endMuonTrajMomentum", &NukeCCQETwoTrack_endMuonTrajMomentum, &b_NukeCCQETwoTrack_endMuonTrajMomentum);
   fChain->SetBranchAddress("NukeCCQETwoTrack_endMuonTrajXPosition", &NukeCCQETwoTrack_endMuonTrajXPosition, &b_NukeCCQETwoTrack_endMuonTrajXPosition);
   fChain->SetBranchAddress("NukeCCQETwoTrack_endMuonTrajYPosition", &NukeCCQETwoTrack_endMuonTrajYPosition, &b_NukeCCQETwoTrack_endMuonTrajYPosition);
   fChain->SetBranchAddress("NukeCCQETwoTrack_endMuonTrajZPosition", &NukeCCQETwoTrack_endMuonTrajZPosition, &b_NukeCCQETwoTrack_endMuonTrajZPosition);
   fChain->SetBranchAddress("NukeCCQETwoTrack_minos_trk_bave", &NukeCCQETwoTrack_minos_trk_bave, &b_NukeCCQETwoTrack_minos_trk_bave);
   fChain->SetBranchAddress("NukeCCQETwoTrack_minos_trk_chi2", &NukeCCQETwoTrack_minos_trk_chi2, &b_NukeCCQETwoTrack_minos_trk_chi2);
   fChain->SetBranchAddress("NukeCCQETwoTrack_minos_trk_end_u", &NukeCCQETwoTrack_minos_trk_end_u, &b_NukeCCQETwoTrack_minos_trk_end_u);
   fChain->SetBranchAddress("NukeCCQETwoTrack_minos_trk_end_v", &NukeCCQETwoTrack_minos_trk_end_v, &b_NukeCCQETwoTrack_minos_trk_end_v);
   fChain->SetBranchAddress("NukeCCQETwoTrack_minos_trk_end_x", &NukeCCQETwoTrack_minos_trk_end_x, &b_NukeCCQETwoTrack_minos_trk_end_x);
   fChain->SetBranchAddress("NukeCCQETwoTrack_minos_trk_end_y", &NukeCCQETwoTrack_minos_trk_end_y, &b_NukeCCQETwoTrack_minos_trk_end_y);
   fChain->SetBranchAddress("NukeCCQETwoTrack_minos_trk_end_z", &NukeCCQETwoTrack_minos_trk_end_z, &b_NukeCCQETwoTrack_minos_trk_end_z);
   fChain->SetBranchAddress("NukeCCQETwoTrack_minos_trk_eqp", &NukeCCQETwoTrack_minos_trk_eqp, &b_NukeCCQETwoTrack_minos_trk_eqp);
   fChain->SetBranchAddress("NukeCCQETwoTrack_minos_trk_eqp_qp", &NukeCCQETwoTrack_minos_trk_eqp_qp, &b_NukeCCQETwoTrack_minos_trk_eqp_qp);
   fChain->SetBranchAddress("NukeCCQETwoTrack_minos_trk_fit_pass", &NukeCCQETwoTrack_minos_trk_fit_pass, &b_NukeCCQETwoTrack_minos_trk_fit_pass);
   fChain->SetBranchAddress("NukeCCQETwoTrack_minos_trk_ndf", &NukeCCQETwoTrack_minos_trk_ndf, &b_NukeCCQETwoTrack_minos_trk_ndf);
   fChain->SetBranchAddress("NukeCCQETwoTrack_minos_trk_p", &NukeCCQETwoTrack_minos_trk_p, &b_NukeCCQETwoTrack_minos_trk_p);
   fChain->SetBranchAddress("NukeCCQETwoTrack_minos_trk_p_curvature", &NukeCCQETwoTrack_minos_trk_p_curvature, &b_NukeCCQETwoTrack_minos_trk_p_curvature);
   fChain->SetBranchAddress("NukeCCQETwoTrack_minos_trk_p_range", &NukeCCQETwoTrack_minos_trk_p_range, &b_NukeCCQETwoTrack_minos_trk_p_range);
   fChain->SetBranchAddress("NukeCCQETwoTrack_minos_trk_qp", &NukeCCQETwoTrack_minos_trk_qp, &b_NukeCCQETwoTrack_minos_trk_qp);
   fChain->SetBranchAddress("NukeCCQETwoTrack_minos_trk_vtx_x", &NukeCCQETwoTrack_minos_trk_vtx_x, &b_NukeCCQETwoTrack_minos_trk_vtx_x);
   fChain->SetBranchAddress("NukeCCQETwoTrack_minos_trk_vtx_y", &NukeCCQETwoTrack_minos_trk_vtx_y, &b_NukeCCQETwoTrack_minos_trk_vtx_y);
   fChain->SetBranchAddress("NukeCCQETwoTrack_minos_trk_vtx_z", &NukeCCQETwoTrack_minos_trk_vtx_z, &b_NukeCCQETwoTrack_minos_trk_vtx_z);
   fChain->SetBranchAddress("NukeCCQETwoTrack_muon_enu", &NukeCCQETwoTrack_muon_enu, &b_NukeCCQETwoTrack_muon_enu);
   fChain->SetBranchAddress("NukeCCQETwoTrack_muon_odClustersAvgTime", &NukeCCQETwoTrack_muon_odClustersAvgTime, &b_NukeCCQETwoTrack_muon_odClustersAvgTime);
   fChain->SetBranchAddress("NukeCCQETwoTrack_muon_odElossMomentum", &NukeCCQETwoTrack_muon_odElossMomentum, &b_NukeCCQETwoTrack_muon_odElossMomentum);
   fChain->SetBranchAddress("NukeCCQETwoTrack_muon_odEndX", &NukeCCQETwoTrack_muon_odEndX, &b_NukeCCQETwoTrack_muon_odEndX);
   fChain->SetBranchAddress("NukeCCQETwoTrack_muon_odEndY", &NukeCCQETwoTrack_muon_odEndY, &b_NukeCCQETwoTrack_muon_odEndY);
   fChain->SetBranchAddress("NukeCCQETwoTrack_muon_odEndZ", &NukeCCQETwoTrack_muon_odEndZ, &b_NukeCCQETwoTrack_muon_odEndZ);
   fChain->SetBranchAddress("NukeCCQETwoTrack_muon_odFaceX", &NukeCCQETwoTrack_muon_odFaceX, &b_NukeCCQETwoTrack_muon_odFaceX);
   fChain->SetBranchAddress("NukeCCQETwoTrack_muon_odFaceY", &NukeCCQETwoTrack_muon_odFaceY, &b_NukeCCQETwoTrack_muon_odFaceY);
   fChain->SetBranchAddress("NukeCCQETwoTrack_muon_odFaceZ", &NukeCCQETwoTrack_muon_odFaceZ, &b_NukeCCQETwoTrack_muon_odFaceZ);
   fChain->SetBranchAddress("NukeCCQETwoTrack_muon_odLastClusZ", &NukeCCQETwoTrack_muon_odLastClusZ, &b_NukeCCQETwoTrack_muon_odLastClusZ);
   fChain->SetBranchAddress("NukeCCQETwoTrack_muon_odStopDistMomentum", &NukeCCQETwoTrack_muon_odStopDistMomentum, &b_NukeCCQETwoTrack_muon_odStopDistMomentum);
   fChain->SetBranchAddress("NukeCCQETwoTrack_muon_odTrackAvgTime", &NukeCCQETwoTrack_muon_odTrackAvgTime, &b_NukeCCQETwoTrack_muon_odTrackAvgTime);
   fChain->SetBranchAddress("NukeCCQETwoTrack_muon_phi", &NukeCCQETwoTrack_muon_phi, &b_NukeCCQETwoTrack_muon_phi);
   fChain->SetBranchAddress("NukeCCQETwoTrack_muon_q2", &NukeCCQETwoTrack_muon_q2, &b_NukeCCQETwoTrack_muon_q2);
   fChain->SetBranchAddress("NukeCCQETwoTrack_muon_score", &NukeCCQETwoTrack_muon_score, &b_NukeCCQETwoTrack_muon_score);
   fChain->SetBranchAddress("NukeCCQETwoTrack_muon_theta", &NukeCCQETwoTrack_muon_theta, &b_NukeCCQETwoTrack_muon_theta);
   fChain->SetBranchAddress("NukeCCQETwoTrack_muon_thetaX", &NukeCCQETwoTrack_muon_thetaX, &b_NukeCCQETwoTrack_muon_thetaX);
   fChain->SetBranchAddress("NukeCCQETwoTrack_muon_thetaY", &NukeCCQETwoTrack_muon_thetaY, &b_NukeCCQETwoTrack_muon_thetaY);
   fChain->SetBranchAddress("NukeCCQETwoTrack_r_minos_trk_bdL", &NukeCCQETwoTrack_r_minos_trk_bdL, &b_NukeCCQETwoTrack_r_minos_trk_bdL);
   fChain->SetBranchAddress("NukeCCQETwoTrack_r_minos_trk_end_dcosx", &NukeCCQETwoTrack_r_minos_trk_end_dcosx, &b_NukeCCQETwoTrack_r_minos_trk_end_dcosx);
   fChain->SetBranchAddress("NukeCCQETwoTrack_r_minos_trk_end_dcosy", &NukeCCQETwoTrack_r_minos_trk_end_dcosy, &b_NukeCCQETwoTrack_r_minos_trk_end_dcosy);
   fChain->SetBranchAddress("NukeCCQETwoTrack_r_minos_trk_end_dcosz", &NukeCCQETwoTrack_r_minos_trk_end_dcosz, &b_NukeCCQETwoTrack_r_minos_trk_end_dcosz);
   fChain->SetBranchAddress("NukeCCQETwoTrack_r_minos_trk_vtx_dcosx", &NukeCCQETwoTrack_r_minos_trk_vtx_dcosx, &b_NukeCCQETwoTrack_r_minos_trk_vtx_dcosx);
   fChain->SetBranchAddress("NukeCCQETwoTrack_r_minos_trk_vtx_dcosy", &NukeCCQETwoTrack_r_minos_trk_vtx_dcosy, &b_NukeCCQETwoTrack_r_minos_trk_vtx_dcosy);
   fChain->SetBranchAddress("NukeCCQETwoTrack_r_minos_trk_vtx_dcosz", &NukeCCQETwoTrack_r_minos_trk_vtx_dcosz, &b_NukeCCQETwoTrack_r_minos_trk_vtx_dcosz);
   fChain->SetBranchAddress("NukeCCQETwoTrack_t_minos_trk_primFSLepMinosInitProjPx", &NukeCCQETwoTrack_t_minos_trk_primFSLepMinosInitProjPx, &b_NukeCCQETwoTrack_t_minos_trk_primFSLepMinosInitProjPx);
   fChain->SetBranchAddress("NukeCCQETwoTrack_t_minos_trk_primFSLepMinosInitProjPy", &NukeCCQETwoTrack_t_minos_trk_primFSLepMinosInitProjPy, &b_NukeCCQETwoTrack_t_minos_trk_primFSLepMinosInitProjPy);
   fChain->SetBranchAddress("NukeCCQETwoTrack_t_minos_trk_primFSLepMinosInitProjPz", &NukeCCQETwoTrack_t_minos_trk_primFSLepMinosInitProjPz, &b_NukeCCQETwoTrack_t_minos_trk_primFSLepMinosInitProjPz);
   fChain->SetBranchAddress("NukeCCQETwoTrack_t_minos_trk_primFSLepMinosInitProjX", &NukeCCQETwoTrack_t_minos_trk_primFSLepMinosInitProjX, &b_NukeCCQETwoTrack_t_minos_trk_primFSLepMinosInitProjX);
   fChain->SetBranchAddress("NukeCCQETwoTrack_t_minos_trk_primFSLepMinosInitProjY", &NukeCCQETwoTrack_t_minos_trk_primFSLepMinosInitProjY, &b_NukeCCQETwoTrack_t_minos_trk_primFSLepMinosInitProjY);
   fChain->SetBranchAddress("NukeCCQETwoTrack_t_minos_trk_primFSLepMinosInitProjZ", &NukeCCQETwoTrack_t_minos_trk_primFSLepMinosInitProjZ, &b_NukeCCQETwoTrack_t_minos_trk_primFSLepMinosInitProjZ);
   fChain->SetBranchAddress("NukeCCQETwoTrack_t_minos_trk_primFSLepMnvFinalPx", &NukeCCQETwoTrack_t_minos_trk_primFSLepMnvFinalPx, &b_NukeCCQETwoTrack_t_minos_trk_primFSLepMnvFinalPx);
   fChain->SetBranchAddress("NukeCCQETwoTrack_t_minos_trk_primFSLepMnvFinalPy", &NukeCCQETwoTrack_t_minos_trk_primFSLepMnvFinalPy, &b_NukeCCQETwoTrack_t_minos_trk_primFSLepMnvFinalPy);
   fChain->SetBranchAddress("NukeCCQETwoTrack_t_minos_trk_primFSLepMnvFinalPz", &NukeCCQETwoTrack_t_minos_trk_primFSLepMnvFinalPz, &b_NukeCCQETwoTrack_t_minos_trk_primFSLepMnvFinalPz);
   fChain->SetBranchAddress("NukeCCQETwoTrack_t_minos_trk_primFSLepMnvFinalX", &NukeCCQETwoTrack_t_minos_trk_primFSLepMnvFinalX, &b_NukeCCQETwoTrack_t_minos_trk_primFSLepMnvFinalX);
   fChain->SetBranchAddress("NukeCCQETwoTrack_t_minos_trk_primFSLepMnvFinalY", &NukeCCQETwoTrack_t_minos_trk_primFSLepMnvFinalY, &b_NukeCCQETwoTrack_t_minos_trk_primFSLepMnvFinalY);
   fChain->SetBranchAddress("NukeCCQETwoTrack_t_minos_trk_primFSLepMnvFinalZ", &NukeCCQETwoTrack_t_minos_trk_primFSLepMnvFinalZ, &b_NukeCCQETwoTrack_t_minos_trk_primFSLepMnvFinalZ);
   fChain->SetBranchAddress("NukeCCQETwoTrack_t_minos_trk_primFSLepMnvInitPx", &NukeCCQETwoTrack_t_minos_trk_primFSLepMnvInitPx, &b_NukeCCQETwoTrack_t_minos_trk_primFSLepMnvInitPx);
   fChain->SetBranchAddress("NukeCCQETwoTrack_t_minos_trk_primFSLepMnvInitPy", &NukeCCQETwoTrack_t_minos_trk_primFSLepMnvInitPy, &b_NukeCCQETwoTrack_t_minos_trk_primFSLepMnvInitPy);
   fChain->SetBranchAddress("NukeCCQETwoTrack_t_minos_trk_primFSLepMnvInitPz", &NukeCCQETwoTrack_t_minos_trk_primFSLepMnvInitPz, &b_NukeCCQETwoTrack_t_minos_trk_primFSLepMnvInitPz);
   fChain->SetBranchAddress("NukeCCQETwoTrack_t_minos_trk_primFSLepMnvInitX", &NukeCCQETwoTrack_t_minos_trk_primFSLepMnvInitX, &b_NukeCCQETwoTrack_t_minos_trk_primFSLepMnvInitX);
   fChain->SetBranchAddress("NukeCCQETwoTrack_t_minos_trk_primFSLepMnvInitY", &NukeCCQETwoTrack_t_minos_trk_primFSLepMnvInitY, &b_NukeCCQETwoTrack_t_minos_trk_primFSLepMnvInitY);
   fChain->SetBranchAddress("NukeCCQETwoTrack_t_minos_trk_primFSLepMnvInitZ", &NukeCCQETwoTrack_t_minos_trk_primFSLepMnvInitZ, &b_NukeCCQETwoTrack_t_minos_trk_primFSLepMnvInitZ);
   fChain->SetBranchAddress("NukeCCQETwoTrack_targetZPos", &NukeCCQETwoTrack_targetZPos, &b_NukeCCQETwoTrack_targetZPos);
   fChain->SetBranchAddress("NukeCCQETwoTrack_trajMuonPhi", &NukeCCQETwoTrack_trajMuonPhi, &b_NukeCCQETwoTrack_trajMuonPhi);
   fChain->SetBranchAddress("NukeCCQETwoTrack_trajMuonProngMomentum", &NukeCCQETwoTrack_trajMuonProngMomentum, &b_NukeCCQETwoTrack_trajMuonProngMomentum);
   fChain->SetBranchAddress("NukeCCQETwoTrack_trajMuonTheta", &NukeCCQETwoTrack_trajMuonTheta, &b_NukeCCQETwoTrack_trajMuonTheta);
   fChain->SetBranchAddress("NukeCCQETwoTrack_isProtonInsideOD", NukeCCQETwoTrack_isProtonInsideOD, &b_NukeCCQETwoTrack_isProtonInsideOD);
   fChain->SetBranchAddress("NukeCCQETwoTrack_ntrajProngProng", NukeCCQETwoTrack_ntrajProngProng, &b_NukeCCQETwoTrack_ntrajProngProng);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_kinked", NukeCCQETwoTrack_proton_kinked, &b_NukeCCQETwoTrack_proton_kinked);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_odMatch", NukeCCQETwoTrack_proton_odMatch, &b_NukeCCQETwoTrack_proton_odMatch);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_trk_pat_history", NukeCCQETwoTrack_proton_trk_pat_history, &b_NukeCCQETwoTrack_proton_trk_pat_history);
   fChain->SetBranchAddress("NukeCCQETwoTrack_trajProtonProngPDG", NukeCCQETwoTrack_trajProtonProngPDG, &b_NukeCCQETwoTrack_trajProtonProngPDG);
   fChain->SetBranchAddress("NukeCCQETwoTrack_trajProtonProngPrimary", NukeCCQETwoTrack_trajProtonProngPrimary, &b_NukeCCQETwoTrack_trajProtonProngPrimary);
   fChain->SetBranchAddress("NukeCCQETwoTrack_coplanarAngle", NukeCCQETwoTrack_coplanarAngle, &b_NukeCCQETwoTrack_coplanarAngle);
   fChain->SetBranchAddress("NukeCCQETwoTrack_endProtonTrajMomentum", NukeCCQETwoTrack_endProtonTrajMomentum, &b_NukeCCQETwoTrack_endProtonTrajMomentum);
   fChain->SetBranchAddress("NukeCCQETwoTrack_endProtonTrajXPosition", NukeCCQETwoTrack_endProtonTrajXPosition, &b_NukeCCQETwoTrack_endProtonTrajXPosition);
   fChain->SetBranchAddress("NukeCCQETwoTrack_endProtonTrajYPosition", NukeCCQETwoTrack_endProtonTrajYPosition, &b_NukeCCQETwoTrack_endProtonTrajYPosition);
   fChain->SetBranchAddress("NukeCCQETwoTrack_endProtonTrajZPosition", NukeCCQETwoTrack_endProtonTrajZPosition, &b_NukeCCQETwoTrack_endProtonTrajZPosition);
   fChain->SetBranchAddress("NukeCCQETwoTrack_muon_endPoint", NukeCCQETwoTrack_muon_endPoint, &b_NukeCCQETwoTrack_muon_endPoint);
   fChain->SetBranchAddress("NukeCCQETwoTrack_muon_startPoint", NukeCCQETwoTrack_muon_startPoint, &b_NukeCCQETwoTrack_muon_startPoint);
   fChain->SetBranchAddress("NukeCCQETwoTrack_open_angle", NukeCCQETwoTrack_open_angle, &b_NukeCCQETwoTrack_open_angle);
   fChain->SetBranchAddress("NukeCCQETwoTrack_pion_chi2_ndf", NukeCCQETwoTrack_pion_chi2_ndf, &b_NukeCCQETwoTrack_pion_chi2_ndf);
   fChain->SetBranchAddress("NukeCCQETwoTrack_pion_score", NukeCCQETwoTrack_pion_score, &b_NukeCCQETwoTrack_pion_score);
   fChain->SetBranchAddress("NukeCCQETwoTrack_pion_score1", NukeCCQETwoTrack_pion_score1, &b_NukeCCQETwoTrack_pion_score1);
   fChain->SetBranchAddress("NukeCCQETwoTrack_pion_score2", NukeCCQETwoTrack_pion_score2, &b_NukeCCQETwoTrack_pion_score2);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_E", NukeCCQETwoTrack_proton_E, &b_NukeCCQETwoTrack_proton_E);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_chi2_ndf", NukeCCQETwoTrack_proton_chi2_ndf, &b_NukeCCQETwoTrack_proton_chi2_ndf);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_ekin", NukeCCQETwoTrack_proton_ekin, &b_NukeCCQETwoTrack_proton_ekin);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_endPointX", NukeCCQETwoTrack_proton_endPointX, &b_NukeCCQETwoTrack_proton_endPointX);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_endPointY", NukeCCQETwoTrack_proton_endPointY, &b_NukeCCQETwoTrack_proton_endPointY);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_endPointZ", NukeCCQETwoTrack_proton_endPointZ, &b_NukeCCQETwoTrack_proton_endPointZ);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_enu", NukeCCQETwoTrack_proton_enu, &b_NukeCCQETwoTrack_proton_enu);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_momentum_shift_BetheBloch_Down", NukeCCQETwoTrack_proton_momentum_shift_BetheBloch_Down, &b_NukeCCQETwoTrack_proton_momentum_shift_BetheBloch_Down);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_momentum_shift_BetheBloch_Up", NukeCCQETwoTrack_proton_momentum_shift_BetheBloch_Up, &b_NukeCCQETwoTrack_proton_momentum_shift_BetheBloch_Up);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_momentum_shift_Birks", NukeCCQETwoTrack_proton_momentum_shift_Birks, &b_NukeCCQETwoTrack_proton_momentum_shift_Birks);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_momentum_shift_MEU_Down", NukeCCQETwoTrack_proton_momentum_shift_MEU_Down, &b_NukeCCQETwoTrack_proton_momentum_shift_MEU_Down);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_momentum_shift_MEU_Up", NukeCCQETwoTrack_proton_momentum_shift_MEU_Up, &b_NukeCCQETwoTrack_proton_momentum_shift_MEU_Up);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_momentum_shift_Mass_Down", NukeCCQETwoTrack_proton_momentum_shift_Mass_Down, &b_NukeCCQETwoTrack_proton_momentum_shift_Mass_Down);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_momentum_shift_Mass_Up", NukeCCQETwoTrack_proton_momentum_shift_Mass_Up, &b_NukeCCQETwoTrack_proton_momentum_shift_Mass_Up);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_p", NukeCCQETwoTrack_proton_p, &b_NukeCCQETwoTrack_proton_p);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_p_calCorrection", NukeCCQETwoTrack_proton_p_calCorrection, &b_NukeCCQETwoTrack_proton_p_calCorrection);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_p_dEdXTool", NukeCCQETwoTrack_proton_p_dEdXTool, &b_NukeCCQETwoTrack_proton_p_dEdXTool);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_p_visEnergy", NukeCCQETwoTrack_proton_p_visEnergy, &b_NukeCCQETwoTrack_proton_p_visEnergy);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_phi", NukeCCQETwoTrack_proton_phi, &b_NukeCCQETwoTrack_proton_phi);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_px", NukeCCQETwoTrack_proton_px, &b_NukeCCQETwoTrack_proton_px);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_py", NukeCCQETwoTrack_proton_py, &b_NukeCCQETwoTrack_proton_py);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_pz", NukeCCQETwoTrack_proton_pz, &b_NukeCCQETwoTrack_proton_pz);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_q2", NukeCCQETwoTrack_proton_q2, &b_NukeCCQETwoTrack_proton_q2);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_score", NukeCCQETwoTrack_proton_score, &b_NukeCCQETwoTrack_proton_score);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_score1", NukeCCQETwoTrack_proton_score1, &b_NukeCCQETwoTrack_proton_score1);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_score1_shift_BetheBloch_Down", NukeCCQETwoTrack_proton_score1_shift_BetheBloch_Down, &b_NukeCCQETwoTrack_proton_score1_shift_BetheBloch_Down);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_score1_shift_BetheBloch_Up", NukeCCQETwoTrack_proton_score1_shift_BetheBloch_Up, &b_NukeCCQETwoTrack_proton_score1_shift_BetheBloch_Up);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_score1_shift_Birks", NukeCCQETwoTrack_proton_score1_shift_Birks, &b_NukeCCQETwoTrack_proton_score1_shift_Birks);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_score1_shift_MEU_Down", NukeCCQETwoTrack_proton_score1_shift_MEU_Down, &b_NukeCCQETwoTrack_proton_score1_shift_MEU_Down);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_score1_shift_MEU_Up", NukeCCQETwoTrack_proton_score1_shift_MEU_Up, &b_NukeCCQETwoTrack_proton_score1_shift_MEU_Up);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_score1_shift_Mass_Down", NukeCCQETwoTrack_proton_score1_shift_Mass_Down, &b_NukeCCQETwoTrack_proton_score1_shift_Mass_Down);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_score1_shift_Mass_Up", NukeCCQETwoTrack_proton_score1_shift_Mass_Up, &b_NukeCCQETwoTrack_proton_score1_shift_Mass_Up);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_score2", NukeCCQETwoTrack_proton_score2, &b_NukeCCQETwoTrack_proton_score2);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_startPointX", NukeCCQETwoTrack_proton_startPointX, &b_NukeCCQETwoTrack_proton_startPointX);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_startPointY", NukeCCQETwoTrack_proton_startPointY, &b_NukeCCQETwoTrack_proton_startPointY);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_startPointZ", NukeCCQETwoTrack_proton_startPointZ, &b_NukeCCQETwoTrack_proton_startPointZ);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_theta", NukeCCQETwoTrack_proton_theta, &b_NukeCCQETwoTrack_proton_theta);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_thetaX", NukeCCQETwoTrack_proton_thetaX, &b_NukeCCQETwoTrack_proton_thetaX);
   fChain->SetBranchAddress("NukeCCQETwoTrack_proton_thetaY", NukeCCQETwoTrack_proton_thetaY, &b_NukeCCQETwoTrack_proton_thetaY);
   fChain->SetBranchAddress("NukeCCQETwoTrack_sys_muon_curve_energy_shift", NukeCCQETwoTrack_sys_muon_curve_energy_shift, &b_NukeCCQETwoTrack_sys_muon_curve_energy_shift);
   fChain->SetBranchAddress("NukeCCQETwoTrack_sys_muon_energy_shift", NukeCCQETwoTrack_sys_muon_energy_shift, &b_NukeCCQETwoTrack_sys_muon_energy_shift);
   fChain->SetBranchAddress("NukeCCQETwoTrack_sys_muon_minerva_energy_shift", NukeCCQETwoTrack_sys_muon_minerva_energy_shift, &b_NukeCCQETwoTrack_sys_muon_minerva_energy_shift);
   fChain->SetBranchAddress("NukeCCQETwoTrack_sys_muon_qSquared_shift", NukeCCQETwoTrack_sys_muon_qSquared_shift, &b_NukeCCQETwoTrack_sys_muon_qSquared_shift);
   fChain->SetBranchAddress("NukeCCQETwoTrack_sys_muon_range_energy_shift", NukeCCQETwoTrack_sys_muon_range_energy_shift, &b_NukeCCQETwoTrack_sys_muon_range_energy_shift);
   fChain->SetBranchAddress("NukeCCQETwoTrack_sys_muon_wSquared_shift", NukeCCQETwoTrack_sys_muon_wSquared_shift, &b_NukeCCQETwoTrack_sys_muon_wSquared_shift);
   fChain->SetBranchAddress("NukeCCQETwoTrack_sys_muon_xbj_shift", NukeCCQETwoTrack_sys_muon_xbj_shift, &b_NukeCCQETwoTrack_sys_muon_xbj_shift);
   fChain->SetBranchAddress("NukeCCQETwoTrack_sys_muon_y_shift", NukeCCQETwoTrack_sys_muon_y_shift, &b_NukeCCQETwoTrack_sys_muon_y_shift);
   fChain->SetBranchAddress("NukeCCQETwoTrack_sys_nu_energy_shift", NukeCCQETwoTrack_sys_nu_energy_shift, &b_NukeCCQETwoTrack_sys_nu_energy_shift);
   fChain->SetBranchAddress("NukeCCQETwoTrack_sys_recoil_energy_shift", NukeCCQETwoTrack_sys_recoil_energy_shift, &b_NukeCCQETwoTrack_sys_recoil_energy_shift);
   fChain->SetBranchAddress("NukeCCQETwoTrack_sys_recoil_qSquared_shift", NukeCCQETwoTrack_sys_recoil_qSquared_shift, &b_NukeCCQETwoTrack_sys_recoil_qSquared_shift);
   fChain->SetBranchAddress("NukeCCQETwoTrack_sys_recoil_wSquared_shift", NukeCCQETwoTrack_sys_recoil_wSquared_shift, &b_NukeCCQETwoTrack_sys_recoil_wSquared_shift);
   fChain->SetBranchAddress("NukeCCQETwoTrack_sys_recoil_xbj_shift", NukeCCQETwoTrack_sys_recoil_xbj_shift, &b_NukeCCQETwoTrack_sys_recoil_xbj_shift);
   fChain->SetBranchAddress("NukeCCQETwoTrack_sys_recoil_y_shift", NukeCCQETwoTrack_sys_recoil_y_shift, &b_NukeCCQETwoTrack_sys_recoil_y_shift);
   fChain->SetBranchAddress("NukeCCQETwoTrack_sys_total_qSquared_shift", NukeCCQETwoTrack_sys_total_qSquared_shift, &b_NukeCCQETwoTrack_sys_total_qSquared_shift);
   fChain->SetBranchAddress("NukeCCQETwoTrack_sys_total_wSquared_shift", NukeCCQETwoTrack_sys_total_wSquared_shift, &b_NukeCCQETwoTrack_sys_total_wSquared_shift);
   fChain->SetBranchAddress("NukeCCQETwoTrack_sys_total_xbj_shift", NukeCCQETwoTrack_sys_total_xbj_shift, &b_NukeCCQETwoTrack_sys_total_xbj_shift);
   fChain->SetBranchAddress("NukeCCQETwoTrack_sys_total_y_shift", NukeCCQETwoTrack_sys_total_y_shift, &b_NukeCCQETwoTrack_sys_total_y_shift);
   fChain->SetBranchAddress("NukeCCQETwoTrack_trajProtonPhi", NukeCCQETwoTrack_trajProtonPhi, &b_NukeCCQETwoTrack_trajProtonPhi);
   fChain->SetBranchAddress("NukeCCQETwoTrack_trajProtonProngMomentum", NukeCCQETwoTrack_trajProtonProngMomentum, &b_NukeCCQETwoTrack_trajProtonProngMomentum);
   fChain->SetBranchAddress("NukeCCQETwoTrack_trajProtonTheta", NukeCCQETwoTrack_trajProtonTheta, &b_NukeCCQETwoTrack_trajProtonTheta);
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
   fChain->SetBranchAddress("mc_wgt_Norm_sz", &mc_wgt_Norm_sz, &b_mc_wgt_Norm_sz);
   fChain->SetBranchAddress("mc_wgt_Norm", &mc_wgt_Norm, &b_mc_wgt_Norm);
   fChain->SetBranchAddress("mc_wgt_ppfx_MIPPNumiYields_sz", &mc_wgt_ppfx_MIPPNumiYields_sz, &b_mc_wgt_ppfx_MIPPNumiYields_sz);
   fChain->SetBranchAddress("mc_wgt_ppfx_MIPPNumiYields", &mc_wgt_ppfx_MIPPNumiYields, &b_mc_wgt_ppfx_MIPPNumiYields);
   fChain->SetBranchAddress("mc_wgt_ppfx_TargetAttenuation_sz", &mc_wgt_ppfx_TargetAttenuation_sz, &b_mc_wgt_ppfx_TargetAttenuation_sz);
   fChain->SetBranchAddress("mc_wgt_ppfx_TargetAttenuation", &mc_wgt_ppfx_TargetAttenuation, &b_mc_wgt_ppfx_TargetAttenuation);
   fChain->SetBranchAddress("mc_wgt_ppfx_NA49_sz", &mc_wgt_ppfx_NA49_sz, &b_mc_wgt_ppfx_NA49_sz);
   fChain->SetBranchAddress("mc_wgt_ppfx_NA49", &mc_wgt_ppfx_NA49, &b_mc_wgt_ppfx_NA49);
   fChain->SetBranchAddress("mc_wgt_ppfx_MIPPKaonsYields_sz", &mc_wgt_ppfx_MIPPKaonsYields_sz, &b_mc_wgt_ppfx_MIPPKaonsYields_sz);
   fChain->SetBranchAddress("mc_wgt_ppfx_MIPPKaonsYields", &mc_wgt_ppfx_MIPPKaonsYields, &b_mc_wgt_ppfx_MIPPKaonsYields);
   fChain->SetBranchAddress("mc_wgt_ppfx_MIPPThinTarget_sz", &mc_wgt_ppfx_MIPPThinTarget_sz, &b_mc_wgt_ppfx_MIPPThinTarget_sz);
   fChain->SetBranchAddress("mc_wgt_ppfx_MIPPThinTarget", &mc_wgt_ppfx_MIPPThinTarget, &b_mc_wgt_ppfx_MIPPThinTarget);
   fChain->SetBranchAddress("mc_wgt_ppfx_Absorption_sz", &mc_wgt_ppfx_Absorption_sz, &b_mc_wgt_ppfx_Absorption_sz);
   fChain->SetBranchAddress("mc_wgt_ppfx_Absorption", &mc_wgt_ppfx_Absorption, &b_mc_wgt_ppfx_Absorption);
   fChain->SetBranchAddress("mc_wgt_ppfx_Others_sz", &mc_wgt_ppfx_Others_sz, &b_mc_wgt_ppfx_Others_sz);
   fChain->SetBranchAddress("mc_wgt_ppfx_Others", &mc_wgt_ppfx_Others, &b_mc_wgt_ppfx_Others);
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

Bool_t NukeCCQETwoTrack_default::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void NukeCCQETwoTrack_default::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t NukeCCQETwoTrack_default::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef NukeCCQETwoTrack_default_cxx
