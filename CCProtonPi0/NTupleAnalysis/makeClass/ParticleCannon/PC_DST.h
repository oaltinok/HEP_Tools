#ifndef PC_DST_h
#define PC_DST_h

// C++ Libraries
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

// ROOT Libraries
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>

// Custom Libraries
#include "../../Libraries/Folder_List.h"
#include "../../Libraries/Data_Functions.h"
#include "../../Libraries/HEP_Functions.h"

using namespace std;

class PC_DST {
    
    public :
        // Custom Functions
        double GetEvis(int min_module, int max_module);
        double CalculateEnergy();
        void InitHistograms();
        void writeHistograms();

        // Histograms
        TH1D* error;

        // Custom Variables
        string rootDir;
        TFile* f;

        // Data
        TTree          *fChain;   //!pointer to the analyzed TTree or TChain
        Int_t           fCurrent; //!current Tree number in a TChain

        // Declaration of leaf types
        Int_t           fmwk_v;
        Int_t           fmwk_r;
        Int_t           fmwk_p;
        Int_t           n_slices;
        Int_t           ev_detector;
        Int_t           ev_det_config;
        Int_t           ev_run;
        Int_t           ev_sub_run;
        Int_t           ev_trigger_type;
        Int_t           ev_cal_settings;
        Int_t           ev_gl_gate;
        Int_t           ev_gate;
        Int_t           ev_gps_time_sec;
        Int_t           ev_gps_time_usec;
        Int_t           ev_readout_time;
        Int_t           ev_errors;
        Int_t           ev_nADC_Frames;
        Int_t           ev_nDisc_Frames;
        Int_t           ev_nFPGA_Frames;
        Int_t           n_febs;
        Int_t           feb_id[509];   //[n_febs]
        Int_t           feb_hv_on[509];   //[n_febs]
        Int_t           feb_hv_targ[509];   //[n_febs]
        Int_t           feb_hv_act[509];   //[n_febs]
        Int_t           feb_hv_per_auto[509];   //[n_febs]
        Double_t        feb_temperature[509];   //[n_febs]
        Int_t           feb_gate_time_stamp[509];   //[n_febs]
        Int_t           n_rawhits;
        Int_t           hit_feb_id[151];   //[n_rawhits]
        UInt_t          hit_flags[151];   //[n_rawhits]
        Int_t           hit_channel_id[151];   //[n_rawhits]
        Int_t           hit_index[151];   //[n_rawhits]
        Int_t           hit_location[151];   //[n_rawhits]
        Int_t           hit_is_mc[151];   //[n_rawhits]
        Int_t           hit_num[151];   //[n_rawhits]
        Int_t           hit_pixel[151];   //[n_rawhits]
        Int_t           hit_board[151];   //[n_rawhits]
        Int_t           hit_chain[151];   //[n_rawhits]
        Int_t           hit_croc[151];   //[n_rawhits]
        Int_t           hit_crate[151];   //[n_rawhits]
        Int_t           hit_link[151];   //[n_rawhits]
        Int_t           hit_disc_fired[151];   //[n_rawhits]
        Int_t           hit_sys_ticks[151];   //[n_rawhits]
        Int_t           hit_delay_ticks[151];   //[n_rawhits]
        Int_t           hit_quarter_ticks[151];   //[n_rawhits]
        Int_t           hit_qlo[151];   //[n_rawhits]
        Int_t           hit_qmed[151];   //[n_rawhits]
        Int_t           hit_qhi[151];   //[n_rawhits]
        Int_t           n_idhits;
        Double_t        hits_id_per_mod;
        Int_t           hit_strip[151];   //[n_rawhits]
        Int_t           hit_plane[151];   //[n_rawhits]
        Int_t           hit_module[151];   //[n_rawhits]
        Int_t           hit_view[151];   //[n_rawhits]
        Int_t           n_odhits;
        Double_t        hits_od_per_mod;
        Int_t           hit_bar[151];   //[n_rawhits]
        Int_t           hit_story[151];   //[n_rawhits]
        Int_t           hit_tower[151];   //[n_rawhits]
        Int_t           hit_frame[151];   //[n_rawhits]
        Int_t           n_vetohits;
        Int_t           hit_wall[151];   //[n_rawhits]
        Int_t           hit_paddle[151];   //[n_rawhits]
        Int_t           hit_pmt[151];   //[n_rawhits]
        Double_t        hit_q[151];   //[n_rawhits]
        Double_t        hit_pe[151];   //[n_rawhits]
        Double_t        hit_norm_energy[151];   //[n_rawhits]
        Double_t        hit_time_raw[151];   //[n_rawhits]
        Double_t        hit_time[151];   //[n_rawhits]
        Int_t           hit_time_slice[151];   //[n_rawhits]
        Double_t        hits_total_pe[151];   //[n_rawhits]
        Int_t           hit_user_color[151];   //[n_rawhits]
        Int_t           n_clusters_id;
        Int_t           clus_id_index[121];   //[n_clusters_id]
        Int_t           clus_id_strip[121];   //[n_clusters_id]
        Int_t           clus_id_plane[121];   //[n_clusters_id]
        Int_t           clus_id_module[121];   //[n_clusters_id]
        Double_t        clus_id_coord[121];   //[n_clusters_id]
        Double_t        clus_id_coordErr[121];   //[n_clusters_id]
        Double_t        clus_id_width[121];   //[n_clusters_id]
        Double_t        clus_id_tpos1[121];   //[n_clusters_id]
        Double_t        clus_id_tpos2[121];   //[n_clusters_id]
        Double_t        clus_id_lpos[121];   //[n_clusters_id]
        Double_t        clus_id_z[121];   //[n_clusters_id]
        Int_t           clus_id_view[121];   //[n_clusters_id]
        Int_t           clus_id_type[121];   //[n_clusters_id]
        Int_t           clus_id_hist[121];   //[n_clusters_id]
        Int_t           clus_id_subdet[121];   //[n_clusters_id]
        Double_t        clus_id_pe[121];   //[n_clusters_id]
        Double_t        clus_id_energy[121];   //[n_clusters_id]
        Double_t        clus_id_time[121];   //[n_clusters_id]
        Int_t           clus_id_time_slice[121];   //[n_clusters_id]
        Int_t           clus_id_size[121];   //[n_clusters_id]
        Int_t           clus_id_hits_idx[121][60];   //[n_clusters_id]
        Int_t           clus_id_usedFor[121];   //[n_clusters_id]
        Int_t           n_clusters_od;
        Int_t           clus_od_index[23];   //[n_clusters_od]
        Double_t        clus_od_z[23];   //[n_clusters_od]
        Int_t           clus_od_frame[23];   //[n_clusters_od]
        Int_t           clus_od_tower[23];   //[n_clusters_od]
        Int_t           clus_od_story[23];   //[n_clusters_od]
        Int_t           clus_od_hist[23];   //[n_clusters_od]
        Double_t        clus_od_pe[23];   //[n_clusters_od]
        Double_t        clus_od_energy[23];   //[n_clusters_od]
        Double_t        clus_od_time[23];   //[n_clusters_od]
        Int_t           clus_od_time_slice[23];   //[n_clusters_od]
        Int_t           clus_od_size[23];   //[n_clusters_od]
        Int_t           clus_od_hits_idx[23][60];   //[n_clusters_od]
        Int_t           n_blobs_id;
        Int_t           blob_id_idx[1];   //[n_blobs_id]
        Int_t           blob_id_subdet[1];   //[n_blobs_id]
        Int_t           blob_id_history[1];   //[n_blobs_id]
        Int_t           blob_id_size[1];   //[n_blobs_id]
        Int_t           blob_id_patrec[1];   //[n_blobs_id]
        Double_t        blob_id_e[1];   //[n_blobs_id]
        Double_t        blob_id_time[1];   //[n_blobs_id]
        Int_t           blob_id_time_slice[1];   //[n_blobs_id]
        Double_t        blob_id_startpoint_x[1];   //[n_blobs_id]
        Double_t        blob_id_startpoint_y[1];   //[n_blobs_id]
        Double_t        blob_id_startpoint_z[1];   //[n_blobs_id]
        Int_t           blob_id_clus_idx[1][1500];   //[n_blobs_id]
        Int_t           n_blobs_od;
        Int_t           blob_od_idx[1];   //[n_blobs_od]
        Int_t           blob_od_history[1];   //[n_blobs_od]
        Int_t           blob_od_size[1];   //[n_blobs_od]
        Int_t           blob_od_patrec[1];   //[n_blobs_od]
        Double_t        blob_od_e[1];   //[n_blobs_od]
        Double_t        blob_od_time[1];   //[n_blobs_od]
        Int_t           blob_od_time_slice[1];   //[n_blobs_od]
        Int_t           blob_od_clus_idx[1][1500];   //[n_blobs_od]
        Int_t           n_tracks;
        Int_t           trk_index[1];   //[n_tracks]
        Int_t           trk_type[1];   //[n_tracks]
        Int_t           trk_patrec[1];   //[n_tracks]
        Int_t           trk_time_slice[1];   //[n_tracks]
        Double_t        trk_vis_energy[1];   //[n_tracks]
        Double_t        trk_theta[1];   //[n_tracks]
        Double_t        trk_phi[1];   //[n_tracks]
        Int_t           trk_hits[1];   //[n_tracks]
        Int_t           trk_dof[1];   //[n_tracks]
        Double_t        trk_chi2perDof[1];   //[n_tracks]
        Double_t        trk_fitMass[1];   //[n_tracks]
        Int_t           trk_nodes[1];   //[n_tracks]
        Double_t        trk_node_X[1][300];   //[n_tracks]
        Double_t        trk_node_Y[1][300];   //[n_tracks]
        Double_t        trk_node_Z[1][300];   //[n_tracks]
        Double_t        trk_node_aX[1][300];   //[n_tracks]
        Double_t        trk_node_aY[1][300];   //[n_tracks]
        Double_t        trk_node_qOverP[1][300];   //[n_tracks]
        Double_t        trk_node_chi2[1][300];   //[n_tracks]
        Int_t           trk_node_cluster_idx[1][300];   //[n_tracks]
        Int_t           trk_usedFor[1];   //[n_tracks]
        Int_t           n_vertices;
        Int_t           vtx_time_slice[2];   //[n_vertices]
        Int_t           vtx_type[2];   //[n_vertices]
        Int_t           vtx_index[2];   //[n_vertices]
        Double_t        vtx_x[2];   //[n_vertices]
        Double_t        vtx_y[2];   //[n_vertices]
        Double_t        vtx_z[2];   //[n_vertices]
        Double_t        vtx_x_err[2];   //[n_vertices]
        Double_t        vtx_y_err[2];   //[n_vertices]
        Double_t        vtx_z_err[2];   //[n_vertices]
        Int_t           vtx_n_tracks[2];   //[n_vertices]
        Int_t           vtx_tracks_idx[2][20];   //[n_vertices]
        Int_t           n_minos_trk;
        Int_t           minos_run;
        Int_t           minos_subrun;
        Int_t           minos_snarl;
        Int_t           minos_trk_idx[1];   //[n_minos_trk]
        Int_t           minos_trk_minervatrk_idx[1];   //[n_minos_trk]
        Int_t           minos_trk_quality[1];   //[n_minos_trk]
        Double_t        minos_trk_pass[1];   //[n_minos_trk]
        Double_t        minos_trk_chi2[1];   //[n_minos_trk]
        Double_t        minos_trk_ndf[1];   //[n_minos_trk]
        Double_t        minos_trk_bave[1];   //[n_minos_trk]
        Double_t        minos_trk_range[1];   //[n_minos_trk]
        Int_t           minos_trk_con[1];   //[n_minos_trk]
        Double_t        minos_trk_p[1];   //[n_minos_trk]
        Double_t        minos_trk_prange[1];   //[n_minos_trk]
        Double_t        minos_trk_qp[1];   //[n_minos_trk]
        Double_t        minos_trk_eqp[1];   //[n_minos_trk]
        Int_t           minos_trk_vtxp[1];   //[n_minos_trk]
        Double_t        minos_trk_vtxu[1];   //[n_minos_trk]
        Double_t        minos_trk_vtxv[1];   //[n_minos_trk]
        Double_t        minos_trk_vtxx[1];   //[n_minos_trk]
        Double_t        minos_trk_vtxy[1];   //[n_minos_trk]
        Double_t        minos_trk_vtxz[1];   //[n_minos_trk]
        Double_t        minos_trk_vtxt[1];   //[n_minos_trk]
        Double_t        minos_trk_mvax[1];   //[n_minos_trk]
        Double_t        minos_trk_mvau[1];   //[n_minos_trk]
        Double_t        minos_trk_mvav[1];   //[n_minos_trk]
        Double_t        minos_trk_vtx_dxdz[1];   //[n_minos_trk]
        Double_t        minos_trk_vtx_dydz[1];   //[n_minos_trk]
        Double_t        minos_trk_vtx_dudz[1];   //[n_minos_trk]
        Double_t        minos_trk_vtx_dvdz[1];   //[n_minos_trk]
        Int_t           minos_trk_endp[1];   //[n_minos_trk]
        Double_t        minos_trk_endu[1];   //[n_minos_trk]
        Double_t        minos_trk_endv[1];   //[n_minos_trk]
        Double_t        minos_trk_endx[1];   //[n_minos_trk]
        Double_t        minos_trk_endy[1];   //[n_minos_trk]
        Double_t        minos_trk_endz[1];   //[n_minos_trk]
        Double_t        minos_trk_endt[1];   //[n_minos_trk]
        Int_t           minos_trk_ns[1];   //[n_minos_trk]
        Int_t           minos_trk_stp_fit[1][564];   //[n_minos_trk]
        Double_t        minos_trk_stp_u[1][564];   //[n_minos_trk]
        Double_t        minos_trk_stp_v[1][564];   //[n_minos_trk]
        Double_t        minos_trk_stp_x[1][564];   //[n_minos_trk]
        Double_t        minos_trk_stp_y[1][564];   //[n_minos_trk]
        Double_t        minos_trk_stp_z[1][564];   //[n_minos_trk]
        Double_t        minos_trk_stp_t[1][564];   //[n_minos_trk]
        Double_t        minos_trk_stp_meu[1][564];   //[n_minos_trk]
        Int_t           n_minos_stp;
        Int_t           minos_stp_plane[1];   //[n_minos_stp]
        Int_t           minos_stp_strip[1];   //[n_minos_stp]
        Int_t           minos_stp_view[1];   //[n_minos_stp]
        Double_t        minos_stp_tpos[1];   //[n_minos_stp]
        Double_t        minos_stp_time[1];   //[n_minos_stp]
        Double_t        minos_stp_z[1];   //[n_minos_stp]
        Double_t        minos_stp_ph[1];   //[n_minos_stp]
        Double_t        minos_stp_pe[1];   //[n_minos_stp]
        Int_t           minos_stp_trkidx[1];   //[n_minos_stp]
        Int_t           minos_sec;
        Int_t           minos_nanosec;
        Double_t        beam_pot;
        Double_t        beam_horncur;
        Double_t        beam_xpos;
        Double_t        beam_ypos;
        Double_t        beam_xwid;
        Double_t        beam_ywid;
        Double_t        beam_dt_nearest;
        Int_t           beam_dt_ok;
        Int_t           beam_pos_ok;
        Int_t           beam_wid_ok;
        Int_t           beam_tor_ok;
        Int_t           beam_horns_ok;
        Int_t           beam_numibeamdb_sec;
        Int_t           beam_numibeamdb_nanosec;
        Int_t           beam_is_good_beam_spill;
        Int_t           beam_is_bad_pot_data_spill;
        Int_t           beam_is_no_beam_spill;
        Int_t           beam_is_bad_data_spill;
        Int_t           beam_is_bad_prof_widx_data_spill;
        Int_t           beam_is_bad_prof_widy_data_spill;
        Int_t           beam_is_bad_xpos_data_spill;
        Int_t           beam_is_bad_ypos_data_spill;
        Int_t           beam_is_bad_horncur_data_spill;
        Int_t           beam_is_bad_nearesttime_spill;
        Int_t           beam_is_bad_beam_spill;
        Int_t           beam_is_bad_pot_spill;
        Int_t           beam_is_bad_xpos_spill;
        Int_t           beam_is_bad_ypos_spill;
        Int_t           beam_is_bad_beamsize_spill;
        Int_t           beam_is_bad_prof_widx_spill;
        Int_t           beam_is_bad_prof_widy_spill;
        Int_t           beam_is_bad_horncur_spill;
        Int_t           beam_is_target_out_spill;
        Int_t           beam_is_bad_beamtype_spill;
        Int_t           beam_is_bad_beam_frac_on_tgt_spill;
        Double_t        beam_tor101;
        Double_t        beam_tr101d;
        Double_t        beam_tortgt;
        Double_t        beam_trtgtd;
        Int_t           mc_run;
        Int_t           mc_subrun;
        Int_t           mc_spill;
        Int_t           n_total_interactions;
        Int_t           mc_MIState;
        Double_t        mc_pot;
        Int_t           mc_beamConfig;
        Int_t           n_interactions;
        Int_t           mc_int_processType[1];   //[n_interactions]
        Int_t           mc_int_nevSpill[1];   //[n_interactions]
        Int_t           mc_int_nevFile[1];   //[n_interactions]
        Int_t           mc_int_channel[1];   //[n_interactions]
        Int_t           mc_int_current[1];   //[n_interactions]
        Int_t           mc_int_charm[1];   //[n_interactions]
        Double_t        mc_int_weight[1];   //[n_interactions]
        Double_t        mc_int_xSection[1];   //[n_interactions]
        Int_t           mc_int_incomingPDG[1];   //[n_interactions]
        Int_t           mc_int_tgtNucleus[1];   //[n_interactions]
        Int_t           mc_int_tgtNucleon[1];   //[n_interactions]
        Int_t           mc_int_targetZ[1];   //[n_interactions]
        Int_t           mc_int_targetA[1];   //[n_interactions]
        Int_t           mc_int_hitQuark[1];   //[n_interactions]
        Int_t           mc_int_seaQuark[1];   //[n_interactions]
        Int_t           mc_int_resID[1];   //[n_interactions]
        Int_t           mc_int_FSLepton[1];   //[n_interactions]
        Double_t        mc_int_incomingE[1];   //[n_interactions]
        Double_t        mc_int_bjorkenX[1];   //[n_interactions]
        Double_t        mc_int_bjorkenY[1];   //[n_interactions]
        Double_t        mc_int_QSquared[1];   //[n_interactions]
        Double_t        mc_int_nucleonT[1];   //[n_interactions]
        Double_t        mc_int_W[1];   //[n_interactions]
        Int_t           mc_int_nFSParticles[1];   //[n_interactions]
        Double_t        mc_int_vtx[1][4];   //[n_interactions]
        Double_t        mc_int_incoming4p[1][4];   //[n_interactions]
        Double_t        mc_int_tgtNucleon4p[1][4];   //[n_interactions]
        Double_t        mc_int_FSLepton4p[1][4];   //[n_interactions]
        Int_t           mc_int_FSPdg[1][30];   //[n_interactions]
        Double_t        mc_int_FSParticlesPx[1][30];   //[n_interactions]
        Double_t        mc_int_FSParticlesPy[1][30];   //[n_interactions]
        Double_t        mc_int_FSParticlesPz[1][30];   //[n_interactions]
        Double_t        mc_int_FSParticlesE[1][30];   //[n_interactions]
        Double_t        mc_flux_proton_P[1][4];   //[n_interactions]
        Double_t        mc_flux_proton_X[1][3];   //[n_interactions]
        Int_t           mc_flux_parent_PDG[1];   //[n_interactions]
        Double_t        mc_flux_parent_prod4P[1][4];   //[n_interactions]
        Double_t        mc_flux_parent_prodPos[1][3];   //[n_interactions]
        Double_t        mc_flux_parent_decay4P[1][4];   //[n_interactions]
        Double_t        mc_flux_parent_decayPos[1][3];   //[n_interactions]
        Int_t           mc_flux_parent_generation[1];   //[n_interactions]
        Int_t           mc_flux_parent_decayMode[1];   //[n_interactions]
        Int_t           mc_flux_secondary_PDG[1];   //[n_interactions]
        Int_t           n_mc_trajectories;
        Bool_t          mc_traj_overflow;
        Double_t        mc_traj_strlength[6];   //[n_mc_trajectories]
        Double_t        mc_traj_curvlength[6];   //[n_mc_trajectories]
        Char_t          mc_traj_leaving[6];   //[n_mc_trajectories]
        Int_t           mc_traj_trkid[6];   //[n_mc_trajectories]
        Int_t           mc_traj_parentid[6];   //[n_mc_trajectories]
        Int_t           mc_traj_pdg[6];   //[n_mc_trajectories]
        Double_t        mc_traj_hit_e[6];   //[n_mc_trajectories]
        Int_t           mc_traj_npoints[6];   //[n_mc_trajectories]
        Double_t        mc_traj_point_x[6][5];   //[n_mc_trajectories]
        Double_t        mc_traj_point_y[6][5];   //[n_mc_trajectories]
        Double_t        mc_traj_point_z[6][5];   //[n_mc_trajectories]
        Double_t        mc_traj_point_t[6][5];   //[n_mc_trajectories]
        Double_t        mc_traj_point_px[6][5];   //[n_mc_trajectories]
        Double_t        mc_traj_point_py[6][5];   //[n_mc_trajectories]
        Double_t        mc_traj_point_pz[6][5];   //[n_mc_trajectories]
        Double_t        mc_traj_point_E[6][5];   //[n_mc_trajectories]
        Int_t           n_mc_id_digits;
        Int_t           mc_id_strip[138];   //[n_mc_id_digits]
        Int_t           mc_id_plane[138];   //[n_mc_id_digits]
        Int_t           mc_id_module[138];   //[n_mc_id_digits]
        Int_t           mc_id_view[138];   //[n_mc_id_digits]
        Double_t        mc_id_pe[138];   //[n_mc_id_digits]
        Double_t        mc_id_time[138];   //[n_mc_id_digits]
        Double_t        mc_id_dE[138];   //[n_mc_id_digits]
        Int_t           mc_id_nmchit[138];   //[n_mc_id_digits]
        Double_t        mc_id_mchit_x[138][2];   //[n_mc_id_digits]
        Double_t        mc_id_mchit_y[138][2];   //[n_mc_id_digits]
        Double_t        mc_id_mchit_z[138][2];   //[n_mc_id_digits]
        Double_t        mc_id_mchit_t[138][2];   //[n_mc_id_digits]
        Int_t           mc_id_mchit_trkid[138][2];   //[n_mc_id_digits]
        Double_t        mc_id_mchit_p[138][2];   //[n_mc_id_digits]
        Double_t        mc_id_mchit_dE[138][2];   //[n_mc_id_digits]
        Double_t        mc_id_mchit_dL[138][2];   //[n_mc_id_digits]
        Int_t           n_mc_od_digits;
        Int_t           mc_od_frame[23];   //[n_mc_od_digits]
        Int_t           mc_od_tower[23];   //[n_mc_od_digits]
        Int_t           mc_od_story[23];   //[n_mc_od_digits]
        Int_t           mc_od_bar[23];   //[n_mc_od_digits]
        Double_t        mc_od_pe[23];   //[n_mc_od_digits]
        Double_t        mc_od_time[23];   //[n_mc_od_digits]
        Double_t        mc_od_dE[23];   //[n_mc_od_digits]
        Int_t           mc_od_nmchit[23];   //[n_mc_od_digits]
        Double_t        mc_od_mchit_x[23][2];   //[n_mc_od_digits]
        Double_t        mc_od_mchit_y[23][2];   //[n_mc_od_digits]
        Double_t        mc_od_mchit_z[23][2];   //[n_mc_od_digits]
        Double_t        mc_od_mchit_t[23][2];   //[n_mc_od_digits]
        Int_t           mc_od_mchit_trkid[23][2];   //[n_mc_od_digits]
        Double_t        mc_od_mchit_p[23][2];   //[n_mc_od_digits]
        Double_t        mc_od_mchit_dE[23][2];   //[n_mc_od_digits]
        Int_t           n_mc_veto_digits;
        Int_t           mc_veto_wall[1];   //[n_mc_veto_digits]
        Int_t           mc_veto_paddle[1];   //[n_mc_veto_digits]
        Double_t        mc_veto_pe[1];   //[n_mc_veto_digits]
        Double_t        mc_veto_time[1];   //[n_mc_veto_digits]
        Double_t        mc_veto_dE[1];   //[n_mc_veto_digits]
        Int_t           mc_veto_nmchit[1];   //[n_mc_veto_digits]
        Double_t        mc_veto_mchit_x[1][2];   //[n_mc_veto_digits]
        Double_t        mc_veto_mchit_y[1][2];   //[n_mc_veto_digits]
        Double_t        mc_veto_mchit_z[1][2];   //[n_mc_veto_digits]
        Double_t        mc_veto_mchit_t[1][2];   //[n_mc_veto_digits]
        Int_t           mc_veto_mchit_trkid[1][2];   //[n_mc_veto_digits]
        Double_t        mc_veto_mchit_p[1][2];   //[n_mc_veto_digits]
        Double_t        mc_veto_mchit_dE[1][2];   //[n_mc_veto_digits]

        // List of branches
        TBranch        *b_fmwk_v;   //!
        TBranch        *b_fwmk_r;   //!
        TBranch        *b_fwmk_p;   //!
        TBranch        *b_n_slices;   //!
        TBranch        *b_ev_detector;   //!
        TBranch        *b_ev_det_config;   //!
        TBranch        *b_ev_run;   //!
        TBranch        *b_ev_sub_run;   //!
        TBranch        *b_ev_trigger_type;   //!
        TBranch        *b_ev_cal_settings;   //!
        TBranch        *b_ev_gl_gate;   //!
        TBranch        *b_ev_gate;   //!
        TBranch        *b_ev_gps_time_sec;   //!
        TBranch        *b_ev_gps_time_usec;   //!
        TBranch        *b_ev_readout_time;   //!
        TBranch        *b_ev_errors;   //!
        TBranch        *b_ev_nADC_Frames;   //!
        TBranch        *b_ev_nDisc_Frames;   //!
        TBranch        *b_ev_nFPGA_Frames;   //!
        TBranch        *b_n_febs;   //!
        TBranch        *b_feb_id;   //!
        TBranch        *b_feb_hv_on;   //!
        TBranch        *b_feb_hv_targ;   //!
        TBranch        *b_feb_hv_act;   //!
        TBranch        *b_feb_hv_per_auto;   //!
        TBranch        *b_feb_temperature;   //!
        TBranch        *b_feb_gate_time_stamp;   //!
        TBranch        *b_n_rawhits;   //!
        TBranch        *b_hit_feb_id;   //!
        TBranch        *b_hit_flags;   //!
        TBranch        *b_hit_channel_id;   //!
        TBranch        *b_hit_index;   //!
        TBranch        *b_hit_location;   //!
        TBranch        *b_hit_is_mc;   //!
        TBranch        *b_hit_num;   //!
        TBranch        *b_hit_pixel;   //!
        TBranch        *b_hit_board;   //!
        TBranch        *b_hit_chain;   //!
        TBranch        *b_hit_croc;   //!
        TBranch        *b_hit_crate;   //!
        TBranch        *b_hit_link;   //!
        TBranch        *b_hit_disc_fired;   //!
        TBranch        *b_hit_sys_ticks;   //!
        TBranch        *b_hit_delay_ticks;   //!
        TBranch        *b_hit_quarter_ticks;   //!
        TBranch        *b_hit_qlo;   //!
        TBranch        *b_hit_qmed;   //!
        TBranch        *b_hit_qhi;   //!
        TBranch        *b_n_idhits;   //!
        TBranch        *b_hits_id_per_mod;   //!
        TBranch        *b_hit_strip;   //!
        TBranch        *b_hit_plane;   //!
        TBranch        *b_hit_module;   //!
        TBranch        *b_hit_view;   //!
        TBranch        *b_n_odhits;   //!
        TBranch        *b_hits_od_per_mod;   //!
        TBranch        *b_hit_bar;   //!
        TBranch        *b_hit_story;   //!
        TBranch        *b_hit_tower;   //!
        TBranch        *b_hit_frame;   //!
        TBranch        *b_n_vetohits;   //!
        TBranch        *b_hit_wall;   //!
        TBranch        *b_hit_paddle;   //!
        TBranch        *b_hit_pmt;   //!
        TBranch        *b_hit_q;   //!
        TBranch        *b_hit_pe;   //!
        TBranch        *b_hit_norm_energy;   //!
        TBranch        *b_hit_time_raw;   //!
        TBranch        *b_hit_time;   //!
        TBranch        *b_hit_time_slice;   //!
        TBranch        *b_hits_total_pe;   //!
        TBranch        *b_hit_user_color;   //!
        TBranch        *b_n_clusters_id;   //!
        TBranch        *b_clus_id_index;   //!
        TBranch        *b_clus_id_strip;   //!
        TBranch        *b_clus_id_plane;   //!
        TBranch        *b_clus_id_module;   //!
        TBranch        *b_clus_id_coord;   //!
        TBranch        *b_clus_id_coordErr;   //!
        TBranch        *b_clus_id_width;   //!
        TBranch        *b_clus_id_tpos1;   //!
        TBranch        *b_clus_id_tpos2;   //!
        TBranch        *b_clus_id_lpos;   //!
        TBranch        *b_clus_id_z;   //!
        TBranch        *b_clus_id_view;   //!
        TBranch        *b_clus_id_type;   //!
        TBranch        *b_clus_id_hist;   //!
        TBranch        *b_clus_id_subdet;   //!
        TBranch        *b_clus_id_pe;   //!
        TBranch        *b_clus_id_energy;   //!
        TBranch        *b_clus_id_time;   //!
        TBranch        *b_clus_id_time_slice;   //!
        TBranch        *b_clus_id_size;   //!
        TBranch        *b_clus_id_hits_idx;   //!
        TBranch        *b_clus_id_usedFor;   //!
        TBranch        *b_n_clusters_od;   //!
        TBranch        *b_clus_od_index;   //!
        TBranch        *b_clus_od_z;   //!
        TBranch        *b_clus_od_frame;   //!
        TBranch        *b_clus_od_tower;   //!
        TBranch        *b_clus_od_story;   //!
        TBranch        *b_clus_od_hist;   //!
        TBranch        *b_clus_od_pe;   //!
        TBranch        *b_clus_od_energy;   //!
        TBranch        *b_clus_od_time;   //!
        TBranch        *b_clus_od_time_slice;   //!
        TBranch        *b_clus_od_size;   //!
        TBranch        *b_clus_od_hits_idx;   //!
        TBranch        *b_n_blobs_id;   //!
        TBranch        *b_blob_id_idx;   //!
        TBranch        *b_blob_id_subdet;   //!
        TBranch        *b_blob_id_history;   //!
        TBranch        *b_blob_id_size;   //!
        TBranch        *b_blob_id_patrec;   //!
        TBranch        *b_blob_id_e;   //!
        TBranch        *b_blob_id_time;   //!
        TBranch        *b_blob_id_time_slice;   //!
        TBranch        *b_blob_id_startpoint_x;   //!
        TBranch        *b_blob_id_startpoint_y;   //!
        TBranch        *b_blob_id_startpoint_z;   //!
        TBranch        *b_blob_id_clus_idx;   //!
        TBranch        *b_n_blobs_od;   //!
        TBranch        *b_blob_od_idx;   //!
        TBranch        *b_blob_od_history;   //!
        TBranch        *b_blob_od_size;   //!
        TBranch        *b_blob_od_patrec;   //!
        TBranch        *b_blob_od_e;   //!
        TBranch        *b_blob_od_time;   //!
        TBranch        *b_blob_od_time_slice;   //!
        TBranch        *b_blob_od_clus_idx;   //!
        TBranch        *b_n_tracks;   //!
        TBranch        *b_trk_index;   //!
        TBranch        *b_trk_type;   //!
        TBranch        *b_trk_patrec;   //!
        TBranch        *b_trk_time_slice;   //!
        TBranch        *b_trk_vis_energy;   //!
        TBranch        *b_trk_theta;   //!
        TBranch        *b_trk_phi;   //!
        TBranch        *b_trk_hits;   //!
        TBranch        *b_trk_dof;   //!
        TBranch        *b_trk_chi2perDof;   //!
        TBranch        *b_trk_fitMass;   //!
        TBranch        *b_trk_nodes;   //!
        TBranch        *b_trk_node_X;   //!
        TBranch        *b_trk_node_Y;   //!
        TBranch        *b_trk_node_Z;   //!
        TBranch        *b_trk_node_aX;   //!
        TBranch        *b_trk_node_aY;   //!
        TBranch        *b_trk_node_qOverP;   //!
        TBranch        *b_trk_node_chi2;   //!
        TBranch        *b_trk_node_cluster_idx;   //!
        TBranch        *b_trk_usedFor;   //!
        TBranch        *b_n_vertices;   //!
        TBranch        *b_vtx_time_slice;   //!
        TBranch        *b_vtx_type;   //!
        TBranch        *b_vtx_index;   //!
        TBranch        *b_vtx_x;   //!
        TBranch        *b_vtx_y;   //!
        TBranch        *b_vtx_z;   //!
        TBranch        *b_vtx_x_err;   //!
        TBranch        *b_vtx_y_err;   //!
        TBranch        *b_vtx_z_err;   //!
        TBranch        *b_vtx_n_tracks;   //!
        TBranch        *b_vtx_tracks_idx;   //!
        TBranch        *b_n_minos_trk;   //!
        TBranch        *b_minos_run;   //!
        TBranch        *b_minos_subrun;   //!
        TBranch        *b_minos_snarl;   //!
        TBranch        *b_minos_trk_idx;   //!
        TBranch        *b_minos_trk_minervatrk_idx;   //!
        TBranch        *b_minos_trk_quality;   //!
        TBranch        *b_minos_trk_pass;   //!
        TBranch        *b_minos_trk_chi2;   //!
        TBranch        *b_minos_trk_ndf;   //!
        TBranch        *b_minos_trk_bave;   //!
        TBranch        *b_minos_trk_range;   //!
        TBranch        *b_minos_trk_con;   //!
        TBranch        *b_minos_trk_p;   //!
        TBranch        *b_minos_trk_prange;   //!
        TBranch        *b_minos_trk_qp;   //!
        TBranch        *b_minos_trk_eqp;   //!
        TBranch        *b_minos_trk_vtxp;   //!
        TBranch        *b_minos_trk_vtxu;   //!
        TBranch        *b_minos_trk_vtxv;   //!
        TBranch        *b_minos_trk_vtxx;   //!
        TBranch        *b_minos_trk_vtxy;   //!
        TBranch        *b_minos_trk_vtxz;   //!
        TBranch        *b_minos_trk_vtxt;   //!
        TBranch        *b_minos_trk_mvax;   //!
        TBranch        *b_minos_trk_mvau;   //!
        TBranch        *b_minos_trk_mvav;   //!
        TBranch        *b_minos_trk_vtx_dxdz;   //!
        TBranch        *b_minos_trk_vtx_dydz;   //!
        TBranch        *b_minos_trk_vtx_dudz;   //!
        TBranch        *b_minos_trk_vtx_dvdz;   //!
        TBranch        *b_minos_trk_endp;   //!
        TBranch        *b_minos_trk_endu;   //!
        TBranch        *b_minos_trk_endv;   //!
        TBranch        *b_minos_trk_endx;   //!
        TBranch        *b_minos_trk_endy;   //!
        TBranch        *b_minos_trk_endz;   //!
        TBranch        *b_minos_trk_endt;   //!
        TBranch        *b_minos_trk_ns;   //!
        TBranch        *b_minos_trk_stp_fit;   //!
        TBranch        *b_minos_trk_stp_u;   //!
        TBranch        *b_minos_trk_stp_v;   //!
        TBranch        *b_minos_trk_stp_x;   //!
        TBranch        *b_minos_trk_stp_y;   //!
        TBranch        *b_minos_trk_stp_z;   //!
        TBranch        *b_minos_trk_stp_t;   //!
        TBranch        *b_minos_trk_stp_meu;   //!
        TBranch        *b_n_minos_stp;   //!
        TBranch        *b_minos_stp_plane;   //!
        TBranch        *b_minos_stp_strip;   //!
        TBranch        *b_minos_stp_view;   //!
        TBranch        *b_minos_stp_tpos;   //!
        TBranch        *b_minos_stp_time;   //!
        TBranch        *b_minos_stp_z;   //!
        TBranch        *b_minos_stp_ph;   //!
        TBranch        *b_minos_stp_pe;   //!
        TBranch        *b_minos_stp_trkidx;   //!
        TBranch        *b_minos_sec;   //!
        TBranch        *b_minos_nanosec;   //!
        TBranch        *b_beam_pot;   //!
        TBranch        *b_beam_horncur;   //!
        TBranch        *b_beamxpos;   //!
        TBranch        *b_beam_ypos;   //!
        TBranch        *b_beam_xwid;   //!
        TBranch        *b_beam_ywid;   //!
        TBranch        *b_beam_dt_nearest;   //!
        TBranch        *b_beam_dt_ok;   //!
        TBranch        *b_beam_pos_ok;   //!
        TBranch        *b_beam_wid_ok;   //!
        TBranch        *b_beam_tor_ok;   //!
        TBranch        *b_beam_horns_ok;   //!
        TBranch        *b_beam_numibeamdb_sec;   //!
        TBranch        *b_beam_numibeamdb_nanosec;   //!
        TBranch        *b_beam_is_good_beam_spill;   //!
        TBranch        *b_beam_is_bad_pot_data_spill;   //!
        TBranch        *b_beam_is_no_beam_spill;   //!
        TBranch        *b_beam_is_bad_data_spill;   //!
        TBranch        *b_beam_is_bad_prof_widx_data_spill;   //!
        TBranch        *b_beam_is_bad_prof_widy_data_spill;   //!
        TBranch        *b_beam_is_bad_xpos_data_spill;   //!
        TBranch        *b_beam_is_bad_ypos_data_spill;   //!
        TBranch        *b_beam_is_bad_horncur_data_spill;   //!
        TBranch        *b_beam_is_bad_nearesttime_spill;   //!
        TBranch        *b_beam_is_bad_beam_spill;   //!
        TBranch        *b_beam_is_bad_pot_spill;   //!
        TBranch        *b_beam_is_bad_xpos_spill;   //!
        TBranch        *b_beam_is_bad_ypos_spill;   //!
        TBranch        *b_beam_is_bad_beamsize_spill;   //!
        TBranch        *b_beam_is_bad_prof_widx_spill;   //!
        TBranch        *b_beam_is_bad_prof_widy_spill;   //!
        TBranch        *b_beam_is_bad_horncur_spill;   //!
        TBranch        *b_beam_is_target_out_spill;   //!
        TBranch        *b_beam_is_bad_beamtype_spill;   //!
        TBranch        *b_beam_is_bad_beam_frac_on_tgt_spill;   //!
        TBranch        *b_beam_tor101;   //!
        TBranch        *b_beam_tr101d;   //!
        TBranch        *b_beam_tortgt;   //!
        TBranch        *b_beam_trtgtd;   //!
        TBranch        *b_mc_run;   //!
        TBranch        *b_mc_subrun;   //!
        TBranch        *b_mc_spill;   //!
        TBranch        *b_n_total_interactions;   //!
        TBranch        *b_mc_MIState;   //!
        TBranch        *b_mc_pot;   //!
        TBranch        *b_mc_beamConfig;   //!
        TBranch        *b_n_interactions;   //!
        TBranch        *b_mc_int_processType;   //!
        TBranch        *b_mc_int_nevSpill;   //!
        TBranch        *b_mc_int_nevFile;   //!
        TBranch        *b_mc_int_channel;   //!
        TBranch        *b_mc_int_current;   //!
        TBranch        *b_mc_int_charm;   //!
        TBranch        *b_mc_int_weight;   //!
        TBranch        *b_mc_int_xSection;   //!
        TBranch        *b_mc_int_incomingPDG;   //!
        TBranch        *b_mc_int_tgtNucleus;   //!
        TBranch        *b_mc_int_tgtNucleon;   //!
        TBranch        *b_mc_int_targetZ;   //!
        TBranch        *b_mc_int_targetA;   //!
        TBranch        *b_mc_int_hitQuark;   //!
        TBranch        *b_mc_int_seaQuark;   //!
        TBranch        *b_mc_int_resID;   //!
        TBranch        *b_mc_int_FSLepton;   //!
        TBranch        *b_mc_int_incomingE;   //!
        TBranch        *b_mc_int_bjorkenX;   //!
        TBranch        *b_mc_int_bjorkenY;   //!
        TBranch        *b_mc_int_QSquared;   //!
        TBranch        *b_mc_int_nucleonT;   //!
        TBranch        *b_mc_int_W;   //!
        TBranch        *b_mc_int_nFSParticles;   //!
        TBranch        *b_mc_int_vtx;   //!
        TBranch        *b_mc_int_incoming4p;   //!
        TBranch        *b_mc_int_tgtNucleon4p;   //!
        TBranch        *b_mc_int_FSLepton4p;   //!
        TBranch        *b_mc_int_FSPdg;   //!
        TBranch        *b_mc_int_FSParticlesPx;   //!
        TBranch        *b_mc_int_FSParticlesPy;   //!
        TBranch        *b_mc_int_FSParticlesPz;   //!
        TBranch        *b_mc_int_FSParticlesE;   //!
        TBranch        *b_mc_flux_proton_P;   //!
        TBranch        *b_mc_flux_proton_X;   //!
        TBranch        *b_mc_flux_parent_PDG;   //!
        TBranch        *b_mc_flux_parent_prod4P;   //!
        TBranch        *b_mc_flux_parent_prodPos;   //!
        TBranch        *b_mc_flux_parent_decay4P;   //!
        TBranch        *b_mc_flux_parent_decayPos;   //!
        TBranch        *b_mc_flux_parent_generation;   //!
        TBranch        *b_mc_flux_parent_decayMode;   //!
        TBranch        *b_mc_flux_secondary_PDG;   //!
        TBranch        *b_n_mc_trajectories;   //!
        TBranch        *b_mc_traj_overflow;   //!
        TBranch        *b_mc_traj_strlength;   //!
        TBranch        *b_mc_traj_curvlength;   //!
        TBranch        *b_mc_traj_leaving;   //!
        TBranch        *b_mc_traj_trkid;   //!
        TBranch        *b_mc_traj_parentid;   //!
        TBranch        *b_mc_traj_pdg;   //!
        TBranch        *b_mc_traj_hit_e;   //!
        TBranch        *b_mc_traj_npoints;   //!
        TBranch        *b_mc_traj_point_x;   //!
        TBranch        *b_mc_traj_point_y;   //!
        TBranch        *b_mc_traj_point_z;   //!
        TBranch        *b_mc_traj_point_t;   //!
        TBranch        *b_mc_traj_point_px;   //!
        TBranch        *b_mc_traj_point_py;   //!
        TBranch        *b_mc_traj_point_pz;   //!
        TBranch        *b_mc_traj_point_E;   //!
        TBranch        *b_n_mc_id_digits;   //!
        TBranch        *b_mc_id_strip;   //!
        TBranch        *b_mc_id_plane;   //!
        TBranch        *b_mc_id_module;   //!
        TBranch        *b_mc_id_view;   //!
        TBranch        *b_mc_id_pe;   //!
        TBranch        *b_mc_id_time;   //!
        TBranch        *b_mc_id_dE;   //!
        TBranch        *b_mc_id_nmchit;   //!
        TBranch        *b_mc_id_mchit_x;   //!
        TBranch        *b_mc_id_mchit_y;   //!
        TBranch        *b_mc_id_mchit_z;   //!
        TBranch        *b_mc_id_mchit_t;   //!
        TBranch        *b_mc_id_mchit_trkid;   //!
        TBranch        *b_mc_id_mchit_p;   //!
        TBranch        *b_mc_id_mchit_dE;   //!
        TBranch        *b_mc_id_mchit_dL;   //!
        TBranch        *b_n_mc_od_digits;   //!
        TBranch        *b_mc_od_frame;   //!
        TBranch        *b_mc_od_tower;   //!
        TBranch        *b_mc_od_story;   //!
        TBranch        *b_mc_od_bar;   //!
        TBranch        *b_mc_od_pe;   //!
        TBranch        *b_mc_od_time;   //!
        TBranch        *b_mc_od_dE;   //!
        TBranch        *b_mc_od_nmchit;   //!
        TBranch        *b_mc_od_mchit_x;   //!
        TBranch        *b_mc_od_mchit_y;   //!
        TBranch        *b_mc_od_mchit_z;   //!
        TBranch        *b_mc_od_mchit_t;   //!
        TBranch        *b_mc_od_mchit_trkid;   //!
        TBranch        *b_mc_od_mchit_p;   //!
        TBranch        *b_mc_od_mchit_dE;   //!
        TBranch        *b_n_mc_veto_digits;   //!
        TBranch        *b_mc_veto_wall;   //!
        TBranch        *b_mc_veto_paddle;   //!
        TBranch        *b_mc_veto_pe;   //!
        TBranch        *b_mc_veto_time;   //!
        TBranch        *b_mc_veto_dE;   //!
        TBranch        *b_mc_veto_nmchit;   //!
        TBranch        *b_mc_veto_mchit_x;   //!
        TBranch        *b_mc_veto_mchit_y;   //!
        TBranch        *b_mc_veto_mchit_z;   //!
        TBranch        *b_mc_veto_mchit_t;   //!
        TBranch        *b_mc_veto_mchit_trkid;   //!
        TBranch        *b_mc_veto_mchit_p;   //!
        TBranch        *b_mc_veto_mchit_dE;   //!

        PC_DST();
        virtual ~PC_DST();
        virtual Int_t    GetEntry(Long64_t entry);
        virtual Long64_t LoadTree(Long64_t entry);
        virtual void     Init(string playlist, TChain *fChain);
        virtual void     Loop(string playlist);
};

#endif
