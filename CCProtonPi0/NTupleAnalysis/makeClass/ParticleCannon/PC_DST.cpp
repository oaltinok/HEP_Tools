#ifndef PC_DST_cpp
#define PC_DST_cpp

#include "PC_DST.h"

using namespace std;

void PC_DST::Loop(string playlist)
{
    // Create Chain and Initialize
    TChain* fChain = new TChain("minerva");
    Init(playlist, fChain);
    if (!fChain) return;
    if (fChain == 0) return;

    Long64_t nbytes = 0, nb = 0;
    Long64_t nentries = fChain->GetEntriesFast();

    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = fChain->GetEntry(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;

        if (ientry == 0) {
            cout<<"\tGetEntry failure "<<jentry<<endl;
            break;
        }

        // Progress Message on Terminal
        if (jentry%25 == 0) cout<<"\tEntry "<<jentry<<endl;

        double reco_E = CalculateEnergy();
        double true_E = mc_int_FSParticlesE[0][0]; 
        double error_E = Data_Functions::getError(true_E, reco_E);

        error->Fill(error_E);
    }

    writeHistograms();
}

double PC_DST::GetEvis(int min_module, int max_module)
{
    double totalE = 0.0;

    for (int i = 0;  i < n_rawhits; i++){
        if (hit_norm_energy[i] < 0.09) continue;
        if (hit_module[i] >= min_module && hit_module[i] <= max_module){
            totalE += hit_norm_energy[i];
        }
    }

    return totalE;
}

double PC_DST::CalculateEnergy()
{
    const double alpha = 1.326;
    const double kE = 2.341;
    const double kH = 9.54;

    double total_evis_target = GetEvis(-5,21);
    double total_evis_trkr = GetEvis(23,84);
    double total_evis_ecal = GetEvis(85,94);
    double total_evis_hcal = GetEvis(95,114);

    double energy = alpha*(total_evis_trkr + kE*total_evis_ecal + kH*total_evis_hcal);

    return energy;
}

void PC_DST::writeHistograms()
{
    std::cout<<">> Writing "<<rootDir<<std::endl;
    f->cd();

    error->Write();

    f->Close();
}

PC_DST::PC_DST()
{
    // Open Output ROOT File
    rootDir = Folder_List::rootOut + Folder_List::ParticleCannon + "PC_Test.root";
 
    cout<<"\tRoot File: "<<rootDir<<endl;
    f = new TFile(rootDir.c_str(),"RECREATE");

    InitHistograms();
}

void PC_DST::InitHistograms()
{
    error = new TH1D("error","Energy Resolution",20,-2,2);
    error->GetXaxis()->SetTitle("(E_{Reco}-E_{True})/E_{True}");
    error->GetYaxis()->SetTitle("N(Events)");
}


PC_DST::~PC_DST()
{
    if (!fChain) return;
    delete fChain->GetCurrentFile();
}

Int_t PC_DST::GetEntry(Long64_t entry)
{
    // Read contents of entry.
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}
Long64_t PC_DST::LoadTree(Long64_t entry)
{
    // Set the environment to read one entry
    if (!fChain) return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0) return centry;
    if (fChain->GetTreeNumber() != fCurrent) {
        fCurrent = fChain->GetTreeNumber();
    }
    return centry;
}

void PC_DST::Init(string playlist, TChain *fChain)
{
    ifstream input_pl(playlist.c_str());
    string filename;

    if( !input_pl.is_open() ){
        std::cerr<<"Cannot open Playlist File!"<<std::endl;
        exit(1);
    }else{
        std::cout<<"Reading Playlist: "<<playlist.c_str()<<std::endl;
    }

    while (input_pl) {
        input_pl>>filename;

        if (!input_pl) break;

        if (filename[0] != '/') break;

        fChain->Add( filename.c_str() );
        //cout<<" Added "<<filename.c_str()<<endl;
    }

    fChain->SetMakeClass(1);
    fChain->SetBranchAddress("fmwk_v", &fmwk_v, &b_fmwk_v);
    fChain->SetBranchAddress("fmwk_r", &fmwk_r, &b_fwmk_r);
    fChain->SetBranchAddress("fmwk_p", &fmwk_p, &b_fwmk_p);
    fChain->SetBranchAddress("n_slices", &n_slices, &b_n_slices);
    fChain->SetBranchAddress("ev_detector", &ev_detector, &b_ev_detector);
    fChain->SetBranchAddress("ev_det_config", &ev_det_config, &b_ev_det_config);
    fChain->SetBranchAddress("ev_run", &ev_run, &b_ev_run);
    fChain->SetBranchAddress("ev_sub_run", &ev_sub_run, &b_ev_sub_run);
    fChain->SetBranchAddress("ev_trigger_type", &ev_trigger_type, &b_ev_trigger_type);
    fChain->SetBranchAddress("ev_cal_settings", &ev_cal_settings, &b_ev_cal_settings);
    fChain->SetBranchAddress("ev_gl_gate", &ev_gl_gate, &b_ev_gl_gate);
    fChain->SetBranchAddress("ev_gate", &ev_gate, &b_ev_gate);
    fChain->SetBranchAddress("ev_gps_time_sec", &ev_gps_time_sec, &b_ev_gps_time_sec);
    fChain->SetBranchAddress("ev_gps_time_usec", &ev_gps_time_usec, &b_ev_gps_time_usec);
    fChain->SetBranchAddress("ev_readout_time", &ev_readout_time, &b_ev_readout_time);
    fChain->SetBranchAddress("ev_errors", &ev_errors, &b_ev_errors);
    fChain->SetBranchAddress("ev_nADC_Frames", &ev_nADC_Frames, &b_ev_nADC_Frames);
    fChain->SetBranchAddress("ev_nDisc_Frames", &ev_nDisc_Frames, &b_ev_nDisc_Frames);
    fChain->SetBranchAddress("ev_nFPGA_Frames", &ev_nFPGA_Frames, &b_ev_nFPGA_Frames);
    fChain->SetBranchAddress("n_febs", &n_febs, &b_n_febs);
    fChain->SetBranchAddress("feb_id", feb_id, &b_feb_id);
    fChain->SetBranchAddress("feb_hv_on", feb_hv_on, &b_feb_hv_on);
    fChain->SetBranchAddress("feb_hv_targ", feb_hv_targ, &b_feb_hv_targ);
    fChain->SetBranchAddress("feb_hv_act", feb_hv_act, &b_feb_hv_act);
    fChain->SetBranchAddress("feb_hv_per_auto", feb_hv_per_auto, &b_feb_hv_per_auto);
    fChain->SetBranchAddress("feb_temperature", feb_temperature, &b_feb_temperature);
    fChain->SetBranchAddress("feb_gate_time_stamp", feb_gate_time_stamp, &b_feb_gate_time_stamp);
    fChain->SetBranchAddress("n_rawhits", &n_rawhits, &b_n_rawhits);
    fChain->SetBranchAddress("hit_feb_id", hit_feb_id, &b_hit_feb_id);
    fChain->SetBranchAddress("hit_flags", hit_flags, &b_hit_flags);
    fChain->SetBranchAddress("hit_channel_id", hit_channel_id, &b_hit_channel_id);
    fChain->SetBranchAddress("hit_index", hit_index, &b_hit_index);
    fChain->SetBranchAddress("hit_location", hit_location, &b_hit_location);
    fChain->SetBranchAddress("hit_is_mc", hit_is_mc, &b_hit_is_mc);
    fChain->SetBranchAddress("hit_num", hit_num, &b_hit_num);
    fChain->SetBranchAddress("hit_pixel", hit_pixel, &b_hit_pixel);
    fChain->SetBranchAddress("hit_board", hit_board, &b_hit_board);
    fChain->SetBranchAddress("hit_chain", hit_chain, &b_hit_chain);
    fChain->SetBranchAddress("hit_croc", hit_croc, &b_hit_croc);
    fChain->SetBranchAddress("hit_crate", hit_crate, &b_hit_crate);
    fChain->SetBranchAddress("hit_link", hit_link, &b_hit_link);
    fChain->SetBranchAddress("hit_disc_fired", hit_disc_fired, &b_hit_disc_fired);
    fChain->SetBranchAddress("hit_sys_ticks", hit_sys_ticks, &b_hit_sys_ticks);
    fChain->SetBranchAddress("hit_delay_ticks", hit_delay_ticks, &b_hit_delay_ticks);
    fChain->SetBranchAddress("hit_quarter_ticks", hit_quarter_ticks, &b_hit_quarter_ticks);
    fChain->SetBranchAddress("hit_qlo", hit_qlo, &b_hit_qlo);
    fChain->SetBranchAddress("hit_qmed", hit_qmed, &b_hit_qmed);
    fChain->SetBranchAddress("hit_qhi", hit_qhi, &b_hit_qhi);
    fChain->SetBranchAddress("n_idhits", &n_idhits, &b_n_idhits);
    fChain->SetBranchAddress("hits_id_per_mod", &hits_id_per_mod, &b_hits_id_per_mod);
    fChain->SetBranchAddress("hit_strip", hit_strip, &b_hit_strip);
    fChain->SetBranchAddress("hit_plane", hit_plane, &b_hit_plane);
    fChain->SetBranchAddress("hit_module", hit_module, &b_hit_module);
    fChain->SetBranchAddress("hit_view", hit_view, &b_hit_view);
    fChain->SetBranchAddress("n_odhits", &n_odhits, &b_n_odhits);
    fChain->SetBranchAddress("hits_od_per_mod", &hits_od_per_mod, &b_hits_od_per_mod);
    fChain->SetBranchAddress("hit_bar", hit_bar, &b_hit_bar);
    fChain->SetBranchAddress("hit_story", hit_story, &b_hit_story);
    fChain->SetBranchAddress("hit_tower", hit_tower, &b_hit_tower);
    fChain->SetBranchAddress("hit_frame", hit_frame, &b_hit_frame);
    fChain->SetBranchAddress("n_vetohits", &n_vetohits, &b_n_vetohits);
    fChain->SetBranchAddress("hit_wall", hit_wall, &b_hit_wall);
    fChain->SetBranchAddress("hit_paddle", hit_paddle, &b_hit_paddle);
    fChain->SetBranchAddress("hit_pmt", hit_pmt, &b_hit_pmt);
    fChain->SetBranchAddress("hit_q", hit_q, &b_hit_q);
    fChain->SetBranchAddress("hit_pe", hit_pe, &b_hit_pe);
    fChain->SetBranchAddress("hit_norm_energy", hit_norm_energy, &b_hit_norm_energy);
    fChain->SetBranchAddress("hit_time_raw", hit_time_raw, &b_hit_time_raw);
    fChain->SetBranchAddress("hit_time", hit_time, &b_hit_time);
    fChain->SetBranchAddress("hit_time_slice", hit_time_slice, &b_hit_time_slice);
    fChain->SetBranchAddress("hits_total_pe", hits_total_pe, &b_hits_total_pe);
    fChain->SetBranchAddress("hit_user_color", hit_user_color, &b_hit_user_color);
    fChain->SetBranchAddress("n_clusters_id", &n_clusters_id, &b_n_clusters_id);
    fChain->SetBranchAddress("clus_id_index", clus_id_index, &b_clus_id_index);
    fChain->SetBranchAddress("clus_id_strip", clus_id_strip, &b_clus_id_strip);
    fChain->SetBranchAddress("clus_id_plane", clus_id_plane, &b_clus_id_plane);
    fChain->SetBranchAddress("clus_id_module", clus_id_module, &b_clus_id_module);
    fChain->SetBranchAddress("clus_id_coord", clus_id_coord, &b_clus_id_coord);
    fChain->SetBranchAddress("clus_id_coordErr", clus_id_coordErr, &b_clus_id_coordErr);
    fChain->SetBranchAddress("clus_id_width", clus_id_width, &b_clus_id_width);
    fChain->SetBranchAddress("clus_id_tpos1", clus_id_tpos1, &b_clus_id_tpos1);
    fChain->SetBranchAddress("clus_id_tpos2", clus_id_tpos2, &b_clus_id_tpos2);
    fChain->SetBranchAddress("clus_id_lpos", clus_id_lpos, &b_clus_id_lpos);
    fChain->SetBranchAddress("clus_id_z", clus_id_z, &b_clus_id_z);
    fChain->SetBranchAddress("clus_id_view", clus_id_view, &b_clus_id_view);
    fChain->SetBranchAddress("clus_id_type", clus_id_type, &b_clus_id_type);
    fChain->SetBranchAddress("clus_id_hist", clus_id_hist, &b_clus_id_hist);
    fChain->SetBranchAddress("clus_id_subdet", clus_id_subdet, &b_clus_id_subdet);
    fChain->SetBranchAddress("clus_id_pe", clus_id_pe, &b_clus_id_pe);
    fChain->SetBranchAddress("clus_id_energy", clus_id_energy, &b_clus_id_energy);
    fChain->SetBranchAddress("clus_id_time", clus_id_time, &b_clus_id_time);
    fChain->SetBranchAddress("clus_id_time_slice", clus_id_time_slice, &b_clus_id_time_slice);
    fChain->SetBranchAddress("clus_id_size", clus_id_size, &b_clus_id_size);
    fChain->SetBranchAddress("clus_id_hits_idx", clus_id_hits_idx, &b_clus_id_hits_idx);
    fChain->SetBranchAddress("clus_id_usedFor", clus_id_usedFor, &b_clus_id_usedFor);
    fChain->SetBranchAddress("n_clusters_od", &n_clusters_od, &b_n_clusters_od);
    fChain->SetBranchAddress("clus_od_index", clus_od_index, &b_clus_od_index);
    fChain->SetBranchAddress("clus_od_z", clus_od_z, &b_clus_od_z);
    fChain->SetBranchAddress("clus_od_frame", clus_od_frame, &b_clus_od_frame);
    fChain->SetBranchAddress("clus_od_tower", clus_od_tower, &b_clus_od_tower);
    fChain->SetBranchAddress("clus_od_story", clus_od_story, &b_clus_od_story);
    fChain->SetBranchAddress("clus_od_hist", clus_od_hist, &b_clus_od_hist);
    fChain->SetBranchAddress("clus_od_pe", clus_od_pe, &b_clus_od_pe);
    fChain->SetBranchAddress("clus_od_energy", clus_od_energy, &b_clus_od_energy);
    fChain->SetBranchAddress("clus_od_time", clus_od_time, &b_clus_od_time);
    fChain->SetBranchAddress("clus_od_time_slice", clus_od_time_slice, &b_clus_od_time_slice);
    fChain->SetBranchAddress("clus_od_size", clus_od_size, &b_clus_od_size);
    fChain->SetBranchAddress("clus_od_hits_idx", clus_od_hits_idx, &b_clus_od_hits_idx);
    fChain->SetBranchAddress("n_blobs_id", &n_blobs_id, &b_n_blobs_id);
    fChain->SetBranchAddress("blob_id_idx", &blob_id_idx, &b_blob_id_idx);
    fChain->SetBranchAddress("blob_id_subdet", &blob_id_subdet, &b_blob_id_subdet);
    fChain->SetBranchAddress("blob_id_history", &blob_id_history, &b_blob_id_history);
    fChain->SetBranchAddress("blob_id_size", &blob_id_size, &b_blob_id_size);
    fChain->SetBranchAddress("blob_id_patrec", &blob_id_patrec, &b_blob_id_patrec);
    fChain->SetBranchAddress("blob_id_e", &blob_id_e, &b_blob_id_e);
    fChain->SetBranchAddress("blob_id_time", &blob_id_time, &b_blob_id_time);
    fChain->SetBranchAddress("blob_id_time_slice", &blob_id_time_slice, &b_blob_id_time_slice);
    fChain->SetBranchAddress("blob_id_startpoint_x", &blob_id_startpoint_x, &b_blob_id_startpoint_x);
    fChain->SetBranchAddress("blob_id_startpoint_y", &blob_id_startpoint_y, &b_blob_id_startpoint_y);
    fChain->SetBranchAddress("blob_id_startpoint_z", &blob_id_startpoint_z, &b_blob_id_startpoint_z);
    fChain->SetBranchAddress("blob_id_clus_idx", &blob_id_clus_idx, &b_blob_id_clus_idx);
    fChain->SetBranchAddress("n_blobs_od", &n_blobs_od, &b_n_blobs_od);
    fChain->SetBranchAddress("blob_od_idx", &blob_od_idx, &b_blob_od_idx);
    fChain->SetBranchAddress("blob_od_history", &blob_od_history, &b_blob_od_history);
    fChain->SetBranchAddress("blob_od_size", &blob_od_size, &b_blob_od_size);
    fChain->SetBranchAddress("blob_od_patrec", &blob_od_patrec, &b_blob_od_patrec);
    fChain->SetBranchAddress("blob_od_e", &blob_od_e, &b_blob_od_e);
    fChain->SetBranchAddress("blob_od_time", &blob_od_time, &b_blob_od_time);
    fChain->SetBranchAddress("blob_od_time_slice", &blob_od_time_slice, &b_blob_od_time_slice);
    fChain->SetBranchAddress("blob_od_clus_idx", &blob_od_clus_idx, &b_blob_od_clus_idx);
    fChain->SetBranchAddress("n_tracks", &n_tracks, &b_n_tracks);
    fChain->SetBranchAddress("trk_index", trk_index, &b_trk_index);
    fChain->SetBranchAddress("trk_type", trk_type, &b_trk_type);
    fChain->SetBranchAddress("trk_patrec", trk_patrec, &b_trk_patrec);
    fChain->SetBranchAddress("trk_time_slice", trk_time_slice, &b_trk_time_slice);
    fChain->SetBranchAddress("trk_vis_energy", trk_vis_energy, &b_trk_vis_energy);
    fChain->SetBranchAddress("trk_theta", trk_theta, &b_trk_theta);
    fChain->SetBranchAddress("trk_phi", trk_phi, &b_trk_phi);
    fChain->SetBranchAddress("trk_hits", trk_hits, &b_trk_hits);
    fChain->SetBranchAddress("trk_dof", trk_dof, &b_trk_dof);
    fChain->SetBranchAddress("trk_chi2perDof", trk_chi2perDof, &b_trk_chi2perDof);
    fChain->SetBranchAddress("trk_fitMass", trk_fitMass, &b_trk_fitMass);
    fChain->SetBranchAddress("trk_nodes", trk_nodes, &b_trk_nodes);
    fChain->SetBranchAddress("trk_node_X", trk_node_X, &b_trk_node_X);
    fChain->SetBranchAddress("trk_node_Y", trk_node_Y, &b_trk_node_Y);
    fChain->SetBranchAddress("trk_node_Z", trk_node_Z, &b_trk_node_Z);
    fChain->SetBranchAddress("trk_node_aX", trk_node_aX, &b_trk_node_aX);
    fChain->SetBranchAddress("trk_node_aY", trk_node_aY, &b_trk_node_aY);
    fChain->SetBranchAddress("trk_node_qOverP", trk_node_qOverP, &b_trk_node_qOverP);
    fChain->SetBranchAddress("trk_node_chi2", trk_node_chi2, &b_trk_node_chi2);
    fChain->SetBranchAddress("trk_node_cluster_idx", trk_node_cluster_idx, &b_trk_node_cluster_idx);
    fChain->SetBranchAddress("trk_usedFor", trk_usedFor, &b_trk_usedFor);
    fChain->SetBranchAddress("n_vertices", &n_vertices, &b_n_vertices);
    fChain->SetBranchAddress("vtx_time_slice", vtx_time_slice, &b_vtx_time_slice);
    fChain->SetBranchAddress("vtx_type", vtx_type, &b_vtx_type);
    fChain->SetBranchAddress("vtx_index", vtx_index, &b_vtx_index);
    fChain->SetBranchAddress("vtx_x", vtx_x, &b_vtx_x);
    fChain->SetBranchAddress("vtx_y", vtx_y, &b_vtx_y);
    fChain->SetBranchAddress("vtx_z", vtx_z, &b_vtx_z);
    fChain->SetBranchAddress("vtx_x_err", vtx_x_err, &b_vtx_x_err);
    fChain->SetBranchAddress("vtx_y_err", vtx_y_err, &b_vtx_y_err);
    fChain->SetBranchAddress("vtx_z_err", vtx_z_err, &b_vtx_z_err);
    fChain->SetBranchAddress("vtx_n_tracks", vtx_n_tracks, &b_vtx_n_tracks);
    fChain->SetBranchAddress("vtx_tracks_idx", vtx_tracks_idx, &b_vtx_tracks_idx);
    fChain->SetBranchAddress("n_minos_trk", &n_minos_trk, &b_n_minos_trk);
    fChain->SetBranchAddress("minos_run", &minos_run, &b_minos_run);
    fChain->SetBranchAddress("minos_subrun", &minos_subrun, &b_minos_subrun);
    fChain->SetBranchAddress("minos_snarl", &minos_snarl, &b_minos_snarl);
    fChain->SetBranchAddress("minos_trk_idx", &minos_trk_idx, &b_minos_trk_idx);
    fChain->SetBranchAddress("minos_trk_minervatrk_idx", &minos_trk_minervatrk_idx, &b_minos_trk_minervatrk_idx);
    fChain->SetBranchAddress("minos_trk_quality", &minos_trk_quality, &b_minos_trk_quality);
    fChain->SetBranchAddress("minos_trk_pass", &minos_trk_pass, &b_minos_trk_pass);
    fChain->SetBranchAddress("minos_trk_chi2", &minos_trk_chi2, &b_minos_trk_chi2);
    fChain->SetBranchAddress("minos_trk_ndf", &minos_trk_ndf, &b_minos_trk_ndf);
    fChain->SetBranchAddress("minos_trk_bave", &minos_trk_bave, &b_minos_trk_bave);
    fChain->SetBranchAddress("minos_trk_range", &minos_trk_range, &b_minos_trk_range);
    fChain->SetBranchAddress("minos_trk_con", &minos_trk_con, &b_minos_trk_con);
    fChain->SetBranchAddress("minos_trk_p", &minos_trk_p, &b_minos_trk_p);
    fChain->SetBranchAddress("minos_trk_prange", &minos_trk_prange, &b_minos_trk_prange);
    fChain->SetBranchAddress("minos_trk_qp", &minos_trk_qp, &b_minos_trk_qp);
    fChain->SetBranchAddress("minos_trk_eqp", &minos_trk_eqp, &b_minos_trk_eqp);
    fChain->SetBranchAddress("minos_trk_vtxp", &minos_trk_vtxp, &b_minos_trk_vtxp);
    fChain->SetBranchAddress("minos_trk_vtxu", &minos_trk_vtxu, &b_minos_trk_vtxu);
    fChain->SetBranchAddress("minos_trk_vtxv", &minos_trk_vtxv, &b_minos_trk_vtxv);
    fChain->SetBranchAddress("minos_trk_vtxx", &minos_trk_vtxx, &b_minos_trk_vtxx);
    fChain->SetBranchAddress("minos_trk_vtxy", &minos_trk_vtxy, &b_minos_trk_vtxy);
    fChain->SetBranchAddress("minos_trk_vtxz", &minos_trk_vtxz, &b_minos_trk_vtxz);
    fChain->SetBranchAddress("minos_trk_vtxt", &minos_trk_vtxt, &b_minos_trk_vtxt);
    fChain->SetBranchAddress("minos_trk_mvax", &minos_trk_mvax, &b_minos_trk_mvax);
    fChain->SetBranchAddress("minos_trk_mvau", &minos_trk_mvau, &b_minos_trk_mvau);
    fChain->SetBranchAddress("minos_trk_mvav", &minos_trk_mvav, &b_minos_trk_mvav);
    fChain->SetBranchAddress("minos_trk_vtx_dxdz", &minos_trk_vtx_dxdz, &b_minos_trk_vtx_dxdz);
    fChain->SetBranchAddress("minos_trk_vtx_dydz", &minos_trk_vtx_dydz, &b_minos_trk_vtx_dydz);
    fChain->SetBranchAddress("minos_trk_vtx_dudz", &minos_trk_vtx_dudz, &b_minos_trk_vtx_dudz);
    fChain->SetBranchAddress("minos_trk_vtx_dvdz", &minos_trk_vtx_dvdz, &b_minos_trk_vtx_dvdz);
    fChain->SetBranchAddress("minos_trk_endp", &minos_trk_endp, &b_minos_trk_endp);
    fChain->SetBranchAddress("minos_trk_endu", &minos_trk_endu, &b_minos_trk_endu);
    fChain->SetBranchAddress("minos_trk_endv", &minos_trk_endv, &b_minos_trk_endv);
    fChain->SetBranchAddress("minos_trk_endx", &minos_trk_endx, &b_minos_trk_endx);
    fChain->SetBranchAddress("minos_trk_endy", &minos_trk_endy, &b_minos_trk_endy);
    fChain->SetBranchAddress("minos_trk_endz", &minos_trk_endz, &b_minos_trk_endz);
    fChain->SetBranchAddress("minos_trk_endt", &minos_trk_endt, &b_minos_trk_endt);
    fChain->SetBranchAddress("minos_trk_ns", &minos_trk_ns, &b_minos_trk_ns);
    fChain->SetBranchAddress("minos_trk_stp_fit", &minos_trk_stp_fit, &b_minos_trk_stp_fit);
    fChain->SetBranchAddress("minos_trk_stp_u", &minos_trk_stp_u, &b_minos_trk_stp_u);
    fChain->SetBranchAddress("minos_trk_stp_v", &minos_trk_stp_v, &b_minos_trk_stp_v);
    fChain->SetBranchAddress("minos_trk_stp_x", &minos_trk_stp_x, &b_minos_trk_stp_x);
    fChain->SetBranchAddress("minos_trk_stp_y", &minos_trk_stp_y, &b_minos_trk_stp_y);
    fChain->SetBranchAddress("minos_trk_stp_z", &minos_trk_stp_z, &b_minos_trk_stp_z);
    fChain->SetBranchAddress("minos_trk_stp_t", &minos_trk_stp_t, &b_minos_trk_stp_t);
    fChain->SetBranchAddress("minos_trk_stp_meu", &minos_trk_stp_meu, &b_minos_trk_stp_meu);
    fChain->SetBranchAddress("n_minos_stp", &n_minos_stp, &b_n_minos_stp);
    fChain->SetBranchAddress("minos_stp_plane", &minos_stp_plane, &b_minos_stp_plane);
    fChain->SetBranchAddress("minos_stp_strip", &minos_stp_strip, &b_minos_stp_strip);
    fChain->SetBranchAddress("minos_stp_view", &minos_stp_view, &b_minos_stp_view);
    fChain->SetBranchAddress("minos_stp_tpos", &minos_stp_tpos, &b_minos_stp_tpos);
    fChain->SetBranchAddress("minos_stp_time", &minos_stp_time, &b_minos_stp_time);
    fChain->SetBranchAddress("minos_stp_z", &minos_stp_z, &b_minos_stp_z);
    fChain->SetBranchAddress("minos_stp_ph", &minos_stp_ph, &b_minos_stp_ph);
    fChain->SetBranchAddress("minos_stp_pe", &minos_stp_pe, &b_minos_stp_pe);
    fChain->SetBranchAddress("minos_stp_trkidx", &minos_stp_trkidx, &b_minos_stp_trkidx);
    fChain->SetBranchAddress("minos_sec", &minos_sec, &b_minos_sec);
    fChain->SetBranchAddress("minos_nanosec", &minos_nanosec, &b_minos_nanosec);
    fChain->SetBranchAddress("beam_pot", &beam_pot, &b_beam_pot);
    fChain->SetBranchAddress("beam_horncur", &beam_horncur, &b_beam_horncur);
    fChain->SetBranchAddress("beam_xpos", &beam_xpos, &b_beamxpos);
    fChain->SetBranchAddress("beam_ypos", &beam_ypos, &b_beam_ypos);
    fChain->SetBranchAddress("beam_xwid", &beam_xwid, &b_beam_xwid);
    fChain->SetBranchAddress("beam_ywid", &beam_ywid, &b_beam_ywid);
    fChain->SetBranchAddress("beam_dt_nearest", &beam_dt_nearest, &b_beam_dt_nearest);
    fChain->SetBranchAddress("beam_dt_ok", &beam_dt_ok, &b_beam_dt_ok);
    fChain->SetBranchAddress("beam_pos_ok", &beam_pos_ok, &b_beam_pos_ok);
    fChain->SetBranchAddress("beam_wid_ok", &beam_wid_ok, &b_beam_wid_ok);
    fChain->SetBranchAddress("beam_tor_ok", &beam_tor_ok, &b_beam_tor_ok);
    fChain->SetBranchAddress("beam_horns_ok", &beam_horns_ok, &b_beam_horns_ok);
    fChain->SetBranchAddress("beam_numibeamdb_sec", &beam_numibeamdb_sec, &b_beam_numibeamdb_sec);
    fChain->SetBranchAddress("beam_numibeamdb_nanosec", &beam_numibeamdb_nanosec, &b_beam_numibeamdb_nanosec);
    fChain->SetBranchAddress("beam_is_good_beam_spill", &beam_is_good_beam_spill, &b_beam_is_good_beam_spill);
    fChain->SetBranchAddress("beam_is_bad_pot_data_spill", &beam_is_bad_pot_data_spill, &b_beam_is_bad_pot_data_spill);
    fChain->SetBranchAddress("beam_is_no_beam_spill", &beam_is_no_beam_spill, &b_beam_is_no_beam_spill);
    fChain->SetBranchAddress("beam_is_bad_data_spill", &beam_is_bad_data_spill, &b_beam_is_bad_data_spill);
    fChain->SetBranchAddress("beam_is_bad_prof_widx_data_spill", &beam_is_bad_prof_widx_data_spill, &b_beam_is_bad_prof_widx_data_spill);
    fChain->SetBranchAddress("beam_is_bad_prof_widy_data_spill", &beam_is_bad_prof_widy_data_spill, &b_beam_is_bad_prof_widy_data_spill);
    fChain->SetBranchAddress("beam_is_bad_xpos_data_spill", &beam_is_bad_xpos_data_spill, &b_beam_is_bad_xpos_data_spill);
    fChain->SetBranchAddress("beam_is_bad_ypos_data_spill", &beam_is_bad_ypos_data_spill, &b_beam_is_bad_ypos_data_spill);
    fChain->SetBranchAddress("beam_is_bad_horncur_data_spill", &beam_is_bad_horncur_data_spill, &b_beam_is_bad_horncur_data_spill);
    fChain->SetBranchAddress("beam_is_bad_nearesttime_spill", &beam_is_bad_nearesttime_spill, &b_beam_is_bad_nearesttime_spill);
    fChain->SetBranchAddress("beam_is_bad_beam_spill", &beam_is_bad_beam_spill, &b_beam_is_bad_beam_spill);
    fChain->SetBranchAddress("beam_is_bad_pot_spill", &beam_is_bad_pot_spill, &b_beam_is_bad_pot_spill);
    fChain->SetBranchAddress("beam_is_bad_xpos_spill", &beam_is_bad_xpos_spill, &b_beam_is_bad_xpos_spill);
    fChain->SetBranchAddress("beam_is_bad_ypos_spill", &beam_is_bad_ypos_spill, &b_beam_is_bad_ypos_spill);
    fChain->SetBranchAddress("beam_is_bad_beamsize_spill", &beam_is_bad_beamsize_spill, &b_beam_is_bad_beamsize_spill);
    fChain->SetBranchAddress("beam_is_bad_prof_widx_spill", &beam_is_bad_prof_widx_spill, &b_beam_is_bad_prof_widx_spill);
    fChain->SetBranchAddress("beam_is_bad_prof_widy_spill", &beam_is_bad_prof_widy_spill, &b_beam_is_bad_prof_widy_spill);
    fChain->SetBranchAddress("beam_is_bad_horncur_spill", &beam_is_bad_horncur_spill, &b_beam_is_bad_horncur_spill);
    fChain->SetBranchAddress("beam_is_target_out_spill", &beam_is_target_out_spill, &b_beam_is_target_out_spill);
    fChain->SetBranchAddress("beam_is_bad_beamtype_spill", &beam_is_bad_beamtype_spill, &b_beam_is_bad_beamtype_spill);
    fChain->SetBranchAddress("beam_is_bad_beam_frac_on_tgt_spill", &beam_is_bad_beam_frac_on_tgt_spill, &b_beam_is_bad_beam_frac_on_tgt_spill);
    fChain->SetBranchAddress("beam_tor101", &beam_tor101, &b_beam_tor101);
    fChain->SetBranchAddress("beam_tr101d", &beam_tr101d, &b_beam_tr101d);
    fChain->SetBranchAddress("beam_tortgt", &beam_tortgt, &b_beam_tortgt);
    fChain->SetBranchAddress("beam_trtgtd", &beam_trtgtd, &b_beam_trtgtd);
    fChain->SetBranchAddress("mc_run", &mc_run, &b_mc_run);
    fChain->SetBranchAddress("mc_subrun", &mc_subrun, &b_mc_subrun);
    fChain->SetBranchAddress("mc_spill", &mc_spill, &b_mc_spill);
    fChain->SetBranchAddress("n_total_interactions", &n_total_interactions, &b_n_total_interactions);
    fChain->SetBranchAddress("mc_MIState", &mc_MIState, &b_mc_MIState);
    fChain->SetBranchAddress("mc_pot", &mc_pot, &b_mc_pot);
    fChain->SetBranchAddress("mc_beamConfig", &mc_beamConfig, &b_mc_beamConfig);
    fChain->SetBranchAddress("n_interactions", &n_interactions, &b_n_interactions);
    fChain->SetBranchAddress("mc_int_processType", mc_int_processType, &b_mc_int_processType);
    fChain->SetBranchAddress("mc_int_nevSpill", mc_int_nevSpill, &b_mc_int_nevSpill);
    fChain->SetBranchAddress("mc_int_nevFile", mc_int_nevFile, &b_mc_int_nevFile);
    fChain->SetBranchAddress("mc_int_channel", mc_int_channel, &b_mc_int_channel);
    fChain->SetBranchAddress("mc_int_current", mc_int_current, &b_mc_int_current);
    fChain->SetBranchAddress("mc_int_charm", mc_int_charm, &b_mc_int_charm);
    fChain->SetBranchAddress("mc_int_weight", mc_int_weight, &b_mc_int_weight);
    fChain->SetBranchAddress("mc_int_xSection", mc_int_xSection, &b_mc_int_xSection);
    fChain->SetBranchAddress("mc_int_incomingPDG", mc_int_incomingPDG, &b_mc_int_incomingPDG);
    fChain->SetBranchAddress("mc_int_tgtNucleus", mc_int_tgtNucleus, &b_mc_int_tgtNucleus);
    fChain->SetBranchAddress("mc_int_tgtNucleon", mc_int_tgtNucleon, &b_mc_int_tgtNucleon);
    fChain->SetBranchAddress("mc_int_targetZ", mc_int_targetZ, &b_mc_int_targetZ);
    fChain->SetBranchAddress("mc_int_targetA", mc_int_targetA, &b_mc_int_targetA);
    fChain->SetBranchAddress("mc_int_hitQuark", mc_int_hitQuark, &b_mc_int_hitQuark);
    fChain->SetBranchAddress("mc_int_seaQuark", mc_int_seaQuark, &b_mc_int_seaQuark);
    fChain->SetBranchAddress("mc_int_resID", mc_int_resID, &b_mc_int_resID);
    fChain->SetBranchAddress("mc_int_FSLepton", mc_int_FSLepton, &b_mc_int_FSLepton);
    fChain->SetBranchAddress("mc_int_incomingE", mc_int_incomingE, &b_mc_int_incomingE);
    fChain->SetBranchAddress("mc_int_bjorkenX", mc_int_bjorkenX, &b_mc_int_bjorkenX);
    fChain->SetBranchAddress("mc_int_bjorkenY", mc_int_bjorkenY, &b_mc_int_bjorkenY);
    fChain->SetBranchAddress("mc_int_QSquared", mc_int_QSquared, &b_mc_int_QSquared);
    fChain->SetBranchAddress("mc_int_nucleonT", mc_int_nucleonT, &b_mc_int_nucleonT);
    fChain->SetBranchAddress("mc_int_W", mc_int_W, &b_mc_int_W);
    fChain->SetBranchAddress("mc_int_nFSParticles", mc_int_nFSParticles, &b_mc_int_nFSParticles);
    fChain->SetBranchAddress("mc_int_vtx", mc_int_vtx, &b_mc_int_vtx);
    fChain->SetBranchAddress("mc_int_incoming4p", mc_int_incoming4p, &b_mc_int_incoming4p);
    fChain->SetBranchAddress("mc_int_tgtNucleon4p", mc_int_tgtNucleon4p, &b_mc_int_tgtNucleon4p);
    fChain->SetBranchAddress("mc_int_FSLepton4p", mc_int_FSLepton4p, &b_mc_int_FSLepton4p);
    fChain->SetBranchAddress("mc_int_FSPdg", mc_int_FSPdg, &b_mc_int_FSPdg);
    fChain->SetBranchAddress("mc_int_FSParticlesPx", mc_int_FSParticlesPx, &b_mc_int_FSParticlesPx);
    fChain->SetBranchAddress("mc_int_FSParticlesPy", mc_int_FSParticlesPy, &b_mc_int_FSParticlesPy);
    fChain->SetBranchAddress("mc_int_FSParticlesPz", mc_int_FSParticlesPz, &b_mc_int_FSParticlesPz);
    fChain->SetBranchAddress("mc_int_FSParticlesE", mc_int_FSParticlesE, &b_mc_int_FSParticlesE);
    fChain->SetBranchAddress("mc_flux_proton_P", mc_flux_proton_P, &b_mc_flux_proton_P);
    fChain->SetBranchAddress("mc_flux_proton_X", mc_flux_proton_X, &b_mc_flux_proton_X);
    fChain->SetBranchAddress("mc_flux_parent_PDG", mc_flux_parent_PDG, &b_mc_flux_parent_PDG);
    fChain->SetBranchAddress("mc_flux_parent_prod4P", mc_flux_parent_prod4P, &b_mc_flux_parent_prod4P);
    fChain->SetBranchAddress("mc_flux_parent_prodPos", mc_flux_parent_prodPos, &b_mc_flux_parent_prodPos);
    fChain->SetBranchAddress("mc_flux_parent_decay4P", mc_flux_parent_decay4P, &b_mc_flux_parent_decay4P);
    fChain->SetBranchAddress("mc_flux_parent_decayPos", mc_flux_parent_decayPos, &b_mc_flux_parent_decayPos);
    fChain->SetBranchAddress("mc_flux_parent_generation", mc_flux_parent_generation, &b_mc_flux_parent_generation);
    fChain->SetBranchAddress("mc_flux_parent_decayMode", mc_flux_parent_decayMode, &b_mc_flux_parent_decayMode);
    fChain->SetBranchAddress("mc_flux_secondary_PDG", mc_flux_secondary_PDG, &b_mc_flux_secondary_PDG);
    fChain->SetBranchAddress("n_mc_trajectories", &n_mc_trajectories, &b_n_mc_trajectories);
    fChain->SetBranchAddress("mc_traj_overflow", &mc_traj_overflow, &b_mc_traj_overflow);
    fChain->SetBranchAddress("mc_traj_strlength", mc_traj_strlength, &b_mc_traj_strlength);
    fChain->SetBranchAddress("mc_traj_curvlength", mc_traj_curvlength, &b_mc_traj_curvlength);
    fChain->SetBranchAddress("mc_traj_leaving", mc_traj_leaving, &b_mc_traj_leaving);
    fChain->SetBranchAddress("mc_traj_trkid", mc_traj_trkid, &b_mc_traj_trkid);
    fChain->SetBranchAddress("mc_traj_parentid", mc_traj_parentid, &b_mc_traj_parentid);
    fChain->SetBranchAddress("mc_traj_pdg", mc_traj_pdg, &b_mc_traj_pdg);
    fChain->SetBranchAddress("mc_traj_hit_e", mc_traj_hit_e, &b_mc_traj_hit_e);
    fChain->SetBranchAddress("mc_traj_npoints", mc_traj_npoints, &b_mc_traj_npoints);
    fChain->SetBranchAddress("mc_traj_point_x", mc_traj_point_x, &b_mc_traj_point_x);
    fChain->SetBranchAddress("mc_traj_point_y", mc_traj_point_y, &b_mc_traj_point_y);
    fChain->SetBranchAddress("mc_traj_point_z", mc_traj_point_z, &b_mc_traj_point_z);
    fChain->SetBranchAddress("mc_traj_point_t", mc_traj_point_t, &b_mc_traj_point_t);
    fChain->SetBranchAddress("mc_traj_point_px", mc_traj_point_px, &b_mc_traj_point_px);
    fChain->SetBranchAddress("mc_traj_point_py", mc_traj_point_py, &b_mc_traj_point_py);
    fChain->SetBranchAddress("mc_traj_point_pz", mc_traj_point_pz, &b_mc_traj_point_pz);
    fChain->SetBranchAddress("mc_traj_point_E", mc_traj_point_E, &b_mc_traj_point_E);
    fChain->SetBranchAddress("n_mc_id_digits", &n_mc_id_digits, &b_n_mc_id_digits);
    fChain->SetBranchAddress("mc_id_strip", mc_id_strip, &b_mc_id_strip);
    fChain->SetBranchAddress("mc_id_plane", mc_id_plane, &b_mc_id_plane);
    fChain->SetBranchAddress("mc_id_module", mc_id_module, &b_mc_id_module);
    fChain->SetBranchAddress("mc_id_view", mc_id_view, &b_mc_id_view);
    fChain->SetBranchAddress("mc_id_pe", mc_id_pe, &b_mc_id_pe);
    fChain->SetBranchAddress("mc_id_time", mc_id_time, &b_mc_id_time);
    fChain->SetBranchAddress("mc_id_dE", mc_id_dE, &b_mc_id_dE);
    fChain->SetBranchAddress("mc_id_nmchit", mc_id_nmchit, &b_mc_id_nmchit);
    fChain->SetBranchAddress("mc_id_mchit_x", mc_id_mchit_x, &b_mc_id_mchit_x);
    fChain->SetBranchAddress("mc_id_mchit_y", mc_id_mchit_y, &b_mc_id_mchit_y);
    fChain->SetBranchAddress("mc_id_mchit_z", mc_id_mchit_z, &b_mc_id_mchit_z);
    fChain->SetBranchAddress("mc_id_mchit_t", mc_id_mchit_t, &b_mc_id_mchit_t);
    fChain->SetBranchAddress("mc_id_mchit_trkid", mc_id_mchit_trkid, &b_mc_id_mchit_trkid);
    fChain->SetBranchAddress("mc_id_mchit_p", mc_id_mchit_p, &b_mc_id_mchit_p);
    fChain->SetBranchAddress("mc_id_mchit_dE", mc_id_mchit_dE, &b_mc_id_mchit_dE);
    fChain->SetBranchAddress("mc_id_mchit_dL", mc_id_mchit_dL, &b_mc_id_mchit_dL);
    fChain->SetBranchAddress("n_mc_od_digits", &n_mc_od_digits, &b_n_mc_od_digits);
    fChain->SetBranchAddress("mc_od_frame", mc_od_frame, &b_mc_od_frame);
    fChain->SetBranchAddress("mc_od_tower", mc_od_tower, &b_mc_od_tower);
    fChain->SetBranchAddress("mc_od_story", mc_od_story, &b_mc_od_story);
    fChain->SetBranchAddress("mc_od_bar", mc_od_bar, &b_mc_od_bar);
    fChain->SetBranchAddress("mc_od_pe", mc_od_pe, &b_mc_od_pe);
    fChain->SetBranchAddress("mc_od_time", mc_od_time, &b_mc_od_time);
    fChain->SetBranchAddress("mc_od_dE", mc_od_dE, &b_mc_od_dE);
    fChain->SetBranchAddress("mc_od_nmchit", mc_od_nmchit, &b_mc_od_nmchit);
    fChain->SetBranchAddress("mc_od_mchit_x", mc_od_mchit_x, &b_mc_od_mchit_x);
    fChain->SetBranchAddress("mc_od_mchit_y", mc_od_mchit_y, &b_mc_od_mchit_y);
    fChain->SetBranchAddress("mc_od_mchit_z", mc_od_mchit_z, &b_mc_od_mchit_z);
    fChain->SetBranchAddress("mc_od_mchit_t", mc_od_mchit_t, &b_mc_od_mchit_t);
    fChain->SetBranchAddress("mc_od_mchit_trkid", mc_od_mchit_trkid, &b_mc_od_mchit_trkid);
    fChain->SetBranchAddress("mc_od_mchit_p", mc_od_mchit_p, &b_mc_od_mchit_p);
    fChain->SetBranchAddress("mc_od_mchit_dE", mc_od_mchit_dE, &b_mc_od_mchit_dE);
    fChain->SetBranchAddress("n_mc_veto_digits", &n_mc_veto_digits, &b_n_mc_veto_digits);
    fChain->SetBranchAddress("mc_veto_wall", mc_veto_wall, &b_mc_veto_wall);
    fChain->SetBranchAddress("mc_veto_paddle", mc_veto_paddle, &b_mc_veto_paddle);
    fChain->SetBranchAddress("mc_veto_pe", mc_veto_pe, &b_mc_veto_pe);
    fChain->SetBranchAddress("mc_veto_time", mc_veto_time, &b_mc_veto_time);
    fChain->SetBranchAddress("mc_veto_dE", mc_veto_dE, &b_mc_veto_dE);
    fChain->SetBranchAddress("mc_veto_nmchit", mc_veto_nmchit, &b_mc_veto_nmchit);
    fChain->SetBranchAddress("mc_veto_mchit_x", mc_veto_mchit_x, &b_mc_veto_mchit_x);
    fChain->SetBranchAddress("mc_veto_mchit_y", mc_veto_mchit_y, &b_mc_veto_mchit_y);
    fChain->SetBranchAddress("mc_veto_mchit_z", mc_veto_mchit_z, &b_mc_veto_mchit_z);
    fChain->SetBranchAddress("mc_veto_mchit_t", mc_veto_mchit_t, &b_mc_veto_mchit_t);
    fChain->SetBranchAddress("mc_veto_mchit_trkid", mc_veto_mchit_trkid, &b_mc_veto_mchit_trkid);
    fChain->SetBranchAddress("mc_veto_mchit_p", mc_veto_mchit_p, &b_mc_veto_mchit_p);
    fChain->SetBranchAddress("mc_veto_mchit_dE", mc_veto_mchit_dE, &b_mc_veto_mchit_dE);
}

#endif
