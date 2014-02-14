/*
================================================================================
Class: CCInclusive
    Core Class for CCNuPionInc Analysis Package
    Includes all data variables for Data Analysis
    Member functions uses member data variables
    
    Outputs .ROOT file filled with 1D or 2D histograms, which is an input file
    for Plotter Class
    
    Uses other classes to define variables or to access specific functions
    
    Main Directory:
        Classes/CCInclusive
        
    Usage:
        > main.cpp declares and controls the class
        > See run function Comments
    
    
    
    Last Revision: 2014_02_14
================================================================================
*/

#ifndef CCInclusive_h
#define CCInclusive_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>

#include <vector>

#include "Libraries/PDG_List.h"
#include "Libraries/Folder_List.h"
#include "Classes/HEP_Misc/HEP_Misc.cpp"
#include "Classes/Particle/Particle.cpp"
#include "Classes/CutNumberList/CutNumberList.cpp"
#include "Classes/BinList/BinList.cpp"


class CCInclusive {
public :
   // -------------------------------------------------------------------------
   //     Specific Functions
   //--------------------------------------------------------------------------
   
    
    
   // -------------------------------------------------------------------------
   //     void run(): Generates a .root file with selected histograms
   //         playlist -> address of the playlist
   //         filename -> file name for the output .root file
   //--------------------------------------------------------------------------
    void run(string playlist, string rootFilename);
    
    
   /* 
   -----------------------------------------------------------------------------
        void initVariables();
        void initHistograms();
         
            Initialize all variables will be used in your analysis
            1) Initiliaze the CutNumbers you want to use
            2) Initiliaze the Histograms that will be created
            
                
   -----------------------------------------------------------------------------
   */
    void initVariables();
    void initHistograms();
    
    void fillParticle(Particle* part);
    
    void fillHistograms();

    void openFiles();
    void closeFiles();
    void writeReadme();
    void writeCutFile();

   // -------------------------------------------------------------------------
   //     Default Functions
   //--------------------------------------------------------------------------
    CCInclusive();
    ~CCInclusive();

    void Init(string playlist, TChain* fChain);
    Int_t GetEntry(Long64_t entry);
    Long64_t LoadTree(Long64_t entry);
    Bool_t Notify();
    void Show(Long64_t entry);
    Int_t Cut(Long64_t entry);
    
   // -------------------------------------------------------------------------
   //     Histograms
   //--------------------------------------------------------------------------
    // True Analysis Variables
    TH1F* beamEnergy;
   
    TH1F* P_muon;
    TH1F* P_proton;
    TH1F* P_pion;
    TH1F* Angle_muon;
    TH1F* Angle_proton;
    TH1F* Angle_pion;
    TH1F* AngleMuon_muon;
    TH1F* AngleMuon_proton;
    TH1F* AngleMuon_pion;
    

   // -------------------------------------------------------------------------
   //     Analysis Variables
   //--------------------------------------------------------------------------
    double maxBeamEnergy;
    HEP_Misc misc;
    TVector3* beam_p3;
    Particle* proton;
    Particle* muon;
    Particle* pion;
    
   // -------------------------------------------------------------------------
   //     Cut Numbers and List of Bins
   //--------------------------------------------------------------------------
    CutNumberList* nCutList;
    BinList* binList;
    
   // -------------------------------------------------------------------------
   //     Files
   //--------------------------------------------------------------------------
    string readmeFile;
    ofstream readme;

   // -------------------------------------------------------------------------
   //     Data
   //--------------------------------------------------------------------------

  // Declaration of leaf types
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
   Int_t           broken_track_most_us_plane;
   Int_t           n_associated_recoil_subprongs;
   Int_t           n_muon_blobs;
   Int_t           n_recoil_iso2Dblobs;
   Int_t           n_recoil_iso3Dblobs;
   Int_t           n_recoil_subprongs;
   Int_t           n_recoil_vtxblobs;
   Int_t           n_tracks;
   Int_t           n_tracks_non_prim;
   Int_t           n_tracks_prim;
   Int_t           n_tracks_prim_forked;
   Int_t           n_tracks_prim_kinked;
   Int_t           n_vertices_startpoint;
   Int_t           pass_canonical_cut;
   Int_t           pass_minosmatch_cut;
   Int_t           pass_nearccinclusive_cut;
   Int_t           phys_energy_in_road_downstream_nplanes;
   Int_t           phys_energy_in_road_upstream_nplanes;
   Int_t           phys_n_dead_discr_pair;
   Int_t           phys_n_dead_discr_pair_in_prim_track_region;
   Int_t           phys_n_dead_discr_pair_two_mod_downstream_prim_track;
   Int_t           phys_n_dead_discr_pair_two_mod_upstream_prim_vtx;
   Int_t           phys_n_dead_discr_pair_upstream_prim_track_proj;
   Int_t           phys_prim_vertex_planeid;
   Int_t           phys_vertex_is_fiducial;
   Double_t        energy_from_mc;
   Double_t        energy_from_mc_fraction;
   Double_t        energy_from_mc_fraction_of_highest;
   Double_t        muon_blob_energy;
   Double_t        muon_phi;
   Double_t        muon_theta;
   Double_t        muon_thetaX;
   Double_t        muon_thetaY;
   Double_t        nu_topological_energy_recoil;
   Double_t        numi_horn_curr;
   Double_t        numi_pot;
   Double_t        numi_x;
   Double_t        numi_x_width;
   Double_t        numi_y;
   Double_t        numi_y_width;
   Double_t        phys_energy_dispersed;
   Double_t        phys_energy_in_road_downstream;
   Double_t        phys_energy_in_road_upstream;
   Double_t        phys_energy_unattached;
   Double_t        prim_vtx_smallest_opening_angle;
   Double_t        primary_track_minerva_energy;
   Double_t        primary_track_minerva_phi;
   Double_t        primary_track_minerva_theta;
   Double_t        vtxprong_energy_cal;
   Double_t        vtxprong_energy_visible;
   Double_t        vtxprong_energy_visible_isoblobs;
   Double_t        vtxprong_energy_visible_vtxblobs;
   Int_t           plane_id_sz;
   Int_t           plane_id[204];   //[plane_id_sz]
   Int_t           n_vtxprong_isoblobs;
   Int_t           vtxprong_isoblob_nclusters[1];   //[n_vtxprong_isoblobs]
   Int_t           n_vtxprong_vtxblobs;
   Int_t           vtxprong_vtxblob_nclusters[2];   //[n_vtxprong_vtxblobs]
   Int_t           clusterU_Angle_sz;
   Double_t        clusterU_Angle[107];   //[clusterU_Angle_sz]
   Int_t           clusterU_Radius_sz;
   Double_t        clusterU_Radius[107];   //[clusterU_Radius_sz]
   Int_t           clusterU_timeDiff_sz;
   Double_t        clusterU_timeDiff[107];   //[clusterU_timeDiff_sz]
   Int_t           clusterU_viewDist_sz;
   Double_t        clusterU_viewDist[107];   //[clusterU_viewDist_sz]
   Int_t           clusterU_visE_binned_sz;
   Double_t        clusterU_visE_binned[10];   //[clusterU_visE_binned_sz]
   Int_t           clusterU_visEnergy_sz;
   Double_t        clusterU_visEnergy[107];   //[clusterU_visEnergy_sz]
   Int_t           clusterU_zDist_sz;
   Double_t        clusterU_zDist[107];   //[clusterU_zDist_sz]
   Int_t           clusterV_Angle_sz;
   Double_t        clusterV_Angle[89];   //[clusterV_Angle_sz]
   Int_t           clusterV_Radius_sz;
   Double_t        clusterV_Radius[89];   //[clusterV_Radius_sz]
   Int_t           clusterV_timeDiff_sz;
   Double_t        clusterV_timeDiff[89];   //[clusterV_timeDiff_sz]
   Int_t           clusterV_viewDist_sz;
   Double_t        clusterV_viewDist[89];   //[clusterV_viewDist_sz]
   Int_t           clusterV_visE_binned_sz;
   Double_t        clusterV_visE_binned[10];   //[clusterV_visE_binned_sz]
   Int_t           clusterV_visEnergy_sz;
   Double_t        clusterV_visEnergy[89];   //[clusterV_visEnergy_sz]
   Int_t           clusterV_zDist_sz;
   Double_t        clusterV_zDist[89];   //[clusterV_zDist_sz]
   Int_t           clusterX_Angle_sz;
   Double_t        clusterX_Angle[179];   //[clusterX_Angle_sz]
   Int_t           clusterX_Radius_sz;
   Double_t        clusterX_Radius[179];   //[clusterX_Radius_sz]
   Int_t           clusterX_timeDiff_sz;
   Double_t        clusterX_timeDiff[179];   //[clusterX_timeDiff_sz]
   Int_t           clusterX_viewDist_sz;
   Double_t        clusterX_viewDist[179];   //[clusterX_viewDist_sz]
   Int_t           clusterX_visE_binned_sz;
   Double_t        clusterX_visE_binned[10];   //[clusterX_visE_binned_sz]
   Int_t           clusterX_visEnergy_sz;
   Double_t        clusterX_visEnergy[179];   //[clusterX_visEnergy_sz]
   Int_t           clusterX_zDist_sz;
   Double_t        clusterX_zDist[179];   //[clusterX_zDist_sz]
   Int_t           plane_visible_energy_sz;
   Double_t        plane_visible_energy[204];   //[plane_visible_energy_sz]
   Double_t        primary_track_minerva_end_position[3];
   Double_t        primary_track_minerva_start_position[3];
   Double_t        vtxprong_isoblob_visenergy[1];   //[n_vtxprong_isoblobs]
   Double_t        vtxprong_vtxblob_visenergy[2];   //[n_vtxprong_vtxblobs]
   Bool_t          truth_has_physics_event;
   Bool_t          truth_pass_CCInclusiveReco;
   Bool_t          truth_pass_plausible;
   Bool_t          truth_pass_fiducial;
   Bool_t          truth_pass_analyzable;
   Bool_t          truth_pass_fiducial_apothem;
   Bool_t          truth_pass_analyzable_apothem;
   Bool_t          truth_pass_ClassifyVertex;
   Bool_t          truth_pass_ClassifyCC;
   Bool_t          truth_pass_MuonEReconstructed;
   Bool_t          truth_pass_MuMinusSignSelection;
   Bool_t          truth_pass_MuPlusSignSelection;
   Bool_t          truth_pass_ClassifyTracks;
   Bool_t          truth_pass_ClassifyRecoil;
   Bool_t          truth_pass_Canonical;
   Bool_t          truth_reco_vertex_fiducial;
   Bool_t          truth_smeared_reco_vertex_fiducial;
   Bool_t          truth_pass_nearInclusive;
   Double_t        truth_muon_fraction_accepted;
   Int_t           genie_wgt_n_shifts;
   Double_t        truth_genie_wgt_AGKYxF1pi[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_AhtBY[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_BhtBY[7];   //[genie_wgt_n_shifts]
   Double_t        truth_genie_wgt_CCQEPauliSupViaFK[7];   //[genie_wgt_n_shifts]
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
   Int_t           CCInclusiveReco_nuFlavor;
   Int_t           CCInclusiveReco_nuHelicity;
   Int_t           CCInclusiveReco_intCurrent;
   Int_t           CCInclusiveReco_intType;
   Double_t        CCInclusiveReco_E;
   Double_t        CCInclusiveReco_Q2;
   Double_t        CCInclusiveReco_x;
   Double_t        CCInclusiveReco_y;
   Double_t        CCInclusiveReco_W;
   Double_t        CCInclusiveReco_score;
   Double_t        CCInclusiveReco_leptonE[4];
   Double_t        CCInclusiveReco_vtx[4];
   Bool_t          CCInclusiveReco_minos_trk_is_contained;
   Bool_t          CCInclusiveReco_minos_trk_is_ok;
   Bool_t          CCInclusiveReco_minos_used_range;
   Bool_t          CCInclusiveReco_minos_used_curvature;
   Int_t           CCInclusiveReco_minos_trk_end_plane;
   Int_t           CCInclusiveReco_minos_trk_quality;
   Int_t           CCInclusiveReco_r_minos_trk_vtx_plane;
   Int_t           CCInclusiveReco_t_minos_trk_numFSMuons;
   Int_t           CCInclusiveReco_t_minos_trk_primFSLeptonPDG;
   Double_t        CCInclusiveReco_Q2_ccqe;
   Double_t        CCInclusiveReco_minos_trk_bave;
   Double_t        CCInclusiveReco_minos_trk_chi2;
   Double_t        CCInclusiveReco_minos_trk_end_u;
   Double_t        CCInclusiveReco_minos_trk_end_v;
   Double_t        CCInclusiveReco_minos_trk_end_x;
   Double_t        CCInclusiveReco_minos_trk_end_y;
   Double_t        CCInclusiveReco_minos_trk_end_z;
   Double_t        CCInclusiveReco_minos_trk_eqp;
   Double_t        CCInclusiveReco_minos_trk_eqp_qp;
   Double_t        CCInclusiveReco_minos_trk_fit_pass;
   Double_t        CCInclusiveReco_minos_trk_ndf;
   Double_t        CCInclusiveReco_minos_trk_p;
   Double_t        CCInclusiveReco_minos_trk_p_curvature;
   Double_t        CCInclusiveReco_minos_trk_p_range;
   Double_t        CCInclusiveReco_minos_trk_qp;
   Double_t        CCInclusiveReco_minos_trk_vtx_x;
   Double_t        CCInclusiveReco_minos_trk_vtx_y;
   Double_t        CCInclusiveReco_minos_trk_vtx_z;
   Double_t        CCInclusiveReco_nu_energy_recoil;
   Double_t        CCInclusiveReco_r_minos_trk_bdL;
   Double_t        CCInclusiveReco_r_minos_trk_end_dcosx;
   Double_t        CCInclusiveReco_r_minos_trk_end_dcosy;
   Double_t        CCInclusiveReco_r_minos_trk_end_dcosz;
   Double_t        CCInclusiveReco_r_minos_trk_vtx_dcosx;
   Double_t        CCInclusiveReco_r_minos_trk_vtx_dcosy;
   Double_t        CCInclusiveReco_r_minos_trk_vtx_dcosz;
   Double_t        CCInclusiveReco_t_minos_trk_primFSLepMinosInitProjPx;
   Double_t        CCInclusiveReco_t_minos_trk_primFSLepMinosInitProjPy;
   Double_t        CCInclusiveReco_t_minos_trk_primFSLepMinosInitProjPz;
   Double_t        CCInclusiveReco_t_minos_trk_primFSLepMinosInitProjX;
   Double_t        CCInclusiveReco_t_minos_trk_primFSLepMinosInitProjY;
   Double_t        CCInclusiveReco_t_minos_trk_primFSLepMinosInitProjZ;
   Double_t        CCInclusiveReco_t_minos_trk_primFSLepMnvFinalPx;
   Double_t        CCInclusiveReco_t_minos_trk_primFSLepMnvFinalPy;
   Double_t        CCInclusiveReco_t_minos_trk_primFSLepMnvFinalPz;
   Double_t        CCInclusiveReco_t_minos_trk_primFSLepMnvFinalX;
   Double_t        CCInclusiveReco_t_minos_trk_primFSLepMnvFinalY;
   Double_t        CCInclusiveReco_t_minos_trk_primFSLepMnvFinalZ;
   Double_t        CCInclusiveReco_t_minos_trk_primFSLepMnvInitPx;
   Double_t        CCInclusiveReco_t_minos_trk_primFSLepMnvInitPy;
   Double_t        CCInclusiveReco_t_minos_trk_primFSLepMnvInitPz;
   Double_t        CCInclusiveReco_t_minos_trk_primFSLepMnvInitX;
   Double_t        CCInclusiveReco_t_minos_trk_primFSLepMnvInitY;
   Double_t        CCInclusiveReco_t_minos_trk_primFSLepMnvInitZ;
   Double_t        CCInclusiveReco_sys_muon_energy_shift[2];
   Double_t        CCInclusiveReco_sys_muon_qSquared_shift[2];
   Double_t        CCInclusiveReco_sys_muon_wSquared_shift[2];
   Double_t        CCInclusiveReco_sys_muon_xbj_shift[2];
   Double_t        CCInclusiveReco_sys_muon_y_shift[2];
   Double_t        CCInclusiveReco_sys_nu_energy_shift[2];
   Double_t        CCInclusiveReco_sys_recoil_energy_shift[2];
   Double_t        CCInclusiveReco_sys_recoil_qSquared_shift[2];
   Double_t        CCInclusiveReco_sys_recoil_wSquared_shift[2];
   Double_t        CCInclusiveReco_sys_recoil_xbj_shift[2];
   Double_t        CCInclusiveReco_sys_recoil_y_shift[2];
   Double_t        CCInclusiveReco_sys_total_qSquared_shift[2];
   Double_t        CCInclusiveReco_sys_total_wSquared_shift[2];
   Double_t        CCInclusiveReco_sys_total_xbj_shift[2];
   Double_t        CCInclusiveReco_sys_total_y_shift[2];
   Int_t           n_prongs;
   Int_t           prong_nParticles[7];   //[n_prongs]
   Int_t           prong_PDGCode[7];   //[n_prongs]
   Double_t        prong_EnergyWeightedMeanClusterTime[7];   //[n_prongs]
   Double_t        prong_VisibleEnergy[7];   //[n_prongs]
   Double_t        prong_nu_energy_dispersed_cal[7];   //[n_prongs]
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
   Double_t        mc_FSPartPx[28];   //[mc_nFSPart]
   Double_t        mc_FSPartPy[28];   //[mc_nFSPart]
   Double_t        mc_FSPartPz[28];   //[mc_nFSPart]
   Double_t        mc_FSPartE[28];   //[mc_nFSPart]
   Int_t           mc_FSPartPDG[28];   //[mc_nFSPart]
   Int_t           mc_er_nPart;
   Int_t           mc_er_ID[50];   //[mc_er_nPart]
   Int_t           mc_er_status[50];   //[mc_er_nPart]
   Double_t        mc_er_posInNucX[50];   //[mc_er_nPart]
   Double_t        mc_er_posInNucY[50];   //[mc_er_nPart]
   Double_t        mc_er_posInNucZ[50];   //[mc_er_nPart]
   Double_t        mc_er_Px[50];   //[mc_er_nPart]
   Double_t        mc_er_Py[50];   //[mc_er_nPart]
   Double_t        mc_er_Pz[50];   //[mc_er_nPart]
   Double_t        mc_er_E[50];   //[mc_er_nPart]
   Int_t           mc_er_FD[50];   //[mc_er_nPart]
   Int_t           mc_er_LD[50];   //[mc_er_nPart]
   Int_t           mc_er_mother[50];   //[mc_er_nPart]
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
   Double_t        mc_wgt_GENIE[1000];   //[mc_wgt_GENIE_sz]
   Int_t           mc_wgt_Flux_Tertiary_sz;
   Double_t        mc_wgt_Flux_Tertiary[1000];   //[mc_wgt_Flux_Tertiary_sz]
   Int_t           mc_wgt_Flux_BeamFocus_sz;
   Double_t        mc_wgt_Flux_BeamFocus[1000];   //[mc_wgt_Flux_BeamFocus_sz]
   Int_t           mc_wgt_Flux_NA49_sz;
   Double_t        mc_wgt_Flux_NA49[1000];   //[mc_wgt_Flux_NA49_sz]
   Double_t        prong_part_score[7];   //[n_prongs]
   Double_t        prong_part_mass[7];   //[n_prongs]
   Int_t           prong_part_charge[7];   //[n_prongs]
   Int_t           prong_part_pid[7];   //[n_prongs]
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
   TBranch        *b_broken_track_most_us_plane;   //!
   TBranch        *b_n_associated_recoil_subprongs;   //!
   TBranch        *b_n_muon_blobs;   //!
   TBranch        *b_n_recoil_iso2Dblobs;   //!
   TBranch        *b_n_recoil_iso3Dblobs;   //!
   TBranch        *b_n_recoil_subprongs;   //!
   TBranch        *b_n_recoil_vtxblobs;   //!
   TBranch        *b_n_tracks;   //!
   TBranch        *b_n_tracks_non_prim;   //!
   TBranch        *b_n_tracks_prim;   //!
   TBranch        *b_n_tracks_prim_forked;   //!
   TBranch        *b_n_tracks_prim_kinked;   //!
   TBranch        *b_n_vertices_startpoint;   //!
   TBranch        *b_pass_canonical_cut;   //!
   TBranch        *b_pass_minosmatch_cut;   //!
   TBranch        *b_pass_nearccinclusive_cut;   //!
   TBranch        *b_phys_energy_in_road_downstream_nplanes;   //!
   TBranch        *b_phys_energy_in_road_upstream_nplanes;   //!
   TBranch        *b_phys_n_dead_discr_pair;   //!
   TBranch        *b_phys_n_dead_discr_pair_in_prim_track_region;   //!
   TBranch        *b_phys_n_dead_discr_pair_two_mod_downstream_prim_track;   //!
   TBranch        *b_phys_n_dead_discr_pair_two_mod_upstream_prim_vtx;   //!
   TBranch        *b_phys_n_dead_discr_pair_upstream_prim_track_proj;   //!
   TBranch        *b_phys_prim_vertex_planeid;   //!
   TBranch        *b_phys_vertex_is_fiducial;   //!
   TBranch        *b_energy_from_mc;   //!
   TBranch        *b_energy_from_mc_fraction;   //!
   TBranch        *b_energy_from_mc_fraction_of_highest;   //!
   TBranch        *b_muon_blob_energy;   //!
   TBranch        *b_muon_phi;   //!
   TBranch        *b_muon_theta;   //!
   TBranch        *b_muon_thetaX;   //!
   TBranch        *b_muon_thetaY;   //!
   TBranch        *b_nu_topological_energy_recoil;   //!
   TBranch        *b_numi_horn_curr;   //!
   TBranch        *b_numi_pot;   //!
   TBranch        *b_numi_x;   //!
   TBranch        *b_numi_x_width;   //!
   TBranch        *b_numi_y;   //!
   TBranch        *b_numi_y_width;   //!
   TBranch        *b_phys_energy_dispersed;   //!
   TBranch        *b_phys_energy_in_road_downstream;   //!
   TBranch        *b_phys_energy_in_road_upstream;   //!
   TBranch        *b_phys_energy_unattached;   //!
   TBranch        *b_prim_vtx_smallest_opening_angle;   //!
   TBranch        *b_primary_track_minerva_energy;   //!
   TBranch        *b_primary_track_minerva_phi;   //!
   TBranch        *b_primary_track_minerva_theta;   //!
   TBranch        *b_vtxprong_energy_cal;   //!
   TBranch        *b_vtxprong_energy_visible;   //!
   TBranch        *b_vtxprong_energy_visible_isoblobs;   //!
   TBranch        *b_vtxprong_energy_visible_vtxblobs;   //!
   TBranch        *b_plane_id_sz;   //!
   TBranch        *b_plane_id;   //!
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
   TBranch        *b_plane_visible_energy_sz;   //!
   TBranch        *b_plane_visible_energy;   //!
   TBranch        *b_primary_track_minerva_end_position;   //!
   TBranch        *b_primary_track_minerva_start_position;   //!
   TBranch        *b_vtxprong_isoblob_visenergy;   //!
   TBranch        *b_vtxprong_vtxblob_visenergy;   //!
   TBranch        *b_truth_has_physics_event;   //!
   TBranch        *b_truth_pass_CCInclusiveReco;   //!
   TBranch        *b_truth_pass_plausible;   //!
   TBranch        *b_truth_pass_fiducial;   //!
   TBranch        *b_truth_pass_analyzable;   //!
   TBranch        *b_truth_pass_fiducial_apothem;   //!
   TBranch        *b_truth_pass_analyzable_apothem;   //!
   TBranch        *b_truth_pass_ClassifyVertex;   //!
   TBranch        *b_truth_pass_ClassifyCC;   //!
   TBranch        *b_truth_pass_MuonEReconstructed;   //!
   TBranch        *b_truth_pass_MuMinusSignSelection;   //!
   TBranch        *b_truth_pass_MuPlusSignSelection;   //!
   TBranch        *b_truth_pass_ClassifyTracks;   //!
   TBranch        *b_truth_pass_ClassifyRecoil;   //!
   TBranch        *b_truth_pass_Canonical;   //!
   TBranch        *b_truth_reco_vertex_fiducial;   //!
   TBranch        *b_truth_smeared_reco_vertex_fiducial;   //!
   TBranch        *b_truth_pass_nearInclusive;   //!
   TBranch        *b_truth_muon_fraction_accepted;   //!
   TBranch        *b_genie_wgt_n_shifts;   //!
   TBranch        *b_truth_genie_wgt_AGKYxF1pi;   //!
   TBranch        *b_truth_genie_wgt_AhtBY;   //!
   TBranch        *b_truth_genie_wgt_BhtBY;   //!
   TBranch        *b_truth_genie_wgt_CCQEPauliSupViaFK;   //!
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
   TBranch        *b_CCInclusiveReco_nuFlavor;   //!
   TBranch        *b_CCInclusiveReco_nuHelicity;   //!
   TBranch        *b_CCInclusiveReco_intCurrent;   //!
   TBranch        *b_CCInclusiveReco_intType;   //!
   TBranch        *b_CCInclusiveReco_E;   //!
   TBranch        *b_CCInclusiveReco_Q2;   //!
   TBranch        *b_CCInclusiveReco_x;   //!
   TBranch        *b_CCInclusiveReco_y;   //!
   TBranch        *b_CCInclusiveReco_W;   //!
   TBranch        *b_CCInclusiveReco_score;   //!
   TBranch        *b_CCInclusiveReco_leptonE;   //!
   TBranch        *b_CCInclusiveReco_vtx;   //!
   TBranch        *b_CCInclusiveReco_minos_trk_is_contained;   //!
   TBranch        *b_CCInclusiveReco_minos_trk_is_ok;   //!
   TBranch        *b_CCInclusiveReco_minos_used_range;   //!
   TBranch        *b_CCInclusiveReco_minos_used_curvature;   //!
   TBranch        *b_CCInclusiveReco_minos_trk_end_plane;   //!
   TBranch        *b_CCInclusiveReco_minos_trk_quality;   //!
   TBranch        *b_CCInclusiveReco_r_minos_trk_vtx_plane;   //!
   TBranch        *b_CCInclusiveReco_t_minos_trk_numFSMuons;   //!
   TBranch        *b_CCInclusiveReco_t_minos_trk_primFSLeptonPDG;   //!
   TBranch        *b_CCInclusiveReco_Q2_ccqe;   //!
   TBranch        *b_CCInclusiveReco_minos_trk_bave;   //!
   TBranch        *b_CCInclusiveReco_minos_trk_chi2;   //!
   TBranch        *b_CCInclusiveReco_minos_trk_end_u;   //!
   TBranch        *b_CCInclusiveReco_minos_trk_end_v;   //!
   TBranch        *b_CCInclusiveReco_minos_trk_end_x;   //!
   TBranch        *b_CCInclusiveReco_minos_trk_end_y;   //!
   TBranch        *b_CCInclusiveReco_minos_trk_end_z;   //!
   TBranch        *b_CCInclusiveReco_minos_trk_eqp;   //!
   TBranch        *b_CCInclusiveReco_minos_trk_eqp_qp;   //!
   TBranch        *b_CCInclusiveReco_minos_trk_fit_pass;   //!
   TBranch        *b_CCInclusiveReco_minos_trk_ndf;   //!
   TBranch        *b_CCInclusiveReco_minos_trk_p;   //!
   TBranch        *b_CCInclusiveReco_minos_trk_p_curvature;   //!
   TBranch        *b_CCInclusiveReco_minos_trk_p_range;   //!
   TBranch        *b_CCInclusiveReco_minos_trk_qp;   //!
   TBranch        *b_CCInclusiveReco_minos_trk_vtx_x;   //!
   TBranch        *b_CCInclusiveReco_minos_trk_vtx_y;   //!
   TBranch        *b_CCInclusiveReco_minos_trk_vtx_z;   //!
   TBranch        *b_CCInclusiveReco_nu_energy_recoil;   //!
   TBranch        *b_CCInclusiveReco_r_minos_trk_bdL;   //!
   TBranch        *b_CCInclusiveReco_r_minos_trk_end_dcosx;   //!
   TBranch        *b_CCInclusiveReco_r_minos_trk_end_dcosy;   //!
   TBranch        *b_CCInclusiveReco_r_minos_trk_end_dcosz;   //!
   TBranch        *b_CCInclusiveReco_r_minos_trk_vtx_dcosx;   //!
   TBranch        *b_CCInclusiveReco_r_minos_trk_vtx_dcosy;   //!
   TBranch        *b_CCInclusiveReco_r_minos_trk_vtx_dcosz;   //!
   TBranch        *b_CCInclusiveReco_t_minos_trk_primFSLepMinosInitProjPx;   //!
   TBranch        *b_CCInclusiveReco_t_minos_trk_primFSLepMinosInitProjPy;   //!
   TBranch        *b_CCInclusiveReco_t_minos_trk_primFSLepMinosInitProjPz;   //!
   TBranch        *b_CCInclusiveReco_t_minos_trk_primFSLepMinosInitProjX;   //!
   TBranch        *b_CCInclusiveReco_t_minos_trk_primFSLepMinosInitProjY;   //!
   TBranch        *b_CCInclusiveReco_t_minos_trk_primFSLepMinosInitProjZ;   //!
   TBranch        *b_CCInclusiveReco_t_minos_trk_primFSLepMnvFinalPx;   //!
   TBranch        *b_CCInclusiveReco_t_minos_trk_primFSLepMnvFinalPy;   //!
   TBranch        *b_CCInclusiveReco_t_minos_trk_primFSLepMnvFinalPz;   //!
   TBranch        *b_CCInclusiveReco_t_minos_trk_primFSLepMnvFinalX;   //!
   TBranch        *b_CCInclusiveReco_t_minos_trk_primFSLepMnvFinalY;   //!
   TBranch        *b_CCInclusiveReco_t_minos_trk_primFSLepMnvFinalZ;   //!
   TBranch        *b_CCInclusiveReco_t_minos_trk_primFSLepMnvInitPx;   //!
   TBranch        *b_CCInclusiveReco_t_minos_trk_primFSLepMnvInitPy;   //!
   TBranch        *b_CCInclusiveReco_t_minos_trk_primFSLepMnvInitPz;   //!
   TBranch        *b_CCInclusiveReco_t_minos_trk_primFSLepMnvInitX;   //!
   TBranch        *b_CCInclusiveReco_t_minos_trk_primFSLepMnvInitY;   //!
   TBranch        *b_CCInclusiveReco_t_minos_trk_primFSLepMnvInitZ;   //!
   TBranch        *b_CCInclusiveReco_sys_muon_energy_shift;   //!
   TBranch        *b_CCInclusiveReco_sys_muon_qSquared_shift;   //!
   TBranch        *b_CCInclusiveReco_sys_muon_wSquared_shift;   //!
   TBranch        *b_CCInclusiveReco_sys_muon_xbj_shift;   //!
   TBranch        *b_CCInclusiveReco_sys_muon_y_shift;   //!
   TBranch        *b_CCInclusiveReco_sys_nu_energy_shift;   //!
   TBranch        *b_CCInclusiveReco_sys_recoil_energy_shift;   //!
   TBranch        *b_CCInclusiveReco_sys_recoil_qSquared_shift;   //!
   TBranch        *b_CCInclusiveReco_sys_recoil_wSquared_shift;   //!
   TBranch        *b_CCInclusiveReco_sys_recoil_xbj_shift;   //!
   TBranch        *b_CCInclusiveReco_sys_recoil_y_shift;   //!
   TBranch        *b_CCInclusiveReco_sys_total_qSquared_shift;   //!
   TBranch        *b_CCInclusiveReco_sys_total_wSquared_shift;   //!
   TBranch        *b_CCInclusiveReco_sys_total_xbj_shift;   //!
   TBranch        *b_CCInclusiveReco_sys_total_y_shift;   //!
   TBranch        *b_n_prongs;   //!
   TBranch        *b_prong_nParticles;   //!
   TBranch        *b_prong_PDGCode;   //!
   TBranch        *b_prong_EnergyWeightedMeanClusterTime;   //!
   TBranch        *b_prong_VisibleEnergy;   //!
   TBranch        *b_prong_nu_energy_dispersed_cal;   //!
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
   TBranch        *b_prong_part_score;   //!
   TBranch        *b_prong_part_mass;   //!
   TBranch        *b_prong_part_charge;   //!
   TBranch        *b_prong_part_pid;   //!
   TBranch        *b_prong_part_E;   //!
   TBranch        *b_prong_part_pos;   //!

   CCInclusive(TTree *tree=0);
   virtual ~CCInclusive();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif
