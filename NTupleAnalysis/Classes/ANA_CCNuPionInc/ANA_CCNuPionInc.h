/*
================================================================================
Class: ANA_CCNuPionInc
    Core Class for CCNuPionInc Analysis Package
    Includes all data variables for Data Analysis
    Member functions uses member data variables
    
    Outputs .ROOT file filled with 1D or 2D histograms, which is an input file
    for Plotter Class
    
    Uses other classes to define variables or to access specific functions
    
    Main Directory:
        Classes/ANA_CCNuPionInc
        
    Usage:
        > main.cpp declares and controls the class
        > See run function Comments
    
    
    
    Last Revision: 2014_02_03
================================================================================
*/

#ifndef ANA_CCNuPionInc_h
#define ANA_CCNuPionInc_h

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


class ANA_CCNuPionInc {
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
    ANA_CCNuPionInc();
    ~ANA_CCNuPionInc();

    void Init(string playlist, TChain* fChain);
    Int_t GetEntry(Long64_t entry);
    Long64_t LoadTree(Long64_t entry);
    Bool_t Notify();
    void Show(Long64_t entry);
    Int_t Cut(Long64_t entry);
    
   // -------------------------------------------------------------------------
   //     Histograms
   //--------------------------------------------------------------------------
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
   Bool_t          isMinosMatchTrack;
   Bool_t          isMinosMatchStub;
   Bool_t          well_fit_vertex;
   Bool_t          isBrokenTrack;
   Int_t           broken_track_most_us_plane;
   Int_t           ddead;
   Int_t           dead;
   Int_t           has_vtx_michel;
   Int_t           n_anchored_long_trk_prongs;
   Int_t           n_anchored_short_trk_prongs;
   Int_t           n_dsp_blob_prongs;
   Int_t           n_iso_blob_prongs;
   Int_t           n_iso_trk_prongs;
   Int_t           n_long_tracks;
   Int_t           n_sepTrk_permutations;
   Int_t           n_short_tracks;
   Int_t           n_startpoint_vertices;
   Int_t           n_twoTrk_permutations;
   Int_t           n_us_muon_clusters;
   Int_t           n_vtx_michel_views;
   Int_t           n_vtx_prongs;
   Int_t           phys_energy_in_road_downstream_nplanes;
   Int_t           phys_energy_in_road_upstream_nplanes;
   Int_t           phys_n_dead_discr_pair;
   Int_t           phys_n_dead_discr_pair_in_prim_track_region;
   Int_t           phys_n_dead_discr_pair_two_mod_downstream_prim_track;
   Int_t           phys_n_dead_discr_pair_two_mod_upstream_prim_vtx;
   Int_t           phys_n_dead_discr_pair_upstream_prim_track_proj;
   Int_t           phys_vertex_is_fiducial;
   Int_t           tdead;
   Int_t           twoTrk_road_planes;
   Int_t           udead;
   Int_t           upstream_plane_num;
   Double_t        dispersedExtraE;
   Double_t        energy_from_mc;
   Double_t        energy_from_mc_fraction;
   Double_t        energy_from_mc_fraction_of_highest;
   Double_t        hadronVisibleE;
   Double_t        muonVisibleE;
   Double_t        muon_phi;
   Double_t        muon_theta;
   Double_t        muon_thetaX;
   Double_t        muon_thetaY;
   Double_t        phys_energy_dispersed;
   Double_t        phys_energy_in_road_downstream;
   Double_t        phys_energy_in_road_upstream;
   Double_t        phys_energy_unattached;
   Double_t        prim_vtx_smallest_opening_angle;
   Double_t        time;
   Double_t        totalIDVisibleE;
   Double_t        totalODVisibleE;
   Double_t        totalVisibleE;
   Double_t        twoTrk_energy_in_road;
   Double_t        unattachedExtraE;
   Double_t        vtxBlobExtraE;
   Double_t        vtx_michel_distance;
   Double_t        well_fit_vertex_angle;
   Int_t           sepTrk_plane[100];
   Int_t           twoTrk_plane[15];
   Double_t        sepTrk_midPtDiff[100];
   Double_t        sepTrk_opening_angle[100];
   Double_t        twoTrk_opening_angle[15];
   Bool_t          truth_has_physics_event;
   Bool_t          truth_isNuSignal;
   Bool_t          truth_isNuBarSignal;
   Bool_t          truth_isTrackablePion;
   Bool_t          truth_isFidVol;
   Bool_t          truth_isPlausible;
   Bool_t          truth_reco_hasGoodObjects;
   Bool_t          truth_reco_isGoodVertex;
   Bool_t          truth_reco_isWellFitVertex;
   Bool_t          truth_reco_isFidVol;
   Bool_t          truth_reco_isFidVol_smeared;
   Bool_t          truth_reco_isMinosMatch;
   Bool_t          truth_reco_isBrokenTrack;
   Int_t           truth_N_gamma;
   Int_t           truth_N_mum;
   Int_t           truth_N_mup;
   Int_t           truth_N_neu;
   Int_t           truth_N_other;
   Int_t           truth_N_pbar;
   Int_t           truth_N_pi0;
   Int_t           truth_N_pim;
   Int_t           truth_N_pip;
   Int_t           truth_N_pro;
   Int_t           truth_Nd_gamma;
   Int_t           truth_Nd_mum;
   Int_t           truth_Nd_mup;
   Int_t           truth_Nd_neu;
   Int_t           truth_Nd_pbar;
   Int_t           truth_Nd_pi0;
   Int_t           truth_Nd_pim;
   Int_t           truth_Nd_pip;
   Int_t           truth_Nd_primaries;
   Int_t           truth_Nd_pro;
   Int_t           truth_mu_charge;
   Int_t           truth_mu_primaryPlanes;
   Int_t           truth_mu_totPlanes;
   Int_t           truth_mu_trackID;
   Int_t           truth_reco_muonCharge;
   Int_t           truth_target_material;
   Int_t           truth_vertex_module;
   Int_t           truth_vertex_plane;
   Double_t        truth_mu_E;
   Double_t        truth_mu_px;
   Double_t        truth_mu_py;
   Double_t        truth_mu_pz;
   Double_t        truth_mu_theta_wrtbeam;
   Int_t           truth_pi_charge[20];
   Int_t           truth_pi_primaryPlanes[20];
   Int_t           truth_pi_totPlanes[20];
   Int_t           truth_pi_trackID[20];
   Int_t           truth_pro_charge[20];
   Int_t           truth_pro_primaryPlanes[20];
   Int_t           truth_pro_totPlanes[20];
   Int_t           truth_pro_trackID[20];
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
   Double_t        truth_pi_E[20];
   Double_t        truth_pi_px[20];
   Double_t        truth_pi_py[20];
   Double_t        truth_pi_pz[20];
   Double_t        truth_pi_theta_wrtbeam[20];
   Double_t        truth_pro_E[20];
   Double_t        truth_pro_px[20];
   Double_t        truth_pro_py[20];
   Double_t        truth_pro_pz[20];
   Double_t        truth_pro_theta_wrtbeam[20];
   Int_t           CCNuPionInc_nuFlavor;
   Int_t           CCNuPionInc_nuHelicity;
   Int_t           CCNuPionInc_intCurrent;
   Int_t           CCNuPionInc_intType;
   Double_t        CCNuPionInc_E;
   Double_t        CCNuPionInc_Q2;
   Double_t        CCNuPionInc_x;
   Double_t        CCNuPionInc_y;
   Double_t        CCNuPionInc_W;
   Double_t        CCNuPionInc_score;
   Double_t        CCNuPionInc_leptonE[4];
   Double_t        CCNuPionInc_vtx[4];
   Bool_t          CCNuPionInc_minos_trk_is_contained;
   Bool_t          CCNuPionInc_minos_trk_is_ok;
   Bool_t          CCNuPionInc_minos_used_range;
   Bool_t          CCNuPionInc_minos_used_curvature;
   Int_t           CCNuPionInc_hadron_number;
   Int_t           CCNuPionInc_minos_trk_end_plane;
   Int_t           CCNuPionInc_minos_trk_quality;
   Int_t           CCNuPionInc_muon_N_minosTracks;
   Int_t           CCNuPionInc_muon_minervaTrack_types;
   Int_t           CCNuPionInc_muon_minosTrackQuality;
   Int_t           CCNuPionInc_muon_roadUpstreamPlanes;
   Int_t           CCNuPionInc_r_minos_trk_vtx_plane;
   Int_t           CCNuPionInc_t_minos_trk_numFSMuons;
   Int_t           CCNuPionInc_t_minos_trk_primFSLeptonPDG;
   Int_t           CCNuPionInc_vtx_module;
   Int_t           CCNuPionInc_vtx_plane;
   Double_t        CCNuPionInc_hadron_recoil;
   Double_t        CCNuPionInc_hadron_recoil_CCInc;
   Double_t        CCNuPionInc_hadron_recoil_CCInc_em;
   Double_t        CCNuPionInc_hadron_recoil_CCInc_kaon;
   Double_t        CCNuPionInc_hadron_recoil_CCInc_neutron;
   Double_t        CCNuPionInc_hadron_recoil_CCInc_other;
   Double_t        CCNuPionInc_hadron_recoil_CCInc_pion;
   Double_t        CCNuPionInc_hadron_recoil_CCInc_proton;
   Double_t        CCNuPionInc_hadron_recoil_CCInc_xtalk;
   Double_t        CCNuPionInc_minos_trk_bave;
   Double_t        CCNuPionInc_minos_trk_chi2;
   Double_t        CCNuPionInc_minos_trk_end_u;
   Double_t        CCNuPionInc_minos_trk_end_v;
   Double_t        CCNuPionInc_minos_trk_end_x;
   Double_t        CCNuPionInc_minos_trk_end_y;
   Double_t        CCNuPionInc_minos_trk_end_z;
   Double_t        CCNuPionInc_minos_trk_eqp;
   Double_t        CCNuPionInc_minos_trk_eqp_qp;
   Double_t        CCNuPionInc_minos_trk_fit_pass;
   Double_t        CCNuPionInc_minos_trk_ndf;
   Double_t        CCNuPionInc_minos_trk_p;
   Double_t        CCNuPionInc_minos_trk_p_curvature;
   Double_t        CCNuPionInc_minos_trk_p_range;
   Double_t        CCNuPionInc_minos_trk_qp;
   Double_t        CCNuPionInc_minos_trk_vtx_x;
   Double_t        CCNuPionInc_minos_trk_vtx_y;
   Double_t        CCNuPionInc_minos_trk_vtx_z;
   Double_t        CCNuPionInc_muon_E;
   Double_t        CCNuPionInc_muon_E_shift;
   Double_t        CCNuPionInc_muon_muScore;
   Double_t        CCNuPionInc_muon_p;
   Double_t        CCNuPionInc_muon_px;
   Double_t        CCNuPionInc_muon_py;
   Double_t        CCNuPionInc_muon_pz;
   Double_t        CCNuPionInc_muon_qp;
   Double_t        CCNuPionInc_muon_qpqpe;
   Double_t        CCNuPionInc_muon_roadUpstreamEnergy;
   Double_t        CCNuPionInc_muon_theta;
   Double_t        CCNuPionInc_muon_theta_biasDown;
   Double_t        CCNuPionInc_muon_theta_biasUp;
   Double_t        CCNuPionInc_r_minos_trk_bdL;
   Double_t        CCNuPionInc_r_minos_trk_end_dcosx;
   Double_t        CCNuPionInc_r_minos_trk_end_dcosy;
   Double_t        CCNuPionInc_r_minos_trk_end_dcosz;
   Double_t        CCNuPionInc_r_minos_trk_vtx_dcosx;
   Double_t        CCNuPionInc_r_minos_trk_vtx_dcosy;
   Double_t        CCNuPionInc_r_minos_trk_vtx_dcosz;
   Double_t        CCNuPionInc_t_minos_trk_primFSLepMinosInitProjPx;
   Double_t        CCNuPionInc_t_minos_trk_primFSLepMinosInitProjPy;
   Double_t        CCNuPionInc_t_minos_trk_primFSLepMinosInitProjPz;
   Double_t        CCNuPionInc_t_minos_trk_primFSLepMinosInitProjX;
   Double_t        CCNuPionInc_t_minos_trk_primFSLepMinosInitProjY;
   Double_t        CCNuPionInc_t_minos_trk_primFSLepMinosInitProjZ;
   Double_t        CCNuPionInc_t_minos_trk_primFSLepMnvFinalPx;
   Double_t        CCNuPionInc_t_minos_trk_primFSLepMnvFinalPy;
   Double_t        CCNuPionInc_t_minos_trk_primFSLepMnvFinalPz;
   Double_t        CCNuPionInc_t_minos_trk_primFSLepMnvFinalX;
   Double_t        CCNuPionInc_t_minos_trk_primFSLepMnvFinalY;
   Double_t        CCNuPionInc_t_minos_trk_primFSLepMnvFinalZ;
   Double_t        CCNuPionInc_t_minos_trk_primFSLepMnvInitPx;
   Double_t        CCNuPionInc_t_minos_trk_primFSLepMnvInitPy;
   Double_t        CCNuPionInc_t_minos_trk_primFSLepMnvInitPz;
   Double_t        CCNuPionInc_t_minos_trk_primFSLepMnvInitX;
   Double_t        CCNuPionInc_t_minos_trk_primFSLepMnvInitY;
   Double_t        CCNuPionInc_t_minos_trk_primFSLepMnvInitZ;
   Double_t        CCNuPionInc_vtx_x;
   Double_t        CCNuPionInc_vtx_y;
   Double_t        CCNuPionInc_vtx_z;
   Int_t           CCNuPionInc_hadron_1stTrackPatRec[10];
   Int_t           CCNuPionInc_hadron_endMichel_bgmodule[10];
   Int_t           CCNuPionInc_hadron_endMichel_category[10];
   Int_t           CCNuPionInc_hadron_endMichel_edmodule[10];
   Int_t           CCNuPionInc_hadron_endMichel_ndigits[10];
   Int_t           CCNuPionInc_hadron_endMichel_nmodules[10];
   Int_t           CCNuPionInc_hadron_endMichel_nplanes[10];
   Int_t           CCNuPionInc_hadron_endMichel_tm_parentpdg[10];
   Int_t           CCNuPionInc_hadron_endMichel_tm_pdg[10];
   Int_t           CCNuPionInc_hadron_endMichel_tm_primarypdg[10];
   Int_t           CCNuPionInc_hadron_hasEndpointMichel[10];
   Int_t           CCNuPionInc_hadron_hasSecondaryMichel[10];
   Int_t           CCNuPionInc_hadron_isDsECAL[10];
   Int_t           CCNuPionInc_hadron_isExiting[10];
   Int_t           CCNuPionInc_hadron_isForked[10];
   Int_t           CCNuPionInc_hadron_isHCAL[10];
   Int_t           CCNuPionInc_hadron_isKinked[10];
   Int_t           CCNuPionInc_hadron_isNuclTargs[10];
   Int_t           CCNuPionInc_hadron_isODMatch[10];
   Int_t           CCNuPionInc_hadron_isSideECAL[10];
   Int_t           CCNuPionInc_hadron_isTracker[10];
   Int_t           CCNuPionInc_hadron_piFit_fails[10];
   Int_t           CCNuPionInc_hadron_piFit_isMM[10];
   Int_t           CCNuPionInc_hadron_piFit_ok[10];
   Int_t           CCNuPionInc_hadron_piFit_range[10];
   Int_t           CCNuPionInc_hadron_proFit_fails[10];
   Int_t           CCNuPionInc_hadron_proFit_isMM[10];
   Int_t           CCNuPionInc_hadron_proFit_ok[10];
   Int_t           CCNuPionInc_hadron_proFit_range[10];
   Int_t           CCNuPionInc_hadron_secMichel_bgmodule[10];
   Int_t           CCNuPionInc_hadron_secMichel_category[10];
   Int_t           CCNuPionInc_hadron_secMichel_edmodule[10];
   Int_t           CCNuPionInc_hadron_secMichel_ndigits[10];
   Int_t           CCNuPionInc_hadron_secMichel_nmodules[10];
   Int_t           CCNuPionInc_hadron_secMichel_nplanes[10];
   Int_t           CCNuPionInc_hadron_tm_PDGCode[10];
   Int_t           CCNuPionInc_hadron_tm_PDGCodeDaughter[10];
   Int_t           CCNuPionInc_hadron_tm_destructCode[10];
   Int_t           CCNuPionInc_hadron_tm_endpt_pdg[10];
   Int_t           CCNuPionInc_hadron_tm_firstInelasticCode[10];
   Int_t           CCNuPionInc_hadron_tm_isTruthMatched[10];
   Int_t           CCNuPionInc_hadron_tm_trackID[10];
   Int_t           CCNuPionInc_hadron_trackNodes[10];
   Double_t        CCNuPionInc_hadron_calE_CCInc[10];
   Double_t        CCNuPionInc_hadron_coneE[10];
   Double_t        CCNuPionInc_hadron_endAvdEdX[10];
   Double_t        CCNuPionInc_hadron_endMichel_bgz[10];
   Double_t        CCNuPionInc_hadron_endMichel_bgzl[10];
   Double_t        CCNuPionInc_hadron_endMichel_distance[10];
   Double_t        CCNuPionInc_hadron_endMichel_edz[10];
   Double_t        CCNuPionInc_hadron_endMichel_edzl[10];
   Double_t        CCNuPionInc_hadron_endMichel_energy[10];
   Double_t        CCNuPionInc_hadron_endMichel_energy_uncorrected[10];
   Double_t        CCNuPionInc_hadron_endMichel_firedFraction[10];
   Double_t        CCNuPionInc_hadron_endMichel_mcfrac[10];
   Double_t        CCNuPionInc_hadron_endMichel_slice_energy[10];
   Double_t        CCNuPionInc_hadron_endMichel_time_diff[10];
   Double_t        CCNuPionInc_hadron_endMichel_tm_otherE[10];
   Double_t        CCNuPionInc_hadron_endpointE[10];
   Double_t        CCNuPionInc_hadron_matRange[10];
   Double_t        CCNuPionInc_hadron_piFit_chi2[10];
   Double_t        CCNuPionInc_hadron_piFit_lastAllScore1[10];
   Double_t        CCNuPionInc_hadron_piFit_lastHalfScore1[10];
   Double_t        CCNuPionInc_hadron_piFit_lastSixScore1[10];
   Double_t        CCNuPionInc_hadron_piFit_score1[10];
   Double_t        CCNuPionInc_hadron_piFit_score1_BetheBloch_biasDown[10];
   Double_t        CCNuPionInc_hadron_piFit_score1_BetheBloch_biasUp[10];
   Double_t        CCNuPionInc_hadron_piFit_score1_Birks[10];
   Double_t        CCNuPionInc_hadron_piFit_score1_MEU_biasDown[10];
   Double_t        CCNuPionInc_hadron_piFit_score1_MEU_biasUp[10];
   Double_t        CCNuPionInc_hadron_piFit_score1_Mass_biasDown[10];
   Double_t        CCNuPionInc_hadron_piFit_score1_Mass_biasUp[10];
   Double_t        CCNuPionInc_hadron_piFit_score2[10];
   Double_t        CCNuPionInc_hadron_piFit_zDiff[10];
   Double_t        CCNuPionInc_hadron_pimFit_SPID[10];
   Double_t        CCNuPionInc_hadron_pimFit_newChi2[10];
   Double_t        CCNuPionInc_hadron_pion_E[10];
   Double_t        CCNuPionInc_hadron_pion_E_BetheBloch_biasDown[10];
   Double_t        CCNuPionInc_hadron_pion_E_BetheBloch_biasUp[10];
   Double_t        CCNuPionInc_hadron_pion_E_Birks[10];
   Double_t        CCNuPionInc_hadron_pion_E_MEU_biasDown[10];
   Double_t        CCNuPionInc_hadron_pion_E_MEU_biasUp[10];
   Double_t        CCNuPionInc_hadron_pion_E_Mass_biasDown[10];
   Double_t        CCNuPionInc_hadron_pion_E_Mass_biasUp[10];
   Double_t        CCNuPionInc_hadron_pion_openingAngle[10];
   Double_t        CCNuPionInc_hadron_pion_px[10];
   Double_t        CCNuPionInc_hadron_pion_py[10];
   Double_t        CCNuPionInc_hadron_pion_pz[10];
   Double_t        CCNuPionInc_hadron_pion_theta[10];
   Double_t        CCNuPionInc_hadron_pion_theta_biasDown[10];
   Double_t        CCNuPionInc_hadron_pion_theta_biasUp[10];
   Double_t        CCNuPionInc_hadron_pipFit_SPID[10];
   Double_t        CCNuPionInc_hadron_pipFit_newChi2[10];
   Double_t        CCNuPionInc_hadron_proFit_SPID[10];
   Double_t        CCNuPionInc_hadron_proFit_chi2[10];
   Double_t        CCNuPionInc_hadron_proFit_newChi2[10];
   Double_t        CCNuPionInc_hadron_proFit_score1[10];
   Double_t        CCNuPionInc_hadron_proFit_score2[10];
   Double_t        CCNuPionInc_hadron_proFit_zDiff[10];
   Double_t        CCNuPionInc_hadron_rawRange[10];
   Double_t        CCNuPionInc_hadron_secMichel_bgz[10];
   Double_t        CCNuPionInc_hadron_secMichel_bgzl[10];
   Double_t        CCNuPionInc_hadron_secMichel_distance[10];
   Double_t        CCNuPionInc_hadron_secMichel_edz[10];
   Double_t        CCNuPionInc_hadron_secMichel_edzl[10];
   Double_t        CCNuPionInc_hadron_secMichel_energy[10];
   Double_t        CCNuPionInc_hadron_secMichel_energy_uncorrected[10];
   Double_t        CCNuPionInc_hadron_secMichel_firedFraction[10];
   Double_t        CCNuPionInc_hadron_secMichel_slice_energy[10];
   Double_t        CCNuPionInc_hadron_secMichel_time_diff[10];
   Double_t        CCNuPionInc_hadron_tm_beginMomentum[10];
   Double_t        CCNuPionInc_hadron_tm_endMomentum[10];
   Double_t        CCNuPionInc_hadron_tm_endpt_fraction[10];
   Double_t        CCNuPionInc_hadron_tm_endpt_otherE[10];
   Double_t        CCNuPionInc_hadron_tm_firstInelasticDistance[10];
   Double_t        CCNuPionInc_hadron_tm_firstInelasticEndP[10];
   Double_t        CCNuPionInc_hadron_tm_fraction[10];
   Double_t        CCNuPionInc_hadron_tm_otherE[10];
   Double_t        CCNuPionInc_hadron_tm_pathlength[10];
   Double_t        CCNuPionInc_hadron_visibleE[10];
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
   Double_t        mc_FSPartPx[19];   //[mc_nFSPart]
   Double_t        mc_FSPartPy[19];   //[mc_nFSPart]
   Double_t        mc_FSPartPz[19];   //[mc_nFSPart]
   Double_t        mc_FSPartE[19];   //[mc_nFSPart]
   Int_t           mc_FSPartPDG[19];   //[mc_nFSPart]
   Int_t           mc_er_nPart;
   Int_t           mc_er_ID[43];   //[mc_er_nPart]
   Int_t           mc_er_status[43];   //[mc_er_nPart]
   Double_t        mc_er_posInNucX[43];   //[mc_er_nPart]
   Double_t        mc_er_posInNucY[43];   //[mc_er_nPart]
   Double_t        mc_er_posInNucZ[43];   //[mc_er_nPart]
   Double_t        mc_er_Px[43];   //[mc_er_nPart]
   Double_t        mc_er_Py[43];   //[mc_er_nPart]
   Double_t        mc_er_Pz[43];   //[mc_er_nPart]
   Double_t        mc_er_E[43];   //[mc_er_nPart]
   Int_t           mc_er_FD[43];   //[mc_er_nPart]
   Int_t           mc_er_LD[43];   //[mc_er_nPart]
   Int_t           mc_er_mother[43];   //[mc_er_nPart]
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
   TBranch        *b_isMinosMatchTrack;   //!
   TBranch        *b_isMinosMatchStub;   //!
   TBranch        *b_well_fit_vertex;   //!
   TBranch        *b_isBrokenTrack;   //!
   TBranch        *b_broken_track_most_us_plane;   //!
   TBranch        *b_ddead;   //!
   TBranch        *b_dead;   //!
   TBranch        *b_has_vtx_michel;   //!
   TBranch        *b_n_anchored_long_trk_prongs;   //!
   TBranch        *b_n_anchored_short_trk_prongs;   //!
   TBranch        *b_n_dsp_blob_prongs;   //!
   TBranch        *b_n_iso_blob_prongs;   //!
   TBranch        *b_n_iso_trk_prongs;   //!
   TBranch        *b_n_long_tracks;   //!
   TBranch        *b_n_sepTrk_permutations;   //!
   TBranch        *b_n_short_tracks;   //!
   TBranch        *b_n_startpoint_vertices;   //!
   TBranch        *b_n_twoTrk_permutations;   //!
   TBranch        *b_n_us_muon_clusters;   //!
   TBranch        *b_n_vtx_michel_views;   //!
   TBranch        *b_n_vtx_prongs;   //!
   TBranch        *b_phys_energy_in_road_downstream_nplanes;   //!
   TBranch        *b_phys_energy_in_road_upstream_nplanes;   //!
   TBranch        *b_phys_n_dead_discr_pair;   //!
   TBranch        *b_phys_n_dead_discr_pair_in_prim_track_region;   //!
   TBranch        *b_phys_n_dead_discr_pair_two_mod_downstream_prim_track;   //!
   TBranch        *b_phys_n_dead_discr_pair_two_mod_upstream_prim_vtx;   //!
   TBranch        *b_phys_n_dead_discr_pair_upstream_prim_track_proj;   //!
   TBranch        *b_phys_vertex_is_fiducial;   //!
   TBranch        *b_tdead;   //!
   TBranch        *b_twoTrk_road_planes;   //!
   TBranch        *b_udead;   //!
   TBranch        *b_upstream_plane_num;   //!
   TBranch        *b_dispersedExtraE;   //!
   TBranch        *b_energy_from_mc;   //!
   TBranch        *b_energy_from_mc_fraction;   //!
   TBranch        *b_energy_from_mc_fraction_of_highest;   //!
   TBranch        *b_hadronVisibleE;   //!
   TBranch        *b_muonVisibleE;   //!
   TBranch        *b_muon_phi;   //!
   TBranch        *b_muon_theta;   //!
   TBranch        *b_muon_thetaX;   //!
   TBranch        *b_muon_thetaY;   //!
   TBranch        *b_phys_energy_dispersed;   //!
   TBranch        *b_phys_energy_in_road_downstream;   //!
   TBranch        *b_phys_energy_in_road_upstream;   //!
   TBranch        *b_phys_energy_unattached;   //!
   TBranch        *b_prim_vtx_smallest_opening_angle;   //!
   TBranch        *b_time;   //!
   TBranch        *b_totalIDVisibleE;   //!
   TBranch        *b_totalODVisibleE;   //!
   TBranch        *b_totalVisibleE;   //!
   TBranch        *b_twoTrk_energy_in_road;   //!
   TBranch        *b_unattachedExtraE;   //!
   TBranch        *b_vtxBlobExtraE;   //!
   TBranch        *b_vtx_michel_distance;   //!
   TBranch        *b_well_fit_vertex_angle;   //!
   TBranch        *b_sepTrk_plane;   //!
   TBranch        *b_twoTrk_plane;   //!
   TBranch        *b_sepTrk_midPtDiff;   //!
   TBranch        *b_sepTrk_opening_angle;   //!
   TBranch        *b_twoTrk_opening_angle;   //!
   TBranch        *b_truth_has_physics_event;   //!
   TBranch        *b_truth_isNuSignal;   //!
   TBranch        *b_truth_isNuBarSignal;   //!
   TBranch        *b_truth_isTrackablePion;   //!
   TBranch        *b_truth_isFidVol;   //!
   TBranch        *b_truth_isPlausible;   //!
   TBranch        *b_truth_reco_hasGoodObjects;   //!
   TBranch        *b_truth_reco_isGoodVertex;   //!
   TBranch        *b_truth_reco_isWellFitVertex;   //!
   TBranch        *b_truth_reco_isFidVol;   //!
   TBranch        *b_truth_reco_isFidVol_smeared;   //!
   TBranch        *b_truth_reco_isMinosMatch;   //!
   TBranch        *b_truth_reco_isBrokenTrack;   //!
   TBranch        *b_truth_N_gamma;   //!
   TBranch        *b_truth_N_mum;   //!
   TBranch        *b_truth_N_mup;   //!
   TBranch        *b_truth_N_neu;   //!
   TBranch        *b_truth_N_other;   //!
   TBranch        *b_truth_N_pbar;   //!
   TBranch        *b_truth_N_pi0;   //!
   TBranch        *b_truth_N_pim;   //!
   TBranch        *b_truth_N_pip;   //!
   TBranch        *b_truth_N_pro;   //!
   TBranch        *b_truth_Nd_gamma;   //!
   TBranch        *b_truth_Nd_mum;   //!
   TBranch        *b_truth_Nd_mup;   //!
   TBranch        *b_truth_Nd_neu;   //!
   TBranch        *b_truth_Nd_pbar;   //!
   TBranch        *b_truth_Nd_pi0;   //!
   TBranch        *b_truth_Nd_pim;   //!
   TBranch        *b_truth_Nd_pip;   //!
   TBranch        *b_truth_Nd_primaries;   //!
   TBranch        *b_truth_Nd_pro;   //!
   TBranch        *b_truth_mu_charge;   //!
   TBranch        *b_truth_mu_primaryPlanes;   //!
   TBranch        *b_truth_mu_totPlanes;   //!
   TBranch        *b_truth_mu_trackID;   //!
   TBranch        *b_truth_reco_muonCharge;   //!
   TBranch        *b_truth_target_material;   //!
   TBranch        *b_truth_vertex_module;   //!
   TBranch        *b_truth_vertex_plane;   //!
   TBranch        *b_truth_mu_E;   //!
   TBranch        *b_truth_mu_px;   //!
   TBranch        *b_truth_mu_py;   //!
   TBranch        *b_truth_mu_pz;   //!
   TBranch        *b_truth_mu_theta_wrtbeam;   //!
   TBranch        *b_truth_pi_charge;   //!
   TBranch        *b_truth_pi_primaryPlanes;   //!
   TBranch        *b_truth_pi_totPlanes;   //!
   TBranch        *b_truth_pi_trackID;   //!
   TBranch        *b_truth_pro_charge;   //!
   TBranch        *b_truth_pro_primaryPlanes;   //!
   TBranch        *b_truth_pro_totPlanes;   //!
   TBranch        *b_truth_pro_trackID;   //!
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
   TBranch        *b_truth_pi_E;   //!
   TBranch        *b_truth_pi_px;   //!
   TBranch        *b_truth_pi_py;   //!
   TBranch        *b_truth_pi_pz;   //!
   TBranch        *b_truth_pi_theta_wrtbeam;   //!
   TBranch        *b_truth_pro_E;   //!
   TBranch        *b_truth_pro_px;   //!
   TBranch        *b_truth_pro_py;   //!
   TBranch        *b_truth_pro_pz;   //!
   TBranch        *b_truth_pro_theta_wrtbeam;   //!
   TBranch        *b_CCNuPionInc_nuFlavor;   //!
   TBranch        *b_CCNuPionInc_nuHelicity;   //!
   TBranch        *b_CCNuPionInc_intCurrent;   //!
   TBranch        *b_CCNuPionInc_intType;   //!
   TBranch        *b_CCNuPionInc_E;   //!
   TBranch        *b_CCNuPionInc_Q2;   //!
   TBranch        *b_CCNuPionInc_x;   //!
   TBranch        *b_CCNuPionInc_y;   //!
   TBranch        *b_CCNuPionInc_W;   //!
   TBranch        *b_CCNuPionInc_score;   //!
   TBranch        *b_CCNuPionInc_leptonE;   //!
   TBranch        *b_CCNuPionInc_vtx;   //!
   TBranch        *b_CCNuPionInc_minos_trk_is_contained;   //!
   TBranch        *b_CCNuPionInc_minos_trk_is_ok;   //!
   TBranch        *b_CCNuPionInc_minos_used_range;   //!
   TBranch        *b_CCNuPionInc_minos_used_curvature;   //!
   TBranch        *b_CCNuPionInc_hadron_number;   //!
   TBranch        *b_CCNuPionInc_minos_trk_end_plane;   //!
   TBranch        *b_CCNuPionInc_minos_trk_quality;   //!
   TBranch        *b_CCNuPionInc_muon_N_minosTracks;   //!
   TBranch        *b_CCNuPionInc_muon_minervaTrack_types;   //!
   TBranch        *b_CCNuPionInc_muon_minosTrackQuality;   //!
   TBranch        *b_CCNuPionInc_muon_roadUpstreamPlanes;   //!
   TBranch        *b_CCNuPionInc_r_minos_trk_vtx_plane;   //!
   TBranch        *b_CCNuPionInc_t_minos_trk_numFSMuons;   //!
   TBranch        *b_CCNuPionInc_t_minos_trk_primFSLeptonPDG;   //!
   TBranch        *b_CCNuPionInc_vtx_module;   //!
   TBranch        *b_CCNuPionInc_vtx_plane;   //!
   TBranch        *b_CCNuPionInc_hadron_recoil;   //!
   TBranch        *b_CCNuPionInc_hadron_recoil_CCInc;   //!
   TBranch        *b_CCNuPionInc_hadron_recoil_CCInc_em;   //!
   TBranch        *b_CCNuPionInc_hadron_recoil_CCInc_kaon;   //!
   TBranch        *b_CCNuPionInc_hadron_recoil_CCInc_neutron;   //!
   TBranch        *b_CCNuPionInc_hadron_recoil_CCInc_other;   //!
   TBranch        *b_CCNuPionInc_hadron_recoil_CCInc_pion;   //!
   TBranch        *b_CCNuPionInc_hadron_recoil_CCInc_proton;   //!
   TBranch        *b_CCNuPionInc_hadron_recoil_CCInc_xtalk;   //!
   TBranch        *b_CCNuPionInc_minos_trk_bave;   //!
   TBranch        *b_CCNuPionInc_minos_trk_chi2;   //!
   TBranch        *b_CCNuPionInc_minos_trk_end_u;   //!
   TBranch        *b_CCNuPionInc_minos_trk_end_v;   //!
   TBranch        *b_CCNuPionInc_minos_trk_end_x;   //!
   TBranch        *b_CCNuPionInc_minos_trk_end_y;   //!
   TBranch        *b_CCNuPionInc_minos_trk_end_z;   //!
   TBranch        *b_CCNuPionInc_minos_trk_eqp;   //!
   TBranch        *b_CCNuPionInc_minos_trk_eqp_qp;   //!
   TBranch        *b_CCNuPionInc_minos_trk_fit_pass;   //!
   TBranch        *b_CCNuPionInc_minos_trk_ndf;   //!
   TBranch        *b_CCNuPionInc_minos_trk_p;   //!
   TBranch        *b_CCNuPionInc_minos_trk_p_curvature;   //!
   TBranch        *b_CCNuPionInc_minos_trk_p_range;   //!
   TBranch        *b_CCNuPionInc_minos_trk_qp;   //!
   TBranch        *b_CCNuPionInc_minos_trk_vtx_x;   //!
   TBranch        *b_CCNuPionInc_minos_trk_vtx_y;   //!
   TBranch        *b_CCNuPionInc_minos_trk_vtx_z;   //!
   TBranch        *b_CCNuPionInc_muon_E;   //!
   TBranch        *b_CCNuPionInc_muon_E_shift;   //!
   TBranch        *b_CCNuPionInc_muon_muScore;   //!
   TBranch        *b_CCNuPionInc_muon_p;   //!
   TBranch        *b_CCNuPionInc_muon_px;   //!
   TBranch        *b_CCNuPionInc_muon_py;   //!
   TBranch        *b_CCNuPionInc_muon_pz;   //!
   TBranch        *b_CCNuPionInc_muon_qp;   //!
   TBranch        *b_CCNuPionInc_muon_qpqpe;   //!
   TBranch        *b_CCNuPionInc_muon_roadUpstreamEnergy;   //!
   TBranch        *b_CCNuPionInc_muon_theta;   //!
   TBranch        *b_CCNuPionInc_muon_theta_biasDown;   //!
   TBranch        *b_CCNuPionInc_muon_theta_biasUp;   //!
   TBranch        *b_CCNuPionInc_r_minos_trk_bdL;   //!
   TBranch        *b_CCNuPionInc_r_minos_trk_end_dcosx;   //!
   TBranch        *b_CCNuPionInc_r_minos_trk_end_dcosy;   //!
   TBranch        *b_CCNuPionInc_r_minos_trk_end_dcosz;   //!
   TBranch        *b_CCNuPionInc_r_minos_trk_vtx_dcosx;   //!
   TBranch        *b_CCNuPionInc_r_minos_trk_vtx_dcosy;   //!
   TBranch        *b_CCNuPionInc_r_minos_trk_vtx_dcosz;   //!
   TBranch        *b_CCNuPionInc_t_minos_trk_primFSLepMinosInitProjPx;   //!
   TBranch        *b_CCNuPionInc_t_minos_trk_primFSLepMinosInitProjPy;   //!
   TBranch        *b_CCNuPionInc_t_minos_trk_primFSLepMinosInitProjPz;   //!
   TBranch        *b_CCNuPionInc_t_minos_trk_primFSLepMinosInitProjX;   //!
   TBranch        *b_CCNuPionInc_t_minos_trk_primFSLepMinosInitProjY;   //!
   TBranch        *b_CCNuPionInc_t_minos_trk_primFSLepMinosInitProjZ;   //!
   TBranch        *b_CCNuPionInc_t_minos_trk_primFSLepMnvFinalPx;   //!
   TBranch        *b_CCNuPionInc_t_minos_trk_primFSLepMnvFinalPy;   //!
   TBranch        *b_CCNuPionInc_t_minos_trk_primFSLepMnvFinalPz;   //!
   TBranch        *b_CCNuPionInc_t_minos_trk_primFSLepMnvFinalX;   //!
   TBranch        *b_CCNuPionInc_t_minos_trk_primFSLepMnvFinalY;   //!
   TBranch        *b_CCNuPionInc_t_minos_trk_primFSLepMnvFinalZ;   //!
   TBranch        *b_CCNuPionInc_t_minos_trk_primFSLepMnvInitPx;   //!
   TBranch        *b_CCNuPionInc_t_minos_trk_primFSLepMnvInitPy;   //!
   TBranch        *b_CCNuPionInc_t_minos_trk_primFSLepMnvInitPz;   //!
   TBranch        *b_CCNuPionInc_t_minos_trk_primFSLepMnvInitX;   //!
   TBranch        *b_CCNuPionInc_t_minos_trk_primFSLepMnvInitY;   //!
   TBranch        *b_CCNuPionInc_t_minos_trk_primFSLepMnvInitZ;   //!
   TBranch        *b_CCNuPionInc_vtx_x;   //!
   TBranch        *b_CCNuPionInc_vtx_y;   //!
   TBranch        *b_CCNuPionInc_vtx_z;   //!
   TBranch        *b_CCNuPionInc_hadron_1stTrackPatRec;   //!
   TBranch        *b_CCNuPionInc_hadron_endMichel_bgmodule;   //!
   TBranch        *b_CCNuPionInc_hadron_endMichel_category;   //!
   TBranch        *b_CCNuPionInc_hadron_endMichel_edmodule;   //!
   TBranch        *b_CCNuPionInc_hadron_endMichel_ndigits;   //!
   TBranch        *b_CCNuPionInc_hadron_endMichel_nmodules;   //!
   TBranch        *b_CCNuPionInc_hadron_endMichel_nplanes;   //!
   TBranch        *b_CCNuPionInc_hadron_endMichel_tm_parentpdg;   //!
   TBranch        *b_CCNuPionInc_hadron_endMichel_tm_pdg;   //!
   TBranch        *b_CCNuPionInc_hadron_endMichel_tm_primarypdg;   //!
   TBranch        *b_CCNuPionInc_hadron_hasEndpointMichel;   //!
   TBranch        *b_CCNuPionInc_hadron_hasSecondaryMichel;   //!
   TBranch        *b_CCNuPionInc_hadron_isDsECAL;   //!
   TBranch        *b_CCNuPionInc_hadron_isExiting;   //!
   TBranch        *b_CCNuPionInc_hadron_isForked;   //!
   TBranch        *b_CCNuPionInc_hadron_isHCAL;   //!
   TBranch        *b_CCNuPionInc_hadron_isKinked;   //!
   TBranch        *b_CCNuPionInc_hadron_isNuclTargs;   //!
   TBranch        *b_CCNuPionInc_hadron_isODMatch;   //!
   TBranch        *b_CCNuPionInc_hadron_isSideECAL;   //!
   TBranch        *b_CCNuPionInc_hadron_isTracker;   //!
   TBranch        *b_CCNuPionInc_hadron_piFit_fails;   //!
   TBranch        *b_CCNuPionInc_hadron_piFit_isMM;   //!
   TBranch        *b_CCNuPionInc_hadron_piFit_ok;   //!
   TBranch        *b_CCNuPionInc_hadron_piFit_range;   //!
   TBranch        *b_CCNuPionInc_hadron_proFit_fails;   //!
   TBranch        *b_CCNuPionInc_hadron_proFit_isMM;   //!
   TBranch        *b_CCNuPionInc_hadron_proFit_ok;   //!
   TBranch        *b_CCNuPionInc_hadron_proFit_range;   //!
   TBranch        *b_CCNuPionInc_hadron_secMichel_bgmodule;   //!
   TBranch        *b_CCNuPionInc_hadron_secMichel_category;   //!
   TBranch        *b_CCNuPionInc_hadron_secMichel_edmodule;   //!
   TBranch        *b_CCNuPionInc_hadron_secMichel_ndigits;   //!
   TBranch        *b_CCNuPionInc_hadron_secMichel_nmodules;   //!
   TBranch        *b_CCNuPionInc_hadron_secMichel_nplanes;   //!
   TBranch        *b_CCNuPionInc_hadron_tm_PDGCode;   //!
   TBranch        *b_CCNuPionInc_hadron_tm_PDGCodeDaughter;   //!
   TBranch        *b_CCNuPionInc_hadron_tm_destructCode;   //!
   TBranch        *b_CCNuPionInc_hadron_tm_endpt_pdg;   //!
   TBranch        *b_CCNuPionInc_hadron_tm_firstInelasticCode;   //!
   TBranch        *b_CCNuPionInc_hadron_tm_isTruthMatched;   //!
   TBranch        *b_CCNuPionInc_hadron_tm_trackID;   //!
   TBranch        *b_CCNuPionInc_hadron_trackNodes;   //!
   TBranch        *b_CCNuPionInc_hadron_calE_CCInc;   //!
   TBranch        *b_CCNuPionInc_hadron_coneE;   //!
   TBranch        *b_CCNuPionInc_hadron_endAvdEdX;   //!
   TBranch        *b_CCNuPionInc_hadron_endMichel_bgz;   //!
   TBranch        *b_CCNuPionInc_hadron_endMichel_bgzl;   //!
   TBranch        *b_CCNuPionInc_hadron_endMichel_distance;   //!
   TBranch        *b_CCNuPionInc_hadron_endMichel_edz;   //!
   TBranch        *b_CCNuPionInc_hadron_endMichel_edzl;   //!
   TBranch        *b_CCNuPionInc_hadron_endMichel_energy;   //!
   TBranch        *b_CCNuPionInc_hadron_endMichel_energy_uncorrected;   //!
   TBranch        *b_CCNuPionInc_hadron_endMichel_firedFraction;   //!
   TBranch        *b_CCNuPionInc_hadron_endMichel_mcfrac;   //!
   TBranch        *b_CCNuPionInc_hadron_endMichel_slice_energy;   //!
   TBranch        *b_CCNuPionInc_hadron_endMichel_time_diff;   //!
   TBranch        *b_CCNuPionInc_hadron_endMichel_tm_otherE;   //!
   TBranch        *b_CCNuPionInc_hadron_endpointE;   //!
   TBranch        *b_CCNuPionInc_hadron_matRange;   //!
   TBranch        *b_CCNuPionInc_hadron_piFit_chi2;   //!
   TBranch        *b_CCNuPionInc_hadron_piFit_lastAllScore1;   //!
   TBranch        *b_CCNuPionInc_hadron_piFit_lastHalfScore1;   //!
   TBranch        *b_CCNuPionInc_hadron_piFit_lastSixScore1;   //!
   TBranch        *b_CCNuPionInc_hadron_piFit_score1;   //!
   TBranch        *b_CCNuPionInc_hadron_piFit_score1_BetheBloch_biasDown;   //!
   TBranch        *b_CCNuPionInc_hadron_piFit_score1_BetheBloch_biasUp;   //!
   TBranch        *b_CCNuPionInc_hadron_piFit_score1_Birks;   //!
   TBranch        *b_CCNuPionInc_hadron_piFit_score1_MEU_biasDown;   //!
   TBranch        *b_CCNuPionInc_hadron_piFit_score1_MEU_biasUp;   //!
   TBranch        *b_CCNuPionInc_hadron_piFit_score1_Mass_biasDown;   //!
   TBranch        *b_CCNuPionInc_hadron_piFit_score1_Mass_biasUp;   //!
   TBranch        *b_CCNuPionInc_hadron_piFit_score2;   //!
   TBranch        *b_CCNuPionInc_hadron_piFit_zDiff;   //!
   TBranch        *b_CCNuPionInc_hadron_pimFit_SPID;   //!
   TBranch        *b_CCNuPionInc_hadron_pimFit_newChi2;   //!
   TBranch        *b_CCNuPionInc_hadron_pion_E;   //!
   TBranch        *b_CCNuPionInc_hadron_pion_E_BetheBloch_biasDown;   //!
   TBranch        *b_CCNuPionInc_hadron_pion_E_BetheBloch_biasUp;   //!
   TBranch        *b_CCNuPionInc_hadron_pion_E_Birks;   //!
   TBranch        *b_CCNuPionInc_hadron_pion_E_MEU_biasDown;   //!
   TBranch        *b_CCNuPionInc_hadron_pion_E_MEU_biasUp;   //!
   TBranch        *b_CCNuPionInc_hadron_pion_E_Mass_biasDown;   //!
   TBranch        *b_CCNuPionInc_hadron_pion_E_Mass_biasUp;   //!
   TBranch        *b_CCNuPionInc_hadron_pion_openingAngle;   //!
   TBranch        *b_CCNuPionInc_hadron_pion_px;   //!
   TBranch        *b_CCNuPionInc_hadron_pion_py;   //!
   TBranch        *b_CCNuPionInc_hadron_pion_pz;   //!
   TBranch        *b_CCNuPionInc_hadron_pion_theta;   //!
   TBranch        *b_CCNuPionInc_hadron_pion_theta_biasDown;   //!
   TBranch        *b_CCNuPionInc_hadron_pion_theta_biasUp;   //!
   TBranch        *b_CCNuPionInc_hadron_pipFit_SPID;   //!
   TBranch        *b_CCNuPionInc_hadron_pipFit_newChi2;   //!
   TBranch        *b_CCNuPionInc_hadron_proFit_SPID;   //!
   TBranch        *b_CCNuPionInc_hadron_proFit_chi2;   //!
   TBranch        *b_CCNuPionInc_hadron_proFit_newChi2;   //!
   TBranch        *b_CCNuPionInc_hadron_proFit_score1;   //!
   TBranch        *b_CCNuPionInc_hadron_proFit_score2;   //!
   TBranch        *b_CCNuPionInc_hadron_proFit_zDiff;   //!
   TBranch        *b_CCNuPionInc_hadron_rawRange;   //!
   TBranch        *b_CCNuPionInc_hadron_secMichel_bgz;   //!
   TBranch        *b_CCNuPionInc_hadron_secMichel_bgzl;   //!
   TBranch        *b_CCNuPionInc_hadron_secMichel_distance;   //!
   TBranch        *b_CCNuPionInc_hadron_secMichel_edz;   //!
   TBranch        *b_CCNuPionInc_hadron_secMichel_edzl;   //!
   TBranch        *b_CCNuPionInc_hadron_secMichel_energy;   //!
   TBranch        *b_CCNuPionInc_hadron_secMichel_energy_uncorrected;   //!
   TBranch        *b_CCNuPionInc_hadron_secMichel_firedFraction;   //!
   TBranch        *b_CCNuPionInc_hadron_secMichel_slice_energy;   //!
   TBranch        *b_CCNuPionInc_hadron_secMichel_time_diff;   //!
   TBranch        *b_CCNuPionInc_hadron_tm_beginMomentum;   //!
   TBranch        *b_CCNuPionInc_hadron_tm_endMomentum;   //!
   TBranch        *b_CCNuPionInc_hadron_tm_endpt_fraction;   //!
   TBranch        *b_CCNuPionInc_hadron_tm_endpt_otherE;   //!
   TBranch        *b_CCNuPionInc_hadron_tm_firstInelasticDistance;   //!
   TBranch        *b_CCNuPionInc_hadron_tm_firstInelasticEndP;   //!
   TBranch        *b_CCNuPionInc_hadron_tm_fraction;   //!
   TBranch        *b_CCNuPionInc_hadron_tm_otherE;   //!
   TBranch        *b_CCNuPionInc_hadron_tm_pathlength;   //!
   TBranch        *b_CCNuPionInc_hadron_visibleE;   //!
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

   ANA_CCNuPionInc(TTree *tree=0);
   virtual ~ANA_CCNuPionInc();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif
