/*
================================================================================
Class: CCProtonPi0
    Core Class for CCProtonPi0 Analysis Package
    Includes all data variables for Data Analysis
    Member functions uses member data variables
    
    Outputs .ROOT file filled with 1D or 2D histograms, which is an input file
    for Plotter Class
    
    Uses other classes to define variables or to access specific functions
    
    Main Directory:
        Classes/CCProtonPi0
        
    Usage:
        > main.cpp declares and controls the class
        > See run function Comments
    
    
    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision:  2014_05_13
================================================================================
*/

#ifndef CCProtonPi0_h
#define CCProtonPi0_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h>

#include <vector>
#include <iostream>
#include <fstream>
#include <string>

// Libraries
#include "../../Libraries/PDG_List.h"
#include "../../Libraries/Data_Functions.h"
#include "../../Libraries/HEP_Functions.h"

// Classes
#include "../BinList/BinList.h"
#include "../Muon/Muon.h"
#include "../Proton/Proton.h"
#include "../Pion/Pion.h"


const double mevSq_to_gevSq = pow(10,6);

class CCProtonPi0 {
public :
   // -------------------------------------------------------------------------
   //     Specific Functions
   //--------------------------------------------------------------------------
   
   // -------------------------------------------------------------------------
   //     void run(): Generates a .root file with selected histograms
   //         playlist -> address of the playlist
   //         filename -> file name for the output .root file
   //--------------------------------------------------------------------------
    void run(std::string playlist);
    
    //--------------------------------------------------------------------------
    //  Runtime and CCProtonPi0 Functions
    //      File: CCProtonPi0.cpp
    //--------------------------------------------------------------------------
    void specifyRunTime();
    void initVariables();
    void initHistograms(); // File: initHistograms.cpp
    void closeFiles();
    void openFiles();
    void fillHistograms();
    void write_RootFile();
    void writeReadme();
    bool isBeamEnergyLow(double maxEnergy);
    int findTrueParticle(int targetPDG);
    int countParticles(int targetPDG, bool applyPCut);
    void get_pID_Stats();
    void fillCCProtonPi0();
    
    //--------------------------------------------------------------------------
    //  Muon Specific Functions
    //      File: Muon_Functions.cpp
    //--------------------------------------------------------------------------
    void fillMuonTrue();
    void fillMuonReco();
    //--------------------------------------------------------------------------
    //  Proton Specific Functions
    //      File: Proton_Functions.cpp
    //--------------------------------------------------------------------------
    void fillProtonTrue();
    void fillProtonReco(int ind);
    int findBestProton();
    bool isProtonShort(int ind);
    //--------------------------------------------------------------------------
    //  Pion Specific Functions
    //      File: Pion_Functions.cpp
    //--------------------------------------------------------------------------
    void fillPionTrue();
    void fillPionReco();
    bool isSinglePion();
    bool isNoMeson();


    //--------------------------------------------------------------------------
    //     Default Functions
    //--------------------------------------------------------------------------
    ~CCProtonPi0();
    void Init(std::string playlist, TChain* fChain);
    Int_t GetEntry(Long64_t entry);
    Long64_t LoadTree(Long64_t entry);
    Bool_t Notify();
    void Show(Long64_t entry);
    Int_t Cut(Long64_t entry);
    
    //--------------------------------------------------------------------------
    //     Histograms
    //--------------------------------------------------------------------------
    TFile* f;
    
    // Analysis Variables
    TH1F* beamEnergy_mc;
    TH1F* beamEnergy_reco;
    TH1F* beamEnergy_error;
    TH2F* beamEnergy_reco_mc;
    
    TH1F* q2_mc;
    TH1F* q2_reco;
    TH1F* q2_error;
    TH2F* q2_reco_mc;
    
    TH1F* vertex_z_true;
    TH1F* vertex_z_reco;
    TH1F* vertex_z_error;
    TH2F* vertex_z_reco_mc;
    
    TH1F* int_channel;
    TH2F* vertex_x_y_true;
    TH2F* vertex_x_y_reco;
    TH1F* n_FSParticles;
    TH1F* n_gammas;
    
    TH1F* pID_purity;
    TH1F* pID_efficiency;
    TH1F* pID_piplus;
    TH1F* pID_piminus;
    TH1F* pID_proton;
    TH1F* pID_other;
    
   // -------------------------------------------------------------------------
   //     Analysis Variables
   //--------------------------------------------------------------------------
    bool isDataAnalysis;
    bool isMC;
    bool applyProtonScore;
    bool is_pID_Studies;
    double maxBeamEnergy;
    int max_nFSPart;
    double minProtonScore;
    TVector3 beam_p3;
    Proton proton;
    Muon muon;
    Pion pion;
    
   // -------------------------------------------------------------------------
   //   List of Bins
   //--------------------------------------------------------------------------
    BinList binList;
    
   // -------------------------------------------------------------------------
   //     Files
   //--------------------------------------------------------------------------
    std::string rootDir;
    std::string plotDir;
    std::string readmeFile;
    std::string cutFile;
    std::string channelTag;
    ofstream readme;
    ofstream cutText;

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
   Int_t           Cut_EndPoint_Michel_Exist;
   Int_t           Cut_Muon_None;
   Int_t           Cut_Muon_Score_Low;
   Int_t           Cut_Proton_None;
   Int_t           Cut_Vertex_Michel_Exist;
   Int_t           Cut_Vertex_None;
   Int_t           Cut_Vertex_Not_Fiducial;
   Int_t           Cut_Vertex_Not_Reconstructable;
   Int_t           Cut_Vertex_Null;
   Int_t           Cut_secEndPoint_Michel_Exist;
   Int_t           broken_track_most_us_plane;
   Int_t           n_anchored_long_trk_prongs;
   Int_t           n_anchored_short_trk_prongs;
   Int_t           n_dsp_blob_prongs;
   Int_t           n_iso_blob_prongs;
   Int_t           n_iso_trk_prongs;
   Int_t           n_long_tracks;
   Int_t           n_short_tracks;
   Int_t           n_startpoint_vertices;
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
   Double_t        unattachedExtraE;
   Double_t        vtxBlobExtraE;
   Double_t        vtx_michel_distance;
   Double_t        well_fit_vertex_angle;
   Bool_t          truth_has_physics_event;
   Bool_t          truth_reco_hasGoodObjects;
   Bool_t          truth_reco_isGoodVertex;
   Bool_t          truth_reco_isWellFitVertex;
   Bool_t          truth_reco_isFidVol;
   Bool_t          truth_reco_isFidVol_smeared;
   Bool_t          truth_reco_isMinosMatch;
   Bool_t          truth_reco_isBrokenTrack;
   Bool_t          truth_isSignal;
   Bool_t          truth_isFidVol;
   Bool_t          truth_isPlausible;
   Int_t           truth_N_deltaplus;
   Int_t           truth_N_gamma;
   Int_t           truth_N_muminus;
   Int_t           truth_N_muplus;
   Int_t           truth_N_neutron;
   Int_t           truth_N_other;
   Int_t           truth_N_pi0;
   Int_t           truth_N_piminus;
   Int_t           truth_N_piplus;
   Int_t           truth_N_proton;
   Int_t           truth_muon_charge;
   Int_t           truth_reco_muonCharge;
   Int_t           truth_target_material;
   Int_t           truth_vertex_module;
   Int_t           truth_vertex_plane;
   Double_t        truth_muon_E;
   Double_t        truth_muon_px;
   Double_t        truth_muon_py;
   Double_t        truth_muon_pz;
   Double_t        truth_muon_theta_wrtbeam;
   Int_t           truth_pi0_trackID[20];
   Int_t           truth_proton_trackID[20];
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
   Double_t        truth_pi0_E[20];
   Double_t        truth_pi0_px[20];
   Double_t        truth_pi0_py[20];
   Double_t        truth_pi0_pz[20];
   Double_t        truth_pi0_theta_wrtbeam[20];
   Double_t        truth_proton_E[20];
   Double_t        truth_proton_px[20];
   Double_t        truth_proton_py[20];
   Double_t        truth_proton_pz[20];
   Double_t        truth_proton_theta_wrtbeam[20];
   Int_t           CCProtonPi0Ana_nuFlavor;
   Int_t           CCProtonPi0Ana_nuHelicity;
   Int_t           CCProtonPi0Ana_intCurrent;
   Int_t           CCProtonPi0Ana_intType;
   Double_t        CCProtonPi0Ana_E;
   Double_t        CCProtonPi0Ana_Q2;
   Double_t        CCProtonPi0Ana_x;
   Double_t        CCProtonPi0Ana_y;
   Double_t        CCProtonPi0Ana_W;
   Double_t        CCProtonPi0Ana_score;
   Double_t        CCProtonPi0Ana_leptonE[4];
   Double_t        CCProtonPi0Ana_vtx[4];
   Bool_t          CCProtonPi0Ana_minos_trk_is_contained;
   Bool_t          CCProtonPi0Ana_minos_trk_is_ok;
   Bool_t          CCProtonPi0Ana_minos_used_range;
   Bool_t          CCProtonPi0Ana_minos_used_curvature;
   Int_t           CCProtonPi0Ana_isMuonInsideOD;
   Int_t           CCProtonPi0Ana_minos_trk_end_plane;
   Int_t           CCProtonPi0Ana_minos_trk_quality;
   Int_t           CCProtonPi0Ana_muon_N_minosTracks;
   Int_t           CCProtonPi0Ana_muon_charge;
   Int_t           CCProtonPi0Ana_muon_minervaTrack_types;
   Int_t           CCProtonPi0Ana_muon_minosTrackQuality;
   Int_t           CCProtonPi0Ana_muon_roadUpstreamPlanes;
   Int_t           CCProtonPi0Ana_ntrajMuonProng;
   Int_t           CCProtonPi0Ana_r_minos_trk_vtx_plane;
   Int_t           CCProtonPi0Ana_t_minos_trk_numFSMuons;
   Int_t           CCProtonPi0Ana_t_minos_trk_primFSLeptonPDG;
   Int_t           CCProtonPi0Ana_trajMuonProngPDG;
   Int_t           CCProtonPi0Ana_trajMuonProngPrimary;
   Int_t           CCProtonPi0Ana_vtx_module;
   Int_t           CCProtonPi0Ana_vtx_plane;
   Double_t        CCProtonPi0Ana_endMuonTrajMomentum;
   Double_t        CCProtonPi0Ana_endMuonTrajXPosition;
   Double_t        CCProtonPi0Ana_endMuonTrajYPosition;
   Double_t        CCProtonPi0Ana_endMuonTrajZPosition;
   Double_t        CCProtonPi0Ana_minos_trk_bave;
   Double_t        CCProtonPi0Ana_minos_trk_chi2;
   Double_t        CCProtonPi0Ana_minos_trk_end_u;
   Double_t        CCProtonPi0Ana_minos_trk_end_v;
   Double_t        CCProtonPi0Ana_minos_trk_end_x;
   Double_t        CCProtonPi0Ana_minos_trk_end_y;
   Double_t        CCProtonPi0Ana_minos_trk_end_z;
   Double_t        CCProtonPi0Ana_minos_trk_eqp;
   Double_t        CCProtonPi0Ana_minos_trk_eqp_qp;
   Double_t        CCProtonPi0Ana_minos_trk_fit_pass;
   Double_t        CCProtonPi0Ana_minos_trk_ndf;
   Double_t        CCProtonPi0Ana_minos_trk_p;
   Double_t        CCProtonPi0Ana_minos_trk_p_curvature;
   Double_t        CCProtonPi0Ana_minos_trk_p_range;
   Double_t        CCProtonPi0Ana_minos_trk_qp;
   Double_t        CCProtonPi0Ana_minos_trk_vtx_x;
   Double_t        CCProtonPi0Ana_minos_trk_vtx_y;
   Double_t        CCProtonPi0Ana_minos_trk_vtx_z;
   Double_t        CCProtonPi0Ana_muon_E;
   Double_t        CCProtonPi0Ana_muon_E_shift;
   Double_t        CCProtonPi0Ana_muon_muScore;
   Double_t        CCProtonPi0Ana_muon_p;
   Double_t        CCProtonPi0Ana_muon_px;
   Double_t        CCProtonPi0Ana_muon_py;
   Double_t        CCProtonPi0Ana_muon_pz;
   Double_t        CCProtonPi0Ana_muon_qp;
   Double_t        CCProtonPi0Ana_muon_qpqpe;
   Double_t        CCProtonPi0Ana_muon_roadUpstreamEnergy;
   Double_t        CCProtonPi0Ana_muon_theta;
   Double_t        CCProtonPi0Ana_muon_theta_biasDown;
   Double_t        CCProtonPi0Ana_muon_theta_biasUp;
   Double_t        CCProtonPi0Ana_r_minos_trk_bdL;
   Double_t        CCProtonPi0Ana_r_minos_trk_end_dcosx;
   Double_t        CCProtonPi0Ana_r_minos_trk_end_dcosy;
   Double_t        CCProtonPi0Ana_r_minos_trk_end_dcosz;
   Double_t        CCProtonPi0Ana_r_minos_trk_vtx_dcosx;
   Double_t        CCProtonPi0Ana_r_minos_trk_vtx_dcosy;
   Double_t        CCProtonPi0Ana_r_minos_trk_vtx_dcosz;
   Double_t        CCProtonPi0Ana_t_minos_trk_primFSLepMinosInitProjPx;
   Double_t        CCProtonPi0Ana_t_minos_trk_primFSLepMinosInitProjPy;
   Double_t        CCProtonPi0Ana_t_minos_trk_primFSLepMinosInitProjPz;
   Double_t        CCProtonPi0Ana_t_minos_trk_primFSLepMinosInitProjX;
   Double_t        CCProtonPi0Ana_t_minos_trk_primFSLepMinosInitProjY;
   Double_t        CCProtonPi0Ana_t_minos_trk_primFSLepMinosInitProjZ;
   Double_t        CCProtonPi0Ana_t_minos_trk_primFSLepMnvFinalPx;
   Double_t        CCProtonPi0Ana_t_minos_trk_primFSLepMnvFinalPy;
   Double_t        CCProtonPi0Ana_t_minos_trk_primFSLepMnvFinalPz;
   Double_t        CCProtonPi0Ana_t_minos_trk_primFSLepMnvFinalX;
   Double_t        CCProtonPi0Ana_t_minos_trk_primFSLepMnvFinalY;
   Double_t        CCProtonPi0Ana_t_minos_trk_primFSLepMnvFinalZ;
   Double_t        CCProtonPi0Ana_t_minos_trk_primFSLepMnvInitPx;
   Double_t        CCProtonPi0Ana_t_minos_trk_primFSLepMnvInitPy;
   Double_t        CCProtonPi0Ana_t_minos_trk_primFSLepMnvInitPz;
   Double_t        CCProtonPi0Ana_t_minos_trk_primFSLepMnvInitX;
   Double_t        CCProtonPi0Ana_t_minos_trk_primFSLepMnvInitY;
   Double_t        CCProtonPi0Ana_t_minos_trk_primFSLepMnvInitZ;
   Double_t        CCProtonPi0Ana_trajMuonPhi;
   Double_t        CCProtonPi0Ana_trajMuonProngEnergy;
   Double_t        CCProtonPi0Ana_trajMuonProngMomentum;
   Double_t        CCProtonPi0Ana_trajMuonProngPx;
   Double_t        CCProtonPi0Ana_trajMuonProngPy;
   Double_t        CCProtonPi0Ana_trajMuonProngPz;
   Double_t        CCProtonPi0Ana_trajMuonTheta;
   Double_t        CCProtonPi0Ana_vtx_x;
   Double_t        CCProtonPi0Ana_vtx_y;
   Double_t        CCProtonPi0Ana_vtx_z;
   Int_t           CCProtonPi0Ana_isProtonInsideOD[10];
   Int_t           CCProtonPi0Ana_ntrajProtonProng[10];
   Int_t           CCProtonPi0Ana_proton_kinked[10];
   Int_t           CCProtonPi0Ana_proton_odMatch[10];
   Int_t           CCProtonPi0Ana_proton_trk_pat_history[10];
   Int_t           CCProtonPi0Ana_trajProtonProngPDG[10];
   Int_t           CCProtonPi0Ana_trajProtonProngPrimary[10];
   Double_t        CCProtonPi0Ana_endProtonTrajMomentum[10];
   Double_t        CCProtonPi0Ana_endProtonTrajXPosition[10];
   Double_t        CCProtonPi0Ana_endProtonTrajYPosition[10];
   Double_t        CCProtonPi0Ana_endProtonTrajZPosition[10];
   Double_t        CCProtonPi0Ana_proton_E[10];
   Double_t        CCProtonPi0Ana_proton_chi2_ndf[10];
   Double_t        CCProtonPi0Ana_proton_ekin[10];
   Double_t        CCProtonPi0Ana_proton_endPointX[10];
   Double_t        CCProtonPi0Ana_proton_endPointY[10];
   Double_t        CCProtonPi0Ana_proton_endPointZ[10];
   Double_t        CCProtonPi0Ana_proton_p[10];
   Double_t        CCProtonPi0Ana_proton_p_calCorrection[10];
   Double_t        CCProtonPi0Ana_proton_p_dEdXTool[10];
   Double_t        CCProtonPi0Ana_proton_p_visEnergy[10];
   Double_t        CCProtonPi0Ana_proton_phi[10];
   Double_t        CCProtonPi0Ana_proton_px[10];
   Double_t        CCProtonPi0Ana_proton_py[10];
   Double_t        CCProtonPi0Ana_proton_pz[10];
   Double_t        CCProtonPi0Ana_proton_score[10];
   Double_t        CCProtonPi0Ana_proton_score1[10];
   Double_t        CCProtonPi0Ana_proton_score2[10];
   Double_t        CCProtonPi0Ana_proton_startPointX[10];
   Double_t        CCProtonPi0Ana_proton_startPointY[10];
   Double_t        CCProtonPi0Ana_proton_startPointZ[10];
   Double_t        CCProtonPi0Ana_proton_theta[10];
   Double_t        CCProtonPi0Ana_proton_thetaX[10];
   Double_t        CCProtonPi0Ana_proton_thetaY[10];
   Double_t        CCProtonPi0Ana_trajProtonPhi[10];
   Double_t        CCProtonPi0Ana_trajProtonProngEnergy[10];
   Double_t        CCProtonPi0Ana_trajProtonProngMomentum[10];
   Double_t        CCProtonPi0Ana_trajProtonProngPx[10];
   Double_t        CCProtonPi0Ana_trajProtonProngPy[10];
   Double_t        CCProtonPi0Ana_trajProtonProngPz[10];
   Double_t        CCProtonPi0Ana_trajProtonTheta[10];
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
   Double_t        mc_FSPartPx[21];   //[mc_nFSPart]
   Double_t        mc_FSPartPy[21];   //[mc_nFSPart]
   Double_t        mc_FSPartPz[21];   //[mc_nFSPart]
   Double_t        mc_FSPartE[21];   //[mc_nFSPart]
   Int_t           mc_FSPartPDG[21];   //[mc_nFSPart]
   Int_t           mc_er_nPart;
   Int_t           mc_er_ID[49];   //[mc_er_nPart]
   Int_t           mc_er_status[49];   //[mc_er_nPart]
   Double_t        mc_er_posInNucX[49];   //[mc_er_nPart]
   Double_t        mc_er_posInNucY[49];   //[mc_er_nPart]
   Double_t        mc_er_posInNucZ[49];   //[mc_er_nPart]
   Double_t        mc_er_Px[49];   //[mc_er_nPart]
   Double_t        mc_er_Py[49];   //[mc_er_nPart]
   Double_t        mc_er_Pz[49];   //[mc_er_nPart]
   Double_t        mc_er_E[49];   //[mc_er_nPart]
   Int_t           mc_er_FD[49];   //[mc_er_nPart]
   Int_t           mc_er_LD[49];   //[mc_er_nPart]
   Int_t           mc_er_mother[49];   //[mc_er_nPart]
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
   Int_t           prong_nParticles[7];   //[n_prongs]
   Double_t        prong_part_score[7];   //[n_prongs]
   Double_t        prong_part_mass[7];   //[n_prongs]
   Int_t           prong_part_charge[7];   //[n_prongs]
   Int_t           prong_part_pid[7];   //[n_prongs]
   std::vector< std::vector<double> > *prong_part_E;
   std::vector< std::vector<double> > *prong_part_pos;

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
   TBranch        *b_Cut_EndPoint_Michel_Exist;   //!
   TBranch        *b_Cut_Muon_None;   //!
   TBranch        *b_Cut_Muon_Score_Low;   //!
   TBranch        *b_Cut_Proton_None;   //!
   TBranch        *b_Cut_Vertex_Michel_Exist;   //!
   TBranch        *b_Cut_Vertex_None;   //!
   TBranch        *b_Cut_Vertex_Not_Fiducial;   //!
   TBranch        *b_Cut_Vertex_Not_Reconstructable;   //!
   TBranch        *b_Cut_Vertex_Null;   //!
   TBranch        *b_Cut_secEndPoint_Michel_Exist;   //!
   TBranch        *b_broken_track_most_us_plane;   //!
   TBranch        *b_n_anchored_long_trk_prongs;   //!
   TBranch        *b_n_anchored_short_trk_prongs;   //!
   TBranch        *b_n_dsp_blob_prongs;   //!
   TBranch        *b_n_iso_blob_prongs;   //!
   TBranch        *b_n_iso_trk_prongs;   //!
   TBranch        *b_n_long_tracks;   //!
   TBranch        *b_n_short_tracks;   //!
   TBranch        *b_n_startpoint_vertices;   //!
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
   TBranch        *b_unattachedExtraE;   //!
   TBranch        *b_vtxBlobExtraE;   //!
   TBranch        *b_vtx_michel_distance;   //!
   TBranch        *b_well_fit_vertex_angle;   //!
   TBranch        *b_truth_has_physics_event;   //!
   TBranch        *b_truth_reco_hasGoodObjects;   //!
   TBranch        *b_truth_reco_isGoodVertex;   //!
   TBranch        *b_truth_reco_isWellFitVertex;   //!
   TBranch        *b_truth_reco_isFidVol;   //!
   TBranch        *b_truth_reco_isFidVol_smeared;   //!
   TBranch        *b_truth_reco_isMinosMatch;   //!
   TBranch        *b_truth_reco_isBrokenTrack;   //!
   TBranch        *b_truth_isSignal;   //!
   TBranch        *b_truth_isFidVol;   //!
   TBranch        *b_truth_isPlausible;   //!
   TBranch        *b_truth_N_deltaplus;   //!
   TBranch        *b_truth_N_gamma;   //!
   TBranch        *b_truth_N_muminus;   //!
   TBranch        *b_truth_N_muplus;   //!
   TBranch        *b_truth_N_neutron;   //!
   TBranch        *b_truth_N_other;   //!
   TBranch        *b_truth_N_pi0;   //!
   TBranch        *b_truth_N_piminus;   //!
   TBranch        *b_truth_N_piplus;   //!
   TBranch        *b_truth_N_proton;   //!
   TBranch        *b_truth_muon_charge;   //!
   TBranch        *b_truth_reco_muonCharge;   //!
   TBranch        *b_truth_target_material;   //!
   TBranch        *b_truth_vertex_module;   //!
   TBranch        *b_truth_vertex_plane;   //!
   TBranch        *b_truth_muon_E;   //!
   TBranch        *b_truth_muon_px;   //!
   TBranch        *b_truth_muon_py;   //!
   TBranch        *b_truth_muon_pz;   //!
   TBranch        *b_truth_muon_theta_wrtbeam;   //!
   TBranch        *b_truth_pi0_trackID;   //!
   TBranch        *b_truth_proton_trackID;   //!
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
   TBranch        *b_truth_pi0_E;   //!
   TBranch        *b_truth_pi0_px;   //!
   TBranch        *b_truth_pi0_py;   //!
   TBranch        *b_truth_pi0_pz;   //!
   TBranch        *b_truth_pi0_theta_wrtbeam;   //!
   TBranch        *b_truth_proton_E;   //!
   TBranch        *b_truth_proton_px;   //!
   TBranch        *b_truth_proton_py;   //!
   TBranch        *b_truth_proton_pz;   //!
   TBranch        *b_truth_proton_theta_wrtbeam;   //!
   TBranch        *b_CCProtonPi0Ana_nuFlavor;   //!
   TBranch        *b_CCProtonPi0Ana_nuHelicity;   //!
   TBranch        *b_CCProtonPi0Ana_intCurrent;   //!
   TBranch        *b_CCProtonPi0Ana_intType;   //!
   TBranch        *b_CCProtonPi0Ana_E;   //!
   TBranch        *b_CCProtonPi0Ana_Q2;   //!
   TBranch        *b_CCProtonPi0Ana_x;   //!
   TBranch        *b_CCProtonPi0Ana_y;   //!
   TBranch        *b_CCProtonPi0Ana_W;   //!
   TBranch        *b_CCProtonPi0Ana_score;   //!
   TBranch        *b_CCProtonPi0Ana_leptonE;   //!
   TBranch        *b_CCProtonPi0Ana_vtx;   //!
   TBranch        *b_CCProtonPi0Ana_minos_trk_is_contained;   //!
   TBranch        *b_CCProtonPi0Ana_minos_trk_is_ok;   //!
   TBranch        *b_CCProtonPi0Ana_minos_used_range;   //!
   TBranch        *b_CCProtonPi0Ana_minos_used_curvature;   //!
   TBranch        *b_CCProtonPi0Ana_isMuonInsideOD;   //!
   TBranch        *b_CCProtonPi0Ana_minos_trk_end_plane;   //!
   TBranch        *b_CCProtonPi0Ana_minos_trk_quality;   //!
   TBranch        *b_CCProtonPi0Ana_muon_N_minosTracks;   //!
   TBranch        *b_CCProtonPi0Ana_muon_charge;   //!
   TBranch        *b_CCProtonPi0Ana_muon_minervaTrack_types;   //!
   TBranch        *b_CCProtonPi0Ana_muon_minosTrackQuality;   //!
   TBranch        *b_CCProtonPi0Ana_muon_roadUpstreamPlanes;   //!
   TBranch        *b_CCProtonPi0Ana_ntrajMuonProng;   //!
   TBranch        *b_CCProtonPi0Ana_r_minos_trk_vtx_plane;   //!
   TBranch        *b_CCProtonPi0Ana_t_minos_trk_numFSMuons;   //!
   TBranch        *b_CCProtonPi0Ana_t_minos_trk_primFSLeptonPDG;   //!
   TBranch        *b_CCProtonPi0Ana_trajMuonProngPDG;   //!
   TBranch        *b_CCProtonPi0Ana_trajMuonProngPrimary;   //!
   TBranch        *b_CCProtonPi0Ana_vtx_module;   //!
   TBranch        *b_CCProtonPi0Ana_vtx_plane;   //!
   TBranch        *b_CCProtonPi0Ana_endMuonTrajMomentum;   //!
   TBranch        *b_CCProtonPi0Ana_endMuonTrajXPosition;   //!
   TBranch        *b_CCProtonPi0Ana_endMuonTrajYPosition;   //!
   TBranch        *b_CCProtonPi0Ana_endMuonTrajZPosition;   //!
   TBranch        *b_CCProtonPi0Ana_minos_trk_bave;   //!
   TBranch        *b_CCProtonPi0Ana_minos_trk_chi2;   //!
   TBranch        *b_CCProtonPi0Ana_minos_trk_end_u;   //!
   TBranch        *b_CCProtonPi0Ana_minos_trk_end_v;   //!
   TBranch        *b_CCProtonPi0Ana_minos_trk_end_x;   //!
   TBranch        *b_CCProtonPi0Ana_minos_trk_end_y;   //!
   TBranch        *b_CCProtonPi0Ana_minos_trk_end_z;   //!
   TBranch        *b_CCProtonPi0Ana_minos_trk_eqp;   //!
   TBranch        *b_CCProtonPi0Ana_minos_trk_eqp_qp;   //!
   TBranch        *b_CCProtonPi0Ana_minos_trk_fit_pass;   //!
   TBranch        *b_CCProtonPi0Ana_minos_trk_ndf;   //!
   TBranch        *b_CCProtonPi0Ana_minos_trk_p;   //!
   TBranch        *b_CCProtonPi0Ana_minos_trk_p_curvature;   //!
   TBranch        *b_CCProtonPi0Ana_minos_trk_p_range;   //!
   TBranch        *b_CCProtonPi0Ana_minos_trk_qp;   //!
   TBranch        *b_CCProtonPi0Ana_minos_trk_vtx_x;   //!
   TBranch        *b_CCProtonPi0Ana_minos_trk_vtx_y;   //!
   TBranch        *b_CCProtonPi0Ana_minos_trk_vtx_z;   //!
   TBranch        *b_CCProtonPi0Ana_muon_E;   //!
   TBranch        *b_CCProtonPi0Ana_muon_E_shift;   //!
   TBranch        *b_CCProtonPi0Ana_muon_muScore;   //!
   TBranch        *b_CCProtonPi0Ana_muon_p;   //!
   TBranch        *b_CCProtonPi0Ana_muon_px;   //!
   TBranch        *b_CCProtonPi0Ana_muon_py;   //!
   TBranch        *b_CCProtonPi0Ana_muon_pz;   //!
   TBranch        *b_CCProtonPi0Ana_muon_qp;   //!
   TBranch        *b_CCProtonPi0Ana_muon_qpqpe;   //!
   TBranch        *b_CCProtonPi0Ana_muon_roadUpstreamEnergy;   //!
   TBranch        *b_CCProtonPi0Ana_muon_theta;   //!
   TBranch        *b_CCProtonPi0Ana_muon_theta_biasDown;   //!
   TBranch        *b_CCProtonPi0Ana_muon_theta_biasUp;   //!
   TBranch        *b_CCProtonPi0Ana_r_minos_trk_bdL;   //!
   TBranch        *b_CCProtonPi0Ana_r_minos_trk_end_dcosx;   //!
   TBranch        *b_CCProtonPi0Ana_r_minos_trk_end_dcosy;   //!
   TBranch        *b_CCProtonPi0Ana_r_minos_trk_end_dcosz;   //!
   TBranch        *b_CCProtonPi0Ana_r_minos_trk_vtx_dcosx;   //!
   TBranch        *b_CCProtonPi0Ana_r_minos_trk_vtx_dcosy;   //!
   TBranch        *b_CCProtonPi0Ana_r_minos_trk_vtx_dcosz;   //!
   TBranch        *b_CCProtonPi0Ana_t_minos_trk_primFSLepMinosInitProjPx;   //!
   TBranch        *b_CCProtonPi0Ana_t_minos_trk_primFSLepMinosInitProjPy;   //!
   TBranch        *b_CCProtonPi0Ana_t_minos_trk_primFSLepMinosInitProjPz;   //!
   TBranch        *b_CCProtonPi0Ana_t_minos_trk_primFSLepMinosInitProjX;   //!
   TBranch        *b_CCProtonPi0Ana_t_minos_trk_primFSLepMinosInitProjY;   //!
   TBranch        *b_CCProtonPi0Ana_t_minos_trk_primFSLepMinosInitProjZ;   //!
   TBranch        *b_CCProtonPi0Ana_t_minos_trk_primFSLepMnvFinalPx;   //!
   TBranch        *b_CCProtonPi0Ana_t_minos_trk_primFSLepMnvFinalPy;   //!
   TBranch        *b_CCProtonPi0Ana_t_minos_trk_primFSLepMnvFinalPz;   //!
   TBranch        *b_CCProtonPi0Ana_t_minos_trk_primFSLepMnvFinalX;   //!
   TBranch        *b_CCProtonPi0Ana_t_minos_trk_primFSLepMnvFinalY;   //!
   TBranch        *b_CCProtonPi0Ana_t_minos_trk_primFSLepMnvFinalZ;   //!
   TBranch        *b_CCProtonPi0Ana_t_minos_trk_primFSLepMnvInitPx;   //!
   TBranch        *b_CCProtonPi0Ana_t_minos_trk_primFSLepMnvInitPy;   //!
   TBranch        *b_CCProtonPi0Ana_t_minos_trk_primFSLepMnvInitPz;   //!
   TBranch        *b_CCProtonPi0Ana_t_minos_trk_primFSLepMnvInitX;   //!
   TBranch        *b_CCProtonPi0Ana_t_minos_trk_primFSLepMnvInitY;   //!
   TBranch        *b_CCProtonPi0Ana_t_minos_trk_primFSLepMnvInitZ;   //!
   TBranch        *b_CCProtonPi0Ana_trajMuonPhi;   //!
   TBranch        *b_CCProtonPi0Ana_trajMuonProngEnergy;   //!
   TBranch        *b_CCProtonPi0Ana_trajMuonProngMomentum;   //!
   TBranch        *b_CCProtonPi0Ana_trajMuonProngPx;   //!
   TBranch        *b_CCProtonPi0Ana_trajMuonProngPy;   //!
   TBranch        *b_CCProtonPi0Ana_trajMuonProngPz;   //!
   TBranch        *b_CCProtonPi0Ana_trajMuonTheta;   //!
   TBranch        *b_CCProtonPi0Ana_vtx_x;   //!
   TBranch        *b_CCProtonPi0Ana_vtx_y;   //!
   TBranch        *b_CCProtonPi0Ana_vtx_z;   //!
   TBranch        *b_CCProtonPi0Ana_isProtonInsideOD;   //!
   TBranch        *b_CCProtonPi0Ana_ntrajProtonProng;   //!
   TBranch        *b_CCProtonPi0Ana_proton_kinked;   //!
   TBranch        *b_CCProtonPi0Ana_proton_odMatch;   //!
   TBranch        *b_CCProtonPi0Ana_proton_trk_pat_history;   //!
   TBranch        *b_CCProtonPi0Ana_trajProtonProngPDG;   //!
   TBranch        *b_CCProtonPi0Ana_trajProtonProngPrimary;   //!
   TBranch        *b_CCProtonPi0Ana_endProtonTrajMomentum;   //!
   TBranch        *b_CCProtonPi0Ana_endProtonTrajXPosition;   //!
   TBranch        *b_CCProtonPi0Ana_endProtonTrajYPosition;   //!
   TBranch        *b_CCProtonPi0Ana_endProtonTrajZPosition;   //!
   TBranch        *b_CCProtonPi0Ana_proton_E;   //!
   TBranch        *b_CCProtonPi0Ana_proton_chi2_ndf;   //!
   TBranch        *b_CCProtonPi0Ana_proton_ekin;   //!
   TBranch        *b_CCProtonPi0Ana_proton_endPointX;   //!
   TBranch        *b_CCProtonPi0Ana_proton_endPointY;   //!
   TBranch        *b_CCProtonPi0Ana_proton_endPointZ;   //!
   TBranch        *b_CCProtonPi0Ana_proton_p;   //!
   TBranch        *b_CCProtonPi0Ana_proton_p_calCorrection;   //!
   TBranch        *b_CCProtonPi0Ana_proton_p_dEdXTool;   //!
   TBranch        *b_CCProtonPi0Ana_proton_p_visEnergy;   //!
   TBranch        *b_CCProtonPi0Ana_proton_phi;   //!
   TBranch        *b_CCProtonPi0Ana_proton_px;   //!
   TBranch        *b_CCProtonPi0Ana_proton_py;   //!
   TBranch        *b_CCProtonPi0Ana_proton_pz;   //!
   TBranch        *b_CCProtonPi0Ana_proton_score;   //!
   TBranch        *b_CCProtonPi0Ana_proton_score1;   //!
   TBranch        *b_CCProtonPi0Ana_proton_score2;   //!
   TBranch        *b_CCProtonPi0Ana_proton_startPointX;   //!
   TBranch        *b_CCProtonPi0Ana_proton_startPointY;   //!
   TBranch        *b_CCProtonPi0Ana_proton_startPointZ;   //!
   TBranch        *b_CCProtonPi0Ana_proton_theta;   //!
   TBranch        *b_CCProtonPi0Ana_proton_thetaX;   //!
   TBranch        *b_CCProtonPi0Ana_proton_thetaY;   //!
   TBranch        *b_CCProtonPi0Ana_trajProtonPhi;   //!
   TBranch        *b_CCProtonPi0Ana_trajProtonProngEnergy;   //!
   TBranch        *b_CCProtonPi0Ana_trajProtonProngMomentum;   //!
   TBranch        *b_CCProtonPi0Ana_trajProtonProngPx;   //!
   TBranch        *b_CCProtonPi0Ana_trajProtonProngPy;   //!
   TBranch        *b_CCProtonPi0Ana_trajProtonProngPz;   //!
   TBranch        *b_CCProtonPi0Ana_trajProtonTheta;   //!
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
  
};

#endif
