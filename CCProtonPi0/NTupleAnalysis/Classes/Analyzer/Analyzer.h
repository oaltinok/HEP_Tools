/*
================================================================================
Class: Analyzer
    Core Class for Analyzer Analysis Package
    Includes all data variables for Data Analysis
    Member functions uses member data variables
    
    Outputs .ROOT file filled with 1D or 2D histograms, which is an input file
    for Plotter Class
    
    Uses other classes to define variables or to access specific functions
    
    Main Directory:
        Classes/Analyzer
        
    Usage:
        > main.cpp declares and controls the class
        > See run function Comments
    
    
    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision:  2015_04_20
================================================================================
*/

#ifndef Analyzer_h
#define Analyzer_h

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
#include <cstdlib>

// Libraries
#include "../../Libraries/PDG_List.h"
#include "../../Libraries/Folder_List.h"
#include "../../Libraries/Data_Functions.h"
#include "../../Libraries/HEP_Functions.h"

// Classes
#include "../NTupleAnalysis/NTupleAnalysis.h"
#include "../BinList/BinList.h"
#include "../CutList/CutList.h"
#include "../PIDTool/PIDTool.h"
#include "../Muon/Muon.h"
#include "../Proton/Proton.h"
#include "../Pion/Pion.h"
#include "../BackgroundTool/BackgroundTool.h"

class Analyzer : public NTupleAnalysis{
    
public :
  
    Analyzer(int nMode);
    ~Analyzer();

   // -------------------------------------------------------------------------
   //     void run(): Generates a .root file with selected histograms
   //         playlist -> address of the playlist
   //---------------------------------------------------------------------------
    void run(std::string playlist);
    
    //--------------------------------------------------------------------------
    //  Initialization Functions
    //      File: initFunctions.cpp
    //--------------------------------------------------------------------------
    void initInteraction();
    void initHistograms(); 
    
    //--------------------------------------------------------------------------
    //  Runtime and Analyzer Functions
    //      File: Analyzer.cpp
    //--------------------------------------------------------------------------
    bool analyzeEvent();
    void fillData();
    void specifyRunTime();
    void closeTextFiles();
    void openTextFiles();
    void fillHistograms();
    void write_RootFile();
    void writeReadme();
    double getCutEfficiency(double nSig, double effBase);
    double getCutPurity(double nSig, double nEvents);
    void fill_mc_w();
    void writeScanFile();
    
    //--------------------------------------------------------------------------
    //  Interaction Specific Functions
    //--------------------------------------------------------------------------
    void getPi0Family();
    bool getCutStatistics();
    void writeFSParticle4P(Long64_t nEntry);
    void fillInteractionTrue();
    void fillInteractionReco();
    double calcDeltaInvariantMass();
    
    //--------------------------------------------------------------------------
    //  Muon Specific Functions
    //      File: MuonFunctions.cpp
    //--------------------------------------------------------------------------
    void fillMuonTrue();
    void fillMuonReco();
    
    //--------------------------------------------------------------------------
    //  Proton Specific Functions
    //      File: ProtonFunctions.cpp
    //--------------------------------------------------------------------------
    void findTrueProton();
    void findRecoProton();
    void fillProtonTrue();
    void fillProtonReco();
    
    //--------------------------------------------------------------------------
    //  Pion Specific Functions
    //      File: PionFunctions.cpp
    //--------------------------------------------------------------------------
    void fillPionTrue();
    void fillPionReco();
    
    //--------------------------------------------------------------------------
    //  Default Functions
    //      File: DefaultFunctions.cpp
    //--------------------------------------------------------------------------
    void Init(string playlist, TChain* fChain);
    Int_t GetEntry(Long64_t entry);
    Long64_t LoadTree(Long64_t entry);
    Bool_t Notify();
    void Show(Long64_t entry);
    Int_t CutEvent(Long64_t entry);
    
    //--------------------------------------------------------------------------
    //     Histograms
    //--------------------------------------------------------------------------
    TFile* f;
    
    // Analysis Variables
    TH1D* beamEnergy_mc;
    TH1D* beamEnergy_reco;
    TH1D* beamEnergy_error;
    TH2D* beamEnergy_reco_mc;
    
    TH1D* beamEnergyCal_mc;
    TH1D* beamEnergyCal_reco;
    TH1D* beamEnergyCal_error;
    TH2D* beamEnergyCal_reco_mc;
    
    TH2D* beamEnergy_beamEnergyCal;
    
    TH1D* deltaInvMass_mc;
    TH1D* deltaInvMass_reco;
    TH1D* deltaInvMass_error;
    TH2D* deltaInvMass_reco_mc;
    
    TH1D* nProngs_hist;
    
    TH1D* pFilter_Status;
    TH1D* pFilter_RejectedEnergy;
    
    TH1D* q2_mc;
    TH1D* q2_reco;
    TH1D* q2_error;
    TH2D* q2_reco_mc;
    
    TH1D* w_mc;
    TH1D* w_reco;
    TH1D* w_error;
    TH2D* w_reco_mc;
    TH1D* wSq_reco;
    
    TH1D* vertex_count;
    
    TH1D* vertex_z_true;
    TH1D* vertex_z_reco;
    TH1D* vertex_z_error;
    TH2D* vertex_z_reco_mc;
   
    TH1D* int_channel;
    TH2D* vertex_x_y_true;
    TH2D* vertex_x_y_reco;
    
    TH1D* mc_w_DIS;
    TH1D* mc_w_RES;
    TH1D* mc_w_CCQE;
    
    TH1D* final_mc_w_DIS;
    TH1D* final_mc_w_RES;
    TH1D* final_mc_w_CCQE;
    
    TH1D* status_Pi0;
    TH1D* status_Pi0_Mother;
    TH1D* status_Pi0_GrandMother;
    
    TH1D* E_Unused_afterReco;
    TH1D* E_Used_afterReco;
    TH1D* time_AllClusters;
    
    TH1D* total_E;
    TH2D* total_E_neutrinoE;
    
    // Michel Study
    bool isMichelStudy;
    double N_trueMichel_before;
    double N_trueMichel_after;
    double N_trueMichel_afterAll;
    double N_noMichel_before;
    double N_noMichel_after;
    TH1D* N_michelElectrons;
    TH1D* michelMuon_P[4];
    TH1D* michelMuon_end_dist_vtx[4];
    TH1D* michelMuon_length[4];
    TH1D* michelMuon_Z[4];
    TH2D* michelMuon_X_Y[4];
    TH1D* michelMuon_Z_vtx[4];
    TH1D* michelPion_P[4];
    TH1D* michelPion_begin_dist_vtx[4];
    TH1D* michelPion_length[4];
    TH1D* michelElectron_E[5];
    TH2D* michelPion_length_dist_vtx[4];
    TH2D* michelMuon_dist_michelPion_length[4];
  
    // -------------------------------------------------------------------------
    // Cut Histograms
    // -------------------------------------------------------------------------
    // Common
    TH1D* hCut_vertexCount;
    TH1D* hCut_nProngs;
    
    // Topology Dependent
    TH1D* hCut_1Prong_Michel;
    TH1D* hCut_2Prong_Michel;
    TH1D* hCut_1Prong_eVis_nuclearTarget;
    TH1D* hCut_2Prong_eVis_nuclearTarget;
    TH1D* hCut_1Prong_eVis_other;
    TH1D* hCut_2Prong_eVis_other;
    TH1D* hCut_1Prong_pi0invMass;
    TH1D* hCut_2Prong_pi0invMass;
    TH1D* hCut_1Prong_gamma1ConvDist;
    TH1D* hCut_2Prong_gamma1ConvDist;
    TH1D* hCut_1Prong_gamma2ConvDist;
    TH1D* hCut_2Prong_gamma2ConvDist;
    TH1D* hCut_1Prong_neutrinoE;
    TH1D* hCut_2Prong_neutrinoE;
    TH1D* hCut_1Prong_UnusedE;
    TH1D* hCut_2Prong_UnusedE;

    
    // 2 Prong Specific
    TH1D* hCut_pIDDiff;
    TH1D* hCut_protonScore_LLR;
    TH1D* hCut_deltaInvMass;

    
   // -------------------------------------------------------------------------
   //     Analysis Variables
   //--------------------------------------------------------------------------
    Muon muon;
    Proton proton;
    Pion pion;
    PIDTool pIDTool;
    BackgroundTool bckgTool;
    CutList cutList;
    bool isDataAnalysis;
    bool analyze_NoProtonEvents;
    bool hasParticleTruthInfo;
    bool is_pID_Studies;
    bool applyProtonScore;
    bool applyPhotonDistance;
    bool applyBeamEnergy;
    bool applyQSq;
    bool applyUnusedE;
    bool writeFSParticleMomentum;
    bool isPassedAllCuts;
    bool applyMaxEvents;
    bool applyDeltaInvMass;
    
    double latest_ScanID;
    int max_nFSPart;
    int indRecoProton;
    int indTrueProton;
    double debug_Counter;
    double pID_KE_Limit;
    double minProtonScore_LLR;
    double minPIDDiff;
    double minPhotonDistance;
    double min_Pi0_invMass;
    double max_Pi0_invMass;
    double min_Delta_invMass;
    double max_Delta_invMass;
    double max_QSq;
    double max_beamEnergy;
    double maxUnusedE;
    double nMaxEvents;
    double SENTINEL;
    TVector3 beam_p3;

   
    string scanFileName;
    std::vector<double> PDG_pi0_Mother;
    std::vector<double> PDG_pi0_GrandMother;
    
   // -------------------------------------------------------------------------
   //   Branch Activation Control
   //--------------------------------------------------------------------------
    bool m_ActivateMC;
    bool m_ActivateInteraction;
    bool m_ActivatePi0; 
    
   // -------------------------------------------------------------------------
   //   List of Bins
   //--------------------------------------------------------------------------
    BinList binList;
    
    
   // -------------------------------------------------------------------------
   //     Files
   //--------------------------------------------------------------------------
    string rootDir;
    string readmeFile;
    string channelTag;
    string failFile[nTopologies];
       
    ofstream readme;
    ofstream failText[nTopologies];

    ofstream roundupText;
    ifstream DSTFileList;

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
   Bool_t          gamma1_isGoodDirection;
   Bool_t          gamma1_isGoodPosition;
   Bool_t          gamma1_isGoodBlob;
   Bool_t          gamma2_isGoodDirection;
   Bool_t          gamma2_isGoodPosition;
   Bool_t          gamma2_isGoodBlob;
   Int_t           Cut_ConeBlobs;
   Int_t           Cut_EndPoint_Michel_Exist;
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
   Int_t           Cut_VtxBlob;
   Int_t           Cut_nProngs;
   Int_t           Cut_secEndPoint_Michel_Exist;
   Int_t           anglescan_ncand;
   Int_t           anglescan_ncandx;
   Int_t           blob_ndof_1;
   Int_t           blob_ndof_2;
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
   Int_t           n_vtx_michel_views;
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
   Double_t        blob_fval_1;
   Double_t        blob_fval_2;
   Double_t        endpoint_michel_distance;
   Double_t        energyUnused_afterReco;
   Double_t        energyUsed_afterReco;
   Double_t        energy_from_mc;
   Double_t        energy_from_mc_fraction;
   Double_t        energy_from_mc_fraction_of_highest;
   Double_t        evis_ECAL;
   Double_t        evis_HCAL;
   Double_t        evis_NuclearTarget;
   Double_t        evis_TotalExceptNuclearTarget;
   Double_t        evis_Tracker;
   Double_t        evis_nearvtx;
   Double_t        evis_total;
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
   Double_t        gamma1_dist_vtx;
   Double_t        gamma1_evis_ecal;
   Double_t        gamma1_evis_hcal;
   Double_t        gamma1_evis_scal;
   Double_t        gamma1_evis_trkr;
   Double_t        gamma1_phi;
   Double_t        gamma1_px;
   Double_t        gamma1_py;
   Double_t        gamma1_pz;
   Double_t        gamma1_score;
   Double_t        gamma1_theta;
   Double_t        gamma1_time;
   Double_t        gamma2_E;
   Double_t        gamma2_dEdx;
   Double_t        gamma2_dist_vtx;
   Double_t        gamma2_evis_ecal;
   Double_t        gamma2_evis_hcal;
   Double_t        gamma2_evis_scal;
   Double_t        gamma2_evis_trkr;
   Double_t        gamma2_phi;
   Double_t        gamma2_px;
   Double_t        gamma2_py;
   Double_t        gamma2_pz;
   Double_t        gamma2_score;
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
   Double_t        prong_showerScore;
   Double_t        reco_eventID;
   Double_t        time;
   Double_t        totalIDVisibleE;
   Double_t        totalODVisibleE;
   Double_t        totalVisibleE;
   Double_t        vtx_michel_distance;
   Int_t           g1dedx_cluster_occupancy_sz;
   Int_t           g1dedx_cluster_occupancy[6];   //[g1dedx_cluster_occupancy_sz]
   Int_t           g2dedx_cluster_occupancy_sz;
   Int_t           g2dedx_cluster_occupancy[6];   //[g2dedx_cluster_occupancy_sz]
   Int_t           Vertex_energy_radii_sz;
   Double_t        Vertex_energy_radii[7];   //[Vertex_energy_radii_sz]
   Int_t           blob_cluster_energy1_sz;
   Double_t        blob_cluster_energy1[6];   //[blob_cluster_energy1_sz]
   Int_t           blob_cluster_energy2_sz;
   Double_t        blob_cluster_energy2[6];   //[blob_cluster_energy2_sz]
   Int_t           g1dedx_cluster_energy_sz;
   Double_t        g1dedx_cluster_energy[6];   //[g1dedx_cluster_energy_sz]
   Int_t           g1dedx_rev_cluster_energy_sz;
   Double_t        g1dedx_rev_cluster_energy[56];   //[g1dedx_rev_cluster_energy_sz]
   Int_t           g2dedx_cluster_energy_sz;
   Double_t        g2dedx_cluster_energy[6];   //[g2dedx_cluster_energy_sz]
   Int_t           g2dedx_rev_cluster_energy_sz;
   Double_t        g2dedx_rev_cluster_energy[22];   //[g2dedx_rev_cluster_energy_sz]
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
   Bool_t          truth_isSignal;
   Bool_t          truth_isSignal_1Pi0;
   Bool_t          truth_isSignal_2Gamma;
   Bool_t          truth_isFidVol;
   Bool_t          truth_AnalyzeEvent;
   Bool_t          truth_isBckg_QELike;
   Bool_t          truth_isBckg_SinglePiPlus;
   Bool_t          truth_isBckg_SinglePiMinus;
   Bool_t          truth_isBckg_MultiPion;
   Bool_t          truth_isBckg_MultiPiZero;
   Bool_t          truth_isBckg_Other;
   Bool_t          truth_isBckg_withAntiMuon;
   Bool_t          truth_isBckg_withMichel;
   Bool_t          truth_isBckg_withPrimaryPi0;
   Bool_t          truth_isBckg_withSecondaryPi0;
   Int_t           truth_N_FSParticles;
   Int_t           truth_N_gamma;
   Int_t           truth_N_pi0;
   Int_t           truth_N_proton;
   Int_t           truth_N_trueMichelElectrons;
   Int_t           truth_muon_charge;
   Int_t           truth_pi0_GrandMother;
   Int_t           truth_pi0_GrandMotherStatus;
   Int_t           truth_pi0_Mother;
   Int_t           truth_pi0_MotherStatus;
   Int_t           truth_pi0_status;
   Int_t           truth_target_material;
   Int_t           truth_vertex_module;
   Int_t           truth_vertex_plane;
   Double_t        truth_eventID;
   Double_t        truth_michelElectron_E;
   Double_t        truth_michelElectron_P;
   Double_t        truth_michelMuon_P;
   Double_t        truth_michelMuon_end_dist_vtx;
   Double_t        truth_michelMuon_length;
   Double_t        truth_michelPion_P;
   Double_t        truth_michelPion_begin_dist_vtx;
   Double_t        truth_michelPion_length;
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
   Double_t        truth_michelMuon_endPoint[3];
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
   Double_t        CCProtonPi0_total_E;
   Double_t        CCProtonPi0_total_px;
   Double_t        CCProtonPi0_total_py;
   Double_t        CCProtonPi0_total_pz;
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
   Int_t           CCProtonPi0_trajProtonProngPDG[10];
   Int_t           CCProtonPi0_trajProtonProngPrimary[10];
   Double_t        CCProtonPi0_endProtonTrajMomentum[10];
   Double_t        CCProtonPi0_endProtonTrajXPosition[10];
   Double_t        CCProtonPi0_endProtonTrajYPosition[10];
   Double_t        CCProtonPi0_endProtonTrajZPosition[10];
   Double_t        CCProtonPi0_pionScore[10];
   Double_t        CCProtonPi0_protonScore[10];
   Double_t        CCProtonPi0_protonScore_LLR[10];
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
   Double_t        mc_FSPartPx[23];   //[mc_nFSPart]
   Double_t        mc_FSPartPy[23];   //[mc_nFSPart]
   Double_t        mc_FSPartPz[23];   //[mc_nFSPart]
   Double_t        mc_FSPartE[23];   //[mc_nFSPart]
   Int_t           mc_FSPartPDG[23];   //[mc_nFSPart]
   Int_t           mc_er_nPart;
   Int_t           mc_er_ID[52];   //[mc_er_nPart]
   Int_t           mc_er_status[52];   //[mc_er_nPart]
   Double_t        mc_er_posInNucX[52];   //[mc_er_nPart]
   Double_t        mc_er_posInNucY[52];   //[mc_er_nPart]
   Double_t        mc_er_posInNucZ[52];   //[mc_er_nPart]
   Double_t        mc_er_Px[52];   //[mc_er_nPart]
   Double_t        mc_er_Py[52];   //[mc_er_nPart]
   Double_t        mc_er_Pz[52];   //[mc_er_nPart]
   Double_t        mc_er_E[52];   //[mc_er_nPart]
   Int_t           mc_er_FD[52];   //[mc_er_nPart]
   Int_t           mc_er_LD[52];   //[mc_er_nPart]
   Int_t           mc_er_mother[52];   //[mc_er_nPart]
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
   TBranch        *b_Cut_VtxBlob;   //!
   TBranch        *b_Cut_nProngs;   //!
   TBranch        *b_Cut_secEndPoint_Michel_Exist;   //!
   TBranch        *b_anglescan_ncand;   //!
   TBranch        *b_anglescan_ncandx;   //!
   TBranch        *b_blob_ndof_1;   //!
   TBranch        *b_blob_ndof_2;   //!
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
   TBranch        *b_n_vtx_michel_views;   //!
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
   TBranch        *b_blob_fval_1;   //!
   TBranch        *b_blob_fval_2;   //!
   TBranch        *b_endpoint_michel_distance;   //!
   TBranch        *b_energyUnused_afterReco;   //!
   TBranch        *b_energyUsed_afterReco;   //!
   TBranch        *b_energy_from_mc;   //!
   TBranch        *b_energy_from_mc_fraction;   //!
   TBranch        *b_energy_from_mc_fraction_of_highest;   //!
   TBranch        *b_evis_ECAL;   //!
   TBranch        *b_evis_HCAL;   //!
   TBranch        *b_evis_NuclearTarget;   //!
   TBranch        *b_evis_TotalExceptNuclearTarget;   //!
   TBranch        *b_evis_Tracker;   //!
   TBranch        *b_evis_nearvtx;   //!
   TBranch        *b_evis_total;   //!
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
   TBranch        *b_gamma1_dist_vtx;   //!
   TBranch        *b_gamma1_evis_ecal;   //!
   TBranch        *b_gamma1_evis_hcal;   //!
   TBranch        *b_gamma1_evis_scal;   //!
   TBranch        *b_gamma1_evis_trkr;   //!
   TBranch        *b_gamma1_phi;   //!
   TBranch        *b_gamma1_px;   //!
   TBranch        *b_gamma1_py;   //!
   TBranch        *b_gamma1_pz;   //!
   TBranch        *b_gamma1_score;   //!
   TBranch        *b_gamma1_theta;   //!
   TBranch        *b_gamma1_time;   //!
   TBranch        *b_gamma2_E;   //!
   TBranch        *b_gamma2_dEdx;   //!
   TBranch        *b_gamma2_dist_vtx;   //!
   TBranch        *b_gamma2_evis_ecal;   //!
   TBranch        *b_gamma2_evis_hcal;   //!
   TBranch        *b_gamma2_evis_scal;   //!
   TBranch        *b_gamma2_evis_trkr;   //!
   TBranch        *b_gamma2_phi;   //!
   TBranch        *b_gamma2_px;   //!
   TBranch        *b_gamma2_py;   //!
   TBranch        *b_gamma2_pz;   //!
   TBranch        *b_gamma2_score;   //!
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
   TBranch        *b_prong_showerScore;   //!
   TBranch        *b_reco_eventID;   //!
   TBranch        *b_time;   //!
   TBranch        *b_totalIDVisibleE;   //!
   TBranch        *b_totalODVisibleE;   //!
   TBranch        *b_totalVisibleE;   //!
   TBranch        *b_vtx_michel_distance;   //!
   TBranch        *b_g1dedx_cluster_occupancy_sz;   //!
   TBranch        *b_g1dedx_cluster_occupancy;   //!
   TBranch        *b_g2dedx_cluster_occupancy_sz;   //!
   TBranch        *b_g2dedx_cluster_occupancy;   //!
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
   TBranch        *b_truth_isSignal_1Pi0;   //!
   TBranch        *b_truth_isSignal_2Gamma;   //!
   TBranch        *b_truth_isFidVol;   //!
   TBranch        *b_truth_AnalyzeEvent;   //!
   TBranch        *b_truth_isBckg_QELike;   //!
   TBranch        *b_truth_isBckg_SinglePiPlus;   //!
   TBranch        *b_truth_isBckg_SinglePiMinus;   //!
   TBranch        *b_truth_isBckg_MultiPion;   //!
   TBranch        *b_truth_isBckg_MultiPiZero;   //!
   TBranch        *b_truth_isBckg_Other;   //!
   TBranch        *b_truth_isBckg_withAntiMuon;   //!
   TBranch        *b_truth_isBckg_withMichel;   //!
   TBranch        *b_truth_isBckg_withPrimaryPi0;   //!
   TBranch        *b_truth_isBckg_withSecondaryPi0;   //!
   TBranch        *b_truth_N_FSParticles;   //!
   TBranch        *b_truth_N_gamma;   //!
   TBranch        *b_truth_N_pi0;   //!
   TBranch        *b_truth_N_proton;   //!
   TBranch        *b_truth_N_trueMichelElectrons;   //!
   TBranch        *b_truth_muon_charge;   //!
   TBranch        *b_truth_pi0_GrandMother;   //!
   TBranch        *b_truth_pi0_GrandMotherStatus;   //!
   TBranch        *b_truth_pi0_Mother;   //!
   TBranch        *b_truth_pi0_MotherStatus;   //!
   TBranch        *b_truth_pi0_status;   //!
   TBranch        *b_truth_target_material;   //!
   TBranch        *b_truth_vertex_module;   //!
   TBranch        *b_truth_vertex_plane;   //!
   TBranch        *b_truth_eventID;   //!
   TBranch        *b_truth_michelElectron_E;   //!
   TBranch        *b_truth_michelElectron_P;   //!
   TBranch        *b_truth_michelMuon_P;   //!
   TBranch        *b_truth_michelMuon_end_dist_vtx;   //!
   TBranch        *b_truth_michelMuon_length;   //!
   TBranch        *b_truth_michelPion_P;   //!
   TBranch        *b_truth_michelPion_begin_dist_vtx;   //!
   TBranch        *b_truth_michelPion_length;   //!
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
   TBranch        *b_truth_michelMuon_endPoint;   //!
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
   TBranch        *b_CCProtonPi0_total_E;   //!
   TBranch        *b_CCProtonPi0_total_px;   //!
   TBranch        *b_CCProtonPi0_total_py;   //!
   TBranch        *b_CCProtonPi0_total_pz;   //!
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
   TBranch        *b_CCProtonPi0_trajProtonProngPDG;   //!
   TBranch        *b_CCProtonPi0_trajProtonProngPrimary;   //!
   TBranch        *b_CCProtonPi0_endProtonTrajMomentum;   //!
   TBranch        *b_CCProtonPi0_endProtonTrajXPosition;   //!
   TBranch        *b_CCProtonPi0_endProtonTrajYPosition;   //!
   TBranch        *b_CCProtonPi0_endProtonTrajZPosition;   //!
   TBranch        *b_CCProtonPi0_pionScore;   //!
   TBranch        *b_CCProtonPi0_protonScore;   //!
   TBranch        *b_CCProtonPi0_protonScore_LLR;   //!
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


};

#endif
