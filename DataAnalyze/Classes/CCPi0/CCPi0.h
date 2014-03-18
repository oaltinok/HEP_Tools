/*
================================================================================
Class: CCPi0
    Core Class for CCPi0 Analysis Package
    Includes all data variables for Data Analysis
    Member functions uses member data variables
    
    Outputs .ROOT file filled with 1D or 2D histograms, which is an input file
    for Plotter Class
    
    Uses other classes to define variables or to access specific functions
    
    Main Directory:
        Classes/CCPi0
        
    Usage:
        > main.cpp declares and controls the class
        > See run function Comments
    
    
    
    Last Revision: 2014_03_18
================================================================================
*/

#ifndef CCPi0_h
#define CCPi0_h

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
#include "Classes/Muon/Muon.cpp"
#include "Classes/Proton/Proton.cpp"
#include "Classes/Pion/Pion.cpp"

const double mevSq_to_gevSq = pow(10,6);

class CCPi0 {
public :
   // -------------------------------------------------------------------------
   //     Specific Functions
   //--------------------------------------------------------------------------
   
    
    
   // -------------------------------------------------------------------------
   //     void run(): Generates a .root file with selected histograms
   //         playlist -> address of the playlist
   //         filename -> file name for the output .root file
   //--------------------------------------------------------------------------
    void run(string playlist);
    
    
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
    
    void fillParticleTrue(Particle& part);
    
    void fillHistograms();

    void openFiles();
    void closeFiles();
    void writeReadme();
    void writeCutFile();

   // -------------------------------------------------------------------------
   //     Default Functions
   //--------------------------------------------------------------------------
    CCPi0();
    ~CCPi0();

    void Init(string playlist, TChain* fChain);
    Int_t GetEntry(Long64_t entry);
    Long64_t LoadTree(Long64_t entry);
    Bool_t Notify();
    void Show(Long64_t entry);
    Int_t Cut(Long64_t entry);
    
   //-------------------------------------------------------------------------
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
    
    TH1F* int_channel;
    TH1F* vertex_z;
    TH1F* n_FSParticles;
    TH1F* n_gammas;
    
    
   // -------------------------------------------------------------------------
   //     Analysis Variables
   //--------------------------------------------------------------------------
    bool isDataAnalysis;
    bool isMC;
    double maxBeamEnergy;
    int max_nFSPart;
    HEP_Misc misc;
    TVector3 beam_p3;
    Proton proton;
    Muon muon;
    Pion pion;
    
   // -------------------------------------------------------------------------
   //     Cut Numbers and List of Bins
   //--------------------------------------------------------------------------
    BinList binList;
    CutNumberList* nCutList;
    
   // -------------------------------------------------------------------------
   //     Files
   //--------------------------------------------------------------------------
    string rootDir;
    string plotDir;
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
   Bool_t          is_GoodDirection1;
   Bool_t          is_GoodPosition1;
   Bool_t          is_GoodDirection2;
   Bool_t          is_GoodPosition2;
   Bool_t          is_GoodBlob1;
   Bool_t          is_GoodBlob2;
   Bool_t          is_anglescan;
   Bool_t          is_anglescan_applied;
   Bool_t          is_houghtransform;
   Bool_t          is_houghtransform_applied;
   Bool_t          is_twoDBlob;
   Int_t           anglescan_ncand;
   Int_t           anglescan_ncandx;
   Int_t           blob_ndof_1;
   Int_t           blob_ndof_2;
   Int_t           broken_track_most_us_plane;
   Int_t           discard_track_count;
   Int_t           dmode;
   Int_t           g1blob_ncluster;
   Int_t           g1blob_ndigit;
   Int_t           g1convidet;
   Int_t           g1mostevispdg;
   Int_t           g2blob_ncluster;
   Int_t           g2blob_ndigit;
   Int_t           g2convidet;
   Int_t           g2mostevispdg;
   Int_t           minos_trk_end_plane;
   Int_t           minos_trk_is_contained;
   Int_t           minos_trk_is_ok;
   Int_t           minos_trk_quality;
   Int_t           minos_trk_used_curvature;
   Int_t           minos_trk_used_range;
   Int_t           nblob_anglescan;
   Int_t           nblob_hough;
   Int_t           nmeson;
   Int_t           nmumcapture;
   Int_t           nmumdecay;
   Int_t           nmupdecay;
   Int_t           nn;
   Int_t           np;
   Int_t           npi0;
   Int_t           npi02;
   Int_t           npim;
   Int_t           npim2;
   Int_t           npimcapture;
   Int_t           npimdecay;
   Int_t           npiminelastic;
   Int_t           npip;
   Int_t           npip2;
   Int_t           npipcapture;
   Int_t           npipdecay;
   Int_t           npipinelastic;
   Int_t           npipm;
   Int_t           od_energeticTower;
   Int_t           phys_energy_in_road_downstream_nplanes;
   Int_t           phys_energy_in_road_upstream_nplanes;
   Int_t           phys_n_dead_discr_pair;
   Int_t           phys_n_dead_discr_pair_in_prim_track_region;
   Int_t           phys_n_dead_discr_pair_two_mod_downstream_prim_track;
   Int_t           phys_n_dead_discr_pair_two_mod_upstream_prim_vtx;
   Int_t           phys_n_dead_discr_pair_upstream_prim_track_proj;
   Int_t           phys_vertex_is_fiducial;
   Int_t           primary_index;
   Int_t           primary_multiplicity;
   Int_t           primary_multiplicity2;
   Int_t           survive_all;
   Int_t           survive_do_muon;
   Int_t           survive_do_vertex;
   Int_t           survive_fiducial;
   Int_t           survive_gammatrack;
   Int_t           survive_has_vertex;
   Int_t           survive_minos_match;
   Int_t           survive_onetrackpervtx;
   Int_t           survive_plausible;
   Int_t           survive_prefilter;
   Int_t           survive_three_vertex;
   Int_t           survive_vtx_blob;
   Int_t           tfiducial;
   Int_t           vertex_count;
   Int_t           vertex_count2;
   Double_t        Dispersed_blob_energy;
   Double_t        Erec;
   Double_t        Filament_Vertex_energy;
   Double_t        Muon_blob_energy;
   Double_t        Q2;
   Double_t        RE_energy_ECAL;
   Double_t        RE_energy_HCAL;
   Double_t        RE_energy_Tracker;
   Double_t        RE_photon_dEdx_1;
   Double_t        RE_photon_dEdx_2;
   Double_t        RE_photon_energy_1;
   Double_t        RE_photon_energy_2;
   Double_t        RE_photon_time_1;
   Double_t        RE_photon_time_2;
   Double_t        RE_scalar;
   Double_t        Rejected_blob_vis_energy;
   Double_t        Sphere_Vertex_energy;
   Double_t        Tn;
   Double_t        Vertex_blob_energy;
   Double_t        W;
   Double_t        W2;
   Double_t        blob_fval_1;
   Double_t        blob_fval_2;
   Double_t        ecalevis;
   Double_t        energy_from_mc;
   Double_t        energy_from_mc_fraction;
   Double_t        energy_from_mc_fraction_of_highest;
   Double_t        g1blob_edge_distance;
   Double_t        g1blob_minsep;
   Double_t        g1blob_vtx_distance;
   Double_t        g1convdist;
   Double_t        g1e;
   Double_t        g1e0;
   Double_t        g1ecaledep;
   Double_t        g1ecalo;
   Double_t        g1g1evis;
   Double_t        g1g2evis;
   Double_t        g1gmevis;
   Double_t        g1hcaledep;
   Double_t        g1idetedep;
   Double_t        g1mostevisfrac;
   Double_t        g1muevis;
   Double_t        g1neutronevis;
   Double_t        g1nukeedep;
   Double_t        g1odetedep;
   Double_t        g1otherevis;
   Double_t        g1othersubdetedep;
   Double_t        g1phi;
   Double_t        g1phi0;
   Double_t        g1pi0evis;
   Double_t        g1pimevis;
   Double_t        g1pipevis;
   Double_t        g1protonevis;
   Double_t        g1sharedevis;
   Double_t        g1sideedep;
   Double_t        g1theta;
   Double_t        g1theta0;
   Double_t        g1totalevis;
   Double_t        g1trkredep;
   Double_t        g2blob_edge_distance;
   Double_t        g2blob_minsep;
   Double_t        g2blob_vtx_distance;
   Double_t        g2convdist;
   Double_t        g2e;
   Double_t        g2e0;
   Double_t        g2ecaledep;
   Double_t        g2ecalo;
   Double_t        g2g1evis;
   Double_t        g2g2evis;
   Double_t        g2gmevis;
   Double_t        g2hcaledep;
   Double_t        g2idetedep;
   Double_t        g2mostevisfrac;
   Double_t        g2muevis;
   Double_t        g2neutronevis;
   Double_t        g2nukeedep;
   Double_t        g2odetedep;
   Double_t        g2otherevis;
   Double_t        g2othersubdetedep;
   Double_t        g2phi;
   Double_t        g2phi0;
   Double_t        g2pi0evis;
   Double_t        g2pimevis;
   Double_t        g2pipevis;
   Double_t        g2protonevis;
   Double_t        g2sharedevis;
   Double_t        g2sideedep;
   Double_t        g2theta;
   Double_t        g2theta0;
   Double_t        g2totalevis;
   Double_t        g2trkredep;
   Double_t        hcalevis;
   Double_t        mgg;
   Double_t        minos_trk_eqp;
   Double_t        minos_trk_fit_pass;
   Double_t        minos_trk_p;
   Double_t        minos_trk_p_curvature;
   Double_t        minos_trk_p_range;
   Double_t        minos_trk_qp;
   Double_t        muon_phi;
   Double_t        muon_theta;
   Double_t        muon_thetaX;
   Double_t        muon_thetaY;
   Double_t        neutronecaledep;
   Double_t        neutronecalo;
   Double_t        neutronhcaledep;
   Double_t        neutronidetedep;
   Double_t        neutronnukeedep;
   Double_t        neutronodetedep;
   Double_t        neutronothersubdetedep;
   Double_t        neutronsideedep;
   Double_t        neutrontrkredep;
   Double_t        nke;
   Double_t        ntgtevis;
   Double_t        oangle;
   Double_t        oangle0;
   Double_t        oangle0x;
   Double_t        od_downstreamFrame;
   Double_t        od_downstreamFrame_z;
   Double_t        od_highStory;
   Double_t        od_highStory_t;
   Double_t        od_lowStory;
   Double_t        od_lowStory_t;
   Double_t        od_maxEnergy;
   Double_t        od_upstreamFrame;
   Double_t        od_upstreamFrame_z;
   Double_t        otherevis;
   Double_t        phys_energy_dispersed;
   Double_t        phys_energy_in_road_downstream;
   Double_t        phys_energy_in_road_upstream;
   Double_t        phys_energy_unattached;
   Double_t        pi0_evis_dispersed_blob;
   Double_t        pi0_evis_muon_blob;
   Double_t        pi0_evis_outtime_blob;
   Double_t        pi0_evis_vtx_blob;
   Double_t        pi0_evisfrac_dispersed_blob;
   Double_t        pi0_evisfrac_muon_blob;
   Double_t        pi0_evisfrac_outtime_blob;
   Double_t        pi0_evisfrac_vtx_blob;
   Double_t        pi0ecaledep;
   Double_t        pi0ecalo;
   Double_t        pi0hcaledep;
   Double_t        pi0idetedep;
   Double_t        pi0nukeedep;
   Double_t        pi0odetedep;
   Double_t        pi0othersubdetedep;
   Double_t        pi0sideedep;
   Double_t        pi0trkredep;
   Double_t        pienergy;
   Double_t        pienergy0;
   Double_t        pimlength;
   Double_t        piphi;
   Double_t        piphi0;
   Double_t        piplength;
   Double_t        pitheta;
   Double_t        pitheta0;
   Double_t        pke;
   Double_t        prim_vtx_smallest_opening_angle;
   Double_t        totalevis;
   Double_t        trkrevis;
   Int_t           anglescan_blob_nc_sz;
   Int_t           anglescan_blob_nc[10];   //[anglescan_blob_nc_sz]
   Int_t           anglescan_blob_ncu_sz;
   Int_t           anglescan_blob_ncu[10];   //[anglescan_blob_ncu_sz]
   Int_t           anglescan_blob_ncv_sz;
   Int_t           anglescan_blob_ncv[10];   //[anglescan_blob_ncv_sz]
   Int_t           anglescan_blob_ncx_sz;
   Int_t           anglescan_blob_ncx[10];   //[anglescan_blob_ncx_sz]
   Int_t           anglescan_blob_nd_sz;
   Int_t           anglescan_blob_nd[10];   //[anglescan_blob_nd_sz]
   Int_t           anglescan_blob_ndu_sz;
   Int_t           anglescan_blob_ndu[10];   //[anglescan_blob_ndu_sz]
   Int_t           anglescan_blob_ndv_sz;
   Int_t           anglescan_blob_ndv[10];   //[anglescan_blob_ndv_sz]
   Int_t           anglescan_blob_ndx_sz;
   Int_t           anglescan_blob_ndx[10];   //[anglescan_blob_ndx_sz]
   Int_t           anglescan_cand_nc_sz;
   Int_t           anglescan_cand_nc[10];   //[anglescan_cand_nc_sz]
   Int_t           anglescan_cand_ncu_sz;
   Int_t           anglescan_cand_ncu[10];   //[anglescan_cand_ncu_sz]
   Int_t           anglescan_cand_ncv_sz;
   Int_t           anglescan_cand_ncv[10];   //[anglescan_cand_ncv_sz]
   Int_t           anglescan_cand_ncx_sz;
   Int_t           anglescan_cand_ncx[10];   //[anglescan_cand_ncx_sz]
   Int_t           anglescan_cand_nd_sz;
   Int_t           anglescan_cand_nd[10];   //[anglescan_cand_nd_sz]
   Int_t           anglescan_cand_ndu_sz;
   Int_t           anglescan_cand_ndu[10];   //[anglescan_cand_ndu_sz]
   Int_t           anglescan_cand_ndv_sz;
   Int_t           anglescan_cand_ndv[10];   //[anglescan_cand_ndv_sz]
   Int_t           anglescan_cand_ndx_sz;
   Int_t           anglescan_cand_ndx[10];   //[anglescan_cand_ndx_sz]
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
   Int_t           multiplicities_sz;
   Int_t           multiplicities[5];   //[multiplicities_sz]
   Int_t           primary_truth_counts_sz;
   Int_t           primary_truth_counts[5];   //[primary_truth_counts_sz]
   Int_t           primary_truth_pdgs1_sz;
   Int_t           primary_truth_pdgs1[5];   //[primary_truth_pdgs1_sz]
   Int_t           primary_truth_pdgs2_sz;
   Int_t           primary_truth_pdgs2[5];   //[primary_truth_pdgs2_sz]
   Int_t           primary_truth_pdgs3_sz;
   Int_t           primary_truth_pdgs3[5];   //[primary_truth_pdgs3_sz]
   Int_t           RE_photon_direction_1_sz;
   Double_t        RE_photon_direction_1[3];   //[RE_photon_direction_1_sz]
   Int_t           RE_photon_direction_2_sz;
   Double_t        RE_photon_direction_2[3];   //[RE_photon_direction_2_sz]
   Int_t           RE_photon_vertex_1_sz;
   Double_t        RE_photon_vertex_1[3];   //[RE_photon_vertex_1_sz]
   Int_t           RE_photon_vertex_2_sz;
   Double_t        RE_photon_vertex_2[3];   //[RE_photon_vertex_2_sz]
   Int_t           deviations_sz;
   Double_t        deviations[5];   //[deviations_sz]
   Int_t           g1convpos_sz;
   Double_t        g1convpos[3];   //[g1convpos_sz]
   Int_t           g1mom_sz;
   Double_t        g1mom[3];   //[g1mom_sz]
   Int_t           g1mom0_sz;
   Double_t        g1mom0[4];   //[g1mom0_sz]
   Int_t           g2convpos_sz;
   Double_t        g2convpos[3];   //[g2convpos_sz]
   Int_t           g2mom_sz;
   Double_t        g2mom[3];   //[g2mom_sz]
   Int_t           g2mom0_sz;
   Double_t        g2mom0[4];   //[g2mom0_sz]
   Int_t           michel_mom_sz;
   Double_t        michel_mom[4];   //[michel_mom_sz]
   Int_t           michel_pos_sz;
   Double_t        michel_pos[4];   //[michel_pos_sz]
   Int_t           mumom_sz;
   Double_t        mumom[5];   //[mumom_sz]
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
   Int_t           pimom_sz;
   Double_t        pimom[4];   //[pimom_sz]
   Int_t           pimom0_sz;
   Double_t        pimom0[4];   //[pimom0_sz]
   Int_t           primary_separations_sz;
   Double_t        primary_separations[5];   //[primary_separations_sz]
   Int_t           primary_trklengths_sz;
   Double_t        primary_trklengths[5];   //[primary_trklengths_sz]
   Int_t           primary_truth_fractions1_sz;
   Double_t        primary_truth_fractions1[5];   //[primary_truth_fractions1_sz]
   Int_t           primary_truth_fractions2_sz;
   Double_t        primary_truth_fractions2[5];   //[primary_truth_fractions2_sz]
   Int_t           primary_truth_fractions3_sz;
   Double_t        primary_truth_fractions3[5];   //[primary_truth_fractions3_sz]
   Int_t           primary_truth_shareds_sz;
   Double_t        primary_truth_shareds[5];   //[primary_truth_shareds_sz]
   Bool_t          truth_has_physics_event;
   Bool_t          truth_reco_minos_match;
   Bool_t          truth_is_fiducial;
   Bool_t          truth_pass_plausible;
   Bool_t          truth_is_ccpi0;
   Bool_t          truth_is_cc1pi0;
   Bool_t          truth_is_ccpi0secondary;
   Bool_t          truth_is_by_pim;
   Bool_t          truth_is_ccpi0x;
   Bool_t          truth_is_other;
   Double_t        truth_MC_photon_energy_1;
   Double_t        truth_MC_photon_energy_2;
   Double_t        truth_MC_pi0_energy;
   Double_t        truth_MC_scalar;
   Double_t        truth_fslepton_E;
   Double_t        truth_fslepton_P;
   Double_t        truth_fslepton_T;
   Double_t        truth_fslepton_phi;
   Double_t        truth_fslepton_theta;
   Double_t        truth_fslepton_theta_x;
   Double_t        truth_fslepton_theta_y;
   Int_t           truth_MC_photon_direction_1_sz;
   Double_t        truth_MC_photon_direction_1[3];   //[truth_MC_photon_direction_1_sz]
   Int_t           truth_MC_photon_direction_2_sz;
   Double_t        truth_MC_photon_direction_2[3];   //[truth_MC_photon_direction_2_sz]
   Int_t           truth_MC_photon_vertex_1_sz;
   Double_t        truth_MC_photon_vertex_1[3];   //[truth_MC_photon_vertex_1_sz]
   Int_t           truth_MC_photon_vertex_2_sz;
   Double_t        truth_MC_photon_vertex_2[3];   //[truth_MC_photon_vertex_2_sz]
   Int_t           truth_MC_pi0_momentum_sz;
   Double_t        truth_MC_pi0_momentum[3];   //[truth_MC_pi0_momentum_sz]
   Int_t           CCPi0AnaTool_nuFlavor;
   Int_t           CCPi0AnaTool_nuHelicity;
   Int_t           CCPi0AnaTool_intCurrent;
   Int_t           CCPi0AnaTool_intType;
   Double_t        CCPi0AnaTool_E;
   Double_t        CCPi0AnaTool_Q2;
   Double_t        CCPi0AnaTool_x;
   Double_t        CCPi0AnaTool_y;
   Double_t        CCPi0AnaTool_W;
   Double_t        CCPi0AnaTool_score;
   Double_t        CCPi0AnaTool_leptonE[4];
   Double_t        CCPi0AnaTool_vtx[4];
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
   Double_t        mc_wgt_GENIE[1000];   //[mc_wgt_GENIE_sz]
   Int_t           mc_wgt_Flux_Tertiary_sz;
   Double_t        mc_wgt_Flux_Tertiary[1000];   //[mc_wgt_Flux_Tertiary_sz]
   Int_t           mc_wgt_Flux_BeamFocus_sz;
   Double_t        mc_wgt_Flux_BeamFocus[1000];   //[mc_wgt_Flux_BeamFocus_sz]
   Int_t           mc_wgt_Flux_NA49_sz;
   Double_t        mc_wgt_Flux_NA49[1000];   //[mc_wgt_Flux_NA49_sz]
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
   TBranch        *b_is_GoodDirection1;   //!
   TBranch        *b_is_GoodPosition1;   //!
   TBranch        *b_is_GoodDirection2;   //!
   TBranch        *b_is_GoodPosition2;   //!
   TBranch        *b_is_GoodBlob1;   //!
   TBranch        *b_is_GoodBlob2;   //!
   TBranch        *b_is_anglescan;   //!
   TBranch        *b_is_anglescan_applied;   //!
   TBranch        *b_is_houghtransform;   //!
   TBranch        *b_is_houghtransform_applied;   //!
   TBranch        *b_is_twoDBlob;   //!
   TBranch        *b_anglescan_ncand;   //!
   TBranch        *b_anglescan_ncandx;   //!
   TBranch        *b_blob_ndof_1;   //!
   TBranch        *b_blob_ndof_2;   //!
   TBranch        *b_broken_track_most_us_plane;   //!
   TBranch        *b_discard_track_count;   //!
   TBranch        *b_dmode;   //!
   TBranch        *b_g1blob_ncluster;   //!
   TBranch        *b_g1blob_ndigit;   //!
   TBranch        *b_g1convidet;   //!
   TBranch        *b_g1mostevispdg;   //!
   TBranch        *b_g2blob_ncluster;   //!
   TBranch        *b_g2blob_ndigit;   //!
   TBranch        *b_g2convidet;   //!
   TBranch        *b_g2mostevispdg;   //!
   TBranch        *b_minos_trk_end_plane;   //!
   TBranch        *b_minos_trk_is_contained;   //!
   TBranch        *b_minos_trk_is_ok;   //!
   TBranch        *b_minos_trk_quality;   //!
   TBranch        *b_minos_trk_used_curvature;   //!
   TBranch        *b_minos_trk_used_range;   //!
   TBranch        *b_nblob_anglescan;   //!
   TBranch        *b_nblob_hough;   //!
   TBranch        *b_nmeson;   //!
   TBranch        *b_nmumcapture;   //!
   TBranch        *b_nmumdecay;   //!
   TBranch        *b_nmupdecay;   //!
   TBranch        *b_nn;   //!
   TBranch        *b_np;   //!
   TBranch        *b_npi0;   //!
   TBranch        *b_npi02;   //!
   TBranch        *b_npim;   //!
   TBranch        *b_npim2;   //!
   TBranch        *b_npimcapture;   //!
   TBranch        *b_npimdecay;   //!
   TBranch        *b_npiminelastic;   //!
   TBranch        *b_npip;   //!
   TBranch        *b_npip2;   //!
   TBranch        *b_npipcapture;   //!
   TBranch        *b_npipdecay;   //!
   TBranch        *b_npipinelastic;   //!
   TBranch        *b_npipm;   //!
   TBranch        *b_od_energeticTower;   //!
   TBranch        *b_phys_energy_in_road_downstream_nplanes;   //!
   TBranch        *b_phys_energy_in_road_upstream_nplanes;   //!
   TBranch        *b_phys_n_dead_discr_pair;   //!
   TBranch        *b_phys_n_dead_discr_pair_in_prim_track_region;   //!
   TBranch        *b_phys_n_dead_discr_pair_two_mod_downstream_prim_track;   //!
   TBranch        *b_phys_n_dead_discr_pair_two_mod_upstream_prim_vtx;   //!
   TBranch        *b_phys_n_dead_discr_pair_upstream_prim_track_proj;   //!
   TBranch        *b_phys_vertex_is_fiducial;   //!
   TBranch        *b_primary_index;   //!
   TBranch        *b_primary_multiplicity;   //!
   TBranch        *b_primary_multiplicity2;   //!
   TBranch        *b_survive_all;   //!
   TBranch        *b_survive_do_muon;   //!
   TBranch        *b_survive_do_vertex;   //!
   TBranch        *b_survive_fiducial;   //!
   TBranch        *b_survive_gammatrack;   //!
   TBranch        *b_survive_has_vertex;   //!
   TBranch        *b_survive_minos_match;   //!
   TBranch        *b_survive_onetrackpervtx;   //!
   TBranch        *b_survive_plausible;   //!
   TBranch        *b_survive_prefilter;   //!
   TBranch        *b_survive_three_vertex;   //!
   TBranch        *b_survive_vtx_blob;   //!
   TBranch        *b_tfiducial;   //!
   TBranch        *b_vertex_count;   //!
   TBranch        *b_vertex_count2;   //!
   TBranch        *b_Dispersed_blob_energy;   //!
   TBranch        *b_Erec;   //!
   TBranch        *b_Filament_Vertex_energy;   //!
   TBranch        *b_Muon_blob_energy;   //!
   TBranch        *b_Q2;   //!
   TBranch        *b_RE_energy_ECAL;   //!
   TBranch        *b_RE_energy_HCAL;   //!
   TBranch        *b_RE_energy_Tracker;   //!
   TBranch        *b_RE_photon_dEdx_1;   //!
   TBranch        *b_RE_photon_dEdx_2;   //!
   TBranch        *b_RE_photon_energy_1;   //!
   TBranch        *b_RE_photon_energy_2;   //!
   TBranch        *b_RE_photon_time_1;   //!
   TBranch        *b_RE_photon_time_2;   //!
   TBranch        *b_RE_scalar;   //!
   TBranch        *b_Rejected_blob_vis_energy;   //!
   TBranch        *b_Sphere_Vertex_energy;   //!
   TBranch        *b_Tn;   //!
   TBranch        *b_Vertex_blob_energy;   //!
   TBranch        *b_W;   //!
   TBranch        *b_W2;   //!
   TBranch        *b_blob_fval_1;   //!
   TBranch        *b_blob_fval_2;   //!
   TBranch        *b_ecalevis;   //!
   TBranch        *b_energy_from_mc;   //!
   TBranch        *b_energy_from_mc_fraction;   //!
   TBranch        *b_energy_from_mc_fraction_of_highest;   //!
   TBranch        *b_g1blob_edge_distance;   //!
   TBranch        *b_g1blob_minsep;   //!
   TBranch        *b_g1blob_vtx_distance;   //!
   TBranch        *b_g1convdist;   //!
   TBranch        *b_g1e;   //!
   TBranch        *b_g1e0;   //!
   TBranch        *b_g1ecaledep;   //!
   TBranch        *b_g1ecalo;   //!
   TBranch        *b_g1g1evis;   //!
   TBranch        *b_g1g2evis;   //!
   TBranch        *b_g1gmevis;   //!
   TBranch        *b_g1hcaledep;   //!
   TBranch        *b_g1idetedep;   //!
   TBranch        *b_g1mostevisfrac;   //!
   TBranch        *b_g1muevis;   //!
   TBranch        *b_g1neutronevis;   //!
   TBranch        *b_g1nukeedep;   //!
   TBranch        *b_g1odetedep;   //!
   TBranch        *b_g1otherevis;   //!
   TBranch        *b_g1othersubdetedep;   //!
   TBranch        *b_g1phi;   //!
   TBranch        *b_g1phi0;   //!
   TBranch        *b_g1pi0evis;   //!
   TBranch        *b_g1pimevis;   //!
   TBranch        *b_g1pipevis;   //!
   TBranch        *b_g1protonevis;   //!
   TBranch        *b_g1sharedevis;   //!
   TBranch        *b_g1sideedep;   //!
   TBranch        *b_g1theta;   //!
   TBranch        *b_g1theta0;   //!
   TBranch        *b_g1totalevis;   //!
   TBranch        *b_g1trkredep;   //!
   TBranch        *b_g2blob_edge_distance;   //!
   TBranch        *b_g2blob_minsep;   //!
   TBranch        *b_g2blob_vtx_distance;   //!
   TBranch        *b_g2convdist;   //!
   TBranch        *b_g2e;   //!
   TBranch        *b_g2e0;   //!
   TBranch        *b_g2ecaledep;   //!
   TBranch        *b_g2ecalo;   //!
   TBranch        *b_g2g1evis;   //!
   TBranch        *b_g2g2evis;   //!
   TBranch        *b_g2gmevis;   //!
   TBranch        *b_g2hcaledep;   //!
   TBranch        *b_g2idetedep;   //!
   TBranch        *b_g2mostevisfrac;   //!
   TBranch        *b_g2muevis;   //!
   TBranch        *b_g2neutronevis;   //!
   TBranch        *b_g2nukeedep;   //!
   TBranch        *b_g2odetedep;   //!
   TBranch        *b_g2otherevis;   //!
   TBranch        *b_g2othersubdetedep;   //!
   TBranch        *b_g2phi;   //!
   TBranch        *b_g2phi0;   //!
   TBranch        *b_g2pi0evis;   //!
   TBranch        *b_g2pimevis;   //!
   TBranch        *b_g2pipevis;   //!
   TBranch        *b_g2protonevis;   //!
   TBranch        *b_g2sharedevis;   //!
   TBranch        *b_g2sideedep;   //!
   TBranch        *b_g2theta;   //!
   TBranch        *b_g2theta0;   //!
   TBranch        *b_g2totalevis;   //!
   TBranch        *b_g2trkredep;   //!
   TBranch        *b_hcalevis;   //!
   TBranch        *b_mgg;   //!
   TBranch        *b_minos_trk_eqp;   //!
   TBranch        *b_minos_trk_fit_pass;   //!
   TBranch        *b_minos_trk_p;   //!
   TBranch        *b_minos_trk_p_curvature;   //!
   TBranch        *b_minos_trk_p_range;   //!
   TBranch        *b_minos_trk_qp;   //!
   TBranch        *b_muon_phi;   //!
   TBranch        *b_muon_theta;   //!
   TBranch        *b_muon_thetaX;   //!
   TBranch        *b_muon_thetaY;   //!
   TBranch        *b_neutronecaledep;   //!
   TBranch        *b_neutronecalo;   //!
   TBranch        *b_neutronhcaledep;   //!
   TBranch        *b_neutronidetedep;   //!
   TBranch        *b_neutronnukeedep;   //!
   TBranch        *b_neutronodetedep;   //!
   TBranch        *b_neutronothersubdetedep;   //!
   TBranch        *b_neutronsideedep;   //!
   TBranch        *b_neutrontrkredep;   //!
   TBranch        *b_nke;   //!
   TBranch        *b_ntgtevis;   //!
   TBranch        *b_oangle;   //!
   TBranch        *b_oangle0;   //!
   TBranch        *b_oangle0x;   //!
   TBranch        *b_od_downstreamFrame;   //!
   TBranch        *b_od_downstreamFrame_z;   //!
   TBranch        *b_od_highStory;   //!
   TBranch        *b_od_highStory_t;   //!
   TBranch        *b_od_lowStory;   //!
   TBranch        *b_od_lowStory_t;   //!
   TBranch        *b_od_maxEnergy;   //!
   TBranch        *b_od_upstreamFrame;   //!
   TBranch        *b_od_upstreamFrame_z;   //!
   TBranch        *b_otherevis;   //!
   TBranch        *b_phys_energy_dispersed;   //!
   TBranch        *b_phys_energy_in_road_downstream;   //!
   TBranch        *b_phys_energy_in_road_upstream;   //!
   TBranch        *b_phys_energy_unattached;   //!
   TBranch        *b_pi0_evis_dispersed_blob;   //!
   TBranch        *b_pi0_evis_muon_blob;   //!
   TBranch        *b_pi0_evis_outtime_blob;   //!
   TBranch        *b_pi0_evis_vtx_blob;   //!
   TBranch        *b_pi0_evisfrac_dispersed_blob;   //!
   TBranch        *b_pi0_evisfrac_muon_blob;   //!
   TBranch        *b_pi0_evisfrac_outtime_blob;   //!
   TBranch        *b_pi0_evisfrac_vtx_blob;   //!
   TBranch        *b_pi0ecaledep;   //!
   TBranch        *b_pi0ecalo;   //!
   TBranch        *b_pi0hcaledep;   //!
   TBranch        *b_pi0idetedep;   //!
   TBranch        *b_pi0nukeedep;   //!
   TBranch        *b_pi0odetedep;   //!
   TBranch        *b_pi0othersubdetedep;   //!
   TBranch        *b_pi0sideedep;   //!
   TBranch        *b_pi0trkredep;   //!
   TBranch        *b_pienergy;   //!
   TBranch        *b_pienergy0;   //!
   TBranch        *b_pimlength;   //!
   TBranch        *b_piphi;   //!
   TBranch        *b_piphi0;   //!
   TBranch        *b_piplength;   //!
   TBranch        *b_pitheta;   //!
   TBranch        *b_pitheta0;   //!
   TBranch        *b_pke;   //!
   TBranch        *b_prim_vtx_smallest_opening_angle;   //!
   TBranch        *b_totalevis;   //!
   TBranch        *b_trkrevis;   //!
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
   TBranch        *b_multiplicities_sz;   //!
   TBranch        *b_multiplicities;   //!
   TBranch        *b_primary_truth_counts_sz;   //!
   TBranch        *b_primary_truth_counts;   //!
   TBranch        *b_primary_truth_pdgs1_sz;   //!
   TBranch        *b_primary_truth_pdgs1;   //!
   TBranch        *b_primary_truth_pdgs2_sz;   //!
   TBranch        *b_primary_truth_pdgs2;   //!
   TBranch        *b_primary_truth_pdgs3_sz;   //!
   TBranch        *b_primary_truth_pdgs3;   //!
   TBranch        *b_RE_photon_direction_1_sz;   //!
   TBranch        *b_RE_photon_direction_1;   //!
   TBranch        *b_RE_photon_direction_2_sz;   //!
   TBranch        *b_RE_photon_direction_2;   //!
   TBranch        *b_RE_photon_vertex_1_sz;   //!
   TBranch        *b_RE_photon_vertex_1;   //!
   TBranch        *b_RE_photon_vertex_2_sz;   //!
   TBranch        *b_RE_photon_vertex_2;   //!
   TBranch        *b_deviations_sz;   //!
   TBranch        *b_deviations;   //!
   TBranch        *b_g1convpos_sz;   //!
   TBranch        *b_g1convpos;   //!
   TBranch        *b_g1mom_sz;   //!
   TBranch        *b_g1mom;   //!
   TBranch        *b_g1mom0_sz;   //!
   TBranch        *b_g1mom0;   //!
   TBranch        *b_g2convpos_sz;   //!
   TBranch        *b_g2convpos;   //!
   TBranch        *b_g2mom_sz;   //!
   TBranch        *b_g2mom;   //!
   TBranch        *b_g2mom0_sz;   //!
   TBranch        *b_g2mom0;   //!
   TBranch        *b_michel_mom_sz;   //!
   TBranch        *b_michel_mom;   //!
   TBranch        *b_michel_pos_sz;   //!
   TBranch        *b_michel_pos;   //!
   TBranch        *b_mumom_sz;   //!
   TBranch        *b_mumom;   //!
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
   TBranch        *b_pimom_sz;   //!
   TBranch        *b_pimom;   //!
   TBranch        *b_pimom0_sz;   //!
   TBranch        *b_pimom0;   //!
   TBranch        *b_primary_separations_sz;   //!
   TBranch        *b_primary_separations;   //!
   TBranch        *b_primary_trklengths_sz;   //!
   TBranch        *b_primary_trklengths;   //!
   TBranch        *b_primary_truth_fractions1_sz;   //!
   TBranch        *b_primary_truth_fractions1;   //!
   TBranch        *b_primary_truth_fractions2_sz;   //!
   TBranch        *b_primary_truth_fractions2;   //!
   TBranch        *b_primary_truth_fractions3_sz;   //!
   TBranch        *b_primary_truth_fractions3;   //!
   TBranch        *b_primary_truth_shareds_sz;   //!
   TBranch        *b_primary_truth_shareds;   //!
   TBranch        *b_truth_has_physics_event;   //!
   TBranch        *b_truth_reco_minos_match;   //!
   TBranch        *b_truth_is_fiducial;   //!
   TBranch        *b_truth_pass_plausible;   //!
   TBranch        *b_truth_is_ccpi0;   //!
   TBranch        *b_truth_is_cc1pi0;   //!
   TBranch        *b_truth_is_ccpi0secondary;   //!
   TBranch        *b_truth_is_by_pim;   //!
   TBranch        *b_truth_is_ccpi0x;   //!
   TBranch        *b_truth_is_other;   //!
   TBranch        *b_truth_MC_photon_energy_1;   //!
   TBranch        *b_truth_MC_photon_energy_2;   //!
   TBranch        *b_truth_MC_pi0_energy;   //!
   TBranch        *b_truth_MC_scalar;   //!
   TBranch        *b_truth_fslepton_E;   //!
   TBranch        *b_truth_fslepton_P;   //!
   TBranch        *b_truth_fslepton_T;   //!
   TBranch        *b_truth_fslepton_phi;   //!
   TBranch        *b_truth_fslepton_theta;   //!
   TBranch        *b_truth_fslepton_theta_x;   //!
   TBranch        *b_truth_fslepton_theta_y;   //!
   TBranch        *b_truth_MC_photon_direction_1_sz;   //!
   TBranch        *b_truth_MC_photon_direction_1;   //!
   TBranch        *b_truth_MC_photon_direction_2_sz;   //!
   TBranch        *b_truth_MC_photon_direction_2;   //!
   TBranch        *b_truth_MC_photon_vertex_1_sz;   //!
   TBranch        *b_truth_MC_photon_vertex_1;   //!
   TBranch        *b_truth_MC_photon_vertex_2_sz;   //!
   TBranch        *b_truth_MC_photon_vertex_2;   //!
   TBranch        *b_truth_MC_pi0_momentum_sz;   //!
   TBranch        *b_truth_MC_pi0_momentum;   //!
   TBranch        *b_CCPi0AnaTool_nuFlavor;   //!
   TBranch        *b_CCPi0AnaTool_nuHelicity;   //!
   TBranch        *b_CCPi0AnaTool_intCurrent;   //!
   TBranch        *b_CCPi0AnaTool_intType;   //!
   TBranch        *b_CCPi0AnaTool_E;   //!
   TBranch        *b_CCPi0AnaTool_Q2;   //!
   TBranch        *b_CCPi0AnaTool_x;   //!
   TBranch        *b_CCPi0AnaTool_y;   //!
   TBranch        *b_CCPi0AnaTool_W;   //!
   TBranch        *b_CCPi0AnaTool_score;   //!
   TBranch        *b_CCPi0AnaTool_leptonE;   //!
   TBranch        *b_CCPi0AnaTool_vtx;   //!
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

   CCPi0(TTree *tree=0);
   virtual ~CCPi0();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif
