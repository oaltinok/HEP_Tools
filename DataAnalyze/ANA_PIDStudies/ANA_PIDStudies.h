//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jul 18 11:41:31 2013 by ROOT version 5.30/00
// from TChain CCInclusiveReco/
//////////////////////////////////////////////////////////

#ifndef ANA_PIDStudies_h
#define ANA_PIDStudies_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

#include "PDG_List.h"
#include "Bin_List.h"

class ANA_PIDStudies {
public :
   // -------------------------------------------------------------------------
   //     Custom Functions
   //--------------------------------------------------------------------------
   
    
    
   // -------------------------------------------------------------------------
   //     void run(): Generates a .root file with selected histograms
   //         playlist -> address of the playlist
   //         filename -> file name for the output .root file
   //         cutFile -> file name for the cutFile that lists the applied analysis cuts
   //         readmeFile -> file name for the readme
   //--------------------------------------------------------------------------
    void run(string playlist, char* filename,string cutFile, string readmeFile);
    double getPercent(double nAll, double nCurrent);
    
    void print();
    
    void openFiles(string cutFile, string readmeFile);
    void closeFiles();
    void writeReadme();

   // -------------------------------------------------------------------------
   //     Default Functions
   //--------------------------------------------------------------------------
    ANA_PIDStudies();
    ~ANA_PIDStudies();

    void Init(string playlist, TChain* fChain);
    Int_t GetEntry(Long64_t entry);
    Long64_t LoadTree(Long64_t entry);
    Bool_t Notify();
    void Show(Long64_t entry);
    Int_t Cut(Long64_t entry);

   // -------------------------------------------------------------------------
   //     Analysis Variables
   //--------------------------------------------------------------------------

   // -------------------------------------------------------------------------
   //     Files
   //--------------------------------------------------------------------------
    ofstream cutText;
    ofstream readme;


   // -------------------------------------------------------------------------
   //     Data
   //--------------------------------------------------------------------------
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

  // Declaration of leaf types
   Int_t           run;
   Int_t           subrun;
   Int_t           gate;
   Int_t           interaction;
   Int_t           n_particles;
   Int_t           part_pdg[2];   //[n_particles]
   Int_t           part_trackId[2];   //[n_particles]
   Int_t           part_fplane[2];   //[n_particles]
   Int_t           part_lplane[2];   //[n_particles]
   Int_t           part_destructCode[2];   //[n_particles]
   Double_t        part_beginPoint[2][3];   //[n_particles]
   Double_t        part_endPoint[2][3];   //[n_particles]
   Double_t        part_theta[2];   //[n_particles]
   Double_t        part_phi[2];   //[n_particles]
   Double_t        part_4p[2][4];   //[n_particles]
   Double_t        part_endTheta[2];   //[n_particles]
   Double_t        part_endPhi[2];   //[n_particles]
   Double_t        part_end4p[2][4];   //[n_particles]
   Double_t        part_lastPlaneDistanceZ[2];   //[n_particles]
   Int_t           n_prongs;
   Int_t           prong_isForked[3];   //[n_prongs]
   Int_t           prong_isKinked[3];   //[n_prongs]
   Int_t           prong_isDSEcal[3];   //[n_prongs]
   Int_t           prong_isSideEcal[3];   //[n_prongs]
   Int_t           prong_isNuclTargs[3];   //[n_prongs]
   Int_t           prong_priTrkPatRec[3];   //[n_prongs]
   Int_t           prong_fplane[3];   //[n_prongs]
   Int_t           prong_lplane[3];   //[n_prongs]
   Int_t           prong_hasSplitClusters[3];   //[n_prongs]
   Int_t           prong_hasEndSplitClusters[3];   //[n_prongs]
   Double_t        prong_beginPoint[3][3];   //[n_prongs]
   Double_t        prong_endPoint[3][3];   //[n_prongs]
   Double_t        prong_beginTheta[3];   //[n_prongs]
   Double_t        prong_beginPhi[3];   //[n_prongs]
   Double_t        prong_endTheta[3];   //[n_prongs]
   Double_t        prong_endPhi[3];   //[n_prongs]
   Double_t        prong_endm1Theta[3];   //[n_prongs]
   Double_t        prong_endm1Phi[3];   //[n_prongs]
   Int_t           tm_partProngMatchPDG[3];   //[n_prongs]
   Int_t           tm_partProngMatchTrackId[3];   //[n_prongs]
   Int_t           tm_partTrackMatchPDG[3];   //[n_prongs]
   Int_t           tm_partTrackMatchTrackId[3];   //[n_prongs]
   Double_t        tm_EFraction[3];   //[n_prongs]
   Double_t        dedx_pion_p[3];   //[n_prongs]
   Double_t        dedx_pion_fitpmin[3];   //[n_prongs]
   Double_t        dedx_pion_fitpmax[3];   //[n_prongs]
   Double_t        dedx_pion_score1[3];   //[n_prongs]
   Double_t        dedx_pion_score2[3];   //[n_prongs]
   Int_t           dedx_pion_fit_by_range[3];   //[n_prongs]
   Int_t           dedx_pion_fit_fails[3];   //[n_prongs]
   Int_t           dedx_pion_nnodes[3];   //[n_prongs]
   Double_t        dedx_pion_chi2[3];   //[n_prongs]
   Double_t        dedx_pion_range_diff[3];   //[n_prongs]
   Double_t        dedx_pion_measuredDEDX[3][100];   //[n_prongs]
   Double_t        dedx_pion_calculatedDEDX[3][100];   //[n_prongs]
   Double_t        dedx_pion_zpositions[3][100];   //[n_prongs]
   Double_t        dedx_proton_p[3];   //[n_prongs]
   Double_t        dedx_proton_score1[3];   //[n_prongs]
   Double_t        dedx_proton_score2[3];   //[n_prongs]
   Int_t           dedx_proton_fit_by_range[3];   //[n_prongs]
   Int_t           dedx_proton_fit_fails[3];   //[n_prongs]
   Int_t           dedx_proton_nnodes[3];   //[n_prongs]
   Double_t        dedx_proton_chi2[3];   //[n_prongs]
   Double_t        dedx_proton_range_diff[3];   //[n_prongs]
   Double_t        dedx_proton_measuredDEDX[3][100];   //[n_prongs]
   Double_t        dedx_proton_calculatedDEDX[3][100];   //[n_prongs]
   Double_t        dedx_proton_zpositions[3][100];   //[n_prongs]
   Double_t        spid_proton_score[3];   //[n_prongs]
   Double_t        spid_piminus_score[3];   //[n_prongs]
   Double_t        spid_piplus_score[3];   //[n_prongs]

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_subrun;   //!
   TBranch        *b_gate;   //!
   TBranch        *b_interaction;   //!
   TBranch        *b_n_particles;   //!
   TBranch        *b_part_pdg;   //!
   TBranch        *b_part_trackId;   //!
   TBranch        *b_part_fplane;   //!
   TBranch        *b_part_lplane;   //!
   TBranch        *b_part_destructCode;   //!
   TBranch        *b_part_beginPoint;   //!
   TBranch        *b_part_endPoint;   //!
   TBranch        *b_part_theta;   //!
   TBranch        *b_part_phi;   //!
   TBranch        *b_part_4p;   //!
   TBranch        *b_part_endTheta;   //!
   TBranch        *b_part_endPhi;   //!
   TBranch        *b_part_end4p;   //!
   TBranch        *b_part_lastPlaneDistanceZ;   //!
   TBranch        *b_n_prongs;   //!
   TBranch        *b_prong_isForked;   //!
   TBranch        *b_prong_isKinked;   //!
   TBranch        *b_prong_isDSEcal;   //!
   TBranch        *b_prong_isSideEcal;   //!
   TBranch        *b_prong_isNuclTargs;   //!
   TBranch        *b_prong_priTrkPatRec;   //!
   TBranch        *b_prong_fplane;   //!
   TBranch        *b_prong_lplane;   //!
   TBranch        *b_prong_hasSplitClusters;   //!
   TBranch        *b_prong_hasEndSplitClusters;   //!
   TBranch        *b_prong_beginPoint;   //!
   TBranch        *b_prong_endPoint;   //!
   TBranch        *b_prong_beginTheta;   //!
   TBranch        *b_prong_beginPhi;   //!
   TBranch        *b_prong_endTheta;   //!
   TBranch        *b_prong_endPhi;   //!
   TBranch        *b_prong_endm1Theta;   //!
   TBranch        *b_prong_endm1Phi;   //!
   TBranch        *b_tm_partProngMatchPDG;   //!
   TBranch        *b_tm_partProngMatchTrackId;   //!
   TBranch        *b_tm_partTrackMatchPDG;   //!
   TBranch        *b_tm_partTrackMatchTrackId;   //!
   TBranch        *b_tm_EFraction;   //!
   TBranch        *b_dedx_pion_p;   //!
   TBranch        *b_dedx_pion_fitpmin;   //!
   TBranch        *b_dedx_pion_fitpmax;   //!
   TBranch        *b_dedx_pion_score1;   //!
   TBranch        *b_dedx_pion_score2;   //!
   TBranch        *b_dedx_pion_fit_by_range;   //!
   TBranch        *b_dedx_pion_fit_fails;   //!
   TBranch        *b_dedx_pion_nnodes;   //!
   TBranch        *b_dedx_pion_chi2;   //!
   TBranch        *b_dedx_pion_range_diff;   //!
   TBranch        *b_dedx_pion_measuredDEDX;   //!
   TBranch        *b_dedx_pion_calculatedDEDX;   //!
   TBranch        *b_dedx_pion_zpositions;   //!
   TBranch        *b_dedx_proton_p;   //!
   TBranch        *b_dedx_proton_score1;   //!
   TBranch        *b_dedx_proton_score2;   //!
   TBranch        *b_dedx_proton_fit_by_range;   //!
   TBranch        *b_dedx_proton_fit_fails;   //!
   TBranch        *b_dedx_proton_nnodes;   //!
   TBranch        *b_dedx_proton_chi2;   //!
   TBranch        *b_dedx_proton_range_diff;   //!
   TBranch        *b_dedx_proton_measuredDEDX;   //!
   TBranch        *b_dedx_proton_calculatedDEDX;   //!
   TBranch        *b_dedx_proton_zpositions;   //!
   TBranch        *b_spid_proton_score;   //!
   TBranch        *b_spid_piminus_score;   //!
   TBranch        *b_spid_piplus_score;   //!

   ANA_PIDStudies(TTree *tree=0);
   virtual ~ANA_PIDStudies();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif
