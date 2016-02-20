#ifndef CCProtonPi0_TruthAnalyzer_h
#define CCProtonPi0_TruthAnalyzer_h

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>

#include <TROOT.h>
#include <TChain.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <PlotUtils/MnvH1D.h>

#include "Cintex/Cintex.h"
#include "../../../Libraries/Folder_List.h"
#include "../../../Classes/NTupleAnalysis/CCProtonPi0_NTupleAnalysis.h"

using namespace PlotUtils;

class CCProtonPi0_TruthAnalyzer : public CCProtonPi0_NTupleAnalysis
{
    public:
        CCProtonPi0_TruthAnalyzer();
        ~CCProtonPi0_TruthAnalyzer();
        void Loop(std::string playlist);

        MnvH1D* muon_P_mc_truth_all_signal;
        MnvH1D* pi0_P_mc_truth_all_signal;
        MnvH1D* muon_theta_mc_truth_all_signal;
        MnvH1D* pi0_theta_mc_truth_all_signal;
        MnvH1D* neutrino_E_mc_truth_all_signal;
        MnvH1D* QSq_mc_truth_all_signal;

    private :
        TFile* f;
        std::string rootDir;

        void FillVertErrorBand_Flux(MnvH1D* h, double var);
        void FillVertErrorBand_Genie(MnvH1D* h, double var);
        void Calc_WeightFromSystematics();
        void AddOtherErrorBands_FillWithCV();
        void AddErrorBands_FillWithCV(MnvH1D* hist);
        double GetPercent(double nAll, double nOther);
        void initHistograms();
        void openTextFiles();
        void resetCounters();
        void writeTextFile();
        void writeHistograms();
        void FillHistogram(MnvH1D *hist, double var);
        void FillSignalHistograms();

        
        // Default Functions    
        void     Init(std::string playlist, TChain* fChain);
        Int_t    GetEntry(Long64_t entry);
        Long64_t LoadTree(Long64_t entry);

        // Class Variables
        std::string file_name;
        ofstream textFile;

        // Counters
        double nAll;
        double nFidVol;
        double nNoFidVol;
        double nSignal;
        // Background With Pi0
        double nBckg_NoPi0;
        double nBckg_SinglePi0;
        double nBckg_MultiPi0;
        // Background Types
        double nBckg_NC;
        double nBckg_AntiNeutrino;
        double nBckg_QELike;
        double nBckg_SingleChargedPion;
        double nBckg_SingleChargedPion_ChargeExchanged;
        double nBckg_DoublePion_WithPi0;
        double nBckg_DoublePion_WithoutPi0;
        double nBckg_MultiPion_WithPi0;
        double nBckg_MultiPion_WithoutPi0;
        double nBckg_Other;

        double cvweight;

        // NTuple Truth Branch
        TTree          *fChain;   //!pointer to the analyzed TTree or TChain
        Int_t           fCurrent; //!current Tree number in a TChain

        // Declaration of leaf types
        Double_t        eventID;
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
        Double_t        mc_FSPartPx[106];   //[mc_nFSPart]
        Double_t        mc_FSPartPy[106];   //[mc_nFSPart]
        Double_t        mc_FSPartPz[106];   //[mc_nFSPart]
        Double_t        mc_FSPartE[106];   //[mc_nFSPart]
        Int_t           mc_FSPartPDG[106];   //[mc_nFSPart]
        Int_t           mc_er_nPart;
        Int_t           mc_er_ID[134];   //[mc_er_nPart]
        Int_t           mc_er_status[134];   //[mc_er_nPart]
        Double_t        mc_er_posInNucX[134];   //[mc_er_nPart]
        Double_t        mc_er_posInNucY[134];   //[mc_er_nPart]
        Double_t        mc_er_posInNucZ[134];   //[mc_er_nPart]
        Double_t        mc_er_Px[134];   //[mc_er_nPart]
        Double_t        mc_er_Py[134];   //[mc_er_nPart]
        Double_t        mc_er_Pz[134];   //[mc_er_nPart]
        Double_t        mc_er_E[134];   //[mc_er_nPart]
        Int_t           mc_er_FD[134];   //[mc_er_nPart]
        Int_t           mc_er_LD[134];   //[mc_er_nPart]
        Int_t           mc_er_mother[134];   //[mc_er_nPart]
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

        // List of branches
        TBranch        *b_eventID;   //!
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

};

#endif
