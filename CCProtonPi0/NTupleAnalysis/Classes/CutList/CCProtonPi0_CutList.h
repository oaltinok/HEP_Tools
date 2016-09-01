/*
================================================================================
Class: CCProtonPi0_CutList
    Member Variables are the CCProtonPi0_Cut Numbers which represent each Selection 
        in the Analysis
    Creates the CCProtonPi0_Cut Table for all topologies in the Analysis
    
    Main Directory:
        Classes/CutList
        
    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
================================================================================
*/
#ifndef CCProtonPi0_CutList_h
#define CCProtonPi0_CutList_h

#include "../NTupleAnalysis/CCProtonPi0_NTupleAnalysis.h"
#include "../Cut/CCProtonPi0_Cut.h"
#include "../BinList/CCProtonPi0_BinList.h"

#include "TVector.h"
#include "TGraph.h"

using namespace PlotUtils;

class CCProtonPi0_CutList : public CCProtonPi0_NTupleAnalysis 
{
    public:
        CCProtonPi0_CutList(bool isModeReduce, bool isMC);
        ~CCProtonPi0_CutList();
        
        void writeCutTable();
        void writeHistograms();
        // -------------------------------------------------------------------------
        // CCProtonPi0_Cut Numbers
        //--------------------------------------------------------------------------
        CCProtonPi0_Cut nCut_All;
        CCProtonPi0_Cut nCut_Vertex_None;
        CCProtonPi0_Cut nCut_Vertex_Not_Reconstructable; 
        CCProtonPi0_Cut nCut_Vertex_Not_Fiducial;
        CCProtonPi0_Cut nCut_Muon_None;              
        CCProtonPi0_Cut nCut_Muon_Angle;
        CCProtonPi0_Cut nCut_Muon_Charge;
        CCProtonPi0_Cut nCut_Vertex_Michel_Exist; 
        CCProtonPi0_Cut nCut_EndPoint_Michel_Exist;
        CCProtonPi0_Cut nCut_secEndPoint_Michel_Exist;
        CCProtonPi0_Cut nCut_Particle_None;
        CCProtonPi0_Cut nCut_Proton_None;            
        CCProtonPi0_Cut nCut_Proton_Bad;            
        CCProtonPi0_Cut nCut_ProtonScore;
        CCProtonPi0_Cut nCut_PreFilter_Pi0;
        CCProtonPi0_Cut nCut_ConeBlobs;
        CCProtonPi0_Cut nCut_BlobDirectionBad;
        CCProtonPi0_Cut nCut_Pi0_Bad;
        CCProtonPi0_Cut nCut_Shower_Michel_Exist;
        CCProtonPi0_Cut nCut_Photon1DistanceLow;
        CCProtonPi0_Cut nCut_Photon2DistanceLow;
        CCProtonPi0_Cut nCut_LowE_SmallAngle;
        CCProtonPi0_Cut nCut_Pi0_invMass;
        CCProtonPi0_Cut nCut_beamEnergy;
        CCProtonPi0_Cut nCut_W;

        // 1 Track Cuts
        CCProtonPi0_Cut nCut_1Track_All;
        CCProtonPi0_Cut nCut_1Track_PreFilter_Pi0;
        CCProtonPi0_Cut nCut_1Track_ConeBlobs;
        CCProtonPi0_Cut nCut_1Track_BlobDirectionBad;
        CCProtonPi0_Cut nCut_1Track_Pi0_Bad;
        CCProtonPi0_Cut nCut_1Track_Photon1DistanceLow;
        CCProtonPi0_Cut nCut_1Track_Photon2DistanceLow;
        CCProtonPi0_Cut nCut_1Track_Pi0_invMass;
        CCProtonPi0_Cut nCut_1Track_beamEnergy;
        CCProtonPi0_Cut nCut_1Track_W;

        // 2 Track Cuts
        CCProtonPi0_Cut nCut_2Track_All;
        CCProtonPi0_Cut nCut_2Track_ProtonScore;
        CCProtonPi0_Cut nCut_2Track_PreFilter_Pi0;
        CCProtonPi0_Cut nCut_2Track_ConeBlobs;
        CCProtonPi0_Cut nCut_2Track_BlobDirectionBad;
        CCProtonPi0_Cut nCut_2Track_Pi0_Bad;
        CCProtonPi0_Cut nCut_2Track_Photon1DistanceLow;
        CCProtonPi0_Cut nCut_2Track_Photon2DistanceLow;
        CCProtonPi0_Cut nCut_2Track_Pi0_invMass;
        CCProtonPi0_Cut nCut_2Track_beamEnergy;
        CCProtonPi0_Cut nCut_2Track_W;

        // -------------------------------------------------------------------------
        // Cut Histograms
        // -------------------------------------------------------------------------
        // Common
        vector<MnvH1D*> hCut_nVertices;
        vector<MnvH1D*> hCut_nTracks;
        vector<MnvH1D*> hCut_nTracks2;
        vector<MnvH1D*> hCut_nTracks_Close;
        vector<MnvH1D*> hCut_nTracks_Far;
        vector<MnvH1D*> hCut_nTracks_Discarded;
        vector<MnvH1D*> hCut_Michel;
        vector<MnvH1D*> hCut_nProtonCandidates;
        vector<MnvH1D*> hCut_nShowerCandidates;
        vector<MnvH1D*> hCut_pi0invMass;
        vector<MnvH1D*> SideBand_muon_P;
        vector<MnvH1D*> SideBand_muon_theta;
        vector<MnvH1D*> SideBand_pi0_P;
        vector<MnvH1D*> SideBand_pi0_KE;
        vector<MnvH1D*> SideBand_pi0_theta;
        vector<MnvH1D*> SideBand_neutrino_E;
        vector<MnvH1D*> SideBand_QSq;
        vector<MnvH1D*> SideBand_W;

        // 1Track
        vector<MnvH1D*> hCut_1Track_nShowerCandidates;
        vector<MnvH1D*> hCut_1Track_eVis_nuclearTarget;
        vector<MnvH1D*> hCut_1Track_eVis_other;
        vector<MnvH1D*> hCut_1Track_pi0invMass;
        vector<MnvH1D*> hCut_1Track_gamma1ConvDist;
        vector<MnvH1D*> hCut_1Track_gamma2ConvDist;
        vector<MnvH1D*> hCut_1Track_neutrinoE;

        // 2Track
        vector<MnvH1D*> hCut_2Track_nShowerCandidates;
        vector<MnvH1D*> hCut_2Track_eVis_nuclearTarget;
        vector<MnvH1D*> hCut_2Track_eVis_other;
        vector<MnvH1D*> hCut_2Track_pi0invMass;
        vector<MnvH1D*> hCut_2Track_gamma1ConvDist;
        vector<MnvH1D*> hCut_2Track_gamma2ConvDist;
        vector<MnvH1D*> hCut_2Track_neutrinoE;
        vector<MnvH1D*> hCut_2Track_protonScore_LLR;
        vector<MnvH1D*> hCut_2Track_deltaInvMass;
 
        // Selected Signal 
        // Signal Q2
        TH1D* mc_Q2_QE;

        TH1D* mc_Q2_RES_1232;
        TH1D* mc_Q2_RES_1535;
        TH1D* mc_Q2_RES_1520;
        TH1D* mc_Q2_RES_Other;

        TH1D* mc_Q2_DIS_1_pi;
        TH1D* mc_Q2_DIS_2_pi;
        TH1D* mc_Q2_DIS_Multi_pi;
        TH1D* mc_Q2_Non_RES;

        // Signal incomingE
        TH1D* mc_incomingE_QE;

        TH1D* mc_incomingE_RES_1232;
        TH1D* mc_incomingE_RES_1535;
        TH1D* mc_incomingE_RES_1520;
        TH1D* mc_incomingE_RES_Other;

        TH1D* mc_incomingE_DIS_1_pi;
        TH1D* mc_incomingE_DIS_2_pi;
        TH1D* mc_incomingE_DIS_Multi_pi;
        TH1D* mc_incomingE_Non_RES;

        // Signal w
        TH1D* mc_w_QE;

        TH1D* mc_w_RES_1232;
        TH1D* mc_w_RES_1535;
        TH1D* mc_w_RES_1520;
        TH1D* mc_w_RES_Other;

        TH1D* mc_w_DIS_1_pi;
        TH1D* mc_w_DIS_2_pi;
        TH1D* mc_w_DIS_Multi_pi;
        TH1D* mc_w_Non_RES;

        // Pi0 Invariant Mass
        TH1D* pi0_invMass_1Track;
        TH1D* pi0_invMass_2Track;
      

        // Michel Electron - Truth Match
        TH1D* michel_piplus_time_diff;
        TH1D* michel_neutron_time_diff;
        TH1D* michel_proton_time_diff;
        TH1D* michel_piminus_time_diff;
        TH1D* michel_other_time_diff;

        TH1D* michel_piplus_energy;
        TH1D* michel_neutron_energy;
        TH1D* michel_proton_energy;
        TH1D* michel_piminus_energy;
        TH1D* michel_other_energy;

        TH1D* michel_piplus_distance;
        TH1D* michel_neutron_distance;
        TH1D* michel_proton_distance;
        TH1D* michel_piminus_distance;
        TH1D* michel_other_distance;

        TH1D* michel_piplus_distance_z;
        TH1D* michel_neutron_distance_z;
        TH1D* michel_proton_distance_z;
        TH1D* michel_piminus_distance_z;
        TH1D* michel_other_distance_z;

        // Pi0 Invariant Mass - Truth Match
        TH1D* signal_invMass_pizero;
        TH1D* signal_invMass_piplus;
        TH1D* signal_invMass_proton;
        TH1D* signal_invMass_neutron;
        TH1D* signal_invMass_other;
 
        TH1D* background_invMass_pizero;
        TH1D* background_invMass_piplus;
        TH1D* background_invMass_proton;
        TH1D* background_invMass_neutron;
        TH1D* background_invMass_other;
      
        // Background Subtraction
        MnvH1D* invMass_all;
        MnvH1D* invMass_mc_reco_all;
        MnvH1D* invMass_mc_reco_signal;
        MnvH1D* invMass_mc_reco_bckg;

        // Studies
        TH2D* signal_gamma_E_cos_openingAngle;
        TH2D* bckg_gamma_E_cos_openingAngle;
        TH2D* bckg_signal_diff_E_cos_openingAngle;
        
        TH3D* signal_E_cosTheta_convLength;
        TH3D* bckg_E_cosTheta_convLength;
        TH3D* bckg_signal_diff_E_cosTheta_convLength;

    private:
        void initHistograms();
        void SetCutNames();
        void OpenTextFiles(bool isMC);
        void formCutVectors();
        void writeAllCuts();
        void write1TrackCuts();
        void write2TrackCuts();
        void writeCutTableHeader(std::ofstream &file);
        void writeCutTableRows(std::ofstream &file, vector<CCProtonPi0_Cut> &nCutVector, bool isAll);
        void writeSingleRow(std::ofstream &file, CCProtonPi0_Cut &currentCut, CCProtonPi0_Cut &eff_base_all, CCProtonPi0_Cut &eff_base_MINOS);
        double getCutEfficiency(CCProtonPi0_Cut &currentCut, CCProtonPi0_Cut &effBase) const;
        double getCutEfficiency(CCProtonPi0_Cut &currentCut, double effBase) const;
        double getCutPurity(CCProtonPi0_Cut &currentCut) const;

        bool m_isMC;
        vector<CCProtonPi0_Cut> nCutVector_All;
        vector<CCProtonPi0_Cut> nCutVector_1Track;
        vector<CCProtonPi0_Cut> nCutVector_2Track;
        
        CCProtonPi0_BinList binList;
        
        std::ofstream cutText_All;
        std::ofstream cutText_1Track;
        std::ofstream cutText_2Track;
        
        TFile* f;
        string rootDir;

        // Number of Signal Events from Truth Info
        double nTrueSignal;
        bool use_nTrueSignal;
    
};




#endif
