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
        //      nTopologies defined in CCProtonPi0_NTupleAnalysis
        //--------------------------------------------------------------------------
        vector<CCProtonPi0_Cut> nCut_All;
        vector<CCProtonPi0_Cut> nCut_Vertex_None;
        vector<CCProtonPi0_Cut> nCut_Vertex_Not_Reconstructable; 
        vector<CCProtonPi0_Cut> nCut_Vertex_Not_Fiducial;
        vector<CCProtonPi0_Cut> nCut_Muon_None;              
        vector<CCProtonPi0_Cut> nCut_Muon_Charge;
        vector<CCProtonPi0_Cut> nCut_Vertex_Michel_Exist; 
        vector<CCProtonPi0_Cut> nCut_EndPoint_Michel_Exist;
        vector<CCProtonPi0_Cut> nCut_secEndPoint_Michel_Exist;
        vector<CCProtonPi0_Cut> nCut_PreFilter_Pi0;
        vector<CCProtonPi0_Cut> nCut_ConeBlobs;
        vector<CCProtonPi0_Cut> nCut_BlobDirectionBad;
        vector<CCProtonPi0_Cut> nCut_BlobsBad;
        vector<CCProtonPi0_Cut> nCut_Photon1DistanceLow;
        vector<CCProtonPi0_Cut> nCut_Photon2DistanceLow;
        vector<CCProtonPi0_Cut> nCut_Pi0_invMass;
        vector<CCProtonPi0_Cut> nCut_Particle_None;
        vector<CCProtonPi0_Cut> nCut_Proton_None;            
        vector<CCProtonPi0_Cut> nCut_ProtonScore;
        vector<CCProtonPi0_Cut> nCut_DeltaInvMass;
        vector<CCProtonPi0_Cut> nCut_beamEnergy;
        vector<CCProtonPi0_Cut> nCut_UnusedE;
        
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
        vector<MnvH1D*> hCut_nProngs;
        vector<MnvH1D*> hCut_nProngs2;
        
        // 1Track
        vector<MnvH1D*> hCut_1Track_Michel;
        vector<MnvH1D*> hCut_1Track_eVis_nuclearTarget;
        vector<MnvH1D*> hCut_1Track_eVis_other;
        vector<MnvH1D*> hCut_1Track_pi0invMass;
        vector<MnvH1D*> hCut_1Track_gamma1ConvDist;
        vector<MnvH1D*> hCut_1Track_gamma2ConvDist;
        vector<MnvH1D*> hCut_1Track_neutrinoE;
        vector<MnvH1D*> hCut_1Track_UnusedE;

        // 2Track
        vector<MnvH1D*> hCut_2Track_Michel;
        vector<MnvH1D*> hCut_2Track_eVis_nuclearTarget;
        vector<MnvH1D*> hCut_2Track_eVis_other;
        vector<MnvH1D*> hCut_2Track_pi0invMass;
        vector<MnvH1D*> hCut_2Track_gamma1ConvDist;
        vector<MnvH1D*> hCut_2Track_gamma2ConvDist;
        vector<MnvH1D*> hCut_2Track_UnusedE;
        vector<MnvH1D*> hCut_2Track_neutrinoE;
        vector<MnvH1D*> hCut_2Track_protonScore_pIDDiff;
        vector<MnvH1D*> hCut_2Track_protonScore_LLR;
        vector<MnvH1D*> hCut_2Track_deltaInvMass;
 
        // MC Only Histograms
        TH1D* mc_w_DIS;
        TH1D* mc_w_RES;
        TH1D* mc_w_CCQE;

        // Pi0 Invariant Mass
        TH1D* pi0_invMass_1Track;
        TH1D* pi0_invMass_2Track;

    private:
        void initHistograms();
        void init_nCutVectors();
        void SetCutNames();
        void OpenTextFiles(bool isMC);
        void formCutVectors();
        void writeCutTableHeader(int t);
        void writeCutTableRows(int t, vector< vector<CCProtonPi0_Cut> > &nCutVector);
        void writeSingleRow(int t, CCProtonPi0_Cut& currentCut);
        double getCutEfficiency(CCProtonPi0_Cut& currentCut, CCProtonPi0_Cut& effBase) const;
        double getCutEfficiency(CCProtonPi0_Cut& currentCut, double effBase) const;
        double getCutPurity(CCProtonPi0_Cut& currentCut) const;

        vector< vector<CCProtonPi0_Cut> > nCutVector_Common;
        vector< vector<CCProtonPi0_Cut> > nCutVector_Topology;
        
        CCProtonPi0_BinList binList;
        
        ofstream cutText[nTopologies];
        string cutFile[nTopologies];
        
        TFile* f;
        string rootDir;

        // Number of Signal Events from Truth Info
        double nTrueSignal;
};




#endif
