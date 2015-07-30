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
        //     CCProtonPi0_Cut Numbers
        //--------------------------------------------------------------------------
        // Common CCProtonPi0_Cut Numbers
        CCProtonPi0_Cut nCut_All;
        CCProtonPi0_Cut nCut_Vertex_None;
        CCProtonPi0_Cut nCut_Vertex_Not_Reconstructable; 
        CCProtonPi0_Cut nCut_Vertex_Not_Fiducial;
        CCProtonPi0_Cut nCut_Muon_None;              
        CCProtonPi0_Cut nCut_Muon_Not_Plausible;
        CCProtonPi0_Cut nCut_Muon_Charge;
        CCProtonPi0_Cut nCut_Vertex_Michel_Exist; 
        CCProtonPi0_Cut nCut_EndPoint_Michel_Exist;
        CCProtonPi0_Cut nCut_secEndPoint_Michel_Exist;
        CCProtonPi0_Cut nCut_PreFilter_Pi0;
        CCProtonPi0_Cut nCut_ConeBlobs;
        CCProtonPi0_Cut nCut_BlobDirectionBad;
        CCProtonPi0_Cut nCut_BlobsBad;
        CCProtonPi0_Cut nCut_Photon1DistanceLow;
        CCProtonPi0_Cut nCut_Photon2DistanceLow;
        CCProtonPi0_Cut nCut_Pi0_invMass;
        
        // nProngs == 1 Cut Numbers (Muon + Pi0 )
        CCProtonPi0_Cut nCut_1Prong_Particle_None;
        CCProtonPi0_Cut nCut_1Prong_Proton_None;            
        CCProtonPi0_Cut nCut_1Prong_ProtonScore;
        CCProtonPi0_Cut nCut_1Prong_DeltaInvMass;
        CCProtonPi0_Cut nCut_1Prong_beamEnergy;
        CCProtonPi0_Cut nCut_1Prong_UnusedE;
        
        // nProngs >= 2 Cut Numbers (Muon + Pi0 + X(No Meson))
        CCProtonPi0_Cut nCut_2Prong_Particle_None;
        CCProtonPi0_Cut nCut_2Prong_Proton_None;            
        CCProtonPi0_Cut nCut_2Prong_ProtonScore;
        CCProtonPi0_Cut nCut_2Prong_DeltaInvMass;
        CCProtonPi0_Cut nCut_2Prong_beamEnergy;
        CCProtonPi0_Cut nCut_2Prong_UnusedE;
 
        // -------------------------------------------------------------------------
        // Cut Histograms
        // -------------------------------------------------------------------------
        // Common
        vector<MnvH1D*> hCut_Michel;
        vector<MnvH1D*> hCut_eVis_nuclearTarget;
        vector<MnvH1D*> hCut_eVis_other;
        vector<MnvH1D*> hCut_pi0invMass;
        vector<MnvH1D*> hCut_gamma1ConvDist;
        vector<MnvH1D*> hCut_gamma2ConvDist;
        
        vector<MnvH1D*> hCut_1Prong_neutrinoE;
        vector<MnvH1D*> hCut_2Prong_neutrinoE;
        vector<MnvH1D*> hCut_1Prong_UnusedE;
        vector<MnvH1D*> hCut_2Prong_UnusedE;

        // 2 Prong Specific
        vector<MnvH1D*> hCut_protonScore_pIDDiff;
        vector<MnvH1D*> hCut_protonScore_LLR;
        vector<MnvH1D*> hCut_deltaInvMass;
 
        // MC Only Histograms
        TH1D* mc_w_DIS;
        TH1D* mc_w_RES;
        TH1D* mc_w_CCQE;
              
    private:
        void initHistograms();
        void SetCutNames();
        void OpenOutputFile();
        void formCutVector();
        void writeCutTableHeader();
        void writeCutTableRows();
        void writeSingleRow(CCProtonPi0_Cut& currentCut);
        void writeSingleRow(CCProtonPi0_Cut& nCut_1Prong, CCProtonPi0_Cut& nCut_2Prong);
        double getCutEfficiency(CCProtonPi0_Cut& currentCut, CCProtonPi0_Cut& effBase) const;
        double getCutEfficiency(CCProtonPi0_Cut& currentCut, double effBase) const;
        double getCutPurity(CCProtonPi0_Cut& currentCut) const;

        vector<CCProtonPi0_Cut> nCutVector;
        
        CCProtonPi0_BinList binList;
        
        ofstream cutText;
        string cutFile;
        
        TFile* f;
        string rootDir;

        // Number of Signal Events from Truth Info
        double nTrueSignal;
};




#endif
