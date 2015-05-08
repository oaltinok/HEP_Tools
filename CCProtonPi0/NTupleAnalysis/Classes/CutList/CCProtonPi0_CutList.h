/*
================================================================================
Class: CCProtonPi0_CutList
    Member Variables are the CCProtonPi0_Cut Numbers which represent each Selection 
        in the Analysis
    Creates the CCProtonPi0_Cut Table for all topologies in the Analysis
    
    Main Directory:
        Classes/CutList
        
    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision:  2015_05_07
================================================================================
*/
#ifndef CCProtonPi0_CutList_h
#define CCProtonPi0_CutList_h

#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <string>

#include "../NTupleAnalysis/CCProtonPi0_NTupleAnalysis.h"
#include "../Cut/CCProtonPi0_Cut.h"

using namespace std;

class CCProtonPi0_CutList : public CCProtonPi0_NTupleAnalysis
{
    public:
        CCProtonPi0_CutList(int nMode);
        ~CCProtonPi0_CutList();
        
        void writeCutTable();
      
        // -------------------------------------------------------------------------
        //     CCProtonPi0_Cut Numbers
        //--------------------------------------------------------------------------
        // Common CCProtonPi0_Cut Numbers
        CCProtonPi0_Cut nCut_All;
        CCProtonPi0_Cut nCut_Vertex_None;
        CCProtonPi0_Cut nCut_Vertex_Not_Reconstructable; 
        CCProtonPi0_Cut nCut_Vertex_Not_Fiducial;
        CCProtonPi0_Cut nCut_Vertex_Count;
        CCProtonPi0_Cut nCut_nProngs;
        
        // nProngs == 1 CCProtonPi0_Cut Numbers
        CCProtonPi0_Cut nCut_1Prong_Muon_None;              
        CCProtonPi0_Cut nCut_1Prong_Muon_Not_Plausible;
        CCProtonPi0_Cut nCut_1Prong_Muon_Score_Low;
        CCProtonPi0_Cut nCut_1Prong_Muon_Charge;
        CCProtonPi0_Cut nCut_1Prong_Vertex_Michel_Exist; 
        CCProtonPi0_Cut nCut_1Prong_EndPoint_Michel_Exist;
        CCProtonPi0_Cut nCut_1Prong_secEndPoint_Michel_Exist;
        CCProtonPi0_Cut nCut_1Prong_PreFilter_Pi0;
        CCProtonPi0_Cut nCut_1Prong_VtxBlob;
        CCProtonPi0_Cut nCut_1Prong_ConeBlobs;
        CCProtonPi0_Cut nCut_1Prong_Pi0_invMass;
        CCProtonPi0_Cut nCut_1Prong_Photon1DistanceLow;
        CCProtonPi0_Cut nCut_1Prong_Photon2DistanceLow;
        CCProtonPi0_Cut nCut_1Prong_beamEnergy;
        CCProtonPi0_Cut nCut_1Prong_UnusedE;
        
        // nProngs == 2 CCProtonPi0_Cut Numbers
        CCProtonPi0_Cut nCut_2Prong_Muon_None;              
        CCProtonPi0_Cut nCut_2Prong_Muon_Not_Plausible;
        CCProtonPi0_Cut nCut_2Prong_Muon_Score_Low;
        CCProtonPi0_Cut nCut_2Prong_Muon_Charge;
        CCProtonPi0_Cut nCut_2Prong_Vertex_Michel_Exist; 
        CCProtonPi0_Cut nCut_2Prong_EndPoint_Michel_Exist;
        CCProtonPi0_Cut nCut_2Prong_secEndPoint_Michel_Exist;
        CCProtonPi0_Cut nCut_2Prong_PreFilter_Pi0;
        CCProtonPi0_Cut nCut_2Prong_VtxBlob;
        CCProtonPi0_Cut nCut_2Prong_ConeBlobs;
        CCProtonPi0_Cut nCut_2Prong_Pi0_invMass;
        CCProtonPi0_Cut nCut_2Prong_Photon1DistanceLow;
        CCProtonPi0_Cut nCut_2Prong_Photon2DistanceLow;
        CCProtonPi0_Cut nCut_2Prong_beamEnergy;
        CCProtonPi0_Cut nCut_2Prong_UnusedE;
        CCProtonPi0_Cut nCut_2Prong_Particle_None;
        CCProtonPi0_Cut nCut_2Prong_Proton_None;            
        CCProtonPi0_Cut nCut_2Prong_ProtonScore;
        CCProtonPi0_Cut nCut_2Prong_DeltaInvMass;
        
    private:
        void SetCutNames();
        void OpenTextFiles();
        void formCutVectors();
        void writeCutTableHeader();
        void writeCutTableRows(vector<CCProtonPi0_Cut> nCutVector, int nProngs, bool isShortList);
        double getCutEfficiency(double nSig, double effBase);
        double getCutPurity(double nSig, double nEvents);
        
        vector<CCProtonPi0_Cut> nCutVector_1Prong;
        vector<CCProtonPi0_Cut> nCutVector_2Prong;
        vector<CCProtonPi0_Cut> nCutVector_1Prong_ShortList;
        vector<CCProtonPi0_Cut> nCutVector_2Prong_ShortList;
        
        ofstream cutText[nTopologies];
        string cutFile[nTopologies];
        
};




#endif
