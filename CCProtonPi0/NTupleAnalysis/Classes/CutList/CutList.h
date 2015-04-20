/*
================================================================================
Class: CutList
    Member Variables are the Cut Numbers which represent each Selection 
        in the Analysis
    Creates the Cut Table for all topologies in the Analysis
    
    Main Directory:
        Classes/CutList
        
    Main Directory:
        Classes/CutList
        
    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision:  2015_04_20
================================================================================
*/
#ifndef CutList_h
#define CutList_h

#include <iostream>
#include <vector>
#include <fstream>
#include <cstdlib>
#include <string>

#include "../Cut/Cut.h"
#include "../../Libraries/Folder_List.h"

using namespace std;

class CutList
{
    public:
        CutList();
        CutList(int nMode);
        ~CutList();
        
        void writeCutTable();
      
        // -------------------------------------------------------------------------
        //     Cut Numbers
        //--------------------------------------------------------------------------
        // Common Cut Numbers
        Cut nCut_All;
        Cut nCut_Vertex_None;
        Cut nCut_Vertex_Not_Reconstructable; 
        Cut nCut_Vertex_Not_Fiducial;
        Cut nCut_Vertex_Count;
        Cut nCut_nProngs;
        
        // nProngs == 1 Cut Numbers
        Cut nCut_1Prong_Muon_None;              
        Cut nCut_1Prong_Muon_Not_Plausible;
        Cut nCut_1Prong_Muon_Score_Low;
        Cut nCut_1Prong_Muon_Charge;
        Cut nCut_1Prong_Vertex_Michel_Exist; 
        Cut nCut_1Prong_EndPoint_Michel_Exist;
        Cut nCut_1Prong_secEndPoint_Michel_Exist;
        Cut nCut_1Prong_PreFilter_Pi0;
        Cut nCut_1Prong_VtxBlob;
        Cut nCut_1Prong_ConeBlobs;
        Cut nCut_1Prong_Pi0_invMass;
        Cut nCut_1Prong_Photon1DistanceLow;
        Cut nCut_1Prong_Photon2DistanceLow;
        Cut nCut_1Prong_beamEnergy;
        Cut nCut_1Prong_UnusedE;
        
        // nProngs == 2 Cut Numbers
        Cut nCut_2Prong_Muon_None;              
        Cut nCut_2Prong_Muon_Not_Plausible;
        Cut nCut_2Prong_Muon_Score_Low;
        Cut nCut_2Prong_Muon_Charge;
        Cut nCut_2Prong_Vertex_Michel_Exist; 
        Cut nCut_2Prong_EndPoint_Michel_Exist;
        Cut nCut_2Prong_secEndPoint_Michel_Exist;
        Cut nCut_2Prong_PreFilter_Pi0;
        Cut nCut_2Prong_VtxBlob;
        Cut nCut_2Prong_ConeBlobs;
        Cut nCut_2Prong_Pi0_invMass;
        Cut nCut_2Prong_Photon1DistanceLow;
        Cut nCut_2Prong_Photon2DistanceLow;
        Cut nCut_2Prong_beamEnergy;
        Cut nCut_2Prong_UnusedE;
        Cut nCut_2Prong_Particle_None;
        Cut nCut_2Prong_Proton_None;            
        Cut nCut_2Prong_ProtonScore;
        Cut nCut_2Prong_DeltaInvMass;
        
    private:
        void SetAnalysisMode(int nMode);
        void SetCutNames();
        void OpenTextFiles();
        void formCutVectors();
        void writeCutTableHeader();
        void writeCutTableRows(vector<Cut> nCutVector, int nProngs, bool isShortList);
        double getCutEfficiency(double nSig, double effBase);
        double getCutPurity(double nSig, double nEvents);
        
        static const int nTopologies = 2;
        
        vector<Cut> nCutVector_1Prong;
        vector<Cut> nCutVector_2Prong;
        vector<Cut> nCutVector_1Prong_ShortList;
        vector<Cut> nCutVector_2Prong_ShortList;
        
        ofstream cutText[nTopologies];
        string cutFile[nTopologies];
        string branchDir;
        
};




#endif
