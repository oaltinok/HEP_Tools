/*
================================================================================
Class: CCProtonPi0_BackgroundTool
    Class Responsible for 
        Calculating Statistics for Background Events
        Creating Background Tables
    
    Main Directory:
        Classes/BackgroundTool
 
    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision:  2015_05_07
================================================================================
*/
#ifndef CCProtonPi0_BackgroundTool_h
#define CCProtonPi0_BackgroundTool_h

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>

#include "../NTupleAnalysis/CCProtonPi0_NTupleAnalysis.h"

using namespace std;

class CCProtonPi0_BackgroundTool : public CCProtonPi0_NTupleAnalysis
{
    public:
        CCProtonPi0_BackgroundTool(int nMode);
        ~CCProtonPi0_BackgroundTool();
        
        void writeBackgroundTable();
        void fillBackgroundBranches(    int nProngs,
                                        bool truth_isBckg_QELike, 
                                        bool truth_isBckg_SinglePiPlus, 
                                        bool truth_isBckg_SinglePiMinus, 
                                        bool truth_isBckg_MultiPion, 
                                        bool truth_isBckg_MultiPiZero, 
                                        bool truth_isBckg_Other,
                                        bool truth_isBckg_withAntiMuon,
                                        bool truth_isBckg_withMichel,
                                        bool truth_isBckg_withPrimaryPi0,
                                        bool truth_isBckg_withSecondaryPi0);

    private:
        static const int nBckgBranch = 6;
        
        string backgroundFile[nTopologies];
        ofstream backgroundText[nTopologies];
        
        vector< vector<double> > bckgVector_1Prong;
        vector< vector<double> > bckgVector_2Prong;
        
        vector<double> bckg_1Prong_QELike;
        vector<double> bckg_1Prong_SinglePiPlus;
        vector<double> bckg_1Prong_SinglePiMinus;
        vector<double> bckg_1Prong_MultiPion;
        vector<double> bckg_1Prong_MultiPiZero;
        vector<double> bckg_1Prong_Other;
        
        vector<double> bckg_2Prong_QELike;
        vector<double> bckg_2Prong_SinglePiPlus;
        vector<double> bckg_2Prong_SinglePiMinus;
        vector<double> bckg_2Prong_MultiPion;
        vector<double> bckg_2Prong_MultiPiZero;
        vector<double> bckg_2Prong_Other;
        
        void setBackgroundBranch(vector<double>& background, bool hasAntiMuon, bool hasMichel, bool hasPrimaryPi0, bool hasSecondaryPi0 );
        void formBackgroundVectors();
        void writeBackgroundTableHeader();
        void writeBackgroundTableRows(vector< vector<double> > bckgVector, int nProngs);
        void OpenTextFiles();

};


#endif

