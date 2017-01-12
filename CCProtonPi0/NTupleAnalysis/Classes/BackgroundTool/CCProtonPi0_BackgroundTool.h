/*
================================================================================
Class: CCProtonPi0_BackgroundTool
    Class Responsible for 
        Calculating Statistics for Background Events
        Creating Background Tables
    
    Main Directory:
        Classes/BackgroundTool
 
    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
================================================================================
*/
#ifndef CCProtonPi0_BackgroundTool_h
#define CCProtonPi0_BackgroundTool_h

#include "../NTupleAnalysis/CCProtonPi0_NTupleAnalysis.h"

using namespace std;

struct Background
{ 
    string name; 
    double nAll;
    double nWithMichel;
};

class CCProtonPi0_BackgroundTool : public CCProtonPi0_NTupleAnalysis
{
    public:
        CCProtonPi0_BackgroundTool(bool isModeReduce);
        ~CCProtonPi0_BackgroundTool();
        
        void set_nTracks(int input);
        void writeBackgroundTable();
        void fillBackgroundWithPi0(bool NoPi0, bool SinglePi0, bool MultiPi0, bool withMichel, double wgt = 1.0);
        void fillBackgroundCompact(bool WithPi0, bool QELike, bool SinglePiPlus, bool Other, double wgt = 1.0);
        
        void fillBackground(bool NC,
                            bool AntiNeutrino,
                            bool QELike,
                            bool SingleChargedPion,
                            bool SingleChargedPion_ChargeExchanged,
                            bool DoublePionWithPi0,
                            bool DoublePionWithoutPi0,
                            bool MultiPionWithPi0,
                            bool MultiPionWithoutPi0,
                            bool Other,
                            bool withMichel,
                            double wgt = 1.0);
    private:
        static const int nTables = 3;
        int nTracks;
        string fileName[nTables];
        ofstream textFile[nTables];
        
        vector< vector<Background> > BackgroundWithPi0Vector;
        vector< vector<Background> > BackgroundCompactVector;
        vector< vector<Background> > BackgroundTypeVector;
       
        // Background with Pi0
        vector<Background> bckg_NoPi0;
        vector<Background> bckg_SinglePi0;
        vector<Background> bckg_MultiPi0;
        vector<Background> bckg_Total_WithPi0;
       
        // Background Types Compact
        vector<Background> bckg_compact_WithPi0;
        vector<Background> bckg_compact_QELike;
        vector<Background> bckg_compact_SinglePiPlus;
        vector<Background> bckg_compact_Other;
        vector<Background> bckg_compact_Total;

        // Background Types 
        vector<Background> bckg_NC;
        vector<Background> bckg_AntiNeutrino;
        vector<Background> bckg_QELike;
        vector<Background> bckg_SingleChargedPion;
        vector<Background> bckg_SingleChargedPion_ChargeExchanged;
        vector<Background> bckg_DoublePionWithPi0;
        vector<Background> bckg_DoublePionWithoutPi0;
        vector<Background> bckg_MultiPionWithPi0;
        vector<Background> bckg_MultiPionWithoutPi0;
        vector<Background> bckg_Other;
        vector<Background> bckg_Total;

        void updateBackground(Background &b, bool withMichel = false, double wgt = 1.0);
        double calcPercent(double nEvents, double nBase);
        void formBackgroundVectors();
        void initBackgrounds();
        void initSingleBackground(vector<Background> &b, string name);
        void writeBackgroundTableHeader();
        void writeBackgroundTableRows(vector< vector<Background> > &bckgVector);
        void OpenTextFiles();

};


#endif

