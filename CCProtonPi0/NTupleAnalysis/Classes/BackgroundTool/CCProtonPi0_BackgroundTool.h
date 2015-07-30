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
        
        void writeBackgroundTable();
        void fillBackgroundWithPi0( int nProngs, bool NoPi0, bool SinglePi0, bool MultiPi0, bool withMichel);
        
        void fillBackground(int nProngs,
                            bool NC,
                            bool AntiNeutrino,
                            bool QELike,
                            bool SinglePion,
                            bool DoublePion,
                            bool MultiPion,
                            bool Other,
                            bool withMichel);
    private:
        string fileName[nTopologies];
        ofstream textFile[nTopologies];
        
        vector< vector<Background> > BackgroundWithPi0Vector;
        vector< vector<Background> > BackgroundTypeVector;
       
        // Background with Pi0
        vector<Background> bckg_NoPi0;
        vector<Background> bckg_SinglePi0;
        vector<Background> bckg_MultiPi0;
        vector<Background> bckg_Total_WithPi0;
        
        // Background Types 
        vector<Background> bckg_NC;
        vector<Background> bckg_AntiNeutrino;
        vector<Background> bckg_QELike;
        vector<Background> bckg_SinglePion;
        vector<Background> bckg_DoublePion;
        vector<Background> bckg_MultiPion;
        vector<Background> bckg_Other;
        vector<Background> bckg_Total;

        void updateBackground(Background &b, bool withMichel);
        int setArrayInd(int nProngs);
        
        double calcPercent(double nEvents, double nBase);
        void formBackgroundVectors();
        void initBackgroundVectors();
        void initSingleBackgroundVector(vector<Background> &b, string name);
        void writeBackgroundTableHeader(int nProngs);
        void writeBackgroundTableRows(vector< vector<Background> > &bckgVector, int nProngs);
        void OpenTextFiles();

};


#endif

