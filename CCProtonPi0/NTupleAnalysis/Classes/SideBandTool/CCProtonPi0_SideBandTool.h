#ifndef CCProtonPi0_SideBandTool_h
#define CCProtonPi0_SideBandTool_h

#include <iomanip>
#include "../NTupleAnalysis/CCProtonPi0_NTupleAnalysis.h"
#include "TObjArray.h"
#include "TFractionFitter.h"

using namespace PlotUtils;

static const int nModels = 6;

struct SideBand
{
    std::string name;
    std::string model_names[nModels];  

    // 0 Ratio from MC 
    // 1 Ratio from Fit
    // 2 Uncertainity
    double fr_signal[3];
    double fr_WithPi0[3];
    double fr_QELike[3];
    double fr_SinglePiPlus[3];
    double fr_Other[3];
    double fr_Total[3];

    TH1D* data;
    TH1D* signal;
    TH1D* WithPi0;
    TH1D* QELike;
    TH1D* SinglePiPlus;
    TH1D* Other;

    TFile* f_mc;
    TFile* f_data;
};

class CCProtonPi0_SideBandTool : public CCProtonPi0_NTupleAnalysis
{
    public:
        CCProtonPi0_SideBandTool();
        ~CCProtonPi0_SideBandTool();
        void Fit();

    private:
        SideBand Michel;
        SideBand pID;
        SideBand LowInvMass;

        void OpenRootFiles();
        void initSideBands();
        void SetNames(SideBand &sb, std::string name);
        void CalcMCRatios(SideBand &sb, bool isLimited = false, int first_bin = 1, int last_bin = 1);
        void GetTH1D(TFile* f, TH1D* &h, std::string var_name);
        void Fit(SideBand &sb, bool isLimitedFit = false, int first_bin = 1, int last_bin = 1);
        void PrintFitResults(SideBand &sb);
        void PrintRatio(double ratio[]);
        void OpenTextFile();
        
        std::string fileName;
        ofstream textFile;
};


#endif

