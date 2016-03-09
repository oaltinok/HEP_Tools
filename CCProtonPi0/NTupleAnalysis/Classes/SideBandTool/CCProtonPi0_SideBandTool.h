#ifndef CCProtonPi0_SideBandTool_h
#define CCProtonPi0_SideBandTool_h

#include <iomanip>
#include "../NTupleAnalysis/CCProtonPi0_NTupleAnalysis.h"
#include "TObjArray.h"
#include "TFractionFitter.h"
#include "TCanvas.h"
#include "TPad.h"
#include "THStack.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"

using namespace PlotUtils;

static const int nModels = 6;

struct SideBand
{
    std::string name;
    std::string model_names[nModels];  

    // 0 Ratio from MC 
    // 1 Ratio from Fit
    // 2 MC/Fit
    // 3 Error 
    double fr_signal[4];
    double fr_WithPi0[4];
    double fr_QELike[4];
    double fr_SinglePiPlus[4];
    double fr_Other[4];
    double fr_Total[4];

    TH1D* data; 
    TH1D* mc_total; 
    TH1D* fit; // Hist for fit

    // 2 Hists for MC Models
    //      ind = 0 for original
    //      ind = 1 for modified according to fit
    TH1D* signal[2];
    TH1D* WithPi0[2];
    TH1D* QELike[2];
    TH1D* SinglePiPlus[2];
    TH1D* Other[2];

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
        double POT_ratio;

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
        void ApplyFitResults(SideBand &sb);

        // Plot Functions
        void ColorHists(SideBand &sb);
        double GetMCScaleRatio(SideBand &sb, bool isArea);
        void Plot(SideBand &sb, int ind, bool isArea);
        std::string fileName;
        ofstream textFile;
};


#endif

