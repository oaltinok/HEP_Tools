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
#include "TLatex.h"
#include "TLine.h"
#include "TMinuit.h"

using namespace PlotUtils;

static const int nModels = 6;

struct SideBand
{
    std::string name;
    std::string model_names[nModels];  

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
        
        void ApplyFitResults(double chisq, double w_WithPi0, double w_QELike, double w_SinglePiPlus);
        void Plot();
        
        SideBand Michel;
        SideBand pID;
        SideBand LowInvMass;

    private:
        double ChiSq;
        double wgt_WithPi0;
        double wgt_QELike;
        double wgt_SinglePiPlus;

        void OpenRootFiles();
        void initSideBands();
        void SetNames(SideBand &sb, std::string name);
        void GetTH1D(TFile* f, TH1D* &h, std::string var_name);
        void ApplyFitResults();
        void ApplyFitResults(SideBand &sb);

        // Plot Functions
        void Plot(SideBand &sb);
        void Plot(SideBand &sb, int ind, bool isArea);
        void ColorHists(SideBand &sb);
        double GetMCScaleRatio(SideBand &sb, bool isArea);
};


#endif

