#ifndef CCProtonPi0_SideBandTool_h
#define CCProtonPi0_SideBandTool_h

#include <iomanip>
#include "../NTupleAnalysis/CCProtonPi0_NTupleAnalysis.h"
#include "TObjArray.h"
#include "TFractionFitter.h"
#include "TCanvas.h"
#include "TPad.h"
#include "THStack.h"
#include "TGaxis.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TLine.h"
#include "TMinuit.h"
#include <PlotUtils/MnvPlotter.h>

using namespace PlotUtils;

static const int nModels = 6;

struct SideBand
{
    std::string name;
    std::string model_names[nModels];  

    MnvH1D* data; 
    MnvH1D* mc_total; 
    TH1D* fit; // Hist for fit

    // 2 Hists for MC Models
    //      ind = 0 for original
    //      ind = 1 for modified according to fit
    MnvH1D* signal[2];
    MnvH1D* WithPi0[2];
    MnvH1D* QELike[2];
    MnvH1D* SinglePiPlus[2];
    MnvH1D* Other[2];

    TFile* f_mc;
    TFile* f_data;
};

class CCProtonPi0_SideBandTool : public CCProtonPi0_NTupleAnalysis
{
    public:
        CCProtonPi0_SideBandTool(std::string var);
        ~CCProtonPi0_SideBandTool();
        
        void ApplyFitResults(double chisq, double par_values[3], double par_errors[3]);
        void Plot();
        
        SideBand Original;
        SideBand Michel;
        SideBand pID;
        SideBand LowInvMass;
        SideBand HighInvMass;

    private:
        double ChiSq;
        double wgt_WithPi0;
        double wgt_QELike;
        double wgt_SinglePiPlus;

        double err_WithPi0;
        double err_QELike;
        double err_SinglePiPlus;

        void OpenRootFiles();
        void OpenRootFiles_Test();
        void initSideBands();
        void SetNames(SideBand &sb, std::string name);
        void GetMnvH1D(TFile* f, MnvH1D* &h, std::string var_name);
        void ApplyFitResults();
        void ApplyFitResults(SideBand &sb);

        // Plot Functions
        void DrawDataMCWithErrorBand(SideBand &sb);
        void Plot(SideBand &sb);
        void Plot(SideBand &sb, int ind, bool isArea);
        void ColorHists(SideBand &sb);
        double GetMCScaleRatio(SideBand &sb, bool isArea);
        double calc_ChiSq(SideBand &sb, int ind, bool isArea);

        std::string var_name;
};


#endif

