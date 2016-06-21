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

struct XSec_Var
{
    std::string name;

    MnvH1D* data; 
    MnvH1D* mc_total; 

    // 2 Hists for MC Models
    //      ind = 0 for original
    //      ind = 1 for modified according to fit
    MnvH1D* signal[2];
    MnvH1D* WithPi0[2];
    MnvH1D* QELike[2];
    MnvH1D* SinglePiPlus[2];
    MnvH1D* Other[2];
};

struct SideBand
{
    std::string name;
    std::string model_names[nModels];  

    // Fit Variable
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

    std::vector<TH1D*> data_all_universes;
    std::vector<TH1D*> mc_total_all_universes;
    std::vector<TH1D*> signal_all_universes;
    std::vector<TH1D*> WithPi0_all_universes;
    std::vector<TH1D*> QELike_all_universes;
    std::vector<TH1D*> SinglePiPlus_all_universes;
    std::vector<TH1D*> Other_all_universes;

    XSec_Var muon_P;
    XSec_Var muon_theta;
    XSec_Var pi0_P;
    XSec_Var pi0_KE;
    XSec_Var pi0_theta;
    XSec_Var neutrino_E;
    XSec_Var QSq;
    XSec_Var W;
    
    // File for Fit Variable
    TFile* f_mc_fit;
    TFile* f_data_fit;

    // File for XSec Variable
    TFile* f_mc_var;
    TFile* f_data_var;

};


class CCProtonPi0_SideBandTool : public CCProtonPi0_NTupleAnalysis
{
    public:
        CCProtonPi0_SideBandTool();
        ~CCProtonPi0_SideBandTool();
        
        void SaveFitResults(double chisq, double par_values[3], double par_errors[3]);
        void ApplyFitResults();
        void Plot();
        
        SideBand Original;
        SideBand Michel;
        SideBand pID;
        SideBand LowInvMass;
        SideBand HighInvMass;
        
        int N_Universes;
        int current_unv;

        std::vector<double> ChiSq;
        std::vector<double> wgt_WithPi0;
        std::vector<double> wgt_QELike;
        std::vector<double> wgt_SinglePiPlus;

        std::vector<double> err_WithPi0;
        std::vector<double> err_QELike;
        std::vector<double> err_SinglePiPlus;

    private:
        void OpenRootFiles();
        void initSideBands();
        void initSideBand_FitHistograms(SideBand &sb);
        void initSideBand_AllUniverses(SideBand& sb);
        void initSideBand_XSecHistograms(SideBand& sb);
        void initSideBand_XSecHistograms(SideBand& sb, XSec_Var& xsec_var, std::string name);
        void SetNames(SideBand &sb, std::string name);
        void ApplyFitResults(SideBand &sb);
        void ApplyFitResults(XSec_Var &xsec_var);

        // Plot Functions
        void DrawDataMCWithErrorBand(SideBand &sb);
        void Plot(SideBand &sb);
        void Plot(SideBand &sb, int ind);
        void Plot(SideBand &sb, XSec_Var &xsec_var, int ind, std::string var_name);
        void Plot(int ind, std::string sb_name, std::string var_name, MnvH1D* data, MnvH1D* mc_total, MnvH1D* signal, MnvH1D* WithPi0, MnvH1D* QELike, MnvH1D* SinglePiPlus, MnvH1D* Other);
        void ColorHists(MnvH1D* data, MnvH1D* signal, MnvH1D* WithPi0, MnvH1D* QELike, MnvH1D* SinglePiPlus, MnvH1D* Other);
        double calc_ChiSq(MnvH1D* data, MnvH1D* signal, MnvH1D* WithPi0, MnvH1D* QELike, MnvH1D* SinglePiPlus, MnvH1D* Other);

        std::string var_name;
};


#endif

