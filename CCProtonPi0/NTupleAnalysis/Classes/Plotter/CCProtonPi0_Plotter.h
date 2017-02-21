/*
   ================================================================================
Class: CCProtonPi0_Plotter
CCProtonPi0_Plotter class includes specific functions for plotting 1D or 2D histograms

Main Directory:
Classes/Plotter/

Usage:
> #include "Classes/Plotter/CCProtonPi0_Plotter.cpp" 
> CCProtonPi0_Plotter p;
> p.plotHistograms(std::string mcFile, std::string plotDir)


Author:        Ozgur Altinok  - ozgur.altinok@tufts.edu
================================================================================
*/

#ifndef CCProtonPi0_Plotter_H
#define CCProtonPi0_Plotter_H
#include <numeric>

#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TColor.h>
#include <TLine.h>
#include <TPad.h>
#include <TPaveStats.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TVectorD.h>
#include <TVirtualFitter.h>
#include <PlotUtils/MnvPlotter.h>
#include <PlotUtils/MnvFluxConstraint.h>
#include <PlotUtils/POTCounter.h>
#include <PlotUtils/TargetUtils.h>
#include <MinervaUnfold/MnvUnfold.h>
#include "Cintex/Cintex.h"

#include "../NTupleAnalysis/CCProtonPi0_NTupleAnalysis.h"
#include "../BinList/CCProtonPi0_BinList.h"
#include "../QSqFitter/CCProtonPi0_QSqFitter.h"

using namespace PlotUtils;

struct rootDir
{
    std::string data;
    std::string mc;
};

class CutArrow
{
    public:
        CutArrow(){ /* Do Nothing */ }
        CutArrow(double cut, std::string arrow_d )
        {
            cut_location = cut;
            arrow_direction = arrow_d;
        }

        double cut_location;
        std::string arrow_direction;
};

class CCProtonPi0_Plotter : public CCProtonPi0_NTupleAnalysis
{
    public:
        CCProtonPi0_Plotter();
        void plotHistograms();

    private:
        bool thesisStyle;

        rootDir rootDir_PC;
        rootDir rootDir_GENIEXSec;
        rootDir rootDir_Truth;
        rootDir rootDir_CrossSection;
        rootDir rootDir_CutHists;
        rootDir rootDir_Interaction;
        rootDir rootDir_Muon;
        rootDir rootDir_Proton;
        rootDir rootDir_Pion;

        // Error Summary Groups
        std::vector<std::string> detGroup;
        std::vector<std::string> genieGroup;
        std::vector<std::string> fsiGroup;
        std::vector<std::string> fluxGroup;
        std::vector<std::string> otherGroup;
        
        // POT Stats
        double POT_Ratio_data_mc;

        void setRootDirs();
        void getPOT_MC();
        void getPOT_Data();

        // Cross Section Plots
        bool plot_muon_P;
        bool plot_muon_theta;
        bool plot_pi0_P;
        bool plot_pi0_KE;
        bool plot_pi0_theta;
        bool plot_QSq;
        bool plot_W;
        bool plot_Enu;

        void PlotXSecVar(std::string var_name, std::string data_var, std::string mc_var, std::string plotDir, std::string plotName);
        void PlotXSecVar_WithMiniBoone(std::string var_name, std::string plotDir);
        void PlotXSecVar_BeforeFSI(std::string var_name, std::string plotDir);
        void PlotXSecVar_FSIType(std::string var_name, std::string plotDir);
        void PlotXSecVar_IntType(std::string var_name, std::string plotDir);
        
        void plotCrossSection();
        void plotOriginalData();
        void plotBackgroundSubtracted();
        void plotBackgroundEstimated();
        void plotEfficiencyCorrected();
        void plotUnfolded();
        void plotFluxIntegrated();
        void plotXSec();
        void plotCrossSection_Check();
        void plotCrossSection_Check(std::string var_name, std::string plotDir);
        void CheckAllUniverses(std::string test_name, MnvH1D* data, MnvH1D* mc);
        void CompareUniversesBinByBin(const std::vector<TH1D*> data_hists, const std::vector<TH1D*> mc_hists, std::string err_name, std::string test_name);

        // Data vs MC
        void plotInteraction_DataMC();
        void plotMuon_DataMC();
        void plotProton_DataMC();
        void plotPion_DataMC();
        void plotCutHistograms_DataMC();

        // MC Only
        void plotTruth();
        void plotInteraction_MCOnly();
        void plotMuon_MCOnly();
        void plotProton_MCOnly();
        void plotPion_MCOnly();
        void plotCutHistograms_MCOnly();

        // True Signal Events
        void plotTruth_Pion();
        void plotTruth_Enu();
        void plotTruth_QSq();
        void plotTruth_W();
        void plotTruth_ShortProton();

        // Other Plots 
        void BckgSubtraction_Studies();
        void QSq_Studies();
        void Studies_2p2h();
        void Draw_QSq_EnuFit(std::string data_dir, std::string mc_dir, std::string var_name, double* pars);
        void Draw_QSq_MaRES_Fit(bool isAreaNorm);
        void Draw_QSq_MaRES_Fit_SB();
        void Draw_QSq_MaRES_Plots();
        void Draw_QSq_MaRES_AreaNorm();
        void Draw_QSq_EnuLimit();
        void Draw_QSq_DeltaSuppression();
        void Draw_QSq_DeltaSuppression_v2(std::string var_name);
        void Draw_QSq_DeltaSuppression_AllPlots();
        void plot_CV_weights();
        void plot_SystematicsInfo();
        void PlotTotalEnuXSec();
        void plotGENIEXSec();
        void plotOtherStudies();
        void plotPC_MINOS_Pi0();
        void plotPC_MINOS_Pi0(std::string rootDir_Steel, std::string rootDir_Carbon, std::string var, std::string label);
        void plot_InvMass_TruthMatch_Stacked(bool isSignal, bool isStacked);
        void plot_Michel_TruthMatch(std::string var);
        void plot_SignalKinematics(std::string var, std::string type, bool isStacked);
        void plot_SignalKinematics();
        void plot_stacked_pi0_P();
        void plot_stacked_pi0_theta();
        void plotStandardHistograms(rootDir &dir, std::string plotDir);

        // Helper Functions
        double Calc_Normalized_NBackground(std::string var_name);
        double Get2DTotalFlow(MnvH2D* h);
        void NormalizeHistogram(MnvH2D* h);
        void NormalizeHistogram(MnvH1D* h);
        void NormalizeHistogram(TH1D* h);
        void ApplyStyle(MnvPlotter* plotter);
        void ApplyStyle_Legend(TLegend* legend);
        void ApplyStyle_Errors(MnvPlotter* plotter, bool groupErrors);
        void ApplyStyle_Thesis();
        void AddCutArrow(MnvPlotter* plotter, CutArrow &cutArrow, double hist_max, double arrow_length);
        void SavePi0InvMassPoints();
        void NormalizeToNormBinWidth(MnvH1D* hist);
        double GetSmallestBinWidth(MnvH1D* hist);
        TH1D* GetBinNormalizedTH1D(MnvH1D* hist, bool WithSystError = false);

        // W Study
        double pars_MC_deltaRES[3];
        double pars_MC_otherRES[3];
        double pars_MC_nonRES_G1[3];
        double pars_MC_nonRES_G2[3];

        void W_Studies();
        void plot_W_FitMinuit(double wgt_DeltaRES, double wgt_OtherRES, double wgt_NonRES);
        void init_W_FitResults();
        void Get_Exponential(double* x, double* y, int nPoints, double a, double b);
        void Get_BreitWigner(double* pars, double* x, double* y, int nPoints);
        void Get_Gaussian(double* pars, double* x, double* y, int nPoints);
        void plot_W_FitResults();
        void W_Fit_Data(std::string fit_name, double* pars_deltaRES, double* pars_otherRES, double* pars_nonRES_G1, double* pars_nonRES_G2, int nPars);
        void W_Fit_MC();
        void W_Fit_MC_DeltaRES();
        void W_Fit_MC_OtherRES();
        void W_Fit_MC_NonRES();
        void printBins_W();
        double Calc_ChiSq_dof(double* data, double* expected, int nPoints, int nPars);
        double Calc_ChiSq(TH1* data, TH1* MC);
        double Calc_ChiSq(TH1* data, TH1* MC, int min_bin, int max_bin);
        void DeltaRes_Studies();
        void expo_fit(const std::string& fileName, const std::string& histName);
        double user_expo(double* x, double* par);

        // Flux Study
        double GetFluxHistContent(MnvH1D* hist, double low1, double low2);
        double Integrate(MnvH1D* hist, int start, int end);
        void GetFlux();
        void PlotFluxHistograms();
        void PlotFluxComparison(std::string plotDir);
        void PlotFluxRatio(std::string plotDir);
        void PlotFluxRebinned(std::string plotDir);
        void PlotFluxRebinned(MnvH1D* original, MnvH1D* rebinned, std::string plotDir);

        // GENIE Tuning Study
        void GENIE_Tuning_Study();
        void XSecVars_GENIE_Tuning_Ratios();
        void GENIE_Tuning_GetHistograms(std::string var_name, std::string data_var, std::string mc_var, TH1D* &data_nominal, TH1D* &data_tuned_v2, TH1D* &data_tuned_v3, TH1D* &mc_nominal, TH1D* &mc_tuned_v2, TH1D* &mc_tuned_v3);
        void GENIE_Tuning_Ratio(std::string var_name, std::string data_var, std::string mc_var);
        void GENIE_Tuning_DataMC_Ratio(std::string var_name, std::string data_var, std::string mc_var);
        void Draw_Comparison_Nominal();
        void Draw_Comparison_Nominal(std::string var_name);
        void Draw_Comparison_DeltaFactor();
        void Draw_Comparison_DeltaFactor(std::string var_name);

        // --------------------------------------------------------------------
        // Plottting Macros - Implemented in CCProtonPi0_Plotter_Macros.cpp
        // --------------------------------------------------------------------
        // Default
        void DrawMnvH1D(rootDir& dir, std::string var_name, std::string plotDir);
        void DrawMnvH1D(MnvH1D* hist1D, std::string var_name, std::string plotDir);
        void DrawMnvH2D(std::string root_dir, std::string var_name, std::string plotDir);
        void DrawMnvH2D(MnvH2D* hist2D, std::string var_name, std::string plotDir, bool isMC);
        void DrawMnvH2D_Signal(rootDir root_dir, std::string var_name, std::string plotDir, double nBckg, bool isMC);
        void Draw1DHist(TH1* hist1D, std::string var_name, std::string plotDir, bool isLogScale = false);

        void Draw1DHist(rootDir &dir, std::string var_name, std::string plotDir, bool isLogScale = false);
        void Draw1DHist_Threshold(rootDir &dir, std::string var_name, std::string plotDir, double threshold = 0, bool isLogScale = false);
        void Draw2DHist(rootDir& dir, std::string var_name, std::string plotDir, double threshold = 0);
        void Draw3DHist(rootDir& dir, std::string var_name, std::string plotDir);

        // MC Only
        void DrawNormalizedMigrationHistogram(rootDir &dir, std::string var_name, std::string plotDir);
        void DrawMCWithErrorBand(rootDir& dir, std::string var_name, std::string plotDir);
        void DrawSignalMC(rootDir &dir, std::string var_name, std::string plotDir, int nCutArrows = 0, CutArrow cutArrow1 = CutArrow(), CutArrow cutArrow2 = CutArrow());
        void DrawErrorBand(MnvH1D* hist, std::string error_name, int error_color, std::string var_name, std::string plotDir);
        void DrawStackedMC(rootDir &dir, std::string var_name, std::string plotDir, int nCutArrows = 0, CutArrow cutArrow1 = CutArrow(), CutArrow cutArrow2 = CutArrow());
        void DrawStackedMC_BckgAll(rootDir &dir, std::string var_name, std::string plotDir, int nCutArrows = 0, CutArrow cutArrow1 = CutArrow(), CutArrow cutArrow2 = CutArrow());
        void DrawStackedMC_BckgType(rootDir &dir, std::string var_name, std::string plotDir, int nCutArrows = 0, CutArrow cutArrow1 = CutArrow(), CutArrow cutArrow2 = CutArrow());

        // Data vs MC
        void DrawDataMC(rootDir& dir, std::string var_name, std::string plotDir);
        void DrawDataMC_Thesis(rootDir& dir, std::string var_name, std::string plotDir);
        void DrawDataMC_Thesis(MnvH1D* data, MnvH1D* mc, std::string var_name, std::string plotDir, bool isXSec = false);
        void DrawDataMC(MnvH1D* data, MnvH1D* mc, std::string var_name, std::string plotDir, bool isXSec = false);
        void DrawDataMC_Signal(rootDir& dir, std::string var_name, std::string plotDir, double nBckg);
        void DrawDataMC_WithRatio(MnvH1D* data, MnvH1D* mc, std::string var_name, std::string plotDir, bool isPOTNorm, bool isXSec = false);
        void DrawDataMC_WithOtherData(MnvH1D* data, MnvH1D* mc, TGraph* otherData, std::string var_name, std::string ext_data_name, std::string plotDir);
        void DrawDataMC_BeforeFSI(MnvH1D* data, MnvH1D* mc, MnvH1D* mc_BeforeFSI, std::string var_name,  std::string plotDir);
        void DrawDataMC_PaperStyle(MnvH1D* data, MnvH1D* mc, MnvH1D* mc_BeforeFSI, std::string var_name,  std::string plotDir);
        void DrawDataMC_FSIType(MnvH1D* data, MnvH1D* mc, std::vector<MnvH1D*> mc_FSIType, std::string var_name,  std::string plotDir);
        void DrawDataMC_IntType(MnvH1D* data, MnvH1D* mc, std::vector<MnvH1D*> mc_IntType, std::string var_name,  std::string plotDir);
        void DrawDataMCRatio(MnvH1D* data, MnvH1D* mc, std::string var_name, std::string plotDir, bool isPOTNorm, bool isXSec = false);
        void DrawDataStackedMC(rootDir &dir, std::string var_name, std::string plotDir, int nCutArrows = 0, CutArrow cutArrow1 = CutArrow(), CutArrow cutArrow2 = CutArrow());
        void DrawDataStackedMC_BckgAll(rootDir &dir, std::string var_name, std::string plotDir, bool isPOTNorm, int nCutArrows = 0, CutArrow cutArrow1 = CutArrow(), CutArrow cutArrow2 = CutArrow());
        void DrawDataStackedMC_BckgType(rootDir &dir, std::string var_name, std::string plotDir, bool isPOTNorm, int nCutArrows = 0, CutArrow cutArrow1 = CutArrow(), CutArrow cutArrow2 = CutArrow());
        void DrawDataStackedMC_WithSignalTypes(rootDir &dir, std::string var_name, std::string plotDir);
        void DrawDataStackedMC_Signal(rootDir &dir, std::string var_name, std::string plotDir, double nBckg);
        void DrawDataMCSignal_Diff(rootDir& dir, std::string var_name, std::string plotDir, double nBckg);
        void DrawDataMCSignal_Diff_2D(rootDir& dir, std::string var_name, std::string plotDir, double nBckg);

        // Other
        MnvH1D* GetBckgSubtractedData(rootDir& dir, std::string var_name, double nBckg);
        MnvH2D* GetBckgSubtractedData_2D(rootDir& dir, std::string var_name, double nBckg);
        void PlotDelta();
        void DrawBackgroundSubtraction(bool isMC);
        void DrawErrorSummary(MnvH1D* hist, std::string var_name, std::string plotDir, bool groupErrors = true);
        void DrawErrorSummary_PaperStyle(MnvH1D* hist, std::string var_name, std::string plotDir, bool groupErrors = true);
        void DrawTGraph(rootDir &dir, std::string var_name, std::string plotDir);
        void DrawEfficiencyCurve(rootDir& dir, std::string var_name, std::string plotDir);
        void DrawStackedMC_GammaEvis(rootDir &dir, int gammaID, std::string plotDir);
        void DrawStackedMC_GammaByPDG(rootDir &dir, std::string var_name, int gammaID, std::string plotDir);

        // Helper
        double FormTObjArray_BckgType(TFile* f_mc, std::string var_name, TObjArray* mc_hists, double &hist_max, double &bin_width); 
        void AddNormBox(MnvPlotter* plotter, bool isPOTNorm, double mc_ratio);
        void SaveRecoRatioPoints(rootDir& dir, std::string var_name, std::string plotDir);
        void Save2DHistPoints(rootDir& dir, std::string var_name, std::string plotDir);     
        double GetMCNormalization(std::string &norm_label, bool isPOTNorm, MnvH1D* data, MnvH1D* mc);

        // Unfolding Study
        void UnfoldingStudy();
        void UnfoldingStudy_Iterations(std::string var_name);
        void UnfoldingStudy_muon_P();
        void UnfoldingStudy_muon_theta();
        void UnfoldingStudy_muon_cos_theta();
        void UnfoldingStudy_pi0_P();
        void UnfoldingStudy_pi0_KE();
        void UnfoldingStudy_pi0_theta();
        void init_UnfoldingHistograms(std::vector<MnvH1D*> &unfolded, std::vector<MnvH1D*> &error, std::vector<MnvH1D*> &diff);
        MnvH1D* CalcUnfoldingError(MnvH1D* diff, MnvH1D* truth);
        MnvH1D*  CalcUnfoldingDiff(MnvH1D* unfolded, MnvH1D* truth);
        void FillUnfoldingHistograms(MnvH1D* &unfolded, MnvH1D* &error, MnvH1D* &diff, MnvH2D* response, MnvH1D* mc_reco, MnvH1D* mc_true, int niter);
        void StyleUnfoldingHistograms(std::vector<MnvH1D*> &hists);
        void PlotUnfolding_StatErrors(const std::vector<MnvH1D*> &hists, const MnvH1D* truth, std::string var_name);
        void PlotUnfolding_Unfolded(const std::vector<MnvH1D*> &hists, const MnvH1D* truth, std::string var_name);
        void PlotUnfolding_Error(const std::vector<MnvH1D*> &hists, std::string var_name);
        void PlotUnfolding_Diff(const std::vector<MnvH1D*> &hists, std::string var_name);
        void PlotUnfolding_TruthComparison();
        void PlotUnfolding_Migration();

        // --------------------------------------------------------------------
        // Systematics
        //      Implementations are in CCProtonPi0_Plotter_Systematics.cpp
        // --------------------------------------------------------------------
        // Error Grouping
        bool IsErrorInGroup(std::string err_name, std::vector<std::string> errGroup);
        void Systematics_SetErrorSummaryGroups();
        void Clear_ErrorSummaryGroups();
        
        // Macros & Plotting
        void Systematics();
        void Systematics_CheckErrorSummary(std::string root_dir, std::string var_name);
        void Systematics_invMass();
        void Systematics_XSec();
        void Systematics_DrawErrorSummary(std::string data_var, std::string mc_var);
        void Systematics_DrawErrorBand_GENIE(std::string mc_var);
        void Systematics_DrawErrorSummary_GENIE(MnvH1D* hist, std::string var_name, std::string error_name, double err_genie_total);
        void Systematics_DrawErrorSummary_Group(std::string data_var, std::string mc_var, std::vector<std::string> errGroup, std::string group_name);
        // Tables
        void Systematics_WriteTables(std::string var_name);
        void Systematics_WriteTable_Fraction(MnvH1D* hist, std::string var_name);
        void Systematics_WriteTable_BinByBin(MnvH1D* hist, std::string var_name);
        TH1D* GetTotalErrorInGroup(MnvH1D* hist, std::vector<std::string> errGroup, bool area_normalized = false);

};

#endif

