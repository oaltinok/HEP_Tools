#ifndef CCProtonPi0_NTupleAnalysis_h
#define CCProtonPi0_NTupleAnalysis_h
/*
   ================================================================================
Class: CCProtonPi0_NTupleAnalysis
Base Class for the CCProtonPi0_NTupleAnalysis Package
Includes Common Member Variables for all Classes

Author:        Ozgur Altinok  - ozgur.altinok@tufts.edu
================================================================================
*/
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>
#include <unordered_map>

#include "Cintex/Cintex.h"
#include <TROOT.h>
#include <TMath.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TFile.h>
#include <PlotUtils/MnvH1D.h>
#include <PlotUtils/MnvH2D.h>
#include <PlotUtils/FluxReweighter.h>

// Libraries
#include "../../Libraries/PDG_List.h"
#include "../../Libraries/Folder_List.h"
#include "../../Libraries/Data_Functions.h"
#include "../../Libraries/HEP_Functions.h"

using namespace PlotUtils;

class CCProtonPi0_NTupleAnalysis
{
    public:
        CCProtonPi0_NTupleAnalysis();

        std::string ana_folder;
        // Constants 
        static const int nFSIType = 6;
        static const int nIntType = 3;
        static const int nHistograms = 11;
        static const int nTopologies = 2;
        static const int n_universe = 100;
        static const int n_lateral_universes = 500;

        static const std::string version;
        static const double EPSILON;
        static const double data_POT;
        static const double mc_POT;
        static const double mc_2p2h_POT;
        static const double POT_ratio;
        static const double POT_ratio_2p2h;
        static const double SENTINEL;
        static const double MeV_to_GeV; 
        static const double MeVSq_to_GeVSq;
        static const double mSq_to_cmSq; 
        static const double mm_to_cm;
        static const double rad_to_deg;
        static const double max_muon_theta;
        static const double min_Enu;
        static const double max_Enu;
        static const double max_W;
        static const double max_muon_angle;
        static const double muon_mass;
        static const double pi0_mass;
        static const double piplus_mass;
        static const double proton_mass;
        static const double neutron_mass;
        static const double beam_theta;
        static const double beam_phi;

        // GENIE Tuning
        static const double genieMaRes;
        static const double genieMaRes1sig;
        // GENIE central value MvRES from electroproduction data fit
        static const double genieMvRes;
        static const double genieMvRes1sig;
        // Reduced MvRES error from electroproduction data fit
        static const double electroProdMvRes1sig;
        // Pion production parameters and errors from deuterium fit
        static const double deuteriumMaRes;
        static const double deuteriumMaRes1sig;
        static const double deuteriumNonResNorm;
        static const double deuteriumNonResNorm1sig;
        static const double deuteriumResNorm;
        static const double deuteriumResNorm1sig;
        // Delta Suppression 
        static const double DeltaFactor_A;
        static const double DeltaFactor_Q0;

        // 2p2h Events
        std::vector<double> fit_2p2h_CV;
        std::vector<double> fit_2p2h_np;
        std::vector<double> fit_2p2h_nn;
        void init2p2hFitResults();
        double Get_2p2h_wgt(Double_t* neutrino_4P, Double_t* muon_4P, std::vector<double> fit_results);
        
        // Event Kinematics 
        bool IsWInRange(double W);
        bool IsEnuInRange(double Enu);

        bool IsEvent2p2h(int type);

        double Calc_Enu_Truth(double muon_E, double proton_E, double pi0_E);
        double Calc_QSq(double Enu, double muon_E, double muon_P, double muon_angle_beam); 
        double Calc_WSq(double Enu, double QSq, double muon_E);

        // Flux Reweighter
        static const bool applyNuEConstraint;
        static const enum FluxReweighter::EFluxVersion new_flux;
        static const enum FluxReweighter::EG4NumiVersion old_flux;

        FluxReweighter* frw;
        bool frw_DefaultInit;
        bool processed_minerva1;
        bool processed_minerva7;
        bool processed_minerva9;
        bool processed_minerva13A;
        bool processed_minerva13B;
        bool processed_minerva13C;
        bool processed_minerva13D;
        bool processed_minerva13E;
        bool processed_2p2h;

        void UpdateFluxReweighter(const int run, int type);
        void ReInitFluxReweighter(enum FluxReweighter::EPlaylist playlist);
        double GetFluxWeight(double Enu, int nuPDG);
        std::vector<double> GetFluxError(double Enu, int nuPDG);

        void OpenTextFile(std::string file_name, std::ofstream &file);
        std::string GetPlaylist(const int run, int type);
        void printBins(const TH1* hist, const std::string var_name, bool useLowEdge = false);
        void printBins(const TH2* hist, const std::string var_name);
        void RunTimeError(std::string message);
        double Average_1DHist(const TH1* hist);

        // Returns a "new" MnvH1D or MnvH2D 
        MnvH1D* GetMnvH1D(TFile* f, std::string var_name);
        MnvH2D* GetMnvH2D(TFile* f, std::string var_name);

        // --------------------------------------------------------------------
        // Errors for Data
        // --------------------------------------------------------------------
        // Vertical Error Bands
        template<class MnvHistoType>
            void AddVertErrorBands_Data(MnvHistoType* h);

        template<class MnvHistoType>
            void AddLeadingErrorBands_Data(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBands_TruthTree(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBandsAndFillWithCV_TruthTree(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBands_FluxHistogram(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBandAndFillWithCV_Flux(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBandAndFillWithCV_Genie(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBandAndFillWithCV_BckgConstraint_WithPi0(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBandAndFillWithCV_BckgConstraint_SinglePiPlus(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBandAndFillWithCV_BckgConstraint_QELike(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBandAndFillWithCV_MichelTrue(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBandAndFillWithCV_MichelFake(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBandAndFillWithCV_TargetMass(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBandAndFillWithCV_2p2h(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBandAndFillWithCV_Unfolding(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBandAndFillWithCV_ProtonTracking(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBandAndFillWithCV_MuonTracking(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBandAndFillWithCV_NeutronResponse(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBandAndFillWithCV_PionResponse(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBandAndFillWithCV_HighMaRES(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBandAndFillWithCV_LowMaRES(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBandAndFillWithCV_DeltaFactor(MnvHistoType* h);

        // Lateral Error Bands
        template<class MnvHistoType>
            void AddLatErrorBands_Data(MnvHistoType* h);

        template<class MnvHistoType>
            void AddLatErrorBandsAndFillWithCV_TruthTree(MnvHistoType* h);

        template<class MnvHistoType>
            void AddLatErrorBands_FluxHistogram(MnvHistoType* h);

        template<class MnvHistoType>
            void AddLatErrorBandAndFillWithCV_ProtonEnergy_MassModel(MnvHistoType* h);

        template<class MnvHistoType>
            void AddLatErrorBandAndFillWithCV_ProtonEnergy_MEU(MnvHistoType* h);

        template<class MnvHistoType>
            void AddLatErrorBandAndFillWithCV_ProtonEnergy_BetheBloch(MnvHistoType* h);

        template<class MnvHistoType>
            void AddLatErrorBandAndFillWithCV_ProtonEnergy_Birks(MnvHistoType* h);

        template<class MnvHistoType>
            void AddLatErrorBandAndFillWithCV_MuonMomentum(MnvHistoType* h);

        template<class MnvHistoType>
            void AddLatErrorBandAndFillWithCV_MuonTheta(MnvHistoType* h);

        template<class MnvHistoType>
            void AddLatErrorBandAndFillWithCV_EM_EnergyScale(MnvHistoType* h);

        // --------------------------------------------------------------------
        // Errors for MC
        // --------------------------------------------------------------------
        // Vertical Error Bands
        template<class MnvHistoType>
            void AddVertErrorBands_MC(MnvHistoType* h);

        template<class MnvHistoType>
            void AddLeadingErrorBands_MC(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBand_Flux(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBand_Genie(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBand_BckgConstraint_WithPi0(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBand_BckgConstraint_SinglePiPlus(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBand_BckgConstraint_QELike(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBand_MichelTrue(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBand_MichelFake(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBand_TargetMass(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBand_2p2h(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBand_Unfolding(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBand_ProtonTracking(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBand_MuonTracking(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBand_NeutronResponse(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBand_PionResponse(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBand_HighMaRES(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBand_LowMaRES(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBand_DeltaFactor(MnvHistoType* h);

        // Lateral Error Bands
        template<class MnvHistoType>
            void AddLatErrorBands_MC(MnvHistoType* h);

        template<class MnvHistoType>
            void AddLatErrorBand_ProtonEnergy_MassModel(MnvHistoType* h);

        template<class MnvHistoType>
            void AddLatErrorBand_ProtonEnergy_MEU(MnvHistoType* h);

        template<class MnvHistoType>
            void AddLatErrorBand_ProtonEnergy_BetheBloch(MnvHistoType* h);

        template<class MnvHistoType>
            void AddLatErrorBand_ProtonEnergy_Birks(MnvHistoType* h);

        template<class MnvHistoType>
            void AddLatErrorBand_MuonMomentum(MnvHistoType* h);

        template<class MnvHistoType>
            void AddLatErrorBand_MuonTheta(MnvHistoType* h);

        template<class MnvHistoType>
            void AddLatErrorBand_EM_EnergyScale(MnvHistoType* h);

        // --------------------------------------------------------------------
        // Systematics - Multi Universe 
        // --------------------------------------------------------------------
        void GetPointersAllUniverses(MnvH1D* mnvh1d_hist, std::vector<TH1D*> &all_universes);
        void GetAllUniverses(MnvH1D* mnvh1d_hist, std::vector<TH1D*> &all_universes, std::vector<std::string> &err_bands, std::vector<int> &hist_ind);
        void ClearAllUniversesVector(std::vector<TH1D*> &all_universes);

    private:
};


#endif
