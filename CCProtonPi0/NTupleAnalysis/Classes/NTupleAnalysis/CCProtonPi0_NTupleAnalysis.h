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
#include <fstream>
#include <cstdlib>
#include <string>
#include <vector>

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
        static const int nHistograms = 7;
        static const int nTopologies = 2;
        static const int n_universe = 100;
        static const int n_lateral_universes = 500;

        static const std::string version;
        static const double EPSILON;
        static const double data_POT;
        static const double mc_POT;
        static const double POT_ratio;
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

        // Event Kinematics 
        bool IsWLow(double true_W);
        bool IsEnuInRange(double Enu);

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
        bool processed_minerva13;

        void UpdateFluxReweighter(const int run);
        void ReInitFluxReweighter(enum FluxReweighter::EPlaylist playlist);
        double GetFluxWeight(double Enu, int nuPDG);
        std::vector<double> GetFluxError(double Enu, int nuPDG);

        void OpenTextFile(std::string file_name, std::ofstream &file);
        std::string GetPlaylist(const int run);
        void printBins(const TH1* hist, const std::string var_name, bool useLowEdge = true);
        void RunTimeError(std::string message);

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
            void AddVertErrorBands_TruthTree(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBands_FluxHistogram(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBandAndFillWithCV_Flux(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBandAndFillWithCV_Genie(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBandAndFillWithCV_BckgConstraint(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBandAndFillWithCV_MichelTrue(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBandAndFillWithCV_MichelFake(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBandAndFillWithCV_TargetMass(MnvHistoType* h);

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

        // Lateral Error Bands
        template<class MnvHistoType>
            void AddLatErrorBands_Data(MnvHistoType* h);

        template<class MnvHistoType>
            void AddLatErrorBands_TruthTree(MnvHistoType* h);

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
            void AddVertErrorBand_Flux(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBand_Genie(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBand_BckgConstraint(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBand_MichelTrue(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBand_MichelFake(MnvHistoType* h);

        template<class MnvHistoType>
            void AddVertErrorBand_TargetMass(MnvHistoType* h);

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
