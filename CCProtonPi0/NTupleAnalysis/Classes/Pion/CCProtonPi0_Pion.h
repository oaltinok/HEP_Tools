/*
================================================================================
Class: CCProtonPi0_Pion -> Derived Class from Particle Base Clas
    CCProtonPi0_Pion Class  inherits Particle Behaviours and 
                extends base class with pion specific parameters

    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
================================================================================
*/

#ifndef CCProtonPi0_Pion_h
#define CCProtonPi0_Pion_h

#include "../Particle/CCProtonPi0_Particle.h"

using namespace PlotUtils;

class CCProtonPi0_Pion : public CCProtonPi0_Particle
{
    public:
        vector<MnvH1D*> invMass;
        vector<MnvH1D*> invMass_Old;
        vector<MnvH1D*> photonEnergy_Asymmetry;

        // Truth Match
        vector<MnvH1D*> evis_frac_true_pi0_reco_all;
        vector<MnvH1D*> evis_frac_reco_pi0_true_pi0;
        vector<MnvH1D*> evis_frac_reco_pi0_reco_all;
        vector<MnvH1D*> evis_frac_reco_nonpi0_reco_all;

        TH1D* g1_evis_proton;
        TH1D* g1_evis_neutron;
        TH1D* g1_evis_pi;
        TH1D* g1_evis_pi0;
        TH1D* g1_evis_muon;

        TH1D* g2_evis_proton;
        TH1D* g2_evis_neutron;
        TH1D* g2_evis_pi;
        TH1D* g2_evis_pi0;
        TH1D* g2_evis_muon;

        TH1D* g3_evis_proton;
        TH1D* g3_evis_neutron;
        TH1D* g3_evis_pi;
        TH1D* g3_evis_pi0;
        TH1D* g3_evis_muon;

        
        // Leading Photon - Energetic Photon
        vector<MnvH1D*> gamma1_E;
        vector<MnvH1D*> gamma1_E_Old;
        vector<MnvH1D*> gamma1_theta;
        vector<MnvH1D*> gamma1_ConvLength;
        TH1D* gamma1_true_E;
        TH1D* gamma1_reco_error_E;
        TH1D* gamma1_reco_Old_error_E;
        TH2D* gamma1_reco_E_true_E;
        TH2D* gamma1_true_E_reco_E_error;

        // Secondary Photon
        vector<MnvH1D*> gamma2_E;
        vector<MnvH1D*> gamma2_E_Old;
        vector<MnvH1D*> gamma2_theta;
        vector<MnvH1D*> gamma2_ConvLength;
        TH1D* gamma2_true_E;
        TH1D* gamma2_reco_error_E;
        TH1D* gamma2_reco_Old_error_E;
        TH2D* gamma2_reco_E_true_E;
        TH2D* gamma2_true_E_reco_E_error;
                
        TH2D* gamma1_E_gamma2_E;
        TH2D* gamma1_convLength_gamma2_convLength;
                
        // Energy Study Histograms -- Temporary
        TH2D* true_E_evis_trkr_ratio;
        TH2D* evis_evis_trkr_ratio;
        TH2D* evis_evis_ecal_ratio;
        TH1D* reco_error_trkr;
        TH1D* calc_error_trkr;
        TH1D* reco_error_ecal;
        TH1D* calc_error_ecal;
        TH1D* reco_error_trkr_ecal;
        TH1D* calc_error_trkr_ecal;

        TH2D* reco_E_true_E;

        TH1D* evis_trkr_ratio_1;
        TH1D* evis_trkr_ratio_2;
        TH1D* evis_trkr_ratio_3;
        TH1D* evis_trkr_ratio_4;
        TH1D* evis_trkr_ratio_5;
        TH1D* evis_trkr_ratio_6;
 
        TH1D* evis_ecal_ratio_1;
        TH1D* evis_ecal_ratio_2;
        TH1D* evis_ecal_ratio_3;
        TH1D* evis_ecal_ratio_4;
        TH1D* evis_ecal_ratio_5;
        
        TH1D* gamma_nPlanes;
        TH1D* gamma1_calc_error;
        TH1D* gamma2_calc_error;
        TH1D* mgg_calc;

        TH1D* Shower_Topology;
        TH1D* energy_frac_trkr;
        TH1D* energy_frac_ecal;
        TH1D* energy_frac_scal;
        TH1D* energy_frac_hcal;
       
        TH1D* energy_trkr;
        TH1D* energy_ecal;
        TH1D* energy_scal;
        TH1D* energy_hcal;
        
        TH2D* evis_trkr_reco_true;
        TH2D* evis_ecal_reco_true;
        TH2D* evis_scal_reco_true;
        TH2D* evis_hcal_reco_true;
        
        TH2D* energy_trkr_reco_true;
        TH2D* energy_ecal_reco_true;
        TH2D* energy_scal_reco_true;
        TH2D* energy_hcal_reco_true;

        // Side ECAL nHits Study
        TH2D* trkr_nHits_reco_correct;
        TH2D* trkr_nHits_reco_true;
        TH2D* scal_nHits_reco_correct;
        TH2D* scal_nHits_reco_true;
        
        CCProtonPi0_Pion(bool isModeReduce, bool isMC, std::string ana_folder);
        void initHistograms();
        void writeHistograms();

    private:
        CCProtonPi0_SingleBin bin_invMass;
        CCProtonPi0_SingleBin bin_photonEnergy_Asymmetry;
        CCProtonPi0_SingleBin bin_photonConvLength;
        CCProtonPi0_SingleBin bin_photonP;
};

#endif

