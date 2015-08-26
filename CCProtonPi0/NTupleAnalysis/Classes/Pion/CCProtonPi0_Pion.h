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
        vector<MnvH1D*> gamma1_theta;
        vector<MnvH1D*> gamma1_ConvLength;
        vector<MnvH1D*> g1_evis_trkr;
        vector<MnvH1D*> g1_evis_scal;
        vector<MnvH1D*> g1_evis_ecal;
        vector<MnvH1D*> g1_evis_hcal;

        // Secondary Photon
        vector<MnvH1D*> gamma2_E;
        vector<MnvH1D*> gamma2_theta;
        vector<MnvH1D*> gamma2_ConvLength;
        vector<MnvH1D*> g2_evis_trkr;
        vector<MnvH1D*> g2_evis_scal;
        vector<MnvH1D*> g2_evis_ecal;
        vector<MnvH1D*> g2_evis_hcal;
        vector<MnvH1D*> g2_evis_frac_scal_trkr;

        // Photon Comparison
        TH1D* gamma1_true_E;
        TH1D* gamma2_true_E;
        TH2D* gamma1_E_gamma2_E;
        TH2D* gamma1_convLength_gamma2_convLength;
        TH2D* gamma1_reco_E_true_E;
        TH2D* gamma2_reco_E_true_E;
        TH1D* gamma2_E_error;
        TH1D* gamma1_E_error;
        TH1D* mgg_reco;
        TH1D* mgg_true;
        TH1D* mgg_error;
        TH2D* mgg_reco_true;
        TH1D* isGamma1_conv_inside;
        TH1D* isGamma2_conv_inside;
        
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

