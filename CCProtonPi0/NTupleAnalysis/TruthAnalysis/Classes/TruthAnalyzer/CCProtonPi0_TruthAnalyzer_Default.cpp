#ifndef CCProtonPi0_TruthAnalyzer_Default_cpp
#define CCProtonPi0_TruthAnalyzer_Default_cpp

#include "CCProtonPi0_TruthAnalyzer.h"

using namespace std;

void CCProtonPi0_TruthAnalyzer::Init(string playlist, TChain* fChain)
{
    // The Init() function is called when the selector needs to initialize
    // a new tree or chain. Typically here the branch addresses and branch
    // pointers of the tree will be set.
    // It is normally not necessary to make changes to the generated
    // code, but the routine can be extended by the user if needed.
    // Init() will be called many times when running on PROOF
    // (once per file to be processed).

    ifstream input_pl(playlist.c_str());
    string filename;

    if( !input_pl.is_open() ){
        cerr<<"Cannot open Playlist File!"<<endl;
        exit(1);
    }else{
        cout<<"Reading Playlist: "<<playlist.c_str()<<endl;
    }

    while (input_pl) {
        input_pl>>filename;

        if (!input_pl) break;

        if (filename[0] != '/') break;

        fChain->Add( filename.c_str() );
        //cout<<" Added "<<filename.c_str()<<endl;
    }

    // Set branch addresses and branch pointers
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("eventID", &eventID, &b_eventID);
    fChain->SetBranchAddress("truth_has_physics_event", &truth_has_physics_event, &b_truth_has_physics_event);
    fChain->SetBranchAddress("truth_isGamma1_conv_inside", &truth_isGamma1_conv_inside, &b_truth_isGamma1_conv_inside);
    fChain->SetBranchAddress("truth_isGamma2_conv_inside", &truth_isGamma2_conv_inside, &b_truth_isGamma2_conv_inside);
    fChain->SetBranchAddress("truth_isSignal", &truth_isSignal, &b_truth_isSignal);
    fChain->SetBranchAddress("truth_isFidVol", &truth_isFidVol, &b_truth_isFidVol);
    fChain->SetBranchAddress("truth_isNC", &truth_isNC, &b_truth_isNC);
    fChain->SetBranchAddress("truth_ReconstructEvent", &truth_ReconstructEvent, &b_truth_ReconstructEvent);
    fChain->SetBranchAddress("truth_isBckg_NoPi0", &truth_isBckg_NoPi0, &b_truth_isBckg_NoPi0);
    fChain->SetBranchAddress("truth_isBckg_SinglePi0", &truth_isBckg_SinglePi0, &b_truth_isBckg_SinglePi0);
    fChain->SetBranchAddress("truth_isBckg_MultiPi0", &truth_isBckg_MultiPi0, &b_truth_isBckg_MultiPi0);
    fChain->SetBranchAddress("truth_isBckg_NC", &truth_isBckg_NC, &b_truth_isBckg_NC);
    fChain->SetBranchAddress("truth_isBckg_AntiNeutrino", &truth_isBckg_AntiNeutrino, &b_truth_isBckg_AntiNeutrino);
    fChain->SetBranchAddress("truth_isBckg_QELike", &truth_isBckg_QELike, &b_truth_isBckg_QELike);
    fChain->SetBranchAddress("truth_isBckg_SingleChargedPion", &truth_isBckg_SingleChargedPion, &b_truth_isBckg_SingleChargedPion);
    fChain->SetBranchAddress("truth_isBckg_SingleChargedPion_ChargeExchanged", &truth_isBckg_SingleChargedPion_ChargeExchanged, &b_truth_isBckg_SingleChargedPion_ChargeExchanged);
    fChain->SetBranchAddress("truth_isBckg_DoublePionWithPi0", &truth_isBckg_DoublePionWithPi0, &b_truth_isBckg_DoublePionWithPi0);
    fChain->SetBranchAddress("truth_isBckg_DoublePionWithoutPi0", &truth_isBckg_DoublePionWithoutPi0, &b_truth_isBckg_DoublePionWithoutPi0);
    fChain->SetBranchAddress("truth_isBckg_MultiPionWithPi0", &truth_isBckg_MultiPionWithPi0, &b_truth_isBckg_MultiPionWithPi0);
    fChain->SetBranchAddress("truth_isBckg_MultiPionWithoutPi0", &truth_isBckg_MultiPionWithoutPi0, &b_truth_isBckg_MultiPionWithoutPi0);
    fChain->SetBranchAddress("truth_isBckg_Other", &truth_isBckg_Other, &b_truth_isBckg_Other);
    fChain->SetBranchAddress("truth_isBckg_withMichel", &truth_isBckg_withMichel, &b_truth_isBckg_withMichel);
    fChain->SetBranchAddress("truth_Bckg_nOther", &truth_Bckg_nOther, &b_truth_Bckg_nOther);
    fChain->SetBranchAddress("truth_Bckg_nPi0_Primary", &truth_Bckg_nPi0_Primary, &b_truth_Bckg_nPi0_Primary);
    fChain->SetBranchAddress("truth_Bckg_nPi0_Secondary", &truth_Bckg_nPi0_Secondary, &b_truth_Bckg_nPi0_Secondary);
    fChain->SetBranchAddress("truth_Bckg_nPi0_Total", &truth_Bckg_nPi0_Total, &b_truth_Bckg_nPi0_Total);
    fChain->SetBranchAddress("truth_Bckg_nPiCharged", &truth_Bckg_nPiCharged, &b_truth_Bckg_nPiCharged);
    fChain->SetBranchAddress("truth_Bckg_nPiCharged_ChargeExchanged", &truth_Bckg_nPiCharged_ChargeExchanged, &b_truth_Bckg_nPiCharged_ChargeExchanged);
    fChain->SetBranchAddress("truth_N_FSParticles", &truth_N_FSParticles, &b_truth_N_FSParticles);
    fChain->SetBranchAddress("truth_N_other", &truth_N_other, &b_truth_N_other);
    fChain->SetBranchAddress("truth_N_pi0", &truth_N_pi0, &b_truth_N_pi0);
    fChain->SetBranchAddress("truth_N_proton", &truth_N_proton, &b_truth_N_proton);
    fChain->SetBranchAddress("truth_N_trueMichelElectrons", &truth_N_trueMichelElectrons, &b_truth_N_trueMichelElectrons);
    fChain->SetBranchAddress("truth_OneShower_evis_most_pdg", &truth_OneShower_evis_most_pdg, &b_truth_OneShower_evis_most_pdg);
    fChain->SetBranchAddress("truth_Rejected_unused_evis_most_pdg", &truth_Rejected_unused_evis_most_pdg, &b_truth_Rejected_unused_evis_most_pdg);
    fChain->SetBranchAddress("truth_ThreeShower_s1_evis_most_pdg", &truth_ThreeShower_s1_evis_most_pdg, &b_truth_ThreeShower_s1_evis_most_pdg);
    fChain->SetBranchAddress("truth_ThreeShower_s2_evis_most_pdg", &truth_ThreeShower_s2_evis_most_pdg, &b_truth_ThreeShower_s2_evis_most_pdg);
    fChain->SetBranchAddress("truth_ThreeShower_s3_evis_most_pdg", &truth_ThreeShower_s3_evis_most_pdg, &b_truth_ThreeShower_s3_evis_most_pdg);
    fChain->SetBranchAddress("truth_blob1_evis_most_pdg", &truth_blob1_evis_most_pdg, &b_truth_blob1_evis_most_pdg);
    fChain->SetBranchAddress("truth_blob2_evis_most_pdg", &truth_blob2_evis_most_pdg, &b_truth_blob2_evis_most_pdg);
    fChain->SetBranchAddress("truth_dispersed_unused_evis_most_pdg", &truth_dispersed_unused_evis_most_pdg, &b_truth_dispersed_unused_evis_most_pdg);
    fChain->SetBranchAddress("truth_pi0_GrandMother", &truth_pi0_GrandMother, &b_truth_pi0_GrandMother);
    fChain->SetBranchAddress("truth_pi0_GrandMotherStatus", &truth_pi0_GrandMotherStatus, &b_truth_pi0_GrandMotherStatus);
    fChain->SetBranchAddress("truth_pi0_Mother", &truth_pi0_Mother, &b_truth_pi0_Mother);
    fChain->SetBranchAddress("truth_pi0_MotherStatus", &truth_pi0_MotherStatus, &b_truth_pi0_MotherStatus);
    fChain->SetBranchAddress("truth_pi0_status", &truth_pi0_status, &b_truth_pi0_status);
    fChain->SetBranchAddress("truth_target_material", &truth_target_material, &b_truth_target_material);
    fChain->SetBranchAddress("truth_vertex_module", &truth_vertex_module, &b_truth_vertex_module);
    fChain->SetBranchAddress("truth_vertex_plane", &truth_vertex_plane, &b_truth_vertex_plane);
    fChain->SetBranchAddress("truth_vertex_unused_evis_most_pdg", &truth_vertex_unused_evis_most_pdg, &b_truth_vertex_unused_evis_most_pdg);
    fChain->SetBranchAddress("truth_OneShower_evis_muon", &truth_OneShower_evis_muon, &b_truth_OneShower_evis_muon);
    fChain->SetBranchAddress("truth_OneShower_evis_neutron", &truth_OneShower_evis_neutron, &b_truth_OneShower_evis_neutron);
    fChain->SetBranchAddress("truth_OneShower_evis_piminus", &truth_OneShower_evis_piminus, &b_truth_OneShower_evis_piminus);
    fChain->SetBranchAddress("truth_OneShower_evis_piplus", &truth_OneShower_evis_piplus, &b_truth_OneShower_evis_piplus);
    fChain->SetBranchAddress("truth_OneShower_evis_pizero", &truth_OneShower_evis_pizero, &b_truth_OneShower_evis_pizero);
    fChain->SetBranchAddress("truth_OneShower_evis_proton", &truth_OneShower_evis_proton, &b_truth_OneShower_evis_proton);
    fChain->SetBranchAddress("truth_OneShower_evis_total_norm", &truth_OneShower_evis_total_norm, &b_truth_OneShower_evis_total_norm);
    fChain->SetBranchAddress("truth_OneShower_evis_total_truth", &truth_OneShower_evis_total_truth, &b_truth_OneShower_evis_total_truth);
    fChain->SetBranchAddress("truth_Rejected_unused_evis_total_norm", &truth_Rejected_unused_evis_total_norm, &b_truth_Rejected_unused_evis_total_norm);
    fChain->SetBranchAddress("truth_Rejected_unused_evis_total_truth", &truth_Rejected_unused_evis_total_truth, &b_truth_Rejected_unused_evis_total_truth);
    fChain->SetBranchAddress("truth_ThreeShower_s1_evis_muon", &truth_ThreeShower_s1_evis_muon, &b_truth_ThreeShower_s1_evis_muon);
    fChain->SetBranchAddress("truth_ThreeShower_s1_evis_neutron", &truth_ThreeShower_s1_evis_neutron, &b_truth_ThreeShower_s1_evis_neutron);
    fChain->SetBranchAddress("truth_ThreeShower_s1_evis_piminus", &truth_ThreeShower_s1_evis_piminus, &b_truth_ThreeShower_s1_evis_piminus);
    fChain->SetBranchAddress("truth_ThreeShower_s1_evis_piplus", &truth_ThreeShower_s1_evis_piplus, &b_truth_ThreeShower_s1_evis_piplus);
    fChain->SetBranchAddress("truth_ThreeShower_s1_evis_pizero", &truth_ThreeShower_s1_evis_pizero, &b_truth_ThreeShower_s1_evis_pizero);
    fChain->SetBranchAddress("truth_ThreeShower_s1_evis_proton", &truth_ThreeShower_s1_evis_proton, &b_truth_ThreeShower_s1_evis_proton);
    fChain->SetBranchAddress("truth_ThreeShower_s1_evis_total_norm", &truth_ThreeShower_s1_evis_total_norm, &b_truth_ThreeShower_s1_evis_total_norm);
    fChain->SetBranchAddress("truth_ThreeShower_s1_evis_total_truth", &truth_ThreeShower_s1_evis_total_truth, &b_truth_ThreeShower_s1_evis_total_truth);
    fChain->SetBranchAddress("truth_ThreeShower_s2_evis_muon", &truth_ThreeShower_s2_evis_muon, &b_truth_ThreeShower_s2_evis_muon);
    fChain->SetBranchAddress("truth_ThreeShower_s2_evis_neutron", &truth_ThreeShower_s2_evis_neutron, &b_truth_ThreeShower_s2_evis_neutron);
    fChain->SetBranchAddress("truth_ThreeShower_s2_evis_piminus", &truth_ThreeShower_s2_evis_piminus, &b_truth_ThreeShower_s2_evis_piminus);
    fChain->SetBranchAddress("truth_ThreeShower_s2_evis_piplus", &truth_ThreeShower_s2_evis_piplus, &b_truth_ThreeShower_s2_evis_piplus);
    fChain->SetBranchAddress("truth_ThreeShower_s2_evis_pizero", &truth_ThreeShower_s2_evis_pizero, &b_truth_ThreeShower_s2_evis_pizero);
    fChain->SetBranchAddress("truth_ThreeShower_s2_evis_proton", &truth_ThreeShower_s2_evis_proton, &b_truth_ThreeShower_s2_evis_proton);
    fChain->SetBranchAddress("truth_ThreeShower_s2_evis_total_norm", &truth_ThreeShower_s2_evis_total_norm, &b_truth_ThreeShower_s2_evis_total_norm);
    fChain->SetBranchAddress("truth_ThreeShower_s2_evis_total_truth", &truth_ThreeShower_s2_evis_total_truth, &b_truth_ThreeShower_s2_evis_total_truth);
    fChain->SetBranchAddress("truth_ThreeShower_s3_evis_muon", &truth_ThreeShower_s3_evis_muon, &b_truth_ThreeShower_s3_evis_muon);
    fChain->SetBranchAddress("truth_ThreeShower_s3_evis_neutron", &truth_ThreeShower_s3_evis_neutron, &b_truth_ThreeShower_s3_evis_neutron);
    fChain->SetBranchAddress("truth_ThreeShower_s3_evis_piminus", &truth_ThreeShower_s3_evis_piminus, &b_truth_ThreeShower_s3_evis_piminus);
    fChain->SetBranchAddress("truth_ThreeShower_s3_evis_piplus", &truth_ThreeShower_s3_evis_piplus, &b_truth_ThreeShower_s3_evis_piplus);
    fChain->SetBranchAddress("truth_ThreeShower_s3_evis_pizero", &truth_ThreeShower_s3_evis_pizero, &b_truth_ThreeShower_s3_evis_pizero);
    fChain->SetBranchAddress("truth_ThreeShower_s3_evis_proton", &truth_ThreeShower_s3_evis_proton, &b_truth_ThreeShower_s3_evis_proton);
    fChain->SetBranchAddress("truth_ThreeShower_s3_evis_total_norm", &truth_ThreeShower_s3_evis_total_norm, &b_truth_ThreeShower_s3_evis_total_norm);
    fChain->SetBranchAddress("truth_ThreeShower_s3_evis_total_truth", &truth_ThreeShower_s3_evis_total_truth, &b_truth_ThreeShower_s3_evis_total_truth);
    fChain->SetBranchAddress("truth_allClusters_evis_pizero", &truth_allClusters_evis_pizero, &b_truth_allClusters_evis_pizero);
    fChain->SetBranchAddress("truth_blob1_evis_muon", &truth_blob1_evis_muon, &b_truth_blob1_evis_muon);
    fChain->SetBranchAddress("truth_blob1_evis_neutron", &truth_blob1_evis_neutron, &b_truth_blob1_evis_neutron);
    fChain->SetBranchAddress("truth_blob1_evis_piminus", &truth_blob1_evis_piminus, &b_truth_blob1_evis_piminus);
    fChain->SetBranchAddress("truth_blob1_evis_piplus", &truth_blob1_evis_piplus, &b_truth_blob1_evis_piplus);
    fChain->SetBranchAddress("truth_blob1_evis_pizero", &truth_blob1_evis_pizero, &b_truth_blob1_evis_pizero);
    fChain->SetBranchAddress("truth_blob1_evis_proton", &truth_blob1_evis_proton, &b_truth_blob1_evis_proton);
    fChain->SetBranchAddress("truth_blob1_evis_total_norm", &truth_blob1_evis_total_norm, &b_truth_blob1_evis_total_norm);
    fChain->SetBranchAddress("truth_blob1_evis_total_truth", &truth_blob1_evis_total_truth, &b_truth_blob1_evis_total_truth);
    fChain->SetBranchAddress("truth_blob2_evis_muon", &truth_blob2_evis_muon, &b_truth_blob2_evis_muon);
    fChain->SetBranchAddress("truth_blob2_evis_neutron", &truth_blob2_evis_neutron, &b_truth_blob2_evis_neutron);
    fChain->SetBranchAddress("truth_blob2_evis_piminus", &truth_blob2_evis_piminus, &b_truth_blob2_evis_piminus);
    fChain->SetBranchAddress("truth_blob2_evis_piplus", &truth_blob2_evis_piplus, &b_truth_blob2_evis_piplus);
    fChain->SetBranchAddress("truth_blob2_evis_pizero", &truth_blob2_evis_pizero, &b_truth_blob2_evis_pizero);
    fChain->SetBranchAddress("truth_blob2_evis_proton", &truth_blob2_evis_proton, &b_truth_blob2_evis_proton);
    fChain->SetBranchAddress("truth_blob2_evis_total_norm", &truth_blob2_evis_total_norm, &b_truth_blob2_evis_total_norm);
    fChain->SetBranchAddress("truth_blob2_evis_total_truth", &truth_blob2_evis_total_truth, &b_truth_blob2_evis_total_truth);
    fChain->SetBranchAddress("truth_dispersed_unused_evis_gamma", &truth_dispersed_unused_evis_gamma, &b_truth_dispersed_unused_evis_gamma);
    fChain->SetBranchAddress("truth_dispersed_unused_evis_muon", &truth_dispersed_unused_evis_muon, &b_truth_dispersed_unused_evis_muon);
    fChain->SetBranchAddress("truth_dispersed_unused_evis_neutron", &truth_dispersed_unused_evis_neutron, &b_truth_dispersed_unused_evis_neutron);
    fChain->SetBranchAddress("truth_dispersed_unused_evis_piminus", &truth_dispersed_unused_evis_piminus, &b_truth_dispersed_unused_evis_piminus);
    fChain->SetBranchAddress("truth_dispersed_unused_evis_piplus", &truth_dispersed_unused_evis_piplus, &b_truth_dispersed_unused_evis_piplus);
    fChain->SetBranchAddress("truth_dispersed_unused_evis_pizero", &truth_dispersed_unused_evis_pizero, &b_truth_dispersed_unused_evis_pizero);
    fChain->SetBranchAddress("truth_dispersed_unused_evis_proton", &truth_dispersed_unused_evis_proton, &b_truth_dispersed_unused_evis_proton);
    fChain->SetBranchAddress("truth_dispersed_unused_evis_total_norm", &truth_dispersed_unused_evis_total_norm, &b_truth_dispersed_unused_evis_total_norm);
    fChain->SetBranchAddress("truth_dispersed_unused_evis_total_truth", &truth_dispersed_unused_evis_total_truth, &b_truth_dispersed_unused_evis_total_truth);
    fChain->SetBranchAddress("truth_ecal_unused_evis_muon", &truth_ecal_unused_evis_muon, &b_truth_ecal_unused_evis_muon);
    fChain->SetBranchAddress("truth_ecal_unused_evis_neutron", &truth_ecal_unused_evis_neutron, &b_truth_ecal_unused_evis_neutron);
    fChain->SetBranchAddress("truth_ecal_unused_evis_piminus", &truth_ecal_unused_evis_piminus, &b_truth_ecal_unused_evis_piminus);
    fChain->SetBranchAddress("truth_ecal_unused_evis_piplus", &truth_ecal_unused_evis_piplus, &b_truth_ecal_unused_evis_piplus);
    fChain->SetBranchAddress("truth_ecal_unused_evis_pizero", &truth_ecal_unused_evis_pizero, &b_truth_ecal_unused_evis_pizero);
    fChain->SetBranchAddress("truth_ecal_unused_evis_proton", &truth_ecal_unused_evis_proton, &b_truth_ecal_unused_evis_proton);
    fChain->SetBranchAddress("truth_ecal_unused_evis_total_norm", &truth_ecal_unused_evis_total_norm, &b_truth_ecal_unused_evis_total_norm);
    fChain->SetBranchAddress("truth_ecal_unused_evis_total_truth", &truth_ecal_unused_evis_total_truth, &b_truth_ecal_unused_evis_total_truth);
    fChain->SetBranchAddress("truth_eventID", &truth_eventID, &b_truth_eventID);
    fChain->SetBranchAddress("truth_hcal_unused_evis_muon", &truth_hcal_unused_evis_muon, &b_truth_hcal_unused_evis_muon);
    fChain->SetBranchAddress("truth_hcal_unused_evis_neutron", &truth_hcal_unused_evis_neutron, &b_truth_hcal_unused_evis_neutron);
    fChain->SetBranchAddress("truth_hcal_unused_evis_piminus", &truth_hcal_unused_evis_piminus, &b_truth_hcal_unused_evis_piminus);
    fChain->SetBranchAddress("truth_hcal_unused_evis_piplus", &truth_hcal_unused_evis_piplus, &b_truth_hcal_unused_evis_piplus);
    fChain->SetBranchAddress("truth_hcal_unused_evis_pizero", &truth_hcal_unused_evis_pizero, &b_truth_hcal_unused_evis_pizero);
    fChain->SetBranchAddress("truth_hcal_unused_evis_proton", &truth_hcal_unused_evis_proton, &b_truth_hcal_unused_evis_proton);
    fChain->SetBranchAddress("truth_hcal_unused_evis_total_norm", &truth_hcal_unused_evis_total_norm, &b_truth_hcal_unused_evis_total_norm);
    fChain->SetBranchAddress("truth_hcal_unused_evis_total_truth", &truth_hcal_unused_evis_total_truth, &b_truth_hcal_unused_evis_total_truth);
    fChain->SetBranchAddress("truth_michelElectron_E", &truth_michelElectron_E, &b_truth_michelElectron_E);
    fChain->SetBranchAddress("truth_michelElectron_P", &truth_michelElectron_P, &b_truth_michelElectron_P);
    fChain->SetBranchAddress("truth_michelMuon_P", &truth_michelMuon_P, &b_truth_michelMuon_P);
    fChain->SetBranchAddress("truth_michelMuon_end_dist_vtx", &truth_michelMuon_end_dist_vtx, &b_truth_michelMuon_end_dist_vtx);
    fChain->SetBranchAddress("truth_michelMuon_length", &truth_michelMuon_length, &b_truth_michelMuon_length);
    fChain->SetBranchAddress("truth_michelPion_P", &truth_michelPion_P, &b_truth_michelPion_P);
    fChain->SetBranchAddress("truth_michelPion_begin_dist_vtx", &truth_michelPion_begin_dist_vtx, &b_truth_michelPion_begin_dist_vtx);
    fChain->SetBranchAddress("truth_michelPion_length", &truth_michelPion_length, &b_truth_michelPion_length);
    fChain->SetBranchAddress("truth_muon_P", &truth_muon_P, &b_truth_muon_P);
    fChain->SetBranchAddress("truth_muon_theta", &truth_muon_theta, &b_truth_muon_theta);
    fChain->SetBranchAddress("truth_other_unused_evis_muon", &truth_other_unused_evis_muon, &b_truth_other_unused_evis_muon);
    fChain->SetBranchAddress("truth_other_unused_evis_neutron", &truth_other_unused_evis_neutron, &b_truth_other_unused_evis_neutron);
    fChain->SetBranchAddress("truth_other_unused_evis_piminus", &truth_other_unused_evis_piminus, &b_truth_other_unused_evis_piminus);
    fChain->SetBranchAddress("truth_other_unused_evis_piplus", &truth_other_unused_evis_piplus, &b_truth_other_unused_evis_piplus);
    fChain->SetBranchAddress("truth_other_unused_evis_pizero", &truth_other_unused_evis_pizero, &b_truth_other_unused_evis_pizero);
    fChain->SetBranchAddress("truth_other_unused_evis_proton", &truth_other_unused_evis_proton, &b_truth_other_unused_evis_proton);
    fChain->SetBranchAddress("truth_other_unused_evis_total_norm", &truth_other_unused_evis_total_norm, &b_truth_other_unused_evis_total_norm);
    fChain->SetBranchAddress("truth_other_unused_evis_total_truth", &truth_other_unused_evis_total_truth, &b_truth_other_unused_evis_total_truth);
    fChain->SetBranchAddress("truth_pi0_P", &truth_pi0_P, &b_truth_pi0_P);
    fChain->SetBranchAddress("truth_pi0_theta", &truth_pi0_theta, &b_truth_pi0_theta);
    fChain->SetBranchAddress("truth_proton_P", &truth_proton_P, &b_truth_proton_P);
    fChain->SetBranchAddress("truth_proton_theta", &truth_proton_theta, &b_truth_proton_theta);
    fChain->SetBranchAddress("truth_total_captured_evis_pizero", &truth_total_captured_evis_pizero, &b_truth_total_captured_evis_pizero);
    fChain->SetBranchAddress("truth_total_captured_evis_total_norm", &truth_total_captured_evis_total_norm, &b_truth_total_captured_evis_total_norm);
    fChain->SetBranchAddress("truth_total_captured_evis_total_truth", &truth_total_captured_evis_total_truth, &b_truth_total_captured_evis_total_truth);
    fChain->SetBranchAddress("truth_vertex_unused_evis_gamma", &truth_vertex_unused_evis_gamma, &b_truth_vertex_unused_evis_gamma);
    fChain->SetBranchAddress("truth_vertex_unused_evis_muon", &truth_vertex_unused_evis_muon, &b_truth_vertex_unused_evis_muon);
    fChain->SetBranchAddress("truth_vertex_unused_evis_neutron", &truth_vertex_unused_evis_neutron, &b_truth_vertex_unused_evis_neutron);
    fChain->SetBranchAddress("truth_vertex_unused_evis_piminus", &truth_vertex_unused_evis_piminus, &b_truth_vertex_unused_evis_piminus);
    fChain->SetBranchAddress("truth_vertex_unused_evis_piplus", &truth_vertex_unused_evis_piplus, &b_truth_vertex_unused_evis_piplus);
    fChain->SetBranchAddress("truth_vertex_unused_evis_pizero", &truth_vertex_unused_evis_pizero, &b_truth_vertex_unused_evis_pizero);
    fChain->SetBranchAddress("truth_vertex_unused_evis_proton", &truth_vertex_unused_evis_proton, &b_truth_vertex_unused_evis_proton);
    fChain->SetBranchAddress("truth_vertex_unused_evis_total_norm", &truth_vertex_unused_evis_total_norm, &b_truth_vertex_unused_evis_total_norm);
    fChain->SetBranchAddress("truth_vertex_unused_evis_total_truth", &truth_vertex_unused_evis_total_truth, &b_truth_vertex_unused_evis_total_truth);
    fChain->SetBranchAddress("truth_gamma1_4P", truth_gamma1_4P, &b_truth_gamma1_4P);
    fChain->SetBranchAddress("truth_gamma1_final_pos", truth_gamma1_final_pos, &b_truth_gamma1_final_pos);
    fChain->SetBranchAddress("truth_gamma1_init_pos", truth_gamma1_init_pos, &b_truth_gamma1_init_pos);
    fChain->SetBranchAddress("truth_gamma2_4P", truth_gamma2_4P, &b_truth_gamma2_4P);
    fChain->SetBranchAddress("truth_gamma2_final_pos", truth_gamma2_final_pos, &b_truth_gamma2_final_pos);
    fChain->SetBranchAddress("truth_gamma2_init_pos", truth_gamma2_init_pos, &b_truth_gamma2_init_pos);
    fChain->SetBranchAddress("genie_wgt_n_shifts", &genie_wgt_n_shifts, &b_genie_wgt_n_shifts);
    fChain->SetBranchAddress("truth_genie_wgt_AGKYxF1pi", truth_genie_wgt_AGKYxF1pi, &b_truth_genie_wgt_AGKYxF1pi);
    fChain->SetBranchAddress("truth_genie_wgt_AhtBY", truth_genie_wgt_AhtBY, &b_truth_genie_wgt_AhtBY);
    fChain->SetBranchAddress("truth_genie_wgt_BhtBY", truth_genie_wgt_BhtBY, &b_truth_genie_wgt_BhtBY);
    fChain->SetBranchAddress("truth_genie_wgt_CCQEPauliSupViaKF", truth_genie_wgt_CCQEPauliSupViaKF, &b_truth_genie_wgt_CCQEPauliSupViaKF);
    fChain->SetBranchAddress("truth_genie_wgt_CV1uBY", truth_genie_wgt_CV1uBY, &b_truth_genie_wgt_CV1uBY);
    fChain->SetBranchAddress("truth_genie_wgt_CV2uBY", truth_genie_wgt_CV2uBY, &b_truth_genie_wgt_CV2uBY);
    fChain->SetBranchAddress("truth_genie_wgt_EtaNCEL", truth_genie_wgt_EtaNCEL, &b_truth_genie_wgt_EtaNCEL);
    fChain->SetBranchAddress("truth_genie_wgt_FrAbs_N", truth_genie_wgt_FrAbs_N, &b_truth_genie_wgt_FrAbs_N);
    fChain->SetBranchAddress("truth_genie_wgt_FrAbs_pi", truth_genie_wgt_FrAbs_pi, &b_truth_genie_wgt_FrAbs_pi);
    fChain->SetBranchAddress("truth_genie_wgt_FrCEx_N", truth_genie_wgt_FrCEx_N, &b_truth_genie_wgt_FrCEx_N);
    fChain->SetBranchAddress("truth_genie_wgt_FrCEx_pi", truth_genie_wgt_FrCEx_pi, &b_truth_genie_wgt_FrCEx_pi);
    fChain->SetBranchAddress("truth_genie_wgt_FrElas_N", truth_genie_wgt_FrElas_N, &b_truth_genie_wgt_FrElas_N);
    fChain->SetBranchAddress("truth_genie_wgt_FrElas_pi", truth_genie_wgt_FrElas_pi, &b_truth_genie_wgt_FrElas_pi);
    fChain->SetBranchAddress("truth_genie_wgt_FrInel_N", truth_genie_wgt_FrInel_N, &b_truth_genie_wgt_FrInel_N);
    fChain->SetBranchAddress("truth_genie_wgt_FrInel_pi", truth_genie_wgt_FrInel_pi, &b_truth_genie_wgt_FrInel_pi);
    fChain->SetBranchAddress("truth_genie_wgt_FrPiProd_N", truth_genie_wgt_FrPiProd_N, &b_truth_genie_wgt_FrPiProd_N);
    fChain->SetBranchAddress("truth_genie_wgt_FrPiProd_pi", truth_genie_wgt_FrPiProd_pi, &b_truth_genie_wgt_FrPiProd_pi);
    fChain->SetBranchAddress("truth_genie_wgt_MFP_N", truth_genie_wgt_MFP_N, &b_truth_genie_wgt_MFP_N);
    fChain->SetBranchAddress("truth_genie_wgt_MFP_pi", truth_genie_wgt_MFP_pi, &b_truth_genie_wgt_MFP_pi);
    fChain->SetBranchAddress("truth_genie_wgt_MaCCQE", truth_genie_wgt_MaCCQE, &b_truth_genie_wgt_MaCCQE);
    fChain->SetBranchAddress("truth_genie_wgt_MaCCQEshape", truth_genie_wgt_MaCCQEshape, &b_truth_genie_wgt_MaCCQEshape);
    fChain->SetBranchAddress("truth_genie_wgt_MaNCEL", truth_genie_wgt_MaNCEL, &b_truth_genie_wgt_MaNCEL);
    fChain->SetBranchAddress("truth_genie_wgt_MaRES", truth_genie_wgt_MaRES, &b_truth_genie_wgt_MaRES);
    fChain->SetBranchAddress("truth_genie_wgt_MvRES", truth_genie_wgt_MvRES, &b_truth_genie_wgt_MvRES);
    fChain->SetBranchAddress("truth_genie_wgt_NormCCQE", truth_genie_wgt_NormCCQE, &b_truth_genie_wgt_NormCCQE);
    fChain->SetBranchAddress("truth_genie_wgt_NormCCRES", truth_genie_wgt_NormCCRES, &b_truth_genie_wgt_NormCCRES);
    fChain->SetBranchAddress("truth_genie_wgt_NormDISCC", truth_genie_wgt_NormDISCC, &b_truth_genie_wgt_NormDISCC);
    fChain->SetBranchAddress("truth_genie_wgt_NormNCRES", truth_genie_wgt_NormNCRES, &b_truth_genie_wgt_NormNCRES);
    fChain->SetBranchAddress("truth_genie_wgt_RDecBR1gamma", truth_genie_wgt_RDecBR1gamma, &b_truth_genie_wgt_RDecBR1gamma);
    fChain->SetBranchAddress("truth_genie_wgt_Rvn1pi", truth_genie_wgt_Rvn1pi, &b_truth_genie_wgt_Rvn1pi);
    fChain->SetBranchAddress("truth_genie_wgt_Rvn2pi", truth_genie_wgt_Rvn2pi, &b_truth_genie_wgt_Rvn2pi);
    fChain->SetBranchAddress("truth_genie_wgt_Rvp1pi", truth_genie_wgt_Rvp1pi, &b_truth_genie_wgt_Rvp1pi);
    fChain->SetBranchAddress("truth_genie_wgt_Rvp2pi", truth_genie_wgt_Rvp2pi, &b_truth_genie_wgt_Rvp2pi);
    fChain->SetBranchAddress("truth_genie_wgt_Theta_Delta2Npi", truth_genie_wgt_Theta_Delta2Npi, &b_truth_genie_wgt_Theta_Delta2Npi);
    fChain->SetBranchAddress("truth_genie_wgt_VecFFCCQEshape", truth_genie_wgt_VecFFCCQEshape, &b_truth_genie_wgt_VecFFCCQEshape);
    fChain->SetBranchAddress("truth_genie_wgt_shifts", truth_genie_wgt_shifts, &b_truth_genie_wgt_shifts);
    fChain->SetBranchAddress("truth_michelMuon_endPoint", truth_michelMuon_endPoint, &b_truth_michelMuon_endPoint);
    fChain->SetBranchAddress("truth_muon_4P", truth_muon_4P, &b_truth_muon_4P);
    fChain->SetBranchAddress("truth_pi0_4P", truth_pi0_4P, &b_truth_pi0_4P);
    fChain->SetBranchAddress("truth_proton_4P", truth_proton_4P, &b_truth_proton_4P);
    fChain->SetBranchAddress("mc_run", &mc_run, &b_mc_run);
    fChain->SetBranchAddress("mc_subrun", &mc_subrun, &b_mc_subrun);
    fChain->SetBranchAddress("mc_nInteractions", &mc_nInteractions, &b_mc_nInteractions);
    fChain->SetBranchAddress("mc_MIState", &mc_MIState, &b_mc_MIState);
    fChain->SetBranchAddress("mc_pot", &mc_pot, &b_mc_pot);
    fChain->SetBranchAddress("mc_beamConfig", &mc_beamConfig, &b_mc_beamConfig);
    fChain->SetBranchAddress("mc_processType", &mc_processType, &b_mc_processType);
    fChain->SetBranchAddress("mc_nthEvtInSpill", &mc_nthEvtInSpill, &b_mc_nthEvtInSpill);
    fChain->SetBranchAddress("mc_nthEvtInFile", &mc_nthEvtInFile, &b_mc_nthEvtInFile);
    fChain->SetBranchAddress("mc_intType", &mc_intType, &b_mc_intType);
    fChain->SetBranchAddress("mc_current", &mc_current, &b_mc_current);
    fChain->SetBranchAddress("mc_charm", &mc_charm, &b_mc_charm);
    fChain->SetBranchAddress("mc_weight", &mc_weight, &b_mc_weight);
    fChain->SetBranchAddress("mc_XSec", &mc_XSec, &b_mc_XSec);
    fChain->SetBranchAddress("mc_diffXSec", &mc_diffXSec, &b_mc_diffXSec);
    fChain->SetBranchAddress("mc_incoming", &mc_incoming, &b_mc_incoming);
    fChain->SetBranchAddress("mc_fluxDriverProb", &mc_fluxDriverProb, &b_mc_fluxDriverProb);
    fChain->SetBranchAddress("mc_targetNucleus", &mc_targetNucleus, &b_mc_targetNucleus);
    fChain->SetBranchAddress("mc_targetZ", &mc_targetZ, &b_mc_targetZ);
    fChain->SetBranchAddress("mc_targetA", &mc_targetA, &b_mc_targetA);
    fChain->SetBranchAddress("mc_targetNucleon", &mc_targetNucleon, &b_mc_targetNucleon);
    fChain->SetBranchAddress("mc_struckQuark", &mc_struckQuark, &b_mc_struckQuark);
    fChain->SetBranchAddress("mc_seaQuark", &mc_seaQuark, &b_mc_seaQuark);
    fChain->SetBranchAddress("mc_resID", &mc_resID, &b_mc_resID);
    fChain->SetBranchAddress("mc_primaryLepton", &mc_primaryLepton, &b_mc_primaryLepton);
    fChain->SetBranchAddress("mc_incomingE", &mc_incomingE, &b_mc_incomingE);
    fChain->SetBranchAddress("mc_Bjorkenx", &mc_Bjorkenx, &b_mc_Bjorkenx);
    fChain->SetBranchAddress("mc_Bjorkeny", &mc_Bjorkeny, &b_mc_Bjorkeny);
    fChain->SetBranchAddress("mc_Q2", &mc_Q2, &b_mc_Q2);
    fChain->SetBranchAddress("mc_nuT", &mc_nuT, &b_mc_nuT);
    fChain->SetBranchAddress("mc_w", &mc_w, &b_mc_w);
    fChain->SetBranchAddress("mc_vtx", mc_vtx, &b_mc_vtx);
    fChain->SetBranchAddress("mc_incomingPartVec", mc_incomingPartVec, &b_mc_incomingPartVec);
    fChain->SetBranchAddress("mc_initNucVec", mc_initNucVec, &b_mc_initNucVec);
    fChain->SetBranchAddress("mc_primFSLepton", mc_primFSLepton, &b_mc_primFSLepton);
    fChain->SetBranchAddress("mc_nFSPart", &mc_nFSPart, &b_mc_nFSPart);
    fChain->SetBranchAddress("mc_FSPartPx", mc_FSPartPx, &b_mc_FSPartPx);
    fChain->SetBranchAddress("mc_FSPartPy", mc_FSPartPy, &b_mc_FSPartPy);
    fChain->SetBranchAddress("mc_FSPartPz", mc_FSPartPz, &b_mc_FSPartPz);
    fChain->SetBranchAddress("mc_FSPartE", mc_FSPartE, &b_mc_FSPartE);
    fChain->SetBranchAddress("mc_FSPartPDG", mc_FSPartPDG, &b_mc_FSPartPDG);
    fChain->SetBranchAddress("mc_er_nPart", &mc_er_nPart, &b_mc_er_nPart);
    fChain->SetBranchAddress("mc_er_ID", mc_er_ID, &b_mc_er_ID);
    fChain->SetBranchAddress("mc_er_status", mc_er_status, &b_mc_er_status);
    fChain->SetBranchAddress("mc_er_posInNucX", mc_er_posInNucX, &b_mc_er_posInNucX);
    fChain->SetBranchAddress("mc_er_posInNucY", mc_er_posInNucY, &b_mc_er_posInNucY);
    fChain->SetBranchAddress("mc_er_posInNucZ", mc_er_posInNucZ, &b_mc_er_posInNucZ);
    fChain->SetBranchAddress("mc_er_Px", mc_er_Px, &b_mc_er_Px);
    fChain->SetBranchAddress("mc_er_Py", mc_er_Py, &b_mc_er_Py);
    fChain->SetBranchAddress("mc_er_Pz", mc_er_Pz, &b_mc_er_Pz);
    fChain->SetBranchAddress("mc_er_E", mc_er_E, &b_mc_er_E);
    fChain->SetBranchAddress("mc_er_FD", mc_er_FD, &b_mc_er_FD);
    fChain->SetBranchAddress("mc_er_LD", mc_er_LD, &b_mc_er_LD);
    fChain->SetBranchAddress("mc_er_mother", mc_er_mother, &b_mc_er_mother);
    fChain->SetBranchAddress("mc_fr_nNuAncestorIDs", &mc_fr_nNuAncestorIDs, &b_mc_fr_nNuAncestorIDs);
    fChain->SetBranchAddress("mc_fr_nuAncestorIDs", mc_fr_nuAncestorIDs, &b_mc_fr_nuAncestorIDs);
    fChain->SetBranchAddress("mc_fr_nuParentID", &mc_fr_nuParentID, &b_mc_fr_nuParentID);
    fChain->SetBranchAddress("mc_fr_decMode", &mc_fr_decMode, &b_mc_fr_decMode);
    fChain->SetBranchAddress("mc_fr_primProtonVtx", mc_fr_primProtonVtx, &b_mc_fr_primProtonVtx);
    fChain->SetBranchAddress("mc_fr_primProtonP", mc_fr_primProtonP, &b_mc_fr_primProtonP);
    fChain->SetBranchAddress("mc_fr_nuParentDecVtx", mc_fr_nuParentDecVtx, &b_mc_fr_nuParentDecVtx);
    fChain->SetBranchAddress("mc_fr_nuParentProdVtx", mc_fr_nuParentProdVtx, &b_mc_fr_nuParentProdVtx);
    fChain->SetBranchAddress("mc_fr_nuParentProdP", mc_fr_nuParentProdP, &b_mc_fr_nuParentProdP);
    fChain->SetBranchAddress("mc_cvweight_total", &mc_cvweight_total, &b_mc_cvweight_total);
    fChain->SetBranchAddress("wgt", &wgt, &b_wgt);
    fChain->SetBranchAddress("mc_cvweight_totalFlux", &mc_cvweight_totalFlux, &b_mc_cvweight_totalFlux);
    fChain->SetBranchAddress("mc_cvweight_totalXsec", &mc_cvweight_totalXsec, &b_mc_cvweight_totalXsec);
    fChain->SetBranchAddress("mc_ppfx1_cvweight", &mc_ppfx1_cvweight, &b_mc_ppfx1_cvweight);
    fChain->SetBranchAddress("mc_gen1_cvweight_total", &mc_gen1_cvweight_total, &b_mc_gen1_cvweight_total);
    fChain->SetBranchAddress("gen1_wgt", &gen1_wgt, &b_gen1_wgt);
    fChain->SetBranchAddress("mc_gen1_cvweight_totalFlux", &mc_gen1_cvweight_totalFlux, &b_mc_gen1_cvweight_totalFlux);
    fChain->SetBranchAddress("mc_gen1_cvweight_NA49", &mc_gen1_cvweight_NA49, &b_mc_gen1_cvweight_NA49);
    fChain->SetBranchAddress("mc_wgt_Flux_BeamFocus_sz", &mc_wgt_Flux_BeamFocus_sz, &b_mc_wgt_Flux_BeamFocus_sz);
    fChain->SetBranchAddress("mc_wgt_Flux_BeamFocus", mc_wgt_Flux_BeamFocus, &b_mc_wgt_Flux_BeamFocus);
    fChain->SetBranchAddress("mc_wgt_gen1_Flux_Tertiary_sz", &mc_wgt_gen1_Flux_Tertiary_sz, &b_mc_wgt_gen1_Flux_Tertiary_sz);
    fChain->SetBranchAddress("mc_wgt_gen1_Flux_Tertiary", mc_wgt_gen1_Flux_Tertiary, &b_mc_wgt_gen1_Flux_Tertiary);
    fChain->SetBranchAddress("mc_wgt_gen1_Flux_NA49_sz", &mc_wgt_gen1_Flux_NA49_sz, &b_mc_wgt_gen1_Flux_NA49_sz);
    fChain->SetBranchAddress("mc_wgt_gen1_Flux_NA49", mc_wgt_gen1_Flux_NA49, &b_mc_wgt_gen1_Flux_NA49);
    fChain->SetBranchAddress("mc_wgt_Norm_sz", &mc_wgt_Norm_sz, &b_mc_wgt_Norm_sz);
    fChain->SetBranchAddress("mc_wgt_Norm", &mc_wgt_Norm, &b_mc_wgt_Norm);
    fChain->SetBranchAddress("mc_wgt_ppfx1_Total_sz", &mc_wgt_ppfx1_Total_sz, &b_mc_wgt_ppfx1_Total_sz);
    fChain->SetBranchAddress("mc_wgt_ppfx1_Total", mc_wgt_ppfx1_Total, &b_mc_wgt_ppfx1_Total);



}

CCProtonPi0_TruthAnalyzer::~CCProtonPi0_TruthAnalyzer()
{
    // Do Nothing!
}

Int_t CCProtonPi0_TruthAnalyzer::GetEntry(Long64_t entry)
{
    // Read contents of entry.
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}

Long64_t CCProtonPi0_TruthAnalyzer::LoadTree(Long64_t entry)
{
    // Set the environment to read one entry
    if (!fChain) return -5;
    Long64_t centry = fChain->LoadTree(entry);
    if (centry < 0) return centry;
    if (fChain->GetTreeNumber() != fCurrent) {
        fCurrent = fChain->GetTreeNumber();
    }
    return centry;
}

#endif

