#ifndef Other_Plots_cpp
#define Other_Plots_cpp

#include "Plotter.h"

using namespace std;

void Plotter::plotCutHistograms()
{
    TFile* f_Root_Signal_Interaction = new TFile(rootDir_Interaction[0].c_str());
    TFile* f_Root_Background_Interaction = new TFile(rootDir_Interaction[1].c_str());
    
    TH1D* h_signal_hCut_vertexCount= (TH1D*)f_Root_Signal_Interaction->Get("hCut_vertexCount");
    TH1D* h_background_hCut_vertexCount = (TH1D*)f_Root_Background_Interaction->Get("hCut_vertexCount");
    plotStacked(h_signal_hCut_vertexCount, h_background_hCut_vertexCount,"Number of Vertices", "hCut_vertexCount.png", otherDir);
    
    TH1D* h_signal_hCut_nProngs= (TH1D*)f_Root_Signal_Interaction->Get("hCut_nProngs");
    TH1D* h_background_hCut_nProngs = (TH1D*)f_Root_Background_Interaction->Get("hCut_nProngs");
    plotStacked(h_signal_hCut_nProngs, h_background_hCut_nProngs,"Number of Prongs at Primary Vertex", "hCut_nProngs.png", otherDir);
    
    TH1D* h_signal_hCut_1Prong_Michel = (TH1D*)f_Root_Signal_Interaction->Get("hCut_1Prong_Michel");
    TH1D* h_background_hCut_1Prong_Michel = (TH1D*)f_Root_Background_Interaction->Get("hCut_1Prong_Michel");
    plotStacked(h_signal_hCut_1Prong_Michel, h_background_hCut_1Prong_Michel,"Does Event has Michel? (0 = No, 1 = Yes) (1Prong)", "hCut_1Prong_Michel.png", otherDir);
    
    TH1D* h_signal_hCut_1Prong_eVis_nuclearTarget = (TH1D*)f_Root_Signal_Interaction->Get("hCut_1Prong_eVis_nuclearTarget");
    TH1D* h_background_hCut_1Prong_eVis_nuclearTarget = (TH1D*)f_Root_Background_Interaction->Get("hCut_1Prong_eVis_nuclearTarget");
    plotStacked(h_signal_hCut_1Prong_eVis_nuclearTarget, h_background_hCut_1Prong_eVis_nuclearTarget,"Visible Energy inside Nuclear Target (1Prong)", "hCut_1Prong_eVis_nuclearTarget.png", otherDir);
    
    TH1D* h_signal_hCut_1Prong_eVis_other = (TH1D*)f_Root_Signal_Interaction->Get("hCut_1Prong_eVis_other");
    TH1D* h_background_hCut_1Prong_eVis_other = (TH1D*)f_Root_Background_Interaction->Get("hCut_1Prong_eVis_other");
    plotStacked(h_signal_hCut_1Prong_eVis_other, h_background_hCut_1Prong_eVis_other,"Visible Energy at Target + ECAL + HCAL (1Prong)", "hCut_1Prong_eVis_other.png", otherDir);
    
    TH1D* h_signal_hCut_1Prong_gamma1ConvDist = (TH1D*)f_Root_Signal_Interaction->Get("hCut_1Prong_gamma1ConvDist");
    TH1D* h_background_hCut_1Prong_gamma1ConvDist = (TH1D*)f_Root_Background_Interaction->Get("hCut_1Prong_gamma1ConvDist");
    plotStacked(h_signal_hCut_1Prong_gamma1ConvDist, h_background_hCut_1Prong_gamma1ConvDist,"Leading Photon Conversion Distance [cm] (1Prong)", "hCut_1Prong_gamma1ConvDist.png", otherDir);
    
    TH1D* h_signal_hCut_1Prong_gamma2ConvDist = (TH1D*)f_Root_Signal_Interaction->Get("hCut_1Prong_gamma2ConvDist");
    TH1D* h_background_hCut_1Prong_gamma2ConvDist = (TH1D*)f_Root_Background_Interaction->Get("hCut_1Prong_gamma2ConvDist");
    plotStacked(h_signal_hCut_1Prong_gamma2ConvDist, h_background_hCut_1Prong_gamma2ConvDist,"Second Photon Conversion Distance [cm] (1Prong)", "hCut_1Prong_gamma2ConvDist.png", otherDir);
    
    TH1D* h_signal_hCut_1Prong_pi0invMass = (TH1D*)f_Root_Signal_Interaction->Get("hCut_1Prong_pi0invMass");
    TH1D* h_background_hCut_1Prong_pi0invMass = (TH1D*)f_Root_Background_Interaction->Get("hCut_1Prong_pi0invMass");
    plotStacked(h_signal_hCut_1Prong_pi0invMass, h_background_hCut_1Prong_pi0invMass,"Pi0 Invariant Mass (1Prong)", "hCut_1Prong_pi0invMass.png", otherDir);
    
    TH1D* h_signal_hCut_1Prong_neutrinoE = (TH1D*)f_Root_Signal_Interaction->Get("hCut_1Prong_neutrinoE");
    TH1D* h_background_hCut_1Prong_neutrinoE = (TH1D*)f_Root_Background_Interaction->Get("hCut_1Prong_neutrinoE");
    plotStacked(h_signal_hCut_1Prong_neutrinoE, h_background_hCut_1Prong_neutrinoE,"Neutrino Energy [GeV] (1Prong)", "hCut_1Prong_neutrinoE.png", otherDir);
    
    TH1D* h_signal_hCut_1Prong_UnusedE = (TH1D*)f_Root_Signal_Interaction->Get("hCut_1Prong_UnusedE");
    TH1D* h_background_hCut_1Prong_UnusedE = (TH1D*)f_Root_Background_Interaction->Get("hCut_1Prong_UnusedE");
    plotStacked(h_signal_hCut_1Prong_UnusedE, h_background_hCut_1Prong_UnusedE,"Unused Cluster Energy at the end of Reconstruction (1Prong)", "hCut_1Prong_UnusedE.png", otherDir);
    
    TH1D* h_signal_hCut_2Prong_Michel = (TH1D*)f_Root_Signal_Interaction->Get("hCut_2Prong_Michel");
    TH1D* h_background_hCut_2Prong_Michel = (TH1D*)f_Root_Background_Interaction->Get("hCut_2Prong_Michel");
    plotStacked(h_signal_hCut_2Prong_Michel, h_background_hCut_2Prong_Michel,"Does Event has Michel? (0 = No, 1 = Yes) (2Prong)", "hCut_2Prong_Michel.png", otherDir);
    
    TH1D* h_signal_hCut_2Prong_eVis_nuclearTarget = (TH1D*)f_Root_Signal_Interaction->Get("hCut_2Prong_eVis_nuclearTarget");
    TH1D* h_background_hCut_2Prong_eVis_nuclearTarget = (TH1D*)f_Root_Background_Interaction->Get("hCut_2Prong_eVis_nuclearTarget");
    plotStacked(h_signal_hCut_2Prong_eVis_nuclearTarget, h_background_hCut_2Prong_eVis_nuclearTarget,"Visible Energy inside Nuclear Target (2Prong)", "hCut_2Prong_eVis_nuclearTarget.png", otherDir);
    
    TH1D* h_signal_hCut_2Prong_eVis_other = (TH1D*)f_Root_Signal_Interaction->Get("hCut_2Prong_eVis_other");
    TH1D* h_background_hCut_2Prong_eVis_other = (TH1D*)f_Root_Background_Interaction->Get("hCut_2Prong_eVis_other");
    plotStacked(h_signal_hCut_2Prong_eVis_other, h_background_hCut_2Prong_eVis_other,"Visible Energy at Target + ECAL + HCAL (2Prong)", "hCut_2Prong_eVis_other.png", otherDir);
    
    TH1D* h_signal_hCut_2Prong_gamma1ConvDist = (TH1D*)f_Root_Signal_Interaction->Get("hCut_2Prong_gamma1ConvDist");
    TH1D* h_background_hCut_2Prong_gamma1ConvDist = (TH1D*)f_Root_Background_Interaction->Get("hCut_2Prong_gamma1ConvDist");
    plotStacked(h_signal_hCut_2Prong_gamma1ConvDist, h_background_hCut_2Prong_gamma1ConvDist,"Leading Photon Conversion Distance [cm] (2Prong)", "hCut_2Prong_gamma1ConvDist.png", otherDir);
    
    TH1D* h_signal_hCut_2Prong_gamma2ConvDist = (TH1D*)f_Root_Signal_Interaction->Get("hCut_2Prong_gamma2ConvDist");
    TH1D* h_background_hCut_2Prong_gamma2ConvDist = (TH1D*)f_Root_Background_Interaction->Get("hCut_2Prong_gamma2ConvDist");
    plotStacked(h_signal_hCut_2Prong_gamma2ConvDist, h_background_hCut_2Prong_gamma2ConvDist,"Second Photon Conversion Distance [cm] (2Prong)", "hCut_2Prong_gamma2ConvDist.png", otherDir);
    
    TH1D* h_signal_hCut_2Prong_pi0invMass = (TH1D*)f_Root_Signal_Interaction->Get("hCut_2Prong_pi0invMass");
    TH1D* h_background_hCut_2Prong_pi0invMass = (TH1D*)f_Root_Background_Interaction->Get("hCut_2Prong_pi0invMass");
    plotStacked(h_signal_hCut_2Prong_pi0invMass, h_background_hCut_2Prong_pi0invMass,"Pi0 Invariant Mass (2Prong)", "hCut_2Prong_pi0invMass.png", otherDir);
    
    TH1D* h_signal_hCut_2Prong_neutrinoE = (TH1D*)f_Root_Signal_Interaction->Get("hCut_2Prong_neutrinoE");
    TH1D* h_background_hCut_2Prong_neutrinoE = (TH1D*)f_Root_Background_Interaction->Get("hCut_2Prong_neutrinoE");
    plotStacked(h_signal_hCut_2Prong_neutrinoE, h_background_hCut_2Prong_neutrinoE,"Neutrino Energy [GeV] (2Prong)", "hCut_2Prong_neutrinoE.png", otherDir);
        
    TH1D* h_signal_hCut_2Prong_UnusedE = (TH1D*)f_Root_Signal_Interaction->Get("hCut_2Prong_UnusedE");
    TH1D* h_background_hCut_2Prong_UnusedE = (TH1D*)f_Root_Background_Interaction->Get("hCut_2Prong_UnusedE");
    plotStacked(h_signal_hCut_2Prong_UnusedE, h_background_hCut_2Prong_UnusedE,"Unused Cluster Energy at the end of Reconstruction (2Prong)", "hCut_2Prong_UnusedE.png", otherDir);
    
    TH1D* h_signal_hCut_protonScore = (TH1D*)f_Root_Signal_Interaction->Get("hCut_protonScore");
    TH1D* h_background_hCut_protonScore = (TH1D*)f_Root_Background_Interaction->Get("hCut_protonScore");
    plotStacked(h_signal_hCut_protonScore, h_background_hCut_protonScore,"Proton Score", "hCut_protonScore.png", otherDir);
    
    TH1D* h_signal_hCut_protonScore_LLR = (TH1D*)f_Root_Signal_Interaction->Get("hCut_protonScore_LLR");
    TH1D* h_background_hCut_protonScore_LLR = (TH1D*)f_Root_Background_Interaction->Get("hCut_protonScore_LLR");
    plotStacked(h_signal_hCut_protonScore_LLR, h_background_hCut_protonScore_LLR,"Proton Score_LLR", "hCut_protonScore_LLR.png", otherDir);
    
    TH1D* h_signal_hCut_pionScore = (TH1D*)f_Root_Signal_Interaction->Get("hCut_pionScore");
    TH1D* h_background_hCut_pionScore = (TH1D*)f_Root_Background_Interaction->Get("hCut_pionScore");
    plotStacked(h_signal_hCut_pionScore, h_background_hCut_pionScore,"Pion Score", "hCut_pionScore.png", otherDir);
    
    TH1D* h_signal_hCut_pIDDiff = (TH1D*)f_Root_Signal_Interaction->Get("hCut_pIDDiff");
    TH1D* h_background_hCut_pIDDiff = (TH1D*)f_Root_Background_Interaction->Get("hCut_pIDDiff");
    plotStacked(h_signal_hCut_pIDDiff, h_background_hCut_pIDDiff,"Proton Score - Pion Score", "hCut_pIDDiff.png", otherDir);
    
    TH1D* h_signal_hCut_deltaInvMass = (TH1D*)f_Root_Signal_Interaction->Get("hCut_deltaInvMass");
    TH1D* h_background_hCut_deltaInvMass = (TH1D*)f_Root_Background_Interaction->Get("hCut_deltaInvMass");
    plotStacked(h_signal_hCut_deltaInvMass, h_background_hCut_deltaInvMass,"Delta+ Invariant Mass", "hCut_deltaInvMass.png", otherDir);
    
    
}

void Plotter::plotSignalBackground()
{
    // Files
    TFile* f_Root_Signal_Interaction = new TFile(rootDir_Interaction[0].c_str());
    TFile* f_Root_Background_Interaction = new TFile(rootDir_Interaction[1].c_str());
    TFile* f_Root_Signal_Muon = new TFile(rootDir_Muon[0].c_str());
    TFile* f_Root_Background_Muon = new TFile(rootDir_Muon[1].c_str());
    TFile* f_Root_Signal_Proton = new TFile(rootDir_Proton[0].c_str());
    TFile* f_Root_Background_Proton = new TFile(rootDir_Proton[1].c_str());
    TFile* f_Root_Signal_Pion = new TFile(rootDir_Pion[0].c_str());
    TFile* f_Root_Background_Pion = new TFile(rootDir_Pion[1].c_str());
    
    // Proton Plots
    TH1D* h_signal_protonScore = (TH1D*)f_Root_Signal_Proton->Get("protonScore");
    TH1D* h_background_protonScore = (TH1D*)f_Root_Background_Proton->Get("protonScore");
    plotStacked(h_signal_protonScore, h_background_protonScore,"Proton Score", "protonScore.png", otherDir);
    
    TH1D* h_signal_pionScore = (TH1D*)f_Root_Signal_Proton->Get("pionScore");
    TH1D* h_background_pionScore = (TH1D*)f_Root_Background_Proton->Get("pionScore");
    plotStacked(h_signal_pionScore, h_background_pionScore,"PionScore", "pionScore.png", otherDir);
        
    TH1D* h_signal_E_Proton = (TH1D*)f_Root_Signal_Proton->Get("E_reco");
    TH1D* h_background_E_Proton = (TH1D*)f_Root_Background_Proton->Get("E_reco");
    plotStacked(h_signal_E_Proton, h_background_E_Proton,"Proton Energy [MeV]", "E_Proton.png", otherDir);
    plotStackedLogScale(h_signal_E_Proton, h_background_E_Proton,"Proton Energy [MeV]", "E_Proton_Log.png", otherDir);
    
    TH1D* h_signal_P_Proton = (TH1D*)f_Root_Signal_Proton->Get("P_reco");
    TH1D* h_background_P_Proton = (TH1D*)f_Root_Background_Proton->Get("P_reco");
    plotStacked(h_signal_P_Proton, h_background_P_Proton,"Proton Momentum [MeV]", "P_Proton.png", otherDir);
    plotStackedLogScale(h_signal_P_Proton, h_background_P_Proton,"Proton Momentum [MeV]", "P_Proton_Log.png", otherDir);
    
    TH1D* h_signal_angleBeam_reco_Proton = (TH1D*)f_Root_Signal_Proton->Get("angleBeam_reco");
    TH1D* h_background_angleBeam_reco_Proton = (TH1D*)f_Root_Background_Proton->Get("angleBeam_reco");
    plotStacked(h_signal_angleBeam_reco_Proton, h_background_angleBeam_reco_Proton,"Proton Angle wrt Beam", "angleBeam_reco_Proton.png", otherDir);
    plotStackedLogScale(h_signal_angleBeam_reco_Proton, h_background_angleBeam_reco_Proton,"Proton Angle wrt Beam", "angleBeam_reco_Proton_Log.png", otherDir);
    
    // Pion Plots
    TH1D* h_signal_gamma1_ConvLength = (TH1D*)f_Root_Signal_Pion->Get("gamma1_ConvLength");
    TH1D* h_background_gamma1_ConvLength = (TH1D*)f_Root_Background_Pion->Get("gamma1_ConvLength");
    plotStacked(h_signal_gamma1_ConvLength, h_background_gamma1_ConvLength,"Leading Photon Conversion Length [cm]", "gamma1_ConvLength.png", otherDir);
    
    TH1D* h_signal_gamma2_ConvLength = (TH1D*)f_Root_Signal_Pion->Get("gamma2_ConvLength");
    TH1D* h_background_gamma2_ConvLength = (TH1D*)f_Root_Background_Pion->Get("gamma2_ConvLength");
    plotStacked(h_signal_gamma2_ConvLength, h_background_gamma2_ConvLength,"Second Photon Conversion Length [cm]", "gamma2_ConvLength.png", otherDir);
    
    TH1D* h_signal_gamma1_nClusters_All = (TH1D*)f_Root_Signal_Pion->Get("gamma1_nClusters_All");
    TH1D* h_background_gamma1_nClusters_All = (TH1D*)f_Root_Background_Pion->Get("gamma1_nClusters_All");
    plotStacked(h_signal_gamma1_nClusters_All, h_background_gamma1_nClusters_All,"Leading Photon N(Clusters)", "gamma1_nClusters_All.png", otherDir);
    
    TH1D* h_signal_gamma2_nClusters_All = (TH1D*)f_Root_Signal_Pion->Get("gamma2_nClusters_All");
    TH1D* h_background_gamma2_nClusters_All = (TH1D*)f_Root_Background_Pion->Get("gamma2_nClusters_All");
    plotStacked(h_signal_gamma2_nClusters_All, h_background_gamma2_nClusters_All,"Second Photon N(Clusters)", "gamma2_nClusters_All.png", otherDir);
    
    TH1D* h_signal_gamma1_nClusters_X = (TH1D*)f_Root_Signal_Pion->Get("gamma1_nClusters_X");
    TH1D* h_background_gamma1_nClusters_X = (TH1D*)f_Root_Background_Pion->Get("gamma1_nClusters_X");
    plotStacked(h_signal_gamma1_nClusters_X, h_background_gamma1_nClusters_X,"Number of X Clusters in the Leading Photon", "gamma1_nClusters_X.png", otherDir);
    
    TH1D* h_signal_gamma2_nClusters_X = (TH1D*)f_Root_Signal_Pion->Get("gamma2_nClusters_X");
    TH1D* h_background_gamma2_nClusters_X = (TH1D*)f_Root_Background_Pion->Get("gamma2_nClusters_X");
    plotStacked(h_signal_gamma2_nClusters_X, h_background_gamma2_nClusters_X,"Number of X Clusters in the Second Photon", "gamma2_nClusters_X.png", otherDir);
    
    TH1D* h_signal_gamma1_Energy= (TH1D*)f_Root_Signal_Pion->Get("gamma1_Energy");
    TH1D* h_background_gamma1_Energy= (TH1D*)f_Root_Background_Pion->Get("gamma1_Energy");
    plotStacked(h_signal_gamma1_Energy, h_background_gamma1_Energy,"Leading Photon Energy [GeV]", "gamma1_Energy.png", otherDir);
    
    TH1D* h_signal_gamma2_Energy= (TH1D*)f_Root_Signal_Pion->Get("gamma2_Energy");
    TH1D* h_background_gamma2_Energy= (TH1D*)f_Root_Background_Pion->Get("gamma2_Energy");
    plotStacked(h_signal_gamma2_Energy, h_background_gamma2_Energy,"Second Photon Energy [GeV]", "gamma2_Energy.png", otherDir);
    
    TH1D* h_signal_photonEnergy_Asymmetry= (TH1D*)f_Root_Signal_Pion->Get("photonEnergy_Asymmetry");
    TH1D* h_background_photonEnergy_Asymmetry= (TH1D*)f_Root_Background_Pion->Get("photonEnergy_Asymmetry");
    plotStacked(h_signal_photonEnergy_Asymmetry, h_background_photonEnergy_Asymmetry,"Photon Energy Asymmetry", "photonEnergy_Asymmetry.png", otherDir);
    
    TH1D* h_signal_E_Pion = (TH1D*)f_Root_Signal_Pion->Get("E_reco");
    TH1D* h_background_E_Pion = (TH1D*)f_Root_Background_Pion->Get("E_reco");
    plotStacked(h_signal_E_Pion, h_background_E_Pion,"Pion Energy [MeV]", "E_Pion.png", otherDir);
    
    TH1D* h_signal_P_Pion = (TH1D*)f_Root_Signal_Pion->Get("P_reco");
    TH1D* h_background_P_Pion = (TH1D*)f_Root_Background_Pion->Get("P_reco");
    plotStacked(h_signal_P_Pion, h_background_P_Pion,"Pion Momentum [MeV]", "P_Pion.png", otherDir);
    
    TH1D* h_signal_angleBeam_reco_Pion = (TH1D*)f_Root_Signal_Pion->Get("angleBeam_reco");
    TH1D* h_background_angleBeam_reco_Pion = (TH1D*)f_Root_Background_Pion->Get("angleBeam_reco");
    plotStacked(h_signal_angleBeam_reco_Pion, h_background_angleBeam_reco_Pion,"Pion Angle wrt Beam", "angleBeam_reco_Pion.png", otherDir);
    
    TH1D* h_signal_invMass = (TH1D*)f_Root_Signal_Pion->Get("invMass");
    TH1D* h_background_invMass = (TH1D*)f_Root_Background_Pion->Get("invMass");
    plotStacked(h_signal_invMass, h_background_invMass,"Pion Invariant Mass [MeV]", "invMass.png", otherDir);
     
    // Muon Plots
    TH1D* h_signal_E_Muon = (TH1D*)f_Root_Signal_Muon->Get("E_reco");
    TH1D* h_background_E_Muon = (TH1D*)f_Root_Background_Muon->Get("E_reco");
    plotStacked(h_signal_E_Muon, h_background_E_Muon,"Muon Energy [GeV]", "E_Muon.png", otherDir);
    
    TH1D* h_signal_P_Muon = (TH1D*)f_Root_Signal_Muon->Get("P_reco");
    TH1D* h_background_P_Muon = (TH1D*)f_Root_Background_Muon->Get("P_reco");
    plotStacked(h_signal_P_Muon, h_background_P_Muon,"Muon Momentum [GeV]", "P_Muon.png", otherDir);
    
    TH1D* h_signal_angleBeam_reco_Muon = (TH1D*)f_Root_Signal_Muon->Get("angleBeam_reco");
    TH1D* h_background_angleBeam_reco_Muon = (TH1D*)f_Root_Background_Muon->Get("angleBeam_reco");
    plotStacked(h_signal_angleBeam_reco_Muon, h_background_angleBeam_reco_Muon,"Muon Angle wrt Beam", "angleBeam_reco_Muon.png", otherDir);
    
    // Interaction Plots 
    TH1D* h_signal_q2_mc = (TH1D*)f_Root_Signal_Interaction->Get("q2_mc");
    TH1D* h_background_q2_mc = (TH1D*)f_Root_Background_Interaction->Get("q2_mc");
    plotStacked(h_signal_q2_mc, h_background_q2_mc,"True Q^{2} [GeV^{2}]", "Q_Square_truth.png", otherDir);
    
    TH1D* h_signal_q2_reco = (TH1D*)f_Root_Signal_Interaction->Get("q2_reco");
    TH1D* h_background_q2_reco = (TH1D*)f_Root_Background_Interaction->Get("q2_reco");
    plotStacked(h_signal_q2_reco, h_background_q2_reco,"Reconstructed Q^{2} [GeV^{2}]", "Q_Square_reco.png", otherDir);
    
    TH1D* h_signal_w_mc = (TH1D*)f_Root_Signal_Interaction->Get("w_mc");
    TH1D* h_background_w_mc = (TH1D*)f_Root_Background_Interaction->Get("w_mc");
    plotStacked(h_signal_w_mc, h_background_w_mc,"True W [GeV]", "W_truth.png", otherDir);
    
    TH1D* h_signal_w_reco = (TH1D*)f_Root_Signal_Interaction->Get("w_reco");
    TH1D* h_background_w_reco = (TH1D*)f_Root_Background_Interaction->Get("w_reco");
    plotStacked(h_signal_w_reco, h_background_w_reco,"Reconstructed W [GeV]", "W_reco.png", otherDir);
    
    TH1D* h_signal_wSq_reco = (TH1D*)f_Root_Signal_Interaction->Get("wSq_reco");
    TH1D* h_background_wSq_reco = (TH1D*)f_Root_Background_Interaction->Get("wSq_reco");
    plotStacked(h_signal_wSq_reco, h_background_wSq_reco,"Reconstructed W^{2} [GeV^{2}]", "wSq_reco.png", otherDir);
    
    TH1D* h_signal_beamEnergy_mc = (TH1D*)f_Root_Signal_Interaction->Get("beamEnergy_mc");
    TH1D* h_background_beamEnergy_mc = (TH1D*)f_Root_Background_Interaction->Get("beamEnergy_mc");
    plotStacked(h_signal_beamEnergy_mc, h_background_beamEnergy_mc,"True Beam Energy [GeV]", "BeamEnergy_truth.png", otherDir);
    
    TH1D* h_signal_beamEnergy_reco = (TH1D*)f_Root_Signal_Interaction->Get("beamEnergy_reco");
    TH1D* h_background_beamEnergy_reco = (TH1D*)f_Root_Background_Interaction->Get("beamEnergy_reco");
    plotStacked(h_signal_beamEnergy_reco, h_background_beamEnergy_reco,"Reconstructed Beam Energy [GeV]", "BeamEnergy_reco.png", otherDir);
      
    TH1D* h_signal_beamEnergyCal_reco = (TH1D*)f_Root_Signal_Interaction->Get("beamEnergyCal_reco");
    TH1D* h_background_beamEnergyCal_reco = (TH1D*)f_Root_Background_Interaction->Get("beamEnergyCal_reco");
    plotStacked(h_signal_beamEnergyCal_reco, h_background_beamEnergyCal_reco,"Reconstructed Calorimetric Beam Energy [GeV]", "BeamEnergyCal_reco.png", otherDir);
    
    TH1D* h_signal_total_E = (TH1D*)f_Root_Signal_Interaction->Get("total_E");
    TH1D* h_background_total_E = (TH1D*)f_Root_Background_Interaction->Get("total_E");
    plotStacked(h_signal_total_E, h_background_total_E,"Total FS Particle Energy [GeV]", "total_E.png", otherDir);
    
    
    delete f_Root_Signal_Interaction;
    delete f_Root_Background_Interaction;
    delete f_Root_Signal_Muon;
    delete f_Root_Background_Muon;
    delete f_Root_Signal_Proton;
    delete f_Root_Background_Proton;
    delete f_Root_Signal_Pion;
    delete f_Root_Background_Pion;
}

 
void Plotter::plot_mc_w_Stacked()
{
    string rootDir = rootDir_Interaction[branchInd];
    string plotDir = plotDir_Interaction[branchInd];
    
    cout<<"\nPlottting Stacked mc_w"<<endl;
    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());
    TCanvas* c1 = new TCanvas();
    THStack *hs = new THStack("hs","TRUE Signal Events (MINOS Matched)");
    TLegend *legend = new TLegend(0.7,0.8,0.9,0.9);  
    
    TH1D* h_mc_w_DIS = (TH1D*)f_Root->Get("mc_w_DIS");
    h_mc_w_DIS->SetFillColor(kRed);
    h_mc_w_DIS->SetMarkerStyle(21);
    h_mc_w_DIS->SetMarkerColor(kRed);
    
    TH1D* h_mc_w_RES = (TH1D*)f_Root->Get("mc_w_RES");
    h_mc_w_RES->SetFillColor(kBlue);
    h_mc_w_RES->SetMarkerStyle(21);
    h_mc_w_RES->SetMarkerColor(kBlue);
      
    legend->AddEntry(h_mc_w_DIS , "DIS", "f");
    legend->AddEntry(h_mc_w_RES, "RES", "f");

    
    hs->Add(h_mc_w_DIS);
    hs->Add(h_mc_w_RES);
    hs->Draw();
    hs->GetXaxis()->SetTitle("mc_w [GeV/c^{2}]");
    hs->GetYaxis()->SetTitle("N(Events)");
    
    legend->Draw();
    
    c1->Print(Form("%s%s",plotDir.c_str(),"mc_w.png"),"png");
    
    delete f_Root;
    delete hs;
    delete legend;
    delete c1;
}

void Plotter::plot_final_mc_w_Stacked()
{
    string rootDir = rootDir_Interaction[branchInd];
    string plotDir = plotDir_Interaction[branchInd];
    
    cout<<"\nPlottting Stacked final_mc_w"<<endl;
    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());
    TCanvas* c1 = new TCanvas();
    THStack *hs = new THStack("hs","TRUE Signal Events (MINOS Matched)");
    TLegend *legend = new TLegend(0.7,0.8,0.9,0.9);  
    
    TH1D* h_mc_w_DIS = (TH1D*)f_Root->Get("final_mc_w_DIS");
    h_mc_w_DIS->SetFillColor(kRed);
    h_mc_w_DIS->SetMarkerStyle(21);
    h_mc_w_DIS->SetMarkerColor(kRed);
    
    TH1D* h_mc_w_RES = (TH1D*)f_Root->Get("final_mc_w_RES");
    h_mc_w_RES->SetFillColor(kBlue);
    h_mc_w_RES->SetMarkerStyle(21);
    h_mc_w_RES->SetMarkerColor(kBlue);
    
    legend->AddEntry(h_mc_w_DIS , "DIS", "f");
    legend->AddEntry(h_mc_w_RES, "RES", "f");
    
    
    hs->Add(h_mc_w_DIS);
    hs->Add(h_mc_w_RES);
    hs->Draw();
    hs->GetXaxis()->SetTitle("mc_w [GeV/c^{2}]");
    hs->GetYaxis()->SetTitle("N(Events)");
    
    legend->Draw();
    
    c1->Print(Form("%s%s",plotDir.c_str(),"final_mc_w.png"),"png");
    
    delete f_Root;
    delete hs;
    delete legend;
    delete c1;
}





#endif
