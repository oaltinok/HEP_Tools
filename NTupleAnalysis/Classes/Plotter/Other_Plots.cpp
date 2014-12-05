#ifndef Other_Plots_cpp
#define Other_Plots_cpp

#include "Plotter.h"

using namespace std;

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
    plotStacked(h_signal_protonScore, h_background_protonScore,"Proton Score", "Stacked_protonScore.png", otherDir);
    
    TH1D* h_signal_pionScore = (TH1D*)f_Root_Signal_Proton->Get("pionScore");
    TH1D* h_background_pionScore = (TH1D*)f_Root_Background_Proton->Get("pionScore");
    plotStacked(h_signal_pionScore, h_background_pionScore,"PionScore", "Stacked_pionScore.png", otherDir);
    
    TH1D* h_signal_E_Proton = (TH1D*)f_Root_Signal_Proton->Get("E_reco");
    TH1D* h_background_E_Proton = (TH1D*)f_Root_Background_Proton->Get("E_reco");
    plotStacked(h_signal_E_Proton, h_background_E_Proton,"Proton Energy [MeV]", "Stacked_E_Proton.png", otherDir);
    
    // Pion Plots
    TH1D* h_signal_photonConvLength = (TH1D*)f_Root_Signal_Pion->Get("photonConvLength");
    TH1D* h_background_photonConvLength = (TH1D*)f_Root_Background_Pion->Get("photonConvLength");
    plotStacked(h_signal_photonConvLength, h_background_photonConvLength,"Photon Conversion Length [cm]", "Stacked_photonConvLength.png", otherDir);
    
    TH1D* h_signal_E_Pion = (TH1D*)f_Root_Signal_Pion->Get("E_reco");
    TH1D* h_background_E_Pion = (TH1D*)f_Root_Background_Pion->Get("E_reco");
    plotStacked(h_signal_E_Pion, h_background_E_Pion,"Pion Energy [MeV]", "Stacked_E_Pion.png", otherDir);
    
    TH1D* h_signal_invMass = (TH1D*)f_Root_Signal_Pion->Get("invMass");
    TH1D* h_background_invMass = (TH1D*)f_Root_Background_Pion->Get("invMass");
    plotStacked(h_signal_invMass, h_background_invMass,"Pion Invariant Mass [MeV]", "Stacked_invMass.png", otherDir);
     
    // Muon Plots
    TH1D* h_signal_E_Muon = (TH1D*)f_Root_Signal_Muon->Get("E_reco");
    TH1D* h_background_E_Muon = (TH1D*)f_Root_Background_Muon->Get("E_reco");
    plotStacked(h_signal_E_Muon, h_background_E_Muon,"Muon Energy [GeV]", "Stacked_E_Muon.png", otherDir);
    
    // Interaction Plots
    TH1D* h_signal_int_channel = (TH1D*)f_Root_Signal_Interaction->Get("int_channel");
    TH1D* h_background_int_channel = (TH1D*)f_Root_Background_Interaction->Get("int_channel");
    plotStacked(h_signal_int_channel, h_background_int_channel,"Channel ID: (1 = QE, 2 = Resonant, 3 = DIS, 4 = Coh pi)", "Stacked_int_channel.png", otherDir);
    
    TH1D* h_signal_nProngs_hist = (TH1D*)f_Root_Signal_Interaction->Get("nProngs_hist");
    TH1D* h_background_nProngs_hist = (TH1D*)f_Root_Background_Interaction->Get("nProngs_hist");
    plotStacked(h_signal_nProngs_hist, h_background_nProngs_hist,"N(Prongs)", "Stacked_nProngs.png", otherDir);
    
    TH1D* h_signal_nProngs_hist_initial = (TH1D*)f_Root_Signal_Interaction->Get("nProngs_hist_initial");
    TH1D* h_background_nProngs_hist_initial = (TH1D*)f_Root_Background_Interaction->Get("nProngs_hist_initial");
    plotStacked(h_signal_nProngs_hist_initial, h_background_nProngs_hist_initial,"N(Prongs)", "Stacked_nProngs_initial.png", otherDir);
    
    TH1D* h_signal_q2_mc = (TH1D*)f_Root_Signal_Interaction->Get("q2_mc");
    TH1D* h_background_q2_mc = (TH1D*)f_Root_Background_Interaction->Get("q2_mc");
    plotStacked(h_signal_q2_mc, h_background_q2_mc,"True Q^{2} [GeV^{2}]", "Stacked_Q_Square_truth.png", otherDir);
    
    TH1D* h_signal_q2_reco = (TH1D*)f_Root_Signal_Interaction->Get("q2_reco");
    TH1D* h_background_q2_reco = (TH1D*)f_Root_Background_Interaction->Get("q2_reco");
    plotStacked(h_signal_q2_reco, h_background_q2_reco,"Reconstructed Q^{2} [GeV^{2}]", "Stacked_Q_Square_reco.png", otherDir);
    
    TH1D* h_signal_w2_mc = (TH1D*)f_Root_Signal_Interaction->Get("w2_mc");
    TH1D* h_background_w2_mc = (TH1D*)f_Root_Background_Interaction->Get("w2_mc");
    plotStacked(h_signal_w2_mc, h_background_w2_mc,"True W [GeV]", "Stacked_W_truth.png", otherDir);
    
    TH1D* h_signal_w2_reco = (TH1D*)f_Root_Signal_Interaction->Get("w2_reco");
    TH1D* h_background_w2_reco = (TH1D*)f_Root_Background_Interaction->Get("w2_reco");
    plotStacked(h_signal_w2_reco, h_background_w2_reco,"Reconstructed W [GeV]", "Stacked_W_reco.png", otherDir);
    
    TH1D* h_signal_beamEnergy_mc = (TH1D*)f_Root_Signal_Interaction->Get("beamEnergy_mc");
    TH1D* h_background_beamEnergy_mc = (TH1D*)f_Root_Background_Interaction->Get("beamEnergy_mc");
    plotStacked(h_signal_beamEnergy_mc, h_background_beamEnergy_mc,"True Beam Energy [GeV]", "Stacked_BeamEnergy_truth.png", otherDir);
    
    TH1D* h_signal_beamEnergy_reco = (TH1D*)f_Root_Signal_Interaction->Get("beamEnergy_reco");
    TH1D* h_background_beamEnergy_reco = (TH1D*)f_Root_Background_Interaction->Get("beamEnergy_reco");
    plotStacked(h_signal_beamEnergy_reco, h_background_beamEnergy_reco,"Reconstructed Beam Energy [GeV]", "Stacked_BeamEnergy_reco.png", otherDir);
      
    TH1D* h_signal_beamEnergyCal_reco = (TH1D*)f_Root_Signal_Interaction->Get("beamEnergyCal_reco");
    TH1D* h_background_beamEnergyCal_reco = (TH1D*)f_Root_Background_Interaction->Get("beamEnergyCal_reco");
    plotStacked(h_signal_beamEnergyCal_reco, h_background_beamEnergyCal_reco,"Reconstructed Calorimetric Beam Energy [GeV]", "Stacked_BeamEnergyCal_reco.png", otherDir);
    
    
    delete f_Root_Signal_Interaction;
    delete f_Root_Background_Interaction;
    delete f_Root_Signal_Muon;
    delete f_Root_Background_Muon;
    delete f_Root_Signal_Proton;
    delete f_Root_Background_Proton;
    delete f_Root_Signal_Pion;
    delete f_Root_Background_Pion;
}


void Plotter::plotDebuggingPlots()
{
    string rootDir = rootDir_Interaction[branchInd];
    string plotDir = plotDir_Interaction[branchInd];
    
    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());
    
    TH1D* h_w2fail_q2 = (TH1D*)f_Root->Get("w2fail_q2");
    plot1D_Hist(h_w2fail_q2,"w2fail_q2.png",plotDir);
    
    TH1D* h_w2fail_Enu = (TH1D*)f_Root->Get("w2fail_Enu");
    plot1D_Hist(h_w2fail_Enu,"w2fail_Enu.png",plotDir);
    
    TH1D* h_w2fail_muon_E= (TH1D*)f_Root->Get("w2fail_muon_E");
    plot1D_Hist(h_w2fail_muon_E,"w2fail_muon_E.png",plotDir);
    
    TH1D* h_w2fail_muon_Pz= (TH1D*)f_Root->Get("w2fail_muon_Pz");
    plot1D_Hist(h_w2fail_muon_Pz,"w2fail_muon_Pz.png",plotDir);
    
    TH1D* h_w2fail_pion_E= (TH1D*)f_Root->Get("w2fail_pion_E");
    plot1D_Hist(h_w2fail_pion_E,"w2fail_pion_E.png",plotDir);
    
    TH1D* h_w2fail_proton_KE= (TH1D*)f_Root->Get("w2fail_proton_KE");
    plot1D_Hist(h_w2fail_proton_KE,"w2fail_proton_KE.png",plotDir);
    
    TH1D* h_w2fail_term1= (TH1D*)f_Root->Get("w2fail_term1");
    plot1D_Hist(h_w2fail_term1,"w2fail_term1.png",plotDir);
    
    TH1D* h_w2fail_w2= (TH1D*)f_Root->Get("w2fail_w2");
    plot1D_Hist(h_w2fail_w2,"w2fail_w2.png",plotDir);
    
    TH2D* h_w2fail_term1_q2= (TH2D*)f_Root->Get("w2fail_term1_q2");
    plot2D_Hist(h_w2fail_term1_q2,"w2fail_term1_q2.png",plotDir);
    
    delete f_Root;
    
}    

void Plotter::plotPID()
{
    pID_proton();
    pID_pion();
    plot_2D_pID();

}

void Plotter::pID_pion()
{
    string rootDir = rootDir_Interaction[branchInd];
    string plotDir = plotDir_Interaction[branchInd];
    
    cout<<"\nPlottting pID for Pion"<<endl;
    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());
    TCanvas* c1 = new TCanvas();
    THStack *hs = new THStack("hs","Pion Score");
    TLegend *legend = new TLegend(0.7,0.8,0.9,0.9); 
    
    TH1D* h_pID_piminus = (TH1D*)f_Root->Get("pID_piminus_pionScore");
    h_pID_piminus->SetFillColor(kYellow);
    h_pID_piminus->SetMarkerStyle(21);
    h_pID_piminus->SetMarkerColor(kYellow);
    
    TH1D* h_pID_piplus = (TH1D*)f_Root->Get("pID_piplus_pionScore");
    h_pID_piplus->SetFillColor(kBlue);
    h_pID_piplus->SetMarkerStyle(21);
    h_pID_piplus->SetMarkerColor(kBlue);
    
    TH1D* h_pID_proton = (TH1D*)f_Root->Get("pID_proton_pionScore");
    h_pID_proton->SetFillColor(kGreen);
    h_pID_proton->SetMarkerStyle(21);
    h_pID_proton->SetMarkerColor(kGreen);
    
    legend->AddEntry(h_pID_proton, "Proton", "f");
    legend->AddEntry(h_pID_piplus, "Pi-Plus", "f");
    legend->AddEntry(h_pID_piminus, "Pi-Minus", "f");
    
    
    hs->Add(h_pID_piminus);
    hs->Add(h_pID_piplus);
    hs->Add(h_pID_proton);
    hs->Draw();
    hs->GetXaxis()->SetTitle("Pion Score");
    hs->GetYaxis()->SetTitle("N(Events)");
    
    legend->Draw();
    
    c1->Print(Form("%s%s",plotDir.c_str(),"stacked_pID_pionScore.png"),"png");
    
    delete f_Root;
    delete hs;
    delete legend;
    delete c1;
    
}

void Plotter::pID_proton()
{
    string rootDir = rootDir_Interaction[branchInd];
    string plotDir = plotDir_Interaction[branchInd];
    
    cout<<"\nPlottting pID for Proton"<<endl;
    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());
    TCanvas* c1 = new TCanvas();
    THStack *hs = new THStack("hs","Proton Score");
    TLegend *legend = new TLegend(0.7,0.8,0.9,0.9); 
    
    TH1D* h_pID_piminus = (TH1D*)f_Root->Get("pID_piminus_protonScore");
    h_pID_piminus->SetFillColor(kYellow);
    h_pID_piminus->SetMarkerStyle(21);
    h_pID_piminus->SetMarkerColor(kYellow);
    
    TH1D* h_pID_piplus = (TH1D*)f_Root->Get("pID_piplus_protonScore");
    h_pID_piplus->SetFillColor(kBlue);
    h_pID_piplus->SetMarkerStyle(21);
    h_pID_piplus->SetMarkerColor(kBlue);
    
    TH1D* h_pID_proton = (TH1D*)f_Root->Get("pID_proton_protonScore");
    h_pID_proton->SetFillColor(kGreen);
    h_pID_proton->SetMarkerStyle(21);
    h_pID_proton->SetMarkerColor(kGreen);
    
    legend->AddEntry(h_pID_proton, "Proton", "f");
    legend->AddEntry(h_pID_piplus, "Pi-Plus", "f");
    legend->AddEntry(h_pID_piminus, "Pi-Minus", "f");
    
    
    hs->Add(h_pID_piminus);
    hs->Add(h_pID_piplus);
    hs->Add(h_pID_proton);
    hs->Draw();
    hs->GetXaxis()->SetTitle("Proton Score");
    hs->GetYaxis()->SetTitle("N(Events)");
    
    legend->Draw();
    
    c1->Print(Form("%s%s",plotDir.c_str(),"stacked_pID_protonScore.png"),"png");
    
    delete f_Root;
    delete hs;
    delete legend;
    delete c1;
    
}

void Plotter::plot_2D_pID()
{
    string rootDir = rootDir_Interaction[branchInd];
    string plotDir = plotDir_Interaction[branchInd];
    
    cout<<"\nPlottting 2D pID Plots"<<endl;
    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());
    
    TH2D* h_pID_proton_pionScore_protonScore = (TH2D*)f_Root->Get("pID_proton_pionScore_protonScore");
    plot2D_Hist(h_pID_proton_pionScore_protonScore,"pID_proton_pionScore_protonScore.png",plotDir);
    
    TH2D* h_pID_piminus_pionScore_protonScore = (TH2D*)f_Root->Get("pID_piminus_pionScore_protonScore");
    plot2D_Hist(h_pID_piminus_pionScore_protonScore,"pID_piminus_pionScore_protonScore.png",plotDir);
    
    TH2D* h_pID_piplus_pionScore_protonScore = (TH2D*)f_Root->Get("pID_piplus_pionScore_protonScore");
    plot2D_Hist(h_pID_piplus_pionScore_protonScore,"pID_piplus_pionScore_protonScore.png",plotDir);
    
    
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
    
    c1->Print(Form("%s%s",plotDir.c_str(),"stacked_mc_w.png"),"png");
    
    delete f_Root;
    delete hs;
    delete legend;
    delete c1;
}





#endif
