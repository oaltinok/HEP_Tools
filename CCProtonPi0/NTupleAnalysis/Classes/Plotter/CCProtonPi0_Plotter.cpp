/*
    See CCProtonPi0_Plotter.h header for Class Information
*/

#ifndef CCProtonPi0_Plotter_cpp
#define CCProtonPi0_Plotter_cpp

#include "CCProtonPi0_Plotter.h"

using namespace std;

void CCProtonPi0_Plotter::plotHistograms()
{
    if (isSignalvsBackground){
//         plotSignalBackground();
//         plotCutHistograms();
    }else{
//         plotInteraction();
        
//         plotMuon();
        
//         plotProton();
        
        plotPion();
        
//        plotPID();
        
//         plot_mc_w_Stacked();
//         plot_final_mc_w_Stacked();
//             plotMichel();
    }
   
    
}

void CCProtonPi0_Plotter::plotSignalRatio(TH1D* h_signal, TH1D* h_background, string fileName, string plotDir, bool isReversed) 
{
    TH1D* h_ratio = new TH1D;
    
    if (isReversed){
        h_background->Copy(*h_ratio);
        
        h_ratio->Divide(h_signal);
        h_ratio->GetYaxis()->SetTitle("Background / Signal");
    }else{
        h_signal->Copy(*h_ratio);
        
        h_ratio->Divide(h_background);
        h_ratio->GetYaxis()->SetTitle("Signal / Background");
    }
    
    TCanvas* c1 = new TCanvas();
    h_ratio->SetLineColor(kBlack);
    h_ratio->SetLineWidth(2);
    
    h_ratio->Draw();
    gPad->Update();
     
    c1->Print(Form("%s%s%s",plotDir.c_str(),"Ratio_",fileName.c_str()),"png");
    
    delete c1;
    delete h_ratio;
}


void CCProtonPi0_Plotter::inform(string rootDir, string plotDir)
{
    cout<<"------------ Plotting ------------"<<endl;
    cout<<"Input File: "<<rootDir<<endl;
    cout<<"Output Folder: "<<plotDir<<endl;
}

void CCProtonPi0_Plotter::plotStackedLogScale(TH1D* h_signal, TH1D* h_background, string plotName, string fileName, string plotDir)
{
    TH1D* h_signalRatio = new TH1D;
    TH1D* h_backgroundRatio = new TH1D;
    
    h_signal->Copy(*h_signalRatio);
    h_background->Copy(*h_backgroundRatio);
    
    TCanvas *c1 = new TCanvas();
    THStack *hs = new THStack("hs",plotName.c_str());
    TLegend *legend = new TLegend(0.7,0.8,0.9,0.9);
    
    c1->SetLogy();
    
    h_signal->SetFillColor(kGreen);
    h_signal->SetLineColor(kGreen);
    h_signal->SetMarkerStyle(21);
    h_signal->SetMarkerColor(kGreen);
    
    h_background->SetFillColor(kRed);
    h_background->SetLineColor(kRed);
    h_background->SetMarkerStyle(21);
    h_background->SetMarkerColor(kRed);
    
    legend->AddEntry(h_signal, "Signal", "f");
    legend->AddEntry(h_background, "Background", "f");
    
    hs->Add(h_background);
    hs->Add(h_signal);
    hs->Draw();
    hs->GetXaxis()->SetTitle(plotName.c_str());
    hs->GetYaxis()->SetTitle("N(Events)");
    
    legend->Draw();
    
    c1->Print(Form("%s%s",plotDir.c_str(),fileName.c_str()),"png");
    
    delete hs;
    delete legend;
    delete c1;
    
    // Plot Signal Ratio
    plotSignalRatio(h_signalRatio,h_backgroundRatio, fileName, plotDir);
    
    
    delete h_signalRatio;
    delete h_backgroundRatio;
    
}

void CCProtonPi0_Plotter::plotStacked(TH1D* h_signal, TH1D* h_background, 
                            string plotName, string fileName, string plotDir, 
                            string signal_label, string background_label,
                            bool isRatioReversed)
{
    TH1D* h_signalRatio = new TH1D;
    TH1D* h_backgroundRatio = new TH1D;
    
    h_signal->Copy(*h_signalRatio);
    h_background->Copy(*h_backgroundRatio);
    
    TCanvas *c1 = new TCanvas();
    THStack *hs = new THStack("hs",plotName.c_str());
    TLegend *legend = new TLegend(0.7,0.8,0.9,0.9);
    
    
    h_signal->SetFillColor(kGreen);
    h_signal->SetLineColor(kGreen);
    h_signal->SetMarkerStyle(21);
    h_signal->SetMarkerColor(kGreen);
    
    h_background->SetFillColor(kRed);
    h_background->SetLineColor(kRed);
    h_background->SetMarkerStyle(21);
    h_background->SetMarkerColor(kRed);
    
    legend->AddEntry(h_signal, signal_label.c_str(), "f");
    legend->AddEntry(h_background, background_label.c_str(), "f");
 
    
    hs->Add(h_background);
    hs->Add(h_signal);
    hs->Draw();
    hs->GetXaxis()->SetTitle(plotName.c_str());
    hs->GetYaxis()->SetTitle("N(Events)");
    
    legend->Draw();
  
    c1->Print(Form("%s%s",plotDir.c_str(),fileName.c_str()),"png");
    
    delete hs;
    delete legend;
    delete c1;
    
    // Plot Signal Ratio
    plotSignalRatio(h_signalRatio,h_backgroundRatio,fileName, plotDir, isRatioReversed);

    delete h_signalRatio;
    delete h_backgroundRatio;
}

void CCProtonPi0_Plotter::plot1D_HistLogScale(TH1D* hist1D, string fileName, string plotDir)
{
    TCanvas* c1 = new TCanvas();
    c1->SetLogy();
    
    hist1D->SetLineColor(kRed);
    hist1D->SetLineWidth(3);
    hist1D->SetFillColor(kRed);
    hist1D->SetFillStyle(3010);
    
    hist1D->Draw();
    gPad->Update();
    
    // Statistics Box
    //     TPaveStats *st = (TPaveStats*)hist1D->FindObject("stats");
    
    c1->Print(Form("%s%s",plotDir.c_str(),fileName.c_str()),"png");
    delete c1;
    
}


void CCProtonPi0_Plotter::plot1D_Hist(TH1D* hist1D, string fileName, string plotDir)
{
    TCanvas* c1 = new TCanvas();
    hist1D->SetLineColor(kRed);
    hist1D->SetLineWidth(3);
    hist1D->SetFillColor(kRed);
    hist1D->SetFillStyle(3010);
    
    hist1D->Draw();
    gPad->Update();
    
    // Statistics Box
//     TPaveStats *st = (TPaveStats*)hist1D->FindObject("stats");
    
    c1->Print(Form("%s%s",plotDir.c_str(),fileName.c_str()),"png");
    delete c1;
    
}

void CCProtonPi0_Plotter::plot2D_Hist(TH2D* hist2D, string fileName, string plotDir)
{
    // Canvas
    Double_t w = 800; 
    Double_t h = 800;
    TCanvas* c1 = new TCanvas("c","c",w,h);
    c1->SetWindowSize(w,h);
    
    // Pad
    TPad *p = new TPad("p","p",0.05,0.05,0.95,0.95);
    p->Draw();
    
    p->cd();
    hist2D->GetYaxis()->SetTitleOffset(1.8);
    hist2D->Draw("colz");
    gPad->Update();
    
    // Statistics Box
    TPaveStats *st = (TPaveStats*)hist2D->FindObject("stats");
    st->SetOptStat(1000000110);
    st->SetX1NDC(0.1); 
    st->SetX2NDC(0.3); 
    st->SetY1NDC(0.8); 
    st->SetY2NDC(0.9); 
   
    c1->Print(Form("%s%s",plotDir.c_str(),fileName.c_str()),"png");
    
    delete p;
    delete c1;
}

//------------------------------------------------------------------------------
// Analysis Mode
//      1) Signal Events
//      2) Background Events
//      Other) All Events
//------------------------------------------------------------------------------
CCProtonPi0_Plotter::CCProtonPi0_Plotter(int nMode)
{
    isSignalvsBackground = false;
    
    if ( nMode == 1) {
        cout<<"----------------------------------------------------------------------"<<endl;
        cout<<"Plot Mode: Signal - Only Signal Events will be Plotted"<<endl;
        branchInd = 0;
    }else if ( nMode == 2){
        cout<<"----------------------------------------------------------------------"<<endl;
        cout<<"Plot Mode: Background - Only Background Events will be Plotted"<<endl;
        branchInd = 1;
    }else if (nMode == 3){
        cout<<"----------------------------------------------------------------------"<<endl;
        cout<<"Plot Mode: All - All Events will be Plotted"<<endl;
        branchInd = 2;
    }else{
        cout<<"----------------------------------------------------------------------"<<endl;
        cout<<"Plot Mode: Signal vs Background"<<endl;
        isSignalvsBackground = true;
    }
    
    cout<<"----------------------------------------------------------------------"<<endl;
    
    setFolders();
}

void CCProtonPi0_Plotter::setFolders()
{
    for ( int i = 0 ; i < nBranches; i++){
        
        // Set Branch: 0 = Signal, 1 = Background, Other = All Events
        if ( i == 0 ) branchDir = Folder_List::signal;
        else if ( i == 1) branchDir = Folder_List::background;
        else branchDir = Folder_List::allEvents;
        
        rootDir_Interaction[i] = Folder_List::output + Folder_List::rootOut + branchDir + "Interaction.root";
        plotDir_Interaction[i] = Folder_List::output + Folder_List::plotOut + branchDir + "Interaction/";
        
        rootDir_PID[i] = Folder_List::output + Folder_List::rootOut + branchDir + "PIDStatistics.root";
        plotDir_PID[i] = Folder_List::output + Folder_List::plotOut + branchDir + "PIDStatistics/";

        rootDir_Muon[i] = Folder_List::output + Folder_List::rootOut + branchDir + "Muon.root";
        plotDir_Muon[i] = Folder_List::output + Folder_List::plotOut + branchDir + "Muon/";

        rootDir_Proton[i] = Folder_List::output + Folder_List::rootOut + branchDir + "Proton.root";
        plotDir_Proton[i] = Folder_List::output + Folder_List::plotOut + branchDir + "Proton/";

        rootDir_Pion[i] = Folder_List::output + Folder_List::rootOut + branchDir + "Pion.root";
        plotDir_Pion[i] = Folder_List::output + Folder_List::plotOut + branchDir + "Pion/";
        
    }
    
    otherDir = Folder_List::output + Folder_List::plotOut + Folder_List::other;
}

void CCProtonPi0_Plotter::plotInteraction()
{
    string rootDir = rootDir_Interaction[branchInd];
    string plotDir = plotDir_Interaction[branchInd];
    
    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());

    // Signal only Plots
    if( branchInd == 0){
        TH1D* h_status_Pi0 = (TH1D*)f_Root->Get("status_Pi0");
        plot1D_Hist(h_status_Pi0,"status_Pi0.png",plotDir);
        
        TH1D* h_status_Pi0_Mother = (TH1D*)f_Root->Get("status_Pi0_Mother");
        plot1D_Hist(h_status_Pi0_Mother,"status_Pi0_Mother.png",plotDir);
        
        TH1D* h_status_Pi0_GrandMother = (TH1D*)f_Root->Get("status_Pi0_GrandMother");
        plot1D_Hist(h_status_Pi0_GrandMother,"status_Pi0_GrandMother.png",plotDir);
    }
     
    TH1D* h_beamEnergy_mc= (TH1D*)f_Root->Get("beamEnergy_mc");
    plot1D_Hist(h_beamEnergy_mc,"beamEnergy_mc.png",plotDir);
    
    TH1D* h_q2_mc= (TH1D*)f_Root->Get("q2_mc");
    plot1D_Hist(h_q2_mc,"q2_mc.png",plotDir);
    
    TH1D* h_w_mc= (TH1D*)f_Root->Get("w_mc");
    plot1D_Hist(h_w_mc,"w_mc.png",plotDir);
    
    TH1D* h_int_channel= (TH1D*)f_Root->Get("int_channel");
    plot1D_Hist(h_int_channel,"int_channel.png",plotDir);
    
    TH1D* h_vertex_z_true= (TH1D*)f_Root->Get("vertex_z_true");
    plot1D_Hist(h_vertex_z_true,"vertex_z_true.png",plotDir);
    
    TH2D* h_vertex_x_y_true= (TH2D*)f_Root->Get("vertex_x_y_true");
    plot2D_Hist(h_vertex_x_y_true,"vertex_x_y_true.png",plotDir);

    TH1D* h_total_E = (TH1D*)f_Root->Get("total_E");
    plot1D_Hist(h_total_E,"total_E.png",plotDir);
    
    TH1D* h_nProngs_hist = (TH1D*)f_Root->Get("nProngs_hist");
    plot1D_Hist(h_nProngs_hist,"nProngs_hist.png",plotDir);
    
    TH1D* h_deltaInvMass_reco = (TH1D*)f_Root->Get("deltaInvMass_reco");
    plot1D_Hist(h_deltaInvMass_reco,"deltaInvMass_reco.png",plotDir);
    
    TH1D* h_vertex_z_reco= (TH1D*)f_Root->Get("vertex_z_reco");
    plot1D_Hist(h_vertex_z_reco,"vertex_z_reco.png",plotDir);
    
    TH2D* h_vertex_x_y_reco= (TH2D*)f_Root->Get("vertex_x_y_reco");
    plot2D_Hist(h_vertex_x_y_reco,"vertex_x_y_reco.png",plotDir);
    
    TH1D* h_pFilter_Status = (TH1D*)f_Root->Get("pFilter_Status");
    plot1D_Hist(h_pFilter_Status,"pFilter_Status.png",plotDir);
    
    TH1D* h_pFilter_RejectedEnergy = (TH1D*)f_Root->Get("pFilter_RejectedEnergy");
    plot1D_Hist(h_pFilter_RejectedEnergy,"pFilter_RejectedEnergy.png",plotDir);
    
    TH1D* h_beamEnergy_reco= (TH1D*)f_Root->Get("beamEnergy_reco");
    plot1D_Hist(h_beamEnergy_reco,"beamEnergy_reco.png",plotDir);
    
    TH1D* h_beamEnergyCal_reco= (TH1D*)f_Root->Get("beamEnergyCal_reco");
    plot1D_Hist(h_beamEnergyCal_reco,"beamEnergyCal_reco.png",plotDir);
    
    TH1D* h_q2_reco= (TH1D*)f_Root->Get("q2_reco");
    plot1D_Hist(h_q2_reco,"q2_reco.png",plotDir);    
    
    TH1D* h_w_reco= (TH1D*)f_Root->Get("w_reco");
    plot1D_Hist(h_w_reco,"w_reco.png",plotDir); 
    
    TH1D* h_wSq_reco= (TH1D*)f_Root->Get("wSq_reco");
    plot1D_Hist(h_wSq_reco,"wSq_reco.png",plotDir); 

    TH2D* h_total_E_neutrinoE = (TH2D*)f_Root->Get("total_E_neutrinoE");
    plot2D_Hist(h_total_E_neutrinoE,"total_E_neutrinoE.png",plotDir);
    
    TH2D* h_deltaInvMass_reco_mc = (TH2D*)f_Root->Get("deltaInvMass_reco_mc");
    plot2D_Hist(h_deltaInvMass_reco_mc,"deltaInvMass_reco_mc.png",plotDir);
    
    TH1D* h_deltaInvMass_error = (TH1D*)f_Root->Get("deltaInvMass_error");
    plot1D_Hist(h_deltaInvMass_error,"deltaInvMass_error.png",plotDir);
        
    TH2D* h_beamEnergy_reco_mc= (TH2D*)f_Root->Get("beamEnergy_reco_mc");
    plot2D_Hist(h_beamEnergy_reco_mc,"beamEnergy_reco_mc.png",plotDir);

    TH1D* h_beamEnergy_error= (TH1D*)f_Root->Get("beamEnergy_error");
    plot1D_Hist(h_beamEnergy_error,"beamEnergy_error.png",plotDir);
    
    TH2D* h_beamEnergyCal_reco_mc= (TH2D*)f_Root->Get("beamEnergyCal_reco_mc");
    plot2D_Hist(h_beamEnergyCal_reco_mc,"beamEnergyCal_reco_mc.png",plotDir);
    
    TH1D* h_beamEnergyCal_error= (TH1D*)f_Root->Get("beamEnergyCal_error");
    plot1D_Hist(h_beamEnergyCal_error,"beamEnergyCal_error.png",plotDir);
    
    TH2D* h_beamEnergy_beamEnergyCal= (TH2D*)f_Root->Get("beamEnergy_beamEnergyCal");
    plot2D_Hist(h_beamEnergy_beamEnergyCal,"beamEnergy_beamEnergyCal.png",plotDir);
    
    TH2D* h_q2_reco_mc= (TH2D*)f_Root->Get("q2_reco_mc");
    plot2D_Hist(h_q2_reco_mc,"q2_reco_mc.png",plotDir);
    
    TH1D* h_q2_error= (TH1D*)f_Root->Get("q2_error");
    plot1D_Hist(h_q2_error,"q2_error.png",plotDir);
    
    TH2D* h_w_reco_mc= (TH2D*)f_Root->Get("w_reco_mc");
    plot2D_Hist(h_w_reco_mc,"w_reco_mc.png",plotDir);
    
    TH1D* h_w_error= (TH1D*)f_Root->Get("w_error");
    plot1D_Hist(h_w_error,"w_error.png",plotDir);

    TH2D* h_vertex_z_reco_mc= (TH2D*)f_Root->Get("vertex_z_reco_mc");
    plot2D_Hist(h_vertex_z_reco_mc,"vertex_z_reco_mc.png",plotDir);
    
    TH1D* h_vertex_z_error= (TH1D*)f_Root->Get("vertex_z_error");
    plot1D_Hist(h_vertex_z_error,"vertex_z_error.png",plotDir);
    
    delete f_Root;
    
}


void CCProtonPi0_Plotter::plotParticleInfo(  string rootDir, string plotDir)
{

    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());

    TH1D* h_E_mc = (TH1D*)f_Root->Get("E_mc");
    plot1D_Hist(h_E_mc,"P_mc.png",plotDir);
    
    TH1D* h_P_mc = (TH1D*)f_Root->Get("P_mc");
    plot1D_Hist(h_P_mc,"P_mc.png",plotDir);
    
    TH1D* h_KE_mc = (TH1D*)f_Root->Get("KE_mc");
    plot1D_Hist(h_KE_mc,"KE_mc.png",plotDir);
    
    TH1D* h_angleBeam_mc = (TH1D*)f_Root->Get("angleBeam_mc");
    plot1D_Hist(h_angleBeam_mc,"angleBeam_mc.png",plotDir);
    
    TH1D* h_angleMuon_mc = (TH1D*)f_Root->Get("angleMuon_mc");
    plot1D_Hist(h_angleMuon_mc,"angleMuon_mc.png",plotDir);

    TH1D* h_partScore = (TH1D*)f_Root->Get("partScore");
    plot1D_Hist(h_partScore,"partScore.png",plotDir);
    
    TH1D* h_E_reco = (TH1D*)f_Root->Get("E_reco");
    plot1D_Hist(h_E_reco,"E_reco.png",plotDir);
    plot1D_HistLogScale(h_E_reco,"E_reco_log.png",plotDir);
    
    TH1D* h_P_reco = (TH1D*)f_Root->Get("P_reco");
    plot1D_Hist(h_P_reco,"P_reco.png",plotDir);
    plot1D_HistLogScale(h_P_reco,"P_reco_log.png",plotDir);
    
    TH1D* h_KE_reco = (TH1D*)f_Root->Get("KE_reco");
    plot1D_Hist(h_KE_reco,"KE_reco.png",plotDir);    
    
    TH1D* h_angleBeam_reco = (TH1D*)f_Root->Get("angleBeam_reco");
    plot1D_Hist(h_angleBeam_reco,"angleBeam_reco.png",plotDir);
    
    TH1D* h_angleMuon_reco = (TH1D*)f_Root->Get("angleMuon_reco");
    plot1D_Hist(h_angleMuon_reco,"angleMuon_reco.png",plotDir);

    TH2D* h_E_reco_mc = (TH2D*)f_Root->Get("E_reco_mc");
    plot2D_Hist(h_E_reco_mc,"E_reco_mc.png",plotDir);
    
    TH2D* h_P_reco_mc = (TH2D*)f_Root->Get("P_reco_mc");
    plot2D_Hist(h_P_reco_mc,"P_reco_mc.png",plotDir);
    
    TH1D* h_E_error = (TH1D*)f_Root->Get("E_error");
    plot1D_Hist(h_E_error,"E_error.png",plotDir);

    TH1D* h_P_error = (TH1D*)f_Root->Get("P_error");
    plot1D_Hist(h_P_error,"P_error.png",plotDir);
    
    TH2D* h_KE_reco_mc = (TH2D*)f_Root->Get("KE_reco_mc");
    plot2D_Hist(h_KE_reco_mc,"KE_reco_mc.png",plotDir);
    
    TH1D* h_KE_error = (TH1D*)f_Root->Get("KE_error");
    plot1D_Hist(h_KE_error,"KE_error.png",plotDir);
    
    TH2D* h_angleBeam_reco_mc = (TH2D*)f_Root->Get("angleBeam_reco_mc");
    plot2D_Hist(h_angleBeam_reco_mc,"angleBeam_reco_mc.png",plotDir);
    
    TH1D* h_angleBeam_error = (TH1D*)f_Root->Get("angleBeam_error");
    plot1D_Hist(h_angleBeam_error,"angleBeam_error.png",plotDir);
    
    TH2D* h_angleMuon_reco_mc = (TH2D*)f_Root->Get("angleMuon_reco_mc");
    plot2D_Hist(h_angleMuon_reco_mc,"angleMuon_reco_mc.png",plotDir);
    
    TH1D* h_angleMuon_error = (TH1D*)f_Root->Get("angleMuon_error");
    plot1D_Hist(h_angleMuon_error,"angleMuon_error.png",plotDir);

    delete f_Root;
}

void CCProtonPi0_Plotter::plotMuon()
{
    plotParticleInfo(rootDir_Muon[branchInd], plotDir_Muon[branchInd]);
}

void CCProtonPi0_Plotter::plotProton()
{    
    string rootDir = rootDir_Proton[branchInd];
    string plotDir = plotDir_Proton[branchInd];
    
    // Standard Plots
    plotParticleInfo(rootDir, plotDir);
    
    // Unique Plots
    TFile* f_Root = new TFile(rootDir.c_str());
    
    TH1D* h_trackLength = (TH1D*)f_Root->Get("trackLength");
    plot1D_Hist(h_trackLength,"trackLength.png",plotDir);
    
    TH1D* h_trackKinked = (TH1D*)f_Root->Get("trackKinked");
    plot1D_Hist(h_trackKinked,"trackKinked.png",plotDir);
    
}

void CCProtonPi0_Plotter::plotPion()
{
    string rootDir = rootDir_Pion[branchInd];
    string plotDir = plotDir_Pion[branchInd];
    
    // Standard Plots
    plotParticleInfo(rootDir, plotDir);
    
    // Unique Plots
    TFile* f_Root = new TFile(rootDir.c_str());
    
    TH1D* h_invMass = (TH1D*)f_Root->Get("invMass");
    plot1D_Hist(h_invMass,"invMass.png",plotDir);
    
    // Photon Conversion Length
    TH1D* h_gamma1_ConvLength = (TH1D*)f_Root->Get("gamma1_ConvLength");
    plot1D_Hist(h_gamma1_ConvLength,"gamma1_ConvLength.png",plotDir);
    
    TH1D* h_gamma2_ConvLength = (TH1D*)f_Root->Get("gamma2_ConvLength");
    plot1D_Hist(h_gamma2_ConvLength,"gamma2_ConvLength.png",plotDir);
    
    TH2D* h_ConvLength_gamma2_gamma1 = (TH2D*)f_Root->Get("ConvLength_gamma2_gamma1");
    plot2D_Hist(h_ConvLength_gamma2_gamma1,"ConvLength_gamma2_gamma1.png",plotDir);
    
    // nClusters
    TH2D* h_nClusters_All_gamma2_gamma1 = (TH2D*)f_Root->Get("nClusters_All_gamma2_gamma1");
    plot2D_Hist(h_nClusters_All_gamma2_gamma1,"nClusters_All_gamma2_gamma1.png",plotDir);
    
    // Gamma Energy
    TH2D* h_Energy_gamma2_gamma1 = (TH2D*)f_Root->Get("Energy_gamma2_gamma1");
    plot2D_Hist(h_Energy_gamma2_gamma1,"Energy_gamma2_gamma1.png",plotDir);
    
    TH1D* h_photonEnergy_Asymmetry = (TH1D*)f_Root->Get("photonEnergy_Asymmetry");
    plot1D_Hist(h_photonEnergy_Asymmetry,"photonEnergy_Asymmetry.png",plotDir);
    
    TH1D* h_photonEnergy_Asymmetry_true = (TH1D*)f_Root->Get("photonEnergy_Asymmetry_true");
    plot1D_Hist(h_photonEnergy_Asymmetry_true,"photonEnergy_Asymmetry_true.png",plotDir);
    
    delete f_Root;
}

void CCProtonPi0_Plotter::plotCutHistograms()
{
    TFile* f_Root_Signal_Interaction = new TFile(rootDir_Interaction[0].c_str());
    TFile* f_Root_Background_Interaction = new TFile(rootDir_Interaction[1].c_str());
    
    TH1D* h_signal_hCut_vertexCount= (TH1D*)f_Root_Signal_Interaction->Get("hCut_vertexCount");
    TH1D* h_background_hCut_vertexCount = (TH1D*)f_Root_Background_Interaction->Get("hCut_vertexCount");
    plotStacked(h_signal_hCut_vertexCount, h_background_hCut_vertexCount,"Number of Vertices", "hCut_vertexCount.png", otherDir);
    
    TH1D* h_signal_hCut_nProngs= (TH1D*)f_Root_Signal_Interaction->Get("hCut_nProngs");
    TH1D* h_background_hCut_nProngs = (TH1D*)f_Root_Background_Interaction->Get("hCut_nProngs");

    
    cout<<"Signal with nProngs 1 = "<<h_signal_hCut_nProngs->GetBinContent(2)<<endl;
    cout<<"Signal with nProngs 2 = "<<h_signal_hCut_nProngs->GetBinContent(3)<<endl;
    cout<<"Signal with nProgns 3+ = "<< h_signal_hCut_nProngs->GetBinContent(4)+
                                        h_signal_hCut_nProngs->GetBinContent(5)+
                                        h_signal_hCut_nProngs->GetBinContent(6)+
                                        h_signal_hCut_nProngs->GetBinContent(7)<<endl;
                                        
    cout<<"Background with nProngs 1 = "<<h_background_hCut_nProngs->GetBinContent(2)<<endl;
    cout<<"Background with nProngs 2 = "<<h_background_hCut_nProngs->GetBinContent(3)<<endl;
    cout<<"Background with nProgns 3+ = "<< h_background_hCut_nProngs->GetBinContent(4)+
                                        h_background_hCut_nProngs->GetBinContent(5)+
                                        h_background_hCut_nProngs->GetBinContent(6)+
                                        h_background_hCut_nProngs->GetBinContent(7)<<endl;
    
    
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
    
    TH1D* h_signal_hCut_protonScore_LLR = (TH1D*)f_Root_Signal_Interaction->Get("hCut_protonScore_LLR");
    TH1D* h_background_hCut_protonScore_LLR = (TH1D*)f_Root_Background_Interaction->Get("hCut_protonScore_LLR");
    plotStacked(h_signal_hCut_protonScore_LLR, h_background_hCut_protonScore_LLR,"Proton Score_LLR", "hCut_protonScore_LLR.png", otherDir);
    
    TH1D* h_signal_hCut_pIDDiff = (TH1D*)f_Root_Signal_Interaction->Get("hCut_pIDDiff");
    TH1D* h_background_hCut_pIDDiff = (TH1D*)f_Root_Background_Interaction->Get("hCut_pIDDiff");
    plotStacked(h_signal_hCut_pIDDiff, h_background_hCut_pIDDiff,"Proton Score - Pion Score", "hCut_pIDDiff.png", otherDir);
    
    TH1D* h_signal_hCut_deltaInvMass = (TH1D*)f_Root_Signal_Interaction->Get("hCut_deltaInvMass");
    TH1D* h_background_hCut_deltaInvMass = (TH1D*)f_Root_Background_Interaction->Get("hCut_deltaInvMass");
    plotStacked(h_signal_hCut_deltaInvMass, h_background_hCut_deltaInvMass,"Delta+ Invariant Mass", "hCut_deltaInvMass.png", otherDir);
    
}


void CCProtonPi0_Plotter::plotMichel()
{
    string rootDir = rootDir_Interaction[branchInd];
    
    inform(rootDir, otherDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());
    
    TH1D* h_N_michelElectrons = (TH1D*)f_Root->Get("N_michelElectrons");
    plot1D_Hist(h_N_michelElectrons,"N_michelElectrons.png",otherDir);
    
    TH1D* h_trueMichel_dist_vtx = (TH1D*)f_Root->Get("trueMichel_dist_vtx");
    TH1D* h_fakeMichel_dist_vtx = (TH1D*)f_Root->Get("fakeMichel_dist_vtx");
    plotStacked(h_fakeMichel_dist_vtx, h_trueMichel_dist_vtx,"Michel Distance to Vertex", "Michel_dist_vtx.png", otherDir,"NO Michel", "TRUE Michel");
//     plot1D_Hist(h_trueMichel_dist_vtx,"trueMichel_dist_vtx.png",otherDir); 
//     plot1D_Hist(h_fakeMichel_dist_vtx,"fakeMichel_dist_vtx.png",otherDir);
    
    TH1D* h_trueMichel_dist_end_point = (TH1D*)f_Root->Get("trueMichel_dist_end_point");
    TH1D* h_fakeMichel_dist_end_point = (TH1D*)f_Root->Get("fakeMichel_dist_end_point");
    plotStacked(h_fakeMichel_dist_end_point, h_trueMichel_dist_end_point,"Michel Distance to End Point", "Michel_dist_end_point.png", otherDir,"NO Michel", "TRUE Michel");
//     plot1D_Hist(h_trueMichel_dist_end_point,"trueMichel_dist_end_pointo.png",otherDir);  
//     plot1D_Hist(h_fakeMichel_dist_end_point,"fakeMichel_dist_end_point.png",otherDir);
    
    TH1D* h_trueMichel_end_Z = (TH1D*)f_Root->Get("trueMichel_end_Z");
    TH1D* h_fakeMichel_end_Z = (TH1D*)f_Root->Get("fakeMichel_end_Z");
    plotStacked(h_fakeMichel_end_Z, h_trueMichel_end_Z,"Michel Prong End Z [mm]", "Michel_end_Z.png", otherDir,"NO Michel", "TRUE Michel");
//     plot1D_Hist(h_trueMichel_end_Z,"trueMichel_end_Z.png",otherDir);  
//     plot1D_Hist(h_fakeMichel_end_Z,"fakeMichel_end_Z.png",otherDir);
    
    TH1D* h_trueMichel_energy = (TH1D*)f_Root->Get("trueMichel_energy");
    TH1D* h_fakeMichel_energy = (TH1D*)f_Root->Get("fakeMichel_energy");
    plotStacked(h_fakeMichel_energy, h_trueMichel_energy,"Michel Prong Energy", "Michel_energy.png", otherDir,"NO Michel", "TRUE Michel");
//     plot1D_Hist(h_trueMichel_energy,"trueMichel_energy.png",otherDir);  
//     plot1D_Hist(h_fakeMichel_energy,"fakeMichel_energy.png",otherDir);
     
    TH1D* h_trueMichel_time_diff = (TH1D*)f_Root->Get("trueMichel_time_diff");
    TH1D* h_fakeMichel_time_diff = (TH1D*)f_Root->Get("fakeMichel_time_diff");
    plotStacked(h_fakeMichel_time_diff, h_trueMichel_time_diff,"Michel Time Difference", "Michel_time_diff.png", otherDir,"NO Michel", "TRUE Michel");  
//     plot1D_Hist(h_trueMichel_time_diff,"trueMichel_time_diff.png",otherDir);  
//     plot1D_Hist(h_fakeMichel_time_diff,"fakeMichel_time_diff.png",otherDir);

    
    TH1D* h_trueMichel_end_Z_vtx_Z = (TH1D*)f_Root->Get("trueMichel_end_Z_vtx_Z");
    TH1D* h_fakeMichel_end_Z_vtx_Z = (TH1D*)f_Root->Get("fakeMichel_end_Z_vtx_Z");
    plotStacked(h_fakeMichel_end_Z_vtx_Z, h_trueMichel_end_Z_vtx_Z,"Michel Z - Vertex Z", "Michel_end_Z_vtx_Z.png", otherDir,"NO Michel", "TRUE Michel");
//     plot1D_Hist(h_trueMichel_end_Z_vtx_Z,"trueMichel_end_Z_vtx_Z.png",otherDir);  
//     plot1D_Hist(h_fakeMichel_end_Z_vtx_Z,"fakeMichel_end_Z_vtx_Z.png",otherDir);
    
 /*   
    TH1D* h_vertex_michelElectron_E = (TH1D*)f_Root->Get("vertex_michelElectron_E");
    TH1D* h_track_michelElectron_E = (TH1D*)f_Root->Get("track_michelElectron_E");
    TH1D* h_track2_michelElectron_E = (TH1D*)f_Root->Get("track2_michelElectron_E");
    TH1D* h_missed_michelElectron_E = (TH1D*)f_Root->Get("missed_michelElectron_E");
    TH1D* h_found_michelElectron_E = (TH1D*)f_Root->Get("found_michelElectron_E");
    MichelTool(h_vertex_michelElectron_E, h_track_michelElectron_E, h_track2_michelElectron_E, h_missed_michelElectron_E,
                "Michel Electron Energy [MeV]", "michelElectron_E.png", otherDir);
    plotStacked(h_found_michelElectron_E, h_missed_michelElectron_E, 
                            "Michel Electron Energy [MeV]", "michelElectron_E.png", otherDir, 
                            "Detected", "Missed");
                
    TH1D* h_vertex_michelMuon_P = (TH1D*)f_Root->Get("vertex_michelMuon_P");
    TH1D* h_track_michelMuon_P = (TH1D*)f_Root->Get("track_michelMuon_P");
    TH1D* h_track2_michelMuon_P = (TH1D*)f_Root->Get("track2_michelMuon_P");
    TH1D* h_missed_michelMuon_P = (TH1D*)f_Root->Get("missed_michelMuon_P");
    MichelTool(h_vertex_michelMuon_P, h_track_michelMuon_P, h_track2_michelMuon_P, h_missed_michelMuon_P,
                "Michel Muon Momentum [MeV]", "michelMuon_P.png", otherDir);

    TH1D* h_vertex_michelMuon_end_dist_vtx = (TH1D*)f_Root->Get("vertex_michelMuon_end_dist_vtx");
    TH1D* h_track_michelMuon_end_dist_vtx = (TH1D*)f_Root->Get("track_michelMuon_end_dist_vtx");
    TH1D* h_track2_michelMuon_end_dist_vtx = (TH1D*)f_Root->Get("track2_michelMuon_end_dist_vtx");
    TH1D* h_missed_michelMuon_end_dist_vtx = (TH1D*)f_Root->Get("missed_michelMuon_end_dist_vtx");
    MichelTool(h_vertex_michelMuon_end_dist_vtx, h_track_michelMuon_end_dist_vtx, h_track2_michelMuon_end_dist_vtx, h_missed_michelMuon_end_dist_vtx,
                "Michel Distance to Vertex [mm]", "michelMuon_end_dist_vtx.png", otherDir);
      
    TH1D* h_vertex_michelMuon_length = (TH1D*)f_Root->Get("vertex_michelMuon_length");
    TH1D* h_track_michelMuon_length = (TH1D*)f_Root->Get("track_michelMuon_length");
    TH1D* h_track2_michelMuon_length = (TH1D*)f_Root->Get("track2_michelMuon_length");
    TH1D* h_missed_michelMuon_length = (TH1D*)f_Root->Get("missed_michelMuon_length");
    MichelTool(h_vertex_michelMuon_length, h_track_michelMuon_length, h_track2_michelMuon_length, h_missed_michelMuon_length,
                "Michel Muon Length [mm]", "michelMuon_length.png", otherDir);
    
    TH1D* h_vertex_michelMuon_Z = (TH1D*)f_Root->Get("vertex_michelMuon_Z");
    TH1D* h_track_michelMuon_Z = (TH1D*)f_Root->Get("track_michelMuon_Z");
    TH1D* h_track2_michelMuon_Z = (TH1D*)f_Root->Get("track2_michelMuon_Z");
    TH1D* h_missed_michelMuon_Z = (TH1D*)f_Root->Get("missed_michelMuon_Z");
    MichelTool(h_vertex_michelMuon_Z, h_track_michelMuon_Z, h_track2_michelMuon_Z, h_missed_michelMuon_Z,
                "Michel Muon End Point Z [mm]", "michelMuon_Z.png", otherDir);
    
    TH1D* h_vertex_michelMuon_Z_vtx = (TH1D*)f_Root->Get("vertex_michelMuon_Z_vtx");
    TH1D* h_track_michelMuon_Z_vtx = (TH1D*)f_Root->Get("track_michelMuon_Z_vtx");
    TH1D* h_track2_michelMuon_Z_vtx = (TH1D*)f_Root->Get("track2_michelMuon_Z_vtx");
    TH1D* h_missed_michelMuon_Z_vtx = (TH1D*)f_Root->Get("missed_michelMuon_Z_vtx");
    MichelTool(h_vertex_michelMuon_Z_vtx, h_track_michelMuon_Z_vtx, h_track2_michelMuon_Z_vtx, h_missed_michelMuon_Z_vtx,
                "Michel Muon End Point Z - Vertex Z [mm]", "michelMuon_Z_vtx.png", otherDir);

    TH1D* h_vertex_michelPion_P = (TH1D*)f_Root->Get("vertex_michelPion_P");
    TH1D* h_track_michelPion_P = (TH1D*)f_Root->Get("track_michelPion_P");
    TH1D* h_track2_michelPion_P = (TH1D*)f_Root->Get("track2_michelPion_P");
    TH1D* h_missed_michelPion_P = (TH1D*)f_Root->Get("missed_michelPion_P");
    MichelTool(h_vertex_michelPion_P, h_track_michelPion_P, h_track2_michelPion_P, h_missed_michelPion_P,
                "Michel Pion Momentum [MeV]", "michelPion_P.png", otherDir);
    
    TH1D* h_vertex_michelPion_begin_dist_vtx = (TH1D*)f_Root->Get("vertex_michelPion_begin_dist_vtx");
    TH1D* h_track_michelPion_begin_dist_vtx = (TH1D*)f_Root->Get("track_michelPion_begin_dist_vtx");
    TH1D* h_track2_michelPion_begin_dist_vtx = (TH1D*)f_Root->Get("track2_michelPion_begin_dist_vtx");
    TH1D* h_missed_michelPion_begin_dist_vtx = (TH1D*)f_Root->Get("missed_michelPion_begin_dist_vtx");
    MichelTool(h_vertex_michelPion_begin_dist_vtx, h_track_michelPion_begin_dist_vtx, h_track2_michelPion_begin_dist_vtx, h_missed_michelPion_begin_dist_vtx,
                "Michel Pion Initial Point Distance to Vertex [mm]", "michelPion_begin_dist_vtx.png", otherDir);
                
    TH1D* h_vertex_michelPion_length = (TH1D*)f_Root->Get("vertex_michelPion_length");
    TH1D* h_track_michelPion_length = (TH1D*)f_Root->Get("track_michelPion_length");
    TH1D* h_track2_michelPion_length = (TH1D*)f_Root->Get("track2_michelPion_length");
    TH1D* h_missed_michelPion_length = (TH1D*)f_Root->Get("missed_michelPion_length");
    MichelTool(h_vertex_michelPion_length, h_track_michelPion_length, h_track2_michelPion_length, h_missed_michelPion_length,
                "Michel Pion Length [mm]", "michelPion_length.png", otherDir);*/
                
//     TH2D* h_vertex_michelMuon_dist_michelPion_length = (TH2D*)f_Root->Get("vertex_michelMuon_dist_michelPion_length");
//     plot2D_Hist(h_vertex_michelMuon_dist_michelPion_length,"vertex_michelMuon_dist_michelPion_length.png",plotDir);
//     
//     TH2D* h_missed_michelMuon_dist_michelPion_length = (TH2D*)f_Root->Get("missed_michelMuon_dist_michelPion_length");
//     plot2D_Hist(h_missed_michelMuon_dist_michelPion_length,"missed_michelMuon_dist_michelPion_length.png",plotDir);
//     
//     TH2D* h_vertex_michelMuon_X_Y = (TH2D*)f_Root->Get("vertex_michelMuon_X_Y");
//     plot2D_Hist(h_vertex_michelMuon_X_Y,"vertex_michelMuon_X_Y.png",plotDir);
//     
//     TH2D* h_missed_michelMuon_X_Y = (TH2D*)f_Root->Get("missed_michelMuon_X_Y");
//     plot2D_Hist(h_missed_michelMuon_X_Y,"missed_michelMuon_X_Y.png",plotDir);
}


void CCProtonPi0_Plotter::MichelTool(TH1D* h_vertex, TH1D* h_track, TH1D* h_track2, TH1D* h_missed,
                         string plotName, string fileName, string plotDir)
{    
    TCanvas* c1 = new TCanvas();
    THStack *hs = new THStack("hs",plotName.c_str());
    TLegend *legend = new TLegend(0.7,0.8,0.9,0.9);  
    
    h_vertex->SetFillColor(kGreen);
    h_vertex->SetMarkerStyle(21);
    h_vertex->SetMarkerColor(kGreen);
    
    h_track->SetFillColor(kBlue);
    h_track->SetMarkerStyle(21);
    h_track->SetMarkerColor(kBlue);
    
    h_track2->SetFillColor(kMagenta);
    h_track2->SetMarkerStyle(21);
    h_track2->SetMarkerColor(kMagenta);
    
    h_missed->SetFillColor(kRed);
    h_missed->SetMarkerStyle(21);
    h_missed->SetMarkerColor(kRed);
      
    legend->AddEntry(h_vertex, "Vertex", "f");
    legend->AddEntry(h_track, "Track", "f");
    legend->AddEntry(h_track2, "2nd Track", "f");
    legend->AddEntry(h_missed, "Missed", "f");
    
    hs->Add(h_missed);
    hs->Add(h_track2);
    hs->Add(h_track);
    hs->Add(h_vertex);
    hs->Draw();
    
    hs->GetXaxis()->SetTitle(plotName.c_str());
    hs->GetYaxis()->SetTitle("N(Events)");
    
    legend->Draw();
    
    c1->Print(Form("%s%s",plotDir.c_str(),fileName.c_str()),"png");
    
    delete hs;
    delete legend;
    delete c1; 
    
}

void CCProtonPi0_Plotter::plotSignalBackground()
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
//     TH1D* h_signal_protonScore = (TH1D*)f_Root_Signal_Proton->Get("protonScore");
//     TH1D* h_background_protonScore = (TH1D*)f_Root_Background_Proton->Get("protonScore");
//     plotStacked(h_signal_protonScore, h_background_protonScore,"Proton Score", "protonScore.png", otherDir);
//     
//     TH1D* h_signal_pionScore = (TH1D*)f_Root_Signal_Proton->Get("pionScore");
//     TH1D* h_background_pionScore = (TH1D*)f_Root_Background_Proton->Get("pionScore");
//     plotStacked(h_signal_pionScore, h_background_pionScore,"PionScore", "pionScore.png", otherDir);
        
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
//     TH1D* h_signal_gamma1_ConvLength = (TH1D*)f_Root_Signal_Pion->Get("gamma1_ConvLength");
//     TH1D* h_background_gamma1_ConvLength = (TH1D*)f_Root_Background_Pion->Get("gamma1_ConvLength");
//     plotStacked(h_signal_gamma1_ConvLength, h_background_gamma1_ConvLength,"Leading Photon Conversion Length [cm]", "gamma1_ConvLength.png", otherDir);
//     
//     TH1D* h_signal_gamma2_ConvLength = (TH1D*)f_Root_Signal_Pion->Get("gamma2_ConvLength");
//     TH1D* h_background_gamma2_ConvLength = (TH1D*)f_Root_Background_Pion->Get("gamma2_ConvLength");
//     plotStacked(h_signal_gamma2_ConvLength, h_background_gamma2_ConvLength,"Second Photon Conversion Length [cm]", "gamma2_ConvLength.png", otherDir);
    
    TH1D* h_signal_gamma1_nClusters_All = (TH1D*)f_Root_Signal_Pion->Get("gamma1_nClusters_All");
    TH1D* h_background_gamma1_nClusters_All = (TH1D*)f_Root_Background_Pion->Get("gamma1_nClusters_All");
    plotStacked(h_signal_gamma1_nClusters_All, h_background_gamma1_nClusters_All,"Leading Photon N(Clusters)", "gamma1_nClusters_All.png", otherDir);
    
    TH1D* h_signal_gamma2_nClusters_All = (TH1D*)f_Root_Signal_Pion->Get("gamma2_nClusters_All");
    TH1D* h_background_gamma2_nClusters_All = (TH1D*)f_Root_Background_Pion->Get("gamma2_nClusters_All");
    plotStacked(h_signal_gamma2_nClusters_All, h_background_gamma2_nClusters_All,"Second Photon N(Clusters)", "gamma2_nClusters_All.png", otherDir);
    
//     TH1D* h_signal_gamma1_nClusters_X = (TH1D*)f_Root_Signal_Pion->Get("gamma1_nClusters_X");
//     TH1D* h_background_gamma1_nClusters_X = (TH1D*)f_Root_Background_Pion->Get("gamma1_nClusters_X");
//     plotStacked(h_signal_gamma1_nClusters_X, h_background_gamma1_nClusters_X,"Number of X Clusters in the Leading Photon", "gamma1_nClusters_X.png", otherDir);
//     
//     TH1D* h_signal_gamma2_nClusters_X = (TH1D*)f_Root_Signal_Pion->Get("gamma2_nClusters_X");
//     TH1D* h_background_gamma2_nClusters_X = (TH1D*)f_Root_Background_Pion->Get("gamma2_nClusters_X");
//     plotStacked(h_signal_gamma2_nClusters_X, h_background_gamma2_nClusters_X,"Number of X Clusters in the Second Photon", "gamma2_nClusters_X.png", otherDir);
    
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
    
   
}

 
void CCProtonPi0_Plotter::plot_mc_w_Stacked()
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

void CCProtonPi0_Plotter::plot_final_mc_w_Stacked()
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

void CCProtonPi0_Plotter::plotPID()
{
    pID_proton();
    pID_proton_LLR();
    plot_2D_pID();
    pIDDiff();
    pIDStats();
    KE();
}

void CCProtonPi0_Plotter::KE()
{
    string rootDir = rootDir_PID[branchInd];
    string plotDir = plotDir_PID[branchInd];
    
    cout<<"\nPlottting Kinetic Energy Plots"<<endl;
    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());
    
    TH1D* h_KE_proton_pIDDiff = (TH1D*)f_Root->Get("KE_proton_pIDDiff");
    TH1D* h_KE_other_pIDDiff = (TH1D*)f_Root->Get("KE_other_pIDDiff");
    plotStacked(h_KE_proton_pIDDiff , h_KE_other_pIDDiff,"KE of True Protons(Green) and Other Particles(Red) using pIDDiff", "KE_proton_pIDDiff.png", plotDir);
    
    TH1D* h_KE_proton_LLR = (TH1D*)f_Root->Get("KE_proton_LLR");
    TH1D* h_KE_other_LLR = (TH1D*)f_Root->Get("KE_other_LLR");
    plotStacked(h_KE_proton_LLR , h_KE_other_LLR,"KE of True Protons(Green) and Other Particles(Red) using LLR", "KE_proton_LLR.png", plotDir);
}


void CCProtonPi0_Plotter::pIDStats()
{
    string rootDir = rootDir_PID[branchInd];
    string plotDir = plotDir_PID[branchInd];
    
    cout<<"\nPlottting pID Statistics"<<endl;
    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());
    
    TH1D* h_purity_LLR = (TH1D*)f_Root->Get("purity_LLR");
    plot1D_Hist(h_purity_LLR,"purity_LLR.png",plotDir);
    
    TH1D* h_efficiency_LLR = (TH1D*)f_Root->Get("efficiency_LLR");
    plot1D_Hist(h_efficiency_LLR,"efficiency_LLR.png",plotDir);
    
    TH1D* h_purityXefficiency_LLR = (TH1D*)f_Root->Get("purityXefficiency_LLR");
    plot1D_Hist(h_purityXefficiency_LLR,"purityXefficiency_LLR.png",plotDir);
    
    TH1D* h_purity_pIDDiff = (TH1D*)f_Root->Get("purity_pIDDiff");
    plot1D_Hist(h_purity_pIDDiff,"purity_pIDDiff.png",plotDir);
    
    TH1D* h_efficiency_pIDDiff = (TH1D*)f_Root->Get("efficiency_pIDDiff");
    plot1D_Hist(h_efficiency_pIDDiff,"efficiency_pIDDiff.png",plotDir);
    
    TH1D* h_purityXefficiency_pIDDiff = (TH1D*)f_Root->Get("purityXefficiency_pIDDiff");
    plot1D_Hist(h_purityXefficiency_pIDDiff,"purityXefficiency_pIDDiff.png",plotDir);
}

void CCProtonPi0_Plotter::pID_proton()
{
    string rootDir = rootDir_PID[branchInd];
    string plotDir = plotDir_PID[branchInd];
    
    cout<<"\nPlottting pID for Proton"<<endl;
    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());
    TCanvas* c1 = new TCanvas();
    THStack *hs = new THStack("hs","Proton Score");
    TLegend *legend = new TLegend(0.7,0.8,0.9,0.9); 
    
//     TH1D* h_pID_other = (TH1D*)f_Root->Get("other_protonScore");
//     h_pID_other->SetFillColor(kRed);
//     h_pID_other->SetMarkerStyle(21);
//     h_pID_other->SetMarkerColor(kRed);
    
    TH1D* h_pID_piminus = (TH1D*)f_Root->Get("piminus_protonScore");
    h_pID_piminus->SetFillColor(kYellow);
    h_pID_piminus->SetMarkerStyle(21);
    h_pID_piminus->SetMarkerColor(kYellow);
    
    TH1D* h_pID_piplus = (TH1D*)f_Root->Get("piplus_protonScore");
    h_pID_piplus->SetFillColor(kBlue);
    h_pID_piplus->SetMarkerStyle(21);
    h_pID_piplus->SetMarkerColor(kBlue);
    
    TH1D* h_pID_proton = (TH1D*)f_Root->Get("proton_protonScore");
    h_pID_proton->SetFillColor(kGreen);
    h_pID_proton->SetMarkerStyle(21);
    h_pID_proton->SetMarkerColor(kGreen);
    
    legend->AddEntry(h_pID_proton, "Proton", "f");
    legend->AddEntry(h_pID_piplus, "Pi-Plus", "f");
    legend->AddEntry(h_pID_piminus, "Pi-Minus", "f");
//     legend->AddEntry(h_pID_other, "Other", "f");
    
    
//     hs->Add(h_pID_other);
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

void CCProtonPi0_Plotter::pIDDiff()
{
    string rootDir = rootDir_PID[branchInd];
    string plotDir = plotDir_PID[branchInd];
    
    cout<<"\nPlottting pIDDiff"<<endl;
    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());
    TCanvas* c1 = new TCanvas();
    THStack *hs = new THStack("hs","Proton Score - Pion Score");
    TLegend *legend = new TLegend(0.7,0.8,0.9,0.9); 
    
    TH1D* h_pID_other = (TH1D*)f_Root->Get("other_pIDDiff");
    h_pID_other->SetFillColor(kRed);
    h_pID_other->SetMarkerStyle(21);
    h_pID_other->SetMarkerColor(kRed);
    
    TH1D* h_pID_piminus = (TH1D*)f_Root->Get("piminus_pIDDiff");
    h_pID_piminus->SetFillColor(kYellow);
    h_pID_piminus->SetMarkerStyle(21);
    h_pID_piminus->SetMarkerColor(kYellow);
    
    TH1D* h_pID_piplus = (TH1D*)f_Root->Get("piplus_pIDDiff");
    h_pID_piplus->SetFillColor(kBlue);
    h_pID_piplus->SetMarkerStyle(21);
    h_pID_piplus->SetMarkerColor(kBlue);
    
    TH1D* h_pID_proton = (TH1D*)f_Root->Get("proton_pIDDiff");
    h_pID_proton->SetFillColor(kGreen);
    h_pID_proton->SetMarkerStyle(21);
    h_pID_proton->SetMarkerColor(kGreen);
    
    legend->AddEntry(h_pID_proton, "Proton", "f");
    legend->AddEntry(h_pID_piplus, "Pi-Plus", "f");
    legend->AddEntry(h_pID_piminus, "Pi-Minus", "f");
    legend->AddEntry(h_pID_other, "Other", "f");
    
    hs->Add(h_pID_other);
    hs->Add(h_pID_piminus);
    hs->Add(h_pID_piplus);
    hs->Add(h_pID_proton);
    hs->Draw();
    hs->GetXaxis()->SetTitle("Proton Score - Pion Score");
    hs->GetYaxis()->SetTitle("N(Events)");
    
    legend->Draw();
    
    c1->Print(Form("%s%s",plotDir.c_str(),"stacked_pIDDiff.png"),"png");
    
    delete f_Root;
    delete hs;
    delete legend;
    delete c1;
    
}

void CCProtonPi0_Plotter::plot_2D_pID()
{
    string rootDir = rootDir_PID[branchInd];
    string plotDir = plotDir_PID[branchInd];
    
    cout<<"\nPlottting 2D pID Plots"<<endl;
    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());
    
    // Pion Score vs Proton Score
    TH2D* h_pID_proton_pionScore_protonScore = (TH2D*)f_Root->Get("proton_pionScore_protonScore");
    plot2D_Hist(h_pID_proton_pionScore_protonScore,"pID_proton_pionScore_protonScore.png",plotDir);
    
    TH2D* h_pID_piminus_pionScore_protonScore = (TH2D*)f_Root->Get("piminus_pionScore_protonScore");
    plot2D_Hist(h_pID_piminus_pionScore_protonScore,"pID_piminus_pionScore_protonScore.png",plotDir);
    
    TH2D* h_pID_piplus_pionScore_protonScore = (TH2D*)f_Root->Get("piplus_pionScore_protonScore");
    plot2D_Hist(h_pID_piplus_pionScore_protonScore,"pID_piplus_pionScore_protonScore.png",plotDir);
    
    // Proton Score vs Proton Score LLR
    TH2D* h_pID_proton_protonScore_protonScore_LLR = (TH2D*)f_Root->Get("proton_protonScore_protonScore_LLR");
    plot2D_Hist(h_pID_proton_protonScore_protonScore_LLR,"pID_proton_protonScore_protonScore_LLR.png",plotDir);
    
    TH2D* h_pID_piplus_protonScore_protonScore_LLR = (TH2D*)f_Root->Get("piplus_protonScore_protonScore_LLR");
    plot2D_Hist(h_pID_piplus_protonScore_protonScore_LLR,"pID_piplus_protonScore_protonScore_LLR.png",plotDir);
    
    TH2D* h_pID_piminus_protonScore_protonScore_LLR = (TH2D*)f_Root->Get("piminus_protonScore_protonScore_LLR");
    plot2D_Hist(h_pID_piminus_protonScore_protonScore_LLR,"pID_piminus_protonScore_protonScore_LLR.png",plotDir);
     
}


void CCProtonPi0_Plotter::pID_proton_LLR()
{
    string rootDir = rootDir_PID[branchInd];
    string plotDir = plotDir_PID[branchInd];
    
    cout<<"\nPlottting pID for Proton"<<endl;
    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());
    TCanvas* c1 = new TCanvas();
    THStack *hs = new THStack("hs","Proton Score Log-Likelihood Ratio");
    TLegend *legend = new TLegend(0.7,0.8,0.9,0.9); 
    
    TH1D* h_pID_other = (TH1D*)f_Root->Get("other_protonScore_LLR");
    h_pID_other->SetFillColor(kRed);
    h_pID_other->SetMarkerStyle(21);
    h_pID_other->SetMarkerColor(kRed);
    
    TH1D* h_pID_piminus = (TH1D*)f_Root->Get("piminus_protonScore_LLR");
    h_pID_piminus->SetFillColor(kYellow);
    h_pID_piminus->SetMarkerStyle(21);
    h_pID_piminus->SetMarkerColor(kYellow);
    
    TH1D* h_pID_piplus = (TH1D*)f_Root->Get("piplus_protonScore_LLR");
    h_pID_piplus->SetFillColor(kBlue);
    h_pID_piplus->SetMarkerStyle(21);
    h_pID_piplus->SetMarkerColor(kBlue);
    
    TH1D* h_pID_proton = (TH1D*)f_Root->Get("proton_protonScore_LLR");
    h_pID_proton->SetFillColor(kGreen);
    h_pID_proton->SetMarkerStyle(21);
    h_pID_proton->SetMarkerColor(kGreen);
    

    legend->AddEntry(h_pID_proton, "Proton", "f");
    legend->AddEntry(h_pID_piplus, "Pi-Plus", "f");
    legend->AddEntry(h_pID_piminus, "Pi-Minus", "f");
    legend->AddEntry(h_pID_other, "Other", "f");
    
    hs->Add(h_pID_other);
    hs->Add(h_pID_piminus);
    hs->Add(h_pID_piplus);
    hs->Add(h_pID_proton);
    hs->Draw();
    hs->GetXaxis()->SetTitle("Proton Score Log-Likelihood Ratio");
    hs->GetYaxis()->SetTitle("N(Events)");
    
    legend->Draw();
    
    c1->Print(Form("%s%s",plotDir.c_str(),"stacked_pID_protonScore_LLR.png"),"png");
    
    delete f_Root;
    delete hs;
    delete legend;
    delete c1;
    
}

#endif





