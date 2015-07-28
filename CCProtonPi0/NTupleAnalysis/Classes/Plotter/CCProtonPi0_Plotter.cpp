/*
    See CCProtonPi0_Plotter.h header for Class Information
*/

#ifndef CCProtonPi0_Plotter_cpp
#define CCProtonPi0_Plotter_cpp

#include "CCProtonPi0_Plotter.h"

using namespace PlotUtils;

void CCProtonPi0_Plotter::plotHistograms()
{
    //plotMuon();
    //plotProton();
    //plotPion();
    //plotInteraction();
    plotCutHistograms();
    
    if (isSignalvsBackground){
         //plotSignalBackground();
    }else{
         //plotInteraction();
               
//        plotPID();
        
//         plot_mc_w_Stacked();
//         plot_final_mc_w_Stacked();
//             plotMichel();
    }
   
}


CCProtonPi0_Plotter::CCProtonPi0_Plotter(int nMode, bool isMC)
{
    isSignalvsBackground = false;

    m_isMC = isMC;
    
    if ( nMode == 1) {
        std::cout<<"----------------------------------------------------------------------"<<std::endl;
        std::cout<<"Plot Mode: Signal - Only Signal Events will be Plotted"<<std::endl;
        branchInd = 0;
    }else if ( nMode == 2){
        std::cout<<"----------------------------------------------------------------------"<<std::endl;
        std::cout<<"Plot Mode: Background - Only Background Events will be Plotted"<<std::endl;
        branchInd = 1;
    }else if (nMode == 3){
        std::cout<<"----------------------------------------------------------------------"<<std::endl;
        std::cout<<"Plot Mode: All - All Events will be Plotted"<<std::endl;
        branchInd = 2;
    }else{
        std::cout<<"----------------------------------------------------------------------"<<std::endl;
        std::cout<<"Plot Mode: Signal vs Background"<<std::endl;
        isSignalvsBackground = true;
    }
    
    std::cout<<"----------------------------------------------------------------------"<<std::endl;
  
    rootDir_CutHists = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + "CutHistograms.root";
    setRootDirs(rootDir_Interaction,"Interaction.root");
    setRootDirs(rootDir_Muon,"Muon.root");
    setRootDirs(rootDir_Proton,"Proton.root");
    setRootDirs(rootDir_Pion,"Pion.root");
    setRootDirs(rootDir_PIDStatistics,"PIDStatistics.root");
    setRootDirs(rootDir_Pi0Blob,"Pi0Blob.root");
    setPlotDirs();
    
}

void CCProtonPi0_Plotter::setRootDirs(rootDir& dirs, std::string fileName )
{
    dirs.mc_signal = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + Folder_List::signal + fileName;
    dirs.mc_background = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + Folder_List::background + fileName;
    dirs.mc_all = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + Folder_List::allEvents + fileName;
    dirs.data = Folder_List::rootOut + Folder_List::Data + Folder_List::analyzed + Folder_List::allEvents + fileName;
}

void CCProtonPi0_Plotter::setPlotDirs()
{
    std::string dataType;

    // MC Data Locations
    for ( int i = 0 ; i < nBranches; i++){
        
        // Set Branch: 0 = Signal, 1 = Background, Other = All Events
        if ( i == 0 ) branchDir = Folder_List::signal;
        else if ( i == 1) branchDir = Folder_List::background;
        else branchDir = Folder_List::allEvents;
       
        if (m_isMC) dataType = Folder_List::MC;
        else dataType = Folder_List::Data;

        plotDir_Interaction[i] = Folder_List::output + Folder_List::plotOut + branchDir + "Interaction/";
        plotDir_PID[i] = Folder_List::output + Folder_List::plotOut + branchDir + "PIDStatistics/";
        plotDir_Muon[i] = Folder_List::output + Folder_List::plotOut + branchDir + "Muon/";
        plotDir_Proton[i] = Folder_List::output + Folder_List::plotOut + branchDir + "Proton/";
        plotDir_Pion[i] = Folder_List::output + Folder_List::plotOut + branchDir + "Pion/";
        plotDir_Pi0Blob[i] = Folder_List::output + Folder_List::plotOut + branchDir + "Pion/";
    }
    
    otherDir = Folder_List::output + Folder_List::plotOut + Folder_List::other;
}

void CCProtonPi0_Plotter::plotInteraction()
{
    std::string plotDir = plotDir_Interaction[branchInd];

    DrawDataMC(rootDir_Interaction,"Enu_1Track",plotDir);
    DrawDataMC(rootDir_Interaction,"Enu_2Track",plotDir);
    DrawDataMC(rootDir_Interaction,"Enu_Cal",plotDir);
    DrawDataMC(rootDir_Interaction,"q2",plotDir);
    DrawDataMC(rootDir_Interaction,"w",plotDir);
    DrawDataMC(rootDir_Interaction,"wSq",plotDir);
    DrawDataMC(rootDir_Interaction,"E_Unused_afterReco",plotDir);
    DrawDataMC(rootDir_Interaction,"E_Used_afterReco",plotDir);
    DrawDataMC(rootDir_Interaction,"deltaInvMass",plotDir);
    DrawDataMC(rootDir_Interaction,"nProngs_hist",plotDir);
}

void CCProtonPi0_Plotter::plotMuon()
{
    std::string plotDir = plotDir_Muon[branchInd];
    plotStandardHistograms(rootDir_Muon, plotDir);
}

void CCProtonPi0_Plotter::plotProton()
{    
    std::string plotDir = plotDir_Proton[branchInd];
    
    //plotStandardHistograms(rootDir_Proton, plotDir);
   
    // Unique Plots
    //DrawDataMC(rootDir_Proton,"trackLength",plotDir);
    //DrawDataMC(rootDir_Proton,"trackKinked",plotDir);
    //DrawDataMC(rootDir_Proton,"partScore",plotDir);

    Draw1DHist(rootDir_Proton,"P_error",plotDir);
    Draw2DHist(rootDir_Proton,"reco_P_true_P",plotDir);
}

void CCProtonPi0_Plotter::plotPion()
{
    std::string plotDir = plotDir_Pion[branchInd];
    
    // Standard Plots
    //plotStandardHistograms(rootDir_Pion, plotDir);
    
    // Unique Plots
    //DrawDataMC(rootDir_Pion,"gamma1_ConvLength",plotDir);
    //DrawDataMC(rootDir_Pion,"gamma2_ConvLength",plotDir);
    //DrawDataMC(rootDir_Pion,"photonEnergy_Asymmetry",plotDir);
    //DrawDataMC(rootDir_Pion,"invMass",plotDir);

    Draw1DHist(rootDir_Pion,"gamma1_P_error",plotDir);
    Draw2DHist(rootDir_Pion,"gamma1_reco_P_true_P",plotDir);
    Draw1DHist(rootDir_Pion,"gamma2_P_error",plotDir);
    Draw2DHist(rootDir_Pion,"gamma2_reco_P_true_P",plotDir);
    Draw2DHist(rootDir_Pion,"gamma1_P_gamma2_P",plotDir);
    Draw2DHist(rootDir_Pion,"gamma1_convLength_gamma2_convLength",plotDir);

}

void CCProtonPi0_Plotter::plotStandardHistograms(rootDir &dir, std::string plotDir)
{
    std::cout<<"Plotting Standard Histograms"<<std::endl;
    
    DrawDataMC(dir, "E", plotDir);
    DrawDataMC(dir, "P", plotDir);
    DrawDataMC(dir, "KE", plotDir);
    DrawDataMC(dir, "theta", plotDir);
    DrawDataMC(dir, "phi", plotDir);
}


void CCProtonPi0_Plotter::plotCutHistograms()
{
    std::string plotDir = plotDir_Interaction[branchInd];
    
    DrawStackedMC(rootDir_CutHists,"hCut_Michel",plotDir);
    DrawStackedMC(rootDir_CutHists,"hCut_eVis_nuclearTarget",plotDir);
    DrawStackedMC(rootDir_CutHists,"hCut_eVis_other",plotDir);
    DrawStackedMC(rootDir_CutHists,"hCut_pi0invMass",plotDir);
    DrawStackedMC(rootDir_CutHists,"hCut_gamma1ConvDist",plotDir);
    DrawStackedMC(rootDir_CutHists,"hCut_gamma2ConvDist",plotDir);
    
    DrawStackedMC(rootDir_CutHists,"hCut_1Prong_neutrinoE",plotDir);
    DrawStackedMC(rootDir_CutHists,"hCut_2Prong_neutrinoE",plotDir);
    DrawStackedMC(rootDir_CutHists,"hCut_1Prong_UnusedE",plotDir);
    DrawStackedMC(rootDir_CutHists,"hCut_2Prong_UnusedE",plotDir);
    
    DrawStackedMC(rootDir_CutHists,"hCut_protonScore_pIDDiff",plotDir);
    DrawStackedMC(rootDir_CutHists,"hCut_protonScore_LLR",plotDir);
    DrawStackedMC(rootDir_CutHists,"hCut_deltaInvMass",plotDir);
}

//
//void CCProtonPi0_Plotter::plotMichel()
//{
////    std::string rootDir = rootDir_Interaction[branchInd];
// //   
// //   inform(rootDir, otherDir);
// //   
// //   TFile* f_Root = new TFile(rootDir.c_str());
// //   
// //   TH1D* h_N_michelElectrons = (TH1D*)f_Root->Get("N_michelElectrons");
// //   Draw1DHist(h_N_michelElectrons,"N_michelElectrons.png",otherDir);
// //   
// //   TH1D* h_trueMichel_dist_vtx = (TH1D*)f_Root->Get("trueMichel_dist_vtx");
// //   TH1D* h_fakeMichel_dist_vtx = (TH1D*)f_Root->Get("fakeMichel_dist_vtx");
// //   plotStacked(h_fakeMichel_dist_vtx, h_trueMichel_dist_vtx,"Michel Distance to Vertex", "Michel_dist_vtx.png", otherDir,"NO Michel", "TRUE Michel");
//////     Draw1DHist(h_trueMichel_dist_vtx,"trueMichel_dist_vtx.png",otherDir); 
//////     Draw1DHist(h_fakeMichel_dist_vtx,"fakeMichel_dist_vtx.png",otherDir);
// //   
// //   TH1D* h_trueMichel_dist_end_point = (TH1D*)f_Root->Get("trueMichel_dist_end_point");
// //   TH1D* h_fakeMichel_dist_end_point = (TH1D*)f_Root->Get("fakeMichel_dist_end_point");
// //   plotStacked(h_fakeMichel_dist_end_point, h_trueMichel_dist_end_point,"Michel Distance to End Point", "Michel_dist_end_point.png", otherDir,"NO Michel", "TRUE Michel");
//////     Draw1DHist(h_trueMichel_dist_end_point,"trueMichel_dist_end_pointo.png",otherDir);  
//////     Draw1DHist(h_fakeMichel_dist_end_point,"fakeMichel_dist_end_point.png",otherDir);
// //   
// //   TH1D* h_trueMichel_end_Z = (TH1D*)f_Root->Get("trueMichel_end_Z");
// //   TH1D* h_fakeMichel_end_Z = (TH1D*)f_Root->Get("fakeMichel_end_Z");
// //   plotStacked(h_fakeMichel_end_Z, h_trueMichel_end_Z,"Michel Prong End Z [mm]", "Michel_end_Z.png", otherDir,"NO Michel", "TRUE Michel");
//////     Draw1DHist(h_trueMichel_end_Z,"trueMichel_end_Z.png",otherDir);  
//////     Draw1DHist(h_fakeMichel_end_Z,"fakeMichel_end_Z.png",otherDir);
// //   
// //   TH1D* h_trueMichel_energy = (TH1D*)f_Root->Get("trueMichel_energy");
// //   TH1D* h_fakeMichel_energy = (TH1D*)f_Root->Get("fakeMichel_energy");
// //   plotStacked(h_fakeMichel_energy, h_trueMichel_energy,"Michel Prong Energy", "Michel_energy.png", otherDir,"NO Michel", "TRUE Michel");
//////     Draw1DHist(h_trueMichel_energy,"trueMichel_energy.png",otherDir);  
//////     Draw1DHist(h_fakeMichel_energy,"fakeMichel_energy.png",otherDir);
// //    
// //   TH1D* h_trueMichel_time_diff = (TH1D*)f_Root->Get("trueMichel_time_diff");
// //   TH1D* h_fakeMichel_time_diff = (TH1D*)f_Root->Get("fakeMichel_time_diff");
// //   plotStacked(h_fakeMichel_time_diff, h_trueMichel_time_diff,"Michel Time Difference", "Michel_time_diff.png", otherDir,"NO Michel", "TRUE Michel");  
//////     Draw1DHist(h_trueMichel_time_diff,"trueMichel_time_diff.png",otherDir);  
//////     Draw1DHist(h_fakeMichel_time_diff,"fakeMichel_time_diff.png",otherDir);
//
// //   
// //   TH1D* h_trueMichel_end_Z_vtx_Z = (TH1D*)f_Root->Get("trueMichel_end_Z_vtx_Z");
// //   TH1D* h_fakeMichel_end_Z_vtx_Z = (TH1D*)f_Root->Get("fakeMichel_end_Z_vtx_Z");
// //   plotStacked(h_fakeMichel_end_Z_vtx_Z, h_trueMichel_end_Z_vtx_Z,"Michel Z - Vertex Z", "Michel_end_Z_vtx_Z.png", otherDir,"NO Michel", "TRUE Michel");
//////     Draw1DHist(h_trueMichel_end_Z_vtx_Z,"trueMichel_end_Z_vtx_Z.png",otherDir);  
//////     Draw1DHist(h_fakeMichel_end_Z_vtx_Z,"fakeMichel_end_Z_vtx_Z.png",otherDir);
// //   
// /*   
//    TH1D* h_vertex_michelElectron_E = (TH1D*)f_Root->Get("vertex_michelElectron_E");
//    TH1D* h_track_michelElectron_E = (TH1D*)f_Root->Get("track_michelElectron_E");
//    TH1D* h_track2_michelElectron_E = (TH1D*)f_Root->Get("track2_michelElectron_E");
//    TH1D* h_missed_michelElectron_E = (TH1D*)f_Root->Get("missed_michelElectron_E");
//    TH1D* h_found_michelElectron_E = (TH1D*)f_Root->Get("found_michelElectron_E");
//    MichelTool(h_vertex_michelElectron_E, h_track_michelElectron_E, h_track2_michelElectron_E, h_missed_michelElectron_E,
//                "Michel Electron Energy [MeV]", "michelElectron_E.png", otherDir);
//    plotStacked(h_found_michelElectron_E, h_missed_michelElectron_E, 
//                            "Michel Electron Energy [MeV]", "michelElectron_E.png", otherDir, 
//                            "Detected", "Missed");
//                
//    TH1D* h_vertex_michelMuon_P = (TH1D*)f_Root->Get("vertex_michelMuon_P");
//    TH1D* h_track_michelMuon_P = (TH1D*)f_Root->Get("track_michelMuon_P");
//    TH1D* h_track2_michelMuon_P = (TH1D*)f_Root->Get("track2_michelMuon_P");
//    TH1D* h_missed_michelMuon_P = (TH1D*)f_Root->Get("missed_michelMuon_P");
//    MichelTool(h_vertex_michelMuon_P, h_track_michelMuon_P, h_track2_michelMuon_P, h_missed_michelMuon_P,
//                "Michel Muon Momentum [MeV]", "michelMuon_P.png", otherDir);
//
//    TH1D* h_vertex_michelMuon_end_dist_vtx = (TH1D*)f_Root->Get("vertex_michelMuon_end_dist_vtx");
//    TH1D* h_track_michelMuon_end_dist_vtx = (TH1D*)f_Root->Get("track_michelMuon_end_dist_vtx");
//    TH1D* h_track2_michelMuon_end_dist_vtx = (TH1D*)f_Root->Get("track2_michelMuon_end_dist_vtx");
//    TH1D* h_missed_michelMuon_end_dist_vtx = (TH1D*)f_Root->Get("missed_michelMuon_end_dist_vtx");
//    MichelTool(h_vertex_michelMuon_end_dist_vtx, h_track_michelMuon_end_dist_vtx, h_track2_michelMuon_end_dist_vtx, h_missed_michelMuon_end_dist_vtx,
//                "Michel Distance to Vertex [mm]", "michelMuon_end_dist_vtx.png", otherDir);
//      
//    TH1D* h_vertex_michelMuon_length = (TH1D*)f_Root->Get("vertex_michelMuon_length");
//    TH1D* h_track_michelMuon_length = (TH1D*)f_Root->Get("track_michelMuon_length");
//    TH1D* h_track2_michelMuon_length = (TH1D*)f_Root->Get("track2_michelMuon_length");
//    TH1D* h_missed_michelMuon_length = (TH1D*)f_Root->Get("missed_michelMuon_length");
//    MichelTool(h_vertex_michelMuon_length, h_track_michelMuon_length, h_track2_michelMuon_length, h_missed_michelMuon_length,
//                "Michel Muon Length [mm]", "michelMuon_length.png", otherDir);
//    
//    TH1D* h_vertex_michelMuon_Z = (TH1D*)f_Root->Get("vertex_michelMuon_Z");
//    TH1D* h_track_michelMuon_Z = (TH1D*)f_Root->Get("track_michelMuon_Z");
//    TH1D* h_track2_michelMuon_Z = (TH1D*)f_Root->Get("track2_michelMuon_Z");
//    TH1D* h_missed_michelMuon_Z = (TH1D*)f_Root->Get("missed_michelMuon_Z");
//    MichelTool(h_vertex_michelMuon_Z, h_track_michelMuon_Z, h_track2_michelMuon_Z, h_missed_michelMuon_Z,
//                "Michel Muon End Point Z [mm]", "michelMuon_Z.png", otherDir);
//    
//    TH1D* h_vertex_michelMuon_Z_vtx = (TH1D*)f_Root->Get("vertex_michelMuon_Z_vtx");
//    TH1D* h_track_michelMuon_Z_vtx = (TH1D*)f_Root->Get("track_michelMuon_Z_vtx");
//    TH1D* h_track2_michelMuon_Z_vtx = (TH1D*)f_Root->Get("track2_michelMuon_Z_vtx");
//    TH1D* h_missed_michelMuon_Z_vtx = (TH1D*)f_Root->Get("missed_michelMuon_Z_vtx");
//    MichelTool(h_vertex_michelMuon_Z_vtx, h_track_michelMuon_Z_vtx, h_track2_michelMuon_Z_vtx, h_missed_michelMuon_Z_vtx,
//                "Michel Muon End Point Z - Vertex Z [mm]", "michelMuon_Z_vtx.png", otherDir);
//
//    TH1D* h_vertex_michelPion_P = (TH1D*)f_Root->Get("vertex_michelPion_P");
//    TH1D* h_track_michelPion_P = (TH1D*)f_Root->Get("track_michelPion_P");
//    TH1D* h_track2_michelPion_P = (TH1D*)f_Root->Get("track2_michelPion_P");
//    TH1D* h_missed_michelPion_P = (TH1D*)f_Root->Get("missed_michelPion_P");
//    MichelTool(h_vertex_michelPion_P, h_track_michelPion_P, h_track2_michelPion_P, h_missed_michelPion_P,
//                "Michel Pion Momentum [MeV]", "michelPion_P.png", otherDir);
//    
//    TH1D* h_vertex_michelPion_begin_dist_vtx = (TH1D*)f_Root->Get("vertex_michelPion_begin_dist_vtx");
//    TH1D* h_track_michelPion_begin_dist_vtx = (TH1D*)f_Root->Get("track_michelPion_begin_dist_vtx");
//    TH1D* h_track2_michelPion_begin_dist_vtx = (TH1D*)f_Root->Get("track2_michelPion_begin_dist_vtx");
//    TH1D* h_missed_michelPion_begin_dist_vtx = (TH1D*)f_Root->Get("missed_michelPion_begin_dist_vtx");
//    MichelTool(h_vertex_michelPion_begin_dist_vtx, h_track_michelPion_begin_dist_vtx, h_track2_michelPion_begin_dist_vtx, h_missed_michelPion_begin_dist_vtx,
//                "Michel Pion Initial Point Distance to Vertex [mm]", "michelPion_begin_dist_vtx.png", otherDir);
//                
//    TH1D* h_vertex_michelPion_length = (TH1D*)f_Root->Get("vertex_michelPion_length");
//    TH1D* h_track_michelPion_length = (TH1D*)f_Root->Get("track_michelPion_length");
//    TH1D* h_track2_michelPion_length = (TH1D*)f_Root->Get("track2_michelPion_length");
//    TH1D* h_missed_michelPion_length = (TH1D*)f_Root->Get("missed_michelPion_length");
//    MichelTool(h_vertex_michelPion_length, h_track_michelPion_length, h_track2_michelPion_length, h_missed_michelPion_length,
//                "Michel Pion Length [mm]", "michelPion_length.png", otherDir);*/
//                
////     TH2D* h_vertex_michelMuon_dist_michelPion_length = (TH2D*)f_Root->Get("vertex_michelMuon_dist_michelPion_length");
////     Draw2DHist(h_vertex_michelMuon_dist_michelPion_length,"vertex_michelMuon_dist_michelPion_length.png",plotDir);
////     
////     TH2D* h_missed_michelMuon_dist_michelPion_length = (TH2D*)f_Root->Get("missed_michelMuon_dist_michelPion_length");
////     Draw2DHist(h_missed_michelMuon_dist_michelPion_length,"missed_michelMuon_dist_michelPion_length.png",plotDir);
////     
////     TH2D* h_vertex_michelMuon_X_Y = (TH2D*)f_Root->Get("vertex_michelMuon_X_Y");
////     Draw2DHist(h_vertex_michelMuon_X_Y,"vertex_michelMuon_X_Y.png",plotDir);
////     
////     TH2D* h_missed_michelMuon_X_Y = (TH2D*)f_Root->Get("missed_michelMuon_X_Y");
////     Draw2DHist(h_missed_michelMuon_X_Y,"missed_michelMuon_X_Y.png",plotDir);
//}
//
//
//void CCProtonPi0_Plotter::MichelTool(TH1D* h_vertex, TH1D* h_track, TH1D* h_track2, TH1D* h_missed,
//                         std::string plotName, std::string fileName, std::string plotDir)
//{    
//    TCanvas* c1 = new TCanvas();
//    THStack *hs = new THStack("hs",plotName.c_str());
//    TLegend *legend = new TLegend(0.7,0.8,0.9,0.9);  
//    
//    h_vertex->SetFillColor(kGreen);
//    h_vertex->SetMarkerStyle(21);
//    h_vertex->SetMarkerColor(kGreen);
//    
//    h_track->SetFillColor(kBlue);
//    h_track->SetMarkerStyle(21);
//    h_track->SetMarkerColor(kBlue);
//    
//    h_track2->SetFillColor(kMagenta);
//    h_track2->SetMarkerStyle(21);
//    h_track2->SetMarkerColor(kMagenta);
//    
//    h_missed->SetFillColor(kRed);
//    h_missed->SetMarkerStyle(21);
//    h_missed->SetMarkerColor(kRed);
//      
//    legend->AddEntry(h_vertex, "Vertex", "f");
//    legend->AddEntry(h_track, "Track", "f");
//    legend->AddEntry(h_track2, "2nd Track", "f");
//    legend->AddEntry(h_missed, "Missed", "f");
//    
//    hs->Add(h_missed);
//    hs->Add(h_track2);
//    hs->Add(h_track);
//    hs->Add(h_vertex);
//    hs->Draw();
//    
//    hs->GetXaxis()->SetTitle(plotName.c_str());
//    hs->GetYaxis()->SetTitle("N(Events)");
//    
//    legend->Draw();
//    
//    c1->Print(Form("%s%s",plotDir.c_str(),fileName.c_str()),"png");
//    
//    delete hs;
//    delete legend;
//    delete c1; 
//    
//}
//
//void CCProtonPi0_Plotter::plotSignalBackground_Pi0Blob()
//{
//    TFile* f_Root_Signal_Pi0Blob = new TFile(rootDir_Pi0Blob[0].c_str());
//    TFile* f_Root_Background_Pi0Blob = new TFile(rootDir_Pi0Blob[1].c_str());
//    
//    TH1D* h_signal_g1_blob_energy = (TH1D*)f_Root_Signal_Pi0Blob->Get("g1_blob_energy");
//    TH1D* h_background_g1_blob_energy = (TH1D*)f_Root_Background_Pi0Blob->Get("g1_blob_energy");
//    plotStacked(h_signal_g1_blob_energy, h_background_g1_blob_energy,"g1_blob_energy", "g1_blob_energy.png", otherDir);
//    
//    TH1D* h_signal_g1_blob_nclusters = (TH1D*)f_Root_Signal_Pi0Blob->Get("g1_blob_nclusters");
//    TH1D* h_background_g1_blob_nclusters = (TH1D*)f_Root_Background_Pi0Blob->Get("g1_blob_nclusters");
//    plotStacked(h_signal_g1_blob_nclusters, h_background_g1_blob_nclusters,"g1_blob_nclusters", "g1_blob_nclusters.png", otherDir);
//
//    TH1D* h_signal_g1_blob_ndigits = (TH1D*)f_Root_Signal_Pi0Blob->Get("g1_blob_ndigits");
//    TH1D* h_background_g1_blob_ndigits = (TH1D*)f_Root_Background_Pi0Blob->Get("g1_blob_ndigits");
//    plotStacked(h_signal_g1_blob_ndigits, h_background_g1_blob_ndigits,"g1_blob_ndigits", "g1_blob_ndigits.png", otherDir);
//
//    TH1D* h_signal_g1_blob_ndof = (TH1D*)f_Root_Signal_Pi0Blob->Get("g1_blob_ndof");
//    TH1D* h_background_g1_blob_ndof = (TH1D*)f_Root_Background_Pi0Blob->Get("g1_blob_ndof");
//    plotStacked(h_signal_g1_blob_ndof, h_background_g1_blob_ndof,"g1_blob_ndof", "g1_blob_ndof.png", otherDir);
//     
//    TH1D* h_signal_g1_blob_fval = (TH1D*)f_Root_Signal_Pi0Blob->Get("g1_blob_fval");
//    TH1D* h_background_g1_blob_fval = (TH1D*)f_Root_Background_Pi0Blob->Get("g1_blob_fval");
//    plotStacked(h_signal_g1_blob_fval, h_background_g1_blob_fval,"g1_blob_fval", "g1_blob_fval.png", otherDir);
//    
//    TH1D* h_signal_g1_blob_dEdx_doublet = (TH1D*)f_Root_Signal_Pi0Blob->Get("g1_blob_dEdx_doublet");
//    TH1D* h_background_g1_blob_dEdx_doublet = (TH1D*)f_Root_Background_Pi0Blob->Get("g1_blob_dEdx_doublet");
//    plotStacked(h_signal_g1_blob_dEdx_doublet, h_background_g1_blob_dEdx_doublet,"g1_blob_dEdx_doublet", "g1_blob_dEdx_doublet.png", otherDir);
// 
//    TH1D* h_signal_g1_blob_dEdx_empty_plane = (TH1D*)f_Root_Signal_Pi0Blob->Get("g1_blob_dEdx_empty_plane");
//    TH1D* h_background_g1_blob_dEdx_empty_plane = (TH1D*)f_Root_Background_Pi0Blob->Get("g1_blob_dEdx_empty_plane");
//    plotStacked(h_signal_g1_blob_dEdx_empty_plane, h_background_g1_blob_dEdx_empty_plane,"g1_blob_dEdx_empty_plane", "g1_blob_dEdx_empty_plane.png", otherDir);
// 
//    TH1D* h_signal_g1_blob_dEdx = (TH1D*)f_Root_Signal_Pi0Blob->Get("g1_blob_dEdx");
//    TH1D* h_background_g1_blob_dEdx = (TH1D*)f_Root_Background_Pi0Blob->Get("g1_blob_dEdx");
//    plotStacked(h_signal_g1_blob_dEdx, h_background_g1_blob_dEdx,"g1_blob_dEdx", "g1_blob_dEdx.png", otherDir);
//
//    TH1D* h_signal_g1_blob_dEdx1 = (TH1D*)f_Root_Signal_Pi0Blob->Get("g1_blob_dEdx1");
//    TH1D* h_background_g1_blob_dEdx1 = (TH1D*)f_Root_Background_Pi0Blob->Get("g1_blob_dEdx1");
//    plotStacked(h_signal_g1_blob_dEdx1, h_background_g1_blob_dEdx1,"g1_blob_dEdx1", "g1_blob_dEdx1.png", otherDir);
//
//    TH1D* h_signal_g1_blob_dEdx_nplane = (TH1D*)f_Root_Signal_Pi0Blob->Get("g1_blob_dEdx_nplane");
//    TH1D* h_background_g1_blob_dEdx_nplane = (TH1D*)f_Root_Background_Pi0Blob->Get("g1_blob_dEdx_nplane");
//    plotStacked(h_signal_g1_blob_dEdx_nplane, h_background_g1_blob_dEdx_nplane,"g1_blob_dEdx_nplane", "g1_blob_dEdx_nplane.png", otherDir);
//
//    // Gamma 2
//    TH1D* h_signal_g2_blob_energy = (TH1D*)f_Root_Signal_Pi0Blob->Get("g2_blob_energy");
//    TH1D* h_background_g2_blob_energy = (TH1D*)f_Root_Background_Pi0Blob->Get("g2_blob_energy");
//    plotStacked(h_signal_g2_blob_energy, h_background_g2_blob_energy,"g2_blob_energy", "g2_blob_energy.png", otherDir);
//    
//    TH1D* h_signal_g2_blob_nclusters = (TH1D*)f_Root_Signal_Pi0Blob->Get("g2_blob_nclusters");
//    TH1D* h_background_g2_blob_nclusters = (TH1D*)f_Root_Background_Pi0Blob->Get("g2_blob_nclusters");
//    plotStacked(h_signal_g2_blob_nclusters, h_background_g2_blob_nclusters,"g2_blob_nclusters", "g2_blob_nclusters.png", otherDir);
//
//    TH1D* h_signal_g2_blob_ndigits = (TH1D*)f_Root_Signal_Pi0Blob->Get("g2_blob_ndigits");
//    TH1D* h_background_g2_blob_ndigits = (TH1D*)f_Root_Background_Pi0Blob->Get("g2_blob_ndigits");
//    plotStacked(h_signal_g2_blob_ndigits, h_background_g2_blob_ndigits,"g2_blob_ndigits", "g2_blob_ndigits.png", otherDir);
//
//}
//
//
//void CCProtonPi0_Plotter::plotSignalBackground()
//{
//    // Files
//    //TFile* f_Root_Signal_Interaction = new TFile(rootDir_Interaction[0].c_str());
//    //TFile* f_Root_Background_Interaction = new TFile(rootDir_Interaction[1].c_str());
////    TFile* f_Root_Signal_Muon = new TFile(rootDir_Muon[0].c_str());
////    TFile* f_Root_Background_Muon = new TFile(rootDir_Muon[1].c_str());
////    TFile* f_Root_Signal_Proton = new TFile(rootDir_Proton[0].c_str());
////    TFile* f_Root_Background_Proton = new TFile(rootDir_Proton[1].c_str());
////    TFile* f_Root_Signal_Pion = new TFile(rootDir_Pion[0].c_str());
////    TFile* f_Root_Background_Pion = new TFile(rootDir_Pion[1].c_str());
////
////    TH1D* h_signal_P_Muon = (TH1D*)f_Root_Signal_Muon->Get("P_reco");
////    TH1D* h_background_P_Muon = (TH1D*)f_Root_Background_Muon->Get("P_reco");
////    plotStacked(h_signal_P_Muon, h_background_P_Muon,"Muon Momentum [GeV]", "P_Muon.png", otherDir);
////    
////    TH1D* h_signal_angleBeam_reco_Muon = (TH1D*)f_Root_Signal_Muon->Get("angleBeam_reco");
////    TH1D* h_background_angleBeam_reco_Muon = (TH1D*)f_Root_Background_Muon->Get("angleBeam_reco");
////    plotStacked(h_signal_angleBeam_reco_Muon, h_background_angleBeam_reco_Muon,"Muon Angle wrt Beam", "angleBeam_reco_Muon.png", otherDir);
////
////    TH1D* h_signal_P_Proton = (TH1D*)f_Root_Signal_Proton->Get("P_reco");
////    TH1D* h_background_P_Proton = (TH1D*)f_Root_Background_Proton->Get("P_reco");
////    plotStacked(h_signal_P_Proton, h_background_P_Proton,"Proton Momentum [MeV]", "P_Proton.png", otherDir);
////    plotStackedLogScale(h_signal_P_Proton, h_background_P_Proton,"Proton Momentum [MeV]", "P_Proton_Log.png", otherDir);
////    
////    TH1D* h_signal_angleBeam_reco_Proton = (TH1D*)f_Root_Signal_Proton->Get("angleBeam_reco");
////    TH1D* h_background_angleBeam_reco_Proton = (TH1D*)f_Root_Background_Proton->Get("angleBeam_reco");
////    plotStacked(h_signal_angleBeam_reco_Proton, h_background_angleBeam_reco_Proton,"Proton Angle wrt Beam", "angleBeam_reco_Proton.png", otherDir);
////    plotStackedLogScale(h_signal_angleBeam_reco_Proton, h_background_angleBeam_reco_Proton,"Proton Angle wrt Beam", "angleBeam_reco_Proton_Log.png", otherDir);
////
////    TH1D* h_signal_P_Pion = (TH1D*)f_Root_Signal_Pion->Get("P_reco");
////    TH1D* h_background_P_Pion = (TH1D*)f_Root_Background_Pion->Get("P_reco");
////    plotStacked(h_signal_P_Pion, h_background_P_Pion,"Pion Momentum [MeV]", "P_Pion.png", otherDir);
////    
////    TH1D* h_signal_angleBeam_reco_Pion = (TH1D*)f_Root_Signal_Pion->Get("angleBeam_reco");
////    TH1D* h_background_angleBeam_reco_Pion = (TH1D*)f_Root_Background_Pion->Get("angleBeam_reco");
////    plotStacked(h_signal_angleBeam_reco_Pion, h_background_angleBeam_reco_Pion,"Pion Angle wrt Beam", "angleBeam_reco_Pion.png", otherDir);
////
//}
//
// 
//void CCProtonPi0_Plotter::plot_mc_w_Stacked()
//{
//   // std::string rootDir = rootDir_Interaction[branchInd];
//   // std::string plotDir = plotDir_Interaction[branchInd];
//   // 
//   // std::cout<<"\nPlottting Stacked mc_w"<<std::endl;
//   // inform(rootDir, plotDir);
//   // 
//   // TFile* f_Root = new TFile(rootDir.c_str());
//   // TCanvas* c1 = new TCanvas();
//   // THStack *hs = new THStack("hs","TRUE Signal Events (MINOS Matched)");
//   // TLegend *legend = new TLegend(0.7,0.8,0.9,0.9);  
//   // 
//   // TH1D* h_mc_w_DIS = (TH1D*)f_Root->Get("mc_w_DIS");
//   // h_mc_w_DIS->SetFillColor(kRed);
//   // h_mc_w_DIS->SetMarkerStyle(21);
//   // h_mc_w_DIS->SetMarkerColor(kRed);
//   // 
//   // TH1D* h_mc_w_RES = (TH1D*)f_Root->Get("mc_w_RES");
//   // h_mc_w_RES->SetFillColor(kBlue);
//   // h_mc_w_RES->SetMarkerStyle(21);
//   // h_mc_w_RES->SetMarkerColor(kBlue);
//   //   
//   // legend->AddEntry(h_mc_w_DIS , "DIS", "f");
//   // legend->AddEntry(h_mc_w_RES, "RES", "f");
//
//   // 
//   // hs->Add(h_mc_w_DIS);
//   // hs->Add(h_mc_w_RES);
//   // hs->Draw();
//   // hs->GetXaxis()->SetTitle("mc_w [GeV/c^{2}]");
//   // hs->GetYaxis()->SetTitle("N(Events)");
//   // 
//   // legend->Draw();
//   // 
//   // c1->Print(Form("%s%s",plotDir.c_str(),"mc_w.png"),"png");
//   // 
//   // delete f_Root;
//   // delete hs;
//   // delete legend;
//   // delete c1;
//}
//
//void CCProtonPi0_Plotter::plot_final_mc_w_Stacked()
//{
//   // std::string rootDir = rootDir_Interaction[branchInd];
//   // std::string plotDir = plotDir_Interaction[branchInd];
//   // 
//   // std::cout<<"\nPlottting Stacked final_mc_w"<<std::endl;
//   // inform(rootDir, plotDir);
//   // 
//   // TFile* f_Root = new TFile(rootDir.c_str());
//   // TCanvas* c1 = new TCanvas();
//   // THStack *hs = new THStack("hs","TRUE Signal Events (MINOS Matched)");
//   // TLegend *legend = new TLegend(0.7,0.8,0.9,0.9);  
//   // 
//   // TH1D* h_mc_w_DIS = (TH1D*)f_Root->Get("final_mc_w_DIS");
//   // h_mc_w_DIS->SetFillColor(kRed);
//   // h_mc_w_DIS->SetMarkerStyle(21);
//   // h_mc_w_DIS->SetMarkerColor(kRed);
//   // 
//   // TH1D* h_mc_w_RES = (TH1D*)f_Root->Get("final_mc_w_RES");
//   // h_mc_w_RES->SetFillColor(kBlue);
//   // h_mc_w_RES->SetMarkerStyle(21);
//   // h_mc_w_RES->SetMarkerColor(kBlue);
//   // 
//   // legend->AddEntry(h_mc_w_DIS , "DIS", "f");
//   // legend->AddEntry(h_mc_w_RES, "RES", "f");
//   // 
//   // 
//   // hs->Add(h_mc_w_DIS);
//   // hs->Add(h_mc_w_RES);
//   // hs->Draw();
//   // hs->GetXaxis()->SetTitle("mc_w [GeV/c^{2}]");
//   // hs->GetYaxis()->SetTitle("N(Events)");
//   // 
//   // legend->Draw();
//   // 
//   // c1->Print(Form("%s%s",plotDir.c_str(),"final_mc_w.png"),"png");
//   // 
//   // delete f_Root;
//   // delete hs;
//   // delete legend;
//   // delete c1;
//}
//
//void CCProtonPi0_Plotter::plotPID()
//{
//    pID_proton();
//    pID_proton_LLR();
//    plot_2D_pID();
//    pIDDiff();
//    pIDStats();
//    KE();
//}
//
//void CCProtonPi0_Plotter::KE()
//{
//    std::string rootDir = rootDir_PID[branchInd];
//    std::string plotDir = plotDir_PID[branchInd];
//    
//    std::cout<<"\nPlottting Kinetic Energy Plots"<<std::endl;
//    inform(rootDir, plotDir);
//    
//    TFile* f_Root = new TFile(rootDir.c_str());
//    
//    TH1D* h_KE_proton_pIDDiff = (TH1D*)f_Root->Get("KE_proton_pIDDiff");
//    TH1D* h_KE_other_pIDDiff = (TH1D*)f_Root->Get("KE_other_pIDDiff");
//    plotStacked(h_KE_proton_pIDDiff , h_KE_other_pIDDiff,"KE of True Protons(Green) and Other Particles(Red) using pIDDiff", "KE_proton_pIDDiff.png", plotDir);
//    
//    TH1D* h_KE_proton_LLR = (TH1D*)f_Root->Get("KE_proton_LLR");
//    TH1D* h_KE_other_LLR = (TH1D*)f_Root->Get("KE_other_LLR");
//    plotStacked(h_KE_proton_LLR , h_KE_other_LLR,"KE of True Protons(Green) and Other Particles(Red) using LLR", "KE_proton_LLR.png", plotDir);
//}
//
//
//void CCProtonPi0_Plotter::pIDStats()
//{
//    //std::string rootDir = rootDir_PID[branchInd];
//    //std::string plotDir = plotDir_PID[branchInd];
//    //
//    //std::cout<<"\nPlottting pID Statistics"<<std::endl;
//    //inform(rootDir, plotDir);
//    //
//    //TFile* f_Root = new TFile(rootDir.c_str());
//    //
//    //TH1D* h_purity_LLR = (TH1D*)f_Root->Get("purity_LLR");
//    //Draw1DHist(h_purity_LLR,"purity_LLR.png",plotDir);
//    //
//    //TH1D* h_efficiency_LLR = (TH1D*)f_Root->Get("efficiency_LLR");
//    //Draw1DHist(h_efficiency_LLR,"efficiency_LLR.png",plotDir);
//    //
//    //TH1D* h_purityXefficiency_LLR = (TH1D*)f_Root->Get("purityXefficiency_LLR");
//    //Draw1DHist(h_purityXefficiency_LLR,"purityXefficiency_LLR.png",plotDir);
//    //
//    //TH1D* h_purity_pIDDiff = (TH1D*)f_Root->Get("purity_pIDDiff");
//    //Draw1DHist(h_purity_pIDDiff,"purity_pIDDiff.png",plotDir);
//    //
//    //TH1D* h_efficiency_pIDDiff = (TH1D*)f_Root->Get("efficiency_pIDDiff");
//    //Draw1DHist(h_efficiency_pIDDiff,"efficiency_pIDDiff.png",plotDir);
//    //
//    //TH1D* h_purityXefficiency_pIDDiff = (TH1D*)f_Root->Get("purityXefficiency_pIDDiff");
//    //Draw1DHist(h_purityXefficiency_pIDDiff,"purityXefficiency_pIDDiff.png",plotDir);
//}
//
//void CCProtonPi0_Plotter::pID_proton()
//{
//    std::string rootDir = rootDir_PID[branchInd];
//    std::string plotDir = plotDir_PID[branchInd];
//    
//    std::cout<<"\nPlottting pID for Proton"<<std::endl;
//    inform(rootDir, plotDir);
//    
//    TFile* f_Root = new TFile(rootDir.c_str());
//    TCanvas* c1 = new TCanvas();
//    THStack *hs = new THStack("hs","Proton Score");
//    TLegend *legend = new TLegend(0.7,0.8,0.9,0.9); 
//    
////     TH1D* h_pID_other = (TH1D*)f_Root->Get("other_protonScore");
////     h_pID_other->SetFillColor(kRed);
////     h_pID_other->SetMarkerStyle(21);
////     h_pID_other->SetMarkerColor(kRed);
//    
//    TH1D* h_pID_piminus = (TH1D*)f_Root->Get("piminus_protonScore");
//    h_pID_piminus->SetFillColor(kYellow);
//    h_pID_piminus->SetMarkerStyle(21);
//    h_pID_piminus->SetMarkerColor(kYellow);
//    
//    TH1D* h_pID_piplus = (TH1D*)f_Root->Get("piplus_protonScore");
//    h_pID_piplus->SetFillColor(kBlue);
//    h_pID_piplus->SetMarkerStyle(21);
//    h_pID_piplus->SetMarkerColor(kBlue);
//    
//    TH1D* h_pID_proton = (TH1D*)f_Root->Get("proton_protonScore");
//    h_pID_proton->SetFillColor(kGreen);
//    h_pID_proton->SetMarkerStyle(21);
//    h_pID_proton->SetMarkerColor(kGreen);
//    
//    legend->AddEntry(h_pID_proton, "Proton", "f");
//    legend->AddEntry(h_pID_piplus, "Pi-Plus", "f");
//    legend->AddEntry(h_pID_piminus, "Pi-Minus", "f");
////     legend->AddEntry(h_pID_other, "Other", "f");
//    
//    
////     hs->Add(h_pID_other);
//    hs->Add(h_pID_piminus);
//    hs->Add(h_pID_piplus);
//    hs->Add(h_pID_proton);
//    hs->Draw();
//    hs->GetXaxis()->SetTitle("Proton Score");
//    hs->GetYaxis()->SetTitle("N(Events)");
//    
//    legend->Draw();
//    
//    c1->Print(Form("%s%s",plotDir.c_str(),"stacked_pID_protonScore.png"),"png");
//    
//    delete f_Root;
//    delete hs;
//    delete legend;
//    delete c1;
//    
//}
//
//void CCProtonPi0_Plotter::pIDDiff()
//{
//    std::string rootDir = rootDir_PID[branchInd];
//    std::string plotDir = plotDir_PID[branchInd];
//    
//    std::cout<<"\nPlottting pIDDiff"<<std::endl;
//    inform(rootDir, plotDir);
//    
//    TFile* f_Root = new TFile(rootDir.c_str());
//    TCanvas* c1 = new TCanvas();
//    THStack *hs = new THStack("hs","Proton Score - Pion Score");
//    TLegend *legend = new TLegend(0.7,0.8,0.9,0.9); 
//    
//    TH1D* h_pID_other = (TH1D*)f_Root->Get("other_pIDDiff");
//    h_pID_other->SetFillColor(kRed);
//    h_pID_other->SetMarkerStyle(21);
//    h_pID_other->SetMarkerColor(kRed);
//    
//    TH1D* h_pID_piminus = (TH1D*)f_Root->Get("piminus_pIDDiff");
//    h_pID_piminus->SetFillColor(kYellow);
//    h_pID_piminus->SetMarkerStyle(21);
//    h_pID_piminus->SetMarkerColor(kYellow);
//    
//    TH1D* h_pID_piplus = (TH1D*)f_Root->Get("piplus_pIDDiff");
//    h_pID_piplus->SetFillColor(kBlue);
//    h_pID_piplus->SetMarkerStyle(21);
//    h_pID_piplus->SetMarkerColor(kBlue);
//    
//    TH1D* h_pID_proton = (TH1D*)f_Root->Get("proton_pIDDiff");
//    h_pID_proton->SetFillColor(kGreen);
//    h_pID_proton->SetMarkerStyle(21);
//    h_pID_proton->SetMarkerColor(kGreen);
//    
//    legend->AddEntry(h_pID_proton, "Proton", "f");
//    legend->AddEntry(h_pID_piplus, "Pi-Plus", "f");
//    legend->AddEntry(h_pID_piminus, "Pi-Minus", "f");
//    legend->AddEntry(h_pID_other, "Other", "f");
//    
//    hs->Add(h_pID_other);
//    hs->Add(h_pID_piminus);
//    hs->Add(h_pID_piplus);
//    hs->Add(h_pID_proton);
//    hs->Draw();
//    hs->GetXaxis()->SetTitle("Proton Score - Pion Score");
//    hs->GetYaxis()->SetTitle("N(Events)");
//    
//    legend->Draw();
//    
//    c1->Print(Form("%s%s",plotDir.c_str(),"stacked_pIDDiff.png"),"png");
//    
//    delete f_Root;
//    delete hs;
//    delete legend;
//    delete c1;
//    
//}
//
//void CCProtonPi0_Plotter::plot_2D_pID()
//{
//    //std::string rootDir = rootDir_PID[branchInd];
//    //std::string plotDir = plotDir_PID[branchInd];
//    //
//    //std::cout<<"\nPlottting 2D pID Plots"<<std::endl;
//    //inform(rootDir, plotDir);
//    //
//    //TFile* f_Root = new TFile(rootDir.c_str());
//    //
//    //// Pion Score vs Proton Score
//    //TH2D* h_pID_proton_pionScore_protonScore = (TH2D*)f_Root->Get("proton_pionScore_protonScore");
//    //Draw2DHist(h_pID_proton_pionScore_protonScore,"pID_proton_pionScore_protonScore.png",plotDir);
//    //
//    //TH2D* h_pID_piminus_pionScore_protonScore = (TH2D*)f_Root->Get("piminus_pionScore_protonScore");
//    //Draw2DHist(h_pID_piminus_pionScore_protonScore,"pID_piminus_pionScore_protonScore.png",plotDir);
//    //
//    //TH2D* h_pID_piplus_pionScore_protonScore = (TH2D*)f_Root->Get("piplus_pionScore_protonScore");
//    //Draw2DHist(h_pID_piplus_pionScore_protonScore,"pID_piplus_pionScore_protonScore.png",plotDir);
//    //
//    //// Proton Score vs Proton Score LLR
//    //TH2D* h_pID_proton_protonScore_protonScore_LLR = (TH2D*)f_Root->Get("proton_protonScore_protonScore_LLR");
//    //Draw2DHist(h_pID_proton_protonScore_protonScore_LLR,"pID_proton_protonScore_protonScore_LLR.png",plotDir);
//    //
//    //TH2D* h_pID_piplus_protonScore_protonScore_LLR = (TH2D*)f_Root->Get("piplus_protonScore_protonScore_LLR");
//    //Draw2DHist(h_pID_piplus_protonScore_protonScore_LLR,"pID_piplus_protonScore_protonScore_LLR.png",plotDir);
//    //
//    //TH2D* h_pID_piminus_protonScore_protonScore_LLR = (TH2D*)f_Root->Get("piminus_protonScore_protonScore_LLR");
//    //Draw2DHist(h_pID_piminus_protonScore_protonScore_LLR,"pID_piminus_protonScore_protonScore_LLR.png",plotDir);
//     
//}
//
//
//void CCProtonPi0_Plotter::pID_proton_LLR()
//{
//    std::string rootDir = rootDir_PID[branchInd];
//    std::string plotDir = plotDir_PID[branchInd];
//    
//    std::cout<<"\nPlottting pID for Proton"<<std::endl;
//    inform(rootDir, plotDir);
//    
//    TFile* f_Root = new TFile(rootDir.c_str());
//    TCanvas* c1 = new TCanvas();
//    THStack *hs = new THStack("hs","Proton Score Log-Likelihood Ratio");
//    TLegend *legend = new TLegend(0.7,0.8,0.9,0.9); 
//    
//    TH1D* h_pID_other = (TH1D*)f_Root->Get("other_protonScore_LLR");
//    h_pID_other->SetFillColor(kRed);
//    h_pID_other->SetMarkerStyle(21);
//    h_pID_other->SetMarkerColor(kRed);
//    
//    TH1D* h_pID_piminus = (TH1D*)f_Root->Get("piminus_protonScore_LLR");
//    h_pID_piminus->SetFillColor(kYellow);
//    h_pID_piminus->SetMarkerStyle(21);
//    h_pID_piminus->SetMarkerColor(kYellow);
//    
//    TH1D* h_pID_piplus = (TH1D*)f_Root->Get("piplus_protonScore_LLR");
//    h_pID_piplus->SetFillColor(kBlue);
//    h_pID_piplus->SetMarkerStyle(21);
//    h_pID_piplus->SetMarkerColor(kBlue);
//    
//    TH1D* h_pID_proton = (TH1D*)f_Root->Get("proton_protonScore_LLR");
//    h_pID_proton->SetFillColor(kGreen);
//    h_pID_proton->SetMarkerStyle(21);
//    h_pID_proton->SetMarkerColor(kGreen);
//    
//
//    legend->AddEntry(h_pID_proton, "Proton", "f");
//    legend->AddEntry(h_pID_piplus, "Pi-Plus", "f");
//    legend->AddEntry(h_pID_piminus, "Pi-Minus", "f");
//    legend->AddEntry(h_pID_other, "Other", "f");
//    
//    hs->Add(h_pID_other);
//    hs->Add(h_pID_piminus);
//    hs->Add(h_pID_piplus);
//    hs->Add(h_pID_proton);
//    hs->Draw();
//    hs->GetXaxis()->SetTitle("Proton Score Log-Likelihood Ratio");
//    hs->GetYaxis()->SetTitle("N(Events)");
//    
//    legend->Draw();
//    
//    c1->Print(Form("%s%s",plotDir.c_str(),"stacked_pID_protonScore_LLR.png"),"png");
//    
//    delete f_Root;
//    delete hs;
//    delete legend;
//    delete c1;
//    
//}
//

//void CCProtonPi0_Plotter::plot_purity_efficiency(TH1D* h_signal, TH1D* h_background, std::string fileName, std::string plotDir, bool keepEventstoRight)
//{
    //// ------------------------------------------------------------------------
    //// Create Histograms
    //// ------------------------------------------------------------------------
    //TH1D* h_purity = new TH1D;
    //TH1D* h_efficiency = new TH1D;
    //TH1D* h_purity_efficiency = new TH1D;

    //// Copy and Reset Histogram Content - Keeps binning and Titles
    //h_signal->Copy(*h_purity);
    //h_signal->Copy(*h_efficiency);
    //h_signal->Copy(*h_purity_efficiency);
    //
    //h_purity->Reset();
    //h_efficiency->Reset();
    //h_purity_efficiency->Reset();
    //
    //// Change Titles for the new histograms
    //if(keepEventstoRight) h_purity->SetTitle("Purity - (Keeping Events to Right)");
    //else h_purity->SetTitle("Purity - (Keeping Events to Left)");
    //h_purity->GetYaxis()->SetTitle("Purity (Captured Signal / Captured Events)");

    //if(keepEventstoRight) h_efficiency->SetTitle("Efficiency - (Keeping Events to Right)");
    //else h_efficiency->SetTitle("Efficiency - (Keeping Events to Left)");
    //h_efficiency->GetYaxis()->SetTitle("Efficiency (Captured Signal / Captured Events)");
    //
    //if(keepEventstoRight) h_purity_efficiency->SetTitle("Purity x Efficiency - (Keeping Events to Right)");
    //else h_purity_efficiency->SetTitle("Purity x Efficiency - (Keeping Events to Left)");
    //h_purity_efficiency->GetYaxis()->SetTitle("Purity x Efficiency");

    //// ------------------------------------------------------------------------
    //// Calculate Statistics
    //// ------------------------------------------------------------------------
    //double nSignal;
    //double nBackground;
    //double nSignalTotal = 0;
    //double nCapturedSignal = 0;
    //double nCapturedEvents = 0;
    //double purity;
    //double efficiency;
    //int nBins = h_signal->GetNbinsX();

    //// Get Total Number of Signal
    //for(int i = nBins; i > 0; i-- ){
    //    nSignalTotal = nSignalTotal +  h_signal->GetBinContent(i); 
    //}

    //// Keep Events to Right - Start Collecting Statistics from Right
    //// Keep Events to Left - Start Collecting Statistics from Left
    //int i;
    //if (keepEventstoRight) i = nBins;
    //else i = 1;
    //while(true){
    //    if(keepEventstoRight && i == 0) break;
    //    if(!keepEventstoRight && i > nBins) break;

    //    // Get Current Bin's Content
    //    nSignal = h_signal->GetBinContent(i);
    //    nBackground = h_background->GetBinContent(i);

    //    // Number of Events up-to this point
    //    nCapturedSignal = nCapturedSignal + nSignal;
    //    nCapturedEvents = nCapturedEvents + nSignal + nBackground;

    //    // Calculate Purity
    //    if( nCapturedEvents == 0 ) purity = 0;
    //    else purity = nCapturedSignal / nCapturedEvents;

    //    // Calculate Efficiency
    //    efficiency =  nCapturedSignal / nSignalTotal;

    //    // Fill Histograms
    //    h_purity->SetBinContent(i,purity);
    //    h_efficiency->SetBinContent(i,efficiency);
    //    h_purity_efficiency->SetBinContent(i,purity*efficiency);
    //
    //    if(keepEventstoRight) i--;
    //    else i++;
    //}

    //// Get Max(purity x efficiency)
    //double value_max = h_purity_efficiency->GetMaximum();
    //int bin_max = h_purity_efficiency->GetMaximumBin();
    //double best_cut = h_purity_efficiency->GetBinLowEdge(bin_max);

    //std::string x_title = h_purity_efficiency->GetXaxis()->GetTitle();
    //
    //h_purity_efficiency->GetXaxis()->SetTitle(Form("%s%s%f",x_title.c_str()," Best Cut = ",best_cut));
    //h_purity_efficiency->GetYaxis()->SetTitle(Form("%s%f","Purity x Efficiency, Max = ",value_max));
    //

    //// ------------------------------------------------------------------------
    //// Plot Histograms
    //// ------------------------------------------------------------------------
    //std::string tag;
    //if (keepEventstoRight) tag = "Right_";
    //else tag = "Left_";
    //std::string purity_fName = tag + "purity_" + fileName;
    //std::string efficiency_fName = tag + "efficiency_" + fileName;
    //std::string purity_efficiency_fName = tag + "purity_efficiency_" + fileName;

    //Draw1DHist(h_purity,purity_fName,plotDir);
    //Draw1DHist(h_efficiency,efficiency_fName,plotDir);
    //Draw1DHist(h_purity_efficiency,purity_efficiency_fName,plotDir);

    //delete h_purity;
    //delete h_efficiency;
    //delete h_purity_efficiency;

//}


#endif





