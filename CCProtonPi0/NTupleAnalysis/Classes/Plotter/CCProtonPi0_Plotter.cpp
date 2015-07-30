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
}


CCProtonPi0_Plotter::CCProtonPi0_Plotter()
{
    data_POT = 9.57E19;
    mc_POT = 9.90E20;
    POT_Ratio_data_mc = data_POT / mc_POT;
    std::cout<<"POT Data = "<<data_POT<<std::endl;
    std::cout<<"POT MC = "<<mc_POT<<std::endl;
    std::cout<<"POT_Data / POT_MC = "<<POT_Ratio_data_mc<<std::endl;
    
    setRootDirs(); 
    setPlotDirs();
}

void CCProtonPi0_Plotter::plotInteraction()
{
    std::cout<<"Plotting Interaction"<<std::endl;
    std::string plotDir = plotDir_Interaction;

    plot_mc_w_Stacked();
    plot_final_mc_w_Stacked();

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

    DrawDataStackedMC(rootDir_Interaction,"Enu_1Track",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"Enu_2Track",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"Enu_Cal",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"q2",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"w",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"wSq",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"E_Unused_afterReco",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"E_Used_afterReco",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"deltaInvMass",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"nProngs_hist",plotDir);
    
    std::cout<<"Plotting Interaction Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plotMuon()
{
    std::cout<<"Plotting Muon"<<std::endl;
    std::string plotDir = plotDir_Muon;
    plotStandardHistograms(rootDir_Muon, plotDir);
    std::cout<<"Plotting Muon Finished!\n"<<std::endl;
}

void CCProtonPi0_Plotter::plotProton()
{    
    std::cout<<"Plotting Proton"<<std::endl;
    std::string plotDir = plotDir_Proton;
    
    plotStandardHistograms(rootDir_Proton, plotDir);
   
    std::cout<<">> Plotting Unique Histograms"<<std::endl;
    // Unique Plots
    DrawDataMC(rootDir_Proton,"trackLength",plotDir);
    DrawDataMC(rootDir_Proton,"trackKinked",plotDir);
    DrawDataMC(rootDir_Proton,"partScore",plotDir);

    DrawDataStackedMC(rootDir_Proton,"trackLength",plotDir);
    DrawDataStackedMC(rootDir_Proton,"trackKinked",plotDir);
    DrawDataStackedMC(rootDir_Proton,"partScore",plotDir);

    Draw1DHist(rootDir_Proton,"P_error",plotDir);
    Draw2DHist(rootDir_Proton,"reco_P_true_P",plotDir);
    
    std::cout<<"Plotting Proton Finished!\n"<<std::endl;
}

void CCProtonPi0_Plotter::plotPion()
{
    std::cout<<"Plotting Pion"<<std::endl;
    std::string plotDir = plotDir_Pion;
    
    // Standard Plots
    plotStandardHistograms(rootDir_Pion, plotDir);
    
    std::cout<<">> Plotting Unique Histograms"<<std::endl;
    // Unique Plots
    DrawDataMC(rootDir_Pion,"gamma1_ConvLength",plotDir);
    DrawDataMC(rootDir_Pion,"gamma2_ConvLength",plotDir);
    DrawDataMC(rootDir_Pion,"photonEnergy_Asymmetry",plotDir);
    DrawDataMC(rootDir_Pion,"invMass",plotDir);

    DrawDataStackedMC(rootDir_Pion,"gamma1_ConvLength",plotDir);
    DrawDataStackedMC(rootDir_Pion,"gamma2_ConvLength",plotDir);
    DrawDataStackedMC(rootDir_Pion,"photonEnergy_Asymmetry",plotDir);
    DrawDataStackedMC(rootDir_Pion,"invMass",plotDir);

    Draw1DHist(rootDir_Pion,"gamma1_P_error",plotDir);
    Draw2DHist(rootDir_Pion,"gamma1_reco_P_true_P",plotDir);
    Draw1DHist(rootDir_Pion,"gamma2_P_error",plotDir);
    Draw2DHist(rootDir_Pion,"gamma2_reco_P_true_P",plotDir);
    Draw2DHist(rootDir_Pion,"gamma1_P_gamma2_P",plotDir);
    Draw2DHist(rootDir_Pion,"gamma1_convLength_gamma2_convLength",plotDir);
    
    std::cout<<"Plotting Pion Finished!\n"<<std::endl;
}

void CCProtonPi0_Plotter::plotStandardHistograms(rootDir &dir, std::string plotDir)
{
    std::cout<<">> Plotting Standard Histograms"<<std::endl;
    
    DrawDataMC(dir, "E", plotDir);
    DrawDataMC(dir, "P", plotDir);
    DrawDataMC(dir, "KE", plotDir);
    DrawDataMC(dir, "theta", plotDir);
    DrawDataMC(dir, "phi", plotDir);

    DrawDataStackedMC(dir, "E", plotDir);
    DrawDataStackedMC(dir, "P", plotDir);
    DrawDataStackedMC(dir, "KE", plotDir);
    DrawDataStackedMC(dir, "theta", plotDir);
    DrawDataStackedMC(dir, "phi", plotDir);
}


void CCProtonPi0_Plotter::plotCutHistograms()
{
    std::cout<<"Plotting CutHistograms"<<std::endl;
    std::string plotDir = plotDir_CutHists;
   
    CutArrow Michel(1,0,60E3,0.1,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_Michel",plotDir,1, Michel);

    CutArrow eVis_nuclearTarget(20,0,50E3,1,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_eVis_nuclearTarget",plotDir, 1, eVis_nuclearTarget);
    
    CutArrow eVis_other_min(50,0,12E3,100,"R"); 
    CutArrow eVis_other_max(2000,0,12E3,100,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_eVis_other",plotDir, 2, eVis_other_min, eVis_other_max);
    
    CutArrow gamma1_ConvDist(15,0,600,5,"R"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_gamma1ConvDist",plotDir, 1, gamma1_ConvDist);
    
    CutArrow gamma2_ConvDist(15,0,300,5,"R"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_gamma2ConvDist",plotDir, 1, gamma2_ConvDist);
    
    CutArrow pi0invMass_min(75,0,300,20,"R"); 
    CutArrow pi0invMass_max(195,0,300,20,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_pi0invMass",plotDir, 2, pi0invMass_min, pi0invMass_max);
    
    CutArrow neutrinoE_1Prong(20,0,190,2,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_1Prong_neutrinoE",plotDir, 1, neutrinoE_1Prong);
    
    CutArrow neutrinoE_2Prong(20,0,170,2,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Prong_neutrinoE",plotDir, 1, neutrinoE_2Prong);
    
    CutArrow unusedE_1Prong(300,0,350,50,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_1Prong_UnusedE",plotDir, 1, unusedE_1Prong);
    
    CutArrow unusedE_2Prong(300,0,325,50,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Prong_UnusedE",plotDir, 1, unusedE_2Prong);
    
    CutArrow protonScore_pIDDiff(0.45,0,25,0.1,"R"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_protonScore_pIDDiff",plotDir, 1, protonScore_pIDDiff);
    
    CutArrow protonScore_LLR(10,0,140,10,"R"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_protonScore_LLR",plotDir, 1, protonScore_LLR);
    
    DrawDataStackedMC(rootDir_CutHists,"hCut_deltaInvMass",plotDir);
    
    std::cout<<"Plotting CutHistograms Finished!"<<std::endl;
}
 
void CCProtonPi0_Plotter::plot_mc_w_Stacked()
{
    // mc_w written during reduce - Its Histogram is with Cut Hists
    std::string root_dir = rootDir_CutHists.mc;
    std::string plotDir = plotDir_Interaction;
    
    std::cout<<"\nPlottting Stacked mc_w"<<std::endl;
    
    TFile* f_Root = new TFile(root_dir.c_str());
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
    std::string root_dir = rootDir_Interaction.mc;
    std::string plotDir = plotDir_Interaction;
    
    std::cout<<"\nPlottting Stacked final_mc_w"<<std::endl;
    
    TFile* f_Root = new TFile(root_dir.c_str());
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

void CCProtonPi0_Plotter::setRootDirs()
{
    // Set MC Root Dir
    rootDir_CutHists.mc = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + "CutHistograms.root";
    rootDir_Interaction.mc = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + "Interaction.root";
    rootDir_Muon.mc = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + "Muon.root";
    rootDir_Proton.mc = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + "Proton.root";
    rootDir_Pion.mc = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + "Pion.root";

    // Set Data Root Dir
    rootDir_CutHists.data = Folder_List::rootOut + Folder_List::Data + Folder_List::analyzed + "CutHistograms.root";
    rootDir_Interaction.data = Folder_List::rootOut + Folder_List::Data + Folder_List::analyzed + "Interaction.root";
    rootDir_Muon.data = Folder_List::rootOut + Folder_List::Data + Folder_List::analyzed + "Muon.root";
    rootDir_Proton.data = Folder_List::rootOut + Folder_List::Data + Folder_List::analyzed + "Proton.root";
    rootDir_Pion.data = Folder_List::rootOut + Folder_List::Data + Folder_List::analyzed + "Pion.root";
}

void CCProtonPi0_Plotter::setPlotDirs()
{
    plotDir_Interaction = Folder_List::output + Folder_List::plotOut + "Interaction/";
    plotDir_CutHists = Folder_List::output + Folder_List::plotOut + "CutHists/";
    plotDir_Muon = Folder_List::output + Folder_List::plotOut + "Muon/";
    plotDir_Proton = Folder_List::output + Folder_List::plotOut + "Proton/";
    plotDir_Pion = Folder_List::output + Folder_List::plotOut + "Pion/";
    otherDir = Folder_List::output + Folder_List::plotOut + Folder_List::other;
}

#endif

