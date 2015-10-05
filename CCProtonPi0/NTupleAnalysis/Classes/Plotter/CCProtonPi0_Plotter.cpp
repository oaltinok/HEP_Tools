/*
    See CCProtonPi0_Plotter.h header for Class Information
*/

#ifndef CCProtonPi0_Plotter_cpp
#define CCProtonPi0_Plotter_cpp

#include "CCProtonPi0_Plotter.h"

using namespace PlotUtils;

void CCProtonPi0_Plotter::plotHistograms()
{
    /*
     *  Plot Functions for CCProtonPi0 Package
     */
    //plotInteraction();
    //plotMuon();
    //plotProton();
    plotPion();
    //plotCutHistograms();
    //plotPi0Blob();
    //SavePi0InvMassPoints();

    /*
     *  Plot Function Reserved for Other Studies
     */
    //plotOtherStudies();
}

// Todo - will revise function
void CCProtonPi0_Plotter::getPOT()
{
    std::string playlist = "Input/Playlists/pl_Data_All.dat"; 
    POTCounter pot_counter;
    double total_pot = pot_counter.getPOTfromPlaylist(playlist);

    std::cout<<"Total POT = "<<total_pot<<std::endl;

    data_POT = 9.61E19 + 2.28E20;
    mc_POT = 9.89E20 + 1.31E21;
    POT_Ratio_data_mc = data_POT / mc_POT;
    std::cout<<"POT Data = "<<data_POT<<std::endl;
    std::cout<<"POT MC = "<<mc_POT<<std::endl;
    std::cout<<"POT_Data / POT_MC = "<<POT_Ratio_data_mc<<std::endl;
 
}

CCProtonPi0_Plotter::CCProtonPi0_Plotter(std::string ana_folder)
{
    setRootDirs(ana_folder); 
    setPlotDirs(ana_folder);
}

void CCProtonPi0_Plotter::plotPi0Blob()
{
    std::cout<<"Plotting Pi0Blob"<<std::endl;
    std::string plotDir = plotDir_Pion; // Plots go under Pion
    
//    // Plot Truth Match - gamma 1
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"g1_evis_most_pdg",plotDir);
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"g1_evis_total_truth",plotDir);
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"g1_evis_frac_pizero",plotDir);
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"g1_evis_frac_piplus",plotDir);
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"g1_evis_frac_piminus",plotDir);
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"g1_evis_frac_proton",plotDir);
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"g1_evis_frac_neutron",plotDir);
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"g1_evis_frac_muon",plotDir);
//
//    // Plot Truth Match - gamma 2
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"g2_evis_most_pdg",plotDir);
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"g2_evis_total_truth",plotDir);
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"g2_evis_frac_pizero",plotDir);
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"g2_evis_frac_piplus",plotDir);
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"g2_evis_frac_piminus",plotDir);
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"g2_evis_frac_proton",plotDir);
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"g2_evis_frac_neutron",plotDir);
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"g2_evis_frac_muon",plotDir);
//    
//    // Total Pi0 Evis
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"captured_evis_frac_all",plotDir);
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"captured_evis_frac_signal",plotDir);


//    // Blob 1 Cluster Z and Strip Numbers
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"g1_cluster_Z",plotDir);
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"g1_max_cluster_Z",plotDir);
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"g1_min_cluster_Z",plotDir);
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"g1_strips",plotDir);
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"g1_max_strip",plotDir);
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"g1_min_strip",plotDir);
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"g1_nPlanes",plotDir);
//    // Blob 2 Cluster Z and Strip Numbers
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"g2_cluster_Z",plotDir);
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"g2_max_cluster_Z",plotDir);
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"g2_min_cluster_Z",plotDir);
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"g2_strips",plotDir);
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"g2_max_strip",plotDir);
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"g2_min_strip",plotDir);
//    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"g2_nPlanes",plotDir);

    std::cout<<"Plotting Pi0Blob Finished!"<<std::endl;
}


void CCProtonPi0_Plotter::plotOtherStudies()
{
    std::cout<<"Plotting Other Studies..."<<std::endl;
    
    std::string plotDir = plotDir_OtherStudies;

    Draw1DHist(rootDir_OtherStudies,"error",plotDir);
    Draw1DHist(rootDir_OtherStudies,"reco_energy",plotDir);
    Draw1DHist(rootDir_OtherStudies,"true_energy",plotDir);
    Draw2DHist(rootDir_OtherStudies,"reco_true_energy",plotDir);
    Draw2DHist(rootDir_OtherStudies,"true_recotrue_energy",plotDir);
    Draw2DHist(rootDir_OtherStudies,"evis_evis_ratio",plotDir,1);
    Draw2DHist(rootDir_OtherStudies,"true_evis_ratio",plotDir,1);
    Draw2DHist(rootDir_OtherStudies,"evis_true",plotDir,1);
    
    std::cout<<"Plotting Other Studies Finished!"<<std::endl;
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

    DrawStackedMC_BckgAll(rootDir_Interaction,"Enu_True",plotDir);
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
 
    Draw1DHist(rootDir_Interaction,"Enu_1Track_Error",plotDir);
    Draw1DHist(rootDir_Interaction,"Enu_2Track_Error",plotDir);
    Draw1DHist(rootDir_Interaction,"Enu_Cal_Error",plotDir);
   
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
    //plotStandardHistograms(rootDir_Pion, plotDir);
    
    //std::cout<<">> Plotting Unique Histograms"<<std::endl;
    // Unique Plots
    //DrawDataMC(rootDir_Pion,"gamma1_E",plotDir);
    //DrawDataMC(rootDir_Pion,"gamma2_E",plotDir);
    //DrawDataMC(rootDir_Pion,"gamma1_theta",plotDir);
    //DrawDataMC(rootDir_Pion,"gamma2_theta",plotDir);
    //DrawDataMC(rootDir_Pion,"gamma1_ConvLength",plotDir);
    //DrawDataMC(rootDir_Pion,"gamma2_ConvLength",plotDir);
    //DrawDataMC(rootDir_Pion,"photonEnergy_Asymmetry",plotDir);
    //DrawDataMC(rootDir_Pion,"invMass",plotDir);

    //DrawDataStackedMC(rootDir_Pion,"gamma1_E",plotDir);
    //DrawDataStackedMC(rootDir_Pion,"gamma2_E",plotDir);
    //DrawDataStackedMC(rootDir_Pion,"gamma1_theta",plotDir);
    //DrawDataStackedMC(rootDir_Pion,"gamma2_theta",plotDir);

  //  DrawStackedMC_BckgAll(rootDir_Pion,"gamma1_E",plotDir);
  //  DrawStackedMC_BckgAll(rootDir_Pion,"gamma2_E",plotDir);
  //  DrawStackedMC_BckgAll(rootDir_Pion,"gamma1_theta",plotDir);
  //  DrawStackedMC_BckgAll(rootDir_Pion,"gamma2_theta",plotDir);
  //  DrawStackedMC_BckgAll(rootDir_Pion,"gamma1_ConvLength",plotDir);
  //  DrawStackedMC_BckgAll(rootDir_Pion,"gamma2_ConvLength",plotDir);
  //  DrawStackedMC_BckgAll(rootDir_Pion,"photonEnergy_Asymmetry",plotDir);
  //  DrawStackedMC_BckgAll(rootDir_Pion,"invMass",plotDir);
    
    // ------------------------------------------------------------------------
    // Other Studies
    // ------------------------------------------------------------------------
    Draw1DHist(rootDir_Pion,"mgg_reco",plotDir);
    Draw1DHist(rootDir_Pion,"mgg_reco_error",plotDir);
    Draw1DHist(rootDir_Pion,"mgg_calc",plotDir);
    Draw1DHist(rootDir_Pion,"mgg_calc_error",plotDir);
    //
    Draw1DHist(rootDir_Pion,"gamma1_true_E",plotDir);
    //Draw1DHist(rootDir_Pion,"gamma1_evis_hcal",plotDir);
    //Draw1DHist(rootDir_Pion,"gamma1_evis_scal",plotDir);
    Draw1DHist(rootDir_Pion,"gamma1_reco_error_E",plotDir);
    Draw1DHist(rootDir_Pion,"gamma1_calc_error_E",plotDir);
    //Draw2DHist(rootDir_Pion,"gamma1_true_E_reco_E_ratio",plotDir);
    //Draw2DHist(rootDir_Pion,"gamma1_true_E_calc_E_ratio",plotDir);
    //Draw2DHist(rootDir_Pion,"gamma1_calc_E_calc_E_ratio",plotDir);
    //Draw2DHist(rootDir_Pion,"gamma1_evis_calc_E_ratio",plotDir);
    Draw2DHist(rootDir_Pion,"gamma1_evis_evis_ratio",plotDir,1);
   // Draw2DHist(rootDir_Pion,"gamma1_reco_E_true_E",plotDir);
   // Draw2DHist(rootDir_Pion,"gamma1_reco_E_calc_E",plotDir);
    Draw2DHist(rootDir_Pion,"gamma1_calc_E_true_E",plotDir);
 
    ofstream text;
    std::string textFile = plotDir_Pion + "evis.txt";
    text.open(textFile.c_str());

    // Get Histogram
    std::string root_dir = rootDir_Pion.mc;
    TFile* f = new TFile(root_dir.c_str());
    TH2D* hist2D = (TH2D*)f->Get("gamma1_evis_evis_ratio");

    // Reset Bins below the threshold
    int nBinsX = hist2D->GetNbinsX();
    int nBinsY = hist2D->GetNbinsY();

    for (int xBin = 1; xBin <= nBinsX; xBin++ ){
        double sum = 0;
        double n = 0;
        for (int yBin = 1; yBin <=nBinsY; yBin++){
            double nEvents = hist2D->GetBinContent(xBin,yBin);
            double y = hist2D->GetYaxis()->GetBinCenter(yBin);
            if ( y <= 2.5 && nEvents > 1){
                sum = sum + y*nEvents;
                n = n + nEvents;
            }
        }
    
        if (n > 0){
            double avg = sum / n;
            double x = hist2D->GetXaxis()->GetBinCenter(xBin);
            text<<x<<" "<<avg<<std::endl;
        }
    }
    
    text.close();
    
    Draw1DHist(rootDir_Pion,"gamma2_true_E",plotDir);
    Draw1DHist(rootDir_Pion,"gamma2_reco_error_E",plotDir);
    Draw1DHist(rootDir_Pion,"gamma2_calc_error_E",plotDir);
    //Draw2DHist(rootDir_Pion,"gamma2_true_E_reco_E_ratio",plotDir);
    //Draw2DHist(rootDir_Pion,"gamma2_true_E_calc_E_ratio",plotDir);
    //Draw2DHist(rootDir_Pion,"gamma2_reco_E_true_E",plotDir);
    Draw2DHist(rootDir_Pion,"gamma2_reco_E_calc_E",plotDir);
    Draw2DHist(rootDir_Pion,"gamma2_calc_E_true_E",plotDir);

    //Draw2DHist(rootDir_Pion,"gamma1_E_gamma2_E",plotDir);
    //Draw2DHist(rootDir_Pion,"gamma1_convLength_gamma2_convLength",plotDir);
   // 
   // Draw1DHist(rootDir_Pion,"isGamma1_conv_inside",plotDir);
   // Draw1DHist(rootDir_Pion,"isGamma2_conv_inside",plotDir);

   // plotGammaEvis();
    //plotPi0TruthMatch();

    std::cout<<"Plotting Pion Finished!\n"<<std::endl;
}

void CCProtonPi0_Plotter::plotPi0TruthMatch()
{
    std::cout<<"Plotting Pi0 Truth Match"<<std::endl;
    std::string plotDir = plotDir_Pion;
    
    // Evis Fraction 
    DrawStackedMC_BckgAll(rootDir_Pion,"evis_frac_reco_pi0_true_pi0",plotDir);
    DrawStackedMC_BckgAll(rootDir_Pion,"evis_frac_true_pi0_reco_all",plotDir);
    DrawStackedMC_BckgAll(rootDir_Pion,"evis_frac_reco_pi0_reco_all",plotDir);
    DrawStackedMC_BckgAll(rootDir_Pion,"evis_frac_reco_nonpi0_reco_all",plotDir);
    
    // Evis By Particle
    DrawStackedMC_GammaByPDG(rootDir_Pion,"evis",1,plotDir);
    DrawStackedMC_GammaByPDG(rootDir_Pion,"evis",2,plotDir);
    DrawStackedMC_GammaByPDG(rootDir_Pion,"evis",3,plotDir);

    // MC Hit Energy by PDG
    Draw1DHist(rootDir_Pion,"g1_digit_E_all",plotDir);
    Draw1DHist(rootDir_Pion,"g1_digit_E_pi0",plotDir);
    Draw1DHist(rootDir_Pion,"g1_digit_E_pi",plotDir);
    Draw1DHist(rootDir_Pion,"g1_digit_E_proton",plotDir);
    Draw1DHist(rootDir_Pion,"g1_digit_E_neutron",plotDir);
    Draw1DHist(rootDir_Pion,"g1_digit_E_muon",plotDir);
 
    Draw1DHist(rootDir_Pion,"g2_digit_E_all",plotDir);
    Draw1DHist(rootDir_Pion,"g2_digit_E_pi0",plotDir);
    Draw1DHist(rootDir_Pion,"g2_digit_E_pi",plotDir);
    Draw1DHist(rootDir_Pion,"g2_digit_E_proton",plotDir);
    Draw1DHist(rootDir_Pion,"g2_digit_E_neutron",plotDir);
    Draw1DHist(rootDir_Pion,"g2_digit_E_muon",plotDir);
   
    DrawStackedMC_GammaByPDG(rootDir_Pion,"digit_E",1,plotDir);
    DrawStackedMC_GammaByPDG(rootDir_Pion,"digit_E",2,plotDir);

    std::cout<<"Plotting Pi0 Truth Match Finished!\n"<<std::endl;
}

void CCProtonPi0_Plotter::plotGammaEvis()
{
    std::cout<<"Plotting Gamma Evis"<<std::endl;
    std::string plotDir = plotDir_Pion;

    DrawStackedMC_BckgAll(rootDir_Pion,"g2_evis_frac_scal_trkr",plotDir);
    
    DrawStackedMC_GammaEvis(rootDir_Pion,1,plotDir);
    DrawStackedMC_GammaEvis(rootDir_Pion,2,plotDir);

    std::cout<<"Plotting Gamma Evis Finished!\n"<<std::endl;
}
void CCProtonPi0_Plotter::plotStandardHistograms(rootDir &dir, std::string plotDir)
{
    std::cout<<">> Plotting Standard Histograms"<<std::endl;
    
    //DrawDataMC(dir, "E", plotDir);
    DrawDataMC(dir, "P", plotDir);
    //DrawDataMC(dir, "KE", plotDir);
    DrawDataMC(dir, "theta", plotDir);
    DrawDataMC(dir, "phi", plotDir);

    //DrawDataStackedMC(dir, "E", plotDir);
    DrawDataStackedMC(dir, "P", plotDir);
    //DrawDataStackedMC(dir, "KE", plotDir);
    DrawDataStackedMC(dir, "theta", plotDir);
    DrawDataStackedMC(dir, "phi", plotDir);
}


void CCProtonPi0_Plotter::plotCutHistograms()
{
    std::cout<<"Plotting CutHistograms"<<std::endl;
    std::string plotDir = plotDir_CutHists;
    
    // ------------------------------------------------------------------------
    // 1 Track
    // ------------------------------------------------------------------------
    CutArrow Michel_1Track(1,0,60E3,0.1,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_Michel",plotDir,1, Michel_1Track);

    CutArrow eVis_nuclearTarget_1Track(20,0,50E3,1,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_eVis_nuclearTarget",plotDir, 1, eVis_nuclearTarget_1Track);
    
    CutArrow eVis_other_min_1Track(50,0,12E3,100,"R"); 
    CutArrow eVis_other_max_1Track(2000,0,12E3,100,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_eVis_other",plotDir, 2, eVis_other_min_1Track, eVis_other_max_1Track);
    
    CutArrow gamma1_ConvDist_1Track(15,0,600,5,"R"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_gamma1ConvDist",plotDir, 1, gamma1_ConvDist_1Track);
    
    CutArrow gamma2_ConvDist_1Track(15,0,300,5,"R"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_gamma2ConvDist",plotDir, 1, gamma2_ConvDist_1Track);
    
    CutArrow pi0invMass_min_1Track(75,0,300,20,"R"); 
    CutArrow pi0invMass_max_1Track(195,0,300,20,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_pi0invMass",plotDir, 2, pi0invMass_min_1Track, pi0invMass_max_1Track);
    
    CutArrow neutrinoE_1Track(20,0,190,2,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_neutrinoE",plotDir, 1, neutrinoE_1Track);
    
    CutArrow unusedE_1Track(300,0,350,50,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_UnusedE",plotDir, 1, unusedE_1Track);
   
    // ------------------------------------------------------------------------
    // 2 Track
    // ------------------------------------------------------------------------
    CutArrow Michel_2Track(1,0,60E3,0.1,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_Michel",plotDir,1, Michel_2Track);

    CutArrow eVis_nuclearTarget_2Track(20,0,50E3,1,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_eVis_nuclearTarget",plotDir, 1, eVis_nuclearTarget_2Track);
    
    CutArrow eVis_other_min_2Track(50,0,12E3,100,"R"); 
    CutArrow eVis_other_max_2Track(2000,0,12E3,100,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_eVis_other",plotDir, 2, eVis_other_min_2Track, eVis_other_max_2Track);
    
    CutArrow gamma1_ConvDist_2Track(15,0,600,5,"R"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_gamma1ConvDist",plotDir, 1, gamma1_ConvDist_2Track);
    
    CutArrow gamma2_ConvDist_2Track(15,0,300,5,"R"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_gamma2ConvDist",plotDir, 1, gamma2_ConvDist_2Track);
    
    CutArrow pi0invMass_min_2Track(75,0,300,20,"R"); 
    CutArrow pi0invMass_max_2Track(195,0,300,20,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_pi0invMass",plotDir, 2, pi0invMass_min_2Track, pi0invMass_max_2Track);

    CutArrow neutrinoE_2Track(20,0,170,2,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_neutrinoE",plotDir, 1, neutrinoE_2Track);
    
    CutArrow unusedE_2Track(300,0,325,50,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_UnusedE",plotDir, 1, unusedE_2Track);
    
    CutArrow protonScore_pIDDiff(0.45,0,25,0.1,"R"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_protonScore_pIDDiff",plotDir, 1, protonScore_pIDDiff);
    
    CutArrow protonScore_LLR(10,0,140,10,"R"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_protonScore_LLR",plotDir, 1, protonScore_LLR);
    
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_deltaInvMass",plotDir);
    
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

void CCProtonPi0_Plotter::setRootDirs(std::string ana_folder)
{
    // Set Other Studies ROOT Dir
    rootDir_OtherStudies.mc = "/minerva/data/users/oaltinok/NTupleAnalysis/ParticleCannon/PC_Test.root"; 
    rootDir_OtherStudies.data = "";
    
    // Set MC Root Dir
    rootDir_CutHists.mc = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + "CutHistograms.root";
    rootDir_Interaction.mc = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + ana_folder + "Interaction.root";
    rootDir_Muon.mc = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + ana_folder + "Muon.root";
    rootDir_Proton.mc = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + ana_folder + "Proton.root";
    rootDir_Pion.mc = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + ana_folder + "Pion.root";
    rootDir_Pi0Blob.mc = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + ana_folder + "Pi0Blob.root";

    // Set Data Root Dir
    rootDir_CutHists.data = Folder_List::rootOut + Folder_List::Data + Folder_List::analyzed + "CutHistograms.root";
    rootDir_Interaction.data = Folder_List::rootOut + Folder_List::Data + Folder_List::analyzed + ana_folder + "Interaction.root";
    rootDir_Muon.data = Folder_List::rootOut + Folder_List::Data + Folder_List::analyzed + ana_folder + "Muon.root";
    rootDir_Proton.data = Folder_List::rootOut + Folder_List::Data + Folder_List::analyzed + ana_folder + "Proton.root";
    rootDir_Pion.data = Folder_List::rootOut + Folder_List::Data + Folder_List::analyzed + ana_folder + "Pion.root";
    rootDir_Pi0Blob.data = Folder_List::rootOut + Folder_List::Data + Folder_List::analyzed + ana_folder + "Pi0Blob.root";
}

void CCProtonPi0_Plotter::setPlotDirs(std::string ana_folder)
{
    plotDir_OtherStudies = Folder_List::output + Folder_List::plotOut + "OtherStudies/";
    plotDir_CutHists = Folder_List::output + Folder_List::plotOut + "CutHists/";
    plotDir_Interaction = Folder_List::output + Folder_List::plotOut + ana_folder + "Interaction/";
    plotDir_Muon = Folder_List::output + Folder_List::plotOut + ana_folder + "Muon/";
    plotDir_Proton = Folder_List::output + Folder_List::plotOut + ana_folder + "Proton/";
    plotDir_Pion = Folder_List::output + Folder_List::plotOut + ana_folder + "Pion/";
    plotDir_Other = Folder_List::output + Folder_List::plotOut + ana_folder + "Other/";
}

void CCProtonPi0_Plotter::SavePi0InvMassPoints()
{
    std::cout<<"Saving Pi0 Invariant Mass Points"<<std::endl;

    std::string rootDir[2];
    std::string textDir[2];
    std::ofstream text[2];

    // Set Dirs
    textDir[0] = Folder_List::output + Folder_List::textOut + "Pi0InvMass_MC_minerva13.txt";
    rootDir[0] = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + "CutHistograms.root";

    textDir[1] = Folder_List::output + Folder_List::textOut + "Pi0InvMass_Data_minerva13.txt";
    rootDir[1] = Folder_List::rootOut + Folder_List::Data + Folder_List::analyzed + "CutHistograms.root";
    
    // Open Output Text Files
    for (int i = 0; i < 2; i++){
        text[i].open( textDir[i].c_str() );
        if( !text[i].is_open() ){
            std::cerr<<"Cannot open output text file: "<<textDir[i]<<std::endl;
            exit(1);
        }else{
            std::cout<<"\t"<<textDir[i]<<std::endl;
        }
    }
   
    // Open Input ROOT Files and Initialize Histograms
    TFile* f_mc = new TFile(rootDir[0].c_str());
    MnvH1D* h_1Track_MC_Signal = (MnvH1D*)f_mc->Get("hCut_1Track_pi0invMass_1");
    MnvH1D* h_2Track_MC_Signal = (MnvH1D*)f_mc->Get("hCut_2Track_pi0invMass_1");
    MnvH1D* h_1Track_MC_Background = (MnvH1D*)f_mc->Get("hCut_1Track_pi0invMass_2");
    MnvH1D* h_2Track_MC_Background = (MnvH1D*)f_mc->Get("hCut_2Track_pi0invMass_2");
    
    TFile* f_data = new TFile(rootDir[1].c_str());
    MnvH1D* h_1Track_Data_All = (MnvH1D*)f_data->Get("hCut_1Track_pi0invMass_0");
    MnvH1D* h_2Track_Data_All = (MnvH1D*)f_data->Get("hCut_2Track_pi0invMass_0");
    
    
    // ------------------------------------------------------------------------
    // Save Histogram Information to Output Text
    // ------------------------------------------------------------------------
    double nBins = h_1Track_MC_Signal->GetNbinsX(); // nBins Are Same for All
       
    // Write MC File
    for (int i = 1; i <= nBins; i++){
        double bin_center = h_1Track_MC_Signal->GetBinCenter(i);
        
        double bin_content_1Track_Signal = h_1Track_MC_Signal->GetBinContent(i);
        double bin_content_1Track_Bckg = h_1Track_MC_Background->GetBinContent(i);
        
        double bin_content_2Track_Signal = h_2Track_MC_Signal->GetBinContent(i);
        double bin_content_2Track_Bckg = h_2Track_MC_Background->GetBinContent(i);

        text[0] <<bin_center<<" "
                <<bin_content_1Track_Signal<<" "
                <<bin_content_1Track_Bckg<<" "
                <<bin_content_2Track_Signal<<" "
                <<bin_content_2Track_Bckg<<std::endl;
    }
    
    // Write Data File
    for (int i = 1; i <= nBins; i++){
        double bin_center = h_1Track_Data_All->GetBinCenter(i);
        
        double bin_content_1Track = h_1Track_Data_All->GetBinContent(i);
        double bin_content_2Track = h_2Track_Data_All->GetBinContent(i);
        
        text[1] <<bin_center<<" "
                <<bin_content_1Track<<" "
                <<bin_content_2Track<<std::endl;
    }

    for (int i = 0; i < 2; i++){
        text[0].close();
    }
    std::cout<<"Done!"<<std::endl;
}

#endif

