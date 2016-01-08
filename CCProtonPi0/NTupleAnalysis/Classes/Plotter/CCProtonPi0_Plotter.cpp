/*
 * See CCProtonPi0_Plotter.h header for Class Information
 */

#ifndef CCProtonPi0_Plotter_cpp
#define CCProtonPi0_Plotter_cpp

#include "CCProtonPi0_Plotter.h"

using namespace PlotUtils;

void CCProtonPi0_Plotter::plotHistograms()
{
    //--------------------------------------------------------------------------
    // Run only once to get the POT
    //--------------------------------------------------------------------------
    //getPOT_MC();
    //getPOT_Data();

    //--------------------------------------------------------------------------
    //  Data vs MC
    //--------------------------------------------------------------------------
    //plotInteraction_DataMC();
    //plotMuon_DataMC();
    //plotProton_DataMC();
    //plotPion_DataMC();
    //plotCutHistograms_DataMC();
    //plotPi0Blob_DataMC();

    //--------------------------------------------------------------------------
    //  MC Only
    //--------------------------------------------------------------------------
    //plotInteraction_MCOnly();
    //plotMuon_MCOnly();
    //plotProton_MCOnly();
    //plotPion_MCOnly();
    plotCutHistograms_MCOnly();
    //plotPi0Blob_MCOnly();
    //plotEfficiencyCurves();

    //--------------------------------------------------------------------------
    //  Plot Function Reserved for Other Studies
    //--------------------------------------------------------------------------
    //SavePi0InvMassPoints();
    //plotOtherStudies();
}

void CCProtonPi0_Plotter::getPOT_MC()
{
    std::string playlist = "Input/Playlists/pl_MC_All.dat"; 
    POTCounter pot_counter;
    double total_pot = pot_counter.getPOTfromPlaylist(playlist);

    std::cout<<"Total MC POT = "<<total_pot<<std::endl;
}

void CCProtonPi0_Plotter::getPOT_Data()
{
    std::string playlist = "Input/Playlists/pl_Data_All.dat"; 
    POTCounter pot_counter;
    double total_pot = pot_counter.getPOTfromPlaylist(playlist);

    std::cout<<"Total Data POT = "<<total_pot<<std::endl;
}

CCProtonPi0_Plotter::CCProtonPi0_Plotter(std::string ana_folder)
{
    //--------------------------------------------------------------------------
    // Set POT -- Run getPOT_MC() and getPOT_Data() Functions once to get POT
    //--------------------------------------------------------------------------
    data_POT = 9.58121e+19; 
    mc_POT = 9.40862e+20;
    POT_Ratio_data_mc = data_POT / mc_POT;
    std::cout<<"POT Data = "<<data_POT<<std::endl;
    std::cout<<"POT MC = "<<mc_POT<<std::endl;
    std::cout<<"POT Ratio = "<<POT_Ratio_data_mc<<std::endl;

    setRootDirs(ana_folder); 
    setPlotDirs(ana_folder);
}

void CCProtonPi0_Plotter::plotPi0Blob_MCOnly()
{
    std::cout<<"Plotting Pi0Blob MC Only"<<std::endl;
    std::string plotDir = plotDir_Pi0Blob; // Plots go under Pion

    // Evis Fraction 
    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"evis_frac_reco_pi0_true_pi0",plotDir);
    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"evis_frac_true_pi0_reco_all",plotDir);
    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"evis_frac_reco_pi0_reco_all",plotDir);
    DrawStackedMC_BckgAll(rootDir_Pi0Blob,"evis_frac_reco_nonpi0_reco_all",plotDir);

    // Evis By Particle
    DrawStackedMC_GammaByPDG(rootDir_Pi0Blob,"evis",1,plotDir);
    DrawStackedMC_GammaByPDG(rootDir_Pi0Blob,"evis",2,plotDir);
    DrawStackedMC_GammaByPDG(rootDir_Pi0Blob,"evis",3,plotDir);

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

    std::cout<<"Plotting Pi0Blob MC Only Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plotPi0Blob_DataMC()
{
    // Do Nothing 
}

void CCProtonPi0_Plotter::plotOtherStudies()
{
    std::cout<<"Plotting Other Studies..."<<std::endl;

    std::string plotDir = plotDir_OtherStudies;

    Draw1DHist(rootDir_OtherStudies,"evis",plotDir);
    Draw1DHist(rootDir_OtherStudies,"true_energy",plotDir);

    Draw2DHist(rootDir_OtherStudies,"evis_evis_ratio",plotDir,1);
    Draw2DHist(rootDir_OtherStudies,"true_evis_ratio",plotDir,1);
    Draw2DHist(rootDir_OtherStudies,"evis_true",plotDir,1);

    Draw1DHist(rootDir_OtherStudies,"error",plotDir);
    Draw1DHist(rootDir_OtherStudies,"reco_energy",plotDir);
    Draw2DHist(rootDir_OtherStudies,"reco_true_energy",plotDir);
    Draw2DHist(rootDir_OtherStudies,"true_recotrue_energy",plotDir);

    std::cout<<"Plotting Other Studies Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plotInteraction_MCOnly()
{
    std::cout<<"Plotting Interaction MC Only"<<std::endl;
    std::string plotDir = plotDir_Interaction;

    //plot_mc_w_Stacked();
    //plot_final_mc_w_Stacked();

    //DrawStackedMC_BckgAll(rootDir_Interaction,"vertex_energy_1Track",plotDir);
    //DrawStackedMC_BckgAll(rootDir_Interaction,"vertex_energy_2Track",plotDir);
    //DrawStackedMC_BckgAll(rootDir_Interaction,"vertex_evis_1Track",plotDir);
    //DrawStackedMC_BckgAll(rootDir_Interaction,"vertex_evis_2Track",plotDir);
    //DrawStackedMC_BckgAll(rootDir_Interaction,"extra_evis_1Track",plotDir);
    //DrawStackedMC_BckgAll(rootDir_Interaction,"extra_evis_2Track",plotDir);

    //Draw1DHist(rootDir_Interaction,"n_ejected_nucleons_1Track",plotDir);
    //Draw1DHist(rootDir_Interaction,"n_ejected_nucleons_2Track",plotDir);

    //Draw1DHist(rootDir_Interaction,"proton_true_P_1Track",plotDir);
    //Draw1DHist(rootDir_Interaction,"proton_true_KE_1Track",plotDir);

    Draw1DHist(rootDir_Interaction,"Enu_True_1Track",plotDir);
    Draw1DHist(rootDir_Interaction,"Enu_True_2Track",plotDir);
    Draw1DHist(rootDir_Interaction,"Enu_1Track_Error",plotDir);
    Draw1DHist(rootDir_Interaction,"Enu_2Track_Error",plotDir);
    //Draw1DHist(rootDir_Interaction,"Enu_1Track_Alt_Error",plotDir);
    
    Draw1DHist(rootDir_Interaction,"Enu_1Track_Diff",plotDir);
    Draw1DHist(rootDir_Interaction,"Enu_2Track_Diff",plotDir);
   
    std::cout<<"Plotting Interaction MC Only Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plotInteraction_DataMC()
{
    std::cout<<"Plotting Interaction Data vs MC"<<std::endl;
    std::string plotDir = plotDir_Interaction;

    DrawDataMC(rootDir_Interaction,"Enu_1Track",plotDir);
    //DrawDataMC(rootDir_Interaction,"Enu_1Track_Alt",plotDir);
    DrawDataMC(rootDir_Interaction,"Enu_2Track",plotDir);
    DrawDataMC(rootDir_Interaction,"Enu",plotDir);
    DrawDataMC(rootDir_Interaction,"QSq",plotDir);
    DrawDataMC(rootDir_Interaction,"WSq",plotDir);
    DrawDataMC(rootDir_Interaction,"W",plotDir);
    DrawDataMC(rootDir_Interaction,"deltaInvMass",plotDir);

    DrawDataStackedMC(rootDir_Interaction,"Enu_1Track",plotDir);
    //DrawDataStackedMC(rootDir_Interaction,"Enu_1Track_Alt",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"Enu_2Track",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"Enu",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"QSq",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"WSq",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"W",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"deltaInvMass",plotDir);

    DrawDataStackedMC(rootDir_Interaction,"vertex_energy_1Track",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"vertex_energy_2Track",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"vertex_evis_1Track",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"vertex_evis_2Track",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"extra_evis_1Track",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"extra_evis_2Track",plotDir);

    std::cout<<"Plotting Interaction Data vs MC Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plotMuon_MCOnly()
{
    std::cout<<"Plotting Muon MC Only"<<std::endl;
    std::string plotDir = plotDir_Muon;

    Draw1DHist(rootDir_Muon,"P_error",plotDir);
    Draw1DHist(rootDir_Muon,"E_error",plotDir);
    Draw1DHist(rootDir_Muon,"E_Diff",plotDir);

    std::cout<<"Plotting Muon MC Only Finished!\n"<<std::endl;

}

void CCProtonPi0_Plotter::plotMuon_DataMC()
{
    std::cout<<"Plotting Muon Data vs MC"<<std::endl;
    std::string plotDir = plotDir_Muon;

    plotStandardHistograms(rootDir_Muon, plotDir);

    std::cout<<"Plotting Muon Data vs MC Finished!\n"<<std::endl;
}

void CCProtonPi0_Plotter::plotProton_MCOnly()
{
    std::cout<<"Plotting Proton MC Only"<<std::endl;
    std::string plotDir = plotDir_Proton;

    DrawStackedMC_BckgAll(rootDir_Proton, "P", plotDir);
    DrawStackedMC_BckgAll(rootDir_Proton, "theta", plotDir);
    DrawStackedMC_BckgAll(rootDir_Proton, "phi", plotDir);
    
    Draw1DHist(rootDir_Proton,"P_error",plotDir);
    Draw1DHist(rootDir_Proton,"E_error",plotDir);
    Draw1DHist(rootDir_Proton,"E_Diff",plotDir);

    std::cout<<"Plotting Proton MC Only Finished!\n"<<std::endl;
}


void CCProtonPi0_Plotter::plotProton_DataMC()
{    
    std::cout<<"Plotting Proton Data vs MC"<<std::endl;
    std::string plotDir = plotDir_Proton;

    plotStandardHistograms(rootDir_Proton, plotDir);

    std::cout<<">> Plotting Unique Histograms"<<std::endl;

    //DrawDataMC(rootDir_Proton,"trackLength",plotDir);
    //DrawDataMC(rootDir_Proton,"trackKinked",plotDir);
    //DrawDataMC(rootDir_Proton,"partScore",plotDir);

    DrawDataStackedMC(rootDir_Proton,"trackLength",plotDir);
    DrawDataStackedMC(rootDir_Proton,"trackKinked",plotDir);
    DrawDataStackedMC(rootDir_Proton,"partScore",plotDir);

    std::cout<<"Plotting Proton Data vs MC Finished!\n"<<std::endl;
}

void CCProtonPi0_Plotter::plotPion_MCOnly()
{
    std::cout<<"Plotting Pion MC Only"<<std::endl;
    std::string plotDir = plotDir_Pion;

    plotPion_True();

    DrawStackedMC_BckgAll(rootDir_Pion,"gamma1_E",plotDir);
    DrawStackedMC_BckgAll(rootDir_Pion,"gamma1_theta",plotDir);
    DrawStackedMC_BckgAll(rootDir_Pion,"gamma1_ConvLength",plotDir);

    DrawStackedMC_BckgAll(rootDir_Pion,"gamma2_E",plotDir);
    DrawStackedMC_BckgAll(rootDir_Pion,"gamma2_theta",plotDir);
    DrawStackedMC_BckgAll(rootDir_Pion,"gamma2_ConvLength",plotDir);

    DrawStackedMC_BckgAll(rootDir_Pion,"photonEnergy_Asymmetry",plotDir);
    DrawStackedMC_BckgAll(rootDir_Pion,"invMass",plotDir);
    DrawStackedMC_BckgAll(rootDir_Pion,"P",plotDir);
    DrawStackedMC_BckgAll(rootDir_Pion,"E",plotDir);
    DrawStackedMC_BckgAll(rootDir_Pion,"theta",plotDir);
    DrawStackedMC_BckgAll(rootDir_Pion,"phi",plotDir);

    std::cout<<"Plotting Pion MC Only Finished!\n"<<std::endl;
}


void CCProtonPi0_Plotter::plotPion_DataMC()
{
    std::cout<<"Plotting Pion Data vs MC"<<std::endl;
    std::string plotDir = plotDir_Pion;

    plotStandardHistograms(rootDir_Pion, plotDir);

    std::cout<<">> Plotting Unique Histograms"<<std::endl;

    //DrawDataMC(rootDir_Pion,"gamma1_E",plotDir);
    //DrawDataMC(rootDir_Pion,"gamma1_theta",plotDir);
    //DrawDataMC(rootDir_Pion,"gamma1_ConvLength",plotDir);

    //DrawDataMC(rootDir_Pion,"gamma2_E",plotDir);
    //DrawDataMC(rootDir_Pion,"gamma2_theta",plotDir);
    //DrawDataMC(rootDir_Pion,"gamma2_ConvLength",plotDir);

    //DrawDataMC(rootDir_Pion,"invMass",plotDir);
    //DrawDataMC(rootDir_Pion,"photonEnergy_Asymmetry",plotDir);

    DrawDataStackedMC(rootDir_Pion,"gamma1_E",plotDir);
    DrawDataStackedMC(rootDir_Pion,"gamma1_theta",plotDir);
    DrawDataStackedMC(rootDir_Pion,"gamma1_ConvLength",plotDir);

    DrawDataStackedMC(rootDir_Pion,"gamma2_E",plotDir);
    DrawDataStackedMC(rootDir_Pion,"gamma2_theta",plotDir);
    DrawDataStackedMC(rootDir_Pion,"gamma2_ConvLength",plotDir);

    DrawDataStackedMC(rootDir_Pion,"invMass",plotDir);
    DrawDataStackedMC(rootDir_Pion,"photonEnergy_Asymmetry",plotDir);

    std::cout<<"Plotting Pion Data vs MC Finished!\n"<<std::endl;
}

void CCProtonPi0_Plotter::plotPion_True()
{
    std::cout<<"Plotting Pion True"<<std::endl;
    std::string plotDir = plotDir_Pion;

    Draw1DHist(rootDir_Pion,"gamma1_true_E",plotDir);
    Draw1DHist(rootDir_Pion,"gamma1_reco_error_E",plotDir);
    Draw2DHist(rootDir_Pion,"gamma1_true_E_reco_E_error",plotDir);
    Draw2DHist(rootDir_Pion,"gamma1_reco_E_true_E",plotDir);

    Draw1DHist(rootDir_Pion,"gamma2_true_E",plotDir);
    Draw1DHist(rootDir_Pion,"gamma2_reco_error_E",plotDir);
    Draw2DHist(rootDir_Pion,"gamma2_true_E_reco_E_error",plotDir);
    Draw2DHist(rootDir_Pion,"gamma2_reco_E_true_E",plotDir);

    Draw1DHist(rootDir_Pion,"P_error",plotDir);
    Draw1DHist(rootDir_Pion,"E_error",plotDir);
    Draw1DHist(rootDir_Pion,"E_Diff",plotDir);

    std::cout<<"Plotting Pion True Finished!\n"<<std::endl;
}

void CCProtonPi0_Plotter::plotStandardHistograms(rootDir &dir, std::string plotDir)
{
    std::cout<<">> Plotting Standard Histograms"<<std::endl;

    //DrawDataMC(dir, "E", plotDir);
    //DrawDataMC(dir, "P", plotDir);
    //DrawDataMC(dir, "KE", plotDir);
    //DrawDataMC(dir, "theta", plotDir);
    //DrawDataMC(dir, "phi", plotDir);

    //DrawDataStackedMC(dir, "E", plotDir);
    DrawDataStackedMC(dir, "P", plotDir);
    //DrawDataStackedMC(dir, "KE", plotDir);
    DrawDataStackedMC(dir, "theta", plotDir);
    DrawDataStackedMC(dir, "phi", plotDir);
}

void CCProtonPi0_Plotter::plotCutHistograms_MCOnly()
{
    std::cout<<"Plotting Cut Histograms MC Only"<<std::endl;
    std::string plotDir = plotDir_CutHists;

   // DrawStackedMC_BckgAll(rootDir_CutHists,"hCut_nTracks",plotDir);
   // DrawStackedMC_BckgAll(rootDir_CutHists,"hCut_nTracks2",plotDir);
   // DrawStackedMC_BckgAll(rootDir_CutHists,"hCut_nTracks_Close",plotDir);
   // DrawStackedMC_BckgAll(rootDir_CutHists,"hCut_nTracks_Far",plotDir);
   // DrawStackedMC_BckgAll(rootDir_CutHists,"hCut_nTracks_Discarded",plotDir);

    // ------------------------------------------------------------------------
    // Common
    // ------------------------------------------------------------------------
    CutArrow Vertex_Count(3,"L"); 
    DrawStackedMC_BckgAll(rootDir_CutHists,"hCut_nVertices",plotDir, 1, Vertex_Count);
       
    CutArrow Proton_Count(3,"L"); 
    DrawStackedMC_BckgAll(rootDir_CutHists,"hCut_nProtonCandidates",plotDir,1 , Proton_Count);

    CutArrow Michel(1,"L"); 
    DrawStackedMC_BckgAll(rootDir_CutHists,"hCut_Michel",plotDir,1, Michel);

    // ------------------------------------------------------------------------
    // 1 Track
    // ------------------------------------------------------------------------
    CutArrow eVis_nuclearTarget_1Track(20,"L"); 
    DrawStackedMC_BckgAll(rootDir_CutHists,"hCut_1Track_eVis_nuclearTarget",plotDir, 1, eVis_nuclearTarget_1Track);

    CutArrow eVis_other_min_1Track(50,"R"); 
    CutArrow eVis_other_max_1Track(2000,"L"); 
    DrawStackedMC_BckgAll(rootDir_CutHists,"hCut_1Track_eVis_other",plotDir, 2, eVis_other_min_1Track, eVis_other_max_1Track);

    CutArrow gamma1_ConvDist_1Track(15,"R"); 
    DrawStackedMC_BckgAll(rootDir_CutHists,"hCut_1Track_gamma1ConvDist",plotDir, 1, gamma1_ConvDist_1Track);

    CutArrow gamma2_ConvDist_1Track(15,"R"); 
    DrawStackedMC_BckgAll(rootDir_CutHists,"hCut_1Track_gamma2ConvDist",plotDir, 1, gamma2_ConvDist_1Track);

    CutArrow pi0invMass_min_1Track(75,"R"); 
    CutArrow pi0invMass_max_1Track(195,"L"); 
    DrawStackedMC_BckgAll(rootDir_CutHists,"hCut_1Track_pi0invMass",plotDir, 2, pi0invMass_min_1Track, pi0invMass_max_1Track);
    DrawStackedMC_BckgAll(rootDir_CutHists,"hCut_1Track_pi0invMass",plotDir, 2, pi0invMass_min_1Track, pi0invMass_max_1Track);

    CutArrow neutrinoE_1Track(20,"L"); 
    DrawStackedMC_BckgAll(rootDir_CutHists,"hCut_1Track_neutrinoE",plotDir, 1, neutrinoE_1Track);

    // ------------------------------------------------------------------------
    // 2 Track
    // ------------------------------------------------------------------------
    CutArrow eVis_nuclearTarget_2Track(20,"L"); 
    DrawStackedMC_BckgAll(rootDir_CutHists,"hCut_2Track_eVis_nuclearTarget",plotDir, 1, eVis_nuclearTarget_2Track);

    CutArrow eVis_other_min_2Track(50,"R"); 
    CutArrow eVis_other_max_2Track(2000,"L"); 
    DrawStackedMC_BckgAll(rootDir_CutHists,"hCut_2Track_eVis_other",plotDir, 2, eVis_other_min_2Track, eVis_other_max_2Track);

    CutArrow gamma1_ConvDist_2Track(15,"R"); 
    DrawStackedMC_BckgAll(rootDir_CutHists,"hCut_2Track_gamma1ConvDist",plotDir, 1, gamma1_ConvDist_2Track);

    CutArrow gamma2_ConvDist_2Track(15,"R"); 
    DrawStackedMC_BckgAll(rootDir_CutHists,"hCut_2Track_gamma2ConvDist",plotDir, 1, gamma2_ConvDist_2Track);

    CutArrow pi0invMass_min_2Track(75,"R"); 
    CutArrow pi0invMass_max_2Track(195,"L"); 
    DrawStackedMC_BckgAll(rootDir_CutHists,"hCut_2Track_pi0invMass",plotDir, 2, pi0invMass_min_2Track, pi0invMass_max_2Track);
    DrawStackedMC_BckgAll(rootDir_CutHists,"hCut_2Track_pi0invMass",plotDir, 2, pi0invMass_min_2Track, pi0invMass_max_2Track);

    CutArrow neutrinoE_2Track(20,"L"); 
    DrawStackedMC_BckgAll(rootDir_CutHists,"hCut_2Track_neutrinoE",plotDir, 1, neutrinoE_2Track);

    CutArrow protonScore_pIDDiff(0.45,"R"); 
    DrawStackedMC_BckgAll(rootDir_CutHists,"hCut_2Track_protonScore_pIDDiff",plotDir, 1, protonScore_pIDDiff);

    CutArrow protonScore_LLR(-10,"R"); 
    DrawStackedMC_BckgAll(rootDir_CutHists,"hCut_2Track_protonScore_LLR",plotDir, 1, protonScore_LLR);

    DrawStackedMC_BckgAll(rootDir_CutHists,"hCut_2Track_deltaInvMass",plotDir);

    std::cout<<"Plotting CutHistograms MC Only Finished!"<<std::endl;
}


void CCProtonPi0_Plotter::plotCutHistograms_DataMC()
{
    std::cout<<"Plotting CutHistograms Data vs MC"<<std::endl;
    std::string plotDir = plotDir_CutHists;

    // ------------------------------------------------------------------------
    // Common
    // ------------------------------------------------------------------------
    DrawDataStackedMC(rootDir_CutHists,"hCut_nVertices",plotDir);
    DrawDataStackedMC(rootDir_CutHists,"hCut_nTracks",plotDir);
    DrawDataStackedMC(rootDir_CutHists,"hCut_nTracks2",plotDir);
    DrawDataStackedMC(rootDir_CutHists,"hCut_nTracks_Close",plotDir);
    DrawDataStackedMC(rootDir_CutHists,"hCut_nTracks_Far",plotDir);
    DrawDataStackedMC(rootDir_CutHists,"hCut_nTracks_Discarded",plotDir);

    CutArrow Michel_1Track(1,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_Michel",plotDir,1, Michel_1Track);

    // ------------------------------------------------------------------------
    // 1 Track
    // ------------------------------------------------------------------------
    CutArrow eVis_nuclearTarget_1Track(20,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_eVis_nuclearTarget",plotDir, 1, eVis_nuclearTarget_1Track);

    CutArrow eVis_other_min_1Track(50,"R"); 
    CutArrow eVis_other_max_1Track(2000,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_eVis_other",plotDir, 2, eVis_other_min_1Track, eVis_other_max_1Track);

    CutArrow gamma1_ConvDist_1Track(15,"R"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_gamma1ConvDist",plotDir, 1, gamma1_ConvDist_1Track);

    CutArrow gamma2_ConvDist_1Track(15,"R"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_gamma2ConvDist",plotDir, 1, gamma2_ConvDist_1Track);

    CutArrow pi0invMass_min_1Track(75,"R"); 
    CutArrow pi0invMass_max_1Track(195,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_pi0invMass",plotDir, 2, pi0invMass_min_1Track, pi0invMass_max_1Track);
    DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_pi0invMass",plotDir, 2, pi0invMass_min_1Track, pi0invMass_max_1Track);

    CutArrow neutrinoE_1Track(20,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_neutrinoE",plotDir, 1, neutrinoE_1Track);

    // ------------------------------------------------------------------------
    // 2 Track
    // ------------------------------------------------------------------------
    CutArrow eVis_nuclearTarget_2Track(20,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_eVis_nuclearTarget",plotDir, 1, eVis_nuclearTarget_2Track);

    CutArrow eVis_other_min_2Track(50,"R"); 
    CutArrow eVis_other_max_2Track(2000,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_eVis_other",plotDir, 2, eVis_other_min_2Track, eVis_other_max_2Track);

    CutArrow gamma1_ConvDist_2Track(15,"R"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_gamma1ConvDist",plotDir, 1, gamma1_ConvDist_2Track);

    CutArrow gamma2_ConvDist_2Track(15,"R"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_gamma2ConvDist",plotDir, 1, gamma2_ConvDist_2Track);

    CutArrow pi0invMass_min_2Track(75,"R"); 
    CutArrow pi0invMass_max_2Track(195,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_pi0invMass",plotDir, 2, pi0invMass_min_2Track, pi0invMass_max_2Track);
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_pi0invMass",plotDir, 2, pi0invMass_min_2Track, pi0invMass_max_2Track);

    CutArrow neutrinoE_2Track(20,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_neutrinoE",plotDir, 1, neutrinoE_2Track);

    CutArrow protonScore_pIDDiff(0.45,"R"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_protonScore_pIDDiff",plotDir, 1, protonScore_pIDDiff);

    CutArrow protonScore_LLR(10,"R"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_protonScore_LLR",plotDir, 1, protonScore_LLR);

    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_deltaInvMass",plotDir);

    std::cout<<"Plotting CutHistograms Data vs MC Finished!"<<std::endl;
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
    rootDir_CutHists.mc = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + "CutHistograms_v2_48a.root";
    rootDir_Interaction.mc = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + ana_folder + "Interaction.root";
    rootDir_Muon.mc = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + ana_folder + "Muon.root";
    rootDir_Proton.mc = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + ana_folder + "Proton.root";
    rootDir_Pion.mc = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + ana_folder + "Pion.root";
    rootDir_Pi0Blob.mc = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + ana_folder + "Pi0Blob.root";

    // Set Data Root Dir
    rootDir_CutHists.data = Folder_List::rootOut + Folder_List::Data + Folder_List::analyzed + "CutHistograms_v2_45a.root";
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
    plotDir_Pi0Blob = Folder_List::output + Folder_List::plotOut + ana_folder + "Pi0Blob/";
    plotDir_Other = Folder_List::output + Folder_List::plotOut + ana_folder + "Other/";
}

void CCProtonPi0_Plotter::SavePi0InvMassPoints()
{
    std::cout<<"Saving Pi0 Invariant Mass Points"<<std::endl;

    std::string rootDir[2];
    std::string textDir[2];
    std::ofstream text[2];

    // Set Dirs
    textDir[0] = Folder_List::output + Folder_List::textOut + "Pi0InvMass_MC_minerva1_v2_38.txt";
    rootDir[0] = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + "CutHistograms.root";

    textDir[1] = Folder_List::output + Folder_List::textOut + "Pi0InvMass_Data_minerva1_v2_38.txt";
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

void CCProtonPi0_Plotter::plotEfficiencyCurves()
{
    // Efficiency Curve for Momentum
    std::string root_dir_all = rootDir_CutHists.mc;
    std::string root_dir_pi0 = rootDir_Pion.mc;
    
    TFile* f_all = new TFile(root_dir_all.c_str());
    TFile* f_pi0 = new TFile(root_dir_pi0.c_str());
    
    TH1D* all_signal_P = (TH1D*)f_all->Get("all_signal_pi0_P");
    TH1D* signal_P = (TH1D*)f_pi0->Get("signal_P");

    DrawEfficiencyCurve("Efficiency_P",plotDir_Pion,all_signal_P,signal_P);
}



#endif

