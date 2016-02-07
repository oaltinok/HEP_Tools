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
    plotInteraction_MCOnly();
    plotMuon_MCOnly();
    //plotProton_MCOnly();
    plotPion_MCOnly();
    //plotCutHistograms_MCOnly();
    //plotPi0Blob_MCOnly();
    plotEfficiencyCurves();

    //--------------------------------------------------------------------------
    //  Plot Function Reserved for Other Studies
    //--------------------------------------------------------------------------
    //SavePi0InvMassPoints();
    //plotOtherStudies();
}

void CCProtonPi0_Plotter::getPOT_MC()
{
    std::string playlist = "Input/Playlists/pl_MC_Merged.dat"; 
    POTCounter pot_counter;
    double total_pot = pot_counter.getPOTfromPlaylist(playlist);

    std::cout<<"Total MC POT = "<<total_pot<<std::endl;
}

void CCProtonPi0_Plotter::getPOT_Data()
{
    std::string playlist = "Input/Playlists/pl_Data_Merged.dat"; 
    POTCounter pot_counter;
    double total_pot = pot_counter.getPOTfromPlaylist(playlist);

    std::cout<<"Total Data POT = "<<total_pot<<std::endl;
}

CCProtonPi0_Plotter::CCProtonPi0_Plotter() : CCProtonPi0_NTupleAnalysis()
{
    //--------------------------------------------------------------------------
    // Set POT -- Run getPOT_MC() and getPOT_Data() Functions once to get POT
    //--------------------------------------------------------------------------
    data_POT = 3.33464e+20; 
    mc_POT = 2.73784e+21;
    POT_Ratio_data_mc = data_POT / mc_POT;
    std::cout<<"POT Data = "<<data_POT<<std::endl;
    std::cout<<"POT MC = "<<mc_POT<<std::endl;
    std::cout<<"POT Ratio = "<<POT_Ratio_data_mc<<std::endl;

    setRootDirs(); 
}

void CCProtonPi0_Plotter::plotPi0Blob_MCOnly()
{
    std::cout<<"Plotting Pi0Blob MC Only"<<std::endl;
    std::string plotDir = Folder_List::plotDir_Pi0Blob; // Plots go under Pion

    // Evis Fraction 
    DrawStackedMC(rootDir_Pi0Blob,"evis_frac_reco_pi0_true_pi0",plotDir);
    DrawStackedMC(rootDir_Pi0Blob,"evis_frac_true_pi0_reco_all",plotDir);
    DrawStackedMC(rootDir_Pi0Blob,"evis_frac_reco_pi0_reco_all",plotDir);
    DrawStackedMC(rootDir_Pi0Blob,"evis_frac_reco_nonpi0_reco_all",plotDir);

    // Evis By Particle
    DrawStackedMC_GammaByPDG(rootDir_Pi0Blob,"evis",1,plotDir);
    DrawStackedMC_GammaByPDG(rootDir_Pi0Blob,"evis",2,plotDir);
    DrawStackedMC_GammaByPDG(rootDir_Pi0Blob,"evis",3,plotDir);

    //    // Plot Truth Match - gamma 1
    //    DrawStackedMC(rootDir_Pi0Blob,"g1_evis_most_pdg",plotDir);
    //    DrawStackedMC(rootDir_Pi0Blob,"g1_evis_total_truth",plotDir);
    //    DrawStackedMC(rootDir_Pi0Blob,"g1_evis_frac_pizero",plotDir);
    //    DrawStackedMC(rootDir_Pi0Blob,"g1_evis_frac_piplus",plotDir);
    //    DrawStackedMC(rootDir_Pi0Blob,"g1_evis_frac_piminus",plotDir);
    //    DrawStackedMC(rootDir_Pi0Blob,"g1_evis_frac_proton",plotDir);
    //    DrawStackedMC(rootDir_Pi0Blob,"g1_evis_frac_neutron",plotDir);
    //    DrawStackedMC(rootDir_Pi0Blob,"g1_evis_frac_muon",plotDir);
    //
    //    // Plot Truth Match - gamma 2
    //    DrawStackedMC(rootDir_Pi0Blob,"g2_evis_most_pdg",plotDir);
    //    DrawStackedMC(rootDir_Pi0Blob,"g2_evis_total_truth",plotDir);
    //    DrawStackedMC(rootDir_Pi0Blob,"g2_evis_frac_pizero",plotDir);
    //    DrawStackedMC(rootDir_Pi0Blob,"g2_evis_frac_piplus",plotDir);
    //    DrawStackedMC(rootDir_Pi0Blob,"g2_evis_frac_piminus",plotDir);
    //    DrawStackedMC(rootDir_Pi0Blob,"g2_evis_frac_proton",plotDir);
    //    DrawStackedMC(rootDir_Pi0Blob,"g2_evis_frac_neutron",plotDir);
    //    DrawStackedMC(rootDir_Pi0Blob,"g2_evis_frac_muon",plotDir);
    //    
    //    // Total Pi0 Evis
    //    DrawStackedMC(rootDir_Pi0Blob,"captured_evis_frac_all",plotDir);
    //    DrawStackedMC(rootDir_Pi0Blob,"captured_evis_frac_signal",plotDir);

    //    // Blob 1 Cluster Z and Strip Numbers
    //    DrawStackedMC(rootDir_Pi0Blob,"g1_cluster_Z",plotDir);
    //    DrawStackedMC(rootDir_Pi0Blob,"g1_max_cluster_Z",plotDir);
    //    DrawStackedMC(rootDir_Pi0Blob,"g1_min_cluster_Z",plotDir);
    //    DrawStackedMC(rootDir_Pi0Blob,"g1_strips",plotDir);
    //    DrawStackedMC(rootDir_Pi0Blob,"g1_max_strip",plotDir);
    //    DrawStackedMC(rootDir_Pi0Blob,"g1_min_strip",plotDir);
    //    DrawStackedMC(rootDir_Pi0Blob,"g1_nPlanes",plotDir);
    //    // Blob 2 Cluster Z and Strip Numbers
    //    DrawStackedMC(rootDir_Pi0Blob,"g2_cluster_Z",plotDir);
    //    DrawStackedMC(rootDir_Pi0Blob,"g2_max_cluster_Z",plotDir);
    //    DrawStackedMC(rootDir_Pi0Blob,"g2_min_cluster_Z",plotDir);
    //    DrawStackedMC(rootDir_Pi0Blob,"g2_strips",plotDir);
    //    DrawStackedMC(rootDir_Pi0Blob,"g2_max_strip",plotDir);
    //    DrawStackedMC(rootDir_Pi0Blob,"g2_min_strip",plotDir);
    //    DrawStackedMC(rootDir_Pi0Blob,"g2_nPlanes",plotDir);

    std::cout<<"Plotting Pi0Blob MC Only Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plotPi0Blob_DataMC()
{
    // Do Nothing 
}

void CCProtonPi0_Plotter::plotOtherStudies()
{
    std::cout<<"Plotting Other Studies..."<<std::endl;
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    Draw2DHist(rootDir_CutHists,"signal_gamma_E_cos_openingAngle",plotDir);
    Draw2DHist(rootDir_CutHists,"bckg_gamma_E_cos_openingAngle",plotDir);
    Draw2DHist(rootDir_CutHists,"bckg_signal_diff_E_cos_openingAngle",plotDir);

    std::cout<<"Plotting Other Studies Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plotInteraction_MCOnly()
{
    std::cout<<"Plotting Interaction MC Only"<<std::endl;
    std::string plotDir = Folder_List::plotDir_Interaction;

    //plot_mc_w_Stacked();
    plot_final_mc_w_Stacked();

    // Response
    DrawNormalizedMigrationHistogram(rootDir_Interaction,"response_neutrino_E",plotDir);
    DrawNormalizedMigrationHistogram(rootDir_Interaction,"response_QSq",plotDir);

    //DrawStackedMC(rootDir_Interaction,"vertex_energy_1Track",plotDir);
    //DrawStackedMC(rootDir_Interaction,"vertex_energy_2Track",plotDir);
    //DrawStackedMC(rootDir_Interaction,"vertex_evis_1Track",plotDir);
    //DrawStackedMC(rootDir_Interaction,"vertex_evis_2Track",plotDir);
    //DrawStackedMC(rootDir_Interaction,"extra_evis_1Track",plotDir);
    //DrawStackedMC(rootDir_Interaction,"extra_evis_2Track",plotDir);

    //Draw1DHist(rootDir_Interaction,"n_ejected_nucleons_1Track",plotDir);
    //Draw1DHist(rootDir_Interaction,"n_ejected_nucleons_2Track",plotDir);

    //Draw1DHist(rootDir_Interaction,"proton_true_P_1Track",plotDir);
    //Draw1DHist(rootDir_Interaction,"proton_true_KE_1Track",plotDir);

    //Draw1DHist(rootDir_Interaction,"Enu_True_1Track",plotDir);
    //Draw1DHist(rootDir_Interaction,"Enu_True_2Track",plotDir);
    //Draw1DHist(rootDir_Interaction,"Enu_1Track_Error",plotDir);
    //Draw1DHist(rootDir_Interaction,"Enu_2Track_Error",plotDir);
    //Draw1DHist(rootDir_Interaction,"Enu_1Track_Alt_Error",plotDir);
    //
    //Draw1DHist(rootDir_Interaction,"Enu_1Track_Diff",plotDir);
    //Draw1DHist(rootDir_Interaction,"Enu_2Track_Diff",plotDir);

    //DrawStackedMC(rootDir_Interaction,"recovered_Pi0_P",plotDir);
    //DrawStackedMC(rootDir_Interaction,"recovered_Pi0_theta",plotDir);

    //Draw1DHist(rootDir_Interaction,"h_extra_muon_energy",plotDir);
    //Draw1DHist(rootDir_Interaction,"h_extra_dispersed_energy",plotDir);
    //Draw1DHist(rootDir_Interaction,"h_extra_rejected_energy",plotDir);

    //DrawStackedMC(rootDir_Interaction,"vertex_energy_1Track",plotDir);
    //DrawStackedMC(rootDir_Interaction,"vertex_evis_1Track",plotDir);
    //DrawStackedMC(rootDir_Interaction,"extra_evis_1Track",plotDir);
    //DrawStackedMC(rootDir_Interaction,"extra_muon_energy_1Track",plotDir);
    //DrawStackedMC(rootDir_Interaction,"extra_dispersed_energy_1Track",plotDir);
    //DrawStackedMC(rootDir_Interaction,"extra_rejected_energy_1Track",plotDir);
    //DrawStackedMC(rootDir_Interaction,"extra_total_energy_1Track",plotDir);
 
    //DrawStackedMC(rootDir_Interaction,"vertex_energy_2Track",plotDir);
    //DrawStackedMC(rootDir_Interaction,"vertex_evis_2Track",plotDir);
    //DrawStackedMC(rootDir_Interaction,"extra_evis_2Track",plotDir);
    //DrawStackedMC(rootDir_Interaction,"extra_muon_energy_2Track",plotDir);
    //DrawStackedMC(rootDir_Interaction,"extra_dispersed_energy_2Track",plotDir);
    //DrawStackedMC(rootDir_Interaction,"extra_rejected_energy_2Track",plotDir);
    //DrawStackedMC(rootDir_Interaction,"extra_total_energy_2Track",plotDir);
 
    //plot_stacked_pi0_P();
    //plot_stacked_pi0_theta();

    std::cout<<"Plotting Interaction MC Only Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plotInteraction_DataMC()
{
    std::cout<<"Plotting Interaction Data vs MC"<<std::endl;
    std::string plotDir = Folder_List::plotDir_Interaction;

    //DrawDataMC(rootDir_Interaction,"Enu_1Track",plotDir);
    //DrawDataMC(rootDir_Interaction,"Enu_1Track_Alt",plotDir);
    //DrawDataMC(rootDir_Interaction,"Enu_2Track",plotDir);
    //DrawDataMC(rootDir_Interaction,"Enu",plotDir);
    //DrawDataMC(rootDir_Interaction,"QSq",plotDir);
    //DrawDataMC(rootDir_Interaction,"WSq",plotDir);
    //DrawDataMC(rootDir_Interaction,"W",plotDir);
    //DrawDataMC(rootDir_Interaction,"deltaInvMass",plotDir);

    DrawDataStackedMC(rootDir_Interaction,"Enu_1Track",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"Enu_1Track_Alt",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"Enu_2Track",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"Enu",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"QSq",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"WSq",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"W",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"W_Calc",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"deltaInvMass",plotDir);

    DrawDataStackedMC(rootDir_Interaction,"vertex_energy_1Track",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"vertex_evis_1Track",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"extra_evis_1Track",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"extra_muon_energy_1Track",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"extra_dispersed_energy_1Track",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"extra_rejected_energy_1Track",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"extra_total_energy_1Track",plotDir);
 
    DrawDataStackedMC(rootDir_Interaction,"vertex_energy_2Track",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"vertex_evis_2Track",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"extra_evis_2Track",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"extra_muon_energy_2Track",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"extra_dispersed_energy_2Track",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"extra_rejected_energy_2Track",plotDir);
    DrawDataStackedMC(rootDir_Interaction,"extra_total_energy_2Track",plotDir);
       
    std::cout<<"Plotting Interaction Data vs MC Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plotMuon_MCOnly()
{
    std::cout<<"Plotting Muon MC Only"<<std::endl;
    std::string plotDir = Folder_List::plotDir_Muon;

    Draw1DHist(rootDir_Muon,"P_error",plotDir);
    Draw1DHist(rootDir_Muon,"E_error",plotDir);
    Draw1DHist(rootDir_Muon,"E_Diff",plotDir);
    
    DrawNormalizedMigrationHistogram(rootDir_Muon,"response_P",plotDir);
    DrawNormalizedMigrationHistogram(rootDir_Muon,"response_theta",plotDir);

    std::cout<<"Plotting Muon MC Only Finished!\n"<<std::endl;

}

void CCProtonPi0_Plotter::plotMuon_DataMC()
{
    std::cout<<"Plotting Muon Data vs MC"<<std::endl;
    std::string plotDir = Folder_List::plotDir_Muon;

    plotStandardHistograms(rootDir_Muon, plotDir);

    std::cout<<"Plotting Muon Data vs MC Finished!\n"<<std::endl;
}

void CCProtonPi0_Plotter::plotProton_MCOnly()
{
    std::cout<<"Plotting Proton MC Only"<<std::endl;
    std::string plotDir = Folder_List::plotDir_Proton;

    DrawStackedMC(rootDir_Proton, "P", plotDir);
    DrawStackedMC(rootDir_Proton, "theta", plotDir);
    DrawStackedMC(rootDir_Proton, "phi", plotDir);
    
    Draw1DHist(rootDir_Proton,"P_error",plotDir);
    Draw1DHist(rootDir_Proton,"E_error",plotDir);
    Draw1DHist(rootDir_Proton,"E_Diff",plotDir);

    std::cout<<"Plotting Proton MC Only Finished!\n"<<std::endl;
}


void CCProtonPi0_Plotter::plotProton_DataMC()
{    
    std::cout<<"Plotting Proton Data vs MC"<<std::endl;
    std::string plotDir = Folder_List::plotDir_Proton;

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
    std::string plotDir = Folder_List::plotDir_Pion;

    plotPion_True();

    DrawStackedMC(rootDir_Pion,"gamma1_E",plotDir);
    DrawStackedMC(rootDir_Pion,"gamma1_theta",plotDir);
    DrawStackedMC(rootDir_Pion,"gamma1_ConvLength",plotDir);

    DrawStackedMC(rootDir_Pion,"gamma2_E",plotDir);
    DrawStackedMC(rootDir_Pion,"gamma2_theta",plotDir);
    DrawStackedMC(rootDir_Pion,"gamma2_ConvLength",plotDir);

    DrawStackedMC(rootDir_Pion,"photonEnergy_Asymmetry",plotDir);
    DrawStackedMC(rootDir_Pion,"invMass",plotDir);
    DrawStackedMC(rootDir_Pion,"cos_openingAngle",plotDir);
    DrawStackedMC(rootDir_Pion,"P",plotDir);
    DrawStackedMC(rootDir_Pion,"E",plotDir);
    DrawStackedMC(rootDir_Pion,"theta",plotDir);
    DrawStackedMC(rootDir_Pion,"phi",plotDir);

    std::cout<<"Plotting Pion MC Only Finished!\n"<<std::endl;
}


void CCProtonPi0_Plotter::plotPion_DataMC()
{
    std::cout<<"Plotting Pion Data vs MC"<<std::endl;
    std::string plotDir = Folder_List::plotDir_Pion;

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
    std::string plotDir = Folder_List::plotDir_Pion;

    Draw1DHist(rootDir_Pion,"gamma1_true_E",plotDir);
    Draw1DHist(rootDir_Pion,"gamma1_reco_error_E",plotDir);
    Draw2DHist(rootDir_Pion,"gamma1_true_E_reco_E_error",plotDir);
    Draw2DHist(rootDir_Pion,"gamma1_reco_E_true_E",plotDir);

    Draw1DHist(rootDir_Pion,"gamma2_true_E",plotDir);
    Draw1DHist(rootDir_Pion,"gamma2_reco_error_E",plotDir);
    Draw2DHist(rootDir_Pion,"gamma2_true_E_reco_E_error",plotDir);
    Draw2DHist(rootDir_Pion,"gamma2_reco_E_true_E",plotDir);
    
    Draw2DHist(rootDir_Pion,"signal_gamma1_E_gamma2_E",plotDir);
    Draw2DHist(rootDir_Pion,"bckg_gamma1_E_gamma2_E",plotDir);
    Draw2DHist(rootDir_Pion,"bckg_signal_diff_E",plotDir);

    Draw2DHist(rootDir_Pion,"signal_gamma1_convLength_gamma2_convLength",plotDir);
    Draw2DHist(rootDir_Pion,"bckg_gamma1_convLength_gamma2_convLength",plotDir);
    Draw2DHist(rootDir_Pion,"bckg_signal_diff_convLength",plotDir);

    Draw1DHist(rootDir_Pion,"P_error",plotDir);
    Draw1DHist(rootDir_Pion,"E_error",plotDir);
    Draw1DHist(rootDir_Pion,"E_Diff",plotDir);
    
    DrawNormalizedMigrationHistogram(rootDir_Pion,"response_P",plotDir);
    DrawNormalizedMigrationHistogram(rootDir_Pion,"response_theta",plotDir);

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
    std::string plotDir = Folder_List::plotDir_CutHists;

   // DrawStackedMC(rootDir_CutHists,"hCut_nTracks",plotDir);
   // DrawStackedMC(rootDir_CutHists,"hCut_nTracks2",plotDir);
   // DrawStackedMC(rootDir_CutHists,"hCut_nTracks_Close",plotDir);
   // DrawStackedMC(rootDir_CutHists,"hCut_nTracks_Far",plotDir);
   // DrawStackedMC(rootDir_CutHists,"hCut_nTracks_Discarded",plotDir);

    // ------------------------------------------------------------------------
    // Common
    // ------------------------------------------------------------------------
    CutArrow Vertex_Count(3,"L"); 
    DrawStackedMC(rootDir_CutHists,"hCut_nVertices",plotDir, 1, Vertex_Count);
       
    CutArrow Proton_Count(3,"L"); 
    DrawStackedMC(rootDir_CutHists,"hCut_nProtonCandidates",plotDir,1 , Proton_Count);
    
    DrawStackedMC(rootDir_CutHists,"hCut_nShowerCandidates",plotDir);
    //DrawStackedMC(rootDir_CutHists,"hCut_1Track_nShowerCandidates",plotDir);
    //DrawStackedMC(rootDir_CutHists,"hCut_2Track_nShowerCandidates",plotDir);

    //CutArrow Michel(1,"L"); 
    //DrawStackedMC(rootDir_CutHists,"hCut_Michel",plotDir,1, Michel);

    CutArrow pi0invMass_min(60,"R"); 
    CutArrow pi0invMass_max(200,"L"); 
    DrawSignalMC(rootDir_CutHists,"hCut_pi0invMass",plotDir, 2, pi0invMass_min, pi0invMass_max);
    DrawSignalMC(rootDir_CutHists,"hCut_pi0invMass_Old",plotDir, 2, pi0invMass_min, pi0invMass_max);
    DrawStackedMC(rootDir_CutHists,"hCut_pi0invMass",plotDir, 2, pi0invMass_min, pi0invMass_max);
    DrawStackedMC(rootDir_CutHists,"hCut_pi0invMass_Old",plotDir, 2, pi0invMass_min, pi0invMass_max);

    // ------------------------------------------------------------------------
    // 1 Track
    // ------------------------------------------------------------------------
    CutArrow eVis_nuclearTarget_1Track(20,"L"); 
    DrawStackedMC(rootDir_CutHists,"hCut_1Track_eVis_nuclearTarget",plotDir, 1, eVis_nuclearTarget_1Track);

    CutArrow eVis_other_min_1Track(50,"R"); 
    CutArrow eVis_other_max_1Track(2500,"L"); 
    DrawStackedMC(rootDir_CutHists,"hCut_1Track_eVis_other",plotDir, 2, eVis_other_min_1Track, eVis_other_max_1Track);

    CutArrow gamma1_ConvDist_1Track(14,"R"); 
    DrawStackedMC(rootDir_CutHists,"hCut_1Track_gamma1ConvDist",plotDir, 1, gamma1_ConvDist_1Track);

    //CutArrow gamma2_ConvDist_1Track(15,"R"); 
    //DrawStackedMC(rootDir_CutHists,"hCut_1Track_gamma2ConvDist",plotDir, 1, gamma2_ConvDist_1Track);
    DrawStackedMC(rootDir_CutHists,"hCut_1Track_gamma2ConvDist",plotDir);

    DrawSignalMC(rootDir_CutHists,"hCut_1Track_pi0invMass",plotDir, 2, pi0invMass_min, pi0invMass_max);
    DrawSignalMC(rootDir_CutHists,"hCut_1Track_pi0invMass_Old",plotDir, 2, pi0invMass_min, pi0invMass_max);
    DrawStackedMC(rootDir_CutHists,"hCut_1Track_pi0invMass",plotDir, 2, pi0invMass_min, pi0invMass_max);
    DrawStackedMC(rootDir_CutHists,"hCut_1Track_pi0invMass_Old",plotDir, 2, pi0invMass_min, pi0invMass_max);

    CutArrow neutrinoE_1Track(20,"L"); 
    DrawStackedMC(rootDir_CutHists,"hCut_1Track_neutrinoE",plotDir, 1, neutrinoE_1Track);

    // ------------------------------------------------------------------------
    // 2 Track
    // ------------------------------------------------------------------------
    CutArrow eVis_nuclearTarget_2Track(20,"L"); 
    DrawStackedMC(rootDir_CutHists,"hCut_2Track_eVis_nuclearTarget",plotDir, 1, eVis_nuclearTarget_2Track);

    CutArrow eVis_other_min_2Track(50,"R"); 
    CutArrow eVis_other_max_2Track(2000,"L"); 
    DrawStackedMC(rootDir_CutHists,"hCut_2Track_eVis_other",plotDir, 2, eVis_other_min_2Track, eVis_other_max_2Track);

    CutArrow gamma1_ConvDist_2Track(14,"R"); 
    DrawStackedMC(rootDir_CutHists,"hCut_2Track_gamma1ConvDist",plotDir, 1, gamma1_ConvDist_2Track);

    //CutArrow gamma2_ConvDist_2Track(15,"R"); 
    //DrawStackedMC(rootDir_CutHists,"hCut_2Track_gamma2ConvDist",plotDir, 1, gamma2_ConvDist_2Track);
    DrawStackedMC(rootDir_CutHists,"hCut_2Track_gamma2ConvDist",plotDir);

    DrawSignalMC(rootDir_CutHists,"hCut_2Track_pi0invMass",plotDir, 2, pi0invMass_min, pi0invMass_max);
    DrawSignalMC(rootDir_CutHists,"hCut_2Track_pi0invMass_Old",plotDir, 2, pi0invMass_min, pi0invMass_max);
    DrawStackedMC(rootDir_CutHists,"hCut_2Track_pi0invMass",plotDir, 2, pi0invMass_min, pi0invMass_max);
    DrawStackedMC(rootDir_CutHists,"hCut_2Track_pi0invMass_Old",plotDir, 2, pi0invMass_min, pi0invMass_max);

    CutArrow neutrinoE_2Track(20,"L"); 
    DrawStackedMC(rootDir_CutHists,"hCut_2Track_neutrinoE",plotDir, 1, neutrinoE_2Track);

    //CutArrow protonScore_pIDDiff(0.45,"R"); 
    //DrawStackedMC(rootDir_CutHists,"hCut_2Track_protonScore_pIDDiff",plotDir, 1, protonScore_pIDDiff);

    CutArrow protonScore_LLR(-10,"R"); 
    DrawStackedMC(rootDir_CutHists,"hCut_2Track_protonScore_LLR",plotDir, 1, protonScore_LLR);

    DrawStackedMC(rootDir_CutHists,"hCut_2Track_deltaInvMass",plotDir);

    std::cout<<"Plotting CutHistograms MC Only Finished!"<<std::endl;
}


void CCProtonPi0_Plotter::plotCutHistograms_DataMC()
{
    std::cout<<"Plotting CutHistograms Data vs MC"<<std::endl;
    std::string plotDir = Folder_List::plotDir_CutHists;

   // DrawDataStackedMC(rootDir_CutHists,"hCut_nTracks",plotDir);
   // DrawDataStackedMC(rootDir_CutHists,"hCut_nTracks2",plotDir);
   // DrawDataStackedMC(rootDir_CutHists,"hCut_nTracks_Close",plotDir);
   // DrawDataStackedMC(rootDir_CutHists,"hCut_nTracks_Far",plotDir);
   // DrawDataStackedMC(rootDir_CutHists,"hCut_nTracks_Discarded",plotDir);

    // ------------------------------------------------------------------------
    // Common
    // ------------------------------------------------------------------------
    CutArrow Vertex_Count(3,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_nVertices",plotDir, 1, Vertex_Count);
       
    CutArrow Proton_Count(3,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_nProtonCandidates",plotDir,1 , Proton_Count);
    
    DrawDataStackedMC(rootDir_CutHists,"hCut_nShowerCandidates",plotDir);
    //DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_nShowerCandidates",plotDir);
    //DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_nShowerCandidates",plotDir);

    //CutArrow Michel(1,"L"); 
    //DrawDataStackedMC(rootDir_CutHists,"hCut_Michel",plotDir,1, Michel);

    CutArrow pi0invMass_min(60,"R"); 
    CutArrow pi0invMass_max(200,"L"); 
    DrawSignalMC(rootDir_CutHists,"hCut_pi0invMass",plotDir, 2, pi0invMass_min, pi0invMass_max);
    DrawSignalMC(rootDir_CutHists,"hCut_pi0invMass_Old",plotDir, 2, pi0invMass_min, pi0invMass_max);
    DrawDataStackedMC(rootDir_CutHists,"hCut_pi0invMass",plotDir, 2, pi0invMass_min, pi0invMass_max);
    DrawDataStackedMC(rootDir_CutHists,"hCut_pi0invMass_Old",plotDir, 2, pi0invMass_min, pi0invMass_max);

    // ------------------------------------------------------------------------
    // 1 Track
    // ------------------------------------------------------------------------
    CutArrow eVis_nuclearTarget_1Track(20,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_eVis_nuclearTarget",plotDir, 1, eVis_nuclearTarget_1Track);

    CutArrow eVis_other_min_1Track(50,"R"); 
    CutArrow eVis_other_max_1Track(2000,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_eVis_other",plotDir, 2, eVis_other_min_1Track, eVis_other_max_1Track);

    CutArrow gamma1_ConvDist_1Track(14,"R"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_gamma1ConvDist",plotDir, 1, gamma1_ConvDist_1Track);

    //CutArrow gamma2_ConvDist_1Track(15,"R"); 
    //DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_gamma2ConvDist",plotDir, 1, gamma2_ConvDist_1Track);
    DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_gamma2ConvDist",plotDir);

    DrawSignalMC(rootDir_CutHists,"hCut_1Track_pi0invMass",plotDir, 2, pi0invMass_min, pi0invMass_max);
    DrawSignalMC(rootDir_CutHists,"hCut_1Track_pi0invMass_Old",plotDir, 2, pi0invMass_min, pi0invMass_max);
    DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_pi0invMass",plotDir, 2, pi0invMass_min, pi0invMass_max);
    DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_pi0invMass_Old",plotDir, 2, pi0invMass_min, pi0invMass_max);

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

    CutArrow gamma1_ConvDist_2Track(14,"R"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_gamma1ConvDist",plotDir, 1, gamma1_ConvDist_2Track);

    //CutArrow gamma2_ConvDist_2Track(15,"R"); 
    //DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_gamma2ConvDist",plotDir, 1, gamma2_ConvDist_2Track);
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_gamma2ConvDist",plotDir);

    DrawSignalMC(rootDir_CutHists,"hCut_2Track_pi0invMass",plotDir, 2, pi0invMass_min, pi0invMass_max);
    DrawSignalMC(rootDir_CutHists,"hCut_2Track_pi0invMass_Old",plotDir, 2, pi0invMass_min, pi0invMass_max);
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_pi0invMass",plotDir, 2, pi0invMass_min, pi0invMass_max);
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_pi0invMass_Old",plotDir, 2, pi0invMass_min, pi0invMass_max);

    CutArrow neutrinoE_2Track(20,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_neutrinoE",plotDir, 1, neutrinoE_2Track);

    //CutArrow protonScore_pIDDiff(0.45,"R"); 
    //DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_protonScore_pIDDiff",plotDir, 1, protonScore_pIDDiff);

    CutArrow protonScore_LLR(-10,"R"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_protonScore_LLR",plotDir, 1, protonScore_LLR);

    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_deltaInvMass",plotDir);

    std::cout<<"Plotting CutHistograms Data vs MC Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plot_mc_w_Stacked()
{
    // mc_w written during reduce - Its Histogram is with Cut Hists
    std::string root_dir = rootDir_CutHists.mc;
    std::string plotDir = Folder_List::plotDir_Interaction;

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

    double nEvents = h_mc_w_DIS->GetEntries() + h_mc_w_RES->GetEntries();
    std::cout<<"nEvents = "<<nEvents<<std::endl;

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
    std::string plotDir = Folder_List::plotDir_Interaction;

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


    double nEvents = h_mc_w_DIS->GetEntries() + h_mc_w_RES->GetEntries();
    std::cout<<"nEvents = "<<nEvents<<std::endl;

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

void CCProtonPi0_Plotter::plot_stacked_pi0_P()
{
    std::string root_dir = rootDir_Interaction.mc;
    std::string plotDir = Folder_List::plotDir_Interaction;

    std::cout<<"\nPlottting Stacked Pi0 Momentum"<<std::endl;

    TFile* f_Root = new TFile(root_dir.c_str());
    TCanvas* c1 = new TCanvas();
    TLegend *legend = new TLegend(0.7,0.8,0.9,0.9);  

    TH1D* h_original = (TH1D*)f_Root->Get("h_original_Pi0_P");
    h_original->SetFillColor(kBlue);
    h_original->SetMarkerStyle(21);
    h_original->SetMarkerColor(kBlue);

    TH1D* h_recovered = (TH1D*)f_Root->Get("h_recovered_Pi0_P");
    h_recovered->SetFillColor(kGreen);
    h_recovered->SetMarkerStyle(21);
    h_recovered->SetMarkerColor(kGreen);

    legend->AddEntry(h_original, "Original", "f");
    legend->AddEntry(h_recovered, "Recovered", "f");

    h_original->Draw();
    h_recovered->Draw("same");

    legend->Draw();

    c1->Print(Form("%s%s",plotDir.c_str(),"compare_pi0_P.png"),"png");

    delete f_Root;
    delete legend;
    delete c1;
}

void CCProtonPi0_Plotter::plot_stacked_pi0_theta()
{
    std::string root_dir = rootDir_Interaction.mc;
    std::string plotDir = Folder_List::plotDir_Interaction;

    std::cout<<"\nPlottting Stacked Pi0 Theta"<<std::endl;

    TFile* f_Root = new TFile(root_dir.c_str());
    TCanvas* c1 = new TCanvas();
    TLegend *legend = new TLegend(0.7,0.8,0.9,0.9);  

    TH1D* h_original = (TH1D*)f_Root->Get("h_original_Pi0_theta");
    h_original->SetFillColor(kBlue);
    h_original->SetMarkerStyle(21);
    h_original->SetMarkerColor(kBlue);

    TH1D* h_recovered = (TH1D*)f_Root->Get("h_recovered_Pi0_theta");
    h_recovered->SetFillColor(kGreen);
    h_recovered->SetMarkerStyle(21);
    h_recovered->SetMarkerColor(kGreen);

    legend->AddEntry(h_original, "Original", "f");
    legend->AddEntry(h_recovered, "Recovered", "f");

    h_original->Draw();
    h_recovered->Draw("same");

    legend->Draw();

    c1->Print(Form("%s%s",plotDir.c_str(),"compare_pi0_theta.png"),"png");

    delete f_Root;
    delete legend;
    delete c1;
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
    std::string plotDir_Efficiency = Folder_List::plotDir_Efficiency;
    std::string root_dir_all = rootDir_Truth.mc;
    std::string root_dir_muon = rootDir_Muon.mc;
    std::string root_dir_pi0 = rootDir_Pion.mc;
    std::string root_dir_interaction = rootDir_Interaction.mc;
    
    TFile* f_all = new TFile(root_dir_all.c_str());
    TFile* f_muon = new TFile(root_dir_muon.c_str());
    TFile* f_pi0 = new TFile(root_dir_pi0.c_str());
    TFile* f_interaction = new TFile(root_dir_interaction.c_str());
    
    // Efficiency Curve for Muon Momentum
    TH1D* all_signal_muon_P = (TH1D*)f_all->Get("all_signal_muon_P");
    TH1D* signal_muon_P = (TH1D*)f_muon->Get("eff_P");
    DrawEfficiencyCurve("Efficiency_Muon_P",plotDir_Efficiency, all_signal_muon_P, signal_muon_P);
 
    // Efficiency Curve for Muon Theta 
    TH1D* all_signal_muon_theta = (TH1D*)f_all->Get("all_signal_muon_theta");
    TH1D* signal_muon_theta = (TH1D*)f_muon->Get("eff_theta");
    DrawEfficiencyCurve("Efficiency_Muon_theta",Folder_List::plotDir_Efficiency, all_signal_muon_theta, signal_muon_theta);

    // Efficiency Curve for Pi0 Momentum
    TH1D* all_signal_pi0_P = (TH1D*)f_all->Get("all_signal_pi0_P");
    TH1D* signal_pi0_P = (TH1D*)f_pi0->Get("eff_P");
    DrawEfficiencyCurve("Efficiency_Pi0_P",plotDir_Efficiency, all_signal_pi0_P, signal_pi0_P);

    // Efficiency Curve for Pi0 Theta 
    TH1D* all_signal_pi0_theta = (TH1D*)f_all->Get("all_signal_pi0_theta");
    TH1D* signal_pi0_theta = (TH1D*)f_pi0->Get("eff_theta");
    DrawEfficiencyCurve("Efficiency_Pi0_theta",plotDir_Efficiency, all_signal_pi0_theta, signal_pi0_theta);
 
    // Efficiency Curve for Neutrino Energy 
    TH1D* all_signal_neutrino_E = (TH1D*)f_all->Get("all_signal_neutrino_E");
    TH1D* signal_neutrino_E = (TH1D*)f_interaction->Get("eff_neutrino_E");
    DrawEfficiencyCurve("Efficiency_Neutrino_E",plotDir_Efficiency, all_signal_neutrino_E, signal_neutrino_E);
  
    // Efficiency Curve for QSq 
    TH1D* all_signal_QSq = (TH1D*)f_all->Get("all_signal_QSq");
    TH1D* signal_QSq = (TH1D*)f_interaction->Get("eff_QSq");
    DrawEfficiencyCurve("Efficiency_QSq",plotDir_Efficiency, all_signal_QSq, signal_QSq);
   
}

void CCProtonPi0_Plotter::setRootDirs()
{
    // Set Other Studies ROOT Dir
    rootDir_Truth.mc = Folder_List::rootDir_Truth_mc;
    rootDir_Truth.data = Folder_List::rootDir_Truth_data;

    rootDir_OtherStudies.mc = Folder_List::rootDir_OtherStudies_mc; 
    rootDir_OtherStudies.data = Folder_List::rootDir_OtherStudies_data;

    // Set MC Root Dir
    rootDir_CutHists.mc = Folder_List::rootDir_CutHists_mc;
    rootDir_Interaction.mc = Folder_List::rootDir_Interaction_mc;
    rootDir_Muon.mc = Folder_List::rootDir_Muon_mc;
    rootDir_Proton.mc = Folder_List::rootDir_Proton_mc;
    rootDir_Pion.mc = Folder_List::rootDir_Pion_mc;
    rootDir_Pi0Blob.mc = Folder_List::rootDir_Pi0Blob_mc;

    // Set Data Root Dir
    rootDir_CutHists.data = Folder_List::rootDir_CutHists_data;
    rootDir_Interaction.data = Folder_List::rootDir_Interaction_data;
    rootDir_Muon.data = Folder_List::rootDir_Muon_data;
    rootDir_Proton.data = Folder_List::rootDir_Proton_data;
    rootDir_Pion.data = Folder_List::rootDir_Pion_data;
    rootDir_Pi0Blob.data = Folder_List::rootDir_Pi0Blob_data;
}

#endif

