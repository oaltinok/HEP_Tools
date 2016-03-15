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
    // Cross Sections
    //--------------------------------------------------------------------------
    //plotErrorSummary();
    //plotOriginalData();
    //plotBackgroundEstimated();
    //plotBackgroundSubtracted();
    //plotUnfolded();
    //plotEfficiencyCorrected();
    plotIntegratedFlux();
    plotCrossSection();
    //plotCrossSection_Check();

    //--------------------------------------------------------------------------
    //  Data vs MC
    //--------------------------------------------------------------------------
    //plotInteraction_DataMC();
    //plotMuon_DataMC();
    //plotProton_DataMC();
    //plotPion_DataMC();
    //plotCutHistograms_DataMC();

    //--------------------------------------------------------------------------
    //  MC Only
    //--------------------------------------------------------------------------
    //plotInteraction_MCOnly();
    //plotMuon_MCOnly();
    //plotProton_MCOnly();
    //plotPion_MCOnly();
    //plotCutHistograms_MCOnly();

    //--------------------------------------------------------------------------
    //  Plot Function Reserved for Other Studies
    //--------------------------------------------------------------------------
    //SavePi0InvMassPoints();
    //plotOtherStudies();
    //plotGENIEXSec();

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
    std::cout<<"POT Data = "<<data_POT<<std::endl;
    std::cout<<"POT MC = "<<mc_POT<<std::endl;
    std::cout<<"POT Ratio = "<<POT_ratio<<std::endl;

    setRootDirs(); 
}

void CCProtonPi0_Plotter::plotCrossSection_Check()
{
    std::cout<<"Plotting Cross Section Check..."<<std::endl;
    std::string plotDir;

    // Check must use MC
    //  For Testing you can turn it off
    bool isMC = true;

    TFile* f_xsec;
    MnvH1D* data;
    MnvH1D* mc;
   
    if(isMC) f_xsec = new TFile(rootDir_CrossSection.mc.c_str());
    else f_xsec = new TFile(rootDir_CrossSection.data.c_str());

    // ------------------------------------------------------------------------
    // Pi0 Momentum
    // ------------------------------------------------------------------------
    plotDir = Folder_List::xsec_pi0_P + Folder_List::plotDir_Check;
    // Background Subtraction Check
    mc = (MnvH1D*)f_xsec->Get("pi0_P_mc_reco_bckg");
    data = (MnvH1D*)f_xsec->Get("pi0_P_bckg_estimated"); 
    if (isMC) DrawDataMC_CrossSection(data,mc,"pi0_P_check_bckg_estimated",plotDir);
    else DrawDataMC(data,mc,"pi0_P_check_bckg_estimated",plotDir);

    mc = (MnvH1D*)f_xsec->Get("pi0_P_mc_reco_signal");
    data = (MnvH1D*)f_xsec->Get("pi0_P_bckg_subtracted"); 
    if (isMC) DrawDataMC_CrossSection(data,mc,"pi0_P_check_bckg_subtracted",plotDir);
    else DrawDataMC(data,mc,"pi0_P_check_bckg_subtracted",plotDir);

    // Unfolding Check
    DrawNormalizedMigrationHistogram(rootDir_CrossSection,"pi0_P_response",plotDir);
    mc = (MnvH1D*)f_xsec->Get("pi0_P_mc_truth_signal");
    data = (MnvH1D*)f_xsec->Get("pi0_P_unfolded"); 
    if (isMC) DrawDataMC_CrossSection(data,mc,"pi0_P_check_unfolding",plotDir);
    else DrawDataMC(data,mc,"pi0_P_check_unfolding",plotDir);

    // Efficiency Check
    mc = (MnvH1D*)f_xsec->Get("pi0_P_mc_truth_all_signal");
    data = (MnvH1D*)f_xsec->Get("pi0_P_efficiency_corrected"); 
    if (isMC) DrawDataMC_CrossSection(data,mc,"pi0_P_check_efficiency",plotDir);
    else DrawDataMC(data,mc,"pi0_P_check_efficiency",plotDir);
    DrawEfficiencyCurve(rootDir_CrossSection,"pi0_P_eff",plotDir);

    // ------------------------------------------------------------------------
    // Muon Momentum
    // ------------------------------------------------------------------------
    plotDir = Folder_List::xsec_muon_P + Folder_List::plotDir_Check;
    // Background Subtraction Check
    mc = (MnvH1D*)f_xsec->Get("muon_P_mc_reco_bckg");
    data = (MnvH1D*)f_xsec->Get("muon_P_bckg_estimated"); 
    if (isMC) DrawDataMC_CrossSection(data,mc,"muon_P_check_bckg_estimated",plotDir);
    else DrawDataMC(data,mc,"muon_P_check_bckg_estimated",plotDir);

    mc = (MnvH1D*)f_xsec->Get("muon_P_mc_reco_signal");
    data = (MnvH1D*)f_xsec->Get("muon_P_bckg_subtracted"); 
    if (isMC) DrawDataMC_CrossSection(data,mc,"muon_P_check_bckg_subtracted",plotDir);
    else DrawDataMC(data,mc,"muon_P_check_bckg_subtracted",plotDir);

    // Unfolding Check
    DrawNormalizedMigrationHistogram(rootDir_CrossSection,"muon_P_response",plotDir);
    mc = (MnvH1D*)f_xsec->Get("muon_P_mc_truth_signal");
    data = (MnvH1D*)f_xsec->Get("muon_P_unfolded"); 
    if (isMC) DrawDataMC_CrossSection(data,mc,"muon_P_check_unfolding",plotDir);
    else DrawDataMC(data,mc,"muon_P_check_unfolding",plotDir);

    // Efficiency Check
    mc = (MnvH1D*)f_xsec->Get("muon_P_mc_truth_all_signal");
    data = (MnvH1D*)f_xsec->Get("muon_P_efficiency_corrected"); 
    if (isMC) DrawDataMC_CrossSection(data,mc,"muon_P_check_efficiency",plotDir);
    else DrawDataMC(data,mc,"muon_P_check_efficiency",plotDir);
    DrawEfficiencyCurve(rootDir_CrossSection,"muon_P_eff",plotDir);

    std::cout<<"Plotting Cross Section Check Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plotOtherStudies()
{
    std::cout<<"Plotting Other Studies..."<<std::endl;
    //std::string plotDir = Folder_List::plotDir_OtherStudies;

    plot_Michel_TruthMatch("time_diff");
    plot_Michel_TruthMatch("energy");
    plot_Michel_TruthMatch("distance");
    
    std::cout<<"Plotting Other Studies Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plotErrorSummary()
{
    std::cout<<"Plotting Error Summary..."<<std::endl;

    TFile* f_xsec_data = new TFile(rootDir_CrossSection.data.c_str());
    TFile* f_xsec_mc = new TFile(rootDir_CrossSection.mc.c_str());
    MnvH1D* hist;

    // Muon Momentum
    std::string plotDir = Folder_List::xsec_muon_P + Folder_List::plotDir_ErrorSummary;

    hist = (MnvH1D*)f_xsec_data->Get("muon_P_all");
    DrawErrorSummary(hist,"muon_P_all",plotDir);

    hist = (MnvH1D*)f_xsec_mc->Get("muon_P_mc_reco_all");
    DrawErrorSummary(hist,"muon_P_mc_reco_all",plotDir);

    hist = (MnvH1D*)f_xsec_mc->Get("muon_P_mc_reco_signal");
    DrawErrorSummary(hist,"muon_P_mc_reco_signal",plotDir);

    hist = (MnvH1D*)f_xsec_mc->Get("muon_P_mc_reco_bckg");
    DrawErrorSummary(hist,"muon_P_mc_reco_bckg",plotDir);

    hist = (MnvH1D*)f_xsec_mc->Get("muon_P_mc_truth_signal");
    DrawErrorSummary(hist,"muon_P_mc_truth_signal",plotDir);

    hist = (MnvH1D*)f_xsec_mc->Get("muon_P_mc_truth_signal_all");
    DrawErrorSummary(hist,"muon_P_mc_truth_signal_all",plotDir);

    TFile* f_data = new TFile(rootDir_Muon.data.c_str());
    TFile* f_mc = new TFile(rootDir_Muon.mc.c_str());
    MnvH1D* data = (MnvH1D*)f_data->Get("muon_P_all");
    MnvH1D* mc = (MnvH1D*)f_mc->Get("muon_P_mc_reco_all");
    DrawDataMC(data,mc,"muon_P_data_MC",plotDir);


    // Pi0 Momentum
    plotDir = Folder_List::xsec_pi0_P + Folder_List::plotDir_ErrorSummary;

    hist = (MnvH1D*)f_xsec_data->Get("pi0_P_all");
    DrawErrorSummary(hist,"pi0_P_all",plotDir);

    hist = (MnvH1D*)f_xsec_mc->Get("pi0_P_mc_reco_all");
    DrawErrorSummary(hist,"pi0_P_mc_reco_all",plotDir);

    hist = (MnvH1D*)f_xsec_mc->Get("pi0_P_mc_reco_signal");
    DrawErrorSummary(hist,"pi0_P_mc_reco_signal",plotDir);

    hist = (MnvH1D*)f_xsec_mc->Get("pi0_P_mc_reco_bckg");
    DrawErrorSummary(hist,"pi0_P_mc_reco_bckg",plotDir);

    hist = (MnvH1D*)f_xsec_mc->Get("pi0_P_mc_truth_signal");
    DrawErrorSummary(hist,"pi0_P_mc_truth_signal",plotDir);

    hist = (MnvH1D*)f_xsec_mc->Get("pi0_P_mc_truth_signal_all");
    DrawErrorSummary(hist,"pi0_P_mc_truth_signal_all",plotDir);

    f_data = new TFile(rootDir_Pion.data.c_str());
    f_mc = new TFile(rootDir_Pion.mc.c_str());
    data = (MnvH1D*)f_data->Get("pi0_P_all");
    mc = (MnvH1D*)f_mc->Get("pi0_P_mc_reco_all");
    DrawDataMC(data,mc,"pi0_P_data_MC",plotDir);

    std::cout<<"Plotting Error Summary Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plotIntegratedFlux()
{
    std::cout<<"Plotting Integrated Flux..."<<std::endl;
    std::string plotDir;
 
    TFile* f_xsec_mc = new TFile(rootDir_CrossSection.mc.c_str());
    TFile* f_xsec_data = new TFile(rootDir_CrossSection.data.c_str());
    MnvH1D* mc;
    MnvH1D* data;

    // QSq 
    plotDir = Folder_List::xsec_QSq + Folder_List::plotDir_IntegratedFlux;
    mc = (MnvH1D*)f_xsec_mc->Get("QSq_integrated_flux");
    data = (MnvH1D*)f_xsec_data->Get("QSq_integrated_flux"); 
    DrawDataMC(data,mc,"QSq_data_MC",plotDir);
    DrawErrorSummary(data,"QSq_integrated_flux",plotDir);
    DrawErrorSummary(mc,"QSq_mc_integrated_flux",plotDir);

    // Pi0 Momentum 
    plotDir = Folder_List::xsec_pi0_P + Folder_List::plotDir_IntegratedFlux;
    mc = (MnvH1D*)f_xsec_mc->Get("pi0_P_integrated_flux");
    data = (MnvH1D*)f_xsec_data->Get("pi0_P_integrated_flux"); 
    DrawDataMC(data,mc,"pi0_P_data_MC",plotDir);
    DrawErrorSummary(data,"pi0_P_integrated_flux",plotDir);
    DrawErrorSummary(mc,"pi0_P_mc_integrated_flux",plotDir);

    // Pi0 Kinetic Energy
    plotDir = Folder_List::xsec_pi0_KE + Folder_List::plotDir_IntegratedFlux;
    mc = (MnvH1D*)f_xsec_mc->Get("pi0_KE_integrated_flux");
    data = (MnvH1D*)f_xsec_data->Get("pi0_KE_integrated_flux"); 
    DrawDataMC(data,mc,"pi0_KE_data_MC",plotDir);
    DrawErrorSummary(data,"pi0_KE_integrated_flux",plotDir);
    DrawErrorSummary(mc,"pi0_KE_mc_integrated_flux",plotDir);

    // Pi0 Theta 
    plotDir = Folder_List::xsec_pi0_theta + Folder_List::plotDir_IntegratedFlux;
    mc = (MnvH1D*)f_xsec_mc->Get("pi0_theta_integrated_flux");
    data = (MnvH1D*)f_xsec_data->Get("pi0_theta_integrated_flux"); 
    DrawDataMC(data,mc,"pi0_theta_data_MC",plotDir);
    DrawErrorSummary(data,"pi0_theta_integrated_flux",plotDir);
    DrawErrorSummary(mc,"pi0_theta_mc_integrated_flux",plotDir);

    // Muon Momentum
    plotDir = Folder_List::xsec_muon_P + Folder_List::plotDir_IntegratedFlux;
    mc = (MnvH1D*)f_xsec_mc->Get("muon_P_integrated_flux");
    data = (MnvH1D*)f_xsec_data->Get("muon_P_integrated_flux"); 
    DrawDataMC(data,mc,"muon_P_data_MC",plotDir);
    DrawErrorSummary(data,"muon_P_integrated_flux",plotDir);
    DrawErrorSummary(mc,"muon_P_mc_integrated_flux",plotDir);

    // Muon Theta 
    plotDir = Folder_List::xsec_muon_theta + Folder_List::plotDir_IntegratedFlux;
    mc = (MnvH1D*)f_xsec_mc->Get("muon_theta_integrated_flux");
    data = (MnvH1D*)f_xsec_data->Get("muon_theta_integrated_flux"); 
    DrawDataMC(data,mc,"muon_theta_data_MC",plotDir);
    DrawErrorSummary(data,"muon_theta_integrated_flux",plotDir);
    DrawErrorSummary(mc,"muon_theta_mc_integrated_flux",plotDir);

    delete f_xsec_mc;
    delete f_xsec_data;
    std::cout<<"Plotting Original Data Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plotOriginalData()
{
    std::cout<<"Plotting Original Data..."<<std::endl;
    std::string plotDir;
 
    TFile* f_xsec_mc = new TFile(rootDir_CrossSection.mc.c_str());
    TFile* f_xsec_data = new TFile(rootDir_CrossSection.data.c_str());
    MnvH1D* mc;
    MnvH1D* data;

    // QSq 
    plotDir = Folder_List::xsec_QSq + Folder_List::plotDir_Original;
    mc = (MnvH1D*)f_xsec_mc->Get("QSq_mc_reco_all");
    data = (MnvH1D*)f_xsec_data->Get("QSq_all"); 
    DrawDataMC(data,mc,"QSq_data_MC",plotDir);
    DrawErrorSummary(data,"QSq_all",plotDir);
    DrawErrorSummary(mc,"QSq_mc_reco_all",plotDir);

    // Pi0 Momentum 
    plotDir = Folder_List::xsec_pi0_P + Folder_List::plotDir_Original;
    mc = (MnvH1D*)f_xsec_mc->Get("pi0_P_mc_reco_all");
    data = (MnvH1D*)f_xsec_data->Get("pi0_P_all"); 
    DrawDataMC(data,mc,"pi0_P_data_MC",plotDir);
    DrawErrorSummary(data,"pi0_P_all",plotDir);
    DrawErrorSummary(mc,"pi0_P_mc_reco_all",plotDir);

    // Pi0 Kinetic Energy
    plotDir = Folder_List::xsec_pi0_KE + Folder_List::plotDir_Original;
    mc = (MnvH1D*)f_xsec_mc->Get("pi0_KE_mc_reco_all");
    data = (MnvH1D*)f_xsec_data->Get("pi0_KE_all"); 
    DrawDataMC(data,mc,"pi0_KE_data_MC",plotDir);
    DrawErrorSummary(data,"pi0_KE_all",plotDir);
    DrawErrorSummary(mc,"pi0_KE_mc_reco_all",plotDir);

    // Pi0 Theta 
    plotDir = Folder_List::xsec_pi0_theta + Folder_List::plotDir_Original;
    mc = (MnvH1D*)f_xsec_mc->Get("pi0_theta_mc_reco_all");
    data = (MnvH1D*)f_xsec_data->Get("pi0_theta_all"); 
    DrawDataMC(data,mc,"pi0_theta_data_MC",plotDir);
    DrawErrorSummary(data,"pi0_theta_all",plotDir);
    DrawErrorSummary(mc,"pi0_theta_mc_reco_all",plotDir);

    // Muon Momentum
    plotDir = Folder_List::xsec_muon_P + Folder_List::plotDir_Original;
    mc = (MnvH1D*)f_xsec_mc->Get("muon_P_mc_reco_all");
    data = (MnvH1D*)f_xsec_data->Get("muon_P_all"); 
    DrawDataMC(data,mc,"muon_P_data_MC",plotDir);
    DrawErrorSummary(data,"muon_P_all",plotDir);
    DrawErrorSummary(mc,"muon_P_mc_reco_all",plotDir);

    // Muon Theta 
    plotDir = Folder_List::xsec_muon_theta + Folder_List::plotDir_Original;
    mc = (MnvH1D*)f_xsec_mc->Get("muon_theta_mc_reco_all");
    data = (MnvH1D*)f_xsec_data->Get("muon_theta_all"); 
    DrawDataMC(data,mc,"muon_theta_data_MC",plotDir);
    DrawErrorSummary(data,"muon_theta_all",plotDir);
    DrawErrorSummary(mc,"muon_theta_mc_reco_all",plotDir);

    delete f_xsec_mc;
    delete f_xsec_data;
    std::cout<<"Plotting Original Data Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plotCrossSection()
{
    std::cout<<"Plotting Cross Sections..."<<std::endl;

    std::string plotDir;
  
    TFile* f_xsec_mc = new TFile(rootDir_CrossSection.mc.c_str());
    TFile* f_xsec_data = new TFile(rootDir_CrossSection.data.c_str());
    MnvH1D* mc;
    MnvH1D* data;

    // QSq
    plotDir = Folder_List::xsec_QSq + Folder_List::plotDir_CrossSection;
    mc = (MnvH1D*)f_xsec_mc->Get("QSq_xsec");
    data = (MnvH1D*)f_xsec_data->Get("QSq_xsec"); 
    DrawDataMC_CrossSection(data,mc,"QSq_xsec_data_MC",plotDir);
    DrawErrorSummary(data,"QSq_xsec_data",plotDir);
    DrawErrorSummary(mc,"QSq_xsec_mc",plotDir);

    // Pi0 Momentum 
    plotDir = Folder_List::xsec_pi0_P + Folder_List::plotDir_CrossSection;
    mc = (MnvH1D*)f_xsec_mc->Get("pi0_P_xsec");
    data = (MnvH1D*)f_xsec_data->Get("pi0_P_xsec"); 
    DrawDataMC_CrossSection(data,mc,"pi0_P_xsec_data_MC",plotDir);
    DrawErrorSummary(data,"pi0_P_xsec_data",plotDir);
    DrawErrorSummary(mc,"pi0_P_xsec_mc",plotDir);

    // Pi0 Kinetic Energy 
    plotDir = Folder_List::xsec_pi0_KE + Folder_List::plotDir_CrossSection;
    mc = (MnvH1D*)f_xsec_mc->Get("pi0_KE_xsec");
    data = (MnvH1D*)f_xsec_data->Get("pi0_KE_xsec"); 
    DrawDataMC_CrossSection(data,mc,"pi0_KE_xsec_data_MC",plotDir);
    DrawErrorSummary(data,"pi0_KE_xsec_data",plotDir);
    DrawErrorSummary(mc,"pi0_KE_xsec_mc",plotDir);

    // Pi0 Theta
    plotDir = Folder_List::xsec_pi0_theta + Folder_List::plotDir_CrossSection;
    mc = (MnvH1D*)f_xsec_mc->Get("pi0_theta_xsec");
    data = (MnvH1D*)f_xsec_data->Get("pi0_theta_xsec"); 
    DrawDataMC_CrossSection(data,mc,"pi0_theta_xsec_data_MC",plotDir);
    DrawErrorSummary(data,"pi0_theta_xsec_data",plotDir);
    DrawErrorSummary(mc,"pi0_theta_xsec_mc",plotDir);

    // Muon Momentum
    plotDir = Folder_List::xsec_muon_P + Folder_List::plotDir_CrossSection;
    mc = (MnvH1D*)f_xsec_mc->Get("muon_P_xsec");
    data = (MnvH1D*)f_xsec_data->Get("muon_P_xsec"); 
    DrawDataMC_CrossSection(data,mc,"muon_P_xsec_data_MC",plotDir);
    DrawErrorSummary(data,"muon_P_xsec_data",plotDir);
    DrawErrorSummary(mc,"muon_P_xsec_mc",plotDir);

    // Muon Theta 
    plotDir = Folder_List::xsec_muon_theta + Folder_List::plotDir_CrossSection;
    mc = (MnvH1D*)f_xsec_mc->Get("muon_theta_xsec");
    data = (MnvH1D*)f_xsec_data->Get("muon_theta_xsec"); 
    DrawDataMC_CrossSection(data,mc,"muon_theta_xsec_data_MC",plotDir);
    DrawErrorSummary(data,"muon_theta_xsec_data",plotDir);
    DrawErrorSummary(mc,"muon_theta_xsec_mc",plotDir);

    delete f_xsec_mc;
    delete f_xsec_data;
    std::cout<<"Plotting Cross Sections Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plotBackgroundEstimated()
{
    std::cout<<"Plotting Background Estimated..."<<std::endl;
    std::string plotDir;

    TFile* f_xsec_mc = new TFile(rootDir_CrossSection.mc.c_str());
    TFile* f_xsec_data = new TFile(rootDir_CrossSection.data.c_str());
    MnvH1D* mc;
    MnvH1D* data;

    // QSq 
    plotDir = Folder_List::xsec_QSq + Folder_List::plotDir_BackgroundEstimated;
    mc = (MnvH1D*)f_xsec_mc->Get("QSq_bckg_estimated");
    data = (MnvH1D*)f_xsec_data->Get("QSq_bckg_estimated"); 
    DrawDataMC(data,mc,"QSq_bckg_estimated_data_MC",plotDir);
    DrawErrorSummary(data,"QSq_bckg_estimated_data",plotDir);
    DrawErrorSummary(mc,"QSq_bckg_estimated_mc",plotDir);

    // Pi0 Momentum 
    plotDir = Folder_List::xsec_pi0_P + Folder_List::plotDir_BackgroundEstimated;
    mc = (MnvH1D*)f_xsec_mc->Get("pi0_P_bckg_estimated");
    data = (MnvH1D*)f_xsec_data->Get("pi0_P_bckg_estimated"); 
    DrawDataMC(data,mc,"pi0_P_bckg_estimated_data_MC",plotDir);
    DrawErrorSummary(data,"pi0_P_bckg_estimated_data",plotDir);
    DrawErrorSummary(mc,"pi0_P_bckg_estimated_mc",plotDir);

    // Pi0 Kinetic Energy 
    plotDir = Folder_List::xsec_pi0_KE + Folder_List::plotDir_BackgroundEstimated;
    mc = (MnvH1D*)f_xsec_mc->Get("pi0_KE_bckg_estimated");
    data = (MnvH1D*)f_xsec_data->Get("pi0_KE_bckg_estimated"); 
    DrawDataMC(data,mc,"pi0_KE_bckg_estimated_data_MC",plotDir);
    DrawErrorSummary(data,"pi0_KE_bckg_estimated_data",plotDir);
    DrawErrorSummary(mc,"pi0_KE_bckg_estimated_mc",plotDir);
    
    // Pi0 Theta 
    plotDir = Folder_List::xsec_pi0_theta + Folder_List::plotDir_BackgroundEstimated;
    mc = (MnvH1D*)f_xsec_mc->Get("pi0_theta_bckg_estimated");
    data = (MnvH1D*)f_xsec_data->Get("pi0_theta_bckg_estimated"); 
    DrawDataMC(data,mc,"pi0_theta_bckg_estimated_data_MC",plotDir);
    DrawErrorSummary(data,"pi0_theta_bckg_estimated_data",plotDir);
    DrawErrorSummary(mc,"pi0_theta_bckg_estimated_mc",plotDir);

    // Muon Momentum
    plotDir = Folder_List::xsec_muon_P + Folder_List::plotDir_BackgroundEstimated;
    mc = (MnvH1D*)f_xsec_mc->Get("muon_P_bckg_estimated");
    data = (MnvH1D*)f_xsec_data->Get("muon_P_bckg_estimated"); 
    DrawDataMC(data,mc,"muon_P_bckg_estimated_data_MC",plotDir);
    DrawErrorSummary(data,"muon_P_bckg_estimated_data",plotDir);
    DrawErrorSummary(mc,"muon_P_bckg_estimated_mc",plotDir);

    // Muon Theta 
    plotDir = Folder_List::xsec_muon_theta + Folder_List::plotDir_BackgroundEstimated;
    mc = (MnvH1D*)f_xsec_mc->Get("muon_theta_bckg_estimated");
    data = (MnvH1D*)f_xsec_data->Get("muon_theta_bckg_estimated"); 
    DrawDataMC(data,mc,"muon_theta_bckg_estimated_data_MC",plotDir);
    DrawErrorSummary(data,"muon_theta_bckg_estimated_data",plotDir);
    DrawErrorSummary(mc,"muon_theta_bckg_estimated_mc",plotDir);

    delete f_xsec_mc;
    delete f_xsec_data;
    std::cout<<"Plotting Background Estimated Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plotBackgroundSubtracted()
{
    std::cout<<"Plotting Background Subtracted..."<<std::endl;
    std::string plotDir;
  
    TFile* f_xsec_mc = new TFile(rootDir_CrossSection.mc.c_str());
    TFile* f_xsec_data = new TFile(rootDir_CrossSection.data.c_str());
    MnvH1D* mc;
    MnvH1D* data;

    // QSq 
    plotDir = Folder_List::xsec_QSq + Folder_List::plotDir_BackgroundSubtracted;
    mc = (MnvH1D*)f_xsec_mc->Get("QSq_bckg_subtracted");
    data = (MnvH1D*)f_xsec_data->Get("QSq_bckg_subtracted"); 
    DrawDataMC(data,mc,"QSq_bckg_subtracted_data_MC",plotDir);
    DrawErrorSummary(data,"QSq_bckg_subtracted_data",plotDir);
    DrawErrorSummary(mc,"QSq_bckg_subtracted_mc",plotDir);

    // Pi0 Momentum 
    plotDir = Folder_List::xsec_pi0_P + Folder_List::plotDir_BackgroundSubtracted;
    mc = (MnvH1D*)f_xsec_mc->Get("pi0_P_bckg_subtracted");
    data = (MnvH1D*)f_xsec_data->Get("pi0_P_bckg_subtracted"); 
    DrawDataMC(data,mc,"pi0_P_bckg_subtracted_data_MC",plotDir);
    DrawErrorSummary(data,"pi0_P_bckg_subtracted_data",plotDir);
    DrawErrorSummary(mc,"pi0_P_bckg_subtracted_mc",plotDir);
 
    // Pi0 Kinetic Energy 
    plotDir = Folder_List::xsec_pi0_KE + Folder_List::plotDir_BackgroundSubtracted;
    mc = (MnvH1D*)f_xsec_mc->Get("pi0_KE_bckg_subtracted");
    data = (MnvH1D*)f_xsec_data->Get("pi0_KE_bckg_subtracted"); 
    DrawDataMC(data,mc,"pi0_KE_bckg_subtracted_data_MC",plotDir);
    DrawErrorSummary(data,"pi0_KE_bckg_subtracted_data",plotDir);
    DrawErrorSummary(mc,"pi0_KE_bckg_subtracted_mc",plotDir);
 
    // Pi0 Theta
    plotDir = Folder_List::xsec_pi0_theta + Folder_List::plotDir_BackgroundSubtracted;
    mc = (MnvH1D*)f_xsec_mc->Get("pi0_theta_bckg_subtracted");
    data = (MnvH1D*)f_xsec_data->Get("pi0_theta_bckg_subtracted"); 
    DrawDataMC(data,mc,"pi0_theta_bckg_subtracted_data_MC",plotDir);
    DrawErrorSummary(data,"pi0_theta_bckg_subtracted_data",plotDir);
    DrawErrorSummary(mc,"pi0_theta_bckg_subtracted_mc",plotDir);
 
    // Muon Momentum
    plotDir = Folder_List::xsec_muon_P + Folder_List::plotDir_BackgroundSubtracted;
    mc = (MnvH1D*)f_xsec_mc->Get("muon_P_bckg_subtracted");
    data = (MnvH1D*)f_xsec_data->Get("muon_P_bckg_subtracted"); 
    DrawDataMC(data,mc,"muon_P_bckg_subtracted_data_MC",plotDir);
    DrawErrorSummary(data,"muon_P_bckg_subtracted_data",plotDir);
    DrawErrorSummary(mc,"muon_P_bckg_subtracted_mc",plotDir);
 
    // Muon Theta 
    plotDir = Folder_List::xsec_muon_theta + Folder_List::plotDir_BackgroundSubtracted;
    mc = (MnvH1D*)f_xsec_mc->Get("muon_theta_bckg_subtracted");
    data = (MnvH1D*)f_xsec_data->Get("muon_theta_bckg_subtracted"); 
    DrawDataMC(data,mc,"muon_theta_bckg_subtracted_data_MC",plotDir);
    DrawErrorSummary(data,"muon_theta_bckg_subtracted_data",plotDir);
    DrawErrorSummary(mc,"muon_theta_bckg_subtracted_mc",plotDir);

    delete f_xsec_mc;
    delete f_xsec_data;
    std::cout<<"Plotting Background Subtracted Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plotUnfolded()
{
    std::cout<<"Plotting Unfolded..."<<std::endl;
    std::string plotDir;
  
    TFile* f_xsec_mc = new TFile(rootDir_CrossSection.mc.c_str());
    TFile* f_xsec_data = new TFile(rootDir_CrossSection.data.c_str());
    MnvH1D* mc;
    MnvH1D* data;
    
    // QSq 
    plotDir = Folder_List::xsec_QSq + Folder_List::plotDir_Unfolding;
    mc = (MnvH1D*)f_xsec_mc->Get("QSq_unfolded");
    data = (MnvH1D*)f_xsec_data->Get("QSq_unfolded"); 
    DrawDataMC(data,mc,"QSq_unfolded_data_MC",plotDir);
    DrawNormalizedMigrationHistogram(rootDir_CrossSection,"QSq_response",plotDir);
    DrawErrorSummary(data,"QSq_unfolded_data",plotDir);
    DrawErrorSummary(mc,"QSq_unfolded_mc",plotDir);

    // Pi0 Momentum 
    plotDir = Folder_List::xsec_pi0_P + Folder_List::plotDir_Unfolding;
    mc = (MnvH1D*)f_xsec_mc->Get("pi0_P_unfolded");
    data = (MnvH1D*)f_xsec_data->Get("pi0_P_unfolded"); 
    DrawDataMC(data,mc,"pi0_P_unfolded_data_MC",plotDir);
    DrawNormalizedMigrationHistogram(rootDir_CrossSection,"pi0_P_response",plotDir);
    DrawErrorSummary(data,"pi0_P_unfolded_data",plotDir);
    DrawErrorSummary(mc,"pi0_P_unfolded_mc",plotDir);
 
    // Pi0 Kinetic Energy
    plotDir = Folder_List::xsec_pi0_KE + Folder_List::plotDir_Unfolding;
    mc = (MnvH1D*)f_xsec_mc->Get("pi0_KE_unfolded");
    data = (MnvH1D*)f_xsec_data->Get("pi0_KE_unfolded"); 
    DrawDataMC(data,mc,"pi0_KE_unfolded_data_MC",plotDir);
    DrawNormalizedMigrationHistogram(rootDir_CrossSection,"pi0_KE_response",plotDir);
    DrawErrorSummary(data,"pi0_KE_unfolded_data",plotDir);
    DrawErrorSummary(mc,"pi0_KE_unfolded_mc",plotDir);
 
    // Pi0 Theta 
    plotDir = Folder_List::xsec_pi0_theta + Folder_List::plotDir_Unfolding;
    mc = (MnvH1D*)f_xsec_mc->Get("pi0_theta_unfolded");
    data = (MnvH1D*)f_xsec_data->Get("pi0_theta_unfolded"); 
    DrawDataMC(data,mc,"pi0_theta_unfolded_data_MC",plotDir);
    DrawNormalizedMigrationHistogram(rootDir_CrossSection,"pi0_theta_response",plotDir);
    DrawErrorSummary(data,"pi0_theta_unfolded_data",plotDir);
    DrawErrorSummary(mc,"pi0_theta_unfolded_mc",plotDir);
  
    // Muon Momentum
    plotDir = Folder_List::xsec_muon_P + Folder_List::plotDir_Unfolding;
    mc = (MnvH1D*)f_xsec_mc->Get("muon_P_unfolded");
    data = (MnvH1D*)f_xsec_data->Get("muon_P_unfolded"); 
    DrawDataMC(data,mc,"muon_P_unfolded_data_MC",plotDir);
    DrawNormalizedMigrationHistogram(rootDir_CrossSection,"muon_P_response",plotDir);
    DrawErrorSummary(data,"muon_P_unfolded_data",plotDir);
    DrawErrorSummary(mc,"muon_P_unfolded_mc",plotDir);
     
    // Muon Theta 
    plotDir = Folder_List::xsec_muon_theta + Folder_List::plotDir_Unfolding;
    mc = (MnvH1D*)f_xsec_mc->Get("muon_theta_unfolded");
    data = (MnvH1D*)f_xsec_data->Get("muon_theta_unfolded"); 
    DrawDataMC(data,mc,"muon_theta_unfolded_data_MC",plotDir);
    DrawNormalizedMigrationHistogram(rootDir_CrossSection,"muon_theta_response",plotDir);
    DrawErrorSummary(data,"muon_theta_unfolded_data",plotDir);
    DrawErrorSummary(mc,"muon_theta_unfolded_mc",plotDir);

    delete f_xsec_mc;
    delete f_xsec_data;
    std::cout<<"Plotting Unfolded Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plotEfficiencyCorrected()
{
    std::cout<<"Plotting Efficiency Corrected..."<<std::endl;
    std::string plotDir;
  
    TFile* f_xsec_mc = new TFile(rootDir_CrossSection.mc.c_str());
    TFile* f_xsec_data = new TFile(rootDir_CrossSection.data.c_str());
    MnvH1D* mc;
    MnvH1D* data;

    // QSq 
    plotDir = Folder_List::xsec_QSq + Folder_List::plotDir_Efficiency;
    mc = (MnvH1D*)f_xsec_mc->Get("QSq_efficiency_corrected");
    data = (MnvH1D*)f_xsec_data->Get("QSq_efficiency_corrected"); 
    DrawDataMC(data,mc,"QSq_efficiency_corrected_data_MC",plotDir);
    DrawEfficiencyCurve(rootDir_CrossSection,"QSq_eff",plotDir);
    DrawErrorSummary(data,"QSq_efficiency_corrected_data",plotDir);
    DrawErrorSummary(mc,"QSq_efficiency_corrected_mc",plotDir);

    // Pi0 Momentum 
    plotDir = Folder_List::xsec_pi0_P + Folder_List::plotDir_Efficiency;
    mc = (MnvH1D*)f_xsec_mc->Get("pi0_P_efficiency_corrected");
    data = (MnvH1D*)f_xsec_data->Get("pi0_P_efficiency_corrected"); 
    DrawDataMC(data,mc,"pi0_P_efficiency_corrected_data_MC",plotDir);
    DrawEfficiencyCurve(rootDir_CrossSection,"pi0_P_eff",plotDir);
    DrawErrorSummary(data,"pi0_P_efficiency_corrected_data",plotDir);
    DrawErrorSummary(mc,"pi0_P_efficiency_corrected_mc",plotDir);

    // Pi0 Kinetic Energy 
    plotDir = Folder_List::xsec_pi0_KE + Folder_List::plotDir_Efficiency;
    mc = (MnvH1D*)f_xsec_mc->Get("pi0_KE_efficiency_corrected");
    data = (MnvH1D*)f_xsec_data->Get("pi0_KE_efficiency_corrected"); 
    DrawDataMC(data,mc,"pi0_KE_efficiency_corrected_data_MC",plotDir);
    DrawEfficiencyCurve(rootDir_CrossSection,"pi0_KE_eff",plotDir);
    DrawErrorSummary(data,"pi0_KE_efficiency_corrected_data",plotDir);
    DrawErrorSummary(mc,"pi0_KE_efficiency_corrected_mc",plotDir);

    // Pi0 Theta
    plotDir = Folder_List::xsec_pi0_theta + Folder_List::plotDir_Efficiency;
    mc = (MnvH1D*)f_xsec_mc->Get("pi0_theta_efficiency_corrected");
    data = (MnvH1D*)f_xsec_data->Get("pi0_theta_efficiency_corrected"); 
    DrawDataMC(data,mc,"pi0_theta_efficiency_corrected_data_MC",plotDir);
    DrawEfficiencyCurve(rootDir_CrossSection,"pi0_theta_eff",plotDir);
    DrawErrorSummary(data,"pi0_theta_efficiency_corrected_data",plotDir);
    DrawErrorSummary(mc,"pi0_theta_efficiency_corrected_mc",plotDir);

    // Muon Momentum
    plotDir = Folder_List::xsec_muon_P + Folder_List::plotDir_Efficiency;
    mc = (MnvH1D*)f_xsec_mc->Get("muon_P_efficiency_corrected");
    data = (MnvH1D*)f_xsec_data->Get("muon_P_efficiency_corrected"); 
    DrawDataMC(data,mc,"muon_P_efficiency_corrected_data_MC",plotDir);
    DrawEfficiencyCurve(rootDir_CrossSection,"muon_P_eff",plotDir);
    DrawErrorSummary(data,"muon_P_efficiency_corrected_data",plotDir);
    DrawErrorSummary(mc,"muon_P_efficiency_corrected_mc",plotDir);

    // Muon Theta 
    plotDir = Folder_List::xsec_muon_theta + Folder_List::plotDir_Efficiency;
    mc = (MnvH1D*)f_xsec_mc->Get("muon_theta_efficiency_corrected");
    data = (MnvH1D*)f_xsec_data->Get("muon_theta_efficiency_corrected"); 
    DrawDataMC(data,mc,"muon_theta_efficiency_corrected_data_MC",plotDir);
    DrawEfficiencyCurve(rootDir_CrossSection,"muon_theta_eff",plotDir);
    DrawErrorSummary(data,"muon_theta_efficiency_corrected_data",plotDir);
    DrawErrorSummary(mc,"muon_theta_efficiency_corrected_mc",plotDir);


    delete f_xsec_mc;
    delete f_xsec_data;
    std::cout<<"Plotting Efficiency Corrected Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plotInteraction_MCOnly()
{
    std::cout<<"Plotting Interaction MC Only"<<std::endl;
    std::string plotDir = Folder_List::plotDir_Interaction;

    Draw1DHist(rootDir_Interaction,"QSq_True",plotDir);
    Draw1DHist(rootDir_Interaction,"QSq_Error",plotDir);
    Draw1DHist(rootDir_Interaction,"QSq_Diff",plotDir);
    
    //plot_mc_w_Stacked();
    //plot_final_mc_w_Stacked();

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

    Draw1DHist(rootDir_Interaction,"Enu_True_1Track",plotDir);
    Draw1DHist(rootDir_Interaction,"Enu_True_2Track",plotDir);
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

   // DrawStackedMC(rootDir_Pion,"gamma1_E",plotDir);
   // DrawStackedMC(rootDir_Pion,"gamma1_theta",plotDir);
   // DrawStackedMC(rootDir_Pion,"gamma1_ConvLength",plotDir);

   // DrawStackedMC(rootDir_Pion,"gamma2_E",plotDir);
   // DrawStackedMC(rootDir_Pion,"gamma2_theta",plotDir);
   // DrawStackedMC(rootDir_Pion,"gamma2_ConvLength",plotDir);

   // DrawStackedMC(rootDir_Pion,"photonEnergy_Asymmetry",plotDir);
   // DrawStackedMC(rootDir_Pion,"invMass",plotDir);
   // DrawStackedMC(rootDir_Pion,"cos_openingAngle",plotDir);
   // DrawStackedMC(rootDir_Pion,"P",plotDir);
   // DrawStackedMC(rootDir_Pion,"E",plotDir);
   // DrawStackedMC(rootDir_Pion,"theta",plotDir);
   // DrawStackedMC(rootDir_Pion,"phi",plotDir);

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

//    Draw1DHist(rootDir_Pion,"gamma1_true_E",plotDir);
//    Draw1DHist(rootDir_Pion,"gamma1_reco_error_E",plotDir);
//    Draw2DHist(rootDir_Pion,"gamma1_true_E_reco_E_error",plotDir);
//    Draw2DHist(rootDir_Pion,"gamma1_reco_E_true_E",plotDir);
//
//    Draw1DHist(rootDir_Pion,"gamma2_true_E",plotDir);
//    Draw1DHist(rootDir_Pion,"gamma2_reco_error_E",plotDir);
//    Draw2DHist(rootDir_Pion,"gamma2_true_E_reco_E_error",plotDir);
//    Draw2DHist(rootDir_Pion,"gamma2_reco_E_true_E",plotDir);
//    
//    Draw2DHist(rootDir_Pion,"signal_gamma1_E_gamma2_E",plotDir);
//    Draw2DHist(rootDir_Pion,"bckg_gamma1_E_gamma2_E",plotDir);
//    Draw2DHist(rootDir_Pion,"bckg_signal_diff_E",plotDir);
//
//    Draw2DHist(rootDir_Pion,"signal_gamma1_convLength_gamma2_convLength",plotDir);
//    Draw2DHist(rootDir_Pion,"bckg_gamma1_convLength_gamma2_convLength",plotDir);
//    Draw2DHist(rootDir_Pion,"bckg_signal_diff_convLength",plotDir);

    // Momentum
    Draw2DHist(rootDir_Pion,"reco_P_true_P",plotDir);
    Draw1DHist(rootDir_Pion,"P_error",plotDir);
    Draw1DHist(rootDir_Pion,"P_Diff",plotDir);
    
    // Energy
    Draw2DHist(rootDir_Pion,"reco_E_true_E",plotDir);
    Draw1DHist(rootDir_Pion,"E_error",plotDir);
    Draw1DHist(rootDir_Pion,"E_Diff",plotDir);
    Draw1DHist(rootDir_Pion,"E_true",plotDir);
    Draw1DHist(rootDir_Pion,"E_reco",plotDir);
 
    // Kinetic Energy
    Draw2DHist(rootDir_Pion,"reco_KE_true_KE",plotDir);
    Draw1DHist(rootDir_Pion,"KE_error",plotDir);
    Draw1DHist(rootDir_Pion,"KE_Diff",plotDir);
    Draw1DHist(rootDir_Pion,"KE_true",plotDir);
    Draw1DHist(rootDir_Pion,"KE_reco",plotDir);

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
    //CutArrow Vertex_Count(3,"L"); 
    //DrawDataStackedMC(rootDir_CutHists,"hCut_nVertices",plotDir, 1, Vertex_Count);
    //   
    //CutArrow Proton_Count(3,"L"); 
    //DrawDataStackedMC(rootDir_CutHists,"hCut_nProtonCandidates",plotDir,1 , Proton_Count);
    //
    //DrawDataStackedMC(rootDir_CutHists,"hCut_nShowerCandidates",plotDir);
    //DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_nShowerCandidates",plotDir);
    //DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_nShowerCandidates",plotDir);

    CutArrow Michel(1,"L"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_Michel",plotDir,1, Michel);

   
    plot_InvMass_TruthMatch_Stacked(true,true);
    plot_InvMass_TruthMatch_Stacked(true,false);
    plot_InvMass_TruthMatch_Stacked(false,true);
    plot_InvMass_TruthMatch_Stacked(false,false);
    CutArrow pi0invMass_min(60,"R"); 
    CutArrow pi0invMass_max(200,"L"); 
    DrawSignalMC(rootDir_CutHists,"hCut_pi0invMass",plotDir, 2, pi0invMass_min, pi0invMass_max);
    DrawDataStackedMC(rootDir_CutHists,"hCut_pi0invMass",plotDir, 2, pi0invMass_min, pi0invMass_max);

    // ------------------------------------------------------------------------
    // 1 Track
    // ------------------------------------------------------------------------
    //CutArrow eVis_nuclearTarget_1Track(20,"L"); 
    //DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_eVis_nuclearTarget",plotDir, 1, eVis_nuclearTarget_1Track);

    //CutArrow eVis_other_min_1Track(50,"R"); 
    //CutArrow eVis_other_max_1Track(2000,"L"); 
    //DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_eVis_other",plotDir, 2, eVis_other_min_1Track, eVis_other_max_1Track);

    //CutArrow gamma1_ConvDist_1Track(14,"R"); 
    //DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_gamma1ConvDist",plotDir, 1, gamma1_ConvDist_1Track);

    //CutArrow gamma2_ConvDist_1Track(15,"R"); 
    //DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_gamma2ConvDist",plotDir, 1, gamma2_ConvDist_1Track);
    //DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_gamma2ConvDist",plotDir);

    DrawSignalMC(rootDir_CutHists,"hCut_1Track_pi0invMass",plotDir, 2, pi0invMass_min, pi0invMass_max);
    DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_pi0invMass",plotDir, 2, pi0invMass_min, pi0invMass_max);

    //CutArrow neutrinoE_1Track(20,"L"); 
    //DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_neutrinoE",plotDir, 1, neutrinoE_1Track);

    // ------------------------------------------------------------------------
    // 2 Track
    // ------------------------------------------------------------------------
    //CutArrow eVis_nuclearTarget_2Track(20,"L"); 
    //DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_eVis_nuclearTarget",plotDir, 1, eVis_nuclearTarget_2Track);

    //CutArrow eVis_other_min_2Track(50,"R"); 
    //CutArrow eVis_other_max_2Track(2000,"L"); 
    //DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_eVis_other",plotDir, 2, eVis_other_min_2Track, eVis_other_max_2Track);

    //CutArrow gamma1_ConvDist_2Track(14,"R"); 
    //DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_gamma1ConvDist",plotDir, 1, gamma1_ConvDist_2Track);

    //CutArrow gamma2_ConvDist_2Track(15,"R"); 
    //DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_gamma2ConvDist",plotDir, 1, gamma2_ConvDist_2Track);
    //DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_gamma2ConvDist",plotDir);

    DrawSignalMC(rootDir_CutHists,"hCut_2Track_pi0invMass",plotDir, 2, pi0invMass_min, pi0invMass_max);
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_pi0invMass",plotDir, 2, pi0invMass_min, pi0invMass_max);

    //CutArrow neutrinoE_2Track(20,"L"); 
    //DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_neutrinoE",plotDir, 1, neutrinoE_2Track);

    //CutArrow protonScore_pIDDiff(0.45,"R"); 
    //DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_protonScore_pIDDiff",plotDir, 1, protonScore_pIDDiff);

    CutArrow protonScore_LLR(-10,"R"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_protonScore_LLR",plotDir, 1, protonScore_LLR);

    //DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_deltaInvMass",plotDir);

    std::cout<<"Plotting CutHistograms Data vs MC Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plot_InvMass_TruthMatch_Stacked(bool isSignal, bool isStacked)
{
    std::cout<<"Plottting InvMass Truth Match"<<std::endl;
    
    if (isSignal) std::cout<<"\tSignal";
    else std::cout<<"\tBackground";
   
    if (isStacked) std::cout<<" - Stacked"<<std::endl;
    else std::cout<<" - Non Stacked"<<std::endl;

    std::string root_dir = rootDir_CutHists.mc;
    std::string plotDir = Folder_List::plotDir_CutHists;
    std::string var;
    std::string type;
    if (isSignal) type = "signal_";
    else type = "background_";
    
    TFile* f_Root = new TFile(root_dir.c_str());
    TCanvas* c1 = new TCanvas("c","c",1280,800);
    THStack *hs = new THStack("hs","m_{#gamma#gamma} with Truth Match");
   
    var = type + "invMass_pizero";
    TH1D* h_pizero = (TH1D*)f_Root->Get(var.c_str());
    h_pizero->SetFillColor(kGreen);
    h_pizero->SetMarkerColor(kGreen);

    var = type + "invMass_piplus";
    TH1D* h_piplus = (TH1D*)f_Root->Get(var.c_str());
    h_piplus->SetFillColor(kRed);
    h_piplus->SetMarkerColor(kRed);

    var = type + "invMass_proton";
    TH1D* h_proton = (TH1D*)f_Root->Get(var.c_str());
    h_proton->SetFillColor(kOrange);
    h_proton->SetMarkerColor(kOrange);

    var = type + "invMass_neutron";
    TH1D* h_neutron = (TH1D*)f_Root->Get(var.c_str());
    h_neutron->SetFillColor(kBlue);
    h_neutron->SetMarkerColor(kBlue);

    var = type + "invMass_other";
    TH1D* h_other = (TH1D*)f_Root->Get(var.c_str());
    h_other->SetFillColor(kGray);
    h_other->SetMarkerColor(kGray);


    hs->Add(h_pizero);
    hs->Add(h_piplus);
    hs->Add(h_proton);
    hs->Add(h_neutron);
    hs->Add(h_other);
    if (isStacked) hs->Draw();
    else hs->Draw("nostack");
    hs->GetXaxis()->SetTitle("m_{#gamma#gamma} [MeV]");
    hs->GetYaxis()->SetTitle("N(Events)");

    // Add Legend
    TLegend *legend = new TLegend(0.7,0.68,0.9,0.9);  
    legend->AddEntry(h_pizero, "#pi^{0}", "f");
    legend->AddEntry(h_piplus, "#pi^{+}", "f");
    legend->AddEntry(h_proton, "proton", "f");
    legend->AddEntry(h_neutron, "neutron", "f");
    legend->AddEntry(h_other, "other", "f");
    legend->SetTextSize(0.05);
    legend->Draw();

    // Add Pi0 InvMass Lines
    double hist_max;
    if (isStacked) hist_max = hs->GetMaximum();
    else hist_max = h_pizero->GetMaximum();
    TLine pi0Mass;
    pi0Mass.SetLineWidth(2);
    pi0Mass.SetLineColor(kBlue);
    pi0Mass.DrawLine(134.98,0,134.98,hist_max);

    TLine pi0Mass_min;
    pi0Mass_min.SetLineWidth(2);
    pi0Mass_min.SetLineColor(kBlack);
    pi0Mass_min.DrawLine(60.0,0,60.0,hist_max);

    TLine pi0Mass_max;
    pi0Mass_max.SetLineWidth(2);
    pi0Mass_max.SetLineColor(kBlack);
    pi0Mass_max.DrawLine(200.0,0,200.0,hist_max);

  
    std::string plot_type;
    std::string out_name;
    if (isStacked) plot_type = "_Stacked";
    else plot_type = "";
    out_name = plotDir + type + "invMass_TruthMatch" + plot_type + ".png"; 

    c1->Print(out_name.c_str(),"png");

    delete f_Root;
    delete hs;
    delete legend;
    delete c1;
}

void CCProtonPi0_Plotter::plot_Michel_TruthMatch(std::string var)
{
    // mc_w written during reduce - Its Histogram is with Cut Hists
    std::string root_dir = rootDir_CutHists.mc;
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    std::cout<<"\nPlottting Stacked Michel Truth Match"<<std::endl;

    TFile* f_Root = new TFile(root_dir.c_str());
    TCanvas* c1 = new TCanvas("c","c",1280,800);
    THStack *hs = new THStack("hs","Michel Electron Truth Match");
    
    std::string var_name = "michel_piplus_" + var;
    TH1D* h_piplus = (TH1D*)f_Root->Get(var_name.c_str());
    h_piplus->SetFillColor(kGreen);
    h_piplus->SetFillStyle(3001);

    var_name = "michel_neutron_" + var;
    TH1D* h_neutron = (TH1D*)f_Root->Get(var_name.c_str());
    h_neutron->SetFillColor(kRed);
    h_neutron->SetFillStyle(3001);

    var_name = "michel_proton_" + var;
    TH1D* h_proton = (TH1D*)f_Root->Get(var_name.c_str());
    h_proton->SetFillColor(kBlue);
    h_proton->SetFillStyle(3001);

    var_name = "michel_piminus_" + var;
    TH1D* h_piminus = (TH1D*)f_Root->Get(var_name.c_str());
    h_piminus->SetFillColor(kMagenta);
    h_piminus->SetFillStyle(3001);

    var_name = "michel_other_" + var;
    TH1D* h_other = (TH1D*)f_Root->Get(var_name.c_str());
    h_other->SetFillColor(kGray);
    h_other->SetFillStyle(3001);

    // Add Legend
    TLegend *legend = new TLegend(0.7,0.68,0.9,0.9);  
    legend->AddEntry(h_other, "other", "f");
    legend->AddEntry(h_piminus, "#pi^{-}", "f");
    legend->AddEntry(h_proton, "proton", "f");
    legend->AddEntry(h_neutron, "neutron", "f");
    legend->AddEntry(h_piplus, "#pi^{+}", "f");

    hs->Add(h_piplus);
    hs->Add(h_neutron);
    hs->Add(h_proton);
    hs->Add(h_piminus);
    hs->Add(h_other);
    hs->Draw();
  
    std::cout<<h_piplus->GetXaxis()->GetTitle()<<std::endl;
    std::cout<<h_piplus->GetYaxis()->GetTitle()<<std::endl;
    
    hs->GetXaxis()->SetTitle(h_piplus->GetXaxis()->GetTitle());
    hs->GetYaxis()->SetTitle(h_piplus->GetYaxis()->GetTitle());

    legend->Draw();

    std::string plot_name = "TruthMatch_michel_" + var_name + ".png";
    c1->Print(Form("%s%s",plotDir.c_str(),plot_name.c_str()),"png");

    delete f_Root;
    delete hs;
    delete legend;
    delete c1;
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


void CCProtonPi0_Plotter::setRootDirs()
{
    // Set Other Studies ROOT Dir
    rootDir_Truth.mc = Folder_List::rootDir_Truth_mc;
    rootDir_Truth.data = Folder_List::rootDir_Truth_data;

    rootDir_OtherStudies.mc = Folder_List::rootDir_OtherStudies_mc; 
    rootDir_OtherStudies.data = Folder_List::rootDir_OtherStudies_data;

    rootDir_CrossSection.mc = Folder_List::rootDir_CrossSection_mc;
    rootDir_CrossSection.data = Folder_List::rootDir_CrossSection_data;

    rootDir_GENIEXSec.mc = Folder_List::rootDir_GENIEXSec;
    rootDir_GENIEXSec.data = Folder_List::rootDir_GENIEXSec;
    
    // Set MC Root Dir
    rootDir_CutHists.mc = Folder_List::rootDir_CutHists_mc;
    rootDir_Interaction.mc = Folder_List::rootDir_Interaction_mc;
    rootDir_Muon.mc = Folder_List::rootDir_Muon_mc;
    rootDir_Proton.mc = Folder_List::rootDir_Proton_mc;
    rootDir_Pion.mc = Folder_List::rootDir_Pion_mc;

    // Set Data Root Dir
    rootDir_CutHists.data = Folder_List::rootDir_CutHists_data;
    rootDir_Interaction.data = Folder_List::rootDir_Interaction_data;
    rootDir_Muon.data = Folder_List::rootDir_Muon_data;
    rootDir_Proton.data = Folder_List::rootDir_Proton_data;
    rootDir_Pion.data = Folder_List::rootDir_Pion_data;
}


void CCProtonPi0_Plotter::plotGENIEXSec()
{
    std::cout<<"Plotting GENIE Cross Sections..."<<std::endl;
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    // Plot 1D Versions
    DrawMnvH1D(rootDir_GENIEXSec,"muon_P_xsec",plotDir);   
    DrawMnvH1D(rootDir_GENIEXSec,"muon_theta_xsec",plotDir);   
    DrawMnvH1D(rootDir_GENIEXSec,"pi0_P_xsec",plotDir);   
    DrawMnvH1D(rootDir_GENIEXSec,"pi0_KE_xsec",plotDir);   
    DrawMnvH1D(rootDir_GENIEXSec,"pi0_theta_xsec",plotDir);   
    DrawMnvH1D(rootDir_GENIEXSec,"QSq_xsec",plotDir);   

    // Plot Comparison
    TFile* f_xsec_mc = new TFile(rootDir_GENIEXSec.mc.c_str());
    TFile* f_xsec_data = new TFile(rootDir_CrossSection.mc.c_str());
    MnvH1D* mc;
    MnvH1D* data;

    mc = (MnvH1D*)f_xsec_mc->Get("muon_P_xsec");
    data = (MnvH1D*)f_xsec_data->Get("muon_P_xsec"); 
    DrawDataMC_CrossSection(data,mc,"muon_P_MC_GENIE",plotDir);
 
    mc = (MnvH1D*)f_xsec_mc->Get("muon_theta_xsec");
    data = (MnvH1D*)f_xsec_data->Get("muon_theta_xsec"); 
    DrawDataMC_CrossSection(data,mc,"muon_theta_MC_GENIE",plotDir);
  
    mc = (MnvH1D*)f_xsec_mc->Get("pi0_P_xsec");
    data = (MnvH1D*)f_xsec_data->Get("pi0_P_xsec"); 
    DrawDataMC_CrossSection(data,mc,"pi0_P_MC_GENIE",plotDir);
   
    mc = (MnvH1D*)f_xsec_mc->Get("pi0_KE_xsec");
    data = (MnvH1D*)f_xsec_data->Get("pi0_KE_xsec"); 
    DrawDataMC_CrossSection(data,mc,"pi0_KE_MC_GENIE",plotDir);
  
    mc = (MnvH1D*)f_xsec_mc->Get("pi0_theta_xsec");
    data = (MnvH1D*)f_xsec_data->Get("pi0_theta_xsec"); 
    DrawDataMC_CrossSection(data,mc,"pi0_theta_MC_GENIE",plotDir);
 
    mc = (MnvH1D*)f_xsec_mc->Get("QSq_xsec");
    data = (MnvH1D*)f_xsec_data->Get("QSq_xsec"); 
    DrawDataMC_CrossSection(data,mc,"QSq_MC_GENIE",plotDir);

    std::cout<<"Done!"<<std::endl;
}
#endif

