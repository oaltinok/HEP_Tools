/*
 * See CCProtonPi0_Plotter.h header for Class Information
 */

#ifndef CCProtonPi0_Plotter_cpp
#define CCProtonPi0_Plotter_cpp

#include "CCProtonPi0_Plotter.h"

using namespace PlotUtils;

void CCProtonPi0_Plotter::plotHistograms()
{
    thesisStyle = true;

    if (thesisStyle){
        ApplyStyle_Thesis();
    }

    //--------------------------------------------------------------------------
    // Run only once to get the POT
    //--------------------------------------------------------------------------
    //getPOT_Data();
    //getPOT_MC();

    //--------------------------------------------------------------------------
    // Cross Sections
    //--------------------------------------------------------------------------
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
    //plotTruth();
    //plotInteraction_MCOnly();
    //plotMuon_MCOnly();
    //plotProton_MCOnly();
    //plotPion_MCOnly();
    //plotCutHistograms_MCOnly();

    //--------------------------------------------------------------------------
    //  Plot Function Reserved for Other Studies
    //--------------------------------------------------------------------------
    //BckgSubtraction_Studies();
    //SavePi0InvMassPoints();
    //plotOtherStudies();
    //plotGENIEXSec();
    //UnfoldingStudy();
    //GENIE_Tuning_Study();
    Systematics();
    //W_Studies();
    //QSq_Studies();
    //Studies_2p2h();
    //DeltaRes_Studies();
    //PlotFluxHistograms();
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
    Systematics_SetErrorSummaryGroups();
}

void CCProtonPi0_Plotter::CheckAllUniverses(std::string test_name, MnvH1D* data, MnvH1D* mc)
{
    
    // Check CV First
    int nBins = data->GetNbinsX();
    for (int i = 0; i < nBins; ++i){
        double data_content = data->GetBinContent(i);
        double mc_content = mc->GetBinContent(i);
        double ratio = data_content/mc_content;
        if (isnan(ratio) || mc_content < EPSILON || data_content < EPSILON) continue;
        else if ( fabs(ratio - 1.0) > EPSILON){
            std::cout<<"-----------------------------------------------------------"<<std::endl;
            std::cout<<"\tClosure Test Failed for "<<test_name<<std::endl;
            std::cout<<"\tCentral Value"<<std::endl;
            std::cout<<"\tBin = "<<i<<std::endl;
            std::cout<<"\tData = "<<data_content<<" MC = "<<mc_content<<std::endl;
            std::cout<<"\tData/MC = "<<data_content/mc_content<<std::endl;
            //exit(1);
        }
    }

    // Vertical Error Bands
    std::vector<std::string> data_vert_errors = data->GetVertErrorBandNames();
    std::vector<std::string> mc_vert_errors = mc->GetVertErrorBandNames();

    for (unsigned int i = 0; i < data_vert_errors.size(); ++i){
        const MnvVertErrorBand* data_err_band = data->GetVertErrorBand(data_vert_errors[i]);
        const MnvVertErrorBand* mc_err_band = mc->GetVertErrorBand(mc_vert_errors[i]);

        const std::vector<TH1D*> data_hists = data_err_band->GetHists();
        const std::vector<TH1D*> mc_hists = mc_err_band->GetHists();
        
        CompareUniversesBinByBin(data_hists, mc_hists, data_vert_errors[i], test_name);
    }

    // Lateral Error Bands
    std::vector<std::string> data_lat_errors = data->GetLatErrorBandNames();
    std::vector<std::string> mc_lat_errors = mc->GetLatErrorBandNames();

    for (unsigned int i = 0; i < data_lat_errors.size(); ++i){
        const MnvLatErrorBand* data_err_band = data->GetLatErrorBand(data_lat_errors[i]);
        const MnvLatErrorBand* mc_err_band = mc->GetLatErrorBand(mc_lat_errors[i]);

        const std::vector<TH1D*> data_hists = data_err_band->GetHists();
        const std::vector<TH1D*> mc_hists = mc_err_band->GetHists();
        
        CompareUniversesBinByBin(data_hists, mc_hists, data_lat_errors[i], test_name);
    }
}

void CCProtonPi0_Plotter::CompareUniversesBinByBin(const std::vector<TH1D*> data_hists, const std::vector<TH1D*> mc_hists, std::string err_name, std::string test_name)
{
    for (unsigned int j = 0; j < data_hists.size(); ++j){
        int nBins = data_hists[j]->GetNbinsX();
        for (int bin = 1; bin <= nBins; ++bin){
            double data_bin_content = data_hists[j]->GetBinContent(bin);
            double mc_bin_content = mc_hists[j]->GetBinContent(bin);
            double ratio = data_bin_content/mc_bin_content;
            if (isnan(ratio) || mc_bin_content == 0 || data_bin_content == 0) continue;
            else if ( fabs(ratio - 1.0) > EPSILON){
                std::cout<<"-----------------------------------------------------------"<<std::endl;
                std::cout<<"\tClosure Test Failed for "<<test_name<<std::endl;
                std::cout<<"\tError Band = "<<err_name<<std::endl;
                std::cout<<"\tUniverse = "<<j<<std::endl;
                std::cout<<"\tBin = "<<bin<<std::endl;
                std::cout<<"\tData = "<<data_bin_content<<" MC = "<<mc_bin_content<<std::endl;
                std::cout<<"\tData/MC = "<<ratio<<std::endl;
            }
        }
    }
}

void CCProtonPi0_Plotter::plotCrossSection_Check(std::string var_name, std::string plotDir)
{
    std::cout<<"Plotting Cross Section Check for "<<var_name<<std::endl;
   
    std::string var;
    
    TFile* f_xsec;
    MnvH1D* data;
    MnvH1D* mc;

    f_xsec = new TFile(rootDir_CrossSection.mc.c_str());

    // Estimated Background 
    var = var_name + "_mc_reco_bckg";
    mc = (MnvH1D*)f_xsec->Get(var.c_str());
    var = var_name + "_bckg_estimated";
    data = (MnvH1D*)f_xsec->Get(var.c_str()); 
    var = var_name + "_check_bckg_estimated";
    DrawDataMC_Thesis(data, mc, var, plotDir, true);
    CheckAllUniverses("Estimated Background",data, mc);

    // Subtracted Background 
    var = var_name + "_mc_reco_signal";
    mc = (MnvH1D*)f_xsec->Get(var.c_str());
    var = var_name + "_bckg_subtracted";
    data = (MnvH1D*)f_xsec->Get(var.c_str()); 
    var = var_name + "_check_bckg_subtracted";
    DrawDataMC_Thesis(data, mc, var, plotDir, true);
    CheckAllUniverses("Subtracted Background",data, mc);

    // Unfolded 
    var = var_name + "_mc_truth_signal";
    mc = (MnvH1D*)f_xsec->Get(var.c_str());
    var = var_name + "_unfolded";
    data = (MnvH1D*)f_xsec->Get(var.c_str()); 
    var = var_name + "_check_unfolded";
    DrawDataMC_Thesis(data, mc, var, plotDir, true);
    CheckAllUniverses("Unfolded",data, mc);
   
    // Efficiency Correction 
    var = var_name + "_mc_truth_all_signal";
    mc = (MnvH1D*)f_xsec->Get(var.c_str());
    var = var_name + "_efficiency_corrected";
    data = (MnvH1D*)f_xsec->Get(var.c_str()); 
    var = var_name + "_check_efficiency_corrected";
    DrawDataMC_Thesis(data, mc, var, plotDir, true);
    CheckAllUniverses("Efficiency Corrected",data, mc);
    
    std::cout<<"Plotting Cross Section Check Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plotCrossSection_Check()
{
    std::string plotDir;

    plotDir = Folder_List::xsec_muon_P + Folder_List::plotDir_Check;
    plotCrossSection_Check("muon_P", plotDir);

    plotDir = Folder_List::xsec_muon_theta + Folder_List::plotDir_Check;
    plotCrossSection_Check("muon_theta", plotDir);

    plotDir = Folder_List::xsec_pi0_P + Folder_List::plotDir_Check;
    plotCrossSection_Check("pi0_P", plotDir);

    plotDir = Folder_List::xsec_pi0_KE + Folder_List::plotDir_Check;
    plotCrossSection_Check("pi0_KE", plotDir);

    plotDir = Folder_List::xsec_pi0_theta + Folder_List::plotDir_Check;
    plotCrossSection_Check("pi0_theta", plotDir);

    plotDir = Folder_List::xsec_QSq + Folder_List::plotDir_Check;
    plotCrossSection_Check("QSq", plotDir);

    plotDir = Folder_List::xsec_W + Folder_List::plotDir_Check;
    plotCrossSection_Check("W", plotDir);

    plotDir = Folder_List::xsec_Enu + Folder_List::plotDir_Check;
    plotCrossSection_Check("Enu", plotDir);
}

void CCProtonPi0_Plotter::plotOtherStudies()
{
    std::cout<<"Plotting Other Studies..."<<std::endl;
    std::string plotDir = Folder_List::plotDir_OtherStudies;
   
    //DrawDataStackedMC_WithSignalTypes(rootDir_CutHists, "SideBand_QSq", plotDir);
    //DrawDataStackedMC(rootDir_CutHists, "SideBand_QSq", plotDir);
   
    TFile* f = new TFile(rootDir_Interaction.mc.c_str());

    MnvH1D* mc =  GetMnvH1D(f, "QSq_MaRES_0");

    std::vector<std::string> vert_errors = mc->GetVertErrorBandNames();
    
    for (unsigned int i = 0; i < vert_errors.size(); ++i){
        std::cout<<vert_errors[i]<<std::endl;
    }


    std::cout<<"Plotting Other Studies Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plotPC_MINOS_Pi0()
{
    std::string rootDir_Steel_1 = "/minerva/data/users/oaltinok/NTupleAnalysis/ParticleCannon/PC_Steel_1.root"; 
    std::string rootDir_Steel_2 = "/minerva/data/users/oaltinok/NTupleAnalysis/ParticleCannon/PC_Steel_2.root"; 
    std::string rootDir_Carbon_1 = "/minerva/data/users/oaltinok/NTupleAnalysis/ParticleCannon/PC_Carbon_1.root"; 
    std::string rootDir_Carbon_2 = "/minerva/data/users/oaltinok/NTupleAnalysis/ParticleCannon/PC_Carbon_2.root"; 
 
    plotPC_MINOS_Pi0(rootDir_Steel_1, rootDir_Carbon_1, "plane_energy", "Energy_1_1");
    plotPC_MINOS_Pi0(rootDir_Steel_1, rootDir_Carbon_1, "nPlanes", "Energy_1_1");
    plotPC_MINOS_Pi0(rootDir_Steel_1, rootDir_Carbon_1, "plane_z", "Energy_1_1");
 
    plotPC_MINOS_Pi0(rootDir_Steel_2, rootDir_Carbon_2, "plane_energy", "Energy_2_0");
    plotPC_MINOS_Pi0(rootDir_Steel_2, rootDir_Carbon_2, "nPlanes", "Energy_2_0");
    plotPC_MINOS_Pi0(rootDir_Steel_2, rootDir_Carbon_2, "plane_z", "Energy_2_0");

}

void CCProtonPi0_Plotter::plotPC_MINOS_Pi0(std::string rootDir_Steel, std::string rootDir_Carbon, std::string var, std::string label)
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;
    
    TFile* f_Root_Steel = new TFile(rootDir_Steel.c_str());
    TFile* f_Root_Carbon = new TFile(rootDir_Carbon.c_str());
    
    // Create Canvas
    TCanvas* c1 = new TCanvas("c","c",1280,800);
    THStack *hs = new THStack("hs",var.c_str());

    TH1D* h_steel = (TH1D*)f_Root_Steel->Get(var.c_str());
    h_steel->SetFillColor(kRed);
    h_steel->SetLineColor(kRed);
    h_steel->SetLineWidth(3);
    h_steel->SetFillStyle(3010);

    TH1D* h_carbon = (TH1D*)f_Root_Carbon->Get(var.c_str());
    h_carbon->SetFillColor(kBlue);
    h_carbon->SetLineColor(kBlue);
    h_carbon->SetLineWidth(3);
    h_carbon->SetFillStyle(3010);

    hs->Add(h_steel);
    hs->Add(h_carbon);
    
    hs->Draw("nostack");
    hs->GetXaxis()->SetTitle(var.c_str());
    hs->GetYaxis()->SetTitle("N(Events)");

    // Add Legend
    TLegend *legend = new TLegend(0.7,0.68,0.9,0.9);  
    legend->AddEntry(h_steel, "Steel", "f");
    legend->AddEntry(h_carbon, "Carbon", "f");
    legend->SetTextSize(0.05);
    legend->Draw();

    std::string out_name = plotDir + var + "_" + label + ".png"; 

    c1->Print(out_name.c_str(),"png");

    delete f_Root_Steel;
    delete f_Root_Carbon;
    delete hs;
    delete legend;
    delete c1;
    
    std::cout<<"Plotting Other Studies Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::PlotTotalEnuXSec()
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    // ROOT Files
    std::string rootDir_data_All = Folder_List::rootOut + Folder_List::Data + Folder_List::analyzed + "CrossSection_v2_81_All.root";
    std::string rootDir_MC_All = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + "CrossSection_v2_81_All.root";
    std::string rootDir_data_W1800 = Folder_List::rootOut + Folder_List::Data + Folder_List::analyzed + "CrossSection_v2_81_W1800.root";
    std::string rootDir_MC_W1800 = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + "CrossSection_v2_81_W1800.root";
    std::string rootDir_data_W1400 = Folder_List::rootOut + Folder_List::Data + Folder_List::analyzed + "CrossSection_v2_81_W1400.root";
    std::string rootDir_MC_W1400 = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + "CrossSection_v2_81_W1400.root";

    // Get Histograms
    TFile* f_data_All = new TFile(rootDir_data_All.c_str());
    TFile* f_MC_All = new TFile(rootDir_MC_All.c_str());

    TFile* f_data_W1800 = new TFile(rootDir_data_W1800.c_str());
    TFile* f_MC_W1800 = new TFile(rootDir_MC_W1800.c_str());

    TFile* f_data_W1400 = new TFile(rootDir_data_W1400.c_str());
    TFile* f_MC_W1400 = new TFile(rootDir_MC_W1400.c_str());

    MnvH1D* data_All = GetMnvH1D(f_data_All, "Enu_xsec");
    MnvH1D* data_W1800 = GetMnvH1D(f_data_W1800, "Enu_xsec");
    MnvH1D* data_W1400 = GetMnvH1D(f_data_W1400, "Enu_xsec");
    MnvH1D* MC_All = GetMnvH1D(f_MC_All, "Enu_xsec");
    MnvH1D* MC_W1800 = GetMnvH1D(f_MC_W1800, "Enu_xsec");
    MnvH1D* MC_W1400 = GetMnvH1D(f_MC_W1400, "Enu_xsec");
   
    // Normalize Histograms
    NormalizeToNormBinWidth(data_All);
    NormalizeToNormBinWidth(MC_All);
    NormalizeToNormBinWidth(data_W1800);
    NormalizeToNormBinWidth(MC_W1800);
    NormalizeToNormBinWidth(data_W1400);
    NormalizeToNormBinWidth(MC_W1400);

    // Style Histograms
    // Data
    data_All->SetMinimum(0.0);
    data_All->SetMaximum(80.0);
    data_All->SetMarkerColor(kBlack);
    data_All->SetMarkerStyle(20);
    data_All->SetMarkerSize(1);
    data_All->SetLineWidth(1);
    data_All->SetLineColor(kBlack);
    data_All->SetFillStyle(0);

    data_W1800->SetMarkerColor(kBlue);
    data_W1800->SetMarkerStyle(20);
    data_W1800->SetMarkerSize(1);
    data_W1800->SetLineWidth(1);
    data_W1800->SetLineColor(kBlue);
    data_W1800->SetFillStyle(0);

    data_W1400->SetMarkerColor(kGreen+3);
    data_W1400->SetMarkerStyle(20);
    data_W1400->SetMarkerSize(1);
    data_W1400->SetLineWidth(1);
    data_W1400->SetLineColor(kGreen+3);
    data_W1400->SetFillStyle(0);

    // MC
    MC_All->SetLineWidth(2);
    MC_All->SetLineColor(kRed);
    MC_All->SetFillStyle(0);

    MC_W1800->SetLineWidth(2);
    MC_W1800->SetLineColor(kMagenta);
    MC_W1800->SetFillStyle(0);

    MC_W1400->SetLineWidth(2);
    MC_W1400->SetLineColor(kTeal);
    MC_W1400->SetFillStyle(0);

    // Plot
    TCanvas* c1 = new TCanvas("c","c",1280,800);

    data_All->Draw("E1 X0");
    MC_All->Draw("SAME HIST");
    
    data_W1800->Draw("SAME E1 X0");
    MC_W1800->Draw("SAME HIST");
    
    data_W1400->Draw("SAME E1 X0");
    MC_W1400->Draw("SAME HIST");

    // Add Legend
    TLegend *legend = new TLegend(0.7,0.68,0.9,0.9);  
    legend->AddEntry(data_All, "Data All");
    legend->AddEntry(MC_All, "MC All");
    legend->AddEntry(data_W1800, "Data (W < 1.8 GeV)");
    legend->AddEntry(MC_W1800, "MC (W < 1.8 GeV)");
    legend->AddEntry(data_W1400, "Data (W < 1.4 GeV)");
    legend->AddEntry(MC_W1400, "MC (W < 1.4 GeV)");
    legend->SetTextSize(0.03);
    legend->Draw();

    // Plot Output
    gStyle->SetOptStat(0); 
    c1->Update();
    std::string out_name;
    out_name = plotDir + "TotalEnuXSec.png"; 

    c1->Print(out_name.c_str(),"png");

    delete legend;
    delete c1;
}


void CCProtonPi0_Plotter::plotOriginalData()
{
    std::cout<<"Plotting Original Data..."<<std::endl;
    std::string plotDir;
    
    if (plot_muon_P){
        plotDir = Folder_List::xsec_muon_P + Folder_List::plotDir_Original;
        PlotXSecVar("muon_P", "all", "mc_reco_all", plotDir, "data_MC" );
    }
  
    if (plot_muon_theta){
        plotDir = Folder_List::xsec_muon_theta + Folder_List::plotDir_Original;
        PlotXSecVar("muon_theta", "all", "mc_reco_all", plotDir, "data_MC" );
    }
  
    if (plot_pi0_P){
        plotDir = Folder_List::xsec_pi0_P + Folder_List::plotDir_Original;
        PlotXSecVar("pi0_P", "all", "mc_reco_all", plotDir, "data_MC" );
    }

    if (plot_pi0_KE){
        plotDir = Folder_List::xsec_pi0_KE + Folder_List::plotDir_Original;
        PlotXSecVar("pi0_KE", "all", "mc_reco_all", plotDir, "data_MC" );
    }
 
    if (plot_pi0_theta){
        plotDir = Folder_List::xsec_pi0_theta + Folder_List::plotDir_Original;
        PlotXSecVar("pi0_theta", "all", "mc_reco_all", plotDir, "data_MC" );
    }

    if (plot_QSq){
        plotDir = Folder_List::xsec_QSq + Folder_List::plotDir_Original;
        PlotXSecVar("QSq", "all", "mc_reco_all", plotDir, "data_MC" );
    }

    if (plot_W){
        plotDir = Folder_List::xsec_W + Folder_List::plotDir_Original;
        PlotXSecVar("W", "all", "mc_reco_all", plotDir, "data_MC" );
    }

    if (plot_Enu){
        plotDir = Folder_List::xsec_Enu + Folder_List::plotDir_Original;
        PlotXSecVar("Enu", "all", "mc_reco_all", plotDir, "data_MC" );
    }
    
    std::cout<<"Plotting Original Data Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plotBackgroundEstimated()
{
    std::cout<<"Plotting Background Estimated..."<<std::endl;
    std::string plotDir;
    
    if (plot_muon_P){
        plotDir = Folder_List::xsec_muon_P + Folder_List::plotDir_BackgroundEstimated;
        PlotXSecVar("muon_P", "bckg_estimated", "bckg_estimated", plotDir, "bckg_estimated_data_MC" );
    }
  
    if (plot_muon_theta){
        plotDir = Folder_List::xsec_muon_theta + Folder_List::plotDir_BackgroundEstimated;
        PlotXSecVar("muon_theta", "bckg_estimated", "bckg_estimated", plotDir, "bckg_estimated_data_MC" );
    }
  
    if (plot_pi0_P){
        plotDir = Folder_List::xsec_pi0_P + Folder_List::plotDir_BackgroundEstimated;
        PlotXSecVar("pi0_P", "bckg_estimated", "bckg_estimated", plotDir, "bckg_estimated_data_MC" );
    }

    if (plot_pi0_KE){
        plotDir = Folder_List::xsec_pi0_KE + Folder_List::plotDir_BackgroundEstimated;
        PlotXSecVar("pi0_KE", "bckg_estimated", "bckg_estimated", plotDir, "bckg_estimated_data_MC" );
    }
 
    if (plot_pi0_theta){
        plotDir = Folder_List::xsec_pi0_theta + Folder_List::plotDir_BackgroundEstimated;
        PlotXSecVar("pi0_theta", "bckg_estimated", "bckg_estimated", plotDir, "bckg_estimated_data_MC" );
    }

    if (plot_QSq){
        plotDir = Folder_List::xsec_QSq + Folder_List::plotDir_BackgroundEstimated;
        PlotXSecVar("QSq", "bckg_estimated", "bckg_estimated", plotDir, "bckg_estimated_data_MC" );
    }

    if (plot_W){
        plotDir = Folder_List::xsec_W + Folder_List::plotDir_BackgroundEstimated;
        PlotXSecVar("W", "bckg_estimated", "bckg_estimated", plotDir, "bckg_estimated_data_MC" );
    }

    if (plot_Enu){
        plotDir = Folder_List::xsec_Enu + Folder_List::plotDir_BackgroundEstimated;
        PlotXSecVar("Enu", "bckg_estimated", "bckg_estimated", plotDir, "bckg_estimated_data_MC" );
    }

    std::cout<<"Plotting Background Estimated Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plotBackgroundSubtracted()
{
    std::cout<<"Plotting Background Subtracted..."<<std::endl;
    std::string plotDir;
 
    if (plot_muon_P){
        plotDir = Folder_List::xsec_muon_P + Folder_List::plotDir_BackgroundSubtracted;
        PlotXSecVar("muon_P", "bckg_subtracted", "bckg_subtracted", plotDir, "bckg_subtracted_data_MC" );
        PlotXSecVar("invMass", "all", "mc_reco_all", plotDir, "invMass_data_MC" );
    }
  
    if (plot_muon_theta){
        plotDir = Folder_List::xsec_muon_theta + Folder_List::plotDir_BackgroundSubtracted;
        PlotXSecVar("muon_theta", "bckg_subtracted", "bckg_subtracted", plotDir, "bckg_subtracted_data_MC" );
        PlotXSecVar("invMass", "all", "mc_reco_all", plotDir, "invMass_data_MC" );
    }
  
    if (plot_pi0_P){
        plotDir = Folder_List::xsec_pi0_P + Folder_List::plotDir_BackgroundSubtracted;
        PlotXSecVar("pi0_P", "bckg_subtracted", "bckg_subtracted", plotDir, "bckg_subtracted_data_MC" );
        PlotXSecVar("invMass", "all", "mc_reco_all", plotDir, "invMass_data_MC" );
    }

    if (plot_pi0_KE){
        plotDir = Folder_List::xsec_pi0_KE + Folder_List::plotDir_BackgroundSubtracted;
        PlotXSecVar("pi0_KE", "bckg_subtracted", "bckg_subtracted", plotDir, "bckg_subtracted_data_MC" );
        PlotXSecVar("invMass", "all", "mc_reco_all", plotDir, "invMass_data_MC" );
    }
 
    if (plot_pi0_theta){
        plotDir = Folder_List::xsec_pi0_theta + Folder_List::plotDir_BackgroundSubtracted;
        PlotXSecVar("pi0_theta", "bckg_subtracted", "bckg_subtracted", plotDir, "bckg_subtracted_data_MC" );
        PlotXSecVar("invMass", "all", "mc_reco_all", plotDir, "invMass_data_MC" );
    }

    if (plot_QSq){
        plotDir = Folder_List::xsec_QSq + Folder_List::plotDir_BackgroundSubtracted;
        PlotXSecVar("QSq", "bckg_subtracted", "bckg_subtracted", plotDir, "bckg_subtracted_data_MC" );
        PlotXSecVar("invMass", "all", "mc_reco_all", plotDir, "invMass_data_MC" );
    }

    if (plot_W){
        plotDir = Folder_List::xsec_W + Folder_List::plotDir_BackgroundSubtracted;
        PlotXSecVar("W", "bckg_subtracted", "bckg_subtracted", plotDir, "bckg_subtracted_data_MC" );
        PlotXSecVar("invMass", "all", "mc_reco_all", plotDir, "invMass_data_MC" );
    }

    if (plot_Enu){
        plotDir = Folder_List::xsec_Enu + Folder_List::plotDir_BackgroundSubtracted;
        PlotXSecVar("Enu", "bckg_subtracted", "bckg_subtracted", plotDir, "bckg_subtracted_data_MC" );
        PlotXSecVar("invMass", "all", "mc_reco_all", plotDir, "invMass_data_MC" );
    }

    std::cout<<"Plotting Background Subtracted Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plotUnfolded()
{
    std::cout<<"Plotting Unfolded..."<<std::endl;
    std::string plotDir;
 
    if (plot_muon_P){
        plotDir = Folder_List::xsec_muon_P + Folder_List::plotDir_Unfolded;
        PlotXSecVar("muon_P", "unfolded", "unfolded", plotDir, "unfolded_data_MC" );
        DrawNormalizedMigrationHistogram(rootDir_CrossSection,"muon_P_response",plotDir);
    }
  
    if (plot_muon_theta){
        plotDir = Folder_List::xsec_muon_theta + Folder_List::plotDir_Unfolded;
        PlotXSecVar("muon_theta", "unfolded", "unfolded", plotDir, "unfolded_data_MC" );
        DrawNormalizedMigrationHistogram(rootDir_CrossSection,"muon_theta_response",plotDir);
    }
  
    if (plot_pi0_P){
        plotDir = Folder_List::xsec_pi0_P + Folder_List::plotDir_Unfolded;
        PlotXSecVar("pi0_P", "unfolded", "unfolded", plotDir, "unfolded_data_MC" );
        DrawNormalizedMigrationHistogram(rootDir_CrossSection,"pi0_P_response",plotDir);
    }

    if (plot_pi0_KE){
        plotDir = Folder_List::xsec_pi0_KE + Folder_List::plotDir_Unfolded;
        PlotXSecVar("pi0_KE", "unfolded", "unfolded", plotDir, "unfolded_data_MC" );
        DrawNormalizedMigrationHistogram(rootDir_CrossSection,"pi0_KE_response",plotDir);
    }
 
    if (plot_pi0_theta){
        plotDir = Folder_List::xsec_pi0_theta + Folder_List::plotDir_Unfolded;
        PlotXSecVar("pi0_theta", "unfolded", "unfolded", plotDir, "unfolded_data_MC" );
        DrawNormalizedMigrationHistogram(rootDir_CrossSection,"pi0_theta_response",plotDir);
    }

    if (plot_QSq){
        plotDir = Folder_List::xsec_QSq + Folder_List::plotDir_Unfolded;
        PlotXSecVar("QSq", "unfolded", "unfolded", plotDir, "unfolded_data_MC" );
        DrawNormalizedMigrationHistogram(rootDir_CrossSection,"QSq_response",plotDir);
    }

    if (plot_W){
        plotDir = Folder_List::xsec_W + Folder_List::plotDir_Unfolded;
        PlotXSecVar("W", "unfolded", "unfolded", plotDir, "unfolded_data_MC" );
        DrawNormalizedMigrationHistogram(rootDir_CrossSection,"W_response",plotDir);
    }

    if (plot_Enu){
        plotDir = Folder_List::xsec_Enu + Folder_List::plotDir_Unfolded;
        PlotXSecVar("Enu", "unfolded", "unfolded", plotDir, "unfolded_data_MC" );
        DrawNormalizedMigrationHistogram(rootDir_CrossSection,"Enu_response",plotDir);
    }

    std::cout<<"Plotting Unfolded Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plotEfficiencyCorrected()
{
    std::cout<<"Plotting Efficiency Corrected..."<<std::endl;
    std::string plotDir;

    if (plot_muon_P){
        plotDir = Folder_List::xsec_muon_P + Folder_List::plotDir_Efficiency;
        PlotXSecVar("muon_P", "efficiency_corrected", "efficiency_corrected", plotDir, "efficiency_corrected_data_MC" );
        DrawEfficiencyCurve(rootDir_CrossSection,"muon_P_eff",plotDir);
        
        DrawMnvH1D(rootDir_CrossSection, "muon_P_mc_truth_all_signal", plotDir);
        DrawMnvH1D(rootDir_CrossSection, "muon_P_mc_truth_signal",plotDir);
    }
  
    if (plot_muon_theta){
        plotDir = Folder_List::xsec_muon_theta + Folder_List::plotDir_Efficiency;
        PlotXSecVar("muon_theta", "efficiency_corrected", "efficiency_corrected", plotDir, "efficiency_corrected_data_MC" );
        DrawEfficiencyCurve(rootDir_CrossSection,"muon_theta_eff",plotDir);
        DrawMnvH1D(rootDir_CrossSection, "muon_theta_mc_truth_all_signal", plotDir);
        DrawMnvH1D(rootDir_CrossSection, "muon_theta_mc_truth_signal",plotDir);
    }
  
    if (plot_pi0_P){
        plotDir = Folder_List::xsec_pi0_P + Folder_List::plotDir_Efficiency;
        PlotXSecVar("pi0_P", "efficiency_corrected", "efficiency_corrected", plotDir, "efficiency_corrected_data_MC" );
        DrawEfficiencyCurve(rootDir_CrossSection,"pi0_P_eff",plotDir);
        DrawMnvH1D(rootDir_CrossSection, "pi0_P_mc_truth_all_signal", plotDir);
        DrawMnvH1D(rootDir_CrossSection, "pi0_P_mc_truth_signal",plotDir);
    }

    if (plot_pi0_KE){
        plotDir = Folder_List::xsec_pi0_KE + Folder_List::plotDir_Efficiency;
        PlotXSecVar("pi0_KE", "efficiency_corrected", "efficiency_corrected", plotDir, "efficiency_corrected_data_MC" );
        DrawEfficiencyCurve(rootDir_CrossSection,"pi0_KE_eff",plotDir);
        DrawMnvH1D(rootDir_CrossSection, "pi0_KE_mc_truth_all_signal", plotDir);
        DrawMnvH1D(rootDir_CrossSection, "pi0_KE_mc_truth_signal",plotDir);
    }
 
    if (plot_pi0_theta){
        plotDir = Folder_List::xsec_pi0_theta + Folder_List::plotDir_Efficiency;
        PlotXSecVar("pi0_theta", "efficiency_corrected", "efficiency_corrected", plotDir, "efficiency_corrected_data_MC" );
        DrawEfficiencyCurve(rootDir_CrossSection,"pi0_theta_eff",plotDir);
        DrawMnvH1D(rootDir_CrossSection, "pi0_theta_mc_truth_all_signal", plotDir);
        DrawMnvH1D(rootDir_CrossSection, "pi0_theta_mc_truth_signal",plotDir);
    }

    if (plot_QSq){
        plotDir = Folder_List::xsec_QSq + Folder_List::plotDir_Efficiency;
        PlotXSecVar("QSq", "efficiency_corrected", "efficiency_corrected", plotDir, "efficiency_corrected_data_MC" );
        DrawEfficiencyCurve(rootDir_CrossSection,"QSq_eff",plotDir);
        DrawMnvH1D(rootDir_CrossSection, "QSq_mc_truth_all_signal", plotDir);
        DrawMnvH1D(rootDir_CrossSection, "QSq_mc_truth_signal",plotDir);
    }

    if (plot_W){
        plotDir = Folder_List::xsec_W + Folder_List::plotDir_Efficiency;
        PlotXSecVar("W", "efficiency_corrected", "efficiency_corrected", plotDir, "efficiency_corrected_data_MC" );
        DrawEfficiencyCurve(rootDir_CrossSection,"W_eff",plotDir);
        DrawMnvH1D(rootDir_CrossSection, "W_mc_truth_all_signal", plotDir);
        DrawMnvH1D(rootDir_CrossSection, "W_mc_truth_signal",plotDir);
    }

    if (plot_Enu){
        plotDir = Folder_List::xsec_Enu + Folder_List::plotDir_Efficiency;
        PlotXSecVar("Enu", "efficiency_corrected", "efficiency_corrected", plotDir, "efficiency_corrected_data_MC" );
        DrawEfficiencyCurve(rootDir_CrossSection,"Enu_eff",plotDir);
        DrawMnvH1D(rootDir_CrossSection, "Enu_mc_truth_all_signal", plotDir);
        DrawMnvH1D(rootDir_CrossSection, "Enu_mc_truth_signal",plotDir);
    }

    std::cout<<"Plotting Efficiency Corrected Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plotFluxIntegrated()
{
    std::cout<<"Plotting Integrated Flux..."<<std::endl;
    std::string plotDir;
 
    if (plot_muon_P){
        plotDir = Folder_List::xsec_muon_P + Folder_List::plotDir_FluxIntegrated;
        PlotXSecVar("muon_P", "flux_integrated", "flux_integrated", plotDir, "flux_integrated_data_MC" );
    }
  
    if (plot_muon_theta){
        plotDir = Folder_List::xsec_muon_theta + Folder_List::plotDir_FluxIntegrated;
        PlotXSecVar("muon_theta", "flux_integrated", "flux_integrated", plotDir, "flux_integrated_data_MC" );
    }
  
    if (plot_pi0_P){
        plotDir = Folder_List::xsec_pi0_P + Folder_List::plotDir_FluxIntegrated;
        PlotXSecVar("pi0_P", "flux_integrated", "flux_integrated", plotDir, "flux_integrated_data_MC" );
    }

    if (plot_pi0_KE){
        plotDir = Folder_List::xsec_pi0_KE + Folder_List::plotDir_FluxIntegrated;
        PlotXSecVar("pi0_KE", "flux_integrated", "flux_integrated", plotDir, "flux_integrated_data_MC" );
    }
 
    if (plot_pi0_theta){
        plotDir = Folder_List::xsec_pi0_theta + Folder_List::plotDir_FluxIntegrated;
        PlotXSecVar("pi0_theta", "flux_integrated", "flux_integrated", plotDir, "flux_integrated_data_MC" );
    }

    if (plot_QSq){
        plotDir = Folder_List::xsec_QSq + Folder_List::plotDir_FluxIntegrated;
        PlotXSecVar("QSq", "flux_integrated", "flux_integrated", plotDir, "flux_integrated_data_MC" );
    }

    if (plot_W){
        plotDir = Folder_List::xsec_W + Folder_List::plotDir_FluxIntegrated;
        PlotXSecVar("W", "flux_integrated", "flux_integrated", plotDir, "flux_integrated_data_MC" );
    }

    if (plot_Enu){
        plotDir = Folder_List::xsec_Enu + Folder_List::plotDir_FluxIntegrated;
        PlotXSecVar("Enu", "flux_integrated", "flux_integrated", plotDir, "flux_integrated_data_MC" );
        
        // Plot Used Flux Histogram
        TFile* f_xsec_mc = new TFile(rootDir_CrossSection.mc.c_str());
        MnvH1D* h_flux_original = GetMnvH1D(f_xsec_mc, "h_flux_minervaLE_FHC");
        MnvH1D* h_flux_rebinned = GetMnvH1D(f_xsec_mc, "h_flux_rebinned");
        PlotFluxRebinned(h_flux_original, h_flux_rebinned, plotDir);
        delete h_flux_original;
        delete h_flux_rebinned;
        delete f_xsec_mc;
    }

    std::cout<<"Plotting Flux Integrated Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plotXSec()
{
    std::cout<<"Plotting Cross Sections..."<<std::endl;
    std::string plotDir;
    
    if (plot_muon_P){
        plotDir = Folder_List::xsec_muon_P + Folder_List::plotDir_CrossSection;
        PlotXSecVar("muon_P", "xsec", "xsec", plotDir, "xsec_data_MC" );
        PlotXSecVar_BeforeFSI("muon_P", plotDir);
        PlotXSecVar_FSIType("muon_P", plotDir);
        PlotXSecVar_IntType("muon_P", plotDir);
    }
  
    if (plot_muon_theta){
        plotDir = Folder_List::xsec_muon_theta + Folder_List::plotDir_CrossSection;
        PlotXSecVar("muon_theta", "xsec", "xsec", plotDir, "xsec_data_MC" );
        PlotXSecVar_BeforeFSI("muon_theta", plotDir);
        PlotXSecVar_FSIType("muon_theta", plotDir);
        PlotXSecVar_IntType("muon_theta", plotDir);
    }
  
    if (plot_pi0_P){
        plotDir = Folder_List::xsec_pi0_P + Folder_List::plotDir_CrossSection;
        PlotXSecVar("pi0_P", "xsec", "xsec", plotDir, "xsec_data_MC" );
        PlotXSecVar_WithMiniBoone("pi0_P", plotDir);
        PlotXSecVar_BeforeFSI("pi0_P", plotDir);
        PlotXSecVar_FSIType("pi0_P", plotDir);
        PlotXSecVar_IntType("pi0_P", plotDir);
    }

    if (plot_pi0_KE){
        plotDir = Folder_List::xsec_pi0_KE + Folder_List::plotDir_CrossSection;
        PlotXSecVar("pi0_KE", "xsec", "xsec", plotDir, "xsec_data_MC" );
        PlotXSecVar_BeforeFSI("pi0_KE", plotDir);
        PlotXSecVar_FSIType("pi0_KE", plotDir);
        PlotXSecVar_IntType("pi0_KE", plotDir);
    }
 
    if (plot_pi0_theta){
        plotDir = Folder_List::xsec_pi0_theta + Folder_List::plotDir_CrossSection;
        PlotXSecVar("pi0_theta", "xsec", "xsec", plotDir, "xsec_data_MC" );
        PlotXSecVar_BeforeFSI("pi0_theta", plotDir);
        PlotXSecVar_FSIType("pi0_theta", plotDir);
        PlotXSecVar_IntType("pi0_theta", plotDir);
    }

    if (plot_QSq){
        plotDir = Folder_List::xsec_QSq + Folder_List::plotDir_CrossSection;
        PlotXSecVar("QSq", "xsec", "xsec", plotDir, "xsec_data_MC" );
        PlotXSecVar_BeforeFSI("QSq", plotDir);
        PlotXSecVar_FSIType("QSq", plotDir);
        PlotXSecVar_IntType("QSq", plotDir);
    }

    if (plot_W){
        plotDir = Folder_List::xsec_W + Folder_List::plotDir_CrossSection;
        PlotXSecVar("W", "xsec", "xsec", plotDir, "xsec_data_MC" );
        PlotXSecVar_BeforeFSI("W", plotDir);
        PlotXSecVar_FSIType("W", plotDir);
        PlotXSecVar_IntType("W", plotDir);
    }

    if (plot_Enu){
        plotDir = Folder_List::xsec_Enu + Folder_List::plotDir_CrossSection;
        PlotXSecVar("Enu", "xsec", "xsec", plotDir, "xsec_data_MC" );
        PlotXSecVar_BeforeFSI("Enu", plotDir);
        PlotXSecVar_FSIType("Enu", plotDir);
        PlotXSecVar_IntType("Enu", plotDir);
    }

    std::cout<<"Plotting Cross Sections Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plotTruth_Enu()
{
    std::string plotDir = Folder_List::plotDir_Interaction;
    //DrawNormalizedMigrationHistogram(rootDir_Interaction, "Enu_All_response", plotDir);
    //DrawNormalizedMigrationHistogram(rootDir_Interaction, "Enu_1Track_response", plotDir);
    //DrawNormalizedMigrationHistogram(rootDir_Interaction, "Enu_2Track_response", plotDir);
  
    //Draw1DHist(rootDir_Interaction,"Enu_Error",plotDir);
    //Draw1DHist(rootDir_Interaction,"Enu_1Track_Error",plotDir);
    //Draw1DHist(rootDir_Interaction,"Enu_2Track_Error",plotDir);
    
    Draw1DHist(rootDir_Interaction,"Enu_Diff",plotDir);
    Draw1DHist(rootDir_Interaction,"Enu_1Track_Diff",plotDir);
    Draw1DHist(rootDir_Interaction,"Enu_2Track_Diff",plotDir);

    //Draw1DHist(rootDir_Truth,"Enu_mc_truth_all_signal", plotDir);
    //Draw1DHist(rootDir_Truth,"QSq_mc_truth_all_signal", plotDir);
}

void CCProtonPi0_Plotter::plotTruth_QSq()
{
    std::string plotDir = Folder_List::plotDir_Interaction;
    //DrawNormalizedMigrationHistogram(rootDir_Interaction, "QSq_All_response", plotDir);
    //DrawNormalizedMigrationHistogram(rootDir_Interaction, "QSq_1Track_response", plotDir);
    //DrawNormalizedMigrationHistogram(rootDir_Interaction, "QSq_2Track_response", plotDir);
    //
    //Draw1DHist(rootDir_Interaction,"QSq_Error",plotDir);
    //Draw1DHist(rootDir_Interaction,"QSq_1Track_Error",plotDir);
    //Draw1DHist(rootDir_Interaction,"QSq_2Track_Error",plotDir);

    Draw1DHist(rootDir_Interaction,"QSq_Diff",plotDir);
    Draw1DHist(rootDir_Interaction,"QSq_1Track_Diff",plotDir);
    Draw1DHist(rootDir_Interaction,"QSq_2Track_Diff",plotDir);
}

void CCProtonPi0_Plotter::plotTruth_W()
{
    std::string plotDir = Folder_List::plotDir_Interaction;
    //DrawNormalizedMigrationHistogram(rootDir_Interaction, "W_response", plotDir);
    //Draw1DHist(rootDir_Interaction,"W_Error",plotDir);
    Draw1DHist(rootDir_Interaction,"W_Diff",plotDir);
}

void CCProtonPi0_Plotter::plotTruth_ShortProton()
{
    std::string plotDir = Folder_List::plotDir_Interaction;
    Draw1DHist(rootDir_Interaction,"proton_true_P_1Track",plotDir);
    Draw1DHist(rootDir_Interaction,"proton_true_KE_1Track",plotDir);
    Draw1DHist(rootDir_Interaction,"proton_true_theta_1Track",plotDir);
}

void CCProtonPi0_Plotter::plotTruth()
{
    std::cout<<"Plotting Truth Histograms"<<std::endl;
    std::string plotDir = Folder_List::plotDir_Truth;

    DrawStackedMC(rootDir_Truth, "CV_weight", plotDir);
    DrawStackedMC(rootDir_Truth, "CV_weight_2p2h", plotDir);
    DrawStackedMC(rootDir_Truth, "CV_weight_Delta", plotDir);
    DrawStackedMC(rootDir_Truth, "CV_weight_CCRES", plotDir);
    DrawStackedMC(rootDir_Truth, "CV_weight_NonRes1pi", plotDir);

    DrawStackedMC(rootDir_Truth, "h_err_2p2h", plotDir);

    DrawStackedMC(rootDir_Truth, "genie_wgt_Theta_Delta2Npi", plotDir);
    DrawStackedMC(rootDir_Truth, "updated_wgt_Theta_Delta2Npi", plotDir);

    DrawStackedMC(rootDir_Truth, "genie_wgt_MaRES", plotDir);
    DrawStackedMC(rootDir_Truth, "updated_wgt_MaRES", plotDir);

    DrawStackedMC(rootDir_Truth, "genie_wgt_MvRES", plotDir);
    DrawStackedMC(rootDir_Truth, "updated_wgt_MvRES", plotDir);

    DrawStackedMC(rootDir_Truth, "genie_wgt_Rvn1pi", plotDir);
    DrawStackedMC(rootDir_Truth, "updated_wgt_Rvn1pi", plotDir);
}

void CCProtonPi0_Plotter::plot_CV_weights()
{
    std::string plotDir = Folder_List::plotDir_Interaction;

    DrawStackedMC(rootDir_Interaction,"CV_weight",plotDir);
    DrawStackedMC(rootDir_Interaction,"CV_weight_Flux",plotDir);
    DrawStackedMC(rootDir_Interaction,"CV_weight_2p2h",plotDir);
    DrawStackedMC(rootDir_Interaction,"CV_weight_Delta",plotDir);
    DrawStackedMC(rootDir_Interaction,"CV_weight_CCRES",plotDir);
    DrawStackedMC(rootDir_Interaction,"CV_weight_NonRes1pi",plotDir);

    DrawStackedMC(rootDir_Interaction,"err_2p2h",plotDir);

    DrawStackedMC(rootDir_Interaction, "genie_wgt_VecFFCCQEshape", plotDir);
    DrawStackedMC(rootDir_Interaction, "genie_wgt_NormDISCC", plotDir);

    DrawStackedMC(rootDir_Interaction,"genie_wgt_Theta_Delta2Npi",plotDir);
    DrawStackedMC(rootDir_Interaction,"updated_wgt_Theta_Delta2Npi",plotDir);

    DrawStackedMC(rootDir_Interaction,"genie_wgt_MaRES",plotDir);
    DrawStackedMC(rootDir_Interaction,"updated_wgt_MaRES",plotDir);

    DrawStackedMC(rootDir_Interaction,"genie_wgt_MvRES",plotDir);
    DrawStackedMC(rootDir_Interaction,"updated_wgt_MvRES",plotDir);

    DrawStackedMC(rootDir_Interaction,"genie_wgt_Rvn1pi",plotDir);
    DrawStackedMC(rootDir_Interaction,"updated_wgt_Rvn1pi",plotDir);
}

void CCProtonPi0_Plotter::plotInteraction_MCOnly()
{
    std::cout<<"Plotting Interaction MC Only"<<std::endl;
    std::string plotDir = Folder_List::plotDir_Interaction;

    //plot_SignalKinematics();
    plot_CV_weights();

    Draw2DHist(rootDir_Interaction,"Enu_flux_wgt",plotDir);
    Draw2DHist(rootDir_Interaction,"Enu_cvweight",plotDir);

    //DrawStackedMC(rootDir_Interaction,"Enu",plotDir);
    //DrawStackedMC(rootDir_Interaction,"QSq",plotDir);
    //DrawStackedMC(rootDir_Interaction,"W",plotDir);
    //DrawStackedMC(rootDir_Interaction,"deltaInvMass",plotDir);

    //Draw1DHist(rootDir_Interaction,"Enu_0",plotDir);
    //plotTruth_Enu();
    //plotTruth_QSq();
    //plotTruth_W();
    //plotTruth_ShortProton();

    //DrawSignalMC(rootDir_Interaction, "W_1Track", plotDir);
    //DrawSignalMC(rootDir_Interaction, "W_2Track", plotDir);
    //DrawSignalMC(rootDir_Interaction, "deltaInvMass", plotDir);

    //Draw1DHist(rootDir_Interaction,"n_ejected_nucleons_1Track",plotDir);
    //Draw1DHist(rootDir_Interaction,"n_ejected_nucleons_2Track",plotDir);
    
    //plot_stacked_pi0_P();
    //plot_stacked_pi0_theta();

    //plot_SystematicsInfo();
    
    std::cout<<"Plotting Interaction MC Only Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plotInteraction_DataMC()
{
    std::cout<<"Plotting Interaction Data vs MC"<<std::endl;
    std::string plotDir = Folder_List::plotDir_Interaction;

    //PlotDelta();

    //DrawDataMC_Thesis(rootDir_Interaction,"Enu_1Track",plotDir);
    //DrawDataMC_Thesis(rootDir_Interaction,"Enu_2Track",plotDir);
    //DrawDataMC_Thesis(rootDir_Interaction,"Enu",plotDir);

    //DrawDataMC_Thesis(rootDir_Interaction,"QSq",plotDir);
    //DrawDataMC_Thesis(rootDir_Interaction,"WSq",plotDir);
    //DrawDataMC_Thesis(rootDir_Interaction,"W",plotDir);
    //DrawDataMC_Thesis(rootDir_Interaction,"deltaInvMass",plotDir);
    //DrawDataMC_Thesis(rootDir_Interaction,"nProtons",plotDir);

    //DrawDataMC_Thesis(rootDir_Interaction,"Enu_1Track",plotDir);
    //DrawDataMC_Thesis(rootDir_Interaction,"Enu_2Track",plotDir);
    //DrawDataMC_Thesis(rootDir_Interaction,"Enu",plotDir);
    //DrawDataMC_Thesis(rootDir_Interaction,"QSq",plotDir);
    //DrawDataMC_Thesis(rootDir_Interaction,"QSq_1Track",plotDir);
    //DrawDataMC_Thesis(rootDir_Interaction,"QSq_2Track",plotDir);
    //DrawDataMC_Thesis(rootDir_Interaction,"WSq",plotDir);
    //DrawDataMC_Thesis(rootDir_Interaction,"W",plotDir);
    //DrawDataMC_Thesis(rootDir_Interaction,"W_1Track",plotDir);
    //DrawDataMC_Thesis(rootDir_Interaction,"W_2Track",plotDir);
    //DrawDataMC_Thesis(rootDir_Interaction,"deltaInvMass",plotDir);

    //DrawDataMC_Thesis(rootDir_Interaction,"vertex_energy_1Track",plotDir);
    //DrawDataMC_Thesis(rootDir_Interaction,"vertex_evis_1Track",plotDir);
    //DrawDataMC_Thesis(rootDir_Interaction,"extra_leftover_energy_1Track",plotDir);
    //DrawDataMC_Thesis(rootDir_Interaction,"extra_muon_energy_1Track",plotDir);
    //DrawDataMC_Thesis(rootDir_Interaction,"extra_rejected_energy_1Track",plotDir);
    //DrawDataMC_Thesis(rootDir_Interaction,"extra_total_energy_1Track",plotDir);

    //DrawDataMC_Thesis(rootDir_Interaction,"vertex_energy_2Track",plotDir);
    //DrawDataMC_Thesis(rootDir_Interaction,"vertex_evis_2Track",plotDir);
    //DrawDataMC_Thesis(rootDir_Interaction,"extra_leftover_energy_2Track",plotDir);
    //DrawDataMC_Thesis(rootDir_Interaction,"extra_muon_energy_2Track",plotDir);
    //DrawDataMC_Thesis(rootDir_Interaction,"extra_rejected_energy_2Track",plotDir);
    //DrawDataMC_Thesis(rootDir_Interaction,"extra_total_energy_2Track",plotDir);

    std::cout<<"Plotting Interaction Data vs MC Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::plotMuon_MCOnly()
{
    std::cout<<"Plotting Muon MC Only"<<std::endl;
    std::string plotDir = Folder_List::plotDir_Muon;

    //Draw1DHist(rootDir_Muon,"E_error",plotDir);
    Draw1DHist(rootDir_Muon,"E_Diff",plotDir);
    //DrawNormalizedMigrationHistogram(rootDir_Muon, "reco_E_true_E", plotDir);

    //Draw1DHist(rootDir_Muon,"muon_P_shift", plotDir);
    //Draw1DHist(rootDir_Truth,"muon_P_mc_truth_all_signal", plotDir);
    //Draw1DHist(rootDir_Truth,"muon_theta_mc_truth_all_signal", plotDir);

    //DrawSignalMC(rootDir_Muon, "P", plotDir);
    //DrawStackedMC(rootDir_Muon, "P", plotDir);
    //DrawSignalMC(rootDir_Muon, "theta", plotDir);
    //DrawStackedMC(rootDir_Muon, "theta", plotDir);
    //DrawSignalMC(rootDir_Muon, "cos_theta", plotDir);
    //DrawStackedMC(rootDir_Muon, "cos_theta", plotDir);

    Draw1DHist(rootDir_Muon,"P_error",plotDir);
    Draw1DHist(rootDir_Muon,"theta_error",plotDir);
    //Draw1DHist(rootDir_Muon,"theta_Diff",plotDir);
    //Draw1DHist(rootDir_Muon,"cos_theta_error",plotDir);

    //DrawNormalizedMigrationHistogram(rootDir_Muon, "muon_P_response", plotDir);
    //DrawNormalizedMigrationHistogram(rootDir_Muon, "muon_theta_response", plotDir);

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

    //DrawStackedMC(rootDir_Proton, "P", plotDir);
    //DrawStackedMC(rootDir_Proton, "theta", plotDir);
    //DrawStackedMC(rootDir_Proton, "phi", plotDir);

    //Draw1DHist(rootDir_Proton,"E_error",plotDir);
    Draw1DHist(rootDir_Proton,"E_Diff",plotDir);
    //DrawNormalizedMigrationHistogram(rootDir_Proton, "reco_E_true_E", plotDir);

    Draw1DHist(rootDir_Proton,"P_error",plotDir);
    //Draw1DHist(rootDir_Proton,"P_Diff",plotDir);
    //DrawNormalizedMigrationHistogram(rootDir_Proton, "proton_P_response", plotDir);
    
    Draw1DHist(rootDir_Proton,"theta_error",plotDir);
    //DrawNormalizedMigrationHistogram(rootDir_Proton, "proton_theta_response", plotDir);
    
    //Draw1DHist(rootDir_Proton,"E_error",plotDir);
    //Draw1DHist(rootDir_Proton,"E_Diff",plotDir);
    
    // Energy Shift Uncertainty
    //Draw1DHist(rootDir_Proton,"energy_shift_BetheBloch",plotDir);
    //Draw1DHist(rootDir_Proton,"energy_shift_Birks",plotDir);
    //Draw1DHist(rootDir_Proton,"energy_shift_MEU",plotDir);
    //Draw1DHist(rootDir_Proton,"energy_shift_Mass",plotDir);
    //Draw1DHist(rootDir_Proton,"energy_shift_Nominal",plotDir);

    std::cout<<"Plotting Proton MC Only Finished!\n"<<std::endl;
}


void CCProtonPi0_Plotter::plotProton_DataMC()
{    
    std::cout<<"Plotting Proton Data vs MC"<<std::endl;
    std::string plotDir = Folder_List::plotDir_Proton;

    plotStandardHistograms(rootDir_Proton, plotDir);

    std::cout<<">> Plotting Unique Histograms"<<std::endl;

    //DrawDataStackedMC(rootDir_Proton,"trackLength",plotDir);
    //DrawDataStackedMC(rootDir_Proton,"trackKinked",plotDir);
    //DrawDataStackedMC(rootDir_Proton,"partScore",plotDir);

    std::cout<<"Plotting Proton Data vs MC Finished!\n"<<std::endl;
}

void CCProtonPi0_Plotter::plotPion_MCOnly()
{
    std::cout<<"Plotting Pion MC Only"<<std::endl;
    std::string plotDir = Folder_List::plotDir_Pion;

    plotTruth_Pion();

    // DrawStackedMC(rootDir_Pion,"gamma1_E",plotDir);
    // DrawStackedMC(rootDir_Pion,"gamma1_theta",plotDir);
    // DrawStackedMC(rootDir_Pion,"gamma1_ConvLength",plotDir);

    // DrawStackedMC(rootDir_Pion,"gamma2_E",plotDir);
    // DrawStackedMC(rootDir_Pion,"gamma2_theta",plotDir);
    // DrawStackedMC(rootDir_Pion,"gamma2_ConvLength",plotDir);

    //DrawStackedMC(rootDir_Pion,"photonEnergy_Asymmetry",plotDir);
    //DrawStackedMC(rootDir_Pion,"invMass",plotDir);
    //DrawStackedMC(rootDir_Pion,"cos_openingAngle",plotDir);
    //DrawStackedMC(rootDir_Pion,"P",plotDir);
    //DrawStackedMC(rootDir_Pion,"KE",plotDir);
    //DrawStackedMC(rootDir_Pion,"theta",plotDir);
    //DrawStackedMC(rootDir_Pion,"phi",plotDir);

    std::cout<<"Plotting Pion MC Only Finished!\n"<<std::endl;
}


void CCProtonPi0_Plotter::plotPion_DataMC()
{
    std::cout<<"Plotting Pion Data vs MC"<<std::endl;
    std::string plotDir = Folder_List::plotDir_Pion;

    plotStandardHistograms(rootDir_Pion, plotDir);

    std::cout<<">> Plotting Unique Histograms"<<std::endl;
    //DrawDataStackedMC(rootDir_Pion,"P_1Track",plotDir);
    //DrawDataStackedMC(rootDir_Pion,"P_2Track",plotDir);
    //DrawDataStackedMC(rootDir_Pion,"theta_1Track",plotDir);
    //DrawDataStackedMC(rootDir_Pion,"theta_2Track",plotDir);

    //DrawDataMC(rootDir_Pion,"gamma1_E",plotDir);
    //DrawDataMC(rootDir_Pion,"gamma1_theta",plotDir);
    //DrawDataMC(rootDir_Pion,"gamma1_ConvLength",plotDir);

    //DrawDataMC(rootDir_Pion,"gamma2_E",plotDir);
    //DrawDataMC(rootDir_Pion,"gamma2_theta",plotDir);
    //DrawDataMC(rootDir_Pion,"gamma2_ConvLength",plotDir);

    //DrawDataMC(rootDir_Pion,"invMass",plotDir);
    //DrawDataMC(rootDir_Pion,"photonEnergy_Asymmetry",plotDir);

    //DrawDataStackedMC(rootDir_Pion,"gamma1_E",plotDir);
    //DrawDataStackedMC(rootDir_Pion,"gamma1_theta",plotDir);
    //DrawDataStackedMC(rootDir_Pion,"gamma1_ConvLength",plotDir);

    //DrawDataStackedMC(rootDir_Pion,"gamma2_E",plotDir);
    //DrawDataStackedMC(rootDir_Pion,"gamma2_theta",plotDir);
    //DrawDataStackedMC(rootDir_Pion,"gamma2_ConvLength",plotDir);

    //DrawDataStackedMC(rootDir_Pion,"invMass",plotDir);
    //DrawDataStackedMC(rootDir_Pion,"photonEnergy_Asymmetry",plotDir);

    std::cout<<"Plotting Pion Data vs MC Finished!\n"<<std::endl;
}

void CCProtonPi0_Plotter::plotTruth_Pion()
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
        //Draw2DHist(rootDir_Pion,"bckg_signal_diff_E",plotDir);
    //
    //    Draw2DHist(rootDir_Pion,"signal_gamma1_convLength_gamma2_convLength",plotDir);
    //    Draw2DHist(rootDir_Pion,"bckg_gamma1_convLength_gamma2_convLength",plotDir);
    //    Draw2DHist(rootDir_Pion,"bckg_signal_diff_convLength",plotDir);

    //Draw1DHist(rootDir_Pion,"E_error",plotDir);
    Draw1DHist(rootDir_Pion,"E_Diff",plotDir);
    //DrawNormalizedMigrationHistogram(rootDir_Pion, "reco_E_true_E", plotDir);
    

    Draw1DHist(rootDir_Pion,"P_error",plotDir);
    //Draw1DHist(rootDir_Pion,"P_Diff",plotDir);
//    DrawNormalizedMigrationHistogram(rootDir_Pion, "pi0_P_response", plotDir);
//    
//    Draw1DHist(rootDir_Pion,"KE_error",plotDir);
//    Draw1DHist(rootDir_Pion,"KE_Diff",plotDir);
//    DrawNormalizedMigrationHistogram(rootDir_Pion, "pi0_KE_response", plotDir);
//
    Draw1DHist(rootDir_Pion,"theta_error",plotDir);
//    DrawNormalizedMigrationHistogram(rootDir_Pion, "pi0_theta_response", plotDir);

    std::cout<<"Plotting Pion True Finished!\n"<<std::endl;
}

void CCProtonPi0_Plotter::plotStandardHistograms(rootDir &dir, std::string plotDir)
{
    std::cout<<">> Plotting Standard Histograms"<<std::endl;

    DrawDataMC_Thesis(dir, "E", plotDir);
    //DrawDataMC_Thesis(dir, "P", plotDir);
    //DrawDataMC_Thesis(dir, "KE", plotDir);
    //DrawDataMC_Thesis(dir, "theta", plotDir);
    //DrawDataMC_Thesis(dir, "phi", plotDir);

    //DrawDataStackedMC(dir, "E", plotDir);
    DrawDataStackedMC(dir, "P", plotDir);
    //DrawDataStackedMC(dir, "KE", plotDir);
    //DrawDataStackedMC(dir, "theta", plotDir);
    //DrawDataStackedMC(dir, "phi", plotDir);
}

void CCProtonPi0_Plotter::plotCutHistograms_MCOnly()
{
    std::cout<<"Plotting Cut Histograms MC Only"<<std::endl;
    std::string plotDir = Folder_List::plotDir_CutHists;

    Draw2DHist(rootDir_CutHists,"bckg_signal_diff_E_cos_openingAngle",plotDir);

    // DrawStackedMC(rootDir_CutHists,"hCut_nTracks",plotDir);
    // DrawStackedMC(rootDir_CutHists,"hCut_nTracks2",plotDir);
    // DrawStackedMC(rootDir_CutHists,"hCut_nTracks_Close",plotDir);
    // DrawStackedMC(rootDir_CutHists,"hCut_nTracks_Far",plotDir);
    // DrawStackedMC(rootDir_CutHists,"hCut_nTracks_Discarded",plotDir);

    // ------------------------------------------------------------------------
    // Common
    // ------------------------------------------------------------------------
    //CutArrow Vertex_Count(3,"L"); 
    //DrawStackedMC(rootDir_CutHists,"hCut_nVertices",plotDir, 1, Vertex_Count);

    //CutArrow Proton_Count(3,"L"); 
    //DrawStackedMC(rootDir_CutHists,"hCut_nProtonCandidates",plotDir,1 , Proton_Count);

    DrawStackedMC(rootDir_CutHists,"hCut_nShowerCandidates",plotDir);
    ////DrawStackedMC(rootDir_CutHists,"hCut_1Track_nShowerCandidates",plotDir);
    ////DrawStackedMC(rootDir_CutHists,"hCut_2Track_nShowerCandidates",plotDir);

    ////CutArrow Michel(1,"L"); 
    ////DrawStackedMC(rootDir_CutHists,"hCut_Michel",plotDir,1, Michel);

    //CutArrow pi0invMass_min(60,"R"); 
    //CutArrow pi0invMass_max(200,"L"); 
    //DrawSignalMC(rootDir_CutHists,"hCut_pi0invMass",plotDir, 2, pi0invMass_min, pi0invMass_max);
    //DrawStackedMC(rootDir_CutHists,"hCut_pi0invMass",plotDir, 2, pi0invMass_min, pi0invMass_max);

    // ------------------------------------------------------------------------
    // 1 Track
    // ------------------------------------------------------------------------
    //CutArrow eVis_nuclearTarget_1Track(20,"L"); 
    //DrawStackedMC(rootDir_CutHists,"hCut_1Track_eVis_nuclearTarget",plotDir, 1, eVis_nuclearTarget_1Track);

    //CutArrow eVis_other_min_1Track(50,"R"); 
    //CutArrow eVis_other_max_1Track(2500,"L"); 
    //DrawStackedMC(rootDir_CutHists,"hCut_1Track_eVis_other",plotDir, 2, eVis_other_min_1Track, eVis_other_max_1Track);

    //CutArrow gamma1_ConvDist_1Track(14,"R"); 
    //DrawStackedMC(rootDir_CutHists,"hCut_1Track_gamma1ConvDist",plotDir, 1, gamma1_ConvDist_1Track);

    ////CutArrow gamma2_ConvDist_1Track(15,"R"); 
    ////DrawStackedMC(rootDir_CutHists,"hCut_1Track_gamma2ConvDist",plotDir, 1, gamma2_ConvDist_1Track);
    //DrawStackedMC(rootDir_CutHists,"hCut_1Track_gamma2ConvDist",plotDir);

    //DrawSignalMC(rootDir_CutHists,"hCut_1Track_pi0invMass",plotDir, 2, pi0invMass_min, pi0invMass_max);
    //DrawStackedMC(rootDir_CutHists,"hCut_1Track_pi0invMass",plotDir, 2, pi0invMass_min, pi0invMass_max);

    //CutArrow neutrinoE_1Track(20,"L"); 
    //DrawStackedMC(rootDir_CutHists,"hCut_1Track_neutrinoE",plotDir, 1, neutrinoE_1Track);

    //// ------------------------------------------------------------------------
    //// 2 Track
    //// ------------------------------------------------------------------------
    //CutArrow eVis_nuclearTarget_2Track(20,"L"); 
    //DrawStackedMC(rootDir_CutHists,"hCut_2Track_eVis_nuclearTarget",plotDir, 1, eVis_nuclearTarget_2Track);

    //CutArrow eVis_other_min_2Track(50,"R"); 
    //CutArrow eVis_other_max_2Track(2000,"L"); 
    //DrawStackedMC(rootDir_CutHists,"hCut_2Track_eVis_other",plotDir, 2, eVis_other_min_2Track, eVis_other_max_2Track);

    //CutArrow gamma1_ConvDist_2Track(14,"R"); 
    //DrawStackedMC(rootDir_CutHists,"hCut_2Track_gamma1ConvDist",plotDir, 1, gamma1_ConvDist_2Track);

    ////CutArrow gamma2_ConvDist_2Track(15,"R"); 
    ////DrawStackedMC(rootDir_CutHists,"hCut_2Track_gamma2ConvDist",plotDir, 1, gamma2_ConvDist_2Track);
    //DrawStackedMC(rootDir_CutHists,"hCut_2Track_gamma2ConvDist",plotDir);

    //DrawSignalMC(rootDir_CutHists,"hCut_2Track_pi0invMass",plotDir, 2, pi0invMass_min, pi0invMass_max);
    //DrawStackedMC(rootDir_CutHists,"hCut_2Track_pi0invMass",plotDir, 2, pi0invMass_min, pi0invMass_max);

    //CutArrow neutrinoE_2Track(20,"L"); 
    //DrawStackedMC(rootDir_CutHists,"hCut_2Track_neutrinoE",plotDir, 1, neutrinoE_2Track);

    ////CutArrow protonScore_pIDDiff(0.45,"R"); 
    ////DrawStackedMC(rootDir_CutHists,"hCut_2Track_protonScore_pIDDiff",plotDir, 1, protonScore_pIDDiff);

    //CutArrow protonScore_LLR(-10,"R"); 
    //DrawStackedMC(rootDir_CutHists,"hCut_2Track_protonScore_LLR",plotDir, 1, protonScore_LLR);

    //DrawStackedMC(rootDir_CutHists,"hCut_2Track_deltaInvMass",plotDir);

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
    DrawDataStackedMC(rootDir_CutHists,"hCut_nShowerCandidates",plotDir);
    //DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_nShowerCandidates",plotDir);
    //DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_nShowerCandidates",plotDir);

    //CutArrow Michel(1,"L"); 
    //DrawDataStackedMC(rootDir_CutHists,"hCut_Michel",plotDir,1, Michel);

    //plot_InvMass_TruthMatch_Stacked(true,true);
    //plot_InvMass_TruthMatch_Stacked(true,false);
    //plot_InvMass_TruthMatch_Stacked(false,true);
    //plot_InvMass_TruthMatch_Stacked(false,false);
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

    //DrawSignalMC(rootDir_CutHists,"hCut_1Track_pi0invMass",plotDir, 2, pi0invMass_min, pi0invMass_max);
    //DrawDataStackedMC(rootDir_CutHists,"hCut_1Track_pi0invMass",plotDir, 2, pi0invMass_min, pi0invMass_max);

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

    CutArrow gamma1_ConvDist_2Track(14,"R"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_gamma1ConvDist",plotDir, 1, gamma1_ConvDist_2Track);

    //CutArrow gamma2_ConvDist_2Track(15,"R"); 
    //DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_gamma2ConvDist",plotDir, 1, gamma2_ConvDist_2Track);
    //DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_gamma2ConvDist",plotDir);

    //DrawSignalMC(rootDir_CutHists,"hCut_2Track_pi0invMass",plotDir, 2, pi0invMass_min, pi0invMass_max);
    //DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_pi0invMass",plotDir, 2, pi0invMass_min, pi0invMass_max);

    //CutArrow neutrinoE_2Track(20,"L"); 
    //DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_neutrinoE",plotDir, 1, neutrinoE_2Track);

    CutArrow protonScore_LLR(-10,"R"); 
    DrawDataStackedMC(rootDir_CutHists,"hCut_2Track_protonScore_LLR",plotDir, 1, protonScore_LLR);

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

void CCProtonPi0_Plotter::plot_SignalKinematics()
{
    plot_SignalKinematics("truth_Enu", "all", true);
    plot_SignalKinematics("truth_QSq", "all", true);
    plot_SignalKinematics("truth_w", "all", true);
 
    plot_SignalKinematics("truth_Enu", "minos", true);
    plot_SignalKinematics("truth_QSq", "minos", true);
    plot_SignalKinematics("truth_w", "minos", true);
 
    plot_SignalKinematics("truth_Enu", "selected", true);
    plot_SignalKinematics("truth_QSq", "selected", true);
    plot_SignalKinematics("truth_w", "selected", true);

    plot_SignalKinematics("mc_incomingE", "selected", true);
    plot_SignalKinematics("mc_Q2", "selected", true);
    plot_SignalKinematics("mc_w", "selected", true);
  
    plot_SignalKinematics("reco_w", "selected", true);
 
    plot_SignalKinematics("reco_bckg_Enu", "selected", true);
    plot_SignalKinematics("reco_bckg_QSq", "selected", true);
    plot_SignalKinematics("reco_bckg_w", "selected", true);
}

void CCProtonPi0_Plotter::plot_SignalKinematics(std::string var, std::string type, bool isStacked)
{
    std::string root_dir;

    if (type.compare("all") == 0){
        std::cout<<"Plotting Signal Kinematics for All Signal Events"<<std::endl;
        root_dir = rootDir_Truth.mc;
    }else if (type.compare("minos") == 0){
        std::cout<<"Plotting Signal Kinematics for MINOS Matched Signal Events"<<std::endl;
        root_dir = rootDir_CutHists.mc;
    }else if (type.compare("selected") == 0){
        std::cout<<"Plotting Signal Kinematics for Selected Signal Events"<<std::endl;
        root_dir = rootDir_Interaction.mc;
    }else{
        std::cout<<"Wrong Signal Type Specified!"<<std::endl;
        std::cout<<"Signal Types: all, minos, selected"<<std::endl;
        return;
    }

    bool isBckg = false;
    std::size_t found = var.find("bckg");
    if (found != std::string::npos){
        isBckg = true;
    }

    std::string plotDir = Folder_List::plotDir_Interaction;
    std::cout<<"\nPlottting "<<var<<std::endl;

    TFile* f_Root = new TFile(root_dir.c_str());
    TCanvas* c1 = new TCanvas("c","c",1280,800);
    THStack *hs = new THStack("hs","Signal Events");
    TLegend *legend = new TLegend(0.75,0.5,0.95,0.9);  
    legend->SetTextFont(42);
    legend->SetTextSize(0.04);

    std::string var_name = var + "_QE";
    TH1D* h_QE = (TH1D*)f_Root->Get(var_name.c_str());
    h_QE->SetFillColor(kRed);
    h_QE->SetLineColor(kRed);
    h_QE->SetLineWidth(2);
    h_QE->SetFillStyle(3001);

    var_name = var + "_RES_1232";
    TH1D* h_RES_1232 = (TH1D*)f_Root->Get(var_name.c_str());
    h_RES_1232->SetFillColor(kGreen);
    h_RES_1232->SetLineColor(kGreen);
    h_RES_1232->SetLineWidth(2);
    h_RES_1232->SetFillStyle(3001);

    var_name = var + "_RES_1535";
    TH1D* h_RES_1535 = (TH1D*)f_Root->Get(var_name.c_str());
    h_RES_1535->SetFillColor(kGreen+1);
    h_RES_1535->SetLineColor(kGreen+1);
    h_RES_1535->SetLineWidth(2);
    h_RES_1535->SetFillStyle(3001);

    var_name = var + "_RES_1520";
    TH1D* h_RES_1520 = (TH1D*)f_Root->Get(var_name.c_str());
    h_RES_1520->SetFillColor(kGreen+2);
    h_RES_1520->SetLineColor(kGreen+2);
    h_RES_1520->SetLineWidth(2);
    h_RES_1520->SetFillStyle(3001);

    var_name = var + "_RES_Other";
    TH1D* h_RES_Other = (TH1D*)f_Root->Get(var_name.c_str());
    h_RES_Other->SetFillColor(kGreen+3);
    h_RES_Other->SetLineColor(kGreen+3);
    h_RES_Other->SetLineWidth(2);
    h_RES_Other->SetFillStyle(3001);

    var_name = var + "_Non_RES";
    TH1D* h_Non_Res = (TH1D*)f_Root->Get(var_name.c_str());
    h_Non_Res->SetFillColor(kMagenta);
    h_Non_Res->SetLineColor(kMagenta);
    h_Non_Res->SetLineWidth(2);
    h_Non_Res->SetFillStyle(3001);

    var_name = var + "_DIS";
    TH1D* h_DIS = (TH1D*)f_Root->Get(var_name.c_str());
    h_DIS->SetFillColor(kAzure);
    h_DIS->SetLineColor(kAzure);
    h_DIS->SetLineWidth(2);
    h_DIS->SetFillStyle(3001);

    var_name = var + "_2p2h";
    TH1D* h_2p2h = (TH1D*)f_Root->Get(var_name.c_str());
    h_2p2h->SetFillColor(kBlack);
    h_2p2h->SetLineColor(kBlack);
    h_2p2h->SetLineWidth(2);
    h_2p2h->SetFillStyle(3001);

    TH1D* h_Coh = NULL;
    if (isBckg){
        var_name = var + "_Coh";
        h_Coh = (TH1D*)f_Root->Get(var_name.c_str());
        h_Coh->SetFillColor(kCyan);
        h_Coh->SetLineColor(kCyan);
        h_Coh->SetLineWidth(2);
        h_Coh->SetFillStyle(3001);
    }

    legend->AddEntry(h_RES_1232, "RES: #Delta(1232)", "f");
    legend->AddEntry(h_RES_1535, "RES: N(1535)", "f");
    legend->AddEntry(h_RES_1520, "RES: N(1520)", "f");
    legend->AddEntry(h_RES_Other, "RES: Other", "f");
    legend->AddEntry(h_Non_Res, "Non Res", "f");
    legend->AddEntry(h_DIS, "DIS", "f");
    legend->AddEntry(h_2p2h, "Valencia 2p2h", "f");
    legend->AddEntry(h_QE, "QE", "f");
    if (isBckg) legend->AddEntry(h_Coh, "Coh", "f");

    if (isBckg) hs->Add(h_Coh);
    hs->Add(h_QE);
    hs->Add(h_2p2h);
    hs->Add(h_DIS);
    hs->Add(h_Non_Res);
    hs->Add(h_RES_Other);
    hs->Add(h_RES_1520);
    hs->Add(h_RES_1535);
    hs->Add(h_RES_1232);

    if (isStacked) hs->Draw();
    else hs->Draw("nostack");

    hs->GetXaxis()->SetTitle(h_RES_1232->GetXaxis()->GetTitle());
    hs->GetYaxis()->SetTitle(h_RES_1232->GetYaxis()->GetTitle());

    legend->Draw();

    std::string plot_type;
    if (isStacked) plot_type = "_" + type + "_Stacked.png";
    else plot_type = "_" + type + ".png";

    if (thesisStyle){
        //hs->GetXaxis()->SetTitle("P_{#mu} [GeV]");
        hs->GetXaxis()->SetTitleFont(62);
        hs->GetXaxis()->SetTitleSize(0.06);
        hs->GetXaxis()->CenterTitle();
        hs->GetXaxis()->SetTitleOffset(1.15);
        hs->GetXaxis()->SetLabelFont(42);
        hs->GetXaxis()->SetLabelSize(0.05);
        hs->GetXaxis()->SetNdivisions(408);

        hs->GetYaxis()->SetTitleFont(62);
        hs->GetYaxis()->SetTitleSize(0.06);
        //hs->GetYaxis()->CenterTitle();
        hs->GetYaxis()->SetTitleOffset(1.2);
        hs->GetYaxis()->SetLabelFont(42);
        hs->GetYaxis()->SetLabelSize(0.05);
        TGaxis::SetMaxDigits(3);
    }

    std::string out_name = plotDir + var + plot_type; 
    c1->Print(out_name.c_str(),"png");

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
    rootDir_PC.mc = Folder_List::rootDir_PC_mc;
    rootDir_PC.data = Folder_List::rootDir_PC_data;

    rootDir_Truth.mc = Folder_List::rootDir_Truth_mc;
    rootDir_Truth.data = Folder_List::rootDir_Truth_data;

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
    //DrawMnvH1D(rootDir_GENIEXSec,"muon_P_xsec",plotDir);   
    //DrawMnvH1D(rootDir_GENIEXSec,"muon_theta_xsec",plotDir);   
    //DrawMnvH1D(rootDir_GENIEXSec,"pi0_P_xsec",plotDir);   
    //DrawMnvH1D(rootDir_GENIEXSec,"pi0_KE_xsec",plotDir);   
    //DrawMnvH1D(rootDir_GENIEXSec,"pi0_theta_xsec",plotDir);   
    //DrawMnvH1D(rootDir_GENIEXSec,"QSq_xsec",plotDir);   

    // Plot Comparison
    TFile* f_xsec_mc = new TFile(rootDir_GENIEXSec.mc.c_str());
    TFile* f_xsec_data = new TFile(rootDir_CrossSection.mc.c_str());
    MnvH1D* mc;
    MnvH1D* data;

    mc = (MnvH1D*)f_xsec_mc->Get("muon_P_xsec");
    data = (MnvH1D*)f_xsec_data->Get("muon_P_xsec"); 
    DrawDataMC(data,mc,"muon_P_MC_GENIE",plotDir,true);

    mc = (MnvH1D*)f_xsec_mc->Get("muon_theta_xsec");
    data = (MnvH1D*)f_xsec_data->Get("muon_theta_xsec"); 
    DrawDataMC(data,mc,"muon_theta_MC_GENIE",plotDir,true);

    mc = (MnvH1D*)f_xsec_mc->Get("pi0_P_xsec");
    data = (MnvH1D*)f_xsec_data->Get("pi0_P_xsec"); 
    DrawDataMC(data,mc,"pi0_P_MC_GENIE",plotDir,true);

    mc = (MnvH1D*)f_xsec_mc->Get("pi0_KE_xsec");
    data = (MnvH1D*)f_xsec_data->Get("pi0_KE_xsec"); 
    DrawDataMC(data,mc,"pi0_KE_MC_GENIE",plotDir,true);

    mc = (MnvH1D*)f_xsec_mc->Get("pi0_theta_xsec");
    data = (MnvH1D*)f_xsec_data->Get("pi0_theta_xsec"); 
    DrawDataMC(data,mc,"pi0_theta_MC_GENIE",plotDir,true);

    mc = (MnvH1D*)f_xsec_mc->Get("QSq_xsec");
    data = (MnvH1D*)f_xsec_data->Get("QSq_xsec"); 
    DrawDataMC(data,mc,"QSq_MC_GENIE",plotDir,true);

    std::cout<<"Done!"<<std::endl;
}

void CCProtonPi0_Plotter::PlotFluxHistograms()
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;
   
    //PlotFluxComparison(plotDir);
    //PlotFluxRatio(plotDir);
    PlotFluxRebinned(plotDir);
}

void CCProtonPi0_Plotter::GetFlux()
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    std::string rootDir_mc = rootDir_Muon.mc;

    TFile* f_mc = new TFile(rootDir_mc.c_str());
   
    MnvH1D* mc = (MnvH1D*)f_mc->Get("muon_theta_mc_reco_all");
    DrawErrorSummary(mc,"muon_theta",plotDir);

    std::string rootDir_mc2 = rootDir_CutHists.mc;

    TFile* f_mc2 = new TFile(rootDir_mc2.c_str());
   
    MnvH1D* mc2 = (MnvH1D*)f_mc2->Get("invMass_mc_reco_all");
    DrawErrorSummary(mc2,"invMass_mc_reco_all",plotDir);

}

void CCProtonPi0_Plotter::PlotFluxComparison(std::string plotDir)
{
    // Get Histograms
    delete frw;
    frw = new FluxReweighter(14, applyNuEConstraint, FluxReweighter::minerva1, new_flux, old_flux);
    MnvH1D* h_flux_minerva1 = new MnvH1D (*(frw->GetFluxReweighted(14)));

    delete frw;
    frw = new FluxReweighter(14, applyNuEConstraint, FluxReweighter::minervaLE_FHC, new_flux, old_flux);
    MnvH1D* h_flux_minervaLE_FHC = new MnvH1D (*(frw->GetFluxReweighted(14)));
    MnvH1D* h_flux_generated = new MnvH1D (*(frw->GetFluxGenerated(14)));

    delete frw;
    frw = new FluxReweighter(14, applyNuEConstraint, FluxReweighter::minerva13, new_flux, old_flux);
    MnvH1D* h_flux_minerva13 = new MnvH1D (*(frw->GetFluxReweighted(14)));


    // Get Integrals in Signal Region
    double area_generated = h_flux_generated->Integral(4,30,"width") * pow(10,6);
    double area_minerva1 = h_flux_minerva1->Integral(4,30,"width") * pow(10,6);
    double area_minervaLE_FHC = h_flux_minervaLE_FHC->Integral(4,30,"width") * pow(10,6);
    double area_minerva13 = h_flux_minerva13->Integral(4,30,"width") * pow(10,6);
   
    // Scale Histograms for Plotting
    double norm_bin_width = h_flux_minerva1->GetNormBinWidth();
    h_flux_minerva1->Scale(norm_bin_width,"width");
    norm_bin_width = h_flux_minervaLE_FHC->GetNormBinWidth();
    h_flux_minervaLE_FHC->Scale(norm_bin_width,"width");
    norm_bin_width = h_flux_generated->GetNormBinWidth();
    h_flux_generated->Scale(norm_bin_width,"width");
    norm_bin_width = h_flux_minerva13->GetNormBinWidth();
    h_flux_minerva13->Scale(norm_bin_width,"width");

    TCanvas* c = new TCanvas("c","c",1280,800);

    // Set Histogram Ranges to Signal Definition
    h_flux_generated->GetXaxis()->SetRangeUser(0,20);
    h_flux_minerva1->GetXaxis()->SetRangeUser(0,20);
    h_flux_minervaLE_FHC->GetXaxis()->SetRangeUser(0,20);
    h_flux_minerva13->GetXaxis()->SetRangeUser(0,20);

    // Style Histograms
    h_flux_generated->SetLineColor(kBlack);
    h_flux_minerva1->SetLineColor(kMagenta);
    h_flux_minervaLE_FHC->SetLineColor(kBlue);
    h_flux_minerva13->SetLineColor(kRed);

    // Draw Histograms
    h_flux_generated->Draw("HIST");
    h_flux_minerva1->Draw("HIST SAME");
    h_flux_minervaLE_FHC->Draw("HIST SAME");
    h_flux_minerva13->Draw("HIST SAME");
   
    // Add Legend
    TLegend *legend = new TLegend(0.7,0.8,0.9,0.9);  
    legend->AddEntry(h_flux_generated, "Old Flux");
    legend->AddEntry(h_flux_minerva1, "minerva1");
    legend->AddEntry(h_flux_minervaLE_FHC, "minervaLE_FHC");
    legend->AddEntry(h_flux_minerva13, "minerva13");
    legend->Draw();

    // Add Line for signal region
    double hist_max = h_flux_generated->GetMaximum();
    TLine line;
    line.SetLineStyle(2);
    line.SetLineWidth(2);
    line.SetLineColor(kBlue);
    line.DrawLine(1.5,0,1.5,hist_max);


    // Add Text
    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.03);
    text.DrawLatex(0.65,0.74,"Signal Region Integrals");
    text.DrawLatex(0.65,0.70,Form("%s%3.2f%s", "Old Flux = ", area_generated," 10^{-6}"));
    text.DrawLatex(0.65,0.66,Form("%s%3.2f%s", "minerva1 = ", area_minerva1," 10^{-6}"));
    text.DrawLatex(0.65,0.62,Form("%s%3.2f%s", "minervaLE_FHC = ", area_minervaLE_FHC," 10^{-6}"));
    text.DrawLatex(0.65,0.58,Form("%s%3.2f%s", "minerva13 = ", area_minerva13," 10^{-6}"));

    std::string out_name = plotDir + "ReweightedFlux.png";
    c->Print(out_name.c_str(), "png");
    
    delete legend;
    delete c;

}


void CCProtonPi0_Plotter::PlotFluxRatio(std::string plotDir)
{
    delete frw;
    frw = new FluxReweighter(14, applyNuEConstraint, FluxReweighter::minervaLE_FHC, new_flux, old_flux);
    MnvH1D* h_flux_minervaLE_FHC = new MnvH1D (*(frw->GetFluxReweighted(14)));
    MnvH1D* h_flux_generated = new MnvH1D (*(frw->GetFluxGenerated(14)));
  
    h_flux_minervaLE_FHC->GetXaxis()->SetRangeUser(0,20);
    h_flux_generated->GetXaxis()->SetRangeUser(0,20);

    // Plot New/Old Flux Ratio
    TCanvas* c = new TCanvas("c","c",1280,800);
    MnvH1D* h_flux_ratio = new MnvH1D (*h_flux_minervaLE_FHC);
    h_flux_ratio->Divide(h_flux_minervaLE_FHC, h_flux_generated);

    h_flux_ratio->SetTitle("New/Old Flux Ratio");
    h_flux_ratio->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
    h_flux_ratio->GetYaxis()->SetTitle("New/Old");
    h_flux_ratio->SetLineColor(kBlack);
    h_flux_ratio->SetMinimum(0.6);
    h_flux_ratio->SetMaximum(1.6);
    h_flux_ratio->Draw("HIST");

    TLine line;
    line.SetLineStyle(2);
    line.SetLineWidth(2);
    line.SetLineColor(kRed);
    line.DrawLine(0,1.0,20,1.0);

    // Add Line for signal region
    line.SetLineColor(kBlue);
    line.DrawLine(1.5,0.6,1.5,1.6);



    std::string out_name = plotDir + "Flux_Ratio.png";
    c->Print(out_name.c_str(), "png");

    delete c;

}

void CCProtonPi0_Plotter::PlotFluxRebinned(std::string plotDir)
{
    delete frw;
    frw = new FluxReweighter(14, applyNuEConstraint, FluxReweighter::minervaLE_FHC, new_flux, old_flux);
    MnvH1D* h_flux_minervaLE_FHC = new MnvH1D (*(frw->GetFluxReweighted(14)));
    
    std::string root_dir = rootDir_Interaction.data;
    TFile* f_Root = new TFile(root_dir.c_str());
    MnvH1D* h_flux_rebinned = GetMnvH1D(f_Root, "Enu_all");
    h_flux_rebinned->Reset(); 
  
    int nBins = h_flux_rebinned->GetNbinsX();
    for (int i = 1; i <= nBins; i++){
        double low = h_flux_rebinned->GetBinLowEdge(i);
        double up = h_flux_rebinned->GetBinLowEdge(i+1);
        double content = GetFluxHistContent(h_flux_minervaLE_FHC,low,up);

        h_flux_rebinned->SetBinContent(i,content);
    }

    std::cout<<"Integrate Original  = "<<Integrate(h_flux_minervaLE_FHC, 1,30)*pow(10,6)<<std::endl;
    std::cout<<"Integrate Rebinned = "<<Integrate(h_flux_rebinned, 1,12)*pow(10,6)<<std::endl;
    //printBins(h_flux_rebinned, "Rebinned", true); 
    //printBins(h_flux_minervaLE_FHC, "Original", true); 
    
    // Scale Histograms for plotting
    double norm_bin_width = h_flux_minervaLE_FHC->GetNormBinWidth();
    h_flux_minervaLE_FHC->Scale(norm_bin_width,"width");
    for (int i = 1; i <= nBins; i++){
        double content = h_flux_rebinned->GetBinContent(i);
        double bin_width = h_flux_rebinned->GetBinWidth(i);
        h_flux_rebinned->SetBinContent(i,content*(0.5/bin_width));
    }

    PlotFluxRebinned(h_flux_minervaLE_FHC, h_flux_rebinned, plotDir);
}

void CCProtonPi0_Plotter::PlotFluxRebinned(MnvH1D* original, MnvH1D* rebinned, std::string plotDir)
{
    TCanvas* c = new TCanvas("c","c",1280,800);
    
    if (thesisStyle){
        rebinned->GetXaxis()->SetTitleFont(62);
        rebinned->GetXaxis()->SetTitleSize(0.06);
        rebinned->GetXaxis()->CenterTitle();
        rebinned->GetXaxis()->SetTitleOffset(1.15);
        rebinned->GetXaxis()->SetLabelFont(42);
        rebinned->GetXaxis()->SetLabelSize(0.05);
        rebinned->GetXaxis()->SetNdivisions(408);

        rebinned->GetYaxis()->SetTitleFont(62);
        rebinned->GetYaxis()->SetTitleSize(0.06);
        //rebinned->GetYaxis()->CenterTitle();
        rebinned->GetYaxis()->SetTitleOffset(1.2);
        rebinned->GetYaxis()->SetLabelFont(42);
        rebinned->GetYaxis()->SetLabelSize(0.05);
        TGaxis::SetMaxDigits(3);
    }
 
    rebinned->GetYaxis()->SetTitle("#nu_{#mu}s/cm^{2}/POT/GeV");
    rebinned->SetLineWidth(2);
    rebinned->SetLineColor(kBlack);
    rebinned->Draw("HIST");

    original->GetXaxis()->SetRangeUser(0,20);
    original->SetLineColor(kRed);
    original->SetLineStyle(2);
    original->SetLineWidth(2);
    original->Draw("HIST SAME");

    //double area_original = original->Integral(1,30,"width") * pow(10,8);
    //double area_rebinned = rebinned->Integral("width") * pow(10,8);

    // Add Text
    //TLatex text;
    //text.SetNDC();
    //text.SetTextSize(0.03);
    //text.DrawLatex(0.65,0.74,"Integrals");
    //text.DrawLatex(0.65,0.70,Form("%s%3.2f%s", "Original Flux = ", area_original," 10^{-8}"));
    //text.DrawLatex(0.65,0.66,Form("%s%3.2f%s", "Rebinned Flux = ", area_rebinned," 10^{-8}"));
   
    // Add Legend
    TLegend *legend = new TLegend(0.6,0.7,0.9,0.9);  
    legend->SetTextFont(42);
    legend->SetTextSize(0.04);
    legend->AddEntry(original, "Original");
    legend->AddEntry(rebinned, "Rebinned");
    legend->Draw();

    gStyle->SetOptStat(0); 
    
    std::string out_name = plotDir + "Flux_Rebinned.png";
    c->Print(out_name.c_str(), "png");

    delete legend;
    delete c;
}
double CCProtonPi0_Plotter::GetFluxHistContent(MnvH1D* hist, double low1, double low2)
{
    std::cout<<"low = "<<low1<<" up = "<<low2<<std::endl;
    double total = 0.0;
    int nBins = hist->GetNbinsX();
    for (int i = 1; i <= nBins; ++i){
        double current_low = hist->GetBinLowEdge(i); 
        if (current_low < low1) continue;
        if (current_low == low2) break;
        total += hist->GetBinContent(i);
        std::cout<<i<<std::endl;
    }

    std::cout<<"Total = "<<total<<std::endl;
    return total;
}

double CCProtonPi0_Plotter::Integrate(MnvH1D* hist, int start, int end)
{
    std::cout<<"Start Bin Low Edge = "<<hist->GetBinLowEdge(start)<<std::endl;
    std::cout<<"End Bin High Edge = "<<hist->GetBinLowEdge(end)+hist->GetBinWidth(end)<<std::endl;

    double total = 0.0;
    int nBins = hist->GetNbinsX();
    for (int i = start; i <= end && i <= nBins; ++i){
        total += hist->GetBinContent(i);
    }

    return total;
}

void CCProtonPi0_Plotter::PlotXSecVar(std::string var_name, std::string data_var, std::string mc_var, std::string plotDir, std::string plotName)
{
    TFile* f_xsec_mc = new TFile(rootDir_CrossSection.mc.c_str());
    TFile* f_xsec_data = new TFile(rootDir_CrossSection.data.c_str());

    bool isXSec;
    if (data_var.compare("xsec") == 0 ) isXSec = true;
    else isXSec = false;

    data_var = var_name + "_" + data_var;
    mc_var = var_name + "_" + mc_var;
    plotName = var_name + "_" + plotName;

    MnvH1D* data = GetMnvH1D(f_xsec_data, data_var);
    MnvH1D* mc = GetMnvH1D(f_xsec_mc, mc_var);

    // Remove Error Bands from MC
    mc->ClearAllErrorBands();
   
    if (thesisStyle){
        DrawDataMC_Thesis(data, mc, plotName, plotDir, isXSec);
    }else{
        DrawDataMC(data, mc, plotName, plotDir, isXSec);
    }

    // Rename for Error Summary
    data_var = "data_" + data_var;
    mc_var = "mc_" + mc_var;
    DrawErrorSummary_PaperStyle(data, data_var, plotDir);
    DrawErrorSummary_PaperStyle(mc, mc_var, plotDir);

    delete data;
    delete mc;
    delete f_xsec_mc;
    delete f_xsec_data;
}

void CCProtonPi0_Plotter::PlotXSecVar_BeforeFSI(std::string var_name, std::string plotDir)
{
    TFile* f_xsec_mc = new TFile(rootDir_CrossSection.mc.c_str());
    TFile* f_xsec_data = new TFile(rootDir_CrossSection.data.c_str());

    std::string data_var = var_name + "_xsec";
    std::string mc_var = var_name + "_xsec_AfterFSI";
    std::string mc_var_BeforeFSI = var_name + "_xsec_BeforeFSI";

    MnvH1D* data = GetMnvH1D(f_xsec_data, data_var);
    MnvH1D* mc = GetMnvH1D(f_xsec_mc, mc_var);
    MnvH1D* mc_BeforeFSI = GetMnvH1D(f_xsec_mc, mc_var_BeforeFSI);

    // Remove Error Bands from MC
    mc->ClearAllErrorBands();
    mc_BeforeFSI->ClearAllErrorBands();

    //DrawDataMC_BeforeFSI(data, mc, mc_BeforeFSI, var_name, plotDir);
    DrawDataMC_PaperStyle(data, mc, mc_BeforeFSI, var_name, plotDir);
  
    delete data;
    delete mc;
    delete mc_BeforeFSI;
    delete f_xsec_mc;
    delete f_xsec_data;
}

void CCProtonPi0_Plotter::PlotXSecVar_IntType(std::string var_name, std::string plotDir)
{
    TFile* f_xsec_mc = new TFile(rootDir_CrossSection.mc.c_str());
    TFile* f_xsec_data = new TFile(rootDir_CrossSection.data.c_str());

    std::string data_var = var_name + "_xsec";
    std::string mc_var = var_name + "_xsec";
    std::string hist_name;

    MnvH1D* data = GetMnvH1D(f_xsec_data, data_var);
    MnvH1D* mc = GetMnvH1D(f_xsec_mc, mc_var);
    mc->ClearAllErrorBands();
   
    MnvH1D* temp = NULL;
    std::vector<MnvH1D*> mc_IntType;
    for (int i = 0; i < nIntType; ++i){
        hist_name = var_name + "_xsec_IntType_" + std::to_string((long long int)i);
        temp = GetMnvH1D(f_xsec_mc, hist_name);
        temp->ClearAllErrorBands();
        mc_IntType.push_back(temp);
    }

    DrawDataMC_IntType(data, mc, mc_IntType, var_name, plotDir);
  
    mc_IntType.clear();

    delete data;
    delete mc;
    delete f_xsec_mc;
    delete f_xsec_data;
}

void CCProtonPi0_Plotter::PlotXSecVar_FSIType(std::string var_name, std::string plotDir)
{
    TFile* f_xsec_mc = new TFile(rootDir_CrossSection.mc.c_str());
    TFile* f_xsec_data = new TFile(rootDir_CrossSection.data.c_str());

    std::string data_var = var_name + "_xsec";
    std::string mc_var = var_name + "_xsec";
    std::string hist_name;

    MnvH1D* data = GetMnvH1D(f_xsec_data, data_var);
    MnvH1D* mc = GetMnvH1D(f_xsec_mc, mc_var);
    mc->ClearAllErrorBands();
   
    MnvH1D* temp = NULL;
    std::vector<MnvH1D*> mc_FSIType;
    for (int i = 0; i < nFSIType; ++i){
        hist_name = var_name + "_xsec_FSIType_" + std::to_string((long long int)i);
        temp = GetMnvH1D(f_xsec_mc, hist_name);
        temp->ClearAllErrorBands();
        mc_FSIType.push_back(temp);
    }

    DrawDataMC_FSIType(data, mc, mc_FSIType, var_name, plotDir);
  
    mc_FSIType.clear();

    delete data;
    delete mc;
    delete f_xsec_mc;
    delete f_xsec_data;
}

void CCProtonPi0_Plotter::PlotXSecVar_WithMiniBoone(std::string var_name, std::string plotDir)
{
    TFile* f_xsec_mc = new TFile(rootDir_CrossSection.mc.c_str());
    TFile* f_xsec_data = new TFile(rootDir_CrossSection.data.c_str());

    std::string data_var = var_name + "_xsec";
    std::string mc_var = var_name + "_xsec";

    MnvH1D* data = GetMnvH1D(f_xsec_data, data_var);
    MnvH1D* mc = GetMnvH1D(f_xsec_mc, mc_var);

    // Remove Error Bands from MC
    mc->ClearAllErrorBands();
   
    // MiniBoone Data
    double x[11] = {0.05, 0.125, 0.175, 0.225, 0.275, 0.35, 0.45, 0.55, 0.70, 0.90, 1.20};
    double y[11] = {3.51, 19.04, 23.50, 20.71, 13.59, 9.75, 5.29, 3.05, 1.36, 0.62, 0.14};
    TGraph* otherData = new TGraph(11,x,y);

    DrawDataMC_WithOtherData(data, mc, otherData, var_name, "MiniBooNE", plotDir);
  
    delete data;
    delete mc;
    delete f_xsec_mc;
    delete f_xsec_data;
}

void CCProtonPi0_Plotter::plotCrossSection()
{
    plot_muon_P = true;
    plot_muon_theta = false;
    plot_pi0_P = false;
    plot_pi0_KE = false;
    plot_pi0_theta = false;
    plot_QSq = false;
    plot_Enu = false;
    plot_W = false;

    //plotOriginalData();
    //plotBackgroundEstimated();
    //plotBackgroundSubtracted();
    //plotUnfolded();
    //plotEfficiencyCorrected();
    //plotFluxIntegrated();
    plotXSec();
}


void CCProtonPi0_Plotter::UnfoldingStudy()
{
    UnfoldingStudy_Iterations("muon_P");
    UnfoldingStudy_Iterations("muon_theta");
    UnfoldingStudy_Iterations("pi0_P");
    UnfoldingStudy_Iterations("pi0_KE");
    UnfoldingStudy_Iterations("pi0_theta");
    UnfoldingStudy_Iterations("QSq");
    UnfoldingStudy_Iterations("Enu");
    //UnfoldingStudy_muon_P();
    //UnfoldingStudy_muon_theta();
    //UnfoldingStudy_muon_cos_theta();
    //UnfoldingStudy_pi0_P();
    //UnfoldingStudy_pi0_KE();
    //UnfoldingStudy_pi0_theta();

    //PlotUnfolding_TruthComparison();
    //PlotUnfolding_Migration();
}

void CCProtonPi0_Plotter::Systematics()
{
    //Systematics_CheckErrorSummary(rootDir_CrossSection.mc, "h_flux_rebinned");
    //Systematics_CheckErrorSummary(rootDir_CutHists.mc, "invMass_mc_reco_all");
    //Systematics_CheckErrorSummary(rootDir_CrossSection.mc, "invMass_mc_reco_all");
    
    //Systematics_XSec();
    //Systematics_invMass();
    
    Systematics_WriteTables("muon_P");
//    Systematics_WriteTables("muon_theta");
//    Systematics_WriteTables("pi0_P");
//    Systematics_WriteTables("pi0_KE");
//    Systematics_WriteTables("pi0_theta");
//    Systematics_WriteTables("QSq");
//    Systematics_WriteTables("Enu");
   
}

void CCProtonPi0_Plotter::PlotDelta()
{
    std::string plotDir = Folder_List::plotDir_Interaction;
    
    TFile* f_data = new TFile(rootDir_Interaction.data.c_str());
    TFile* f_mc = new TFile(rootDir_Interaction.mc.c_str());
 
    MnvH1D* data = GetMnvH1D(f_data, "DeltaTransverse_data");
    MnvH1D* mc = GetMnvH1D(f_mc, "DeltaTransverse_mc");
    DrawDataMC(data, mc, "DeltaTransverseMomentum", plotDir, false);

    Draw2DHist(rootDir_Interaction, "DeltaTransverse_mc_res", plotDir);

}

void CCProtonPi0_Plotter::plot_SystematicsInfo()
{
    std::string plotDir = Folder_List::plotDir_Interaction;
    
    Draw1DHist(rootDir_Interaction,"normal_rand_numbers",plotDir);
    Draw1DHist(rootDir_Interaction,"em_shift_rand_numbers",plotDir);
    Draw1DHist(rootDir_Interaction,"muonP_shift_rand_numbers",plotDir);
    Draw1DHist(rootDir_Interaction,"muon_theta_shift_rand_numbers",plotDir);
    Draw1DHist(rootDir_Interaction,"Birks_shift_rand_numbers",plotDir);

    Draw1DHist(rootDir_Interaction,"Err_NeutronResponse",plotDir);
    Draw1DHist(rootDir_Interaction,"Err_PionResponse",plotDir);
    Draw1DHist(rootDir_Interaction,"Err_MuonTracking",plotDir);
}


void CCProtonPi0_Plotter::W_Studies()
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    // W Distributions
    //DrawDataStackedMC(rootDir_Interaction,"W",plotDir);
    //DrawDataStackedMC(rootDir_Interaction,"W_1Track",plotDir);
    //DrawDataStackedMC(rootDir_Interaction,"W_2Track",plotDir);
    //DrawDataStackedMC(rootDir_Interaction,"deltaInvMass",plotDir);

    // W Quality
    //Draw1DHist(rootDir_Interaction,"W_Error",plotDir);
    //Draw1DHist(rootDir_Interaction,"W_Diff",plotDir);
    //Draw2DHist(rootDir_Interaction, "W_response", plotDir);

    //DrawDataStackedMC_WithSignalTypes(rootDir_Interaction, "W_p_pi0", plotDir);
    //DrawDataStackedMC_WithSignalTypes(rootDir_Interaction, "W_All", plotDir);
    //DrawDataStackedMC_WithSignalTypes(rootDir_Interaction, "W_1", plotDir);
    //DrawDataStackedMC_WithSignalTypes(rootDir_Interaction, "W_2", plotDir);

    //printBins_W();
    //init_W_FitResults();
    //plot_W_FitResults();
    plot_W_FitMinuit(1.0,1.0,1.0);
    plot_W_FitMinuit(0.931897,0.01,1.78231);
}

void CCProtonPi0_Plotter::BckgSubtraction_Studies()
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;
    //DrawDataStackedMC_WithSignalTypes(rootDir_Interaction, "pi0_invMass_All", plotDir);
    //DrawDataStackedMC_WithSignalTypes(rootDir_Interaction, "pi0_invMass_1Track", plotDir);
    //DrawDataStackedMC_WithSignalTypes(rootDir_Interaction, "pi0_invMass_2Track", plotDir);
    //DrawDataStackedMC_WithSignalTypes(rootDir_Interaction, "pi0_invMass_DeltaRES", plotDir);
 
    //DrawDataMC(rootDir_Interaction, "pi0_invMass_All_0", plotDir);
    //DrawDataMC(rootDir_Interaction, "pi0_invMass_1Track_0", plotDir);
    //DrawDataMC(rootDir_Interaction, "pi0_invMass_2Track_0", plotDir);
    DrawDataMC(rootDir_Interaction, "pi0_invMass_DeltaRES_0", plotDir);
   
    //Calc_Normalized_NBackground("pi0_invMass_All");
    //Calc_Normalized_NBackground("pi0_invMass_1Track");
    //Calc_Normalized_NBackground("pi0_invMass_2Track");
    Calc_Normalized_NBackground("pi0_invMass_DeltaRES");
}

double CCProtonPi0_Plotter::Calc_Normalized_NBackground(std::string var_name)
{
    std::string var;

    TFile* f_data = new TFile(Folder_List::rootDir_Interaction_data.c_str());
    TFile* f_mc = new TFile(Folder_List::rootDir_Interaction_mc.c_str());

    var = var_name + "_0"; 
    MnvH1D* h_data = GetMnvH1D(f_data, var);
    double nAll = h_data->Integral();

    var = var_name + "_2"; 
    MnvH1D* h_mc_bckg = GetMnvH1D(f_mc,var);
    h_mc_bckg->Scale(POT_ratio);
    double nBckg = h_mc_bckg->Integral();

    h_data->Add(h_mc_bckg, -1);
    double nSignal = h_data->Integral();

    delete h_data;
    delete h_mc_bckg;
    delete f_data;
    delete f_mc;

    std::cout<<"Background Calculated for "<<var_name<<std::endl;
    std::cout<<"\tAll = "<<nAll<<std::endl;
    std::cout<<"\tSignal = "<<nSignal<<std::endl;
    std::cout<<"\tBckg = "<<nBckg<<std::endl;

    return nBckg;
}

void CCProtonPi0_Plotter::QSq_Studies()
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    // QSq Distributions
    //DrawDataStackedMC(rootDir_Interaction,"QSq",plotDir);
    //DrawDataStackedMC(rootDir_Interaction,"QSq_1Track",plotDir);
    //DrawDataStackedMC(rootDir_Interaction,"QSq_2Track",plotDir);

    // QSq Quality
    //Draw1DHist(rootDir_Interaction,"QSq_Error",plotDir);
    //Draw1DHist(rootDir_Interaction,"QSq_Diff",plotDir);
    //Draw2DHist(rootDir_Interaction, "QSq_response", plotDir);
    //Save2DHistPoints(rootDir_Interaction, "QSq_response", plotDir);

    // Neutrino Energy Limit
    //Draw_QSq_EnuLimit();
    
    // MaRES Fit
    //Draw_QSq_MaRES_Plots();
    Draw_QSq_MaRES_Fit(false); // POT
    //Draw_QSq_MaRES_Fit(true); // Area
    //Draw_QSq_MaRES_AreaNorm();
    
    // Delta Suppression 
    //Draw_QSq_DeltaSuppression();
    //Draw_QSq_DeltaSuppression_AllPlots();
    //Draw_QSq_DeltaSuppression_v2("muon_P_xsec");
    //Draw_QSq_DeltaSuppression_v2("muon_theta_xsec");
    //Draw_QSq_DeltaSuppression_v2("pi0_P_xsec");
    //Draw_QSq_DeltaSuppression_v2("pi0_KE_xsec");
    //Draw_QSq_DeltaSuppression_v2("pi0_theta_xsec");
    //Draw_QSq_DeltaSuppression_v2("QSq_xsec");
    //Draw_QSq_DeltaSuppression_v2("Enu_xsec");
    //Draw_QSq_DeltaSuppression_v2("W_xsec");
}

double CCProtonPi0_Plotter::Calc_ChiSq_dof(double* data, double* expected, int nPoints, int nPars)
{
    double ChiSq = 0.0;
    
    for (int i = 0; i < nPoints; ++i){
        ChiSq += std::pow((data[i]-expected[i]),2) / expected[i];
    }

    double dof = nPoints - nPars;
    double ChiSq_dof = ChiSq / dof;

    return ChiSq_dof;
}

double CCProtonPi0_Plotter::Calc_ChiSq(TH1* data, TH1* MC, int min_bin, int max_bin)
{
    double ChiSq = 0.0;
    double nPoints = 0.0;
  
    for (int i = min_bin; i <= max_bin; ++i){
        double nData = data->GetBinContent(i);
        double nMC = MC->GetBinContent(i);
       
        ChiSq += std::pow((nData-nMC),2) / nMC;
        nPoints++;
    }

    std::cout<<"nPoints Used in ChiSq = "<<nPoints<<std::endl;

    return ChiSq;
}

double CCProtonPi0_Plotter::Calc_ChiSq(TH1* data, TH1* MC)
{
    int nBins = data->GetNbinsX();
    double ChiSq = 0.0;
    double nPoints = 0.0;
   
    for (int i = 1; i <= nBins; ++i){
        double nData = data->GetBinContent(i);
        double nMC = MC->GetBinContent(i);
       
        ChiSq += std::pow((nData-nMC),2) / nMC;
        nPoints++;
    }

    std::cout<<"nPoints Used in ChiSq = "<<nPoints<<std::endl;

    return ChiSq;
}

void CCProtonPi0_Plotter::DeltaRes_Studies()
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;
    
    Draw1DHist(rootDir_Interaction, "resID", plotDir);
    Draw1DHist(rootDir_Interaction, "resID_theta", plotDir);
//
//    DrawDataMC_Signal(rootDir_Interaction, "deltaInvMass", plotDir, 883.973);
//    DrawDataMC_Signal(rootDir_Interaction, "Delta_pi_P", plotDir, 883.973);
//    DrawDataMC_Signal(rootDir_Interaction, "Delta_pi_theta", plotDir, 883.973);
//    DrawDataMC_Signal(rootDir_Interaction, "Delta_pi_phi", plotDir, 883.973);
//
//    DrawMnvH2D_Signal(rootDir_Interaction, "Delta_pi_phi_theta", plotDir, 883.973, false);
//    DrawMnvH2D_Signal(rootDir_Interaction, "Delta_pi_phi_theta", plotDir, 883.973, true);
//
//    DrawMnvH2D_Signal(rootDir_Interaction, "Delta_pi_P_theta", plotDir, 883.973, false);
//    DrawMnvH2D_Signal(rootDir_Interaction, "Delta_pi_P_theta", plotDir, 883.973, true);
}

void CCProtonPi0_Plotter::Studies_2p2h()
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;
    
    DrawDataStackedMC_Signal(rootDir_Interaction, "vertex_energy_1Track", plotDir, 1288.07);
    DrawDataStackedMC_Signal(rootDir_Interaction, "vertex_energy_2Track", plotDir, 1709.03);
    DrawDataStackedMC_Signal(rootDir_Interaction, "vertex_energy_All", plotDir, 2997.1);

    DrawDataStackedMC_Signal(rootDir_Interaction, "vertex_evis_1Track", plotDir, 1288.07);
    DrawDataStackedMC_Signal(rootDir_Interaction, "vertex_evis_2Track", plotDir, 1709.03);
    DrawDataStackedMC_Signal(rootDir_Interaction, "vertex_evis_All", plotDir, 2997.1);

    DrawDataMC_Signal(rootDir_Interaction, "vertex_energy_1Track", plotDir, 1288.07);
    DrawDataMC_Signal(rootDir_Interaction, "vertex_energy_2Track", plotDir, 1709.03);
    DrawDataMC_Signal(rootDir_Interaction, "vertex_energy_All", plotDir, 2997.1);

    DrawDataMC_Signal(rootDir_Interaction, "vertex_evis_1Track", plotDir, 1288.07);
    DrawDataMC_Signal(rootDir_Interaction, "vertex_evis_2Track", plotDir, 1709.03);
    DrawDataMC_Signal(rootDir_Interaction, "vertex_evis_All", plotDir, 2997.1);

    DrawDataMCSignal_Diff(rootDir_Interaction, "vertex_energy_1Track", plotDir, 1288.07);
    DrawDataMCSignal_Diff(rootDir_Interaction, "vertex_energy_2Track", plotDir, 1709.03);
    DrawDataMCSignal_Diff(rootDir_Interaction, "vertex_energy_All", plotDir, 2997.1);

    DrawDataMCSignal_Diff(rootDir_Interaction, "vertex_evis_1Track", plotDir, 1288.07);
    DrawDataMCSignal_Diff(rootDir_Interaction, "vertex_evis_2Track", plotDir, 1709.03);
    DrawDataMCSignal_Diff(rootDir_Interaction, "vertex_evis_All", plotDir, 2997.1);

    DrawMnvH2D(rootDir_Interaction.data, "q3_q0_All_0", plotDir);
    DrawMnvH2D(rootDir_Interaction.mc, "q3_q0_All_0", plotDir);

    DrawMnvH2D(rootDir_Interaction.data, "W_QSq_All_0", plotDir);
    DrawMnvH2D(rootDir_Interaction.mc, "W_QSq_All_0", plotDir);

    DrawDataMCSignal_Diff_2D(rootDir_Interaction, "q3_q0_1Track", plotDir, 1288.07);
    DrawDataMCSignal_Diff_2D(rootDir_Interaction, "q3_q0_2Track", plotDir, 1709.03);
    DrawDataMCSignal_Diff_2D(rootDir_Interaction, "q3_q0_All", plotDir, 2997.1);

    DrawDataMCSignal_Diff_2D(rootDir_Interaction, "W_QSq_1Track", plotDir, 1288.07);
    DrawDataMCSignal_Diff_2D(rootDir_Interaction, "W_QSq_2Track", plotDir, 1709.03);
    DrawDataMCSignal_Diff_2D(rootDir_Interaction, "W_QSq_All", plotDir, 2997.1);
}

#endif

