#include "CCProtonPi0_Plotter.h"

using namespace PlotUtils;

void CCProtonPi0_Plotter::GENIE_Tuning_Study()
{
    XSecVars_GENIE_Tuning_Ratios();
    //XSecVars_CVWeight();
}

void CCProtonPi0_Plotter::GENIE_Tuning_Ratio(std::string var_name)
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    std::string data_dir = "/minerva/data/users/oaltinok/NTupleAnalysis_After_GENIE_Tuning/Data/Analyzed/CrossSection.root";
    std::string mc_dir_before = "/minerva/data/users/oaltinok/NTupleAnalysis_Before_GENIE_Tuning/MC/Analyzed/CrossSection.root";
    std::string mc_dir_after = "/minerva/data/users/oaltinok/NTupleAnalysis_After_GENIE_Tuning/MC/Analyzed/CrossSection.root";

    TFile* f_data = new TFile(data_dir.c_str());
    TFile* f_mc_before = new TFile(mc_dir_before.c_str());
    TFile* f_mc_after = new TFile(mc_dir_after.c_str());

    // Get Data Histogram
    MnvH1D* temp = GetMnvH1D(f_data, var_name + "_all");
    TH1D* data_base = GetBinNormalizedTH1D(temp);
    delete temp;

    // Get Histogram Before Tuning
    temp = GetMnvH1D(f_mc_before, var_name + "_mc_reco_all");
    TH1D* mc_base = GetBinNormalizedTH1D(temp);
    delete temp;

    // Get Histogram After Tuning
    temp = GetMnvH1D(f_mc_after, var_name + "_mc_reco_all");
    TH1D* ratio_data = GetBinNormalizedTH1D(temp);
    TH1D* ratio_mc = GetBinNormalizedTH1D(temp);
    delete temp;

    // Scale MC
    ratio_data->Scale(POT_ratio);

    // Divide Histograms
    ratio_data->Divide(data_base);
    ratio_mc->Divide(mc_base);
    delete data_base;
    delete mc_base;

    // Create Canvas
    TCanvas* c = new TCanvas("c","c",1280,800);

    // Style Histograms
    ratio_data->GetYaxis()->SetTitle("Ratio");
    ratio_data->SetMinimum(0.5);
    ratio_data->SetMaximum(1.5);
    ratio_data->SetLineColor(kBlack);
    ratio_data->SetLineWidth(3);
    ratio_data->SetFillColor(kWhite);

    ratio_mc->SetMinimum(0.5);
    ratio_mc->SetMaximum(1.5);
    ratio_mc->SetLineColor(kRed);
    ratio_mc->SetLineWidth(3);
    ratio_mc->SetFillColor(kWhite);

    // Plot
    ratio_data->Draw("HIST");
    ratio_mc->Draw("HIST SAME");
    gPad->Update();
    gStyle->SetOptStat(0); 
    c->Update();

    // Add Ratio = 1 Line 
    TLine ratio_1;
    ratio_1.SetLineWidth(2);
    ratio_1.SetLineStyle(7);
    ratio_1.SetLineColor(kBlue);
    double line_min = ratio_data->GetBinLowEdge(1);
    double line_max = ratio_data->GetBinLowEdge(ratio_data->GetNbinsX()+1);
    ratio_1.DrawLine(line_min,1,line_max,1);

    TLegend *legend = new TLegend(0.6,0.8,0.9,0.9);  
    legend->AddEntry(ratio_data, "Data", "l" );
    legend->AddEntry(ratio_mc, "MC w/o GENIE Tuning","l");
    legend->SetTextSize(0.04);
    legend->Draw();
 
    var_name = var_name + "_ratio";
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".png"), "png");

    delete ratio_data;
    delete ratio_mc;
    delete legend;
    delete c;
    delete f_data;
    delete f_mc_before;
    delete f_mc_after;
}

void CCProtonPi0_Plotter::XSecVars_GENIE_Tuning_Ratios()
{
    GENIE_Tuning_Ratio("muon_P");
    GENIE_Tuning_Ratio("muon_theta");
    GENIE_Tuning_Ratio("pi0_P");
    GENIE_Tuning_Ratio("pi0_KE");
    GENIE_Tuning_Ratio("pi0_theta");
    GENIE_Tuning_Ratio("QSq");
    GENIE_Tuning_Ratio("Enu");
}

void CCProtonPi0_Plotter::XSecVars_CVWeight()
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;
//
//    DrawDataStackedMC(rootDir_Interaction,"muon_P_GENIE_Tuning",plotDir);
//    DrawDataStackedMC(rootDir_Interaction,"muon_P_NO_GENIE_Tuning",plotDir);
//    DrawDataStackedMC(rootDir_Interaction,"muon_P_Delta",plotDir);
//    DrawDataStackedMC(rootDir_Interaction,"muon_P_CCRES",plotDir);
//    DrawDataStackedMC(rootDir_Interaction,"muon_P_NonRes1pi",plotDir);
    
    //DrawDataStackedMC(rootDir_Interaction,"muon_theta_GENIE_Tuning",plotDir);
    //DrawDataStackedMC(rootDir_Interaction,"muon_theta_Delta",plotDir);
    //DrawDataStackedMC(rootDir_Interaction,"muon_theta_CCRES",plotDir);
    //DrawDataStackedMC(rootDir_Interaction,"muon_theta_NonRes1pi",plotDir);
 
//    DrawDataStackedMC(rootDir_Interaction,"pi0_P_GENIE_Tuning",plotDir);
//    DrawDataStackedMC(rootDir_Interaction,"pi0_P_NO_GENIE_Tuning",plotDir);
//    DrawDataStackedMC(rootDir_Interaction,"pi0_P_Delta",plotDir);
//    DrawDataStackedMC(rootDir_Interaction,"pi0_P_CCRES",plotDir);
//    DrawDataStackedMC(rootDir_Interaction,"pi0_P_NonRes1pi",plotDir);
//// 
  //  DrawDataStackedMC(rootDir_Interaction,"pi0_KE_GENIE_Tuning",plotDir);
  //  DrawDataStackedMC(rootDir_Interaction,"pi0_KE_NO_GENIE_Tuning",plotDir);
  //  DrawDataStackedMC(rootDir_Interaction,"pi0_KE_Delta",plotDir);
  //  DrawDataStackedMC(rootDir_Interaction,"pi0_KE_CCRES",plotDir);
  //  DrawDataStackedMC(rootDir_Interaction,"pi0_KE_NonRes1pi",plotDir);
 
//    DrawDataStackedMC(rootDir_Interaction,"pi0_theta_GENIE_Tuning",plotDir);
//    DrawDataStackedMC(rootDir_Interaction,"pi0_theta_NO_GENIE_Tuning",plotDir);
//    DrawDataStackedMC(rootDir_Interaction,"pi0_theta_Delta",plotDir);
//    DrawDataStackedMC(rootDir_Interaction,"pi0_theta_CCRES",plotDir);
//    DrawDataStackedMC(rootDir_Interaction,"pi0_theta_NonRes1pi",plotDir);
// 
//    DrawDataStackedMC(rootDir_Interaction,"QSq_GENIE_Tuning",plotDir);
//    DrawDataStackedMC(rootDir_Interaction,"QSq_NO_GENIE_Tuning",plotDir);
//    DrawDataStackedMC(rootDir_Interaction,"QSq_Delta",plotDir);
//    DrawDataStackedMC(rootDir_Interaction,"QSq_CCRES",plotDir);
//    DrawDataStackedMC(rootDir_Interaction,"QSq_NonRes1pi",plotDir);
// 
//    DrawDataStackedMC(rootDir_Interaction,"Enu_GENIE_Tuning",plotDir);
//    DrawDataStackedMC(rootDir_Interaction,"Enu_NO_GENIE_Tuning",plotDir);
//    DrawDataStackedMC(rootDir_Interaction,"Enu_Delta",plotDir);
//    DrawDataStackedMC(rootDir_Interaction,"Enu_CCRES",plotDir);
//    DrawDataStackedMC(rootDir_Interaction,"Enu_NonRes1pi",plotDir);
}
