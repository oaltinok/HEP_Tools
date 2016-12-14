#include "CCProtonPi0_Plotter.h"

using namespace PlotUtils;

void CCProtonPi0_Plotter::GENIE_Tuning_Study()
{
    XSecVars_GENIE_Tuning_Ratios();
}

// Returns "new" TH1D's
void CCProtonPi0_Plotter::GENIE_Tuning_GetHistograms(std::string var_name, std::string data_var, std::string mc_var, TH1D* &data_nominal, TH1D* &data_tuned_v3, TH1D* &data_tuned_v4, TH1D* &mc_nominal, TH1D* &mc_tuned_v3, TH1D* &mc_tuned_v4)
{
    std::string data_dir_nominal = "/minerva/data/users/oaltinok/NTupleAnalysis_GENIE_Nominal/Data/Analyzed/CrossSection.root";
    std::string data_dir_tuned_v3 = "/minerva/data/users/oaltinok/NTupleAnalysis_GENIE_Tuning_v3/Data/Analyzed/CrossSection.root";
    std::string data_dir_tuned_v4 = "/minerva/data/users/oaltinok/NTupleAnalysis_GENIE_Tuning_v4/Data/Analyzed/CrossSection.root";

    std::string mc_dir_nominal = "/minerva/data/users/oaltinok/NTupleAnalysis_GENIE_Nominal/MC/Analyzed/CrossSection.root";
    std::string mc_dir_tuned_v3 = "/minerva/data/users/oaltinok/NTupleAnalysis_GENIE_Tuning_v3/MC/Analyzed/CrossSection.root";
    std::string mc_dir_tuned_v4 = "/minerva/data/users/oaltinok/NTupleAnalysis_GENIE_Tuning_v4/MC/Analyzed/CrossSection.root";

    TFile* f_data_nominal = new TFile(data_dir_nominal.c_str());
    TFile* f_data_tuned_v3 = new TFile(data_dir_tuned_v3.c_str());
    TFile* f_data_tuned_v4 = new TFile(data_dir_tuned_v4.c_str());

    TFile* f_mc_nominal = new TFile(mc_dir_nominal.c_str());
    TFile* f_mc_tuned_v3 = new TFile(mc_dir_tuned_v3.c_str());
    TFile* f_mc_tuned_v4 = new TFile(mc_dir_tuned_v4.c_str());

    MnvH1D* temp = NULL;

    // Get Data Histogram Before Tuning
    temp = GetMnvH1D(f_data_nominal, var_name + "_" + data_var);
    data_nominal = GetBinNormalizedTH1D(temp);
    delete temp;

    // Get Data Histogram After Tuning
    temp = GetMnvH1D(f_data_tuned_v3, var_name + "_" + data_var);
    data_tuned_v3 = GetBinNormalizedTH1D(temp);
    delete temp;

    // Get Data Histogram After Tuning
    temp = GetMnvH1D(f_data_tuned_v4, var_name + "_" + data_var);
    data_tuned_v4 = GetBinNormalizedTH1D(temp);
    delete temp;

    // Get MC Histogram Before Tuning
    temp = GetMnvH1D(f_mc_nominal, var_name + "_" + mc_var);
    mc_nominal = GetBinNormalizedTH1D(temp);
    delete temp;

    // Get MC Histogram After Tuning
    temp = GetMnvH1D(f_mc_tuned_v3, var_name + "_" + mc_var);
    mc_tuned_v3 = GetBinNormalizedTH1D(temp);
    delete temp;

    // Get MC Histogram After Tuning
    temp = GetMnvH1D(f_mc_tuned_v4, var_name + "_" + mc_var);
    mc_tuned_v4 = GetBinNormalizedTH1D(temp);
    delete temp;

    delete f_data_nominal;
    delete f_data_tuned_v3;
    delete f_data_tuned_v4;
    delete f_mc_nominal;
    delete f_mc_tuned_v3;
    delete f_mc_tuned_v4;
}

void CCProtonPi0_Plotter::GENIE_Tuning_DataMC_Ratio(std::string var_name, std::string data_var, std::string mc_var)
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    TH1D* data_nominal = NULL;
    TH1D* data_tuned_v3 = NULL;
    TH1D* data_tuned_v4 = NULL;

    TH1D* mc_nominal = NULL;
    TH1D* mc_tuned_v3 = NULL;
    TH1D* mc_tuned_v4 = NULL;

    GENIE_Tuning_GetHistograms(var_name, data_var, mc_var, data_nominal, data_tuned_v3, data_tuned_v4, mc_nominal, mc_tuned_v3, mc_tuned_v4);

    // Get Ratios
    TH1D* ratio_nominal = (TH1D*) data_nominal->Clone();
    TH1D* ratio_v3 = (TH1D*) data_tuned_v3->Clone();
    TH1D* ratio_v4 = (TH1D*) data_tuned_v4->Clone();

    // Divide Histograms
    ratio_nominal->Divide(mc_nominal);
    ratio_v3->Divide(mc_tuned_v3);
    ratio_v4->Divide(mc_tuned_v4);
    
    delete data_nominal;
    delete data_tuned_v3;
    delete data_tuned_v4;
    delete mc_nominal;
    delete mc_tuned_v3;
    delete mc_tuned_v4;

    // Create Canvas
    TCanvas* c = new TCanvas("c","c",1280,800);
    TPad* pad = (TPad*)c->cd();
    pad->cd();
    pad->SetTopMargin(0.30);

    // Style Histograms
    ratio_nominal->GetYaxis()->SetTitle("Ratio");
    ratio_nominal->SetMinimum(0.5);
    ratio_nominal->SetMaximum(1.5);
    ratio_nominal->SetLineColor(kBlack);
    ratio_nominal->SetLineWidth(3);
    ratio_nominal->SetFillColor(kWhite);

    ratio_v3->GetYaxis()->SetTitle("Ratio");
    ratio_v3->SetMinimum(0.5);
    ratio_v3->SetMaximum(1.5);
    ratio_v3->SetLineColor(kRed);
    ratio_v3->SetLineWidth(3);
    ratio_v3->SetFillColor(kWhite);

    ratio_v4->GetYaxis()->SetTitle("Ratio");
    ratio_v4->SetMinimum(0.5);
    ratio_v4->SetMaximum(1.5);
    ratio_v4->SetLineColor(kBlue);
    ratio_v4->SetLineWidth(3);
    ratio_v4->SetFillColor(kWhite);


    // Plot
    ratio_nominal->Draw("HIST");
    ratio_v3->Draw("HIST SAME");
    ratio_v4->Draw("HIST SAME");
    gPad->Update();
    gStyle->SetOptStat(0); 
    c->Update();

    // Add Ratio = 1 Line 
    TLine ratio_1;
    ratio_1.SetLineWidth(2);
    ratio_1.SetLineStyle(7);
    ratio_1.SetLineColor(kGreen+3);
    double line_min = ratio_v3->GetBinLowEdge(1);
    double line_max = ratio_v3->GetBinLowEdge(ratio_v3->GetNbinsX()+1);
    ratio_1.DrawLine(line_min,1,line_max,1);

    TLegend *legend = new TLegend(0.4,0.7,0.9,0.9);  
    legend->AddEntry(ratio_nominal, "Nominal(Data/GENIE)", "l" );
    legend->AddEntry(ratio_v3, "Tuned+2p2h(Data/GENIE)", "l" );
    legend->AddEntry(ratio_v4, "Only 2p2h(Data/GENIE)", "l" );
    legend->SetTextSize(0.04);
    legend->Draw();
 
    var_name = var_name + "_DataMC_ratio";
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".png"), "png");

    delete legend;
    delete c;
    delete ratio_nominal;
    delete ratio_v3;
    delete ratio_v4;
}

void CCProtonPi0_Plotter::GENIE_Tuning_Ratio(std::string var_name, std::string data_var, std::string mc_var)
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    TH1D* data_nominal = NULL;
    TH1D* data_tuned_v3 = NULL;
    TH1D* data_tuned_v4 = NULL;

    TH1D* mc_nominal = NULL;
    TH1D* mc_tuned_v3 = NULL;
    TH1D* mc_tuned_v4 = NULL;


    GENIE_Tuning_GetHistograms(var_name, data_var, mc_var, data_nominal, data_tuned_v3, data_tuned_v4, mc_nominal, mc_tuned_v3, mc_tuned_v4);

    // Get Ratios
    TH1D* mc_ratio_v3 = (TH1D*) mc_tuned_v3->Clone();
    TH1D* mc_ratio_v4 = (TH1D*) mc_tuned_v4->Clone();

    TH1D* data_ratio_v3 = (TH1D*) data_tuned_v3->Clone();
    TH1D* data_ratio_v4 = (TH1D*) data_tuned_v4->Clone();

    // Divide Histograms
    mc_ratio_v3->Divide(mc_nominal);
    mc_ratio_v4->Divide(mc_nominal);
    data_ratio_v3->Divide(data_nominal);
    data_ratio_v4->Divide(data_nominal);
    delete data_nominal;
    delete data_tuned_v3;
    delete data_tuned_v4;
    delete mc_nominal;
    delete mc_tuned_v3;
    delete mc_tuned_v4;

    // Create Canvas
    TCanvas* c = new TCanvas("c","c",1280,800);
    TPad* pad = (TPad*)c->cd();
    pad->cd();
    pad->SetTopMargin(0.30);

    // Style Histograms
    mc_ratio_v3->GetYaxis()->SetTitle("Ratio");
    mc_ratio_v3->SetMinimum(0.5);
    mc_ratio_v3->SetMaximum(1.5);
    mc_ratio_v3->SetLineColor(kRed);
    mc_ratio_v3->SetLineWidth(3);
    mc_ratio_v3->SetFillColor(kWhite);

    mc_ratio_v4->SetMinimum(0.5);
    mc_ratio_v4->SetMaximum(1.5);
    mc_ratio_v4->SetLineColor(kRed);
    mc_ratio_v4->SetLineWidth(3);
    mc_ratio_v4->SetLineStyle(2);
    mc_ratio_v4->SetFillColor(kWhite);

    data_ratio_v3->SetMinimum(0.5);
    data_ratio_v3->SetMaximum(1.5);
    data_ratio_v3->SetLineColor(kBlack);
    data_ratio_v3->SetLineWidth(3);
    data_ratio_v3->SetFillColor(kWhite);

    data_ratio_v4->SetMinimum(0.5);
    data_ratio_v4->SetMaximum(1.5);
    data_ratio_v4->SetLineColor(kBlack);
    data_ratio_v4->SetLineWidth(3);
    data_ratio_v4->SetLineStyle(2);
    data_ratio_v4->SetFillColor(kWhite);

    gStyle->SetEndErrorSize(6);

    // Plot
    mc_ratio_v3->Draw("HIST ");
    mc_ratio_v4->Draw("HIST SAME");
    data_ratio_v3->Draw("HIST SAME");
    data_ratio_v4->Draw("HIST SAME");
    gPad->Update();
    gStyle->SetOptStat(0); 
    c->Update();

    // Add Ratio = 1 Line 
    TLine ratio_1;
    ratio_1.SetLineWidth(2);
    ratio_1.SetLineStyle(7);
    ratio_1.SetLineColor(kGreen+3);
    double line_min = mc_ratio_v3->GetBinLowEdge(1);
    double line_max = mc_ratio_v3->GetBinLowEdge(mc_ratio_v3->GetNbinsX()+1);
    ratio_1.DrawLine(line_min,1,line_max,1);

    TLegend *legend = new TLegend(0.4,0.7,0.9,0.9);  
    legend->AddEntry(mc_ratio_v3, "GENIE(Tuned+2p2h) / GENIE", "l");
    legend->AddEntry(mc_ratio_v4, "GENIE(Only 2p2h) / GENIE", "l");
    legend->AddEntry(data_ratio_v3, "Data(Tuned+2p2h) / Data", "l");
    legend->AddEntry(data_ratio_v4, "Data(Only 2p2h) / Data", "l");
    legend->SetTextSize(0.04);
    legend->Draw();
 
    var_name = var_name + "_" + data_var + "_ratio";
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".png"), "png");

    delete c;
    delete legend;
    delete mc_ratio_v3;
    delete mc_ratio_v4;
    delete data_ratio_v3;
    delete data_ratio_v4;
}

void CCProtonPi0_Plotter::XSecVars_GENIE_Tuning_Ratios()
{
//    GENIE_Tuning_Ratio("muon_P", "all", "mc_reco_all");
    GENIE_Tuning_Ratio("muon_theta", "all", "mc_reco_all");
//    GENIE_Tuning_Ratio("pi0_P", "all", "mc_reco_all");
//    GENIE_Tuning_Ratio("pi0_KE", "all", "mc_reco_all");
//    GENIE_Tuning_Ratio("pi0_theta", "all", "mc_reco_all");
    //GENIE_Tuning_Ratio("QSq", "all", "mc_reco_all");
//    GENIE_Tuning_Ratio("W", "all", "mc_reco_all");
//    GENIE_Tuning_Ratio("Enu", "all", "mc_reco_all");
//
//    GENIE_Tuning_Ratio("muon_P", "bckg_subtracted", "bckg_subtracted");
    GENIE_Tuning_Ratio("muon_theta", "bckg_subtracted", "bckg_subtracted");
//    GENIE_Tuning_Ratio("pi0_P", "bckg_subtracted", "bckg_subtracted");
//    GENIE_Tuning_Ratio("pi0_KE", "bckg_subtracted", "bckg_subtracted");
//    GENIE_Tuning_Ratio("pi0_theta", "bckg_subtracted", "bckg_subtracted");
    //GENIE_Tuning_Ratio("QSq", "bckg_subtracted", "bckg_subtracted");
//    GENIE_Tuning_Ratio("W", "bckg_subtracted", "bckg_subtracted");
//    GENIE_Tuning_Ratio("Enu", "bckg_subtracted", "bckg_subtracted");

//    GENIE_Tuning_Ratio("muon_P", "bckg_estimated", "bckg_estimated");
    GENIE_Tuning_Ratio("muon_theta", "bckg_estimated", "bckg_estimated");
//    GENIE_Tuning_Ratio("pi0_P", "bckg_estimated", "bckg_estimated");
//    GENIE_Tuning_Ratio("pi0_KE", "bckg_estimated", "bckg_estimated");
//    GENIE_Tuning_Ratio("pi0_theta", "bckg_estimated", "bckg_estimated");
    //GENIE_Tuning_Ratio("QSq", "bckg_estimated", "bckg_estimated");
//    GENIE_Tuning_Ratio("W", "bckg_estimated", "bckg_estimated");
//    GENIE_Tuning_Ratio("Enu", "bckg_estimated", "bckg_estimated");

//    GENIE_Tuning_Ratio("muon_P", "unfolded", "unfolded");
    GENIE_Tuning_Ratio("muon_theta", "unfolded", "unfolded");
//    GENIE_Tuning_Ratio("pi0_P", "unfolded", "unfolded");
//    GENIE_Tuning_Ratio("pi0_KE", "unfolded", "unfolded");
//    GENIE_Tuning_Ratio("pi0_theta", "unfolded", "unfolded");
    //GENIE_Tuning_Ratio("QSq", "unfolded", "unfolded");
//    GENIE_Tuning_Ratio("W", "unfolded", "unfolded");
//    GENIE_Tuning_Ratio("Enu", "unfolded", "unfolded");
//
//    GENIE_Tuning_Ratio("muon_P", "efficiency_corrected", "efficiency_corrected");
    GENIE_Tuning_Ratio("muon_theta", "efficiency_corrected", "efficiency_corrected");
//    GENIE_Tuning_Ratio("pi0_P", "efficiency_corrected", "efficiency_corrected");
//    GENIE_Tuning_Ratio("pi0_KE", "efficiency_corrected", "efficiency_corrected");
//    GENIE_Tuning_Ratio("pi0_theta", "efficiency_corrected", "efficiency_corrected");
    //GENIE_Tuning_Ratio("QSq", "efficiency_corrected", "efficiency_corrected");
//    GENIE_Tuning_Ratio("W", "efficiency_corrected", "efficiency_corrected");
//    GENIE_Tuning_Ratio("Enu", "efficiency_corrected", "efficiency_corrected");

    GENIE_Tuning_Ratio("muon_P", "xsec", "xsec");
    GENIE_Tuning_Ratio("muon_theta", "xsec", "xsec");
    GENIE_Tuning_Ratio("pi0_P", "xsec", "xsec");
    GENIE_Tuning_Ratio("pi0_KE", "xsec", "xsec");
    GENIE_Tuning_Ratio("pi0_theta", "xsec", "xsec");
    GENIE_Tuning_Ratio("QSq", "xsec", "xsec");
    GENIE_Tuning_Ratio("W", "xsec", "xsec");
    GENIE_Tuning_Ratio("Enu", "xsec", "xsec");

    GENIE_Tuning_DataMC_Ratio("muon_P", "xsec", "xsec");
    GENIE_Tuning_DataMC_Ratio("muon_theta", "xsec", "xsec");
    GENIE_Tuning_DataMC_Ratio("pi0_P", "xsec", "xsec");
    GENIE_Tuning_DataMC_Ratio("pi0_KE", "xsec", "xsec");
    GENIE_Tuning_DataMC_Ratio("pi0_theta", "xsec", "xsec");
    GENIE_Tuning_DataMC_Ratio("QSq", "xsec", "xsec");
    GENIE_Tuning_DataMC_Ratio("W", "xsec", "xsec");
    GENIE_Tuning_DataMC_Ratio("Enu", "xsec", "xsec");
}

