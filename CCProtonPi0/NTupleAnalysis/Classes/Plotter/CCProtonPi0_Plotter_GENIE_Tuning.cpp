#include "CCProtonPi0_Plotter.h"

using namespace PlotUtils;

void CCProtonPi0_Plotter::GENIE_Tuning_Study()
{
    XSecVars_GENIE_Tuning_Ratios();
}

// Returns "new" TH1D's
void CCProtonPi0_Plotter::GENIE_Tuning_GetHistograms(std::string var_name, std::string data_var, std::string mc_var, TH1D* &data_nominal, TH1D* &data_tuned, TH1D* &mc_nominal, TH1D* &mc_tuned )
{
    std::string data_dir_before = "/minerva/data/users/oaltinok/NTupleAnalysis_GENIE_Nominal/Data/Analyzed/CrossSection.root";
    std::string data_dir_after = "/minerva/data/users/oaltinok/NTupleAnalysis_GENIE_Tuning_v2/Data/Analyzed/CrossSection.root";
    std::string mc_dir_before = "/minerva/data/users/oaltinok/NTupleAnalysis_GENIE_Nominal/MC/Analyzed/CrossSection.root";
    std::string mc_dir_after = "/minerva/data/users/oaltinok/NTupleAnalysis_GENIE_Tuning_v2/MC/Analyzed/CrossSection.root";

    TFile* f_data_before = new TFile(data_dir_before.c_str());
    TFile* f_data_after = new TFile(data_dir_after.c_str());
    TFile* f_mc_before = new TFile(mc_dir_before.c_str());
    TFile* f_mc_after = new TFile(mc_dir_after.c_str());

    // Get Data Histogram Before Tuning
    MnvH1D* temp = GetMnvH1D(f_data_before, var_name + "_" + data_var);
    data_nominal = GetBinNormalizedTH1D(temp);
    delete temp;

    // Get Data Histogram After Tuning
    temp = GetMnvH1D(f_data_after, var_name + "_" + data_var);
    data_tuned = GetBinNormalizedTH1D(temp);
    delete temp;

    // Get MC Histogram Before Tuning
    temp = GetMnvH1D(f_mc_before, var_name + "_" + mc_var);
    mc_nominal = GetBinNormalizedTH1D(temp);
    delete temp;

    // Get MC Histogram After Tuning
    temp = GetMnvH1D(f_mc_after, var_name + "_" + mc_var);
    mc_tuned = GetBinNormalizedTH1D(temp);
    delete temp;

    delete f_data_before;
    delete f_data_after;
    delete f_mc_before;
    delete f_mc_after;
}

void CCProtonPi0_Plotter::GENIE_Tuning_DataMC_Ratio(std::string var_name, std::string data_var, std::string mc_var)
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    TH1D* data_nominal = NULL;
    TH1D* data_tuned = NULL;
    TH1D* mc_nominal = NULL;
    TH1D* mc_tuned = NULL;

    GENIE_Tuning_GetHistograms(var_name, data_var, mc_var, data_nominal, data_tuned, mc_nominal, mc_tuned);

    // Get Ratios
    TH1D* tuned_data_mc = (TH1D*) data_tuned->Clone();
    TH1D* nominal_data_mc = (TH1D*) data_nominal->Clone();

    // Divide Histograms
    tuned_data_mc->Divide(mc_tuned);
    nominal_data_mc->Divide(mc_nominal);
    delete data_nominal;
    delete data_tuned;
    delete mc_nominal;
    delete mc_tuned;

    // Create Canvas
    TCanvas* c = new TCanvas("c","c",1280,800);

    // Style Histograms
    tuned_data_mc->GetYaxis()->SetTitle("Ratio");
    tuned_data_mc->SetMinimum(0.0);
    tuned_data_mc->SetMaximum(2.0);
    tuned_data_mc->SetLineColor(kRed);
    tuned_data_mc->SetLineWidth(3);
    tuned_data_mc->SetFillColor(kWhite);

    nominal_data_mc->SetMinimum(0.0);
    nominal_data_mc->SetMaximum(2.0);
    nominal_data_mc->SetLineColor(kBlack);
    nominal_data_mc->SetLineWidth(3);
    nominal_data_mc->SetFillColor(kWhite);

    // Plot
    tuned_data_mc->Draw("HIST");
    nominal_data_mc->Draw("HIST SAME");
    gPad->Update();
    gStyle->SetOptStat(0); 
    c->Update();

    // Add Ratio = 1 Line 
    TLine ratio_1;
    ratio_1.SetLineWidth(2);
    ratio_1.SetLineStyle(7);
    ratio_1.SetLineColor(kGreen+3);
    double line_min = tuned_data_mc->GetBinLowEdge(1);
    double line_max = tuned_data_mc->GetBinLowEdge(tuned_data_mc->GetNbinsX()+1);
    ratio_1.DrawLine(line_min,1,line_max,1);

    TLegend *legend = new TLegend(0.6,0.7,0.9,0.9);  
    legend->AddEntry(tuned_data_mc, "Data^{#diamond} / GENIE^{#diamond}", "l" );
    legend->AddEntry(nominal_data_mc, "Data / GENIE", "l" );
    legend->SetTextSize(0.04);
    legend->Draw();
 
    var_name = var_name + "_DataMC_ratio";
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".png"), "png");

    delete legend;
    delete c;
    delete tuned_data_mc;
    delete nominal_data_mc;
}

void CCProtonPi0_Plotter::GENIE_Tuning_Ratio(std::string var_name, std::string data_var, std::string mc_var)
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    TH1D* data_nominal = NULL;
    TH1D* data_tuned = NULL;
    TH1D* mc_nominal = NULL;
    TH1D* mc_tuned = NULL;

    GENIE_Tuning_GetHistograms(var_name, data_var, mc_var, data_nominal, data_tuned, mc_nominal, mc_tuned);

    //std::cout<<"Ratios"<<std::endl;
    //std::cout<<var_name + "_" + data_var + " Nominal = "<<data_nominal->Integral()<<std::endl;
    //std::cout<<var_name + "_" + data_var + " Tuned = "<<data_tuned->Integral()<<std::endl;
    //std::cout<<var_name + "_" + mc_var + " Nominal = "<<mc_nominal->Integral()<<std::endl;
    //std::cout<<var_name + "_" + mc_var + " Tuned = "<<mc_tuned->Integral()<<std::endl;

    // Get Ratios
    TH1D* tuned_mc_ratio_nominal_mc = (TH1D*) mc_tuned->Clone();
    TH1D* tuned_data_ratio_nominal_data = (TH1D*) data_tuned->Clone();

    // Divide Histograms
    tuned_mc_ratio_nominal_mc->Divide(mc_nominal);
    tuned_data_ratio_nominal_data->Divide(data_nominal);
    delete data_nominal;
    delete data_tuned;
    delete mc_nominal;
    delete mc_tuned;

    // Create Canvas
    TCanvas* c = new TCanvas("c","c",1280,800);

    // Style Histograms
    tuned_mc_ratio_nominal_mc->GetYaxis()->SetTitle("Ratio");
    tuned_mc_ratio_nominal_mc->SetMinimum(0.0);
    tuned_mc_ratio_nominal_mc->SetMaximum(2.0);
    tuned_mc_ratio_nominal_mc->SetLineColor(kRed);
    tuned_mc_ratio_nominal_mc->SetLineWidth(3);
    tuned_mc_ratio_nominal_mc->SetFillColor(kWhite);

    tuned_data_ratio_nominal_data->SetMinimum(0.0);
    tuned_data_ratio_nominal_data->SetMaximum(2.0);
    tuned_data_ratio_nominal_data->SetLineColor(kBlack);
    tuned_data_ratio_nominal_data->SetLineWidth(3);
    tuned_data_ratio_nominal_data->SetFillColor(kWhite);

    // Plot
    tuned_mc_ratio_nominal_mc->Draw("HIST");
    tuned_data_ratio_nominal_data->Draw("HIST SAME");
    gPad->Update();
    gStyle->SetOptStat(0); 
    c->Update();

    // Add Ratio = 1 Line 
    TLine ratio_1;
    ratio_1.SetLineWidth(2);
    ratio_1.SetLineStyle(7);
    ratio_1.SetLineColor(kGreen+3);
    double line_min = tuned_mc_ratio_nominal_mc->GetBinLowEdge(1);
    double line_max = tuned_mc_ratio_nominal_mc->GetBinLowEdge(tuned_mc_ratio_nominal_mc->GetNbinsX()+1);
    ratio_1.DrawLine(line_min,1,line_max,1);

    TLegend *legend = new TLegend(0.6,0.7,0.9,0.9);  
    legend->AddEntry(tuned_mc_ratio_nominal_mc, "GENIE^{#diamond} / GENIE", "l");
    legend->AddEntry(tuned_data_ratio_nominal_data, "Data^{#diamond} / Data", "l");
    legend->SetTextSize(0.04);
    legend->Draw();
 
    var_name = var_name + "_" + data_var + "_ratio";
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".png"), "png");

    delete legend;
    delete c;
}

void CCProtonPi0_Plotter::XSecVars_GENIE_Tuning_Ratios()
{
    GENIE_Tuning_Ratio("muon_P", "all", "mc_reco_all");
    GENIE_Tuning_Ratio("muon_theta", "all", "mc_reco_all");
    GENIE_Tuning_Ratio("pi0_P", "all", "mc_reco_all");
    GENIE_Tuning_Ratio("pi0_KE", "all", "mc_reco_all");
    GENIE_Tuning_Ratio("pi0_theta", "all", "mc_reco_all");
    GENIE_Tuning_Ratio("QSq", "all", "mc_reco_all");
    GENIE_Tuning_Ratio("W", "all", "mc_reco_all");
    GENIE_Tuning_Ratio("Enu", "all", "mc_reco_all");

    GENIE_Tuning_Ratio("muon_P", "bckg_subtracted", "bckg_subtracted");
    GENIE_Tuning_Ratio("muon_theta", "bckg_subtracted", "bckg_subtracted");
    GENIE_Tuning_Ratio("pi0_P", "bckg_subtracted", "bckg_subtracted");
    GENIE_Tuning_Ratio("pi0_KE", "bckg_subtracted", "bckg_subtracted");
    GENIE_Tuning_Ratio("pi0_theta", "bckg_subtracted", "bckg_subtracted");
    GENIE_Tuning_Ratio("QSq", "bckg_subtracted", "bckg_subtracted");
    GENIE_Tuning_Ratio("W", "bckg_subtracted", "bckg_subtracted");
    GENIE_Tuning_Ratio("Enu", "bckg_subtracted", "bckg_subtracted");

    GENIE_Tuning_Ratio("muon_P", "unfolded", "unfolded");
    GENIE_Tuning_Ratio("muon_theta", "unfolded", "unfolded");
    GENIE_Tuning_Ratio("pi0_P", "unfolded", "unfolded");
    GENIE_Tuning_Ratio("pi0_KE", "unfolded", "unfolded");
    GENIE_Tuning_Ratio("pi0_theta", "unfolded", "unfolded");
    GENIE_Tuning_Ratio("QSq", "unfolded", "unfolded");
    GENIE_Tuning_Ratio("W", "unfolded", "unfolded");
    GENIE_Tuning_Ratio("Enu", "unfolded", "unfolded");

    GENIE_Tuning_Ratio("muon_P", "efficiency_corrected", "efficiency_corrected");
    GENIE_Tuning_Ratio("muon_theta", "efficiency_corrected", "efficiency_corrected");
    GENIE_Tuning_Ratio("pi0_P", "efficiency_corrected", "efficiency_corrected");
    GENIE_Tuning_Ratio("pi0_KE", "efficiency_corrected", "efficiency_corrected");
    GENIE_Tuning_Ratio("pi0_theta", "efficiency_corrected", "efficiency_corrected");
    GENIE_Tuning_Ratio("QSq", "efficiency_corrected", "efficiency_corrected");
    GENIE_Tuning_Ratio("W", "efficiency_corrected", "efficiency_corrected");
    GENIE_Tuning_Ratio("Enu", "efficiency_corrected", "efficiency_corrected");

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

