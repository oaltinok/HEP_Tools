#include "CCProtonPi0_Plotter.h"

using namespace PlotUtils;

void CCProtonPi0_Plotter::GENIE_Tuning_Study()
{
    XSecVars_GENIE_Tuning_Ratios();
}

// Returns "new" TH1D's
void CCProtonPi0_Plotter::GENIE_Tuning_GetHistograms(std::string var_name, std::string data_var, std::string mc_var, TH1D* &data_nominal, TH1D* &data_tuned_test, TH1D* &data_tuned_best, TH1D* &mc_nominal, TH1D* &mc_tuned_test, TH1D* &mc_tuned_best)
{
    std::string data_dir_nominal = "/minerva/data/users/oaltinok/NTupleAnalysis_Nominal/Data/Analyzed/CrossSection.root";
    std::string data_dir_tuned_test = "/minerva/data/users/oaltinok/NTupleAnalysis_AllNonRES/Data/Analyzed/CrossSection.root";
    std::string data_dir_tuned_best = "/minerva/data/users/oaltinok/NTupleAnalysis_Best/Data/Analyzed/CrossSection.root";

    std::string mc_dir_nominal = "/minerva/data/users/oaltinok/NTupleAnalysis_Nominal/MC/Analyzed/CrossSection.root";
    std::string mc_dir_tuned_test = "/minerva/data/users/oaltinok/NTupleAnalysis_AllNonRES/MC/Analyzed/CrossSection.root";
    std::string mc_dir_tuned_best = "/minerva/data/users/oaltinok/NTupleAnalysis_Best/MC/Analyzed/CrossSection.root";

    TFile* f_data_nominal = new TFile(data_dir_nominal.c_str());
    TFile* f_data_tuned_test = new TFile(data_dir_tuned_test.c_str());
    TFile* f_data_tuned_best = new TFile(data_dir_tuned_best.c_str());

    TFile* f_mc_nominal = new TFile(mc_dir_nominal.c_str());
    TFile* f_mc_tuned_test = new TFile(mc_dir_tuned_test.c_str());
    TFile* f_mc_tuned_best = new TFile(mc_dir_tuned_best.c_str());

    MnvH1D* temp = NULL;

    // Get Data Histogram Before Tuning
    temp = GetMnvH1D(f_data_nominal, var_name + "_" + data_var);
    data_nominal = GetBinNormalizedTH1D(temp);
    delete temp;

    // Get Data Histogram After Tuning
    temp = GetMnvH1D(f_data_tuned_test, var_name + "_" + data_var);
    data_tuned_test = GetBinNormalizedTH1D(temp);
    delete temp;

    // Get Data Histogram After Tuning
    temp = GetMnvH1D(f_data_tuned_best, var_name + "_" + data_var);
    data_tuned_best = GetBinNormalizedTH1D(temp);
    delete temp;

    // Get MC Histogram Before Tuning
    temp = GetMnvH1D(f_mc_nominal, var_name + "_" + mc_var);
    mc_nominal = GetBinNormalizedTH1D(temp);
    delete temp;

    // Get MC Histogram After Tuning
    temp = GetMnvH1D(f_mc_tuned_test, var_name + "_" + mc_var);
    mc_tuned_test = GetBinNormalizedTH1D(temp);
    delete temp;

    // Get MC Histogram After Tuning
    temp = GetMnvH1D(f_mc_tuned_best, var_name + "_" + mc_var);
    mc_tuned_best = GetBinNormalizedTH1D(temp);
    delete temp;

    delete f_data_nominal;
    delete f_data_tuned_test;
    delete f_data_tuned_best;
    delete f_mc_nominal;
    delete f_mc_tuned_test;
    delete f_mc_tuned_best;
}

void CCProtonPi0_Plotter::GENIE_Tuning_DataMC_Ratio(std::string var_name, std::string data_var, std::string mc_var)
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    TH1D* data_nominal = NULL;
    TH1D* data_tuned_test = NULL;
    TH1D* data_tuned_best = NULL;

    TH1D* mc_nominal = NULL;
    TH1D* mc_tuned_test = NULL;
    TH1D* mc_tuned_best = NULL;

    GENIE_Tuning_GetHistograms(var_name, data_var, mc_var, data_nominal, data_tuned_test, data_tuned_best, mc_nominal, mc_tuned_test, mc_tuned_best);

    // Get Ratios
    TH1D* ratio_nominal = (TH1D*) data_nominal->Clone();
    TH1D* ratio_test = (TH1D*) data_tuned_test->Clone();
    TH1D* ratio_best = (TH1D*) data_tuned_best->Clone();

    // Divide Histograms
    ratio_nominal->Divide(mc_nominal);
    ratio_test->Divide(mc_tuned_test);
    ratio_best->Divide(mc_tuned_best);
    
    delete data_nominal;
    delete data_tuned_test;
    delete data_tuned_best;
    delete mc_nominal;
    delete mc_tuned_test;
    delete mc_tuned_best;

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

    ratio_test->GetYaxis()->SetTitle("Ratio");
    ratio_test->SetMinimum(0.5);
    ratio_test->SetMaximum(1.5);
    ratio_test->SetLineColor(kRed);
    ratio_test->SetLineWidth(3);
    ratio_test->SetFillColor(kWhite);

    ratio_best->GetYaxis()->SetTitle("Ratio");
    ratio_best->SetMinimum(0.5);
    ratio_best->SetMaximum(1.5);
    ratio_best->SetLineColor(kBlue);
    ratio_best->SetLineWidth(3);
    ratio_best->SetFillColor(kWhite);


    // Plot
    ratio_nominal->Draw("HIST");
    ratio_test->Draw("HIST SAME");
    ratio_best->Draw("HIST SAME");
    gPad->Update();
    gStyle->SetOptStat(0); 
    c->Update();

    // Add Ratio = 1 Line 
    TLine ratio_1;
    ratio_1.SetLineWidth(2);
    ratio_1.SetLineStyle(7);
    ratio_1.SetLineColor(kGreen+3);
    double line_min = ratio_test->GetBinLowEdge(1);
    double line_max = ratio_test->GetBinLowEdge(ratio_test->GetNbinsX()+1);
    ratio_1.DrawLine(line_min,1,line_max,1);

    TLegend *legend = new TLegend(0.4,0.7,0.9,0.9);  
    legend->AddEntry(ratio_nominal, "Nominal(Data/GENIE)", "l" );
    legend->AddEntry(ratio_test, "All NonRES1pi(Data/GENIE)", "l" );
    legend->AddEntry(ratio_best, "Best(Data/GENIE)", "l" );
    legend->SetTextSize(0.04);
    legend->Draw();
 
    var_name = var_name + "_DataMC_ratio";
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".png"), "png");

    delete legend;
    delete c;
    delete ratio_nominal;
    delete ratio_test;
    delete ratio_best;
}

void CCProtonPi0_Plotter::GENIE_Tuning_Ratio(std::string var_name, std::string data_var, std::string mc_var)
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    TH1D* data_nominal = NULL;
    TH1D* data_tuned_test = NULL;
    TH1D* data_tuned_best = NULL;

    TH1D* mc_nominal = NULL;
    TH1D* mc_tuned_test = NULL;
    TH1D* mc_tuned_best = NULL;


    GENIE_Tuning_GetHistograms(var_name, data_var, mc_var, data_nominal, data_tuned_test, data_tuned_best, mc_nominal, mc_tuned_test, mc_tuned_best);

    // Get Ratios
    TH1D* mc_ratio_test = (TH1D*) mc_tuned_test->Clone();
    TH1D* mc_ratio_best = (TH1D*) mc_tuned_best->Clone();

    TH1D* data_ratio_test = (TH1D*) data_tuned_test->Clone();
    TH1D* data_ratio_best = (TH1D*) data_tuned_best->Clone();

    // Divide Histograms
    mc_ratio_test->Divide(mc_nominal);
    mc_ratio_best->Divide(mc_nominal);
    data_ratio_test->Divide(data_nominal);
    data_ratio_best->Divide(data_nominal);
    delete data_nominal;
    delete data_tuned_test;
    delete data_tuned_best;
    delete mc_nominal;
    delete mc_tuned_test;
    delete mc_tuned_best;

    // Create Canvas
    TCanvas* c = new TCanvas("c","c",1280,800);
    TPad* pad = (TPad*)c->cd();
    pad->cd();
    pad->SetTopMargin(0.30);

    // Style Histograms
    mc_ratio_test->GetYaxis()->SetTitle("Ratio");
    mc_ratio_test->SetMinimum(0.5);
    mc_ratio_test->SetMaximum(1.5);
    mc_ratio_test->SetLineColor(kRed);
    mc_ratio_test->SetLineWidth(3);
    mc_ratio_test->SetFillColor(kWhite);

    mc_ratio_best->SetMinimum(0.5);
    mc_ratio_best->SetMaximum(1.5);
    mc_ratio_best->SetLineColor(kRed);
    mc_ratio_best->SetLineWidth(3);
    mc_ratio_best->SetLineStyle(2);
    mc_ratio_best->SetFillColor(kWhite);

    data_ratio_test->SetMinimum(0.5);
    data_ratio_test->SetMaximum(1.5);
    data_ratio_test->SetLineColor(kBlack);
    data_ratio_test->SetLineWidth(3);
    data_ratio_test->SetFillColor(kWhite);

    data_ratio_best->SetMinimum(0.5);
    data_ratio_best->SetMaximum(1.5);
    data_ratio_best->SetLineColor(kBlack);
    data_ratio_best->SetLineWidth(3);
    data_ratio_best->SetLineStyle(2);
    data_ratio_best->SetFillColor(kWhite);

    gStyle->SetEndErrorSize(6);

    // Plot
    mc_ratio_test->Draw("HIST ");
    mc_ratio_best->Draw("HIST SAME");
    data_ratio_test->Draw("HIST SAME");
    data_ratio_best->Draw("HIST SAME");
    gPad->Update();
    gStyle->SetOptStat(0); 
    c->Update();

    // Add Ratio = 1 Line 
    TLine ratio_1;
    ratio_1.SetLineWidth(2);
    ratio_1.SetLineStyle(7);
    ratio_1.SetLineColor(kGreen+3);
    double line_min = mc_ratio_test->GetBinLowEdge(1);
    double line_max = mc_ratio_test->GetBinLowEdge(mc_ratio_test->GetNbinsX()+1);
    ratio_1.DrawLine(line_min,1,line_max,1);

    TLegend *legend = new TLegend(0.4,0.7,0.9,0.9);  
    legend->AddEntry(mc_ratio_test, "GENIE(All NonRES1pi) / GENIE", "l");
    legend->AddEntry(mc_ratio_best, "GENIE(Best) / GENIE", "l");
    legend->AddEntry(data_ratio_test, "Data(All NonRES1pi) / Data", "l");
    legend->AddEntry(data_ratio_best, "Data(Best) / Data", "l");
    legend->SetTextSize(0.04);
    legend->Draw();
 
    var_name = var_name + "_" + data_var + "_ratio";
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".png"), "png");

    delete c;
    delete legend;
    delete mc_ratio_test;
    delete mc_ratio_best;
    delete data_ratio_test;
    delete data_ratio_best;
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

