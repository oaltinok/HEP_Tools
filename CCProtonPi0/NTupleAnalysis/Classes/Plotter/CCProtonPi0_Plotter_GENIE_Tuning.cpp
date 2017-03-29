#include "CCProtonPi0_Plotter.h"

using namespace PlotUtils;

void CCProtonPi0_Plotter::GENIE_Tuning_Study()
{
    //XSecVars_GENIE_Tuning_Ratios();

    Draw_Comparison_DeltaFactor();
    //Draw_Comparison_Nominal();
}

void CCProtonPi0_Plotter::Draw_Comparison_DeltaFactor()
{
    Draw_Comparison_DeltaFactor("muon_P_xsec");
    Draw_Comparison_DeltaFactor("muon_theta_xsec");
    Draw_Comparison_DeltaFactor("pi0_P_xsec");
    Draw_Comparison_DeltaFactor("pi0_KE_xsec");
    Draw_Comparison_DeltaFactor("pi0_theta_xsec");
    Draw_Comparison_DeltaFactor("QSq_xsec");
    Draw_Comparison_DeltaFactor("Enu_xsec");
}

void CCProtonPi0_Plotter::Draw_Comparison_Nominal()
{
    Draw_Comparison_Nominal("muon_P_xsec");
    Draw_Comparison_Nominal("muon_theta_xsec");
    Draw_Comparison_Nominal("pi0_P_xsec");
    Draw_Comparison_Nominal("pi0_KE_xsec");
    Draw_Comparison_Nominal("pi0_theta_xsec");
    Draw_Comparison_Nominal("QSq_xsec");
    Draw_Comparison_Nominal("Enu_xsec");
}

void CCProtonPi0_Plotter::Draw_Comparison_DeltaFactor(std::string var_name)
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    std::string final_data_dir = "/minerva/data/users/oaltinok/NTupleAnalysis_Final_XSecs/Data/Analyzed/CrossSection.root";
    std::string final_mc_dir = "/minerva/data/users/oaltinok/NTupleAnalysis_Final_XSecs/MC/Analyzed/CrossSection.root";

    std::string deltafactor_data_dir = "/minerva/data/users/oaltinok/NTupleAnalysis_DeltaFactor_Applied/Data/Analyzed/CrossSection.root";
    std::string deltafactor_mc_dir = "/minerva/data/users/oaltinok/NTupleAnalysis_DeltaFactor_Applied/MC/Analyzed/CrossSection.root";

    TFile* f_data_final = new TFile(final_data_dir.c_str());
    TFile* f_mc_final = new TFile(final_mc_dir.c_str());

    TFile* f_data_deltafactor = new TFile(deltafactor_data_dir.c_str());
    TFile* f_mc_deltafactor = new TFile(deltafactor_mc_dir.c_str());

    // --------------------------------------------------------------------
    // Get Histograms 
    // --------------------------------------------------------------------
    MnvH1D* data_deltafactor = GetMnvH1D(f_data_deltafactor, var_name);
    MnvH1D* mc_deltafactor = GetMnvH1D(f_mc_deltafactor, var_name);
  
    TH1D* h_data_deltafactor = GetBinNormalizedTH1D(data_deltafactor, true);
    TH1D* h_mc_deltafactor = GetBinNormalizedTH1D(mc_deltafactor, true);

    MnvH1D* data_final = GetMnvH1D(f_data_final, var_name);
    MnvH1D* mc_final = GetMnvH1D(f_mc_final, var_name);
  
    TH1D* h_data_final = GetBinNormalizedTH1D(data_final, true);
    TH1D* h_mc_final = GetBinNormalizedTH1D(mc_final, true);

    h_data_deltafactor->SetMarkerStyle(22);
    h_data_deltafactor->SetMarkerSize(1.5);
    h_data_deltafactor->SetMarkerColor(kMagenta);
    h_data_deltafactor->SetLineWidth(2);
    h_data_deltafactor->SetLineColor(kMagenta);

    h_mc_deltafactor->SetLineWidth(3);
    h_mc_deltafactor->SetLineColor(kBlue);
    h_mc_deltafactor->SetFillColor(kWhite);

    h_data_final->SetMarkerStyle(20);
    h_data_final->SetLineStyle(1);
    h_data_final->SetMarkerSize(1);
    h_data_final->SetMarkerColor(kBlack);
    h_data_final->SetLineWidth(2);
    h_data_final->SetLineColor(kBlack);

    h_mc_final->SetLineWidth(3);
    h_mc_final->SetLineStyle(1);
    h_mc_final->SetLineColor(kRed);
    h_mc_final->SetFillColor(kWhite);
    
    TCanvas* c = new TCanvas("c","c",800,800);

    if (thesisStyle){
        gStyle->SetEndErrorSize(6);

        h_data_deltafactor->GetXaxis()->SetTitleFont(62);
        h_data_deltafactor->GetXaxis()->SetTitleSize(0.06);
        h_data_deltafactor->GetXaxis()->CenterTitle();
        h_data_deltafactor->GetXaxis()->SetTitleOffset(1.15);
        h_data_deltafactor->GetXaxis()->SetLabelFont(42);
        h_data_deltafactor->GetXaxis()->SetLabelSize(0.05);
        h_data_deltafactor->GetXaxis()->SetNdivisions(408);

        h_data_deltafactor->GetYaxis()->SetTitleFont(62);
        h_data_deltafactor->GetYaxis()->SetTitleSize(0.06);
        //h_data_deltafactor->GetYaxis()->CenterTitle();
        h_data_deltafactor->GetYaxis()->SetTitleOffset(1.1);
        h_data_deltafactor->GetYaxis()->SetLabelFont(42);
        h_data_deltafactor->GetYaxis()->SetLabelSize(0.05);
        TGaxis::SetMaxDigits(3);
    }

    h_data_deltafactor->SetMaximum(h_data_deltafactor->GetMaximum()*1.75);
    h_data_deltafactor->SetMinimum(0.0);
    h_data_deltafactor->Draw("E1 X0");
    h_mc_deltafactor->Draw("HIST SAME");
    h_data_final->Draw("E1 X0 SAME");
    h_mc_final->Draw("HIST SAME");

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");
    text.DrawLatex(0.30,0.83,"Data POT: 3.33E+20");

    // TLegend
    TLegend *legend = new TLegend(0.43,0.70,0.9,0.9);  
    ApplyStyle_Legend(legend);
    legend->SetTextSize(0.03);
    legend->AddEntry(h_data_final, "Central Value Data", "lep");
    legend->AddEntry(h_mc_final, "Central Value GENIE", "l");
    legend->AddEntry(h_data_deltafactor, "CC-RES Supp. Data", "lep");
    legend->AddEntry(h_mc_deltafactor, "CC-RES Supp. GENIE", "l");
    legend->Draw();

    // Save Plot 
    gStyle->SetOptStat(0); 
    c->Update();
    std::string plot_name = var_name + "_Comp_DeltaFactor.png";
    c->Print(Form("%s%s",plotDir.c_str(),plot_name.c_str()), "png");

    delete h_data_deltafactor;
    delete data_deltafactor;
    delete h_data_final;
    delete h_mc_deltafactor;
    delete mc_deltafactor;
    delete h_mc_final;
    delete c;
    delete f_data_deltafactor;
    delete f_data_final;
    delete f_mc_deltafactor;
    delete f_mc_final;
}

void CCProtonPi0_Plotter::Draw_Comparison_Nominal(std::string var_name)
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    //std::string final_data_dir = "/minerva/data/users/oaltinok/NTupleAnalysis_Best/Data/Analyzed/CrossSection.root";
    //std::string final_mc_dir = "/minerva/data/users/oaltinok/NTupleAnalysis_Best/MC/Analyzed/CrossSection.root";

    std::string final_data_dir = "/minerva/data/users/oaltinok/NTupleAnalysis/Data/Analyzed/CrossSection.root";
    std::string final_mc_dir = "/minerva/data/users/oaltinok/NTupleAnalysis/MC/Analyzed/CrossSection.root";

    std::string nominal_data_dir = "/minerva/data/users/oaltinok/NTupleAnalysis_Nominal/Data/Analyzed/CrossSection.root";
    std::string nominal_mc_dir = "/minerva/data/users/oaltinok/NTupleAnalysis_Nominal/MC/Analyzed/CrossSection.root";

    TFile* f_data_final = new TFile(final_data_dir.c_str());
    TFile* f_mc_final = new TFile(final_mc_dir.c_str());

    TFile* f_data_nominal = new TFile(nominal_data_dir.c_str());
    TFile* f_mc_nominal = new TFile(nominal_mc_dir.c_str());
    
    // --------------------------------------------------------------------
    // Get Histograms 
    // --------------------------------------------------------------------
    MnvH1D* data_nominal = GetMnvH1D(f_data_nominal, var_name);
    MnvH1D* mc_nominal = GetMnvH1D(f_mc_nominal, var_name);
  
    TH1D* h_data_nominal = GetBinNormalizedTH1D(data_nominal, true);
    TH1D* h_mc_nominal = GetBinNormalizedTH1D(mc_nominal, true);

    MnvH1D* data_final = GetMnvH1D(f_data_final, var_name);
    MnvH1D* mc_final = GetMnvH1D(f_mc_final, var_name);
  
    TH1D* h_data_final = GetBinNormalizedTH1D(data_final, true);
    TH1D* h_mc_final = GetBinNormalizedTH1D(mc_final, true);

    h_data_nominal->SetMarkerStyle(20);
    h_data_nominal->SetMarkerSize(1);
    h_data_nominal->SetMarkerColor(kBlue);
    h_data_nominal->SetLineWidth(2);
    h_data_nominal->SetLineColor(kBlue);

    h_mc_nominal->SetLineWidth(3);
    h_mc_nominal->SetLineColor(kBlue);
    h_mc_nominal->SetFillColor(kWhite);

    h_data_final->SetMarkerStyle(20);
    h_data_final->SetLineStyle(1);
    h_data_final->SetMarkerSize(1);
    h_data_final->SetMarkerColor(kBlack);
    h_data_final->SetLineWidth(2);
    h_data_final->SetLineColor(kBlack);

    h_mc_final->SetLineWidth(3);
    h_mc_final->SetLineStyle(1);
    h_mc_final->SetLineColor(kRed);
    h_mc_final->SetFillColor(kWhite);
    
    TCanvas* c = new TCanvas("c","c",800,800);

    if (thesisStyle){
        gStyle->SetEndErrorSize(6);

        h_data_nominal->GetXaxis()->SetTitleFont(62);
        h_data_nominal->GetXaxis()->SetTitleSize(0.06);
        h_data_nominal->GetXaxis()->CenterTitle();
        h_data_nominal->GetXaxis()->SetTitleOffset(1.15);
        h_data_nominal->GetXaxis()->SetLabelFont(42);
        h_data_nominal->GetXaxis()->SetLabelSize(0.05);
        h_data_nominal->GetXaxis()->SetNdivisions(408);

        h_data_nominal->GetYaxis()->SetTitleFont(62);
        h_data_nominal->GetYaxis()->SetTitleSize(0.06);
        //h_data_nominal->GetYaxis()->CenterTitle();
        h_data_nominal->GetYaxis()->SetTitleOffset(1.1);
        h_data_nominal->GetYaxis()->SetLabelFont(42);
        h_data_nominal->GetYaxis()->SetLabelSize(0.05);
        TGaxis::SetMaxDigits(3);
    }

    h_data_nominal->SetMaximum(h_data_nominal->GetMaximum()*1.75);
    h_data_nominal->SetMinimum(0.0);
    h_data_nominal->Draw("E1 X0");
    h_mc_nominal->Draw("HIST SAME");
    h_data_final->Draw("E1 X0 SAME");
    h_mc_final->Draw("HIST SAME");

    // Add Normalization Labels
    TLatex text;
    text.SetNDC();
    text.SetTextColor(kBlue);
    text.SetTextSize(0.03);
    text.SetTextAlign(22);
    text.DrawLatex(0.30,0.87,"POT Normalized");
    text.DrawLatex(0.30,0.83,"Data POT: 3.33E+20");

    // TLegend
    TLegend *legend = new TLegend(0.45,0.70,0.9,0.9);  
    ApplyStyle_Legend(legend);
    legend->AddEntry(h_data_final, "Tuned Data", "lep");
    legend->AddEntry(h_mc_final, "Tuned GENIE", "l");
    legend->AddEntry(h_data_nominal, "Nominal Data", "lep");
    legend->AddEntry(h_mc_nominal, "Nominal GENIE", "l");
    legend->Draw();

    // Add Text
//    int nBins = h_data_nominal->GetNbinsX();
//    TLatex text;
//    text.SetNDC();
//    text.SetTextSize(0.03);
//    text.DrawLatex(0.17,0.85,Form("%s%3.2f", "GENIE CV #chi^{2} = ", Calc_ChiSq(h_data_nominal, h_mc_nominal,1,nBins)));
//    text.DrawLatex(0.17,0.81,Form("%s%3.2f", "#Delta Supp. #chi^{2} = ", Calc_ChiSq(h_data_final, h_mc_final,1,nBins)));

    // Save Plot 
    gStyle->SetOptStat(0); 
    c->Update();
    std::string plot_name = var_name + "_Comp_Nominal.png";
    c->Print(Form("%s%s",plotDir.c_str(),plot_name.c_str()), "png");

    delete h_data_nominal;
    delete data_nominal;
    delete h_data_final;
    delete h_mc_nominal;
    delete mc_nominal;
    delete h_mc_final;
    delete c;
    delete f_data_nominal;
    delete f_data_final;
    delete f_mc_nominal;
    delete f_mc_final;
}

// Returns "new" TH1D's
void CCProtonPi0_Plotter::GENIE_Tuning_GetHistograms(std::string var_name, std::string data_var, std::string mc_var, TH1D* &data_nominal, TH1D* &data_tuned_test, TH1D* &data_tuned_best, TH1D* &mc_nominal, TH1D* &mc_tuned_test, TH1D* &mc_tuned_best)
{
    std::string data_dir_nominal = "/minerva/data/users/oaltinok/NTupleAnalysis_Nominal/Data/Analyzed/CrossSection.root";
    std::string data_dir_tuned_test = "/minerva/data/users/oaltinok/NTupleAnalysis_DeltaRES/Data/Analyzed/CrossSection.root";
    //std::string data_dir_tuned_best = "/minerva/data/users/oaltinok/NTupleAnalysis_Only2p2h/Data/Analyzed/CrossSection.root";
    std::string data_dir_tuned_best = "/minerva/data/users/oaltinok/NTupleAnalysis_Best/Data/Analyzed/CrossSection.root";

    std::string mc_dir_nominal = "/minerva/data/users/oaltinok/NTupleAnalysis_Nominal/MC/Analyzed/CrossSection.root";
    std::string mc_dir_tuned_test = "/minerva/data/users/oaltinok/NTupleAnalysis_DeltaRES/MC/Analyzed/CrossSection.root";
    //std::string mc_dir_tuned_best = "/minerva/data/users/oaltinok/NTupleAnalysis_Only2p2h/MC/Analyzed/CrossSection.root";
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
    TCanvas* c = new TCanvas("c","c",800,800);

    // Style Histograms
    ratio_nominal->GetYaxis()->SetTitle("Ratio");
    ratio_nominal->SetMinimum(0.0);
    ratio_nominal->SetMaximum(2);
    ratio_nominal->SetLineColor(kBlue);
    ratio_nominal->SetLineWidth(3);
    ratio_nominal->SetFillColor(kWhite);

    ratio_test->GetYaxis()->SetTitle("Ratio");
    ratio_test->SetLineColor(kRed);
    ratio_test->SetLineWidth(3);
    ratio_test->SetFillColor(kWhite);

    ratio_best->GetYaxis()->SetTitle("Ratio");
    ratio_best->SetLineColor(kBlack);
    ratio_best->SetLineWidth(3);
    ratio_best->SetFillColor(kWhite);

    if (thesisStyle){
        gStyle->SetEndErrorSize(6);

        ratio_nominal->GetXaxis()->SetTitleFont(62);
        ratio_nominal->GetXaxis()->SetTitleSize(0.06);
        ratio_nominal->GetXaxis()->CenterTitle();
        ratio_nominal->GetXaxis()->SetTitleOffset(1.15);
        ratio_nominal->GetXaxis()->SetLabelFont(42);
        ratio_nominal->GetXaxis()->SetLabelSize(0.05);
        ratio_nominal->GetXaxis()->SetNdivisions(408);

        ratio_nominal->GetYaxis()->SetTitleFont(62);
        ratio_nominal->GetYaxis()->SetTitleSize(0.06);
        //ratio_nominal->GetYaxis()->CenterTitle();
        ratio_nominal->GetYaxis()->SetTitleOffset(1.1);
        ratio_nominal->GetYaxis()->SetLabelFont(42);
        ratio_nominal->GetYaxis()->SetLabelSize(0.05);
        TGaxis::SetMaxDigits(3);
    }

    // Plot
    ratio_nominal->Draw("HIST");
    ratio_best->Draw("HIST SAME");
    //ratio_test->Draw("HIST SAME");
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

    TLegend *legend = new TLegend(0.35,0.8,0.9,0.9);  
    ApplyStyle_Legend(legend);
    legend->AddEntry(ratio_best, "Tuned (Data/GENIE)", "l" );
    legend->AddEntry(ratio_nominal, "Nominal (Data/GENIE)", "l" );
    //legend->AddEntry(ratio_test, "#Delta^{++} Decay(Data/GENIE)", "l" );
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
    TCanvas* c = new TCanvas("c","c",800,800);

    // Style Histograms
    mc_ratio_best->GetYaxis()->SetTitle("Ratio");
    mc_ratio_best->SetMaximum(2.0);
    mc_ratio_best->SetMinimum(0.0);
    mc_ratio_best->SetLineColor(kRed);
    mc_ratio_best->SetLineWidth(3);
    mc_ratio_best->SetFillColor(kWhite);

    data_ratio_best->SetLineColor(kBlack);
    data_ratio_best->SetLineWidth(3);
    data_ratio_best->SetFillColor(kWhite);

    mc_ratio_test->SetLineColor(kRed);
    mc_ratio_test->SetLineWidth(3);
    mc_ratio_test->SetLineStyle(2);
    mc_ratio_test->SetFillColor(kWhite);

    data_ratio_test->SetLineColor(kBlack);
    data_ratio_test->SetLineWidth(3);
    data_ratio_test->SetLineStyle(2);
    data_ratio_test->SetFillColor(kWhite);

    if (thesisStyle){
        gStyle->SetEndErrorSize(6);

        mc_ratio_best->GetXaxis()->SetTitleFont(62);
        mc_ratio_best->GetXaxis()->SetTitleSize(0.06);
        mc_ratio_best->GetXaxis()->CenterTitle();
        mc_ratio_best->GetXaxis()->SetTitleOffset(1.15);
        mc_ratio_best->GetXaxis()->SetLabelFont(42);
        mc_ratio_best->GetXaxis()->SetLabelSize(0.05);
        mc_ratio_best->GetXaxis()->SetNdivisions(408);

        mc_ratio_best->GetYaxis()->SetTitleFont(62);
        mc_ratio_best->GetYaxis()->SetTitleSize(0.06);
        //mc_ratio_best->GetYaxis()->CenterTitle();
        mc_ratio_best->GetYaxis()->SetTitleOffset(1.1);
        mc_ratio_best->GetYaxis()->SetLabelFont(42);
        mc_ratio_best->GetYaxis()->SetLabelSize(0.05);
        TGaxis::SetMaxDigits(3);
    }

    // Plot
    mc_ratio_best->Draw("HIST");
    data_ratio_best->Draw("HIST SAME");
    //mc_ratio_test->Draw("HIST SAME");
    //data_ratio_test->Draw("HIST SAME");
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

    TLegend *legend = new TLegend(0.30,0.8,0.9,0.9);  
    ApplyStyle_Legend(legend);
    legend->AddEntry(data_ratio_best, "Data(Tuned/Nominal)", "l");
    legend->AddEntry(mc_ratio_best, "GENIE(Tuned/Nominal)", "l");
    //legend->AddEntry(data_ratio_test, "Data(#Delta{++} Decay) / Data", "l");
    //legend->AddEntry(mc_ratio_test, "GENIE(#Delta{++} Decay) / GENIE", "l");
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
    //GENIE_Tuning_Ratio("muon_theta", "all", "mc_reco_all");
//    GENIE_Tuning_Ratio("pi0_P", "all", "mc_reco_all");
//    GENIE_Tuning_Ratio("pi0_KE", "all", "mc_reco_all");
//    GENIE_Tuning_Ratio("pi0_theta", "all", "mc_reco_all");
    //GENIE_Tuning_Ratio("QSq", "all", "mc_reco_all");
//    GENIE_Tuning_Ratio("W", "all", "mc_reco_all");
//    GENIE_Tuning_Ratio("Enu", "all", "mc_reco_all");
//
//    GENIE_Tuning_Ratio("muon_P", "bckg_subtracted", "bckg_subtracted");
    //GENIE_Tuning_Ratio("muon_theta", "bckg_subtracted", "bckg_subtracted");
//    GENIE_Tuning_Ratio("pi0_P", "bckg_subtracted", "bckg_subtracted");
//    GENIE_Tuning_Ratio("pi0_KE", "bckg_subtracted", "bckg_subtracted");
//    GENIE_Tuning_Ratio("pi0_theta", "bckg_subtracted", "bckg_subtracted");
    //GENIE_Tuning_Ratio("QSq", "bckg_subtracted", "bckg_subtracted");
//    GENIE_Tuning_Ratio("W", "bckg_subtracted", "bckg_subtracted");
//    GENIE_Tuning_Ratio("Enu", "bckg_subtracted", "bckg_subtracted");

//    GENIE_Tuning_Ratio("muon_P", "bckg_estimated", "bckg_estimated");
    //GENIE_Tuning_Ratio("muon_theta", "bckg_estimated", "bckg_estimated");
//    GENIE_Tuning_Ratio("pi0_P", "bckg_estimated", "bckg_estimated");
//    GENIE_Tuning_Ratio("pi0_KE", "bckg_estimated", "bckg_estimated");
//    GENIE_Tuning_Ratio("pi0_theta", "bckg_estimated", "bckg_estimated");
    //GENIE_Tuning_Ratio("QSq", "bckg_estimated", "bckg_estimated");
//    GENIE_Tuning_Ratio("W", "bckg_estimated", "bckg_estimated");
//    GENIE_Tuning_Ratio("Enu", "bckg_estimated", "bckg_estimated");

//    GENIE_Tuning_Ratio("muon_P", "unfolded", "unfolded");
    //GENIE_Tuning_Ratio("muon_theta", "unfolded", "unfolded");
//    GENIE_Tuning_Ratio("pi0_P", "unfolded", "unfolded");
//    GENIE_Tuning_Ratio("pi0_KE", "unfolded", "unfolded");
//    GENIE_Tuning_Ratio("pi0_theta", "unfolded", "unfolded");
    //GENIE_Tuning_Ratio("QSq", "unfolded", "unfolded");
//    GENIE_Tuning_Ratio("W", "unfolded", "unfolded");
//    GENIE_Tuning_Ratio("Enu", "unfolded", "unfolded");
//
//    GENIE_Tuning_Ratio("muon_P", "efficiency_corrected", "efficiency_corrected");
    //GENIE_Tuning_Ratio("muon_theta", "efficiency_corrected", "efficiency_corrected");
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

