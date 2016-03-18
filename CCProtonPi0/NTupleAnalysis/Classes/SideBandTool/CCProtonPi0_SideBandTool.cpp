#include "CCProtonPi0_SideBandTool.h"

using namespace PlotUtils;

CCProtonPi0_SideBandTool::CCProtonPi0_SideBandTool() : CCProtonPi0_NTupleAnalysis()
{
    OpenRootFiles();
    initSideBands();
}

void CCProtonPi0_SideBandTool::OpenRootFiles()
{
    std::string rootDir;

    rootDir = Folder_List::rootDir_sideBand_Original_mc;
    Original.f_mc = new TFile(rootDir.c_str());

    rootDir = Folder_List::rootDir_sideBand_Original_data;
    Original.f_data = new TFile(rootDir.c_str());

    rootDir = Folder_List::rootDir_sideBand_Michel_mc;
    Michel.f_mc = new TFile(rootDir.c_str());

    rootDir = Folder_List::rootDir_sideBand_Michel_data;
    Michel.f_data = new TFile(rootDir.c_str());

    rootDir = Folder_List::rootDir_sideBand_pID_mc;
    pID.f_mc = new TFile(rootDir.c_str());

    rootDir = Folder_List::rootDir_sideBand_pID_data;
    pID.f_data = new TFile(rootDir.c_str());

    //rootDir = Folder_List::rootDir_sideBand_LowInvMass_mc;
    rootDir = Folder_List::rootDir_sideBand_HighInvMass_mc;
    LowInvMass.f_mc = new TFile(rootDir.c_str());

    //rootDir = Folder_List::rootDir_sideBand_LowInvMass_data;
    rootDir = Folder_List::rootDir_sideBand_HighInvMass_data;
    LowInvMass.f_data = new TFile(rootDir.c_str());

    rootDir = Folder_List::rootDir_sideBand_HighInvMass_mc;
    HighInvMass.f_mc = new TFile(rootDir.c_str());

    rootDir = Folder_List::rootDir_sideBand_HighInvMass_data;
    HighInvMass.f_data = new TFile(rootDir.c_str());
}

void CCProtonPi0_SideBandTool::initSideBands()
{
    SetNames(Original, "Original");
    GetTH1D(Original.f_data, Original.data, "hCut_pi0invMass_0");   
    GetTH1D(Original.f_mc, Original.mc_total, "hCut_pi0invMass_0");   
    GetTH1D(Original.f_mc, Original.signal[0], "hCut_pi0invMass_1");   
    GetTH1D(Original.f_mc, Original.WithPi0[0], "hCut_pi0invMass_3");   
    GetTH1D(Original.f_mc, Original.QELike[0], "hCut_pi0invMass_4");   
    GetTH1D(Original.f_mc, Original.SinglePiPlus[0], "hCut_pi0invMass_5");   
    GetTH1D(Original.f_mc, Original.Other[0], "hCut_pi0invMass_6");   

    SetNames(Michel, "Michel");
    GetTH1D(Michel.f_data, Michel.data, "hCut_pi0invMass_0");   
    GetTH1D(Michel.f_mc, Michel.mc_total, "hCut_pi0invMass_0");   
    GetTH1D(Michel.f_mc, Michel.signal[0], "hCut_pi0invMass_1");   
    GetTH1D(Michel.f_mc, Michel.WithPi0[0], "hCut_pi0invMass_3");   
    GetTH1D(Michel.f_mc, Michel.QELike[0], "hCut_pi0invMass_4");   
    GetTH1D(Michel.f_mc, Michel.SinglePiPlus[0], "hCut_pi0invMass_5");   
    GetTH1D(Michel.f_mc, Michel.Other[0], "hCut_pi0invMass_6");   

    SetNames(pID, "pID");
    GetTH1D(pID.f_data, pID.data, "hCut_pi0invMass_0");   
    GetTH1D(pID.f_mc, pID.mc_total, "hCut_pi0invMass_0");   
    GetTH1D(pID.f_mc, pID.signal[0], "hCut_pi0invMass_1");   
    GetTH1D(pID.f_mc, pID.WithPi0[0], "hCut_pi0invMass_3");   
    GetTH1D(pID.f_mc, pID.QELike[0], "hCut_pi0invMass_4");   
    GetTH1D(pID.f_mc, pID.SinglePiPlus[0], "hCut_pi0invMass_5");   
    GetTH1D(pID.f_mc, pID.Other[0], "hCut_pi0invMass_6");   

    SetNames(LowInvMass, "LowInvMass");
    GetTH1D(LowInvMass.f_data, LowInvMass.data, "hCut_pi0invMass_0");   
    GetTH1D(LowInvMass.f_mc, LowInvMass.mc_total, "hCut_pi0invMass_0");   
    GetTH1D(LowInvMass.f_mc, LowInvMass.signal[0], "hCut_pi0invMass_1");   
    GetTH1D(LowInvMass.f_mc, LowInvMass.WithPi0[0], "hCut_pi0invMass_3");   
    GetTH1D(LowInvMass.f_mc, LowInvMass.QELike[0], "hCut_pi0invMass_4");   
    GetTH1D(LowInvMass.f_mc, LowInvMass.SinglePiPlus[0], "hCut_pi0invMass_5");   
    GetTH1D(LowInvMass.f_mc, LowInvMass.Other[0], "hCut_pi0invMass_6");   

    SetNames(HighInvMass, "HighInvMass");
    GetTH1D(HighInvMass.f_data, HighInvMass.data, "hCut_pi0invMass_0");   
    GetTH1D(HighInvMass.f_mc, HighInvMass.mc_total, "hCut_pi0invMass_0");   
    GetTH1D(HighInvMass.f_mc, HighInvMass.signal[0], "hCut_pi0invMass_1");   
    GetTH1D(HighInvMass.f_mc, HighInvMass.WithPi0[0], "hCut_pi0invMass_3");   
    GetTH1D(HighInvMass.f_mc, HighInvMass.QELike[0], "hCut_pi0invMass_4");   
    GetTH1D(HighInvMass.f_mc, HighInvMass.SinglePiPlus[0], "hCut_pi0invMass_5");   
    GetTH1D(HighInvMass.f_mc, HighInvMass.Other[0], "hCut_pi0invMass_6");   
}

void CCProtonPi0_SideBandTool::GetTH1D(TFile* f, TH1D* &h, std::string var_name)
{
    MnvH1D* temp;
    temp = new MnvH1D(*(MnvH1D*)f->Get(var_name.c_str()));
    h = new TH1D(temp->GetCVHistoWithStatError()); 
    delete temp;
}

void CCProtonPi0_SideBandTool::SetNames(SideBand &sb, std::string name)
{
    sb.name = name;
    sb.model_names[0] = "Signal";
    sb.model_names[1] = "WithPi0";
    sb.model_names[2] = "QELike";
    sb.model_names[3] = "SinglePiPlus";
    sb.model_names[4] = "Other";
    sb.model_names[5] = "Total";
}

CCProtonPi0_SideBandTool::~CCProtonPi0_SideBandTool()
{

}

void CCProtonPi0_SideBandTool::Plot()
{
    Plot(Original);
    Plot(Michel);
    Plot(pID);
    Plot(LowInvMass);
    Plot(HighInvMass);
}

void CCProtonPi0_SideBandTool::Plot(SideBand &sb)
{
    std::cout<<"Plotting "<<sb.name<<std::endl;
    // Original 
    Plot(sb, 0, true);
    Plot(sb, 0, false);
    
    // Modified
    Plot(sb, 1, true);
    Plot(sb, 1, false);
}

void CCProtonPi0_SideBandTool::Plot(SideBand &sb, int ind, bool isArea)
{
    std::string type;
    if (ind == 0) type = "Original";
    else type = "Modified";

    std::string norm;
    if (isArea) norm = "Area";
    else norm = "POT";
    std::string plot_title = "Side Band: " + sb.name + " " + type + " " + norm + " Normalized";

    TCanvas* c = new TCanvas("c","c",1280,800);
    THStack* hs = new THStack("hs",plot_title.c_str());

    // Get & Scale MC Models
    ColorHists(sb);
    double mc_ratio = GetMCScaleRatio(sb, isArea);
    TH1D* h_signal = new TH1D (*(sb.signal[ind]));
    TH1D* h_WithPi0 = new TH1D (*(sb.WithPi0[ind]));
    TH1D* h_QELike = new TH1D (*(sb.QELike[ind]));
    TH1D* h_SinglePiPlus = new TH1D (*(sb.SinglePiPlus[ind]));
    TH1D* h_Other = new TH1D (*(sb.Other[ind]));

    h_signal->Scale(mc_ratio);
    h_WithPi0->Scale(mc_ratio);
    h_QELike->Scale(mc_ratio);
    h_SinglePiPlus->Scale(mc_ratio);
    h_Other->Scale(mc_ratio);

    // Get MC Total Area
    double area_mc = 0.0;
    area_mc += h_signal->Integral();
    area_mc += h_WithPi0->Integral();
    area_mc += h_QELike->Integral();
    area_mc += h_SinglePiPlus->Integral();
    area_mc += h_Other->Integral();

    // Plot MC Models
    hs->Add(h_WithPi0);  
    hs->Add(h_QELike);  
    hs->Add(h_SinglePiPlus);  
    hs->Add(h_Other);  
    hs->Add(h_signal);  

    hs->Draw("HIST");
    hs->GetXaxis()->SetTitle("#pi^{0} Invariant Mass [MeV]");
    hs->GetYaxis()->SetTitle("N(Events)");

    // Plot Data
    sb.data->Draw("SAME E1 X0");
    double area_data = sb.data->Integral();

    // Add Legend
    TLegend *legend = new TLegend(0.7,0.68,0.9,0.9);  
    legend->AddEntry(sb.data, "Data");
    legend->AddEntry(h_signal, "Signal", "f");
    legend->AddEntry(h_Other, "Bckg: Other", "f");
    legend->AddEntry(h_SinglePiPlus, "Bckg: 1 #pi^{+}", "f");
    legend->AddEntry(h_QELike, "Bckg: QE Like", "f");
    legend->AddEntry(h_WithPi0, "Bckg: #pi^{0} + X", "f");
    legend->SetTextSize(0.03);
    legend->Draw();

    // Add Pi0 InvMass Lines
    double hist_max = sb.data->GetMaximum();
    hs->SetMaximum(hist_max * 1.2);
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

    // Add Weights as Text to Modified Plot 
    if (ind != 0){
        TLatex* text = new TLatex;
        text->SetTextSize(0.03);
        text->SetNDC();
        text->DrawLatex(0.7, 0.64, "Fit Results");
        text->DrawLatex(0.7, 0.60, Form("%s%3.2f", "Fit #chi^{2} = ", ChiSq));
        text->DrawLatex(0.7, 0.56, Form("#color[4]{%s%3.2f}", "wgt(SinglePiPlus) = ", wgt_SinglePiPlus));
        text->DrawLatex(0.7, 0.52, Form("%s%3.2f", "wgt(QELike) = ", wgt_QELike));
        text->DrawLatex(0.7, 0.48, Form("#color[2]{%s%3.2f}", "wgt(WithPi0) = ", wgt_WithPi0));
        delete text;
    }

    // Add Plot-Area and Plot-ChiSq
    TLatex* areaText = new TLatex;
    areaText->SetNDC();
    areaText->SetTextSize(0.03);
    areaText->SetTextColor(kBlue);
    areaText->DrawLatex(0.15, 0.85,Form("Area(Data)/Area(MC) = %3.2f",area_data/area_mc));
    double plot_chisq = calc_ChiSq(sb, ind, isArea);
    areaText->DrawLatex(0.15, 0.81, Form("Plot #chi^{2} = %3.2f", plot_chisq));
    delete areaText;

    // Plot Output
    std::string plotDir = Folder_List::plotDir_SideBand;
    std::string out_name;
    out_name = plotDir + sb.name + "_" + type + "_" + norm + ".png"; 

    c->Print(out_name.c_str(),"png");

    delete legend;
    delete hs;
    delete c;

    // Create Another Canvas for Data MC Difference Plot
    TCanvas* c2 = new TCanvas("c2","c2",1280,800);
    
    TH1D* h_data_mc_diff = new TH1D( "h_data_mc_diff","Data - MC for each bin",100,-100,100);
    h_data_mc_diff->GetXaxis()->SetTitle("N(Data)-N(MC)");
    h_data_mc_diff->GetYaxis()->SetTitle("N(Events)");
 
    for (int i = 1; i <= sb.data->GetNbinsX(); ++i){
        double diff = sb.data->GetBinContent(i) - sb.mc_total->GetBinContent(i)*mc_ratio;
        h_data_mc_diff->Fill(diff);
    }
    
    h_data_mc_diff->SetLineColor(kRed);
    h_data_mc_diff->SetLineWidth(3);
    h_data_mc_diff->SetFillColor(kRed);
    h_data_mc_diff->SetFillStyle(3010);

    h_data_mc_diff->Draw("HIST");
    gPad->Update();
    gStyle->SetOptStat(111111); 

    // Plot Output
    out_name = plotDir + sb.name + "_" + type + "_" + norm + "Difference.png"; 

    c2->Print(out_name.c_str(),"png");

    delete c2;
}

void CCProtonPi0_SideBandTool::ColorHists(SideBand &sb)
{
    // MC Models
    for (int i = 0; i < 2; ++i){
        sb.signal[i]->SetFillColor(kGreen);
        sb.signal[i]->SetFillStyle(3001);

        sb.WithPi0[i]->SetFillColor(kRed);
        sb.WithPi0[i]->SetFillStyle(3001);

        sb.QELike[i]->SetFillColor(kOrange);
        sb.QELike[i]->SetFillStyle(3001);

        sb.SinglePiPlus[i]->SetFillColor(kBlue);
        sb.SinglePiPlus[i]->SetFillStyle(3001);

        sb.Other[i]->SetFillColor(kGray);
        sb.Other[i]->SetFillStyle(3001);
    }

    // Data
    sb.data->SetMarkerColor(kBlack);
    sb.data->SetMarkerStyle(20);
    sb.data->SetMarkerSize(1);
    sb.data->SetLineWidth(1);
    sb.data->SetLineColor(kBlack);
    sb.data->SetFillStyle(0);
}

double CCProtonPi0_SideBandTool::GetMCScaleRatio(SideBand &sb, bool isArea)
{
    double mc_ratio;
    if (isArea){
        double data_area = sb.data->Integral();
        double mc_area = sb.mc_total->Integral();
        mc_ratio = data_area / mc_area;
    }else{
        mc_ratio = POT_ratio;
    }

    return mc_ratio;
}

void CCProtonPi0_SideBandTool::ApplyFitResults(double chisq, double w_WithPi0, double w_QELike, double w_SinglePiPlus )
{
    ChiSq = chisq;
    wgt_WithPi0 = w_WithPi0;
    wgt_QELike = w_QELike;
    wgt_SinglePiPlus = w_SinglePiPlus;

    ApplyFitResults();
}

void CCProtonPi0_SideBandTool::ApplyFitResults()
{
    ApplyFitResults(Original);
    ApplyFitResults(Michel);
    ApplyFitResults(pID);
    ApplyFitResults(LowInvMass);
    ApplyFitResults(HighInvMass);
}

void CCProtonPi0_SideBandTool::ApplyFitResults(SideBand &sb)
{
    std::cout<<"Applying Fit Result to "<<sb.name<<std::endl;
    // Clone Original Histograms
    sb.signal[1] = new TH1D (*sb.signal[0]);
    sb.WithPi0[1] = new TH1D (*sb.WithPi0[0]);
    sb.QELike[1] = new TH1D (*sb.QELike[0]);
    sb.SinglePiPlus[1] = new TH1D (*sb.SinglePiPlus[0]);
    sb.Other[1] = new TH1D (*sb.Other[0]);

    // Scale 
    sb.WithPi0[1]->Scale(wgt_WithPi0);
    sb.QELike[1]->Scale(wgt_QELike);
    sb.SinglePiPlus[1]->Scale(wgt_SinglePiPlus);
}


double CCProtonPi0_SideBandTool::calc_ChiSq(SideBand &sb, int ind, bool isArea)
{
    double mc_ratio = GetMCScaleRatio(sb, isArea);
    double chi_sq = 0.0;

    for (int i = 1; i <= sb.data->GetNbinsX(); ++i){
        // Get N(Events) in Single Bin
        double data = sb.data->GetBinContent(i);
        double signal = sb.signal[ind]->GetBinContent(i);
        double WithPi0 = sb.WithPi0[ind]->GetBinContent(i);
        double QELike = sb.QELike[ind]->GetBinContent(i);
        double SinglePiPlus = sb.SinglePiPlus[ind]->GetBinContent(i);
        double Other = sb.Other[ind]->GetBinContent(i);

        // Add All MC and scale them
        double MC_total = (signal + WithPi0 + QELike + SinglePiPlus + Other) * mc_ratio;
        
        double delta = pow((data-MC_total),2) / data;
        chi_sq += delta;
    }
    
    return chi_sq;
}

