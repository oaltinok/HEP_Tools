#include "CCProtonPi0_SideBandTool.h"

using namespace PlotUtils;

void CCProtonPi0_SideBandTool::Fit()
{
    Fit(Michel);
    Fit(pID);
    Fit(LowInvMass, true, 1, 3);
}

void CCProtonPi0_SideBandTool::Fit(SideBand &sb, bool isLimitedFit, int first_bin, int last_bin)
{
    std::cout<<"Fitting Side Band: "<<sb.name<<std::endl;

    // Crate TObjArray
    TObjArray* mc_models = new TObjArray(5);
    mc_models->Add(sb.signal[0]);
    mc_models->Add(sb.WithPi0[0]);
    mc_models->Add(sb.QELike[0]);
    mc_models->Add(sb.SinglePiPlus[0]);
    mc_models->Add(sb.Other[0]);

    // Fit data and mc_models
    TFractionFitter* fitter = new TFractionFitter(sb.data, mc_models, "Q");
    // Constrain Fractions [0,1]
    fitter->Constrain(0,0.0,1.0);
    fitter->Constrain(1,0.0,1.0);
    fitter->Constrain(2,0.0,1.0);
    fitter->Constrain(3,0.0,1.0);
    fitter->Constrain(4,0.0,1.0);

    if (isLimitedFit){
        fitter->SetRangeX(first_bin, last_bin);
    }

    Int_t status = fitter->Fit();

    if (status != 0) {
        std::cout<<"\tFit Error!"<<std::endl;
        return;
    }

    sb.fit = new TH1D (*(TH1D*)fitter->GetPlot());
    
    // Fit Results
    double total = 0;
    double total_err = 0;
    Double_t f;                         // Fraction
    Double_t f_err;                     // Fraction Error
    fitter->GetResult(0, f, f_err);     // access fitted value and uncertainty for parameter 1
    sb.fr_signal[1] = f;
    sb.fr_signal[2] = sb.fr_signal[1]/sb.fr_signal[0];
    sb.fr_signal[3] = f_err;
    total += f;
    total_err += f_err;

    fitter->GetResult(1, f, f_err);
    sb.fr_WithPi0[1] = f;
    sb.fr_WithPi0[2] = sb.fr_WithPi0[1]/sb.fr_WithPi0[0];
    sb.fr_WithPi0[3] = f_err; 
    total += f;
    total_err += f_err;
   
    fitter->GetResult(2, f, f_err); 
    sb.fr_QELike[1] = f;
    sb.fr_QELike[2] = sb.fr_QELike[1]/sb.fr_QELike[0];
    sb.fr_QELike[3] = f_err;  
    total += f;
    total_err += f_err;
    
    fitter->GetResult(3, f, f_err); 
    sb.fr_SinglePiPlus[1] = f;
    sb.fr_SinglePiPlus[2] = sb.fr_SinglePiPlus[1]/sb.fr_SinglePiPlus[0];
    sb.fr_SinglePiPlus[3] = f_err;  
    total += f;
    total_err += f_err;
    
    fitter->GetResult(4, f, f_err);
    sb.fr_Other[1] = f;
    sb.fr_Other[2] = sb.fr_Other[1]/sb.fr_Other[0];
    sb.fr_Other[3] = f_err;  
    total += f;
    total_err += f_err;
  
    sb.fr_Total[1] = total;
    sb.fr_Total[2] = 0;
    sb.fr_Total[3] = total_err;

    PrintFitResults(sb);
  
    ApplyFitResults(sb);
    
    // Plot Fit Results
    for (int i = 0; i < 2; ++i){
        Plot(sb, i, true);
        Plot(sb, i, false);
    }
  
    delete fitter;
    std::cout<<"Done!"<<std::endl;
}

CCProtonPi0_SideBandTool::CCProtonPi0_SideBandTool() : CCProtonPi0_NTupleAnalysis()
{
    std::cout<<"Initializing CCProtonPi0_SideBandTool"<<std::endl;

    double data_POT = 3.33009e+20; 
    double mc_POT = 2.73881e+21;
    POT_ratio = data_POT/mc_POT;

    OpenTextFile();
    OpenRootFiles();
    initSideBands();

    std::cout<<"Done!"<<std::endl;
    std::cout<<std::endl;
}

void CCProtonPi0_SideBandTool::OpenRootFiles()
{
    std::string rootDir;

    rootDir = Folder_List::rootDir_sideBand_Michel_mc;
    Michel.f_mc = new TFile(rootDir.c_str());

    rootDir = Folder_List::rootDir_sideBand_Michel_data;
    Michel.f_data = new TFile(rootDir.c_str());

    rootDir = Folder_List::rootDir_sideBand_pID_mc;
    pID.f_mc = new TFile(rootDir.c_str());

    rootDir = Folder_List::rootDir_sideBand_pID_data;
    pID.f_data = new TFile(rootDir.c_str());

    rootDir = Folder_List::rootDir_sideBand_LowInvMass_mc;
    LowInvMass.f_mc = new TFile(rootDir.c_str());

    rootDir = Folder_List::rootDir_sideBand_LowInvMass_data;
    LowInvMass.f_data = new TFile(rootDir.c_str());
}

void CCProtonPi0_SideBandTool::initSideBands()
{
    SetNames(Michel, "Michel");
    GetTH1D(Michel.f_data, Michel.data, "hCut_pi0invMass_0");   
    GetTH1D(Michel.f_mc, Michel.mc_total, "hCut_pi0invMass_0");   
    GetTH1D(Michel.f_mc, Michel.signal[0], "hCut_pi0invMass_1");   
    GetTH1D(Michel.f_mc, Michel.WithPi0[0], "hCut_pi0invMass_3");   
    GetTH1D(Michel.f_mc, Michel.QELike[0], "hCut_pi0invMass_4");   
    GetTH1D(Michel.f_mc, Michel.SinglePiPlus[0], "hCut_pi0invMass_5");   
    GetTH1D(Michel.f_mc, Michel.Other[0], "hCut_pi0invMass_6");   
    CalcMCRatios(Michel);

    SetNames(pID, "pID");
    GetTH1D(pID.f_data, pID.data, "hCut_pi0invMass_0");   
    GetTH1D(pID.f_mc, pID.mc_total, "hCut_pi0invMass_0");   
    GetTH1D(pID.f_mc, pID.signal[0], "hCut_pi0invMass_1");   
    GetTH1D(pID.f_mc, pID.WithPi0[0], "hCut_pi0invMass_3");   
    GetTH1D(pID.f_mc, pID.QELike[0], "hCut_pi0invMass_4");   
    GetTH1D(pID.f_mc, pID.SinglePiPlus[0], "hCut_pi0invMass_5");   
    GetTH1D(pID.f_mc, pID.Other[0], "hCut_pi0invMass_6");   
    CalcMCRatios(pID);

    SetNames(LowInvMass, "LowInvMass");
    GetTH1D(LowInvMass.f_data, LowInvMass.data, "hCut_pi0invMass_0");   
    GetTH1D(LowInvMass.f_mc, LowInvMass.mc_total, "hCut_pi0invMass_0");   
    GetTH1D(LowInvMass.f_mc, LowInvMass.signal[0], "hCut_pi0invMass_1");   
    GetTH1D(LowInvMass.f_mc, LowInvMass.WithPi0[0], "hCut_pi0invMass_3");   
    GetTH1D(LowInvMass.f_mc, LowInvMass.QELike[0], "hCut_pi0invMass_4");   
    GetTH1D(LowInvMass.f_mc, LowInvMass.SinglePiPlus[0], "hCut_pi0invMass_5");   
    GetTH1D(LowInvMass.f_mc, LowInvMass.Other[0], "hCut_pi0invMass_6");   
    CalcMCRatios(LowInvMass, true, 1, 3);
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

void CCProtonPi0_SideBandTool::CalcMCRatios(SideBand &sb, bool isLimited, int first_bin, int last_bin)
{
    // If the Integral is NOT Limited then integrate over ALL bins
    if (!isLimited){
        first_bin = 1;
        last_bin = sb.signal[0]->GetNbinsX();
    }

    double signal = sb.signal[0]->Integral(first_bin, last_bin);
    double WithPi0 = sb.WithPi0[0]->Integral(first_bin, last_bin);
    double QELike = sb.QELike[0]->Integral(first_bin, last_bin);
    double SinglePiPlus = sb.SinglePiPlus[0]->Integral(first_bin, last_bin);
    double Other = sb.Other[0]->Integral(first_bin, last_bin);
    
    double total = signal + WithPi0 + QELike + SinglePiPlus + Other;

    double r_signal = signal / total;
    double r_WithPi0 = WithPi0 / total;
    double r_QELike = QELike / total;
    double r_SinglePiPlus = SinglePiPlus / total;
    double r_Other = Other / total;

    sb.fr_signal[0] = r_signal;
    sb.fr_WithPi0[0] = r_WithPi0;
    sb.fr_QELike[0] = r_QELike;
    sb.fr_SinglePiPlus[0] = r_SinglePiPlus;
    sb.fr_Other[0] = r_Other;
    sb.fr_Total[0] = 1.0;
}

void CCProtonPi0_SideBandTool::PrintFitResults(SideBand &sb)
{
    using namespace std;
    textFile<<left;
    textFile<<setw(16)<<"Type"; 
    textFile<<setw(16)<<"MC Ratio"; 
    textFile<<setw(16)<<"Fit Ratio"; 
    textFile<<setw(16)<<"Fit/MC"; 
    textFile<<setw(16)<<"Error"; 
    textFile<<endl; 

    textFile<<setw(16)<<sb.model_names[0]; PrintRatio(sb.fr_signal);
    textFile<<setw(16)<<sb.model_names[1]; PrintRatio(sb.fr_WithPi0);
    textFile<<setw(16)<<sb.model_names[2]; PrintRatio(sb.fr_QELike);
    textFile<<setw(16)<<sb.model_names[3]; PrintRatio(sb.fr_SinglePiPlus);
    textFile<<setw(16)<<sb.model_names[4]; PrintRatio(sb.fr_Other);
    textFile<<setw(16)<<sb.model_names[5]; PrintRatio(sb.fr_Total);
    
    textFile<<endl;
}

void CCProtonPi0_SideBandTool::PrintRatio(double ratio[])
{
    using namespace std;
    textFile<<left;
    textFile<<setw(16)<<setprecision(3)<<ratio[0];
    textFile<<setw(16)<<setprecision(3)<<ratio[1];
    textFile<<setw(16)<<setprecision(3)<<ratio[2];
    textFile<<setw(16)<<setprecision(3)<<ratio[3];
    textFile<<endl;
}

void CCProtonPi0_SideBandTool::OpenTextFile()
{
    fileName = Folder_List::output + Folder_List::textOut + "SideBandTables.txt";
    CCProtonPi0_NTupleAnalysis::OpenTextFile(fileName, textFile);
}

CCProtonPi0_SideBandTool::~CCProtonPi0_SideBandTool()
{
    textFile.close();
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

    // Plot Fit
    sb.fit->Draw("SAME E1");

    // Add Legend
    TLegend *legend = new TLegend(0.7,0.68,0.9,0.9);  
    legend->AddEntry(sb.data, "Data");
    legend->AddEntry(sb.fit, "Fit Result");
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

    // Plot Output
    std::string plotDir = Folder_List::plotDir_SideBand;
    std::string out_name;
    out_name = plotDir + sb.name + "_" + type + "_" + norm + ".png"; 

    c->Print(out_name.c_str(),"png");
   
    delete legend;
    delete hs;
    delete c;
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
    
    // Fit
    sb.fit->SetLineColor(kMagenta);
    sb.fit->SetLineWidth(3);
    sb.fit->SetFillStyle(0);
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


void CCProtonPi0_SideBandTool::ApplyFitResults(SideBand &sb)
{
    // Clone Original Histograms
    sb.signal[1] = new TH1D (*sb.signal[0]);
    sb.WithPi0[1] = new TH1D (*sb.WithPi0[0]);
    sb.QELike[1] = new TH1D (*sb.QELike[0]);
    sb.SinglePiPlus[1] = new TH1D (*sb.SinglePiPlus[0]);
    sb.Other[1] = new TH1D (*sb.Other[0]);
    
    // Scale 
    sb.signal[1]->Scale(sb.fr_signal[1]/sb.fr_signal[0]);
    sb.WithPi0[1]->Scale(sb.fr_WithPi0[1]/sb.fr_WithPi0[0]);
    sb.QELike[1]->Scale(sb.fr_QELike[1]/sb.fr_QELike[0]);
    sb.SinglePiPlus[1]->Scale(sb.fr_SinglePiPlus[1]/sb.fr_SinglePiPlus[0]);
    sb.Other[1]->Scale(sb.fr_Other[1]/sb.fr_Other[0]);
}



