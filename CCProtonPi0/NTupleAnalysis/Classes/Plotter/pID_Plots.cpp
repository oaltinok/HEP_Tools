#ifndef pID_Plots_cpp
#define pID_Plots_cpp

#include "Plotter.h"

void Plotter::plotPID()
{
    pID_proton();
    pID_proton_LLR();
    plot_2D_pID();
    pIDDiff();
    pIDStats();
    KE();
}

void Plotter::KE()
{
    string rootDir = rootDir_PID[branchInd];
    string plotDir = plotDir_PID[branchInd];
    
    cout<<"\nPlottting Kinetic Energy Plots"<<endl;
    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());
    
    TH1D* h_KE_proton_pIDDiff = (TH1D*)f_Root->Get("KE_proton_pIDDiff");
    TH1D* h_KE_other_pIDDiff = (TH1D*)f_Root->Get("KE_other_pIDDiff");
    plotStacked(h_KE_proton_pIDDiff , h_KE_other_pIDDiff,"KE of True Protons(Green) and Other Particles(Red) using pIDDiff", "KE_proton_pIDDiff.png", plotDir);
    
    TH1D* h_KE_proton_LLR = (TH1D*)f_Root->Get("KE_proton_LLR");
    TH1D* h_KE_other_LLR = (TH1D*)f_Root->Get("KE_other_LLR");
    plotStacked(h_KE_proton_LLR , h_KE_other_LLR,"KE of True Protons(Green) and Other Particles(Red) using LLR", "KE_proton_LLR.png", plotDir);
}


void Plotter::pIDStats()
{
    string rootDir = rootDir_PID[branchInd];
    string plotDir = plotDir_PID[branchInd];
    
    cout<<"\nPlottting pID Statistics"<<endl;
    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());
    
    TH1D* h_purity_LLR = (TH1D*)f_Root->Get("purity_LLR");
    plot1D_Hist(h_purity_LLR,"purity_LLR.png",plotDir);
    
    TH1D* h_efficiency_LLR = (TH1D*)f_Root->Get("efficiency_LLR");
    plot1D_Hist(h_efficiency_LLR,"efficiency_LLR.png",plotDir);
    
    TH1D* h_purityXefficiency_LLR = (TH1D*)f_Root->Get("purityXefficiency_LLR");
    plot1D_Hist(h_purityXefficiency_LLR,"purityXefficiency_LLR.png",plotDir);
    
    TH1D* h_purity_pIDDiff = (TH1D*)f_Root->Get("purity_pIDDiff");
    plot1D_Hist(h_purity_pIDDiff,"purity_pIDDiff.png",plotDir);
    
    TH1D* h_efficiency_pIDDiff = (TH1D*)f_Root->Get("efficiency_pIDDiff");
    plot1D_Hist(h_efficiency_pIDDiff,"efficiency_pIDDiff.png",plotDir);
    
    TH1D* h_purityXefficiency_pIDDiff = (TH1D*)f_Root->Get("purityXefficiency_pIDDiff");
    plot1D_Hist(h_purityXefficiency_pIDDiff,"purityXefficiency_pIDDiff.png",plotDir);
}

void Plotter::pID_proton()
{
    string rootDir = rootDir_PID[branchInd];
    string plotDir = plotDir_PID[branchInd];
    
    cout<<"\nPlottting pID for Proton"<<endl;
    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());
    TCanvas* c1 = new TCanvas();
    THStack *hs = new THStack("hs","Proton Score");
    TLegend *legend = new TLegend(0.7,0.8,0.9,0.9); 
    
//     TH1D* h_pID_other = (TH1D*)f_Root->Get("other_protonScore");
//     h_pID_other->SetFillColor(kRed);
//     h_pID_other->SetMarkerStyle(21);
//     h_pID_other->SetMarkerColor(kRed);
    
    TH1D* h_pID_piminus = (TH1D*)f_Root->Get("piminus_protonScore");
    h_pID_piminus->SetFillColor(kYellow);
    h_pID_piminus->SetMarkerStyle(21);
    h_pID_piminus->SetMarkerColor(kYellow);
    
    TH1D* h_pID_piplus = (TH1D*)f_Root->Get("piplus_protonScore");
    h_pID_piplus->SetFillColor(kBlue);
    h_pID_piplus->SetMarkerStyle(21);
    h_pID_piplus->SetMarkerColor(kBlue);
    
    TH1D* h_pID_proton = (TH1D*)f_Root->Get("proton_protonScore");
    h_pID_proton->SetFillColor(kGreen);
    h_pID_proton->SetMarkerStyle(21);
    h_pID_proton->SetMarkerColor(kGreen);
    
    legend->AddEntry(h_pID_proton, "Proton", "f");
    legend->AddEntry(h_pID_piplus, "Pi-Plus", "f");
    legend->AddEntry(h_pID_piminus, "Pi-Minus", "f");
//     legend->AddEntry(h_pID_other, "Other", "f");
    
    
//     hs->Add(h_pID_other);
    hs->Add(h_pID_piminus);
    hs->Add(h_pID_piplus);
    hs->Add(h_pID_proton);
    hs->Draw();
    hs->GetXaxis()->SetTitle("Proton Score");
    hs->GetYaxis()->SetTitle("N(Events)");
    
    legend->Draw();
    
    c1->Print(Form("%s%s",plotDir.c_str(),"stacked_pID_protonScore.png"),"png");
    
    delete f_Root;
    delete hs;
    delete legend;
    delete c1;
    
}

void Plotter::pIDDiff()
{
    string rootDir = rootDir_PID[branchInd];
    string plotDir = plotDir_PID[branchInd];
    
    cout<<"\nPlottting pIDDiff"<<endl;
    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());
    TCanvas* c1 = new TCanvas();
    THStack *hs = new THStack("hs","Proton Score - Pion Score");
    TLegend *legend = new TLegend(0.7,0.8,0.9,0.9); 
    
    TH1D* h_pID_other = (TH1D*)f_Root->Get("other_pIDDiff");
    h_pID_other->SetFillColor(kRed);
    h_pID_other->SetMarkerStyle(21);
    h_pID_other->SetMarkerColor(kRed);
    
    TH1D* h_pID_piminus = (TH1D*)f_Root->Get("piminus_pIDDiff");
    h_pID_piminus->SetFillColor(kYellow);
    h_pID_piminus->SetMarkerStyle(21);
    h_pID_piminus->SetMarkerColor(kYellow);
    
    TH1D* h_pID_piplus = (TH1D*)f_Root->Get("piplus_pIDDiff");
    h_pID_piplus->SetFillColor(kBlue);
    h_pID_piplus->SetMarkerStyle(21);
    h_pID_piplus->SetMarkerColor(kBlue);
    
    TH1D* h_pID_proton = (TH1D*)f_Root->Get("proton_pIDDiff");
    h_pID_proton->SetFillColor(kGreen);
    h_pID_proton->SetMarkerStyle(21);
    h_pID_proton->SetMarkerColor(kGreen);
    
    legend->AddEntry(h_pID_proton, "Proton", "f");
    legend->AddEntry(h_pID_piplus, "Pi-Plus", "f");
    legend->AddEntry(h_pID_piminus, "Pi-Minus", "f");
    legend->AddEntry(h_pID_other, "Other", "f");
    
    hs->Add(h_pID_other);
    hs->Add(h_pID_piminus);
    hs->Add(h_pID_piplus);
    hs->Add(h_pID_proton);
    hs->Draw();
    hs->GetXaxis()->SetTitle("Proton Score - Pion Score");
    hs->GetYaxis()->SetTitle("N(Events)");
    
    legend->Draw();
    
    c1->Print(Form("%s%s",plotDir.c_str(),"stacked_pIDDiff.png"),"png");
    
    delete f_Root;
    delete hs;
    delete legend;
    delete c1;
    
}

void Plotter::plot_2D_pID()
{
    string rootDir = rootDir_PID[branchInd];
    string plotDir = plotDir_PID[branchInd];
    
    cout<<"\nPlottting 2D pID Plots"<<endl;
    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());
    
    // Pion Score vs Proton Score
    TH2D* h_pID_proton_pionScore_protonScore = (TH2D*)f_Root->Get("proton_pionScore_protonScore");
    plot2D_Hist(h_pID_proton_pionScore_protonScore,"pID_proton_pionScore_protonScore.png",plotDir);
    
    TH2D* h_pID_piminus_pionScore_protonScore = (TH2D*)f_Root->Get("piminus_pionScore_protonScore");
    plot2D_Hist(h_pID_piminus_pionScore_protonScore,"pID_piminus_pionScore_protonScore.png",plotDir);
    
    TH2D* h_pID_piplus_pionScore_protonScore = (TH2D*)f_Root->Get("piplus_pionScore_protonScore");
    plot2D_Hist(h_pID_piplus_pionScore_protonScore,"pID_piplus_pionScore_protonScore.png",plotDir);
    
    // Proton Score vs Proton Score LLR
    TH2D* h_pID_proton_protonScore_protonScore_LLR = (TH2D*)f_Root->Get("proton_protonScore_protonScore_LLR");
    plot2D_Hist(h_pID_proton_protonScore_protonScore_LLR,"pID_proton_protonScore_protonScore_LLR.png",plotDir);
    
    TH2D* h_pID_piplus_protonScore_protonScore_LLR = (TH2D*)f_Root->Get("piplus_protonScore_protonScore_LLR");
    plot2D_Hist(h_pID_piplus_protonScore_protonScore_LLR,"pID_piplus_protonScore_protonScore_LLR.png",plotDir);
    
    TH2D* h_pID_piminus_protonScore_protonScore_LLR = (TH2D*)f_Root->Get("piminus_protonScore_protonScore_LLR");
    plot2D_Hist(h_pID_piminus_protonScore_protonScore_LLR,"pID_piminus_protonScore_protonScore_LLR.png",plotDir);
     
}


void Plotter::pID_proton_LLR()
{
    string rootDir = rootDir_PID[branchInd];
    string plotDir = plotDir_PID[branchInd];
    
    cout<<"\nPlottting pID for Proton"<<endl;
    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());
    TCanvas* c1 = new TCanvas();
    THStack *hs = new THStack("hs","Proton Score Log-Likelihood Ratio");
    TLegend *legend = new TLegend(0.7,0.8,0.9,0.9); 
    
    TH1D* h_pID_other = (TH1D*)f_Root->Get("other_protonScore_LLR");
    h_pID_other->SetFillColor(kRed);
    h_pID_other->SetMarkerStyle(21);
    h_pID_other->SetMarkerColor(kRed);
    
    TH1D* h_pID_piminus = (TH1D*)f_Root->Get("piminus_protonScore_LLR");
    h_pID_piminus->SetFillColor(kYellow);
    h_pID_piminus->SetMarkerStyle(21);
    h_pID_piminus->SetMarkerColor(kYellow);
    
    TH1D* h_pID_piplus = (TH1D*)f_Root->Get("piplus_protonScore_LLR");
    h_pID_piplus->SetFillColor(kBlue);
    h_pID_piplus->SetMarkerStyle(21);
    h_pID_piplus->SetMarkerColor(kBlue);
    
    TH1D* h_pID_proton = (TH1D*)f_Root->Get("proton_protonScore_LLR");
    h_pID_proton->SetFillColor(kGreen);
    h_pID_proton->SetMarkerStyle(21);
    h_pID_proton->SetMarkerColor(kGreen);
    

    
    legend->AddEntry(h_pID_proton, "Proton", "f");
    legend->AddEntry(h_pID_piplus, "Pi-Plus", "f");
    legend->AddEntry(h_pID_piminus, "Pi-Minus", "f");
    legend->AddEntry(h_pID_other, "Other", "f");
    
    hs->Add(h_pID_other);
    hs->Add(h_pID_piminus);
    hs->Add(h_pID_piplus);
    hs->Add(h_pID_proton);
    hs->Draw();
    hs->GetXaxis()->SetTitle("Proton Score Log-Likelihood Ratio");
    hs->GetYaxis()->SetTitle("N(Events)");
    
    legend->Draw();
    
    c1->Print(Form("%s%s",plotDir.c_str(),"stacked_pID_protonScore_LLR.png"),"png");
    
    delete f_Root;
    delete hs;
    delete legend;
    delete c1;
    
}



#endif

