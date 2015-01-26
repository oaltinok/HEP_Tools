#ifndef pID_Plots_cpp
#define pID_Plots_cpp

#include "Plotter.h"

void Plotter::plotPID()
{
    pID_proton();
    pID_pion();
    pID_proton_LLR();
    pID_pion_LLR();
    plot_2D_pID();
    pIDSum();
    pIDDiff();
}

void Plotter::pID_pion()
{
    string rootDir = rootDir_Interaction[branchInd];
    string plotDir = plotDir_Interaction[branchInd];
    
    cout<<"\nPlottting pID for Pion"<<endl;
    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());
    TCanvas* c1 = new TCanvas();
    THStack *hs = new THStack("hs","Pion Score");
    TLegend *legend = new TLegend(0.7,0.8,0.9,0.9); 
    
    TH1D* h_pID_piminus = (TH1D*)f_Root->Get("pID_piminus_pionScore");
    h_pID_piminus->SetFillColor(kYellow);
    h_pID_piminus->SetMarkerStyle(21);
    h_pID_piminus->SetMarkerColor(kYellow);
    
    TH1D* h_pID_piplus = (TH1D*)f_Root->Get("pID_piplus_pionScore");
    h_pID_piplus->SetFillColor(kBlue);
    h_pID_piplus->SetMarkerStyle(21);
    h_pID_piplus->SetMarkerColor(kBlue);
    
    TH1D* h_pID_proton = (TH1D*)f_Root->Get("pID_proton_pionScore");
    h_pID_proton->SetFillColor(kGreen);
    h_pID_proton->SetMarkerStyle(21);
    h_pID_proton->SetMarkerColor(kGreen);
    
    legend->AddEntry(h_pID_proton, "Proton", "f");
    legend->AddEntry(h_pID_piplus, "Pi-Plus", "f");
    legend->AddEntry(h_pID_piminus, "Pi-Minus", "f");
    
    
    hs->Add(h_pID_piminus);
    hs->Add(h_pID_piplus);
    hs->Add(h_pID_proton);
    hs->Draw();
    hs->GetXaxis()->SetTitle("Pion Score");
    hs->GetYaxis()->SetTitle("N(Events)");
    
    legend->Draw();
    
    c1->Print(Form("%s%s",plotDir.c_str(),"stacked_pID_pionScore.png"),"png");
    
    delete f_Root;
    delete hs;
    delete legend;
    delete c1;
    
}

void Plotter::pID_proton()
{
    string rootDir = rootDir_Interaction[branchInd];
    string plotDir = plotDir_Interaction[branchInd];
    
    cout<<"\nPlottting pID for Proton"<<endl;
    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());
    TCanvas* c1 = new TCanvas();
    THStack *hs = new THStack("hs","Proton Score");
    TLegend *legend = new TLegend(0.7,0.8,0.9,0.9); 
    
    TH1D* h_pID_piminus = (TH1D*)f_Root->Get("pID_piminus_protonScore");
    h_pID_piminus->SetFillColor(kYellow);
    h_pID_piminus->SetMarkerStyle(21);
    h_pID_piminus->SetMarkerColor(kYellow);
    
    TH1D* h_pID_piplus = (TH1D*)f_Root->Get("pID_piplus_protonScore");
    h_pID_piplus->SetFillColor(kBlue);
    h_pID_piplus->SetMarkerStyle(21);
    h_pID_piplus->SetMarkerColor(kBlue);
    
    TH1D* h_pID_proton = (TH1D*)f_Root->Get("pID_proton_protonScore");
    h_pID_proton->SetFillColor(kGreen);
    h_pID_proton->SetMarkerStyle(21);
    h_pID_proton->SetMarkerColor(kGreen);
    
    legend->AddEntry(h_pID_proton, "Proton", "f");
    legend->AddEntry(h_pID_piplus, "Pi-Plus", "f");
    legend->AddEntry(h_pID_piminus, "Pi-Minus", "f");
    
    
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

void Plotter::pIDSum()
{
    string rootDir = rootDir_Interaction[branchInd];
    string plotDir = plotDir_Interaction[branchInd];
    
    cout<<"\nPlottting pIDSum"<<endl;
    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());
    TCanvas* c1 = new TCanvas();
    THStack *hs = new THStack("hs","Proton Score + Pion Score");
    TLegend *legend = new TLegend(0.7,0.8,0.9,0.9); 
    
    TH1D* h_pID_piminus = (TH1D*)f_Root->Get("pID_piminus_pIDSum");
    h_pID_piminus->SetFillColor(kYellow);
    h_pID_piminus->SetMarkerStyle(21);
    h_pID_piminus->SetMarkerColor(kYellow);
    
    TH1D* h_pID_piplus = (TH1D*)f_Root->Get("pID_piplus_pIDSum");
    h_pID_piplus->SetFillColor(kBlue);
    h_pID_piplus->SetMarkerStyle(21);
    h_pID_piplus->SetMarkerColor(kBlue);
    
    TH1D* h_pID_proton = (TH1D*)f_Root->Get("pID_proton_pIDSum");
    h_pID_proton->SetFillColor(kGreen);
    h_pID_proton->SetMarkerStyle(21);
    h_pID_proton->SetMarkerColor(kGreen);
    
    legend->AddEntry(h_pID_proton, "Proton", "f");
    legend->AddEntry(h_pID_piplus, "Pi-Plus", "f");
    legend->AddEntry(h_pID_piminus, "Pi-Minus", "f");
    
    
    hs->Add(h_pID_piminus);
    hs->Add(h_pID_piplus);
    hs->Add(h_pID_proton);
    hs->Draw();
    hs->GetXaxis()->SetTitle("Proton Score + Pion Score");
    hs->GetYaxis()->SetTitle("N(Events)");
    
    legend->Draw();
    
    c1->Print(Form("%s%s",plotDir.c_str(),"stacked_pIDSum.png"),"png");
    
    delete f_Root;
    delete hs;
    delete legend;
    delete c1;
    
}

void Plotter::pIDDiff()
{
    string rootDir = rootDir_Interaction[branchInd];
    string plotDir = plotDir_Interaction[branchInd];
    
    cout<<"\nPlottting pIDDiff"<<endl;
    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());
    TCanvas* c1 = new TCanvas();
    THStack *hs = new THStack("hs","Proton Score - Pion Score");
    TLegend *legend = new TLegend(0.7,0.8,0.9,0.9); 
    
    TH1D* h_pID_piminus = (TH1D*)f_Root->Get("pID_piminus_pIDDiff");
    h_pID_piminus->SetFillColor(kYellow);
    h_pID_piminus->SetMarkerStyle(21);
    h_pID_piminus->SetMarkerColor(kYellow);
    
    TH1D* h_pID_piplus = (TH1D*)f_Root->Get("pID_piplus_pIDDiff");
    h_pID_piplus->SetFillColor(kBlue);
    h_pID_piplus->SetMarkerStyle(21);
    h_pID_piplus->SetMarkerColor(kBlue);
    
    TH1D* h_pID_proton = (TH1D*)f_Root->Get("pID_proton_pIDDiff");
    h_pID_proton->SetFillColor(kGreen);
    h_pID_proton->SetMarkerStyle(21);
    h_pID_proton->SetMarkerColor(kGreen);
    
    legend->AddEntry(h_pID_proton, "Proton", "f");
    legend->AddEntry(h_pID_piplus, "Pi-Plus", "f");
    legend->AddEntry(h_pID_piminus, "Pi-Minus", "f");
    
    
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
    string rootDir = rootDir_Interaction[branchInd];
    string plotDir = plotDir_Interaction[branchInd];
    
    cout<<"\nPlottting 2D pID Plots"<<endl;
    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());
    
    // Pion Score vs Proton Score
    TH2D* h_pID_proton_pionScore_protonScore = (TH2D*)f_Root->Get("pID_proton_pionScore_protonScore");
    plot2D_Hist(h_pID_proton_pionScore_protonScore,"pID_proton_pionScore_protonScore.png",plotDir);
    
    TH2D* h_pID_piminus_pionScore_protonScore = (TH2D*)f_Root->Get("pID_piminus_pionScore_protonScore");
    plot2D_Hist(h_pID_piminus_pionScore_protonScore,"pID_piminus_pionScore_protonScore.png",plotDir);
    
    TH2D* h_pID_piplus_pionScore_protonScore = (TH2D*)f_Root->Get("pID_piplus_pionScore_protonScore");
    plot2D_Hist(h_pID_piplus_pionScore_protonScore,"pID_piplus_pionScore_protonScore.png",plotDir);
    
    // Proton Score vs Proton Score LLR
    TH2D* h_pID_proton_protonScore_protonScore_LLR = (TH2D*)f_Root->Get("pID_proton_protonScore_protonScore_LLR");
    plot2D_Hist(h_pID_proton_protonScore_protonScore_LLR,"pID_proton_protonScore_protonScore_LLR.png",plotDir);
    
    TH2D* h_pID_piplus_protonScore_protonScore_LLR = (TH2D*)f_Root->Get("pID_piplus_protonScore_protonScore_LLR");
    plot2D_Hist(h_pID_piplus_protonScore_protonScore_LLR,"pID_piplus_protonScore_protonScore_LLR.png",plotDir);
    
    TH2D* h_pID_piminus_protonScore_protonScore_LLR = (TH2D*)f_Root->Get("pID_piminus_protonScore_protonScore_LLR");
    plot2D_Hist(h_pID_piminus_protonScore_protonScore_LLR,"pID_piminus_protonScore_protonScore_LLR.png",plotDir);
    
    // Pion Score vs Pion Score LLR
    TH2D* h_pID_proton_pionScore_pionScore_LLR = (TH2D*)f_Root->Get("pID_proton_pionScore_pionScore_LLR");
    plot2D_Hist(h_pID_proton_pionScore_pionScore_LLR,"pID_proton_pionScore_pionScore_LLR.png",plotDir);
    
    TH2D* h_pID_piplus_pionScore_pionScore_LLR = (TH2D*)f_Root->Get("pID_piplus_pionScore_pionScore_LLR");
    plot2D_Hist(h_pID_piplus_pionScore_pionScore_LLR,"pID_piplus_pionScore_pionScore_LLR.png",plotDir);
    
    TH2D* h_pID_piminus_pionScore_pionScore_LLR = (TH2D*)f_Root->Get("pID_piminus_pionScore_pionScore_LLR");
    plot2D_Hist(h_pID_piminus_pionScore_pionScore_LLR,"pID_piminus_pionScore_pionScore_LLR.png",plotDir);
     
}

void Plotter::pID_pion_LLR()
{
    string rootDir = rootDir_Interaction[branchInd];
    string plotDir = plotDir_Interaction[branchInd];
    
    cout<<"\nPlottting pID for Pion"<<endl;
    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());
    TCanvas* c1 = new TCanvas();
    THStack *hs = new THStack("hs","Pion Score Log-Likelihood Ratio");
    TLegend *legend = new TLegend(0.7,0.8,0.9,0.9); 
    
    TH1D* h_pID_piminus = (TH1D*)f_Root->Get("pID_piminus_pionScore_LLR");
    h_pID_piminus->SetFillColor(kYellow);
    h_pID_piminus->SetMarkerStyle(21);
    h_pID_piminus->SetMarkerColor(kYellow);
    
    TH1D* h_pID_piplus = (TH1D*)f_Root->Get("pID_piplus_pionScore_LLR");
    h_pID_piplus->SetFillColor(kBlue);
    h_pID_piplus->SetMarkerStyle(21);
    h_pID_piplus->SetMarkerColor(kBlue);
    
    TH1D* h_pID_proton = (TH1D*)f_Root->Get("pID_proton_pionScore_LLR");
    h_pID_proton->SetFillColor(kGreen);
    h_pID_proton->SetMarkerStyle(21);
    h_pID_proton->SetMarkerColor(kGreen);
    
    legend->AddEntry(h_pID_proton, "Proton", "f");
    legend->AddEntry(h_pID_piplus, "Pi-Plus", "f");
    legend->AddEntry(h_pID_piminus, "Pi-Minus", "f");
    
    
    hs->Add(h_pID_piminus);
    hs->Add(h_pID_piplus);
    hs->Add(h_pID_proton);
    hs->Draw();
    hs->GetXaxis()->SetTitle("Pion Score Log-Likelihood Ratio");
    hs->GetYaxis()->SetTitle("N(Events)");
    
    legend->Draw();
    
    c1->Print(Form("%s%s",plotDir.c_str(),"stacked_pID_pionScore_LLR.png"),"png");
    
    delete f_Root;
    delete hs;
    delete legend;
    delete c1;
    
}

void Plotter::pID_proton_LLR()
{
    string rootDir = rootDir_Interaction[branchInd];
    string plotDir = plotDir_Interaction[branchInd];
    
    cout<<"\nPlottting pID for Proton"<<endl;
    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());
    TCanvas* c1 = new TCanvas();
    THStack *hs = new THStack("hs","Proton Score Log-Likelihood Ratio");
    TLegend *legend = new TLegend(0.7,0.8,0.9,0.9); 
    
    TH1D* h_pID_piminus = (TH1D*)f_Root->Get("pID_piminus_protonScore_LLR");
    h_pID_piminus->SetFillColor(kYellow);
    h_pID_piminus->SetMarkerStyle(21);
    h_pID_piminus->SetMarkerColor(kYellow);
    
    TH1D* h_pID_piplus = (TH1D*)f_Root->Get("pID_piplus_protonScore_LLR");
    h_pID_piplus->SetFillColor(kBlue);
    h_pID_piplus->SetMarkerStyle(21);
    h_pID_piplus->SetMarkerColor(kBlue);
    
    TH1D* h_pID_proton = (TH1D*)f_Root->Get("pID_proton_protonScore_LLR");
    h_pID_proton->SetFillColor(kGreen);
    h_pID_proton->SetMarkerStyle(21);
    h_pID_proton->SetMarkerColor(kGreen);
    
    legend->AddEntry(h_pID_proton, "Proton", "f");
    legend->AddEntry(h_pID_piplus, "Pi-Plus", "f");
    legend->AddEntry(h_pID_piminus, "Pi-Minus", "f");
    
    
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

