/*
================================================================================
Study: pID_Studies
    Stand Alone ROOT Function that reads a specific .ROOT file and generates
    specific output related with the study.
    
    pID_Studies reads the ROOT file generated for CCDeltaPlus Class
    
    
    Main Directory:
        Studies/
        
    Usage:
        > root -l
        > .L pID_Studies
        > pID_Studies()
    
    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision:  2014_05_02
================================================================================
*/

#include <TH1.h>
#include <TH2.h>
#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h>

#include "../Libraries/Folder_List.h" // Contains File and Folder Paths

void pID_Studies() 
{
    string rootDir = getFileLocation("../",Folder_List::f_Root_CCDeltaPlus);
    string plotDir = getFileLocation("../",Folder_List::f_Plot_CCDeltaPlus);
    
    inform(rootDir, plotDir);
    
    TFile* f_Root = new TFile(rootDir.c_str());
    
    
    get_pID_Stats(f_Root);
    get_pID_Plots(f_Root, plotDir);

}

void get_pID_Plots(TFile* f_Root, string plotDir)
{
    THStack *hs = new THStack("hs","Proton Score");
    
    TH1F* h_pID_other = f_Root->Get("pID_other");
    h_pID_other->SetFillColor(kRed);
    h_pID_other->SetMarkerStyle(21);
    h_pID_other->SetMarkerColor(kRed);
    
    TH1F* h_pID_piminus = f_Root->Get("pID_piminus");
    h_pID_piminus->SetFillColor(kYellow);
    h_pID_piminus->SetMarkerStyle(21);
    h_pID_piminus->SetMarkerColor(kYellow);
    
    TH1F* h_pID_piplus = f_Root->Get("pID_piplus");
    h_pID_piplus->SetFillColor(kBlue);
    h_pID_piplus->SetMarkerStyle(21);
    h_pID_piplus->SetMarkerColor(kBlue);
    
    TH1F* h_pID_proton = f_Root->Get("pID_proton");
    h_pID_proton->SetFillColor(kGreen);
    h_pID_proton->SetMarkerStyle(21);
    h_pID_proton->SetMarkerColor(kGreen);
    
    TCanvas* c1 = new TCanvas();

    hs->Add(h_pID_other);
    hs->Add(h_pID_piminus);
    hs->Add(h_pID_piplus);
    hs->Add(h_pID_proton);
    hs->Draw();
    hs->GetXaxis()->SetTitle("Proton Score");
    hs->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",0.05));
    
    
    c1->Print(Form("%s%s",plotDir.c_str(),"stacked_pID"),"png");
}


void get_pID_Stats(TFile* f_Root)
{
    
    TH1F* pID_proton  = f_Root->Get("pID_proton");
    TH1F* pID_piplus  = f_Root->Get("pID_piplus");
    TH1F* pID_piminus = f_Root->Get("pID_piminus");
    TH1F* pID_other   = f_Root->Get("pID_other");

    double nProton = 0;
    double nTotalProton = 0;
    double nCapturedEvents = 0;
    double nEvents = 0;
    double purity;
    double efficiency;
    int nBins = 20;
    
    //Get Total Proton
    for(int i = nBins; i >= 1; i--){
        
    }
    cout<<"nTotalProton = "<<nTotalProton<<endl;
    
    //Get Total Proton
    for(int i = nBins; i >= 1; i--){
        nTotalProton = nTotalProton + pID_proton->GetBinContent(i);
    }
    cout<<"nTotalProton = "<<nTotalProton<<endl;
    
    
    for(int i = nBins; i >= 1; i--){
        nProton = nProton + pID_proton->GetBinContent(i);
        nCapturedEvents =   nCapturedEvents+
                            pID_proton->GetBinContent(i) +
                            pID_piplus->GetBinContent(i) +
                            pID_piminus->GetBinContent(i) +
                            pID_other->GetBinContent(i);
        nEvents =           pID_proton->GetBinContent(i) +
                            pID_piplus->GetBinContent(i) +
                            pID_piminus->GetBinContent(i) +
                            pID_other->GetBinContent(i);
                            
        purity = nProton / nCapturedEvents;
        efficiency = nProton / nTotalProton;
//         cout<<"pID = "<<pID_proton->GetBinLowEdge(i)<<" Purity = "<<purity<<" Efficiency = "<<efficiency<<endl;
        cout<<pID_proton->GetBinLowEdge(i)<<" "<<purity<<" "<<efficiency<<endl;
    }
}

void inform(string rootDir, string plotDir)
{
    cout<<"------------ Plotting ------------"<<endl;
    cout<<"Input File: "<<rootDir<<endl;
    cout<<"Output Folder: "<<plotDir<<endl;

}