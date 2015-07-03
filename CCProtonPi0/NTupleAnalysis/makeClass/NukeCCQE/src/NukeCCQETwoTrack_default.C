#ifndef NukeCCQETwoTrack_default_cxx
#define NukeCCQETwoTrack_default_cxx

#include "NukeCCQETwoTrack_default.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#include <iostream>

using namespace std;

void NukeCCQETwoTrack_default::Loop()
{
    if (fChain == 0) return;

    TFile* f = new TFile("NukeCCQE_default.root","RECREATE");


    TH1D* open_angle = new TH1D("open_angle", "Openning Angle", 100 , 0,3.15);
    TH1D* proton_p_error = new TH1D("proton_p_error", "(Reco - True)/True", 100 , -2,2);
    TH1D* proton_p_reco = new TH1D("proton_p_reco", "Reconstructed Proton Momentum [MeV]", 41, -150.0,4000);
    TH1D* proton_p_mc = new TH1D("proton_p_mc", "True Proton Momentum [MeV]", 41, -150.0,4000);
    TH2D* proton_p_reco_mc = new TH2D("proton_p_reco_mc", "True vs Reconstructed Proton Momentum [MeV]", 41, -150.0,4000, 41, -150.0, 4000);

    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;
  
    double error;
    double p_reco;
    double p_mc;
    double NBad = 0;
    double NGood = 0;
    double NKinked = 0;
    double NNoKink = 0;
    double NKinkedAll = 0;
    double NNoKinkAll = 0;

    for (Long64_t jentry=0; jentry<nentries;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry);   nbytes += nb;
        // if (Cut(ientry) < 0) continue;

        if ( (jentry%10000)==0 ) cout<<mc_run<<" "<<mc_subrun<<endl;
       
        for ( int i = 0; i < 10; i++){
           
            if (NukeCCQETwoTrack_trajProtonProngPDG[i] == 2212){
                
                if (NukeCCQETwoTrack_proton_p[i] > 0 ){
                    p_reco = NukeCCQETwoTrack_proton_p[i]; 
                    NGood++;
                }else{
                    p_reco = -100.0;
                    NBad++;
                }
                
                p_mc = NukeCCQETwoTrack_trajProtonProngMomentum[i];            
                error = (p_reco - p_mc)/p_mc;

                proton_p_reco->Fill(p_reco);
                proton_p_mc->Fill(p_mc);
                proton_p_reco_mc->Fill(p_reco,p_mc);
                proton_p_error->Fill(error);
  
                if(NukeCCQETwoTrack_proton_kinked[i] == 1) NKinkedAll++;
                else NNoKinkAll++;

                if(error < -0.2 && error > -1){
                    open_angle->Fill(NukeCCQETwoTrack_open_angle[i]);
                    if(NukeCCQETwoTrack_proton_kinked[i] == 1) NKinked++;
                    else NNoKink++;
                }
            }
        }
    }

    double total = NGood + NBad;
    double percent = NGood / total * 100;

    cout<<"NGood = "<<NGood<<" NBad = "<<NBad<<endl;
    cout<<"Good Percent over Total = "<<percent<<endl;
    cout<<"NKinked = "<<NKinked<<" NNoKink = "<<NNoKink<<endl;
    cout<<"NKinkedAll = "<<NKinkedAll<<" NNoKinkAll = "<<NNoKinkAll<<endl;

    //plot1D_HistLogScale(proton_p_error,"proton_p_error.png","Plots/");
    //plot1D_Hist(proton_p_error,"modified_proton_p_error.png","Plots/");
    plot1D_Hist(open_angle,"modified_openning_angle.png","Plots/");
    //plot1D_Hist(proton_p_reco,"modified_proton_p_reco.png","Plots/");
    //plot1D_HistLogScale(proton_p_reco,"modified_proton_p_reco_log.png","Plots/");
    //
    //plot1D_Hist(proton_p_mc,"modified_proton_p_mc.png","Plots/");
    //
    //plot2D_Hist(proton_p_reco_mc,"modified_proton_p_reco_mc.png","Plots/");

    f->Write();

}

void NukeCCQETwoTrack_default::plot1D_HistLogScale(TH1D* hist1D, string fileName, string plotDir)
{
    TCanvas* c1 = new TCanvas();
    c1->SetLogy();
    
    hist1D->SetLineColor(kRed);
    hist1D->SetLineWidth(3);
    hist1D->SetFillColor(kRed);
    hist1D->SetFillStyle(3010);
    
    hist1D->Draw();
    gPad->Update();
    
    // Statistics Box
    //     TPaveStats *st = (TPaveStats*)hist1D->FindObject("stats");
    
    c1->Print(Form("%s%s",plotDir.c_str(),fileName.c_str()),"png");
    delete c1;
    
}

void NukeCCQETwoTrack_default::plot1D_Hist(TH1D* hist1D, string fileName, string plotDir)
{
    TCanvas* c1 = new TCanvas();
    hist1D->SetLineColor(kRed);
    hist1D->SetLineWidth(3);
    hist1D->SetFillColor(kRed);
    hist1D->SetFillStyle(3010);
    
    hist1D->Draw();
    gPad->Update();
    
    // Statistics Box
//     TPaveStats *st = (TPaveStats*)hist1D->FindObject("stats");
    
    c1->Print(Form("%s%s",plotDir.c_str(),fileName.c_str()),"png");
    delete c1;
    
}

void NukeCCQETwoTrack_default::plot2D_Hist(TH2D* hist2D, string fileName, string plotDir)
{
    // Canvas
    Double_t w = 800; 
    Double_t h = 800;
    TCanvas* c1 = new TCanvas("c","c",w,h);
    c1->SetWindowSize(w,h);
    
    // Pad
    TPad *p = new TPad("p","p",0.05,0.05,0.95,0.95);
    p->Draw();
    
    p->cd();
    hist2D->GetYaxis()->SetTitleOffset(1.8);
    hist2D->Draw("colz");
    gPad->Update();
    
    // Statistics Box
    TPaveStats *st = (TPaveStats*)hist2D->FindObject("stats");
    st->SetOptStat(1000000110);
    st->SetX1NDC(0.1); 
    st->SetX2NDC(0.3); 
    st->SetY1NDC(0.8); 
    st->SetY2NDC(0.9); 
   
    c1->Print(Form("%s%s",plotDir.c_str(),fileName.c_str()),"png");
    
    delete p;
    delete c1;
}

#endif

