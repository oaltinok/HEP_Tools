#define Other_Sample_cxx
#include "Other_Sample.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void Other_Sample::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L Other_Sample.C
//      Root > Other_Sample t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
    
    reco_P = new TH1D( "reco_P","Reconstructed Proton Momentum",20,0.0,2000 );
    reco_P->GetXaxis()->SetTitle("Reconstructed P_{p} [MeV]");
    reco_P->GetYaxis()->SetTitle("Events");
    
    reco_P_true_P = new TH2D( "reco_P_true_P","True vs Reconstructed Proton Momentum",20,0.0,2000,20,0.0,2000);
    reco_P_true_P->GetXaxis()->SetTitle("Reconstructed P_{p} [MeV]");
    reco_P_true_P->GetYaxis()->SetTitle("True P_{p} [MeV]");

    P_error = new TH1D( "P_error","Error on Proton Momentum",400,-2,2 );
    P_error->GetXaxis()->SetTitle("(P_{Reco}-P_{True})/P_{True}");
    P_error->GetYaxis()->SetTitle("Events");

    double reco_p;
    double true_p;
    double error;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
        
        for(int i = 0; i < 10; i++){
            if (CCProtonPi0_trajProtonProngPDG[i] == 2212){
                reco_p = CCProtonPi0_all_protons_P[i];
                true_p = CCProtonPi0_trajProtonProngMomentum[i];
                error = (reco_p-true_p)/true_p;

                reco_P->Fill(reco_p);
                reco_P_true_P->Fill(reco_p,true_p);
                P_error->Fill(error);
            }
        } 
   }


    plot1D_Hist(P_error, "binary_P_error.png","Plots/");
    plot1D_Hist(reco_P, "binary_reco_P.png","Plots/");
    plot2D_Hist(reco_P_true_P, "binary_reco_P_true_P.png","Plots/");


}

void Other_Sample::plot1D_Hist(TH1D* hist1D, string fileName, string plotDir)
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

void Other_Sample::plot2D_Hist(TH2D* hist2D, string fileName, string plotDir)
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
