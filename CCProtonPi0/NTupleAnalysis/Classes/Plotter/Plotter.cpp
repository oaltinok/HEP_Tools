/*
    See Plotter.h header for Class Information
*/

#ifndef PLOTTER_CPP
#define PLOTTER_CPP
#include "Plotter.h"

using namespace std;


void Plotter::plotHistograms()
{
    if (isSignalvsBackground){
//         plotSignalBackground();
//         plotCutHistograms();
    }else{
//         plotInteraction();
        
//         plotMuon();
        
//         plotProton();
        
//         plotPion();
        
//        plotPID();
        
//         plot_mc_w_Stacked();
//         plot_final_mc_w_Stacked();
            plotMichel();
    }
   
    
}

void Plotter::plotSignalRatio(TH1D* h_signal, TH1D* h_background, string fileName, string plotDir) 
{
    TH1D* h_ratio = new TH1D;
    
    h_signal->Copy(*h_ratio);
    
    h_ratio->Divide(h_background);
    h_ratio->GetYaxis()->SetTitle("Signal / Background");
    
    TCanvas* c1 = new TCanvas();
    h_ratio->SetLineColor(kBlack);
    h_ratio->SetLineWidth(2);
    
    h_ratio->Draw();
    gPad->Update();
     
    c1->Print(Form("%s%s%s",plotDir.c_str(),"Ratio_",fileName.c_str()),"png");
    
    delete c1;
    delete h_ratio;
}


void Plotter::inform(string rootDir, string plotDir)
{
    cout<<"------------ Plotting ------------"<<endl;
    cout<<"Input File: "<<rootDir<<endl;
    cout<<"Output Folder: "<<plotDir<<endl;
}

void Plotter::plotStackedLogScale(TH1D* h_signal, TH1D* h_background, string plotName, string fileName, string plotDir)
{
    TH1D* h_signalRatio = new TH1D;
    TH1D* h_backgroundRatio = new TH1D;
    
    h_signal->Copy(*h_signalRatio);
    h_background->Copy(*h_backgroundRatio);
    
    TCanvas *c1 = new TCanvas();
    THStack *hs = new THStack("hs",plotName.c_str());
    TLegend *legend = new TLegend(0.7,0.8,0.9,0.9);
    
    c1->SetLogy();
    
    h_signal->SetFillColor(kGreen);
    h_signal->SetLineColor(kGreen);
    h_signal->SetMarkerStyle(21);
    h_signal->SetMarkerColor(kGreen);
    
    h_background->SetFillColor(kRed);
    h_background->SetLineColor(kRed);
    h_background->SetMarkerStyle(21);
    h_background->SetMarkerColor(kRed);
    
    legend->AddEntry(h_signal, "Signal", "f");
    legend->AddEntry(h_background, "Background", "f");
    
    hs->Add(h_background);
    hs->Add(h_signal);
    hs->Draw();
    hs->GetXaxis()->SetTitle(plotName.c_str());
    hs->GetYaxis()->SetTitle("N(Events)");
    
    legend->Draw();
    
    c1->Print(Form("%s%s",plotDir.c_str(),fileName.c_str()),"png");
    
    delete hs;
    delete legend;
    delete c1;
    
    // Plot Signal Ratio
    plotSignalRatio(h_signalRatio,h_backgroundRatio, fileName, plotDir);
    
    
    delete h_signalRatio;
    delete h_backgroundRatio;
    
}

void Plotter::plotStacked(TH1D* h_signal, TH1D* h_background, 
                            string plotName, string fileName, string plotDir, 
                            string signal_label, string background_label)
{
    TH1D* h_signalRatio = new TH1D;
    TH1D* h_backgroundRatio = new TH1D;
    
    h_signal->Copy(*h_signalRatio);
    h_background->Copy(*h_backgroundRatio);
    
    TCanvas *c1 = new TCanvas();
    THStack *hs = new THStack("hs",plotName.c_str());
    TLegend *legend = new TLegend(0.7,0.8,0.9,0.9);
    
    
    h_signal->SetFillColor(kGreen);
    h_signal->SetLineColor(kGreen);
    h_signal->SetMarkerStyle(21);
    h_signal->SetMarkerColor(kGreen);
    
    h_background->SetFillColor(kRed);
    h_background->SetLineColor(kRed);
    h_background->SetMarkerStyle(21);
    h_background->SetMarkerColor(kRed);
    
    legend->AddEntry(h_signal, signal_label.c_str(), "f");
    legend->AddEntry(h_background, background_label.c_str(), "f");
 
    
    hs->Add(h_background);
    hs->Add(h_signal);
    hs->Draw();
    hs->GetXaxis()->SetTitle(plotName.c_str());
    hs->GetYaxis()->SetTitle("N(Events)");
    
    legend->Draw();
  
    c1->Print(Form("%s%s",plotDir.c_str(),fileName.c_str()),"png");
    
    delete hs;
    delete legend;
    delete c1;
    
    // Plot Signal Ratio
    plotSignalRatio(h_signalRatio,h_backgroundRatio,fileName, plotDir);

    delete h_signalRatio;
    delete h_backgroundRatio;
}

void Plotter::plot1D_HistLogScale(TH1D* hist1D, string fileName, string plotDir)
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


void Plotter::plot1D_Hist(TH1D* hist1D, string fileName, string plotDir)
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

void Plotter::plot2D_Hist(TH2D* hist2D, string fileName, string plotDir)
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

//------------------------------------------------------------------------------
// Analysis Mode
//      1) Signal Events
//      2) Background Events
//      Other) All Events
//------------------------------------------------------------------------------
Plotter::Plotter(int nMode)
{
    isSignalvsBackground = false;
    
    if ( nMode == 1) {
        cout<<"----------------------------------------------------------------------"<<endl;
        cout<<"Plot Mode: Signal - Only Signal Events will be Plotted"<<endl;
        branchInd = 0;
    }else if ( nMode == 2){
        cout<<"----------------------------------------------------------------------"<<endl;
        cout<<"Plot Mode: Background - Only Background Events will be Plotted"<<endl;
        branchInd = 1;
    }else if (nMode == 3){
        cout<<"----------------------------------------------------------------------"<<endl;
        cout<<"Plot Mode: All - All Events will be Plotted"<<endl;
        branchInd = 2;
    }else{
        cout<<"----------------------------------------------------------------------"<<endl;
        cout<<"Plot Mode: Signal vs Background"<<endl;
        isSignalvsBackground = true;
    }
    
    cout<<"----------------------------------------------------------------------"<<endl;
    
    setFolders();
}

void Plotter::setFolders()
{
    for ( int i = 0 ; i < nBranches; i++){
        
        // Set Branch: 0 = Signal, 1 = Background, Other = All Events
        if ( i == 0 ) branchDir = Folder_List::signal;
        else if ( i == 1) branchDir = Folder_List::background;
        else branchDir = Folder_List::allEvents;
        
        rootDir_Interaction[i] = Folder_List::output + Folder_List::rootOut + branchDir + "Interaction.root";
        plotDir_Interaction[i] = Folder_List::output + Folder_List::plotOut + branchDir + "Interaction/";
        
        rootDir_PID[i] = Folder_List::output + Folder_List::rootOut + branchDir + "PIDStatistics.root";
        plotDir_PID[i] = Folder_List::output + Folder_List::plotOut + branchDir + "PIDStatistics/";

        rootDir_Muon[i] = Folder_List::output + Folder_List::rootOut + branchDir + "Muon.root";
        plotDir_Muon[i] = Folder_List::output + Folder_List::plotOut + branchDir + "Muon/";

        rootDir_Proton[i] = Folder_List::output + Folder_List::rootOut + branchDir + "Proton.root";
        plotDir_Proton[i] = Folder_List::output + Folder_List::plotOut + branchDir + "Proton/";

        rootDir_Pion[i] = Folder_List::output + Folder_List::rootOut + branchDir + "Pion.root";
        plotDir_Pion[i] = Folder_List::output + Folder_List::plotOut + branchDir + "Pion/";
        
    }
    
    otherDir = Folder_List::output + Folder_List::plotOut + Folder_List::other;
}

#endif





