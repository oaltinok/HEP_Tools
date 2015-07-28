/*
    See CCProtonPi0_Plotter.h header for Class Information
*/

#ifndef CCProtonPi0_Plotter_Macros_cpp
#define CCProtonPi0_Plotter_Macros_cpp

#include "CCProtonPi0_Plotter.h"

using namespace PlotUtils;

void CCProtonPi0_Plotter::DrawStackedMC(std::string root_dir, std::string var_name, std::string plotDir)
{
   DrawStackedMC_BckgAll(root_dir,var_name,plotDir); 
   DrawStackedMC_BckgWithPi0(root_dir,var_name,plotDir); 
   DrawStackedMC_BckgType(root_dir,var_name,plotDir); 
}

void CCProtonPi0_Plotter::DrawStackedMC_BckgAll(std::string root_dir, std::string var_name, std::string plotDir)
{
    TFile* f_mc = new TFile(root_dir.c_str());

    // ------------------------------------------------------------------------
    // Fill TObjArray
    // ------------------------------------------------------------------------
    TObjArray* mc_hists = new TObjArray;
    MnvH1D* temp;
    std::string var;
   
    // Get Signal
    var = Form("%s_%d",var_name.c_str(),0);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Signal");
    mc_hists->Add(temp);
    
    // Get All Background
    var = Form("%s_%d",var_name.c_str(),1);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Background");
    mc_hists->Add(temp);
    
    
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas();
    plotter->DrawStackedMC(mc_hists,1.0 ,"TR");
    
    // Print Plot
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),"_bckg_all.png"), "png");

    delete c;
}

void CCProtonPi0_Plotter::DrawStackedMC_BckgType(std::string root_dir, std::string var_name, std::string plotDir)
{
    TFile* f_mc = new TFile(root_dir.c_str());

    // ------------------------------------------------------------------------
    // Fill TObjArray
    // ------------------------------------------------------------------------
    TObjArray* mc_hists = new TObjArray;
    MnvH1D* temp;
    std::string var;
   
    // Get Signal
    var = Form("%s_%d",var_name.c_str(),0);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Signal");
    mc_hists->Add(temp);
    
    // Get Bckg: NC
    var = Form("%s_%d",var_name.c_str(),5);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: NC");
    mc_hists->Add(temp);
    
    // Get Bckg: AntiNeutrino
    var = Form("%s_%d",var_name.c_str(),6);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: AntiNeutrino");
    mc_hists->Add(temp);

    // Get Bckg: QELike
    var = Form("%s_%d",var_name.c_str(),7);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: QELike");
    mc_hists->Add(temp);
   
    // Get Bckg: SinglePion
    var = Form("%s_%d",var_name.c_str(),8);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: SinglePion");
    mc_hists->Add(temp);

    // Get Bckg: DoublePion
    var = Form("%s_%d",var_name.c_str(),9);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: DoublePion");
    mc_hists->Add(temp);

    // Get Bckg: MultiPion
    var = Form("%s_%d",var_name.c_str(),10);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: MultiPion");
    mc_hists->Add(temp);
 
    // Get Bckg: Other
    var = Form("%s_%d",var_name.c_str(),11);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: Other");
    mc_hists->Add(temp);
    
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas("c","c",1280,800);
    plotter->DrawStackedMC(mc_hists,1.0 ,"TR");
    plotter->AddCutArrow(1.0, 0.0, 400, 0.4, "R"); 
    PlotUtils::t_PlotStyle style = kCompactStyle;    
    plotter->ApplyStyle(style); 
    
    // Print Plot
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),"_bckg_Type.png"), "png");

    delete c;
}
void CCProtonPi0_Plotter::DrawStackedMC_BckgWithPi0(std::string root_dir, std::string var_name, std::string plotDir)
{
    TFile* f_mc = new TFile(root_dir.c_str());

    // ------------------------------------------------------------------------
    // Fill TObjArray
    // ------------------------------------------------------------------------
    TObjArray* mc_hists = new TObjArray;
    MnvH1D* temp;
    std::string var;
   
    // Get Signal
    var = Form("%s_%d",var_name.c_str(),0);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Signal");
    mc_hists->Add(temp);
    
    // Get Bckg: NoPi0
    var = Form("%s_%d",var_name.c_str(),2);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: NoPi0");
    mc_hists->Add(temp);
    
    // Get Bckg: SinglePi0
    var = Form("%s_%d",var_name.c_str(),3);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: SinglePi0");
    mc_hists->Add(temp);

    // Get Bckg: MultiPi0
    var = Form("%s_%d",var_name.c_str(),4);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: MultiPi0");
    mc_hists->Add(temp);
    
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas();
    plotter->DrawStackedMC(mc_hists,1.0 ,"TR");
    t_PlotStyle style = kCCCohStyle;    
    plotter->ApplyStyle(style); 
    
    // Print Plot
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),"_bckg_Pi0.png"), "png");

    delete c;
}
void CCProtonPi0_Plotter::DrawDataMCRatio(rootDir& dir, std::string var_name, std::string plotDir)
{
    std::string rootDir_mc = dir.mc_all;
    std::string rootDir_data = dir.data;
    
    TFile* f_mc = new TFile(rootDir_mc.c_str());
    TFile* f_data = new TFile(rootDir_data.c_str());
   
    MnvH1D* mc = (MnvH1D*)f_mc->Get(var_name.c_str());
    MnvH1D* data = (MnvH1D*)f_data->Get(var_name.c_str()); 

    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas();
 
    // POT Normalize
    const int POT_exp = 19;
    const double POT_mc = 98.94; // E+19
    const double POT_data = 9.56; // E+19
    double weight = POT_data/POT_mc;
   
    plotter->DrawDataMCRatio(data, mc, weight);
    
    // Add Plot Labels
    const double y_pos = 0.88;
    const double y_diff = 0.033;
    //plotter->AddPlotLabel("MINERvA Preliminary",0.3,y_pos,y_diff,kBlue);
    plotter->AddPlotLabel("Playlist: minerva1",0.3,y_pos,y_diff,kBlue);
    plotter->AddPlotLabel("POT Normalized",0.3,y_pos-2*y_diff, y_diff); 
    plotter->AddPlotLabel(Form("%3.2fE%d Data POT",POT_data,POT_exp),0.3,y_pos-3*y_diff,y_diff);

    // Print Plot
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),"_ratio.png"), "png");

    delete c;
}

void CCProtonPi0_Plotter::DrawDataMC(rootDir& dir, std::string var_name, std::string plotDir)
{
    std::string rootDir_mc = dir.mc_all;
    std::string rootDir_data = dir.data;
    
    TFile* f_mc = new TFile(rootDir_mc.c_str());
    TFile* f_data = new TFile(rootDir_data.c_str());
   
    MnvH1D* mc = (MnvH1D*)f_mc->Get(var_name.c_str());
    MnvH1D* data = (MnvH1D*)f_data->Get(var_name.c_str()); 

    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas();
 
    // POT Normalize
    const int POT_exp = 19;
    const double POT_mc = 98.94; // E+19
    const double POT_data = 9.56; // E+19
    double weight = POT_data/POT_mc;
   
    plotter->DrawDataMC(data, mc, weight, "TR", false);
    
    // Add Plot Labels
    const double y_pos = 0.88;
    const double y_diff = 0.033;
    //plotter->AddPlotLabel("MINERvA Preliminary",0.3,y_pos,y_diff,kBlue);
    plotter->AddPlotLabel("Playlist: minerva1",0.3,y_pos,y_diff,kBlue);
    plotter->AddPlotLabel("POT Normalized",0.3,y_pos-2*y_diff, y_diff); 
    plotter->AddPlotLabel(Form("%3.2fE%d Data POT",POT_data,POT_exp),0.3,y_pos-3*y_diff,y_diff);

    // Print Plot
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".png"), "png");

    delete c;
    
    // Now Plot the Ratio
    DrawDataMCRatio(dir, var_name, plotDir);
}

void CCProtonPi0_Plotter::Draw1DHist(rootDir& dir, std::string var_name, std::string plotDir, bool isLogScale)
{
    std::string root_dir;

    if (branchInd == 0) root_dir = dir.mc_signal;
    else if (branchInd == 1) root_dir = dir.mc_background;
    else root_dir = dir.mc_all;
   
    // Get Histogram
    TFile* f = new TFile(root_dir.c_str());
    TH1D* hist1D = (TH1D*)f->Get(var_name.c_str());

    // Create Canvas
    TCanvas* c = new TCanvas();
    if(isLogScale) c->SetLogy();

    // Plot Options
    hist1D->SetLineColor(kRed);
    hist1D->SetLineWidth(3);
    hist1D->SetFillColor(kRed);
    hist1D->SetFillStyle(3010);
    
    hist1D->Draw();
    gPad->Update();
    
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".png"), "png");
    
    delete c;
    
}

void CCProtonPi0_Plotter::Draw2DHist(rootDir& dir, std::string var_name, std::string plotDir)
{
    std::string root_dir;

    if (branchInd == 0) root_dir = dir.mc_signal;
    else if (branchInd == 1) root_dir = dir.mc_background;
    else root_dir = dir.mc_all;
   
    // Get Histogram
    TFile* f = new TFile(root_dir.c_str());
    TH2D* hist2D = (TH2D*)f->Get(var_name.c_str());
    
    // Canvas
    Double_t w = 800; 
    Double_t h = 800;
    TCanvas* c = new TCanvas("c","c",w,h);
    c->SetWindowSize(w,h);
    
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
   
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".png"), "png");
    
    delete p;
    delete c;
}

//
//void CCProtonPi0_Plotter::plotSignalRatio(TH1D* h_signal, TH1D* h_background, std::string fileName, std::string plotDir, bool isReversed) 
//{
//    TH1D* h_ratio = new TH1D;
//    
//    if (isReversed){
//        h_background->Copy(*h_ratio);
//        
//        h_ratio->Divide(h_signal);
//        h_ratio->GetYaxis()->SetTitle("Background / Signal");
//    }else{
//        h_signal->Copy(*h_ratio);
//        
//        h_ratio->Divide(h_background);
//        h_ratio->GetYaxis()->SetTitle("Signal / Background");
//    }
//    
//    TCanvas* c1 = new TCanvas();
//    h_ratio->SetLineColor(kBlack);
//    h_ratio->SetLineWidth(2);
//    
//    h_ratio->Draw();
//    gPad->Update();
//     
//    c1->Print(Form("%s%s%s",plotDir.c_str(),"Ratio_",fileName.c_str()),"png");
//    
//    delete c1;
//    delete h_ratio;
//}
//
//void CCProtonPi0_Plotter::plotStackedLogScale(TH1D* h_signal, TH1D* h_background, std::string plotName, std::string fileName, std::string plotDir)
//{
//    TH1D* h_signalRatio = new TH1D;
//    TH1D* h_backgroundRatio = new TH1D;
//    
//    h_signal->Copy(*h_signalRatio);
//    h_background->Copy(*h_backgroundRatio);
//    
//    TCanvas *c1 = new TCanvas();
//    THStack *hs = new THStack("hs",plotName.c_str());
//    TLegend *legend = new TLegend(0.7,0.8,0.9,0.9);
//    
//    c1->SetLogy();
//    
//    h_signal->SetFillColor(kGreen);
//    h_signal->SetLineColor(kGreen);
//    h_signal->SetMarkerStyle(21);
//    h_signal->SetMarkerColor(kGreen);
//    
//    h_background->SetFillColor(kRed);
//    h_background->SetLineColor(kRed);
//    h_background->SetMarkerStyle(21);
//    h_background->SetMarkerColor(kRed);
//    
//    legend->AddEntry(h_signal, "Signal", "f");
//    legend->AddEntry(h_background, "Background", "f");
//    
//    hs->Add(h_background);
//    hs->Add(h_signal);
//    hs->Draw();
//    hs->GetXaxis()->SetTitle(plotName.c_str());
//    hs->GetYaxis()->SetTitle("N(Events)");
//    
//    legend->Draw();
//    
//    c1->Print(Form("%s%s",plotDir.c_str(),fileName.c_str()),"png");
//    
//    delete hs;
//    delete legend;
//    delete c1;
//    
//    // Plot Signal Ratio
//    plotSignalRatio(h_signalRatio,h_backgroundRatio, fileName, plotDir);
//    
//    delete h_signalRatio;
//    delete h_backgroundRatio;
//    
//}
//
//void CCProtonPi0_Plotter::plotStacked(TH1D* h_signal, TH1D* h_background, 
//                            std::string plotName, std::string fileName, std::string plotDir, 
//                            std::string signal_label, std::string background_label,
//                            bool isRatioReversed)
//{
//  
//    TH1D* h_signalRatio = new TH1D;
//    TH1D* h_backgroundRatio = new TH1D;
// 
//    h_signal->Copy(*h_signalRatio);
//    h_background->Copy(*h_backgroundRatio);
//
//    // Create Plot Objects
//    TCanvas *c1 = new TCanvas();
//    THStack *hs = new THStack("hs",plotName.c_str());
//    TLegend *legend = new TLegend(0.7,0.8,0.9,0.9);
//    
//    h_signal->SetFillColor(kGreen);
//    h_signal->SetLineColor(kGreen);
//    h_signal->SetMarkerStyle(21);
//    h_signal->SetMarkerColor(kGreen);
//    
//    h_background->SetFillColor(kRed);
//    h_background->SetLineColor(kRed);
//    h_background->SetMarkerStyle(21);
//    h_background->SetMarkerColor(kRed);
//    
//    legend->AddEntry(h_signal, signal_label.c_str(), "f");
//    legend->AddEntry(h_background, background_label.c_str(), "f");
// 
//    
//    hs->Add(h_background);
//    hs->Add(h_signal);
//    hs->Draw();
//    hs->GetXaxis()->SetTitle(plotName.c_str());
//    hs->GetYaxis()->SetTitle("N(Events)");
//    
//    legend->Draw();
//  
//    c1->Print(Form("%s%s",plotDir.c_str(),fileName.c_str()),"png");
//    
//    delete hs;
//    delete legend;
//    delete c1;
//    
//    // Plot Signal Ratio
//    plotSignalRatio(h_signalRatio,h_backgroundRatio,fileName, plotDir, isRatioReversed);
//    
//    // Plot Statistics
//    plot_purity_efficiency(h_signal,h_background,fileName,plotDir,true);
//    plot_purity_efficiency(h_signal,h_background,fileName,plotDir,false);
//
//    delete h_signalRatio;
//    delete h_backgroundRatio;
//}
//
//



#endif


