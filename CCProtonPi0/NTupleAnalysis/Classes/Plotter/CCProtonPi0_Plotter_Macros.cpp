/*
    See CCProtonPi0_Plotter.h header for Class Information
*/

#ifndef CCProtonPi0_Plotter_Macros_cpp
#define CCProtonPi0_Plotter_Macros_cpp

#include "CCProtonPi0_Plotter.h"

using namespace PlotUtils;

void CCProtonPi0_Plotter::DrawDataStackedMC(rootDir &dir, std::string var_name, std::string plotDir, int nCutArrows, CutArrow cutArrow1, CutArrow cutArrow2)
{
   DrawDataStackedMC_BckgAll(dir,var_name,plotDir, nCutArrows, cutArrow1, cutArrow2);
   //DrawDataStackedMC_BckgWithPi0(dir,var_name,plotDir, nCutArrows, cutArrow1, cutArrow2);
   //DrawDataStackedMC_BckgType(dir,var_name,plotDir, nCutArrows, cutArrow1, cutArrow2);
}

void CCProtonPi0_Plotter::DrawDataStackedMC_BckgAll(rootDir &dir, std::string var_name, std::string plotDir, int nCutArrows, CutArrow cutArrow1, CutArrow cutArrow2)
{
    std::string rootDir_data = dir.data;
    std::string rootDir_mc = dir.mc;

    TFile* f_data = new TFile(rootDir_data.c_str());
    TFile* f_mc = new TFile(rootDir_mc.c_str());
  
    std::string var;
    
    // ------------------------------------------------------------------------
    // Get Data Histogram
    // ------------------------------------------------------------------------
    var = Form("%s_%d",var_name.c_str(),0);
    MnvH1D* data = (MnvH1D*)f_data->Get(var.c_str()); 

    // ------------------------------------------------------------------------
    // Fill TObjArray - For MC Histograms
    // ------------------------------------------------------------------------
    TObjArray* mc_hists = new TObjArray;
    MnvH1D* temp;
   
    // Get All Background
    var = Form("%s_%d",var_name.c_str(),2);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Background");
    mc_hists->Add(temp);
    
    // Get Signal
    var = Form("%s_%d",var_name.c_str(),1);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Signal");
    mc_hists->Add(temp);

    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas("c","c",1280,800);
    ApplyStyle(plotter);
    plotter->DrawDataStackedMC(data,mc_hists,POT_Ratio_data_mc ,"TR","Data",2,1);
   
    // Add Plot Labels
    const double y_pos = 0.88;
    const double y_diff = 0.033;
    plotter->AddPlotLabel("Playlist: minerva1",0.3,y_pos,y_diff,kBlue);
    plotter->AddPOTNormBox(data_POT,mc_POT,0.3,y_pos-y_diff);
    
    // If Cut Histogram - Add Cut Arrows
    if (nCutArrows == 1){
        AddCutArrow(plotter, cutArrow1);
    }else if (nCutArrows == 2){
        AddCutArrow(plotter,cutArrow1);
        AddCutArrow(plotter,cutArrow2);
    }

    if (var_name.compare("hCut_1Track_pi0invMass") == 0 || var_name.compare("hCut_2Track_pi0invMass") == 0 ){
        TLine pi0Mass;
        pi0Mass.SetLineWidth(2);
        pi0Mass.SetLineColor(kBlue);
        pi0Mass.DrawLine(134.98,0,134.98,300);
    }
    // Print Plot
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),"_bckg_all.png"), "png");

    delete c;
}

void CCProtonPi0_Plotter::DrawDataStackedMC_BckgType(rootDir &dir, std::string var_name, std::string plotDir, int nCutArrows, CutArrow cutArrow1, CutArrow cutArrow2)
{
    std::string rootDir_data = dir.data;
    std::string rootDir_mc = dir.mc;

    TFile* f_data = new TFile(rootDir_data.c_str());
    TFile* f_mc = new TFile(rootDir_mc.c_str());
 
    std::string var;
    
    // ------------------------------------------------------------------------
    // Get Data Histogram
    // ------------------------------------------------------------------------
    var = Form("%s_%d",var_name.c_str(),0);
    MnvH1D* data = (MnvH1D*)f_data->Get(var.c_str()); 

    // ------------------------------------------------------------------------
    // Fill TObjArray - For MC Histograms
    // ------------------------------------------------------------------------
    TObjArray* mc_hists = new TObjArray;
    MnvH1D* temp;
   
    // Get Signal
    var = Form("%s_%d",var_name.c_str(),1);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Signal");
    mc_hists->Add(temp);
    
    // Get Bckg: NC
    var = Form("%s_%d",var_name.c_str(),6);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: NC");
    mc_hists->Add(temp);
    
    // Get Bckg: AntiNeutrino
    var = Form("%s_%d",var_name.c_str(),7);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: AntiNeutrino");
    mc_hists->Add(temp);

    // Get Bckg: QELike
    var = Form("%s_%d",var_name.c_str(),8);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: QELike");
    mc_hists->Add(temp);
   
    // Get Bckg: SinglePion
    var = Form("%s_%d",var_name.c_str(),9);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: SinglePion");
    mc_hists->Add(temp);

    // Get Bckg: DoublePion
    var = Form("%s_%d",var_name.c_str(),10);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: DoublePion");
    mc_hists->Add(temp);

    // Get Bckg: MultiPion
    var = Form("%s_%d",var_name.c_str(),11);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: MultiPion");
    mc_hists->Add(temp);
 
    // Get Bckg: Other
    var = Form("%s_%d",var_name.c_str(),12);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: Other");
    mc_hists->Add(temp);
    
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas("c","c",1280,800);
    ApplyStyle(plotter);
    plotter->DrawDataStackedMC(data, mc_hists,POT_Ratio_data_mc ,"TR","Data",3,2);
 
    // Add Plot Labels
    const double y_pos = 0.88;
    const double y_diff = 0.033;
    plotter->AddPlotLabel("Playlist: minerva1",0.3,y_pos,y_diff,kBlue);
    plotter->AddPOTNormBox(data_POT,mc_POT,0.3,y_pos-y_diff);   
    
    // If Cut Histogram - Add Cut Arrows
    if (nCutArrows == 1){
        AddCutArrow(plotter, cutArrow1);
    }else if (nCutArrows == 2){
        AddCutArrow(plotter,cutArrow1);
        AddCutArrow(plotter,cutArrow2);
    }

    // Print Plot
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),"_bckg_Type.png"), "png");

    delete c;
}

void CCProtonPi0_Plotter::DrawDataStackedMC_BckgWithPi0(rootDir &dir, std::string var_name, std::string plotDir, int nCutArrows, CutArrow cutArrow1, CutArrow cutArrow2)
{
    std::string rootDir_data = dir.data;
    std::string rootDir_mc = dir.mc;

    TFile* f_data = new TFile(rootDir_data.c_str());
    TFile* f_mc = new TFile(rootDir_mc.c_str());
    
    std::string var;

    // ------------------------------------------------------------------------
    // Get Data Histogram
    // ------------------------------------------------------------------------
    var = Form("%s_%d",var_name.c_str(),0);
    MnvH1D* data = (MnvH1D*)f_data->Get(var.c_str()); 

    // ------------------------------------------------------------------------
    // Fill TObjArray - For MC Histograms
    // ------------------------------------------------------------------------
    TObjArray* mc_hists = new TObjArray;
    MnvH1D* temp;
   
    // Get Signal
    var = Form("%s_%d",var_name.c_str(),1);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Signal");
    mc_hists->Add(temp);
    
    // Get Bckg: NoPi0
    var = Form("%s_%d",var_name.c_str(),3);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: NoPi0");
    mc_hists->Add(temp);
    
    // Get Bckg: SinglePi0
    var = Form("%s_%d",var_name.c_str(),4);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: SinglePi0");
    mc_hists->Add(temp);

    // Get Bckg: MultiPi0
    var = Form("%s_%d",var_name.c_str(),5);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Bckg: MultiPi0");
    mc_hists->Add(temp);
    
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    TCanvas* c = new TCanvas("c","c",1280,800);
    MnvPlotter* plotter = new MnvPlotter();
    ApplyStyle(plotter);
    plotter->DrawDataStackedMC(data,mc_hists,POT_Ratio_data_mc ,"TR","Data",3,2);

    // If Cut Histogram - Add Cut Arrows
    if (nCutArrows == 1){
        AddCutArrow(plotter, cutArrow1);
    }else if (nCutArrows == 2){
        AddCutArrow(plotter,cutArrow1);
        AddCutArrow(plotter,cutArrow2);
    }
    
    // Add Plot Labels
    const double y_pos = 0.88;
    const double y_diff = 0.033;
    plotter->AddPlotLabel("Playlist: minerva1",0.3,y_pos,y_diff,kBlue);
    plotter->AddPOTNormBox(data_POT,mc_POT,0.3,y_pos-y_diff);
    
    // Print Plot
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),"_bckg_Pi0.png"), "png");

    delete c;
}

void CCProtonPi0_Plotter::ApplyStyle(MnvPlotter* plotter)
{
    plotter->hist_min_zero = false;
    plotter->mc_line_width = 0;
    plotter->axis_title_size_x = 0.04;
    plotter->axis_title_size_y = 0.04;
    plotter->axis_label_size = 0.03;
    plotter->legend_text_size = 0.02;
}


void CCProtonPi0_Plotter::DrawDataMCRatio(rootDir& dir, std::string var_name, std::string plotDir)
{
    std::string rootDir_mc = dir.mc;
    std::string rootDir_data = dir.data;
    
    TFile* f_mc = new TFile(rootDir_mc.c_str());
    TFile* f_data = new TFile(rootDir_data.c_str());
   
    std::string var = Form("%s_%d",var_name.c_str(),0);
    MnvH1D* mc = (MnvH1D*)f_mc->Get(var.c_str());
    MnvH1D* data = (MnvH1D*)f_data->Get(var.c_str()); 

    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas("c","c",1280,800);
 
    // Plot
    plotter->DrawDataMCRatio(data, mc, POT_Ratio_data_mc);
    
    // Add Plot Labels
    const double y_pos = 0.88;
    const double y_diff = 0.033;
    plotter->AddPlotLabel("Playlist: minerva1",0.3,y_pos,y_diff,kBlue);
    plotter->AddPOTNormBox(data_POT,mc_POT,0.3,y_pos-y_diff);

    // Print Plot
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),"_ratio.png"), "png");

    delete c;
}

void CCProtonPi0_Plotter::DrawMCWithErrorBand(rootDir& dir, std::string var_name, std::string plotDir)
{
    std::string rootDir_mc = dir.mc;
    
    TFile* f_mc = new TFile(rootDir_mc.c_str());
   
    TH1D* mc = (TH1D*)f_mc->Get(var_name.c_str());

    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas("c","c",1280,800);

    // Plot
    plotter->DrawMCWithErrorBand(mc);

    // Add Plot Labels
    const double y_pos = 0.88;
    const double y_diff = 0.033;
    plotter->AddPlotLabel("Playlist: minerva1",0.3,y_pos,y_diff,kBlue);
    
    // Print Plot
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".png"), "png");

    delete c;
}

void CCProtonPi0_Plotter::DrawDataMC(rootDir& dir, std::string var_name, std::string plotDir)
{
    std::string rootDir_mc = dir.mc;
    std::string rootDir_data = dir.data;
    
    TFile* f_mc = new TFile(rootDir_mc.c_str());
    TFile* f_data = new TFile(rootDir_data.c_str());
   
    std::string var = Form("%s_%d",var_name.c_str(),0);
    MnvH1D* mc = (MnvH1D*)f_mc->Get(var.c_str());
    MnvH1D* data = (MnvH1D*)f_data->Get(var.c_str()); 

    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas("c","c",1280,800);

    // Plot
    plotter->DrawDataMC(data, mc, POT_Ratio_data_mc, "TR", false);

    // Add Plot Labels
    const double y_pos = 0.88;
    const double y_diff = 0.033;
    plotter->AddPlotLabel("Playlist: minerva1",0.3,y_pos,y_diff,kBlue);
    plotter->AddPOTNormBox(data_POT,mc_POT,0.3,y_pos-y_diff);
    
    // Print Plot
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".png"), "png");

    delete c;
    
    // Now Plot the Ratio
    DrawDataMCRatio(dir, var_name, plotDir);
}
void CCProtonPi0_Plotter::Draw1DHist(rootDir& dir, std::string var_name, std::string plotDir, bool isLogScale)
{
    std::string root_dir = dir.mc;
   
    // Get Histogram
    TFile* f = new TFile(root_dir.c_str());
    TH1D* hist1D = (TH1D*)f->Get(var_name.c_str());

    // Create Canvas
    TCanvas* c = new TCanvas("c","c",1280,800);
    if(isLogScale) c->SetLogy();

    // Plot Options
    hist1D->SetLineColor(kRed);
    hist1D->SetLineWidth(3);
    hist1D->SetFillColor(kRed);
    hist1D->SetFillStyle(3010);
   
    hist1D->Draw();
    gPad->Update();
    gStyle->SetOptStat("nemr"); 
    
    int max_bin = hist1D->GetMaximumBin();
    double max_bin_value = hist1D->GetBinCenter(max_bin);
    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.03);
    text.DrawLatex(0.78,0.7,Form("%s%3.2f", "Peak at ",max_bin_value));
    
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".png"), "png");
    delete c;
    
}

void CCProtonPi0_Plotter::Draw2DHist(rootDir& dir, std::string var_name, std::string plotDir, double threshold)
{
    std::string root_dir = dir.mc;
   
    // Get Histogram
    TFile* f = new TFile(root_dir.c_str());
    TH2D* hist2D = (TH2D*)f->Get(var_name.c_str());

    // Canvas
    Double_t w = 800; 
    Double_t h = 800;
    TCanvas* c = new TCanvas("c","c",w,h);
    c->SetWindowSize(w,h);
   
    // Reset Bins below the threshold
    int nBinsX = hist2D->GetNbinsX();
    int nBinsY = hist2D->GetNbinsY();
    for (int xBin = 1; xBin <= nBinsX; xBin++ ){
        for (int yBin = 1; yBin <=nBinsY; yBin++){
            int nEvents = hist2D->GetBinContent(xBin,yBin);
            if (nEvents <= threshold){
                hist2D->SetBinContent(xBin,yBin,0);
            }
        }
    }

    // Pad
    TPad *p = new TPad("p","p",0.05,0.05,0.95,0.95);
    p->Draw();
    
    p->cd();
    hist2D->GetYaxis()->SetTitleOffset(1.8);
    hist2D->Draw("colz");
    gPad->Update();
   
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".png"), "png");
    
    delete p;
    delete c;
}

void CCProtonPi0_Plotter::AddCutArrow(MnvPlotter* plotter, CutArrow &cutArrow)
{
    double cut_location = cutArrow.cut_location;
    double ymin = cutArrow.ymin;
    double ymax = cutArrow.ymax;
    double arrow_length = cutArrow.arrow_length;
    std::string arrow_direction = cutArrow.arrow_direction;
    
    plotter->AddCutArrow(cut_location, ymin, ymax, arrow_length, arrow_direction); 
}

void CCProtonPi0_Plotter::DrawStackedMC_BckgAll(rootDir &dir, std::string var_name, std::string plotDir, int nCutArrows, CutArrow cutArrow1, CutArrow cutArrow2)
{
    std::string rootDir_mc = dir.mc;

    TFile* f_mc = new TFile(rootDir_mc.c_str());
  
    std::string var;
    
    // ------------------------------------------------------------------------
    // Fill TObjArray - For MC Histograms
    // ------------------------------------------------------------------------
    TObjArray* mc_hists = new TObjArray;
    MnvH1D* temp;
   
    // Get All Background
    var = Form("%s_%d",var_name.c_str(),2);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Background");
    double nBckg = temp->GetEntries(); 
    mc_hists->Add(temp);
    
    // Get Signal
    var = Form("%s_%d",var_name.c_str(),1);
    temp = (MnvH1D*)f_mc->Get(var.c_str());
    temp->SetTitle("Signal");
    double nSignal = temp->GetEntries(); 
    mc_hists->Add(temp);

    // ------------------------------------------------------------------------
    // Plot  - If you want A Log Plot, axis_minimum = 0.1
    // ------------------------------------------------------------------------
    TCanvas* c = new TCanvas(var_name.c_str(),var_name.c_str(),1280,800);
    //gPad->SetLogy(); 
    //gPad->Update(); 
    MnvPlotter* plotter = new MnvPlotter();
    ApplyStyle(plotter);
    //plotter->axis_minimum = 0.1;
    plotter->DrawStackedMC(mc_hists,1,"TR");
   
    // Add Plot Labels
    const double y_pos = 0.88;
    const double y_diff = 0.033;
    plotter->AddPlotLabel(Form("nSignal = %3.0f, nBckg = %3.0f",nSignal,nBckg),0.3,y_pos,y_diff,kBlue); 
    //plotter->AddPOTNormBox(data_POT,mc_POT,0.3,y_pos-y_diff);
    
    // If Cut Histogram - Add Cut Arrows
    if (nCutArrows == 1){
        AddCutArrow(plotter, cutArrow1);
    }else if (nCutArrows == 2){
        AddCutArrow(plotter,cutArrow1);
        AddCutArrow(plotter,cutArrow2);
    }

    // Print Plot
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),"_mc_bckg_all.png"), "png");

    delete c;
}


void CCProtonPi0_Plotter::DrawStackedMC_GammaEvis(rootDir &dir, int gammaID, std::string plotDir)
{
    std::string rootDir_mc = dir.mc;

    TFile* f_mc = new TFile(rootDir_mc.c_str());
  
    // Use All MC Events -- Indice = 0
    std::string trkr = Form("g%d_evis_trkr_0",gammaID);
    std::string scal = Form("g%d_evis_scal_0",gammaID);
    std::string ecal = Form("g%d_evis_ecal_0",gammaID);
    std::string hcal = Form("g%d_evis_hcal_0",gammaID);
    
    // ------------------------------------------------------------------------
    // Fill TObjArray - For MC Histograms
    // ------------------------------------------------------------------------
    TObjArray* mc_hists = new TObjArray;
    MnvH1D* temp;
   
    // Get Tracker Evis
    temp = (MnvH1D*)f_mc->Get(trkr.c_str());
    temp->SetTitle("Tracker");
    double nTrkr = temp->GetEntries(); 
    mc_hists->Add(temp);
    
    // Get SideECAL Evis
    temp = (MnvH1D*)f_mc->Get(scal.c_str());
    temp->SetTitle("SideECAL");
    mc_hists->Add(temp);

    // Get ECAL Evis
    temp = (MnvH1D*)f_mc->Get(ecal.c_str());
    temp->SetTitle("ECAL");
    mc_hists->Add(temp);

    // Get HCAL Evis
    temp = (MnvH1D*)f_mc->Get(hcal.c_str());
    temp->SetTitle("HCAL");
    mc_hists->Add(temp);
    
    // ------------------------------------------------------------------------
    // Plot 
    // ------------------------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas("c","c",1280,800);
    ApplyStyle(plotter);
    plotter->DrawStackedMC(mc_hists,1,"TR");
    

    // Add Plot Labels
    const double y_pos = 0.88;
    const double y_diff = 0.033;
    plotter->AddPlotLabel(Form("nEvents = %3.0f",nTrkr),0.3,y_pos,y_diff,kBlue); 
    plotter->AddHistoTitle(Form("Visible Energy for Gamma %d",gammaID));


    // Print Plot
    c->Print(Form("%s%s%d%s",plotDir.c_str(),"gamma_",gammaID,"_mc_evis_stacked.png"), "png");

    delete c;
}

void CCProtonPi0_Plotter::DrawStackedMC_GammaByPDG(rootDir &dir, std::string var_name, int gammaID, std::string plotDir)
{
    // ------------------------------------------------------------------------- 
    // Get Histograms from File 
    // ------------------------------------------------------------------------- 
    std::string rootDir_mc = dir.mc;

    TFile* f_mc = new TFile(rootDir_mc.c_str());
  
    std::string pi0 = Form("g%d_%s_pi0",gammaID,var_name.c_str());
    std::string pi = Form("g%d_%s_pi",gammaID,var_name.c_str());
    std::string proton = Form("g%d_%s_proton",gammaID,var_name.c_str());
    std::string neutron = Form("g%d_%s_neutron",gammaID,var_name.c_str());
    std::string muon = Form("g%d_%s_muon",gammaID,var_name.c_str());
   
    TH1D* h_pi0 = new TH1D;   
    TH1D* h_pi = new TH1D;   
    TH1D* h_proton = new TH1D;   
    TH1D* h_neutron = new TH1D;   
    TH1D* h_muon = new TH1D;   

    h_pi0 = (TH1D*)f_mc->Get(pi0.c_str());
    h_pi = (TH1D*)f_mc->Get(pi.c_str());
    h_proton = (TH1D*)f_mc->Get(proton.c_str());
    h_neutron = (TH1D*)f_mc->Get(neutron.c_str());
    h_muon = (TH1D*)f_mc->Get(muon.c_str());

    // Set Histogram Style
    h_pi0->SetFillColor(kGreen);
    h_pi0->SetLineColor(kGreen);
    h_pi0->SetMarkerStyle(21);
    h_pi0->SetMarkerColor(kGreen);
    
    h_pi->SetFillColor(kRed);
    h_pi->SetLineColor(kRed);
    h_pi->SetMarkerStyle(21);
    h_pi->SetMarkerColor(kRed);
   
    h_proton->SetFillColor(kBlue);
    h_proton->SetLineColor(kBlue);
    h_proton->SetMarkerStyle(21);
    h_proton->SetMarkerColor(kBlue);

    h_neutron->SetFillColor(kCyan);
    h_neutron->SetLineColor(kCyan);
    h_neutron->SetMarkerStyle(21);
    h_neutron->SetMarkerColor(kCyan);

    h_muon->SetFillColor(kMagenta);
    h_muon->SetLineColor(kMagenta);
    h_muon->SetMarkerStyle(21);
    h_muon->SetMarkerColor(kMagenta);

    // Create Canvas and Form THStack
    TCanvas* c = new TCanvas("c","c",1280,800);
    THStack* hs = new THStack("hs",Form("%s by Particle for Gamma %d",var_name.c_str(),gammaID));
    TLegend* legend = new TLegend(0.7,0.8,0.9,0.9);

    //c->SetLogy();
    // Legend
    legend->AddEntry(h_pi0, "#pi^{0}", "f");
    legend->AddEntry(h_pi, "#pi^{#pm}", "f");
    legend->AddEntry(h_proton, "proton", "f");
    legend->AddEntry(h_neutron, "neutron", "f");
    legend->AddEntry(h_muon, "#mu^{-}", "f");
    
    hs->Add(h_pi0);
    hs->Add(h_pi);
    hs->Add(h_proton);
    hs->Add(h_neutron);
    hs->Add(h_muon);
    hs->Draw();
    hs->GetXaxis()->SetTitle(Form("%s by Particle for Gamma %d",var_name.c_str(),gammaID));
    hs->GetYaxis()->SetTitle("N(Events)");
    
    legend->Draw();
    
    c->Print(Form("%s%s%d%s%s%s",plotDir.c_str(),"gamma_",gammaID,"_mc_",var_name.c_str(),"_by_particle.png"), "png");

    delete hs;
    delete c;
}

#endif


