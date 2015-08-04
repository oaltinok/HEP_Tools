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
   DrawDataStackedMC_BckgWithPi0(dir,var_name,plotDir, nCutArrows, cutArrow1, cutArrow2);
   DrawDataStackedMC_BckgType(dir,var_name,plotDir, nCutArrows, cutArrow1, cutArrow2);
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
    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas("c","c",1280,800);
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
    
    c->Print(Form("%s%s%s",plotDir.c_str(),var_name.c_str(),".png"), "png");
    
    delete c;
    
}

void CCProtonPi0_Plotter::Draw2DHist(rootDir& dir, std::string var_name, std::string plotDir)
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
    
    // Pad
    TPad *p = new TPad("p","p",0.05,0.05,0.95,0.95);
    p->Draw();
    
    p->cd();
    hist2D->GetYaxis()->SetTitleOffset(1.8);
    hist2D->Draw("colz");
    gPad->Update();
    
    // Statistics Box
    //TPaveStats *ps = (TPaveStats*)se->GetPrimitive("stats");    
    //st->SetOptStat(1000000110);
    //st->SetX1NDC(0.1); 
    //st->SetX2NDC(0.3); 
    //st->SetY1NDC(0.8); 
    //st->SetY2NDC(0.9); 
   
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



#endif


