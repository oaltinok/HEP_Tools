#include "CCProtonPi0_Plotter.h"

const int max_iter = 4;

using namespace PlotUtils;

void CCProtonPi0_Plotter::UnfoldingStudy_muon_cos_theta()
{
    std::cout<<"Unfolding Study for muon_cos_theta"<<std::endl;

    std::string plotDir = Folder_List::plotDir_OtherStudies;

    // Use Different MC Samples to do the study
    // Train -- MC for Response Matrix
    // Sample -- MC for reco and truth Values
    rootDir Train;
    rootDir Sample;
    Train.mc = Folder_List::rootDir_Muon_Train;
    Sample.mc = Folder_List::rootDir_Muon_Sample;

    TFile* f_Train = new TFile(Train.mc.c_str());
    TFile* f_Sample = new TFile(Sample.mc.c_str());

    // Init Histograms Based on Pi0 Momentum
    MnvH2D* response = new MnvH2D(*(MnvH2D*)f_Train->Get("muon_cos_theta_response")); 
    MnvH1D* mc_reco = new MnvH1D(*(MnvH1D*)f_Sample->Get("muon_cos_theta_mc_reco_signal")); 
    MnvH1D* mc_true = new MnvH1D(*(MnvH1D*)f_Sample->Get("muon_cos_theta_mc_truth_signal")); 

    std::vector<MnvH1D*> unfolded;
    std::vector<MnvH1D*> error;
    std::vector<MnvH1D*> diff;
    init_UnfoldingHistograms(unfolded, error, diff);

    // 0 Iteration -- Original
    unfolded[0] = mc_reco;
    diff[0] = CalcUnfoldingDiff(unfolded[0], mc_true);
    error[0] = CalcUnfoldingError(unfolded[0], mc_true);

    // Fill Histograms with Different N(Iterations)
    for(int i = 1; i <= max_iter; ++i){
        FillUnfoldingHistograms(unfolded[i], error[i], diff[i], response, mc_reco, mc_true, i);
    }

    StyleUnfoldingHistograms(unfolded);
    StyleUnfoldingHistograms(error);
    StyleUnfoldingHistograms(diff);

    PlotUnfolding_Unfolded(unfolded, mc_true, "cos_muon_theta");
    PlotUnfolding_StatErrors(unfolded, mc_true, "cos_muon_theta");
    PlotUnfolding_Error(error, "cos_muon_theta");
    PlotUnfolding_Diff(diff, "cos_muon_theta");

    delete response;
    delete mc_reco;
    delete mc_true;
    delete f_Train;
    delete f_Sample;

    std::cout<<"Unfolding Study Finished!\n"<<std::endl;
}
void CCProtonPi0_Plotter::UnfoldingStudy_muon_theta()
{
    std::cout<<"Unfolding Study for muon_theta"<<std::endl;

    std::string plotDir = Folder_List::plotDir_OtherStudies;

    // Use Different MC Samples to do the study
    // Train -- MC for Response Matrix
    // Sample -- MC for reco and truth Values
    rootDir Train;
    rootDir Sample;
    Train.mc = Folder_List::rootDir_Muon_Train;
    Sample.mc = Folder_List::rootDir_Muon_Sample;

    TFile* f_Train = new TFile(Train.mc.c_str());
    TFile* f_Sample = new TFile(Sample.mc.c_str());

    // Init Histograms Based on Pi0 Momentum
    MnvH2D* response = new MnvH2D(*(MnvH2D*)f_Train->Get("muon_theta_response")); 
    MnvH1D* mc_reco = new MnvH1D(*(MnvH1D*)f_Sample->Get("muon_theta_mc_reco_signal")); 
    MnvH1D* mc_true = new MnvH1D(*(MnvH1D*)f_Sample->Get("muon_theta_mc_truth_signal")); 

    std::vector<MnvH1D*> unfolded;
    std::vector<MnvH1D*> error;
    std::vector<MnvH1D*> diff;
    init_UnfoldingHistograms(unfolded, error, diff);

    // 0 Iteration -- Original
    unfolded[0] = mc_reco;
    diff[0] = CalcUnfoldingDiff(unfolded[0], mc_true);
    error[0] = CalcUnfoldingError(unfolded[0], mc_true);

    // Fill Histograms with Different N(Iterations)
    for(int i = 1; i <= max_iter; ++i){
        FillUnfoldingHistograms(unfolded[i], error[i], diff[i], response, mc_reco, mc_true, i);
    }

    StyleUnfoldingHistograms(unfolded);
    StyleUnfoldingHistograms(error);
    StyleUnfoldingHistograms(diff);

    PlotUnfolding_Unfolded(unfolded, mc_true, "muon_theta");
    PlotUnfolding_StatErrors(unfolded, mc_true, "muon_theta");
    PlotUnfolding_Error(error, "muon_theta");
    PlotUnfolding_Diff(diff, "muon_theta");

    delete response;
    delete mc_reco;
    delete mc_true;
    delete f_Train;
    delete f_Sample;

    std::cout<<"Unfolding Study Finished!\n"<<std::endl;
}

void CCProtonPi0_Plotter::UnfoldingStudy_muon_P()
{
    std::cout<<"Unfolding Study for muon_P"<<std::endl;

    std::string plotDir = Folder_List::plotDir_OtherStudies;

    // Use Different MC Samples to do the study
    // Train -- MC for Response Matrix
    // Sample -- MC for reco and truth Values
    rootDir Train;
    rootDir Sample;
    Train.mc = Folder_List::rootDir_Muon_Train;
    Sample.mc = Folder_List::rootDir_Muon_Sample;

    TFile* f_Train = new TFile(Train.mc.c_str());
    TFile* f_Sample = new TFile(Sample.mc.c_str());

    // Init Histograms Based on Pi0 Momentum
    MnvH2D* response = new MnvH2D(*(MnvH2D*)f_Train->Get("muon_P_response")); 
    MnvH1D* mc_reco = new MnvH1D(*(MnvH1D*)f_Sample->Get("muon_P_mc_reco_signal")); 
    MnvH1D* mc_true = new MnvH1D(*(MnvH1D*)f_Sample->Get("muon_P_mc_truth_signal")); 

    std::vector<MnvH1D*> unfolded;
    std::vector<MnvH1D*> error;
    std::vector<MnvH1D*> diff;
    init_UnfoldingHistograms(unfolded, error, diff);

    // 0 Iteration -- Original
    unfolded[0] = new MnvH1D(*mc_reco);
    diff[0] = CalcUnfoldingDiff(unfolded[0], mc_true);
    error[0] = CalcUnfoldingError(unfolded[0], mc_true);

    // Fill Histograms with Different N(Iterations)
    for(int i = 1; i <= max_iter; ++i){
        FillUnfoldingHistograms(unfolded[i], error[i], diff[i], response, mc_reco, mc_true, i);
    }

    StyleUnfoldingHistograms(unfolded);
    StyleUnfoldingHistograms(error);
    StyleUnfoldingHistograms(diff);

    PlotUnfolding_Unfolded(unfolded, mc_true, "muon_P");
    PlotUnfolding_StatErrors(unfolded, mc_true, "muon_P");
    PlotUnfolding_Error(error, "muon_P");
    PlotUnfolding_Diff(diff, "muon_P");

    delete response;
    delete mc_reco;
    delete mc_true;
    delete f_Train;
    delete f_Sample;

    std::cout<<"Unfolding Study Finished!\n"<<std::endl;
}

void CCProtonPi0_Plotter::UnfoldingStudy_pi0_KE()
{
    std::cout<<"Unfolding Study for pi0_KE"<<std::endl;

    std::string plotDir = Folder_List::plotDir_OtherStudies;

    // Use Different MC Samples to do the study
    // Train -- MC for Response Matrix
    // Sample -- MC for reco and truth Values
    rootDir Train;
    rootDir Sample;
    Train.mc = Folder_List::rootDir_Pion_Train;
    Sample.mc = Folder_List::rootDir_Pion_Sample;

    TFile* f_Train = new TFile(Train.mc.c_str());
    TFile* f_Sample = new TFile(Sample.mc.c_str());

    // Init Histograms Based on Pi0 Momentum
    MnvH2D* response = new MnvH2D(*(MnvH2D*)f_Train->Get("pi0_KE_response")); 
    MnvH1D* mc_reco = new MnvH1D(*(MnvH1D*)f_Sample->Get("pi0_KE_mc_reco_signal")); 
    MnvH1D* mc_true = new MnvH1D(*(MnvH1D*)f_Sample->Get("pi0_KE_mc_truth_signal")); 

    std::vector<MnvH1D*> unfolded;
    std::vector<MnvH1D*> error;
    std::vector<MnvH1D*> diff;
    init_UnfoldingHistograms(unfolded, error, diff);

    // 0 Iteration -- Original
    unfolded[0] = mc_reco;
    diff[0] = CalcUnfoldingDiff(unfolded[0], mc_true);
    error[0] = CalcUnfoldingError(unfolded[0], mc_true);

    // Fill Histograms with Different N(Iterations)
    for(int i = 1; i <= max_iter; ++i){
        FillUnfoldingHistograms(unfolded[i], error[i], diff[i], response, mc_reco, mc_true, i);
    }

    StyleUnfoldingHistograms(unfolded);
    StyleUnfoldingHistograms(error);
    StyleUnfoldingHistograms(diff);

    PlotUnfolding_Unfolded(unfolded, mc_true, "pi0_KE");
    PlotUnfolding_StatErrors(unfolded, mc_true, "pi0_KE");
    PlotUnfolding_Error(error, "pi0_KE");
    PlotUnfolding_Diff(diff, "pi0_KE");

    delete response;
    delete mc_reco;
    delete mc_true;
    delete f_Train;
    delete f_Sample;

    std::cout<<"Unfolding Study Finished!\n"<<std::endl;
}

void CCProtonPi0_Plotter::UnfoldingStudy_pi0_P()
{
    std::cout<<"Unfolding Study for pi0_P"<<std::endl;

    std::string plotDir = Folder_List::plotDir_OtherStudies;

    // Use Different MC Samples to do the study
    // Train -- MC for Response Matrix
    // Sample -- MC for reco and truth Values
    rootDir Train;
    rootDir Sample;
    Train.mc = Folder_List::rootDir_Pion_Train;
    Sample.mc = Folder_List::rootDir_Pion_Sample;

    TFile* f_Train = new TFile(Train.mc.c_str());
    TFile* f_Sample = new TFile(Sample.mc.c_str());

    // Init Histograms Based on Pi0 Momentum
    MnvH2D* response = new MnvH2D(*(MnvH2D*)f_Train->Get("pi0_P_response")); 
    MnvH1D* mc_reco = new MnvH1D(*(MnvH1D*)f_Sample->Get("pi0_P_mc_reco_signal")); 
    MnvH1D* mc_true = new MnvH1D(*(MnvH1D*)f_Sample->Get("pi0_P_mc_truth_signal")); 

    std::vector<MnvH1D*> unfolded;
    std::vector<MnvH1D*> error;
    std::vector<MnvH1D*> diff;
    init_UnfoldingHistograms(unfolded, error, diff);

    // 0 Iteration -- Original
    unfolded[0] = mc_reco;
    diff[0] = CalcUnfoldingDiff(unfolded[0], mc_true);
    error[0] = CalcUnfoldingError(unfolded[0], mc_true);

    // Fill Histograms with Different N(Iterations)
    for(int i = 1; i <= max_iter; ++i){
        FillUnfoldingHistograms(unfolded[i], error[i], diff[i], response, mc_reco, mc_true, i);
    }

    StyleUnfoldingHistograms(unfolded);
    StyleUnfoldingHistograms(error);
    StyleUnfoldingHistograms(diff);

    PlotUnfolding_Unfolded(unfolded, mc_true, "pi0_P");
    PlotUnfolding_StatErrors(unfolded, mc_true, "pi0_P");
    PlotUnfolding_Error(error, "pi0_P");
    PlotUnfolding_Diff(diff, "pi0_P");

    delete response;
    delete mc_reco;
    delete mc_true;
    delete f_Train;
    delete f_Sample;

    std::cout<<"Unfolding Study Finished!\n"<<std::endl;
}

void CCProtonPi0_Plotter::UnfoldingStudy_pi0_theta()
{
    std::cout<<"Unfolding Study for pi0_theta"<<std::endl;

    std::string plotDir = Folder_List::plotDir_OtherStudies;

    // Use Different MC Samples to do the study
    // Train -- MC for Response Matrix
    // Sample -- MC for reco and truth Values
    rootDir Train;
    rootDir Sample;
    Train.mc = Folder_List::rootDir_Pion_Train;
    Sample.mc = Folder_List::rootDir_Pion_Sample;

    TFile* f_Train = new TFile(Train.mc.c_str());
    TFile* f_Sample = new TFile(Sample.mc.c_str());

    // Init Histograms Based on Pi0 Momentum
    MnvH2D* response = new MnvH2D(*(MnvH2D*)f_Train->Get("pi0_theta_response")); 
    MnvH1D* mc_reco = new MnvH1D(*(MnvH1D*)f_Sample->Get("pi0_theta_mc_reco_signal")); 
    MnvH1D* mc_true = new MnvH1D(*(MnvH1D*)f_Sample->Get("pi0_theta_mc_truth_signal")); 

    std::vector<MnvH1D*> unfolded;
    std::vector<MnvH1D*> error;
    std::vector<MnvH1D*> diff;
    init_UnfoldingHistograms(unfolded, error, diff);

    // 0 Iteration -- Original
    unfolded[0] = mc_reco;
    diff[0] = CalcUnfoldingDiff(unfolded[0], mc_true);
    error[0] = CalcUnfoldingError(unfolded[0], mc_true);

    // Fill Histograms with Different N(Iterations)
    for(int i = 1; i <= max_iter; ++i){
        FillUnfoldingHistograms(unfolded[i], error[i], diff[i], response, mc_reco, mc_true, i);
    }

    StyleUnfoldingHistograms(unfolded);
    StyleUnfoldingHistograms(error);
    StyleUnfoldingHistograms(diff);

    PlotUnfolding_Unfolded(unfolded, mc_true, "pi0_theta");
    PlotUnfolding_StatErrors(unfolded, mc_true, "pi0_theta");
    PlotUnfolding_Error(error, "pi0_theta");
    PlotUnfolding_Diff(diff, "pi0_theta");

    delete response;
    delete mc_reco;
    delete mc_true;
    delete f_Train;
    delete f_Sample;

    std::cout<<"Unfolding Study Finished!\n"<<std::endl;
}

void CCProtonPi0_Plotter::init_UnfoldingHistograms(std::vector<MnvH1D*> &unfolded, std::vector<MnvH1D*> &error, std::vector<MnvH1D*> &diff)
{
    // Fill with NULL Pointers
    MnvH1D* temp = NULL;

    for(int i = 0; i <= max_iter; ++i){
        unfolded.push_back(temp);
        error.push_back(temp);
        diff.push_back(temp);
    }
}

MnvH1D* CCProtonPi0_Plotter::CalcUnfoldingDiff(MnvH1D* unfolded, MnvH1D* truth)
{
    MnvH1D* diff = new MnvH1D(*unfolded);
    diff->Add(truth,-1); // Subtract Truth
    diff->GetYaxis()->SetTitle("abs(Reco-True)");

    // Take Absolute Value
    int nBins = diff->GetNbinsX();
    for (int i = 1; i <= nBins; ++i){
        double temp = diff->GetBinContent(i);
        if (temp < 0) diff->SetBinContent(i,-1*temp);
    }

    return diff;
}

MnvH1D* CCProtonPi0_Plotter::CalcUnfoldingError(MnvH1D* unfolded, MnvH1D* truth)
{
    // Calculate Reco-Truth
    MnvH1D* temp_diff = new MnvH1D(*unfolded);  // Start with Unfolded
    temp_diff->Add(truth, -1);                  // Subtract Truth

    MnvH1D* error = new MnvH1D(*temp_diff);
    error->Reset(); // Clear Bins
    error->Divide(temp_diff,truth);                  // Divide Difference with Truth
    error->GetYaxis()->SetTitle("(Reco-True)/True");

    delete temp_diff;
    return error;
}

void CCProtonPi0_Plotter::FillUnfoldingHistograms(MnvH1D* &unfolded, MnvH1D* &error, MnvH1D* &diff, MnvH2D* response, MnvH1D* mc_reco, MnvH1D* mc_true, int niter)
{
    unfolded = NULL;
    MinervaUnfold::MnvUnfold::Get().UnfoldHisto(unfolded, response, mc_reco, RooUnfold::kBayes, niter, true);
    diff = CalcUnfoldingDiff(unfolded, mc_true);
    error = CalcUnfoldingError(unfolded, mc_true);
}

void CCProtonPi0_Plotter::PlotUnfolding_StatErrors(const std::vector<MnvH1D*> &hists, const MnvH1D* truth, std::string var_name)
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    // Labels for THStack
    std::string plot_title = hists[1]->GetTitle();
    std::string X_title = hists[1]->GetXaxis()->GetTitle();
    std::string Y_title = hists[1]->GetYaxis()->GetTitle();

    TCanvas* c1 = new TCanvas("c","c",1280,800);
    THStack *hs = new THStack("hs",plot_title.c_str());
    TLegend *legend = new TLegend(0.7,0.68,0.9,0.9);  

    // Get Stat Errors for Truth
    TH1D* errTruth = new TH1D(truth->GetStatError(true));
    errTruth->SetFillStyle(0);
    errTruth->SetLineWidth(3);
    errTruth->SetLineColor(1);
    hs->Add(errTruth);
    legend->AddEntry(errTruth, "Truth", "L");

    // Get Stat Errors for All Iterations
    std::vector<TH1D*> errIter;
    for ( int i = 0; i <= max_iter; ++i){
        TH1D* temp = new TH1D(hists[i]->GetStatError(true));
        errIter.push_back(temp);
    }

    // Add Marker to No Unfolding
    errIter[0]->SetMarkerColor(kRed);
    legend->AddEntry(errIter[0],"No Unfolding","L");
    hs->Add(errIter[0]);

    for(int i = 1; i <= max_iter; ++i){
        errIter[i]->SetLineWidth(2);
        errIter[i]->SetLineStyle(7);
        legend->AddEntry(errIter[i], Form("%s%d","N(Iterations) = ",i),"L");
        hs->Add(errIter[i]);
    }

    hs->Draw("nostack");
    legend->Draw();

    // Add Axis Labels AFTER drawing hs
    hs->GetXaxis()->SetTitle(X_title.c_str());
    hs->GetYaxis()->SetTitle("Statistical Fractional Error");

    var_name = var_name + "_StatErrors.png";
    c1->Print(Form("%s%s",plotDir.c_str(),var_name.c_str()),"png");

    delete legend;
    delete hs;
    delete c1;
}

void CCProtonPi0_Plotter::PlotUnfolding_Unfolded(const std::vector<MnvH1D*> &hists, const MnvH1D* truth, std::string var_name)
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    // Labels for THStack
    std::string plot_title = hists[1]->GetTitle();
    std::string X_title = hists[1]->GetXaxis()->GetTitle();
    std::string Y_title = hists[1]->GetYaxis()->GetTitle();

    TCanvas* c1 = new TCanvas("c","c",1280,800);
    THStack *hs = new THStack("hs",plot_title.c_str());
    TLegend *legend = new TLegend(0.7,0.68,0.9,0.9);  

    // Create Temp Objects, we do not change the actual objects
    MnvH1D* tempTruth = new MnvH1D(*truth);
    std::vector<MnvH1D*> tempHists;
    for (int i = 0; i <= max_iter; ++i){
        MnvH1D* temp = new MnvH1D(*hists[i]);
        // Normalize to Bin Width
        double norm_bin_width = temp->GetNormBinWidth();
        temp->Scale(norm_bin_width, "width");
        tempHists.push_back(temp);
    }

    // Add Truth Value to Unfolded Only 
    tempTruth->SetMarkerStyle(21);
    tempTruth->SetMarkerSize(2);
    tempTruth->SetFillStyle(0);
    tempTruth->SetLineWidth(3);
    tempTruth->SetLineColor(1);
    double norm_bin_width = tempTruth->GetNormBinWidth();
    tempTruth->Scale(norm_bin_width, "width");
    hs->Add(tempTruth);
    legend->AddEntry(tempTruth, "Truth", "LP");

    // Add Marker to No Unfolding
    tempHists[0]->SetMarkerStyle(20);
    tempHists[0]->SetMarkerSize(2);
    tempHists[0]->SetMarkerColor(kRed);
    legend->AddEntry(tempHists[0],"No Unfolding","LP");
    hs->Add(tempHists[0]);

    for(int i = 1; i <= max_iter; ++i){
        legend->AddEntry(tempHists[i], Form("%s%d","N(Iterations) = ",i),"L");
        hs->Add(tempHists[i]);
    }

    hs->Draw("nostack");
    legend->Draw();

    // Add Axis Labels AFTER drawing hs
    hs->GetXaxis()->SetTitle(X_title.c_str());
    hs->GetYaxis()->SetTitle(Y_title.c_str());

    var_name = var_name + "_Unfolded.png";
    c1->Print(Form("%s%s",plotDir.c_str(),var_name.c_str()),"png");


    // Root does not like deleting histograms
    delete legend;
    delete hs;
    delete c1;
}

void CCProtonPi0_Plotter::PlotUnfolding_Error(const std::vector<MnvH1D*> &hists, std::string var_name)
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    // Labels for THStack
    std::string plot_title = hists[1]->GetTitle();
    std::string X_title = hists[1]->GetXaxis()->GetTitle();
    std::string Y_title = hists[1]->GetYaxis()->GetTitle();

    TCanvas* c1 = new TCanvas("c","c",1280,800);
    TLegend *legend = new TLegend(0.7,0.68,0.9,0.9);  
    THStack *hs = new THStack("hs",plot_title.c_str());

    // Create Temp Objects, we do not change the actual objects
    std::vector<MnvH1D*> tempHists;
    for (int i = 0; i <= max_iter; ++i){
        MnvH1D* temp = new MnvH1D(*hists[i]);
        // Normalize to Bin Width
        double norm_bin_width = temp->GetNormBinWidth();
        temp->Scale(norm_bin_width, "width");
        tempHists.push_back(temp);
    }

    legend->AddEntry(tempHists[0],"No Unfolding","L");
    hs->Add(tempHists[0]);

    for(int i = 1; i <= max_iter; ++i){
        legend->AddEntry(tempHists[i], Form("%s%d","N(Iterations) = ",i),"L");
        hs->Add(tempHists[i]);
    }

    hs->SetMinimum(-0.7);
    hs->SetMaximum(0.7);
    hs->Draw("nostack HIST");
    legend->Draw();

    // Add Axis Labels AFTER drawing hs
    hs->GetXaxis()->SetTitle(X_title.c_str());
    hs->GetYaxis()->SetTitle(Y_title.c_str());

    var_name = var_name + "_Error.png";
    c1->Print(Form("%s%s",plotDir.c_str(),var_name.c_str()),"png");

    // Root does not like deleting histograms
    delete legend;
    delete hs;
    delete c1;
}

void CCProtonPi0_Plotter::PlotUnfolding_Diff(const std::vector<MnvH1D*> &hists, std::string var_name)
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    // Labels for THStack
    std::string plot_title = hists[1]->GetTitle();
    std::string X_title = hists[1]->GetXaxis()->GetTitle();
    std::string Y_title = hists[1]->GetYaxis()->GetTitle();

    TCanvas* c1 = new TCanvas("c","c",1280,800);
    TLegend *legend = new TLegend(0.7,0.68,0.9,0.9);  
    THStack *hs = new THStack("hs",plot_title.c_str());

    // Create Temp Objects, we do not change the actual objects
    std::vector<MnvH1D*> tempHists;
    for (int i = 0; i <= max_iter; ++i){
        MnvH1D* temp = new MnvH1D(*hists[i]);
        // Normalize to Bin Width
        double norm_bin_width = temp->GetNormBinWidth();
        temp->Scale(norm_bin_width, "width");
        tempHists.push_back(temp);
    }

    legend->AddEntry(tempHists[0],"No Unfolding","L");
    hs->Add(tempHists[0]);

    for(int i = 1; i <= max_iter; ++i){
        legend->AddEntry(tempHists[i], Form("%s%d","N(Iterations) = ",i),"L");
        hs->Add(tempHists[i]);
    }

    hs->Draw("nostack HIST");
    legend->Draw();

    // Add Axis Labels AFTER drawing hs
    hs->GetXaxis()->SetTitle(X_title.c_str());
    hs->GetYaxis()->SetTitle(Y_title.c_str());

    // Add Integrals as Text 
    TLatex* text = new TLatex;
    text->SetTextSize(0.03);
    text->SetNDC();
    text->DrawLatex(0.7, 0.64, "Area of the Curves");
    text->DrawLatex(0.7, 0.60, Form("%s%3.2f", "Area(0 Iter) = ", tempHists[0]->Integral()));
    text->DrawLatex(0.7, 0.56, Form("%s%3.2f", "Area(1 Iter) = ", tempHists[1]->Integral()));
    text->DrawLatex(0.7, 0.52, Form("%s%3.2f", "Area(2 Iter) = ", tempHists[2]->Integral()));
    text->DrawLatex(0.7, 0.48, Form("%s%3.2f", "Area(3 Iter) = ", tempHists[3]->Integral()));
    text->DrawLatex(0.7, 0.44, Form("%s%3.2f", "Area(4 Iter) = ", tempHists[4]->Integral()));
    delete text;

    var_name = var_name + "_Diff.png";
    c1->Print(Form("%s%s",plotDir.c_str(),var_name.c_str()),"png");

    
    // Root does not like deleting histograms
    delete legend;
    delete hs;
    delete c1;
}

void CCProtonPi0_Plotter::PlotUnfolding_Migration()
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;
    // Sample is Data
    // Train is MC
    rootDir pi0;
    rootDir muon;
    pi0.mc = Folder_List::rootDir_Pion_Train;
    muon.mc = Folder_List::rootDir_Muon_Train;

    // Plot Migration Histograms
    DrawNormalizedMigrationHistogram(muon, "muon_P_response", plotDir);
    DrawNormalizedMigrationHistogram(muon, "muon_theta_response", plotDir);
    DrawNormalizedMigrationHistogram(muon, "muon_cos_theta_response", plotDir);
    DrawNormalizedMigrationHistogram(pi0, "pi0_P_response", plotDir);
    DrawNormalizedMigrationHistogram(pi0, "pi0_KE_response", plotDir);
    DrawNormalizedMigrationHistogram(pi0, "pi0_theta_response", plotDir);
}

void CCProtonPi0_Plotter::PlotUnfolding_TruthComparison()
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    // Sample is Data
    // Train is MC
    rootDir pi0;
    rootDir muon;
    pi0.data = Folder_List::rootDir_Pion_Sample;
    pi0.mc = Folder_List::rootDir_Pion_Train;
    muon.data = Folder_List::rootDir_Muon_Sample;
    muon.mc = Folder_List::rootDir_Muon_Train;

    TFile* f_muon_mc = new TFile(muon.mc.c_str());
    TFile* f_muon_data = new TFile(muon.data.c_str());
    TFile* f_pi0_mc = new TFile(pi0.mc.c_str());
    TFile* f_pi0_data = new TFile(pi0.data.c_str());

    MnvH1D* data = (MnvH1D*)f_muon_data->Get("muon_P_mc_truth_signal");
    MnvH1D* mc= (MnvH1D*)f_muon_mc->Get("muon_P_mc_truth_signal");
    DrawDataMC(data, mc, "muon_P_truth_comparison", plotDir);

    data = (MnvH1D*)f_muon_data->Get("muon_theta_mc_truth_signal");
    mc= (MnvH1D*)f_muon_mc->Get("muon_theta_mc_truth_signal");
    DrawDataMC(data, mc, "muon_theta_truth_comparison", plotDir);

    data = (MnvH1D*)f_muon_data->Get("muon_cos_theta_mc_truth_signal");
    mc= (MnvH1D*)f_muon_mc->Get("muon_cos_theta_mc_truth_signal");
    DrawDataMC(data, mc, "muon_cos_theta_truth_comparison", plotDir);

    data = (MnvH1D*)f_pi0_data->Get("pi0_P_mc_truth_signal");
    mc= (MnvH1D*)f_pi0_mc->Get("pi0_P_mc_truth_signal");
    DrawDataMC(data, mc, "pi0_P_truth_comparison", plotDir);

    data = (MnvH1D*)f_pi0_data->Get("pi0_KE_mc_truth_signal");
    mc= (MnvH1D*)f_pi0_mc->Get("pi0_KE_mc_truth_signal");
    DrawDataMC(data, mc, "pi0_KE_truth_comparison", plotDir);

    data = (MnvH1D*)f_pi0_data->Get("pi0_theta_mc_truth_signal");
    mc= (MnvH1D*)f_pi0_mc->Get("pi0_theta_mc_truth_signal");
    DrawDataMC(data, mc, "pi0_theta_truth_comparison", plotDir);

    delete f_muon_mc;
    delete f_muon_data;
    delete f_pi0_mc;
    delete f_pi0_data;
}

void CCProtonPi0_Plotter::StyleUnfoldingHistograms(std::vector<MnvH1D*> &hists)
{
    // Common Style for All Histograms
    for (int i = 0; i <=max_iter; ++i){
        hists[i]->SetFillStyle(0);            
        hists[i]->SetLineWidth(2);
    }

    // Unique Style 
    hists[0]->SetLineWidth(3);          // No Iterations
    hists[0]->SetLineColor(kRed);       // No Iterations

    hists[1]->SetLineColor(kBlue); 
    hists[2]->SetLineColor(kCyan);  
    hists[3]->SetLineColor(kGreen+4);
    hists[4]->SetLineColor(kMagenta);
}

