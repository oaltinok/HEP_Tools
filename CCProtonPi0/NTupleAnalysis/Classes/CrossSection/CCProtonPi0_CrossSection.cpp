#ifndef CCProtonPi0_CrossSection_cpp
#define CCProtonPi0_CrossSection_cpp

#include "CCProtonPi0_CrossSection.h"

using namespace PlotUtils;

CCProtonPi0_CrossSection::CCProtonPi0_CrossSection(bool isMC) : CCProtonPi0_NTupleAnalysis()
{
    std::cout<<"Initializing CCProtonPi0_CrossSection"<<std::endl;

    iteration = 1; // Unfolding Iteration
    min_invMass = 60; // MeV
    max_invMass = 200; // MeV
    
    m_isMC = isMC;

    OpenRootFiles();
    initHistograms();
    initXSecs();
    std::cout<<"Done!"<<std::endl;
}

void CCProtonPi0_CrossSection::Calc_CrossSections()
{
    Calc_Normalized_NBackground();
    
    Calc_CrossSection(muon_P);   
    Calc_CrossSection(muon_theta);   
    Calc_CrossSection(pi0_P);   
    Calc_CrossSection(pi0_KE);   
    Calc_CrossSection(pi0_theta);   
    Calc_CrossSection(QSq);   
    
    writeHistograms();
}

void CCProtonPi0_CrossSection::Calc_CrossSection(XSec &var)
{
    std::cout<<"\n-----------------------------------------------------------------------"<<std::endl;
    std::cout<<"Calculating Cross Section for "<<var.name<<std::endl;

    // Data Correction
    var.bckg_subtracted = Subtract_Background(var.all, var.mc_reco_bckg, var.bckg_estimated,var.name);
    var.unfolded = Unfold_Data(var.bckg_subtracted, var.response, var.name);   
    var.efficiency_corrected = Efficiency_Divide(var.unfolded, var.eff, var.name);   

    // Integrate Flux
    var.integrated_flux = Integrate_Flux(var.efficiency_corrected, var.name, var.isEv);
    
    // Calculate Final Cross Section
    var.xsec= Calc_FinalCrossSection(var.efficiency_corrected, var.integrated_flux, var.name); 
 
    Style_XSec(var);
    std::cout<<"Done!"<<std::endl;
    std::cout<<"-----------------------------------------------------------------------"<<std::endl;
}

void CCProtonPi0_CrossSection::Style_XSec(XSec &var)
{
    // Add Labels
    var.xsec->SetTitle(var.plot_title.c_str());
    var.xsec->GetXaxis()->SetTitle(var.plot_xlabel.c_str());
    var.xsec->GetYaxis()->SetTitle(var.plot_ylabel.c_str());
    
    // Scale Cross Section Result to match with Label
    var.xsec->Scale(var.scale);
}

void CCProtonPi0_CrossSection::Calc_Normalized_NBackground()
{
    std::cout<<"\nCalculating N(Background) in Data"<<std::endl;
   
    // ------------------------------------------------------------------------
    // Get Histograms
    // ------------------------------------------------------------------------
  
    // Get TH1D from MnvH1D
    TH1D* data_invMass = new TH1D(invMass_all->GetCVHistoWithStatError());
    TH1D* signal_invMass = new TH1D(invMass_mc_reco_signal->GetCVHistoWithStatError());
    TH1D* bckg_invMass = new TH1D(invMass_mc_reco_bckg->GetCVHistoWithStatError());
    
    invMass_fit_signal = new TH1D(*data_invMass); // Start as data hist (will subtract bckg to get correct signal data)
    invMass_fit_signal->SetName("invMass_fit_signal");
    
    invMass_fit_bckg = new TH1D(*bckg_invMass);
    invMass_fit_bckg->SetName("invMass_bckg");

    // ------------------------------------------------------------------------
    // Estimate Background using TFractionFitter
    // ------------------------------------------------------------------------
    TObjArray* mc_models = new TObjArray(2);
    mc_models->Add(signal_invMass);
    mc_models->Add(bckg_invMass);

    TFractionFitter* fitter = new TFractionFitter(data_invMass, mc_models, "Q");
    fitter->Constrain(1, 0.0, 1.0);
    
    std::cout<<"\tFitting Background Shape..."<<std::endl;
    Int_t status = fitter->Fit();

    if (status != 0) {
        std::cerr<<"\t\tFit Error!"<<std::endl;
    }
    
    invMass_fit_result  = new TH1D (*(TH1D*)fitter->GetPlot());
    invMass_fit_result->SetName("data_fit_result"); 

    Double_t N_Data = data_invMass->Integral();
    Double_t f0 =     0.0; // fittted signal fraction
    Double_t f1 =     0.0; // fitted background fraction
    Double_t f0_err = 0.0; // fitted signal fraction uncertainty
    Double_t f1_err = 0.0; // fitted background fraction uncertainty
    
    fitter->GetResult(0,f0,f0_err); // access fitted value and uncertainty for parameter 1
    fitter->GetResult(1,f1,f1_err); // access fitted value and uncertainty for parameter 2

    double rf0 = 1.e2 * f0_err/f0;
    double rf1 = 1.e2 * f1_err/f1;
    Uncertainity_Bckg = rf1;
    std::cout<<"\tFit Results:"<<std::endl;
    std::cout<<"\t\tSignal Fraction = "<<f0<<" Error = "<<f0_err<<" Uncertainity = "<<rf0<<std::endl;
    std::cout<<"\t\tBackground Fraction = "<<f1<<" Error = "<<f1_err<<" Uncertainity = "<<rf1<<std::endl;

    if (status == 0) {
        double area = invMass_fit_bckg->Integral();
        double total_fitted_background = f1 * N_Data; // Total N(Bckg) in Data
        std::cout<<"\tAll Events in Data (Whole Spectrum) = "<<N_Data<<std::endl;
        std::cout<<"\tBackground in Data (Whole Spectrum) = "<<total_fitted_background<<std::endl;
        invMass_fit_bckg->Scale(total_fitted_background/area); 
        // Subtract scaled bckg hist from data hist
        invMass_fit_signal->Add(invMass_fit_bckg,-1); 

            // this calculation overestimates the number of background events
            // in the lower and upper bins
            // nbkg = fitted_background->Integral(__lower_bin, __upper_bin);

            // integrate the number of background events correctly
            // notice the mass range, not bin range
        N_Background_Data = Integrate_SignalRegion(invMass_fit_bckg);
        double N_Signal_Data = Integrate_SignalRegion(invMass_fit_signal);
        double N_All_Data = Integrate_SignalRegion(data_invMass);

        std::cout<<std::endl;
        std::cout<<"\tAll Events in Data (Signal Region) = "<<N_All_Data<<std::endl; 
        std::cout<<"\tBackground in Data (Signal Region) = "<<N_Background_Data<<std::endl; 
        std::cout<<"\tSignal in Data (Signal Region) = "<<N_Signal_Data<<std::endl; 
    }

    delete mc_models;
    delete fitter;
}

double CCProtonPi0_CrossSection::Integrate_SignalRegion(TH1D* h)
{
    int lower_bin = h->FindBin(min_invMass);
    int upper_bin = h->FindBin(max_invMass);
    
    double nbkg = h->Integral(lower_bin, upper_bin);
    
        // assume uniform event distribution within the bin,
        // calculate the over counting for the lower bin
    double low_edge = h->GetBinLowEdge(lower_bin);
        //std::cout << "lower bound " << lower_bound << std::endl;
        //std::cout << "low_edge    " << low_edge    << std::endl;
    double dm1 = min_invMass - low_edge;
    
    double bin_width1   = h->GetBinWidth(lower_bin);
    double bin_content1 = h->GetBinContent(lower_bin);
    double overcount1   = (dm1/bin_width1) * bin_content1;
    
        // assume uniform event distribution within the bin,
        // calculate the over counting for the upper bin
    double high_edge = h->GetBinLowEdge(upper_bin) + h->GetBinWidth(upper_bin);
        //std::cout << "upper bound " << upper_bound << std::endl;
        //std::cout << "high_edge    " << high_edge    << std::endl;
    double dm2 = high_edge - max_invMass;
    
    double bin_width2   = h->GetBinWidth(upper_bin);
    double bin_content2 = h->GetBinContent(upper_bin);
    double overcount2   = (dm2/bin_width2) * bin_content2;
    
        // correct for over-counting
    double corrected_nbkg = nbkg - overcount1 - overcount2;

        /*
    std::cout << " uncorrected integral, corrected integral: "
              << nbkg << "," << corrected_nbkg
              << std::endl;
        */
    
    return corrected_nbkg;
}

void CCProtonPi0_CrossSection::NormalizeHistogram(MnvH1D* h)
{
    std::cout<<"\tNormalizing Background Shape"<<std::endl;
    int NBins = h->GetNbinsX();
    double area = h->Integral();
    double nOverFlow = h->GetBinContent(NBins+1);
    std::cout<<"\t\tBefore Norm = "<<area<<std::endl;
    h->Scale(1/(area+nOverFlow));
    std::cout<<"\t\tAfter Norm = "<<h->Integral()<<std::endl;
    std::cout<<"\tDone!"<<std::endl;
}

MnvH1D* CCProtonPi0_CrossSection::Subtract_Background(MnvH1D* data, MnvH1D* mc_bckg, MnvH1D* &bckg_estimated, std::string var_name)
{
    std::cout<<"Subtracting Background for "<<var_name<<std::endl;
    // Init Histogram
    MnvH1D* bckg_subtracted = new MnvH1D(*data); 
    std::string hist_name = var_name + "_bckg_subtracted";
    bckg_subtracted->SetName(hist_name.c_str());
    bckg_subtracted->SetTitle("Background Subtracted Data");
   
    // Get Shape of Background
    //      Unit area -- Including Overflow
    //      After Normalization, Area MUST be less than 1
    NormalizeHistogram(mc_bckg);

    // Estimate N Background in Variable
    mc_bckg->Scale(N_Background_Data);

    // Subtracted Background
    bckg_subtracted->Add(mc_bckg,-1);
  
    // Estimated Background
    bckg_estimated = new MnvH1D(*mc_bckg);
    hist_name = var_name + "_bckg_estimated";
    bckg_estimated->SetName(hist_name.c_str());
    bckg_estimated->SetTitle("Estimated Background in Data"); 

    std::cout<<"\tTotal Data Area = "<<data->Integral()<<std::endl;
    std::cout<<"\tEstimated Background in Data Area = "<<bckg_estimated->Integral()<<std::endl;
    std::cout<<"\tBackground Subtracted Data Area = "<<bckg_subtracted->Integral()<<std::endl;

    std::cout<<"Done!"<<std::endl;

    return bckg_subtracted;
}

MnvH1D* CCProtonPi0_CrossSection::Unfold_Data(MnvH1D* bckg_subtracted, MnvH2D* response, std::string var_name)
{
    std::cout<<"Unfolding Data for "<<var_name<<std::endl;
    // Init Histogram
    MnvH1D* unfolded = 0;

    std::cout<<"\tNumber of iteratons = "<<iteration<<std::endl;
    // Use MnvUnfold to Unfold Data
    MinervaUnfold::MnvUnfold::Get().UnfoldHisto(unfolded, response, bckg_subtracted, RooUnfold::kBayes, iteration, true);

    // Set Name of the Histogram
    std::string hist_name = var_name + "_unfolded";
    unfolded->SetName(hist_name.c_str());
    unfolded->SetTitle("Unfolded Data");

    std::cout<<"Done!"<<std::endl;

    return unfolded;
}

MnvH1D* CCProtonPi0_CrossSection::Efficiency_Divide(MnvH1D* unfolded, MnvH1D* eff, std::string var_name)
{
    std::cout<<"Efficiency Correction for "<<var_name<<std::endl;
  
    // Init Histogram
    MnvH1D* efficiency_corrected = new MnvH1D(*unfolded);
    std::string hist_name = var_name + "_efficiency_corrected";
    efficiency_corrected->SetName(hist_name.c_str());
    efficiency_corrected->SetTitle("Efficiency Corrected Data");

    // Divide by Efficiency
    std::cout<<"\tArea Before = "<<efficiency_corrected->Integral()<<std::endl;
    efficiency_corrected->Divide(unfolded, eff);
    std::cout<<"\tArea After = "<<efficiency_corrected->Integral()<<std::endl;
    
    std::cout<<"Done!"<<std::endl;

    return efficiency_corrected;
}

MnvH1D* CCProtonPi0_CrossSection::Integrate_Flux(MnvH1D* data_efficiency_corrected, std::string var_name, bool isEv)
{
    std::cout<<"Integrating Flux for "<<var_name<<std::endl;
    
    MnvH1D* integrated_flux = calc_flux(data_efficiency_corrected, rootDir_flux, isEv);
    
    std::string hist_name = var_name + "_integrated_flux";
    integrated_flux->SetName(hist_name.c_str());

    std::cout<<"Done!"<<std::endl;

    return integrated_flux;
}

MnvH1D* CCProtonPi0_CrossSection::Calc_FinalCrossSection(MnvH1D* data_efficiency_corrected, MnvH1D* integrated_flux, std::string var_name)
{
    MnvH1D* h_xs = new MnvH1D(*data_efficiency_corrected);
    std::string hist_name = var_name + "_xsec";
    h_xs->SetName(hist_name.c_str());

    // Divide the efficiency corrected distribution by the integrated flux
    h_xs->Divide(data_efficiency_corrected, integrated_flux);

    // Set POT
    double pot;
    if (m_isMC) pot = mc_POT;
    else pot = data_POT;

    //fiducial is modules 27-79 inclusive
    int nplanes = 2 * ( 79 - 27 + 1 );
    double n_atoms    = TargetUtils::Get().GetTrackerNCarbonAtoms(nplanes, m_isMC, 850.0);
    double n_nucleons = TargetUtils::Get().GetTrackerNNucleons(nplanes, m_isMC, 850.0);

    std::cout << "Number of C atoms   (1e30): " << n_atoms/1e30 << std::endl;
    std::cout << "Number of neutrons: (1e30): " << n_nucleons/1e30 << std::endl;

    h_xs->Scale(1/pot);
    h_xs->Scale(1/n_nucleons);

    h_xs->Scale(1e40); // to quote xs in 1-e40, 1e-42 for theta;

    return h_xs;
}


void CCProtonPi0_CrossSection::initXSecs()
{
    // See CCProtonPi0_CrossSection_Style.cpp for Definitions 
    init_muon_P();
    init_muon_theta();
    init_pi0_P();
    init_pi0_KE();
    init_pi0_theta();
    init_QSq();
}

void CCProtonPi0_CrossSection::OpenRootFiles()
{
    // Flux File
    rootDir_flux = Folder_List::rootDir_Flux_new;

    // Output File
    if (m_isMC) rootDir_out = Folder_List::rootDir_CrossSection_mc;
    else rootDir_out = Folder_List::rootDir_CrossSection_data;
    std::cout<<"\tRoot File: "<<rootDir_out<<std::endl;
    f_out = new TFile(rootDir_out.c_str(), "RECREATE");

    // Truth File
    std::string rootDir = Folder_List::rootDir_Truth_mc;
    f_truth = new TFile(rootDir.c_str());

    // Cut Histograms -- For Invariant Mass 
    rootDir = Folder_List::rootDir_CutHists_data;
    f_data_cutHists = new TFile(rootDir.c_str());
 
    rootDir = Folder_List::rootDir_CutHists_mc;
    f_mc_cutHists = new TFile(rootDir.c_str());
}

void CCProtonPi0_CrossSection::initHistograms()
{
    // Pi0 Invariant Mass
    if (m_isMC) invMass_all = new MnvH1D(*(MnvH1D*)f_mc_cutHists->Get("invMass_mc_reco_all"));
    else invMass_all = new MnvH1D(*(MnvH1D*)f_data_cutHists->Get("invMass_all"));

    invMass_mc_reco_signal = new MnvH1D(*(MnvH1D*)f_mc_cutHists->Get("invMass_mc_reco_signal")); 
    invMass_mc_reco_bckg  = new MnvH1D(*(MnvH1D*)f_mc_cutHists->Get("invMass_mc_reco_bckg")); 
}

void CCProtonPi0_CrossSection::initHistograms(XSec &var)
{
    std::string hist_name;
    
    if (m_isMC){ 
        hist_name = var.name + "_mc_reco_all";
        var.all = new MnvH1D(*(MnvH1D*)var.f_mc->Get(hist_name.c_str())); 
    }else{ 
        hist_name = var.name + "_all";
        var.all = new MnvH1D(*(MnvH1D*)var.f_data->Get(hist_name.c_str())); 
    }    

    hist_name = var.name + "_mc_reco_signal";
    var.mc_reco_signal = new MnvH1D(*(MnvH1D*)var.f_mc->Get(hist_name.c_str())); 

    hist_name = var.name + "_mc_reco_bckg";
    var.mc_reco_bckg = new MnvH1D(*(MnvH1D*)var.f_mc->Get(hist_name.c_str())); 

    // All Signal is on Truth File
    hist_name = var.name + "_mc_truth_all_signal";
    var.mc_truth_all_signal = new MnvH1D(*(MnvH1D*)f_truth->Get(hist_name.c_str())); 
    
    hist_name = var.name + "_mc_truth_signal";
    var.mc_truth_signal = new MnvH1D(*(MnvH1D*)var.f_mc->Get(hist_name.c_str())); 

    hist_name = var.name + "_eff"; 
    var.eff = new MnvH1D(*var.mc_truth_signal);
    var.eff->SetName(hist_name.c_str());
    var.eff->Divide(var.mc_truth_signal, var.mc_truth_all_signal);

    hist_name = var.name + "_response";
    var.response = new MnvH2D(*(MnvH2D*)var.f_mc->Get(hist_name.c_str())); 
}

void CCProtonPi0_CrossSection::writeHistograms(XSec &var)
{
    // Data Histograms
    var.xsec->Write();
    var.all->Write();
    var.bckg_subtracted->Write();
    var.bckg_estimated->Write();
    var.unfolded->Write();
    var.efficiency_corrected->Write();
    var.integrated_flux->Write();
   
    // MC Truth Histograms
    var.mc_truth_all_signal->Write();
    var.mc_truth_signal->Write();
    var.mc_reco_signal->Write();
    var.mc_reco_bckg->Write();
    var.eff->Write();
    var.response->Write();
}

void CCProtonPi0_CrossSection::writeHistograms()
{
    std::cout<<">> Writing "<<rootDir_out<<std::endl;
   
    f_out->cd();
    
    // Pi0 Invariant Mass
    invMass_fit_result->Write();
    invMass_fit_bckg->Write();
    invMass_fit_signal->Write();
    invMass_all->Write();
    invMass_mc_reco_signal->Write();
    invMass_mc_reco_bckg->Write();
    
    writeHistograms(muon_P);
    writeHistograms(muon_theta);
    writeHistograms(pi0_P);
    writeHistograms(pi0_KE);
    writeHistograms(pi0_theta);
    writeHistograms(QSq);

    f_out->Close();
}

#endif

