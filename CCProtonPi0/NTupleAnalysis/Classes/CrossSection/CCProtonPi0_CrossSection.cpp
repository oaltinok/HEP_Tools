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
    
    data_POT = 9.44854e+19;   
    mc_POT = 9.0182e+20;
    
    m_isMC = isMC;

    OpenRootFiles();
    initHistograms();
    std::cout<<"Done!"<<std::endl;
}

void CCProtonPi0_CrossSection::Calc_CrossSections()
{
    Calc_Normalized_NBackground();
    
    Calc_CrossSection_muon_P();   
    Calc_CrossSection_muon_theta();   
    
    Calc_CrossSection_pi0_P();   

    // Re Init Default Histograms Again before Writing to File
    // We Scaled some of them during the data corrections
    initHistograms();
    writeHistograms();
}

void CCProtonPi0_CrossSection::Calc_CrossSection_muon_P()
{
    std::cout<<"\n-----------------------------------------------------------------------"<<std::endl;
    std::cout<<"Calculating Cross Section for Muon Momentum"<<std::endl;
    bool isEv = false;

    // Data Correction
    muon_P_bckg_subtracted = Subtract_Background(muon_P_all, muon_P_mc_reco_bckg, muon_P_bckg_estimated,"muon_P");
    muon_P_unfolded = Unfold_Data(muon_P_bckg_subtracted, muon_P_response, "muon_P");   
    muon_P_efficiency_corrected = Efficiency_Divide(muon_P_unfolded, muon_P_eff, "muon_P");   

    // Integrate Flux
    muon_P_integrated_flux = Integrate_Flux(muon_P_efficiency_corrected, "muon_P", isEv);
    
    // Calculate Final Cross Section
    muon_P_xsec = Calc_FinalCrossSection(muon_P_efficiency_corrected, muon_P_integrated_flux, "muon_P");
  
    Style_muon_P();
    std::cout<<"Done!"<<std::endl;
    std::cout<<"-----------------------------------------------------------------------"<<std::endl;
}

void CCProtonPi0_CrossSection::Style_muon_P()
{
    // Add Labels
    muon_P_xsec->SetTitle("Differential Cross Section for P_{#mu}");
    muon_P_xsec->GetXaxis()->SetTitle("Muon Momentum [GeV]");
    muon_P_xsec->GetYaxis()->SetTitle("d#sigma/d_{P_{#mu}} (10^{-40} cm^{2}/nucleon/GeV)");
    
    // Scale Cross Section Result to match with Label
    double bin_width = 1.0; // GeV
    double quote_width = 1.0; // per GeV
    double scale = quote_width / bin_width; 
    muon_P_xsec->Scale(scale);
}

void CCProtonPi0_CrossSection::Calc_CrossSection_muon_theta()
{
    std::cout<<"\n-----------------------------------------------------------------------"<<std::endl;
    std::cout<<"Calculating Cross Section for Muon Theta"<<std::endl;
    bool isEv = false;

    // Data Correction
    muon_theta_bckg_subtracted = Subtract_Background(muon_theta_all, muon_theta_mc_reco_bckg, muon_theta_bckg_estimated,"muon_theta");
    muon_theta_unfolded = Unfold_Data(muon_theta_bckg_subtracted, muon_theta_response, "muon_theta");   
    muon_theta_efficiency_corrected = Efficiency_Divide(muon_theta_unfolded, muon_theta_eff, "muon_theta");   

    // Integrate Flux
    muon_theta_integrated_flux = Integrate_Flux(muon_theta_efficiency_corrected, "muon_theta", isEv);
    
    // Calculate Final Cross Section
    muon_theta_xsec = Calc_FinalCrossSection(muon_theta_efficiency_corrected, muon_theta_integrated_flux, "muon_theta");
  
    Style_muon_theta();

    std::cout<<"Done!"<<std::endl;
    std::cout<<"-----------------------------------------------------------------------"<<std::endl;
}

void CCProtonPi0_CrossSection::Style_muon_theta()
{
    // Add Labels
    muon_theta_xsec->SetTitle("Differential Cross Section for #theta_{#mu}");
    muon_theta_xsec->GetXaxis()->SetTitle("Muon Theta [degree]");
    muon_theta_xsec->GetYaxis()->SetTitle("d#sigma/d_{#theta_{#mu}} (10^{-40} cm^{2}/nucleon/degree)");
    
    // Scale Cross Section Result to match with Label
    double bin_width = 2.08; // degree 
    double quote_width = 1.0; // per degree
    double scale = quote_width / bin_width; 
    muon_theta_xsec->Scale(scale);
}

void CCProtonPi0_CrossSection::Calc_CrossSection_pi0_P()
{
    std::cout<<"\n-----------------------------------------------------------------------"<<std::endl;
    std::cout<<"Calculating Cross Section for Pi0 Momentum Data"<<std::endl;
    bool isEv = false;

    // Data Correction
    pi0_P_bckg_subtracted = Subtract_Background(pi0_P_all, pi0_P_mc_reco_bckg, pi0_P_bckg_estimated, "pi0_P");
    pi0_P_unfolded = Unfold_Data(pi0_P_bckg_subtracted, pi0_P_response, "pi0_P");   
    pi0_P_efficiency_corrected = Efficiency_Divide(pi0_P_unfolded, pi0_P_eff, "pi0_P");   

    // Integrate Flux
    pi0_P_integrated_flux = Integrate_Flux(pi0_P_efficiency_corrected, "pi0_P", isEv);

    // Calculate Final Cross Section
    pi0_P_xsec = Calc_FinalCrossSection(pi0_P_efficiency_corrected, pi0_P_integrated_flux, "pi0_P");
    pi0_P_xsec->SetNormBinWidth(0.1);

    Style_pi0_P();

    std::cout<<"Done!"<<std::endl;
    std::cout<<"-----------------------------------------------------------------------"<<std::endl;
}

void CCProtonPi0_CrossSection::Style_pi0_P()
{
    // Add Labels
    pi0_P_xsec->SetTitle("Differential Cross Section for P_{#pi^{0}}");
    pi0_P_xsec->GetXaxis()->SetTitle("Pion Momentum [GeV]");
    pi0_P_xsec->GetYaxis()->SetTitle("d#sigma/d_{P_{#pi^{0}}} (10^{-40} cm^{2}/nucleon/GeV)");

    // Style Cross Section Result to match with Label
    double bin_width = 0.1; // GeV
    double quote_width = 1.0; // per GeV
    double scale = quote_width / bin_width; 
    pi0_P_xsec->Scale(scale);

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
    
    // ------------------------------------------------------------------------
    // Estimate Background using TFractionFitter
    // ------------------------------------------------------------------------
    TObjArray* mc_models = new TObjArray(2);
    mc_models->Add(signal_invMass);
    mc_models->Add(bckg_invMass);

    TFractionFitter* fitter = new TFractionFitter(data_invMass, mc_models, "Q");
    fitter->Constrain(0, 0.0, 1.0);
    fitter->Constrain(1, 0.0, 1.0);
    
    std::cout<<"\tFitting Background Shape..."<<std::endl;
    Int_t status = fitter->Fit();

    if (status != 0) {
        std::cerr<<"\t\tFit Error!"<<std::endl;
    }
    
    fit_result  = new TH1F (*(TH1F*)fitter->GetPlot());
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
        TH1D* fitted_background = (TH1D*) bckg_invMass->Clone("fitted_bkg");
        double area = fitted_background->Integral();
        double total_fitted_background = f1 * N_Data; // Total N(Bckg) in Data
        std::cout<<"\tAll Events in Data (Whole Spectrum) = "<<N_Data<<std::endl;
        std::cout<<"\tBackground in Data (Whole Spectrum) = "<<total_fitted_background<<std::endl;
        fitted_background->Scale(total_fitted_background/area); 

            // this calculation overestimates the number of background events
            // in the lower and upper bins
            // nbkg = fitted_background->Integral(__lower_bin, __upper_bin);

            // integrate the number of background events correctly
            // notice the mass range, not bin range
        double N_All_Data = Integrate_SignalRegion(data_invMass);
        N_Background_Data = Integrate_SignalRegion(fitted_background);

        std::cout<<std::endl;
        std::cout<<"\tAll Events in Data (Signal Region) = "<<N_All_Data<<std::endl; 
        std::cout<<"\tBackground in Data (Signal Region) = "<<N_Background_Data<<std::endl; 
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
    
    bool apply_Enu_cut = false;
    MnvH1D* integrated_flux = calc_flux(data_efficiency_corrected, rootDir_flux, apply_Enu_cut, isEv);
    
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
    double n_atoms    = TargetUtils::Get().GetTrackerNCarbonAtoms( nplanes, m_isMC, 850.0 );
    double n_nucleons = TargetUtils::Get().GetTrackerNNucleons(nplanes, m_isMC, 850.0);

    std::cout << "Number of C atoms   (1e30): " << n_atoms/1e30 << std::endl;
    std::cout << "Number of neutrons: (1e30): " << n_nucleons/1e30 << std::endl;

    h_xs->Scale(1/pot);
    h_xs->Scale(1/n_nucleons);

    h_xs->Scale(1e40); // to quote xs in 1-e40, 1e-42 for theta;

    return h_xs;
}


void CCProtonPi0_CrossSection::OpenRootFiles()
{
    std::cout<<"Opening Root Files..."<<std::endl;
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

    // Data Files
    rootDir = Folder_List::rootDir_CutHists_data;
    f_data_cutHists = new TFile(rootDir.c_str());
    rootDir = Folder_List::rootDir_Muon_data;
    f_data_muon = new TFile(rootDir.c_str());
    rootDir = Folder_List::rootDir_Pion_data;
    f_data_pi0 = new TFile(rootDir.c_str());

    // MC Files
    rootDir = Folder_List::rootDir_CutHists_mc;
    f_mc_cutHists = new TFile(rootDir.c_str());
    rootDir = Folder_List::rootDir_Muon_mc;
    f_mc_muon = new TFile(rootDir.c_str());
    rootDir = Folder_List::rootDir_Pion_mc;
    f_mc_pi0 = new TFile(rootDir.c_str());

    std::cout<<"Done!"<<std::endl;
}

void CCProtonPi0_CrossSection::initHistograms()
{
    std::cout<<"Initializing Histograms..."<<std::endl;
    // ------------------------------------------------------------------------
    // Pi0 Invariant Mass
    // ------------------------------------------------------------------------
    if (m_isMC) invMass_all = new MnvH1D(*(MnvH1D*)f_mc_cutHists->Get("invMass_mc_reco_all"));
    else invMass_all = new MnvH1D(*(MnvH1D*)f_data_cutHists->Get("invMass_all"));

    invMass_mc_reco_signal = new MnvH1D(*(MnvH1D*)f_mc_cutHists->Get("invMass_mc_reco_signal")); 
    invMass_mc_reco_bckg  = new MnvH1D(*(MnvH1D*)f_mc_cutHists->Get("invMass_mc_reco_bckg")); 

    // ------------------------------------------------------------------------
    // Muon Momentum
    // ------------------------------------------------------------------------
    if (m_isMC) muon_P_all = new MnvH1D(*(MnvH1D*)f_mc_muon->Get("muon_P_mc_reco_all")); 
    else muon_P_all = new MnvH1D(*(MnvH1D*)f_data_muon->Get("muon_P_all")); 
    
    muon_P_mc_reco_signal = new MnvH1D(*(MnvH1D*)f_mc_muon->Get("muon_P_mc_reco_signal")); 
    muon_P_mc_reco_bckg = new MnvH1D(*(MnvH1D*)f_mc_muon->Get("muon_P_mc_reco_bckg")); 

    muon_P_mc_truth_all_signal = new MnvH1D(*(MnvH1D*)f_truth->Get("muon_P_mc_truth_all_signal")); 
    muon_P_mc_truth_signal = new MnvH1D(*(MnvH1D*)f_mc_muon->Get("muon_P_mc_truth_signal")); 
 
    muon_P_eff = new MnvH1D(*muon_P_mc_truth_signal);
    muon_P_eff->SetName("muon_P_eff");
    muon_P_eff->Divide(muon_P_mc_truth_signal, muon_P_mc_truth_all_signal);

    muon_P_response = new MnvH2D(*(MnvH2D*)f_mc_muon->Get("muon_P_response")); 
 
    // ------------------------------------------------------------------------
    // Muon Theta 
    // ------------------------------------------------------------------------
    if (m_isMC) muon_theta_all = new MnvH1D(*(MnvH1D*)f_mc_muon->Get("muon_theta_mc_reco_all")); 
    else muon_theta_all = new MnvH1D(*(MnvH1D*)f_data_muon->Get("muon_theta_all")); 
    
    muon_theta_mc_reco_signal = new MnvH1D(*(MnvH1D*)f_mc_muon->Get("muon_theta_mc_reco_signal")); 
    muon_theta_mc_reco_bckg = new MnvH1D(*(MnvH1D*)f_mc_muon->Get("muon_theta_mc_reco_bckg")); 

    muon_theta_mc_truth_all_signal = new MnvH1D(*(MnvH1D*)f_truth->Get("muon_theta_mc_truth_all_signal")); 
    muon_theta_mc_truth_signal = new MnvH1D(*(MnvH1D*)f_mc_muon->Get("muon_theta_mc_truth_signal")); 
 
    muon_theta_eff = new MnvH1D(*muon_theta_mc_truth_signal);
    muon_theta_eff->SetName("muon_theta_eff");
    muon_theta_eff->Divide(muon_theta_mc_truth_signal, muon_theta_mc_truth_all_signal);

    muon_theta_response = new MnvH2D(*(MnvH2D*)f_mc_muon->Get("muon_theta_response")); 
     
    // ------------------------------------------------------------------------
    // Pi0 Momentum
    // ------------------------------------------------------------------------
    if (m_isMC) pi0_P_all = new MnvH1D(*(MnvH1D*)f_mc_pi0->Get("pi0_P_mc_reco_all")); 
    else pi0_P_all = new MnvH1D(*(MnvH1D*)f_data_pi0->Get("pi0_P_all")); 
 
    pi0_P_mc_reco_signal = new MnvH1D(*(MnvH1D*)f_mc_pi0->Get("pi0_P_mc_reco_signal")); 
    pi0_P_mc_reco_bckg = new MnvH1D(*(MnvH1D*)f_mc_pi0->Get("pi0_P_mc_reco_bckg")); 
   
    pi0_P_mc_truth_all_signal = new MnvH1D(*(MnvH1D*)f_truth->Get("pi0_P_mc_truth_all_signal")); 
    pi0_P_mc_truth_signal = new MnvH1D(*(MnvH1D*)f_mc_pi0->Get("pi0_P_mc_truth_signal")); 
   
    pi0_P_eff = new MnvH1D(*pi0_P_mc_truth_signal);
    pi0_P_eff->SetName("pi0_P_eff");
    pi0_P_eff->Divide(pi0_P_mc_truth_signal, pi0_P_mc_truth_all_signal);

    pi0_P_response = new MnvH2D(*(MnvH2D*)f_mc_pi0->Get("pi0_P_response")); 
 
    std::cout<<"Done!"<<std::endl;
}

void CCProtonPi0_CrossSection::writeHistograms()
{
    std::cout<<">> Writing "<<rootDir_out<<std::endl;
   
    f_out->cd();
    
    // Pi0 Invariant Mass
    fit_result->Write();
    invMass_all->Write();
    invMass_mc_reco_signal->Write();
    invMass_mc_reco_bckg->Write();

    // Muon Momentum - Data
    muon_P_xsec->Write();
    muon_P_all->Write();
    muon_P_bckg_subtracted->Write();
    muon_P_bckg_estimated->Write();
    muon_P_unfolded->Write();
    muon_P_efficiency_corrected->Write();
    muon_P_integrated_flux->Write();
   
    // Muon Momentum - MC Truth
    muon_P_mc_truth_all_signal->Write();
    muon_P_mc_truth_signal->Write();
    muon_P_mc_reco_signal->Write();
    muon_P_mc_reco_bckg->Write();
    muon_P_eff->Write();
    muon_P_response->Write();
 
    // Muon Momentum - Data
    muon_theta_xsec->Write();
    muon_theta_all->Write();
    muon_theta_bckg_subtracted->Write();
    muon_theta_bckg_estimated->Write();
    muon_theta_unfolded->Write();
    muon_theta_efficiency_corrected->Write();
    muon_theta_integrated_flux->Write();
   
    // Muon Momentum - MC Truth
    muon_theta_mc_truth_all_signal->Write();
    muon_theta_mc_truth_signal->Write();
    muon_theta_mc_reco_signal->Write();
    muon_theta_mc_reco_bckg->Write();
    muon_theta_eff->Write();
    muon_theta_response->Write();
   
    // Pi0 Momentum - Data
    pi0_P_xsec->Write();
    pi0_P_all->Write();
    pi0_P_bckg_subtracted->Write();
    pi0_P_bckg_estimated->Write();
    pi0_P_unfolded->Write();
    pi0_P_efficiency_corrected->Write();
    pi0_P_integrated_flux->Write();

    // Pi0 Momentum - MC Truth
    pi0_P_mc_truth_all_signal->Write();
    pi0_P_mc_truth_signal->Write();
    pi0_P_mc_reco_signal->Write();
    pi0_P_mc_reco_bckg->Write();
    pi0_P_eff->Write();
    pi0_P_response->Write();

    f_out->Close();
}

#endif

