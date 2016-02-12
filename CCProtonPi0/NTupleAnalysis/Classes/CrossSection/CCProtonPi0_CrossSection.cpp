#ifndef CCProtonPi0_CrossSection_cpp
#define CCProtonPi0_CrossSection_cpp

#include "CCProtonPi0_CrossSection.h"

using namespace PlotUtils;

CCProtonPi0_CrossSection::CCProtonPi0_CrossSection() : CCProtonPi0_NTupleAnalysis()
{
    std::cout<<"Initializing CCProtonPi0_CrossSection"<<std::endl;

    iteration = 1; // Unfolding Iteration
    min_invMass = 60; // MeV
    max_invMass = 200; // MeV
    
    data_POT = 9.44854e+19;   
    mc_POT = 9.0182e+20;
    
    OpenRootFiles();
    initHistograms();
    std::cout<<"Done!"<<std::endl;
}

void CCProtonPi0_CrossSection::Calc_CrossSections()
{
    Calc_Normalized_NBackground();
    Calc_CrossSection_Muon_P();   
    Calc_CrossSection_Pi0_P();   

    Label_CrossSection_Hists();

    // Re Init Default Histograms Again before Writing to File
    // We Scaled some of them during the data corrections
    initHistograms();
    writeHistograms();
}

void CCProtonPi0_CrossSection::Calc_CrossSection_Muon_P()
{
    std::cout<<"\n-----------------------------------------------------------------------"<<std::endl;
    std::cout<<"Calculating Cross Section for Muon Momentum"<<std::endl;

    // Data Correction
    data_bckg_subtracted_muon_P = Subtract_Background(data_all_muon_P, mc_reco_bckg_muon_P, data_bckg_estimated_muon_P,  "muon_P");
    data_unfolded_muon_P = Unfold_Data(data_bckg_subtracted_muon_P, response_muon_P, "muon_P");   
    data_efficiency_corrected_muon_P = Efficiency_Divide(data_unfolded_muon_P, eff_muon_P, "muon_P");   

    // Integrate Flux
    data_integrated_flux_muon_P = Integrate_Flux(data_efficiency_corrected_muon_P, "muon_P",false);
    mc_truth_integrated_flux_muon_P = Integrate_Flux(mc_truth_all_signal_muon_P, "muon_P",false);
    
    // Calculate Cross Section for Data
    data_xsec_muon_P = Calc_CrossSection(data_efficiency_corrected_muon_P, data_integrated_flux_muon_P, data_POT, false);
    data_xsec_muon_P->SetName("data_xsec_muon_P");
    data_xsec_muon_P->SetNormBinWidth(1.0);
    
    // Calculate Cross Section for MC Truth All 
    mc_truth_xsec_muon_P = Calc_CrossSection(mc_truth_all_signal_muon_P, mc_truth_integrated_flux_muon_P, mc_POT, true);
    mc_truth_xsec_muon_P->SetName("mc_truth_xsec_muon_P");
    mc_truth_xsec_muon_P->SetNormBinWidth(1.0);
    std::cout<<"Done!"<<std::endl;
    std::cout<<"-----------------------------------------------------------------------"<<std::endl;
}


void CCProtonPi0_CrossSection::Calc_CrossSection_Pi0_P()
{
    std::cout<<"\n-----------------------------------------------------------------------"<<std::endl;
    std::cout<<"Calculating Cross Section for Pi0 Momentum"<<std::endl;

    // Data Correction
    data_bckg_subtracted_pi0_P = Subtract_Background(data_all_pi0_P, mc_reco_bckg_pi0_P, data_bckg_estimated_pi0_P, "pi0_P");
    data_unfolded_pi0_P = Unfold_Data(data_bckg_subtracted_pi0_P, response_pi0_P, "pi0_P");   
    data_efficiency_corrected_pi0_P = Efficiency_Divide(data_unfolded_pi0_P, eff_pi0_P, "pi0_P");   

    // Integrate Flux
    data_integrated_flux_pi0_P = Integrate_Flux(data_efficiency_corrected_pi0_P, "pi0_P",false);
    mc_truth_integrated_flux_pi0_P = Integrate_Flux(mc_truth_all_signal_pi0_P, "pi0_P",false);

    // Calculate Cross Section for Data
    data_xsec_pi0_P = Calc_CrossSection(data_efficiency_corrected_pi0_P, data_integrated_flux_pi0_P, data_POT, false);
    data_xsec_pi0_P->SetName("data_xsec_pi0_P");
    data_xsec_pi0_P->SetNormBinWidth(0.1);

    // Calculate Cross Section for MC Truth All 
    mc_truth_xsec_pi0_P = Calc_CrossSection(mc_truth_all_signal_pi0_P, mc_truth_integrated_flux_pi0_P, mc_POT, true);
    mc_truth_xsec_pi0_P->SetName("mc_truth_xsec_pi0_P");
    mc_truth_xsec_pi0_P->SetNormBinWidth(0.1);

    std::cout<<"Done!"<<std::endl;
    std::cout<<"-----------------------------------------------------------------------"<<std::endl;
}

void CCProtonPi0_CrossSection::Calc_Normalized_NBackground()
{
    std::cout<<"\nCalculating N(Background) in Data"<<std::endl;
    
    // Get Histograms
    std::string rootDir_data = Folder_List::rootDir_CutHists_data;
    std::string rootDir_mc = Folder_List::rootDir_CutHists_mc;

    TFile* f_data = new TFile(rootDir_data.c_str());
    TFile* f_mc = new TFile(rootDir_mc.c_str());
    
    TH1D* data_invMass = (TH1D*)f_data->Get("data_invMass_All"); 
    TH1D* signal_invMass = (TH1D*)f_mc->Get("signal_invMass_All"); 
    TH1D* bckg_invMass = (TH1D*)f_mc->Get("bckg_invMass_All"); 
    
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
        std::cout<<"\tBackground in Data (Whole Spectrum) = "<<total_fitted_background<<std::endl;
        fitted_background->Scale(total_fitted_background/area); 

            // this calculation overestimates the number of background events
            // in the lower and upper bins
            // nbkg = fitted_background->Integral(__lower_bin, __upper_bin);

            // integrate the number of background events correctly
            // notice the mass range, not bin range
        N_Background_Data = Integrate_SignalRegion(fitted_background);

        std::cout<<"\tBackground in Data (Signal Region) = "<<N_Background_Data<<std::endl; 
    }

    delete mc_models;
    delete fitter;
    delete f_mc;
    delete f_data;
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
    std::string hist_name = "data_bckg_subtracted_" + var_name;
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
    hist_name = "data_bckg_estimated_" + var_name;
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
    std::string hist_name = "data_unfolded_" + var_name;
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
    std::string hist_name = "data_efficiency_corrected_" + var_name;
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

    std::string hist_name = "integrated_flux_" + var_name;
    integrated_flux->SetName(hist_name.c_str());

    std::cout<<"Done!"<<std::endl;

    return integrated_flux;
}

void CCProtonPi0_CrossSection::OpenRootFiles()
{
    // Flux File
    rootDir_flux = Folder_List::rootDir_Flux_new;

    // Output File
    rootDir_out = Folder_List::rootDir_CrossSection;
    std::cout<<"\tRoot File: "<<rootDir_out<<std::endl;
    f_out = new TFile(rootDir_out.c_str(), "RECREATE");

    // Truth File
    std::string rootDir = Folder_List::rootDir_Truth_mc;
    f_truth = new TFile(rootDir.c_str());

    // Data Files
    rootDir = Folder_List::rootDir_Muon_data;
    f_data_muon = new TFile(rootDir.c_str());
    rootDir = Folder_List::rootDir_Pion_data;
    f_data_pi0 = new TFile(rootDir.c_str());

    // MC Files
    rootDir = Folder_List::rootDir_Muon_mc;
    f_mc_muon = new TFile(rootDir.c_str());
    rootDir = Folder_List::rootDir_Pion_mc;
    f_mc_pi0 = new TFile(rootDir.c_str());
}

void CCProtonPi0_CrossSection::initHistograms()
{
    // ------------------------------------------------------------------------
    // Muon Momentum
    // ------------------------------------------------------------------------
    data_all_muon_P = new MnvH1D(*(MnvH1D*)f_data_muon->Get("data_all_muon_P")); 
    data_all_muon_P->SetTitle("Original Data P_{#mu}");

    mc_reco_all_muon_P = new MnvH1D(*(MnvH1D*)f_mc_muon->Get("mc_reco_all_muon_P")); 
    mc_reco_all_muon_P->SetTitle("MC Reconstructed P_{#mu}");
    mc_reco_all_muon_P->SetName("mc_reco_all_muon_P");

    mc_reco_signal_muon_P = new MnvH1D(*(MnvH1D*)f_mc_muon->Get("mc_reco_signal_muon_P")); 
    mc_reco_signal_muon_P->SetTitle("MC Reconstructed Signal P_{#mu}");

    mc_reco_bckg_muon_P = new MnvH1D(*(MnvH1D*)f_mc_muon->Get("mc_reco_bckg_muon_P")); 
    mc_reco_bckg_muon_P->SetTitle("MC Background P_{#mu}");

    mc_truth_all_signal_muon_P = new MnvH1D(*(MnvH1D*)f_truth->Get("mc_truth_all_signal_muon_P")); 
    mc_truth_all_signal_muon_P->SetTitle("MC Truth All Signal P_{#mu}");
    mc_truth_signal_muon_P = new MnvH1D(*(MnvH1D*)f_mc_muon->Get("mc_truth_signal_muon_P")); 
    mc_truth_signal_muon_P->SetTitle("MC Truth Signal P_{#mu}");

    response_muon_P = new MnvH2D(*(MnvH2D*)f_mc_muon->Get("response_P")); 
    response_muon_P->SetName("response_muon_P");

    eff_muon_P = new MnvH1D(*mc_truth_signal_muon_P);
    eff_muon_P->SetName("eff_muon_P");
    eff_muon_P->Divide(mc_truth_signal_muon_P, mc_truth_all_signal_muon_P);
    
    // ------------------------------------------------------------------------
    // Pi0 Momentum
    // ------------------------------------------------------------------------
    data_all_pi0_P = new MnvH1D(*(MnvH1D*)f_data_pi0->Get("data_all_pi0_P")); 
    data_all_pi0_P->SetTitle("Original Data P_{#pi^{0}}");
 
    mc_reco_all_pi0_P = new MnvH1D(*(MnvH1D*)f_mc_pi0->Get("mc_reco_all_pi0_P")); 
    mc_reco_all_pi0_P->SetTitle("MC Reconstructed P_{#mu}");
    mc_reco_all_pi0_P->SetName("mc_reco_all_pi0_P");
   
    mc_reco_signal_pi0_P = new MnvH1D(*(MnvH1D*)f_mc_pi0->Get("mc_reco_signal_pi0_P")); 
    mc_reco_signal_pi0_P->SetTitle("MC Reconstructed Signal P_{#pi^{0}}");
 
    mc_reco_bckg_pi0_P = new MnvH1D(*(MnvH1D*)f_mc_pi0->Get("mc_reco_bckg_pi0_P")); 
    mc_reco_bckg_pi0_P->SetTitle("MC Reconstructed Background P_{#pi^{0}}");
 
    mc_truth_all_signal_pi0_P = new MnvH1D(*(MnvH1D*)f_truth->Get("mc_truth_all_signal_pi0_P")); 
    mc_truth_all_signal_pi0_P->SetTitle("MC Truth All Signal P_{#pi^{0}}");
   
    mc_truth_signal_pi0_P = new MnvH1D(*(MnvH1D*)f_mc_pi0->Get("mc_truth_signal_pi0_P")); 
    mc_truth_signal_pi0_P->SetTitle("MC Truth Signal P_{#pi^{0}}");
    
    response_pi0_P = new MnvH2D(*(MnvH2D*)f_mc_pi0->Get("response_P")); 
    response_pi0_P->SetName("response_pi0_P");
 
    eff_pi0_P = new MnvH1D(*mc_truth_signal_pi0_P);
    eff_pi0_P->SetName("eff_pi0_P");
    eff_pi0_P->Divide(mc_truth_signal_pi0_P, mc_truth_all_signal_pi0_P);
    
}

MnvH1D* CCProtonPi0_CrossSection::Calc_CrossSection(MnvH1D* data_efficiency_corrected, MnvH1D* integrated_flux, double pot, bool isMC)
{
    MnvH1D* h_xs = new MnvH1D(*data_efficiency_corrected);

    // Divide the efficiency corrected distribution by the integrated flux
    h_xs->Divide(data_efficiency_corrected, integrated_flux);

    //fiducial is modules 27-80 inclusive
    int nplanes = 2 * ( 80 - 27 + 1 );
    double n_atoms    = TargetUtils::Get().GetTrackerNCarbonAtoms( nplanes, isMC, 850.0 );
    double n_nucleons = TargetUtils::Get().GetTrackerNNucleons(nplanes, isMC, 850.0);

    std::cout << "Number of C atoms   (1e30): " << n_atoms/1e30 << std::endl;
    std::cout << "Number of neutrons: (1e30): " << n_nucleons/1e30 << std::endl;

    h_xs->Scale(1/pot);
    h_xs->Scale(1/n_nucleons);

    h_xs->Scale(1e40); // to quote xs in 1-e40, 1e-42 for theta;

    return h_xs;
}

void CCProtonPi0_CrossSection::Label_CrossSection_Hists()
{
    data_xsec_muon_P->SetTitle("Differential Cross Section for P_{#mu}");
    data_xsec_muon_P->GetXaxis()->SetTitle("Muon Momentum [GeV]");
    data_xsec_muon_P->GetYaxis()->SetTitle("d#sigma/d_{P_{#mu}} (10^{-40} cm^{2}/nucleon/GeV)");

    mc_truth_xsec_muon_P->SetTitle("Differential Cross Section for P_{#mu}");
    mc_truth_xsec_muon_P->GetXaxis()->SetTitle("Muon Momentum [GeV]");
    mc_truth_xsec_muon_P->GetYaxis()->SetTitle("d#sigma/d_{P_{#mu}} (10^{-40} cm^{2}/nucleon/GeV)");

    data_xsec_pi0_P->SetTitle("Differential Cross Section for P_{#pi^{0}}");
    data_xsec_pi0_P->GetXaxis()->SetTitle("Pion Momentum [GeV]");
    data_xsec_pi0_P->GetYaxis()->SetTitle("d#sigma/d_{P_{#pi^{0}}} (10^{-40} cm^{2}/nucleon/GeV)");

    mc_truth_xsec_pi0_P->SetTitle("Differential Cross Section for P_{#pi^{0}}");
    mc_truth_xsec_pi0_P->GetXaxis()->SetTitle("Pion Momentum [GeV]");
    mc_truth_xsec_pi0_P->GetYaxis()->SetTitle("d#sigma/d_{P_{#pi^{0}}} (10^{-40} cm^{2}/nucleon/GeV)");
}

void CCProtonPi0_CrossSection::writeHistograms()
{
    std::cout<<">> Writing "<<rootDir_out<<std::endl;
   
    f_out->cd();
  
    // Muon Momentum
    data_xsec_muon_P->Write();
    data_all_muon_P->Write();
    data_bckg_subtracted_muon_P->Write();
    data_bckg_estimated_muon_P->Write();
    data_unfolded_muon_P->Write();
    data_efficiency_corrected_muon_P->Write();
    data_integrated_flux_muon_P->Write();
    
    mc_reco_all_muon_P->Write();
    mc_reco_signal_muon_P->Write();
    mc_reco_bckg_muon_P->Write();
    mc_truth_xsec_muon_P->Write();
    mc_truth_all_signal_muon_P->Write();
    mc_truth_signal_muon_P->Write();
    mc_truth_integrated_flux_muon_P->Write();
    
    response_muon_P->Write();
    eff_muon_P->Write();
    

    // Pi0 Momentum
    data_xsec_pi0_P->Write();
    data_all_pi0_P->Write();
    data_bckg_subtracted_pi0_P->Write();
    data_bckg_estimated_pi0_P->Write();
    data_unfolded_pi0_P->Write();
    data_efficiency_corrected_pi0_P->Write();
    data_integrated_flux_pi0_P->Write();
    
    mc_reco_all_pi0_P->Write();
    mc_reco_signal_pi0_P->Write();
    mc_reco_bckg_pi0_P->Write();
    mc_truth_xsec_pi0_P->Write();
    mc_truth_all_signal_pi0_P->Write();
    mc_truth_signal_pi0_P->Write();
    mc_truth_integrated_flux_pi0_P->Write();
    
    response_pi0_P->Write();
    eff_pi0_P->Write();


    f_out->Close();
}

#endif

