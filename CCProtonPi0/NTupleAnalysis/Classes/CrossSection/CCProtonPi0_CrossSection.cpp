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

    // Open output Log File
    if (m_isMC){
        text_out_name = Folder_List::output + Folder_List::textOut + "CrossSection_Log_MC.txt";
    }else{
        text_out_name = Folder_List::output + Folder_List::textOut + "CrossSection_Log_Data.txt";
    }
    OpenTextFile(text_out_name, text_out);

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
    Calc_CrossSection(W);   
    Calc_CrossSection(Enu);   
    
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
    var.flux_integrated = Integrate_Flux(var.efficiency_corrected, var.name);
    
    // Calculate Final Cross Section
    var.xsec = Calc_FinalCrossSection(var.flux_integrated, var.name); 
 
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
    var.xsec->SetNormBinWidth(var.smallest_bin_width);
    double norm_bin_width = var.xsec->GetNormBinWidth();
    double quote_width = 1.0; // per Unit (degree, GeV, GeV^2 etc...) 
    double scale = quote_width/norm_bin_width; 
    var.xsec->Scale(scale);
}

void CCProtonPi0_CrossSection::Calc_Normalized_NBackground()
{
    text_out<<"\nCalculating N(Background) in Data"<<std::endl;
   
    // ------------------------------------------------------------------------
    // Get Histograms for All Universes
    // ------------------------------------------------------------------------
    std::vector<TH1D*> data_invMass;
    std::vector<TH1D*> bckg_invMass;
    std::vector<TH1D*> signal_invMass;
    
    std::vector<std::string> dummy_err_bands;
    std::vector<int> dummy_hist_ind;

    GetAllVertUniverses(invMass_all, data_invMass, dummy_err_bands, dummy_hist_ind);
    // Clear Dummy Vectors
    dummy_err_bands.clear();
    dummy_hist_ind.clear();
    GetAllVertUniverses(invMass_mc_reco_bckg, bckg_invMass, dummy_err_bands, dummy_hist_ind);
    // Clear Dummy Vectors
    dummy_err_bands.clear();
    dummy_hist_ind.clear();

    // Get signal_invMass (Background Subtracted)
    //      Start as data hist (will subtract bckg to get correct signal in data)
    for (unsigned int i = 0; i < data_invMass.size(); ++i){
        TH1D* temp = new TH1D(*data_invMass[i]);
        temp->SetName("invMass_fit_signal");
        signal_invMass.push_back(temp);
    }
    
    // POT Normalize MC Background Histogram if we are using Data
    if (!m_isMC){
        for (unsigned int i = 0; i < bckg_invMass.size(); ++i){
            bckg_invMass[i]->Scale(POT_ratio);    
        }
    }

    // ------------------------------------------------------------------------
    // Estimate N(Background) in Data in each universe
    // ------------------------------------------------------------------------
    for (unsigned int i = 0; i < data_invMass.size(); ++i){
        text_out<<"Estimating Backgroung in Universe: "<<i<<std::endl;
        // Get N(Events)
        double total_data = data_invMass[i]->Integral();
        double total_bckg = bckg_invMass[i]->Integral();

        text_out<<"\tAll Events in Data (Whole Spectrum) = "<<total_data<<std::endl;
        text_out<<"\tBackground in Data (Whole Spectrum) = "<<total_bckg<<std::endl;

        // Subtract bckg hist from data hist
        signal_invMass[i]->Add(bckg_invMass[i],-1); 

        // this calculation overestimates the number of background events
        // in the lower and upper bins
        // nbkg = fitted_background->Integral(__lower_bin, __upper_bin);

        // integrate the number of background events correctly
        // notice the mass range, not bin range
        N_Background_Data.push_back(Integrate_SignalRegion(bckg_invMass[i]));
        double N_Signal_Data = Integrate_SignalRegion(signal_invMass[i]);
        double N_All_Data = Integrate_SignalRegion(data_invMass[i]);

        double percent_signal = N_Signal_Data/N_All_Data*100;
        double percent_bckg = N_Background_Data[i]/N_All_Data*100.0;

        text_out<<std::endl;
        text_out<<"\tAll Events in Data (Signal Region) = "<<N_All_Data<<std::endl; 
        text_out<<"\tBackground in Data (Signal Region) = "<<N_Background_Data[i]<<" "<<percent_bckg<<std::endl; 
        text_out<<"\tSignal in Data (Signal Region) = "<<N_Signal_Data<<" "<<percent_signal<<std::endl; 
        text_out<<std::endl;
        text_out<<std::endl;
    }

    // Clear Allocated Memory
    ClearAllUniversesVector(data_invMass);
    ClearAllUniversesVector(bckg_invMass);
    ClearAllUniversesVector(signal_invMass);
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

void CCProtonPi0_CrossSection::NormalizeHistogram(TH1D* h)
{
    text_out<<"\tNormalizing Background Shape on TH1D"<<std::endl;
    int NBins = h->GetNbinsX();
    double area = h->Integral();
    double nOverFlow = h->GetBinContent(NBins+1);
    text_out<<"\t\tBefore Norm = "<<area<<std::endl;
    h->Scale(1/(area+nOverFlow));
    text_out<<"\t\tAfter Norm = "<<h->Integral()<<std::endl;
    text_out<<"\tDone!"<<std::endl;
}

void CCProtonPi0_CrossSection::NormalizeHistogram(MnvH1D* h)
{
    text_out<<"\tNormalizing Background Shape on MnvH1D"<<std::endl;
    int NBins = h->GetNbinsX();
    double area = h->Integral();
    double nOverFlow = h->GetBinContent(NBins+1);
    text_out<<"\t\tBefore Norm = "<<area<<std::endl;
    h->Scale(1/(area+nOverFlow),"",false); // Scale only on CentralValue
    text_out<<"\t\tAfter Norm = "<<h->Integral()<<std::endl;
    text_out<<"\tDone!"<<std::endl;
}

MnvH1D* CCProtonPi0_CrossSection::Subtract_Background(MnvH1D* data, MnvH1D* mc_bckg, MnvH1D* &bckg_estimated, std::string var_name)
{
    text_out<<"Subtracting Background for "<<var_name<<std::endl;

    // Init MnvH1D
    MnvH1D* bckg_subtracted = new MnvH1D(*data); 
    std::string hist_name = var_name + "_bckg_subtracted";
    bckg_subtracted->SetName(hist_name.c_str());
    bckg_subtracted->SetTitle("Background Subtracted Data");

    // ------------------------------------------------------------------------
    //  Subtract Background from Central Value -- No Universe Subtraction
    // ------------------------------------------------------------------------
    text_out<<"\tSubtracting Background in Central Value"<<std::endl;

    // Get Shape of Background
    //      Unit area -- Including Overflow
    //      After Normalization, Area MUST be less than 1
    NormalizeHistogram(mc_bckg);

    // Estimate N Background in Variable
    // [0] is for CV Value
    mc_bckg->Scale(N_Background_Data[0],"",false);

    // Subtracted Background -- Use TH1D::Add to add only to Central Value Histogram
    bckg_subtracted->TH1D::Add(mc_bckg,-1);

    text_out<<std::endl;
    text_out<<"\tTotal Data Area = "<<data->Integral()<<std::endl;
    text_out<<"\tBackground Subtracted Data Area = "<<bckg_subtracted->Integral()<<std::endl;
    text_out<<std::endl;
 
    // ------------------------------------------------------------------------
    // Get Pointers to All Universes
    // ------------------------------------------------------------------------
    // >> These are "Pointers" NOT "new" histograms
    // >> This vector does "NOT" include Central Value Histogram
    std::vector<TH1D*> data_all_universes;
    std::vector<TH1D*> bckg_subtracted_all_universes;
    std::vector<TH1D*> mc_bckg_all_universes;
    
    GetPointersAllVertUniverses(data, data_all_universes);
    GetPointersAllVertUniverses(bckg_subtracted, bckg_subtracted_all_universes);
    GetPointersAllVertUniverses(mc_bckg, mc_bckg_all_universes);

    text_out<<"N(Universes) = "<<data_all_universes.size()<<std::endl;
   
    // Sanity Check
    if ( bckg_subtracted_all_universes.size() != mc_bckg_all_universes.size()){
        text_out<<"WARNING! - Subtract Background N(Universes) NOT Same!"<<std::endl;
        exit(1);
    }
    
    // Loop Over All Universes
    for (unsigned int i = 0; i < bckg_subtracted_all_universes.size(); ++i){
        text_out<<"\tSubtracting Background in Universe = "<<i<<std::endl;
       
        // Get Shape of Background
        //      Unit area -- Including Overflow
        //      After Normalization, Area MUST be less than 1
        NormalizeHistogram(mc_bckg_all_universes[i]);

        // Estimate N Background in Variable
        mc_bckg_all_universes[i]->Scale(N_Background_Data[i+1]);

        text_out<<"\tBackground Subtracted Data = "<<bckg_subtracted_all_universes[i]->Integral()<<std::endl;
        // Subtracted Background
        bckg_subtracted_all_universes[i]->Add(mc_bckg_all_universes[i],-1);
       
        text_out<<"\tTotal Data Area = "<<data_all_universes[i]->Integral()<<std::endl;
        text_out<<"\tEstimated Background in Data = "<<mc_bckg_all_universes[i]->Integral()<<std::endl;
        text_out<<"\tBackground Subtracted Data = "<<bckg_subtracted_all_universes[i]->Integral()<<std::endl;
        text_out<<std::endl;
    }

    // Estimated Background
    bckg_estimated = new MnvH1D(*mc_bckg);
    hist_name = var_name + "_bckg_estimated";
    bckg_estimated->SetName(hist_name.c_str());
    bckg_estimated->SetTitle("Estimated Background in Data"); 
   
    text_out<<"Done!"<<std::endl;

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

MnvH1D* CCProtonPi0_CrossSection::Integrate_Flux(MnvH1D* data_efficiency_corrected, std::string var_name)
{
    std::cout<<"Integrating Flux for "<<var_name<<std::endl;

    MnvH1D* flux_integrated = new MnvH1D (*data_efficiency_corrected); 
   
    std::string hist_name = var_name + "_flux_integrated";
    flux_integrated->SetName(hist_name.c_str());
    flux_integrated->SetTitle("Flux Integrated Data");

    // Divide by total integral except Enu
    // Enu divided bin by bin
    if (var_name.compare("Enu") == 0 ){
        std::cout<<"Flux Integration for Neutrino Energy!"<<std::endl; 
        flux_integrated->Divide(data_efficiency_corrected, h_flux_rebinned);
    }else{
        flux_integrated->Scale(1/flux_integral);
    }
    std::cout<<"Done!"<<std::endl;

    return flux_integrated;
}

MnvH1D* CCProtonPi0_CrossSection::Calc_FinalCrossSection(MnvH1D* flux_integrated, std::string var_name)
{
    MnvH1D* h_xs = new MnvH1D(*flux_integrated);
    std::string hist_name = var_name + "_xsec";
    h_xs->SetName(hist_name.c_str());

    // Set POT
    double pot;
    if (m_isMC) pot = mc_POT;
    else pot = data_POT;

    //fiducial is modules 27-79 inclusive
    double det_mass     = TargetUtils::Get().GetTrackerMass(5991.37, 8363.92, m_isMC, 850.0);
    double n_atoms      = TargetUtils::Get().GetTrackerNCarbonAtoms(5991.37, 8363.92, m_isMC, 850.0);
    double n_nucleons   = TargetUtils::Get().GetTrackerNNucleons(5991.37, 8363.92, m_isMC, 850.0);

    std::cout<<"Detector Mass         = "<<det_mass<<std::endl;
    std::cout<<"Number of C atoms     = "<<n_atoms/1e30<<std::endl;
    std::cout<<"Number of Nucleons    = "<<n_nucleons/1e30<<std::endl;

    h_xs->Scale(1/pot);
    h_xs->Scale(1/n_nucleons);

    h_xs->Scale(1e40); // to quote xs in 1-e40;

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
    init_W();
    init_Enu();
}

void CCProtonPi0_CrossSection::OpenRootFiles()
{
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
    initFluxHistograms();
    
    // Pi0 Invariant Mass
    if (m_isMC) invMass_all = GetMnvH1D(f_mc_cutHists, "invMass_mc_reco_all");
    else invMass_all = GetMnvH1D(f_data_cutHists, "invMass_all");

    invMass_mc_reco_signal = GetMnvH1D(f_mc_cutHists, "invMass_mc_reco_signal");
    invMass_mc_reco_bckg = GetMnvH1D(f_mc_cutHists, "invMass_mc_reco_bckg");
}

void CCProtonPi0_CrossSection::initFluxHistograms()
{
    // Get Reweighted Flux Histogram
    // Use minervaLE_FHC because I compared others and all are almost same
    // Integral of this is the average
    delete frw;
    frw = new FluxReweighter(14, applyNuEConstraint, FluxReweighter::minervaLE_FHC, new_flux, old_flux);
    h_flux_minervaLE_FHC = new MnvH1D (*(frw->GetFluxReweighted(14)));
    h_flux_minervaLE_FHC->SetName("h_flux_minervaLE_FHC");
    h_flux_minervaLE_FHC->Scale(1/mSq_to_cmSq); // Our measurement scale is cm2
    flux_integral = h_flux_minervaLE_FHC->Integral(4,30,"width");
    std::cout<<"Signal Region Flux Integral = "<<flux_integral<<std::endl; 
   
    // Get Rebinned Flux for Neutrino Energy
    h_flux_rebinned = new MnvH1D( "h_flux_rebinned","Flux (rebinned)", binList.size_Enu, binList.a_Enu);
    h_flux_rebinned->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    h_flux_rebinned->GetYaxis()->SetTitle("#nu_{#mu} s/cm^{2}/P.O.T./GeV");
    AddVertErrorBands_MC(h_flux_rebinned);

    int nBins = h_flux_rebinned->GetNbinsX();
    for (int i = 1; i <= nBins; i++){
        double low = h_flux_rebinned->GetBinLowEdge(i);
        double up = h_flux_rebinned->GetBinLowEdge(i+1);
        double content = GetFluxHistContent(h_flux_minervaLE_FHC,low,up);
        double bin_width = h_flux_rebinned->GetBinWidth(i);
        h_flux_rebinned->SetBinContent(i,content/bin_width);
    }
}

void CCProtonPi0_CrossSection::initHistograms(XSec &var)
{
    std::string hist_name;
    
    if (m_isMC){ 
        hist_name = var.name + "_mc_reco_all";
        var.all = GetMnvH1D(var.f_mc, hist_name);
    }else{ 
        hist_name = var.name + "_all";
        var.all = GetMnvH1D(var.f_data, hist_name);
    }    

    hist_name = var.name + "_mc_reco_signal";
    var.mc_reco_signal = GetMnvH1D(var.f_mc, hist_name);

    hist_name = var.name + "_mc_reco_bckg";
    var.mc_reco_bckg = GetMnvH1D(var.f_mc, hist_name);

    // All Signal is on Truth File
    hist_name = var.name + "_mc_truth_all_signal";
    var.mc_truth_all_signal = GetMnvH1D(f_truth, hist_name);
    
    hist_name = var.name + "_mc_truth_signal";
    var.mc_truth_signal = GetMnvH1D(var.f_mc, hist_name);

    hist_name = var.name + "_eff"; 
    var.eff = new MnvH1D(*var.mc_truth_signal);
    var.eff->SetName(hist_name.c_str());
    var.eff->Divide(var.mc_truth_signal, var.mc_truth_all_signal);

    hist_name = var.name + "_response";
    var.response = GetMnvH2D(var.f_mc, hist_name);
}

double CCProtonPi0_CrossSection::GetFluxHistContent(MnvH1D* hist, double low1, double low2)
{
    double total = 0.0;
    int nBins = hist->GetNbinsX();
    for (int i = 1; i <= nBins; ++i){
        double current_low = hist->GetBinLowEdge(i); 
        if (current_low < low1) continue;
        if (current_low == low2) break;
        total += hist->GetBinContent(i)*hist->GetBinWidth(i);
    }
    return total;
}

double CCProtonPi0_CrossSection::GetSmallestBinWidth(MnvH1D* hist)
{
    double smallest = 99999999;
    int nBins = hist->GetNbinsX();
    for (int i = 0; i <= nBins; ++i){
        double current = hist->GetBinWidth(i);
        if (current < smallest) smallest = current;
    }

    return smallest;
}

void CCProtonPi0_CrossSection::writeHistograms(XSec &var)
{
    // Data Histograms
    var.xsec->Write();
    
    var.all->SetNormBinWidth(var.smallest_bin_width);
    var.all->Write();
    
    var.bckg_subtracted->SetNormBinWidth(var.smallest_bin_width);
    var.bckg_subtracted->Write();

    var.bckg_estimated->SetNormBinWidth(var.smallest_bin_width);
    var.bckg_estimated->Write();

    var.unfolded->SetNormBinWidth(var.smallest_bin_width);
    var.unfolded->Write();

    var.efficiency_corrected->SetNormBinWidth(var.smallest_bin_width);
    var.efficiency_corrected->Write();

    var.flux_integrated->SetNormBinWidth(var.smallest_bin_width);
    var.flux_integrated->Write();
   
    // MC Truth Histograms
    var.mc_truth_all_signal->SetNormBinWidth(var.smallest_bin_width);
    var.mc_truth_all_signal->Write();

    var.mc_truth_signal->SetNormBinWidth(var.smallest_bin_width);
    var.mc_truth_signal->Write();

    var.mc_reco_signal->SetNormBinWidth(var.smallest_bin_width);
    var.mc_reco_signal->Write();
    
    var.mc_reco_bckg->SetNormBinWidth(var.smallest_bin_width);
    var.mc_reco_bckg->Write();

    var.eff->SetNormBinWidth(var.smallest_bin_width);
    var.eff->Write();

    var.response->Write();
}

void CCProtonPi0_CrossSection::writeHistograms()
{
    std::cout<<">> Writing "<<rootDir_out<<std::endl;
   
    f_out->cd();

    h_flux_minervaLE_FHC->Write();
    h_flux_rebinned->Write();

    // Pi0 Invariant Mass
    invMass_all->Write();
    invMass_mc_reco_signal->Write();
    invMass_mc_reco_bckg->Write();
    
    writeHistograms(muon_P);
    writeHistograms(muon_theta);
    writeHistograms(pi0_P);
    writeHistograms(pi0_KE);
    writeHistograms(pi0_theta);
    writeHistograms(QSq);
    writeHistograms(W);
    writeHistograms(Enu);

    f_out->Close();
}

#endif

