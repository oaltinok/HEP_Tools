#ifndef CCProtonPi0_CrossSection_cpp
#define CCProtonPi0_CrossSection_cpp

#include "CCProtonPi0_CrossSection.h"

using namespace PlotUtils;

CCProtonPi0_CrossSection::CCProtonPi0_CrossSection(bool isMC) : CCProtonPi0_NTupleAnalysis()
{
    std::cout<<"Initializing CCProtonPi0_CrossSection"<<std::endl;

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

    // Results w/o systematics
    RemoveErrorBands = false;

    OpenRootFiles();
    initHistograms();
    initXSecs();

    std::cout<<"Done!"<<std::endl;
}

void CCProtonPi0_CrossSection::Calc_CrossSections()
{
    Calc_Normalized_NBackground();

    // Regular Cross Section Calculation
    Calc_CrossSection(muon_P);   
    Calc_CrossSection(muon_theta);   
    Calc_CrossSection(pi0_P);   
    Calc_CrossSection(pi0_KE);   
    Calc_CrossSection(pi0_theta);   
    Calc_CrossSection(QSq);   
    Calc_CrossSection(W);   
    Calc_CrossSection(Enu);   
 
    // Before FSI Cross Section Calculation
    Calc_CrossSection_BeforeFSI(muon_P);   
    Calc_CrossSection_BeforeFSI(muon_theta);   
    Calc_CrossSection_BeforeFSI(pi0_P);   
    Calc_CrossSection_BeforeFSI(pi0_KE);   
    Calc_CrossSection_BeforeFSI(pi0_theta);   
    Calc_CrossSection_BeforeFSI(QSq);   
    Calc_CrossSection_BeforeFSI(W);   
    Calc_CrossSection_BeforeFSI(Enu);   
  
    // FSI Type Cross Section Calculation
    Calc_CrossSection_FSIType(muon_P);   
    Calc_CrossSection_FSIType(muon_theta);   
    Calc_CrossSection_FSIType(pi0_P);   
    Calc_CrossSection_FSIType(pi0_KE);   
    Calc_CrossSection_FSIType(pi0_theta);   
    Calc_CrossSection_FSIType(QSq);   
    Calc_CrossSection_FSIType(W);   
    Calc_CrossSection_FSIType(Enu);   
 
    // Int Type Cross Section Calculation
    Calc_CrossSection_IntType(muon_P);   
    Calc_CrossSection_IntType(muon_theta);   
    Calc_CrossSection_IntType(pi0_P);   
    Calc_CrossSection_IntType(pi0_KE);   
    Calc_CrossSection_IntType(pi0_theta);   
    Calc_CrossSection_IntType(QSq);   
    Calc_CrossSection_IntType(W);   
    Calc_CrossSection_IntType(Enu);   

    writeHistograms();
}

void CCProtonPi0_CrossSection::Calc_CrossSection(XSec &var)
{
    std::cout<<"\n-----------------------------------------------------------------------"<<std::endl;
    std::cout<<"Calculating Cross Section for "<<var.name<<std::endl;

    // Data Correction
    var.bckg_subtracted = Subtract_Background(var.all, var.mc_reco_bckg, var.bckg_estimated,var.name);
    var.unfolded = Unfold_Data(var.bckg_subtracted, var.response, var.name, var.nIterations);   
    var.efficiency_corrected = Efficiency_Divide(var.unfolded, var.eff, var.name);   

    // Integrate Flux
    var.flux_integrated = Integrate_Flux(var.efficiency_corrected, var.name, "_flux_integrated");
    
    // Calculate Final Cross Section
    var.xsec = Calc_FinalCrossSection(var.flux_integrated, var.name, "_xsec"); 
 
    AddLabels_XSecHist(var);
    Scale_XSecHist(var.xsec, var.smallest_bin_width);

    std::cout<<"Done!"<<std::endl;
    std::cout<<"-----------------------------------------------------------------------"<<std::endl;
}

void CCProtonPi0_CrossSection::Calc_CrossSection_BeforeFSI(XSec &var)
{
    std::cout<<"\n-----------------------------------------------------------------------"<<std::endl;
    std::cout<<"Calculating Cross Section for "<<var.name<<" Before FSI"<<std::endl;

    // Integrate Flux
    var.flux_integrated_BeforeFSI = Integrate_Flux(var.efficiency_corrected_BeforeFSI, var.name, "_flux_integrated_BeforeFSI");
    
    // Calculate Final Cross Section
    var.xsec_BeforeFSI = Calc_FinalCrossSection(var.flux_integrated_BeforeFSI, var.name, "_xsec_BeforeFSI"); 
 
    Scale_XSecHist(var.xsec_BeforeFSI, var.smallest_bin_width);

    std::cout<<"Done!"<<std::endl;
    std::cout<<"-----------------------------------------------------------------------"<<std::endl;
}

void CCProtonPi0_CrossSection::Calc_CrossSection_FSIType(XSec &var)
{
    std::cout<<"\n-----------------------------------------------------------------------"<<std::endl;

    for (int i = 0; i < nFSIType; ++i){
        std::cout<<"Calculating Cross Section for "<<var.name<<" FSI Type "<<i<<std::endl;

        // Integrate Flux
        var.flux_integrated_FSIType[i] = Integrate_Flux(var.efficiency_corrected_FSIType[i], var.name, Form("%s_%d", "_flux_integrated_FSIType",i));

        // Calculate Final Cross Section
        var.xsec_FSIType[i] = Calc_FinalCrossSection(var.flux_integrated_FSIType[i], var.name, Form("%s_%d","_xsec_FSIType",i)); 

        Scale_XSecHist(var.xsec_FSIType[i], var.smallest_bin_width);
    }
    std::cout<<"Done!"<<std::endl;
    std::cout<<"-----------------------------------------------------------------------"<<std::endl;
}

void CCProtonPi0_CrossSection::Calc_CrossSection_IntType(XSec &var)
{
    std::cout<<"\n-----------------------------------------------------------------------"<<std::endl;

    for (int i = 0; i < nIntType; ++i){
        std::cout<<"Calculating Cross Section for "<<var.name<<" Int Type "<<i<<std::endl;

        // Integrate Flux
        var.flux_integrated_IntType[i] = Integrate_Flux(var.efficiency_corrected_IntType[i], var.name, Form("%s_%d", "_flux_integrated_IntType",i));

        // Calculate Final Cross Section
        var.xsec_IntType[i] = Calc_FinalCrossSection(var.flux_integrated_IntType[i], var.name, Form("%s_%d","_xsec_IntType",i)); 

        Scale_XSecHist(var.xsec_IntType[i], var.smallest_bin_width);
    }
    std::cout<<"Done!"<<std::endl;
    std::cout<<"-----------------------------------------------------------------------"<<std::endl;
}

void CCProtonPi0_CrossSection::AddLabels_XSecHist(XSec &var)
{
    // Add Labels
    var.xsec->SetTitle(var.plot_title.c_str());
    var.xsec->GetXaxis()->SetTitle(var.plot_xlabel.c_str());
    var.xsec->GetYaxis()->SetTitle(var.plot_ylabel.c_str());
}

void CCProtonPi0_CrossSection::Scale_XSecHist(MnvH1D* xsec_hist, double smallest_bin_width)
{
    // Scale Cross Section Result to match with Label  
    xsec_hist->SetNormBinWidth(smallest_bin_width);
    double norm_bin_width = xsec_hist->GetNormBinWidth();
    double quote_width = 1.0; // per Unit (degree, GeV, GeV^2 etc...) 
    double scale = quote_width/norm_bin_width; 
    xsec_hist->Scale(scale);
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
    
    std::vector<std::string> all_err_bands;
    std::vector<std::string> bckg_err_bands;
    std::vector<int> all_hist_ind;
    std::vector<int> bckg_hist_ind;

    GetAllUniverses(invMass_all, data_invMass, all_err_bands, all_hist_ind);
    GetAllUniverses(invMass_mc_reco_bckg, bckg_invMass, bckg_err_bands, bckg_hist_ind);

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
        text_out<<"Estimating Background in Error Band: "<<all_err_bands[i]<<" Universe = "<<all_hist_ind[i]<<std::endl;

        // Subtract bckg hist from data hist
        signal_invMass[i]->Add(bckg_invMass[i],-1); 

        N_Background_Data.push_back(bckg_invMass[i]->Integral());
        
        double N_Signal_Data = signal_invMass[i]->Integral();
        double N_All_Data = data_invMass[i]->Integral();

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

void CCProtonPi0_CrossSection::NormalizeHistogram(TH1D* h)
{
    text_out<<"\tNormalizing Background Shape on TH1D"<<std::endl;
    int NBins = h->GetNbinsX();
    double area = h->Integral();
    double nOverFlow = h->GetBinContent(NBins+1);
    double nUnderFlow = h->GetBinContent(0);
    text_out<<"\t\tBefore Norm = "<<area<<std::endl;
    h->Scale(1/(area+nOverFlow+nUnderFlow));
    text_out<<"\t\tAfter Norm = "<<h->Integral()<<std::endl;
    text_out<<"\tDone!"<<std::endl;
}

void CCProtonPi0_CrossSection::NormalizeHistogram(MnvH1D* h)
{
    text_out<<"\tNormalizing Background Shape on MnvH1D"<<std::endl;
    int NBins = h->GetNbinsX();
    double area = h->Integral();
    double nOverFlow = h->GetBinContent(NBins+1);
    double nUnderFlow = h->GetBinContent(0);
    text_out<<"\t\tBefore Norm = "<<area<<std::endl;
    h->Scale(1/(area+nOverFlow+nUnderFlow),"",false); // Scale only on CentralValue
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
    
    GetPointersAllUniverses(data, data_all_universes);
    GetPointersAllUniverses(bckg_subtracted, bckg_subtracted_all_universes);
    GetPointersAllUniverses(mc_bckg, mc_bckg_all_universes);

    text_out<<"N(Universes) = "<<data_all_universes.size()<<std::endl;
   
    // Sanity Check
    if ( bckg_subtracted_all_universes.size() != mc_bckg_all_universes.size()){
        text_out<<"WARNING! - Subtract Background N(Universes) NOT Same!"<<std::endl;
        exit(1);
    }
    
    // Get Universes Names  -- For Logging purposes only
    std::vector<TH1D*> dummy_hists;
    std::vector<std::string> err_bands;
    std::vector<int> hist_ind;

    GetAllUniverses(data, dummy_hists, err_bands, hist_ind);
    ClearAllUniversesVector(dummy_hists);

    // Loop Over All Universes
    for (unsigned int i = 0; i < bckg_subtracted_all_universes.size(); ++i){
        text_out<<"\tSubtracting Background in Error Band: "<<err_bands[i+1]<<" Universe = "<<hist_ind[i+1]<<std::endl;
       
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

MnvH1D* CCProtonPi0_CrossSection::Unfold_Data(MnvH1D* bckg_subtracted, MnvH2D* response, std::string var_name, int nIter)
{
    std::cout<<"Unfolding Data for "<<var_name<<std::endl;
    // Init Histogram
    MnvH1D* unfolded = 0;

    std::cout<<"\tNumber of iterations = "<<nIter<<std::endl;

    // Use MnvUnfold to Unfold Data
    MinervaUnfold::MnvUnfold::Get().UnfoldHisto(unfolded, response, bckg_subtracted, RooUnfold::kBayes, nIter, true);

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

MnvH1D* CCProtonPi0_CrossSection::Integrate_Flux(MnvH1D* data_efficiency_corrected, std::string var_name, std::string hist_name)
{
    std::cout<<"Integrating Flux for "<<var_name<<std::endl;

    MnvH1D* flux_integrated = new MnvH1D (*data_efficiency_corrected); 
   
    hist_name = var_name + hist_name;
    flux_integrated->SetName(hist_name.c_str());
    flux_integrated->SetTitle("Flux Integrated Data");

    // Divide by total integral except Enu
    // Enu divided bin by bin
    if (var_name.compare("Enu") == 0 ){
        std::cout<<"Flux Integration for Neutrino Energy!"<<std::endl; 
        flux_integrated->Divide(data_efficiency_corrected, h_flux_rebinned);
    }else{
        // Scale All Universes with Central Value
        flux_integrated->Scale(1/cv_flux_integral);

        // Rescale Flux Universes with their own integrals
        MnvVertErrorBand* flux_err = flux_integrated->GetVertErrorBand("Flux");
        std::vector<TH1D*> flux_unv = flux_err->GetHists();
        for (unsigned int i = 0; i < flux_unv.size(); ++i){
            flux_unv[i]->Scale(cv_flux_integral);
            flux_unv[i]->Scale(1/unv_flux_integral[i]);
        }

    }
    std::cout<<"Done!"<<std::endl;

    return flux_integrated;
}

MnvH1D* CCProtonPi0_CrossSection::Calc_FinalCrossSection(MnvH1D* flux_integrated, std::string var_name, std::string hist_name)
{
    MnvH1D* h_xs = new MnvH1D(*flux_integrated);
    hist_name = var_name + hist_name;
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

    if (RemoveErrorBands){
        invMass_all->ClearAllErrorBands();
        invMass_mc_reco_signal->ClearAllErrorBands();
        invMass_mc_reco_bckg->ClearAllErrorBands();
    }
}


void CCProtonPi0_CrossSection::initFluxHistograms()
{
    // Get Reweighted Flux Histogram
    delete frw;
    frw = new FluxReweighter(14, applyNuEConstraint, FluxReweighter::minervaLE_FHC, new_flux, old_flux);
    h_flux_minervaLE_FHC = new MnvH1D (*(frw->GetFluxReweighted(14)));
    h_flux_minervaLE_FHC->SetName("h_flux_minervaLE_FHC");
    h_flux_minervaLE_FHC->Scale(1/mSq_to_cmSq); // Our measurement scale is cm2

    // Flux Universes have different integrals -- We need to get them
    IntegrateAllFluxUniverses();

    // Rebin the Flux Histogram for Neutrino Energy Divide() Operation
    RebinFluxHistogram();

    // Add Missing Error Bands after rebinning
    AddErrorBands_FluxHistogram();
}

void CCProtonPi0_CrossSection::initHistograms(XSec &var)
{
    MnvH1D* temp = NULL;
    for (int i = 0; i < nFSIType; ++i){
        var.efficiency_corrected_FSIType.push_back(temp);
        var.flux_integrated_FSIType.push_back(temp);
        var.xsec_FSIType.push_back(temp);
    }
    for (int i = 0; i < nIntType; ++i){
        var.efficiency_corrected_IntType.push_back(temp);
        var.flux_integrated_IntType.push_back(temp);
        var.xsec_IntType.push_back(temp);
    }
    
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

    // All Signal Before FSI is on Truth File
    hist_name = var.name + "_mc_truth_all_signal_BeforeFSI";
    var.efficiency_corrected_BeforeFSI = GetMnvH1D(f_truth, hist_name);
    for (int i = 0; i < nFSIType; ++i){
        hist_name = var.name + "_mc_truth_all_signal_FSIType_" + std::to_string((long long int)i);
        var.efficiency_corrected_FSIType[i] = GetMnvH1D(f_truth, hist_name);
    }
 
    for (int i = 0; i < nIntType; ++i){
        hist_name = var.name + "_mc_truth_all_signal_IntType_" + std::to_string((long long int)i);
        var.efficiency_corrected_IntType[i] = GetMnvH1D(f_truth, hist_name);
    }

    hist_name = var.name + "_response";
    var.response = GetMnvH2D(var.f_mc, hist_name);

    // For
    if (RemoveErrorBands){
        var.all->ClearAllErrorBands();
        var.mc_reco_signal->ClearAllErrorBands();
        var.mc_reco_bckg->ClearAllErrorBands();
        var.mc_truth_all_signal->ClearAllErrorBands();
        var.mc_truth_signal->ClearAllErrorBands();
        var.eff->ClearAllErrorBands();
        var.efficiency_corrected_BeforeFSI->ClearAllErrorBands();
        for (int i = 0; i < nFSIType; ++i){
            var.efficiency_corrected_FSIType[i]->ClearAllErrorBands();
        }
        for (int i = 0; i < nIntType; ++i){
            var.efficiency_corrected_IntType[i]->ClearAllErrorBands();
        }
        var.response->ClearAllErrorBands();
    }
}

double CCProtonPi0_CrossSection::GetFluxHistContent(TH1* hist, double low1, double low2)
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

void CCProtonPi0_CrossSection::AddErrorBands_FluxHistogram()
{
    // Flux Histogram have only Flux Error Bands -- Other Error Bands are added and filled with CV
    AddVertErrorBands_FluxHistogram(h_flux_rebinned);
    AddLatErrorBands_FluxHistogram(h_flux_rebinned);
}

void CCProtonPi0_CrossSection::IntegrateAllFluxUniverses()
{
    // Integrate for Central Value
    cv_flux_integral = h_flux_minervaLE_FHC->Integral(4,30,"width"); // Whole Range 1.5-20 GeV
    //cv_flux_integral = h_flux_minervaLE_FHC->Integral(4,8,"width"); // Low Enu 1.5-4.0 GeV
    //cv_flux_integral = h_flux_minervaLE_FHC->Integral(8,20,"width"); // High Enu 4.0-10 GeV
    text_out<<"Signal Region Flux Integrals"<<std::endl;    
    text_out<<"\tCentral Value = "<<cv_flux_integral<<std::endl;    
    
    // Integrate Flux Error Band Universes
    MnvVertErrorBand* flux_err_band = h_flux_minervaLE_FHC->GetVertErrorBand("Flux");
    const std::vector<TH1D*> flux_err_band_universes = flux_err_band->GetHists();

    for (unsigned int i = 0; i < flux_err_band_universes.size(); ++i){
        double temp_integral = flux_err_band_universes[i]->Integral(4,30,"width"); // Whole Range 1.5-20 GeV
        //double temp_integral = flux_err_band_universes[i]->Integral(4,8,"width"); // Low Enu 1.5-4.0 GeV
        //double temp_integral = flux_err_band_universes[i]->Integral(8,20,"width"); // High Enu 4.0-10 GeV
        unv_flux_integral.push_back(temp_integral);
        text_out<<"\tFlux Unv "<<i<<" = "<<temp_integral<<std::endl;
    }
}

void CCProtonPi0_CrossSection::RebinFluxHistogram()
{
    // Create an empty MnvH1D and Fill it with Central Value 
    h_flux_rebinned = new MnvH1D( "h_flux_rebinned","Flux (rebinned)", binList.size_Enu, binList.a_Enu);
    h_flux_rebinned->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    h_flux_rebinned->GetYaxis()->SetTitle("#nu_{#mu} s/cm^{2}/P.O.T./GeV");

    RebinFluxHistogram(h_flux_rebinned, h_flux_minervaLE_FHC);

    // ------------------------------------------------------------------------
    // Get All Universes in Flux Error Band and Rebin Them
    // ------------------------------------------------------------------------
    // 1) Get All Universes from h_flux_minervaLE_FHC as reference
    std::vector<TH1D*> reference_flux_err_unv = h_flux_minervaLE_FHC->GetVertErrorBand("Flux")->GetHists();

    // 2) Add Empty Error Band to flux_rebinned
    h_flux_rebinned->AddVertErrorBand("Flux", reference_flux_err_unv.size());

    // 3) Get All Universes from h_flux_rebinned
    std::vector<TH1D*> rebinned_flux_err_unv = h_flux_rebinned->GetVertErrorBand("Flux")->GetHists();

    // 4) Fill Empty Universes on flux_rebinned by using flux_minerva_LE_FHC
    if (rebinned_flux_err_unv.size() != reference_flux_err_unv.size()){
        RunTimeError("h_flux_rebinned have wrong N(Universes) in Flux Error Band");
    }
    for (unsigned int unv = 0; unv < reference_flux_err_unv.size(); ++unv){
        RebinFluxHistogram(rebinned_flux_err_unv[unv], reference_flux_err_unv[unv]); 
    }

    // Check Rebinned Flux Histogram Integrals -- For Testing
    text_out<<"Testing! -- Rebinned Flux Integrals"<<std::endl;
    // Integrate for Central Value
    text_out<<"Signal Region Flux Integrals"<<std::endl;    
    text_out<<"\tCentral Value = "<<h_flux_rebinned->Integral(2,12,"width")<<std::endl;    
    // Integrate Flux Error Band Universes
    MnvVertErrorBand* flux_err_band = h_flux_rebinned->GetVertErrorBand("Flux");
    const std::vector<TH1D*> flux_err_band_universes = flux_err_band->GetHists();

    for (unsigned int i = 0; i < flux_err_band_universes.size(); ++i){
        double temp_integral = flux_err_band_universes[i]->Integral(2,12,"width");
        text_out<<"\tFlux Unv "<<i<<" = "<<temp_integral<<std::endl;
    }

}

void CCProtonPi0_CrossSection::RebinFluxHistogram(TH1* rebinned, TH1* reference)
{
    int nBins = rebinned->GetNbinsX();
    for (int i = 1; i <= nBins; i++){
        double low = rebinned->GetBinLowEdge(i);
        double up = rebinned->GetBinLowEdge(i+1);
        double content = GetFluxHistContent(reference,low,up);
        double bin_width = rebinned->GetBinWidth(i);
        rebinned->SetBinContent(i,content/bin_width);
    }
}

void CCProtonPi0_CrossSection::writeHistograms(XSec &var)
{
    // Data Histograms
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

    var.efficiency_corrected_BeforeFSI->SetNormBinWidth(var.smallest_bin_width);
    var.efficiency_corrected_BeforeFSI->Write();
    
    for (int i = 0; i < nFSIType; ++i){
        var.efficiency_corrected_FSIType[i]->SetNormBinWidth(var.smallest_bin_width);
        var.efficiency_corrected_FSIType[i]->Write();
    }
    for (int i = 0; i < nIntType; ++i){
        var.efficiency_corrected_IntType[i]->SetNormBinWidth(var.smallest_bin_width);
        var.efficiency_corrected_IntType[i]->Write();
    }

    var.flux_integrated->SetNormBinWidth(var.smallest_bin_width);
    var.flux_integrated->Write();
 
    var.flux_integrated_BeforeFSI->SetNormBinWidth(var.smallest_bin_width);
    var.flux_integrated_BeforeFSI->Write();
 
    for (int i = 0; i < nFSIType; ++i){
        var.flux_integrated_FSIType[i]->SetNormBinWidth(var.smallest_bin_width);
        var.flux_integrated_FSIType[i]->Write();
    }
    for (int i = 0; i < nIntType; ++i){
        var.flux_integrated_IntType[i]->SetNormBinWidth(var.smallest_bin_width);
        var.flux_integrated_IntType[i]->Write();
    }
 
    var.xsec->Write();
    var.xsec_BeforeFSI->Write();

    for (int i = 0; i < nFSIType; ++i){
        var.xsec_FSIType[i]->Write();
    }
    for (int i = 0; i < nIntType; ++i){
        var.xsec_IntType[i]->Write();
    }

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

