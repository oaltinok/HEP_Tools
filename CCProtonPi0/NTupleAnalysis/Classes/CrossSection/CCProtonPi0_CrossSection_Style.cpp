#include "CCProtonPi0_CrossSection.h"

void CCProtonPi0_CrossSection::init_muon_P()
{
    muon_P.name = "muon_P";
    muon_P.isEv = false;

    muon_P.plot_title = "Differential Cross Section for P_{#mu}";
    muon_P.plot_xlabel = "Muon Momentum [GeV]";
    muon_P.plot_ylabel = "d#sigma/dp_{#mu} (10^{-40} cm^{2}/nucleon/GeV)";

    // Scale Cross Section Result to match with Label
    muon_P.bin_width = 0.5; // GeV
    muon_P.quote_width = 1.0; // per GeV
    muon_P.scale = muon_P.quote_width / muon_P.bin_width; 

    // ROOT Files
    std::string rootDir = Folder_List::rootDir_Muon_data;
    muon_P.f_data = new TFile(rootDir.c_str());
    
    rootDir = Folder_List::rootDir_Muon_mc;
    muon_P.f_mc = new TFile(rootDir.c_str());

    initHistograms(muon_P);
}

void CCProtonPi0_CrossSection::init_muon_theta()
{
    muon_theta.name = "muon_theta";
    muon_theta.isEv = false;
 
    muon_theta.plot_title = "Differential Cross Section for #theta_{#mu}";
    muon_theta.plot_xlabel = "Muon Angle [degree]";
    muon_theta.plot_ylabel = "d#sigma/d#theta_{#mu} (10^{-40} cm^{2}/nucleon/degree)";

    // Scale Cross Section Result to match with Label
    muon_theta.bin_width = 2.08; // degree
    muon_theta.quote_width = 1.0; // per GeV
    muon_theta.scale = muon_theta.quote_width / muon_theta.bin_width; 

    // ROOT Files
    std::string rootDir = Folder_List::rootDir_Muon_data;
    muon_theta.f_data = new TFile(rootDir.c_str());
    
    rootDir = Folder_List::rootDir_Muon_mc;
    muon_theta.f_mc = new TFile(rootDir.c_str());

    initHistograms(muon_theta);
}

void CCProtonPi0_CrossSection::init_pi0_P()
{
    pi0_P.name = "pi0_P";
    pi0_P.isEv = false;

    pi0_P.plot_title = "Differential Cross Section for P_{#pi^{0}}";
    pi0_P.plot_xlabel = "Pion Momentum [GeV]";
    pi0_P.plot_ylabel = "d#sigma/dp_{#pi^{0}} (10^{-40} cm^{2}/nucleon/GeV)";

    // Scale Cross Section Result to match with Label
    pi0_P.bin_width = 0.1; // GeV
    pi0_P.quote_width = 1.0; // per GeV
    pi0_P.scale = pi0_P.quote_width / pi0_P.bin_width; 

    // ROOT Files
    std::string rootDir = Folder_List::rootDir_Pion_data;
    pi0_P.f_data = new TFile(rootDir.c_str());
    
    rootDir = Folder_List::rootDir_Pion_mc;
    pi0_P.f_mc = new TFile(rootDir.c_str());

    initHistograms(pi0_P);
}

void CCProtonPi0_CrossSection::init_pi0_KE()
{
    pi0_KE.name = "pi0_KE";
    pi0_KE.isEv = false;

    pi0_KE.plot_title = "Differential Cross Section for T_{#pi^{0}}";
    pi0_KE.plot_xlabel = "Pion Kinetic Energy [GeV]";
    pi0_KE.plot_ylabel = "d#sigma/dT_{#pi^{0}} (10^{-40} cm^{2}/nucleon/GeV)";

    // Scale Cross Section Result to match with Label
    pi0_KE.bin_width = 0.1; // GeV
    pi0_KE.quote_width = 1.0; // per GeV
    pi0_KE.scale = pi0_KE.quote_width / pi0_KE.bin_width; 

    // ROOT Files
    std::string rootDir = Folder_List::rootDir_Pion_data;
    pi0_KE.f_data = new TFile(rootDir.c_str());
    
    rootDir = Folder_List::rootDir_Pion_mc;
    pi0_KE.f_mc = new TFile(rootDir.c_str());

    initHistograms(pi0_KE);
}

void CCProtonPi0_CrossSection::init_pi0_theta()
{
    pi0_theta.name = "pi0_theta";
    pi0_theta.isEv = false;

    pi0_theta.plot_title = "Differential Cross Section for #theta_{#pi^{0}}";
    pi0_theta.plot_xlabel = "Pion theta [GeV]";
    pi0_theta.plot_ylabel = "d#sigma/d#theta_{#pi^{0}} (10^{-40} cm^{2}/nucleon/degree)";

    // Scale Cross Section Result to match with Label
    pi0_theta.bin_width = 10; // degree 
    pi0_theta.quote_width = 1.0; // per degree
    pi0_theta.scale = pi0_theta.quote_width / pi0_theta.bin_width; 

    // ROOT Files
    std::string rootDir = Folder_List::rootDir_Pion_data;
    pi0_theta.f_data = new TFile(rootDir.c_str());
    
    rootDir = Folder_List::rootDir_Pion_mc;
    pi0_theta.f_mc = new TFile(rootDir.c_str());

    initHistograms(pi0_theta);
}

void CCProtonPi0_CrossSection::init_QSq()
{
    QSq.name = "QSq";
    QSq.isEv = false;

    QSq.plot_title = "Differential Cross Section for Q^{2}";
    QSq.plot_xlabel = "Q^{2} [GeV^{2}]";
    QSq.plot_ylabel = "d#sigma/dQ^{2} (10^{-40} cm^{2}/nucleon/GeV^{2})";

    // Scale Cross Section Result to match with Label
    QSq.bin_width = 0.1; // GeV^2
    QSq.quote_width = 1.0; // per GeV^2
    QSq.scale = QSq.quote_width / QSq.bin_width; 

    // ROOT Files
    std::string rootDir = Folder_List::rootDir_Interaction_data;
    QSq.f_data = new TFile(rootDir.c_str());
    
    rootDir = Folder_List::rootDir_Interaction_mc;
    QSq.f_mc = new TFile(rootDir.c_str());

    initHistograms(QSq);
}


