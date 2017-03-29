#include "CCProtonPi0_CrossSection.h"

void CCProtonPi0_CrossSection::init_muon_P()
{
    muon_P.name = "muon_P";
    muon_P.isEv = false;

    muon_P.plot_title = "Differential Cross Section for P_{#mu}";
    muon_P.plot_xlabel = "Muon Momentum (GeV/c)";
    muon_P.plot_ylabel = "d#sigma/dp_{#mu} (10^{-40} cm^{2}/nucleon/GeV/c)";

    muon_P.nIterations = 2;

    // ROOT Files
    std::string rootDir = Folder_List::rootDir_Muon_data;
    muon_P.f_data = new TFile(rootDir.c_str());
    
    rootDir = Folder_List::rootDir_Muon_mc;
    muon_P.f_mc = new TFile(rootDir.c_str());

    initHistograms(muon_P);

    muon_P.smallest_bin_width = GetSmallestBinWidth(muon_P.all);
}

void CCProtonPi0_CrossSection::init_muon_theta()
{
    muon_theta.name = "muon_theta";
    muon_theta.isEv = false;
 
    muon_theta.plot_title = "Differential Cross Section for #theta_{#mu}";
    muon_theta.plot_xlabel = "Muon Angle (deg)";
    muon_theta.plot_ylabel = "d#sigma/d#theta_{#mu} (10^{-40} cm^{2}/nucleon/deg)";

    muon_theta.nIterations = 1;
    
    // ROOT Files
    std::string rootDir = Folder_List::rootDir_Muon_data;
    muon_theta.f_data = new TFile(rootDir.c_str());
    
    rootDir = Folder_List::rootDir_Muon_mc;
    muon_theta.f_mc = new TFile(rootDir.c_str());

    initHistograms(muon_theta);
    
    muon_theta.smallest_bin_width = GetSmallestBinWidth(muon_theta.all);
}

void CCProtonPi0_CrossSection::init_pi0_P()
{
    pi0_P.name = "pi0_P";
    pi0_P.isEv = false;

    pi0_P.plot_title = "Differential Cross Section for P_{#pi^{0}}";
    pi0_P.plot_xlabel = "Pion Momentum (GeV/c)";
    pi0_P.plot_ylabel = "d#sigma/dp_{#pi^{0}} (10^{-40} cm^{2}/nucleon/GeV/c)";

    pi0_P.nIterations = 5;
    
    // ROOT Files
    std::string rootDir = Folder_List::rootDir_Pion_data;
    pi0_P.f_data = new TFile(rootDir.c_str());
    
    rootDir = Folder_List::rootDir_Pion_mc;
    pi0_P.f_mc = new TFile(rootDir.c_str());

    initHistograms(pi0_P);
    
    pi0_P.smallest_bin_width = GetSmallestBinWidth(pi0_P.all);
}

void CCProtonPi0_CrossSection::init_pi0_KE()
{
    pi0_KE.name = "pi0_KE";
    pi0_KE.isEv = false;

    pi0_KE.plot_title = "Differential Cross Section for T_{#pi^{0}}";
    pi0_KE.plot_xlabel = "Pion Kinetic Energy (GeV)";
    pi0_KE.plot_ylabel = "d#sigma/dT_{#pi^{0}} (10^{-40} cm^{2}/nucleon/GeV)";

    pi0_KE.nIterations = 4;
    
    // ROOT Files
    std::string rootDir = Folder_List::rootDir_Pion_data;
    pi0_KE.f_data = new TFile(rootDir.c_str());
    
    rootDir = Folder_List::rootDir_Pion_mc;
    pi0_KE.f_mc = new TFile(rootDir.c_str());

    initHistograms(pi0_KE);

    pi0_KE.smallest_bin_width = GetSmallestBinWidth(pi0_KE.all);
}

void CCProtonPi0_CrossSection::init_pi0_theta()
{
    pi0_theta.name = "pi0_theta";
    pi0_theta.isEv = false;

    pi0_theta.plot_title = "Differential Cross Section for #theta_{#pi^{0}}";
    pi0_theta.plot_xlabel = "Pion Angle (deg)";
    pi0_theta.plot_ylabel = "d#sigma/d#theta_{#pi^{0}} (10^{-40} cm^{2}/nucleon/deg)";

    pi0_theta.nIterations = 4;

    // ROOT Files
    std::string rootDir = Folder_List::rootDir_Pion_data;
    pi0_theta.f_data = new TFile(rootDir.c_str());
    
    rootDir = Folder_List::rootDir_Pion_mc;
    pi0_theta.f_mc = new TFile(rootDir.c_str());

    initHistograms(pi0_theta);

    pi0_theta.smallest_bin_width = GetSmallestBinWidth(pi0_theta.all);
}

void CCProtonPi0_CrossSection::init_W()
{
    W.name = "W";
    W.isEv = false;

    W.plot_title = "Differential Cross Section for W";
    W.plot_xlabel = "W (GeV)";
    W.plot_ylabel = "d#sigma/dW (10^{-40} cm^{2}/nucleon/GeV)";

    W.nIterations = 4;
    
    // ROOT Files
    std::string rootDir = Folder_List::rootDir_Interaction_data;
    W.f_data = new TFile(rootDir.c_str());
    
    rootDir = Folder_List::rootDir_Interaction_mc;
    W.f_mc = new TFile(rootDir.c_str());

    initHistograms(W);

    W.smallest_bin_width = GetSmallestBinWidth(W.all);
}

void CCProtonPi0_CrossSection::init_QSq()
{
    QSq.name = "QSq";
    QSq.isEv = false;

    QSq.plot_title = "Differential Cross Section for Q^{2}";
    QSq.plot_xlabel = "Q^{2} (GeV^{2})";
    QSq.plot_ylabel = "d#sigma/dQ^{2} (10^{-40} cm^{2}/nucleon/GeV^{2})";

    QSq.nIterations = 4;
    
    // ROOT Files
    std::string rootDir = Folder_List::rootDir_Interaction_data;
    QSq.f_data = new TFile(rootDir.c_str());
    
    rootDir = Folder_List::rootDir_Interaction_mc;
    QSq.f_mc = new TFile(rootDir.c_str());

    initHistograms(QSq);

    QSq.smallest_bin_width = GetSmallestBinWidth(QSq.all);
}

void CCProtonPi0_CrossSection::init_Enu()
{
    Enu.name = "Enu";
    Enu.isEv = true;

    Enu.plot_title = "Total Cross Section for E_{#nu}";
    Enu.plot_xlabel = "Neutrino Energy (GeV)";
    Enu.plot_ylabel = "#sigma(E_{#nu}) (10^{-40} cm^{2}/nucleon)";

    Enu.nIterations = 4;

    // ROOT Files
    std::string rootDir = Folder_List::rootDir_Interaction_data;
    Enu.f_data = new TFile(rootDir.c_str());
    
    rootDir = Folder_List::rootDir_Interaction_mc;
    Enu.f_mc = new TFile(rootDir.c_str());

    initHistograms(Enu);

    Enu.smallest_bin_width = GetSmallestBinWidth(Enu.all);
}


