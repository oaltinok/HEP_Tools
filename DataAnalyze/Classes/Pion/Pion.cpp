/*
    See Pion.h header for Class Information
*/

#include "Pion.h"

Pion::Pion()
{

    cout<<"Initializing Pion Particle"<<endl;
    
    // File Locations
    rootDir =   Folder_List::f_Root_Pion;
    plotDir =   Folder_List::f_Plot_Pion;
    
    cout<<"\tRoot File: "<<rootDir<<endl;
    cout<<"\tPlot Output Folder: "<<plotDir<<endl;
    
    // Create Root File 
    f = new TFile(rootDir.c_str(),"RECREATE");

    // Initialize Bins
    bin_P.setBin(30, 0.0, 3000.0);
    bin_KE.setBin(30, 0.0, 3000.0);

    cout<<"\tInitializing Histograms "<<endl;
    P_mc = new TH1F( "P_mc","True Pion Momentum",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    P_mc->GetXaxis()->SetTitle("True Pion Momentum MeV");
    P_mc->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_P.get_width()));
    
    P_reco = new TH1F( "P_reco","Reconstructed Pion Momentum",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    P_reco->GetXaxis()->SetTitle("Reconstructed Pion Momentum MeV");
    P_reco->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_P.get_width()));
    
    P_error = new TH1F( "P_error","Error on Pion Momentum",bin_error.get_nBins(), bin_error.get_min(), bin_error.get_max() );
    P_error->GetXaxis()->SetTitle("(True - Reco) / True");
    P_error->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_error.get_width()));
    
    P_reco_mc = new TH2F( "P_reco_mc","True vs Reconstructed Pion Momentum",
                                bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max(),
                                bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max());
    P_reco_mc->GetXaxis()->SetTitle("Reconstructed Pion Momentum MeV");
    P_reco_mc->GetYaxis()->SetTitle("True Pion Momentum MeV");
    
    KE_mc = new TH1F( "KE_mc","True Pion Kinetic Energy",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    KE_mc->GetXaxis()->SetTitle("True Pion Kinetic Energy MeV");
    KE_mc->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_P.get_width()));
    
    KE_reco = new TH1F( "KE_reco","Reconstructed Pion Kinetic Energy",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    KE_reco->GetXaxis()->SetTitle("Reconstructed Pion Kinetic Energy MeV");
    KE_reco->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_P.get_width()));
    
    KE_error = new TH1F( "KE_error","Error on Pion Kinetic Energy",bin_error.get_nBins(), bin_error.get_min(), bin_error.get_max() );
    KE_error->GetXaxis()->SetTitle("(True - Reco) / True");
    KE_error->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_error.get_width()));
    
    KE_reco_mc = new TH2F( "KE_reco_mc","True vs Reconstructed Pion Kinetic Energy",
                                bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max(),
                                bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max());
    KE_reco_mc->GetXaxis()->SetTitle("Reconstructed Pion Kinetic Energy MeV");
    KE_reco_mc->GetYaxis()->SetTitle("True Pion Kinetic Energy MeV");
    
    cout<<"Initialization Complete! "<<endl;
}

void Pion::set_kineticEnergy(bool isMC)
{
    int type = getDataType(isMC);
    
    kineticEnergy[type] = p4[type].Energy() - restMass;
}

