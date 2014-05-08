/*
    See Proton.h header for Class Information
*/

#ifndef Proton_cpp
#define Proton_cpp
#include "Proton.h"

using namespace std;

Proton::Proton()
{
    cout<<"Initializing Proton Particle"<<endl;
    
    // File Locations
    rootDir =   "Output/RootFiles/Proton.root";
    plotDir =   "Output/Plots/Proton/";
    
    cout<<"\tRoot File: "<<rootDir<<endl;
    cout<<"\tPlot Output Folder: "<<plotDir<<endl;
    
    // Create Root File 
    f = new TFile(rootDir.c_str(),"RECREATE");

    // Initialize Bins
    bin_P.setBin(40, 0.0, 4000.0);
    bin_KE.setBin(40, 0.0, 4000.0);

    cout<<"\tInitializing Histograms "<<endl;
    partScore = new TH1F( "partScore","Proton Particle Score",bin_partScore.get_nBins(), bin_partScore.get_min(), bin_partScore.get_max() );
    partScore->GetXaxis()->SetTitle("Particle Score");
    partScore->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_P.get_width()));
    
    P_mc = new TH1F( "P_mc","True Proton Momentum",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    P_mc->GetXaxis()->SetTitle("True Proton Momentum MeV");
    P_mc->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_P.get_width()));
    
    P_reco = new TH1F( "P_reco","Reconstructed Proton Momentum",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    P_reco->GetXaxis()->SetTitle("Reconstructed Proton Momentum MeV");
    P_reco->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_P.get_width()));
    
    P_error = new TH1F( "P_error","Error on Proton Momentum",bin_error.get_nBins(), bin_error.get_min(), bin_error.get_max() );
    P_error->GetXaxis()->SetTitle("(Reco- True) / True");
    P_error->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_error.get_width()));
    
    P_reco_mc = new TH2F( "P_reco_mc","True vs Reconstructed Proton Momentum",
                                bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max(),
                                bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max());
    P_reco_mc->GetXaxis()->SetTitle("Reconstructed Proton Momentum MeV");
    P_reco_mc->GetYaxis()->SetTitle("True Proton Momentum MeV");
    
    KE_mc = new TH1F( "KE_mc","True Proton Kinetic Energy",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    KE_mc->GetXaxis()->SetTitle("True Proton Kinetic Energy MeV");
    KE_mc->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_P.get_width()));
    
    KE_reco = new TH1F( "KE_reco","Reconstructed Proton Kinetic Energy",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    KE_reco->GetXaxis()->SetTitle("Reconstructed Proton Kinetic Energy MeV");
    KE_reco->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_P.get_width()));
    
    KE_error = new TH1F( "KE_error","Error on Proton Kinetic Energy",bin_error.get_nBins(), bin_error.get_min(), bin_error.get_max() );
    KE_error->GetXaxis()->SetTitle("(Reco- True) / True");
    KE_error->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_error.get_width()));
    
    KE_reco_mc = new TH2F( "KE_reco_mc","True vs Reconstructed Proton Kinetic Energy",
                                bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max(),
                                bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max());
    KE_reco_mc->GetXaxis()->SetTitle("Reconstructed Proton Kinetic Energy MeV");
    KE_reco_mc->GetYaxis()->SetTitle("True Proton Kinetic Energy MeV");
    
    angleMuon_mc = new TH1F( "angleMuon_mc","True Proton Angle wrt. Muon",bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max() );
    angleMuon_mc->GetXaxis()->SetTitle("True Proton Angle wrt. Muon");
    angleMuon_mc->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_angle.get_width()));
    
    angleMuon_reco = new TH1F( "angleMuon_reco","Reconstructed Proton Angle wrt. Muon",bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max() );
    angleMuon_reco->GetXaxis()->SetTitle("Reconstructed Proton Angle wrt. Muon");
    angleMuon_reco->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_angle.get_width()));
    
    angleMuon_error = new TH1F( "angleMuon_error","Error on Proton Angle wrt. Muon",bin_error.get_nBins(), bin_error.get_min(), bin_error.get_max() );
    angleMuon_error->GetXaxis()->SetTitle("(Reco- True) / True");
    angleMuon_error->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_error.get_width()));
    
    angleMuon_reco_mc = new TH2F( "angleMuon_reco_mc","True vs Reconstructed Proton Angle wrt. Muon",
                                bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max(),
                                bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max());
    angleMuon_reco_mc->GetXaxis()->SetTitle("Reconstructed Proton Angle wrt. Muon");
    angleMuon_reco_mc->GetYaxis()->SetTitle("True Proton Angle wrt. Muon");
    
    angleBeam_mc = new TH1F( "angleBeam_mc","True Proton Angle wrt. Beam",bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max() );
    angleBeam_mc->GetXaxis()->SetTitle("True Proton Angle wrt. Beam");
    angleBeam_mc->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_angle.get_width()));
    
    angleBeam_reco = new TH1F( "angleBeam_reco","Reconstructed Proton Angle wrt. Beam",bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max() );
    angleBeam_reco->GetXaxis()->SetTitle("Reconstructed Proton Angle wrt. Beam");
    angleBeam_reco->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_angle.get_width()));
    
    angleBeam_error = new TH1F( "angleBeam_error","Error on Proton Angle wrt. Beam",bin_error.get_nBins(), bin_error.get_min(), bin_error.get_max() );
    angleBeam_error->GetXaxis()->SetTitle("(Reco- True) / True");
    angleBeam_error->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_error.get_width()));
    
    angleBeam_reco_mc = new TH2F( "angleBeam_reco_mc","True vs Reconstructed Proton Angle wrt. Beam",
                                bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max(),
                                bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max());
    angleBeam_reco_mc->GetXaxis()->SetTitle("Reconstructed Proton Angle wrt. Beam");
    angleBeam_reco_mc->GetYaxis()->SetTitle("True Proton Angle wrt. Beam");
    
    cout<<"Initialization Complete! "<<endl;
}

void Proton::set_kineticEnergy(bool isMC)
{
    int type = getDataType(isMC);
    
    kineticEnergy[type] = p4[type].Energy() - restMass;
}

#endif


