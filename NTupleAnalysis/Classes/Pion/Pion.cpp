/*
    See Pion.h header for Class Information
*/
#ifndef Pion_cpp
#define Pion_cpp

#include "Pion.h"

using namespace std;

Pion::Pion()
{

    cout<<"Initializing Pion Particle"<<endl;
    
    // File Locations
    rootDir =   "Output/RootFiles/Pion.root";
    plotDir =   "Output/Plots/Pion/";
    
    cout<<"\tRoot File: "<<rootDir<<endl;
    cout<<"\tPlot Output Folder: "<<plotDir<<endl;
    
    // Create Root File 
    f = new TFile(rootDir.c_str(),"RECREATE");

    // Initialize Bins
    bin_P.setBin(30, 0.0, 3000.0);
    bin_KE.setBin(30, 0.0, 3000.0);
    bin_invMass.setBin(50,0.0,500.0);

    cout<<"\tInitializing Histograms "<<endl;
    //--------------------------------------------------------------------------
    // Invariant Mass
    //--------------------------------------------------------------------------
    invMass = new TH1D( "invMass","Reconstructed Pi0 Invariant Mass",bin_invMass.get_nBins(), bin_invMass.get_min(), bin_invMass.get_max() );
    invMass->GetXaxis()->SetTitle("Reconstructed Pi0 Invariant Mass [MeV]");
    invMass->GetYaxis()->SetTitle(Form("Candidates / %3.2f [MeV]",bin_invMass.get_width()));

    invMass_0Pi0 = new TH1D( "invMass_0Pi0","True N(Pi0) = 0",bin_invMass.get_nBins(), bin_invMass.get_min(), bin_invMass.get_max() );
    invMass_0Pi0->GetXaxis()->SetTitle("Reconstructed Pion Invariant Mass [MeV]");
    invMass_0Pi0->GetYaxis()->SetTitle(Form("Number of Pions / %3.1f [MeV]",bin_invMass.get_width()));
    
    invMass_1Pi0 = new TH1D( "invMass_1Pi0","True N(Pi0) = 1",bin_invMass.get_nBins(), bin_invMass.get_min(), bin_invMass.get_max() );
    invMass_1Pi0->GetXaxis()->SetTitle("Reconstructed Pion Invariant Mass [MeV]");
    invMass_1Pi0->GetYaxis()->SetTitle(Form("Number of Pions / %3.1f [MeV]",bin_invMass.get_width()));
    
    invMass_MultPi0 = new TH1D( "invMass_MultPi0","True N(Pi0) = 2+",bin_invMass.get_nBins(), bin_invMass.get_min(), bin_invMass.get_max() );
    invMass_MultPi0->GetXaxis()->SetTitle("Reconstructed Pion Invariant Mass [MeV]");
    invMass_MultPi0->GetYaxis()->SetTitle(Form("Number of Pions / %3.1f [MeV]",bin_invMass.get_width()));
    
    //--------------------------------------------------------------------------
    // Momentum
    //--------------------------------------------------------------------------
    P_mc = new TH1D( "P_mc","True Pion Momentum",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    P_mc->GetXaxis()->SetTitle("True Pion Momentum [MeV]");
    P_mc->GetYaxis()->SetTitle(Form("Number of Pions / %3.1f [MeV]",bin_P.get_width()));
    
    P_reco = new TH1D( "P_reco","Reconstructed Pion Momentum",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    P_reco->GetXaxis()->SetTitle("Reconstructed Pion Momentum [MeV]");
    P_reco->GetYaxis()->SetTitle(Form("Number of Pions / %3.1f [MeV]",bin_P.get_width()));
    
    P_error = new TH1D( "P_error","Error on Pion Momentum",bin_error.get_nBins(), bin_error.get_min(), bin_error.get_max() );
    P_error->GetXaxis()->SetTitle("(Reco- True) / True");
    P_error->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",bin_error.get_width()));
    
    P_reco_mc = new TH2D( "P_reco_mc","True vs Reconstructed Pion Momentum",
                                bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max(),
                                bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max());
    P_reco_mc->GetXaxis()->SetTitle("Reconstructed Pion Momentum [MeV]");
    P_reco_mc->GetYaxis()->SetTitle("True Pion Momentum [MeV]");
    
    // Momentum Studies
    P_reco_0Pi0 = new TH1D( "P_reco_0Pi0","True N(Pi0) = 0",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    P_reco_0Pi0->GetXaxis()->SetTitle("Reconstructed Pion Momentum [MeV]");
    P_reco_0Pi0->GetYaxis()->SetTitle(Form("Number of Pions / %3.1f [MeV]",bin_P.get_width()));
    
    P_reco_1Pi0 = new TH1D( "P_reco_1Pi0","True N(Pi0) = 1",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    P_reco_1Pi0->GetXaxis()->SetTitle("Reconstructed Pion Momentum [MeV]");
    P_reco_1Pi0->GetYaxis()->SetTitle(Form("Number of Pions / %3.1f [MeV]",bin_P.get_width()));
    
    P_reco_MultPi0 = new TH1D( "P_reco_MultPi0","True N(Pi0) = 2+",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    P_reco_MultPi0->GetXaxis()->SetTitle("Reconstructed Pion Momentum [MeV]");
    P_reco_MultPi0->GetYaxis()->SetTitle(Form("Number of Pions / %3.1f [MeV]",bin_P.get_width()));
    
    P_reco_mc_1Pi0 = new TH2D( "P_reco_mc_1Pi0","True N(Pi0) = 1",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max(),bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    P_reco_mc_1Pi0->GetXaxis()->SetTitle("Reconstructed Pion Momentum [MeV]");
    P_reco_mc_1Pi0->GetYaxis()->SetTitle("True Pion Momentum [MeV]");
    
    P_error_1Pi0 = new TH1D( "P_error_1Pi0","Error on Pion Momentum True N(Pi0) = 1",bin_error.get_nBins(), bin_error.get_min(), bin_error.get_max() );
    P_error_1Pi0->GetXaxis()->SetTitle("(Reco- True) / True");
    P_error_1Pi0->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",bin_error.get_width()));
    
    P_reco_mc_MultPi0 = new TH2D( "P_reco_mc_MultPi0","True N(Pi0) = 2+",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max(),bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    P_reco_mc_MultPi0->GetXaxis()->SetTitle("Reconstructed Pion Momentum [MeV]");
    P_reco_mc_MultPi0->GetYaxis()->SetTitle("True Pion Momentum [MeV]");
    
    P_error_MultPi0 = new TH1D( "P_error_MultPi0","Error on Pion Momentum True N(Pi0) = 2+",bin_error.get_nBins(), bin_error.get_min(), bin_error.get_max() );
    P_error_MultPi0->GetXaxis()->SetTitle("(Reco- True) / True");
    P_error_MultPi0->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",bin_error.get_width()));
    
    //--------------------------------------------------------------------------
    // Kinetic Energy
    //--------------------------------------------------------------------------
    KE_mc = new TH1D( "KE_mc","True Pion Kinetic Energy",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    KE_mc->GetXaxis()->SetTitle("True Pion Kinetic Energy [MeV]");
    KE_mc->GetYaxis()->SetTitle(Form("Number of Pions / %3.1f [MeV]",bin_P.get_width()));
    
    KE_reco = new TH1D( "KE_reco","Reconstructed Pion Kinetic Energy",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    KE_reco->GetXaxis()->SetTitle("Reconstructed Pion Kinetic Energy [MeV]");
    KE_reco->GetYaxis()->SetTitle(Form("Number of Pions / %3.1f [MeV]",bin_P.get_width()));
    
    KE_error = new TH1D( "KE_error","Error on Pion Kinetic Energy",bin_error.get_nBins(), bin_error.get_min(), bin_error.get_max() );
    KE_error->GetXaxis()->SetTitle("(Reco- True) / True");
    KE_error->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",bin_error.get_width()));
    
    KE_reco_mc = new TH2D( "KE_reco_mc","True vs Reconstructed Pion Kinetic Energy",
                                bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max(),
                                bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max());
    KE_reco_mc->GetXaxis()->SetTitle("Reconstructed Pion Kinetic Energy [MeV]");
    KE_reco_mc->GetYaxis()->SetTitle("True Pion Kinetic Energy [MeV]");
    
    //--------------------------------------------------------------------------
    // Angle wrt Beam
    //--------------------------------------------------------------------------
    angleBeam_mc = new TH1D( "angleBeam_mc","True Pion Angle wrt. Beam",bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max() );
    angleBeam_mc->GetXaxis()->SetTitle("True Pion Angle wrt. Beam [Degree]");
    angleBeam_mc->GetYaxis()->SetTitle(Form("Number of Pions / %3.1f [Degree]",bin_angle.get_width()));
    
    angleBeam_reco = new TH1D( "angleBeam_reco","Reconstructed Pion Angle wrt. Beam",bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max() );
    angleBeam_reco->GetXaxis()->SetTitle("Reconstructed Pion Angle wrt. Beam [Degree]");
    angleBeam_reco->GetYaxis()->SetTitle(Form("Number of Pions / %3.1f [Degree]",bin_angle.get_width()));
    
    angleBeam_error = new TH1D( "angleBeam_error","Error on Pion Angle wrt. Beam",bin_error.get_nBins(), bin_error.get_min(), bin_error.get_max() );
    angleBeam_error->GetXaxis()->SetTitle("(Reco- True) / True");
    angleBeam_error->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",bin_error.get_width()));
    
    angleBeam_reco_mc = new TH2D( "angleBeam_reco_mc","True vs Reconstructed Pion Angle wrt. Beam",
                                bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max(),
                                bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max());
    angleBeam_reco_mc->GetXaxis()->SetTitle("Reconstructed Pion Angle wrt. Beam [Degree]");
    angleBeam_reco_mc->GetYaxis()->SetTitle("True Pion Angle wrt. Beam [Degree]");
    
    //--------------------------------------------------------------------------
    // Angle wrt Muon
    //--------------------------------------------------------------------------
    angleMuon_mc = new TH1D( "angleMuon_mc","True Pion Angle wrt. Muon",bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max() );
    angleMuon_mc->GetXaxis()->SetTitle("True Pion Angle wrt. Muon [Degree]");
    angleMuon_mc->GetYaxis()->SetTitle(Form("Number of Pions / %3.1f [Degree]",bin_angle.get_width()));
    
    angleMuon_reco = new TH1D( "angleMuon_reco","Reconstructed Pion Angle wrt. Muon",bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max() );
    angleMuon_reco->GetXaxis()->SetTitle("Reconstructed Pion Angle wrt. Muon [Degree]");
    angleMuon_reco->GetYaxis()->SetTitle(Form("Number of Pions / %3.1f [Degree]",bin_angle.get_width()));
    
    angleMuon_error = new TH1D( "angleMuon_error","Error on Pion Angle wrt. Muon",bin_error.get_nBins(), bin_error.get_min(), bin_error.get_max() );
    angleMuon_error->GetXaxis()->SetTitle("(Reco- True) / True");
    angleMuon_error->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",bin_error.get_width()));
    
    angleMuon_reco_mc = new TH2D( "angleMuon_reco_mc","True vs Reconstructed Pion Angle wrt. Muon",
                                bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max(),
                                bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max());
    angleMuon_reco_mc->GetXaxis()->SetTitle("Reconstructed Pion Angle wrt. Muon [Degree]");
    angleMuon_reco_mc->GetYaxis()->SetTitle("True Pion Angle wrt. Muon ");
    
    //--------------------------------------------------------------------------
    // Particle Score
    //--------------------------------------------------------------------------
    partScore = new TH1D( "partScore","Pion Particle Score",bin_partScore.get_nBins(), bin_partScore.get_min(), bin_partScore.get_max() );
    partScore->GetXaxis()->SetTitle("Particle Score");
    partScore->GetYaxis()->SetTitle(Form("Number of Pions / %3.1f ",bin_P.get_width()));
    
    
    cout<<"Initialization Complete! "<<endl;
}

void Pion::set_kineticEnergy(bool isMC)
{
    int type = getDataType(isMC);
    
    kineticEnergy[type] = p4[type].Energy() - restMass;
}

#endif
