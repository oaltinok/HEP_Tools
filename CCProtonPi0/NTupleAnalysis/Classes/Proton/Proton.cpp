/*
    See Proton.h header for Class Information
*/

#ifndef Proton_cpp
#define Proton_cpp
#include "Proton.h"

using namespace std;

Proton::Proton()
{
    cout<<"Wrong usage of Proton! Must include Analysis Mode"<<endl;
    exit(EXIT_FAILURE);
}

Proton::Proton(int nMode)
{
    cout<<"Initializing Proton Particle"<<endl;
    
    setAnalysisMode(nMode);
    
    // File Locations
    rootDir =   Folder_List::output + Folder_List::rootOut + branchDir + "Proton.root";
    
    cout<<"\tRoot File: "<<rootDir<<endl;
    
    // Create Root File 
    f = new TFile(rootDir.c_str(),"RECREATE");

    // Initialize Bins
    bin_E.setBin(30,0.0,3000.0);
    bin_P.setBin(20, 0.0, 2000.0);
    bin_KE.setBin(20, 0.0, 2000.0);

    cout<<"\tInitializing Histograms "<<endl;
       
    // Default Histograms
    partScore = new TH1D( "partScore","Proton Particle Score",bin_partScore.get_nBins(), bin_partScore.get_min(), bin_partScore.get_max() );
    partScore->GetXaxis()->SetTitle("Particle Score");
    partScore->GetYaxis()->SetTitle(Form("Number of Protons / %3.1f ",bin_partScore.get_width()));
    
    E_mc = new TH1D( "E_mc","True Proton Energy",bin_E.get_nBins(), bin_E.get_min(), bin_E.get_max() );
    E_mc->GetXaxis()->SetTitle("True Proton Energy [MeV]");
    E_mc->GetYaxis()->SetTitle(Form("Number of Protons / %3.1f [MeV] ",bin_E.get_width()));
    
    E_reco = new TH1D( "E_reco","Reconstructed Proton Energy",bin_E.get_nBins(), bin_E.get_min(), bin_E.get_max() );
    E_reco->GetXaxis()->SetTitle("Reconstructed Proton Energy [MeV]");
    E_reco->GetYaxis()->SetTitle(Form("Number of Protons / %3.1f [MeV] ",bin_E.get_width()));
    
    E_error = new TH1D( "E_error","Error on Proton Energy",bin_error.get_nBins(), bin_error.get_min(), bin_error.get_max() );
    E_error->GetXaxis()->SetTitle("(Reco- True) / True");
    E_error->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",bin_error.get_width()));
    
    E_reco_mc = new TH2D( "E_reco_mc","True vs Reconstructed Proton Energy",
                          bin_E.get_nBins(), bin_E.get_min(), bin_E.get_max(),
                          bin_E.get_nBins(), bin_E.get_min(), bin_E.get_max());
                          E_reco_mc->GetXaxis()->SetTitle("Reconstructed Proton Energy [MeV]");
                          E_reco_mc->GetYaxis()->SetTitle("True Proton Energy [MeV]");
       
    P_mc = new TH1D( "P_mc","True Proton Momentum",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    P_mc->GetXaxis()->SetTitle("True Proton Momentum [MeV]");
    P_mc->GetYaxis()->SetTitle(Form("Number of Protons / %3.1f [MeV] ",bin_P.get_width()));
    
    P_reco = new TH1D( "P_reco","Reconstructed Proton Momentum",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    P_reco->GetXaxis()->SetTitle("Reconstructed Proton Momentum [MeV]");
    P_reco->GetYaxis()->SetTitle(Form("Number of Protons / %3.1f [MeV] ",bin_P.get_width()));
    
    P_error = new TH1D( "P_error","Error on Proton Momentum",bin_error.get_nBins(), bin_error.get_min(), bin_error.get_max() );
    P_error->GetXaxis()->SetTitle("(Reco- True) / True");
    P_error->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",bin_error.get_width()));
    
    P_reco_mc = new TH2D( "P_reco_mc","True vs Reconstructed Proton Momentum",
                                bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max(),
                                bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max());
    P_reco_mc->GetXaxis()->SetTitle("Reconstructed Proton Momentum [MeV]");
    P_reco_mc->GetYaxis()->SetTitle("True Proton Momentum [MeV]");
    
    KE_mc = new TH1D( "KE_mc","True Proton Kinetic Energy",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    KE_mc->GetXaxis()->SetTitle("True Proton Kinetic Energy [MeV]");
    KE_mc->GetYaxis()->SetTitle(Form("Number of Protons / %3.1f [MeV] ",bin_P.get_width()));
    
    KE_reco = new TH1D( "KE_reco","Reconstructed Proton Kinetic Energy",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    KE_reco->GetXaxis()->SetTitle("Reconstructed Proton Kinetic Energy [MeV]");
    KE_reco->GetYaxis()->SetTitle(Form("Number of Protons / %3.1f [MeV] ",bin_P.get_width()));
    
    KE_error = new TH1D( "KE_error","Error on Proton Kinetic Energy",bin_error.get_nBins(), bin_error.get_min(), bin_error.get_max() );
    KE_error->GetXaxis()->SetTitle("(Reco- True) / True");
    KE_error->GetYaxis()->SetTitle(Form("Candidates / %3.2f",bin_error.get_width()));
    
    KE_reco_mc = new TH2D( "KE_reco_mc","True vs Reconstructed Proton Kinetic Energy",
                                bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max(),
                                bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max());
    KE_reco_mc->GetXaxis()->SetTitle("Reconstructed Proton Kinetic Energy [MeV]");
    KE_reco_mc->GetYaxis()->SetTitle("True Proton Kinetic Energy [MeV]");
    
    angleMuon_mc = new TH1D( "angleMuon_mc","True Proton Angle wrt. Muon",bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max() );
    angleMuon_mc->GetXaxis()->SetTitle("True Proton Angle wrt. Muon [Degree]");
    angleMuon_mc->GetYaxis()->SetTitle(Form("Number of Protons / %3.1f [Degree]",bin_angle.get_width()));
    
    angleMuon_reco = new TH1D( "angleMuon_reco","Reconstructed Proton Angle wrt. Muon",bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max() );
    angleMuon_reco->GetXaxis()->SetTitle("Reconstructed Proton Angle wrt. Muon [Degree]");
    angleMuon_reco->GetYaxis()->SetTitle(Form("Number of Protons / %3.1f [Degree]",bin_angle.get_width()));
    
    angleMuon_error = new TH1D( "angleMuon_error","Error on Proton Angle wrt. Muon",bin_error.get_nBins(), bin_error.get_min(), bin_error.get_max() );
    angleMuon_error->GetXaxis()->SetTitle("(Reco- True) / True");
    angleMuon_error->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",bin_error.get_width()));
    
    angleMuon_reco_mc = new TH2D( "angleMuon_reco_mc","True vs Reconstructed Proton Angle wrt. Muon",
                                bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max(),
                                bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max());
    angleMuon_reco_mc->GetXaxis()->SetTitle("Reconstructed Proton Angle wrt. Muon [Degree]");
    angleMuon_reco_mc->GetYaxis()->SetTitle("True Proton Angle wrt. Muon [Degree]");
    
    angleBeam_mc = new TH1D( "angleBeam_mc","True Proton Angle wrt. Beam",bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max() );
    angleBeam_mc->GetXaxis()->SetTitle("True Proton Angle wrt. Beam [Degree]");
    angleBeam_mc->GetYaxis()->SetTitle(Form("Number of Protons / %3.1f [Degree]",bin_angle.get_width()));
    
    angleBeam_reco = new TH1D( "angleBeam_reco","Reconstructed Proton Angle wrt. Beam",bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max() );
    angleBeam_reco->GetXaxis()->SetTitle("Reconstructed Proton Angle wrt. Beam [Degree]");
    angleBeam_reco->GetYaxis()->SetTitle(Form("Number of Protons / %3.1f [Degree]",bin_angle.get_width()));
    
    angleBeam_error = new TH1D( "angleBeam_error","Error on Proton Angle wrt. Beam",bin_error.get_nBins(), bin_error.get_min(), bin_error.get_max() );
    angleBeam_error->GetXaxis()->SetTitle("(Reco- True) / True");
    angleBeam_error->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",bin_error.get_width()));
    
    angleBeam_reco_mc = new TH2D( "angleBeam_reco_mc","True vs Reconstructed Proton Angle wrt. Beam",
                                bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max(),
                                bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max());
    angleBeam_reco_mc->GetXaxis()->SetTitle("Reconstructed Proton Angle wrt. Beam [Degree]");
    angleBeam_reco_mc->GetYaxis()->SetTitle("True Proton Angle wrt. Beam [Degree]");
    
    cout<<"Done!"<<endl;
}

void Proton::set_kineticEnergy(bool isMC)
{
    int type = getDataType(isMC);
    
    kineticEnergy[type] = p4[type].Energy() - restMass;
}

#endif


