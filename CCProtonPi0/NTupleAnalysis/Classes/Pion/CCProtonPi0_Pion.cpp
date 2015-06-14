/*
    See CCProtonPi0_Pion.h header for Class Information
*/
#ifndef CCProtonPi0_Pion_cpp
#define CCProtonPi0_Pion_cpp

#include "CCProtonPi0_Pion.h"

using namespace std;

CCProtonPi0_Pion::CCProtonPi0_Pion(int nMode) : CCProtonPi0_Particle(nMode)
{
    cout<<"Initializing CCProtonPi0_Pion"<<endl;    
    
    // File Locations
    rootDir =   Folder_List::output + Folder_List::rootOut + branchDir + "Pion.root";
    
    cout<<"\tRoot File: "<<rootDir<<endl;
    
    // Create Root File 
    f = new TFile(rootDir.c_str(),"RECREATE");

    // Initialize Bins
    bin_P.setBin(17, 0.0, 1700.0);
    bin_KE.setBin(30, 0.0, 3000.0);
    bin_invMass.setBin(60,0.0,600.0);
    bin_photonConvLength.setBin(50,0.0,100.0);
    bin_blob_energy.setBin(100,0.0,1000);
    
    initHistograms();   
    
    cout<<"Done!"<<endl;
}

void CCProtonPi0_Pion::initHistograms()
{
    gamma1_ConvLength = new TH1D( "gamma1_ConvLength","Leading Photon Conversion Length",bin_photonConvLength.get_nBins(), bin_photonConvLength.get_min(), bin_photonConvLength.get_max() );
    gamma1_ConvLength->GetXaxis()->SetTitle("Leading Photon Distance from Vertex [cm]");
    gamma1_ConvLength->GetYaxis()->SetTitle("N(Events)");
    
    gamma2_ConvLength = new TH1D( "gamma2_ConvLength","Secondary Photon Conversion Length",bin_photonConvLength.get_nBins(), bin_photonConvLength.get_min(), bin_photonConvLength.get_max() );
    gamma2_ConvLength->GetXaxis()->SetTitle("Second Photon Distance from Vertex [cm]");
    gamma2_ConvLength->GetYaxis()->SetTitle("N(Events)");
    
    ConvLength_gamma2_gamma1= new TH2D( "ConvLength_gamma2_gamma1","Leading vs Second Photon Conversion Length",bin_photonConvLength.get_nBins(), bin_photonConvLength.get_min(), bin_photonConvLength.get_max(),bin_photonConvLength.get_nBins(), bin_photonConvLength.get_min(), bin_photonConvLength.get_max() );
    ConvLength_gamma2_gamma1->GetXaxis()->SetTitle("Second Photon Distance from Vertex [cm]");
    ConvLength_gamma2_gamma1->GetYaxis()->SetTitle("Leading Photon Distance from Vertex [cm]");
     
    photonEnergy_Asymmetry = new TH1D( "photonEnergy_Asymmetry","Photon Energy Asymmetry",bin_blob_energy.get_nBins(), bin_blob_energy.get_min(), bin_blob_energy.get_max());
    photonEnergy_Asymmetry->GetXaxis()->SetTitle("Photon Energy Asymmetry - E(G2)/E(G1)");
    photonEnergy_Asymmetry->GetYaxis()->SetTitle("N(Events)");
    
    photonEnergy_Asymmetry_true = new TH1D( "photonEnergy_Asymmetry_true","True Photon Energy Asymmetry",bin_blob_energy.get_nBins(), bin_blob_energy.get_min(), bin_blob_energy.get_max());
    photonEnergy_Asymmetry_true->GetXaxis()->SetTitle("True Photon Energy Asymmetry - E(G2)/E(G1)");
    photonEnergy_Asymmetry_true->GetYaxis()->SetTitle("N(Events)");    
    
    //--------------------------------------------------------------------------
    // Invariant Mass
    //--------------------------------------------------------------------------
    invMass = new TH1D( "invMass","Reconstructed Pi0 Invariant Mass",bin_invMass.get_nBins(), bin_invMass.get_min(), bin_invMass.get_max() );
    invMass->GetXaxis()->SetTitle("Reconstructed Pi0 Invariant Mass [MeV]");
    invMass->GetYaxis()->SetTitle(Form("Candidates / %3.2f [MeV]",bin_invMass.get_width()));   

    //--------------------------------------------------------------------------
    // Momentum and Energy
    //--------------------------------------------------------------------------
    E_mc = new TH1D( "E_mc","True Pion Energy",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    E_mc->GetXaxis()->SetTitle("True Pion Energy [MeV]");
    E_mc->GetYaxis()->SetTitle(Form("Number of Pions / %3.1f [MeV]",bin_P.get_width()));
    
    E_reco = new TH1D( "E_reco","Reconstructed Pion Energy",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    E_reco->GetXaxis()->SetTitle("Reconstructed Pion Energy [MeV]");
    E_reco->GetYaxis()->SetTitle(Form("Number of Pions / %3.1f [MeV]",bin_P.get_width()));
    
    E_error = new TH1D( "E_error","Error on Pion Energy",bin_error.get_nBins(), bin_error.get_min(), bin_error.get_max() );
    E_error->GetXaxis()->SetTitle("(Reco- True) / True");
    E_error->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",bin_error.get_width()));
    
    E_reco_mc = new TH2D( "E_reco_mc","True vs Reconstructed Pion Energy",
                          bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max(),
                          bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max());
                          E_reco_mc->GetXaxis()->SetTitle("Reconstructed Pion Energy [MeV]");
                          E_reco_mc->GetYaxis()->SetTitle("True Pion Energy [MeV]");
    
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
    

}
CCProtonPi0_Pion::~CCProtonPi0_Pion()
{    
    
}

void CCProtonPi0_Pion::set_kineticEnergy(bool isMC)
{
    int type = getDataType(isMC);
    
    kineticEnergy[type] = p4[type].Energy() - restMass;
}

#endif
