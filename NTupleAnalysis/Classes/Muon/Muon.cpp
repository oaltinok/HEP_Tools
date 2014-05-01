/*
    See Muon.h header for Class Information
*/
#include "Muon.h"

using namespace std;

Muon::Muon()
{
    cout<<"Initializing Muon Particle"<<endl;
    
    // File Locations
    rootDir =   Folder_List::f_Root_Muon;
    plotDir =   Folder_List::f_Plot_Muon;
    
    cout<<"\tRoot File: "<<rootDir<<endl;
    cout<<"\tPlot Output Folder: "<<plotDir<<endl;
    
    // Create Root File 
    f = new TFile(rootDir.c_str(),"RECREATE");

    // Initialize Bins
    bin_P.setBin(100,0.0,10000.0);
    bin_KE.setBin(100,0.0,10000.0);

    cout<<"\tInitializing Histograms "<<endl;
    partScore = new TH1F( "partScore","Muon Particle Score",bin_partScore.get_nBins(), bin_partScore.get_min(), bin_partScore.get_max() );
    partScore->GetXaxis()->SetTitle("Particle Score");
    partScore->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_partScore.get_width()));
    
    P_mc = new TH1F( "P_mc","True Muon Momentum",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    P_mc->GetXaxis()->SetTitle("True Muon Momentum MeV");
    P_mc->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_P.get_width()));
    
    P_reco = new TH1F( "P_reco","Reconstructed Muon Momentum",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    P_reco->GetXaxis()->SetTitle("Reconstructed Muon Momentum MeV");
    P_reco->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_P.get_width()));
    
    P_error = new TH1F( "P_error","Error on Muon Momentum",bin_error.get_nBins(), bin_error.get_min(), bin_error.get_max() );
    P_error->GetXaxis()->SetTitle("(True - Reco) / True");
    P_error->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_error.get_width()));
    
    P_reco_mc = new TH2F( "P_reco_mc","True vs Reconstructed Muon Momentum",
                                bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max(),
                                bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max());
    P_reco_mc->GetXaxis()->SetTitle("Reconstructed Muon Momentum MeV");
    P_reco_mc->GetYaxis()->SetTitle("True Muon Momentum MeV");
    
    KE_mc = new TH1F( "KE_mc","True Muon Kinetic Energy",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    KE_mc->GetXaxis()->SetTitle("True Muon Kinetic Energy MeV");
    KE_mc->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_P.get_width()));
    
    KE_reco = new TH1F( "KE_reco","Reconstructed Muon Kinetic Energy",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    KE_reco->GetXaxis()->SetTitle("Reconstructed Muon Kinetic Energy MeV");
    KE_reco->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_P.get_width()));
    
    KE_error = new TH1F( "KE_error","Error on Muon Kinetic Energy",bin_error.get_nBins(), bin_error.get_min(), bin_error.get_max() );
    KE_error->GetXaxis()->SetTitle("(True - Reco) / True");
    KE_error->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_error.get_width()));
    
    KE_reco_mc = new TH2F( "KE_reco_mc","True vs Reconstructed Muon Kinetic Energy",
                                bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max(),
                                bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max());
    KE_reco_mc->GetXaxis()->SetTitle("Reconstructed Muon Kinetic Energy MeV");
    KE_reco_mc->GetYaxis()->SetTitle("True Muon Kinetic Energy MeV");
    
    angleMuon_mc = new TH1F( "angleMuon_mc","True Muon Angle wrt. Muon",bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max() );
    angleMuon_mc->GetXaxis()->SetTitle("True Muon Angle wrt. Muon");
    angleMuon_mc->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_angle.get_width()));
    
    angleMuon_reco = new TH1F( "angleMuon_reco","Reconstructed Muon Angle wrt. Muon",bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max() );
    angleMuon_reco->GetXaxis()->SetTitle("Reconstructed Muon Angle wrt. Muon");
    angleMuon_reco->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_angle.get_width()));
    
    angleMuon_error = new TH1F( "angleMuon_error","Error on Muon Angle wrt. Muon",bin_error.get_nBins(), bin_error.get_min(), bin_error.get_max() );
    angleMuon_error->GetXaxis()->SetTitle("(True - Reco) / True");
    angleMuon_error->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_error.get_width()));
    
    angleMuon_reco_mc = new TH2F( "angleMuon_reco_mc","True vs Reconstructed Muon Angle wrt. Muon",
                                bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max(),
                                bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max());
    angleMuon_reco_mc->GetXaxis()->SetTitle("Reconstructed Muon Angle wrt. Muon");
    angleMuon_reco_mc->GetYaxis()->SetTitle("True Muon Angle wrt. Muon");
    
    angleBeam_mc = new TH1F( "angleBeam_mc","True Muon Angle wrt. Beam",bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max() );
    angleBeam_mc->GetXaxis()->SetTitle("True Muon Angle wrt. Beam");
    angleBeam_mc->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_angle.get_width()));
    
    angleBeam_reco = new TH1F( "angleBeam_reco","Reconstructed Muon Angle wrt. Beam",bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max() );
    angleBeam_reco->GetXaxis()->SetTitle("Reconstructed Muon Angle wrt. Beam");
    angleBeam_reco->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_angle.get_width()));
    
    angleBeam_error = new TH1F( "angleBeam_error","Error on Muon Angle wrt. Beam",bin_error.get_nBins(), bin_error.get_min(), bin_error.get_max() );
    angleBeam_error->GetXaxis()->SetTitle("(True - Reco) / True");
    angleBeam_error->GetYaxis()->SetTitle(Form("Candidates / %3.1f ",bin_error.get_width()));
    
    angleBeam_reco_mc = new TH2F( "angleBeam_reco_mc","True vs Reconstructed Muon Angle wrt. Beam",
                                bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max(),
                                bin_angle.get_nBins(), bin_angle.get_min(), bin_angle.get_max());
    angleBeam_reco_mc->GetXaxis()->SetTitle("Reconstructed Muon Angle wrt. Beam");
    angleBeam_reco_mc->GetYaxis()->SetTitle("True Muon Angle wrt. Beam");
    
    cout<<"Initialization Complete! "<<endl;
}

bool Muon::get_isMinosMatched()
{
    return isMinosMatched;
}

void Muon::set_isMinosMatched(bool input)
{
    isMinosMatched = input;
}

void Muon::set_angleMuon(Particle &mu, bool isMC)
{
    // There is only 1 muon in the interaction
    // Do not need to calculate the angle wrt itself
    // Set the angle to zero
    for( int i = 0; i < N_DATA_TYPE; i++){
        angleMuon[i] = 0.0;
    }
}

void Muon::set_kineticEnergy(bool isMC)
{
    int type = getDataType(isMC);
    
    kineticEnergy[type] = p4[type].Energy() - restMass;
}




