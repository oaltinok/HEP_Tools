/*
    See CCProtonPi0_Muon.h header for Class Information
*/
#ifndef CCProtonPi0_Muon_cpp
#define CCProtonPi0_Muon_cpp

#include "CCProtonPi0_Muon.h"

using namespace std;

CCProtonPi0_Muon::CCProtonPi0_Muon(int nMode) : CCProtonPi0_Particle(nMode)
{
    cout<<"Initializing CCProtonPi0_Muon"<<endl;
      
    if(nMode == 0){
        cout<<"\tNTuple Reduce Mode -- Will not create ROOT Files"<<endl;
    }else{
        // File Locations
        rootDir = Folder_List::rootOut_analyzed + branchDir + "Muon.root";

        cout<<"\tRoot File: "<<rootDir<<endl;
        
        // Create Root File 
        f = new TFile(rootDir.c_str(),"RECREATE");

        // Initialize Bins
        bin_P.setBin(10,0.0,10.0);
        bin_KE.setBin(100,0.0,10.0);
        bin_AngleBeam.setBin(90,0.0,90.0);
        
        initHistograms();
    }
    
    cout<<"Done!"<<endl;
}


void CCProtonPi0_Muon::initHistograms()
{
    partScore = new TH1D( "partScore","Muon Particle Score",bin_partScore.get_nBins(), bin_partScore.get_min(), bin_partScore.get_max() );
    partScore->GetXaxis()->SetTitle("Particle Score");
    partScore->GetYaxis()->SetTitle(Form("Number of Muons / %3.1f ",bin_partScore.get_width()));
    
    E_mc = new TH1D( "E_mc","True Muon Energy",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    E_mc->GetXaxis()->SetTitle("True Muon Energy [GeV]");
    E_mc->GetYaxis()->SetTitle(Form("Number of Muons / %3.1f [GeV] ",bin_P.get_width()));
    
    E_reco = new TH1D( "E_reco","Reconstructed Muon Energy",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    E_reco->GetXaxis()->SetTitle("Reconstructed Muon Energy [GeV]");
    E_reco->GetYaxis()->SetTitle(Form("Number of Muons / %3.1f [GeV] ",bin_P.get_width()));
    
    E_error = new TH1D( "E_error","Error on Muon Energy",bin_error.get_nBins(), bin_error.get_min(), bin_error.get_max() );
    E_error->GetXaxis()->SetTitle("(Reco- True) / True");
    E_error->GetYaxis()->SetTitle(Form("Candidates  / %3.2f",bin_error.get_width()));
    
    E_reco_mc = new TH2D( "E_reco_mc","True vs Reconstructed Muon Energy",
                          bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max(),
                          bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max());
                          E_reco_mc->GetXaxis()->SetTitle("Reconstructed Muon Energy [GeV]");
                          E_reco_mc->GetYaxis()->SetTitle("True Muon Energy [GeV]");
    
    P_mc = new TH1D( "P_mc","True Muon Momentum",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    P_mc->GetXaxis()->SetTitle("True Muon Momentum [GeV]");
    P_mc->GetYaxis()->SetTitle(Form("Number of Muons / %3.1f [GeV] ",bin_P.get_width()));
    
    P_reco = new TH1D( "P_reco","Reconstructed Muon Momentum",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    P_reco->GetXaxis()->SetTitle("Reconstructed Muon Momentum [GeV]");
    P_reco->GetYaxis()->SetTitle(Form("Number of Muons / %3.1f [GeV] ",bin_P.get_width()));
    
    P_error = new TH1D( "P_error","Error on Muon Momentum",bin_error.get_nBins(), bin_error.get_min(), bin_error.get_max() );
    P_error->GetXaxis()->SetTitle("(Reco- True) / True");
    P_error->GetYaxis()->SetTitle(Form("Candidates  / %3.2f",bin_error.get_width()));
    
    P_reco_mc = new TH2D( "P_reco_mc","True vs Reconstructed Muon Momentum",
                                bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max(),
                                bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max());
    P_reco_mc->GetXaxis()->SetTitle("Reconstructed Muon Momentum [GeV]");
    P_reco_mc->GetYaxis()->SetTitle("True Muon Momentum [GeV]");
    
    KE_mc = new TH1D( "KE_mc","True Muon Kinetic Energy",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    KE_mc->GetXaxis()->SetTitle("True Muon Kinetic Energy [GeV]");
    KE_mc->GetYaxis()->SetTitle(Form("Number of Muons / %3.1f [GeV]",bin_P.get_width()));
    
    KE_reco = new TH1D( "KE_reco","Reconstructed Muon Kinetic Energy",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    KE_reco->GetXaxis()->SetTitle("Reconstructed Muon Kinetic Energy [GeV]");
    KE_reco->GetYaxis()->SetTitle(Form("Number of Muons / %3.1f [GeV]",bin_P.get_width()));
    
    KE_error = new TH1D( "KE_error","Error on Muon Kinetic Energy",bin_error.get_nBins(), bin_error.get_min(), bin_error.get_max() );
    KE_error->GetXaxis()->SetTitle("(Reco- True) / True");
    KE_error->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",bin_error.get_width()));
    
    KE_reco_mc = new TH2D( "KE_reco_mc","True vs Reconstructed Muon Kinetic Energy",
                                bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max(),
                                bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max());
    KE_reco_mc->GetXaxis()->SetTitle("Reconstructed Muon Kinetic Energy [GeV]");
    KE_reco_mc->GetYaxis()->SetTitle("True Muon Kinetic Energy [GeV]");
    
    angleMuon_mc = new TH1D( "angleMuon_mc","True Muon Angle wrt. Muon",bin_AngleBeam.get_nBins(), bin_AngleBeam.get_min(), bin_AngleBeam.get_max() );
    angleMuon_mc->GetXaxis()->SetTitle("True Muon Angle wrt. Muon [Degree]");
    angleMuon_mc->GetYaxis()->SetTitle(Form("Number of Muons / %3.1f [Degree] ",bin_AngleBeam.get_width()));
    
    angleMuon_reco = new TH1D( "angleMuon_reco","Reconstructed Muon Angle wrt. Muon",bin_AngleBeam.get_nBins(), bin_AngleBeam.get_min(), bin_AngleBeam.get_max() );
    angleMuon_reco->GetXaxis()->SetTitle("Reconstructed Muon Angle wrt. Muon [Degree]");
    angleMuon_reco->GetYaxis()->SetTitle(Form("Number of Muons / %3.1f [Degree] ",bin_AngleBeam.get_width()));
    
    angleMuon_error = new TH1D( "angleMuon_error","Error on Muon Angle wrt. Muon",bin_error.get_nBins(), bin_error.get_min(), bin_error.get_max() );
    angleMuon_error->GetXaxis()->SetTitle("(Reco- True) / True");
    angleMuon_error->GetYaxis()->SetTitle(Form("Candidates  / %3.2f ",bin_error.get_width()));
    
    angleMuon_reco_mc = new TH2D( "angleMuon_reco_mc","True vs Reconstructed Muon Angle wrt. Muon",
                                bin_AngleBeam.get_nBins(), bin_AngleBeam.get_min(), bin_AngleBeam.get_max(),
                                bin_AngleBeam.get_nBins(), bin_AngleBeam.get_min(), bin_AngleBeam.get_max());
    angleMuon_reco_mc->GetXaxis()->SetTitle("Reconstructed Muon Angle wrt. Muon [Degree]");
    angleMuon_reco_mc->GetYaxis()->SetTitle("True Muon Angle wrt. Muon [Degree]");
    
    angleBeam_mc = new TH1D( "angleBeam_mc","True Muon Angle wrt. Beam",bin_AngleBeam.get_nBins(), bin_AngleBeam.get_min(), bin_AngleBeam.get_max() );
    angleBeam_mc->GetXaxis()->SetTitle("True Muon Angle wrt. Beam [Degree]" );
    angleBeam_mc->GetYaxis()->SetTitle(Form("Number of Muons / %3.1f [Degree]",bin_AngleBeam.get_width()));
    
    angleBeam_reco = new TH1D( "angleBeam_reco","Reconstructed Muon Angle wrt. Beam",bin_AngleBeam.get_nBins(), bin_AngleBeam.get_min(), bin_AngleBeam.get_max() );
    angleBeam_reco->GetXaxis()->SetTitle("Reconstructed Muon Angle wrt. Beam [Degree]");
    angleBeam_reco->GetYaxis()->SetTitle(Form("Number of Muons / %3.1f [Degree]",bin_AngleBeam.get_width()));
    
    angleBeam_error = new TH1D( "angleBeam_error","Error on Muon Angle wrt. Beam",bin_error.get_nBins(), bin_error.get_min(), bin_error.get_max() );
    angleBeam_error->GetXaxis()->SetTitle("(Reco- True) / True");
    angleBeam_error->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",bin_AngleBeam.get_width()));
    
    angleBeam_reco_mc = new TH2D( "angleBeam_reco_mc","True vs Reconstructed Muon Angle wrt. Beam",
                                bin_AngleBeam.get_nBins(), bin_AngleBeam.get_min(), bin_AngleBeam.get_max(),
                                bin_AngleBeam.get_nBins(), bin_AngleBeam.get_min(), bin_AngleBeam.get_max());
    angleBeam_reco_mc->GetXaxis()->SetTitle("Reconstructed Muon Angle wrt. Beam [Degree]");
    angleBeam_reco_mc->GetYaxis()->SetTitle("True Muon Angle wrt. Beam [Degree]");

}

bool CCProtonPi0_Muon::get_isMinosMatched()
{
    return isMinosMatched;
}

void CCProtonPi0_Muon::set_isMinosMatched(bool input)
{
    isMinosMatched = input;
}

void CCProtonPi0_Muon::set_angleMuon(CCProtonPi0_Particle &mu, bool isMC)
{
    // There is only 1 CCProtonPi0_Muon in the interaction
    // Do not need to calculate the angle wrt itself
    // Set the angle to zero
    for( int i = 0; i < N_DATA_TYPE; i++){
        angleMuon[i] = 0.0;
    }
    
    // Silence Compiler Warning
    (void)mu;
    (void)isMC;   
}

void CCProtonPi0_Muon::set_kineticEnergy(bool isMC)
{
    int type = getDataType(isMC);
    
    kineticEnergy[type] = p4[type].Energy() - restMass;
}

#endif



