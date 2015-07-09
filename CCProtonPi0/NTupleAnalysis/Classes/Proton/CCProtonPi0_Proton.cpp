/*
    See CCProtonPi0_Proton.h header for Class Information
*/

#ifndef CCProtonPi0_Proton_cpp
#define CCProtonPi0_Proton_cpp

#include "CCProtonPi0_Proton.h"

using namespace std;

CCProtonPi0_Proton::CCProtonPi0_Proton(int nMode) : CCProtonPi0_Particle(nMode)
{
    cout<<"Initializing CCProtonPi0_Proton"<<endl;
        
    if(nMode == 0){
        cout<<"\tNTuple Reduce Mode -- Will not create ROOT Files"<<endl;
    }else{
        // File Locations
        rootDir = Folder_List::rootOut_analyzed + branchDir + "Proton.root";
        
        cout<<"\tRoot File: "<<rootDir<<endl;
        
        // Create Root File 
        f = new TFile(rootDir.c_str(),"RECREATE");

        // Initialize Bins
        bin_E.setBin(60,0.0,3000.0);
        bin_P.setBin(40, 0.0, 2000.0);
        bin_KE.setBin(40, 0.0, 2000.0);
        bin_trackLength.setBin(250,0.0,2500.0);
        bin_trackKinked.setBin(2,0.0,2.0);
        
        initHistograms();        
    }
    cout<<"Done!"<<endl;
}


void CCProtonPi0_Proton::initHistograms()
{
    // Unique Histograms
    trackLength = new TH1D( "trackLength","Proton Track Length",bin_trackLength.get_nBins(), bin_trackLength.get_min(), bin_trackLength.get_max() );
    trackLength->GetXaxis()->SetTitle("Proton Track Length [mm]");
    trackLength->GetYaxis()->SetTitle("N(Events)");
    
    trackKinked = new TH1D( "trackKinked","Proton Track Kinked or NOT",bin_trackKinked.get_nBins(), bin_trackKinked.get_min(), bin_trackKinked.get_max() );
    trackKinked->GetXaxis()->SetTitle("Proton Track Kinked or NOT");
    trackKinked->GetYaxis()->SetTitle("N(Events)");

    // Default Histograms
    partScore = new TH1D( "partScore","Proton Particle Score",binList.particleScore_LLR.get_nBins(), binList.particleScore_LLR.get_min(), binList.particleScore_LLR.get_max() );
    partScore->GetXaxis()->SetTitle("Particle Score");
    partScore->GetYaxis()->SetTitle(Form("Number of Protons / %3.1f ",binList.particleScore_LLR.get_width()));
    
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

}
void CCProtonPi0_Proton::set_kineticEnergy(bool isMC)
{
    int type = getDataType(isMC);
    
    kineticEnergy[type] = p4[type].Energy() - restMass;
}

#endif


