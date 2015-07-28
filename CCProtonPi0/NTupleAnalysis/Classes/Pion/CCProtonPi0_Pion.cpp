/*
    See CCProtonPi0_Pion.h header for Class Information
*/
#ifndef CCProtonPi0_Pion_cpp
#define CCProtonPi0_Pion_cpp

#include "CCProtonPi0_Pion.h"

using namespace PlotUtils;

CCProtonPi0_Pion::CCProtonPi0_Pion(int nMode, bool isMC) : CCProtonPi0_Particle(nMode)
{
    std::cout<<"Initializing CCProtonPi0_Pion"<<std::endl;    
    
    if(nMode == 0){
        std::cout<<"\tNTuple Reduce Mode -- Will not create ROOT Files"<<std::endl;
    }else{
        // File Locations
        if (isMC) rootDir = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + branchDir + "Pion.root";
        else rootDir = Folder_List::rootOut + Folder_List::Data + Folder_List::analyzed + branchDir + "Pion.root";      
        
        std::cout<<"\tRoot File: "<<rootDir<<std::endl;
        
        // Create Root File 
        f = new TFile(rootDir.c_str(),"RECREATE");

        // Initialize Bins
        bin_P.setBin(17, 0.0, 1.7);
        bin_KE.setBin(30, 0.0, 3.0);
        bin_invMass.setBin(60,0.0,600.0);
        bin_photonConvLength.setBin(50,0.0,100.0);
        bin_photonP.setBin(20,0.0,1.0);
        bin_photonEnergy_Asymmetry.setBin(100,0.0,1.0);
        
        initHistograms();   
    }
    
    std::cout<<"Done!"<<std::endl;
}

void CCProtonPi0_Pion::initHistograms()
{
    // Unique Histograms
    // Leading Photon - Energetic Photon
    gamma1_ConvLength = new MnvH1D( "gamma1_ConvLength","Leading Photon Conversion Length",bin_photonConvLength.get_nBins(), bin_photonConvLength.get_min(), bin_photonConvLength.get_max() );
    gamma1_ConvLength->GetXaxis()->SetTitle("Photon Conversion Length [cm]");
    gamma1_ConvLength->GetYaxis()->SetTitle(Form("Events / %3.2f [cm]",bin_photonConvLength.get_width()));
    
    gamma1_P = new MnvH1D( "gamma1_P","Leading Photon Momentum",bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max() );
    gamma1_P->GetXaxis()->SetTitle("P_{#gamma_{1}} [GeV]");
    gamma1_P->GetYaxis()->SetTitle(Form("Events / %3.2f [GeV]",bin_photonP.get_width()));

    gamma1_theta = new MnvH1D( "gamma1_theta","Reconstructed Leading Photon Theta",binList.angle.get_nBins(), binList.angle.get_min(), binList.angle.get_max() );
    gamma1_theta->GetXaxis()->SetTitle("Reconstructed #theta_{#gamma_{1}} [Degree]");
    gamma1_theta->GetYaxis()->SetTitle(Form("Events / %3.1f [Degree]",binList.angle.get_width()));

    gamma1_reco_P_true_P = new TH2D( "gamma1_reco_P_true_P","Leading Photon True vs Reconstructed Momentum",bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max(),bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max());
    gamma1_reco_P_true_P->GetXaxis()->SetTitle("Reco P_{#gamma_{1}} [GeV]");
    gamma1_reco_P_true_P->GetYaxis()->SetTitle("True P_{#gamma_{1}} [GeV]");

    gamma1_P_error = new TH1D( "gamma1_P_error","Error on Leading Photon Momentum",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    gamma1_P_error->GetXaxis()->SetTitle("(P_{Reco}-P_{True})/P_{True}");
    gamma1_P_error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));

    gamma2_ConvLength = new MnvH1D( "gamma2_ConvLength","Secondary Photon Conversion Length",bin_photonConvLength.get_nBins(), bin_photonConvLength.get_min(), bin_photonConvLength.get_max() );
    gamma2_ConvLength->GetXaxis()->SetTitle("Photon Conversion Length [cm]");
    gamma2_ConvLength->GetYaxis()->SetTitle(Form("Events / %3.2f [cm]",bin_photonConvLength.get_width()));
    
    gamma2_P = new MnvH1D( "gamma2_P","Secondary Photon Momentum",bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max() );
    gamma2_P->GetXaxis()->SetTitle("P_{#gamma_{1}} [GeV]");
    gamma2_P->GetYaxis()->SetTitle(Form("Events / %3.2f [GeV]",bin_photonP.get_width()));

    gamma2_theta = new MnvH1D( "gamma2_theta","Reconstructed Secondary Photon Theta",binList.angle.get_nBins(), binList.angle.get_min(), binList.angle.get_max() );
    gamma2_theta->GetXaxis()->SetTitle("Reconstructed #theta_{#gamma_{1}} [Degree]");
    gamma2_theta->GetYaxis()->SetTitle(Form("Events / %3.1f [Degree]",binList.angle.get_width()));

    gamma2_reco_P_true_P = new TH2D( "gamma2_reco_P_true_P","Secondary Photon True vs Reconstructed Momentum",bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max(),bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max() );
    gamma2_reco_P_true_P->GetXaxis()->SetTitle("Reco P_{#gamma_{1}} [GeV]");
    gamma2_reco_P_true_P->GetYaxis()->SetTitle("True P_{#gamma_{1}} [GeV]");
 
    gamma2_P_error = new TH1D( "gamma2_P_error","Error on Secondary Photon Momentum",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    gamma2_P_error->GetXaxis()->SetTitle("(P_{Reco}-P_{True})/P_{True}");
    gamma2_P_error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));
   
    gamma1_convLength_gamma2_convLength= new TH2D( "gamma1_convLength_gamma2_convLength","Leading vs Second Photon Conversion Length",bin_photonConvLength.get_nBins(), bin_photonConvLength.get_min(), bin_photonConvLength.get_max(),bin_photonConvLength.get_nBins(), bin_photonConvLength.get_min(), bin_photonConvLength.get_max() );
    gamma1_convLength_gamma2_convLength->GetXaxis()->SetTitle("Leading Photon Distance from Vertex [cm]");
    gamma1_convLength_gamma2_convLength->GetYaxis()->SetTitle("Second Photon Distance from Vertex [cm]");
     
    gamma1_P_gamma2_P = new TH2D( "gamma1_P_gamma2_P","Leading Photon vs Secondary Photon Momentum",bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max(),bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max() );
    gamma1_P_gamma2_P->GetXaxis()->SetTitle("Reconstructed P_{#gamma_{1}} [GeV]");
    gamma1_P_gamma2_P->GetYaxis()->SetTitle("Reconstructed P_{#gamma_{2}} [GeV]");
    
    photonEnergy_Asymmetry = new MnvH1D( "photonEnergy_Asymmetry","Photon Energy Asymmetry",bin_photonEnergy_Asymmetry.get_nBins(), bin_photonEnergy_Asymmetry.get_min(), bin_photonEnergy_Asymmetry.get_max());
    photonEnergy_Asymmetry->GetXaxis()->SetTitle("Photon Energy Asymmetry - E(G2)/E(G1)");
    photonEnergy_Asymmetry->GetYaxis()->SetTitle("N(Events)");
   
    invMass = new MnvH1D( "invMass","Reconstructed Pi0 Invariant Mass",bin_invMass.get_nBins(), bin_invMass.get_min(), bin_invMass.get_max() );
    invMass->GetXaxis()->SetTitle("Reconstructed m_{#gamma#gamma} [MeV]");
    invMass->GetYaxis()->SetTitle(Form("Events / %3.2f [MeV]",bin_invMass.get_width()));   

    // Standard Histograms 
    E = new MnvH1D( "E","Reconstructed Pion Energy",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    E->GetXaxis()->SetTitle("Reconstructed E_{#pi^{0}} [GeV]");
    E->GetYaxis()->SetTitle(Form("Pions / %3.1f [GeV]",bin_P.get_width()));
      
    P = new MnvH1D( "P","Reconstructed Pion Momentum",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    P->GetXaxis()->SetTitle("Reconstructed P_{#pi^{0}} [GeV]");
    P->GetYaxis()->SetTitle(Form("Pions / %3.1f [GeV]",bin_P.get_width()));
    
    KE = new MnvH1D( "KE","Reconstructed Pion Kinetic Energy",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    KE->GetXaxis()->SetTitle("Reconstructed T_{#pi^{0}} [GeV]");
    KE->GetYaxis()->SetTitle(Form("Pions / %3.1f [GeV]",bin_P.get_width()));
    
    theta = new MnvH1D( "theta","Reconstructed Pion Theta",binList.angle.get_nBins(), binList.angle.get_min(), binList.angle.get_max() );
    theta->GetXaxis()->SetTitle("Reconstructed #theta_{#pi^{0}} [Degree]");
    theta->GetYaxis()->SetTitle(Form("Pions / %3.1f [Degree]",binList.angle.get_width()));
    
    phi = new MnvH1D( "phi","Reconstructed Pion Phi",binList.angle.get_nBins(), binList.angle.get_min(), binList.angle.get_max() );
    phi->GetXaxis()->SetTitle("Reconstructed #phi_{#pi^{0}} [Degree]");
    phi->GetYaxis()->SetTitle(Form("Pions / %3.1f [Degree]",binList.angle.get_width()));

}

void CCProtonPi0_Pion::writeHistograms()
{
    std::cout<<">> Writing "<<rootDir<<std::endl;
    f->cd();

    // Unique Histograms
    invMass->Write();
    
    // Leading Photon
    gamma1_ConvLength->Write();
    gamma1_P->Write();
    gamma1_theta->Write();
    gamma1_reco_P_true_P->Write();
    gamma1_P_error->Write();

    // Secondary Photon
    gamma2_ConvLength->Write();
    gamma2_P->Write();
    gamma2_theta->Write();
    gamma2_reco_P_true_P->Write();
    gamma2_P_error->Write();
 
    // Photon Comparsion
    gamma1_P_gamma2_P->Write();
    gamma1_convLength_gamma2_convLength->Write();
    photonEnergy_Asymmetry->Write();
    
    // Standard Histograms
    E->Write();
    P->Write();
    KE->Write();
    theta->Write();
    phi->Write();
}
#endif
