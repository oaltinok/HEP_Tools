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
        bin_photonEnergy_Asymmetry.setBin(100,0.0,1.0);
        
        initHistograms();   
    }
    
    std::cout<<"Done!"<<std::endl;
}

void CCProtonPi0_Pion::initHistograms()
{
    // Unique Histograms
    gamma1_ConvLength = new MnvH1D( "gamma1_ConvLength","Leading Photon Conversion Length",bin_photonConvLength.get_nBins(), bin_photonConvLength.get_min(), bin_photonConvLength.get_max() );
    gamma1_ConvLength->GetXaxis()->SetTitle("Photon Conversion Length [cm]");
    gamma1_ConvLength->GetYaxis()->SetTitle(Form("Events / %3.2f [cm]",bin_photonConvLength.get_width()));
    
    gamma2_ConvLength = new MnvH1D( "gamma2_ConvLength","Secondary Photon Conversion Length",bin_photonConvLength.get_nBins(), bin_photonConvLength.get_min(), bin_photonConvLength.get_max() );
    gamma2_ConvLength->GetXaxis()->SetTitle("Photon Conversion Length [cm]");
    gamma2_ConvLength->GetYaxis()->SetTitle(Form("Events / %3.2f [cm]",bin_photonConvLength.get_width()));
    
    ConvLength_gamma2_gamma1= new MnvH2D( "ConvLength_gamma2_gamma1","Leading vs Second Photon Conversion Length",bin_photonConvLength.get_nBins(), bin_photonConvLength.get_min(), bin_photonConvLength.get_max(),bin_photonConvLength.get_nBins(), bin_photonConvLength.get_min(), bin_photonConvLength.get_max() );
    ConvLength_gamma2_gamma1->GetXaxis()->SetTitle("Second Photon Distance from Vertex [cm]");
    ConvLength_gamma2_gamma1->GetYaxis()->SetTitle("Leading Photon Distance from Vertex [cm]");
     
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
    gamma1_ConvLength->Write();
    gamma2_ConvLength->Write();
    ConvLength_gamma2_gamma1->Write();
    photonEnergy_Asymmetry->Write();
    invMass->Write();
    E->Write();
    P->Write();
    KE->Write();
    theta->Write();
    phi->Write();
}
#endif
