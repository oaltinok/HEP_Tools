/*
    See CCProtonPi0_Pion.h header for Class Information
*/
#ifndef CCProtonPi0_Pion_cpp
#define CCProtonPi0_Pion_cpp

#include "CCProtonPi0_Pion.h"

using namespace PlotUtils;

CCProtonPi0_Pion::CCProtonPi0_Pion(bool isModeReduce, bool isMC, std::string ana_folder) : CCProtonPi0_Particle()
{
    std::cout<<"Initializing CCProtonPi0_Pion"<<std::endl;    
    
    if (isModeReduce){
        std::cout<<"\tNTuple Reduce Mode -- Will not create ROOT Files"<<std::endl;
    }else{
        // File Locations
        if (isMC) rootDir = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + ana_folder + "Pion.root";
        else rootDir = Folder_List::rootOut + Folder_List::Data + Folder_List::analyzed + ana_folder + "Pion.root";      
        
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
    MnvH1D* temp = NULL;

    for (int i = 0; i < nHistograms; i++){
        // Unique Histograms
        // Leading Photon - Energetic Photon
        temp = new MnvH1D( Form("%s_%d","gamma1_ConvLength",i),"Leading Photon Conversion Length",bin_photonConvLength.get_nBins(), bin_photonConvLength.get_min(), bin_photonConvLength.get_max() );
        temp->GetXaxis()->SetTitle("Photon Conversion Length [cm]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f [cm]",bin_photonConvLength.get_width()));
        gamma1_ConvLength.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","gamma1_P",i),"Leading Photon Momentum",bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max() );
        temp->GetXaxis()->SetTitle("P_{#gamma_{1}} [GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f [GeV]",bin_photonP.get_width()));
        gamma1_P.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","gamma1_theta",i),"Reconstructed Leading Photon Theta",binList.angle.get_nBins(), binList.angle.get_min(), binList.angle.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed #theta_{#gamma_{1}} [Degree]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.1f [Degree]",binList.angle.get_width()));
        gamma1_theta.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","gamma2_ConvLength",i),"Secondary Photon Conversion Length",bin_photonConvLength.get_nBins(), bin_photonConvLength.get_min(), bin_photonConvLength.get_max() );
        temp->GetXaxis()->SetTitle("Photon Conversion Length [cm]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f [cm]",bin_photonConvLength.get_width()));
        gamma2_ConvLength.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","gamma2_P",i),"Secondary Photon Momentum",bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max() );
        temp->GetXaxis()->SetTitle("P_{#gamma_{1}} [GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f [GeV]",bin_photonP.get_width()));
        gamma2_P.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","gamma2_theta",i),"Reconstructed Secondary Photon Theta",binList.angle.get_nBins(), binList.angle.get_min(), binList.angle.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed #theta_{#gamma_{1}} [Degree]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.1f [Degree]",binList.angle.get_width()));
        gamma2_theta.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","photonEnergy_Asymmetry",i),"Photon Energy Asymmetry",bin_photonEnergy_Asymmetry.get_nBins(), bin_photonEnergy_Asymmetry.get_min(), bin_photonEnergy_Asymmetry.get_max());
        temp->GetXaxis()->SetTitle("Photon Energy Asymmetry - E(G2)/E(G1)");
        temp->GetYaxis()->SetTitle("N(Events)");
        photonEnergy_Asymmetry.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","invMass",i),"Reconstructed Pi0 Invariant Mass",bin_invMass.get_nBins(), bin_invMass.get_min(), bin_invMass.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed m_{#gamma#gamma} [MeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f [MeV]",bin_invMass.get_width()));   
        invMass.push_back(temp);

        // Standard Histograms 
        temp = new MnvH1D( Form("%s_%d","E",i),"Reconstructed Pion Energy",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed E_{#pi^{0}} [GeV]");
        temp->GetYaxis()->SetTitle(Form("Pions / %3.1f [GeV]",bin_P.get_width()));
        E.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","P",i),"Reconstructed Pion Momentum",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed P_{#pi^{0}} [GeV]");
        temp->GetYaxis()->SetTitle(Form("Pions / %3.1f [GeV]",bin_P.get_width()));
        P.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","KE",i),"Reconstructed Pion Kinetic Energy",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed T_{#pi^{0}} [GeV]");
        temp->GetYaxis()->SetTitle(Form("Pions / %3.1f [GeV]",bin_P.get_width()));
        KE.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","theta",i),"Reconstructed Pion Theta",binList.angle.get_nBins(), binList.angle.get_min(), binList.angle.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed #theta_{#pi^{0}} [Degree]");
        temp->GetYaxis()->SetTitle(Form("Pions / %3.1f [Degree]",binList.angle.get_width()));
        theta.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","phi",i),"Reconstructed Pion Phi",binList.angle.get_nBins(), binList.angle.get_min(), binList.angle.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed #phi_{#pi^{0}} [Degree]");
        temp->GetYaxis()->SetTitle(Form("Pions / %3.1f [Degree]",binList.angle.get_width()));
        phi.push_back(temp);
    }
    // Other
    gamma1_reco_P_true_P = new TH2D( "gamma1_reco_P_true_P","Leading Photon True vs Reconstructed Momentum",bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max(),bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max());
    gamma1_reco_P_true_P->GetXaxis()->SetTitle("Reco P_{#gamma_{1}} [GeV]");
    gamma1_reco_P_true_P->GetYaxis()->SetTitle("True P_{#gamma_{1}} [GeV]");

    gamma1_P_error = new TH1D( "gamma1_P_error","Error on Leading Photon Momentum",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    gamma1_P_error->GetXaxis()->SetTitle("(P_{Reco}-P_{True})/P_{True}");
    gamma1_P_error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));

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

    mgg = new TH1D( "mgg","Corrected Pi0 Invariant Mass",bin_invMass.get_nBins(), bin_invMass.get_min(), bin_invMass.get_max() );
    mgg->GetXaxis()->SetTitle("Corrected m_{#gamma#gamma} [MeV]");
    mgg->GetYaxis()->SetTitle(Form("Events / %3.2f [MeV]",bin_invMass.get_width()));   
}

void CCProtonPi0_Pion::writeHistograms()
{
    std::cout<<">> Writing "<<rootDir<<std::endl;
    f->cd();

    for (int i = 0; i < nHistograms; i++){
        // Unique Histograms
        invMass[i]->Write();
        photonEnergy_Asymmetry[i]->Write();
        
        // Leading Photon
        gamma1_ConvLength[i]->Write();
        gamma1_P[i]->Write();
        gamma1_theta[i]->Write();
        // Secondary Photon
        gamma2_ConvLength[i]->Write();
        gamma2_P[i]->Write();
        gamma2_theta[i]->Write();


        // Standard Histograms
        E[i]->Write();
        P[i]->Write();
        KE[i]->Write();
        theta[i]->Write();
        phi[i]->Write();
    }
    
    // Photon Comparsion
    mgg->Write();
    gamma1_reco_P_true_P->Write();
    gamma1_P_error->Write();
    gamma2_reco_P_true_P->Write();
    gamma2_P_error->Write();
 
    gamma1_P_gamma2_P->Write();
    gamma1_convLength_gamma2_convLength->Write();

    f->Close();

}
#endif
