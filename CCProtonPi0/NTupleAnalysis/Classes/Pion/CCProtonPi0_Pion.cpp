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
        bin_E.setBin(17, 0.0, 1.7);
        bin_E_Diff.setBin(100, -0.5, 0.5);
        bin_KE.setBin(30, 0.0, 3.0);
        bin_invMass.setBin(28,60,200.0);
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
        // --------------------------------------------------------------------
        // Unique Histograms
        // --------------------------------------------------------------------
        
        // Leading Photon - Energetic Photon
        temp = new MnvH1D( Form("%s_%d","gamma1_ConvLength",i),"Leading Photon Conversion Length",bin_photonConvLength.get_nBins(), bin_photonConvLength.get_min(), bin_photonConvLength.get_max() );
        temp->GetXaxis()->SetTitle("#gamma_{1} Conversion Length [cm]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f [cm]",bin_photonConvLength.get_width()));
        gamma1_ConvLength.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","gamma1_E",i),"Leading Photon Energy",bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max() );
        temp->GetXaxis()->SetTitle("E_{#gamma_{1}} [GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f [GeV]",bin_photonP.get_width()));
        gamma1_E.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","gamma1_theta",i),"Reconstructed Leading Photon Theta",binList.angle.get_nBins(), binList.angle.get_min(), binList.angle.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed #theta_{#gamma_{1}} [Degree]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.1f [Degree]",binList.angle.get_width()));
        gamma1_theta.push_back(temp);

        // Secondary Photon
        temp = new MnvH1D( Form("%s_%d","gamma2_ConvLength",i),"Secondary Photon Conversion Length",bin_photonConvLength.get_nBins(), bin_photonConvLength.get_min(), bin_photonConvLength.get_max() );
        temp->GetXaxis()->SetTitle("#gamma_{2} Conversion Length [cm]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f [cm]",bin_photonConvLength.get_width()));
        gamma2_ConvLength.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","gamma2_E",i),"Secondary Photon Energy",bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max() );
        temp->GetXaxis()->SetTitle("E_{#gamma_{2}} [GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f [GeV]",bin_photonP.get_width()));
        gamma2_E.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","gamma2_theta",i),"Reconstructed Secondary Photon Theta",binList.angle.get_nBins(), binList.angle.get_min(), binList.angle.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed #theta_{#gamma_{2}} [Degree]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.1f [Degree]",binList.angle.get_width()));
        gamma2_theta.push_back(temp);

        // Other
        temp = new MnvH1D( Form("%s_%d","photonEnergy_Asymmetry",i),"Photon Energy Asymmetry",bin_photonEnergy_Asymmetry.get_nBins(), bin_photonEnergy_Asymmetry.get_min(), bin_photonEnergy_Asymmetry.get_max());
        temp->GetXaxis()->SetTitle("Photon Energy Asymmetry - E_{#gamma_{1}}/E_{$gamma_{2}}");
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

        double binsP[17] = {0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 1.0, 1.5, 1.7};
        temp = new MnvH1D( Form("%s_%d","P",i),"Reconstructed Pion Momentum",16,binsP);
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

    double binsP[11] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0, 1.5};
    signal_P = new TH1D( "signal_P","Momentum for Signal Events",10,binsP);
    signal_P->GetXaxis()->SetTitle("Momentum [GeV]");
    signal_P->GetYaxis()->SetTitle("N(Events)");
 
    signal_theta = new TH1D( "signal_theta","Theta for Signal Events",binList.angle.get_nBins(), binList.angle.get_min(), binList.angle.get_max());
    signal_theta->GetXaxis()->SetTitle("Theta");
    signal_theta->GetYaxis()->SetTitle("N(Events)");
     

    // Truth Energy - Gamma 1
    gamma1_true_E = new TH1D( "gamma1_true_E","Leading Photon True Energy",bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max());
    gamma1_true_E->GetXaxis()->SetTitle("True E_{#gamma_{1}} [GeV]");
    gamma1_true_E->GetYaxis()->SetTitle("N(Events)");
 
    gamma1_reco_error_E = new TH1D( "gamma1_reco_error_E","Error on Reconstructed Energy",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    gamma1_reco_error_E->GetXaxis()->SetTitle("(E_{Reco}-E_{True})/E_{True}");
    gamma1_reco_error_E->GetYaxis()->SetTitle("N(Events)");
 
    gamma1_reco_E_true_E = new TH2D( "gamma1_reco_E_true_E","Leading Photon True vs Reconstructed Energy",bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max(),bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max());
    gamma1_reco_E_true_E->GetXaxis()->SetTitle("Reco E_{#gamma_{1}} [GeV]");
    gamma1_reco_E_true_E->GetYaxis()->SetTitle("True E_{#gamma_{1}} [GeV]");

    gamma1_true_E_reco_E_error = new TH2D( "gamma1_true_E_reco_E_error","(E_{Reco}-E_{True})/E_{True} vs E_{True}",bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max(),binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    gamma1_true_E_reco_E_error->GetXaxis()->SetTitle("E_{True}");
    gamma1_true_E_reco_E_error->GetYaxis()->SetTitle("(E_{Reco}-E_{True})/E_{True}");

    // Truth Energy - Gamma 2
    gamma2_true_E = new TH1D( "gamma2_true_E","Secondary Photon True Energy",bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max());
    gamma2_true_E->GetXaxis()->SetTitle("True E_{#gamma_{1}} [GeV]");
    gamma2_true_E->GetYaxis()->SetTitle("N(Events)");
 
    gamma2_reco_error_E = new TH1D( "gamma2_reco_error_E","Error on Reconstructed Energy",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    gamma2_reco_error_E->GetXaxis()->SetTitle("(E_{Reco}-E_{True})/E_{True}");
    gamma2_reco_error_E->GetYaxis()->SetTitle("N(Events)");
 
    gamma2_reco_E_true_E = new TH2D( "gamma2_reco_E_true_E","Secondary Photon True vs Reconstructed Energy",bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max(),bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max());
    gamma2_reco_E_true_E->GetXaxis()->SetTitle("Reco E_{#gamma_{1}} [GeV]");
    gamma2_reco_E_true_E->GetYaxis()->SetTitle("True E_{#gamma_{1}} [GeV]");

    gamma2_true_E_reco_E_error = new TH2D( "gamma2_true_E_reco_E_error","(E_{Reco}-E_{True})/E_{True} vs E_{True}",bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max(),binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    gamma2_true_E_reco_E_error->GetXaxis()->SetTitle("E_{True}");
    gamma2_true_E_reco_E_error->GetYaxis()->SetTitle("(E_{Reco}-E_{True})/E_{True}");

    // Other
    gamma1_convLength_gamma2_convLength= new TH2D( "gamma1_convLength_gamma2_convLength","Leading vs Second Photon Conversion Length",bin_photonConvLength.get_nBins(), bin_photonConvLength.get_min(), bin_photonConvLength.get_max(),bin_photonConvLength.get_nBins(), bin_photonConvLength.get_min(), bin_photonConvLength.get_max() );
    gamma1_convLength_gamma2_convLength->GetXaxis()->SetTitle("Leading Photon Distance from Vertex [cm]");
    gamma1_convLength_gamma2_convLength->GetYaxis()->SetTitle("Second Photon Distance from Vertex [cm]");
     
    signal_gamma1_E_gamma2_E = new TH2D( "signal_gamma1_E_gamma2_E","Signal: Gamma1 Energy vs Gamma2 Energy",20,0.0,200.0,20,0.0,200.0);
    signal_gamma1_E_gamma2_E->GetXaxis()->SetTitle("Reconstructed E_{#gamma_{1}} [MeV]");
    signal_gamma1_E_gamma2_E->GetYaxis()->SetTitle("Reconstructed E_{#gamma_{2}} [MeV]");
 
    bckg_gamma1_E_gamma2_E = new TH2D( "bckg_gamma1_E_gamma2_E","Background: Gamma1 Energy vs Gamma2 Energy",20,0.0,200.0,20,0.0,200.0);
    bckg_gamma1_E_gamma2_E->GetXaxis()->SetTitle("Reconstructed E_{#gamma_{1}} [MeV]");
    bckg_gamma1_E_gamma2_E->GetYaxis()->SetTitle("Reconstructed E_{#gamma_{2}} [MeV]");

    bckg_signal_diff = new TH2D( "bckg_signal_diff","Background - Signal: Gamma1 Energy vs Gamma2 Energy",20,0.0,200.0,20,0.0,200.0);
    bckg_signal_diff->GetXaxis()->SetTitle("Reconstructed E_{#gamma_{1}} [MeV]");
    bckg_signal_diff->GetYaxis()->SetTitle("Reconstructed E_{#gamma_{2}} [MeV]");

    reco_P_true_P = new TH2D( "reco_P_true_P","True vs Reconstructed #pi^0 Momentum",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max(), bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max());
    reco_P_true_P->GetXaxis()->SetTitle("Reconstructed P_{#pi^0} [GeV]");
    reco_P_true_P->GetYaxis()->SetTitle("True P_{#pi^0} [GeV]");

    P_error = new TH1D( "P_error","Error on #pi^0 Momentum",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    P_error->GetXaxis()->SetTitle("(P_{Reco}-P_{True})/P_{True}");
    P_error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));

    reco_E_true_E = new TH2D( "reco_E_true_E","True vs Reconstructed #pi^0 Energy",bin_E.get_nBins(), bin_E.get_min(), bin_E.get_max(), bin_E.get_nBins(), bin_E.get_min(), bin_E.get_max());
    reco_E_true_E->GetXaxis()->SetTitle("Reconstructed E_{#pi^0} [GeV]");
    reco_E_true_E->GetYaxis()->SetTitle("True E_{#pi^0} [GeV]");

    E_error = new TH1D( "E_error","Error on #pi^0 Energy",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    E_error->GetXaxis()->SetTitle("(E_{Reco}-E_{True})/E_{True}");
    E_error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));

    E_Diff = new TH1D( "E_Diff","Error on #pi^0 Energy",bin_E_Diff.get_nBins(), bin_E_Diff.get_min(), bin_E_Diff.get_max() );
    E_Diff->GetXaxis()->SetTitle("(E_{Reco}-E_{True})/E_{True}");
    E_Diff->GetYaxis()->SetTitle(Form("Events / %3.2f ",bin_E_Diff.get_width()));
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
        gamma1_E[i]->Write();
        gamma1_theta[i]->Write();
       
        // Secondary Photon
        gamma2_ConvLength[i]->Write();
        gamma2_E[i]->Write();
        gamma2_theta[i]->Write();
       
        // Standard Histograms
        E[i]->Write();
        P[i]->Write();
        KE[i]->Write();
        theta[i]->Write();
        phi[i]->Write();
    }
    signal_P->Write();
    signal_theta->Write();

    // Photon Comparison
    bckg_signal_diff->Add(signal_gamma1_E_gamma2_E, -1);
    bckg_signal_diff->Add(bckg_gamma1_E_gamma2_E, +1);
    signal_gamma1_E_gamma2_E->Write();
    bckg_gamma1_E_gamma2_E->Write();
    bckg_signal_diff->Write();

    // Gamma1 True Energy
    gamma1_true_E->Write();
    gamma1_reco_error_E->Write();
    gamma1_reco_E_true_E->Write();
    gamma1_true_E_reco_E_error->Write();
  
    // Gamma2 True Energy
    gamma2_true_E->Write();
    gamma2_reco_error_E->Write();
    gamma2_reco_E_true_E->Write();
    gamma2_true_E_reco_E_error->Write();

    reco_P_true_P->Write();
    P_error->Write();

    reco_E_true_E->Write();
    E_error->Write();
    
    E_Diff->Write();

    f->Close();
}
#endif
