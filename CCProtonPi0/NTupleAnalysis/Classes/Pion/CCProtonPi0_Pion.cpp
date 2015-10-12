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
        bin_photonP.setBin(80,0.0,1.0);
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

        // Truth Match
        temp = new MnvH1D( Form("%s_%d","evis_frac_true_pi0_reco_all",i),"Visible Energy #pi^0 Fraction",binList.fraction.get_nBins(), binList.fraction.get_min(), binList.fraction.get_max() );
        temp->GetXaxis()->SetTitle("True E_{vis}^{#pi^{0}} / Reco E_{vis}^{Total}");
        temp->GetYaxis()->SetTitle("N(Events)");
        evis_frac_true_pi0_reco_all.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","evis_frac_reco_pi0_true_pi0",i),"Visible Energy #pi^0 Fraction",binList.fraction.get_nBins(), binList.fraction.get_min(), binList.fraction.get_max() );
        temp->GetXaxis()->SetTitle("Reco E_{vis}^{#pi^{0}} / True E_{vis}^{#pi^{0}}");
        temp->GetYaxis()->SetTitle("N(Events)");
        evis_frac_reco_pi0_true_pi0.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","evis_frac_reco_pi0_reco_all",i),"Visible Energy #pi^0 Fraction",binList.fraction.get_nBins(), binList.fraction.get_min(), binList.fraction.get_max() );
        temp->GetXaxis()->SetTitle("Reco E_{vis}^{#pi^{0}} / Reco E_{vis}^{Total}");
        temp->GetYaxis()->SetTitle("N(Events)");
        evis_frac_reco_pi0_reco_all.push_back(temp);
 
        temp = new MnvH1D( Form("%s_%d","evis_frac_reco_nonpi0_reco_all",i),"Visible Energy Non-#pi^0 Fraction",binList.fraction.get_nBins(), binList.fraction.get_min(), binList.fraction.get_max() );
        temp->GetXaxis()->SetTitle("Reco E_{vis}^{Non-#pi^{0}} / Reco E_{vis}^{Total}");
        temp->GetYaxis()->SetTitle("N(Events)");
        evis_frac_reco_nonpi0_reco_all.push_back(temp);

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

        temp = new MnvH1D( Form("%s_%d","gamma2_ConvLength",i),"Secondary Photon Conversion Length",bin_photonConvLength.get_nBins(), bin_photonConvLength.get_min(), bin_photonConvLength.get_max() );
        temp->GetXaxis()->SetTitle("#gamma_{2} Conversion Length [cm]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f [cm]",bin_photonConvLength.get_width()));
        gamma2_ConvLength.push_back(temp);

        // Secondary Photon
        temp = new MnvH1D( Form("%s_%d","gamma2_E",i),"Secondary Photon Energy",bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max() );
        temp->GetXaxis()->SetTitle("E_{#gamma_{2}} [GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f [GeV]",bin_photonP.get_width()));
        gamma2_E.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","gamma2_theta",i),"Reconstructed Secondary Photon Theta",binList.angle.get_nBins(), binList.angle.get_min(), binList.angle.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed #theta_{#gamma_{2}} [Degree]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.1f [Degree]",binList.angle.get_width()));
        gamma2_theta.push_back(temp);

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

    gamma1_true_E_reco_E_ratio = new TH2D( "gamma1_true_E_reco_E_ratio","E_{Reco}/E_{True} vs E_{True}",bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max(),binList.ratio.get_nBins(), binList.ratio.get_min(), binList.ratio.get_max() );
    gamma1_true_E_reco_E_ratio->GetXaxis()->SetTitle("E_{True}");
    gamma1_true_E_reco_E_ratio->GetYaxis()->SetTitle("E_{Reco}/E_{True}");

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

   gamma2_true_E_reco_E_ratio = new TH2D( "gamma2_true_E_reco_E_ratio","E_{Reco}/E_{True} vs E_{True}",bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max(),binList.ratio.get_nBins(), binList.ratio.get_min(), binList.ratio.get_max() );
    gamma2_true_E_reco_E_ratio->GetXaxis()->SetTitle("E_{True}");
    gamma2_true_E_reco_E_ratio->GetYaxis()->SetTitle("E_{Reco}/E_{True}");

    // Other
    gamma1_convLength_gamma2_convLength= new TH2D( "gamma1_convLength_gamma2_convLength","Leading vs Second Photon Conversion Length",bin_photonConvLength.get_nBins(), bin_photonConvLength.get_min(), bin_photonConvLength.get_max(),bin_photonConvLength.get_nBins(), bin_photonConvLength.get_min(), bin_photonConvLength.get_max() );
    gamma1_convLength_gamma2_convLength->GetXaxis()->SetTitle("Leading Photon Distance from Vertex [cm]");
    gamma1_convLength_gamma2_convLength->GetYaxis()->SetTitle("Second Photon Distance from Vertex [cm]");
     
    gamma1_E_gamma2_E = new TH2D( "gamma1_E_gamma2_E","Leading Photon vs Secondary Photon Energy",bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max(),bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max() );
    gamma1_E_gamma2_E->GetXaxis()->SetTitle("Reconstructed E_{#gamma_{1}} [GeV]");
    gamma1_E_gamma2_E->GetYaxis()->SetTitle("Reconstructed E_{#gamma_{2}} [GeV]");

    // ------------------------------------------------------------------------
    // Truth Match
    // ------------------------------------------------------------------------
    g1_evis_proton = new TH1D("g1_evis_proton","Visible Energy",binList.gamma_evis_pdg.get_nBins(), binList.gamma_evis_pdg.get_min(), binList.gamma_evis_pdg.get_max() );
    g1_evis_proton->GetXaxis()->SetTitle("Evis");
    g1_evis_proton->GetYaxis()->SetTitle("N(Events)");

    g1_evis_neutron = new TH1D("g1_evis_neutron","Visible Energy",binList.gamma_evis_pdg.get_nBins(), binList.gamma_evis_pdg.get_min(), binList.gamma_evis_pdg.get_max() );
    g1_evis_neutron->GetXaxis()->SetTitle("Evis");
    g1_evis_neutron->GetYaxis()->SetTitle("N(Events)");

    g1_evis_pi = new TH1D("g1_evis_pi","Visible Energy",binList.gamma_evis_pdg.get_nBins(), binList.gamma_evis_pdg.get_min(), binList.gamma_evis_pdg.get_max() );
    g1_evis_pi->GetXaxis()->SetTitle("Evis");
    g1_evis_pi->GetYaxis()->SetTitle("N(Events)");

    g1_evis_pi0 = new TH1D("g1_evis_pi0","Visible Energy",binList.gamma_evis_pdg.get_nBins(), binList.gamma_evis_pdg.get_min(), binList.gamma_evis_pdg.get_max() );
    g1_evis_pi0->GetXaxis()->SetTitle("Evis");
    g1_evis_pi0->GetYaxis()->SetTitle("N(Events)");

    g1_evis_muon = new TH1D("g1_evis_muon","Visible Energy",binList.gamma_evis_pdg.get_nBins(), binList.gamma_evis_pdg.get_min(), binList.gamma_evis_pdg.get_max() );
    g1_evis_muon->GetXaxis()->SetTitle("Evis");
    g1_evis_muon->GetYaxis()->SetTitle("N(Events)");

    g2_evis_proton = new TH1D("g2_evis_proton","Visible Energy",binList.gamma_evis_pdg.get_nBins(), binList.gamma_evis_pdg.get_min(), binList.gamma_evis_pdg.get_max() );
    g2_evis_proton->GetXaxis()->SetTitle("Evis");
    g2_evis_proton->GetYaxis()->SetTitle("N(Events)");

    g2_evis_neutron = new TH1D("g2_evis_neutron","Visible Energy",binList.gamma_evis_pdg.get_nBins(), binList.gamma_evis_pdg.get_min(), binList.gamma_evis_pdg.get_max() );
    g2_evis_neutron->GetXaxis()->SetTitle("Evis");
    g2_evis_neutron->GetYaxis()->SetTitle("N(Events)");

    g2_evis_pi = new TH1D("g2_evis_pi","Visible Energy",binList.gamma_evis_pdg.get_nBins(), binList.gamma_evis_pdg.get_min(), binList.gamma_evis_pdg.get_max() );
    g2_evis_pi->GetXaxis()->SetTitle("Evis");
    g2_evis_pi->GetYaxis()->SetTitle("N(Events)");

    g2_evis_pi0 = new TH1D("g2_evis_pi0","Visible Energy",binList.gamma_evis_pdg.get_nBins(), binList.gamma_evis_pdg.get_min(), binList.gamma_evis_pdg.get_max() );
    g2_evis_pi0->GetXaxis()->SetTitle("Evis");
    g2_evis_pi0->GetYaxis()->SetTitle("N(Events)");

    g2_evis_muon = new TH1D("g2_evis_muon","Visible Energy",binList.gamma_evis_pdg.get_nBins(), binList.gamma_evis_pdg.get_min(), binList.gamma_evis_pdg.get_max() );
    g2_evis_muon->GetXaxis()->SetTitle("Evis");
    g2_evis_muon->GetYaxis()->SetTitle("N(Events)");
   
    g3_evis_proton = new TH1D("g3_evis_proton","Visible Energy",binList.pi0_evis_pdg.get_nBins(), binList.pi0_evis_pdg.get_min(), binList.pi0_evis_pdg.get_max() );
    g3_evis_proton->GetXaxis()->SetTitle("Evis");
    g3_evis_proton->GetYaxis()->SetTitle("N(Events)");

    g3_evis_neutron = new TH1D("g3_evis_neutron","Visible Energy",binList.pi0_evis_pdg.get_nBins(), binList.pi0_evis_pdg.get_min(), binList.pi0_evis_pdg.get_max() );
    g3_evis_neutron->GetXaxis()->SetTitle("Evis");
    g3_evis_neutron->GetYaxis()->SetTitle("N(Events)");

    g3_evis_pi = new TH1D("g3_evis_pi","Visible Energy",binList.pi0_evis_pdg.get_nBins(), binList.pi0_evis_pdg.get_min(), binList.pi0_evis_pdg.get_max() );
    g3_evis_pi->GetXaxis()->SetTitle("Evis");
    g3_evis_pi->GetYaxis()->SetTitle("N(Events)");

    g3_evis_pi0 = new TH1D("g3_evis_pi0","Visible Energy",binList.pi0_evis_pdg.get_nBins(), binList.pi0_evis_pdg.get_min(), binList.pi0_evis_pdg.get_max() );
    g3_evis_pi0->GetXaxis()->SetTitle("Evis");
    g3_evis_pi0->GetYaxis()->SetTitle("N(Events)");

    g3_evis_muon = new TH1D("g3_evis_muon","Visible Energy",binList.pi0_evis_pdg.get_nBins(), binList.pi0_evis_pdg.get_min(), binList.pi0_evis_pdg.get_max() );
    g3_evis_muon->GetXaxis()->SetTitle("Evis");
    g3_evis_muon->GetYaxis()->SetTitle("N(Events)");
}

void CCProtonPi0_Pion::writeHistograms()
{
    std::cout<<">> Writing "<<rootDir<<std::endl;
    f->cd();

    for (int i = 0; i < nHistograms; i++){
        // Unique Histograms
        invMass[i]->Write();
        photonEnergy_Asymmetry[i]->Write();
        
        evis_frac_reco_pi0_true_pi0[i]->Write();
        evis_frac_true_pi0_reco_all[i]->Write();
        evis_frac_reco_pi0_reco_all[i]->Write();
        evis_frac_reco_nonpi0_reco_all[i]->Write();

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

    // Photon Comparison
    gamma1_E_gamma2_E->Write();
    gamma1_convLength_gamma2_convLength->Write();

    // Truth Match
    g1_evis_proton->Write();
    g1_evis_neutron->Write();
    g1_evis_pi->Write();
    g1_evis_pi0->Write();
    g1_evis_muon->Write();

    g2_evis_proton->Write();
    g2_evis_neutron->Write();
    g2_evis_pi->Write();
    g2_evis_pi0->Write();
    g2_evis_muon->Write();

    g3_evis_proton->Write();
    g3_evis_neutron->Write();
    g3_evis_pi->Write();
    g3_evis_pi0->Write();
    g3_evis_muon->Write();

    // Gamma1 True Energy
    gamma1_true_E->Write();
    gamma1_reco_error_E->Write();
    gamma1_reco_E_true_E->Write();
    gamma1_true_E_reco_E_ratio->Write();
  
    // Gamma2 True Energy
    gamma2_true_E->Write();
    gamma2_reco_error_E->Write();
    gamma2_reco_E_true_E->Write();
    gamma2_true_E_reco_E_ratio->Write();
    
    f->Close();
}
#endif
