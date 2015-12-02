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
        bin_invMass.setBin(20,50,250.0);
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

        temp = new MnvH1D( Form("%s_%d","gamma1_E_Old",i),"Leading Photon Energy Old",bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max() );
        temp->GetXaxis()->SetTitle("E_{#gamma_{1}} [GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f [GeV]",bin_photonP.get_width()));
        gamma1_E_Old.push_back(temp);

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

        temp = new MnvH1D( Form("%s_%d","gamma2_E_Old",i),"Secondary Photon Energy Old",bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max() );
        temp->GetXaxis()->SetTitle("E_{#gamma_{2}} [GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f [GeV]",bin_photonP.get_width()));
        gamma2_E_Old.push_back(temp);

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

        temp = new MnvH1D( Form("%s_%d","invMass_Old",i),"Reconstructed Pi0 Invariant Mass Old",bin_invMass.get_nBins(), bin_invMass.get_min(), bin_invMass.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed m_{#gamma#gamma} [MeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f [MeV]",bin_invMass.get_width()));   
        invMass_Old.push_back(temp);

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
 
    gamma1_reco_Old_error_E = new TH1D( "gamma1_reco_Old_error_E","Error on Reconstructed Energy Old",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    gamma1_reco_Old_error_E->GetXaxis()->SetTitle("(E_{Reco}-E_{True})/E_{True}");
    gamma1_reco_Old_error_E->GetYaxis()->SetTitle("N(Events)");

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
  
    gamma2_reco_Old_error_E = new TH1D( "gamma2_reco_Old_error_E","Error on Reconstructed Energy Old",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    gamma2_reco_Old_error_E->GetXaxis()->SetTitle("(E_{Reco}-E_{True})/E_{True}");
    gamma2_reco_Old_error_E->GetYaxis()->SetTitle("N(Events)");
 
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
     
    gamma1_E_gamma2_E = new TH2D( "gamma1_E_gamma2_E","Leading Photon vs Secondary Photon Energy",bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max(),bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max() );
    gamma1_E_gamma2_E->GetXaxis()->SetTitle("Reconstructed E_{#gamma_{1}} [GeV]");
    gamma1_E_gamma2_E->GetYaxis()->SetTitle("Reconstructed E_{#gamma_{2}} [GeV]");

    reco_E_true_E = new TH2D( "reco_E_true_E","Photon True vs Reconstructed Energy",bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max(),bin_photonP.get_nBins(), bin_photonP.get_min(), bin_photonP.get_max());
    reco_E_true_E->GetXaxis()->SetTitle("Reconstructed Energy [GeV]");
    reco_E_true_E->GetYaxis()->SetTitle("True Energy [GeV]");

    true_E_evis_trkr_ratio = new TH2D( "true_E_evis_trkr_ratio","Showers Contained inside Tracker",40,0,600,40,1,2.0);
    true_E_evis_trkr_ratio->GetXaxis()->SetTitle("E_{True} [MeV]");
    true_E_evis_trkr_ratio->GetYaxis()->SetTitle("E_{True}/E_{Visible}");

    double bins[7] = {0,50,100,150,200,300,700};
    evis_evis_trkr_ratio = new TH2D( "evis_evis_trkr_ratio","Showers Contained inside Tracker",6,bins,10,1,2);
    evis_evis_trkr_ratio->GetXaxis()->SetTitle("E_{Visible} [MeV]");
    evis_evis_trkr_ratio->GetYaxis()->SetTitle("E_{True}/E_{Visible}");

    evis_evis_ecal_ratio = new TH2D( "evis_evis_ecal_ratio","Showers Contained inside ECAL",5,bins,10,1,6.0);
    evis_evis_ecal_ratio->GetXaxis()->SetTitle("E_{Visible} [MeV]");
    evis_evis_ecal_ratio->GetYaxis()->SetTitle("E_{True}/E_{Visible}");

    evis_trkr_ratio_1 = new TH1D("evis_trkr_ratio_1","0 MeV < Evis < 50 MeV",10,1,2);
    evis_trkr_ratio_1->GetXaxis()->SetTitle("E_{True}/E_{Visible}");
    evis_trkr_ratio_1->GetYaxis()->SetTitle("N(Events)");

    evis_trkr_ratio_2 = new TH1D("evis_trkr_ratio_2","50 MeV < Evis < 100 MeV",10,1,2);
    evis_trkr_ratio_2->GetXaxis()->SetTitle("E_{True}/E_{Visible}");
    evis_trkr_ratio_2->GetYaxis()->SetTitle("N(Events)");

    evis_trkr_ratio_3 = new TH1D("evis_trkr_ratio_3","100 MeV < Evis < 150 MeV",10,1,2);
    evis_trkr_ratio_3->GetXaxis()->SetTitle("E_{True}/E_{Visible}");
    evis_trkr_ratio_3->GetYaxis()->SetTitle("N(Events)");

    evis_trkr_ratio_4 = new TH1D("evis_trkr_ratio_4","150 MeV < Evis < 200 MeV",10,1,2);
    evis_trkr_ratio_4->GetXaxis()->SetTitle("E_{True}/E_{Visible}");
    evis_trkr_ratio_4->GetYaxis()->SetTitle("N(Events)");

    evis_trkr_ratio_5 = new TH1D("evis_trkr_ratio_5","200 MeV < Evis < 300 MeV",10,1,2);
    evis_trkr_ratio_5->GetXaxis()->SetTitle("E_{True}/E_{Visible}");
    evis_trkr_ratio_5->GetYaxis()->SetTitle("N(Events)");

    evis_trkr_ratio_6 = new TH1D("evis_trkr_ratio_6","300 MeV < Evis < 700 MeV",10,1,2);
    evis_trkr_ratio_6->GetXaxis()->SetTitle("E_{True}/E_{Visible}");
    evis_trkr_ratio_6->GetYaxis()->SetTitle("N(Events)");


    evis_ecal_ratio_1 = new TH1D("evis_ecal_ratio_1","0 MeV < Evis < 50 MeV",10,1,6);
    evis_ecal_ratio_1->GetXaxis()->SetTitle("E_{True}/E_{Visible}");
    evis_ecal_ratio_1->GetYaxis()->SetTitle("N(Events)");

    evis_ecal_ratio_2 = new TH1D("evis_ecal_ratio_2","50 MeV < Evis < 100 MeV",10,1,6);
    evis_ecal_ratio_2->GetXaxis()->SetTitle("E_{True}/E_{Visible}");
    evis_ecal_ratio_2->GetYaxis()->SetTitle("N(Events)");

    evis_ecal_ratio_3 = new TH1D("evis_ecal_ratio_3","100 MeV < Evis < 150 MeV",10,1,6);
    evis_ecal_ratio_3->GetXaxis()->SetTitle("E_{True}/E_{Visible}");
    evis_ecal_ratio_3->GetYaxis()->SetTitle("N(Events)");

    evis_ecal_ratio_4 = new TH1D("evis_ecal_ratio_4","150 MeV < Evis < 200 MeV",10,1,6);
    evis_ecal_ratio_4->GetXaxis()->SetTitle("E_{True}/E_{Visible}");
    evis_ecal_ratio_4->GetYaxis()->SetTitle("N(Events)");

    evis_ecal_ratio_5 = new TH1D("evis_ecal_ratio_5","200 MeV < Evis < 300 MeV",10,1,6);
    evis_ecal_ratio_5->GetXaxis()->SetTitle("E_{True}/E_{Visible}");
    evis_ecal_ratio_5->GetYaxis()->SetTitle("N(Events)");

    reco_error_trkr = new TH1D("reco_error_trkr","Error on Reconstructed Energy (Tracker Contained)",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    reco_error_trkr->GetXaxis()->SetTitle("(E_{Reco}-E_{True})/E_{True}");
    reco_error_trkr->GetYaxis()->SetTitle("N(Events)");
 
    reco_error_ecal = new TH1D("reco_error_ecal","Error on Reconstructed Energy (ECAL Contained)",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    reco_error_ecal->GetXaxis()->SetTitle("(E_{Reco}-E_{True})/E_{True}");
    reco_error_ecal->GetYaxis()->SetTitle("N(Events)");

    reco_error_trkr_ecal = new TH1D("reco_error_trkr_ecal","Error on Reconstructed Energy (Tracker + ECAL)",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    reco_error_trkr_ecal->GetXaxis()->SetTitle("(E_{Reco}-E_{True})/E_{True}");
    reco_error_trkr_ecal->GetYaxis()->SetTitle("N(Events)");
 
    calc_error_trkr = new TH1D("calc_error_trkr","Error on Corrected Energy (Tracker Contained)",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    calc_error_trkr->GetXaxis()->SetTitle("(E_{Reco}-E_{True})/E_{True}");
    calc_error_trkr->GetYaxis()->SetTitle("N(Events)");
 
    calc_error_ecal = new TH1D("calc_error_ecal","Error on Corrected Energy (ECAL Contained)",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    calc_error_ecal->GetXaxis()->SetTitle("(E_{Reco}-E_{True})/E_{True}");
    calc_error_ecal->GetYaxis()->SetTitle("N(Events)");

    calc_error_trkr_ecal = new TH1D("calc_error_trkr_ecal","Error on Reconstructed Energy (Tracker + ECAL)",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    calc_error_trkr_ecal->GetXaxis()->SetTitle("(E_{Reco}-E_{True})/E_{True}");
    calc_error_trkr_ecal->GetYaxis()->SetTitle("N(Events)");

    gamma1_calc_error = new TH1D( "gamma1_calc_error","Error on Corrected Energy Gamma 1",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    gamma1_calc_error->GetXaxis()->SetTitle("(E_{Reco}-E_{True})/E_{True}");
    gamma1_calc_error->GetYaxis()->SetTitle("N(Events)");
 
    gamma2_calc_error = new TH1D( "gamma2_calc_error","Error on Corrected Energy Gamma 2",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    gamma2_calc_error->GetXaxis()->SetTitle("(E_{Reco}-E_{True})/E_{True}");
    gamma2_calc_error->GetYaxis()->SetTitle("N(Events)");

    mgg_calc = new TH1D("mgg_calc","Corrected Pi0 Invariant Mass",bin_invMass.get_nBins(), bin_invMass.get_min(), bin_invMass.get_max() );
    mgg_calc->GetXaxis()->SetTitle("Reconstructed m_{#gamma#gamma} [MeV]");
    mgg_calc->GetYaxis()->SetTitle(Form("Events / %3.2f [MeV]",bin_invMass.get_width()));   

    gamma_nPlanes = new TH1D("gamma_nPlanes","Shower Contained inside ECAL",20,0,20);
    gamma_nPlanes->GetXaxis()->SetTitle("N(Planes)");
    gamma_nPlanes->GetYaxis()->SetTitle("Events");   

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
 
    Shower_Topology = new TH1D( "Shower_Topology","Shower Topology",binList.multiplicity.get_nBins(),binList.multiplicity.get_min(),binList.multiplicity.get_max());
    Shower_Topology->GetXaxis()->SetTitle("0:TRKR, 1:ECAL, 2:SCAL, 3:TRKR+SCAL, 4:TRKR+ECAL, 5:ECAL+SCAL, 6:Other");
    Shower_Topology->GetYaxis()->SetTitle("N(Showers)");
 
    energy_trkr = new TH1D("energy_trkr","Reconstructed Energy in Tracker",20,0,500);
    energy_trkr->GetXaxis()->SetTitle("Reconstructed Energy [MeV]");
    energy_trkr->GetYaxis()->SetTitle("N(Events)");

    energy_ecal = new TH1D("energy_ecal","Reconstructed Energy in ECAL",20,0,500);
    energy_ecal->GetXaxis()->SetTitle("Reconstructed Energy [MeV]");
    energy_ecal->GetYaxis()->SetTitle("N(Events)");

    energy_hcal = new TH1D("energy_hcal","Reconstructed Energy in HCAL",20,0,500);
    energy_hcal->GetXaxis()->SetTitle("Reconstructed Energy [MeV]");
    energy_hcal->GetYaxis()->SetTitle("N(Events)");

    energy_scal = new TH1D("energy_scal","Reconstructed Energy in Side ECAL",20,0,500);
    energy_scal->GetXaxis()->SetTitle("Reconstructed Energy [MeV]");
    energy_scal->GetYaxis()->SetTitle("N(Events)");

    energy_frac_trkr = new TH1D("energy_frac_trkr","Reconstructed Energy Fraction in Tracker",22,0,1.1);
    energy_frac_trkr->GetXaxis()->SetTitle("E(Tracker)/E(Total)");
    energy_frac_trkr->GetYaxis()->SetTitle("N(Events)");
 
    energy_frac_ecal = new TH1D("energy_frac_ecal","Reconstructed Energy Fraction in ECAL",22,0,1.1);
    energy_frac_ecal->GetXaxis()->SetTitle("E(ECAL)/E(Total)");
    energy_frac_ecal->GetYaxis()->SetTitle("N(Events)");
 
    energy_frac_hcal = new TH1D("energy_frac_hcal","Reconstructed Energy Fraction in HCAL",22,0,1.1);
    energy_frac_hcal->GetXaxis()->SetTitle("E(HCAL)/E(Total)");
    energy_frac_hcal->GetYaxis()->SetTitle("N(Events)");
 
    energy_frac_scal = new TH1D("energy_frac_scal","Reconstructed Energy Fraction in Side ECAL",22,0,1.1);
    energy_frac_scal->GetXaxis()->SetTitle("E(Side ECAL)/E(Total)");
    energy_frac_scal->GetYaxis()->SetTitle("N(Events)");

    evis_trkr_reco_true = new TH2D( "evis_trkr_reco_true","True vs Reco Tracker Visible Energy",20,5,150,20,5,150);
    evis_trkr_reco_true->GetXaxis()->SetTitle("E_{Visible}^{Reco} [MeV]");
    evis_trkr_reco_true->GetYaxis()->SetTitle("E_{Visible}^{True} [MeV]");

    evis_ecal_reco_true = new TH2D( "evis_ecal_reco_true","True vs Reco ECAL Visible Energy",20,5,50,20,5,50);
    evis_ecal_reco_true->GetXaxis()->SetTitle("E_{Visible}^{Reco} [MeV]");
    evis_ecal_reco_true->GetYaxis()->SetTitle("E_{Visible}^{True} [MeV]");

    evis_scal_reco_true = new TH2D( "evis_scal_reco_true","True vs Reco Side ECAL Visible Energy",20,5,50,20,5,50);
    evis_scal_reco_true->GetXaxis()->SetTitle("E_{Visible}^{Reco} [MeV]");
    evis_scal_reco_true->GetYaxis()->SetTitle("E_{Visible}^{True} [MeV]");

    evis_hcal_reco_true = new TH2D( "evis_hcal_reco_true","True vs Reco HCAL Visible Energy",20,5,50,20,5,50);
    evis_hcal_reco_true->GetXaxis()->SetTitle("E_{Visible}^{Reco} [MeV]");
    evis_hcal_reco_true->GetYaxis()->SetTitle("E_{Visible}^{True} [MeV]");
 
    energy_trkr_reco_true = new TH2D( "energy_trkr_reco_true","True vs Reco Tracker Energy",20,5,250,20,5,250);
    energy_trkr_reco_true->GetXaxis()->SetTitle("E_{Reco} [MeV]");
    energy_trkr_reco_true->GetYaxis()->SetTitle("E_{True} [MeV]");

    energy_ecal_reco_true = new TH2D( "energy_ecal_reco_true","True vs Reco ECAL Energy",20,5,50,20,5,50);
    energy_ecal_reco_true->GetXaxis()->SetTitle("E_{Reco} [MeV]");
    energy_ecal_reco_true->GetYaxis()->SetTitle("E_{True} [MeV]");

    energy_scal_reco_true = new TH2D( "energy_scal_reco_true","True vs Reco Side ECAL Energy",20,5,50,20,5,50);
    energy_scal_reco_true->GetXaxis()->SetTitle("E_{Reco} [MeV]");
    energy_scal_reco_true->GetYaxis()->SetTitle("E_{True} [MeV]");

    energy_hcal_reco_true = new TH2D( "energy_hcal_reco_true","True vs Reco HCAL Energy",20,5,50,20,5,50);
    energy_hcal_reco_true->GetXaxis()->SetTitle("E_{Reco} [MeV]");
    energy_hcal_reco_true->GetYaxis()->SetTitle("E_{True} [MeV]");
 
    // Side ECAL nHits Study
    trkr_nHits_reco_correct = new TH2D( "trkr_nHits_reco_correct","Correct vs All Tracker N(Hits)",100,0,100,100,0,100);
    trkr_nHits_reco_correct->GetXaxis()->SetTitle("Tracker N(Hits) based on Strip Number");
    trkr_nHits_reco_correct->GetYaxis()->SetTitle("Correctly assigned Tracker N(Hits)");
 
    trkr_nHits_reco_true = new TH2D( "trkr_nHits_reco_true","True vs Reco Tracker N(Hits)",100,0,100,100,0,100);
    trkr_nHits_reco_true->GetXaxis()->SetTitle("Tracker N(Hits) based on Strip Number");
    trkr_nHits_reco_true->GetYaxis()->SetTitle("Tracker N(Hits) based on TRUE Hit Position");
 
    scal_nHits_reco_correct = new TH2D( "scal_nHits_reco_correct","Correct vs All SCAL N(Hits)",100,0,100,100,0,100);
    scal_nHits_reco_correct->GetXaxis()->SetTitle("SCAL N(Hits) based on Strip Number");
    scal_nHits_reco_correct->GetYaxis()->SetTitle("Correctly assigned SCAL N(Hits)");
 
    scal_nHits_reco_true = new TH2D( "scal_nHits_reco_true","True vs Reco SCAL N(Hits)",100,0,100,100,0,100);
    scal_nHits_reco_true->GetXaxis()->SetTitle("SCAL N(Hits) based on Strip Number");
    scal_nHits_reco_true->GetYaxis()->SetTitle("SCAL N(Hits) based on TRUE Hit Position");
 
    trkr_nHits_2_reco_true = new TH2D( "trkr_nHits_2_reco_true","True vs Reco Tracker N(Hits) with Geometrical Method",100,0,100,100,0,100);
    trkr_nHits_2_reco_true->GetXaxis()->SetTitle("Tracker N(Hits) based Geometrical Assumption");
    trkr_nHits_2_reco_true->GetYaxis()->SetTitle("True Tracker N(Hits)");
  
    scal_nHits_2_reco_true = new TH2D( "scal_nHits_2_reco_true","True vs Reco SCAL N(Hits) with Geometrical Method",100,0,100,100,0,100);
    scal_nHits_2_reco_true->GetXaxis()->SetTitle("SCAL N(Hits) based Geometrical Assumption");
    scal_nHits_2_reco_true->GetYaxis()->SetTitle("True SCAL N(Hits)");
 
    old_model_trkr_nHits_error = new TH1D( "old_model_trkr_nHits_error","Old Model Tracker nHits Error",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    old_model_trkr_nHits_error->GetXaxis()->SetTitle("(N_{Reco}-N_{True})/N_{True}");
    old_model_trkr_nHits_error->GetYaxis()->SetTitle("N(Events)");
 
    new_model_trkr_nHits_error = new TH1D( "new_model_trkr_nHits_error","New Model Tracker nHits Error",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    new_model_trkr_nHits_error->GetXaxis()->SetTitle("(N_{Reco}-N_{True})/N_{True}");
    new_model_trkr_nHits_error->GetYaxis()->SetTitle("N(Events)");

    improved_model_trkr_nHits_error = new TH1D( "improved_model_trkr_nHits_error","Improved Model Tracker nHits Error",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    improved_model_trkr_nHits_error->GetXaxis()->SetTitle("(N_{Reco}-N_{True})/N_{True}");
    improved_model_trkr_nHits_error->GetYaxis()->SetTitle("N(Events)");

    old_model_scal_nHits_error = new TH1D( "old_model_scal_nHits_error","Old Model Side ECAL nHits Error",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    old_model_scal_nHits_error->GetXaxis()->SetTitle("(N_{Reco}-N_{True})/N_{True}");
    old_model_scal_nHits_error->GetYaxis()->SetTitle("N(Events)");
 
    new_model_scal_nHits_error = new TH1D( "new_model_scal_nHits_error","New Model Side ECAL nHits Error",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    new_model_scal_nHits_error->GetXaxis()->SetTitle("(N_{Reco}-N_{True})/N_{True}");
    new_model_scal_nHits_error->GetYaxis()->SetTitle("N(Events)");
 
    improved_model_scal_nHits_error = new TH1D( "improved_model_scal_nHits_error","Improved Model Side ECAL nHits Error",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    improved_model_scal_nHits_error->GetXaxis()->SetTitle("(N_{Reco}-N_{True})/N_{True}");
    improved_model_scal_nHits_error->GetYaxis()->SetTitle("N(Events)");

    scal_minZ_evis = new TH1D( "scal_minZ_evis","Visible Cluster Energy at min Z",30, 0, 30 );
    scal_minZ_evis->GetXaxis()->SetTitle("E_{Visible} [MeV]");
    scal_minZ_evis->GetYaxis()->SetTitle("N(Events)");

    scal_minZ_nDigits = new TH1D( "scal_minZ_nDigits","N(Digits) at min Z",10, 0, 10 );
    scal_minZ_nDigits->GetXaxis()->SetTitle("N(Digits)");
    scal_minZ_nDigits->GetYaxis()->SetTitle("N(Events)");


}

void CCProtonPi0_Pion::writeHistograms()
{
    std::cout<<">> Writing "<<rootDir<<std::endl;
    f->cd();

    for (int i = 0; i < nHistograms; i++){
        // Unique Histograms
        invMass[i]->Write();
        invMass_Old[i]->Write();
        photonEnergy_Asymmetry[i]->Write();
        
        evis_frac_reco_pi0_true_pi0[i]->Write();
        evis_frac_true_pi0_reco_all[i]->Write();
        evis_frac_reco_pi0_reco_all[i]->Write();
        evis_frac_reco_nonpi0_reco_all[i]->Write();

        // Leading Photon
        gamma1_ConvLength[i]->Write();
        gamma1_E[i]->Write();
        gamma1_E_Old[i]->Write();
        gamma1_theta[i]->Write();
       
        // Secondary Photon
        gamma2_ConvLength[i]->Write();
        gamma2_E[i]->Write();
        gamma2_E_Old[i]->Write();
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
    gamma1_reco_Old_error_E->Write();
    gamma1_reco_E_true_E->Write();
    gamma1_true_E_reco_E_error->Write();
  
    // Gamma2 True Energy
    gamma2_true_E->Write();
    gamma2_reco_Old_error_E->Write();
    gamma2_reco_error_E->Write();
    gamma2_reco_E_true_E->Write();
    gamma2_true_E_reco_E_error->Write();

    reco_E_true_E->Write();
    
    true_E_evis_trkr_ratio->Write();
    evis_evis_trkr_ratio->Write();
    evis_evis_ecal_ratio->Write();
    reco_error_trkr->Write();
    reco_error_ecal->Write();
    reco_error_trkr_ecal->Write();
    calc_error_trkr->Write();
    calc_error_ecal->Write();
    calc_error_trkr_ecal->Write();

    evis_trkr_ratio_1->Write();
    evis_trkr_ratio_2->Write();
    evis_trkr_ratio_3->Write();
    evis_trkr_ratio_4->Write();
    evis_trkr_ratio_5->Write();
    evis_trkr_ratio_6->Write();
 
    evis_ecal_ratio_1->Write();
    evis_ecal_ratio_2->Write();
    evis_ecal_ratio_3->Write();
    evis_ecal_ratio_4->Write();
    evis_ecal_ratio_5->Write();
   
    gamma1_calc_error->Write();
    gamma2_calc_error->Write();
    mgg_calc->Write();
    gamma_nPlanes->Write();

    Shower_Topology->Write();
    
    energy_trkr->Write();
    energy_ecal->Write();
    energy_hcal->Write();
    energy_scal->Write();

    energy_frac_trkr->Write();
    energy_frac_ecal->Write();
    energy_frac_hcal->Write();
    energy_frac_scal->Write();
    
    evis_trkr_reco_true->Write();
    evis_ecal_reco_true->Write();
    evis_scal_reco_true->Write();
    evis_hcal_reco_true->Write();
 
    energy_trkr_reco_true->Write();
    energy_ecal_reco_true->Write();
    energy_scal_reco_true->Write();
    energy_hcal_reco_true->Write();
   
    trkr_nHits_reco_correct->Write();
    trkr_nHits_reco_true->Write();
    scal_nHits_reco_correct->Write();
    scal_nHits_reco_true->Write();

    trkr_nHits_2_reco_true->Write();
    scal_nHits_2_reco_true->Write();

    old_model_trkr_nHits_error->Write();
    new_model_trkr_nHits_error->Write();
    improved_model_trkr_nHits_error->Write();

    old_model_scal_nHits_error->Write();
    new_model_scal_nHits_error->Write();
    improved_model_scal_nHits_error->Write();

    scal_minZ_evis->Write();
    scal_minZ_nDigits->Write();

    f->Close();
}
#endif
