/*
    See CCProtonPi0_Interaction.h header for Class Information
*/
#ifndef CCProtonPi0_Interaction_cpp
#define CCProtonPi0_Interaction_cpp

#include "CCProtonPi0_Interaction.h"

using namespace PlotUtils;

CCProtonPi0_Interaction::CCProtonPi0_Interaction(bool isModeReduce, bool isMC, std::string ana_folder) : CCProtonPi0_NTupleAnalysis()
{
    std::cout<<"Initializing CCProtonPi0_Interaction"<<std::endl;
    
    if(isModeReduce){
        std::cout<<"\tNTuple Reduce Mode -- Will not create ROOT Files"<<std::endl;
    }else{
        if (isMC) rootDir = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + ana_folder + "Interaction.root";
        else rootDir = Folder_List::rootOut + Folder_List::Data + Folder_List::analyzed  + ana_folder + "Interaction.root";
        
        std::cout<<"\tRoot File: "<<rootDir<<std::endl;
 
        // Create Root File 
        f = new TFile(rootDir.c_str(),"RECREATE");
        
        initHistograms();
    }
    
    std::cout<<"Done!"<<std::endl;
}


void CCProtonPi0_Interaction::initHistograms()
{
    MnvH1D* temp = NULL;

    for (int i = 0; i < nHistograms; i++){
        temp = new MnvH1D( Form("%s_%d","Enu_True",i),"True Beam Energy",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
        temp->GetXaxis()->SetTitle("True E_{#nu}[GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.beamE.get_width()));
        Enu_True.push_back(temp);
        
        temp = new MnvH1D( Form("%s_%d","Enu_1Track",i),"Reconstructed Beam Energy - 1 Track",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed E_{#nu} - 1Track [GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.beamE.get_width()));
        Enu_1Track.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","Enu_1Track_Alt",i),"Reconstructed Beam Energy - 1 Track Alternative",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed E_{#nu} - 1Track Alternative [GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.beamE.get_width()));
        Enu_1Track_Alt.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","Enu_2Track",i),"Reconstructed Beam Energy - 2 Track",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed E_{#nu} - 2 Track [GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.beamE.get_width()));
        Enu_2Track.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","QSq",i),"Reconstructed Q^{2}",binList.q2.get_nBins(), binList.q2.get_min(), binList.q2.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed Q^{2} [GeV^{2}]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.q2.get_width()));
        QSq.push_back(temp);    

        temp = new MnvH1D( Form("%s_%d","WSq",i),"Reconstructed W^{2}",binList.wSq.get_nBins(), binList.wSq.get_min(), binList.wSq.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed W^{2} [GeV^{2}]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.wSq.get_width()));
        WSq.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","W",i),"Reconstructed W",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed W [GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.w.get_width()));
        W.push_back(temp);
 
        temp = new MnvH1D( Form("%s_%d","vertex_energy_1Track",i),"Vertex Blob Energy (r = 90mm)",binList.vertex_energy.get_nBins(), binList.vertex_energy.get_min(), binList.vertex_energy.get_max() );
        temp->GetXaxis()->SetTitle("Vertex Blob Energy [MeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.vertex_energy.get_width()));
        vertex_energy_1Track.push_back(temp);
  
        temp = new MnvH1D( Form("%s_%d","vertex_energy_2Track",i),"Vertex Blob Energy (r = 90mm)",binList.vertex_energy.get_nBins(), binList.vertex_energy.get_min(), binList.vertex_energy.get_max() );
        temp->GetXaxis()->SetTitle("Vertex Blob Energy [MeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.vertex_energy.get_width()));
        vertex_energy_2Track.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","vertex_energy_Corrected_1Track",i),"Vertex Blob Energy (r = 90mm)",binList.vertex_energy.get_nBins(), binList.vertex_energy.get_min(), binList.vertex_energy.get_max() );
        temp->GetXaxis()->SetTitle("Vertex Blob Energy [MeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.vertex_energy.get_width()));
        vertex_energy_Corrected_1Track.push_back(temp);
  
        temp = new MnvH1D( Form("%s_%d","vertex_energy_Corrected_2Track",i),"Vertex Blob Energy (r = 90mm)",binList.vertex_energy.get_nBins(), binList.vertex_energy.get_min(), binList.vertex_energy.get_max() );
        temp->GetXaxis()->SetTitle("Vertex Blob Energy [MeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.vertex_energy.get_width()));
        vertex_energy_Corrected_2Track.push_back(temp);


        temp = new MnvH1D( Form("%s_%d","vertex_evis_1Track",i),"Vertex Blob Visible Energy (r = 90mm)",binList.vertex_evis.get_nBins(), binList.vertex_evis.get_min(), binList.vertex_evis.get_max() );
        temp->GetXaxis()->SetTitle("Vertex Blob Visible Energy [MeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.vertex_evis.get_width()));
        vertex_evis_1Track.push_back(temp);
  
        temp = new MnvH1D( Form("%s_%d","vertex_evis_2Track",i),"Vertex Blob Visible Energy (r = 90mm)",binList.vertex_evis.get_nBins(), binList.vertex_evis.get_min(), binList.vertex_evis.get_max() );
        temp->GetXaxis()->SetTitle("Vertex Blob Visible Energy [MeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.vertex_evis.get_width()));
        vertex_evis_2Track.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","extra_energy_1Track",i),"Extra Energy (r = 300mm)",binList.extra_energy.get_nBins(), binList.extra_energy.get_min(), binList.extra_energy.get_max() );
        temp->GetXaxis()->SetTitle("Extra Energy [MeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.extra_energy.get_width()));
        extra_energy_1Track.push_back(temp);
 
        temp = new MnvH1D( Form("%s_%d","extra_energy_2Track",i),"Extra Energy (r = 300mm)",binList.extra_energy.get_nBins(), binList.extra_energy.get_min(), binList.extra_energy.get_max() );
        temp->GetXaxis()->SetTitle("Extra Energy [MeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.extra_energy.get_width()));
        extra_energy_2Track.push_back(temp);
      
        temp = new MnvH1D( Form("%s_%d","deltaInvMass",i),"Reconstructed Delta+ Invariant Mass",binList.deltaInvMass.get_nBins(), binList.deltaInvMass.get_min(), binList.deltaInvMass.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed #Delta^{+} Inv. Mass [MeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f [MeV] ",binList.deltaInvMass.get_width()));
        deltaInvMass.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","E_Unused_afterReco",i),"Unused Cluster Energy after Reconstruction",binList.UnusedE.get_nBins(), binList.UnusedE.get_min(), binList.UnusedE.get_max() );
        temp->GetXaxis()->SetTitle("Unused E_{Visible} after Reconstruction [MeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.UnusedE.get_width()));
        E_Unused_afterReco.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","E_Used_afterReco",i),"Used Cluster Energy after Reconstruction",binList.UsedE.get_nBins(), binList.UsedE.get_min(), binList.UsedE.get_max() );
        temp->GetXaxis()->SetTitle("Used E_{Visible} after Reconstruction [MeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.UsedE.get_width()));
        E_Used_afterReco.push_back(temp);

        // Vertex
        temp = new MnvH1D( Form("%s_%d","vertex_z",i),"Reconstructed Interaction Vertex",binList.vertex_z.get_nBins(), binList.vertex_z.get_min(), binList.vertex_z.get_max() );
        temp->GetXaxis()->SetTitle("z = 4293 Target, #bf{z = 5810 Interaction Region}, z = 8614 ECAL, z = 9088 HCAL");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.vertex_z.get_width()));
        vertex_z.push_back(temp);

        // --------------------------------------------------------------------
        // Visible Energy -- Neutrino Energy Study
        // --------------------------------------------------------------------
        // 1 Track
        temp = new MnvH1D( Form("%s_%d","evis_total_1Track",i),"Total Visible Energy -- 1Track",100,0.0,2000.0 );
        temp->GetXaxis()->SetTitle("Total Visible Energy [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        evis_total_1Track.push_back(temp);
 
        temp = new MnvH1D( Form("%s_%d","evis_muon_1Track",i),"Muon Visible Energy -- 1Track",100,0.0,1000.0 );
        temp->GetXaxis()->SetTitle("Muon Visible Energy [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        evis_muon_1Track.push_back(temp);
  
        temp = new MnvH1D( Form("%s_%d","evis_pi0_1Track",i),"Pi0 Visible Energy -- 1Track",100,0.0,1000.0 );
        temp->GetXaxis()->SetTitle("Pi0 Visible Energy [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        evis_pi0_1Track.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","evis_hadron_1Track",i),"Hadronic Visible Energy -- 1Track",100,0.0,1500.0 );
        temp->GetXaxis()->SetTitle("E_{Visible}^{Total} - E_{Visible}^{#mu} [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        evis_hadron_1Track.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","energy_hadron_true_1Track",i),"Hadronic True Energy -- 1Track",100,0.0,5000.0 );
        temp->GetXaxis()->SetTitle("E_{True}^{#nu} - E_{True}^{#mu} [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        energy_hadron_true_1Track.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","energy_hadron_reco_1Track",i),"Hadronic Reco Energy -- 1Track",100,0.0,5000.0 );
        temp->GetXaxis()->SetTitle("E_{Reco}^{#nu} - E_{Reco}^{#mu} [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        energy_hadron_reco_1Track.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","evis_extra_1Track",i),"Extra Visible Energy -- 1Track",100,0.0,1000.0 );
        temp->GetXaxis()->SetTitle("E_{Visible}^{Total} - (E_{Visible}^{#mu} + E_{Visible}^{#pi^{0}}) [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        evis_extra_1Track.push_back(temp);
  
        temp = new MnvH1D( Form("%s_%d","energy_extra_true_1Track",i),"Extra True Energy -- 1Track",50,0.0,500.0 );
        temp->GetXaxis()->SetTitle("(E_{True}^{#nu} + m_{n}) - (E_{True}^{#mu} + E_{True}^{#pi^{0}}) [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        energy_extra_true_1Track.push_back(temp);

        // 2 Track
        temp = new MnvH1D( Form("%s_%d","evis_total_2Track",i),"Total Visible Energy -- 2Track",100,0.0,2000.0 );
        temp->GetXaxis()->SetTitle("Total Visible Energy [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        evis_total_2Track.push_back(temp);
 
        temp = new MnvH1D( Form("%s_%d","evis_muon_2Track",i),"Muon Visible Energy -- 2Track",100,0.0,1000.0 );
        temp->GetXaxis()->SetTitle("Muon Visible Energy [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        evis_muon_2Track.push_back(temp);
  
        temp = new MnvH1D( Form("%s_%d","evis_pi0_2Track",i),"Pi0 Visible Energy -- 2Track",100,0.0,1000.0 );
        temp->GetXaxis()->SetTitle("Pi0 Visible Energy [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        evis_pi0_2Track.push_back(temp);
 
        temp = new MnvH1D( Form("%s_%d","evis_proton_2Track",i),"Proton Visible Energy -- 2Track",100,0.0,1000.0 );
        temp->GetXaxis()->SetTitle("Proton Visible Energy [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        evis_proton_2Track.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","evis_hadron_2Track",i),"Hadronic Visible Energy -- 2Track",100,0.0,1500.0 );
        temp->GetXaxis()->SetTitle("E_{Visible}^{Total} - E_{Visible}^{#mu} [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        evis_hadron_2Track.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","energy_hadron_true_2Track",i),"Hadronic True Energy -- 2Track",100,0.0,5000.0 );
        temp->GetXaxis()->SetTitle("E_{True}^{#nu} - E_{True}^{#mu} [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        energy_hadron_true_2Track.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","energy_hadron_reco_2Track",i),"Hadronic Reco Energy -- 2Track",100,0.0,5000.0 );
        temp->GetXaxis()->SetTitle("E_{Reco}^{#nu} - E_{Reco}^{#mu} [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        energy_hadron_reco_2Track.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","evis_hadron_nopi0_2Track",i),"Hadronic Visible Energy (#pi^0 E_{Visible} Subtracted) -- 2Track",100,0.0,1500.0 );
        temp->GetXaxis()->SetTitle("E_{Visible}^{Total} - (E_{Visible}^{#mu} + E_{Visible}^{#pi^{0}}) [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        evis_hadron_nopi0_2Track.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","energy_hadron_nopi0_true_2Track",i),"Hadronic True Energy (#pi^0 E_{True} Subtracted)-- 2Track",100,0.0,5000.0 );
        temp->GetXaxis()->SetTitle("E_{True}^{#nu} - (E_{True}^{#mu} + E_{True}^{#pi^{0}}) [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        energy_hadron_nopi0_true_2Track.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","evis_extra_2Track",i),"Extra Visible Energy -- 2Track",100,0.0,1000.0 );
        temp->GetXaxis()->SetTitle("E_{Visible}^{Total} - (E_{Visible}^{#mu} + E_{Visible}^{#pi^{0}} + E_{Visible}^{p}) [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        evis_extra_2Track.push_back(temp);
  
        temp = new MnvH1D( Form("%s_%d","energy_extra_true_2Track",i),"Extra True Energy -- 2Track",50,0.0,500.0 );
        temp->GetXaxis()->SetTitle("(E_{True}^{#nu} + m_{n}) - (E_{True}^{#mu} + E_{True}^{#pi^{0}} + E_{True}^{p}) [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        energy_extra_true_2Track.push_back(temp);
    }
    
    // MC Only Histograms
    final_mc_w_DIS = new TH1D( "final_mc_w_DIS","True W for DIS",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max() );
    final_mc_w_DIS->GetXaxis()->SetTitle("True W for DIS [GeV]");
    final_mc_w_DIS->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.mc_w.get_width()));
    
    final_mc_w_RES = new TH1D( "final_mc_w_RES","True W for RES",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max() );
    final_mc_w_RES->GetXaxis()->SetTitle("True W for RES [GeV]");
    final_mc_w_RES->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.mc_w.get_width()));
    
    final_mc_w_CCQE = new TH1D( "final_mc_w_CCQE","True W for CCQE",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max() );
    final_mc_w_CCQE->GetXaxis()->SetTitle("True W for CCQE [GeV]");
    final_mc_w_CCQE->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.mc_w.get_width()));

    // Short Proton True Information
    proton_true_P_1Track = new TH1D( "proton_true_P_1Track","Short Proton True Momentum",binList.short_proton_P.get_nBins(), binList.short_proton_P.get_min(), binList.short_proton_P.get_max() );
    proton_true_P_1Track->GetXaxis()->SetTitle("Short Proton P_{True} [MeV]");
    proton_true_P_1Track->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.short_proton_P.get_width()));

    proton_true_KE_1Track = new TH1D( "proton_true_KE_1Track","Short Proton True Momentum",binList.short_proton_KE.get_nBins(), binList.short_proton_KE.get_min(), binList.short_proton_KE.get_max() );
    proton_true_KE_1Track->GetXaxis()->SetTitle("Short Proton P_{True} [MeV]");
    proton_true_KE_1Track->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.short_proton_KE.get_width()));

    // Neutrino Energy Error
    Enu_1Track_Error = new TH1D("Enu_1Track_Error","Neutrino Energy Error - 1 Track",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    Enu_1Track_Error->GetXaxis()->SetTitle("(E_{#nu}^{Reco}-E_{#nu}^{True})/E_{#nu}^{True}");
    Enu_1Track_Error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));

    Enu_1Track_Alt_Error = new TH1D("Enu_1Track_Alt_Error","Neutrino Energy Error - 1 Track Alternative",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    Enu_1Track_Alt_Error->GetXaxis()->SetTitle("(E_{#nu}^{Reco}-E_{#nu}^{True})/E_{#nu}^{True}");
    Enu_1Track_Alt_Error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));

    Enu_2Track_Error = new TH1D("Enu_2Track_Error","Neutrino Energy Error - 2 Track",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    Enu_2Track_Error->GetXaxis()->SetTitle("(E_{#nu}^{Reco}-E_{#nu}^{True})/E_{#nu}^{True}");
    Enu_2Track_Error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));

    Enu_1Track_Corrected_Error = new TH1D("Enu_1Track_Corrected_Error","Neutrino Energy Error - 1 Track Corrected",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    Enu_1Track_Corrected_Error->GetXaxis()->SetTitle("(E_{#nu}^{Reco}-E_{#nu}^{True})/E_{#nu}^{True}");
    Enu_1Track_Corrected_Error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));

    Enu_2Track_Corrected_Error = new TH1D("Enu_2Track_Corrected_Error","Neutrino Energy Error - 2 Track Corrected",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    Enu_2Track_Corrected_Error->GetXaxis()->SetTitle("(E_{#nu}^{Reco}-E_{#nu}^{True})/E_{#nu}^{True}");
    Enu_2Track_Corrected_Error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));

    Enu_1Track_True = new TH2D("Enu_1Track_True","Neutrino Energy True vs 1 Track",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max(),binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    Enu_1Track_True->GetXaxis()->SetTitle("E_{#nu}^{Reco} [GeV]");
    Enu_1Track_True->GetYaxis()->SetTitle("E_{#nu}^{True} [GeV]");

    Enu_1Track_Alt_True = new TH2D("Enu_1Track_Alt_True","Neutrino Energy True vs 1 Track Alternative",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max(),binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    Enu_1Track_Alt_True->GetXaxis()->SetTitle("E_{#nu}^{Reco} [GeV]");
    Enu_1Track_Alt_True->GetYaxis()->SetTitle("E_{#nu}^{True} [GeV]");

    Enu_1Track_1Track_Alt = new TH2D("Enu_1Track_1Track_Alt","Neutrino Energy 1Track vs 1 Track Alternative",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max(),binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    Enu_1Track_1Track_Alt->GetXaxis()->SetTitle("E_{#nu}^{Reco} [GeV]");
    Enu_1Track_1Track_Alt->GetYaxis()->SetTitle("E_{#nu}^{Reco} Alternative [GeV]");

    Enu_2Track_True = new TH2D("Enu_2Track_True","Neutrino Energy True vs 2 Track",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max(),binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    Enu_2Track_True->GetXaxis()->SetTitle("E_{#nu}^{Reco} [GeV]");
    Enu_2Track_True->GetYaxis()->SetTitle("E_{#nu}^{True} [GeV]");

    // Neutrino Energy Study
    evis_extra_energy_extra_true_1Track = new TH2D("evis_extra_energy_extra_true_1Track","Extra Energy True vs Visible Energy -- 1Track",100,0,1000,100,0,1000);
    evis_extra_energy_extra_true_1Track->GetXaxis()->SetTitle("E_{Visible}^{Total} - (E_{Visible}^{#mu} + E_{Visible}^{#pi^{0}} + E_{Visible}^{p}) [MeV]");
    evis_extra_energy_extra_true_1Track->GetYaxis()->SetTitle("E_{True}^{Total} - (E_{True}^{#mu} + E_{True}^{#pi^{0}} + E_{True}^{p}) [MeV]");

    evis_extra_energy_extra_true_2Track = new TH2D("evis_extra_energy_extra_true_2Track","Extra Energy True vs Visible Energy -- 2Track",100,0,1000,100,0,1000);
    evis_extra_energy_extra_true_2Track->GetXaxis()->SetTitle("E_{Visible}^{Total} - (E_{Visible}^{#mu} + E_{Visible}^{#pi^{0}} + E_{Visible}^{p}) [MeV]");
    evis_extra_energy_extra_true_2Track->GetYaxis()->SetTitle("E_{True}^{Total} - (E_{True}^{#mu} + E_{True}^{#pi^{0}} + E_{True}^{p}) [MeV]");

    energy_hadron_reco_true_1Track = new TH2D("energy_hadron_reco_true_1Track","Hadronic Energy Truve vs Reco -- 1Track",100,0,5000,100,0,5000);
    energy_hadron_reco_true_1Track->GetXaxis()->SetTitle("E_{Reco}^{Total} - E_{Reco}^{#mu}");
    energy_hadron_reco_true_1Track->GetYaxis()->SetTitle("E_{True}^{Total} - E_{True}^{#mu}");

    energy_hadron_reco_true_2Track = new TH2D("energy_hadron_reco_true_2Track","Hadronic Energy Truve vs Reco -- 2Track",100,0,5000,100,0,5000);
    energy_hadron_reco_true_2Track->GetXaxis()->SetTitle("E_{Reco}^{Total} - E_{Reco}^{#mu}");
    energy_hadron_reco_true_2Track->GetYaxis()->SetTitle("E_{True}^{Total} - E_{True}^{#mu}");

    vertex_evis_true_proton_KE = new TH2D("vertex_evis_true_proton_KE","True Proton KE vs Vertex Blob Evis (Signal, 1Track)",50,0.0,500.0,50,0.0,500.0);
    vertex_evis_true_proton_KE->GetXaxis()->SetTitle("E_{Visible} [MeV]");
    vertex_evis_true_proton_KE->GetYaxis()->SetTitle("T_{Proton}^{True} [MeV]");

    vertex_evis_vertex_evis_ratio = new TH2D("vertex_evis_vertex_evis_ratio","True Proton KE/Vertex Blob Evis vs Vertex Blob Evis (Signal, 1Track)",10,0.0,500.0,10,0.0,2.0);
    vertex_evis_vertex_evis_ratio->GetXaxis()->SetTitle("E_{Visible} [MeV]");
    vertex_evis_vertex_evis_ratio->GetYaxis()->SetTitle("T_{Proton}^{True}/E_{Visible} [MeV]");

    // Vertex Sphere Energy
    extra_energy_50_1Track = new TH1D("extra_energy_50_1Track","Extra Energy radius = 50 mm -- 1Track",50,0.0,500.0);
    extra_energy_50_1Track->GetXaxis()->SetTitle("Extra Energy [MeV]");
    extra_energy_50_1Track->GetYaxis()->SetTitle("N(Events)");

    extra_energy_100_1Track = new TH1D("extra_energy_100_1Track","Extra Energy radius = 100 mm -- 1Track",50,0.0,500.0);
    extra_energy_100_1Track->GetXaxis()->SetTitle("Extra Energy [MeV]");
    extra_energy_100_1Track->GetYaxis()->SetTitle("N(Events)");

    extra_energy_150_1Track = new TH1D("extra_energy_150_1Track","Extra Energy radius = 150 mm -- 1Track",50,0.0,500.0);
    extra_energy_150_1Track->GetXaxis()->SetTitle("Extra Energy [MeV]");
    extra_energy_150_1Track->GetYaxis()->SetTitle("N(Events)");

    extra_energy_200_1Track = new TH1D("extra_energy_200_1Track","Extra Energy radius = 200 mm -- 1Track",50,0.0,500.0);
    extra_energy_200_1Track->GetXaxis()->SetTitle("Extra Energy [MeV]");
    extra_energy_200_1Track->GetYaxis()->SetTitle("N(Events)");

    extra_energy_300_1Track = new TH1D("extra_energy_300_1Track","Extra Energy radius = 300 mm -- 1Track",50,0.0,500.0);
    extra_energy_300_1Track->GetXaxis()->SetTitle("Extra Energy [MeV]");
    extra_energy_300_1Track->GetYaxis()->SetTitle("N(Events)");

    extra_energy_500_1Track = new TH1D("extra_energy_500_1Track","Extra Energy radius = 500 mm -- 1Track",50,0.0,500.0);
    extra_energy_500_1Track->GetXaxis()->SetTitle("Extra Energy [MeV]");
    extra_energy_500_1Track->GetYaxis()->SetTitle("N(Events)");

    extra_energy_50_2Track = new TH1D("extra_energy_50_2Track","Extra Energy radius = 50 mm -- 2Track",50,0.0,500.0);
    extra_energy_50_2Track->GetXaxis()->SetTitle("Extra Energy [MeV]");
    extra_energy_50_2Track->GetYaxis()->SetTitle("N(Events)");

    extra_energy_100_2Track = new TH1D("extra_energy_100_2Track","Extra Energy radius = 100 mm -- 2Track",50,0.0,500.0);
    extra_energy_100_2Track->GetXaxis()->SetTitle("Extra Energy [MeV]");
    extra_energy_100_2Track->GetYaxis()->SetTitle("N(Events)");

    extra_energy_150_2Track = new TH1D("extra_energy_150_2Track","Extra Energy radius = 150 mm -- 2Track",50,0.0,500.0);
    extra_energy_150_2Track->GetXaxis()->SetTitle("Extra Energy [MeV]");
    extra_energy_150_2Track->GetYaxis()->SetTitle("N(Events)");

    extra_energy_200_2Track = new TH1D("extra_energy_200_2Track","Extra Energy radius = 200 mm -- 2Track",50,0.0,500.0);
    extra_energy_200_2Track->GetXaxis()->SetTitle("Extra Energy [MeV]");
    extra_energy_200_2Track->GetYaxis()->SetTitle("N(Events)");

    extra_energy_300_2Track = new TH1D("extra_energy_300_2Track","Extra Energy radius = 300 mm -- 2Track",50,0.0,500.0);
    extra_energy_300_2Track->GetXaxis()->SetTitle("Extra Energy [MeV]");
    extra_energy_300_2Track->GetYaxis()->SetTitle("N(Events)");

    extra_energy_500_2Track = new TH1D("extra_energy_500_2Track","Extra Energy radius = 500 mm -- 2Track",50,0.0,500.0);
    extra_energy_500_2Track->GetXaxis()->SetTitle("Extra Energy [MeV]");
    extra_energy_500_2Track->GetYaxis()->SetTitle("N(Events)");

}

void CCProtonPi0_Interaction::writeHistograms()
{
    std::cout<<">> Writing "<<rootDir<<std::endl;
   
    f->cd();

    for (int i = 0; i < nHistograms; i++){
        // Event Kinematics
        Enu_True[i]->Write();
        Enu_1Track[i]->Write();
        Enu_1Track_Alt[i]->Write();
        Enu_2Track[i]->Write();
        QSq[i]->Write();
        WSq[i]->Write();
        W[i]->Write();
       
        // Vertex & Extra Energy
        vertex_energy_1Track[i]->Write();
        vertex_energy_2Track[i]->Write();
        vertex_energy_Corrected_1Track[i]->Write();
        vertex_energy_Corrected_2Track[i]->Write();
        vertex_evis_1Track[i]->Write();
        vertex_evis_2Track[i]->Write();
        extra_energy_1Track[i]->Write();
        extra_energy_2Track[i]->Write();

        // Reconstruction 
        E_Unused_afterReco[i]->Write();
        E_Used_afterReco[i]->Write();

        // Other Event Parameters 
        deltaInvMass[i]->Write();
        
        // --------------------------------------------------------------------
        // Visible Energy -- Neutrino Energy Study
        // --------------------------------------------------------------------
        // 1 Track
        evis_total_1Track[i]->Write();
        evis_muon_1Track[i]->Write();
        evis_pi0_1Track[i]->Write();

        evis_hadron_1Track[i]->Write();
        energy_hadron_true_1Track[i]->Write();
        energy_hadron_reco_1Track[i]->Write();

        evis_extra_1Track[i]->Write();
        energy_extra_true_1Track[i]->Write();

        // 2 Track
        evis_total_2Track[i]->Write();
        evis_muon_2Track[i]->Write();
        evis_pi0_2Track[i]->Write();
        evis_proton_2Track[i]->Write();

        evis_hadron_2Track[i]->Write();
        energy_hadron_true_2Track[i]->Write();
        energy_hadron_reco_2Track[i]->Write();

        evis_hadron_nopi0_2Track[i]->Write();
        energy_hadron_nopi0_true_2Track[i]->Write();
        
        evis_extra_2Track[i]->Write();
        energy_extra_true_2Track[i]->Write();
        // --------------------------------------------------------------------
    }
    
    // MC Only Histograms
    final_mc_w_DIS->Write();
    final_mc_w_RES->Write();
    final_mc_w_CCQE->Write();

    proton_true_P_1Track->Write();
    proton_true_KE_1Track->Write();

    Enu_1Track_Error->Write();
    Enu_1Track_Alt_Error->Write();
    Enu_2Track_Error->Write();
    Enu_1Track_Corrected_Error->Write();
    Enu_2Track_Corrected_Error->Write();

    Enu_1Track_True->Write();
    Enu_1Track_Alt_True->Write();
    Enu_1Track_1Track_Alt->Write();
    Enu_2Track_True->Write();
  
    energy_hadron_reco_true_1Track->Write();
    energy_hadron_reco_true_2Track->Write();
    
    evis_extra_energy_extra_true_1Track->Write();
    evis_extra_energy_extra_true_2Track->Write();

    vertex_evis_true_proton_KE->Write();
    vertex_evis_vertex_evis_ratio->Write();

    extra_energy_50_1Track->Write();
    extra_energy_100_1Track->Write();
    extra_energy_150_1Track->Write();
    extra_energy_200_1Track->Write();
    extra_energy_300_1Track->Write();
    extra_energy_500_1Track->Write();

    extra_energy_50_2Track->Write();
    extra_energy_100_2Track->Write();
    extra_energy_150_2Track->Write();
    extra_energy_200_2Track->Write();
    extra_energy_300_2Track->Write();
    extra_energy_500_2Track->Write();

    f->Close();
}



#endif
