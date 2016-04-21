/*
    See CCProtonPi0_Interaction.h header for Class Information
*/
#ifndef CCProtonPi0_Interaction_cpp
#define CCProtonPi0_Interaction_cpp

#include "CCProtonPi0_Interaction.h"

using namespace PlotUtils;

CCProtonPi0_Interaction::CCProtonPi0_Interaction(bool isModeReduce, bool isMC) {
    std::cout<<"Initializing CCProtonPi0_Interaction"<<std::endl;
    
    if(isModeReduce){
        std::cout<<"\tNTuple Reduce Mode -- Will not create ROOT Files"<<std::endl;
    }else{
        if (isMC) rootDir = Folder_List::rootDir_Interaction_mc;
        else rootDir = Folder_List::rootDir_Interaction_data;
        
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

        temp = new MnvH1D( Form("%s_%d","Enu",i),"Reconstructed Beam Energy - All Events",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed E_{#nu} - All Events [GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.beamE.get_width()));
        Enu.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","QSq",i),"Reconstructed Q^{2}",binList.Q2.get_nBins(), binList.Q2.get_min(), binList.Q2.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed Q^{2} [GeV^{2}]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.Q2.get_width()));
        QSq.push_back(temp);    

        temp = new MnvH1D( Form("%s_%d","WSq",i),"Reconstructed W^{2}",binList.wSq.get_nBins(), binList.wSq.get_min(), binList.wSq.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed W^{2} [GeV^{2}]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.wSq.get_width()));
        WSq.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","W",i),"Reconstructed W",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed W [GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.w.get_width()));
        W.push_back(temp);
 
        temp = new MnvH1D( Form("%s_%d","W_Calc",i),"Reconstructed W",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed W [GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.w.get_width()));
        W_Calc.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","vertex_energy_1Track",i),"Vertex Blob Energy (r = 90mm)",binList.vertex_energy.get_nBins(), binList.vertex_energy.get_min(), binList.vertex_energy.get_max() );
        temp->GetXaxis()->SetTitle("Vertex Blob Energy [MeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.vertex_energy.get_width()));
        vertex_energy_1Track.push_back(temp);
  
        temp = new MnvH1D( Form("%s_%d","vertex_energy_2Track",i),"Vertex Blob Energy (r = 90mm)",binList.vertex_energy.get_nBins(), binList.vertex_energy.get_min(), binList.vertex_energy.get_max() );
        temp->GetXaxis()->SetTitle("Vertex Blob Energy [MeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.vertex_energy.get_width()));
        vertex_energy_2Track.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","vertex_evis_1Track",i),"Vertex Blob Visible Energy (r = 90mm)",binList.vertex_evis.get_nBins(), binList.vertex_evis.get_min(), binList.vertex_evis.get_max() );
        temp->GetXaxis()->SetTitle("Vertex Blob Visible Energy [MeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.vertex_evis.get_width()));
        vertex_evis_1Track.push_back(temp);
  
        temp = new MnvH1D( Form("%s_%d","vertex_evis_2Track",i),"Vertex Blob Visible Energy (r = 90mm)",binList.vertex_evis.get_nBins(), binList.vertex_evis.get_min(), binList.vertex_evis.get_max() );
        temp->GetXaxis()->SetTitle("Vertex Blob Visible Energy [MeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.vertex_evis.get_width()));
        vertex_evis_2Track.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","extra_evis_1Track",i),"Extra Leftover Visible Energy (r = 300mm)",binList.extra_energy.get_nBins(), binList.extra_energy.get_min(), binList.extra_energy.get_max() );
        temp->GetXaxis()->SetTitle("Extra Leftover Energy [MeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.extra_energy.get_width()));
        extra_evis_1Track.push_back(temp);
 
        temp = new MnvH1D( Form("%s_%d","extra_evis_2Track",i),"Extra Leftover Visible Energy (r = 300mm)",binList.extra_energy.get_nBins(), binList.extra_energy.get_min(), binList.extra_energy.get_max() );
        temp->GetXaxis()->SetTitle("Extra Leftover Energy [MeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.extra_energy.get_width()));
        extra_evis_2Track.push_back(temp);
      
        temp = new MnvH1D( Form("%s_%d","deltaInvMass",i),"Reconstructed #Delta^{+} Invariant Mass",binList.deltaInvMass.get_nBins(), binList.deltaInvMass.get_min(), binList.deltaInvMass.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed #Delta^{+} Inv. Mass [GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f [MeV] ",binList.deltaInvMass.get_width()));
        deltaInvMass.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","vertex_z",i),"Reconstructed Interaction Vertex",binList.vertex_z.get_nBins(), binList.vertex_z.get_min(), binList.vertex_z.get_max() );
        temp->GetXaxis()->SetTitle("z = 4293 Target, #bf{z = 5810 Interaction Region}, z = 8614 ECAL, z = 9088 HCAL");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.vertex_z.get_width()));
        vertex_z.push_back(temp);

        // Extra Energy
        temp = new MnvH1D( Form("%s_%d","extra_dispersed_energy_1Track",i),"Extra Dispersed Energy 1 Track", 20, 0.0, 500.0 );
        temp->GetXaxis()->SetTitle("Extra Dispersed Energy [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        extra_dispersed_energy_1Track.push_back(temp);
        
        temp = new MnvH1D( Form("%s_%d","extra_muon_energy_1Track",i),"Extra Muon Energy 1 Track", 20, 0.0, 500.0 );
        temp->GetXaxis()->SetTitle("Extra Muon Energy [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        extra_muon_energy_1Track.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","extra_rejected_energy_1Track",i),"Extra Rejected Visible Energy 1 Track", 20, 0.0, 500.0 );
        temp->GetXaxis()->SetTitle("Extra Rejected Visible Energy [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        extra_rejected_energy_1Track.push_back(temp);
    
        temp = new MnvH1D( Form("%s_%d","extra_total_energy_1Track",i),"Total Extra Energy 1 Track", 20, 0.0, 500.0 );
        temp->GetXaxis()->SetTitle("Total Extra Energy [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        extra_total_energy_1Track.push_back(temp);
        
        temp = new MnvH1D( Form("%s_%d","extra_dispersed_energy_2Track",i),"Extra Dispersed Energy 2 Track", 20, 0.0, 500.0 );
        temp->GetXaxis()->SetTitle("Extra Dispersed Energy [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        extra_dispersed_energy_2Track.push_back(temp);
        
        temp = new MnvH1D( Form("%s_%d","extra_muon_energy_2Track",i),"Extra Muon Energy 2 Track", 20, 0.0, 500.0 );
        temp->GetXaxis()->SetTitle("Extra Muon Energy [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        extra_muon_energy_2Track.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","extra_rejected_energy_2Track",i),"Extra Rejected Visible Energy 2 Track", 20, 0.0, 500.0 );
        temp->GetXaxis()->SetTitle("Extra Rejected Visible Energy [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        extra_rejected_energy_2Track.push_back(temp);
    
        temp = new MnvH1D( Form("%s_%d","extra_total_energy_2Track",i),"Total Extra Energy 2 Track", 20, 0.0, 500.0 );
        temp->GetXaxis()->SetTitle("Total Extra Energy [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        extra_total_energy_2Track.push_back(temp);

    }

    // Cross Section Variables
    QSq_all = new MnvH1D( "QSq_all","Data All Q^{2}", binList.size_QSq, binList.a_QSq);
    QSq_all->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    QSq_all->GetYaxis()->SetTitle("N(Events)");

    QSq_mc_truth_signal = new MnvH1D( "QSq_mc_truth_signal","MC Truth Signal Q^{2}", binList.size_QSq, binList.a_QSq);
    QSq_mc_truth_signal->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    QSq_mc_truth_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_MC(QSq_mc_truth_signal);

    QSq_mc_reco_all = new MnvH1D( "QSq_mc_reco_all","MC All Reconstructed Q^{2}", binList.size_QSq, binList.a_QSq);
    QSq_mc_reco_all->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    QSq_mc_reco_all->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_MC(QSq_mc_reco_all);

    QSq_mc_reco_signal = new MnvH1D( "QSq_mc_reco_signal","MC Reconstructed Signal Q^{2}", binList.size_QSq, binList.a_QSq);
    QSq_mc_reco_signal->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    QSq_mc_reco_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_MC(QSq_mc_reco_signal);

    QSq_mc_reco_bckg = new MnvH1D( "QSq_mc_reco_bckg","MC Reconstructed Background Q^{2}", binList.size_QSq, binList.a_QSq);
    QSq_mc_reco_bckg->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    QSq_mc_reco_bckg->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_MC(QSq_mc_reco_bckg);

    QSq_response = new MnvH2D( "QSq_response","Signal Q^{2}", binList.size_QSq, binList.a_QSq, binList.size_QSq, binList.a_QSq);
    QSq_response->GetXaxis()->SetTitle("Reconstructed Q^{2} [GeV^{2}]");
    QSq_response->GetYaxis()->SetTitle("True Q^{2} [GeV^{2}]");
    AddVertErrorBands_MC(QSq_response);

    Enu_all = new MnvH1D( "Enu_all","Data All E_{#nu}", binList.size_Enu, binList.a_Enu);
    Enu_all->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    Enu_all->GetYaxis()->SetTitle("N(Events)");

    Enu_mc_truth_signal = new MnvH1D( "Enu_mc_truth_signal","MC Truth Signal E_{#nu}", binList.size_Enu, binList.a_Enu);
    Enu_mc_truth_signal->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    Enu_mc_truth_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_MC(Enu_mc_truth_signal);

    Enu_mc_reco_all = new MnvH1D( "Enu_mc_reco_all","MC All Reconstructed E_{#nu}", binList.size_Enu, binList.a_Enu);
    Enu_mc_reco_all->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    Enu_mc_reco_all->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_MC(Enu_mc_reco_all);

    Enu_mc_reco_signal = new MnvH1D( "Enu_mc_reco_signal","MC Reconstructed Signal E_{#nu}", binList.size_Enu, binList.a_Enu);
    Enu_mc_reco_signal->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    Enu_mc_reco_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_MC(Enu_mc_reco_signal);

    Enu_mc_reco_bckg = new MnvH1D( "Enu_mc_reco_bckg","MC Reconstructed Background E_{#nu}", binList.size_Enu, binList.a_Enu);
    Enu_mc_reco_bckg->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    Enu_mc_reco_bckg->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_MC(Enu_mc_reco_bckg);

    Enu_response = new MnvH2D( "Enu_response","Signal E_{#nu}", binList.size_Enu, binList.a_Enu, binList.size_Enu, binList.a_Enu);
    Enu_response->GetXaxis()->SetTitle("Reconstructed E_{#nu} [GeV]");
    Enu_response->GetYaxis()->SetTitle("True E_{#nu} [GeV]");
    AddVertErrorBands_MC(Enu_response);

    
    // MC Only Histograms
    final_mc_w_DIS = new TH1D( "final_mc_w_DIS","True W for DIS",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max() );
    final_mc_w_DIS->GetXaxis()->SetTitle("True W for DIS [GeV]");
    final_mc_w_DIS->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.mc_w.get_width()));
    
    final_mc_w_RES = new TH1D( "final_mc_w_RES","True W for RES",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max() );
    final_mc_w_RES->GetXaxis()->SetTitle("True W for RES [GeV]");
    final_mc_w_RES->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.mc_w.get_width()));
   
    final_mc_Q2_DIS = new TH1D( "final_mc_Q2_DIS","True Q^{2} for DIS",binList.mc_Q2.get_nBins(), binList.mc_Q2.get_min(), binList.mc_Q2.get_max() );
    final_mc_Q2_DIS->GetXaxis()->SetTitle("True Q^{2} for DIS [GeV^{2}]");
    final_mc_Q2_DIS->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.mc_Q2.get_width()));
    
    final_mc_Q2_RES = new TH1D( "final_mc_Q2_RES","True Q^{2} for RES",binList.mc_Q2.get_nBins(), binList.mc_Q2.get_min(), binList.mc_Q2.get_max() );
    final_mc_Q2_RES->GetXaxis()->SetTitle("True Q^{2} for RES [GeV^{2}]");
    final_mc_Q2_RES->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.mc_Q2.get_width()));

    // QSq: True, Error, Diff
    QSq_True = new MnvH1D( "QSq_True","True Q^{2}",binList.QSq.get_nBins(), binList.QSq.get_min(), binList.QSq.get_max() );
    QSq_True->GetXaxis()->SetTitle("Q^{2} [Gev^{2}]");
    QSq_True->GetYaxis()->SetTitle("N(Events)");

    QSq_Error = new MnvH1D( "QSq_Error","Q^{2} Error",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    QSq_Error->GetXaxis()->SetTitle("(Q^{2}_{Reco} - Q^{2}_{True})/Q^{2}_{True}");
    QSq_Error->GetYaxis()->SetTitle("N(Events)");

    QSq_Diff = new MnvH1D( "QSq_Diff","Q^{2} Difference ",binList.beamE_Diff.get_nBins(), binList.beamE_Diff.get_min(), binList.beamE_Diff.get_max() );
    QSq_Diff->GetXaxis()->SetTitle("Q^{2}_{Reco} - Q^{2}_{True}");
    QSq_Diff->GetYaxis()->SetTitle("N(Events)");

    // Short Proton True Information
    proton_true_P_1Track = new TH1D( "proton_true_P_1Track","Short Proton True Momentum",binList.short_proton_P.get_nBins(), binList.short_proton_P.get_min(), binList.short_proton_P.get_max() );
    proton_true_P_1Track->GetXaxis()->SetTitle("Short Proton P_{True} [MeV]");
    proton_true_P_1Track->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.short_proton_P.get_width()));

    proton_true_KE_1Track = new TH1D( "proton_true_KE_1Track","Short Proton True Momentum",binList.short_proton_KE.get_nBins(), binList.short_proton_KE.get_min(), binList.short_proton_KE.get_max() );
    proton_true_KE_1Track->GetXaxis()->SetTitle("Short Proton P_{True} [MeV]");
    proton_true_KE_1Track->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.short_proton_KE.get_width()));

    // Ejected Nucleons
    n_ejected_nucleons_1Track = new TH1D( "n_ejected_nucleons_1Track","N(Nucleons) out of Nucleus",binList.multiplicity.get_nBins(), binList.multiplicity.get_min(), binList.multiplicity.get_max() );
    n_ejected_nucleons_1Track->GetXaxis()->SetTitle("N(Nucleons) out of Nucleus");
    n_ejected_nucleons_1Track->GetYaxis()->SetTitle("N(Events)");

    n_ejected_nucleons_2Track = new TH1D( "n_ejected_nucleons_2Track","N(Nucleons) out of Nucleus",binList.multiplicity.get_nBins(), binList.multiplicity.get_min(), binList.multiplicity.get_max() );
    n_ejected_nucleons_2Track->GetXaxis()->SetTitle("N(Nucleons) out of Nucleus");
    n_ejected_nucleons_2Track->GetYaxis()->SetTitle("N(Events)");

    // ------------------------------------------------------------------------
    // Neutrino Energy: Truth, Error, Difference
    // ------------------------------------------------------------------------
    Enu_True_1Track = new TH1D("Enu_True_1Track","True Neutrino Energy - 1 Track",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    Enu_True_1Track->GetXaxis()->SetTitle("(E_{#nu}^{Reco}-E_{#nu}^{True})/E_{#nu}^{True}");
    Enu_True_1Track->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.beamE.get_width()));
  
    Enu_True_2Track = new TH1D("Enu_True_2Track","True Neutrino Energy - 1 Track",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    Enu_True_2Track->GetXaxis()->SetTitle("(E_{#nu}^{Reco}-E_{#nu}^{True})/E_{#nu}^{True}");
    Enu_True_2Track->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.beamE.get_width()));
    
    Enu_1Track_Error = new TH1D("Enu_1Track_Error","Neutrino Energy Error - 1 Track",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    Enu_1Track_Error->GetXaxis()->SetTitle("(E_{#nu}^{Reco}-E_{#nu}^{True})/E_{#nu}^{True}");
    Enu_1Track_Error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));

    Enu_1Track_Alt_Error = new TH1D("Enu_1Track_Alt_Error","Neutrino Energy Error - 1 Track Alternative",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    Enu_1Track_Alt_Error->GetXaxis()->SetTitle("(E_{#nu}^{Reco}-E_{#nu}^{True})/E_{#nu}^{True}");
    Enu_1Track_Alt_Error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));

    Enu_2Track_Error = new TH1D("Enu_2Track_Error","Neutrino Energy Error - 2 Track",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    Enu_2Track_Error->GetXaxis()->SetTitle("(E_{#nu}^{Reco}-E_{#nu}^{True})/E_{#nu}^{True}");
    Enu_2Track_Error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));

    Enu_1Track_Diff = new TH1D("Enu_1Track_Diff","Neutrino Energy Difference - 1 Track",binList.beamE_Diff.get_nBins(), binList.beamE_Diff.get_min(), binList.beamE_Diff.get_max() );
    Enu_1Track_Diff->GetXaxis()->SetTitle("E_{#nu}^{Reco}-E_{#nu}^{True} [GeV]");
    Enu_1Track_Diff->GetYaxis()->SetTitle("N(Events)");

    Enu_2Track_Diff = new TH1D("Enu_2Track_Diff","Neutrino Energy Difference - 2 Track",binList.beamE_Diff.get_nBins(), binList.beamE_Diff.get_min(), binList.beamE_Diff.get_max() );
    Enu_2Track_Diff->GetXaxis()->SetTitle("E_{#nu}^{Reco}-E_{#nu}^{True} [GeV]");
    Enu_2Track_Diff->GetYaxis()->SetTitle("N(Events)");

    // Extra Energy
    h_extra_muon_energy = new TH1D("h_extra_muon_energy","Extra Muon Energy",20,0.0,1000);
    h_extra_muon_energy->GetXaxis()->SetTitle("Extra Muon Energy [MeV]");
    h_extra_muon_energy->GetYaxis()->SetTitle("N(Events)");

    h_extra_dispersed_energy = new TH1D("h_extra_dispersed_energy","Extra Dispersed Energy",20,0.0,1000);
    h_extra_dispersed_energy->GetXaxis()->SetTitle("Extra Dispersed Energy [MeV]");
    h_extra_dispersed_energy->GetYaxis()->SetTitle("N(Events)");
 
    h_extra_rejected_energy = new TH1D("h_extra_rejected_energy","Extra Rejected Energy",20,0.0,1000);
    h_extra_rejected_energy->GetXaxis()->SetTitle("Extra Rejected Energy [MeV]");
    h_extra_rejected_energy->GetYaxis()->SetTitle("N(Events)");
}

void CCProtonPi0_Interaction::writeHistograms()
{
    std::cout<<">> Writing "<<rootDir<<std::endl;
   
    f->cd();

    for (int i = 0; i < nHistograms; i++){
        // Event Kinematics
        Enu_1Track[i]->Write();
        Enu_1Track_Alt[i]->Write();
        Enu_2Track[i]->Write();
        Enu[i]->Write();
        QSq[i]->Write();
        WSq[i]->Write();
        W[i]->Write();
        W_Calc[i]->Write();
       
        // Vertex & Extra Energy
        vertex_energy_1Track[i]->Write();
        vertex_energy_2Track[i]->Write();
        vertex_evis_1Track[i]->Write();
        vertex_evis_2Track[i]->Write();
        extra_evis_1Track[i]->Write();
        extra_evis_2Track[i]->Write();

        // Other Event Parameters 
        deltaInvMass[i]->Write();

        // Extra Energy
        extra_muon_energy_1Track[i]->Write();
        extra_dispersed_energy_1Track[i]->Write();
        extra_rejected_energy_1Track[i]->Write();
        extra_total_energy_1Track[i]->Write();
        extra_muon_energy_2Track[i]->Write();
        extra_dispersed_energy_2Track[i]->Write();
        extra_rejected_energy_2Track[i]->Write();
        extra_total_energy_2Track[i]->Write();
    }
    
    QSq_all->Write();
    QSq_mc_truth_signal->Write();
    QSq_mc_reco_all->Write();
    QSq_mc_reco_signal->Write();
    QSq_mc_reco_bckg->Write();
    QSq_response->Write();
 
    Enu_all->Write();
    Enu_mc_truth_signal->Write();
    Enu_mc_reco_all->Write();
    Enu_mc_reco_signal->Write();
    Enu_mc_reco_bckg->Write();
    Enu_response->Write();
  
    // MC Only Histograms
    final_mc_w_DIS->Write();
    final_mc_w_RES->Write();
 
    final_mc_Q2_DIS->Write();
    final_mc_Q2_RES->Write();
    
    // Short Proton
    proton_true_P_1Track->Write();
    proton_true_KE_1Track->Write();
    
    // Ejected Nucleons
    n_ejected_nucleons_1Track->Write();
    n_ejected_nucleons_2Track->Write();

    //Extra Energy
    h_extra_muon_energy->Write();
    h_extra_dispersed_energy->Write();
    h_extra_rejected_energy->Write();
    
    // QSq Truth, Error, Difference
    QSq_True->Write();
    QSq_Error->Write();
    QSq_Diff->Write();
    
    // Neutrino Energy: Truth, Error, Difference
    Enu_True_1Track->Write();
    Enu_True_2Track->Write();
    
    Enu_1Track_Error->Write();
    Enu_1Track_Alt_Error->Write();
    Enu_2Track_Error->Write();

    Enu_1Track_Diff->Write();
    Enu_2Track_Diff->Write();

    f->Close();
}



#endif
