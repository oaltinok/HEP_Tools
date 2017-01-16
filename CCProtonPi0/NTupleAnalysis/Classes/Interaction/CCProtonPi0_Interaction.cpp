/*
    See CCProtonPi0_Interaction.h header for Class Information
*/
#ifndef CCProtonPi0_Interaction_cpp
#define CCProtonPi0_Interaction_cpp

#include "CCProtonPi0_Interaction.h"

using namespace PlotUtils;

CCProtonPi0_Interaction::CCProtonPi0_Interaction(bool isModeReduce, bool isMC) 
{
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
        temp = new MnvH1D( Form("%s_%d","CV_weight_2p2h",i),"Central Value Weight -- 2p2h", 80,0.0,2.0);
        temp->GetXaxis()->SetTitle("Central Value Weight -- 2p2h");
        temp->GetYaxis()->SetTitle("N(Events)");
        CV_weight_2p2h.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","CV_weight_Delta",i),"Central Value Weight -- Delta", 80,0.0,2.0);
        temp->GetXaxis()->SetTitle("Central Value Weight -- Delta");
        temp->GetYaxis()->SetTitle("N(Events)");
        CV_weight_Delta.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","CV_weight_CCRES",i),"Central Value Weight -- CCRES", 80,0.0,2.0);
        temp->GetXaxis()->SetTitle("Central Value Weight -- CCRES");
        temp->GetYaxis()->SetTitle("N(Events)");
        CV_weight_CCRES.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","CV_weight_NonRes1pi",i),"Central Value Weight -- NonRes1pi", 80,0.0,2.0);
        temp->GetXaxis()->SetTitle("Central Value Weight -- NonRes1pi");
        temp->GetYaxis()->SetTitle("N(Events)");
        CV_weight_NonRes1pi.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","CV_weight",i),"Central Value Weight", 80,0.0,2.0);
        temp->GetXaxis()->SetTitle("Central Value Weight");
        temp->GetYaxis()->SetTitle("N(Events)");
        CV_weight.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","err_2p2h",i),"2p2h 1#sigma Weights", 80,0.0,2.0);
        temp->GetXaxis()->SetTitle("2p2h 1#sigma Weights");
        temp->GetYaxis()->SetTitle("N(Events)");
        err_2p2h.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","genie_wgt_Theta_Delta2Npi",i),"Theta_Delta2Npi GENIE 1#sigma Weights", 80,0.0,2.0);
        temp->GetXaxis()->SetTitle("Theta_Delta2Npi GENIE 1#sigma Weights");
        temp->GetYaxis()->SetTitle("N(Events)");
        genie_wgt_Theta_Delta2Npi.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","updated_wgt_Theta_Delta2Npi",i),"Theta_Delta2Npi Updated 1#sigma Weights", 80,0.0,2.0);
        temp->GetXaxis()->SetTitle("Theta_Delta2Npi Updated 1#sigma Weights");
        temp->GetYaxis()->SetTitle("N(Events)");
        updated_wgt_Theta_Delta2Npi.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","genie_wgt_MaRES",i),"MaRES GENIE 1#sigma Weights", 80,0.0,2.0);
        temp->GetXaxis()->SetTitle("MaRES GENIE 1#sigma Weights");
        temp->GetYaxis()->SetTitle("N(Events)");
        genie_wgt_MaRES.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","updated_wgt_MaRES",i),"MaRES Updated 1#sigma Weights", 80,0.0,2.0);
        temp->GetXaxis()->SetTitle("MaRES Updated 1#sigma Weights");
        temp->GetYaxis()->SetTitle("N(Events)");
        updated_wgt_MaRES.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","genie_wgt_MvRES",i),"MvRES GENIE 1#sigma Weights", 80,0.0,2.0);
        temp->GetXaxis()->SetTitle("MvRES GENIE 1#sigma Weights");
        temp->GetYaxis()->SetTitle("N(Events)");
        genie_wgt_MvRES.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","updated_wgt_MvRES",i),"MvRES Updated 1#sigma Weights", 80,0.0,2.0);
        temp->GetXaxis()->SetTitle("MvRES Updated 1#sigma Weights");
        temp->GetYaxis()->SetTitle("N(Events)");
        updated_wgt_MvRES.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","genie_wgt_Rvn1pi",i),"Rvn1pi GENIE 1#sigma Weights", 80,0.0,2.0);
        temp->GetXaxis()->SetTitle("Rvn1pi GENIE 1#sigma Weights");
        temp->GetYaxis()->SetTitle("N(Events)");
        genie_wgt_Rvn1pi.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","updated_wgt_Rvn1pi",i),"Rvn1pi Updated 1#sigma Weights", 80,0.0,2.0);
        temp->GetXaxis()->SetTitle("Rvn1pi Updated 1#sigma Weights");
        temp->GetYaxis()->SetTitle("N(Events)");
        updated_wgt_Rvn1pi.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","Enu",i),"Reconstructed Beam Energy - All Events", binList.size_Enu, binList.a_Enu);
        temp->GetXaxis()->SetTitle("Reconstructed E_{#nu} - All Events [GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.beamE.get_width()));
        Enu.push_back(temp);
 
        temp = new MnvH1D( Form("%s_%d","Enu_1Track",i),"Reconstructed Beam Energy - 1 Track", binList.size_Enu, binList.a_Enu);
        temp->GetXaxis()->SetTitle("Reconstructed E_{#nu} - 1Track [GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.beamE.get_width()));
        Enu_1Track.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","Enu_2Track",i),"Reconstructed Beam Energy - 2 Track", binList.size_Enu, binList.a_Enu);
        temp->GetXaxis()->SetTitle("Reconstructed E_{#nu} - 2 Track [GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.beamE.get_width()));
        Enu_2Track.push_back(temp);
       
        temp = new MnvH1D( Form("%s_%d","QSq",i),"Reconstructed Q^{2}", binList.size_QSq, binList.a_QSq);
        temp->GetXaxis()->SetTitle("Reconstructed Q^{2} [GeV^{2}]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.Q2.get_width()));
        QSq.push_back(temp);    

        temp = new MnvH1D( Form("%s_%d","QSq_1Track",i),"Reconstructed Q^{2} 1Track", binList.size_QSq, binList.a_QSq);
        temp->GetXaxis()->SetTitle("Reconstructed Q^{2} [GeV^{2}]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.Q2.get_width()));
        QSq_1Track.push_back(temp);    

        temp = new MnvH1D( Form("%s_%d","QSq_2Track",i),"Reconstructed Q^{2} 2Track", binList.size_QSq, binList.a_QSq);
        temp->GetXaxis()->SetTitle("Reconstructed Q^{2} [GeV^{2}]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.Q2.get_width()));
        QSq_2Track.push_back(temp);    

        temp = new MnvH1D( Form("%s_%d","WSq",i),"Reconstructed W^{2}",binList.wSq.get_nBins(), binList.wSq.get_min(), binList.wSq.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed W^{2} [GeV^{2}]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.wSq.get_width()));
        WSq.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","WSq_1Track",i),"Reconstructed W^{2} 1Track",binList.wSq.get_nBins(), binList.wSq.get_min(), binList.wSq.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed W^{2} [GeV^{2}]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.wSq.get_width()));
        WSq_1Track.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","WSq_2Track",i),"Reconstructed W^{2} 2Track",binList.wSq.get_nBins(), binList.wSq.get_min(), binList.wSq.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed W^{2} [GeV^{2}]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.wSq.get_width()));
        WSq_2Track.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","W",i),"Reconstructed W",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed W [GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.w.get_width()));
        W.push_back(temp);
 
        temp = new MnvH1D( Form("%s_%d","W_1Track",i),"Reconstructed W 1Track",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed W [GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.w.get_width()));
        W_1Track.push_back(temp);
 
        temp = new MnvH1D( Form("%s_%d","W_2Track",i),"Reconstructed W 2Track",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed W [GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.w.get_width()));
        W_2Track.push_back(temp);

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

         
        temp = new MnvH1D( Form("%s_%d","deltaInvMass",i),"Reconstructed p#pi^{0} Invariant Mass",binList.deltaInvMass.get_nBins(), binList.deltaInvMass.get_min(), binList.deltaInvMass.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed p#pi^{0} Inv. Mass [GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f [MeV] ",binList.deltaInvMass.get_width()));
        deltaInvMass.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","vertex_z",i),"Reconstructed Interaction Vertex",binList.vertex_z.get_nBins(), binList.vertex_z.get_min(), binList.vertex_z.get_max() );
        temp->GetXaxis()->SetTitle("z = 4293 Target, #bf{z = 5810 Interaction Region}, z = 8614 ECAL, z = 9088 HCAL");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.vertex_z.get_width()));
        vertex_z.push_back(temp);

        // Extra Energy
         temp = new MnvH1D( Form("%s_%d","extra_leftover_energy_1Track",i),"Extra Leftover Visible Energy (r = 300mm)",binList.extra_energy.get_nBins(), binList.extra_energy.get_min(), binList.extra_energy.get_max() );
        temp->GetXaxis()->SetTitle("Extra Leftover Energy [MeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.extra_energy.get_width()));
        extra_leftover_energy_1Track.push_back(temp);
       
        temp = new MnvH1D( Form("%s_%d","extra_muon_energy_1Track",i),"Extra Muon Energy 1 Track",binList.extra_energy.get_nBins(), binList.extra_energy.get_min(), binList.extra_energy.get_max() );

        temp->GetXaxis()->SetTitle("Extra Muon Energy [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        extra_muon_energy_1Track.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","extra_rejected_energy_1Track",i),"Extra Rejected Visible Energy 1 Track",binList.extra_energy.get_nBins(), binList.extra_energy.get_min(), binList.extra_energy.get_max() );

        temp->GetXaxis()->SetTitle("Extra Rejected Visible Energy [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        extra_rejected_energy_1Track.push_back(temp);
    
        temp = new MnvH1D( Form("%s_%d","extra_total_energy_1Track",i),"Total Extra Energy 1 Track",binList.extra_energy.get_nBins(), binList.extra_energy.get_min(), binList.extra_energy.get_max() );

        temp->GetXaxis()->SetTitle("Total Extra Energy [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        extra_total_energy_1Track.push_back(temp);
     
        temp = new MnvH1D( Form("%s_%d","extra_leftover_energy_2Track",i),"Extra Leftover Visible Energy (r = 300mm)",binList.extra_energy.get_nBins(), binList.extra_energy.get_min(), binList.extra_energy.get_max() );
        temp->GetXaxis()->SetTitle("Extra Leftover Energy [MeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.extra_energy.get_width()));
        extra_leftover_energy_2Track.push_back(temp);
       
        temp = new MnvH1D( Form("%s_%d","extra_muon_energy_2Track",i),"Extra Muon Energy 2 Track",binList.extra_energy.get_nBins(), binList.extra_energy.get_min(), binList.extra_energy.get_max() );

        temp->GetXaxis()->SetTitle("Extra Muon Energy [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        extra_muon_energy_2Track.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","extra_rejected_energy_2Track",i),"Extra Rejected Visible Energy 2 Track",binList.extra_energy.get_nBins(), binList.extra_energy.get_min(), binList.extra_energy.get_max() );

        temp->GetXaxis()->SetTitle("Extra Rejected Visible Energy [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        extra_rejected_energy_2Track.push_back(temp);
    
        temp = new MnvH1D( Form("%s_%d","extra_total_energy_2Track",i),"Total Extra Energy 2 Track",binList.extra_energy.get_nBins(), binList.extra_energy.get_min(), binList.extra_energy.get_max() );

        temp->GetXaxis()->SetTitle("Total Extra Energy [MeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        extra_total_energy_2Track.push_back(temp);

        temp = new MnvH1D(Form("%s_%d","W_p_pi0",i),"Reconstructed p#pi^{0} Invariant Mass",30,0.5,2.0);
        temp->GetXaxis()->SetTitle("Reconstructed p#pi^{0} Inv. Mass [GeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        W_p_pi0.push_back(temp);
    
        temp = new MnvH1D(Form("%s_%d","W_All",i),"Reconstructed W All Events",30,0.5,2.0);
        temp->GetXaxis()->SetTitle("Reconstructed W[GeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        W_All.push_back(temp);
     
        temp = new MnvH1D(Form("%s_%d","W_1",i),"Reconstructed W 1 Track Events",30,0.5,2.0);
        temp->GetXaxis()->SetTitle("Reconstructed W[GeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        W_1.push_back(temp);
 
        temp = new MnvH1D(Form("%s_%d","W_2",i),"Reconstructed W 2 Track Events",30,0.5,2.0);
        temp->GetXaxis()->SetTitle("Reconstructed W[GeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        W_2.push_back(temp);

        temp = new MnvH1D(Form("%s_%d","QSq_All",i),"Q^{2} for All Events", binList.size_QSq, binList.a_QSq);
        temp->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
        temp->GetYaxis()->SetTitle("N(Events)");
        QSq_All.push_back(temp);
 
        temp = new MnvH1D(Form("%s_%d","QSq_LowEnu",i),"Q^{2} for E_{#nu} < 4 GeV", binList.size_QSq, binList.a_QSq);
        temp->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
        temp->GetYaxis()->SetTitle("N(Events)");
        QSq_LowEnu.push_back(temp);
 
        temp = new MnvH1D(Form("%s_%d","QSq_HighEnu",i),"Q^{2} for E_{#nu} > 4 GeV", binList.size_QSq, binList.a_QSq);
        temp->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
        temp->GetYaxis()->SetTitle("N(Events)");
        QSq_HighEnu.push_back(temp);
    }

    TH1D* h_temp;
    for (int i = 0; i <= 100; ++i){
        h_temp = new TH1D(Form("%s_%d","QSq_LowMaRES",i),"Q^{2}", 20, 0.0, 2.0);
        h_temp->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
        h_temp->GetYaxis()->SetTitle("N(Events)");
        QSq_LowMaRES.push_back(h_temp);

        h_temp = new TH1D(Form("%s_%d","QSq_HighMaRES",i),"Q^{2}", 20, 0.0, 2.0);
        h_temp->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
        h_temp->GetYaxis()->SetTitle("N(Events)");
        QSq_HighMaRES.push_back(h_temp);
    }

    // Cross Section Variables
    QSq_all = new MnvH1D( "QSq_all","Data All Q^{2}", binList.size_QSq, binList.a_QSq);
    QSq_all->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    QSq_all->GetYaxis()->SetTitle("N(Events)");

    QSq_mc_truth_signal = new MnvH1D( "QSq_mc_truth_signal","MC Truth Signal Q^{2}", binList.size_QSq, binList.a_QSq);
    QSq_mc_truth_signal->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    QSq_mc_truth_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_MC(QSq_mc_truth_signal);
    AddLatErrorBands_MC(QSq_mc_truth_signal);

    QSq_mc_reco_all = new MnvH1D( "QSq_mc_reco_all","MC All Reconstructed Q^{2}", binList.size_QSq, binList.a_QSq);
    QSq_mc_reco_all->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    QSq_mc_reco_all->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_MC(QSq_mc_reco_all);
    AddLatErrorBands_MC(QSq_mc_reco_all);

    QSq_mc_reco_signal = new MnvH1D( "QSq_mc_reco_signal","MC Reconstructed Signal Q^{2}", binList.size_QSq, binList.a_QSq);
    QSq_mc_reco_signal->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    QSq_mc_reco_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_MC(QSq_mc_reco_signal);
    AddLatErrorBands_MC(QSq_mc_reco_signal);

    QSq_mc_reco_bckg = new MnvH1D( "QSq_mc_reco_bckg","MC Reconstructed Background Q^{2}", binList.size_QSq, binList.a_QSq);
    QSq_mc_reco_bckg->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    QSq_mc_reco_bckg->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_MC(QSq_mc_reco_bckg);
    AddLatErrorBands_MC(QSq_mc_reco_bckg);

    QSq_response = new MnvH2D( "QSq_response","Signal Q^{2}", binList.size_QSq, binList.a_QSq, binList.size_QSq, binList.a_QSq);
    QSq_response->GetXaxis()->SetTitle("Reconstructed Q^{2} [GeV^{2}]");
    QSq_response->GetYaxis()->SetTitle("True Q^{2} [GeV^{2}]");
    AddVertErrorBands_MC(QSq_response);
    AddLatErrorBands_MC(QSq_response);

    Enu_all = new MnvH1D( "Enu_all","Data All E_{#nu}", binList.size_Enu, binList.a_Enu);
    Enu_all->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    Enu_all->GetYaxis()->SetTitle("N(Events)");

    Enu_mc_truth_signal = new MnvH1D( "Enu_mc_truth_signal","MC Truth Signal E_{#nu}", binList.size_Enu, binList.a_Enu);
    Enu_mc_truth_signal->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    Enu_mc_truth_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_MC(Enu_mc_truth_signal);
    AddLatErrorBands_MC(Enu_mc_truth_signal);

    Enu_mc_reco_all = new MnvH1D( "Enu_mc_reco_all","MC All Reconstructed E_{#nu}", binList.size_Enu, binList.a_Enu);
    Enu_mc_reco_all->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    Enu_mc_reco_all->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_MC(Enu_mc_reco_all);
    AddLatErrorBands_MC(Enu_mc_reco_all);

    Enu_mc_reco_signal = new MnvH1D( "Enu_mc_reco_signal","MC Reconstructed Signal E_{#nu}", binList.size_Enu, binList.a_Enu);
    Enu_mc_reco_signal->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    Enu_mc_reco_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_MC(Enu_mc_reco_signal);
    AddLatErrorBands_MC(Enu_mc_reco_signal);

    Enu_mc_reco_bckg = new MnvH1D( "Enu_mc_reco_bckg","MC Reconstructed Background E_{#nu}", binList.size_Enu, binList.a_Enu);
    Enu_mc_reco_bckg->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    Enu_mc_reco_bckg->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_MC(Enu_mc_reco_bckg);
    AddLatErrorBands_MC(Enu_mc_reco_bckg);

    Enu_response = new MnvH2D( "Enu_response","Signal E_{#nu}", binList.size_Enu, binList.a_Enu, binList.size_Enu, binList.a_Enu);
    Enu_response->GetXaxis()->SetTitle("Reconstructed E_{#nu} [GeV]");
    Enu_response->GetYaxis()->SetTitle("True E_{#nu} [GeV]");
    AddVertErrorBands_MC(Enu_response);
    AddLatErrorBands_MC(Enu_response);

    W_all = new MnvH1D( "W_all","Data All W",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max() );
    W_all->GetXaxis()->SetTitle("W [GeV]");
    W_all->GetYaxis()->SetTitle("N(Events)");

    W_mc_truth_signal = new MnvH1D( "W_mc_truth_signal","MC Truth Signal W",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max() );
    W_mc_truth_signal->GetXaxis()->SetTitle("W [GeV]");
    W_mc_truth_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_MC(W_mc_truth_signal);
    AddLatErrorBands_MC(W_mc_truth_signal);

    W_mc_reco_all = new MnvH1D( "W_mc_reco_all","MC All Reconstructed W",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max() );
    W_mc_reco_all->GetXaxis()->SetTitle("W [GeV]");
    W_mc_reco_all->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_MC(W_mc_reco_all);
    AddLatErrorBands_MC(W_mc_reco_all);

    W_mc_reco_signal = new MnvH1D( "W_mc_reco_signal","MC Reconstructed Signal W",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max() );
    W_mc_reco_signal->GetXaxis()->SetTitle("W [GeV]");
    W_mc_reco_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_MC(W_mc_reco_signal);
    AddLatErrorBands_MC(W_mc_reco_signal);

    W_mc_reco_bckg = new MnvH1D( "W_mc_reco_bckg","MC Reconstructed Background W",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max() );
    W_mc_reco_bckg->GetXaxis()->SetTitle("W [GeV]");
    W_mc_reco_bckg->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_MC(W_mc_reco_bckg);
    AddLatErrorBands_MC(W_mc_reco_bckg);

    W_response = new MnvH2D( "W_response","Signal W",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max(),binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max() );
    W_response->GetXaxis()->SetTitle("Reconstructed W [GeV]");
    W_response->GetYaxis()->SetTitle("True W [GeV]");
    AddVertErrorBands_MC(W_response);
    AddLatErrorBands_MC(W_response);

    // QSq Error, Diff
    WSq_QSq_Diff = new TH2D( "WSq_QSq_Diff","Q^{2} Difference vs W^{2}",binList.wSq.get_nBins(), binList.wSq.get_min(), binList.wSq.get_max(),binList.QSq_Diff.get_nBins(), binList.QSq_Diff.get_min(), binList.QSq_Diff.get_max() );
    WSq_QSq_Diff->GetXaxis()->SetTitle("W^{2} [GeV^{2}]");
    WSq_QSq_Diff->GetYaxis()->SetTitle("Q^{2}_{Reco} - Q^{2}_{True} [GeV^{2}]");

    QSq_All_response = new MnvH2D( "QSq_All_response","Signal Q^{2} All", binList.size_QSq, binList.a_QSq, binList.size_QSq, binList.a_QSq);
    QSq_All_response->GetXaxis()->SetTitle("Reconstructed Q^{2} [GeV^{2}]");
    QSq_All_response->GetYaxis()->SetTitle("True Q^{2} [GeV^{2}]");

    QSq_1Track_response = new MnvH2D( "QSq_1Track_response","Signal Q^{2} 1Track", binList.size_QSq, binList.a_QSq, binList.size_QSq, binList.a_QSq);
    QSq_1Track_response->GetXaxis()->SetTitle("Reconstructed Q^{2} [GeV^{2}]");
    QSq_1Track_response->GetYaxis()->SetTitle("True Q^{2} [GeV^{2}]");

    QSq_2Track_response = new MnvH2D( "QSq_2Track_response","Signal Q^{2} 1Track", binList.size_QSq, binList.a_QSq, binList.size_QSq, binList.a_QSq);
    QSq_2Track_response->GetXaxis()->SetTitle("Reconstructed Q^{2} [GeV^{2}]");
    QSq_2Track_response->GetYaxis()->SetTitle("True Q^{2} [GeV^{2}]");

    QSq_Error = new MnvH1D( "QSq_Error","Q^{2} Error",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    QSq_Error->GetXaxis()->SetTitle("(Q^{2}_{Reco} - Q^{2}_{True})/Q^{2}_{True}");
    QSq_Error->GetYaxis()->SetTitle("N(Events)");

    QSq_1Track_Error = new MnvH1D( "QSq_1Track_Error","Q^{2} Error - 1Track",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    QSq_1Track_Error->GetXaxis()->SetTitle("(Q^{2}_{Reco} - Q^{2}_{True})/Q^{2}_{True}");
    QSq_1Track_Error->GetYaxis()->SetTitle("N(Events)");

    QSq_2Track_Error = new MnvH1D( "QSq_2Track_Error","Q^{2} Error - 2Track",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    QSq_2Track_Error->GetXaxis()->SetTitle("(Q^{2}_{Reco} - Q^{2}_{True})/Q^{2}_{True}");
    QSq_2Track_Error->GetYaxis()->SetTitle("N(Events)");

    QSq_Diff = new MnvH1D( "QSq_Diff","Q^{2} Difference ",binList.QSq_Diff.get_nBins(), binList.QSq_Diff.get_min(), binList.QSq_Diff.get_max() );
    QSq_Diff->GetXaxis()->SetTitle("Q^{2}_{Reco} - Q^{2}_{True}");
    QSq_Diff->GetYaxis()->SetTitle("N(Events)");

    QSq_1Track_Diff = new MnvH1D( "QSq_1Track_Diff","Q^{2} Difference - 1Track",binList.QSq_Diff.get_nBins(), binList.QSq_Diff.get_min(), binList.QSq_Diff.get_max() );
    QSq_1Track_Diff->GetXaxis()->SetTitle("Q^{2}_{Reco} - Q^{2}_{True}");
    QSq_1Track_Diff->GetYaxis()->SetTitle("N(Events)");

    QSq_2Track_Diff = new MnvH1D( "QSq_2Track_Diff","Q^{2} Difference - 2Track",binList.QSq_Diff.get_nBins(), binList.QSq_Diff.get_min(), binList.QSq_Diff.get_max() );
    QSq_2Track_Diff->GetXaxis()->SetTitle("Q^{2}_{Reco} - Q^{2}_{True}");
    QSq_2Track_Diff->GetYaxis()->SetTitle("N(Events)");

    // Short Proton True Information
    nProtons = new MnvH1D( "nProtons","Number of Tracked Protons",5,0.0,5.0);
    nProtons->GetXaxis()->SetTitle("N(Protons)");
    nProtons->GetYaxis()->SetTitle("N(Events)");

    proton_true_P_1Track = new TH1D( "proton_true_P_1Track","Short Proton True Momentum",binList.short_proton_P.get_nBins(), binList.short_proton_P.get_min(), binList.short_proton_P.get_max() );
    proton_true_P_1Track->GetXaxis()->SetTitle("Short Proton P_{True} [MeV]");
    proton_true_P_1Track->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.short_proton_P.get_width()));

    proton_true_KE_1Track = new TH1D( "proton_true_KE_1Track","Short Proton True Momentum",binList.short_proton_KE.get_nBins(), binList.short_proton_KE.get_min(), binList.short_proton_KE.get_max() );
    proton_true_KE_1Track->GetXaxis()->SetTitle("Short Proton P_{True} [MeV]");
    proton_true_KE_1Track->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.short_proton_KE.get_width()));

    proton_true_theta_1Track = new TH1D( "proton_true_theta_1Track","Short Proton True #theta_{p}",binList.angle.get_nBins(), binList.angle.get_min(), binList.angle.get_max() );
    proton_true_theta_1Track->GetXaxis()->SetTitle("True #theta_{p} [Degree]");
    proton_true_theta_1Track->GetYaxis()->SetTitle(Form("Protons / %3.1f [Degree]",binList.angle.get_width()));

    Polarization_data = new MnvH1D( "Polarization_data","Pion Polarization in #Delta Rest Frame",10,-1,1);
    Polarization_data->GetXaxis()->SetTitle("Polarization");
    Polarization_data->GetYaxis()->SetTitle("N(Events)");

    Polarization_mc = new MnvH1D( "Polarization_mc","Pion Polarization in #Delta Rest Frame",10,-1,1);
    Polarization_mc->GetXaxis()->SetTitle("Polarization");
    Polarization_mc->GetYaxis()->SetTitle("N(Events)");

    DeltaTransverse_data = new MnvH1D( "DeltaTransverse_data","Delta Transverse Momentum (Double Transverse)",21,-500,500);
    DeltaTransverse_data->GetXaxis()->SetTitle("Delta Transverse Momentum [MeV]");
    DeltaTransverse_data->GetYaxis()->SetTitle("N(Events)");

    DeltaTransverse_mc = new MnvH1D( "DeltaTransverse_mc","Delta Transverse Momentum (Double Transverse)",21,-500,500);
    DeltaTransverse_mc->GetXaxis()->SetTitle("Delta Transverse Momentum [MeV]");
    DeltaTransverse_mc->GetYaxis()->SetTitle("N(Events)");

    DeltaTransverse_mc_res = new MnvH2D( "DeltaTransverse_mc_res","Double Transverse Momentum Residual vs Truth",21,-500,500,21,-50,50);
    DeltaTransverse_mc_res->GetXaxis()->SetTitle("Truth Double Transverse Momentum [MeV]");
    DeltaTransverse_mc_res->GetYaxis()->SetTitle("Double Transverse Momentum Residual [MeV]");

    // Ejected Nucleons
    n_ejected_nucleons_1Track = new TH1D( "n_ejected_nucleons_1Track","N(Nucleons) out of Nucleus",binList.multiplicity.get_nBins(), binList.multiplicity.get_min(), binList.multiplicity.get_max() );
    n_ejected_nucleons_1Track->GetXaxis()->SetTitle("N(Nucleons) out of Nucleus");
    n_ejected_nucleons_1Track->GetYaxis()->SetTitle("N(Events)");

    n_ejected_nucleons_2Track = new TH1D( "n_ejected_nucleons_2Track","N(Nucleons) out of Nucleus",binList.multiplicity.get_nBins(), binList.multiplicity.get_min(), binList.multiplicity.get_max() );
    n_ejected_nucleons_2Track->GetXaxis()->SetTitle("N(Nucleons) out of Nucleus");
    n_ejected_nucleons_2Track->GetYaxis()->SetTitle("N(Events)");

    // ------------------------------------------------------------------------
    // W: Truth, Error, Difference
    // ------------------------------------------------------------------------
    W_Error = new TH1D("W_Error","W Error",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    W_Error->GetXaxis()->SetTitle("(W_{Reco}-W_{True})/W_{True}");
    W_Error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));
 
    W_Diff = new TH1D("W_Diff","W Difference",binList.W_Diff.get_nBins(), binList.W_Diff.get_min(), binList.W_Diff.get_max() );
    W_Diff->GetXaxis()->SetTitle("W_{Reco}-W_{True} [GeV]");
    W_Diff->GetYaxis()->SetTitle("N(Events)");

    // ------------------------------------------------------------------------
    // Neutrino Energy: Truth, Error, Difference
    // ------------------------------------------------------------------------
    Enu_All_response = new MnvH2D( "Enu_All_response","Signal E_{#nu} All", binList.size_Enu, binList.a_Enu, binList.size_Enu, binList.a_Enu);
    Enu_All_response->GetXaxis()->SetTitle("Reconstructed E_{#nu} [GeV]");
    Enu_All_response->GetYaxis()->SetTitle("True E_{#nu} [GeV]");

    Enu_1Track_response = new MnvH2D( "Enu_1Track_response","Signal E_{#nu} 1Track", binList.size_Enu, binList.a_Enu, binList.size_Enu, binList.a_Enu);
    Enu_1Track_response->GetXaxis()->SetTitle("Reconstructed E_{#nu} [GeV]");
    Enu_1Track_response->GetYaxis()->SetTitle("True E_{#nu} [GeV]");

    Enu_2Track_response = new MnvH2D( "Enu_2Track_response","Signal E_{#nu} 2Track", binList.size_Enu, binList.a_Enu, binList.size_Enu, binList.a_Enu);
    Enu_2Track_response->GetXaxis()->SetTitle("Reconstructed E_{#nu} [GeV]");
    Enu_2Track_response->GetYaxis()->SetTitle("True E_{#nu} [GeV]");

    Enu_Error = new TH1D("Enu_Error","Neutrino Energy Error",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    Enu_Error->GetXaxis()->SetTitle("(E_{#nu}^{Reco}-E_{#nu}^{True})/E_{#nu}^{True}");
    Enu_Error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));
  
    Enu_1Track_Error = new TH1D("Enu_1Track_Error","Neutrino Energy Error - 1 Track",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    Enu_1Track_Error->GetXaxis()->SetTitle("(E_{#nu}^{Reco}-E_{#nu}^{True})/E_{#nu}^{True}");
    Enu_1Track_Error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));

    Enu_2Track_Error = new TH1D("Enu_2Track_Error","Neutrino Energy Error - 2 Track",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    Enu_2Track_Error->GetXaxis()->SetTitle("(E_{#nu}^{Reco}-E_{#nu}^{True})/E_{#nu}^{True}");
    Enu_2Track_Error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));

    Enu_Diff = new TH1D("Enu_Diff","Neutrino Energy Difference",binList.beamE_Diff.get_nBins(), binList.beamE_Diff.get_min(), binList.beamE_Diff.get_max() );
    Enu_Diff->GetXaxis()->SetTitle("E_{#nu}^{Reco}-E_{#nu}^{True} [GeV]");
    Enu_Diff->GetYaxis()->SetTitle("N(Events)");

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

    h_extra_leftover_energy = new TH1D("h_extra_leftover_energy","Extra Leftover Energy",20,0.0,1000);
    h_extra_leftover_energy->GetXaxis()->SetTitle("Extra Leftover Energy [MeV]");
    h_extra_leftover_energy->GetYaxis()->SetTitle("N(Events)");
 
    h_extra_rejected_energy = new TH1D("h_extra_rejected_energy","Extra Rejected Energy",20,0.0,1000);
    h_extra_rejected_energy->GetXaxis()->SetTitle("Extra Rejected Energy [MeV]");
    h_extra_rejected_energy->GetYaxis()->SetTitle("N(Events)");

    // ------------------------------------------------------------------------
    // Signal Q2
    // ------------------------------------------------------------------------
    mc_Q2_QE = new TH1D("mc_Q2_QE","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    mc_Q2_QE->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    mc_Q2_QE->GetYaxis()->SetTitle("N(Events)");

    mc_Q2_RES_1232 = new TH1D("mc_Q2_RES_1232","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    mc_Q2_RES_1232->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    mc_Q2_RES_1232->GetYaxis()->SetTitle("N(Events)");

    mc_Q2_RES_1535 = new TH1D("mc_Q2_RES_1535","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    mc_Q2_RES_1535->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    mc_Q2_RES_1535->GetYaxis()->SetTitle("N(Events)");

    mc_Q2_RES_1520 = new TH1D("mc_Q2_RES_1520","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    mc_Q2_RES_1520->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    mc_Q2_RES_1520->GetYaxis()->SetTitle("N(Events)");

    mc_Q2_RES_Other = new TH1D("mc_Q2_RES_Other","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    mc_Q2_RES_Other->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    mc_Q2_RES_Other->GetYaxis()->SetTitle("N(Events)");

    mc_Q2_DIS = new TH1D("mc_Q2_DIS","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    mc_Q2_DIS->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    mc_Q2_DIS->GetYaxis()->SetTitle("N(Events)");

    mc_Q2_2p2h = new TH1D("mc_Q2_2p2h","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    mc_Q2_2p2h->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    mc_Q2_2p2h->GetYaxis()->SetTitle("N(Events)");

    mc_Q2_Non_RES = new TH1D("mc_Q2_Non_RES","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    mc_Q2_Non_RES->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    mc_Q2_Non_RES->GetYaxis()->SetTitle("N(Events)");

    // ------------------------------------------------------------------------
    // Signal Truth Q2
    // ------------------------------------------------------------------------
    truth_QSq_QE = new TH1D("truth_QSq_QE","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    truth_QSq_QE->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    truth_QSq_QE->GetYaxis()->SetTitle("N(Events)");

    truth_QSq_RES_1232 = new TH1D("truth_QSq_RES_1232","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    truth_QSq_RES_1232->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    truth_QSq_RES_1232->GetYaxis()->SetTitle("N(Events)");

    truth_QSq_RES_1535 = new TH1D("truth_QSq_RES_1535","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    truth_QSq_RES_1535->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    truth_QSq_RES_1535->GetYaxis()->SetTitle("N(Events)");

    truth_QSq_RES_1520 = new TH1D("truth_QSq_RES_1520","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    truth_QSq_RES_1520->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    truth_QSq_RES_1520->GetYaxis()->SetTitle("N(Events)");

    truth_QSq_RES_Other = new TH1D("truth_QSq_RES_Other","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    truth_QSq_RES_Other->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    truth_QSq_RES_Other->GetYaxis()->SetTitle("N(Events)");

    truth_QSq_DIS = new TH1D("truth_QSq_DIS","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    truth_QSq_DIS->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    truth_QSq_DIS->GetYaxis()->SetTitle("N(Events)");

    truth_QSq_2p2h = new TH1D("truth_QSq_2p2h","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    truth_QSq_2p2h->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    truth_QSq_2p2h->GetYaxis()->SetTitle("N(Events)");

    truth_QSq_Non_RES = new TH1D("truth_QSq_Non_RES","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    truth_QSq_Non_RES->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    truth_QSq_Non_RES->GetYaxis()->SetTitle("N(Events)");

    // ------------------------------------------------------------------------
    // Background Truth Q2
    // ------------------------------------------------------------------------
    reco_bckg_QSq_QE = new TH1D("reco_bckg_QSq_QE","Q^{2} for Background Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    reco_bckg_QSq_QE->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    reco_bckg_QSq_QE->GetYaxis()->SetTitle("N(Events)");

    reco_bckg_QSq_RES_1232 = new TH1D("reco_bckg_QSq_RES_1232","Q^{2} for Background Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    reco_bckg_QSq_RES_1232->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    reco_bckg_QSq_RES_1232->GetYaxis()->SetTitle("N(Events)");

    reco_bckg_QSq_RES_1535 = new TH1D("reco_bckg_QSq_RES_1535","Q^{2} for Background Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    reco_bckg_QSq_RES_1535->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    reco_bckg_QSq_RES_1535->GetYaxis()->SetTitle("N(Events)");

    reco_bckg_QSq_RES_1520 = new TH1D("reco_bckg_QSq_RES_1520","Q^{2} for Background Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    reco_bckg_QSq_RES_1520->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    reco_bckg_QSq_RES_1520->GetYaxis()->SetTitle("N(Events)");

    reco_bckg_QSq_RES_Other = new TH1D("reco_bckg_QSq_RES_Other","Q^{2} for Background Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    reco_bckg_QSq_RES_Other->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    reco_bckg_QSq_RES_Other->GetYaxis()->SetTitle("N(Events)");

    reco_bckg_QSq_DIS = new TH1D("reco_bckg_QSq_DIS","Q^{2} for Background Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    reco_bckg_QSq_DIS->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    reco_bckg_QSq_DIS->GetYaxis()->SetTitle("N(Events)");

    reco_bckg_QSq_2p2h = new TH1D("reco_bckg_QSq_2p2h","Q^{2} for Background Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    reco_bckg_QSq_2p2h->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    reco_bckg_QSq_2p2h->GetYaxis()->SetTitle("N(Events)");

    reco_bckg_QSq_Non_RES = new TH1D("reco_bckg_QSq_Non_RES","Q^{2} for Background Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    reco_bckg_QSq_Non_RES->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    reco_bckg_QSq_Non_RES->GetYaxis()->SetTitle("N(Events)");

    reco_bckg_QSq_Coh = new TH1D("reco_bckg_QSq_Coh","Q^{2} for Background Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    reco_bckg_QSq_Coh->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    reco_bckg_QSq_Coh->GetYaxis()->SetTitle("N(Events)");

    // ------------------------------------------------------------------------
    // Signal incomingE
    // ------------------------------------------------------------------------
    mc_incomingE_QE = new TH1D("mc_incomingE_QE","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    mc_incomingE_QE->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    mc_incomingE_QE->GetYaxis()->SetTitle("N(Events)");

    mc_incomingE_RES_1232 = new TH1D("mc_incomingE_RES_1232","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    mc_incomingE_RES_1232->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    mc_incomingE_RES_1232->GetYaxis()->SetTitle("N(Events)");

    mc_incomingE_RES_1535 = new TH1D("mc_incomingE_RES_1535","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    mc_incomingE_RES_1535->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    mc_incomingE_RES_1535->GetYaxis()->SetTitle("N(Events)");

    mc_incomingE_RES_1520 = new TH1D("mc_incomingE_RES_1520","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    mc_incomingE_RES_1520->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    mc_incomingE_RES_1520->GetYaxis()->SetTitle("N(Events)");

    mc_incomingE_RES_Other = new TH1D("mc_incomingE_RES_Other","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    mc_incomingE_RES_Other->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    mc_incomingE_RES_Other->GetYaxis()->SetTitle("N(Events)");

    mc_incomingE_DIS = new TH1D("mc_incomingE_DIS","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    mc_incomingE_DIS->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    mc_incomingE_DIS->GetYaxis()->SetTitle("N(Events)");

    mc_incomingE_2p2h = new TH1D("mc_incomingE_2p2h","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    mc_incomingE_2p2h->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    mc_incomingE_2p2h->GetYaxis()->SetTitle("N(Events)");

    mc_incomingE_Non_RES = new TH1D("mc_incomingE_Non_RES","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    mc_incomingE_Non_RES->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    mc_incomingE_Non_RES->GetYaxis()->SetTitle("N(Events)");

    // ------------------------------------------------------------------------
    // Signal Truth incomingE
    // ------------------------------------------------------------------------
    truth_Enu_QE = new TH1D("truth_Enu_QE","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    truth_Enu_QE->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    truth_Enu_QE->GetYaxis()->SetTitle("N(Events)");

    truth_Enu_RES_1232 = new TH1D("truth_Enu_RES_1232","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    truth_Enu_RES_1232->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    truth_Enu_RES_1232->GetYaxis()->SetTitle("N(Events)");

    truth_Enu_RES_1535 = new TH1D("truth_Enu_RES_1535","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    truth_Enu_RES_1535->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    truth_Enu_RES_1535->GetYaxis()->SetTitle("N(Events)");

    truth_Enu_RES_1520 = new TH1D("truth_Enu_RES_1520","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    truth_Enu_RES_1520->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    truth_Enu_RES_1520->GetYaxis()->SetTitle("N(Events)");

    truth_Enu_RES_Other = new TH1D("truth_Enu_RES_Other","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    truth_Enu_RES_Other->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    truth_Enu_RES_Other->GetYaxis()->SetTitle("N(Events)");

    truth_Enu_DIS = new TH1D("truth_Enu_DIS","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    truth_Enu_DIS->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    truth_Enu_DIS->GetYaxis()->SetTitle("N(Events)");

    truth_Enu_2p2h = new TH1D("truth_Enu_2p2h","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    truth_Enu_2p2h->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    truth_Enu_2p2h->GetYaxis()->SetTitle("N(Events)");

    truth_Enu_Non_RES = new TH1D("truth_Enu_Non_RES","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    truth_Enu_Non_RES->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    truth_Enu_Non_RES->GetYaxis()->SetTitle("N(Events)");

    // ------------------------------------------------------------------------
    // Background Truth incomingE
    // ------------------------------------------------------------------------
    reco_bckg_Enu_QE = new TH1D("reco_bckg_Enu_QE","E_{#nu} for Background Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    reco_bckg_Enu_QE->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    reco_bckg_Enu_QE->GetYaxis()->SetTitle("N(Events)");

    reco_bckg_Enu_RES_1232 = new TH1D("reco_bckg_Enu_RES_1232","E_{#nu} for Background Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    reco_bckg_Enu_RES_1232->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    reco_bckg_Enu_RES_1232->GetYaxis()->SetTitle("N(Events)");

    reco_bckg_Enu_RES_1535 = new TH1D("reco_bckg_Enu_RES_1535","E_{#nu} for Background Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    reco_bckg_Enu_RES_1535->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    reco_bckg_Enu_RES_1535->GetYaxis()->SetTitle("N(Events)");

    reco_bckg_Enu_RES_1520 = new TH1D("reco_bckg_Enu_RES_1520","E_{#nu} for Background Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    reco_bckg_Enu_RES_1520->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    reco_bckg_Enu_RES_1520->GetYaxis()->SetTitle("N(Events)");

    reco_bckg_Enu_RES_Other = new TH1D("reco_bckg_Enu_RES_Other","E_{#nu} for Background Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    reco_bckg_Enu_RES_Other->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    reco_bckg_Enu_RES_Other->GetYaxis()->SetTitle("N(Events)");

    reco_bckg_Enu_DIS = new TH1D("reco_bckg_Enu_DIS","E_{#nu} for Background Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    reco_bckg_Enu_DIS->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    reco_bckg_Enu_DIS->GetYaxis()->SetTitle("N(Events)");

    reco_bckg_Enu_2p2h = new TH1D("reco_bckg_Enu_2p2h","E_{#nu} for Background Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    reco_bckg_Enu_2p2h->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    reco_bckg_Enu_2p2h->GetYaxis()->SetTitle("N(Events)");

    reco_bckg_Enu_Non_RES = new TH1D("reco_bckg_Enu_Non_RES","E_{#nu} for Background Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    reco_bckg_Enu_Non_RES->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    reco_bckg_Enu_Non_RES->GetYaxis()->SetTitle("N(Events)");

    reco_bckg_Enu_Coh = new TH1D("reco_bckg_Enu_Coh","E_{#nu} for Background Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    reco_bckg_Enu_Coh->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    reco_bckg_Enu_Coh->GetYaxis()->SetTitle("N(Events)");

    // ------------------------------------------------------------------------
    // Signal w
    // ------------------------------------------------------------------------
    mc_w_QE = new TH1D("mc_w_QE","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    mc_w_QE->GetXaxis()->SetTitle("W [GeV]");
    mc_w_QE->GetYaxis()->SetTitle("N(Events)");
 
    mc_w_RES_1232 = new TH1D("mc_w_RES_1232","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    mc_w_RES_1232->GetXaxis()->SetTitle("W [GeV]");
    mc_w_RES_1232->GetYaxis()->SetTitle("N(Events)");

    mc_w_RES_1535 = new TH1D("mc_w_RES_1535","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    mc_w_RES_1535->GetXaxis()->SetTitle("W [GeV]");
    mc_w_RES_1535->GetYaxis()->SetTitle("N(Events)");

    mc_w_RES_1520 = new TH1D("mc_w_RES_1520","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    mc_w_RES_1520->GetXaxis()->SetTitle("W [GeV]");
    mc_w_RES_1520->GetYaxis()->SetTitle("N(Events)");

    mc_w_RES_Other = new TH1D("mc_w_RES_Other","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    mc_w_RES_Other->GetXaxis()->SetTitle("W [GeV]");
    mc_w_RES_Other->GetYaxis()->SetTitle("N(Events)");

    mc_w_DIS = new TH1D("mc_w_DIS","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    mc_w_DIS->GetXaxis()->SetTitle("W [GeV]");
    mc_w_DIS->GetYaxis()->SetTitle("N(Events)");

    mc_w_2p2h = new TH1D("mc_w_2p2h","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    mc_w_2p2h->GetXaxis()->SetTitle("W [GeV]");
    mc_w_2p2h->GetYaxis()->SetTitle("N(Events)");

    mc_w_Non_RES = new TH1D("mc_w_Non_RES","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    mc_w_Non_RES->GetXaxis()->SetTitle("W [GeV]");
    mc_w_Non_RES->GetYaxis()->SetTitle("N(Events)");

    // ------------------------------------------------------------------------
    // Signal Truth w
    // ------------------------------------------------------------------------
    truth_w_QE = new TH1D("truth_w_QE","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    truth_w_QE->GetXaxis()->SetTitle("W [GeV]");
    truth_w_QE->GetYaxis()->SetTitle("N(Events)");
 
    truth_w_RES_1232 = new TH1D("truth_w_RES_1232","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    truth_w_RES_1232->GetXaxis()->SetTitle("W [GeV]");
    truth_w_RES_1232->GetYaxis()->SetTitle("N(Events)");

    truth_w_RES_1535 = new TH1D("truth_w_RES_1535","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    truth_w_RES_1535->GetXaxis()->SetTitle("W [GeV]");
    truth_w_RES_1535->GetYaxis()->SetTitle("N(Events)");

    truth_w_RES_1520 = new TH1D("truth_w_RES_1520","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    truth_w_RES_1520->GetXaxis()->SetTitle("W [GeV]");
    truth_w_RES_1520->GetYaxis()->SetTitle("N(Events)");

    truth_w_RES_Other = new TH1D("truth_w_RES_Other","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    truth_w_RES_Other->GetXaxis()->SetTitle("W [GeV]");
    truth_w_RES_Other->GetYaxis()->SetTitle("N(Events)");

    truth_w_DIS = new TH1D("truth_w_DIS","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    truth_w_DIS->GetXaxis()->SetTitle("W [GeV]");
    truth_w_DIS->GetYaxis()->SetTitle("N(Events)");

    truth_w_2p2h = new TH1D("truth_w_2p2h","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    truth_w_2p2h->GetXaxis()->SetTitle("W [GeV]");
    truth_w_2p2h->GetYaxis()->SetTitle("N(Events)");

    truth_w_Non_RES = new TH1D("truth_w_Non_RES","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    truth_w_Non_RES->GetXaxis()->SetTitle("W [GeV]");
    truth_w_Non_RES->GetYaxis()->SetTitle("N(Events)");

    // ------------------------------------------------------------------------
    // Background Truth w
    // ------------------------------------------------------------------------
    reco_bckg_w_QE = new TH1D("reco_bckg_w_QE","W for Background Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    reco_bckg_w_QE->GetXaxis()->SetTitle("W [GeV]");
    reco_bckg_w_QE->GetYaxis()->SetTitle("N(Events)");
 
    reco_bckg_w_RES_1232 = new TH1D("reco_bckg_w_RES_1232","W for Background Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    reco_bckg_w_RES_1232->GetXaxis()->SetTitle("W [GeV]");
    reco_bckg_w_RES_1232->GetYaxis()->SetTitle("N(Events)");

    reco_bckg_w_RES_1535 = new TH1D("reco_bckg_w_RES_1535","W for Background Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    reco_bckg_w_RES_1535->GetXaxis()->SetTitle("W [GeV]");
    reco_bckg_w_RES_1535->GetYaxis()->SetTitle("N(Events)");

    reco_bckg_w_RES_1520 = new TH1D("reco_bckg_w_RES_1520","W for Background Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    reco_bckg_w_RES_1520->GetXaxis()->SetTitle("W [GeV]");
    reco_bckg_w_RES_1520->GetYaxis()->SetTitle("N(Events)");

    reco_bckg_w_RES_Other = new TH1D("reco_bckg_w_RES_Other","W for Background Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    reco_bckg_w_RES_Other->GetXaxis()->SetTitle("W [GeV]");
    reco_bckg_w_RES_Other->GetYaxis()->SetTitle("N(Events)");

    reco_bckg_w_DIS = new TH1D("reco_bckg_w_DIS","W for Background Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    reco_bckg_w_DIS->GetXaxis()->SetTitle("W [GeV]");
    reco_bckg_w_DIS->GetYaxis()->SetTitle("N(Events)");

    reco_bckg_w_2p2h = new TH1D("reco_bckg_w_2p2h","W for Background Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    reco_bckg_w_2p2h->GetXaxis()->SetTitle("W [GeV]");
    reco_bckg_w_2p2h->GetYaxis()->SetTitle("N(Events)");

    reco_bckg_w_Non_RES = new TH1D("reco_bckg_w_Non_RES","W for Background Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    reco_bckg_w_Non_RES->GetXaxis()->SetTitle("W [GeV]");
    reco_bckg_w_Non_RES->GetYaxis()->SetTitle("N(Events)");

    reco_bckg_w_Coh = new TH1D("reco_bckg_w_Coh","W for Background Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    reco_bckg_w_Coh->GetXaxis()->SetTitle("W [GeV]");
    reco_bckg_w_Coh->GetYaxis()->SetTitle("N(Events)");

    // ------------------------------------------------------------------------
    // Signal w
    // ------------------------------------------------------------------------
    reco_w_QE = new TH1D("reco_w_QE","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    reco_w_QE->GetXaxis()->SetTitle("W [GeV]");
    reco_w_QE->GetYaxis()->SetTitle("N(Events)");
 
    reco_w_RES_1232 = new TH1D("reco_w_RES_1232","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    reco_w_RES_1232->GetXaxis()->SetTitle("W [GeV]");
    reco_w_RES_1232->GetYaxis()->SetTitle("N(Events)");

    reco_w_RES_1535 = new TH1D("reco_w_RES_1535","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    reco_w_RES_1535->GetXaxis()->SetTitle("W [GeV]");
    reco_w_RES_1535->GetYaxis()->SetTitle("N(Events)");

    reco_w_RES_1520 = new TH1D("reco_w_RES_1520","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    reco_w_RES_1520->GetXaxis()->SetTitle("W [GeV]");
    reco_w_RES_1520->GetYaxis()->SetTitle("N(Events)");

    reco_w_RES_Other = new TH1D("reco_w_RES_Other","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    reco_w_RES_Other->GetXaxis()->SetTitle("W [GeV]");
    reco_w_RES_Other->GetYaxis()->SetTitle("N(Events)");

    reco_w_DIS = new TH1D("reco_w_DIS","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    reco_w_DIS->GetXaxis()->SetTitle("W [GeV]");
    reco_w_DIS->GetYaxis()->SetTitle("N(Events)");

    reco_w_2p2h = new TH1D("reco_w_2p2h","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    reco_w_2p2h->GetXaxis()->SetTitle("W [GeV]");
    reco_w_2p2h->GetYaxis()->SetTitle("N(Events)");

    reco_w_Non_RES = new TH1D("reco_w_Non_RES","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    reco_w_Non_RES->GetXaxis()->SetTitle("W [GeV]");
    reco_w_Non_RES->GetYaxis()->SetTitle("N(Events)");

    // Random Number Generator Histograms 
    normal_rand_numbers = new TH1D( "normal_rand_numbers","Normal Random Numbers",50,-3.0,3.0);
    normal_rand_numbers->GetXaxis()->SetTitle("Normal Random Numbers");
    normal_rand_numbers->GetYaxis()->SetTitle("N(Numbers)");

    em_shift_rand_numbers = new TH1D( "em_shift_rand_numbers","EM Energy Scale Shift Random Numbers",50,-0.1,0.1);
    em_shift_rand_numbers->GetXaxis()->SetTitle("EM Energy Scale Shift");
    em_shift_rand_numbers->GetYaxis()->SetTitle("N(Numbers)");

    muonP_shift_rand_numbers = new TH1D( "muonP_shift_rand_numbers","Muon Momentum Shift Random Numbers",50,-0.1,0.1);
    muonP_shift_rand_numbers->GetXaxis()->SetTitle("Muon Momentum Shift");
    muonP_shift_rand_numbers->GetYaxis()->SetTitle("N(Numbers)");

    muon_theta_shift_rand_numbers = new TH1D( "muon_theta_shift_rand_numbers","Muon Theta Shift Random Numbers",50,-0.02,0.02);
    muon_theta_shift_rand_numbers->GetXaxis()->SetTitle("Muon Theta Shift");
    muon_theta_shift_rand_numbers->GetYaxis()->SetTitle("N(Numbers)");

    Birks_shift_rand_numbers = new TH1D( "Birks_shift_rand_numbers","Proton Energy Birks Shift Random Numbers",50,-0.02,0.02);
    Birks_shift_rand_numbers->GetXaxis()->SetTitle("Proton Energy Birks Shift");
    Birks_shift_rand_numbers->GetYaxis()->SetTitle("N(Numbers)");

    Err_NeutronResponse = new TH1D( "Err_NeutronResponse","Neutron Response Error",50,-1,1);
    Err_NeutronResponse->GetXaxis()->SetTitle("Error used as (wgt = 1 #pm error)");
    Err_NeutronResponse->GetYaxis()->SetTitle("N(Events)");

    Err_PionResponse = new TH1D( "Err_PionResponse","Pion Response Error",50,-1,1);
    Err_PionResponse->GetXaxis()->SetTitle("Error used as (wgt = 1 #pm error)");
    Err_PionResponse->GetYaxis()->SetTitle("N(Events)");

    Err_MuonTracking = new TH1D( "Err_MuonTracking","Muon Tracking Error",50,-1,1);
    Err_MuonTracking->GetXaxis()->SetTitle("Error used as (wgt = 1 #pm error)");
    Err_MuonTracking->GetYaxis()->SetTitle("N(Events)");


}

void CCProtonPi0_Interaction::writeHistograms()
{
    std::cout<<">> Writing "<<rootDir<<std::endl;
   
    f->cd();

    for (int i = 0; i < nHistograms; i++){
        CV_weight[i]->Write();
        CV_weight_2p2h[i]->Write();
        CV_weight_Delta[i]->Write();
        CV_weight_CCRES[i]->Write();
        CV_weight_NonRes1pi[i]->Write();
        err_2p2h[i]->Write();
        genie_wgt_Theta_Delta2Npi[i]->Write();
        updated_wgt_Theta_Delta2Npi[i]->Write();
        genie_wgt_MaRES[i]->Write();
        updated_wgt_MaRES[i]->Write();
        genie_wgt_MvRES[i]->Write();
        updated_wgt_MvRES[i]->Write();
        genie_wgt_Rvn1pi[i]->Write();
        updated_wgt_Rvn1pi[i]->Write();

        // Event Kinematics
        Enu[i]->Write();
        Enu_1Track[i]->Write();
        Enu_2Track[i]->Write();
        QSq[i]->Write();
        QSq_1Track[i]->Write();
        QSq_2Track[i]->Write();
        WSq[i]->Write();
        WSq_1Track[i]->Write();
        WSq_2Track[i]->Write();
        W[i]->Write();
        W_1Track[i]->Write();
        W_2Track[i]->Write();
       
        // Vertex Energy
        vertex_energy_1Track[i]->Write();
        vertex_energy_2Track[i]->Write();
        vertex_evis_1Track[i]->Write();
        vertex_evis_2Track[i]->Write();

        // Other Event Parameters 
        deltaInvMass[i]->Write();

        // Extra Energy
        extra_leftover_energy_1Track[i]->Write();
        extra_muon_energy_1Track[i]->Write();
        extra_rejected_energy_1Track[i]->Write();
        extra_total_energy_1Track[i]->Write();

        extra_leftover_energy_2Track[i]->Write();
        extra_muon_energy_2Track[i]->Write();
        extra_rejected_energy_2Track[i]->Write();
        extra_total_energy_2Track[i]->Write();
    }
   
    for (int i = 0; i < 10; ++i){
        W_p_pi0[i]->Write();
        W_All[i]->Write();
        W_1[i]->Write();
        W_2[i]->Write();
        QSq_All[i]->Write();
        QSq_LowEnu[i]->Write();
        QSq_HighEnu[i]->Write();
    }

    for (int i = 0; i <=100; ++i){
        QSq_LowMaRES[i]->Write();
        QSq_HighMaRES[i]->Write();
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
 
    W_all->Write();
    W_mc_truth_signal->Write();
    W_mc_reco_all->Write();
    W_mc_reco_signal->Write();
    W_mc_reco_bckg->Write();
    W_response->Write();
    
    // Short Proton
    nProtons->Write();
    proton_true_P_1Track->Write();
    proton_true_KE_1Track->Write();
    proton_true_theta_1Track->Write();
    
    // Ejected Nucleons
    n_ejected_nucleons_1Track->Write();
    n_ejected_nucleons_2Track->Write();

    //Extra Energy
    h_extra_muon_energy->Write();
    h_extra_rejected_energy->Write();
    
    // QSq Truth, Error, Difference
    WSq_QSq_Diff->Write();
    QSq_All_response->Write();
    QSq_1Track_response->Write();
    QSq_2Track_response->Write();
    
    QSq_Error->Write();
    QSq_1Track_Error->Write();
    QSq_2Track_Error->Write();

    QSq_Diff->Write();
    QSq_1Track_Diff->Write();
    QSq_2Track_Diff->Write();
   
    // W Error, Difference
    W_Diff->Write();
    W_Error->Write();

    // Neutrino Energy: Error, Difference
    Enu_All_response->Write();
    Enu_1Track_response->Write();
    Enu_2Track_response->Write();
    
    Enu_Error->Write();
    Enu_1Track_Error->Write();
    Enu_2Track_Error->Write();

    Enu_Diff->Write();
    Enu_1Track_Diff->Write();
    Enu_2Track_Diff->Write();

    // Signal Q2
    mc_Q2_QE->Write();
    mc_Q2_RES_1232->Write();
    mc_Q2_RES_1535->Write();
    mc_Q2_RES_1520->Write();
    mc_Q2_RES_Other->Write();
   
    mc_Q2_DIS->Write();
    mc_Q2_2p2h->Write();
    mc_Q2_Non_RES->Write();
 
    truth_QSq_QE->Write();
    truth_QSq_RES_1232->Write();
    truth_QSq_RES_1535->Write();
    truth_QSq_RES_1520->Write();
    truth_QSq_RES_Other->Write();
   
    truth_QSq_DIS->Write();
    truth_QSq_2p2h->Write();
    truth_QSq_Non_RES->Write();
 
    reco_bckg_QSq_QE->Write();
    reco_bckg_QSq_RES_1232->Write();
    reco_bckg_QSq_RES_1535->Write();
    reco_bckg_QSq_RES_1520->Write();
    reco_bckg_QSq_RES_Other->Write();
   
    reco_bckg_QSq_DIS->Write();
    reco_bckg_QSq_2p2h->Write();
    reco_bckg_QSq_Non_RES->Write();
    reco_bckg_QSq_Coh->Write();

    // Signal incomingE
    mc_incomingE_QE->Write();
    mc_incomingE_RES_1232->Write();
    mc_incomingE_RES_1535->Write();
    mc_incomingE_RES_1520->Write();
    mc_incomingE_RES_Other->Write();
   
    mc_incomingE_DIS->Write();
    mc_incomingE_2p2h->Write();
    mc_incomingE_Non_RES->Write();
 
    truth_Enu_QE->Write();
    truth_Enu_RES_1232->Write();
    truth_Enu_RES_1535->Write();
    truth_Enu_RES_1520->Write();
    truth_Enu_RES_Other->Write();
   
    truth_Enu_DIS->Write();
    truth_Enu_2p2h->Write();
    truth_Enu_Non_RES->Write();
 
    reco_bckg_Enu_QE->Write();
    reco_bckg_Enu_RES_1232->Write();
    reco_bckg_Enu_RES_1535->Write();
    reco_bckg_Enu_RES_1520->Write();
    reco_bckg_Enu_RES_Other->Write();
   
    reco_bckg_Enu_DIS->Write();
    reco_bckg_Enu_2p2h->Write();
    reco_bckg_Enu_Non_RES->Write();
    reco_bckg_Enu_Coh->Write();

    // Signal w
    mc_w_QE->Write();
    mc_w_RES_1232->Write();
    mc_w_RES_1535->Write();
    mc_w_RES_1520->Write();
    mc_w_RES_Other->Write();
   
    mc_w_DIS->Write();
    mc_w_2p2h->Write();
    mc_w_Non_RES->Write();

    truth_w_QE->Write();
    truth_w_RES_1232->Write();
    truth_w_RES_1535->Write();
    truth_w_RES_1520->Write();
    truth_w_RES_Other->Write();
   
    truth_w_DIS->Write();
    truth_w_2p2h->Write();
    truth_w_Non_RES->Write();

    reco_bckg_w_QE->Write();
    reco_bckg_w_RES_1232->Write();
    reco_bckg_w_RES_1535->Write();
    reco_bckg_w_RES_1520->Write();
    reco_bckg_w_RES_Other->Write();
   
    reco_bckg_w_DIS->Write();
    reco_bckg_w_2p2h->Write();
    reco_bckg_w_Non_RES->Write();
    reco_bckg_w_Coh->Write();

    // Signal reco w
    reco_w_QE->Write();
    reco_w_RES_1232->Write();
    reco_w_RES_1535->Write();
    reco_w_RES_1520->Write();
    reco_w_RES_Other->Write();
   
    reco_w_DIS->Write();
    reco_w_2p2h->Write();
    reco_w_Non_RES->Write();

    Polarization_data->Write();
    Polarization_mc->Write();
    DeltaTransverse_data->Write();
    DeltaTransverse_mc->Write();
    DeltaTransverse_mc_res->Write();

    normal_rand_numbers->Write();
    em_shift_rand_numbers->Write();
    muonP_shift_rand_numbers->Write();
    muon_theta_shift_rand_numbers->Write();
    Birks_shift_rand_numbers->Write();

    Err_NeutronResponse->Write();
    Err_PionResponse->Write();
    Err_MuonTracking->Write();

    f->Close();
}



#endif
