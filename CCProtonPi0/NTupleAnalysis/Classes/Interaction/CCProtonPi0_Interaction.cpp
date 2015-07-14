/*
    See CCProtonPi0_Interaction.h header for Class Information
*/
#ifndef CCProtonPi0_Interaction_cpp
#define CCProtonPi0_Interaction_cpp

#include "CCProtonPi0_Interaction.h"

using namespace std;

CCProtonPi0_Interaction::CCProtonPi0_Interaction(int nMode, bool isMC) : CCProtonPi0_NTupleAnalysis(nMode)
{
    cout<<"Initializing CCProtonPi0_Interaction"<<endl;
    
    if(nMode == 0){
        cout<<"\tNTuple Reduce Mode -- Will not create ROOT Files"<<endl;
    }else{
        if (isMC) rootDir = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed +  branchDir + "Interaction.root";
        else rootDir = Folder_List::rootOut + Folder_List::Data + Folder_List::analyzed + branchDir + "Interaction.root";
        
        cout<<"\tRoot File: "<<rootDir<<endl;
 
        // Create Root File 
        f = new TFile(rootDir.c_str(),"RECREATE");
        
        initHistograms();
    }
    
    cout<<"Done!"<<endl;
}


void CCProtonPi0_Interaction::initHistograms()
{
 
    proton_p = new TH1D( "proton_p","Proton Momentum Original",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    proton_p->GetXaxis()->SetTitle("(Reco-True)/True");
    proton_p->GetYaxis()->SetTitle("Reco NOT Modified");
    
    proton_p_shifted = new TH1D( "proton_p_shifted","Proton Momentum Shifted",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    proton_p_shifted->GetXaxis()->SetTitle("(Reco-True)/True");
    proton_p_shifted->GetYaxis()->SetTitle("Reco Modified by a factor 1.050");
 
    status_Pi0 = new TH1D( "status_Pi0","Pi0 Status",binList.particleStatus.get_nBins(), binList.particleStatus.get_min(), binList.particleStatus.get_max() );
    status_Pi0->GetXaxis()->SetTitle("GENIE Status Codes");
    status_Pi0->GetYaxis()->SetTitle("N(Events)");
    
    status_Pi0_Mother = new TH1D( "status_Pi0_Mother","Pi0 Mother Status",binList.particleStatus.get_nBins(), binList.particleStatus.get_min(), binList.particleStatus.get_max() );
    status_Pi0_Mother->GetXaxis()->SetTitle("GENIE Status Codes");
    status_Pi0_Mother->GetYaxis()->SetTitle("N(Events)");
    
    status_Pi0_GrandMother = new TH1D( "status_Pi0_GrandMother","Pi0 GrandMother Status",binList.particleStatus.get_nBins(), binList.particleStatus.get_min(), binList.particleStatus.get_max() );
    status_Pi0_GrandMother->GetXaxis()->SetTitle("GENIE Status Codes");
    status_Pi0_GrandMother->GetYaxis()->SetTitle("N(Events)");
    
    Enu_1Track_mc = new TH1D( "Enu_1Track_mc","True Beam Energy - 1 Track",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    Enu_1Track_mc->GetXaxis()->SetTitle("True Beam Energy [GeV]");
    Enu_1Track_mc->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.beamE.get_width()));
    
    Enu_1Track_reco = new TH1D( "Enu_1Track_reco","Reconstructed Beam Energy - 1 Track",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    Enu_1Track_reco->GetXaxis()->SetTitle("Reconstructed Beam Energy [GeV]");
    Enu_1Track_reco->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.beamE.get_width()));
    
    Enu_1Track_error = new TH1D( "Enu_1Track_error","Error on Beam Energy - 1 Track",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    Enu_1Track_error->GetXaxis()->SetTitle("(True-Reco) / True");
    Enu_1Track_error->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.error.get_width()));
    
    Enu_1Track_reco_mc = new TH2D( "Enu_1Track_reco_mc","True vs Reconstructed Beam Energy - 1 Track",
    binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max(),
    binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max());
    Enu_1Track_reco_mc->GetXaxis()->SetTitle("Reconstructed Beam Energy [GeV]");
    Enu_1Track_reco_mc->GetYaxis()->SetTitle("True Beam Energy [GeV]");
    
    Enu_2Track_mc = new TH1D( "Enu_2Track_mc","True Beam Energy - 2 Track",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    Enu_2Track_mc->GetXaxis()->SetTitle("True Beam Energy [GeV]");
    Enu_2Track_mc->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.beamE.get_width()));
    
    Enu_2Track_reco = new TH1D( "Enu_2Track_reco","Reconstructed Beam Energy - 2 Track",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    Enu_2Track_reco->GetXaxis()->SetTitle("Reconstructed Beam Energy [GeV]");
    Enu_2Track_reco->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.beamE.get_width()));
    
    Enu_2Track_error = new TH1D( "Enu_2Track_error","Error on Beam Energy - 2 Track",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    Enu_2Track_error->GetXaxis()->SetTitle("(True-Reco) / True");
    Enu_2Track_error->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.error.get_width()));
    
    Enu_2Track_reco_mc = new TH2D( "Enu_2Track_reco_mc","True vs Reconstructed Beam Energy - 2 Track",
    binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max(),
    binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max());
    Enu_2Track_reco_mc->GetXaxis()->SetTitle("Reconstructed Beam Energy [GeV]");
    Enu_2Track_reco_mc->GetYaxis()->SetTitle("True Beam Energy [GeV]");

    Enu_Cal_mc = new TH1D( "Enu_Cal_mc","True Beam Energy",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    Enu_Cal_mc->GetXaxis()->SetTitle("True Beam Energy [GeV]");
    Enu_Cal_mc->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.beamE.get_width()));
    
    Enu_Cal_reco = new TH1D( "Enu_Cal_reco","Reconstructed Calorimetric Beam Energy",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    Enu_Cal_reco->GetXaxis()->SetTitle("Reconstructed Calorimetric  Beam Energy [GeV]");
    Enu_Cal_reco->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.beamE.get_width()));
    
    Enu_Cal_error = new TH1D( "Enu_Cal_error","Error on Calorimetric Beam Energy",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    Enu_Cal_error->GetXaxis()->SetTitle("(True-Reco) / True");
    Enu_Cal_error->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.error.get_width()));
    
    Enu_Cal_reco_mc = new TH2D( "Enu_Cal_reco_mc","True vs Reconstructed Calorimetric Beam Energy",
                                   binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max(),
                                   binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max());
    Enu_Cal_reco_mc->GetXaxis()->SetTitle("Reconstructed Calorimetric Beam Energy [GeV]");
    Enu_Cal_reco_mc->GetYaxis()->SetTitle("True Beam Energy [GeV]");
    
    Enu_1Track_Enu_Cal = new TH2D( "Enu_1Track_Enu_Cal","Calorimetric Beam Energy vs Beam Energy from Tp",
                                      binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max(),
                                      binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max());
    Enu_1Track_Enu_Cal->GetXaxis()->SetTitle("Reconstructed Beam Energy from Tp [GeV]");
    Enu_1Track_Enu_Cal->GetYaxis()->SetTitle("Reconstructed Calorimetric Beam Energy [GeV]");
    
    q2_mc = new TH1D( "q2_mc","True Q^{2}",binList.q2.get_nBins(), binList.q2.get_min(), binList.q2.get_max() );
    q2_mc->GetXaxis()->SetTitle("True Q^{2} [GeV^{2}]");
    q2_mc->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.q2.get_width()));
    
    q2_reco = new TH1D( "q2_reco","Reconstructed Q^{2}",binList.q2.get_nBins(), binList.q2.get_min(), binList.q2.get_max() );
    q2_reco->GetXaxis()->SetTitle("Reconstructed Q^{2} [GeV^{2}]");
    q2_reco->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.q2.get_width()));
    
    q2_error = new TH1D( "q2_error","Error on Q^{2}",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    q2_error->GetXaxis()->SetTitle("(True-Reco) / True");
    q2_error->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.error.get_width()));
    
    q2_reco_mc = new TH2D( "q2_reco_mc","True vs Reconstructed Q^{2}",
                            binList.q2.get_nBins(), binList.q2.get_min(), binList.q2.get_max(),
                            binList.q2.get_nBins(), binList.q2.get_min(), binList.q2.get_max());
    q2_reco_mc->GetXaxis()->SetTitle("Reconstructed Q^{2} [GeV^{2}]");
    q2_reco_mc->GetYaxis()->SetTitle("True Q^{2} [GeV^{2}]");
    
    w_mc = new TH1D( "w_mc","True W",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max() );
    w_mc->GetXaxis()->SetTitle("True W [GeV]");
    w_mc->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.w.get_width()));
    
    w_reco = new TH1D( "w_reco","Reconstructed W",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max() );
    w_reco->GetXaxis()->SetTitle("Reconstructed W [GeV]");
    w_reco->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.w.get_width()));
    
    wSq_reco = new TH1D( "wSq_reco","Reconstructed W^{2}",binList.wSq.get_nBins(), binList.wSq.get_min(), binList.wSq.get_max() );
    wSq_reco->GetXaxis()->SetTitle("Reconstructed W^{2} [GeV^{2}]");
    wSq_reco->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.wSq.get_width()));
    
    w_error = new TH1D( "w_error","Error on W",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    w_error->GetXaxis()->SetTitle("(True-Reco) / True");
    w_error->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.error.get_width()));
    
    w_reco_mc = new TH2D( "w_reco_mc","True vs Reconstructed W",
                           binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max(),
                           binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max());
    w_reco_mc->GetXaxis()->SetTitle("Reconstructed W [GeV]");
    w_reco_mc->GetYaxis()->SetTitle("True W [GeV]");
                           
    mc_w_DIS = new TH1D( "mc_w_DIS","True W for DIS",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max() );
    mc_w_DIS->GetXaxis()->SetTitle("True W for DIS [GeV]");
    mc_w_DIS->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.w.get_width()));
    
    mc_w_RES = new TH1D( "mc_w_RES","True W for RES",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max() );
    mc_w_RES->GetXaxis()->SetTitle("True W for RES [GeV]");
    mc_w_RES->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.w.get_width()));
    
    mc_w_CCQE = new TH1D( "mc_w_CCQE","True W for CCQE",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max() );
    mc_w_CCQE->GetXaxis()->SetTitle("True W for CCQE [GeV]");
    mc_w_CCQE->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.w.get_width()));
    
    final_mc_w_DIS = new TH1D( "final_mc_w_DIS","True W for DIS",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max() );
    final_mc_w_DIS->GetXaxis()->SetTitle("True W for DIS [GeV]");
    final_mc_w_DIS->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.w.get_width()));
    
    final_mc_w_RES = new TH1D( "final_mc_w_RES","True W for RES",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max() );
    final_mc_w_RES->GetXaxis()->SetTitle("True W for RES [GeV]");
    final_mc_w_RES->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.w.get_width()));
    
    final_mc_w_CCQE = new TH1D( "final_mc_w_CCQE","True W for CCQE",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max() );
    final_mc_w_CCQE->GetXaxis()->SetTitle("True W for CCQE [GeV]");
    final_mc_w_CCQE->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.w.get_width()));
    
       
    // Interaction Type
    int_channel = new TH1D( "int_channel","Interaction Channel",binList.int_channel.get_nBins(), binList.int_channel.get_min(), binList.int_channel.get_max() );
    int_channel->GetXaxis()->SetTitle("1 = QE, 2 = Resonant, 3 = DIS, 4 = Coh pi");
    int_channel->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.int_channel.get_width()));
   
    nProngs_hist = new TH1D( "nProngs_hist","Number of Prongs",binList.multiplicity.get_nBins(), binList.multiplicity.get_min(), binList.multiplicity.get_max() );
    nProngs_hist->GetXaxis()->SetTitle("Number of Prongs");
    nProngs_hist->GetYaxis()->SetTitle("Candidates");
    
    // Vertex
    vertex_count = new TH1D( "vertex_count","Number of Vertices",binList.objectCount.get_nBins(), binList.objectCount.get_min(), binList.objectCount.get_max() );
    vertex_count->GetXaxis()->SetTitle("Number of Vertices");
    vertex_count->GetYaxis()->SetTitle("N(Events)");
    
    vertex_x_y_true = new TH2D( "vertex_x_y_true","True Vertex X vs Y",   binList.vertex_x_y.get_nBins(), binList.vertex_x_y.get_min(), binList.vertex_x_y.get_max(),
    binList.vertex_x_y.get_nBins(), binList.vertex_x_y.get_min(), binList.vertex_x_y.get_max());
    vertex_x_y_true->GetXaxis()->SetTitle("True Vertex X [mm]");
    vertex_x_y_true->GetYaxis()->SetTitle("True Vertex Y [mm]");
    
    vertex_x_y_reco = new TH2D( "vertex_x_y_reco","Reconstructed Vertex X vs Y",   binList.vertex_x_y.get_nBins(), binList.vertex_x_y.get_min(), binList.vertex_x_y.get_max(),
    binList.vertex_x_y.get_nBins(), binList.vertex_x_y.get_min(), binList.vertex_x_y.get_max());
    vertex_x_y_reco->GetXaxis()->SetTitle("Reconstructed Vertex X [mm]");
    vertex_x_y_reco->GetYaxis()->SetTitle("Reconstructed Vertex Y [mm]");
    
    vertex_z_true = new TH1D( "vertex_z_true","True Vertex Z",binList.vertex_z.get_nBins(), binList.vertex_z.get_min(), binList.vertex_z.get_max() );
    vertex_z_true->GetXaxis()->SetTitle("z = 4293 Target, #bf{z = 5810 Interaction Region}, z = 8614 ECAL, z = 9088 HCAL");
    vertex_z_true->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.vertex_z.get_width()));
    
    vertex_z_reco = new TH1D( "vertex_z_reco","Reconstructed Vertex Z",binList.vertex_z.get_nBins(), binList.vertex_z.get_min(), binList.vertex_z.get_max() );
    vertex_z_reco->GetXaxis()->SetTitle("z = 4293 Target, #bf{z = 5810 Interaction Region}, z = 8614 ECAL, z = 9088 HCAL");
    vertex_z_reco->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.vertex_z.get_width()));
    
    vertex_z_error = new TH1D( "vertex_z_error","Error on Vertex Z",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    vertex_z_error->GetXaxis()->SetTitle("(True-Reco) / True");
    vertex_z_error->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.error.get_width()));
    
    vertex_z_reco_mc = new TH2D( "vertex_z_reco_mc","True vs Reconstructed Vertex Z",
    binList.vertex_z.get_nBins(), binList.vertex_z.get_min(), binList.vertex_z.get_max(),
    binList.vertex_z.get_nBins(), binList.vertex_z.get_min(), binList.vertex_z.get_max());
    vertex_z_reco_mc->GetXaxis()->SetTitle("Reconstructed Vertex Z [mm]");
    vertex_z_reco_mc->GetYaxis()->SetTitle("True Vertex Z [mm]");
    
    deltaInvMass_reco = new TH1D( "deltaInvMass_reco","Reconstructed Delta+ Invariant Mass",binList.deltaInvMass.get_nBins(), binList.deltaInvMass.get_min(), binList.deltaInvMass.get_max() );
    deltaInvMass_reco->GetXaxis()->SetTitle("Reconstructed Delta+ Invariant Mass [MeV]");
    deltaInvMass_reco->GetYaxis()->SetTitle(Form("Candidates / %3.2f [MeV] ",binList.deltaInvMass.get_width()));
    
    deltaInvMass_mc= new TH1D( "deltaInvMass_mc","True Delta+ Invariant Mass",binList.deltaInvMass.get_nBins(), binList.deltaInvMass.get_min(), binList.deltaInvMass.get_max() );
    deltaInvMass_mc->GetXaxis()->SetTitle("True Delta+ Invariant Mass [MeV]");
    deltaInvMass_mc->GetYaxis()->SetTitle(Form("Candidates / %3.2f [MeV] ",binList.deltaInvMass.get_width()));
    
    deltaInvMass_error = new TH1D( "deltaInvMass_error","Delta+ Invariant Mass Error",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    deltaInvMass_error->GetXaxis()->SetTitle("(True-Reco) / True");
    deltaInvMass_error->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.error.get_width()));
    
    deltaInvMass_reco_mc = new TH2D( "deltaInvMass_reco_mc","True vs Reconstructed Delta+ Invariant Mass",
    binList.deltaInvMass.get_nBins(), binList.deltaInvMass.get_min(), binList.deltaInvMass.get_max(),
    binList.deltaInvMass.get_nBins(), binList.deltaInvMass.get_min(), binList.deltaInvMass.get_max());
    deltaInvMass_reco_mc->GetXaxis()->SetTitle("Reconstructed Delta+ Invariant Mass [MeV]");
    deltaInvMass_reco_mc->GetYaxis()->SetTitle("True Delta+ Invariant Mass [MeV]");
    
    pFilter_Status = new TH1D( "pFilter_Status","Prefilter() Result",binList.preFilter_Status.get_nBins(), binList.preFilter_Status.get_min(), binList.preFilter_Status.get_max() );
    pFilter_Status->GetXaxis()->SetTitle("0 = Passes Filter, 1 = Target Filter, 2 = Max Other, 3 = Min Other");
    pFilter_Status->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.preFilter_Status.get_width()));
    
    pFilter_RejectedEnergy = new TH1D( "pFilter_RejectedEnergy","Rejected Energy by preFilter()",binList.preFilter_RejectedEnergy.get_nBins(), binList.preFilter_RejectedEnergy.get_min(), binList.preFilter_RejectedEnergy.get_max() );
    pFilter_RejectedEnergy->GetXaxis()->SetTitle("Rejected Energy by preFilter()");
    pFilter_RejectedEnergy->GetYaxis()->SetTitle(Form("Candidates / %3.2f [MeV]",binList.preFilter_RejectedEnergy.get_width()));
    
    // Cluster Energy after Pi0 Reconstruction
    E_Unused_afterReco = new TH1D( "E_Unused_afterReco","Unused Cluster Energy after Reconstruction",binList.UnusedE.get_nBins(), binList.UnusedE.get_min(), binList.UnusedE.get_max() );
    E_Unused_afterReco->GetXaxis()->SetTitle("Unused Cluster Energy after Reconstruction [MeV]");
    E_Unused_afterReco->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.UnusedE.get_width()));
    
    E_Used_afterReco = new TH1D( "E_Used_afterReco","Used Cluster Energy after Reconstruction",binList.UsedE.get_nBins(), binList.UsedE.get_min(), binList.UsedE.get_max() );
    E_Used_afterReco->GetXaxis()->SetTitle("Used Cluster Energy after Reconstruction [MeV]");
    E_Used_afterReco->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.UsedE.get_width()));
    
    // Cluster Timing
    time_AllClusters = new TH1D( "time_AllClusters","Cluster Time",binList.time.get_nBins(), binList.time.get_min(), binList.time.get_max() );
    time_AllClusters->GetXaxis()->SetTitle("Cluster Time");
    time_AllClusters->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.time.get_width()));
    
    // Total Final State Particle Energy
    total_E = new TH1D( "total_E","Total FS Particle Energy",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    total_E->GetXaxis()->SetTitle("Total FS Particle Energy [GeV]");
    total_E->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.beamE.get_width()));
    
    total_E_neutrinoE = new TH2D( "total_E_neutrinoE","Neutrino Energy vs Total FS Particle Energy",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max(),binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    total_E_neutrinoE->GetXaxis()->SetTitle("Total FS Particle Energy [GeV]");
    total_E_neutrinoE->GetYaxis()->SetTitle("Neutrino Energy [GeV]");
}

void CCProtonPi0_Interaction::write_RootFile()
{
    cout<<">> Writing "<<rootDir<<endl;
    f->Write();
}



#endif
