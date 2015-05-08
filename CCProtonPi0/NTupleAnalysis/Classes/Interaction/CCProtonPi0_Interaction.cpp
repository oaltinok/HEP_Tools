/*
    See CCProtonPi0_Interaction.h header for Class Information
*/
#ifndef CCProtonPi0_Interaction_cpp
#define CCProtonPi0_Interaction_cpp

#include "CCProtonPi0_Interaction.h"

using namespace std;

CCProtonPi0_Interaction::CCProtonPi0_Interaction(int nMode) : CCProtonPi0_NTupleAnalysis(nMode)
{
    cout<<"Initializing CCProtonPi0_Interaction"<<endl;
    
    rootDir = Folder_List::output + Folder_List::rootOut + branchDir + "Interaction.root";
    
    cout<<"\tRoot File: "<<rootDir<<endl;
 
    // Create Root File 
    f = new TFile(rootDir.c_str(),"RECREATE");
    
    initHistograms();
  
    cout<<"Done!"<<endl;
}

void CCProtonPi0_Interaction::initHistograms()
{
    status_Pi0 = new TH1D( "status_Pi0","Pi0 Status",binList.particleStatus.get_nBins(), binList.particleStatus.get_min(), binList.particleStatus.get_max() );
    status_Pi0->GetXaxis()->SetTitle("GENIE Status Codes");
    status_Pi0->GetYaxis()->SetTitle("N(Events)");
    
    status_Pi0_Mother = new TH1D( "status_Pi0_Mother","Pi0 Mother Status",binList.particleStatus.get_nBins(), binList.particleStatus.get_min(), binList.particleStatus.get_max() );
    status_Pi0_Mother->GetXaxis()->SetTitle("GENIE Status Codes");
    status_Pi0_Mother->GetYaxis()->SetTitle("N(Events)");
    
    status_Pi0_GrandMother = new TH1D( "status_Pi0_GrandMother","Pi0 GrandMother Status",binList.particleStatus.get_nBins(), binList.particleStatus.get_min(), binList.particleStatus.get_max() );
    status_Pi0_GrandMother->GetXaxis()->SetTitle("GENIE Status Codes");
    status_Pi0_GrandMother->GetYaxis()->SetTitle("N(Events)");
    
    beamEnergy_mc = new TH1D( "beamEnergy_mc","True Beam Energy",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    beamEnergy_mc->GetXaxis()->SetTitle("True Beam Energy [GeV]");
    beamEnergy_mc->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.beamE.get_width()));
    
    beamEnergy_reco = new TH1D( "beamEnergy_reco","Reconstructed Beam Energy",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    beamEnergy_reco->GetXaxis()->SetTitle("Reconstructed Beam Energy [GeV]");
    beamEnergy_reco->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.beamE.get_width()));
    
    beamEnergy_error = new TH1D( "beamEnergy_error","Error on Beam Energy",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    beamEnergy_error->GetXaxis()->SetTitle("(True-Reco) / True");
    beamEnergy_error->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.error.get_width()));
    
    beamEnergy_reco_mc = new TH2D( "beamEnergy_reco_mc","True vs Reconstructed Beam Energy",
    binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max(),
    binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max());
    beamEnergy_reco_mc->GetXaxis()->SetTitle("Reconstructed Beam Energy [GeV]");
    beamEnergy_reco_mc->GetYaxis()->SetTitle("True Beam Energy [GeV]");
    
    beamEnergyCal_mc = new TH1D( "beamEnergyCal_mc","True Beam Energy",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    beamEnergyCal_mc->GetXaxis()->SetTitle("True Beam Energy [GeV]");
    beamEnergyCal_mc->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.beamE.get_width()));
    
    beamEnergyCal_reco = new TH1D( "beamEnergyCal_reco","Reconstructed Calorimetric Beam Energy",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    beamEnergyCal_reco->GetXaxis()->SetTitle("Reconstructed Calorimetric  Beam Energy [GeV]");
    beamEnergyCal_reco->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.beamE.get_width()));
    
    beamEnergyCal_error = new TH1D( "beamEnergyCal_error","Error on Calorimetric Beam Energy",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    beamEnergyCal_error->GetXaxis()->SetTitle("(True-Reco) / True");
    beamEnergyCal_error->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.error.get_width()));
    
    beamEnergyCal_reco_mc = new TH2D( "beamEnergyCal_reco_mc","True vs Reconstructed Calorimetric Beam Energy",
                                   binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max(),
                                   binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max());
    beamEnergyCal_reco_mc->GetXaxis()->SetTitle("Reconstructed Calorimetric Beam Energy [GeV]");
    beamEnergyCal_reco_mc->GetYaxis()->SetTitle("True Beam Energy [GeV]");
    
    beamEnergy_beamEnergyCal = new TH2D( "beamEnergy_beamEnergyCal","Calorimetric Beam Energy vs Beam Energy from Tp",
                                      binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max(),
                                      binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max());
    beamEnergy_beamEnergyCal->GetXaxis()->SetTitle("Reconstructed Beam Energy from Tp [GeV]");
    beamEnergy_beamEnergyCal->GetYaxis()->SetTitle("Reconstructed Calorimetric Beam Energy [GeV]");
    
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
    
    nProngs_hist = new TH1D( "nProngs_hist","Number of Prongs",binList.multiplicity.get_nBins(), binList.multiplicity.get_min(), binList.multiplicity.get_max() );
    nProngs_hist->GetXaxis()->SetTitle("Number of Prongs");
    nProngs_hist->GetYaxis()->SetTitle("Candidates"); 

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
    
    //! ------------------------------------------------------------------------
    //! Cut Histograms
    //! ------------------------------------------------------------------------
    hCut_vertexCount = new TH1D( "hCut_vertexCount","Number of Vertices",binList.objectCount.get_nBins(), binList.objectCount.get_min(), binList.objectCount.get_max() );
    hCut_vertexCount->GetXaxis()->SetTitle("Number of Vertices");
    hCut_vertexCount->GetYaxis()->SetTitle("N(Events)");
    
    hCut_nProngs = new TH1D( "hCut_nProngs","Number of Prongs before CUTS",binList.multiplicity.get_nBins(), binList.multiplicity.get_min(), binList.multiplicity.get_max() );
    hCut_nProngs->GetXaxis()->SetTitle("Number of Prongs");
    hCut_nProngs->GetYaxis()->SetTitle("N(Events)");
    
    hCut_1Prong_Michel = new TH1D( "hCut_1Prong_Michel","Event Has Michel?",binList.michelID.get_nBins(), binList.michelID.get_min(), binList.michelID.get_max() );
    hCut_1Prong_Michel->GetXaxis()->SetTitle("0 = No Michel, 1 = Michel");
    hCut_1Prong_Michel->GetYaxis()->SetTitle("N(Events)");

    hCut_1Prong_eVis_nuclearTarget = new TH1D( "hCut_1Prong_eVis_nuclearTarget","Visible Energy in Nuclear Target",binList.eVis_nuclearTarget.get_nBins(), binList.eVis_nuclearTarget.get_min(), binList.eVis_nuclearTarget.get_max() );
    hCut_1Prong_eVis_nuclearTarget->GetXaxis()->SetTitle("Visible Energy in Nuclear Target [MeV]");
    hCut_1Prong_eVis_nuclearTarget->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.eVis_nuclearTarget.get_width()));
    
    hCut_1Prong_eVis_other = new TH1D( "hCut_1Prong_eVis_other","Visible Energy in Tracker + ECAL + HCAL",binList.eVis_other.get_nBins(), binList.eVis_other.get_min(), binList.eVis_other.get_max() );
    hCut_1Prong_eVis_other->GetXaxis()->SetTitle("Visible Energy in Tracker + ECAL + HCAL [MeV]");
    hCut_1Prong_eVis_other->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.eVis_other.get_width()));
    
    hCut_1Prong_gamma1ConvDist = new TH1D( "hCut_1Prong_gamma1ConvDist","Leading Photon Conversion Distance",binList.bin_photonConvLength.get_nBins(), binList.bin_photonConvLength.get_min(), binList.bin_photonConvLength.get_max() );
    hCut_1Prong_gamma1ConvDist->GetXaxis()->SetTitle("Leading Photon Conversion Distance");
    hCut_1Prong_gamma1ConvDist->GetYaxis()->SetTitle(Form("Candidates / %3.2f [MeV]",binList.bin_photonConvLength.get_width()));
    
    hCut_1Prong_gamma2ConvDist = new TH1D( "hCut_1Prong_gamma2ConvDist","Second Photon Conversion Distance",binList.bin_photonConvLength.get_nBins(), binList.bin_photonConvLength.get_min(), binList.bin_photonConvLength.get_max() );
    hCut_1Prong_gamma2ConvDist->GetXaxis()->SetTitle("Second Photon Conversion Distance");
    hCut_1Prong_gamma2ConvDist->GetYaxis()->SetTitle(Form("Candidates / %3.2f [MeV]",binList.bin_photonConvLength.get_width()));
    
    hCut_1Prong_pi0invMass = new TH1D( "hCut_1Prong_pi0invMass","Reconstructed Pi0 Invariant Mass",binList.pi0_invMass.get_nBins(), binList.pi0_invMass.get_min(), binList.pi0_invMass.get_max() );
    hCut_1Prong_pi0invMass->GetXaxis()->SetTitle("Reconstructed Pi0 Invariant Mass [MeV]");
    hCut_1Prong_pi0invMass->GetYaxis()->SetTitle(Form("Candidates / %3.2f [MeV]",binList.pi0_invMass.get_width()));
    
    hCut_1Prong_neutrinoE = new TH1D( "hCut_1Prong_neutrinoE","Reconstructed Beam Energy",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    hCut_1Prong_neutrinoE->GetXaxis()->SetTitle("Reconstructed Beam Energy [GeV]");
    hCut_1Prong_neutrinoE->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.beamE.get_width()));
    
    hCut_1Prong_UnusedE = new TH1D( "hCut_1Prong_UnusedE","Unused Cluster Energy after Pi0 Reconstruction",binList.UnusedE.get_nBins(), binList.UnusedE.get_min(), binList.UnusedE.get_max() );
    hCut_1Prong_UnusedE->GetXaxis()->SetTitle("Unused Cluster Energy after Pi0 Reconstruction [MeV]");
    hCut_1Prong_UnusedE->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.UnusedE.get_width()));
    
    hCut_2Prong_Michel = new TH1D( "hCut_2Prong_Michel","Event Has Michel?",binList.michelID.get_nBins(), binList.michelID.get_min(), binList.michelID.get_max() );
    hCut_2Prong_Michel->GetXaxis()->SetTitle("0 = No Michel, 1 = Michel");
    hCut_2Prong_Michel->GetYaxis()->SetTitle("N(Events)");

    hCut_2Prong_eVis_nuclearTarget = new TH1D( "hCut_2Prong_eVis_nuclearTarget","Visible Energy in Nuclear Target",binList.eVis_nuclearTarget.get_nBins(), binList.eVis_nuclearTarget.get_min(), binList.eVis_nuclearTarget.get_max() );
    hCut_2Prong_eVis_nuclearTarget->GetXaxis()->SetTitle("Visible Energy in Nuclear Target [MeV]");
    hCut_2Prong_eVis_nuclearTarget->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.eVis_nuclearTarget.get_width()));
    
    hCut_2Prong_eVis_other = new TH1D( "hCut_2Prong_eVis_other","Visible Energy in Tracker + ECAL + HCAL",binList.eVis_other.get_nBins(), binList.eVis_other.get_min(), binList.eVis_other.get_max() );
    hCut_2Prong_eVis_other->GetXaxis()->SetTitle("Visible Energy in Tracker + ECAL + HCAL [MeV]");
    hCut_2Prong_eVis_other->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.eVis_other.get_width()));
    
    hCut_2Prong_gamma1ConvDist = new TH1D( "hCut_2Prong_gamma1ConvDist","Leading Photon Conversion Distance",binList.bin_photonConvLength.get_nBins(), binList.bin_photonConvLength.get_min(), binList.bin_photonConvLength.get_max() );
    hCut_2Prong_gamma1ConvDist->GetXaxis()->SetTitle("Leading Photon Conversion Distance");
    hCut_2Prong_gamma1ConvDist->GetYaxis()->SetTitle(Form("Candidates / %3.2f [MeV]",binList.bin_photonConvLength.get_width()));
    
    hCut_2Prong_gamma2ConvDist = new TH1D( "hCut_2Prong_gamma2ConvDist","Second Photon Conversion Distance",binList.bin_photonConvLength.get_nBins(), binList.bin_photonConvLength.get_min(), binList.bin_photonConvLength.get_max() );
    hCut_2Prong_gamma2ConvDist->GetXaxis()->SetTitle("Second Photon Conversion Distance");
    hCut_2Prong_gamma2ConvDist->GetYaxis()->SetTitle(Form("Candidates / %3.2f [MeV]",binList.bin_photonConvLength.get_width()));
    
    hCut_2Prong_pi0invMass = new TH1D( "hCut_2Prong_pi0invMass","Reconstructed Pi0 Invariant Mass",binList.pi0_invMass.get_nBins(), binList.pi0_invMass.get_min(), binList.pi0_invMass.get_max() );
    hCut_2Prong_pi0invMass->GetXaxis()->SetTitle("Reconstructed Pi0 Invariant Mass [MeV]");
    hCut_2Prong_pi0invMass->GetYaxis()->SetTitle(Form("Candidates / %3.2f [MeV]",binList.pi0_invMass.get_width()));
    
    hCut_2Prong_neutrinoE = new TH1D( "hCut_2Prong_neutrinoE","Reconstructed Beam Energy",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    hCut_2Prong_neutrinoE->GetXaxis()->SetTitle("Reconstructed Beam Energy [GeV]");
    hCut_2Prong_neutrinoE->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.beamE.get_width()));
    
    hCut_2Prong_UnusedE = new TH1D( "hCut_2Prong_UnusedE","Unused Cluster Energy after Pi0 Reconstruction",binList.UnusedE.get_nBins(), binList.UnusedE.get_min(), binList.UnusedE.get_max() );
    hCut_2Prong_UnusedE->GetXaxis()->SetTitle("Unused Cluster Energy after Pi0 Reconstruction [MeV]");
    hCut_2Prong_UnusedE->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.UnusedE.get_width()));
    
    hCut_pIDDiff = new TH1D( "hCut_pIDDiff","Proton Score - Pion Score",binList.particleScoreDiff.get_nBins(), binList.particleScoreDiff.get_min(), binList.particleScoreDiff.get_max() );
    hCut_pIDDiff->GetXaxis()->SetTitle("Proton Score - Pion Score");
    hCut_pIDDiff->GetYaxis()->SetTitle("N(Events)");
    
    hCut_protonScore_LLR = new TH1D( "hCut_protonScore_LLR","proton_protonScore_LLR",binList.particleScore_LLR.get_nBins(), binList.particleScore_LLR.get_min(), binList.particleScore_LLR.get_max() );
    hCut_protonScore_LLR->GetXaxis()->SetTitle("proton_protonScore_LLR");
    hCut_protonScore_LLR->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScore_LLR.get_width()));
    
    hCut_deltaInvMass = new TH1D( "hCut_deltaInvMass","deltaInvMass",binList.deltaInvMass.get_nBins(), binList.deltaInvMass.get_min(), binList.deltaInvMass.get_max() );
    hCut_deltaInvMass->GetXaxis()->SetTitle("hCut_deltaInvMass");
    hCut_deltaInvMass->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.deltaInvMass.get_width()));
    
}

void CCProtonPi0_Interaction::write_RootFile()
{
    cout<<">> Writing "<<rootDir<<endl;
    f->Write();
}



#endif
