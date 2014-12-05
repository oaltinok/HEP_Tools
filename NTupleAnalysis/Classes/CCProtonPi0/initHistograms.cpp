#ifndef initHistograms_cpp
#define initHistograms_cpp

#include "CCProtonPi0.h"
using namespace std;

void CCProtonPi0::initHistograms()
{
    cout<<"\tInitializing CCProtonPi0 Histograms"<<endl;
    
    //--------------------------------------------------------------------------
    // PID Histograms
    
    pID_purity = new TH1D( "pID_purity","Proton Purity",binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max() );
    pID_purity->GetXaxis()->SetTitle("Proton Purity = Captured Proton / Captured Total Events");
    pID_purity->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScore.get_width()));
    
    pID_efficiency = new TH1D( "pID_efficiency","Proton Efficiency",binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max() );
    pID_efficiency->GetXaxis()->SetTitle("Proton Efficiency = Captured Proton / Total Protons");
    pID_efficiency->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScore.get_width()));
    
    // Proton Score 
    pID_piplus_protonScore = new TH1D( "pID_piplus_protonScore","piplus_protonScore",binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max() );
    pID_piplus_protonScore->GetXaxis()->SetTitle("piplus_protonScore");
    pID_piplus_protonScore->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScore.get_width()));
    
    pID_piminus_protonScore = new TH1D( "pID_piminus_protonScore","piminus_protonScore",binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max() );
    pID_piminus_protonScore->GetXaxis()->SetTitle("piminus_protonScore");
    pID_piminus_protonScore->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScore.get_width()));
    
    pID_proton_protonScore = new TH1D( "pID_proton_protonScore","proton_protonScore",binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max() );
    pID_proton_protonScore->GetXaxis()->SetTitle("proton_protonScore");
    pID_proton_protonScore->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScore.get_width()));
    
    pID_other_protonScore = new TH1D( "pID_other_protonScore","other_protonScore",binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max() );
    pID_other_protonScore->GetXaxis()->SetTitle("other_protonScore");
    pID_other_protonScore->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScore.get_width()));
    
    // Pion Score
    pID_piplus_pionScore = new TH1D( "pID_piplus_pionScore","piplus_pionScore",binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max() );
    pID_piplus_pionScore->GetXaxis()->SetTitle("piplus_pionScore");
    pID_piplus_pionScore->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScore.get_width()));
    
    pID_piminus_pionScore = new TH1D( "pID_piminus_pionScore","piminus_pionScore",binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max() );
    pID_piminus_pionScore->GetXaxis()->SetTitle("piminus_pionScore");
    pID_piminus_pionScore->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScore.get_width()));
    
    pID_proton_pionScore = new TH1D( "pID_proton_pionScore","proton_pionScore",binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max() );
    pID_proton_pionScore->GetXaxis()->SetTitle("proton_pionScore");
    pID_proton_pionScore->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScore.get_width()));
    
    pID_other_pionScore = new TH1D( "pID_other_pionScore","other_pionScore",binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max() );
    pID_other_pionScore->GetXaxis()->SetTitle("other_pionScore");
    pID_other_pionScore->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScore.get_width()));
    
    // 2D Particle Score
    pID_proton_pionScore_protonScore = new TH2D( "pID_proton_pionScore_protonScore","pID_proton_pionScore_protonScore",
                                                 binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max(),
                                                 binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max());
    pID_proton_pionScore_protonScore->GetXaxis()->SetTitle("Pion Score");
    pID_proton_pionScore_protonScore->GetYaxis()->SetTitle("Proton Score");
    
    pID_piplus_pionScore_protonScore = new TH2D( "pID_piplus_pionScore_protonScore","pID_piplus_pionScore_protonScore",
                                                 binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max(),
                                                 binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max());
    pID_piplus_pionScore_protonScore->GetXaxis()->SetTitle("Pion Score");
    pID_piplus_pionScore_protonScore->GetYaxis()->SetTitle("Proton Score");
                                                 
    pID_piminus_pionScore_protonScore = new TH2D( "pID_piminus_pionScore_protonScore","pID_piminus_pionScore_protonScore",
                                                 binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max(),
                                                 binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max());
    pID_piminus_pionScore_protonScore->GetXaxis()->SetTitle("Pion Score");
    pID_piminus_pionScore_protonScore->GetYaxis()->SetTitle("Proton Score");
    
    //--------------------------------------------------------------------------
    
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
    
    w2_mc = new TH1D( "w2_mc","True W",binList.w2.get_nBins(), binList.w2.get_min(), binList.w2.get_max() );
    w2_mc->GetXaxis()->SetTitle("True W [GeV]");
    w2_mc->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.w2.get_width()));
    
    w2_reco = new TH1D( "w2_reco","Reconstructed W",binList.w2.get_nBins(), binList.w2.get_min(), binList.w2.get_max() );
    w2_reco->GetXaxis()->SetTitle("Reconstructed W [GeV]");
    w2_reco->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.w2.get_width()));
    
    w2_error = new TH1D( "w2_error","Error on W",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    w2_error->GetXaxis()->SetTitle("(True-Reco) / True");
    w2_error->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.error.get_width()));
    
    w2_reco_mc = new TH2D( "w2_reco_mc","True vs Reconstructed W",
                           binList.w2.get_nBins(), binList.w2.get_min(), binList.w2.get_max(),
                           binList.w2.get_nBins(), binList.w2.get_min(), binList.w2.get_max());
                           w2_reco_mc->GetXaxis()->SetTitle("Reconstructed W [GeV]");
                           w2_reco_mc->GetYaxis()->SetTitle("True W [GeV]");
                           
    mc_w_DIS = new TH1D( "mc_w_DIS","True W for DIS",binList.w2.get_nBins(), binList.w2.get_min(), binList.w2.get_max() );
    mc_w_DIS->GetXaxis()->SetTitle("True W for DIS [GeV]");
    mc_w_DIS->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.w2.get_width()));
    
    mc_w_RES = new TH1D( "mc_w_RES","True W for RES",binList.w2.get_nBins(), binList.w2.get_min(), binList.w2.get_max() );
    mc_w_RES->GetXaxis()->SetTitle("True W for RES [GeV]");
    mc_w_RES->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.w2.get_width()));
    
    mc_w_CCQE = new TH1D( "mc_w_CCQE","True W for CCQE",binList.w2.get_nBins(), binList.w2.get_min(), binList.w2.get_max() );
    mc_w_CCQE->GetXaxis()->SetTitle("True W for CCQE [GeV]");
    mc_w_CCQE->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.w2.get_width()));
    
    int_channel = new TH1D( "int_channel","Interaction Channel",binList.int_channel.get_nBins(), binList.int_channel.get_min(), binList.int_channel.get_max() );
    int_channel->GetXaxis()->SetTitle("1 = QE, 2 = Resonant, 3 = DIS, 4 = Coh pi");
    int_channel->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.int_channel.get_width()));
    
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
    
    nProngs_hist_initial = new TH1D( "nProngs_hist_initial","Number of Prongs before CUTS",binList.multiplicity.get_nBins(), binList.multiplicity.get_min(), binList.multiplicity.get_max() );
    nProngs_hist_initial->GetXaxis()->SetTitle("Number of Prongs");
    nProngs_hist_initial->GetYaxis()->SetTitle("Candidates"); 
    
    // Debugging Histograms
    w2fail_w2 = new TH1D( "w2fail_w2","mN^{2}+2mN(Epi + Tp) - Q^{2}",binList.w2fail.get_nBins(), binList.w2fail.get_min(), binList.w2fail.get_max() );
    w2fail_w2->GetXaxis()->SetTitle("mN^{2}+2mN(Epi + Tp) - Q^{2} [GeV^{2}]");
    w2fail_w2->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.w2fail.get_width()));
    
    w2fail_q2 = new TH1D( "w2fail_q2","Reconstructed Q^{2}",binList.q2.get_nBins(), binList.q2.get_min(), binList.q2.get_max() );
    w2fail_q2->GetXaxis()->SetTitle("Reconstructed Q^{2} [GeV^{2}]");
    w2fail_q2->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.q2.get_width()));
    
    w2fail_term1 = new TH1D( "w2fail_term1","mN^{2}+2mN(Epi + Tp)",binList.q2.get_nBins(), binList.q2.get_min(), binList.q2.get_max() );
    w2fail_term1->GetXaxis()->SetTitle("mN^{2}+2mN(Epi + Tp) [GeV^{2}]");
    w2fail_term1->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.q2.get_width()));
    
    w2fail_term1_q2 = new TH2D( "w2fail_term1_q2","Q^{2} vs (mN^{2}+2mN(Epi + Tp))",  binList.q2.get_nBins(), binList.q2.get_min(), binList.q2.get_max(),
                                                                        binList.q2.get_nBins(), binList.q2.get_min(), binList.q2.get_max() );
    w2fail_term1_q2->GetXaxis()->SetTitle("mN^{2}+2mN(Epi + Tp) [GeV^{2}]");
    w2fail_term1_q2->GetYaxis()->SetTitle("Reconstructed Q^{2} [GeV^{2}]");
    
    w2fail_Enu = new TH1D( "w2fail_Enu","Reconstructed Beam Energy",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    w2fail_Enu->GetXaxis()->SetTitle("Reconstructed Beam Energy [GeV]");
    w2fail_Enu->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.beamE.get_width()));
    
    w2fail_muon_E = new TH1D( "w2fail_muon_E","Reconstructed Muon Energy",binList.muonE.get_nBins(), binList.muonE.get_min(), binList.muonE.get_max() );
    w2fail_muon_E->GetXaxis()->SetTitle("Reconstructed Muon Energy [GeV]");
    w2fail_muon_E->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.muonE.get_width()));
    
    w2fail_muon_Pz = new TH1D( "w2fail_muon_Pz","Reconstructed Muon Momentum(z)",binList.muonE.get_nBins(), binList.muonE.get_min(), binList.muonE.get_max() );
    w2fail_muon_Pz->GetXaxis()->SetTitle("Reconstructed Muon Momentum(z) [GeV]");
    w2fail_muon_Pz->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.muonE.get_width()));
    
    w2fail_pion_E = new TH1D( "w2fail_pion_E","Reconstructed Pion Energy",binList.pionE.get_nBins(), binList.pionE.get_min(), binList.pionE.get_max() );
    w2fail_pion_E->GetXaxis()->SetTitle("Reconstructed Pion Energy [GeV]");
    w2fail_pion_E->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.pionE.get_width()));
    
    w2fail_proton_KE = new TH1D( "w2fail_proton_KE","Reconstructed Proton Kinetic Energy",binList.protonKE.get_nBins(), binList.protonKE.get_min(), binList.protonKE.get_max() );
    w2fail_proton_KE->GetXaxis()->SetTitle("Reconstructed Proton Kinetic Energy [GeV]");
    w2fail_proton_KE->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.protonKE.get_width()));
    
    
}

#endif
