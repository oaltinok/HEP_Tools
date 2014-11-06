#ifndef initHistograms_cpp
#define initHistograms_cpp

#include "CCProtonPi0.h"

void CCProtonPi0::initHistograms()
{
    pID_purity = new TH1D( "pID_purity","Proton Purity",binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max() );
    pID_purity->GetXaxis()->SetTitle("Proton Purity = Captured Proton / Captured Total Events");
    pID_purity->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScore.get_width()));
    
    pID_efficiency = new TH1D( "pID_efficiency","Proton Efficiency",binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max() );
    pID_efficiency->GetXaxis()->SetTitle("Proton Efficiency = Captured Proton / Total Protons");
    pID_efficiency->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScore.get_width()));
    
    pID_piplus = new TH1D( "pID_piplus","Pi Plus",binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max() );
    pID_piplus->GetXaxis()->SetTitle("Pi Plus");
    pID_piplus->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScore.get_width()));
    
    pID_piminus = new TH1D( "pID_piminus","Pi Minus",binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max() );
    pID_piminus->GetXaxis()->SetTitle("Pi Minus");
    pID_piminus->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScore.get_width()));
    
    pID_proton = new TH1D( "pID_proton","Proton",binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max() );
    pID_proton->GetXaxis()->SetTitle("Proton");
    pID_proton->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScore.get_width()));
    
    pID_other = new TH1D( "pID_other","Other",binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max() );
    pID_other->GetXaxis()->SetTitle("Other");
    pID_other->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScore.get_width()));
    
    beamEnergy_mc = new TH1D( "beamEnergy_mc","True Beam Energy",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    beamEnergy_mc->GetXaxis()->SetTitle("True Beam Energy [MeV]");
    beamEnergy_mc->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.beamE.get_width()));
    
    beamEnergy_reco = new TH1D( "beamEnergy_reco","Reconstructed Beam Energy",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    beamEnergy_reco->GetXaxis()->SetTitle("Reconstructed Beam Energy [MeV]");
    beamEnergy_reco->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.beamE.get_width()));
    
    beamEnergy_error = new TH1D( "beamEnergy_error","Error on Beam Energy",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    beamEnergy_error->GetXaxis()->SetTitle("(Reco- True) / True");
    beamEnergy_error->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.error.get_width()));
    
    beamEnergy_reco_mc = new TH2D( "beamEnergy_reco_mc","True vs Reconstructed Beam Energy",
                                binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max(),
                                binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max());
    beamEnergy_reco_mc->GetXaxis()->SetTitle("Reconstructed Beam Energy [MeV]");
    beamEnergy_reco_mc->GetYaxis()->SetTitle("True Beam Energy [MeV]");
    
    q2_mc = new TH1D( "q2_mc","True Q^{2}",binList.q2.get_nBins(), binList.q2.get_min(), binList.q2.get_max() );
    q2_mc->GetXaxis()->SetTitle("True Q^{2} [MeV]");
    q2_mc->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.q2.get_width()));
    
    q2_reco = new TH1D( "q2_reco","Reconstructed Q^{2}",binList.q2.get_nBins(), binList.q2.get_min(), binList.q2.get_max() );
    q2_reco->GetXaxis()->SetTitle("Reconstructed Q^{2} [MeV]");
    q2_reco->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.q2.get_width()));
    
    q2_error = new TH1D( "q2_error","Error on Q^{2}",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    q2_error->GetXaxis()->SetTitle("(Reco- True) / True");
    q2_error->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.error.get_width()));
    
    q2_reco_mc = new TH2D( "q2_reco_mc","True vs Reconstructed Q^{2}",
                                binList.q2.get_nBins(), binList.q2.get_min(), binList.q2.get_max(),
                                binList.q2.get_nBins(), binList.q2.get_min(), binList.q2.get_max());
    q2_reco_mc->GetXaxis()->SetTitle("Reconstructed Q^{2} [MeV]");
    q2_reco_mc->GetYaxis()->SetTitle("True Q^{2} [MeV]");
    
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
    vertex_z_error->GetXaxis()->SetTitle("(Reco- True) / True");
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
    deltaInvMass_error->GetXaxis()->SetTitle("(Reco- True) / True");
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
}

#endif
