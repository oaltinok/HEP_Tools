/*
    See CCProtonPi0_Muon.h header for Class Information
*/
#ifndef CCProtonPi0_Muon_cpp
#define CCProtonPi0_Muon_cpp

#include "CCProtonPi0_Muon.h"

using namespace PlotUtils;

CCProtonPi0_Muon::CCProtonPi0_Muon(bool isModeReduce, bool isMC) : CCProtonPi0_Particle()
{
    std::cout<<"Initializing CCProtonPi0_Muon"<<std::endl;
    
    if(isModeReduce){
        std::cout<<"\tNTuple Reduce Mode -- Will not create ROOT Files"<<std::endl;
    }else{
        // File Locations
        if (isMC) rootDir = Folder_List::rootDir_Muon_mc;
        else rootDir = Folder_List::rootDir_Muon_data;

        std::cout<<"\tRoot File: "<<rootDir<<std::endl;
        
        // Create Root File 
        f = new TFile(rootDir.c_str(),"RECREATE");

        // Initialize Bins
        bin_E.setBin(10,0.0,10.0);
        bin_E_Diff.setBin(100,-1.0,1.0);
        bin_P.setBin(10,0.0,10.0);
        bin_KE.setBin(100,0.0,10.0);
        bin_muonTheta.setBin(12,0.0,25.0);
        
        initHistograms();
    }
    
    std::cout<<"Done!"<<std::endl;
}


void CCProtonPi0_Muon::initHistograms()
{
    MnvH1D* temp = NULL;

    for (int i = 0; i < nHistograms; i++){
        temp = new MnvH1D( Form("%s_%d","E",i),"Reconstructed Muon Energy",bin_E.get_nBins(), bin_E.get_min(), bin_E.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed E_{#mu} [GeV]");
        temp->GetYaxis()->SetTitle(Form("Muons / %3.1f [GeV] ",bin_E.get_width()));
        E.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","P",i),"Reconstructed Muon Momentum",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed P_{#mu} [GeV]");
        temp->GetYaxis()->SetTitle(Form("Muons / %3.1f [GeV] ",bin_P.get_width()));
        P.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","KE",i),"Reconstructed Muon Kinetic Energy",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed T_{#mu} [GeV]");
        temp->GetYaxis()->SetTitle(Form("Muons / %3.1f [GeV]",bin_P.get_width()));
        KE.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","theta",i),"Reconstructed Muon Theta",bin_muonTheta.get_nBins(), bin_muonTheta.get_min(), bin_muonTheta.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed #theta_{#mu} [Degree]");
        temp->GetYaxis()->SetTitle(Form("Muons / %3.1f [Degree] ",bin_muonTheta.get_width()));
        theta.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","phi",i),"Reconstructed Muon Phi",binList.angle.get_nBins(), binList.angle.get_min(), binList.angle.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed #phi_{#mu}[Degree]");
        temp->GetYaxis()->SetTitle(Form("Muons / %3.1f [Degree]",binList.angle.get_width()));
        phi.push_back(temp);
    }

        
    reco_P_true_P = new TH2D( "reco_P_true_P","True vs Reconstructed Muon Momentum",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max(), bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max());
    reco_P_true_P->GetXaxis()->SetTitle("Reconstructed P_{#mu} [GeV]");
    reco_P_true_P->GetYaxis()->SetTitle("True P_{#mu} [GeV]");

    P_error = new TH1D( "P_error","Error on Muon Momentum",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    P_error->GetXaxis()->SetTitle("(P_{Reco}-P_{True})/P_{True}");
    P_error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));
 
    reco_E_true_E = new TH2D( "reco_E_true_E","True vs Reconstructed Muon Energy",bin_E.get_nBins(), bin_E.get_min(), bin_E.get_max(), bin_E.get_nBins(), bin_E.get_min(), bin_E.get_max());
    reco_E_true_E->GetXaxis()->SetTitle("Reconstructed E_{#mu} [GeV]");
    reco_E_true_E->GetYaxis()->SetTitle("True E_{#mu} [GeV]");

    E_error = new TH1D( "E_error","Error on Muon Energy",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    E_error->GetXaxis()->SetTitle("(E_{Reco}-E_{True})/E_{True}");
    E_error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));

    E_Diff = new TH1D( "E_Diff","Difference on Muon Energy",bin_E_Diff.get_nBins(), bin_E_Diff.get_min(), bin_E_Diff.get_max() );
    E_Diff->GetXaxis()->SetTitle("E_{Reco}-E_{True} [GeV]");
    E_Diff->GetYaxis()->SetTitle(Form("Events / %3.2f ",bin_E_Diff.get_width()));

    // Cross Section Variables 
    muon_P_all = new MnvH1D( "muon_P_all","Data All P_{#mu}",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    muon_P_all->GetXaxis()->SetTitle("P_{#mu} [GeV]");
    muon_P_all->GetYaxis()->SetTitle("N(Events)");

    muon_P_mc_truth_signal = new MnvH1D( "muon_P_mc_truth_signal","MC Truth Signal P_{#mu}",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    muon_P_mc_truth_signal->GetXaxis()->SetTitle("P_{#mu} [GeV]");
    muon_P_mc_truth_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_MC(muon_P_mc_truth_signal);

    muon_P_mc_reco_all = new MnvH1D( "muon_P_mc_reco_all","MC All Reconstructed P_{#mu}",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    muon_P_mc_reco_all->GetXaxis()->SetTitle("P_{#mu} [GeV]");
    muon_P_mc_reco_all->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_MC(muon_P_mc_reco_all);

    muon_P_mc_reco_signal = new MnvH1D( "muon_P_mc_reco_signal","MC Reconstructed Signal P_{#mu}",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    muon_P_mc_reco_signal->GetXaxis()->SetTitle("P_{#mu} [GeV]");
    muon_P_mc_reco_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_MC(muon_P_mc_reco_signal);

    muon_P_mc_reco_bckg = new MnvH1D( "muon_P_mc_reco_bckg","MC Reconstructed Background P_{#mu}",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    muon_P_mc_reco_bckg->GetXaxis()->SetTitle("P_{#mu} [GeV]");
    muon_P_mc_reco_bckg->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_MC(muon_P_mc_reco_bckg);

    muon_P_response = new MnvH2D( "muon_P_response","Signal Muon Momentum",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max(),bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    muon_P_response->GetXaxis()->SetTitle("Reconstructed P_{#mu} [GeV]");
    muon_P_response->GetYaxis()->SetTitle("True P_{#mu} [GeV]");
    AddVertErrorBands_MC(muon_P_response);

    // Muon Theta
    muon_theta_all = new MnvH1D( "muon_theta_all","Data All #theta_{#mu}",bin_muonTheta.get_nBins(), bin_muonTheta.get_min(), bin_muonTheta.get_max() );
    muon_theta_all->GetXaxis()->SetTitle("#theta_{#mu} [degree]");
    muon_theta_all->GetYaxis()->SetTitle("N(Events)");

    muon_theta_mc_truth_signal = new MnvH1D( "muon_theta_mc_truth_signal","MC Truth Signal #theta_{#mu}",bin_muonTheta.get_nBins(), bin_muonTheta.get_min(), bin_muonTheta.get_max() );
    muon_theta_mc_truth_signal->GetXaxis()->SetTitle("#theta_{#mu} [degree]");
    muon_theta_mc_truth_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_MC(muon_theta_mc_truth_signal);

    muon_theta_mc_reco_all = new MnvH1D( "muon_theta_mc_reco_all","MC All Reconstructed #theta_{#mu}",bin_muonTheta.get_nBins(), bin_muonTheta.get_min(), bin_muonTheta.get_max() );
    muon_theta_mc_reco_all->GetXaxis()->SetTitle("#theta_{#mu} [degree]");
    muon_theta_mc_reco_all->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_MC(muon_theta_mc_reco_all);

    muon_theta_mc_reco_signal = new MnvH1D( "muon_theta_mc_reco_signal","MC Reconstructed Signal #theta_{#mu}",bin_muonTheta.get_nBins(), bin_muonTheta.get_min(), bin_muonTheta.get_max() );
    muon_theta_mc_reco_signal->GetXaxis()->SetTitle("#theta_{#mu} [degree]");
    muon_theta_mc_reco_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_MC(muon_theta_mc_reco_signal);

    muon_theta_mc_reco_bckg = new MnvH1D( "muon_theta_mc_reco_bckg","MC Reconstructed Background #theta_{#mu}",bin_muonTheta.get_nBins(), bin_muonTheta.get_min(), bin_muonTheta.get_max() );
    muon_theta_mc_reco_bckg->GetXaxis()->SetTitle("#theta_{#mu} [degree]");
    muon_theta_mc_reco_bckg->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_MC(muon_theta_mc_reco_bckg);

    muon_theta_response = new MnvH2D( "muon_theta_response","Signal Muon Momentum",bin_muonTheta.get_nBins(), bin_muonTheta.get_min(), bin_muonTheta.get_max(),bin_muonTheta.get_nBins(), bin_muonTheta.get_min(), bin_muonTheta.get_max() );
    muon_theta_response->GetXaxis()->SetTitle("Reconstructed #theta_{#mu} [degree]");
    muon_theta_response->GetYaxis()->SetTitle("True #theta_{#mu} [degree]");
    AddVertErrorBands_MC(muon_theta_response);
}


void CCProtonPi0_Muon::writeHistograms()
{
    std::cout<<">> Writing "<<rootDir<<std::endl;
    f->cd();

    for (int i = 0; i < nHistograms; i++){
        E[i]->Write();
        P[i]->Write();
        KE[i]->Write();
        theta[i]->Write();
        phi[i]->Write();
    }

    P_error->Write();
    reco_P_true_P->Write();

    E_error->Write();
    reco_E_true_E->Write();

    E_Diff->Write();

    muon_P_all->Write();
    muon_P_mc_truth_signal->Write();
    muon_P_mc_reco_all->Write();
    muon_P_mc_reco_signal->Write();
    muon_P_mc_reco_bckg->Write();
    muon_P_response->Write();
 
    muon_theta_all->Write();
    muon_theta_mc_truth_signal->Write();
    muon_theta_mc_reco_all->Write();
    muon_theta_mc_reco_signal->Write();
    muon_theta_mc_reco_bckg->Write();
    muon_theta_response->Write();
   
    f->Close();
}


#endif



