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
        bin_theta_Diff.setBin(40,-5.0,5.0);
        bin_KE.setBin(100,0.0,10.0);
        
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

        temp = new MnvH1D( Form("%s_%d","P",i),"Reconstructed Muon Momentum",binList.size_muon_P, binList.a_muon_P);
        temp->GetXaxis()->SetTitle("Reconstructed P_{#mu} [GeV]");
        temp->GetYaxis()->SetTitle("Events/Bin");
        P.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","KE",i),"Reconstructed Muon Kinetic Energy",bin_KE.get_nBins(), bin_KE.get_min(), bin_KE.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed T_{#mu} [GeV]");
        temp->GetYaxis()->SetTitle(Form("Muons / %3.1f [GeV]",bin_KE.get_width()));
        KE.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","theta",i),"Reconstructed Muon Theta",binList.size_muon_theta, binList.a_muon_theta);
        temp->GetXaxis()->SetTitle("Reconstructed #theta_{#mu} [Degree]");
        temp->GetYaxis()->SetTitle("N(events)");
        theta.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","cos_theta",i),"Reconstructed Muon cos(#theta)",binList.muon_cos_theta.get_nBins(), binList.muon_cos_theta.get_min(), binList.muon_cos_theta.get_max());
        temp->GetXaxis()->SetTitle("Reconstructed cos(#theta_{#mu})");
        temp->GetYaxis()->SetTitle("N(events)");
        cos_theta.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","phi",i),"Reconstructed Muon Phi",binList.angle.get_nBins(), binList.angle.get_min(), binList.angle.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed #phi_{#mu}[Degree]");
        temp->GetYaxis()->SetTitle(Form("Muons / %3.1f [Degree]",binList.angle.get_width()));
        phi.push_back(temp);
    }

    // Error Histograms
    P_error = new TH1D( "P_error","Error on Muon Momentum",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    P_error->GetXaxis()->SetTitle("(P_{Reco}-P_{True})/P_{True}");
    P_error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));
 
    E_error = new TH1D( "E_error","Error on Muon Energy",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    E_error->GetXaxis()->SetTitle("(E_{Reco}-E_{True})/E_{True}");
    E_error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));
 
    KE_error = new TH1D( "KE_error","Error on Muon Kinetic Energy",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    KE_error->GetXaxis()->SetTitle("(T_{Reco}-T_{True})/T_{True}");
    KE_error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));
 
    theta_error = new TH1D( "theta_error","Error on Muon theta",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    theta_error->GetXaxis()->SetTitle("(#theta_{Reco}-#theta_{True})/#theta_{True}");
    theta_error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));
    
    cos_theta_error = new TH1D( "cos_theta_error","Error on Muon cos(#theta)",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    cos_theta_error->GetXaxis()->SetTitle("(cos(#theta_{Reco})-cos(#theta_{True}))/cos(#theta_{True})");
    cos_theta_error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));

    theta_Diff = new TH1D( "theta_Diff","Difference on Muon Theta",bin_theta_Diff.get_nBins(), bin_theta_Diff.get_min(), bin_theta_Diff.get_max() );
    theta_Diff->GetXaxis()->SetTitle("#theta_{Reco}-#theta_{True} [degree]");
    theta_Diff->GetYaxis()->SetTitle(Form("Events / %3.2f ",bin_theta_Diff.get_width()));

    reco_P_true_P = new TH2D( "reco_P_true_P","True vs Reconstructed Muon Momentum",binList.size_muon_P, binList.a_muon_P,binList.size_muon_P, binList.a_muon_P);
    reco_P_true_P->GetXaxis()->SetTitle("Reconstructed P_{#mu} [GeV]");
    reco_P_true_P->GetYaxis()->SetTitle("True P_{#mu} [GeV]");

    reco_E_true_E = new TH2D( "reco_E_true_E","True vs Reconstructed Muon Energy",bin_E.get_nBins(), bin_E.get_min(), bin_E.get_max(), bin_E.get_nBins(), bin_E.get_min(), bin_E.get_max());
    reco_E_true_E->GetXaxis()->SetTitle("Reconstructed E_{#mu} [GeV]");
    reco_E_true_E->GetYaxis()->SetTitle("True E_{#mu} [GeV]");

    E_Diff = new TH1D( "E_Diff","Difference on Muon Energy",bin_E_Diff.get_nBins(), bin_E_Diff.get_min(), bin_E_Diff.get_max() );
    E_Diff->GetXaxis()->SetTitle("E_{Reco}-E_{True} [GeV]");
    E_Diff->GetYaxis()->SetTitle(Form("Events / %3.2f ",bin_E_Diff.get_width()));

    // Cross Section Variables 
    muon_P_all = new MnvH1D( "muon_P_all","Data All P_{#mu}",binList.size_muon_P, binList.a_muon_P);
    muon_P_all->GetXaxis()->SetTitle("Muon Momentum (GeV/c)");
    muon_P_all->GetYaxis()->SetTitle("Events/Bin");

    muon_P_mc_truth_signal = new MnvH1D( "muon_P_mc_truth_signal","MC Truth Signal P_{#mu}",binList.size_muon_P, binList.a_muon_P);
    muon_P_mc_truth_signal->GetXaxis()->SetTitle("Muon Momentum (GeV/c)");
    muon_P_mc_truth_signal->GetYaxis()->SetTitle("Events/Bin");
    AddVertErrorBands_MC(muon_P_mc_truth_signal);
    AddLatErrorBands_MC(muon_P_mc_truth_signal);

    muon_P_mc_reco_all = new MnvH1D( "muon_P_mc_reco_all","MC All Reconstructed P_{#mu}",binList.size_muon_P, binList.a_muon_P);
    muon_P_mc_reco_all->GetXaxis()->SetTitle("Muon Momentum (GeV/c)");
    muon_P_mc_reco_all->GetYaxis()->SetTitle("Events/Bin");
    AddVertErrorBands_MC(muon_P_mc_reco_all);
    AddLatErrorBands_MC(muon_P_mc_reco_all);

    muon_P_mc_reco_signal = new MnvH1D( "muon_P_mc_reco_signal","MC Reconstructed Signal P_{#mu}",binList.size_muon_P, binList.a_muon_P);
    muon_P_mc_reco_signal->GetXaxis()->SetTitle("Muon Momentum (GeV/c)");
    muon_P_mc_reco_signal->GetYaxis()->SetTitle("Events/Bin");
    AddVertErrorBands_MC(muon_P_mc_reco_signal);
    AddLatErrorBands_MC(muon_P_mc_reco_signal);

    muon_P_mc_reco_bckg = new MnvH1D( "muon_P_mc_reco_bckg","MC Reconstructed Background P_{#mu}",binList.size_muon_P, binList.a_muon_P);
    muon_P_mc_reco_bckg->GetXaxis()->SetTitle("Muon Momentum (GeV/c)");
    muon_P_mc_reco_bckg->GetYaxis()->SetTitle("Events/Bin");
    AddVertErrorBands_MC(muon_P_mc_reco_bckg);
    AddLatErrorBands_MC(muon_P_mc_reco_bckg);

    muon_P_response = new MnvH2D( "muon_P_response","Signal Muon Momentum",binList.size_muon_P, binList.a_muon_P,binList.size_muon_P, binList.a_muon_P);
    muon_P_response->GetXaxis()->SetTitle("Reconstructed Muon Momentum (GeV/c)");
    muon_P_response->GetYaxis()->SetTitle("True Muon Momentum (GeV/c)");
    AddVertErrorBands_MC(muon_P_response);
    AddLatErrorBands_MC(muon_P_response);

    // Muon Theta
    muon_theta_all = new MnvH1D( "muon_theta_all","Data All #theta_{#mu}",binList.size_muon_theta, binList.a_muon_theta);
    muon_theta_all->GetXaxis()->SetTitle("Muon Angle (deg)");
    muon_theta_all->GetYaxis()->SetTitle("Events/Bin");

    muon_theta_mc_truth_signal = new MnvH1D( "muon_theta_mc_truth_signal","MC Truth Signal #theta_{#mu}",binList.size_muon_theta, binList.a_muon_theta);
    muon_theta_mc_truth_signal->GetXaxis()->SetTitle("Muon Angle (deg)");
    muon_theta_mc_truth_signal->GetYaxis()->SetTitle("Events/Bin");
    AddVertErrorBands_MC(muon_theta_mc_truth_signal);
    AddLatErrorBands_MC(muon_theta_mc_truth_signal);

    muon_theta_mc_reco_all = new MnvH1D( "muon_theta_mc_reco_all","MC All Reconstructed #theta_{#mu}",binList.size_muon_theta, binList.a_muon_theta);
    muon_theta_mc_reco_all->GetXaxis()->SetTitle("Muon Angle (deg)");
    muon_theta_mc_reco_all->GetYaxis()->SetTitle("Events/Bin");
    AddVertErrorBands_MC(muon_theta_mc_reco_all);
    AddLatErrorBands_MC(muon_theta_mc_reco_all);

    muon_theta_mc_reco_signal = new MnvH1D( "muon_theta_mc_reco_signal","MC Reconstructed Signal #theta_{#mu}",binList.size_muon_theta, binList.a_muon_theta);
    muon_theta_mc_reco_signal->GetXaxis()->SetTitle("Muon Angle (deg)");
    muon_theta_mc_reco_signal->GetYaxis()->SetTitle("Events/Bin");
    AddVertErrorBands_MC(muon_theta_mc_reco_signal);
    AddLatErrorBands_MC(muon_theta_mc_reco_signal);

    muon_theta_mc_reco_bckg = new MnvH1D( "muon_theta_mc_reco_bckg","MC Reconstructed Background #theta_{#mu}",binList.size_muon_theta, binList.a_muon_theta);
    muon_theta_mc_reco_bckg->GetXaxis()->SetTitle("Muon Angle (deg)");
    muon_theta_mc_reco_bckg->GetYaxis()->SetTitle("Events/Bin");
    AddVertErrorBands_MC(muon_theta_mc_reco_bckg);
    AddLatErrorBands_MC(muon_theta_mc_reco_bckg);

    muon_theta_response = new MnvH2D( "muon_theta_response","Signal Muon Angle",binList.size_muon_theta, binList.a_muon_theta,binList.size_muon_theta, binList.a_muon_theta);
    muon_theta_response->GetXaxis()->SetTitle("Reconstructed Muon Angle (deg)");
    muon_theta_response->GetYaxis()->SetTitle("True Muon Angle (deg)");
    AddVertErrorBands_MC(muon_theta_response);
    AddLatErrorBands_MC(muon_theta_response);

    muon_P_shift = new MnvH1D( "muon_P_shift","Muon Momentum Shuft",50,-100,500);
    muon_P_shift->GetXaxis()->SetTitle("Muon Momentum Shift[MeV]");
    muon_P_shift->GetYaxis()->SetTitle("Events/Bin");
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
        cos_theta[i]->Write();
        phi[i]->Write();
    }

    P_error->Write();
    E_error->Write();
    KE_error->Write();
    theta_error->Write();
    cos_theta_error->Write();
    
    reco_P_true_P->Write();
    reco_E_true_E->Write();

    E_Diff->Write();
    theta_Diff->Write();

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

    muon_P_shift->Write();

    f->Close();
}


#endif



