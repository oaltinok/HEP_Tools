/*
    See CCProtonPi0_Interaction.h header for Class Information
*/
#ifndef CCProtonPi0_Interaction_cpp
#define CCProtonPi0_Interaction_cpp

#include "CCProtonPi0_Interaction.h"

using namespace PlotUtils;

CCProtonPi0_Interaction::CCProtonPi0_Interaction(int nMode, bool isMC) : CCProtonPi0_NTupleAnalysis(nMode)
{
    std::cout<<"Initializing CCProtonPi0_Interaction"<<std::endl;
    
    if(nMode == 0){
        std::cout<<"\tNTuple Reduce Mode -- Will not create ROOT Files"<<std::endl;
    }else{
        if (isMC) rootDir = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed +  branchDir + "Interaction.root";
        else rootDir = Folder_List::rootOut + Folder_List::Data + Folder_List::analyzed + branchDir + "Interaction.root";
        
        std::cout<<"\tRoot File: "<<rootDir<<std::endl;
 
        // Create Root File 
        f = new TFile(rootDir.c_str(),"RECREATE");
        
        initHistograms();
    }
    
    std::cout<<"Done!"<<std::endl;
}


void CCProtonPi0_Interaction::initHistograms()
{
    Enu_1Track = new MnvH1D( "Enu_1Track","Reconstructed Beam Energy - 1 Track",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    Enu_1Track->GetXaxis()->SetTitle("Reconstructed E_{#nu} - 1Track [GeV]");
    Enu_1Track->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.beamE.get_width()));
    
    Enu_2Track = new MnvH1D( "Enu_2Track","Reconstructed Beam Energy - 2 Track",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    Enu_2Track->GetXaxis()->SetTitle("Reconstructed E_{#nu} - 2 Track [GeV]");
    Enu_2Track->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.beamE.get_width()));
    
    Enu_Cal = new MnvH1D( "Enu_Cal","Reconstructed Calorimetric Beam Energy",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    Enu_Cal->GetXaxis()->SetTitle("Reconstructed E_{#nu} - Calorimetry [GeV]");
    Enu_Cal->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.beamE.get_width()));
    
    q2 = new MnvH1D( "q2","Reconstructed Q^{2}",binList.q2.get_nBins(), binList.q2.get_min(), binList.q2.get_max() );
    q2->GetXaxis()->SetTitle("Reconstructed Q^{2} [GeV^{2}]");
    q2->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.q2.get_width()));
    
    w = new MnvH1D( "w","Reconstructed W",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max() );
    w->GetXaxis()->SetTitle("Reconstructed W [GeV]");
    w->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.w.get_width()));
    
    wSq = new MnvH1D( "wSq","Reconstructed W^{2}",binList.wSq.get_nBins(), binList.wSq.get_min(), binList.wSq.get_max() );
    wSq->GetXaxis()->SetTitle("Reconstructed W^{2} [GeV^{2}]");
    wSq->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.wSq.get_width()));
    
    nProngs_hist = new MnvH1D( "nProngs_hist","Number of Prongs",binList.multiplicity.get_nBins(), binList.multiplicity.get_min(), binList.multiplicity.get_max() );
    nProngs_hist->GetXaxis()->SetTitle("Number of Tracks");
    nProngs_hist->GetYaxis()->SetTitle("Events");
    
    deltaInvMass = new MnvH1D( "deltaInvMass","Reconstructed Delta+ Invariant Mass",binList.deltaInvMass.get_nBins(), binList.deltaInvMass.get_min(), binList.deltaInvMass.get_max() );
    deltaInvMass->GetXaxis()->SetTitle("Reconstructed #Delta^{+} Inv. Mass [MeV]");
    deltaInvMass->GetYaxis()->SetTitle(Form("Events / %3.2f [MeV] ",binList.deltaInvMass.get_width()));
    
    E_Unused_afterReco = new MnvH1D( "E_Unused_afterReco","Unused Cluster Energy after Reconstruction",binList.UnusedE.get_nBins(), binList.UnusedE.get_min(), binList.UnusedE.get_max() );
    E_Unused_afterReco->GetXaxis()->SetTitle("Unused E_{Visible} after Reconstruction [MeV]");
    E_Unused_afterReco->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.UnusedE.get_width()));
    
    E_Used_afterReco = new MnvH1D( "E_Used_afterReco","Used Cluster Energy after Reconstruction",binList.UsedE.get_nBins(), binList.UsedE.get_min(), binList.UsedE.get_max() );
    E_Used_afterReco->GetXaxis()->SetTitle("Used E_{Visible} after Reconstruction [MeV]");
    E_Used_afterReco->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.UsedE.get_width()));
    
    // MC Only Histograms
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

}

void CCProtonPi0_Interaction::writeHistograms()
{
    std::cout<<">> Writing "<<rootDir<<std::endl;
   
    f->cd();

    // Event Kinematics
    Enu_1Track->Write();
    Enu_2Track->Write();
    Enu_Cal->Write();
    q2->Write();
    w->Write();
    wSq->Write();
   
    // Reconstruction 
    E_Unused_afterReco->Write();
    E_Used_afterReco->Write();
    
    // Other Event Parameters 
    deltaInvMass->Write();
    nProngs_hist->Write();
   
    // MC Only Histograms
    mc_w_DIS->Write();
    mc_w_RES->Write();
    mc_w_CCQE->Write();
    
    final_mc_w_DIS->Write();
    final_mc_w_RES->Write();
    final_mc_w_CCQE->Write();
}



#endif
