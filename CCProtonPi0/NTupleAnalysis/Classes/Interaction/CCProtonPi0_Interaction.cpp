/*
    See CCProtonPi0_Interaction.h header for Class Information
*/
#ifndef CCProtonPi0_Interaction_cpp
#define CCProtonPi0_Interaction_cpp

#include "CCProtonPi0_Interaction.h"

using namespace PlotUtils;

CCProtonPi0_Interaction::CCProtonPi0_Interaction(bool isModeReduce, bool isMC) : CCProtonPi0_NTupleAnalysis()
{
    std::cout<<"Initializing CCProtonPi0_Interaction"<<std::endl;
    
    if(isModeReduce){
        std::cout<<"\tNTuple Reduce Mode -- Will not create ROOT Files"<<std::endl;
    }else{
        if (isMC) rootDir = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + "Interaction.root";
        else rootDir = Folder_List::rootOut + Folder_List::Data + Folder_List::analyzed  + "Interaction.root";
        
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

        temp = new MnvH1D( Form("%s_%d","Enu_2Track",i),"Reconstructed Beam Energy - 2 Track",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed E_{#nu} - 2 Track [GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.beamE.get_width()));
        Enu_2Track.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","Enu_Cal",i),"Reconstructed Calorimetric Beam Energy",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed E_{#nu} - Calorimetry [GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.beamE.get_width()));
        Enu_Cal.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","q2",i),"Reconstructed Q^{2}",binList.q2.get_nBins(), binList.q2.get_min(), binList.q2.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed Q^{2} [GeV^{2}]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.q2.get_width()));
        q2.push_back(temp);    

        temp = new MnvH1D( Form("%s_%d","w",i),"Reconstructed W",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed W [GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.w.get_width()));
        w.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","wSq",i),"Reconstructed W^{2}",binList.wSq.get_nBins(), binList.wSq.get_min(), binList.wSq.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed W^{2} [GeV^{2}]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.wSq.get_width()));
        wSq.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","nProngs_hist",i),"Number of Prongs",binList.multiplicity.get_nBins(), binList.multiplicity.get_min(), binList.multiplicity.get_max() );
        temp->GetXaxis()->SetTitle("Number of Tracks");
        temp->GetYaxis()->SetTitle("Events");
        nProngs_hist.push_back(temp);

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

    }
    
    // MC Only Histograms
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

    for (int i = 0; i < nHistograms; i++){
        // Event Kinematics
        Enu_1Track[i]->Write();
        Enu_2Track[i]->Write();
        Enu_Cal[i]->Write();
        q2[i]->Write();
        w[i]->Write();
        wSq[i]->Write();

        // Reconstruction 
        E_Unused_afterReco[i]->Write();
        E_Used_afterReco[i]->Write();

        // Other Event Parameters 
        deltaInvMass[i]->Write();
        nProngs_hist[i]->Write();
    }
    // MC Only Histograms
    final_mc_w_DIS->Write();
    final_mc_w_RES->Write();
    final_mc_w_CCQE->Write();
    
    f->Close();
}



#endif
