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
        if (isMC) rootDir = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + "Muon.root";
        else rootDir = Folder_List::rootOut + Folder_List::Data + Folder_List::analyzed + "Muon.root";

        std::cout<<"\tRoot File: "<<rootDir<<std::endl;
        
        // Create Root File 
        f = new TFile(rootDir.c_str(),"RECREATE");

        // Initialize Bins
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
        temp = new MnvH1D( Form("%s_%d","E",i),"Reconstructed Muon Energy",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed E_{#mu} [GeV]");
        temp->GetYaxis()->SetTitle(Form("Muons / %3.1f [GeV] ",bin_P.get_width()));
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

    f->Close();
}


#endif



