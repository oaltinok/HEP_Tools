/*
    See CCProtonPi0_Muon.h header for Class Information
*/
#ifndef CCProtonPi0_Muon_cpp
#define CCProtonPi0_Muon_cpp

#include "CCProtonPi0_Muon.h"

using namespace PlotUtils;

CCProtonPi0_Muon::CCProtonPi0_Muon(int nMode, bool isMC) : CCProtonPi0_Particle(nMode)
{
    std::cout<<"Initializing CCProtonPi0_Muon"<<std::endl;
    
    if(nMode == 0){
        std::cout<<"\tNTuple Reduce Mode -- Will not create ROOT Files"<<std::endl;
    }else{
        // File Locations
        if (isMC) rootDir = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + branchDir + "Muon.root";
        else rootDir = Folder_List::rootOut + Folder_List::Data + Folder_List::analyzed + branchDir + "Muon.root";

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
    E = new MnvH1D( "E","Reconstructed Muon Energy",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    E->GetXaxis()->SetTitle("Reconstructed Muon Energy [GeV]");
    E->GetYaxis()->SetTitle(Form("Number of Muons / %3.1f [GeV] ",bin_P.get_width()));
    
    P = new MnvH1D( "P","Reconstructed Muon Momentum",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    P->GetXaxis()->SetTitle("Reconstructed Muon Momentum [GeV]");
    P->GetYaxis()->SetTitle(Form("Number of Muons / %3.1f [GeV] ",bin_P.get_width()));
    
    KE = new MnvH1D( "KE","Reconstructed Muon Kinetic Energy",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    KE->GetXaxis()->SetTitle("Reconstructed Muon Kinetic Energy [GeV]");
    KE->GetYaxis()->SetTitle(Form("Number of Muons / %3.1f [GeV]",bin_P.get_width()));
   
    theta = new MnvH1D( "theta","Reconstructed Muon Theta",bin_muonTheta.get_nBins(), bin_muonTheta.get_min(), bin_muonTheta.get_max() );
    theta->GetXaxis()->SetTitle("Theta [Degree]");
    theta->GetYaxis()->SetTitle(Form("Number of Muons / %3.1f [Degree] ",bin_muonTheta.get_width()));

    phi = new MnvH1D( "phi","Reconstructed Muon Phi",binList.angle.get_nBins(), binList.angle.get_min(), binList.angle.get_max() );
    phi->GetXaxis()->SetTitle("Phi [Degree]");
    phi->GetYaxis()->SetTitle(Form("Number of Muons / %3.1f [Degree]",binList.angle.get_width()));
}

void CCProtonPi0_Muon::writeHistograms()
{
    std::cout<<">> Writing "<<rootDir<<std::endl;
    f->cd();
    E->Write();
    P->Write();
    KE->Write();
    theta->Write();
    phi->Write();
}


#endif



