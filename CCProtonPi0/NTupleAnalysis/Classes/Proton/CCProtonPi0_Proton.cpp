/*
    See CCProtonPi0_Proton.h header for Class Information
*/

#ifndef CCProtonPi0_Proton_cpp
#define CCProtonPi0_Proton_cpp

#include "CCProtonPi0_Proton.h"

using namespace PlotUtils;

CCProtonPi0_Proton::CCProtonPi0_Proton(int nMode, bool isMC) : CCProtonPi0_Particle(nMode)
{
    std::cout<<"Initializing CCProtonPi0_Proton"<<std::endl;
        
    if(nMode == 0){
        std::cout<<"\tNTuple Reduce Mode -- Will not create ROOT Files"<<std::endl;
    }else{
        // File Locations
        if (isMC) rootDir = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + branchDir + "Proton.root";
        else rootDir = Folder_List::rootOut + Folder_List::Data + Folder_List::analyzed + branchDir + "Proton.root";
        
        std::cout<<"\tRoot File: "<<rootDir<<std::endl;
        
        // Create Root File 
        f = new TFile(rootDir.c_str(),"RECREATE");

        // Initialize Bins
        bin_E.setBin(25, 0.5 ,3.0);
        bin_P.setBin(20, 0.0, 2.0);
        bin_KE.setBin(20, 0.0, 2.0);
        bin_trackLength.setBin(25,0.0,250.0);
        bin_trackKinked.setBin(2,0.0,2.0);
        
        initHistograms();        
    }
    std::cout<<"Done!"<<std::endl;
}

void CCProtonPi0_Proton::initHistograms()
{
    // Unique Histograms
    trackLength = new MnvH1D( "trackLength","Proton Track Length",bin_trackLength.get_nBins(), bin_trackLength.get_min(), bin_trackLength.get_max() );
    trackLength->GetXaxis()->SetTitle("Proton Track Length [cm]");
    trackLength->GetYaxis()->SetTitle(Form("Protons / %3.1f cm ",bin_trackLength.get_width()));
    
    trackKinked = new MnvH1D( "trackKinked","Proton Track Kinked or NOT",bin_trackKinked.get_nBins(), bin_trackKinked.get_min(), bin_trackKinked.get_max() );
    trackKinked->GetXaxis()->SetTitle("Proton Track Kinked or NOT");
    trackKinked->GetYaxis()->SetTitle("N(Events)");

    partScore = new MnvH1D( "partScore","Proton Particle Score (LLR)",binList.particleScore_LLR.get_nBins(), binList.particleScore_LLR.get_min(), binList.particleScore_LLR.get_max() );
    partScore->GetXaxis()->SetTitle("Particle Score");
    partScore->GetYaxis()->SetTitle(Form("Protons / %3.1f ",binList.particleScore_LLR.get_width()));
   
    // Standard Histograms
    E = new MnvH1D( "E","Reconstructed Proton Energy",bin_E.get_nBins(), bin_E.get_min(), bin_E.get_max() );
    E->GetXaxis()->SetTitle("Reconstructed E_{p} [GeV]");
    E->GetYaxis()->SetTitle(Form("Protons / %3.1f [GeV] ",bin_E.get_width()));
    
    P = new MnvH1D( "P","Reconstructed Proton Momentum",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    P->GetXaxis()->SetTitle("Reconstructed P_{p} [GeV]");
    P->GetYaxis()->SetTitle(Form("Protons / %3.1f [GeV] ",bin_P.get_width()));
       
    KE = new MnvH1D( "KE","Reconstructed Proton Kinetic Energy",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
    KE->GetXaxis()->SetTitle("Reconstructed T_{p} [GeV]");
    KE->GetYaxis()->SetTitle(Form("Protons / %3.1f [GeV] ",bin_P.get_width()));
    
    theta = new MnvH1D( "theta","Reconstructed #theta_{p}",binList.angle.get_nBins(), binList.angle.get_min(), binList.angle.get_max() );
    theta->GetXaxis()->SetTitle("Reconstructed #theta_{p} [Degree]");
    theta->GetYaxis()->SetTitle(Form("Protons / %3.1f [Degree]",binList.angle.get_width()));
    
    phi = new MnvH1D( "phi","Reconstructed #phi_{p}",binList.angle.get_nBins(), binList.angle.get_min(), binList.angle.get_max() );
    phi->GetXaxis()->SetTitle("Reconstructed #phi_{p} [Degree]");
    phi->GetYaxis()->SetTitle(Form("Protons / %3.1f [Degree]",binList.angle.get_width()));

    reco_P_true_P = new TH2D( "reco_P_true_P","True vs Reconstructed Proton Momentum",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max(), bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max());
    reco_P_true_P->GetXaxis()->SetTitle("Reconstructed P_{p} [GeV]");
    reco_P_true_P->GetYaxis()->SetTitle("True P_{p} [GeV]");

    P_error = new TH1D( "P_error","Error on Proton Momentum",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    P_error->GetXaxis()->SetTitle("(P_{Reco}-P_{True})/P_{True}");
    P_error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));
}

void CCProtonPi0_Proton::writeHistograms()
{
    std::cout<<">> Writing "<<rootDir<<std::endl;
    f->cd();
    trackLength->Write();
    trackKinked->Write();
    partScore->Write();
    E->Write();
    P->Write();
    KE->Write();
    theta->Write();
    phi->Write();
    reco_P_true_P->Write();
    P_error->Write();
}

#endif


