/*
    See CCProtonPi0_Proton.h header for Class Information
*/

#ifndef CCProtonPi0_Proton_cpp
#define CCProtonPi0_Proton_cpp

#include "CCProtonPi0_Proton.h"

using namespace PlotUtils;

CCProtonPi0_Proton::CCProtonPi0_Proton(bool isModeReduce, bool isMC) : CCProtonPi0_Particle()
{
    std::cout<<"Initializing CCProtonPi0_Proton"<<std::endl;
        
    if (isModeReduce){
        std::cout<<"\tNTuple Reduce Mode -- Will not create ROOT Files"<<std::endl;
    }else{
        // File Locations
        if (isMC) rootDir = Folder_List::rootDir_Proton_mc;
        else rootDir = Folder_List::rootDir_Proton_data;
        
        std::cout<<"\tRoot File: "<<rootDir<<std::endl;
        
        // Create Root File 
        f = new TFile(rootDir.c_str(),"RECREATE");

        // Initialize Bins
        bin_E.setBin(25, 0.5 ,3.0);
        bin_E_Diff.setBin(100, -0.1,0.1);
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
    MnvH1D* temp = NULL;

    for(int i = 0; i < nHistograms; i++){
        // Unique Histograms
        temp = new MnvH1D( Form("%s_%d","trackLength",i),"Proton Track Length",bin_trackLength.get_nBins(), bin_trackLength.get_min(), bin_trackLength.get_max() );
        temp->GetXaxis()->SetTitle("Proton Track Length [cm]");
        temp->GetYaxis()->SetTitle(Form("Protons / %3.1f cm ",bin_trackLength.get_width()));
        trackLength.push_back(temp);    

        temp = new MnvH1D( Form("%s_%d","trackKinked",i),"Proton Track Kinked or NOT",bin_trackKinked.get_nBins(), bin_trackKinked.get_min(), bin_trackKinked.get_max() );
        temp->GetXaxis()->SetTitle("Proton Track Kinked or NOT");
        temp->GetYaxis()->SetTitle("N(Events)");
        trackKinked.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","partScore",i),"Proton Particle Score (LLR)",binList.particleScore_LLR.get_nBins(), binList.particleScore_LLR.get_min(), binList.particleScore_LLR.get_max() );
        temp->GetXaxis()->SetTitle("Particle Score");
        temp->GetYaxis()->SetTitle(Form("Protons / %3.1f ",binList.particleScore_LLR.get_width()));
        partScore.push_back(temp);

        // Standard Histograms
        temp = new MnvH1D( Form("%s_%d","E",i),"Reconstructed Proton Energy",bin_E.get_nBins(), bin_E.get_min(), bin_E.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed E_{p} [GeV]");
        temp->GetYaxis()->SetTitle(Form("Protons / %3.1f [GeV] ",bin_E.get_width()));
        E.push_back(temp); 

        temp = new MnvH1D( Form("%s_%d","P",i),"Reconstructed Proton Momentum",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed P_{p} [GeV]");
        temp->GetYaxis()->SetTitle(Form("Protons / %3.1f [GeV] ",bin_P.get_width()));
        P.push_back(temp); 

        temp = new MnvH1D( Form("%s_%d","KE",i),"Reconstructed Proton Kinetic Energy",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed T_{p} [GeV]");
        temp->GetYaxis()->SetTitle(Form("Protons / %3.1f [GeV] ",bin_P.get_width()));
        KE.push_back(temp); 

        temp = new MnvH1D( Form("%s_%d","theta",i),"Reconstructed #theta_{p}",binList.angle.get_nBins(), binList.angle.get_min(), binList.angle.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed #theta_{p} [Degree]");
        temp->GetYaxis()->SetTitle(Form("Protons / %3.1f [Degree]",binList.angle.get_width()));
        theta.push_back(temp); 

        temp = new MnvH1D( Form("%s_%d","phi",i),"Reconstructed #phi_{p}",binList.angle.get_nBins(), binList.angle.get_min(), binList.angle.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed #phi_{p} [Degree]");
        temp->GetYaxis()->SetTitle(Form("Protons / %3.1f [Degree]",binList.angle.get_width()));
        phi.push_back(temp); 
    }
    
    reco_P_true_P = new TH2D( "reco_P_true_P","True vs Reconstructed Proton Momentum",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max(), bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max());
    reco_P_true_P->GetXaxis()->SetTitle("Reconstructed P_{p} [GeV]");
    reco_P_true_P->GetYaxis()->SetTitle("True P_{p} [GeV]");

    P_error = new TH1D( "P_error","Error on Proton Momentum",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    P_error->GetXaxis()->SetTitle("(P_{Reco}-P_{True})/P_{True}");
    P_error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));

    reco_E_true_E = new TH2D( "reco_E_true_E","True vs Reconstructed Proton Energy",bin_E.get_nBins(), bin_E.get_min(), bin_E.get_max(), bin_E.get_nBins(), bin_E.get_min(), bin_E.get_max());
    reco_E_true_E->GetXaxis()->SetTitle("Reconstructed E_{p} [GeV]");
    reco_E_true_E->GetYaxis()->SetTitle("True E_{p} [GeV]");

    E_error = new TH1D( "E_error","Error on Proton Energy",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    E_error->GetXaxis()->SetTitle("(E_{Reco}-E_{True})/E_{True}");
    E_error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));

    E_Diff = new TH1D( "E_Diff","Difference on Proton Energy",bin_E_Diff.get_nBins(), bin_E_Diff.get_min(), bin_E_Diff.get_max() );
    E_Diff->GetXaxis()->SetTitle("E_{Reco}-E_{True} [GeV]");
    E_Diff->GetYaxis()->SetTitle(Form("Events / %3.2f ",bin_E_Diff.get_width()));

}

void CCProtonPi0_Proton::writeHistograms()
{
    std::cout<<">> Writing "<<rootDir<<std::endl;
    f->cd();
    for (int i = 0; i < nHistograms; i++){
        trackLength[i]->Write();
        trackKinked[i]->Write();
        partScore[i]->Write();
        E[i]->Write();
        P[i]->Write();
        KE[i]->Write();
        theta[i]->Write();
        phi[i]->Write();
    }
    
    reco_P_true_P->Write();
    P_error->Write();

    reco_E_true_E->Write();
    E_error->Write();

    E_Diff->Write();

    f->Close();
}

#endif


