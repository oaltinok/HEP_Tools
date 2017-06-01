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
        bin_P_Diff.setBin(100, -0.1,0.1);
        bin_KE.setBin(20, 0.0, 2.0);
        bin_trackLength.setBin(25,0.0,250.0);
        bin_trackKinked.setBin(2,0.0,2.0);
        bin_theta_Diff.setBin(40,-5.0,5.0);
        
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
        temp->GetYaxis()->SetTitle("Events/Bin");
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

    proton_theta_response = new MnvH2D( "proton_theta_response","Signal Proton Angle", binList.angle.get_nBins(), binList.angle.get_min(), binList.angle.get_max(),binList.angle.get_nBins(), binList.angle.get_min(), binList.angle.get_max() );
    proton_theta_response->GetXaxis()->SetTitle("Reconstructed #theta_{#mu} [degree]");
    proton_theta_response->GetYaxis()->SetTitle("True #theta_{#mu} [degree]");
    AddVertErrorBands_MC(proton_theta_response);
 
    theta_error = new TH1D( "theta_error","Error on cos(theta)",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    theta_error->GetXaxis()->SetTitle("(#theta_{Reco}-#theta_{True})/#theta_{True}");
    theta_error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));
 
    theta_Diff = new TH1D( "theta_Diff","Difference on Proton Theta",bin_theta_Diff.get_nBins(), bin_theta_Diff.get_min(), bin_theta_Diff.get_max() );
    theta_Diff->GetXaxis()->SetTitle("#theta_{Reco}-#theta_{True} [degree]");
    theta_Diff->GetYaxis()->SetTitle(Form("Events / %3.2f ",bin_theta_Diff.get_width()));

    proton_P_response = new MnvH2D( "proton_P_response","True vs Reconstructed Proton Momentum",bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max(), bin_P.get_nBins(), bin_P.get_min(), bin_P.get_max());
    proton_P_response->GetXaxis()->SetTitle("Reconstructed P_{p} [GeV]");
    proton_P_response->GetYaxis()->SetTitle("True P_{p} [GeV]");
    AddVertErrorBands_MC(proton_P_response);

    P_error = new TH1D( "P_error","Error on Proton Momentum",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    P_error->GetXaxis()->SetTitle("(P_{Reco}-P_{True})/P_{True}");
    P_error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));

    P_Diff = new TH1D( "P_Diff","Difference on Proton Energy",bin_P_Diff.get_nBins(), bin_P_Diff.get_min(), bin_P_Diff.get_max() );
    P_Diff->GetXaxis()->SetTitle("E_{Reco}-E_{True} [GeV]");
    P_Diff->GetYaxis()->SetTitle(Form("Events / %3.2f ",bin_P_Diff.get_width()));

    reco_E_true_E = new TH2D( "reco_E_true_E","True vs Reconstructed Proton Energy",bin_E.get_nBins(), bin_E.get_min(), bin_E.get_max(), bin_E.get_nBins(), bin_E.get_min(), bin_E.get_max());
    reco_E_true_E->GetXaxis()->SetTitle("Reconstructed E_{p} [GeV]");
    reco_E_true_E->GetYaxis()->SetTitle("True E_{p} [GeV]");

    E_error = new TH1D( "E_error","Error on Proton Energy",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    E_error->GetXaxis()->SetTitle("(E_{Reco}-E_{True})/E_{True}");
    E_error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));

    E_Diff = new TH1D( "E_Diff","Difference on Proton Energy",bin_E_Diff.get_nBins(), bin_E_Diff.get_min(), bin_E_Diff.get_max() );
    E_Diff->GetXaxis()->SetTitle("E_{Reco}-E_{True} [GeV]");
    E_Diff->GetYaxis()->SetTitle(Form("Events / %3.2f ",bin_E_Diff.get_width()));

    // dEdX Uncertainties 
    energy_shift_BetheBloch = new TH1D( "energy_shift_BetheBloch","Proton Energy Shift by BetheBloch", 100,-50,50);
    energy_shift_BetheBloch->GetXaxis()->SetTitle("Energy Shift [MeV]");
    energy_shift_BetheBloch->GetYaxis()->SetTitle("Events/Bin");

    energy_shift_Birks = new TH1D( "energy_shift_Birks","Proton Energy Shift by Birks", 100,-50,50);
    energy_shift_Birks->GetXaxis()->SetTitle("Energy Shift [MeV]");
    energy_shift_Birks->GetYaxis()->SetTitle("Events/Bin");

    energy_shift_MEU = new TH1D( "energy_shift_MEU","Proton Energy Shift by MEU", 100,-50,50);
    energy_shift_MEU->GetXaxis()->SetTitle("Energy Shift [MeV]");
    energy_shift_MEU->GetYaxis()->SetTitle("Events/Bin");

    energy_shift_Mass = new TH1D( "energy_shift_Mass","Proton Energy Shift by Mass", 100,-50,50);
    energy_shift_Mass->GetXaxis()->SetTitle("Energy Shift [MeV]");
    energy_shift_Mass->GetYaxis()->SetTitle("Events/Bin");

    energy_shift_Nominal = new TH1D( "energy_shift_Nominal","Proton Energy Shift by Nominal", 100,-50,50);
    energy_shift_Nominal->GetXaxis()->SetTitle("Energy Shift [MeV]");
    energy_shift_Nominal->GetYaxis()->SetTitle("Events/Bin");

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

    proton_theta_response->Write();
    theta_error->Write();
    theta_Diff->Write();

    proton_P_response->Write();
    P_error->Write();

    reco_E_true_E->Write();
    E_error->Write();

    E_Diff->Write();
    P_Diff->Write();

    energy_shift_BetheBloch->Write();
    energy_shift_Birks->Write();
    energy_shift_MEU->Write();
    energy_shift_Mass->Write();
    energy_shift_Nominal->Write();
    
    f->Close();
}

#endif


