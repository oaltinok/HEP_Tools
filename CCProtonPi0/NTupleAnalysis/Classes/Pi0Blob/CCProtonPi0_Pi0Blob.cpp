/*
    See CCProtonPi0_Pi0Blob.h header for Class Information
*/
#ifndef CCProtonPi0_Pi0Blob_cpp
#define CCProtonPi0_Pi0Blob_cpp

#include "CCProtonPi0_Pi0Blob.h"

using namespace std;

CCProtonPi0_Pi0Blob::CCProtonPi0_Pi0Blob(int nMode) : CCProtonPi0_NTupleAnalysis(nMode)
{
    cout<<"Initializing CCProtonPi0_Pi0Blob"<<endl;

    rootDir = Folder_List::output + Folder_List::rootOut + branchDir + "Pi0Blob.root";
    
    cout<<"\tRoot File: "<<rootDir<<endl;
 
    // Create Root File 
    f = new TFile(rootDir.c_str(),"RECREATE");

    // Initialize Bins
    bin_blob_nclusters.setBin(40,0.0,40.0);
    bin_blob_minsep.setBin(200,0.0,2000.0);
    bin_blob_energy.setBin(100,0.0,1000);
    bin_blob_ndigits.setBin(100,0,100);
   
    initHistograms();
  
    cout<<"Done!"<<endl;
}

void CCProtonPi0_Pi0Blob::initHistograms()
{
    // Gamma1 
    g1_blob_ndigits = new TH1D( "g1_blob_ndigits","Photon N(Digits)",bin_blob_ndigits.get_nBins(), bin_blob_ndigits.get_min(), bin_blob_ndigits.get_max() );
    g1_blob_ndigits->GetXaxis()->SetTitle("Photon N(Digits)");
    g1_blob_ndigits->GetYaxis()->SetTitle("N(Events)");
    
    g1_blob_nclusters = new TH1D( "g1_blob_nclusters","Leading Photon N(Clusters)",bin_blob_nclusters.get_nBins(), bin_blob_nclusters.get_min(), bin_blob_nclusters.get_max() );
    g1_blob_nclusters->GetXaxis()->SetTitle("Leading Photon N(Clusters)");
    g1_blob_nclusters->GetYaxis()->SetTitle("N(Events)");
    
    g1_blob_energy = new TH1D( "g1_blob_energy","Leading Photon Energy",bin_blob_energy.get_nBins(), bin_blob_energy.get_min(), bin_blob_energy.get_max() );
    g1_blob_energy->GetXaxis()->SetTitle("Blob Energy [MeV]");
    g1_blob_energy->GetYaxis()->SetTitle("N(Events)");

    g1_blob_minsep = new TH1D( "g1_blob_minsep","Blob Min Separation from Vertex",bin_blob_minsep.get_nBins(), bin_blob_minsep.get_min(), bin_blob_minsep.get_max() );
    g1_blob_minsep->GetXaxis()->SetTitle("Blob Min Separation from Vertex");
    g1_blob_minsep->GetYaxis()->SetTitle("N(Events)");
 
     // Gamma 2
    g2_blob_ndigits = new TH1D( "g2_blob_ndigits","Photon N(Digits)",bin_blob_ndigits.get_nBins(), bin_blob_ndigits.get_min(), bin_blob_ndigits.get_max() );
    g2_blob_ndigits->GetXaxis()->SetTitle("Photon N(Digits)");
    g2_blob_ndigits->GetYaxis()->SetTitle("N(Events)");
    
    g2_blob_nclusters = new TH1D( "g2_blob_nclusters","Second Photon N(Clusters)",bin_blob_nclusters.get_nBins(), bin_blob_nclusters.get_min(), bin_blob_nclusters.get_max() );
    g2_blob_nclusters->GetXaxis()->SetTitle("Second Photon N(Clusters)");
    g2_blob_nclusters->GetYaxis()->SetTitle("N(Events)");
  
    g2_blob_energy = new TH1D( "g2_blob_energy","Second Photon Energy",bin_blob_energy.get_nBins(), bin_blob_energy.get_min(), bin_blob_energy.get_max() );
    g2_blob_energy->GetXaxis()->SetTitle("Blob Energy [MeV]");
    g2_blob_energy->GetYaxis()->SetTitle("N(Events)");

    g2_blob_minsep = new TH1D( "g2_blob_minsep","Blob Min Separation from Vertex",bin_blob_minsep.get_nBins(), bin_blob_minsep.get_min(), bin_blob_minsep.get_max() );
    g2_blob_minsep->GetXaxis()->SetTitle("Blob Min Separation from Vertex");
    g2_blob_minsep->GetYaxis()->SetTitle("N(Events)");

}

void CCProtonPi0_Pi0Blob::write_RootFile()
{
    cout<<">> Writing "<<rootDir<<endl;
    f->Write();
}


#endif

