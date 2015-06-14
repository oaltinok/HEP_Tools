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
    bin_blob_ndof.setBin(70,0,70);
    bin_blob_fval.setBin(100,0,5000);
    bin_blob_dEdx_doublet.setBin(2,0,2);
    bin_blob_dEdx.setBin(50,0,100);
    bin_blob_dEdx_nplane.setBin(40,0,40);
    bin_blob_dEdx_cluster_energy.setBin(100,0,500);
   
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


    g1_blob_ndof = new TH1D( "g1_blob_ndof","Blob Fit: ndof",bin_blob_ndof.get_nBins(), bin_blob_ndof.get_min(), bin_blob_ndof.get_max() );
    g1_blob_ndof->GetXaxis()->SetTitle("Blob Fit: ndof");
    g1_blob_ndof->GetYaxis()->SetTitle("N(Events)");

    g1_blob_fval = new TH1D( "g1_blob_fval","Blob Fit: fval",bin_blob_fval.get_nBins(), bin_blob_fval.get_min(), bin_blob_fval.get_max() );
    g1_blob_fval->GetXaxis()->SetTitle("Blob Fit: fval");
    g1_blob_fval->GetYaxis()->SetTitle("N(Events)");
    
    g1_blob_dEdx_doublet = new TH1D( "g1_blob_dEdx_doublet","blob_dEdx_doublet",bin_blob_dEdx_doublet.get_nBins(), bin_blob_dEdx_doublet.get_min(), bin_blob_dEdx_doublet.get_max() );
    g1_blob_dEdx_doublet->GetXaxis()->SetTitle("blob_dEdx_doublet");
    g1_blob_dEdx_doublet->GetYaxis()->SetTitle("N(Events)");

    g1_blob_dEdx_empty_plane = new TH1D( "g1_blob_dEdx_empty_plane","blob_dEdx_empty_plane",bin_blob_dEdx_doublet.get_nBins(), bin_blob_dEdx_doublet.get_min(), bin_blob_dEdx_doublet.get_max() );
    g1_blob_dEdx_empty_plane->GetXaxis()->SetTitle("blob_dEdx_empty_plane");
    g1_blob_dEdx_empty_plane->GetYaxis()->SetTitle("N(Events)");

    g1_blob_dEdx = new TH1D( "g1_blob_dEdx","blob_dEdx",bin_blob_dEdx.get_nBins(), bin_blob_dEdx.get_min(), bin_blob_dEdx.get_max() );
    g1_blob_dEdx->GetXaxis()->SetTitle("blob_dEdx");
    g1_blob_dEdx->GetYaxis()->SetTitle("N(Events)");
 
    g1_blob_dEdx1 = new TH1D( "g1_blob_dEdx1","blob_dEdx1",bin_blob_dEdx.get_nBins(), bin_blob_dEdx.get_min(), bin_blob_dEdx.get_max() );
    g1_blob_dEdx1->GetXaxis()->SetTitle("blob_dEdx1");
    g1_blob_dEdx1->GetYaxis()->SetTitle("N(Events)");

    g1_blob_dEdx_nplane = new TH1D( "g1_blob_dEdx_nplane","blob_dEdx_nplane",bin_blob_dEdx_nplane.get_nBins(), bin_blob_dEdx_nplane.get_min(), bin_blob_dEdx_nplane.get_max() );
    g1_blob_dEdx_nplane->GetXaxis()->SetTitle("blob_dEdx_nplane");
    g1_blob_dEdx_nplane->GetYaxis()->SetTitle("N(Events)");

    g1_blob_dEdx_cluster_energy = new TH1D( "g1_blob_dEdx_cluster_energy","blob_dEdx_cluster_energy",bin_blob_dEdx_cluster_energy.get_nBins(), bin_blob_dEdx_cluster_energy.get_min(), bin_blob_dEdx_cluster_energy.get_max() );
    g1_blob_dEdx_cluster_energy->GetXaxis()->SetTitle("blob_dEdx_cluster_energy");
    g1_blob_dEdx_cluster_energy->GetYaxis()->SetTitle("N(Events)");

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

