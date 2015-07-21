#ifndef CCProtonPi0_MichelTool_cpp
#define CCProtonPi0_MichelTool_cpp

#include "CCProtonPi0_MichelTool.h"

using namespace std;

CCProtonPi0_MichelTool::CCProtonPi0_MichelTool(int nMode, bool isMC) : CCProtonPi0_NTupleAnalysis(nMode)
{
    cout<<"Initializing CCProtonPi0_MichelTool"<<endl;
    
    if(nMode == 0){
        cout<<"\tNTuple Reduce Mode -- Will not create ROOT Files"<<endl;
    }else{
        // File Locations
        if (isMC) rootDir = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + branchDir + "MichelTool.root";
        else rootDir = Folder_List::rootOut + Folder_List::Data + Folder_List::analyzed + branchDir + "MichelTool.root";      

        cout<<"\tRoot File: "<<rootDir<<endl;
        
        // Create Root File 
        f = new TFile(rootDir.c_str(),"RECREATE");
        
        initHistograms();
        
        N_trueMichel_before = 0.0;
        N_trueMichel_after = 0.0;
        N_trueMichel_afterAll = 0.0;
        N_noMichel_before = 0.0;
        N_noMichel_after = 0.0;
        N_detectedMichel_true = 0.0;
        N_detectedMichel_fake = 0.0;
        N_missedMichel_true = 0.0;
        N_missedMichel_fake = 0.0;
        N_detectedMichel_true_signal = 0.0;
        N_detectedMichel_fake_signal = 0.0;
        N_missedMichel_true_signal = 0.0;
        N_missedMichel_fake_signal = 0.0;
        N_selected_detectedMichel_true = 0.0;
        N_selected_detectedMichel_fake = 0.0;
        N_selected_missedMichel_true = 0.0;
        N_selected_missedMichel_fake = 0.0;
    }
    
    cout<<"Done!"<<endl;
}


void CCProtonPi0_MichelTool::initHistograms()
{
    trueMichel_dist_vtx = new TH1D( "trueMichel_dist_vtx","TRUE Michel Distance to Vertex (Reconstructed)",binList.michelMuon_end_dist_vtx.get_nBins(), binList.michelMuon_end_dist_vtx.get_min(), binList.michelMuon_end_dist_vtx.get_max() );
    trueMichel_dist_vtx->GetXaxis()->SetTitle("TRUE Michel Distance to Vertex (Reconstructed) [mm]");
    trueMichel_dist_vtx->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelMuon_end_dist_vtx.get_width()));
    
    fakeMichel_dist_vtx = new TH1D( "fakeMichel_dist_vtx","FAKE Michel Distance to Vertex (Reconstructed)",binList.michelMuon_end_dist_vtx.get_nBins(), binList.michelMuon_end_dist_vtx.get_min(), binList.michelMuon_end_dist_vtx.get_max() );
    fakeMichel_dist_vtx->GetXaxis()->SetTitle("FAKE Michel Distance to Vertex (Reconstructed) [mm]");
    fakeMichel_dist_vtx->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelMuon_end_dist_vtx.get_width()));

    trueMichel_dist_end_point = new TH1D( "trueMichel_dist_end_point","TRUE Michel Distance to Track End Point (Reconstructed)",binList.michelMuon_end_dist_vtx.get_nBins(), binList.michelMuon_end_dist_vtx.get_min(), binList.michelMuon_end_dist_vtx.get_max() );
    trueMichel_dist_end_point->GetXaxis()->SetTitle("TRUE Michel Distance to Track End Point (Reconstructed) [mm]");
    trueMichel_dist_end_point->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelMuon_end_dist_vtx.get_width()));
    
    fakeMichel_dist_end_point = new TH1D( "fakeMichel_dist_end_point","FAKE Michel Distance to Track End Point (Reconstructed)",binList.michelMuon_end_dist_vtx.get_nBins(), binList.michelMuon_end_dist_vtx.get_min(), binList.michelMuon_end_dist_vtx.get_max() );
    fakeMichel_dist_end_point->GetXaxis()->SetTitle("FAKE Michel Distance to Track End Point (Reconstructed) [mm]");
    fakeMichel_dist_end_point->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelMuon_end_dist_vtx.get_width()));
    
    trueMichel_energy = new TH1D( "trueMichel_energy","TRUE Michel Electron Energy",binList.michelMuon_P.get_nBins(), binList.michelMuon_P.get_min(), binList.michelMuon_P.get_max() );
    trueMichel_energy->GetXaxis()->SetTitle("TRUE Michel Electron Energy [MeV]");
    trueMichel_energy->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelMuon_P.get_width()));
    
    fakeMichel_energy = new TH1D( "fakeMichel_energy","FAKE Michel Electron Energy",binList.michelMuon_P.get_nBins(), binList.michelMuon_P.get_min(), binList.michelMuon_P.get_max() );
    fakeMichel_energy->GetXaxis()->SetTitle("FAKE Michel Electron Energy [MeV]");
    fakeMichel_energy->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelMuon_P.get_width()));
    
    trueMichel_end_Z = new TH1D( "trueMichel_end_Z","TRUE Michel Prong end Z",binList.vertex_z.get_nBins(), binList.vertex_z.get_min(), binList.vertex_z.get_max() );
    trueMichel_end_Z->GetXaxis()->SetTitle("z = 4293 Target, #bf{z = 5810 Interaction Region}, z = 8614 ECAL, z = 9088 HCAL");
    trueMichel_end_Z->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.vertex_z.get_width()));
    
    fakeMichel_end_Z = new TH1D( "fakeMichel_end_Z","FAKE Michel Prong end Z",binList.vertex_z.get_nBins(), binList.vertex_z.get_min(), binList.vertex_z.get_max() );
    fakeMichel_end_Z->GetXaxis()->SetTitle("z = 4293 Target, #bf{z = 5810 Interaction Region}, z = 8614 ECAL, z = 9088 HCAL");
    fakeMichel_end_Z->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.vertex_z.get_width()));
    
    trueMichel_time_diff = new TH1D( "trueMichel_time_diff","TRUE Michel Time Difference",binList.michel_time_diff.get_nBins(), binList.michel_time_diff.get_min(), binList.michel_time_diff.get_max() );
    trueMichel_time_diff->GetXaxis()->SetTitle("TRUE Michel Time Difference");
    trueMichel_time_diff->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michel_time_diff.get_width()));    

    fakeMichel_time_diff = new TH1D( "fakeMichel_time_diff","FAKE Michel Time Difference",binList.michel_time_diff.get_nBins(), binList.michel_time_diff.get_min(), binList.michel_time_diff.get_max() );
    fakeMichel_time_diff->GetXaxis()->SetTitle("TRUE Michel Time Difference");
    fakeMichel_time_diff->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michel_time_diff.get_width()));
   
    trueMichel_end_Z_vtx_Z = new TH1D( "trueMichel_end_Z_vtx_Z","TRUE Michel Muon End Point Z - Vertex Z",binList.michelMuon_Z_vtx.get_nBins(), binList.michelMuon_Z_vtx.get_min(), binList.michelMuon_Z_vtx.get_max() );
    trueMichel_end_Z_vtx_Z->GetXaxis()->SetTitle("Michel Muon End Point Z - Vertex Z");
    trueMichel_end_Z_vtx_Z->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelMuon_Z_vtx.get_width()));
    
    fakeMichel_end_Z_vtx_Z = new TH1D( "fakeMichel_end_Z_vtx_Z","FAKE Michel Muon End Point Z - Vertex Z",binList.michelMuon_Z_vtx.get_nBins(), binList.michelMuon_Z_vtx.get_min(), binList.michelMuon_Z_vtx.get_max() );
    fakeMichel_end_Z_vtx_Z->GetXaxis()->SetTitle("Michel Muon End Point Z - Vertex Z");
    fakeMichel_end_Z_vtx_Z->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelMuon_Z_vtx.get_width()));    
    
    N_michelElectrons = new TH1D( "N_michelElectrons","Number of Michel Electrons in an Event",binList.multiplicity.get_nBins(), binList.multiplicity.get_min(), binList.multiplicity.get_max() );
    N_michelElectrons->GetXaxis()->SetTitle("Number of Michel Electrons");
    N_michelElectrons->GetYaxis()->SetTitle("N(Events)");
    
    michelElectron_E[0] = new TH1D( "vertex_michelElectron_E","Michel Electron Energy",binList.michelMuon_P.get_nBins(), binList.michelMuon_P.get_min(), binList.michelMuon_P.get_max() );
    michelElectron_E[0]->GetXaxis()->SetTitle("Michel Electron Energy [MeV]");
    michelElectron_E[0]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelMuon_P.get_width()));
    
    michelMuon_P[0] = new TH1D( "vertex_michelMuon_P","Michel Muon Momentum",binList.michelMuon_P.get_nBins(), binList.michelMuon_P.get_min(), binList.michelMuon_P.get_max() );
    michelMuon_P[0]->GetXaxis()->SetTitle("Michel Muon Momentum [MeV]");
    michelMuon_P[0]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelMuon_P.get_width()));
    
    michelMuon_end_dist_vtx[0] = new TH1D( "vertex_michelMuon_end_dist_vtx","Michel Muon End Point Distance to Vertex",binList.michelMuon_end_dist_vtx.get_nBins(), binList.michelMuon_end_dist_vtx.get_min(), binList.michelMuon_end_dist_vtx.get_max() );
    michelMuon_end_dist_vtx[0]->GetXaxis()->SetTitle("Michel Muon End Point Distance to Vertex [mm]");
    michelMuon_end_dist_vtx[0]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelMuon_end_dist_vtx.get_width()));

    michelMuon_length[0] = new TH1D( "vertex_michelMuon_length","Michel Muon Track Length",binList.michelMuon_length.get_nBins(), binList.michelMuon_length.get_min(), binList.michelMuon_length.get_max() );
    michelMuon_length[0]->GetXaxis()->SetTitle("Michel Muon Track Length [mm]");
    michelMuon_length[0]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelMuon_length.get_width()));
    
    michelMuon_X_Y[0] = new TH2D( "vertex_michelMuon_X_Y","Michel Muon End Point X vs Y",   binList.vertex_x_y.get_nBins(), binList.vertex_x_y.get_min(), binList.vertex_x_y.get_max(),
    binList.vertex_x_y.get_nBins(), binList.vertex_x_y.get_min(), binList.vertex_x_y.get_max());
    michelMuon_X_Y[0]->GetXaxis()->SetTitle("Michel Muon End Point X [mm]");
    michelMuon_X_Y[0]->GetYaxis()->SetTitle("Michel Muon End Point [mm]");
    
    michelMuon_Z[0] = new TH1D( "vertex_michelMuon_Z","Michel Muon End Point Z",binList.vertex_z.get_nBins(), binList.vertex_z.get_min(), binList.vertex_z.get_max() );
    michelMuon_Z[0]->GetXaxis()->SetTitle("z = 4293 Target, #bf{z = 5810 Interaction Region}, z = 8614 ECAL, z = 9088 HCAL");
    michelMuon_Z[0]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.vertex_z.get_width()));
    
    michelMuon_Z_vtx[0] = new TH1D( "vertex_michelMuon_Z_vtx","Michel Muon End Point Z - Vertex Z",binList.michelMuon_Z_vtx.get_nBins(), binList.michelMuon_Z_vtx.get_min(), binList.michelMuon_Z_vtx.get_max() );
    michelMuon_Z_vtx[0]->GetXaxis()->SetTitle("Michel Muon End Point Z - Vertex Z");
    michelMuon_Z_vtx[0]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelMuon_Z_vtx.get_width()));
    
    michelPion_P[0] = new TH1D( "vertex_michelPion_P","Michel Pion Momentum",binList.michelPion_P.get_nBins(), binList.michelPion_P.get_min(), binList.michelPion_P.get_max() );
    michelPion_P[0]->GetXaxis()->SetTitle("Michel Pion Momentum [MeV]");
    michelPion_P[0]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelPion_P.get_width()));
    
    michelPion_begin_dist_vtx[0] = new TH1D( "vertex_michelPion_begin_dist_vtx","Michel Pion Inital Point Distance to Vertex",binList.michelPion_begin_dist_vtx.get_nBins(), binList.michelPion_begin_dist_vtx.get_min(), binList.michelPion_begin_dist_vtx.get_max() );
    michelPion_begin_dist_vtx[0]->GetXaxis()->SetTitle("Michel Pion Inital Point Distance to Vertex [mm]");
    michelPion_begin_dist_vtx[0]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelPion_begin_dist_vtx.get_width()));

    michelPion_length[0] = new TH1D( "vertex_michelPion_length","Michel Pion Track Length",binList.michelPion_length.get_nBins(), binList.michelPion_length.get_min(), binList.michelPion_length.get_max() );
    michelPion_length[0]->GetXaxis()->SetTitle("Michel Pion Track Length [mm]");
    michelPion_length[0]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelPion_length.get_width()));
    
    michelMuon_dist_michelPion_length[0] = new TH2D( "vertex_michelMuon_dist_michelPion_length","Michel Pion Length vs Michel Muon Distance to Vertex",binList.michelMuon_end_dist_vtx.get_nBins(), binList.michelMuon_end_dist_vtx.get_min(), binList.michelMuon_end_dist_vtx.get_max(),binList.michelPion_begin_dist_vtx.get_nBins(), binList.michelPion_begin_dist_vtx.get_min(), binList.michelPion_begin_dist_vtx.get_max());
    michelMuon_dist_michelPion_length[0]->GetXaxis()->SetTitle("Michel Muon End Point Distance to Vertex [mm]");
    michelMuon_dist_michelPion_length[0]->GetYaxis()->SetTitle("Michel Pion Length [mm]");
    
    // Michel Study - Found Events by Track End Point Michel
    michelElectron_E[1] = new TH1D( "track_michelElectron_E","Michel Electron Energy",binList.michelMuon_P.get_nBins(), binList.michelMuon_P.get_min(), binList.michelMuon_P.get_max() );
    michelElectron_E[1]->GetXaxis()->SetTitle("Michel Electron Energy [MeV]");
    michelElectron_E[1]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelMuon_P.get_width()));
    
    michelMuon_P[1] = new TH1D( "track_michelMuon_P","Michel Muon Momentum",binList.michelMuon_P.get_nBins(), binList.michelMuon_P.get_min(), binList.michelMuon_P.get_max() );
    michelMuon_P[1]->GetXaxis()->SetTitle("Michel Muon Momentum [MeV]");
    michelMuon_P[1]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelMuon_P.get_width()));
    
    michelMuon_end_dist_vtx[1] = new TH1D( "track_michelMuon_end_dist_vtx","Michel Muon End Point Distance to Vertex",binList.michelMuon_end_dist_vtx.get_nBins(), binList.michelMuon_end_dist_vtx.get_min(), binList.michelMuon_end_dist_vtx.get_max() );
    michelMuon_end_dist_vtx[1]->GetXaxis()->SetTitle("Michel Muon End Point Distance to Vertex [mm]");
    michelMuon_end_dist_vtx[1]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelMuon_end_dist_vtx.get_width()));

    michelMuon_length[1] = new TH1D( "track_michelMuon_length","Michel Muon Track Length",binList.michelMuon_length.get_nBins(), binList.michelMuon_length.get_min(), binList.michelMuon_length.get_max() );
    michelMuon_length[1]->GetXaxis()->SetTitle("Michel Muon Track Length [mm]");
    michelMuon_length[1]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelMuon_length.get_width()));
    
    michelMuon_X_Y[1] = new TH2D( "track_michelMuon_X_Y","Michel Muon End Point X vs Y",   binList.vertex_x_y.get_nBins(), binList.vertex_x_y.get_min(), binList.vertex_x_y.get_max(),
    binList.vertex_x_y.get_nBins(), binList.vertex_x_y.get_min(), binList.vertex_x_y.get_max());
    michelMuon_X_Y[1]->GetXaxis()->SetTitle("Michel Muon End Point X [mm]");
    michelMuon_X_Y[1]->GetYaxis()->SetTitle("Michel Muon End Point [mm]");
    
    michelMuon_Z[1] = new TH1D( "track_michelMuon_Z","Michel Muon End Point Z",binList.vertex_z.get_nBins(), binList.vertex_z.get_min(), binList.vertex_z.get_max() );
    michelMuon_Z[1]->GetXaxis()->SetTitle("z = 4293 Target, #bf{z = 5810 Interaction Region}, z = 8614 ECAL, z = 9088 HCAL");
    michelMuon_Z[1]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.vertex_z.get_width()));
    
    michelMuon_Z_vtx[1] = new TH1D( "track_michelMuon_Z_vtx","Michel Muon End Point Z - Vertex Z",binList.michelMuon_Z_vtx.get_nBins(), binList.michelMuon_Z_vtx.get_min(), binList.michelMuon_Z_vtx.get_max() );
    michelMuon_Z_vtx[1]->GetXaxis()->SetTitle("Michel Muon End Point Z - Vertex Z");
    michelMuon_Z_vtx[1]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelMuon_Z_vtx.get_width()));
    
    michelPion_P[1] = new TH1D( "track_michelPion_P","Michel Pion Momentum",binList.michelPion_P.get_nBins(), binList.michelPion_P.get_min(), binList.michelPion_P.get_max() );
    michelPion_P[1]->GetXaxis()->SetTitle("Michel Pion Momentum [MeV]");
    michelPion_P[1]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelPion_P.get_width()));
    
    michelPion_begin_dist_vtx[1] = new TH1D( "track_michelPion_begin_dist_vtx","Michel Pion Inital Point Distance to Vertex",binList.michelPion_begin_dist_vtx.get_nBins(), binList.michelPion_begin_dist_vtx.get_min(), binList.michelPion_begin_dist_vtx.get_max() );
    michelPion_begin_dist_vtx[1]->GetXaxis()->SetTitle("Michel Pion Inital Point Distance to Vertex [mm]");
    michelPion_begin_dist_vtx[1]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelPion_begin_dist_vtx.get_width()));

    michelPion_length[1] = new TH1D( "track_michelPion_length","Michel Pion Track Length",binList.michelPion_length.get_nBins(), binList.michelPion_length.get_min(), binList.michelPion_length.get_max() );
    michelPion_length[1]->GetXaxis()->SetTitle("Michel Pion Track Length [mm]");
    michelPion_length[1]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelPion_length.get_width()));
    
    michelMuon_dist_michelPion_length[1] = new TH2D( "track_michelMuon_dist_michelPion_length","Michel Pion Length vs Michel Muon Distance to Vertex",binList.michelMuon_end_dist_vtx.get_nBins(), binList.michelMuon_end_dist_vtx.get_min(), binList.michelMuon_end_dist_vtx.get_max(),binList.michelPion_begin_dist_vtx.get_nBins(), binList.michelPion_begin_dist_vtx.get_min(), binList.michelPion_begin_dist_vtx.get_max());
    michelMuon_dist_michelPion_length[1]->GetXaxis()->SetTitle("Michel Muon End Point Distance to Vertex [mm]");
    michelMuon_dist_michelPion_length[1]->GetYaxis()->SetTitle("Michel Pion Length [mm]");
    
    // Michel Study - Found Events by Secondary Track End Point Michel
    michelElectron_E[2] = new TH1D( "track2_michelElectron_E","Michel Electron Energy",binList.michelMuon_P.get_nBins(), binList.michelMuon_P.get_min(), binList.michelMuon_P.get_max() );
    michelElectron_E[2]->GetXaxis()->SetTitle("Michel Electron Energy [MeV]");
    michelElectron_E[2]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelMuon_P.get_width()));
    
    michelMuon_P[2] = new TH1D( "track2_michelMuon_P","Michel Muon Momentum",binList.michelMuon_P.get_nBins(), binList.michelMuon_P.get_min(), binList.michelMuon_P.get_max() );
    michelMuon_P[2]->GetXaxis()->SetTitle("Michel Muon Momentum [MeV]");
    michelMuon_P[2]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelMuon_P.get_width()));
    
    michelMuon_end_dist_vtx[2] = new TH1D( "track2_michelMuon_end_dist_vtx","Michel Muon End Point Distance to Vertex",binList.michelMuon_end_dist_vtx.get_nBins(), binList.michelMuon_end_dist_vtx.get_min(), binList.michelMuon_end_dist_vtx.get_max() );
    michelMuon_end_dist_vtx[2]->GetXaxis()->SetTitle("Michel Muon End Point Distance to Vertex [mm]");
    michelMuon_end_dist_vtx[2]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelMuon_end_dist_vtx.get_width()));

    michelMuon_length[2] = new TH1D( "track2_michelMuon_length","Michel Muon Track Length",binList.michelMuon_length.get_nBins(), binList.michelMuon_length.get_min(), binList.michelMuon_length.get_max() );
    michelMuon_length[2]->GetXaxis()->SetTitle("Michel Muon Track Length [mm]");
    michelMuon_length[2]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelMuon_length.get_width()));
    
    michelMuon_X_Y[2] = new TH2D( "track2_michelMuon_X_Y","Michel Muon End Point X vs Y",   binList.vertex_x_y.get_nBins(), binList.vertex_x_y.get_min(), binList.vertex_x_y.get_max(),
    binList.vertex_x_y.get_nBins(), binList.vertex_x_y.get_min(), binList.vertex_x_y.get_max());
    michelMuon_X_Y[2]->GetXaxis()->SetTitle("Michel Muon End Point X [mm]");
    michelMuon_X_Y[2]->GetYaxis()->SetTitle("Michel Muon End Point [mm]");
    
    michelMuon_Z[2] = new TH1D( "track2_michelMuon_Z","Michel Muon End Point Z",binList.vertex_z.get_nBins(), binList.vertex_z.get_min(), binList.vertex_z.get_max() );
    michelMuon_Z[2]->GetXaxis()->SetTitle("z = 4293 Target, #bf{z = 5810 Interaction Region}, z = 8614 ECAL, z = 9088 HCAL");
    michelMuon_Z[2]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.vertex_z.get_width()));
    
    michelMuon_Z_vtx[2] = new TH1D( "track2_michelMuon_Z_vtx","Michel Muon End Point Z - Vertex Z",binList.michelMuon_Z_vtx.get_nBins(), binList.michelMuon_Z_vtx.get_min(), binList.michelMuon_Z_vtx.get_max() );
    michelMuon_Z_vtx[2]->GetXaxis()->SetTitle("Michel Muon End Point Z - Vertex Z");
    michelMuon_Z_vtx[2]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelMuon_Z_vtx.get_width()));
    
    michelPion_P[2] = new TH1D( "track2_michelPion_P","Michel Pion Momentum",binList.michelPion_P.get_nBins(), binList.michelPion_P.get_min(), binList.michelPion_P.get_max() );
    michelPion_P[2]->GetXaxis()->SetTitle("Michel Pion Momentum [MeV]");
    michelPion_P[2]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelPion_P.get_width()));
    
    michelPion_begin_dist_vtx[2] = new TH1D( "track2_michelPion_begin_dist_vtx","Michel Pion Inital Point Distance to Vertex",binList.michelPion_begin_dist_vtx.get_nBins(), binList.michelPion_begin_dist_vtx.get_min(), binList.michelPion_begin_dist_vtx.get_max() );
    michelPion_begin_dist_vtx[2]->GetXaxis()->SetTitle("Michel Pion Inital Point Distance to Vertex [mm]");
    michelPion_begin_dist_vtx[2]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelPion_begin_dist_vtx.get_width()));

    michelPion_length[2] = new TH1D( "track2_michelPion_length","Michel Pion Track Length",binList.michelPion_length.get_nBins(), binList.michelPion_length.get_min(), binList.michelPion_length.get_max() );
    michelPion_length[2]->GetXaxis()->SetTitle("Michel Pion Track Length [mm]");
    michelPion_length[2]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelPion_length.get_width()));
    
    michelMuon_dist_michelPion_length[2] = new TH2D( "track2_michelMuon_dist_michelPion_length","Michel Pion Length vs Michel Muon Distance to Vertex",binList.michelMuon_end_dist_vtx.get_nBins(), binList.michelMuon_end_dist_vtx.get_min(), binList.michelMuon_end_dist_vtx.get_max(),binList.michelPion_begin_dist_vtx.get_nBins(), binList.michelPion_begin_dist_vtx.get_min(), binList.michelPion_begin_dist_vtx.get_max());
    michelMuon_dist_michelPion_length[2]->GetXaxis()->SetTitle("Michel Muon End Point Distance to Vertex [mm]");
    michelMuon_dist_michelPion_length[2]->GetYaxis()->SetTitle("Michel Pion Length [mm]");
    
    // Michel Study - Missed Events
    michelElectron_E[3] = new TH1D( "missed_michelElectron_E","Michel Electron Energy",binList.michelMuon_P.get_nBins(), binList.michelMuon_P.get_min(), binList.michelMuon_P.get_max() );
    michelElectron_E[3]->GetXaxis()->SetTitle("Michel Electron Energy [MeV]");
    michelElectron_E[3]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelMuon_P.get_width()));
    
    michelMuon_P[3] = new TH1D( "missed_michelMuon_P","Michel Muon Momentum",binList.michelMuon_P.get_nBins(), binList.michelMuon_P.get_min(), binList.michelMuon_P.get_max() );
    michelMuon_P[3]->GetXaxis()->SetTitle("Michel Muon Momentum [MeV]");
    michelMuon_P[3]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelMuon_P.get_width()));
    
    michelMuon_end_dist_vtx[3] = new TH1D( "missed_michelMuon_end_dist_vtx","Michel Muon End Point Distance to Vertex",binList.michelMuon_end_dist_vtx.get_nBins(), binList.michelMuon_end_dist_vtx.get_min(), binList.michelMuon_end_dist_vtx.get_max() );
    michelMuon_end_dist_vtx[3]->GetXaxis()->SetTitle("Michel Muon End Point Distance to Vertex [mm]");
    michelMuon_end_dist_vtx[3]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelMuon_end_dist_vtx.get_width()));

    michelMuon_length[3] = new TH1D( "missed_michelMuon_length","Michel Muon Track Length",binList.michelMuon_length.get_nBins(), binList.michelMuon_length.get_min(), binList.michelMuon_length.get_max() );
    michelMuon_length[3]->GetXaxis()->SetTitle("Michel Muon Track Length [mm]");
    michelMuon_length[3]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelMuon_length.get_width()));
    
    michelMuon_X_Y[3] = new TH2D( "missed_michelMuon_X_Y","Michel Muon End Point X vs Y",   binList.vertex_x_y.get_nBins(), binList.vertex_x_y.get_min(), binList.vertex_x_y.get_max(),
    binList.vertex_x_y.get_nBins(), binList.vertex_x_y.get_min(), binList.vertex_x_y.get_max());
    michelMuon_X_Y[3]->GetXaxis()->SetTitle("Michel Muon End Point X [mm]");
    michelMuon_X_Y[3]->GetYaxis()->SetTitle("Michel Muon End Point [mm]");
    
    michelMuon_Z[3] = new TH1D( "missed_michelMuon_Z","Michel Muon End Point Z",binList.vertex_z.get_nBins(), binList.vertex_z.get_min(), binList.vertex_z.get_max() );
    michelMuon_Z[3]->GetXaxis()->SetTitle("z = 4293 Target, #bf{z = 5810 Interaction Region}, z = 8614 ECAL, z = 9088 HCAL");
    michelMuon_Z[3]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.vertex_z.get_width()));
    
    michelMuon_Z_vtx[3] = new TH1D( "missed_michelMuon_Z_vtx","Michel Muon End Point Z - Vertex Z",binList.michelMuon_Z_vtx.get_nBins(), binList.michelMuon_Z_vtx.get_min(), binList.michelMuon_Z_vtx.get_max() );
    michelMuon_Z_vtx[3]->GetXaxis()->SetTitle("Michel Muon End Point Z - Vertex Z");
    michelMuon_Z_vtx[3]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelMuon_Z_vtx.get_width()));
    
    michelPion_P[3] = new TH1D( "missed_michelPion_P","Michel Pion Momentum",binList.michelPion_P.get_nBins(), binList.michelPion_P.get_min(), binList.michelPion_P.get_max() );
    michelPion_P[3]->GetXaxis()->SetTitle("Michel Pion Momentum [MeV]");
    michelPion_P[3]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelPion_P.get_width()));
    
    michelPion_begin_dist_vtx[3] = new TH1D( "missed_michelPion_begin_dist_vtx","Michel Pion Inital Point Distance to Vertex",binList.michelPion_begin_dist_vtx.get_nBins(), binList.michelPion_begin_dist_vtx.get_min(), binList.michelPion_begin_dist_vtx.get_max() );
    michelPion_begin_dist_vtx[3]->GetXaxis()->SetTitle("Michel Pion Inital Point Distance to Vertex [mm]");
    michelPion_begin_dist_vtx[3]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelPion_begin_dist_vtx.get_width()));

    michelPion_length[3] = new TH1D( "missed_michelPion_length","Michel Pion Track Length",binList.michelPion_length.get_nBins(), binList.michelPion_length.get_min(), binList.michelPion_length.get_max() );
    michelPion_length[3]->GetXaxis()->SetTitle("Michel Pion Track Length [mm]");
    michelPion_length[3]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelPion_length.get_width()));
    
    michelMuon_dist_michelPion_length[3] = new TH2D( "missed_michelMuon_dist_michelPion_length","Michel Pion Length vs Michel Muon Distance to Vertex",binList.michelMuon_end_dist_vtx.get_nBins(), binList.michelMuon_end_dist_vtx.get_min(), binList.michelMuon_end_dist_vtx.get_max(),binList.michelPion_begin_dist_vtx.get_nBins(), binList.michelPion_begin_dist_vtx.get_min(), binList.michelPion_begin_dist_vtx.get_max());
    michelMuon_dist_michelPion_length[3]->GetXaxis()->SetTitle("Michel Muon End Point Distance to Vertex [mm]");
    michelMuon_dist_michelPion_length[3]->GetYaxis()->SetTitle("Michel Pion Length [mm]");
    
    // All Found Michels
    michelElectron_E[4] = new TH1D( "found_michelElectron_E","Michel Electron Energy",binList.michelMuon_P.get_nBins(), binList.michelMuon_P.get_min(), binList.michelMuon_P.get_max() );
    michelElectron_E[4]->GetXaxis()->SetTitle("Michel Electron Energy [MeV]");
    michelElectron_E[4]->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.michelMuon_P.get_width()));
}

void CCProtonPi0_MichelTool::write_RootFile()
{
    cout<<">> Writing "<<rootDir<<endl;
    f->Write();    
}





#endif
