#define ANA_PIDStudies_cxx
#include "ANA_PIDStudies.h"

#include <TH1.h>
#include <TH2.h>
#include <TMath.h>

ANA_PIDStudies::ANA_PIDStudies
{
    // Do Nothing For Now
}

void ANA_PIDStudies::run(string playlist, char* filename, string cutFile, string readmeFile)
{

    //------------------------------------------------------------------------
    // Open files for writing
    //------------------------------------------------------------------------

    openFiles(cutFile,readmeFile);

    TFile* f = new TFile(filename,"RECREATE");

    //------------------------------------------------------------------------
    // Create chain
    //------------------------------------------------------------------------

    TChain* fChain = new TChain("pidStudies");
    Init(playlist, fChain);

    if (!fChain) return;
    if (fChain == 0) return;

    //------------------------------------------------------------------------
    // Select only required Branches for Performance
    //------------------------------------------------------------------------
//     fChain->SetBranchStatus("*",0);  // disable all branches

    //------------------------------------------------------------------------
    //  1D Histograms
    //------------------------------------------------------------------------

    // Measured dEdX Pions
    TH1F* dEdX_pion = new TH1F("measured_pion_dEdX","Measured Pion dEdX",NBINS_dEdX , MIN_dEdX , MAX_dEdX );
    dEdX_pion->GetXaxis()->SetTitle("Measured dEdX for Pions");
    dEdX_pion->GetYaxis()->SetTitle( Form("Candidates / %4.2f ",WIDTH_dEdX) );
    
    TH1F* dEdX_difference = new TH1F("dEdX_difference","dEdX Difference",NBINS_delta_dEdX , MIN_delta_dEdX , MAX_delta_dEdX );
    dEdX_difference->GetXaxis()->SetTitle("dEdX Difference");
    dEdX_difference->GetYaxis()->SetTitle( Form("Candidates / %4.2f ",WIDTH_delta_dEdX) );
    
    TH1F* pion_score = new TH1F("pion_score","Pion Score",NBINS_pID , MIN_pID , MAX_pID);
    pion_score->GetXaxis()->SetTitle("Pion Score");
    pion_score->GetYaxis()->SetTitle( Form("Candidates / %4.2f ",WIDTH_pID) );
    
    TH1F* proton_score = new TH1F("proton_score","Proton Score",NBINS_pID , MIN_pID , MAX_pID);
    proton_score->GetXaxis()->SetTitle("Proton Score");
    proton_score->GetYaxis()->SetTitle( Form("Candidates / %4.2f ",WIDTH_pID) );
    
    TH2F* dEdX_Nodes = new TH2F("dEdX_Nodes","dEdX vs Node #",NBINS_Nodes , MIN_Nodes , MAX_Nodes , NBINS_dEdX , MIN_dEdX , MAX_dEdX );
    dEdX_Nodes->GetXaxis()->SetTitle("Node #");
    dEdX_Nodes->GetYaxis()->SetTitle("dEdX");
    
    TH2F* pion_score_proton_mom = new TH2F("pion_score_proton_mom","Proton Momentum vs Pion Score",NBINS_pID , MIN_pID , MAX_pID, NBINS_Pp , MIN_Pp , MAX_Pp );
    pion_score_proton_mom->GetXaxis()->SetTitle("Pion Score");
    pion_score_proton_mom->GetYaxis()->SetTitle("True Proton Momentum [MeV]");
    

    //------------------------------------------------------------------------
    // Loop over Chain
    //------------------------------------------------------------------------
        
    Long64_t nentries = fChain->GetEntriesFast();
    cout<<"There are "<<nentries<<" entries!"<<endl;
    
    double P;
    
    double nAll, nLow, nMiddle, nHigh;
    
    nAll = 0; nLow = 0; nMiddle = 0; nHigh = 0;

    Long64_t nbytes = 0, nb = 0;
    for (Long64_t jentry=0; jentry<nentries;jentry++) {

        Long64_t ientry = fChain->GetEntry(jentry);
        if (ientry == 0) {
            cout<<"GetEntry failure "<<jentry<<endl;
            break;
        }
        nb = fChain->GetEntry(jentry);   nbytes += nb;    

        if (jentry%5000 == 0){
            cout<<" Entry "<<jentry<<endl;
        }
        

        for( int i = 0; (i < n_prongs) && (i < 3); i++){
            pion_score->Fill(dedx_pion_score1[i]);
            proton_score->Fill(dedx_proton_score1[i]);
            P = calcMomentum(part_4p[0][1],part_4p[0][2],part_4p[0][3]);
            
            // Do not count pion_score = -1
            if (dedx_pion_score1[i] >= 0.0){
                nAll++;
                if(dedx_pion_score1[i] <= 0.2){
                    nLow++;
                }else if(dedx_pion_score1[i] < 0.7){
                    nMiddle++;
                }else{
                    nHigh++;
                }
            }
            
            pion_score_proton_mom->Fill(dedx_pion_score1[i], P);
        }
    
    }
    
    //Write on Cut File
    cutText<<"All Events            = "<<getPercent(nAll,nAll)<<endl;
    cutText<<"PionScore Low         = "<<getPercent(nAll,nLow)<<endl;
    cutText<<"PionScore Middle      = "<<getPercent(nAll,nMiddle)<<endl;
    cutText<<"PionScore High        = "<<getPercent(nAll,nHigh)<<endl;

    closeFiles();

    f->Write();

}

double ANA_PIDStudies::calcMomentum(double px, double py, double pz)
{
    double P;
    
    P = sqrt( px*px + py*py + pz*pz);
    
    return P;

}


void ANA_PIDStudies::print()
{
    cout<<ev_run<<" "<<ev_sub_run<<" "<<ev_gate<<endl;
}

double ANA_PIDStudies::getPercent(double nAll, double nCurrent)
{
    double percent;

    percent = (nCurrent / nAll) * 100;

    return percent;
}


void ANA_PIDStudies::closeFiles()
{
    cutText.close();
    readme.close();
}

void ANA_PIDStudies::openFiles(string cutFile, string readmeFile)
{
    cutText.open( cutFile.c_str() );
    if( !cutText.is_open() ){
        cerr<<"Cannot open output text file: "<<cutFile<<endl;
        exit(1);
    }

    readme.open( readmeFile.c_str() );
    if( !readme.is_open() ){
        cerr<<"Cannot open output text file: "<<readmeFile<<endl;
        exit(1);
    }

//     writeReadme();


}

void ANA_PIDStudies::writeReadme()
{
  
        
}

// -------------------------------------------------------------------------
//     Default Functions
//--------------------------------------------------------------------------

#ifdef ANA_PIDStudies_cxx
ANA_PIDStudies::ANA_PIDStudies()
{
    
}

ANA_PIDStudies::~ANA_PIDStudies()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t ANA_PIDStudies::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t ANA_PIDStudies::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void ANA_PIDStudies::Init(string playlist, TChain* fChain)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   ifstream input_pl(playlist.c_str());
   string filename;

    if( !input_pl.is_open() ){
        cerr<<"Cannot open Playlist File!"<<endl;
        exit(1);
    }else{
        cout<<"Playlist: "<<playlist.c_str()<<endl;
    }


   while (input_pl) {
     input_pl>>filename;
     
     if (!input_pl) break;
    
     if (filename[0] != '/') break;
    
     fChain->Add( filename.c_str() );
//      cout<<" Added "<<filename.c_str()<<endl;
   }

   // Set branch addresses and branch pointers
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("subrun", &subrun, &b_subrun);
   fChain->SetBranchAddress("gate", &gate, &b_gate);
   fChain->SetBranchAddress("interaction", &interaction, &b_interaction);
   fChain->SetBranchAddress("n_particles", &n_particles, &b_n_particles);
   fChain->SetBranchAddress("part_pdg", part_pdg, &b_part_pdg);
   fChain->SetBranchAddress("part_trackId", part_trackId, &b_part_trackId);
   fChain->SetBranchAddress("part_fplane", part_fplane, &b_part_fplane);
   fChain->SetBranchAddress("part_lplane", part_lplane, &b_part_lplane);
   fChain->SetBranchAddress("part_destructCode", part_destructCode, &b_part_destructCode);
   fChain->SetBranchAddress("part_beginPoint", part_beginPoint, &b_part_beginPoint);
   fChain->SetBranchAddress("part_endPoint", part_endPoint, &b_part_endPoint);
   fChain->SetBranchAddress("part_theta", part_theta, &b_part_theta);
   fChain->SetBranchAddress("part_phi", part_phi, &b_part_phi);
   fChain->SetBranchAddress("part_4p", part_4p, &b_part_4p);
   fChain->SetBranchAddress("part_endTheta", part_endTheta, &b_part_endTheta);
   fChain->SetBranchAddress("part_endPhi", part_endPhi, &b_part_endPhi);
   fChain->SetBranchAddress("part_end4p", part_end4p, &b_part_end4p);
   fChain->SetBranchAddress("part_lastPlaneDistanceZ", part_lastPlaneDistanceZ, &b_part_lastPlaneDistanceZ);
   fChain->SetBranchAddress("n_prongs", &n_prongs, &b_n_prongs);
   fChain->SetBranchAddress("prong_isForked", prong_isForked, &b_prong_isForked);
   fChain->SetBranchAddress("prong_isKinked", prong_isKinked, &b_prong_isKinked);
   fChain->SetBranchAddress("prong_isDSEcal", prong_isDSEcal, &b_prong_isDSEcal);
   fChain->SetBranchAddress("prong_isSideEcal", prong_isSideEcal, &b_prong_isSideEcal);
   fChain->SetBranchAddress("prong_isNuclTargs", prong_isNuclTargs, &b_prong_isNuclTargs);
   fChain->SetBranchAddress("prong_priTrkPatRec", prong_priTrkPatRec, &b_prong_priTrkPatRec);
   fChain->SetBranchAddress("prong_fplane", prong_fplane, &b_prong_fplane);
   fChain->SetBranchAddress("prong_lplane", prong_lplane, &b_prong_lplane);
   fChain->SetBranchAddress("prong_hasSplitClusters", prong_hasSplitClusters, &b_prong_hasSplitClusters);
   fChain->SetBranchAddress("prong_hasEndSplitClusters", prong_hasEndSplitClusters, &b_prong_hasEndSplitClusters);
   fChain->SetBranchAddress("prong_beginPoint", prong_beginPoint, &b_prong_beginPoint);
   fChain->SetBranchAddress("prong_endPoint", prong_endPoint, &b_prong_endPoint);
   fChain->SetBranchAddress("prong_beginTheta", prong_beginTheta, &b_prong_beginTheta);
   fChain->SetBranchAddress("prong_beginPhi", prong_beginPhi, &b_prong_beginPhi);
   fChain->SetBranchAddress("prong_endTheta", prong_endTheta, &b_prong_endTheta);
   fChain->SetBranchAddress("prong_endPhi", prong_endPhi, &b_prong_endPhi);
   fChain->SetBranchAddress("prong_endm1Theta", prong_endm1Theta, &b_prong_endm1Theta);
   fChain->SetBranchAddress("prong_endm1Phi", prong_endm1Phi, &b_prong_endm1Phi);
   fChain->SetBranchAddress("tm_partProngMatchPDG", tm_partProngMatchPDG, &b_tm_partProngMatchPDG);
   fChain->SetBranchAddress("tm_partProngMatchTrackId", tm_partProngMatchTrackId, &b_tm_partProngMatchTrackId);
   fChain->SetBranchAddress("tm_partTrackMatchPDG", tm_partTrackMatchPDG, &b_tm_partTrackMatchPDG);
   fChain->SetBranchAddress("tm_partTrackMatchTrackId", tm_partTrackMatchTrackId, &b_tm_partTrackMatchTrackId);
   fChain->SetBranchAddress("tm_EFraction", tm_EFraction, &b_tm_EFraction);
   fChain->SetBranchAddress("dedx_pion_p", dedx_pion_p, &b_dedx_pion_p);
   fChain->SetBranchAddress("dedx_pion_fitpmin", dedx_pion_fitpmin, &b_dedx_pion_fitpmin);
   fChain->SetBranchAddress("dedx_pion_fitpmax", dedx_pion_fitpmax, &b_dedx_pion_fitpmax);
   fChain->SetBranchAddress("dedx_pion_score1", dedx_pion_score1, &b_dedx_pion_score1);
   fChain->SetBranchAddress("dedx_pion_score2", dedx_pion_score2, &b_dedx_pion_score2);
   fChain->SetBranchAddress("dedx_pion_fit_by_range", dedx_pion_fit_by_range, &b_dedx_pion_fit_by_range);
   fChain->SetBranchAddress("dedx_pion_fit_fails", dedx_pion_fit_fails, &b_dedx_pion_fit_fails);
   fChain->SetBranchAddress("dedx_pion_nnodes", dedx_pion_nnodes, &b_dedx_pion_nnodes);
   fChain->SetBranchAddress("dedx_pion_chi2", dedx_pion_chi2, &b_dedx_pion_chi2);
   fChain->SetBranchAddress("dedx_pion_range_diff", dedx_pion_range_diff, &b_dedx_pion_range_diff);
   fChain->SetBranchAddress("dedx_pion_measuredDEDX", dedx_pion_measuredDEDX, &b_dedx_pion_measuredDEDX);
   fChain->SetBranchAddress("dedx_pion_calculatedDEDX", dedx_pion_calculatedDEDX, &b_dedx_pion_calculatedDEDX);
   fChain->SetBranchAddress("dedx_pion_zpositions", dedx_pion_zpositions, &b_dedx_pion_zpositions);
   fChain->SetBranchAddress("dedx_proton_p", dedx_proton_p, &b_dedx_proton_p);
   fChain->SetBranchAddress("dedx_proton_score1", dedx_proton_score1, &b_dedx_proton_score1);
   fChain->SetBranchAddress("dedx_proton_score2", dedx_proton_score2, &b_dedx_proton_score2);
   fChain->SetBranchAddress("dedx_proton_fit_by_range", dedx_proton_fit_by_range, &b_dedx_proton_fit_by_range);
   fChain->SetBranchAddress("dedx_proton_fit_fails", dedx_proton_fit_fails, &b_dedx_proton_fit_fails);
   fChain->SetBranchAddress("dedx_proton_nnodes", dedx_proton_nnodes, &b_dedx_proton_nnodes);
   fChain->SetBranchAddress("dedx_proton_chi2", dedx_proton_chi2, &b_dedx_proton_chi2);
   fChain->SetBranchAddress("dedx_proton_range_diff", dedx_proton_range_diff, &b_dedx_proton_range_diff);
   fChain->SetBranchAddress("dedx_proton_measuredDEDX", dedx_proton_measuredDEDX, &b_dedx_proton_measuredDEDX);
   fChain->SetBranchAddress("dedx_proton_calculatedDEDX", dedx_proton_calculatedDEDX, &b_dedx_proton_calculatedDEDX);
   fChain->SetBranchAddress("dedx_proton_zpositions", dedx_proton_zpositions, &b_dedx_proton_zpositions);
   fChain->SetBranchAddress("spid_proton_score", spid_proton_score, &b_spid_proton_score);
   fChain->SetBranchAddress("spid_piminus_score", spid_piminus_score, &b_spid_piminus_score);
   fChain->SetBranchAddress("spid_piplus_score", spid_piplus_score, &b_spid_piplus_score);
   Notify();
}

Bool_t ANA_PIDStudies::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

    cout<<"Initialization is OK!"<<endl;
   return kTRUE;
}

void ANA_PIDStudies::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

Int_t ANA_PIDStudies::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef ANA_PIDStudies_cxx

