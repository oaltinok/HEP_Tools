/*
    See PIDTool.h header for Class Information
*/
#ifndef PIDTool_cpp
#define PIDTool_cpp

#include "PIDTool.h"

using namespace std;

void PIDTool::get_pID_Stats()
{
    get_pID_Stats_LLR();
    get_pID_Stats_pIDDiff();
}

void PIDTool::get_pID_Stats_LLR()
{
    cout<<"=== Calculating pID Statistics for LLR ==="<<endl;

    double nProton = 0;
    double nTotalProton = 0;
    double nCapturedEvents = 0;
    double nEvents = 0;
    double purity;
    double efficiency;
    int nBins = binList.particleScore_LLR.get_nBins();
    
    //Get Total Proton
    for(int i = nBins; i >= 1; i--){
        nTotalProton = nTotalProton + proton_protonScore_LLR->GetBinContent(i);
    }
    cout<<"nTotalProton = "<<nTotalProton<<endl;
    
    
    for(int i = nBins; i >= 1; i--){
        nProton = nProton + proton_protonScore_LLR->GetBinContent(i);
        
        nCapturedEvents =   nCapturedEvents +
                            proton_protonScore_LLR->GetBinContent(i) +
                            piplus_protonScore_LLR->GetBinContent(i) +
                            piminus_protonScore_LLR->GetBinContent(i) +
                            other_protonScore_LLR->GetBinContent(i);
                            
        nEvents =           proton_protonScore_LLR->GetBinContent(i) +
                            piplus_protonScore_LLR->GetBinContent(i) +
                            piminus_protonScore_LLR->GetBinContent(i) +
                            other_protonScore_LLR->GetBinContent(i);
                            
        if (nCapturedEvents == 0 ) purity = 0;
        else purity = nProton / nCapturedEvents;
        
        efficiency = nProton / nTotalProton;
        purity_LLR->SetBinContent(i,purity);
        efficiency_LLR->SetBinContent(i,efficiency);
        purityXefficiency_LLR->SetBinContent(i,purity*efficiency);
//         cout<<"pID = "<<proton_protonScore_LLR->GetBinLowEdge(i)<<" Purity = "<<purity<<" Efficiency = "<<efficiency<<endl;
    }
    
    double maxPurityXEfficiency = purityXefficiency_LLR->GetMaximum();
    int bin_max = purityXefficiency_LLR->GetMaximumBin();
    double pIDDiff = purityXefficiency_LLR->GetBinLowEdge(bin_max);
    cout<<"Maximum = "<<maxPurityXEfficiency<<" Bin = "<<pIDDiff<<endl;
}

void PIDTool::get_pID_Stats_pIDDiff()
{
    cout<<"=== Calculating pID Statistics for pID Diff ==="<<endl;

    double nProton = 0;
    double nTotalProton = 0;
    double nCapturedEvents = 0;
    double nEvents = 0;
    double purity;
    double efficiency;
    int nBins = binList.particleScoreDiff.get_nBins();
    
    //Get Total Proton
    for(int i = nBins; i >= 1; i--){
        nTotalProton = nTotalProton + proton_pIDDiff->GetBinContent(i);
    }
    cout<<"nTotalProton = "<<nTotalProton<<endl;
    
    
    for(int i = nBins; i >= 1; i--){
        nProton = nProton + proton_pIDDiff->GetBinContent(i);
        
        nCapturedEvents =   nCapturedEvents +
                            proton_pIDDiff->GetBinContent(i) +
                            piplus_pIDDiff->GetBinContent(i) +
                            piminus_pIDDiff->GetBinContent(i) +
                            other_pIDDiff->GetBinContent(i);
                            
        nEvents =           proton_pIDDiff->GetBinContent(i) +
                            piplus_pIDDiff->GetBinContent(i) +
                            piminus_pIDDiff->GetBinContent(i) +
                            other_pIDDiff->GetBinContent(i);
                            
        if (nCapturedEvents == 0 ) purity = 0;
        else purity = nProton / nCapturedEvents;
        
        efficiency = nProton / nTotalProton;
        purity_pIDDiff->SetBinContent(i,purity);
        efficiency_pIDDiff->SetBinContent(i,efficiency);
        purityXefficiency_pIDDiff->SetBinContent(i,purity*efficiency);
//         cout<<"pID Diff = "<<proton_pIDDiff->GetBinLowEdge(i)<<" Purity = "<<purity<<" Efficiency = "<<efficiency<<endl;
    }
    
    double maxPurityXEfficiency = purityXefficiency_pIDDiff->GetMaximum();
    int bin_max = purityXefficiency_pIDDiff->GetMaximumBin();
    double pIDDiff = purityXefficiency_pIDDiff->GetBinLowEdge(bin_max);
    cout<<"Maximum = "<<maxPurityXEfficiency<<" Bin = "<<pIDDiff<<endl;
    
}



void PIDTool::FillHistograms(   double protonScore_LLR, double protonScore, double pionScore,
                                int truthPDG, double prongE)
{
    if(truthPDG == PDG_List::proton){
        proton_protonScore_LLR->Fill(protonScore_LLR);
        proton_protonScore->Fill(protonScore);
        proton_pionScore_protonScore->Fill(pionScore,protonScore);
        proton_pIDDiff->Fill(protonScore - pionScore);
        proton_protonScore_protonScore_LLR->Fill(protonScore,protonScore_LLR);
        KE_proton_pIDDiff->Fill(prongE - 938.27);
        KE_proton_LLR->Fill(prongE - 938.27);
    }else if(truthPDG == PDG_List::pi_plus){
        piplus_protonScore_LLR->Fill(protonScore_LLR);
        piplus_protonScore->Fill(protonScore);
        piplus_pionScore_protonScore->Fill(pionScore,protonScore);
        piplus_pIDDiff->Fill(protonScore - pionScore);
        piplus_protonScore_protonScore_LLR->Fill(protonScore,protonScore_LLR);
        KE_other_pIDDiff->Fill(prongE - 938.27);
        KE_other_LLR->Fill(prongE - 938.27);
    }else if(truthPDG == PDG_List::pi_minus){
        piminus_protonScore_LLR->Fill(protonScore_LLR);
        piminus_protonScore->Fill(protonScore);
        piminus_pionScore_protonScore->Fill(pionScore,protonScore);
        piminus_pIDDiff->Fill(protonScore - pionScore);
        piminus_protonScore_protonScore_LLR->Fill(protonScore,protonScore_LLR);
        KE_other_pIDDiff->Fill(prongE - 938.27);
        KE_other_LLR->Fill(prongE - 938.27);
    }else{
        other_protonScore_LLR->Fill(protonScore_LLR);
        other_pIDDiff->Fill(protonScore - pionScore);
        KE_other_pIDDiff->Fill(prongE - 938.27);
        KE_other_LLR->Fill(prongE - 938.27);
    }
    
}


PIDTool::PIDTool()
{
    cout<<"Wrong usage of PIDTool! Must include Analysis Mode"<<endl;
    exit(EXIT_FAILURE);
}

PIDTool::PIDTool(int nMode)
{
    cout<<"Initializing PIDTool"<<endl;
    
    SetAnalysisMode(nMode);
    
    // File Locations
    rootDir =   Folder_List::output + Folder_List::rootOut + branchDir + "PIDStatistics.root";

    cout<<"\tRoot File: "<<rootDir<<endl;
    
    // Create Root File 
    f = new TFile(rootDir.c_str(),"RECREATE");
    
    // Initialize Histograms
    purity_LLR = new TH1D( "purity_LLR","Proton Purity",binList.particleScore_LLR.get_nBins(), binList.particleScore_LLR.get_min(), binList.particleScore_LLR.get_max()  );
    purity_LLR->GetXaxis()->SetTitle("Proton Purity = Captured Proton / Captured Total Events");
    purity_LLR->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScore.get_width()));
    
    efficiency_LLR = new TH1D( "efficiency_LLR","Proton Efficiency",binList.particleScore_LLR.get_nBins(), binList.particleScore_LLR.get_min(), binList.particleScore_LLR.get_max() );
    efficiency_LLR->GetXaxis()->SetTitle("Proton Efficiency = Captured Proton / Total Protons");
    efficiency_LLR->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScore.get_width()));
    
    purityXefficiency_LLR = new TH1D( "purityXefficiency_LLR","Proton Purity X Efficiency",binList.particleScore_LLR.get_nBins(), binList.particleScore_LLR.get_min(), binList.particleScore_LLR.get_max()  );
    purityXefficiency_LLR->GetXaxis()->SetTitle("Purity X Efficiency");
    purityXefficiency_LLR->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScore.get_width()));
    
    purity_pIDDiff = new TH1D( "purity_pIDDiff","Proton Purity using pID Diff",binList.particleScoreDiff.get_nBins(), binList.particleScoreDiff.get_min(), binList.particleScoreDiff.get_max()  );
    purity_pIDDiff->GetXaxis()->SetTitle("Proton Purity = Captured Proton / Captured Total Events");
    purity_pIDDiff->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScoreDiff.get_width()));
    
    efficiency_pIDDiff = new TH1D( "efficiency_pIDDiff","Proton Efficiency using pID Diff",binList.particleScoreDiff.get_nBins(), binList.particleScoreDiff.get_min(), binList.particleScoreDiff.get_max() );
    efficiency_pIDDiff->GetXaxis()->SetTitle("Proton Efficiency = Captured Proton / Total Protons");
    efficiency_pIDDiff->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScoreDiff.get_width()));
    
    purityXefficiency_pIDDiff = new TH1D( "purityXefficiency_pIDDiff","Proton Purity X Efficiency using pID Diff",binList.particleScoreDiff.get_nBins(), binList.particleScoreDiff.get_min(), binList.particleScoreDiff.get_max()  );
    purityXefficiency_pIDDiff->GetXaxis()->SetTitle("Purity X Efficiency");
    purityXefficiency_pIDDiff->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScoreDiff.get_width()));
    
    // Proton Score LLR
    piplus_protonScore_LLR = new TH1D( "piplus_protonScore_LLR","piplus_protonScore_LLR",binList.particleScore_LLR.get_nBins(), binList.particleScore_LLR.get_min(), binList.particleScore_LLR.get_max() );
    piplus_protonScore_LLR->GetXaxis()->SetTitle("piplus_protonScore_LLR");
    piplus_protonScore_LLR->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScore_LLR.get_width()));
    
    piminus_protonScore_LLR = new TH1D( "piminus_protonScore_LLR","piminus_protonScore_LLR",binList.particleScore_LLR.get_nBins(), binList.particleScore_LLR.get_min(), binList.particleScore_LLR.get_max() );
    piminus_protonScore_LLR->GetXaxis()->SetTitle("piminus_protonScore_LLR");
    piminus_protonScore_LLR->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScore_LLR.get_width()));
    
    proton_protonScore_LLR = new TH1D( "proton_protonScore_LLR","proton_protonScore_LLR",binList.particleScore_LLR.get_nBins(), binList.particleScore_LLR.get_min(), binList.particleScore_LLR.get_max() );
    proton_protonScore_LLR->GetXaxis()->SetTitle("proton_protonScore_LLR");
    proton_protonScore_LLR->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScore_LLR.get_width()));
    
    other_protonScore_LLR = new TH1D( "other_protonScore_LLR","other_protonScore_LLR",binList.particleScore_LLR.get_nBins(), binList.particleScore_LLR.get_min(), binList.particleScore_LLR.get_max() );
    other_protonScore_LLR->GetXaxis()->SetTitle("other_protonScore_LLR");
    other_protonScore_LLR->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScore_LLR.get_width()));
    
    // Proton Score 
    piplus_protonScore = new TH1D( "piplus_protonScore","piplus_protonScore",binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max() );
    piplus_protonScore->GetXaxis()->SetTitle("piplus_protonScore");
    piplus_protonScore->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScore.get_width()));
    
    piminus_protonScore = new TH1D( "piminus_protonScore","piminus_protonScore",binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max() );
    piminus_protonScore->GetXaxis()->SetTitle("piminus_protonScore");
    piminus_protonScore->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScore.get_width()));
    
    proton_protonScore = new TH1D( "proton_protonScore","proton_protonScore",binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max() );
    proton_protonScore->GetXaxis()->SetTitle("proton_protonScore");
    proton_protonScore->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScore.get_width()));
    
    
    
    // 2D Particle Score
    proton_pionScore_protonScore = new TH2D( "proton_pionScore_protonScore","proton_pionScore_protonScore",
                                                 binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max(),
                                                 binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max());
    proton_pionScore_protonScore->GetXaxis()->SetTitle("Pion Score");
    proton_pionScore_protonScore->GetYaxis()->SetTitle("Proton Score");
    
    piplus_pionScore_protonScore = new TH2D( "piplus_pionScore_protonScore","piplus_pionScore_protonScore",
                                                 binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max(),
                                                 binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max());
    piplus_pionScore_protonScore->GetXaxis()->SetTitle("Pion Score");
    piplus_pionScore_protonScore->GetYaxis()->SetTitle("Proton Score");
                                                 
    piminus_pionScore_protonScore = new TH2D( "piminus_pionScore_protonScore","piminus_pionScore_protonScore",
                                                 binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max(),
                                                 binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max());
    piminus_pionScore_protonScore->GetXaxis()->SetTitle("Pion Score");
    piminus_pionScore_protonScore->GetYaxis()->SetTitle("Proton Score");
    
    // 2D Proton Score vs ProtonScore_LLR
    proton_protonScore_protonScore_LLR = new TH2D( "proton_protonScore_protonScore_LLR","proton_protonScore_protonScore_LLR",
                                                 binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max(),
                                                 binList.particleScore_LLR.get_nBins(), binList.particleScore_LLR.get_min(), binList.particleScore_LLR.get_max());
    proton_protonScore_protonScore_LLR->GetXaxis()->SetTitle("Proton Score");
    proton_protonScore_protonScore_LLR->GetYaxis()->SetTitle("Proton Score LLR");
    
    piplus_protonScore_protonScore_LLR = new TH2D( "piplus_protonScore_protonScore_LLR","piplus_protonScore_protonScore_LLR",
                                                 binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max(),
                                                 binList.particleScore_LLR.get_nBins(), binList.particleScore_LLR.get_min(), binList.particleScore_LLR.get_max());
    piplus_protonScore_protonScore_LLR->GetXaxis()->SetTitle("Proton Score");
    piplus_protonScore_protonScore_LLR->GetYaxis()->SetTitle("Proton Score LLR");
                                                 
    piminus_protonScore_protonScore_LLR = new TH2D( "piminus_protonScore_protonScore_LLR","piminus_protonScore_protonScore_LLR",
                                                 binList.particleScore.get_nBins(), binList.particleScore.get_min(), binList.particleScore.get_max(),
                                                 binList.particleScore_LLR.get_nBins(), binList.particleScore_LLR.get_min(), binList.particleScore_LLR.get_max());
    piminus_protonScore_protonScore_LLR->GetXaxis()->SetTitle("Proton Score");
    piminus_protonScore_protonScore_LLR->GetYaxis()->SetTitle("Proton Score LLR");
    
    // Particle Score Difference
    
    proton_pIDDiff = new TH1D( "proton_pIDDiff","Proton Score - Pion Score for TRUE Protons",binList.particleScoreDiff.get_nBins(), binList.particleScoreDiff.get_min(), binList.particleScoreDiff.get_max() );
    proton_pIDDiff->GetXaxis()->SetTitle("Proton Score - Pion Score");
    proton_pIDDiff->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScoreDiff.get_width()));
    
    piplus_pIDDiff = new TH1D( "piplus_pIDDiff","Proton Score - Pion Score for TRUE Pi-Plus",binList.particleScoreDiff.get_nBins(), binList.particleScoreDiff.get_min(), binList.particleScoreDiff.get_max() );
    piplus_pIDDiff->GetXaxis()->SetTitle("Proton Score - Pion Score");
    piplus_pIDDiff->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScoreDiff.get_width()));
    
    piminus_pIDDiff = new TH1D( "piminus_pIDDiff","Proton Score - Pion Score for TRUE Pi-Minus",binList.particleScoreDiff.get_nBins(), binList.particleScoreDiff.get_min(), binList.particleScoreDiff.get_max() );
    piminus_pIDDiff->GetXaxis()->SetTitle("Proton Score - Pion Score");
    piminus_pIDDiff->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScoreDiff.get_width()));
    
    other_pIDDiff = new TH1D( "other_pIDDiff","Proton Score - Pion Score for TRUE Other",binList.particleScoreDiff.get_nBins(), binList.particleScoreDiff.get_min(), binList.particleScoreDiff.get_max() );
    other_pIDDiff->GetXaxis()->SetTitle("Proton Score - Pion Score");
    other_pIDDiff->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScoreDiff.get_width()));
    
    //--------------------------------------------------------------------------
    
    
    //! ------------------------------------------------------------------------
    //! Kinetic Energy
    //! ------------------------------------------------------------------------
    
    SingleBin bin_KE;
    bin_KE.setBin(20, 0.0, 2000.0);
    KE_proton_pIDDiff = new TH1D( "KE_proton_pIDDiff","Kinetic Energy of True Protons pID Diff",bin_KE.get_nBins(), bin_KE.get_min(), bin_KE.get_max() );
    KE_proton_pIDDiff->GetXaxis()->SetTitle("Proton Kinetic Energy [MeV]");
    KE_proton_pIDDiff->GetYaxis()->SetTitle(Form("Number of Protons / %3.1f [MeV] ",bin_KE.get_width()));
    
    KE_other_pIDDiff = new TH1D( "KE_other_pIDDiff","Kinetic Energy of Fake Protons pID Diff",bin_KE.get_nBins(), bin_KE.get_min(), bin_KE.get_max() );
    KE_other_pIDDiff->GetXaxis()->SetTitle("Fake Proton Kinetic Energy [MeV]");
    KE_other_pIDDiff->GetYaxis()->SetTitle(Form("Number of Protons / %3.1f [MeV] ",bin_KE.get_width()));
    
    KE_proton_LLR = new TH1D( "KE_proton_LLR","Kinetic Energy of True Protons LLR",bin_KE.get_nBins(), bin_KE.get_min(), bin_KE.get_max() );
    KE_proton_LLR->GetXaxis()->SetTitle("Proton Kinetic Energy [MeV]");
    KE_proton_LLR->GetYaxis()->SetTitle(Form("Number of Protons / %3.1f [MeV] ",bin_KE.get_width()));
    
    KE_other_LLR = new TH1D( "KE_other_LLR","Kinetic Energy of Fake Protons LLR",bin_KE.get_nBins(), bin_KE.get_min(), bin_KE.get_max() );
    KE_other_LLR->GetXaxis()->SetTitle("Fake Proton Kinetic Energy [MeV]");
    KE_other_LLR->GetYaxis()->SetTitle(Form("Number of Protons / %3.1f [MeV] ",bin_KE.get_width()));
    
    cout<<"Done!"<<endl;

}

void PIDTool::SetAnalysisMode(int nMode)
{
    if ( nMode == 1) {
        cout<<"\tAnalysis Mode: Signal - Only Signal Events will be Analysed"<<endl;
        branchDir = Folder_List::signal;
    }else if ( nMode == 2){
        cout<<"\tAnalysis Mode: Background - Only Background Events will be Analysed"<<endl;
        branchDir = Folder_List::background;
    }else{
        cout<<"\tAnalysis Mode: All - All Events will be Analysed"<<endl;
        branchDir = Folder_List::allEvents;
    }
}

void PIDTool::write_RootFile()
{
    cout<<">> Writing "<<rootDir<<endl;
    f->Write();
}


#endif

