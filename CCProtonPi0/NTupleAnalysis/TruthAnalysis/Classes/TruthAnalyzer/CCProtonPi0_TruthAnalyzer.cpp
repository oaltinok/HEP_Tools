#ifndef CCProtonPi0_TruthAnalyzer_cpp
#define CCProtonPi0_TruthAnalyzer_cpp


#include "../../../Classes/Analyzer/new_flux.h"

#include "CCProtonPi0_TruthAnalyzer.h"

using namespace PlotUtils;

void CCProtonPi0_TruthAnalyzer::Loop(std::string playlist)
{
    // Control Flow
    bool applyMaxEvents = false;
    double nMaxEvents = 100000;

    //------------------------------------------------------------------------
    // Create chain
    //------------------------------------------------------------------------
    TChain* fChain = new TChain("Truth");
    Init(playlist, fChain);

    if (!fChain) return;
    if (fChain == 0) return;
    
    // Read Flux & Calc Weights 
    new_flux::get().read_oldflux_histogram(Folder_List::rootDir_Flux_old);
    new_flux::get().read_newflux_histogram(Folder_List::rootDir_Flux_new);
    new_flux::get().calc_weights();

    // Disable Branches for Performance
    fChain->SetBranchStatus("*", false);
    fChain->SetBranchStatus("truth_is*", true);
    fChain->SetBranchStatus("truth_pi0_*", true);
    fChain->SetBranchStatus("truth_muon_*", true);
    fChain->SetBranchStatus("*genie_wgt_*", true);
    fChain->SetBranchStatus("mc_Q2", true);
    fChain->SetBranchStatus("mc_incomingE", true);

    //------------------------------------------------------------------------
    // Loop over Chain
    //------------------------------------------------------------------------
    std::cout<<"Looping over all entries"<<std::endl;

    Long64_t nentries = fChain->GetEntriesFast();
    Long64_t nbytes = 0, nb = 0;

    for (Long64_t jentry=0; jentry<nentries;jentry++) {

        nb = fChain->GetEntry(jentry);   nbytes += nb;    
        Long64_t ientry = fChain->GetEntry(jentry);

        if (ientry == 0) {
            std::cout<<"\tGetEntry failure "<<jentry<<std::endl;
            break;
        }

        // Progress Message on Terminal
        if (jentry % 1000000 == 0 ) std::cout<<"Entry = "<<jentry<<std::endl;

        if (applyMaxEvents && jentry == nMaxEvents){
            std::cout<<"\tReached Event Limit!"<<std::endl;
            break;
        }

        nAll++;
        if (truth_isFidVol) nFidVol++;
        else{
            nNoFidVol++;
            continue;
        }

        // Count Signal and Background
        if (truth_isSignal){
            FillSignalHistograms();
            nSignal++;
        }

        // Background With Pi0
        if (truth_isBckg_NoPi0) nBckg_NoPi0++;
        else if (truth_isBckg_SinglePi0) nBckg_SinglePi0++;
        else if (truth_isBckg_MultiPi0) nBckg_MultiPi0++;

        // Background Types
        if (truth_isBckg_NC) nBckg_NC++;
        else if (truth_isBckg_AntiNeutrino) nBckg_AntiNeutrino++;
        else if (truth_isBckg_QELike) nBckg_QELike++;
        else if (truth_isBckg_SingleChargedPion) nBckg_SingleChargedPion++;
        else if (truth_isBckg_SingleChargedPion_ChargeExchanged) nBckg_SingleChargedPion_ChargeExchanged++;
        else if (truth_isBckg_DoublePionWithPi0) nBckg_DoublePion_WithPi0++;
        else if (truth_isBckg_DoublePionWithoutPi0) nBckg_DoublePion_WithoutPi0++;
        else if (truth_isBckg_MultiPionWithPi0) nBckg_MultiPion_WithPi0++;
        else if (truth_isBckg_MultiPionWithoutPi0) nBckg_MultiPion_WithoutPi0++;
        else if (truth_isBckg_Other) nBckg_Other++;
        else if (truth_isFidVol && !truth_isSignal) std::cout<<"WARNING! No Background Type"<<std::endl;
    }

    // Add Other Error Bands and Fill With CV
    AddOtherErrorBands_FillWithCV();

    writeTextFile();
    writeHistograms();
}

void CCProtonPi0_TruthAnalyzer::writeTextFile() 
{
    double totalBackgroundWithPi0 = nBckg_NoPi0 + nBckg_SinglePi0 + nBckg_MultiPi0;
    double totalBackground = nBckg_NC + nBckg_AntiNeutrino + nBckg_QELike + nBckg_SingleChargedPion + nBckg_SingleChargedPion_ChargeExchanged + nBckg_DoublePion_WithPi0 + nBckg_DoublePion_WithoutPi0 + nBckg_MultiPion_WithPi0 + nBckg_MultiPion_WithoutPi0 + nBckg_Other;
    double totalEvents1 = nSignal + totalBackgroundWithPi0;
    double totalEvents2 = nSignal + totalBackground;

    // Formatting for Text Output
    textFile<<std::fixed;
    textFile<<std::setprecision(2);

    // Write to Text File
    textFile<<"nAll = "<<nAll<<std::endl;
    textFile<<"nFidVol = "<<nFidVol<<" "<<GetPercent(nAll,nFidVol)<<"%"<<std::endl;
    textFile<<"nNoFidVol = "<<nNoFidVol<<" "<<GetPercent(nAll,nNoFidVol)<<"%"<<std::endl;
    textFile<<"================================================================"<<std::endl;
    textFile<<"nSignal = "<<nSignal<<" "<<GetPercent(nFidVol,nSignal)<<"%"<<std::endl;
    textFile<<"----------------------------------------------------------------"<<std::endl;
    textFile<<"nBckg_NoPi0 = "<<nBckg_NoPi0<<" "<<GetPercent(nFidVol,nBckg_NoPi0)<<"%"<<std::endl;
    textFile<<"nBckg_SinglePi0 = "<<nBckg_SinglePi0<<" "<<GetPercent(nFidVol,nBckg_SinglePi0)<<"%"<<std::endl;
    textFile<<"nBckg_MultiPi0 = "<<nBckg_MultiPi0<<" "<<GetPercent(nFidVol,nBckg_MultiPi0)<<"%"<<std::endl;
    textFile<<"----------------------------------------------------------------"<<std::endl;
    textFile<<"nBckg_NC = "<<nBckg_NC<<" "<<GetPercent(nFidVol,nBckg_NC)<<"%"<<std::endl;
    textFile<<"nBckg_AntiNeutrino = "<<nBckg_AntiNeutrino<<" "<<GetPercent(nFidVol,nBckg_AntiNeutrino)<<"%"<<std::endl;
    textFile<<"nBckg_QELike = "<<nBckg_QELike<<" "<<GetPercent(nFidVol,nBckg_QELike)<<"%"<<std::endl;
    textFile<<"nBckg_SingleChargedPion = "<<nBckg_SingleChargedPion<<" "<<GetPercent(nFidVol,nBckg_SingleChargedPion)<<"%"<<std::endl;
    textFile<<"nBckg_SingleChargedPion_ChargeExchanged = "<<nBckg_SingleChargedPion_ChargeExchanged<<" "<<GetPercent(nFidVol,nBckg_SingleChargedPion_ChargeExchanged)<<"%"<<std::endl;
    textFile<<"nBckg_DoublePion_WithPi0 = "<<nBckg_DoublePion_WithPi0<<" "<<GetPercent(nFidVol,nBckg_DoublePion_WithPi0)<<"%"<<std::endl;
    textFile<<"nBckg_DoublePion_WithoutPi0 = "<<nBckg_DoublePion_WithoutPi0<<" "<<GetPercent(nFidVol,nBckg_DoublePion_WithoutPi0)<<"%"<<std::endl;
    textFile<<"nBckg_MultiPion_WithPi0 = "<<nBckg_MultiPion_WithPi0<<" "<<GetPercent(nFidVol,nBckg_MultiPion_WithPi0)<<"%"<<std::endl;
    textFile<<"nBckg_MultiPion_WithoutPi0 = "<<nBckg_MultiPion_WithoutPi0<<" "<<GetPercent(nFidVol,nBckg_MultiPion_WithoutPi0)<<"%"<<std::endl;
    textFile<<"nBckg_Other = "<<nBckg_Other<<" "<<GetPercent(nFidVol,nBckg_Other)<<"%"<<std::endl;
    textFile<<"----------------------------------------------------------------"<<std::endl;

    textFile<<"Total Signal & Background With Pi0 = "<<totalEvents1<<" "<<GetPercent(nFidVol,totalEvents1)<<"%"<<std::endl;
    textFile<<"Total Signal & Background = "<<totalEvents2<<" "<<GetPercent(nFidVol,totalEvents2)<<"%"<<std::endl;

    textFile.close();
}

double CCProtonPi0_TruthAnalyzer::GetPercent(double nAll, double nOther)
{
    double percent = (nOther/nAll) * 100;
    return percent;
}

CCProtonPi0_TruthAnalyzer::CCProtonPi0_TruthAnalyzer() : CCProtonPi0_NTupleAnalysis()

{
    std::cout<<"Initializing TruthAnalyzer!"<<std::endl;

    // Required for MINERvA Framework Classes
    ROOT::Cintex::Cintex::Enable();

    // Open ROOT File
    rootDir = Folder_List::rootDir_Truth_mc;

    std::cout<<"\tRoot File: "<<rootDir<<std::endl;
    f = new TFile(rootDir.c_str(),"RECREATE");

    initHistograms();
    
    openTextFiles();

    resetCounters();

    std::cout<<"Finished"<<std::endl;
}

void CCProtonPi0_TruthAnalyzer::resetCounters() 
{
    nAll = 0.0;
    nFidVol = 0.0;
    nNoFidVol = 0.0;
    nSignal = 0.0;

    // Background With Pi0
    nBckg_NoPi0 = 0.0;
    nBckg_SinglePi0 = 0.0;
    nBckg_MultiPi0 = 0.0;

    // Background Types
    nBckg_NC = 0.0;
    nBckg_AntiNeutrino = 0.0;
    nBckg_QELike = 0.0;
    nBckg_SingleChargedPion = 0.0;
    nBckg_SingleChargedPion_ChargeExchanged = 0.0;
    nBckg_DoublePion_WithPi0 = 0.0;
    nBckg_DoublePion_WithoutPi0 = 0.0;
    nBckg_MultiPion_WithPi0 = 0.0;
    nBckg_MultiPion_WithoutPi0 = 0.0;
    nBckg_Other = 0.0;
}


void CCProtonPi0_TruthAnalyzer::openTextFiles() 
{
    // Open TextFiles
    file_name = "TruthInfo.dat";
    textFile.open(file_name.c_str());
    if (!textFile.is_open()){
        std::cerr<<"Cannot open output text file: "<<file_name<<std::endl;
        exit(1);
    }else{
        std::cout<<"\t"<<file_name<<std::endl;
    }
}

void CCProtonPi0_TruthAnalyzer::initHistograms()
{
    // ------------------------------------------------------------------------
    // Muon Variables
    // ------------------------------------------------------------------------
    int nBins_muon_P = 10;
    double min_muon_P = 0.0;
    double max_muon_P = 10.0;
    muon_P_mc_truth_all_signal = new MnvH1D( "muon_P_mc_truth_all_signal","Muon Momentum for Signal Events",nBins_muon_P, min_muon_P, max_muon_P);
    muon_P_mc_truth_all_signal->GetXaxis()->SetTitle("Momentum [GeV]");
    muon_P_mc_truth_all_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBand_Flux(muon_P_mc_truth_all_signal);
    AddVertErrorBand_Genie(muon_P_mc_truth_all_signal);

    int nBins_muon_theta = 12;
    double min_muon_theta = 0.0;
    double max_muon_theta = 25.0;
    muon_theta_mc_truth_all_signal = new MnvH1D( "muon_theta_mc_truth_all_signal","Pi0 Muon Theta for Signal Events",nBins_muon_theta,min_muon_theta,max_muon_theta);
    muon_theta_mc_truth_all_signal->GetXaxis()->SetTitle("Theta");
    muon_theta_mc_truth_all_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBand_Flux(muon_theta_mc_truth_all_signal);
    AddVertErrorBand_Genie(muon_theta_mc_truth_all_signal);

    // ------------------------------------------------------------------------
    // Pi0 Variables
    // ------------------------------------------------------------------------
    int nBins_pi0_P = 17;
    double min_pi0_P = 0.0;
    double max_pi0_P = 1.7;
    pi0_P_mc_truth_all_signal = new MnvH1D( "pi0_P_mc_truth_all_signal","Pi0 Momentum for Signal Events",nBins_pi0_P, min_pi0_P, max_pi0_P);
    pi0_P_mc_truth_all_signal->GetXaxis()->SetTitle("Momentum [GeV]");
    pi0_P_mc_truth_all_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBand_Flux(pi0_P_mc_truth_all_signal);
    AddVertErrorBand_Genie(pi0_P_mc_truth_all_signal);

    int nBins_pi0_theta = 18;
    double min_pi0_theta = 0.0;
    double max_pi0_theta = 180.0;
    pi0_theta_mc_truth_all_signal = new MnvH1D( "pi0_theta_mc_truth_all_signal","Theta for Signal Events",nBins_pi0_theta, min_pi0_theta, max_pi0_theta);
    pi0_theta_mc_truth_all_signal->GetXaxis()->SetTitle("Theta");
    pi0_theta_mc_truth_all_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBand_Flux(pi0_theta_mc_truth_all_signal);
    AddVertErrorBand_Genie(pi0_theta_mc_truth_all_signal);

    // ------------------------------------------------------------------------
    // Neutrino Energy & Q2
    // ------------------------------------------------------------------------
    int nBins_neutrino_E = 20;
    double min_neutrino_E = 0.0;
    double max_neutrino_E = 20.0;
    neutrino_E_mc_truth_all_signal = new MnvH1D( "neutrino_E_mc_truth_all_signal","Neutrino Energy for Signal Events",nBins_neutrino_E, min_neutrino_E, max_neutrino_E);
    neutrino_E_mc_truth_all_signal->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
    neutrino_E_mc_truth_all_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBand_Flux(neutrino_E_mc_truth_all_signal);
    AddVertErrorBand_Genie(neutrino_E_mc_truth_all_signal);

    int nBins_QSq = 40;
    double min_QSq = 0.0;
    double max_QSq = 4.0;
    QSq_mc_truth_all_signal = new MnvH1D( "QSq_mc_truth_all_signal","Q^{2} for Signal Events",nBins_QSq,min_QSq,max_QSq);
    QSq_mc_truth_all_signal->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    QSq_mc_truth_all_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBand_Flux(QSq_mc_truth_all_signal);
    AddVertErrorBand_Genie(QSq_mc_truth_all_signal);
}

void CCProtonPi0_TruthAnalyzer::FillHistogram(MnvH1D* hist, double var)
{
    hist->Fill(var, cvweight);
    FillVertErrorBand_Flux(hist, var);
    FillVertErrorBand_Genie(hist, var);
}

void CCProtonPi0_TruthAnalyzer::FillVertErrorBand_Flux(MnvH1D* h, double var)
{
    double enu0 = mc_incomingE/1.e3;
    std::vector<double> random_weights = new_flux::get().get_random_weights(enu0);
    h->FillVertErrorBand("Flux",  var, &random_weights[0],  cvweight, 1.0);
}

void CCProtonPi0_TruthAnalyzer::FillVertErrorBand_Genie(MnvH1D* h, double var)
{
    h->FillVertErrorBand("GENIE_AGKYxF1pi"         ,var, truth_genie_wgt_AGKYxF1pi[2]        , truth_genie_wgt_AGKYxF1pi[4]        , cvweight);
    h->FillVertErrorBand("GENIE_AhtBY"             ,var, truth_genie_wgt_AhtBY[2]            , truth_genie_wgt_AhtBY[4]            , cvweight);
    h->FillVertErrorBand("GENIE_BhtBY"             ,var, truth_genie_wgt_BhtBY[2]            , truth_genie_wgt_BhtBY[4]            , cvweight);
    h->FillVertErrorBand("GENIE_CCQEPauliSupViaKF" ,var, truth_genie_wgt_CCQEPauliSupViaKF[2], truth_genie_wgt_CCQEPauliSupViaKF[4], cvweight);
    h->FillVertErrorBand("GENIE_CV1uBY"            ,var, truth_genie_wgt_CV1uBY[2]           , truth_genie_wgt_CV1uBY[4]           , cvweight);
    h->FillVertErrorBand("GENIE_CV2uBY"            ,var, truth_genie_wgt_CV2uBY[2]           , truth_genie_wgt_CV2uBY[4]           , cvweight);
    h->FillVertErrorBand("GENIE_EtaNCEL"           ,var, truth_genie_wgt_EtaNCEL[2]          , truth_genie_wgt_EtaNCEL[4]          , cvweight);
    h->FillVertErrorBand("GENIE_FrAbs_N"           ,var, truth_genie_wgt_FrAbs_N[2]          , truth_genie_wgt_FrAbs_N[4]          , cvweight);
    h->FillVertErrorBand("GENIE_FrAbs_pi"          ,var, truth_genie_wgt_FrAbs_pi[2]         , truth_genie_wgt_FrAbs_pi[4]         , cvweight);
    h->FillVertErrorBand("GENIE_FrCEx_N"           ,var, truth_genie_wgt_FrCEx_N[2]          , truth_genie_wgt_FrCEx_N[4]          , cvweight);
    h->FillVertErrorBand("GENIE_FrCEx_pi"          ,var, truth_genie_wgt_FrCEx_pi[2]         , truth_genie_wgt_FrCEx_pi[4]         , cvweight);
    h->FillVertErrorBand("GENIE_FrElas_N"          ,var, truth_genie_wgt_FrElas_N[2]         , truth_genie_wgt_FrElas_N[4]         , cvweight);
    h->FillVertErrorBand("GENIE_FrElas_pi"         ,var, truth_genie_wgt_FrElas_pi[2]        , truth_genie_wgt_FrElas_pi[4]        , cvweight);
    h->FillVertErrorBand("GENIE_FrInel_N"          ,var, truth_genie_wgt_FrInel_N[2]         , truth_genie_wgt_FrInel_N[4]         , cvweight);
    h->FillVertErrorBand("GENIE_FrInel_pi"         ,var, truth_genie_wgt_FrInel_pi[2]        , truth_genie_wgt_FrInel_pi[4]        , cvweight);
    h->FillVertErrorBand("GENIE_FrPiProd_N"        ,var, truth_genie_wgt_FrPiProd_N[2]       , truth_genie_wgt_FrPiProd_N[4]       , cvweight);
    h->FillVertErrorBand("GENIE_FrPiProd_pi"       ,var, truth_genie_wgt_FrPiProd_pi[2]      , truth_genie_wgt_FrPiProd_pi[4]      , cvweight);
    h->FillVertErrorBand("GENIE_MFP_N"             ,var, truth_genie_wgt_MFP_N[2]            , truth_genie_wgt_MFP_N[4]            , cvweight);
    h->FillVertErrorBand("GENIE_MFP_pi"            ,var, truth_genie_wgt_MFP_pi[2]           , truth_genie_wgt_MFP_pi[4]           , cvweight);
    h->FillVertErrorBand("GENIE_MaCCQE"            ,var, truth_genie_wgt_MaCCQE[2]           , truth_genie_wgt_MaCCQE[4]           , cvweight);
    h->FillVertErrorBand("GENIE_MaCCQEshape"       ,var, truth_genie_wgt_MaCCQEshape[2]      , truth_genie_wgt_MaCCQEshape[4]      , cvweight);
    h->FillVertErrorBand("GENIE_MaNCEL"            ,var, truth_genie_wgt_MaNCEL[2]           , truth_genie_wgt_MaNCEL[4]           , cvweight);
    h->FillVertErrorBand("GENIE_MaRES"             ,var, truth_genie_wgt_MaRES[2]            , truth_genie_wgt_MaRES[4]            , cvweight);
    h->FillVertErrorBand("GENIE_MvRES"             ,var, truth_genie_wgt_MvRES[2]            , truth_genie_wgt_MvRES[4]            , cvweight);
    h->FillVertErrorBand("GENIE_NormCCQE"          ,var, truth_genie_wgt_NormCCQE[2]         , truth_genie_wgt_NormCCQE[4]         , cvweight);
    h->FillVertErrorBand("GENIE_NormCCRES"         ,var, truth_genie_wgt_NormCCRES[2]        , truth_genie_wgt_NormCCRES[4]        , cvweight);
    h->FillVertErrorBand("GENIE_NormDISCC"         ,var, truth_genie_wgt_NormDISCC[2]        , truth_genie_wgt_NormDISCC[4]        , cvweight);
    h->FillVertErrorBand("GENIE_NormNCRES"         ,var, truth_genie_wgt_NormNCRES[2]        , truth_genie_wgt_NormNCRES[4]        , cvweight);
    h->FillVertErrorBand("GENIE_RDecBR1gamma"      ,var, truth_genie_wgt_RDecBR1gamma[2]     , truth_genie_wgt_RDecBR1gamma[4]     , cvweight);
    h->FillVertErrorBand("GENIE_Rvn1pi"            ,var, truth_genie_wgt_Rvn1pi[2]           , truth_genie_wgt_Rvn1pi[4]           , cvweight);
    h->FillVertErrorBand("GENIE_Rvn2pi"            ,var, truth_genie_wgt_Rvn2pi[2]           , truth_genie_wgt_Rvn2pi[4]           , cvweight);
    h->FillVertErrorBand("GENIE_Rvp1pi"            ,var, truth_genie_wgt_Rvp1pi[2]           , truth_genie_wgt_Rvp1pi[4]           , cvweight);
    h->FillVertErrorBand("GENIE_Rvp2pi"            ,var, truth_genie_wgt_Rvp2pi[2]           , truth_genie_wgt_Rvp2pi[4]           , cvweight);
    h->FillVertErrorBand("GENIE_Theta_Delta2Npi"   ,var, truth_genie_wgt_Theta_Delta2Npi[2]  , truth_genie_wgt_Theta_Delta2Npi[4]  , cvweight);
    h->FillVertErrorBand("GENIE_VecFFCCQEshape"    ,var, truth_genie_wgt_VecFFCCQEshape[2]   , truth_genie_wgt_VecFFCCQEshape[4]   , cvweight);
}

void CCProtonPi0_TruthAnalyzer::Calc_WeightFromSystematics()
{
    // Update cvweight with Flux
    double enu0 = mc_incomingE * MeV_to_GeV;
    cvweight = new_flux::get().get_cvweight(enu0);
    //std::cout<<"cvweight = "<<cvweight<<std::endl;
}

void CCProtonPi0_TruthAnalyzer::AddOtherErrorBands_FillWithCV()
{
    AddErrorBands_FillWithCV(muon_P_mc_truth_all_signal);
    AddErrorBands_FillWithCV(muon_theta_mc_truth_all_signal);
    AddErrorBands_FillWithCV(pi0_P_mc_truth_all_signal);
    AddErrorBands_FillWithCV(pi0_theta_mc_truth_all_signal);
    AddErrorBands_FillWithCV(neutrino_E_mc_truth_all_signal);
    AddErrorBands_FillWithCV(QSq_mc_truth_all_signal);
}

void CCProtonPi0_TruthAnalyzer::AddErrorBands_FillWithCV(MnvH1D* hist)
{
    AddVertErrorBandAndFillWithCV_MuonTracking(hist);
}

void CCProtonPi0_TruthAnalyzer::FillSignalHistograms()
{
    Calc_WeightFromSystematics();
    FillHistogram(muon_P_mc_truth_all_signal, truth_muon_P * MeV_to_GeV);
    FillHistogram(muon_theta_mc_truth_all_signal, truth_muon_theta * rad_to_deg);
    FillHistogram(pi0_P_mc_truth_all_signal, truth_pi0_P * MeV_to_GeV);
    FillHistogram(pi0_theta_mc_truth_all_signal, truth_pi0_theta * rad_to_deg);
    FillHistogram(neutrino_E_mc_truth_all_signal, mc_incomingE * MeV_to_GeV);
    FillHistogram(QSq_mc_truth_all_signal, mc_Q2 * MeVSq_to_GeVSq);
}

void CCProtonPi0_TruthAnalyzer::writeHistograms()
{
    f->cd();

    muon_P_mc_truth_all_signal->Write();
    muon_theta_mc_truth_all_signal->Write();
    pi0_P_mc_truth_all_signal->Write();
    pi0_theta_mc_truth_all_signal->Write();
    neutrino_E_mc_truth_all_signal->Write();
    QSq_mc_truth_all_signal->Write();

    f->Close();
}

#endif


