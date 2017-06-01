#ifndef CCProtonPi0_TruthAnalyzer_cpp
#define CCProtonPi0_TruthAnalyzer_cpp

#include "CCProtonPi0_TruthAnalyzer.h"

using namespace PlotUtils;

void CCProtonPi0_TruthAnalyzer::Loop(std::string playlist)
{
    // Control Flow
    bool applyMaxEvents = false;
    double nMaxEvents = 100000000;

    fillErrors_ByHand = true;
    
    applyGENIETuning_Complete = true;
    applyGENIETuning_DeltaSuppression = false;

    //------------------------------------------------------------------------
    // Create chain
    //------------------------------------------------------------------------
    TChain* fChain = new TChain("Truth");
    Init(playlist, fChain);

    if (!fChain) return;
    if (fChain == 0) return;
    
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

        CalcEventWeight();

        nAll.increment(cvweight);

        if (truth_isFidVol) nFidVol.increment(cvweight);
        else nFidVol_Out.increment(cvweight);

        // Count Events Only inside the Fidicual Volume
        if (truth_isFidVol){

            // QSq Study
            //if ( !IsEnuInRange(mc_incomingE)) continue;

            FillGENIE_Tuning();

            // Count Signal Out due to Acceptance or Kinematics
            if (truth_isSignalOut_Acceptance) nSignalOut_Acceptance.increment(cvweight);
            if (truth_isSignalOut_Kinematics) nSignalOut_Kinematics.increment(cvweight);

            // Count Signal and Background inside Fiducial Volume
            if (truth_isSignal) nSignal.increment(cvweight);
            else nBckg.increment(cvweight);

            // Count if Muon Angle is Small or Large

            TVector3 muon_3P(mc_primFSLepton[0], mc_primFSLepton[1], mc_primFSLepton[2]);
            TVector3 unit_muon_3P = muon_3P.Unit();
            TVector3 unit_z_axis(0.0, 0.0, 1.0);

            double muon_angle = acos(unit_muon_3P.Dot(unit_z_axis));
            muon_angle = muon_angle * rad_to_deg;
            if (muon_angle > 25) nMuonAngle_Large.increment(cvweight);
            else nMuonAngle_Small.increment(cvweight);
            
            if (truth_isSignalOut_Kinematics){
                FillSignalOut_Kinematics();    
            }

            if (truth_isSignal_BeforeFSI){
                nSignal_BeforeFSI.increment(cvweight);

                if (truth_W_exp_BeforeFSI > 0 && truth_W_exp_BeforeFSI < 1400){
                    GetDelta_pi_angles_BeforeFSI(); // Also Calculates the deltaInvMass_BeforeFSI
                }

                FillSignal_XSec_Variables_BeforeFSI();
            }

            CountSignalFeed();
         
            if (truth_isSignal){

                CalcDeltaRichKinematics();

                // Detailed Interaction breakout
                FillSignal_InteractionType(); // Also counts the N(Signal) in each type

                // Default Signal Histograms
                FillSignal_XSec_Variables();

                // FSI Type
                int FSIType = GetFSIType();
                CountFSIType(FSIType);
                FillSignal_XSec_Variables_FSIType(FSIType);

                // Int Type
                int IntType = GetIntType();
                FillSignal_XSec_Variables_IntType(IntType);
            }
        }
    }

    // Add Other Error Bands and Fill With CV
    AddOtherErrorBands_FillWithCV();

    writeTextFile();
    writeHistograms();
}

bool CCProtonPi0_TruthAnalyzer::isDeltaRichSample()
{
    return (truth_W_exp > 0.0 && truth_W_exp < 1400);
}

void CCProtonPi0_TruthAnalyzer::CalcDeltaRichKinematics()
{
    m_deltaInvMass = calcDeltaInvariantMass() * MeV_to_GeV;
    GetDelta_pi_angles();
}

void CCProtonPi0_TruthAnalyzer::CountSignalFeed()
{
    if (truth_isSignal_BeforeFSI && !truth_isSignal){
        nFSI_FeedOut.increment(cvweight);
        CountFSIType_FeedOut(); 
    }

    if ( !truth_isSignal_BeforeFSI && truth_isSignal){
        nFSI_FeedIn.increment(cvweight);
    }
}

void CCProtonPi0_TruthAnalyzer::CountFSIType_FeedOut()
{
    int nPi0 = 0;
    int nChargedPi = 0;

    // Count Number of Pions
    for (int i = 0; i < mc_er_nPart; ++i){
        if (mc_er_status[i] != 1) continue;
        if (mc_er_ID[i] == 111) nPi0++;    
        else if (std::abs(mc_er_ID[i]) == 211 ) nChargedPi++;
    }

    if ( (nPi0 + nChargedPi) == 0){
        nFSI_FeedOut_Abs.increment(cvweight);
    }else if ( (nPi0 + nChargedPi) > 1){
        nFSI_FeedOut_MultiPi.increment(cvweight);
    }else if (nPi0 == 0 && nChargedPi == 1){
        nFSI_FeedOut_Cex.increment(cvweight);
    }else if (nPi0 == 1 && nChargedPi == 0){
        nFSI_FeedOut_Other.increment(cvweight);
    }else{
        std::cout<<"WARNING! CountFSIType_FeedOut() Different FSI Type"<<std::endl;
    }
}


void CCProtonPi0_TruthAnalyzer::FillSignalOut_Kinematics()
{
    // Fill Enu even if Enu < 1.5 GeV
    if (isSignal_EnuLow()){
        FillHistogram(Enu_mc_truth_all_signal_AfterFSI, mc_incomingE * MeV_to_GeV);
        FillHistogram(Enu_mc_truth_all_signal_BeforeFSI, mc_incomingE * MeV_to_GeV);
    }

    // Fill W even if W > 1.8 GeV
    if (isSignal_WHigh()){
        if (mc_intType == 1){
            FillHistogram(truth_w_QE, truth_W_exp * MeV_to_GeV);
        }else if (mc_intType == 2){
            if (mc_resID == 0){ 
                FillHistogram(truth_w_RES_1232, truth_W_exp * MeV_to_GeV);
            }else if (mc_resID == 1){
                FillHistogram(truth_w_RES_1535, truth_W_exp * MeV_to_GeV);
            }else if (mc_resID == 2){
                FillHistogram(truth_w_RES_1520, truth_W_exp * MeV_to_GeV);
            }else{
                FillHistogram(truth_w_RES_Other, truth_W_exp * MeV_to_GeV);
            }
        }else if ( IsGenieNonRES() ){
            FillHistogram(truth_w_Non_RES, truth_W_exp * MeV_to_GeV);
        }else if (mc_intType == 3){
            FillHistogram(truth_w_DIS, truth_W_exp * MeV_to_GeV);
        }else if ( IsEvent2p2h(mc_intType) ){
            FillHistogram(truth_w_2p2h, truth_W_exp * MeV_to_GeV);
        }else{
            std::cout<<"WARNING! Signal Event with different interaction Type!"<<std::endl;
        }
    }
}

void CCProtonPi0_TruthAnalyzer::CountFSIType(int type)
{
    switch (type) {
        case 0: nFSI_NonInteracting.increment(cvweight);
                break;
        case 1: nFSI_Elastic.increment(cvweight);
                break;
        case 2: nFSI_Inelastic.increment(cvweight);
                break;
        case 3: nFSI_ChargeExchange.increment(cvweight);
                break;
        case 4: nFSI_MultiPi.increment(cvweight);
                break;
        case 5: nFSI_NucleonToPi.increment(cvweight);
                break;
        default: RunTimeError("Wrong FSI Type");
    }
}

int CCProtonPi0_TruthAnalyzer::GetIntType()
{
    int IntType = -1;
    
    if (mc_intType == 2){
        if (mc_resID == 0){ 
            // Delta Resonance
            IntType = 0;
        }else{
            // Other Resonance
            IntType = 1;
        }
    }else{
        // Non-Resonance 
        IntType = 2;
    }

    return IntType;
}

int CCProtonPi0_TruthAnalyzer::GetFSIType()
{
    const double eps = 0.99;
    //-------------------------------------------------------------------------
    // FSI Based on NPions
    //-------------------------------------------------------------------------
    int nPi0 = 0;
    int nChargedPi = 0;

    // Count Number of Pions
    for (int i = 0; i < mc_er_nPart; ++i){
        if (mc_er_status[i] != 14) continue;
        if (mc_er_ID[i] == 111) nPi0++;    
        else if (std::abs(mc_er_ID[i]) == 211 ) nChargedPi++;
    }

    if ((nChargedPi + nPi0) > 1){
        return 4;
    }else if (nPi0 == 0 && nChargedPi == 1){
        return 3;
    }else if (nPi0 == 0){
        return 5;
    }

    //-------------------------------------------------------------------------
    // FSI Based on Pion Scattering 
    //-------------------------------------------------------------------------
    if (nPi0 == 1){
        double E[3];
        double KE[3];
        double Px[3];
        double Py[3];
        double Pz[3];
        // Loop again to get Pi0 Initial and Final Energy
        for (int i = 0; i < mc_er_nPart; ++i){
            if (mc_er_ID[i] != 111) continue;
            
            if (mc_er_status[i] == 14){
                E[0] = mc_er_E[i];
                KE[0] = E[0]-pi0_mass;
                Px[0] = mc_er_Px[i];
                Py[0] = mc_er_Py[i];
                Pz[0] = mc_er_Pz[i];
            }else if (mc_er_status[i] == 1){
                E[1] = mc_er_E[i];
                KE[1] = E[1]-pi0_mass;
                Px[1] = mc_er_Px[i];
                Py[1] = mc_er_Py[i];
                Pz[1] = mc_er_Pz[i];
            }
        }

        // Calculate the Initial and Final Difference
        E[2] = std::abs(E[0]-E[1]);
        KE[2] = std::abs(KE[0]-KE[1]);
        Px[2] = std::abs(Px[0]-Px[1]);
        Py[2] = std::abs(Py[0]-Py[1]);
        Pz[2] = std::abs(Pz[0]-Pz[1]);

        bool isKE_Changed = KE[2] > eps;
        bool isDirection_Changed = (Px[2] > eps || Py[2] > eps || Pz[2] > eps);

        if (isKE_Changed){
            // Kinetic Energy is Changed
            return 2;
        }else if (isDirection_Changed){
            // KE Conserved, direction changed
            return 1;
        }else{
            // Nothing Changed
            return 0;
        }
    } 
    
    std::cout<<"WARNING! In this stage there must be single Pi0"<<std::endl;
    return -1;
    
}

void CCProtonPi0_TruthAnalyzer::writeTextFile() 
{
    // Formatting for Text Output
    textFile<<std::fixed;
    textFile<<std::setprecision(1);

    // All Events
    WriteCounter(nAll, nAll);
    textFile<<std::endl;
    
    WriteCounter(nFidVol, nAll);
    WriteCounter(nFidVol_Out, nAll);
    textFile<<std::endl;

    // Events inside Fiducial Volume
    WriteCounter(nSignal, nFidVol);
    WriteCounter(nSignal_BeforeFSI, nFidVol);
    WriteCounter(nSignalOut_Acceptance, nFidVol);
    WriteCounter(nSignalOut_Kinematics, nFidVol);
    WriteCounter(nBckg, nFidVol);
    WriteCounter(nMuonAngle_Small, nFidVol);
    WriteCounter(nMuonAngle_Large, nFidVol);
    textFile<<std::endl;

    // FSI Types
    WriteCounter(nFSI_FeedOut, nBckg);
    WriteCounter(nFSI_FeedOut_Abs, nBckg);
    WriteCounter(nFSI_FeedOut_MultiPi, nBckg);
    WriteCounter(nFSI_FeedOut_Cex, nBckg);
    WriteCounter(nFSI_FeedOut_Other, nBckg);

    textFile<<std::endl;
    WriteCounter(nFSI_FeedIn, nSignal);
    WriteCounter(nFSI_NonInteracting, nSignal);
    WriteCounter(nFSI_Elastic, nSignal);
    WriteCounter(nFSI_Inelastic, nSignal);
    WriteCounter(nFSI_ChargeExchange, nSignal);
    WriteCounter(nFSI_MultiPi, nSignal);
    WriteCounter(nFSI_NucleonToPi, nSignal);
    textFile<<std::endl;

    // Signal Types
    WriteCounter(nQE, nSignal);
    textFile<<std::endl;

    WriteCounter(nRES_1232, nSignal);
    WriteCounter(nRES_1535, nSignal);
    WriteCounter(nRES_1520, nSignal);
    WriteCounter(nRES_Other, nSignal);
    textFile<<std::endl;
   
    WriteCounter(nDIS, nSignal);
    WriteCounter(nNon_RES, nSignal);
    WriteCounter(n2p2h, nSignal);
    textFile<<std::endl;
  
    WriteCounter(nTempCounter1, nAll);
    WriteCounter(nTempCounter2, nAll);
    textFile<<std::endl;

    WriteCounter(DeltaStudy_nAll, nSignal);
    WriteCounter(DeltaStudy_nDeltaRES, DeltaStudy_nAll);
    WriteCounter(DeltaStudy_nOtherRES, DeltaStudy_nAll);
    WriteCounter(DeltaStudy_nNonRES, DeltaStudy_nAll);
    textFile<<std::endl;

    WriteCounter(nShowers, nShowers);
    WriteCounter(nShowers_Out, nShowers);
    textFile<<std::endl;

    textFile.close();
}

void CCProtonPi0_TruthAnalyzer::WriteCounter(CCProtonPi0_Counter Counter, CCProtonPi0_Counter PercentBase)
{
    textFile<<Counter.getName()<<"\t"<<Counter.getCount()<<"\t"<<GetPercent(PercentBase,Counter)<<std::endl;
}

double CCProtonPi0_TruthAnalyzer::GetPercent(CCProtonPi0_Counter nAll, CCProtonPi0_Counter nOther)
{
    double percent = (nOther.getCount()/nAll.getCount()) * 100;
    return percent;
}

CCProtonPi0_TruthAnalyzer::CCProtonPi0_TruthAnalyzer() : 
    CCProtonPi0_NTupleAnalysis(),
    BckgConstrainer(Folder_List::BckgConstraints_TruthAnalysis),
    QSqFitter()
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
    nAll.setName("nAll");

    nTempCounter1.setName("Muon");
    nTempCounter2.setName("AntiMuon");

    nFidVol.setName("nFidVol");
    nFidVol_Out.setName("nFidVol_Out");

    nSignal.setName("nSignal");
    nSignal_BeforeFSI.setName("nSignal_BeforeFSI");
    nSignalOut_Acceptance.setName("nSignalOut_Acceptance");
    nSignalOut_Kinematics.setName("nSignalOut_Kinematics");
    nBckg.setName("nBckg");
    nMuonAngle_Small.setName("nMuonAngle_Small");
    nMuonAngle_Large.setName("nMuonAngle_Large");

    // FSI 
    nFSI_FeedOut.setName("nFSI_FeedOut");
    nFSI_FeedOut_Abs.setName("nFSI_FeedOut_Abs");
    nFSI_FeedOut_MultiPi.setName("nFSI_FeedOut_MultiPi");
    nFSI_FeedOut_Cex.setName("nFSI_FeedOut_Cex");
    nFSI_FeedOut_Other.setName("nFSI_FeedOut_Other");


    nFSI_FeedIn.setName("nFSI_FeedIn");
    nFSI_NonInteracting.setName("nFSI_NonInteracting");
    nFSI_Elastic.setName("nFSI_Elastic");
    nFSI_Inelastic.setName("nFSI_Inelastic");
    nFSI_ChargeExchange.setName("nFSI_ChargeExchange");
    nFSI_MultiPi.setName("nFSI_MultiPi");
    nFSI_NucleonToPi.setName("nFSI_NucleonToPi");

    // Signal Type
    nQE.setName("nQuasi_Elastic");
    
    nRES_1232.setName("nSignal_RES_Delta");
    nRES_1535.setName("nSignal_RES_1535");
    nRES_1520.setName("nSignal_RES_1520");
    nRES_Other.setName("nSignal_RES_Other");
    
    nDIS.setName("nSignal_DIS");
    nNon_RES.setName("nSignal_NonRES");
    n2p2h.setName("nSignal_2p2h");
    
    // Delta Study
    DeltaStudy_nAll.setName("DeltaStudy_nAll");
    DeltaStudy_nDeltaRES.setName("DeltaStudy_nDeltaRES");
    DeltaStudy_nOtherRES.setName("DeltaStudy_nOtherRES");
    DeltaStudy_nNonRES.setName("DeltaStudy_nNonRES");

    nShowers.setName("nShowers");
    nShowers_Out.setName("nShowers_Out");
}

void CCProtonPi0_TruthAnalyzer::openTextFiles() 
{
    // Open TextFiles
    OpenTextFile("TruthInfo.txt", textFile);
    OpenTextFile("EventRecord.txt", logFile);

    logFile<<std::left;
    logFile.width(8); logFile<<"PDG"<<" "; 
    logFile.width(8); logFile<<"Status"<<" "; 
    logFile.width(8); logFile<<"Mother"<<" ";    
    logFile.width(12); logFile<<"ER E"<<" ";      
    logFile.width(12); logFile<<"Traj E"<<" ";      
    logFile<<endl;

}

void CCProtonPi0_TruthAnalyzer::initHistograms()
{
    MnvH1D* temp = NULL;
    for (int i = 0; i < nHistograms; i++){
        temp = new MnvH1D( Form("%s_%d","CV_weight_2p2h",i),"Central Value Weight -- 2p2h", 80,0.0,2.0);
        temp->GetXaxis()->SetTitle("Central Value Weight -- 2p2h");
        temp->GetYaxis()->SetTitle("N(Events)");
        CV_weight_2p2h.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","CV_weight_Delta",i),"Central Value Weight -- Delta", 80,0.0,2.0);
        temp->GetXaxis()->SetTitle("Central Value Weight -- Delta");
        temp->GetYaxis()->SetTitle("N(Events)");
        CV_weight_Delta.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","CV_weight_CCRES",i),"Central Value Weight -- CCRES", 80,0.0,2.0);
        temp->GetXaxis()->SetTitle("Central Value Weight -- CCRES");
        temp->GetYaxis()->SetTitle("N(Events)");
        CV_weight_CCRES.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","CV_weight_NonRes1pi",i),"Central Value Weight -- NonRes1pi", 80,0.0,2.0);
        temp->GetXaxis()->SetTitle("Central Value Weight -- NonRes1pi");
        temp->GetYaxis()->SetTitle("N(Events)");
        CV_weight_NonRes1pi.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","CV_weight",i),"Central Value Weight", 80,0.0,2.0);
        temp->GetXaxis()->SetTitle("Central Value Weight");
        temp->GetYaxis()->SetTitle("N(Events)");
        CV_weight.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","h_err_2p2h",i),"2p2h 1#sigma Weights", 80,0.0,2.0);
        temp->GetXaxis()->SetTitle("2p2h 1#sigma Weights");
        temp->GetYaxis()->SetTitle("N(Events)");
        h_err_2p2h.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","genie_wgt_Theta_Delta2Npi",i),"Theta_Delta2Npi GENIE 1#sigma Weights", 80,0.0,2.0);
        temp->GetXaxis()->SetTitle("Theta_Delta2Npi GENIE 1#sigma Weights");
        temp->GetYaxis()->SetTitle("N(Events)");
        genie_wgt_Theta_Delta2Npi.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","updated_wgt_Theta_Delta2Npi",i),"Theta_Delta2Npi Updated 1#sigma Weights", 80,0.0,2.0);
        temp->GetXaxis()->SetTitle("Theta_Delta2Npi Updated 1#sigma Weights");
        temp->GetYaxis()->SetTitle("N(Events)");
        updated_wgt_Theta_Delta2Npi.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","genie_wgt_MaRES",i),"MaRES GENIE 1#sigma Weights", 80,0.0,2.0);
        temp->GetXaxis()->SetTitle("MaRES GENIE 1#sigma Weights");
        temp->GetYaxis()->SetTitle("N(Events)");
        genie_wgt_MaRES.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","updated_wgt_MaRES",i),"MaRES Updated 1#sigma Weights", 80,0.0,2.0);
        temp->GetXaxis()->SetTitle("MaRES Updated 1#sigma Weights");
        temp->GetYaxis()->SetTitle("N(Events)");
        updated_wgt_MaRES.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","genie_wgt_MvRES",i),"MvRES GENIE 1#sigma Weights", 80,0.0,2.0);
        temp->GetXaxis()->SetTitle("MvRES GENIE 1#sigma Weights");
        temp->GetYaxis()->SetTitle("N(Events)");
        genie_wgt_MvRES.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","updated_wgt_MvRES",i),"MvRES Updated 1#sigma Weights", 80,0.0,2.0);
        temp->GetXaxis()->SetTitle("MvRES Updated 1#sigma Weights");
        temp->GetYaxis()->SetTitle("N(Events)");
        updated_wgt_MvRES.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","genie_wgt_Rvn1pi",i),"Rvn1pi GENIE 1#sigma Weights", 80,0.0,2.0);
        temp->GetXaxis()->SetTitle("Rvn1pi GENIE 1#sigma Weights");
        temp->GetYaxis()->SetTitle("N(Events)");
        genie_wgt_Rvn1pi.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","updated_wgt_Rvn1pi",i),"Rvn1pi Updated 1#sigma Weights", 80,0.0,2.0);
        temp->GetXaxis()->SetTitle("Rvn1pi Updated 1#sigma Weights");
        temp->GetYaxis()->SetTitle("N(Events)");
        updated_wgt_Rvn1pi.push_back(temp);
    }

    // ------------------------------------------------------------------------
    // Test Variables
    // ------------------------------------------------------------------------
    Test_pi0_P = new TH1D("Test_pi0_P","Pi0 Momentum Check",40,-100,100);
    Test_pi0_P->GetXaxis()->SetTitle("P_{Trajectory} - P_{EventRecord} [Mev]");
    Test_pi0_P->GetYaxis()->SetTitle("NEvents()");

    Test_pi0_theta = new TH1D("Test_pi0_theta","Pi0 Angle Check",40,-10,10);
    Test_pi0_theta->GetXaxis()->SetTitle("#theta_{Trajectory} - #theta_{EventRecord} [degree]");
    Test_pi0_theta->GetYaxis()->SetTitle("NEvents()");

    // ------------------------------------------------------------------------
    // Muon Variables
    // ------------------------------------------------------------------------
    muon_P_mc_truth_all_signal = new MnvH1D( "muon_P_mc_truth_all_signal","Muon Momentum for Signal Events",binList.size_muon_P, binList.a_muon_P);
    muon_P_mc_truth_all_signal->GetXaxis()->SetTitle("P_{#mu} [GeV]");
    muon_P_mc_truth_all_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(muon_P_mc_truth_all_signal);

    muon_theta_mc_truth_all_signal = new MnvH1D( "muon_theta_mc_truth_all_signal","Muon Theta for Signal Events",binList.size_muon_theta, binList.a_muon_theta);
    muon_theta_mc_truth_all_signal->GetXaxis()->SetTitle("#theta_{#mu} [deg]");
    muon_theta_mc_truth_all_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(muon_theta_mc_truth_all_signal);

    // ------------------------------------------------------------------------
    // Pi0 Variables
    // ------------------------------------------------------------------------
    pi0_P_mc_truth_all_signal = new MnvH1D( "pi0_P_mc_truth_all_signal","Pi0 Momentum for Signal Events",binList.size_pi0_P, binList.a_pi0_P);
    pi0_P_mc_truth_all_signal->GetXaxis()->SetTitle("P_{#pi^{0}} [GeV]");
    pi0_P_mc_truth_all_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(pi0_P_mc_truth_all_signal);

    pi0_KE_mc_truth_all_signal = new MnvH1D( "pi0_KE_mc_truth_all_signal","Pi0 Kinetic Energy for Signal Events",binList.size_pi0_KE, binList.a_pi0_KE);
    pi0_KE_mc_truth_all_signal->GetXaxis()->SetTitle("T_{#pi^{0}} [GeV]");
    pi0_KE_mc_truth_all_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(pi0_KE_mc_truth_all_signal);
    
    pi0_theta_mc_truth_all_signal = new MnvH1D( "pi0_theta_mc_truth_all_signal","Pi0 Theta for Signal Events",binList.size_pi0_theta, binList.a_pi0_theta);
    pi0_theta_mc_truth_all_signal->GetXaxis()->SetTitle("#theta_{#pi^{0}}[deg]");
    pi0_theta_mc_truth_all_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(pi0_theta_mc_truth_all_signal);

    // ------------------------------------------------------------------------
    // Neutrino Energy, Q2 & W
    // ------------------------------------------------------------------------
    Enu_mc_truth_all_signal = new MnvH1D( "Enu_mc_truth_all_signal","Neutrino Energy for Signal Events",binList.size_Enu, binList.a_Enu);
    Enu_mc_truth_all_signal->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
    Enu_mc_truth_all_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(Enu_mc_truth_all_signal);

    QSq_mc_truth_all_signal = new MnvH1D( "QSq_mc_truth_all_signal","Q^{2} for Signal Events",binList.size_QSq, binList.a_QSq);
    QSq_mc_truth_all_signal->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    QSq_mc_truth_all_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(QSq_mc_truth_all_signal);

    W_mc_truth_all_signal = new MnvH1D( "W_mc_truth_all_signal","W for Signal Events",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max());
    W_mc_truth_all_signal->GetXaxis()->SetTitle("W_{exp} (GeV)");
    W_mc_truth_all_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(W_mc_truth_all_signal);

    deltaInvMass_mc_truth_all_signal = new MnvH1D( "deltaInvMass_mc_truth_all_signal","deltaInvMass for Signal Events",binList.deltaInvMass.get_nBins(), binList.deltaInvMass.get_min(), binList.deltaInvMass.get_max());
    deltaInvMass_mc_truth_all_signal->GetXaxis()->SetTitle("deltaInvMass");
    deltaInvMass_mc_truth_all_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(deltaInvMass_mc_truth_all_signal);

    Delta_pi_theta_mc_truth_all_signal = new MnvH1D( "Delta_pi_theta_mc_truth_all_signal","Delta_pi_theta for Signal Events",binList.Delta_pi_theta.get_nBins(), binList.Delta_pi_theta.get_min(), binList.Delta_pi_theta.get_max());
    Delta_pi_theta_mc_truth_all_signal->GetXaxis()->SetTitle("Delta_pi_theta");
    Delta_pi_theta_mc_truth_all_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(Delta_pi_theta_mc_truth_all_signal);

    Delta_pi_phi_mc_truth_all_signal = new MnvH1D( "Delta_pi_phi_mc_truth_all_signal","Delta_pi_phi for Signal Events",binList.Delta_pi_phi.get_nBins(), binList.Delta_pi_phi.get_min(), binList.Delta_pi_phi.get_max());
    Delta_pi_phi_mc_truth_all_signal->GetXaxis()->SetTitle("Delta_pi_phi");
    Delta_pi_phi_mc_truth_all_signal->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(Delta_pi_phi_mc_truth_all_signal);

    // ------------------------------------------------------------------------
    // Muon Variables -- Before FSI
    // ------------------------------------------------------------------------
    muon_P_mc_truth_all_signal_BeforeFSI = new MnvH1D( "muon_P_mc_truth_all_signal_BeforeFSI","Muon Momentum for Signal Events",40, 0.0, 10.0);
    muon_P_mc_truth_all_signal_BeforeFSI->GetXaxis()->SetTitle("Momentum [GeV]");
    muon_P_mc_truth_all_signal_BeforeFSI->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(muon_P_mc_truth_all_signal_BeforeFSI);

    muon_theta_mc_truth_all_signal_BeforeFSI = new MnvH1D( "muon_theta_mc_truth_all_signal_BeforeFSI","Muon Theta for Signal Events",40, 0.0, 25.0);
    muon_theta_mc_truth_all_signal_BeforeFSI->GetXaxis()->SetTitle("Theta");
    muon_theta_mc_truth_all_signal_BeforeFSI->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(muon_theta_mc_truth_all_signal_BeforeFSI);

    // ------------------------------------------------------------------------
    // Pi0 Variables -- Before FSI
    // ------------------------------------------------------------------------
    pi0_P_mc_truth_all_signal_BeforeFSI = new MnvH1D( "pi0_P_mc_truth_all_signal_BeforeFSI","Pi0 Momentum for Signal Events",40, 0.0, 1.50);
    pi0_P_mc_truth_all_signal_BeforeFSI->GetXaxis()->SetTitle("Momentum [GeV]");
    pi0_P_mc_truth_all_signal_BeforeFSI->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(pi0_P_mc_truth_all_signal_BeforeFSI);

    pi0_KE_mc_truth_all_signal_BeforeFSI = new MnvH1D( "pi0_KE_mc_truth_all_signal_BeforeFSI","Pi0 Kinetic Energy for Signal Events",40, 0.0, 1.0);
    pi0_KE_mc_truth_all_signal_BeforeFSI->GetXaxis()->SetTitle("Kinetic Energy [GeV]");
    pi0_KE_mc_truth_all_signal_BeforeFSI->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(pi0_KE_mc_truth_all_signal_BeforeFSI);

    pi0_theta_mc_truth_all_signal_BeforeFSI = new MnvH1D( "pi0_theta_mc_truth_all_signal_BeforeFSI","Pi0 Theta for Signal Events",40, 0.0, 180.0);
    pi0_theta_mc_truth_all_signal_BeforeFSI->GetXaxis()->SetTitle("Theta");
    pi0_theta_mc_truth_all_signal_BeforeFSI->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(pi0_theta_mc_truth_all_signal_BeforeFSI);

    // ------------------------------------------------------------------------
    // Neutrino Energy, Q2 & W -- Before FSI
    // ------------------------------------------------------------------------
    Enu_mc_truth_all_signal_BeforeFSI = new MnvH1D( "Enu_mc_truth_all_signal_BeforeFSI","Neutrino Energy for Signal Events",binList.size_Enu_Fine, binList.a_Enu_Fine);
    Enu_mc_truth_all_signal_BeforeFSI->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
    Enu_mc_truth_all_signal_BeforeFSI->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(Enu_mc_truth_all_signal_BeforeFSI);

    QSq_mc_truth_all_signal_BeforeFSI = new MnvH1D( "QSq_mc_truth_all_signal_BeforeFSI","Q^{2} for Signal Events",40, 0.0, 2.0);
    QSq_mc_truth_all_signal_BeforeFSI->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    QSq_mc_truth_all_signal_BeforeFSI->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(QSq_mc_truth_all_signal_BeforeFSI);

    W_mc_truth_all_signal_BeforeFSI = new MnvH1D( "W_mc_truth_all_signal_BeforeFSI","W for Signal Events",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max());
    W_mc_truth_all_signal_BeforeFSI->GetXaxis()->SetTitle("W_{exp} (GeV)");
    W_mc_truth_all_signal_BeforeFSI->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(W_mc_truth_all_signal_BeforeFSI);

    deltaInvMass_mc_truth_all_signal_BeforeFSI = new MnvH1D( "deltaInvMass_mc_truth_all_signal_BeforeFSI","",40, 1.0, 1.6);
    deltaInvMass_mc_truth_all_signal_BeforeFSI->GetXaxis()->SetTitle("deltaInvMass");
    deltaInvMass_mc_truth_all_signal_BeforeFSI->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(deltaInvMass_mc_truth_all_signal_BeforeFSI);

    Delta_pi_theta_mc_truth_all_signal_BeforeFSI = new MnvH1D( "Delta_pi_theta_mc_truth_all_signal_BeforeFSI","",40, -1.0, 1.0);
    Delta_pi_theta_mc_truth_all_signal_BeforeFSI->GetXaxis()->SetTitle("Delta_pi_theta");
    Delta_pi_theta_mc_truth_all_signal_BeforeFSI->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(Delta_pi_theta_mc_truth_all_signal_BeforeFSI);

    Delta_pi_phi_mc_truth_all_signal_BeforeFSI = new MnvH1D( "Delta_pi_phi_mc_truth_all_signal_BeforeFSI","",40, 0.0, 360.0);
    Delta_pi_phi_mc_truth_all_signal_BeforeFSI->GetXaxis()->SetTitle("Delta_pi_phi");
    Delta_pi_phi_mc_truth_all_signal_BeforeFSI->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(Delta_pi_phi_mc_truth_all_signal_BeforeFSI);

    // ------------------------------------------------------------------------
    // Muon Variables -- After FSI
    // ------------------------------------------------------------------------
    muon_P_mc_truth_all_signal_AfterFSI = new MnvH1D( "muon_P_mc_truth_all_signal_AfterFSI","Muon Momentum for Signal Events",40, 0.0, 10.0);
    muon_P_mc_truth_all_signal_AfterFSI->GetXaxis()->SetTitle("Momentum [GeV]");
    muon_P_mc_truth_all_signal_AfterFSI->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(muon_P_mc_truth_all_signal_AfterFSI);

    muon_theta_mc_truth_all_signal_AfterFSI = new MnvH1D( "muon_theta_mc_truth_all_signal_AfterFSI","Muon Theta for Signal Events",40, 0.0, 25.0);
    muon_theta_mc_truth_all_signal_AfterFSI->GetXaxis()->SetTitle("Theta");
    muon_theta_mc_truth_all_signal_AfterFSI->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(muon_theta_mc_truth_all_signal_AfterFSI);

    // ------------------------------------------------------------------------
    // Pi0 Variables -- After FSI
    // ------------------------------------------------------------------------
    pi0_P_mc_truth_all_signal_AfterFSI = new MnvH1D( "pi0_P_mc_truth_all_signal_AfterFSI","Pi0 Momentum for Signal Events",40, 0.0, 1.50);
    pi0_P_mc_truth_all_signal_AfterFSI->GetXaxis()->SetTitle("Momentum [GeV]");
    pi0_P_mc_truth_all_signal_AfterFSI->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(pi0_P_mc_truth_all_signal_AfterFSI);

    pi0_KE_mc_truth_all_signal_AfterFSI = new MnvH1D( "pi0_KE_mc_truth_all_signal_AfterFSI","Pi0 Kinetic Energy for Signal Events",40, 0.0, 1.0);
    pi0_KE_mc_truth_all_signal_AfterFSI->GetXaxis()->SetTitle("Kinetic Energy [GeV]");
    pi0_KE_mc_truth_all_signal_AfterFSI->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(pi0_KE_mc_truth_all_signal_AfterFSI);

    pi0_theta_mc_truth_all_signal_AfterFSI = new MnvH1D( "pi0_theta_mc_truth_all_signal_AfterFSI","Pi0 Theta for Signal Events",40, 0.0, 180.0);
    pi0_theta_mc_truth_all_signal_AfterFSI->GetXaxis()->SetTitle("Theta");
    pi0_theta_mc_truth_all_signal_AfterFSI->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(pi0_theta_mc_truth_all_signal_AfterFSI);

    // ------------------------------------------------------------------------
    // Neutrino Energy, Q2 & W -- After FSI
    // ------------------------------------------------------------------------
    Enu_mc_truth_all_signal_AfterFSI = new MnvH1D( "Enu_mc_truth_all_signal_AfterFSI","Neutrino Energy for Signal Events",binList.size_Enu_Fine, binList.a_Enu_Fine);
    Enu_mc_truth_all_signal_AfterFSI->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
    Enu_mc_truth_all_signal_AfterFSI->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(Enu_mc_truth_all_signal_AfterFSI);

    QSq_mc_truth_all_signal_AfterFSI = new MnvH1D( "QSq_mc_truth_all_signal_AfterFSI","Q^{2} for Signal Events",40, 0.0, 2.0);
    QSq_mc_truth_all_signal_AfterFSI->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    QSq_mc_truth_all_signal_AfterFSI->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(QSq_mc_truth_all_signal_AfterFSI);

    W_mc_truth_all_signal_AfterFSI = new MnvH1D( "W_mc_truth_all_signal_AfterFSI","W for Signal Events",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max());
    W_mc_truth_all_signal_AfterFSI->GetXaxis()->SetTitle("W_{exp} (GeV)");
    W_mc_truth_all_signal_AfterFSI->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(W_mc_truth_all_signal_AfterFSI);

    deltaInvMass_mc_truth_all_signal_AfterFSI = new MnvH1D( "deltaInvMass_mc_truth_all_signal_AfterFSI","",40, 1.0, 1.6);
    deltaInvMass_mc_truth_all_signal_AfterFSI->GetXaxis()->SetTitle("deltaInvMass");
    deltaInvMass_mc_truth_all_signal_AfterFSI->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(deltaInvMass_mc_truth_all_signal_AfterFSI);

    Delta_pi_theta_mc_truth_all_signal_AfterFSI = new MnvH1D( "Delta_pi_theta_mc_truth_all_signal_AfterFSI","",40, -1.0, 1.0);
    Delta_pi_theta_mc_truth_all_signal_AfterFSI->GetXaxis()->SetTitle("Delta_pi_theta");
    Delta_pi_theta_mc_truth_all_signal_AfterFSI->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(Delta_pi_theta_mc_truth_all_signal_AfterFSI);

    Delta_pi_phi_mc_truth_all_signal_AfterFSI = new MnvH1D( "Delta_pi_phi_mc_truth_all_signal_AfterFSI","",40, 0.0, 360.0);
    Delta_pi_phi_mc_truth_all_signal_AfterFSI->GetXaxis()->SetTitle("Delta_pi_phi");
    Delta_pi_phi_mc_truth_all_signal_AfterFSI->GetYaxis()->SetTitle("N(Events)");
    AddVertErrorBands_TruthTree(Delta_pi_phi_mc_truth_all_signal_AfterFSI);

    // ------------------------------------------------------------------------
    // FSI Type 
    // ------------------------------------------------------------------------
    for (int i = 0; i < nFSIType; ++i){
        temp = new MnvH1D( Form("%s_%d","muon_P_mc_truth_all_signal_FSIType",i),"Muon Momentum for Signal Events",binList.size_muon_P, binList.a_muon_P);
        temp->GetXaxis()->SetTitle("Momentum [GeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBands_TruthTree(temp);
        muon_P_mc_truth_all_signal_FSIType.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","muon_theta_mc_truth_all_signal_FSIType",i),"Muon Theta for Signal Events",binList.size_muon_theta, binList.a_muon_theta);
        temp->GetXaxis()->SetTitle("Theta");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBands_TruthTree(temp);
        muon_theta_mc_truth_all_signal_FSIType.push_back(temp);
        
        temp = new MnvH1D( Form("%s_%d","pi0_P_mc_truth_all_signal_FSIType",i),"Pi0 Momentum for Signal Events",binList.size_pi0_P, binList.a_pi0_P);
        temp->GetXaxis()->SetTitle("Momentum [GeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBands_TruthTree(temp);
        pi0_P_mc_truth_all_signal_FSIType.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","pi0_KE_mc_truth_all_signal_FSIType",i),"Pi0 Kinetic Energy for Signal Events",binList.size_pi0_KE, binList.a_pi0_KE);
        temp->GetXaxis()->SetTitle("Kinetic Energy [GeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBands_TruthTree(temp);
        pi0_KE_mc_truth_all_signal_FSIType.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","pi0_theta_mc_truth_all_signal_FSIType",i), "Pi0 Theta for Signal Events",binList.size_pi0_theta, binList.a_pi0_theta);
        temp->GetXaxis()->SetTitle("Theta");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBands_TruthTree(temp);
        pi0_theta_mc_truth_all_signal_FSIType.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","Enu_mc_truth_all_signal_FSIType",i),"Neutrino Energy for Signal Events",binList.size_Enu, binList.a_Enu);
        temp->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBands_TruthTree(temp);
        Enu_mc_truth_all_signal_FSIType.push_back(temp);

        temp = new MnvH1D( Form("%s_%d", "QSq_mc_truth_all_signal_FSIType",i),"Q^{2} for Signal Events",binList.size_QSq, binList.a_QSq);
        temp->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBands_TruthTree(temp);
        QSq_mc_truth_all_signal_FSIType.push_back(temp);

        temp = new MnvH1D( Form("%s_%d", "W_mc_truth_all_signal_FSIType",i), "W for Signal Events",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max());
        temp->GetXaxis()->SetTitle("W_{exp} (GeV)");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBands_TruthTree(temp);
        W_mc_truth_all_signal_FSIType.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","deltaInvMass_mc_truth_all_signal_FSIType",i),"",binList.deltaInvMass.get_nBins(), binList.deltaInvMass.get_min(), binList.deltaInvMass.get_max());
        temp->GetXaxis()->SetTitle("deltaInvMass");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBands_TruthTree(temp);
        deltaInvMass_mc_truth_all_signal_FSIType.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","Delta_pi_theta_mc_truth_all_signal_FSIType",i),"",binList.Delta_pi_theta.get_nBins(), binList.Delta_pi_theta.get_min(), binList.Delta_pi_theta.get_max());
        temp->GetXaxis()->SetTitle("Delta_pi_theta");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBands_TruthTree(temp);
        Delta_pi_theta_mc_truth_all_signal_FSIType.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","Delta_pi_phi_mc_truth_all_signal_FSIType",i),"",binList.Delta_pi_phi.get_nBins(), binList.Delta_pi_phi.get_min(), binList.Delta_pi_phi.get_max());
        temp->GetXaxis()->SetTitle("Delta_pi_phi");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBands_TruthTree(temp);
        Delta_pi_phi_mc_truth_all_signal_FSIType.push_back(temp);
    }

    // ------------------------------------------------------------------------
    // Interaction Type 
    // ------------------------------------------------------------------------
    for (int i = 0; i < nIntType; i++){
        temp = new MnvH1D( Form("%s_%d","muon_P_mc_truth_all_signal_IntType",i),"Muon Momentum for Signal Events",binList.size_muon_P, binList.a_muon_P);
        temp->GetXaxis()->SetTitle("Momentum [GeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBands_TruthTree(temp);
        muon_P_mc_truth_all_signal_IntType.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","muon_theta_mc_truth_all_signal_IntType",i),"Muon Theta for Signal Events",binList.size_muon_theta, binList.a_muon_theta);
        temp->GetXaxis()->SetTitle("Theta");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBands_TruthTree(temp);
        muon_theta_mc_truth_all_signal_IntType.push_back(temp);
        
        temp = new MnvH1D( Form("%s_%d","pi0_P_mc_truth_all_signal_IntType",i),"Pi0 Momentum for Signal Events",binList.size_pi0_P, binList.a_pi0_P);
        temp->GetXaxis()->SetTitle("Momentum [GeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBands_TruthTree(temp);
        pi0_P_mc_truth_all_signal_IntType.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","pi0_KE_mc_truth_all_signal_IntType",i),"Pi0 Kinetic Energy for Signal Events",binList.size_pi0_KE, binList.a_pi0_KE);
        temp->GetXaxis()->SetTitle("Kinetic Energy [GeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBands_TruthTree(temp);
        pi0_KE_mc_truth_all_signal_IntType.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","pi0_theta_mc_truth_all_signal_IntType",i), "Pi0 Theta for Signal Events",binList.size_pi0_theta, binList.a_pi0_theta);
        temp->GetXaxis()->SetTitle("Theta");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBands_TruthTree(temp);
        pi0_theta_mc_truth_all_signal_IntType.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","Enu_mc_truth_all_signal_IntType",i),"Neutrino Energy for Signal Events",binList.size_Enu, binList.a_Enu);
        temp->GetXaxis()->SetTitle("Neutrino Energy [GeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBands_TruthTree(temp);
        Enu_mc_truth_all_signal_IntType.push_back(temp);

        temp = new MnvH1D( Form("%s_%d", "QSq_mc_truth_all_signal_IntType",i),"Q^{2} for Signal Events",binList.size_QSq, binList.a_QSq);
        temp->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBands_TruthTree(temp);
        QSq_mc_truth_all_signal_IntType.push_back(temp);

        temp = new MnvH1D( Form("%s_%d", "W_mc_truth_all_signal_IntType",i), "W for Signal Events",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max());
        temp->GetXaxis()->SetTitle("W_{exp} (GeV)");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBands_TruthTree(temp);
        W_mc_truth_all_signal_IntType.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","deltaInvMass_mc_truth_all_signal_IntType",i),"",binList.deltaInvMass.get_nBins(), binList.deltaInvMass.get_min(), binList.deltaInvMass.get_max());
        temp->GetXaxis()->SetTitle("deltaInvMass");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBands_TruthTree(temp);
        deltaInvMass_mc_truth_all_signal_IntType.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","Delta_pi_theta_mc_truth_all_signal_IntType",i),"",binList.Delta_pi_theta.get_nBins(), binList.Delta_pi_theta.get_min(), binList.Delta_pi_theta.get_max());
        temp->GetXaxis()->SetTitle("Delta_pi_theta");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBands_TruthTree(temp);
        Delta_pi_theta_mc_truth_all_signal_IntType.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","Delta_pi_phi_mc_truth_all_signal_IntType",i),"",binList.Delta_pi_phi.get_nBins(), binList.Delta_pi_phi.get_min(), binList.Delta_pi_phi.get_max());
        temp->GetXaxis()->SetTitle("Delta_pi_phi");
        temp->GetYaxis()->SetTitle("N(Events)");
        AddVertErrorBands_TruthTree(temp);
        Delta_pi_phi_mc_truth_all_signal_IntType.push_back(temp);
    }

    // ------------------------------------------------------------------------
    // Signal Q2
    // ------------------------------------------------------------------------
    truth_QSq_QE = new TH1D("truth_QSq_QE","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    truth_QSq_QE->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    truth_QSq_QE->GetYaxis()->SetTitle("N(Events)");

    truth_QSq_RES_1232 = new TH1D("truth_QSq_RES_1232","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    truth_QSq_RES_1232->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    truth_QSq_RES_1232->GetYaxis()->SetTitle("N(Events)");

    truth_QSq_RES_1535 = new TH1D("truth_QSq_RES_1535","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    truth_QSq_RES_1535->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    truth_QSq_RES_1535->GetYaxis()->SetTitle("N(Events)");

    truth_QSq_RES_1520 = new TH1D("truth_QSq_RES_1520","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    truth_QSq_RES_1520->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    truth_QSq_RES_1520->GetYaxis()->SetTitle("N(Events)");

    truth_QSq_RES_Other = new TH1D("truth_QSq_RES_Other","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    truth_QSq_RES_Other->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    truth_QSq_RES_Other->GetYaxis()->SetTitle("N(Events)");

    truth_QSq_DIS = new TH1D("truth_QSq_DIS","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    truth_QSq_DIS->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    truth_QSq_DIS->GetYaxis()->SetTitle("N(Events)");

    truth_QSq_Non_RES = new TH1D("truth_QSq_Non_RES","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    truth_QSq_Non_RES->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    truth_QSq_Non_RES->GetYaxis()->SetTitle("N(Events)");

    truth_QSq_2p2h = new TH1D("truth_QSq_2p2h","Q^{2} for Signal Events",binList.mc_QSq.get_nBins(), binList.mc_QSq.get_min(), binList.mc_QSq.get_max());
    truth_QSq_2p2h->GetXaxis()->SetTitle("Q^{2} [GeV^{2}]");
    truth_QSq_2p2h->GetYaxis()->SetTitle("N(Events)");

    // ------------------------------------------------------------------------
    // Signal incomingE
    // ------------------------------------------------------------------------
    truth_Enu_QE = new TH1D("truth_Enu_QE","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    truth_Enu_QE->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    truth_Enu_QE->GetYaxis()->SetTitle("N(Events)");

    truth_Enu_RES_1232 = new TH1D("truth_Enu_RES_1232","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    truth_Enu_RES_1232->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    truth_Enu_RES_1232->GetYaxis()->SetTitle("N(Events)");

    truth_Enu_RES_1535 = new TH1D("truth_Enu_RES_1535","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    truth_Enu_RES_1535->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    truth_Enu_RES_1535->GetYaxis()->SetTitle("N(Events)");

    truth_Enu_RES_1520 = new TH1D("truth_Enu_RES_1520","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    truth_Enu_RES_1520->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    truth_Enu_RES_1520->GetYaxis()->SetTitle("N(Events)");

    truth_Enu_RES_Other = new TH1D("truth_Enu_RES_Other","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    truth_Enu_RES_Other->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    truth_Enu_RES_Other->GetYaxis()->SetTitle("N(Events)");

    truth_Enu_DIS = new TH1D("truth_Enu_DIS","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    truth_Enu_DIS->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    truth_Enu_DIS->GetYaxis()->SetTitle("N(Events)");

    truth_Enu_Non_RES = new TH1D("truth_Enu_Non_RES","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    truth_Enu_Non_RES->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    truth_Enu_Non_RES->GetYaxis()->SetTitle("N(Events)");

    truth_Enu_2p2h = new TH1D("truth_Enu_2p2h","E_{#nu} for Signal Events",binList.mc_incomingE.get_nBins(), binList.mc_incomingE.get_min(), binList.mc_incomingE.get_max());
    truth_Enu_2p2h->GetXaxis()->SetTitle("E_{#nu} [GeV]");
    truth_Enu_2p2h->GetYaxis()->SetTitle("N(Events)");

    // ------------------------------------------------------------------------
    // Signal w
    // ------------------------------------------------------------------------
    truth_w_QE = new TH1D("truth_w_QE","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    truth_w_QE->GetXaxis()->SetTitle("W_{exp} (GeV)");
    truth_w_QE->GetYaxis()->SetTitle("N(Events)");
 
    truth_w_RES_1232 = new TH1D("truth_w_RES_1232","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    truth_w_RES_1232->GetXaxis()->SetTitle("W_{exp} (GeV)");
    truth_w_RES_1232->GetYaxis()->SetTitle("N(Events)");

    truth_w_RES_1535 = new TH1D("truth_w_RES_1535","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    truth_w_RES_1535->GetXaxis()->SetTitle("W_{exp} (GeV)");
    truth_w_RES_1535->GetYaxis()->SetTitle("N(Events)");

    truth_w_RES_1520 = new TH1D("truth_w_RES_1520","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    truth_w_RES_1520->GetXaxis()->SetTitle("W_{exp} (GeV)");
    truth_w_RES_1520->GetYaxis()->SetTitle("N(Events)");

    truth_w_RES_Other = new TH1D("truth_w_RES_Other","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    truth_w_RES_Other->GetXaxis()->SetTitle("W_{exp} (GeV)");
    truth_w_RES_Other->GetYaxis()->SetTitle("N(Events)");

    truth_w_DIS = new TH1D("truth_w_DIS","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    truth_w_DIS->GetXaxis()->SetTitle("W_{exp} (GeV)");
    truth_w_DIS->GetYaxis()->SetTitle("N(Events)");

    truth_w_Non_RES = new TH1D("truth_w_Non_RES","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    truth_w_Non_RES->GetXaxis()->SetTitle("W_{exp} (GeV)");
    truth_w_Non_RES->GetYaxis()->SetTitle("N(Events)");

    truth_w_2p2h = new TH1D("truth_w_2p2h","W for Signal Events",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max());
    truth_w_2p2h->GetXaxis()->SetTitle("W_{exp} (GeV)");
    truth_w_2p2h->GetYaxis()->SetTitle("N(Events)");

    // Delta RES Study
    delta_anisotropy = new MnvH1D( "delta_anisotropy","Theta_Delta2Npi GENIE 1#sigma Weights", 80,0.0,2.0);
    delta_anisotropy->GetXaxis()->SetTitle("Theta_Delta2Npi GENIE 1#sigma Weights");
    delta_anisotropy->GetYaxis()->SetTitle("N(Events)");
}

int CCProtonPi0_TruthAnalyzer::GetBackgroundTypeInd()
{
    // Check For Signal
    if (truth_isSignal){
        cout<<"WARNING! Signal Event requested Background Ind! - Returning -1"<<endl;
        return -1;
    }

    if (truth_isBckg_Compact_WithPi0) return 3;
    else if (truth_isBckg_Compact_QELike) return 4;
    else if (truth_isBckg_Compact_SinglePiPlus) return 5;
    else if (truth_isBckg_Compact_Other) return 6;
    else{
        cout<<"WARNING! No Background Type Found - Returning -1"<<endl;
        return -1;
    }
}

void CCProtonPi0_TruthAnalyzer::FillHistogram(vector<MnvH1D*> &hist, double var)
{
    // Always Fill hist[0]
    hist[0]->Fill(var, cvweight);

    // Fill Signal
    if (truth_isSignal){
        hist[1]->Fill(var, cvweight);
    }else{
        // Fill Background
        hist[2]->Fill(var, cvweight); // Always Fill ind == 2 -- All Background

        // Fill Background Type
        int ind = GetBackgroundTypeInd();
        hist[ind]->Fill(var, cvweight);
    }
}

void CCProtonPi0_TruthAnalyzer::FillHistogram(MnvH1D* hist, double var)
{
    hist->Fill(var, cvweight);

    if (fillErrors_ByHand){
        FillVertErrorBand_Flux_ByHand(hist, var);
        FillVertErrorBand_Genie_ByHand(hist, var);
        FillVertErrorBand_2p2h_ByHand(hist, var);

        //FillVertErrorBand_HighMaRES_ByHand(hist, var);
        //FillVertErrorBand_LowMaRES_ByHand(hist, var);
        //FillVertErrorBand_DeltaFactor_ByHand(hist, var);
    }else{
        FillVertErrorBand_Flux(hist, var);
        FillVertErrorBand_Genie(hist, var);
        FillVertErrorBand_2p2h(hist, var);
    }
}

void CCProtonPi0_TruthAnalyzer::FillHistogram(TH1D* hist, double var)
{
    hist->Fill(var, cvweight);
}

void CCProtonPi0_TruthAnalyzer::FillVertErrorBand_HighMaRES_ByHand(MnvH1D* h, double var)
{
    std::vector<double> wgts_up = QSqFitter.GetWeights(truth_genie_wgt_MaRES[4], truth_genie_wgt_MaRES[5]);
    FillVertErrorBand_ByHand(h, var, "HighMaRES", wgts_up);
}

void CCProtonPi0_TruthAnalyzer::FillVertErrorBand_LowMaRES_ByHand(MnvH1D* h, double var)
{
    std::vector<double> wgts_dn = QSqFitter.GetWeights(truth_genie_wgt_MaRES[2], truth_genie_wgt_MaRES[1]);
    FillVertErrorBand_ByHand(h, var, "LowMaRES", wgts_dn);
}

void CCProtonPi0_TruthAnalyzer::FillVertErrorBand_DeltaFactor_ByHand(MnvH1D* h, double var)
{
    std::vector<double> wgts = GetDeltaFactorWeights();
    FillVertErrorBand_ByHand(h, var, "DeltaFactor", wgts);
}

std::vector<double> CCProtonPi0_TruthAnalyzer::GetDeltaFactorWeights()
{
    const double A = 1.0;
    const double A_MINOS = 1.010;
    const double Q0_MINOS = 0.156;
    
    // Vary Q0
    std::vector<double> Q0_Vector;
    for (double q0 = 0.050; q0 <= 0.155; q0 += 0.001){
        Q0_Vector.push_back(q0);
    }

    // Fill Delta Factor Weights
    //      First Element is MINOS Factor
    std::vector<double> weights;
    double QSq = truth_QSq_exp * MeVSq_to_GeVSq;
    double deltaFactor = GetDeltaFactor(QSq, A_MINOS, Q0_MINOS);
    weights.push_back(deltaFactor);

    for (unsigned int i = 0; i < Q0_Vector.size(); ++i){
        deltaFactor = GetDeltaFactor(QSq, A, Q0_Vector[i]);
        weights.push_back(deltaFactor);
    }


    // Debugging
//    std::cout<<"Delta Weights Size = "<<weights.size()<<std::endl;
//    for (unsigned int i = 0; i < weights.size(); ++i){
//        std::cout<<weights[i]<<std::endl;
//    }

    return weights;
}

void CCProtonPi0_TruthAnalyzer::FillVertErrorBand_Flux(MnvH1D* h, double var)
{
    std::vector<double> flux_errors = GetFluxError(mc_incomingE * MeV_to_GeV, mc_incoming);
    h->FillVertErrorBand("Flux",  var, &flux_errors[0],  cvweight);
}

void CCProtonPi0_TruthAnalyzer::FillVertErrorBand_Flux_ByHand(MnvH1D* h, double var)
{
    std::vector<double> flux_errors = GetFluxError(mc_incomingE * MeV_to_GeV, mc_incoming);
    FillVertErrorBand_ByHand(h, var, "Flux", flux_errors);
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
    //h->FillVertErrorBand("GENIE_MaCCQEshape"       ,var, truth_genie_wgt_MaCCQEshape[2]      , truth_genie_wgt_MaCCQEshape[4]      , cvweight);
    h->FillVertErrorBand("GENIE_MaNCEL"            ,var, truth_genie_wgt_MaNCEL[2]           , truth_genie_wgt_MaNCEL[4]           , cvweight);
    //h->FillVertErrorBand("GENIE_NormCCQE"          ,var, truth_genie_wgt_NormCCQE[2]         , truth_genie_wgt_NormCCQE[4]         , cvweight);
    h->FillVertErrorBand("GENIE_NormDISCC"         ,var, truth_genie_wgt_NormDISCC[2]        , truth_genie_wgt_NormDISCC[4]        , cvweight);
    h->FillVertErrorBand("GENIE_NormNCRES"         ,var, truth_genie_wgt_NormNCRES[2]        , truth_genie_wgt_NormNCRES[4]        , cvweight);
    h->FillVertErrorBand("GENIE_RDecBR1gamma"      ,var, truth_genie_wgt_RDecBR1gamma[2]     , truth_genie_wgt_RDecBR1gamma[4]     , cvweight);
    h->FillVertErrorBand("GENIE_Rvn2pi"            ,var, truth_genie_wgt_Rvn2pi[2]           , truth_genie_wgt_Rvn2pi[4]           , cvweight);
    h->FillVertErrorBand("GENIE_Rvp2pi"            ,var, truth_genie_wgt_Rvp2pi[2]           , truth_genie_wgt_Rvp2pi[4]           , cvweight);
    h->FillVertErrorBand("GENIE_VecFFCCQEshape"    ,var, truth_genie_wgt_VecFFCCQEshape[2]   , truth_genie_wgt_VecFFCCQEshape[4]   , cvweight);

    h->FillVertErrorBand("GENIE_NormCCRES"         ,var, updated_genie_wgt_NormCCRES[2]        , updated_genie_wgt_NormCCRES[4]        , cvweight);
    h->FillVertErrorBand("GENIE_Theta_Delta2Npi"   ,var, updated_genie_wgt_Theta_Delta2Npi[2]  , updated_genie_wgt_Theta_Delta2Npi[4]  , cvweight);
    h->FillVertErrorBand("GENIE_MaRES"             ,var, updated_genie_wgt_MaRES[2]            , updated_genie_wgt_MaRES[4]            , cvweight);
    h->FillVertErrorBand("GENIE_MvRES"             ,var, updated_genie_wgt_MvRES[2]            , updated_genie_wgt_MvRES[4]            , cvweight);
    h->FillVertErrorBand("GENIE_Rvn1pi"            ,var, updated_genie_wgt_Rvn1pi[2]           , updated_genie_wgt_Rvn1pi[4]           , cvweight);
    h->FillVertErrorBand("GENIE_Rvp1pi"            ,var, updated_genie_wgt_Rvp1pi[2]           , updated_genie_wgt_Rvp1pi[4]           , cvweight);
}

void CCProtonPi0_TruthAnalyzer::FillVertErrorBand_Genie_ByHand(MnvH1D* h, double var)
{
    FillVertErrorBand_ByHand(h, var, "GENIE_AGKYxF1pi"         , truth_genie_wgt_AGKYxF1pi[2]        , truth_genie_wgt_AGKYxF1pi[4]        );
    FillVertErrorBand_ByHand(h, var, "GENIE_AhtBY"             , truth_genie_wgt_AhtBY[2]            , truth_genie_wgt_AhtBY[4]            );
    FillVertErrorBand_ByHand(h, var, "GENIE_BhtBY"             , truth_genie_wgt_BhtBY[2]            , truth_genie_wgt_BhtBY[4]            );
    FillVertErrorBand_ByHand(h, var, "GENIE_CCQEPauliSupViaKF" , truth_genie_wgt_CCQEPauliSupViaKF[2], truth_genie_wgt_CCQEPauliSupViaKF[4]);
    FillVertErrorBand_ByHand(h, var, "GENIE_CV1uBY"            , truth_genie_wgt_CV1uBY[2]           , truth_genie_wgt_CV1uBY[4]           );
    FillVertErrorBand_ByHand(h, var, "GENIE_CV2uBY"            , truth_genie_wgt_CV2uBY[2]           , truth_genie_wgt_CV2uBY[4]           );
    FillVertErrorBand_ByHand(h, var, "GENIE_EtaNCEL"           , truth_genie_wgt_EtaNCEL[2]          , truth_genie_wgt_EtaNCEL[4]          );
    FillVertErrorBand_ByHand(h, var, "GENIE_FrAbs_N"           , truth_genie_wgt_FrAbs_N[2]          , truth_genie_wgt_FrAbs_N[4]          );
    FillVertErrorBand_ByHand(h, var, "GENIE_FrAbs_pi"          , truth_genie_wgt_FrAbs_pi[2]         , truth_genie_wgt_FrAbs_pi[4]         );
    FillVertErrorBand_ByHand(h, var, "GENIE_FrCEx_N"           , truth_genie_wgt_FrCEx_N[2]          , truth_genie_wgt_FrCEx_N[4]          );
    FillVertErrorBand_ByHand(h, var, "GENIE_FrCEx_pi"          , truth_genie_wgt_FrCEx_pi[2]         , truth_genie_wgt_FrCEx_pi[4]         );
    FillVertErrorBand_ByHand(h, var, "GENIE_FrElas_N"          , truth_genie_wgt_FrElas_N[2]         , truth_genie_wgt_FrElas_N[4]         );
    FillVertErrorBand_ByHand(h, var, "GENIE_FrElas_pi"         , truth_genie_wgt_FrElas_pi[2]        , truth_genie_wgt_FrElas_pi[4]        );
    FillVertErrorBand_ByHand(h, var, "GENIE_FrInel_N"          , truth_genie_wgt_FrInel_N[2]         , truth_genie_wgt_FrInel_N[4]         );
    FillVertErrorBand_ByHand(h, var, "GENIE_FrInel_pi"         , truth_genie_wgt_FrInel_pi[2]        , truth_genie_wgt_FrInel_pi[4]        );
    FillVertErrorBand_ByHand(h, var, "GENIE_FrPiProd_N"        , truth_genie_wgt_FrPiProd_N[2]       , truth_genie_wgt_FrPiProd_N[4]       );
    FillVertErrorBand_ByHand(h, var, "GENIE_FrPiProd_pi"       , truth_genie_wgt_FrPiProd_pi[2]      , truth_genie_wgt_FrPiProd_pi[4]      );
    FillVertErrorBand_ByHand(h, var, "GENIE_MFP_N"             , truth_genie_wgt_MFP_N[2]            , truth_genie_wgt_MFP_N[4]            );
    FillVertErrorBand_ByHand(h, var, "GENIE_MFP_pi"            , truth_genie_wgt_MFP_pi[2]           , truth_genie_wgt_MFP_pi[4]           );
    FillVertErrorBand_ByHand(h, var, "GENIE_MaCCQE"            , truth_genie_wgt_MaCCQE[2]           , truth_genie_wgt_MaCCQE[4]           );
    //FillVertErrorBand_ByHand(h, var, "GENIE_MaCCQEshape"       , truth_genie_wgt_MaCCQEshape[2]      , truth_genie_wgt_MaCCQEshape[4]      );
    FillVertErrorBand_ByHand(h, var, "GENIE_MaNCEL"            , truth_genie_wgt_MaNCEL[2]           , truth_genie_wgt_MaNCEL[4]           );
    //FillVertErrorBand_ByHand(h, var, "GENIE_NormCCQE"          , truth_genie_wgt_NormCCQE[2]         , truth_genie_wgt_NormCCQE[4]         );
    FillVertErrorBand_ByHand(h, var, "GENIE_NormDISCC"         , truth_genie_wgt_NormDISCC[2]        , truth_genie_wgt_NormDISCC[4]        );
    FillVertErrorBand_ByHand(h, var, "GENIE_NormNCRES"         , truth_genie_wgt_NormNCRES[2]        , truth_genie_wgt_NormNCRES[4]        );
    FillVertErrorBand_ByHand(h, var, "GENIE_RDecBR1gamma"      , truth_genie_wgt_RDecBR1gamma[2]     , truth_genie_wgt_RDecBR1gamma[4]     );
    FillVertErrorBand_ByHand(h, var, "GENIE_Rvn2pi"            , truth_genie_wgt_Rvn2pi[2]           , truth_genie_wgt_Rvn2pi[4]           );
    FillVertErrorBand_ByHand(h, var, "GENIE_Rvp2pi"            , truth_genie_wgt_Rvp2pi[2]           , truth_genie_wgt_Rvp2pi[4]           );
    FillVertErrorBand_ByHand(h, var, "GENIE_VecFFCCQEshape"    , truth_genie_wgt_VecFFCCQEshape[2]   , truth_genie_wgt_VecFFCCQEshape[4]   );

    FillVertErrorBand_ByHand(h, var, "GENIE_NormCCRES"         , updated_genie_wgt_NormCCRES[2]        , updated_genie_wgt_NormCCRES[4]        );
    FillVertErrorBand_ByHand(h, var, "GENIE_Theta_Delta2Npi"   , updated_genie_wgt_Theta_Delta2Npi[2]  , updated_genie_wgt_Theta_Delta2Npi[4]  );
    FillVertErrorBand_ByHand(h, var, "GENIE_MaRES"             , updated_genie_wgt_MaRES[2]            , updated_genie_wgt_MaRES[4]            );
    FillVertErrorBand_ByHand(h, var, "GENIE_MvRES"             , updated_genie_wgt_MvRES[2]            , updated_genie_wgt_MvRES[4]            );
    FillVertErrorBand_ByHand(h, var, "GENIE_Rvn1pi"            , updated_genie_wgt_Rvn1pi[2]           , updated_genie_wgt_Rvn1pi[4]           );
    FillVertErrorBand_ByHand(h, var, "GENIE_Rvp1pi"            , updated_genie_wgt_Rvp1pi[2]           , updated_genie_wgt_Rvp1pi[4]           );
}

void CCProtonPi0_TruthAnalyzer::FillVertErrorBand_ByHand(MnvH1D* h, double var, std::string error_name, double err_down, double err_up)
{
    std::vector<double> errors;
    errors.push_back(err_down);
    errors.push_back(err_up);

    FillVertErrorBand_ByHand(h, var, error_name, errors);
}

void CCProtonPi0_TruthAnalyzer::FillVertErrorBand_ByHand(MnvH1D* h, double var, std::string error_name, std::vector<double> errors)
{
    // Get a Pointer to Error Band
    MnvVertErrorBand* err_band =  h->GetVertErrorBand(error_name);

    // Get a Pointer to Histograms 
    std::vector<TH1D*> err_hists = err_band->GetHists();
   
    // Sanity Check
    if (err_hists.size() != errors.size()) {
        std::cout<<"WARNING! Can not Fill Vertical Error Band: "<<error_name<<std::endl;
        exit(1);
    }

    // Fill Error Band Base Histogram with Default cvweight
    int cvbin = err_band->TH1D::Fill(var, cvweight);

    // Update Errors
    if (cvbin == -1){
        cvbin = err_band->FindBin(var);
    }

    for (unsigned int i = 0; i < err_hists.size(); ++i ) {
        
        // wgt_bckg is  1.0 -- there is no Background
        double wgt_bckg = 1.0;
      
        const double applyWeight = cvweight * wgt_bckg;
        const double wgtU = errors[i]*applyWeight;
        err_hists[i]->AddBinContent( cvbin, wgtU );

        const double err = err_hists[i]->GetBinError(cvbin);
        const double newerr2 = err*err + wgtU*wgtU; 
        const double newerr = (0.<newerr2) ? sqrt(newerr2) : 0.;
        err_hists[i]->SetBinError( cvbin, newerr );
    }
}

void CCProtonPi0_TruthAnalyzer::CalcEventWeight()
{
    UpdateFluxReweighter(mc_run, mc_intType); 

    // Reset cvweight
    cvweight = 1.0;
    cvweight_2p2h = 1.0;
    cvweight_Delta = 1.0;
    cvweight_CCRES = 1.0;
    cvweight_NonRes1pi = 1.0;

    // Replace cvweight with Flux Weight
    cvweight = GetFluxWeight(mc_incomingE * MeV_to_GeV, mc_incoming);

    // 2p2h Event Weights
    if ( IsEvent2p2h(mc_intType) ){
        cvweight_2p2h = Get_2p2h_wgt(mc_incomingPartVec, mc_primFSLepton, fit_2p2h_CV) * POT_ratio_2p2h;
        cvweight *= cvweight_2p2h;
        // Calc 2p2h Uncertainty 
        Get2p2hErr();
    }else{
        // 2p2h Uncertainty is 0.0 for other events 
        err_2p2h = 0.0;
    }

    if (applyGENIETuning_Complete){

        // Delta decay anisotropy weighting per DocDB 9850.  Weight is 1.0 for non-delta resonance interactions.
        cvweight_Delta *= ( 1.0 + truth_genie_wgt_Theta_Delta2Npi[4] ) / 2.0;
        cvweight *= cvweight_Delta;

        if ( IsGenieRvn1pi() || IsGenieRvp1pi() ){
            cvweight_NonRes1pi *= deuteriumNonResNorm;
            cvweight *= cvweight_NonRes1pi;
        }

        if ( IsGenieCCRes() ){
            cvweight_CCRES *= deuteriumResNorm; 
            cvweight_CCRES *= GetMaResWeight(deuteriumResNorm);
            cvweight *= cvweight_CCRES;
        }
    }

    if ( applyGENIETuning_DeltaSuppression && IsGenieCCRes() ){
        // Delta Suppression Factor
        double QSq = truth_QSq_exp * MeVSq_to_GeVSq;
        cvweight_CCRES *= GetDeltaFactor(QSq, DeltaFactor_A, DeltaFactor_Q0);
        cvweight *= cvweight_CCRES;
    }

    UpdateGENIESystematics();
}

void CCProtonPi0_TruthAnalyzer::AddOtherErrorBands_FillWithCV()
{
    // Add Vertical Error Bands
    AddVertErrorBandsAndFillWithCV_TruthTree(muon_P_mc_truth_all_signal);
    AddVertErrorBandsAndFillWithCV_TruthTree(muon_theta_mc_truth_all_signal);
    AddVertErrorBandsAndFillWithCV_TruthTree(pi0_P_mc_truth_all_signal);
    AddVertErrorBandsAndFillWithCV_TruthTree(pi0_KE_mc_truth_all_signal);
    AddVertErrorBandsAndFillWithCV_TruthTree(pi0_theta_mc_truth_all_signal);
    AddVertErrorBandsAndFillWithCV_TruthTree(Enu_mc_truth_all_signal);
    AddVertErrorBandsAndFillWithCV_TruthTree(QSq_mc_truth_all_signal);
    AddVertErrorBandsAndFillWithCV_TruthTree(W_mc_truth_all_signal);
    AddVertErrorBandsAndFillWithCV_TruthTree(deltaInvMass_mc_truth_all_signal);
    AddVertErrorBandsAndFillWithCV_TruthTree(Delta_pi_theta_mc_truth_all_signal);
    AddVertErrorBandsAndFillWithCV_TruthTree(Delta_pi_phi_mc_truth_all_signal);

    // Add Lateral Error Bands
    AddLatErrorBandsAndFillWithCV_TruthTree(muon_P_mc_truth_all_signal);
    AddLatErrorBandsAndFillWithCV_TruthTree(muon_theta_mc_truth_all_signal);
    AddLatErrorBandsAndFillWithCV_TruthTree(pi0_P_mc_truth_all_signal);
    AddLatErrorBandsAndFillWithCV_TruthTree(pi0_KE_mc_truth_all_signal);
    AddLatErrorBandsAndFillWithCV_TruthTree(pi0_theta_mc_truth_all_signal);
    AddLatErrorBandsAndFillWithCV_TruthTree(Enu_mc_truth_all_signal);
    AddLatErrorBandsAndFillWithCV_TruthTree(QSq_mc_truth_all_signal);
    AddLatErrorBandsAndFillWithCV_TruthTree(W_mc_truth_all_signal);
    AddLatErrorBandsAndFillWithCV_TruthTree(deltaInvMass_mc_truth_all_signal);
    AddLatErrorBandsAndFillWithCV_TruthTree(Delta_pi_theta_mc_truth_all_signal);
    AddLatErrorBandsAndFillWithCV_TruthTree(Delta_pi_phi_mc_truth_all_signal);

    // Add Vertical Error Bands
    AddVertErrorBandsAndFillWithCV_TruthTree(muon_P_mc_truth_all_signal_BeforeFSI);
    AddVertErrorBandsAndFillWithCV_TruthTree(muon_theta_mc_truth_all_signal_BeforeFSI);
    AddVertErrorBandsAndFillWithCV_TruthTree(pi0_P_mc_truth_all_signal_BeforeFSI);
    AddVertErrorBandsAndFillWithCV_TruthTree(pi0_KE_mc_truth_all_signal_BeforeFSI);
    AddVertErrorBandsAndFillWithCV_TruthTree(pi0_theta_mc_truth_all_signal_BeforeFSI);
    AddVertErrorBandsAndFillWithCV_TruthTree(Enu_mc_truth_all_signal_BeforeFSI);
    AddVertErrorBandsAndFillWithCV_TruthTree(QSq_mc_truth_all_signal_BeforeFSI);
    AddVertErrorBandsAndFillWithCV_TruthTree(W_mc_truth_all_signal_BeforeFSI);
    AddVertErrorBandsAndFillWithCV_TruthTree(deltaInvMass_mc_truth_all_signal_BeforeFSI);
    AddVertErrorBandsAndFillWithCV_TruthTree(Delta_pi_theta_mc_truth_all_signal_BeforeFSI);
    AddVertErrorBandsAndFillWithCV_TruthTree(Delta_pi_phi_mc_truth_all_signal_BeforeFSI);

    // Add Lateral Error Bands
    AddLatErrorBandsAndFillWithCV_TruthTree(muon_P_mc_truth_all_signal_BeforeFSI);
    AddLatErrorBandsAndFillWithCV_TruthTree(muon_theta_mc_truth_all_signal_BeforeFSI);
    AddLatErrorBandsAndFillWithCV_TruthTree(pi0_P_mc_truth_all_signal_BeforeFSI);
    AddLatErrorBandsAndFillWithCV_TruthTree(pi0_KE_mc_truth_all_signal_BeforeFSI);
    AddLatErrorBandsAndFillWithCV_TruthTree(pi0_theta_mc_truth_all_signal_BeforeFSI);
    AddLatErrorBandsAndFillWithCV_TruthTree(Enu_mc_truth_all_signal_BeforeFSI);
    AddLatErrorBandsAndFillWithCV_TruthTree(QSq_mc_truth_all_signal_BeforeFSI);
    AddLatErrorBandsAndFillWithCV_TruthTree(W_mc_truth_all_signal_BeforeFSI);
    AddLatErrorBandsAndFillWithCV_TruthTree(deltaInvMass_mc_truth_all_signal_BeforeFSI);
    AddLatErrorBandsAndFillWithCV_TruthTree(Delta_pi_theta_mc_truth_all_signal_BeforeFSI);
    AddLatErrorBandsAndFillWithCV_TruthTree(Delta_pi_phi_mc_truth_all_signal_BeforeFSI);

    // Add Vertical Error Bands
    AddVertErrorBandsAndFillWithCV_TruthTree(muon_P_mc_truth_all_signal_AfterFSI);
    AddVertErrorBandsAndFillWithCV_TruthTree(muon_theta_mc_truth_all_signal_AfterFSI);
    AddVertErrorBandsAndFillWithCV_TruthTree(pi0_P_mc_truth_all_signal_AfterFSI);
    AddVertErrorBandsAndFillWithCV_TruthTree(pi0_KE_mc_truth_all_signal_AfterFSI);
    AddVertErrorBandsAndFillWithCV_TruthTree(pi0_theta_mc_truth_all_signal_AfterFSI);
    AddVertErrorBandsAndFillWithCV_TruthTree(Enu_mc_truth_all_signal_AfterFSI);
    AddVertErrorBandsAndFillWithCV_TruthTree(QSq_mc_truth_all_signal_AfterFSI);
    AddVertErrorBandsAndFillWithCV_TruthTree(W_mc_truth_all_signal_AfterFSI);
    AddVertErrorBandsAndFillWithCV_TruthTree(deltaInvMass_mc_truth_all_signal_AfterFSI);
    AddVertErrorBandsAndFillWithCV_TruthTree(Delta_pi_theta_mc_truth_all_signal_AfterFSI);
    AddVertErrorBandsAndFillWithCV_TruthTree(Delta_pi_phi_mc_truth_all_signal_AfterFSI);

    // Add Lateral Error Bands
    AddLatErrorBandsAndFillWithCV_TruthTree(muon_P_mc_truth_all_signal_AfterFSI);
    AddLatErrorBandsAndFillWithCV_TruthTree(muon_theta_mc_truth_all_signal_AfterFSI);
    AddLatErrorBandsAndFillWithCV_TruthTree(pi0_P_mc_truth_all_signal_AfterFSI);
    AddLatErrorBandsAndFillWithCV_TruthTree(pi0_KE_mc_truth_all_signal_AfterFSI);
    AddLatErrorBandsAndFillWithCV_TruthTree(pi0_theta_mc_truth_all_signal_AfterFSI);
    AddLatErrorBandsAndFillWithCV_TruthTree(Enu_mc_truth_all_signal_AfterFSI);
    AddLatErrorBandsAndFillWithCV_TruthTree(QSq_mc_truth_all_signal_AfterFSI);
    AddLatErrorBandsAndFillWithCV_TruthTree(W_mc_truth_all_signal_AfterFSI);
    AddLatErrorBandsAndFillWithCV_TruthTree(deltaInvMass_mc_truth_all_signal_AfterFSI);
    AddLatErrorBandsAndFillWithCV_TruthTree(Delta_pi_theta_mc_truth_all_signal_AfterFSI);
    AddLatErrorBandsAndFillWithCV_TruthTree(Delta_pi_phi_mc_truth_all_signal_AfterFSI);

    for (int i = 0; i < nFSIType; ++i){
        // Add Vertical Error Bands
        AddVertErrorBandsAndFillWithCV_TruthTree(muon_P_mc_truth_all_signal_FSIType[i]);
        AddVertErrorBandsAndFillWithCV_TruthTree(muon_theta_mc_truth_all_signal_FSIType[i]);
        AddVertErrorBandsAndFillWithCV_TruthTree(pi0_P_mc_truth_all_signal_FSIType[i]);
        AddVertErrorBandsAndFillWithCV_TruthTree(pi0_KE_mc_truth_all_signal_FSIType[i]);
        AddVertErrorBandsAndFillWithCV_TruthTree(pi0_theta_mc_truth_all_signal_FSIType[i]);
        AddVertErrorBandsAndFillWithCV_TruthTree(Enu_mc_truth_all_signal_FSIType[i]);
        AddVertErrorBandsAndFillWithCV_TruthTree(QSq_mc_truth_all_signal_FSIType[i]);
        AddVertErrorBandsAndFillWithCV_TruthTree(W_mc_truth_all_signal_FSIType[i]);
        AddVertErrorBandsAndFillWithCV_TruthTree(deltaInvMass_mc_truth_all_signal_FSIType[i]);
        AddVertErrorBandsAndFillWithCV_TruthTree(Delta_pi_theta_mc_truth_all_signal_FSIType[i]);
        AddVertErrorBandsAndFillWithCV_TruthTree(Delta_pi_phi_mc_truth_all_signal_FSIType[i]);

        // Add Lateral Error Bands
        AddLatErrorBandsAndFillWithCV_TruthTree(muon_P_mc_truth_all_signal_FSIType[i]);
        AddLatErrorBandsAndFillWithCV_TruthTree(muon_theta_mc_truth_all_signal_FSIType[i]);
        AddLatErrorBandsAndFillWithCV_TruthTree(pi0_P_mc_truth_all_signal_FSIType[i]);
        AddLatErrorBandsAndFillWithCV_TruthTree(pi0_KE_mc_truth_all_signal_FSIType[i]);
        AddLatErrorBandsAndFillWithCV_TruthTree(pi0_theta_mc_truth_all_signal_FSIType[i]);
        AddLatErrorBandsAndFillWithCV_TruthTree(Enu_mc_truth_all_signal_FSIType[i]);
        AddLatErrorBandsAndFillWithCV_TruthTree(QSq_mc_truth_all_signal_FSIType[i]);
        AddLatErrorBandsAndFillWithCV_TruthTree(W_mc_truth_all_signal_FSIType[i]);
        AddLatErrorBandsAndFillWithCV_TruthTree(deltaInvMass_mc_truth_all_signal_FSIType[i]);
        AddLatErrorBandsAndFillWithCV_TruthTree(Delta_pi_theta_mc_truth_all_signal_FSIType[i]);
        AddLatErrorBandsAndFillWithCV_TruthTree(Delta_pi_phi_mc_truth_all_signal_FSIType[i]);
    }

    for (int i = 0; i < nIntType; ++i){
        // Add Vertical Error Bands
        AddVertErrorBandsAndFillWithCV_TruthTree(muon_P_mc_truth_all_signal_IntType[i]);
        AddVertErrorBandsAndFillWithCV_TruthTree(muon_theta_mc_truth_all_signal_IntType[i]);
        AddVertErrorBandsAndFillWithCV_TruthTree(pi0_P_mc_truth_all_signal_IntType[i]);
        AddVertErrorBandsAndFillWithCV_TruthTree(pi0_KE_mc_truth_all_signal_IntType[i]);
        AddVertErrorBandsAndFillWithCV_TruthTree(pi0_theta_mc_truth_all_signal_IntType[i]);
        AddVertErrorBandsAndFillWithCV_TruthTree(Enu_mc_truth_all_signal_IntType[i]);
        AddVertErrorBandsAndFillWithCV_TruthTree(QSq_mc_truth_all_signal_IntType[i]);
        AddVertErrorBandsAndFillWithCV_TruthTree(W_mc_truth_all_signal_IntType[i]);
        AddVertErrorBandsAndFillWithCV_TruthTree(deltaInvMass_mc_truth_all_signal_IntType[i]);
        AddVertErrorBandsAndFillWithCV_TruthTree(Delta_pi_theta_mc_truth_all_signal_IntType[i]);
        AddVertErrorBandsAndFillWithCV_TruthTree(Delta_pi_phi_mc_truth_all_signal_IntType[i]);

        // Add Lateral Error Bands
        AddLatErrorBandsAndFillWithCV_TruthTree(muon_P_mc_truth_all_signal_IntType[i]);
        AddLatErrorBandsAndFillWithCV_TruthTree(muon_theta_mc_truth_all_signal_IntType[i]);
        AddLatErrorBandsAndFillWithCV_TruthTree(pi0_P_mc_truth_all_signal_IntType[i]);
        AddLatErrorBandsAndFillWithCV_TruthTree(pi0_KE_mc_truth_all_signal_IntType[i]);
        AddLatErrorBandsAndFillWithCV_TruthTree(pi0_theta_mc_truth_all_signal_IntType[i]);
        AddLatErrorBandsAndFillWithCV_TruthTree(Enu_mc_truth_all_signal_IntType[i]);
        AddLatErrorBandsAndFillWithCV_TruthTree(QSq_mc_truth_all_signal_IntType[i]);
        AddLatErrorBandsAndFillWithCV_TruthTree(W_mc_truth_all_signal_IntType[i]);
        AddLatErrorBandsAndFillWithCV_TruthTree(deltaInvMass_mc_truth_all_signal_IntType[i]);
        AddLatErrorBandsAndFillWithCV_TruthTree(Delta_pi_theta_mc_truth_all_signal_IntType[i]);
        AddLatErrorBandsAndFillWithCV_TruthTree(Delta_pi_phi_mc_truth_all_signal_IntType[i]);
    }

}

void CCProtonPi0_TruthAnalyzer::FillSignal_Test()
{
    for (int i = 0; i < mc_er_nPart; ++i){
        if (mc_er_ID[i] != 111) continue;
        logFile<<std::left;
        logFile.width(8); logFile<<mc_er_ID[i]<<" "; 
        logFile.width(8); logFile<<mc_er_status[i]<<" "; 
        logFile.width(8); logFile<<mc_er_mother[i]<<" ";    
        logFile.width(12); logFile<<mc_er_E[i]<<" ";      
        logFile.width(12); logFile<<truth_pi0_4P[3]<<" ";      
        logFile<<std::endl;
    }
    logFile<<std::endl;
    logFile<<std::endl;
}

void CCProtonPi0_TruthAnalyzer::FillSignal_XSec_Variables_IntType(int type)
{
    FillHistogram(muon_P_mc_truth_all_signal_IntType[type], truth_muon_P * MeV_to_GeV);
    FillHistogram(muon_theta_mc_truth_all_signal_IntType[type], truth_muon_theta_beam * rad_to_deg);
    FillHistogram(pi0_P_mc_truth_all_signal_IntType[type], truth_pi0_P * MeV_to_GeV);
    FillHistogram(pi0_KE_mc_truth_all_signal_IntType[type], truth_pi0_KE * MeV_to_GeV);
    FillHistogram(pi0_theta_mc_truth_all_signal_IntType[type], truth_pi0_theta_beam * rad_to_deg);
    FillHistogram(Enu_mc_truth_all_signal_IntType[type], mc_incomingE * MeV_to_GeV);
    FillHistogram(QSq_mc_truth_all_signal_IntType[type], truth_QSq_exp * MeVSq_to_GeVSq);
    FillHistogram(W_mc_truth_all_signal_IntType[type], truth_W_exp * MeV_to_GeV);
    if (isDeltaRichSample()){
        FillHistogram(deltaInvMass_mc_truth_all_signal_IntType[type], m_deltaInvMass);
        FillHistogram(Delta_pi_theta_mc_truth_all_signal_IntType[type], m_Delta_pi_theta);
        FillHistogram(Delta_pi_phi_mc_truth_all_signal_IntType[type], m_Delta_pi_phi);
    }
}

void CCProtonPi0_TruthAnalyzer::FillSignal_XSec_Variables_FSIType(int type)
{
    FillHistogram(muon_P_mc_truth_all_signal_FSIType[type], truth_muon_P * MeV_to_GeV);
    FillHistogram(muon_theta_mc_truth_all_signal_FSIType[type], truth_muon_theta_beam * rad_to_deg);
    FillHistogram(pi0_P_mc_truth_all_signal_FSIType[type], truth_pi0_P * MeV_to_GeV);
    FillHistogram(pi0_KE_mc_truth_all_signal_FSIType[type], truth_pi0_KE * MeV_to_GeV);
    FillHistogram(pi0_theta_mc_truth_all_signal_FSIType[type], truth_pi0_theta_beam * rad_to_deg);
    FillHistogram(Enu_mc_truth_all_signal_FSIType[type], mc_incomingE * MeV_to_GeV);
    FillHistogram(QSq_mc_truth_all_signal_FSIType[type], truth_QSq_exp * MeVSq_to_GeVSq);
    FillHistogram(W_mc_truth_all_signal_FSIType[type], truth_W_exp * MeV_to_GeV);
    if (isDeltaRichSample()){
        FillHistogram(deltaInvMass_mc_truth_all_signal_FSIType[type], m_deltaInvMass);
        FillHistogram(Delta_pi_theta_mc_truth_all_signal_FSIType[type], m_Delta_pi_theta);
        FillHistogram(Delta_pi_phi_mc_truth_all_signal_FSIType[type], m_Delta_pi_phi);
    }
}

void CCProtonPi0_TruthAnalyzer::FillSignal_XSec_Variables_BeforeFSI()
{
    FillHistogram(muon_P_mc_truth_all_signal_BeforeFSI, truth_muon_P_BeforeFSI * MeV_to_GeV);
    FillHistogram(muon_theta_mc_truth_all_signal_BeforeFSI, truth_muon_theta_beam_BeforeFSI * rad_to_deg);
    FillHistogram(pi0_P_mc_truth_all_signal_BeforeFSI, truth_pi0_P_BeforeFSI * MeV_to_GeV);
    FillHistogram(pi0_KE_mc_truth_all_signal_BeforeFSI, truth_pi0_KE_BeforeFSI * MeV_to_GeV);
    FillHistogram(pi0_theta_mc_truth_all_signal_BeforeFSI, truth_pi0_theta_beam_BeforeFSI * rad_to_deg);
    FillHistogram(QSq_mc_truth_all_signal_BeforeFSI, truth_QSq_exp_BeforeFSI * MeVSq_to_GeVSq);
    FillHistogram(W_mc_truth_all_signal_BeforeFSI, truth_W_exp_BeforeFSI * MeV_to_GeV);
    FillHistogram(Enu_mc_truth_all_signal_BeforeFSI, truth_Enu_BeforeFSI * MeV_to_GeV);
    if (truth_W_exp_BeforeFSI > 0 && truth_W_exp_BeforeFSI < 1400){
        FillHistogram(deltaInvMass_mc_truth_all_signal_BeforeFSI, m_deltaInvMass_BeforeFSI);
        FillHistogram(Delta_pi_theta_mc_truth_all_signal_BeforeFSI, m_Delta_pi_theta_BeforeFSI);
        FillHistogram(Delta_pi_phi_mc_truth_all_signal_BeforeFSI, m_Delta_pi_phi_BeforeFSI);
    }
}

void CCProtonPi0_TruthAnalyzer::FillSignal_XSec_Variables()
{
    FillHistogram(muon_P_mc_truth_all_signal, truth_muon_P * MeV_to_GeV);
    FillHistogram(muon_theta_mc_truth_all_signal, truth_muon_theta_beam * rad_to_deg);
    FillHistogram(pi0_P_mc_truth_all_signal, truth_pi0_P * MeV_to_GeV);
    FillHistogram(pi0_KE_mc_truth_all_signal, truth_pi0_KE * MeV_to_GeV);
    FillHistogram(pi0_theta_mc_truth_all_signal, truth_pi0_theta_beam * rad_to_deg);
    FillHistogram(Enu_mc_truth_all_signal, mc_incomingE * MeV_to_GeV);
    FillHistogram(QSq_mc_truth_all_signal, truth_QSq_exp * MeVSq_to_GeVSq);
    FillHistogram(W_mc_truth_all_signal, truth_W_exp * MeV_to_GeV);
    if (isDeltaRichSample()){
        FillHistogram(deltaInvMass_mc_truth_all_signal, m_deltaInvMass);
        FillHistogram(Delta_pi_theta_mc_truth_all_signal, m_Delta_pi_theta);
        FillHistogram(Delta_pi_phi_mc_truth_all_signal, m_Delta_pi_phi);
    }

    // After FSI -- Same XSec with Finer Binning
    FillHistogram(muon_P_mc_truth_all_signal_AfterFSI, truth_muon_P * MeV_to_GeV);
    FillHistogram(muon_theta_mc_truth_all_signal_AfterFSI, truth_muon_theta_beam * rad_to_deg);
    FillHistogram(pi0_P_mc_truth_all_signal_AfterFSI, truth_pi0_P * MeV_to_GeV);
    FillHistogram(pi0_KE_mc_truth_all_signal_AfterFSI, truth_pi0_KE * MeV_to_GeV);
    FillHistogram(pi0_theta_mc_truth_all_signal_AfterFSI, truth_pi0_theta_beam * rad_to_deg);
    FillHistogram(Enu_mc_truth_all_signal_AfterFSI, mc_incomingE * MeV_to_GeV);
    FillHistogram(QSq_mc_truth_all_signal_AfterFSI, truth_QSq_exp * MeVSq_to_GeVSq);
    FillHistogram(W_mc_truth_all_signal_AfterFSI, truth_W_exp * MeV_to_GeV);
    if (isDeltaRichSample()){
        FillHistogram(deltaInvMass_mc_truth_all_signal_AfterFSI, m_deltaInvMass);
        FillHistogram(Delta_pi_theta_mc_truth_all_signal_AfterFSI, m_Delta_pi_theta);
        FillHistogram(Delta_pi_phi_mc_truth_all_signal_AfterFSI, m_Delta_pi_phi);
    }
}

void CCProtonPi0_TruthAnalyzer::FillSignal_InteractionType()
{
    if (mc_intType == 1){
        nQE.increment(cvweight);
        FillHistogram(truth_QSq_QE, truth_QSq_exp * MeVSq_to_GeVSq);
        FillHistogram(truth_Enu_QE, mc_incomingE * MeV_to_GeV);
        FillHistogram(truth_w_QE, truth_W_exp * MeV_to_GeV);
    }else if (mc_intType == 2){
        if (mc_resID == 0){ 
            nRES_1232.increment(cvweight); 
            FillHistogram(truth_QSq_RES_1232, truth_QSq_exp * MeVSq_to_GeVSq);
            FillHistogram(truth_Enu_RES_1232, mc_incomingE * MeV_to_GeV);
            FillHistogram(truth_w_RES_1232, truth_W_exp * MeV_to_GeV);
        }else if (mc_resID == 1){
            nRES_1535.increment(cvweight); 
            FillHistogram(truth_QSq_RES_1535, truth_QSq_exp * MeVSq_to_GeVSq);
            FillHistogram(truth_Enu_RES_1535, mc_incomingE * MeV_to_GeV);
            FillHistogram(truth_w_RES_1535, truth_W_exp * MeV_to_GeV);
        }else if (mc_resID == 2){
            nRES_1520.increment(cvweight); 
            FillHistogram(truth_QSq_RES_1520, truth_QSq_exp * MeVSq_to_GeVSq);
            FillHistogram(truth_Enu_RES_1520, mc_incomingE * MeV_to_GeV);
            FillHistogram(truth_w_RES_1520, truth_W_exp * MeV_to_GeV);
        }else{
            nRES_Other.increment(cvweight); 
            FillHistogram(truth_QSq_RES_Other, truth_QSq_exp * MeVSq_to_GeVSq);
            FillHistogram(truth_Enu_RES_Other, mc_incomingE * MeV_to_GeV);
            FillHistogram(truth_w_RES_Other, truth_W_exp * MeV_to_GeV);
        }
    }else if ( IsGenieNonRES() ){
        nNon_RES.increment(cvweight);
        FillHistogram(truth_QSq_Non_RES, truth_QSq_exp * MeVSq_to_GeVSq);
        FillHistogram(truth_Enu_Non_RES, mc_incomingE * MeV_to_GeV);
        FillHistogram(truth_w_Non_RES, truth_W_exp * MeV_to_GeV);
    }else if (mc_intType == 3){
        nDIS.increment(cvweight);
        FillHistogram(truth_QSq_DIS, truth_QSq_exp * MeVSq_to_GeVSq);
        FillHistogram(truth_Enu_DIS, mc_incomingE * MeV_to_GeV);
        FillHistogram(truth_w_DIS, truth_W_exp * MeV_to_GeV);
    }else if ( IsEvent2p2h(mc_intType) ){
        n2p2h.increment(cvweight);
        FillHistogram(truth_QSq_2p2h, truth_QSq_exp * MeVSq_to_GeVSq);
        FillHistogram(truth_Enu_2p2h, mc_incomingE * MeV_to_GeV);
        FillHistogram(truth_w_2p2h, truth_W_exp * MeV_to_GeV);
    }else{
        std::cout<<"WARNING! Signal Event with different interaction Type!"<<std::endl;
    }
}

int CCProtonPi0_TruthAnalyzer::Get_nFS_pions()
{
    int nFS_pions = 0;

    for (int i = 0; i < mc_er_nPart; ++i){
        if (std::abs(mc_er_ID[i]) == 211 || mc_er_ID[i] == 111){
            if (isMother_DIS_Fragment(i)) nFS_pions++;
        }
    }

    return nFS_pions;
}

bool CCProtonPi0_TruthAnalyzer::isMother_DIS_Fragment(int ind)
{
    int mother_ind = mc_er_mother[ind];

    if (mc_er_status[mother_ind] == 12) return true;
    else return false;
}

void CCProtonPi0_TruthAnalyzer::PrintEventRecord()
{
    std::cout<<std::endl;
    std::cout<<std::left;
    std::cout.width(4); std::cout<<"ID"<<" "; 
    std::cout.width(12); std::cout<<"PDG"<<" "; 
    std::cout.width(8); std::cout<<"Status"<<" ";    
    std::cout.width(8); std::cout<<"Mother"<<" ";      
    std::cout.width(12); std::cout<<"Energy"<<" ";
    std::cout<<std::endl;
    for (int i = 0; i < mc_er_nPart; ++i){
        std::cout.width(4); std::cout<<i<<" "; 
        std::cout.width(12); std::cout<<mc_er_ID[i]<<" "; 
        std::cout.width(8); std::cout<<mc_er_status[i]<<" ";    
        std::cout.width(8); std::cout<<mc_er_mother[i]<<" ";      
        std::cout.width(12); std::cout<<mc_er_E[i]<<" ";
        std::cout<<std::endl;
    }
}

double CCProtonPi0_TruthAnalyzer::GetMaResWeight( double newMaRes )
{
    // Indices 2 and 4 correspond to -1 and +1 sigma shifts, respectively
    int weightIndex = newMaRes < genieMaRes ? 2 : 4;
    return 1.0 + fabs( newMaRes - genieMaRes ) * ( truth_genie_wgt_MaRES[weightIndex] - 1.0 ) / genieMaRes1sig;
}

double CCProtonPi0_TruthAnalyzer::GetMvResWeight(double newMvRes )
{
  // Indices 2 and 4 correspond to -1 and +1 sigma shifts, respectively
  int weightIndex = newMvRes < genieMvRes ? 2 : 4;
  return 1.0 + fabs( newMvRes - genieMvRes ) * ( truth_genie_wgt_MvRES[weightIndex] - 1.0 ) / genieMvRes1sig;
}

bool CCProtonPi0_TruthAnalyzer::IsGenieCCRes()
{
    return mc_current == 1 && mc_intType == 2;
}

bool CCProtonPi0_TruthAnalyzer::IsGenieNonRES()
{
    // NonRES 1pi 
    if ( IsGenieRvn1pi() || IsGenieRvp1pi() ) return true;

    // NonRES 2pi
    if ( truth_genie_wgt_Rvn2pi[2] != 1 || truth_genie_wgt_Rvp2pi[2] != 1) return true;

    return false;
}

bool CCProtonPi0_TruthAnalyzer::IsGenieRvn1pi()
{
    // Weight is NOT equal to 1 for Rvn1pi events
    return truth_genie_wgt_Rvn1pi[2] != 1;
}

bool CCProtonPi0_TruthAnalyzer::IsGenieRvp1pi()
{
    // Weight is NOT equal to 1 for Rvp1pi events
    return truth_genie_wgt_Rvp1pi[2] != 1;
}

bool CCProtonPi0_TruthAnalyzer::IsGenie_NonRES_n_piplus()
{
    if ( IsGenieRvn1pi() || IsGenieRvp1pi() ){
        std::vector<int> prim_part = GetPrimaryParticles();

        if (prim_part.empty()) return false;

        if (prim_part[0] == 2112 && prim_part[1] == 211) return true;
        else return false;
    }else{
        return false;
    }
}

std::vector<int> CCProtonPi0_TruthAnalyzer::GetPrimaryParticles()
{
    std::vector<int> prim_part;

    for (int i = 0; i < mc_er_nPart; ++i){
        if (mc_er_status[i] == 14){ 
            prim_part.push_back(mc_er_ID[i]);
        }
    }

    return prim_part;
}

void CCProtonPi0_TruthAnalyzer::initUpdatedGenieWeights()
{
    for (int i = 0; i < 7; ++i){
        updated_genie_wgt_Theta_Delta2Npi[i] = truth_genie_wgt_Theta_Delta2Npi[i];
        updated_genie_wgt_NormCCRES[i] = truth_genie_wgt_NormCCRES[i];
        updated_genie_wgt_MaRES[i] = truth_genie_wgt_MaRES[i];
        updated_genie_wgt_MvRES[i] = truth_genie_wgt_MvRES[i];
        updated_genie_wgt_Rvn1pi[i] = truth_genie_wgt_Rvn1pi[i];
        updated_genie_wgt_Rvp1pi[i] = truth_genie_wgt_Rvp1pi[i];
    }
}

void CCProtonPi0_TruthAnalyzer::UpdateGENIESystematics()
{
    // ------------------------------------------------------------------------
    // First init updated_genie_wgts with Nominal Values  
    // ------------------------------------------------------------------------
    initUpdatedGenieWeights();
 
    // ------------------------------------------------------------------------
    // Now Overwrite them with updated values
    // ------------------------------------------------------------------------
    if (applyGENIETuning_Complete){

        // Delta decay non-isotropy weights per DocDB 9850
        updated_genie_wgt_Theta_Delta2Npi[2] = 2.0 / (1.0 + truth_genie_wgt_Theta_Delta2Npi[4]);
        updated_genie_wgt_Theta_Delta2Npi[4] = 2.0 * truth_genie_wgt_Theta_Delta2Npi[4] / ( 1.0 + truth_genie_wgt_Theta_Delta2Npi[4]); 
        if (IsGenieCCRes()){
            updated_genie_wgt_MaRES[2] = GetMaResWeight(deuteriumMaRes - deuteriumMaRes1sig) / GetMaResWeight(deuteriumMaRes);
            updated_genie_wgt_MaRES[4] = GetMaResWeight(deuteriumMaRes + deuteriumMaRes1sig) / GetMaResWeight(deuteriumMaRes);

            updated_genie_wgt_MvRES[2] = GetMvResWeight(genieMvRes - electroProdMvRes1sig); // GENIE MvRES not changed
            updated_genie_wgt_MvRES[4] = GetMvResWeight(genieMvRes + electroProdMvRes1sig); // GENIE MvRES not changed
        
            updated_genie_wgt_NormCCRES[2] = 1.0 - deuteriumResNorm1sig; 
            updated_genie_wgt_NormCCRES[4] = 1.0 + deuteriumResNorm1sig; 
        }

        if (IsGenieRvn1pi() || IsGenieRvp1pi()){
            updated_genie_wgt_Rvn1pi[2] = 1.0 - deuteriumNonResNorm1sig;
            updated_genie_wgt_Rvn1pi[4] = 1.0 + deuteriumNonResNorm1sig;
 
            updated_genie_wgt_Rvp1pi[2] = 1.0 - deuteriumNonResNorm1sig;
            updated_genie_wgt_Rvp1pi[4] = 1.0 + deuteriumNonResNorm1sig;
        }
    }

    bool isDebug = false;
    if (isDebug){
        std::cout<<std::endl;
        std::cout<<"Updated GENIE Systematics"<<std::endl;
        std::cout<<"Theta[2] = "<<truth_genie_wgt_Theta_Delta2Npi[2]<<" "<<updated_genie_wgt_Theta_Delta2Npi[2]<<std::endl;
        std::cout<<"Theta[4] = "<<truth_genie_wgt_Theta_Delta2Npi[4]<<" "<<updated_genie_wgt_Theta_Delta2Npi[4]<<std::endl;
        std::cout<<"NormCCRES[2] = "<<truth_genie_wgt_NormCCRES[2]<<" "<<updated_genie_wgt_NormCCRES[2]<<std::endl;
        std::cout<<"NormCCRES[4] = "<<truth_genie_wgt_NormCCRES[4]<<" "<<updated_genie_wgt_NormCCRES[4]<<std::endl;
        std::cout<<"MaRES[2] = "<<truth_genie_wgt_MaRES[2]<<" "<<updated_genie_wgt_MaRES[2]<<std::endl;
        std::cout<<"MaRES[4] = "<<truth_genie_wgt_MaRES[4]<<" "<<updated_genie_wgt_MaRES[4]<<std::endl;
        std::cout<<"MvRES[2] = "<<truth_genie_wgt_MvRES[2]<<" "<<updated_genie_wgt_MvRES[2]<<std::endl;
        std::cout<<"MvRES[4] = "<<truth_genie_wgt_MvRES[4]<<" "<<updated_genie_wgt_MvRES[4]<<std::endl;
        std::cout<<"Rvn1pi[2] = "<<truth_genie_wgt_Rvn1pi[2]<<" "<<updated_genie_wgt_Rvn1pi[2]<<std::endl;
        std::cout<<"Rvn1pi[4] = "<<truth_genie_wgt_Rvn1pi[4]<<" "<<updated_genie_wgt_Rvn1pi[4]<<std::endl;
        std::cout<<"Rvp1pi[2] = "<<truth_genie_wgt_Rvp1pi[2]<<" "<<updated_genie_wgt_Rvp1pi[2]<<std::endl;
        std::cout<<"Rvp1pi[4] = "<<truth_genie_wgt_Rvp1pi[4]<<" "<<updated_genie_wgt_Rvp1pi[4]<<std::endl;
    }

}

void CCProtonPi0_TruthAnalyzer::FillGENIE_Tuning()
{
    FillHistogram(CV_weight, cvweight);
    FillHistogram(CV_weight_2p2h, cvweight_2p2h);
    FillHistogram(CV_weight_Delta, cvweight_Delta);
    FillHistogram(CV_weight_CCRES, cvweight_CCRES);
    FillHistogram(CV_weight_NonRes1pi, cvweight_NonRes1pi);

    FillHistogram(h_err_2p2h, 1-err_2p2h);
    FillHistogram(h_err_2p2h, 1+err_2p2h);
  
    FillHistogram(genie_wgt_Theta_Delta2Npi, truth_genie_wgt_Theta_Delta2Npi[2]);
    FillHistogram(genie_wgt_Theta_Delta2Npi, truth_genie_wgt_Theta_Delta2Npi[4]);

    FillHistogram(updated_wgt_Theta_Delta2Npi, updated_genie_wgt_Theta_Delta2Npi[2]);
    FillHistogram(updated_wgt_Theta_Delta2Npi, updated_genie_wgt_Theta_Delta2Npi[4]);

    FillHistogram(genie_wgt_MaRES, truth_genie_wgt_MaRES[2]);
    FillHistogram(genie_wgt_MaRES, truth_genie_wgt_MaRES[4]);

    FillHistogram(updated_wgt_MaRES, updated_genie_wgt_MaRES[2]);
    FillHistogram(updated_wgt_MaRES, updated_genie_wgt_MaRES[4]);

    FillHistogram(genie_wgt_MvRES, truth_genie_wgt_MvRES[2]);
    FillHistogram(genie_wgt_MvRES, truth_genie_wgt_MvRES[4]);

    FillHistogram(updated_wgt_MvRES, updated_genie_wgt_MvRES[2]);
    FillHistogram(updated_wgt_MvRES, updated_genie_wgt_MvRES[4]);

    FillHistogram(genie_wgt_Rvn1pi, truth_genie_wgt_Rvn1pi[2]);
    FillHistogram(genie_wgt_Rvn1pi, truth_genie_wgt_Rvn1pi[4]);

    FillHistogram(updated_wgt_Rvn1pi, updated_genie_wgt_Rvn1pi[2]);
    FillHistogram(updated_wgt_Rvn1pi, updated_genie_wgt_Rvn1pi[4]);
}

void CCProtonPi0_TruthAnalyzer::Get2p2hErr()
{
    // Get weights using different fits
    double wgt_cv = Get_2p2h_wgt(mc_incomingPartVec, mc_primFSLepton, fit_2p2h_CV);  
    double wgt_np = Get_2p2h_wgt(mc_incomingPartVec, mc_primFSLepton, fit_2p2h_np);  
    double wgt_nn = Get_2p2h_wgt(mc_incomingPartVec, mc_primFSLepton, fit_2p2h_nn);  

    // Calculate differences
    double err_np = (wgt_np - wgt_cv) / wgt_cv;
    double err_nn = (wgt_nn - wgt_cv) / wgt_cv;

    // Make sure 0 is 0
    err_np = err_np < EPSILON ? 0.0 : err_np;
    err_nn = err_nn < EPSILON ? 0.0 : err_np;

    // Get the average difference
    err_2p2h = (err_np + err_nn) / 2.0;
}

void CCProtonPi0_TruthAnalyzer::FillVertErrorBand_2p2h(MnvH1D* h, double var)
{
    double correctionErr = err_2p2h;
    h->FillVertErrorBand("2p2h", var, 1-correctionErr, 1+correctionErr, cvweight);
}

void CCProtonPi0_TruthAnalyzer::FillVertErrorBand_2p2h_ByHand(MnvH1D* h, double var)
{
    double correctionErr = err_2p2h;
    FillVertErrorBand_ByHand(h, var, "2p2h", 1-correctionErr, 1+correctionErr);
}

double CCProtonPi0_TruthAnalyzer::GetDeltaFactor(double QSq, double A, double Q0)
{
    double factor;

    if (IsGenieCCRes()){
        factor = A / (1 + exp(1-(sqrt(QSq)/Q0)));
    }else{
        factor = 1.0;
    }

    return factor;
}

double CCProtonPi0_TruthAnalyzer::calcDeltaInvariantMass()
{
    TLorentzVector delta_4P = GetDelta4P();
    double M_delta = delta_4P.M(); 

    return M_delta;
}

TLorentzVector CCProtonPi0_TruthAnalyzer::GetDelta4P()
{
    double delta_4P[4];

    for(int i = 0; i < 4; i++){
        delta_4P[i] = truth_proton_4P[i] + truth_pi0_4P[i];
    }

    TLorentzVector delta_lorentz(delta_4P); 

    return delta_lorentz;
}

TLorentzVector CCProtonPi0_TruthAnalyzer::Get_Neutrino_4P(const double Enu)
{
    const double theta = -0.05887;

    double Px = 0.0;
    double Py = Enu*sin(theta);
    double Pz = Enu*cos(theta); 

    TLorentzVector beam_4P(Px,Py,Pz,Enu);

    return beam_4P;
}

void CCProtonPi0_TruthAnalyzer::GetDelta_pi_angles_BeforeFSI()
{
    bool isDebug = false;

    using namespace std;

    // ------------------------------------------------------------------------
    // Get Proton and Pi0 Indice
    // ------------------------------------------------------------------------
    int proton_i = -1;
    int pi0_i = -1;

    for (int i = 0; i < mc_er_nPart; ++i){
        if (mc_er_status[i] != 14) continue;
   
        if (mc_er_ID[i] == 111) pi0_i = i;    
        if (mc_er_ID[i] == 2212) proton_i = i;
    }

    if (proton_i == -1 || pi0_i == -1){
        RunTimeError("Not Signal Before FSI");
    }

    // --------------------------------------------------------------------
    // Form Proton-Pion System
    // --------------------------------------------------------------------
    double delta_4P_v[4];
    delta_4P_v[0] = mc_er_Px[proton_i] + mc_er_Px[pi0_i];
    delta_4P_v[1] = mc_er_Py[proton_i] + mc_er_Py[pi0_i];
    delta_4P_v[2] = mc_er_Pz[proton_i] + mc_er_Pz[pi0_i];
    delta_4P_v[3] = mc_er_E[proton_i] + mc_er_E[pi0_i];

    TLorentzVector delta_4P(delta_4P_v); 
    
    // --------------------------------------------------------------------
    // Find Boost from LAB Frame to Delta Rest Frame
    // --------------------------------------------------------------------
    double M_delta = delta_4P.M(); 
    m_deltaInvMass_BeforeFSI = M_delta * MeV_to_GeV;

    if (isDebug) cout<<"Delta Mass = "<<M_delta * MeV_to_GeV<<endl;

    double gamma = delta_4P.Gamma();
    double boost_x = -(delta_4P.Px()/ (gamma*M_delta));
    double boost_y = -(delta_4P.Py() / (gamma*M_delta));
    double boost_z = -(delta_4P.Pz() / (gamma*M_delta));

    if (isDebug){
        cout<<"Delta P (Before) = "<<delta_4P.Vect().Mag()<<endl;
        cout<<"Delta P (Before) = "<<delta_4P.Vect().x();
        cout<<" "<<delta_4P.Vect().y();
        cout<<" "<<delta_4P.Vect().z()<<endl;
    }

    delta_4P.Boost(boost_x,boost_y,boost_z);

    if (isDebug){
        cout<<"Delta P (After) = "<<delta_4P.Vect().Mag()<<std::endl;
        cout<<"Delta P (After) = "<<delta_4P.Vect().x();
        cout<<" "<<delta_4P.Vect().y();
        cout<<" "<<delta_4P.Vect().z()<<endl;
        cout<<"Delta Mass (After) = "<<delta_4P.M() * MeV_to_GeV<<endl;
    }

    // --------------------------------------------------------------------
    // Boost All 4-Momentums to Delta Rest Frame
    // --------------------------------------------------------------------
    TLorentzVector beam_4P = Get_Neutrino_4P(mc_incomingE);

    TLorentzVector muon_4P(truth_muon_4P[0], truth_muon_4P[1], truth_muon_4P[2], truth_muon_4P[3]); 
    TLorentzVector pi0_4P(mc_er_Px[pi0_i], mc_er_Py[pi0_i], mc_er_Pz[pi0_i], mc_er_E[pi0_i]); 

    if (isDebug){
        cout<<"Before Boost"<<endl;
        cout<<"beam 4P = ("<<beam_4P.Px()<<", "<<beam_4P.Py()<<", "<<beam_4P.Pz()<<", "<<beam_4P.E()<<")"<<endl;
        cout<<"muon 4P = ("<<muon_4P.Px()<<", "<<muon_4P.Py()<<", "<<muon_4P.Pz()<<", "<<muon_4P.E()<<")"<<endl;
        cout<<"Muon Mass = "<<muon_4P.M()<<endl;
        cout<<"pi0 4P = ("<<pi0_4P.Px()<<", "<<pi0_4P.Py()<<", "<<pi0_4P.Pz()<<", "<<pi0_4P.E()<<")"<<endl;
        cout<<"pi0 Mass = "<<pi0_4P.M()<<endl;
    }

    beam_4P.Boost(boost_x, boost_y, boost_z);
    muon_4P.Boost(boost_x, boost_y, boost_z);
    pi0_4P.Boost(boost_x, boost_y, boost_z);

    if (isDebug){
        cout<<"After Boost"<<endl;
        cout<<"beam 4P = ("<<beam_4P.Px()<<", "<<beam_4P.Py()<<", "<<beam_4P.Pz()<<", "<<beam_4P.E()<<")"<<endl;
        cout<<"muon 4P = ("<<muon_4P.Px()<<", "<<muon_4P.Py()<<", "<<muon_4P.Pz()<<", "<<muon_4P.E()<<")"<<endl;
        cout<<"Muon Mass = "<<muon_4P.M()<<endl;
        cout<<"pi0 4P = ("<<pi0_4P.Px()<<", "<<pi0_4P.Py()<<", "<<pi0_4P.Pz()<<", "<<pi0_4P.E()<<")"<<endl;
        cout<<"pi0 Mass = "<<pi0_4P.M()<<endl;
    }

    // --------------------------------------------------------------------
    // Form Axes
    // --------------------------------------------------------------------
    TVector3 beam_3P = beam_4P.Vect();
    TVector3 muon_3P = muon_4P.Vect();

    // Z Axis -- Momentum Transfer Axis
    TVector3 axis_z(beam_3P.Px() - muon_3P.Px(), beam_3P.Py() - muon_3P.Py(), beam_3P.Pz() - muon_3P.Pz());
    TVector3 unit_axis_z = axis_z.Unit();

    // Y Axis -- Beam x Muon
    TVector3 axis_y = beam_3P.Cross(muon_3P);
    TVector3 unit_axis_y = axis_y.Unit();

    // X Axis -- Right Handed Coordinate System y_axis x z_axis = x_axis
    TVector3 axis_x = axis_y.Cross(axis_z);
    TVector3 unit_axis_x = axis_x.Unit();

    // --------------------------------------------------------------------
    // Calculate Angles
    // --------------------------------------------------------------------
    TVector3 pi0_3P = pi0_4P.Vect();
    TVector3 unit_pi0_3P = pi0_3P.Unit();
    double pion_P = pi0_3P.Mag();

    double cos_theta_z = unit_pi0_3P.Dot(unit_axis_z);

    double theta_x = acos(unit_pi0_3P.Dot(unit_axis_x)); 
    double theta_y = acos(unit_pi0_3P.Dot(unit_axis_y)); 
    double theta_z = acos(cos_theta_z); 

    double Px = pion_P*cos(theta_x);
    double Py = pion_P*cos(theta_y);

    double phi_x = acos(Px/(pion_P*sin(theta_z)));
    double phi = phi_x * TMath::RadToDeg();    

    // Convert to 360
    if (Py < 0 ){
        phi = 360 - phi;
    }

    m_Delta_pi_theta_BeforeFSI = cos_theta_z;
    m_Delta_pi_phi_BeforeFSI = phi;

    if (isDebug){
        cout<<"pion_theta = "<<cos_theta_z<<endl;
        cout<<"pion_phi = "<<phi<<endl;
        cout<<endl;
    }
}

void CCProtonPi0_TruthAnalyzer::GetDelta_pi_angles()
{
    bool isDebug = false;

    using namespace std;

    // --------------------------------------------------------------------
    // Find Boost from LAB Frame to Delta Rest Frame
    // --------------------------------------------------------------------
    TLorentzVector delta_4P = GetDelta4P();
    double M_delta = delta_4P.M(); 
    if (isDebug) cout<<"Delta Mass = "<<M_delta * MeV_to_GeV<<endl;

    double gamma = delta_4P.Gamma();
    double boost_x = -(delta_4P.Px()/ (gamma*M_delta));
    double boost_y = -(delta_4P.Py() / (gamma*M_delta));
    double boost_z = -(delta_4P.Pz() / (gamma*M_delta));

    if (isDebug){
        cout<<"Delta P (Before) = "<<delta_4P.Vect().Mag()<<endl;
        cout<<"Delta P (Before) = "<<delta_4P.Vect().x();
        cout<<" "<<delta_4P.Vect().y();
        cout<<" "<<delta_4P.Vect().z()<<endl;
    }

    delta_4P.Boost(boost_x,boost_y,boost_z);

    if (isDebug){
        cout<<"Delta P (After) = "<<delta_4P.Vect().Mag()<<std::endl;
        cout<<"Delta P (After) = "<<delta_4P.Vect().x();
        cout<<" "<<delta_4P.Vect().y();
        cout<<" "<<delta_4P.Vect().z()<<endl;
        cout<<"Delta Mass (After) = "<<delta_4P.M() * MeV_to_GeV<<endl;
    }

    // --------------------------------------------------------------------
    // Boost All 4-Momentums to Delta Rest Frame
    // --------------------------------------------------------------------
    TLorentzVector beam_4P = Get_Neutrino_4P(mc_incomingE);

    TLorentzVector muon_4P(truth_muon_4P[0], truth_muon_4P[1], truth_muon_4P[2], truth_muon_4P[3]); 
    TLorentzVector pi0_4P(truth_pi0_4P[0], truth_pi0_4P[1], truth_pi0_4P[2], truth_pi0_4P[3]); 

    if (isDebug){
        cout<<"Before Boost"<<endl;
        cout<<"beam 4P = ("<<beam_4P.Px()<<", "<<beam_4P.Py()<<", "<<beam_4P.Pz()<<", "<<beam_4P.E()<<")"<<endl;
        cout<<"muon 4P = ("<<muon_4P.Px()<<", "<<muon_4P.Py()<<", "<<muon_4P.Pz()<<", "<<muon_4P.E()<<")"<<endl;
        cout<<"pi0 4P = ("<<pi0_4P.Px()<<", "<<pi0_4P.Py()<<", "<<pi0_4P.Pz()<<", "<<pi0_4P.E()<<")"<<endl;
    }

    beam_4P.Boost(boost_x,boost_y,boost_z);
    muon_4P.Boost(boost_x,boost_y,boost_z);
    pi0_4P.Boost(boost_x,boost_y,boost_z);

    if (isDebug){
        cout<<"After Boost"<<endl;
        cout<<"beam 4P = ("<<beam_4P.Px()<<", "<<beam_4P.Py()<<", "<<beam_4P.Pz()<<", "<<beam_4P.E()<<")"<<endl;
        cout<<"muon 4P = ("<<muon_4P.Px()<<", "<<muon_4P.Py()<<", "<<muon_4P.Pz()<<", "<<muon_4P.E()<<")"<<endl;
        cout<<"pi0 4P = ("<<pi0_4P.Px()<<", "<<pi0_4P.Py()<<", "<<pi0_4P.Pz()<<", "<<pi0_4P.E()<<")"<<endl;
    }

    // --------------------------------------------------------------------
    // Form Axes
    // --------------------------------------------------------------------
    TVector3 beam_3P = beam_4P.Vect();
    TVector3 muon_3P = muon_4P.Vect();

    // Z Axis -- Momentum Transfer Axis
    TVector3 axis_z(beam_3P.Px() - muon_3P.Px(), beam_3P.Py() - muon_3P.Py(), beam_3P.Pz() - muon_3P.Pz());
    TVector3 unit_axis_z = axis_z.Unit();

    // Y Axis -- Beam x Muon
    TVector3 axis_y = beam_3P.Cross(muon_3P);
    TVector3 unit_axis_y = axis_y.Unit();

    // X Axis -- Right Handed Coordinate System y_axis x z_axis = x_axis
    TVector3 axis_x = axis_y.Cross(axis_z);
    TVector3 unit_axis_x = axis_x.Unit();

    // --------------------------------------------------------------------
    // Calculate Angles
    // --------------------------------------------------------------------
    TVector3 pi0_3P = pi0_4P.Vect();
    TVector3 unit_pi0_3P = pi0_3P.Unit();
    double pion_P = pi0_3P.Mag();

    double cos_theta_z = unit_pi0_3P.Dot(unit_axis_z);

    double theta_x = acos(unit_pi0_3P.Dot(unit_axis_x)); 
    double theta_y = acos(unit_pi0_3P.Dot(unit_axis_y)); 
    double theta_z = acos(cos_theta_z); 

    double Px = pion_P*cos(theta_x);
    double Py = pion_P*cos(theta_y);

    double phi_x = acos(Px/(pion_P*sin(theta_z)));
    double phi = phi_x * TMath::RadToDeg();    

    // Convert to 360
    if (Py < 0 ){
        phi = 360 - phi;
    }

    m_Delta_pi_theta = cos_theta_z;
    m_Delta_pi_phi = phi;

    if (isDebug){
        cout<<"pion_theta = "<<m_Delta_pi_theta<<endl;
        cout<<"pion_phi = "<<m_Delta_pi_phi<<endl;
        cout<<endl;
    }
}

bool CCProtonPi0_TruthAnalyzer::isSignal_EnuLow()
{
    if (truth_W_exp > max_W|| mc_incomingE > max_Enu){
        return false;
    }else if (mc_incomingE < min_Enu){
        return true;
    }
    
    std::cout<<"mc_incomingE = "<<mc_incomingE<<" W = "<<truth_W_exp<<std::endl;
    RunTimeError("SignalOut_Kinematics other than W or Enu");
    return false;
    
}

bool CCProtonPi0_TruthAnalyzer::isSignal_WHigh()
{
    if (mc_incomingE < min_Enu || mc_incomingE > max_Enu){
        return false;
    }else if ( truth_W_exp >= max_W){
        return true;
    }
    
    std::cout<<"mc_incomingE = "<<mc_incomingE<<" W = "<<truth_W_exp<<std::endl;
    RunTimeError("SignalOut_Kinematics other than W or Enu");
    return false;
}

void CCProtonPi0_TruthAnalyzer::writeHistograms()
{
    f->cd();

    // Test Variables
    Test_pi0_P->Write();
    Test_pi0_theta->Write();

    // Cross Section Variables
    muon_P_mc_truth_all_signal->Write();
    muon_theta_mc_truth_all_signal->Write();
    pi0_P_mc_truth_all_signal->Write();
    pi0_KE_mc_truth_all_signal->Write();
    pi0_theta_mc_truth_all_signal->Write();
    Enu_mc_truth_all_signal->Write();
    QSq_mc_truth_all_signal->Write();
    W_mc_truth_all_signal->Write();
    deltaInvMass_mc_truth_all_signal->Write();
    Delta_pi_theta_mc_truth_all_signal->Write();
    Delta_pi_phi_mc_truth_all_signal->Write();

    // Cross Section Variables -- Before FSI
    muon_P_mc_truth_all_signal_BeforeFSI->Write();
    muon_theta_mc_truth_all_signal_BeforeFSI->Write();
    pi0_P_mc_truth_all_signal_BeforeFSI->Write();
    pi0_KE_mc_truth_all_signal_BeforeFSI->Write();
    pi0_theta_mc_truth_all_signal_BeforeFSI->Write();
    Enu_mc_truth_all_signal_BeforeFSI->Write();
    QSq_mc_truth_all_signal_BeforeFSI->Write();
    W_mc_truth_all_signal_BeforeFSI->Write();
    deltaInvMass_mc_truth_all_signal_BeforeFSI->Write();
    Delta_pi_theta_mc_truth_all_signal_BeforeFSI->Write();
    Delta_pi_phi_mc_truth_all_signal_BeforeFSI->Write();

    // Cross Section Variables -- After FSI
    muon_P_mc_truth_all_signal_AfterFSI->Write();
    muon_theta_mc_truth_all_signal_AfterFSI->Write();
    pi0_P_mc_truth_all_signal_AfterFSI->Write();
    pi0_KE_mc_truth_all_signal_AfterFSI->Write();
    pi0_theta_mc_truth_all_signal_AfterFSI->Write();
    Enu_mc_truth_all_signal_AfterFSI->Write();
    QSq_mc_truth_all_signal_AfterFSI->Write();
    W_mc_truth_all_signal_AfterFSI->Write();
    deltaInvMass_mc_truth_all_signal_AfterFSI->Write();
    Delta_pi_theta_mc_truth_all_signal_AfterFSI->Write();
    Delta_pi_phi_mc_truth_all_signal_AfterFSI->Write();

    // FSI Type
    for (int i = 0; i < nFSIType; ++i){
        muon_P_mc_truth_all_signal_FSIType[i]->Write();
        muon_theta_mc_truth_all_signal_FSIType[i]->Write();
        pi0_P_mc_truth_all_signal_FSIType[i]->Write();
        pi0_KE_mc_truth_all_signal_FSIType[i]->Write();
        pi0_theta_mc_truth_all_signal_FSIType[i]->Write();
        QSq_mc_truth_all_signal_FSIType[i]->Write();
        W_mc_truth_all_signal_FSIType[i]->Write();
        Enu_mc_truth_all_signal_FSIType[i]->Write();
        deltaInvMass_mc_truth_all_signal_FSIType[i]->Write();
        Delta_pi_theta_mc_truth_all_signal_FSIType[i]->Write();
        Delta_pi_phi_mc_truth_all_signal_FSIType[i]->Write();
    }

    // Interaction Type
    for (int i = 0; i < nIntType; ++i){
        muon_P_mc_truth_all_signal_IntType[i]->Write();
        muon_theta_mc_truth_all_signal_IntType[i]->Write();
        pi0_P_mc_truth_all_signal_IntType[i]->Write();
        pi0_KE_mc_truth_all_signal_IntType[i]->Write();
        pi0_theta_mc_truth_all_signal_IntType[i]->Write();
        QSq_mc_truth_all_signal_IntType[i]->Write();
        W_mc_truth_all_signal_IntType[i]->Write();
        Enu_mc_truth_all_signal_IntType[i]->Write();
        deltaInvMass_mc_truth_all_signal_IntType[i]->Write();
        Delta_pi_theta_mc_truth_all_signal_IntType[i]->Write();
        Delta_pi_phi_mc_truth_all_signal_IntType[i]->Write();
    }

    for (int i = 0; i < nHistograms; i++){
        CV_weight[i]->Write();
        CV_weight_2p2h[i]->Write();
        CV_weight_Delta[i]->Write();
        CV_weight_CCRES[i]->Write();
        CV_weight_NonRes1pi[i]->Write();
        h_err_2p2h[i]->Write();
        genie_wgt_Theta_Delta2Npi[i]->Write();
        updated_wgt_Theta_Delta2Npi[i]->Write();
        genie_wgt_MaRES[i]->Write();
        updated_wgt_MaRES[i]->Write();
        genie_wgt_MvRES[i]->Write();
        updated_wgt_MvRES[i]->Write();
        genie_wgt_Rvn1pi[i]->Write();
        updated_wgt_Rvn1pi[i]->Write();
    }

    // Signal Q2
    truth_QSq_QE->Write();
    truth_QSq_RES_1232->Write();
    truth_QSq_RES_1535->Write();
    truth_QSq_RES_1520->Write();
    truth_QSq_RES_Other->Write();
   
    truth_QSq_DIS->Write();
    truth_QSq_Non_RES->Write();
    truth_QSq_2p2h->Write();
 
    // Signal incomingE
    truth_Enu_QE->Write();
    truth_Enu_RES_1232->Write();
    truth_Enu_RES_1535->Write();
    truth_Enu_RES_1520->Write();
    truth_Enu_RES_Other->Write();
   
    truth_Enu_DIS->Write();
    truth_Enu_Non_RES->Write();
    truth_Enu_2p2h->Write();
 
    // Signal w
    truth_w_QE->Write();
    truth_w_RES_1232->Write();
    truth_w_RES_1535->Write();
    truth_w_RES_1520->Write();
    truth_w_RES_Other->Write();
   
    truth_w_DIS->Write();
    truth_w_Non_RES->Write();
    truth_w_2p2h->Write();

    delta_anisotropy->Write();

    f->Close();
}

#endif


