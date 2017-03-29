#ifndef CCProtonPi0_NTupleAnalysis_cpp
#define CCProtonPi0_NTupleanalysis_cpp

#include "CCProtonPi0_NTupleAnalysis.h"

using namespace PlotUtils;

// Initialize Constants
const std::string CCProtonPi0_NTupleAnalysis::version = "DeltaSuppressed";

const double CCProtonPi0_NTupleAnalysis::EPSILON = 1.0e-3; 

const double CCProtonPi0_NTupleAnalysis::data_POT = 3.33153e+20;
const double CCProtonPi0_NTupleAnalysis::mc_POT = 2.21867e+21; 
const double CCProtonPi0_NTupleAnalysis::mc_2p2h_POT = 2.56925e+21;
const double CCProtonPi0_NTupleAnalysis::POT_ratio = data_POT/mc_POT;
const double CCProtonPi0_NTupleAnalysis::POT_ratio_2p2h = mc_POT/mc_2p2h_POT;

const double CCProtonPi0_NTupleAnalysis::max_muon_theta = 25; // degree
const double CCProtonPi0_NTupleAnalysis::min_Enu = 1500; // MeV
//const double CCProtonPi0_NTupleAnalysis::max_Enu = 4000; // MeV -- QSq Enu Study
const double CCProtonPi0_NTupleAnalysis::max_Enu = 20000; // MeV
//const double CCProtonPi0_NTupleAnalysis::min_Enu = 4000; // MeV -- QSq Enu Study
const double CCProtonPi0_NTupleAnalysis::max_W = 1800; // MeV


const double CCProtonPi0_NTupleAnalysis::SENTINEL = -9.9;
const double CCProtonPi0_NTupleAnalysis::MeV_to_GeV = pow(10,-3);
const double CCProtonPi0_NTupleAnalysis::MeVSq_to_GeVSq = pow(10,-6);
const double CCProtonPi0_NTupleAnalysis::mSq_to_cmSq = pow(10,4);
const double CCProtonPi0_NTupleAnalysis::mm_to_cm = pow(10,-1);
const double CCProtonPi0_NTupleAnalysis::rad_to_deg = 180.0/M_PI;
const double CCProtonPi0_NTupleAnalysis::beam_theta = -0.05887;  //rad 
const double CCProtonPi0_NTupleAnalysis::beam_phi = M_PI_2;  //rad 

const double CCProtonPi0_NTupleAnalysis::muon_mass = 105.6583;      // MeV
const double CCProtonPi0_NTupleAnalysis::pi0_mass = 134.9766;       // MeV
const double CCProtonPi0_NTupleAnalysis::piplus_mass = 139.5701;    // MeV
const double CCProtonPi0_NTupleAnalysis::proton_mass = 938.272013;  // MeV
const double CCProtonPi0_NTupleAnalysis::neutron_mass = 939.56536;  // MeV

// GENIE Tuning
const double CCProtonPi0_NTupleAnalysis::genieMaRes              = 1.12;
const double CCProtonPi0_NTupleAnalysis::genieMaRes1sig          = 0.2 * genieMaRes;
// GENIE central value MvRES from electroproduction data fit
const double CCProtonPi0_NTupleAnalysis::genieMvRes              = 0.84;
const double CCProtonPi0_NTupleAnalysis::genieMvRes1sig          = 0.1 * genieMvRes;
// Reduced MvRES error from electroproduction data fit
const double CCProtonPi0_NTupleAnalysis::electroProdMvRes1sig    = 0.03 * genieMvRes;
// Pion production parameters and errors from deuterium fit
const double CCProtonPi0_NTupleAnalysis::deuteriumMaRes          = 0.94;
const double CCProtonPi0_NTupleAnalysis::deuteriumMaRes1sig      = 0.05;
const double CCProtonPi0_NTupleAnalysis::deuteriumNonResNorm     = 0.46;
const double CCProtonPi0_NTupleAnalysis::deuteriumNonResNorm1sig = 0.04;
const double CCProtonPi0_NTupleAnalysis::deuteriumResNorm        = 1.15;
const double CCProtonPi0_NTupleAnalysis::deuteriumResNorm1sig    = 0.07;
// Delta Suppression
const double CCProtonPi0_NTupleAnalysis::DeltaFactor_A = 1.0;
const double CCProtonPi0_NTupleAnalysis::DeltaFactor_Q0 = 0.116;

// Flux Correction
const bool CCProtonPi0_NTupleAnalysis::applyNuEConstraint = true;
const FluxReweighter::EFluxVersion CCProtonPi0_NTupleAnalysis::new_flux = FluxReweighter::gen2thin;
const FluxReweighter::EG4NumiVersion CCProtonPi0_NTupleAnalysis::old_flux = FluxReweighter::g4numiv5;

CCProtonPi0_NTupleAnalysis::CCProtonPi0_NTupleAnalysis()
{
    // Required for MINERvA Framework Classes
    ROOT::Cintex::Cintex::Enable();

    // For Adding Error Bars, Initialize FluxReweighter with minerva1
    // It will be updated according to playlist of the run during Analysis Stage
    std::cout<<"Initializing FluxReweighter for Flux Error Band"<<std::endl;
    frw_DefaultInit = true;    
    frw = new FluxReweighter(14, applyNuEConstraint, FluxReweighter::minerva1, new_flux, old_flux);
    processed_minerva1 = false;
    processed_minerva7 = false;
    processed_minerva9 = false;
    processed_minerva13A = false;
    processed_minerva13B = false;
    processed_minerva13C = false;
    processed_minerva13D = false;
    processed_minerva13E = false;
    processed_2p2h = false;

    init2p2hFitResults();
}

void CCProtonPi0_NTupleAnalysis::OpenTextFile(std::string file_name, std::ofstream &file)
{
    file.open(file_name.c_str());
    if (!file.is_open()){
        std::cerr<<"Cannot open output text file: "<<file_name<<std::endl;
        exit(1);
    }else{
        std::cout<<"\t"<<file_name<<std::endl;
    }
}

MnvH1D* CCProtonPi0_NTupleAnalysis::GetMnvH1D(TFile* f, std::string var_name)
{
    MnvH1D* h = new MnvH1D( * dynamic_cast<MnvH1D*>(f->Get(var_name.c_str())) );
    h->SetDirectory(NULL);
    return h;
}

MnvH2D* CCProtonPi0_NTupleAnalysis::GetMnvH2D(TFile* f, std::string var_name)
{
    MnvH2D* h = new MnvH2D( * dynamic_cast<MnvH2D*>(f->Get(var_name.c_str())) );
    h->SetDirectory(NULL);
    return h;
}

// --------------------------------------------------------------------
// Errors for Data
// --------------------------------------------------------------------
    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBands_Data(MnvHistoType* h)
{
    AddVertErrorBandAndFillWithCV_Genie(h);
    AddVertErrorBandAndFillWithCV_Flux(h);
    AddVertErrorBandAndFillWithCV_BckgConstraint_WithPi0(h);
    AddVertErrorBandAndFillWithCV_BckgConstraint_SinglePiPlus(h);
    AddVertErrorBandAndFillWithCV_BckgConstraint_QELike(h);
    AddVertErrorBandAndFillWithCV_MichelTrue(h);
    AddVertErrorBandAndFillWithCV_MichelFake(h);
    AddVertErrorBandAndFillWithCV_TargetMass(h);
    AddVertErrorBandAndFillWithCV_2p2h(h);
    AddVertErrorBandAndFillWithCV_Unfolding(h);
    AddVertErrorBandAndFillWithCV_MuonTracking(h);
    AddVertErrorBandAndFillWithCV_ProtonTracking(h);
    AddVertErrorBandAndFillWithCV_NeutronResponse(h);
    AddVertErrorBandAndFillWithCV_PionResponse(h);

    // QSq Study
    //AddVertErrorBandAndFillWithCV_HighMaRES(h);
    //AddVertErrorBandAndFillWithCV_LowMaRES(h);
    //AddVertErrorBandAndFillWithCV_DeltaFactor(h);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBands_Data<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBands_Data<MnvH2D>(MnvH2D* h);

// Truth Tree Only have GENIE and Flux Errors, others are handled as Data
    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBands_TruthTree(MnvHistoType* h)
{
    AddVertErrorBand_Flux(h);
    AddVertErrorBand_Genie(h);
    AddVertErrorBand_2p2h(h);

    //AddVertErrorBand_HighMaRES(h);
    //AddVertErrorBand_LowMaRES(h);
    //AddVertErrorBand_DeltaFactor(h);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBands_TruthTree<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBands_TruthTree<MnvH2D>(MnvH2D* h);
    
    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBandsAndFillWithCV_TruthTree(MnvHistoType* h)
{
    AddVertErrorBandAndFillWithCV_BckgConstraint_WithPi0(h);
    AddVertErrorBandAndFillWithCV_BckgConstraint_SinglePiPlus(h);
    AddVertErrorBandAndFillWithCV_BckgConstraint_QELike(h);
    AddVertErrorBandAndFillWithCV_MichelTrue(h);
    AddVertErrorBandAndFillWithCV_MichelFake(h);
    AddVertErrorBandAndFillWithCV_TargetMass(h);
    AddVertErrorBandAndFillWithCV_Unfolding(h);
    AddVertErrorBandAndFillWithCV_MuonTracking(h);
    AddVertErrorBandAndFillWithCV_ProtonTracking(h);
    AddVertErrorBandAndFillWithCV_NeutronResponse(h);
    AddVertErrorBandAndFillWithCV_PionResponse(h);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandsAndFillWithCV_TruthTree<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandsAndFillWithCV_TruthTree<MnvH2D>(MnvH2D* h);

// Flux Histogram have Flux Errors, others are handled as Data
    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBands_FluxHistogram(MnvHistoType* h)
{
    AddVertErrorBandAndFillWithCV_Genie(h);
    AddVertErrorBandAndFillWithCV_BckgConstraint_WithPi0(h);
    AddVertErrorBandAndFillWithCV_BckgConstraint_SinglePiPlus(h);
    AddVertErrorBandAndFillWithCV_BckgConstraint_QELike(h);
    AddVertErrorBandAndFillWithCV_MichelTrue(h);
    AddVertErrorBandAndFillWithCV_MichelFake(h);
    AddVertErrorBandAndFillWithCV_TargetMass(h);
    AddVertErrorBandAndFillWithCV_2p2h(h);
    AddVertErrorBandAndFillWithCV_Unfolding(h);
    AddVertErrorBandAndFillWithCV_MuonTracking(h);
    AddVertErrorBandAndFillWithCV_ProtonTracking(h);
    AddVertErrorBandAndFillWithCV_NeutronResponse(h);
    AddVertErrorBandAndFillWithCV_PionResponse(h);

    // QSq Study
    //AddVertErrorBandAndFillWithCV_HighMaRES(h);
    //AddVertErrorBandAndFillWithCV_LowMaRES(h);
    //AddVertErrorBandAndFillWithCV_DeltaFactor(h);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBands_FluxHistogram<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBands_FluxHistogram<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_Flux(MnvHistoType* h)
{
    h->AddVertErrorBandAndFillWithCV("Flux",  100);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_Flux<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_Flux<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_Genie(MnvHistoType* h)
{
    h->AddVertErrorBandAndFillWithCV("GENIE_AGKYxF1pi"         ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_AhtBY"             ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_BhtBY"             ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_CCQEPauliSupViaKF" ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_CV1uBY"            ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_CV2uBY"            ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_EtaNCEL"           ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_FrAbs_N"           ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_FrAbs_pi"          ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_FrCEx_N"           ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_FrCEx_pi"          ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_FrElas_N"          ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_FrElas_pi"         ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_FrInel_N"          ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_FrInel_pi"         ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_FrPiProd_N"        ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_FrPiProd_pi"       ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_MFP_N"             ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_MFP_pi"            ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_MaCCQE"            ,2);
    //h->AddVertErrorBandAndFillWithCV("GENIE_MaCCQEshape"       ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_MaNCEL"            ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_MaRES"             ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_MvRES"             ,2);
    //h->AddVertErrorBandAndFillWithCV("GENIE_NormCCQE"          ,2);
    //h->AddVertErrorBandAndFillWithCV("GENIE_NormCCRES"         ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_NormDISCC"         ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_NormNCRES"         ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_RDecBR1gamma"      ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_Rvn1pi"            ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_Rvn2pi"            ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_Rvp1pi"            ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_Rvp2pi"            ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_Theta_Delta2Npi"   ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_VecFFCCQEshape"    ,2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_Genie<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_Genie<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_BckgConstraint_WithPi0(MnvHistoType* h)
{
    h->AddVertErrorBandAndFillWithCV("BckgConstraint_WithPi0", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_BckgConstraint_WithPi0<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_BckgConstraint_WithPi0<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_BckgConstraint_SinglePiPlus(MnvHistoType* h)
{
    h->AddVertErrorBandAndFillWithCV("BckgConstraint_SinglePiPlus", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_BckgConstraint_SinglePiPlus<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_BckgConstraint_SinglePiPlus<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_BckgConstraint_QELike(MnvHistoType* h)
{
    h->AddVertErrorBandAndFillWithCV("BckgConstraint_QELike", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_BckgConstraint_QELike<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_BckgConstraint_QELike<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_MichelFake(MnvHistoType* h)
{
    h->AddVertErrorBandAndFillWithCV("MichelFake", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_MichelFake<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_MichelFake<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_MichelTrue(MnvHistoType* h)
{
    h->AddVertErrorBandAndFillWithCV("MichelTrue", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_MichelTrue<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_MichelTrue<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_TargetMass(MnvHistoType* h)
{
    h->AddVertErrorBandAndFillWithCV("TargetMass", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_TargetMass<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_TargetMass<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_2p2h(MnvHistoType* h)
{
    h->AddVertErrorBandAndFillWithCV("2p2h", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_2p2h<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_2p2h<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_Unfolding(MnvHistoType* h)
{
    h->AddVertErrorBandAndFillWithCV("Unfolding", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_Unfolding<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_Unfolding<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_ProtonTracking(MnvHistoType* h)
{
    h->AddVertErrorBandAndFillWithCV("ProtonTracking", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_ProtonTracking<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_ProtonTracking<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_MuonTracking(MnvHistoType* h)
{
    h->AddVertErrorBandAndFillWithCV("MuonTracking", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_MuonTracking<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_MuonTracking<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_NeutronResponse(MnvHistoType* h)
{
    h->AddVertErrorBandAndFillWithCV("NeutronResponse", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_NeutronResponse<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_NeutronResponse<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_PionResponse(MnvHistoType* h)
{
    h->AddVertErrorBandAndFillWithCV("PionResponse", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_PionResponse<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_PionResponse<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_HighMaRES(MnvHistoType* h)
{
    h->AddVertErrorBandAndFillWithCV("HighMaRES", 201);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_HighMaRES<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_HighMaRES<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_LowMaRES(MnvHistoType* h)
{
    h->AddVertErrorBandAndFillWithCV("LowMaRES", 201);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_LowMaRES<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_LowMaRES<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_DeltaFactor(MnvHistoType* h)
{
    h->AddVertErrorBandAndFillWithCV("DeltaFactor", 106);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_DeltaFactor<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_DeltaFactor<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddLatErrorBands_Data(MnvHistoType* h)
{
    AddLatErrorBandAndFillWithCV_ProtonEnergy_MassModel(h);
    AddLatErrorBandAndFillWithCV_ProtonEnergy_MEU(h);
    AddLatErrorBandAndFillWithCV_ProtonEnergy_BetheBloch(h);
    AddLatErrorBandAndFillWithCV_ProtonEnergy_Birks(h);
    AddLatErrorBandAndFillWithCV_MuonMomentum(h);
    AddLatErrorBandAndFillWithCV_MuonTheta(h);
    AddLatErrorBandAndFillWithCV_EM_EnergyScale(h);
}
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBands_Data<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBands_Data<MnvH2D>(MnvH2D* h);

// All Truth Tree Lateral Error Bands are Handled as Data
    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddLatErrorBandsAndFillWithCV_TruthTree(MnvHistoType* h)
{
    AddLatErrorBands_Data(h);
}
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBandsAndFillWithCV_TruthTree<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBandsAndFillWithCV_TruthTree<MnvH2D>(MnvH2D* h);

// All Flux Histogram Lateral Error Bands are Handled as Data
    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddLatErrorBands_FluxHistogram(MnvHistoType* h)
{
    AddLatErrorBands_Data(h);
}
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBands_FluxHistogram<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBands_FluxHistogram<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddLatErrorBandAndFillWithCV_ProtonEnergy_MassModel(MnvHistoType* h)
{
    h->AddLatErrorBandAndFillWithCV("ProtonEnergy_MassModel", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBandAndFillWithCV_ProtonEnergy_MassModel<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBandAndFillWithCV_ProtonEnergy_MassModel<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddLatErrorBandAndFillWithCV_ProtonEnergy_MEU(MnvHistoType* h)
{
    h->AddLatErrorBandAndFillWithCV("ProtonEnergy_MEU", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBandAndFillWithCV_ProtonEnergy_MEU<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBandAndFillWithCV_ProtonEnergy_MEU<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddLatErrorBandAndFillWithCV_ProtonEnergy_BetheBloch(MnvHistoType* h)
{
    h->AddLatErrorBandAndFillWithCV("ProtonEnergy_BetheBloch", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBandAndFillWithCV_ProtonEnergy_BetheBloch<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBandAndFillWithCV_ProtonEnergy_BetheBloch<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddLatErrorBandAndFillWithCV_ProtonEnergy_Birks(MnvHistoType* h)
{
    h->AddLatErrorBandAndFillWithCV("ProtonEnergy_Birks", n_lateral_universes);
}
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBandAndFillWithCV_ProtonEnergy_Birks<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBandAndFillWithCV_ProtonEnergy_Birks<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddLatErrorBandAndFillWithCV_MuonMomentum(MnvHistoType* h)
{
    h->AddLatErrorBandAndFillWithCV("MuonMomentum", n_lateral_universes);
}
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBandAndFillWithCV_MuonMomentum<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBandAndFillWithCV_MuonMomentum<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddLatErrorBandAndFillWithCV_MuonTheta(MnvHistoType* h)
{
    h->AddLatErrorBandAndFillWithCV("MuonTheta", n_lateral_universes);
}
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBandAndFillWithCV_MuonTheta<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBandAndFillWithCV_MuonTheta<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddLatErrorBandAndFillWithCV_EM_EnergyScale(MnvHistoType* h)
{
    h->AddLatErrorBandAndFillWithCV("EM_EnergyScale", n_lateral_universes);
}
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBandAndFillWithCV_EM_EnergyScale<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBandAndFillWithCV_EM_EnergyScale<MnvH2D>(MnvH2D* h);


// --------------------------------------------------------------------
//  Errors for MC
// --------------------------------------------------------------------
    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBands_MC(MnvHistoType* h)
{
    AddVertErrorBand_Genie(h);
    AddVertErrorBand_Flux(h);
    AddVertErrorBand_BckgConstraint_WithPi0(h);
    AddVertErrorBand_BckgConstraint_SinglePiPlus(h);
    AddVertErrorBand_BckgConstraint_QELike(h);
    AddVertErrorBand_MichelTrue(h);
    AddVertErrorBand_MichelFake(h);
    AddVertErrorBand_TargetMass(h);
    AddVertErrorBand_2p2h(h);
    AddVertErrorBand_Unfolding(h);
    AddVertErrorBand_MuonTracking(h);
    AddVertErrorBand_ProtonTracking(h);
    AddVertErrorBand_NeutronResponse(h);
    AddVertErrorBand_PionResponse(h);

    // QSq Study
    //AddVertErrorBand_HighMaRES(h);
    //AddVertErrorBand_LowMaRES(h);
    //AddVertErrorBand_DeltaFactor(h);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBands_MC<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBands_MC<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_Flux(MnvHistoType* h)
{
    frw->AddFluxErrorBand(h);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_Flux<MnvH1D>( MnvH1D* h );
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_Flux<MnvH2D>( MnvH2D* h );

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_Genie(MnvHistoType* h)
{
    h->AddVertErrorBand("GENIE_AGKYxF1pi"         ,2);
    h->AddVertErrorBand("GENIE_AhtBY"             ,2);
    h->AddVertErrorBand("GENIE_BhtBY"             ,2);
    h->AddVertErrorBand("GENIE_CCQEPauliSupViaKF" ,2);
    h->AddVertErrorBand("GENIE_CV1uBY"            ,2);
    h->AddVertErrorBand("GENIE_CV2uBY"            ,2);
    h->AddVertErrorBand("GENIE_EtaNCEL"           ,2);
    h->AddVertErrorBand("GENIE_FrAbs_N"           ,2);
    h->AddVertErrorBand("GENIE_FrAbs_pi"          ,2);
    h->AddVertErrorBand("GENIE_FrCEx_N"           ,2);
    h->AddVertErrorBand("GENIE_FrCEx_pi"          ,2);
    h->AddVertErrorBand("GENIE_FrElas_N"          ,2);
    h->AddVertErrorBand("GENIE_FrElas_pi"         ,2);
    h->AddVertErrorBand("GENIE_FrInel_N"          ,2);
    h->AddVertErrorBand("GENIE_FrInel_pi"         ,2);
    h->AddVertErrorBand("GENIE_FrPiProd_N"        ,2);
    h->AddVertErrorBand("GENIE_FrPiProd_pi"       ,2);
    h->AddVertErrorBand("GENIE_MFP_N"             ,2);
    h->AddVertErrorBand("GENIE_MFP_pi"            ,2);
    h->AddVertErrorBand("GENIE_MaCCQE"            ,2);
    //h->AddVertErrorBand("GENIE_MaCCQEshape"       ,2);
    h->AddVertErrorBand("GENIE_MaNCEL"            ,2);
    h->AddVertErrorBand("GENIE_MaRES"             ,2);
    h->AddVertErrorBand("GENIE_MvRES"             ,2);
    //h->AddVertErrorBand("GENIE_NormCCQE"          ,2);
    //h->AddVertErrorBand("GENIE_NormCCRES"         ,2);
    h->AddVertErrorBand("GENIE_NormDISCC"         ,2);
    h->AddVertErrorBand("GENIE_NormNCRES"         ,2);
    h->AddVertErrorBand("GENIE_RDecBR1gamma"      ,2);
    h->AddVertErrorBand("GENIE_Rvn1pi"            ,2);
    h->AddVertErrorBand("GENIE_Rvn2pi"            ,2);
    h->AddVertErrorBand("GENIE_Rvp1pi"            ,2);
    h->AddVertErrorBand("GENIE_Rvp2pi"            ,2);
    h->AddVertErrorBand("GENIE_Theta_Delta2Npi"   ,2);
    h->AddVertErrorBand("GENIE_VecFFCCQEshape"    ,2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_Genie<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_Genie<MnvH2D>(MnvH2D* h);

     template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_BckgConstraint_WithPi0(MnvHistoType* h)
{
    h->AddVertErrorBand("BckgConstraint_WithPi0", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_BckgConstraint_WithPi0<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_BckgConstraint_WithPi0<MnvH2D>(MnvH2D* h);

     template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_BckgConstraint_SinglePiPlus(MnvHistoType* h)
{
    h->AddVertErrorBand("BckgConstraint_SinglePiPlus", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_BckgConstraint_SinglePiPlus<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_BckgConstraint_SinglePiPlus<MnvH2D>(MnvH2D* h);

     template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_BckgConstraint_QELike(MnvHistoType* h)
{
    h->AddVertErrorBand("BckgConstraint_QELike", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_BckgConstraint_QELike<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_BckgConstraint_QELike<MnvH2D>(MnvH2D* h);

     template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_MichelFake(MnvHistoType* h)
{
    h->AddVertErrorBand("MichelFake", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_MichelFake<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_MichelFake<MnvH2D>(MnvH2D* h);

     template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_MichelTrue(MnvHistoType* h)
{
    h->AddVertErrorBand("MichelTrue", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_MichelTrue<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_MichelTrue<MnvH2D>(MnvH2D* h);

     template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_TargetMass(MnvHistoType* h)
{
    h->AddVertErrorBand("TargetMass", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_TargetMass<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_TargetMass<MnvH2D>(MnvH2D* h);

     template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_2p2h(MnvHistoType* h)
{
    h->AddVertErrorBand("2p2h", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_2p2h<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_2p2h<MnvH2D>(MnvH2D* h);

     template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_Unfolding(MnvHistoType* h)
{
    h->AddVertErrorBand("Unfolding", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_Unfolding<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_Unfolding<MnvH2D>(MnvH2D* h);

   template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_ProtonTracking(MnvHistoType* h)
{
    h->AddVertErrorBand("ProtonTracking", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_ProtonTracking<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_ProtonTracking<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_MuonTracking(MnvHistoType* h)
{
    h->AddVertErrorBand("MuonTracking", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_MuonTracking<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_MuonTracking<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_NeutronResponse(MnvHistoType* h)
{
    h->AddVertErrorBand("NeutronResponse", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_NeutronResponse<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_NeutronResponse<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_PionResponse(MnvHistoType* h)
{
    h->AddVertErrorBand("PionResponse", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_PionResponse<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_PionResponse<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_HighMaRES(MnvHistoType* h)
{
    h->AddVertErrorBand("HighMaRES", 201);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_HighMaRES<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_HighMaRES<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_LowMaRES(MnvHistoType* h)
{
    h->AddVertErrorBand("LowMaRES", 201);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_LowMaRES<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_LowMaRES<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_DeltaFactor(MnvHistoType* h)
{
    h->AddVertErrorBand("DeltaFactor", 106);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_DeltaFactor<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_DeltaFactor<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddLatErrorBands_MC(MnvHistoType* h)
{
    AddLatErrorBand_ProtonEnergy_MassModel(h);
    AddLatErrorBand_ProtonEnergy_MEU(h);
    AddLatErrorBand_ProtonEnergy_BetheBloch(h);
    AddLatErrorBand_ProtonEnergy_Birks(h);
    AddLatErrorBand_MuonMomentum(h);
    AddLatErrorBand_MuonTheta(h);
    AddLatErrorBand_EM_EnergyScale(h);
}
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBands_MC<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBands_MC<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddLatErrorBand_ProtonEnergy_MassModel(MnvHistoType* h)
{
    h->AddLatErrorBand("ProtonEnergy_MassModel", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBand_ProtonEnergy_MassModel<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBand_ProtonEnergy_MassModel<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddLatErrorBand_ProtonEnergy_MEU(MnvHistoType* h)
{
    h->AddLatErrorBand("ProtonEnergy_MEU", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBand_ProtonEnergy_MEU<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBand_ProtonEnergy_MEU<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddLatErrorBand_ProtonEnergy_BetheBloch(MnvHistoType* h)
{
    h->AddLatErrorBand("ProtonEnergy_BetheBloch", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBand_ProtonEnergy_BetheBloch<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBand_ProtonEnergy_BetheBloch<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddLatErrorBand_MuonMomentum(MnvHistoType* h)
{
    h->AddLatErrorBand("MuonMomentum", n_lateral_universes);
}
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBand_MuonMomentum<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBand_MuonMomentum<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddLatErrorBand_MuonTheta(MnvHistoType* h)
{
    h->AddLatErrorBand("MuonTheta", n_lateral_universes);
}
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBand_MuonTheta<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBand_MuonTheta<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddLatErrorBand_ProtonEnergy_Birks(MnvHistoType* h)
{
    h->AddLatErrorBand("ProtonEnergy_Birks", n_lateral_universes);
}
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBand_ProtonEnergy_Birks<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBand_ProtonEnergy_Birks<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddLatErrorBand_EM_EnergyScale(MnvHistoType* h)
{
    h->AddLatErrorBand("EM_EnergyScale", n_lateral_universes);
}
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBand_EM_EnergyScale<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBand_EM_EnergyScale<MnvH2D>(MnvH2D* h);


std::string CCProtonPi0_NTupleAnalysis::GetPlaylist(const int run, int type)
{
    std::string playlist;

    if ( IsEvent2p2h(type) ) playlist = "minerva_2p2h";
    else if (run >= 10200 && run <= 10249) playlist = "minerva1";
    else if (run >= 10250 && run <= 10254) playlist = "minerva7";
    else if (run >= 10255 && run <= 10259) playlist = "minerva9";
    else if (run >= 12200 && run <= 12209) playlist = "minerva13A";
    else if (run >= 12210 && run <= 12219) playlist = "minerva13B";
    else if (run >= 13200 && run <= 13299) playlist = "minerva13C";
    else if (run >= 14201 && run <= 14209) playlist = "minerva13D";
    else if (run >= 14210 && run <= 14229) playlist = "minerva13E";
    else RunTimeError("ERROR: NO Playlist Found!");

    return playlist;
}

void CCProtonPi0_NTupleAnalysis::UpdateFluxReweighter(int run, int type)
{
    std::string playlist = GetPlaylist(run, type);

    if (!processed_minerva1 && playlist.compare("minerva1") == 0){
        std::cout<<"Playlist: minerva1"<<std::endl;
        ReInitFluxReweighter(FluxReweighter::minerva1);
        processed_minerva1 = true;
    }else if (!processed_minerva7 && playlist.compare("minerva7") == 0){
        std::cout<<"Playlist: minerva7"<<std::endl;
        ReInitFluxReweighter(FluxReweighter::minervaLE_FHC);
        processed_minerva7 = true;
    }else if (!processed_minerva9 && playlist.compare("minerva9") == 0){
        std::cout<<"Playlist: minerva9"<<std::endl;
        ReInitFluxReweighter(FluxReweighter::minervaLE_FHC);
        processed_minerva9 = true;
    }else if (!processed_minerva13A && playlist.compare("minerva13A") == 0){
        std::cout<<"Playlist: minerva13A"<<std::endl;
        ReInitFluxReweighter(FluxReweighter::minerva13);
        processed_minerva13A = true;
    }else if (!processed_minerva13B && playlist.compare("minerva13B") == 0){
        std::cout<<"Playlist: minerva13B"<<std::endl;
        ReInitFluxReweighter(FluxReweighter::minerva13);
        processed_minerva13B = true;
    }else if (!processed_minerva13C && playlist.compare("minerva13C") == 0){
        std::cout<<"Playlist: minerva13C"<<std::endl;
        ReInitFluxReweighter(FluxReweighter::minerva13);
        processed_minerva13C = true;
    }else if (!processed_minerva13D && playlist.compare("minerva13D") == 0){
        std::cout<<"Playlist: minerva13D"<<std::endl;
        ReInitFluxReweighter(FluxReweighter::minerva13);
        processed_minerva13D = true;
    }else if (!processed_minerva13E && playlist.compare("minerva13E") == 0){
        std::cout<<"Playlist: minerva13E"<<std::endl;
        ReInitFluxReweighter(FluxReweighter::minerva13);
        processed_minerva13E = true;
    }else if ( processed_minerva1 && processed_minerva13C && processed_minerva7 && processed_minerva9 ){
        if (!processed_2p2h && playlist.compare("minerva_2p2h") == 0){
            std::cout<<"Playlist: 2p2h"<<std::endl;
            ReInitFluxReweighter(FluxReweighter::minerva13);
            processed_2p2h = true;
        }
    }
}

void CCProtonPi0_NTupleAnalysis::ReInitFluxReweighter(enum FluxReweighter::EPlaylist playlist)
{
    if(frw != NULL || frw_DefaultInit){
        frw_DefaultInit = false; 
        delete frw;
    }
    frw = new FluxReweighter(14, applyNuEConstraint, playlist, new_flux, old_flux);
}

void CCProtonPi0_NTupleAnalysis::printBins(const TH2* hist, const std::string var_name)
{
    std::cout<<"Printing Bin Content for "<<var_name<<std::endl;

    int nBinsX = hist->GetNbinsX();
    int nBinsY = hist->GetNbinsY();

    for (int y = nBinsY; y >= 1; --y){
        for (int x = 1; x <= nBinsX; ++x){
            double content = hist->GetBinContent(x,y);
            std::cout<<content<<" ";
        }
        std::cout<<std::endl;
    }
}

void CCProtonPi0_NTupleAnalysis::printBins(const TH1* hist, const std::string var_name, bool useLowEdge)
{
    std::cout<<std::left;
    std::cout<<"Printing Bin Content for "<<var_name<<std::endl;
    if (useLowEdge){
        std::cout.width(12); std::cout<<"BinLowEdge";
    }else{
        std::cout.width(12); std::cout<<"BinCenter";
    }
    std::cout.width(12); std::cout<<"Content";
    std::cout.width(12); std::cout<<"Percent"<<std::endl;

    double nEntries = hist->Integral();
    double total = 0;
    int nBins = hist->GetNbinsX();
    for (int i = 1; i <= nBins; ++i){
        double content = hist->GetBinContent(i);
        total += content;
        if (useLowEdge){
            std::cout.width(12); std::cout<<hist->GetBinLowEdge(i);
        }else{
            std::cout.width(12); std::cout<<hist->GetBinCenter(i);
        }
        std::cout.width(12); std::cout<<std::setprecision(3)<<content;
        std::cout.width(12); std::cout<<std::setprecision(3)<<content/nEntries*100;
        std::cout<<std::endl;
    }

    std::cout.width(12); std::cout<<"Total";
    std::cout.width(12); std::cout<<std::setprecision(3)<<total;
    std::cout.width(12); std::cout<<std::setprecision(3)<<total/nEntries*100;
    std::cout<<std::endl;

    std::cout<<"UnderFlow = "<<hist->GetBinContent(0)<<std::endl;
    std::cout<<"OverFlow = "<<hist->GetBinContent(nBins+1)<<std::endl;

}

bool CCProtonPi0_NTupleAnalysis::IsWInRange(double W)
{
    if (W <= max_W) return true;
    else return false;
}

bool CCProtonPi0_NTupleAnalysis::IsEnuInRange(double Enu)
{
    if (Enu >= min_Enu && Enu <= max_Enu) return true;
    //if (Enu >= min_Enu && Enu < max_Enu) return true; // QSq Study
    else return false;
}

double CCProtonPi0_NTupleAnalysis::Calc_Enu_Truth(double muon_E, double proton_E, double pi0_E)
{
    double Enu;
    double proton_KE = proton_E - proton_mass;

    Enu = muon_E + proton_KE + pi0_E; 

    return Enu;
}

double CCProtonPi0_NTupleAnalysis::Calc_QSq(double Enu, double muon_E, double muon_P, double muon_angle_beam) 
{
    double QSq = 2*Enu*(muon_E - muon_P*cos(muon_angle_beam))-(muon_mass*muon_mass);

    return QSq;
}

double CCProtonPi0_NTupleAnalysis::Calc_WSq(double Enu, double QSq, double muon_E)
{
    // Calculate WSq - Use eq. in Research Logbook page 31
    double WSq = neutron_mass*neutron_mass + 2*neutron_mass*(Enu - muon_E) - QSq; 

    return WSq;
}

void CCProtonPi0_NTupleAnalysis::GetAllUniverses(MnvH1D* mnvh1d_hist, std::vector<TH1D*> &all_universes, std::vector<std::string> &err_bands, std::vector<int> &hist_ind)
{
    // Check for input vector
    if (!all_universes.empty()){
        std::cout<<"WARNING! input vector<TH1D*> all_universes is NOT empty"<<std::endl;
        std::cout<<"Returning without change!"<<std::endl;
        return;
    }

    // ------------------------------------------------------------------------
    // Add CV Histogram as first element -- all_universes[0]
    // ------------------------------------------------------------------------
    TH1D* cv_hist =  new TH1D(mnvh1d_hist->GetCVHistoWithStatError());
    all_universes.push_back(cv_hist);
    err_bands.push_back("CentralValue");
    hist_ind.push_back(0);

    // ------------------------------------------------------------------------
    // Add Other Universes from Error Bands
    // ------------------------------------------------------------------------
    // Get Vert Error Band Names
    std::vector<std::string> vert_err_names = mnvh1d_hist->GetVertErrorBandNames();

    // Loop over all Vertical Error Bands
    for (unsigned int i = 0; i < vert_err_names.size(); ++i){
        MnvVertErrorBand* err_band =  mnvh1d_hist->GetVertErrorBand(vert_err_names[i]);
        // Get All Histograms from it
        std::vector<TH1D*> err_hists = err_band->GetHists();
        for (unsigned int j = 0; j < err_hists.size(); ++j){
            TH1D* temp = new TH1D(*err_hists[j]);
            all_universes.push_back(temp);
            err_bands.push_back(vert_err_names[i]);
            hist_ind.push_back(j);
        }
    }

    // Get Lat Error Band Names
    std::vector<std::string> lat_err_names = mnvh1d_hist->GetLatErrorBandNames();

    // Loop over all Lateral Error Bands
    for (unsigned int i = 0; i < lat_err_names.size(); ++i){
        MnvLatErrorBand* err_band =  mnvh1d_hist->GetLatErrorBand(lat_err_names[i]);
        // Get All Histograms from it
        std::vector<TH1D*> err_hists = err_band->GetHists();
        for (unsigned int j = 0; j < err_hists.size(); ++j){
            TH1D* temp = new TH1D(*err_hists[j]);
            all_universes.push_back(temp);
            err_bands.push_back(lat_err_names[i]);
            hist_ind.push_back(j);
        }
    }
}

void CCProtonPi0_NTupleAnalysis::GetPointersAllUniverses(MnvH1D* mnvh1d_hist, std::vector<TH1D*> &all_universes)
{
    // ------------------------------------------------------------------------
    // Add All Universes from All Error Bands
    // ------------------------------------------------------------------------
    // Get Vert Error Band Names
    std::vector<std::string> vert_err_names = mnvh1d_hist->GetVertErrorBandNames();

    // Loop over all Vertical Error Bands
    for (unsigned int i = 0; i < vert_err_names.size(); ++i){
        MnvVertErrorBand* err_band =  mnvh1d_hist->GetVertErrorBand(vert_err_names[i]);
        // Get All Histograms from it
        std::vector<TH1D*> err_hists = err_band->GetHists();
        for (unsigned int j = 0; j < err_hists.size(); ++j){
            all_universes.push_back(err_hists[j]);
        }
    }

    // Get Lat Error Band Names
    std::vector<std::string> lat_err_names = mnvh1d_hist->GetLatErrorBandNames();

    // Loop over all Lateral Error Bands
    for (unsigned int i = 0; i < lat_err_names.size(); ++i){
        MnvLatErrorBand* err_band =  mnvh1d_hist->GetLatErrorBand(lat_err_names[i]);
        // Get All Histograms from it
        std::vector<TH1D*> err_hists = err_band->GetHists();
        for (unsigned int j = 0; j < err_hists.size(); ++j){
            all_universes.push_back(err_hists[j]);
        }
    }
}

void CCProtonPi0_NTupleAnalysis::ClearAllUniversesVector(std::vector<TH1D*> &all_universes)
{
    // Delete all TH1Ds First
    for (unsigned int i = 0; i < all_universes.size(); ++i){
        delete all_universes[i];
    }

    // Remove All Elements from Vector
    all_universes.clear();

}

void CCProtonPi0_NTupleAnalysis::RunTimeError(std::string message)
{
    std::cout<<">> Error: "<<message<<std::endl;
    exit(1); 
}

double CCProtonPi0_NTupleAnalysis::GetFluxWeight(double Enu, int nuPDG)
{
    double flux_weight = frw->GetFluxCVWeight(Enu, nuPDG);
    return flux_weight;
}

std::vector<double> CCProtonPi0_NTupleAnalysis::GetFluxError(double Enu, int nuPDG)
{
    std::vector<double> flux_error = frw->GetFluxErrorWeights(Enu, nuPDG);
    return flux_error;
}

bool CCProtonPi0_NTupleAnalysis::IsEvent2p2h(int type)
{
    if ( type == 8 ) return true;
    else return false;
}

// This is the 2D Gaussian weight function for 2p2h Weights 
// Central Value= fit_2p2h_CV;
// +-1sigma = fit_2p2h_np, fit_2p2h_nn;
double CCProtonPi0_NTupleAnalysis::Get_2p2h_wgt(Double_t* neutrino_4P, Double_t* muon_4P, std::vector<double> fit_results)
{
    // Calculate q0 & q3
    double q0 = neutrino_4P[3] - muon_4P[3]; 
    double q3_x = neutrino_4P[0] - muon_4P[0]; 
    double q3_y = neutrino_4P[1] - muon_4P[1]; 
    double q3_z = neutrino_4P[2] - muon_4P[2]; 
    double q3 = HEP_Functions::calcMomentum(q3_x, q3_y, q3_z);  
    
    // Convert to GeV
    q0 = q0 * MeV_to_GeV;
    q3 = q3 * MeV_to_GeV;

    // 2D Gaussian Parameters
    double norm = fit_results[0];
    double meanq0 = fit_results[1];
    double meanq3 = fit_results[2];
    double sigmaq0 = fit_results[3];
    double sigmaq3 = fit_results[4];
    double corr = fit_results[5];

    double z =  std::pow((q0 - meanq0),2) / std::pow(sigmaq0,2) + 
                std::pow((q3 - meanq3),2) / std::pow(sigmaq3,2) -
                2*corr*(q0-meanq0)*(q3-meanq3)/ (sigmaq0 * sigmaq3);

    double ret = norm*exp( -0.5 * z / (1 - std::pow(corr,2)) );

    // The 2d gaussian +1: the function used for the "no scale down" fits
    ret = 1+ ret;
    
    //std::cout<<"2p2h Weight for q0 = "<<q0<<" q3 = "<<q3<<std::endl;
    //std::cout<<"\t"<<ret<<std::endl;
    return ret;
}

void CCProtonPi0_NTupleAnalysis::init2p2hFitResults()
{
    fit_2p2h_CV.push_back(18.5896);
    fit_2p2h_CV.push_back(0.256895);
    fit_2p2h_CV.push_back(0.509694);
    fit_2p2h_CV.push_back(0.0528243);
    fit_2p2h_CV.push_back(0.125003);
    fit_2p2h_CV.push_back(0.949053);
    
    fit_2p2h_np.push_back(19.2035);
    fit_2p2h_np.push_back(0.238834);
    fit_2p2h_np.push_back(0.494791);
    fit_2p2h_np.push_back(0.0557671);
    fit_2p2h_np.push_back(0.127039);
    fit_2p2h_np.push_back(0.889579);

    fit_2p2h_nn.push_back(78.0529);
    fit_2p2h_nn.push_back(0.284221);
    fit_2p2h_nn.push_back(0.52907);
    fit_2p2h_nn.push_back(0.0551371);
    fit_2p2h_nn.push_back(0.121867);
    fit_2p2h_nn.push_back(0.981287);
}

#endif

