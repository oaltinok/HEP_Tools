#ifndef CCProtonPi0_NTupleAnalysis_cpp
#define CCProtonPi0_NTupleanalysis_cpp

#include "CCProtonPi0_NTupleAnalysis.h"

using namespace PlotUtils;

// Initialize Constants
const std::string CCProtonPi0_NTupleAnalysis::version = "v2_93";

const double CCProtonPi0_NTupleAnalysis::EPSILON = 1.0e-3; 

const double CCProtonPi0_NTupleAnalysis::data_POT = 3.33534e+20;
const double CCProtonPi0_NTupleAnalysis::mc_POT = 2.21773e+21; 
const double CCProtonPi0_NTupleAnalysis::POT_ratio = data_POT/mc_POT;

const double CCProtonPi0_NTupleAnalysis::min_Enu = 1500; // MeV
const double CCProtonPi0_NTupleAnalysis::max_Enu = 20000; // MeV
const double CCProtonPi0_NTupleAnalysis::max_W = 140000; // MeV

const double CCProtonPi0_NTupleAnalysis::SENTINEL = -9.9;
const double CCProtonPi0_NTupleAnalysis::MeV_to_GeV = pow(10,-3);
const double CCProtonPi0_NTupleAnalysis::MeVSq_to_GeVSq = pow(10,-6);
const double CCProtonPi0_NTupleAnalysis::mSq_to_cmSq = pow(10,4);
const double CCProtonPi0_NTupleAnalysis::mm_to_cm = pow(10,-1);
const double CCProtonPi0_NTupleAnalysis::rad_to_deg = 180.0/M_PI;

const double CCProtonPi0_NTupleAnalysis::muon_mass = 105.66;    // MeV
const double CCProtonPi0_NTupleAnalysis::pi0_mass = 134.98;     // MeV
const double CCProtonPi0_NTupleAnalysis::piplus_mass = 139.57;  // MeV
const double CCProtonPi0_NTupleAnalysis::proton_mass = 938.27;  // MeV
const double CCProtonPi0_NTupleAnalysis::neutron_mass = 939.57; // MeV

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
    processed_minerva13 = false;
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
    AddVertErrorBandAndFillWithCV_BckgConstraint(h);
    AddVertErrorBandAndFillWithCV_TargetMass(h);
    AddVertErrorBandAndFillWithCV_MuonTracking(h);
    AddVertErrorBandAndFillWithCV_ProtonTracking(h);
    AddVertErrorBandAndFillWithCV_NeutronResponse(h);
    AddVertErrorBandAndFillWithCV_PionResponse(h);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBands_Data<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBands_Data<MnvH2D>(MnvH2D* h);

// Truth Tree Only have GENIE and Flux Errors, others are handled as Data
    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBands_TruthTree(MnvHistoType* h)
{
    AddVertErrorBandAndFillWithCV_BckgConstraint(h);
    AddVertErrorBandAndFillWithCV_TargetMass(h);
    AddVertErrorBandAndFillWithCV_MuonTracking(h);
    AddVertErrorBandAndFillWithCV_ProtonTracking(h);
    AddVertErrorBandAndFillWithCV_NeutronResponse(h);
    AddVertErrorBandAndFillWithCV_PionResponse(h);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBands_TruthTree<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBands_TruthTree<MnvH2D>(MnvH2D* h);

// Flux Histogram have Flux Errors, others are handled as Data
    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBands_FluxHistogram(MnvHistoType* h)
{
    AddVertErrorBandAndFillWithCV_Genie(h);
    AddVertErrorBandAndFillWithCV_BckgConstraint(h);
    AddVertErrorBandAndFillWithCV_TargetMass(h);
    AddVertErrorBandAndFillWithCV_MuonTracking(h);
    AddVertErrorBandAndFillWithCV_ProtonTracking(h);
    AddVertErrorBandAndFillWithCV_NeutronResponse(h);
    AddVertErrorBandAndFillWithCV_PionResponse(h);
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
void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_BckgConstraint(MnvHistoType* h)
{
    h->AddVertErrorBandAndFillWithCV("BckgConstraint", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_BckgConstraint<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_BckgConstraint<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_TargetMass(MnvHistoType* h)
{
    h->AddVertErrorBandAndFillWithCV("TargetMass", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_TargetMass<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_TargetMass<MnvH2D>(MnvH2D* h);

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
void CCProtonPi0_NTupleAnalysis::AddLatErrorBands_TruthTree(MnvHistoType* h)
{
    AddLatErrorBands_Data(h);
}
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBands_TruthTree<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddLatErrorBands_TruthTree<MnvH2D>(MnvH2D* h);

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
    AddVertErrorBand_BckgConstraint(h);
    AddVertErrorBand_TargetMass(h);
    AddVertErrorBand_MuonTracking(h);
    AddVertErrorBand_ProtonTracking(h);
    AddVertErrorBand_NeutronResponse(h);
    AddVertErrorBand_PionResponse(h);
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
void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_BckgConstraint(MnvHistoType* h)
{
    h->AddVertErrorBand("BckgConstraint", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_BckgConstraint<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_BckgConstraint<MnvH2D>(MnvH2D* h);

     template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_TargetMass(MnvHistoType* h)
{
    h->AddVertErrorBand("TargetMass", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_TargetMass<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_TargetMass<MnvH2D>(MnvH2D* h);

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


std::string CCProtonPi0_NTupleAnalysis::GetPlaylist(const int run)
{
    std::string playlist;
    if (run >= 10200 && run <= 10249) playlist = "minerva1";
    else if (run >= 10250 && run <= 10254) playlist = "minerva7";
    else if (run >= 10255 && run <= 10259) playlist = "minerva9";
    else if (run >= 12200 && run <= 12209) playlist = "minerva13A";
    else if (run >= 12210 && run <= 12219) playlist = "minerva13B";
    else if (run >= 13200 && run <= 13299) playlist = "minerva13C";
    else if (run >= 14201 && run <= 14209) playlist = "minerva13D";
    else if (run >= 14210 && run <= 14229) playlist = "minerva13E";
    else std::cout<<"WARNING: NO Playlist Found!"<<std::endl;

    return playlist;
}

void CCProtonPi0_NTupleAnalysis::UpdateFluxReweighter(int run)
{
    std::string playlist = GetPlaylist(run);

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
    }else if (!processed_minerva13 && playlist.find("minerva13") != std::string::npos){
        std::cout<<"Playlist: minerva13"<<std::endl;
        ReInitFluxReweighter(FluxReweighter::minerva13);
        processed_minerva13 = true;
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
        std::cout.width(12); std::cout<<content;
        std::cout.width(12); std::cout<<content/nEntries*100;
        std::cout<<std::endl;
    }

    std::cout.width(12); std::cout<<"Total";
    std::cout.width(12); std::cout<<total;
    std::cout.width(12); std::cout<<total/nEntries*100;
    std::cout<<std::endl;

}

bool CCProtonPi0_NTupleAnalysis::IsWLow(double true_W)
{
    if (true_W <= max_W) return true;
    else return false;
}

bool CCProtonPi0_NTupleAnalysis::IsEnuInRange(double Enu)
{
    if (Enu >= min_Enu && Enu <= max_Enu) return true;
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



#endif

