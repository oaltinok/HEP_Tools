#ifndef CCProtonPi0_NTupleAnalysis_cpp
#define CCProtonPi0_NTupleanalysis_cpp

#include "CCProtonPi0_NTupleAnalysis.h"

using namespace PlotUtils;

// Initialize Constants
const std::string CCProtonPi0_NTupleAnalysis::version = "v2_80";

const double CCProtonPi0_NTupleAnalysis::data_POT = 3.33458e+20;
const double CCProtonPi0_NTupleAnalysis::mc_POT = 2.2185e+21;
const double CCProtonPi0_NTupleAnalysis::POT_ratio = data_POT/mc_POT;

const double CCProtonPi0_NTupleAnalysis::min_Enu = 1500; // MeV
const double CCProtonPi0_NTupleAnalysis::max_Enu = 20000; // MeV

const double CCProtonPi0_NTupleAnalysis::SENTINEL = -9.9;
const double CCProtonPi0_NTupleAnalysis::MeV_to_GeV = pow(10,-3);
const double CCProtonPi0_NTupleAnalysis::MeVSq_to_GeVSq = pow(10,-6);
const double CCProtonPi0_NTupleAnalysis::mSq_to_cmSq = pow(10,4);
const double CCProtonPi0_NTupleAnalysis::mm_to_cm = pow(10,-1);
const double CCProtonPi0_NTupleAnalysis::rad_to_deg = 180.0/M_PI;

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

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBands_Data(MnvHistoType* h)
{
    AddVertErrorBandAndFillWithCV_Flux(h);
    AddVertErrorBandAndFillWithCV_Genie(h);
    AddVertErrorBandAndFillWithCV_MuonTracking(h);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBands_Data<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBands_Data<MnvH2D>(MnvH2D* h);

    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_Flux(MnvHistoType* h)
{
    h->AddVertErrorBandAndFillWithCV("Flux",  n_universe);
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
    h->AddVertErrorBandAndFillWithCV("GENIE_MaCCQEshape"       ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_MaNCEL"            ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_MaRES"             ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_MvRES"             ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_NormCCQE"          ,2);
    h->AddVertErrorBandAndFillWithCV("GENIE_NormCCRES"         ,2);
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
void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_MuonTracking(MnvHistoType* h)
{
    h->AddVertErrorBandAndFillWithCV("MuonTracking", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_MuonTracking<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_MuonTracking<MnvH2D>(MnvH2D* h);

// MC Error Bars
    template<class MnvHistoType>
void CCProtonPi0_NTupleAnalysis::AddVertErrorBands_MC(MnvHistoType* h)
{
    AddVertErrorBand_Flux(h);
    AddVertErrorBand_Genie(h);
    AddVertErrorBand_MuonTracking(h);
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
    h->AddVertErrorBand("GENIE_MaCCQEshape"       ,2);
    h->AddVertErrorBand("GENIE_MaNCEL"            ,2);
    h->AddVertErrorBand("GENIE_MaRES"             ,2);
    h->AddVertErrorBand("GENIE_MvRES"             ,2);
    h->AddVertErrorBand("GENIE_NormCCQE"          ,2);
    h->AddVertErrorBand("GENIE_NormCCRES"         ,2);
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
void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_MuonTracking(MnvHistoType* h)
{
    h->AddVertErrorBand("MuonTracking", 2);
}
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_MuonTracking<MnvH1D>(MnvH1D* h);
template void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_MuonTracking<MnvH2D>(MnvH2D* h);

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


#endif

