#ifndef CCProtonPi0_NTupleAnalysis_cpp
#define CCProtonPi0_NTupleanalysis_cpp

#include "CCProtonPi0_NTupleAnalysis.h"

using namespace PlotUtils;

// Initialize Constants
const std::string CCProtonPi0_NTupleAnalysis::version = "v2_60";
const double CCProtonPi0_NTupleAnalysis::SENTINEL = -9.9;
const double CCProtonPi0_NTupleAnalysis::MeV_to_GeV = pow(10,-3);
const double CCProtonPi0_NTupleAnalysis::MeVSq_to_GeVSq = pow(10,-6);
const double CCProtonPi0_NTupleAnalysis::mm_to_cm = pow(10,-1);
const double CCProtonPi0_NTupleAnalysis::rad_to_deg = 180.0/M_PI;

CCProtonPi0_NTupleAnalysis::CCProtonPi0_NTupleAnalysis()
{
    // Required for MINERvA Framework Classes
    ROOT::Cintex::Cintex::Enable();
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

void CCProtonPi0_NTupleAnalysis::AddVertErrorBands_Data(MnvH1D* h)
{
    AddVertErrorBandAndFillWithCV_Flux(h);
    AddVertErrorBandAndFillWithCV_Genie(h);
    AddVertErrorBandAndFillWithCV_MuonTracking(h);
}

void CCProtonPi0_NTupleAnalysis::AddVertErrorBands_Data(MnvH2D* h)
{
    AddVertErrorBandAndFillWithCV_Flux(h);
    AddVertErrorBandAndFillWithCV_Genie(h);
    AddVertErrorBandAndFillWithCV_MuonTracking(h);
}

void CCProtonPi0_NTupleAnalysis::AddVertErrorBands_MC(MnvH1D* h)
{
    AddVertErrorBand_Flux(h);
    AddVertErrorBand_Genie(h);
    AddVertErrorBand_MuonTracking(h);
}

void CCProtonPi0_NTupleAnalysis::AddVertErrorBands_MC(MnvH2D* h)
{
    AddVertErrorBand_Flux(h);
    AddVertErrorBand_Genie(h);
    AddVertErrorBand_MuonTracking(h);
}

void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_Flux(MnvH1D* h)
{
    h->AddVertErrorBandAndFillWithCV("Flux",  n_universe);
}

void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_Flux(MnvH2D* h)
{
    h->AddVertErrorBandAndFillWithCV("Flux",  n_universe);
}

void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_Genie(MnvH1D* h)
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

void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_Genie(MnvH2D* h)
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

void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_Flux(MnvH1D* h)
{
    h->AddVertErrorBand("Flux",  n_universe);
}

void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_Flux(MnvH2D* h)
{
    h->AddVertErrorBand("Flux",  n_universe);
}

void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_Genie(MnvH1D* h)
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

void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_Genie(MnvH2D* h)
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

void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_MuonTracking(MnvH1D* h)
{
    h->AddVertErrorBandAndFillWithCV("MuonTracking", 2);
}

void CCProtonPi0_NTupleAnalysis::AddVertErrorBandAndFillWithCV_MuonTracking(MnvH2D* h)
{
    h->AddVertErrorBandAndFillWithCV("MuonTracking", 2);
}

void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_MuonTracking(MnvH1D* h)
{
    h->AddVertErrorBand("MuonTracking", 2);
}

void CCProtonPi0_NTupleAnalysis::AddVertErrorBand_MuonTracking(MnvH2D* h)
{
    h->AddVertErrorBand("MuonTracking", 2);
}


#endif

