#ifndef CCProtonPi0_Analyzer_Systematics_cpp
#define CCProtonPi0_Analyzer_Systematics_cpp

#include "CCProtonPi0_Analyzer.h"

using namespace std;

void CCProtonPi0_Analyzer::AddErrorBands_Data()
{
    AddVertErrorBands_Data(pi0.pi0_P_all);
    AddVertErrorBands_Data(pi0.pi0_KE_all);
    AddVertErrorBands_Data(pi0.pi0_theta_all);
    AddVertErrorBands_Data(muon.muon_P_all);
    AddVertErrorBands_Data(muon.muon_theta_all);
    AddVertErrorBands_Data(interaction.QSq_all);
    AddVertErrorBands_Data(interaction.Enu_all);
    AddVertErrorBands_Data(interaction.W_all);

    AddLatErrorBands_Data(pi0.pi0_P_all);
    AddLatErrorBands_Data(pi0.pi0_KE_all);
    AddLatErrorBands_Data(pi0.pi0_theta_all);
    AddLatErrorBands_Data(muon.muon_P_all);
    AddLatErrorBands_Data(muon.muon_theta_all);
    AddLatErrorBands_Data(interaction.QSq_all);
    AddLatErrorBands_Data(interaction.Enu_all);
    AddLatErrorBands_Data(interaction.W_all);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_Genie(MnvH1D* h, double var)
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
    h->FillVertErrorBand("GENIE_MaRES"             ,var, truth_genie_wgt_MaRES[2]            , truth_genie_wgt_MaRES[4]            , cvweight);
    h->FillVertErrorBand("GENIE_MvRES"             ,var, truth_genie_wgt_MvRES[2]            , truth_genie_wgt_MvRES[4]            , cvweight);
    //h->FillVertErrorBand("GENIE_NormCCQE"          ,var, truth_genie_wgt_NormCCQE[2]         , truth_genie_wgt_NormCCQE[4]         , cvweight);
    //h->FillVertErrorBand("GENIE_NormCCRES"         ,var, truth_genie_wgt_NormCCRES[2]        , truth_genie_wgt_NormCCRES[4]        , cvweight);
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

void CCProtonPi0_Analyzer::FillVertErrorBand_Genie_ByHand(MnvH1D* h, double var)
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
    FillVertErrorBand_ByHand(h, var, "GENIE_MaRES"             , truth_genie_wgt_MaRES[2]            , truth_genie_wgt_MaRES[4]            );
    FillVertErrorBand_ByHand(h, var, "GENIE_MvRES"             , truth_genie_wgt_MvRES[2]            , truth_genie_wgt_MvRES[4]            );
    //FillVertErrorBand_ByHand(h, var, "GENIE_NormCCQE"          , truth_genie_wgt_NormCCQE[2]         , truth_genie_wgt_NormCCQE[4]         );
    //FillVertErrorBand_ByHand(h, var, "GENIE_NormCCRES"         , truth_genie_wgt_NormCCRES[2]        , truth_genie_wgt_NormCCRES[4]        );
    FillVertErrorBand_ByHand(h, var, "GENIE_NormDISCC"         , truth_genie_wgt_NormDISCC[2]        , truth_genie_wgt_NormDISCC[4]        );
    FillVertErrorBand_ByHand(h, var, "GENIE_NormNCRES"         , truth_genie_wgt_NormNCRES[2]        , truth_genie_wgt_NormNCRES[4]        );
    FillVertErrorBand_ByHand(h, var, "GENIE_RDecBR1gamma"      , truth_genie_wgt_RDecBR1gamma[2]     , truth_genie_wgt_RDecBR1gamma[4]     );
    FillVertErrorBand_ByHand(h, var, "GENIE_Rvn1pi"            , truth_genie_wgt_Rvn1pi[2]           , truth_genie_wgt_Rvn1pi[4]           );
    FillVertErrorBand_ByHand(h, var, "GENIE_Rvn2pi"            , truth_genie_wgt_Rvn2pi[2]           , truth_genie_wgt_Rvn2pi[4]           );
    FillVertErrorBand_ByHand(h, var, "GENIE_Rvp1pi"            , truth_genie_wgt_Rvp1pi[2]           , truth_genie_wgt_Rvp1pi[4]           );
    FillVertErrorBand_ByHand(h, var, "GENIE_Rvp2pi"            , truth_genie_wgt_Rvp2pi[2]           , truth_genie_wgt_Rvp2pi[4]           );
    FillVertErrorBand_ByHand(h, var, "GENIE_Theta_Delta2Npi"   , truth_genie_wgt_Theta_Delta2Npi[2]  , truth_genie_wgt_Theta_Delta2Npi[4]  );
    FillVertErrorBand_ByHand(h, var, "GENIE_VecFFCCQEshape"    , truth_genie_wgt_VecFFCCQEshape[2]   , truth_genie_wgt_VecFFCCQEshape[4]   );
}

void CCProtonPi0_Analyzer::FillVertErrorBand_Genie(MnvH2D* h, double xval, double yval)
{
    h->FillVertErrorBand("GENIE_AGKYxF1pi"         ,xval, yval, truth_genie_wgt_AGKYxF1pi[2]        , truth_genie_wgt_AGKYxF1pi[4]        , cvweight);
    h->FillVertErrorBand("GENIE_AhtBY"             ,xval, yval, truth_genie_wgt_AhtBY[2]            , truth_genie_wgt_AhtBY[4]            , cvweight);
    h->FillVertErrorBand("GENIE_BhtBY"             ,xval, yval, truth_genie_wgt_BhtBY[2]            , truth_genie_wgt_BhtBY[4]            , cvweight);
    h->FillVertErrorBand("GENIE_CCQEPauliSupViaKF" ,xval, yval, truth_genie_wgt_CCQEPauliSupViaKF[2], truth_genie_wgt_CCQEPauliSupViaKF[4], cvweight);
    h->FillVertErrorBand("GENIE_CV1uBY"            ,xval, yval, truth_genie_wgt_CV1uBY[2]           , truth_genie_wgt_CV1uBY[4]           , cvweight);
    h->FillVertErrorBand("GENIE_CV2uBY"            ,xval, yval, truth_genie_wgt_CV2uBY[2]           , truth_genie_wgt_CV2uBY[4]           , cvweight);
    h->FillVertErrorBand("GENIE_EtaNCEL"           ,xval, yval, truth_genie_wgt_EtaNCEL[2]          , truth_genie_wgt_EtaNCEL[4]          , cvweight);
    h->FillVertErrorBand("GENIE_FrAbs_N"           ,xval, yval, truth_genie_wgt_FrAbs_N[2]          , truth_genie_wgt_FrAbs_N[4]          , cvweight);
    h->FillVertErrorBand("GENIE_FrAbs_pi"          ,xval, yval, truth_genie_wgt_FrAbs_pi[2]         , truth_genie_wgt_FrAbs_pi[4]         , cvweight);
    h->FillVertErrorBand("GENIE_FrCEx_N"           ,xval, yval, truth_genie_wgt_FrCEx_N[2]          , truth_genie_wgt_FrCEx_N[4]          , cvweight);
    h->FillVertErrorBand("GENIE_FrCEx_pi"          ,xval, yval, truth_genie_wgt_FrCEx_pi[2]         , truth_genie_wgt_FrCEx_pi[4]         , cvweight);
    h->FillVertErrorBand("GENIE_FrElas_N"          ,xval, yval, truth_genie_wgt_FrElas_N[2]         , truth_genie_wgt_FrElas_N[4]         , cvweight);
    h->FillVertErrorBand("GENIE_FrElas_pi"         ,xval, yval, truth_genie_wgt_FrElas_pi[2]        , truth_genie_wgt_FrElas_pi[4]        , cvweight);
    h->FillVertErrorBand("GENIE_FrInel_N"          ,xval, yval, truth_genie_wgt_FrInel_N[2]         , truth_genie_wgt_FrInel_N[4]         , cvweight);
    h->FillVertErrorBand("GENIE_FrInel_pi"         ,xval, yval, truth_genie_wgt_FrInel_pi[2]        , truth_genie_wgt_FrInel_pi[4]        , cvweight);
    h->FillVertErrorBand("GENIE_FrPiProd_N"        ,xval, yval, truth_genie_wgt_FrPiProd_N[2]       , truth_genie_wgt_FrPiProd_N[4]       , cvweight);
    h->FillVertErrorBand("GENIE_FrPiProd_pi"       ,xval, yval, truth_genie_wgt_FrPiProd_pi[2]      , truth_genie_wgt_FrPiProd_pi[4]      , cvweight);
    h->FillVertErrorBand("GENIE_MFP_N"             ,xval, yval, truth_genie_wgt_MFP_N[2]            , truth_genie_wgt_MFP_N[4]            , cvweight);
    h->FillVertErrorBand("GENIE_MFP_pi"            ,xval, yval, truth_genie_wgt_MFP_pi[2]           , truth_genie_wgt_MFP_pi[4]           , cvweight);
    h->FillVertErrorBand("GENIE_MaCCQE"            ,xval, yval, truth_genie_wgt_MaCCQE[2]           , truth_genie_wgt_MaCCQE[4]           , cvweight);
    //h->FillVertErrorBand("GENIE_MaCCQEshape"       ,xval, yval, truth_genie_wgt_MaCCQEshape[2]      , truth_genie_wgt_MaCCQEshape[4]      , cvweight);
    h->FillVertErrorBand("GENIE_MaNCEL"            ,xval, yval, truth_genie_wgt_MaNCEL[2]           , truth_genie_wgt_MaNCEL[4]           , cvweight);
    h->FillVertErrorBand("GENIE_MaRES"             ,xval, yval, truth_genie_wgt_MaRES[2]            , truth_genie_wgt_MaRES[4]            , cvweight);
    h->FillVertErrorBand("GENIE_MvRES"             ,xval, yval, truth_genie_wgt_MvRES[2]            , truth_genie_wgt_MvRES[4]            , cvweight);
    //h->FillVertErrorBand("GENIE_NormCCQE"          ,xval, yval, truth_genie_wgt_NormCCQE[2]         , truth_genie_wgt_NormCCQE[4]         , cvweight);
    //h->FillVertErrorBand("GENIE_NormCCRES"         ,xval, yval, truth_genie_wgt_NormCCRES[2]        , truth_genie_wgt_NormCCRES[4]        , cvweight);
    h->FillVertErrorBand("GENIE_NormDISCC"         ,xval, yval, truth_genie_wgt_NormDISCC[2]        , truth_genie_wgt_NormDISCC[4]        , cvweight);
    h->FillVertErrorBand("GENIE_NormNCRES"         ,xval, yval, truth_genie_wgt_NormNCRES[2]        , truth_genie_wgt_NormNCRES[4]        , cvweight);
    h->FillVertErrorBand("GENIE_RDecBR1gamma"      ,xval, yval, truth_genie_wgt_RDecBR1gamma[2]     , truth_genie_wgt_RDecBR1gamma[4]     , cvweight);
    h->FillVertErrorBand("GENIE_Rvn1pi"            ,xval, yval, truth_genie_wgt_Rvn1pi[2]           , truth_genie_wgt_Rvn1pi[4]           , cvweight);
    h->FillVertErrorBand("GENIE_Rvn2pi"            ,xval, yval, truth_genie_wgt_Rvn2pi[2]           , truth_genie_wgt_Rvn2pi[4]           , cvweight);
    h->FillVertErrorBand("GENIE_Rvp1pi"            ,xval, yval, truth_genie_wgt_Rvp1pi[2]           , truth_genie_wgt_Rvp1pi[4]           , cvweight);
    h->FillVertErrorBand("GENIE_Rvp2pi"            ,xval, yval, truth_genie_wgt_Rvp2pi[2]           , truth_genie_wgt_Rvp2pi[4]           , cvweight);
    h->FillVertErrorBand("GENIE_Theta_Delta2Npi"   ,xval, yval, truth_genie_wgt_Theta_Delta2Npi[2]  , truth_genie_wgt_Theta_Delta2Npi[4]  , cvweight);
    h->FillVertErrorBand("GENIE_VecFFCCQEshape"    ,xval, yval, truth_genie_wgt_VecFFCCQEshape[2]   , truth_genie_wgt_VecFFCCQEshape[4]   , cvweight);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_Genie_ByHand(MnvH2D* h, double xval, double yval)
{
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_AGKYxF1pi"         , truth_genie_wgt_AGKYxF1pi[2]        , truth_genie_wgt_AGKYxF1pi[4]        );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_AhtBY"             , truth_genie_wgt_AhtBY[2]            , truth_genie_wgt_AhtBY[4]            );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_BhtBY"             , truth_genie_wgt_BhtBY[2]            , truth_genie_wgt_BhtBY[4]            );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_CCQEPauliSupViaKF" , truth_genie_wgt_CCQEPauliSupViaKF[2], truth_genie_wgt_CCQEPauliSupViaKF[4]);
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_CV1uBY"            , truth_genie_wgt_CV1uBY[2]           , truth_genie_wgt_CV1uBY[4]           );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_CV2uBY"            , truth_genie_wgt_CV2uBY[2]           , truth_genie_wgt_CV2uBY[4]           );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_EtaNCEL"           , truth_genie_wgt_EtaNCEL[2]          , truth_genie_wgt_EtaNCEL[4]          );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_FrAbs_N"           , truth_genie_wgt_FrAbs_N[2]          , truth_genie_wgt_FrAbs_N[4]          );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_FrAbs_pi"          , truth_genie_wgt_FrAbs_pi[2]         , truth_genie_wgt_FrAbs_pi[4]         );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_FrCEx_N"           , truth_genie_wgt_FrCEx_N[2]          , truth_genie_wgt_FrCEx_N[4]          );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_FrCEx_pi"          , truth_genie_wgt_FrCEx_pi[2]         , truth_genie_wgt_FrCEx_pi[4]         );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_FrElas_N"          , truth_genie_wgt_FrElas_N[2]         , truth_genie_wgt_FrElas_N[4]         );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_FrElas_pi"         , truth_genie_wgt_FrElas_pi[2]        , truth_genie_wgt_FrElas_pi[4]        );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_FrInel_N"          , truth_genie_wgt_FrInel_N[2]         , truth_genie_wgt_FrInel_N[4]         );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_FrInel_pi"         , truth_genie_wgt_FrInel_pi[2]        , truth_genie_wgt_FrInel_pi[4]        );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_FrPiProd_N"        , truth_genie_wgt_FrPiProd_N[2]       , truth_genie_wgt_FrPiProd_N[4]       );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_FrPiProd_pi"       , truth_genie_wgt_FrPiProd_pi[2]      , truth_genie_wgt_FrPiProd_pi[4]      );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_MFP_N"             , truth_genie_wgt_MFP_N[2]            , truth_genie_wgt_MFP_N[4]            );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_MFP_pi"            , truth_genie_wgt_MFP_pi[2]           , truth_genie_wgt_MFP_pi[4]           );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_MaCCQE"            , truth_genie_wgt_MaCCQE[2]           , truth_genie_wgt_MaCCQE[4]           );
    //FillVertErrorBand_ByHand(h, xval, yval, "GENIE_MaCCQEshape"       , truth_genie_wgt_MaCCQEshape[2]      , truth_genie_wgt_MaCCQEshape[4]      );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_MaNCEL"            , truth_genie_wgt_MaNCEL[2]           , truth_genie_wgt_MaNCEL[4]           );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_MaRES"             , truth_genie_wgt_MaRES[2]            , truth_genie_wgt_MaRES[4]            );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_MvRES"             , truth_genie_wgt_MvRES[2]            , truth_genie_wgt_MvRES[4]            );
    //FillVertErrorBand_ByHand(h, xval, yval, "GENIE_NormCCQE"          , truth_genie_wgt_NormCCQE[2]         , truth_genie_wgt_NormCCQE[4]         );
    //FillVertErrorBand_ByHand(h, xval, yval, "GENIE_NormCCRES"         , truth_genie_wgt_NormCCRES[2]        , truth_genie_wgt_NormCCRES[4]        );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_NormDISCC"         , truth_genie_wgt_NormDISCC[2]        , truth_genie_wgt_NormDISCC[4]        );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_NormNCRES"         , truth_genie_wgt_NormNCRES[2]        , truth_genie_wgt_NormNCRES[4]        );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_RDecBR1gamma"      , truth_genie_wgt_RDecBR1gamma[2]     , truth_genie_wgt_RDecBR1gamma[4]     );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_Rvn1pi"            , truth_genie_wgt_Rvn1pi[2]           , truth_genie_wgt_Rvn1pi[4]           );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_Rvn2pi"            , truth_genie_wgt_Rvn2pi[2]           , truth_genie_wgt_Rvn2pi[4]           );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_Rvp1pi"            , truth_genie_wgt_Rvp1pi[2]           , truth_genie_wgt_Rvp1pi[4]           );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_Rvp2pi"            , truth_genie_wgt_Rvp2pi[2]           , truth_genie_wgt_Rvp2pi[4]           );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_Theta_Delta2Npi"   , truth_genie_wgt_Theta_Delta2Npi[2]  , truth_genie_wgt_Theta_Delta2Npi[4]  );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_VecFFCCQEshape"    , truth_genie_wgt_VecFFCCQEshape[2]   , truth_genie_wgt_VecFFCCQEshape[4]   );
}


void CCProtonPi0_Analyzer::FillHistogramWithVertErrors(vector<MnvH1D*> &hist, double var)
{
    if (m_isMC){
        // Always Fill hist[0]
        FillHistogramWithVertErrors(hist[0], var);

        // Fill Signal
        if (truth_isSignal){
            FillHistogramWithVertErrors(hist[1], var);
        }else{
            // Fill Background
            FillHistogramWithVertErrors(hist[2], var);

            // Fill Background Type
            int ind = GetBackgroundTypeInd();
            FillHistogramWithVertErrors(hist[ind], var);
        }
    }else{
        FillHistogram(hist[0], var);
    }
}

void CCProtonPi0_Analyzer::FillVertErrorBand_ByHand(MnvH1D* h, double var, std::string error_name, std::vector<double> errors)
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
        
        // wgt_bckg is universe_wgt / cv_wgt 
        double wgt_bckg = applyBckgConstraints_Unv ? GetBckgConstraint(error_name, i) : 1.0;
      
        const double applyWeight = cvweight * wgt_bckg;
        const double wgtU = errors[i]*applyWeight;
        err_hists[i]->AddBinContent( cvbin, wgtU );

        const double err = err_hists[i]->GetBinError(cvbin);
        const double newerr2 = err*err + wgtU*wgtU; 
        const double newerr = (0.<newerr2) ? sqrt(newerr2) : 0.;
        err_hists[i]->SetBinError( cvbin, newerr );
    }
}

void CCProtonPi0_Analyzer::FillVertErrorBand_ByHand(MnvH1D* h, double var, std::string error_name, double err_down, double err_up)
{
    std::vector<double> errors;
    errors.push_back(err_down);
    errors.push_back(err_up);

    FillVertErrorBand_ByHand(h, var, error_name, errors);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_ByHand(MnvH2D* h, double xval, double yval, std::string error_name, std::vector<double> errors)
{
    // Get a Pointer to Error Band
    MnvVertErrorBand2D* err_band =  h->GetVertErrorBand(error_name);

    // Get a Pointer to Histograms 
    std::vector<TH2D*> err_hists = err_band->GetHists();
   
    // Sanity Check
    if (err_hists.size() != errors.size()) {
        std::cout<<"WARNING! Can not Fill Vertical Error Band: "<<error_name<<std::endl;
        exit(1);
    }

    int cvbin = err_band->TH2D::Fill( xval, yval, cvweight );

    if( cvbin == -1 ){
        cvbin = err_band->FindBin( xval, yval );
    }

    for( unsigned int i = 0; i < err_hists.size(); ++i ){
        
        // wgt_bckg is universe_wgt / cv_wgt
        //double wgt_bckg = applyBckgConstraints_Unv ? GetBckgConstraint(error_name, i) : 1.0;

        double wgt_bckg = 1.0; // We fill MnvH2D only for Signal Events: No need to search for Bckg Constraints

        const double applyWeight = cvweight * wgt_bckg;
        const double wgtU = errors[i]*applyWeight;
        err_hists[i]->AddBinContent( cvbin, wgtU );

        const double err = err_hists[i]->GetBinError(cvbin);
        const double newerr2 = err*err + wgtU*wgtU; 
        const double newerr = (0.<newerr2) ? sqrt(newerr2) : 0.;
        err_hists[i]->SetBinError( cvbin, newerr );
    }
}

void CCProtonPi0_Analyzer::FillVertErrorBand_ByHand(MnvH2D* h, double xval, double yval, std::string error_name, double err_down, double err_up)
{
    std::vector<double> errors;
    errors.push_back(err_down);
    errors.push_back(err_up);

    FillVertErrorBand_ByHand(h, xval, yval, error_name, errors);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_NeutronResponse(MnvH1D* h, double var)
{
    double correctionErr = GetNeutronResponseErr();

    if(!counter2.isCounted && correctionErr != 0){
        counter2.count++; 
        counter2.isCounted = true;
    }
    h->FillVertErrorBand("NeutronResponse", var, 1-correctionErr, 1+correctionErr, cvweight);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_NeutronResponse_ByHand(MnvH1D* h, double var)
{
    double correctionErr = GetNeutronResponseErr();
    FillVertErrorBand_ByHand(h, var, "NeutronResponse", 1-correctionErr, 1+correctionErr);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_NeutronResponse(MnvH2D* h, double xval, double yval)
{
    double correctionErr = GetNeutronResponseErr();

    h->FillVertErrorBand("NeutronResponse", xval, yval, 1-correctionErr, 1+correctionErr, cvweight);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_NeutronResponse_ByHand(MnvH2D* h, double xval, double yval)
{
    double correctionErr = GetNeutronResponseErr();

    FillVertErrorBand_ByHand(h, xval, yval, "NeutronResponse", 1-correctionErr, 1+correctionErr);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_PionResponse(MnvH1D* h, double var)
{
    double correctionErr = 0.0;

    if (truth_isBckg_SingleChargedPion_ChargeExchanged){
        correctionErr = 0.5;
        if(!counter1.isCounted){
            counter1.count++; 
            counter1.isCounted = true;
        }
    }

    h->FillVertErrorBand("PionResponse", var, 1-correctionErr, 1+correctionErr, cvweight);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_PionResponse_ByHand(MnvH1D* h, double var)
{
    double correctionErr = 0.0;

    if (truth_isBckg_SingleChargedPion_ChargeExchanged){
        correctionErr = 0.5;
    }
   
    FillVertErrorBand_ByHand(h, var, "PionResponse", 1-correctionErr, 1+correctionErr);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_PionResponse(MnvH2D* h, double xval, double yval)
{
    double correctionErr = 0.0;

    if (truth_isBckg_SingleChargedPion_ChargeExchanged){
        correctionErr = 0.5;
    }

    h->FillVertErrorBand("PionResponse", xval, yval, 1-correctionErr, 1+correctionErr, cvweight);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_PionResponse_ByHand(MnvH2D* h, double xval, double yval)
{
    double correctionErr = 0.0;

    if (truth_isBckg_SingleChargedPion_ChargeExchanged){
        correctionErr = 0.5;
    }
    
    FillVertErrorBand_ByHand(h, xval, yval, "PionResponse", 1-correctionErr, 1+correctionErr);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_MuonTracking(MnvH1D* h, double var)
{
    double correctionErr = GetMINOSCorrectionErr();
    h->FillVertErrorBand("MuonTracking", var, 1-correctionErr, 1+correctionErr, cvweight);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_MuonTracking_ByHand(MnvH1D* h, double var)
{
    double correctionErr = GetMINOSCorrectionErr();
    FillVertErrorBand_ByHand(h, var, "MuonTracking", 1-correctionErr, 1+correctionErr);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_MuonTracking(MnvH2D* h, double xval, double yval)
{
    double correctionErr = GetMINOSCorrectionErr();
    h->FillVertErrorBand("MuonTracking", xval, yval, 1-correctionErr, 1+correctionErr, cvweight);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_MuonTracking_ByHand(MnvH2D* h, double xval, double yval)
{
    double correctionErr = GetMINOSCorrectionErr();
 
    FillVertErrorBand_ByHand(h, xval, yval, "MuonTracking", 1-correctionErr, 1+correctionErr);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_Flux(MnvH1D* h, double var)
{
    std::vector<double> flux_errors = GetFluxError();
    h->FillVertErrorBand("Flux",  var, &flux_errors[0],  cvweight, 1.0);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_Flux_ByHand(MnvH1D* h, double var)
{
    std::vector<double> flux_errors = GetFluxError();

    FillVertErrorBand_ByHand(h, var, "Flux", flux_errors);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_Flux(MnvH2D* h, double xval, double yval)
{
    std::vector<double> flux_errors = GetFluxError();
    h->FillVertErrorBand("Flux",  xval, yval,  &flux_errors[0],  cvweight, 1.0);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_Flux_ByHand(MnvH2D* h, double xval, double yval)
{
    std::vector<double> flux_errors = GetFluxError();

    FillVertErrorBand_ByHand(h, xval, yval, "Flux", flux_errors);
}

void CCProtonPi0_Analyzer::FillHistogramWithVertErrors(MnvH1D* hist, double var)
{
    // Fill CV Value
    hist->Fill(var, cvweight);

    // Fill Vertical Error Bands
    if (fillErrors_ByHand){
        FillVertErrorBand_Genie_ByHand(hist, var);
        FillVertErrorBand_Flux_ByHand(hist, var);
        FillVertErrorBand_MuonTracking_ByHand(hist, var);
        //FillVertErrorBand_PionResponse_ByHand(hist, var);
        //FillVertErrorBand_NeutronResponse_ByHand(hist, var);
        FillVertErrorBand_PionResponse(hist, var);
        FillVertErrorBand_NeutronResponse(hist, var);
    }else{
        FillVertErrorBand_Genie(hist, var);
        FillVertErrorBand_Flux(hist, var);
        FillVertErrorBand_MuonTracking(hist, var);
        FillVertErrorBand_PionResponse(hist, var);
        FillVertErrorBand_NeutronResponse(hist, var);
    }
}

void CCProtonPi0_Analyzer::FillHistogramWithVertErrors(MnvH2D* hist, double xval, double yval)
{
    // Fill CV Value
    hist->Fill(xval,yval, cvweight);

    // Fill Vertical Error Bands
    if (fillErrors_ByHand){
        FillVertErrorBand_Genie_ByHand(hist, xval, yval);
        FillVertErrorBand_Flux_ByHand(hist, xval, yval);
        FillVertErrorBand_MuonTracking_ByHand(hist, xval, yval);
        //FillVertErrorBand_PionResponse_ByHand(hist, xval, yval);
        //FillVertErrorBand_NeutronResponse_ByHand(hist, xval, yval);
        FillVertErrorBand_PionResponse(hist, xval, yval);
        FillVertErrorBand_NeutronResponse(hist, xval, yval);
    }else{
        FillVertErrorBand_Genie(hist, xval, yval);
        FillVertErrorBand_Flux(hist, xval, yval);
        FillVertErrorBand_MuonTracking(hist, xval, yval);
        FillVertErrorBand_PionResponse(hist, xval, yval);
        FillVertErrorBand_NeutronResponse(hist, xval, yval);
    }
}

void CCProtonPi0_Analyzer::initLateralErrorBandShifts(bool isModeReduce)
{
    // Fill No Random Shift Vector
    for (int i = 0; i < n_lateral_universes; ++i){
        no_random_shifts.push_back(0.0);
    }   

    // Fill EM Energy Scale Shifts
    double mc_1sigma = 0.013;
    double data_1sigma = 0.019; 
    double em_uncertainty = sqrt(mc_1sigma*mc_1sigma + data_1sigma*data_1sigma);
    em_random_shifts = RandNumGenerator.GetRandomShifts(em_uncertainty); // ~ Gaussian(0.0, em_uncertainty)

    if (!isModeReduce){
        // Fill Histograms for Normal Random Numbers & EM Random Shift
        //  Muon Momentum Shift changes event by event
        std::vector<double> normal_rand_numbers = RandNumGenerator.GetNormalRandomVector();
        for (unsigned int i = 0; i < normal_rand_numbers.size(); ++i){
            interaction.normal_rand_numbers->Fill(normal_rand_numbers[i]);
            interaction.em_shift_rand_numbers->Fill(em_random_shifts[i]);
        }
        ismuonP_shifts_filled = false;
    }else{
        ismuonP_shifts_filled = true;
    }

}

void CCProtonPi0_Analyzer::Calc_muonP_random_shifts()
{
    double muonP_uncertainty = muon_E_shift/muon_P;
    muonP_random_shifts.clear();
    muonP_random_shifts = RandNumGenerator.GetRandomShifts(muonP_uncertainty); // ~ Gaussian(0.0, muonP_uncertainty)

    // Fill MuonP shifts only for 1 event
    if (!ismuonP_shifts_filled){
        for (unsigned int i = 0; i < muonP_random_shifts.size(); ++i){
            interaction.muonP_shift_rand_numbers->Fill(muonP_random_shifts[i]);
        }
        ismuonP_shifts_filled = true;
    }
}

void CCProtonPi0_Analyzer::FillLatErrorBands_ByHand()
{
    FillLatErrorBand_EM_EnergyScale();
    FillLatErrorBand_MuonMomentum();
}

void CCProtonPi0_Analyzer::FillLatErrorBand_EM_EnergyScale()
{
    std::string err_name = "EM_EnergyScale";
    double reco_muon_theta = GetCorrectedMuonTheta();
    
    for (int i = 0; i < n_lateral_universes; ++i){

        // Pi0 Variables
        double pi0_invMass_i = (1.0 + em_random_shifts[i]) * pi0_invMass;
        double pi0_P_i = (1.0 + em_random_shifts[i]) * pi0_P;
        double pi0_E_i = sqrt(pi0_P_i*pi0_P_i + pi0_mass*pi0_mass);
        double pi0_KE_i = pi0_E_i - pi0_mass;
        double Enu_i = Calc_Enu_shifted(muon_E, pi0_E_i); // Use actual muon energy and shifted pi0 energy
        double QSq_i = Calc_QSq(Enu_i, muon_E, muon_P, reco_muon_theta); // Use actual muon energy
        double WSq_i = Calc_WSq(Enu_i, QSq_i, muon_E); // Use actual muon energy
        double W_i = (WSq_i > 0) ? sqrt(WSq_i) : -1; 

        double pi0_P_shift = (pi0_P_i - pi0_P) * MeV_to_GeV;
        double pi0_KE_shift = (pi0_KE_i - pi0_KE) * MeV_to_GeV;
        double Enu_shift = (Enu_i - m_Enu) * MeV_to_GeV;
        double QSq_shift = (QSq_i - m_QSq) * MeVSq_to_GeVSq;
        double W_shift = (W_i - m_W) * MeV_to_GeV;

        bool PassedCuts = IsEnuInRange(Enu_i) && IsInvMassInRange(pi0_invMass_i);
        //PassedCuts = true; // For Comparing with Auto Fill
        
        if (PassedCuts){
            //-----------------------------------------------------------------
            // pi0_P
            //-----------------------------------------------------------------
            FillLatErrorBand_SingleUniverse(pi0.pi0_P_mc_reco_all, err_name, i, pi0_P * MeV_to_GeV, pi0_P_shift);
            if (truth_isSignal){
                FillLatErrorBand_SingleUniverse(pi0.pi0_P_mc_truth_signal, err_name, i, truth_pi0_P * MeV_to_GeV, 0.0);
                FillLatErrorBand_SingleUniverse(pi0.pi0_P_mc_reco_signal, err_name, i, pi0_P * MeV_to_GeV, pi0_P_shift);
                FillLatErrorBand_SingleUniverse(pi0.pi0_P_response, err_name, i, pi0_P * MeV_to_GeV, truth_pi0_P * MeV_to_GeV, pi0_P_shift, 0.0);
            }else{
                FillLatErrorBand_SingleUniverse(pi0.pi0_P_mc_reco_bckg, err_name, i, pi0_P * MeV_to_GeV, pi0_P_shift);
            }
            //-----------------------------------------------------------------
            // pi0_KE
            //-----------------------------------------------------------------
            FillLatErrorBand_SingleUniverse(pi0.pi0_KE_mc_reco_all, err_name, i, pi0_KE * MeV_to_GeV, pi0_KE_shift);
            if (truth_isSignal){
                FillLatErrorBand_SingleUniverse(pi0.pi0_KE_mc_truth_signal, err_name, i, truth_pi0_KE * MeV_to_GeV, 0.0);
                FillLatErrorBand_SingleUniverse(pi0.pi0_KE_mc_reco_signal, err_name, i, pi0_KE * MeV_to_GeV, pi0_KE_shift);
                FillLatErrorBand_SingleUniverse(pi0.pi0_KE_response, err_name, i, pi0_KE * MeV_to_GeV, truth_pi0_KE * MeV_to_GeV, pi0_KE_shift, 0.0);
            }else{
                FillLatErrorBand_SingleUniverse(pi0.pi0_KE_mc_reco_bckg, err_name, i, pi0_KE * MeV_to_GeV, pi0_KE_shift);
            }
            //-----------------------------------------------------------------
            // pi0_theta
            //-----------------------------------------------------------------
            FillLatErrorBand_SingleUniverse(pi0.pi0_theta_mc_reco_all, err_name, i, pi0_theta * TMath::RadToDeg(), 0.0);
            if (truth_isSignal){
                FillLatErrorBand_SingleUniverse(pi0.pi0_theta_mc_truth_signal, err_name, i, truth_pi0_theta * TMath::RadToDeg(), 0.0);
                FillLatErrorBand_SingleUniverse(pi0.pi0_theta_mc_reco_signal, err_name, i, pi0_theta * TMath::RadToDeg(), 0.0);
                FillLatErrorBand_SingleUniverse(pi0.pi0_theta_response, err_name, i, pi0_theta * TMath::RadToDeg(), truth_pi0_theta * TMath::RadToDeg(), 0.0, 0.0);
            }else{
                FillLatErrorBand_SingleUniverse(pi0.pi0_theta_mc_reco_bckg, err_name, i, pi0_theta * TMath::RadToDeg(), 0.0);
            }
            //-----------------------------------------------------------------
            // muon_P
            //-----------------------------------------------------------------
            FillLatErrorBand_SingleUniverse(muon.muon_P_mc_reco_all, err_name, i, muon_P * MeV_to_GeV, 0.0);
            if (truth_isSignal){
                FillLatErrorBand_SingleUniverse(muon.muon_P_mc_truth_signal, err_name, i, truth_muon_P * MeV_to_GeV, 0.0);
                FillLatErrorBand_SingleUniverse(muon.muon_P_mc_reco_signal, err_name, i, muon_P * MeV_to_GeV, 0.0);
                FillLatErrorBand_SingleUniverse(muon.muon_P_response, err_name, i, muon_P * MeV_to_GeV, truth_muon_P * MeV_to_GeV, 0.0, 0.0);
            }else{
                FillLatErrorBand_SingleUniverse(muon.muon_P_mc_reco_bckg, err_name, i, muon_P * MeV_to_GeV, 0.0);
            }
            //-----------------------------------------------------------------
            // muon_theta
            //-----------------------------------------------------------------
            FillLatErrorBand_SingleUniverse(muon.muon_theta_mc_reco_all, err_name, i, muon_theta * TMath::RadToDeg(), 0.0);
            if (truth_isSignal){
                FillLatErrorBand_SingleUniverse(muon.muon_theta_mc_truth_signal, err_name, i, truth_muon_theta * TMath::RadToDeg(), 0.0);
                FillLatErrorBand_SingleUniverse(muon.muon_theta_mc_reco_signal, err_name, i, muon_theta * TMath::RadToDeg(), 0.0);
                FillLatErrorBand_SingleUniverse(muon.muon_theta_response, err_name, i, muon_theta * TMath::RadToDeg(), truth_muon_theta * TMath::RadToDeg(), 0.0, 0.0);
            }else{
                FillLatErrorBand_SingleUniverse(muon.muon_theta_mc_reco_bckg, err_name, i, muon_theta * TMath::RadToDeg(), 0.0);
            }
            //-----------------------------------------------------------------
            // QSq
            //-----------------------------------------------------------------
            FillLatErrorBand_SingleUniverse(interaction.QSq_mc_reco_all, err_name, i, m_QSq * MeVSq_to_GeVSq, QSq_shift);
            if (truth_isSignal){
                FillLatErrorBand_SingleUniverse(interaction.QSq_mc_truth_signal, err_name, i, m_QSq_Truth * MeVSq_to_GeVSq, 0.0);
                FillLatErrorBand_SingleUniverse(interaction.QSq_mc_reco_signal, err_name, i, m_QSq * MeVSq_to_GeVSq, QSq_shift);
                FillLatErrorBand_SingleUniverse(interaction.QSq_response, err_name, i, m_QSq * MeVSq_to_GeVSq, m_QSq_Truth * MeVSq_to_GeVSq, QSq_shift, 0.0);
            }else{
                FillLatErrorBand_SingleUniverse(interaction.QSq_mc_reco_bckg, err_name, i, m_QSq * MeVSq_to_GeVSq, QSq_shift);
            }
            //-----------------------------------------------------------------
            // W
            //-----------------------------------------------------------------
            FillLatErrorBand_SingleUniverse(interaction.W_mc_reco_all, err_name, i, m_W * MeV_to_GeV, W_shift);
            if (truth_isSignal){
                FillLatErrorBand_SingleUniverse(interaction.W_mc_truth_signal, err_name, i, m_W_Truth * MeV_to_GeV, 0.0);
                FillLatErrorBand_SingleUniverse(interaction.W_mc_reco_signal, err_name, i, m_W * MeV_to_GeV, W_shift);
                FillLatErrorBand_SingleUniverse(interaction.W_response, err_name, i, m_W * MeV_to_GeV, m_W_Truth * MeV_to_GeV, W_shift, 0.0);
            }else{
                FillLatErrorBand_SingleUniverse(interaction.W_mc_reco_bckg, err_name, i, m_W * MeV_to_GeV, W_shift);
            }
            //-----------------------------------------------------------------
            // Enu
            //-----------------------------------------------------------------
            FillLatErrorBand_SingleUniverse(interaction.Enu_mc_reco_all, err_name, i, m_Enu * MeV_to_GeV, Enu_shift);
            if (truth_isSignal){
                FillLatErrorBand_SingleUniverse(interaction.Enu_mc_truth_signal, err_name, i, m_Enu_Truth * MeV_to_GeV, 0.0);
                FillLatErrorBand_SingleUniverse(interaction.Enu_mc_reco_signal, err_name, i, m_Enu * MeV_to_GeV, Enu_shift);
                FillLatErrorBand_SingleUniverse(interaction.Enu_response, err_name, i, m_Enu * MeV_to_GeV, m_Enu_Truth * MeV_to_GeV, Enu_shift, 0.0);
            }else{
                FillLatErrorBand_SingleUniverse(interaction.Enu_mc_reco_bckg, err_name, i, m_Enu * MeV_to_GeV, Enu_shift);
            }
        }else{
           //counter1.count++;
        }
    }
}

void CCProtonPi0_Analyzer::FillLatErrorBand_MuonMomentum()
{
    std::string err_name = "MuonMomentum";
    double reco_muon_theta = GetCorrectedMuonTheta();
   
    // For each event muonP_random_shifts is different
    Calc_muonP_random_shifts();
    for (int i = 0; i < n_lateral_universes; ++i){

        // Muon Variables
        double muon_E_i = (1.0 + muonP_random_shifts[i]) * muon_E;
        double muon_P_i = (1.0 + muonP_random_shifts[i]) * muon_P;
        double Enu_i = Calc_Enu_shifted(muon_E_i, pi0_E); // Use shifted muon energy and actual pi0 energy
        double QSq_i = Calc_QSq(Enu_i, muon_E_i, muon_P_i, reco_muon_theta); // Use shifted muon energy
        double WSq_i = Calc_WSq(Enu_i, QSq_i, muon_E_i); // Use shifted muon energy
        double W_i = (WSq_i > 0) ? sqrt(WSq_i) : -1; 

        double muon_P_shift = (muon_P_i - muon_P) * MeV_to_GeV;
        double Enu_shift = (Enu_i - m_Enu) * MeV_to_GeV;
        double QSq_shift = (QSq_i - m_QSq) * MeVSq_to_GeVSq;
        double W_shift = (W_i - m_W) * MeV_to_GeV;

        bool PassedCuts = IsEnuInRange(Enu_i);

        if (PassedCuts){
            //-----------------------------------------------------------------
            // muon_P
            //-----------------------------------------------------------------
            FillLatErrorBand_SingleUniverse(muon.muon_P_mc_reco_all, err_name, i, muon_P * MeV_to_GeV, muon_P_shift);
            if (truth_isSignal){
                FillLatErrorBand_SingleUniverse(muon.muon_P_mc_truth_signal, err_name, i, truth_muon_P * MeV_to_GeV, 0.0);
                FillLatErrorBand_SingleUniverse(muon.muon_P_mc_reco_signal, err_name, i, muon_P * MeV_to_GeV, muon_P_shift);
                FillLatErrorBand_SingleUniverse(muon.muon_P_response, err_name, i, muon_P * MeV_to_GeV, truth_muon_P * MeV_to_GeV, muon_P_shift, 0.0);
            }else{
                FillLatErrorBand_SingleUniverse(muon.muon_P_mc_reco_bckg, err_name, i, muon_P * MeV_to_GeV, muon_P_shift);
            }
            //-----------------------------------------------------------------
            // muon_theta
            //-----------------------------------------------------------------
            FillLatErrorBand_SingleUniverse(muon.muon_theta_mc_reco_all, err_name, i, reco_muon_theta * TMath::RadToDeg(), 0.0);
            if (truth_isSignal){
                FillLatErrorBand_SingleUniverse(muon.muon_theta_mc_truth_signal, err_name, i, truth_muon_theta * TMath::RadToDeg(), 0.0);
                FillLatErrorBand_SingleUniverse(muon.muon_theta_mc_reco_signal, err_name, i, reco_muon_theta * TMath::RadToDeg(), 0.0);
                FillLatErrorBand_SingleUniverse(muon.muon_theta_response, err_name, i, reco_muon_theta * TMath::RadToDeg(), truth_muon_theta * TMath::RadToDeg(), 0.0, 0.0);
            }else{
                FillLatErrorBand_SingleUniverse(muon.muon_theta_mc_reco_bckg, err_name, i, reco_muon_theta * TMath::RadToDeg(), 0.0);
            }
            //-----------------------------------------------------------------
            // pi0_P
            //-----------------------------------------------------------------
            FillLatErrorBand_SingleUniverse(pi0.pi0_P_mc_reco_all, err_name, i, pi0_P * MeV_to_GeV, 0.0);
            if (truth_isSignal){
                FillLatErrorBand_SingleUniverse(pi0.pi0_P_mc_truth_signal, err_name, i, truth_pi0_P * MeV_to_GeV, 0.0);
                FillLatErrorBand_SingleUniverse(pi0.pi0_P_mc_reco_signal, err_name, i, pi0_P * MeV_to_GeV, 0.0);
                FillLatErrorBand_SingleUniverse(pi0.pi0_P_response, err_name, i, pi0_P * MeV_to_GeV, truth_pi0_P * MeV_to_GeV, 0.0, 0.0);
            }else{
                FillLatErrorBand_SingleUniverse(pi0.pi0_P_mc_reco_bckg, err_name, i, pi0_P * MeV_to_GeV, 0.0);
            }
            //-----------------------------------------------------------------
            // pi0_KE
            //-----------------------------------------------------------------
            FillLatErrorBand_SingleUniverse(pi0.pi0_KE_mc_reco_all, err_name, i, pi0_KE * MeV_to_GeV, 0.0);
            if (truth_isSignal){
                FillLatErrorBand_SingleUniverse(pi0.pi0_KE_mc_truth_signal, err_name, i, truth_pi0_KE * MeV_to_GeV, 0.0);
                FillLatErrorBand_SingleUniverse(pi0.pi0_KE_mc_reco_signal, err_name, i, pi0_KE * MeV_to_GeV, 0.0);
                FillLatErrorBand_SingleUniverse(pi0.pi0_KE_response, err_name, i, pi0_KE * MeV_to_GeV, truth_pi0_KE * MeV_to_GeV, 0.0, 0.0);
            }else{
                FillLatErrorBand_SingleUniverse(pi0.pi0_KE_mc_reco_bckg, err_name, i, pi0_KE * MeV_to_GeV, 0.0);
            }
            //-----------------------------------------------------------------
            // pi0_theta
            //-----------------------------------------------------------------
            FillLatErrorBand_SingleUniverse(pi0.pi0_theta_mc_reco_all, err_name, i, pi0_theta * TMath::RadToDeg(), 0.0);
            if (truth_isSignal){
                FillLatErrorBand_SingleUniverse(pi0.pi0_theta_mc_truth_signal, err_name, i, truth_pi0_theta * TMath::RadToDeg(), 0.0);
                FillLatErrorBand_SingleUniverse(pi0.pi0_theta_mc_reco_signal, err_name, i, pi0_theta * TMath::RadToDeg(), 0.0);
                FillLatErrorBand_SingleUniverse(pi0.pi0_theta_response, err_name, i, pi0_theta * TMath::RadToDeg(), truth_pi0_theta * TMath::RadToDeg(), 0.0, 0.0);
            }else{
                FillLatErrorBand_SingleUniverse(pi0.pi0_theta_mc_reco_bckg, err_name, i, pi0_theta * TMath::RadToDeg(), 0.0);
            }
            //-----------------------------------------------------------------
            // QSq
            //-----------------------------------------------------------------
            FillLatErrorBand_SingleUniverse(interaction.QSq_mc_reco_all, err_name, i, m_QSq * MeVSq_to_GeVSq, QSq_shift);
            if (truth_isSignal){
                FillLatErrorBand_SingleUniverse(interaction.QSq_mc_truth_signal, err_name, i, m_QSq_Truth * MeVSq_to_GeVSq, 0.0);
                FillLatErrorBand_SingleUniverse(interaction.QSq_mc_reco_signal, err_name, i, m_QSq * MeVSq_to_GeVSq, QSq_shift);
                FillLatErrorBand_SingleUniverse(interaction.QSq_response, err_name, i, m_QSq * MeVSq_to_GeVSq, m_QSq_Truth * MeVSq_to_GeVSq, QSq_shift, 0.0);
            }else{
                FillLatErrorBand_SingleUniverse(interaction.QSq_mc_reco_bckg, err_name, i, m_QSq * MeVSq_to_GeVSq, QSq_shift);
            }
            //-----------------------------------------------------------------
            // W
            //-----------------------------------------------------------------
            FillLatErrorBand_SingleUniverse(interaction.W_mc_reco_all, err_name, i, m_W * MeV_to_GeV, W_shift);
            if (truth_isSignal){
                FillLatErrorBand_SingleUniverse(interaction.W_mc_truth_signal, err_name, i, m_W_Truth * MeV_to_GeV, 0.0);
                FillLatErrorBand_SingleUniverse(interaction.W_mc_reco_signal, err_name, i, m_W * MeV_to_GeV, W_shift);
                FillLatErrorBand_SingleUniverse(interaction.W_response, err_name, i, m_W * MeV_to_GeV, m_W_Truth * MeV_to_GeV, W_shift, 0.0);
            }else{
                FillLatErrorBand_SingleUniverse(interaction.W_mc_reco_bckg, err_name, i, m_W * MeV_to_GeV, W_shift);
            }
            //-----------------------------------------------------------------
            // Enu
            //-----------------------------------------------------------------
            FillLatErrorBand_SingleUniverse(interaction.Enu_mc_reco_all, err_name, i, m_Enu * MeV_to_GeV, Enu_shift);
            if (truth_isSignal){
                FillLatErrorBand_SingleUniverse(interaction.Enu_mc_truth_signal, err_name, i, m_Enu_Truth * MeV_to_GeV, 0.0);
                FillLatErrorBand_SingleUniverse(interaction.Enu_mc_reco_signal, err_name, i, m_Enu * MeV_to_GeV, Enu_shift);
                FillLatErrorBand_SingleUniverse(interaction.Enu_response, err_name, i, m_Enu * MeV_to_GeV, m_Enu_Truth * MeV_to_GeV, Enu_shift, 0.0);
            }else{
                FillLatErrorBand_SingleUniverse(interaction.Enu_mc_reco_bckg, err_name, i, m_Enu * MeV_to_GeV, Enu_shift);
            }
        }else{
           //counter2.count++;
        }
    }
}

void CCProtonPi0_Analyzer::FillLatErrorBand_EM_EnergyScale_invMass()
{
    std::string err_name = "EM_EnergyScale";
    
    for (int i = 0; i < n_lateral_universes; ++i){

        // EM Energy Dependent Variables
        double gamma1_E_i = (1.0 + em_random_shifts[i]) * gamma1_E;
        double gamma2_E_i = (1.0 + em_random_shifts[i]) * gamma2_E;
        double pi0_P_i = (1.0 + em_random_shifts[i]) * pi0_P;
        double pi0_E_i = sqrt(pi0_P_i*pi0_P_i + pi0_mass*pi0_mass);
        double Enu_i = Calc_Enu_shifted(muon_E, pi0_E_i); // Use actual muon energy and shifted pi0 energy

        bool PassedCuts = IsEnuInRange(Enu_i) && !IsOpeningAngleSmallAndEnergyLow(gamma1_E_i, gamma2_E_i);
        //PassedCuts = true; // For testing with Auto Fill

        if (PassedCuts){
            double pi0_invMass_i = (1.0 + em_random_shifts[i]) * pi0_invMass;
            double pi0_invMass_shift = pi0_invMass_i - pi0_invMass;
            FillLatErrorBand_SingleUniverse(cutList.invMass_mc_reco_all, err_name, i, pi0_invMass, pi0_invMass_shift);
            if (truth_isSignal){
                FillLatErrorBand_SingleUniverse(cutList.invMass_mc_reco_signal, err_name, i, pi0_invMass, pi0_invMass_shift);
            }else{
                FillLatErrorBand_SingleUniverse(cutList.invMass_mc_reco_bckg, err_name, i, pi0_invMass, pi0_invMass_shift);
            }
        }else{
           //counter1.count++;
        }
    }
}

void CCProtonPi0_Analyzer::FillLatErrorBand_MuonMomentum_invMass()
{
    std::string err_name = "MuonMomentum";
    
    // For each event muonP_random_shifts is different
    Calc_muonP_random_shifts();
 
    for (int i = 0; i < n_lateral_universes; ++i){

        // Muon Variables
        double muon_E_i = (1.0 + muonP_random_shifts[i]) * muon_E;
        double Enu_i = Calc_Enu_shifted(muon_E_i, pi0_E); // Use shifted muon energy and actual pi0 energy

        bool PassedCuts = IsEnuInRange(Enu_i);

        if (PassedCuts){
            FillLatErrorBand_SingleUniverse(cutList.invMass_mc_reco_all, err_name, i, pi0_invMass, 0.0);
            if (truth_isSignal){
                FillLatErrorBand_SingleUniverse(cutList.invMass_mc_reco_signal, err_name, i, pi0_invMass, 0.0);
            }else{
                FillLatErrorBand_SingleUniverse(cutList.invMass_mc_reco_bckg, err_name, i, pi0_invMass, 0.0);
            }
        }else{
           //counter2.count++;
        }
    }
}

void CCProtonPi0_Analyzer::FillLatErrorBand_EM_EnergyScale_SideBand_invMass()
{
    std::string err_name = "EM_EnergyScale";
    
    for (int i = 0; i < n_lateral_universes; ++i){

        // EM Energy Dependent Variables
        double gamma1_E_i = (1.0 + em_random_shifts[i]) * gamma1_E;
        double gamma2_E_i = (1.0 + em_random_shifts[i]) * gamma2_E;
        double pi0_P_i = (1.0 + em_random_shifts[i]) * pi0_P;
        double pi0_E_i = sqrt(pi0_P_i*pi0_P_i + pi0_mass*pi0_mass);
        double Enu_i = Calc_Enu_shifted(muon_E, pi0_E_i); // Use actual muon energy and shifted pi0 energy

        bool PassedCuts = IsEnuInRange(Enu_i) && !IsOpeningAngleSmallAndEnergyLow(gamma1_E_i, gamma2_E_i);

        if (PassedCuts){
            double pi0_invMass_i = (1.0 + em_random_shifts[i]) * pi0_invMass;
            double pi0_invMass_shift = pi0_invMass_i - pi0_invMass;
            
            // Fill All Events on [0]
            FillLatErrorBand_SingleUniverse(cutList.hCut_pi0invMass[0], err_name, i, pi0_invMass, pi0_invMass_shift);
            if (truth_isSignal){
                // Fill All Signal on [1]
                FillLatErrorBand_SingleUniverse(cutList.hCut_pi0invMass[1], err_name, i, pi0_invMass, pi0_invMass_shift);
            }else{
                // Fill All Background on [2]
                FillLatErrorBand_SingleUniverse(cutList.hCut_pi0invMass[2], err_name, i, pi0_invMass, pi0_invMass_shift);

                // Fill Background Type on [ind]
                int ind = GetBackgroundTypeInd();
                FillLatErrorBand_SingleUniverse(cutList.hCut_pi0invMass[ind], err_name, i, pi0_invMass, pi0_invMass_shift);
            }
        }
    }
}

void CCProtonPi0_Analyzer::FillLatErrorBand_MuonMomentum_SideBand_invMass()
{
    std::string err_name = "MuonMomentum";
    
    // For each event muonP_random_shifts is different
    Calc_muonP_random_shifts();
 
    for (int i = 0; i < n_lateral_universes; ++i){

        // Muon Variables
        double muon_E_i = (1.0 + muonP_random_shifts[i]) * muon_E;
        double Enu_i = Calc_Enu_shifted(muon_E_i, pi0_E); // Use shifted muon energy and actual pi0 energy

        bool PassedCuts = IsEnuInRange(Enu_i);

        if (PassedCuts){
            // Fill All Events on [0]
            FillLatErrorBand_SingleUniverse(cutList.hCut_pi0invMass[0], err_name, i, pi0_invMass, 0.0);
            if (truth_isSignal){
                // Fill Signal Events on [1]
                FillLatErrorBand_SingleUniverse(cutList.hCut_pi0invMass[1], err_name, i, pi0_invMass, 0.0);
            }else{
                // Fill All Background on [2]
                FillLatErrorBand_SingleUniverse(cutList.hCut_pi0invMass[2], err_name, i, pi0_invMass, 0.0);
            
                // Fill Background Type on [ind]
                int ind = GetBackgroundTypeInd();
                FillLatErrorBand_SingleUniverse(cutList.hCut_pi0invMass[ind], err_name, i, pi0_invMass, 0.0);
            }
        }
    }
}

void CCProtonPi0_Analyzer::FillLatErrorBand_SingleUniverse(MnvH1D* hist, std::string err_name, int unv, double var, double shift)
{
    // Get a Pointer to Error Band
    MnvLatErrorBand* err_band =  hist->GetLatErrorBand(err_name);

    // Get a Pointer to Histograms 
    std::vector<TH1D*> err_hists = err_band->GetHists();
  
    // Fill Error Band Base Histogram with Default cvweight
    // Fill Only once with Universe 0
    if (unv == 0){ 
        err_band->TH1D::Fill(var, cvweight);
    }

    int cvbin = err_band->FindBin(var);

    const int nbins = err_band->GetNbinsX();
    const double cvBinLowEdge  = (cvbin < err_band->GetNbinsX()+1) ? err_band->GetBinLowEdge( cvbin )  : err_band->GetBinLowEdge(cvbin-1) + err_band->GetBinWidth(cvbin-1);
    const double cvBinHighEdge = (cvbin < err_band->GetNbinsX()+1) ? err_band->GetBinLowEdge( cvbin+1) : cvBinLowEdge + err_band->GetBinWidth(cvbin-1);
    
    // Do not Fill Error Band if Shift is NOT Physical
    if( MnvHist::IsNotPhysicalShift(shift) ) return;

    const double shiftVal = var + shift;
    int bin = cvbin;
    if (shiftVal < var && shiftVal < cvBinLowEdge && bin>0)  {
        for( ; bin != 0; --bin )
            if( err_band->GetBinLowEdge( bin ) < shiftVal )
                break;
    }  
    else if (shiftVal > var && shiftVal > cvBinHighEdge && bin < err_band->GetNbinsX()+1) {
        for( ; bin != nbins+1; ++bin )
            if( shiftVal <  (err_band->GetBinLowEdge( bin ) + err_band->GetBinWidth( bin ) ) )
                break;
    }      

    // wgt_bckg is universe_wgt / cv_wgt
    double wgt_bckg = applyBckgConstraints_Unv ? GetBckgConstraint(err_name, unv) : 1.0;
    double wgtU = cvweight * wgt_bckg;
    err_hists[unv]->AddBinContent( bin, wgtU );

    const double err = err_hists[unv]->GetBinError(bin);
    const double newerr2 = err*err + wgtU*wgtU;
    const double newerr = (0.<newerr2) ? sqrt(newerr2) : 0.;
    err_hists[unv]->SetBinError( bin, newerr );

}

void CCProtonPi0_Analyzer::FillLatErrorBand_SingleUniverse(MnvH2D* hist, std::string err_name, int unv, double xval, double yval, double x_shift, double y_shift)
{
    // Get a Pointer to Error Band
    MnvLatErrorBand2D* err_band =  hist->GetLatErrorBand(err_name);

    // Get a Pointer to Histograms 
    std::vector<TH2D*> err_hists = err_band->GetHists();
  
    // Fill Error Band Base Histogram with Default cvweight
    // Fill Only once with Universe 0
    if (unv == 0){ 
        err_band->TH2D::Fill(xval, yval, cvweight);
    }

    // Do not Fill Error Band if Shift is NOT Physical
    if( MnvHist::IsNotPhysicalShift(x_shift) || MnvHist::IsNotPhysicalShift(y_shift) ) return;

    const double x_shiftVal = xval + x_shift;
    const double y_shiftVal = yval + y_shift;
    int bin = hist->FindBin( x_shiftVal, y_shiftVal );
  
    // wgt_bckg is universe_wgt / cv_wgt
    //double wgt_bckg = applyBckgConstraints_Unv ? GetBckgConstraint(err_name, unv) : 1.0;

    double wgt_bckg = 1.0; // We fill MnvH2D only for Signal Events: No need to search for Bckg Constraints
    double wgtU = cvweight * wgt_bckg;
    err_hists[unv]->AddBinContent( bin, wgtU );

    const double err = err_hists[unv]->GetBinError(bin);
    const double newerr2 = err*err + wgtU*wgtU;
    const double newerr = (0.<newerr2) ? sqrt(newerr2) : 0.;
    err_hists[unv]->SetBinError( bin, newerr );

}

double CCProtonPi0_Analyzer::GetBckgConstraint(std::string error_name, int hist_ind)
{
    // Find the Bckg Constraint if it is one of the constrained events
    if (truth_isBckg_Compact_SinglePiPlus || truth_isBckg_Compact_QELike || truth_isBckg_Compact_WithPi0){
        int begin = 0;
        int end = 0;
        GetSearchRange(error_name, begin, end); 
        for (int i = begin; i <=end; ++i){
            if (hist_ind == error_hist_inds[i] && error_name.compare(error_names[i]) == 0){
                if (truth_isBckg_Compact_SinglePiPlus) return error_wgt_SinglePiPlus[i] / cv_wgt_SinglePiPlus;
                else if (truth_isBckg_Compact_QELike) return error_wgt_QELike[i] / cv_wgt_QELike;
                else if (truth_isBckg_Compact_WithPi0) return error_wgt_WithPi0[i] / cv_wgt_WithPi0;
                else{
                    std::cout<<"WARNING! You tried to get a Background Constraint for a Non-Constraint Event"<<std::endl;
                    exit(1);
                }
            }
        }
    //  If event is Signal or Bckg_Other weight is 1.0
    }else{
        return 1.0;
    }

    std::cout<<"WARNING! Can not find BckgConstraint for "<<error_name<<" Hist = "<<hist_ind<<std::endl;
    exit(1);
}

void CCProtonPi0_Analyzer::GetSearchRange(std::string err_name, int &begin, int &end)
{
    if (err_name.compare("Flux") == 0){
        begin = 0;
        end = 99;
    }else if (err_name.find("GENIE") != std::string::npos){
        begin = 100;
        end = 169;
    }else if (err_name.compare("MuonTracking") == 0){
        begin = 170;
        end = 171;
    }else if (err_name.compare("EM_EnergyScale") == 0){
        begin = 172;
        end = 671;
    }else if (err_name.compare("MuonMomentum") == 0){
        begin = 672;
        end = 1172;
    }else{
        std::cout<<"WARNING! Can not find Search Rangefor "<<err_name<<std::endl;
        exit(1);
    }
}

void CCProtonPi0_Analyzer::ReadBckgConstraints()
{
    ifstream file;
    file.open(Folder_List::BckgConstraints.c_str());

    if (!file.is_open()){
        std::cout<<"WARNING! Cannot open input file "<<Folder_List::BckgConstraints<<std::endl;
        exit(1);
    }

    // Read Header and Discard (Don't use it)
    std::string line;
    getline(file, line);

    std::string error_name;
    int hist_ind;
    double dummy;
    double wgt_SinglePiPlus;
    double wgt_QELike;
    double wgt_WithPi0;

    // Read CV Result and Save CV Weights
    getline(file, line);
    std::stringstream cv_line_ss(line);
    cv_line_ss>>error_name>>hist_ind>>dummy>>dummy>>wgt_SinglePiPlus>>wgt_QELike>>wgt_WithPi0>>dummy>>dummy>>dummy;
   
    if (error_name.compare("CentralValue") == 0){
        cv_wgt_SinglePiPlus = wgt_SinglePiPlus;
        cv_wgt_QELike = wgt_QELike;
        cv_wgt_WithPi0 = wgt_WithPi0;
    }else{
        std::cout<<"WARNING! Cannot Read CentralValue Bckg Constraints!"<<std::endl;
        exit(1);
    }

    // Read Other Lines
    while(!file.eof()){
        getline(file,line);
        std::stringstream line_ss(line);
        line_ss>>error_name>>hist_ind>>dummy>>dummy>>wgt_SinglePiPlus>>wgt_QELike>>wgt_WithPi0>>dummy>>dummy>>dummy;
      
        error_names.push_back(error_name);
        error_hist_inds.push_back(hist_ind);
        error_wgt_SinglePiPlus.push_back(wgt_SinglePiPlus);
        error_wgt_QELike.push_back(wgt_QELike);
        error_wgt_WithPi0.push_back(wgt_WithPi0);
    }

    file.close();

//    for (unsigned int i = 0; i < error_names.size(); ++i){
//        std::cout<<i<<" ";
//        std::cout<<error_names[i]<<" ";
//        std::cout<<error_hist_inds[i]<<" ";
//        std::cout<<error_wgt_SinglePiPlus[i]<<" ";
//        std::cout<<error_wgt_QELike[i]<<" ";
//        std::cout<<error_wgt_WithPi0[i]<<std::endl;
//    }

}

void CCProtonPi0_Analyzer::FillLatErrorBands_Auto()
{
    // Test Function for Comparing Fill_ByHand
    std::vector<double> pi0_P_random_shifts;
    for (unsigned int i = 0; i < em_random_shifts.size(); ++i ){
        double pi0_P_i = (1.0 + em_random_shifts[i]) * pi0_P;
        double shift = pi0_P_i - pi0_P;
        pi0_P_random_shifts.push_back(shift * MeV_to_GeV);
    }

    // Fill for pi0 -- Shifted Distribution
    pi0.pi0_P_mc_reco_all->FillLatErrorBand("EM_EnergyScale", pi0_P * MeV_to_GeV, pi0_P_random_shifts, cvweight);    
    pi0.pi0_P_response->FillLatErrorBand("EM_EnergyScale", pi0_P * MeV_to_GeV, truth_pi0_P * MeV_to_GeV, pi0_P_random_shifts, no_random_shifts, cvweight);    

    // Fill for muon -- No Shift 
    muon.muon_P_mc_reco_all->FillLatErrorBand("EM_EnergyScale", muon_P * MeV_to_GeV, no_random_shifts, cvweight);    
    muon.muon_P_response->FillLatErrorBand("EM_EnergyScale", muon_P * MeV_to_GeV, truth_muon_P * MeV_to_GeV, no_random_shifts, no_random_shifts, cvweight);    
}

void CCProtonPi0_Analyzer::FillLatErrorBands_invMass_Auto()
{
    // Test Function for Comparing Fill_ByHand
    std::vector<double> pi0_invMass_random_shifts;
    for (unsigned int i = 0; i < em_random_shifts.size(); ++i ){
        double pi0_invMass_i = (1.0 + em_random_shifts[i]) * pi0_invMass;
        double shift = pi0_invMass_i - pi0_invMass;
        pi0_invMass_random_shifts.push_back(shift);
    }

    // Fill for EM_EnergyScale -- Shifted Distribution
    cutList.invMass_mc_reco_all->FillLatErrorBand("EM_EnergyScale", pi0_invMass, pi0_invMass_random_shifts, cvweight);    
    cutList.invMass_mc_reco_signal->FillLatErrorBand("EM_EnergyScale", pi0_invMass, pi0_invMass_random_shifts, cvweight);    
    cutList.invMass_mc_reco_bckg->FillLatErrorBand("EM_EnergyScale", pi0_invMass, pi0_invMass_random_shifts, cvweight);    

    // Fill for Muon Momentum -- No Shift
    cutList.invMass_mc_reco_all->FillLatErrorBand("MuonMomentum", pi0_invMass, no_random_shifts, cvweight);    
    cutList.invMass_mc_reco_signal->FillLatErrorBand("MuonMomentum", pi0_invMass, no_random_shifts, cvweight);    
    cutList.invMass_mc_reco_bckg->FillLatErrorBand("MuonMomentum", pi0_invMass, no_random_shifts, cvweight);    
}


double CCProtonPi0_Analyzer::GetNeutronResponseErr()
{
    int nNeutrons = CountFSParticles(2112, 150);

    if (nNeutrons < 1) return 0.0;

    int ind = GetLeadingNeutronInd();
    if (ind == -1) return 0.0; 

    double neutron_path_length = CalcNeutronPathLength(ind);
    bool isInelastic = isNeutronInelastic(ind); 
    
    bool reweightable = (neutron_path_length > 0.0) && isInelastic; 
   
    double coeff = 0.0;
    if (reweightable) {

        double neutron_P = HEP_Functions::calcMomentum(detmc_traj_px0[ind], detmc_traj_py0[ind], detmc_traj_pz0[ind]);

        // From Kevin, taken from coherent pion analysis
        const double c0 = 0.0010528;
        const double c1 = -0.3269;
        const double c2 = 159.55;
        double l0 = 1/(c0 + c1/neutron_P + c2/(neutron_P * neutron_P));
        double l  = neutron_path_length;
        bool isContained = isPointContained(detmc_traj_xf[ind], detmc_traj_yf[ind], detmc_traj_zf[ind]);
        if (isContained) coeff = (-l/l0) / (exp(l/l0) - 1);
        else coeff = l/l0;
    }

    double neutron_KE = detmc_traj_E0[ind] - neutron_mass; 
    double f = Calc_f(neutron_KE);

    return coeff * f;
}

double CCProtonPi0_Analyzer::Calc_f(double neutron_KE)
{
    double var_unc = std::max(0.0,std::min(150.0,-150 + neutron_KE)/1500.)
        + std::min(0.25,0.1 + (3*std::max(0.0,50 - neutron_KE))/500.);

    return var_unc;
}

int CCProtonPi0_Analyzer::GetLeadingNeutronInd()
{
    double Emax = -1.0;
    int ind = -1;

    for (int i = 0; i < detmc_ntrajectory2; ++i){
            int parent = detmc_traj_mother[i];
            if (parent > 0) continue;
            
            int pdg = detmc_traj_pdg[i];
            if (pdg != 2112) continue;

            double E = detmc_traj_E0[i];
            
            if (E > Emax) ind = i;
    }
    
    return ind;
}

// contained in the detector (apothem = 1100, zmin = 5430.0, zmax = 10000.0)
// x,y,z in mm
bool CCProtonPi0_Analyzer::isPointContained(double x, double y, double z)
{
    if ((5430.0 < z && z < 10000.0) && inside_hexagon(x,y,1100.0)) return true;

    return false;
}

// inside hexagon test
// x,y, apothem should have the same units
bool CCProtonPi0_Analyzer::inside_hexagon(double x, double y, double apothem)
{
    double max_bound =  apothem;
    double min_bound = -apothem;

    double delta_phi = 60.0 * TMath::DegToRad();
    double phi = 0.0;
    for (int i = 0; i < 3; ++i) {
        double xtmp = cos(phi)*x + sin(phi)*y;
        phi += delta_phi;

        if (xtmp < min_bound || xtmp > max_bound) return false;
    }

    return true;

}

void CCProtonPi0_Analyzer::PrintLeadingNeutron(int i)
{
    double neutron_P = HEP_Functions::calcMomentum(detmc_traj_px0[i], detmc_traj_py0[i], detmc_traj_pz0[i]);
    double neutron_KE = detmc_traj_E0[i] - neutron_mass; 
    double neutron_path_length = CalcNeutronPathLength(i);

    std::cout<<"neutron_PDG = "<<detmc_traj_pdg[i]<<" "<<"Mother = "<<detmc_traj_mother[i]<<std::endl;
    std::cout<<"neutron_P = "<<neutron_P<<std::endl;
    std::cout<<"neutron_KE = "<<neutron_KE<<std::endl;
    std::cout<<"neutron_path_length = "<<neutron_path_length<<std::endl;
    std::cout<<"--------------------"<<std::endl;
}

double CCProtonPi0_Analyzer::CalcNeutronPathLength(int i)
{
    double path_length = 0.0;
    bool isContained = isPointContained(detmc_traj_xf[i], detmc_traj_yf[i], detmc_traj_zf[i]);

    if (isContained){
        path_length = HEP_Functions::calcDistance(detmc_traj_x0[i], detmc_traj_y0[i], detmc_traj_z0[i], detmc_traj_xf[i], detmc_traj_yf[i], detmc_traj_zf[i]); 
    }else{
            // calculate contained length here
            // step along neutron initial direction (should be from the last contained
            // point, but we don't have that info for now) in step of 5cm
        double P = HEP_Functions::calcMomentum(detmc_traj_px0[i], detmc_traj_py0[i], detmc_traj_pz0[i]);
        double sx = detmc_traj_px0[i]/P;
        double sy = detmc_traj_py0[i]/P;
        double sz = detmc_traj_pz0[i]/P;

        double length = 0.0;
        const double step_size = 50.0; 
        for (;; length += step_size) {
            double new_x = detmc_traj_px0[i] + length * sx;
            double new_y = detmc_traj_py0[i] + length * sy;
            double new_z = detmc_traj_pz0[i] + length * sz;

            if (!isPointContained(new_x,new_y,new_z)) break;
        }
    
        path_length = length;
    }

    return path_length;
}

bool CCProtonPi0_Analyzer::isNeutronInelastic(int ind)
{
    const char* processes[] = {"hadElastic","NeutronInelastic",
        "ProtonInelastic","PionMinusInelastic",
        "PionPlusInelastic","CHIPS",
        "Decay",
        "conv","compt"};

    bool verbose = false;
    if (verbose) {
        std::cout<<"\tthe neutron is: "<<setw(6)<<detmc_traj_id[ind];
        std::cout<<setw(15)<<detmc_traj_xf[ind];
        std::cout<<setw(15)<<detmc_traj_yf[ind];
        std::cout<<setw(15)<<detmc_traj_zf[ind];
        std::cout<<std::endl;
    }

    for (int i = 0; i < detmc_ntrajectory2; ++i) {
        int parent = detmc_traj_mother[i];

        if (parent == 0) continue;

        if (parent != detmc_traj_id[ind]) continue;

        int pdg  = detmc_traj_pdg[i];
        int proc = detmc_traj_proc[i];

        double x = detmc_traj_x0[i];
        double y = detmc_traj_y0[i];
        double z = detmc_traj_z0[i];

        if (verbose) {
            std::cout<<"\tdaughter: "<<"parent = "<<setw(6)<<parent;
            std::cout<<" pdg = "<<setw(8)<<pdg;
            std::cout<<setw(15)<<x<<setw(15)<<y<<setw(15)<<z<< std::endl;
        }

        double d = HEP_Functions::calcDistance(x,y,z, detmc_traj_xf[ind], detmc_traj_yf[ind], detmc_traj_zf[ind]);

        if (d < 1.e-3) {
            if (verbose) {
                std::cout<<"\tend point daughter: " <<setw(6)<<"parent = "<<parent;
                std::cout<<setw(15)<<x<<setw(15)<<y<<setw(15)<<z<<setw(20)<<processes[proc]<<" "<<proc<<std::endl;
            }

            if (proc == 1) return true;
        }
    }

    return false;
}

#endif //CCProtonPi0_Analyzer_Systematics_cpp
