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
    AddVertErrorBands_Data(interaction.QSq_MaRES[0]);

    AddLatErrorBands_Data(pi0.pi0_P_all);
    AddLatErrorBands_Data(pi0.pi0_KE_all);
    AddLatErrorBands_Data(pi0.pi0_theta_all);
    AddLatErrorBands_Data(muon.muon_P_all);
    AddLatErrorBands_Data(muon.muon_theta_all);
    AddLatErrorBands_Data(interaction.QSq_all);
    AddLatErrorBands_Data(interaction.Enu_all);
    AddLatErrorBands_Data(interaction.W_all);

    // Delta(1232) enriched Sample
    AddLeadingErrorBands_Data(interaction.deltaInvMass_all);
    AddLeadingErrorBands_Data(interaction.Delta_pi_theta_all);
    AddLeadingErrorBands_Data(interaction.Delta_pi_phi_all);
    
    // 2p2h Search
    AddLeadingErrorBands_Data(interaction.muon_theta_muon_KE_all);
    AddLeadingErrorBands_Data(interaction.q3_q0_all);

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
    //h->FillVertErrorBand("GENIE_NormCCQE"          ,xval, yval, truth_genie_wgt_NormCCQE[2]         , truth_genie_wgt_NormCCQE[4]         , cvweight);
    h->FillVertErrorBand("GENIE_NormDISCC"         ,xval, yval, truth_genie_wgt_NormDISCC[2]        , truth_genie_wgt_NormDISCC[4]        , cvweight);
    h->FillVertErrorBand("GENIE_NormNCRES"         ,xval, yval, truth_genie_wgt_NormNCRES[2]        , truth_genie_wgt_NormNCRES[4]        , cvweight);
    h->FillVertErrorBand("GENIE_RDecBR1gamma"      ,xval, yval, truth_genie_wgt_RDecBR1gamma[2]     , truth_genie_wgt_RDecBR1gamma[4]     , cvweight);
    h->FillVertErrorBand("GENIE_Rvn2pi"            ,xval, yval, truth_genie_wgt_Rvn2pi[2]           , truth_genie_wgt_Rvn2pi[4]           , cvweight);
    h->FillVertErrorBand("GENIE_Rvp2pi"            ,xval, yval, truth_genie_wgt_Rvp2pi[2]           , truth_genie_wgt_Rvp2pi[4]           , cvweight);
    h->FillVertErrorBand("GENIE_VecFFCCQEshape"    ,xval, yval, truth_genie_wgt_VecFFCCQEshape[2]   , truth_genie_wgt_VecFFCCQEshape[4]   , cvweight);

    h->FillVertErrorBand("GENIE_NormCCRES"         ,xval, yval, updated_genie_wgt_NormCCRES[2]        , updated_genie_wgt_NormCCRES[4]        , cvweight);
    h->FillVertErrorBand("GENIE_Theta_Delta2Npi"   ,xval, yval, updated_genie_wgt_Theta_Delta2Npi[2]  , updated_genie_wgt_Theta_Delta2Npi[4]  , cvweight);
    h->FillVertErrorBand("GENIE_MaRES"             ,xval, yval, updated_genie_wgt_MaRES[2]            , updated_genie_wgt_MaRES[4]            , cvweight);
    h->FillVertErrorBand("GENIE_MvRES"             ,xval, yval, updated_genie_wgt_MvRES[2]            , updated_genie_wgt_MvRES[4]            , cvweight);
    h->FillVertErrorBand("GENIE_Rvn1pi"            ,xval, yval, updated_genie_wgt_Rvn1pi[2]           , updated_genie_wgt_Rvn1pi[4]           , cvweight);
    h->FillVertErrorBand("GENIE_Rvp1pi"            ,xval, yval, updated_genie_wgt_Rvp1pi[2]           , updated_genie_wgt_Rvp1pi[4]           , cvweight);
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
    //FillVertErrorBand_ByHand(h, xval, yval, "GENIE_NormCCQE"          , truth_genie_wgt_NormCCQE[2]         , truth_genie_wgt_NormCCQE[4]         );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_NormDISCC"         , truth_genie_wgt_NormDISCC[2]        , truth_genie_wgt_NormDISCC[4]        );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_NormNCRES"         , truth_genie_wgt_NormNCRES[2]        , truth_genie_wgt_NormNCRES[4]        );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_RDecBR1gamma"      , truth_genie_wgt_RDecBR1gamma[2]     , truth_genie_wgt_RDecBR1gamma[4]     );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_Rvn2pi"            , truth_genie_wgt_Rvn2pi[2]           , truth_genie_wgt_Rvn2pi[4]           );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_Rvp2pi"            , truth_genie_wgt_Rvp2pi[2]           , truth_genie_wgt_Rvp2pi[4]           );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_VecFFCCQEshape"    , truth_genie_wgt_VecFFCCQEshape[2]   , truth_genie_wgt_VecFFCCQEshape[4]   );

    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_NormCCRES"         , updated_genie_wgt_NormCCRES[2]        , updated_genie_wgt_NormCCRES[4]        );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_Theta_Delta2Npi"   , updated_genie_wgt_Theta_Delta2Npi[2]  , updated_genie_wgt_Theta_Delta2Npi[4]  );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_MaRES"             , updated_genie_wgt_MaRES[2]            , updated_genie_wgt_MaRES[4]            );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_MvRES"             , updated_genie_wgt_MvRES[2]            , updated_genie_wgt_MvRES[4]            );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_Rvn1pi"            , updated_genie_wgt_Rvn1pi[2]           , updated_genie_wgt_Rvn1pi[4]           );
    FillVertErrorBand_ByHand(h, xval, yval, "GENIE_Rvp1pi"            , updated_genie_wgt_Rvp1pi[2]           , updated_genie_wgt_Rvp1pi[4]           );
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
        double wgt_bckg = applyBckgConstraints_Unv ? GetBckgConstraint(error_name, i) : 1.0;
        const double applyWeight = cvweight * wgt_bckg;
        const double wgtU = errors[i]*applyWeight;
        err_hists[i]->AddBinContent( cvbin, wgtU );

        //const double err = err_hists[i]->GetBinError(cvbin);
        //const double newerr2 = err*err + wgtU*wgtU; 
        //const double newerr = (0.<newerr2) ? sqrt(newerr2) : 0.;
        //err_hists[i]->SetBinError( cvbin, newerr );
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
    double correctionErr = GetPionResponseErr();
    h->FillVertErrorBand("PionResponse", var, 1-correctionErr, 1+correctionErr, cvweight);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_PionResponse_ByHand(MnvH1D* h, double var)
{
    double correctionErr = GetPionResponseErr();
    FillVertErrorBand_ByHand(h, var, "PionResponse", 1-correctionErr, 1+correctionErr);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_PionResponse(MnvH2D* h, double xval, double yval)
{
    double correctionErr = GetPionResponseErr();
    h->FillVertErrorBand("PionResponse", xval, yval, 1-correctionErr, 1+correctionErr, cvweight);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_PionResponse_ByHand(MnvH2D* h, double xval, double yval)
{
    double correctionErr = GetPionResponseErr();
    FillVertErrorBand_ByHand(h, xval, yval, "PionResponse", 1-correctionErr, 1+correctionErr);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_BckgConstraint_WithPi0(MnvH1D* h, double var)
{
    double correctionErr = GetBckgConstraint_WithPi0Err();
    h->FillVertErrorBand("BckgConstraint_WithPi0", var, 1-correctionErr, 1+correctionErr, cvweight);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_BckgConstraint_WithPi0_ByHand(MnvH1D* h, double var)
{
    double correctionErr = GetBckgConstraint_WithPi0Err();
    FillVertErrorBand_ByHand(h, var, "BckgConstraint_WithPi0", 1-correctionErr, 1+correctionErr);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_BckgConstraint_WithPi0(MnvH2D* h, double xval, double yval)
{
    double correctionErr = GetBckgConstraint_WithPi0Err();
    h->FillVertErrorBand("BckgConstraint_WithPi0", xval, yval, 1-correctionErr, 1+correctionErr, cvweight);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_BckgConstraint_WithPi0_ByHand(MnvH2D* h, double xval, double yval)
{
    double correctionErr = GetBckgConstraint_WithPi0Err();
    FillVertErrorBand_ByHand(h, xval, yval, "BckgConstraint_WithPi0", 1-correctionErr, 1+correctionErr);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_BckgConstraint_SinglePiPlus(MnvH1D* h, double var)
{
    double correctionErr = GetBckgConstraint_SinglePiPlusErr();
    h->FillVertErrorBand("BckgConstraint_SinglePiPlus", var, 1-correctionErr, 1+correctionErr, cvweight);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_BckgConstraint_SinglePiPlus_ByHand(MnvH1D* h, double var)
{
    double correctionErr = GetBckgConstraint_SinglePiPlusErr();
    FillVertErrorBand_ByHand(h, var, "BckgConstraint_SinglePiPlus", 1-correctionErr, 1+correctionErr);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_BckgConstraint_SinglePiPlus(MnvH2D* h, double xval, double yval)
{
    double correctionErr = GetBckgConstraint_SinglePiPlusErr();
    h->FillVertErrorBand("BckgConstraint_SinglePiPlus", xval, yval, 1-correctionErr, 1+correctionErr, cvweight);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_BckgConstraint_SinglePiPlus_ByHand(MnvH2D* h, double xval, double yval)
{
    double correctionErr = GetBckgConstraint_SinglePiPlusErr();
    FillVertErrorBand_ByHand(h, xval, yval, "BckgConstraint_SinglePiPlus", 1-correctionErr, 1+correctionErr);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_BckgConstraint_QELike(MnvH1D* h, double var)
{
    double correctionErr = GetBckgConstraint_QELikeErr();
    h->FillVertErrorBand("BckgConstraint_QELike", var, 1-correctionErr, 1+correctionErr, cvweight);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_BckgConstraint_QELike_ByHand(MnvH1D* h, double var)
{
    double correctionErr = GetBckgConstraint_QELikeErr();
    FillVertErrorBand_ByHand(h, var, "BckgConstraint_QELike", 1-correctionErr, 1+correctionErr);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_BckgConstraint_QELike(MnvH2D* h, double xval, double yval)
{
    double correctionErr = GetBckgConstraint_QELikeErr();
    h->FillVertErrorBand("BckgConstraint_QELike", xval, yval, 1-correctionErr, 1+correctionErr, cvweight);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_BckgConstraint_QELike_ByHand(MnvH2D* h, double xval, double yval)
{
    double correctionErr = GetBckgConstraint_QELikeErr();
    FillVertErrorBand_ByHand(h, xval, yval, "BckgConstraint_QELike", 1-correctionErr, 1+correctionErr);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_MichelTrue(MnvH1D* h, double var)
{
    double correctionErr = GetMichelTrueErr();
    h->FillVertErrorBand("MichelTrue", var, 1-correctionErr, 1+correctionErr, cvweight);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_MichelTrue_ByHand(MnvH1D* h, double var)
{
    double correctionErr = GetMichelTrueErr();
    FillVertErrorBand_ByHand(h, var, "MichelTrue", 1-correctionErr, 1+correctionErr);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_MichelTrue(MnvH2D* h, double xval, double yval)
{
    double correctionErr = GetMichelTrueErr();
    h->FillVertErrorBand("MichelTrue", xval, yval, 1-correctionErr, 1+correctionErr, cvweight);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_MichelTrue_ByHand(MnvH2D* h, double xval, double yval)
{
    double correctionErr = GetMichelTrueErr();
    FillVertErrorBand_ByHand(h, xval, yval, "MichelTrue", 1-correctionErr, 1+correctionErr);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_MichelFake(MnvH1D* h, double var)
{
    double correctionErr = GetMichelFakeErr();
    h->FillVertErrorBand("MichelFake", var, 1-correctionErr, 1+correctionErr, cvweight);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_MichelFake_ByHand(MnvH1D* h, double var)
{
    double correctionErr = GetMichelFakeErr();
    FillVertErrorBand_ByHand(h, var, "MichelFake", 1-correctionErr, 1+correctionErr);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_MichelFake(MnvH2D* h, double xval, double yval)
{
    double correctionErr = GetMichelFakeErr();
    h->FillVertErrorBand("MichelFake", xval, yval, 1-correctionErr, 1+correctionErr, cvweight);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_MichelFake_ByHand(MnvH2D* h, double xval, double yval)
{
    double correctionErr = GetMichelFakeErr();
    FillVertErrorBand_ByHand(h, xval, yval, "MichelFake", 1-correctionErr, 1+correctionErr);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_TargetMass(MnvH1D* h, double var)
{
    double correctionErr = GetTargetMassErr();
    h->FillVertErrorBand("TargetMass", var, 1-correctionErr, 1+correctionErr, cvweight);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_TargetMass_ByHand(MnvH1D* h, double var)
{
    double correctionErr = GetTargetMassErr();
    FillVertErrorBand_ByHand(h, var, "TargetMass", 1-correctionErr, 1+correctionErr);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_TargetMass(MnvH2D* h, double xval, double yval)
{
    double correctionErr = GetTargetMassErr();
    h->FillVertErrorBand("TargetMass", xval, yval, 1-correctionErr, 1+correctionErr, cvweight);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_TargetMass_ByHand(MnvH2D* h, double xval, double yval)
{
    double correctionErr = GetTargetMassErr();
    FillVertErrorBand_ByHand(h, xval, yval, "TargetMass", 1-correctionErr, 1+correctionErr);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_2p2h(MnvH1D* h, double var)
{
    double correctionErr = err_2p2h;
    h->FillVertErrorBand("2p2h", var, 1-correctionErr, 1+correctionErr, cvweight);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_2p2h_ByHand(MnvH1D* h, double var)
{
    double correctionErr = err_2p2h;
    FillVertErrorBand_ByHand(h, var, "2p2h", 1-correctionErr, 1+correctionErr);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_2p2h(MnvH2D* h, double xval, double yval)
{
    double correctionErr = err_2p2h;
    h->FillVertErrorBand("2p2h", xval, yval, 1-correctionErr, 1+correctionErr, cvweight);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_2p2h_ByHand(MnvH2D* h, double xval, double yval)
{
    double correctionErr = err_2p2h;
    FillVertErrorBand_ByHand(h, xval, yval, "2p2h", 1-correctionErr, 1+correctionErr);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_Unfolding(MnvH1D* h, double var)
{
    double correctionErr = GetUnfoldingErr();
    h->FillVertErrorBand("Unfolding", var, 1-correctionErr, 1+correctionErr, cvweight);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_Unfolding_ByHand(MnvH1D* h, double var)
{
    double correctionErr = GetUnfoldingErr();
    FillVertErrorBand_ByHand(h, var, "Unfolding", 1-correctionErr, 1+correctionErr);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_Unfolding(MnvH2D* h, double xval, double yval)
{
    double correctionErr = GetUnfoldingErr();
    h->FillVertErrorBand("Unfolding", xval, yval, 1-correctionErr, 1+correctionErr, cvweight);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_Unfolding_ByHand(MnvH2D* h, double xval, double yval)
{
    double correctionErr = GetUnfoldingErr();
    FillVertErrorBand_ByHand(h, xval, yval, "Unfolding", 1-correctionErr, 1+correctionErr);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_ProtonTracking(MnvH1D* h, double var)
{
    double correctionErr = GetProtonTrackingErr();
    h->FillVertErrorBand("ProtonTracking", var, 1-correctionErr, 1+correctionErr, cvweight);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_ProtonTracking_ByHand(MnvH1D* h, double var)
{
    double correctionErr = GetProtonTrackingErr();
    FillVertErrorBand_ByHand(h, var, "ProtonTracking", 1-correctionErr, 1+correctionErr);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_ProtonTracking(MnvH2D* h, double xval, double yval)
{
    double correctionErr = GetProtonTrackingErr();
    h->FillVertErrorBand("ProtonTracking", xval, yval, 1-correctionErr, 1+correctionErr, cvweight);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_ProtonTracking_ByHand(MnvH2D* h, double xval, double yval)
{
    double correctionErr = GetProtonTrackingErr();
    FillVertErrorBand_ByHand(h, xval, yval, "ProtonTracking", 1-correctionErr, 1+correctionErr);
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
    std::vector<double> flux_errors = GetFluxError(mc_incomingE * MeV_to_GeV, mc_incoming);
    h->FillVertErrorBand("Flux",  var, &flux_errors[0],  cvweight, 1.0);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_Flux_ByHand(MnvH1D* h, double var)
{
    std::vector<double> flux_errors = GetFluxError(mc_incomingE * MeV_to_GeV, mc_incoming);
    FillVertErrorBand_ByHand(h, var, "Flux", flux_errors);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_Flux(MnvH2D* h, double xval, double yval)
{
    std::vector<double> flux_errors = GetFluxError(mc_incomingE * MeV_to_GeV, mc_incoming);
    h->FillVertErrorBand("Flux",  xval, yval,  &flux_errors[0],  cvweight, 1.0);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_Flux_ByHand(MnvH2D* h, double xval, double yval)
{
    std::vector<double> flux_errors = GetFluxError(mc_incomingE * MeV_to_GeV, mc_incoming);
    FillVertErrorBand_ByHand(h, xval, yval, "Flux", flux_errors);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_HighMaRES_ByHand(MnvH1D* h, double var)
{
    std::vector<double> wgts_up = QSqFitter.GetWeights(truth_genie_wgt_MaRES[4], truth_genie_wgt_MaRES[5]);
    FillVertErrorBand_ByHand(h, var, "HighMaRES", wgts_up);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_HighMaRES_ByHand(MnvH2D* h, double xval, double yval)
{
    std::vector<double> wgts_up = QSqFitter.GetWeights(truth_genie_wgt_MaRES[4], truth_genie_wgt_MaRES[5]);
    FillVertErrorBand_ByHand(h, xval, yval, "HighMaRES", wgts_up);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_LowMaRES_ByHand(MnvH1D* h, double var)
{
    std::vector<double> wgts_dn = QSqFitter.GetWeights(truth_genie_wgt_MaRES[2], truth_genie_wgt_MaRES[1]);
    FillVertErrorBand_ByHand(h, var, "LowMaRES", wgts_dn);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_LowMaRES_ByHand(MnvH2D* h, double xval, double yval)
{
    std::vector<double> wgts_dn = QSqFitter.GetWeights(truth_genie_wgt_MaRES[2], truth_genie_wgt_MaRES[1]);
    FillVertErrorBand_ByHand(h, xval, yval, "LowMaRES", wgts_dn);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_DeltaFactor_ByHand(MnvH1D* h, double var)
{
    std::vector<double> wgts = GetDeltaFactorWeights();
    FillVertErrorBand_ByHand(h, var, "DeltaFactor", wgts);
}

void CCProtonPi0_Analyzer::FillVertErrorBand_DeltaFactor_ByHand(MnvH2D* h, double xval, double yval)
{
    std::vector<double> wgts = GetDeltaFactorWeights();
    FillVertErrorBand_ByHand(h, xval, yval, "DeltaFactor", wgts);
}

std::vector<double> CCProtonPi0_Analyzer::GetDeltaFactorWeights()
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

void CCProtonPi0_Analyzer::FillHistogramWithVertErrors(MnvH1D* hist, double var)
{
    // Fill CV Value
    hist->Fill(var, cvweight);

    // Fill Vertical Error Bands
    if (fillErrors_ByHand){
        FillVertErrorBand_Genie_ByHand(hist, var);
        FillVertErrorBand_Flux_ByHand(hist, var);
        FillVertErrorBand_BckgConstraint_WithPi0_ByHand(hist, var);
        FillVertErrorBand_BckgConstraint_SinglePiPlus_ByHand(hist, var);
        FillVertErrorBand_BckgConstraint_QELike_ByHand(hist, var);
        FillVertErrorBand_MichelTrue_ByHand(hist, var);
        FillVertErrorBand_MichelFake_ByHand(hist, var);
        FillVertErrorBand_TargetMass_ByHand(hist, var);
        FillVertErrorBand_Unfolding_ByHand(hist, var);
        FillVertErrorBand_2p2h_ByHand(hist, var);
        FillVertErrorBand_ProtonTracking_ByHand(hist, var);
        FillVertErrorBand_MuonTracking_ByHand(hist, var);
        FillVertErrorBand_PionResponse_ByHand(hist, var);
        FillVertErrorBand_NeutronResponse_ByHand(hist, var);
    
        // QSq Study
        //FillVertErrorBand_HighMaRES_ByHand(hist, var);
        //FillVertErrorBand_LowMaRES_ByHand(hist, var);
        //FillVertErrorBand_DeltaFactor_ByHand(hist, var);
    }else{
        FillVertErrorBand_Genie(hist, var);
        FillVertErrorBand_Flux(hist, var);
        FillVertErrorBand_BckgConstraint_WithPi0(hist, var);
        FillVertErrorBand_BckgConstraint_SinglePiPlus(hist, var);
        FillVertErrorBand_BckgConstraint_QELike(hist, var);
        FillVertErrorBand_MichelTrue(hist, var);
        FillVertErrorBand_MichelFake(hist, var);
        FillVertErrorBand_TargetMass(hist, var);
        FillVertErrorBand_Unfolding(hist, var);
        FillVertErrorBand_2p2h(hist, var);
        FillVertErrorBand_ProtonTracking(hist, var);
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
        FillVertErrorBand_BckgConstraint_WithPi0_ByHand(hist, xval, yval);
        FillVertErrorBand_BckgConstraint_SinglePiPlus_ByHand(hist, xval, yval);
        FillVertErrorBand_BckgConstraint_QELike_ByHand(hist, xval, yval);
        FillVertErrorBand_MichelTrue_ByHand(hist, xval, yval);
        FillVertErrorBand_MichelFake_ByHand(hist, xval, yval);
        FillVertErrorBand_TargetMass_ByHand(hist, xval, yval);
        FillVertErrorBand_Unfolding_ByHand(hist, xval, yval);
        FillVertErrorBand_2p2h_ByHand(hist, xval, yval);
        FillVertErrorBand_ProtonTracking_ByHand(hist, xval, yval);
        FillVertErrorBand_MuonTracking_ByHand(hist, xval, yval);
        FillVertErrorBand_PionResponse_ByHand(hist, xval, yval);
        FillVertErrorBand_NeutronResponse_ByHand(hist, xval, yval);

        // QSq Study
        //FillVertErrorBand_HighMaRES_ByHand(hist, xval, yval);
        //FillVertErrorBand_LowMaRES_ByHand(hist, xval, yval);
        //FillVertErrorBand_DeltaFactor_ByHand(hist, xval, yval);
    }else{
        FillVertErrorBand_Genie(hist, xval, yval);
        FillVertErrorBand_Flux(hist, xval, yval);
        FillVertErrorBand_BckgConstraint_WithPi0(hist, xval, yval);
        FillVertErrorBand_BckgConstraint_SinglePiPlus(hist, xval, yval);
        FillVertErrorBand_BckgConstraint_QELike(hist, xval, yval);
        FillVertErrorBand_MichelTrue(hist, xval, yval);
        FillVertErrorBand_MichelFake(hist, xval, yval);
        FillVertErrorBand_TargetMass(hist, xval, yval);
        FillVertErrorBand_Unfolding(hist, xval, yval);
        FillVertErrorBand_2p2h(hist, xval, yval);
        FillVertErrorBand_ProtonTracking(hist, xval, yval);
        FillVertErrorBand_MuonTracking(hist, xval, yval);
        FillVertErrorBand_PionResponse(hist, xval, yval);
        FillVertErrorBand_NeutronResponse(hist, xval, yval);
    }
}

void CCProtonPi0_Analyzer::FillHistogramWithLeadingErrors(MnvH1D* hist, double var)
{
    // Fill CV Value
    hist->Fill(var, cvweight);
    
    FillVertErrorBand_Genie_ByHand(hist, var);
    FillVertErrorBand_Flux_ByHand(hist, var);
    FillVertErrorBand_2p2h_ByHand(hist, var);
}

void CCProtonPi0_Analyzer::FillHistogramWithLeadingErrors(MnvH2D* hist, double xval, double yval)
{
    // Fill CV Value
    hist->Fill(xval, yval, cvweight);

    FillVertErrorBand_Genie_ByHand(hist, xval, yval);
    FillVertErrorBand_Flux_ByHand(hist, xval, yval);
    FillVertErrorBand_2p2h_ByHand(hist, xval, yval);
}

void CCProtonPi0_Analyzer::initLateralErrorBandShifts(bool isModeReduce)
{
    Calc_no_random_shifts();
    Calc_em_energy_random_shifts();
    Calc_muon_theta_random_shifts();
 
    // Fill Random Shift Histograms if it is Analysis Mode
    if (isModeReduce){
        ismuonP_shifts_filled = true;
        isBirks_shifts_filled = true;
    }else{
        Fill_RandomShiftHistograms();
        ismuonP_shifts_filled = false;
        isBirks_shifts_filled = false;
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

void CCProtonPi0_Analyzer::Calc_Birks_random_shifts()
{
    if (nProtonCandidates == 0) return;

    Birks_random_shifts2D.clear();

    for (int i = 0; i < nProtonCandidates; ++i){
        double Birks_uncertainty = all_protons_energy_shift_Birks[i]/all_protons_E[i];
        std::vector<double> Birks_random_shifts = RandNumGenerator.GetRandomShifts(Birks_uncertainty); // ~ Gaussian(0.0, Birks_uncertainty)
        // Fill Birks shifts only for 1 event
        if (!isBirks_shifts_filled && Birks_uncertainty != 0.0){
            for (unsigned int i = 0; i < Birks_random_shifts.size(); ++i){
                interaction.Birks_shift_rand_numbers->Fill(Birks_random_shifts[i]);
            }
            isBirks_shifts_filled = true;
        }
        
        Birks_random_shifts2D.push_back(Birks_random_shifts);
    }
}

void CCProtonPi0_Analyzer::FillLatErrorBands_ByHand()
{
    FillLatErrorBand_EM_EnergyScale();
    FillLatErrorBand_MuonMomentum();
    FillLatErrorBand_MuonTheta();
    FillLatErrorBand_ProtonEnergy("ProtonEnergy_MassModel");
    FillLatErrorBand_ProtonEnergy("ProtonEnergy_MEU");
    FillLatErrorBand_ProtonEnergy("ProtonEnergy_BetheBloch");
    FillLatErrorBand_ProtonEnergy_Birks();
}

void CCProtonPi0_Analyzer::FillLatErrorBand_ProtonEnergy(std::string err_name)
{
    double reco_muon_theta = GetCorrectedMuonTheta();

    bool PassedCuts;
    double Enu_i;
    double QSq_i;
    double WSq_i;
    double W_i;

    for (int i = 0; i < 2; ++i){

        if (nProtonCandidates > 0){
            std::vector<double> proton_energy_shifts = GetProtonEnergyShifts(err_name, i);

            // Debugging
            //std::cout<<err_name<<" Universe = "<<i<<std::endl;
            //for (unsigned int j = 0; j < proton_energy_shifts.size(); ++j){
            //    std::cout<<"\t"<<proton_energy_shifts[j]<<std::endl;;
            //}

            // Calculate Event Kinematics with Shifted Proton Energy 
            double total_proton_KE_i = 0.0;
            for (int p = 0; p < nProtonCandidates; ++p ){
                total_proton_KE_i += all_protons_KE[p] + proton_energy_shifts[p]; 
            }
            Enu_i = Calc_Enu_shifted(muon_E, pi0_E, total_proton_KE_i); // Use actual muon energy and pi0 energy but shifted proton kinetic energy
            QSq_i = Calc_QSq(Enu_i, muon_E, muon_P, reco_muon_theta); // Use actual muon energy
            WSq_i = Calc_WSq(Enu_i, QSq_i, muon_E); // Use actual muon energy
            W_i = (WSq_i > 0) ? sqrt(WSq_i) : -1; 

            PassedCuts = IsEnuInRange(Enu_i) && IsWInRange(W_i);

        }else{
            // For No Proton Events, event passes all cuts
            PassedCuts = true;
            Enu_i = m_Enu;
            QSq_i = m_QSq;
            W_i = m_W;
        }

        if (PassedCuts){
            double Enu_shift = (Enu_i - m_Enu) * MeV_to_GeV;
            double QSq_shift = (QSq_i - m_QSq) * MeVSq_to_GeVSq;
            double W_shift = (W_i - m_W) * MeV_to_GeV;

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
            FillLatErrorBand_SingleUniverse(muon.muon_theta_mc_reco_all, err_name, i, reco_muon_theta * TMath::RadToDeg(), 0.0);
            if (truth_isSignal){
                FillLatErrorBand_SingleUniverse(muon.muon_theta_mc_truth_signal, err_name, i, truth_muon_theta * TMath::RadToDeg(), 0.0);
                FillLatErrorBand_SingleUniverse(muon.muon_theta_mc_reco_signal, err_name, i, reco_muon_theta * TMath::RadToDeg(), 0.0);
                FillLatErrorBand_SingleUniverse(muon.muon_theta_response, err_name, i, reco_muon_theta * TMath::RadToDeg(), truth_muon_theta * TMath::RadToDeg(), 0.0, 0.0);
            }else{
                FillLatErrorBand_SingleUniverse(muon.muon_theta_mc_reco_bckg, err_name, i, reco_muon_theta * TMath::RadToDeg(), 0.0);
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
            // Did not passed cuts
        }
    }
}

void CCProtonPi0_Analyzer::FillLatErrorBand_ProtonEnergy_Birks()
{
    std::string err_name = "ProtonEnergy_Birks";
    double reco_muon_theta = GetCorrectedMuonTheta();
  
    bool PassedCuts;
    double Enu_i;
    double QSq_i;
    double WSq_i;
    double W_i;

    // For each event Birks_random_shifts is different
    Calc_Birks_random_shifts();

    for (int i = 0; i < n_lateral_universes; ++i){

        if (nProtonCandidates > 0){
            // Calculate Event Kinematics with Shifted Proton Energy 
            double total_proton_KE_i = 0.0;
            for (int p = 0; p < nProtonCandidates; ++p ){
                total_proton_KE_i += (1.0 + Birks_random_shifts2D[p][i]) * all_protons_KE[p]; 
            }
            Enu_i = Calc_Enu_shifted(muon_E, pi0_E, total_proton_KE_i); // Use actual muon energy and pi0 energy but shifted proton kinetic energy
            QSq_i = Calc_QSq(Enu_i, muon_E, muon_P, reco_muon_theta); // Use actual muon energy
            WSq_i = Calc_WSq(Enu_i, QSq_i, muon_E); // Use actual muon energy
            W_i = (WSq_i > 0) ? sqrt(WSq_i) : -1; 

            PassedCuts = IsEnuInRange(Enu_i) && IsWInRange(W_i);
        }else{
            // For No Proton Events, event passes all cuts
            PassedCuts = true;
            Enu_i = m_Enu;
            QSq_i = m_QSq;
            W_i = m_W;
        }

        if (PassedCuts){
            double Enu_shift = (Enu_i - m_Enu) * MeV_to_GeV;
            double QSq_shift = (QSq_i - m_QSq) * MeVSq_to_GeVSq;
            double W_shift = (W_i - m_W) * MeV_to_GeV;

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
            FillLatErrorBand_SingleUniverse(muon.muon_theta_mc_reco_all, err_name, i, reco_muon_theta * TMath::RadToDeg(), 0.0);
            if (truth_isSignal){
                FillLatErrorBand_SingleUniverse(muon.muon_theta_mc_truth_signal, err_name, i, truth_muon_theta * TMath::RadToDeg(), 0.0);
                FillLatErrorBand_SingleUniverse(muon.muon_theta_mc_reco_signal, err_name, i, reco_muon_theta * TMath::RadToDeg(), 0.0);
                FillLatErrorBand_SingleUniverse(muon.muon_theta_response, err_name, i, reco_muon_theta * TMath::RadToDeg(), truth_muon_theta * TMath::RadToDeg(), 0.0, 0.0);
            }else{
                FillLatErrorBand_SingleUniverse(muon.muon_theta_mc_reco_bckg, err_name, i, reco_muon_theta * TMath::RadToDeg(), 0.0);
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
            // Did not passed cuts
        }
    }
}

void CCProtonPi0_Analyzer::FillLatErrorBand_EM_EnergyScale()
{
    std::string err_name = "EM_EnergyScale";
    double reco_muon_theta = GetCorrectedMuonTheta();
    
    for (int i = 0; i < n_lateral_universes; ++i){

        // EM Energy Dependent Variables
        double gamma1_E_i = (1.0 + em_energy_random_shifts[i]) * gamma1_E;
        double gamma2_E_i = (1.0 + em_energy_random_shifts[i]) * gamma2_E;
        double pi0_invMass_i = std::sqrt(2*gamma1_E_i*gamma2_E_i*(1-pi0_cos_openingAngle));
        double pi0_P_i = (1.0 + em_energy_random_shifts[i]) * pi0_P;
        double pi0_E_i = sqrt(pi0_P_i*pi0_P_i + pi0_mass*pi0_mass);
        double pi0_KE_i = pi0_E_i - pi0_mass;
        double Enu_i = Calc_Enu_shifted(muon_E, pi0_E_i, m_total_proton_KE); // Use actual muon energy and shifted pi0 energy
        double QSq_i = Calc_QSq(Enu_i, muon_E, muon_P, reco_muon_theta); // Use actual muon energy
        double WSq_i = Calc_WSq(Enu_i, QSq_i, muon_E); // Use actual muon energy
        double W_i = (WSq_i > 0) ? sqrt(WSq_i) : -1; 

        bool PassedCuts = IsEnuInRange(Enu_i) && IsWInRange(W_i) && IsInvMassInRange(pi0_invMass_i) && !IsOpeningAngleSmallAndEnergyLow(gamma1_E_i, gamma2_E_i);

        if (PassedCuts){
            double pi0_P_shift = (pi0_P_i - pi0_P) * MeV_to_GeV;
            double pi0_KE_shift = (pi0_KE_i - pi0_KE) * MeV_to_GeV;
            double Enu_shift = (Enu_i - m_Enu) * MeV_to_GeV;
            double QSq_shift = (QSq_i - m_QSq) * MeVSq_to_GeVSq;
            double W_shift = (W_i - m_W) * MeV_to_GeV;

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
            FillLatErrorBand_SingleUniverse(muon.muon_theta_mc_reco_all, err_name, i, reco_muon_theta * TMath::RadToDeg(), 0.0);
            if (truth_isSignal){
                FillLatErrorBand_SingleUniverse(muon.muon_theta_mc_truth_signal, err_name, i, truth_muon_theta * TMath::RadToDeg(), 0.0);
                FillLatErrorBand_SingleUniverse(muon.muon_theta_mc_reco_signal, err_name, i, reco_muon_theta * TMath::RadToDeg(), 0.0);
                FillLatErrorBand_SingleUniverse(muon.muon_theta_response, err_name, i, reco_muon_theta * TMath::RadToDeg(), truth_muon_theta * TMath::RadToDeg(), 0.0, 0.0);
            }else{
                FillLatErrorBand_SingleUniverse(muon.muon_theta_mc_reco_bckg, err_name, i, reco_muon_theta * TMath::RadToDeg(), 0.0);
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
            // Does not satisfy cuts
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
        double muon_P_i = (1.0 + muonP_random_shifts[i]) * muon_P;
        double muon_E_i = sqrt(muon_P_i * muon_P_i + muon_mass * muon_mass);
        double Enu_i = Calc_Enu_shifted(muon_E_i, pi0_E, m_total_proton_KE); // Use shifted muon energy and actual pi0 energy
        double QSq_i = Calc_QSq(Enu_i, muon_E_i, muon_P_i, reco_muon_theta); // Use shifted muon energy
        double WSq_i = Calc_WSq(Enu_i, QSq_i, muon_E_i); // Use shifted muon energy
        double W_i = (WSq_i > 0) ? sqrt(WSq_i) : -1; 

        bool PassedCuts = IsEnuInRange(Enu_i) && IsWInRange(W_i);

        if (PassedCuts){
            double muon_P_shift = (muon_P_i - muon_P) * MeV_to_GeV;
            double Enu_shift = (Enu_i - m_Enu) * MeV_to_GeV;
            double QSq_shift = (QSq_i - m_QSq) * MeVSq_to_GeVSq;
            double W_shift = (W_i - m_W) * MeV_to_GeV;

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
            // Does not satisfy cuts
        }
    }
}

void CCProtonPi0_Analyzer::FillLatErrorBand_MuonTheta()
{
    std::string err_name = "MuonTheta";
    double reco_muon_theta = GetCorrectedMuonTheta();

    for (int i = 0; i < n_lateral_universes; ++i){
        
        double muon_theta_i = (1.0 + muon_theta_random_shifts[i]) * reco_muon_theta;
        double QSq_i = Calc_QSq(m_Enu, muon_E, muon_P, muon_theta_i); // Use actual neutrino energy and shifted muon theta
        double WSq_i = Calc_WSq(m_Enu, QSq_i, muon_E); // Use shifted QSq 
        double W_i = (WSq_i > 0) ? sqrt(WSq_i) : -1; 
        
        bool PassedCuts =  IsWInRange(W_i);

        if (PassedCuts){

            double muon_theta_shift = (muon_theta_i - reco_muon_theta) * TMath::RadToDeg();
            double QSq_shift = (QSq_i - m_QSq) * MeVSq_to_GeVSq;
            double W_shift = (W_i - m_W) * MeV_to_GeV;

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
            FillLatErrorBand_SingleUniverse(muon.muon_theta_mc_reco_all, err_name, i, reco_muon_theta * TMath::RadToDeg(), muon_theta_shift);
            if (truth_isSignal){
                FillLatErrorBand_SingleUniverse(muon.muon_theta_mc_truth_signal, err_name, i, truth_muon_theta * TMath::RadToDeg(), 0.0);
                FillLatErrorBand_SingleUniverse(muon.muon_theta_mc_reco_signal, err_name, i, reco_muon_theta * TMath::RadToDeg(), muon_theta_shift);
                FillLatErrorBand_SingleUniverse(muon.muon_theta_response, err_name, i, reco_muon_theta * TMath::RadToDeg(), truth_muon_theta * TMath::RadToDeg(), muon_theta_shift, 0.0);
            }else{
                FillLatErrorBand_SingleUniverse(muon.muon_theta_mc_reco_bckg, err_name, i, reco_muon_theta * TMath::RadToDeg(), muon_theta_shift);
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
            FillLatErrorBand_SingleUniverse(interaction.Enu_mc_reco_all, err_name, i, m_Enu * MeV_to_GeV, 0.0);
            if (truth_isSignal){
                FillLatErrorBand_SingleUniverse(interaction.Enu_mc_truth_signal, err_name, i, m_Enu_Truth * MeV_to_GeV, 0.0);
                FillLatErrorBand_SingleUniverse(interaction.Enu_mc_reco_signal, err_name, i, m_Enu * MeV_to_GeV, 0.0);
                FillLatErrorBand_SingleUniverse(interaction.Enu_response, err_name, i, m_Enu * MeV_to_GeV, m_Enu_Truth * MeV_to_GeV, 0.0, 0.0);
            }else{
                FillLatErrorBand_SingleUniverse(interaction.Enu_mc_reco_bckg, err_name, i, m_Enu * MeV_to_GeV, 0.0);
            }
        }else{
            // Did NOT passed Cuts!
        }
    }
}

void CCProtonPi0_Analyzer::FillLatErrorBand_ProtonEnergy_invMass(std::string err_name)
{
    bool PassedCuts;
    double reco_muon_theta = GetCorrectedMuonTheta();

    for (int i = 0; i < 2; ++i){
        if (nProtonCandidates > 0){
            std::vector<double> proton_energy_shifts = GetProtonEnergyShifts(err_name, i);

            // Calculate Event Kinematics with Shifted Proton Energy 
            double total_proton_KE_i = 0.0;
            for (int p = 0; p < nProtonCandidates; ++p ){
                total_proton_KE_i += all_protons_KE[p] + proton_energy_shifts[p]; 
            }
            double Enu_i = Calc_Enu_shifted(muon_E, pi0_E, total_proton_KE_i); // Use actual muon energy and pi0 energy but shifted proton kinetic energy
            double QSq_i = Calc_QSq(Enu_i, muon_E, muon_P, reco_muon_theta); // Use actual muon energy
            double WSq_i = Calc_WSq(Enu_i, QSq_i, muon_E); // Use actual muon energy
            double W_i = (WSq_i > 0) ? sqrt(WSq_i) : -1; 

            PassedCuts = IsEnuInRange(Enu_i) && IsWInRange(W_i);

        }else{
            // If there are no protons, event passes cuts
            PassedCuts = true;
        }
        
        if (PassedCuts){
            FillLatErrorBand_SingleUniverse(cutList.invMass_mc_reco_all, err_name, i, pi0_invMass, 0.0);
            if (truth_isSignal){
                FillLatErrorBand_SingleUniverse(cutList.invMass_mc_reco_signal, err_name, i, pi0_invMass, 0.0);
            }else{
                FillLatErrorBand_SingleUniverse(cutList.invMass_mc_reco_bckg, err_name, i, pi0_invMass, 0.0);
            }
        }else{
            // Did Not passed cuts
        }
    }
}

void CCProtonPi0_Analyzer::FillLatErrorBand_ProtonEnergy_Birks_invMass()
{
    std::string err_name = "ProtonEnergy_Birks";
  
    bool PassedCuts;

    // For each event Birks_random_shifts is different
    Calc_Birks_random_shifts();
    double reco_muon_theta = GetCorrectedMuonTheta();

    for (int i = 0; i < n_lateral_universes; ++i){

        if (nProtonCandidates > 0){
        // Calculate Event Kinematics with Shifted Proton Energy 
        double total_proton_KE_i = 0.0;
        for (int p = 0; p < nProtonCandidates; ++p ){
            total_proton_KE_i += (1.0 + Birks_random_shifts2D[p][i]) * all_protons_KE[p]; 
        }
        double Enu_i = Calc_Enu_shifted(muon_E, pi0_E, total_proton_KE_i); // Use actual muon energy and pi0 energy but shifted proton kinetic energy
        double QSq_i = Calc_QSq(Enu_i, muon_E, muon_P, reco_muon_theta); // Use actual muon energy
        double WSq_i = Calc_WSq(Enu_i, QSq_i, muon_E); // Use actual muon energy
        double W_i = (WSq_i > 0) ? sqrt(WSq_i) : -1; 

        PassedCuts = IsEnuInRange(Enu_i) && IsWInRange(W_i);

        }else{
            // If there are no protons, event passes cuts
            PassedCuts = true;
        }
        
        if (PassedCuts){
            FillLatErrorBand_SingleUniverse(cutList.invMass_mc_reco_all, err_name, i, pi0_invMass, 0.0);
            if (truth_isSignal){
                FillLatErrorBand_SingleUniverse(cutList.invMass_mc_reco_signal, err_name, i, pi0_invMass, 0.0);
            }else{
                FillLatErrorBand_SingleUniverse(cutList.invMass_mc_reco_bckg, err_name, i, pi0_invMass, 0.0);
            }
        }else{
            // Did Not passed cuts
        }
    }
}

void CCProtonPi0_Analyzer::FillLatErrorBand_EM_EnergyScale_invMass()
{
    std::string err_name = "EM_EnergyScale";
    
    double reco_muon_theta = GetCorrectedMuonTheta();
    for (int i = 0; i < n_lateral_universes; ++i){

        double gamma1_E_i = (1.0 + em_energy_random_shifts[i]) * gamma1_E;
        double gamma2_E_i = (1.0 + em_energy_random_shifts[i]) * gamma2_E;
        double pi0_invMass_i = std::sqrt(2*gamma1_E_i*gamma2_E_i*(1-pi0_cos_openingAngle));
        double pi0_P_i = (1.0 + em_energy_random_shifts[i]) * pi0_P;
        double pi0_E_i = sqrt(pi0_P_i*pi0_P_i + pi0_mass*pi0_mass);
        double Enu_i = Calc_Enu_shifted(muon_E, pi0_E_i, m_total_proton_KE); // Use actual muon energy and shifted pi0 energy
        double QSq_i = Calc_QSq(Enu_i, muon_E, muon_P, reco_muon_theta); // Use actual muon energy
        double WSq_i = Calc_WSq(Enu_i, QSq_i, muon_E); // Use actual muon energy
        double W_i = (WSq_i > 0) ? sqrt(WSq_i) : -1; 

        bool PassedCuts = IsEnuInRange(Enu_i) && IsWInRange(W_i) && IsInvMassInRange(pi0_invMass_i) && !IsOpeningAngleSmallAndEnergyLow(gamma1_E_i, gamma2_E_i);

        if (PassedCuts){
            double pi0_invMass_shift = pi0_invMass_i - pi0_invMass;
            FillLatErrorBand_SingleUniverse(cutList.invMass_mc_reco_all, err_name, i, pi0_invMass, pi0_invMass_shift);
            if (truth_isSignal){
                FillLatErrorBand_SingleUniverse(cutList.invMass_mc_reco_signal, err_name, i, pi0_invMass, pi0_invMass_shift);
            }else{
                FillLatErrorBand_SingleUniverse(cutList.invMass_mc_reco_bckg, err_name, i, pi0_invMass, pi0_invMass_shift);
            }
        }else{
            // Does not satisfy cuts
        }
    }
}

void CCProtonPi0_Analyzer::FillLatErrorBand_MuonMomentum_invMass()
{
    std::string err_name = "MuonMomentum";
    
    // For each event muonP_random_shifts is different
    Calc_muonP_random_shifts();
    double reco_muon_theta = GetCorrectedMuonTheta();
 
    for (int i = 0; i < n_lateral_universes; ++i){

        // Muon Variables
        double muon_P_i = (1.0 + muonP_random_shifts[i]) * muon_P;
        double muon_E_i = sqrt(muon_P_i * muon_P_i + muon_mass * muon_mass);
        double Enu_i = Calc_Enu_shifted(muon_E_i, pi0_E, m_total_proton_KE); // Use shifted muon energy and actual pi0 energy
        double QSq_i = Calc_QSq(Enu_i, muon_E_i, muon_P_i, reco_muon_theta); // Use shifted muon energy
        double WSq_i = Calc_WSq(Enu_i, QSq_i, muon_E_i); // Use shifted muon energy
        double W_i = (WSq_i > 0) ? sqrt(WSq_i) : -1; 

        bool PassedCuts = IsEnuInRange(Enu_i) && IsWInRange(W_i);

        if (PassedCuts){
            FillLatErrorBand_SingleUniverse(cutList.invMass_mc_reco_all, err_name, i, pi0_invMass, 0.0);
            if (truth_isSignal){
                FillLatErrorBand_SingleUniverse(cutList.invMass_mc_reco_signal, err_name, i, pi0_invMass, 0.0);
            }else{
                FillLatErrorBand_SingleUniverse(cutList.invMass_mc_reco_bckg, err_name, i, pi0_invMass, 0.0);
            }
        }else{
            // Does not satisfy cuts
        }
    }
}

void CCProtonPi0_Analyzer::FillLatErrorBand_MuonTheta_invMass()
{
    std::string err_name = "MuonTheta";
    double reco_muon_theta = GetCorrectedMuonTheta();
    
    for (int i = 0; i < n_lateral_universes; ++i){

        double muon_theta_i = (1.0 + muon_theta_random_shifts[i]) * reco_muon_theta;
        double QSq_i = Calc_QSq(m_Enu, muon_E, muon_P, muon_theta_i); // Use actual neutrino energy and shifted muon theta
        double WSq_i = Calc_WSq(m_Enu, QSq_i, muon_E); // Use shifted QSq 
        double W_i = (WSq_i > 0) ? sqrt(WSq_i) : -1; 

        bool PassedCuts = IsWInRange(W_i);

        if (PassedCuts){
            FillLatErrorBand_SingleUniverse(cutList.invMass_mc_reco_all, err_name, i, pi0_invMass, 0.0);
            if (truth_isSignal){
                FillLatErrorBand_SingleUniverse(cutList.invMass_mc_reco_signal, err_name, i, pi0_invMass, 0.0);
            }else{
                FillLatErrorBand_SingleUniverse(cutList.invMass_mc_reco_bckg, err_name, i, pi0_invMass, 0.0);
            }
        }else{
            // Not satisfied cuts
        }
    }
}


void CCProtonPi0_Analyzer::FillLatErrorBand_ProtonEnergy_SideBand_invMass(std::string err_name)
{
    bool PassedCuts;
    double reco_muon_theta = GetCorrectedMuonTheta();

    for (int i = 0; i < 2; ++i){
        if (nProtonCandidates > 0){
            std::vector<double> proton_energy_shifts = GetProtonEnergyShifts(err_name, i);

            // Calculate Event Kinematics with Shifted Proton Energy 
            double total_proton_KE_i = 0.0;
            for (int p = 0; p < nProtonCandidates; ++p ){
                total_proton_KE_i += all_protons_KE[p] + proton_energy_shifts[p]; 
            }
            double Enu_i = Calc_Enu_shifted(muon_E, pi0_E, total_proton_KE_i); // Use actual muon energy and pi0 energy but shifted proton kinetic energy
            double QSq_i = Calc_QSq(Enu_i, muon_E, muon_P, reco_muon_theta); // Use actual muon energy
            double WSq_i = Calc_WSq(Enu_i, QSq_i, muon_E); // Use actual muon energy
            double W_i = (WSq_i > 0) ? sqrt(WSq_i) : -1; 

            PassedCuts = IsEnuInRange(Enu_i) && IsWInRange(W_i);

        }else{
            // For No Proton Events, event passes all cuts
            PassedCuts = true;
        }

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

void CCProtonPi0_Analyzer::FillLatErrorBand_ProtonEnergy_Birks_SideBand_invMass()
{
    std::string err_name = "ProtonEnergy_Birks";

    bool PassedCuts;
    double reco_muon_theta = GetCorrectedMuonTheta();
    
    // For each event Birks_random_shifts is different
    Calc_Birks_random_shifts();

    for (int i = 0; i < n_lateral_universes; ++i){
        
        if (nProtonCandidates > 0){
            // Calculate Event Kinematics with Shifted Proton Energy 
            double total_proton_KE_i = 0.0;
            for (int p = 0; p < nProtonCandidates; ++p ){
                total_proton_KE_i += (1.0 + Birks_random_shifts2D[p][i]) * all_protons_KE[p]; 
            }
            double Enu_i = Calc_Enu_shifted(muon_E, pi0_E, total_proton_KE_i); // Use actual muon energy and pi0 energy but shifted proton kinetic energy
            double QSq_i = Calc_QSq(Enu_i, muon_E, muon_P, reco_muon_theta); // Use actual muon energy
            double WSq_i = Calc_WSq(Enu_i, QSq_i, muon_E); // Use actual muon energy
            double W_i = (WSq_i > 0) ? sqrt(WSq_i) : -1; 

            PassedCuts = IsEnuInRange(Enu_i) && IsWInRange(W_i);

        }else{
            // For No Proton Events, event passes all cuts
            PassedCuts = true;
        }

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

void CCProtonPi0_Analyzer::FillLatErrorBand_EM_EnergyScale_SideBand_invMass()
{
    std::string err_name = "EM_EnergyScale";
    double reco_muon_theta = GetCorrectedMuonTheta();
    
    for (int i = 0; i < n_lateral_universes; ++i){

        // EM Energy Dependent Variables
        double gamma1_E_i = (1.0 + em_energy_random_shifts[i]) * gamma1_E;
        double gamma2_E_i = (1.0 + em_energy_random_shifts[i]) * gamma2_E;
        double pi0_invMass_i = std::sqrt(2*gamma1_E_i*gamma2_E_i*(1-pi0_cos_openingAngle));
        double pi0_P_i = (1.0 + em_energy_random_shifts[i]) * pi0_P;
        double pi0_E_i = sqrt(pi0_P_i*pi0_P_i + pi0_mass*pi0_mass);
        double Enu_i = Calc_Enu_shifted(muon_E, pi0_E_i, m_total_proton_KE); // Use actual muon energy and shifted pi0 energy
        double QSq_i = Calc_QSq(Enu_i, muon_E, muon_P, reco_muon_theta); // Use actual muon energy
        double WSq_i = Calc_WSq(Enu_i, QSq_i, muon_E); // Use actual muon energy
        double W_i = (WSq_i > 0) ? sqrt(WSq_i) : -1; 

        bool PassedCuts = IsEnuInRange(Enu_i) && IsWInRange(W_i) && !IsOpeningAngleSmallAndEnergyLow(gamma1_E_i, gamma2_E_i);

        if (PassedCuts){
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
    double reco_muon_theta = GetCorrectedMuonTheta();
 
    for (int i = 0; i < n_lateral_universes; ++i){

        double muon_P_i = (1.0 + muonP_random_shifts[i]) * muon_P;
        double muon_E_i = sqrt(muon_P_i * muon_P_i + muon_mass * muon_mass);
        double Enu_i = Calc_Enu_shifted(muon_E_i, pi0_E, m_total_proton_KE); // Use shifted muon energy and actual pi0 energy
        double QSq_i = Calc_QSq(Enu_i, muon_E_i, muon_P_i, reco_muon_theta); // Use shifted muon energy
        double WSq_i = Calc_WSq(Enu_i, QSq_i, muon_E_i); // Use shifted muon energy
        double W_i = (WSq_i > 0) ? sqrt(WSq_i) : -1; 

        bool PassedCuts = IsEnuInRange(Enu_i) && IsWInRange(W_i);

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

void CCProtonPi0_Analyzer::FillLatErrorBand_MuonTheta_SideBand_invMass()
{
    std::string err_name = "MuonTheta";
    
    double reco_muon_theta = GetCorrectedMuonTheta();
    
    for (int i = 0; i < n_lateral_universes; ++i){

        double muon_theta_i = (1.0 + muon_theta_random_shifts[i]) * reco_muon_theta;
        double QSq_i = Calc_QSq(m_Enu, muon_E, muon_P, muon_theta_i); // Use actual neutrino energy and shifted muon theta
        double WSq_i = Calc_WSq(m_Enu, QSq_i, muon_E); // Use shifted QSq 
        double W_i = (WSq_i > 0) ? sqrt(WSq_i) : -1; 

        bool PassedCuts = IsWInRange(W_i);

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
        }else{
            // Does NOT satisfy cuts
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

    // Get Weight for the specific universe
    // wgt_bckg is universe_wgt / cv_wgt
    double wgt_bckg = applyBckgConstraints_Unv ? GetBckgConstraint(err_name, unv) : 1.0;
    double wgtU = cvweight * wgt_bckg;
    err_hists[unv]->AddBinContent( bin, wgtU );

    // Update Bin Error
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
    double wgt_bckg = applyBckgConstraints_Unv ? GetBckgConstraint(err_name, unv) : 1.0;
    double wgtU = cvweight * wgt_bckg;
    err_hists[unv]->AddBinContent( bin, wgtU );

//    const double err = err_hists[unv]->GetBinError(bin);
//    const double newerr2 = err*err + wgtU*wgtU;
//    const double newerr = (0.<newerr2) ? sqrt(newerr2) : 0.;
//    err_hists[unv]->SetBinError( bin, newerr );

}

double CCProtonPi0_Analyzer::GetBckgConstraint(std::string error_name, int hist_ind)
{
    // Find the Bckg Constraint if it is one of the constrained events
    if (truth_isBckg_Compact_SinglePiPlus){
        double bckg_wgt_SinglePiPlus = BckgConstrainer.GetBckgConstraint(error_name, hist_ind, "SinglePiPlus");
        return bckg_wgt_SinglePiPlus / cv_wgt_SinglePiPlus;
    }else if (truth_isBckg_Compact_QELike){
        double bckg_wgt_QELike = BckgConstrainer.GetBckgConstraint(error_name, hist_ind, "QELike");
        return bckg_wgt_QELike / cv_wgt_QELike;
    }else if (truth_isBckg_Compact_WithPi0){
        double bckg_wgt_WithPi0 = BckgConstrainer.GetBckgConstraint(error_name, hist_ind, "WithPi0");
        return bckg_wgt_WithPi0 / cv_wgt_WithPi0;
    }else{
        return 1.0;
    }
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

double CCProtonPi0_Analyzer::GetProtonTrackingErr()
{
    double track_eff_error = 0.0;
    if (nProtonCandidates > 0){
        track_eff_error = proton_length < 8*17.0 ? 0.046 : 0.003; //Ref: TN048
    }
    return track_eff_error;
}

std::vector<double> CCProtonPi0_Analyzer::GetProtonEnergyShifts(std::string err_name, int unv)
{
    std::vector<double> energy_shifts;

    if (err_name.compare("ProtonEnergy_MassModel") == 0){
        if (unv == 0) FillProtonEnergyShiftVector(energy_shifts, all_protons_energy_shift_Mass_Up);
        else if (unv == 1) FillProtonEnergyShiftVector(energy_shifts, all_protons_energy_shift_Mass_Down);
        else RunTimeError("Proton Systematics, Wrong Number of Universes!");
    }else if (err_name.compare("ProtonEnergy_MEU") == 0){
        if (unv == 0) FillProtonEnergyShiftVector(energy_shifts, all_protons_energy_shift_MEU_Up);
        else if (unv == 1) FillProtonEnergyShiftVector(energy_shifts, all_protons_energy_shift_MEU_Down);
        else RunTimeError("Proton Systematics, Wrong Number of Universes!");
    }else if (err_name.compare("ProtonEnergy_BetheBloch") == 0){
        if (unv == 0) FillProtonEnergyShiftVector(energy_shifts, all_protons_energy_shift_BetheBloch_Up);
        else if (unv == 1) FillProtonEnergyShiftVector(energy_shifts, all_protons_energy_shift_BetheBloch_Down);
        else RunTimeError("Proton Systematics, Wrong Number of Universes!");
    }else{
        RunTimeError("Proton Systematics, Wrong Error Band Name!");
    }

    return energy_shifts;
}

void CCProtonPi0_Analyzer::FillProtonEnergyShiftVector(std::vector<double> &energy_shifts, Double_t shifts[10])
{
    for (int i = 0; i < nProtonCandidates; ++i){
        energy_shifts.push_back(shifts[i]);
    }
}

double CCProtonPi0_Analyzer::GetBckgConstraint_WithPi0Err()
{
    double correctionErr;

    if (truth_isBckg_Compact_WithPi0) correctionErr = cv_err_WithPi0;
    else correctionErr = 0.0;

    return correctionErr;
}

double CCProtonPi0_Analyzer::GetBckgConstraint_SinglePiPlusErr()
{
    double correctionErr;

    if (truth_isBckg_Compact_SinglePiPlus) correctionErr = cv_err_SinglePiPlus;
    else correctionErr = 0.0;

    return correctionErr;
}

double CCProtonPi0_Analyzer::GetBckgConstraint_QELikeErr()
{
    double correctionErr;

    if (truth_isBckg_Compact_QELike) correctionErr = cv_err_QELike;
    else correctionErr = 0.0;

    return correctionErr;
}

double CCProtonPi0_Analyzer::GetMichelFakeErr()
{
    // Uncertainty is 0.5% hence weights 1-0.005 AND 1+0.005
    // docDB 11443, slide 38 
    
    // Michel Fake Error Applied only to Events without True Michel
    double correctionErr = 0.0;
    if (!truth_isBckg_withMichel) correctionErr = 0.005;
    return correctionErr;
}

double CCProtonPi0_Analyzer::GetMichelTrueErr()
{
    // Uncertainty is 1.1% hence weights 1-0.011 AND 1+0.011
    // docDB 11443, slide 38 
    
    // Michel True Error Applied only to Events with True Michel
    double correctionErr = 0.0;
    if (truth_isBckg_withMichel) correctionErr = 0.011;
    return correctionErr;
}

double CCProtonPi0_Analyzer::GetTargetMassErr()
{
    // Uncertainty is 1.4% hence weights 1-0.014 AND 1+0.014
    return 0.014;
}

void CCProtonPi0_Analyzer::Get2p2hErr()
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

//    if (err_2p2h != 0.0){
//        std::cout<<" "<<std::endl;
//        std::cout<<"err_np = "<<err_np<<std::endl;
//        std::cout<<"err_nn = "<<err_nn<<std::endl; 
//        std::cout<<"err_2p2h = "<<err_2p2h<<std::endl;
//    }
}

double CCProtonPi0_Analyzer::GetUnfoldingErr()
{
    // Uncertainty is 1% hence weights 1-0.01 AND 1+0.01
    // docDB 12426, slide 22
    return 0.01;
}

double CCProtonPi0_Analyzer::GetPionResponseErr()
{
    double correctionErr = 0.0;

    if (truth_isBckg_SingleChargedPion_ChargeExchanged){
        correctionErr = 0.5;
    }

    return correctionErr;
}

void CCProtonPi0_Analyzer::Calc_muon_theta_random_shifts()
{
    // Calculate Muon Theta Error from ThetaX and ThetaY Error 
    const double muonThetaX_Err = 0.001;  //rad (error on beam angleX correction from Edgar Valencia's DocDB 11550) 
    const double muonThetaY_Err = 0.0009; //rad (error on beam angleY correction from Edgar Valencia's DocDB 11550)
    
    // This is merely: sec^2(theta_X) + sec^2(theta_Y) = sec^2(theta) +1 
    double cos_muonTheta_Err = 1/(sqrt(1.0/pow(cos(muonThetaX_Err), 2) + 1.0/pow(cos(muonThetaY_Err), 2) -1) ); 
    double muonTheta_Err = acos(cos_muonTheta_Err);
    muon_theta_random_shifts = RandNumGenerator.GetRandomShifts(muonTheta_Err); // ~ Gaussian(0.0, muonTheta_Err)
}

void CCProtonPi0_Analyzer::Calc_em_energy_random_shifts()
{
    // Fill EM Energy Scale Shifts
    double mc_1sigma = 0.013;
    double data_1sigma = 0.019; 
    double em_uncertainty = sqrt(mc_1sigma*mc_1sigma + data_1sigma*data_1sigma);
    em_energy_random_shifts = RandNumGenerator.GetRandomShifts(em_uncertainty); // ~ Gaussian(0.0, em_uncertainty)
}

void CCProtonPi0_Analyzer::Calc_no_random_shifts()
{
    for (int i = 0; i < n_lateral_universes; ++i){
        no_random_shifts.push_back(0.0);
    }   
}

void CCProtonPi0_Analyzer::Fill_RandomShiftHistograms()
{
    // Fill Histograms for Constant Shifts
    std::vector<double> normal_rand_numbers = RandNumGenerator.GetNormalRandomVector();
    for (unsigned int i = 0; i < normal_rand_numbers.size(); ++i){
        interaction.normal_rand_numbers->Fill(normal_rand_numbers[i]);
        interaction.em_shift_rand_numbers->Fill(em_energy_random_shifts[i]);
        interaction.muon_theta_shift_rand_numbers->Fill(muon_theta_random_shifts[i]);
    }
}

void CCProtonPi0_Analyzer::initUpdatedGenieWeights()
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

void CCProtonPi0_Analyzer::UpdateGENIESystematics()
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

#endif //CCProtonPi0_Analyzer_Systematics_cpp
