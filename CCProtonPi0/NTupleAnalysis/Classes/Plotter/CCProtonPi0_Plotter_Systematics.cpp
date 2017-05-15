#include "CCProtonPi0_Plotter.h"

using namespace PlotUtils;

void CCProtonPi0_Plotter::Systematics_CheckErrorSummary(std::string root_dir, std::string var_name)
{
    std::string plotDir = Folder_List::plotDir_Systematics_Summary;
    TFile* f_mc = new TFile(root_dir.c_str());

    MnvH1D* mc = GetMnvH1D(f_mc, var_name);
    
    DrawErrorSummary(mc, var_name, plotDir);

    Systematics_WriteTable_Fraction(mc, var_name);
    Systematics_WriteTable_BinByBin(mc, var_name);

    delete mc;
}

void CCProtonPi0_Plotter::Systematics_invMass()
{
    Systematics_DrawErrorSummary("invMass_all", "invMass_mc_reco_all");
}

void CCProtonPi0_Plotter::Systematics_ErrorEvolution(std::string var_name, std::string err_name, int unv)
{
    Systematics_SingleUniverse(var_name, "bckg_subtracted", err_name, unv); 
    Systematics_SingleUniverse(var_name, "unfolded", err_name, unv); 
    Systematics_SingleUniverse(var_name, "efficiency_corrected", err_name, unv); 
    Systematics_SingleUniverse(var_name, "xsec", err_name, unv); 
}

void CCProtonPi0_Plotter::Systematics_SingleUniverse(std::string var_name, std::string data_type, std::string err_name, int unv)
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;

    TFile* f = new TFile(rootDir_CrossSection.data.c_str());
    std::string var = var_name + "_" + data_type;

    MnvH1D* hist = GetMnvH1D(f, var);
    double norm_bin_width = hist->GetNormBinWidth();

    TH1D* h = NULL;

    if (err_name.compare("CV") == 0){
        h = GetBinNormalizedTH1D(hist); 
    }else{
        MnvVertErrorBand* err = hist->GetVertErrorBand(err_name);
        h = err->GetHist(unv);
        h->Scale(norm_bin_width, "width");
    }

    std::string plot_name = var + "_" + err_name;
    Draw1DHist(h, plot_name, plotDir);

}

void CCProtonPi0_Plotter::Systematics_XSec()
{
    Systematics_DrawErrorSummary("muon_P_xsec", "muon_P_xsec");
    Systematics_DrawErrorSummary("muon_theta_xsec", "muon_theta_xsec");
    Systematics_DrawErrorSummary("pi0_P_xsec", "pi0_P_xsec");
    Systematics_DrawErrorSummary("pi0_KE_xsec", "pi0_KE_xsec");
    Systematics_DrawErrorSummary("pi0_theta_xsec", "pi0_theta_xsec");
    Systematics_DrawErrorSummary("QSq_xsec", "QSq_xsec");
    Systematics_DrawErrorSummary("Enu_xsec", "Enu_xsec");
    //Systematics_DrawErrorSummary("W_xsec", "W_xsec");

    // Group I -- Detector
    Systematics_DrawErrorSummary_Group("muon_P_xsec", "muon_P_xsec", detGroup, "Detector");
    Systematics_DrawErrorSummary_Group("muon_theta_xsec", "muon_theta_xsec", detGroup, "Detector");
    Systematics_DrawErrorSummary_Group("pi0_P_xsec", "pi0_P_xsec", detGroup, "Detector");
    Systematics_DrawErrorSummary_Group("pi0_KE_xsec", "pi0_KE_xsec", detGroup, "Detector");
    Systematics_DrawErrorSummary_Group("pi0_theta_xsec", "pi0_theta_xsec", detGroup, "Detector");
    Systematics_DrawErrorSummary_Group("QSq_xsec", "QSq_xsec", detGroup, "Detector");
    Systematics_DrawErrorSummary_Group("Enu_xsec", "Enu_xsec", detGroup, "Detector");
    //Systematics_DrawErrorSummary_Group("W_xsec", "W_xsec", detGroup, "Detector");

    // Group II -- GENIE Cross Section 
    Systematics_DrawErrorSummary_Group("muon_P_xsec", "muon_P_xsec", genieGroup, "GENIE_XSec");
    Systematics_DrawErrorSummary_Group("muon_theta_xsec", "muon_theta_xsec", genieGroup, "GENIE_XSec");
    Systematics_DrawErrorSummary_Group("pi0_P_xsec", "pi0_P_xsec", genieGroup, "GENIE_XSec");
    Systematics_DrawErrorSummary_Group("pi0_KE_xsec", "pi0_KE_xsec", genieGroup, "GENIE_XSec");
    Systematics_DrawErrorSummary_Group("pi0_theta_xsec", "pi0_theta_xsec", genieGroup, "GENIE_XSec");
    Systematics_DrawErrorSummary_Group("QSq_xsec", "QSq_xsec", genieGroup, "GENIE_XSec");
    Systematics_DrawErrorSummary_Group("Enu_xsec", "Enu_xsec", genieGroup, "GENIE_XSec");
    //Systematics_DrawErrorSummary_Group("W_xsec", "W_xsec", genieGroup, "GENIE_XSec");

    // Group III -- GENIE FSI 
    Systematics_DrawErrorSummary_Group("muon_P_xsec", "muon_P_xsec", fsiGroup, "GENIE_FSI");
    Systematics_DrawErrorSummary_Group("muon_theta_xsec", "muon_theta_xsec", fsiGroup, "GENIE_FSI");
    Systematics_DrawErrorSummary_Group("pi0_P_xsec", "pi0_P_xsec", fsiGroup, "GENIE_FSI");
    Systematics_DrawErrorSummary_Group("pi0_KE_xsec", "pi0_KE_xsec", fsiGroup, "GENIE_FSI");
    Systematics_DrawErrorSummary_Group("pi0_theta_xsec", "pi0_theta_xsec", fsiGroup, "GENIE_FSI");
    Systematics_DrawErrorSummary_Group("QSq_xsec", "QSq_xsec", fsiGroup, "GENIE_FSI");
    Systematics_DrawErrorSummary_Group("Enu_xsec", "Enu_xsec", fsiGroup, "GENIE_FSI");
    //Systematics_DrawErrorSummary_Group("W_xsec", "W_xsec", fsiGroup, "GENIE_FSI");

    // Group IV -- Flux 
    Systematics_DrawErrorSummary_Group("muon_P_xsec", "muon_P_xsec", fluxGroup, "Flux");
    Systematics_DrawErrorSummary_Group("muon_theta_xsec", "muon_theta_xsec", fluxGroup, "Flux");
    Systematics_DrawErrorSummary_Group("pi0_P_xsec", "pi0_P_xsec", fluxGroup, "Flux");
    Systematics_DrawErrorSummary_Group("pi0_KE_xsec", "pi0_KE_xsec", fluxGroup, "Flux");
    Systematics_DrawErrorSummary_Group("pi0_theta_xsec", "pi0_theta_xsec", fluxGroup, "Flux");
    Systematics_DrawErrorSummary_Group("QSq_xsec", "QSq_xsec", fluxGroup, "Flux");
    Systematics_DrawErrorSummary_Group("Enu_xsec", "Enu_xsec", fluxGroup, "Flux");
    //Systematics_DrawErrorSummary_Group("W_xsec", "W_xsec", fluxGroup, "Flux");

    // Group V -- Other 
    Systematics_DrawErrorSummary_Group("muon_P_xsec", "muon_P_xsec", otherGroup, "Other");
    Systematics_DrawErrorSummary_Group("muon_theta_xsec", "muon_theta_xsec", otherGroup, "Other");
    Systematics_DrawErrorSummary_Group("pi0_P_xsec", "pi0_P_xsec", otherGroup, "Other");
    Systematics_DrawErrorSummary_Group("pi0_KE_xsec", "pi0_KE_xsec", otherGroup, "Other");
    Systematics_DrawErrorSummary_Group("pi0_theta_xsec", "pi0_theta_xsec", otherGroup, "Other");
    Systematics_DrawErrorSummary_Group("QSq_xsec", "QSq_xsec", otherGroup, "Other");
    Systematics_DrawErrorSummary_Group("Enu_xsec", "Enu_xsec", otherGroup, "Other");
    //Systematics_DrawErrorSummary_Group("W_xsec", "W_xsec", otherGroup, "Other");

    //Systematics_DrawErrorBand_GENIE("muon_P_xsec");
    //Systematics_DrawErrorBand_GENIE("muon_theta_xsec");
    //Systematics_DrawErrorBand_GENIE("pi0_P_xsec");
    //Systematics_DrawErrorBand_GENIE("pi0_KE_xsec");
    //Systematics_DrawErrorBand_GENIE("pi0_theta_xsec");
    //Systematics_DrawErrorBand_GENIE("QSq_xsec");
    ////Systematics_DrawErrorBand_GENIE("W_xsec");
    //Systematics_DrawErrorBand_GENIE("Enu_xsec");
}

void CCProtonPi0_Plotter::Systematics_DrawErrorSummary(std::string data_var, std::string mc_var)
{
    // Silence compiler warning
    (void) mc_var;
    
    std::string plotDir = Folder_List::plotDir_Systematics_Summary;
    TFile* f_xsec_data = new TFile(rootDir_CrossSection.data.c_str());
    //TFile* f_xsec_mc = new TFile(rootDir_CrossSection.mc.c_str());

    MnvH1D* data = GetMnvH1D(f_xsec_data, data_var);
    //MnvH1D* mc = GetMnvH1D(f_xsec_mc, mc_var);
    
    data_var = "data_" + data_var;

    DrawErrorSummary(data, data_var, plotDir);
    //DrawErrorSummary(mc, mc_var, plotDir);

    delete data;
    //delete mc;
}

void CCProtonPi0_Plotter::Systematics_DrawErrorSummary_Group(std::string data_var, std::string mc_var, std::vector<std::string> errGroup, std::string group_name)
{
    // Silence compiler warning
    (void) mc_var;
 
    std::string plotDir = Folder_List::plotDir_Systematics_Summary;
    TFile* f_xsec_data = new TFile(rootDir_CrossSection.data.c_str());

    MnvH1D* data = GetMnvH1D(f_xsec_data, data_var);
   
    // Remove Vertical Error Bands not in the Group
    std::vector<std::string> vert_errors = data->GetVertErrorBandNames();
    
    for (unsigned int i = 0; i < vert_errors.size(); ++i){
        if (!IsErrorInGroup(vert_errors[i], errGroup)){ 
            MnvVertErrorBand* non_group = data->PopVertErrorBand(vert_errors[i]);
            delete non_group;
        }
    }

    // Remove Lateral Error Bands not in the Group
    std::vector<std::string> lat_errors = data->GetLatErrorBandNames();
    
    for (unsigned int i = 0; i < lat_errors.size(); ++i){
        if (!IsErrorInGroup(lat_errors[i], errGroup)){ 
            MnvLatErrorBand* non_group = data->PopLatErrorBand(lat_errors[i]);
            delete non_group;
        }
    }

    data_var = data_var + "_" + group_name;

    DrawErrorSummary_PaperStyle(data, data_var, plotDir, false);

    delete data;
    delete f_xsec_data;
}

void CCProtonPi0_Plotter::Systematics_DrawErrorBand_GENIE(std::string data_var)
{
    TFile* f_xsec_data = new TFile(rootDir_CrossSection.data.c_str());

    MnvH1D* data = GetMnvH1D(f_xsec_data, data_var);

    // First Remove Non-GENIE Error Bands
    std::vector<std::string> vert_errors = data->GetVertErrorBandNames();
    for (unsigned int i = 0; i < vert_errors.size(); ++i){
        std::size_t found = vert_errors[i].find("GENIE");
        if (found == std::string::npos){
            MnvVertErrorBand* non_genie = data->PopVertErrorBand(vert_errors[i]);
            delete non_genie;
        }
    }
    
    std::vector<std::string> lat_errors = data->GetLatErrorBandNames();
    for (unsigned int i = 0; i < lat_errors.size(); ++i){
            MnvLatErrorBand* non_genie = data->PopLatErrorBand(lat_errors[i]);
            delete non_genie;
    }

    std::vector<std::string> genie_errors = data->GetVertErrorBandNames();
    TH1D* err_total = new TH1D(data->GetTotalError(false,true));
    double err_genie_total = err_total->Integral();
    delete err_total;

    // Now Remove GENIE Error Bands one by one and plot them without grouping
    int n_errors = genie_errors.size();
    for (int i = n_errors-1; i >= 0; --i){
        std::string temp_name = data_var + "_" + std::to_string((long long int)i) + "_" + genie_errors[i]; 
        Systematics_DrawErrorSummary_GENIE(data, temp_name, genie_errors[i], err_genie_total);
        MnvVertErrorBand* genie = data->PopVertErrorBand(genie_errors[i]);
        delete genie;        
    }

    delete data;
}

void CCProtonPi0_Plotter::Systematics_DrawErrorSummary_GENIE(MnvH1D* hist, std::string var_name, std::string error_name, double err_genie_total)
{
    std::string plotDir = Folder_List::plotDir_Systematics_GENIE;
    
    MnvPlotter* plotter = new MnvPlotter();
    TCanvas* c = new TCanvas("c","c",1280,800);

    // Do not Group GENIE Errors
    ApplyStyle_Errors(plotter, false);

    plotter->DrawErrorSummary(hist,"TR", false);

    // Add Plot Labels
    plotter->AddHistoTitle(hist->GetTitle());

    // Current Total 
    TH1D* h_err_total = new TH1D(hist->GetTotalError(false,true));

    // Single Fractional
    std::vector<std::string> errGroup;
    errGroup.push_back(error_name);
    TH1D* h_err_single = GetTotalErrorInGroup(hist, errGroup); 
    double err_total = h_err_total->Integral();
    double err_single = h_err_single->Integral();

    delete h_err_total;
    delete h_err_single;
 
    err_genie_total *= err_genie_total;
    err_total *= err_total;
    err_single *= err_single;  

    // Add Text
    TLatex text;
    text.SetNDC();
    text.SetTextSize(0.03);
    text.SetTextColor(kBlue);
    text.DrawLatex(0.17,0.87,Form("%s = %.2f%%", "GENIE Total ",  err_genie_total/err_genie_total*100));
    text.DrawLatex(0.17,0.84,Form("%s = %.2f%%", "Current Total ",  err_total/err_genie_total*100));
    text.DrawLatex(0.17,0.81,Form("%s = %.2f%%", error_name.c_str(),  err_single/err_genie_total*100));

    // Print Plot
    std::string out_name = plotDir + var_name + ".png";
    c->Print(out_name.c_str(), "png");

    delete c;
    delete plotter;
}

void CCProtonPi0_Plotter::Systematics_WriteTables(std::string var_name)
{
    std::string plotDir = Folder_List::plotDir_OtherStudies;
    std::string root_dir = "/minerva/data/users/oaltinok/NTupleAnalysis_Final_XSecs/Data/Analyzed/CrossSection.root";
    TFile* f_data = new TFile(rootDir_CrossSection.data.c_str());

    std::string var;
    MnvH1D* hist = NULL;

    // Original
    var = var_name + "_" + "all";
    hist = GetMnvH1D(f_data, var);
    Systematics_WriteTable_Fraction(hist, var);
    Systematics_WriteTable_BinByBin(hist, var);
    delete hist;
 
    // Background Subtracted 
    var = var_name + "_" + "bckg_subtracted";
    hist = GetMnvH1D(f_data, var);
    Systematics_WriteTable_Fraction(hist, var);
    Systematics_WriteTable_BinByBin(hist, var);
    delete hist;
 
    // Unfolded 
    var = var_name + "_" + "unfolded";
    hist = GetMnvH1D(f_data, var);
    Systematics_WriteTable_Fraction(hist, var);
    Systematics_WriteTable_BinByBin(hist, var);
    delete hist;
 
    // Efficiency Corrected 
    var = var_name + "_" + "efficiency_corrected";
    hist = GetMnvH1D(f_data, var);
    Systematics_WriteTable_Fraction(hist, var);
    Systematics_WriteTable_BinByBin(hist, var);
    delete hist;
 
    // Flux Integrated 
    var = var_name + "_" + "flux_integrated";
    hist = GetMnvH1D(f_data, var);
    Systematics_WriteTable_Fraction(hist, var);
    Systematics_WriteTable_BinByBin(hist, var);
    delete hist;
 
    // Cross Section 
    var = var_name + "_" + "xsec";
    hist = GetMnvH1D(f_data, var);
    Systematics_WriteTable_Fraction(hist, var);
    Systematics_WriteTable_BinByBin(hist, var);
    delete hist;
    
    delete f_data;
}

void CCProtonPi0_Plotter::Systematics_WriteTable_Fraction(MnvH1D* hist, std::string var_name)
{
    // Open Text File
    std::string file_name = Folder_List::plotDir_Systematics_Summary + "Table_Fractional_" + var_name + ".txt";
    ofstream file;
    OpenTextFile(file_name, file);

    // Write Header
    file<<std::left;
    file.width(36); file<<"ErrorName"<<" "; 
    file.width(12); file<<"Fraction";    
    file<<std::endl;

    std::vector<std::string> errGroup;
    // Loop Over Vertical Errors
    std::vector<std::string> vert_errors = hist->GetVertErrorBandNames();
    for (unsigned int i = 0; i < vert_errors.size(); ++i){
        errGroup.clear();
        errGroup.push_back(vert_errors[i]);
        TH1D* h_err_single = GetTotalErrorInGroup(hist, errGroup); 
        int nBins = h_err_single->GetNbinsX();
        double err_single = h_err_single->Integral()/nBins;
        
        file.width(36); file<<vert_errors[i]<<" ";
        file.width(12); file<<err_single;
        file<<std::endl;
        delete h_err_single;
    }

    // Loop Over Lateral Errors
    std::vector<std::string> lat_errors = hist->GetLatErrorBandNames();
    for (unsigned int i = 0; i < lat_errors.size(); ++i){
        errGroup.clear();
        errGroup.push_back(lat_errors[i]);
        TH1D* h_err_single = GetTotalErrorInGroup(hist, errGroup); 
        int nBins = h_err_single->GetNbinsX();
        double err_single = h_err_single->Integral()/nBins;
        
        file.width(36); file<<lat_errors[i]<<" ";
        file.width(12); file<<err_single<<" ";
        file<<std::endl;
        delete h_err_single;
    }

    file.close();
}

void CCProtonPi0_Plotter::Systematics_WriteTable_BinByBin(MnvH1D* hist, std::string var_name)
{
    bool forMATLAB = false;

    // Open Text File
    std::string label = forMATLAB ? "MATLAB" : "Excel";
    std::string file_name = Folder_List::plotDir_Systematics_Summary + "Table_BinByBin_" + var_name + "_" + label + ".txt";
    ofstream file;
    OpenTextFile(file_name, file);

    if (!forMATLAB){
        // Write Header
        file<<std::left;
        file.width(16); file<<"BinRange"<<" ";    
        file.width(12); file<<"(I)Detector"<<" "; 
        file.width(12); file<<"(II)GENIE"<<" "; 
        file.width(12); file<<"(III)FSI"<<" "; 
        file.width(12); file<<"(IV)Flux"<<" "; 
        file.width(12); file<<"(V)Other"<<" "; 
        file.width(12); file<<"Total"; 
        file<<std::endl;
    }

    // Get Errors
    TH1D* h_err_total = new TH1D(hist->GetTotalError(false,true));
    TH1D* h_err_det = GetTotalErrorInGroup(hist, detGroup); // Returns new TH1D -- delete at the end 
    TH1D* h_err_genie = GetTotalErrorInGroup(hist, genieGroup); 
    TH1D* h_err_fsi = GetTotalErrorInGroup(hist, fsiGroup); 
    TH1D* h_err_flux = GetTotalErrorInGroup(hist, fluxGroup); 
    TH1D* h_err_other = GetTotalErrorInGroup(hist, otherGroup); 

    int nBins = h_err_total->GetNbinsX();
    for (int i = 1; i <= nBins; ++i){
        double bin_min = h_err_total->GetBinLowEdge(i);
        double bin_width = h_err_total->GetBinWidth(i);

        double detector = h_err_det->GetBinContent(i);
        double genie = h_err_genie->GetBinContent(i);
        double fsi = h_err_fsi->GetBinContent(i);
        double flux = h_err_flux->GetBinContent(i);
        double other = h_err_other->GetBinContent(i);
        double total = h_err_total->GetBinContent(i);

        if (!forMATLAB){
            file.width(16); file<<Form("%2.1f - %2.1f",bin_min,bin_min+bin_width)<<" ";    
        }
        file.width(12); file<<Form("%3.2f",detector*100)<<" ";
        file.width(12); file<<Form("%3.2f",genie*100)<<" ";
        file.width(12); file<<Form("%3.2f",fsi*100)<<" ";
        file.width(12); file<<Form("%3.2f",flux*100)<<" ";
        file.width(12); file<<Form("%3.2f",other*100)<<" ";
        file.width(12); file<<Form("%3.2f",total*100);
        file<<std::endl;
    }
    
    file.close();

    delete h_err_total;
    delete h_err_det;
    delete h_err_genie;
    delete h_err_fsi;
    delete h_err_flux;
    delete h_err_other;
}

void CCProtonPi0_Plotter::Systematics_SetErrorSummaryGroups()
{
    //---------------------------------------------------------------------
    // (I) Detector Response
    //---------------------------------------------------------------------
    detGroup.push_back("NeutronResponse");
    detGroup.push_back("PionResponse");
    detGroup.push_back("MuonTracking");
    detGroup.push_back("EM_EnergyScale");
    detGroup.push_back("MuonMomentum");
    detGroup.push_back("MuonTheta");
    detGroup.push_back("MichelTrue");
    detGroup.push_back("MichelFake");
    detGroup.push_back("TargetMass");
    detGroup.push_back("ProtonTracking");
    detGroup.push_back("ProtonEnergy_Birks");
    detGroup.push_back("ProtonEnergy_BetheBloch");
    detGroup.push_back("ProtonEnergy_MassModel");
    detGroup.push_back("ProtonEnergy_MEU");
   
    //---------------------------------------------------------------------
    // (II) GENIE Cross Section
    //---------------------------------------------------------------------
    genieGroup.push_back("GENIE_AhtBY"             );
    genieGroup.push_back("GENIE_BhtBY"             );
    genieGroup.push_back("GENIE_CCQEPauliSupViaKF" );
    genieGroup.push_back("GENIE_CV1uBY"            );
    genieGroup.push_back("GENIE_CV2uBY"            );
    genieGroup.push_back("GENIE_EtaNCEL"           );
    genieGroup.push_back("GENIE_MaCCQE"            );
    //genieGroup.push_back("GENIE_MaCCQEshape"       );
    genieGroup.push_back("GENIE_MaNCEL"            );
    genieGroup.push_back("GENIE_MaRES"             );
    genieGroup.push_back("GENIE_MvRES"             );
    //genieGroup.push_back("GENIE_NormCCQE"          );
    genieGroup.push_back("GENIE_NormCCRES"         );
    genieGroup.push_back("GENIE_NormDISCC"         );
    genieGroup.push_back("GENIE_NormNCRES"         );
    genieGroup.push_back("GENIE_Rvn1pi"            );
    genieGroup.push_back("GENIE_Rvn2pi"            );
    genieGroup.push_back("GENIE_Rvp1pi"            );
    genieGroup.push_back("GENIE_Rvp2pi"            );
    genieGroup.push_back("GENIE_VecFFCCQEshape"    );

    //---------------------------------------------------------------------
    // (III) GENIE FSI
    //---------------------------------------------------------------------
    fsiGroup.push_back("GENIE_AGKYxF1pi"         );
    fsiGroup.push_back("GENIE_FrAbs_N"           );
    fsiGroup.push_back("GENIE_FrAbs_pi"          );
    fsiGroup.push_back("GENIE_FrCEx_N"           );
    fsiGroup.push_back("GENIE_FrCEx_pi"          );
    fsiGroup.push_back("GENIE_FrElas_N"          );
    fsiGroup.push_back("GENIE_FrElas_pi"         );
    fsiGroup.push_back("GENIE_FrInel_N"          );
    fsiGroup.push_back("GENIE_FrInel_pi"         );
    fsiGroup.push_back("GENIE_FrPiProd_N"        );
    fsiGroup.push_back("GENIE_FrPiProd_pi"       );
    fsiGroup.push_back("GENIE_MFP_N"             );
    fsiGroup.push_back("GENIE_MFP_pi"            );
    fsiGroup.push_back("GENIE_RDecBR1gamma"      );
    fsiGroup.push_back("GENIE_Theta_Delta2Npi"   );
     
    //---------------------------------------------------------------------
    // (IV) Flux 
    //---------------------------------------------------------------------
    fluxGroup.push_back("Flux");

    //---------------------------------------------------------------------
    // (V) Other
    //---------------------------------------------------------------------
    otherGroup.push_back("BckgConstraint_WithPi0");
    otherGroup.push_back("BckgConstraint_QELike");
    otherGroup.push_back("BckgConstraint_SinglePiPlus");
    otherGroup.push_back("Unfolding");
    otherGroup.push_back("2p2h");
}

bool CCProtonPi0_Plotter::IsErrorInGroup(std::string err_name, std::vector<std::string> errGroup)
{
    for (unsigned int i = 0; i < errGroup.size(); ++i){
        if (err_name.compare(errGroup[i]) == 0) return true;
    }

    return false;
}

TH1D* CCProtonPi0_Plotter::GetTotalErrorInGroup(MnvH1D* hist, std::vector<std::string> errGroup, bool area_normalized)
{
    // Get a Copy of the Input Histogram
    MnvH1D* tempHist = new MnvH1D(*hist);

    // Remove Vertical Error Bands not in the Group
    std::vector<std::string> vert_errors = tempHist->GetVertErrorBandNames();
    
    for (unsigned int i = 0; i < vert_errors.size(); ++i){
        if (!IsErrorInGroup(vert_errors[i], errGroup)){ 
            MnvVertErrorBand* non_group = tempHist->PopVertErrorBand(vert_errors[i]);
            delete non_group;
        }
    }

    // Remove Lateral Error Bands not in the Group
    std::vector<std::string> lat_errors = tempHist->GetLatErrorBandNames();
    
    for (unsigned int i = 0; i < lat_errors.size(); ++i){
        if (!IsErrorInGroup(lat_errors[i], errGroup)){ 
            MnvLatErrorBand* non_group = tempHist->PopLatErrorBand(lat_errors[i]);
            delete non_group;
        }
    }

    // Debugging: Print What is left?
//    std::vector<std::string> all_errors = tempHist->GetErrorBandNames();
//    std::cout<<"What is Left?"<<std::endl;
//    for (unsigned int i = 0; i < all_errors.size(); ++i){
//        std::cout<<"\t"<<all_errors[i]<<std::endl;
//    }

    // Get Total Error Histogram
    TH1D* total_error = new TH1D(tempHist->GetTotalError(false,true,area_normalized));

    delete tempHist;
    
    return total_error;
}

void CCProtonPi0_Plotter::Clear_ErrorSummaryGroups()
{
    detGroup.clear();
    genieGroup.clear();
    fsiGroup.clear();
    fluxGroup.clear();
    otherGroup.clear();
}
