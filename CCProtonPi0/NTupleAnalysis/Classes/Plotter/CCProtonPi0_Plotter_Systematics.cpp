#include "CCProtonPi0_Plotter.h"

using namespace PlotUtils;

void CCProtonPi0_Plotter::Systematics_RawData()
{
    std::cout<<"Plotting Error Summary for Raw Data..."<<std::endl;

    Systematics_DrawErrorSummary("muon_P_all", "muon_P_mc_reco_all");
    Systematics_DrawErrorSummary("muon_theta_all", "muon_theta_mc_reco_all");
    Systematics_DrawErrorSummary("pi0_P_all", "pi0_P_mc_reco_all");
    Systematics_DrawErrorSummary("pi0_KE_all", "pi0_KE_mc_reco_all");
    Systematics_DrawErrorSummary("pi0_theta_all", "pi0_theta_mc_reco_all");
    Systematics_DrawErrorSummary("QSq_all", "QSq_mc_reco_all");
    Systematics_DrawErrorSummary("W_all", "W_mc_reco_all");
    Systematics_DrawErrorSummary("Enu_all", "Enu_mc_reco_all");

    Systematics_DrawErrorBand_GENIE("muon_P_mc_reco_all");
//    Systematics_DrawErrorBand_GENIE("muon_theta_mc_reco_all");
//    Systematics_DrawErrorBand_GENIE("pi0_P_mc_reco_all");
//    Systematics_DrawErrorBand_GENIE("pi0_KE_mc_reco_all");
//    Systematics_DrawErrorBand_GENIE("pi0_theta_mc_reco_all");
//    Systematics_DrawErrorBand_GENIE("QSq_mc_reco_all");
//    Systematics_DrawErrorBand_GENIE("W_mc_reco_all");
//    Systematics_DrawErrorBand_GENIE("Enu_mc_reco_all");

    std::cout<<"Plotting Error Summary for Raw Data Finished!"<<std::endl;
}

void CCProtonPi0_Plotter::Systematics_DrawErrorSummary(std::string data_var, std::string mc_var)
{
    std::string plotDir = Folder_List::plotDir_Systematics_Summary;
    TFile* f_xsec_data = new TFile(rootDir_CrossSection.data.c_str());
    TFile* f_xsec_mc = new TFile(rootDir_CrossSection.mc.c_str());

    MnvH1D* data = GetMnvH1D(f_xsec_data, data_var);
    MnvH1D* mc = GetMnvH1D(f_xsec_mc, mc_var);

    data_var = "data_" + data_var;
    mc_var = "mc_" + mc_var;

    DrawErrorSummary(data, data_var, plotDir);
    DrawErrorSummary(mc, mc_var, plotDir);

    Systematics_WriteTable(mc, mc_var);
    
    delete data;
    delete mc;
}

void CCProtonPi0_Plotter::Systematics_DrawErrorBand_GENIE(std::string mc_var)
{
    TFile* f_xsec_mc = new TFile(rootDir_CrossSection.mc.c_str());

    MnvH1D* mc = GetMnvH1D(f_xsec_mc, mc_var);

    std::vector<std::string> all_errors = mc->GetVertErrorBandNames();
    
    // First Remove Non-GENIE Error Bands
    for (unsigned int i = 0; i < all_errors.size(); ++i){
        std::size_t found = all_errors[i].find("GENIE");
        if (found == std::string::npos){
            MnvVertErrorBand* non_genie = mc->PopVertErrorBand(all_errors[i]);
            delete non_genie;
        }
    }

    
    std::vector<std::string> genie_errors = mc->GetVertErrorBandNames();
    TH1D* err_total = new TH1D(mc->GetTotalError(false));
    double err_genie_total = err_total->Integral();
    delete err_total;

    // Now Remove GENIE Error Bands one by one and plot them without grouping
    int n_errors = genie_errors.size();
    for (int i = n_errors-1; i >= 0; --i){
        std::string temp_name = mc_var + "_" + std::to_string((long long int)i) + "_" + genie_errors[i]; 
        Systematics_DrawErrorSummary_GENIE(mc, temp_name, genie_errors[i], err_genie_total);
        MnvVertErrorBand* genie = mc->PopVertErrorBand(genie_errors[i]);
        delete genie;        
    }

    delete mc;
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

    TH1D* h_err_total = new TH1D(hist->GetTotalError(false));
    MnvVertErrorBand* error_band = hist->GetVertErrorBand(error_name);
    TH1D* h_err_single = new TH1D(error_band->GetErrorBand()); 
    double err_total = h_err_total->Integral();
    double err_single = h_err_single->Integral();

    err_genie_total = err_genie_total*err_genie_total;
    err_total = err_total * err_total;
    err_single = err_single * err_single;
    delete h_err_total;
    delete h_err_single;
   
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

void CCProtonPi0_Plotter::Systematics_WriteTable(MnvH1D* hist, std::string var_name)
{
    // Open Text File
    std::string file_name = Folder_List::plotDir_Systematics_Summary + "Table_" + var_name + ".txt";
    ofstream file;
    OpenTextFile(file_name, file);

    // Write Header
    file<<std::left;
    file.width(36); file<<"Error Name"<<" "; 
    file.width(12); file<<"Percent";    
    file<<std::endl;

    // Get Err Total
    TH1D* h_err_total = new TH1D(hist->GetTotalError(false));
    double err_total = h_err_total->Integral();
    err_total *= err_total; // Use Square of the Error
    delete h_err_total;

    // Loop Over All Errors
    std::vector<std::string> all_errors = hist->GetVertErrorBandNames();
    double total_check = 0.0;
    for (unsigned int i = 0; i < all_errors.size(); ++i){
        MnvVertErrorBand* error_band = hist->GetVertErrorBand(all_errors[i]);
        TH1D* h_err_single = new TH1D(error_band->GetErrorBand()); 
        double err_single = h_err_single->Integral();
        err_single *= err_single; // Use Square of the Error
        total_check += err_single;
        double percent = err_single/err_total*100.0;
        
        file.width(36); file<<all_errors[i]<<" ";
        file.width(12); file<<percent<<std::endl;
        delete h_err_single;
    }

    // Write Total
    double total_check_percent = total_check/err_total*100.0;
    file.width(36); file<<"Total"<<" ";
    file.width(12); file<<total_check_percent<<std::endl;

    file.close();
}

