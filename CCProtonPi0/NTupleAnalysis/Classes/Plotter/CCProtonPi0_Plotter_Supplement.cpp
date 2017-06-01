#ifndef CCProtonPi0_Plotter_Supplement_cpp
#define CCProtonPi0_Plotter_Supplement_cpp

#include "CCProtonPi0_Plotter.h"

using namespace PlotUtils;

void CCProtonPi0_Plotter::Supplement_Tables()
{
    Supplement_XSec();
    Supplement_Errors();
    Supplement_Flux();
    Supplement_Correlation();
}

void CCProtonPi0_Plotter::Supplement_Correlation()
{
    std::string file_name = Folder_List::plotDir_Paper + "Correlation_Tables.txt";
    ofstream file;
    OpenTextFile(file_name, file);

    Supplement_Correlation("muon_P_xsec", file, "%2.1f - %2.1f", "Bins (GeV/c)");
    Supplement_Correlation("muon_theta_xsec", file, "%2.0f - %2.0f", "Bins (deg)");
    Supplement_Correlation("pi0_KE_xsec", file, "%2.2f - %2.2f", "Bins (GeV)");
    Supplement_Correlation("pi0_theta_xsec", file, "%2.0f - %2.0f", "Bins (deg)");
    Supplement_Correlation("Enu_xsec", file, "%2.1f - %2.1f", "Bins (GeV)");
    Supplement_Correlation("QSq_xsec", file, "%2.2f - %2.2f", "Bins (GeV$^{2}$)");
    Supplement_Correlation("W_xsec", file, "%2.1f - %2.1f", "Bins (GeV)");
    Supplement_Correlation("deltaInvMass_xsec", file, "%2.2f - %2.2f", "Bins (GeV)");
    Supplement_Correlation("Delta_pi_theta_xsec", file, "%2.1f - %2.1f", "Bins (Bins)");
    Supplement_Correlation("Delta_pi_phi_xsec", file, "%2.0f - %2.0f", "Bins (deg)");
    file.close();
}

void CCProtonPi0_Plotter::Supplement_Flux()
{
    std::string file_name = Folder_List::plotDir_Paper + "Flux_Table.txt";
    ofstream file;
    OpenTextFile(file_name, file);

    TFile* f_xsec_mc = new TFile(rootDir_CrossSection.mc.c_str());
    MnvH1D* hist = GetMnvH1D(f_xsec_mc, "h_flux_minervaLE_FHC");

    int nBins = hist->GetNbinsX();
    for (int i = 1; i <= nBins; ++i){
        double bin_min = hist->GetBinLowEdge(i);
        double bin_width = hist->GetBinWidth(i);
        double flux = hist->GetBinContent(i);
        if (bin_min >= 20.0) break;
       
        file.precision(2);
        file.width(14); file<<Form("%3.1f - %3.1f",bin_min,bin_min+bin_width)<<" & ";   
        file<<std::scientific; file.width(6); file<<flux<<" \\\\";
        file<<std::endl;
    }

    file.close();
}

void CCProtonPi0_Plotter::Supplement_XSec()
{
    std::string file_name = Folder_List::plotDir_Paper + "XSec_Tables.txt";
    ofstream file;
    OpenTextFile(file_name, file);

    Supplement_XSec("muon_P_xsec", file, "%2.1f - %2.1f");
    Supplement_XSec("muon_theta_xsec", file, "%2.0f - %2.0f");
    Supplement_XSec("pi0_KE_xsec", file, "%2.2f - %2.2f");
    Supplement_XSec("pi0_theta_xsec", file, "%2.0f - %2.0f");
    Supplement_XSec("Enu_xsec", file, "%2.1f - %2.1f");
    Supplement_XSec("QSq_xsec", file, "%2.2f - %2.2f");
    Supplement_XSec("W_xsec", file, "%2.1f - %2.1f");
    Supplement_XSec("deltaInvMass_xsec", file, "%2.2f - %2.2f");
    Supplement_XSec("Delta_pi_theta_xsec", file, "%2.1f - %2.1f");
    Supplement_XSec("Delta_pi_phi_xsec", file, "%2.0f - %2.0f");
    file.close();
}

void CCProtonPi0_Plotter::Supplement_Errors()
{
    std::string file_name = Folder_List::plotDir_Paper + "Systematics_Tables.txt";
    ofstream file;
    OpenTextFile(file_name, file);
    
    Supplement_Errors("muon_P_xsec", file, "%2.1f - %2.1f");
    Supplement_Errors("muon_theta_xsec", file, "%2.0f - %2.0f");
    Supplement_Errors("pi0_KE_xsec", file, "%2.2f - %2.2f");
    Supplement_Errors("pi0_theta_xsec", file, "%2.0f - %2.0f");
    Supplement_Errors("Enu_xsec", file, "%2.1f - %2.1f");
    Supplement_Errors("QSq_xsec", file, "%2.2f - %2.2f");
    Supplement_Errors("W_xsec", file, "%2.1f - %2.1f");
    Supplement_Errors("deltaInvMass_xsec", file, "%2.2f - %2.2f");
    Supplement_Errors("Delta_pi_theta_xsec", file, "%2.1f - %2.1f");
    Supplement_Errors("Delta_pi_phi_xsec", file, "%2.0f - %2.0f");
    file.close();
}

void CCProtonPi0_Plotter::Supplement_Correlation(std::string var_name, std::ofstream& file, std::string bin_format, std::string bin_name)
{
    TFile* f_xsec_data = new TFile(rootDir_CrossSection.data.c_str());

    MnvH1D* h = GetMnvH1D(f_xsec_data, var_name);
    MnvH1D* hist = new MnvH1D(h->GetBinNormalizedCopy());
    TMatrixD matrix = hist->GetTotalCorrelationMatrix();

    file<<var_name<<std::endl;
    int nBins = hist->GetNbinsX();
    // Print First Line
    file.width(14); file<<bin_name<<" & ";
    for(int i = 1; i <= nBins; ++i){
        double bin_min = hist->GetBinLowEdge(i);
        double bin_width = hist->GetBinWidth(i);

        // Handle Last Column    
        if (i == nBins){
            file.width(12); file<<Form(bin_format.c_str(),bin_min,bin_min+bin_width)<<" \\\\ "<<std::endl;;   
        }else{
            file.width(12); file<<Form(bin_format.c_str(),bin_min,bin_min+bin_width)<<" & ";   
        }
    }
    file<<"\\hline"<<std::endl;

    // Print Table
    for (int row = 1; row <= nBins; ++row){
        double bin_min = hist->GetBinLowEdge(row);
        double bin_width = hist->GetBinWidth(row);
        file.width(14); file<<Form(bin_format.c_str(),bin_min,bin_min+bin_width)<<" & ";   
        TMatrixDRow current_row = matrix[row];
        for (int col = 1; col <= nBins; ++col){
            double val = current_row[col];
            if ( col == nBins){
                file.width(12); file<<Form("%6.4f", val)<<" \\\\ "<<std::endl;
            }else{
                file.width(12); file<<Form("%6.4f", val)<<" & ";
            }
        }
    }
    file<<std::endl;
    file<<"---------------------------------------------"<<std::endl;
    file<<std::endl;
   
    delete h;
    delete hist;
    delete f_xsec_data;
}


void CCProtonPi0_Plotter::Supplement_XSec(std::string var_name, std::ofstream& file, std::string bin_format)
{
    TFile* f_xsec_data = new TFile(rootDir_CrossSection.data.c_str());

    MnvH1D* h = GetMnvH1D(f_xsec_data, var_name);
    MnvH1D* hist = new MnvH1D(h->GetBinNormalizedCopy());
    TH1D* h_err_stat = new TH1D(hist->GetStatError(true));
    TH1D* h_err_syst = new TH1D(hist->GetTotalError(false,true));
    TH1D* h_err_total = new TH1D(hist->GetTotalError(true,true));

    file<<var_name<<std::endl;
    int nBins = hist->GetNbinsX();
    for (int i = 1; i <= nBins; ++i){
        double bin_min = hist->GetBinLowEdge(i);
        double bin_width = hist->GetBinWidth(i);
        double xsec = hist->GetBinContent(i);
        double err_stat = h_err_stat->GetBinContent(i);
        double err_syst  = h_err_syst->GetBinContent(i);
        double err_total = h_err_total->GetBinContent(i);

        file.width(14); file<<Form(bin_format.c_str(),bin_min,bin_min+bin_width)<<" & ";   
        file.width(6); file<<Form("%3.2f", xsec)<<" & ";
        file.width(4); file<<Form("%3.0f", err_stat*100)<<" & ";
        file.width(4); file<<Form("%3.0f", err_syst*100)<<" & ";
        file.width(4); file<<Form("%3.0f", err_total*100)<<" \\\\";
        file<<std::endl;
    }
    file<<std::endl;
    file<<"---------------------------------------------"<<std::endl;
    file<<std::endl;
   
    delete h;
    delete hist;
    delete h_err_stat;
    delete h_err_syst;
    delete h_err_total;
    delete f_xsec_data;
}

void CCProtonPi0_Plotter::Supplement_Errors(std::string var_name, std::ofstream& file, std::string bin_format)
{
    TFile* f_xsec_data = new TFile(rootDir_CrossSection.data.c_str());

    MnvH1D* h = GetMnvH1D(f_xsec_data, var_name);
    MnvH1D* hist = new MnvH1D(h->GetBinNormalizedCopy());

    // Get Errors
    TH1D* h_err_total = new TH1D(hist->GetTotalError(false,true));
    TH1D* h_err_det = GetTotalErrorInGroup(hist, detGroup); // Returns new TH1D -- delete at the end 
    TH1D* h_err_genie = GetTotalErrorInGroup(hist, genieGroup); 
    TH1D* h_err_fsi = GetTotalErrorInGroup(hist, fsiGroup); 
    TH1D* h_err_flux = GetTotalErrorInGroup(hist, fluxGroup); 
    TH1D* h_err_other = GetTotalErrorInGroup(hist, otherGroup); 

    file<<var_name<<std::endl;
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

        file.width(14); file<<Form(bin_format.c_str(),bin_min,bin_min+bin_width)<<" & ";   
        file.width(4); file<<Form("%2.1f",detector*100)<<" & ";
        file.width(4); file<<Form("%2.1f",genie*100)<<" & ";
        file.width(4); file<<Form("%2.1f",fsi*100)<<" & ";
        file.width(4); file<<Form("%2.1f",flux*100)<<" & ";
        file.width(4); file<<Form("%2.1f",other*100)<<" & ";
        file.width(4); file<<Form("%2.1f",total*100)<<" \\\\";
        file<<std::endl;
    }
    file<<std::endl;
    file<<"---------------------------------------------"<<std::endl;
    file<<std::endl;
   
    delete h_err_total;
    delete h_err_det;
    delete h_err_genie;
    delete h_err_fsi;
    delete h_err_flux;
    delete h_err_other;
    delete hist;
    delete h;
    delete f_xsec_data;

}

#endif

