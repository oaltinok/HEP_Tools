#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <map>
#include <string>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <sys/stat.h>
#include <getopt.h>
#include <memory>

#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <Cintex/Cintex.h>
#include <TMath.h>
#include <TStyle.h>
#include <TMatrixD.h>
#include <TColor.h>
#include <TPaveText.h>

#include <PlotUtils/MnvH1D.h>
#include <PlotUtils/MnvH2D.h>
#include <PlotUtils/MnvPlotter.h>

#include "calc_flux.h"
#include "calc_xs.h"
#include "plot_attr.h"
#include "error_group.h"
#include "error_color.h"
#include "add_error_band.h"
#include "draw_data_mcvariation.h"
#include "draw_data_mcvariation_ratio.h"
#include "draw_data_datavariation.h"
#include "draw_data_mcbreakdown.h"
#include "draw_data_mcstacked.h"
#include "data_correction.h"
#include "print_error_group.h"
#include "draw_systematics.h"
#include "get_h1d.h"
#include "remove_bins.h"

using std::setw;
using std::sqrt;

using PlotUtils::MnvH1D;
using PlotUtils::MnvVertErrorBand;
using PlotUtils::MnvH2D;
using PlotUtils::MnvPlotter;
using MinervaUnfold::MnvUnfold;



namespace {

        // Must be called AFTER SetRootEnv to override its behaviors
    void customize_plotter(MnvPlotter* plotter)
    {
        gStyle->SetCanvasDefW(900);
        gStyle->SetCanvasDefH(675); // 4x3 aspect ratio
        gStyle->SetPadRightMargin(0.05);

        gStyle->SetStripDecimals(false);
        plotter->legend_text_size  = 0.045;
        plotter->legend_text_font = 42; // default 62 (bold)
        plotter->data_marker_size = 1.3;
        plotter->axis_title_font_x   = 42;
        plotter->axis_title_size_x   = 0.06;
        plotter->axis_title_offset_x = 1.1;
        plotter->axis_title_font_y   = 42;
        plotter->axis_title_size_y   = 0.06;
        plotter->axis_title_offset_y = 1.1;
        plotter->axis_label_size = 0.05;
        plotter->axis_label_font = 42;

            //gStyle->SetAxisColor(kGray + 1,  "XY");
            //gStyle->SetLabelColor(kGray + 1, "XY");
            //gStyle->SetFrameLineColor(kGray + 1);
        
    }

        // Must be called AFTER SetRootEnv to override its behaviors
    void customize_plotter_carrie(MnvPlotter* plotter)
    {
        gStyle->SetCanvasDefW(900);
        gStyle->SetCanvasDefH(750);
        gStyle->SetPadRightMargin(0.15);

        gStyle->SetStripDecimals(false);
        plotter->legend_offset_y = -0.05;
        plotter->legend_offset_x = +0.05;
    }

    void customize_plotter_uncertainty(MnvPlotter* plotter)
    {
        plotter->height_nspaces_per_hist = 1.3;
        plotter->width_xspace_per_letter = 0.4;
        plotter->legend_text_size        = 0.035;
    }

    
    void write_hist_title(MnvPlotter* plotter)
    {
            //const char* plot_title = "#scale[0.9]{#color[1]{#bar{#nu}_{#mu}Tracker #rightarrow #mu^{+}1#pi^{0}X (W < 1.8 GeV)}}";
            const char* plot_title = "";
            //const char* plot_title = "#scale[0.9]{#color[1]{#bar{#nu}_{#mu}Tracker #rightarrow #mu^{+}1#pi^{0}X (X has no mesons)}}";    //const char* plot_title = "#scale[0.9]{#color[1]{#bar{#nu}_{#mu} + CH #rightarrow #mu^{+} + 1#pi^{0} + nucleons}}";
        plotter->AddHistoTitle(plot_title,0.05);
    }

    void write_label(MnvPlotter* plotter)
    {
        TPaveText* text = new TPaveText(0.18, 0.85, 0.28, 0.90, "NDC");
        text->AddText("b)  ");
        text->SetTextColor(kBlue);
        text->SetFillColor(kWhite);
        text->SetTextFont(42);
        text->Draw("SAME");

    }
    
    void write_pot_exposure(MnvPlotter* plotter, double pot)
    {
            //plotter->WriteNorm("POT normalized", 0.32, 0.84,0.035, 0.909e+20);
            //plotter->WriteNorm("POT normalized", 0.32, 0.84,0.035, pot);
            //plotter->WriteNorm("POT Normalized", 0.29, 0.84,0.035, pot);
        TPaveText* text = new TPaveText(0.18, 0.80, 0.54, 0.90, "NDC");
        text->AddText("b)  #scale[1.15]{#bar{#nu}_{#mu} + CH #rightarrow #mu^{+} + #pi^{0} + X}");
        text->AddText("#scale[0.85]{POT Normalized}");
        text->SetTextColor(kBlue);
        text->SetFillColor(kWhite);
        text->SetTextFont(42);
        text->Draw("SAME");

            //std::cout << "Top margin:  " << gStyle->GetPadTopMargin() << std::endl;
            //std::cout << "Left margin: " << gStyle->GetPadLeftMargin() << std::endl;

    }

    void write_area_exposure(MnvPlotter* plotter, double pot)
    {
            //plotter->WriteNorm("Area normalized", 0.32, 0.84,0.035, 0.909e+20);
            //plotter->WriteNorm("Area normalized", 0.32, 0.84,0.035, pot);
            //plotter->WriteNorm("Area Normalized", 0.29, 0.84,0.035, pot);
        TPaveText* text = new TPaveText(0.18, 0.80, 0.54, 0.90, "NDC");
        text->AddText("b)  #scale[1.15]{#bar{#nu}_{#mu} + CH #rightarrow #mu^{+} + #pi^{0} + X}");
        text->AddText("#scale[0.85]{Area Normalized}");
        text->SetTextColor(kBlue);
        text->SetFillColor(kWhite);
        text->SetTextFont(42);
        text->Draw("SAME");


    }

void dump(TH1D* cv, TH1D* stat = NULL, TH1D* syst = NULL)
    {
        std::cout.setf(std::ios_base::fixed);
        std::cout.precision(3);
        std::cout << "Dumping histogram with name: " << cv->GetName() << std::endl;

        if (stat && syst) {
            std::cout << "       bin   low        high      value      tot      tot_frac     stat      stat_frac      sys      sys_frac" << std::endl;
            for (int i = 0; i <= cv->GetNbinsX() + 1; ++i) {
                double content = cv->GetBinContent(i);
                std::cout << setw(8)  << i
                          << setw(10) << cv->GetBinLowEdge(i)
                          << setw(10) << cv->GetBinLowEdge(i) + cv->GetBinWidth(i)
                          << setw(10) << cv->GetBinContent(i)
                          << setw(10) << cv->GetBinError(i)
                          << setw(10) << (content > 0 ? 100*cv->GetBinError(i)/content : 0.0) << " (%) "
                          << setw(10) << stat->GetBinContent(i)
                          << setw(10) << (content > 0 ? 100*stat->GetBinContent(i)/content : 0.0) << " (%) "
                          << setw(10) << syst->GetBinContent(i)
                          << setw(10) << (content > 0 ? 100*syst->GetBinContent(i)/content : 0.0) << " (%) "
                          << std::endl;
            }
        }


        if (stat && !syst) {
            for (int i = 0; i <= cv->GetNbinsX() + 1; ++i) {
                double content = cv->GetBinContent(i);
                std::cout << setw(8)  << i
                          << setw(10) << cv->GetBinLowEdge(i)
                          << setw(10) << cv->GetBinLowEdge(i) + cv->GetBinWidth(i)
                          << setw(10) << cv->GetBinContent(i)
                          << setw(10) << cv->GetBinError(i)
                          << setw(10) << (content > 0 ? 100*cv->GetBinError(i)/content : 0.0) << " (%) "
                          << setw(10) << stat->GetBinContent(i)
                          << setw(10) << (content > 0 ? 100*stat->GetBinContent(i)/content : 0.0) << " (%) "
                          << std::endl;
            }
        }

        if (!stat && syst) {
            for (int i = 0; i <= cv->GetNbinsX() + 1; ++i) {
                double content = cv->GetBinContent(i);
                std::cout << setw(8)  << i
                          << setw(10) << cv->GetBinLowEdge(i)
                          << setw(10) << cv->GetBinLowEdge(i) + cv->GetBinWidth(i)
                          << setw(10) << cv->GetBinContent(i)
                          << setw(10) << cv->GetBinError(i)
                          << setw(10) << (content > 0 ? 100*cv->GetBinError(i)/content : 0.0) << " (%) "
                          << setw(10) << syst->GetBinContent(i)
                          << setw(10) << (content > 0 ? 100*syst->GetBinContent(i)/content : 0.0) << " (%) "
                          << std::endl;
            }
        }

        if (!stat && !syst) {
            for (int i = 0; i <= cv->GetNbinsX() + 1; ++i) {
                double content = cv->GetBinContent(i);
                std::cout << setw(8)  << i
                          << setw(10) << cv->GetBinCenter(i)
                          << setw(10) << cv->GetBinContent(i)
                          << setw(10) << cv->GetBinError(i)
                          << setw(10) << (content > 0 ? 100*cv->GetBinError(i)/content : 0.0) << " (%) "
                          << std::endl;
            }
        }
        
    }


    void area_normalize_universes(MnvH1D* h)
    {
        double cvArea = h->Integral();
        std::vector<std::string> vertNames = h->GetVertErrorBandNames();
        for (std::vector<std::string>::iterator name = vertNames.begin();
             name != vertNames.end(); ++name) {
            MnvVertErrorBand* errorBand = h->GetVertErrorBand(*name);
            assert(errorBand);

            std::vector<TH1D*> universe_hists = errorBand->GetHists();
            assert(!universe_hists.empty());
            for (std::vector<TH1D*>::iterator hvar = universe_hists.begin();
                 hvar != universe_hists.end(); ++hvar) {
                (*hvar)->Scale(cvArea/(*hvar)->Integral());
            }
        }
    }
    

    double integral(MnvH1D* h)
    {
        double area = 0.0;
        double norm_width = h->GetNormBinWidth();
        for (int i = 0; i <= h->GetNbinsX(); ++i) {
            double content = h->GetBinContent(i);
            double width   = h->GetBinWidth(i);

            area += content * width/norm_width;
        }

        return area;
    }
    
        
    double calc_area_ratio(MnvH1D* mnvh1d_data, MnvH1D* mnvh1d_mc)
    {


        double area_data  = mnvh1d_data->Integral();
        assert(area_data > 0.0);
        
        double area_mc  = mnvh1d_mc->Integral();
        assert(area_mc > 0.0);

            //std::cout << "Title: " << mnvh1d_mc->GetTitle() << std::endl;
            //std::cout << "\tData area: " << area_data << std::endl;
            //std::cout << "\tMC area:   " << area_mc << std::endl;
        
        return area_data/area_mc;
	
	
    }

        // type = 1: not bin width normalized: for neutrino energy xs
        // type = 0: return the usual bin normalized copy
        // type = 2: reserved
    MnvH1D mnvh1d_adapter(MnvH1D* input, int type)
    {
        if (type == 1) return MnvH1D(*input);
        if (type == 0) return input->GetBinNormalizedCopy();

        assert(false);

        return input->GetBinNormalizedCopy();
    }

    
    MnvH1D* mnvh1d_adapter_ptr(MnvH1D* input, int type)
    {
        MnvH1D* result(NULL);
        if (type == 1) result = new MnvH1D(*input);
        if (type == 0) result = new MnvH1D(input->GetBinNormalizedCopy());
        
        assert(result);
        mem_mgr::get().manage(result);
        
        return result;
        
    }
    
    
    const TH2D* reformat(const TH2D* migration)
    {
        int Nx = migration->GetNbinsX();
        int Ny = migration->GetNbinsY();


        TH2D* reformatted_hist = new TH2D("migration_formatted", "Migration matrix",
                                          Nx, 0.0, (double) Nx,
                                          Ny, 0.0, (double) Ny);
        
        for (int i = 1; i <= Nx; ++i) {
            for (int j = 1; j <= Ny; ++j) {
                double content = migration->GetBinContent(i,j);

                reformatted_hist->SetBinContent(i,j,content);
            }
        }

        return reformatted_hist;
    }


    double* convert_to_array(const std::vector<double>& v)
    {
            // note the square bracket!
        double* a = new double[v.size()];
        for (std::size_t i = 0; i < v.size(); ++i) {
            a[i] = v[i];
        }
        
        return a;
    }

    
    MnvH1D* neut_to_mnvh1d(const std::string& neut_xsec, const std::string& variable)
    {
        std::cout << "Read NEUT cross section " << std::endl;
        
        MnvH1D* output(NULL);

            // add w here since I want to report the w without the W < 1.8 cut.
            // For the pion variables, there was no W cut, so I'm adding W here
            // Add fake histogram to the nuwro-cross-section.root and neut-cross-section.root so that
            // the program runs
        if (variable == "pimom" || variable == "theta" || variable == "kinetic" || variable == "w") {

            if (neut_xsec.find(".root") == std::string::npos) {
                std::cerr << "Not a root file, expect root file" << std::endl;
                exit(1);
            }
                    
            TFile* neut_file = new TFile(neut_xsec.c_str());
            TH1D*  hneut_xs  = (TH1D*) neut_file->Get(Form("neut_%s", variable.c_str()));
            
            output = new MnvH1D(*hneut_xs);

        } else {

            if (neut_xsec.find(".root") != std::string::npos) {
                std::cerr << "Got a root file, expect text file" << std::endl;
                exit(1);
                
            }
                // read text file
                // convert text file to TH1D*
                // convert TH1D* to MnvH1D*
            std::vector<double> xvec;
            std::vector<double> yvec;

            std::ifstream ifs(neut_xsec.c_str());
            while (!ifs.eof()) {
                double xval;
                double yval;

                ifs >> xval >> yval;

                if (ifs.eof()) break;

                xvec.push_back(xval);
                yvec.push_back(yval * 10);
                
            }

                // create and fill the histogram
            TH1D* hneut_xs = new TH1D(Form("neut-%s", variable.c_str()), "NEUT xsec", xvec.size()-1, convert_to_array(xvec));
            for (unsigned int i = 0; i < xvec.size(); ++i) {
                hneut_xs->SetBinContent(i+1, yvec[i]);
            }
            
            output = new MnvH1D(*hneut_xs);
            
        }
        
        return output;

    }

        // nuwro_xsec: NuWro cross section .dist file name
        // variable: variable (used only for pion momentum and angle
        // col = 1: without FSI
        // col = 2: with FSI
    MnvH1D* nuwro_to_mnvh1d(const std::string& nuwro_xsec, const std::string& variable, int col)
    {
        std::cout << "Read NuWro cross section " << std::endl;
        
        MnvH1D* output(NULL);
            // see comment in neut_to_mnvh1d
        if (variable == "pimom" || variable == "theta" || variable == "kinetic" || variable == "w") {

            if (nuwro_xsec.find(".root") == std::string::npos) {
                std::cerr << "Not a root file, expect root file" << std::endl;
                exit(1);
            }
                    
            TFile* nuwro_file = new TFile(nuwro_xsec.c_str());
            TH1D*  hnuwro_xs  = (TH1D*) nuwro_file->Get(Form("nuwro_%s", variable.c_str()));            
            
            output = new MnvH1D(*hnuwro_xs);

        } else {

            if (nuwro_xsec.find(".root") != std::string::npos) {
                std::cerr << "Got a root file, expect text file" << std::endl;
                exit(1);
                
            }
                // read text file
                // convert text file to TH1D*
                // convert TH1D* to MnvH1D*
            std::vector<double> xvec;
            std::vector<double> yvec;
            
            std::ifstream ifs(nuwro_xsec.c_str());
            while (!ifs.eof()) {
                double xval;
                double yval1;
                double yval2;
                
                ifs >> xval >> yval1 >> yval2;

                if (ifs.eof()) break;

                xvec.push_back(xval);
                if (col == 1) yvec.push_back(yval1 * 1e40);
                if (col == 2) yvec.push_back(yval2 * 1e40);
                
            }

                // create and fill the histogram
            TH1D* hnuwro_xs = new TH1D(Form("nuwro-%s", variable.c_str()), "NuWro xsec", xvec.size()-1, convert_to_array(xvec));
            for (unsigned int i = 0; i <= xvec.size(); ++i) {
                hnuwro_xs->SetBinContent(i+1, yvec[i]);
            }
            
            output = new MnvH1D(*hnuwro_xs);
            
        }
        
        return output;

    }

        // pass the second histogram for binning info
    void print_matrix_stdcout(const TMatrixD& mat, MnvH1D* h)
    {
        int N = mat.GetNrows();
        for (int i = 1; i < N; ++i) {
            printf("%8d %12.2f %12.2f   ", i,
                   h->GetBinLowEdge(i),
                   h->GetBinLowEdge(i) + h->GetBinWidth(i));
            
            for (int j = 1; j < N; ++j) {
                printf(" %10.4f ", mat(i,j));
            }
            
            std::cout << std::endl;
        }
    }


    void print_matrix_ofstream(const TMatrixD& mat, MnvH1D* h, std::ofstream& ofs)
    {
        int N = mat.GetNrows();
        for (int i = 1; i < N-1; ++i) {
            ofs << Form("%6.2f - %6.2f   ",
                        h->GetBinLowEdge(i),
                        h->GetBinLowEdge(i) + h->GetBinWidth(i));
            
            for (int j = 1; j < N-1; ++j) {
                ofs << Form(" & %10.4f ", mat(i,j));
            }
            
            ofs << "\\\\ \n";
        }

    }


    void print_correlation_matrix_from_covariance(const TMatrixD& mat, MnvH1D* h)
    {
        int N = mat.GetNrows();
        for (int i = 1; i < N; ++i) {
            printf("%8d %12.2f %12.2f   ", i,
                   h->GetBinLowEdge(i),
                   h->GetBinLowEdge(i) + h->GetBinWidth(i));
            
            for (int j = 1; j < N; ++j) {
                double A = sqrt(mat(i,i) * mat(j,j));
                double corr = (A > 0.0 ? mat(i,j)/A : 0.0);
                printf(" %10.4f ", corr);
            }
            
            std::cout << std::endl;
        }
   
    }

}


int main(int argc, char* argv[])
{

    ROOT::Cintex::Cintex::Enable();

    TH1::AddDirectory(false);

    
    bool __reweight_flux = false;
    bool __enu_cut       = false;
    bool __carrie_style  = false;
    bool __legacy_code   = false;
    bool __draw_genie_band = false;
    
    int    __iteration = 1;
    double __nsigma   =   2.0;
    double __sigma    =  30.0;
    double __reweight_min = 0.0;
    double __reweight_max = 0.0;
    double __reweight_amount = 0.0;

    
    std::string __minerva_histogram_file;
    std::string __minerva_efficiency_file;
    std::string __minerva_response_file;
    std::string __downstream_histogram_file;
    std::string __downstream_efficiency_file;
    std::string __downstream_response_file;

    std::string __flux_file;
    std::string __xsec_file;
    std::string __variable;
    std::string __figure_format;
    std::string __neut_file;
    std::string __nuwro_file;

    std::string __export_xsec_data;
        
    for (;;) {

        int option_index = 0;
        static struct option long_options[] = {
            {"with-enu-cut",          no_argument,       0,  0 },
            {"reweight-flux",         no_argument,       0,  0 },
            {"carrie-style",          no_argument,       0,  0 },
            {"legacy-code",           no_argument,       0,  0 },
            {"draw-genie-band",       no_argument,       0,  0 },
            {"iteration",             required_argument, 0,  0 },
            {"sigma",                 required_argument, 0,  0 },
            {"nsigma",                required_argument, 0,  0 },
            {"reweight-min",          required_argument, 0,  0 },
            {"reweight-max",          required_argument, 0,  0 },
            {"reweight-amount",       required_argument, 0,  0 },
            {"minerva-histogram",     required_argument, 0,  0 },
            {"minerva-efficiency",    required_argument, 0,  0 },
            {"minerva-response",      required_argument, 0,  0 },
            {"downstream-histogram",  required_argument, 0,  0 },
            {"downstream-efficiency", required_argument, 0,  0 },
            {"downstream-response",   required_argument, 0,  0 },
            {"xsec",                  required_argument, 0,  0 },
            {"flux",                  required_argument, 0,  0 },
            {"neut",                  required_argument, 0,  0 },
            {"nuwro",                 required_argument, 0,  0 },
            {"variable",              required_argument, 0,  0 },
            {"figure-format",         required_argument, 0,  0 },
            {"export",                required_argument, 0,  0 },
            {0,         0,                 0,  0 }
        };
        
        int c = getopt_long(argc, argv, "",
                            long_options, &option_index);
        
        if (c != 0 ) break;
        
        std::string opt_name(long_options[option_index].name);

        
        if      (opt_name == "with-enu-cut")          __enu_cut         = true;
        else if (opt_name == "reweight-flux")         __reweight_flux   = true;
        else if (opt_name == "carrie-style")          __carrie_style    = true;
        else if (opt_name == "legacy-code")           __legacy_code     = true;
        else if (opt_name == "draw-genie-band")       __draw_genie_band = true;
        else if (opt_name == "iteration")             __iteration       = atoi(optarg);
        else if (opt_name == "sigma")                 __sigma           = atof(optarg);
        else if (opt_name == "nsigma")                __nsigma          = atof(optarg);
        else if (opt_name == "minerva-histogram")     __minerva_histogram_file     = optarg;
        else if (opt_name == "minerva-efficiency")    __minerva_efficiency_file    = optarg;
        else if (opt_name == "minerva-response")      __minerva_response_file      = optarg;
        else if (opt_name == "downstream-histogram")  __downstream_histogram_file  = optarg;
        else if (opt_name == "downstream-efficiency") __downstream_efficiency_file = optarg;
        else if (opt_name == "downstream-response")   __downstream_response_file   = optarg;
        else if (opt_name == "xsec")                  __xsec_file        = optarg;
        else if (opt_name == "flux")                  __flux_file        = optarg;
        else if (opt_name == "neut")                  __neut_file        = optarg;
        else if (opt_name == "nuwro")                 __nuwro_file       = optarg;
        else if (opt_name == "variable")              __variable         = optarg;
        else if (opt_name == "figure-format")         __figure_format    = optarg;
        else if (opt_name == "export")                __export_xsec_data = optarg;
        else if (opt_name == "reweight-min")          __reweight_min     = atof(optarg);
        else if (opt_name == "reweight-max")          __reweight_max     = atof(optarg);
        else if (opt_name == "reweight-amount")       __reweight_amount  = atof(optarg);
        else {}
    }
    
    std::cout << "__draw_genie_band  " << __draw_genie_band  << std::endl;
    std::cout << "__reweight_flux    " << __reweight_flux    << std::endl;
    std::cout << "__carrie_style     " << __carrie_style     << std::endl;
    std::cout << "__legacy_code      " << __legacy_code      << std::endl;
    std::cout << "__enu_cut          " << __enu_cut          << std::endl;
    std::cout << "__iteration        " << __iteration        << std::endl;
    std::cout << "__sigma            " << __sigma            << std::endl;
    std::cout << "__nsigma           " << __nsigma           << std::endl;
    std::cout << "__reweight_min     " << __reweight_min     << std::endl;
    std::cout << "__reweight_max     " << __reweight_max     << std::endl;
    std::cout << "__reweight_amount  " << __reweight_amount  << std::endl;
    std::cout << "__minerva_histogram_file     " << __minerva_histogram_file      << std::endl;
    std::cout << "__minerva_efficiency_file    " << __minerva_efficiency_file     << std::endl;
    std::cout << "__minerva_response_file      " << __minerva_response_file       << std::endl;
    std::cout << "__downstream_histogram_file  " << __downstream_histogram_file   << std::endl;
    std::cout << "__downstream_efficiency_file " << __downstream_efficiency_file  << std::endl;
    std::cout << "__downstream_response_file   " << __downstream_response_file    << std::endl;
    std::cout << "__xsec_file                  " << __xsec_file                   << std::endl;
    std::cout << "__flux_file                  " << __flux_file                   << std::endl;
    std::cout << "__neut_file                  " << __neut_file                   << std::endl;
    std::cout << "__nuwro_file                 " << __nuwro_file                  << std::endl;
    std::cout << "__variable                   " << __variable                    << std::endl;
    std::cout << "__figure_format              " << __figure_format               << std::endl;
    std::cout << "__export_xsec_data           " << __export_xsec_data            << std::endl;
    
        // check if file actually exists;
    std::map<std::string, std::string> filename_map;
    filename_map.insert(std::make_pair("__minerva_histogram_file",     __minerva_histogram_file));
    filename_map.insert(std::make_pair("__minerva_efficiency_file",    __minerva_efficiency_file));
    filename_map.insert(std::make_pair("__minerva_response_file",      __minerva_response_file));
    filename_map.insert(std::make_pair("__downstream_histogram_file",  __downstream_histogram_file));
    filename_map.insert(std::make_pair("__downstream_efficiency_file", __downstream_efficiency_file));
    filename_map.insert(std::make_pair("__downstream_response_file",   __downstream_response_file));
    filename_map.insert(std::make_pair("__xsec_file",  __xsec_file));
    filename_map.insert(std::make_pair("__flux_file",  __flux_file));
        //filename_map.insert(std::make_pair("__neut_file",  __neut_file));
        //filename_map.insert(std::make_pair("__nuwro_file", __nuwro_file));
    for (std::map<std::string, std::string>::iterator n = filename_map.begin();
         n != filename_map.end(); ++n) {
        std::string name = n->second;
        std::string args = n->first;
        if (name.empty()) {

            std::cerr << "args not set: " << args << std::endl;
            exit(1);
        }
        
        struct stat buffer;   
        if (stat(name.c_str(), &buffer) != 0) {
            std::cerr << "File not found: " << name << std::endl;
            exit(1);
        }
        
    }



    
    const char* variable_name = __variable.c_str();



    const double m0 = 134.9766;

    const double lower_peak = m0 - __nsigma * __sigma;
    const double upper_peak = m0 + __nsigma * __sigma;
   
 
    
        //------------------------------------------------------
    MnvPlotter* plotter = new MnvPlotter();
    plotter->error_summary_group_map = create_error_group_for_table();
    plotter->error_color_map = get_error_colors();
    plotter->SetRootEnv();
    if  (__carrie_style) customize_plotter_carrie(plotter);
    else customize_plotter(plotter);
    gStyle->SetEndErrorSize(4);
    
    std::cout << "MnvPlotter::axis_title_offset_y: " << plotter->axis_title_offset_y << std::endl;
    std::cout << "MnvPlotter::axis_title_size_y:   " << plotter->axis_title_size_y << std::endl;
    std::cout << "MnvPlotter::axis_label_size:     " << plotter->axis_label_size << std::endl;

    
    // const double data_pot  =  0.909e20; old processing
    const int minerva_nmodules    = 54; // default tracker fiducial volume
    const int downstream_nmodules = 27; // for default fiducial volume cut at 7200 mm
    const double vol_ratio = (1.0 * downstream_nmodules) / minerva_nmodules;
    
    const double mc_pot    = 10.000e20;
    const double data_minerva_pot  =  1.061e20;
    const double data_downstream_pot =  0.9458e20; 
    const double recorded_pot = data_minerva_pot + data_downstream_pot;
    const double equivalent_pot = data_downstream_pot * vol_ratio; // downstream exposure is equivalent to
                                                                   // the full detector with this exposure
    
    const double tot_pot = data_minerva_pot + equivalent_pot;
    const double mnv_scale_pot = data_minerva_pot/mc_pot;
    const double tot_scale_pot = tot_pot/mc_pot;

    std::cout << "Minerva 5 POT (e20)  " << data_minerva_pot/1.e20    << std::endl;
    std::cout << "Downstream POT (e20) " << data_downstream_pot/1.e20 << std::endl;
    std::cout << "Recorded POT (e20)   " << recorded_pot/1.e20      << std::endl;
    std::cout << "Equivalent POT (e20) " << equivalent_pot/1.e20      << std::endl;
    std::cout << "Total POT (e20)      " << tot_pot/1.e20             << std::endl;
    

    
    std::map<std::string, plot_attr> plot_attr_map = create_plot_attr_map();
    plot_attr plot_attr_variable = plot_attr_map[std::string(variable_name)];
    plot_attr plot_attr_mgg      = plot_attr_map["mgg"];
    double width = plot_attr_variable.binwidth;
    
    bool is_enu = __variable == "enu";

    
    TFile* minerva_hist_file       = new TFile(__minerva_histogram_file.c_str());
    TFile* minerva_response_file   = new TFile(__minerva_response_file.c_str());
    TFile* minerva_efficiency_file = new TFile(__minerva_efficiency_file.c_str());

    TFile* downstream_hist_file       = new TFile(__downstream_histogram_file.c_str());
    TFile* downstream_response_file   = new TFile(__downstream_response_file.c_str());
    TFile* downstream_efficiency_file = new TFile(__downstream_efficiency_file.c_str());

    TFile* xsec_file = new TFile(__xsec_file.c_str());
    
    assert(minerva_hist_file);
    assert(minerva_response_file);
    assert(minerva_efficiency_file);

    assert(downstream_hist_file);
    assert(downstream_response_file);
    assert(downstream_efficiency_file);
    
    

        // minerva5
    data_correction handler_minerva5(minerva_hist_file,
                                     minerva_response_file,
                                     minerva_efficiency_file,
                                     xsec_file,
                                     "minerva",
                                     __variable);

    handler_minerva5.background_subtract(lower_peak, upper_peak);
    handler_minerva5.unfold(__iteration);
    handler_minerva5.efficiency_divide();

    handler_minerva5.set_norm_bin_width(width);

    MnvH1D* mnvh1d_minerva_mcsig_raw          = handler_minerva5.__mnvh1d_mcsig_raw;
    MnvH1D* mnvh1d_minerva_mcsig_produced     = handler_minerva5.__mnvh1d_gnsig_nocut;
    MnvH1D* mnvh1d_minerva_mcsig_selected     = handler_minerva5.__mnvh1d_gnsig_allcut;
    MnvH1D* mnvh1d_minerva_efficiency         = handler_minerva5.__mnvh1d_efficiency;
    MnvH1D* mnvh1d_minerva_data_bkgsub        = handler_minerva5.__mnvh1d_data_bkgsub;
    MnvH1D* mnvh1d_minerva_data_bkg           = handler_minerva5.__mnvh1d_estimated_background;
    MnvH1D* mnvh1d_minerva_data_unfolded      = handler_minerva5.__mnvh1d_data_unfolded;
    MnvH1D* mnvh1d_minerva_data_effcorrected  = handler_minerva5.__mnvh1d_data_eff_corrected;
        
    
        // downstream
    data_correction handler_downstream(downstream_hist_file,
                                       downstream_response_file,
                                       downstream_efficiency_file,
                                       xsec_file,      // always minerva5 for true xs
                                       "downstream",
                                       __variable);
    
    handler_downstream.background_subtract(lower_peak, upper_peak);
    handler_downstream.unfold(__iteration);
    handler_downstream.efficiency_divide();
    
    handler_downstream.set_norm_bin_width(width);
    

    MnvH1D* mnvh1d_downstream_efficiency         = handler_downstream.__mnvh1d_efficiency;
    MnvH1D* mnvh1d_downstream_data_bkgsub        = handler_downstream.__mnvh1d_data_bkgsub;
    MnvH1D* mnvh1d_downstream_data_bkg           = handler_downstream.__mnvh1d_estimated_background;
    MnvH1D* mnvh1d_downstream_data_unfolded      = handler_downstream.__mnvh1d_data_unfolded;
    MnvH1D* mnvh1d_downstream_data_effcorrected  = handler_downstream.__mnvh1d_data_eff_corrected;


        // integrated flux, input histogram is used as template for multi-universe
    MnvH1D* mnvh1d_integrated_flux = calc_flux(mnvh1d_minerva_data_effcorrected,
                                               __flux_file,
                                               __enu_cut,is_enu,
                                               __reweight_flux,
                                               __reweight_min,
                                               __reweight_max,
                                               __reweight_amount);
    assert(mnvh1d_integrated_flux);
    
        // final cross section using minerva5 data
    MnvH1D* mnvh1d_minerva_data_multiflux_xs = calc_xs(mnvh1d_minerva_data_effcorrected, mnvh1d_integrated_flux,
                                                       data_minerva_pot,false);
    assert(mnvh1d_minerva_data_multiflux_xs);
    mnvh1d_minerva_data_multiflux_xs->Scale(plot_attr_variable.scale);
    mnvh1d_minerva_data_multiflux_xs->GetXaxis()->SetTitle(plot_attr_variable.xs_xlabel.c_str());
    mnvh1d_minerva_data_multiflux_xs->GetYaxis()->SetTitle(plot_attr_variable.xs_ylabel.c_str());
    
        // final cross section using downstream data
    MnvH1D* mnvh1d_downstream_data_multiflux_xs = calc_xs(mnvh1d_downstream_data_effcorrected, mnvh1d_integrated_flux,
                                                          data_downstream_pot/2.0,false); // half detector
    assert(mnvh1d_downstream_data_multiflux_xs);
    mnvh1d_downstream_data_multiflux_xs->Scale(plot_attr_variable.scale);   
    mnvh1d_downstream_data_multiflux_xs->GetXaxis()->SetTitle(plot_attr_variable.xs_xlabel.c_str());    
    mnvh1d_downstream_data_multiflux_xs->GetYaxis()->SetTitle(plot_attr_variable.xs_ylabel.c_str());
    
        // add two datasets and then calculate the cross sections
    MnvH1D* mnvh1d_tot_data_effcorrected = new MnvH1D(*mnvh1d_minerva_data_effcorrected);
    mnvh1d_tot_data_effcorrected->Add(mnvh1d_downstream_data_effcorrected);

    MnvH1D* mnvh1d_tot_data_multiflux_xs = calc_xs(mnvh1d_tot_data_effcorrected, mnvh1d_integrated_flux, tot_pot, false);
    mnvh1d_tot_data_multiflux_xs->Scale(plot_attr_variable.scale);   
    mnvh1d_tot_data_multiflux_xs->GetXaxis()->SetTitle(plot_attr_variable.xs_xlabel.c_str());
    mnvh1d_tot_data_multiflux_xs->GetXaxis()->SetNdivisions(plot_attr_variable.xndivision);

    mnvh1d_tot_data_multiflux_xs->GetXaxis()->SetLimits(plot_attr_variable.xs_xmin,
                                                        plot_attr_variable.xs_xmax);

    mnvh1d_tot_data_multiflux_xs->GetYaxis()->CenterTitle();
    mnvh1d_tot_data_multiflux_xs->GetYaxis()->SetTitle(plot_attr_variable.xs_ylabel.c_str());
    mnvh1d_tot_data_multiflux_xs->GetYaxis()->SetNdivisions(plot_attr_variable.yndivision);
    

        // pot x target weighted average of the sigma_5 and sigma_ds as suggested by Kevin
        // (sigma_5 x target5 x pot5 + sigma_ds x target_ds x pot_ds)/(target5 x pot5 + target_ds x pot_ds)
    MnvH1D* mnvh1d_avg1_data_multiflux_xs = new MnvH1D(*mnvh1d_minerva_data_multiflux_xs);
    mnvh1d_avg1_data_multiflux_xs->Scale(54.0 * data_minerva_pot/1e20);
    mnvh1d_avg1_data_multiflux_xs->Add(mnvh1d_downstream_data_multiflux_xs, 27.0 * data_downstream_pot/1e20);
    mnvh1d_avg1_data_multiflux_xs->Scale(1.0/(54.0 * data_minerva_pot/1e20 + 27.0 * data_downstream_pot/1e20));

        // pot weighted average that I came up with
        //MnvH1D* mnvh1d_avg2_data_multiflux_xs = new MnvH1D(*mnvh1d_minerva_data_multiflux_xs);
        //mnvh1d_avg2_data_multiflux_xs->Scale(data_minerva_pot/1e20);
        //mnvh1d_avg2_data_multiflux_xs->Add(mnvh1d_downstream_data_multiflux_xs, data_downstream_pot/1e20);
        //mnvh1d_avg2_data_multiflux_xs->Scale(1.0/(data_minerva_pot/1e20 + data_downstream_pot/1e20));
    MnvH1D* mnvh1d_tot_data_effcorrected2 = new MnvH1D(*mnvh1d_minerva_data_effcorrected);
    mnvh1d_tot_data_effcorrected2->Add(mnvh1d_downstream_data_effcorrected, 2);
    MnvH1D* mnvh1d_avg2_data_multiflux_xs = calc_xs(mnvh1d_tot_data_effcorrected2, mnvh1d_integrated_flux, recorded_pot, false);
    
        // true cross sections by GENIE
    MnvH1D* mnvh1d_integrated_flux_smooth = calc_flux(handler_minerva5.__mnvh1d_gnsig_xsec_smooth,
                                                      __flux_file,
                                                      __enu_cut,is_enu,
                                                      __reweight_flux,
                                                      __reweight_min,
                                                      __reweight_max,
                                                      __reweight_amount);
    
    assert(mnvh1d_integrated_flux_smooth);
    
    MnvH1D* mnvh1d_gnsig_multiflux_xs = calc_xs(handler_minerva5.__mnvh1d_gnsig_xsec_smooth,
                                                mnvh1d_integrated_flux_smooth, mc_pot, true);
    
    assert(mnvh1d_gnsig_multiflux_xs);
    mnvh1d_gnsig_multiflux_xs->Scale(plot_attr_variable.scale);  
    mnvh1d_gnsig_multiflux_xs->SetNormBinWidth(width);
    mnvh1d_gnsig_multiflux_xs->GetYaxis()->SetTitle(plot_attr_variable.xs_ylabel.c_str());
    mnvh1d_gnsig_multiflux_xs->GetXaxis()->SetTitle(plot_attr_variable.xs_xlabel.c_str());    
    //if (__variable == "pimom") mnvh1d_gnsig_multiflux_xs->Scale(plot_attr_variable.scale*2);
    
    
    MnvH1D* mnvh1d_denominator_nofsi = (MnvH1D*) minerva_efficiency_file->Get(Form("mc_cut00_primary_%s",variable_name));

    assert(mnvh1d_denominator_nofsi);
    mnvh1d_denominator_nofsi->SetNormBinWidth(width);
    mnvh1d_denominator_nofsi->GetYaxis()->CenterTitle();
    MnvH1D* mnvh1d_gnsig_multiflux_xs_nofsi = calc_xs(mnvh1d_denominator_nofsi,
                                                      mnvh1d_integrated_flux_smooth, mc_pot, true);
    assert(mnvh1d_gnsig_multiflux_xs_nofsi);
    mnvh1d_gnsig_multiflux_xs_nofsi->Scale(plot_attr_variable.scale);  
    
    const bool as_frac = true;
    const bool not_as_frac = false;
    const bool cov_area_norm = true;
    const bool not_cov_area_norm = false;
    const bool include_stat = true;
        //const bool not_include_stat = false;
  

        // print cross sections and uncertainties
    std::cout << "Final cross sections (minerva5 only): " << __variable << std::endl;
        //dump(&mnvh1d_adapter(mnvh1d_minerva_data_multiflux_xs, is_enu).GetCVHistoWithError(include_stat, not_cov_area_norm),  
        // &mnvh1d_adapter(mnvh1d_minerva_data_multiflux_xs, is_enu).GetStatError(not_as_frac),                    
        // &mnvh1d_adapter(mnvh1d_minerva_data_multiflux_xs, is_enu).GetTotalError(not_include_stat, not_as_frac, not_cov_area_norm)); 

       // print cross sections and uncertainties
    std::cout << "Final cross sections (downstream only): " << __variable << std::endl;
        //dump(&mnvh1d_adapter(mnvh1d_downstream_data_multiflux_xs, is_enu).GetCVHistoWithError(include_stat, not_cov_area_norm),
        // &mnvh1d_adapter(mnvh1d_downstream_data_multiflux_xs, is_enu).GetStatError(not_as_frac),                
        // &mnvh1d_adapter(mnvh1d_downstream_data_multiflux_xs, is_enu).GetTotalError(not_include_stat, not_as_frac, not_cov_area_norm));

    std::cout << "Final cross sections (minerva5 + downstream): " << __variable << std::endl;
        //dump(&mnvh1d_adapter(mnvh1d_tot_data_multiflux_xs, is_enu).GetCVHistoWithError(include_stat, not_cov_area_norm),
        // &mnvh1d_adapter(mnvh1d_tot_data_multiflux_xs, is_enu).GetStatError(not_as_frac), 
        // &mnvh1d_adapter(mnvh1d_tot_data_multiflux_xs, is_enu).GetTotalError(not_include_stat, not_as_frac, not_cov_area_norm));

    std::cout << "Total cross section: "
              << mnvh1d_adapter(mnvh1d_tot_data_multiflux_xs, is_enu).GetCVHistoWithError(include_stat, not_cov_area_norm).Integral("width")
              << std::endl;


        /*
    draw_systematics draw_flux(mnvh1d_integrated_flux);
    draw_flux.draw_all(Form("%s-mc-flux-uncertainty", variable_name));
    
    draw_systematics draw_produced(mnvh1d_minerva_mcsig_produced);
    draw_produced.draw_all(Form("%s-mc-produced-uncertainty", variable_name));

    draw_systematics draw_selected(mnvh1d_minerva_mcsig_selected);
    draw_selected.draw_all(Form("%s-mc-selected-uncertainty", variable_name));

    draw_systematics draw_eff(mnvh1d_minerva_efficiency);
    draw_eff.draw_all(Form("%s-mc-efficiency-uncertainty", variable_name));

    draw_systematics draw_raw(mnvh1d_minerva_mcsig_raw);
    draw_raw.draw_all(Form("%s-mc-raw-uncertainty", variable_name));
    
    draw_systematics draw_bkg(mnvh1d_minerva_data_bkg);
    draw_bkg.draw_all(Form("%s-data-background-uncertainty", variable_name));

    draw_systematics draw_bkgsub(mnvh1d_minerva_data_bkgsub);
    draw_bkgsub.draw_all(Form("%s-data-bkgsub-uncertainty", variable_name));

    draw_systematics draw_unfolded(mnvh1d_minerva_data_unfolded);
    draw_unfolded.draw_all(Form("%s-data-unfolded-uncertainty", variable_name));

    draw_systematics draw_effcorr(mnvh1d_minerva_data_effcorrected);
    draw_effcorr.draw_all(Form("%s-data-effcorr-uncertainty", variable_name));

    draw_systematics draw_xsec(mnvh1d_minerva_data_multiflux_xs);
    draw_xsec.draw_all(Form("%s-data-xsec-uncertainty", variable_name));
        */
    
    TCanvas* c = new TCanvas("c1");
    
    
        // plots cross sections by minerva5 data
    TGaxis::SetMaxDigits(5);
    plotter->axis_minimum = 0.001;
    plotter->axis_maximum = plot_attr_variable.xs_area_ymax;
    



    c = new TCanvas("c2");
    plotter->DrawDataMC(get_h1d_with_total_pot_err(mnvh1d_minerva_data_multiflux_xs, is_enu),
                        mnvh1d_adapter_ptr(mnvh1d_gnsig_multiflux_xs,       is_enu),
                        1.0,
                        plot_attr_variable.legend_position.c_str(), false);
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    write_pot_exposure(plotter, data_minerva_pot);
    c->Print(Form("%s-mnv-data-mc-xs-multiflux-pot%s",variable_name, __figure_format.c_str()));
    delete c;
    
    c = new TCanvas("c3");
    plotter->DrawDataMC(get_h1d_with_total_area_err(mnvh1d_minerva_data_multiflux_xs, is_enu),
                        mnvh1d_adapter_ptr(mnvh1d_gnsig_multiflux_xs,       is_enu),
                        calc_area_ratio(mnvh1d_minerva_data_multiflux_xs, mnvh1d_gnsig_multiflux_xs),
                        plot_attr_variable.legend_position.c_str(), false);
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    write_area_exposure(plotter, data_minerva_pot);
    c->Print(Form("%s-mnv-data-mc-xs-multiflux-area%s",variable_name, __figure_format.c_str()));
    delete c;

        // plot cross sections by downstream data
    c = new TCanvas("c4");
    plotter->DrawDataMC(get_h1d_with_total_pot_err(mnvh1d_downstream_data_multiflux_xs, is_enu),
                        mnvh1d_adapter_ptr(mnvh1d_gnsig_multiflux_xs,          is_enu),
                        1.0,
                        plot_attr_variable.legend_position.c_str(), false);
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    write_pot_exposure(plotter, equivalent_pot);
    c->Print(Form("%s-frz-data-mc-xs-multiflux-pot%s",variable_name, __figure_format.c_str()));
    delete c;
    
    c = new TCanvas("c5");
    plotter->DrawDataMC(get_h1d_with_total_area_err(mnvh1d_downstream_data_multiflux_xs, is_enu),
                        mnvh1d_adapter_ptr(mnvh1d_gnsig_multiflux_xs,             is_enu),
                        calc_area_ratio(mnvh1d_downstream_data_multiflux_xs, mnvh1d_gnsig_multiflux_xs),
                        plot_attr_variable.legend_position.c_str(), false);
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    write_area_exposure(plotter, equivalent_pot);
    c->Print(Form("%s-frz-data-mc-xs-multiflux-area%s",variable_name, __figure_format.c_str()));
    delete c;

    
        // cross sections using combined data
    plotter->axis_maximum = plot_attr_variable.xs_area_ymax;
    c = new TCanvas("c6");
    plotter->DrawDataMC(get_h1d_with_total_pot_err(mnvh1d_tot_data_multiflux_xs, is_enu),
                        mnvh1d_adapter_ptr(mnvh1d_gnsig_multiflux_xs,   is_enu),
                        1.0,
                        plot_attr_variable.legend_position.c_str(), false);
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    write_pot_exposure(plotter, recorded_pot);
    c->Print(Form("%s-tot-data-mc-xs-multiflux-pot%s",variable_name, __figure_format.c_str()));
    delete c;
    
    c = new TCanvas("c7");
    plotter->DrawDataMC(get_h1d_with_total_area_err(mnvh1d_tot_data_multiflux_xs, is_enu),
                        mnvh1d_adapter_ptr(mnvh1d_gnsig_multiflux_xs,      is_enu),
                        calc_area_ratio(mnvh1d_tot_data_multiflux_xs, mnvh1d_gnsig_multiflux_xs),
                        plot_attr_variable.legend_position.c_str(), false);
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    write_area_exposure(plotter, recorded_pot);
    c->Print(Form("%s-tot-data-mc-xs-multiflux-area%s",variable_name, __figure_format.c_str()));
    delete c;
    
    c = new TCanvas("c8");
    plotter->DrawDataMCRatio(get_h1d_with_total_pot_err(mnvh1d_tot_data_multiflux_xs, is_enu),
                             mnvh1d_adapter_ptr(mnvh1d_gnsig_multiflux_xs,  0),
                             1.0,
                             true, // draw one line
                             0.0,  // plotMin
                             2.0   // plotMax
                             ); 
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    write_pot_exposure(plotter, recorded_pot);
    c->Print(Form("%s-tot-data-mc-xs-multiflux-pot-ratio%s",variable_name, __figure_format.c_str()));
    delete c;

    
    
    // draw minerva5 and downstream comparison
    mnvh1d_minerva_data_multiflux_xs->SetTitle("minerva5");
    mnvh1d_downstream_data_multiflux_xs->SetTitle("downstream");
    mnvh1d_tot_data_multiflux_xs->SetTitle("combined");
    mnvh1d_avg1_data_multiflux_xs->SetTitle("averaged1");
    mnvh1d_avg2_data_multiflux_xs->SetTitle("averaged2");


    c = new TCanvas("c9");
    draw_data_datavariation(plotter,
                            get_h1d_with_stat_err(mnvh1d_minerva_data_multiflux_xs,    is_enu),
                            get_h1d_with_stat_err(mnvh1d_downstream_data_multiflux_xs, is_enu),
                            NULL,
                            1.0, 1.0, 
                            kSolid, kDashed,
                            kRed, kRed);

    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    write_pot_exposure(plotter, data_minerva_pot);
    c->Print(Form("%s-dataset1-data-mc-xs-multiflux-pot%s",variable_name, __figure_format.c_str()));
    delete c;

        // draw minerva5 and combined comparison
    c = new TCanvas("c10");
    draw_data_datavariation(plotter,
                            get_h1d_with_stat_err(mnvh1d_minerva_data_multiflux_xs,  is_enu),
                            get_h1d_with_stat_err(mnvh1d_tot_data_multiflux_xs,      is_enu),
                            NULL,
                            1.0, 1.0, 
                            kSolid, kDashed,
                            kRed, kRed);

    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    write_pot_exposure(plotter, data_minerva_pot);
    c->Print(Form("%s-dataset2-data-mc-xs-multiflux-pot%s",variable_name, __figure_format.c_str()));
    delete c;

        // draw combined and averaged 1 comparison
    c = new TCanvas("c11");
    draw_data_datavariation(plotter,
                            get_h1d_with_stat_err(mnvh1d_avg1_data_multiflux_xs, is_enu),
                            get_h1d_with_stat_err(mnvh1d_tot_data_multiflux_xs,  is_enu),
                            NULL,
                            1.0, 1.0, 
                            kSolid, kDashed,
                            kRed, kRed);
    
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    write_pot_exposure(plotter, recorded_pot);
    c->Print(Form("%s-dataset3-data-mc-xs-multiflux-pot%s",variable_name, __figure_format.c_str()));
    delete c;


        // draw combined and averaged 2 comparison
    c = new TCanvas("c12");
    draw_data_datavariation(plotter,
                            get_h1d_with_stat_err(mnvh1d_avg2_data_multiflux_xs, is_enu),
                            get_h1d_with_stat_err(mnvh1d_tot_data_multiflux_xs,  is_enu),
                            NULL,
                            1.0, 1.0, 
                            kSolid, kDashed,
                            kRed, kRed);
    
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    write_pot_exposure(plotter, recorded_pot);
    c->Print(Form("%s-dataset4-data-mc-xs-multiflux-pot%s",variable_name, __figure_format.c_str()));
    delete c;

    
        // draw data vs GENIE with and without FSI
    mnvh1d_gnsig_multiflux_xs->SetTitle("GENIE w/ FSI");
    mnvh1d_gnsig_multiflux_xs_nofsi->SetTitle("GENIE w/o FSI");

    c = new TCanvas("c13");
    draw_data_mcvariation(plotter,
                          get_h1d_with_total_pot_err(mnvh1d_tot_data_multiflux_xs,  is_enu),
                          get_h1d_with_stat_err(mnvh1d_tot_data_multiflux_xs,   is_enu),
                          mnvh1d_adapter_ptr(mnvh1d_gnsig_multiflux_xs,         is_enu),
                          mnvh1d_adapter_ptr(mnvh1d_gnsig_multiflux_xs_nofsi,   is_enu),
                          NULL,
                          NULL,
                          NULL,
                          1.0, 1.0, 1.0, 1.0, 1.0,
                          kSolid, kDashed, kSolid, kSolid, kSolid,
                          kRed, kRed, kBlue, kBlue, kBlue,
                          plot_attr_variable.legend_position.c_str(), false);
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    write_pot_exposure(plotter, recorded_pot);
    c->Print(Form("%s-data-mc-xs-multiflux-variation-pot%s",variable_name, __figure_format.c_str()));
    delete c;
    
    c = new TCanvas("c14");
    draw_data_mcvariation(plotter,
                          get_h1d_with_total_area_err(mnvh1d_tot_data_multiflux_xs, is_enu),
                          get_h1d_with_stat_err(mnvh1d_tot_data_multiflux_xs, is_enu),
                          mnvh1d_adapter_ptr(mnvh1d_gnsig_multiflux_xs,         is_enu),
                          mnvh1d_adapter_ptr(mnvh1d_gnsig_multiflux_xs_nofsi,   is_enu),
                          NULL,
                          NULL,
                          NULL,
                          calc_area_ratio(mnvh1d_tot_data_multiflux_xs, mnvh1d_gnsig_multiflux_xs),
                          calc_area_ratio(mnvh1d_tot_data_multiflux_xs, mnvh1d_gnsig_multiflux_xs_nofsi),
                          1.0, 1.0, 1.0,
                          kSolid, kDashed, kSolid, kSolid, kSolid,
                          kRed, kRed, kBlue, kBlue, kBlue,
                          plot_attr_variable.legend_position.c_str(), false);
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    write_area_exposure(plotter, recorded_pot);
    c->Print(Form("%s-data-mc-xs-multiflux-variation-area%s",variable_name, __figure_format.c_str()));
    delete c;



    
        // draw data vs nuwro prediction
        //TFile* nuwro_xs_file = new TFile("nuwro-cross-section.root", "read");
        //TH1D* hnuwro_multiflux_xs = (TH1D*) nuwro_xs_file->Get(Form("nuwro_%s", __variable.c_str()));
        //if (__variable == "pimom" || __variable == "theta") assert(hnuwro_multiflux_xs);

        // draw data vs NEUT prediction
        //TFile* neut_xs_file = new TFile("neut-cross-section.root", "read");
        //TH1D* hneut_multiflux_xs = (TH1D*) neut_xs_file->Get(Form("neut_%s", __variable.c_str()));
        //if (__variable == "pimom" || __variable == "theta") assert(hneut_multiflux_xs);

#ifdef model_comparison
    MnvH1D* mnvh1d_neut_xs(NULL);
    MnvH1D* mnvh1d_nuwro_xs(NULL);
    
    if (!__neut_file.empty())  mnvh1d_neut_xs  = neut_to_mnvh1d(__neut_file, __variable);
    if (!__nuwro_file.empty()) mnvh1d_nuwro_xs = nuwro_to_mnvh1d(__nuwro_file, __variable, 2); 

        //assert(mnvh1d_neut_xs);
        //assert(mnvh1d_nuwro_xs);
    
    if ((mnvh1d_neut_xs && mnvh1d_nuwro_xs)) {
        
            //MnvH1D*  mnvh1d_nuwro_multiflux_xs = new MnvH1D(*hnuwro_multiflux_xs);
            //mnvh1d_nuwro_multiflux_xs->SetNormBinWidth(width);
                
            //MnvH1D* mnvh1d_neut_multiflux_xs = new MnvH1D(*hneut_multiflux_xs);
            //double new_width = get_min_binwidth(mnvh1d_neut_multiflux_xs);
            //mnvh1d_neut_multiflux_xs->SetNormBinWidth(new_width);
            //std::cout << "new_width: " << new_width << std::endl;

            // Unfortunately the pion variables and the muon variables have difference format:
            // the (old) pion-variable predictions need binwidth normalized while the (new) muon-variable
            // predictions are simply curves, which don't need binwidth normalized
        if (__variable == "pimom" || __variable == "theta") {
            mnvh1d_neut_xs->Scale(width, "width");
            mnvh1d_nuwro_xs->Scale(width, "width");
        }
        
        mnvh1d_nuwro_xs->SetTitle("NuWro");
        mnvh1d_neut_xs->SetTitle("NEUT");
        gStyle->SetNdivisions(plot_attr_variable.xndivision, "X");
        gStyle->SetNdivisions(plot_attr_variable.yndivision, "Y");
        

            // generator comparison, no GENIE w/o FSI
        c = new TCanvas("c15");
        draw_data_mcvariation(plotter,
                              get_h1d_with_total_pot_err(mnvh1d_tot_data_multiflux_xs, is_enu),
                              get_h1d_with_stat_err(mnvh1d_tot_data_multiflux_xs,  is_enu),
                              mnvh1d_adapter_ptr(mnvh1d_gnsig_multiflux_xs,        is_enu),
                              NULL,
                              mnvh1d_nuwro_xs,
                              mnvh1d_neut_xs,
                              NULL,
                              1.0, 1.0, 1.0, 1.0, 1.0,
//                              kSolid, kSolid, kSolid, kSolid, kSolid,  // for presentation
                              kSolid, kDashed, 7, 5, kSolid,          // for paper
                              kRed, kRed, kOrange-6, kGreen+2, kBlue,
                              plot_attr_variable.legend_position.c_str(),
                              false, // use hist title
                              is_enu);
        plotter->WritePreliminary("TL");
        write_hist_title(plotter);
        write_pot_exposure(plotter, recorded_pot);
        c->Print(Form("%s-tot-data-mc-xs-multiflux-models1-pot%s",variable_name, __figure_format.c_str()));
        delete c;


        c = new TCanvas("c16");
        draw_data_mcvariation(plotter,
                              get_h1d_with_total_area_err(mnvh1d_tot_data_multiflux_xs, is_enu),
                              get_h1d_with_stat_err(mnvh1d_tot_data_multiflux_xs,       is_enu),
                              mnvh1d_adapter_ptr(mnvh1d_gnsig_multiflux_xs,             is_enu),
                              NULL,
                              mnvh1d_nuwro_xs,
                              mnvh1d_neut_xs, 
                              NULL,
                              calc_area_ratio(mnvh1d_adapter_ptr(mnvh1d_tot_data_multiflux_xs, is_enu),
                                              mnvh1d_adapter_ptr(mnvh1d_gnsig_multiflux_xs,    is_enu)),    
                              1.0,
                              calc_area_ratio(mnvh1d_adapter_ptr(mnvh1d_tot_data_multiflux_xs, is_enu), mnvh1d_nuwro_xs),
                              calc_area_ratio(mnvh1d_adapter_ptr(mnvh1d_tot_data_multiflux_xs, is_enu), mnvh1d_neut_xs),  
                              1.0,
//                              kSolid, kSolid, kSolid, kSolid, kSolid, // for presentation
                              kSolid, kDashed, 7, 5, kSolid,         // for paper
                              kRed, kRed, kOrange-6, kGreen+2, kBlue,
                              plot_attr_variable.legend_position.c_str(),
                              false, // use hist title
                              is_enu);
        plotter->WritePreliminary("TL");
        write_hist_title(plotter);
        write_area_exposure(plotter, recorded_pot);
        c->Print(Form("%s-tot-data-mc-xs-multiflux-models1-area%s",variable_name, __figure_format.c_str()));
        delete c;


            // generator comparison, with GENIE w/o FSI
        c = new TCanvas("c17");
        draw_data_mcvariation(plotter,
                              get_h1d_with_total_pot_err(mnvh1d_tot_data_multiflux_xs,  is_enu),
                              get_h1d_with_stat_err(mnvh1d_tot_data_multiflux_xs,   is_enu),
                              mnvh1d_adapter_ptr(mnvh1d_gnsig_multiflux_xs,         is_enu),
                              mnvh1d_adapter_ptr(mnvh1d_gnsig_multiflux_xs_nofsi,   is_enu),
                              mnvh1d_nuwro_xs,
                              mnvh1d_neut_xs,
                              NULL,
                              1.0, 1.0, 1.0, 1.0, 1.0,
//                              kSolid, kDashed, kSolid, kSolid, kSolid,  // for presentation
                              kSolid, kDashed, 7, 5, kSolid,          // for paper
                              kRed, kRed, kOrange-6, kGreen+2, kBlue,
                              plot_attr_variable.legend_position.c_str(),
                              false, // use hist title
                              is_enu);
        plotter->WritePreliminary("TL");
        write_hist_title(plotter);
        write_pot_exposure(plotter, recorded_pot);
        c->Print(Form("%s-tot-data-mc-xs-multiflux-models2-pot%s",variable_name, __figure_format.c_str()));
        delete c;

        c = new TCanvas("c18");
        draw_data_mcvariation_ratio(plotter,
                                    get_h1d_with_total_pot_err(mnvh1d_tot_data_multiflux_xs,  is_enu),
                                    get_h1d_with_stat_err(mnvh1d_tot_data_multiflux_xs,   is_enu),
                                    mnvh1d_adapter_ptr(mnvh1d_gnsig_multiflux_xs,         is_enu),
                                    mnvh1d_adapter_ptr(mnvh1d_gnsig_multiflux_xs_nofsi,   is_enu),
                                    NULL,
                                    NULL,
                                    NULL,
                                    1.0, 1.0, 1.0, 1.0, 1.0,
//                              kSolid, kDashed, kSolid, kSolid, kSolid,  // for presentation
                                    kSolid, kDashed, 7, 5, kSolid,          // for paper
                                    kRed, kRed, kOrange-6, kGreen+2, kBlue,
                                    plot_attr_variable.legend_position.c_str(),
                                    false, // use hist title
                                    is_enu);
        plotter->WritePreliminary("TL");
        write_hist_title(plotter);
        write_pot_exposure(plotter, recorded_pot);
        c->Print(Form("%s-tot-data-mc-xs-multiflux-models2-pot-ratio%s",variable_name, __figure_format.c_str()));
        delete c;

        
        c = new TCanvas("c19");
        draw_data_mcvariation(plotter,
                              get_h1d_with_total_area_err(mnvh1d_tot_data_multiflux_xs, is_enu),
                              get_h1d_with_stat_err(mnvh1d_tot_data_multiflux_xs,       is_enu),
                              mnvh1d_adapter_ptr(mnvh1d_gnsig_multiflux_xs,             is_enu),
                              mnvh1d_adapter_ptr(mnvh1d_gnsig_multiflux_xs_nofsi,       is_enu),
                              mnvh1d_nuwro_xs,
                              mnvh1d_neut_xs, 
                              NULL,
                              calc_area_ratio(mnvh1d_adapter_ptr(mnvh1d_tot_data_multiflux_xs,     is_enu),
                                              mnvh1d_adapter_ptr(mnvh1d_gnsig_multiflux_xs,        is_enu)), 
                              calc_area_ratio(mnvh1d_adapter_ptr(mnvh1d_tot_data_multiflux_xs,     is_enu),
                                              mnvh1d_adapter_ptr(mnvh1d_gnsig_multiflux_xs_nofsi,  is_enu)), 
                              calc_area_ratio(mnvh1d_adapter_ptr(mnvh1d_tot_data_multiflux_xs,     is_enu),  mnvh1d_nuwro_xs), 
                              calc_area_ratio(mnvh1d_adapter_ptr(mnvh1d_tot_data_multiflux_xs,     is_enu),  mnvh1d_neut_xs),  
                              1.0, 
//                              kSolid, kDashed, kSolid, kSolid, kSolid,     // for presentation
                              kSolid, kDashed, 7, 5, kSolid,            // for paper
                              kRed, kRed, kOrange-6, kGreen+2, kBlue, 
                              plot_attr_variable.legend_position.c_str(),
                              false, // use hist title
                              is_enu);
        plotter->WritePreliminary("TL");
        write_hist_title(plotter);
        write_area_exposure(plotter, recorded_pot);
        c->Print(Form("%s-tot-data-mc-xs-multiflux-models2-area%s",variable_name, __figure_format.c_str()));
        delete c;


            // Calculate Chi2 for different comparisons" 
        int ndof_genie       = 0;
        int ndof_genie_nofsi = 0;
        int ndof_nuwro       = 0;
        int ndof_neut        = 0;
        
        int ndof_genie_shape       = 0;
        int ndof_genie_nofsi_shape = 0;
        int ndof_nuwro_shape       = 0;
        int ndof_neut_shape        = 0;
        
    
        bool shapeOnly = true;
        bool absolute = !shapeOnly;
        bool useDataErrorMatrix = true;
        
            // Absolute chi2
        double chi2_genie
            = plotter->Chi2DataMC(mnvh1d_tot_data_multiflux_xs,
                                  mnvh1d_gnsig_multiflux_xs,
                                  ndof_genie,
                                  1.0,
                                  useDataErrorMatrix,
                                  absolute);
    
        double chi2_genie_nofsi
            = plotter->Chi2DataMC(mnvh1d_tot_data_multiflux_xs,
                                  mnvh1d_gnsig_multiflux_xs_nofsi,
                                  ndof_genie_nofsi,
                                  1.0,
                                  useDataErrorMatrix,
                                  absolute);
        
        double chi2_nuwro
            = plotter->Chi2DataMC(mnvh1d_tot_data_multiflux_xs,
                                  mnvh1d_nuwro_xs,
                                  ndof_nuwro,
                                  1.0,
                                  useDataErrorMatrix,
                                  absolute);

        double chi2_neut
            = plotter->Chi2DataMC(mnvh1d_tot_data_multiflux_xs,
                                  mnvh1d_neut_xs,
                                  ndof_neut,
                                  1.0,
                                  useDataErrorMatrix,
                                  absolute);
        
        
            // shape Chi2
        double chi2_genie_shape
            = plotter->Chi2DataMC(mnvh1d_tot_data_multiflux_xs,
                                  mnvh1d_gnsig_multiflux_xs,
                                  ndof_genie_shape,
                                  calc_area_ratio(mnvh1d_tot_data_multiflux_xs, mnvh1d_gnsig_multiflux_xs),
                                  useDataErrorMatrix,
                                  shapeOnly);
        
        double chi2_genie_nofsi_shape
            = plotter->Chi2DataMC(mnvh1d_tot_data_multiflux_xs,
                                  mnvh1d_gnsig_multiflux_xs_nofsi,
                                  ndof_genie_nofsi_shape,
                                  calc_area_ratio(mnvh1d_tot_data_multiflux_xs, mnvh1d_gnsig_multiflux_xs_nofsi),
                                  useDataErrorMatrix,
                                  shapeOnly);

        double chi2_nuwro_shape
            = plotter->Chi2DataMC(mnvh1d_tot_data_multiflux_xs,
                                  mnvh1d_nuwro_xs,
                                  ndof_nuwro_shape,
                                  calc_area_ratio(mnvh1d_tot_data_multiflux_xs, mnvh1d_nuwro_xs),
                                  useDataErrorMatrix,
                                  shapeOnly);

        double chi2_neut_shape
            = plotter->Chi2DataMC(mnvh1d_tot_data_multiflux_xs,
                                  mnvh1d_neut_xs,
                                  ndof_neut_shape,
                                  calc_area_ratio(mnvh1d_tot_data_multiflux_xs, mnvh1d_neut_xs),
                                  useDataErrorMatrix,
                                  shapeOnly);
        

        
        std::cout << "Calculate chi2 for " << __variable << std::endl;
        std::cout << "Absolute chi2  " << std::endl;
        std::cout << "GENIE          " << chi2_genie       << " " << ndof_genie       << std::endl;
        std::cout << "GENIE (No FSI) " << chi2_genie_nofsi << " " << ndof_genie_nofsi << std::endl;
        std::cout << "NuWro          " << chi2_nuwro       << " " << ndof_nuwro       << std::endl;
        std::cout << "NEUT           " << chi2_neut        << " " << ndof_neut        << std::endl;        
    
        std::cout << "Shape only chi2" << std::endl;
        std::cout << "GENIE          " << chi2_genie_shape         << " " << ndof_genie_shape       << std::endl;
        std::cout << "GENIE (No FSI) " << chi2_genie_nofsi_shape   << " " << ndof_genie_nofsi_shape << std::endl;
        std::cout << "NuWro          " << chi2_nuwro_shape         << " " << ndof_nuwro_shape       << std::endl;
        std::cout << "NEUT           " << chi2_neut_shape          << " " << ndof_neut_shape        << std::endl;

        
    }

#endif
            // follow pi0 signal in the nucleus
    MnvH1D* mnvh1d_modified_signal   = static_cast<MnvH1D*>( minerva_efficiency_file->Get(Form("mc_cut00_modified_%s",   variable_name)));
    MnvH1D* mnvh1d_nonint_signal     = static_cast<MnvH1D*>( minerva_efficiency_file->Get(Form("mc_cut00_nonint_%s",     variable_name)));
    MnvH1D* mnvh1d_elastic_signal    = static_cast<MnvH1D*>( minerva_efficiency_file->Get(Form("mc_cut00_elastic_%s",    variable_name)));
    MnvH1D* mnvh1d_inelastic2_signal = static_cast<MnvH1D*>( minerva_efficiency_file->Get(Form("mc_cut00_inelastic2_%s", variable_name)));
    MnvH1D* mnvh1d_cex_signal        = static_cast<MnvH1D*>( minerva_efficiency_file->Get(Form("mc_cut00_cex_%s",        variable_name)));
    MnvH1D* mnvh1d_multipion_signal  = static_cast<MnvH1D*>( minerva_efficiency_file->Get(Form("mc_cut00_multipion_%s",  variable_name)));
    
    assert(mnvh1d_modified_signal);
    assert(mnvh1d_nonint_signal);
    assert(mnvh1d_elastic_signal);
    assert(mnvh1d_inelastic2_signal);
    assert(mnvh1d_cex_signal);
    assert(mnvh1d_multipion_signal);

        
    MnvH1D* mnvh1d_mcsig_multiflux_xs_modified   = calc_xs(mnvh1d_modified_signal,   mnvh1d_integrated_flux, mc_pot, true);
    MnvH1D* mnvh1d_mcsig_multiflux_xs_nonint     = calc_xs(mnvh1d_nonint_signal,     mnvh1d_integrated_flux, mc_pot, true);
    MnvH1D* mnvh1d_mcsig_multiflux_xs_elastic    = calc_xs(mnvh1d_elastic_signal,    mnvh1d_integrated_flux, mc_pot, true);
    MnvH1D* mnvh1d_mcsig_multiflux_xs_inelastic2 = calc_xs(mnvh1d_inelastic2_signal, mnvh1d_integrated_flux, mc_pot, true);
    MnvH1D* mnvh1d_mcsig_multiflux_xs_cex        = calc_xs(mnvh1d_cex_signal,        mnvh1d_integrated_flux, mc_pot, true);
    MnvH1D* mnvh1d_mcsig_multiflux_xs_multipion  = calc_xs(mnvh1d_multipion_signal,  mnvh1d_integrated_flux, mc_pot, true);
        
 

    assert(mnvh1d_mcsig_multiflux_xs_modified);
    assert(mnvh1d_mcsig_multiflux_xs_nonint);
    assert(mnvh1d_mcsig_multiflux_xs_elastic);
    assert(mnvh1d_mcsig_multiflux_xs_inelastic2);
    assert(mnvh1d_mcsig_multiflux_xs_cex);
    assert(mnvh1d_mcsig_multiflux_xs_multipion);
    
        
    mnvh1d_mcsig_multiflux_xs_modified->Scale(plot_attr_variable.scale);
    mnvh1d_mcsig_multiflux_xs_nonint->Scale(plot_attr_variable.scale);
    mnvh1d_mcsig_multiflux_xs_elastic->Scale(plot_attr_variable.scale);
    mnvh1d_mcsig_multiflux_xs_inelastic2->Scale(plot_attr_variable.scale);
    mnvh1d_mcsig_multiflux_xs_cex->Scale(plot_attr_variable.scale);
    mnvh1d_mcsig_multiflux_xs_multipion->Scale(plot_attr_variable.scale);
    
    mnvh1d_mcsig_multiflux_xs_modified->SetNormBinWidth(width);
    mnvh1d_mcsig_multiflux_xs_nonint->SetNormBinWidth(width);
    mnvh1d_mcsig_multiflux_xs_elastic->SetNormBinWidth(width);
    mnvh1d_mcsig_multiflux_xs_inelastic2->SetNormBinWidth(width);
    mnvh1d_mcsig_multiflux_xs_cex->SetNormBinWidth(width);
    mnvh1d_mcsig_multiflux_xs_multipion->SetNormBinWidth(width);

    
    c = new TCanvas("c20");
    TObjArray* objarray_xs_multiflux_origin_breakdown = new TObjArray;
    
    mnvh1d_mcsig_multiflux_xs_nonint->SetTitle("#pi^{0} Non-interacting");
    mnvh1d_mcsig_multiflux_xs_elastic->SetTitle("#pi^{0} Elastic");
    mnvh1d_mcsig_multiflux_xs_inelastic2->SetTitle("#pi^{0} Inelastic");
    mnvh1d_mcsig_multiflux_xs_cex->SetTitle("#pi^{-} #rightarrow #pi^{0}");
    mnvh1d_mcsig_multiflux_xs_multipion->SetTitle("Multi-#pi #rightarrow #pi^{0}");

    
    if (true) {

            // decorations
        mnvh1d_mcsig_multiflux_xs_multipion->SetLineColor(kCyan-8);
        mnvh1d_mcsig_multiflux_xs_multipion->SetFillColor(kCyan-8);
        mnvh1d_mcsig_multiflux_xs_multipion->SetFillStyle(1001);   
        
        mnvh1d_mcsig_multiflux_xs_cex->SetLineColor(kGreen-2);
        mnvh1d_mcsig_multiflux_xs_cex->SetFillColor(kGreen-2);
        mnvh1d_mcsig_multiflux_xs_cex->SetFillStyle(1001);         
        
        mnvh1d_mcsig_multiflux_xs_inelastic2->SetLineColor(kGray+1);
        mnvh1d_mcsig_multiflux_xs_inelastic2->SetFillColor(kGray+1);
        mnvh1d_mcsig_multiflux_xs_inelastic2->SetFillStyle(1001);

        mnvh1d_mcsig_multiflux_xs_elastic->SetLineColor(kOrange+6);
        mnvh1d_mcsig_multiflux_xs_elastic->SetFillColor(kOrange+6);
        mnvh1d_mcsig_multiflux_xs_elastic->SetFillStyle(1001);     
        

        mnvh1d_mcsig_multiflux_xs_nonint->SetLineColor(kRed+1);
        mnvh1d_mcsig_multiflux_xs_nonint->SetFillColor(kRed+1);
        mnvh1d_mcsig_multiflux_xs_nonint->SetFillStyle(1001); 
   
    } else {
        
            // decorations
        mnvh1d_mcsig_multiflux_xs_multipion->SetLineColor(17);
        mnvh1d_mcsig_multiflux_xs_multipion->SetFillColor(17);
        mnvh1d_mcsig_multiflux_xs_multipion->SetFillStyle(1001);   
        
        mnvh1d_mcsig_multiflux_xs_cex->SetLineColor(kMagenta);     // old color 17
        mnvh1d_mcsig_multiflux_xs_cex->SetFillColor(kMagenta);     // old color 17
        mnvh1d_mcsig_multiflux_xs_cex->SetFillStyle(3013);         
        

        mnvh1d_mcsig_multiflux_xs_inelastic2->SetLineColor(kRed);
        mnvh1d_mcsig_multiflux_xs_inelastic2->SetFillColor(kRed);
        mnvh1d_mcsig_multiflux_xs_inelastic2->SetFillStyle(1001);

        mnvh1d_mcsig_multiflux_xs_elastic->SetLineColor(kGreen);
        mnvh1d_mcsig_multiflux_xs_elastic->SetFillColor(kGreen);
        mnvh1d_mcsig_multiflux_xs_elastic->SetFillStyle(3007);     
        

        mnvh1d_mcsig_multiflux_xs_nonint->SetLineColor(kBlue);
        mnvh1d_mcsig_multiflux_xs_nonint->SetFillColor(kBlue);
        mnvh1d_mcsig_multiflux_xs_nonint->SetFillStyle(3010);      // old fill 3006
        
    }
    
    objarray_xs_multiflux_origin_breakdown->Add(mnvh1d_mcsig_multiflux_xs_nonint);
    objarray_xs_multiflux_origin_breakdown->Add(mnvh1d_mcsig_multiflux_xs_elastic);
    objarray_xs_multiflux_origin_breakdown->Add(mnvh1d_mcsig_multiflux_xs_inelastic2);
    objarray_xs_multiflux_origin_breakdown->Add(mnvh1d_mcsig_multiflux_xs_cex);
    objarray_xs_multiflux_origin_breakdown->Add(mnvh1d_mcsig_multiflux_xs_multipion);

    int old_line_width = plotter->mc_line_width;
    plotter->mc_line_width = 1;
    

    if (is_enu) plotter->draw_normalized_to_bin_width = false;
        // redraw using the customized version

    int old_xndivision = gStyle->GetNdivisions("X");
    int old_yndivision = gStyle->GetNdivisions("Y");
    gStyle->SetNdivisions(plot_attr_variable.xndivision, "X");
    gStyle->SetNdivisions(plot_attr_variable.yndivision, "Y");
    
    c = new TCanvas("c21");
    draw_data_mcstacked(plotter,
                        mnvh1d_tot_data_multiflux_xs,
                        objarray_xs_multiflux_origin_breakdown,
                        1.0,
                        plot_attr_variable.legend_position.c_str(),
                            //"Data", 2, 1, 3001,
                        "Data", -1, 1, -1,
                        plot_attr_variable.xs_xlabel.c_str(),
                        plot_attr_variable.xs_ylabel.c_str(),
                        not_cov_area_norm,
                        is_enu);
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    write_pot_exposure(plotter, recorded_pot);
    c->Print(Form("%s-data-mc-xs-multiflux-origin-breakdown-pot%s", variable_name, __figure_format.c_str()));
    delete c;
    
        // area normalized
    c = new TCanvas("c22");
    draw_data_mcstacked(plotter,
                        mnvh1d_tot_data_multiflux_xs,
                        objarray_xs_multiflux_origin_breakdown,
                        calc_area_ratio(mnvh1d_tot_data_multiflux_xs, mnvh1d_gnsig_multiflux_xs),
                        plot_attr_variable.legend_position.c_str(),
                            //"Data", 2, 1, 3001,
                        "Data", -1, 1, -1,
                        plot_attr_variable.xs_xlabel.c_str(),
                        plot_attr_variable.xs_ylabel.c_str(),
                        cov_area_norm,
                        is_enu); 
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    write_area_exposure(plotter, recorded_pot);
    c->Print(Form("%s-data-mc-xs-multiflux-origin-breakdown-area%s", variable_name, __figure_format.c_str()));
    delete c;


        // old histogram files do not have the process breakdown histograms
    if (!__legacy_code) {
        
            // type3 signal breakdown
        MnvH1D* mnvh1d_type3_cat1 = static_cast<MnvH1D*>(xsec_file->Get(Form("mc_cut00_type3_%s_1", variable_name)));
        MnvH1D* mnvh1d_type3_cat2 = static_cast<MnvH1D*>(xsec_file->Get(Form("mc_cut00_type3_%s_2", variable_name)));
        MnvH1D* mnvh1d_type3_cat3 = static_cast<MnvH1D*>(xsec_file->Get(Form("mc_cut00_type3_%s_3", variable_name)));
        MnvH1D* mnvh1d_type3_cat4 = static_cast<MnvH1D*>(xsec_file->Get(Form("mc_cut00_type3_%s_4", variable_name)));
        MnvH1D* mnvh1d_type3_cat5 = static_cast<MnvH1D*>(xsec_file->Get(Form("mc_cut00_type3_%s_5", variable_name)));
        
        assert(mnvh1d_type3_cat1);
        assert(mnvh1d_type3_cat2);
        assert(mnvh1d_type3_cat3);
        assert(mnvh1d_type3_cat4);
        assert(mnvh1d_type3_cat5);
        
        MnvH1D* mnvh1d_mcsig_multiflux_xs_type3_cat1 = calc_xs(mnvh1d_type3_cat1, mnvh1d_integrated_flux, mc_pot, true);
        MnvH1D* mnvh1d_mcsig_multiflux_xs_type3_cat2 = calc_xs(mnvh1d_type3_cat2, mnvh1d_integrated_flux, mc_pot, true);
        MnvH1D* mnvh1d_mcsig_multiflux_xs_type3_cat3 = calc_xs(mnvh1d_type3_cat3, mnvh1d_integrated_flux, mc_pot, true);
        MnvH1D* mnvh1d_mcsig_multiflux_xs_type3_cat4 = calc_xs(mnvh1d_type3_cat4, mnvh1d_integrated_flux, mc_pot, true);
        MnvH1D* mnvh1d_mcsig_multiflux_xs_type3_cat5 = calc_xs(mnvh1d_type3_cat5, mnvh1d_integrated_flux, mc_pot, true);
        
        assert(mnvh1d_mcsig_multiflux_xs_type3_cat1);
        assert(mnvh1d_mcsig_multiflux_xs_type3_cat2);
        assert(mnvh1d_mcsig_multiflux_xs_type3_cat3);
        assert(mnvh1d_mcsig_multiflux_xs_type3_cat4);
        assert(mnvh1d_mcsig_multiflux_xs_type3_cat5);
        
        mnvh1d_mcsig_multiflux_xs_type3_cat1->Scale(plot_attr_variable.scale);
        mnvh1d_mcsig_multiflux_xs_type3_cat2->Scale(plot_attr_variable.scale);
        mnvh1d_mcsig_multiflux_xs_type3_cat3->Scale(plot_attr_variable.scale);
        mnvh1d_mcsig_multiflux_xs_type3_cat4->Scale(plot_attr_variable.scale);
        mnvh1d_mcsig_multiflux_xs_type3_cat5->Scale(plot_attr_variable.scale);
        
        mnvh1d_mcsig_multiflux_xs_type3_cat1->SetNormBinWidth(width);
        mnvh1d_mcsig_multiflux_xs_type3_cat2->SetNormBinWidth(width);
        mnvh1d_mcsig_multiflux_xs_type3_cat3->SetNormBinWidth(width);
        mnvh1d_mcsig_multiflux_xs_type3_cat4->SetNormBinWidth(width);
        mnvh1d_mcsig_multiflux_xs_type3_cat5->SetNormBinWidth(width);

        
        c = new TCanvas("c23");
        TObjArray* objarray_xs_multiflux_process_breakdown = new TObjArray;
        
        mnvh1d_mcsig_multiflux_xs_type3_cat1->SetTitle("Quasi-Elastic");
        mnvh1d_mcsig_multiflux_xs_type3_cat2->SetTitle("Delta resonance");
        mnvh1d_mcsig_multiflux_xs_type3_cat3->SetTitle("Other resonances");
        mnvh1d_mcsig_multiflux_xs_type3_cat4->SetTitle("Non-Resonant");
        mnvh1d_mcsig_multiflux_xs_type3_cat5->SetTitle("DIS");
        
        if (__carrie_style) {
            
                // decorations
            mnvh1d_mcsig_multiflux_xs_type3_cat1->SetLineColor(kRed+1);
            mnvh1d_mcsig_multiflux_xs_type3_cat1->SetFillColor(kRed+1);
            mnvh1d_mcsig_multiflux_xs_type3_cat1->SetFillStyle(1001);
            
            mnvh1d_mcsig_multiflux_xs_type3_cat2->SetLineColor(kOrange + 6);     
            mnvh1d_mcsig_multiflux_xs_type3_cat2->SetFillColor(kOrange + 6);     
            mnvh1d_mcsig_multiflux_xs_type3_cat2->SetFillStyle(1001);         
            
            mnvh1d_mcsig_multiflux_xs_type3_cat3->SetLineColor(kGray + 1);
            mnvh1d_mcsig_multiflux_xs_type3_cat3->SetFillColor(kGray + 1);
            mnvh1d_mcsig_multiflux_xs_type3_cat3->SetFillStyle(1001);     
            
            mnvh1d_mcsig_multiflux_xs_type3_cat4->SetLineColor(kGreen-2);
            mnvh1d_mcsig_multiflux_xs_type3_cat4->SetFillColor(kGreen-2);
            mnvh1d_mcsig_multiflux_xs_type3_cat4->SetFillStyle(1001);     
            
            
            mnvh1d_mcsig_multiflux_xs_type3_cat5->SetLineColor(kCyan + 1);
            mnvh1d_mcsig_multiflux_xs_type3_cat5->SetFillColor(kCyan + 1);
            mnvh1d_mcsig_multiflux_xs_type3_cat5->SetFillStyle(1001);     
            
        
        } else {

            mnvh1d_mcsig_multiflux_xs_type3_cat1->SetLineColor(17);
            mnvh1d_mcsig_multiflux_xs_type3_cat1->SetFillColor(17);
            mnvh1d_mcsig_multiflux_xs_type3_cat1->SetFillStyle(1001);
            
            mnvh1d_mcsig_multiflux_xs_type3_cat2->SetLineColor(17);
            mnvh1d_mcsig_multiflux_xs_type3_cat2->SetFillColor(17);
            mnvh1d_mcsig_multiflux_xs_type3_cat2->SetFillStyle(1001); // 3001        
            
            mnvh1d_mcsig_multiflux_xs_type3_cat3->SetLineColor(kGreen+2);
            mnvh1d_mcsig_multiflux_xs_type3_cat3->SetFillColor(kGreen+2);
            mnvh1d_mcsig_multiflux_xs_type3_cat3->SetFillStyle(1001); // 3003
                //mnvh1d_mcsig_multiflux_xs_type3_cat3->SetLineStyle(kDashed);
            
            mnvh1d_mcsig_multiflux_xs_type3_cat4->SetLineColor(kRed+2);
            mnvh1d_mcsig_multiflux_xs_type3_cat4->SetFillColor(kRed+2);
            mnvh1d_mcsig_multiflux_xs_type3_cat4->SetFillStyle(1001); // 3005     
                //mnvh1d_mcsig_multiflux_xs_type3_cat4->SetLineStyle(9); // 3005     
            
            mnvh1d_mcsig_multiflux_xs_type3_cat5->SetLineColor(kBlue);
            mnvh1d_mcsig_multiflux_xs_type3_cat5->SetFillColor(kBlue);
            mnvh1d_mcsig_multiflux_xs_type3_cat5->SetFillStyle(3010);     
            
        }
    
            //objarray_xs_multiflux_process_breakdown->Add(mnvh1d_mcsig_multiflux_xs_type3_cat1);
        objarray_xs_multiflux_process_breakdown->Add(mnvh1d_mcsig_multiflux_xs_type3_cat4);
        objarray_xs_multiflux_process_breakdown->Add(mnvh1d_mcsig_multiflux_xs_type3_cat3);
        objarray_xs_multiflux_process_breakdown->Add(mnvh1d_mcsig_multiflux_xs_type3_cat2);
            //objarray_xs_multiflux_process_breakdown->Add(mnvh1d_mcsig_multiflux_xs_type3_cat5);
        
    

        c = new TCanvas("c24");
        draw_data_mcstacked(plotter,
                            mnvh1d_tot_data_multiflux_xs,
                            objarray_xs_multiflux_process_breakdown,
                            1.0,
                            plot_attr_variable.legend_position.c_str(),
                                //"Data", 2, 1, 3001,
                            "Data", -1, 1, -1,
                            plot_attr_variable.xs_xlabel.c_str(),
                            plot_attr_variable.xs_ylabel.c_str(),
                            not_cov_area_norm,
                            is_enu);
        plotter->WritePreliminary("TL");
        write_hist_title(plotter);
        write_pot_exposure(plotter, recorded_pot);
        c->Print(Form("%s-data-mc-xs-multiflux-process-breakdown-pot%s", variable_name, __figure_format.c_str()));
        delete c;
    
            // area normalized
        c = new TCanvas("c25");
        draw_data_mcstacked(plotter,
                            mnvh1d_tot_data_multiflux_xs,
                            objarray_xs_multiflux_process_breakdown,
                            calc_area_ratio(mnvh1d_tot_data_multiflux_xs, mnvh1d_gnsig_multiflux_xs),
                            plot_attr_variable.legend_position.c_str(),
                                //"Data", 2, 1, 3001,
                            "Data", -1, 1, -1,
                            plot_attr_variable.xs_xlabel.c_str(),
                            plot_attr_variable.xs_ylabel.c_str(),
                            cov_area_norm,
                            is_enu);
        plotter->WritePreliminary("TL");
        write_hist_title(plotter);
        write_area_exposure(plotter, recorded_pot);
        c->Print(Form("%s-data-mc-xs-multiflux-process-breakdown-area%s", variable_name, __figure_format.c_str()));
        delete c;
        
        gStyle->SetNdivisions(old_xndivision, "X");
        gStyle->SetNdivisions(old_yndivision, "Y");
        
        plotter->mc_line_width = old_line_width;
        if (is_enu) plotter->draw_normalized_to_bin_width = true;


        c = new TCanvas("c26");
        draw_data_mcbreakdown(plotter,
                              get_h1d_with_total_pot_err(mnvh1d_tot_data_multiflux_xs,     is_enu),
                              mnvh1d_adapter_ptr(mnvh1d_gnsig_multiflux_xs,            is_enu),
                              mnvh1d_adapter_ptr(mnvh1d_mcsig_multiflux_xs_type3_cat1, is_enu),
                              mnvh1d_adapter_ptr(mnvh1d_mcsig_multiflux_xs_type3_cat2, is_enu),
                              mnvh1d_adapter_ptr(mnvh1d_mcsig_multiflux_xs_type3_cat3, is_enu),
                              mnvh1d_adapter_ptr(mnvh1d_mcsig_multiflux_xs_type3_cat4, is_enu),
                              mnvh1d_adapter_ptr(mnvh1d_mcsig_multiflux_xs_type3_cat5, is_enu),
                              1.0,
                              plot_attr_variable.legend_position.c_str(), false);
        plotter->WritePreliminary("TL");
        write_hist_title(plotter);
        write_pot_exposure(plotter, recorded_pot);
        c->Print(Form("%s-data-mc-xs-multiflux-component-type3-pot%s",variable_name, __figure_format.c_str()));
        delete c;
        
        c = new TCanvas("c27");
        draw_data_mcbreakdown(plotter,
                              get_h1d_with_total_area_err(mnvh1d_tot_data_multiflux_xs, is_enu),
                              mnvh1d_adapter_ptr(mnvh1d_gnsig_multiflux_xs,             is_enu),
                              mnvh1d_adapter_ptr(mnvh1d_mcsig_multiflux_xs_type3_cat1,  is_enu),
                              mnvh1d_adapter_ptr(mnvh1d_mcsig_multiflux_xs_type3_cat2,  is_enu),
                              mnvh1d_adapter_ptr(mnvh1d_mcsig_multiflux_xs_type3_cat3,  is_enu),
                              mnvh1d_adapter_ptr(mnvh1d_mcsig_multiflux_xs_type3_cat4,  is_enu),
                              mnvh1d_adapter_ptr(mnvh1d_mcsig_multiflux_xs_type3_cat5,  is_enu),
                              calc_area_ratio(mnvh1d_tot_data_multiflux_xs, mnvh1d_gnsig_multiflux_xs),
                              plot_attr_variable.legend_position.c_str(), false);
        plotter->WritePreliminary("TL");
        write_hist_title(plotter);
        write_area_exposure(plotter, recorded_pot);
        c->Print(Form("%s-data-mc-xs-multiflux-component-type3-area%s",variable_name, __figure_format.c_str()));
        delete c;
        
    } // legacy_code

    
    
    std::ofstream fstream_chi2(Form("chi2-%s.txt", variable_name));
    
    std::vector<std::string> vertNames = mnvh1d_gnsig_multiflux_xs->GetVertErrorBandNames();
    for (std::vector<std::string>::iterator s = vertNames.begin();
         s != vertNames.end(); ++s) {

        if (!__draw_genie_band) break;
    
        std::string name = *s;
        
        if (name.find("GENIE") == std::string::npos) continue;
        
        std::cout << "\tknob: " << name << std::endl;   
        
        std::vector<TH1D*> variation_hists = mnvh1d_gnsig_multiflux_xs->GetVertErrorBand(name)->GetHists();
    
        assert(variation_hists.size() == 2);
    
        TH1D  h_cv      = mnvh1d_gnsig_multiflux_xs->GetCVHistoWithStatError();
        TH1D* h_shiftdn = variation_hists.at(0);
        TH1D* h_shiftup = variation_hists.at(1);
        
        
        MnvH1D* mnvh1d_shiftdn = new MnvH1D(*h_shiftdn);
        MnvH1D* mnvh1d_shiftup = new MnvH1D(*h_shiftup);

            // remove GENIE_ from name
        name.replace(0,6,"");
        
        mnvh1d_shiftup->SetTitle(Form("%s +1#sigma", name.c_str()));
        mnvh1d_shiftdn->SetTitle(Form("%s -1#sigma", name.c_str()));
        
        mnvh1d_shiftdn->SetNormBinWidth(width);
        mnvh1d_shiftup->SetNormBinWidth(width);
        
        c = new TCanvas("c27");
        draw_data_mcvariation(plotter,
                              get_h1d_with_total_pot_err(mnvh1d_tot_data_multiflux_xs, is_enu),
                              get_h1d_with_stat_err(mnvh1d_tot_data_multiflux_xs,      is_enu),
                              mnvh1d_adapter_ptr(mnvh1d_gnsig_multiflux_xs,            is_enu),
                              mnvh1d_adapter_ptr(mnvh1d_shiftup,                       is_enu),
                              mnvh1d_adapter_ptr(mnvh1d_shiftdn,                       is_enu),
                              NULL, NULL,
                              1.0, 1.0, 1.0, 1.0, 1.0,
                              kSolid, kDashed, kDashed, kSolid, kSolid,
                              kRed, kWhite, kBlue, kBlue, kBlue,
                              plot_attr_variable.legend_position.c_str(), false, is_enu);
        plotter->WritePreliminary("TL");
        write_hist_title(plotter);
        write_pot_exposure(plotter, recorded_pot);
        c->Print(Form("optimize/%s-data-mc-xs-multiflux-pot-%s%s",variable_name, name.c_str(), __figure_format.c_str()));
        delete c;
  
        c = new TCanvas("c29");
        draw_data_mcvariation(plotter,
                              get_h1d_with_total_area_err(mnvh1d_tot_data_multiflux_xs,  is_enu),
                              get_h1d_with_stat_err(mnvh1d_tot_data_multiflux_xs,        is_enu),
                              mnvh1d_adapter_ptr(mnvh1d_gnsig_multiflux_xs,              is_enu),
                              mnvh1d_adapter_ptr(mnvh1d_shiftup,                         is_enu),
                              mnvh1d_adapter_ptr(mnvh1d_shiftdn,                         is_enu),
                              NULL, NULL,
                              calc_area_ratio(mnvh1d_tot_data_multiflux_xs, mnvh1d_gnsig_multiflux_xs),
                              calc_area_ratio(mnvh1d_tot_data_multiflux_xs, mnvh1d_shiftup),
                              calc_area_ratio(mnvh1d_tot_data_multiflux_xs, mnvh1d_shiftdn),
                              1.0, 1.0,
                              kSolid, kDashed, kDashed, kSolid, kSolid,
                              kRed, kWhite, kBlue, kBlue, kBlue,
                              plot_attr_variable.legend_position.c_str(), false, is_enu);
        plotter->WritePreliminary("TL");
        write_hist_title(plotter);
        write_area_exposure(plotter, recorded_pot);
        c->Print(Form("optimize/%s-data-mc-xs-multiflux-area-%s%s",variable_name, name.c_str(), __figure_format.c_str()));
        delete c;
    
        // calculate chi2
        bool useDataErrorMatrix = true;
        bool shapeOnly = true;
        bool absolute = !shapeOnly;
        int  ndof = 0;
   
        // absolute
        double abs_chi2_cv = plotter->Chi2DataMC(mnvh1d_tot_data_multiflux_xs,
                                             mnvh1d_gnsig_multiflux_xs,
                                             ndof,
                                             1.0,
                                             useDataErrorMatrix,
                                             absolute);
   
        double abs_chi2_up = plotter->Chi2DataMC(mnvh1d_tot_data_multiflux_xs,
                                             mnvh1d_shiftup,
                                             ndof,
                                             1.0,
                                             useDataErrorMatrix,
                                             absolute);
   
        double abs_chi2_dn = plotter->Chi2DataMC(mnvh1d_tot_data_multiflux_xs,
                                             mnvh1d_shiftdn,
                                             ndof,
                                             1.0,
                                             useDataErrorMatrix,
                                             absolute);
   
        // shape only
        double shape_chi2_cv = plotter->Chi2DataMC(mnvh1d_tot_data_multiflux_xs,
                                             mnvh1d_gnsig_multiflux_xs,
                                             ndof,
                                             calc_area_ratio(mnvh1d_tot_data_multiflux_xs, mnvh1d_gnsig_multiflux_xs),
                                             useDataErrorMatrix,
                                             shapeOnly);
   
        double shape_chi2_up = plotter->Chi2DataMC(mnvh1d_tot_data_multiflux_xs,
                                             mnvh1d_shiftup,
                                             ndof,
                                             calc_area_ratio(mnvh1d_tot_data_multiflux_xs, mnvh1d_shiftup),
                                             useDataErrorMatrix,
                                             shapeOnly);
   
        double shape_chi2_dn = plotter->Chi2DataMC(mnvh1d_tot_data_multiflux_xs,
                                             mnvh1d_shiftdn,
                                             ndof,
                                             calc_area_ratio(mnvh1d_tot_data_multiflux_xs, mnvh1d_shiftdn),
                                             useDataErrorMatrix,
                                             shapeOnly);
   
        fstream_chi2 << ndof << ","
                     << abs_chi2_dn   << "," << abs_chi2_cv   << "," << abs_chi2_up    << ","
                     << shape_chi2_dn << "," << shape_chi2_cv << "," << shape_chi2_up  << ","
                     << name
                     << std::endl;
        
    } // loop over knobs
    
    fstream_chi2.close();

    
        //------------------------------------------------------------------------
        // DRAW UNCERTAINTIES
        //------------------------------------------------------------------------
        //plotter = new MnvPlotter(kCCNuPionIncStyle);
    customize_plotter_uncertainty(plotter);
        // draw cross sections uncertainties
        // minerva5 data only
    plotter->axis_maximum = 0.5;
    plotter->legend_n_columns = 2;
    
    double epsilon = 0.001;






  c = new TCanvas("c30");
    mnvh1d_minerva_data_bkg->GetXaxis()->SetTitle(plot_attr_variable.xs_xlabel.c_str());
    mnvh1d_minerva_data_bkg->GetYaxis()->CenterTitle();
    plotter->DrawErrorSummary(mnvh1d_minerva_data_bkg,
                              plot_attr_variable.legend_position.c_str(),
                              true,
                              true,
                              epsilon,
                              not_cov_area_norm,
                              "",
                              as_frac,
                              plot_attr_variable.xlabel.c_str()); 
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    c->Print(Form("%s-mnv-data-bkg-pot-error%s",variable_name, __figure_format.c_str()));
    delete c;


    c = new TCanvas("c31");
    plotter->DrawErrorSummary(mnvh1d_minerva_data_bkg,
                              plot_attr_variable.legend_position.c_str(),
                              true,
                              true,
                              epsilon,
                              cov_area_norm,
                              "",
                              as_frac,
                              plot_attr_variable.xlabel.c_str()); 
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    c->Print(Form("%s-mnv-data-bkg-area-error%s",variable_name, __figure_format.c_str()));
    delete c;
    
        //---------------------------------------------------------------------


    c = new TCanvas("c32");
    mnvh1d_minerva_data_bkgsub->GetXaxis()->SetTitle(plot_attr_variable.xs_xlabel.c_str());
    mnvh1d_minerva_data_bkgsub->GetYaxis()->CenterTitle();
    plotter->DrawErrorSummary(mnvh1d_minerva_data_bkgsub,
                              plot_attr_variable.legend_position.c_str(),
                              true,
                              true,
                              epsilon,
                              not_cov_area_norm,
                              "",
                              as_frac,
                              plot_attr_variable.xlabel.c_str()); 
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    c->Print(Form("%s-mnv-data-bkgsub-pot-error%s",variable_name, __figure_format.c_str()));
    delete c;


    c = new TCanvas("c33");
    plotter->DrawErrorSummary(mnvh1d_minerva_data_bkgsub,
                              plot_attr_variable.legend_position.c_str(),
                              true,
                              true,
                              epsilon,
                              cov_area_norm,
                              "",
                              as_frac,
                              plot_attr_variable.xlabel.c_str()); 
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    c->Print(Form("%s-mnv-data-bkgsub-area-error%s",variable_name, __figure_format.c_str()));
    delete c;
    
        //-------------------------------------------------------------------------
    c = new TCanvas("c34");
    mnvh1d_minerva_data_unfolded->GetXaxis()->SetTitle(plot_attr_variable.xs_xlabel.c_str());
    mnvh1d_minerva_data_unfolded->GetYaxis()->CenterTitle();
    plotter->DrawErrorSummary(mnvh1d_minerva_data_unfolded,
                              plot_attr_variable.legend_position.c_str(),
                              true,
                              true,
                              epsilon,
                              not_cov_area_norm,
                              "",
                              as_frac,
                              plot_attr_variable.xlabel.c_str()); 
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    c->Print(Form("%s-mnv-data-unfolded-pot-error%s",variable_name, __figure_format.c_str()));
    delete c;


    c = new TCanvas("c35");
    plotter->DrawErrorSummary(mnvh1d_minerva_data_unfolded,
                              plot_attr_variable.legend_position.c_str(),
                              true,
                              true,
                              epsilon,
                              cov_area_norm,
                              "",
                              as_frac,
                              plot_attr_variable.xlabel.c_str()); 
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    c->Print(Form("%s-mnv-data-unfolded-area-error%s",variable_name, __figure_format.c_str()));
    delete c;

        //-----------------------------------------------------------------------
    c = new TCanvas("c36");
    mnvh1d_minerva_mcsig_produced->GetXaxis()->SetTitle(plot_attr_variable.xs_xlabel.c_str());
    mnvh1d_minerva_mcsig_produced->GetYaxis()->CenterTitle();
    plotter->DrawErrorSummary(mnvh1d_minerva_mcsig_produced,
                              plot_attr_variable.legend_position.c_str(),
                              true,
                              true,
                              epsilon,
                              not_cov_area_norm,
                              "",
                              as_frac,
                              plot_attr_variable.xlabel.c_str()); 
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    c->Print(Form("%s-mnv-mc-produced-pot-error%s",variable_name, __figure_format.c_str()));
    delete c;


    c = new TCanvas("c37");
    plotter->DrawErrorSummary(mnvh1d_minerva_mcsig_produced,
                              plot_attr_variable.legend_position.c_str(),
                              true,
                              true,
                              epsilon,
                              cov_area_norm,
                              "",
                              as_frac,
                              plot_attr_variable.xlabel.c_str()); 
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    c->Print(Form("%s-mnv-mc-produced-area-error%s",variable_name, __figure_format.c_str()));
    delete c;

        //------------------------------------------------------------------------
    c = new TCanvas("c38");
    mnvh1d_minerva_mcsig_selected->GetXaxis()->SetTitle(plot_attr_variable.xs_xlabel.c_str());
    mnvh1d_minerva_mcsig_selected->GetYaxis()->CenterTitle();
    plotter->DrawErrorSummary(mnvh1d_minerva_mcsig_selected,
                              plot_attr_variable.legend_position.c_str(),
                              true,
                              true,
                              epsilon,
                              not_cov_area_norm,
                              "",
                              as_frac,
                              plot_attr_variable.xlabel.c_str()); 
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    c->Print(Form("%s-mnv-mc-selected-pot-error%s",variable_name, __figure_format.c_str()));
    delete c;


    c = new TCanvas("c39");
    plotter->DrawErrorSummary(mnvh1d_minerva_mcsig_selected,
                              plot_attr_variable.legend_position.c_str(),
                              true,
                              true,
                              epsilon,
                              cov_area_norm,
                              "",
                              as_frac,
                              plot_attr_variable.xlabel.c_str()); 
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    c->Print(Form("%s-mnv-mc-selected-area-error%s",variable_name, __figure_format.c_str()));
    delete c;
        //-----------------------------------------------------------------------
    c = new TCanvas("c40");
    mnvh1d_minerva_efficiency->GetXaxis()->SetTitle(plot_attr_variable.xs_xlabel.c_str());
    mnvh1d_minerva_efficiency->GetYaxis()->CenterTitle();
    plotter->DrawErrorSummary(mnvh1d_minerva_efficiency,
                              plot_attr_variable.legend_position.c_str(),
                              true,
                              true,
                              epsilon,
                              not_cov_area_norm,
                              "",
                              as_frac,
                              plot_attr_variable.xlabel.c_str()); 
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    c->Print(Form("%s-mnv-mc-efficiency-pot-error%s",variable_name, __figure_format.c_str()));
    delete c;


    c = new TCanvas("c41");
    plotter->DrawErrorSummary(mnvh1d_minerva_efficiency,
                              plot_attr_variable.legend_position.c_str(),
                              true,
                              true,
                              epsilon,
                              cov_area_norm,
                              "",
                              as_frac,
                              plot_attr_variable.xlabel.c_str()); 
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    c->Print(Form("%s-mnv-mc-efficiency-area-error%s",variable_name, __figure_format.c_str()));
    delete c;

        //------------------------------------------------------------------------
    c = new TCanvas("c42");
    mnvh1d_minerva_data_effcorrected->GetXaxis()->SetTitle(plot_attr_variable.xs_xlabel.c_str());
    mnvh1d_minerva_data_effcorrected->GetYaxis()->CenterTitle();
    plotter->DrawErrorSummary(mnvh1d_minerva_data_effcorrected,
                              plot_attr_variable.legend_position.c_str(),
                              true,
                              true,
                              epsilon,
                              not_cov_area_norm,
                              "",
                              as_frac,
                              plot_attr_variable.xlabel.c_str()); 
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    c->Print(Form("%s-mnv-data-effcorrected-pot-error%s",variable_name, __figure_format.c_str()));
    delete c;


    c = new TCanvas("c43");
    plotter->DrawErrorSummary(mnvh1d_minerva_data_effcorrected,
                              plot_attr_variable.legend_position.c_str(),
                              true,
                              true,
                              epsilon,
                              cov_area_norm,
                              "",
                              as_frac,
                              plot_attr_variable.xlabel.c_str()); 
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    c->Print(Form("%s-mnv-data-effcorrected-area-error%s",variable_name, __figure_format.c_str()));
    delete c;

    
        //-------------------------------------------------------------------------
    c = new TCanvas("c44");
    mnvh1d_minerva_data_multiflux_xs->GetXaxis()->SetTitle(plot_attr_variable.xs_xlabel.c_str());
    mnvh1d_minerva_data_multiflux_xs->GetYaxis()->CenterTitle();
    plotter->DrawErrorSummary(mnvh1d_minerva_data_multiflux_xs,
                              plot_attr_variable.legend_position.c_str(),
                              true,
                              true,
                              epsilon,
                              not_cov_area_norm,
                              "",
                              as_frac,
                              plot_attr_variable.xlabel.c_str()); 
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
        //write_pot_exposure(plotter, data_minerva_pot);
    c->Print(Form("%s-mnv-data-mc-xs-multiflux-pot-error%s",variable_name, __figure_format.c_str()));
    delete c;


    c = new TCanvas("c45");
    plotter->DrawErrorSummary(mnvh1d_minerva_data_multiflux_xs,
                              plot_attr_variable.legend_position.c_str(),
                              true,
                              true,
                              epsilon,
                              cov_area_norm,
                              "",
                              as_frac,
                              plot_attr_variable.xlabel.c_str()); 
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
        //write_area_exposure(plotter, data_minerva_pot);
    c->Print(Form("%s-mnv-data-mc-xs-multiflux-area-error%s",variable_name, __figure_format.c_str()));
    delete c;

    
        // downstream data only
    mnvh1d_downstream_data_multiflux_xs->GetXaxis()->SetTitle(plot_attr_variable.xs_xlabel.c_str());
    mnvh1d_downstream_data_multiflux_xs->GetYaxis()->CenterTitle();
    c = new TCanvas("c46");
    plotter->DrawErrorSummary(mnvh1d_downstream_data_multiflux_xs,
                              plot_attr_variable.legend_position.c_str(),
                              true,
                              true,
                              epsilon,
                              not_cov_area_norm,
                              "",
                              as_frac,
                              plot_attr_variable.xlabel.c_str()); 
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
        //write_pot_exposure(plotter, equivalent_pot);
    c->Print(Form("%s-frz-data-mc-xs-multiflux-pot-error%s",variable_name, __figure_format.c_str()));
    delete c;

    c = new TCanvas("c35");
    plotter->DrawErrorSummary(mnvh1d_downstream_data_multiflux_xs,
                              plot_attr_variable.legend_position.c_str(),
                              true,
                              true,
                              epsilon,
                              cov_area_norm,
                              "",
                              as_frac,
                              plot_attr_variable.xlabel.c_str()); 
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
        //write_area_exposure(plotter, equivalent_pot);
    c->Print(Form("%s-frz-data-mc-xs-multiflux-area-error%s",variable_name, __figure_format.c_str()));
    delete c;


        //--------------------------- draw and print the combined xsec
    if (__variable == "kinetic") mnvh1d_tot_data_multiflux_xs = remove_right_bin(mnvh1d_tot_data_multiflux_xs);

    
        // minerva5 + downstream
    c = new TCanvas("c47");
    plotter->DrawErrorSummary(mnvh1d_tot_data_multiflux_xs,
                              plot_attr_variable.legend_position.c_str(),
                              true,
                              true,
                              epsilon,
                              not_cov_area_norm,
                              "",
                              as_frac,
                              plot_attr_variable.xlabel.c_str()); 
        //plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    write_label(plotter);
        //write_pot_exposure(plotter, tot_pot);
    c->Print(Form("%s-tot-data-mc-xs-multiflux-pot-error%s",variable_name, __figure_format.c_str()));
    delete c;

    
    c = new TCanvas("c48");
    plotter->DrawErrorSummary(mnvh1d_tot_data_multiflux_xs,
                              plot_attr_variable.legend_position.c_str(),
                              true,
                              true,
                              epsilon,
                              cov_area_norm,
                              "",
                              as_frac,
                              plot_attr_variable.xlabel.c_str()); 
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
        //write_area_exposure(plotter, tot_pot);
    c->Print(Form("%s-tot-data-mc-xs-multiflux-area-error%s",variable_name, __figure_format.c_str()));
    delete c;


#define draw_flux_group    
#ifdef draw_flux_group
    std::pair<std::map<std::string, TH1D*>, std::map<std::string, TH1D*> > result
        = add_group_error(mnvh1d_tot_data_multiflux_xs);
    std::map<std::string, TH1D*> error_group_map = result.first;

    
    c = new TCanvas("c49");

    std::string groupName("Flux");
    error_group_map[groupName]->Draw("HIST");
    error_group_map[groupName]->SetLineWidth(2);
    error_group_map[groupName]->SetLineColor(kYellow-3);
    error_group_map[groupName]->GetYaxis()->SetRangeUser(0.0001, 0.3);
    error_group_map[groupName]->GetYaxis()->SetTitle("Fractional Uncertainty");
    error_group_map[groupName]->GetYaxis()->CenterTitle();
    error_group_map[groupName]->GetYaxis()->SetLabelSize(0.05);
    error_group_map[groupName]->GetYaxis()->SetTitleSize(0.06);
    error_group_map[groupName]->GetYaxis()->SetTitleOffset(1.1);
    error_group_map[groupName]->GetXaxis()->SetLabelSize(0.05);
    error_group_map[groupName]->GetXaxis()->SetTitleSize(0.06);
    error_group_map[groupName]->GetXaxis()->SetTitleOffset(1.1);
    error_group_map[groupName]->GetXaxis()->SetNdivisions(plot_attr_variable.xndivision);
    
    plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    write_pot_exposure(plotter, recorded_pot);
    c->Print(Form("%s-tot-data-mc-xs-multiflux-pot-flux-error%s",variable_name, __figure_format.c_str()));
    delete c;

#endif
    draw_systematics drawObj(mnvh1d_tot_data_multiflux_xs);
    drawObj.draw_all_group(Form("%s-data-xsec-uncertainty", variable_name));
        //drawObj.draw_one("test_prefix", "Flux");


    
    std::cout << "Correlation matrix " << std::endl;
    TMatrixD corr_mat_pot = mnvh1d_tot_data_multiflux_xs->GetTotalCorrelationMatrix(not_cov_area_norm, 
                                                                                    include_stat);  
        //corr_mat.Print();
    print_matrix_stdcout(corr_mat_pot, mnvh1d_tot_data_multiflux_xs);
    

    std::cout << "Error matrix " << std::endl;
    TMatrixD error_mat_pot = mnvh1d_tot_data_multiflux_xs->GetTotalErrorMatrix(include_stat, 
                                                                               not_as_frac,  
                                                                               not_cov_area_norm); 
        //error_mat.Print();
    print_matrix_stdcout(error_mat_pot, mnvh1d_tot_data_multiflux_xs);

    TMatrixD corr_mat_area = mnvh1d_tot_data_multiflux_xs->GetTotalCorrelationMatrix(cov_area_norm, 
                                                                                     include_stat);  
  
    std::ofstream correlation_matrix_pot_ofs(Form("%s-correlation-matrix-pot.txt", variable_name));
    print_matrix_ofstream(corr_mat_pot,
                          mnvh1d_tot_data_multiflux_xs,
                          correlation_matrix_pot_ofs);
    
    std::ofstream correlation_matrix_area_ofs(Form("%s-correlation-matrix-area.txt", variable_name));
    print_matrix_ofstream(corr_mat_area,
                          mnvh1d_tot_data_multiflux_xs,
                          correlation_matrix_area_ofs);
    
    
        //std::cout << "Check divide" << std::endl;
        //print_correlation_matrix_from_covariance(error_mat_pot, mnvh1d_tot_data_multiflux_xs);

    std::cout << "xsec-data for variable: " << __variable << std::endl;
    print_error_group(mnvh1d_tot_data_multiflux_xs);

    
    c = new TCanvas("c47");
    plotter->DrawErrorSummary(mnvh1d_tot_data_multiflux_xs,
                              plot_attr_variable.legend_position.c_str(),
                              true,
                              true,
                              epsilon,
                              not_cov_area_norm,
                              "",
                              as_frac,
                              plot_attr_variable.xlabel.c_str()); 
        //plotter->WritePreliminary("TL");
    write_hist_title(plotter);
    write_label(plotter);
        //write_pot_exposure(plotter, tot_pot);
    c->Print(Form("%s-tot-data-mc-xs-multiflux-pot-error-v2%s",variable_name, __figure_format.c_str()));
    delete c;

    
    if (!__export_xsec_data.empty()) {
        TFile f(__export_xsec_data.c_str(), "update");

        mnvh1d_tot_data_multiflux_xs->SetName(Form("%s-xsec-data", variable_name));
        mnvh1d_tot_data_multiflux_xs->Write();

        f.Close();
    }
}


