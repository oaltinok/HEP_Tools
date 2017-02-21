#include <iostream>
#include <iomanip>
#include <map>
#include <string>
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <sys/stat.h>
#include <getopt.h>


#include <TFile.h>
#include <TString.h>
#include <Cintex/Cintex.h>
#include <TMath.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TPaveText.h>

#include <PlotUtils/MnvH1D.h>
#include <PlotUtils/MnvPlotter.h>

#include "plot_attr.h"
#include "draw_data_and_models.h"
#include "draw_data_and_models_ratio.h"
#include "read_data_and_convert_to_graph.h"
#include "get_h1d.h"


using std::setw;
using std::sqrt;

using PlotUtils::MnvH1D;
using PlotUtils::MnvVertErrorBand;
using PlotUtils::MnvPlotter;



namespace {

  
        // Must be called AFTER SetRootEnv to override its behaviors
    void customize_plotter(MnvPlotter* plotter)
    {
        gStyle->SetCanvasDefW(900);
        gStyle->SetCanvasDefH(675); 
        gStyle->SetPadRightMargin(0.05);
        
        gStyle->SetStripDecimals(false);
        plotter->legend_text_size  = 0.045;
        plotter->legend_text_font = 42; 
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
            //gStyle->SetFrameLineWidth(1);
            //gStyle->SetEndErrorSize(4);
        
    }

    void write_hist_title(MnvPlotter* plotter)
    {
        const char* plot_title = "";
        plotter->AddHistoTitle(plot_title,0.05);
    }

    void write_pot_exposure(MnvPlotter* plotter, double pot)
    {
        TPaveText* text = new TPaveText(0.18, 0.80, 0.54, 0.90, "NDC");
        text->AddText("#scale[1.15]{b)  #bar{#nu}_{#mu} + CH #rightarrow #mu^{+} + #pi^{0} + X}");
        text->AddText("#scale[0.85]{POT Normalized}");
        text->SetTextColor(kBlue);
        text->SetFillColor(kWhite);
        text->SetTextFont(42);
        text->Draw("SAME");
    }

    TH1D* remove_last_bin(const TH1D* h1d)
    {
        unsigned int N = h1d->GetNbinsX();
        const double* edges = h1d->GetXaxis()->GetXbins()->GetArray();
        double new_edges[N];
        for (unsigned int i = 0; i <= N; ++i) {
            new_edges[i] = edges[i];
            printf("%5d %8.2f\n", i, edges[i]);
        }
        
        TH1D* result = new TH1D("result", "Trimmed histogram", N-1, new_edges);
        for (unsigned int bin = 0; bin <= N; ++bin) {
            double val = h1d->GetBinContent(bin);
            double err = h1d->GetBinError(bin);
            
            result->SetBinContent(bin, val);
            result->SetBinError(bin, err);
        }

        result->GetXaxis()->ImportAttributes(h1d->GetXaxis());
        result->GetYaxis()->ImportAttributes(h1d->GetYaxis());
        
        return result;
    }
}


// As the name implies, compare data xsec to different generator predictions:
// GENIE with and without FSI, NuWro, and NEUT
int main(int argc, char* argv[])
{

    ROOT::Cintex::Cintex::Enable();

    TH1::AddDirectory(false);

    std::string __xsec_data;
    std::string __xsec_genie_fsi;
    std::string __xsec_genie_nofsi;
    std::string __xsec_nuwro;
    std::string __xsec_neut;
     
    std::string __variable;
    std::string __figure_format(".png");

    
    for (;;) {

        int option_index = 0;
        static struct option long_options[] = {
            {"variable",      required_argument, 0,  0 },
            {"figure-format", required_argument, 0,  0 },
            {"data",          required_argument, 0,  0 },
            {"genie-fsi",     required_argument, 0,  0 },
            {"genie-nofsi",   required_argument, 0,  0 },
            {"nuwro",         required_argument, 0,  0 },
            {"neut",          required_argument, 0,  0 },
            {0,               0,                 0,  0 }
        };
        
        int c = getopt_long(argc, argv, "",
                            long_options, &option_index);
        
        if (c != 0 ) break;
        
        std::string opt_name(long_options[option_index].name);

        
        if (opt_name == "figure-format")     __figure_format = optarg;
        else if (opt_name == "variable")     __variable = optarg;
        else if (opt_name == "data")         __xsec_data  = optarg;
        else if (opt_name == "genie-fsi")    __xsec_genie_fsi  = optarg;
        else if (opt_name == "genie-nofsi")  __xsec_genie_nofsi  = optarg;
        else if (opt_name == "nuwro")        __xsec_nuwro  = optarg;
        else if (opt_name == "neut")         __xsec_neut  = optarg;
        else {}
        
    }
    

    std::cout << "__figure_format     "       << __figure_format     << std::endl;
    std::cout << "__variable          "       << __variable          << std::endl;
    std::cout << "__xsec_data         "       << __xsec_data         << std::endl;
    std::cout << "__xsec_genie_fsi    "       << __xsec_genie_fsi    << std::endl;
    std::cout << "__xsec_genie_nofsi  "       << __xsec_genie_nofsi  << std::endl;
    std::cout << "__xsec_nuwro        "       << __xsec_nuwro        << std::endl;
    std::cout << "__xsec_neut         "       << __xsec_neut         << std::endl;
    
        // check if file actually exists;
    std::map<std::string, std::string> filename_map;
    filename_map.insert(std::make_pair("__xsec_data",        __xsec_data));
    filename_map.insert(std::make_pair("__xsec_genie_fsi",   __xsec_genie_fsi));
    filename_map.insert(std::make_pair("__xsec_genie_nofsi", __xsec_genie_nofsi));
    filename_map.insert(std::make_pair("__xsec_nuwro",       __xsec_nuwro));
    filename_map.insert(std::make_pair("__xsec_neut",        __xsec_neut));
    
    
    for (std::map<std::string, std::string>::iterator n = filename_map.begin();
         n != filename_map.end(); ++n) {
        std::string name = n->second;
        std::string args = n->first;
        if (name.empty()) {
            std::cerr << "args not set: " << args << std::endl;
            continue;
        }
        
        struct stat buffer;   
        if (stat(name.c_str(), &buffer) != 0) {
            std::cerr << "File not found: " << name << std::endl;
            exit(1);
        }
        
    }

    MnvPlotter* plotter = new MnvPlotter();
    plotter->SetRootEnv();
    customize_plotter(plotter);
    const double data_minerva_pot  =  1.061e20;
    const double data_downstream_pot =  0.9458e20; 
    const double recorded_pot = data_minerva_pot + data_downstream_pot;
        
    std::map<std::string, plot_attr> plot_attr_map = create_plot_attr_map();
    plot_attr plot_attr_variable = plot_attr_map[__variable];    

    
        //--------------------------------------------------------------------------
    TFile* data_file = new TFile(__xsec_data.c_str(), "READ");
    MnvH1D* mnvh1d_data = static_cast<MnvH1D*>(data_file->Get(Form("%s-xsec-data", __variable.c_str())));
    assert(mnvh1d_data);

    TH1D* data_with_total_err(NULL);
    TH1D* data_with_stat_err(NULL);
    if (__variable != "enu") {
        data_with_total_err = get_var_h1d_with_total_pot_err(mnvh1d_data);
        data_with_stat_err  = get_var_h1d_with_stat_err(mnvh1d_data);
        
    } else {
        data_with_total_err = get_enu_h1d_with_total_pot_err(mnvh1d_data);
        data_with_stat_err  = get_enu_h1d_with_stat_err(mnvh1d_data);
    }

    if (__variable != "pimom" && __variable != "theta") {
        data_with_total_err = remove_last_bin(data_with_total_err);
        data_with_stat_err = remove_last_bin(data_with_stat_err);
    }
    
    data_with_total_err->SetMarkerStyle(plotter->data_marker);
    data_with_total_err->SetMarkerSize(plotter->data_marker_size);
    data_with_total_err->SetMarkerColor(plotter->data_color);
    data_with_total_err->SetLineWidth(plotter->data_line_width);
    data_with_total_err->SetLineStyle(plotter->data_line_style);
    data_with_total_err->SetLineColor(plotter->data_color);
    data_with_stat_err->SetLineWidth(plotter->data_line_width);
    data_with_stat_err->SetLineStyle(plotter->data_line_style);
    data_with_stat_err->SetLineColor(plotter->data_color);

    data_with_total_err->GetXaxis()->SetTitle(plot_attr_variable.xs_xlabel.c_str());
    data_with_total_err->GetXaxis()->SetTitleFont(plotter->axis_title_font_x);
    data_with_total_err->GetXaxis()->SetTitleSize(plotter->axis_title_size_x);
    data_with_total_err->GetXaxis()->SetTitleOffset(plotter->axis_title_offset_x);
    data_with_total_err->GetXaxis()->SetLabelFont(plotter->axis_label_font);
    data_with_total_err->GetXaxis()->SetLabelSize(plotter->axis_label_size);
    data_with_total_err->GetXaxis()->SetNdivisions(plot_attr_variable.xndivision);
    data_with_total_err->GetXaxis()->CenterTitle(kTRUE);
    data_with_total_err->GetXaxis()->SetLimits(plot_attr_variable.xs_xmin, plot_attr_variable.xs_xmax);

    data_with_total_err->GetYaxis()->SetTitle(plot_attr_variable.xs_ylabel.c_str());
    data_with_total_err->GetYaxis()->SetTitleFont(plotter->axis_title_font_y);
    data_with_total_err->GetYaxis()->SetTitleSize(plotter->axis_title_size_y);
    data_with_total_err->GetYaxis()->SetTitleOffset(plotter->axis_title_offset_y);
    data_with_total_err->GetYaxis()->SetLabelFont(plotter->axis_label_font);
    data_with_total_err->GetYaxis()->SetLabelSize(plotter->axis_label_size);
    data_with_total_err->GetYaxis()->SetNdivisions(plot_attr_variable.yndivision);
    data_with_total_err->GetYaxis()->CenterTitle(kTRUE);
    data_with_total_err->GetYaxis()->SetTitleOffset(1.1);
    
    std::cout << "Title offset: " << data_with_total_err->GetYaxis()->GetTitleOffset() << std::endl;
    
    std::vector<TGraph*> models;
    
    if (!__xsec_genie_fsi.empty()) {
        TGraph* graph = read_data_and_convert_to_graph(__xsec_genie_fsi);
        graph->SetTitle("GENIE w/ FSI");
        graph->SetLineColor(kRed);
        graph->SetLineStyle(kSolid);
        graph->SetLineWidth(plotter->mc_line_width);
        models.push_back(graph);
    }
    
    if (!__xsec_genie_nofsi.empty()) {
        TGraph* graph = read_data_and_convert_to_graph(__xsec_genie_nofsi);
        graph->SetTitle("GENIE w/o FSI");
        graph->SetLineColor(kRed);
        graph->SetLineStyle(kDashed);
        graph->SetLineWidth(plotter->mc_line_width);
        models.push_back(graph);
    }

    if (!__xsec_nuwro.empty()) {
        TGraph* graph = read_data_and_convert_to_graph(__xsec_nuwro);
        graph->SetTitle("NuWro");
        graph->SetLineColor(kOrange-6);
        graph->SetLineStyle(7);
        graph->SetLineWidth(plotter->mc_line_width);
        models.push_back(graph);
    }

    if (!__xsec_neut.empty()) {
        TGraph* graph = read_data_and_convert_to_graph(__xsec_neut);
        graph->SetTitle("NEUT");
        graph->SetLineColor(kGreen+2);
        graph->SetLineStyle(5);
        graph->SetLineWidth(plotter->mc_line_width);
        models.push_back(graph);
    }



    
    TCanvas* c1 = new TCanvas("c1");
    plotter->axis_maximum = plot_attr_variable.xs_area_ymax;
    draw_data_and_models(plotter, data_with_total_err, data_with_stat_err, models, __variable == "enu");
    write_hist_title(plotter);
    write_pot_exposure(plotter, recorded_pot);
    c1->Print(Form("%s-data-models-pot%s", __variable.c_str(), __figure_format.c_str()));
    delete c1;


        //c1 = new TCanvas("c1");
        //draw_data_and_models_ratio(plotter, data_with_total_err, data_with_stat_err, models, __variable == "enu");
        //write_hist_title(plotter);
        //write_pot_exposure(plotter, recorded_pot);
        //c1->Print(Form("%s-data-models-pot-ratio%s", __variable.c_str(), __figure_format.c_str()));
        //delete c1;

        // Data range is settable
        // Things we might want to do with different models:
        //  - chi^2 comparison
        //  - Shape comparison with chi^2
    

}

