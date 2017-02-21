#ifndef draw_systematics_h
#define draw_systematics_h

#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <string>

#include <TPaveText.h>
#include <TLegend.h>
#include <TCanvas.h>
#include <TColor.h>

#include <PlotUtils/MnvH1D.h>

#include "mem_mgr.h"
#include "error_group.h"

using PlotUtils::MnvH1D;

namespace {

    void replace_char(std::string& s, char old_char, char new_char)
    {
        for (unsigned int i = 0; i < s.length(); ++i) {
            if (s[i] == old_char) s[i] = new_char;
        }
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

class draw_systematics {
  public:
    explicit draw_systematics(MnvH1D* h);

    void draw_one(const std::string& figure_prefix, const std::string& name = "Muon_momentum") const;
    void draw_individual(const std::string& figure_prefix) const;
    void draw_each_group(const std::string& figure_prefix) const;
    void draw_all_group(const std::string& figure_name) const;
    void draw_total(const std::string& figure_prefix) const;
    void draw_all(const std::string& figure_prefix) const;
    
  private:
    void do_draw(TH1D* h, const std::string& figure_name, const std::string& legend = "") const;
    
    std::vector<TH1D*> __err_hist_list;             // vector of all systematic uncertainties
    std::map<std::string, TH1D*> __err_hist_map;    // map of all systematic uncertainties indexed by errName
    std::map<std::string, TH1D*> __err_group_map;   // map of systematic uncertainty group indexed by groupName

    TH1D* __h1d_total_err;
    TH1D* __h1d_stat_err;
    
};

draw_systematics::draw_systematics(MnvH1D* mnvh1d)
:   __h1d_total_err(0),
    __h1d_stat_err(0)
{

    bool include_stat = true;
    bool as_frac = true;
    bool not_cov_area_norm = false;
    std::vector<std::string> vertNames = mnvh1d->GetVertErrorBandNames();
    for (std::vector<std::string>::iterator name = vertNames.begin();
         name != vertNames.end(); ++name) {
        
        std::string errName = *name;

        TH1D h1d_absolute_err =  mnvh1d->GetVertErrorBand(errName)->GetErrorBand(as_frac, not_cov_area_norm);
        TH1D* h1d_ptr_absolute_err = static_cast<TH1D*>(h1d_absolute_err.Clone(Form("tmp_abs_%s", errName.c_str())));

        mem_mgr::get().manage(h1d_ptr_absolute_err);

        __err_hist_map.insert(std::make_pair(errName,h1d_ptr_absolute_err));
        __err_hist_list.push_back(h1d_ptr_absolute_err);
    }


    std::vector<std::string> latNames = mnvh1d->GetLatErrorBandNames();
    for (std::vector<std::string>::iterator name = latNames.begin();
         name != latNames.end(); ++name) {
        
        std::string errName = *name;
        
        TH1D h1d_absolute_err =  mnvh1d->GetLatErrorBand(errName)->GetErrorBand(as_frac, not_cov_area_norm);
        TH1D* h1d_ptr_absolute_err = static_cast<TH1D*>(h1d_absolute_err.Clone(Form("tmp_abs_%s", errName.c_str())));

        mem_mgr::get().manage(h1d_ptr_absolute_err);

        __err_hist_map.insert(std::make_pair(errName,h1d_ptr_absolute_err));
        __err_hist_list.push_back(h1d_ptr_absolute_err);
    }
    
    std::pair<std::map<std::string, TH1D*>, std::map<std::string, TH1D*> > result = add_group_error(mnvh1d);
    __err_group_map = result.first;

    __h1d_total_err = static_cast<TH1D*>(mnvh1d->GetTotalError(include_stat, as_frac, not_cov_area_norm).Clone());
    __h1d_stat_err  = static_cast<TH1D*>(mnvh1d->GetStatError(as_frac).Clone());

    assert(__h1d_total_err);
    assert(__h1d_stat_err);
    

}

void draw_systematics::draw_one(const std::string& figure_prefix,
                                const std::string& name) const
{
    for (std::map<std::string, TH1D*>::const_iterator i = __err_hist_map.begin();
         i != __err_hist_map.end(); ++i) {
        std::string errName(i->first);
        
        if (errName == name) {
            TH1D* hist = i->second;
            std::string figure_name = figure_prefix + '_' + errName;
            do_draw(hist, figure_name, errName);
            std::cout << "Prefix: " << figure_prefix << std::endl;
            for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
                printf("\t%2d %8.2f\n", bin, hist->GetBinContent(bin));
            }
        }
    }
}

void draw_systematics::draw_individual(const std::string& figure_prefix) const
{
    for (std::map<std::string, TH1D*>::const_iterator i = __err_hist_map.begin();
         i != __err_hist_map.end(); ++i) {
        std::string errName(i->first);
        replace_char(errName, ' ', '_');
        TH1D* hist = i->second;
        std::string figure_name = figure_prefix + '_' + errName;
        do_draw(hist, figure_name, errName);
    }
}

void draw_systematics::draw_each_group(const std::string& figure_prefix) const
{

    for (std::map<std::string, TH1D*>::const_iterator i = __err_group_map.begin();
         i != __err_group_map.end(); ++i) {
        std::string groupName(i->first);
        replace_char(groupName, ' ', '_');
        TH1D* hist = i->second;
        std::string figure_name = figure_prefix + '_' + groupName;
        do_draw(hist, figure_name, groupName);
    }
    
}

void draw_systematics::draw_all_group(const std::string& figure_name) const 
{

    std::map<std::string, int> lines;
    lines["Detector"]            = 2;
    lines["Cross Section Model"] = 2;
    lines["FSI Model" ]          = 2;
    lines["Flux"]                = 2;
    lines["Other"]               = 2;

    std::map<std::string, int> styles;
    styles["Detector"]            = 3;
    styles["Cross Section Model"] = 5;
    styles["FSI Model" ]          = 9;
    styles["Flux"]                = 6;
    styles["Other"]               = 7;
    
    std::map<std::string, int> colors;
    colors["Detector"]            = kBlue+2;
    colors["Cross Section Model"] = kRed+2;
    colors["FSI Model" ]          = kMagenta+2;
    colors["Flux"]                = TColor::GetColor("#8b4513");
    colors["Other"]               = kGreen+3;

    std::vector<std::string> group_ordered_list;
    group_ordered_list.push_back("Detector");
    group_ordered_list.push_back("Cross Section Model");
    group_ordered_list.push_back("FSI Model");
    group_ordered_list.push_back("Flux");
    group_ordered_list.push_back("Other");
    
    std::map<std::string, std::string> keys;
    keys["Detector"]            = "Detector";
    keys["Cross Section Model"] = "X-Sec Model";
    keys["FSI Model" ]          = "FSI Model";
    keys["Flux"]                = "Flux";
    keys["Other"]               = "Other";

    
    TH1D* __h1d_total_err_tmp = remove_last_bin(__h1d_total_err);

    
    __h1d_total_err_tmp->SetLineWidth(2);
    __h1d_total_err_tmp->SetLineColor(kBlack);
    __h1d_total_err_tmp->SetLineStyle(kSolid);
    __h1d_total_err_tmp->GetYaxis()->SetRangeUser(0.0001, 0.5);
    __h1d_total_err_tmp->GetYaxis()->SetTitle("Fractional Uncertainty");
    __h1d_total_err_tmp->GetYaxis()->CenterTitle();
    __h1d_total_err_tmp->GetYaxis()->SetLabelSize(0.05);
    __h1d_total_err_tmp->GetYaxis()->SetTitleSize(0.06);
    __h1d_total_err_tmp->GetYaxis()->SetTitleOffset(1.1);
    __h1d_total_err_tmp->GetXaxis()->SetLabelSize(0.05);
    __h1d_total_err_tmp->GetXaxis()->SetTitleSize(0.06);
    __h1d_total_err_tmp->GetXaxis()->SetTitleOffset(1.1);
    __h1d_total_err_tmp->GetXaxis()->CenterTitle();
    __h1d_total_err_tmp->GetXaxis()->SetRangeUser(1.0, 9.5);

    __h1d_stat_err->SetLineWidth(2);
    __h1d_stat_err->SetLineColor(kBlack);
    __h1d_stat_err->SetLineStyle(kDashed);

    TCanvas* c = new TCanvas("c1");
    
    __h1d_total_err_tmp->Draw("HIST");
    __h1d_stat_err->Draw("SAME");
    

    TLegend* legend = new TLegend(0.52,0.74,0.94,0.9);

    legend->SetNColumns(2);
    legend->SetTextSize(0.038);
    legend->SetTextFont(42);
    legend->SetFillColor(0);
    legend->SetLineColor(0);
    legend->AddEntry(__h1d_total_err, "Total Error", "l");
    legend->AddEntry(__h1d_stat_err,  "Statistical", "l");
    std::cout << "Legend border:         " << legend->GetBorderSize() << std::endl;
    std::cout << "Legend col separation: " << legend->GetColumnSeparation() << std::endl;
    std::cout << "Legend row separation: " << legend->GetEntrySeparation() << std::endl;
    legend->SetColumnSeparation(0.06);
    legend->SetEntrySeparation(0.12);
    for (std::vector<std::string>::iterator g = group_ordered_list.begin();
         g != group_ordered_list.end(); ++g) {
        std::string groupName = *g;
        TH1D* hist = __err_group_map.at(groupName);
        hist->SetLineWidth(lines[groupName]);
        hist->SetLineStyle(styles[groupName]);
        hist->SetLineColor(colors[groupName]);
        hist->Draw("HISTSAME");
        legend->AddEntry(hist, keys[groupName].c_str(), "l");
    }
    
    legend->Draw("SAME");

    TPaveText* text = new TPaveText(0.16, 0.80, 0.51, 0.90, "NDC");
    text->AddText("b)  #scale[1.15]{#bar{#nu}_{#mu} + CH #rightarrow #mu^{+} + #pi^{0} + X}");
    text->AddText("#scale[0.85]{POT Normalized}");
    text->SetTextColor(kBlue);
    text->SetFillColor(kWhite);
    text->SetTextFont(42);
    text->Draw("SAME");

    std::string figure_format("pdf");
    c->Print(Form("%s.%s", figure_name.c_str(), figure_format.c_str()));

    delete legend;
    delete c;
}

void draw_systematics::draw_total(const std::string& figure_prefix) const
{
    TH1D* h1d_total_syst = add_in_quadrature(__err_hist_list, "total_syst");
    std::string figure_name = figure_prefix + "_total_syst";
    do_draw(h1d_total_syst, figure_name, "Total syst");
}

void draw_systematics::draw_all(const std::string& figure_prefix) const
{
    draw_individual(figure_prefix);
    draw_each_group(figure_prefix);
    draw_total(figure_prefix);
}

void draw_systematics::do_draw(TH1D* h, const std::string& figure_name, const std::string& legend) const 
{
    TCanvas* c = new TCanvas("c1");
    
    h->Draw("HIST");
    h->SetLineWidth(2);
    h->SetLineColor(kYellow-3);
    h->GetYaxis()->SetRangeUser(0.0001, 0.3);
    h->GetYaxis()->SetTitle("Fractional Uncertainty");
    h->GetYaxis()->CenterTitle();
    h->GetYaxis()->SetLabelSize(0.05);
    h->GetYaxis()->SetTitleSize(0.06);
    h->GetYaxis()->SetTitleOffset(1.1);
    h->GetXaxis()->SetLabelSize(0.05);
    h->GetXaxis()->SetTitleSize(0.06);
    h->GetXaxis()->SetTitleOffset(1.1);

    TLegend* l = new TLegend(0.4,0.7,0.7,0.9);
    l->SetTextSize(0.05);
    l->SetTextFont(42);
    l->SetFillColor(0);
    l->SetLineColor(0);
    l->AddEntry(h, legend.c_str(), "l");
    l->Draw("SAME");
    
    std::string figure_format("png");
    c->Print(Form("%s.%s", figure_name.c_str(), figure_format.c_str()));

    delete l;
    delete c;
}

#endif

