#ifndef new_flux_h_seen
#define new_flux_h_seen

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <cstdlib>
#include <cassert>

#include <TH1D.h>
#include <TFile.h>
#include <TRandom1.h>
#include <TString.h>

#include "../../Libraries/Folder_List.h"
#include <PlotUtils/MnvH1D.h>

using PlotUtils::MnvH1D;

/// A simple uniform histogram object
class TinyHist {
  public:
    TinyHist() {}
    explicit TinyHist(const TH1D* hist);

    int GetNbinsX() const;
    int FindBin(double val) const;
    double GetBinContent(int bin) const;
    double GetBinWidth(int bin) const;
    
    void SetBinContent(int bin, double val);
    void Reset();
    void Print();
    
  private:
    std::map<int, double> __hist_data;
    double __binwidth;
};

TinyHist::TinyHist(const TH1D* hist)
{
    __binwidth = hist->GetBinWidth(1);

    for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
        double content = hist->GetBinContent(bin);
        __hist_data.insert(std::make_pair(bin, content));
    }
}

int TinyHist::GetNbinsX() const
{
    return __hist_data.size();
}

double TinyHist::GetBinContent(int bin) const
{
    std::map<int, double>::const_iterator target = __hist_data.find(bin);
    if (target != __hist_data.end()) return target->second;

    return 0.0;
}

int TinyHist::FindBin(double val) const
{
    return int(val/__binwidth) + 1;
}

double TinyHist::GetBinWidth(int) const
{
    return __binwidth;
}

void TinyHist::Reset()
{
    for (std::map<int, double>::iterator i = __hist_data.begin();
         i != __hist_data.end(); ++i) {
        i->second = 0.0;
    }
}

void TinyHist::SetBinContent(int bin, double val)
{
    __hist_data[bin] = val;
}

void TinyHist::Print()
{
    for (std::map<int, double>::iterator i = __hist_data.begin();
         i != __hist_data.end(); ++i) {
        printf("%5d %8.4e\n", i->first, i->second);
    }
    
}

class new_flux {
  public:
    static new_flux& get();
    
    void read_oldflux_histogram(const std::string& filename);
    void read_newflux_histogram(const std::string& filename);
    void calc_weights(bool verbose = false);
    
    double get_cvweight(double enu0);
    std::vector<double> get_random_weights(double enu0);
    
  private:
    new_flux();
    new_flux(const new_flux&) {}
    new_flux& operator=(const new_flux&) { return *this;}

    
    ~new_flux();

    TinyHist __h1d_oldflux_cv;
    TinyHist __h1d_newflux_cv;
    TinyHist __h1d_flux_ratio;
    
    TinyHist __mnvh1d_newflux;
    
    std::vector<TinyHist>  __newflux_universe_histograms;
    std::vector<double>    __overflow_random_weights;
    std::map<int, std::vector<double> > __random_weight_map;
    
};

new_flux::new_flux()
{}

new_flux::~new_flux()
{}


new_flux& new_flux::get()
{
    static new_flux singleton;

    return singleton;
}


/// 1) read the MnvH1D histogram for the new flux
/// 2) get the flux central-value from the MnvH1D. This will be
///      used below to obtain the central-value reweighting function
void new_flux::read_oldflux_histogram(const std::string& filename)
{
    std::cout << "Reading old flux file: " << filename << std::endl;
    
    TFile* f = new TFile(filename.c_str(), "READ");

    MnvH1D* mnvh1d_flux = static_cast<MnvH1D*>(f->Get("flux_E_unweighted"));
    TH1D* h1d_tmp = static_cast<TH1D*>(mnvh1d_flux->GetCVHistoWithError().Clone("oldflux_cv"));

    assert(h1d_tmp);

    h1d_tmp->Rebin(5);   // new bin width 0.5 vs 0.1 GeV
    h1d_tmp->Scale(0.2); // to preserve the bin normalization
    
    __h1d_oldflux_cv = TinyHist(h1d_tmp);

    //std::cout << "Old flux histogram" << std::endl;
    //__h1d_oldflux_cv.Print();

    delete mnvh1d_flux;
    delete h1d_tmp;
    
    f->Close();
    delete f;
    
}


/// 2015/11/22: the flux file contains TH1D, not MnvH1D
/// 1) read the MnvH1D histogram for the new flux
/// 2) get the flux central-value from the MnvH1D
/// 3) get the flux systematic uncertainty from the MnvH1D
/// 4) Divide the flux central-value by the old flux to get
///      the flux ratio that is used as reweighting function
///
/// 2015/12/24: use the covariance matrix from Kevin
/// 1) Generate the random weights using the covariance matrix
/// 2) The random weights in turn are used to generate universe
///    histograms
/// 3) Use the central-value flux (also from Kevin) and the universe
///    histograms to construct an MnvH1D object
void new_flux::read_newflux_histogram(const std::string& filename)
{

    std::cout << "Reading new flux file: " << filename << std::endl;
    
    TFile* f = new TFile(filename.c_str(), "READ");
    
    MnvH1D* mnvh1d_tmp = static_cast<MnvH1D*>(f->Get("flux_E_cvweighted"));
    TH1D* h1d_tmp      = static_cast<TH1D*>(mnvh1d_tmp->GetCVHistoWithError().Clone("newflux_cv"));

    __h1d_newflux_cv = TinyHist(h1d_tmp);

    //std::cout << "New flux histogram" << std::endl;
    //__h1d_newflux_cv.Print();

    std::vector<TH1D*> universe_tmp_histograms = mnvh1d_tmp->GetVertErrorBand("Flux")->GetHists();
    for (std::vector<TH1D*>::iterator h = universe_tmp_histograms.begin();
         h != universe_tmp_histograms.end(); ++h) {
        TinyHist t(*h);
        __newflux_universe_histograms.push_back(t);
    }

    assert(!__newflux_universe_histograms.empty());
    
    __h1d_flux_ratio = __h1d_newflux_cv;
    __h1d_flux_ratio.Reset();
    
    f->Close();
    delete f;
    
}

void new_flux::calc_weights(bool verbose)
{

    std::cout << "__h1d_oldflux_cv: " << &__h1d_oldflux_cv << std::endl;
    std::cout << "__h1d_newflux_cv: " << &__h1d_newflux_cv << std::endl;
        // cv weight
        // new flux has different range up to only 20 GeV instead of 100 GeV
    for (int bin = 1; bin < __h1d_oldflux_cv.GetNbinsX(); ++bin) {
        double old_width = __h1d_oldflux_cv.GetBinWidth(bin);
        double new_width = __h1d_newflux_cv.GetBinWidth(bin);
        assert(old_width == new_width);
        
        double old_flux = __h1d_oldflux_cv.GetBinContent(bin); 
        double new_flux = __h1d_newflux_cv.GetBinContent(bin); 
        
        double ratio = (0.0 < old_flux) ? new_flux/old_flux : 1.0;
        __h1d_flux_ratio.SetBinContent(bin, ratio);
        
        if (verbose) {
            printf("Bin %02d cvweight %.4f\n", bin, ratio);
        }
        
    }

    std::cout << "New flux / Old flux ratio" << std::endl;
    //__h1d_flux_ratio.Print();
    
        // universe weights
        // random weights are specific the each bin to handle
        // bin-to-bin correlation
    for (int bin = 1; bin <= __h1d_newflux_cv.GetNbinsX(); ++bin) {
        int i = 0;
        for (std::vector<TinyHist>::const_iterator h = __newflux_universe_histograms.begin();
             h != __newflux_universe_histograms.end(); ++h) {
            double oldflux = __h1d_oldflux_cv.GetBinContent(bin);
            double universe_newflux    = h->GetBinContent(bin);

            double weight = (0.0 < oldflux) ? universe_newflux/oldflux : 1.0;
            __random_weight_map[bin].push_back(weight);

            if (verbose) {
                printf("\t Bin %2d universe %03d weight %.4f\n", bin, i, weight);
            }

            ++i;
        }
    }

        // for enu0 > 20 GeV (the range of the covariance matrix
        // from Kevin), assume 10% uncertainty
    TRandom1 generator(123456);
    for (unsigned int i = 0; i < __newflux_universe_histograms.size(); ++i) {
        double last_bin_cv = __h1d_flux_ratio.GetBinContent(__h1d_flux_ratio.GetNbinsX());
        double random_weight = generator.Gaus(last_bin_cv, 0.1);
        __overflow_random_weights.push_back(random_weight);
    }
    
}

double new_flux::get_cvweight(double enu0)
{
  
    double bin    = __h1d_flux_ratio.FindBin(enu0);
    double weight = __h1d_flux_ratio.GetBinContent(bin);
    if (bin > __h1d_newflux_cv.GetNbinsX()) {
        
        return __h1d_flux_ratio.GetBinContent(__h1d_flux_ratio.GetNbinsX());

    }
    return weight;
}


std::vector<double> new_flux::get_random_weights(double enu0)
{

    int bin = __h1d_newflux_cv.FindBin(enu0);

    if (bin > __h1d_newflux_cv.GetNbinsX()) {
        
        return __overflow_random_weights;

    }

    std::vector<double> random_weights = __random_weight_map[bin];
    assert(!random_weights.empty());
    
    return random_weights;

}

/// New flux as posted by Leo 11/20/2015:
/// 1) does not have MnvH1D object
/// 2) has different enu range 120 vs 100
/// 3) used different units for the Y axis
/// 4) only has statistical uncertainties


#endif
