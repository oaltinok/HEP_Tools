#include <cassert>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <map>
#include <cstdlib>
#include <sstream>
#include <unistd.h>
#include <sys/stat.h>
#include <getopt.h>

#include <TChain.h>
#include <TFile.h>
#include <Cintex/Cintex.h>
#include <TMath.h>
#include <TParameter.h>

#include <PlotUtils/MnvH1D.h>
#include <PlotUtils/MnvNormalization.h>

#include "hist_factory.h"
#include "make_chain.h"
#include "mc_struct.h"
#include "data_struct.h"

#include "mnvh1d_collection.h"
#include "add_error_band.h"
#include "fill_error_band.h"
#include "truth_kincalc.h"
#include "kincalc.h"
#include "random_shifts.h"
#include "constants.h"
#include "piecewise_delta_flux.h"
#include "smooth_delta_flux.h"
#include "new_flux.h"

using PlotUtils::MnvH1D;
using PlotUtils::MnvNormalizer;
using std::setw;
using std::sqrt;
using std::pow;
using std::cos;
using std::sin;


namespace {

    double calc_distance(double* vtx, double* other_vtx)
    {
        double distance = sqrt(pow(other_vtx[0]-vtx[0],2) +
                               pow(other_vtx[1]-vtx[1],2) +
                               pow(other_vtx[2]-vtx[2],2));
        
        return distance;
    }
    
    void do_count_mc(int cut, 
                     bool signal,
                     bool background,
                     std::map<int, int>& mcsig_counter,
                     std::map<int, int>& mcbkg_counter,
                     std::map<int, int>& mctot_counter)
    {
        mctot_counter[cut] += 1;
        
        if (signal)     mcsig_counter[cut] += 1;
        if (background) mcbkg_counter[cut] += 1;
    }
    
    void do_count_data(int cut,
                       std::map<int, int>& data_counter)
    {
        data_counter[cut] += 1;
    }

    const double rad2deg = TMath::RadToDeg();
    
    const double MeV = 1.0;
    const double GeV = 1.e3;
    
}

int main(int argc, char* argv[])
{

    ROOT::Cintex::Cintex::Enable();

    TH1::AddDirectory(false);

    long int __max_nevent = 1e10;
    
    bool __cvweight             = false;
    bool __minos_correction     = false;
    bool __generated            = false;
    bool __calorimetry          = false;
    bool __resonant_formula     = false;
    bool __enu_cut              = false;
    bool __w_cut                = false;
    bool __closure_test         = false;
    bool __unfold_test          = false;
    bool __reweight_flux        = false;
    bool __flux_smooth_curve    = false;
    bool __flux_piecewise_curve = false;
    bool __is_newflux           = false;

    double __emin =  1.5;
    double __emax = 10.0;
    double __wmin =  0.0;
    double __wmax =  1.8;
   
    double __nsigma   =   2.0;
    double __nsigmasb =   2.0;
    double __offset   =   0.0;
    double __sigma    =  30.0;
    double __mfp      = 150.0;
    double __reweight_min = 0.0;
    double __reweight_max = 0.0;
    double __reweight_amount = 0.0;
    
    std::string __detector;
    std::string __mc_file;
    std::string __data_file;
    std::string __histogram_file;
    std::string __label;
    std::string __newflux_file;
    std::string __oldflux_file;
    
    for (;;) {

        int option_index = 0;
        static struct option long_options[] = {
            {"cvweight",              no_argument,       0,  0 },
            {"minos-correction",      no_argument,       0,  0 },
            {"generated",             no_argument,       0,  0 },
            {"calorimetry",           no_argument,       0,  0 },
            {"resonant-formula",      no_argument,       0,  0 },
            {"with-enu-cut",          no_argument,       0,  0 },
            {"with-w-cut",            no_argument,       0,  0 },
            {"closure-test",          no_argument,       0,  0 },
            {"unfold-test",           no_argument,       0,  0 },
            {"reweight-flux",         no_argument,       0,  0 },
            {"flux-smooth-curve",     no_argument,       0,  0 },
            {"flux-piecewise-curve",  no_argument,       0,  0 },
            {"sigma",                 required_argument, 0,  0 },
            {"nsigma",                required_argument, 0,  0 },
            {"nsigmasb",              required_argument, 0,  0 },
            {"offset",                required_argument, 0,  0 },
            {"mfp",                   required_argument, 0,  0 },
            {"reweight-min",          required_argument, 0,  0 },
            {"reweight-max",          required_argument, 0,  0 },
            {"reweight-amount",       required_argument, 0,  0 },
            {"detector",              required_argument, 0,  0 },
            {"mc-file",               required_argument, 0,  0 },
            {"data-file",             required_argument, 0,  0 },
            {"histogram-file",        required_argument, 0,  0 },
            {"label",                 required_argument, 0,  0 },
            {"current-flux-file",     required_argument, 0,  0 },
            {"new-flux-file",         required_argument, 0,  0 },
            {"max-nevent",            required_argument, 0,  0 },
            {0,         0,                 0,  0 }
        };
        
        int c = getopt_long(argc, argv, "",
                            long_options, &option_index);
        
        if (c != 0 ) break;
        
        std::string opt_name(long_options[option_index].name);

        
        if      (opt_name == "cvweight")             __cvweight         = true;
        else if (opt_name == "minos-correction")     __minos_correction = true;
        else if (opt_name == "generated")            __generated        = true;
        else if (opt_name == "calorimetry")          __calorimetry      = true;
        else if (opt_name == "resonant-formula")     __resonant_formula = true;
        else if (opt_name == "with-enu-cut")         __enu_cut          = true;
        else if (opt_name == "with-w-cut")           __w_cut            = true;
        else if (opt_name == "closure-test")         __closure_test     = true;
        else if (opt_name == "unfold-test")          __unfold_test      = true;
        else if (opt_name == "reweight-flux")        __reweight_flux    = true;
        else if (opt_name == "flux-smooth-curve")    __flux_smooth_curve    = true;
        else if (opt_name == "flux-piecewise-curve") __flux_piecewise_curve = true;
        else if (opt_name == "sigma")                __sigma    = atof(optarg);
        else if (opt_name == "nsigma")               __nsigma   = atof(optarg);
        else if (opt_name == "nsigmasb")             __nsigmasb = atof(optarg);
        else if (opt_name == "offset")               __offset   = atof(optarg);
        else if (opt_name == "mfp")                  __mfp      = atof(optarg);
        else if (opt_name == "detector")             __detector          = optarg;
        else if (opt_name == "mc-file")              __mc_file           = optarg;
        else if (opt_name == "data-file")            __data_file         = optarg;
        else if (opt_name == "histogram-file")       __histogram_file    = optarg;
        else if (opt_name == "label")                __label             = optarg;
        else if (opt_name == "new-flux-file")        __newflux_file      = optarg;
        else if (opt_name == "current-flux-file")    __oldflux_file      = optarg;
        else if (opt_name == "reweight-min")         __reweight_min      = atof(optarg);
        else if (opt_name == "reweight-max")         __reweight_max      = atof(optarg);
        else if (opt_name == "reweight-amount")      __reweight_amount   = atof(optarg);
        else if (opt_name == "max-nevent")           __max_nevent        = atoi(optarg);
        else {}
    }

    __is_newflux = !__newflux_file.empty();
   
    std::cout << "__cvweight              " << __cvweight         << std::endl;
    std::cout << "__minos_correction      " << __minos_correction << std::endl;
    std::cout << "__generated             " << __generated        << std::endl;
    std::cout << "__calorimetry           " << __calorimetry      << std::endl;
    std::cout << "__resonant_formula      " << __resonant_formula << std::endl;
    std::cout << "__enu_cut               " << __enu_cut          << std::endl;
    std::cout << "__w_cut                 " << __w_cut            << std::endl;
    std::cout << "__closure_test          " << __closure_test     << std::endl;
    std::cout << "__unfold_test           " << __unfold_test      << std::endl;
    std::cout << "__reweight_flux         " << __reweight_flux    << std::endl;
    std::cout << "__flux_smooth_curve     " << __flux_smooth_curve    << std::endl;
    std::cout << "__flux_piecewise_curve  " << __flux_piecewise_curve << std::endl;
    std::cout << "__sigma                 " << __sigma            << std::endl;
    std::cout << "__nsigma                " << __nsigma           << std::endl;
    std::cout << "__nsigmasb              " << __nsigmasb         << std::endl;
    std::cout << "__offset                " << __offset           << std::endl;
    std::cout << "__mfp                   " << __mfp              << std::endl;
    std::cout << "__reweight_min          " << __reweight_min     << std::endl;
    std::cout << "__reweight_max          " << __reweight_max     << std::endl;
    std::cout << "__reweight_amount       " << __reweight_amount  << std::endl;
    std::cout << "__detector              " << __detector         << std::endl;
    std::cout << "__mc_file               " << __mc_file          << std::endl;
    std::cout << "__data_file             " << __data_file        << std::endl;
    std::cout << "__histogram_file        " << __histogram_file   << std::endl;
    std::cout << "__label                 " << __label            << std::endl;
    std::cout << "__newflux_file          " << __newflux_file     << std::endl;
    std::cout << "__oldflux_file          " << __oldflux_file     << std::endl;
    std::cout << "__max_nevent            " << __max_nevent       << std::endl;
    
    
    // Make sure detector names are always correct
    if (__detector != "minerva" && __detector != "downstream") assert(false);

    // Avoid the situation where the detector name is 'downstream', but the event files
    // are 'minerva' by mistake
    if (__detector == "downstream" &&
        __mc_file.find("downstream") == std::string::npos) assert(false);

        // Make sure the files exists
    struct stat buffer;
    if (stat(__mc_file.c_str(), &buffer) != 0) {
        std::cerr << "File not found: " << __mc_file << std::endl;
        exit(1);
    }

    if (stat(__data_file.c_str(), &buffer) != 0) {
        std::cerr << "File not found: " << __data_file << std::endl;
        exit(1);
    }
    
    
    
    // Make sure that only one method is chose
    if (( __calorimetry &&  __resonant_formula) ||
        (!__calorimetry && !__resonant_formula)) assert(false);
  
    
    const double m0 = constants().m0;

    const double lower_peak = m0 - __nsigma * __sigma;
    const double upper_peak = m0 + __nsigma * __sigma;
    const double sblow_lower = std::max(0.0, m0 - (__nsigma + __offset + __nsigmasb) * __sigma);
    const double sblow_upper  = m0 - (__nsigma + __offset) * __sigma;
    const double sbhigh_lower = m0 + (__nsigma + __offset) * __sigma;
    const double sbhigh_upper = m0 + (__nsigma + __offset + __nsigmasb) * __sigma;

    std::cout << "lower_peak:   " << lower_peak << std::endl;
    std::cout << "upper_peak:   " << upper_peak << std::endl;
    std::cout << "sblow_lower:  " << sblow_lower << std::endl;
    std::cout << "sblow_upper:  " << sblow_upper << std::endl;
    std::cout << "sbhigh_lower: " << sbhigh_lower << std::endl;
    std::cout << "sbhigh_upper: " << sbhigh_upper << std::endl;


    MnvNormalizer normalizer("minerva5");

    std::cout << "Creating histograms" << std::endl;

    MnvH1D* mnvh1d_mc_mgg    = new MnvH1D(*(make_hist_wrapper_1d(0, Form("%s-mc_mgg",    __detector.c_str()))));
    MnvH1D* mnvh1d_mcsig_mgg = new MnvH1D(*(make_hist_wrapper_1d(0, Form("%s-mcsig_mgg", __detector.c_str()))));
    MnvH1D* mnvh1d_mcbkg_mgg = new MnvH1D(*(make_hist_wrapper_1d(0, Form("%s-mcbkg_mgg", __detector.c_str()))));

        // before energy scale correction
    MnvH1D* mnvh1d_mc_mgg2    = new MnvH1D(*(make_hist_wrapper_1d(0, Form("%s-mc_mgg2",    __detector.c_str()))));
    MnvH1D* mnvh1d_mcsig_mgg2 = new MnvH1D(*(make_hist_wrapper_1d(0, Form("%s-mcsig_mgg2", __detector.c_str()))));
    MnvH1D* mnvh1d_mcbkg_mgg2 = new MnvH1D(*(make_hist_wrapper_1d(0, Form("%s-mcbkg_mgg2", __detector.c_str()))));

    
        // costheta
    MnvH1D* mnvh1d_mc_costheta    = new MnvH1D(*make_hist_wrapper_1d(1, Form("%s-mc_costheta_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_mcsig_costheta = new MnvH1D(*make_hist_wrapper_1d(1, Form("%s-mcsig_costheta_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_mcbkg_costheta = new MnvH1D(*make_hist_wrapper_1d(1, Form("%s-mcbkg_costheta_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gn_costheta    = new MnvH1D(*make_hist_wrapper_1d(1, Form("%s-gn_costheta_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_gnsig_costheta = new MnvH1D(*make_hist_wrapper_1d(1, Form("%s-gnsig_costheta_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gnbkg_costheta = new MnvH1D(*make_hist_wrapper_1d(1, Form("%s-gnbkg_costheta_masspeak", __detector.c_str())));

        // theta
    MnvH1D* mnvh1d_mc_theta    = new MnvH1D(*make_hist_wrapper_1d(5, Form("%s-mc_theta_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_mcsig_theta = new MnvH1D(*make_hist_wrapper_1d(5, Form("%s-mcsig_theta_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_mcbkg_theta = new MnvH1D(*make_hist_wrapper_1d(5, Form("%s-mcbkg_theta_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gn_theta    = new MnvH1D(*make_hist_wrapper_1d(5, Form("%s-gn_theta_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_gnsig_theta = new MnvH1D(*make_hist_wrapper_1d(5, Form("%s-gnsig_theta_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gnbkg_theta = new MnvH1D(*make_hist_wrapper_1d(5, Form("%s-gnbkg_theta_masspeak", __detector.c_str())));
    
        // pi0 total energy
    MnvH1D* mnvh1d_mc_pienergy    = new MnvH1D(*make_hist_wrapper_1d(10, Form("%s-mc_pienergy_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_mcsig_pienergy = new MnvH1D(*make_hist_wrapper_1d(10, Form("%s-mcsig_pienergy_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_mcbkg_pienergy = new MnvH1D(*make_hist_wrapper_1d(10, Form("%s-mcbkg_pienergy_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gn_pienergy    = new MnvH1D(*make_hist_wrapper_1d(10, Form("%s-gn_pienergy_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_gnsig_pienergy = new MnvH1D(*make_hist_wrapper_1d(10, Form("%s-gnsig_pienergy_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gnbkg_pienergy = new MnvH1D(*make_hist_wrapper_1d(10, Form("%s-gnbkg_pienergy_masspeak", __detector.c_str())));

        // pi0 momentum
    MnvH1D* mnvh1d_mc_pimom    = new MnvH1D(*make_hist_wrapper_1d(15, Form("%s-mc_pimom_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_mcsig_pimom = new MnvH1D(*make_hist_wrapper_1d(15, Form("%s-mcsig_pimom_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_mcbkg_pimom = new MnvH1D(*make_hist_wrapper_1d(15, Form("%s-mcbkg_pimom_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gn_pimom    = new MnvH1D(*make_hist_wrapper_1d(15, Form("%s-gn_pimom_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_gnsig_pimom = new MnvH1D(*make_hist_wrapper_1d(15, Form("%s-gnsig_pimom_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gnbkg_pimom = new MnvH1D(*make_hist_wrapper_1d(15, Form("%s-gnbkg_pimom_masspeak", __detector.c_str())));

        // pi0 kinetic energy
    MnvH1D* mnvh1d_mc_kinetic    = new MnvH1D(*make_hist_wrapper_1d(20, Form("%s-mc_kinetic_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_mcsig_kinetic = new MnvH1D(*make_hist_wrapper_1d(20, Form("%s-mcsig_kinetic_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_mcbkg_kinetic = new MnvH1D(*make_hist_wrapper_1d(20, Form("%s-mcbkg_kinetic_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gn_kinetic    = new MnvH1D(*make_hist_wrapper_1d(20, Form("%s-gn_kinetic_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_gnsig_kinetic = new MnvH1D(*make_hist_wrapper_1d(20, Form("%s-gnsig_kinetic_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gnbkg_kinetic = new MnvH1D(*make_hist_wrapper_1d(20, Form("%s-gnbkg_kinetic_masspeak", __detector.c_str())));

        // muon momentum
    MnvH1D* mnvh1d_mc_mmom    = new MnvH1D(*make_hist_wrapper_1d(25, Form("%s-mc_mmom_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_mcsig_mmom = new MnvH1D(*make_hist_wrapper_1d(25, Form("%s-mcsig_mmom_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_mcbkg_mmom = new MnvH1D(*make_hist_wrapper_1d(25, Form("%s-mcbkg_mmom_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gn_mmom    = new MnvH1D(*make_hist_wrapper_1d(25, Form("%s-gn_mmom_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_gnsig_mmom = new MnvH1D(*make_hist_wrapper_1d(25, Form("%s-gnsig_mmom_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gnbkg_mmom = new MnvH1D(*make_hist_wrapper_1d(25, Form("%s-gnbkg_mmom_masspeak", __detector.c_str())));
    
        // muon longitudinal momentum
    MnvH1D* mnvh1d_mc_pz    = new MnvH1D(*make_hist_wrapper_1d(50, Form("%s-mc_pz_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_mcsig_pz = new MnvH1D(*make_hist_wrapper_1d(50, Form("%s-mcsig_pz_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_mcbkg_pz = new MnvH1D(*make_hist_wrapper_1d(50, Form("%s-mcbkg_pz_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gn_pz    = new MnvH1D(*make_hist_wrapper_1d(50, Form("%s-gn_pz_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_gnsig_pz = new MnvH1D(*make_hist_wrapper_1d(50, Form("%s-gnsig_pz_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gnbkg_pz = new MnvH1D(*make_hist_wrapper_1d(50, Form("%s-gnbkg_pz_masspeak", __detector.c_str())));
    
    
        // muon transverse momentum
    MnvH1D* mnvh1d_mc_pT    = new MnvH1D(*make_hist_wrapper_1d(50, Form("%s-mc_pT_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_mcsig_pT = new MnvH1D(*make_hist_wrapper_1d(50, Form("%s-mcsig_pT_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_mcbkg_pT = new MnvH1D(*make_hist_wrapper_1d(50, Form("%s-mcbkg_pT_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gn_pT    = new MnvH1D(*make_hist_wrapper_1d(50, Form("%s-gn_pT_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_gnsig_pT = new MnvH1D(*make_hist_wrapper_1d(50, Form("%s-gnsig_pT_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gnbkg_pT = new MnvH1D(*make_hist_wrapper_1d(50, Form("%s-gnbkg_pT_masspeak", __detector.c_str())));

        // muon theta
    MnvH1D* mnvh1d_mc_mutheta    = new MnvH1D(*make_hist_wrapper_1d(30, Form("%s-mc_mutheta_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_mcsig_mutheta = new MnvH1D(*make_hist_wrapper_1d(30, Form("%s-mcsig_mutheta_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_mcbkg_mutheta = new MnvH1D(*make_hist_wrapper_1d(30, Form("%s-mcbkg_mutheta_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gn_mutheta    = new MnvH1D(*make_hist_wrapper_1d(30, Form("%s-gn_mutheta_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_gnsig_mutheta = new MnvH1D(*make_hist_wrapper_1d(30, Form("%s-gnsig_mutheta_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gnbkg_mutheta = new MnvH1D(*make_hist_wrapper_1d(30, Form("%s-gnbkg_mutheta_masspeak", __detector.c_str())));

        // Neutrino energy
    MnvH1D* mnvh1d_mc_enu    = new MnvH1D(*make_hist_wrapper_1d(35, Form("%s-mc_enu_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_mcsig_enu = new MnvH1D(*make_hist_wrapper_1d(35, Form("%s-mcsig_enu_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_mcbkg_enu = new MnvH1D(*make_hist_wrapper_1d(35, Form("%s-mcbkg_enu_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gn_enu    = new MnvH1D(*make_hist_wrapper_1d(35, Form("%s-gn_enu_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_gnsig_enu = new MnvH1D(*make_hist_wrapper_1d(35, Form("%s-gnsig_enu_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gnbkg_enu = new MnvH1D(*make_hist_wrapper_1d(35, Form("%s-gnbkg_enu_masspeak", __detector.c_str())));
    
        // Q^2
    MnvH1D* mnvh1d_mc_q2    = new MnvH1D(*make_hist_wrapper_1d(40, Form("%s-mc_q2_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_mcsig_q2 = new MnvH1D(*make_hist_wrapper_1d(40, Form("%s-mcsig_q2_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_mcbkg_q2 = new MnvH1D(*make_hist_wrapper_1d(40, Form("%s-mcbkg_q2_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gn_q2    = new MnvH1D(*make_hist_wrapper_1d(40, Form("%s-gn_q2_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_gnsig_q2 = new MnvH1D(*make_hist_wrapper_1d(40, Form("%s-gnsig_q2_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gnbkg_q2 = new MnvH1D(*make_hist_wrapper_1d(40, Form("%s-gnbkg_q2_masspeak", __detector.c_str())));
    
        // Hadronic invariant mass
    MnvH1D* mnvh1d_mc_w    = new MnvH1D(*make_hist_wrapper_1d(45, Form("%s-mc_w_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_mcsig_w = new MnvH1D(*make_hist_wrapper_1d(45, Form("%s-mcsig_w_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_mcbkg_w = new MnvH1D(*make_hist_wrapper_1d(45, Form("%s-mcbkg_w_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gn_w    = new MnvH1D(*make_hist_wrapper_1d(45, Form("%s-gn_w_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_gnsig_w = new MnvH1D(*make_hist_wrapper_1d(45, Form("%s-gnsig_w_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gnbkg_w = new MnvH1D(*make_hist_wrapper_1d(45, Form("%s-gnbkg_w_masspeak", __detector.c_str())));
    
    

    std::cout << "Add vertical error band" << std::endl;

        // invariant mass
    add_error_bands<0>(mnvh1d_mc_mgg,     __is_newflux);
    add_error_bands<0>(mnvh1d_mcsig_mgg,  __is_newflux);
    add_error_bands<0>(mnvh1d_mcbkg_mgg,  __is_newflux);

    add_error_bands<0>(mnvh1d_mc_mgg2,     __is_newflux);
    add_error_bands<0>(mnvh1d_mcsig_mgg2,  __is_newflux);
    add_error_bands<0>(mnvh1d_mcbkg_mgg2,  __is_newflux);

        // costheta
    add_error_bands<1>(mnvh1d_mc_costheta,     __is_newflux);
    add_error_bands<1>(mnvh1d_mcsig_costheta,  __is_newflux);
    add_error_bands<1>(mnvh1d_mcbkg_costheta,  __is_newflux);
    add_error_bands<1>(mnvh1d_gn_costheta,     __is_newflux);
    add_error_bands<1>(mnvh1d_gnsig_costheta,  __is_newflux);
    add_error_bands<1>(mnvh1d_gnbkg_costheta,  __is_newflux);
    
        // theta
    add_error_bands<5>(mnvh1d_mc_theta,     __is_newflux);
    add_error_bands<5>(mnvh1d_mcsig_theta,  __is_newflux);
    add_error_bands<5>(mnvh1d_mcbkg_theta,  __is_newflux);
    add_error_bands<5>(mnvh1d_gn_theta,     __is_newflux);
    add_error_bands<5>(mnvh1d_gnsig_theta,  __is_newflux);
    add_error_bands<5>(mnvh1d_gnbkg_theta,  __is_newflux);

        // pienergy
    add_error_bands<10>(mnvh1d_mc_pienergy,     __is_newflux);
    add_error_bands<10>(mnvh1d_mcsig_pienergy,  __is_newflux);
    add_error_bands<10>(mnvh1d_mcbkg_pienergy,  __is_newflux);
    add_error_bands<10>(mnvh1d_gn_pienergy,     __is_newflux);
    add_error_bands<10>(mnvh1d_gnsig_pienergy,  __is_newflux);
    add_error_bands<10>(mnvh1d_gnbkg_pienergy,  __is_newflux);

        // pimom
    add_error_bands<15>(mnvh1d_mc_pimom,     __is_newflux);
    add_error_bands<15>(mnvh1d_mcsig_pimom,  __is_newflux);
    add_error_bands<15>(mnvh1d_mcbkg_pimom,  __is_newflux);
    add_error_bands<15>(mnvh1d_gn_pimom,     __is_newflux);
    add_error_bands<15>(mnvh1d_gnsig_pimom,  __is_newflux);
    add_error_bands<15>(mnvh1d_gnbkg_pimom,  __is_newflux);

        // kinetic
    add_error_bands<20>(mnvh1d_mc_kinetic,     __is_newflux);
    add_error_bands<20>(mnvh1d_mcsig_kinetic,  __is_newflux);
    add_error_bands<20>(mnvh1d_mcbkg_kinetic,  __is_newflux);
    add_error_bands<20>(mnvh1d_gn_kinetic,     __is_newflux);
    add_error_bands<20>(mnvh1d_gnsig_kinetic,  __is_newflux);
    add_error_bands<20>(mnvh1d_gnbkg_kinetic,  __is_newflux);

        // muon momentum
    add_error_bands<25>(mnvh1d_mc_mmom,     __is_newflux);
    add_error_bands<25>(mnvh1d_mcsig_mmom,  __is_newflux);
    add_error_bands<25>(mnvh1d_mcbkg_mmom,  __is_newflux);
    add_error_bands<25>(mnvh1d_gn_mmom,     __is_newflux);
    add_error_bands<25>(mnvh1d_gnsig_mmom,  __is_newflux);
    add_error_bands<25>(mnvh1d_gnbkg_mmom,  __is_newflux);


        // muon longitudinal momentum
    add_error_bands<50>(mnvh1d_mc_pz,     __is_newflux);
    add_error_bands<50>(mnvh1d_mcsig_pz,  __is_newflux);
    add_error_bands<50>(mnvh1d_mcbkg_pz,  __is_newflux);
    add_error_bands<50>(mnvh1d_gn_pz,     __is_newflux);
    add_error_bands<50>(mnvh1d_gnsig_pz,  __is_newflux);
    add_error_bands<50>(mnvh1d_gnbkg_pz,  __is_newflux);

        // muon transverse momentum
    add_error_bands<55>(mnvh1d_mc_pT,     __is_newflux);
    add_error_bands<55>(mnvh1d_mcsig_pT,  __is_newflux);
    add_error_bands<55>(mnvh1d_mcbkg_pT,  __is_newflux);
    add_error_bands<55>(mnvh1d_gn_pT,     __is_newflux);
    add_error_bands<55>(mnvh1d_gnsig_pT,  __is_newflux);
    add_error_bands<55>(mnvh1d_gnbkg_pT,  __is_newflux);

    
        // muon theta
    add_error_bands<30>(mnvh1d_mc_mutheta,     __is_newflux);
    add_error_bands<30>(mnvh1d_mcsig_mutheta,  __is_newflux);
    add_error_bands<30>(mnvh1d_mcbkg_mutheta,  __is_newflux);
    add_error_bands<30>(mnvh1d_gn_mutheta,     __is_newflux);
    add_error_bands<30>(mnvh1d_gnsig_mutheta,  __is_newflux);
    add_error_bands<30>(mnvh1d_gnbkg_mutheta,  __is_newflux);


        // Neutrino energy
    add_error_bands<35>(mnvh1d_mc_enu,     __is_newflux);
    add_error_bands<35>(mnvh1d_mcsig_enu,  __is_newflux);
    add_error_bands<35>(mnvh1d_mcbkg_enu,  __is_newflux);
    add_error_bands<35>(mnvh1d_gn_enu,     __is_newflux);
    add_error_bands<35>(mnvh1d_gnsig_enu,  __is_newflux);
    add_error_bands<35>(mnvh1d_gnbkg_enu,  __is_newflux);


        // Q2
    add_error_bands<40>(mnvh1d_mc_q2,     __is_newflux);
    add_error_bands<40>(mnvh1d_mcsig_q2,  __is_newflux);
    add_error_bands<40>(mnvh1d_mcbkg_q2,  __is_newflux);
    add_error_bands<40>(mnvh1d_gn_q2,     __is_newflux);
    add_error_bands<40>(mnvh1d_gnsig_q2,  __is_newflux);
    add_error_bands<40>(mnvh1d_gnbkg_q2,  __is_newflux);


        // Hadronic invariant mass
    add_error_bands<45>(mnvh1d_mc_w,     __is_newflux);
    add_error_bands<45>(mnvh1d_mcsig_w,  __is_newflux);
    add_error_bands<45>(mnvh1d_mcbkg_w,  __is_newflux);
    add_error_bands<45>(mnvh1d_gn_w,     __is_newflux);
    add_error_bands<45>(mnvh1d_gnsig_w,  __is_newflux);
    add_error_bands<45>(mnvh1d_gnbkg_w,  __is_newflux);

    if (__is_newflux) {
        new_flux::get().read_oldflux_histogram(__oldflux_file); 
        new_flux::get().read_newflux_histogram(__newflux_file);
        new_flux::get().calc_weights(true);
    }
    
    std::cout << "Creating MC chain... " << std::endl;
    TChain* mc_chain = make_chain(__mc_file.c_str());
    
    std::cout << "Total MC entries: " << mc_chain->GetEntries() << std::endl;

    mc_struct* mc_object = new mc_struct;
    mc_object->Init(mc_chain);

    std::cout.setf(std::ios_base::fixed);
    std::cout.precision(2);

    
    std::map<int, int> mcsig_counter;
    std::map<int, int> mcbkg_counter;
    std::map<int, int> mctot_counter;
    std::map<int, int> data_counter;
    
    bool is_downstream = (__detector == "downstream");
    const int n_mc = mc_chain->GetEntriesFast();
    int readEntry = 0;
    int selectedEvent = 0;
    
    std::cout << "Fill histograms..." << std::endl;
    for (int entry = 0; entry < n_mc; ++entry) { 
        
        mc_object->GetEntry(entry);

        ++readEntry;
        
        if (readEntry > __max_nevent) break;
        
        if (mc_object->mgg > 1e4) continue;

        if (readEntry%1000 == 0) std::cout << "Processed " << (int)(100.0*readEntry/n_mc) << " (%) "<< std::endl;
                      
        const double costheta0  = cos(mc_object->pithetab0 * deg2rad);
        double theta0    = mc_object->pithetab0;
        double pimom0    = mc_object->pimom0[3]/GeV;
        double pienergy0 = mc_object->pienergy0/GeV;
        double kinetic0  = (mc_object->pienergy0-m0)/GeV;
        double mmom0     = mc_object->truth_fslepton_P; // already in GeV
        double pz0       = mmom0 * cos(mc_object->truth_fslepton_theta);
        double pT0       = mmom0 * sin(mc_object->truth_fslepton_theta);
        double mutheta0   = mc_object->truth_fslepton_theta * rad2deg;
        double enu0      = mc_object->mc_incomingE/GeV;
        double q20       = mc_object->mc_Q2/(GeV * GeV);
        double w0        = mc_object->mc_w/GeV;
        
        if (!__generated) {
            
            cc1pi0::truth_kincalc calc(mc_object->mc_incomingE,
                                       mc_object->truth_fslepton_P * GeV,
                                       mc_object->truth_fslepton_theta);
                                       
            q20        = calc.q2()/(GeV * GeV);
            w0         = calc.w()/GeV;
        }
           
        int mc_current = mc_object->mc_current;
        int mc_incoming = mc_object->mc_incoming;
        bool truth_is_cc1pi0 = mc_object->truth_is_cc1pi0;
        bool truth_is_fiducial = mc_object->truth_is_fiducial;
        bool enu_cut = (enu0 > __emin) && (enu0 < __emax);
        bool w_cut   = (w0   > __wmin) && (w0   < __wmax);
        
        bool signal = truth_is_fiducial && truth_is_cc1pi0 && mc_incoming == -14 && mc_current==1;
        if (__enu_cut) signal = signal && enu_cut;
        if (__w_cut)   signal = signal && w_cut;
        
        if (is_downstream) signal = signal && (mc_object->mc_vtx[2] > 7200.0);
        bool background = !signal;
                
        double cvweight = 1.0;
        if (__is_newflux) {
            
            cvweight = new_flux::get().get_cvweight(enu0);

        } else {
        
            if (__cvweight) cvweight *= mc_object->wgt;

        }
        
        if (__minos_correction) {
            const double minos_eff_correction = normalizer.GetCorrection(mc_object->minos_trk_p);
            cvweight *= minos_eff_correction;
        }

        if (__reweight_flux) {
            
            
            if (__flux_smooth_curve) {
                
                double amount = smooth_delta_flux::get().calculate(enu0);
                cvweight *= (1.0 + amount);
                
            } else if (__flux_piecewise_curve) {
                
                double amount = piecewise_delta_flux::get().calculate(enu0);
                cvweight *= (1.0 + amount);
                
            } else {

                if (__reweight_min < enu0 && enu0 < __reweight_max) cvweight *= ( 1.0 + __reweight_amount);
                
            }

        }
        
 
        do_count_mc(0, signal, background, mcsig_counter, mcbkg_counter, mctot_counter);
        if (mc_object->survive_minos_match != 1) continue;        
        
        do_count_mc(1, signal, background, mcsig_counter, mcbkg_counter, mctot_counter);
        if (mc_object->minos_trk_qp <= 0) continue;        

        do_count_mc(2, signal, background, mcsig_counter, mcbkg_counter, mctot_counter);
        if (mc_object->survive_fiducial != 1) continue;
        if (is_downstream && mc_object->vtx[2] < 7200.0) continue;

        
        bool oneTrack1 = mc_object->primary_multiplicity2 == 0;
        bool oneTrack2 = mc_object->max_ntrack < 2;
        bool gammaLike = mc_object->max_deviation < 50;

        do_count_mc(3, signal, background, mcsig_counter, mcbkg_counter, mctot_counter);
        if (!oneTrack1) continue;

        do_count_mc(4, signal, background, mcsig_counter, mcbkg_counter, mctot_counter);
        if (!oneTrack2) continue;

        do_count_mc(5, signal, background, mcsig_counter, mcbkg_counter, mctot_counter);
        if (!gammaLike) continue;

        do_count_mc(6, signal, background, mcsig_counter, mcbkg_counter, mctot_counter);
        if (mc_object->ntgtevis >= 20.0) continue;        

        do_count_mc(7, signal, background, mcsig_counter, mcbkg_counter, mctot_counter);
        if (mc_object->otherevis <= 80.0) continue;        

        do_count_mc(8, signal, background, mcsig_counter, mcbkg_counter, mctot_counter);
        if (mc_object->otherevis >= 2000.0) continue;        

        do_count_mc(9, signal, background, mcsig_counter, mcbkg_counter, mctot_counter);
        bool goodGamma = mc_object->is_GoodBlob1 && mc_object->is_GoodBlob2;
        if (!goodGamma) continue;        

        do_count_mc(10, signal, background, mcsig_counter, mcbkg_counter, mctot_counter);
        if (mc_object->final_blob_ncx[1]<2) continue;

        do_count_mc(11, signal, background, mcsig_counter, mcbkg_counter, mctot_counter);
        double g1mfp = calc_distance(mc_object->vtx, mc_object->RE_photon_vertex_1);
        double g2mfp = calc_distance(mc_object->vtx, mc_object->RE_photon_vertex_2);
        
        if (g1mfp < __mfp) continue;
        if (g2mfp < __mfp) continue;    

        do_count_mc(12, signal, background, mcsig_counter, mcbkg_counter, mctot_counter);

        TLorentzVector muon_p4(mc_object->muon_px,
                               mc_object->muon_py,
                               mc_object->muon_pz,
                               mc_object->muon_E);

        TLorentzVector muon_p4_rot = muon_p4;
        muon_p4_rot.RotateX(-3.3 * deg2rad);
      
        
        TLorentzVector pion_p4(mc_object->pimom[0],
                               mc_object->pimom[1],
                               mc_object->pimom[2],
                               mc_object->pienergy); 
        
        const double alpha = m0/115.0;// from fit MC to gaussian
        
        const double costheta  = cos(mc_object->pithetab * deg2rad);
        double theta     = mc_object->pithetab;
        double pimom     = alpha*mc_object->pimom[3]/GeV;
        double pienergy  = alpha*mc_object->pienergy/GeV;
        double kinetic   = (alpha*mc_object->pienergy-m0)/GeV;
        double mmom      = mc_object->muon_p/GeV;
        double pz        = muon_p4_rot.Pz()/GeV;
        double pT        = muon_p4_rot.Pt()/GeV;
        double mutheta    = mc_object->muon_theta * rad2deg;
        const double mgg = alpha*mc_object->mgg; 

        
        double enu       = 0.0;
        double q2        = 0.0;
        double w         = 0.0;
       
        if (__calorimetry) {
            
            cc1pi0::kincalc calc(muon_p4, pion_p4, alpha,
                                    true,
                                    mc_object->Vertex_blob_energy,
                                    mc_object->Dispersed_blob_energy); 
        
            enu       = calc.enu()/GeV; 
            q2        = calc.q2()/(GeV * GeV);
            w         = calc.w()/GeV;
       
        }
        
        if (__resonant_formula) {
            
            cc1pi0::kincalc calc(muon_p4, pion_p4, alpha);
        
            enu       = calc.enu()/GeV; 
            q2        = calc.q2()/(GeV * GeV);
            w         = calc.w()/GeV;
       
        }

        bool not_within_enu_range_cv = (enu < __emin) || (__emax < enu);
        bool not_within_w_range_cv   = (w < __wmin)   || (__wmax < w);
        
        bool pass_enu_cut_cv = true;
        if (__enu_cut && not_within_enu_range_cv) pass_enu_cut_cv = false;

        bool pass_w_cut_cv = true;
        if (__w_cut && not_within_w_range_cv) pass_w_cut_cv = false;
        
        bool pass_mgg_cut_cv = (lower_peak < mgg) && (mgg < upper_peak);

        if (pass_enu_cut_cv)  do_count_mc(13, signal, background, mcsig_counter, mcbkg_counter, mctot_counter);
        if (pass_w_cut_cv)    do_count_mc(14, signal, background, mcsig_counter, mcbkg_counter, mctot_counter);


            // handle mgg distribution differently :(, i.e., without the invariant mass cut
        if (pass_enu_cut_cv && pass_w_cut_cv) {

            mnvh1d_mc_mgg->Fill(mgg, cvweight);
            mnvh1d_mc_mgg2->Fill(mc_object->mgg, cvweight);

            fill_vert_error_bands(mnvh1d_mc_mgg,  mgg,            mc_object, cvweight,  __is_newflux);
            fill_vert_error_bands(mnvh1d_mc_mgg2, mc_object->mgg, mc_object, cvweight,  __is_newflux);

            if (signal) {

                mnvh1d_mcsig_mgg->Fill(mgg, cvweight);
                mnvh1d_mcsig_mgg2->Fill(mc_object->mgg, cvweight);

                fill_vert_error_bands(mnvh1d_mcsig_mgg,  mgg,              mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_mcsig_mgg2, mc_object->mgg,   mc_object, cvweight,  __is_newflux);

            }

            if (background) {

                mnvh1d_mcbkg_mgg->Fill(mgg, cvweight);
                mnvh1d_mcbkg_mgg2->Fill(mc_object->mgg, cvweight);
                fill_vert_error_bands(mnvh1d_mcbkg_mgg,  mgg,            mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_mcbkg_mgg2, mc_object->mgg, mc_object, cvweight,  __is_newflux);
                
            }

        }
        
            // handle cv: total mc, signal, background
        if (pass_enu_cut_cv && pass_w_cut_cv && pass_mgg_cut_cv) {
            ++selectedEvent;            

            do_count_mc(15, signal, background, mcsig_counter, mcbkg_counter, mctot_counter);
        
            mnvh1d_mc_costheta->Fill(costheta, cvweight);
            mnvh1d_mc_theta->Fill(theta,       cvweight);
            mnvh1d_mc_pienergy->Fill(pienergy, cvweight);
            mnvh1d_mc_pimom->Fill(pimom,       cvweight);
            mnvh1d_mc_kinetic->Fill(kinetic,   cvweight);
            mnvh1d_mc_mmom->Fill(mmom,         cvweight);
            mnvh1d_mc_pz->Fill(pz,             cvweight);
            mnvh1d_mc_pT->Fill(pT,             cvweight);
            mnvh1d_mc_mutheta->Fill(mutheta,   cvweight);
            mnvh1d_mc_enu->Fill(enu,           cvweight);
            mnvh1d_mc_q2->Fill(q2,             cvweight);
            mnvh1d_mc_w->Fill(w,               cvweight);
        
            mnvh1d_gn_costheta->Fill(costheta0, cvweight);
            mnvh1d_gn_theta->Fill(theta0,       cvweight);
            mnvh1d_gn_pienergy->Fill(pienergy0, cvweight);
            mnvh1d_gn_pimom->Fill(pimom0,       cvweight);
            mnvh1d_gn_kinetic->Fill(kinetic0,   cvweight);
            mnvh1d_gn_mmom->Fill(mmom0,         cvweight);
            mnvh1d_gn_pz->Fill(pz0,             cvweight);
            mnvh1d_gn_pT->Fill(pT0,             cvweight);
            mnvh1d_gn_mutheta->Fill(mutheta0,   cvweight);
            mnvh1d_gn_enu->Fill(enu0,           cvweight);
            mnvh1d_gn_q2->Fill(q20,             cvweight);
            mnvh1d_gn_w->Fill(w0,               cvweight);

                        
            fill_vert_error_bands(mnvh1d_mc_costheta, costheta, mc_object, cvweight,  __is_newflux);
            fill_vert_error_bands(mnvh1d_mc_theta,    theta,    mc_object, cvweight,  __is_newflux);
            fill_vert_error_bands(mnvh1d_mc_pienergy, pienergy, mc_object, cvweight,  __is_newflux);
            fill_vert_error_bands(mnvh1d_mc_pimom,    pimom,    mc_object, cvweight,  __is_newflux);
            fill_vert_error_bands(mnvh1d_mc_kinetic,  kinetic,  mc_object, cvweight,  __is_newflux);
            fill_vert_error_bands(mnvh1d_mc_mmom,     mmom,     mc_object, cvweight,  __is_newflux);
            fill_vert_error_bands(mnvh1d_mc_pz,       pz,       mc_object, cvweight,  __is_newflux);
            fill_vert_error_bands(mnvh1d_mc_pT,       pT,       mc_object, cvweight,  __is_newflux);
            fill_vert_error_bands(mnvh1d_mc_mutheta,  mutheta,  mc_object, cvweight,  __is_newflux);
            fill_vert_error_bands(mnvh1d_mc_enu,      enu,      mc_object, cvweight,  __is_newflux);
            fill_vert_error_bands(mnvh1d_mc_q2,       q2,       mc_object, cvweight,  __is_newflux);
            fill_vert_error_bands(mnvh1d_mc_w,        w,        mc_object, cvweight,  __is_newflux);
            
            fill_vert_error_bands(mnvh1d_gn_costheta, costheta0, mc_object, cvweight,  __is_newflux);
            fill_vert_error_bands(mnvh1d_gn_theta,    theta0,    mc_object, cvweight,  __is_newflux);
            fill_vert_error_bands(mnvh1d_gn_pienergy, pienergy0, mc_object, cvweight,  __is_newflux);
            fill_vert_error_bands(mnvh1d_gn_pimom,    pimom0,    mc_object, cvweight,  __is_newflux);
            fill_vert_error_bands(mnvh1d_gn_kinetic,  kinetic0,  mc_object, cvweight,  __is_newflux);
            fill_vert_error_bands(mnvh1d_gn_mmom,     mmom0,     mc_object, cvweight,  __is_newflux);
            fill_vert_error_bands(mnvh1d_gn_pz,       pz0,       mc_object, cvweight,  __is_newflux);
            fill_vert_error_bands(mnvh1d_gn_pT,       pT0,       mc_object, cvweight,  __is_newflux);
            fill_vert_error_bands(mnvh1d_gn_mutheta,  mutheta0,  mc_object, cvweight,  __is_newflux);
            fill_vert_error_bands(mnvh1d_gn_enu,      enu0 ,     mc_object, cvweight,  __is_newflux);
            fill_vert_error_bands(mnvh1d_gn_q2,       q20,       mc_object, cvweight,  __is_newflux);
            fill_vert_error_bands(mnvh1d_gn_w,        w0,        mc_object, cvweight,  __is_newflux);
            
            if (signal) {
                            
                mnvh1d_mcsig_costheta->Fill(costheta, cvweight);
                mnvh1d_mcsig_theta->Fill(theta,       cvweight);
                mnvh1d_mcsig_pienergy->Fill(pienergy, cvweight);
                mnvh1d_mcsig_pimom->Fill(pimom,       cvweight);
                mnvh1d_mcsig_kinetic->Fill(kinetic,   cvweight);
                mnvh1d_mcsig_mmom->Fill(mmom,         cvweight);
                mnvh1d_mcsig_pz->Fill(pz,             cvweight);
                mnvh1d_mcsig_pT->Fill(pT,             cvweight);
                mnvh1d_mcsig_mutheta->Fill(mutheta,   cvweight);
                mnvh1d_mcsig_enu->Fill(enu,           cvweight);
                mnvh1d_mcsig_q2->Fill(q2,             cvweight);
                mnvh1d_mcsig_w->Fill(w,               cvweight);
                
                mnvh1d_gnsig_costheta->Fill(costheta0, cvweight);
                mnvh1d_gnsig_theta->Fill(theta0,       cvweight);
                mnvh1d_gnsig_pienergy->Fill(pienergy0, cvweight);
                mnvh1d_gnsig_pimom->Fill(pimom0,       cvweight);
                mnvh1d_gnsig_kinetic->Fill(kinetic0,   cvweight);
                mnvh1d_gnsig_mmom->Fill(mmom0,         cvweight);
                mnvh1d_gnsig_pz->Fill(pz0,             cvweight);
                mnvh1d_gnsig_pT->Fill(pT0,             cvweight);
                mnvh1d_gnsig_mutheta->Fill(mutheta0,   cvweight);
                mnvh1d_gnsig_enu->Fill(enu0,           cvweight);
                mnvh1d_gnsig_q2->Fill(q20,             cvweight);
                mnvh1d_gnsig_w->Fill(w0,               cvweight);

                
                fill_vert_error_bands(mnvh1d_mcsig_costheta, costheta, mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_mcsig_theta,    theta,    mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_mcsig_pienergy, pienergy, mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_mcsig_pimom,    pimom,    mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_mcsig_kinetic,  kinetic,  mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_mcsig_mmom,     mmom,     mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_mcsig_pz,       pz,       mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_mcsig_pT,       pT,       mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_mcsig_mutheta,  mutheta,  mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_mcsig_enu,      enu,      mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_mcsig_q2,       q2,       mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_mcsig_w,        w,        mc_object, cvweight,  __is_newflux);
                
                fill_vert_error_bands(mnvh1d_gnsig_costheta, costheta0, mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_gnsig_theta,    theta0,    mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_gnsig_pienergy, pienergy0, mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_gnsig_pimom,    pimom0,    mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_gnsig_kinetic,  kinetic0,  mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_gnsig_mmom,     mmom0,     mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_gnsig_pz,       pz0,       mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_gnsig_pT,       pT0,       mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_gnsig_mutheta,  mutheta0,  mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_gnsig_enu,      enu0,      mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_gnsig_q2,       q20,       mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_gnsig_w,        w0,        mc_object, cvweight,  __is_newflux);
                
            }

            if (background) {
                
                mnvh1d_mcbkg_costheta->Fill(costheta, cvweight);
                mnvh1d_mcbkg_theta->Fill(theta,       cvweight);
                mnvh1d_mcbkg_pienergy->Fill(pienergy, cvweight);
                mnvh1d_mcbkg_pimom->Fill(pimom,       cvweight);
                mnvh1d_mcbkg_kinetic->Fill(kinetic,   cvweight);
                mnvh1d_mcbkg_mmom->Fill(mmom,         cvweight);
                mnvh1d_mcbkg_pz->Fill(pz,             cvweight);
                mnvh1d_mcbkg_pT->Fill(pT,             cvweight);
                mnvh1d_mcbkg_mutheta->Fill(mutheta,   cvweight);
                mnvh1d_mcbkg_enu->Fill(enu,           cvweight);
                mnvh1d_mcbkg_q2->Fill(q2,             cvweight);
                mnvh1d_mcbkg_w->Fill(w,               cvweight);
            
                mnvh1d_gnbkg_costheta->Fill(costheta0, cvweight);
                mnvh1d_gnbkg_theta->Fill(theta0,       cvweight);
                mnvh1d_gnbkg_pienergy->Fill(pienergy0, cvweight);
                mnvh1d_gnbkg_pimom->Fill(pimom0,       cvweight);
                mnvh1d_gnbkg_kinetic->Fill(kinetic0,   cvweight);
                mnvh1d_gnbkg_mmom->Fill(mmom0,         cvweight);
                mnvh1d_gnbkg_pz->Fill(pz0,             cvweight);
                mnvh1d_gnbkg_pT->Fill(pT0,             cvweight);
                mnvh1d_gnbkg_mutheta->Fill(mutheta0,   cvweight);
                mnvh1d_gnbkg_enu->Fill(enu0,           cvweight);
                mnvh1d_gnbkg_q2->Fill(q20,             cvweight);
                mnvh1d_gnbkg_w->Fill(w0,               cvweight);

                    
                fill_vert_error_bands(mnvh1d_mcbkg_costheta, costheta, mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_mcbkg_theta,    theta,    mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_mcbkg_pienergy, pienergy, mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_mcbkg_pimom,    pimom,    mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_mcbkg_kinetic,  kinetic,  mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_mcbkg_mmom,     mmom,     mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_mcbkg_pz,       pz,       mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_mcbkg_pT,       pT,       mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_mcbkg_mutheta,  mutheta,  mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_mcbkg_enu,      enu,      mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_mcbkg_q2,       q2,       mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_mcbkg_w,        w,        mc_object, cvweight,  __is_newflux);
                            
                fill_vert_error_bands(mnvh1d_gnbkg_costheta, costheta0, mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_gnbkg_theta,    theta0,    mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_gnbkg_pienergy, pienergy0, mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_gnbkg_pimom,    pimom0,    mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_gnbkg_kinetic,  kinetic0,  mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_gnbkg_mmom,     mmom0,     mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_gnbkg_pz,       pz0,       mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_gnbkg_pT,       pT0,       mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_gnbkg_mutheta,  mutheta0,  mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_gnbkg_enu,      enu0,      mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_gnbkg_q2,       q20,       mc_object, cvweight,  __is_newflux);
                fill_vert_error_bands(mnvh1d_gnbkg_w,        w0,        mc_object, cvweight,  __is_newflux);
                
            }
            
        }

            // loop over lateral-error band universes
        for (int i = 0; i < constants().n_lateral_universe; ++i) {
            double em_uncertainty = 0.022;                                // em-response uncertainty
            double mr_uncertainty = mc_object->muon_dp/mc_object->muon_p; // muon-response uncertainty
            std::vector<double> em_random_shifts = random_shifts::Get().get_random_shifts<15>(em_uncertainty); // ~ Gaussian(0.0, em_uncertainty)
            std::vector<double> mr_random_shifts = random_shifts::Get().get_random_shifts<25>(mr_uncertainty); // ~ Gaussian(0.0, mr_uncertainty)
            
            double enu_i[3] = { 0.0 }; // actually the [1] and [2] elements are not used in this program
            double q2_i[3]  = { 0.0 }; // same
            double w_i[3]   = { 0.0 }; // same
            
            double em_shift = em_random_shifts[i];
            double mr_shift = mr_random_shifts[i];

            TLorentzVector pion_p4_i = (1.0 + em_shift) * pion_p4;
            TLorentzVector muon_p4_i = (1.0 + mr_shift) * muon_p4;
        
            if (__calorimetry) {
                cc1pi0::kincalc calc1(muon_p4_i, pion_p4_i, alpha,   // Vary both muon and pion 4-momenta.
                                      true,                          // The results are used to check if the event passes
                                      mc_object->Vertex_blob_energy, // the enu and w cuts in the universe, but not used to fill histograms
                                      mc_object->Dispersed_blob_energy); 
                
                enu_i[0] = calc1.enu()/GeV;
                q2_i[0]  = calc1.q2()/(GeV*GeV);
                w_i[0]   = calc1.w()/GeV;
            
                cc1pi0::kincalc calc2(muon_p4_i, pion_p4, alpha,   // Vary only the muon 4-momentum.
                                      true,                        // The results has the effect of the shifted muon 4-momentum
                                      mc_object->Vertex_blob_energy,
                                      mc_object->Dispersed_blob_energy); 
                
                enu_i[1] = calc2.enu()/GeV;
                q2_i[1]  = calc2.q2()/(GeV*GeV);
                w_i[1]   = calc2.w()/GeV;
            
                cc1pi0::kincalc calc3(muon_p4, pion_p4_i, alpha,   // Vary only the pion 4-momentum
                                      true,                        // The results has the effect of the shifted pion 4-momentum
                                      mc_object->Vertex_blob_energy,
                                      mc_object->Dispersed_blob_energy); 
                
                enu_i[2] = calc3.enu()/GeV;
                q2_i[2]  = calc3.q2()/(GeV*GeV);
                w_i[2]   = calc3.w()/GeV;
                
            
            }

            if (__resonant_formula) {
                cc1pi0::kincalc calc1(muon_p4_i, pion_p4_i, alpha); // see comment above in __calorimetry
                
                enu_i[0] = calc1.enu()/GeV;
                q2_i[0]  = calc1.q2()/(GeV*GeV);
                w_i[0]   = calc1.w()/GeV;
            
                cc1pi0::kincalc calc2(muon_p4_i, pion_p4, alpha);   // see comment above in __calorimetry
                
                enu_i[1] = calc2.enu()/GeV;
                q2_i[1]  = calc2.q2()/(GeV*GeV);
                w_i[1]   = calc2.w()/GeV;
                
                cc1pi0::kincalc calc3(muon_p4, pion_p4_i, alpha);   // see comment above in __calorimetry
                
                enu_i[2] = calc3.enu()/GeV;
                q2_i[2]  = calc3.q2()/(GeV*GeV);
                w_i[2]   = calc3.w()/GeV;
         
            }

                   // Muon response systematics
            bool not_within_enu_range_mr_i = (enu_i[1] < __emin) || (__emax < enu_i[1]);
            bool not_within_w_range_mr_i   = (w_i[1] < __wmin)   || (__wmax < w_i[1]);
            
            bool pass_enu_cut_mr_i = true;
            if (__enu_cut && not_within_enu_range_mr_i) pass_enu_cut_mr_i = false;

            bool pass_w_cut_mr_i = true;
            if (__w_cut && not_within_w_range_mr_i) pass_w_cut_mr_i = false;


                // EM energy scale systematic
            bool not_within_enu_range_em_i = (enu_i[2] < __emin) || (__emax < enu_i[2]);
            bool not_within_w_range_em_i   = (w_i[2] < __wmin)   || (__wmax < w_i[2]);
            
            bool pass_enu_cut_em_i = true;
            if (__enu_cut && not_within_enu_range_em_i) pass_enu_cut_em_i = false;

            bool pass_w_cut_em_i = true;
            if (__w_cut && not_within_w_range_em_i) pass_w_cut_em_i = false;

            double mgg_i = (1.0 + em_shift) * mgg;
            bool pass_mgg_cut_em_i = (lower_peak < mgg_i) && (mgg_i < upper_peak);
            
            if (pass_enu_cut_mr_i && pass_w_cut_mr_i) {
                    // Only the event population change
                fill_lateral_universe_histogram(mnvh1d_mc_mgg,  "Muon_momentum",   i, mgg,            cvweight);
                fill_lateral_universe_histogram(mnvh1d_mc_mgg2, "Muon_momentum",   i, mc_object->mgg, cvweight);
                
                if (signal) {

                    fill_lateral_universe_histogram(mnvh1d_mcsig_mgg,  "Muon_momentum",   i, mgg,            cvweight);
                    fill_lateral_universe_histogram(mnvh1d_mcsig_mgg2, "Muon_momentum",   i, mc_object->mgg, cvweight);
                }

                if (background) {

                    fill_lateral_universe_histogram(mnvh1d_mcbkg_mgg,  "Muon_momentum",   i, mgg,            cvweight);
                    fill_lateral_universe_histogram(mnvh1d_mcbkg_mgg2, "Muon_momentum",   i, mc_object->mgg, cvweight);
                }

            }

            if (pass_enu_cut_em_i && pass_w_cut_em_i) {
                    // Both the event population and value change
                fill_lateral_universe_histogram(mnvh1d_mc_mgg,  "EM_energy_scale", i, mgg_i,          cvweight);
                fill_lateral_universe_histogram(mnvh1d_mc_mgg2, "EM_energy_scale", i, mc_object->mgg, cvweight); // no need to handle correctly this histogram
                
                if (signal) {

                    fill_lateral_universe_histogram(mnvh1d_mcsig_mgg,  "EM_energy_scale", i, mgg_i,          cvweight);
                    fill_lateral_universe_histogram(mnvh1d_mcsig_mgg2, "EM_energy_scale", i, mc_object->mgg, cvweight); // same
                }

                if (background) {

                    fill_lateral_universe_histogram(mnvh1d_mcbkg_mgg,  "EM_energy_scale", i, mgg_i,          cvweight);
                    fill_lateral_universe_histogram(mnvh1d_mcbkg_mgg2, "EM_energy_scale", i, mc_object->mgg, cvweight); // same
                
                }

            }

                // Muon momentum systematics
            if (pass_enu_cut_mr_i && pass_w_cut_mr_i && pass_mgg_cut_cv) {

                double mmom_i = (1.0 + mr_shift) * mmom;
                double pz_i   = (1.0 + mr_shift) * pz;
                double pT_i   = (1.0 + mr_shift) * pT;
                
                fill_lateral_universe_histogram(mnvh1d_mc_costheta, "Muon_momentum", i, costheta,  cvweight); // - for variables that does not depend on the 
                fill_lateral_universe_histogram(mnvh1d_mc_theta,    "Muon_momentum", i, theta,     cvweight); //    muon momentum, fill the cv value
                fill_lateral_universe_histogram(mnvh1d_mc_pienergy, "Muon_momentum", i, pienergy,  cvweight);
                fill_lateral_universe_histogram(mnvh1d_mc_pimom,    "Muon_momentum", i, pimom,     cvweight);
                fill_lateral_universe_histogram(mnvh1d_mc_kinetic,  "Muon_momentum", i, kinetic,   cvweight);
                fill_lateral_universe_histogram(mnvh1d_mc_mmom,     "Muon_momentum", i, mmom_i,    cvweight); // - for variables that depend on muon momentum
                fill_lateral_universe_histogram(mnvh1d_mc_pz,       "Muon_momentum", i, pz_i,      cvweight); //    explicity, fill the shifted momentum
                fill_lateral_universe_histogram(mnvh1d_mc_pT,       "Muon_momentum", i, pT_i,      cvweight);
                fill_lateral_universe_histogram(mnvh1d_mc_mutheta,  "Muon_momentum", i, mutheta,   cvweight);
                fill_lateral_universe_histogram(mnvh1d_mc_enu,      "Muon_momentum", i, enu_i[1],  cvweight); // - for variables that depend on the muon momentum
                fill_lateral_universe_histogram(mnvh1d_mc_q2,       "Muon_momentum", i, q2_i[1],   cvweight); //    implicitly, fill the shifted value due to 
                fill_lateral_universe_histogram(mnvh1d_mc_w,        "Muon_momentum", i, w_i[1],    cvweight); //    the shift in the muon momentum

                if (signal) {

                    fill_lateral_universe_histogram(mnvh1d_mcsig_costheta, "Muon_momentum", i, costheta,  cvweight); // - for variables that does not depend on the 
                    fill_lateral_universe_histogram(mnvh1d_mcsig_theta,    "Muon_momentum", i, theta,     cvweight); //    muon momentum, fill the cv value
                    fill_lateral_universe_histogram(mnvh1d_mcsig_pienergy, "Muon_momentum", i, pienergy,  cvweight);
                    fill_lateral_universe_histogram(mnvh1d_mcsig_pimom,    "Muon_momentum", i, pimom,     cvweight);
                    fill_lateral_universe_histogram(mnvh1d_mcsig_kinetic,  "Muon_momentum", i, kinetic,   cvweight);
                    fill_lateral_universe_histogram(mnvh1d_mcsig_mmom,     "Muon_momentum", i, mmom_i,    cvweight); // - for variables that depend on muon momentum
                    fill_lateral_universe_histogram(mnvh1d_mcsig_pz,       "Muon_momentum", i, pz_i,      cvweight); //    explicity, fill the shifted momentum
                    fill_lateral_universe_histogram(mnvh1d_mcsig_pT,       "Muon_momentum", i, pT_i,      cvweight);
                    fill_lateral_universe_histogram(mnvh1d_mcsig_mutheta,  "Muon_momentum", i, mutheta,   cvweight);
                    fill_lateral_universe_histogram(mnvh1d_mcsig_enu,      "Muon_momentum", i, enu_i[1],  cvweight); // - for variables that depend on the muon momentum
                    fill_lateral_universe_histogram(mnvh1d_mcsig_q2,       "Muon_momentum", i, q2_i[1],   cvweight); //    implicitly, fill the shifted value due to 
                    fill_lateral_universe_histogram(mnvh1d_mcsig_w,        "Muon_momentum", i, w_i[1],    cvweight); //    the shift in the muon momentum
                    
                }

                if (background) {

                    fill_lateral_universe_histogram(mnvh1d_mcbkg_costheta, "Muon_momentum", i, costheta,  cvweight); // - for variables that does not depend on the 
                    fill_lateral_universe_histogram(mnvh1d_mcbkg_theta,    "Muon_momentum", i, theta,     cvweight); //    muon momentum, fill the cv value
                    fill_lateral_universe_histogram(mnvh1d_mcbkg_pienergy, "Muon_momentum", i, pienergy,  cvweight);
                    fill_lateral_universe_histogram(mnvh1d_mcbkg_pimom,    "Muon_momentum", i, pimom,     cvweight);
                    fill_lateral_universe_histogram(mnvh1d_mcbkg_kinetic,  "Muon_momentum", i, kinetic,   cvweight);
                    fill_lateral_universe_histogram(mnvh1d_mcbkg_mmom,     "Muon_momentum", i, mmom_i,    cvweight); // - for variables that depend on muon momentum
                    fill_lateral_universe_histogram(mnvh1d_mcbkg_pz,       "Muon_momentum", i, pz_i,      cvweight); //    explicity, fill the shifted momentum
                    fill_lateral_universe_histogram(mnvh1d_mcbkg_pT,       "Muon_momentum", i, pT_i,      cvweight);
                    fill_lateral_universe_histogram(mnvh1d_mcbkg_mutheta,  "Muon_momentum", i, mutheta,   cvweight);
                    fill_lateral_universe_histogram(mnvh1d_mcbkg_enu,      "Muon_momentum", i, enu_i[1],  cvweight); // - for variables that depend on the muon momentum
                    fill_lateral_universe_histogram(mnvh1d_mcbkg_q2,       "Muon_momentum", i, q2_i[1],   cvweight); //    implicitly, fill the shifted value due to 
                    fill_lateral_universe_histogram(mnvh1d_mcbkg_w,        "Muon_momentum", i, w_i[1],    cvweight); //    the shift in the muon momentum
                               
                }


            }

                // EM energy scale systematics
            if (pass_enu_cut_em_i && pass_w_cut_em_i && pass_mgg_cut_em_i) {

                double pienergy_i = (1.0 + em_shift) * pienergy;
                double pimom_i    = (1.0 + em_shift) * pimom;
                double kinetic_i  = (1.0 + em_shift) * kinetic;
                
                fill_lateral_universe_histogram(mnvh1d_mc_costheta, "EM_energy_scale", i, costheta,    cvweight); // - for variables that does not depend on the 
                fill_lateral_universe_histogram(mnvh1d_mc_theta,    "EM_energy_scale", i, theta,       cvweight); //    em energy scale, fill the cv value
                fill_lateral_universe_histogram(mnvh1d_mc_pienergy, "EM_energy_scale", i, pienergy_i,  cvweight); // - for variables that depend on em energy scale
                fill_lateral_universe_histogram(mnvh1d_mc_pimom,    "EM_energy_scale", i, pimom_i,     cvweight); //    explicity, fill the shifted value
                fill_lateral_universe_histogram(mnvh1d_mc_kinetic,  "EM_energy_scale", i, kinetic_i,   cvweight);
                fill_lateral_universe_histogram(mnvh1d_mc_mmom,     "EM_energy_scale", i, mmom,        cvweight);
                fill_lateral_universe_histogram(mnvh1d_mc_pz,       "EM_energy_scale", i, pz,          cvweight);
                fill_lateral_universe_histogram(mnvh1d_mc_pT,       "EM_energy_scale", i, pT,          cvweight);
                fill_lateral_universe_histogram(mnvh1d_mc_mutheta,  "EM_energy_scale", i, mutheta,     cvweight);
                fill_lateral_universe_histogram(mnvh1d_mc_enu,      "EM_energy_scale", i, enu_i[2],    cvweight); // - for variables that depend on the muon momentum
                fill_lateral_universe_histogram(mnvh1d_mc_q2,       "EM_energy_scale", i, q2_i[2],     cvweight); //   implicitly, fill the shifted value due to 
                fill_lateral_universe_histogram(mnvh1d_mc_w,        "EM_energy_scale", i, w_i[2],      cvweight); //   the shift in the muon momentum
        

                if (signal) {

                    fill_lateral_universe_histogram(mnvh1d_mcsig_costheta, "EM_energy_scale", i, costheta,    cvweight); // - for variables that does not depend on the 
                    fill_lateral_universe_histogram(mnvh1d_mcsig_theta,    "EM_energy_scale", i, theta,       cvweight); //    em energy scale, fill the cv value
                    fill_lateral_universe_histogram(mnvh1d_mcsig_pienergy, "EM_energy_scale", i, pienergy_i,  cvweight); // - for variables that depend on em energy scale
                    fill_lateral_universe_histogram(mnvh1d_mcsig_pimom,    "EM_energy_scale", i, pimom_i,     cvweight); //    explicity, fill the shifted value
                    fill_lateral_universe_histogram(mnvh1d_mcsig_kinetic,  "EM_energy_scale", i, kinetic_i,   cvweight);
                    fill_lateral_universe_histogram(mnvh1d_mcsig_mmom,     "EM_energy_scale", i, mmom,        cvweight);
                    fill_lateral_universe_histogram(mnvh1d_mcsig_pz,       "EM_energy_scale", i, pz,          cvweight);
                    fill_lateral_universe_histogram(mnvh1d_mcsig_pT,       "EM_energy_scale", i, pT,          cvweight);
                    fill_lateral_universe_histogram(mnvh1d_mcsig_mutheta,  "EM_energy_scale", i, mutheta,     cvweight);
                    fill_lateral_universe_histogram(mnvh1d_mcsig_enu,      "EM_energy_scale", i, enu_i[2],    cvweight); // - for variables that depend on the muon momentum
                    fill_lateral_universe_histogram(mnvh1d_mcsig_q2,       "EM_energy_scale", i, q2_i[2],     cvweight); //   implicitly, fill the shifted value due to 
                    fill_lateral_universe_histogram(mnvh1d_mcsig_w,        "EM_energy_scale", i, w_i[2],      cvweight); //   the shift in the muon momentum
                                     
                }

                if (background) {

                    fill_lateral_universe_histogram(mnvh1d_mcbkg_costheta, "EM_energy_scale", i, costheta,    cvweight); // - for variables that does not depend on the 
                    fill_lateral_universe_histogram(mnvh1d_mcbkg_theta,    "EM_energy_scale", i, theta,       cvweight); //    em energy scale, fill the cv value
                    fill_lateral_universe_histogram(mnvh1d_mcbkg_pienergy, "EM_energy_scale", i, pienergy_i,  cvweight); // - for variables that depend on em energy scale
                    fill_lateral_universe_histogram(mnvh1d_mcbkg_pimom,    "EM_energy_scale", i, pimom_i,     cvweight); //    explicity, fill the shifted value
                    fill_lateral_universe_histogram(mnvh1d_mcbkg_kinetic,  "EM_energy_scale", i, kinetic_i,   cvweight);
                    fill_lateral_universe_histogram(mnvh1d_mcbkg_mmom,     "EM_energy_scale", i, mmom,        cvweight);
                    fill_lateral_universe_histogram(mnvh1d_mcbkg_pz,       "EM_energy_scale", i, pz,          cvweight);
                    fill_lateral_universe_histogram(mnvh1d_mcbkg_pT,       "EM_energy_scale", i, pT,          cvweight);
                    fill_lateral_universe_histogram(mnvh1d_mcbkg_mutheta,  "EM_energy_scale", i, mutheta,     cvweight);
                    fill_lateral_universe_histogram(mnvh1d_mcbkg_enu,      "EM_energy_scale", i, enu_i[2],    cvweight); // - for variables that depend on the muon momentum
                    fill_lateral_universe_histogram(mnvh1d_mcbkg_q2,       "EM_energy_scale", i, q2_i[2],     cvweight); //   implicitly, fill the shifted value due to 
                    fill_lateral_universe_histogram(mnvh1d_mcbkg_w,        "EM_energy_scale", i, w_i[2],      cvweight); //   the shift in the muon momentum
                               
                }


            }

        
        } // loop over universes

        
    }

    std::cout << "Selected events: " << selectedEvent << std::endl;
    
    delete mc_object; mc_object = 0;
    
    MnvH1D* mnvh1d_data_mgg  = new MnvH1D(*(make_hist_wrapper_1d(0, Form("%s-data_mgg", __detector.c_str()))));

        // costheta
    MnvH1D* mnvh1d_data_costheta  = new MnvH1D(*make_hist_wrapper_1d(1, Form("%s-data_costheta_masspeak",  __detector.c_str())));
    MnvH1D* mnvh1d_mdsig_costheta = new MnvH1D(*make_hist_wrapper_1d(1, Form("%s-mdsig_costheta_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_mdbkg_costheta = new MnvH1D(*make_hist_wrapper_1d(1, Form("%s-mdbkg_costheta_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gd_costheta    = new MnvH1D(*make_hist_wrapper_1d(1, Form("%s-gd_costheta_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_gdsig_costheta = new MnvH1D(*make_hist_wrapper_1d(1, Form("%s-gdsig_costheta_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gdbkg_costheta = new MnvH1D(*make_hist_wrapper_1d(1, Form("%s-gdbkg_costheta_masspeak", __detector.c_str())));
    
        // theta
    MnvH1D* mnvh1d_data_theta  = new MnvH1D(*make_hist_wrapper_1d(5, Form("%s-data_theta_masspeak",  __detector.c_str())));
    MnvH1D* mnvh1d_mdsig_theta = new MnvH1D(*make_hist_wrapper_1d(5, Form("%s-mdsig_theta_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_mdbkg_theta = new MnvH1D(*make_hist_wrapper_1d(5, Form("%s-mdbkg_theta_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gd_theta    = new MnvH1D(*make_hist_wrapper_1d(5, Form("%s-gd_theta_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_gdsig_theta = new MnvH1D(*make_hist_wrapper_1d(5, Form("%s-gdsig_theta_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gdbkg_theta = new MnvH1D(*make_hist_wrapper_1d(5, Form("%s-gdbkg_theta_masspeak", __detector.c_str())));

        // pi0 total energy
    MnvH1D* mnvh1d_data_pienergy  = new MnvH1D(*make_hist_wrapper_1d(10, Form("%s-data_pienergy_masspeak",  __detector.c_str())));
    MnvH1D* mnvh1d_mdsig_pienergy = new MnvH1D(*make_hist_wrapper_1d(10, Form("%s-mdsig_pienergy_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_mdbkg_pienergy = new MnvH1D(*make_hist_wrapper_1d(10, Form("%s-mdbkg_pienergy_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gd_pienergy    = new MnvH1D(*make_hist_wrapper_1d(10, Form("%s-gd_pienergy_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_gdsig_pienergy = new MnvH1D(*make_hist_wrapper_1d(10, Form("%s-gdsig_pienergy_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gdbkg_pienergy = new MnvH1D(*make_hist_wrapper_1d(10, Form("%s-gdbkg_pienergy_masspeak", __detector.c_str())));

        // pi0 momentum
    MnvH1D* mnvh1d_data_pimom  = new MnvH1D(*make_hist_wrapper_1d(15, Form("%s-data_pimom_masspeak",  __detector.c_str())));
    MnvH1D* mnvh1d_mdsig_pimom = new MnvH1D(*make_hist_wrapper_1d(15, Form("%s-mdsig_pimom_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_mdbkg_pimom = new MnvH1D(*make_hist_wrapper_1d(15, Form("%s-mdbkg_pimom_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gd_pimom    = new MnvH1D(*make_hist_wrapper_1d(15, Form("%s-gd_pimom_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_gdsig_pimom = new MnvH1D(*make_hist_wrapper_1d(15, Form("%s-gdsig_pimom_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gdbkg_pimom = new MnvH1D(*make_hist_wrapper_1d(15, Form("%s-gdbkg_pimom_masspeak", __detector.c_str())));

        // pi0 kinetic energy
    MnvH1D* mnvh1d_data_kinetic  = new MnvH1D(*make_hist_wrapper_1d(20, Form("%s-data_kinetic_masspeak",  __detector.c_str())));
    MnvH1D* mnvh1d_mdsig_kinetic = new MnvH1D(*make_hist_wrapper_1d(20, Form("%s-mdsig_kinetic_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_mdbkg_kinetic = new MnvH1D(*make_hist_wrapper_1d(20, Form("%s-mdbkg_kinetic_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gd_kinetic    = new MnvH1D(*make_hist_wrapper_1d(20, Form("%s-gd_kinetic_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_gdsig_kinetic = new MnvH1D(*make_hist_wrapper_1d(20, Form("%s-gdsig_kinetic_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gdbkg_kinetic = new MnvH1D(*make_hist_wrapper_1d(20, Form("%s-gdbkg_kinetic_masspeak", __detector.c_str())));

        // muon momentum
    MnvH1D* mnvh1d_data_mmom  = new MnvH1D(*make_hist_wrapper_1d(25, Form("%s-data_mmom_masspeak",  __detector.c_str())));
    MnvH1D* mnvh1d_mdsig_mmom = new MnvH1D(*make_hist_wrapper_1d(25, Form("%s-mdsig_mmom_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_mdbkg_mmom = new MnvH1D(*make_hist_wrapper_1d(25, Form("%s-mdbkg_mmom_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gd_mmom    = new MnvH1D(*make_hist_wrapper_1d(25, Form("%s-gd_mmom_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_gdsig_mmom = new MnvH1D(*make_hist_wrapper_1d(25, Form("%s-gdsig_mmom_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gdbkg_mmom = new MnvH1D(*make_hist_wrapper_1d(25, Form("%s-gdbkg_mmom_masspeak", __detector.c_str())));

        // muon longitudinal momentum
    MnvH1D* mnvh1d_data_pz  = new MnvH1D(*make_hist_wrapper_1d(50, Form("%s-data_pz_masspeak",  __detector.c_str())));
    MnvH1D* mnvh1d_mdsig_pz = new MnvH1D(*make_hist_wrapper_1d(50, Form("%s-mdsig_pz_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_mdbkg_pz = new MnvH1D(*make_hist_wrapper_1d(50, Form("%s-mdbkg_pz_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gd_pz    = new MnvH1D(*make_hist_wrapper_1d(50, Form("%s-gd_pz_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_gdsig_pz = new MnvH1D(*make_hist_wrapper_1d(50, Form("%s-gdsig_pz_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gdbkg_pz = new MnvH1D(*make_hist_wrapper_1d(50, Form("%s-gdbkg_pz_masspeak", __detector.c_str())));

        // muon transverse momentum
    MnvH1D* mnvh1d_data_pT  = new MnvH1D(*make_hist_wrapper_1d(55, Form("%s-data_pT_masspeak",  __detector.c_str())));
    MnvH1D* mnvh1d_mdsig_pT = new MnvH1D(*make_hist_wrapper_1d(55, Form("%s-mdsig_pT_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_mdbkg_pT = new MnvH1D(*make_hist_wrapper_1d(55, Form("%s-mdbkg_pT_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gd_pT    = new MnvH1D(*make_hist_wrapper_1d(55, Form("%s-gd_pT_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_gdsig_pT = new MnvH1D(*make_hist_wrapper_1d(55, Form("%s-gdsig_pT_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gdbkg_pT = new MnvH1D(*make_hist_wrapper_1d(55, Form("%s-gdbkg_pT_masspeak", __detector.c_str())));

        // muon theta
    MnvH1D* mnvh1d_data_mutheta  = new MnvH1D(*make_hist_wrapper_1d(30, Form("%s-data_mutheta_masspeak",  __detector.c_str())));
    MnvH1D* mnvh1d_mdsig_mutheta = new MnvH1D(*make_hist_wrapper_1d(30, Form("%s-mdsig_mutheta_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_mdbkg_mutheta = new MnvH1D(*make_hist_wrapper_1d(30, Form("%s-mdbkg_mutheta_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gd_mutheta    = new MnvH1D(*make_hist_wrapper_1d(30, Form("%s-gd_mutheta_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_gdsig_mutheta = new MnvH1D(*make_hist_wrapper_1d(30, Form("%s-gdsig_mutheta_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gdbkg_mutheta = new MnvH1D(*make_hist_wrapper_1d(30, Form("%s-gdbkg_mutheta_masspeak", __detector.c_str())));

        // Neutrino energy
    MnvH1D* mnvh1d_data_enu  = new MnvH1D(*make_hist_wrapper_1d(35, Form("%s-data_enu_masspeak",  __detector.c_str())));
    MnvH1D* mnvh1d_mdsig_enu = new MnvH1D(*make_hist_wrapper_1d(35, Form("%s-mdsig_enu_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_mdbkg_enu = new MnvH1D(*make_hist_wrapper_1d(35, Form("%s-mdbkg_enu_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gd_enu    = new MnvH1D(*make_hist_wrapper_1d(35, Form("%s-gd_enu_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_gdsig_enu = new MnvH1D(*make_hist_wrapper_1d(35, Form("%s-gdsig_enu_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gdbkg_enu = new MnvH1D(*make_hist_wrapper_1d(35, Form("%s-gdbkg_enu_masspeak", __detector.c_str())));

        // Q^2
    MnvH1D* mnvh1d_data_q2  = new MnvH1D(*make_hist_wrapper_1d(40, Form("%s-data_q2_masspeak",  __detector.c_str())));
    MnvH1D* mnvh1d_mdsig_q2 = new MnvH1D(*make_hist_wrapper_1d(40, Form("%s-mdsig_q2_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_mdbkg_q2 = new MnvH1D(*make_hist_wrapper_1d(40, Form("%s-mdbkg_q2_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gd_q2    = new MnvH1D(*make_hist_wrapper_1d(40, Form("%s-gd_q2_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_gdsig_q2 = new MnvH1D(*make_hist_wrapper_1d(40, Form("%s-gdsig_q2_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gdbkg_q2 = new MnvH1D(*make_hist_wrapper_1d(40, Form("%s-gdbkg_q2_masspeak", __detector.c_str())));

        // Hadronic invariant mass
    MnvH1D* mnvh1d_data_w  = new MnvH1D(*make_hist_wrapper_1d(45, Form("%s-data_w_masspeak",  __detector.c_str())));
    MnvH1D* mnvh1d_mdsig_w = new MnvH1D(*make_hist_wrapper_1d(45, Form("%s-mdsig_w_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_mdbkg_w = new MnvH1D(*make_hist_wrapper_1d(45, Form("%s-mdbkg_w_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gd_w    = new MnvH1D(*make_hist_wrapper_1d(45, Form("%s-gd_w_masspeak",    __detector.c_str())));
    MnvH1D* mnvh1d_gdsig_w = new MnvH1D(*make_hist_wrapper_1d(45, Form("%s-gdsig_w_masspeak", __detector.c_str())));
    MnvH1D* mnvh1d_gdbkg_w = new MnvH1D(*make_hist_wrapper_1d(45, Form("%s-gdbkg_w_masspeak", __detector.c_str())));


        //---------------------- PROCESS DATA --------------------------------------------------------------

    std::cout << "Creating DATA chain " << std::endl;
    TChain* data_chain = make_chain(__data_file.c_str());
    const int n_data = data_chain->GetEntries();

    std::cout << "Total data entries: " << n_data << std::endl;
    
    mc_struct* data_object = new mc_struct; // important change so that cv weight can be applied to
    data_object->Init(data_chain);          // the fake data

    std::ofstream ofs("data-signal-candidates-file.dat");

    TH1D* h_warping_function_pimom(NULL);
    TH1D* h_warping_function_theta(NULL);
    
    TFile warping_file("warping-function.root");
    
    if (__unfold_test) {
        
        h_warping_function_pimom = (TH1D*) warping_file.Get("pimom");
        h_warping_function_theta = (TH1D*) warping_file.Get("theta");
        
        assert(h_warping_function_pimom);
        assert(h_warping_function_theta);
        
    }
        
    std::cout << "Fill data histograms..." << std::endl;
    readEntry = 0;
    selectedEvent = 0;
    for (int entry = 0; entry < n_data; ++entry) {
        
        data_object->GetEntry(entry);

        ++readEntry;

        if (readEntry > __max_nevent) break;

        if (data_object->mgg > 1e4) continue;

        if (readEntry%1000 == 0) std::cout << "Processed " << (int)(100.0*readEntry/n_data) << " (%) "<< std::endl;
        
        
        do_count_data(0, data_counter);        
        if (data_object->survive_minos_match != 1) continue;

        do_count_data(1, data_counter);        
        if (data_object->minos_trk_qp <= 0) continue;        

        do_count_data(2, data_counter);        
        if (data_object->survive_fiducial != 1) continue;
        if (is_downstream && data_object->vtx[2] < 7200.0) continue;

        bool oneTrack1 = data_object->primary_multiplicity2 == 0;
        bool oneTrack2 = data_object->max_ntrack < 2;
        bool gammaLike = data_object->max_deviation < 50;

                
        do_count_data(3, data_counter);        
        if (!oneTrack1) continue;

        do_count_data(4, data_counter);        
        if (!oneTrack2) continue;

        do_count_data(5, data_counter);        
        if (!gammaLike) continue;

        do_count_data(6, data_counter);        
        if (data_object->ntgtevis >= 20.0) continue;

        do_count_data(7, data_counter);        
        if (data_object->otherevis <= 80.0) continue;

        do_count_data(8, data_counter);        
        if (data_object->otherevis >= 2000.0) continue;        

        do_count_data(9, data_counter);        

        bool goodGamma = data_object->is_GoodBlob1 && data_object->is_GoodBlob2;
        if (!goodGamma) continue;        
        do_count_data(10, data_counter);        


        if (data_object->final_blob_ncx[1]<2) continue;
        do_count_data(11, data_counter);        

        double g1mfp = calc_distance(data_object->vtx, data_object->RE_photon_vertex_1);
        double g2mfp = calc_distance(data_object->vtx, data_object->RE_photon_vertex_2);
        
        if (g1mfp < __mfp) continue;
        if (g2mfp < __mfp) continue;    
        do_count_data(12, data_counter);        

        
        TLorentzVector muon_p4(data_object->muon_px,
                               data_object->muon_py,
                               data_object->muon_pz,
                               data_object->muon_E);
        TLorentzVector muon_p4_rot = muon_p4;
        muon_p4_rot.RotateX(-3.3 * deg2rad);

        TLorentzVector pion_p4(data_object->pimom[0],
                               data_object->pimom[1],
                               data_object->pimom[2],
                               data_object->pienergy); 



        double alpha = m0/121.0; // same calibration as the EML fit
        if (__closure_test || __unfold_test) alpha = m0/115.0;
        
            // reconstructed variables
        const double costheta  = cos(data_object->pithetab * deg2rad);
        double theta     = data_object->pithetab;
        double pienergy  = alpha*data_object->pienergy/GeV;
        double pimom     = alpha*data_object->pimom[3]/GeV;
        double kinetic   = (alpha*data_object->pienergy - m0)/GeV;
        double mmom      = data_object->muon_p/GeV;
        double pz        = muon_p4_rot.Pz()/GeV;
        double pT        = muon_p4_rot.Pt()/GeV;
        double mutheta   = data_object->muon_theta * rad2deg;
        const double mgg = alpha*data_object->mgg; 

                      
        double enu       = 0.0;
        double q2        = 0.0;
        double w         = 0.0;
       
        if (__calorimetry) {
            
            cc1pi0::kincalc calc(muon_p4, pion_p4, alpha,
                                    true,
                                    data_object->Vertex_blob_energy,
                                    data_object->Dispersed_blob_energy); 
        
            enu       = calc.enu()/GeV; 
            q2        = calc.q2()/(GeV * GeV);
            w         = calc.w()/GeV;
       
        }
        
        if (__resonant_formula) {
            
            cc1pi0::kincalc calc(muon_p4, pion_p4, alpha);
        
            enu       = calc.enu()/GeV; 
            q2        = calc.q2()/(GeV * GeV);
            w         = calc.w()/GeV;
       
        }
      
      
        if (__enu_cut && (enu < __emin)) continue;
        if (__enu_cut && (enu > __emax)) continue;
        do_count_data(13, data_counter);        

        if (__w_cut && (w < __wmin)) continue;
        if (__w_cut && (w > __wmax)) continue;

        do_count_data(14, data_counter);        
        
        mnvh1d_data_mgg->Fill(mgg);        
        
        if (mgg < lower_peak || mgg > upper_peak) continue;
        
        ++selectedEvent;


        mnvh1d_data_costheta->Fill(costheta);
        mnvh1d_data_theta->Fill(theta);
        mnvh1d_data_pienergy->Fill(pienergy);
        mnvh1d_data_pimom->Fill(pimom);
        mnvh1d_data_kinetic->Fill(kinetic);
        mnvh1d_data_mmom->Fill(mmom);
        mnvh1d_data_pz->Fill(pz);
        mnvh1d_data_pT->Fill(pT);
        mnvh1d_data_mutheta->Fill(mutheta);
        mnvh1d_data_enu->Fill(enu);
        mnvh1d_data_q2->Fill(q2);
        mnvh1d_data_w->Fill(w);
        

            // For mockup data 
        const double costheta0  = cos(data_object->pithetab0 * deg2rad);
        double theta0    = data_object->pithetab0;
        double pimom0    = data_object->pimom0[3]/GeV;
        double pienergy0 = data_object->pienergy0/GeV;
        double kinetic0  = (data_object->pienergy0-m0)/GeV;
        double mmom0     = data_object->truth_fslepton_P; // already in GeV
        double pz0       = mmom0 * cos(data_object->truth_fslepton_theta);
        double pT0       = mmom0 * sin(data_object->truth_fslepton_theta);
        double mutheta0  = data_object->truth_fslepton_theta * rad2deg;
        double enu0      = data_object->mc_incomingE/GeV;
        double q20       = data_object->mc_Q2/(GeV * GeV);
        double w0        = data_object->mc_w/GeV;;
 


        
        mnvh1d_gd_costheta->Fill(costheta0);
        mnvh1d_gd_theta->Fill(theta0);
        mnvh1d_gd_pienergy->Fill(pienergy0);
        mnvh1d_gd_pimom->Fill(pimom0);
        mnvh1d_gd_kinetic->Fill(kinetic0);
        mnvh1d_gd_mmom->Fill(mmom0);
        mnvh1d_gd_pz->Fill(pz0);
        mnvh1d_gd_pT->Fill(pT0);
        mnvh1d_gd_mutheta->Fill(mutheta0);
        mnvh1d_gd_enu->Fill(enu0);
        mnvh1d_gd_q2->Fill(q20);
        mnvh1d_gd_w->Fill(w0);


            // applicable only for mock data, separate into true signal and background
        int mc_current = data_object->mc_current;
        int mc_incoming = data_object->mc_incoming;
        bool truth_is_cc1pi0 = data_object->truth_is_cc1pi0;
        bool truth_is_fiducial = data_object->truth_is_fiducial;
        bool enu_cut = (enu0 > __emin) && (enu0 < __emax);
        bool w_cut   = (w0   > __wmin) && (w0   < __wmax);
        
        bool signal = truth_is_fiducial && truth_is_cc1pi0 && mc_incoming == -14 && mc_current==1;
        if (__enu_cut) signal = signal && enu_cut;
        if (__w_cut)   signal = signal && w_cut;
        
        if (is_downstream) signal = signal && (data_object->mc_vtx[2] > 7200.0);
        bool background = !signal;

        
        if (signal) {
            
            mnvh1d_mdsig_costheta->Fill(costheta);
            mnvh1d_mdsig_theta->Fill(theta);
            mnvh1d_mdsig_pienergy->Fill(pienergy);
            mnvh1d_mdsig_pimom->Fill(pimom);
            mnvh1d_mdsig_kinetic->Fill(kinetic);
            mnvh1d_mdsig_mmom->Fill(mmom);
            mnvh1d_mdsig_pz->Fill(pz);
            mnvh1d_mdsig_pT->Fill(pT);
            mnvh1d_mdsig_mutheta->Fill(mutheta);
            mnvh1d_mdsig_enu->Fill(enu);
            mnvh1d_mdsig_q2->Fill(q2);
            mnvh1d_mdsig_w->Fill(w);
            
            
            mnvh1d_gdsig_costheta->Fill(costheta0);
            mnvh1d_gdsig_theta->Fill(theta0);
            mnvh1d_gdsig_pienergy->Fill(pienergy0);
            mnvh1d_gdsig_pimom->Fill(pimom0);
            mnvh1d_gdsig_kinetic->Fill(kinetic0);
            mnvh1d_gdsig_mmom->Fill(mmom0);
            mnvh1d_gdsig_pz->Fill(pz0);
            mnvh1d_gdsig_pT->Fill(pT0);
            mnvh1d_gdsig_mutheta->Fill(mutheta0);
            mnvh1d_gdsig_enu->Fill(enu0);
            mnvh1d_gdsig_q2->Fill(q20);
            mnvh1d_gdsig_w->Fill(w0);

        }

        if (background) {
            
            mnvh1d_mdbkg_costheta->Fill(costheta);
            mnvh1d_mdbkg_theta->Fill(theta);
            mnvh1d_mdbkg_pienergy->Fill(pienergy);
            mnvh1d_mdbkg_pimom->Fill(pimom);
            mnvh1d_mdbkg_kinetic->Fill(kinetic);
            mnvh1d_mdbkg_mmom->Fill(mmom);
            mnvh1d_mdbkg_pz->Fill(pz);
            mnvh1d_mdbkg_pT->Fill(pT);
            mnvh1d_mdbkg_mutheta->Fill(mutheta);
            mnvh1d_mdbkg_enu->Fill(enu);
            mnvh1d_mdbkg_q2->Fill(q2);
            mnvh1d_mdbkg_w->Fill(w);


            mnvh1d_gdbkg_costheta->Fill(costheta0);
            mnvh1d_gdbkg_theta->Fill(theta0);
            mnvh1d_gdbkg_pienergy->Fill(pienergy0);
            mnvh1d_gdbkg_pimom->Fill(pimom0);
            mnvh1d_gdbkg_kinetic->Fill(kinetic0);
            mnvh1d_gdbkg_mmom->Fill(mmom0);
            mnvh1d_gdbkg_pz->Fill(pz0);
            mnvh1d_gdbkg_pT->Fill(pT0);
            mnvh1d_gdbkg_mutheta->Fill(mutheta0);
            mnvh1d_gdbkg_enu->Fill(enu0);
            mnvh1d_gdbkg_q2->Fill(q20);
            mnvh1d_gdbkg_w->Fill(w0);


        }
           
    }

    warping_file.Close();
    
    ofs.close();

    delete data_object; data_object = 0;
    
    std::cout << "Selected events: " << selectedEvent << std::endl;
    
        // Add[]ErrorBandAndFillCV should be done after the event loop, of course!
    add_error_bands_and_fill_with_cv<0>(mnvh1d_data_mgg,  __is_newflux);

        // costheta
    add_error_bands_and_fill_with_cv <1>(mnvh1d_data_costheta,    __is_newflux);
    add_error_bands_and_fill_with_cv <1>(mnvh1d_mdsig_costheta,   __is_newflux);
    add_error_bands_and_fill_with_cv <1>(mnvh1d_mdbkg_costheta,   __is_newflux);
    add_error_bands_and_fill_with_cv <1>(mnvh1d_gd_costheta,      __is_newflux);
    add_error_bands_and_fill_with_cv <1>(mnvh1d_gdsig_costheta,   __is_newflux);
    add_error_bands_and_fill_with_cv <1>(mnvh1d_gdbkg_costheta,   __is_newflux);

        // theta
    add_error_bands_and_fill_with_cv <5>(mnvh1d_data_theta,       __is_newflux);
    add_error_bands_and_fill_with_cv <5>(mnvh1d_mdsig_theta,      __is_newflux);
    add_error_bands_and_fill_with_cv <5>(mnvh1d_mdbkg_theta,      __is_newflux);
    add_error_bands_and_fill_with_cv <5>(mnvh1d_gd_theta,         __is_newflux);
    add_error_bands_and_fill_with_cv <5>(mnvh1d_gdsig_theta,      __is_newflux);
    add_error_bands_and_fill_with_cv <5>(mnvh1d_gdbkg_theta,      __is_newflux);

        // pienergy
    add_error_bands_and_fill_with_cv<10>(mnvh1d_data_pienergy,    __is_newflux);
    add_error_bands_and_fill_with_cv<10>(mnvh1d_mdsig_pienergy,   __is_newflux);
    add_error_bands_and_fill_with_cv<10>(mnvh1d_mdbkg_pienergy,   __is_newflux);
    add_error_bands_and_fill_with_cv<10>(mnvh1d_gd_pienergy,      __is_newflux);
    add_error_bands_and_fill_with_cv<10>(mnvh1d_gdsig_pienergy,   __is_newflux);
    add_error_bands_and_fill_with_cv<10>(mnvh1d_gdbkg_pienergy,   __is_newflux);

        // pimom
    add_error_bands_and_fill_with_cv<15>(mnvh1d_data_pimom,       __is_newflux);
    add_error_bands_and_fill_with_cv<15>(mnvh1d_mdsig_pimom,      __is_newflux);
    add_error_bands_and_fill_with_cv<15>(mnvh1d_mdbkg_pimom,      __is_newflux);
    add_error_bands_and_fill_with_cv<15>(mnvh1d_gd_pimom,         __is_newflux);
    add_error_bands_and_fill_with_cv<15>(mnvh1d_gdsig_pimom,      __is_newflux);
    add_error_bands_and_fill_with_cv<15>(mnvh1d_gdbkg_pimom,      __is_newflux);
    
        // kinetic
    add_error_bands_and_fill_with_cv<20>(mnvh1d_data_kinetic,     __is_newflux);
    add_error_bands_and_fill_with_cv<20>(mnvh1d_mdsig_kinetic,    __is_newflux);
    add_error_bands_and_fill_with_cv<20>(mnvh1d_mdbkg_kinetic,    __is_newflux);
    add_error_bands_and_fill_with_cv<20>(mnvh1d_gd_kinetic,       __is_newflux);
    add_error_bands_and_fill_with_cv<20>(mnvh1d_gdsig_kinetic,    __is_newflux);
    add_error_bands_and_fill_with_cv<20>(mnvh1d_gdbkg_kinetic,    __is_newflux);

        // mmom
    add_error_bands_and_fill_with_cv<25>(mnvh1d_data_mmom,        __is_newflux);
    add_error_bands_and_fill_with_cv<25>(mnvh1d_mdsig_mmom,       __is_newflux);
    add_error_bands_and_fill_with_cv<25>(mnvh1d_mdbkg_mmom,       __is_newflux);
    add_error_bands_and_fill_with_cv<25>(mnvh1d_gd_mmom,          __is_newflux);
    add_error_bands_and_fill_with_cv<25>(mnvh1d_gdsig_mmom,       __is_newflux);
    add_error_bands_and_fill_with_cv<25>(mnvh1d_gdbkg_mmom,       __is_newflux);


        // longitudinal mmom
    add_error_bands_and_fill_with_cv<50>(mnvh1d_data_pz,        __is_newflux);
    add_error_bands_and_fill_with_cv<50>(mnvh1d_mdsig_pz,       __is_newflux);
    add_error_bands_and_fill_with_cv<50>(mnvh1d_mdbkg_pz,       __is_newflux);
    add_error_bands_and_fill_with_cv<50>(mnvh1d_gd_pz,          __is_newflux);
    add_error_bands_and_fill_with_cv<50>(mnvh1d_gdsig_pz,       __is_newflux);
    add_error_bands_and_fill_with_cv<50>(mnvh1d_gdbkg_pz,       __is_newflux);


        // transverse mmom
    add_error_bands_and_fill_with_cv<55>(mnvh1d_data_pT,        __is_newflux);
    add_error_bands_and_fill_with_cv<55>(mnvh1d_mdsig_pT,       __is_newflux);
    add_error_bands_and_fill_with_cv<55>(mnvh1d_mdbkg_pT,       __is_newflux);
    add_error_bands_and_fill_with_cv<55>(mnvh1d_gd_pT,          __is_newflux);
    add_error_bands_and_fill_with_cv<55>(mnvh1d_gdsig_pT,       __is_newflux);
    add_error_bands_and_fill_with_cv<55>(mnvh1d_gdbkg_pT,       __is_newflux);

        // mutheta
    add_error_bands_and_fill_with_cv<30>(mnvh1d_data_mutheta,     __is_newflux);
    add_error_bands_and_fill_with_cv<30>(mnvh1d_mdsig_mutheta,    __is_newflux);
    add_error_bands_and_fill_with_cv<30>(mnvh1d_mdbkg_mutheta,    __is_newflux);
    add_error_bands_and_fill_with_cv<30>(mnvh1d_gd_mutheta,       __is_newflux);
    add_error_bands_and_fill_with_cv<30>(mnvh1d_gdsig_mutheta,    __is_newflux);
    add_error_bands_and_fill_with_cv<30>(mnvh1d_gdbkg_mutheta,    __is_newflux);

        // enu
    add_error_bands_and_fill_with_cv<35>(mnvh1d_data_enu,         __is_newflux);
    add_error_bands_and_fill_with_cv<35>(mnvh1d_mdsig_enu,        __is_newflux);
    add_error_bands_and_fill_with_cv<35>(mnvh1d_mdbkg_enu,        __is_newflux);
    add_error_bands_and_fill_with_cv<35>(mnvh1d_gd_enu,           __is_newflux);
    add_error_bands_and_fill_with_cv<35>(mnvh1d_gdsig_enu,        __is_newflux);
    add_error_bands_and_fill_with_cv<35>(mnvh1d_gdbkg_enu,        __is_newflux);

        // q2
    add_error_bands_and_fill_with_cv<40>(mnvh1d_data_q2,          __is_newflux);
    add_error_bands_and_fill_with_cv<40>(mnvh1d_mdsig_q2,         __is_newflux);
    add_error_bands_and_fill_with_cv<40>(mnvh1d_mdbkg_q2,         __is_newflux);
    add_error_bands_and_fill_with_cv<40>(mnvh1d_gd_q2,            __is_newflux);
    add_error_bands_and_fill_with_cv<40>(mnvh1d_gdsig_q2,         __is_newflux);
    add_error_bands_and_fill_with_cv<40>(mnvh1d_gdbkg_q2,         __is_newflux);

        // w
    add_error_bands_and_fill_with_cv<45>(mnvh1d_data_w,           __is_newflux);
    add_error_bands_and_fill_with_cv<45>(mnvh1d_mdsig_w,          __is_newflux);
    add_error_bands_and_fill_with_cv<45>(mnvh1d_mdbkg_w,          __is_newflux);
    add_error_bands_and_fill_with_cv<45>(mnvh1d_gd_w,             __is_newflux);
    add_error_bands_and_fill_with_cv<45>(mnvh1d_gdsig_w,          __is_newflux);
    add_error_bands_and_fill_with_cv<45>(mnvh1d_gdbkg_w,          __is_newflux);
    

    TFile f(__histogram_file.c_str(), "recreate");

    std::cout << "Write histograms to file..." << std::endl;

        // Save MC histograms
    mnvh1d_mc_mgg->Write();
    mnvh1d_mcsig_mgg->Write();
    mnvh1d_mcbkg_mgg->Write();

    mnvh1d_mc_mgg2->Write();
    mnvh1d_mcsig_mgg2->Write();
    mnvh1d_mcbkg_mgg2->Write();

        // pi0 costheta
    mnvh1d_mc_costheta->Write();
    mnvh1d_mcsig_costheta->Write();
    mnvh1d_mcbkg_costheta->Write();
    
    mnvh1d_gn_costheta->Write();
    mnvh1d_gnsig_costheta->Write();
    mnvh1d_gnbkg_costheta->Write();

    
    mnvh1d_data_costheta->Write();
    mnvh1d_mdsig_costheta->Write();
    mnvh1d_mdbkg_costheta->Write();
    
    mnvh1d_gd_costheta->Write();
    mnvh1d_gdsig_costheta->Write();
    mnvh1d_gdbkg_costheta->Write();

    
        // pi0 theta
    mnvh1d_mc_theta->Write();
    mnvh1d_mcsig_theta->Write();
    mnvh1d_mcbkg_theta->Write();

    mnvh1d_gn_theta->Write();
    mnvh1d_gnsig_theta->Write();
    mnvh1d_gnbkg_theta->Write();


    mnvh1d_data_theta->Write();
    mnvh1d_mdsig_theta->Write();
    mnvh1d_mdbkg_theta->Write();

    mnvh1d_gd_theta->Write();
    mnvh1d_gdsig_theta->Write();
    mnvh1d_gdbkg_theta->Write();

    
        // pi0 total energy
    mnvh1d_mc_pienergy->Write();
    mnvh1d_mcsig_pienergy->Write();
    mnvh1d_mcbkg_pienergy->Write();
    
    mnvh1d_gn_pienergy->Write();
    mnvh1d_gnsig_pienergy->Write();
    mnvh1d_gnbkg_pienergy->Write();


    mnvh1d_data_pienergy->Write();
    mnvh1d_mdsig_pienergy->Write();
    mnvh1d_mdbkg_pienergy->Write();
    
    mnvh1d_gd_pienergy->Write();
    mnvh1d_gdsig_pienergy->Write();
    mnvh1d_gdbkg_pienergy->Write();

    
        // pi0 momentum
    mnvh1d_mc_pimom->Write();
    mnvh1d_mcsig_pimom->Write();
    mnvh1d_mcbkg_pimom->Write();

    mnvh1d_gn_pimom->Write();
    mnvh1d_gnsig_pimom->Write();
    mnvh1d_gnbkg_pimom->Write();


    mnvh1d_data_pimom->Write();
    mnvh1d_mdsig_pimom->Write();
    mnvh1d_mdbkg_pimom->Write();

    mnvh1d_gd_pimom->Write();
    mnvh1d_gdsig_pimom->Write();
    mnvh1d_gdbkg_pimom->Write();

    
        // pi0 kinetic energy
    mnvh1d_mc_kinetic->Write();
    mnvh1d_mcsig_kinetic->Write();
    mnvh1d_mcbkg_kinetic->Write();

    mnvh1d_gn_kinetic->Write();
    mnvh1d_gnsig_kinetic->Write();
    mnvh1d_gnbkg_kinetic->Write();


    mnvh1d_data_kinetic->Write();
    mnvh1d_mdsig_kinetic->Write();
    mnvh1d_mdbkg_kinetic->Write();

    mnvh1d_gd_kinetic->Write();
    mnvh1d_gdsig_kinetic->Write();
    mnvh1d_gdbkg_kinetic->Write();

    
        // muon momentum
    mnvh1d_mc_mmom->Write();
    mnvh1d_mcsig_mmom->Write();
    mnvh1d_mcbkg_mmom->Write();

    mnvh1d_gn_mmom->Write();
    mnvh1d_gnsig_mmom->Write();
    mnvh1d_gnbkg_mmom->Write();


    mnvh1d_data_mmom->Write();
    mnvh1d_mdsig_mmom->Write();
    mnvh1d_mdbkg_mmom->Write();

    mnvh1d_gd_mmom->Write();
    mnvh1d_gdsig_mmom->Write();
    mnvh1d_gdbkg_mmom->Write();

        // muon longitudinal momentum
    mnvh1d_mc_pz->Write();
    mnvh1d_mcsig_pz->Write();
    mnvh1d_mcbkg_pz->Write();

    mnvh1d_gn_pz->Write();
    mnvh1d_gnsig_pz->Write();
    mnvh1d_gnbkg_pz->Write();


    mnvh1d_data_pz->Write();
    mnvh1d_mdsig_pz->Write();
    mnvh1d_mdbkg_pz->Write();

    mnvh1d_gd_pz->Write();
    mnvh1d_gdsig_pz->Write();
    mnvh1d_gdbkg_pz->Write();


        // muon transverse momentum
    mnvh1d_mc_pT->Write();
    mnvh1d_mcsig_pT->Write();
    mnvh1d_mcbkg_pT->Write();

    mnvh1d_gn_pT->Write();
    mnvh1d_gnsig_pT->Write();
    mnvh1d_gnbkg_pT->Write();


    mnvh1d_data_pT->Write();
    mnvh1d_mdsig_pT->Write();
    mnvh1d_mdbkg_pT->Write();

    mnvh1d_gd_pT->Write();
    mnvh1d_gdsig_pT->Write();
    mnvh1d_gdbkg_pT->Write();

    
        // muon angle
    mnvh1d_mc_mutheta->Write();
    mnvh1d_mcsig_mutheta->Write();
    mnvh1d_mcbkg_mutheta->Write();

    mnvh1d_gn_mutheta->Write();
    mnvh1d_gnsig_mutheta->Write();
    mnvh1d_gnbkg_mutheta->Write();


    mnvh1d_data_mutheta->Write();
    mnvh1d_mdsig_mutheta->Write();
    mnvh1d_mdbkg_mutheta->Write();

    mnvh1d_gd_mutheta->Write();
    mnvh1d_gdsig_mutheta->Write();
    mnvh1d_gdbkg_mutheta->Write();


        // Neutrino energy
    mnvh1d_mc_enu->Write();
    mnvh1d_mcsig_enu->Write();
    mnvh1d_mcbkg_enu->Write();

    mnvh1d_gn_enu->Write();
    mnvh1d_gnsig_enu->Write();
    mnvh1d_gnbkg_enu->Write();


    mnvh1d_data_enu->Write();
    mnvh1d_mdsig_enu->Write();
    mnvh1d_mdbkg_enu->Write();

    mnvh1d_gd_enu->Write();
    mnvh1d_gdsig_enu->Write();
    mnvh1d_gdbkg_enu->Write();

    
        // Q2
    mnvh1d_mc_q2->Write();
    mnvh1d_mcsig_q2->Write();
    mnvh1d_mcbkg_q2->Write();

    mnvh1d_gn_q2->Write();
    mnvh1d_gnsig_q2->Write();
    mnvh1d_gnbkg_q2->Write();


    mnvh1d_data_q2->Write();
    mnvh1d_mdsig_q2->Write();
    mnvh1d_mdbkg_q2->Write();

    mnvh1d_gd_q2->Write();
    mnvh1d_gdsig_q2->Write();
    mnvh1d_gdbkg_q2->Write();
    

        // Hadronic invariant mass
    mnvh1d_mc_w->Write();
    mnvh1d_mcsig_w->Write();
    mnvh1d_mcbkg_w->Write();

    mnvh1d_gn_w->Write();
    mnvh1d_gnsig_w->Write();
    mnvh1d_gnbkg_w->Write();


    mnvh1d_data_w->Write();
    mnvh1d_mdsig_w->Write();
    mnvh1d_mdbkg_w->Write();

    mnvh1d_gd_w->Write();
    mnvh1d_gdsig_w->Write();
    mnvh1d_gdbkg_w->Write();

    
        // Save data histograms
    mnvh1d_data_mgg->Write();


    TParameter<double> par_emin("__emin", __emin);
    TParameter<double> par_emax("__emax", __emax);
    TParameter<double> par_wmin("__wmin", __wmin);
    TParameter<double> par_wmax("__wmax", __wmax);

    TParameter<double> par_mfp("__mfp", __mfp);
    TParameter<double> par_sigma("__sigma", __sigma);
    TParameter<double> par_nsigma("__nsigma", __nsigma);

    par_emin.Write();
    par_emax.Write();
    par_wmin.Write();
    par_wmax.Write();

    par_mfp.Write();
    par_sigma.Write();
    par_nsigma.Write();


    
    f.Close();
    
    std::cout << "Event tally" << std::endl;    
    std::cout << "  cut     signal     background      mc     data " << std::endl;

    std::ofstream evt_ofs(Form("%s-%s-event-tally.txt",__detector.c_str(), __label.c_str()));
    
    for (std::map<int, int>::iterator c = mctot_counter.begin();
         c != mctot_counter.end(); ++c) {
        std::cout << setw(6) << c->first
                  << setw(10) << mcsig_counter[c->first] 
                  << setw(12) << mcbkg_counter[c->first]
                  << setw(12) << mctot_counter[c->first]
                  << setw(12) << data_counter[c->first]
                  <<  std::endl;
        
        evt_ofs << setw(6) << c->first << ","
                << setw(10) << mcsig_counter[c->first] << ","
                << setw(12) << mcbkg_counter[c->first] << ","
                << setw(12) << mctot_counter[c->first] << ","
                << setw(12) << data_counter[c->first]  << ","
                <<  "\n";
        
    }
    
    evt_ofs.close();

    return 0;
    
}

// it is not a bug: not implement new flux for mockup data since the test
// with the nominal flux should be enough
