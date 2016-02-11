#ifndef calc_flux_h
#define calc_flux_h

#include <PlotUtils/TargetUtils.h>

#include "CCProtonPi0_CrossSection.h"
#include "constants.h"
#include "integrate.h"
#include "piecewise_delta_flux.h"
#include "smooth_delta_flux.h"

using std::setw;

using PlotUtils::MnvH1D;
using PlotUtils::TargetUtils;

namespace {

        // this version has problem with the bin-edge effect, but
        // negligible for flux-integrated cross sections (not true for enu)
    double root_integrate(TH1D* flux_histogram, double e_min, double e_max)
    {
            // integrate flux CV
        int b_min = flux_histogram->FindBin(e_min);
        int b_max = flux_histogram->FindBin(e_max);
    
        double integrated_flux = flux_histogram->Integral( b_min, b_max, "width" );    
        
        return integrated_flux;
    }

}

MnvH1D* CCProtonPi0_CrossSection::calc_flux( MnvH1D* mnvh1d_template,          // Template histogram to copy the binning for the flux histogram 
                   const std::string& flux_filename, // The flux file name
                   bool enu_cut,                     // whether to apply the enu cut
                   bool isEv,                        // if this variable is Enu since it is treated differently
                   bool __reweight_flux ,     // whether to do the flux reweight study
                   double __reweight_emin,     // lower bound of the reweighted region 
                   double __reweight_emax,     // upper bound of the reweighted region
                   double __reweight_amount,   // amount to reweight
                   bool flux_smooth_curve,   // reweight the flux over full neutrino energy range
                   bool flux_piecewise_curve)
{
    
        // Effectively the whole neutrino spectrum
    double E_min = 0.0;
    double E_max = 20.0;

        // Apply the neutrino energy cut if 'enu_cut' is TRUE
    if (enu_cut) {
        E_min = constants().emin;  // GeV
        E_max = constants().emax;  // GeV
    }

    const int n_universe = constants().n_universe;
 
    std::cout << " Using flux file: " << flux_filename << std::endl;
    TFile* flux_file = new TFile(flux_filename.c_str(), "READ" );
    assert(flux_file);
    
    MnvH1D* mnvh1d_flux_original = dynamic_cast<MnvH1D*>(flux_file->Get( "flux_E_cvweighted"));
    if (!mnvh1d_flux_original) {
            // not an MnvH1D object: can happen if just try to calculate the true GENIE xsec
        std::cout << "Note that the flux is not an MnvH1D object" << std::endl;
        TH1D* tmpflux = static_cast<TH1D*>(flux_file->Get( "flux_E_cvweighted"));
        mnvh1d_flux_original = new MnvH1D(*tmpflux);
    }
    
    assert(mnvh1d_flux_original);
    std::cout << "Print the flux (0-15 GeV)" << std::endl;
    for (int i = 1; ; ++i) { break;
        double low_edge  = mnvh1d_flux_original->GetBinLowEdge(i);
        double high_edge = mnvh1d_flux_original->GetBinLowEdge(i) + mnvh1d_flux_original->GetBinWidth(i);
        double content = mnvh1d_flux_original->GetBinContent(i);
        std::cout << "\t" << setw(8)
                  << Form("%2.2f", low_edge)
                  << " - "
                  << Form("%2.2f", high_edge)
                  << setw(15) << Form("%.4e",content)
                  << std::endl;

        if (high_edge > 15.0) break;
        
    }
    
    if (mnvh1d_flux_original->GetXaxis()->GetXmin() > E_min ||
        mnvh1d_flux_original->GetXaxis()->GetXmax() < E_max) {
        
        std::cerr << "Flux range too small: ("
                  << mnvh1d_flux_original->GetXaxis()->GetXmin() << ","
                  << mnvh1d_flux_original->GetXaxis()->GetXmax() << ")"
                  << std::endl;
        exit(1);   
        
    }
        
        // make an integrated flux histogram with the same bins as the eff. corrected distribution
    MnvH1D* mnvh1d_flux_integrated = (MnvH1D*) mnvh1d_template->Clone( "flux_integrated" );
    mnvh1d_flux_integrated->Reset();
    mnvh1d_flux_integrated->ClearAllErrorBands();

        // change the flux if flux reweighting is being applied.
    if (__reweight_flux) {

            // reweight the central-value
        for (int i = 1;  i <= mnvh1d_flux_original->GetNbinsX(); ++i) {
            double enu0    = mnvh1d_flux_original->GetBinCenter(i);
            double content = mnvh1d_flux_original->GetBinContent(i);
            
            if (flux_smooth_curve) {
                
                double amount = smooth_delta_flux::get().calculate(enu0);
                content *= (1.0 + amount);
                
            } else if (flux_piecewise_curve) {
                
                double amount = piecewise_delta_flux::get().calculate(enu0);
                content *= (1.0 + amount);
                
            } else {
                
                if (__reweight_emin <= enu0 && enu0 <= __reweight_emax) content *= (1.0 + __reweight_amount);
                
            }
            
            mnvh1d_flux_original->SetBinContent(i, content);

        }

            // reweight universe histograms
        std::vector<std::string> error_names = mnvh1d_flux_original->GetVertErrorBandNames();
        for (std::vector<std::string>::iterator error_name = error_names.begin();
             error_name != error_names.end(); ++error_name) {

            for (int universe = 0; universe < n_universe; ++universe ){

                TH1D* h_flux_original_var = mnvh1d_flux_original->GetVertErrorBand(*error_name)->GetHist(universe);
                
                for( int i = 1; i <= h_flux_original_var->GetNbinsX(); ++i ){  
                    double enu0    = h_flux_original_var->GetBinCenter(i);
                    double content = h_flux_original_var->GetBinContent(i);
                    
                    if (flux_smooth_curve) {
                        
                        double amount = smooth_delta_flux::get().calculate(enu0);
                        content *= (1.0 + amount);
                        
                    } else if (flux_piecewise_curve) {
                        
                        double amount = piecewise_delta_flux::get().calculate(enu0);
                        content *= (1.0 + amount);
                        
                    } else {
                        
                        if (__reweight_emin <= enu0 && enu0 <= __reweight_emax) content *= (1.0 + __reweight_amount);
                        
                    }

                    
                    h_flux_original_var->SetBinContent(i, content);
                    
                }
                
            }
            
        }
        
    } // reweight_flux
    

    if ( isEv ){
            // Note that in general the flux is binned much finer than the data binning
            // Loop over these tiny bins for each data bin
        std::cout << "Enu variable, integrate the flux over each data bin" << std::endl;
            // integrate flux CV
        for( int i=1; i <= mnvh1d_flux_integrated->GetNbinsX(); ++i ){  
            double e_min = mnvh1d_template->GetBinLowEdge(i); 
            double e_max = mnvh1d_template->GetBinLowEdge(i+1); 
            if (e_max > E_max) break;
            int b_min = mnvh1d_flux_original->FindBin(e_min);
            int b_max = mnvh1d_flux_original->FindBin(e_max) - 1;
            double flux_cv_fixed = integrate(mnvh1d_flux_original, e_min, e_max);
            double width = mnvh1d_flux_original->GetBinWidth(1);             // flux histogram has uniform bins
            mnvh1d_flux_integrated->SetBinContent(i, flux_cv_fixed * width); // 'width' is the flux bin width, similar to
                                                                             // "width" option in TH1::Integral()
            
            std::cout << "\tthere are " << setw(4) << (b_max-b_min) << " flux bins for this data bin: "
                      << setw(4) << i << setw(8) << Form("%2.2f", e_min) << " - " << Form("%2.2f", e_max)
                    //<< setw(12) << Form("%.4e", flux_cv)
                      << setw(12) << Form("%.4e", flux_cv_fixed)
                      << std::endl;

    
        }
        
            // integrate flux universes 
        std::vector<std::string> error_names = mnvh1d_flux_original->GetVertErrorBandNames();
        for (std::vector<std::string>::iterator error_name=error_names.begin();
             error_name!=error_names.end(); ++error_name) {

            if (n_universe > (int) mnvh1d_flux_original->GetVertErrorBand(*error_name)->GetNHists()) {
                std::cout << "Not enough universes from flux histogram" << std::endl;
                exit(1);
            }

            std::vector<TH1D*> histos;
            for (int universe = 0; universe < n_universe; ++universe) {
                TH1D* histo = new TH1D( (TH1D)*mnvh1d_flux_integrated);
                histo->SetName(Form("tmp_universe_%i", universe));
                histo->Reset();
                for (int i=1; i < mnvh1d_flux_integrated->GetNbinsX()+2; ++i){  
                    double e_min = mnvh1d_template->GetBinLowEdge(i); 
                    double e_max = mnvh1d_template->GetBinLowEdge(i+1); 
                    TH1D* h_flux_original_var = mnvh1d_flux_original->GetVertErrorBand(*error_name)->GetHist(universe);
                    double flux_var_fixed = integrate(h_flux_original_var, e_min, e_max);
                    double width = mnvh1d_flux_original->GetBinWidth(1);
                    histo->SetBinContent( i, flux_var_fixed * width);
                }
                
                histos.push_back(histo);
            }
            
            mnvh1d_flux_integrated->AddVertErrorBand(*error_name, histos);
            
            for(std::vector<TH1D*>::iterator h = histos.begin();
                h!=histos.end(); ++h){
                delete (*h);
            }
            
        }
    }
    
    else {  // other variables: integrate over the whole enu range [1.5,20] (GeV)
            // integrate flux CV
        std::cout << "Non-Enu variable, integrate the flux over neutrino energy" << std::endl;
        std::cout << "Integrating range: " << E_min << " < Ev < " << E_max << std::endl;
        double integrated_flux_cv = root_integrate(mnvh1d_flux_original, E_min, E_max);
        for (int i = 1; i < mnvh1d_flux_integrated->GetNbinsX()+2; ++i) {  
            mnvh1d_flux_integrated->SetBinContent(i, integrated_flux_cv);
        }

            // integrate flux universes
            // create integrated flux for all universes
        std::vector<std::string> error_names = mnvh1d_flux_original->GetVertErrorBandNames();
        for (std::vector<std::string>::iterator error_name=error_names.begin();
             error_name!=error_names.end(); ++error_name){

            std::cout << "\t error_band: " << *error_name << std::endl;
            if (n_universe > (int)mnvh1d_flux_original->GetVertErrorBand(*error_name)->GetNHists()) {
                std::cout << "Not enough universes from flux histogram" << std::endl;
                exit(1);
            }

            std::vector<TH1D*> histos;
            for (int universe = 0; universe < n_universe; ++universe) {
                TH1D* histo = new TH1D( (TH1D)*mnvh1d_flux_integrated );
                histo->SetName( Form( "tmp_universe_%i", universe ) );
                histo->Reset();
                TH1D* h_flux_original_var = mnvh1d_flux_original->GetVertErrorBand(*error_name)->GetHist(universe);
                double integrated_flux_var = root_integrate(h_flux_original_var, E_min, E_max);
                    //std::cout << "\tuniverse " << universe << " integrated flux " << integrated_flux_var << std::endl;
                for (int i=1; i < mnvh1d_flux_integrated->GetNbinsX()+2; ++i){  
                    histo->SetBinContent( i, integrated_flux_var );
                }

                histos.push_back( histo );
            }

            mnvh1d_flux_integrated->AddVertErrorBand(*error_name, histos);
            
            for (std::vector<TH1D*>::iterator h = histos.begin();
                h!=histos.end(); ++h) {
                delete (*h);
            }
    
        }
    }
        // original flux is no needed anymore
    delete mnvh1d_flux_original;
    
    mnvh1d_flux_integrated->AddMissingErrorBandsAndFillWithCV(*mnvh1d_template);
    
    
        // integrated flux units are /m^2/POT, convert it to /cm^2/POT 
    mnvh1d_flux_integrated->Scale(1.0e-4 );
    std::cout << "Integrated flux (10^-8 v/cm^2/POT)" << std::endl;
    std::cout.precision(5);
    for (int i = 1; i <= mnvh1d_flux_integrated->GetNbinsX(); ++i) {
        //std::cout << setw(8) << i << setw(20) << mnvh1d_flux_integrated->GetBinContent(i)/1e-8
                  //<< std::endl;
    }

    
    return mnvh1d_flux_integrated;

}


#endif
