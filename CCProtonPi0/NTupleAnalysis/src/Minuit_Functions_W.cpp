#include <iostream>
#include <cstdlib>
#include <cmath>
#include <string>
#include <TFile.h>
#include <TMinuit.h>
#include <PlotUtils/MnvH1D.h>
#include "../Libraries/Folder_List.h"

using namespace PlotUtils;

void NormalizeHistogram(MnvH1D* h)
{
    int NBins = h->GetNbinsX();
    double area = h->Integral();
    double nOverFlow = h->GetBinContent(NBins+1);
    double nUnderFlow = h->GetBinContent(0);
    h->Scale(1/(area+nOverFlow+nUnderFlow),"",false); // Scale only on CentralValue
}

void calc_ChiSq_W(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    // Silence Unused Variable Warning
    (void) npar;
    (void) gin;
    (void) iflag;
    
    double data_POT = 3.33153e+20;
    double mc_POT = 2.21867e+21; 
    double POT_ratio = data_POT/mc_POT;

    double ChiSq = 0.0;

    // Get Background Subtracted Data
    TFile* f_data = new TFile(Folder_List::rootDir_Interaction_data.c_str());
    TFile* f_mc = new TFile(Folder_List::rootDir_Interaction_mc.c_str());
    MnvH1D* data = new MnvH1D( * dynamic_cast<MnvH1D*>(f_data->Get("W_All_0")) );
    MnvH1D* mc_bckg = new MnvH1D( * dynamic_cast<MnvH1D*>(f_mc->Get("W_All_2")) );

    NormalizeHistogram(mc_bckg);
    
    const double nBckg = 2997.1; // Got number from Cross Section Calculation Log File
    mc_bckg->Scale(nBckg);

    data->Add(mc_bckg, -1); 
   
    // Get MC Signal Types
    MnvH1D* mc_DeltaRES = new MnvH1D( * dynamic_cast<MnvH1D*>(f_mc->Get("W_All_7")) );
    MnvH1D* mc_OtherRES = new MnvH1D( * dynamic_cast<MnvH1D*>(f_mc->Get("W_All_8")) );
    MnvH1D* mc_NonRES = new MnvH1D( * dynamic_cast<MnvH1D*>(f_mc->Get("W_All_9")) );

    int nBins = data->GetNbinsX(); 

    for (int i = 1; i <= nBins ; ++i) {
        double nData = data->GetBinContent(i);
        if (nData == 0) continue;

        // par[] will be the weights associated with that background
        double nDeltaRES = par[0] * mc_DeltaRES->GetBinContent(i) * POT_ratio;
        double nOtherRES = par[1] * mc_OtherRES->GetBinContent(i) * POT_ratio;
        double nNonRES = par[2] * mc_NonRES->GetBinContent(i) * POT_ratio;
        
        double nTotalMC = nDeltaRES + nOtherRES + nNonRES;

        double delta  = std::pow((nData - nTotalMC),2)/nData;
        ChiSq += delta;
    }

    delete data;
    delete mc_bckg;
    delete mc_DeltaRES; 
    delete mc_OtherRES; 
    delete mc_NonRES; 
    delete f_data;
    delete f_mc;

    f = ChiSq;
    return;
}

void FitMinuit_W()
{
    TMinuit *ptMinuit = new TMinuit(3);  //initialize TMinuit with a maximum of 3 params
    //
    //  select verbose level:
    //    default :     (58 lines in this test)
    //    -1 : minimum  ( 4 lines in this test)
    //     0 : low      (31 lines)
    //     1 : medium   (61 lines)
    //     2 : high     (89 lines)
    //     3 : maximum (199 lines in this test)
    //
    ptMinuit->SetPrintLevel();
    // set the user function that calculates chi_square (the value to minimize)
    ptMinuit->SetFCN(calc_ChiSq_W);

    Double_t arglist[10];
    Int_t ierflg = 0;

    arglist[0] = 1;
    ptMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

    // Set starting values and step sizes for parameters
    //      Start from 1.0 and step 0.01
    ptMinuit->mnparm(0, "DeltaRES", 1.0, 0.01, 0.5, 2.0, ierflg);
    ptMinuit->mnparm(1, "OtherRES", 1.0, 0.01, 0.01, 2.0, ierflg);
    ptMinuit->mnparm(2, "NonRES", 1.0, 0.01, 0.5, 2.0, ierflg);

    // Now ready for minimization step
    arglist[0] = 500;
    arglist[1] = 1.;
    ptMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

    // Print results
    std::cout<<"\nPrint results from minuit\n";
    double fParamVal[3];
    double fParamErr[3];
    ptMinuit->GetParameter(0, fParamVal[0], fParamErr[0]);
    std::cout << "DeltaRES = " << fParamVal[0] << "\n";
    ptMinuit->GetParameter(1, fParamVal[1], fParamErr[1]);
    std::cout << "OtherRES = " << fParamVal[1] << "\n";
    ptMinuit->GetParameter(2, fParamVal[2], fParamErr[2]);
    std::cout << "NonRES = " << fParamVal[2] << "\n";
    
    // if you want to access to these parameters, use:
    Double_t min_ChiSq,edm,errdef;
    Int_t nvpar,nparx,icstat;
    ptMinuit->mnstat(min_ChiSq,edm,errdef,nvpar,nparx,icstat);
    //void mnstat(Double_t &fmin, Double_t &fedm, Double_t &errdef, Int_t &npari, Int_t &nparx, Int_t &istat) 
    //*-*-*-*-*Returns concerning the current status of the minimization*-*-*-*-*
    //*-*      =========================================================
    //*-*       User-called
    //*-*          Namely, it returns:
    //*-*        FMIN: the best function value found so far
    //*-*        FEDM: the estimated vertical distance remaining to minimum
    //*-*        ERRDEF: the value of UP defining parameter uncertainties
    //*-*        NPARI: the number of currently variable parameters
    //*-*        NPARX: the highest (external) parameter number defined by user
    //*-*        ISTAT: a status integer indicating how good is the covariance
    //*-*           matrix:  0= not calculated at all
    //*-*                    1= approximation only, not accurate
    //*-*                    2= full matrix, but forced positive-definite
    //*-*                    3= full accurate covariance matrix
    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
    
    
    std::cout << "\n";
    std::cout << " Minimum chi square = " << min_ChiSq<< "\n";
    std::cout << " Estimated vert. distance to min. = " << edm << "\n";
    std::cout << " Number of variable parameters = " << nvpar << "\n";
    std::cout << " Highest number of parameters defined by user = " << nparx << "\n";
    std::cout << " Status of covariance matrix = " << icstat << "\n";

    std::cout << "\n";
    ptMinuit->mnprin(1,min_ChiSq);
    //*-*-*-*Prints the values of the parameters at the time of the call*-*-*-*-*
    //*-*    ===========================================================
    //*-*        also prints other relevant information such as function value,
    //*-*        estimated distance to minimum, parameter errors, step sizes.
    //*-*
    //*-*         According to the value of IKODE, the printout is:
    //*-*    IKODE=INKODE= 0    only info about function value
    //*-*                  1    parameter values, errors, limits
    //*-*                  2    values, errors, step sizes, internal values
    //*-*                  3    values, errors, step sizes, first derivs.
    //*-*                  4    values, parabolic errors, MINOS errors
    //*-*    when INKODE=5, MNPRIN chooses IKODE=1,2, or 3, according to ISW(2)
    //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

    delete ptMinuit;
}

