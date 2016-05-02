#include "../Classes/SideBandTool/CCProtonPi0_SideBandTool.h"

using namespace std;

CCProtonPi0_SideBandTool sbtool;

double calc_ChiSq_SideBand(SideBand &sb, Double_t *par, bool isPartial = false, int min_bin = 1, int max_bin = 1)
{
    if (!isPartial){
        min_bin = 1;
        max_bin = sb.data->GetNbinsX();
    }
    
    if (min_bin == max_bin){
        cout<<"Wrong Range for Fit"<<endl;
        exit(EXIT_FAILURE);
    }

    double ChiSq = 0.0;

    for (int i = 1; i <= max_bin; ++i) {
        double nData = sb.data->GetBinContent(i);
        
        // Do not use Signal and Other in Fit
        double nSignal = sb.signal[0]->GetBinContent(i) * sbtool.POT_ratio;
        double nOther = sb.Other[0]->GetBinContent(i) * sbtool.POT_ratio;

        // par[] will be the weights associated with that background
        double nWithPi0 = par[0] * sb.WithPi0[0]->GetBinContent(i) * sbtool.POT_ratio;
        double nQELike = par[1] * sb.QELike[0]->GetBinContent(i) * sbtool.POT_ratio;
        double nSinglePiPlus = par[2] * sb.SinglePiPlus[0]->GetBinContent(i) * sbtool.POT_ratio;
        
        double nTotalMC = nSignal + nWithPi0 + nQELike + nSinglePiPlus + nOther;

        double delta  = std::pow((nData - nTotalMC),2)/nData;
        ChiSq += delta;
    }

    return ChiSq;
}

void calc_ChiSq(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    // Silence Unused Variable Warning
    (void) npar;
    (void) gin;
    (void) iflag;

    double ChiSq = 0;

    // Calculate ChiSq for Michel for ALL Bins
    ChiSq += calc_ChiSq_SideBand(sbtool.Michel, par);
   
    // Calculate ChiSq for pID for ALL Bins
    ChiSq += calc_ChiSq_SideBand(sbtool.pID, par);

    // Calculate ChiSq for Low Inv Mass Region for first 6 bins
    ChiSq += calc_ChiSq_SideBand(sbtool.LowInvMass, par, true, 1, 6);
 
    // Calculate ChiSq for High Inv Mass Region for last 30 bins
    ChiSq += calc_ChiSq_SideBand(sbtool.HighInvMass, par, true, 21, 50);
    
    f = ChiSq;
    return;
}

void FitMinuit()
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
    ptMinuit->SetFCN(calc_ChiSq);

    Double_t arglist[10];
    Int_t ierflg = 0;

    arglist[0] = 1;
    ptMinuit->mnexcm("SET ERR", arglist ,1,ierflg);

    // Set starting values and step sizes for parameters
    //      Start from 1.0 and step 0.01
    ptMinuit->mnparm(0, "WithPi0", 1.0, 0.01, 0.5, 2.0 , ierflg);
    ptMinuit->mnparm(1, "QELike", 1.0, 0.01, 0.5, 2.0, ierflg);
    ptMinuit->mnparm(2, "SinglePiPlus", 1.0, 0.01, 0.5, 2.0, ierflg);

    // Now ready for minimization step
    arglist[0] = 500;
    arglist[1] = 1.;
    ptMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

    // Print results
    std::cout<<"\nPrint results from minuit\n";
    double fParamVal[3];
    double fParamErr[3];
    ptMinuit->GetParameter(0, fParamVal[0], fParamErr[0]);
    std::cout << "WithPi0 = " << fParamVal[0] << "\n";
    ptMinuit->GetParameter(1, fParamVal[1], fParamErr[1]);
    std::cout << "QELike = " << fParamVal[1] << "\n";
    ptMinuit->GetParameter(2, fParamVal[2], fParamErr[2]);
    std::cout << "SinglePiPlus = " << fParamVal[2] << "\n";
    
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
    
    sbtool.ApplyFitResults(min_ChiSq, fParamVal[0],fParamVal[1],fParamVal[2]);
    sbtool.Plot();
    
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

}

void FitSideBands()
{
    cout<<"======================================================================"<<endl;
    cout<<"Fitting Side Bands..."<<endl;
    cout<<"======================================================================"<<endl;
    FitMinuit();
}


