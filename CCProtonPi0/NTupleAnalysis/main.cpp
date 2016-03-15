/*
================================================================================
main.cpp
    Main Function that controls the Analysis Package
    
    Classes Used:
        CCProtonPi0_Analyzer core class for the package
        CCProtonPi0_Plotter class includes functions specific to generate Plots
    
    Input and Output folders assigned in this function
    Edit isTest variable to run a Test job or all playlist
    
    Build Package using "make"
    
    Usage:
    Reduce NTuple to a Single File
            > ./main.exe reduce
    Analyze Reduced File
            > ./main.exe run
    Plot Analyzed File
            > ./main.exe plot

    Author:        Ozgur Altinok  - ozgur.altinok@tufts.edu
================================================================================
*/

// Include Required Classes
#include "Classes/Analyzer/CCProtonPi0_Analyzer.h"
#include "Classes/CrossSection/CCProtonPi0_CrossSection.h"
#include "Classes/SideBandTool/CCProtonPi0_SideBandTool.h"
#include "Classes/Plotter/CCProtonPi0_Plotter.h"
#include "Cintex/Cintex.h"

#include <string>
#include <ctime>

using namespace std;

const string runOption_Run = "run";
const string runOption_Plot = "plot";
const string runOption_Reduce = "reduce";
const string runOption_CrossSection = "calc";
const string runOption_FitSideBand = "fit";

const string typeOption_mc = "mc";
const string typeOption_data = "data";

int GetMode(int argc, char* argv[]);
void showInputError(char *argv[]);
void Plot();
void FitSideBands();
void Reduce(string playlist, bool isMC);
void Analyze(string playlist, bool isMC);
void Calculate_CrossSection(bool isMC);
void FitMinuit();
double calc_ChiSq_SideBand(SideBand &sb, Double_t *par, bool isPartial = false, int min_bin = 1, int max_bin = 1);
void calc_ChiSq(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

CCProtonPi0_SideBandTool sbtool;
int main(int argc, char *argv[] )
{
    time_t timeStart; time(&timeStart);
    time_t timeEnd;
    double timeDiff;
    int timeDiff_m;
    int timeDiff_s;
    string pl_reduce;
    string pl_analyze;
    bool isMC;

    TH1::AddDirectory(false);
    
    // Check User Command
    int nMode = GetMode(argc,argv);
    if (nMode == 0){
        showInputError(argv);
        return 0;
    }

    if ( nMode < 0) isMC = true;
    else isMC = false;

    if (nMode != 10 || nMode != 20){
        if (isMC){
            cout<<"MC Playlists Selected!\n"<<endl;
            pl_reduce = "Input/Playlists/pl_MC_Merged.dat"; 
            pl_analyze = "Input/Playlists/pl_MC_Reduced.dat"; 
        }else{
            cout<<"Data Playlists Selected!\n"<<endl;
            pl_reduce = "Input/Playlists/pl_Data_Merged.dat"; 
            pl_analyze = "Input/Playlists/pl_Data_Reduced.dat"; 
        }
    }

    ROOT::Cintex::Cintex::Enable();

    if ( abs(nMode) == 1) Reduce(pl_reduce, isMC);
    else if ( abs(nMode) == 2) Analyze(pl_analyze, isMC);
    else if ( abs(nMode) == 3) Calculate_CrossSection(isMC);
    else if ( nMode == 10) Plot();
    else if ( nMode == 20) FitSideBands();
    else{
        cout<<"Problem on Mode!, Returning"<<endl;
        return 0;
    }

    time(&timeEnd);
    timeDiff = ( timeEnd - timeStart );
    
    timeDiff_m = (int)timeDiff / 60;
    timeDiff_s = (int)timeDiff % 60;
    
    cout<<"Time to Complete = "<<timeDiff_m<<"\' "<<timeDiff_s<<"\" or "<<timeDiff<<" seconds"<<endl;
    return 0;
}

void Reduce(string playlist, bool isMC)
{
    bool isModeReduce = true;
    cout<<"\n"<<endl;
    cout<<"======================================================================"<<endl;
    cout<<"Reducing NTuples..."<<endl;
    cout<<"======================================================================"<<endl;
    CCProtonPi0_Analyzer t(isModeReduce, isMC);
    t.reduce(playlist);
}

void Analyze(string playlist, bool isMC)
{
    bool isModeReduce = false;
    cout<<"\n"<<endl;
    cout<<"======================================================================"<<endl;
    cout<<"Analyzing NTuples, Creating Histograms..."<<endl;
    cout<<"======================================================================"<<endl;
    CCProtonPi0_Analyzer analyzer(isModeReduce, isMC);
    analyzer.analyze(playlist);
}

void Calculate_CrossSection(bool isMC)
{
    cout<<"\n"<<endl;
    cout<<"======================================================================"<<endl;
    cout<<"Calculating Cross Section..."<<endl;
    cout<<"======================================================================"<<endl;
    CCProtonPi0_CrossSection crossSection(isMC);
    crossSection.Calc_CrossSections();
}


void Plot()
{
    cout<<"======================================================================"<<endl;
    cout<<"Plotting Histograms..."<<endl;
    cout<<"======================================================================"<<endl;
    CCProtonPi0_Plotter plotter;
    plotter.plotHistograms();
}

void FitSideBands()
{
    cout<<"======================================================================"<<endl;
    cout<<"Fitting Side Bands..."<<endl;
    cout<<"======================================================================"<<endl;
    FitMinuit();

}

/*
 *  1   reduce
 *  2   run
 *  3   calculate cross section
 *  10  plot
 *  20  fit side band
 *   
*/
int GetMode(int argc, char* argv[])
{
    // argc can only be 2 or 3
    if (argc != 2 && argc != 3) return 0;

    std::string runSelect = argv[1];
    if (argc == 2){
        if (runSelect.compare(runOption_Plot) == 0) return 10;
        else if (runSelect.compare(runOption_FitSideBand) == 0) return 20;
        else return 0;
    }
     
    std::string typeSelect = argv[2];
    // First check for ERROR
    if (runSelect.compare(runOption_Reduce) != 0 && runSelect.compare(runOption_Run) != 0 && runSelect.compare(runOption_CrossSection) != 0) return 0;
    if (typeSelect.compare(typeOption_mc) != 0 && typeSelect.compare(typeOption_data) != 0) return 0;

    // Passed ERROR Check - Valid Input    
    if (runSelect.compare(runOption_Reduce) == 0){
        if (typeSelect.compare(typeOption_mc) == 0) return -1;
        else if (typeSelect.compare(typeOption_data) == 0) return 1;
        else return 0;
    }

    if (runSelect.compare(runOption_Run) == 0){
        if (typeSelect.compare(typeOption_mc) == 0) return -2;
        else if (typeSelect.compare(typeOption_data) == 0) return 2;
        else return 0;
    }

    if (runSelect.compare(runOption_CrossSection) == 0){
        if (typeSelect.compare(typeOption_mc) == 0) return -3;
        else if (typeSelect.compare(typeOption_data) == 0) return 3;
        else return 0;
    }

    return 0;
}

void showInputError(char *argv[])
{
    cout<<std::left;
    cout<<"----------------------------------------------------------------------"<<endl;
    cout<<"Not a valid syntax!"<<endl;
    cout<<"----------------------------------------------------------------------"<<endl;
    cout<<"Correct Syntax for NTuple Reduce"<<endl;
    cout<<"\t"<<argv[0]<<" "<<runOption_Reduce<<" "<<typeOption_mc<<endl;
    cout<<"\t"<<argv[0]<<" "<<runOption_Reduce<<" "<<typeOption_data<<"\n"<<endl;
    cout<<"Correct Syntax for NTuple Analysis"<<endl;
    cout<<"\t"<<argv[0]<<" "<<runOption_Run<<" "<<typeOption_mc<<endl;
    cout<<"\t"<<argv[0]<<" "<<runOption_Run<<" "<<typeOption_data<<"\n"<<endl;
    cout<<"Correct Syntax for Calculating Cross Section"<<endl;
    cout<<"\t"<<argv[0]<<" "<<runOption_CrossSection<<" "<<typeOption_mc<<endl;
    cout<<"\t"<<argv[0]<<" "<<runOption_CrossSection<<" "<<typeOption_data<<"\n"<<endl;
    cout<<"Correct Syntax for Plotting"<<endl;
    cout<<"\t"<<argv[0]<<" "<<runOption_Plot<<"\n"<<endl;
    cout<<"Correct Syntax for Fitting SideBands"<<endl;
    cout<<"\t"<<argv[0]<<" "<<runOption_FitSideBand<<"\n"<<endl;
    cout<<"----------------------------------------------------------------------"<<endl;
}

double calc_ChiSq_SideBand(SideBand &sb, Double_t *par, bool isPartial, int min_bin, int max_bin)
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
    double ChiSq = 0;
    
    // Calculate ChiSq for pID for ALL Bins
    ChiSq += calc_ChiSq_SideBand(sbtool.pID, par);

    // Calculate ChiSq for Low Inv Mass Region for first 3 bins
    ChiSq += calc_ChiSq_SideBand(sbtool.LowInvMass, par, true, 1, 3);
    
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



