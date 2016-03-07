#include "CCProtonPi0_SideBandTool.h"

using namespace PlotUtils;

void CCProtonPi0_SideBandTool::Fit()
{
    Fit(Michel);
    Fit(pID);
    Fit(LowInvMass, true, 1, 3);
}

void CCProtonPi0_SideBandTool::Fit(SideBand &sb, bool isLimitedFit, int first_bin, int last_bin)
{
    std::cout<<"Fitting Side Band: "<<sb.name<<std::endl;

    // Crate TObjArray
    TObjArray* mc_models = new TObjArray(5);
    mc_models->Add(sb.signal);
    mc_models->Add(sb.WithPi0);
    mc_models->Add(sb.QELike);
    mc_models->Add(sb.SinglePiPlus);
    mc_models->Add(sb.Other);

    // Fit data and mc_models
    TFractionFitter* fitter = new TFractionFitter(sb.data, mc_models, "Q");
    // Constrain Fractions [0,1]
    fitter->Constrain(0, 0.0, 1.0);
    fitter->Constrain(1, 0.0, 1.0);
    fitter->Constrain(2, 0.0, 1.0);
    fitter->Constrain(3, 0.0, 1.0);
    fitter->Constrain(4, 0.0, 1.0);

    if (isLimitedFit){
        fitter->SetRangeX(first_bin, last_bin);
    }

    Int_t status = fitter->Fit();

    if (status != 0) {
        std::cout<<"\tFit Error!"<<std::endl;
        return;
    }

    // Fit Results
    double total = 0;
    double total_err = 0;
    Double_t f;                         // Fraction
    Double_t f_err;                     // Fraction Error
    fitter->GetResult(0, f, f_err);     // access fitted value and uncertainty for parameter 1
    sb.fr_signal[1] = f;
    sb.fr_signal[2] = f_err;
    total += f;
    total_err += f_err;

    fitter->GetResult(1, f, f_err);
    sb.fr_WithPi0[1] = f;
    sb.fr_WithPi0[2] = f_err; 
    total += f;
    total_err += f_err;
   
    fitter->GetResult(2, f, f_err); 
    sb.fr_QELike[1] = f;
    sb.fr_QELike[2] = f_err;  
    total += f;
    total_err += f_err;
    
    fitter->GetResult(3, f, f_err); 
    sb.fr_SinglePiPlus[1] = f;
    sb.fr_SinglePiPlus[2] = f_err;  
    total += f;
    total_err += f_err;
    
    fitter->GetResult(4, f, f_err);
    sb.fr_Other[1] = f;
    sb.fr_Other[2] = f_err;  
    total += f;
    total_err += f_err;
  
    sb.fr_Total[1] = total;
    sb.fr_Total[2] = total_err;

    PrintFitResults(sb);

    std::cout<<"Done!"<<std::endl;
}

CCProtonPi0_SideBandTool::CCProtonPi0_SideBandTool() : CCProtonPi0_NTupleAnalysis()
{
    std::cout<<"Initializing CCProtonPi0_SideBandTool"<<std::endl;

    OpenTextFile();
    OpenRootFiles();
    initSideBands();

    std::cout<<"Done!"<<std::endl;
    std::cout<<std::endl;
}

void CCProtonPi0_SideBandTool::OpenRootFiles()
{
    std::string rootDir;

    rootDir = Folder_List::rootDir_sideBand_Michel_mc;
    Michel.f_mc = new TFile(rootDir.c_str());

    rootDir = Folder_List::rootDir_sideBand_Michel_data;
    Michel.f_data = new TFile(rootDir.c_str());

    rootDir = Folder_List::rootDir_sideBand_pID_mc;
    pID.f_mc = new TFile(rootDir.c_str());

    rootDir = Folder_List::rootDir_sideBand_pID_data;
    pID.f_data = new TFile(rootDir.c_str());

    rootDir = Folder_List::rootDir_sideBand_LowInvMass_mc;
    LowInvMass.f_mc = new TFile(rootDir.c_str());

    rootDir = Folder_List::rootDir_sideBand_LowInvMass_data;
    LowInvMass.f_data = new TFile(rootDir.c_str());
}

void CCProtonPi0_SideBandTool::initSideBands()
{
    SetNames(Michel, "Michel");
    GetTH1D(Michel.f_data, Michel.data, "hCut_pi0invMass_0");   
    GetTH1D(Michel.f_mc, Michel.signal, "hCut_pi0invMass_1");   
    GetTH1D(Michel.f_mc, Michel.WithPi0, "hCut_pi0invMass_3");   
    GetTH1D(Michel.f_mc, Michel.QELike, "hCut_pi0invMass_4");   
    GetTH1D(Michel.f_mc, Michel.SinglePiPlus, "hCut_pi0invMass_5");   
    GetTH1D(Michel.f_mc, Michel.Other, "hCut_pi0invMass_6");   
    CalcMCRatios(Michel);

    SetNames(pID, "pID");
    GetTH1D(pID.f_data, pID.data, "hCut_pi0invMass_0");   
    GetTH1D(pID.f_mc, pID.signal, "hCut_pi0invMass_1");   
    GetTH1D(pID.f_mc, pID.WithPi0, "hCut_pi0invMass_3");   
    GetTH1D(pID.f_mc, pID.QELike, "hCut_pi0invMass_4");   
    GetTH1D(pID.f_mc, pID.SinglePiPlus, "hCut_pi0invMass_5");   
    GetTH1D(pID.f_mc, pID.Other, "hCut_pi0invMass_6");   
    CalcMCRatios(pID);

    SetNames(LowInvMass, "LowInvMass");
    GetTH1D(LowInvMass.f_data, LowInvMass.data, "hCut_pi0invMass_0");   
    GetTH1D(LowInvMass.f_mc, LowInvMass.signal, "hCut_pi0invMass_1");   
    GetTH1D(LowInvMass.f_mc, LowInvMass.WithPi0, "hCut_pi0invMass_3");   
    GetTH1D(LowInvMass.f_mc, LowInvMass.QELike, "hCut_pi0invMass_4");   
    GetTH1D(LowInvMass.f_mc, LowInvMass.SinglePiPlus, "hCut_pi0invMass_5");   
    GetTH1D(LowInvMass.f_mc, LowInvMass.Other, "hCut_pi0invMass_6");   
    CalcMCRatios(LowInvMass, true, 1, 3);
}

void CCProtonPi0_SideBandTool::GetTH1D(TFile* f, TH1D* &h, std::string var_name)
{
    MnvH1D* temp;
    temp = new MnvH1D(*(MnvH1D*)f->Get(var_name.c_str()));
    h = new TH1D(temp->GetCVHistoWithStatError()); 
    delete temp;
}

void CCProtonPi0_SideBandTool::SetNames(SideBand &sb, std::string name)
{
    sb.name = name;
    sb.model_names[0] = "Signal";
    sb.model_names[1] = "WithPi0";
    sb.model_names[2] = "QELike";
    sb.model_names[3] = "SinglePiPlus";
    sb.model_names[4] = "Other";
    sb.model_names[5] = "Total";
}

void CCProtonPi0_SideBandTool::CalcMCRatios(SideBand &sb, bool isLimited, int first_bin, int last_bin)
{
    // If the Integral is NOT Limited then integrate over ALL bins
    if (!isLimited){
        first_bin = 1;
        last_bin = sb.signal->GetNbinsX();
    }

    double signal = sb.signal->Integral(first_bin, last_bin);
    double WithPi0 = sb.WithPi0->Integral(first_bin, last_bin);
    double QELike = sb.QELike->Integral(first_bin, last_bin);
    double SinglePiPlus = sb.SinglePiPlus->Integral(first_bin, last_bin);
    double Other = sb.Other->Integral(first_bin, last_bin);
    
    double total = signal + WithPi0 + QELike + SinglePiPlus + Other;

    sb.fr_signal[0] = signal / total;
    sb.fr_WithPi0[0] = WithPi0 / total;
    sb.fr_QELike[0] = QELike / total;
    sb.fr_SinglePiPlus[0] = SinglePiPlus / total;
    sb.fr_Other[0] = Other / total;
    sb.fr_Total[0] = total / total;
}

void CCProtonPi0_SideBandTool::PrintFitResults(SideBand &sb)
{
    using namespace std;
    textFile<<left;
    textFile<<setw(16)<<"Type"; 
    textFile<<setw(16)<<"MC Ratio"; 
    textFile<<setw(16)<<"Fit Ratio"; 
    textFile<<setw(16)<<"Error"; 
    textFile<<endl; 

    textFile<<setw(16)<<sb.model_names[0]; PrintRatio(sb.fr_signal);
    textFile<<setw(16)<<sb.model_names[1]; PrintRatio(sb.fr_WithPi0);
    textFile<<setw(16)<<sb.model_names[2]; PrintRatio(sb.fr_QELike);
    textFile<<setw(16)<<sb.model_names[3]; PrintRatio(sb.fr_SinglePiPlus);
    textFile<<setw(16)<<sb.model_names[4]; PrintRatio(sb.fr_Other);
    textFile<<setw(16)<<sb.model_names[5]; PrintRatio(sb.fr_Total);
    
    textFile<<endl;
}

void CCProtonPi0_SideBandTool::PrintRatio(double ratio[])
{
    using namespace std;
    textFile<<left;
    textFile<<setw(16)<<setprecision(3)<<ratio[0];
    textFile<<setw(16)<<setprecision(3)<<ratio[1];
    textFile<<setw(16)<<setprecision(3)<<ratio[2];
    textFile<<endl;
}

void CCProtonPi0_SideBandTool::OpenTextFile()
{
    fileName = Folder_List::output + Folder_List::textOut + "SideBandTables.txt";
    CCProtonPi0_NTupleAnalysis::OpenTextFile(fileName, textFile);
}

CCProtonPi0_SideBandTool::~CCProtonPi0_SideBandTool()
{
    textFile.close();
}

