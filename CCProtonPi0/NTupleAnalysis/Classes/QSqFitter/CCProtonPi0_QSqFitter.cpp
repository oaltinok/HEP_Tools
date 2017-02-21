#ifndef CCProtonPi0_QSqFitter_cpp
#define CCProtonPi0_QSqFitter_cpp

#include "CCProtonPi0_QSqFitter.h"

const double CCProtonPi0_QSqFitter::GENIE_MaRES = 1.12; // GeV
const double CCProtonPi0_QSqFitter::x1_1sigma = 0.0;
const double CCProtonPi0_QSqFitter::x2_1sigma = 100.0;
const double CCProtonPi0_QSqFitter::x1_2sigma = 100.0;
const double CCProtonPi0_QSqFitter::x2_2sigma = 200.0;

CCProtonPi0_QSqFitter::CCProtonPi0_QSqFitter()
{
    isDebug = false;
    FillMaRESVector(MaRESVector_dn, 0.8, 0.6);
    FillMaRESVector(MaRESVector_up, 1.2, 1.4);
}

void CCProtonPi0_QSqFitter::FillMaRESVector(std::vector<double> &MaRESVector, double one_sigma, double two_sigma)
{  
    double y1_1sigma = GENIE_MaRES; // CV 
    double y2_1sigma = GENIE_MaRES * one_sigma; // 20% Change 
    double y1_2sigma = GENIE_MaRES * one_sigma;
    double y2_2sigma = GENIE_MaRES * two_sigma; // 40% Change

    // Get Eq. of Line for 1 Sigma 
    double m_1sigma = Calc_Slope(x1_1sigma, x2_1sigma, y1_1sigma, y2_1sigma);
    double c_1sigma = Calc_Constant(m_1sigma, x1_1sigma, y1_1sigma); 

    // Fill Weights for 1 Sigma
    for (int i = 0; i <= 100; i++){
        double current_MaRES = Calc_Weight(m_1sigma, c_1sigma, i);
        MaRESVector.push_back(current_MaRES);
        if (isDebug){
            std::cout<<"\t"<<i<<" "<<current_MaRES<<std::endl;
        }
    }

    // Get Eq. of Line for 2 Sigma 
    double m_2sigma = Calc_Slope(x1_2sigma, x2_2sigma, y1_2sigma, y2_2sigma);
    double c_2sigma = Calc_Constant(m_2sigma, x1_2sigma, y1_2sigma); 

    for (int i = 101; i <= 200; i++){
        double current_MaRES = Calc_Weight(m_2sigma, c_2sigma, i);
        MaRESVector.push_back(current_MaRES);
        if (isDebug){
            std::cout<<"\t"<<i<<" "<<current_MaRES<<std::endl;
        }
    }
}

std::vector<double> CCProtonPi0_QSqFitter::GetWeights(double genie_1sigma, double genie_2sigma)
{
    std::vector<double> wgts;
   
    double y1_1sigma = 1.0; // CV Weight
    double y2_1sigma = genie_1sigma;
    double y1_2sigma = genie_1sigma;
    double y2_2sigma = genie_2sigma;

    if (isDebug){ 
        std::cout<<"Genie 1Sigma = "<<genie_1sigma<<std::endl;
        std::cout<<"Genie 2Sigma = "<<genie_2sigma<<std::endl;
        std::cout<<"Partitioned Weights"<<std::endl;
    }

    // Get Eq. of Line for 1 Sigma 
    double m_1sigma = Calc_Slope(x1_1sigma, x2_1sigma, y1_1sigma, y2_1sigma);
    double c_1sigma = Calc_Constant(m_1sigma, x1_1sigma, y1_1sigma); 

    // Fill Weights for 1 Sigma
    for (int i = 0; i <= 100; i++){
        double wgt = Calc_Weight(m_1sigma, c_1sigma, i);
        wgts.push_back(wgt);
        if (isDebug){
            std::cout<<"\t"<<i<<" "<<wgt<<std::endl;
        }
    }

    // Get Eq. of Line for 2 Sigma 
    double m_2sigma = Calc_Slope(x1_2sigma, x2_2sigma, y1_2sigma, y2_2sigma);
    double c_2sigma = Calc_Constant(m_2sigma, x1_2sigma, y1_2sigma); 

    for (int i = 101; i <= 200; i++){
        double wgt = Calc_Weight(m_2sigma, c_2sigma, i);
        wgts.push_back(wgt);
        if (isDebug){
            std::cout<<"\t"<<i<<" "<<wgt<<std::endl;
        }
    }

    return wgts;
}

double CCProtonPi0_QSqFitter::Calc_Slope(double x1, double x2, double y1, double y2)
{
    double m = (y2 - y1) / (x2 - x1);
    return m;
}

double CCProtonPi0_QSqFitter::Calc_Constant(double m, double x, double y)
{
    // y = mx + c
    double c = y - m*x;
    return c;
}

double CCProtonPi0_QSqFitter::Calc_Weight(double m, double c, double x)
{
    double y = m*x + c;
    return y;
}

int CCProtonPi0_QSqFitter::GetMinChiSq_DeltaFactor()
{
    FillChiSqVector_DeltaFactor();
    int min = FindMinChiSq(ChiSqVector_DeltaFactor);

    std::cout<<"min = "<<min<<std::endl;
    return min;
}

int CCProtonPi0_QSqFitter::GetMinChiSq(bool isAreaNorm)
{
    if (isAreaNorm){
        FillChiSqVector_AreaNorm(ChiSqVector_dn, false);
        FillChiSqVector_AreaNorm(ChiSqVector_up, true);
        int min_dn = FindMinChiSq(ChiSqVector_dn);
        int min_up = FindMinChiSq(ChiSqVector_up);

        std::cout<<"Minimum ChiSqs for Area Normalized"<<std::endl;
        std::cout<<"min_dn = "<<min_dn<<" min_up = "<<min_up<<std::endl;
        return min_up;
    }else{
        FillChiSqVector(ChiSqVector_dn, false);
        FillChiSqVector(ChiSqVector_up, true);
        int min_dn = FindMinChiSq(ChiSqVector_dn);
        int min_up = FindMinChiSq(ChiSqVector_up);

        std::cout<<"min_dn = "<<min_dn<<" min_up = "<<min_up<<std::endl;
        return min_up;
    }
}

int CCProtonPi0_QSqFitter::FindMinChiSq(std::vector<double> &ChiSqVector)
{
    double current_min = ChiSqVector[0];
    int current_min_ind = 0;
    
    for (unsigned int  i = 1; i < ChiSqVector.size(); ++i){
        if (ChiSqVector[i] < current_min){
            current_min = ChiSqVector[i];
            current_min_ind = i;
        }
    }

    return current_min_ind;
}

void CCProtonPi0_QSqFitter::NormalizeHistogram(TH1D* h)
{
    int NBins = h->GetNbinsX();
    double area = h->Integral();
    double nOverFlow = h->GetBinContent(NBins+1);
    double nUnderFlow = h->GetBinContent(0);
    h->Scale(1/(area+nOverFlow+nUnderFlow)); // Scale only on CentralValue
}

void CCProtonPi0_QSqFitter::FillChiSqVector_SB(std::vector<double> &ChiSqVector, bool isUpShift)
{
    // ------------------------------------------------------------------------
    // Read Files
    // ------------------------------------------------------------------------
    std::string data_dir_SB_Michel = "/minerva/data/users/oaltinok/NTupleAnalysis/Data/Analyzed/CutHistograms_Michel.root";
    std::string data_dir_SB_pID = "/minerva/data/users/oaltinok/NTupleAnalysis/Data/Analyzed/CutHistograms_pID.root";
    std::string data_dir_SB_LowInvMass = "/minerva/data/users/oaltinok/NTupleAnalysis/Data/Analyzed/CutHistograms_LowInvMass.root";
    std::string data_dir_SB_HighInvMass = "/minerva/data/users/oaltinok/NTupleAnalysis/Data/Analyzed/CutHistograms_HighInvMass.root";
    std::string data_dir = "/minerva/data/users/oaltinok/NTupleAnalysis/Data/Analyzed/Interaction.root";
 
    std::string mc_dir_SB_Michel = "/minerva/data/users/oaltinok/NTupleAnalysis/MC/Analyzed/CutHistograms_Michel.root";
    std::string mc_dir_SB_pID = "/minerva/data/users/oaltinok/NTupleAnalysis/MC/Analyzed/CutHistograms_pID.root";
    std::string mc_dir_SB_LowInvMass = "/minerva/data/users/oaltinok/NTupleAnalysis/MC/Analyzed/CutHistograms_LowInvMass.root";
    std::string mc_dir_SB_HighInvMass = "/minerva/data/users/oaltinok/NTupleAnalysis/MC/Analyzed/CutHistograms_HighInvMass.root";
    std::string mc_dir = "/minerva/data/users/oaltinok/NTupleAnalysis/MC/Analyzed/Interaction.root";
   
    TFile* f_data_SB_Michel = new TFile(data_dir_SB_Michel.c_str());
    TFile* f_data_SB_pID = new TFile(data_dir_SB_pID.c_str());
    TFile* f_data_SB_LowInvMass = new TFile(data_dir_SB_LowInvMass.c_str());
    TFile* f_data_SB_HighInvMass = new TFile(data_dir_SB_HighInvMass.c_str());
    TFile* f_data = new TFile(data_dir.c_str());
 
    TFile* f_mc_SB_Michel = new TFile(mc_dir_SB_Michel.c_str());
    TFile* f_mc_SB_pID = new TFile(mc_dir_SB_pID.c_str());
    TFile* f_mc_SB_LowInvMass = new TFile(mc_dir_SB_LowInvMass.c_str());
    TFile* f_mc_SB_HighInvMass = new TFile(mc_dir_SB_HighInvMass.c_str());
    TFile* f_mc = new TFile(mc_dir.c_str());

    // ------------------------------------------------------------------------
    // Get Histograms
    // ------------------------------------------------------------------------
    MnvH1D* h_data_SB_Michel = GetMnvH1D(f_data_SB_Michel, "SideBand_QSq_0");
    MnvH1D* h_data_SB_pID = GetMnvH1D(f_data_SB_pID, "SideBand_QSq_0");
    MnvH1D* h_data_SB_LowInvMass = GetMnvH1D(f_data_SB_LowInvMass, "SideBand_QSq_0");
    MnvH1D* h_data_SB_HighInvMass = GetMnvH1D(f_data_SB_HighInvMass, "SideBand_QSq_0");
    MnvH1D* h_data = GetMnvH1D(f_data, "QSq_MaRES_0");

    MnvH1D* h_mc_SB_Michel = GetMnvH1D(f_mc_SB_Michel, "SideBand_QSq_0");
    MnvH1D* h_mc_SB_pID = GetMnvH1D(f_mc_SB_pID, "SideBand_QSq_0");
    MnvH1D* h_mc_SB_LowInvMass = GetMnvH1D(f_mc_SB_LowInvMass, "SideBand_QSq_0");
    MnvH1D* h_mc_SB_HighInvMass = GetMnvH1D(f_mc_SB_HighInvMass, "SideBand_QSq_0");
    MnvH1D* h_mc = GetMnvH1D(f_mc, "QSq_MaRES_0");

    // ------------------------------------------------------------------------
    // Get Universes
    // ------------------------------------------------------------------------
    std::string err_name = isUpShift ? "HighMaRES" : "LowMaRES";

    MnvVertErrorBand* err_data_SB_Michel = h_data_SB_Michel->GetVertErrorBand(err_name);
    MnvVertErrorBand* err_data_SB_pID = h_data_SB_pID->GetVertErrorBand(err_name);
    MnvVertErrorBand* err_data_SB_LowInvMass = h_data_SB_LowInvMass->GetVertErrorBand(err_name);
    MnvVertErrorBand* err_data_SB_HighInvMass = h_data_SB_HighInvMass->GetVertErrorBand(err_name);
    MnvVertErrorBand* err_data = h_data->GetVertErrorBand(err_name);

    MnvVertErrorBand* err_mc_SB_Michel = h_mc_SB_Michel->GetVertErrorBand(err_name);
    MnvVertErrorBand* err_mc_SB_pID = h_mc_SB_pID->GetVertErrorBand(err_name);
    MnvVertErrorBand* err_mc_SB_LowInvMass = h_mc_SB_LowInvMass->GetVertErrorBand(err_name);
    MnvVertErrorBand* err_mc_SB_HighInvMass = h_mc_SB_HighInvMass->GetVertErrorBand(err_name);
    MnvVertErrorBand* err_mc = h_mc->GetVertErrorBand(err_name);

    std::vector<TH1D*> unv_data_SB_Michel = err_data_SB_Michel->GetHists();
    std::vector<TH1D*> unv_data_SB_pID = err_data_SB_pID->GetHists();
    std::vector<TH1D*> unv_data_SB_LowInvMass = err_data_SB_LowInvMass->GetHists();
    std::vector<TH1D*> unv_data_SB_HighInvMass = err_data_SB_HighInvMass->GetHists();
    std::vector<TH1D*> unv_data = err_data->GetHists();
 
    std::vector<TH1D*> unv_mc_SB_Michel = err_mc_SB_Michel->GetHists();
    std::vector<TH1D*> unv_mc_SB_pID = err_mc_SB_pID->GetHists();
    std::vector<TH1D*> unv_mc_SB_LowInvMass = err_mc_SB_LowInvMass->GetHists();
    std::vector<TH1D*> unv_mc_SB_HighInvMass = err_mc_SB_HighInvMass->GetHists();
    std::vector<TH1D*> unv_mc = err_mc->GetHists();
  
    // ------------------------------------------------------------------------
    // Calculate Global Chi Squre  
    // ------------------------------------------------------------------------
    for (int i = 0; i < 201; ++i){
        double Global_ChiSq = 0.0;
        //Global_ChiSq += Calc_ChiSq(unv_data_SB_Michel[i], unv_mc_SB_Michel[i]);
        //Global_ChiSq += Calc_ChiSq(unv_data_SB_pID[i], unv_mc_SB_pID[i]);
        //Global_ChiSq += Calc_ChiSq(unv_data_SB_LowInvMass[i], unv_mc_SB_LowInvMass[i]);
        //Global_ChiSq += Calc_ChiSq(unv_data_SB_HighInvMass[i], unv_mc_SB_HighInvMass[i]);
        Global_ChiSq += Calc_ChiSq(unv_data[i], unv_mc[i]);
        ChiSqVector.push_back(Global_ChiSq);
    }

    if (isDebug){
        std::cout<<"ChiSq Vector for "<<err_name<<std::endl;
        for (unsigned int i = 0; i < ChiSqVector.size(); ++i){
            std::cout<<"\t"<<i<<" "<<ChiSqVector[i]<<std::endl;
        }
    }

    // Clean Memory
    delete h_data_SB_Michel;
    delete h_data_SB_pID;
    delete h_data_SB_LowInvMass;
    delete h_data_SB_HighInvMass;
    delete h_data;

    delete h_mc_SB_Michel;
    delete h_mc_SB_pID;
    delete h_mc_SB_LowInvMass;
    delete h_mc_SB_HighInvMass;
    delete h_mc;

    delete f_data_SB_Michel;
    delete f_data_SB_pID;
    delete f_data_SB_LowInvMass;
    delete f_data_SB_HighInvMass;
    delete f_data;

    delete f_mc_SB_Michel;
    delete f_mc_SB_pID;
    delete f_mc_SB_LowInvMass;
    delete f_mc_SB_HighInvMass;
    delete f_mc;
}

void CCProtonPi0_QSqFitter::FillChiSqVector_DeltaFactor()
{
    // ------------------------------------------------------------------------
    // Read Files
    // ------------------------------------------------------------------------
    std::string data_dir = Folder_List::rootDir_CrossSection_data;
    std::string mc_dir = Folder_List::rootDir_CrossSection_mc;
   
    TFile* f_data = new TFile(data_dir.c_str());
    TFile* f_mc = new TFile(mc_dir.c_str());

    // ------------------------------------------------------------------------
    // Get Histograms
    // ------------------------------------------------------------------------
    MnvH1D* h_data = GetMnvH1D(f_data, "QSq_xsec");
    MnvH1D* h_mc = GetMnvH1D(f_mc, "QSq_xsec");
 
    // ------------------------------------------------------------------------
    // Get Universes
    // ------------------------------------------------------------------------
    std::string err_name = "DeltaFactor";
    MnvVertErrorBand* err_data = h_data->GetVertErrorBand(err_name);
    MnvVertErrorBand* err_mc = h_mc->GetVertErrorBand(err_name);

    std::vector<TH1D*> unv_data = err_data->GetHists();
    std::vector<TH1D*> unv_mc = err_mc->GetHists();
  
    // ------------------------------------------------------------------------
    // Calculate Chi Squre  
    // ------------------------------------------------------------------------
    for (unsigned int i = 0; i < unv_mc.size(); ++i){
        double ChiSq = Calc_ChiSq_Delta(unv_data[i], unv_mc[i]);
        ChiSqVector_DeltaFactor.push_back(ChiSq);
    }

    if (isDebug){
        std::cout<<"ChiSq Vector for "<<err_name<<std::endl;
        for (unsigned int i = 0; i < ChiSqVector_DeltaFactor.size(); ++i){
            std::cout<<"\t"<<i<<" "<<ChiSqVector_DeltaFactor[i]<<std::endl;
        }
    }

    // Clean Memory
    delete h_data;
    delete h_mc;
    delete f_data;
    delete f_mc;
}

double CCProtonPi0_QSqFitter::GetSmallestBinWidth(MnvH1D* hist)
{
    double smallest = 99999999;
    int nBins = hist->GetNbinsX();
    for (int i = 0; i <= nBins; ++i){
        double current = hist->GetBinWidth(i);
        if (current < smallest) smallest = current;
    }

    return smallest;
}

void CCProtonPi0_QSqFitter::FillChiSqVector(std::vector<double> &ChiSqVector, bool isUpShift)
{
    // ------------------------------------------------------------------------
    // Read Files
    // ------------------------------------------------------------------------
    std::string data_dir = "/minerva/data/users/oaltinok/NTupleAnalysis_MaRES_Fit_XSecs/Data/Analyzed/CrossSection.root";
    std::string mc_dir = "/minerva/data/users/oaltinok/NTupleAnalysis_MaRES_Fit_XSecs/MC/Analyzed/CrossSection.root";

    TFile* f_data = new TFile(data_dir.c_str());
    TFile* f_mc = new TFile(mc_dir.c_str());

    // ------------------------------------------------------------------------
    // Get Histograms
    // ------------------------------------------------------------------------
    MnvH1D* h_data = GetMnvH1D(f_data, "QSq_xsec");
    MnvH1D* h_mc = GetMnvH1D(f_mc, "QSq_xsec");

    // ------------------------------------------------------------------------
    // Get Universes
    // ------------------------------------------------------------------------
    std::string err_name = isUpShift ? "HighMaRES" : "LowMaRES";

    MnvVertErrorBand* err_data = h_data->GetVertErrorBand(err_name);
    MnvVertErrorBand* err_mc = h_mc->GetVertErrorBand(err_name);

    std::vector<TH1D*> unv_data = err_data->GetHists();
    std::vector<TH1D*> unv_mc = err_mc->GetHists();
  
    // ------------------------------------------------------------------------
    // Calculate Chi Squre  
    // ------------------------------------------------------------------------
    for (int i = 0; i < 201; ++i){
        double ChiSq = Calc_ChiSq(unv_data[i], unv_mc[i]);
        ChiSqVector.push_back(ChiSq);
    }

    if (isDebug){
        std::cout<<"ChiSq Vector for "<<err_name<<std::endl;
        for (unsigned int i = 0; i < ChiSqVector.size(); ++i){
            std::cout<<"\t"<<i<<" "<<ChiSqVector[i]<<std::endl;
        }
    }

    // Clean Memory
    delete h_data;
    delete h_mc;
    delete f_data;
    delete f_mc;

}

void CCProtonPi0_QSqFitter::FillChiSqVector_AreaNorm(std::vector<double> &ChiSqVector, bool isUpShift)
{
    // ------------------------------------------------------------------------
    // Read Files
    // ------------------------------------------------------------------------
    std::string data_dir = "/minerva/data/users/oaltinok/NTupleAnalysis_MaRES_Fit_XSecs/Data/Analyzed/CrossSection.root";
    std::string mc_dir = "/minerva/data/users/oaltinok/NTupleAnalysis_MaRES_Fit_XSecs/MC/Analyzed/CrossSection.root";

    TFile* f_data = new TFile(data_dir.c_str());
    TFile* f_mc = new TFile(mc_dir.c_str());

    // ------------------------------------------------------------------------
    // Get Histograms
    // ------------------------------------------------------------------------
    MnvH1D* h_data = GetMnvH1D(f_data, "QSq_xsec");
    MnvH1D* h_mc = GetMnvH1D(f_mc, "QSq_xsec");

    // ------------------------------------------------------------------------
    // Get Universes
    // ------------------------------------------------------------------------
    std::string err_name = isUpShift ? "HighMaRES" : "LowMaRES";

    MnvVertErrorBand* err_data = h_data->GetVertErrorBand(err_name);
    MnvVertErrorBand* err_mc = h_mc->GetVertErrorBand(err_name);

    std::vector<TH1D*> unv_data = err_data->GetHists();
    std::vector<TH1D*> unv_mc = err_mc->GetHists();
  
    // ------------------------------------------------------------------------
    // Calculate Chi Squre  
    // ------------------------------------------------------------------------
    for (int i = 0; i < 201; ++i){
        AreaNormalize(unv_data[i], unv_mc[i]);
        double ChiSq = Calc_ChiSq(unv_data[i], unv_mc[i]);
        ChiSqVector.push_back(ChiSq);
    }

    if (isDebug){
        std::cout<<"ChiSq Vector for "<<err_name<<std::endl;
        for (unsigned int i = 0; i < ChiSqVector.size(); ++i){
            std::cout<<"\t"<<i<<" "<<ChiSqVector[i]<<std::endl;
        }
    }

    // Clean Memory
    delete h_data;
    delete h_mc;
    delete f_data;
    delete f_mc;

}

double CCProtonPi0_QSqFitter::Calc_ChiSq(TH1D* data, TH1D* MC)
{
    // Do not Scale MC for XSec
    //MC->Scale(POT_ratio);

    int nBins = data->GetNbinsX();
    double ChiSq = 0.0;

    // Do not use first 2 bins
    for (int i = 3; i <= nBins; ++i){
        double nData = data->GetBinContent(i);
        double nMC = MC->GetBinContent(i);
        
        ChiSq += std::pow((nData-nMC),2) / nMC;
    }

    return ChiSq;
}

void CCProtonPi0_QSqFitter::AreaNormalize(TH1D* data, TH1D* MC)
{
    // Do not Scale MC for XSec
    //MC->Scale(POT_ratio);

    int nBins = data->GetNbinsX();
    double area_data = data->Integral(3,nBins,"width");
    double area_mc = MC->Integral(3,nBins,"width");
    double ratio = area_data/area_mc;
    MC->Scale(ratio);
}

double CCProtonPi0_QSqFitter::Calc_ChiSq_Delta(TH1D* data, TH1D* MC)
{
    // Do not Scale MC for XSec
    //MC->Scale(POT_ratio);

    //int nBins = data->GetNbinsX();
    double ChiSq = 0.0;

    for (int i = 1; i <= 3; ++i){
        double nData = data->GetBinContent(i);
        double nMC = MC->GetBinContent(i);
        
        ChiSq += std::pow((nData-nMC),2) / nMC;
    }

    return ChiSq;
}
#endif

