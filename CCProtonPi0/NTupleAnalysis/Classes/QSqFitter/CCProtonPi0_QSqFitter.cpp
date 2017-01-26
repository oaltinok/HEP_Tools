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

int CCProtonPi0_QSqFitter::GetMinChiSq()
{
    FillChiSqVector(ChiSqVector_dn, false);
    FillChiSqVector(ChiSqVector_up, true);
    int min_dn = FindMinChiSq(ChiSqVector_dn);
    int min_up = FindMinChiSq(ChiSqVector_up);

    std::cout<<"min_up = "<<min_up<<" min_dn = "<<min_dn<<std::endl;
    return min_up;
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

void CCProtonPi0_QSqFitter::FillChiSqVector(std::vector<double> &ChiSqVector, bool isUpShift)
{
    std::string var_type = isUpShift ? "HighMaRES" : "LowMaRES";
   
    std::string var_name = "QSq_" + var_type;
    std::string bckg_name = "pi0_invMass_" + var_type;

    TFile* f_data = new TFile(Folder_List::rootDir_Interaction_data.c_str());
    TFile* f_mc = new TFile(Folder_List::rootDir_Interaction_mc.c_str());

    double nData_All;
    double nData_Signal;
    double nData_Bckg;

    for (int i = 0; i <= 200; ++i){
        // --------------------------------------------------------------------
        // Calculate Total N(Background) in this Universe
        // --------------------------------------------------------------------
        // Data is Same in Each Universe
        TH1D* invMass_data = new TH1D( * dynamic_cast<TH1D*>(f_data->Get(Form("%s_%d", "pi0_invMass_All", 0))));
        nData_All = invMass_data->Integral(); 

        TH1D* invMass_mc = new TH1D( * dynamic_cast<TH1D*>(f_mc->Get(Form("%s_%d", bckg_name.c_str(), i))));
        invMass_mc->Scale(POT_ratio);
        nData_Bckg = invMass_mc->Integral();

        invMass_data->Add(invMass_mc, -1);
        nData_Signal = invMass_data->Integral();
        delete invMass_mc;
        delete invMass_data;
        
        //std::cout<<"nData_All = "<<nData_All<<std::endl;
        //std::cout<<"nData_Signal = "<<nData_Signal<<std::endl;
        //std::cout<<"nData_Bckg = "<<nData_Bckg<<std::endl;

        // --------------------------------------------------------------------
        // Background Subtract in this Universe 
        // --------------------------------------------------------------------
        // Get Background Subtracted Data
        TH1D* h_data = new TH1D( * dynamic_cast<TH1D*>(f_data->Get(Form("%s_%d", var_name.c_str(), 0))));
        TH1D* h_mc_bckg = new TH1D( * dynamic_cast<TH1D*>(f_mc->Get(Form("%s_%s_%d", var_name.c_str(), "Bckg",i))));
        NormalizeHistogram(h_mc_bckg);
        h_mc_bckg->Scale(nData_Bckg);
        h_data->Add(h_mc_bckg, -1);

        // Get MC Signal
        TH1D* h_mc_signal = new TH1D( * dynamic_cast<TH1D*>(f_mc->Get(Form("%s_%d", var_name.c_str(), i))));
        h_mc_signal->Scale(POT_ratio);
        double ChiSq = Calc_ChiSq(h_data, h_mc_signal);
        ChiSqVector.push_back(ChiSq);

        delete h_data;
        delete h_mc_signal;
        delete h_mc_bckg;
    }

    if (isDebug){
        std::cout<<"ChiSq Vector"<<std::endl;
        for (unsigned int i = 0; i < ChiSqVector.size(); ++i){
            std::cout<<"\t"<<i<<" "<<ChiSqVector[i]<<std::endl;
        }
    }
}

double CCProtonPi0_QSqFitter::Calc_ChiSq(TH1D* data, TH1D* MC)
{
    int nBins = data->GetNbinsX();
    double ChiSq = 0.0;

    for (int i = 3; i < nBins-5; ++i){
        double nData = data->GetBinContent(i);
        double nMC = MC->GetBinContent(i);
        
        ChiSq += std::pow((nData-nMC),2) / nMC;
    }

    // Combine Last 5 Bins
    double nData_Last5 = 0.0;
    double nMC_Last5 = 0.0;
    for (int i = nBins - 5; i <= nBins; ++i){
        nData_Last5 += data->GetBinContent(i);
        nMC_Last5 += MC->GetBinContent(i);
    }

    ChiSq += std::pow((nData_Last5-nMC_Last5),2) / nMC_Last5;

    return ChiSq;
}

#endif

