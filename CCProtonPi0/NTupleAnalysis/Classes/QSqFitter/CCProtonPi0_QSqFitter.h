#ifndef CCProtonPi0_QSqFitter_h
#define CCProtonPi0_QSqFitter_h

#include "../NTupleAnalysis/CCProtonPi0_NTupleAnalysis.h"
#include "../../Libraries/Folder_List.h"

class CCProtonPi0_QSqFitter: public CCProtonPi0_NTupleAnalysis
{
    public:
        CCProtonPi0_QSqFitter();
        std::vector<double> GetWeights(double genie_1sigma, double genie_2sigma);

        int GetMinChiSq();
        std::vector<double> MaRESVector_up;
        std::vector<double> ChiSqVector_up;
        std::vector<double> MaRESVector_dn;
        std::vector<double> ChiSqVector_dn;

    private:
        static const double GENIE_MaRES;
        static const double x1_1sigma; // First Element -- 0 (CV)
        static const double x2_1sigma; // 1Sigma Element -- 50
        static const double x1_2sigma; // 1Sigma Element -- 50;
        static const double x2_2sigma; // 2Sigma Element -- 100;

        bool isDebug;

        int FindMinChiSq(std::vector<double> &ChiSqVector);
        void FillMaRESVector(std::vector<double> &MaRESVector, double one_sigma, double two_sigma);
        void FillChiSqVector(std::vector<double> &ChiSqVector, bool isUpShift);
        void NormalizeHistogram(TH1D* h);
        double Calc_ChiSq(TH1D* data, TH1D* MC);
        double Calc_Slope(double x1, double x2, double y1, double y2);
        double Calc_Constant(double m, double x, double y);
        double Calc_Weight(double m, double c, double x);
};


#endif

