#ifndef CCProtonPi0_RandNumGenerator_h
#define CCProtonPi0_RandNumGenerator_h

#include <TRandom1.h>

#include "../NTupleAnalysis/CCProtonPi0_NTupleAnalysis.h"

class CCProtonPi0_RandNumGenerator : public CCProtonPi0_NTupleAnalysis 
{
    public:
        CCProtonPi0_RandNumGenerator();
        ~CCProtonPi0_RandNumGenerator();

        std::vector<double> GetNormalRandomVector();
        std::vector<double> GetRandomShifts(double sigma);
        
        void Print_NormalRandomVector();

    private:
        TRandom1 generator;
        
        void FillNormalRandomVector();
        std::vector<double> normal_random_vector;

};

#endif
