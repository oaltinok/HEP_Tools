#ifndef CCProtonPi0_RandNumGenerator_cpp
#define CCProtonPi0_RandNumGenerator_cpp

#include "CCProtonPi0_RandNumGenerator.h"

CCProtonPi0_RandNumGenerator::CCProtonPi0_RandNumGenerator() : 
    CCProtonPi0_NTupleAnalysis(), 
    generator(1987)
{
    FillNormalRandomVector();
}

CCProtonPi0_RandNumGenerator::~CCProtonPi0_RandNumGenerator()
{
    // Do Nothing!
}

void CCProtonPi0_RandNumGenerator::FillNormalRandomVector()
{
    const int n_universe = n_lateral_universes;
    
    for (int i = 0; i < n_universe; ++i) {    
        double val = generator.Gaus(0.0, 1.0);
        normal_random_vector.push_back(val);
    }
}

std::vector<double> CCProtonPi0_RandNumGenerator::GetRandomShifts(double sigma)
{

    std::vector<double> random_shifts;

    for (unsigned int i = 0; i < normal_random_vector.size(); ++i){
        double temp = sigma * normal_random_vector[i];
        random_shifts.push_back(temp);
    }
    
    return random_shifts;
}

void CCProtonPi0_RandNumGenerator::Print_NormalRandomVector()
{
    std::cout<<"RandNumGenerator Normal Random Vector"<<std::endl;
    for (unsigned int i = 0; i < normal_random_vector.size(); ++i){
        std::cout<<normal_random_vector[i]<<std::endl;
    }
}

std::vector<double> CCProtonPi0_RandNumGenerator::GetNormalRandomVector()
{
    return normal_random_vector;
}

#endif
