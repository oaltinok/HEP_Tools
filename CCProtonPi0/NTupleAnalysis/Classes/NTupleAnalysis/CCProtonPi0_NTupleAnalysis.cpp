#ifndef CCProtonPi0_NTupleAnalysis_cpp
#define CCProtonPi0_NTupleanalysis_cpp

#include "CCProtonPi0_NTupleAnalysis.h"

using namespace std;

// Initialize Constants
const string CCProtonPi0_NTupleAnalysis::version = "v2_49";
const double CCProtonPi0_NTupleAnalysis::SENTINEL = -9.9;
const double CCProtonPi0_NTupleAnalysis::MeV_to_GeV = pow(10,-3);
const double CCProtonPi0_NTupleAnalysis::MeVSq_to_GeVSq = pow(10,-6);
const double CCProtonPi0_NTupleAnalysis::mm_to_cm = pow(10,-1);

CCProtonPi0_NTupleAnalysis::CCProtonPi0_NTupleAnalysis()
{
    // Required for MINERvA Framework Classes
    ROOT::Cintex::Cintex::Enable();
}

void CCProtonPi0_NTupleAnalysis::OpenTextFile(string file_name, std::ofstream &file)
{
    file.open(file_name.c_str());
    if (!file.is_open()){
        cerr<<"Cannot open output text file: "<<file_name<<endl;
        exit(1);
    }else{
        cout<<"\t"<<file_name<<endl;
    }
}

#endif

