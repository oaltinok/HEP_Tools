#ifndef Plotter_cpp
#define Plotter_cpp


#include "Plotter.h"

using namespace std;

void Plotter::SetRootFile(TFile* input)
{
    f_Root = input;
}

Plotter::Plotter()
{
    isDebugging = false;    
}


#endif
