/*
    See CCProtonPi0_SingleBin.h header for Class Information
*/
#ifndef CCProtonPi0_SingleBin_cpp
#define CCProtonPi0_SingleBin_cpp

#include "CCProtonPi0_SingleBin.h"

using namespace std;

CCProtonPi0_SingleBin::CCProtonPi0_SingleBin()
{

}

void CCProtonPi0_SingleBin::setBin(int in_nBins, double in_min, double in_max)
{
    set_nBins(in_nBins);
    set_min(in_min);
    set_max(in_max);

}

int CCProtonPi0_SingleBin::get_nBins()
{
    return nBins;
}

double CCProtonPi0_SingleBin::get_max()
{
    return max;
}

double CCProtonPi0_SingleBin::get_min()
{
    return min;
}

double CCProtonPi0_SingleBin::get_width()
{
    width = (max - min) / ((double)nBins);

    return width;
}

void CCProtonPi0_SingleBin::set_nBins(int input)
{
    nBins = input;
}

void CCProtonPi0_SingleBin::set_max(double input)
{
    max = input;
}

void CCProtonPi0_SingleBin::set_min(double input)
{
    min = input;
}

#endif
