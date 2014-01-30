/*
    See SingleBin.h header for Class Information
*/
#include "SingleBin.h"

using namespace std;

SingleBin::SingleBin()
{

}

void SingleBin::setBin(int in_nBins, double in_min, double in_max)
{
    set_nBins(in_nBins);
    set_min(in_min);
    set_max(in_max);

}

int SingleBin::get_nBins()
{
    return nBins;
}

double SingleBin::get_max()
{
    return max;
}

double SingleBin::get_min()
{
    return min;
}

double SingleBin::get_width()
{
    width = (max - min) / ((double)nBins);

    return width;
}

void SingleBin::set_nBins(int input)
{
    nBins = input;
}

void SingleBin::set_max(double input)
{
    max = input;
}

void SingleBin::set_min(double input)
{
    min = input;
}
