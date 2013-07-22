#ifndef Bin_List.h
#define Bin_List.h

using namespace std;

double getBinWidth(double min,double max,int nbins)
{
    double width;

    width = (max - min) / ((double)nbins);

    return width;

}

const double mev_to_gev = 0.001;
const double mevSq_to_gevSq = 1e-6;

// Neutrino Energy [GeV]
const int NBINS_Ev = 40;
const double MIN_Ev = 1;
const double MAX_Ev = 5;
const double WIDTH_Ev = getBinWidth(MIN_Ev,MAX_Ev,NBINS_Ev);



// W [GeV]
const int NBINS_W = 40;
const double MIN_W = 0;
const double MAX_W = 4;
const double WIDTH_W = getBinWidth(MIN_W,MAX_W,NBINS_W);



#endif