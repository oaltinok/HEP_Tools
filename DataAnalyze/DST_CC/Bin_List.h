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
const int NBINS_energyNeutrino = 40;
const double MIN_energyNeutrino = 1;
const double MAX_energyNeutrino = 5;
const double WIDTH_energyNeutrino = getBinWidth(MIN_energyNeutrino,MAX_energyNeutrino,NBINS_energyNeutrino);



// W [GeV]
const int NBINS_W = 40;
const double MIN_W = 0;
const double MAX_W = 4;



#endif