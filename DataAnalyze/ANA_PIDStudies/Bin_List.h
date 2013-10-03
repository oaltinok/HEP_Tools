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


// Pion Score
const int NBINS_pID = 50;
const double MIN_pID = 0;
const double MAX_pID = 1;
const double WIDTH_pID  = getBinWidth(MIN_pID ,MAX_pID ,NBINS_pID );

// Node #
const int NBINS_Nodes = 41;
const double MIN_Nodes  = 0;
const double MAX_Nodes  = 40;
const double WIDTH_Nodes  = 1;

// dEdX
const int NBINS_dEdX = 40;
const double MIN_dEdX  = 0;
const double MAX_dEdX  = 2;
const double WIDTH_dEdX  = getBinWidth(MIN_dEdX ,MAX_dEdX ,NBINS_dEdX );

// delta_dEdX
const int NBINS_delta_dEdX = 80;
const double MIN_delta_dEdX  = -2;
const double MAX_delta_dEdX  = 2;
const double WIDTH_delta_dEdX  = getBinWidth(MIN_delta_dEdX ,MAX_delta_dEdX ,NBINS_delta_dEdX );

// Proton Energy [MeV]
const int NBINS_Pp = 50;
const double MIN_Pp = 200;
const double MAX_Pp = 800;
const double WIDTH_Pp = getBinWidth(MIN_Pp,MAX_Pp,NBINS_Pp);

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