/*
================================================================================
Library: Bin_List
    Binning information for ROOT Histograms
    List of Commonly used bins for the analysis 
    
    Main Directory:
        Libraries/
    
    Usage:
        > #include "Libraries/Bin_List.h" 
        > Bin_List::NBINS_Error
            
    
    Last Revision: 2014_01_27
================================================================================
*/
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

// Error
const int NBINS_Error = 400;
const double MIN_Error= -2;
const double MAX_Error  = 2;
const double WIDTH_Error  = getBinWidth(MIN_Error ,MAX_Error ,NBINS_Error);

// Particle Multiplicity
const int NBINS_Multiplicity = 11;
const double MIN_Multiplicity = 0;
const double MAX_Multiplicity  = 10;
const double WIDTH_Multiplicity  = getBinWidth(MIN_Multiplicity ,MAX_Multiplicity ,NBINS_Multiplicity);

// Vertex Z Position
const int NBINS_PosZ = 5507;
const double MIN_PosZ = 4480;
const double MAX_PosZ  = 9987;
const double WIDTH_PosZ  = getBinWidth(MIN_PosZ ,MAX_PosZ ,NBINS_PosZ);

// Proton Px [MeV]
const int NBINS_Proton_Px = 24;
const double MIN_Proton_Px = -1200;
const double MAX_Proton_Px  = 1200;
const double WIDTH_Proton_Px  = getBinWidth(MIN_Proton_Px ,MAX_Proton_Px ,NBINS_Proton_Px);

// Proton Py [MeV]
const int NBINS_Proton_Py = 24;
const double MIN_Proton_Py = -1200;
const double MAX_Proton_Py  = 1200;
const double WIDTH_Proton_Py  = getBinWidth(MIN_Proton_Py ,MAX_Proton_Py ,NBINS_Proton_Py);

// Proton Pz [MeV]
const int NBINS_Proton_Pz = 20;
const double MIN_Proton_Pz = -500;
const double MAX_Proton_Pz  = 1500;
const double WIDTH_Proton_Pz  = getBinWidth(MIN_Proton_Pz ,MAX_Proton_Pz ,NBINS_Proton_Pz);

// Proton PTotal [MeV]
const int NBINS_Proton_PTotal = 40;
const double MIN_Proton_PTotal = 0;
const double MAX_Proton_PTotal  = 4000;
const double WIDTH_Proton_PTotal  = getBinWidth(MIN_Proton_PTotal ,MAX_Proton_PTotal ,NBINS_Proton_PTotal);


// Nucleon PTotal [MeV]
const int NBINS_Nucleon_PTotal = 40;
const double MIN_Nucleon_PTotal = 0;
const double MAX_Nucleon_PTotal  = 4000;
const double WIDTH_Nucleon_PTotal  = getBinWidth(MIN_Nucleon_PTotal ,MAX_Nucleon_PTotal ,NBINS_Nucleon_PTotal);

// Pion PTotal [MeV]
const int NBINS_Pion_PTotal = 40;
const double MIN_Pion_PTotal = 0;
const double MAX_Pion_PTotal  = 4000;
const double WIDTH_Pion_PTotal  = getBinWidth(MIN_Pion_PTotal ,MAX_Pion_PTotal ,NBINS_Pion_PTotal);

// Proton Energy [MeV]
const int NBINS_Proton_Energy = 21;
const double MIN_Proton_Energy = 900;
const double MAX_Proton_Energy  = 3000;
const double WIDTH_Proton_Energy  = getBinWidth(MIN_Proton_Energy ,MAX_Proton_Energy ,NBINS_Proton_Energy);

// Muon Px [MeV]
const int NBINS_Muon_Px = 30;
const double MIN_Muon_Px = -1500;
const double MAX_Muon_Px  = 1500;
const double WIDTH_Muon_Px  = getBinWidth(MIN_Muon_Px ,MAX_Muon_Px ,NBINS_Muon_Px);

// Muon Py [MeV]
const int NBINS_Muon_Py = 30;
const double MIN_Muon_Py = -1500;
const double MAX_Muon_Py  = 1500;
const double WIDTH_Muon_Py  = getBinWidth(MIN_Muon_Py ,MAX_Muon_Py ,NBINS_Muon_Py);

// Muon Pz [MeV]
const int NBINS_Muon_Pz = 160;
const double MIN_Muon_Pz = -1000;
const double MAX_Muon_Pz  = 15000;
const double WIDTH_Muon_Pz  = getBinWidth(MIN_Muon_Pz ,MAX_Muon_Pz ,NBINS_Muon_Pz);

// Muon PTotal [MeV]
const int NBINS_Muon_PTotal = 100;
const double MIN_Muon_PTotal = 0;
const double MAX_Muon_PTotal  = 10000;
const double WIDTH_Muon_PTotal  = getBinWidth(MIN_Muon_PTotal ,MAX_Muon_PTotal ,NBINS_Muon_PTotal);

// Muon Energy [MeV]
const int NBINS_Muon_Energy = 50;
const double MIN_Muon_Energy = 0;
const double MAX_Muon_Energy  = 5000;
const double WIDTH_Muon_Energy  = getBinWidth(MIN_Muon_Energy ,MAX_Muon_Energy ,NBINS_Muon_Energy);

// Neutrino Energy [GeV]
const int NBINS_Ev = 50;
const double MIN_Ev = 0;
const double MAX_Ev = 5000;
const double WIDTH_Ev = getBinWidth(MIN_Ev,MAX_Ev,NBINS_Ev);

// Neutrino Energy [GeV]
const int NBINS_neutrino_energy_diff = 80;
const double MIN_neutrino_energy_diff = -4000;
const double MAX_neutrino_energy_diff = 4000;
const double WIDTH_neutrino_energy_diff = getBinWidth(MIN_neutrino_energy_diff,MAX_neutrino_energy_diff,NBINS_neutrino_energy_diff);

// W [GeV]
const int NBINS_W = 40;
const double MIN_W = 0;
const double MAX_W = 4;
const double WIDTH_W = getBinWidth(MIN_W,MAX_W,NBINS_W);

// Angles(Theta = angle wrt Z)
const int NBINS_Angle = 360;
const double MIN_Angle = -3.14;
const double MAX_Angle = 3.14;
const double WIDTH_Angle = getBinWidth(MIN_Angle,MAX_Angle,NBINS_Angle);

// Angles(Theta = angle wrt Z)
const int NBINS_Angle_Beam_Muon = 180;
const double MIN_Angle_Beam_Muon = 0;
const double MAX_Angle_Beam_Muon = 180;
const double WIDTH_Angle_Beam_Muon = getBinWidth(MIN_Angle_Beam_Muon,MAX_Angle_Beam_Muon,NBINS_Angle_Beam_Muon);

// Coplanar Angle
const int NBINS_coplanar = 360;
const double MIN_coplanar = -1.0;
const double MAX_coplanar = 1.0;
const double WIDTH_coplanar = getBinWidth(MIN_coplanar,MAX_coplanar,NBINS_coplanar);

// Proton Score
const int NBINS_proton_score = 40;
const double MIN_proton_score = 0.0;
const double MAX_proton_score = 1.0;
const double WIDTH_proton_score = getBinWidth(MIN_proton_score,MAX_proton_score,NBINS_proton_score);

// Proton Track Length
const int NBINS_proton_length = 200;
const double MIN_proton_length = 0.0;
const double MAX_proton_length = 2000;
const double WIDTH_proton_length = getBinWidth(MIN_proton_length,MAX_proton_length,NBINS_proton_length);


#endif


