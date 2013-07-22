#ifndef PDG_List.h
#define PDG_List.h

using namespace std;

// -----------------------------------------------------------------------------
//      Leptons
//------------------------------------------------------------------------------
const int electron = 11;
const int positron = -electron;

const int mu_minus = 13;
const int mu_plus = -mu_minus;

const int tau_minus = 15;
const int tau_plus = -tau_minus;

const int neutrino_electron = 12;
const int neutrino_muon = 14;
const int neutrino_tau = 16;


// -----------------------------------------------------------------------------
//      Light Baryons
//------------------------------------------------------------------------------
const int proton = 2212;
const int neutron = 2112;


// -----------------------------------------------------------------------------
//      Strange Mesons
//------------------------------------------------------------------------------
const int kaon_zero_L = 130;
const int kaon_zero_S = 310;
const int kaon_zero = 311;
const int kaon_plus = 321;
const int kaon_minus = -kaon_plus;



// -----------------------------------------------------------------------------
//      Light Mesons
//------------------------------------------------------------------------------
const int pi_zero = 111;
const int pi_plus = 211;
const int pi_minus = -pi_plus;


#endif

