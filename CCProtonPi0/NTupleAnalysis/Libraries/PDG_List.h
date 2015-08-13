/*
================================================================================
Namespace: PDG_List
    Particle Data Group Particle Identification List
    List of Commonly used particles PDGs
    
    Please contact ozgur.altinok@tufts.edu if you find a mistake in the list
    
    Main Directory:
        Libraries/
    
    Usage:
        > #include "Libraries/PDG_List.h" 
        > PDG_List::proton      // Returns proton PDG (2212)
            
    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
================================================================================
*/

#ifndef PDG_List_h
#define PDG_List_h

namespace PDG_List
{

const int gamma = 22;

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

// -----------------------------------------------------------------------------
//      Baryon Resonance
//------------------------------------------------------------------------------
const double Delta_p_1232 = 2214;
const double Delta_pp_1232 = 2224;

const double Delta_p_1620 = 2122;
const double Delta_pp_1620 = 2222;

const double Delta_p_1700 = 12214;
const double Delta_pp_1700 = 12224;

const double N_p_1440 = 12212;
const double N_p_1520 = 2124;
const double N_p_1535 = 22212;
const double N_p_1650 = 32212;
const double N_p_1675 = 2216;
const double N_p_1680 = 12216;
}

#endif

