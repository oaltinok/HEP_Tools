/*
================================================================================
Function List: Pion_Functions.cpp
    Contains Pion Specific Functions
    
    Functions Defined under CCDeltaPlus Namespace

    Usage:
        > #include "CCDeltaPlus.cpp"
        > #ifdef CCDeltaPlus_cxx 
    
    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision:  2014_05_01
================================================================================
*/

#include "CCDeltaPlus.cpp"
#ifdef CCDeltaPlus_cxx

void CCDeltaPlus::fillPionReco()
{
    // Fill 4-Momentum
    pion.set_p4(    pimom[0],
                    pimom[1],
                    pimom[2],
                    pienergy,
                    false);
    
    // set Angle wrt Beam
    pion.set_angleBeam(beam_p3, false);
    
    // set Angle wrt Muon
    pion.set_angleMuon(muon, false);
}

void CCDeltaPlus::fillPionTrue()
{
    
}

/*
--------------------------------------------------------------------------------
 findPion:
    Returns the single pi_plus indice
--------------------------------------------------------------------------------
*/
int CCDeltaPlus::findPion()
{
  
    // Return a valid indice only if there is only single pi+ in Final State
    if(isSinglePion()){
        return findParticle(PDG_List::pi_zero);
    }else{
        return -1;
    }
}

bool CCDeltaPlus::isSinglePion()
{
    int nPlus;
    int nMinus;
    int nZero;
    
    nPlus = countParticles(PDG_List::pi_plus,false);
    nMinus = countParticles(PDG_List::pi_minus,false);
    nZero = countParticles(PDG_List::pi_zero,true);
    
    if(nMinus == 0 && nZero == 1 && nPlus == 0){
        return true;
    }else{
        return false;
    }

}

bool CCDeltaPlus::isNoMeson()
{
    int nPlus;
    int nMinus;
    int nZero;
    int nKaons;
    int nGamma;
    int nPions;

    
    nPlus = countParticles(PDG_List::pi_plus,false);
    nMinus = countParticles(PDG_List::pi_minus,false);
    nZero = countParticles(PDG_List::pi_zero,false);
    nKaons = countParticles(PDG_List::kaon_zero_L,false);
    nKaons = nKaons +  countParticles(PDG_List::kaon_zero_S,false);
    nKaons = nKaons +  countParticles(PDG_List::kaon_zero,false);
    nKaons = nKaons +  countParticles(PDG_List::kaon_plus,false);
    nKaons = nKaons +  countParticles(PDG_List::kaon_minus,false);
    nGamma = countParticles(PDG_List::gamma,false);
    
    nPions = nMinus + nZero + nPlus;
    
    if(nPions == 0 && nKaons == 0 & nGamma == 0 ){
        return true;
    }else{
        return false;
    }

}


#endif // #ifdef CCDeltaPlus_cxx