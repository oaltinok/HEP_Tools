/*
================================================================================
Function List: Cuts.cpp
    Contains all Data Selection(Cut) Functions
    
    Functions Defined under CCPi0 Class

    Usage:
        > #include "CCPi0.cpp"
        > #ifdef CCPi0_cxx 
    
    Last Revision: 2014_03_17
================================================================================
*/

#include "CCPi0.cpp"
#ifdef CCPi0_cxx

/*
--------------------------------------------------------------------------------
 Volume Cut: isVertexContained()
    Vertex of the interaction must be inside Detector Interaction Region
--------------------------------------------------------------------------------
*/
bool CCPi0::isVertexContained()
{
    double rMax = 600*sqrt(2);
    double zMin = 5900;
    double zMax = 8000;
    double r;
    
    if(mc_vtx[2] > zMin && mc_vtx[2] < zMax){
        r = sqrt(mc_vtx[0]*mc_vtx[0] + mc_vtx[1]*mc_vtx[1]);
        if(r < rMax){
            return true;
        }
    }
    
    return false;
}

/*
--------------------------------------------------------------------------------
 MinosCut: isMinosMatched
    Returns True if the muon matches with MINOS detector
    Returns False if the muon does NOT matches with MINOS detector
--------------------------------------------------------------------------------
*/
bool CCPi0::isMinosMatched()
{
    if(minos_track_match){
        return true;
    }else{
        return false;
    }

}


/*
--------------------------------------------------------------------------------
 Volume Cut: isisBeamEnergyLow(double maxEnergy)
    Incoming Beam Energy must be lower than a maximum Energy
--------------------------------------------------------------------------------
*/
bool CCPi0::isBeamEnergyLow(double maxEnergy)
{
    if(mc_incomingE > maxEnergy){
        return false;
    }else{
        return true;
    }

}

/*
--------------------------------------------------------------------------------
 findParticle:
    Returns the array indice of the Final State Particle for given PDG
--------------------------------------------------------------------------------
*/
int CCPi0::findParticle(int targetPDG)
{
    int ind = -1;
    
    for(int i = 0; i < mc_nFSPart && 1 < max_nFSPart; i++){
        if(  mc_FSPartPDG[i] == targetPDG ){
            ind = i;
            break;
        }
    }

    return ind;

}


/*
--------------------------------------------------------------------------------
 findProton:
    Returns the interaction proton indice
--------------------------------------------------------------------------------
*/
int CCPi0::findProton()
{
    int currentProtonind = -1;
    double currentProtonP = 0;
    
    TVector3 p3;
    
    for(int i = 0; i < mc_nFSPart && 1 < max_nFSPart; i++ ){
    
        // Momentum of Proton (No Proton at rest)
        // Update Proton Energy until you find the highest energy
        if( mc_FSPartPDG[i] == PDG_List::proton){
            p3.SetXYZ(mc_FSPartPx[i],mc_FSPartPy[i],mc_FSPartPz[i]);
            if (currentProtonP < p3.Mag()){
                currentProtonP = p3.Mag();
                currentProtonind = i;
            }
        }
    }
    
    if (currentProtonP > 0){
        return currentProtonind;
    }else{
        return -1;
    }

}

bool CCPi0::isProtonShort(int ind)
{
    const double protonMass = 938; 
    const double minProtonKE = 120;
    double protonKE;
    double protonE;
    
    protonE = mc_FSPartE[ind];
    protonKE =  protonE - protonMass;
    if( protonKE < minProtonKE){
        return true;
    }else{
        return false;
    }
}

/*
--------------------------------------------------------------------------------
 findPion:
    Returns the single pi_plus indice
--------------------------------------------------------------------------------
*/
int CCPi0::findPion()
{
  
    // Return a valid indice only if there is only single pi+ in Final State
    if(isSinglePion()){
        return findParticle(PDG_List::pi_zero);
    }else{
        return -1;
    }
}

bool CCPi0::isSinglePion()
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

bool CCPi0::isNoMeson()
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

/*
--------------------------------------------------------------------------------
 countParticles:
    Returns the number of particles in the Final State
        Input 
            int targetPDG
            bool applyPCut - Variable for selecting particles with momentum 
                                (no particle at rest)
--------------------------------------------------------------------------------
*/
int CCPi0::countParticles(int targetPDG, bool applyPCut)
{
    int count = 0;
    TVector3 p3;
    
    for(int i = 0; i < mc_nFSPart && 1 < max_nFSPart; i++ ){
        if( mc_FSPartPDG[i] == targetPDG){
            if(applyPCut){
                p3.SetXYZ(mc_FSPartPx[i],mc_FSPartPy[i],mc_FSPartPz[i]);
                if(p3.Mag() > 0){
                    count++;
                }
            }
            else{
                count++;
            }
        }
    }
    
    return count;

}


#endif
