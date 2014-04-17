/*
================================================================================
Function List: Cuts.cpp
    Contains all Data Selection(Cut) Functions
    
    Functions Defined under CCDeltaPlus Class

    Usage:
        > #include "CCDeltaPlus.cpp"
        > #ifdef CCDeltaPlus_cxx 
    
    Last Revision: 2014_04_17
================================================================================
*/

#include "CCDeltaPlus.cpp"
#ifdef CCDeltaPlus_cxx

/*
--------------------------------------------------------------------------------
 Volume Cut: isVertexContained()
    Vertex of the interaction must be inside Detector Interaction Region
--------------------------------------------------------------------------------
*/
bool CCDeltaPlus::isVertexContained()
{

    const double fiducialApothem        =  850.0; 
    const double fiducialUpstreamZ      = 5990.0; 
    const double fiducialDownstreamZ    = 8340.0; 
    
    double r;
    
    if(mc_vtx[2] > fiducialUpstreamZ && mc_vtx[2] < fiducialDownstreamZ ){
        r = sqrt(mc_vtx[0]*mc_vtx[0] + mc_vtx[1]*mc_vtx[1]);
        if(r < fiducialApothem){
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
// bool CCDeltaPlus::isMinosMatched()
// {
//     if(minos_track_match){
//         return true;
//     }else{
//         return false;
//     }
// 
// }


/*
--------------------------------------------------------------------------------
 Volume Cut: isisBeamEnergyLow(double maxEnergy)
    Incoming Beam Energy must be lower than a maximum Energy
--------------------------------------------------------------------------------
*/
bool CCDeltaPlus::isBeamEnergyLow(double maxEnergy)
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
int CCDeltaPlus::findParticle(int targetPDG)
{
    int ind = -1;
    
    for(int i = 0; i < mc_nFSPart && i < max_nFSPart; i++){
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
int CCDeltaPlus::findProton()
{
    int currentProtonind = -1;
    double currentProtonP = 0;
    
    TVector3 p3;
    
    for(int i = 0; i < mc_nFSPart && i < max_nFSPart; i++ ){
    
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

bool CCDeltaPlus::isProtonShort(int ind)
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
int CCDeltaPlus::countParticles(int targetPDG, bool applyPCut)
{
    int count = 0;
    TVector3 p3;
    
    for(int i = 0; i < mc_nFSPart && i < max_nFSPart; i++ ){
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
