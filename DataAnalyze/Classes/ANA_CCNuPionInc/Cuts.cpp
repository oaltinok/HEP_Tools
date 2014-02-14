/*
================================================================================
Function List: Cuts.cpp
    Contains all Data Selection(Cut) Functions
    
    Functions Defined under ANA_CCNuPionInc Class

    Usage:
        > #include "ANA_CCNuPionInc.cpp"
        > #ifdef ANA_CCNuPionInc_cxx 
    
    Last Revision: 2014_01_28
================================================================================
*/

#include "ANA_CCNuPionInc.cpp"
#ifdef ANA_CCNuPionInc_cxx

/*
--------------------------------------------------------------------------------
 Volume Cut: isVertexContained()
    Vertex of the interaction must be inside Detector Interaction Region
--------------------------------------------------------------------------------
*/
bool ANA_CCNuPionInc::isVertexContained()
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
 findParticle:
    Returns the array indice of the Final State Particle for given PDG
--------------------------------------------------------------------------------
*/
int ANA_CCNuPionInc::findParticle(int targetPDG)
{
    int ind = -1;
    
    for(int i = 0; i < mc_nFSPart; i++){
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
int ANA_CCNuPionInc::findProton()
{
    int currentProtonind = -1;
    double currentProtonE = -1;
    
    for(int i = 0; i < mc_nFSPart; i++ ){
    
        // Energy of Proton (Only the fastest proton)
        // Update Proton Energy until you find the highest energy
        if( mc_FSPartPDG[i] == PDG_List::proton){
            if (currentProtonE < mc_FSPartE[i]){
                currentProtonE = mc_FSPartE[i];
                currentProtonind = i;
            }
        }
    }
    
    if (currentProtonE > 0){
        return currentProtonind;
    }else{
        return -1;
    }

}

/*
--------------------------------------------------------------------------------
 findPion:
    Returns the single pi_plus indice
--------------------------------------------------------------------------------
*/
int ANA_CCNuPionInc::findPion()
{
  
    // Return a valid indice only if there is only single pi+ in Final State
    if(isSinglePion()){
        return findParticle(PDG_List::pi_plus);
    }else{
        return -1;
    }
}

bool ANA_CCNuPionInc::isSinglePion()
{
    int nPlus;
    int nMinus;
    int nZero;
    
    nPlus = countParticles(PDG_List::pi_plus,true);
    nMinus = countParticles(PDG_List::pi_minus,false);
    nZero = countParticles(PDG_List::pi_zero,false);
    
    if(nMinus == 0 && nZero == 0 && nPlus == 1){
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
int ANA_CCNuPionInc::countParticles(int targetPDG, bool applyPCut)
{
    int count = 0;
    TVector3 p3;
    
    for(int i = 0; i < mc_nFSPart; i++ ){
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
