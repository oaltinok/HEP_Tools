/*
    See Muon.h header for Class Information
*/
#include "Muon.h"

using namespace std;

Muon::Muon()
{
    
}

bool Muon::get_isMinosMatched()
{
    return isMinosMatched;
}

void Muon::set_isMinosMatched(bool input)
{
    isMinosMatched = input;
}

void Muon::set_angleMuon(Particle &mu, bool isMC)
{
    // There is only 1 muon in the interaction
    // Do not need to calculate the angle wrt itself
    // Set the angle to zero
    for( int i = 0; i < N_DATA_TYPE; i++){
        angleMuon[i] = 0.0;
    }
}

void Muon::set_kineticEnergy(bool isMC)
{
    int type = getDataType(isMC);
    
    kineticEnergy[type] = p4[type].Energy() - restMass;
}




