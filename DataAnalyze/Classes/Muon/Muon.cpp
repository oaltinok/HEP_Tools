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

void Muon::set_angleMuon(TVector3 muonp3)
{
    angleMuon[1] = 0;
}




