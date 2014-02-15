/*
    See Particle.h header for Class Information
*/
#include "Particle.h"

using namespace std;


Particle::Particle()
{
    p4 = new TLorentzVector[N_CHANNELS];
    angleBeam = new double[N_CHANNELS];
    angleMuon = new double[N_CHANNELS];
}

virtual void Particle::set_angleMuon(TVector3 muonp3)
{
    angleMuon[1] = p4[1].Angle(muonp3);

}

void Particle::set_angleBeam(TVector3 beamp3)
{
    angleBeam[1] = p4[1].Angle(beamp3);
}


