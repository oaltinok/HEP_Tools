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


