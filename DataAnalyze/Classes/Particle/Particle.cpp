/*
    See Particle.h header for Class Information
*/
#include "Particle.h"

using namespace std;

Particle::Particle()
{
    bin_error.setBin(400, -2.0, 2.0);
    bin_angle.setBin(180, 0.0, 180.0);
}

void Particle::fill_Histograms()
{
    P_reco->Fill(momentum[0]);
    P_mc->Fill(momentum[1]);
    P_error->Fill(momentum[2]);
    P_reco_mc->Fill(momentum[0],momentum[1]);
    
    KE_reco->Fill(kineticEnergy[0]);
    KE_mc->Fill(kineticEnergy[1]);
    KE_error->Fill(kineticEnergy[2]);
    KE_reco_mc->Fill(kineticEnergy[0],kineticEnergy[1]);
    
    angleMuon_reco->Fill(angleMuon[0] * TMath::RadToDeg());
    angleMuon_mc->Fill(angleMuon[1] * TMath::RadToDeg());
    angleMuon_error->Fill(angleMuon[2]);
    angleMuon_reco_mc->Fill(angleMuon[0] * TMath::RadToDeg(),angleMuon[1] * TMath::RadToDeg());
    
    angleBeam_reco->Fill(angleBeam[0] * TMath::RadToDeg());
    angleBeam_mc->Fill(angleBeam[1] * TMath::RadToDeg());
    angleBeam_error->Fill(angleBeam[2]);
    angleBeam_reco_mc->Fill(angleBeam[0] * TMath::RadToDeg(),angleBeam[1] * TMath::RadToDeg());

}

void Particle::write_RootFile()
{
    cout<<">> Writing "<<rootDir<<endl;
    f->Write();
}

virtual void Particle::set_angleMuon(Particle &mu, bool isMC)
{
    int type = getDataType(isMC);
    TVector3 muonp3 = mu.p4[type].Vect(); // Get Muon 3-Momentum
    
    angleMuon[type] = p4[type].Angle( muonp3 );
}

void Particle::set_angleBeam(TVector3 beamp3, bool isMC)
{
    int type = getDataType(isMC);
    angleBeam[type] = p4[type].Angle(beamp3);
}

void Particle::set_p4(double px, double py, double pz, double E, bool isMC)
{
    int type = getDataType(isMC);
    
    p4[type].SetPxPyPzE(px,py,pz,E);    // Set 4-Momentum
    set_momentum(isMC);                 // Set 3-Momentum
    set_kineticEnergy(isMC);            // Set Kinetic Energy
}

void Particle::set_errors()
{
    momentum[2] = calc_error(momentum[1], momentum[0]);
    kineticEnergy[2] = calc_error(kineticEnergy[1], kineticEnergy[0]);
    angleMuon[2] = calc_error(angleMuon[1] * TMath::RadToDeg(), angleMuon[0] * TMath::RadToDeg());
    angleBeam[2] = calc_error(angleBeam[1] * TMath::RadToDeg(), angleBeam[0] * TMath::RadToDeg());
}

void Particle::set_momentum(bool isMC)
{
    int type = getDataType(isMC);
    
    momentum[type] = p4[type].P();
}


double Particle::calc_error(double trueValue, double recoValue)
{
    double error;
    
    error = (trueValue - recoValue) / trueValue;
    
    return error;
}

int Particle::getDataType(bool isMC)
{
    int type; // Indice for data type ( reco = 0, mc = 1 )
    
    if ( isMC ){
        type = 1;
    }else{
        type = 0;
    }
    
    return type;

}

Particle::~Particle()
{

}


