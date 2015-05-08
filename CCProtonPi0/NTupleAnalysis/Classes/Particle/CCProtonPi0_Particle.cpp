/*
    See CCProtonPi0_Particle.h header for Class Information
*/

#ifndef CCProtonPi0_Particle_cpp
#define CCProtonPi0_Particle_cpp

#include "CCProtonPi0_Particle.h"

using namespace std;


CCProtonPi0_Particle::CCProtonPi0_Particle(int nMode) : CCProtonPi0_NTupleAnalysis(nMode)
{
    bin_partScore.setBin(100,0.0,1.0);
    bin_error.setBin(400, -2.0, 2.0);
    bin_angle.setBin(180, 0.0, 180.0);
}

void CCProtonPi0_Particle::fill_Histograms()
{
    partScore->Fill(particleScore);
    
    E_reco->Fill(energy[0]);
    E_mc->Fill(energy[1]);
    E_error->Fill(energy[2]);
    E_reco_mc->Fill(energy[0],energy[1]);
   
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

void CCProtonPi0_Particle::write_RootFile()
{
    cout<<">> Writing "<<rootDir<<endl;
    f->Write();
}

void CCProtonPi0_Particle::set_angleMuon(CCProtonPi0_Particle &mu, bool isMC)
{
    int type = getDataType(isMC);
    TVector3 muonp3 = mu.p4[type].Vect(); // Get Muon 3-Momentum
    
    angleMuon[type] = p4[type].Angle( muonp3 );
}

void CCProtonPi0_Particle::set_kineticEnergy(bool isMC)
{
    int type = getDataType(isMC);
    cout<<"Warning using CCProtonPi0_Particle::set_kineticEnergy!"<<endl;
    kineticEnergy[type] = -1.0;
}

void CCProtonPi0_Particle::set_angleBeam(TVector3 beamp3, bool isMC)
{
    int type = getDataType(isMC);
    angleBeam[type] = p4[type].Angle(beamp3);
}

void CCProtonPi0_Particle::set_p4(double px, double py, double pz, double E, bool isMC)
{
    int type = getDataType(isMC);
    
    p4[type].SetPxPyPzE(px,py,pz,E);    // Set 4-Momentum
    momentum[type] = p4[type].P();      // Set 3-Momentum
    energy[type] = p4[type][3];         // Set Total Energy
    set_kineticEnergy(isMC);            // Set Kinetic Energy
}

void CCProtonPi0_Particle::set_errors()
{
    energy[2] = Data_Functions::getError(energy[1], energy[0]);
    momentum[2] = Data_Functions::getError(momentum[1], momentum[0]);
    kineticEnergy[2] = Data_Functions::getError(kineticEnergy[1], kineticEnergy[0]);
    angleMuon[2] = Data_Functions::getError(angleMuon[1] * TMath::RadToDeg(), angleMuon[0] * TMath::RadToDeg());
    angleBeam[2] = Data_Functions::getError(angleBeam[1] * TMath::RadToDeg(), angleBeam[0] * TMath::RadToDeg());
}

int CCProtonPi0_Particle::getDataType(bool isMC)
{
    int type; // Indice for data type ( reco = 0, mc = 1 )
    
    if ( isMC ){
        type = 1;
    }else{
        type = 0;
    }
    
    return type;

}

CCProtonPi0_Particle::~CCProtonPi0_Particle()
{    
    delete partScore;
    
    delete E_mc;
    delete E_reco;
    delete E_error;
    delete E_reco_mc;
    
    delete P_mc;
    delete P_reco;
    delete P_error;
    delete P_reco_mc;
    
    delete KE_mc;
    delete KE_reco;
    delete KE_error;
    delete KE_reco_mc;
    
    delete angleBeam_mc;
    delete angleBeam_reco;
    delete angleBeam_error;
    delete angleBeam_reco_mc;
    
    delete angleMuon_mc;
    delete angleMuon_reco;
    delete angleMuon_error;
    delete angleMuon_reco_mc;   
    
    delete f;
}
#endif


