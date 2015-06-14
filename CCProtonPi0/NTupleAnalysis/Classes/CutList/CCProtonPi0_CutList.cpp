/*
    See CutList.h header for Class Information
*/
#ifndef CCProtonPi0_CutList_cpp
#define CCProtonPi0_CutList_cpp

#include "CCProtonPi0_CutList.h"

using namespace std;

CCProtonPi0_CutList::CCProtonPi0_CutList(int nMode) : CCProtonPi0_NTupleAnalysis(nMode)
{
    cout<<"Initializing CCProtonPi0_CutList"<<endl;
    
    SetCutNames();
    OpenOutputFile();

    cout<<"Done!"<<endl;
}

void CCProtonPi0_CutList::SetCutNames()
{
    // Common Cut Numbers
    nCut_All.set_Name("All");
    nCut_Vertex_None.set_Name("Vertex_None");
    nCut_Vertex_Not_Reconstructable.set_Name("Vertex_Not_Reconstructable"); 
    nCut_Vertex_Not_Fiducial.set_Name("Vertex_Not_Fiducial");
    nCut_Vertex_Count.set_Name("Vertex_Count");  
    nCut_Muon_None.set_Name("Muon_None");              
    nCut_Muon_Not_Plausible.set_Name("Muon_Not_Plausible");
    nCut_Muon_Score_Low.set_Name("Muon_Score_Low");
    nCut_Muon_Charge.set_Name("Muon_Charge");
    nCut_Vertex_Michel_Exist.set_Name("Vertex_Michel_Exist"); 
    nCut_EndPoint_Michel_Exist.set_Name("EndPoint_Michel_Exist");
    nCut_secEndPoint_Michel_Exist.set_Name("secEndPoint_Michel_Exist");
    nCut_PreFilter_Pi0.set_Name("PreFilter_Pi0");
    nCut_ConeBlobs.set_Name("ConeBlobs");
    nCut_BlobsBad.set_Name("BlobsBad");
    nCut_Pi0BlobCuts.set_Name("Pi0BlobCuts");
    nCut_Photon1DistanceLow.set_Name("Photon1DistanceLow");
    nCut_Photon2DistanceLow.set_Name("Photon2DistanceLow");
    nCut_Pi0_invMass.set_Name("Pi0_invMass");

    // nProngs == 1 Cut Numbers (Muon + Pi0)
    nCut_1Prong_Particle_None.set_Name("Particle_None");
    nCut_1Prong_Proton_None.set_Name("Proton_None");            
    nCut_1Prong_ProtonScore.set_Name("Proton_Score");
    nCut_1Prong_DeltaInvMass.set_Name("Delta_invMass");
    nCut_1Prong_beamEnergy.set_Name("beamEnergy");
    nCut_1Prong_UnusedE.set_Name("UnusedE");

    // nProngs >= 2 Cut Numbers (Muon + Pi0 + X(No Meson))
    nCut_2Prong_Particle_None.set_Name("Particle_None");
    nCut_2Prong_Proton_None.set_Name("Proton_None");            
    nCut_2Prong_ProtonScore.set_Name("Proton_Score");
    nCut_2Prong_DeltaInvMass.set_Name("Delta_invMass");
    nCut_2Prong_beamEnergy.set_Name("beamEnergy");
    nCut_2Prong_UnusedE.set_Name("UnusedE");

}

void CCProtonPi0_CutList::OpenOutputFile()
{
    // Open Cut Files
    cutFile = Folder_List::output + Folder_List::textOut + branchDir + "CutTable.txt";
    
    cutText.open( cutFile.c_str() );
    if( !cutText.is_open() ){
        cerr<<"Cannot open output text file: "<<cutFile<<endl;
        exit(1);
    }else{
        cout<<"\t"<<cutFile<<endl;
    }
     
}

void CCProtonPi0_CutList::writeCutTableHeader()
{
    cutText<<std::left;
    
    cutText.width(35); cutText<<"Cut"<<" "; 
    
    cutText.width(12); cutText<<"N(Events)"<<" ";    
    cutText.width(12); cutText<<"N(Signal)"<<" ";      
    cutText.width(12); cutText<<"Eff(FidVol)"<<" ";      
    cutText.width(12); cutText<<"Eff(MINOS)"<<" ";      
    cutText.width(12); cutText<<"Purity"<<" ";
    cutText.width(12); cutText<<"N(Study1)"<<" "; 
    cutText.width(12); cutText<<"N(Study2)"<<" "; 
    cutText<<endl;
}

double CCProtonPi0_CutList::getCutEfficiency(CCProtonPi0_Cut& currentCut, CCProtonPi0_Cut& effBase) const
{
    double nSignal_current = currentCut.nSignal.getCount();
    double nSignal_effBase = effBase.nSignal.getCount(); 
    double eff = (nSignal_current / nSignal_effBase) * 100.0;
    return eff;   
}

double CCProtonPi0_CutList::getCutPurity(CCProtonPi0_Cut& currentCut) const
{
    double nSignal_current = currentCut.nSignal.getCount();
    double nEvents_current = currentCut.nEvent.getCount();
    double purity  = (nSignal_current / nEvents_current) * 100.0;
    return purity;   
}

void CCProtonPi0_CutList::formCutVector()
{   
    nCutVector.push_back(nCut_All);
    nCutVector.push_back(nCut_Vertex_None);
    nCutVector.push_back(nCut_Vertex_Not_Reconstructable); 
    nCutVector.push_back(nCut_Vertex_Not_Fiducial);
    nCutVector.push_back(nCut_Vertex_Count); 
    nCutVector.push_back(nCut_Muon_None);              
    nCutVector.push_back(nCut_Muon_Not_Plausible);
    nCutVector.push_back(nCut_Muon_Score_Low);
    nCutVector.push_back(nCut_Muon_Charge);
    nCutVector.push_back(nCut_Vertex_Michel_Exist); 
    nCutVector.push_back(nCut_EndPoint_Michel_Exist);
    nCutVector.push_back(nCut_secEndPoint_Michel_Exist);
    nCutVector.push_back(nCut_PreFilter_Pi0);
    nCutVector.push_back(nCut_ConeBlobs);
    nCutVector.push_back(nCut_BlobsBad);
    nCutVector.push_back(nCut_Pi0BlobCuts);
    nCutVector.push_back(nCut_Photon1DistanceLow);
    nCutVector.push_back(nCut_Photon2DistanceLow);
    nCutVector.push_back(nCut_Pi0_invMass);
    
}

void CCProtonPi0_CutList::writeCutTable()
{
    formCutVector();
    
    cout<<">> Writing "<<cutFile<<endl;
    
    writeCutTableHeader();
    writeCutTableRows();
}

void CCProtonPi0_CutList::writeCutTableRows()
{
    // Write General Cuts upto Proton Reconstruction    
    for( unsigned int i = 0; i < nCutVector.size(); i++){
        writeSingleRow(nCutVector[i]);    
    }

    cutText<<endl;

    //Write Separate Cuts;
    writeSingleRow(nCut_1Prong_Particle_None, nCut_2Prong_Particle_None);
    writeSingleRow(nCut_1Prong_Proton_None, nCut_2Prong_Proton_None);
    writeSingleRow(nCut_1Prong_ProtonScore, nCut_2Prong_ProtonScore);
    writeSingleRow(nCut_1Prong_DeltaInvMass, nCut_2Prong_DeltaInvMass);
    writeSingleRow(nCut_1Prong_beamEnergy, nCut_2Prong_beamEnergy);
    writeSingleRow(nCut_1Prong_UnusedE, nCut_2Prong_UnusedE);
}

void CCProtonPi0_CutList::writeSingleRow(CCProtonPi0_Cut& currentCut)
{
    double eff_FidVolume;
    double eff_MINOS;
    double purity;    

    eff_FidVolume = getCutEfficiency(currentCut,nCut_Vertex_Not_Fiducial);
    eff_MINOS = getCutEfficiency(currentCut,nCut_Muon_None);
    purity = getCutPurity(currentCut);
            
    cutText.width(35); cutText<<currentCut.get_Name()<<" ";
    cutText.width(12); cutText<<currentCut.nEvent.getCount()<<" ";
    
    // Total Signal
    cutText.width(12); cutText<<currentCut.nSignal.getCount()<<" ";

    // Efficiency
    if ( eff_FidVolume <= 100){
        cutText.width(12); cutText<<eff_FidVolume<<" ";
    }else{
        cutText.width(12); cutText<<"N/A"<<" ";    
    }

    if ( eff_MINOS <= 100){
        cutText.width(12); cutText<<eff_MINOS<<" ";
    }else{
        cutText.width(12); cutText<<"N/A"<<" ";    
    }    
    
    // Purity
    cutText.width(12); cutText<<purity<<" ";

    // Number of Events which are studied
    cutText.width(12); cutText<<currentCut.nStudy1.getCount()<<" ";
    cutText.width(12); cutText<<currentCut.nStudy2.getCount()<<" ";
    
    cutText<<endl;
}



void CCProtonPi0_CutList::writeSingleRow(CCProtonPi0_Cut& nCut_1Prong, CCProtonPi0_Cut& nCut_2Prong)
{
    double eff_FidVolume;
    double eff_MINOS;
    double purity;    

    // ------------------------------------------------------------------------
    // Write nCut_1Prong  Statistics
    // ------------------------------------------------------------------------
    eff_FidVolume = getCutEfficiency(nCut_1Prong,nCut_Vertex_Not_Fiducial);
    eff_MINOS = getCutEfficiency(nCut_1Prong,nCut_Muon_None);
    purity = getCutPurity(nCut_1Prong);
            
    cutText.unsetf( std::ios::floatfield ); 
    cutText.width(35); cutText<<nCut_1Prong.get_Name()<<" ";
    cutText.width(12); cutText<<nCut_1Prong.nEvent.getCount()<<" ";
    
    // Total Signal
    cutText.width(12); cutText<<nCut_1Prong.nSignal.getCount()<<" ";

    cutText.precision(4); 

    // Efficiency
    if ( eff_FidVolume <= 100){
        cutText.width(12); cutText<<eff_FidVolume<<" ";
    }else{
        cutText.width(12); cutText<<"N/A"<<" ";    
    }

    if ( eff_MINOS <= 100){
        cutText.width(12); cutText<<eff_MINOS<<" ";
    }else{
        cutText.width(12); cutText<<"N/A"<<" ";    
    }    
    
    // Purity
    cutText.width(12); cutText<<purity<<" ";

    // Number of Events which are studied
    cutText.width(12); cutText<<nCut_1Prong.nStudy1.getCount()<<" ";
    cutText.width(12); cutText<<nCut_1Prong.nStudy2.getCount()<<" ";

    // ------------------------------------------------------------------------
    // Write nCut_2Prong Statistics
    // ------------------------------------------------------------------------
    eff_FidVolume = getCutEfficiency(nCut_2Prong,nCut_Vertex_Not_Fiducial);
    eff_MINOS = getCutEfficiency(nCut_2Prong,nCut_Muon_None);
    purity = getCutPurity(nCut_2Prong);
            
    cutText.width(12); cutText<<nCut_2Prong.nEvent.getCount()<<" ";
    
    // Total Signal
    cutText.width(12); cutText<<nCut_2Prong.nSignal.getCount()<<" ";

    cutText.precision(4); 

    // Efficiency
    if ( eff_FidVolume <= 100){
        cutText.width(12); cutText<<eff_FidVolume<<" ";
    }else{
        cutText.width(12); cutText<<"N/A"<<" ";    
    }

    if ( eff_MINOS <= 100){
        cutText.width(12); cutText<<eff_MINOS<<" ";
    }else{
        cutText.width(12); cutText<<"N/A"<<" ";    
    }    
    
    // Purity
    cutText.width(12); cutText<<purity<<" ";

    // Number of Events which are studied
    cutText.width(12); cutText<<nCut_2Prong.nStudy1.getCount()<<" ";
    cutText.width(12); cutText<<nCut_2Prong.nStudy2.getCount()<<" ";

    cutText<<endl;
}

CCProtonPi0_CutList::~CCProtonPi0_CutList()
{
    cutText.close(); 
}


#endif
