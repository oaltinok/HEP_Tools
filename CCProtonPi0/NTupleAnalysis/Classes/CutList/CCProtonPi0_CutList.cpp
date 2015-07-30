/*
    See CutList.h header for Class Information
*/
#ifndef CCProtonPi0_CutList_cpp
#define CCProtonPi0_CutList_cpp

#include "CCProtonPi0_CutList.h"

using namespace PlotUtils;

CCProtonPi0_CutList::CCProtonPi0_CutList(bool isModeReduce, bool isMC) : CCProtonPi0_NTupleAnalysis()
{
    cout<<"Initializing CCProtonPi0_CutList"<<endl;
    
    if(isModeReduce){
        nTrueSignal = 240237;
        SetCutNames();
        OpenOutputFile();
        
        // File Locations
        if (isMC) rootDir = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + "CutHistograms.root";
        else rootDir = Folder_List::rootOut + Folder_List::Data + Folder_List::analyzed + "CutHistograms.root";
        
        cout<<"\tRoot File: "<<rootDir<<endl;
 
        // Create Root File 
        f = new TFile(rootDir.c_str(),"RECREATE");
        
        initHistograms();
    }else{
        cout<<"\tNTuple Analysis Mode -- Will not create ROOT & Text Files"<<endl;
    }

    cout<<"Done!"<<endl;
}

/*
 *  Histograms
 *      See Following Page for Histogram Indices
 *          https://cdcvs.fnal.gov/redmine/projects/minerva/wiki/Ozgur's_scratch_page
 * */
void CCProtonPi0_CutList::initHistograms()
{
    MnvH1D* temp = NULL;

    for (int i = 0; i < nHistograms; i++){
        temp = new MnvH1D( Form("%s_%d","hCut_Michel",i),Form("%d",i),binList.michelID.get_nBins(), binList.michelID.get_min(), binList.michelID.get_max() );
        temp->GetXaxis()->SetTitle("0 = No Michel, 1 = Michel");
        temp->GetYaxis()->SetTitle("N(Events)");
        hCut_Michel.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","hCut_eVis_nuclearTarget",i),"Visible Energy in Nuclear Target",binList.eVis_nuclearTarget.get_nBins(), binList.eVis_nuclearTarget.get_min(), binList.eVis_nuclearTarget.get_max() );
        temp->GetXaxis()->SetTitle("Visible Energy in Nuclear Target [MeV]");
        temp->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.eVis_nuclearTarget.get_width()));
        hCut_eVis_nuclearTarget.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","hCut_eVis_other",i),"Visible Energy in Tracker + ECAL + HCAL",binList.eVis_other.get_nBins(), binList.eVis_other.get_min(), binList.eVis_other.get_max() );
        temp->GetXaxis()->SetTitle("Visible Energy in Tracker + ECAL + HCAL [MeV]");
        temp->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.eVis_other.get_width()));
        hCut_eVis_other.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","hCut_gamma1ConvDist",i),"Leading Photon Conversion Distance",binList.bin_photonConvLength.get_nBins(), binList.bin_photonConvLength.get_min(), binList.bin_photonConvLength.get_max() );
        temp->GetXaxis()->SetTitle("Leading Photon Conversion Distance");
        temp->GetYaxis()->SetTitle(Form("Candidates / %3.2f [cm]",binList.bin_photonConvLength.get_width()));
        hCut_gamma1ConvDist.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","hCut_gamma2ConvDist",i),"Second Photon Conversion Distance",binList.bin_photonConvLength.get_nBins(), binList.bin_photonConvLength.get_min(), binList.bin_photonConvLength.get_max() );
        temp->GetXaxis()->SetTitle("Second Photon Conversion Distance");
        temp->GetYaxis()->SetTitle(Form("Candidates / %3.2f [cm]",binList.bin_photonConvLength.get_width()));
        hCut_gamma2ConvDist.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","hCut_pi0invMass",i),"Reconstructed Pi0 Invariant Mass",binList.pi0_invMass.get_nBins(), binList.pi0_invMass.get_min(), binList.pi0_invMass.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed Pi0 Invariant Mass [MeV]");
        temp->GetYaxis()->SetTitle(Form("Candidates / %3.2f [MeV]",binList.pi0_invMass.get_width()));
        hCut_pi0invMass.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","hCut_1Prong_neutrinoE",i),"Reconstructed Beam Energy",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed Beam Energy [GeV]");
        temp->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.beamE.get_width()));
        hCut_1Prong_neutrinoE.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","hCut_1Prong_UnusedE",i),"Unused Cluster Energy after Pi0 Reconstruction",binList.UnusedE.get_nBins(), binList.UnusedE.get_min(), binList.UnusedE.get_max() );
        temp->GetXaxis()->SetTitle("Unused Cluster Energy after Pi0 Reconstruction [MeV]");
        temp->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.UnusedE.get_width()));
        hCut_1Prong_UnusedE.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","hCut_2Prong_neutrinoE",i),"Reconstructed Beam Energy",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed Beam Energy [GeV]");
        temp->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.beamE.get_width()));
        hCut_2Prong_neutrinoE.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","hCut_2Prong_UnusedE",i),"Unused Cluster Energy after Pi0 Reconstruction",binList.UnusedE.get_nBins(), binList.UnusedE.get_min(), binList.UnusedE.get_max() );
        temp->GetXaxis()->SetTitle("Unused Cluster Energy after Pi0 Reconstruction [MeV]");
        temp->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.UnusedE.get_width()));
        hCut_2Prong_UnusedE.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","hCut_protonScore_pIDDiff",i),"Proton Score - Pion Score",binList.particleScoreDiff.get_nBins(), binList.particleScoreDiff.get_min(), binList.particleScoreDiff.get_max() );
        temp->GetXaxis()->SetTitle("Proton Score - Pion Score");
        temp->GetYaxis()->SetTitle("N(Events)");
        hCut_protonScore_pIDDiff.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","hCut_protonScore_LLR",i),"proton_protonScore_LLR",binList.particleScore_LLR.get_nBins(), binList.particleScore_LLR.get_min(), binList.particleScore_LLR.get_max() );
        temp->GetXaxis()->SetTitle("proton_protonScore_LLR");
        temp->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScore_LLR.get_width()));
        hCut_protonScore_LLR.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","hCut_deltaInvMass",i),"deltaInvMass",binList.deltaInvMass.get_nBins(), binList.deltaInvMass.get_min(), binList.deltaInvMass.get_max() );
        temp->GetXaxis()->SetTitle("hCut_deltaInvMass");
        temp->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.deltaInvMass.get_width()));
        hCut_deltaInvMass.push_back(temp);
    }

    // MC Only Histograms
    mc_w_DIS = new TH1D( "mc_w_DIS","True W for DIS",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max() );
    mc_w_DIS->GetXaxis()->SetTitle("True W for DIS [GeV]");
    mc_w_DIS->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.w.get_width()));
    
    mc_w_RES = new TH1D( "mc_w_RES","True W for RES",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max() );
    mc_w_RES->GetXaxis()->SetTitle("True W for RES [GeV]");
    mc_w_RES->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.w.get_width()));
    
    mc_w_CCQE = new TH1D( "mc_w_CCQE","True W for CCQE",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max() );
    mc_w_CCQE->GetXaxis()->SetTitle("True W for CCQE [GeV]");
    mc_w_CCQE->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.w.get_width()));
}

void CCProtonPi0_CutList::SetCutNames()
{
    // Common Cut Numbers
    nCut_All.set_Name("All");
    nCut_Vertex_None.set_Name("Vertex_None");
    nCut_Vertex_Not_Reconstructable.set_Name("Vertex_Not_Reconstructable"); 
    nCut_Vertex_Not_Fiducial.set_Name("Vertex_Not_Fiducial");
    nCut_Muon_None.set_Name("Muon_None");              
    nCut_Muon_Not_Plausible.set_Name("Muon_Not_Plausible");
    nCut_Muon_Charge.set_Name("Muon_Charge");
    nCut_Vertex_Michel_Exist.set_Name("Vertex_Michel_Exist"); 
    nCut_EndPoint_Michel_Exist.set_Name("EndPoint_Michel_Exist");
    nCut_secEndPoint_Michel_Exist.set_Name("secEndPoint_Michel_Exist");
    nCut_PreFilter_Pi0.set_Name("PreFilter_Pi0");
    nCut_ConeBlobs.set_Name("ConeBlobs");
    nCut_BlobDirectionBad.set_Name("BlobDirectionBad");
    nCut_BlobsBad.set_Name("BlobsBad");
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
    cutFile = Folder_List::output + Folder_List::textOut + "CutTable.txt";
    
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
    cutText.width(12); cutText<<"Eff(AllSignal)"<<" ";      
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

double CCProtonPi0_CutList::getCutEfficiency(CCProtonPi0_Cut& currentCut, double effBase) const
{
    double nSignal_current = currentCut.nSignal.getCount();
    double nSignal_effBase = effBase; 
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
    nCutVector.push_back(nCut_Muon_None);              
    nCutVector.push_back(nCut_Muon_Not_Plausible);
    nCutVector.push_back(nCut_Muon_Charge);
    nCutVector.push_back(nCut_Vertex_Michel_Exist); 
    nCutVector.push_back(nCut_EndPoint_Michel_Exist);
    nCutVector.push_back(nCut_secEndPoint_Michel_Exist);
    nCutVector.push_back(nCut_PreFilter_Pi0);
    nCutVector.push_back(nCut_ConeBlobs);
    nCutVector.push_back(nCut_BlobDirectionBad);
    nCutVector.push_back(nCut_BlobsBad);
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
    double eff_AllSignal;
    double eff_MINOS;
    double purity;    

    eff_AllSignal = getCutEfficiency(currentCut,nTrueSignal);
    eff_MINOS = getCutEfficiency(currentCut,nCut_Muon_None);
    purity = getCutPurity(currentCut);
            
    cutText.width(35); cutText<<currentCut.get_Name()<<" ";
    cutText.width(12); cutText<<currentCut.nEvent.getCount()<<" ";
    
    // Total Signal
    cutText.width(12); cutText<<currentCut.nSignal.getCount()<<" ";

    // Efficiency
    if ( eff_AllSignal <= 100){
        cutText.width(12); cutText<<eff_AllSignal<<" ";
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
    double eff_AllSignal;
    double eff_MINOS;
    double purity;    

    // ------------------------------------------------------------------------
    // Write nCut_1Prong  Statistics
    // ------------------------------------------------------------------------
    eff_AllSignal = getCutEfficiency(nCut_1Prong,nTrueSignal);
    eff_MINOS = getCutEfficiency(nCut_1Prong,nCut_Muon_None);
    purity = getCutPurity(nCut_1Prong);
            
    cutText.unsetf( std::ios::floatfield ); 
    cutText.width(35); cutText<<nCut_1Prong.get_Name()<<" ";
    cutText.width(12); cutText<<nCut_1Prong.nEvent.getCount()<<" ";
    
    // Total Signal
    cutText.width(12); cutText<<nCut_1Prong.nSignal.getCount()<<" ";

    cutText.precision(4); 

    // Efficiency
    if ( eff_AllSignal <= 100){
        cutText.width(12); cutText<<eff_AllSignal<<" ";
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
    eff_AllSignal = getCutEfficiency(nCut_2Prong,nCut_Vertex_Not_Fiducial);
    eff_MINOS = getCutEfficiency(nCut_2Prong,nCut_Muon_None);
    purity = getCutPurity(nCut_2Prong);
            
    cutText.width(12); cutText<<nCut_2Prong.nEvent.getCount()<<" ";
    
    // Total Signal
    cutText.width(12); cutText<<nCut_2Prong.nSignal.getCount()<<" ";

    cutText.precision(4); 

    // Efficiency
    if ( eff_AllSignal <= 100){
        cutText.width(12); cutText<<eff_AllSignal<<" ";
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

void CCProtonPi0_CutList::writeHistograms()
{
    std::cout<<">> Writing "<<rootDir<<std::endl;
    f->cd();
    for (int i = 0; i < nHistograms; i++){
        // Common
        hCut_Michel[i]->Write();
        hCut_eVis_nuclearTarget[i]->Write();
        hCut_eVis_other[i]->Write();
        hCut_pi0invMass[i]->Write();
        hCut_gamma1ConvDist[i]->Write();
        hCut_gamma2ConvDist[i]->Write();
        
        hCut_1Prong_neutrinoE[i]->Write();
        hCut_2Prong_neutrinoE[i]->Write();
        hCut_1Prong_UnusedE[i]->Write();
        hCut_2Prong_UnusedE[i]->Write();

        // 2 Prong Specific
        hCut_protonScore_pIDDiff[i]->Write();
        hCut_protonScore_LLR[i]->Write();
        hCut_deltaInvMass[i]->Write();
    }

    // MC Only
    mc_w_DIS->Write();
    mc_w_RES->Write();
    mc_w_CCQE->Write();

    f->Close();
}

#endif

