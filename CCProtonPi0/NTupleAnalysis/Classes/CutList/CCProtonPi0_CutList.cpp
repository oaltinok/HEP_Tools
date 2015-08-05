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
        init_nCutVectors();
        SetCutNames();
        OpenTextFiles(isMC);
        
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
        // --------------------------------------------------------------------
        // 1 Track
        // --------------------------------------------------------------------
        temp = new MnvH1D( Form("%s_%d","hCut_1Track_Michel",i),Form("%d",i),binList.michelID.get_nBins(), binList.michelID.get_min(), binList.michelID.get_max() );
        temp->GetXaxis()->SetTitle("0 = No Michel, 1 = Michel");
        temp->GetYaxis()->SetTitle("N(Events)");
        hCut_1Track_Michel.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","hCut_1Track_eVis_nuclearTarget",i),"Visible Energy in Nuclear Target",binList.eVis_nuclearTarget.get_nBins(), binList.eVis_nuclearTarget.get_min(), binList.eVis_nuclearTarget.get_max() );
        temp->GetXaxis()->SetTitle("Visible Energy in Nuclear Target [MeV]");
        temp->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.eVis_nuclearTarget.get_width()));
        hCut_1Track_eVis_nuclearTarget.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","hCut_1Track_eVis_other",i),"Visible Energy in Tracker + ECAL + HCAL",binList.eVis_other.get_nBins(), binList.eVis_other.get_min(), binList.eVis_other.get_max() );
        temp->GetXaxis()->SetTitle("Visible Energy in Tracker + ECAL + HCAL [MeV]");
        temp->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.eVis_other.get_width()));
        hCut_1Track_eVis_other.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","hCut_1Track_gamma1ConvDist",i),"Leading Photon Conversion Distance",binList.bin_photonConvLength.get_nBins(), binList.bin_photonConvLength.get_min(), binList.bin_photonConvLength.get_max() );
        temp->GetXaxis()->SetTitle("Leading Photon Conversion Distance");
        temp->GetYaxis()->SetTitle(Form("Candidates / %3.2f [cm]",binList.bin_photonConvLength.get_width()));
        hCut_1Track_gamma1ConvDist.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","hCut_1Track_gamma2ConvDist",i),"Second Photon Conversion Distance",binList.bin_photonConvLength.get_nBins(), binList.bin_photonConvLength.get_min(), binList.bin_photonConvLength.get_max() );
        temp->GetXaxis()->SetTitle("Second Photon Conversion Distance");
        temp->GetYaxis()->SetTitle(Form("Candidates / %3.2f [cm]",binList.bin_photonConvLength.get_width()));
        hCut_1Track_gamma2ConvDist.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","hCut_1Track_pi0invMass",i),"Reconstructed Pi0 Invariant Mass",binList.pi0_invMass.get_nBins(), binList.pi0_invMass.get_min(), binList.pi0_invMass.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed Pi0 Invariant Mass [MeV]");
        temp->GetYaxis()->SetTitle(Form("Candidates / %3.2f [MeV]",binList.pi0_invMass.get_width()));
        hCut_1Track_pi0invMass.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","hCut_1Track_neutrinoE",i),"Reconstructed Beam Energy",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed Beam Energy [GeV]");
        temp->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.beamE.get_width()));
        hCut_1Track_neutrinoE.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","hCut_1Track_UnusedE",i),"Unused Cluster Energy after Pi0 Reconstruction",binList.UnusedE.get_nBins(), binList.UnusedE.get_min(), binList.UnusedE.get_max() );
        temp->GetXaxis()->SetTitle("Unused Cluster Energy after Pi0 Reconstruction [MeV]");
        temp->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.UnusedE.get_width()));
        hCut_1Track_UnusedE.push_back(temp);

        // --------------------------------------------------------------------
        // 2 Track
        // --------------------------------------------------------------------
        temp = new MnvH1D( Form("%s_%d","hCut_2Track_Michel",i),Form("%d",i),binList.michelID.get_nBins(), binList.michelID.get_min(), binList.michelID.get_max() );
        temp->GetXaxis()->SetTitle("0 = No Michel, 1 = Michel");
        temp->GetYaxis()->SetTitle("N(Events)");
        hCut_2Track_Michel.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","hCut_2Track_eVis_nuclearTarget",i),"Visible Energy in Nuclear Target",binList.eVis_nuclearTarget.get_nBins(), binList.eVis_nuclearTarget.get_min(), binList.eVis_nuclearTarget.get_max() );
        temp->GetXaxis()->SetTitle("Visible Energy in Nuclear Target [MeV]");
        temp->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.eVis_nuclearTarget.get_width()));
        hCut_2Track_eVis_nuclearTarget.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","hCut_2Track_eVis_other",i),"Visible Energy in Tracker + ECAL + HCAL",binList.eVis_other.get_nBins(), binList.eVis_other.get_min(), binList.eVis_other.get_max() );
        temp->GetXaxis()->SetTitle("Visible Energy in Tracker + ECAL + HCAL [MeV]");
        temp->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.eVis_other.get_width()));
        hCut_2Track_eVis_other.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","hCut_2Track_gamma1ConvDist",i),"Leading Photon Conversion Distance",binList.bin_photonConvLength.get_nBins(), binList.bin_photonConvLength.get_min(), binList.bin_photonConvLength.get_max() );
        temp->GetXaxis()->SetTitle("Leading Photon Conversion Distance");
        temp->GetYaxis()->SetTitle(Form("Candidates / %3.2f [cm]",binList.bin_photonConvLength.get_width()));
        hCut_2Track_gamma1ConvDist.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","hCut_2Track_gamma2ConvDist",i),"Second Photon Conversion Distance",binList.bin_photonConvLength.get_nBins(), binList.bin_photonConvLength.get_min(), binList.bin_photonConvLength.get_max() );
        temp->GetXaxis()->SetTitle("Second Photon Conversion Distance");
        temp->GetYaxis()->SetTitle(Form("Candidates / %3.2f [cm]",binList.bin_photonConvLength.get_width()));
        hCut_2Track_gamma2ConvDist.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","hCut_2Track_pi0invMass",i),"Reconstructed Pi0 Invariant Mass",binList.pi0_invMass.get_nBins(), binList.pi0_invMass.get_min(), binList.pi0_invMass.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed Pi0 Invariant Mass [MeV]");
        temp->GetYaxis()->SetTitle(Form("Candidates / %3.2f [MeV]",binList.pi0_invMass.get_width()));
        hCut_2Track_pi0invMass.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","hCut_2Track_neutrinoE",i),"Reconstructed Beam Energy",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed Beam Energy [GeV]");
        temp->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.beamE.get_width()));
        hCut_2Track_neutrinoE.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","hCut_2Track_UnusedE",i),"Unused Cluster Energy after Pi0 Reconstruction",binList.UnusedE.get_nBins(), binList.UnusedE.get_min(), binList.UnusedE.get_max() );
        temp->GetXaxis()->SetTitle("Unused Cluster Energy after Pi0 Reconstruction [MeV]");
        temp->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.UnusedE.get_width()));
        hCut_2Track_UnusedE.push_back(temp);
        
        temp = new MnvH1D( Form("%s_%d","hCut_2Track_protonScore_pIDDiff",i),"Proton Score - Pion Score",binList.particleScoreDiff.get_nBins(), binList.particleScoreDiff.get_min(), binList.particleScoreDiff.get_max() );
        temp->GetXaxis()->SetTitle("Proton Score - Pion Score");
        temp->GetYaxis()->SetTitle("N(Events)");
        hCut_2Track_protonScore_pIDDiff.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","hCut_2Track_protonScore_LLR",i),"proton_protonScore_LLR",binList.particleScore_LLR.get_nBins(), binList.particleScore_LLR.get_min(), binList.particleScore_LLR.get_max() );
        temp->GetXaxis()->SetTitle("proton_protonScore_LLR");
        temp->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.particleScore_LLR.get_width()));
        hCut_2Track_protonScore_LLR.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","hCut_2Track_deltaInvMass",i),"deltaInvMass",binList.deltaInvMass.get_nBins(), binList.deltaInvMass.get_min(), binList.deltaInvMass.get_max() );
        temp->GetXaxis()->SetTitle("hCut_2Track_deltaInvMass");
        temp->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.deltaInvMass.get_width()));
        hCut_2Track_deltaInvMass.push_back(temp);
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

void CCProtonPi0_CutList::init_nCutVectors()
{
    for (int i = 0; i < nTopologies; i++){
        nCut_All.push_back(CCProtonPi0_Cut());
        nCut_Vertex_None.push_back(CCProtonPi0_Cut());
        nCut_Vertex_Not_Reconstructable.push_back(CCProtonPi0_Cut()); 
        nCut_Vertex_Not_Fiducial.push_back(CCProtonPi0_Cut());
        nCut_Muon_None.push_back(CCProtonPi0_Cut());              
        nCut_Muon_Not_Plausible.push_back(CCProtonPi0_Cut());
        nCut_Muon_Charge.push_back(CCProtonPi0_Cut());
        nCut_Vertex_Michel_Exist.push_back(CCProtonPi0_Cut()); 
        nCut_EndPoint_Michel_Exist.push_back(CCProtonPi0_Cut());
        nCut_secEndPoint_Michel_Exist.push_back(CCProtonPi0_Cut());
        nCut_PreFilter_Pi0.push_back(CCProtonPi0_Cut());
        nCut_ConeBlobs.push_back(CCProtonPi0_Cut());
        nCut_BlobDirectionBad.push_back(CCProtonPi0_Cut());
        nCut_BlobsBad.push_back(CCProtonPi0_Cut());
        nCut_Photon1DistanceLow.push_back(CCProtonPi0_Cut());
        nCut_Photon2DistanceLow.push_back(CCProtonPi0_Cut());
        nCut_Pi0_invMass.push_back(CCProtonPi0_Cut());
        nCut_Particle_None.push_back(CCProtonPi0_Cut());
        nCut_Proton_None.push_back(CCProtonPi0_Cut());            
        nCut_ProtonScore.push_back(CCProtonPi0_Cut());
        nCut_DeltaInvMass.push_back(CCProtonPi0_Cut());
        nCut_beamEnergy.push_back(CCProtonPi0_Cut());
        nCut_UnusedE.push_back(CCProtonPi0_Cut());
    }
}

void CCProtonPi0_CutList::SetCutNames()
{
    for (int i = 0; i < nTopologies; i++){
        nCut_All[i].set_Name("All");
        nCut_Vertex_None[i].set_Name("Vertex_None");
        nCut_Vertex_Not_Reconstructable[i].set_Name("Vertex_Not_Reconstructable"); 
        nCut_Vertex_Not_Fiducial[i].set_Name("Vertex_Not_Fiducial");
        nCut_Muon_None[i].set_Name("Muon_None");              
        nCut_Muon_Not_Plausible[i].set_Name("Muon_Not_Plausible");
        nCut_Muon_Charge[i].set_Name("Muon_Charge");
        nCut_Vertex_Michel_Exist[i].set_Name("Vertex_Michel_Exist"); 
        nCut_EndPoint_Michel_Exist[i].set_Name("EndPoint_Michel_Exist");
        nCut_secEndPoint_Michel_Exist[i].set_Name("secEndPoint_Michel_Exist");
        nCut_PreFilter_Pi0[i].set_Name("PreFilter_Pi0");
        nCut_ConeBlobs[i].set_Name("ConeBlobs");
        nCut_BlobDirectionBad[i].set_Name("BlobDirectionBad");
        nCut_BlobsBad[i].set_Name("BlobsBad");
        nCut_Photon1DistanceLow[i].set_Name("Photon1DistanceLow");
        nCut_Photon2DistanceLow[i].set_Name("Photon2DistanceLow");
        nCut_Pi0_invMass[i].set_Name("Pi0_invMass");
        nCut_Particle_None[i].set_Name("Particle_None");
        nCut_Proton_None[i].set_Name("Proton_None");            
        nCut_ProtonScore[i].set_Name("Proton_Score");
        nCut_DeltaInvMass[i].set_Name("Delta_invMass");
        nCut_beamEnergy[i].set_Name("beamEnergy");
        nCut_UnusedE[i].set_Name("UnusedE");
    }
}

void CCProtonPi0_CutList::OpenTextFiles(bool isMC)
{
    std::string type;

    if (isMC) type = "CutTable_MC_";
    else type = "CutTable_Data_";

    // Open Cut Files
    cutFile[0] = Folder_List::output + Folder_List::textOut + type + "1Track.txt";
    cutFile[1] = Folder_List::output + Folder_List::textOut + type + "2Track.txt";
    
    for (int i = 0; i < nTopologies; i++){
        cutText[i].open( cutFile[i].c_str() );
        if( !cutText[i].is_open() ){
            cerr<<"Cannot open output text file: "<<cutFile[i]<<endl;
            exit(1);
        }else{
            cout<<"\t"<<cutFile[i]<<endl;
        }
    }
     
}

void CCProtonPi0_CutList::writeCutTableHeader(int t)
{
    cutText[t]<<std::left;
    
    cutText[t].width(35); cutText[t]<<"Cut"<<" "; 
    
    cutText[t].width(12); cutText[t]<<"N(Events)"<<" ";    
    cutText[t].width(12); cutText[t]<<"N(Signal)"<<" ";      
    cutText[t].width(12); cutText[t]<<"Eff(All)"<<" ";      
    cutText[t].width(12); cutText[t]<<"Eff(MINOS)"<<" ";      
    cutText[t].width(12); cutText[t]<<"Purity"<<" ";
    cutText[t].width(12); cutText[t]<<"N(Study1)"<<" "; 
    cutText[t].width(12); cutText[t]<<"N(Study2)"<<" "; 
    cutText[t]<<endl;
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

void CCProtonPi0_CutList::formCutVectors()
{   

    nCutVector_Common.push_back(nCut_All);
    nCutVector_Common.push_back(nCut_Vertex_None);
    nCutVector_Common.push_back(nCut_Vertex_Not_Reconstructable); 
    nCutVector_Common.push_back(nCut_Vertex_Not_Fiducial);
    nCutVector_Topology.push_back(nCut_Muon_None);              
    nCutVector_Topology.push_back(nCut_Muon_Not_Plausible);
    nCutVector_Topology.push_back(nCut_Muon_Charge);
    nCutVector_Topology.push_back(nCut_Vertex_Michel_Exist); 
    nCutVector_Topology.push_back(nCut_EndPoint_Michel_Exist);
    nCutVector_Topology.push_back(nCut_secEndPoint_Michel_Exist);
    nCutVector_Topology.push_back(nCut_PreFilter_Pi0);
    nCutVector_Topology.push_back(nCut_ConeBlobs);
    nCutVector_Topology.push_back(nCut_BlobDirectionBad);
    nCutVector_Topology.push_back(nCut_BlobsBad);
    nCutVector_Topology.push_back(nCut_Photon1DistanceLow);
    nCutVector_Topology.push_back(nCut_Photon2DistanceLow);
    nCutVector_Topology.push_back(nCut_Pi0_invMass);
    nCutVector_Topology.push_back(nCut_Particle_None);
    nCutVector_Topology.push_back(nCut_Proton_None);
    nCutVector_Topology.push_back(nCut_ProtonScore);
    nCutVector_Topology.push_back(nCut_DeltaInvMass);
    nCutVector_Topology.push_back(nCut_beamEnergy);
    nCutVector_Topology.push_back(nCut_UnusedE);
}

void CCProtonPi0_CutList::writeCutTable()
{
    formCutVectors();
    
    for (int t = 0; t < nTopologies; t++){
        cout<<">> Writing "<<cutFile[t]<<endl;
    
        writeCutTableHeader(t);
        writeCutTableRows(t,nCutVector_Common);
        cutText[t]<<endl;
        writeCutTableRows(t,nCutVector_Topology);
    }
}

void CCProtonPi0_CutList::writeCutTableRows(int t, vector< vector<CCProtonPi0_Cut> > &nCutVector)
{
    // Write General Cuts upto Proton Reconstruction    
    for( unsigned int i = 0; i < nCutVector.size(); i++){
        writeSingleRow(t, nCutVector[i][t]);    
    }

    cutText[t]<<endl;
}

void CCProtonPi0_CutList::writeSingleRow(int t, CCProtonPi0_Cut& currentCut)
{
    double eff_AllSignal;
    double eff_MINOS;
    double purity;    

    eff_AllSignal = getCutEfficiency(currentCut,nCut_All[t]);
    eff_MINOS = getCutEfficiency(currentCut,nCut_Muon_None[t]);
    purity = getCutPurity(currentCut);
            
    cutText[t].width(35); cutText[t]<<currentCut.get_Name()<<" ";
    cutText[t].width(12); cutText[t]<<currentCut.nEvent.getCount()<<" ";
    
    // Total Signal
    cutText[t].width(12); cutText[t]<<currentCut.nSignal.getCount()<<" ";

    // Efficiency
    if ( eff_AllSignal <= 100){
        cutText[t].width(12); cutText[t]<<eff_AllSignal<<" ";
    }else{
        cutText[t].width(12); cutText[t]<<"N/A"<<" ";    
    }

    if ( eff_MINOS <= 100){
        cutText[t].width(12); cutText[t]<<eff_MINOS<<" ";
    }else{
        cutText[t].width(12); cutText[t]<<"N/A"<<" ";    
    }    
    
    // Purity
    cutText[t].width(12); cutText[t]<<purity<<" ";

    // Number of Events which are studied
    cutText[t].width(12); cutText[t]<<currentCut.nStudy1.getCount()<<" ";
    cutText[t].width(12); cutText[t]<<currentCut.nStudy2.getCount()<<" ";
    
    cutText[t]<<endl;
}

CCProtonPi0_CutList::~CCProtonPi0_CutList()
{
    for (int i = 0; i < nTopologies; i++){
        cutText[i].close(); 
    }
}

void CCProtonPi0_CutList::writeHistograms()
{
    std::cout<<">> Writing "<<rootDir<<std::endl;
    f->cd();
    for (int i = 0; i < nHistograms; i++){
        // 1 Track
        hCut_1Track_Michel[i]->Write();
        hCut_1Track_eVis_nuclearTarget[i]->Write();
        hCut_1Track_eVis_other[i]->Write();
        hCut_1Track_pi0invMass[i]->Write();
        hCut_1Track_gamma1ConvDist[i]->Write();
        hCut_1Track_gamma2ConvDist[i]->Write();
        hCut_1Track_neutrinoE[i]->Write();
        hCut_1Track_UnusedE[i]->Write();
       
        // 2 Track 
        hCut_2Track_Michel[i]->Write();
        hCut_2Track_eVis_nuclearTarget[i]->Write();
        hCut_2Track_eVis_other[i]->Write();
        hCut_2Track_pi0invMass[i]->Write();
        hCut_2Track_gamma1ConvDist[i]->Write();
        hCut_2Track_gamma2ConvDist[i]->Write();
        hCut_2Track_neutrinoE[i]->Write();
        hCut_2Track_UnusedE[i]->Write();
        hCut_2Track_protonScore_pIDDiff[i]->Write();
        hCut_2Track_protonScore_LLR[i]->Write();
        hCut_2Track_deltaInvMass[i]->Write();
    }

    // MC Only
    mc_w_DIS->Write();
    mc_w_RES->Write();
    mc_w_CCQE->Write();

    f->Close();
}

#endif

