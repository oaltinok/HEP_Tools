/*
    See CCProtonPi0_Interaction.h header for Class Information
*/
#ifndef CCProtonPi0_Interaction_cpp
#define CCProtonPi0_Interaction_cpp

#include "CCProtonPi0_Interaction.h"

using namespace PlotUtils;

CCProtonPi0_Interaction::CCProtonPi0_Interaction(bool isModeReduce, bool isMC, std::string ana_folder) : CCProtonPi0_NTupleAnalysis()
{
    std::cout<<"Initializing CCProtonPi0_Interaction"<<std::endl;
    
    if(isModeReduce){
        std::cout<<"\tNTuple Reduce Mode -- Will not create ROOT Files"<<std::endl;
    }else{
        if (isMC) rootDir = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + ana_folder + "Interaction.root";
        else rootDir = Folder_List::rootOut + Folder_List::Data + Folder_List::analyzed  + ana_folder + "Interaction.root";
        
        std::cout<<"\tRoot File: "<<rootDir<<std::endl;
 
        // Create Root File 
        f = new TFile(rootDir.c_str(),"RECREATE");
        
        initHistograms();
    }
    
    std::cout<<"Done!"<<std::endl;
}


void CCProtonPi0_Interaction::initHistograms()
{
    MnvH1D* temp = NULL;

    for (int i = 0; i < nHistograms; i++){
        temp = new MnvH1D( Form("%s_%d","Enu_True",i),"True Beam Energy",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
        temp->GetXaxis()->SetTitle("True E_{#nu}[GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.beamE.get_width()));
        Enu_True.push_back(temp);
        
        temp = new MnvH1D( Form("%s_%d","Enu_1Track",i),"Reconstructed Beam Energy - 1 Track",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed E_{#nu} - 1Track [GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.beamE.get_width()));
        Enu_1Track.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","Enu_1Track_Alt",i),"Reconstructed Beam Energy - 1 Track Alternative",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed E_{#nu} - 1Track Alternative [GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.beamE.get_width()));
        Enu_1Track_Alt.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","Enu_2Track",i),"Reconstructed Beam Energy - 2 Track",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed E_{#nu} - 2 Track [GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.beamE.get_width()));
        Enu_2Track.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","QSq",i),"Reconstructed Q^{2}",binList.q2.get_nBins(), binList.q2.get_min(), binList.q2.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed Q^{2} [GeV^{2}]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.q2.get_width()));
        QSq.push_back(temp);    

        temp = new MnvH1D( Form("%s_%d","WSq",i),"Reconstructed W^{2}",binList.wSq.get_nBins(), binList.wSq.get_min(), binList.wSq.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed W^{2} [GeV^{2}]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.wSq.get_width()));
        WSq.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","W",i),"Reconstructed W",binList.w.get_nBins(), binList.w.get_min(), binList.w.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed W [GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.w.get_width()));
        W.push_back(temp);
 
        temp = new MnvH1D( Form("%s_%d","vertex_energy_1Track",i),"Vertex Blob Energy (r = 90mm)",binList.vertex_energy.get_nBins(), binList.vertex_energy.get_min(), binList.vertex_energy.get_max() );
        temp->GetXaxis()->SetTitle("Vertex Blob Energy [MeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.vertex_energy.get_width()));
        vertex_energy_1Track.push_back(temp);
  
        temp = new MnvH1D( Form("%s_%d","vertex_energy_2Track",i),"Vertex Blob Energy (r = 90mm)",binList.vertex_energy.get_nBins(), binList.vertex_energy.get_min(), binList.vertex_energy.get_max() );
        temp->GetXaxis()->SetTitle("Vertex Blob Energy [MeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.vertex_energy.get_width()));
        vertex_energy_2Track.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","vertex_evis_1Track",i),"Vertex Blob Visible Energy (r = 90mm)",binList.vertex_evis.get_nBins(), binList.vertex_evis.get_min(), binList.vertex_evis.get_max() );
        temp->GetXaxis()->SetTitle("Vertex Blob Visible Energy [MeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.vertex_evis.get_width()));
        vertex_evis_1Track.push_back(temp);
  
        temp = new MnvH1D( Form("%s_%d","vertex_evis_2Track",i),"Vertex Blob Visible Energy (r = 90mm)",binList.vertex_evis.get_nBins(), binList.vertex_evis.get_min(), binList.vertex_evis.get_max() );
        temp->GetXaxis()->SetTitle("Vertex Blob Visible Energy [MeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.vertex_evis.get_width()));
        vertex_evis_2Track.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","extra_energy_1Track",i),"Extra Energy (r = 300mm)",binList.extra_energy.get_nBins(), binList.extra_energy.get_min(), binList.extra_energy.get_max() );
        temp->GetXaxis()->SetTitle("Extra Energy [MeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.extra_energy.get_width()));
        extra_energy_1Track.push_back(temp);
 
        temp = new MnvH1D( Form("%s_%d","extra_energy_2Track",i),"Extra Energy (r = 300mm)",binList.extra_energy.get_nBins(), binList.extra_energy.get_min(), binList.extra_energy.get_max() );
        temp->GetXaxis()->SetTitle("Extra Energy [MeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.extra_energy.get_width()));
        extra_energy_2Track.push_back(temp);
      
        temp = new MnvH1D( Form("%s_%d","deltaInvMass",i),"Reconstructed #Delta^{+} Invariant Mass",binList.deltaInvMass.get_nBins(), binList.deltaInvMass.get_min(), binList.deltaInvMass.get_max() );
        temp->GetXaxis()->SetTitle("Reconstructed #Delta^{+} Inv. Mass [GeV]");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f [MeV] ",binList.deltaInvMass.get_width()));
        deltaInvMass.push_back(temp);

        temp = new MnvH1D( Form("%s_%d","vertex_z",i),"Reconstructed Interaction Vertex",binList.vertex_z.get_nBins(), binList.vertex_z.get_min(), binList.vertex_z.get_max() );
        temp->GetXaxis()->SetTitle("z = 4293 Target, #bf{z = 5810 Interaction Region}, z = 8614 ECAL, z = 9088 HCAL");
        temp->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.vertex_z.get_width()));
        vertex_z.push_back(temp);
}
    
    // MC Only Histograms
    final_mc_w_DIS = new TH1D( "final_mc_w_DIS","True W for DIS",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max() );
    final_mc_w_DIS->GetXaxis()->SetTitle("True W for DIS [GeV]");
    final_mc_w_DIS->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.mc_w.get_width()));
    
    final_mc_w_RES = new TH1D( "final_mc_w_RES","True W for RES",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max() );
    final_mc_w_RES->GetXaxis()->SetTitle("True W for RES [GeV]");
    final_mc_w_RES->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.mc_w.get_width()));
    
    final_mc_w_CCQE = new TH1D( "final_mc_w_CCQE","True W for CCQE",binList.mc_w.get_nBins(), binList.mc_w.get_min(), binList.mc_w.get_max() );
    final_mc_w_CCQE->GetXaxis()->SetTitle("True W for CCQE [GeV]");
    final_mc_w_CCQE->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.mc_w.get_width()));

    // Short Proton True Information
    proton_true_P_1Track = new TH1D( "proton_true_P_1Track","Short Proton True Momentum",binList.short_proton_P.get_nBins(), binList.short_proton_P.get_min(), binList.short_proton_P.get_max() );
    proton_true_P_1Track->GetXaxis()->SetTitle("Short Proton P_{True} [MeV]");
    proton_true_P_1Track->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.short_proton_P.get_width()));

    proton_true_KE_1Track = new TH1D( "proton_true_KE_1Track","Short Proton True Momentum",binList.short_proton_KE.get_nBins(), binList.short_proton_KE.get_min(), binList.short_proton_KE.get_max() );
    proton_true_KE_1Track->GetXaxis()->SetTitle("Short Proton P_{True} [MeV]");
    proton_true_KE_1Track->GetYaxis()->SetTitle(Form("Candidates / %3.2f ",binList.short_proton_KE.get_width()));

    n_ejected_nucleons_1Track = new TH1D( "n_ejected_nucleons_1Track","N(Nucleons) out of Nucleus",binList.multiplicity.get_nBins(), binList.multiplicity.get_min(), binList.multiplicity.get_max() );
    n_ejected_nucleons_1Track->GetXaxis()->SetTitle("N(Nucleons) out of Nucleus");
    n_ejected_nucleons_1Track->GetYaxis()->SetTitle("N(Events)");

    n_ejected_nucleons_2Track = new TH1D( "n_ejected_nucleons_2Track","N(Nucleons) out of Nucleus",binList.multiplicity.get_nBins(), binList.multiplicity.get_min(), binList.multiplicity.get_max() );
    n_ejected_nucleons_2Track->GetXaxis()->SetTitle("N(Nucleons) out of Nucleus");
    n_ejected_nucleons_2Track->GetYaxis()->SetTitle("N(Events)");

    // Neutrino Energy Error
    Enu_1Track_Error = new TH1D("Enu_1Track_Error","Neutrino Energy Error - 1 Track",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    Enu_1Track_Error->GetXaxis()->SetTitle("(E_{#nu}^{Reco}-E_{#nu}^{True})/E_{#nu}^{True}");
    Enu_1Track_Error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));

    Enu_1Track_Alt_Error = new TH1D("Enu_1Track_Alt_Error","Neutrino Energy Error - 1 Track Alternative",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    Enu_1Track_Alt_Error->GetXaxis()->SetTitle("(E_{#nu}^{Reco}-E_{#nu}^{True})/E_{#nu}^{True}");
    Enu_1Track_Alt_Error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));

    Enu_2Track_Error = new TH1D("Enu_2Track_Error","Neutrino Energy Error - 2 Track",binList.error.get_nBins(), binList.error.get_min(), binList.error.get_max() );
    Enu_2Track_Error->GetXaxis()->SetTitle("(E_{#nu}^{Reco}-E_{#nu}^{True})/E_{#nu}^{True}");
    Enu_2Track_Error->GetYaxis()->SetTitle(Form("Events / %3.2f ",binList.error.get_width()));

    Enu_1Track_True = new TH2D("Enu_1Track_True","Neutrino Energy True vs 1 Track",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max(),binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    Enu_1Track_True->GetXaxis()->SetTitle("E_{#nu}^{Reco} [GeV]");
    Enu_1Track_True->GetYaxis()->SetTitle("E_{#nu}^{True} [GeV]");

    Enu_1Track_Alt_True = new TH2D("Enu_1Track_Alt_True","Neutrino Energy True vs 1 Track Alternative",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max(),binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    Enu_1Track_Alt_True->GetXaxis()->SetTitle("E_{#nu}^{Reco} [GeV]");
    Enu_1Track_Alt_True->GetYaxis()->SetTitle("E_{#nu}^{True} [GeV]");

    Enu_1Track_1Track_Alt = new TH2D("Enu_1Track_1Track_Alt","Neutrino Energy 1Track vs 1 Track Alternative",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max(),binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    Enu_1Track_1Track_Alt->GetXaxis()->SetTitle("E_{#nu}^{Reco} [GeV]");
    Enu_1Track_1Track_Alt->GetYaxis()->SetTitle("E_{#nu}^{Reco} Alternative [GeV]");

    Enu_2Track_True = new TH2D("Enu_2Track_True","Neutrino Energy True vs 2 Track",binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max(),binList.beamE.get_nBins(), binList.beamE.get_min(), binList.beamE.get_max() );
    Enu_2Track_True->GetXaxis()->SetTitle("E_{#nu}^{Reco} [GeV]");
    Enu_2Track_True->GetYaxis()->SetTitle("E_{#nu}^{True} [GeV]");
}

void CCProtonPi0_Interaction::writeHistograms()
{
    std::cout<<">> Writing "<<rootDir<<std::endl;
   
    f->cd();

    for (int i = 0; i < nHistograms; i++){
        // Event Kinematics
        Enu_True[i]->Write();
        Enu_1Track[i]->Write();
        Enu_1Track_Alt[i]->Write();
        Enu_2Track[i]->Write();
        QSq[i]->Write();
        WSq[i]->Write();
        W[i]->Write();
       
        // Vertex & Extra Energy
        vertex_energy_1Track[i]->Write();
        vertex_energy_2Track[i]->Write();
        vertex_evis_1Track[i]->Write();
        vertex_evis_2Track[i]->Write();
        extra_energy_1Track[i]->Write();
        extra_energy_2Track[i]->Write();

        // Other Event Parameters 
        deltaInvMass[i]->Write();
    }
    
    // MC Only Histograms
    final_mc_w_DIS->Write();
    final_mc_w_RES->Write();
    final_mc_w_CCQE->Write();

    proton_true_P_1Track->Write();
    proton_true_KE_1Track->Write();
    
    n_ejected_nucleons_1Track->Write();
    n_ejected_nucleons_2Track->Write();

    Enu_1Track_Error->Write();
    Enu_1Track_Alt_Error->Write();
    Enu_2Track_Error->Write();

    Enu_1Track_True->Write();
    Enu_1Track_Alt_True->Write();
    Enu_1Track_1Track_Alt->Write();
    Enu_2Track_True->Write();

    f->Close();
}



#endif
