/*
    See CCProtonPi0_Pi0Blob.h header for Class Information
*/
#ifndef CCProtonPi0_Pi0Blob_cpp
#define CCProtonPi0_Pi0Blob_cpp

#include "CCProtonPi0_Pi0Blob.h"

using namespace PlotUtils;

CCProtonPi0_Pi0Blob::CCProtonPi0_Pi0Blob(bool isModeReduce, bool isMC, std::string ana_folder) : CCProtonPi0_NTupleAnalysis()
{
    std::cout<<"Initializing CCProtonPi0_Pi0Blob"<<std::endl;

    if(isModeReduce){
        std::cout<<"\tNTuple Reduce Mode -- Will not create ROOT Files"<<std::endl;
    }else{
        // File Locations
        if (isMC) rootDir = Folder_List::rootOut + Folder_List::MC + Folder_List::analyzed + ana_folder + "Pi0Blob.root";
        else rootDir = Folder_List::rootOut + Folder_List::Data + Folder_List::analyzed + ana_folder + "Pi0Blob.root";      

        std::cout<<"\tRoot File: "<<rootDir<<std::endl;
     
        // Create Root File 
        f = new TFile(rootDir.c_str(),"RECREATE");

        initBins();
        initHistograms();
    }  
    std::cout<<"Done!"<<std::endl;
}

void CCProtonPi0_Pi0Blob::initHistograms()
{
    MnvH1D* temp = NULL;

    for (int i = 0; i < nHistograms; i++){
        temp = new MnvH1D(Form("%s_%d","g1_evis_most_pdg",i),Form("%s %d","Most Evis PDGi for Gamma", 1), binList.multiplicity.get_nBins(), binList.multiplicity.get_min(), binList.multiplicity.get_max());
        temp->GetXaxis()->SetTitle("0: #pi^{0}, 1: #pi^{+}, 2: #pi^{-}, 3: n, 4: p, 5: #mu^{-}, 6:Other");
        temp->GetYaxis()->SetTitle("N(Events)");
        g1_evis_most_pdg.push_back(temp);

        temp = new MnvH1D(Form("%s_%d","g1_evis_total_truth",i),Form("%s %d","Total Visible Energy for Gamma",1),bin_blob_energy.get_nBins(), bin_blob_energy.get_min(), bin_blob_energy.get_max() );
        temp->GetXaxis()->SetTitle("Total Visible Energy [GeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        g1_evis_total_truth.push_back(temp);

        temp = new MnvH1D(Form("%s_%d","g1_evis_frac_pizero",i),Form("%s %d","pizero Visible Energy Fraction for Gamma",1),binList.fraction.get_nBins(), binList.fraction.get_min(), binList.fraction.get_max() );
        temp->GetXaxis()->SetTitle("Evis pizero / Evis Total");
        temp->GetYaxis()->SetTitle("N(Events)");
        g1_evis_frac_pizero.push_back(temp);

        temp = new MnvH1D(Form("%s_%d","g1_evis_frac_piplus",i),Form("%s %d","piplus Visible Energy Fraction for Gamma",1),binList.fraction.get_nBins(), binList.fraction.get_min(), binList.fraction.get_max() );
        temp->GetXaxis()->SetTitle("Evis piplus / Evis Total");
        temp->GetYaxis()->SetTitle("N(Events)");
        g1_evis_frac_piplus.push_back(temp);

        temp = new MnvH1D(Form("%s_%d","g1_evis_frac_piminus",i),Form("%s %d","piminus Visible Energy Fraction for Gamma",1),binList.fraction.get_nBins(), binList.fraction.get_min(), binList.fraction.get_max() );
        temp->GetXaxis()->SetTitle("Evis piminus / Evis Total");
        temp->GetYaxis()->SetTitle("N(Events)");
        g1_evis_frac_piminus.push_back(temp);

        temp = new MnvH1D(Form("%s_%d","g1_evis_frac_proton",i),Form("%s %d","proton Visible Energy Fraction for Gamma",1),binList.fraction.get_nBins(), binList.fraction.get_min(), binList.fraction.get_max() );
        temp->GetXaxis()->SetTitle("Evis proton / Evis Total");
        temp->GetYaxis()->SetTitle("N(Events)");
        g1_evis_frac_proton.push_back(temp);

        temp = new MnvH1D(Form("%s_%d","g1_evis_frac_neutron",i),Form("%s %d","neutron Visible Energy Fraction for Gamma",1),binList.fraction.get_nBins(), binList.fraction.get_min(), binList.fraction.get_max() );
        temp->GetXaxis()->SetTitle("Evis neutron / Evis Total");
        temp->GetYaxis()->SetTitle("N(Events)");
        g1_evis_frac_neutron.push_back(temp);

        temp = new MnvH1D(Form("%s_%d","g1_evis_frac_muon",i),Form("%s %d","muon Visible Energy Fraction for Gamma",1),binList.fraction.get_nBins(), binList.fraction.get_min(), binList.fraction.get_max() );
        temp->GetXaxis()->SetTitle("Evis muon / Evis Total");
        temp->GetYaxis()->SetTitle("N(Events)");
        g1_evis_frac_muon.push_back(temp);

        // Gamma 2
        temp = new MnvH1D(Form("%s_%d","g2_evis_most_pdg",i),Form("%s %d","Most Evis PDGi for Gamma", 2), binList.multiplicity.get_nBins(), binList.multiplicity.get_min(), binList.multiplicity.get_max());
        temp->GetXaxis()->SetTitle("0: #pi^{0}, 1: #pi^{+}, 2: #pi^{-}, 3: n, 4: p, 5: #mu^{-}, 6:Other");
        temp->GetYaxis()->SetTitle("N(Events)");
        g2_evis_most_pdg.push_back(temp);

        temp = new MnvH1D(Form("%s_%d","g2_evis_total_truth",i),Form("%s %d","Total Visible Energy for Gamma",2),bin_blob_energy.get_nBins(), bin_blob_energy.get_min(), bin_blob_energy.get_max() );
        temp->GetXaxis()->SetTitle("Total Visible Energy [GeV]");
        temp->GetYaxis()->SetTitle("N(Events)");
        g2_evis_total_truth.push_back(temp);

        temp = new MnvH1D(Form("%s_%d","g2_evis_frac_pizero",i),Form("%s %d","pizero Visible Energy Fraction for Gamma",2),binList.fraction.get_nBins(), binList.fraction.get_min(), binList.fraction.get_max() );
        temp->GetXaxis()->SetTitle("Evis pizero / Evis Total");
        temp->GetYaxis()->SetTitle("N(Events)");
        g2_evis_frac_pizero.push_back(temp);

        temp = new MnvH1D(Form("%s_%d","g2_evis_frac_piplus",i),Form("%s %d","piplus Visible Energy Fraction for Gamma",2),binList.fraction.get_nBins(), binList.fraction.get_min(), binList.fraction.get_max() );
        temp->GetXaxis()->SetTitle("Evis piplus / Evis Total");
        temp->GetYaxis()->SetTitle("N(Events)");
        g2_evis_frac_piplus.push_back(temp);

        temp = new MnvH1D(Form("%s_%d","g2_evis_frac_piminus",i),Form("%s %d","piminus Visible Energy Fraction for Gamma",2),binList.fraction.get_nBins(), binList.fraction.get_min(), binList.fraction.get_max() );
        temp->GetXaxis()->SetTitle("Evis piminus / Evis Total");
        temp->GetYaxis()->SetTitle("N(Events)");
        g2_evis_frac_piminus.push_back(temp);

        temp = new MnvH1D(Form("%s_%d","g2_evis_frac_proton",i),Form("%s %d","proton Visible Energy Fraction for Gamma",2),binList.fraction.get_nBins(), binList.fraction.get_min(), binList.fraction.get_max() );
        temp->GetXaxis()->SetTitle("Evis proton / Evis Total");
        temp->GetYaxis()->SetTitle("N(Events)");
        g2_evis_frac_proton.push_back(temp);

        temp = new MnvH1D(Form("%s_%d","g2_evis_frac_neutron",i),Form("%s %d","neutron Visible Energy Fraction for Gamma",2),binList.fraction.get_nBins(), binList.fraction.get_min(), binList.fraction.get_max() );
        temp->GetXaxis()->SetTitle("Evis neutron / Evis Total");
        temp->GetYaxis()->SetTitle("N(Events)");
        g2_evis_frac_neutron.push_back(temp);

        temp = new MnvH1D(Form("%s_%d","g2_evis_frac_muon",i),Form("%s %d","muon Visible Energy Fraction for Gamma",2),binList.fraction.get_nBins(), binList.fraction.get_min(), binList.fraction.get_max() );
        temp->GetXaxis()->SetTitle("Evis muon / Evis Total");
        temp->GetYaxis()->SetTitle("N(Events)");
        g2_evis_frac_muon.push_back(temp);

        temp = new MnvH1D(Form("%s_%d","captured_evis_frac_all",i),"All Events: Pi0 Evis Capture Fraction",binList.fraction2.get_nBins(), binList.fraction2.get_min(), binList.fraction2.get_max() );
        temp->GetXaxis()->SetTitle("Pi0 Evis Captured / Pi0 Evis Total");
        temp->GetYaxis()->SetTitle("N(Events)");
        captured_evis_frac_all.push_back(temp);

        temp = new MnvH1D(Form("%s_%d","captured_evis_frac_signal",i),"Signal Events: Pi0 Evis Capture Fraction",binList.fraction2.get_nBins(), binList.fraction2.get_min(), binList.fraction2.get_max() );
        temp->GetXaxis()->SetTitle("Pi0 Evis Captured / Pi0 Evis Total");
        temp->GetYaxis()->SetTitle("N(Events)");
        captured_evis_frac_signal.push_back(temp);

        // Reco Values
        temp = new MnvH1D(Form("%s_%d","g1_nPlanes",i),Form("%s %d","N(Planes) for Gamma ",1),binList.shower_length.get_nBins(), binList.shower_length.get_min(), binList.shower_length.get_max() );
        temp->GetXaxis()->SetTitle("N(Planes)");
        temp->GetYaxis()->SetTitle("N(Events)");
        g1_nPlanes.push_back(temp);

        temp = new MnvH1D(Form("%s_%d","g2_nPlanes",i),Form("%s %d","N(Planes) for Gamma ",2),binList.shower_length.get_nBins(), binList.shower_length.get_min(), binList.shower_length.get_max() );
        temp->GetXaxis()->SetTitle("N(Planes)");
        temp->GetYaxis()->SetTitle("N(Events)");
        g2_nPlanes.push_back(temp);
    } //end loop all Histograms 
}

void CCProtonPi0_Pi0Blob::initBins()
{
    // Initialize Bins
    bin_blob_energy.setBin(100,0.0,0.7);
}

void CCProtonPi0_Pi0Blob::writeHistograms()
{
    std::cout<<">> Writing "<<rootDir<<std::endl;
    f->cd();
    
    // Write Truth Match Histograms
    for (int i = 0; i < nHistograms; i++){
        g1_evis_most_pdg[i]->Write();
        g1_evis_total_truth[i]->Write();
        g1_evis_frac_pizero[i]->Write();
        g1_evis_frac_piplus[i]->Write();
        g1_evis_frac_piminus[i]->Write();
        g1_evis_frac_proton[i]->Write();
        g1_evis_frac_neutron[i]->Write();
        g1_evis_frac_muon[i]->Write();
        
        g2_evis_most_pdg[i]->Write();
        g2_evis_total_truth[i]->Write();
        g2_evis_frac_pizero[i]->Write();
        g2_evis_frac_piplus[i]->Write();
        g2_evis_frac_piminus[i]->Write();
        g2_evis_frac_proton[i]->Write();
        g2_evis_frac_neutron[i]->Write();
        g2_evis_frac_muon[i]->Write();
        
        captured_evis_frac_all[i]->Write();
        captured_evis_frac_signal[i]->Write();

        g1_nPlanes[i]->Write();
        g2_nPlanes[i]->Write();
    }

    f->Close();
}


#endif

