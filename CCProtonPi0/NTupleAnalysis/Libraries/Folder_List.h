/*
================================================================================
Namespace: Folder_List
    Folder List for easier folder management



Main Directory:
Libraries/

Usage:
> #include "Libraries/Folder_List.h" 
> PDG_List::output      // Returns default Output Directory

Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
================================================================================
*/

#ifndef Folder_List_h
#define Folder_List_h

#include <string>

namespace Folder_List
{
    //-------------------------------------------------------------------------
    // Data files on /minerva/data/ Disk
    //-------------------------------------------------------------------------
    const std::string rootOut = "/minerva/data/users/oaltinok/NTupleAnalysis/";
    const std::string MC = "MC/";
    const std::string Data = "Data/";
    const std::string ParticleCannon = "ParticleCannon/";
    
    const std::string analyzed = "Analyzed/";
    const std::string reduced = "Reduced/"; 

    //-------------------------------------------------------------------------
    // NTuple Analysis Output
    //-------------------------------------------------------------------------
    // Default Folder for Output
    const std::string output = "Output/";
    
    // SubFolders for Output
    const std::string textOut = "TextFiles/";
    const std::string plotOut = "Plots/";
    
    // Folder Branches
    const std::string signal = "Signal/";
    const std::string background = "Background/";
    const std::string allEvents = "AllEvents/";
    const std::string other = "Other/";

    //-------------------------------------------------------------------------
    // ROOT Dirs
    //-------------------------------------------------------------------------
    // Set Truth ROOT Dir
    const std::string rootDir_Truth_mc = rootOut + MC + analyzed + "TruthHistograms.root";
    const std::string rootDir_Truth_data = "";

    // Set Other Studies ROOT Dir
    const std::string rootDir_OtherStudies_mc = "/minerva/data/users/oaltinok/NTupleAnalysis/ParticleCannon/PC_Test.root"; 
    const std::string rootDir_OtherStudies_data = "";

    const std::string rootDir_CrossSection = rootOut + Data + analyzed + "CrossSection.root";
    
    // Set MC Root Dir;
    const std::string rootDir_CutHists_mc = rootOut + MC + analyzed + "CutHistograms.root";
    const std::string rootDir_Interaction_mc = rootOut + MC + analyzed + "Interaction.root";
    const std::string rootDir_Muon_mc = rootOut + MC + analyzed + "Muon.root";
    const std::string rootDir_Proton_mc = rootOut + MC + analyzed + "Proton.root";
    const std::string rootDir_Pion_mc = rootOut + MC + analyzed + "Pion.root";
    const std::string rootDir_Pi0Blob_mc = rootOut + MC + analyzed + "Pi0Blob.root";

    // Set Data Root Dir
    const std::string rootDir_CutHists_data = rootOut + Data + analyzed + "CutHistograms.root";
    const std::string rootDir_Interaction_data = rootOut + Data + analyzed + "Interaction.root";
    const std::string rootDir_Muon_data = rootOut + Data + analyzed + "Muon.root";
    const std::string rootDir_Proton_data = rootOut + Data + analyzed + "Proton.root";
    const std::string rootDir_Pion_data = rootOut + Data + analyzed + "Pion.root";
    const std::string rootDir_Pi0Blob_data = rootOut + Data + analyzed + "Pi0Blob.root";

    //-------------------------------------------------------------------------
    // Plot Dirs
    //-------------------------------------------------------------------------
    const std::string plotDir_Efficiency = output + plotOut + "Efficiency/";
    const std::string plotDir_OtherStudies = output + plotOut + "OtherStudies/";
    const std::string plotDir_CutHists = output + plotOut + "CutHists/";
    const std::string plotDir_Interaction = output + plotOut + "Interaction/";
    const std::string plotDir_Muon = output + plotOut + "Muon/";
    const std::string plotDir_Proton = output + plotOut + "Proton/";
    const std::string plotDir_Pion = output + plotOut + "Pion/";
    const std::string plotDir_Pi0Blob = output + plotOut + "Pi0Blob/";
    const std::string plotDir_CrossSection = output + plotOut +"CrossSection/";
}

#endif 
