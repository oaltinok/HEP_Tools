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
    
}

#endif 
