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
Last Revision:  2014_11_11
================================================================================
*/

#ifndef Folder_List_h
#define Folder_List_h

#include <string>

namespace Folder_List
{

    // Default Folder for Output
    const std::string output = "Output/";
    
    // SubFolders for Output
    const std::string rootOut = "RootFiles/";
    const std::string textOut = "TextFiles/";
    const std::string plotOut = "Plots/";
    
    // Folder Branches
    const std::string signal = "Signal/";
    const std::string background = "Background/";
    const std::string allEvents = "AllEvents/";
    
    // Other
    const std::string other = "Other/";

}

#endif 
