/*
================================================================================
Library: Folder_List
    Folder list for Input, Output Files
    
    Main Directory:
        Libraries/
    
    Usage:
        > #include "Libraries/Folder_List.h" 
        > Folder_List::IN + Folder_List::PLAYLIST + Folder_List::F_PL_TEST
            // Returns the Location of Input Test Playlist

    
    Last Revision: 2014_03_18
================================================================================
*/

#ifndef Folder_List_h
#define Folder_List_h

#include<string>

using namespace std;

/*
--------------------------------------------------------------------------------
Methods for getting Complete File Location
--------------------------------------------------------------------------------
*/
string getFileLocation(string s1, string s2)
{
    return s1+s2;
}

string getFileLocation(string s1, string s2, string s3)
{
    return s1+s2+s3;
}


/*
--------------------------------------------------------------------------------
 Channel Tag
--------------------------------------------------------------------------------
*/

const string channelTag = "Test";

/*
--------------------------------------------------------------------------------
 File Names --Starts with F --
--------------------------------------------------------------------------------
*/

// Playlists
const string F_PL_TEST = "pl_MC_Test.dat";
const string F_PL_FULL = "pl_MC_minerva1.dat";


// ROOT Files
const string F_ROOT_MC = "MC_minerva1.root";
const string F_ROOT_DATA = "DATA_minerva1.root";
const string F_ROOT_MUON = "Muon.root";
const string F_ROOT_PION = "Pion.root";
const string F_ROOT_PROTON = "Proton.root";
const string F_ROOT_CCDELTAPLUS = "CCDeltaPlus.root";

// TEXT Files
const string F_TEXT_CUT = "CutTable";
const string F_TEXT_README = "readme.txt";


/*
--------------------------------------------------------------------------------
 Folder Locations
--------------------------------------------------------------------------------
*/

// Output Folder Locations
const string OUT = "Output/";
const string ROOT = "RootFiles/";
const string TEXT = "TextFiles/";
const string PLOT = "Plots/";
const string PLOT_MUON = "Plots/Muon/";
const string PLOT_PROTON = "Plots/Proton/";
const string PLOT_PION = "Plots/Pion/";
const string PLOT_CCDELTAPLUS= "Plots/CCDeltaPlus/";

// Input Folder Locations    
const string IN = "Input/";
const string PLAYLIST = "Playlists/";

/*
--------------------------------------------------------------------------------
 Complete Address
--------------------------------------------------------------------------------
*/

// Particle Specific Root Files - Generated Automatically with CCPi0::run()
string f_Root_Muon =    getFileLocation(Folder_List::OUT, Folder_List::ROOT,Folder_List::F_ROOT_MUON);
string f_Root_Proton =  getFileLocation(Folder_List::OUT, Folder_List::ROOT,Folder_List::F_ROOT_PROTON);
string f_Root_Pion =    getFileLocation(Folder_List::OUT, Folder_List::ROOT,Folder_List::F_ROOT_PION);

// Particle Specific Plot Folders - Folders Created with setup.sh
string f_Plot_Muon =    getFileLocation(Folder_List::OUT, Folder_List::PLOT_MUON );
string f_Plot_Proton =  getFileLocation(Folder_List::OUT, Folder_List::PLOT_PROTON);
string f_Plot_Pion =    getFileLocation(Folder_List::OUT, Folder_List::PLOT_PION);

// Default Files and Folders
string f_Root_CCDeltaPlus =   getFileLocation(Folder_List::OUT, Folder_List::ROOT,Folder_List::F_ROOT_CCDELTAPLUS);
string f_Plot_CCDeltaPlus =   getFileLocation(Folder_List::OUT, Folder_List::PLOT_CCDELTAPLUS);

// Playlist File Addresses
string f_PL_Test =      getFileLocation(Folder_List::IN, Folder_List::PLAYLIST, Folder_List::F_PL_TEST);
string f_PL_Full =      getFileLocation(Folder_List::IN, Folder_List::PLAYLIST, Folder_List::F_PL_FULL);




#endif