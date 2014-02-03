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

    
    Last Revision: 2014_02_03
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
 File Names --Starts with F --
--------------------------------------------------------------------------------
*/

// Playlists
const string F_PL_TEST = "pl_MC_Test.dat";
const string F_PL_FULL = "pl_MC_Brandon.dat";

// ROOT Files
const string F_ROOT_MC = "MC_minerva1.root";
const string F_ROOT_DATA = "DATA_minerva1.root";

// TEXT Files
const string F_TEXT_CUT = "CutTable.txt";
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

// Input Folder Locations    
const string IN = "Input/";
const string PLAYLIST = "Playlists/";


#endif