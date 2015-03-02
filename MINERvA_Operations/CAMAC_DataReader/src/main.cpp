/*
================================================================================
main.cpp
    Main Function that controls reading CAMAC readout file
    
    Classes Used:
        CAMAC_Data - See Class Header for more Information
    
    Compile:
        > make
        
    Usage:
        > ./CAMACDataReader <CAMAC_readout_file>
    
    Author:        Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision: 2015_03_02
================================================================================
*/

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>

#include "../CAMAC_DataReader/CAMAC_Data.h"

using namespace std;

void ErrorUsage();
string GetLatestDataFile();


int main(int argc, char* argv[])
{
    CAMAC_Data Data;
    string latestDataDir;
    bool isNewDataFile = false;
    
    // Run First Time
    latestDataDir = GetLatestDataFile();
    Data.SetDataFileDir(latestDataDir);
    Data.ProcessData();
   
//     while(true){
//         
//     }

    // Check for new data file and re-run again if needed
    for (int i = 0; i < 3; i++){
        // Get Latest Data Dir
        latestDataDir = GetLatestDataFile();
        
        // Check whether the Data File Changed
        isNewDataFile = Data.ChangeDataFileDir(latestDataDir);
        
        // Process Data if the data file is NEW
        if(isNewDataFile){
            Data.ProcessData();
        }
    }

    return 0;
}

string GetLatestDataFile()
{
    ifstream commandFile;
    string dir = "DataFiles/";
    string fileName;
    
    // Use system command and write output to command.dat
    system("ls -t DataFiles/ >> command.dat");
    
    // Open Command.dat File
    commandFile.open("command.dat");
    
    // Read Command.dat File
    //      First Line is the latest modified file
    getline(commandFile,fileName);
    
    // Remove command.dat
    system("rm -f command.dat");
    
    // Append fileName to the Dir
    dir.append(fileName);
    
    return dir;
}






