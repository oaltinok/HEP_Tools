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
    Last Revision: 2015_03_03
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
   
    // Check for new data file and re-run again if needed
    while(true){
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

/*
	GetLatestDataFile()
		Gets latest camac readout file from /work/data
*/
string GetLatestDataFile()
{
    ifstream commandFile;
    string fileName;
    
    // Use system command and write output to command.dat
    system("ls -t /work/data/TB*camac.dat >> command.dat");
    
    // Open Command.dat File
    commandFile.open("command.dat");
    
    // Read Command.dat File
    //      First Line is the latest modified file
    getline(commandFile,fileName);
    
    // Remove command.dat
    system("rm -f command.dat");
    
    return fileName;
}






