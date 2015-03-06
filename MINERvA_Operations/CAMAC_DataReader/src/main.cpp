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

void AutoMode();
void ManualMode(char* argv[]);
void ErrorUsage(char* argv[]);
string GetLatestDataFile();
string GetSingleDataFile(string run, string subrun);
void CheckUsage(int argc, char* argv[]);

int main(int argc, char* argv[])
{	
	CheckUsage(argc,argv);

    // Auto Mode - Automatically reads the latest camacreadout file
    if (argc == 1) AutoMode();
    else ManualMode(argv);

    return 0;
}

void CheckUsage(int argc, char* argv[])
{	
	string keyRun = "-r";
	string keySubrun = "-s";
	
	// Check Number of Arguments
	if( argc != 1 && argc != 5) ErrorUsage(argv);
	
	// Check Keys
	if(argc == 5){
		if(keyRun.compare(argv[1]) != 0) ErrorUsage(argv);
		if(keySubrun.compare(argv[3]) != 0) ErrorUsage(argv);
	}	
}


void ManualMode(char* argv[])
{
	CAMAC_Data Data;
	string dataDir;
	
	dataDir = GetSingleDataFile(argv[2],argv[4]);
	
	if( dataDir != ""){ 
		Data.SetDataFileDir(dataDir);
		Data.ProcessData();
	}else{
		cout<<"ERROR: No File for run = "<<argv[2]<<" subrun = "<<argv[4]<<endl;
		exit(EXIT_FAILURE);
	}
}

void AutoMode()
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
        if(isNewDataFile) Data.ProcessData();
    }
}


void ErrorUsage(char* argv[])
{
	cout<<"\n\n";
	cout<<"Error on Usage! Correct Usage as follows:"<<endl;
	cout<<"Auto Mode:"<<endl;
	cout<<"\t"<<argv[0]<<endl;
	cout<<"Manual Mode:"<<endl;
	cout<<"\t"<<argv[0]<<" -r <run_number> -s <subrun_number>"<<endl;
	cout<<"\n\n";
	exit(EXIT_FAILURE);
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

string GetSingleDataFile(string run, string subrun)
{
    ifstream commandFile;
    string fileName;
    string command = "ls -t /work/data/TB*" + run + "_*" + subrun + "_cosmc" + "*_camac.dat" + " >>command.dat";
    
    // Use system command and write output to command.dat
    system(command.c_str());
    
    // Open Command.dat File
    commandFile.open("command.dat");
    
    // Read Command.dat File
    //      First Line is the requested file
    getline(commandFile,fileName);
    
    // Remove command.dat
    system("rm -f command.dat");
    
    return fileName;
}







