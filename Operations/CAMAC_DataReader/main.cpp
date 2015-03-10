/*
================================================================================
main.cpp
    Main Function for package
    
    Classes Used:
        CAMAC_DataReader - See Class Header for more Information
    
    Compile:
        > make
        
    Usage:
        > ./CAMACDataReader
        or
        > ./CAMACDataReader -r <run_number> -s <run_number>
    
    Author:        Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision: 2015_03_10
================================================================================
*/
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <string>

// Include All Classes
#include "Classes/CAMAC_DataReader/CAMAC_DataReader.h"

using namespace std;

void AutoMode();
void ManualMode(char* argv[]);
void ErrorUsage(char* argv[]);
void CheckUsage(int argc, char* argv[]);

int main(int argc, char* argv[])
{
    CheckUsage(argc,argv);
    
    if (argc == 1) AutoMode();
    else ManualMode(argv);

    return 0;
}


void ManualMode(char* argv[])
{
    string run = argv[2];
    string subrun = argv[4];
        
    CAMAC_DataReader CDR("Manual", run, subrun);
    
    CDR.RunManual();
}

void AutoMode()
{
	cout<<" Auto Mode is not working yet!"<<endl;
	exit(EXIT_FAILURE);
    CAMAC_DataReader CDR("Auto");

    while(true){
        CDR.RunAuto();
    }   
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
