/*
    See CAMAC_DataReader.h header or Class Information
*/
#ifndef CAMAC_DataReader_cpp
#define CAMAC_DataReader_cpp

#include "../CAMAC_DataReader/CAMAC_DataReader.h"

using namespace std;

CAMAC_DataReader::CAMAC_DataReader()
{
    //     rootDir = "DataFiles/Test.root";
    rootDir = "DataFiles/NearlineCurrentHistos.root";
    dataDir_TOF = "DataFiles/testreadoutcamacv6_camac.dat";
    dataDir_Veto = "DataFiles/testreadoutcamacv6_camac.dat";
    branchName_TOF = "TOF_Data";
    branchName_Veto = "Veto_Data"; 
    
    TOF.SetRootFileDir(rootDir);
    TOF.SetDataFileDir(dataDir_TOF);
    TOF.SetBranchName(branchName_TOF);
    
    Veto.SetRootFileDir(rootDir);
    Veto.SetDataFileDir(dataDir_Veto);
    Veto.SetBranchName(branchName_Veto);
        
}


void CAMAC_DataReader::ReadData()
{
    cout<<"Reading CAMAC TOF Data"<<endl;
    TOF.ProcessData();
    
    cout<<"Reading CAMAC Veto Data"<<endl;
    Veto.ProcessData();
}


CAMAC_DataReader::~CAMAC_DataReader()
{

}

#endif
