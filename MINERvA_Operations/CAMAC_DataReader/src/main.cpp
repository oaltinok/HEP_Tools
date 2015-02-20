/*
================================================================================
main.cpp
    Main Function that controls the TOF and Veto Scalars
    
    Classes Used:
        CAMAC_DataReader
    
    CAMAC_DataReader uses other objects to read & process CAMAC Data
        See CAMAC_DataReader Header for more information
   
    Compile:
        > make
        
    Usage:
        > ./CAMACDataReader
    
    Author:        Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision: 2015_02_20
================================================================================
*/
#include "../CAMAC_DataReader/CAMAC_DataReader.h"

using namespace std;

int main(int argc, char* argv[])
{
    
    CAMAC_DataReader dataReader;
    
    dataReader.ReadData();
   
    return 0;
}






