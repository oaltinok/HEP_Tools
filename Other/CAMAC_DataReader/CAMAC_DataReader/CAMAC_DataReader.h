/*
================================================================================
Class: CAMAC_DataReader
    Controls CAMAC Data Objects:
        TOF
        Veto
    Assigns Data Locations for the CAMAC Data Objects
    Process Data for each CAMAC Data Object
            
    
    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision:  2015_02_20
================================================================================
*/

#ifndef CAMAC_DataReader_h
#define CAMAC_DataReader_h

#include <iostream>

#include "CAMAC_Data.h"

class CAMAC_DataReader
{
    public:
        CAMAC_DataReader();
        ~CAMAC_DataReader();
        
        void ReadData();

    private:
        CAMAC_Data TOF;
        CAMAC_Data Veto;
        
        // I/O Files
        std::string rootDir;
        std::string dataDir_TOF;
        std::string dataDir_Veto;
        std::string branchName_TOF;
        std::string branchName_Veto;   
};


#endif
