#ifndef CCProtonPi0_Pi0BlobTool_cpp
#define CCProtonPi0_Pi0BlobTool_cpp

// #include <Event/Vertex.h>
#include <Event/IDCluster.h>
#include <Event/IDBlob.h>

#include "Pi0BlobTool.h"

using namespace std;

Pi0BlobTool::Pi0BlobTool()
{
    // Do Nothing!
}

bool Pi0BlobTool::isBlob1_Good( Minerva::IDBlob* blob )
{
    // Sanity Check
    if (blob == NULL) return false;

    bool isBlobGood =   isNDigitsHigh(blob, minNDigits_blob1 ) &&
                        isEnergyHigh(blob, minEnergy_blob1);

    return isBlobGood;
}

bool Pi0BlobTool::isBlob2_Good( Minerva::IDBlob* blob )
{
    // Sanity Check
    if (blob == NULL) return false;

    bool isBlobGood = isNDigitsHigh(blob, minNDigits_blob2 );

    return isBlobGood;
}

bool Pi0BlobTool::isNDigitsHigh( Minerva::IDBlob* blob, const double minNDigits)
{
    double blob_ndigits = blob->getAllDigits().size();

    if (blob_ndigits > minNDigits ) return true;
    else return false;
}

bool Pi0BlobTool::isEnergyHigh( Minerva::IDBlob* blob, const double minEnergy )
{
    double blob_energy = blob->energy();

    if ( blob_energy > minEnergy) return true;
    else return false;
}

#endif

