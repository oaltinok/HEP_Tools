/*
    See ClusterVectorInfo.h for Class Information
*/
#include <GaudiKernel/SmartRef.h>

#include <Event/IDCluster.h>

#include "ClusterVectorInfo.h"


ClusterVectorInfo::ClusterVectorInfo(const SmartRefVector<Minerva::IDCluster>& clusters,
                                     bool unusedOnly,
                                     bool noLowActivity,
                                     bool noLowLight)
    : fAllClusters(clusters),
      fDigitCountX(0),
      fDigitCountU(0),
      fDigitCountV(0),
      fUnusedOnly(unusedOnly),
      fNoLowActivity(noLowActivity),
      fNoLowLight(noLowLight)
      
{
    
    for (   SmartRefVector<Minerva::IDCluster>::const_iterator cluster = fAllClusters.begin();
            cluster != fAllClusters.end(); 
            ++cluster )
    {
        // Checks
        if ( fUnusedOnly    && (*cluster)->history() != Minerva::IDCluster::Unused ) continue;
        if ( fNoLowActivity && (*cluster)->type() == Minerva::IDCluster::LowActivity  ) continue;
        if ( fNoLowLight    && (*cluster)->pe()/(*cluster)->iddigs() <= 3.0 ) continue;
        
        // Set Clusters in X View
        if ( (*cluster)->view() == Minerva::IDCluster::X ) {
            fDigitCountX += (*cluster)->iddigs();
            fXClusters.push_back(*cluster);
        }
        // Set Clusters in U View
        if ( (*cluster)->view() == Minerva::IDCluster::U ) {
            fDigitCountU += (*cluster)->iddigs();
            fUClusters.push_back(*cluster);
        }
        // Set Clusters in V View
        if ( (*cluster)->view() == Minerva::IDCluster::V ) {
            fDigitCountV += (*cluster)->iddigs();
            fVClusters.push_back(*cluster);
        }
    }
}

//==============================================================================
// Get Functions for Number of Clusters
//==============================================================================
// All Clusters
unsigned int ClusterVectorInfo::GetN() const 
{
    return GetNx() + GetNu() + GetNv();
}

// X View
unsigned int ClusterVectorInfo::GetNx() const 
{
    return fXClusters.size();
}

// U View
unsigned int ClusterVectorInfo::GetNu() const 
{
    return fUClusters.size();
}

// V View
unsigned int ClusterVectorInfo::GetNv() const 
{
    return fVClusters.size();
}

//==============================================================================
// Get Functions for Digit Count
//==============================================================================
// All Digits
unsigned int ClusterVectorInfo::GetDigitCount() const 
{
    return fDigitCountX + fDigitCountU + fDigitCountV;
}

// X View
unsigned int ClusterVectorInfo::GetDigitCountX() const 
{
    return fDigitCountX;
}

// U View
unsigned int ClusterVectorInfo::GetDigitCountU() const 
{
    return fDigitCountU;
}

// V View
unsigned int ClusterVectorInfo::GetDigitCountV() const 
{
    return fDigitCountV;
}

//==============================================================================
// Get Functions for Cluster Vectors in Specific Views
//==============================================================================
// X View
const SmartRefVector<Minerva::IDCluster>& ClusterVectorInfo::GetXClusters() const 
{
    return fXClusters;
}

// U View
const SmartRefVector<Minerva::IDCluster>& ClusterVectorInfo::GetUClusters() const 
{
    return fUClusters;
}

// V View
const SmartRefVector<Minerva::IDCluster>& ClusterVectorInfo::GetVClusters() const 
{
    return fVClusters;
}
