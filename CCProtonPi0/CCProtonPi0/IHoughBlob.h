/*
    IHoughBlob.h Duplicated from CCPi0 Package on 2014-05-24
        Purpose: Make CCProtonPi0 Package independent of CCPi0 Package
        Future: Common Tools and Functions will be combined under AnaUtils or
                PionUtils
                
    Original Author:    Trung Le
    Author:             Ozgur Altinok  - ozgur.altinok@tufts.edu
    Date:               2014_05_24
    Last Revision:      2014_05_24
*/
#ifndef IHOUGHBLOB_H 
#define IHOUGHBLOB_H 1

// Include files
// -------------
// inheritance
#include "GaudiKernel/IAlgTool.h"

// forwards
#include "Event/VectorTypeDefs.h" //includes MinervaEventFwd.h

static const InterfaceID IID_IHoughBlob( "IHoughBlob", 1, 0 );

class IHoughBlob : virtual public IAlgTool {

  public:

    /// Return the interface ID
    static const InterfaceID& interfaceID() { return IID_IHoughBlob; }

    virtual StatusCode initialize() = 0;
    virtual StatusCode finalize()   = 0;

    virtual StatusCode GetViewClustersAboveThreshold( SmartRefVector<Minerva::IDCluster> &idClusterVec,
                                                      SmartRefVector<Minerva::IDCluster> &idClusterView, 
                                                      Minerva::IDCluster::View view, double pecut ) const = 0;

    virtual StatusCode PseudoCone(SmartRefVector<Minerva::IDCluster> &Seed, SmartRefVector<Minerva::IDCluster> &ClusVectorX,
                                  Gaudi::XYZVector direction, Gaudi::XYZPoint vert ) const =0;

    virtual StatusCode XUVMatch(SmartRefVector<Minerva::IDCluster> &Seed, SmartRefVector<Minerva::IDCluster> &ClusVectorU,
                                SmartRefVector<Minerva::IDCluster> &ClusVectorV, double match ) const = 0;
    
    virtual StatusCode XUVMatch(SmartRef<Minerva::IDCluster> Cluster, SmartRefVector<Minerva::IDCluster> &Seed, 
                                SmartRefVector<Minerva::IDCluster> &ClusVectorU,	SmartRefVector<Minerva::IDCluster> &ClusVectorV, 
                                double zmin, double zmax, double match) const = 0;
    
    virtual StatusCode AddClusterInsideCone(SmartRef<Minerva::IDCluster> UnuCluster, std::vector<Minerva::IDBlob*> &idBlobs,
                                            Gaudi::XYZPoint vert) const = 0;
    
    virtual StatusCode Angle( SmartRef<Minerva::IDCluster> Cluster, Gaudi::XYZVector direction, Gaudi::XYZPoint vert, double &radius ) const = 0;
    
    virtual StatusCode Create2dHTSeed( SmartRefVector<Minerva::IDCluster> &idClusterView, SmartRefVector<Minerva::IDCluster> &HT2dClusters, 
                                       double r, double theta, Gaudi::XYZPoint ref, double &spX, double &spZ ) const = 0;	
    
    virtual StatusCode GetDirection( Minerva::IDBlob *idBlob ) const = 0;
    
    virtual bool  GetDirection( Minerva::IDBlob *idBlob, Gaudi::XYZPoint vertex ) const = 0;
    
    virtual bool GetStartPosition( Minerva::IDBlob *idBlob, Gaudi::XYZPoint vert, bool is_vertex = false ) const = 0;
    
    virtual StatusCode idBlobdEdx( Minerva::IDBlob *idblob, double &dEdx ) const = 0;
    
    virtual StatusCode isPhoton( SmartRefVector<Minerva::IDCluster> Seed, Gaudi::XYZPoint vtX ) const = 0;
    
    virtual StatusCode getBlobEnergyTime( Minerva::IDBlob *idblob, double &energy,
                                          double& tracker_evis,
                                          double& ecal_evis,
                                          double& hcal_evis,
                                          double& scal_evis) const = 0;
    
};
#endif
