#ifndef HTBLOB_H
#define HTBLOB_H 1

// inheritance
#include "MinervaUtils/MinervaHistoTool.h"
#include "CCDeltaPlus/IHoughBlob.h"

// forwards
class IMinervaMathTool;
class TH1D;

class HTBlob: public MinervaHistoTool, virtual public IHoughBlob {

  public:

    // constructor
    HTBlob( const std::string& type, const std::string& name, const IInterface* parent );

    // destructor
    virtual ~HTBlob();

    virtual StatusCode initialize();    //< Algorithm initialization
    virtual StatusCode finalize();    //< Algorithm finalization

    StatusCode GetViewClustersAboveThreshold( SmartRefVector<Minerva::IDCluster> &idClusterVec,
                                              SmartRefVector<Minerva::IDCluster> &idClusterView, 
                                              Minerva::IDCluster::View view, double pe_cut ) const;

    StatusCode PseudoCone(SmartRefVector<Minerva::IDCluster> &Seed, SmartRefVector<Minerva::IDCluster> &ClusVectorX,
                          Gaudi::XYZVector direction, Gaudi::XYZPoint vert ) const;
    
    StatusCode XUVMatch(SmartRefVector<Minerva::IDCluster> &Seed, SmartRefVector<Minerva::IDCluster> &ClusVectorU,
                        SmartRefVector<Minerva::IDCluster> &ClusVectorV, double match ) const;
    
    StatusCode XUVMatch(SmartRef<Minerva::IDCluster> Cluster, SmartRefVector<Minerva::IDCluster> &Seed, 
                        SmartRefVector<Minerva::IDCluster> &ClusVectorU,	SmartRefVector<Minerva::IDCluster> &ClusVectorV,
                        double zmin, double zmax, double match) const;
    
    StatusCode AddClusterInsideCone(SmartRef<Minerva::IDCluster> UnuCluster, std::vector<Minerva::IDBlob*> &idBlobs, 
                                    Gaudi::XYZPoint vert) const;
    
    StatusCode Angle( SmartRef<Minerva::IDCluster> Cluster, Gaudi::XYZVector direction, Gaudi::XYZPoint vert, double &radius ) const;
    
    StatusCode Create2dHTSeed( SmartRefVector<Minerva::IDCluster> &idClusterView, SmartRefVector<Minerva::IDCluster> &HT2dClusters, 
                               double r, double theta, Gaudi::XYZPoint ref, double &spX, double &spZ) const;
    
    StatusCode GetDirection( Minerva::IDBlob *idBlob ) const;
    
    bool GetDirection( Minerva::IDBlob *idBlob, Gaudi::XYZPoint vertex ) const;
    
    bool GetStartPosition( Minerva::IDBlob *idBlob, Gaudi::XYZPoint vert, bool is_vertex ) const;

    StatusCode idBlobdEdx( Minerva::IDBlob *idblob, double &dEdx ) const;
    
    StatusCode isPhoton( SmartRefVector<Minerva::IDCluster> Seed, Gaudi::XYZPoint vtX ) const;
    
    StatusCode getBlobEnergyTime( Minerva::IDBlob *idblob, double &energy ) const;
    
    StatusCode invariantMass( Minerva::IDBlob* idblob1, Minerva::IDBlob* idblob2, double &mass, Gaudi::XYZPoint vert ) const;

    StatusCode energyCorrection( double &energy, double z, double radius ) const;

  private:

    IMinervaMathTool* m_mathTool;
	
    int m_planesdEdx;
    
    //vertex blob to event vertex
    double m_minDistanceStripPhoton;
    double m_minDistanceModulePhoton;
    
        // calibration constants
    double  m_scalefactor;
    double  m_calibrationTracker;
    double  m_calibrationECal;
    double  m_calibrationHCal;
    
    StatusCode FinddEdxPlanes(TH1D *h, int &index, double &dEdx) const;
    
};

#endif // HBLOB_H
