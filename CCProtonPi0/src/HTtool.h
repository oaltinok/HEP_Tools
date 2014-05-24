/*
    HTtool.h Duplicated from CCPi0 Package on 2014-05-24
        Purpose: Make CCProtonPi0 Package independent of CCPi0 Package
        Future: Common Tools and Functions will be combined under AnaUtils or
                PionUtils
                
    Original Author:    Trung Le
    Author:             Ozgur Altinok  - ozgur.altinok@tufts.edu
    Date:               2014_05_24
    Last Revision:      2014_05_24
*/
#ifndef HTTOOL_H 
#define HTTOOL_H 1

// inheritance
#include "MinervaUtils/MinervaHistoTool.h"
#include "CCProtonPi0/IHoughTool.h"


//forwards
class TH2D;


class HTtool: public MinervaHistoTool, virtual public IHoughTool {

  public:

    HTtool( const std::string& type, const std::string& name, const IInterface* parent ); ///< Standard constructor
    virtual ~HTtool(); ///< Destructor

    StatusCode initialize();
    StatusCode finalize();
    
    StatusCode GetReference( SmartRefVector<Minerva::IDCluster> idClusterView, Gaudi::XYZPoint &ref ) const;
    StatusCode GetReference( BlobSeeds::seedCandVector Seeds, Gaudi::XYZPoint &ref ) const;

    StatusCode Hough2D( SmartRefVector<Minerva::IDCluster> idClusterVec, double &r, double &theta, const Gaudi::XYZPoint& ref ) const;
    StatusCode Hough2D( SmartRefVector<Minerva::IDCluster> idClusterVec, double &r, double &theta, const Gaudi::XYZPoint& ref, const Gaudi::XYZPoint& vert ) const;
    StatusCode Hough2D( BlobSeeds::seedCandVector Seeds, double &r, double &theta, const Gaudi::XYZPoint& ref ) const;
    StatusCode Hough2D( BlobSeeds::seedCandVector Seeds, double &r, double &theta, const Gaudi::XYZPoint& ref, const Gaudi::XYZPoint& vert ) const;
    
    StatusCode FillHough1Cluster( SmartRef<Minerva::IDCluster> idCluster, TH2D *h, const Gaudi::XYZPoint& ref) const;
    StatusCode FillHough1Cluster( BlobSeeds::seedCandidate *Seed, TH2D *h, const Gaudi::XYZPoint& ref) const;
    
  private:
			
    double a;

};
#endif
