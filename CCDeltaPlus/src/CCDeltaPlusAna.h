#ifndef CCDELTAPLUSANA_H 
#define CCDELTAPLUSANA_H 1

#include <utility>

// ineritance
#include "AnaUtils/MinervaAnalysisTool.h"


// forwards
class IMuonUtils;
class IExtraEnergyTool;
class IIDAnchoredBlobCreator;
class ITrackLinearPropagator;
class IHitTaggerTool;
class IHoughBlob;
class IHoughTool;

class IClusterUtilsTool;
class IIDBlobCreator;
class IODBlobCreator;
class IBlobCreatorUtils;
class IIDIsolatedIDBlobCreator;
class IIDBlobSeedingTool;
class IGeomUtilSvc;
class IMinervaCoordSysTool;
class IMinervaMathTool;

namespace Minerva {
  class DeDetector;
  class DeOuterDetector;
}

//! This class is for Reconstruct Pi0 using muon match vertex
class CCDeltaPlusAna : public MinervaAnalysisTool
{
 public:

  //! Standard constructor
  CCDeltaPlusAna( const std::string& type, const std::string& name, const IInterface* parent );

  //! Destructor (mandatory for inheritance)
  virtual ~CCDeltaPlusAna(){};

  //! Initialize (mandatory for inheritance)
  StatusCode initialize();
  //! Finalize (mandatory for inheritance)
  StatusCode finalize();

  //! Reconstruct the event (mandatory for inheritance)
  StatusCode reconstructEvent( Minerva::PhysicsEvent* event, Minerva::GenMinInteraction* truth = NULL ) const;

  //! Attach an interpretations to the event (mandatory for inheritance)
  StatusCode interpretEvent( const Minerva::PhysicsEvent* event, const Minerva::GenMinInteraction* truth,
                             std::vector<Minerva::NeutrinoInt*>& nuInts ) const;

  StatusCode tagTruth( Minerva::GenMinInteraction* truth ) const;


 private:

  IExtraEnergyTool*       m_extraEnergyTool;
  IIDAnchoredBlobCreator* m_idConeScanBlob;
  ITrackLinearPropagator* m_trackPropagator;
  IBlobCreatorUtils*      m_blobUtils;
  IHitTaggerTool* 	  m_colorTag;
  IHoughBlob* 		  m_idHoughBlob;
  IHoughTool* 		  m_idHoughTool;
  IIDBlobSeedingTool*     m_blobSeedingTool;
  IMinervaCoordSysTool*   m_minervaCoordSysTool;
  IMinervaMathTool*       m_mathTool;
  
  const Minerva::DeDetector*      m_idDet;                ///< Inner detector
  const Minerva::DeOuterDetector* m_odDet;                ///< Outer detector
  IGeomUtilSvc*                   m_GeomUtilSvc;          ///< GeomUtilSvc
  
  double m_fiducialApothem;      ///< Vertex Cut (mm) : Is this event signal?
  double m_fiducialUpstreamZ;    ///< Vertex Cut (mm) : Is this event signal?
  double m_fiducialDownstreamZ;  ///< Vertex Cut (mm) : Is this event signal?
  double m_qOverpChargeCut;      ///< q/p charge cut
  double m_energyHoughlimit;     ///< Energy limit to start using Hough T.
  double m_rejectedClustersTime; ///< window time  to allow clusters
  bool   new_impl_;
  bool   fProcessOneEvent;
  int    fRunNumber;
  int    fSubrunNumber;
  int    fGateNumber;

  bool   fAllowUVMatchWithMoreTolerance;
  double fUVMatchTolerance;
  double fUVMatchMoreTolerance;

  
      //dEdx
  int    m_planesdEdx;
  double m_distance_dEdx;
  
      // vertex blob
  bool 	 m_sphereVertex;
  double m_maxSearchD;
  double m_maxStartingD;
  double m_maxSearchGap;		
  bool	 m_filamentVertex;
  double m_maxSearchDFila;
  double m_maxStartingDFila;
  double m_maxSearchGapFila;
  bool   m_filterClusterTypes;
  
  double m_extraEnergyCylinderRadius; 		///< Cylinder Cut (mm) 
  double m_extraEnergyCylinderUpstreamLength;  	///< Cylinder Cut (mm) 
  double m_extraEnergyCylinderDownstreamLength;	///< Cylinder Cut (mm) 
  double m_extraEnergyLowerTimeWindow;    	///< Cylinder Cut (ns) 
  double m_extraEnergyUpperTimeWindow; 		///< Cylinder Cut (ns) 
  double m_extraEnergyPECut; 			///< Cylinder Cut (MeV)?

  mutable std::map<int, Minerva::TG4Trajectory*> fTrajectoryMap;
  mutable const Minerva::TG4Trajectory* fPizero;
  mutable const Minerva::TG4Trajectory* fGamma1;
  mutable const Minerva::TG4Trajectory* fGamma2;

  StatusCode 	fillPi0Branches(Minerva::PhysicsEvent *event,
                                Minerva::IDBlob* idblob1, Minerva::IDBlob* idblob2,
                                const SmartRef<Minerva::Vertex>& vertex ) const;
  
  StatusCode 	getMCPi0Info(Minerva::GenMinInteraction* truth) const;
  bool 	 	ispi0( Minerva::GenMinInteraction* truth, int charge) const;
  bool          is1pi0 (const Minerva::GenMinInteraction* v) const;
  bool 		ispi0x( Minerva::GenMinInteraction* truth) const;
  bool          shouldAnalyzeMC( const Minerva::GenMinInteraction *truth ) const;
  std::pair<bool,bool> ispi0secondary() const; /// the second member is true if it is produced by pi-
  
  
  StatusCode 	PreFilter(Minerva::PhysicsEvent *event ) const;
  
  StatusCode 	DoMuon(Minerva::PhysicsEvent *event, Minerva::GenMinInteraction* truth,
                       SmartRef<Minerva::Vertex>& vertex ) const;
  
  StatusCode 	DoVertex(Minerva::PhysicsEvent *event,  const SmartRef<Minerva::Vertex>& vertex) const;
  
  StatusCode 	VtxBlob(Minerva::PhysicsEvent *event,   const SmartRef<Minerva::Vertex>& vertex ) const;
  
  StatusCode 	ConeBlobs(Minerva::PhysicsEvent *event, const SmartRef<Minerva::Vertex>& vertex ) const;
  
  StatusCode 	DispersedBlob(Minerva::PhysicsEvent *event ) const;
  
  StatusCode  AngleScanBlob(SmartRefVector<Minerva::IDCluster> idClusters,
                             const SmartRef<Minerva::Vertex>& vertex,
                             std::vector<Minerva::IDBlob*>& outBlobs) const;
  
  StatusCode  HoughBlob(SmartRefVector<Minerva::IDCluster> idClusters,
                         const SmartRef<Minerva::Vertex>& vertex,
                         std::vector<Minerva::IDBlob*>& outBlobs) const;
  
  StatusCode  processBlobs( Minerva::PhysicsEvent *event, std::vector<Minerva::IDBlob*> idBlobs) const;

  StatusCode  ODActivity( Minerva::PhysicsEvent *event, std::vector<Minerva::IDBlob*> idBlobs ) const;

  void CheckBlobDigitAssignment(Minerva::PhysicsEvent* event, const Minerva::IDBlob* blob, int index) const;

  void FillTrajectoryMap() const;
  void SummarizeTruth(Minerva::PhysicsEvent* event) const;
  void SummarizeHadronTruth(Minerva::PhysicsEvent* event) const;
  void GammaEdepInfo(Minerva::PhysicsEvent* event) const;

  double CalcTrackLength(const SmartRef<Minerva::Track>& track) const;
  double CalcMinBlobSeparation(const Minerva::IDBlob* blob, const SmartRef<Minerva::Vertex>& vertex) const;
  
  SmartRefVector<Minerva::IDCluster> FilterInSphereClusters(const SmartRef<Minerva::Vertex>& vertex,
                                                            const SmartRefVector<Minerva::IDCluster>& clusters,
                                                            const double radius) const;
  
  std::pair<int,double>
      OneParLineFitBlob(const Minerva::IDBlob* blob, 
                        const SmartRef<Minerva::Vertex>& vertex,
                        const SmartRef<Minerva::Track>& muonTrack) const;
  
  double CalcDistanceFromBlobAxisToVertex(const Minerva::IDBlob* blob,
                                          const SmartRef<Minerva::Vertex>& vertex) const;

  double CalcDistanceFromVertexToExiting(const Minerva::IDBlob* blob,
                                         const SmartRef<Minerva::Vertex>& vertex) const;

  void PrintDigitVectorInfo(const SmartRefVector<Minerva::IDDigit>& digits) const;
  
};

#endif

