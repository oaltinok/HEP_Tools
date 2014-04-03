/*
================================================================================
CCDeltaPlusAna

    Reconstruction Package:
        Exclusive Channel for muon,proton,pi0 on the final state
        Contains Clone of CCPi0 for Pi0 Reconstruction 
        Contains Clone of NukeCCQE for Short Track Reconstruction 
        Both Clones may be inherited in future, for now they are just copied
    
    Main Package:
        AnalysisFramework/Ana/CCDeltaPlus/
        
    Setup:
        > getpack -u Ana/CCDeltaPlus
        > cmt config
        > cmt make
        
    Usage:
        There is a corresponding Options file under Tools/SystemTests
        Use Tools/ProductionScripts/ana_scripts/ProcessAna.py

    
    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision:  2014_04_03
    
================================================================================
*/
#ifndef CCDELTAPLUSANA_H 
#define CCDELTAPLUSANA_H 1

#include <utility>

// ineritance
#include "AnaUtils/MinervaAnalysisTool.h"

//-- Forward Declarations
#include "Event/MinervaEventFwd.h"


//=============================================================================
// NukeCCQE Classes
//=============================================================================
class IConeUtilsTool;
class Cone;
class ParticleExtraDataDefs;
class IEnergyLoss;
class IEnergyCorrectionTool;
class IAbsorberStacker;
class IMinervaCoordSysTool;
class INuclearTargetTool;
class NuclearTarget;
class IAnchoredTrackFormation;
class IPrimaryBlobProngTool;
class ITrackPropagator;
class IMinosMinervaTransformer;
class IVertexFitter;
class IProngClassificationTool;
class IODProngClassificationTool;
class IParticleTool;
class IParticleMakerTool; 
class IBlobCreatorUtils; 
class IVertexEnergyStudyTool;
class IMichelTool; 

class IHitTaggerTool;
class MinervaObjectSort;

class IProtonUtils;
class IMCTrackTool;
class MCTrack;

//=============================================================================
// CCPi0 Classes
//=============================================================================
class IMuonUtils;
class IExtraEnergyTool;
class IIDAnchoredBlobCreator;
class ITrackLinearPropagator;
class IHoughBlob;
class IHoughTool;

class IClusterUtilsTool;
class IIDBlobCreator;
class IODBlobCreator;
class IIDIsolatedIDBlobCreator;
class IIDBlobSeedingTool;
class IGeomUtilSvc;
class IMinervaMathTool;



namespace Minerva {
  class TimeSlice;
  class DePlane;
  class DeDetector;
  class DeOuterDetector;
  class NuclearTarget;
}

//! This class is for Reconstruct Pi0 using muon match vertex
class CCDeltaPlusAna : public MinervaAnalysisTool
{
 private:
   typedef std::vector<Minerva::NeutrinoInt*> NeutrinoVect;

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
  StatusCode interpretEvent( const Minerva::PhysicsEvent* event, const Minerva::GenMinInteraction* truth, NeutrinoVect& nuInts ) const;
  StatusCode tammy_interpretEvent( const Minerva::PhysicsEvent *event, const Minerva::GenMinInteraction *truth, NeutrinoVect& nuInts ) const;

  StatusCode tagTruth( Minerva::GenMinInteraction* truth ) const;
  


 private:
 
//=============================================================================
// NukeCCQE Private Functions
//=============================================================================
 
  StatusCode interpretFailEvent( Minerva::PhysicsEvent* event ) const;
  StatusCode makeConeVtxBlobProng( SmartRef<Minerva::Vertex>& vertex, std::vector<Minerva::Prong*>* vtxBlobProngVect ) const;

  Minerva::IDClusterVect      getClusters( Minerva::PhysicsEvent* event, Minerva::Track::PatRecHistory pattern ) const;
  Minerva::Vertex*            getStopPointVertex( Minerva::PhysicsEvent* event, SmartRef<Minerva::Track>& track ) const;
  SmartRef<Minerva::Particle> findMuonParticleCandidate( SmartRef<Minerva::Prong>& prong ) const;

  std::string odMuonProng( SmartRef<Minerva::Prong>& prong, SmartRef<Minerva::Particle>& particle ) const;

  bool isInsideTargetFiducial( SmartRef<Minerva::Vertex>& vertex, std::string& name, int& passVertexZCut ) const;
  
  bool protonProng( Minerva::ProngVect& prongs, SmartRef<Minerva::Prong>& prong, SmartRef<Minerva::Particle>& particle ) const;
  bool muonProng( Minerva::ProngVect& prongs, SmartRef<Minerva::Prong>& prong, SmartRef<Minerva::Particle>& particle, 
                  std::string& event_type ) const;

  bool createdAnchoredShortTracks( Minerva::PhysicsEvent* event, Minerva::Vertex* vertex ) const;
  bool createdKinkedShortTracks( Minerva::PhysicsEvent* event ) const;

  bool createdTrackedProngs( Minerva::PhysicsEvent* event ) const;
  bool createdTrackedParticles( Minerva::ProngVect& prongs ) const;

  bool createdIDBlobProng( Minerva::PhysicsEvent* event ) const;
  bool createdODBlobProng( Minerva::PhysicsEvent* event ) const;

  bool classifyProngs( Minerva::PhysicsEvent* event ) const;

  void correctProtonProngEnergy( SmartRef<Minerva::Prong>& protonProng, double& p_calCorrection, double& p_visEnergyCorrection ) const;
  void correctProtonProngEnergy( SmartRef<Minerva::Prong>& protonProng, SmartRef<Minerva::Particle> particle ) const {
       double p_calCorrection = 0, p_visEnergyCorrection = 0;
       correctProtonProngEnergy(protonProng,p_calCorrection,p_visEnergyCorrection);
       double E  = particle->momentumVec().E() + p_calCorrection;
       double p  = sqrt( E*E - particle->mass()*particle->mass() );
       double px = p*sin((protonProng->minervaTracks().front())->theta())*cos((protonProng->minervaTracks().front())->phi());
       double py = p*sin((protonProng->minervaTracks().front())->theta())*sin((protonProng->minervaTracks().front())->phi());
       double pz = p*cos((protonProng->minervaTracks().front())->theta());
       Gaudi::XYZTVector updateFourMomentum(px,py,pz,E);
       particle->setMomentumVec(updateFourMomentum);
       return;
  }

  void linkProngToPrimaryProng( Minerva::Prong* prong, Minerva::ProngVect& primaryProngs ) const;
  void setTrackDirection( Minerva::Track* track, Minerva::Vertex* vertex ) const;
  void setMuonID( Minerva::NeutrinoInt* neutrino, SmartRef<Minerva::Prong>& prong ) const;
  void setODMuonParticleData( Minerva::NeutrinoInt* neutrino, SmartRef<Minerva::Particle>& particle ) const;
  void setTrackProngTruth( Minerva::NeutrinoInt* neutrino, SmartRef<Minerva::Prong>& prong ) const;
  void setBlobProngTruth( Minerva::NeutrinoInt* neutrino, Minerva::ProngVect& prongs ) const;
  void setProtonParticleData( Minerva::NeutrinoInt* nu, SmartRef<Minerva::Prong>& prong, SmartRef<Minerva::Particle>& particle,
                              double vertex_z ) const;
  void setPionParticleData( Minerva::NeutrinoInt* neutrino, SmartRef<Minerva::Prong>& prong ) const;
  void setGenMinTruthInformation( Minerva::PhysicsEvent* event, Minerva::GenMinInteraction* truth ) const;
  void fillBlobEventData( Minerva::PhysicsEvent* event ) const;
  void fillMuonMomentumFromProton( Minerva::NeutrinoInt* nu, Minerva::TrackVect& tracks, double proton_theta, double proton_p ) const;

  bool tagMichelElectrons(Minerva::PhysicsEvent* event) const;

  IConeUtilsTool*             m_coneUtilsTool;
  IEnergyCorrectionTool*      m_energyCorrectionTool;
  IEnergyLoss*                m_energyLoss;
  IAbsorberStacker*           m_absorberStacker;
  IMinervaCoordSysTool*       m_coordSysTool;
  ITrackPropagator*           m_propagateToMinos;
  IMinosMinervaTransformer*   m_mmt;
  IMCTrackTool*               m_mcTrackTool;
  IHitTaggerTool*             m_hitTagger;
  Minerva::DeOuterDetector*   m_odDet;
  Minerva::DeDetector*        m_idDet;

  INuclearTargetTool*         m_nuclearTargetTool;
  std::string                 m_nuclearTargetToolAlias;

  IProtonUtils*               m_protonUtils;
  std::string                 m_protonUtilsAlias;

  IPrimaryBlobProngTool*      m_primaryBlobProngTool;
  std::string                 m_primaryBlobProngAlias;

  IAnchoredTrackFormation*    m_anchoredTracker;
  std::string                 m_anchorShortTrackerAlias;

  IVertexFitter*              m_vertexFitter;
  std::string                 m_vertexFitterAlias;

  IProngClassificationTool*   m_prongIntersection;
  std::string                 m_prongIntersectionAlias;

  IODProngClassificationTool* m_odMatchTool;
  std::string                 m_odMatchAlias;

  IParticleMakerTool*         m_particleMaker;
  std::string                 m_particleMakerAlias; 

  IBlobCreatorUtils*          m_blobCreatorUtils; 
  std::string                 m_blobCreatorUtilsAlias; 

  IVertexEnergyStudyTool*     m_vertexEnergyStudyTool; 
  std::string                 m_vtxEngStudyToolAlias;

  IMichelTool*                m_michelTool;
  std::string                 m_michelToolAlias;

  double                      m_fvApothem;
  double                      m_min_mva_score;
  double                      m_CCQEBindingEnergyMeV;
  double                      m_clusMinTimeWindow;
  double                      m_clusMaxTimeWindow;
  double                      m_fuzzRadius;

  int                         m_numSearchRadii;
  double                      m_searchStepSize;
  double                      m_maxSearchDistance;
  double 		      m_maxStartingDistance;
  double 		      m_maxAllowedSearchGap;
  double 		      m_maxSeparationBlobVertex;

  bool                        m_applyCalCorrection;

  double                      m_tar1VertexZUSCut;
  double                      m_tar1VertexZDSCut;

  double                      m_tar2VertexZUSCut;
  double                      m_tar2VertexZDSCut;

  double 		      m_tar3VertexZUSCut;
  double                      m_tar3VertexZDSCut;

  double                      m_tar4VertexZUSCut;
  double                      m_tar4VertexZDSCut;

  double                      m_tar5VertexZUSCut;
  double                      m_tar5VertexZDSCut;

  double                      m_maxSearchDistance_VESTool;
  double                      m_maxStartingDistance_VESTool;
  double                      m_maxAllowedSearchGap_VESTool;

  double                      m_michel_downstreamZ;
  double                      m_michel_upstreamZ;

  std::vector<double>                   m_fvUpstreamZ;
  std::vector<double>                   m_fvDownstreamZ;
  std::vector<std::string>              m_targetNames;
  std::vector<Minerva::NuclearTarget* > m_targets;

  int     m_muonProngColor; 
  int     m_protonProngColor; 
  int     m_primaryVertexColor; 
  int     m_secondaryVertexColor; 
  int     m_endPointVertexColor; 
  int     m_unattachedProngColor;
 

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

