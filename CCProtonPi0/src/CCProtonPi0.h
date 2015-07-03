/*
================================================================================
CCProtonPi0

    Reconstruction Package:
        Exclusive Channel for muon,proton,pi0 on the final state
        Uses Nightly Build
    
    Main Package:
        AnalysisFramework/Ana/CCProtonPi0/
        
    Setup:
        > getpack -u Ana/CCProtonPi0
        > cmt config
        > cmt make
        
    Usage:
        There is a corresponding Options file under Tools/SystemTests
        Use Tools/ProductionScripts/ana_scripts/ProcessAna.py

    
    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
    Version:        v2_13
    
================================================================================
*/
#ifndef CCPROTONPI0_H 
#define CCPROTONPI0_H 1

#include <utility>
#include <fstream>

// ineritance
#include "AnaUtils/MinervaAnalysisTool.h"

//-- Forward Declarations
#include "Event/MinervaEventFwd.h"

class IAnchoredTrackFormation;
class IBlobCreatorUtils;
class ICCPionIncUtils;
class ICalorimetryUtils;
class IConeUtilsTool;
class IEnergyCorrectionTool;
class IExtraEnergyTool;
class IGeomUtilSvc;
class IGetCalAttenuation;
class IGetDeadTime;
class IGiGaGeomCnvSvc;
class IHitTaggerTool;
class IHoughBlob;
class IHoughTool;
class IIDAnchoredBlobCreator;
class IMCTrackTool;
class IMichelTool;
class IMinervaCoordSysTool;
class IMinervaMathTool;
class IMinervaObjectAssociator;
class IODProngClassificationTool;
class IParticleMakerTool;
class IParticleTool;
class IProngClassificationTool;
class IProtonUtils;
class IRecoObjectTimeTool;
class IVertexEnergyStudyTool;
class IVertexFitter;
class TRandom3;

namespace Minerva {
  class DeDetector;
  class DeOuterDetector;
}

//! This class is for Reconstruct Pi0 using muon match vertex
class CCProtonPi0 : public MinervaAnalysisTool
{
    private:
        typedef std::vector<Minerva::NeutrinoInt*> NeutrinoVect;
        
    public:
        
        //! Standard constructor
        CCProtonPi0(const std::string& type, 
                    const std::string& name, 
                    const IInterface* parent );
        
        //! Destructor (mandatory for inheritance)
        virtual ~CCProtonPi0(){};
        
        StatusCode initialize();
        StatusCode finalize();
        
        //! Reconstruct the event (mandatory for inheritance)
        StatusCode reconstructEvent(Minerva::PhysicsEvent* event, 
                                    Minerva::GenMinInteraction* truthEvent = NULL ) const;
        
        //! Attach an interpretations to the event (mandatory for inheritance)
        StatusCode interpretEvent(  const Minerva::PhysicsEvent* event, 
                                    const Minerva::GenMinInteraction* truthEvent, 
                                    NeutrinoVect& interaction_hyp ) const;
        
        StatusCode tagTruth( Minerva::GenMinInteraction* truthEvent ) const;
        
        bool truthIsPlausible( const Minerva::PhysicsEvent*) const;
        
    private:
        std::vector<std::string>   m_dedx_uncertainties;
        
        // Fiducial Volume
        double m_fidHexApothem;
        double m_fidUpStreamZ;
        double m_fidDownStreamZ;
        
        // Analyzable Volume
        double m_recoHexApothem;
        double m_recoUpStreamZ;
        double m_recoDownStreamZ;
        
        double m_beamAngleBias;
        
        double m_detectableGammaE;
        double m_detectablePi0KE;
        double m_detectableProtonKE;
        
        double m_minMuonScore;
        double m_minProtonScore;
        
        // Algorihm Flow
        bool m_writeFSParticle_Table;
        bool m_store_all_events;
        bool m_makeShortTracks;
        bool m_DoPlausibilityCuts;
        bool m_DoTruthMatch;
        bool m_reconstruct_NoProtonEvents;
        bool m_applyExtraMichelCuts;
        
        // Tool Names
        std::string m_particleToolName;
        std::string m_particleToolAlias;
        std::string m_protonUtilsAlias;
        std::string m_michelTrkToolAlias;
        std::string m_michelVtxToolAlias;
        
        // Prong and Cluster Colors
        int m_Color_muonProng;
        int m_Color_protonProng;
        int m_Color_primaryVertex;
        int m_Color_secondaryVertex;
        int m_Color_endPointVertex;
        int m_Color_unattachedProng;
        int m_Color_Gamma1Prong;
        int m_Color_Gamma2Prong;
        int m_Color_GammaOtherProng;
        int m_Color_clusterUnused;
        int m_Color_clusterUsed;
        int m_Color_VertexFila;
        int m_Color_VertexSphere;
        int m_Color_RejectedBlob;
             
        // VtxBlob
        bool 	 m_sphereVertex;
        double  m_maxSearchD;
        double  m_maxStartingD;
        double  m_maxSearchGap;		
        bool	 m_filamentVertex;
        double  m_maxSearchDFila;
        double  m_maxStartingDFila;
        double  m_maxSearchGapFila;
        bool    m_filterClusterTypes;
        bool    fSkipLowEnergyClusterVtxEnergy;
        bool    fThresholdVertexEnergy;
        
        // ConeBlobs
        double m_qOverpChargeCut;      ///< q/p charge cut
        double m_energyHoughlimit;     ///< Energy limit to start using Hough T.
        double m_rejectedClustersTime; ///< window time  to allow clusters
        
        bool m_TrytoRecoverBlobReco;
        bool m_ApplyAttenuationCorrection;
        
        bool   m_AllowUVMatchWithMoreTolerance;
        double m_UVMatchTolerance;
        double m_UVMatchMoreTolerance;
  
        double m_extraEnergyCylinderRadius;           ///< Cylinder Cut (mm) 
        double m_extraEnergyCylinderUpstreamLength;   ///< Cylinder Cut (mm) 
        double m_extraEnergyCylinderDownstreamLength; ///< Cylinder Cut (mm) 
        double m_extraEnergyLowerTimeWindow;          ///< Cylinder Cut (ns) 
        double m_extraEnergyUpperTimeWindow;          ///< Cylinder Cut (ns) 
        double m_extraEnergyPECut;                    ///< Cylinder Cut (MeV)?
        
        double m_maxSeedLongTrackChi2;
        double m_maxSeedShortTrackChi2;
        
        TRandom3*                 m_randomGen;
        unsigned long int         m_randomSeed;
        
        Minerva::DeDetector*        m_InnerDetector;
        Minerva::DeOuterDetector*   m_OuterDetector;
        
        // Duplicated variables -- will fix in future
        const Minerva::DeDetector*      m_idDet;                ///< Inner detector
        const Minerva::DeOuterDetector* m_odDet;                ///< Outer detector
        
        IAnchoredTrackFormation* m_anchoredTracker;
        IBlobCreatorUtils* m_blobUtils;
        ICCPionIncUtils* m_ccPionIncUtils;
        ICalorimetryUtils* m_caloUtils;
        IConeUtilsTool* m_coneUtilsTool;
        IEnergyCorrectionTool* m_energyCorrectionTool;
        IExtraEnergyTool* m_extraEnergyTool;
        IGeomUtilSvc* m_GeomUtilSvc;
        IGetCalAttenuation* m_AttenuationCorrectionTool;
        IGetDeadTime* m_getDeadTimeTool;
        IGiGaGeomCnvSvc* m_gigaCnvSvc;
        IHitTaggerTool* m_hitTagger;
        IHoughBlob* m_idHoughBlob;
        IHoughTool* m_idHoughTool;
        IIDAnchoredBlobCreator* m_idConeScanBlob;
        IIDAnchoredBlobCreator* m_stopPointBlobTool;
        IMCTrackTool* m_MCTrackTool;
        IMichelTool* m_michelTrkTool;
        IMichelTool* m_michelVtxTool;
        IMinervaCoordSysTool* m_coordSysTool;
        IMinervaMathTool* m_mathTool;
        IMinervaObjectAssociator* m_objectAssociator;
        IODProngClassificationTool* m_odMatchTool;
        IParticleMakerTool* m_particleMaker;
        IParticleTool* m_LikelihoodPIDTool;
        IParticleTool* m_particleTool;
        IProngClassificationTool* m_prongIntersection;
        IProtonUtils* m_protonUtils;
        IRecoObjectTimeTool* m_recoTimeTool;
        IVertexEnergyStudyTool* m_vertexEnergyStudyTool;
        IVertexFitter* m_vertexFitter;

        //! Private Functions
        SmartRefVector<Minerva::IDCluster> FilterInSphereClusters(  const SmartRefVector<Minerva::IDCluster>& clusters, const double sphereRadius, std::vector<double>& radii) const;
        StatusCode HoughBlob( SmartRefVector<Minerva::IDCluster> idClusters, std::vector<Minerva::IDBlob*>& outBlobs) const;
        StatusCode ODActivity( Minerva::PhysicsEvent *event, std::vector<Minerva::IDBlob*> idBlobs ) const;
        StatusCode getNearestPlane( double z, int & module_return, int & plane_return) const;
        StatusCode interpretFailEvent( Minerva::PhysicsEvent* event ) const;
        bool ConeBlobs( Minerva::PhysicsEvent *event ) const;
        bool PreFilterPi0( Minerva::PhysicsEvent *event ) const;
        bool AreBlobsGood() const;
        void VtxBlob( Minerva::PhysicsEvent *event ) const;
        bool checkMichel(Minerva::GenMinInteraction* truthEvent) const;
        bool checkPionAbsorption(Minerva::GenMinInteraction* truthEvent) const;
        void correctProtonProngEnergy(  SmartRef<Minerva::Prong>& protonProng, double& p_calCorrection, double& p_visEnergyCorrection ) const;
        bool createTrackedParticles(Minerva::ProngVect& prongs ) const;
        bool getProtonProng(    Minerva::ProngVect& primaryProngs ) const;
        bool isBackgroundQELike(Minerva::GenMinInteraction* truthEvent) const;
        bool isMichelProngGood(Minerva::Prong &michelProng) const;
        bool isMotherPrimary(std::vector<int>& motherList, int mother ) const;
        bool isSinglePi0( Minerva::GenMinInteraction* truthEvent, int nPi0, int nGamma) const;
        bool setMuonData( Minerva::NeutrinoInt* nuInt ) const;
        bool setPi0Data( Minerva::NeutrinoInt* nuInt ) const;
        bool setProtonData( Minerva::NeutrinoInt* nuInt ) const;
        bool tagSignal( Minerva::GenMinInteraction* truthEvent) const;
        double Calc_QSq(double Enu) const;
        double Calc_WSq(double Enu, double QSq) const;
        double Calc_Enu_Cal(double hadronEnergy) const;
        double Calc_Enu_1Track() const;
        double Calc_Enu_2Track() const;
        double Calc_Longitudinal_Momentum(Gaudi::LorentzVector particle_4P) const;
        double Calc_Px_wrt_Beam(Gaudi::LorentzVector particle_4P) const;
        double Calc_Py_wrt_Beam(Gaudi::LorentzVector particle_4P) const;
        double CalcMinBlobSeparation( const Minerva::IDBlob* blob) const;
        double calcDistance( double x1, double y1, double z1,double x2, double y2, double z2) const;
        double getClusterEnergy( Minerva::PhysicsEvent* event, std::string input_clusterType) const;
        int getMichelPion(std::vector<int>& piList, int ID ) const;
        std::pair<int,double> OneParLineFitBlob(const Minerva::IDBlob* blob) const;
        std::vector<double> GetBlobClusterEnergy( const Minerva::IDBlob* blob ) const;
        void ApplyAttenuationCorrection(Minerva::IDBlob* blob) const;
        void Calculate_dEdx( const Minerva::IDBlob* blob, Minerva::PhysicsEvent* event, unsigned int blob_number) const;
        void SetSignalKinematics(Minerva::GenMinInteraction* truthEvent) const;
        void processBlobs( Minerva::PhysicsEvent *event, std::vector<Minerva::IDBlob*> idBlobs) const;
        void saveMichelElectron(Minerva::GenMinInteraction* truthEvent, int muon_ID) const;
        void saveMichelProngToNTuple(Minerva::PhysicsEvent* event, Minerva::Prong &michelProng) const;
        void setBlobData(Minerva::PhysicsEvent* event) const;
        void setEventKinematics(Minerva::NeutrinoInt* nuInt, double hadronVisibleEnergy) const;
        void setPi0GenieRecord(Minerva::GenMinInteraction* truthEvent) const;
        void setTargetMaterial(Minerva::GenMinInteraction* truthEvent) const;
        void setTrackProngTruth( Minerva::NeutrinoInt* neutrino, Minerva::ProngVect& prongs ) const;
        void setVertexData( Minerva::NeutrinoInt* nuInt, const Minerva::PhysicsEvent* event ) const;
        void tagBackground(Minerva::GenMinInteraction* truthEvent) const;
        void writeBackgroundType(Minerva::GenMinInteraction* truthEvent) const;
        void writeEventRecord(Minerva::GenMinInteraction* truthEvent, bool isSignal) const;
        void writeFSParticleTable(bool isSignal) const;
        void setTrackDirection( Minerva::Track* track, Minerva::Vertex* vertex ) const;
        Minerva::IDClusterVect getClusters( Minerva::PhysicsEvent* event ) const;
        bool createdAnchoredShortTracks( Minerva::PhysicsEvent* event, Minerva::Vertex* vertex, bool make_primary_short_tracks ) const;
};

#endif // CCPROTONPI0_H 

