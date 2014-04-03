#undef NDEBUG
#include <cassert>
#include <algorithm>
#include <iostream>
#include <iomanip>

#include <TString.h>
#include <TDatabasePDG.h>
#include <TVector3.h>
#include <TMath.h>

#include "CCDeltaPlusAna.h"
#include "CCDeltaPlus/IHoughBlob.h"
#include "CCDeltaPlus/IHoughTool.h"

#include "Event/DAQHeader.h"
#include "Event/GenMinHeader.h"
#include "Event/MCIDDigit.h"
#include "Event/TG4Trajectory.h"
#include "Event/IDCluster.h"
#include "Event/VectorTypeDefs.h"
#include "Event/Vertex.h"
#include "Event/Track.h"
#include "Event/TimeSlice.h"

#include "GeoUtils/IMinervaCoordSysTool.h"
#include "GeoUtils/INuclearTargetTool.h"
#include "GeoUtils/NuclearTarget.h"
#include "GeoUtils/IAbsorberStacker.h"
#include "GeoUtils/AbsorberStackerIDContainer.h"

#include "MinervaUtils/IMinervaMathTool.h"
#include "MinervaUtils/IHitTaggerTool.h"
#include "MinervaUtils/MinervaObjectSort.h"

#include "MinervaDet/IGeomUtilSvc.h"
#include "MinervaDet/DeSubdet.h"
#include "MinervaDet/DeDetector.h"
#include "MinervaDet/DePlane.h"
#include "ODDet/DeOuterDetector.h"

#include "EnergyRecTools/IExtraEnergyTool.h"
#include "EnergyRecTools/IEnergyCorrectionTool.h"
#include "EnergyRecTools/IEnergyLoss.h"

#include "RecUtils/IClusterUtilsTool.h"
#include "RecUtils/IConeUtilsTool.h"
#include "RecUtils/Cone.h"
#include "RecUtils/ParticleExtraDataDefs.h"

#include "RecInterfaces/IIDBlobCreator.h"
#include "RecInterfaces/IODBlobCreator.h"
#include "RecInterfaces/ITrackLinearPropagator.h"
#include "RecInterfaces/IFiducialPointTool.h"
#include "RecInterfaces/IAnchoredTrackFormation.h"
#include "RecInterfaces/ITrackPropagator.h"
#include "RecoStudies/IVertexEnergyStudyTool.h"

#include "BlobFormation/IIDIsolatedIDBlobCreator.h"
#include "BlobFormation/IIDBlobSeedingTool.h"
#include "BlobFormation/IIDAnchoredBlobCreator.h"
#include "BlobFormation/IBlobCreatorUtils.h" 

#include "AnaUtils/IMuonUtils.h"
#include "AnaUtils/IPhysicsCalculator.h"
#include "AnaUtils/IProtonUtils.h"
#include "AnaUtils/IMCTrackTool.h"
#include "AnaUtils/MCTrack.h"

#include "ProngMaker/IProngClassificationTool.h"
#include "ProngMaker/IODProngClassificationTool.h"
#include "ProngMaker/IMichelTool.h"

#include "ParticleMaker/IParticleTool.h"
#include "ParticleMaker/IParticleMakerTool.h"

#include "DetDesc/Material.h"
#include "MinosInterface/IMinosMinervaTransformer.h"
#include "EventRecInterfaces/IPrimaryBlobProngTool.h"
#include "VertexCreation/IVertexFitter.h"

#include "TraverseHistory.h"
#include "TrackTruthInfo.h"
#include "DigitVectorTruthInfo.h"
#include "HitVectorTruthInfo.h"
#include "OneParLineFit.h"
#include "Par.h"
#include "AngleScan.h"
#include "ClusterVectorInfo.h"


/* 
    Global variables, instead of class data members so that they can be
    assigned in 'const' methods. In general, we can use 'mutable' members,
    but it does not seem to work in the framework 
*/
SmartRef<Minerva::Track>    m_MuonTrack;
SmartRef<Minerva::Prong>    m_MuonProng;
SmartRef<Minerva::Particle> m_MuonParticle;
SmartRef<Minerva::Vertex>   m_PrimaryVertex;

const double double_epsilon = std::numeric_limits<double>::epsilon();

namespace {
    bool InsideHexagon(double x, double y, double w)
    {
        double max = w/2;
        double min = -w/2;
        
        double deltaphi = 60.0*TMath::DegToRad();
        double phi = 0.0;
        for (int i = 0; i < 3; ++i) {
            double xtmp = std::cos(phi)*x + std::sin(phi)*y;
            phi += deltaphi;
            
            if (xtmp < min || xtmp > max) return false;
        }
        
        return true;
        
    }

    const double rad2deg = TMath::RadToDeg();
    const double mp = 938.272046;
    const double mn = 939.565379;
    const double mmuon  = 105.65;
    const double theta_b = 3.3*TMath::DegToRad();
    const double sintheta_b = std::sin(theta_b);
    const double costheta_b = std::cos(theta_b);

    class greaterTotalEnergy : public std::binary_function
    <
        const Minerva::TG4Trajectory*,
        const Minerva::TG4Trajectory*,
        bool
        > {
    public:
        bool operator() (const Minerva::TG4Trajectory* lhs,
                        const Minerva::TG4Trajectory* rhs) {
            return lhs->GetInitialMomentum().E() > rhs->GetInitialMomentum().E();
        }
    
    };


    class PdgEqual : public std::unary_function
    <
        const Minerva::TG4Trajectory*,
        bool
        > {
    public:
        PdgEqual(int pdg) : fPdg(pdg) {}
        bool operator()(const Minerva::TG4Trajectory* lhs) {
            return lhs->GetPDGCode() == fPdg;
        }
        
    private:
        int fPdg;
    };
    
}
    
using std::setw;


DECLARE_TOOL_FACTORY( CCDeltaPlusAna );
    
    
//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
CCDeltaPlusAna::CCDeltaPlusAna(const std::string& type, const std::string& name, const IInterface* parent ) :
MinervaAnalysisTool( type, name, parent ) 
{
    info() <<"Enter CCDeltaPlusAna::CCDeltaPlusAna() -- Default Constructor" << endmsg;
    info() <<"=============================================================================" <<endmsg;
    declareInterface<IInteractionHypothesis>(this);
    
    //! mandatory declaration of analysis signature: CCDeltaPlusReco
    m_anaSignature = "CCDeltaPlusAna";
    
    // Protected properties from IInteractionHypothesis.
    m_hypMeths.push_back( m_anaSignature );
    declareProperty("HypothesisMethods", m_hypMeths);
    
    debug() << " CCDeltaPlusAna Hypothesis added " << endmsg;
    
    // Vertex Cut : Is this event signal?
    declareProperty( "FiducialApothem",      m_fiducialApothem      =  850.0* CLHEP::millimeter ); 
    declareProperty( "FiducialUpstreamZ",    m_fiducialUpstreamZ    = 5912.0* CLHEP::millimeter ); 
    declareProperty( "FiducialDownstreamZ",  m_fiducialDownstreamZ  = 8445.0* CLHEP::millimeter ); 
    
    // Vertex blob
    declareProperty( "UseSphereVertex", 		m_sphereVertex	= true );
    declareProperty( "SphereMaxSearchDistance", 	m_maxSearchD  	=  90.0 * CLHEP::millimeter ); 
    declareProperty( "SphereMaxStartingDistance", m_maxStartingD	=  90.0 * CLHEP::millimeter ); 
    declareProperty( "SphereMaxAllowedSearchGap", m_maxSearchGap	= 180.0 * CLHEP::millimeter ); 
    declareProperty( "UseFilamentVertex", 	m_filamentVertex = false );
    declareProperty( "FilamentMaxSearchDistance",   m_maxSearchDFila    = 500.0 * CLHEP::millimeter );
    declareProperty( "FilamentMaxStartingDistance", m_maxStartingDFila  = 68.0 * CLHEP::millimeter ); 
    declareProperty( "FilamentMaxAllowedSearchGap", m_maxSearchGapFila  = 51.0 * CLHEP::millimeter ); 
    declareProperty( "FilterClusterTypes",   m_filterClusterTypes    = true );
    
    //dEdx
    declareProperty( "NumberPlanesdEdX", 	m_planesdEdx 	= 4 ); 
    m_distance_dEdx = m_planesdEdx*17 * CLHEP::millimeter; 
    
    // Ana tools
    declareProperty( "qOverpChargeCut",       m_qOverpChargeCut       = 0   );
    declareProperty( "HoughEnergyLimit",      m_energyHoughlimit      = 900 * CLHEP::MeV ); 
    declareProperty( "RejectedClustersTime",  m_rejectedClustersTime  = 25 * CLHEP::ns );
    
    /* Number from Cesar Sotelo's studies.
        For forward-going track, the cylinder starts 20cm DOWNSTREAM of the first track node
    */
    declareProperty( "ExtraEnergyCylinderUpstreamOffset",  m_extraEnergyCylinderUpstreamLength = -200.0*CLHEP::mm);
    declareProperty( "ExtraEnergyCylinderDownstreamOffset",m_extraEnergyCylinderDownstreamLength = 50.0*CLHEP::mm);
    declareProperty( "ExtraEnergyCylinderRadius",          m_extraEnergyCylinderRadius = 50*CLHEP::mm);
    declareProperty( "ExtraEnergyLowerTimeWindow",         m_extraEnergyLowerTimeWindow = 25.0*CLHEP::ns);
    declareProperty( "ExtraEnergyUpperTimeWindow",         m_extraEnergyUpperTimeWindow = 25.0*CLHEP::ns);
    declareProperty( "ExtraEnergyPhotoElectronCut",        m_extraEnergyPECut = 0.0);
    
    declareProperty( "new_impl", new_impl_ = true);
    declareProperty( "ProcessOneEvent", fProcessOneEvent = false);
    declareProperty( "Run", fRunNumber = -1);
    declareProperty( "Subrun", fSubrunNumber = -1);
    declareProperty( "Gate", fGateNumber = -1);
    
    declareProperty( "UVMatchTolerance", fUVMatchTolerance = 10.0*CLHEP::mm);
    declareProperty( "UVMatchMoreTolerance", fUVMatchMoreTolerance = 100.0*CLHEP::mm);
    declareProperty( "AllowUVMatchWithMoreTolerance", fAllowUVMatchWithMoreTolerance = true);
    
    
   declareProperty("ApplyCalCorrection",          m_applyCalCorrection      = true);
   declareProperty("MinimumMuonMVAPIDScore",      m_min_mva_score           = 0.0);
   declareProperty("CCQEBindingEnergyMeV",        m_CCQEBindingEnergyMeV    = 34.0*CLHEP::MeV);

   declareProperty( "MaxSearchDistance_VESTool",      m_maxSearchDistance_VESTool     = 300.0*CLHEP::mm );
   declareProperty( "MaxStartingDistance_VESTool",    m_maxStartingDistance_VESTool   = 300.0*CLHEP::mm );
   declareProperty( "MaxAllowedSearchGap_VESTool",    m_maxAllowedSearchGap_VESTool   = 1e6*CLHEP::mm );

   declareProperty("ClusterVtxMinTimeWindow",     m_clusMinTimeWindow       = -20*CLHEP::ns);
   declareProperty("ClusterVtxMaxTimeWindow",     m_clusMaxTimeWindow       = 35*CLHEP::ns);
   declareProperty("FuzzRadius",                  m_fuzzRadius              = 50.0*CLHEP::mm);

   declareProperty("VtxProngNumSearchRadii",      m_numSearchRadii          = 1);
   declareProperty("VtxProngSearchStepSize",      m_searchStepSize          = 100*CLHEP::mm);
   declareProperty("VtxProngMaxSearchDistance",   m_maxSearchDistance       = 100*CLHEP::mm);
   declareProperty("VtxProngMaxStartDistance",    m_maxStartingDistance     = 100*CLHEP::mm);
   declareProperty("VtxProngMaxSearchGap",        m_maxAllowedSearchGap     = 1e6*CLHEP::mm);
   declareProperty("VtxProngMaxSepBlobVertex",    m_maxSeparationBlobVertex = 150*CLHEP::mm);

   declareProperty("NuclearTargetToolAlias",      m_nuclearTargetToolAlias  = "CCDeltaPlusAnaTargetTool");
   declareProperty("ProtonUtilsAlias",            m_protonUtilsAlias        = "CCDeltaPlusAnaProtonUtils");

   declareProperty("PrimaryBlobProngAlias",       m_primaryBlobProngAlias   = "CCDeltaPlusAnaPrimaryBlobProng");

   declareProperty("AnchorShortTrackerAlias",     m_anchorShortTrackerAlias = "CCDeltaPlusAnaShortTracker");
   declareProperty("VertexFitterAlias",           m_vertexFitterAlias       = "CCDeltaPlusAnaVertexFitter");

   declareProperty("ODMatchProngToolAlias",       m_odMatchAlias            = "CCDeltaPlusAnaODMatchTool");
   declareProperty("ProngIntersectionAlias",      m_prongIntersectionAlias  = "CCDeltaPlusAnaProngIntersect");
   declareProperty("ParticleMakerAlias",          m_particleMakerAlias      = "CCDeltaPlusAnaParticleMaker");
   declareProperty("MichelToolAlias",             m_michelToolAlias         = "CCDeltaPlusAnaMichelTool"); 

   declareProperty("BlobCreatorUtilsAlias",       m_blobCreatorUtilsAlias   = "CCDeltaPlusAnaBlobCreator"); 
   declareProperty("VertexEnergyStudyToolAlias",  m_vtxEngStudyToolAlias    = "CCDeltaPlusAnaVESTool"); 

   declareProperty("FidVolApothem",               m_fvApothem = 860.0*CLHEP::mm);
   declareProperty("FidVolUpstreamZ",             m_fvUpstreamZ);
   declareProperty("FidVolDownstreamZ",           m_fvDownstreamZ);
   declareProperty("TargetNames",                 m_targetNames);

   declareProperty("Target1VertexZUSCut",         m_tar1VertexZUSCut     = 4467.19*CLHEP::mm);
   declareProperty("Target1VertexZDSCut",         m_tar1VertexZDSCut     = 4510.09*CLHEP::mm);

   declareProperty("Target2VertexZUSCut",         m_tar2VertexZUSCut     = 4688.26*CLHEP::mm);
   declareProperty("Target2VertexZDSCut",         m_tar2VertexZDSCut     = 4731.16*CLHEP::mm);

   declareProperty("Target3VertexZUSCut",         m_tar3VertexZUSCut     = 4892.41*CLHEP::mm);
   declareProperty("Target3VertexZDSCut",         m_tar3VertexZDSCut     = 4989.68*CLHEP::mm);

   declareProperty("Target4VertexZUSCut",         m_tar4VertexZUSCut     = 5630.80*CLHEP::mm);
   declareProperty("Target4VertexZDSCut",         m_tar4VertexZDSCut     = 5664.80*CLHEP::mm);

   declareProperty("Target5VertexZUSCut",         m_tar5VertexZUSCut     = 5756.71*CLHEP::mm);
   declareProperty("Target5VertexZDSCut",         m_tar5VertexZDSCut     = 5801.24*CLHEP::mm);

   declareProperty("MuonProngColor",              m_muonProngColor       = 0x228B22); //-- green
   declareProperty("ProtonProngColor",            m_protonProngColor     = 0x9932CC); //-- purple
   declareProperty("PrimaryVertexProngColor",     m_primaryVertexColor   = 0xFF0000); //-- red
   declareProperty("SecondaryVertexProngColor",   m_secondaryVertexColor = 0xFFA500); //-- orange
   declareProperty("TrackEndProngColor",          m_endPointVertexColor  = 0xFF1493); //-- pink
   declareProperty("UnattachedProngColor",        m_unattachedProngColor = 0x0000FF); //-- blue

   //! Loose fiducial volume requirement for truth michel
   declareProperty("MichelDownStreamZ",          m_michel_downstreamZ     = 8524.19*CLHEP::mm); //module 83?
   declareProperty("MichelUpStreamZ",            m_michel_upstreamZ       = 5900.91*CLHEP::mm); //module 25?

   //! fiducial volume for the nuclear targets region
   m_fvUpstreamZ.push_back(4284.46*CLHEP::mm);
   m_fvDownstreamZ.push_back(5975.04*CLHEP::mm);
   m_targetNames.push_back("NuclearTargets");

   //! values corresponds the zposition at the face of the upstream module
   m_fvUpstreamZ.push_back(4417.10*CLHEP::mm); 
   m_fvUpstreamZ.push_back(4638.18*CLHEP::mm); 
   m_fvUpstreamZ.push_back(4859.33*CLHEP::mm);
   m_fvUpstreamZ.push_back(5580.80*CLHEP::mm);
   m_fvUpstreamZ.push_back(5713.53*CLHEP::mm);

   //! values corresponds to the zposition at the back of the downstream module
   m_fvDownstreamZ.push_back(4543.01*CLHEP::mm); 
   m_fvDownstreamZ.push_back(4764.09*CLHEP::mm);
   m_fvDownstreamZ.push_back(5029.38*CLHEP::mm);
   m_fvDownstreamZ.push_back(5706.71*CLHEP::mm);
   m_fvDownstreamZ.push_back(5839.35*CLHEP::mm);
   
   //! target names
   m_targetNames.push_back("Target001");
   m_targetNames.push_back("Target002");
   m_targetNames.push_back("Target003");
   m_targetNames.push_back("Target004");
   m_targetNames.push_back("Target005");
    
    info() <<"Exit CCDeltaPlusAna::CCDeltaPlusAna() -- Default Constructor" << endmsg;
    
}
    
//=============================================================================
// Initialize
//=============================================================================
StatusCode CCDeltaPlusAna::initialize()
{
    info() <<"Enter CCDeltaPlusAna::initialize()" << endmsg;
    info() <<"=============================================================================" <<endmsg;
    
    //! Initialize the base class.  This will fail if you did not define m_anaSignature.
    StatusCode sc = this->MinervaAnalysisTool::initialize();
    if( sc.isFailure() ) { 
        return Error( "Failed to initialize!", sc ); 
    
    }
    info()<<"   initialized MinervaAnalysisTool"<<endmsg;
    
    /*
    ============================================================================
        Initializing Analysis Tools
    ============================================================================
    */

    try{ 
        m_nuclearTargetTool = tool<INuclearTargetTool>("NuclearTargetTool", m_nuclearTargetToolAlias); 
            m_nuclearTargetTool->m_locked = false;
            m_nuclearTargetTool->addAllPassiveNuclearTargets();
        
            Minerva::NuclearTarget* refTarget     = m_nuclearTargetTool->addNuclearTarget(5);
            Minerva::NuclearTarget* plasticTarget = m_nuclearTargetTool->addNuclearTarget(27,60,refTarget->getTarget());
    
            m_targetNames.push_back( plasticTarget->getName() ); 
            m_fvUpstreamZ.push_back( plasticTarget->getZStart() );
            m_fvDownstreamZ.push_back( plasticTarget->getZEnd() );
    
            m_nuclearTargetTool->lock();
            m_targets = m_nuclearTargetTool->getNuclearTargets();
    } 
    catch( GaudiException& e ) {
        error() << "Could not obtain NuclearTargetTool: " << m_nuclearTargetToolAlias << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{ 
        m_protonUtils = tool<IProtonUtils>("ProtonUtils", m_protonUtilsAlias); 
    }
    catch( GaudiException& e ){
        error() << "Could not obtain ProtonUtils: " << m_protonUtilsAlias << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{ 
        m_coneUtilsTool = tool<IConeUtilsTool>("ConeUtilsTool"); 
    }
    catch( GaudiException& e ){
        error() << "Could not obtain ConeUtilsTool!" << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{ 
        m_energyCorrectionTool = tool<IEnergyCorrectionTool>("EnergyCorrectionTool"); 
    }
    catch( GaudiException& e ){
        error() << "Could not obtain EnergyCorrectionTool!" << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{ 
        m_energyLoss = tool<IEnergyLoss>("EnergyLoss"); 
    } 
    catch( GaudiException& e ) {
        error() << "Could not obtain EnergyLoss!" << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{ 
        m_absorberStacker = tool<IAbsorberStacker>("AbsorberStacker"); 
    }
    catch( GaudiException& e ) {
        error() << "Could not obtain AbsorberStacker!" << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{ 
        m_coordSysTool = tool<IMinervaCoordSysTool>("MinervaCoordSysTool"); 
    }
    catch( GaudiException& e ) {
        error() << "Could not obtain MinervaCoordSysTool!" << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{ 
        m_propagateToMinos = tool<ITrackPropagator>("TrackPropagator"); 
    }
    catch( GaudiException& e ) {
        error() << "Could not obtain TrackPropagator!" << endmsg;
        return StatusCode::FAILURE;
    } 
    
    try{ 
        m_mmt = tool<IMinosMinervaTransformer>("MinosMinervaTransformer"); 
    }
    catch( GaudiException& e ) {
        error() << "Could not obtain MinosMinervaTransformer!" << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{ 
        m_primaryBlobProngTool = tool<IPrimaryBlobProngTool>("PrimaryBlobProngTool",m_primaryBlobProngAlias); 
    }
    catch( GaudiException& e ){
        error() << "Could not obtain PrimaryBlobProngTool: " << m_primaryBlobProngAlias << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{ 
        m_anchoredTracker = tool<IAnchoredTrackFormation>("AnchoredShortTracker",m_anchorShortTrackerAlias); 
    }
    catch( GaudiException& e ){
        error() << "Could not obtain AnchoredShortTracker: " << m_anchorShortTrackerAlias << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{ 
        m_vertexFitter = tool<IVertexFitter>("VertexFitterKalman",m_vertexFitterAlias); 
    }
    catch( GaudiException& e ){
        error() << "Could not obtain VertexFitterKalman: " << m_vertexFitterAlias << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{ 
        m_prongIntersection = tool<IProngClassificationTool>("ProngIntersectionTool", m_prongIntersectionAlias); 
    }
    catch( GaudiException& e ){
        error() << "Could not obtain ProngClassificationTool: " << m_prongIntersectionAlias << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{ 
        m_odMatchTool = tool<IODProngClassificationTool>("ODTrackMatchTool", m_odMatchAlias); 
    }
    catch( GaudiException& e ){
        error() << "Could not obtain ODTrackMatchTool: " << m_odMatchAlias << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{ 
        m_particleMaker = tool<IParticleMakerTool>("ParticleMakerTool", m_particleMakerAlias); 
    }
    catch( GaudiException& e ){
        error() << "Could not obtain ParticleMakerTool: " << m_particleMakerAlias << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{ 
        m_mcTrackTool = tool<IMCTrackTool>("MCTrackTool"); 
    }
    catch( GaudiException& e){
        error() << "Could not obtain MCTrackTool! " <<endmsg;
        return StatusCode::FAILURE;
    }
   
    try{ 
        m_hitTagger = tool<IHitTaggerTool>( "HitTaggerTool" ); 
    }
    catch( GaudiException& e ) {
        error() << "Could not obtain tool: HitTaggerTool!" << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{ 
        m_michelTool = tool<IMichelTool>("MichelTool",m_michelToolAlias); 
    }
    catch( GaudiException& e ) {
        error() << "Could not obtain tool: MichelTool!" << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{ 
        m_blobCreatorUtils = tool<IBlobCreatorUtils>("BlobCreatorUtils", m_blobCreatorUtilsAlias); 
    }
    catch( GaudiException& e ){
        error() << "Could not obtain BlobCreatorUtils: " << m_blobCreatorUtilsAlias << endmsg;
        return StatusCode::FAILURE;
    } 
    
    try{ 
        m_vertexEnergyStudyTool = tool<IVertexEnergyStudyTool>("VertexEnergyStudyTool", m_vtxEngStudyToolAlias); 
    }
    catch( GaudiException& e ){
        error() << "Could not obtain VertexEnergyStudyTool: " << m_vtxEngStudyToolAlias << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{
        m_extraEnergyTool = tool<IExtraEnergyTool>("ExtraEnergyTool");
    }
    catch(GaudiException& e){
        error()<<"Could not obtain tool: ExtraEnergyTool<oaltinok_version>" << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{
        m_idConeScanBlob = tool<IIDAnchoredBlobCreator>("ConeScanIDBlobCreator");
    } 
    catch (GaudiException& e){
        error() << " Could not obtain tool: ConeScanIDBlobCreator " << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{
        m_trackPropagator = tool<ITrackLinearPropagator>("TrackLinearPropagator");
    } 
    catch(GaudiException& e){
        error()<<"Could not obtain tool: TrackLinearPropagator<oaltinok_version>" << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{
        m_idHoughBlob = tool<IHoughBlob>("HTBlob");
    }
    catch (GaudiException& e){
        error() << " Could not obtain tool: HBlob " << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{
        m_idHoughTool = tool<IHoughTool>("HTtool");
    }
    catch (GaudiException& e){
        error() << " Could not obtain tool: HTtool " << endmsg;
        return StatusCode::FAILURE;
    }
    
    try {
        m_blobSeedingTool = tool<IIDBlobSeedingTool>("BlobSeedingTool"); 
    }
    catch( GaudiException& e ) { 
        error() << "Could not obtain tool: BlobSeedingTool !" << endmsg;
        return StatusCode::FAILURE;
    }
    
    try {
        m_mathTool = tool<IMinervaMathTool>("MinervaMathTool");
    }
    catch( GaudiException& e ) {
        error() << "Could not obtain tool: MinervaMathTool!" << endmsg;
        return StatusCode::FAILURE;
    }
    
    /*
    ============================================================================
        Initializing Detector
    ============================================================================
    */
    debug()<<"Initializing Inner and Outer Detector"<<endmsg;
    try{ 
        m_idDet = getDet<Minerva::DeDetector>("/dd/Structure/Minerva/Detector/InnerDetector"); 
    }
    catch( GaudiException& e ){
        error() << "Could not retrieve the Inner Detector!" << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{ 
        m_odDet = getDet<Minerva::DeOuterDetector>("/dd/Structure/Minerva/Detector/OuterDetector"); 
    }
    catch( GaudiException& e ) {
        error() << "Could not retrieve the Outer Detector!" << endmsg;
        return StatusCode::FAILURE;
    }
    
    
//     service("GeomUtilSvc", m_GeomUtilSvc, true);
//     m_idDet = m_GeomUtilSvc->getIDDet();
//     m_odDet = m_GeomUtilSvc->getODDet();
    
    Constants().k_hcal = 0.088;
    Constants().k_ecal = 0.360;
    Constants().k_trkr = 0.760;
    Constants().coangle = 15.0;
    Constants().clength = 3000.0;
    Constants().xminevis = 20.0;
    Constants().uminevis = 10.0;
    Constants().vminevis = 10.0;
    
    //! declare common branches
    declareCommonPhysicsAnaBranches();
    declareMinosMuonBranches();
    declareGenieWeightBranches();
    
    /*
    ============================================================================
        Select the branches you want in your AnaTuple
    ============================================================================
    */
    
    /*
    ============================================================================
        NukeCCQE Variables - Tammy's Package
    ============================================================================
    */
    
    debug() << "Initializing Tammy's Variables Started" << endmsg;

    declareIntTruthBranch("tammy_has_michel_electron", 0);
    declareContainerDoubleTruthBranch("tammy_has_michel_from_tammy_pion_plus_momentum");
    declareContainerDoubleTruthBranch("tammy_has_michel_from_tammy_pion_minus_momentum");
    
    //-- Event Michel variables
    declareContainerIntEventBranch("tammy_has_michel_category");
    declareContainerIntEventBranch("tammy_has_michel_vertex_type");
    declareContainerIntEventBranch("tammy_has_michel_in_vertex_point");
    declareContainerDoubleEventBranch("tammy_has_michel_distance");
    declareContainerDoubleEventBranch("tammy_has_michel_energy");
    declareContainerIntEventBranch("tammy_has_michel_ndigits");
    declareContainerDoubleEventBranch("tammy_has_michel_time_diff");
    
    declareIntEventBranch("tammy_genie_n_neutrinos",     0);
    declareIntEventBranch("tammy_genie_n_muons",         0);
    declareIntEventBranch("tammy_genie_n_mesons",        0);
    declareIntEventBranch("tammy_genie_n_heavy_baryons", 0);
    declareIntEventBranch("tammy_genie_n_photons",       0);
    declareIntEventBranch("tammy_genie_n_protons",       0);
    declareIntEventBranch("tammy_genie_n_neutrons",      0);
    declareIntEventBranch("tammy_genie_n_pions",         0);
    declareIntEventBranch("tammy_genie_n_pi_zeros",      0);
    declareIntEventBranch("tammy_genie_n_charms",        0);
    declareIntEventBranch("tammy_genie_n_kaons",         0);
    declareIntEventBranch("tammy_genie_n_others",        0);
    declareIntEventBranch("tammy_genie_n_particles",     0);
    
    declareIntEventBranch("tammy_intraNukeNParticles",                 0);
    declareIntEventBranch("tammy_intraNukeDeltaPlusPlusDecay",         0);
    declareIntEventBranch("tammy_intraNukeNeutronQuasiElasticScatter", 0);
    declareIntEventBranch("tammy_intraNukeOtherProcess",               0);
    declareDoubleEventBranch("tammy_intraNukeProtonMomentum",          0);
    declareContainerDoubleEventBranch("tammy_intraNukeProtonMomentumVec", 4, -1);
    
    declareIntEventBranch("tammy_NoInteractionVertex",        0);
    declareIntEventBranch("tammy_NullVertex",                 0);
    declareIntEventBranch("tammy_FailFidVolume",              0);
    declareIntEventBranch("tammy_FailOutTracks",              0);
    declareIntEventBranch("tammy_UnattachedProngsWithTracks", 0);
    declareIntEventBranch("tammy_CreatedShortTracks",         0);
    declareIntEventBranch("tammy_FailShortOutTrack",          0);
    declareIntEventBranch("tammy_FailRefitFidVolume",         0);
    declareIntEventBranch("tammy_FailExitingProng",           0);
    declareIntEventBranch("tammy_FailContainedProng",         0);
    
    declareIntEventBranch("tammy_muon_enters_front",   0);
    declareIntEventBranch("tammy_proton_enters_front", 0);
    
    declareIntEventBranch("tammy_vtx_fit_converged",         -1);
    declareDoubleEventBranch("tammy_vtx_fit_chi2",           -1);
    declareContainerDoubleEventBranch("tammy_fit_vtx",    3, -1);
    
    declareIntEventBranch("tammy_passVertexZCut",               -1);
    declareIntEventBranch("tammy_timeSlice",                    -1);
    declareIntEventBranch("tammy_classification",               -1);
    declareIntEventBranch("tammy_n_odClusters",                 -1);
    declareIntEventBranch("tammy_n_odClustersWithTimeCut",      -1);
    declareDoubleEventBranch("tammy_muonFuzzEnergy",             0);
    declareDoubleEventBranch("tammy_protonFuzzEnergy",           0);
    declareDoubleEventBranch("tammy_primaryVertexEnergy",        0);
    declareDoubleEventBranch("tammy_secondaryVertexEnergy",      0);
    declareDoubleEventBranch("tammy_endPointEnergy",             0);
    declareDoubleEventBranch("tammy_isolatedEnergy",             0);
    declareDoubleEventBranch("tammy_isolatedEnergy_tracker",     0);
    declareDoubleEventBranch("tammy_isolatedEnergy_ecal",        0);
    declareDoubleEventBranch("tammy_isolatedEnergy_hcal",        0);
    declareDoubleEventBranch("tammy_isolatedEnergy_targets",     0);
    declareDoubleEventBranch("tammy_hadronic_energy",            0);
    declareDoubleEventBranch("tammy_odEnergy",                   0);
    declareDoubleEventBranch("tammy_odEnergyWithTimeCut",        0);
    
    declareIntBranch(m_hypMeths, "tammy_targetID",                  -1);
    declareIntBranch(m_hypMeths, "tammy_targetZ",                   -1);
    declareIntBranch(m_hypMeths, "tammy_muon_minos_track",           0);
    declareIntBranch(m_hypMeths, "tammy_muon_minos_stub",            0);
    declareIntBranch(m_hypMeths, "tammy_muon_od_track",              0);
    declareIntBranch(m_hypMeths, "tammy_muon_side_ecal",             0);
    declareIntBranch(m_hypMeths, "tammy_muon_down_hcal",             0);
    declareIntBranch(m_hypMeths, "tammy_muon_trk_pat_history",      -1);
    declareIntBranch(m_hypMeths, "tammy_proton_trk_pat_history",    -1);
    declareIntBranch(m_hypMeths, "tammy_muon_n_values",              0);
    declareIntBranch(m_hypMeths, "tammy_proton_n_values",            0);
    declareIntBranch(m_hypMeths, "tammy_proton_kinked",             -1);
    declareIntBranch(m_hypMeths, "tammy_proton_odMatch",            -1);
    declareIntBranch(m_hypMeths, "tammy_pOK",                        0);
    declareIntBranch(m_hypMeths, "tammy_inside_minos_partial_plane", 0);
    
    declareDoubleBranch(m_hypMeths, "tammy_targetZPos",             -1);
    declareDoubleBranch(m_hypMeths, "tammy_muon_score",             -1);
    declareDoubleBranch(m_hypMeths, "tammy_muon_theta",             -1);
    declareDoubleBranch(m_hypMeths, "tammy_muon_thetaX",            -1);
    declareDoubleBranch(m_hypMeths, "tammy_muon_thetaY",            -1);
    declareDoubleBranch(m_hypMeths, "tammy_muon_phi",               -1);
    declareDoubleBranch(m_hypMeths, "tammy_muon_enu",               -1);
    declareDoubleBranch(m_hypMeths, "tammy_muon_q2",                -1);
    declareDoubleBranch(m_hypMeths, "tammy_calc_muon_p",            -1);
    declareDoubleBranch(m_hypMeths, "tammy_exit_muon_p",            -1);
    declareDoubleBranch(m_hypMeths, "tammy_proton_score",           -1);
    declareDoubleBranch(m_hypMeths, "tammy_proton_score1",          -1);
    declareDoubleBranch(m_hypMeths, "tammy_proton_score2",          -1);
    declareDoubleBranch(m_hypMeths, "tammy_proton_theta",           -1);
    declareDoubleBranch(m_hypMeths, "tammy_proton_thetaX",          -1);
    declareDoubleBranch(m_hypMeths, "tammy_proton_thetaY",          -1);
    declareDoubleBranch(m_hypMeths, "tammy_proton_phi",             -1);
    declareDoubleBranch(m_hypMeths, "tammy_proton_enu",             -1);
    declareDoubleBranch(m_hypMeths, "tammy_proton_q2",              -1);
    declareDoubleBranch(m_hypMeths, "tammy_proton_chi2_ndf",        -1);
    declareDoubleBranch(m_hypMeths, "tammy_proton_p_calCorrection", -1);
    declareDoubleBranch(m_hypMeths, "tammy_proton_p_visEnergy",     -1);
    declareDoubleBranch(m_hypMeths, "tammy_pion_score",             -1);
    declareDoubleBranch(m_hypMeths, "tammy_pion_score1",            -1);
    declareDoubleBranch(m_hypMeths, "tammy_pion_score2",            -1);
    declareDoubleBranch(m_hypMeths, "tammy_pion_chi2_ndf",          -1);
    declareDoubleBranch(m_hypMeths, "tammy_open_angle",             -1);
    declareDoubleBranch(m_hypMeths, "tammy_coplanarAngle",          -1);
    
    declareIntBranch(m_hypMeths,    "tammy_muon_odLastFrame",          -1);
    declareIntBranch(m_hypMeths,    "tammy_muon_odLastStory",          -1);
    declareDoubleBranch(m_hypMeths, "tammy_muon_odFaceX",              -9999);
    declareDoubleBranch(m_hypMeths, "tammy_muon_odFaceY",              -9999); 
    declareDoubleBranch(m_hypMeths, "tammy_muon_odFaceZ",              -9999);
    declareDoubleBranch(m_hypMeths, "tammy_muon_odEndX",               -9999);            
    declareDoubleBranch(m_hypMeths, "tammy_muon_odEndY",               -9999);
    declareDoubleBranch(m_hypMeths, "tammy_muon_odEndZ",               -9999);            
    declareDoubleBranch(m_hypMeths, "tammy_muon_odLastClusZ",          -9999);
    declareDoubleBranch(m_hypMeths, "tammy_muon_odStopDistMomentum",   -9999);            
    declareDoubleBranch(m_hypMeths, "tammy_muon_odElossMomentum",      -9999); 
    declareDoubleBranch(m_hypMeths, "tammy_muon_odTrackAvgTime",       -9999);
    declareDoubleBranch(m_hypMeths, "tammy_muotammy_n_odClustersAvgTime",    -9999);          
    
    declareIntBranch(m_hypMeths,    "tammy_isProtonInsideOD",        -1);
    declareIntBranch(m_hypMeths,    "tammy_isMuonInsideOD",          -1);
    declareIntBranch(m_hypMeths,    "tammy_ntrajProngProng",         -1);
    declareIntBranch(m_hypMeths,    "tammy_ntrajMuonProng",          -1);
    declareIntBranch(m_hypMeths,    "tammy_trajProtonProngPrimary",  -1);
    declareIntBranch(m_hypMeths,    "tammy_trajMuonProngPrimary",    -1);
    declareIntBranch(m_hypMeths,    "tammy_trajProtonProngPDG",      -1);
    declareIntBranch(m_hypMeths,    "tammy_trajMuonProngPDG",        -1);
    declareDoubleBranch(m_hypMeths, "tammy_trajMuonProngMomentum",   -1);
    declareDoubleBranch(m_hypMeths, "tammy_trajProtonProngMomentum", -1);
    declareDoubleBranch(m_hypMeths, "tammy_trajProtonTheta",         -1);
    declareDoubleBranch(m_hypMeths, "tammy_trajProtonPhi",           -1);
    declareDoubleBranch(m_hypMeths, "tammy_trajMuonTheta",           -1);
    declareDoubleBranch(m_hypMeths, "tammy_trajMuonPhi",             -1);
    declareDoubleBranch(m_hypMeths, "tammy_endProtonTrajMomentum",   -9999);
    declareDoubleBranch(m_hypMeths, "tammy_endMuonTrajMomentum",     -9999);
    declareDoubleBranch(m_hypMeths, "tammy_endProtonTrajXPosition",  -9999);
    declareDoubleBranch(m_hypMeths, "tammy_endProtonTrajYPosition",  -9999);
    declareDoubleBranch(m_hypMeths, "tammy_endProtonTrajZPosition",  -9999);
    declareDoubleBranch(m_hypMeths, "tammy_endMuonTrajXPosition",    -9999);
    declareDoubleBranch(m_hypMeths, "tammy_endMuonTrajYPosition",    -9999);
    declareDoubleBranch(m_hypMeths, "tammy_endMuonTrajZPosition",    -9999);
    
    declareContainerDoubleBranch(m_hypMeths, "tammy_minos_uv",             2, -1);
    declareContainerDoubleBranch(m_hypMeths, "tammy_muon_startPoint",      3, -1);
    declareContainerDoubleBranch(m_hypMeths, "tammy_muon_endPoint",        3, -1);
    declareContainerDoubleBranch(m_hypMeths, "tammy_proton_startPoint",    3, -1);
    declareContainerDoubleBranch(m_hypMeths, "tammy_proton_endPoint",      3, -1);
    declareContainerDoubleBranch(m_hypMeths, "tammy_proton_4p",            4, -1);
    declareContainerDoubleBranch(m_hypMeths, "tammy_proton_p_values",    200, -1);
    declareContainerDoubleBranch(m_hypMeths, "tammy_proton_chi2_values", 200, -1);
    declareContainerDoubleBranch(m_hypMeths, "tammy_muon_p_values",      200, -1);
    declareContainerDoubleBranch(m_hypMeths, "tammy_muon_chi2_values",   200, -1);
    
    declareContainerIntBranch(m_hypMeths, "tammy_primaryVtxPDG",   200, -1);
    declareContainerIntBranch(m_hypMeths, "tammy_secondaryVtxPDG", 200, -1);
    declareContainerIntBranch(m_hypMeths, "tammy_endPointVtxPDG",  200, -1);
    declareContainerIntBranch(m_hypMeths, "tammy_isolatedPDG",     200, -1);
    declareContainerIntBranch(m_hypMeths, "tammy_muonFuzzPDG",     200, -1);
    declareContainerIntBranch(m_hypMeths, "tammy_protonFuzzPDG",   200, -1);
    
    declareContainerIntBranch(m_hypMeths, "tammy_primaryVtxParentId",   200, -1);
    declareContainerIntBranch(m_hypMeths, "tammy_secondaryVtxParentId", 200, -1);
    declareContainerIntBranch(m_hypMeths, "tammy_endPointVtxParentId",  200, -1);
    declareContainerIntBranch(m_hypMeths, "tammy_isolatedParentId",     200, -1);
    declareContainerIntBranch(m_hypMeths, "tammy_muonFuzzParentId",     200, -1);
    declareContainerIntBranch(m_hypMeths, "tammy_protonFuzzParentId",   200, -1);
    
    declareContainerDoubleBranch(m_hypMeths, "tammy_primaryVtxEnergy",   200, -1);
    declareContainerDoubleBranch(m_hypMeths, "tammy_secondaryVtxEnergy", 200, -1);
    declareContainerDoubleBranch(m_hypMeths, "tammy_endPointVtxEnergy",  200, -1);
    declareContainerDoubleBranch(m_hypMeths, "tammy_isolatedEnergy",     200, -1);
    declareContainerDoubleBranch(m_hypMeths, "tammy_muonFuzzEnergy",     200, -1);
    declareContainerDoubleBranch(m_hypMeths, "tammy_protonFuzzEnergy",   200, -1);
    
    declareContainerDoubleBranch(m_hypMeths, "tammy_primaryVtxTrueEnergy",   200, -1);
    declareContainerDoubleBranch(m_hypMeths, "tammy_secondaryVtxTrueEnergy", 200, -1);
    declareContainerDoubleBranch(m_hypMeths, "tammy_endPointVtxTrueEnergy",  200, -1);
    declareContainerDoubleBranch(m_hypMeths, "tammy_isolatedTrueEnergy",     200, -1);
    declareContainerDoubleBranch(m_hypMeths, "tammy_muonFuzzTrueEnergy",     200, -1);
    declareContainerDoubleBranch(m_hypMeths, "tammy_protonFuzzTrueEnergy",   200, -1);
    
    declareDoubleEventBranch("tammy_isolatedEnergy_tracker",     0);
    declareDoubleEventBranch("tammy_isolatedEnergy_ecal",        0);
    declareDoubleEventBranch("tammy_isolatedEnergy_hcal",        0);
    declareDoubleEventBranch("tammy_isolatedEnergy_targets",     0);
    
    debug() << "Initializing Tammy's Variables Finished!" << endmsg;
    
    /*
    ============================================================================
        CCPi0 Variables - Trung's Package
    ============================================================================
    */
    
    debug() << "Initializing Trungs's Variables Started" << endmsg;
    
    //! Declare Truth Branches
    declareDoubleTruthBranch( "MC_pi0_energy", -9999 );
    declareContainerDoubleTruthBranch( "MC_pi0_momentum" );
    
    declareDoubleTruthBranch( "MC_scalar", -9999);
    declareDoubleTruthBranch( "MC_photon_energy_1", -9999 );
    declareDoubleTruthBranch( "MC_photon_energy_2", -9999 );
    declareContainerDoubleTruthBranch( "MC_photon_direction_1");
    declareContainerDoubleTruthBranch( "MC_photon_direction_2");
    declareContainerDoubleTruthBranch( "MC_photon_vertex_1" );
    declareContainerDoubleTruthBranch( "MC_photon_vertex_2" );
    
    
    // Declare Event Branches
    declareDoubleEventBranch( "Vertex_blob_energy", -9999 );
    declareDoubleEventBranch( "Filament_Vertex_energy", -9999 );
    declareDoubleEventBranch( "Sphere_Vertex_energy", -9999 );
    declareDoubleEventBranch( "Dispersed_blob_energy", -9999 );
    declareDoubleEventBranch( "Rejected_blob_vis_energy", -9999 );
    declareDoubleEventBranch( "Muon_blob_energy", -9999 );
    
    declareDoubleEventBranch( "RE_scalar", -9999);
    declareDoubleEventBranch( "RE_photon_energy_1", -9999 );
    declareDoubleEventBranch( "RE_photon_energy_2", -9999 );
    declareDoubleEventBranch( "RE_photon_dEdx_1", -9999 );
    declareDoubleEventBranch( "RE_photon_dEdx_2", -9999 );
    declareDoubleEventBranch( "RE_photon_time_1", -9999 );
    declareDoubleEventBranch( "RE_photon_time_2", -9999 );
    declareBoolEventBranch("is_GoodDirection1");
    declareBoolEventBranch("is_GoodPosition1");
    declareBoolEventBranch("is_GoodDirection2");
    declareBoolEventBranch("is_GoodPosition2");
    
    declareBoolEventBranch("is_GoodBlob1" );
    declareBoolEventBranch("is_GoodBlob2" );
    declareContainerDoubleEventBranch( "RE_photon_direction_1");
    declareContainerDoubleEventBranch( "RE_photon_direction_2");
    declareContainerDoubleEventBranch( "RE_photon_vertex_1");
    declareContainerDoubleEventBranch( "RE_photon_vertex_2");
    
    declareDoubleEventBranch( "RE_energy_Tracker", -9999 );
    declareDoubleEventBranch( "RE_energy_ECAL", -9999 );
    declareDoubleEventBranch( "RE_energy_HCAL", -9999 );
    
    
    declareIntEventBranch("minos_trk_is_contained",  -1);
    declareIntEventBranch("minos_trk_quality",       -1);
    declareIntEventBranch("minos_trk_is_ok",         -1);
    declareDoubleEventBranch("minos_trk_fit_pass",   -1);
    declareDoubleEventBranch("minos_trk_qp",         -1);
    declareDoubleEventBranch("minos_trk_eqp",        -1);
    declareDoubleEventBranch("minos_trk_p",          -1);
    declareDoubleEventBranch("minos_trk_p_range",    -1);
    declareDoubleEventBranch("minos_trk_p_curvature",-1);
    declareIntEventBranch("minos_trk_used_range",    -1);
    declareIntEventBranch("minos_trk_used_curvature",-1);
    declareIntEventBranch("minos_trk_end_plane",     -1);
    
    
    declareBoolEventBranch("is_anglescan" );
    declareBoolEventBranch("is_anglescan_applied");
    declareBoolEventBranch("is_houghtransform" );
    declareBoolEventBranch("is_houghtransform_applied");
    declareBoolEventBranch("is_twoDBlob" );
    
        // OD info
    declareDoubleEventBranch( "od_upstreamFrame", -9999 );
    declareDoubleEventBranch( "od_downstreamFrame", -9999 );
    declareDoubleEventBranch( "od_upstreamFrame_z", -9999 );
    declareDoubleEventBranch( "od_downstreamFrame_z", -9999 );
    
    declareDoubleEventBranch( "od_highStory", -9999 );
    declareDoubleEventBranch( "od_lowStory", -9999 );
    declareDoubleEventBranch( "od_highStory_t", -9999 );
    declareDoubleEventBranch( "od_lowStory_t", -9999 );
    
    declareDoubleEventBranch( "od_maxEnergy", -9999 );
    declareIntEventBranch( "od_energeticTower", -9999 );
    
    declareContainerDoubleEventBranch( "od_distanceBlobTower");
    declareContainerDoubleEventBranch( "od_towerEnergy" );
    declareContainerDoubleEventBranch( "od_towerNClusters" );
    declareContainerDoubleEventBranch( "od_towerTime" );
    
    declareContainerDoubleEventBranch( "od_idBlobTime" );
    declareContainerDoubleEventBranch( "od_towerTimeBlobOD" );
    declareContainerDoubleEventBranch( "od_towerTimeBlobMuon" );
    
    // Truth branch
    declareBoolTruthBranch( "reco_minos_match" );
    declareBoolTruthBranch( "is_fiducial" );
    declareBoolTruthBranch( "pass_plausible" );
    
    declareDoubleTruthBranch( "fslepton_E", -9999 );
    declareDoubleTruthBranch( "fslepton_P", -9999 );
    declareDoubleTruthBranch( "fslepton_T", -9999 );
    declareDoubleTruthBranch( "fslepton_theta", -9999 );
    declareDoubleTruthBranch( "fslepton_theta_x", -9999 );
    declareDoubleTruthBranch( "fslepton_theta_y", -9999 );
    declareDoubleTruthBranch( "fslepton_phi", -9999 );
    
    declareBoolTruthBranch( "is_ccpi0" );
    declareBoolTruthBranch( "is_cc1pi0");
    declareBoolTruthBranch( "is_ccpi0secondary" );
    declareBoolTruthBranch( "is_by_pim");
    declareBoolTruthBranch( "is_ccpi0x" );
    declareBoolTruthBranch( "is_other" );
    
    declareIntEventBranch( "survive_prefilter",   -1);
    declareIntEventBranch( "survive_do_muon",     -1);
    declareIntEventBranch( "survive_do_vertex",   -1);
    declareIntEventBranch( "survive_vtx_blob",    -1);
    declareIntEventBranch( "survive_all",         -1);
    declareIntEventBranch( "survive_minos_match", -1);
    declareIntEventBranch( "survive_has_vertex",  -1);
    declareIntEventBranch( "survive_fiducial",    -1);
    declareIntEventBranch( "tfiducial",           -1);
    
    declareIntEventBranch( "survive_plausible",-1);
    declareIntEventBranch( "survive_three_vertex",-1);
    declareIntEventBranch( "survive_onetrackpervtx",-1);
    declareIntEventBranch( "survive_gammatrack",-1);
    
    declareIntEventBranch( "primary_index", -1);
    declareIntEventBranch( "primary_multiplicity",-1);
    declareIntEventBranch( "primary_multiplicity2",-1); 
    declareIntEventBranch( "vertex_count",-1);
    declareIntEventBranch( "vertex_count2", -1);
    declareIntEventBranch( "discard_track_count",-1);
    declareContainerIntEventBranch( "multiplicities" );
    declareContainerDoubleEventBranch( "deviations" );
    
    declareContainerDoubleEventBranch("primary_separations");
    declareContainerDoubleEventBranch("primary_trklengths");
    declareContainerIntEventBranch("primary_truth_counts");
    declareContainerDoubleEventBranch("primary_truth_shareds");
    declareContainerDoubleEventBranch("primary_truth_fractions1");
    declareContainerDoubleEventBranch("primary_truth_fractions2");
    declareContainerDoubleEventBranch("primary_truth_fractions3");
    declareContainerIntEventBranch("primary_truth_pdgs1");
    declareContainerIntEventBranch("primary_truth_pdgs2");
    declareContainerIntEventBranch("primary_truth_pdgs3");
    
    declareDoubleEventBranch("pi0_evis_muon_blob");
    declareDoubleEventBranch("pi0_evis_vtx_blob");
    declareDoubleEventBranch("pi0_evis_outtime_blob");
    declareDoubleEventBranch("pi0_evis_dispersed_blob");
    declareDoubleEventBranch("pi0_evisfrac_muon_blob");
    declareDoubleEventBranch("pi0_evisfrac_vtx_blob");
    declareDoubleEventBranch("pi0_evisfrac_outtime_blob");
    declareDoubleEventBranch("pi0_evisfrac_dispersed_blob");
    
    declareIntEventBranch("blob_ndof_1");
    declareIntEventBranch("blob_ndof_2");
    declareDoubleEventBranch("blob_fval_1");
    declareDoubleEventBranch("blob_fval_2");
    
    
    declareIntEventBranch("npi0",-1);
    declareIntEventBranch("npip",-1);
    declareIntEventBranch("npim",-1);
    declareIntEventBranch("npipm",-1);
    declareIntEventBranch("nmeson",-1);
    declareIntEventBranch("npi02",-1);
    declareIntEventBranch("np",-1);
    declareIntEventBranch("nn",-1);
    declareDoubleEventBranch("pke",-1.);
    declareDoubleEventBranch("nke",-1.);
    
    declareIntEventBranch("dmode",-1);
    declareDoubleEventBranch("pienergy0",-1.);
    declareDoubleEventBranch("pitheta0", -1.);
    declareDoubleEventBranch("piphi0",   -1.);
    declareDoubleEventBranch("oangle0",  -1.);
    declareDoubleEventBranch("oangle0x", -1.);
    
    declareDoubleEventBranch("g1e0",     -1.);
    declareDoubleEventBranch("g1theta0", -1.);
    declareDoubleEventBranch("g1phi0",   -1.);
    declareDoubleEventBranch("g2e0",     -1.);
    declareDoubleEventBranch("g2theta0", -1.);
    declareDoubleEventBranch("g2phi0",   -1.);
    declareContainerDoubleEventBranch("pimom0");
    declareContainerDoubleEventBranch("g1mom0");
    declareContainerDoubleEventBranch("g2mom0");
    declareIntEventBranch("g1convidet",-1);       // Whether the gamma's converts inside the  
    declareIntEventBranch("g2convidet",-1);       // inner detector
    declareDoubleEventBranch("g1convdist",-1.0);
    declareDoubleEventBranch("g2convdist",-1.0);
    declareContainerDoubleEventBranch("g1convpos");
    declareContainerDoubleEventBranch("g2convpos");
    
    declareIntEventBranch("npim2", -1);          /* pi- variables */
    declareIntEventBranch("npimcapture", -1);
    declareIntEventBranch("npiminelastic", -1);
    declareIntEventBranch("npimdecay", -1);
    declareIntEventBranch("npip2", -1);
    declareIntEventBranch("npipcapture", -1);    /* pi+ variables */
    declareIntEventBranch("npipinelastic", -1);
    declareIntEventBranch("npipdecay", -1);
    declareDoubleEventBranch("pimlength", -1);
    declareDoubleEventBranch("piplength", -1);
    
    declareIntEventBranch("nmupdecay", -1);
    declareIntEventBranch("nmumdecay", -1);
    declareIntEventBranch("nmumcapture", -1);
    declareContainerDoubleEventBranch("michel_pos");
    declareContainerDoubleEventBranch("michel_mom");
    
    
    declareDoubleEventBranch("g1nukeedep", -1.);
    declareDoubleEventBranch("g2nukeedep", -1.);
    declareDoubleEventBranch("g1trkredep", -1.);
    declareDoubleEventBranch("g2trkredep", -1.);
    declareDoubleEventBranch("g1ecaledep", -1.);
    declareDoubleEventBranch("g2ecaledep", -1.);
    declareDoubleEventBranch("g1hcaledep", -1.);
    declareDoubleEventBranch("g2hcaledep", -1.);
    declareDoubleEventBranch("g1sideedep", -1.);
    declareDoubleEventBranch("g2sideedep", -1.);
    declareDoubleEventBranch("g1idetedep", -1.);
    declareDoubleEventBranch("g2idetedep", -1.);
    declareDoubleEventBranch("g1odetedep", -1.);
    declareDoubleEventBranch("g2odetedep", -1.);
    declareDoubleEventBranch("g1othersubdetedep", -1.);
    declareDoubleEventBranch("g2othersubdetedep", -1.);
    declareDoubleEventBranch("g1ecalo", -1.0);
    declareDoubleEventBranch("g2ecalo", -1.0);
    
    declareDoubleEventBranch("pi0nukeedep", -1.);
    declareDoubleEventBranch("pi0trkredep", -1.);
    declareDoubleEventBranch("pi0ecaledep", -1.);
    declareDoubleEventBranch("pi0hcaledep", -1.);
    declareDoubleEventBranch("pi0sideedep", -1.);
    declareDoubleEventBranch("pi0idetedep", -1.);
    declareDoubleEventBranch("pi0odetedep", -1.);
    declareDoubleEventBranch("pi0othersubdetedep", -1.);
    declareDoubleEventBranch("pi0ecalo",-1.0);
    
    declareDoubleEventBranch("neutronnukeedep", -1.);
    declareDoubleEventBranch("neutrontrkredep", -1.);
    declareDoubleEventBranch("neutronecaledep", -1.);
    declareDoubleEventBranch("neutronhcaledep", -1.);
    declareDoubleEventBranch("neutronsideedep", -1.);
    declareDoubleEventBranch("neutronidetedep", -1.);
    declareDoubleEventBranch("neutronodetedep", -1.);
    declareDoubleEventBranch("neutronothersubdetedep", -1.);
    declareDoubleEventBranch("neutronecalo",-1.0);
    
    declareDoubleEventBranch("totalevis", -1.0);
    declareDoubleEventBranch("ntgtevis",  -1.0);
    declareDoubleEventBranch("trkrevis",  -1.0);
    declareDoubleEventBranch("ecalevis",  -1.0);
    declareDoubleEventBranch("hcalevis",  -1.0);
    declareDoubleEventBranch("otherevis", -1.0);
    
    
    declareIntEventBranch("nblob_anglescan",-1);
    declareIntEventBranch("nblob_hough", -1);
    declareIntEventBranch("anglescan_ncandx", -1);
    declareIntEventBranch("anglescan_ncand", -1);
    
    
    declareContainerIntEventBranch("anglescan_candx_nc");
    declareContainerIntEventBranch("anglescan_candx_nd");
    
    declareContainerIntEventBranch("anglescan_cand_nc");
    declareContainerIntEventBranch("anglescan_cand_ncx");
    declareContainerIntEventBranch("anglescan_cand_ncu");
    declareContainerIntEventBranch("anglescan_cand_ncv");
    declareContainerIntEventBranch("anglescan_cand_nd");
    declareContainerIntEventBranch("anglescan_cand_ndx");
    declareContainerIntEventBranch("anglescan_cand_ndu");
    declareContainerIntEventBranch("anglescan_cand_ndv");
    
    declareContainerIntEventBranch("anglescan_blob_nc");
    declareContainerIntEventBranch("anglescan_blob_ncx");
    declareContainerIntEventBranch("anglescan_blob_ncu");
    declareContainerIntEventBranch("anglescan_blob_ncv");
    declareContainerIntEventBranch("anglescan_blob_nd");
    declareContainerIntEventBranch("anglescan_blob_ndx");
    declareContainerIntEventBranch("anglescan_blob_ndu");
    declareContainerIntEventBranch("anglescan_blob_ndv");
    
    declareContainerIntEventBranch("hough_blob_nc");
    declareContainerIntEventBranch("hough_blob_ncx");
    declareContainerIntEventBranch("hough_blob_ncu");
    declareContainerIntEventBranch("hough_blob_ncv");
    declareContainerIntEventBranch("hough_blob_nd");
    declareContainerIntEventBranch("hough_blob_ndx");
    declareContainerIntEventBranch("hough_blob_ndu");
    declareContainerIntEventBranch("hough_blob_ndv");
    
    declareContainerIntEventBranch("final_blob_nc");
    declareContainerIntEventBranch("final_blob_ncx");
    declareContainerIntEventBranch("final_blob_ncu");
    declareContainerIntEventBranch("final_blob_ncv");
    declareContainerIntEventBranch("final_blob_nd");
    declareContainerIntEventBranch("final_blob_ndx");
    declareContainerIntEventBranch("final_blob_ndu");
    declareContainerIntEventBranch("final_blob_ndv");
    
    
    
    declareIntEventBranch("g1blob_ndigit",-1);
    declareIntEventBranch("g1blob_ncluster",-1);
    declareDoubleEventBranch("g1blob_minsep", -1.0);
    declareDoubleEventBranch("g1blob_vtx_distance",-1.0);
    declareDoubleEventBranch("g1blob_edge_distance",-1.0);
    
    declareIntEventBranch("g2blob_ndigit",-1);
    declareIntEventBranch("g2blob_ncluster",-1);
    declareDoubleEventBranch("g2blob_minsep", -1.0);
    declareDoubleEventBranch("g2blob_vtx_distance",-1.0);
    declareDoubleEventBranch("g2blob_edge_distance",-1.0);
    

    declareDoubleEventBranch("pienergy",1.e6);
    declareDoubleEventBranch("pitheta", 1.e6);
    declareDoubleEventBranch("piphi",   1.e6);
    declareDoubleEventBranch("oangle",  1.e6);
    
    declareDoubleEventBranch("g1e",     1.e6);
    declareDoubleEventBranch("g1theta", 1.e6);
    declareDoubleEventBranch("g1phi",   1.e6);
    declareDoubleEventBranch("g2e",     1.e6);
    declareDoubleEventBranch("g2theta", 1.e6);
    declareDoubleEventBranch("g2phi",   1.e6);
    declareContainerDoubleEventBranch("pimom");
    declareContainerDoubleEventBranch("g1mom");
    declareContainerDoubleEventBranch("g2mom");
    
    declareDoubleEventBranch("mgg", 1.e6);
    
    declareContainerDoubleEventBranch("mumom");
    declareDoubleEventBranch("Tn", -1.0);
    declareDoubleEventBranch("Tn2", -1.0);
    declareDoubleEventBranch("Erec", -1.0);
    declareDoubleEventBranch("Erec2", -1.0);
    declareDoubleEventBranch("Q2", -1.0);
    declareDoubleEventBranch("W2", -1.0);
    declareDoubleEventBranch("W",  -1.0);
    
    // Truth info of reconstructed EM showers
    declareDoubleEventBranch("g1pi0evis", -1.0);
    declareDoubleEventBranch("g1g1evis", -1.0);
    declareDoubleEventBranch("g1g2evis", -1.0);
    declareDoubleEventBranch("g1neutronevis", -1.0);
    declareDoubleEventBranch("g1protonevis", -1.0);
    declareDoubleEventBranch("g1pipevis", -1.0);
    declareDoubleEventBranch("g1pimevis", -1.0);
    declareDoubleEventBranch("g1gmevis", -1.0);
    declareDoubleEventBranch("g1muevis", -1.0);
    declareDoubleEventBranch("g1otherevis", -1.0);
    declareDoubleEventBranch("g1totalevis", -1.0);
    declareDoubleEventBranch("g1sharedevis", -1.0);
    declareDoubleEventBranch("g1mostevisfrac", -1.0);
    declareIntEventBranch("g1mostevispdg",-1);
    
    declareDoubleEventBranch("g2pi0evis", -1.0);
    declareDoubleEventBranch("g2g1evis", -1.0);
    declareDoubleEventBranch("g2g2evis", -1.0);
    declareDoubleEventBranch("g2neutronevis", -1.0);
    declareDoubleEventBranch("g2protonevis", -1.0);
    declareDoubleEventBranch("g2pipevis", -1.0);
    declareDoubleEventBranch("g2pimevis", -1.0);
    declareDoubleEventBranch("g2gmevis", -1.0);
    declareDoubleEventBranch("g2muevis", -1.0);
    declareDoubleEventBranch("g2otherevis", -1.0);
    declareDoubleEventBranch("g2totalevis", -1.0);
    declareDoubleEventBranch("g2sharedevis", -1.0);
    declareDoubleEventBranch("g2mostevisfrac", -1.0);   /* fraction of evis */
    declareIntEventBranch("g2mostevispdg",-1);          /* pdg of the particle contributing most */
    
    // NeutrinoInt branches 
    //declareDoubleBranch( m_hypMeths, "muon_P", -9999 );
    //declareDoubleBranch( m_hypMeths, "muon_E", -9999 );
    //declareDoubleBranch( m_hypMeths, "muon_T", -9999 );
    //declareDoubleBranch( m_hypMeths, "muon_Pscale", -9999 );
    
    debug() << "Initializing Trungs's Variables Finished!" << endmsg;
    
    info() <<"Exit CCDeltaPlusAna::initialize()" << endmsg;
    
    return sc;
}
    
//=============================================================================
// reconstructEvent() --
//=============================================================================
StatusCode CCDeltaPlusAna::reconstructEvent( Minerva::PhysicsEvent *event, Minerva::GenMinInteraction* truth ) const
{

    debug() <<"Enter CCDeltaPlusAna::reconstructEvent()" << endmsg;
    debug() <<"=============================================================================" <<endmsg;
    debug() <<" "<<endmsg;
    
    Minerva::NeutrinoInt *nuInt = new Minerva::NeutrinoInt( "CCDeltaPlusAna" );
    NeutrinoVect defaultHyp;
    defaultHyp.push_back(nuInt);
    addInteractionHyp(event,defaultHyp);
    
    
    //-- check if this a plausible event ( MC only )
    if( truth && !truthIsPlausible(truth) ) {
        debug() << "   This is not a plausible MC event! returning!" << endmsg;
        return StatusCode::SUCCESS;
    }
    
    //-- fill in some truth information ( at the moment, prefer not to move information to the Truth ntuple )
//     if( truth ) setGenMinTruthInformation(event,truth);
    
    //-- interaction vertex check
    if( !event->hasInteractionVertex() ) {
        debug() << "  The event does not have an interaction vertex!" << endmsg;
        event->setIntData("tammy_NoInteractionVertex",1);
        return interpretFailEvent(event);
    } 
    
    //-- get the interaction vertex 
    SmartRef<Minerva::Vertex> vertex = event->interactionVertex();
    if( !vertex ) { 
        bool pass = true; 
        std::string tag = "BadObject";
        event->filtertaglist()->addFilterTag(tag,pass);
        event->setIntData("tammy_NullVertex",1);
        error() << "This vertex is NULL! Flag this event as bad!" << endmsg;
        return interpretFailEvent(event);
    }
    
    //-- check if vertex is inside passive or active targets fiducial volume 
//     int tammy_passVertexZCut = 0;
//     std::string name   = "";
//     if( !isInsideTargetFiducial(vertex,name,tammy_passVertexZCut) ) { 
//         debug() << "  The vertex is not within any of the targets fiducial volume!" << endmsg; 
//         event->setIntData("tammy_FailFidVolume",1);
//         return interpretFailEvent(event);
//     }else{
//         debug() << "  The vertex is inside " << name << " fiducial volume with position = " << vertex->position() << endmsg;
//     }
//     event->setIntData(name,1);
    
    
    //-- check the number of outgoing tracks 
    
    debug() << "  The number of outgoing tracks is = " << vertex->getOutgoingTracks().size() << endmsg;
    if( vertex->getOutgoingTracks().size() < 2 ) {
        debug() << "  This vertex has " << vertex->getOutgoingTracks().size() << " outgoing tracks!" << endmsg;
        event->setIntData("tammy_FailOutTracks",1);
        return interpretFailEvent(event);
    } else{
        debug() << "  The number of outgoing tracks is = " << vertex->getOutgoingTracks().size() << endmsg;
    }
    
    
    //-- check the number of unattached prongs in the event with tracks
    int ntracks = 0;
    Minerva::ProngVect unattachedProngs = event->select<Minerva::Prong>("Used:Unused","!IsSubProng");
    for(unsigned int i = 0; i < unattachedProngs.size(); i++){
        if( !unattachedProngs[i]->minervaTracks().empty() ){
            ntracks++;
        }
    }
    if( ntracks != 0 ) {
        debug() << "  This event has unattached prongs with tracks!" << endmsg;
        event->setIntData("tammy_UnattachedProngsWithTracks",1);
//         return interpretFailEvent(event);
    } else{
        debug() << "  The number of unattached prong is = " << unattachedProngs.size() << endmsg;
    }
    
    //-- make tracks using the anchor short tracker and refit the vertex
    bool createdTracks = createdAnchoredShortTracks(event,vertex);
    if( createdTracks ) {
        event->setIntData("tammy_CreatedShortTracks", 1);
        debug() << "  Created short anchored track(s)!" << endmsg;
    } else{
        debug() << "  Did not create any short tracks." << endmsg;
    }
    
    //-- re-check the number of outgoing tracks
    debug() << "  After createdAnchoredShortTracks() the number of outgoing primary tracks: " << vertex->getOutgoingTracks().size() << endmsg;
//     if( vertex->getOutgoingTracks().size() != 2 ) {
//         debug() << "  The number of outgoing tracks after creating short tracks = " << vertex->getOutgoingTracks().size() << endmsg;
//         event->setIntData("tammy_FailShortOutTrack",1);
//         return interpretFailEvent(event);
//     } else{
//         debug() << "  The vertex has 2 outgoing primary tracks!" << endmsg;
//     }
    
    //-- set vertex fit information
    double fit_chi2 = 0;
    if( vertex->hasDoubleData("fit_chi2") ) {
        fit_chi2 = vertex->getDoubleData("fit_chi2"); 
        event->setDoubleData("tammy_vtx_fit_chi2",fit_chi2);
        debug() << "   The vertex fit chi2 = " << fit_chi2 << endmsg;  
    }
    
    int fit_converged = 0;
    if( vertex->hasIntData("fit_converged") ) {
        fit_converged = vertex->getIntData("fit_converged");
        event->setIntData("tammy_vtx_fit_converged",fit_converged);
        debug() << "   The vertex fit convergence = " << fit_converged << endmsg;
    }
    
    //-- re-set the vertex position for vertex fit failures for long-long tracks combination
    SmartRefVector<Minerva::Track> trks = vertex->getOutgoingTracks();
    if( fit_converged == 0 ) {
        if( (trks[0]->patRecHistory() == Minerva::Track::LongPatRec3View && trks[1]->patRecHistory() == Minerva::Track::LongPatRec3View) ||
            (trks[0]->patRecHistory() == Minerva::Track::LongPatRec3View && trks[1]->patRecHistory() == Minerva::Track::LongPatRec2View) ||
            (trks[0]->patRecHistory() == Minerva::Track::LongPatRec2View && trks[1]->patRecHistory() == Minerva::Track::LongPatRec3View) ) {
            
            m_vertexFitter->fit(vertex);
        
            std::vector<double> tammy_fit_vtx;
            tammy_fit_vtx.push_back( vertex->position().x() );
            tammy_fit_vtx.push_back( vertex->position().y() );
            tammy_fit_vtx.push_back( vertex->position().z() );
        
            event->setContainerDoubleData("fit_vtx",tammy_fit_vtx);
            event->setInteractionVertex(vertex);
            Minerva::ProngVect prongs = event->primaryProngs();
            
            for(unsigned int p = 0; p < prongs.size(); p++){
                prongs[p]->setSourceVertex(vertex);
            }
        }
    }
    
    //-- make sure the update vertex is inside the targets fiducial volume
//     if( createdTracks ) {
//         name.clear();
//         tammy_passVertexZCut = 0;
//         vertex = event->interactionVertex();
//         
//         if( !isInsideTargetFiducial(vertex,name,tammy_passVertexZCut) ) {
//             debug() << "  The update vertex with position = " << vertex->position() << " is not within fiducial volume!" << endmsg;
//             event->setIntData("tammy_FailRefitFidVolume",1);
//             return interpretFailEvent(event);
//         }else{
//             debug() << "  The update vertex with position = " << vertex->position() << " is inside the fiducial volume!" << endmsg;
//         }
//         event->setIntData(name,1);
//     }
    
    //-- make new prongs to include the new created primary track
    bool makeProngs = false;
    if( createdTracks ){
        makeProngs = createdTrackedProngs(event);
    }
    debug() << "  Was the creation of a new prong successful? " << makeProngs << endmsg;
    
    //-- check if kinked tracks can be made; add them and linked vertices to the prong
    bool createdKinkedTracks = createdKinkedShortTracks(event);
    debug() << "  Was the creation of a kinked track(s) successful? " << createdKinkedTracks << endmsg;
    
    //-- (re) classify track-based prongs
    if( makeProngs || createdKinkedTracks ) {
        bool classify = classifyProngs(event);
        debug() << "  Was the tammy_classification of the track-based prongs successful? " << classify << endmsg;
    }
    
    //-- get all of the primary prongs in the event
    Minerva::ProngVect primaryProngs = event->primaryProngs();
    debug() << "  The number of primary prongs = " << primaryProngs.size() << endmsg;
    
    //-- create new particles ( non-Minos-bit positive particles )
    bool makeParticles = createdTrackedParticles(primaryProngs);
    debug() << "  Was the creation of particles successful? " << makeParticles << endmsg;
    
    //-- check if the other prong is not contained and has a muon particle association
    std::string eventType = "";
    SmartRef<Minerva::Prong>    muon;
    SmartRef<Minerva::Particle> muonParticle;
    if( !muonProng(primaryProngs,muon,muonParticle,eventType) ) {
        debug() << "  Didn't find any not exiting in the inner detector bit-positive prong with a muon particle!" << endmsg;
        event->setIntData("tammy_FailExitingProng",1);
        return interpretFailEvent(event);
    } else{
        debug() << "  Tag the muon's prong with bit-field = " << muon->typeBitsToString() << endmsg;
    }
    
    //-- check if one of the primary prong is contained and has a dEdX proton particle 
    SmartRef<Minerva::Prong>    proton;
    SmartRef<Minerva::Particle> protonParticle;
    primaryProngs = event->primaryProngs();
    if( !protonProng(primaryProngs,proton,protonParticle) ) {
        debug() << "  Didn't find any contained in the tracker bit-positive prong with a proton particle!" << endmsg;
        event->setIntData("tammy_FailContainedProng",1);
        return interpretFailEvent(event);
    } else{
        debug() << "  Tag the proton's prong with bit-field = " << proton->typeBitsToString() << endmsg;
    }
    
    //-- check event tammy_classification
//     if( eventType == "golden" ){
//         event->setIntData("tammy_classification",1);
//     }else if( eventType == "silver" ){
//         event->setIntData("tammy_classification",2);
//     }else { 
//         error() << "The event was not classified!" << endmsg; 
//         return interpretFailEvent(event); 
//     }
    
    //-- tag if the muon and proton candidates enter front of the detector
//     event->setIntData("tammy_muon_enters_front",int(muon->EnterFrontID()));
//     event->setIntData("tammy_proton_enters_front",int(proton->EnterFrontID()));
    
    //-- tag the tammy_timeSlice
    event->setIntData("tammy_timeSlice",int((muon->minervaTracks()[0])->timeSlice()));
    
    //-- tag if the reconstructed vertex passes the vtx cut
//     event->setIntData("tammy_passVertexZCut",tammy_passVertexZCut);
    
    //-- make idblob based prongs
    if( createdIDBlobProng(event) ) {
        classifyProngs(event);
        debug() << "  The creation of idblob based prong(s) was successful? " << endmsg; 
    } else{
        debug() << "  No idblob based prongs were created. " << endmsg;
    }
    
    //-- make odblob based prongs
    if( createdODBlobProng(event) ) {
        classifyProngs(event);
        debug() << "  The creation of odblob based prong(s) was successful? " << endmsg;
    } else{
        debug() << "  No odblob based prong were created." << endmsg; 
    }
    
    //-- look for michels
    if( tagMichelElectrons(event) ) {
        debug() << "  A Michel electron was found!" << endmsg;
    } else{
        debug() << "  Did not find a Michel electron." << endmsg;
    }
    
    //-- fill the event anatuple data
    fillBlobEventData(event);
    
    //-- call the interpret event function;
    NeutrinoVect nuInts;
    tammy_interpretEvent(event,truth,nuInts); 
    
    //-- mark the event
//     markEvent(event);
    
    //-- add neutrino interactions to the physics event
    addInteractionHyp(event,nuInts);
    

    
    debug() << "Exit NukeCCQE::reconstructEvent" << endmsg;


    Minerva::GenMinHeader* header(NULL);
    if (exist<Minerva::GenMinHeader>(Minerva::GenMinHeaderLocation::Default)) {
        header = get<Minerva::GenMinHeader>(Minerva::GenMinHeaderLocation::Default); 
    }

    if (header) {
        int run    = header->RunNumber();
        int subrun = header->SubRunNumber();
        int gate   = header->SpillNumber();
        if (fProcessOneEvent && ((run != fRunNumber) || (subrun != fSubrunNumber) || (gate != fGateNumber))) {
            return StatusCode::SUCCESS;
        }
    }
    
    fTrajectoryMap.clear();

    FillTrajectoryMap();

    SummarizeTruth(event);
    SummarizeHadronTruth(event);
    
    if (fPizero) GammaEdepInfo(event);
    
    
    
    debug() << gateData() << endmsg;
    
    markEvent(event);
    
    
    if ( !DoMuon(event, truth, vertex) )   return StatusCode::SUCCESS;
    event->setIntData("survive_do_muon", 1);
    
    m_PrimaryVertex = vertex;
    
    if ( !DoVertex(event, vertex) ) return StatusCode::SUCCESS;
    event->setIntData("survive_do_vertex", 1);
    
    if ( !PreFilter(event) ) return StatusCode::SUCCESS;
    event->setIntData("survive_prefilter", 1);
    
    if ( !VtxBlob(event, vertex) ) return StatusCode::SUCCESS;
    event->setIntData("survive_vtx_blob", 1);
    
    if ( !ConeBlobs(event, vertex) ) return StatusCode::SUCCESS;
    
    if ( !DispersedBlob(event) ) return StatusCode::SUCCESS;
    event->setIntData("survive_all",1);
    
    

        //markEvent( event );
    std::cout << "Passing all cuts " << std::endl;
    // Now interpret the event and add NeutrinoInts
    NeutrinoVect interactions;
    interpretEvent( event, truth, interactions );
    
    // Add the newly create NeutrinoInts to this PhysicsEvent
    StatusCode sc = addInteractionHyp( event, interactions );
    
    
    
    //-- fill AnaNTuple
    fillCommonPhysicsAnaBranches(event);
    
    //-- fill NuMI branches
    fillNuMIBranches(event);
    
    //-- fill vertex activity branches
    fillVertexActivityStudyBranches(event);
    
    //-- fill particle response branches
    if( haveNeutrinoMC() ) {
        unattachedProngs.clear();
        unattachedProngs = event->select<Minerva::Prong>("Used:Unused","All");
        for(unsigned int p = 0; p < unattachedProngs.size(); p++) {
            if( unattachedProngs[p]->filtertaglist()->filterTagExists("UnattachedEnergy") ) {
                Minerva::IDClusterVect clusters = unattachedProngs[p]->getAllIDClusters();
                fillParticleResponseBranches(event,clusters);
            }
        }
    }
    
    
    
    return StatusCode::SUCCESS;
    
}
    
    
StatusCode CCDeltaPlusAna::tammy_interpretEvent( const Minerva::PhysicsEvent *event, const Minerva::GenMinInteraction *truth, NeutrinoVect& nuInts ) const
{
    debug() << "Enter tammy_interpretEvent()" << endmsg;
    StatusCode sc;
    
    //-- initialize
    if( !truth ) info() << "  This is a MC event." << endmsg;
    
    //-- get the interaction vertex
    SmartRef<Minerva::Vertex> vertex = event->interactionVertex();
    
    //-- get the prongs
    Minerva::ProngVect primaryProngs    = event->primaryProngs();
    Minerva::ProngVect unattachedProngs = event->select<Minerva::Prong>("Used:Unused","All");
    debug() << "  The number of primary and unattached prongs = ( " 
            <<  primaryProngs.size() << ", " << unattachedProngs.size() << " )" << endmsg;
    
    //-- data storage
    SmartRef<Minerva::Prong> muonProng;
    SmartRef<Minerva::Prong> protonProng;
    
    SmartRef<Minerva::Particle> muonParticle;
    SmartRef<Minerva::Particle> protonParticle;
    
    Minerva::ProngVect primaryVertexProngs;
    Minerva::ProngVect secondaryVertexProngs;
    Minerva::ProngVect endPointProngs;
    Minerva::ProngVect isolatedProngs;
    
    //-- loop over primaryProngs and retrieve best particle
    for(unsigned int prong = 0; prong < primaryProngs.size(); prong++) {
    
        //-- get prong blobs
        Minerva::IDBlobVect blobs = primaryProngs[prong]->idblobs();
    
        //-- proton prong
        if( primaryProngs[prong]->filtertaglist()->filterTagExists("PrimaryProton") ) {
            protonProng    = primaryProngs[prong];
            protonParticle = primaryProngs[prong]->bestParticle();
            m_hitTagger->applyColorTag(protonProng,m_protonProngColor);
            for(unsigned int b = 0; b < blobs.size(); b++) {
            if( blobs[b]->history() == Minerva::IDBlob::Hidden ) m_hitTagger->applyColorTag(blobs[b],m_protonProngColor);
            }
        }
    
        //-- muon prong
        if( primaryProngs[prong]->filtertaglist()->filterTagExists("PrimaryMuon") ) {
        muonProng    = primaryProngs[prong];
        muonParticle = primaryProngs[prong]->bestParticle(); 
        m_hitTagger->applyColorTag(muonProng,m_muonProngColor);
        for(unsigned int b = 0; b < blobs.size(); b++) {
            if( blobs[b]->history() == Minerva::IDBlob::Hidden ) m_hitTagger->applyColorTag(blobs[b],m_muonProngColor);
        }
        }
    
        //-- vertex energy prong
        if( primaryProngs[prong]->filtertaglist()->filterTagExists("PrimaryVertexEnergy") ) {
        primaryVertexProngs.push_back( primaryProngs[prong] );
        m_hitTagger->applyColorTag(primaryProngs[prong],m_primaryVertexColor);
        }
    }
    
    //-- loop over unattached prongs
    for(unsigned int prong = 0; prong < unattachedProngs.size(); prong++) {
        if( unattachedProngs[prong]->filtertaglist()->filterTagExists("UnattachedEnergy") ) {
        isolatedProngs.push_back( unattachedProngs[prong] );
        m_hitTagger->applyColorTag(unattachedProngs[prong],m_unattachedProngColor);
        } else if( unattachedProngs[prong]->filtertaglist()->filterTagExists("SecondaryVertexEnergy") ) {
        secondaryVertexProngs.push_back( unattachedProngs[prong] );
        m_hitTagger->applyColorTag(unattachedProngs[prong],m_secondaryVertexColor);
        } else if( unattachedProngs[prong]->filtertaglist()->filterTagExists("EndPointVertexEnergy") ) {
        endPointProngs.push_back( unattachedProngs[prong] );
        m_hitTagger->applyColorTag(unattachedProngs[prong],m_endPointVertexColor);
        } 
    }
    
    //-- constructed the nuke two-track cc-qe neutrino hypothesis
    Minerva::NeutrinoInt* ccqeHyp = new Minerva::NeutrinoInt(m_anaSignature);
    ccqeHyp->setNeutrinoFlavor(Minerva::NeutrinoInt::MuonFlavor);
    ccqeHyp->setInteractionCurrent(Minerva::NeutrinoInt::ChargedCurrent);
    ccqeHyp->setInteractionType(Minerva::NeutrinoInt::QuasiElastic);
    
    //-- vertex information
    Gaudi::XYZTVector vtx;
    vtx.SetCoordinates( vertex->position().x(), vertex->position().y(), vertex->position().z(), muonParticle->startPos().t() );
    ccqeHyp->setVertex(vtx);
    double vertexZ = vertex->position().z();
    
    //-- track information
    Minerva::TrackVect muonTracks = muonProng->minervaTracks();
    Minerva::Track::PatRecHistory muonPatRecHistory = muonTracks[0]->patRecHistory();
    ccqeHyp->setIntData("tammy_muon_trk_pat_history",int(muonPatRecHistory));
    debug() << "  muon's track pattern recognition history = " << muonPatRecHistory << endmsg;
    
    Minerva::TrackVect protonTracks = protonProng->minervaTracks();
    Minerva::Track::PatRecHistory protonPatRecHistory = protonTracks[0]->patRecHistory();
    ccqeHyp->setIntData("tammy_proton_trk_pat_history",int(protonPatRecHistory));
    debug() << "  proton's track pattern recognition history = " << protonPatRecHistory << endmsg; 
    
    //-- muon and proton open angle
    double term1 = sin(muonTracks[0]->theta())*sin(protonTracks[0]->theta())*cos(muonTracks[0]->phi()-protonTracks[0]->phi());
    double term2 = cos(muonTracks[0]->theta())*cos(protonTracks[0]->theta());
    double tammy_open_angle =  acos( term1 + term2 );
    ccqeHyp->setDoubleData("open_angle",tammy_open_angle);
    debug() << "  open angle = " << tammy_open_angle*180/3.14159 
            << ", muon theta = " << muonTracks[0]->theta()*180/3.14159
            << ", proton theta = " << protonTracks[0]->theta()*180/3.14159 << endmsg;
    
    //-- passive target/materials
    int tammy_targetID  = -1;
    int materialZ = -1;
    
    Minerva::NuclearTarget* target = NULL; 
    for(unsigned int i = 0; i < m_targetNames.size(); i++) {
        if( m_targetNames[i] == "NuclearTargets" ) continue;
        if( !event->hasIntData(m_targetNames[i]) ) continue;
        target = m_nuclearTargetTool->getNuclearTarget(m_targetNames[i]); 
        tammy_targetID = target->getTargetID();
        debug() << "  retrieve the target name = " << target->getName() << ", and id = " << tammy_targetID << endmsg;  
    }
    
    if( target ) {
        const Material* material =  m_nuclearTargetTool->getSectionMaterial(target,vertex->position().x(),vertex->position().y());
        if( material ) {
        materialZ = (int)material->Z();
        debug() << "  retrieve the target's section material name = " << material->name() << ", and Z = " << material->Z() << endmsg;
        }
    }
    
    ccqeHyp->setIntData( "targetID",tammy_targetID );
    ccqeHyp->setIntData( "tammy_targetZ",materialZ );
    
    double tammy_targetZ = -1;
    if( tammy_targetID >= 1 && tammy_targetID <= 5 ) 
        tammy_targetZ = m_nuclearTargetTool->getSectionZCenter(target,vertex->position().x(),vertex->position().y());
    
    ccqeHyp->setDoubleData( "tammy_targetZPos",tammy_targetZ );
    debug() << "  the target's section z position = " << tammy_targetZ << endmsg;
    
    //-- prongs' truth information
    if( haveNeutrinoMC() ) {
        Minerva::ProngVect muonProngVect, protonProngVect;
        muonProngVect.push_back( muonProng );
        protonProngVect.push_back( protonProng );
    
        setTrackProngTruth(ccqeHyp,muonProng);
        debug() << "	set the muon prong track truth information" << endmsg;
    
        setTrackProngTruth(ccqeHyp,protonProng);
        debug() << "	set the proton prong track truth information" << endmsg;
    
        setBlobProngTruth(ccqeHyp,muonProngVect);
        debug() << "	set the muon prong blob truth information" << endmsg;
    
        setBlobProngTruth(ccqeHyp,protonProngVect);
        debug() << "	set the proton prong blob truth information" << endmsg;
    
        setBlobProngTruth(ccqeHyp,primaryVertexProngs);
        debug() << "	set the primary vertex blob truth information" << endmsg;
    
        setBlobProngTruth(ccqeHyp,secondaryVertexProngs);
        debug() << "	set the secondary vertex blob truth information" << endmsg;
    
        setBlobProngTruth(ccqeHyp,endPointProngs);
        debug() << "	set the track endpoint blob truth information" << endmsg;
    
        setBlobProngTruth(ccqeHyp,isolatedProngs);
        debug() << "	set the unattached energy blob truth information" << endmsg;
    }
    
    //! fill the muon minos branches
    fillMinosMuonBranches(ccqeHyp,muonProng);
    
    //-- kinematics from muon and proton
    int charge = -1;
    MuonUtils->muonCharge(muonProng,charge); 
    
    //-- neutrino helicity
    if( charge > 0 ) ccqeHyp->setNeutrinoHelicity( Minerva::NeutrinoInt::Antineutrino );
    else if( charge < 0 ) ccqeHyp->setNeutrinoHelicity( Minerva::NeutrinoInt::Neutrino );
    else { ccqeHyp->setNeutrinoHelicity( Minerva::NeutrinoInt::UnknownHelicity ); charge = -1.0; }
    
    //-- rotated muon's angle w.r.t beam
    double tammy_muon_theta  = m_coordSysTool->thetaWRTBeam(muonParticle->momentumVec());
    double tammy_muon_thetaX = m_coordSysTool->thetaXWRTBeam(muonParticle->momentumVec());
    double tammy_muon_thetaY = m_coordSysTool->thetaYWRTBeam(muonParticle->momentumVec());
    double tammy_muon_phi    = m_coordSysTool->phiWRTBeam(muonParticle->momentumVec());
    
    //-- corrected muon energy loss in target
    Gaudi::LorentzVector muonfourVec;
    sc = m_energyCorrectionTool->getCorrectedEnergy(muonProng,muonParticle,vertexZ,muonfourVec);  
    if( !sc ) muonfourVec = muonParticle->momentumVec();
    else muonParticle->setMomentumVec(muonfourVec);
    
    //-- calculate kinematics
    double p_muon   = sqrt( muonfourVec.E()*muonfourVec.E() - MinervaUnits::M_mu*MinervaUnits::M_mu );  
    double E_muon   = muonfourVec.E();
    double enu_muon = PhysicsCalculator->nuEnergyCCQE(E_muon,p_muon,tammy_muon_theta,charge,m_CCQEBindingEnergyMeV);
    double Q2_muon  = PhysicsCalculator->qSquaredCCQE(E_muon,p_muon,tammy_muon_theta,charge,m_CCQEBindingEnergyMeV);
    
    //-- stored in extra containers
    Gaudi::LorentzVector lepton;
    lepton.SetCoordinates(muonfourVec.px(),muonfourVec.py(), muonfourVec.pz(), muonfourVec.E());
    ccqeHyp->setLeptonEnergy(lepton);
    
    std::vector<double> tammy_muon_startPoint, tammy_muon_endPoint;
    tammy_muon_endPoint.push_back( (muonProng->minervaTracks().back())->lastState().x() );
    tammy_muon_endPoint.push_back( (muonProng->minervaTracks().back())->lastState().y() );
    tammy_muon_endPoint.push_back( (muonProng->minervaTracks().back())->lastState().z() );
    
    tammy_muon_startPoint.push_back( (muonProng->minervaTracks().front())->firstState().x() );
    tammy_muon_startPoint.push_back( (muonProng->minervaTracks().front())->firstState().y() );
    tammy_muon_startPoint.push_back( (muonProng->minervaTracks().front())->firstState().z() );
    
    //-- set muon data
    ccqeHyp->setDoubleData("tammy_muon_score",muonParticle->score());
    ccqeHyp->setDoubleData("tammy_muon_theta",tammy_muon_theta);
    ccqeHyp->setDoubleData("tammy_muon_thetaX",tammy_muon_thetaX);
    ccqeHyp->setDoubleData("tammy_muon_thetaY",tammy_muon_thetaY);
    ccqeHyp->setDoubleData("tammy_muon_phi",tammy_muon_phi);
    ccqeHyp->setDoubleData("tammy_muon_enu",enu_muon);
    ccqeHyp->setDoubleData("tammy_muon_q2",Q2_muon);
    ccqeHyp->setContainerDoubleData("tammy_muon_endPoint",tammy_muon_endPoint);
    ccqeHyp->setContainerDoubleData("tammy_muon_startPoint",tammy_muon_startPoint);
    
    //-- set proton data
    setProtonParticleData(ccqeHyp,protonProng,protonParticle,vertexZ);
    
    //-- rotated proton's angle w.r.t beam
    double tammy_proton_theta  = m_coordSysTool->thetaWRTBeam(protonParticle->momentumVec());
    double tammy_proton_phi    = m_coordSysTool->phiWRTBeam(protonParticle->momentumVec());
    
    //-- coplanarity angle
    Gaudi::XYZVector p_nu(0.0,0.0,1.0);
    Gaudi::XYZVector p_mu(sin(tammy_muon_theta)*cos(tammy_muon_phi),sin(tammy_muon_theta)*sin(tammy_muon_phi),cos(tammy_muon_theta));
    Gaudi::XYZVector p_proton(sin(tammy_proton_theta)*cos(tammy_proton_phi),sin(tammy_proton_theta)*sin(tammy_proton_phi),cos(tammy_proton_theta));
    
    Gaudi::XYZVector nvec = p_nu.Cross(p_mu);
    Gaudi::XYZVector mvec = p_nu.Cross(p_proton);
    
    double nlength = sqrt( nvec.x()*nvec.x() + nvec.y()*nvec.y() + nvec.z()*nvec.z() );
    double mlength = sqrt( mvec.x()*mvec.x() + mvec.y()*mvec.y() + mvec.z()*mvec.z() );
    double coplane = acos((nvec.Dot(mvec))/(nlength*mlength));
    ccqeHyp->setDoubleData("tammy_coplanarAngle",coplane);
    debug() << "   coplanar angle = " << coplane*180/3.14159 << endmsg;
    
    //-- set some pion information
    setPionParticleData(ccqeHyp,protonProng);
    debug() << "   set pion particle information" << endmsg;
    
    //-- set od muon particle data
    setODMuonParticleData(ccqeHyp,muonParticle);
    debug() << "   set the od muon particle data information" << endmsg;
    
    setMuonID(ccqeHyp,muonProng);
    debug() << "   set muon ID information" << endmsg;
    
    //-- fill some additional muon information
    //double tammy_proton_p = sqrt( protonfourVec.E()*protonfourVec.E() - MinervaUnits::M_proton*MinervaUnits::M_proton );
    //fillMuonMomentumFromProton(ccqeHyp,muonTracks,tammy_proton_theta,tammy_proton_p);
    //debug() << "   fill muon momentum from proton information" << endmsg;
    
    //-- fill common branches
    fillSystematicShiftsBranches(ccqeHyp,event);
    
    //-- store ccqe hypothesis
    nuInts.push_back( ccqeHyp );
    
    debug() << "Exit tammy_interpretEvent()" << endmsg;
    
    return StatusCode::SUCCESS;
    
}
//=============================================================================
// interpretEvent() --
//=============================================================================
StatusCode CCDeltaPlusAna::interpretEvent( const Minerva::PhysicsEvent *event, const Minerva::GenMinInteraction *truth, NeutrinoVect& nuInts ) const
{
    debug() << " == CCDeltaPlusAna::interpretEvent ==" << endmsg;
    if( event ) debug() << " working on Interpret Event and pass all cuts " << endmsg;
    if( truth ) debug() << "This event has a matched MC interaction<oaltinok_version>" << endmsg;
    
        // get charge of muon if available
    int muon_charge = 0;
    MuonUtils->muonCharge(m_MuonProng, muon_charge, m_qOverpChargeCut );
    
        // get origin of muon
    Gaudi::XYZTVector muon_position = m_MuonParticle->startPos();
    
        //----------------------------------------------------------------------------
        // Calculate stuff 
        //----------------------------------------------------------------------------
    double muon_P = m_MuonParticle->momentumVec().P();
    double muon_E = m_MuonParticle->momentumVec().E();
    double muon_T = muon_E - MinervaUnits::M_mu;
    
        // get muon momenum +/- shifts and re-calculate stuff 
    double muon_P_scale = 1.0;
    if( muon_P ) muon_P_scale = MuonUtils->calculateMomentumCorrection(m_MuonProng )/muon_P;


        // change units
    muon_P *= CLHEP::MeV/CLHEP::GeV;
    muon_E *= CLHEP::MeV/CLHEP::GeV;
    muon_T *= CLHEP::MeV/CLHEP::GeV;
    
        //! If you decide you want to interpret the event, create a new NeutrinoInt.
    Minerva::NeutrinoInt *nuInt = new Minerva::NeutrinoInt( "CCDeltaPlusAna" );
    
        //! Add the NeutrinoInt to the vector in return value
    nuInts.push_back( nuInt );
    
    nuInt->setNeutrinoFlavor( Minerva::NeutrinoInt::MuonFlavor );
    nuInt->setInteractionCurrent( Minerva::NeutrinoInt::ChargedCurrent );

    nuInt->setNeutrinoHelicity( PhysicsCalculator->getHelicity( muon_charge ) );
    nuInt->setLeptonEnergy( m_MuonParticle->momentumVec() );
    nuInt->setVertex( muon_position );
    
    nuInt->setDoubleData( "muon_P", muon_P );
    nuInt->setDoubleData( "muon_E", muon_E );
    nuInt->setDoubleData( "muon_T", muon_T );
    nuInt->setDoubleData( "muon_Pscale", muon_P_scale );


    return StatusCode::SUCCESS;
    
}
    
    
//=======================================================================
//  PreFilter
//=======================================================================
StatusCode CCDeltaPlusAna::PreFilter(Minerva::PhysicsEvent *event ) const
{
        //if ( event->processType() ==  Minerva::PhysicsEvent::RockParticle )  {
        //info() << " Jumping rock Muon " << endmsg;
        //return StatusCode::FAILURE;
        //}
    
    SmartRefVector<Minerva::IDCluster> anaClusters
        = event->select<Minerva::IDCluster>("Unused","!LowActivity&!XTalkCandidate");
    SmartRefVector< Minerva::IDCluster >::iterator itClus = anaClusters.begin();
    
    double ntgtEvis = 0.0;
    double trkrEvis = 0.0;
    double ecalEvis = 0.0;
    double hcalEvis = 0.0;
    double total = 0.0;
    for ( ; itClus != anaClusters.end(); itClus++ ){
        const double energy = (*itClus)->energy();
        Minerva::IDCluster::Subdet subdet = (*itClus)->subdet();
        total += energy;
        if      (subdet == Minerva::IDCluster::NuclTargs) ntgtEvis += energy;
        else if (subdet == Minerva::IDCluster::Tracker)   trkrEvis += energy;
        else if (subdet == Minerva::IDCluster::ECAL)      ecalEvis += energy;
        else if (subdet == Minerva::IDCluster::HCAL)      hcalEvis += energy;
        else {}
    }
    
    const double otherevis = trkrEvis + ecalEvis + hcalEvis;
    
    event->setDoubleData("totalevis", total);
    event->setDoubleData("ntgtevis",  ntgtEvis);
    event->setDoubleData("trkrevis",  trkrEvis);
    event->setDoubleData("ecalevis",  ecalEvis);
    event->setDoubleData("hcalevis",  hcalEvis);
    event->setDoubleData("otherevis", otherevis);
    
    if      ( ntgtEvis > 20)    return StatusCode::FAILURE;
    else if ( otherevis > 2000) return StatusCode::FAILURE; // energy bigger than 1.7*1.2 GeV must be ignored
    else if ( otherevis < 80 )  return StatusCode::FAILURE; // energy smaller than 80*1.2 Mev must be ignored
    else return StatusCode::SUCCESS;
    
    return StatusCode::SUCCESS;
    }
    
//=======================================================================
//  DoMuon
//=======================================================================
StatusCode CCDeltaPlusAna::DoMuon(Minerva::PhysicsEvent *event, Minerva::GenMinInteraction* truth,
                                SmartRef<Minerva::Vertex>& vertex ) const
{
    // - Muon exist
    bool has_muon = MuonUtils->findMuonProng( event, m_MuonProng, m_MuonParticle );
    
    if( truth ) truth->filtertaglist()->setOrAddFilterTag( "reco_minos_match", has_muon );
    if( !has_muon ){
        debug() << "Did not find a muon prong! This cannot be a CCDeltaPlus event." << endmsg;
        return StatusCode::FAILURE; // We didn't crash.
    }

    event->setIntData("survive_minos_match", 1);

    if (has_muon) {

        int is_contained = -1;    // 1)
        int quality      = -1;    // 2)
        int is_ok        = -1;    // 3)
        double fit_pass  = -1;    // 4)
        double qp        = 0.0;   // 5)
        double eqp       = 0.0;   // 6)
        double p         = -1.0;  // 7)
        double p_range   = -1.0;  // 8)
        double p_curve   = -1.0;  // 9)
        int used_range   = -1;    // 10)
        int used_curve   = -1;    // 11)
        int end_plane    = -1;    // 12)
        
        m_MuonProng->getIntData(    ProngExtraDataDefs::MinosTrackContained(),  is_contained);
        m_MuonProng->getIntData(    ProngExtraDataDefs::MinosTrackQuality(),    quality);
        m_MuonProng->getIntData(    ProngExtraDataDefs::MinosOK(),              is_ok);
        m_MuonProng->getDoubleData( ProngExtraDataDefs::MinosTrackFitPass(),    fit_pass);
        m_MuonProng->getDoubleData( ProngExtraDataDefs::MinosTrackQP(),         qp);
        m_MuonProng->getDoubleData( ProngExtraDataDefs::MinosTrackSigmaQP(),    eqp);
        m_MuonProng->getDoubleData( ProngExtraDataDefs::MinosMomentum(),        p);
        m_MuonProng->getDoubleData( ProngExtraDataDefs::MinosTrackPRange(),     p_range);
        m_MuonProng->getIntData(    ProngExtraDataDefs::MinosRecoByRange(),     used_range );
        m_MuonProng->getIntData(    ProngExtraDataDefs::MinosRecoByCurvature(), used_curve);
        m_MuonProng->getIntData(    ProngExtraDataDefs::MinosTrackEndPlane(),   end_plane);

        p_curve = qp != 0.0 ? std::abs(1./qp) : -1;

        event->setIntData(   "minos_trk_is_contained",  is_contained);
        event->setIntData(   "minos_trk_quality",       quality);
        event->setIntData(   "minos_trk_is_ok",         is_ok);
        event->setDoubleData("minos_trk_fit_pass",      fit_pass);
        event->setDoubleData("minos_trk_qp",            qp);
        event->setDoubleData("minos_trk_eqp",           eqp);
        event->setDoubleData("minos_trk_p",             p);
        event->setDoubleData("minos_trk_prange",        p_range);
        event->setDoubleData("minos_trk_pcurvature",    p_curve);
        event->setIntData(   "minos_trk_used_range",    used_range);
        event->setIntData(   "minos_trk_used_curvature",used_curve);
        event->setIntData(   "minos_trk_end_plane",     end_plane);
        
    }
    
    if( !event->hasInteractionVertex() ) {
        debug() << " The interaction vertex for this event is not set - this tool cannot do reconstruction without one."
                << endmsg;
    return StatusCode::FAILURE; // We didn't crash.
    }

    event->setIntData("survive_has_vertex",1);
    
    vertex = event->interactionVertex();
    
    Gaudi::XYZPoint position = vertex->position();
    
    if( !FiducialPointTool->isFiducial( position, m_fiducialApothem, m_fiducialUpstreamZ, m_fiducialDownstreamZ ) ){
        debug() << " Interaction Vertex is not fiducial!" << endmsg;
        return StatusCode::FAILURE; // We didn't crash.
    }

    event->setIntData("survive_fiducial",1);

    bool plausible = shouldAnalyzeMC(truth);
    event->setIntData("survive_plausible", plausible);
    
    bool transverse_fiducial = InsideHexagon(position.X(),position.Y(),2*850.0);
    event->setIntData("tfiducial", transverse_fiducial);
    
// make it primary for fill fillCommonPhysicsAnaBranches
    m_MuonProng->filtertaglist()->setOrAddFilterTag( AnaFilterTags::PrimaryMuon(), true ); 
    
        //  -- Muon Blob
    Minerva::IDBlob *muonBlob = new Minerva::IDBlob;	
//SmartRef<Minerva::IDBlob> muonBlob;	
    double muon_blob_energy = 0;

    SmartRefVector<Minerva::IDCluster> extraMuonClusters
        = m_extraEnergyTool->getExtraIDClustersCylinder( m_MuonProng,
                                                        m_extraEnergyCylinderUpstreamLength,
                                                        m_extraEnergyCylinderDownstreamLength,
                                                        m_extraEnergyCylinderRadius,
                                                        IExtraEnergyTool::k_all,
                                                        IExtraEnergyTool::k_inside,
                                                        m_extraEnergyLowerTimeWindow,
                                                        m_extraEnergyUpperTimeWindow,
                                                        m_extraEnergyPECut,
                                                        IExtraEnergyTool::k_clusters_unused,
                                                        -100 ); /* -100 means clusters in all modules! */

    debug() << " Extra muon contain " << extraMuonClusters.size() << " clusters<oaltinok_version>" << endmsg;
    std::vector<BlobSeeds::seedCandVector> VectorSeeds;
    m_blobSeedingTool->makeViewSortedBlobSeeds( extraMuonClusters, VectorSeeds );
    std::vector<BlobSeeds::seedCandVector>::iterator itSeeds = VectorSeeds.begin();
    for ( ; itSeeds != VectorSeeds.end(); itSeeds++ ){
        BlobSeeds::seedCandVector::iterator itSeed = (*itSeeds).begin();
        for ( ; itSeed != (*itSeeds).end(); itSeed++ ){
            debug() << " Seed with Energy " << (*itSeed)->seedenergy << " " << ((*itSeed)->IdClustVec).size()
                    <<  endmsg;
            if ( (*itSeed)->seedenergy < 30 )	{ 
                muonBlob->add( (*itSeed)->IdClustVec );
                debug() << " Blob size " << muonBlob->nclusters() << endmsg; 
            }
        }
    }

    if ( muonBlob->nclusters() ) {
        muonBlob->setHistory( Minerva::IDBlob::Unused );
        muonBlob->setPatRecHistory( Minerva::IDBlob::IsolatedIDBlobPatRec );
        m_idHoughBlob->getBlobEnergyTime( muonBlob, muon_blob_energy);
        debug()<< "Adding Muon blob with " << muonBlob->nclusters() << " clusters; energy = "  << muon_blob_energy
            << endmsg;
        addObject( event, muonBlob );
    }
    
    event->setDoubleData( "Muon_blob_energy", muon_blob_energy );

    DigitVectorTruthInfo info;
    info.ParseTruth(muonBlob->getAllDigits(),fTrajectoryMap);
    double evis = info.GetEdepByPdg(111);
    event->setDoubleData("pi0_evis_muon_blob", evis);
    
    
    return StatusCode::SUCCESS;

}
    
//=======================================================================
//  DoVertex
//=======================================================================
StatusCode CCDeltaPlusAna::DoVertex(Minerva::PhysicsEvent *event,
                                const SmartRef<Minerva::Vertex>& vertex) const
{
    debug() << " CCDeltaPlusAna::DoVertex " << endmsg;

    const Gaudi::XYZPoint& pos = vertex->position();
    
    SmartRefVector<Minerva::Vertex> allVertices  = event->select<Minerva::Vertex>( "All","StartPoint" );

    event->setIntData("vertex_count",allVertices.size());

    SmartRefVector<Minerva::Vertex>::iterator primaryVertexIterator = allVertices.end();
    for (SmartRefVector<Minerva::Vertex>::iterator vtx = allVertices.begin() ;
        vtx != allVertices.end(); ++vtx ){

        if (*vtx == vertex) {
            primaryVertexIterator = vtx;
            break;
        }
    }

    event->setIntData("primary_index", std::distance(allVertices.begin(),primaryVertexIterator));
    
    if (primaryVertexIterator != allVertices.end()) allVertices.erase(primaryVertexIterator);

    event->setIntData("vertex_count2", allVertices.size()); /* Number of secondary vertices */

    int primaryVertexTrackCount = vertex->getNTracks();
    event->setIntData("primary_multiplicity", primaryVertexTrackCount);

    SmartRef<Minerva::Track> muonTrack = m_MuonProng->minervaTracks().front();
    assert(muonTrack);
    m_MuonTrack = muonTrack;
    
    const int maxTrack = 5;
    SmartRefVector<Minerva::Track> tracks = vertex->getTracks();
    std::vector<double> primary_separations(maxTrack,1e6);
    std::vector<double> primary_trklengths(maxTrack,1e6);

    std::vector<int>    primary_truth_counts(maxTrack,-1);
    std::vector<double> primary_truth_shareds(maxTrack,0);
    std::vector<int>    primary_truth_pdgs1(maxTrack,-1);
    std::vector<int>    primary_truth_pdgs2(maxTrack,-1);
    std::vector<int>    primary_truth_pdgs3(maxTrack,-1);
    std::vector<double> primary_truth_fractions1(maxTrack,0);
    std::vector<double> primary_truth_fractions2(maxTrack,0);
    std::vector<double> primary_truth_fractions3(maxTrack,0);
    std::size_t index = 0;
    std::size_t trackMultiplicity = 0;
    SmartRefVector<Minerva::Track> farTracks;
    for (SmartRefVector<Minerva::Track>::iterator t = tracks.begin();
        t != tracks.end(); ++t) {
        if (muonTrack && (*t == muonTrack)) {
            std::cout << "\t\tSkipping the muon track" << std::endl;
            continue;
        }

            /* Not sure which end of the track is closest to the primary vertex.
            Take the shorter of the distances to both track ends */
        Gaudi::XYZPoint firstPos = (*t)->firstNode()->position();
        Gaudi::XYZPoint lastPos  = (*t)->lastNode()->position();

        const double distance = std::min((firstPos-vertex->position()).R(),
                                        (lastPos-vertex->position()).R());
        std::cout << "\t\tDistance to primary vertex: " << distance << std::endl;

            /* Count tracks close to the primary vertex and use it
            as an updated track multiplicity */
        if (distance < 50.0) ++trackMultiplicity;
        if (distance > 50.0) farTracks.push_back(*t);

        const double trackLength = CalcTrackLength(*t);
        std::cout << "\t\tTrack length: " << trackLength << std::endl;

        primary_separations[index] = distance;
        primary_trklengths[index]  = trackLength;

        TrackTruthInfo truthInfo;
        truthInfo.ParseTruth(*t,fTrajectoryMap);

        primary_truth_counts[index]  = truthInfo.GetTruthParticleCount();
        primary_truth_shareds[index] = truthInfo.GetSharedFraction();
        
        int pdg1 = truthInfo.GetTruthPdg(0);

        std::cout << "\t\tTruth particle: " << pdg1
                << " " << truthInfo.GetEdepFraction()
                << " " << truthInfo.GetSharedFraction()
                << std::endl;

        primary_truth_pdgs1[index]      = pdg1;
        primary_truth_fractions1[index] = truthInfo.GetEdepFraction();
        
        int pdg2 = truthInfo.GetTruthPdg(1);
        std::cout << "\t\tTruth particle: " << pdg2
                << " " << truthInfo.GetEdepFraction()
                << " " << truthInfo.GetSharedFraction()
                << std::endl;
        primary_truth_pdgs2[index]      = pdg2;
        primary_truth_fractions2[index] = truthInfo.GetEdepFraction();
        
        int pdg3 = truthInfo.GetTruthPdg(2);
        std::cout << "\t\tTruth particle: " << pdg3
                << " " << truthInfo.GetEdepFraction()
                << " " << truthInfo.GetSharedFraction()
                << std::endl;

        primary_truth_pdgs3[index]      = pdg3;
        primary_truth_fractions3[index] = truthInfo.GetEdepFraction();

        ++index;

        if (index > maxTrack-1) break;
    }

    event->setIntData("primary_multiplicity2", trackMultiplicity);
    
    event->setContainerDoubleData("primary_separations",primary_separations);
    event->setContainerDoubleData("primary_trklengths",primary_trklengths);
    event->setContainerIntData("primary_truth_counts",   primary_truth_counts);
    event->setContainerDoubleData("primary_truth_shareds",primary_truth_shareds);
    event->setContainerIntData("primary_truth_pdgs1",   primary_truth_pdgs1);
    event->setContainerIntData("primary_truth_pdgs2",   primary_truth_pdgs2);
    event->setContainerIntData("primary_truth_pdgs3",   primary_truth_pdgs3);
    event->setContainerDoubleData("primary_truth_fractions1",primary_truth_fractions1);
    event->setContainerDoubleData("primary_truth_fractions2",primary_truth_fractions2);
    event->setContainerDoubleData("primary_truth_fractions3",primary_truth_fractions3);
    

    int multiTrackVertexCount = 0;
    SmartRefVector<Minerva::Track> tracksToClear;

    std::vector<int>    multiplicities;
    std::vector<double> deviations;
    for (SmartRefVector<Minerva::Vertex>::iterator vtx = allVertices.begin() ;
        vtx != allVertices.end(); ++vtx ){
        
        SmartRefVector<Minerva::Track> tracks = (*vtx)->getTracks();
        if ( tracks.size() != 1 ) ++multiTrackVertexCount;

            /* For unclear reason, Jose projects and discards later only
            the first track (outgoing from secondary vertices) */
            
        Gaudi::XYZPoint projectedVertex;
        SmartRef<Minerva::Track> firstTrack = tracks.front();
        m_trackPropagator->propagate(firstTrack, pos.z(), projectedVertex);
        double distance = (pos-projectedVertex).R();
        if ( distance < 50 ) tracksToClear.push_back(firstTrack);

        multiplicities.push_back(tracks.size());
        deviations.push_back(distance);
    }

        /* Put them in the same container before deleting */
    std::copy(tracksToClear.begin(),tracksToClear.end(),
            std::back_inserter(farTracks));
    
    if (!farTracks.empty()) {
        debug() << " Discarding " <<  farTracks.size() << " tracks, to get photons " << endmsg;
        std::size_t oldsize = farTracks.size();
        Minerva::EventMgr *mgr = getEventMgr(allVertices.front());
        discardObject(mgr,farTracks);
        std::size_t newsize = farTracks.size();
        assert(newsize == oldsize);
        
    }
    event->setIntData("discard_track_count", tracksToClear.size());
    event->setContainerIntData("multiplicities", multiplicities);
    event->setContainerDoubleData("deviations", deviations);

    
    return StatusCode::SUCCESS;

}

//=======================================================================
//  VtxBlob
//=======================================================================
StatusCode CCDeltaPlusAna::VtxBlob(Minerva::PhysicsEvent *event, const SmartRef<Minerva::Vertex>& vertex ) const
{

debug() << " CCDeltaPlusAna::VtxBlob " << endmsg;

    //  --  Vertex blob
double vertex_energy = 0;
double vertex_energy_sphere = 0;
double vertex_energy_filament= 0;

const Gaudi::XYZPoint& pos = vertex->position();

if ( m_filamentVertex ){
    
    SmartRefVector<Minerva::IDCluster> analyFilaClusters
        = event->select<Minerva::IDCluster>("Unused","!LowActivity&!XTalkCandidate");

    SmartRefVector<Minerva::IDCluster> filaClusters;
    m_blobUtils->fillProximateClusterVec(pos,analyFilaClusters,filaClusters,
                                        m_maxSearchDFila,m_maxStartingDFila,m_maxSearchGapFila);

    if (!filaClusters.empty()) {
        
        Minerva::IDBlob* vtxFilaBlob = new Minerva::IDBlob();
        m_blobUtils->insertIDClusters( filaClusters, vtxFilaBlob, Minerva::IDBlob::VertexBlobPatRec );
        m_idHoughBlob->getBlobEnergyTime(vtxFilaBlob,vertex_energy_filament);
        debug() << " Adding Filament vertex " << vtxFilaBlob->nclusters()
                << " clusters; energy = "  << vertex_energy_filament << endmsg;
        addObject( event, vtxFilaBlob );
        m_colorTag->applyColorTag( vtxFilaBlob, 0x9900FF );
    }

}

if ( m_sphereVertex ){

    SmartRefVector<Minerva::IDCluster> unusedClusters
        = event->select<Minerva::IDCluster>("Unused","!LowActivity&!XTalkCandidate");
    
    SmartRefVector<Minerva::IDCluster> sphereClusters = FilterInSphereClusters(vertex,unusedClusters,m_maxSearchD);
    
    if (!sphereClusters.empty()) {

        Minerva::IDBlob* vtxSphereBlob = new Minerva::IDBlob();
        m_blobUtils->insertIDClusters( sphereClusters, vtxSphereBlob, Minerva::IDBlob::VertexBlobPatRec );
        m_idHoughBlob->getBlobEnergyTime(vtxSphereBlob,vertex_energy_sphere);
        debug() << " Adding Sphere vertex " << vtxSphereBlob->nclusters()
                << " clusters; energy = "  <<  vertex_energy_sphere << endmsg;
        addObject( event, vtxSphereBlob );
        m_colorTag->applyColorTag( vtxSphereBlob, 0x9900FF );

        DigitVectorTruthInfo info;
        info.ParseTruth(vtxSphereBlob->getAllDigits(),fTrajectoryMap);
        double evis = info.GetEdepByPdg(111);
        event->setDoubleData("pi0_evis_vtx_blob", evis);
    }
    
}	

vertex_energy = vertex_energy_sphere + vertex_energy_filament;
event->setDoubleData( "Vertex_blob_energy", vertex_energy );
event->setDoubleData( "Filament_Vertex_energy", vertex_energy_filament );
event->setDoubleData( "Sphere_Vertex_energy", vertex_energy_sphere );

return StatusCode::SUCCESS;

}

//=======================================================================
//  ConeBlobs
//=======================================================================
StatusCode CCDeltaPlusAna::ConeBlobs(Minerva::PhysicsEvent *event,
                                const SmartRef<Minerva::Vertex>& vertex) const
{

    debug() << " CCDeltaPlusAna::ConeBlobs " << endmsg;

        // variables to clean up the clusters
    SmartRefVector<Minerva::IDCluster> preidClusters
        = event->select<Minerva::IDCluster>("Unused","!LowActivity&!XTalkCandidate");
    SmartRefVector<Minerva::IDCluster> usableClusters;
    SmartRefVector<Minerva::IDCluster> outTimeClusters;
    
        // get origin of muon
    Gaudi::XYZTVector muon_position = m_MuonParticle->startPos();
    
    double energyTracker = 0, energyECAL = 0, energyHCAL = 0;
    bool isAngleScan = false;
    bool isAngleScanApplied = false;
    bool isHough = false;
    bool isHoughApplied = false;
    
    SmartRefVector<Minerva::IDCluster>::iterator it_clus;
    for ( it_clus = preidClusters.begin(); it_clus != preidClusters.end(); it_clus++){
        
        if ( (*it_clus)->pe()/(*it_clus)->iddigs() <= 3 ) continue;
        
        if ( (*it_clus)->subdet() ==  Minerva::IDCluster::Tracker ) energyTracker += (*it_clus)->energy();
        if ( (*it_clus)->subdet() ==  Minerva::IDCluster::ECAL )    energyECAL += (*it_clus)->energy();
        if ( (*it_clus)->subdet() ==  Minerva::IDCluster::HCAL )    energyHCAL += (*it_clus)->energy();
        
            // must include clusters close to MUON vertex time < 25 ns
        if ( std::abs( (*it_clus)->time() -  muon_position.T() ) < m_rejectedClustersTime ) {
            usableClusters.push_back(*it_clus);
        } else {
            outTimeClusters.push_back(*it_clus);
        }
    }

    event->setDoubleData( "RE_energy_Tracker", energyTracker );
    event->setDoubleData( "RE_energy_ECAL", energyECAL );
    event->setDoubleData( "RE_energy_HCAL", energyHCAL );

    info() << " Energy to analyze " << energyTracker + energyECAL + energyHCAL << endmsg;

    std::vector<Minerva::IDBlob*> foundBlobs;
    unsigned int nblob_anglescan = 0;
    unsigned int nblob_hough     = 0;
    if ( ( energyTracker + energyECAL + energyHCAL ) < m_energyHoughlimit ) {

        if (new_impl_) {

            AngleScan angleScanAlg(usableClusters,vertex->position());
            angleScanAlg.AllowUVMatchWithMoreTolerance(fAllowUVMatchWithMoreTolerance);
            angleScanAlg.SetUVMatchTolerance(fUVMatchTolerance);
            angleScanAlg.DoReco();

            event->setIntData("anglescan_ncandx", angleScanAlg.GetNxCandidate());
            event->setIntData("anglescan_ncand",  angleScanAlg.GetNCandidate());

            std::vector<int> candx_nc;
            std::vector<int> candx_nd;
            const std::vector<SmartRefVector<Minerva::IDCluster> >& xshowerCand = angleScanAlg.GetXShowerCandVector();
            for (std::vector<SmartRefVector<Minerva::IDCluster> >::const_iterator s
                    = xshowerCand.begin();
                s != xshowerCand.end(); ++s) {
                ClusterVectorInfo clusterInfo(*s);
                candx_nc.push_back(clusterInfo.GetN());
                candx_nd.push_back(clusterInfo.GetDigitCount());
            }
            
            event->setContainerIntData("anglescan_candx_nc", candx_nc);
            event->setContainerIntData("anglescan_candx_nd", candx_nd);
                        
            
            std::vector<int> cand_nc;
            std::vector<int> cand_ncx;
            std::vector<int> cand_ncu;
            std::vector<int> cand_ncv;

            std::vector<int> cand_nd;
            std::vector<int> cand_ndx;
            std::vector<int> cand_ndu;
            std::vector<int> cand_ndv;
            const std::vector<SmartRefVector<Minerva::IDCluster> >& showerCand = angleScanAlg.GetShowerCandVector();
            for (std::vector<SmartRefVector<Minerva::IDCluster> >::const_iterator s
                    = showerCand.begin();
                s != showerCand.end(); ++s) {
                ClusterVectorInfo clusterInfo(*s);
                cand_nc.push_back(clusterInfo.GetN());
                cand_ncx.push_back(clusterInfo.GetNx());
                cand_ncu.push_back(clusterInfo.GetNu());
                cand_ncv.push_back(clusterInfo.GetNv());

                cand_nd.push_back(clusterInfo.GetDigitCount());
                cand_ndx.push_back(clusterInfo.GetDigitCountX());
                cand_ndu.push_back(clusterInfo.GetDigitCountU());
                cand_ndv.push_back(clusterInfo.GetDigitCountV());
            }

            event->setContainerIntData("anglescan_cand_nc",   cand_nc);
            event->setContainerIntData("anglescan_cand_ncx",  cand_ncx);
            event->setContainerIntData("anglescan_cand_ncu",  cand_ncu);
            event->setContainerIntData("anglescan_cand_ncv",  cand_ncv);
            
            event->setContainerIntData("anglescan_cand_nd",   cand_nd);
            event->setContainerIntData("anglescan_cand_ndx",  cand_ndx);
            event->setContainerIntData("anglescan_cand_ndu",  cand_ndu);
            event->setContainerIntData("anglescan_cand_ndv",  cand_ndv);

            std::vector<Minerva::IDBlob*> angleScanBlobs = angleScanAlg.GetShowers();
            for (std::vector<Minerva::IDBlob*>::const_iterator b = angleScanBlobs.begin();
                b != angleScanBlobs.end(); ++b) {
                if ( !m_idHoughBlob->isPhoton((*b)->clusters(),vertex->position())) continue;
                foundBlobs.push_back(*b);
            }

            std::vector<int> blob_nc;
            std::vector<int> blob_ncx;
            std::vector<int> blob_ncu;
            std::vector<int> blob_ncv;

            std::vector<int> blob_nd;
            std::vector<int> blob_ndx;
            std::vector<int> blob_ndu;
            std::vector<int> blob_ndv;        
            for (std::vector<Minerva::IDBlob*>::const_iterator b = foundBlobs.begin();
                b != foundBlobs.end(); ++b) {
                ClusterVectorInfo clusterInfo((*b)->clusters());
                blob_nc.push_back(clusterInfo.GetN());
                blob_ncx.push_back(clusterInfo.GetNx());
                blob_ncu.push_back(clusterInfo.GetNu());
                blob_ncv.push_back(clusterInfo.GetNv());

                blob_nd.push_back(clusterInfo.GetDigitCount());
                blob_ndx.push_back(clusterInfo.GetDigitCountX());
                blob_ndu.push_back(clusterInfo.GetDigitCountU());
                blob_ndv.push_back(clusterInfo.GetDigitCountV());
            }

            event->setContainerIntData("anglescan_blob_nc", blob_nc);
            event->setContainerIntData("anglescan_blob_ncx",blob_ncx);
            event->setContainerIntData("anglescan_blob_ncu",blob_ncu);
            event->setContainerIntData("anglescan_blob_ncv",blob_ncv);

            event->setContainerIntData("anglescan_blob_nd", blob_nd);
            event->setContainerIntData("anglescan_blob_ndx",blob_ndx);
            event->setContainerIntData("anglescan_blob_ndu",blob_ndu);
            event->setContainerIntData("anglescan_blob_ndv",blob_ndv);
            
        } // end if new_impl
        else {
            
            StatusCode sc = AngleScanBlob(usableClusters, vertex, foundBlobs );
        
        }

        nblob_anglescan = foundBlobs.size();
        isAngleScan = foundBlobs.size() == 2;
        isAngleScanApplied = true;
        
    }
    if ( !isAngleScan) { 

            /* Blobs found by the AngleScan are not managed, delete them here and
            empty the container */
        for (std::vector<Minerva::IDBlob*>::iterator b = foundBlobs.begin();
            b != foundBlobs.end(); ++b) {
            delete *b;
        }
        foundBlobs.clear();
        
        StatusCode sc = HoughBlob(usableClusters, vertex, foundBlobs );
        nblob_hough = foundBlobs.size();
        
        if ( sc && foundBlobs.size() == 2) isHough = true;
        isHoughApplied = true;

        std::vector<int> blob_nc;
        std::vector<int> blob_ncx;
        std::vector<int> blob_ncu;
        std::vector<int> blob_ncv;

        std::vector<int> blob_nd;
        std::vector<int> blob_ndx;
        std::vector<int> blob_ndu;
        std::vector<int> blob_ndv;
        for (std::vector<Minerva::IDBlob*>::const_iterator b = foundBlobs.begin();
            b != foundBlobs.end(); ++b) {
            ClusterVectorInfo clusterInfo((*b)->clusters());
            blob_nc.push_back(clusterInfo.GetN());
            blob_ncx.push_back(clusterInfo.GetNx());
            blob_ncu.push_back(clusterInfo.GetNu());
            blob_ncv.push_back(clusterInfo.GetNv());

            blob_nd.push_back(clusterInfo.GetDigitCount());
            blob_ndx.push_back(clusterInfo.GetDigitCountX());
            blob_ndu.push_back(clusterInfo.GetDigitCountU());
            blob_ndv.push_back(clusterInfo.GetDigitCountV());
        }

        event->setContainerIntData("hough_blob_nc", blob_nc);
        event->setContainerIntData("hough_blob_ncx",blob_ncx);
        event->setContainerIntData("hough_blob_ncu",blob_ncu);
        event->setContainerIntData("hough_blob_ncv",blob_ncv);

        event->setContainerIntData("hough_blob_nd", blob_nd);
        event->setContainerIntData("hough_blob_ndx",blob_ndx);
        event->setContainerIntData("hough_blob_ndu",blob_ndu);
        event->setContainerIntData("hough_blob_ndv",blob_ndv);
        
    }
        /* if either isAngleScan or isHough is true, the blobs are finally stored
        in the TES, which are managed. Otherwise, delete them */
    double minBlobSep1 = -1.;
    double minBlobSep2 = -1.;
    if (isAngleScan || isHough) {
        processBlobs(event,foundBlobs);
        fillPi0Branches(event,foundBlobs[0],foundBlobs[1],vertex);

        minBlobSep1 = CalcMinBlobSeparation(foundBlobs[0],vertex);
        minBlobSep2 = CalcMinBlobSeparation(foundBlobs[1],vertex);

        if (foundBlobs[0]->energy() < foundBlobs[1]->energy()) std::swap(minBlobSep1,minBlobSep2);
        
        std::cout << "\tmin_sep: " << minBlobSep1 << " " << minBlobSep2 << std::endl;

        std::vector<int> blob_nc;
        std::vector<int> blob_ncx;
        std::vector<int> blob_ncu;
        std::vector<int> blob_ncv;

        std::vector<int> blob_nd;
        std::vector<int> blob_ndx;
        std::vector<int> blob_ndu;
        std::vector<int> blob_ndv;
        for (std::vector<Minerva::IDBlob*>::const_iterator b = foundBlobs.begin();
            b != foundBlobs.end(); ++b) {
            ClusterVectorInfo clusterInfo((*b)->clusters());
            blob_nc.push_back(clusterInfo.GetN());
            blob_ncx.push_back(clusterInfo.GetNx());
            blob_ncu.push_back(clusterInfo.GetNu());
            blob_ncv.push_back(clusterInfo.GetNv());
            
            blob_nd.push_back(clusterInfo.GetDigitCount());
            blob_ndx.push_back(clusterInfo.GetDigitCountX());
            blob_ndu.push_back(clusterInfo.GetDigitCountU());
            blob_ndv.push_back(clusterInfo.GetDigitCountV());
        }

        event->setContainerIntData("final_blob_nc", blob_nc);
        event->setContainerIntData("final_blob_ncx",blob_ncx);
        event->setContainerIntData("final_blob_ncu",blob_ncu);
        event->setContainerIntData("final_blob_ncv",blob_ncv);

        event->setContainerIntData("final_blob_nd", blob_nd);
        event->setContainerIntData("final_blob_ndx",blob_ndx);
        event->setContainerIntData("final_blob_ndu",blob_ndu);
        event->setContainerIntData("final_blob_ndv",blob_ndv);

        
        ODActivity(event,foundBlobs);
            
    } else {
        for (std::vector<Minerva::IDBlob*>::iterator b = foundBlobs.begin();
            b != foundBlobs.end(); ++b) { 
            delete *b;
        }
        foundBlobs.clear();
    }
    
    event->filtertaglist()->setOrAddFilterTag( "is_anglescan", isAngleScan );
    event->filtertaglist()->setOrAddFilterTag( "is_anglescan_applied",isAngleScanApplied);
    event->filtertaglist()->setOrAddFilterTag( "is_houghtransform", isHough );
    event->filtertaglist()->setOrAddFilterTag( "is_houghtransform_applied",isHoughApplied);
    
    event->setIntData("nblob_anglescan", nblob_anglescan);
    event->setIntData("nblob_hough",     nblob_hough);
    event->setDoubleData("g1blob_minsep", minBlobSep1);
    event->setDoubleData("g2blob_minsep", minBlobSep2);
    
        //storing rejected id Clusters
    Minerva::IDBlob* rejectedBlob = new Minerva::IDBlob();
    
    if (!outTimeClusters.empty()) {
        m_blobUtils->insertIDClusters(outTimeClusters, rejectedBlob, Minerva::IDBlob::DispersedIDBlobPatRec );
        debug()<< "Adding rejected blob with vis energy = "  << rejectedBlob->energy() << endmsg;
        addObject( event, rejectedBlob );
        m_colorTag->applyColorTag( rejectedBlob, 0x000000 );
    }

    event->setDoubleData( "Rejected_blob_vis_energy", rejectedBlob->energy() );

    DigitVectorTruthInfo info;
    info.ParseTruth(rejectedBlob->getAllDigits(),fTrajectoryMap);
    double evis = info.GetEdepByPdg(111);
    event->setDoubleData("pi0_evis_outtime_blob", evis);
    
    return StatusCode::SUCCESS;

}

//=======================================================================
//  AngleScan
//=======================================================================
StatusCode CCDeltaPlusAna::AngleScanBlob(SmartRefVector<Minerva::IDCluster> idClusters,
                                    const SmartRef<Minerva::Vertex>& vertex,
                                    std::vector<Minerva::IDBlob*>& outBlobs) const
{
    debug() << " CCDeltaPlusAna::AngleScanBlob " << endmsg;

    const Gaudi::XYZPoint& pos = vertex->position();
    
        // main interface to create blob, this tool already exclude lowactivity clusters.
    std::vector<Minerva::IDBlob*>* preidBlobs = new std::vector<Minerva::IDBlob*>;
    m_idConeScanBlob->createIDBlobs(idClusters,preidBlobs,pos);
    
    std::vector<Minerva::IDBlob*>::iterator itBlob = preidBlobs->begin();
    for ( ; itBlob != preidBlobs->end(); itBlob ++ ){
        if ( !m_idHoughBlob->isPhoton( (*itBlob)->clusters(), pos ) ) continue;
        outBlobs.push_back(*itBlob);
    }

    info() << " Angle Scan is done " << endmsg;
    
    return StatusCode::SUCCESS;
}

//=======================================================================
//  HoughBlob
//=======================================================================
StatusCode CCDeltaPlusAna::HoughBlob(SmartRefVector<Minerva::IDCluster> idClusters,
                                const SmartRef<Minerva::Vertex>& vertex,
                                std::vector<Minerva::IDBlob*>& outBlobs) const
{
    debug() << " CCDeltaPlusAna::HoughBlob " << endmsg;

    const Gaudi::XYZPoint& pos = vertex->position();
    
    SmartRefVector<Minerva::IDCluster> idClusViewX, idClusViewU, idClusViewV;
    SmartRefVector<Minerva::IDCluster>::iterator it_clus;

    if ( !m_idHoughBlob->GetViewClustersAboveThreshold( idClusters, idClusViewX, Minerva::IDCluster::X, 0 )) return StatusCode::FAILURE;
    m_idHoughBlob->GetViewClustersAboveThreshold( idClusters, idClusViewU, Minerva::IDCluster::U, 0 );
    m_idHoughBlob->GetViewClustersAboveThreshold( idClusters, idClusViewV, Minerva::IDCluster::V, 0 );
    
    double r, theta, spX, spZ;
    std::vector<Minerva::IDBlob*> idBlobs;
    
    for ( int i = 0; i < 3; i++){

        Minerva::IDBlob *idBlob = new Minerva::IDBlob;
        SmartRefVector<Minerva::IDCluster> idHTSeed;
        Gaudi::XYZPoint ref(0,0,0);// energetic cluster position
        
        if ( !m_idHoughTool->GetReference( idClusViewX, ref ) ) break;
        debug() << " Pass Get Reference<oaltinok_version>" << endmsg;	

        if ( !m_idHoughTool->Hough2D( idClusViewX, r, theta, ref, pos ) ) break;
        debug() << " Pass Get Hough2d<oaltinok_version>" << endmsg;	
        
        if ( !m_idHoughBlob->Create2dHTSeed( idClusViewX, idHTSeed, r, theta, ref, spX, spZ ) ) break;
        debug() << " Pass Get Create2dHT " << endmsg;	
        
        Gaudi::XYZPoint startpoint(spX, 0, spZ);
        Gaudi::XYZVector direction(-1/tan(theta*CLHEP::pi/180),0,1);

        if ( !m_idHoughBlob->PseudoCone(idHTSeed, idClusViewX, direction, pos ) ) continue;
        debug() << " Pass PseudoCone " << endmsg;
        
        if ( !m_idHoughBlob->XUVMatch(idHTSeed, idClusViewU, idClusViewV, 50 ) ) continue;
        debug() << " Pass XUV Match - First " << endmsg;	
        
        if ( !m_idHoughBlob->isPhoton(idHTSeed, pos) ) continue;
        
        idBlob->add( idHTSeed );
        debug() << " Added idBlob with size " << idHTSeed.size() << endmsg;

        idBlobs.push_back(idBlob);
        
    }

    if ( idBlobs.empty() ) return StatusCode::FAILURE;
    
    std::vector<Minerva::IDBlob*>::iterator itBlob;
    for (itBlob = idBlobs.begin(); itBlob != idBlobs.end(); ++itBlob ){
        Minerva::IDBlob *blobTemp = new Minerva::IDBlob;
        SmartRefVector<Minerva::IDCluster> idSeed = (*itBlob)->clusters();
        if ( !m_idHoughBlob->XUVMatch( idSeed, idClusViewU, idClusViewV, 120 ) ) continue;
        blobTemp->add(idSeed);
        double Uclusters = blobTemp->nclusters(Minerva::IDCluster::U);
        double Vclusters = blobTemp->nclusters(Minerva::IDCluster::V);
        debug() << *itBlob << " " << Uclusters << " U clusters; " << Vclusters << " V clusters " << endmsg;
        (*itBlob)->clear();
        if ( Uclusters != .0 || Vclusters != .0 ) {
            debug() << *itBlob << " inserting to FinalBlobs " << endmsg;
            (*itBlob)->add(idSeed);
            m_idHoughBlob->GetStartPosition( *itBlob, pos, true );
            m_idHoughBlob->GetDirection( *itBlob, pos );
            outBlobs.push_back(*itBlob);
        }
    }

    for ( it_clus = idClusViewX.begin(); it_clus != idClusViewX.end(); it_clus++ )
            m_idHoughBlob->AddClusterInsideCone( *it_clus, outBlobs, pos );

    for ( it_clus = idClusViewU.begin(); it_clus != idClusViewU.end(); it_clus++ )
            m_idHoughBlob->AddClusterInsideCone( *it_clus, outBlobs, pos );

    for ( it_clus = idClusViewV.begin(); it_clus != idClusViewV.end(); it_clus++ )
            m_idHoughBlob->AddClusterInsideCone( *it_clus, outBlobs, pos );


    info() << " Hough Transform is done! " << endmsg; 

    
    return StatusCode::SUCCESS;
}

//=======================================================================
//  processBlobs
//=======================================================================
StatusCode CCDeltaPlusAna::processBlobs( Minerva::PhysicsEvent *event, std::vector<Minerva::IDBlob*> idBlobs) const
{

    int count = 0;
    std::vector<Minerva::IDBlob*>::iterator itBlob = idBlobs.begin();    
    for (; itBlob != idBlobs.end(); ++itBlob){
        (*itBlob)->setHistory( Minerva::IDBlob::Unused );
        (*itBlob)->setPatRecHistory( Minerva::IDBlob::Pi0IDBlobPatRec );
        if ( (*itBlob)->nclusters( Minerva::IDCluster::X ) > 0 )  count++;
        if ( (*itBlob)->nclusters( Minerva::IDCluster::U ) > 0 )  count++;
        if ( (*itBlob)->nclusters( Minerva::IDCluster::V ) > 0 )  count++;
        debug() << " Storing Photon blob with " << (*itBlob)->nclusters() << " clusters "
                << endmsg;
        addObject( event, *itBlob );
    }


    if ( count != 6 ) event->filtertaglist()->setOrAddFilterTag( "is_twoDBlob", true );
    m_colorTag->applyColorTag( (idBlobs)[0], 0xFF0000 );
    m_colorTag->applyColorTag( (idBlobs)[1], 0x0000FF );

    debug() << "pi0 candidate<oaltinok_version>" << endmsg;

    std::pair<int,double> result1 = OneParLineFitBlob(idBlobs[0],m_PrimaryVertex,m_MuonTrack);
    std::pair<int,double> result2 = OneParLineFitBlob(idBlobs[1],m_PrimaryVertex,m_MuonTrack);

    event->setIntData("blob_ndof_1",result1.first);
    event->setIntData("blob_ndof_2",result2.first);
    event->setDoubleData("blob_fval_1", result1.second);
    event->setDoubleData("blob_fval_2", result2.second);
    
    return StatusCode::SUCCESS;
}

//=======================================================================
//  DispersedBlob
//=======================================================================
StatusCode CCDeltaPlusAna::DispersedBlob( Minerva::PhysicsEvent *event ) const
{

//  --  Dispersed blob
    double dispersed_energy = 0;
    
    Minerva::IDBlob* dispersedBlob = new Minerva::IDBlob();
    SmartRefVector<Minerva::IDCluster> preidClusters
        = event->select<Minerva::IDCluster>("Unused","!LowActivity&!XTalkCandidate");
    SmartRefVector<Minerva::IDCluster> dispersedClusters;
    SmartRefVector<Minerva::IDCluster>::iterator it_clus;

    for ( it_clus = preidClusters.begin(); it_clus != preidClusters.end(); it_clus++){
        if ( (*it_clus)->pe()/(*it_clus)->iddigs() <= 3 ) continue; 
        dispersedClusters.push_back(*it_clus);
    }
    
    if (!dispersedClusters.empty()) {
        m_blobUtils->insertIDClusters( dispersedClusters, dispersedBlob, Minerva::IDBlob::DispersedIDBlobPatRec );
        m_idHoughBlob->getBlobEnergyTime(dispersedBlob,dispersed_energy);
        debug()<< "Adding dispersed blob with " << dispersedBlob->nclusters()
            << " clusters; Dispersed energy = "  << dispersed_energy << endmsg;
        addObject( event, dispersedBlob );
    }

    event->setDoubleData( "Dispersed_blob_energy", dispersed_energy );

    info() << "Dispersed blob energy: " << dispersed_energy << endmsg;
    PrintDigitVectorInfo(dispersedBlob->getAllDigits());
    
    DigitVectorTruthInfo info;
    info.ParseTruth(dispersedBlob->getAllDigits(),fTrajectoryMap);
    double evis = info.GetEdepByPdg(111);
    event->setDoubleData("pi0_evis_dispersed_blob",evis);
    
    return StatusCode::SUCCESS;

}

//=======================================================================
//  invariantMass -  will calculate direction and startpoint
//=======================================================================
StatusCode CCDeltaPlusAna::fillPi0Branches(Minerva::PhysicsEvent *event, Minerva::IDBlob* idblob1, Minerva::IDBlob* idblob2,
                                        const SmartRef<Minerva::Vertex>& vertex ) const
{
    debug() << "CCDeltaPlusAna::fillPi0Branches "  << endmsg;

    const Gaudi::XYZPoint& pos = vertex->position();
            
    bool goodPosition1  = m_idHoughBlob->GetStartPosition(idblob1, pos, true );
    bool goodDirection1 = m_idHoughBlob->GetDirection(idblob1, pos );
    bool isGoodBlob1 = false;
    if (goodPosition1 && goodDirection1) isGoodBlob1 = true;
    
    bool goodPosition2  = m_idHoughBlob->GetStartPosition(idblob2, pos, true );
    bool goodDirection2 = m_idHoughBlob->GetDirection(idblob2, pos );
    bool isGoodBlob2 = false;
    if (goodPosition2 && goodDirection2) isGoodBlob2 = true;
    if (!isGoodBlob1 || !isGoodBlob2) return StatusCode::SUCCESS;
    double energy1 = 0.0;
    double energy2 = 0.0;
    m_idHoughBlob->getBlobEnergyTime( idblob1, energy1 );
    m_idHoughBlob->getBlobEnergyTime( idblob2, energy2 );

        /* Make sure gamma1 is the more energetic one */
    if (energy2 > energy1) {
        std::swap(goodPosition1, goodPosition2);  /* Swap variables already assigned */
        std::swap(goodDirection1,goodDirection2);
        std::swap(isGoodBlob1,   isGoodBlob2);
        std::swap(energy1,energy2);
        std::swap(idblob1,idblob2);               /* And the blobs themselves */
    }
    
    double time1 = idblob1->time();
    double time2 = idblob2->time();

    double dEdx1 = 0.0;
    double dEdx2 = 0.0;
    m_idHoughBlob->idBlobdEdx( idblob1, dEdx1 );
    m_idHoughBlob->idBlobdEdx( idblob2, dEdx2 );

    Gaudi::XYZVector direction1 = idblob1->direction(); 
    Gaudi::XYZVector direction2	= idblob2->direction();


    std::vector<double> position1;
    std::vector<double> position2;
    position1.push_back(idblob1->startPoint().x());
    position1.push_back(idblob1->startPoint().y());
    position1.push_back(idblob1->startPoint().z());
    position2.push_back(idblob2->startPoint().x());
    position2.push_back(idblob2->startPoint().y());
    position2.push_back(idblob2->startPoint().z());

    
    double cos_oangle  = direction1.Dot(direction2);
    
    std::vector<double> direc_1, direc_2;
    direc_1.push_back(direction1.x());
    direc_1.push_back(direction1.y());
    direc_1.push_back(direction1.z());
    
    direc_2.push_back(direction2.x());
    direc_2.push_back(direction2.y());
    direc_2.push_back(direction2.z());
    std::cout << "Direction: " << direction1 << " " << direction2
            << cos_oangle << std::endl;
    event->setDoubleData( "RE_scalar", cos_oangle );
    event->setDoubleData( "RE_photon_energy_1", energy1 );
    event->setDoubleData( "RE_photon_energy_2", energy2 );
    event->setDoubleData( "RE_photon_dEdx_1", dEdx1 );
    event->setDoubleData( "RE_photon_dEdx_2", dEdx2 );
    event->setDoubleData( "RE_photon_time_1", time1 );
    event->setDoubleData( "RE_photon_time_2", time2 );
    event->filtertaglist()->setOrAddFilterTag("is_GoodDirection1", goodDirection1);
    event->filtertaglist()->setOrAddFilterTag("is_GoodPosition1",  goodPosition1);
    event->filtertaglist()->setOrAddFilterTag("is_GoodDirection2", goodDirection2);
    event->filtertaglist()->setOrAddFilterTag("is_GoodPosition2",  goodPosition2);
    event->filtertaglist()->setOrAddFilterTag("is_GoodBlob1", isGoodBlob1);
    event->filtertaglist()->setOrAddFilterTag("is_GoodBlob2", isGoodBlob2);
    event->setContainerDoubleData( "RE_photon_direction_1", direc_1 );
    event->setContainerDoubleData( "RE_photon_direction_2", direc_2 );
    event->setContainerDoubleData( "RE_photon_vertex_1", position1 );
    event->setContainerDoubleData( "RE_photon_vertex_2", position2 );

    event->setIntData("g1blob_ndigit",  idblob1->getAllDigits().size());
    event->setIntData("g2blob_ndigit",  idblob2->getAllDigits().size());
    event->setIntData("g1blob_ncluster",idblob1->nclusters());
    event->setIntData("g2blob_ncluster",idblob2->nclusters());

    
    TVector3 g1mom(direction1.x(),direction1.y(),direction1.z());
    TVector3 g2mom(direction2.x(),direction2.y(),direction2.z());
    g1mom *= energy1;
    g2mom *= energy2;

    TVector3 pimom = g1mom + g2mom;
    
    const double oangle = g1mom.Angle(g2mom)*TMath::RadToDeg();
    event->setDoubleData("oangle",oangle);
    const double mgg = std::sqrt(2*energy1*energy2*(1-cos_oangle));
    event->setDoubleData("mgg", mgg);

    
    event->setDoubleData("pienergy", energy1+energy2);
    event->setDoubleData("pitheta",  pimom.Theta()*TMath::RadToDeg());
    event->setDoubleData("piphi",    pimom.Phi()*TMath::RadToDeg());

    event->setDoubleData("g1e",    energy1);
    event->setDoubleData("g1theta",g1mom.Theta()*TMath::RadToDeg());
    event->setDoubleData("g1phi",  g1mom.Phi()*TMath::RadToDeg());
    event->setDoubleData("g2e",    energy2);
    event->setDoubleData("g2theta",g2mom.Theta()*TMath::RadToDeg());
    event->setDoubleData("g2phi",  g2mom.Phi()*TMath::RadToDeg());

    std::vector<double> pimomVec;
    std::vector<double> g1momVec;
    std::vector<double> g2momVec;

    pimomVec.push_back(pimom.x());
    pimomVec.push_back(pimom.y());
    pimomVec.push_back(pimom.z());
    pimomVec.push_back(pimom.Mag());

    g1momVec.push_back(g1mom.x());
    g1momVec.push_back(g1mom.y());
    g1momVec.push_back(g1mom.z());
    g2momVec.push_back(g2mom.x());
    g2momVec.push_back(g2mom.y());
    g2momVec.push_back(g2mom.z());
    
    event->setContainerDoubleData("pimom",pimomVec);
    event->setContainerDoubleData("g1mom",g1momVec);
    event->setContainerDoubleData("g2mom",g2momVec);


    const double dvtx1 = CalcDistanceFromBlobAxisToVertex(idblob1,vertex);
    const double dvtx2 = CalcDistanceFromBlobAxisToVertex(idblob2,vertex);
    event->setDoubleData("g1blob_vtx_distance",dvtx1);
    event->setDoubleData("g2blob_vtx_distance",dvtx2);

    const double dmax1 = CalcDistanceFromVertexToExiting(idblob1,vertex);
    const double dmax2 = CalcDistanceFromVertexToExiting(idblob2,vertex);
    event->setDoubleData("g1blob_edge_distance",dmax1);
    event->setDoubleData("g2blob_edge_distance",dmax2);
    
    
    CheckBlobDigitAssignment(event,idblob1,1);
    CheckBlobDigitAssignment(event,idblob2,2);

        // Convenient variables for analysis
        // 1) Muon and pi0 momentum in the detector coordinates
    const double pmu_x = m_MuonParticle->momentumVec().Px();
    const double pmu_y = m_MuonParticle->momentumVec().Py();
    const double pmu_z = m_MuonParticle->momentumVec().Pz();
    const double Emu    = m_MuonParticle->momentumVec().E();
    
    const double ppi_x = pimom.x();
    const double ppi_y = pimom.y();
    const double ppi_z = pimom.z();
    const double Epi    = energy1 + energy2;
    
        // 2) transform into beam coordinates
    const double pmu_y2 =  costheta_b*pmu_y + sintheta_b*pmu_z;
    const double pmu_z2 = -sintheta_b*pmu_y + costheta_b*pmu_z;
    
    const double ppi_y2 =  costheta_b*ppi_y + sintheta_b*ppi_z;
    const double ppi_z2 = -sintheta_b*ppi_y + costheta_b*ppi_z;
    
        // 3) Calculate the neutrino energy
    info() << "Ozgur>> Using new Neutrino Energy Estimation"<< endmsg;
    const double pT   = pmu_x + pmu_y2 + ppi_x + ppi_y2;
    const double pTx = pmu_x + ppi_x;
    const double pTy = pmu_y2 + ppi_y2;
    const double pmuL = Emu - pmu_z2;
    const double ppiL = Epi - ppi_z2;
    const double pL   = pmuL + ppiL;
    const double Tn   = 0.5*(pL*pL + pT*pT)/(mn - pL);
    const double Tn2   = 0.5*(pL*pL + pTx*pTx + pTy*pTy)/(mn - pL);
    const double Erec = Emu + Epi + Tn;
    const double Erec2 = Emu + Epi + Tn2;

    const double Q2 = 2*Erec*(Emu-pmu_z2) - mmuon*mmuon;
    const double W2 = mp*mp + 2*mp*(Erec-Emu) - Q2;
    const double W  = std::sqrt(std::max(0.,W2));

    std::vector<double> mumom;
    mumom.push_back(pmu_x);
    mumom.push_back(pmu_y);
    mumom.push_back(pmu_z);
    mumom.push_back(m_MuonParticle->momentumVec().R());
    mumom.push_back(Emu);
    
    event->setContainerDoubleData("mumom", mumom);
    event->setDoubleData("Tn", Tn);
    event->setDoubleData("Tn2", Tn2);
    event->setDoubleData("Erec", Erec);
    event->setDoubleData("Erec2", Erec2);
    event->setDoubleData("Q2", Q2);
    event->setDoubleData("W2", W2);
    event->setDoubleData("W", W);
        
    return StatusCode::SUCCESS;

}
    
//=============================================================================
// tagTruth 
//=============================================================================
StatusCode CCDeltaPlusAna::tagTruth( Minerva::GenMinInteraction* truth ) const 
{
    debug() <<"Enter CCDeltaPlusAna::tagTruth()" << endmsg;
    debug() <<"=============================================================================" <<endmsg;
    debug() <<" "<<endmsg;
    
    
    if( !truth ) {
        warning() << "The GenMinInteraction is NULL!" << endmsg;
        return StatusCode::SUCCESS;
    }
    
//     fillGenieWeightBranches(truth);
//     
//     Minerva::PhysicsEvents* events = getPhysicsEvents();
//     for(Minerva::PhysicsEvents::iterator event = events->begin(); event != events->end(); event++) {
//         if( !(*event)->hasInteractionVertex() ) continue;
//         SmartRef<Minerva::Vertex> vertex = (*event)->interactionVertex();
//         if( !vertex ) continue;     
//         std::string name = "";
//         int tammy_passVertexZCut = 0;
//         if( !isInsideTargetFiducial(vertex,name,tammy_passVertexZCut) ) continue;     
//         if( tammy_passVertexZCut == 0 ) continue;
//         if( vertex->getOutgoingTracks().size() > 2 ) continue;
//         int ntracks = 0;
//         Minerva::ProngVect unattachedProngs = (*event)->select<Minerva::Prong>("Used:Unused","!IsSubProng");    
//         for(unsigned int i = 0; i < unattachedProngs.size(); i++) if( !unattachedProngs[i]->minervaTracks().empty() ) ntracks++; 
//         if( ntracks != 0 ) continue;
//         
//         if( truth->processType() == 1 && truth->current() == 1 && truth->incoming() == 14 ) {
//         int nprotons = 0;
//         const std::vector<int> fsParticlePDG = truth->fSpdg();
//         for(unsigned int i = 0; i < fsParticlePDG.size(); i++) { if( fsParticlePDG[i] == 2212 ) nprotons++; }
//     
//         if( nprotons != 0 ) truth->filtertaglist()->addFilterTag( "CCDeltaPlusAnaSignal", true );
//         }  
//     }
//     
//     //-- tag truth michel electrons
//     verbose() << "  Tagging truth michel electrons." << endmsg;
//     
//     Gaudi::LorentzVector positron_pos;
//     Gaudi::LorentzVector electron_pos;
//     
//     std::map<int,double> tammy_pion_plus_ids_momentum, tammy_pion_minus_ids_momentum;
//     std::map<int,int> tammy_muon_plus_id_parent, tammy_muon_minus_id_parent;
//     std::map<int,Gaudi::LorentzVector> electron_parent_pos, positron_parent_pos;
//     std::vector<double> tammy_pion_plus_momentum, tammy_pion_minus_momentum;
//     
//     Minerva::TG4Trajectories* alltrajects = get<Minerva::TG4Trajectories>( "MC/TG4Trajectories" );
//     for(Minerva::TG4Trajectories::iterator it = alltrajects->begin(); it != alltrajects->end(); ++it) {
//         Minerva::TG4Trajectory* traject = *it;
//     
//         //! look for charged pions
//         if( traject->GetPDGCode()==211 ){
//             tammy_pion_plus_ids_momentum[traject->GetTrackId()] = (traject->GetInitialMomentum()).P();
//         }
//         if( traject->GetPDGCode()==-211 ){
//             tammy_pion_minus_ids_momentum[traject->GetTrackId()] = (traject->GetInitialMomentum()).P();
//         }
//     
//         //! look for muons from decays
//         if( traject->GetProcessName() == "Decay" && traject->GetPDGCode() == -13 ) {
//             tammy_muon_plus_id_parent[traject->GetTrackId()] = traject->GetParentId();
//         }
//     
//         if( traject->GetProcessName() == "Decay" && traject->GetPDGCode() == 13 ) {
//             tammy_muon_minus_id_parent[traject->GetTrackId()] = traject->GetParentId();
//         }
//     
//         //! look for electrons from decays
//         if( traject->GetProcessName() == "Decay" && traject->GetPDGCode() == -11 ) {
//             positron_parent_pos[traject->GetParentId()] = traject->GetInitialPosition();
//         }
//     
//         if( traject->GetProcessName() == "muMinusCaptureAtRest" && traject->GetPDGCode() == 11 ) {
//             electron_parent_pos[traject->GetParentId()] = traject->GetInitialPosition();
//         }
//     
//     } // end of TG4Trajectories loop
//     
//     //! loop over electrons
//     for( std::map<int,Gaudi::LorentzVector>::iterator electron = electron_parent_pos.begin(); electron != electron_parent_pos.end(); ++electron ){
//         // check if electron parent is a muon-
//         if( tammy_muon_minus_id_parent.count( electron->first ) ){
//         //check if electron is inside fiducial region
//         Gaudi::LorentzVector electron_pos = electron->second;
//         Gaudi::XYZPoint p( electron_pos.x(), electron_pos.y(), electron_pos.z() );
//         if( FiducialPointTool->isFiducial( p, m_analyzeApothem, m_michel_upstreamZ, m_michel_downstreamZ ) ){
//             truth->setIntData("tammy_has_michel_electron",1);
//             // check if muon parent is a pion-
//             if( tammy_pion_minus_ids_momentum.count( tammy_muon_minus_id_parent[electron->first] ) ){
//                 tammy_pion_minus_momentum.push_back( tammy_pion_minus_ids_momentum[ tammy_muon_minus_id_parent[electron->first] ]  );
//             }
//         }
//         }
//     }
//     
//     //! loop over positrons
//     for( std::map<int,Gaudi::LorentzVector>::iterator positron = positron_parent_pos.begin(); positron != positron_parent_pos.end(); ++positron ){
//         // check if positron parent is a muon+
//         if( tammy_muon_plus_id_parent.count( positron->first ) ){
//             //check if electron is inside fiducial region
//             Gaudi::LorentzVector positron_pos = positron->second;
//             Gaudi::XYZPoint p( positron_pos.x(), positron_pos.y(), positron_pos.z() );
//             if( FiducialPointTool->isFiducial( p, m_analyzeApothem, m_michel_upstreamZ, m_michel_downstreamZ ) ){
//                 truth->setIntData("tammy_has_michel_electron",1);
//                 // check if muon parent is a pion+ 
//                 if( tammy_pion_plus_ids_momentum.count( tammy_muon_plus_id_parent[positron->first] ) ){
//                     tammy_pion_plus_momentum.push_back( tammy_pion_plus_ids_momentum[ tammy_muon_plus_id_parent[positron->first] ]  );
//                 }
//             }
//         }
//     }
//     
//     truth->setContainerDoubleData("tammy_has_michel_from_tammy_pion_plus_momentum", tammy_pion_plus_momentum);
//     truth->setContainerDoubleData("tammy_has_michel_from_tammy_pion_minus_momentum", tammy_pion_minus_momentum);
//     verbose() << "   End tagging truth michel electrons." << endmsg;
//     verbose() << "Exit NukeCCQE::tagTruth" << endmsg;
//         
//         
//     //! Check to see if the MC event is plausibly analyzable - important for Data overlay.
//     bool isPlausible = truthIsPlausible( truth );
//     truth->filtertaglist()->setOrAddFilterTag( "pass_plausible", isPlausible );
//         
//     // is interaction in fiducial volume?
//     bool is_fiducial = FiducialPointTool->isFiducial( truth, m_fiducialApothem, m_fiducialUpstreamZ, m_fiducialDownstreamZ );
//     truth->filtertaglist()->setOrAddFilterTag( "is_fiducial", is_fiducial );
//     
//     debug() << " Pass fiducial point " << is_fiducial<< endmsg;
//     
//     // get primary final state lepton stuff
//     double fslepton_P = truth->PrimFSLepton().P(); 
//     double fslepton_E = truth->PrimFSLepton().E();
//     double fslepton_T = fslepton_E - MinervaUnits::M_mu;
//     double fslepton_theta = m_minervaCoordSysTool->thetaWRTBeam( truth->PrimFSLepton() );
//     double fslepton_theta_x = m_minervaCoordSysTool->thetaXWRTBeam( truth->PrimFSLepton() );
//     double fslepton_theta_y = m_minervaCoordSysTool->thetaYWRTBeam( truth->PrimFSLepton() );
//     double fslepton_phi = m_minervaCoordSysTool->phiWRTBeam( truth->PrimFSLepton() );
//     
//     truth->setDoubleData( "fslepton_E", fslepton_E*CLHEP::MeV/CLHEP::GeV );
//     truth->setDoubleData( "fslepton_P", fslepton_P*CLHEP::MeV/CLHEP::GeV );
//     truth->setDoubleData( "fslepton_T", fslepton_T*CLHEP::MeV/CLHEP::GeV );
//     truth->setDoubleData( "fslepton_theta", fslepton_theta );
//     truth->setDoubleData( "fslepton_theta_x", fslepton_theta_x );
//     truth->setDoubleData( "fslepton_theta_y", fslepton_theta_y );
//     truth->setDoubleData( "fslepton_phi", fslepton_phi );
//     
//     debug() << " Filled truth " << endmsg;
//         // clasify by ccpi0 ccpi0x other
//     bool is_ccpi0          = false;
//     bool is_cc1pi0         = false;
//     bool is_ccpi0secondary = false;
//     bool is_by_pim         = false;
//     bool is_ccpi0x         = false;
//     bool is_other          = false;
//     
//     int charge = 0;
//     if(truth->primaryLepton()==13 ) charge = -1;
//     else if(truth->primaryLepton()==-13 ) charge = 1;
//     
//     if ( is1pi0(truth)) {
//         is_cc1pi0 = true;
//         getMCPi0Info(truth);
//     }
//     
//     if ( ispi0( truth, charge ) )	{  is_ccpi0 = true; }
//     else if ( ispi0x( truth ) ) is_ccpi0x = true;
//     else is_other = true;
//     
//     std::pair<bool,bool> secondary_pi0_info = ispi0secondary();
//     if (secondary_pi0_info.first)  is_ccpi0secondary = true;
//     if (secondary_pi0_info.second) is_by_pim = true;
//     
//     truth->filtertaglist()->setOrAddFilterTag( "is_ccpi0",          is_ccpi0 );
//     truth->filtertaglist()->setOrAddFilterTag( "is_cc1pi0",         is_cc1pi0 );
//     truth->filtertaglist()->setOrAddFilterTag( "is_ccpi0secondary", is_ccpi0secondary );
//     truth->filtertaglist()->setOrAddFilterTag( "is_by_pim",         is_by_pim);
//     truth->filtertaglist()->setOrAddFilterTag( "is_ccpi0x",         is_ccpi0x );
//     truth->filtertaglist()->setOrAddFilterTag( "is_other",          is_other );
    
    
    debug() <<"Do Nothing!" << endmsg;
    debug() <<"Exit CCDeltaPlusAna::tagTruth()" << endmsg;
    
    return StatusCode::SUCCESS;
    
}
    
//=============================================================================
// ispi0
//=============================================================================
bool CCDeltaPlusAna::ispi0( Minerva::GenMinInteraction* truth, int charge) const
{
    if ( truth->nParticlesFS() == 3 && charge == 1 && truth->fSpdg()[1] == 2112 && truth->fSpdg()[2] == 111 ) {
        return true;
    }
    else if ( truth->nParticlesFS() == 3 && charge == -1 && truth->fSpdg()[1] == 2212 && truth->fSpdg()[2] == 111 ) {
        return true;
    }
    else return false;
}
    
bool CCDeltaPlusAna::is1pi0(const Minerva::GenMinInteraction* v) const
{
    int pizeroCount = 0;
    int pionCount = 0;
    int othermesonCount = 0;

    const std::vector<int>& fsPdgs = v->fSpdg();
    for (std::vector<int>::const_iterator i = fsPdgs.begin();
        i != fsPdgs.end(); ++i) {
        const int pdg = *i;
        const int abspdg = std::abs(pdg);
        
        if (pdg == 111) ++pizeroCount;

        if (abspdg == 211) ++pionCount;
        
        const bool notpion = (pdg != 111) && (abspdg != 211);
        if (notpion) {
            TParticlePDG* particle = TDatabasePDG::Instance()->GetParticle(pdg);
                /* particle could be a nucleus fragment which makes 'particle' invalid,
                do not assert(particle) */
            if (particle && std::string(particle->ParticleClass()) == "Meson") ++othermesonCount;
            
        }
        
    }

    return pizeroCount == 1 && pionCount == 0 && othermesonCount == 0;
}

//=============================================================================
// ispi0secondary
//=============================================================================
std::pair<bool,bool> CCDeltaPlusAna::ispi0secondary( ) const
{
        /* Is there a secondary pi0, is it by pi- charge exchange? */
    bool secondary_pi0 = false;
    bool by_pim = false;
    for (std::map<int, Minerva::TG4Trajectory*>::const_iterator traj = fTrajectoryMap.begin();
        traj != fTrajectoryMap.end(); ++traj) {
        const Minerva::TG4Trajectory* trajectory = traj->second;

        if (trajectory->GetParentId() == 0) continue; /* Skip primary particle */

        if (trajectory->GetPDGCode() == 111) {
            secondary_pi0 = true;
            std::map<int,Minerva::TG4Trajectory*>::const_iterator target
                = fTrajectoryMap.find(trajectory->GetParentId());
            if (target != fTrajectoryMap.end() && target->second->GetPDGCode() == -211) by_pim = true;
        
            break;
        }
        
    }

    return std::make_pair(secondary_pi0,by_pim);
        
}
    
//=============================================================================
// ispi0x
//=============================================================================
bool CCDeltaPlusAna::ispi0x( Minerva::GenMinInteraction* truth) const
{
// - Pi0 form primary particles?
    for( unsigned int i = 0; i < truth->nParticlesFS(); i++ ){
        if( truth->fSpdg()[i] == 111 && truth->fsParticlesE()[i] > 50 ) return true;
    }
    
    return false;
}

//=======================================================================
//  getMCPi0info
//=======================================================================
StatusCode CCDeltaPlusAna::getMCPi0Info(Minerva::GenMinInteraction* truth) const
{


    if( !exist<Minerva::TG4Trajectories>( "MC/TG4Trajectories" ) )
    {
        warning() << "No TG4Trajectories found in this gate. Skipping to next one." << endmsg;
        return StatusCode::SUCCESS;
    }

    
    Minerva::TG4Trajectories* trajectories = get<Minerva::TG4Trajectories>(Minerva::TG4TrajectoryLocation::Default);
    Minerva::TG4Trajectories::iterator traj_it = trajectories->begin();
    
    int trkid = 0;
    Gaudi::LorentzVector a[2], mc_pi0, b[2], vtx[2], vtemp[2];

    Minerva::TG4Trajectory* pizero(NULL);
    std::vector<Minerva::TG4Trajectory*> pizeroDaughters;
    for( ; traj_it != trajectories->end(); traj_it++ ) {
        if( (*traj_it)->GetPDGCode()==111 && (*traj_it)->GetParentId()== 0 ) {
            trkid = (*traj_it)->GetTrackId();
            mc_pi0 = (*traj_it)->GetInitialMomentum();
            pizero = *traj_it;
            break;
        }
    }

    int count = 0;	
    for ( traj_it = trajectories->begin(); traj_it != trajectories->end(); traj_it++ ) {
        if ( (*traj_it)->GetPDGCode()==22 && (*traj_it)->GetParentId()== trkid ) {
            a[count]   = (*traj_it)->GetInitialMomentum();
            vtx[count] = (*traj_it)->GetFinalPosition();
            count++;
            pizeroDaughters.push_back(*traj_it);
        }
        if ( count == 2 ) break;
    } 

//need order photon1 most energetic
    if (a[0].e() < a[1].e()){
        b[0] = a[1]; b[1] = a[0];
        a[0] = b[0]; a[1] = b[1];
        
        vtemp[0] = vtx[1]; vtemp[1] = vtx[0];
        vtx[0] = vtemp[0]; vtx[1] = vtemp[1];
    }

        /* This assert statement was triggered, probably because of a pi0
        outside tracking volume. Return StatusCode::SUCCESS if empty
        for now */
        //assert(!pizeroDaughters.empty());

    if (pizeroDaughters.empty()) return StatusCode::SUCCESS;
    
    Minerva::TG4Trajectory* g1(NULL);
    Minerva::TG4Trajectory* g2(NULL);
    if (pizeroDaughters.front()->GetInitialMomentum().E() >=
        pizeroDaughters.back()->GetInitialMomentum().E()) {
        g1 = pizeroDaughters.front();
        g2 = pizeroDaughters.back();
    } else {
        g1 = pizeroDaughters.back();
        g2 = pizeroDaughters.front();
    }
    
    std::vector<double> direc_1;
    direc_1.push_back(a[0].px());
    direc_1.push_back(a[0].py());
    direc_1.push_back(a[0].pz());
    std::vector<double> direc_2;
    direc_2.push_back(a[1].px());
    direc_2.push_back(a[1].py());
    direc_2.push_back(a[1].pz());

    std::vector<double> position_1;
    position_1.push_back(vtx[0].px());
    position_1.push_back(vtx[0].py());
    position_1.push_back(vtx[0].pz());
    std::vector<double> position_2;
    position_2.push_back(vtx[1].px());
    position_2.push_back(vtx[1].py());
    position_2.push_back(vtx[1].pz());
    
    std::vector<double> direc_pi0;
    direc_pi0.push_back(mc_pi0.px());
    direc_pi0.push_back(mc_pi0.py());
    direc_pi0.push_back(mc_pi0.pz());
    
    double scal = a[0].px()*a[1].px() + a[0].py()*a[1].py() + a[0].pz()*a[1].pz(); 
    scal = scal/sqrt(pow(a[0].px(),2)+pow(a[0].py(),2)+pow(a[0].pz(),2));
    scal = scal/sqrt(pow(a[1].px(),2)+pow(a[1].py(),2)+pow(a[1].pz(),2));

    truth->setDoubleData( "MC_pi0_energy",mc_pi0.e() );
    truth->setDoubleData( "MC_photon_energy_1",a[0].e() );
    truth->setDoubleData( "MC_photon_energy_2",a[1].e() );
    truth->setDoubleData( "MC_scalar", scal );
    
    truth->setContainerDoubleData( "MC_pi0_momentum", direc_pi0 );
    truth->setContainerDoubleData( "MC_photon_direction_1",  direc_1 );
    truth->setContainerDoubleData( "MC_photon_direction_2",  direc_2 );
    truth->setContainerDoubleData( "MC_photon_vertex_1",   position_1 );
    truth->setContainerDoubleData( "MC_photon_vertex_2",   position_2 );


    truth->setDoubleData("pienergy0", pizero->GetInitialMomentum().E());
    truth->setDoubleData("pitheta0",  pizero->GetInitialMomentum().Theta()*TMath::RadToDeg());
    truth->setDoubleData("piphi0",    pizero->GetInitialMomentum().Phi()*TMath::RadToDeg());

    truth->setDoubleData("g1e0",    g1->GetInitialMomentum().E());
    truth->setDoubleData("g1theta0",g1->GetInitialMomentum().Theta()*TMath::RadToDeg());
    truth->setDoubleData("g1phi0",  g1->GetInitialMomentum().Phi()*TMath::RadToDeg());
    truth->setDoubleData("g2e0",    g2->GetInitialMomentum().E());
    truth->setDoubleData("g2theta0",g2->GetInitialMomentum().Theta()*TMath::RadToDeg());
    truth->setDoubleData("g2phi0",  g2->GetInitialMomentum().Phi()*TMath::RadToDeg());

    truth->setDoubleData("g1conv_distance",(g1->GetInitialPosition().Vect()-truth->Vtx().Vect()).R());
    truth->setDoubleData("g2conv_distance",(g2->GetInitialPosition().Vect()-truth->Vtx().Vect()).R());

    std::vector<double> pimomVec;
    std::vector<double> g1momVec;
    std::vector<double> g2momVec;
    pimomVec.push_back(pizero->GetInitialMomentum().Px());
    pimomVec.push_back(pizero->GetInitialMomentum().Py());
    pimomVec.push_back(pizero->GetInitialMomentum().Pz());
    pimomVec.push_back(pizero->GetInitialMomentum().P());
    
    g1momVec.push_back(g1->GetInitialMomentum().Px());
    g1momVec.push_back(g1->GetInitialMomentum().Py());
    g1momVec.push_back(g1->GetInitialMomentum().Pz());
    
    g2momVec.push_back(g2->GetInitialMomentum().Px());
    g2momVec.push_back(g2->GetInitialMomentum().Py());
    g2momVec.push_back(g2->GetInitialMomentum().Pz());

    truth->setContainerDoubleData("pimom0",pimomVec);
    truth->setContainerDoubleData("g1mom0",g1momVec);
    truth->setContainerDoubleData("g2mom0",g2momVec);

    TVector3 g1mom(g1->GetInitialMomentum().Px(),
                g1->GetInitialMomentum().Py(),
                g1->GetInitialMomentum().Pz());
    TVector3 g2mom(g2->GetInitialMomentum().Px(),
                g2->GetInitialMomentum().Py(),
                g2->GetInitialMomentum().Pz());
    
    const double oangle0 = g1mom.Angle(g2mom)*TMath::RadToDeg();
    truth->setDoubleData("oangle0",oangle0);
    
    return StatusCode::SUCCESS;
    
}
    

//=======================================================================
//  shouldAnalyzeMC - Gabe P.
//=======================================================================
bool CCDeltaPlusAna::shouldAnalyzeMC( const Minerva::GenMinInteraction *truth ) const
{
    debug() << " CCDeltaPlusAna::shouldAnalyzeMC " << endmsg;

    if( getDAQHeader()->isMCTrigger() && NULL == truth ) return false;
    if( truth ) {
        if( !muonIsPlausible ( m_MuonProng ) ){
            debug() << "This MC muon in not plausible because it has a MC energy fraction below threshold."
                    << endmsg; 
            
            return false;
        }
    }
    
    return true;

}

void CCDeltaPlusAna::CheckBlobDigitAssignment(Minerva::PhysicsEvent* event, const Minerva::IDBlob* blob, int index) const
{
    std::cout << "Blob energy: " << blob->energy() << std::endl;
    
    FillTrajectoryMap();
    
    SmartRefVector<Minerva::IDDigit> digits = blob->getAllDigits();

    DigitVectorTruthInfo digitVectorInfo;
    digitVectorInfo.ParseTruth(digits,fTrajectoryMap);

    const double pi0evis  = digitVectorInfo.GetEdepByPdg(111);
    double g1evis = 0.0;
    double g2evis = 0.0;
    if (fGamma1) g1evis   = digitVectorInfo.GetEdepByTrackId(fGamma1->GetTrackId());
    if (fGamma2) g2evis   = digitVectorInfo.GetEdepByTrackId(fGamma2->GetTrackId());
    const double nevis    = digitVectorInfo.GetEdepByPdg(2112);
    const double pevis    = digitVectorInfo.GetEdepByPdg(2212);
    const double pipevis  = digitVectorInfo.GetEdepByPdg(+211);
    const double pimevis  = digitVectorInfo.GetEdepByPdg(-211);
    const double pipmevis = pipevis + pimevis;
    const double gmevis   = digitVectorInfo.GetEdepByPdg(22);
    const double muevis   = digitVectorInfo.GetEdepByPdg(+13)+digitVectorInfo.GetEdepByPdg(-13);
    const double total    = digitVectorInfo.GetTotalTruthEnergy();
    const double other    = total -(pi0evis+nevis+pevis+pipmevis+gmevis+muevis);
    const double shared   = digitVectorInfo.GetSharedTruthEnergy();
    const int pdg         = digitVectorInfo.GetMostEvisPdg();
    const double frac     = digitVectorInfo.GetEdepByPdg(pdg)/total;

    event->setDoubleData(Form("g%1dpi0evis",index),     pi0evis);
    event->setDoubleData(Form("g%1dg1evis",index),      g1evis);
    event->setDoubleData(Form("g%1dg2evis",index),      g2evis);
    event->setDoubleData(Form("g%1dneutronevis",index), nevis);
    event->setDoubleData(Form("g%1dprotonevis",index),  pevis);
    event->setDoubleData(Form("g%1dpipevis",index),     pipevis);
    event->setDoubleData(Form("g%1dpimevis",index),     pimevis);
    event->setDoubleData(Form("g%1dgmevis",index),      gmevis);
    event->setDoubleData(Form("g%1dmuevis",index),      muevis);
    event->setDoubleData(Form("g%1dotherevis",index),   other);
    event->setDoubleData(Form("g%1dtotalevis",index),   total);
    event->setDoubleData(Form("g%1dsharedevis",index),  shared);
    event->setDoubleData(Form("g%1dmostevisfrac",index),frac);
    event->setIntData(Form("g%1dmostevispdg",index),pdg);
}

void CCDeltaPlusAna::FillTrajectoryMap() const
{
        /* The trajectory map must be emptied for every new physics event */
    if (fTrajectoryMap.empty()) {
        Minerva::TG4Trajectories* trajectories(NULL);
        if (exist<Minerva::TG4Trajectories>(Minerva::TG4TrajectoryLocation::Default)) {
                trajectories = get<Minerva::TG4Trajectories>(Minerva::TG4TrajectoryLocation::Default);
        }

        if (!trajectories) return;
        
        for (Minerva::TG4Trajectories::iterator t = trajectories->begin();
            t != trajectories->end(); ++t) {
            Minerva::TG4Trajectory* traj = *t;
            debug() << "\t TrackId: " << traj->GetTrackId() << " parentId: " << traj->GetParentId() << endmsg;
            fTrajectoryMap[traj->GetTrackId()] = traj;
        }
    }
    
    if (fTrajectoryMap.empty()) {
        warning() << "TrajectoryMap empty " << endmsg;
        return;
    }
    
    debug()<< "Trajectory map: " << fTrajectoryMap.size() << endmsg;
}


void CCDeltaPlusAna::SummarizeTruth(Minerva::PhysicsEvent* event) const
{
    fPizero = NULL;
    fGamma1 = NULL;
    fGamma2 = NULL;
    
    int pizeroCount = 0;
    int npip = 0;
    int npim = 0;
    int pionCount = 0;
    int othermesonCount = 0;
    int protonCount = 0;
    int neutronCount = 0;
    double totalProtonKE = 0.0;
    double totalNeutronKE = 0.0;

    std::vector<const Minerva::TG4Trajectory*> fsPizeros;
    std::set<int> secondaryParticles;
    for (std::map<int,Minerva::TG4Trajectory*>::const_iterator t = fTrajectoryMap.begin();
        t != fTrajectoryMap.end(); ++t) {
        const Minerva::TG4Trajectory* traj = t->second;
        const int pdg = traj->GetPDGCode();
        const int abspdg = std::abs(pdg);
        const double E = traj->GetInitialMomentum().E();
        
        if (traj->GetParentId() > 0) {
            secondaryParticles.insert(pdg);
            continue;
        }

        if (abspdg == 11) continue;
        
        if (pdg == 22) continue;
        if (pdg == 2112) {
            ++neutronCount;
            totalNeutronKE += (E-mn);
            continue;
        }
        if (pdg == 2212) {
            ++protonCount;
            totalProtonKE += (E-mp);
            continue;
        }

        if (pdg == 111) {
            ++pizeroCount;
            fsPizeros.push_back(traj);
        }

        if (pdg == 211) ++npip;
        if (pdg ==-211) ++npim;
        if (abspdg == 211) ++pionCount;
        
        const bool not_pion = (pdg != 111) && (abspdg != 211);
        if (not_pion) {
            TParticlePDG* particle = TDatabasePDG::Instance()->GetParticle(pdg);
            assert(particle);
            if (std::string(particle->ParticleClass()) == "Meson") ++othermesonCount;
            
        }
                
    }

    const Minerva::TG4Trajectory* pizero(NULL);
    if (!fsPizeros.empty()) {
        std::sort(fsPizeros.begin(), fsPizeros.end(), greaterTotalEnergy());
        pizero = fsPizeros.front();
    }
    
    fPizero = pizero;
    
        // Another loop to get the kinematics of daughter gamma's from the pi0 decay
    std::vector<const Minerva::TG4Trajectory*> pizeroDaughters;
    for (std::map<int,Minerva::TG4Trajectory*>::const_iterator t = fTrajectoryMap.begin();
        t != fTrajectoryMap.end(); ++t) {
        const TG4Trajectory* traj = t->second;
            /* Skip final-state particles */
        if (traj->GetParentId() == 0) continue;
        const TG4Trajectory* parent = fTrajectoryMap[traj->GetParentId()];
        if (parent && parent == pizero) pizeroDaughters.push_back(traj);

    }

    event->setIntData("npi0", pizeroCount);
    event->setIntData("npip", npip);
    event->setIntData("npim", npim);
    event->setIntData("npipm",pionCount);
    event->setIntData("nmeson",othermesonCount);
    event->setIntData("npi02",std::count(secondaryParticles.begin(),
                                        secondaryParticles.end(), 111));
    event->setIntData("np",protonCount);
    event->setIntData("nn",neutronCount);
    event->setDoubleData("pke",totalProtonKE);
    event->setDoubleData("nke",totalNeutronKE);

    event->setIntData("dmode",pizeroDaughters.size()== 2 ? 0 : 1);

    if (pizero) {
        event->setDoubleData("pitheta0",pizero->GetInitialMomentum().Theta()*rad2deg);
        event->setDoubleData("piphi0",  pizero->GetInitialMomentum().Phi()*rad2deg);
        event->setDoubleData("pienergy0",pizero->GetInitialMomentum().E());
        
        std::vector<double> pimom0;
        pimom0.push_back(pizero->GetInitialMomentum().Px());
        pimom0.push_back(pizero->GetInitialMomentum().Py());
        pimom0.push_back(pizero->GetInitialMomentum().Pz());
        pimom0.push_back(pizero->GetInitialMomentum().P());

        event->setContainerDoubleData("pimom0",pimom0);
    }

        /* Guard against decay outside the detector! */
    if (pizero && !pizeroDaughters.empty()) {
        std::sort(pizeroDaughters.begin(),pizeroDaughters.end(), greaterTotalEnergy());
        fGamma1 = pizeroDaughters.front();
        fGamma2 = pizeroDaughters.back();
        const Minerva::TG4Trajectory* g1 = pizeroDaughters.front(); /* Variables to use locally */
        const Minerva::TG4Trajectory* g2 = pizeroDaughters.back();
        const Gaudi::LorentzVector g1p4 = g1->GetInitialMomentum();
        const Gaudi::LorentzVector g2p4 = g2->GetInitialMomentum();
        assert(g1p4.E()>= g2p4.E());
        const TVector3 g1mom(g1p4.Px(),g1p4.Py(),g1p4.Pz());
        const TVector3 g2mom(g2p4.Px(),g2p4.Py(),g2p4.Pz());

        event->setDoubleData("g1theta0",g1p4.Theta()*rad2deg);
        event->setDoubleData("g1phi0",  g1p4.Phi()*rad2deg);
        event->setDoubleData("g1e0",    g1p4.E());
        event->setDoubleData("g2theta0",g2p4.Theta()*rad2deg);
        event->setDoubleData("g2phi0",  g2p4.Phi()*rad2deg);
        event->setDoubleData("g2e0",    g2p4.E());

        event->setDoubleData("oangle0",g1mom.Angle(g2mom)*rad2deg);

        const double s1x = g1mom.Unit().X();
        const double s1z = g1mom.Unit().Z();
        const double s2x = g2mom.Unit().X();
        const double s2z = g2mom.Unit().Z();

        const double g1theta_x = std::atan2(s1x,s1z)*rad2deg;
        const double g2theta_x = std::atan2(s2x,s2z)*rad2deg;
        double delta = std::max(g1theta_x,g2theta_x)-std::min(g1theta_x,g2theta_x);
        if (delta > 180.0) delta = 360-delta; 
        event->setDoubleData("oangle0x", delta);
        
        std::vector<double> g1momVec;
        g1momVec.push_back(g1p4.Px());
        g1momVec.push_back(g1p4.Py());
        g1momVec.push_back(g1p4.Pz());
        g1momVec.push_back(g1p4.P());

        std::vector<double> g2momVec;
        g2momVec.push_back(g2p4.Px());
        g2momVec.push_back(g2p4.Py());
        g2momVec.push_back(g2p4.Pz());
        g2momVec.push_back(g2p4.P());

        event->setContainerDoubleData("g1mom0",g1momVec);
        event->setContainerDoubleData("g2mom0",g2momVec);

        std::vector<double> g1convPos;
        g1convPos.push_back(g1->GetFinalPosition().X());
        g1convPos.push_back(g1->GetFinalPosition().Y());
        g1convPos.push_back(g1->GetFinalPosition().Z());

        std::vector<double> g2convPos;
        g2convPos.push_back(g2->GetFinalPosition().X());
        g2convPos.push_back(g2->GetFinalPosition().Y());
        g2convPos.push_back(g2->GetFinalPosition().Z());

        event->setContainerDoubleData("g1convpos",g1convPos);
        event->setContainerDoubleData("g2convpos",g2convPos);

            /* Type conversion to be used with isInside() */
        Gaudi::XYZPoint g1Stop(g1->GetFinalPosition().X(),
                            g1->GetFinalPosition().Y(),
                            g1->GetFinalPosition().Z());
        Gaudi::XYZPoint g2Stop(g2->GetFinalPosition().X(),
                            g2->GetFinalPosition().Y(),
                            g2->GetFinalPosition().Z());
        
        int g1convidet = m_idDet->isInside(g1Stop) ? 1 : 0;
        int g2convidet = m_idDet->isInside(g2Stop) ? 1 : 0;
        event->setIntData("g1convidet",g1convidet);
        event->setIntData("g2convidet",g2convidet);
        
        double g1convdist = (g1->GetFinalPosition().Vect() - g1->GetInitialPosition().Vect()).R();
        double g2convdist = (g2->GetFinalPosition().Vect() - g2->GetInitialPosition().Vect()).R();
        event->setDoubleData("g1convdist",g1convdist);
        event->setDoubleData("g2convdist",g2convdist);
        
    }
}

void CCDeltaPlusAna::SummarizeHadronTruth(Minerva::PhysicsEvent* event) const
{ 
        // Find Michel electron from pi+/- --> mu+/- decays
        // 1) Find secondary muons
    debug() << "Find mu+/- decays<oaltinok_version>" << endmsg;
    std::vector<const Minerva::TG4Trajectory*> secondaryMuons;
    for (std::map<int,Minerva::TG4Trajectory*>::const_iterator t = fTrajectoryMap.begin();
        t != fTrajectoryMap.end(); ++t) {
        const TG4Trajectory* traj = t->second;
        const int abspdg = std::abs(traj->GetPDGCode());
        if (traj->GetParentId() == 0) continue;          /* Skip final-state particles */
        if (abspdg == 22) continue;                      /* Skip gamma to speed up */
        if (abspdg == 13) secondaryMuons.push_back(traj);

    }
        // 2) Find the secondary muon daughters
    std::map<const TG4Trajectory*,std::vector<const TG4Trajectory*> > muonDaughterMap;
    for (std::vector<const Minerva::TG4Trajectory*>::const_iterator muon = secondaryMuons.begin();
        muon != secondaryMuons.end(); ++muon) {
        std::vector<const Minerva::TG4Trajectory*> decayDaughters;
        for (std::map<int,Minerva::TG4Trajectory*>::const_iterator t = fTrajectoryMap.begin();
            t != fTrajectoryMap.end(); ++t) {
            const TG4Trajectory* traj = t->second;
            if (traj->GetParentId() == 0) continue;
            if (traj->GetPDGCode() == 22) continue;
            const Minerva::TG4Trajectory* parent = fTrajectoryMap[traj->GetParentId()];
            if (parent && parent == *muon &&
                (traj->GetProcessName() == "muMinusCaptureAtRest" ||
                traj->GetProcessName() == "Decay")) {
                decayDaughters.push_back(traj);
            }
        }

        muonDaughterMap.insert(std::make_pair(*muon,decayDaughters));
    }
        // 3) Check if the daughers are consistent with decay or capture
    unsigned int nmupdecay   = 0;
    unsigned int nmumdecay   = 0;
    unsigned int nmumcapture = 0;
    debug() << "Total decays: " << muonDaughterMap.size() << endmsg;
    for (std::map<const TG4Trajectory*, std::vector<const TG4Trajectory*> >::iterator m
            = muonDaughterMap.begin();
        m != muonDaughterMap.end(); ++m) {
        const TG4Trajectory* parent = m->first;
        const std::vector<const TG4Trajectory*>& daughters = m->second;
        const TG4Trajectory* positron(NULL);
        const TG4Trajectory* electron(NULL);
        const TG4Trajectory* nue(NULL);
        const TG4Trajectory* antinue(NULL);
        const TG4Trajectory* numu(NULL);
        const TG4Trajectory* antinumu(NULL);
        debug() << "\tParent: " << parent->GetPDGCode() << endmsg;
        for (std::vector<const TG4Trajectory*>::const_iterator d = daughters.begin();
            d != daughters.end(); ++d) {
            const int pdg = (*d)->GetPDGCode();
            debug() << "\t\tDaughter: " << pdg << endmsg;
            if (pdg == +11) electron = *d;
            if (pdg == -11) positron = *d;
            if (pdg == +12) nue      = *d;
            if (pdg == -12) antinue  = *d;
            if (pdg == +14) numu     = *d;
            if (pdg == -14) antinumu = *d;
        }

        if (parent->GetPDGCode() == 13 && electron && antinue && numu) {
            ++nmumdecay;
            debug() << "\t\t ==> mu- decay<oaltinok_version>" << endmsg;
            std::vector<double> pos;
            pos.push_back(electron->GetInitialPosition().X());
            pos.push_back(electron->GetInitialPosition().Y());
            pos.push_back(electron->GetInitialPosition().Z());
            pos.push_back(electron->GetInitialPosition().T());
            std::vector<double> mom;
            mom.push_back(electron->GetInitialMomentum().Px());
            mom.push_back(electron->GetInitialMomentum().Py());
            mom.push_back(electron->GetInitialMomentum().Pz());
            mom.push_back(electron->GetInitialMomentum().E());

            event->setContainerDoubleData("michel_pos",pos);
            event->setContainerDoubleData("michel_mom",mom);
        }

        if (parent->GetPDGCode() == 13 && !electron && numu) {
            ++nmumcapture;
            debug() << "\t\t ==> mu- capture<oaltinok_version>" << endmsg;
        }
        
        if (parent->GetPDGCode() ==-13 && positron && nue && antinumu) {
            ++nmupdecay;
            debug() << "\t\t ==> mu+ decay<oaltinok_version>" << endmsg;
        }
    }

    event->setIntData("nmupdecay",nmupdecay);
    event->setIntData("nmumdecay",nmumdecay);
    event->setIntData("nmumcapture",nmumcapture);
    
    debug() << "Find pi+/- reactions<oaltinok_version>" << endmsg;
    std::vector<const TG4Trajectory*> chargedPions;
    std::vector<const TG4Trajectory*> charged2Pions;
    for (std::map<int,Minerva::TG4Trajectory*>::const_iterator t = fTrajectoryMap.begin();
        t != fTrajectoryMap.end(); ++t) {
        const TG4Trajectory* traj = t->second;
        const int pdg = traj->GetPDGCode();
        if (pdg == 22) continue;
        if (pdg == 11) continue;
        if (pdg ==-11) continue;
        const int abspdg = std::abs(pdg);
        if (abspdg == 211) {      
                /* Final-state pi- */
            if (traj->GetParentId() == 0) { 
                chargedPions.push_back(traj);
                debug() << "\tFound final-state pi: " << traj->GetTrackId() << endmsg;
                continue;
            }
                /* Secondary pi-, but reject scatterd one */
            else {
                
                const TG4Trajectory* parent(NULL);
                if (traj->GetParentId() > 0) parent = fTrajectoryMap[traj->GetParentId()];
                if (parent && std::abs(parent->GetPDGCode()) != 211) {
                    debug() << "\tFound secondary pi:   " << traj->GetTrackId() << " " << traj->GetParentId()
                            << endmsg;
                    charged2Pions.push_back(traj);
                }
            }
        }
    }

    unsigned int npipcapture   =-1;
    unsigned int npipinelastic = 0;
    unsigned int npipdecay     = 0;
    unsigned int npimcapture   = 0;
    unsigned int npiminelastic = 0;
    unsigned int npimdecay     = 0;
    for (std::vector<const TG4Trajectory*>::const_iterator p = chargedPions.begin();
        p != chargedPions.end(); ++p) {
        const int pdg = (*p)->GetPDGCode();
        if (pdg == 211) debug() << "\tpi+ daughters<oaltinok_version>" << endmsg;
        if (pdg ==-211) debug() << "\tpi- daughters<oaltinok_version>" << endmsg;
        std::string proc("");
        for (std::map<int,Minerva::TG4Trajectory*>::const_iterator t = fTrajectoryMap.begin();
            t != fTrajectoryMap.end(); ++t) {
            const TG4Trajectory* traj = t->second;
            if (traj->GetParentId() == 0) continue;
            const TG4Trajectory* parent = fTrajectoryMap[traj->GetParentId()];
            if (parent == *p) {
                const double deltaT = traj->GetInitialPosition().T() - parent->GetInitialPosition().T();
                debug() << "\t\t" << setw(10) << traj->GetPDGCode() << " " << traj->GetProcessName()
                        << setw(10) << deltaT
                        << endmsg;
                proc += traj->GetProcessName();
            }
        }
        
        if      (proc.find("CHIPSNuclearCaptureAtRest") != std::string::npos) ++npimcapture;
        else if (proc.find("PionMinusInelastic") != std::string::npos)        ++npiminelastic;
        else if (proc.find("Decay") != std::string::npos && pdg == -211)      ++npimdecay;
        else if (proc.find("PionPlusInelastic") != std::string::npos)         ++npipinelastic;
        else if (proc.find("Decay") != std::string::npos && pdg == +211)      ++npipdecay;
        else {}
    }

    event->setIntData("npimcapture",   npimcapture);
    event->setIntData("npiminelastic", npiminelastic);
    event->setIntData("npimdecay",     npimdecay);
    event->setIntData("npipcapture",   npipcapture);
    event->setIntData("npipinelastic", npipinelastic);
    event->setIntData("npipdecay",     npipdecay);

    event->setIntData("npim2",std::count_if(chargedPions.begin(),chargedPions.end(),PdgEqual(-211)));
    event->setIntData("npip2",std::count_if(chargedPions.begin(),chargedPions.end(),PdgEqual(+211)));
        

    double pimTrackLength = -1.0;
    double pipTrackLength = -1.0;
    
    std::vector<const TG4Trajectory*>::const_iterator pimIter
        = std::find_if(chargedPions.begin(), chargedPions.end(),PdgEqual(-211));
    std::vector<const TG4Trajectory*>::const_iterator pipIter
        = std::find_if(chargedPions.begin(), chargedPions.end(),PdgEqual(+211));
    
    if (pimIter != chargedPions.end()) {
        const TG4Trajectory* pim = *pimIter;
        pimTrackLength = (pim->GetFinalPosition().Vect() - pim->GetInitialPosition().Vect()).R();
    }
    if (pipIter != chargedPions.end()) {
        const TG4Trajectory* pip = *pipIter;
        pipTrackLength = (pip->GetFinalPosition().Vect() - pip->GetInitialPosition().Vect()).R();
    }
    
    event->setDoubleData("pimlength",pimTrackLength);
    event->setDoubleData("piplength",pipTrackLength);
    

}

void CCDeltaPlusAna::GammaEdepInfo(Minerva::PhysicsEvent* event) const
{
    Minerva::MCHits* idhits(NULL);
    if (exist<Minerva::MCHits>(Minerva::MCHitLocation::ID)) {
        idhits = get<Minerva::MCHits>(Minerva::MCHitLocation::ID);
    }

    Minerva::MCHits* odhits(NULL);
    if (exist<Minerva::MCHits>(Minerva::MCHitLocation::OD)) {
        odhits = get<Minerva::MCHits>(Minerva::MCHitLocation::OD);
    }
    
    SmartRefVector<Minerva::MCHit> idhitvec;
    SmartRefVector<Minerva::MCHit> odhitvec;
    SmartRefVector<Minerva::MCHit> nukehits;
    SmartRefVector<Minerva::MCHit> trkrhits;
    SmartRefVector<Minerva::MCHit> ecalhits;
    SmartRefVector<Minerva::MCHit> hcalhits;
    SmartRefVector<Minerva::MCHit> sidehits;
    SmartRefVector<Minerva::MCHit> other_subdet_hits;

    for (Minerva::MCHits::iterator h = idhits->begin();
        h != idhits->end(); ++h) {
        idhitvec.push_back(*h);
        const Minerva::MCHit* hit = *h;
        const double x = (hit->StartX() + hit->StopX())/2.0;
        const double y = (hit->StartY() + hit->StopY())/2.0;
        const double z = (hit->StartZ() + hit->StopZ())/2.0;

        Gaudi::XYZPoint thePoint(x,y,z);
        const Minerva::DeSubdet* subdet = m_idDet->getDeSubdet(thePoint);
        Minerva::SubdetID::Subdet subdetId = subdet->getSubdetID().subdet();
        if      (subdetId == Minerva::SubdetID::Tracker)   trkrhits.push_back(*h);
        else if (subdetId == Minerva::SubdetID::ECAL)      ecalhits.push_back(*h);
        else if (subdetId == Minerva::SubdetID::HCAL)      hcalhits.push_back(*h);
        else if (subdetId == Minerva::SubdetID::NuclTargs) nukehits.push_back(*h);
        else other_subdet_hits.push_back(*h);

            /* Ten strips near the detector edges in the Tracker and NTGT regions
            belong to the side ECALs. Store them in 'sidehits' and remove them from
            'trkrhits' or 'nukehits' to avoid double counting
            */
        if (subdetId == Minerva::SubdetID::Tracker ||
            subdetId == Minerva::SubdetID::NuclTargs) {
            Minerva::StripID stripId = subdet->getStripID(thePoint);
            unsigned int strip = stripId.strip();
            if (strip < 10 || strip > 117) {
                sidehits.push_back(*h);
                if (subdetId == Minerva::SubdetID::Tracker)   trkrhits.pop_back();
                if (subdetId == Minerva::SubdetID::NuclTargs) nukehits.pop_back();
            }
        }

    }

    std::copy(odhits->begin(), odhits->end(), std::back_inserter(odhitvec));
    
    std::cout << "HitVector: nuke " << nukehits.size() << std::endl;
    std::cout << "HitVector: trkr " << trkrhits.size() << std::endl;
    std::cout << "HitVector: ecal " << ecalhits.size() << std::endl;
    std::cout << "HitVector: hcal " << hcalhits.size() << std::endl;
    std::cout << "HitVector: side " << sidehits.size() << std::endl;
    std::cout << "HitVector: id   " << idhitvec.size() << std::endl;
    std::cout << "HitVector: od   " << odhitvec.size() << std::endl;

    if (fPizero) {
        std::cout << "Traversing tracking history..." << std::endl;
        HitVectorTruthInfo nukeHitInfo;
        nukeHitInfo.ParseTruth(nukehits,fTrajectoryMap);

        HitVectorTruthInfo trkrHitInfo;
        trkrHitInfo.ParseTruth(trkrhits,fTrajectoryMap);

        HitVectorTruthInfo ecalHitInfo;
        ecalHitInfo.ParseTruth(ecalhits,fTrajectoryMap);

        HitVectorTruthInfo hcalHitInfo;
        hcalHitInfo.ParseTruth(hcalhits,fTrajectoryMap);

        HitVectorTruthInfo sideHitInfo;
        sideHitInfo.ParseTruth(sidehits,fTrajectoryMap);

        HitVectorTruthInfo odetHitInfo;
        odetHitInfo.ParseTruth(odhitvec,fTrajectoryMap);

        HitVectorTruthInfo otherSubdetHitInfo;
        otherSubdetHitInfo.ParseTruth(other_subdet_hits,fTrajectoryMap);

        event->setDoubleData("pi0nukeedep", nukeHitInfo.GetEdepByPdg(111));
        event->setDoubleData("pi0trkredep", trkrHitInfo.GetEdepByPdg(111));
        event->setDoubleData("pi0ecaledep", ecalHitInfo.GetEdepByPdg(111));
        event->setDoubleData("pi0hcaledep", hcalHitInfo.GetEdepByPdg(111));
        event->setDoubleData("pi0sideedep", sideHitInfo.GetEdepByPdg(111));
        event->setDoubleData("pi0idetedep",
                            nukeHitInfo.GetEdepByPdg(111) +
                            trkrHitInfo.GetEdepByPdg(111) +
                            ecalHitInfo.GetEdepByPdg(111) +
                            hcalHitInfo.GetEdepByPdg(111) +
                            sideHitInfo.GetEdepByPdg(111));
        
        event->setDoubleData("pi0odetedep", odetHitInfo.GetEdepByPdg(111));
        event->setDoubleData("pi0othersubdetedep",otherSubdetHitInfo.GetEdepByPdg(111));

        const double pi0ecalo
            = nukeHitInfo.GetEdepByPdg(111)/Constants().k_trkr  /* Treat NTGT as Tracker */
            + trkrHitInfo.GetEdepByPdg(111)/Constants().k_trkr
            + ecalHitInfo.GetEdepByPdg(111)/Constants().k_ecal
            + hcalHitInfo.GetEdepByPdg(111)/Constants().k_hcal
            + sideHitInfo.GetEdepByPdg(111)/Constants().k_ecal; /* Treat side ECAL as ECAL */
        event->setDoubleData("pi0ecalo", pi0ecalo);
            

            /* pi0 --> gamma gamma  decay */
        if (fGamma1 && fGamma2) {
            const int trackId1 = fGamma1->GetTrackId();
            const int trackId2 = fGamma2->GetTrackId();
            event->setDoubleData("g1nukeedep",nukeHitInfo.GetEdepByTrackId(trackId1));
            event->setDoubleData("g2nukeedep",nukeHitInfo.GetEdepByTrackId(trackId2));

            event->setDoubleData("g1trkredep",trkrHitInfo.GetEdepByTrackId(trackId1));
            event->setDoubleData("g2trkredep",trkrHitInfo.GetEdepByTrackId(trackId2));

            event->setDoubleData("g1ecaledep",ecalHitInfo.GetEdepByTrackId(trackId1));
            event->setDoubleData("g2ecaledep",ecalHitInfo.GetEdepByTrackId(trackId2));
            
            event->setDoubleData("g1hcaledep",hcalHitInfo.GetEdepByTrackId(trackId1));
            event->setDoubleData("g2hcaledep",hcalHitInfo.GetEdepByTrackId(trackId2));
            
            event->setDoubleData("g1sideedep",sideHitInfo.GetEdepByTrackId(trackId1));
            event->setDoubleData("g2sideedep",sideHitInfo.GetEdepByTrackId(trackId2));

            event->setDoubleData("g1idetedep",
                                nukeHitInfo.GetEdepByTrackId(trackId1) +
                                trkrHitInfo.GetEdepByTrackId(trackId1) +
                                ecalHitInfo.GetEdepByTrackId(trackId1) +
                                hcalHitInfo.GetEdepByTrackId(trackId1) +
                                sideHitInfo.GetEdepByTrackId(trackId1)
                                );
            event->setDoubleData("g2idetedep",
                                nukeHitInfo.GetEdepByTrackId(trackId2) +
                                trkrHitInfo.GetEdepByTrackId(trackId2) +
                                ecalHitInfo.GetEdepByTrackId(trackId2) +
                                hcalHitInfo.GetEdepByTrackId(trackId2) +
                                sideHitInfo.GetEdepByTrackId(trackId2)
                                );
            
            event->setDoubleData("g1odetedep",odetHitInfo.GetEdepByTrackId(trackId1));
            event->setDoubleData("g2odetedep",odetHitInfo.GetEdepByTrackId(trackId2));

            event->setDoubleData("g1othersubdetedep",otherSubdetHitInfo.GetEdepByTrackId(trackId1));
            event->setDoubleData("g2othersubdetedep",otherSubdetHitInfo.GetEdepByTrackId(trackId2));

            const double g1ecalo
                = nukeHitInfo.GetEdepByTrackId(trackId1)/Constants().k_trkr
                + trkrHitInfo.GetEdepByTrackId(trackId1)/Constants().k_trkr
                + ecalHitInfo.GetEdepByTrackId(trackId1)/Constants().k_ecal
                + hcalHitInfo.GetEdepByTrackId(trackId1)/Constants().k_hcal
                + sideHitInfo.GetEdepByTrackId(trackId1)/Constants().k_ecal;
            const double g2ecalo
                = nukeHitInfo.GetEdepByTrackId(trackId2)/Constants().k_trkr
                + trkrHitInfo.GetEdepByTrackId(trackId2)/Constants().k_trkr
                + ecalHitInfo.GetEdepByTrackId(trackId2)/Constants().k_ecal
                + hcalHitInfo.GetEdepByTrackId(trackId2)/Constants().k_hcal
                + sideHitInfo.GetEdepByTrackId(trackId2)/Constants().k_ecal; 

            event->setDoubleData("g1ecalo", g1ecalo);
            event->setDoubleData("g2ecalo", g2ecalo);
        }

            /* While we are at it, get the energy deposit by primary neutrons */
        event->setDoubleData("neutronnukeedep", nukeHitInfo.GetEdepByPdg(2112));
        event->setDoubleData("neutrontrkredep", trkrHitInfo.GetEdepByPdg(2112));
        event->setDoubleData("neutronecaledep", ecalHitInfo.GetEdepByPdg(2112));
        event->setDoubleData("neutronhcaledep", hcalHitInfo.GetEdepByPdg(2112));
        event->setDoubleData("neutronsideedep", sideHitInfo.GetEdepByPdg(2112));
        event->setDoubleData("neutronidetedep",
                            nukeHitInfo.GetEdepByPdg(2112) +
                            trkrHitInfo.GetEdepByPdg(2112) +
                            ecalHitInfo.GetEdepByPdg(2112) +
                            hcalHitInfo.GetEdepByPdg(2112) +
                            sideHitInfo.GetEdepByPdg(2112));
        
        event->setDoubleData("neutronodetedep", odetHitInfo.GetEdepByPdg(2112));
        event->setDoubleData("neutronothersubdetedep",otherSubdetHitInfo.GetEdepByPdg(2112));


    }
    
}

double CCDeltaPlusAna::CalcTrackLength(const SmartRef<Minerva::Track>& track) const 
{
    const Minerva::Track::NodeContainer& nodes = track->nodes();

    double trackLength = 0.0;
    Minerva::Track::NodeContainer::const_iterator nextToEnd = nodes.end();
    --nextToEnd;
    for (Minerva::Track::NodeContainer::const_iterator node = nodes.begin();
        node != nodes.end(); ++node) {
            /* Handle the last node */
        if (node != nextToEnd) {
            Minerva::Track::NodeContainer::const_iterator next = node; ++next;
            TVector3 currentPos((*node)->state().x(),(*node)->state().y(), (*node)->z());
            TVector3 nextPos((*next)->state().x(),(*next)->state().y(),(*next)->z());
            trackLength += (nextPos-currentPos).Mag();
        }
        
    }

    return trackLength;
}

double CCDeltaPlusAna::CalcMinBlobSeparation(const Minerva::IDBlob* blob,
                                        const SmartRef<Minerva::Vertex>& vertex) const
{
    const Gaudi::XYZPoint& pos = vertex->position();
    const double x0 = pos.X();
    const double y0 = pos.Y();
    const double z0 = pos.Z();
    double u0 = m_mathTool->calcUfromXY(x0,y0);
    double v0 = m_mathTool->calcVfromXY(x0,y0);

    const SmartRefVector<Minerva::IDCluster>& clusters = blob->clusters();
    double dmin = 1e6;
    for (SmartRefVector<Minerva::IDCluster>::const_iterator c = clusters.begin();
        c != clusters.end(); ++c) {
        const double energy = (*c)->energy();
        const double z      = (*c)->z();
        const double t      = (*c)->position();
            /* Skip low-energy clusters, normal MIP cluster has 1.7*2 MeV = 3.4 MeV */
        if (energy < 1.0) continue;

            /* Tranverse position: X,U,or V */
        double t0  = 0.0; 
        switch ((*c)->view()) {

            case Minerva::IDCluster::X:
                t0 = x0;
                break;

            case Minerva::IDCluster::U:
                t0 = u0;
                break;

            case Minerva::IDCluster::V:
                t0 = v0;
                break;

            default:
                warning() << "Invalid view<oaltinok_version>" << endmsg;
                exit(1);
                
        }
            /* Distance from this cluster to the fitted vertex */
        double distance = sqrt(pow(z-z0,2) + pow(t-t0,2));
        dmin = std::min(dmin,distance);
    }

    return dmin;
}

SmartRefVector<Minerva::IDCluster>
CCDeltaPlusAna::FilterInSphereClusters(const SmartRef<Minerva::Vertex>& vertex,
                                    const SmartRefVector<Minerva::IDCluster>& clusters,
                                    const double sphereRadius) const
{
    const Gaudi::XYZPoint& pos = vertex->position();
    const double x0 = pos.X();
    const double y0 = pos.Y();
    const double z0 = pos.Z();
    double u0 = m_mathTool->calcUfromXY(x0,y0);
    double v0 = m_mathTool->calcVfromXY(x0,y0);

    SmartRefVector<Minerva::IDCluster> sphereClusters;
    for (SmartRefVector<Minerva::IDCluster>::const_iterator c = clusters.begin();
        c != clusters.end(); ++c) {
        const double energy = (*c)->energy();
        const double z      = (*c)->z();
        const double t      = (*c)->position();
            /* Skip low-energy clusters, normal MIP cluster has 1.7*2 MeV = 3.4 MeV */
        if (energy < 1.0) continue;
        
            /* Tranverse position: X,U,or V */
        double t0  = 0.0; 
        switch ((*c)->view()) {
            
            case Minerva::IDCluster::X:
                t0 = x0;
                break;
                
            case Minerva::IDCluster::U:
                t0 = u0;
                break;
                
            case Minerva::IDCluster::V:
                t0 = v0;
                break;
                
            default:
                warning() << "Invalid view<oaltinok_version>" << endmsg;
                exit(1);
                
        }
            /* Distance from this cluster to the primary vertex */
        double radius = sqrt(pow(z-z0,2) + pow(t-t0,2));
        if (radius < sphereRadius) sphereClusters.push_back(*c);
    }

    return sphereClusters;
    
}

/* Do a line fit for blob digits in the X view. Force the fitted line to
go through the primary vertex
*/

std::pair<int,double>
CCDeltaPlusAna::OneParLineFitBlob(const Minerva::IDBlob* blob, 
                            const SmartRef<Minerva::Vertex>& vertex,
                            const SmartRef<Minerva::Track>& muonTrack) const
{
        /* Get all X blob digits via clusters since clusters have 'view' info */
    SmartRefVector<Minerva::IDDigit> digits;
    const SmartRefVector<Minerva::IDCluster>& clusters = blob->clusters();
    for (SmartRefVector<Minerva::IDCluster>::const_iterator c = clusters.begin();
        c != clusters.end(); ++c) {
        if ((*c)->view() != Minerva::IDCluster::X) continue;
        const SmartRefVector<Minerva::IDDigit>& cdigits = (*c)->digits();
        std::copy(cdigits.begin(), cdigits.end(), std::back_inserter(digits));
    }
        // Filter digits in the X view here
    std::vector<LineFit::Point> data;
    double total = 0.0;
    for (SmartRefVector<Minerva::IDDigit>::const_iterator d = digits.begin();
        d != digits.end(); ++d) {

        if ((*d)->pe() < 3.0) continue;
        
        const Minerva::DePlane* plane = m_idDet->getDePlane((*d)->stripid());
        assert(plane);
        if (!plane) continue;

        total += (*d)->normEnergy();
        
        const double z = plane->getZCenter();
        const double x = plane->getTPos((*d)->stripid());
        const double w = (*d)->normEnergy();
        data.push_back(LineFit::Point(z,x,0.0,w));
    
    }
    
    const Gaudi::XYZPoint& pos = vertex->position();
    const double x0 = pos.X();
    const double z0 = pos.Z();

        /* Track slope near the vertex to calculate x(z) */
    const double tkx = muonTrack->nearestNode(z0)->state().ax();

    const double zmin  = -20.0;      /* Scan range in Z +/- 2 cm around the vertex */
    const double zmax  = +20.0;
    const double zstep = 5.0;        /* Scan step */

    std::vector<double> fval_vector; /* fit chi2 of the 9 scanned points */
    double fval_min = 1.e12;         /* The minimum chi2 from the scanned points */
    unsigned int ndof = digits.size();
    for (double dz = zmin; dz <= zmax; dz += zstep) {
        const double z_i = z0 + dz;
        const double x_i = tkx*(z_i - z0) + x0;
            /*
        std::cout << "\tScan i: " << (int)(dz-zmin)/zstep << " "
                << z_i << " " << x_i
                << std::endl;
            */
        TFitterMinuit* minuit = new TFitterMinuit;
        OneParLineFit fittingFunction(data);
        fittingFunction.SetFixedPoint(z_i,x_i);
        minuit->SetMinuitFCN(&fittingFunction);

            /* Set paramter (index,name,init,err,range) */
        minuit->SetParameter(0,"kx",1.0,0.1,-1000,1000);

        minuit->SetPrintLevel(-1);
        minuit->CreateMinimizer();

        int errno = minuit->Minimize();
        if (errno > 0) {
            std::cerr << "Fitting error: " << errno << std::endl;
            fval_vector.push_back(1.e9);
            continue;
        }

            //const double kfit = minuit->GetParameter(0);
        const double fval = fittingFunction.GetFCN();
        
        fval_min = std::min(fval,fval_min);
    }
    std::cout << "\tLine fitting: chi2/ndof/total: " << fval_min << "/" << ndof
            << "/" << total 
            << std::endl;
    
    return std::make_pair<int,double>(ndof,fval_min/total);
}

double CCDeltaPlusAna::CalcDistanceFromBlobAxisToVertex(const Minerva::IDBlob* blob,
                                                    const SmartRef<Minerva::Vertex>& vertex) const
{
    const Gaudi::XYZPoint& x0 = vertex->position();
    const Gaudi::XYZPoint& x1 = blob->startPoint();

    const Gaudi::XYZVector& s = blob->direction();
    
    const double d0 = 1000.0;
        /* Second point on the blob axis 1m away from the x1 */
    Gaudi::XYZPoint x2(x1.X() + d0*s.X(),
                    x1.Y() + d0*s.Y(),
                    x1.Z() + d0*s.Z()); 

        /* Calculate the distance from the vertex to the blob axis defined by x1 and x2
        The line is defined by x1 and x2.
        The point of interest is x0
        Then the distance is:
        d = |(x0-x1)x(x0-x2)|/|x2-x1| */
    Gaudi::XYZVector x01 = x0-x1;
    Gaudi::XYZVector x02 = x0-x2;
    const double dvtx = (x01.Cross(x02)).R()/(x2-x1).R();
    std::cout << "\tDistance_to_vertex: " << dvtx << std::endl;

        /* Self-test using another point on the line 2m away from the x1*/
    const double d1 = 2000.0;
    Gaudi::XYZPoint x3(x1.X() + d1*s.X(),
                    x1.Y() + d1*s.Y(),
                    x1.Z() + d1*s.Z());
    
    Gaudi::XYZVector x31 = x3-x1;
    Gaudi::XYZVector x32 = x3-x2;
    const double dvtx2 = (x31.Cross(x32)).R()/(x2-x1).R();
    std::cout << "\t Self-test distance: " << dvtx2 << std::endl;
    assert(dvtx2 < 1.e-6);

    return dvtx;
}


double CCDeltaPlusAna::CalcDistanceFromVertexToExiting(const Minerva::IDBlob* blob,
                                                    const SmartRef<Minerva::Vertex>& vertex) const
{
    const Gaudi::XYZPoint& pos = vertex->position();
    const Gaudi::XYZVector& s = blob->direction();
    
    double d_step = 50.0; /* Accuracy on dmax = dstep/2 */
    double dmax = 0.0;    /* Distance from vertex to where the blob would exit */
    unsigned int nstep = 0;
    for (double d = 0.0; ; d += d_step) {
        dmax = d;
        ++nstep;
        Gaudi::XYZPoint r = pos;
        r += d*s;

            /* Not side exiting blob */
        if (r.Z() < 5000.0 || r.Z() > 9000.0) break;

        if (!InsideHexagon(r.X(),r.Y(),2200.0)) break;

    }

    std::cout << "\t nstep before exiting:    " << nstep << std::endl;
    std::cout << "\t Distance before exiting: " << dmax << std::endl;
    
    return dmax;
}

//=======================================================================
//  Finalize
//=======================================================================
StatusCode CCDeltaPlusAna::finalize()
{
    debug() << "CCDeltaPlusAna::finalize()" << endmsg;
    
        // finalize the base class.
    StatusCode sc = this->MinervaAnalysisTool::finalize();
    if( sc.isFailure() ) return Error( "Failed to finalize!", sc );
    
    return sc;

}

void CCDeltaPlusAna::PrintDigitVectorInfo(const SmartRefVector<Minerva::IDDigit>& digits) const
{
    info() << "Total blob digits: " << digits.size() << endmsg;
    for (SmartRefVector<Minerva::IDDigit>::const_iterator d = digits.begin();
        d != digits.end(); ++d) {
        info() << "(module,strip,energy) : " << (*d)->module() << " "
            << setw(5) << (*d)->strip() << " "
            << setw(6) << (*d)->normEnergy()
            << endmsg;
            
    }
}

//=======================================================================
//  ODActivity
//=======================================================================
StatusCode CCDeltaPlusAna::ODActivity( Minerva::PhysicsEvent *event, std::vector<Minerva::IDBlob*> idBlobs ) const
{
    debug() << " Starting ODActivity::EnergeticTower " << endmsg;
    
        // get origin of muon
    Gaudi::XYZTVector muon_position = m_MuonParticle->startPos();
    
    double od_energy[6] = {0,0,0,0,0,0}; 
    int od_nclusters[6] = {0,0,0,0,0,0}; 
    double od_time[6] = {0,0,0,0,0,0};
    double od_pe[6] = {0,0,0,0,0,0};
    int frame_low = 1000;
    int frame_high = -1000;
    Minerva::BarID barId_z_low;
    Minerva::BarID barId_z_high;
    
    int story_low = 1000;
    int story_high = -1000;
    Minerva::BarID barId_t_low;
    Minerva::BarID barId_t_high;
    
    Gaudi::XYZPoint bar_pos;
    
    SmartRefVector<Minerva::ODCluster> odClusters = event->select<Minerva::ODCluster>( "All", "All" );
    if ( odClusters.empty() ) return StatusCode::FAILURE;
    
    SmartRefVector<Minerva::ODCluster>::iterator itOD = odClusters.begin();

    for ( ; itOD != odClusters.end(); itOD++ ){
        const Minerva::BarID barID = (*itOD)->barid();
        int frame = barID.frame();
        
        if (frame < frame_low) {
            frame_low = frame;
            barId_z_low = barID;
        }
        if (frame > frame_high) {
            frame_high = frame;
            barId_z_high = barID;
        }
        
        int story = barID.story();
        
        if (story < story_low) {
            story_low = story;
            barId_t_low = barID;
        }
        if (story > story_high) {
            story_high = story;
            barId_t_high = barID;
        }
        
        double energy = (*itOD)->energy();
        int tower = barID.tower();
        if (tower>=1 && tower<=6) {
            od_nclusters[tower-1]++;
            od_energy[tower-1] += energy;	
            od_time[tower-1] += (*itOD)->time()*(*itOD)->pe(); // od time charge weigthed 
            od_pe[tower-1] += (*itOD)->pe();
        } else {
            warning() << "invalid tower number: " << tower << endmsg;
        }

    }

    debug() << " Pass loop over ODClusters<oaltinok_version>" << endmsg;
    
    double maxODE = -1;
    int od_mostEnergeticTower = 0;
    std::vector<double> energyVector, timeVector, nclustersVector, distanceVector;
    std::vector<double> timeBlobODVector, timeBlobMuonVector, timeblobVector;

    for (int i=0; i<6; i++) {
        if ( od_pe[i] != 0 ) od_time[i] = od_time[i]/od_pe[i];
        energyVector.push_back(od_energy[i]);
        timeVector.push_back(od_time[i]);
        nclustersVector.push_back(od_nclusters[i]);
        if (od_energy[i]>maxODE) {
            maxODE = od_energy[i];
            od_mostEnergeticTower = i+1;
        }
    }

        //IDBlobs and ODTower
    Minerva::IDCluster::View view;
    if ( od_mostEnergeticTower == 2 || od_mostEnergeticTower == 5 ) view =  Minerva::IDCluster::X;
    else if ( od_mostEnergeticTower == 3 || od_mostEnergeticTower == 6 ) view =  Minerva::IDCluster::U;
    else if ( od_mostEnergeticTower == 1 || od_mostEnergeticTower == 4 ) view =  Minerva::IDCluster::V;
    else {
        warning() << " invalid tower number on IDBlobTower " << od_mostEnergeticTower << endmsg;
        return StatusCode::FAILURE;
    }
    
    std::vector<Minerva::IDBlob*>::iterator itBlob = idBlobs.begin();
    
    for ( ; itBlob != idBlobs.end(); itBlob++ ){
        
        double min_position = 10000, max_position = -10000;
        SmartRefVector<Minerva::IDCluster> idClusters =  (*itBlob)->clusters();
        SmartRefVector<Minerva::IDCluster>::iterator itClus = idClusters.begin();
        
        for ( ; itClus !=  idClusters.end(); itClus++ ){

            if ( (*itClus)->view() != view ) continue;
            if ( od_mostEnergeticTower <= 3 ) {
                if ( (*itClus)->position() < min_position )
                        min_position = (*itClus)->position();
            }
            else {
                if ( (*itClus)->position() > max_position )
                        max_position = (*itClus)->position();
            }
        }
        
        if ( od_mostEnergeticTower < 3 ) distanceVector.push_back( fabs(-1060-min_position) );
        else distanceVector.push_back( fabs(1060+max_position) );
        debug() << " ID Blob time " << (*itBlob)->time() << "; Tower time " << od_time[od_mostEnergeticTower-1]
                << " Muon time " << muon_position.T() << endmsg; 
        timeblobVector.push_back( (*itBlob)->time() );
        timeBlobODVector.push_back( (*itBlob)->time()-od_time[od_mostEnergeticTower-1] );
        timeBlobMuonVector.push_back( (*itBlob)->time()- muon_position.T());
    }
    
    event->setDoubleData( "od_upstreamFrame", frame_low );
    event->setDoubleData( "od_downstreamFrame", frame_high );
    event->setDoubleData( "od_upstreamFrame_z", m_odDet->getDeODFrame(barId_z_low)->getZCenter() );
    event->setDoubleData( "od_downstreamFrame_z", m_odDet->getDeODFrame(barId_z_high)->getZCenter() );

    event->setDoubleData( "od_highStory", story_high );
    event->setDoubleData( "od_lowStory", story_low );
    bar_pos = m_odDet->getDeODFrame(barId_t_high)->getBarPos(barId_t_high);
    event->setDoubleData( "od_highStory_t", sqrt( bar_pos.x()*bar_pos.x() + bar_pos.y()*bar_pos.y()) );
    bar_pos = m_odDet->getDeODFrame(barId_t_low)->getBarPos(barId_t_low);
    event->setDoubleData( "od_lowStory_t", sqrt( bar_pos.x()*bar_pos.x() + bar_pos.y()*bar_pos.y()) );
    
    event->setDoubleData( "od_maxEnergy", maxODE );
    event->setIntData( "od_energeticTower", od_mostEnergeticTower );
    
    event->setContainerDoubleData( "od_distanceBlobTower", distanceVector );
    event->setContainerDoubleData( "od_towerEnergy", energyVector );
    event->setContainerDoubleData( "od_towerNClusters", nclustersVector );
    event->setContainerDoubleData( "od_towerTime", timeVector );
    event->setContainerDoubleData( "od_idBlobTime", timeblobVector );
    event->setContainerDoubleData( "od_towerTimeBlobOD", timeBlobODVector );
    event->setContainerDoubleData( "od_towerTimeBlobMuon", timeBlobMuonVector );
    
    debug() << " Ending ODActivity::EnergeticTower " << endmsg;
    
    return StatusCode::SUCCESS;
}
    
    
//--------------------------------------------------------------------------------
// NukeCCQE Private Functions
//--------------------------------------------------------------------------------
    

//----------------------------------------------------------------------------------------
// interpret Events which fails the reconstructor cuts
//----------------------------------------------------------------------------------------
StatusCode CCDeltaPlusAna::interpretFailEvent( Minerva::PhysicsEvent* event ) const
{  
   debug() << "Enter CCDeltaPlusAna::interpretFailEvent" << endmsg;

   Minerva::NeutrinoInt* ccqeHyp = new Minerva::NeutrinoInt(m_anaSignature);
   NeutrinoVect nuInts;
   nuInts.push_back( ccqeHyp );
   markEvent(event);
   addInteractionHyp(event,nuInts);

   debug() << "Exit CCDeltaPlusAna::interpretFailEvent" << endmsg;
   return StatusCode::SUCCESS;
}

// //------------------------------------------------------------------------------------------
// // tag truth information - called independently from the above algorithms
// //------------------------------------------------------------------------------------------
// StatusCode CCDeltaPlusAna::tagTruth( Minerva::GenMinInteraction* truth) const
// {
//    verbose() << "Enter CCDeltaPlusAna::tagTruth" << endmsg;
//    if( !truth ) {
//      warning() << "The GenMinInteraction is NULL!" << endmsg;
//      return StatusCode::SUCCESS;
//    }
// 
//    fillGenieWeightBranches(truth);
//   
//    Minerva::PhysicsEvents* events = getPhysicsEvents();
//    for(Minerva::PhysicsEvents::iterator event = events->begin(); event != events->end(); event++) {
//      if( !(*event)->hasInteractionVertex() ) continue;
//      SmartRef<Minerva::Vertex> vertex = (*event)->interactionVertex();
//      if( !vertex ) continue;     
//      std::string name = "";
//      int tammy_passVertexZCut = 0;
//      if( !isInsideTargetFiducial(vertex,name,tammy_passVertexZCut) ) continue;     
//      if( tammy_passVertexZCut == 0 ) continue;
//      if( vertex->getOutgoingTracks().size() > 2 ) continue;
//      int ntracks = 0;
//      Minerva::ProngVect unattachedProngs = (*event)->select<Minerva::Prong>("Used:Unused","!IsSubProng");    
//      for(unsigned int i = 0; i < unattachedProngs.size(); i++) if( !unattachedProngs[i]->minervaTracks().empty() ) ntracks++; 
//      if( ntracks != 0 ) continue;
//      
//      if( truth->processType() == 1 && truth->current() == 1 && truth->incoming() == 14 ) {
//        int nprotons = 0;
//        const std::vector<int> fsParticlePDG = truth->fSpdg();
//        for(unsigned int i = 0; i < fsParticlePDG.size(); i++) { if( fsParticlePDG[i] == 2212 ) nprotons++; }
// 
//        if( nprotons != 0 ) truth->filtertaglist()->addFilterTag( "CCDeltaPlusAnaSignal", true );
//      }  
//    }
// 
//    //-- tag truth michel electrons
//    verbose() << "  Tagging truth michel electrons." << endmsg;
// 
//    Gaudi::LorentzVector positron_pos;
//    Gaudi::LorentzVector electron_pos;
// 
//    std::map<int,double> tammy_pion_plus_ids_momentum, tammy_pion_minus_ids_momentum;
//    std::map<int,int> tammy_muon_plus_id_parent, tammy_muon_minus_id_parent;
//    std::map<int,Gaudi::LorentzVector> electron_parent_pos, positron_parent_pos;
//    std::vector<double> tammy_pion_plus_momentum, tammy_pion_minus_momentum;
// 
//    Minerva::TG4Trajectories* alltrajects = get<Minerva::TG4Trajectories>( "MC/TG4Trajectories" );
//    for(Minerva::TG4Trajectories::iterator it = alltrajects->begin(); it != alltrajects->end(); ++it) {
//      Minerva::TG4Trajectory* traject = *it;
// 
//      //! look for charged pions
//      if( traject->GetPDGCode()==211 ){
//        tammy_pion_plus_ids_momentum[traject->GetTrackId()] = (traject->GetInitialMomentum()).P();
//      }
//      if( traject->GetPDGCode()==-211 ){
//        tammy_pion_minus_ids_momentum[traject->GetTrackId()] = (traject->GetInitialMomentum()).P();
//      }
// 
//      //! look for muons from decays
//      if( traject->GetProcessName() == "Decay" && traject->GetPDGCode() == -13 ) {
//        tammy_muon_plus_id_parent[traject->GetTrackId()] = traject->GetParentId();
//      }
// 
//      if( traject->GetProcessName() == "Decay" && traject->GetPDGCode() == 13 ) {
//        tammy_muon_minus_id_parent[traject->GetTrackId()] = traject->GetParentId();
//      }
// 
//      //! look for electrons from decays
//      if( traject->GetProcessName() == "Decay" && traject->GetPDGCode() == -11 ) {
//        positron_parent_pos[traject->GetParentId()] = traject->GetInitialPosition();
//      }
// 
//      if( traject->GetProcessName() == "muMinusCaptureAtRest" && traject->GetPDGCode() == 11 ) {
//        electron_parent_pos[traject->GetParentId()] = traject->GetInitialPosition();
//      }
// 
//    } // end of TG4Trajectories loop
// 
//    //! loop over electrons
//    for( std::map<int,Gaudi::LorentzVector>::iterator electron = electron_parent_pos.begin(); electron != electron_parent_pos.end(); ++electron ){
//     // check if electron parent is a muon-
//      if( tammy_muon_minus_id_parent.count( electron->first ) ){
//        //check if electron is inside fiducial region
//        Gaudi::LorentzVector electron_pos = electron->second;
//        Gaudi::XYZPoint p( electron_pos.x(), electron_pos.y(), electron_pos.z() );
//        if( FiducialPointTool->isFiducial( p, m_analyzeApothem, m_michel_upstreamZ, m_michel_downstreamZ ) ){
//          truth->setIntData("tammy_has_michel_electron",1);
//          // check if muon parent is a pion-
//          if( tammy_pion_minus_ids_momentum.count( tammy_muon_minus_id_parent[electron->first] ) ){
//            tammy_pion_minus_momentum.push_back( tammy_pion_minus_ids_momentum[ tammy_muon_minus_id_parent[electron->first] ]  );
//          }
//        }
//      }
//    }
// 
//    //! loop over positrons
//    for( std::map<int,Gaudi::LorentzVector>::iterator positron = positron_parent_pos.begin(); positron != positron_parent_pos.end(); ++positron ){
//     // check if positron parent is a muon+
//      if( tammy_muon_plus_id_parent.count( positron->first ) ){
//        //check if electron is inside fiducial region
//        Gaudi::LorentzVector positron_pos = positron->second;
//        Gaudi::XYZPoint p( positron_pos.x(), positron_pos.y(), positron_pos.z() );
//        if( FiducialPointTool->isFiducial( p, m_analyzeApothem, m_michel_upstreamZ, m_michel_downstreamZ ) ){
//          truth->setIntData("tammy_has_michel_electron",1);
//          // check if muon parent is a pion+ 
//          if( tammy_pion_plus_ids_momentum.count( tammy_muon_plus_id_parent[positron->first] ) ){
//            tammy_pion_plus_momentum.push_back( tammy_pion_plus_ids_momentum[ tammy_muon_plus_id_parent[positron->first] ]  );
//          }
//        }
//      }
//    }
// 
//    truth->setContainerDoubleData("tammy_has_michel_from_tammy_pion_plus_momentum", tammy_pion_plus_momentum);
//    truth->setContainerDoubleData("tammy_has_michel_from_tammy_pion_minus_momentum", tammy_pion_minus_momentum);
//    verbose() << "   End tagging truth michel electrons." << endmsg;
//    verbose() << "Exit CCDeltaPlusAna::tagTruth" << endmsg;
//    return StatusCode::SUCCESS;
// }

//--------------------------------------------------------------------------------------------------
// created odblob based prongs
//--------------------------------------------------------------------------------------------------
bool CCDeltaPlusAna::createdODBlobProng( Minerva::PhysicsEvent* event ) const
{
   verbose() << "Enter CCDeltaPlusAna::createdODBlobProng" << endmsg;
   bool make = false;

   //-- interaction vertex
   SmartRef<Minerva::Vertex> vertex = event->interactionVertex();

   //-- get all of the odblobs int the event
   Minerva::ODClusterVect odclusters = event->select<Minerva::ODCluster>("Unused","!LowActivity&!XTalkCandidate");
   verbose() << "  The number of unused now low activity ODClusters = " << odclusters.size() << endmsg;

   double tammy_odEnergy = 0, tammy_odEnergyWithTimeCut = 0;
   int tammy_n_odClustersWithTimeCut = 0;
   
   for(unsigned int clus = 0; clus < odclusters.size(); clus++) {
     tammy_odEnergy += odclusters[clus]->energy();
     double timeDifference = vertex->getTime() - odclusters[clus]->time();

     if( timeDifference >= m_clusMinTimeWindow && timeDifference <= m_clusMaxTimeWindow ) continue;
     tammy_odEnergyWithTimeCut += odclusters[clus]->energy();
     tammy_n_odClustersWithTimeCut++;
   }

   event->setIntData("tammy_n_odClusters",int(odclusters.size()));
   event->setIntData("tammy_n_odClustersWithTimeCut",tammy_n_odClustersWithTimeCut);
   event->setDoubleData("odEnergy",tammy_odEnergy);
   event->setDoubleData("tammy_odEnergyWithTimeCut",tammy_odEnergyWithTimeCut);

   verbose() << "Exit CCDeltaPlusAna::createdODBlobProng" << endmsg;
   return make;
}

//--------------------------------------------------------------------------------------------------
// created idblobs for non low activity id clusters
//--------------------------------------------------------------------------------------------------
bool CCDeltaPlusAna::createdIDBlobProng( Minerva::PhysicsEvent* event ) const
{
   verbose() << "Enter CCDeltaPlusAna::createdIDBlobProng" << endmsg;
   StatusCode sc;
   bool make = false;

   //-- get unused id clusters
   Minerva::IDClusterVect idclusters = event->select<Minerva::IDCluster>("Unused","!LowActivity&!XTalkCandidate");
   verbose() << "  The number of unused non low activity IDClusters = " << idclusters.size() << endmsg;

   //-- get vertices 
   Minerva::VertexVect vertices = event->select<Minerva::Vertex>("All");
   verbose() << "  The number of vertices in the event = " << vertices.size() << endmsg;

   //-- get all of the primary prongs
   Minerva::ProngVect prongs = event->primaryProngs();
   verbose() << "  The number of primary prongs in the event = " << prongs.size() << endmsg;

   //-- hide id clusters that do not fall in the the vertex time window
   SmartRef<Minerva::Vertex> vertex = event->interactionVertex();
   for(unsigned int clus = 0; clus < idclusters.size(); clus++) {
     double timeDifference = vertex->getTime() - idclusters[clus]->time();

     if( timeDifference >= m_clusMinTimeWindow && timeDifference <= m_clusMaxTimeWindow ) continue;
     idclusters[clus]->setHistory(Minerva::IDCluster::Hidden);
     idclusters[clus]->setIntData("isOutSideVertexTimeWindow",1);
   }

   //-- tag for vertex prongs
   bool pass = true;
   std::string tag = "PrimaryVertexEnergy";

   //-- get event manager
   Minerva::EventMgr* mgr = getEventMgr(vertex);

   //! make vertex prongs at all others vertices in the event
   for(unsigned int v = 0; v < vertices.size(); v++) {

     verbose() << "   The vertex type = " << vertices[v]->type() << endmsg;
     verbose() << "   The vertex position = " << vertices[v]->position() << endmsg;

     //-- create a container of vertex prongs
     std::vector<Minerva::Prong*>* vtxBlobProngVect = new std::vector<Minerva::Prong*>; 

     //-- make vertex prongs
     if( vertices[v]->type() == Minerva::Vertex::StartPoint ) {
       sc = m_primaryBlobProngTool->makeMultipleRadiiStyleVtxBlobProngs(mgr,idclusters,vertices[v],m_searchStepSize,m_numSearchRadii,
                                                                        m_maxSeparationBlobVertex,vtxBlobProngVect);  
     } else if( vertices[v]->type() == Minerva::Vertex::StopPoint ) {
       sc = makeConeVtxBlobProng(vertices[v],vtxBlobProngVect);

     } else {
       sc = m_primaryBlobProngTool->makeFilamentStyleVtxBlobProngs(mgr,idclusters,vertices[v],m_maxSearchDistance,
                                              m_maxStartingDistance,m_maxAllowedSearchGap,m_maxSeparationBlobVertex,vtxBlobProngVect);
     }

     verbose() << "  The number of created vtxBlobProngs = " << vtxBlobProngVect->size() << endmsg;

     //-- linked blob prongs the vertex and primary prongs
     if( !sc || vtxBlobProngVect->empty() ) { delete vtxBlobProngVect; continue; }
     make = true;

     std::vector<Minerva::Prong*>::iterator pit = vtxBlobProngVect->begin();
     for( ; pit != vtxBlobProngVect->end(); pit++) {
       if( vertices[v]->type() == Minerva::Vertex::StartPoint ) {
         event->promoteProngToPrimary(SmartRef<Minerva::Prong>(*pit));
         (*pit)->setSourceVertex(vertices[v]);
         (*pit)->setHistory(Minerva::Prong::Used);
         (*pit)->filtertaglist()->addFilterTag(tag,pass);
       } else {
         (*pit)->addLinkVertex(vertices[v]);
         (*pit)->setHistory(Minerva::Prong::Used);
         linkProngToPrimaryProng((*pit),prongs);
         if( vertices[v]->type() == Minerva::Vertex::Kinked ) {
           (*pit)->filtertaglist()->addFilterTag("SecondaryVertexEnergy",pass);
         } else if( vertices[v]->type() == Minerva::Vertex::StopPoint ) {
           (*pit)->filtertaglist()->addFilterTag("EndPointVertexEnergy",pass);
         }
       }
     }
   }

   //-- make a vector isolated blobs that should be associated with a the inner part of the primary prongs
   for(unsigned int p = 0; p < prongs.size(); p++) {
     verbose() << "   prong bit field = " << prongs[p]->typeBitsToString() << endmsg;

     //-- check if there are isolated blobs (3D) that should be associated with the inner part of the track
     std::vector<Minerva::IDBlob*>* trackBlobs = new std::vector<Minerva::IDBlob*>;
     sc = MuonUtils->findMuonEMBlobs(prongs[p],idclusters,m_maxSearchDistance,trackBlobs);
     verbose() << "  The number of created hidden track blobs = " << trackBlobs->size() << endmsg;

     //-- add blobs to prong and the store
     if( sc && !trackBlobs->empty() ) {
       std::vector<Minerva::IDBlob*>::iterator bit = trackBlobs->begin();
       for( ; bit != trackBlobs->end(); bit++) {
         addObject(event,(*bit));
         prongs[p]->add(SmartRef<Minerva::IDBlob>(*bit));
         (*bit)->setIntData("fuzzAttachedTrack",1);
         (*bit)->setHistory(Minerva::IDBlob::Used);
       }
     } else delete trackBlobs;

     //-- check if there are clusters that should be associated with the inner part of the track
     Minerva::IDBlob *fuzzBlob = new Minerva::IDBlob();
     sc = MuonUtils->findMuonFuzz(prongs[p],m_maxSearchDistance,m_fuzzRadius,fuzzBlob);
     verbose() << "  The number of clusters found in the fuzz blob = " << fuzzBlob->nclusters() << endmsg;

     if( sc && fuzzBlob->nclusters() ) {
       addObject(event,fuzzBlob);
       prongs[p]->add(SmartRef<Minerva::IDBlob>(fuzzBlob));
       fuzzBlob->setIntData("fuzzAttachedTrack",1);
       fuzzBlob->setHistory(Minerva::IDBlob::Used);
     } else delete fuzzBlob;
   }  

   //-- get remaining clusters after making vertex and fuzz blobs
   idclusters = event->select<Minerva::IDCluster>("Unused","!LowActivity&!XTalkCandidate");
   verbose() << "  The remaining number of clusters after making vertex and em blobs = " << idclusters.size() << endmsg;

   if( !idclusters.empty() ) {

     //-- store the remaining energy in a isolated prong
     Minerva::Prong* newProng = new Minerva::Prong("CCDeltaPlusAnaUnattachedProng");

     //-- add clusters to prong
     for(unsigned int clus = 0; clus < idclusters.size(); clus++) {
       idclusters[clus]->setHistory(Minerva::IDCluster::Used);
       newProng->add(idclusters[clus]);
     }

     //-- add prong to the store
     newProng->filtertaglist()->addFilterTag("UnattachedEnergy",pass);
     addObject(event,newProng);
     make = true;
   }

   verbose() << "Exit CCDeltaPlusAna::createdIDBlobProng" << endmsg;
   return make;
}

//---------------------------------------------------------------------------------------
// created particles for negative bit Minos prongs
//---------------------------------------------------------------------------------------
bool CCDeltaPlusAna::createdTrackedParticles( Minerva::ProngVect& prongs ) const
{
   verbose() << "Enter CCDeltaPlusAna::createdParticles" << endmsg;
   bool makeParticles = false;

   //-- check if the prongs are odMatch
   Minerva::EventMgr* mgr = getEventMgr(prongs[0]);
   Minerva::ODClusterVect odClusters = mgr->select<Minerva::ODCluster>("Unused","!LowActivity&!XTalkCandidate");
   m_odMatchTool->classifyProngs(prongs,odClusters);

   //-- loop over prongs
   for(unsigned int p = 0; p < prongs.size(); p++) {
     verbose() << "  The prong of bit-field = " << prongs[p]->typeBitsToString() << endmsg;

     //-- make sure the prong is not a bad object
     bool pass = true; std::string tag = "BadObject";
     if( prongs[p]->filtertaglist()->filterTagExists(tag) ) {
       if( prongs[p]->filtertaglist()->checkFilterTag(tag,pass) ) {
         error() << "This prong = " << prongs[p]->typeBitsToString() << " has been flag as a \"BadObject\", skipping!" << endmsg;
         continue;
       }
     }
     
     //-- skipped prongs with particles
     verbose() << "  The prong has n particles = " << prongs[p]->particles().size() << endmsg;
     if( !prongs[p]->particles().empty() ) continue;

     //! select the particle hypotheses candidates
     std::vector< Minerva::Particle::ID > hypotheses;
     hypotheses.push_back( Minerva::Particle::Muon );

     if( !prongs[p]->OdMatch() ) {
       hypotheses.push_back( Minerva::Particle::Proton );
       hypotheses.push_back( Minerva::Particle::Pion );
     }

     //! select the particle tools
     IParticleMakerTool::NameAliasListType toolsToUse;
     toolsToUse.push_back( std::make_pair("MuonMVAPIDTool","MuonMVAPIDTool") );
         
     if( !prongs[p]->OdMatch() ) toolsToUse.push_back( std::make_pair("dEdXTool","dEdXTool") );

     //-- make particles
     Minerva::Prong* prong = prongs[p];

     if( m_particleMaker->makeParticles(prong,hypotheses,toolsToUse).isSuccess() ) {
       verbose() << "  The prong of bit-field = " << prongs[p]->typeBitsToString() 
                 << " has " << prong->particles().size() << " number of particles attached." << endmsg; 
       makeParticles = true;
     } else verbose() << "  Did not make particles for the prong type = " << prongs[p]->typeBitsToString() << endmsg;
   }

   verbose() << "Exit CCDeltaPlusAna::createdParticles" << endmsg;
   return makeParticles;
}

//---------------------------------------------------------------------------------------
// classify or re-classify negative Minos-bit prongs
//---------------------------------------------------------------------------------------
bool CCDeltaPlusAna::classifyProngs( Minerva::PhysicsEvent* event ) const
{
   verbose() << "Enter CCDeltaPlusAna::classifyProngs" << endmsg;

   Minerva::ProngVect primaryProngs    = event->primaryProngs();
   Minerva::ProngVect unattachedProngs = event->select<Minerva::Prong>("Used:Unused");

   Minerva::ProngVect prongs;
   for(unsigned int p = 0; p < primaryProngs.size(); p++)    prongs.push_back( primaryProngs[p] );
   for(unsigned int p = 0; p < unattachedProngs.size(); p++) prongs.push_back( unattachedProngs[p] );

   verbose() << "  The number of prongs = " << prongs.size() << endmsg;
   bool classify = false;  

   //-- loop over prongs
   for(unsigned int p = 0; p < prongs.size(); p++) {

     verbose() << "  The prong bit = " << prongs[p]->typeBitsToString() 
               << " number of tracks = " << prongs[p]->minervaTracks().size() 
               << ", and the number of idblobs = " << prongs[p]->idblobs().size() << endmsg;

     //-- classify prongs
     if( prongs[p]->typeBitsToString().size() > 8 ) continue;
     if( m_prongIntersection->classifyProng(prongs[p]).isSuccess() ) {
       verbose() << "  The prong bit after tammy_classification = " << prongs[p]->typeBitsToString() << endmsg;
       classify = true;
     }
   }

   verbose() << "Exit CCDeltaPlusAna::classifyProngs" << endmsg;
   return classify;
}

//---------------------------------------------------------------------------------------
// rebuild prongs to include the new create tracks
//---------------------------------------------------------------------------------------
bool CCDeltaPlusAna::createdTrackedProngs( Minerva::PhysicsEvent* event ) const
{
   verbose() << "Enter CCDeltaPlusAna::createdProngs" << endmsg;

   //-- get interaction vertex and outgoing tracks
   SmartRef<Minerva::Vertex> vertex = event->interactionVertex();
   Minerva::TrackVect tracks = vertex->getOutgoingTracks();
   verbose() << "  The number of outgoing tracks = " << tracks.size() << endmsg;

   //-- get the created new primary track
   SmartRef<Minerva::Track> track;
   for(unsigned int trk = 0; trk < tracks.size(); trk++) {
     if( tracks[trk]->hasIntData("createdTrack") ) track = tracks[trk];
   }

   //-- check if the created track was found
   if( !track ) return false;

   //-- make new primary track prong 
   Minerva::Prong* newProng = new Minerva::Prong("CCDeltaPlusAnaPrimaryTrackProng");
   verbose() << "  Made a new Prong with creation signature : " << newProng->creationSignature() << endmsg;

   //-- add track and source vertex to the prong
   track->setHistory(Minerva::Track::Used);
   newProng->add(track);
   newProng->setSourceVertex(vertex);

   //-- add prong to the store and promote to a primary prong
   addObject(event,newProng);
   event->promoteProngToPrimary(SmartRef<Minerva::Prong>(newProng));

   verbose() << "Exit CCDeltaPlusAna::createdProngs" << endmsg;
   return true;
}

//---------------------------------------------------------------------------------------
// created kinked anchored tracks
//---------------------------------------------------------------------------------------
bool CCDeltaPlusAna::createdKinkedShortTracks( Minerva::PhysicsEvent* event ) const
{
   verbose() << "Enter CCDeltaPlusAna::createdKinkedShortTracks" << endmsg;
   bool createdKinked = false;
   StatusCode sc;

   //-- get all of the primary prongs in the event
   Minerva::ProngVect prongs = event->primaryProngs();
   verbose() << "  The number of primary prongs = " << prongs.size() << endmsg;

   //-- loop over primary prongs
   for(unsigned int p = 0; p < prongs.size(); p++) {

     //-- get the track
     SmartRef<Minerva::Track> track = prongs[p]->minervaTracks().back();

     //-- get stop point vertex
     Minerva::Vertex* vertex = getStopPointVertex(event,track);
     if( !vertex || vertex->type() != Minerva::Vertex::StopPoint ) continue;
     verbose() << "  The vertex's position = " << vertex->position() << endmsg;

     //-- continue making tracks
     bool can_make_kinked_track = true;
     while( can_make_kinked_track ) {

       //-- make short tracks
       bool make_short_tracks = createdAnchoredShortTracks(event,vertex); 
       if( !make_short_tracks ) { can_make_kinked_track = false; break; }

       //! set the flag
       createdKinked = true;

       //! add the tracks to the prong
       Minerva::TrackVect tracks = vertex->getOutgoingTracks();
       for(unsigned int t = 0; t < tracks.size(); ++t) {
         setTrackDirection(tracks[t],vertex);
         tracks[t]->setHistory(Minerva::Track::Used);
         prongs[p]->add(tracks[t]);
       }

       //! set prong and vertex bit types
       if( vertex->getNOutgoingTracks() == 1 ) {
         prongs[p]->setKinked(true);
         prongs[p]->addLinkVertex(vertex);
         vertex->setVertexFlags(Minerva::Vertex::Kinked);
       } else if( vertex->getNOutgoingTracks() > 1 ) {
         prongs[p]->setForked(true);
         prongs[p]->addLinkVertex(vertex);
         vertex->setVertexFlags(Minerva::Vertex::Forked);
         break;
       }

       //! get the track
       track = prongs[p]->minervaTracks().back();

       //! get vertex
       vertex = getStopPointVertex(event,track);
       if( !vertex ) break;

     } //-- end while loop over track making 
   } //-- end loop over primary prongs         

   verbose() << "Exit CCDeltaPlusAna::createdKinkedShortTracks" << endmsg;
   return createdKinked;
}

//---------------------------------------------------------------------------------------
// create short anchored tracks
//---------------------------------------------------------------------------------------
bool CCDeltaPlusAna::createdAnchoredShortTracks( Minerva::PhysicsEvent* event, Minerva::Vertex* vertex ) const
{
   verbose() << "Enter CCDeltaPlusAna::createdAnchoredShortTracks" << endmsg;
   StatusCode sc;
   bool createdTracks = false;

   //-- get the time slice
   Minerva::TimeSlice* tammy_timeSlice = getTimeSlice( vertex );
   short slice = tammy_timeSlice ? tammy_timeSlice->sliceNumber() : -1;
   verbose() << "   The slice associated with the event is = " << slice << endmsg;

   //-- anchored short track parameters
   double searchDist = 500.*CLHEP::mm;   //!-- maximum cluster search distance in a view  
   double startDist  = 500.*CLHEP::mm;   //!-- maximum cluster starting distance in a view  
   double gap        = 500.*CLHEP::mm;   //!-- maximum cluster allowed search gap in a view 
   double range3D    = 250.*CLHEP::mm;   //!-- maximum 3D distance between track endpoint and cluster 3D position
   double range2D    = 200.*CLHEP::mm;   //!-- maximum 2D distance between track endpoint and cluster 2D position
   bool   rawClus    = false;  //!-- used the default cluster list method for retreiving clusters
   bool   applyMS    = false;  //!-- apply the multiple scattering correction when fitting the track
   bool   force3D    = true;   //!-- enforced the x,y,z distance cuts

   //-- create a container to store tracks
   std::vector< Minerva::Track* > trackCont;
   trackCont.clear();

   //-- run the anchored short tracker
   bool can_make_anchor_tracks = true;
   while( can_make_anchor_tracks ) {
     Minerva::IDClusterVect idclusters = getClusters(event,Minerva::Track::FourHitPatRec);
     verbose() << "  The number of IDClusters = " << idclusters.size() << endmsg;

     if( idclusters.empty() ) { can_make_anchor_tracks = false; break; }

     //-- creates short tracks
     std::vector<Minerva::Track*> *tracks = new std::vector<Minerva::Track*>; 
     sc = m_anchoredTracker->createAnchoredTracks(idclusters,vertex,tracks,rawClus,searchDist,startDist,gap,
                                                  applyMS,force3D,range3D,range2D);
     verbose() << "  The number of anchor tracks created = " << tracks->size() << endmsg; 

     if( !sc || tracks->empty() ) { delete tracks; can_make_anchor_tracks = false; break; }
  
     //-- store tracks in the container
     for(std::vector<Minerva::Track*>::iterator trk = tracks->begin(); trk != tracks->end(); trk++) { 
       (*trk)->setPatRecHistory( Minerva::Track::FourHitPatRec );
       trackCont.push_back( (*trk) );
     }
   } //-- end loop for making anchored tracks 

   //-- run the vertex energy short tracker
   bool can_make_vertex_tracks = true;
   while( can_make_vertex_tracks ) {
     Minerva::IDClusterVect idclusters = getClusters(event,Minerva::Track::NotAssigned1);
     verbose() << "  The number of IDClusters = " << idclusters.size() << endmsg;
 
     if( idclusters.empty() ) { can_make_vertex_tracks = false; break; }
     
     std::vector<Minerva::Track*> *deShortTrackVector = new std::vector<Minerva::Track*>; 
     Gaudi::XYZTVector zDirectionVector(0, 0, 1, 0); 
     Gaudi::XYZPoint vertexPosition = vertex->position();
     sc = m_vertexEnergyStudyTool->reconstructVertexActivity(event,vertexPosition,zDirectionVector,idclusters,deShortTrackVector); 
     verbose() << "  The number of vertex short tracks created = " << deShortTrackVector->size() << endmsg;
 
     if( !sc || deShortTrackVector->empty() ) { delete deShortTrackVector; can_make_vertex_tracks = false; break; }

     //-- store tracks in the container
     for(std::vector<Minerva::Track*>::iterator trk = deShortTrackVector->begin(); trk != deShortTrackVector->end(); trk++) {
       (*trk)->setPatRecHistory( Minerva::Track::NotAssigned1 );
       trackCont.push_back( (*trk) );
     }
   } //-- end loop for making vertex short tracks

   //-- add tracks to store and vertex object
   for(std::vector<Minerva::Track*>::iterator track = trackCont.begin(); track != trackCont.end(); track++) {

     //-- add information to the track, vertex and event
     SmartRef<Minerva::Track> inTrack(*track);
     addObject(event,(*track));
     setTrackDirection(*track,vertex);
     (*track)->setTimeSlice(slice);
     (*track)->setIntData( "createdTrack",1 );
     vertex->addOutgoingTrack(inTrack);

     //-- create a stopPoint Vertex
     Minerva::Vertex* stopPointVertex = new Minerva::Vertex();
     stopPointVertex->setIncomingTrack(inTrack);
     stopPointVertex->setPosition((*track)->lastState().position());
     stopPointVertex->setTimeSlice(slice);
     stopPointVertex->setVertexFlags(Minerva::Vertex::StopPoint);
     addObject(event,stopPointVertex);
     createdTracks = true;
     verbose() << "  The new stop vertex's position = " << stopPointVertex->position() << endmsg;
   }

   //-- refit tracks
   Gaudi::XYZPoint  originalPosition = vertex->position();
   Gaudi::XYZVector originalPosError = vertex->positionErr();
   verbose() << "  The original vertex's position = " << vertex->position() << endmsg; 
    
   //-- refit vertex
   m_vertexFitter->fit(vertex); 

   std::vector<double> tammy_fit_vtx;
   tammy_fit_vtx.push_back( vertex->position().x() );
   tammy_fit_vtx.push_back( vertex->position().y() );
   tammy_fit_vtx.push_back( vertex->position().z() ); 

   event->setContainerDoubleData("fit_vtx",tammy_fit_vtx);

   SmartRef<Minerva::Vertex> vtx(vertex);
   event->setInteractionVertex( vtx );
   Minerva::ProngVect prongs = event->primaryProngs();
   for(unsigned int p = 0; p < prongs.size(); p++) prongs[p]->setSourceVertex(vertex);

   verbose() << "  The return vertex's position = " << vertex->position() << endmsg;   
   verbose() << "Exit CCDeltaPlusAna::createdAnchoredShortTracks" << endmsg;
   return createdTracks;
}

//---------------------------------------------------------------------------------------
// return the momentum analyzable contained ( proton candidate ) prong/particle
//---------------------------------------------------------------------------------------
bool CCDeltaPlusAna::protonProng( Minerva::ProngVect& prongs, SmartRef<Minerva::Prong>& prong, SmartRef<Minerva::Particle>& particle ) const
{
   verbose() << "Enter CCDeltaPlusAna::protonProng" << endmsg;
 
   bool pass = true; 
   std::string tag = "PrimaryProton";

   //-- get prong candidate
   Minerva::ProngVect tmpProngs;
   for(unsigned int p = 0; p < prongs.size(); p++) {
     if( prongs[p]->filtertaglist()->filterTagExists("PrimaryMuon") ) continue;
     tmpProngs.push_back( prongs[p] );
   }

   //-- find proton particle using ProtonUtils
   bool isProton = m_protonUtils->findProtonProng(tmpProngs,prong,particle);
   if( isProton ) { 
     prong->filtertaglist()->addFilterTag(tag,pass);
     prong->updateBestParticle(particle);
   }

   verbose() << "  found proton prong = " << isProton << endmsg;
   verbose() << "Exit CCDeltaPlusAna::protonProng" << endmsg;
   return isProton;
}

//----------------------------------------------------------------------------------------------
// return the momentum analyzable non-contained ( muon candidate ) prong/particle
//----------------------------------------------------------------------------------------------
bool CCDeltaPlusAna::muonProng( Minerva::ProngVect& prongs, SmartRef<Minerva::Prong>& prong, SmartRef<Minerva::Particle>& particle, 
                                  std::string& type ) const
{
   verbose() << "Enter CCDeltaPlusAna::muonProng" << endmsg;

   type = "unknown";
   bool pass = true;
   std::string tag = "PrimaryMuon";

   //-- containers
   Minerva::ProngVect odtmp, hcaltmp, ecaltmp, exittmp;

   //-- find muon particle 
   if( MuonUtils->findMuonProng(prongs,prong,particle) ) type = "golden"; 
   else {

     //-- loop over prongs 
     for(unsigned int p = 0; p < prongs.size(); p++) {
       if( prongs[p]->MinosTrack() ||  prongs[p]->MinosStub() ) { 
         type     = "silver"; 
         prong    = prongs[p]; 
         particle = findMuonParticleCandidate(prongs[p]); 
         break; 
       }
       else if( !prongs[p]->MinosTrack() && !prongs[p]->MinosStub() ) {
         if( prongs[p]->OdMatch() ) odtmp.push_back( prongs[p] );
         else if( !prongs[p]->OdMatch() && prongs[p]->DsHcal() )   hcaltmp.push_back( prongs[p] );
         else if( !prongs[p]->OdMatch() && prongs[p]->SideEcal() ) ecaltmp.push_back( prongs[p] );
       }
     }
    
     if( type == "unknown" ) {
       if( odtmp.size() == 1 ) { type = odMuonProng(odtmp[0],particle); prong = odtmp[0]; }
       else if( hcaltmp.size() == 1 ) { type = "silver"; prong = hcaltmp[0]; particle = findMuonParticleCandidate(hcaltmp[0]); }
       else if( ecaltmp.size() == 1 ) { type = "silver"; prong = ecaltmp[0]; particle = findMuonParticleCandidate(ecaltmp[0]); }
       else if( odtmp.size()   == 2 || hcaltmp.size() == 2 || ecaltmp.size() == 2 ) {
              
          if( !odtmp.empty() )        for(unsigned int p = 0; p < odtmp.size(); p++)   exittmp.push_back( odtmp[p] );
          else if( !hcaltmp.empty() ) for(unsigned int p = 0; p < hcaltmp.size(); p++) exittmp.push_back( hcaltmp[p] );
          else if( !ecaltmp.empty() ) for(unsigned int p = 0; p < ecaltmp.size(); p++) exittmp.push_back( ecaltmp[p] );

          Gaudi::XYZPoint trk1_pt1 = (exittmp[0]->minervaTracks().front())->firstState().position();
          Gaudi::XYZPoint trk2_pt1 = (exittmp[1]->minervaTracks().front())->firstState().position();

          Gaudi::XYZPoint trk1_pt2 = (exittmp[0]->minervaTracks().back())->lastState().position();
          Gaudi::XYZPoint trk2_pt2 = (exittmp[1]->minervaTracks().back())->lastState().position();

          double l1 = sqrt( pow(trk1_pt1.x()-trk1_pt2.x(),2) + pow(trk1_pt1.y()-trk1_pt2.y(),2) + pow(trk1_pt1.z()-trk1_pt2.z(),2) );
          double l2 = sqrt( pow(trk2_pt1.x()-trk2_pt2.x(),2) + pow(trk2_pt1.y()-trk2_pt2.y(),2) + pow(trk2_pt1.z()-trk2_pt2.z(),2) );

          if( l1 > l2 ) {
            if( !odtmp.empty() ) { type = odMuonProng(exittmp[0],particle); prong = exittmp[0]; }
            else { type = "silver"; prong = exittmp[0]; particle = findMuonParticleCandidate(exittmp[0]); }
          } else {
            if( !odtmp.empty() ) { type = odMuonProng(exittmp[1],particle); prong = exittmp[1]; }
            else { type = "silver"; prong = exittmp[1]; particle = findMuonParticleCandidate(exittmp[1]); }
          }
       }
     }
   } //-- if statement 

   //-- tag the prong
   bool isMuon = false;
   if( (type == "golden" || type == "silver") && particle ) { 
     prong->filtertaglist()->addFilterTag(tag,pass); 
     prong->updateBestParticle(particle);
     isMuon = true; 
   }

   verbose() << "  found muon prong = " << isMuon << " and the event type = " << type << endmsg;
   verbose() << "Exit CCDeltaPlusAna::muonProng" << endmsg;
   return isMuon;
}

//---------------------------------------------------------------------------------
// return the momentum analyzable and good particle id for odMatch prongs
//---------------------------------------------------------------------------------
std::string CCDeltaPlusAna::odMuonProng( SmartRef<Minerva::Prong>& prong, SmartRef<Minerva::Particle>& particle ) const
{
   verbose() << "Enter CCDeltaPlusAna::odMuonProng" << endmsg;
  
   //-- check if this is an upstream or downstream od prong 
   std::string type = "unknown"; 
   verbose() << "  Is OdMatchUS = " << prong->OdMatchUS() << ", OdMatchDS = " << prong->OdMatchDS() << endmsg;

   //-- loop over particles
   Minerva::ParticleVect particles = prong->particles();
   for(unsigned int part = 0; part < particles.size(); part++) {
     
     if( particles[part]->methodSignature().find("MuonMVAPID") == std::string::npos ) continue;

     verbose() << "  particle method of signature = " << particles[part]->methodSignature() 
               << ", score = " << particles[part]->score() << endmsg;
         
     if( particles[part]->idcode() == Minerva::Particle::Muon && particles[part]->score() > m_min_mva_score ) {
       type = "silver";
       particle = particles[part];  
     }   
   } 

   verbose() << "  the event type = " << type << endmsg;
   verbose() << "Exit CCDeltaPlusAna::odMuonProng" << endmsg;
   return type;
}

//--------------------------------------------------------------------------------
// targets fiducial volume
//--------------------------------------------------------------------------------
bool CCDeltaPlusAna::isInsideTargetFiducial( SmartRef<Minerva::Vertex>& vertex, std::string& tar_name, int& tammy_passVertexZCut ) const
{
   verbose() << "Enter CCDeltaPlusAna::isInsideTargetFiducial" << endmsg;

   tammy_passVertexZCut = 0;
   Gaudi::XYZPoint point   = vertex->position();
   bool insideTargetRegion = false; 

   for(unsigned int tar = 0; tar < m_targetNames.size(); tar++) {
     if( !m_coordSysTool->inFiducial(point,m_fvApothem,m_fvUpstreamZ[tar],m_fvDownstreamZ[tar]) ) continue; 
     tar_name = m_targetNames[tar]; 

     if( tar_name == "Target005" ) {
       if( m_coordSysTool->inFiducial(point,m_fvApothem,m_tar5VertexZUSCut,m_tar5VertexZDSCut) ) tammy_passVertexZCut = 1;
     } else if( tar_name == "Target004" ) {
       if( m_coordSysTool->inFiducial(point,m_fvApothem,m_tar4VertexZUSCut,m_tar4VertexZDSCut) ) tammy_passVertexZCut = 1;
     } else if( tar_name == "Target003" ) {
       if( m_coordSysTool->inFiducial(point,m_fvApothem,m_tar3VertexZUSCut,m_tar3VertexZDSCut) ) tammy_passVertexZCut = 1;
     } else if( tar_name == "Target002" ) {
       if( m_coordSysTool->inFiducial(point,m_fvApothem,m_tar2VertexZUSCut,m_tar2VertexZDSCut) ) tammy_passVertexZCut = 1;
     } else if( tar_name == "Target001" ) {
       if( m_coordSysTool->inFiducial(point,m_fvApothem,m_tar1VertexZUSCut,m_tar1VertexZDSCut) ) tammy_passVertexZCut = 1;
     }
     
     if( tar_name.find("NuclearTargets") != std::string::npos ) insideTargetRegion = true;
     else return true;
   }

   verbose() << "Exit CCDeltaPlusAna::isInsideTargetFiducial" << endmsg;
   if( insideTargetRegion ) return true;
   return false;
}

//------------------------------------------------------------------------------------
//  grab a list of clusters for creating a short track
//------------------------------------------------------------------------------------
Minerva::IDClusterVect CCDeltaPlusAna::getClusters( Minerva::PhysicsEvent* event, Minerva::Track::PatRecHistory pattern ) const
{
   verbose() << "Enter CCDeltaPlusAna::getClusters" << endmsg;
   Minerva::IDClusterVect retClusters; 
   StatusCode sc;
 
   //-- get interaction vertex
   SmartRef<Minerva::Vertex> vertex  = event->interactionVertex();   

   //-- get the clusters
   const char* clus_type = pattern == Minerva::Track::FourHitPatRec ? "Trackable:HeavyIonizing" : "!LowActivity&!XTalkCandidate";
   Minerva::IDClusterVect idClusters = event->select<Minerva::IDCluster>("Unused",clus_type);
   if( idClusters.size() > 200 ) return retClusters;
   verbose() << "  the number of idClusters before anchoring = " << idClusters.size() << endmsg;
   
   //-- create new anchored ID tracks to retrieve an axis for the cone
   std::vector<Minerva::Track*>* tracks = new std::vector<Minerva::Track*>;

   if( pattern == Minerva::Track::FourHitPatRec ) {
     sc = m_anchoredTracker->createAnchoredTracks(idClusters,vertex,tracks,true,6000,800,1000,false,false,50,850);
   } else {

     Minerva::IDClusterVect clustersNearVtxVector;
     Gaudi::XYZPoint vertexPosition = vertex->position();
     sc = m_blobCreatorUtils->fillProximateClusterVec(vertexPosition,idClusters,clustersNearVtxVector,m_maxSearchDistance_VESTool,
                                                      m_maxStartingDistance_VESTool,m_maxAllowedSearchGap_VESTool);
     
     if( !sc || clustersNearVtxVector.empty() ) { 
       Gaudi::XYZTVector zDirectionVector(0, 0, 1, 0);
       sc = m_vertexEnergyStudyTool->reconstructVertexActivity(event,vertexPosition,zDirectionVector,clustersNearVtxVector,tracks);
     }
   }

   //-- sanity checks for the number of tracks
   if( !sc || tracks->empty() ) { delete tracks; return retClusters; }
   verbose() << "   the number of tracks = " << tracks->size() << endmsg;

   //-- set track direction
   std::vector<Minerva::Track*>::iterator trk;
   for(trk = tracks->begin(); trk != tracks->end(); trk++) setTrackDirection(*trk,vertex);

   //-- retreive the track with the closest first state position to the vertex
   Minerva::Track* track = *tracks->begin();
   if( tracks->size() > 1 ) {
     double min_separation = 9999.9;
     for(trk = tracks->begin(); trk != tracks->end(); trk++) {
       double separation = MathTool->distanceBetween(vertex->position(),(*trk)->firstState().position());
       if( separation < min_separation ) { track = *trk; min_separation = separation; }
     }
   }

   //-- store the track's clusters 
   Minerva::IDClusterVect clusTmp;
   Minerva::Track::NodeContainer::iterator node;
   for(node = track->nodes().begin(); node != track->nodes().end(); node++) clusTmp.push_back( (*node)->idcluster() );

   //-- cone angle
   double coneAngle = 30*CLHEP::degree;

   //-- create cone
   Gaudi::XYZPoint firstState = track->firstState().position();
   double theta = track->theta();
   double phi   = track->phi();

   if( track->direction() == Minerva::Track::Backward ) {
     theta = Gaudi::Units::pi - theta;
     phi   = Gaudi::Units::twopi - phi;
   }

   Gaudi::Polar3DVector axis(1,theta,phi);
   Cone cone(firstState,axis,coneAngle);
   verbose() << "    cone's firstState position = " << firstState << ", axis = " << axis << endmsg;  

   //-- clean up tracks
   for( trk = tracks->begin(); trk != tracks->end(); trk++) { (*trk)->reset(); delete *trk;  }
   delete tracks;

   //-- get fill container with clusters inside of cone
   Minerva::IDClusterVect clusters = event->select<Minerva::IDCluster>("Unused","!LowActivity&!XTalkCandidate");
   verbose() << "   the number of clusters after anchoring = " << clusters.size() << endmsg;

   for(unsigned int clus = 0; clus < clusters.size(); clus++) {
     if( m_coneUtilsTool->isInsideCone(clusters[clus],cone) ) retClusters.push_back( clusters[clus] ); 
     else { 
       for(unsigned int tmp = 0; tmp < clusTmp.size(); tmp++) 
         if( clusTmp[tmp] == clusters[clus] ) retClusters.push_back( clusters[clus] );
     }
   }

   verbose() << "   the number of clusters return from coning = " << retClusters.size() << endmsg;
   verbose() << "Exit CCDeltaPlusAna::getClusters" << endmsg;
   return retClusters;
}

//---------------------------------------------------------------------------------
// get the stop vertex associated with the track
//---------------------------------------------------------------------------------
Minerva::Vertex* CCDeltaPlusAna::getStopPointVertex( Minerva::PhysicsEvent* event, SmartRef<Minerva::Track>& track ) const
{
   verbose() << "Enter CCDeltaPlusAna::getVertex" << endmsg;
   Minerva::Vertex* vertex = (Minerva::Vertex*)NULL;
   Minerva::VertexVect vertices = event->select<Minerva::Vertex>("Unused:Used","All");
   for(unsigned int v = 0; v < vertices.size(); v++) {
     if( vertices[v]->getNIncomingTracks() == 0 )   continue;
     if( vertices[v]->getIncomingTrack() == track ) vertex = vertices[v];
   }
   verbose() << "Exit CCDeltaPlusAna::getVertex" << endmsg;
   return vertex;
}

//-------------------------------------------------------------------------------
// corrected for the proton energy
//-------------------------------------------------------------------------------
void CCDeltaPlusAna::correctProtonProngEnergy( SmartRef<Minerva::Prong>& protonProng, double& p_calCorrection, 
                                                 double& p_visEnergyCorrection ) const
{
   verbose() << "Enter CCDeltaPlusAna::correctProtonProngEnergy" << endmsg;

   //-- initialize
   double calEnergy = 0.0;
   double visEnergy = 0.0;
    
   //-- get the proton prong subprongs
   Minerva::ProngVect subProngs = protonProng->subProngs();
   if( !subProngs.empty() ) {
     verbose() << "  The number of subProngs = " << subProngs.size() << endmsg;

     //-- loop over sub prongs
     for(unsigned int p = 0; p < subProngs.size(); p++) {

       //-- get the subprong link vertex
       Minerva::Vertex* vertex      = (Minerva::Vertex*)NULL;
       Minerva::VertexVect vertices = subProngs[p]->linkVertices();
       
       //-- skip subprongs attached to secondary vertex -- add back in energy at the dEdX stage 
       for(unsigned int v = 0; v < vertices.size(); ++v) {
         if( vertices[v]->type() == Minerva::Vertex::Kinked ) vertex = vertices[v];
       }
       if( vertex ) continue;

       //-- apply polyline calorimetric energy correction to the subProngs
       calEnergy += CaloUtils->applyCalConsts(subProngs[p]);
       verbose() << "  The calorimetric energy = " << calEnergy << " MeV" << endmsg;

       //-- get the total visible energy of the sub prongs
       Minerva::IDClusterVect clusters = subProngs[p]->getAllIDClusters();
       for(unsigned int clus = 0; clus < clusters.size(); clus++) visEnergy += clusters[clus]->energy();
       verbose() << "  The visible energy = " << visEnergy << " MeV" << endmsg;
     }
   }

   //-- get the proton prong fuzz energy
   Minerva::IDBlobVect fuzzBlobs = protonProng->getAllIDBlobs();
   for(unsigned int blob = 0; blob < fuzzBlobs.size(); blob++) {
     if( !fuzzBlobs[blob]->hasIntData("fuzzAttachedTrack") ) continue;
    
     //-- apply polyline calorimetric energy correction to the fuzz energy  
     calEnergy += CaloUtils->applyCalConsts(fuzzBlobs[blob]);

     //-- get the total visible energy of the fuzz blobs
     Minerva::IDClusterVect clusters = fuzzBlobs[blob]->clusters();
     for(unsigned int clus = 0; clus < clusters.size(); clus++) visEnergy += clusters[clus]->energy();
   }

   //-- get the prong first track
   SmartRef<Minerva::Track> track = (*protonProng->minervaTracks().begin());

   //-- update the particles four momentum
   Minerva::ParticleVect particles = protonProng->particles();
   for(unsigned int part = 0; part < particles.size(); part++) {
     if( particles[part]->idcode() != Minerva::Particle::Proton ) continue;

     Gaudi::XYZTVector fourMomentum = particles[part]->momentumVec();
     verbose() << "   particle idcode = " << particles[part]->idcode() << ", and four momomentum = " << fourMomentum << endmsg;

     double E_calCorrection = fourMomentum.E() + calEnergy;
     p_calCorrection = sqrt( pow(E_calCorrection,2) - pow(particles[part]->mass(),2) );
     verbose() << "  update energy using calorimetric correction = " << E_calCorrection << endmsg;

     double E_visEnergyCorrection = fourMomentum.E() + visEnergy;
     p_visEnergyCorrection = sqrt( pow(E_visEnergyCorrection,2) - pow(particles[part]->mass(),2) );     
     verbose() << " update energy using visible energy correction = " << E_visEnergyCorrection << endmsg;
   }     
        
   verbose() << "Exit CCDeltaPlusAna::correctProtonProngEnergy" << endmsg;
   return;
}


//--------------------------------------------------------------------------------
// make cone style vertex blobs 
//--------------------------------------------------------------------------------
StatusCode CCDeltaPlusAna::makeConeVtxBlobProng( SmartRef<Minerva::Vertex>& vertex, std::vector<Minerva::Prong*>* vtxBlobProngVect ) const
{
   verbose() << "Enter CCDeltaPlusAna::makeConeVtxBlobProng" << endmsg;
   Minerva::EventMgr* mgr = getEventMgr(vertex);

   //-- store all unused clusters in a container
   Minerva::IDClusterVect clusters = mgr->select<Minerva::IDCluster>("Unused","!LowActivity&!XTalkCandidate");
   verbose() << "  The number of remaining unused non low activity IDClusters = " << clusters.size() << endmsg;

   //! get the vertex incoming track
   SmartRef<Minerva::Track> track = vertex->getIncomingTrack();
   if( !track ) return StatusCode::SUCCESS;

   //-- cone angle
   double coneAngle = 30*CLHEP::degree;
   
   //-- create cone
   Gaudi::XYZPoint lastState = track->lastState().position();
   double theta = track->theta();
   double phi   = track->phi();

   if( track->direction() == Minerva::Track::Backward ) {
     theta = Gaudi::Units::pi - theta;
     phi   = Gaudi::Units::twopi - phi;
   }

   Gaudi::Polar3DVector axis(1,theta,phi);
   Cone cone(lastState,axis,coneAngle);
   verbose() << "    cone's last state position = " << lastState << ", axis = " << axis << endmsg;  

   //-- get clusters inside the cone
   Minerva::IDClusterVect coneClusters;
   for(unsigned int clus = 0; clus < clusters.size(); clus++) {
     if( m_coneUtilsTool->isInsideCone(clusters[clus],cone) ) coneClusters.push_back( clusters[clus] ); 
   }

   verbose() << "    the number of clusters inside cone = " << coneClusters.size() << endmsg;
   if( coneClusters.empty() ) return StatusCode::SUCCESS;

   //-- created a vertex blob and prong
   Minerva::IDBlob* newBlob = new Minerva::IDBlob( Minerva::IDBlob::UnknownSubdet );
   addObject(mgr,newBlob);
   newBlob->add(coneClusters);
   newBlob->setPatRecHistory(Minerva::IDBlob::NotAssigned2);
   newBlob->setHistory(Minerva::IDBlob::Used);

   Minerva::Prong* newProng = new Minerva::Prong("CCDeltaPlusAnaVtxProng");
   addObject(mgr,newProng);
   newProng->add( SmartRef<Minerva::IDBlob>(newBlob) );
   vertex->addIDBlob( SmartRef<Minerva::IDBlob>(newBlob) );

   //-- set clusters history
   for(unsigned int clus = 0; clus < coneClusters.size(); clus++) coneClusters[clus]->setHistory(Minerva::IDCluster::Used);

   //-- return prong to user
   vtxBlobProngVect->push_back( newProng );

   verbose() << "Exit CCDeltaPlusAna::makeConeVtxBlobProng" << endmsg;
   return StatusCode::SUCCESS;
}

//--------------------------------------------------------------------------------
// return the prong's muon particle
//--------------------------------------------------------------------------------
SmartRef<Minerva::Particle> CCDeltaPlusAna::findMuonParticleCandidate( SmartRef<Minerva::Prong>& prong ) const
{
   SmartRef<Minerva::Particle> muon;
   Minerva::ParticleVect particles = prong->particles();
   for(unsigned int part = 0; part < particles.size(); part++) {
     if( particles[part]->idcode() == Minerva::Particle::Muon ) muon = particles[part];
   }
   return muon;
}

//---------------------------------------------------------------------------------
// link a prong to a primary prong using the vertex
//---------------------------------------------------------------------------------
void CCDeltaPlusAna::linkProngToPrimaryProng( Minerva::Prong* prong, Minerva::ProngVect& primaryProngs ) const
{
   Minerva::VertexVect prongVertices = prong->linkVertices();
   for(unsigned int p = 0; p < primaryProngs.size(); p++) {
     Minerva::VertexVect primaryProngVertices = primaryProngs[p]->linkVertices();
     for(unsigned int vtx = 0; vtx < primaryProngVertices.size(); vtx++) {
       for(unsigned int v = 0; v < prongVertices.size(); v++) {
         if( prongVertices[v] == primaryProngVertices[vtx] ) primaryProngs[p]->add(SmartRef<Minerva::Prong>(prong));
       }
     }         
   }      
   return;
}

//----------------------------------------------------------------------------------
// set the created track direction
//----------------------------------------------------------------------------------
void CCDeltaPlusAna::setTrackDirection( Minerva::Track* track, Minerva::Vertex* vertex ) const
{
   int backward = 0, forward = 0;
   Minerva::Track::NodeContainer nodes = track->nodes();
   for(unsigned int node = 0; node < nodes.size(); node++) {
     if( nodes[node]->z() < vertex->position().z() )      backward++;
     else if( nodes[node]->z() > vertex->position().z() ) forward++;
   }
   if( backward > forward )      track->setDirection(Minerva::Track::Backward);
   else if( forward > backward ) track->setDirection(Minerva::Track::Forward);
   return;
}

//---------------------------------------------------------------------------------
// fill blob event data
//---------------------------------------------------------------------------------
void CCDeltaPlusAna::fillBlobEventData( Minerva::PhysicsEvent* event ) const 
{
    SmartRef<Minerva::Vertex> vertex    = event->interactionVertex();
    Minerva::ProngVect primaryProngs    = event->primaryProngs();
    Minerva::ProngVect unattachedProngs = event->select<Minerva::Prong>("Used:Unused","All");
    
    Minerva::ProngVect prongs;
    for(unsigned int p = 0; p < primaryProngs.size(); p++)    prongs.push_back( primaryProngs[p] );
    for(unsigned int p = 0; p < unattachedProngs.size(); p++) prongs.push_back( unattachedProngs[p] );
    
    double tammy_primaryVertexEnergy = 0, tammy_secondaryVertexEnergy = 0, tammy_endPointEnergy = 0, tammy_isolatedEnergy = 0;
    double tammy_isolatedEnergy_tracker = 0, tammy_isolatedEnergy_ecal = 0, tammy_isolatedEnergy_hcal = 0;
    double tammy_isolatedEnergy_targets = 0, tammy_isolatedEnergy_calorimetry = 0;
    double tammy_hadronic_energy = 0;
    
    for(unsigned int p = 0; p < prongs.size(); p++) {
        if( prongs[p]->filtertaglist()->filterTagExists("PrimaryMuon") ) continue;
        Minerva::IDClusterVect clusters = prongs[p]->getAllIDClusters();
        tammy_hadronic_energy += CaloUtils->applyCalConsts(clusters);
    } 
    
    for(unsigned int p = 0; p < prongs.size(); p++) {
        if( prongs[p]->filtertaglist()->filterTagExists("PrimaryVertexEnergy") ) {
            tammy_primaryVertexEnergy = prongs[p]->minervaVisibleEnergySum();
        } else if( prongs[p]->filtertaglist()->filterTagExists("SecondaryVertexEnergy") ) {
            tammy_secondaryVertexEnergy = prongs[p]->minervaVisibleEnergySum();
        } else if( prongs[p]->filtertaglist()->filterTagExists("EndPointVertexEnergy") ) { 
            tammy_endPointEnergy = prongs[p]->minervaVisibleEnergySum();
        } else if( prongs[p]->filtertaglist()->filterTagExists("UnattachedEnergy") ) {     
            tammy_isolatedEnergy = prongs[p]->minervaVisibleEnergySum();
            tammy_isolatedEnergy_calorimetry = CaloUtils->applyCalConsts(prongs[p]);
    
            Minerva::IDClusterVect clusters = prongs[p]->getAllIDClusters();
            for(unsigned int c = 0; c < clusters.size(); c++) {
                if( clusters[c]->subdet() == Minerva::IDCluster::NuclTargs ) {
                    tammy_isolatedEnergy_targets += clusters[c]->energy();       
                } else if( clusters[c]->subdet() == Minerva::IDCluster::Tracker ) {
                    tammy_isolatedEnergy_tracker += clusters[c]->energy();
                } else if( clusters[c]->subdet() == Minerva::IDCluster::ECAL ) {
                    tammy_isolatedEnergy_ecal += clusters[c]->energy();
                } else if( clusters[c]->subdet() == Minerva::IDCluster::HCAL ) {
                    tammy_isolatedEnergy_hcal += clusters[c]->energy();
                }
            }
        } 
    }
    
    double tammy_muonFuzzEnergy = 0, tammy_protonFuzzEnergy = 0;
    for(unsigned int p = 0; p < prongs.size(); p++) {
        Minerva::IDBlobVect blobs = prongs[p]->getAllIDBlobs();
        for(unsigned int b = 0; b < blobs.size(); b++) {
            if( !blobs[b]->hasIntData("fuzzAttachedTrack") )      continue;
            if( !blobs[b]->history() != Minerva::IDBlob::Hidden ) continue;
            if( prongs[p]->filtertaglist()->filterTagExists("PrimaryProton") )    tammy_protonFuzzEnergy += blobs[b]->energy();
            else if( prongs[p]->filtertaglist()->filterTagExists("PrimaryMuon") ) tammy_muonFuzzEnergy   += blobs[b]->energy();
        }
    }
    
    event->setDoubleData("muonFuzzEnergy",tammy_muonFuzzEnergy);
    event->setDoubleData("protonFuzzEnergy",tammy_protonFuzzEnergy);
    event->setDoubleData("primaryVertexEnergy",tammy_primaryVertexEnergy);
    event->setDoubleData("secondaryVertexEnergy",tammy_secondaryVertexEnergy);
    event->setDoubleData("endPointEnergy",tammy_endPointEnergy);
    event->setDoubleData("isolatedEnergy",tammy_isolatedEnergy);
    event->setDoubleData("tammy_isolatedEnergy_tracker",tammy_isolatedEnergy_tracker);   
    event->setDoubleData("tammy_isolatedEnergy_ecal",tammy_isolatedEnergy_ecal);
    event->setDoubleData("tammy_isolatedEnergy_hcal",tammy_isolatedEnergy_hcal);
    event->setDoubleData("tammy_isolatedEnergy_targets",tammy_isolatedEnergy_targets);
    event->setDoubleData("hadronic_energy",tammy_hadronic_energy);
    
    return;
}

//---------------------------------------------------------------------------------
// set the blobs Geant4 truth information
//---------------------------------------------------------------------------------
void CCDeltaPlusAna::setBlobProngTruth( Minerva::NeutrinoInt* neutrino, Minerva::ProngVect& prongs ) const
{
   if( prongs.empty() ) return;
   Minerva::TG4Trajectories* trajectories = get<Minerva::TG4Trajectories>( "MC/TG4Trajectories" ); 

   std::string name = "";
   if( prongs[0]->filtertaglist()->filterTagExists("PrimaryVertexEnergy") )        name = "tammy_primaryVtx";
   else if( prongs[0]->filtertaglist()->filterTagExists("SecondaryVertexEnergy") ) name = "tammy_secondaryVtx";
   else if( prongs[0]->filtertaglist()->filterTagExists("EndPointVertexEnergy") )  name = "tammy_endPointVtx";
   else if( prongs[0]->filtertaglist()->filterTagExists("UnattachedEnergy") )      name = "isolated";
   else if( prongs[0]->filtertaglist()->filterTagExists("PrimaryMuon") )           name = "muonFuzz";
   else if( prongs[0]->filtertaglist()->filterTagExists("PrimaryProton") )         name = "protonFuzz";
   else return;

   std::string pdgName        = name + "PDG";
   std::string primaryName    = name + "ParentId";
   std::string energyName     = name + "Energy";
   std::string trueEnergyName = name + "TrueEnergy";

   std::vector<int> pdg(200,-1);
   std::vector<int> primary(200,-1);
   std::vector<double> energy(200,-1);
   std::vector<double> trueEnergy(200,-1);

   int i = 0;
   for(unsigned int p = 0; p < prongs.size(); p++) {
     std::vector< Minerva::IDCluster* > clusters;

     if( name.find("secondary") != std::string::npos || name.find("endPoint") != std::string::npos ) {
       for(unsigned int b = 0; b < prongs[p]->getAllIDBlobs().size(); b++) {
         IDClusterVect tmp = prongs[p]->getAllIDBlobs()[b]->clusters();
         for(unsigned int c = 0; c < tmp.size(); c++) clusters.push_back( tmp[c] );
       }
     } else {
       for(unsigned int c = 0; c < prongs[p]->getAllIDClusters().size(); c++) clusters.push_back( prongs[p]->getAllIDClusters()[c] );
     }

     for(unsigned int clus = 0; clus < clusters.size(); clus++) {
       if( name.find("Fuzz") != std::string::npos && clusters[clus]->history() != Minerva::IDCluster::Hidden ) continue;

       std::vector<Minerva::IDCluster*> clusVect;
       clusVect.push_back( clusters[clus] );

       double other_energy = 0.0;
       std::map<const Minerva::TG4Trajectory*,double> trajMap = TruthMatcher->getTG4Trajectories(clusVect, other_energy);
       if( trajMap.empty() ) continue;
      
       const Minerva::TG4Trajectory* traj = NULL;
       double frac_energy = 0.0;
       std::map<const Minerva::TG4Trajectory*,double>::iterator tj = trajMap.begin();
       for( ; tj != trajMap.end(); tj++) {
         if( (*tj).second < frac_energy ) continue;
         frac_energy = (*tj).second;
         traj        = (*tj).first;
       }  
    
       if( i >= 200 ) continue; 

       pdg.at(i)        = traj->GetPDGCode();
       trueEnergy.at(i) = traj->GetInitialMomentum().E();
       energy.at(i)     = clusters[clus]->energy();

       if( traj->GetParentId() == 0 ) primary.at(i) = traj->GetPDGCode();
       else { 

           bool foundPrimary = false;
           while( !foundPrimary ) {

             bool foundMatch = false;
             Gaudi::XYZPoint start( traj->GetInitialPosition().x(), traj->GetInitialPosition().y(), traj->GetInitialPosition().z() );
 
             for(Minerva::TG4Trajectories::iterator it = trajectories->begin(); it != trajectories->end(); it++) {
               Gaudi::XYZPoint end( (*it)->GetFinalPosition().x(), (*it)->GetFinalPosition().y(), (*it)->GetFinalPosition().z() );
               if( end == start ) { traj = (*it); foundMatch = true; break; }
             }

             if( !foundMatch ) break;
             if( traj->GetParentId() == 0 ) { primary.at(i) = traj->GetPDGCode(); foundPrimary = true; break; }
           }           
       }  

       i++;
     } 
   }

   neutrino->setContainerIntData(pdgName,pdg);
   neutrino->setContainerIntData(primaryName,primary);
   neutrino->setContainerDoubleData(energyName,energy);
   neutrino->setContainerDoubleData(trueEnergyName,trueEnergy);
   return;
}

//---------------------------------------------------------------------------------
// set the track prong Geant4 truth information
//---------------------------------------------------------------------------------
void CCDeltaPlusAna::setTrackProngTruth( Minerva::NeutrinoInt* neutrino, SmartRef<Minerva::Prong>& prong ) const
{
   std::vector<const Minerva::TG4Trajectory*> parent;
   Minerva::TG4Trajectories* mc_trajectories = get<Minerva::TG4Trajectories>( "MC/TG4Trajectories" );
   for(Minerva::TG4Trajectories::iterator it = mc_trajectories->begin(); it != mc_trajectories->end(); ++it) {
     if( (*it)->GetParentId() == 0 ) parent.push_back( (*it) );
   }

   std::vector< std::pair<const Minerva::TG4Trajectory*,double> > trajectories;
   Minerva::TrackVect tracks = prong->minervaTracks();

   for(unsigned int trk = 0; trk < tracks.size(); trk++) {
     double other_energy = 0.0;
     std::map<const Minerva::TG4Trajectory*,double> trajMap = TruthMatcher->getTG4Trajectories(tracks[trk], other_energy);
     std::map<const Minerva::TG4Trajectory*,double>::iterator it = trajMap.begin();
     if( it == trajMap.end() ) continue;
     if( trajMap.empty() )     continue;

     std::vector< std::pair<const Minerva::TG4Trajectory*,double> > parentVect;
     for( ; it != trajMap.end(); it++) {

       const Minerva::TG4Trajectory* tmp = (*it).first;    
       double fraction = (*it).second;

       if( tmp->GetParentId() == 0 ) parentVect.push_back( std::make_pair(tmp,fraction) );
       else {

          int  attempts = 0;
          bool foundPrimary = false;
          while( !foundPrimary ) {
           
            bool foundMatch = false;
            Gaudi::XYZPoint start( tmp->GetInitialPosition().x(), tmp->GetInitialPosition().y(), tmp->GetInitialPosition().z() );

            for(Minerva::TG4Trajectories::iterator tj = mc_trajectories->begin(); tj != mc_trajectories->end(); ++tj) {
              Gaudi::XYZPoint end( (*tj)->GetFinalPosition().x(), (*tj)->GetFinalPosition().y(), (*tj)->GetFinalPosition().z() );
              if( end == start ) { tmp = (*tj); foundMatch = true; break; }
            }

            if( !foundMatch ) break;
            else attempts++;

            if( foundMatch && attempts == 10 ) break;
            if( tmp->GetParentId() == 0 ) { parentVect.push_back( std::make_pair(tmp,fraction) ); foundPrimary = true; break; }
          }
       }
     }

     std::vector< std::pair<const Minerva::TG4Trajectory*,double> > updateParentVect;
     for(unsigned int i = 0; i < parent.size(); i++) {
       double energy = 0.0;

       for(unsigned int j = 0; j < parentVect.size(); j++) {
         if( parentVect[j].first->GetTrackId() == parent[i]->GetTrackId() ) energy += parentVect[j].second; 
       }

       updateParentVect.push_back( std::make_pair(parent[i],energy) );
     }
       
     const Minerva::TG4Trajectory* tmp_traj = NULL;
     double frac_energy = 0.0; 
     for(unsigned int i = 0; i < updateParentVect.size(); i++) {
       if( updateParentVect[i].second < frac_energy ) continue;
       frac_energy = updateParentVect[i].second; 
       tmp_traj    = updateParentVect[i].first; 
     }

     trajectories.push_back( std::make_pair(tmp_traj,frac_energy) );
   }

   if( trajectories.empty() ) return;
   const Minerva::TG4Trajectory* traj = NULL;
   double max = 0;
   for(unsigned int i = 0; i < trajectories.size(); i++) {
     if( trajectories[i].second < max ) continue;
     max  = trajectories[i].second;
     traj = trajectories[i].first;
   }

   if( !traj ) return;

   Gaudi::LorentzVector traj_4p = traj->GetInitialMomentum();
   double tj_p     = sqrt( traj_4p.px()*traj_4p.px() + traj_4p.py()*traj_4p.py() + traj_4p.pz()*traj_4p.pz() );
   double tj_theta = traj->GetInitialMomentum().Theta();
   double tj_phi   = traj->GetInitialMomentum().Phi();

   int primary = 0;
   if( traj->GetProcessName().find("Primary") != std::string::npos ) primary = 1; 

   if( prong->filtertaglist()->filterTagExists("PrimaryProton") ) {
     neutrino->setIntData("ntrajProtonProng",int(trajectories.size()));
     neutrino->setIntData("tammy_trajProtonProngPrimary",primary);
     neutrino->setIntData("tammy_trajProtonProngPDG",traj->GetPDGCode());
     neutrino->setDoubleData("tammy_trajProtonProngMomentum",tj_p);
     neutrino->setDoubleData("tammy_trajProtonTheta",tj_theta);
     neutrino->setDoubleData("tammy_trajProtonPhi",tj_phi);
   } else if( prong->filtertaglist()->filterTagExists("PrimaryMuon") ) {
     neutrino->setIntData("tammy_ntrajMuonProng",int(trajectories.size()));
     neutrino->setIntData("tammy_trajMuonProngPrimary",primary);
     neutrino->setIntData("tammy_trajMuonProngPDG",traj->GetPDGCode());
     neutrino->setDoubleData("tammy_trajMuonProngMomentum",tj_p);
     neutrino->setDoubleData("tammy_trajMuonTheta",tj_theta);
     neutrino->setDoubleData("tammy_trajMuonPhi",tj_phi);
   }

   const Minerva::TG4Trajectory* endTraj = trajectories[trajectories.size()-1].first;
   if( !endTraj ) return;

   const Minerva::TG4TrajectoryPoints& points = endTraj->GetTrajectoryPoints();
   
   if( !points.empty() ) {
     Minerva::TG4TrajectoryPoint* endPoint = *points.rbegin();
     if( endPoint ) {
       double trajPx = endPoint->GetMomentum().px();
       double trajPy = endPoint->GetMomentum().py();
       double trajPz = endPoint->GetMomentum().pz();
       double trajMomentum  = sqrt( trajPx*trajPx + trajPy*trajPy + trajPz*trajPz );

       double trajX = endPoint->GetPosition().x();
       double trajY = endPoint->GetPosition().y();
       double trajZ = endPoint->GetPosition().z();
       Gaudi::XYZPoint point(trajX,trajY,trajZ);

       if( prong->filtertaglist()->filterTagExists("PrimaryProton") ) {
         neutrino->setIntData("tammy_isProtonInsideOD",int(m_odDet->isInside(point)));
         neutrino->setDoubleData("tammy_tammy_endProtonTrajMomentum",trajMomentum);
         neutrino->setDoubleData("tammy_endProtonTrajXPosition",trajX);
         neutrino->setDoubleData("tammy_endProtonTrajYPosition",trajY);
         neutrino->setDoubleData("tammy_endProtonTrajZPosition",trajZ);
       } else if( prong->filtertaglist()->filterTagExists("PrimaryMuon") ) {
         neutrino->setIntData("tammy_isMuonInsideOD",int(m_odDet->isInside(point)));
         neutrino->setDoubleData("tammy_tammy_endMuonTrajMomentum",trajMomentum);
         neutrino->setDoubleData("tammy_endMuonTrajXPosition",trajX);
         neutrino->setDoubleData("tammy_endMuonTrajYPosition",trajY);
         neutrino->setDoubleData("tammy_endMuonTrajZPosition",trajZ);
       }
     } 
   }

   return;
}

//----------------------------------------------------------------------------------
// set the intranuke and genie truth information
//----------------------------------------------------------------------------------
void CCDeltaPlusAna::setGenMinTruthInformation( Minerva::PhysicsEvent* event, Minerva::GenMinInteraction* truth ) const
{
    verbose() << "Enter CCDeltaPlusAna::setGenMinTruthInformation" << endmsg;
    
    int tammy_genie_n_neutrinos     = 0;
    int tammy_genie_n_muons         = 0;
    int tammy_genie_n_mesons        = 0;
    int tammy_genie_n_heavy_baryons = 0;
    int tammy_genie_n_photons       = 0;
    int tammy_genie_n_protons       = 0;
    int tammy_genie_n_neutrons      = 0;
    int tammy_genie_n_pions         = 0;
    int tammy_genie_n_pi_zeros      = 0;
    int tammy_genie_n_charms        = 0;
    int tammy_genie_n_kaons         = 0;
    int tammy_genie_n_others        = 0;
    int tammy_genie_n_particles     = 0;
    
    std::vector< std::pair<int,double> > nucleonVect;
    
    std::vector<double> pxVec = truth->fsParticlesPx();
    std::vector<double> pyVec = truth->fsParticlesPy();
    std::vector<double> pzVec = truth->fsParticlesPz();
    std::vector<int>    pdg   = truth->fSpdg();
    
    tammy_genie_n_particles = (int)pdg.size();
    
    for(unsigned int p = 0; p < pdg.size(); p++) {
        double ptot = sqrt( pxVec[p]*pxVec[p] + pyVec[p]*pyVec[p] + pzVec[p]*pzVec[p] );
    
        if( fabs( pdg[p] ) == 13 ) tammy_genie_n_muons++;
        else if( fabs( pdg[p] ) == 14 || fabs( pdg[p] ) == 15 || fabs( pdg[p] ) == 16 ) tammy_genie_n_neutrinos++;
        else if( fabs( pdg[p] ) == 22 )   tammy_genie_n_photons++;
        else if( fabs( pdg[p] ) == 2212 ) { tammy_genie_n_protons++;  nucleonVect.push_back( std::make_pair(pdg[p],ptot) ); }
        else if( fabs( pdg[p] ) == 2112 ) { tammy_genie_n_neutrons++; nucleonVect.push_back( std::make_pair(pdg[p],ptot) ); }
        else if( fabs( pdg[p] ) == 211 )  { tammy_genie_n_pions++;    tammy_genie_n_mesons++; }
        else if( fabs( pdg[p] ) == 111 )  { tammy_genie_n_pi_zeros++; tammy_genie_n_mesons++; }
        else if( fabs( pdg[p] ) == 411 || fabs( pdg[p] ) == 421 || fabs( pdg[p] ) == 431 ) { tammy_genie_n_charms++; tammy_genie_n_mesons++; }
        else if( fabs( pdg[p] ) == 321 || fabs( pdg[p] ) == 311 || fabs( pdg[p] ) == 310 || fabs( pdg[p] ) == 130 ) {
        tammy_genie_n_kaons++; tammy_genie_n_mesons++;
        }
        else if( fabs( pdg[p] ) > 3000 && fabs( pdg[p] ) < 5000 ) tammy_genie_n_heavy_baryons++;
        else tammy_genie_n_others++;
    }
    
    event->setIntData("tammy_genie_n_neutrinos",tammy_genie_n_neutrinos);
    event->setIntData("tammy_genie_n_muons",tammy_genie_n_muons);
    event->setIntData("tammy_genie_n_mesons",tammy_genie_n_mesons);
    event->setIntData("tammy_genie_n_heavy_baryons",tammy_genie_n_heavy_baryons);
    event->setIntData("tammy_genie_n_photons",tammy_genie_n_photons);
    event->setIntData("tammy_genie_n_protons",tammy_genie_n_protons);
    event->setIntData("tammy_genie_n_neutrons",tammy_genie_n_neutrons);
    event->setIntData("tammy_genie_n_pions",tammy_genie_n_pions);
    event->setIntData("tammy_genie_n_pi_zeros",tammy_genie_n_pi_zeros);
    event->setIntData("tammy_genie_n_charms",tammy_genie_n_charms);
    event->setIntData("tammy_genie_n_kaons",tammy_genie_n_kaons);
    event->setIntData("tammy_genie_n_others",tammy_genie_n_others);
    event->setIntData("tammy_genie_n_particles",tammy_genie_n_particles);
    
    bool isCCQElike = false;
    if( truth->current()         == Minerva::GenMinInteraction::kChargedCurrent &&
        truth->interactionType() == Minerva::GenMinInteraction::kQEL && truth->incoming() ==  14 ){
        
            isCCQElike = true;
        }
    else if( tammy_genie_n_muons == 1 && tammy_genie_n_protons != 0 && tammy_genie_n_mesons == 0 &&
                tammy_genie_n_heavy_baryons == 0 && tammy_genie_n_photons == 0 ){ 
                
            isCCQElike = true;
        }
    
    double momentum = -1;
    std::vector< double > momentumVec(4,-1);  
    
    int index   = -1;
    int scatter =  0;
    int delta   =  0; 
    int other   =  0;
    
    Minerva::GenMinEventRecord* eventRecord = truth->eventRecord();
    int n_fs_parts = eventRecord->eRecNParticleGen();
    
    double leading_tammy_proton_p = 0;
    for(unsigned int i = 0; i < nucleonVect.size(); i++) {
        if( nucleonVect[i].first != 2212 ) continue;
        if( nucleonVect[i].second > leading_tammy_proton_p ) leading_tammy_proton_p = nucleonVect[i].second;
    }
    
    if( isCCQElike ) {
        std::vector< std::pair<int,int> > list;
        list.clear();
    
        for(int i = 0; i < n_fs_parts; i++) {
        if( eventRecord->eRecPartID(i) != 2212 ) continue;
        
        double px = eventRecord->eRecMomentumPx(i);
        double py = eventRecord->eRecMomentumPy(i);
        double pz = eventRecord->eRecMomentumPz(i);
        double p  = sqrt( px*px + py*py + pz*pz );
    
        if( p != leading_tammy_proton_p ) continue;
        list.push_back( std::make_pair(i,eventRecord->eRecPartID(i)) );
    
        int mother_index = eventRecord->eRecMother(i);
        bool foundMother = false;
    
        while( !foundMother ) {
            if( mother_index == -1 ) foundMother = true;
            else {
                list.push_back( std::make_pair(mother_index,eventRecord->eRecPartID(mother_index)) );
                mother_index = eventRecord->eRecMother( mother_index );
            }
        }
    
        break;
        }
        
        std::vector< std::pair<int,int> >::reverse_iterator rit;
        for(rit = list.rbegin(); rit != list.rend(); rit++) {
        if( eventRecord->eRecPartStatus( (*rit).first ) == 0 ) continue;       
        
        if( eventRecord->eRecPartID( (*rit).first )      == 2112 ) scatter = 1;
        else if( eventRecord->eRecPartID( (*rit).first ) == 2224 ) delta   = 1;
    
        if( scatter == 1 || delta == 1 ) {
            rit++;
            if( (*rit).second == 2212 ) index = (*rit).first;
            else { scatter = 0; delta = 0; }
        }
        
        break;
        }
    }
    
    if( index != -1 ) {
        double px = eventRecord->eRecMomentumPx(index);
        double py = eventRecord->eRecMomentumPy(index);
        double pz = eventRecord->eRecMomentumPz(index);
        double p  = sqrt( px*px + py*py + pz*pz );
        double e  = eventRecord->eRecEnergy(index);
    
        momentumVec.at(0) = px;
        momentumVec.at(1) = py;
        momentumVec.at(2) = pz;
        momentumVec.at(3) = e;
        momentum = p;
    } else other = 1;  
    
    event->setDoubleData("tammy_intraNukeProtonMomentum",momentum);
    event->setContainerDoubleData("tammy_intraNukeProtonMomentumVec",momentumVec);
    event->setIntData("tammy_intraNukeNParticles",n_fs_parts);
    event->setIntData("tammy_intraNukeDeltaPlusPlusDecay",delta);
    event->setIntData("tammy_intraNukeNeutronQuasiElasticScatter",scatter);
    event->setIntData("tammy_intraNukeOtherProcess",other);
    
    verbose() << "Exit CCDeltaPlusAna::setGenMinTruthInformation" << endmsg;
    return;
}

//---------------------------------------------------------------------------------
// set the muon id
//---------------------------------------------------------------------------------
void CCDeltaPlusAna::setMuonID( Minerva::NeutrinoInt* neutrino, SmartRef<Minerva::Prong>& prong ) const
{
   if( prong->MinosTrack() )     neutrino->setIntData( "tammy_muon_minos_track", 1 );
   else if( prong->MinosStub() ) neutrino->setIntData( "tammy_muon_minos_stub", 1 );
   else if( !prong->MinosTrack() && !prong->MinosStub() ) {
     if( prong->OdMatch() ) neutrino->setIntData( "tammy_muon_od_track", 1 );
     else if( !prong->OdMatch() ) {
       if( prong->SideEcal() )    neutrino->setIntData( "tammy_muon_side_ecal", 1 );
       else if( prong->DsHcal() ) neutrino->setIntData( "tammy_muon_down_hcal", 1 );
     } 
   }
   return;
}

//----------------------------------------------------------------------------------
// set proton particle data
//----------------------------------------------------------------------------------
void CCDeltaPlusAna::setProtonParticleData( Minerva::NeutrinoInt* neutrino, SmartRef<Minerva::Prong>& prong,
                                              SmartRef<Minerva::Particle>& particle, double vertexZ ) const 
{
   double tammy_proton_theta  = m_coordSysTool->thetaWRTBeam(particle->momentumVec());
   double tammy_proton_thetaX = m_coordSysTool->thetaXWRTBeam(particle->momentumVec());
   double tammy_proton_thetaY = m_coordSysTool->thetaYWRTBeam(particle->momentumVec());
   double tammy_proton_phi    = m_coordSysTool->phiWRTBeam(particle->momentumVec());

   Gaudi::LorentzVector protonfourVec;
   StatusCode sc = m_energyCorrectionTool->getCorrectedEnergy(prong,particle,vertexZ,protonfourVec);
   if( !sc ) protonfourVec = particle->momentumVec();
   else particle->setMomentumVec(protonfourVec);

   double p_calCorrection = -1, p_visEnergyCorrection = -1;
   correctProtonProngEnergy(prong,p_calCorrection,p_visEnergyCorrection);

   neutrino->setDoubleData("tammy_proton_p_calCorrection",p_calCorrection);
   neutrino->setDoubleData("tammy_proton_p_visEnergy",p_visEnergyCorrection);

   double T_proton   = protonfourVec.E() - MinervaUnits::M_p;
   double enu_proton = m_protonUtils->nuEnergyCCQE(tammy_proton_theta,T_proton,m_CCQEBindingEnergyMeV);
   double Q2_proton  = m_protonUtils->qSquaredCCQE(T_proton,m_CCQEBindingEnergyMeV);

   std::vector<double> proton4p;
   proton4p.push_back( protonfourVec.px() );
   proton4p.push_back( protonfourVec.py() );
   proton4p.push_back( protonfourVec.pz() );
   proton4p.push_back( protonfourVec.E() );

   std::vector<double> tammy_proton_startPoint, tammy_proton_endPoint;
   tammy_proton_endPoint.push_back( (prong->minervaTracks().back())->lastState().x() );
   tammy_proton_endPoint.push_back( (prong->minervaTracks().back())->lastState().y() );
   tammy_proton_endPoint.push_back( (prong->minervaTracks().back())->lastState().z() );

   tammy_proton_startPoint.push_back( (prong->minervaTracks().front())->firstState().x() );
   tammy_proton_startPoint.push_back( (prong->minervaTracks().front())->firstState().y() );
   tammy_proton_startPoint.push_back( (prong->minervaTracks().front())->firstState().z() );

   neutrino->setIntData("tammy_proton_kinked",(int)prong->Kinked());
   neutrino->setIntData("tammy_proton_odMatch",(int)prong->OdMatch());
   neutrino->setDoubleData("tammy_proton_score",particle->score());
   neutrino->setDoubleData("tammy_proton_score1",particle->getDoubleData("score1"));
   neutrino->setDoubleData("tammy_proton_score2",particle->getDoubleData("score2"));
   neutrino->setDoubleData("tammy_proton_chi2_ndf",particle->getDoubleData("chi2_ndf"));
   neutrino->setDoubleData("tammy_proton_theta",tammy_proton_theta);
   neutrino->setDoubleData("tammy_proton_thetaX",tammy_proton_thetaX);
   neutrino->setDoubleData("tammy_proton_thetaY",tammy_proton_thetaY);
   neutrino->setDoubleData("tammy_proton_phi",tammy_proton_phi);
   neutrino->setDoubleData("tammy_proton_enu",enu_proton);
   neutrino->setDoubleData("tammy_proton_q2",Q2_proton);
   neutrino->setContainerDoubleData("tammy_proton_endPoint",tammy_proton_endPoint);
   neutrino->setContainerDoubleData("tammy_proton_startPoint",tammy_proton_startPoint);
   neutrino->setContainerDoubleData("tammy_proton_4p",proton4p);

   return;
}

//-----------------------------------------------------------------------------------
// set some pion variables
//-----------------------------------------------------------------------------------
void CCDeltaPlusAna::setPionParticleData( Minerva::NeutrinoInt* neutrino, SmartRef<Minerva::Prong>& prong ) const
{
   Minerva::ParticleVect particles = prong->particles();
   for(unsigned int p = 0; p < particles.size(); p++) {
     if( particles[p]->idcode() == Minerva::Particle::Pion && !particles[p]->isMultiMass() ) {
       neutrino->setDoubleData("tammy_pion_score",particles[p]->score());
       neutrino->setDoubleData("tammy_pion_score1",particles[p]->getDoubleData(ParticleExtraDataDefs::dEdXScore1()));
       neutrino->setDoubleData("tammy_pion_score2",particles[p]->getDoubleData(ParticleExtraDataDefs::dEdXScore2()));
       neutrino->setDoubleData("tammy_pion_chi2_ndf",particles[p]->getDoubleData(ParticleExtraDataDefs::dEdXReducedChiSquared()));
     }
   }
   return;
}

//--------------------------------------------------------------------------------
// set od muon particle data
//--------------------------------------------------------------------------------
void CCDeltaPlusAna::setODMuonParticleData( Minerva::NeutrinoInt* nu, SmartRef<Minerva::Particle>& particle ) const
{
   if( particle->hasDoubleData(ParticleExtraDataDefs::ODFaceXPosition()) ) 
     nu->setDoubleData( "tammy_muon_odFaceX",particle->getDoubleData(ParticleExtraDataDefs::ODFaceXPosition()) );

   if( particle->hasDoubleData(ParticleExtraDataDefs::ODFaceYPosition()) )  
     nu->setDoubleData( "tammy_muon_odFaceY",particle->getDoubleData(ParticleExtraDataDefs::ODFaceYPosition()) );

   if( particle->hasDoubleData(ParticleExtraDataDefs::ODFaceZPosition()) ) 
     nu->setDoubleData( "tammy_muon_odFaceZ",particle->getDoubleData(ParticleExtraDataDefs::ODFaceZPosition()) );

   if( particle->hasDoubleData(ParticleExtraDataDefs::ODEndXPosition()) ) 
     nu->setDoubleData( "tammy_muon_odEndX",particle->getDoubleData(ParticleExtraDataDefs::ODEndXPosition()) );

   if( particle->hasDoubleData(ParticleExtraDataDefs::ODEndYPosition()) ) 
     nu->setDoubleData( "tammy_muon_odEndY",particle->getDoubleData(ParticleExtraDataDefs::ODEndYPosition()) );

   if( particle->hasDoubleData(ParticleExtraDataDefs::ODEndZPosition()) ) 
     nu->setDoubleData( "tammy_muon_odEndZ",particle->getDoubleData(ParticleExtraDataDefs::ODEndZPosition()) );

   if( particle->hasDoubleData(ParticleExtraDataDefs::ODLastClusZ()) )
     nu->setDoubleData( "tammy_muon_odLastClusZ",particle->getDoubleData(ParticleExtraDataDefs::ODLastClusZ()) );

   if( particle->hasIntData(ParticleExtraDataDefs::ODLastStory()) )
     nu->setIntData( "tammy_muon_odLastStory",particle->getIntData(ParticleExtraDataDefs::ODLastStory()) );

   if( particle->hasIntData(ParticleExtraDataDefs::ODLastFrame()) )
     nu->setIntData( "tammy_muon_odLastFrame",particle->getIntData(ParticleExtraDataDefs::ODLastFrame()) );

   if( particle->hasDoubleData(ParticleExtraDataDefs::ODStopDistMomentum()) )
     nu->setDoubleData( "tammy_muon_odStopDistMomentum",particle->getDoubleData(ParticleExtraDataDefs::ODStopDistMomentum()) );

   if( particle->hasDoubleData(ParticleExtraDataDefs::ODdEdXMomentum()) )
     nu->setDoubleData( "tammy_muon_odElossMomentum",particle->getDoubleData(ParticleExtraDataDefs::ODdEdXMomentum()) );

   if( particle->hasDoubleData(ParticleExtraDataDefs::ODTrackAvgTime()) )
     nu->setDoubleData( "tammy_muon_odTrackAvgTime",particle->getDoubleData(ParticleExtraDataDefs::ODTrackAvgTime()) );

   if( particle->hasDoubleData(ParticleExtraDataDefs::ODClustersAvgTime()) )
     nu->setDoubleData( "tammy_muotammy_n_odClustersAvgTime",particle->getDoubleData(ParticleExtraDataDefs::ODClustersAvgTime()) );
   return;
}

//-----------------------------------------------------------
// calculate the muon momentum using proton kinematics
//-----------------------------------------------------------
void CCDeltaPlusAna::fillMuonMomentumFromProton( Minerva::NeutrinoInt* nu, Minerva::TrackVect& tracks, 
                                                   double tammy_proton_theta, double tammy_proton_p ) const
{
   //! check if the minerva track projects inside of minos
   Minerva::State  stateAtMinos;
   Gaudi::XYZPoint extrapToMinos;
   Gaudi::XYZPoint extrapToMinerva;

   double minos_front_face = m_mmt->getZShift();
   m_propagateToMinos->propagate( tracks[tracks.size()-1], minos_front_face, stateAtMinos );

   double minos_u, minos_v, minos_z;
   m_mmt->minervaToMinos( stateAtMinos.x(), stateAtMinos.y(), stateAtMinos.z(), minos_u, minos_v, minos_z);

   bool hit_partial = m_mmt->insideMinosPartialPlane(minos_u/1000.0,minos_v/1000.0,0.1);

   std::vector<double> tammy_minos_uv;
   tammy_minos_uv.push_back( minos_u );
   tammy_minos_uv.push_back( minos_v );

   //! get calculate muon energy
   double tammy_muon_theta = tracks[0]->theta();
   double tammy_muon_p     = tammy_proton_p * sin(tammy_proton_theta) / sin(tammy_muon_theta);

   //! get muon track absorbers
   std::vector<Minerva::AbsorberIDPtr> id_absorbers;
   for(unsigned int trk = 0; trk < tracks.size(); trk++) {

     bool projectToDetectorEnd = false;
     if(trk == tracks.size()-1) projectToDetectorEnd = true;

     std::vector<Minerva::AbsorberIDPtr> tmp_absorber;
     if( m_absorberStacker->AbsorberStackerGenerator(tracks[trk],tmp_absorber,projectToDetectorEnd).isSuccess() ) {
       for(std::vector<Minerva::AbsorberIDPtr>::iterator itAb = tmp_absorber.begin(); itAb != tmp_absorber.end(); ++itAb) {
         id_absorbers.push_back( *itAb );
       }
     }
   }
    
   if( id_absorbers.empty() ) return;

   //! energy conserved
   bool energyOk = true;
   if( std::isnan(tammy_muon_p) ) energyOk = false;

   //! calculate the energy loss per absorber
   double ptot = tammy_muon_p;
   double etot = sqrt( pow(ptot,2) + pow(MinervaUnits::M_muon,2) );
   double ekin = etot - MinervaUnits::M_muon; 

   if( energyOk ) { 
     std::vector<Minerva::AbsorberIDPtr>::iterator it = id_absorbers.begin();
     for( ; it != id_absorbers.end(); it++) {
       double par[4] = { etot, -1, -1, -1 }; //! energy, sigma, dE, pathLength
       double path   = (*it)->getPathLength();

       ILVolume::Intersection intersection = (*it)->getIntersection();
       if( !m_energyLoss->EnergyLossCalculator(intersection,path,MinervaUnits::M_muon,par).isSuccess() ) continue;
       if( par[2] < 0 ) { energyOk = false; break; }

       ekin = ekin - par[2];
       etot = ekin + MinervaUnits::M_muon;
     }
   }

   double p = sqrt( pow(etot,2) - pow(MinervaUnits::M_muon,2) );

   nu->setIntData("tammy_pOK",(int)energyOk);
   nu->setIntData("tammy_inside_minos_partial_plane",(int)hit_partial);

   nu->setDoubleData("tammy_calc_muon_p",tammy_muon_p);
   nu->setDoubleData("tammy_exit_muon_p",p);

   nu->setContainerDoubleData("minos_uv",tammy_minos_uv);
   return;
}

//--------------------------------------------------------------------
//! search for michels at vertices in the event
//--------------------------------------------------------------------
bool CCDeltaPlusAna::tagMichelElectrons( Minerva::PhysicsEvent* event ) const
{
   verbose() << "Enter CCDeltaPlusAna::tagMichelElectrons" << endmsg;
   bool foundMichel = false;
  
   //-- get the tagged muon prong
   Minerva::ProngVect primaryProngs = event->primaryProngs();
   SmartRef<Minerva::Prong> muonProng; 

   for(ProngVect::iterator itProngs = primaryProngs.begin(); itProngs != primaryProngs.end(); ++itProngs) {
     if( (*itProngs)->filtertaglist()->filterTagExists( "PrimaryMuon" ) ) {
       muonProng = *itProngs;
       break;
     }
   }

   std::vector< int > tammy_has_michel_category;
   std::vector< int > tammy_has_michel_vertex_type; // 0: int vtx, 1: proton's sec vtx, 2: proton's trk endpt
   std::vector< int > tammy_has_michel_in_vertex_point; //vertex number
   std::vector< double > tammy_has_michel_distance;
   std::vector< double > tammy_has_michel_energy;
   std::vector< int >    tammy_has_michel_ndigits;
   std::vector< double > tammy_has_michel_time_diff;

   //-- look for michel electrons in vertices
   Minerva::Prong tammy_vtx_michelProng;
   const SmartRefVector<Minerva::Vertex> vertices = event->select<Minerva::Vertex>();

   verbose() << "   Looking for michels in "  << vertices.size()-1 << " vertex points" << endmsg;
   for(SmartRefVector<Minerva::Vertex>::const_iterator itVtx = vertices.begin(); itVtx != vertices.end(); ++itVtx) {
    
     //-- no need to look at the endpoint of the muon track
     if( (*itVtx)->getIncomingTrack() == muonProng->minervaTracks()[0] ) continue;

     bool michel = true;
     Gaudi::XYZPoint pos = (*itVtx)->position();
     int x = pos.x(); int y = pos.y(); int z = pos.z();
     double t = (*itVtx)->getTime();
     // Follow Dun's procedure to loop over vertices with different times
     while( michel )
     {
       michel = m_michelTool->findMichel( x, y, z, t, tammy_vtx_michelProng );
       if (!michel) break;

       int tammy_vtx_type = -1;

       if( (*itVtx)->type() == Minerva::Vertex::StartPoint )
         tammy_vtx_type = 0;
       else if ( (*itVtx)->type() == Minerva::Vertex::Kinked )
         tammy_vtx_type = 1;
       else if ( (*itVtx)->type() == Minerva::Vertex::StopPoint )
         tammy_vtx_type = 2;
       else
         debug()<<"Could not determine vertex point!"<<endmsg;

       tammy_has_michel_vertex_type.push_back( tammy_vtx_type);
       tammy_has_michel_category.push_back( tammy_vtx_michelProng.getIntData("category") );
       tammy_has_michel_in_vertex_point.push_back( itVtx - vertices.begin() );
       tammy_has_michel_energy.push_back( tammy_vtx_michelProng.getDoubleData("energy") );
       tammy_has_michel_distance.push_back( tammy_vtx_michelProng.getDoubleData("distance") );
       tammy_has_michel_ndigits.push_back( tammy_vtx_michelProng.getIntData("ndigits") );
       tammy_has_michel_time_diff.push_back( tammy_vtx_michelProng.getDoubleData("time_diff") );

       verbose() << "    found a michel " << endmsg;
       t += tammy_vtx_michelProng.getDoubleData("meantime_diff") + 100;
     }

   }

   //-- store values in ntuple
   event->setContainerIntData("tammy_has_michel_category",tammy_has_michel_category);
   event->setContainerIntData("has_michel_vertex_type",tammy_has_michel_vertex_type);
   event->setContainerIntData("has_michel_in_vertex_point",tammy_has_michel_in_vertex_point);
   event->setContainerDoubleData("has_michel_distance",tammy_has_michel_distance);
   event->setContainerDoubleData("has_michel_energy",tammy_has_michel_energy);
   event->setContainerIntData("has_michel_ndigits",tammy_has_michel_ndigits);
   event->setContainerDoubleData("has_michel_time_diff",tammy_has_michel_time_diff);
   
   if (tammy_has_michel_category.size()>0) foundMichel = true;

   verbose() << "Exit CCDeltaPlusAna::tagMichelElectrons" << endmsg;
   return foundMichel;
}

