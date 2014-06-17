/*
    See CCProtonPi0.h header for Class Information
*/
#include <cassert>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <numeric>

// Local Libraries
#include "CCProtonPi0.h"
#include "AngleScan.h"
#include "ClusterVectorInfo.h"
#include "DigitVectorTruthInfo.h"
#include "OneParLineFit.h"
#include "TwoParLineFit.h"

#include "CCProtonPi0/IHoughBlob.h"
#include "CCProtonPi0/IHoughTool.h"

// ROOT Libraries
#include "TString.h"
#include "TDatabasePDG.h"
#include "TVector3.h"
#include "TMath.h"
#include "TRandom3.h"

// Minerva Analysis Framework
#include "RecInterfaces/IFiducialPointTool.h"
#include "RecInterfaces/IRecoObjectTimeTool.h"

#include "EnergyRecTools/IEnergyCorrectionTool.h"
#include "EnergyRecTools/IExtraEnergyTool.h"

#include "ProngMaker/IMichelTool.h"
#include "ProngMaker/IProngClassificationTool.h"
#include "ProngMaker/IODProngClassificationTool.h"

#include "ParticleMaker/IParticleMakerTool.h"
#include "ParticleMaker/IParticleTool.h"

#include "AnaUtils/IProtonUtils.h"
#include "AnaUtils/ICCPionIncUtils.h"
#include "AnaUtils/IMCTrackTool.h"
#include "AnaUtils/MCTrack.h"

#include "MinervaUtils/IMinervaMathTool.h"
#include "MinervaUtils/IHitTaggerTool.h"
#include "MinervaUtils/IMinervaObjectAssociator.h"

#include "MinervaDet/IGeomUtilSvc.h"
#include "MinervaDet/DeDetector.h"
#include "ODDet/DeOuterDetector.h"

#include "BlobFormation/IIDAnchoredBlobCreator.h"
#include "BlobFormation/IBlobCreatorUtils.h"

#include "CalTools/IGetCalAttenuation.h"

#include "GeoUtils/IMinervaCoordSysTool.h"

#include "BadChannels/IGetDeadTime.h"

#include "GiGaCnv/IGiGaGeomCnvSvc.h"
#include "G4Navigator.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Material.hh"

//==============================================================================
// Global variables
//==============================================================================
/*
   Global Variables used instead of class data members so that they can be
   assigned in 'const' methods. In general, we can use 'mutable' members,
   but it does not seem to work in the framework 
*/
SmartRef<Minerva::Track>    m_MuonTrack;
SmartRef<Minerva::Prong>    m_MuonProng;
SmartRef<Minerva::Particle> m_MuonParticle;
SmartRef<Minerva::Vertex>   m_PrimaryVertex;
std::map<int, Minerva::TG4Trajectory*> fTrajectoryMap;

namespace
{
    class increasingDistanceFromVertex : public std::binary_function
    <
        SmartRef<Minerva::IDCluster>,
        SmartRef<Minerva::IDCluster>,
        bool
        > {
    public:
        explicit increasingDistanceFromVertex(Gaudi::XYZPoint pos) : fVertex(pos) {}
        
        bool operator() (const SmartRef<Minerva::IDCluster>& lhs,
                         const SmartRef<Minerva::IDCluster>& rhs) const {

            
            const double z0 = fVertex.Z();
            std::cout<<z0<<std::endl;

            return std::abs(lhs->z()-z0) < std::abs(rhs->z()-z0);
        }

    private:
        Gaudi::XYZPoint fVertex;
    };

}


DECLARE_TOOL_FACTORY( CCProtonPi0 );

using namespace Minerva;
    
//==============================================================================
// Standard constructor 
//==============================================================================
CCProtonPi0::CCProtonPi0(const std::string& type, const std::string& name, const IInterface* parent ) :
MinervaAnalysisTool( type, name, parent ) 
{
    info() <<"--------------------------------------------------------------------------"<<endmsg;
    info() <<"Enter CCProtonPi0::CCProtonPi0() -- Default Constructor" << endmsg;
    
    declareInterface<IInteractionHypothesis>(this);
    //! mandatory declaration of analysis signature: CCProtonPi0
    m_anaSignature = "CCProtonPi0";
    
    // Private Properties
    declareProperty("StoreAllEvents",      m_store_all_events =    true);
    declareProperty("DoPlausibilityCuts",  m_doPlausibilityCuts =  true);
    declareProperty("MakeShortTracks",     m_makeShortTracks =     true);
    
    declareProperty("BeamAngleBias",       m_beamAngleBias = 0.006*CLHEP::radian );
    
    declareProperty("MinMuonScore",        m_minMuonScore = 0.9 );
    declareProperty("MinProtonScore",      m_minProtonScore = 0.0 );
    
    declareProperty("MuonProngColor",              m_muonProngColor       = 0x228B22); //-- green
    declareProperty("ProtonProngColor",            m_protonProngColor     = 0x9932CC); //-- purple
    declareProperty("PrimaryVertexProngColor",     m_primaryVertexColor   = 0xFF0000); //-- red
    declareProperty("SecondaryVertexProngColor",   m_secondaryVertexColor = 0xFFA500); //-- orange
    declareProperty("TrackEndProngColor",          m_endPointVertexColor  = 0xFF1493); //-- pink
    declareProperty("UnattachedProngColor",        m_unattachedProngColor = 0x0000FF); //-- blue
    
    
    declareProperty("MichelTrkToolAlias",   m_michelTrkToolAlias    = "CCProtonPi0_MichelTrackTool");
    declareProperty("MichelVtxToolAlias",   m_michelVtxToolAlias    = "CCProtonPi0_MichelVertexTool");
    declareProperty("ProtonUtilsAlias",     m_protonUtilsAlias      = "CCProtonPi0_ProtonUtils");
    declareProperty("ParticleToolName",     m_particleToolName      = "dEdXTool" );
    declareProperty("ParticleToolAlias",    m_particleToolAlias     = "dEdXTool" );
    
    // Vertex blob
    declareProperty( "SkipLowEnergyClusterVtxEnergy", fSkipLowEnergyClusterVtxEnergy = true);
    declareProperty( "ThresholdVertexEnergy",         fThresholdVertexEnergy = 1.0*CLHEP::MeV);
    declareProperty( "UseSphereVertex", 		m_sphereVertex	= true );
    declareProperty( "SphereMaxSearchDistance", 	m_maxSearchD  	=  90.0 * CLHEP::millimeter ); 
    declareProperty( "SphereMaxStartingDistance", m_maxStartingD	=  90.0 * CLHEP::millimeter ); 
    declareProperty( "SphereMaxAllowedSearchGap", m_maxSearchGap	= 180.0 * CLHEP::millimeter ); 
    declareProperty( "UseFilamentVertex", 	m_filamentVertex = false );
    declareProperty( "FilamentMaxSearchDistance",   m_maxSearchDFila    = 500.0 * CLHEP::millimeter );
    declareProperty( "FilamentMaxStartingDistance", m_maxStartingDFila  = 68.0 * CLHEP::millimeter ); 
    declareProperty( "FilamentMaxAllowedSearchGap", m_maxSearchGapFila  = 51.0 * CLHEP::millimeter ); 
    declareProperty( "FilterClusterTypes",   m_filterClusterTypes    = true );
    
    // Ana tools
    declareProperty( "qOverpChargeCut",       m_qOverpChargeCut       = 0   );
    declareProperty( "HoughEnergyLimit",      m_energyHoughlimit      = 900 * CLHEP::MeV ); 
    declareProperty( "RejectedClustersTime",  m_rejectedClustersTime  = 25 * CLHEP::ns );
    
    
    /* 
        Number from Cesar Sotelo's studies.
        For forward-going track, the cylinder starts 20cm DOWNSTREAM of the first track node
    */
    declareProperty( "ExtraEnergyCylinderUpstreamOffset",  m_extraEnergyCylinderUpstreamLength = -200.0*CLHEP::mm);
    declareProperty( "ExtraEnergyCylinderDownstreamOffset",m_extraEnergyCylinderDownstreamLength = 50.0*CLHEP::mm);
    declareProperty( "ExtraEnergyCylinderRadius",          m_extraEnergyCylinderRadius = 50*CLHEP::mm);
    declareProperty( "ExtraEnergyLowerTimeWindow",         m_extraEnergyLowerTimeWindow = 25.0*CLHEP::ns);
    declareProperty( "ExtraEnergyUpperTimeWindow",         m_extraEnergyUpperTimeWindow = 25.0*CLHEP::ns);
    declareProperty( "ExtraEnergyPhotoElectronCut",        m_extraEnergyPECut = 0.0);
    
    declareProperty( "new_impl", new_impl_ = true);
    declareProperty( "try_to_recover", try_to_recover_ = false);
    declareProperty( "attenuation_correction", attenuation_correction_ = true);
    
    declareProperty( "UVMatchTolerance", fUVMatchTolerance = 10.0*CLHEP::mm);
    declareProperty( "UVMatchMoreTolerance", fUVMatchMoreTolerance = 100.0*CLHEP::mm);
    declareProperty( "AllowUVMatchWithMoreTolerance", fAllowUVMatchWithMoreTolerance = true);
    

    // Protected properties from IInteractionHypothesis.
    m_hypMeths.push_back( m_anaSignature );
    declareProperty("HypothesisMethods", m_hypMeths);
    info() << " CCProtonPi0 Hypothesis added " << endmsg;
    
    
    info() <<"Exit CCProtonPi0::CCProtonPi0() -- Default Constructor" << endmsg;
    info() <<"--------------------------------------------------------------------------"<<endmsg;
}
    
//==============================================================================
// Initialize
//==============================================================================
StatusCode CCProtonPi0::initialize()
{
    info() <<"--------------------------------------------------------------------------"<<endmsg;
    info() <<"Enter CCProtonPi0::initialize()" << endmsg;
    
    //! Initialize the base class.  This will fail if you did not define m_anaSignature.
    StatusCode sc = this->MinervaAnalysisTool::initialize();
    if( sc.isFailure() ) { 
        return Error( "Failed to initialize!", sc ); 
    }
    info()<<"   initialized MinervaAnalysisTool"<<endmsg;
    
    //Seed the RNG
    m_randomGen = new TRandom3( m_randomSeed );   
    if ( m_randomGen == NULL ) return StatusCode::FAILURE;
    
    
    m_detectableGammaE      = 0;        // MeV
    m_detectablePi0KE       = 0;        // MeV
    m_detectableProtonKE    = 120;      // MeV
    
    // Fiducial Volume
    m_fidHexApothem  = 850.0*CLHEP::mm;
    m_fidUpStreamZ   = 5990.0*CLHEP::mm;   // ~middle of module 27, plane 1
    m_fidDownStreamZ = 8340.0*CLHEP::mm;   // ~middle of module 79, plane 1
    
    // Reconstructable Volume
    m_recoHexApothem  = 1000.0*CLHEP::mm; 
    m_recoUpStreamZ   = 5750.0*CLHEP::mm;
    m_recoDownStreamZ = 8700.0*CLHEP::mm;
    
    
    //! Initializing Analysis Tools
    try {
        m_InnerDetector = getDet<Minerva::DeDetector>("/dd/Structure/Minerva/Detector/InnerDetector");
    } catch(GaudiException& e) {
        error() << "Could not obtain detector: Inner Detector" << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{ m_OuterDetector = getDet<Minerva::DeOuterDetector>("/dd/Structure/Minerva/Detector/OuterDetector"); }
    catch( GaudiException& e ) {
        error() << "Could not retrieve the Outer Detector!" << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{
        m_coordSysTool = tool<IMinervaCoordSysTool>("MinervaCoordSysTool");
    }catch(GaudiException& e){
        error() << "Could not obtain tool: MinervaCoordSysTool" << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{ 
        m_protonUtils = tool<IProtonUtils>("ProtonUtils", m_protonUtilsAlias); 
    }catch( GaudiException& e ){
        error() << "Could not obtain tool: ProtonUtils" << m_protonUtilsAlias<< endmsg;
        return StatusCode::FAILURE;
    }
    
    try{ 
        m_prongIntersection = tool<IProngClassificationTool>("ProngIntersectionTool"); 
    }catch( GaudiException& e ){
        error() << "Could not obtain ProngClassificationTool: " << m_prongIntersection << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{ 
        m_particleMaker = tool<IParticleMakerTool>("ParticleMakerTool"); 
    }catch( GaudiException& e ){
        error() << "Could not obtain ParticleMakerTool: " << m_particleMaker << endmsg;
        return StatusCode::FAILURE;
    }
      
    try {
        m_particleTool = tool<IParticleTool>(m_particleToolName, m_particleToolAlias);
    } catch(GaudiException& e) {
        error() << "Could not obtain particleTool" <<endmsg;
        return StatusCode::FAILURE;
    }
    
    try{ 
        m_odMatchTool = tool<IODProngClassificationTool>("ODTrackMatchTool"); 
    }catch( GaudiException& e ){
        error() << "Could not obtain ProngClassificationTool: " << m_odMatchTool << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{ 
        m_energyCorrectionTool = tool<IEnergyCorrectionTool>("EnergyCorrectionTool"); 
    }catch( GaudiException& e ){
        error() << "Could not obtain tool: EnergyCorrectionTool" << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{ 
        m_hitTagger = tool<IHitTaggerTool>( "HitTaggerTool" ); 
    }
    catch( GaudiException& e ) {
        error() << "Could not obtain tool: HitTaggerTool" << endmsg;
        return StatusCode::FAILURE;
    }
   
    try {
        m_michelTrkTool = tool<IMichelTool>("MichelTool",m_michelTrkToolAlias);
    } catch(GaudiException& e){
        error()<<"Could not obtain tool: MichelTool" << endmsg;
        return StatusCode::FAILURE;
    } 
    
    try {
        m_michelVtxTool = tool<IMichelTool>("MichelTool",m_michelVtxToolAlias);
    } catch(GaudiException& e){
        error()<<"Could not obtain tool: MichelTool" << endmsg;
        return StatusCode::FAILURE;
    }
     
    try {
        m_ccPionIncUtils = tool<ICCPionIncUtils>("CCPionIncUtils");
    } catch( GaudiException& e ) {
        error() << "Could not obtain tool: CCPionIncUtils" << endmsg;
        return StatusCode::FAILURE;
    }
    
    try {
        m_objectAssociator = tool<IMinervaObjectAssociator>("MinervaObjectAssociator");
    } catch( GaudiException& e ) {
        error() << "Could not obtain tool: MinervaObjectAssociator" << endmsg;
        return StatusCode::FAILURE;
    }
    
    try {
        m_caloUtils = tool<ICalorimetryUtils>("CalorimetryUtils");  
    } catch( GaudiException& e ) {    
        error() << "Could not obtain tool: CalorimetryUtils!" << endmsg;    
        return StatusCode::FAILURE;  
    }
      
    try {
        m_stopPointBlobTool = tool<IIDAnchoredBlobCreator>("VertexBlobCreator", "CCDeltaPlusVtxBlobCreator");
    } catch( GaudiException& e ){
        error() <<"Could not obtain VertexBlobCreator"<<endmsg;
        return StatusCode::FAILURE;
    }
    
    try {
        m_extraEnergyTool = tool<IExtraEnergyTool>("ExtraEnergyTool");
    } catch( GaudiException& e ) {
        error() << "Could not obtain tool: ExtraEnergyTool" << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{
        m_getDeadTimeTool = tool<IGetDeadTime>( "GetDeadTime" );
    } catch(GaudiException& e){
        error()<<"Could not obtain tool: GetDeadTime" << endmsg;
        return StatusCode::FAILURE;
    }
    
    try {
        m_recoTimeTool = tool<IRecoObjectTimeTool>( "RecoObjectTimeTool" );
    }catch (GaudiException& e) {
        error() << "Could not obtain tool: RecoObjectTimeTool" << endmsg;
        return StatusCode::FAILURE;
    } 
    
    try {
        m_MCTrackTool = tool<IMCTrackTool>("MCTrackTool");
    } catch (GaudiException& e) {
        error() << "Could not obtain tool: MCTrackTool" << endmsg;
        return StatusCode::FAILURE;
    }
    
    try {
        m_mathTool = tool<IMinervaMathTool>("MinervaMathTool");
    }
    catch( GaudiException& e ) {
        error() << "Could not obtain tool: MinervaMathTool!" << endmsg;
        return StatusCode::FAILURE;
    }
    
    //! Blob Tools
    try{
        m_blobUtils = tool<IBlobCreatorUtils>("BlobCreatorUtils");
    } 
    catch( GaudiException& e ){
        error() << "Could not obtain tool: BlobCreatorUtils! " << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{
        m_idHoughBlob = tool<IHoughBlob>("HTBlob");
    }
    catch (GaudiException& e){
        error() << " Could not obtain tool: HTBlob " << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{
        m_idHoughTool = tool<IHoughTool>("HTtool");
    }
    catch (GaudiException& e){
        error() << " Could not obtain tool: HTtool " << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{
        m_idConeScanBlob = tool<IIDAnchoredBlobCreator>("ConeScanIDBlobCreator");
    } 
    catch (GaudiException& e){
        error() << " Could not obtain tool: ConeScanIDBlobCreator " << endmsg;
        return StatusCode::FAILURE;
    }
    
    try {
        m_AttenuationCorrectionTool = tool<IGetCalAttenuation>("GetCalAttenuation", "AttenuationCorrectionTool");
    } catch( GaudiException& e ){
        error() << "Could not obtain tool: CC1pi0AttenuationCorrection tool" << endmsg; 
        return StatusCode::FAILURE;
    }
  
    
    debug() <<" Obtained all tools!"<<endmsg;
    
    //! Services
    try{
        m_gigaCnvSvc = svc<IGiGaGeomCnvSvc>("GiGaGeo", true);
    } catch( GaudiException& e ){
        error() <<"Could not obtain GiGaGeo"<<endmsg;
        return StatusCode::FAILURE;
    }
    
    service("GeomUtilSvc", m_GeomUtilSvc, true);
    m_idDet = m_GeomUtilSvc->getIDDet();
    m_odDet = m_GeomUtilSvc->getODDet();
    
    fTrajectoryMap.clear();
    
    
    //! declare common branches
    declareCommonPhysicsAnaBranches();
    declareMinosMuonBranches();
    declareGenieWeightBranches();
    
    //--------------------------------------------------------------------------
    //! Select the branches you want in your AnaTuple
    //--------------------------------------------------------------------------
    
    //! truth  
    declareBoolTruthBranch( "reco_hasGoodObjects");
    declareBoolTruthBranch( "reco_isGoodVertex");  
    declareBoolTruthBranch( "reco_isWellFitVertex");
    declareBoolTruthBranch( "reco_isFidVol");
    declareBoolTruthBranch( "reco_isFidVol_smeared");
    declareBoolTruthBranch( "reco_isMinosMatch");
    declareBoolTruthBranch( "reco_isBrokenTrack");
    declareIntTruthBranch(  "reco_muonCharge", 0);
    
    declareBoolTruthBranch( "isSignal");
    declareBoolTruthBranch( "isFidVol");
    declareBoolTruthBranch( "isPlausible");
    
    declareIntTruthBranch( "vertex_module", 500);
    declareIntTruthBranch( "vertex_plane", 0);
    declareIntTruthBranch( "target_material", -1);
    
    declareIntTruthBranch("N_proton", -1 );
    declareIntTruthBranch("N_neutron", -1 );
    declareIntTruthBranch("N_muminus", -1 );
    declareIntTruthBranch("N_muplus", -1 );
    declareIntTruthBranch("N_pi0", -1 );
    declareIntTruthBranch("N_piplus", -1 );
    declareIntTruthBranch("N_piminus", -1 );
    declareIntTruthBranch("N_deltaplus", -1 );
    declareIntTruthBranch("N_gamma", -1 );
    declareIntTruthBranch("N_other", -1 );

    declareDoubleTruthBranch("muon_px", 0 );
    declareDoubleTruthBranch("muon_py", 0 );
    declareDoubleTruthBranch("muon_pz", 0 );
    declareDoubleTruthBranch("muon_E", -9.0 );
    declareDoubleTruthBranch("muon_theta_wrtbeam", -9.0 );
    declareIntTruthBranch("muon_charge", 0 );
    
    declareContainerDoubleTruthBranch("proton_px", 20, 0 );
    declareContainerDoubleTruthBranch("proton_py", 20, 0 );
    declareContainerDoubleTruthBranch("proton_pz", 20, 0 );
    declareContainerDoubleTruthBranch("proton_E",  20, -9.0 );
    declareContainerDoubleTruthBranch("proton_vtx_x", 20, 0 );
    declareContainerDoubleTruthBranch("proton_vtx_y", 20, 0 );
    declareContainerDoubleTruthBranch("proton_vtx_z", 20, 0 );
    declareContainerDoubleTruthBranch("proton_theta_wrtbeam",  20, -9.0 );
    declareContainerIntTruthBranch("proton_trackID", 20, -1 );    
    declareContainerIntTruthBranch("proton_parentID", 20, -1 );
    
    declareContainerDoubleTruthBranch("pi0_px", 20, 0 );
    declareContainerDoubleTruthBranch("pi0_py", 20, 0 );
    declareContainerDoubleTruthBranch("pi0_pz", 20, 0 );
    declareContainerDoubleTruthBranch("pi0_E",  20, -9.0 );
    declareContainerDoubleTruthBranch("pi0_vtx_x", 20, 0 );
    declareContainerDoubleTruthBranch("pi0_vtx_y", 20, 0 );
    declareContainerDoubleTruthBranch("pi0_vtx_z", 20, 0 );    
    declareContainerDoubleTruthBranch("pi0_theta_wrtbeam",  20, -9.0 );
    declareContainerIntTruthBranch("pi0_trackID", 20, -1 );
    declareContainerIntTruthBranch("pi0_parentID", 20, -1 );
    
    declareContainerDoubleTruthBranch("gamma_px", 20, 0 );
    declareContainerDoubleTruthBranch("gamma_py", 20, 0 );
    declareContainerDoubleTruthBranch("gamma_pz", 20, 0 );
    declareContainerDoubleTruthBranch("gamma_E",  20, -9.0 );
    declareContainerDoubleTruthBranch("gamma_vtx_x", 20, 0 );
    declareContainerDoubleTruthBranch("gamma_vtx_y", 20, 0 );
    declareContainerDoubleTruthBranch("gamma_vtx_z", 20, 0 );    
    declareContainerDoubleTruthBranch("gamma_theta_wrtbeam",  20, -9.0 );
    declareContainerIntTruthBranch("gamma_trackID", 20, -1 );
    declareContainerIntTruthBranch("gamma_parentID", 20, -1 );
    
    //! Truth Match for Prongs
    declareIntBranch(m_hypMeths,    "isMuonInsideOD",        -1);
    declareIntBranch(m_hypMeths,    "ntrajMuonProng",        -1);
    declareIntBranch(m_hypMeths,    "trajMuonProngPrimary",  -1);
    declareIntBranch(m_hypMeths,    "trajMuonProngPDG",      -1);
    declareDoubleBranch(m_hypMeths, "trajMuonProngPx", -1);
    declareDoubleBranch(m_hypMeths, "trajMuonProngPy", -1);
    declareDoubleBranch(m_hypMeths, "trajMuonProngPz", -1);
    declareDoubleBranch(m_hypMeths, "trajMuonProngEnergy", -1);
    declareDoubleBranch(m_hypMeths, "trajMuonProngMomentum", -1);
    declareDoubleBranch(m_hypMeths, "trajMuonProngPSelf", -1);
    declareDoubleBranch(m_hypMeths, "trajMuonTheta",         -1);
    declareDoubleBranch(m_hypMeths, "trajMuonPhi",           -1);
    declareDoubleBranch(m_hypMeths, "endMuonTrajMomentum",   -9999);
    declareDoubleBranch(m_hypMeths, "endMuonTrajXPosition",  -9999);
    declareDoubleBranch(m_hypMeths, "endMuonTrajYPosition",  -9999);
    declareDoubleBranch(m_hypMeths, "endMuonTrajZPosition",  -9999);
    
    declareContainerIntBranch(m_hypMeths,    "isProtonInsideOD",        10, -1);
    declareContainerIntBranch(m_hypMeths,    "ntrajProtonProng",         10, -1);
    declareContainerIntBranch(m_hypMeths,    "trajProtonProngPrimary",  10, -1);
    declareContainerIntBranch(m_hypMeths,    "trajProtonProngPDG",      10, -1);
    declareContainerDoubleBranch(m_hypMeths, "trajProtonProngPx", 10, -1);
    declareContainerDoubleBranch(m_hypMeths, "trajProtonProngPy", 10, -1);
    declareContainerDoubleBranch(m_hypMeths, "trajProtonProngPz", 10, -1);
    declareContainerDoubleBranch(m_hypMeths, "trajProtonProngEnergy", 10, -1);
    declareContainerDoubleBranch(m_hypMeths, "trajProtonProngMomentum", 10, -1);
    declareContainerDoubleBranch(m_hypMeths, "trajProtonProngPSelf", 10, -1);
    declareContainerDoubleBranch(m_hypMeths, "trajProtonTheta",         10, -1);
    declareContainerDoubleBranch(m_hypMeths, "trajProtonPhi",           10, -1);
    declareContainerDoubleBranch(m_hypMeths, "endProtonTrajMomentum",   10, -9999);
    declareContainerDoubleBranch(m_hypMeths, "endProtonTrajXPosition",  10, -9999);
    declareContainerDoubleBranch(m_hypMeths, "endProtonTrajYPosition",  10, -9999);
    declareContainerDoubleBranch(m_hypMeths, "endProtonTrajZPosition",  10, -9999);
    
    
    //! Event - Cut Results
    declareIntEventBranch( "Cut_Vertex_None", -1 );
    declareIntEventBranch( "Cut_Vertex_Null", -1 );
    declareIntEventBranch( "Cut_Vertex_Not_Reconstructable", -1 );
    declareIntEventBranch( "Cut_Vertex_Not_Fiducial", -1 );
    declareIntEventBranch( "Cut_Muon_None",-1);
    declareIntEventBranch( "Cut_Muon_Score_Low",-1);
    declareIntEventBranch( "Cut_Muon_Charge",-1);
    declareIntEventBranch( "Cut_Vertex_Michel_Exist", -1 );
    declareIntEventBranch( "Cut_EndPoint_Michel_Exist", -1 );
    declareIntEventBranch( "Cut_secEndPoint_Michel_Exist", -1 );
    declareIntEventBranch( "Cut_Proton_None", -1 );
    declareIntEventBranch( "Cut_PreFilter_Pi0", -1 );
    
    //! Event - General reco
    declareIntEventBranch( "n_startpoint_vertices", -1 );
    declareIntEventBranch( "n_long_tracks", -1);
    declareIntEventBranch( "n_short_tracks", -1);
    declareIntEventBranch( "n_anchored_long_trk_prongs", -1 );
    declareIntEventBranch( "n_anchored_short_trk_prongs", -1 );
    declareIntEventBranch( "n_vtx_prongs", -1 );
    declareIntEventBranch( "n_iso_trk_prongs", -1);
    declareIntEventBranch( "n_iso_blob_prongs", -1);  
    declareIntEventBranch( "n_dsp_blob_prongs", -1);
    declareIntEventBranch( "n_us_muon_clusters", 0);
    declareDoubleEventBranch( "time", -1.0 );
  
    //! Event - Michel 
    declareIntEventBranch( "n_vtx_michel_views", 0 );
    declareDoubleEventBranch( "vtx_michel_distance", -1.0);
    
    //! Event - Energy
    declareDoubleEventBranch( "muonVisibleE", -1.0 );
    declareDoubleEventBranch( "hadronVisibleE", -1.0 );
    declareDoubleEventBranch( "totalVisibleE", -1.0 );
    declareDoubleEventBranch( "totalIDVisibleE", -1.0 );
    declareDoubleEventBranch( "totalODVisibleE", -1.0 );
    declareDoubleEventBranch( "vtxBlobExtraE",  -1.0 );
    declareDoubleEventBranch( "unattachedExtraE", -1.0 );
    declareDoubleEventBranch( "dispersedExtraE", -1.0 );
     
    declareBoolEventBranch( "isMinosMatchTrack");
    declareBoolEventBranch( "isMinosMatchStub");
    declareBoolEventBranch( "well_fit_vertex");
    declareDoubleEventBranch( "well_fit_vertex_angle", -99.9 );
    declareBoolEventBranch( "isBrokenTrack");
    
    //! Event - Pi0 
    // May convert to NuInt type variable
    // Additional Branches CCProtonPi0 Styles
    declareDoubleEventBranch("pi0_px",-1.0);
    declareDoubleEventBranch("pi0_py",-1.0);
    declareDoubleEventBranch("pi0_pz",-1.0);
    declareDoubleEventBranch("pi0_E",-1.0);
    
    declareDoubleEventBranch("gamma1_px",-1.0);
    declareDoubleEventBranch("gamma1_py",-1.0);
    declareDoubleEventBranch("gamma1_pz",-1.0);
    declareDoubleEventBranch("gamma1_E",-1.0);
    
    declareDoubleEventBranch("gamma2_px",-1.0);
    declareDoubleEventBranch("gamma2_py",-1.0);
    declareDoubleEventBranch("gamma2_pz",-1.0);
    declareDoubleEventBranch("gamma2_E",-1.0);

    // CCPi0 Default Branches
    declareDoubleEventBranch("pienergy",1.e6);
    declareDoubleEventBranch("pitheta", 1.e6);
    declareDoubleEventBranch("piphi",   1.e6);
    declareDoubleEventBranch("pithetax", 1.e6);
    declareDoubleEventBranch("pithetay",   1.e6);
    declareDoubleEventBranch("oangle",  1.e6);
    
    declareDoubleEventBranch("mgg", 1.e6);
    
    declareDoubleEventBranch("pithetab", 1.e6);
    declareDoubleEventBranch("piphib",   1.e6);
    declareDoubleEventBranch("pithetaxb", 1.e6);
    declareDoubleEventBranch("pithetayb",   1.e6);
    
    declareDoubleEventBranch("g1e",     1.e6);
    declareDoubleEventBranch("g1theta", 1.e6);
    declareDoubleEventBranch("g1phi",   1.e6);
    declareDoubleEventBranch("g2e",     1.e6);
    declareDoubleEventBranch("g2theta", 1.e6);
    declareDoubleEventBranch("g2phi",   1.e6);
    
    declareContainerDoubleEventBranch("pimom",4,-1);
    declareContainerDoubleEventBranch("g1mom",3,-1);
    declareContainerDoubleEventBranch("g2mom",3,-1);
   
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
    
    //! Event - PreFilterPi0()
    declareDoubleEventBranch("evis_nearvtx_total", -1.0);
    declareDoubleEventBranch("evis_nearvtx_x", -1.0);
    declareDoubleEventBranch("evis_nearvtx_u", -1.0);
    declareDoubleEventBranch("evis_nearvtx_v", -1.0);
    
    declareDoubleEventBranch("evis_total_x", -1.0);
    declareDoubleEventBranch("evis_ntgt_x",  -1.0);
    declareDoubleEventBranch("evis_trkr_x",  -1.0);
    declareDoubleEventBranch("evis_ecal_x",  -1.0);
    declareDoubleEventBranch("evis_hcal_x",  -1.0);
    
    declareDoubleEventBranch("evis_total_u", -1.0);
    declareDoubleEventBranch("evis_ntgt_u",  -1.0);
    declareDoubleEventBranch("evis_trkr_u",  -1.0);
    declareDoubleEventBranch("evis_ecal_u",  -1.0);
    declareDoubleEventBranch("evis_hcal_u",  -1.0);
    
    declareDoubleEventBranch("evis_total_v", -1.0);
    declareDoubleEventBranch("evis_ntgt_v",  -1.0);
    declareDoubleEventBranch("evis_trkr_v",  -1.0);
    declareDoubleEventBranch("evis_ecal_v",  -1.0);
    declareDoubleEventBranch("evis_hcal_v",  -1.0);
    
    declareDoubleEventBranch("evis_total", -1.0);
    declareDoubleEventBranch("evis_ntgt",  -1.0);
    declareDoubleEventBranch("evis_trkr",  -1.0);
    declareDoubleEventBranch("evis_ecal",  -1.0);
    declareDoubleEventBranch("evis_hcal",  -1.0);
    declareDoubleEventBranch("evis_other", -1.0);
    
    //! Event - VtxBlob()
    declareContainerDoubleEventBranch("Vertex_energy_radii");
    declareDoubleEventBranch( "Vertex_blob_energy", -9999 );
    declareDoubleEventBranch( "Filament_Vertex_energy", -9999 );
    declareDoubleEventBranch( "Sphere_Vertex_energy", -9999 );
    
    //! Event - ConeBlobs()
    declareDoubleEventBranch( "RE_energy_Tracker", -9999 );
    declareDoubleEventBranch( "RE_energy_ECAL", -9999 );
    declareDoubleEventBranch( "RE_energy_HCAL", -9999 );
    
    declareIntEventBranch("nblob_anglescan",-1);
    declareIntEventBranch("nblob_hough", -1);
    declareIntEventBranch("anglescan_ncandx", -1);
    declareIntEventBranch("anglescan_ncand", -1);
    
    declareContainerDoubleEventBranch("mgg_vector");
    declareContainerDoubleEventBranch("good_mgg_vector");
    
    //! Event - Gamma dEdX Information
    declareIntEventBranch("g1dedx_nplane", -1);
    declareIntEventBranch("g1dedx_doublet", -1);
    declareIntEventBranch("g1dedx_empty_plane", -1);
    declareDoubleEventBranch("g1dedx_total", -1.0);
    declareDoubleEventBranch("g1dedx", -1.0);
    declareDoubleEventBranch("g1dedx_total1", -1.0);
    declareDoubleEventBranch("g1dedx1", -1.0);
    declareContainerIntEventBranch("g1dedx_cluster_occupancy");
    declareContainerDoubleEventBranch("g1dedx_cluster_energy");
    declareContainerDoubleEventBranch("g1dedx_rev_cluster_energy");
    
    declareIntEventBranch("g2dedx_nplane", -1);
    declareIntEventBranch("g2dedx_doublet", -1);
    declareIntEventBranch("g2dedx_empty_plane", -1);
    declareDoubleEventBranch("g2dedx_total", -1.0);
    declareDoubleEventBranch("g2dedx", -1.0);
    declareDoubleEventBranch("g2dedx_total1", -1.0);
    declareDoubleEventBranch("g2dedx1", -1.0);
    declareContainerIntEventBranch("g2dedx_cluster_occupancy");
    declareContainerDoubleEventBranch("g2dedx_cluster_energy");
    declareContainerDoubleEventBranch("g2dedx_rev_cluster_energy");
    
    //! Event -  Truth info of reconstructed EM showers
    // Will be modified in future
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
    
    //! Event - Mostly for debugging purposes will be trimmed in future
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
    
    declareContainerDoubleEventBranch("blob_cluster_energy1");
    declareContainerDoubleEventBranch("blob_cluster_energy2");
    
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
    
    declareDoubleEventBranch("g1recotrkrevis", -1.0);
    declareDoubleEventBranch("g1recoecalevis", -1.0);
    declareDoubleEventBranch("g1recohcalevis", -1.0);
    declareDoubleEventBranch("g1recoscalevis", -1.0);
    
    declareDoubleEventBranch("g2recotrkrevis", -1.0);
    declareDoubleEventBranch("g2recoecalevis", -1.0);
    declareDoubleEventBranch("g2recohcalevis", -1.0);
    declareDoubleEventBranch("g2recoscalevis", -1.0);
    
    //! Event - OD info
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
    
    //! Event - Vertex Information -- Filled in reconstructEvent()
    declareIntEventBranch("vtx_total_count",-1);
    declareIntEventBranch("vtx_secondary_count",-1);
    declareIntEventBranch("vtx_primary_index",-1);
    declareIntEventBranch("vtx_primary_multiplicity",-1);
    
    //! NeutrinoInt - Primary Vertex -- Filled in interpretEvent()
    declareIntBranch( m_hypMeths, "vtx_module", -99);
    declareIntBranch( m_hypMeths, "vtx_plane",-1);
    declareDoubleBranch( m_hypMeths, "vtx_x",0.0);
    declareDoubleBranch( m_hypMeths, "vtx_y",0.0);
    declareDoubleBranch( m_hypMeths, "vtx_z",0.0);
    
    //! NeutrinoInt - Muon -- Filled in setMuonParticleData()
    declareIntBranch( m_hypMeths, "muon_minervaTrack_types", -1);
    declareIntBranch( m_hypMeths, "muon_N_minosTracks", -1);
    declareIntBranch( m_hypMeths, "muon_minosTrackQuality", -1);
    declareIntBranch( m_hypMeths, "muon_roadUpstreamPlanes", -1);
    declareIntBranch( m_hypMeths, "muon_charge", -99);
    declareDoubleBranch( m_hypMeths, "muon_roadUpstreamEnergy", 0.0);
    declareDoubleBranch( m_hypMeths, "muon_E", 0.0);
    declareDoubleBranch( m_hypMeths, "muon_p", 0.0);
    declareDoubleBranch( m_hypMeths, "muon_px", 0.0);
    declareDoubleBranch( m_hypMeths, "muon_py", 0.0);
    declareDoubleBranch( m_hypMeths, "muon_pz", 0.0);  
    declareDoubleBranch( m_hypMeths, "muon_theta", 0.0);
    declareDoubleBranch( m_hypMeths, "muon_theta_biasUp", 0.0);
    declareDoubleBranch( m_hypMeths, "muon_theta_biasDown", 0.0); 
    declareDoubleBranch( m_hypMeths, "muon_muScore", -1.0);
    declareDoubleBranch( m_hypMeths, "muon_qp", 99.0);
    declareDoubleBranch( m_hypMeths, "muon_qpqpe", 99.0);
    declareDoubleBranch( m_hypMeths, "muon_E_shift", 0.0);
    
    //! NeutrinoInt - Proton -- Filled in setProtonParticleData()
    declareContainerIntBranch(m_hypMeths,    "proton_trk_pat_history", 10, -1);
    declareContainerIntBranch(m_hypMeths,    "proton_kinked",      10, -1);
    declareContainerIntBranch(m_hypMeths,    "proton_odMatch",     10, -1);
    declareContainerDoubleBranch(m_hypMeths, "proton_startPointX", 10, -1);
    declareContainerDoubleBranch(m_hypMeths, "proton_startPointY", 10, -1);
    declareContainerDoubleBranch(m_hypMeths, "proton_startPointZ", 10, -1);
    declareContainerDoubleBranch(m_hypMeths, "proton_endPointX",   10, -1);
    declareContainerDoubleBranch(m_hypMeths, "proton_endPointY",   10, -1);
    declareContainerDoubleBranch(m_hypMeths, "proton_endPointZ",   10, -1);
    declareContainerDoubleBranch(m_hypMeths, "proton_score",    10, -1);
    declareContainerDoubleBranch(m_hypMeths, "proton_score1",   10, -1);
    declareContainerDoubleBranch(m_hypMeths, "proton_score2",   10, -1);
    declareContainerDoubleBranch(m_hypMeths, "proton_chi2_ndf", 10, -1);
    declareContainerDoubleBranch(m_hypMeths, "proton_theta",    10, -1);
    declareContainerDoubleBranch(m_hypMeths, "proton_thetaX",   10, -1);
    declareContainerDoubleBranch(m_hypMeths, "proton_thetaY",   10, -1);
    declareContainerDoubleBranch(m_hypMeths, "proton_phi",      10, -1);
    declareContainerDoubleBranch(m_hypMeths, "proton_ekin",     10, -1);
    declareContainerDoubleBranch(m_hypMeths, "proton_E",        10, -1);
    declareContainerDoubleBranch(m_hypMeths, "proton_p",        10, -1);
    declareContainerDoubleBranch(m_hypMeths, "proton_px",       10, -1);
    declareContainerDoubleBranch(m_hypMeths, "proton_py",       10, -1);
    declareContainerDoubleBranch(m_hypMeths, "proton_pz",       10, -1);
    declareContainerDoubleBranch(m_hypMeths, "proton_p_calCorrection", 10, -1);
    declareContainerDoubleBranch(m_hypMeths, "proton_p_visEnergy",     10, -1);
    declareContainerDoubleBranch(m_hypMeths, "proton_p_dEdXTool",      10, -1);
    
   
    
    info() <<"Exit CCProtonPi0::initialize()" << endmsg;
    info() <<"--------------------------------------------------------------------------"<<endmsg;
    
    return sc;
}
    
//==============================================================================
//
// reconstructEvent() --
//
//==============================================================================
StatusCode CCProtonPi0::reconstructEvent( Minerva::PhysicsEvent *event, Minerva::GenMinInteraction* truthEvent ) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() <<"Enter CCProtonPi0::reconstructEvent()" << endmsg;
    
    //--------------------------------------------------------------------------
    //! Initialize truthEvent reco booleans
    //--------------------------------------------------------------------------
    if (truthEvent) {
        truthEvent->filtertaglist()->setOrAddFilterTag( "reco_hasGoodObjects", false );
        truthEvent->filtertaglist()->setOrAddFilterTag( "reco_isGoodVertex", false );
        truthEvent->filtertaglist()->setOrAddFilterTag( "reco_isWellFitVertex", false );
        truthEvent->filtertaglist()->setOrAddFilterTag( "reco_isFidVol", false );
        truthEvent->filtertaglist()->setOrAddFilterTag( "reco_isMinosMatch", false );
        truthEvent->filtertaglist()->setOrAddFilterTag( "reco_isBrokenTrack", false );
        truthEvent->filtertaglist()->setOrAddFilterTag( "reco_isFidVol_smeared", false );
    }
    
    //--------------------------------------------------------------------------
    //! Initialize reco booleans
    //--------------------------------------------------------------------------
    event->filtertaglist()->setOrAddFilterTag( "well_fit_vertex", false );
    event->filtertaglist()->setOrAddFilterTag( "isMinosMatchTrack", false ); 
    event->filtertaglist()->setOrAddFilterTag( "isMinosMatchStub", false );
    event->filtertaglist()->setOrAddFilterTag( "isBrokenTrack", false ); 
    
    if( truthEvent ){
        info() << "This is a MC event." << endmsg;
    }
    
    //--------------------------------------------------------------------------
    //! Check if this a plausible event ( MC only )
    //--------------------------------------------------------------------------
    // default plausibility means that underlying truth is likely something to be rejected by standard Minos Match CC analysis
    if( truthEvent && m_doPlausibilityCuts && !truthIsPlausible(truthEvent) ) {
        debug() << "This is not a plausible MC event! returning!" << endmsg;
        return StatusCode::SUCCESS;
    }
    
    
    //--------------------------------------------------------------------------
    //! Does this event carry a reco object with a BadObject flag?
    //--------------------------------------------------------------------------
    if( event->filtertaglist()->isFilterTagTrue( AnaFilterTags::BadObject() ) ) { 
        error() << "Found an event flagged with a BadObject! Refusing to analyze..." << endmsg;
        counter("REFUSED_A_BADOBJECT") += 1;
        return StatusCode::SUCCESS; // Things are bad, but we didn't crash.
    }
    counter("REFUSED_A_BADOBJECT") += 0;
    if (truthEvent){
        truthEvent->filtertaglist()->setOrAddFilterTag( "reco_hasGoodObjects", true); 
    }
    
    //==========================================================================
    //
    // Vertex Reconstruction
    //
    //==========================================================================
    debug() << "START: Vertex Reconstruction..." << endmsg;
    
    //--------------------------------------------------------------------------
    //! MAKE CUT - if NO Interaction Vertex
    //--------------------------------------------------------------------------
    if( !(event->hasInteractionVertex()) ) {
        debug() << "The event does not have an interaction vertex!" << endmsg;
        event->setIntData("Cut_Vertex_None",1);
        if( m_store_all_events ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    } 
       
    //--------------------------------------------------------------------------
    //! MAKE CUT - if NULL Interaction Vertex
    //--------------------------------------------------------------------------
    

    if( event->interactionVertex() == NULL ) { 
        bool pass = true; 
        std::string tag = "BadObject";
        event->filtertaglist()->addFilterTag(tag,pass);
        event->setIntData("Cut_Vertex_Null",1);
        error() << "This vertex is NULL! Flag this event as bad!" << endmsg;
        if( m_store_all_events ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    }
    SmartRef<Minerva::Vertex> m_PrimaryVertex = event->interactionVertex();
    
    //--------------------------------------------------------------------------
    //! MAKE CUT - if Interaction Vertex is NOT in Reconstructable Volume
    //--------------------------------------------------------------------------
    // "Vertex is NOT in Reconstructable Volume" means we can not run vertex-anchored
    // short tracker with a meaningful result outside of that volume.
    // Only the events that pass this CUT are used in vertex-anchored short tracker 
    if ( !FiducialPointTool->isFiducial(m_PrimaryVertex->position(), m_recoHexApothem, m_recoUpStreamZ, m_recoDownStreamZ ) ) {
        Gaudi::XYZPoint vtx_position = m_PrimaryVertex->position();
        debug() <<"Interaction Vertex is NOT in reconstructable volume = ("<<vtx_position.x()<<","<<vtx_position.y()<<","<<vtx_position.z()<<")"<< endmsg;
        event->setIntData("Cut_Vertex_Not_Reconstructable",1);
        if( m_store_all_events ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    }
    
    //--------------------------------------------------------------------------
    //! Create vertex-anchored short track Prongs, refit vertex
    //--------------------------------------------------------------------------
    int n_anchored_long_trk_prongs = event->primaryProngs().size() - 1;
    if (m_makeShortTracks) {
        debug()<<"Making short tracks"<<endmsg;
        m_ccPionIncUtils->makeShortTracks(event); 
        debug()<<"Finished making short tracks"<<endmsg;
    }  
    int n_anchored_short_trk_prongs = event->primaryProngs().size() - n_anchored_long_trk_prongs - 1;
    int n_iso_trk_prongs = (event->select<Prong>("Used:Unused","All")).size() - event->primaryProngs().size();
    
    Minerva::TrackVect all_tracks = event->select<Track>("Used:Unused","All");
    for (Minerva::TrackVect::iterator itTrk = all_tracks.begin(); itTrk != all_tracks.end(); ++itTrk) {
        Minerva::Track::NodeContainer nodes = (*itTrk)->nodes();
        Minerva::Track::NodeContainer::const_iterator node = nodes.begin();
        for (  ;node != nodes.end() -1; ++node) {
            if ( (*node)->z() == (*(node+1))->z()) {
                warning()<<"Duplicate node at "<<(*node)->z()<<" on track with type: "<<(*itTrk)->patRecHistory()<<endmsg;
                warning()<<"  Energies: "<<(*node)->idcluster()->energy()<<" "<<(*(node+1))->idcluster()->energy()<<endmsg;
            }
        }
    }
    
    //--------------------------------------------------------------------------
    //! MAKE CUT - if Interaction Vertex is NOT in Fiducial Volume
    //-------------------------------------------------------------------------- 
    debug()<<"Making fiducial volume cut"<<endmsg;
    
    // check for NaN first
    Gaudi::XYZPoint vtx_position = m_PrimaryVertex->position();
    if (vtx_position.z() != vtx_position.z()) {
        warning()<<"NaN vertex!"<<endmsg;
        if( m_store_all_events ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    }
    
    //smear the vertex position by 1mm 1mm 10mm, then record whether the smeared vertex passes the fiducial volume cut
    if (truthEvent) {
        Gaudi::XYZPoint vtx_smear(  vtx_position.x() + m_randomGen->Gaus(0.0,0.91*CLHEP::mm), 
                                    vtx_position.y() + m_randomGen->Gaus(0.0,1.25*CLHEP::mm), 
                                    vtx_position.z() + m_randomGen->Gaus(0.0, 10.0*CLHEP::mm));
        if ( FiducialPointTool->isFiducial(vtx_smear, m_fidHexApothem, m_fidUpStreamZ, m_fidDownStreamZ ) ) {
            truthEvent->filtertaglist()->setOrAddFilterTag( "reco_isFidVol_smeared", true );
        }
        else{
            truthEvent->filtertaglist()->setOrAddFilterTag( "reco_isFidVol_smeared", false );
        }
    }  
    
    if( !FiducialPointTool->isFiducial( vtx_position, m_fidHexApothem, m_fidUpStreamZ, m_fidDownStreamZ ) ){
        debug() <<"Interaction Vertex is not in fiducial volume = ("<<vtx_position.x()<<","<<vtx_position.y()<<","<<vtx_position.z()<<")"<< endmsg;
        event->setIntData("Cut_Vertex_Not_Fiducial",1);
        if( m_store_all_events ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    }
    
    //--------------------------------------------------------------------------
    //! VERTEX Passed all CUTS - Extend Vertex Reconstruction
    //--------------------------------------------------------------------------
    // Get Vertex Count
    SmartRefVector<Minerva::Vertex> allVertices  = event->select<Minerva::Vertex>( "All","StartPoint" );
    debug()<<"N(vertices) = "<<allVertices.size()<<endmsg;
    event->setIntData("vtx_total_count",allVertices.size());
    
    SmartRefVector<Minerva::Vertex>::iterator iter_primary_vtx = allVertices.end();
    SmartRefVector<Minerva::Vertex>::iterator iter_vtx;
    for ( iter_vtx = allVertices.begin(); iter_vtx != allVertices.end(); ++iter_vtx ){
        if (*iter_vtx == m_PrimaryVertex) {
            iter_primary_vtx = iter_vtx;
            break;
        }
    }
    event->setIntData("vtx_primary_index", std::distance(allVertices.begin(),iter_primary_vtx));
    if (iter_primary_vtx != allVertices.end()) allVertices.erase(iter_primary_vtx);
    event->setIntData("vtx_secondary_count", allVertices.size()); /* Number of secondary vertices */
    
    debug() << "FINISH: Vertex Reconstruction!" << endmsg;

    //==========================================================================
    //
    // Muon Reconstruction
    //
    //==========================================================================

    debug() << "START: Muon Reconstruction..." << endmsg;

    bool is_minos_track = false, is_minos_stub = false;
    bool foundMuon = MuonUtils->findMuonProng( event, m_MuonProng, m_MuonParticle );
    debug() << "m_MuonProng and m_MuonParticle is Saved!" << endmsg;
    
    //--------------------------------------------------------------------------
    //! MAKE CUT - if NO GOOD MUON (MINOS Matched)
    //--------------------------------------------------------------------------  
    if( !foundMuon ){
        debug() << "Did not find a muon prong!" << endmsg;
        event->setIntData("Cut_Muon_None",1);
        if( m_store_all_events ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    }
    
    if ( !m_MuonProng ) { 
        warning() << "Identified a muon Prong, but it is NULL!" << endmsg;
        return StatusCode::FAILURE; // We sort of did crash... 
    }
    
    //--------------------------------------------------------------------------
    //! Check if this a plausible Muon ( MC only )
    //--------------------------------------------------------------------------
    double mc_frac = -1.0;
    if ( m_doPlausibilityCuts && !muonIsPlausible( m_MuonProng, mc_frac) ) {
        debug()<<"Muon is not plausible"<<endmsg;
        if( m_store_all_events ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    }
    
    //--------------------------------------------------------------------------
    //! MAKE CUT - if MUON Score is LOW
    //--------------------------------------------------------------------------  
    debug() << "Muon Particle Score: " << m_MuonParticle->score() << endmsg;
    if(m_MuonParticle->score() < m_minMuonScore){
        debug()<<"Muon prong does not pass score cut"<<endmsg;
        event->setIntData("Cut_Muon_Score_Low",1);
        if( m_store_all_events ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    }
    
    //--------------------------------------------------------------------------
    //! MAKE CUT - if Muon has positive charge (AntiMuon)
    //--------------------------------------------------------------------------
    int charge = -99;
    MuonUtils->muonCharge(m_MuonProng,charge);
    if(charge == 1){
        debug()<<"AntiMuon Contamination"<<endmsg;
        event->setIntData("Cut_Muon_Charge",1);
        if( m_store_all_events ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    }
        
    //--------------------------------------------------------------------------
    //! MUON Passed All Cuts -- tag it as "PrimaryMuon"
    //--------------------------------------------------------------------------
    m_MuonProng->filtertaglist()->setOrAddFilterTag( "PrimaryMuon", true );
    m_MuonParticle->filtertaglist()->setOrAddFilterTag( "PrimaryMuon", true );
    
    // Save Muon Track info
    SmartRef<Minerva::Track> muonTrack = m_MuonProng->minervaTracks().front();
    m_MuonTrack = muonTrack;
    debug() << "m_MuonTrack is Saved!" << endmsg;
     
    if (m_MuonProng->MinosTrack()) is_minos_track = true;
    if (m_MuonProng->MinosStub()) is_minos_stub = true;
    if (is_minos_stub && is_minos_track) counter("MuonHasMinosStubAndTrack")++;
    else counter("MuonHasMinosStubAndTrack")+=0;
    if (!is_minos_stub && !is_minos_track) counter("MuonIsNotMinosMatched")++;
    else counter("MuonIsNotMinosMatched")+=0;
    
    event->filtertaglist()->setOrAddFilterTag("isMinosMatchTrack", is_minos_track );
    event->filtertaglist()->setOrAddFilterTag("isMinosMatchStub", is_minos_stub );
    if (truthEvent) truthEvent->filtertaglist()->setOrAddFilterTag( "reco_isMinosMatch", true );
    

    // Get Muon Energy
    double muon_visible_energy = m_MuonProng->minervaVisibleEnergySum();
    m_hitTagger->applyColorTag(m_MuonProng, m_muonProngColor);
    event->setTime( m_recoTimeTool->prongBestTime(m_MuonProng) );
    SmartRefVector<Track> muonTracks = m_MuonProng->minervaTracks(); 
    
    debug() << "FINISH: Muon Reconstruction" << endmsg;
    
    //==========================================================================
    //
    // Michel Electrons
    //
    //==========================================================================
    
    //--------------------------------------------------------------------------
    //! MAKE CUT - Vertex Michels
    //--------------------------------------------------------------------------
    debug() << "START: Vertex Michels" << endmsg;
    Minerva::Prong vtx_michel_prong;
    bool foundMichel = m_michelVtxTool->findMichel( m_PrimaryVertex, vtx_michel_prong );
    if (foundMichel) {
        debug()<<"Found a Michel Electron!"<<endmsg;
        event->setIntData("Cut_Vertex_Michel_Exist",1);
        event->setIntData("n_vtx_michel_views",vtx_michel_prong.getIntData("category"));
        event->setDoubleData("vtx_michel_distance",vtx_michel_prong.getDoubleData("distance"));
        if( m_store_all_events ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    }else{
        debug()<<"There are NO Vertex Michel Electrons in the event!"<<endmsg;
    }
    
    debug() << "FINISH: Vertex Michel" << endmsg;
    
    //--------------------------------------------------------------------------
    //! MAKE CUT - End Point Michels
    //--------------------------------------------------------------------------
    debug()<<"START: End Point Michel"<<endmsg;
  
    // Get All Prongs to search for end point michels
    ProngVect primaryProngs = event->primaryProngs();

    for (ProngVect::iterator iterProngs = primaryProngs.begin(); iterProngs != primaryProngs.end(); ++iterProngs) {
    
        // Skip Muon Prong
        if ( (*iterProngs) == m_MuonProng ) {
            debug()<<"Muon prong skipping!"<<endmsg;
            continue;
        }
            
        //Make an endpoint vertex blob for the hadron prong - currently this depends on the prong order!!!
        debug()<<"Make endpoint blob"<<endmsg;
        SmartRef<Minerva::Vertex> endpoint_vtx;
        m_objectAssociator->getVertex_fromTrackBack( endpoint_vtx, (*iterProngs)->minervaTracks().back() );
        if (!endpoint_vtx) {
            warning()<<"Could not find a back vertex for this prong!"<<endmsg;
            continue;
        }        
        
        IDClusterVect StopPointClusters = event->select<Minerva::IDCluster>("Unused","!XTalkCandidate");
        if (!StopPointClusters.empty()) {
            Minerva::IDBlob* StopPointBlob = new Minerva::IDBlob;
            debug()<<"CreateIDBlobs"<<endmsg;
            m_stopPointBlobTool->createIDBlobs(StopPointClusters, StopPointBlob, endpoint_vtx);
            m_hitTagger->applyColorTag(StopPointBlob, m_endPointVertexColor);
            addObject(event, StopPointBlob);        
            endpoint_vtx->addIDBlob(StopPointBlob);
            (*iterProngs)->add(StopPointBlob);
            (*iterProngs)->setDoubleData("endpointE", m_caloUtils->applyCalConsts(StopPointBlob, "Default") );
        }
        
        //Search for endpoint michels
        Minerva::Prong michelProng;
        bool foundMichel = m_michelTrkTool->findMichel( endpoint_vtx, michelProng );

        if (foundMichel){
            debug()<<"Found an End Point Michel Electron!"<<endmsg;
            event->setIntData("Cut_EndPoint_Michel_Exist",1);
            if( m_store_all_events ) return interpretFailEvent(event); 
            else return StatusCode::SUCCESS; 
        }
        
        // Search for secondary Michels
        if ( (*iterProngs)->minervaTracks().size() > 1) {
            debug()<<"Searching for secondary end point michels..."<<endmsg;
            for (SmartRefVector<Minerva::Track>::iterator iterTracks = (*iterProngs)->minervaTracks().begin(); iterTracks != (*iterProngs)->minervaTracks().end() - 1; ++iterTracks) {
                SmartRef<Minerva::Vertex> current_end_vtx;
                m_objectAssociator->getVertex_fromTrackBack( current_end_vtx, *iterTracks );
                foundMichel = m_michelTrkTool->findMichel( current_end_vtx, michelProng);
                if (foundMichel) {
                    debug()<<"Found a Secondary End Point Michel Electron!"<<endmsg;
                    event->setIntData("Cut_secEndPoint_Michel_Exist",1);
                    if( m_store_all_events ) return interpretFailEvent(event); 
                    else return StatusCode::SUCCESS; 
                }
            }
        }
    }//end loop over primary prongs
    
    debug()<<"FINISH: End Point Michel"<<endmsg;
    
    //==========================================================================
    //
    // Proton Reconstruction
    //
    //==========================================================================
    
    //--------------------------------------------------------------------------
    //! MAKE CUT - if NO GOOD Proton
    //--------------------------------------------------------------------------
    debug() << "START: Proton Reconstruction" << endmsg;
    //-- get all of the primary prongs in the event
    primaryProngs = event->primaryProngs();
    
    //-- create new particles
    debug() << "Creating particles with Proton and Pion hypotheses" <<endmsg;
    bool makeParticles = createTrackedParticles(primaryProngs);
    if (makeParticles){
        debug() << "Creation of Particles are SUCCESFUL!"<< endmsg;
    }else{
        warning() << "Creation of Particles are FAILED!"<< endmsg;
    }
    
    //-- check if one of the primary prong is contained and has a dEdX proton particle 
    Minerva::ProngVect protonProngs;
    Minerva::ParticleVect protonParticles;
    bool foundProton = getProtonProng(primaryProngs, protonProngs, protonParticles);
    if( foundProton ) {
        for(unsigned int i = 0; i < protonProngs.size(); i++) {
            debug() << "Tag the proton's prong with bit-field = " << protonProngs[i]->typeBitsToString() << endmsg;
        }
    } 
    else {
        debug() << "Didn't find any contained in the tracker bit-positive prong with a proton particle!" << endmsg;
        event->setIntData("Cut_Proton_None",1);
        if( m_store_all_events ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    }
    
    debug()<<"FINISH: Proton Reconstruction"<<endmsg;
    
    
    //==========================================================================
    //
    // Pi0 Reconstruction
    //
    //==========================================================================
    
    //--------------------------------------------------------------------------
    //! MAKE CUT - if fails PreFilterPi0()
    //--------------------------------------------------------------------------
    if ( PreFilterPi0(event) ){
        debug()<<"Passed PreFilterPi0()"<<endmsg;
    }else{
        event->setIntData("Cut_PreFilter_Pi0",1);
        if( m_store_all_events ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS;  
    }
    
    //--------------------------------------------------------------------------
    //! MAKE CUT - if fails VtxBlob()
    //--------------------------------------------------------------------------
    debug() << "START: Vertex Blob Reconstruction" << endmsg;
    if ( VtxBlob(event, m_PrimaryVertex) ){
        debug()<<"Succesful VtxBlob Reconstruction!"<<endmsg;
    }else{
        if( m_store_all_events ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS;  
    }
    
    debug()<<"FINISH: Vertex Blob Reconstruction"<<endmsg;
    
    //--------------------------------------------------------------------------
    //! MAKE CUT - if fails ConeBlobs()
    //--------------------------------------------------------------------------
    debug() << "START: Cone Blob Reconstruction" << endmsg;
    if ( ConeBlobs(event, m_PrimaryVertex) ){
        debug()<<"Succesful ConeBlobs Reconstruction!"<<endmsg;
    }else{
        if( m_store_all_events ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS;  
    }
    
    debug()<<"FINISH: Cone Blob Reconstruction"<<endmsg;
    
    
    //--------------------------------------------------------------------------
    //! @todo Determine if vertex has broken track
    //--------------------------------------------------------------------------
//     unsigned int broken_US_plane = 0;  
//     bool has_broken_track = vertexHasBrokenTrack(m_PrimaryVertex, broken_US_plane);
//     event->filtertaglist()->setOrAddFilterTag("isBrokenTrack",has_broken_track);  
//     if (truthEvent) truthEvent->filtertaglist()->setOrAddFilterTag( "reco_isBrokenTrack", has_broken_track);
//     
//     unsigned int upstream_plane_num = getMostUpstreamPlane( m_PrimaryVertex );
//     event->setIntData("upstream_plane_num", upstream_plane_num);
    
    
    //--------------------------------------------------------------
    //! Get total visible energy in PhysicsEvent
    //--------------------------------------------------------------
    SmartRefVector<IDCluster> idClusters = event->select<IDCluster>("Used:Unused","!XTalkCandidate");
    SmartRefVector<ODCluster> odClusters = event->select<ODCluster>("Used:Unused","!XTalkCandidate");  
    double odVisibleEnergy = m_extraEnergyTool->getODEnergy(odClusters, -1.0);
    double idVisibleEnergy = m_extraEnergyTool->getIDEnergy(idClusters, -1.0);
    double totalVisibleEnergy = odVisibleEnergy + idVisibleEnergy;  


    //--------------------------------------------------------------
    //! count other reconstructed quantities
    //--------------------------------------------------------------
    int n_long_tracks = (event->select<Track>("Used:Unused","LongPatRec3View:LongPatRec2View")).size();
    int n_short_tracks = (event->select<Track>("Used:Unused","FourHitPatRec")).size();
    
    //--------------------------------------------------------------
    //! Calculate dead time
    //--------------------------------------------------------------
    m_getDeadTimeTool->setDeadTimeTool( getDeadTimeTable() );
    
    // get dead time data
    std::map<const char*, int> intDeadData;
    std::map<const char*, double> doubleDeadData;
    DeadTimeTool->totalDead( event, intDeadData, doubleDeadData );
    
    debug() << "Size of dead time data: " << intDeadData.size() << endmsg;
    int dead = 0, udead = 0, ddead = 0, tdead = 0;
    for( std::map<const char*, int>::iterator it_map = intDeadData.begin(); it_map != intDeadData.end(); it_map++ ) {
        debug() << "  " << (*it_map).first << ": " << (*it_map).second << endmsg;
        if((std::string)(*it_map).first=="UDEAD")
        udead += (*it_map).second;
        if((std::string)(*it_map).first=="DDEAD")
        ddead += (*it_map).second;
        if((std::string)(*it_map).first=="TDEAD")
        tdead += (*it_map).second;
        if((std::string)(*it_map).first=="DEAD")
        dead += (*it_map).second; 
    }
    
    //--------------------------------------------------------------
    //! event time
    //--------------------------------------------------------------  
    double event_time = event->time();

    //--------------------------------------------------------------------------
    //! Finish filling event portion of ntuple 
    //--------------------------------------------------------------------------
    event->setDoubleData("time", event_time);
    
    event->setIntData( "dead", dead );
    event->setIntData( "udead", udead );
    event->setIntData( "ddead", ddead );
    event->setIntData( "tdead", tdead );
    
//     event->setIntData("n_startpoint_vertices", n_startpoint_vertices);
    event->setIntData("n_long_tracks", n_long_tracks);
    event->setIntData("n_short_tracks", n_short_tracks);
    event->setIntData("n_anchored_long_trk_prongs", n_anchored_long_trk_prongs);
    event->setIntData("n_anchored_short_trk_prongs", n_anchored_short_trk_prongs);
//     event->setIntData("n_vtx_blob_prongs", n_vtx_blob_prongs);
    event->setIntData("n_iso_trk_prongs", n_iso_trk_prongs);
//     event->setIntData("n_iso_blob_prongs", n_iso_blob_prongs);
//     event->setIntData("n_dsp_blob_prongs", n_dsp_blob_prongs);
//     
    event->setDoubleData("muonVisibleE", muon_visible_energy );
//     event->setDoubleData("hadronVisibleE", hadron_visible_energy );
    
    event->setDoubleData( "totalVisibleE",   totalVisibleEnergy );
    event->setDoubleData( "totalIDVisibleE", idVisibleEnergy );
    event->setDoubleData( "totalODVisibleE", odVisibleEnergy );
//     //event->setDoubleData( "vtxTrackExtraE",  vtx_track_extra_E);
//     event->setDoubleData( "vtxBlobExtraE",   vtx_blob_energy );
//     event->setDoubleData( "unattachedExtraE",unattached_extra_energy );
//     event->setDoubleData( "dispersedExtraE", dispersed_energy );  
    

    fillCommonPhysicsAnaBranches( event );

    //--------------------------------------------------------------------------
    //! Call the interpretEvent function.
    //--------------------------------------------------------------------------
    NeutrinoVect interactions;
    StatusCode interpret = this->interpretEvent( event, truthEvent, interactions );
    
    //! If there were any neutrino interactions reconstructed, mark the event.
    if( interactions.size() == 1 ){
        markEvent( event );
    }
    else {
        warning()<<"interpretEvent() returns "<<interactions.size()<<" interactions!"<<endmsg;
        return StatusCode::SUCCESS; // We didn't crash.
    }
    
    //! Add the interactions to the event. Use MinervaHistoTool::addInteractionHyp.
    StatusCode sc = addInteractionHyp( event, interactions );
    

    debug() <<"Exit CCProtonPi0::reconstructEvent()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    return sc;
    
}
    
//==============================================================================
//
// interpretEvent()
//
//==============================================================================
StatusCode CCProtonPi0::interpretEvent( const Minerva::PhysicsEvent *event, const Minerva::GenMinInteraction *truthEvent, NeutrinoVect& interaction_hyp ) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() <<"Enter CCProtonPi0::interpretEvent()" << endmsg;
    
    if( truthEvent ){
        debug() << "This is a MC event." << endmsg;
    }
    
    if( !event ){
        debug() << "NULL Event" << endmsg;
        return StatusCode::FAILURE;
    }
    
    //--------------------------------------------------------------------------
    //! Create interaction hypothesis
    //--------------------------------------------------------------------------
    Minerva::NeutrinoInt* nuInt = new Minerva::NeutrinoInt(m_anaSignature);
    interaction_hyp.push_back( nuInt );

    //--------------------------------------------------------------------------
    //! Identify and Store - Muon and Proton Prongs
    //--------------------------------------------------------------------------
    ProngVect primaryProngs = event->primaryProngs();
    
    SmartRef<Prong> muonProng = (Prong*)NULL;
    ProngVect protonProngs;
    ParticleVect protonParticles;
    
    for (ProngVect::iterator itProng = primaryProngs.begin(); itProng != primaryProngs.end(); ++itProng) {

        bool isPrimaryMuon = false;
        bool isPrimaryProton = false;
    
        if ( (*itProng)->filtertaglist()->filterTagExists("PrimaryMuon") ) {
            (*itProng)->filtertaglist()->checkFilterTag( "PrimaryMuon", isPrimaryMuon );
        }
        
        if ( (*itProng)->filtertaglist()->filterTagExists("PrimaryProton") ) {
            (*itProng)->filtertaglist()->checkFilterTag( "PrimaryProton", isPrimaryProton );
        }    
 
        if (isPrimaryMuon && !isPrimaryProton){
            debug()<<"muonProng Identified!"<<endmsg;
            muonProng = *itProng;
        }
        
        if (isPrimaryProton && !isPrimaryMuon){
            debug()<<"protonProng Identified!"<<endmsg;
            protonProngs.push_back(*itProng);
            protonParticles.push_back((*itProng)->bestParticle());
        }
        
        if (isPrimaryMuon && isPrimaryProton ) {
            warning()<<"Prong is two primary particles!"<<endmsg;
        }
    }
    
    
    
    //--------------------------------------------------------------
    //! Get MCTracks vector
    //--------------------------------------------------------------
//     std::vector<MCTrack> MCTracksVector;
//     if (truthEvent) {
//         MCTracksVector = m_MCTrackTool->getMCTracks(truthEvent);
//     }
    
    //--------------------------------------------------------------------------  
    //! Calculate and Set Muon Kinematics
    //--------------------------------------------------------------------------
    setMuonParticleData( nuInt, m_MuonProng);

    //--------------------------------------------------------------------------  
    //! Calculate and Set Proton Kinematics
    //--------------------------------------------------------------------------
    
    // VertexZ Required for setProtonParticleData()
    double vertexZ = event->interactionVertex()->position().z();
    
    setProtonParticleData(nuInt,protonProngs,protonParticles,vertexZ);
    

    //--------------------------------------------------------------------------
    //! Calculate and Set Vertex Parameters
    //--------------------------------------------------------------------------
    Gaudi::XYZTVector vtx_position( event->interactionVertex()->position().x(), 
                                    event->interactionVertex()->position().y(),
                                    event->interactionVertex()->position().z(),
                                    event->time() );     
    nuInt->setVertex( vtx_position );
    nuInt->setScore( 1.0 );
    
    // Set primaryVertexData
    int vtx_module, vtx_plane;
    debug()<<"Calling getNearestPlane, vtx is "<<vtx_position.z()<<endmsg;
    getNearestPlane(vtx_position.z(), vtx_module, vtx_plane); 
    nuInt->setIntData("vtx_module", vtx_module);
    nuInt->setIntData("vtx_plane", vtx_plane);
    nuInt->setDoubleData("vtx_x", vtx_position.x() );
    nuInt->setDoubleData("vtx_y", vtx_position.y() );
    nuInt->setDoubleData("vtx_z", vtx_position.z() );
    

    

    
    
    //--------------------------------------------------------------------------
    //! Interaction Parameters
    //--------------------------------------------------------------------------
    //nuInt->setNeutrinoHelicity( getHelicity( mu_charge ) );
    nuInt->setInteractionCurrent( Minerva::NeutrinoInt::ChargedCurrent );
    nuInt->setInteractionType( Minerva::NeutrinoInt::UnknownInt );
    
    //--------------------------------------------------------------------------
    //! Truth Matching for Muon and Proton Prongs
    //--------------------------------------------------------------------------
    if( haveNeutrinoMC() ) {
        Minerva::ProngVect muonProngs;
        muonProngs.push_back( muonProng );
        
        setTrackProngTruth(nuInt,muonProngs);
        debug() << "       set the muon prong track truth information" << endmsg;
    
        setTrackProngTruth(nuInt,protonProngs);
        debug() << "       set the proton prong track truth information" << endmsg;
    }
    

    debug() <<"Exit CCProtonPi0::interpretEvent()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    return StatusCode::SUCCESS;
    
}
   

//==============================================================================
// tagTruth()
//==============================================================================
StatusCode CCProtonPi0::tagTruth( Minerva::GenMinInteraction* truthEvent ) const 
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter: CCProtonPi0::tagTruth()" << endmsg;

    if (!truthEvent) {
        warning()<<"Passed a null truthEvent to tagTruth()!"<<endmsg;
        return StatusCode::SUCCESS;
    }
    
    //--------------------------------------------------------------------------
    //! Fill GENIE Weight Branches
    //--------------------------------------------------------------------------
    StatusCode sc = fillGenieWeightBranches( truthEvent );
    if (sc.isFailure() ) {
        warning()<<"Genie weight branch filling failed!"<<endmsg;
        return sc;
    }
    
    
    //--------------------------------------------------------------------------
    //! MC Plausibility
    //--------------------------------------------------------------------------
    bool isPlausible = truthIsPlausible(truthEvent);
    truthEvent->filtertaglist()->setOrAddFilterTag( "isPlausible", isPlausible); 
    
    
    //--------------------------------------------------------------------------
    //! is Fiducial volume
    //--------------------------------------------------------------------------
    bool isFidVol = FiducialPointTool->isFiducial(truthEvent, m_fidHexApothem, m_fidUpStreamZ, m_fidDownStreamZ );
    truthEvent->filtertaglist()->setOrAddFilterTag( "isFidVol", isFidVol );
    debug()<<"Finished determining true_isFidVol: "<<isFidVol<<endmsg;
    
    
    //--------------------------------------------------------------------------
    //! get TG4Trajectories
    //--------------------------------------------------------------------------
    const SmartRefVector<Minerva::TG4Trajectory> pri_trajectories = truthEvent->trajectories();
    SmartRefVector<Minerva::TG4Trajectory>::const_iterator it_mcpart;
    if (pri_trajectories.size() != truthEvent->fSpdg().size() && pri_trajectories.size() + 1 != truthEvent->fSpdg().size() ) {
        warning()<<"Number of GENIE FS particles "<<truthEvent->fSpdg().size()<<
                " does not match number of Geant4 primary trajectories "<<pri_trajectories.size()<<"!"<<endmsg;
                
        for ( SmartRefVector<Minerva::TG4Trajectory>::const_iterator itTraj = pri_trajectories.begin(); itTraj != pri_trajectories.end(); ++itTraj ) {
        warning()<<"    G4 particle: "<<(*itTraj)->GetPDGCode()<<endmsg;
        }
        for ( std::vector<int>::const_iterator itGENIE = truthEvent->fSpdg().begin(); itGENIE != truthEvent->fSpdg().end(); ++itGENIE ) {
        warning()<<"    GENIE part:  "<<(*itGENIE)<<endmsg;
        }
    }
      
    //--------------------------------------------------------------------------
    //! Print the List of FS Particles
    //--------------------------------------------------------------------------
    debug()<<"Final State Particle List"<<endmsg;
    debug()<<"Parent\t|\tParticle"<<endmsg;
    for (it_mcpart = pri_trajectories.begin(); it_mcpart != pri_trajectories.end(); ++it_mcpart) {
        debug()<<(*it_mcpart)->GetParentId()<<"\t|\t"<<(*it_mcpart)->GetPDGCode()<<endmsg;
    }
    
    //--------------------------------------------------------------------------
    //! Count the Number of FS Particles
    //--------------------------------------------------------------------------
    int N_proton    = 0;
    int N_neutron   = 0;
    int N_muminus   = 0;
    int N_muplus    = 0;
    int N_pi0       = 0;
    int N_piplus    = 0;
    int N_piminus   = 0;
    int N_deltaplus = 0;
    int N_gamma     = 0;
    int N_other     = 0;

    for (it_mcpart = pri_trajectories.begin(); it_mcpart != pri_trajectories.end(); ++it_mcpart) {
    
        if ( ( (*it_mcpart)->GetPDGCode() ) == 2212 ) N_proton++;
        else if( ( (*it_mcpart)->GetPDGCode() ) == 2112 ) N_neutron++;
        else if( ( (*it_mcpart)->GetPDGCode() ) == 13 ) N_muminus++;
        else if( ( (*it_mcpart)->GetPDGCode() ) == -13 ) N_muplus++;
        else if( ( (*it_mcpart)->GetPDGCode() ) == 111 ) N_pi0++;
        else if( ( (*it_mcpart)->GetPDGCode() ) == 211 ) N_piplus++;
        else if( ( (*it_mcpart)->GetPDGCode() ) == -211 ) N_piminus++;
        else if( ( (*it_mcpart)->GetPDGCode() ) == 2214 ) N_deltaplus++;
        else if( ( (*it_mcpart)->GetPDGCode() ) == 22 ) N_gamma++;
        else N_other++;
        
    }
    
    //--------------------------------------------------------------------------
    //! Find Signal -- CC Neutrino Interaction with FS Particles: muon, proton, pi0
    //--------------------------------------------------------------------------
    bool isSignal = false;

    int t_current = truthEvent->current();
    int t_neutrinoPDG = truthEvent->incoming();
    
    // CC Neutrino Interaction
    if ( t_current==1 && t_neutrinoPDG == 14 ) {
        // Atleast 1 proton, 1 pi0, 0 pi+, 0 pi-
        if(N_proton > 0 && N_pi0 == 1 && N_piplus == 0 && N_piminus == 0){
            isSignal = true;
        }
    }
    truthEvent->filtertaglist()->setOrAddFilterTag( "isSignal", isSignal );
   
    
    //------------------------------------------------------------
    //! Target Material
    //------------------------------------------------------------
    Gaudi::LorentzVector t_vtx = truthEvent->Vtx(); 
    int t_targetMaterial = -1; 
    
    G4ThreeVector G4vtx(t_vtx.x(), t_vtx.y(), t_vtx.z());
    G4ThreeVector G4vec(0, 0, 1); // The next function needs a direction vector as argument. Here's a dummy.
    G4Navigator* geantNav = new G4Navigator;
    geantNav->SetWorldVolume(m_gigaCnvSvc->world());
    G4VPhysicalVolume* physVol=geantNav->LocateGlobalPointAndSetup(G4vtx, &G4vec, false, true);
    G4Material* material=physVol->GetLogicalVolume()->GetMaterial(); 
    delete geantNav;
    
    if ( material->GetName() == "/dd/Materials/Minerva/PlasticScint" ) {
        t_targetMaterial = 0;
    }
    else if ( material->GetName() == "/dd/Materials/Minerva/PSTitaniumDioxide" ) {
        t_targetMaterial = 1;
    }
    else if ( material->GetName() == "/dd/Materials/Minerva/GreyEpoxy" ) {
        t_targetMaterial = 2;
    }
    else if ( material->GetName() == "/dd/Materials/Minerva/Lexan" ) {
        t_targetMaterial = 3;
    }
    else if ( material->GetName() == "/dd/Materials/Minerva/GreenFiber" ) {
        t_targetMaterial = 4;
    }
    else if ( material->GetName() == "/dd/Materials/Minerva/PureLead" ) {
        t_targetMaterial = 5;
    }
    else if ( material->GetName() == "/dd/Materials/Minerva/PureCarbon" ) {
        t_targetMaterial = 6;
    }  
    else if ( material->GetName() == "/dd/Materials/Minerva/PureAluminum" ) {
        t_targetMaterial = 7;
    }   
    else if ( material->GetName() == "/dd/Materials/Minerva/StainlessSteel" ) {
        t_targetMaterial = 8;
    }   
    else if ( material->GetName() == "/dd/Materials/Air" ) {
        t_targetMaterial = 9;
    }
    else t_targetMaterial = 10; 
    
    //--------------------------------------------------------------
    //! Vertex module and plane
    //--------------------------------------------------------------   
    debug()<<"Calling getNearestPlane, t_vtx is "<<t_vtx.z()<<endmsg;    
    int vertex_module;
    int vertex_plane;
    getNearestPlane(t_vtx.z(), vertex_module, vertex_plane);  
    plot1D(2*vertex_module + vertex_plane + 0.5,"nearest_plane","Nearest Plane to point in z",-40,260, 300);
    plot2D(t_vtx.z(), 2*vertex_module + vertex_plane +0.5,"nearest_plane_v_z","",4000.0,10000.0, 120, -40, 260, 300);
    

    
    //!Get all trajectories
    Minerva::TG4Trajectories* all_trajectories = NULL;
    if (exist<Minerva::TG4Trajectories>(Minerva::TG4TrajectoryLocation::Default)) {
        all_trajectories = get<Minerva::TG4Trajectories>(Minerva::TG4TrajectoryLocation::Default);
        debug() <<"There are "<<all_trajectories->size()<<" TG4Trajectories"<<endmsg; 
    }
    if (!all_trajectories) {
        warning()<<"Could not find trajectories in TES"<<endmsg;
    }
    
    //--------------------------------------------------------------------------
    //! Get Muon Kinematics
    //--------------------------------------------------------------------------
    Gaudi::LorentzVector t_mu4p = truthEvent->PrimFSLepton();
    double t_muon_px = t_mu4p.px();
    double t_muon_py = t_mu4p.py();
    double t_muon_pz = t_mu4p.pz();
    double t_muon_E = t_mu4p.E();
    double t_muon_theta=m_coordSysTool->thetaWRTBeam(t_mu4p);
    int t_muon_charge = 0;
    
    if (t_current==1 && t_neutrinoPDG==14) t_muon_charge = -1;
    else if ( t_current==1 && t_neutrinoPDG==-14) t_muon_charge = 1;
    
    //--------------------------------------------------------------------------
    //! Get Proton, Pi0 and Gamma Kinematics
    //--------------------------------------------------------------------------
    int nMaxParticles = 20;
    
    std::vector<double> t_proton_px;
    std::vector<double> t_proton_py;
    std::vector<double> t_proton_pz;
    std::vector<double> t_proton_E;
    std::vector<double> t_proton_vtx_x;
    std::vector<double> t_proton_vtx_y;
    std::vector<double> t_proton_vtx_z;
    std::vector<double> t_proton_theta;
    std::vector<int> t_proton_trackID;
    std::vector<int> t_proton_parentID;
    
    std::vector<double> t_pi0_px;
    std::vector<double> t_pi0_py;
    std::vector<double> t_pi0_pz;
    std::vector<double> t_pi0_E;
    std::vector<double> t_pi0_vtx_x;
    std::vector<double> t_pi0_vtx_y;
    std::vector<double> t_pi0_vtx_z;
    std::vector<double> t_pi0_theta;
    std::vector<int> t_pi0_trackID;
    std::vector<int> t_pi0_parentID;
    
    std::vector<double> t_gamma_px;
    std::vector<double> t_gamma_py;
    std::vector<double> t_gamma_pz;
    std::vector<double> t_gamma_E;
    std::vector<double> t_gamma_vtx_x;
    std::vector<double> t_gamma_vtx_y;
    std::vector<double> t_gamma_vtx_z;
    std::vector<double> t_gamma_theta;
    std::vector<int> t_gamma_trackID;
    std::vector<int> t_gamma_parentID;
    
    // Initialize Vectors
    for (int i = 0; i < nMaxParticles; i++){
        t_proton_px.push_back(-1.0);
        t_proton_py.push_back(-1.0);
        t_proton_pz.push_back(-1.0);
        t_proton_E.push_back(-1.0);
        t_proton_vtx_x.push_back(-1.0);
        t_proton_vtx_y.push_back(-1.0);
        t_proton_vtx_z.push_back(-1.0);
        t_proton_theta.push_back(-9.0);
        t_proton_trackID.push_back(-1);
        t_proton_parentID.push_back(-1);
        
        t_pi0_px.push_back(-1.0);
        t_pi0_py.push_back(-1.0);
        t_pi0_pz.push_back(-1.0);
        t_pi0_E.push_back(-1.0);
        t_pi0_vtx_x.push_back(-1.0);
        t_pi0_vtx_y.push_back(-1.0);
        t_pi0_vtx_z.push_back(-1.0);
        t_pi0_theta.push_back(-9.0);
        t_pi0_trackID.push_back(-1);
        t_pi0_parentID.push_back(-1);
        
        t_gamma_px.push_back(-1.0);
        t_gamma_py.push_back(-1.0);
        t_gamma_pz.push_back(-1.0);
        t_gamma_E.push_back(-1.0);
        t_gamma_vtx_x.push_back(-1.0);
        t_gamma_vtx_y.push_back(-1.0);
        t_gamma_vtx_z.push_back(-1.0);
        t_gamma_theta.push_back(-9.0);
        t_gamma_trackID.push_back(-1);
        t_gamma_parentID.push_back(-1);
    }
    
    int nPi0 = 0;
    int nProton = 0;
    int nGamma = 0;
    
    for (it_mcpart = pri_trajectories.begin(); it_mcpart != pri_trajectories.end(); ++it_mcpart) {
        
        Gaudi::LorentzVector temp_4p = (*it_mcpart)->GetInitialMomentum();
        Gaudi::LorentzVector temp_vtx = (*it_mcpart)->GetFinalPosition();
        int partPDG = (*it_mcpart)->GetPDGCode();
        
        if ( (partPDG == 2212) && nProton < nMaxParticles){
            t_proton_px[nProton] = temp_4p.px();
            t_proton_py[nProton] = temp_4p.py();
            t_proton_pz[nProton] = temp_4p.pz();
            t_proton_E[nProton]  = temp_4p.E();
            t_proton_vtx_x[nProton] = temp_vtx.x();
            t_proton_vtx_y[nProton] = temp_vtx.y();
            t_proton_vtx_z[nProton] = temp_vtx.z();
            t_proton_theta[nProton] = m_coordSysTool->thetaWRTBeam(temp_4p); 
            t_proton_trackID[nProton] = (*it_mcpart)->GetTrackId();
            nProton++;
        }else if ( (partPDG == 111) && nPi0 < nMaxParticles){
            t_pi0_px[nPi0] = temp_4p.px();
            t_pi0_py[nPi0] = temp_4p.py();
            t_pi0_pz[nPi0] = temp_4p.pz();
            t_pi0_E[nPi0]  = temp_4p.E();
            t_pi0_vtx_x[nPi0] = temp_vtx.x();
            t_pi0_vtx_y[nPi0] = temp_vtx.y();
            t_pi0_vtx_z[nPi0] = temp_vtx.z();
            t_pi0_theta[nPi0] = m_coordSysTool->thetaWRTBeam(temp_4p); 
            t_pi0_trackID[nPi0] = (*it_mcpart)->GetTrackId();
            nPi0++;
        }else if ( (partPDG == 22) && nGamma < nMaxParticles){
            t_gamma_px[nGamma] = temp_4p.px();
            t_gamma_py[nGamma] = temp_4p.py();
            t_gamma_pz[nGamma] = temp_4p.pz();
            t_gamma_E[nGamma]  = temp_4p.E();
            t_gamma_vtx_x[nGamma] = temp_vtx.x();
            t_gamma_vtx_y[nGamma] = temp_vtx.y();
            t_gamma_vtx_z[nGamma] = temp_vtx.z();
            t_gamma_theta[nGamma] = m_coordSysTool->thetaWRTBeam(temp_4p); 
            t_gamma_trackID[nGamma] = (*it_mcpart)->GetTrackId();
            t_gamma_parentID[nGamma] = (*it_mcpart)->GetParentId();
            nGamma++;
        }
    }
    
    //--------------------------------------------------------------------------
    //! fill the truthEvent info
    //--------------------------------------------------------------------------
    debug()<<"Filling ntuple for truth information!"<<endmsg;
    
    truthEvent->setIntData("vertex_module", vertex_module);
    truthEvent->setIntData("vertex_plane", vertex_plane);  
    truthEvent->setIntData("target_material", t_targetMaterial);
    
    truthEvent->setDoubleData("muon_px", t_muon_px);
    truthEvent->setDoubleData("muon_py", t_muon_py);
    truthEvent->setDoubleData("muon_pz", t_muon_pz);
    truthEvent->setDoubleData("muon_E",  t_muon_E);
    truthEvent->setDoubleData("muon_theta_wrtbeam",  t_muon_theta);
    truthEvent->setIntData("muon_charge", t_muon_charge);
    
    truthEvent->setContainerDoubleData("proton_px", t_proton_px);
    truthEvent->setContainerDoubleData("proton_py", t_proton_py);
    truthEvent->setContainerDoubleData("proton_pz", t_proton_pz);
    truthEvent->setContainerDoubleData("proton_E",  t_proton_E);
    truthEvent->setContainerDoubleData("proton_vtx_x", t_proton_vtx_x);
    truthEvent->setContainerDoubleData("proton_vtx_y", t_proton_vtx_y);
    truthEvent->setContainerDoubleData("proton_vtx_z", t_proton_vtx_z);
    truthEvent->setContainerDoubleData("proton_theta_wrtbeam",  t_proton_theta);
    truthEvent->setContainerIntData("proton_trackID", t_proton_trackID);
    truthEvent->setContainerIntData("proton_parentID", t_proton_parentID);
    
    truthEvent->setContainerDoubleData("pi0_px", t_pi0_px);
    truthEvent->setContainerDoubleData("pi0_py", t_pi0_py);
    truthEvent->setContainerDoubleData("pi0_pz", t_pi0_pz);
    truthEvent->setContainerDoubleData("pi0_E",  t_pi0_E);
    truthEvent->setContainerDoubleData("pi0_vtx_x", t_pi0_vtx_x);
    truthEvent->setContainerDoubleData("pi0_vtx_y", t_pi0_vtx_y);
    truthEvent->setContainerDoubleData("pi0_vtx_z", t_pi0_vtx_z);
    truthEvent->setContainerDoubleData("pi0_theta_wrtbeam",  t_pi0_theta);
    truthEvent->setContainerIntData("pi0_trackID", t_pi0_trackID);
    truthEvent->setContainerIntData("pi0_parentID", t_pi0_parentID);
    
    truthEvent->setContainerDoubleData("gamma_px", t_gamma_px);
    truthEvent->setContainerDoubleData("gamma_py", t_gamma_py);
    truthEvent->setContainerDoubleData("gamma_pz", t_gamma_pz);
    truthEvent->setContainerDoubleData("gamma_E",  t_gamma_E);
    truthEvent->setContainerDoubleData("gamma_vtx_x", t_gamma_vtx_x);
    truthEvent->setContainerDoubleData("gamma_vtx_y", t_gamma_vtx_y);
    truthEvent->setContainerDoubleData("gamma_vtx_z", t_gamma_vtx_z);
    truthEvent->setContainerDoubleData("gamma_theta_wrtbeam",  t_gamma_theta);
    truthEvent->setContainerIntData("gamma_trackID", t_gamma_trackID);
    truthEvent->setContainerIntData("gamma_parentID", t_gamma_parentID);
    
    truthEvent->setIntData("N_proton",  N_proton);
    truthEvent->setIntData("N_neutron",  N_neutron);
    truthEvent->setIntData("N_muminus",  N_muminus);
    truthEvent->setIntData("N_muplus",  N_muplus);
    truthEvent->setIntData("N_pi0",  N_pi0);
    truthEvent->setIntData("N_piplus",  N_piplus);
    truthEvent->setIntData("N_piminus",  N_piminus);
    truthEvent->setIntData("N_deltaplus",N_deltaplus);
    truthEvent->setIntData("N_gamma", N_gamma);
    truthEvent->setIntData("N_other", N_other );
    
    debug() <<"Exit CCProtonPi0::tagTruth()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    return StatusCode::SUCCESS;
    
}

//==============================================================================
//  Finalize
//==============================================================================
StatusCode CCProtonPi0::finalize()
{
    info() <<"--------------------------------------------------------------------------"<<endmsg;
    info() << "Enter: finalize()" << endmsg;
    
    
    // finalize the base class.
    StatusCode sc = this->MinervaAnalysisTool::finalize();
    if( sc.isFailure() ) return Error( "Failed to finalize!", sc );
    
    
    info() <<"Exit CCProtonPi0::CCProtonPi0::finalize()" << endmsg;
    info() <<"--------------------------------------------------------------------------"<<endmsg;
    return sc;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
//
//  Private Functions
//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
    
//----------------------------------------------------------------------------------------
// interpret Events which fails the reconstructor cuts
//----------------------------------------------------------------------------------------
StatusCode CCProtonPi0::interpretFailEvent( Minerva::PhysicsEvent* event ) const
{  

    debug() <<"Exit CCProtonPi0::reconstructEvent() through interpretFailEvent()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter CCProtonPi0::interpretFailEvent()" << endmsg;
    
    Minerva::NeutrinoInt* nuInt = new Minerva::NeutrinoInt(m_anaSignature);
    NeutrinoVect nuInts;
    nuInts.push_back( nuInt );
    markEvent(event);
    addInteractionHyp(event,nuInts);
    fillCommonPhysicsAnaBranches(event);
    fillNuMIBranches(event);
    
    debug() << "Exit CCProtonPi0::interpretFailEvent()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    return StatusCode::SUCCESS;
}

//---------------------------------------------------------------------------------
// set the track prong Geant4 truth information
//---------------------------------------------------------------------------------
void CCProtonPi0::setTrackProngTruth( Minerva::NeutrinoInt* neutrino, Minerva::ProngVect& prongs ) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter CCProtonPi0::setTrackProngTruth()" << endmsg;
    
    std::vector<const Minerva::TG4Trajectory*> parent;
    Minerva::TG4Trajectories* mc_trajectories = get<Minerva::TG4Trajectories>( "MC/TG4Trajectories" );
    for(Minerva::TG4Trajectories::iterator iter = mc_trajectories->begin(); iter != mc_trajectories->end(); ++iter) {
        if( (*iter)->GetParentId() == 0 ) parent.push_back( (*iter) );
    }
    
    std::vector< std::pair<const Minerva::TG4Trajectory*,double> > trajectories;
    std::vector<int> ntraj(10,-1);
    std::vector<int> nprim(10,-1);
    std::vector<int> pdg(10,-1);
    std::vector<int> inside_od(10,-1);
    std::vector<double> traj_Px(10,-1);
    std::vector<double> traj_Py(10,-1);
    std::vector<double> traj_Pz(10,-1);
    std::vector<double> traj_E(10,-1);
    std::vector<double> traj_P(10,-1);
    std::vector<double> traj_PSelf(10,-1);
    std::vector<double> traj_theta(10,-1);
    std::vector<double> traj_phi(10,-1);
    std::vector<double> end_trajMomentum(10,-1);
    std::vector<double> end_trajX(10,-1);
    std::vector<double> end_trajY(10,-1);
    std::vector<double> end_trajZ(10,-1);
    
    for(unsigned int p = 0; p < prongs.size(); p++) {
        SmartRef<Minerva::Prong> prong = prongs[p];
        Minerva::TrackVect tracks = prong->minervaTracks();
    
        for(unsigned int trk = 0; trk < tracks.size(); trk++) {
            double other_energy = 0.0;
            std::map<const Minerva::TG4Trajectory*,double> trajMap = TruthMatcher->getTG4Trajectories(tracks[trk],other_energy);
            std::map<const Minerva::TG4Trajectory*,double>::iterator it;
        
            if( trajMap.empty() ) continue;
        
            std::vector< std::pair<const Minerva::TG4Trajectory*,double> > parentVect;
            for(it = trajMap.begin(); it != trajMap.end(); it++) {
        
                const Minerva::TG4Trajectory* tmp = (*it).first;    
                double fraction = (*it).second;
        
                if( tmp->GetParentId() == 0 ){
                    parentVect.push_back( std::make_pair(tmp,fraction) );
                }else {
                    int  attempts = 0;
                    bool foundPrimary = false;
                    while( !foundPrimary ) {
                    
                        bool foundMatch = false;
                        Gaudi::XYZPoint start( tmp->GetInitialPosition().x(), tmp->GetInitialPosition().y(), tmp->GetInitialPosition().z() );
            
                        for(Minerva::TG4Trajectories::iterator tj = mc_trajectories->begin(); tj != mc_trajectories->end(); ++tj) {
                            Gaudi::XYZPoint end( (*tj)->GetFinalPosition().x(), (*tj)->GetFinalPosition().y(), (*tj)->GetFinalPosition().z() );
                            if( end == start ) { 
                                tmp = (*tj); foundMatch = true; 
                                break; 
                            }
                        }
            
                        if( !foundMatch ) break;
                        else attempts++;
            
                        if( foundMatch && attempts == 10 ) break;
                        if( tmp->GetParentId() == 0 ) { 
                            parentVect.push_back( std::make_pair(tmp,fraction) ); 
                            foundPrimary = true; 
                            break; 
                        }
                    } //! end of while loop
                } //! end of else
            } //! end loop over traj 
    
            std::vector< std::pair<const Minerva::TG4Trajectory*,double> > updateParentVect;
            for(unsigned int i = 0; i < parent.size(); i++) {
                double energy = 0.0;
                for(unsigned int j = 0; j < parentVect.size(); j++) {
                    if( parentVect[j].first->GetTrackId() == parent[i]->GetTrackId() ){
                        energy += parentVect[j].second; 
                    }
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
    
        if( !traj ) continue;
    
        Gaudi::LorentzVector traj_4p = traj->GetInitialMomentum();
        double PSelf = sqrt(traj_4p.Px()*traj_4p.Px() +
                            traj_4p.Py()*traj_4p.Py() +
                            traj_4p.Pz()*traj_4p.Pz());
        
        int primary = 0;
        if( traj->GetProcessName().find("Primary") != std::string::npos ) primary = 1; 
    
        ntraj[p]      = int(trajectories.size());
        nprim[p]      = primary;
        pdg[p]        = traj->GetPDGCode();
        traj_Px[p]    = traj_4p.Px();
        traj_Py[p]    = traj_4p.Py();
        traj_Pz[p]    = traj_4p.Pz();
        traj_E[p]     = traj_4p.E();
        traj_P[p]     = traj_4p.P();
        traj_PSelf[p] = PSelf;
        traj_theta[p] = traj_4p.Theta();
        traj_phi[p]   = traj_4p.Phi();
    
        const Minerva::TG4Trajectory* endTraj = trajectories[trajectories.size()-1].first;
        if( !endTraj ) continue;
    
        const Minerva::TG4TrajectoryPoints& points = endTraj->GetTrajectoryPoints();
    
        if( !points.empty() ) {
            Minerva::TG4TrajectoryPoint* endPoint = *points.rbegin();
            if( endPoint ) {
                double end_trajPx = endPoint->GetMomentum().px();
                double end_trajPy = endPoint->GetMomentum().py();
                double end_trajPz = endPoint->GetMomentum().pz();
                end_trajMomentum[p]  = sqrt(    end_trajPx*end_trajPx + 
                                                end_trajPy*end_trajPy + 
                                                end_trajPz*end_trajPz );
        
                end_trajX[p] = endPoint->GetPosition().x();
                end_trajY[p] = endPoint->GetPosition().y();
                end_trajZ[p] = endPoint->GetPosition().z();
                Gaudi::XYZPoint point(end_trajX[p],end_trajY[p],end_trajZ[p]);
                inside_od[p] = int(m_OuterDetector->isInside(point));
            } 
        }
    } //! end loop over prongs
    
    //-- fill
    if( prongs[0]->filtertaglist()->filterTagExists("PrimaryProton") ) {
        neutrino->setContainerIntData("ntrajProtonProng",ntraj);
        neutrino->setContainerIntData("trajProtonProngPrimary",nprim);
        neutrino->setContainerIntData("trajProtonProngPDG",pdg);
        neutrino->setContainerDoubleData("trajProtonProngPx",traj_Px);
        neutrino->setContainerDoubleData("trajProtonProngPy",traj_Py);
        neutrino->setContainerDoubleData("trajProtonProngPz",traj_Pz);
        neutrino->setContainerDoubleData("trajProtonProngEnergy",traj_E);
        neutrino->setContainerDoubleData("trajProtonProngMomentum",traj_P);
        neutrino->setContainerDoubleData("trajProtonProngPSelf",traj_PSelf);
        neutrino->setContainerDoubleData("trajProtonTheta",traj_theta);
        neutrino->setContainerDoubleData("trajProtonPhi",traj_phi);
        neutrino->setContainerIntData("isProtonInsideOD",inside_od);
        neutrino->setContainerDoubleData("endProtonTrajMomentum",end_trajMomentum);
        neutrino->setContainerDoubleData("endProtonTrajXPosition",end_trajX);
        neutrino->setContainerDoubleData("endProtonTrajYPosition",end_trajY);
        neutrino->setContainerDoubleData("endProtonTrajZPosition",end_trajZ);
    } else if( prongs[0]->filtertaglist()->filterTagExists("PrimaryMuon") ) {
        neutrino->setIntData("ntrajMuonProng",ntraj[0]);
        neutrino->setIntData("trajMuonProngPrimary",nprim[0]);
        neutrino->setIntData("trajMuonProngPDG",pdg[0]);
        neutrino->setDoubleData("trajMuonProngPx",traj_Px[0]);
        neutrino->setDoubleData("trajMuonProngPy",traj_Py[0]);
        neutrino->setDoubleData("trajMuonProngPz",traj_Pz[0]);
        neutrino->setDoubleData("trajMuonProngEnergy",traj_E[0]);
        neutrino->setDoubleData("trajMuonProngMomentum",traj_P[0]);
        neutrino->setDoubleData("trajMuonProngPSelf",traj_PSelf[0]);
        neutrino->setDoubleData("trajMuonTheta",traj_theta[0]);
        neutrino->setDoubleData("trajMuonPhi",traj_phi[0]);
        neutrino->setIntData("isMuonInsideOD",inside_od[0]);
        neutrino->setDoubleData("endMuonTrajMomentum",end_trajMomentum[0]);
        neutrino->setDoubleData("endMuonTrajXPosition",end_trajX[0]);
        neutrino->setDoubleData("endMuonTrajYPosition",end_trajY[0]);
        neutrino->setDoubleData("endMuonTrajZPosition",end_trajZ[0]);
    }
    
    
    debug() << "Exit CCProtonPi0::setTrackProngTruth()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    
    return;
}




//==============================================================================
// Set Muon particle data
//==============================================================================
StatusCode CCProtonPi0::setMuonParticleData(   Minerva::NeutrinoInt* nuInt, SmartRef<Minerva::Prong>& muonProng) const 
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter CCProtonPi0::setMuonParticleData()" << endmsg;
    
    Gaudi::LorentzVector muon_4p;
    
    SmartRefVector<Track>::iterator iterTrk;
    SmartRef<Minerva::Particle> muonPart = (Minerva::Particle*)NULL;
    int muon_minervaTrack_types = 0;
    int muon_N_minosTracks = 0;
    if (!muonProng) {
        fatal()<<"NO Muon!!"<<endmsg;
        return StatusCode::FAILURE;
    }
    
    SmartRefVector<Minerva::Particle> particle_hypotheses = muonProng->particles();
    SmartRefVector<Minerva::Particle>::iterator iterPart;
    for( iterPart = particle_hypotheses.begin(); iterPart != particle_hypotheses.end(); ++iterPart ) {  
        bool is_muon_particle = false;
        (*iterPart)->filtertaglist()->checkFilterTag( "PrimaryMuon", is_muon_particle );
        if (is_muon_particle) {
            muonPart = *iterPart;
            break;
        }
    }
    
    SmartRefVector<Track> muonTracks = muonProng->minervaTracks();
    int nlong = 0, nshort = 0;
    for (iterTrk = muonTracks.begin(); iterTrk != muonTracks.end(); ++iterTrk) {
        if ((*iterTrk)->type() == Track::Long) nlong++;
        if ((*iterTrk)->type() == Track::Short) nshort++;
    }
    
    if (nlong > 0 && nshort == 0) muon_minervaTrack_types = 1;
    else if (nlong==0 && nshort > 0) muon_minervaTrack_types = 2;
    else if (nlong>0 && nshort>0) muon_minervaTrack_types = 3;
    
    muon_N_minosTracks = muonProng->minosTracks().size();   
    
    //! muon energy in road upstream info
    int muon_roadUpstreamPlanes = -1;
    int charge = -99;
    double muon_roadUpstreamEnergy = AnaToolUtils->energyInTheRoadUpstream(muonProng, muon_roadUpstreamPlanes);  
    
    //! Fill muon particle hypothesis info
    //! @todo - cleanup redundant minos track information
    double muon_E = -9.0, muon_p = 0.0, muon_px = 0.0, muon_py = 0.0, muon_pz = 0.0;
    double muon_theta = -9.0, muon_theta_biasUp = -9.0, muon_theta_biasDown = -9.0;
    double muon_muScore = -1.0, muon_qp = -9999.9, muon_qpqpe = -9999.9, muon_E_shift;
    int muon_minosTrackQuality = 0;
    if (muonPart) {
        muon_4p    = muonPart->momentumVec();
        muon_px    = muon_4p.px();
        muon_py    = muon_4p.py();
        muon_pz    = muon_4p.pz();    
        muon_E     = muon_4p.E();
        muon_p     = muon_4p.P();
        muon_theta = m_coordSysTool->thetaWRTBeam(muon_4p);
        muon_theta_biasUp = m_coordSysTool->thetaWRTBeam(muon_4p,m_beamAngleBias) - muon_theta;
        muon_theta_biasDown = m_coordSysTool->thetaWRTBeam(muon_4p, -1.0*m_beamAngleBias) - muon_theta;
        
        debug()<<"Muon:  "<<muon_px<<" "<<muon_py<<" "<<muon_pz<<" "<<muon_E<<endmsg;
        
        muon_muScore = muonPart->score();
        
        muon_qpqpe = MuonUtils->minosQPQPE(muonProng);
        if (muonProng->MinosTrack()) {
            SmartRef<MinosRecoTrack> minosTrack = muonProng->minosTracks()[0];
            muon_qp = minosTrack->qp();
            muon_minosTrackQuality = minosTrack->quality();      
        }
        else if ( !muonProng->MinosStub() ) {
            warning()<<"No MINOS track or stub!  How is this a muon?"<<endmsg;
        }
        
        //the shifts are so small compared to MINOS-matched momentum that we can approximate the momentum shift as the energy shift (with ~0.2% at 1.5 GeV/c)
        muon_E_shift = MuonUtils->calculateMomentumCorrection(muonProng);
        
        //Get Muon Charge
        MuonUtils->muonCharge(muonProng,charge); 
    }
    
    debug()<<"Filling Muon Ntuple Variables"<<endmsg;
    
    nuInt->setLeptonEnergy( muon_4p );
    
    //! ntuple Muon Variables
    nuInt->setIntData("muon_minervaTrack_types", muon_minervaTrack_types);
    nuInt->setIntData("muon_N_minosTracks", muon_N_minosTracks);
    nuInt->setIntData("muon_minosTrackQuality", muon_minosTrackQuality);
    nuInt->setIntData("muon_roadUpstreamPlanes", muon_roadUpstreamPlanes);
    nuInt->setIntData("muon_charge",charge);
    nuInt->setDoubleData("muon_roadUpstreamEnergy", muon_roadUpstreamEnergy);
    nuInt->setDoubleData("muon_E",muon_E);
    nuInt->setDoubleData("muon_p",muon_p);
    nuInt->setDoubleData("muon_px",muon_px);
    nuInt->setDoubleData("muon_py",muon_py);
    nuInt->setDoubleData("muon_pz",muon_pz);
    nuInt->setDoubleData("muon_theta",muon_theta);
    nuInt->setDoubleData("muon_theta_biasUp",muon_theta_biasUp);
    nuInt->setDoubleData("muon_theta_biasDown",muon_theta_biasDown);
    nuInt->setDoubleData("muon_muScore", muon_muScore);
    nuInt->setDoubleData("muon_qp",muon_qp );
    nuInt->setDoubleData("muon_qpqpe",muon_qpqpe);
    nuInt->setDoubleData("muon_E_shift",muon_E_shift);

    
    fillMinosMuonBranches(nuInt, muonProng);

    debug() << "Exit CCProtonPi0::setMuonParticleData()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    return StatusCode::SUCCESS;
}

//==============================================================================
// Find the plane nearest to a point
//==============================================================================
StatusCode CCProtonPi0::getNearestPlane(double z, int & module_return, int & plane_return) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter CCProtonPi0::getNearestPlane()" << endmsg;
    
    // Got From Brandon's CCNuPionInc - 2014_04_14
    // testing new MinervaDet routine
    // works as advertised, but a couple problems:
    // 1) Does not go to more downstream plane (upstream if backwards track) if point is in passive material
    
    DePlane const * pPlane = m_InnerDetector->getClosestDePlane(z);
    if (!pPlane) {
        error()<<"getClosestDePlane failed on z "<<z<<endmsg;
        return StatusCode::FAILURE;
    }
    
    PlaneID planeid = pPlane->getPlaneID();
    module_return = planeid.module();
    plane_return = planeid.plane();

    debug() << "Exit CCProtonPi0::getNearestPlane()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    return StatusCode::SUCCESS;
  
}

//==============================================================================
// Created particles for negative bit Minos prongs
//==============================================================================
bool CCProtonPi0::createTrackedParticles( Minerva::ProngVect& prongs ) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() <<"Enter CCProtonPi0::createTrackedParticles()" << endmsg;
    bool makeParticles = false;
    
    //-- check if the prongs are odMatch
    Minerva::EventMgr* mgr = getEventMgr(prongs[0]);
    Minerva::ODClusterVect odClusters = mgr->select<Minerva::ODCluster>("Unused","!LowActivity&!XTalkCandidate");
    m_odMatchTool->classifyProngs(prongs,odClusters);
    
    //-- loop over prongs
    for(unsigned int p = 0; p < prongs.size(); p++) {
        debug() << "The prong of bit-field = " << prongs[p]->typeBitsToString() << endmsg;
    
        //-- make sure the prong is not a bad object
        bool pass = true; std::string tag = "BadObject";
        if( prongs[p]->filtertaglist()->filterTagExists(tag) ) {
            if( prongs[p]->filtertaglist()->checkFilterTag(tag,pass) ) {
                error() << "This prong = " << prongs[p]->typeBitsToString() << " has been flag as a \"BadObject\", skipping!" << endmsg;
                continue;
            }
        }
        
        //-- skipped prongs with particles
        debug() << "The prong has n particles = " << prongs[p]->particles().size() << endmsg;
        if( !prongs[p]->particles().empty() ) continue;
    
        //! select the particle hypotheses candidates and particle tools
        std::vector< Minerva::Particle::ID > hypotheses;
        IParticleMakerTool::NameAliasListType toolsToUse;

        if( !prongs[p]->OdMatch() ) {
            hypotheses.push_back( Minerva::Particle::Proton );
            hypotheses.push_back( Minerva::Particle::Pion );
            toolsToUse.push_back( std::make_pair("dEdXTool","dEdXTool") );
        }

        //-- make particles
        Minerva::Prong* prong = prongs[p];
        bool isParticleCreated = m_particleMaker->makeParticles(prong,hypotheses,toolsToUse);
    
        if( isParticleCreated ) {
            debug() << "The prong of bit-field = " << prongs[p]->typeBitsToString() 
                    << " has " << prong->particles().size() << " number of particles attached." << endmsg; 
            makeParticles = true;
        } else{
            debug() << "Did not make particles for the prong type = " << prongs[p]->typeBitsToString() << endmsg;
        }
    
        hypotheses.clear();
    }
    
    debug() << "Exit CCProtonPi0::createTrackedParticles()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    return makeParticles;
}

//==============================================================================
// Return the momentum analyzable contained ( proton candidate ) prong/particle
//==============================================================================
bool CCProtonPi0::getProtonProng(    Minerva::ProngVect& hadronProngs, 
                                        Minerva::ProngVect& protonProngs,
                                        Minerva::ParticleVect& protonParticles ) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter CCProtonPi0::getProtonProng()" << endmsg;
  
    debug() << "N(hadronProngs) =  " << hadronProngs.size() << endmsg;
    
    // Initialize
    bool isProton = false;
    int nProtons = 0;
    int nGoodProtons = 0; 
    
    // Get All Proton Candidates
    for(unsigned int p = 0; p < hadronProngs.size(); p++) {
        debug() <<"Checking prong "<<p<<endmsg;
        
        // Sanity Check
        if( hadronProngs[p]->filtertaglist()->filterTagExists("PrimaryMuon") ){
            debug() << "Muon Prong, skipping this prong! "<< endmsg;
            continue;
        }
    
        // Temp Storage for ProngVect, Single Prong and Single Particle
        Minerva::ProngVect tmpProngs;
        SmartRef<Minerva::Prong> prong       = (Minerva::Prong*)NULL;
        SmartRef<Minerva::Particle> particle = (Minerva::Particle*)NULL;
        
        // Push current Hadron Prong to temp ProngVect
        tmpProngs.push_back( hadronProngs[p] );
        
        // Find Proton using m_protonUtils
        isProton = m_protonUtils->findProtonProng(tmpProngs,prong,particle); 
        
        if( isProton ) {
            nProtons++;
            debug() <<"Found a proton prong!"<< endmsg;
            debug() <<"Proton Particle Score: " << particle->score() << endmsg;
            prong->filtertaglist()->setOrAddFilterTag("PrimaryProton",true);
            particle->filtertaglist()->setOrAddFilterTag("PrimaryProton",true);
            prong->updateBestParticle(particle);
            protonProngs.push_back( prong );
            protonParticles.push_back( particle );
        }
    }

    debug() << "Found "<<nProtons<<" Proton Candidates!"<<endmsg;
    debug() << "Applying minProtonScore Cut with = "<<m_minProtonScore<<endmsg;
    
    // Apply Score Cut
    for(unsigned int p = 0; p < protonParticles.size(); p++) {
        if (protonParticles[p]->score() > m_minProtonScore){
            nGoodProtons++;
        }else{
            protonProngs.erase(protonProngs.begin()+p);
            protonParticles.erase(protonParticles.begin()+p);
        }
    }
   
    
    debug() << "Found "<<nGoodProtons<<" Good Protons!"<<endmsg;
    debug() << "N(protonProngs) =  " << protonProngs.size() << endmsg;
    debug() << "N(protonParticles) =  " << protonParticles.size() << endmsg;
   
    
    debug() << "Exit CCProtonPi0::getProtonProng()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    
    return isProton;
}

//==============================================================================
// Set proton particle data
//==============================================================================
void CCProtonPi0::setProtonParticleData(     Minerva::NeutrinoInt* nuInt, 
                                                Minerva::ProngVect& protonProngs,
                                                Minerva::ParticleVect& protonParticles, 
                                                double vertexZ ) const 
{
    // Got from Tammy's NukeCCQETwoTrack - 2014_04_15
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter CCProtonPi0::setProtonParticleData()" << endmsg;
    
    std::vector<double> p_calCorrection(10,-1);
    std::vector<double> p_visEnergyCorrection(10,-1);
    std::vector<double> p_dedx(10,-1);
    
    std::vector<double> proton_theta(10,-1);
    std::vector<double> proton_thetaX(10,-1);
    std::vector<double> proton_thetaY(10,-1);
    std::vector<double> proton_phi(10,-1);
    
    std::vector<double> E(10,-1);
    std::vector<double> px(10,-1);
    std::vector<double> py(10,-1);
    std::vector<double> pz(10,-1);
    std::vector<double> p(10,-1);
    
    std::vector<double> ekin(10,-1);
    std::vector<double> enu(10,-1);
    std::vector<double> Q2(10,-1);
    
    std::vector<double> proton_end_x(10,-1);
    std::vector<double> proton_end_y(10,-1);
    std::vector<double> proton_end_z(10,-1);
    
    std::vector<double> proton_start_x(10,-1);
    std::vector<double> proton_start_y(10,-1);
    std::vector<double> proton_start_z(10,-1);
    
    std::vector<int> kinked(10,-1);
    std::vector<int> odMatch(10,-1);
    
    std::vector<double> score(10,-1);
    std::vector<double> score1(10,-1);
    std::vector<double> score2(10,-1);
    std::vector<double> chi2(10,-1);
    
//     std::map< std::string, std::vector<double> > dedxMapMomentum;
//     std::map< std::string, std::vector<double> > dedxMapScore;
    
//     for(unsigned int i = 0; i < m_dedx_uncertainties.size(); i++) {
//         std::string name = m_dedx_uncertainties[i];
//         std::vector<double> tmp(10,-1);
//         dedxMapMomentum.insert( std::pair< std::string,std::vector<double> >(name,tmp) );
//         dedxMapScore.insert( std::pair< std::string,std::vector<double> >(name,tmp) );
//     }
    
    for(unsigned int i = 0; i < protonProngs.size(); i++) {
        double theta = protonProngs[i]->minervaTracks().front()->theta();
        double phi   = protonProngs[i]->minervaTracks().front()->phi();
    
        proton_theta[i]  = m_coordSysTool->thetaWRTBeam(protonParticles[i]->momentumVec());
        proton_thetaX[i] = m_coordSysTool->thetaXWRTBeam(protonParticles[i]->momentumVec());
        proton_thetaY[i] = m_coordSysTool->thetaYWRTBeam(protonParticles[i]->momentumVec());
        proton_phi[i]    = m_coordSysTool->phiWRTBeam(protonParticles[i]->momentumVec());
//         double p_original = sqrt( pow(protonParticles[i]->momentumVec().E(),2) - MinervaUnits::M_p*MinervaUnits::M_p );
    
        Gaudi::LorentzVector dEdXprotonfourVec;
        StatusCode sc = m_energyCorrectionTool->getCorrectedEnergy(protonProngs[i],protonParticles[i],vertexZ,dEdXprotonfourVec);
        if( !sc ) dEdXprotonfourVec = protonParticles[i]->momentumVec();
        else protonParticles[i]->setMomentumVec(dEdXprotonfourVec);
    
        p_dedx[i] = sqrt( dEdXprotonfourVec.E()*dEdXprotonfourVec.E() - MinervaUnits::M_proton*MinervaUnits::M_proton );
        correctProtonProngEnergy(protonProngs[i],p_calCorrection[i],p_visEnergyCorrection[i]);
    
        E[i]  = sqrt( p_calCorrection[i]*p_calCorrection[i] + MinervaUnits::M_proton*MinervaUnits::M_proton );
        px[i] = p_calCorrection[i]*sin(theta)*cos(phi);
        py[i] = p_calCorrection[i]*sin(theta)*sin(phi);
        pz[i] = p_calCorrection[i]*cos(theta);
        p[i]  = sqrt( px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i] );
    
        if( protonProngs[i]->minervaTracks().front()->direction() == Minerva::Track::Backward && pz[i] > 0. ) pz[i] *= -1.0;
        Gaudi::LorentzVector protonfourVec(px[i],py[i],pz[i],E[i]);
    
        ekin[i] = protonfourVec.E() - MinervaUnits::M_p;
    
        proton_end_x[i]   = (protonProngs[i]->minervaTracks().back())->lastState().x();
        proton_end_y[i]   = (protonProngs[i]->minervaTracks().back())->lastState().y();
        proton_end_z[i]   = (protonProngs[i]->minervaTracks().back())->lastState().z();
    
        proton_start_x[i] = (protonProngs[i]->minervaTracks().back())->firstState().x();
        proton_start_y[i] = (protonProngs[i]->minervaTracks().back())->firstState().y();
        proton_start_z[i] = (protonProngs[i]->minervaTracks().back())->firstState().z();
    
        kinked[i]  = (int)protonProngs[i]->Kinked();
        odMatch[i] = (int)protonProngs[i]->OdMatch();
        score[i]   = protonParticles[i]->score();
        score1[i]  = protonParticles[i]->getDoubleData("score1");
        score2[i]  = protonParticles[i]->getDoubleData("score2");
        chi2[i]    = protonParticles[i]->getDoubleData("chi2_ndf");
    
//         for(unsigned int j = 0; j < m_dedx_uncertainties.size(); j++) {
//         std::string name = m_dedx_uncertainties[j];
//     
//         double score_varied  = protonParticles[i]->hasDoubleData("score1_" + name) ? protonParticles[i]->getDoubleData("score1_" + name) : 0.0;
//         double energy_varied = protonParticles[i]->hasDoubleData("FitE_"   + name) ? protonParticles[i]->getDoubleData("FitE_"   + name) : 0.0;
//         double p_varied      = sqrt( energy_varied*energy_varied - MinervaUnits::M_proton*MinervaUnits::M_proton );
//     
//         double score_shift   = protonParticles[i]->getDoubleData("score1") - score_varied;
//         double p_shift       = p_original - p_varied;
//         
//         (dedxMapMomentum.find(name))->second.at(i) = p_shift;
//         (dedxMapScore.find(name))->second.at(i)    = score_shift;
//         }
    }
    
    nuInt->setContainerDoubleData("proton_p_calCorrection",p_calCorrection);
    nuInt->setContainerDoubleData("proton_p_visEnergy",p_visEnergyCorrection);
    nuInt->setContainerDoubleData("proton_p_dEdXTool",p_dedx);
    
    nuInt->setContainerDoubleData("proton_endPointX",proton_end_x);
    nuInt->setContainerDoubleData("proton_endPointY",proton_end_y);
    nuInt->setContainerDoubleData("proton_endPointZ",proton_end_z);
    
    nuInt->setContainerDoubleData("proton_startPointX",proton_start_x);
    nuInt->setContainerDoubleData("proton_startPointY",proton_start_y);
    nuInt->setContainerDoubleData("proton_startPointZ",proton_start_z);
    
    nuInt->setContainerDoubleData("proton_p",p);
    nuInt->setContainerDoubleData("proton_px",px);
    nuInt->setContainerDoubleData("proton_py",py);
    nuInt->setContainerDoubleData("proton_pz",pz);
    nuInt->setContainerDoubleData("proton_E",E);
    
    nuInt->setContainerIntData("proton_kinked",kinked);
    nuInt->setContainerIntData("proton_odMatch",odMatch);
    
    nuInt->setContainerDoubleData("proton_score",score);
    nuInt->setContainerDoubleData("proton_score1",score1);
    nuInt->setContainerDoubleData("proton_score2",score2);
    nuInt->setContainerDoubleData("proton_chi2_ndf",chi2);
    
    nuInt->setContainerDoubleData("proton_theta",proton_theta);
    nuInt->setContainerDoubleData("proton_thetaX",proton_thetaX);
    nuInt->setContainerDoubleData("proton_thetaY",proton_thetaY);
    nuInt->setContainerDoubleData("proton_phi",proton_phi);
    
    nuInt->setContainerDoubleData("proton_ekin",ekin);
    
//     for(unsigned int j = 0; j < m_dedx_uncertainties.size(); j++) {
//         std::string name = m_dedx_uncertainties[j];
//     
//         std::vector<double> p_shift_vect     = dedxMapMomentum.find(name)->second;
//         std::vector<double> score_shift_vect = dedxMapScore.find(name)->second;
//     
//         nuInt->setContainerDoubleData("proton_momentum_shift_" + name, p_shift_vect);
//         nuInt->setContainerDoubleData("proton_score1_shift_" + name, score_shift_vect);
//     }
    
    debug() << "Exit CCProtonPi0::setProtonParticleData()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    
    return;
}

//==============================================================================
// Correct Proton Energy
//==============================================================================
void CCProtonPi0::correctProtonProngEnergy(  SmartRef<Minerva::Prong>& protonProng, 
                                                double& p_calCorrection, 
                                                double& p_visEnergyCorrection ) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter CCProtonPi0::correctProtonProngEnergy()" << endmsg;
    
    //-- initialize
    double calEnergy = 0.0;
    double visEnergy = 0.0;
        
    //-- get the proton prong subprongs
    Minerva::ProngVect subProngs = protonProng->subProngs();
    if( !subProngs.empty() ) {
        debug() << "The number of subProngs = " << subProngs.size() << endmsg;
    
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
        debug() << "The calorimetric energy = " << calEnergy << " MeV" << endmsg;
    
        //-- get the total visible energy of the sub prongs
        Minerva::IDClusterVect clusters = subProngs[p]->getAllIDClusters();
        for(unsigned int clus = 0; clus < clusters.size(); clus++) visEnergy += clusters[clus]->energy();
        debug() << "The visible energy = " << visEnergy << " MeV" << endmsg;
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
        debug() << "particle idcode = " << particles[part]->idcode() << ", and four momentum = " << fourMomentum << endmsg;
    
        double E_calCorrection = fourMomentum.E() + calEnergy;
        p_calCorrection = sqrt( pow(E_calCorrection,2) - pow(particles[part]->mass(),2) );
        debug() << "update energy using calorimetric correction = " << E_calCorrection << endmsg;
    
        double E_visEnergyCorrection = fourMomentum.E() + visEnergy;
        p_visEnergyCorrection = sqrt( pow(E_visEnergyCorrection,2) - pow(particles[part]->mass(),2) );     
        debug() << "update energy using visible energy correction = " << E_visEnergyCorrection << endmsg;
    }     
            
    debug() << "Exit CCProtonPi0::correctProtonProngEnergy()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    return;
}

//==============================================================================
// setPi0ParticleData
//==============================================================================

StatusCode CCProtonPi0::setPi0ParticleData( Minerva::PhysicsEvent *event, 
                                            Minerva::IDBlob* idblob1, 
                                            Minerva::IDBlob* idblob2,
                                            const SmartRef<Minerva::Vertex>& vertex ) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "CCProtonPi0::setPi0ParticleData()" << endmsg;

    const Gaudi::XYZPoint& pos = vertex->position();
            
    bool goodPosition1  = m_idHoughBlob->GetStartPosition(idblob1, pos, true );
    bool goodDirection1 = m_idHoughBlob->GetDirection(idblob1, pos );
    bool isGoodBlob1 = false;
    if (goodPosition1 && goodDirection1) isGoodBlob1 = true;
    
    bool goodPosition2  = m_idHoughBlob->GetStartPosition(idblob2, pos, true );
    bool goodDirection2 = m_idHoughBlob->GetDirection(idblob2, pos );
    bool isGoodBlob2 = false;
    if (goodPosition2 && goodDirection2) isGoodBlob2 = true;
    // Process only if we have 2 Bood Blobs
    if (!isGoodBlob1 || !isGoodBlob2){
        debug() << "Don't have 2 Good Blobs! -- CANNOT Set Pi0 Particle Data!" <<endmsg; 
        return StatusCode::SUCCESS;
    }

    if (attenuation_correction_) {
        ApplyAttenuationCorrection(idblob1);
        ApplyAttenuationCorrection(idblob2);
    }
    
    double g1energy   = 0.0;
    double g1trkrevis = 0.0;
    double g1ecalevis = 0.0;
    double g1hcalevis = 0.0;
    double g1scalevis = 0.0;
    
    double g2energy   = 0.0;
    double g2trkrevis = 0.0;
    double g2ecalevis = 0.0;
    double g2hcalevis = 0.0;
    double g2scalevis = 0.0;

    
    m_idHoughBlob->getBlobEnergyTime( idblob1, g1energy, g1trkrevis, g1ecalevis, g1hcalevis, g1scalevis );
    m_idHoughBlob->getBlobEnergyTime( idblob2, g2energy, g2trkrevis, g2ecalevis, g2hcalevis, g2scalevis );

    /* Make sure gamma1 is the more energetic one */
    if (g2energy > g1energy) {
        std::swap(goodPosition1, goodPosition2);  /* Swap variables already assigned */
        std::swap(goodDirection1,goodDirection2);
        std::swap(isGoodBlob1,   isGoodBlob2);
        std::swap(g1energy,g2energy);
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
    
    debug() << "Direction: " << direction1 << " " << direction2
              << cos_oangle << endmsg;

    event->setDoubleData( "RE_scalar", cos_oangle );
    event->setDoubleData( "RE_photon_energy_1", g1energy );
    event->setDoubleData( "RE_photon_energy_2", g2energy );
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
    g1mom *= g1energy;
    g2mom *= g2energy;

    TVector3 pimom = g1mom + g2mom;
    
    const double oangle = g1mom.Angle(g2mom)*TMath::RadToDeg();
    event->setDoubleData("oangle", oangle);
    
    const double mgg = std::sqrt(2*g1energy*g2energy*(1-cos_oangle));
    event->setDoubleData("mgg", mgg);

    // CCProtonPi0 Style Branches
    event->setDoubleData("pi0_px",pimom.x());
    event->setDoubleData("pi0_py",pimom.y());
    event->setDoubleData("pi0_pz",pimom.z());
    event->setDoubleData("pi0_E",g1energy+g2energy);
    
    event->setDoubleData("gamma1_px",g1mom.x());
    event->setDoubleData("gamma1_py",g1mom.y());
    event->setDoubleData("gamma1_pz",g1mom.z());
    event->setDoubleData("gamma1_E",g1energy);
    
    event->setDoubleData("gamma2_px",g2mom.x());
    event->setDoubleData("gamma2_py",g2mom.y());
    event->setDoubleData("gamma2_pz",g2mom.z());
    event->setDoubleData("gamma2_E",g2energy);

    // CCPi0 Default Branches
    event->setDoubleData("pienergy", g1energy+g2energy);
    event->setDoubleData("pitheta",  pimom.Theta()*TMath::RadToDeg());
    event->setDoubleData("piphi",    pimom.Phi()*TMath::RadToDeg());
    event->setDoubleData("pithetax", std::atan2(pimom.X(),pimom.Z())*TMath::RadToDeg());
    event->setDoubleData("pithetay", std::atan2(pimom.Y(),pimom.Z())*TMath::RadToDeg());

    event->setDoubleData("g1e",    g1energy);
    event->setDoubleData("g1theta",g1mom.Theta()*TMath::RadToDeg());
    event->setDoubleData("g1phi",  g1mom.Phi()*TMath::RadToDeg());
    event->setDoubleData("g2e",    g2energy);
    event->setDoubleData("g2theta",g2mom.Theta()*TMath::RadToDeg());
    event->setDoubleData("g2phi",  g2mom.Phi()*TMath::RadToDeg());

    event->setDoubleData("g1recotrkrevis", g1trkrevis);
    event->setDoubleData("g1recoecalevis", g1ecalevis);
    event->setDoubleData("g1recohcalevis", g1hcalevis);
    event->setDoubleData("g1recoscalevis", g1scalevis);

    event->setDoubleData("g2recotrkrevis", g2trkrevis);
    event->setDoubleData("g2recoecalevis", g2ecalevis);
    event->setDoubleData("g2recohcalevis", g2hcalevis);
    event->setDoubleData("g2recoscalevis", g2scalevis);

    
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

    CalculatedEdx(idblob1,event,1,vertex);
    CalculatedEdx(idblob2,event,2,vertex);
    
//     CheckBlobDigitAssignment(event,idblob1,1);
//     CheckBlobDigitAssignment(event,idblob2,2);

//         // Convenient variables for analysis
//         // 1) Muon and pi0 momentum in the detector coordinates
//     const double pmu_x = m_MuonParticle->momentumVec().Px();
//     const double pmu_y = m_MuonParticle->momentumVec().Py();
//     const double pmu_z = m_MuonParticle->momentumVec().Pz();
//     const double Emu    = m_MuonParticle->momentumVec().E();
//     
//     const double ppi_x = pimom.x();
//     const double ppi_y = pimom.y();
//     const double ppi_z = pimom.z();
//     const double Epi    = g1energy + g2energy;
//     
//         // 2) transform into beam coordinates
//     const double pmu_y2 =  costheta_b*pmu_y + sintheta_b*pmu_z;
//     const double pmu_z2 = -sintheta_b*pmu_y + costheta_b*pmu_z;
//     
//     const double ppi_y2 =  costheta_b*ppi_y + sintheta_b*ppi_z;
//     const double ppi_z2 = -sintheta_b*ppi_y + costheta_b*ppi_z;
// 
//         // 3) Calculate the neutrino energy
//     const double pTx = pmu_x + ppi_x;
//     const double pTy = pmu_y2 + ppi_y2;
//     
//     const double pmuL = Emu - pmu_z2;
//     const double ppiL = Epi - ppi_z2;
//     const double pL   = pmuL + ppiL;
//     const double Tn   = 0.5*(pL*pL + pTx*pTx + pTy*pTy)/(mn - pL);
//     const double Erec = Emu + Epi + Tn;
// 
//     const double Q2 = 2*Erec*(Emu-pmu_z2) - mmuon*mmuon;
//     const double W2 = mp*mp + 2*mp*(Erec-Emu) - Q2;
//     const double W  = std::sqrt(std::max(0.,W2));
// 
//         // Re-calculate these kinematic variables using 'better' pi0 energy
//     const double m0 = 134.9764;
//     const double Erec_es = Emu + m0/mgg*Epi + Tn;
//     const double Q2_es = 2*Erec_es*(Emu-pmu_z2) - mmuon*mmuon;
//     const double W2_es = mp*mp + 2*mp*(Erec_es-Emu) - Q2_es;
//     const double W_es  = std::sqrt(std::max(0.,W2_es));

    
//     event->setDoubleData("Tn", Tn);
//     event->setDoubleData("Erec", Erec);
//     event->setDoubleData("Q2", Q2);
//     event->setDoubleData("W2", W2);
//     event->setDoubleData("W", W);
//     
//     event->setDoubleData("Erec_es", Erec_es);
//     event->setDoubleData("Q2_es", Q2_es);
//     event->setDoubleData("W2_es", W2_es);
//     event->setDoubleData("W_es", W_es);


    
//     std::vector<double> mumom;
//     mumom.push_back(pmu_x);
//     mumom.push_back(pmu_y);
//     mumom.push_back(pmu_z);
//     mumom.push_back(m_MuonParticle->momentumVec().R());
//     mumom.push_back(Emu);
//     
//     event->setContainerDoubleData("mumom", mumom);

//     TVector3 pimom2(pimom.X(),ppi_y2,ppi_z2);
//     event->setDoubleData("pithetab",  pimom2.Theta()*TMath::RadToDeg());
//     event->setDoubleData("piphib",    pimom2.Phi()*TMath::RadToDeg());
//     event->setDoubleData("pithetaxb", std::atan2(pimom2.X(),pimom2.Z())*TMath::RadToDeg());
//     event->setDoubleData("pithetayb", std::atan2(pimom2.Y(),pimom2.Z())*TMath::RadToDeg());

    
    debug() << "Exit CCProtonPi0::setPi0ParticleData()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;

    return StatusCode::SUCCESS;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
//
//  CCPi0 Functions -- Blob Tools
//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


//==============================================================================
//  VtxBlob
//==============================================================================
StatusCode CCProtonPi0::VtxBlob(Minerva::PhysicsEvent *event, const SmartRef<Minerva::Vertex>& vertex ) const
{
    
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter CCProtonPi0::VtxBlob()" << endmsg;
    
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
            double dummy_trkrevis = 0.0;
            double dummy_ecalevis = 0.0;
            double dummy_hcalevis = 0.0;
            double dummy_scalevis = 0.0;
            m_idHoughBlob->getBlobEnergyTime(vtxFilaBlob,vertex_energy_filament, 
                                            dummy_trkrevis,dummy_ecalevis, dummy_hcalevis, dummy_scalevis);
    
            debug() << " Adding Filament vertex " << vtxFilaBlob->nclusters()
                    << " clusters; energy = "  << vertex_energy_filament << endmsg;
            addObject( event, vtxFilaBlob );
            m_hitTagger->applyColorTag( vtxFilaBlob, 0x9900FF );
        }
    } // if (m_filamentVertex)
    
    
    if ( m_sphereVertex ){
    
        SmartRefVector<Minerva::IDCluster> unusedClusters
            = event->select<Minerva::IDCluster>("Unused","!LowActivity&!XTalkCandidate");
    
        std::vector<double> radii;
        radii.push_back(50.0);
        radii.push_back(100.0);
        radii.push_back(150.0);
        radii.push_back(200.0);
        radii.push_back(250.0);
        radii.push_back(300.0);
        radii.push_back(500.0);
    
        // Note: return the vertex energies within these 'shells' in the parameter 'radii'
        SmartRefVector<Minerva::IDCluster> sphereClusters = FilterInSphereClusters(vertex,unusedClusters,m_maxSearchD,radii);
    
        event->setContainerDoubleData("Vertex_energy_radii", radii);
        
        if (!sphereClusters.empty()) {
    
            Minerva::IDBlob* vtxSphereBlob = new Minerva::IDBlob();
            m_blobUtils->insertIDClusters( sphereClusters, vtxSphereBlob, Minerva::IDBlob::VertexBlobPatRec );
            double dummy_trkrevis = 0.0;
            double dummy_ecalevis = 0.0;
            double dummy_hcalevis = 0.0;
            double dummy_scalevis = 0.0;
            m_idHoughBlob->getBlobEnergyTime(vtxSphereBlob,vertex_energy_sphere,
                                            dummy_trkrevis, dummy_ecalevis, dummy_hcalevis, dummy_scalevis);
    
            debug() << " Adding Sphere vertex " << vtxSphereBlob->nclusters()
                    << " clusters; energy = "  <<  vertex_energy_sphere << endmsg;
            addObject( event, vtxSphereBlob );
            m_hitTagger->applyColorTag( vtxSphereBlob, 0x9900FF );
        }
    } // if ( m_sphereVertex )
    
    vertex_energy = vertex_energy_sphere + vertex_energy_filament;
    event->setDoubleData( "Vertex_blob_energy", vertex_energy );
    event->setDoubleData( "Filament_Vertex_energy", vertex_energy_filament );
    event->setDoubleData( "Sphere_Vertex_energy", vertex_energy_sphere );
  
    debug() << "Exit CCProtonPi0::VtxBlob()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;

  return StatusCode::SUCCESS;

}

//==============================================================================
//  FilterInSphereClusters()
//==============================================================================
SmartRefVector<Minerva::IDCluster> CCProtonPi0::FilterInSphereClusters(const SmartRef<Minerva::Vertex>& vertex,
                                     const SmartRefVector<Minerva::IDCluster>& clusters,
                                     const double sphereRadius,
                                     std::vector<double>& radii) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter CCProtonPi0::FilterInSphereClusters()" << endmsg;
    
    const Gaudi::XYZPoint& pos = vertex->position();
    const double x0 = pos.X();
    const double y0 = pos.Y();
    const double z0 = pos.Z();
    double u0 = m_mathTool->calcUfromXY(x0,y0);
    double v0 = m_mathTool->calcVfromXY(x0,y0);

    std::vector<double> vertexEnergyVector(radii.size(),0.0);
    SmartRefVector<Minerva::IDCluster> sphereClusters;
    
    for (SmartRefVector<Minerva::IDCluster>::const_iterator c = clusters.begin();
            c != clusters.end(); ++c) {
        const double energy = (*c)->energy();
        const double z      = (*c)->z();
        const double t      = (*c)->position();
            /* Skip low-energy clusters, normal MIP cluster has 1.7*2 MeV = 3.4 MeV */
        if (fSkipLowEnergyClusterVtxEnergy && energy < fThresholdVertexEnergy) continue;
        
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
                warning() << "Invalid view" << endmsg;
                exit(1);
                
        }
            /* Distance from this cluster to the primary vertex */
        double radius = sqrt(pow(z-z0,2) + pow(t-t0,2));
        if (radius < sphereRadius) sphereClusters.push_back(*c);


        for (std::size_t i = 0; i < radii.size(); ++i) {
            if (radius < radii[i]) vertexEnergyVector[i] += energy;
        }
                
    }

    vertexEnergyVector.swap(radii);
    
    debug() << "Exit CCProtonPi0::FilterInSphereClusters()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    
    return sphereClusters;
    
}

//==============================================================================
//  ConeBlobs
//==============================================================================
StatusCode CCProtonPi0::ConeBlobs(Minerva::PhysicsEvent *event,
                                   const SmartRef<Minerva::Vertex>& vertex) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter CCProtonPi0::ConeBlobs()" << endmsg;

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
    bool additional = false;
    bool isHough = false;
    bool isHoughApplied = false;
    
    // Get the Energy to analyze on all Clusters: Tracker, ECAL, HCAL
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
    // Save Energy Information
    event->setDoubleData( "RE_energy_Tracker", energyTracker );
    event->setDoubleData( "RE_energy_ECAL", energyECAL );
    event->setDoubleData( "RE_energy_HCAL", energyHCAL );

    debug() << " Energy to analyze " << energyTracker + energyECAL + energyHCAL << endmsg;

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


        if (try_to_recover_ && foundBlobs.size() > 2) {

            debug() << "Finding best blob combination" << endmsg;
            std::vector<double> masses;
            std::vector<double> good_masses;
            std::vector<Minerva::IDBlob*> bestBlobs;
            double dM = 1.e6;
            for (std::vector<Minerva::IDBlob*>::iterator b1 = foundBlobs.begin();
                 b1 != foundBlobs.end(); ++b1) {
                std::vector<Minerva::IDBlob*>::iterator b2 = b1+1;

                if (b2 == foundBlobs.end()) break;

                ClusterVectorInfo clusterInfo1((*b1)->clusters());
                if (clusterInfo1.GetNx() < 2) continue;
                
                double e1 = 0.0;
                double dummy1 = 0.0;
                double dummy2 = 0.0;
                double dummy3 = 0.0;
                double dummy4 = 0.0;
                m_idHoughBlob->getBlobEnergyTime((*b1),e1,dummy1,dummy2,dummy3,dummy4);
                // double t1 = (*b1)->time();

                bool goodPosition1  = m_idHoughBlob->GetStartPosition((*b1),vertex->position(),true);
                bool goodDirection1 = m_idHoughBlob->GetDirection((*b1),vertex->position());
                Gaudi::XYZVector direction1 = (*b1)->direction();
                
                if (!goodPosition1 || !goodDirection1) break;
                
                for (; b2 != foundBlobs.end(); ++b2) {

                    ClusterVectorInfo clusterInfo2((*b2)->clusters());
                    if (clusterInfo2.GetNx() < 2) continue;
                    
                    double e2 = 0.0;
                    m_idHoughBlob->getBlobEnergyTime((*b2),e2,dummy1,dummy2,dummy3,dummy4);
                    // double t2 = (*b2)->time();

                    bool goodPosition2  = m_idHoughBlob->GetStartPosition((*b2),vertex->position(),true);
                    bool goodDirection2 = m_idHoughBlob->GetDirection((*b2),vertex->position());

                    if (!goodPosition2 || !goodDirection2) continue;
                    
                    Gaudi::XYZVector direction2 = (*b2)->direction();

                    const double cos_oangle = direction1.Dot(direction2);
                    const double mgg = std::sqrt(2*e1*e2*(1-cos_oangle));
                    debug() << "\tCombination (" << std::distance(foundBlobs.begin(),b1) << ","
                              << std::distance(foundBlobs.begin(),b2) << "): " << mgg
                              << endmsg;

                    masses.push_back(mgg);

                    if (std::abs(mgg-134.9) < dM) {
                        debug() << "\t replace dM: " << dM << " " << " with " << std::abs(mgg-134.9) << endmsg;
                        dM = std::abs(mgg-134.9);
                        bestBlobs.clear();
                        bestBlobs.push_back(*b1);
                        bestBlobs.push_back(*b2);
                    }

                        /*
                    if (60.0 < mgg && mgg < 180) {
                        good_masses.push_back(mgg);

                        bestBlobs.push_back(*b1);
                        bestBlobs.push_back(*b2);
                    }
                        */
                }
                     
                
            } 

                //additional = good_masses.size() == 1;
                //if (good_masses.size() == 1) foundBlobs.swap(bestBlobs);
            additional = ( bestBlobs.size() == 2 );
            if (bestBlobs.size() == 2) foundBlobs.swap(bestBlobs);
            
            event->setContainerDoubleData("mgg_vector", masses);
            event->setContainerDoubleData("good_mgg_vector", good_masses);
        }
        
    }

    if ( !isAngleScan && !additional) { 

        /* 
            Blobs found by the AngleScan are not managed, delete them here and
            empty the container 
        */
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
    
    /* 
        if either isAngleScan or isHough is true, the blobs are finally stored
        in the TES, which are managed. Otherwise, delete them 
    */
    double minBlobSep1 = -1.;
    double minBlobSep2 = -1.;
    if (isAngleScan || isHough || additional) {
        processBlobs(event,foundBlobs);
        setPi0ParticleData(event,foundBlobs[0],foundBlobs[1],vertex);

        minBlobSep1 = CalcMinBlobSeparation(foundBlobs[0],vertex);
        minBlobSep2 = CalcMinBlobSeparation(foundBlobs[1],vertex);

        debug() << "\tmin_sep: " << minBlobSep1 << " " << minBlobSep2 << endmsg;

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
        
        std::vector<double> blob_cluster_energy1 = GetBlobClusterEnergy(foundBlobs[0],vertex);
        std::vector<double> blob_cluster_energy2 = GetBlobClusterEnergy(foundBlobs[1],vertex);
        event->setContainerDoubleData("blob_cluster_energy1", blob_cluster_energy1);
        event->setContainerDoubleData("blob_cluster_energy2", blob_cluster_energy2);
            
        
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
        m_hitTagger->applyColorTag( rejectedBlob, 0x000000 );
    }

    event->setDoubleData( "Rejected_blob_vis_energy", rejectedBlob->energy() );

//     DigitVectorTruthInfo info;
//     info.ParseTruth(rejectedBlob->getAllDigits(),fTrajectoryMap);
//     double evis = info.GetEdepByPdg(111);
//     event->setDoubleData("pi0_evis_outtime_blob", evis);
    
    debug() << "Exit CCProtonPi0::ConeBlobs()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    
    return StatusCode::SUCCESS;

}

//==============================================================================
//  AngleScanBlob()
//==============================================================================
StatusCode CCProtonPi0::AngleScanBlob(SmartRefVector<Minerva::IDCluster> idClusters,
                                       const SmartRef<Minerva::Vertex>& vertex,
                                       std::vector<Minerva::IDBlob*>& outBlobs) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter CCProtonPi0::AngleScanBlob()" << endmsg;

    const Gaudi::XYZPoint& pos = vertex->position();
    
        // main interface to create blob, this tool already exclude lowactivity clusters.
    std::vector<Minerva::IDBlob*>* preidBlobs = new std::vector<Minerva::IDBlob*>;
    m_idConeScanBlob->createIDBlobs(idClusters,preidBlobs,pos);
    
    std::vector<Minerva::IDBlob*>::iterator itBlob = preidBlobs->begin();
    for ( ; itBlob != preidBlobs->end(); itBlob ++ ){
        if ( !m_idHoughBlob->isPhoton( (*itBlob)->clusters(), pos ) ) continue;
        outBlobs.push_back(*itBlob);
    }

    debug() << "Exit CCProtonPi0::AngleScanBlob()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    
    return StatusCode::SUCCESS;
}

//==============================================================================
//  HoughBlob
//==============================================================================
StatusCode CCProtonPi0::HoughBlob(SmartRefVector<Minerva::IDCluster> idClusters,
                                   const SmartRef<Minerva::Vertex>& vertex,
                                   std::vector<Minerva::IDBlob*>& outBlobs) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter CCProtonPi0::HoughBlob()" << endmsg;

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
        debug() << " Pass Get Reference" << endmsg;	
	
        if ( !m_idHoughTool->Hough2D( idClusViewX, r, theta, ref, pos ) ) break;
        debug() << " Pass Get Hough2d" << endmsg;	
        
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
    
    debug() << "Exit CCProtonPi0::HoughBlob()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;

    return StatusCode::SUCCESS;
}

//==============================================================================
//  processBlobs
//==============================================================================
StatusCode CCProtonPi0::processBlobs( Minerva::PhysicsEvent *event, std::vector<Minerva::IDBlob*> idBlobs) const
{

    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter CCProtonPi0::processBlobs()" << endmsg;

    int count = 0;
    std::vector<Minerva::IDBlob*>::iterator itBlob = idBlobs.begin();    
    for (; itBlob != idBlobs.end(); ++itBlob){
        (*itBlob)->setHistory( Minerva::IDBlob::Unused );
        (*itBlob)->setPatRecHistory( Minerva::IDBlob::Pi0IDBlobPatRec );
        if ( (*itBlob)->nclusters( Minerva::IDCluster::X ) > 0 )  count++;
        if ( (*itBlob)->nclusters( Minerva::IDCluster::U ) > 0 )  count++;
        if ( (*itBlob)->nclusters( Minerva::IDCluster::V ) > 0 )  count++;
        debug() << " Storing Photon blob (key: " << (*itBlob)->key() << ") with " << (*itBlob)->nclusters() << " clusters "
                << endmsg;
        addObject( event, *itBlob );
    }

    if ( count != 6 ) event->filtertaglist()->setOrAddFilterTag( "is_twoDBlob", true );
    m_hitTagger->applyColorTag( (idBlobs)[0], 0xFF0000 );
    m_hitTagger->applyColorTag( (idBlobs)[1], 0x0000FF );

    debug() << "pi0 candidate" << endmsg;
    debug() << " photon 1 is blob: " << idBlobs[0]->key() << endmsg;
    debug() << " photon 2 is blob: " << idBlobs[1]->key() << endmsg;
    idBlobs[0]->setIntData("Photon1", true);
    idBlobs[1]->setIntData("Photon2", true);
//    debug() << "  photon 1 blob info: " << *(idBlobs[0]) << endmsg;

//     std::pair<int,double> result1 = OneParLineFitBlob(idBlobs[0],m_PrimaryVertex,m_MuonTrack);
//     std::pair<int,double> result2 = OneParLineFitBlob(idBlobs[1],m_PrimaryVertex,m_MuonTrack);
// 
//     double g1distance = TwoParLineFitBlobVtxDistance(idBlobs[0],m_PrimaryVertex);
//     double g2distance = TwoParLineFitBlobVtxDistance(idBlobs[1],m_PrimaryVertex);
// 
//     event->setIntData("blob_ndof_1",result1.first);
//     event->setIntData("blob_ndof_2",result2.first);
//     event->setDoubleData("blob_fval_1", result1.second);
//     event->setDoubleData("blob_fval_2", result2.second);
//     event->setDoubleData("g1blob_vtx_distance", g1distance);
//     event->setDoubleData("g2blob_vtx_distance", g2distance);
    
    
    debug() << "Exit CCProtonPi0::processBlobs()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    
    
    return StatusCode::SUCCESS;
}

//==============================================================================
//  OneParLineFitBlob
//==============================================================================
/* 
    Do a line fit for blob digits in the X view. Force the fitted line to
        go through the primary vertex
 */
std::pair<int,double> CCProtonPi0::OneParLineFitBlob(const Minerva::IDBlob* blob, 
                             const SmartRef<Minerva::Vertex>& vertex,
                             const SmartRef<Minerva::Track>& muonTrack) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter CCProtonPi0::OneParLineFitBlob" << endmsg;
    
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
          
//         debug() << "\tScan i: " << (int)(dz-zmin)/zstep << " "
//                   << z_i << " " << x_i
//                   << endmsg;

        TFitterMinuit* minuit = new TFitterMinuit;
        OneParLineFit fittingFunction(data);
        fittingFunction.SetFixedPoint(z_i,x_i);
        minuit->SetMinuitFCN(&fittingFunction);

        /* Set parameter (index,name,init,err,range) */
        minuit->SetParameter(0,"kx",1.0,0.1,-1000,1000);

        minuit->SetPrintLevel(-1);
        minuit->CreateMinimizer();

        int errno = minuit->Minimize();
        if (errno > 0) {
            std::cerr << "Fitting error: " << errno << endmsg;
            fval_vector.push_back(1.e9);
            continue;
        }

        //const double kfit = minuit->GetParameter(0);
        const double fval = fittingFunction.GetFCN();
        
        fval_min = std::min(fval,fval_min);
    }
    debug() << "\tLine fitting: chi2/ndof/total: " << fval_min << "/" << ndof
              << "/" << total 
              << endmsg;
    
    
    debug() << "Exit CCProtonPi0::OneParLineFitBlob()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    
    return std::make_pair<int,double>(ndof,fval_min/total);
}

//==============================================================================
//  TwoParLineFitBlobVtxDistance
//==============================================================================
//  Require input blob with starting position and direction
double CCProtonPi0::TwoParLineFitBlobVtxDistance(const Minerva::IDBlob* blob,
                                                  const SmartRef<Minerva::Vertex>& vertex) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter CCProtonPi0::TwoParLineFitBlobVtxDistance()" << endmsg;
    
    const double x0 = vertex->position().Z();
    const double y0 = vertex->position().X();
    
    /* Get all X blob digits via clusters since clusters have 'view' info */
    std::vector<LineFit::Point> data;
    const SmartRefVector<Minerva::IDCluster>& clusters = blob->clusters();
    for (SmartRefVector<Minerva::IDCluster>::const_iterator c = clusters.begin();
         c != clusters.end(); ++c) {

        if ((*c)->view() != Minerva::IDCluster::X) continue;

        // use cluster (x,z) position for now. Could try to use digits later
        const double x = (*c)->z();         // not a typo: z is x in 2D
        const double y = (*c)->position();  // tranverse position is y in 2D
        const double w = (*c)->energy();    // use cluster energy as weight

        data.push_back(LineFit::Point(x,y,0.0,w));

    }

    TwoParLineFit fittingFunc(data);
    TFitterMinuit* minuit = new TFitterMinuit;
    minuit->SetMinuitFCN(&fittingFunc);
    minuit->SetParameter(0,"k",  1.0,  5.0,-1000,1000);
    minuit->SetParameter(1,"b",100.0,500.0,-1.e6,1.e6);
    

    minuit->SetPrintLevel(-1);
    minuit->CreateMinimizer();

    int errno = minuit->Minimize();

    if (errno > 0) {
        std::cerr << "TwoParLineFit:: fitting error: " << errno << endmsg;
    }

    // Parameters from the fit in the form a*x + b*y + c = 0
    const double a = minuit->GetParameter(0);
    const double b = -1.0;
    const double c = minuit->GetParameter(1);

    // distance from the fitted line to the muon vertex |a*x0 + b*y0 + c|/sqrt(a*a+b*b)
    const double distance = std::abs(a*x0 + b*y0 + c)/std::sqrt(a*a + b*b);
    debug() << "TwoParLineFit: distance to vtx: " << distance << endmsg;

    debug() << "Exit CCProtonPi0::TwoParLineFitBlobVtxDistance()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    
    return distance;
}

//==============================================================================
//  ApplyAttenuationCorrection
//==============================================================================
void CCProtonPi0::ApplyAttenuationCorrection(Minerva::IDBlob* blob) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter CCProtonPi0::ApplyAttenuationCorrection()" << endmsg;

    info() << "AttenuationCorrection (before): " << blob->energy() << endmsg;

    const double x0 = blob->startPoint().x();
    const double y0 = blob->startPoint().y();
    const double z0 = blob->startPoint().z();

    const double sx = blob->direction().x();
    const double sy = blob->direction().y();
    const double sz = blob->direction().z();
    
    if (std::abs(sz) < std::cos(89.0*TMath::DegToRad())) {
        debug() << "Cannot make attenuation correction for high-angle gamma" << endmsg;
        return;
    }
    
    SmartRefVector<Minerva::IDCluster> clusters = blob->clusters();
    for (SmartRefVector<Minerva::IDCluster>::iterator c = clusters.begin();
         c != clusters.end(); ++c) {
        const double z = (*c)->z();
        const double t = (z-z0)/sz;
        const double x = x0 + sx*t;
        const double y = y0 + sy*t;
            
        m_AttenuationCorrectionTool->calibrateCluster(*c,Gaudi::XYZPoint(x,y,z));
    }

    blob->updateEnergy();
    
    info() << "AttenuationCorrection (after):  " << blob->energy() << endmsg;
    
    debug() << "Exit CCProtonPi0::ApplyAttenuationCorrection()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
}

//==============================================================================
//  CalcMinBlobSeparation
//==============================================================================
double CCProtonPi0::CalcMinBlobSeparation(const Minerva::IDBlob* blob,
                                           const SmartRef<Minerva::Vertex>& vertex) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter CCProtonPi0::CalcMinBlobSeparation()" << endmsg;
    
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
                warning() << "Invalid view" << endmsg;
                exit(1);
                
        }
        /* Distance from this cluster to the fitted vertex */
        double distance = sqrt(pow(z-z0,2) + pow(t-t0,2));
        dmin = std::min(dmin,distance);
    }
    
    debug() << "Exit CCProtonPi0::CalcMinBlobSeparation()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;

    return dmin;
}

//==============================================================================
//  GetBlobClusterEnergy
//==============================================================================
std::vector<double> CCProtonPi0::GetBlobClusterEnergy(const Minerva::IDBlob* blob, const SmartRef<Minerva::Vertex>& vertex) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() <<"Enter CCProtonPi0::GetBlobClusterEnergy()" << endmsg;

    std::vector<double> clusterEnergies;
    
    SmartRefVector<Minerva::IDCluster> clusters = blob->clusters();

    Gaudi::XYZPoint pos = vertex->position();
    debug()<<"Vertex to increasingDistanceFromVertex = ("<<pos.X()<<","<<pos.Y()<<","<<pos.Z()<<")"<<endmsg;
    
    std::sort(clusters.begin(), clusters.end(), increasingDistanceFromVertex(pos));

    for (SmartRefVector<Minerva::IDCluster>::iterator c = clusters.begin();
         c != clusters.end(); ++c) {
        clusterEnergies.push_back((*c)->energy());
        if (std::distance(clusters.begin(),c) > 4) break;
    }

    debug() << "Exit CCProtonPi0::GetBlobClusterEnergy()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    
    return clusterEnergies;
}

//==============================================================================
//  ODActivity
//==============================================================================
StatusCode CCProtonPi0::ODActivity( Minerva::PhysicsEvent *event, std::vector<Minerva::IDBlob*> idBlobs ) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() <<"Enter CCProtonPi0::ODActivity()" << endmsg;
    
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

    debug() << " Pass loop over ODClusters" << endmsg;
    
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
    
    debug() <<"Exit CCProtonPi0::ODActivity()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    
    return StatusCode::SUCCESS;
}

//==============================================================================
//  CalcDistanceFromBlobAxisToVertex
//==============================================================================
double CCProtonPi0::CalcDistanceFromBlobAxisToVertex(const Minerva::IDBlob* blob,
                                                      const SmartRef<Minerva::Vertex>& vertex) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() <<"Enter CCProtonPi0::CalcDistanceFromBlobAxisToVertex()" << endmsg;
    
    const Gaudi::XYZPoint& x0 = vertex->position();
    const Gaudi::XYZPoint& x1 = blob->startPoint();

    const Gaudi::XYZVector& s = blob->direction();
    
    const double d0 = 1000.0;
    /* Second point on the blob axis 1m away from the x1 */
    Gaudi::XYZPoint x2(x1.X() + d0*s.X(),
                       x1.Y() + d0*s.Y(),
                       x1.Z() + d0*s.Z()); 

    /* 
        Calculate the distance from the vertex to the blob axis defined by x1 and x2
            The line is defined by x1 and x2.
            The point of interest is x0
            Then the distance is:
            d = |(x0-x1)x(x0-x2)|/|x2-x1| 
    */
    Gaudi::XYZVector x01 = x0-x1;
    Gaudi::XYZVector x02 = x0-x2;
    const double dvtx = (x01.Cross(x02)).R()/(x2-x1).R();
    debug() << "\tDistance_to_vertex: " << dvtx << endmsg;

    /* Self-test using another point on the line 2m away from the x1*/
    const double d1 = 2000.0;
    Gaudi::XYZPoint x3(x1.X() + d1*s.X(),
                       x1.Y() + d1*s.Y(),
                       x1.Z() + d1*s.Z());
    
    Gaudi::XYZVector x31 = x3-x1;
    Gaudi::XYZVector x32 = x3-x2;
    const double dvtx2 = (x31.Cross(x32)).R()/(x2-x1).R();
    debug() << "\t Self-test distance: " << dvtx2 << endmsg;
    assert(dvtx2 < 1.e-6);
    
    debug() <<"Exit CCProtonPi0::CalcDistanceFromBlobAxisToVertex()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;

    return dvtx;
}

//==============================================================================
//  CalcDistanceFromVertexToExiting
//==============================================================================
double CCProtonPi0::CalcDistanceFromVertexToExiting(const Minerva::IDBlob* blob,
                                                     const SmartRef<Minerva::Vertex>& vertex) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() <<"Enter CCProtonPi0::CalcDistanceFromVertexToExiting()" << endmsg;
    
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

    debug() << "\t nstep before exiting:    " << nstep << endmsg;
    debug() << "\t Distance before exiting: " << dmax << endmsg;
    
    debug() <<"Exit CCProtonPi0::CalcDistanceFromVertexToExiting()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    
    return dmax;
}

//==============================================================================
//  CalculatedEdx
//==============================================================================
void CCProtonPi0::CalculatedEdx(const Minerva::IDBlob* blob,
                                  Minerva::PhysicsEvent* event, unsigned int blob_number,
                                  const SmartRef<Minerva::Vertex>& vertex) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() <<"Enter CCProtonPi0::CalculatedEdx()" << endmsg;

    debug() << "Calculate dEdx for blob: " << blob_number <<  " (forward,backward) ("
              << (blob->direction().Z() > 0) << "," << (blob->direction().Z() < 0) << ")"
              << endmsg;
    std::map<int, std::vector<double> > planeClusterEnergyMap;
    SmartRefVector<Minerva::IDCluster> clusters = blob->clusters();
    Gaudi::XYZPoint pos = vertex->position();
    debug()<<"Vertex to increasingDistanceFromVertex = ("<<pos.X()<<","<<pos.Y()<<","<<pos.Z()<<")"<<endmsg;
    
    std::sort(clusters.begin(), clusters.end(), increasingDistanceFromVertex(pos));
    
    for (SmartRefVector<Minerva::IDCluster>::iterator c = clusters.begin();
         c != clusters.end(); ++c) {
        const Minerva::PlaneID& planeId = (*c)->planeid();
        unsigned int plane  = planeId.plane();
        unsigned int module = planeId.module();
        const unsigned int planeNo = 2*module + plane;
        debug() << "\t(module,plane,planeNo): " << module << " " << plane << " "
                  << planeNo << endmsg;

        planeClusterEnergyMap[planeNo].push_back((*c)->energy());
    }
    
    
    std::vector<int>    blob_cluster_occupancy;
    std::vector<double> blob_cluster_energy;
    
    double total     = 0.0;
    double total1    = 0.0; // exclude the first plane
    bool doublet     = false;
    bool empty_plane = false;
    for (std::map<int, std::vector<double> >::iterator pl = planeClusterEnergyMap.begin();
         pl != planeClusterEnergyMap.end(); ++pl) {
        const int plane = pl->first;
        const std::vector<double>& clusterEnergies = pl->second;
        
        std::map<int, std::vector<double> >::iterator next_pl = pl; ++next_pl;
        const int next_plane = next_pl->first;

        const double total_plane_energy
            = std::accumulate(clusterEnergies.begin(), clusterEnergies.end(), 0.0);
        blob_cluster_occupancy.push_back(clusterEnergies.size());
        blob_cluster_energy.push_back(total_plane_energy);

        total += total_plane_energy;
        if (pl != planeClusterEnergyMap.begin()) total1 += total_plane_energy;
        
        if (clusterEnergies.size() > 1)     doublet = true;

        if (next_pl != planeClusterEnergyMap.end() &&
            std::abs(next_plane-plane) > 1) empty_plane = true;
        
        if (std::distance(planeClusterEnergyMap.begin(),pl) > 4) break;
    }

    // Less complicated logic for cluster energies from the shower back since
    // I'm not making plot for the NIM paper!
    std::vector<double> blob_rev_cluster_energy; // cluster energies from the shower back (rev. for reverse)
    for (std::map<int, std::vector<double> >::reverse_iterator pl = planeClusterEnergyMap.rbegin();
         pl != planeClusterEnergyMap.rend(); ++pl) {
        const std::vector<double>& clusterEnergies = pl->second;
        const double total_plane_energy
            = std::accumulate(clusterEnergies.begin(), clusterEnergies.end(), 0.0);
        blob_rev_cluster_energy.push_back(total_plane_energy);
    }
    

    if (blob_number == 1) {
        event->setIntData("g1dedx_nplane", planeClusterEnergyMap.size());
        event->setIntData("g1dedx_doublet", (int) doublet);
        event->setIntData("g1dedx_empty_plane", (int) empty_plane);
        event->setDoubleData("g1dedx_total", total);
        event->setDoubleData("g1dedx", total/blob_cluster_energy.size());
        event->setDoubleData("g1dedx_total1", total1);
        event->setDoubleData("g1dedx1", total1/std::max((int)blob_cluster_energy.size()-1,1));
        event->setContainerIntData("g1dedx_cluster_occupancy", blob_cluster_occupancy);
        event->setContainerDoubleData("g1dedx_cluster_energy", blob_cluster_energy);

        event->setContainerDoubleData("g1dedx_rev_cluster_energy", blob_rev_cluster_energy);

    }

    if (blob_number == 2) {
        event->setIntData("g2dedx_nplane", planeClusterEnergyMap.size());
        event->setIntData("g2dedx_doublet", (int) doublet);
        event->setIntData("g2dedx_empty_plane", (int) empty_plane);
        event->setDoubleData("g2dedx_total", total);
        event->setDoubleData("g2dedx", total/blob_cluster_energy.size());
        event->setDoubleData("g2dedx_total1", total1);
        event->setDoubleData("g2dedx1", total1/std::max((int)blob_cluster_energy.size()-1,1));
        event->setContainerIntData("g2dedx_cluster_occupancy", blob_cluster_occupancy);
        event->setContainerDoubleData("g2dedx_cluster_energy", blob_cluster_energy);

        event->setContainerDoubleData("g2dedx_rev_cluster_energy", blob_rev_cluster_energy);
    }
    
    debug() <<"Exit CCProtonPi0::CalculatedEdx()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;

}

//==============================================================================
//  CheckBlobDigitAssignment
//==============================================================================
// void CCProtonPi0::CheckBlobDigitAssignment(Minerva::PhysicsEvent* event, const Minerva::IDBlob* blob, int index) const
// {
//     debug() <<"--------------------------------------------------------------------------"<<endmsg;
//     debug() <<"Enter CCProtonPi0::CheckBlobDigitAssignment()" << endmsg;
//     
//     debug() << "Blob energy: " << blob->energy() << endmsg;
//     
//     FillTrajectoryMap();
//     
//     SmartRefVector<Minerva::IDDigit> digits = blob->getAllDigits();
// 
//     DigitVectorTruthInfo digitVectorInfo;
//     digitVectorInfo.ParseTruth(digits,fTrajectoryMap);
// 
//     const double pi0evis  = digitVectorInfo.GetEdepByPdg(111);
//     double g1evis = 0.0;
//     double g2evis = 0.0;
//     if (fGamma1) g1evis   = digitVectorInfo.GetEdepByTrackId(fGamma1->GetTrackId());
//     if (fGamma2) g2evis   = digitVectorInfo.GetEdepByTrackId(fGamma2->GetTrackId());
//     const double nevis    = digitVectorInfo.GetEdepByPdg(2112);
//     const double pevis    = digitVectorInfo.GetEdepByPdg(2212);
//     const double pipevis  = digitVectorInfo.GetEdepByPdg(+211);
//     const double pimevis  = digitVectorInfo.GetEdepByPdg(-211);
//     const double pipmevis = pipevis + pimevis;
//     const double gmevis   = digitVectorInfo.GetEdepByPdg(22);
//     const double muevis   = digitVectorInfo.GetEdepByPdg(+13)+digitVectorInfo.GetEdepByPdg(-13);
//     const double total    = digitVectorInfo.GetTotalTruthEnergy();
//     const double other    = total -(pi0evis+nevis+pevis+pipmevis+gmevis+muevis);
//     const double shared   = digitVectorInfo.GetSharedTruthEnergy();
//     const int pdg         = digitVectorInfo.GetMostEvisPdg();
//     const double frac     = digitVectorInfo.GetEdepByPdg(pdg)/total;
// 
//     event->setDoubleData(Form("g%1dpi0evis",index),     pi0evis);
//     event->setDoubleData(Form("g%1dg1evis",index),      g1evis);
//     event->setDoubleData(Form("g%1dg2evis",index),      g2evis);
//     event->setDoubleData(Form("g%1dneutronevis",index), nevis);
//     event->setDoubleData(Form("g%1dprotonevis",index),  pevis);
//     event->setDoubleData(Form("g%1dpipevis",index),     pipevis);
//     event->setDoubleData(Form("g%1dpimevis",index),     pimevis);
//     event->setDoubleData(Form("g%1dgmevis",index),      gmevis);
//     event->setDoubleData(Form("g%1dmuevis",index),      muevis);
//     event->setDoubleData(Form("g%1dotherevis",index),   other);
//     event->setDoubleData(Form("g%1dtotalevis",index),   total);
//     event->setDoubleData(Form("g%1dsharedevis",index),  shared);
//     event->setDoubleData(Form("g%1dmostevisfrac",index),frac);
//     event->setIntData(Form("g%1dmostevispdg",index),pdg);
//     
//     debug() <<"Exit CCProtonPi0::CheckBlobDigitAssignment()" << endmsg;
//     debug() <<"--------------------------------------------------------------------------"<<endmsg;
// }


//==============================================================================
//  InsideHexagon
//==============================================================================
bool CCProtonPi0::InsideHexagon(double x, double y, double w) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() <<"Enter CCProtonPi0::InsideHexagon()" << endmsg;
    
    double max = w/2;
    double min = -w/2;
    
    double deltaphi = 60.0*TMath::DegToRad();
    double phi = 0.0;
    for (int i = 0; i < 3; ++i) {
        double xtmp = std::cos(phi)*x + std::sin(phi)*y;
        phi += deltaphi;
        
        if (xtmp < min || xtmp > max) return false;
    }
    
    debug() <<"Exit CCProtonPi0::InsideHexagon()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    
    return true;
    
}

//==============================================================================
//  PreFilterPi0
//==============================================================================
StatusCode CCProtonPi0::PreFilterPi0(Minerva::PhysicsEvent *event ) const
{

    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() <<"Enter CCProtonPi0::PreFilterPi0()" << endmsg;
        //if ( event->processType() ==  Minerva::PhysicsEvent::RockParticle )  {
        //info() << " Jumping rock Muon " << endmsg;
        //return StatusCode::FAILURE;
        //}

    SmartRefVector<Minerva::IDCluster> unusedClusters
        = event->select<Minerva::IDCluster>("Unused","!LowActivity&!XTalkCandidate");

    const Gaudi::XYZPoint& vertex = event->interactionVertex()->position();
    const double inter_mod_z   = 45.0; // From geometry spreadsheet
    const double inter_plane_z = inter_mod_z/2.0;
    const double z0 = vertex.Z();

    SmartRefVector<Minerva::IDCluster> nearVertexClusters;
    for (SmartRefVector<Minerva::IDCluster>::iterator c = unusedClusters.begin();
         c != unusedClusters.end(); c++ ){
        const double c_z = (*c)->z();
        if (std::abs(c_z-z0) < 5*inter_plane_z) nearVertexClusters.push_back(*c);
    }

    std::cout << "Total near-vertex cluster: " << nearVertexClusters.size() << std::endl;

    double nearvtx_total  = 0.0;
    double nearvtx_evisX = 0.0;
    double nearvtx_evisU = 0.0;
    double nearvtx_evisV = 0.0;
    
    for (SmartRefVector<Minerva::IDCluster>::iterator c = nearVertexClusters.begin();
         c != nearVertexClusters.end(); c++ ){
        const double energy = (*c)->energy();
        nearvtx_total += energy;
        if ((*c)->view() == Minerva::IDCluster::X) nearvtx_evisX += energy;
        if ((*c)->view() == Minerva::IDCluster::U) nearvtx_evisU += energy;
        if ((*c)->view() == Minerva::IDCluster::V) nearvtx_evisV += energy;
    }

    event->setDoubleData("evis_nearvtx_total", nearvtx_total);
    event->setDoubleData("evis_nearvtx_x", nearvtx_evisX);
    event->setDoubleData("evis_nearvtx_u", nearvtx_evisU);
    event->setDoubleData("evis_nearvtx_u", nearvtx_evisV);

    
    double ntgtEvis = 0.0;
    double trkrEvis = 0.0;
    double ecalEvis = 0.0;
    double hcalEvis = 0.0;
    double total = 0.0;
    
    double ntgtEvisX = 0.0;
    double trkrEvisX = 0.0;
    double ecalEvisX = 0.0;
    double hcalEvisX = 0.0;
    double totalX = 0.0;
    
    double ntgtEvisU = 0.0;
    double trkrEvisU = 0.0;
    double ecalEvisU = 0.0;
    double hcalEvisU = 0.0;
    double totalU = 0.0;
    
    double ntgtEvisV = 0.0;
    double trkrEvisV = 0.0;
    double ecalEvisV = 0.0;
    double hcalEvisV = 0.0;
    double totalV = 0.0;

    for (SmartRefVector<Minerva::IDCluster>::iterator c = unusedClusters.begin();
         c != unusedClusters.end(); c++ ){
        const double energy = (*c)->energy();
        Minerva::IDCluster::Subdet subdet = (*c)->subdet();
        total += energy;
        if      (subdet == Minerva::IDCluster::NuclTargs) {
            ntgtEvis += energy;
            if ((*c)->view() == Minerva::IDCluster::X) ntgtEvisX += energy;
            if ((*c)->view() == Minerva::IDCluster::U) ntgtEvisU += energy;
            if ((*c)->view() == Minerva::IDCluster::V) ntgtEvisV += energy;
        }
        else if (subdet == Minerva::IDCluster::Tracker) {
            trkrEvis += energy;
            if ((*c)->view() == Minerva::IDCluster::X) trkrEvisX += energy;
            if ((*c)->view() == Minerva::IDCluster::U) trkrEvisU += energy;
            if ((*c)->view() == Minerva::IDCluster::V) trkrEvisV += energy;
        }
        else if (subdet == Minerva::IDCluster::ECAL) {
            ecalEvis += energy;
            if ((*c)->view() == Minerva::IDCluster::X) ecalEvisX += energy;
            if ((*c)->view() == Minerva::IDCluster::U) ecalEvisU += energy;
            if ((*c)->view() == Minerva::IDCluster::V) ecalEvisV += energy;
        }
        else if (subdet == Minerva::IDCluster::HCAL) {
            hcalEvis += energy;
            if ((*c)->view() == Minerva::IDCluster::X) hcalEvisX += energy;
            if ((*c)->view() == Minerva::IDCluster::U) hcalEvisU += energy;
            if ((*c)->view() == Minerva::IDCluster::V) hcalEvisV += energy;
        }
        
        else {}
    }
    
    const double otherevis = trkrEvis + ecalEvis + hcalEvis;
    
    event->setDoubleData("evis_total_x", totalX);
    event->setDoubleData("evis_ntgt_x",  ntgtEvisX);
    event->setDoubleData("evis_trkr_x",  trkrEvisX);
    event->setDoubleData("evis_ecal_x",  ecalEvisX);
    event->setDoubleData("evis_hcal_x",  hcalEvisX);
    
    event->setDoubleData("evis_total_u", totalU);
    event->setDoubleData("evis_ntgt_u",  ntgtEvisU);
    event->setDoubleData("evis_trkr_u",  trkrEvisU);
    event->setDoubleData("evis_ecal_u",  ecalEvisU);
    event->setDoubleData("evis_hcal_u",  hcalEvisU);
    
    event->setDoubleData("evis_total_v", totalV);
    event->setDoubleData("evis_ntgt_v",  ntgtEvisV);
    event->setDoubleData("evis_trkr_v",  trkrEvisV);
    event->setDoubleData("evis_ecal_v",  ecalEvisV);
    event->setDoubleData("evis_hcal_v",  hcalEvisV);
    
    event->setDoubleData("evis_total", total);
    event->setDoubleData("evis_ntgt",  ntgtEvis);
    event->setDoubleData("evis_trkr",  trkrEvis);
    event->setDoubleData("evis_ecal",  ecalEvis);
    event->setDoubleData("evis_hcal",  hcalEvis);
    event->setDoubleData("evis_other", otherevis);
    
    debug() <<"Exit CCProtonPi0::PreFilterPi0()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;

    if      ( ntgtEvis > 20)    return StatusCode::FAILURE;
    else if ( otherevis > 2000) return StatusCode::FAILURE; // energy bigger than 1.7*1.2 GeV must be ignored
    else if ( otherevis < 80 )  return StatusCode::FAILURE; // energy smaller than 80*1.2 Mev must be ignored
    else return StatusCode::SUCCESS;
    
    return StatusCode::SUCCESS;
    
}










 
