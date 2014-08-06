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
SmartRef<Minerva::Vertex>       m_PrimaryVertex;
SmartRef<Minerva::Track>        m_MuonTrack;
SmartRef<Minerva::Prong>        m_MuonProng;
SmartRef<Minerva::Particle>     m_MuonParticle;
Minerva::ProngVect    m_ProtonProngs;
Minerva::ParticleVect m_ProtonParticles;

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
//             std::cout<<z0<<std::endl;

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
    declareProperty("WriteFSParticleTable", m_writeFSParticle_Table =   false);
    declareProperty("StoreAllEvents",       m_store_all_events      =   true);
    declareProperty("DoPlausibilityCuts",   m_DoPlausibilityCuts    =   true);
    declareProperty("DoTruthMatch",         m_DoTruthMatch          =   true);
    declareProperty("MakeShortTracks",      m_makeShortTracks       =   true);
    
    declareProperty("BeamAngleBias",       m_beamAngleBias = 0.006*CLHEP::radian );
    
    declareProperty("MinMuonScore",        m_minMuonScore = 0.9 );
    declareProperty("MinProtonScore",      m_minProtonScore = 0.0 );
    
    //--------------------------------------------------------------------------
    // Arachne Colors
    //--------------------------------------------------------------------------
    // Prongs and Tracks
    declareProperty("MuonProngColor",           m_Color_muonProng       = 0x228B22); //-- green
    declareProperty("ProtonProngColor",         m_Color_protonProng     = 0x0000FF); //-- blue
    declareProperty("TrackEndProngColor",       m_Color_endPointVertex  = 0xFFA500); //-- orange
    declareProperty("UnattachedProngColor",     m_Color_unattachedProng = 0xFFFF00); //-- yellow
    
    // Blobs and Clusters
    declareProperty("BlobRejectedColor",        m_Color_RejectedBlob    = 0xFF0000); //-- red
    declareProperty("ClusterUnusedColor",       m_Color_clusterUnused   = 0xFF00FF); //-- magenta
    declareProperty("ClusterUsedColor",         m_Color_clusterUsed     = 0x008080); //-- cyan
    declareProperty("Gamma1ProngColor",         m_Color_Gamma1Prong     = 0x00FF00); //-- light green (lime)
    declareProperty("Gamma2ProngColor",         m_Color_Gamma2Prong     = 0x00BFFF); //-- light blue 
    
    // Vertex
    declareProperty("PrimaryVertexProngColor",  m_Color_primaryVertex   = 0x000000); //-- black
    declareProperty("SecondaryVertexProngColor",m_Color_secondaryVertex = 0x000000); //-- black
    declareProperty("VertexFilaColor",          m_Color_VertexFila      = 0x800080); //-- purple
    declareProperty("VertexSphereColor",        m_Color_VertexSphere    = 0xFA8072); //-- salmon 
    
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
    
    declareProperty( "TrytoRecoverBlobReco", m_TrytoRecoverBlobReco = false);
    declareProperty( "ApplyAttenuationCorrection", m_ApplyAttenuationCorrection = true);
    
    declareProperty( "UVMatchTolerance", m_UVMatchTolerance = 10.0*CLHEP::mm);
    declareProperty( "UVMatchMoreTolerance", m_UVMatchMoreTolerance = 100.0*CLHEP::mm);
    declareProperty( "AllowUVMatchWithMoreTolerance", m_AllowUVMatchWithMoreTolerance = true);
    
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
    
    //! declare common branches
    declareCommonPhysicsAnaBranches();
    declareMinosMuonBranches();
    declareGenieWeightBranches();
    
    //--------------------------------------------------------------------------
    //! Select the branches you want in your AnaTuple
    //--------------------------------------------------------------------------
    double SENTINEL = -9.9;
    
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
    declareIntEventBranch( "Cut_Event_Not_Plausible", -1 );
    declareIntEventBranch( "Cut_Event_Has_BadObject", -1 );
    declareIntEventBranch( "Cut_Vertex_None", -1 );
    declareIntEventBranch( "Cut_Vertex_Null", -1 );
    declareIntEventBranch( "Cut_Vertex_Not_Reconstructable", -1 );
    declareIntEventBranch( "Cut_Vertex_Not_Fiducial", -1 );
    declareIntEventBranch( "Cut_Muon_None",-1);
    declareIntEventBranch( "Cut_Muon_Not_Plausible",-1);
    declareIntEventBranch( "Cut_Muon_Score_Low",-1);
    declareIntEventBranch( "Cut_Muon_Charge",-1);
    declareIntEventBranch( "Cut_Vertex_Michel_Exist", -1 );
    declareIntEventBranch( "Cut_EndPoint_Michel_Exist", -1 );
    declareIntEventBranch( "Cut_secEndPoint_Michel_Exist", -1 );
    declareIntEventBranch( "Cut_Particle_None", -1 );
    declareIntEventBranch( "Cut_Proton_None", -1 );
    declareIntEventBranch( "Cut_PreFilter_Pi0", -1 );
    declareIntEventBranch( "Cut_VtxBlob", -1 );
    declareIntEventBranch( "Cut_ConeBlobs", -1 );
    
    //! Event - General reco
    declareIntEventBranch( "n_long_tracks", -1);
    declareIntEventBranch( "n_short_tracks", -1);
    declareIntEventBranch( "n_anchored_long_trk_prongs", -1 );
    declareIntEventBranch( "n_anchored_short_trk_prongs", -1 );
    declareIntEventBranch( "n_iso_trk_prongs", -1);
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
    
    //! Event - PreFilterPi0()
    declareIntEventBranch("preFilter_Result", -1);
    declareDoubleEventBranch("preFilter_rejectedEnergy", -1.0);
    
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
    declareDoubleEventBranch("Vertex_blob_energy", -9999 );
    declareDoubleEventBranch("Filament_Vertex_energy", -9999 );
    declareDoubleEventBranch("Sphere_Vertex_energy", -9999 );
    
    //! Event - ConeBlobs()
    declareDoubleEventBranch("RE_energy_Tracker", -9999 );
    declareDoubleEventBranch("RE_energy_ECAL", -9999 );
    declareDoubleEventBranch("RE_energy_HCAL", -9999 );
    
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

    declareDoubleEventBranch("g1blob_minsep", -1.0);
    declareDoubleEventBranch("g2blob_minsep", -1.0);

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
    
    //! NeutrinoInt - Muon -- Filled in setMuonData()
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
    
    //! NeutrinoInt - Proton -- Filled in setProtonData()
    declareContainerIntBranch(m_hypMeths,    "proton_trk_pat_history", 10, -1);
    declareContainerIntBranch(m_hypMeths,    "proton_kinked",      10, -1);
    declareContainerIntBranch(m_hypMeths,    "proton_odMatch",     10, -1);
    declareContainerIntBranch(m_hypMeths,    "proton_isRecoGood",     10, -1);
    declareContainerDoubleBranch(m_hypMeths, "proton_startPointX", 10, SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_startPointY", 10, SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_startPointZ", 10, SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_endPointX",   10, SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_endPointY",   10, SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_endPointZ",   10, SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_score",    10, SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_score1",   10, SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_score2",   10, SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_chi2_ndf", 10, SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_theta",    10, SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_thetaX",   10, SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_thetaY",   10, SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_phi",      10, SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_ekin",     10, SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_E",        10, SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_p",        10, SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_px",       10, SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_py",       10, SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_pz",       10, SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_p_calCorrection", 10, SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_p_visEnergy",     10, SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_p_dEdXTool",      10, SENTINEL);
    
    //! Event - Pi0 -- Filles in setPi0Data()
    // May convert to NuInt type variable
    
    declareDoubleEventBranch("pi0_px",SENTINEL);
    declareDoubleEventBranch("pi0_py",SENTINEL);
    declareDoubleEventBranch("pi0_pz",SENTINEL);
    declareDoubleEventBranch("pi0_E",SENTINEL);
    declareDoubleEventBranch("pi0_invMass", SENTINEL);
    declareDoubleEventBranch("pi0_theta", SENTINEL);
    declareDoubleEventBranch("pi0_phi",   SENTINEL);
    declareDoubleEventBranch("pi0_thetaX", SENTINEL);
    declareDoubleEventBranch("pi0_thetaY",   SENTINEL);
    declareDoubleEventBranch("pi0_openingAngle",  SENTINEL);
    declareDoubleEventBranch("pi0_cos_openingAngle", SENTINEL);
    
    declareDoubleEventBranch("gamma1_px",SENTINEL);
    declareDoubleEventBranch("gamma1_py",SENTINEL);
    declareDoubleEventBranch("gamma1_pz",SENTINEL);
    declareDoubleEventBranch("gamma1_E",SENTINEL);
    declareDoubleEventBranch("gamma1_theta",SENTINEL);
    declareDoubleEventBranch("gamma1_phi",SENTINEL);
    declareDoubleEventBranch("gamma1_dEdx",SENTINEL);
    declareDoubleEventBranch("gamma1_time",SENTINEL);
    declareDoubleEventBranch("gamma1_dist_vtx",-1.0);
    declareDoubleEventBranch("gamma1_dist_exit",-1.0);
    declareContainerDoubleEventBranch("gamma1_direction",3,SENTINEL);
    declareContainerDoubleEventBranch("gamma1_vertex",3,SENTINEL);
    declareIntEventBranch("gamma1_blob_ndigits",-1);
    declareIntEventBranch("gamma1_blob_nclusters",-1);
    declareDoubleEventBranch("gamma1_evis_trkr", SENTINEL);
    declareDoubleEventBranch("gamma1_evis_ecal", SENTINEL);
    declareDoubleEventBranch("gamma1_evis_hcal", SENTINEL);
    declareDoubleEventBranch("gamma1_evis_scal", SENTINEL);
    declareBoolEventBranch("gamma1_isGoodDirection");
    declareBoolEventBranch("gamma1_isGoodPosition");
    declareBoolEventBranch("gamma1_isGoodBlob");
    
    declareDoubleEventBranch("gamma2_px",SENTINEL);
    declareDoubleEventBranch("gamma2_py",SENTINEL);
    declareDoubleEventBranch("gamma2_pz",SENTINEL);
    declareDoubleEventBranch("gamma2_E",SENTINEL);
    declareDoubleEventBranch("gamma2_theta",SENTINEL);
    declareDoubleEventBranch("gamma2_phi",SENTINEL);
    declareDoubleEventBranch("gamma2_dEdx",SENTINEL);
    declareDoubleEventBranch("gamma2_time",SENTINEL);
    declareDoubleEventBranch("gamma2_dist_vtx",-1.0);
    declareDoubleEventBranch("gamma2_dist_exit",-1.0);
    declareContainerDoubleEventBranch("gamma2_direction",3,SENTINEL);
    declareContainerDoubleEventBranch("gamma2_vertex",3,SENTINEL);
    declareIntEventBranch("gamma2_blob_ndigits",-1);
    declareIntEventBranch("gamma2_blob_nclusters",-1);
    declareDoubleEventBranch("gamma2_evis_trkr", SENTINEL);
    declareDoubleEventBranch("gamma2_evis_ecal", SENTINEL);
    declareDoubleEventBranch("gamma2_evis_hcal", SENTINEL);
    declareDoubleEventBranch("gamma2_evis_scal", SENTINEL);
    declareBoolEventBranch("gamma2_isGoodDirection");
    declareBoolEventBranch("gamma2_isGoodPosition");
    declareBoolEventBranch("gamma2_isGoodBlob");
    
    
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
    //! Initialize Global Variables
    //--------------------------------------------------------------------------
    // SmartRefs
    m_PrimaryVertex = NULL;
    m_MuonTrack     = NULL;
    m_MuonProng     = NULL;
    m_MuonParticle  = NULL;
    // Vectors
    m_ProtonProngs.clear();
    m_ProtonParticles.clear();
    
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
    //! MAKE CUT - IF Event is NOT Plausible (MC Only)
    //--------------------------------------------------------------------------
    // default plausibility means that underlying truth is likely something to be rejected by standard Minos Match CC analysis
    if( truthEvent && m_DoPlausibilityCuts && !truthIsPlausible(truthEvent) ) {
        debug() << "This is not a plausible MC event! returning!" << endmsg;
        event->setIntData("Cut_Event_Not_Plausible",1);
        if( m_store_all_events ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    }
    
    //--------------------------------------------------------------------------
    //! MAKE CUT - IF Event has a Bad Object
    //--------------------------------------------------------------------------
    if( event->filtertaglist()->isFilterTagTrue( AnaFilterTags::BadObject() ) ) { 
        error() << "Found an event flagged with a BadObject! Refusing to analyze..." << endmsg;
        event->setIntData("Cut_Event_Has_BadObject",1);
        if( m_store_all_events ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    }
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
    m_PrimaryVertex = event->interactionVertex();
    
    //--------------------------------------------------------------------------
    //! MAKE CUT - if Interaction Vertex is NOT in Reconstructable Volume
    //--------------------------------------------------------------------------
    // "Vertex is NOT in Reconstructable Volume" means we can not run vertex-anchored
    // short tracker with a meaningful result outside of that volume.
    // Only the events that pass this CUT are used in vertex-anchored short tracker
    Gaudi::XYZPoint vtx_position = m_PrimaryVertex->position();
    if ( !FiducialPointTool->isFiducial(vtx_position, m_recoHexApothem, m_recoUpStreamZ, m_recoDownStreamZ ) ) {
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
    vtx_position = m_PrimaryVertex->position();
    if (vtx_position.z() != vtx_position.z()) {
        warning()<<"NaN vertex!"<<endmsg;
        event->setIntData("Cut_Vertex_Null",1);
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
    
    // Sanity Checks
    if(m_MuonProng == NULL){
        warning() << "m_MuonProng == NULL" << endmsg;
    }else if( m_MuonParticle == NULL){
        warning() << "m_MuonParticle == NULL" << endmsg;
    }
    
    //--------------------------------------------------------------------------
    //! MAKE CUT - IF Muon is NOT Plausible (MC Only)
    //--------------------------------------------------------------------------
    double mc_frac = -1.0;
    if ( m_DoPlausibilityCuts && !muonIsPlausible( m_MuonProng, mc_frac) ) {
        debug()<<"Muon is not plausible"<<endmsg;
        event->setIntData("Cut_Muon_Not_Plausible",1);
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
    m_hitTagger->applyColorTag(m_MuonProng, m_Color_muonProng);
    event->setTime( m_recoTimeTool->prongBestTime(m_MuonProng) );
    
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
            m_hitTagger->applyColorTag(StopPointBlob, m_Color_endPointVertex);
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

    debug() << "START: Proton Reconstruction" << endmsg;
    
    //--------------------------------------------------------------------------
    //! MAKE CUT - if Particle Creation Fails (Proton and Pion Hypotheses)
    //--------------------------------------------------------------------------
    // Get all of the primary prongs in the event
    primaryProngs = event->primaryProngs();
    
    // Create new particles with Proton and Pion hypotheses
    debug() << "Creating particles with Proton and Pion hypotheses" <<endmsg;
    bool makeParticles = createTrackedParticles(primaryProngs);
    if (!makeParticles){
        debug() << "Creation of Particles are FAILED!"<< endmsg;
        event->setIntData("Cut_Particle_None",1);
        if( m_store_all_events ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    }
    
    //--------------------------------------------------------------------------
    //! MAKE CUT - if NO Prong contains a Proton Particle
    //--------------------------------------------------------------------------
    bool foundProton = getProtonProng(primaryProngs);
    if( foundProton ) {
        for(unsigned int i = 0; i < m_ProtonProngs.size(); i++) {
            debug() << "Tag the proton's prong with bit-field = " << m_ProtonProngs[i]->typeBitsToString() << endmsg;
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
    if ( VtxBlob(event) ){
        debug()<<"Succesful VtxBlob Reconstruction!"<<endmsg;
    }else{
        event->setIntData("Cut_VtxBlob",1);
        if( m_store_all_events ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS;  
    }
    
    debug()<<"FINISH: Vertex Blob Reconstruction"<<endmsg;
    
    //--------------------------------------------------------------------------
    //! MAKE CUT - if fails ConeBlobs() -- Main Function for Pi0 Reconstruction
    //--------------------------------------------------------------------------
    debug() << "START: Cone Blob Reconstruction" << endmsg;
    if ( ConeBlobs(event) ){
        debug()<<"Succesful ConeBlobs Reconstruction!"<<endmsg;
    }else{
        event->setIntData("Cut_ConeBlobs",1);
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


    //--------------------------------------------------------------------------
    //! Color Inner Detector Unused Clusters
    //--------------------------------------------------------------------------
    for (   SmartRefVector<Minerva::IDCluster>::iterator iter_c = idClusters.begin();
            iter_c != idClusters.end(); 
            ++iter_c) 
    {
        if ( (int)(*iter_c)->history() == Minerva::IDCluster::Unused ){
            m_hitTagger->applyColorTag((*iter_c), m_Color_clusterUnused);
        }
    }
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //! count other reconstructed quantities
    //--------------------------------------------------------------------------
    int n_long_tracks = (event->select<Track>("Used:Unused","LongPatRec3View:LongPatRec2View")).size();
    int n_short_tracks = (event->select<Track>("Used:Unused","FourHitPatRec")).size();
    
    //--------------------------------------------------------------------------
    //! Calculate dead time
    //--------------------------------------------------------------------------
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
    
    event->setIntData("n_long_tracks", n_long_tracks);
    event->setIntData("n_short_tracks", n_short_tracks);
    event->setIntData("n_anchored_long_trk_prongs", n_anchored_long_trk_prongs);
    event->setIntData("n_anchored_short_trk_prongs", n_anchored_short_trk_prongs);
    event->setIntData("n_iso_trk_prongs", n_iso_trk_prongs);
    
    event->setDoubleData("muonVisibleE", muon_visible_energy );
    
    event->setDoubleData( "totalVisibleE",   totalVisibleEnergy );
    event->setDoubleData( "totalIDVisibleE", idVisibleEnergy );
    event->setDoubleData( "totalODVisibleE", odVisibleEnergy );
    
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
    bool muonFilled = setMuonData( nuInt );
    if( !muonFilled ) {
        warning()<<"Muon NTuple Branches did not filled!"<<endmsg;
        return StatusCode::FAILURE;
    }

    //--------------------------------------------------------------------------  
    //! Calculate and Set Proton Kinematics
    //--------------------------------------------------------------------------
    bool protonFilled = setProtonData(nuInt);
    if( !protonFilled ) {
        warning()<<"Proton NTuple Branches did not filled!"<<endmsg;
        return StatusCode::FAILURE;
    }
   
    //--------------------------------------------------------------------------
    //! Calculate and Set Vertex Parameters
    //--------------------------------------------------------------------------
    setVertexData(nuInt,event);
    

    //--------------------------------------------------------------------------
    //! Interaction Parameters
    //--------------------------------------------------------------------------
    //nuInt->setNeutrinoHelicity( getHelicity( mu_charge ) );
    nuInt->setInteractionCurrent( Minerva::NeutrinoInt::ChargedCurrent );
    nuInt->setInteractionType( Minerva::NeutrinoInt::UnknownInt );
    
    //--------------------------------------------------------------------------
    //! Truth Matching for Muon and Proton Prongs
    //--------------------------------------------------------------------------
    if( haveNeutrinoMC() && m_DoTruthMatch ) {
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
    if ( t_current == 1 && t_neutrinoPDG == 14 ) {
        // Atleast 1 proton, 1 pi0, 0 pi+, 0 pi-
        if(N_proton > 0 && N_pi0 == 1 && N_piplus == 0 && N_piminus == 0 && N_other == 0){
            isSignal = true;
        }
    }
    truthEvent->filtertaglist()->setOrAddFilterTag( "isSignal", isSignal );
    
    
    //--------------------------------------------------------------------------
    //! Create Final State Particle Table
    //--------------------------------------------------------------------------
    if (m_writeFSParticle_Table){
        // Using Reverse Iterator to have the table in correct order
        SmartRefVector<Minerva::TG4Trajectory>::const_reverse_iterator rit_mcpart;
        
        info() <<"--------------------------------------------------------------"<<endmsg;
        info() <<"Final State Particle Table"<<endmsg;
        if(isSignal){
            info() <<">> Marked as Signal! <<"<<endmsg;
        }
        info() <<"ID\tPDG\tParent\tPx\tPy\tPz\tE"<<endmsg;
        for (rit_mcpart = pri_trajectories.rbegin(); rit_mcpart != pri_trajectories.rend(); ++rit_mcpart) {
            Gaudi::LorentzVector temp_4p = (*rit_mcpart)->GetInitialMomentum();
            info()  <<(*rit_mcpart)->GetTrackId()<<"\t"
                    <<(*rit_mcpart)->GetPDGCode()<<"\t"
                    <<(*rit_mcpart)->GetParentId()<<"\t"
                    <<temp_4p.px()<<"\t"
                    <<temp_4p.py()<<"\t"
                    <<temp_4p.pz()<<"\t"
                    <<temp_4p.E()<<"\t"
                    <<endmsg;
        }
        info() <<"--------------------------------------------------------------"<<endmsg;
    }
   
    
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
    
//------------------------------------------------------------------------------
// interpret Events which fails the reconstructor cuts
//------------------------------------------------------------------------------
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
    debug()<<"Debug Point = 0"<<endmsg;
    for(unsigned int p = 0; p < prongs.size(); p++) {
        SmartRef<Minerva::Prong> prong = prongs[p];
        Minerva::TrackVect tracks = prong->minervaTracks();
        debug()<<"Debug Point = 1 "<<prongs.size()<<endmsg;
        for(unsigned int trk = 0; trk < tracks.size(); trk++) {
            double other_energy = 0.0;
            debug()<<"Debug Point = 2 "<<tracks.size()<<endmsg;
            std::map<const Minerva::TG4Trajectory*,double> trajMap = TruthMatcher->getTG4Trajectories(tracks[trk],other_energy);
            std::map<const Minerva::TG4Trajectory*,double>::iterator it;
            debug()<<"Debug Point = 3 "<<tracks.size()<<endmsg;
            if( trajMap.empty() ) continue;
        
            std::vector< std::pair<const Minerva::TG4Trajectory*,double> > parentVect;
            for(it = trajMap.begin(); it != trajMap.end(); it++) {
                debug()<<"Debug Point = 4 "<<endmsg;
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
            debug()<<"Debug Point = 5 "<<endmsg;
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
// Set Vertex Data 
//==============================================================================
void CCProtonPi0::setVertexData( Minerva::NeutrinoInt* nuInt, const Minerva::PhysicsEvent* event ) const 
{
    Gaudi::XYZTVector vtx_position( m_PrimaryVertex->position().x(), 
                                    m_PrimaryVertex->position().y(),
                                    m_PrimaryVertex->position().z(),
                                    event->time() );     
    nuInt->setVertex( vtx_position );
    nuInt->setScore( 1.0 );
    nuInt->setDoubleData("vtx_x", vtx_position.x() );
    nuInt->setDoubleData("vtx_y", vtx_position.y() );
    nuInt->setDoubleData("vtx_z", vtx_position.z() );
    
    // Get Vertex Module and Planes
    int vtx_module, vtx_plane;
    debug()<<"Calling getNearestPlane, vtx is "<<vtx_position.z()<<endmsg;
    getNearestPlane(vtx_position.z(), vtx_module, vtx_plane); 
    nuInt->setIntData("vtx_module", vtx_module);
    nuInt->setIntData("vtx_plane", vtx_plane);
}

//==============================================================================
// Set Muon particle data
//==============================================================================
bool CCProtonPi0::setMuonData( Minerva::NeutrinoInt* nuInt ) const 
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter CCProtonPi0::setMuonData()" << endmsg;
    
    // Sanity Check
    if (!m_MuonProng || !m_MuonParticle) {
        fatal()<< "NO Muon!!" <<endmsg;
        return false;
    }
    
    //--------------------------------------------------------------------------
    //! Get Muon Track Information
    //--------------------------------------------------------------------------
    SmartRefVector<Track>::iterator iterTrk;
    SmartRefVector<Track> muonTracks = m_MuonProng->minervaTracks();
    
    int muon_minervaTrack_types = 0;
    int muon_N_minosTracks = 0;
    int nlong = 0;
    int nshort = 0;
    
    //! Number of Short and Long Tracks
    for (iterTrk = muonTracks.begin(); iterTrk != muonTracks.end(); ++iterTrk) {
        if ((*iterTrk)->type() == Track::Long) nlong++;
        if ((*iterTrk)->type() == Track::Short) nshort++;
    }
    
    if (nlong > 0 && nshort == 0) muon_minervaTrack_types = 1;
    else if (nlong == 0 && nshort > 0) muon_minervaTrack_types = 2;
    else if (nlong >0 && nshort >0) muon_minervaTrack_types = 3;
    
    //! Get number of MINOS Tracks
    muon_N_minosTracks = m_MuonProng->minosTracks().size();   
    
    //--------------------------------------------------------------------------
    //! Get Muon Kinematics
    //--------------------------------------------------------------------------
    Gaudi::LorentzVector muon_4p;
    double SENTINEL = -9.9;
    
    int muon_charge = -99;
    int muon_roadUpstreamPlanes = -99;
    int muon_minosTrackQuality = -99;
    double muon_px = SENTINEL;
    double muon_py = SENTINEL;
    double muon_pz = SENTINEL;
    double muon_E = SENTINEL;
    double muon_p = SENTINEL;
    double muon_theta = SENTINEL;
    double muon_theta_biasUp = SENTINEL;
    double muon_theta_biasDown = SENTINEL;
    double muon_muScore = SENTINEL; 
    double muon_qp = SENTINEL; 
    double muon_qpqpe = SENTINEL;
    double muon_E_shift = SENTINEL;
    double muon_roadUpstreamEnergy = SENTINEL;
    
    // 4-Momentum and Angle Information
    muon_4p    = m_MuonParticle->momentumVec();
    muon_px    = muon_4p.px();
    muon_py    = muon_4p.py();
    muon_pz    = muon_4p.pz();    
    muon_E     = muon_4p.E();
    muon_p     = muon_4p.P();
    muon_theta = m_coordSysTool->thetaWRTBeam(muon_4p);
    muon_theta_biasUp = m_coordSysTool->thetaWRTBeam(muon_4p,m_beamAngleBias) - muon_theta;
    muon_theta_biasDown = m_coordSysTool->thetaWRTBeam(muon_4p, -1.0*m_beamAngleBias) - muon_theta;
    
    // Muon Score
    muon_muScore = m_MuonParticle->score();
    
    // Muon Q/P / Q/P-Error
    muon_qpqpe = MuonUtils->minosQPQPE(m_MuonProng);
    if (m_MuonProng->MinosTrack()) {
        SmartRef<MinosRecoTrack> minosTrack = m_MuonProng->minosTracks()[0];
        muon_qp = minosTrack->qp();
        muon_minosTrackQuality = minosTrack->quality();      
    }
    else if ( !m_MuonProng->MinosStub() ) {
        warning()<<"No MINOS track or stub!  How is this a muon?"<<endmsg;
    }
    
    // Muon energy in road upstream info
    muon_roadUpstreamEnergy = AnaToolUtils->energyInTheRoadUpstream(m_MuonProng, muon_roadUpstreamPlanes);  
    
    // Muon Energy Shift
    //      the shifts are so small compared to MINOS-matched momentum that we can 
    //      approximate the momentum shift as the energy shift (with ~0.2% at 1.5 GeV/c)
    muon_E_shift = MuonUtils->calculateMomentumCorrection(m_MuonProng);
    
    //Get Muon Charge
    MuonUtils->muonCharge(m_MuonProng,muon_charge); 
    
    //--------------------------------------------------------------------------
    // Debugging: Check values
        debug()<<"P4(Muon) = ( "
        <<muon_4p.px()<<", "
        <<muon_4p.py()<<", "
        <<muon_4p.pz()<<", "
        <<muon_4p.E()<<" )"
        <<endmsg;
    //--------------------------------------------------------------------------
    
    debug()<<"Filling Muon Ntuple Variables"<<endmsg;
    
    //--------------------------------------------------------------------------
    //! Fill Muon Branches
    //--------------------------------------------------------------------------
    nuInt->setLeptonEnergy( muon_4p );
    
    nuInt->setIntData("muon_minervaTrack_types", muon_minervaTrack_types);
    nuInt->setIntData("muon_N_minosTracks", muon_N_minosTracks);
    nuInt->setIntData("muon_minosTrackQuality", muon_minosTrackQuality);
    nuInt->setIntData("muon_roadUpstreamPlanes", muon_roadUpstreamPlanes);
    nuInt->setIntData("muon_charge",muon_charge);
    nuInt->setDoubleData("muon_roadUpstreamEnergy", muon_roadUpstreamEnergy);
    nuInt->setDoubleData("muon_px",muon_px);
    nuInt->setDoubleData("muon_py",muon_py);
    nuInt->setDoubleData("muon_pz",muon_pz);
    nuInt->setDoubleData("muon_E",muon_E);
    nuInt->setDoubleData("muon_p",muon_p);
    nuInt->setDoubleData("muon_theta",muon_theta);
    nuInt->setDoubleData("muon_theta_biasUp",muon_theta_biasUp);
    nuInt->setDoubleData("muon_theta_biasDown",muon_theta_biasDown);
    nuInt->setDoubleData("muon_muScore", muon_muScore);
    nuInt->setDoubleData("muon_qp",muon_qp );
    nuInt->setDoubleData("muon_qpqpe",muon_qpqpe);
    nuInt->setDoubleData("muon_E_shift",muon_E_shift);

    fillMinosMuonBranches(nuInt, m_MuonProng);

    debug() << "Exit CCProtonPi0::setMuonData()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    
    return true;
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
    
    // Check if the prongs are odMatch
    Minerva::EventMgr* mgr = getEventMgr(prongs[0]);
    Minerva::ODClusterVect odClusters = mgr->select<Minerva::ODCluster>("Unused","!LowActivity&!XTalkCandidate");
    m_odMatchTool->classifyProngs(prongs,odClusters);
    
    // Loop over prongs
    for(unsigned int p = 0; p < prongs.size(); p++) {
        debug() << "The prong of bit-field = " << prongs[p]->typeBitsToString() << endmsg;
    
        // Make sure the prong is not a bad object
        bool pass = true; std::string tag = "BadObject";
        if( prongs[p]->filtertaglist()->filterTagExists(tag) ) {
            if( prongs[p]->filtertaglist()->checkFilterTag(tag,pass) ) {
                error() << "This prong = " << prongs[p]->typeBitsToString() << " has been flagged as a \"BadObject\", skipping!" << endmsg;
                continue;
            }
        }
        
        // Skipped prongs with particles
        debug() << "The prong has n particles = " << prongs[p]->particles().size() << endmsg;
        if( !prongs[p]->particles().empty() ) continue;
    
        // Select the particle hypotheses candidates and particle tools
        std::vector< Minerva::Particle::ID > hypotheses;
        IParticleMakerTool::NameAliasListType toolsToUse;

        // Do NOT use prongs which have OD Match
        if( !prongs[p]->OdMatch() ) {
            hypotheses.push_back( Minerva::Particle::Proton );
            hypotheses.push_back( Minerva::Particle::Pion );
            toolsToUse.push_back( std::make_pair("dEdXTool","dEdXTool") );
        }

        // Make particles
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
// Uses Global Variables: m_ProtonProngs and m_ProtonParticles
//==============================================================================
bool CCProtonPi0::getProtonProng(   Minerva::ProngVect& primaryProngs) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter CCProtonPi0::getProtonProng()" << endmsg;
  
    debug() << "N(primaryProngs) =  " << primaryProngs.size() << endmsg;
    
    // Initialize
    bool isProtonExist = false;
    
    // Get All Proton Candidates
    for(unsigned int p = 0; p < primaryProngs.size(); p++) {
        debug() <<"Checking prong "<<p<<endmsg;
        // Skip Muon Prong
        if( primaryProngs[p]->filtertaglist()->filterTagExists("PrimaryMuon") ){
            debug() << "Muon Prong, skipping this prong! "<< endmsg;
            continue;
        }
        
        // Temp Storage for ProngVect, Single Prong
        Minerva::ProngVect tmpProngs;
        SmartRef<Minerva::Prong> prong       = (Minerva::Prong*)NULL;
        SmartRef<Minerva::Particle> particle = (Minerva::Particle*)NULL;
        
        // Push current Prong to temp ProngVect
        tmpProngs.push_back( primaryProngs[p] );
        
        // Find Proton using m_protonUtils
        bool isProton = m_protonUtils->findProtonProng(tmpProngs,prong,particle); 
        
        if( isProton ) {
            debug() <<"Found a proton prong!"<< endmsg;
            debug() <<"Proton Particle Score: " << particle->score() << endmsg;
            prong->filtertaglist()->setOrAddFilterTag("PrimaryProton",true);
            particle->filtertaglist()->setOrAddFilterTag("PrimaryProton",true);
            prong->updateBestParticle(particle);
            m_hitTagger->applyColorTag(prong, m_Color_protonProng);
            m_ProtonProngs.push_back( prong );
            m_ProtonParticles.push_back( particle );
            isProtonExist = true;
        }
    }

    debug() << "Found "<<m_ProtonParticles.size()<<" Proton Candidates!"<<endmsg;
   
    debug() << "Exit CCProtonPi0::getProtonProng()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    
    return isProtonExist;
}

//==============================================================================
// Set proton particle data
//==============================================================================
bool CCProtonPi0::setProtonData( Minerva::NeutrinoInt* nuInt ) const 
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter CCProtonPi0::setProtonData()" << endmsg;
    
    debug() <<"Inside setProtonData()"<<endmsg;
    debug() <<"m_ProtonProngs.size() = "<<m_ProtonProngs.size()
            <<" m_ProtonParticles.size() = "<<m_ProtonParticles.size()<<endmsg;
            
    // Sanity Check
    if ( m_ProtonProngs.size() == 0 ) {
        warning()<< "m_ProtonProngs is empty" <<endmsg;
        return false;
    }else if ( m_ProtonParticles.size() == 0){
        warning()<< "m_ProtonParticles is empty" <<endmsg;
        return false;
    }else if ( m_ProtonProngs.size() != m_ProtonParticles.size()){
        warning()<< "m_ProtonProngs and m_ProtonParticles NOT SAME SIZE" <<endmsg;
        return false;
    }
    
    // Declare Variables
    double SENTINEL = -9.9;
    std::vector<double> p_calCorrection(10,SENTINEL);
    std::vector<double> p_visEnergyCorrection(10,SENTINEL);
    std::vector<double> p_dedx(10,SENTINEL);
    
    std::vector<double> proton_theta(10,SENTINEL);
    std::vector<double> proton_thetaX(10,SENTINEL);
    std::vector<double> proton_thetaY(10,SENTINEL);
    std::vector<double> proton_phi(10,SENTINEL);
    
    std::vector<double> E(10,SENTINEL);
    std::vector<double> px(10,SENTINEL);
    std::vector<double> py(10,SENTINEL);
    std::vector<double> pz(10,SENTINEL);
    std::vector<double> p(10,SENTINEL);
    
    std::vector<double> ekin(10,SENTINEL);
    std::vector<double> enu(10,SENTINEL);
    std::vector<double> Q2(10,SENTINEL);
    
    std::vector<double> proton_end_x(10,SENTINEL);
    std::vector<double> proton_end_y(10,SENTINEL);
    std::vector<double> proton_end_z(10,SENTINEL);
    
    std::vector<double> proton_start_x(10,SENTINEL);
    std::vector<double> proton_start_y(10,SENTINEL);
    std::vector<double> proton_start_z(10,SENTINEL);
    
    std::vector<int> kinked(10,-1);
    std::vector<int> odMatch(10,-1);
    std::vector<int> isRecoGood(10,-1);
    
    std::vector<double> score(10,SENTINEL);
    std::vector<double> score1(10,SENTINEL);
    std::vector<double> score2(10,SENTINEL);
    std::vector<double> chi2(10,SENTINEL);
    
    // Get Vertex Z Position
    double vertexZ = m_PrimaryVertex->position().z();
    
    // Loop over all Proton Candidates
    for(unsigned int i = 0; i < m_ProtonProngs.size(); i++) {
        SmartRef<Minerva::Prong> prong = m_ProtonProngs[i];
        SmartRef<Minerva::Particle> particle = m_ProtonParticles[i];
        
        double theta = prong->minervaTracks().front()->theta();
        double phi   = prong->minervaTracks().front()->phi();
        
        // Mark Prong if Particle Energy == 0
        if ( particle->momentumVec().E() == 0.0 ){
            debug()<<"particle energy from dEdX = 0";
            isRecoGood[i] = -1;
        }else{
            isRecoGood[i] = 1;
        }
    
        proton_theta[i]  = m_coordSysTool->thetaWRTBeam(particle->momentumVec());
        proton_thetaX[i] = m_coordSysTool->thetaXWRTBeam(particle->momentumVec());
        proton_thetaY[i] = m_coordSysTool->thetaYWRTBeam(particle->momentumVec());
        proton_phi[i]    = m_coordSysTool->phiWRTBeam(particle->momentumVec());
        
        //----------------------------------------------------------------------
        // Debugging: Check values
            debug()<<"Theta(Proton) = "<<theta
            <<" Phi(Proton) = "<<phi
            <<endmsg;
        //----------------------------------------------------------------------
        
        //----------------------------------------------------------------------
        // Debugging: Check values
            debug()<<"Before m_energyCorrectionTool->getCorrectedEnergy"<<endmsg;
            debug()<<"P4(Proton) = ( "
            <<particle->momentumVec().Px()<<", "
            <<particle->momentumVec().Py()<<", "
            <<particle->momentumVec().Pz()<<", "
            <<particle->momentumVec().E()<<" )"
            <<endmsg;
        //----------------------------------------------------------------------
    
        debug()<<"Entering m_energyCorrectionTool->getCorrectedEnergy..."<<endmsg;
        Gaudi::LorentzVector dEdXprotonfourVec;
        StatusCode sc = m_energyCorrectionTool->getCorrectedEnergy(prong,particle,vertexZ,dEdXprotonfourVec);
        
        // if m_energyCorrectionTool->getCorrectedEnergy returns SUCCESS
        if( sc ){
            debug()<<"m_energyCorrectionTool->getCorrectedEnergy return SUCCESS"<<endmsg;
            debug()<<"particle->setMomentumVec(dEdXprotonfourVec)"<<endmsg;
            particle->setMomentumVec(dEdXprotonfourVec);
        } else{
            debug()<<"m_energyCorrectionTool->getCorrectedEnergy return FAILURE"<<endmsg;
            debug()<<"dEdXprotonfourVec = particle->momentumVec()"<<endmsg;
            dEdXprotonfourVec = particle->momentumVec();
        }
        
        p_dedx[i] = sqrt( dEdXprotonfourVec.E()*dEdXprotonfourVec.E() - MinervaUnits::M_proton*MinervaUnits::M_proton );
        
        //----------------------------------------------------------------------
        // Debugging: Check values
            debug()<<"Before correctProtonProngEnergy"<<endmsg;
            debug()<<"P4(Proton) = ( "
            <<particle->momentumVec().Px()<<", "
            <<particle->momentumVec().Py()<<", "
            <<particle->momentumVec().Pz()<<", "
            <<particle->momentumVec().E()<<" )"
            <<endmsg;
        //----------------------------------------------------------------------
        
        bool isCorrected = correctProtonProngEnergy(prong,p_calCorrection[i],p_visEnergyCorrection[i]);
        
        // if correctProtonProngEnergy() returns false ==> use default particle 4-Momentum
        if ( !isCorrected){
            warning()<<"correctProtonProngEnergy() returned false!"<<endmsg;
            p_calCorrection[i] = sqrt(  particle->momentumVec().Px() * particle->momentumVec().Px() +
                                        particle->momentumVec().Py() * particle->momentumVec().Py() +
                                        particle->momentumVec().Pz() * particle->momentumVec().Pz() );
        }
        //----------------------------------------------------------------------
        // Debugging: Check values
            debug()<<"p_calCorrection[i] = "<<p_calCorrection[i]<<endmsg;
            debug()<<"p_visEnergyCorrection[i] = "<<p_visEnergyCorrection[i]<<endmsg;
        //----------------------------------------------------------------------
        
        E[i]  = sqrt( p_calCorrection[i]*p_calCorrection[i] + MinervaUnits::M_proton*MinervaUnits::M_proton );
        px[i] = p_calCorrection[i]*sin(theta)*cos(phi);
        py[i] = p_calCorrection[i]*sin(theta)*sin(phi);
        pz[i] = p_calCorrection[i]*cos(theta);
        p[i]  = sqrt( px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i] );
        
        //----------------------------------------------------------------------
        // Debugging: Check values
            debug()<<"After correctProtonProngEnergy"<<endmsg;
            debug()<<"P4(Proton) = ( "
            <<px[i]<<", "
            <<py[i]<<", "
            <<pz[i]<<", "
            <<E[i]<<" )"
            <<endmsg;
        //----------------------------------------------------------------------
    
        if( prong->minervaTracks().front()->direction() == Minerva::Track::Backward && pz[i] > 0. ) pz[i] *= -1.0;
        Gaudi::LorentzVector protonfourVec(px[i],py[i],pz[i],E[i]);
    
        ekin[i] = protonfourVec.E() - MinervaUnits::M_p;
    
        proton_end_x[i]   = (prong->minervaTracks().back())->lastState().x();
        proton_end_y[i]   = (prong->minervaTracks().back())->lastState().y();
        proton_end_z[i]   = (prong->minervaTracks().back())->lastState().z();
    
        proton_start_x[i] = (prong->minervaTracks().back())->firstState().x();
        proton_start_y[i] = (prong->minervaTracks().back())->firstState().y();
        proton_start_z[i] = (prong->minervaTracks().back())->firstState().z();
    
        kinked[i]  = (int)prong->Kinked();
        odMatch[i] = (int)prong->OdMatch();
        score[i]   = particle->score();
        score1[i]  = particle->getDoubleData("score1");
        score2[i]  = particle->getDoubleData("score2");
        chi2[i]    = particle->getDoubleData("chi2_ndf");
        
        //----------------------------------------------------------------------
        // Debugging: Check values
            debug()<<"P4(Proton) = ( "
            <<px[i]<<", "
            <<py[i]<<", "
            <<pz[i]<<", "
            <<E[i]<<" )"
            <<" score = "
            <<score[i]
            <<endmsg;
        //----------------------------------------------------------------------
        
        //----------------------------------------------------------------------
        // Debugging: Check values
            debug()<<"ekin[i] = "<<ekin[i]<<endmsg;
            debug()<<"kinked[i] = "<<kinked[i]<<endmsg;
            debug()<<"odMatch[i] = "<<odMatch[i]<<endmsg;
            debug()<<"isRecoGood[i] = "<<isRecoGood[i]<<endmsg;
        //----------------------------------------------------------------------
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
    
    nuInt->setContainerDoubleData("proton_px",px);
    nuInt->setContainerDoubleData("proton_py",py);
    nuInt->setContainerDoubleData("proton_pz",pz);
    nuInt->setContainerDoubleData("proton_E",E);
    nuInt->setContainerDoubleData("proton_p",p);
    
    nuInt->setContainerIntData("proton_kinked",kinked);
    nuInt->setContainerIntData("proton_odMatch",odMatch);
    nuInt->setContainerIntData("proton_isRecoGood",isRecoGood);
    
    nuInt->setContainerDoubleData("proton_score",score);
    nuInt->setContainerDoubleData("proton_score1",score1);
    nuInt->setContainerDoubleData("proton_score2",score2);
    nuInt->setContainerDoubleData("proton_chi2_ndf",chi2);
    
    nuInt->setContainerDoubleData("proton_theta",proton_theta);
    nuInt->setContainerDoubleData("proton_thetaX",proton_thetaX);
    nuInt->setContainerDoubleData("proton_thetaY",proton_thetaY);
    nuInt->setContainerDoubleData("proton_phi",proton_phi);
    
    nuInt->setContainerDoubleData("proton_ekin",ekin);
    
    debug() << "Exit CCProtonPi0::setProtonData()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    
    return true;
}

//==============================================================================
// Correct Proton Energy
//==============================================================================
bool CCProtonPi0::correctProtonProngEnergy(  SmartRef<Minerva::Prong>& protonProng, 
                                                double& p_calCorrection, 
                                                double& p_visEnergyCorrection ) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter CCProtonPi0::correctProtonProngEnergy()" << endmsg;
    
    //-- initialize
    bool isCorrected = true;
    double E_calCorrection = 0.0;
    double E_visEnergyCorrection = 0.0;
    double calEnergy = 0.0;
    double visEnergy = 0.0;
    p_calCorrection         = 0.0;
    p_visEnergyCorrection   = 0.0;
    
        
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
        for(unsigned int clus = 0; clus < clusters.size(); clus++){
            visEnergy += clusters[clus]->energy(); 
        }
    }
    //----------------------------------------------------------------------
    // Debugging: Check values
        debug() <<"The fuzz blobs energies!"<<endmsg;
        debug() <<"calEnergy = "<<calEnergy<<endmsg;
        debug() <<"visEnergy = "<<visEnergy<<endmsg;
    //----------------------------------------------------------------------
    
    //-- get the prong first track
    SmartRef<Minerva::Track> track = (*protonProng->minervaTracks().begin());
    
    //-- update the particles four momentum
    Minerva::ParticleVect particles = protonProng->particles();
    debug() <<"Prong N(Particles) = "<<particles.size()<<endmsg;
    for(unsigned int part = 0; part < particles.size(); part++) {
        if( particles[part]->idcode() != Minerva::Particle::Proton ) continue;
    
        Gaudi::XYZTVector fourMomentum = particles[part]->momentumVec();
        debug() << "particle idcode = " << particles[part]->idcode() << ", and four momentum = " << fourMomentum << endmsg;
    
        E_calCorrection = fourMomentum.E() + calEnergy;
        E_visEnergyCorrection = fourMomentum.E() + visEnergy;
        
        // Break if Energy is lower than particle rest mass
        if( E_calCorrection < particles[part]->mass()){
            warning() << "E_calCorrection < Particle Mass "<<E_calCorrection<<" < "<<particles[part]->mass()<<endmsg;
            isCorrected = false;
            break;
        }else if( E_visEnergyCorrection < particles[part]->mass()){
            warning() << " E_visEnergyCorrection < Particle Mass "<<E_visEnergyCorrection<<" < "<<particles[part]->mass()<<endmsg;
            isCorrected = false;
            break;
        }
        
        p_calCorrection += sqrt( pow(E_calCorrection,2) - pow(particles[part]->mass(),2) );
        debug() << "update energy using calorimetric correction = " << E_calCorrection << endmsg;
        debug() << "update momentum using calorimetric correction = " << p_calCorrection << endmsg;
        
        p_visEnergyCorrection += sqrt( pow(E_visEnergyCorrection,2) - pow(particles[part]->mass(),2) );     
        debug() << "update energy using visible energy correction = " << E_visEnergyCorrection << endmsg;
        debug() << "update momentum using visible correction = " << p_visEnergyCorrection << endmsg;
    }     
    
    debug() << "Exit CCProtonPi0::correctProtonProngEnergy()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    
    // p_visEnergyCorrection == 0 leads to a 0 4-Momentum
    return isCorrected;
}

//==============================================================================
// setPi0Data
//==============================================================================

bool CCProtonPi0::setPi0Data(   Minerva::PhysicsEvent *event, 
                                        Minerva::IDBlob* idblob1, 
                                        Minerva::IDBlob* idblob2) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "CCProtonPi0::setPi0Data()" << endmsg;
    
    const Gaudi::XYZPoint& vtx_position = m_PrimaryVertex->position();
    
    // Check Gamma1 Quality
    bool goodPosition1  = m_idHoughBlob->GetStartPosition(idblob1, vtx_position, true );
    bool goodDirection1 = m_idHoughBlob->GetDirection(idblob1, vtx_position );
    bool isGoodBlob1 = false;
    if (goodPosition1 && goodDirection1) isGoodBlob1 = true;
    
    // Check Gamma2 Quality
    bool goodPosition2  = m_idHoughBlob->GetStartPosition(idblob2, vtx_position, true );
    bool goodDirection2 = m_idHoughBlob->GetDirection(idblob2, vtx_position );
    bool isGoodBlob2 = false;
    if (goodPosition2 && goodDirection2) isGoodBlob2 = true;
    
    event->filtertaglist()->setOrAddFilterTag("gamma1_isGoodDirection", goodDirection1);
    event->filtertaglist()->setOrAddFilterTag("gamma1_isGoodPosition",  goodPosition1);
    event->filtertaglist()->setOrAddFilterTag("gamma1_isGoodBlob", isGoodBlob1);
    
    event->filtertaglist()->setOrAddFilterTag("gamma2_isGoodBlob", isGoodBlob2);
    event->filtertaglist()->setOrAddFilterTag("gamma2_isGoodDirection", goodDirection2);
    event->filtertaglist()->setOrAddFilterTag("gamma2_isGoodPosition",  goodPosition2);
    
    // Process only if we have 2 Good Blobs
    if (!isGoodBlob1 || !isGoodBlob2){
        debug() << "Don't have 2 Good Blobs! -- CANNOT Set Pi0 Particle Data!" <<endmsg;
        return false;
    }

    if (m_ApplyAttenuationCorrection) {
        ApplyAttenuationCorrection(idblob1);
        ApplyAttenuationCorrection(idblob2);
    }
    
    // Get Blob Energy and Time
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

    //--------------------------------------------------------------------------
    //    Make sure Gamma1 is the more energetic one 
    //--------------------------------------------------------------------------
    if (g2energy > g1energy) {
        std::swap(goodPosition1, goodPosition2);  /* Swap variables already assigned */
        std::swap(goodDirection1,goodDirection2);
        std::swap(isGoodBlob1,   isGoodBlob2);
        std::swap(g1energy,g2energy);
        std::swap(idblob1,idblob2);               /* And the blobs themselves */
    }
    
    // Get Blob Time
    double time1 = idblob1->time();
    double time2 = idblob2->time();

    // Get Blob dEdx
    double dEdx1 = 0.0;
    double dEdx2 = 0.0;
    m_idHoughBlob->idBlobdEdx( idblob1, dEdx1 );
    m_idHoughBlob->idBlobdEdx( idblob2, dEdx2 );

    // Get Blob Directions
    Gaudi::XYZVector direction1 = idblob1->direction(); 
    Gaudi::XYZVector direction2	= idblob2->direction();

    // Save Gamma1 and Gamma2 Starting Position (blob vertex)
    std::vector<double> position1;
    std::vector<double> position2;
    position1.push_back(idblob1->startPoint().x());
    position1.push_back(idblob1->startPoint().y());
    position1.push_back(idblob1->startPoint().z());
    
    position2.push_back(idblob2->startPoint().x());
    position2.push_back(idblob2->startPoint().y());
    position2.push_back(idblob2->startPoint().z());

    // Save Gamma1 and Gamma2 Direction
    std::vector<double> direc_1;
    std::vector<double> direc_2;
    direc_1.push_back(direction1.x());
    direc_1.push_back(direction1.y());
    direc_1.push_back(direction1.z());
    
    direc_2.push_back(direction2.x());
    direc_2.push_back(direction2.y());
    direc_2.push_back(direction2.z());
    
    // Calculate Gamma1 and Gamma2 3-Momentums
    TVector3 g1mom(direction1.x(),direction1.y(),direction1.z());
    TVector3 g2mom(direction2.x(),direction2.y(),direction2.z());
    g1mom *= g1energy;
    g2mom *= g2energy;
    
    // Calculate Distance: Vertex - Blob
    const double gamma1_dist_vtx = CalcDistanceFromBlobAxisToVertex(idblob1);
    const double gamma2_dist_vtx = CalcDistanceFromBlobAxisToVertex(idblob2);
    
    // Calculate Distance: Blob - Exit
    const double gamma1_dist_exit = CalcDistanceFromVertexToExiting(idblob1);
    const double gamma2_dist_exit = CalcDistanceFromVertexToExiting(idblob2);
    
    // Calculate Pi0 3-Momentum
    TVector3 pimom = g1mom + g2mom;
    
    // Calculate Opening Angle
    const double openingAngle       = (g1mom.Angle(g2mom))*TMath::RadToDeg();
    const double cos_openingAngle   = direction1.Dot(direction2);
    
    // Calculate invariant Mass of Pi0
    const double invMass = std::sqrt(2*g1energy*g2energy*(1-cos_openingAngle));
    
    // Calculate dEdX for Blobs
    CalculatedEdx(idblob1,event,1);
    CalculatedEdx(idblob2,event,2);
    
    //--------------------------------------------------------------------------
    //    Fill Branches
    //--------------------------------------------------------------------------
    
    // Pi0 Information
    event->setDoubleData("pi0_openingAngle", openingAngle);
    event->setDoubleData("pi0_cos_openingAngle", cos_openingAngle );
    event->setDoubleData("pi0_invMass", invMass);
    event->setDoubleData("pi0_px",pimom.x());
    event->setDoubleData("pi0_py",pimom.y());
    event->setDoubleData("pi0_pz",pimom.z());
    event->setDoubleData("pi0_E",g1energy+g2energy);
    event->setDoubleData("pi0_theta",  pimom.Theta()*TMath::RadToDeg());
    event->setDoubleData("pi0_phi",    pimom.Phi()*TMath::RadToDeg());
    event->setDoubleData("pi0_thetaX", std::atan2(pimom.X(),pimom.Z())*TMath::RadToDeg());
    event->setDoubleData("pi0_thetaY", std::atan2(pimom.Y(),pimom.Z())*TMath::RadToDeg());
    
    // Gamma1 Information
    event->setDoubleData("gamma1_px",g1mom.x());
    event->setDoubleData("gamma1_py",g1mom.y());
    event->setDoubleData("gamma1_pz",g1mom.z());
    event->setDoubleData("gamma1_E",g1energy);
    event->setDoubleData("gamma1_theta",g1mom.Theta()*TMath::RadToDeg());
    event->setDoubleData("gamma1_phi",  g1mom.Phi()*TMath::RadToDeg());
    event->setDoubleData("gamma1_dEdx", dEdx1 );
    event->setDoubleData("gamma1_time", time1 );
    event->setDoubleData("gamma1_dist_vtx",gamma1_dist_vtx);
    event->setDoubleData("gamma1_dist_exit",gamma1_dist_exit);
    event->setContainerDoubleData("gamma1_direction", direc_1 );
    event->setContainerDoubleData("gamma1_vertex", position1 );
    event->setIntData("gamma1_blob_ndigits",  idblob1->getAllDigits().size());
    event->setIntData("gamma1_blob_nclusters",idblob1->nclusters());
    event->setDoubleData("gamma1_evis_trkr", g1trkrevis);
    event->setDoubleData("gamma1_evis_ecal", g1ecalevis);
    event->setDoubleData("gamma1_evis_hcal", g1hcalevis);
    event->setDoubleData("gamma1_evis_scal", g1scalevis);
    
    // Gamma2 Information
    event->setDoubleData("gamma2_px",g2mom.x());
    event->setDoubleData("gamma2_py",g2mom.y());
    event->setDoubleData("gamma2_pz",g2mom.z());
    event->setDoubleData("gamma2_E",g2energy);
    event->setDoubleData("gamma2_theta",g2mom.Theta()*TMath::RadToDeg());
    event->setDoubleData("gamma2_phi",  g2mom.Phi()*TMath::RadToDeg());
    event->setDoubleData("gamma2_dEdx", dEdx2 );
    event->setDoubleData("gamma2_time", time2 );
    event->setDoubleData("gamma2_dist_vtx",gamma2_dist_vtx);
    event->setDoubleData("gamma2_dist_exit",gamma2_dist_exit);
    event->setContainerDoubleData("gamma2_direction", direc_2 );
    event->setContainerDoubleData("gamma2_vertex", position2 );
    event->setIntData("gamma2_blob_nclusters",idblob2->nclusters());
    event->setIntData("gamma2_blob_ndigits",  idblob2->getAllDigits().size());
    event->setDoubleData("gamma2_evis_trkr", g2trkrevis);
    event->setDoubleData("gamma2_evis_ecal", g2ecalevis);
    event->setDoubleData("gamma2_evis_hcal", g2hcalevis);
    event->setDoubleData("gamma2_evis_scal", g2scalevis);

    debug() << "Exit CCProtonPi0::setPi0Data()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;

    return true;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
//
//  CCPi0 Functions -- Blob Tools
//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>


//==============================================================================
//  VtxBlob
//==============================================================================
bool CCProtonPi0::VtxBlob(Minerva::PhysicsEvent *event) const
{
    
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter CCProtonPi0::VtxBlob()" << endmsg;
    
    const Gaudi::XYZPoint& vtx_position = m_PrimaryVertex->position();

    //  --  Vertex blob
    double vertex_energy = 0;
    double vertex_energy_sphere = 0;
    double vertex_energy_filament= 0;

    if ( m_filamentVertex ){
        
        SmartRefVector<Minerva::IDCluster> analyFilaClusters 
            = event->select<Minerva::IDCluster>("Unused","!LowActivity&!XTalkCandidate");
    
        SmartRefVector<Minerva::IDCluster> filaClusters;
        m_blobUtils->fillProximateClusterVec(vtx_position,analyFilaClusters,filaClusters,
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
            m_hitTagger->applyColorTag( vtxFilaBlob, m_Color_VertexFila );
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
        SmartRefVector<Minerva::IDCluster> sphereClusters = FilterInSphereClusters(unusedClusters,m_maxSearchD,radii);
    
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
            m_hitTagger->applyColorTag( vtxSphereBlob, m_Color_VertexSphere );
        }
    } // if ( m_sphereVertex )
    
    vertex_energy = vertex_energy_sphere + vertex_energy_filament;
    event->setDoubleData( "Vertex_blob_energy", vertex_energy );
    event->setDoubleData( "Filament_Vertex_energy", vertex_energy_filament );
    event->setDoubleData( "Sphere_Vertex_energy", vertex_energy_sphere );
  
    debug() << "Exit CCProtonPi0::VtxBlob()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;

  return true;

}

//==============================================================================
//  FilterInSphereClusters()
//==============================================================================
SmartRefVector<Minerva::IDCluster> CCProtonPi0::FilterInSphereClusters( const SmartRefVector<Minerva::IDCluster>& clusters,
                                     const double sphereRadius,
                                     std::vector<double>& radii) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter CCProtonPi0::FilterInSphereClusters()" << endmsg;
    
    //--------------------------------------------------------------------------
    // Debugging: Check values
    const Gaudi::XYZPoint& vtx_position = m_PrimaryVertex->position();
    debug() <<"m_PrimaryVertex->position() = "<<endmsg;
    debug() <<"("<<vtx_position.x()<<","<<vtx_position.y()<<","<<vtx_position.z()<<")"<< endmsg;
    //--------------------------------------------------------------------------
    
    const double x0 = vtx_position.X();
    const double y0 = vtx_position.Y();
    const double z0 = vtx_position.Z();
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
bool CCProtonPi0::ConeBlobs( Minerva::PhysicsEvent *event ) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter CCProtonPi0::ConeBlobs()" << endmsg;
    
    // Initialize Bool Variables
    bool isAngleScan        = false;
    bool isAngleScanApplied = false;
    bool additional         = false;
    bool isHough            = false;
    bool isHoughApplied     = false;
    bool isPi0DataSet       = false;
    
    // Get Vertex Position
    const Gaudi::XYZPoint& vtx_position = m_PrimaryVertex->position();
    
    // variables to clean up the clusters
    SmartRefVector<Minerva::IDCluster> preidClusters 
            = event->select<Minerva::IDCluster>("Unused","!LowActivity&!XTalkCandidate");
    SmartRefVector<Minerva::IDCluster> usableClusters;
    SmartRefVector<Minerva::IDCluster> outTimeClusters;
    
    // ---------------------------------------------------------------------
    // Debugging
    debug()<<"Cluster History in ConeBlobs()"<<endmsg;
    // Check Cluster History
    // Loop over Clusters
    for (   SmartRefVector<Minerva::IDCluster>::iterator iter_c = preidClusters.begin();
            iter_c != preidClusters.end(); 
            ++iter_c) 
    {
        if( (int)(*iter_c)->history() == Minerva::IDCluster::Used ){
            debug()<<"Used"<<endmsg;
        }else if ( (int)(*iter_c)->history() == Minerva::IDCluster::Unused ){
            debug()<<"Unused"<<endmsg;
        }
    }
    // ---------------------------------------------------------------------    
    
    // get origin of muon
    Gaudi::XYZTVector muon_position = m_MuonParticle->startPos();
    
    //--------------------------------------------------------------------------
    // Get the Energy to Analyze on all Clusters: Tracker, ECAL, HCAL and
    // Fill usableClusters
    //--------------------------------------------------------------------------
    double energyTracker    = 0.0;
    double energyECAL       = 0.0;
    double energyHCAL       = 0.0;
    double EnergyToAnalyze  = 0.0;
    
    // Loop over all UNUSED Clusters
    SmartRefVector<Minerva::IDCluster>::iterator it_clus;
    for ( it_clus = preidClusters.begin(); it_clus != preidClusters.end(); ++it_clus){
        
        if ( (*it_clus)->pe()/(*it_clus)->iddigs() <= 3 ) continue;
        
        if ( (*it_clus)->subdet() ==  Minerva::IDCluster::Tracker ) energyTracker += (*it_clus)->energy();
        if ( (*it_clus)->subdet() ==  Minerva::IDCluster::ECAL )    energyECAL += (*it_clus)->energy();
        if ( (*it_clus)->subdet() ==  Minerva::IDCluster::HCAL )    energyHCAL += (*it_clus)->energy();
        
        // Fill usableClusters 
        // Include clusters close to MUON vertex time < 25 ns
        if ( std::abs( (*it_clus)->time() -  muon_position.T() ) < m_rejectedClustersTime ) {
            usableClusters.push_back(*it_clus);
        } else {
            outTimeClusters.push_back(*it_clus);
        }
    }
    EnergyToAnalyze = energyTracker + energyECAL + energyHCAL;
    debug() << "Energy to analyze = "<< EnergyToAnalyze << endmsg;
    
    // Return False if the Energy is Higher than Limit
    if ( EnergyToAnalyze >= m_energyHoughlimit ){
        debug() <<"Energy to analyze is higher than the limit! "
                <<EnergyToAnalyze<<" >= "<<m_energyHoughlimit<<endmsg;
        debug() <<"Exit CCProtonPi0::ConeBlobs()" << endmsg;
        debug() <<"--------------------------------------------------------------------------"<<endmsg;
        return false;
    }
    
    // Save Energy Information to NTuples
    event->setDoubleData("RE_energy_Tracker", energyTracker );
    event->setDoubleData("RE_energy_ECAL", energyECAL );
    event->setDoubleData("RE_energy_HCAL", energyHCAL );


    //--------------------------------------------------------------------------
    // Analyze usableClusters using AngleScan Class
    //--------------------------------------------------------------------------
    std::vector<Minerva::IDBlob*> foundBlobs;
    unsigned int nblob_anglescan = 0;
    unsigned int nblob_hough     = 0;
    
    // Create AngleScan Object
    AngleScan angleScanAlg(usableClusters,vtx_position);
    angleScanAlg.AllowUVMatchWithMoreTolerance(m_AllowUVMatchWithMoreTolerance);
    angleScanAlg.SetUVMatchTolerance(m_UVMatchTolerance);
    angleScanAlg.DoReco();

    event->setIntData("anglescan_ncandx", angleScanAlg.GetNxCandidate());   // Number of Shower Candidates in X
    event->setIntData("anglescan_ncand",  angleScanAlg.GetNCandidate());    // Number of Shower Candidates

    std::vector<int> candx_nc;
    std::vector<int> candx_nd;
    // Get Shower Candidates in X View
    const std::vector<SmartRefVector<Minerva::IDCluster> >& xshowerCand = angleScanAlg.GetXShowerCandVector();
    // Loop over Shower Candidates in X View
    for (   std::vector<SmartRefVector<Minerva::IDCluster> >::const_iterator s = xshowerCand.begin();
            s != xshowerCand.end(); 
            ++s) 
    {
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
    for (   std::vector<SmartRefVector<Minerva::IDCluster> >::const_iterator s = showerCand.begin();
            s != showerCand.end(); 
            ++s) 
    {
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
    for (   std::vector<Minerva::IDBlob*>::const_iterator b = angleScanBlobs.begin();
            b != angleScanBlobs.end(); 
            ++b) 
    {
        if ( !m_idHoughBlob->isPhoton((*b)->clusters(),vtx_position)) continue;
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
    for (   std::vector<Minerva::IDBlob*>::const_iterator b = foundBlobs.begin();
            b != foundBlobs.end(); 
            ++b) 
    {
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
        

    nblob_anglescan = foundBlobs.size();
    isAngleScan = foundBlobs.size() == 2;
    isAngleScanApplied = true;


    if (m_TrytoRecoverBlobReco && foundBlobs.size() > 2) {

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

            bool goodPosition1  = m_idHoughBlob->GetStartPosition((*b1),vtx_position,true);
            bool goodDirection1 = m_idHoughBlob->GetDirection((*b1),vtx_position);
            Gaudi::XYZVector direction1 = (*b1)->direction();
            
            if (!goodPosition1 || !goodDirection1) break;
            
            for (; b2 != foundBlobs.end(); ++b2) {

                ClusterVectorInfo clusterInfo2((*b2)->clusters());
                if (clusterInfo2.GetNx() < 2) continue;
                
                double e2 = 0.0;
                m_idHoughBlob->getBlobEnergyTime((*b2),e2,dummy1,dummy2,dummy3,dummy4);
                // double t2 = (*b2)->time();

                bool goodPosition2  = m_idHoughBlob->GetStartPosition((*b2),vtx_position,true);
                bool goodDirection2 = m_idHoughBlob->GetDirection((*b2),vtx_position);

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
        
        StatusCode sc = HoughBlob(usableClusters, foundBlobs );
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
        isPi0DataSet = setPi0Data(event,foundBlobs[0],foundBlobs[1]);

        //----------------------------------------------------------------------
        // Return if setPi0Data() Fails
        // Returning false ensures no execution of interpretEvent()
        // ---------------------------------------------------------------------
        if (!isPi0DataSet){
            debug() << "Pi0 Data was not Set, Returning!..."<<endmsg;
            debug() << "Exit CCProtonPi0::ConeBlobs()" << endmsg;
            debug() <<"--------------------------------------------------------------------------"<<endmsg;
            return false;
        }
        
        minBlobSep1 = CalcMinBlobSeparation(foundBlobs[0]);
        minBlobSep2 = CalcMinBlobSeparation(foundBlobs[1]);

        debug() << "\tmin_sep: " << minBlobSep1 << " " << minBlobSep2 << endmsg;

        std::vector<int> blob_nc;
        std::vector<int> blob_ncx;
        std::vector<int> blob_ncu;
        std::vector<int> blob_ncv;

        std::vector<int> blob_nd;
        std::vector<int> blob_ndx;
        std::vector<int> blob_ndu;
        std::vector<int> blob_ndv;
        for (   std::vector<Minerva::IDBlob*>::const_iterator b = foundBlobs.begin();
                b != foundBlobs.end(); 
                ++b) 
        {
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
        
        std::vector<double> blob_cluster_energy1 = GetBlobClusterEnergy(foundBlobs[0]);
        std::vector<double> blob_cluster_energy2 = GetBlobClusterEnergy(foundBlobs[1]);
        event->setContainerDoubleData("blob_cluster_energy1", blob_cluster_energy1);
        event->setContainerDoubleData("blob_cluster_energy2", blob_cluster_energy2);
        
        ODActivity(event,foundBlobs);
            
    } else {
        for (   std::vector<Minerva::IDBlob*>::iterator b = foundBlobs.begin();
                b != foundBlobs.end(); 
                ++b) 
        { 
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
        m_hitTagger->applyColorTag( rejectedBlob, m_Color_RejectedBlob ); // red
    }

    event->setDoubleData( "Rejected_blob_vis_energy", rejectedBlob->energy() );
    
    debug() << "Exit CCProtonPi0::ConeBlobs()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    
    return isPi0DataSet;

}


//==============================================================================
//  HoughBlob
//==============================================================================
StatusCode CCProtonPi0::HoughBlob(SmartRefVector<Minerva::IDCluster> idClusters,
                                   std::vector<Minerva::IDBlob*>& outBlobs) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter CCProtonPi0::HoughBlob()" << endmsg;

    //--------------------------------------------------------------------------
    // Debugging: Check values
    const Gaudi::XYZPoint& vtx_position = m_PrimaryVertex->position();
    debug() <<"m_PrimaryVertex->position() = "<<endmsg;
    debug() <<"("<<vtx_position.x()<<","<<vtx_position.y()<<","<<vtx_position.z()<<")"<< endmsg;
    //--------------------------------------------------------------------------
    
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
	
        if ( !m_idHoughTool->Hough2D( idClusViewX, r, theta, ref, vtx_position ) ) break;
        debug() << " Pass Get Hough2d" << endmsg;	
        
        if ( !m_idHoughBlob->Create2dHTSeed( idClusViewX, idHTSeed, r, theta, ref, spX, spZ ) ) break;
        debug() << " Pass Get Create2dHT " << endmsg;	
        
        Gaudi::XYZPoint startpoint(spX, 0, spZ);
        Gaudi::XYZVector direction(-1/tan(theta*CLHEP::pi/180),0,1);

        if ( !m_idHoughBlob->PseudoCone(idHTSeed, idClusViewX, direction, vtx_position ) ) continue;
        debug() << " Pass PseudoCone " << endmsg;
        
        if ( !m_idHoughBlob->XUVMatch(idHTSeed, idClusViewU, idClusViewV, 50 ) ) continue;
        debug() << " Pass XUV Match - First " << endmsg;	
        
        if ( !m_idHoughBlob->isPhoton(idHTSeed, vtx_position) ) continue;
        
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
            m_idHoughBlob->GetStartPosition( *itBlob, vtx_position, true );
            m_idHoughBlob->GetDirection( *itBlob, vtx_position );
            outBlobs.push_back(*itBlob);
        }
    }

    for ( it_clus = idClusViewX.begin(); it_clus != idClusViewX.end(); it_clus++ )
            m_idHoughBlob->AddClusterInsideCone( *it_clus, outBlobs, vtx_position );

    for ( it_clus = idClusViewU.begin(); it_clus != idClusViewU.end(); it_clus++ )
            m_idHoughBlob->AddClusterInsideCone( *it_clus, outBlobs, vtx_position );

    for ( it_clus = idClusViewV.begin(); it_clus != idClusViewV.end(); it_clus++ )
            m_idHoughBlob->AddClusterInsideCone( *it_clus, outBlobs, vtx_position );


    info() << " Hough Transform is done! " << endmsg;
    
    debug() << "Exit CCProtonPi0::HoughBlob()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;

    return StatusCode::SUCCESS;
}

//==============================================================================
//  processBlobs
//==============================================================================
void CCProtonPi0::processBlobs( Minerva::PhysicsEvent *event, std::vector<Minerva::IDBlob*> idBlobs) const
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

    if ( count != 6 ){ 
        event->filtertaglist()->setOrAddFilterTag( "is_twoDBlob", true );
    }
    m_hitTagger->applyColorTag( (idBlobs)[0], m_Color_Gamma1Prong ); // light green
    m_hitTagger->applyColorTag( (idBlobs)[1], m_Color_Gamma2Prong ); // light blue

    debug() << "pi0 candidate" << endmsg;
    debug() << " photon 1 is blob: " << idBlobs[0]->key() << endmsg;
    debug() << " photon 2 is blob: " << idBlobs[1]->key() << endmsg;
    idBlobs[0]->setIntData("Photon1", true);
    idBlobs[1]->setIntData("Photon2", true);

    debug() << "Exit CCProtonPi0::processBlobs()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
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
double CCProtonPi0::CalcMinBlobSeparation(const Minerva::IDBlob* blob) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter CCProtonPi0::CalcMinBlobSeparation()" << endmsg;

    //--------------------------------------------------------------------------
    // Debugging: Check values
    const Gaudi::XYZPoint& vtx_position = m_PrimaryVertex->position();
    debug() <<"m_PrimaryVertex->position() = "<<endmsg;
    debug() <<"("<<vtx_position.x()<<","<<vtx_position.y()<<","<<vtx_position.z()<<")"<< endmsg;
    //--------------------------------------------------------------------------
    
    const double x0 = vtx_position.X();
    const double y0 = vtx_position.Y();
    const double z0 = vtx_position.Z();
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
std::vector<double> CCProtonPi0::GetBlobClusterEnergy(const Minerva::IDBlob* blob) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() <<"Enter CCProtonPi0::GetBlobClusterEnergy()" << endmsg;
    
    //--------------------------------------------------------------------------
    // Debugging: Check values
    Gaudi::XYZPoint vtx_position = m_PrimaryVertex->position();
    debug() <<"m_PrimaryVertex->position() = "<<endmsg;
    debug() <<"("<<vtx_position.x()<<","<<vtx_position.y()<<","<<vtx_position.z()<<")"<< endmsg;
    //--------------------------------------------------------------------------

    std::vector<double> clusterEnergies;
    
    SmartRefVector<Minerva::IDCluster> clusters = blob->clusters();

    debug()<<"Vertex to increasingDistanceFromVertex = ("<<vtx_position.X()<<","<<vtx_position.Y()<<","<<vtx_position.Z()<<")"<<endmsg;
    
    std::sort(clusters.begin(), clusters.end(), increasingDistanceFromVertex(vtx_position));

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
double CCProtonPi0::CalcDistanceFromBlobAxisToVertex( const Minerva::IDBlob* blob ) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() <<"Enter CCProtonPi0::CalcDistanceFromBlobAxisToVertex()" << endmsg;
    
    //--------------------------------------------------------------------------
    // Debugging: Check values
    Gaudi::XYZPoint vtx_position = m_PrimaryVertex->position();
    debug() <<"m_PrimaryVertex->position() = "<<endmsg;
    debug() <<"("<<vtx_position.x()<<","<<vtx_position.y()<<","<<vtx_position.z()<<")"<< endmsg;
    //--------------------------------------------------------------------------
    
    const Gaudi::XYZPoint& x0 = m_PrimaryVertex->position();
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
double CCProtonPi0::CalcDistanceFromVertexToExiting(const Minerva::IDBlob* blob) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() <<"Enter CCProtonPi0::CalcDistanceFromVertexToExiting()" << endmsg;
    
    //--------------------------------------------------------------------------
    // Debugging: Check values
    const Gaudi::XYZPoint& vtx_position = m_PrimaryVertex->position();
    debug() <<"m_PrimaryVertex->position() = "<<endmsg;
    debug() <<"("<<vtx_position.x()<<","<<vtx_position.y()<<","<<vtx_position.z()<<")"<< endmsg;
    //--------------------------------------------------------------------------
    
    const Gaudi::XYZVector& s = blob->direction();
    
    double d_step = 50.0; /* Accuracy on dmax = dstep/2 */
    double dmax = 0.0;    /* Distance from vertex to where the blob would exit */
    unsigned int nstep = 0;
    for (double d = 0.0; ; d += d_step) {
        dmax = d;
        ++nstep;
        Gaudi::XYZPoint r = vtx_position;
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
void CCProtonPi0::CalculatedEdx(    const Minerva::IDBlob* blob,
                                    Minerva::PhysicsEvent* event, 
                                    unsigned int blob_number) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() <<"Enter CCProtonPi0::CalculatedEdx()" << endmsg;
    
    //--------------------------------------------------------------------------
    // Debugging: Check values
    Gaudi::XYZPoint vtx_position = m_PrimaryVertex->position();
    debug() <<"m_PrimaryVertex->position() = "<<endmsg;
    debug() <<"("<<vtx_position.x()<<","<<vtx_position.y()<<","<<vtx_position.z()<<")"<< endmsg;
    //--------------------------------------------------------------------------

    debug() << "Calculate dEdx for blob: " << blob_number <<  " (forward,backward) ("
              << (blob->direction().Z() > 0) << "," << (blob->direction().Z() < 0) << ")"
              << endmsg;
    std::map<int, std::vector<double> > planeClusterEnergyMap;
    SmartRefVector<Minerva::IDCluster> clusters = blob->clusters();

    debug()<<"Vertex to increasingDistanceFromVertex = ("<<vtx_position.X()<<","<<vtx_position.Y()<<","<<vtx_position.Z()<<")"<<endmsg;
    
    std::sort(clusters.begin(), clusters.end(), increasingDistanceFromVertex(vtx_position));
    
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
bool CCProtonPi0::PreFilterPi0(Minerva::PhysicsEvent *event) const
{

    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() <<"Enter CCProtonPi0::PreFilterPi0()" << endmsg;
    
    //--------------------------------------------------------------------------
    // Debugging: Check values
    const Gaudi::XYZPoint& vtx_position = m_PrimaryVertex->position();
    debug() <<"m_PrimaryVertex->position() = "<<endmsg;
    debug() <<"("<<vtx_position.x()<<","<<vtx_position.y()<<","<<vtx_position.z()<<")"<< endmsg;
    //--------------------------------------------------------------------------

    // Visible Energy Limits for the PreFilter [MeV]
    // Other = Tracker + ECAL + HCAL
    const double max_Evis_Target = 20;
    const double min_Evis_Other = 80;   // energy smaller than 80*1.2 Mev must be ignored
    const double max_Evis_Other = 2000; // energy bigger than 1.7*1.2 GeV must be ignored
    
    const double inter_mod_z   = 45.0; // From geometry spreadsheet
    const double inter_plane_z = inter_mod_z/2.0;
    const double z0 = vtx_position.Z();
    
    // Use only UNUSED clusters
    SmartRefVector<Minerva::IDCluster> unusedClusters;
    unusedClusters = event->select<Minerva::IDCluster>("Unused","!LowActivity&!XTalkCandidate");

    //--------------------------------------------------------------------------
    // Find Near-Vertex Clusters
    //--------------------------------------------------------------------------
    double nearvtx_total = 0.0;
    double nearvtx_evisX = 0.0;
    double nearvtx_evisU = 0.0;
    double nearvtx_evisV = 0.0;
    
    // Save all near-vertex clusters
    SmartRefVector<Minerva::IDCluster> nearVertexClusters;
    for (   SmartRefVector<Minerva::IDCluster>::iterator iter_cluster = unusedClusters.begin();
            iter_cluster != unusedClusters.end();
            ++iter_cluster ){
            
        const double cluster_z = (*iter_cluster)->z();
        if (std::abs(cluster_z-z0) < 5*inter_plane_z){ 
            nearVertexClusters.push_back(*iter_cluster);
        }
    }

    debug() << "Total Number of near-vertex clusters: " << nearVertexClusters.size() << endmsg;

    // Save Visible Energy of near-vertex clusters for each view
    for (   SmartRefVector<Minerva::IDCluster>::iterator iter_cluster = nearVertexClusters.begin();
            iter_cluster != nearVertexClusters.end(); 
            ++iter_cluster ){
            
        const double energy = (*iter_cluster)->energy();
        Minerva::IDCluster::View view = (*iter_cluster)->view();
        
        nearvtx_total += energy;
        if (view == Minerva::IDCluster::X) nearvtx_evisX += energy;
        if (view == Minerva::IDCluster::U) nearvtx_evisU += energy;
        if (view == Minerva::IDCluster::V) nearvtx_evisV += energy;
    }

    event->setDoubleData("evis_nearvtx_total", nearvtx_total);
    event->setDoubleData("evis_nearvtx_x", nearvtx_evisX);
    event->setDoubleData("evis_nearvtx_u", nearvtx_evisU);
    event->setDoubleData("evis_nearvtx_u", nearvtx_evisV);


    //--------------------------------------------------------------------------
    // Save Visible Energy of all UNUSED Hits for each sub detector and view
    // Target, Tracker, ECAL, HCAL, TOTAL
    //--------------------------------------------------------------------------
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

    for (   SmartRefVector<Minerva::IDCluster>::iterator iter_cluster = unusedClusters.begin();
            iter_cluster != unusedClusters.end(); 
            ++iter_cluster ){
            
        const double energy = (*iter_cluster)->energy();
        Minerva::IDCluster::Subdet subdet = (*iter_cluster)->subdet();
        Minerva::IDCluster::View view = (*iter_cluster)->view();
//         if ((*iter_cluster)->view() == Minerva::IDCluster::X) ntgtEvisX += energy;
        
        // Total Visible Energy
        total += energy;
        // Visible Energy in SubDetector = Target
        if (subdet == Minerva::IDCluster::NuclTargs) {
            ntgtEvis += energy;
            if (view == Minerva::IDCluster::X) ntgtEvisX += energy;
            if (view == Minerva::IDCluster::U) ntgtEvisU += energy;
            if (view == Minerva::IDCluster::V) ntgtEvisV += energy;
        }
        // Visible Energy in SubDetector = Tracker
        else if (subdet == Minerva::IDCluster::Tracker) {
            trkrEvis += energy;
            if (view == Minerva::IDCluster::X) trkrEvisX += energy;
            if (view == Minerva::IDCluster::U) trkrEvisU += energy;
            if (view == Minerva::IDCluster::V) trkrEvisV += energy;
        }
        // Visible Energy in SubDetector = ECAL
        else if (subdet == Minerva::IDCluster::ECAL) {
            ecalEvis += energy;
            if (view == Minerva::IDCluster::X) ecalEvisX += energy;
            if (view == Minerva::IDCluster::U) ecalEvisU += energy;
            if (view == Minerva::IDCluster::V) ecalEvisV += energy;
        }
        // Visible Energy in SubDetector = HCAL
        else if (subdet == Minerva::IDCluster::HCAL) {
            hcalEvis += energy;
            if (view == Minerva::IDCluster::X) hcalEvisX += energy;
            if (view == Minerva::IDCluster::U) hcalEvisU += energy;
            if (view == Minerva::IDCluster::V) hcalEvisV += energy;
        }
        else{
            debug() <<"Cluster SubDetector does not found!"<<endmsg;
        }
    }
    
    // Visible Energy except Target Region
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


    
    //--------------------------------------------------------------------------
    // preFilter_Result
    // 0 = Passes Filter
    // 1 = Rejected by Target Filter
    // 2 = Rejected by Max Other
    // 3 = Rejected by Min Other
    //--------------------------------------------------------------------------
    
    // Check Target Region
    if ( ntgtEvis > max_Evis_Target){
        debug() <<"Evis(target) bigger than max value!"<<endmsg;
        debug() <<"Evis(target) = "<<ntgtEvis<<" max =  "<<max_Evis_Target<<endmsg;
        event->setIntData("preFilter_Result",1);
        return false;
    }
    // Check Other Regions (Tracker + ECAL + HCAL)
    if ( otherevis > max_Evis_Other){
        debug() <<"Evis(other) bigger than max value!"<<endmsg;
        debug() <<"Evis(other) = "<<otherevis<<" max =  "<<max_Evis_Other<<endmsg;
        event->setIntData("preFilter_Result",2);
        event->setDoubleData("preFilter_rejectedEnergy",otherevis);
        return false;
    }else if( otherevis < min_Evis_Other ){
        debug() <<"Evis(other) lower than min value!"<<endmsg;
        debug() <<"Evis(other) = "<<otherevis<<" min =  "<<min_Evis_Other<<endmsg;
        event->setIntData("preFilter_Result",3);
        event->setDoubleData("preFilter_rejectedEnergy",otherevis);
        return false;
    }

    // Passed all Tests
    debug() <<"Evis(target) = "<<ntgtEvis<<" Evis(other) = "<<otherevis<<endmsg;
    debug() <<"max_Target = "<<max_Evis_Target<<" min_Other =  "<<min_Evis_Other<<" max_Other =  "<<max_Evis_Other<<endmsg;
    event->setIntData("preFilter_Result",0);
    
    return true;
    
}










 
