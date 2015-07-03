/*
    See CCProtonPi0.h header for Class Information
*/
#ifndef CCProtonPi0_cpp 
#define CCProtonPi0_cpp 1

#include <cassert>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <numeric>

// Local Libraries
#include "CCProtonPi0.h"
#include "Pi0Reco/AngleScan.h"
#include "Pi0Reco/ClusterVectorInfo.h"
#include "Pi0Reco/DigitVectorTruthInfo.h"
#include "Pi0Reco/OneParLineFit.h"
#include "Pi0Reco/TwoParLineFit.h"

#include "CCProtonPi0/IHoughBlob.h"
#include "CCProtonPi0/IHoughTool.h"

// ROOT Libraries
#include "TString.h"
#include "TDatabasePDG.h"
#include "TVector3.h"
#include "TMath.h"
#include "TRandom3.h"

// Gaudi
#include "GaudiKernel/PhysicalConstants.h"

// Minerva Analysis Framework
#include "AnaUtils/ICCPionIncUtils.h"
#include "AnaUtils/IMCTrackTool.h"
#include "AnaUtils/IProtonUtils.h"
#include "AnaUtils/MCTrack.h"
#include "BadChannels/IGetDeadTime.h"
#include "BlobFormation/IBlobCreatorUtils.h"
#include "BlobFormation/IIDAnchoredBlobCreator.h"
#include "CalTools/IGetCalAttenuation.h"
#include "EnergyRecTools/IEnergyCorrectionTool.h"
#include "EnergyRecTools/IExtraEnergyTool.h"
#include "Event/TimeSlice.h"
#include "G4Material.hh"
#include "G4Navigator.hh"
#include "G4ThreeVector.hh"
#include "G4VPhysicalVolume.hh"
#include "GeoUtils/IMinervaCoordSysTool.h"
#include "GiGaCnv/IGiGaGeomCnvSvc.h"
#include "MinervaDet/DeDetector.h"
#include "MinervaDet/IGeomUtilSvc.h"
#include "MinervaUtils/IHitTaggerTool.h"
#include "MinervaUtils/IMinervaMathTool.h"
#include "MinervaUtils/IMinervaObjectAssociator.h"
#include "ODDet/DeOuterDetector.h"
#include "ParticleMaker/IParticleMakerTool.h"
#include "ParticleMaker/IParticleTool.h"
#include "ProngMaker/IMichelTool.h"
#include "ProngMaker/IODProngClassificationTool.h"
#include "ProngMaker/IProngClassificationTool.h"
#include "RecInterfaces/IAnchoredTrackFormation.h"
#include "RecInterfaces/IFiducialPointTool.h"
#include "RecInterfaces/IRecoObjectTimeTool.h"
#include "RecUtils/Cone.h"
#include "RecUtils/IConeUtilsTool.h"
#include "RecoStudies/PDG.h"
#include "RecoStudies/IVertexEnergyStudyTool.h"
#include "VertexCreation/IVertexFitter.h"

//==============================================================================
// Global variables
//==============================================================================
/*
   Global Variables used instead of class data members so that they can be
   assigned in 'const' methods. In general, we can use 'mutable' members,
   but it does not seem to work in the framework 
*/
Gaudi::LorentzVector m_muon_4P;
Gaudi::LorentzVector m_proton_4P;
Gaudi::LorentzVector m_pi0_4P;
SmartRef<Minerva::Vertex>   m_PrimaryVertex;
SmartRef<Minerva::Track>    m_MuonTrack;
SmartRef<Minerva::Prong>    m_MuonProng;
SmartRef<Minerva::Particle> m_MuonParticle;
SmartRef<Minerva::IDBlob>   m_Pi0Blob1;
SmartRef<Minerva::IDBlob>   m_Pi0Blob2;
Minerva::ProngVect    m_ProtonProngs;
Minerva::ParticleVect m_ProtonParticles;



// Counters for Functions - Debugging Purposes
double N_tagTruth;
double N_reconstructEvent;


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
    debug() <<"Enter CCProtonPi0::CCProtonPi0() -- Default Constructor" << endmsg;
    
    declareInterface<IInteractionHypothesis>(this);
    // Mandatory declaration of analysis signature: CCProtonPi0
    m_anaSignature = "CCProtonPi0";
    
	 // Private Properties
    declareProperty("WriteFSParticleTable", m_writeFSParticle_Table =   false);
    declareProperty("StoreAllEvents",       m_store_all_events      =   true);
    declareProperty("DoPlausibilityCuts",   m_DoPlausibilityCuts    =   true);
    declareProperty("DoTruthMatch",         m_DoTruthMatch          =   true);
    declareProperty("MakeShortTracks",      m_makeShortTracks       =   true);
    declareProperty("RecoNoProtonEvents",   m_reconstruct_NoProtonEvents = true);
    declareProperty("ApplyExtraMichelCuts", m_applyExtraMichelCuts = false);
    
    declareProperty("BeamAngleBias",       m_beamAngleBias = 0.006*CLHEP::radian );
    
    declareProperty("MinMuonScore",        m_minMuonScore = 0.9 );
    declareProperty("MinProtonScore",      m_minProtonScore = 0.0 );
	
    declareProperty("MaxSeedLongTrackChi2",  m_maxSeedLongTrackChi2 = 10 );
    declareProperty("MaxSeedShortTrackChi2", m_maxSeedShortTrackChi2 = 50 );  // short tracks can be constructed from > 1 cluster per plane, so chi^2s are usually bigger
	 
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
    declareProperty("GammaOtherProngColor",     m_Color_GammaOtherProng = 0xFFA500); //-- orange
    
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

    // dedx uncertainties
    m_dedx_uncertainties.push_back("Mass_Up");
    m_dedx_uncertainties.push_back("Mass_Down");
    m_dedx_uncertainties.push_back("MEU_Up");
    m_dedx_uncertainties.push_back("MEU_Down");
    m_dedx_uncertainties.push_back("BetheBloch_Up");
    m_dedx_uncertainties.push_back("BetheBloch_Down");
    m_dedx_uncertainties.push_back("Birks");

    // Protected properties from IInteractionHypothesis.
    m_hypMeths.push_back( m_anaSignature );
    declareProperty("HypothesisMethods", m_hypMeths);
    info() << " CCProtonPi0 Hypothesis added " << endmsg;
   
}
    
//==============================================================================
// Initialize
//==============================================================================
StatusCode CCProtonPi0::initialize()
{
    debug() <<"Enter CCProtonPi0::initialize()" << endmsg;
    
    // Initialize the base class.  This will fail if you did not define m_anaSignature.
    StatusCode sc = this->MinervaAnalysisTool::initialize();
    if( sc.isFailure() ) { 
        return Error( "Failed to initialize!", sc ); 
    }
    info()<<"   initialized MinervaAnalysisTool"<<endmsg;
    
    //Seed the RNG
    m_randomGen = new TRandom3( m_randomSeed );   
    if ( m_randomGen == NULL ) return StatusCode::FAILURE;
    
    
    // Reset Counters for Functions - Debugging Purposes
    // They can be used to track events in DST and Ana
    N_tagTruth = 0;
    N_reconstructEvent = 0;
    
    m_detectableGammaE      = 0;        // MeV
    m_detectablePi0KE       = 0;        // MeV
    m_detectableProtonKE    = 120;      // MeV
    
    // Fiducial Volume
    m_fidHexApothem  = 850.0*CLHEP::mm;
    m_fidUpStreamZ   = 5990.0*CLHEP::mm;   // ~middle of module 27, plane 1
    m_fidDownStreamZ = 8340.0*CLHEP::mm;   // ~middle of module 79, plane 1
    
    // Reconstructable Volume
    m_recoHexApothem  = 1000.0*CLHEP::mm; 
    m_recoUpStreamZ   = 5810.0*CLHEP::mm;
    m_recoDownStreamZ = 8600.0*CLHEP::mm;
//     m_recoUpStreamZ   = 5750.0*CLHEP::mm;
//     m_recoDownStreamZ = 8700.0*CLHEP::mm;
    
    
    // Initializing Analysis Tools
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
     
    try {
        m_LikelihoodPIDTool = tool<IParticleTool>( "LikelihoodPIDTool" );
    }catch( GaudiException& e ) {
        error() << "Could not obtain tool: LikelihoodPIDTool!" << endmsg;
        return StatusCode::FAILURE;
    }
   
    try{ 
        m_vertexFitter = tool<IVertexFitter>("VertexFitterKalman"); 
    } catch( GaudiException& e ){
        error() << "Could not obtain VertexFitterKalman"<<endmsg;
        return StatusCode::FAILURE;
    }

    try{ 
        m_vertexEnergyStudyTool = tool<IVertexEnergyStudyTool>("VertexEnergyStudyTool"); 
    }catch( GaudiException& e ){
        error() << "Could not obtain VertexEnergyStudyTool"<< endmsg;
        return StatusCode::FAILURE;
    }

    try { 
        m_anchoredTracker = tool<IAnchoredTrackFormation>("AnchoredShortTracker"); 
    }catch( GaudiException& e ){
        error() << "Could not obtain AnchoredShortTracker"<< endmsg;
        return StatusCode::FAILURE;
    }

    try { 
        m_coneUtilsTool = tool<IConeUtilsTool>("ConeUtilsTool"); 
    } catch( GaudiException& e ){
        error() << "Could not obtain ConeUtilsTool!" << endmsg;
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

    // Blob Tools
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
    
    // Services
    try{
        m_gigaCnvSvc = svc<IGiGaGeomCnvSvc>("GiGaGeo", true);
    } catch( GaudiException& e ){
        error() <<"Could not obtain GiGaGeo"<<endmsg;
        return StatusCode::FAILURE;
    }
    
    service("GeomUtilSvc", m_GeomUtilSvc, true);
    m_idDet = m_GeomUtilSvc->getIDDet();
    m_odDet = m_GeomUtilSvc->getODDet();
    
    // declare common branches
    declareCommonPhysicsAnaBranches();
    declareMinosMuonBranches();
    declareGenieWeightBranches();
    
    //--------------------------------------------------------------------------
    // Select the branches you want in your AnaTuple
    //--------------------------------------------------------------------------
    const double SENTINEL = -9.9; 
    // ------------------------------------------------------------------------
    // Truth Branches
    // ------------------------------------------------------------------------
    declareDoubleTruthBranch("eventID",0.0);
    
    declareIntTruthBranch("N_FSParticles", -1 );
    declareIntTruthBranch("N_proton", -1 );
    declareIntTruthBranch("N_pi0", -1 );
    declareIntTruthBranch("N_gamma", -1 );
    
    declareIntTruthBranch( "vertex_module", 500);
    declareIntTruthBranch( "vertex_plane", 0);
    declareIntTruthBranch( "target_material", -1); 
    
    declareContainerDoubleTruthBranch("muon_4P", 4, SENTINEL );
    declareContainerDoubleTruthBranch("proton_4P", 4, SENTINEL );
    declareContainerDoubleTruthBranch("pi0_4P", 4, SENTINEL );
    declareContainerDoubleTruthBranch("gamma1_4P", 4, SENTINEL );
    declareContainerDoubleTruthBranch("gamma2_4P", 4, SENTINEL );
    
    declareIntTruthBranch("pi0_status", -9 );
    declareIntTruthBranch("pi0_Mother", -9 );
    declareIntTruthBranch("pi0_MotherStatus", -9 );   
    declareIntTruthBranch("pi0_GrandMother", -9 );
    declareIntTruthBranch("pi0_GrandMotherStatus", -9 );
    
    // Signal and Background
    declareBoolTruthBranch("isSignal");
    declareBoolTruthBranch("isSignal_1Pi0");
    declareBoolTruthBranch("isSignal_2Gamma");
    declareBoolTruthBranch("isFidVol");
    declareBoolTruthBranch("AnalyzeEvent");
    
    // Background Types
    declareBoolTruthBranch("isBckg_QELike");
    declareBoolTruthBranch("isBckg_SinglePiPlus");
    declareBoolTruthBranch("isBckg_SinglePiMinus");
    declareBoolTruthBranch("isBckg_MultiPion");
    declareBoolTruthBranch("isBckg_MultiPiZero");
    declareBoolTruthBranch("isBckg_Other");
    
    // Background Branching - For Each Background Type
    declareBoolTruthBranch("isBckg_withAntiMuon");
    declareBoolTruthBranch("isBckg_withMichel");
    declareBoolTruthBranch("isBckg_withPrimaryPi0");
    declareBoolTruthBranch("isBckg_withSecondaryPi0");
    
    // Truth Michel Information
    declareIntTruthBranch("N_trueMichelElectrons", -1 );
    declareDoubleTruthBranch("michelElectron_P", SENTINEL );
    declareDoubleTruthBranch("michelElectron_E", SENTINEL );
    declareContainerDoubleTruthBranch("michelMuon_endPoint", 3, SENTINEL ); 
    declareDoubleTruthBranch("michelMuon_P", SENTINEL );
    declareDoubleTruthBranch("michelMuon_length", SENTINEL );
    declareDoubleTruthBranch("michelMuon_end_dist_vtx", SENTINEL );
    declareDoubleTruthBranch("michelPion_P", SENTINEL );
    declareDoubleTruthBranch("michelPion_length", SENTINEL );
    declareDoubleTruthBranch("michelPion_begin_dist_vtx", SENTINEL );

    // ------------------------------------------------------------------------       
    // Event Branches
    // ------------------------------------------------------------------------       
    // Cut Results
    declareIntEventBranch( "Cut_Vertex_None", -1 );
    declareIntEventBranch( "Cut_Vertex_Not_Reconstructable", -1 );
    declareIntEventBranch( "Cut_Vertex_Not_Fiducial", -1 );
    declareIntEventBranch( "Cut_UnattachedProngsWithTracks", -1 );
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
    declareIntEventBranch( "Cut_ConeBlobs", -1 );
    declareIntEventBranch( "Cut_BlobsBad", -1 );
    
    // General Reco
    declareIntEventBranch( "nProngs", -1);
    declareIntEventBranch( "nTracks", -1);
    declareDoubleEventBranch( "time", -1.0 );
    declareDoubleEventBranch( "reco_eventID", -1.0 );
   
    // Vertex Information 
    declareIntEventBranch("vtx_total_count",-1);
    declareIntEventBranch("vtx_secondary_count",-1);
    declareIntEventBranch("vtx_primary_index",-1);
    declareIntEventBranch("vtx_primary_multiplicity",-1);
    declareIntEventBranch("vtx_fit_converged",-1);
    declareDoubleEventBranch("vtx_fit_chi2", -1);
    declareContainerDoubleEventBranch("fit_vtx", 3, -1);

    // Michel Prong  
    declareDoubleEventBranch( "michelProng_distance", -1.0);
    declareDoubleEventBranch( "michelProng_energy", -1.0);
    declareDoubleEventBranch( "michelProng_time_diff", -1.0);
    declareDoubleEventBranch( "michelProng_end_Z", -1.0);
    declareDoubleEventBranch( "michelProng_begin_Z", -1.0);
    
    // Visible Energy inside Detector
    declareDoubleEventBranch( "muonVisibleE", -1.0 );
    declareDoubleEventBranch( "hadronVisibleE", -1.0 );
    declareDoubleEventBranch( "totalVisibleE", -1.0 );
    declareDoubleEventBranch( "totalIDVisibleE", -1.0 );
    declareDoubleEventBranch( "totalODVisibleE", -1.0 );
    declareDoubleEventBranch( "energyUnused_afterReco", -1.0 );
    declareDoubleEventBranch( "energyUsed_afterReco", -1.0 );
      
    // PreFilterPi0()
    declareIntEventBranch("preFilter_Result", -1);
    declareDoubleEventBranch("preFilter_rejectedEnergy", -1.0);
    declareDoubleEventBranch("evis_nearvtx", -1.0);
    declareDoubleEventBranch("evis_total", -1.0);
    declareDoubleEventBranch("evis_NuclearTarget",  -1.0);
    declareDoubleEventBranch("evis_Tracker",  -1.0);
    declareDoubleEventBranch("evis_ECAL",  -1.0);
    declareDoubleEventBranch("evis_HCAL",  -1.0);
    declareDoubleEventBranch("evis_TotalExceptNuclearTarget", -1.0);
    
    // VtxBlob()
    declareContainerDoubleEventBranch("Vertex_energy_radii");
    declareDoubleEventBranch("Vertex_blob_energy", -9999 );
    declareDoubleEventBranch("Filament_Vertex_energy", -9999 );
    declareDoubleEventBranch("Sphere_Vertex_energy", -9999 );
    
    // ConeBlobs()
    declareDoubleEventBranch("RE_energy_Tracker", -9999 );
    declareDoubleEventBranch("RE_energy_ECAL", -9999 );
    declareDoubleEventBranch("RE_energy_HCAL", -9999 );
    declareIntEventBranch("anglescan_ncandx", -1);
    declareIntEventBranch("anglescan_ncand", -1);
    
    // processBlobs() -- Called from ConeBlobs()
    declareIntEventBranch("blob_ndof_1");
    declareIntEventBranch("blob_ndof_2");
    declareDoubleEventBranch("blob_fval_1");
    declareDoubleEventBranch("blob_fval_2");    

    // ODActivity() -- Called from ConeBlobs()
    declareIntEventBranch( "od_energeticTower", -9999 );
    declareDoubleEventBranch( "od_upstreamFrame", -9999 );
    declareDoubleEventBranch( "od_downstreamFrame", -9999 );
    declareDoubleEventBranch( "od_upstreamFrame_z", -9999 );
    declareDoubleEventBranch( "od_downstreamFrame_z", -9999 );
    declareDoubleEventBranch( "od_highStory", -9999 );
    declareDoubleEventBranch( "od_lowStory", -9999 );
    declareDoubleEventBranch( "od_highStory_t", -9999 );
    declareDoubleEventBranch( "od_lowStory_t", -9999 );
    declareDoubleEventBranch( "od_maxEnergy", -9999 );
    declareContainerDoubleEventBranch( "od_distanceBlobTower");
    declareContainerDoubleEventBranch( "od_towerEnergy" );
    declareContainerDoubleEventBranch( "od_towerNClusters" );
    declareContainerDoubleEventBranch( "od_towerTime" );
    declareContainerDoubleEventBranch( "od_idBlobTime" );
    declareContainerDoubleEventBranch( "od_towerTimeBlobOD" );
    declareContainerDoubleEventBranch( "od_towerTimeBlobMuon" );
 
    // Blob Data -- Filled in setBlobData()
    declareIntEventBranch("gamma1_blob_ndigits",-1);
    declareIntEventBranch("gamma1_blob_nclusters",-1);
    declareDoubleEventBranch("gamma1_blob_energy",-1);
    declareDoubleEventBranch("gamma1_blob_minsep", -1.0);
    
    declareIntEventBranch("gamma2_blob_ndigits",-1);
    declareIntEventBranch("gamma2_blob_nclusters",-1);
    declareDoubleEventBranch("gamma2_blob_energy",-1);
    declareDoubleEventBranch("gamma2_blob_minsep", -1.0);

    // Calculate_dEdX() -- Called from setBlobData()
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

    //-------------------------------------------------------------------------
    // NeutrinoInt Branches 
    //-------------------------------------------------------------------------
    // Primary Vertex
    declareIntBranch( m_hypMeths, "vtx_module", -99);
    declareIntBranch( m_hypMeths, "vtx_plane",-1);
    declareDoubleBranch( m_hypMeths, "vtx_x",0.0);
    declareDoubleBranch( m_hypMeths, "vtx_y",0.0);
    declareDoubleBranch( m_hypMeths, "vtx_z",0.0);
    
    // Event Kinematics -- Filled in setEventKinematics()
    declareDoubleBranch( m_hypMeths, "neutrino_E_1Track", -9.9);
    declareDoubleBranch( m_hypMeths, "neutrino_E_2Track", -9.9);
    declareDoubleBranch( m_hypMeths, "neutrino_E_Cal", -9.9);
    declareDoubleBranch( m_hypMeths, "QSq_1Track", -9.9);
    declareDoubleBranch( m_hypMeths, "QSq_Cal", -9.9);
    declareDoubleBranch( m_hypMeths, "WSq_1Track", -9.9);
    declareDoubleBranch( m_hypMeths, "WSq_Cal", -9.9);
    declareDoubleBranch( m_hypMeths, "W_1Track", -9.9); 
    declareDoubleBranch( m_hypMeths, "W_Cal", -9.9); 
    
    // Muon Kinematics  -- Filled in setMuonData()
    declareIntBranch( m_hypMeths, "muon_hasMinosMatchTrack", -1);
    declareIntBranch( m_hypMeths, "muon_hasMinosMatchStub", -1);
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
    
    // Proton Kinematics -- Filled in setProtonData()
    declareContainerIntBranch(m_hypMeths,    "proton_kinked", 10,    -1);
    declareContainerIntBranch(m_hypMeths,    "proton_odMatch", 10,  -1);
    declareContainerDoubleBranch(m_hypMeths, "proton_length",10,  SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_startPointX",10,  SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_startPointY",10,  SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_startPointZ", 10, SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_endPointX", 10,   SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_endPointY", 10,   SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_endPointZ", 10,   SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "protonScore",   10,  SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "pionScore",   10,  SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "protonScore_LLR",  10,   SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_chi2_ndf", 10, SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_theta",   10,  SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_thetaX",   10, SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_thetaY",   10, SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_phi",       10,SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_ekin",      10,SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_E",         10,SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_p",         10,SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_px",        10,SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_py",        10,SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_pz",        10,SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_p_calCorrection", 10, SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_p_visEnergy",      10,SENTINEL);
    declareContainerDoubleBranch(m_hypMeths, "proton_p_dEdXTool",       10,SENTINEL);
    
    // Pi0 & Gamma1,2 Kinematics -- Filled in setPi0Data()    
    declareDoubleBranch(m_hypMeths,"pi0_px",SENTINEL);
    declareDoubleBranch(m_hypMeths,"pi0_py",SENTINEL);
    declareDoubleBranch(m_hypMeths,"pi0_pz",SENTINEL);
    declareDoubleBranch(m_hypMeths,"pi0_E",SENTINEL);
    declareDoubleBranch(m_hypMeths,"pi0_invMass", SENTINEL);
    declareDoubleBranch(m_hypMeths,"pi0_theta", SENTINEL);
    declareDoubleBranch(m_hypMeths,"pi0_phi",   SENTINEL);
    declareDoubleBranch(m_hypMeths,"pi0_thetaX", SENTINEL);
    declareDoubleBranch(m_hypMeths,"pi0_thetaY",   SENTINEL);
    declareDoubleBranch(m_hypMeths,"pi0_openingAngle",  SENTINEL);
    declareDoubleBranch(m_hypMeths,"pi0_cos_openingAngle", SENTINEL);
    
    declareDoubleBranch(m_hypMeths,"gamma1_px",SENTINEL);
    declareDoubleBranch(m_hypMeths,"gamma1_py",SENTINEL);
    declareDoubleBranch(m_hypMeths,"gamma1_pz",SENTINEL);
    declareDoubleBranch(m_hypMeths,"gamma1_E",SENTINEL);
    declareDoubleBranch(m_hypMeths,"gamma1_theta",SENTINEL);
    declareDoubleBranch(m_hypMeths,"gamma1_phi",SENTINEL);
    declareDoubleBranch(m_hypMeths,"gamma1_dEdx",SENTINEL);
    declareDoubleBranch(m_hypMeths,"gamma1_time",SENTINEL);
    declareDoubleBranch(m_hypMeths,"gamma1_dist_vtx",SENTINEL);
    declareContainerDoubleBranch(m_hypMeths,"gamma1_direction",3,SENTINEL);
    declareContainerDoubleBranch(m_hypMeths,"gamma1_vertex",3,SENTINEL);
    declareDoubleBranch(m_hypMeths,"gamma1_evis_trkr", SENTINEL);
    declareDoubleBranch(m_hypMeths,"gamma1_evis_ecal", SENTINEL);
    declareDoubleBranch(m_hypMeths,"gamma1_evis_hcal", SENTINEL);
    declareDoubleBranch(m_hypMeths,"gamma1_evis_scal", SENTINEL);
       
    declareDoubleBranch(m_hypMeths,"gamma2_px",SENTINEL);
    declareDoubleBranch(m_hypMeths,"gamma2_py",SENTINEL);
    declareDoubleBranch(m_hypMeths,"gamma2_pz",SENTINEL);
    declareDoubleBranch(m_hypMeths,"gamma2_E",SENTINEL);
    declareDoubleBranch(m_hypMeths,"gamma2_theta",SENTINEL);
    declareDoubleBranch(m_hypMeths,"gamma2_phi",SENTINEL);
    declareDoubleBranch(m_hypMeths,"gamma2_dEdx",SENTINEL);
    declareDoubleBranch(m_hypMeths,"gamma2_time",SENTINEL);
    declareDoubleBranch(m_hypMeths,"gamma2_dist_vtx",SENTINEL);
    declareContainerDoubleBranch(m_hypMeths,"gamma2_direction",3,SENTINEL);
    declareContainerDoubleBranch(m_hypMeths,"gamma2_vertex",3,SENTINEL);
    declareDoubleBranch(m_hypMeths,"gamma2_evis_trkr", SENTINEL);
    declareDoubleBranch(m_hypMeths,"gamma2_evis_ecal", SENTINEL);
    declareDoubleBranch(m_hypMeths,"gamma2_evis_hcal", SENTINEL);
    declareDoubleBranch(m_hypMeths,"gamma2_evis_scal", SENTINEL);
    
    // Truth Match for Prongs
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
 
    return sc;
}
    
//==============================================================================
//
// reconstructEvent() --
//
//==============================================================================
StatusCode CCProtonPi0::reconstructEvent( Minerva::PhysicsEvent *event, Minerva::GenMinInteraction* truthEvent ) const
{
    debug()<<"Enter CCProtonPi0::reconstructEvent()" << endmsg;
    debug()<<"Version: Modified while the job test2 is Running!"<<endmsg; 
    //--------------------------------------------------------------------------
    // Count Events Enters reconstructEvent()
    //      N_reconstructEvent is used as the reco_eventID
    //--------------------------------------------------------------------------
    debug()<<"reconstructEvent() Count = "<<N_reconstructEvent<<" tagTruth() Count = "<<N_tagTruth-1<<endmsg;
    event->setDoubleData("reco_eventID", N_reconstructEvent);
    N_reconstructEvent++;
    
    //--------------------------------------------------------------------------
    // Initialize Global Variables
    //--------------------------------------------------------------------------
    // SmartRefs
    m_PrimaryVertex = NULL;
    m_MuonTrack = NULL;
    m_MuonProng = NULL;
    m_MuonParticle = NULL;
    m_Pi0Blob1 = NULL;
    m_Pi0Blob2 = NULL;
    // Vectors
    m_ProtonProngs.clear();
    m_ProtonParticles.clear();
    
    //--------------------------------------------------------------------------
    // Do NOT Analyze the Event -
    //      if the TRUE Vertex is NOT Fiducial
    //      if the AnalyzeEvent is False - TagTruth() Decides This
    //--------------------------------------------------------------------------    
    if( truthEvent ){
        info() << "This is a MC event." << endmsg;
        if ( !(truthEvent->filtertaglist()->isFilterTagTrue("isFidVol")) ){
            info() << "True Vertex is NOT Fiducial! Skipping Event! " <<endmsg;
            return StatusCode::SUCCESS; 
        }
        
        if ( !(truthEvent->filtertaglist()->isFilterTagTrue("AnalyzeEvent")) ){
            info() << "TagTruth() Marked Event as Do NOT Analyze! Skipping Event! " <<endmsg;
            return StatusCode::SUCCESS; 
        }
    }
 
    //--------------------------------------------------------------------------
    // Do NOT Analyze the Event - if Event has a Bad Object
    //--------------------------------------------------------------------------
    if( event->filtertaglist()->isFilterTagTrue( AnaFilterTags::BadObject() ) ) { 
        warning() << "Found an event flagged with a BadObject! Refusing to analyze..." << endmsg;
        return StatusCode::SUCCESS; 
    }
    
    
    //==========================================================================
    //
    // Vertex Reconstruction
    //
    //==========================================================================
    debug() << "START: Vertex Reconstruction..." << endmsg;
    
    //--------------------------------------------------------------------------
    // MAKE CUT - if NO or NULL Interaction Vertex
    //--------------------------------------------------------------------------
    if( !(event->hasInteractionVertex()) || event->interactionVertex() == NULL ){
        debug() << "The event has NO or NULL vertex!" << endmsg;
        event->setIntData("Cut_Vertex_None",1);
        if( m_store_all_events ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    }
    
    // Save Primary Vertex
    m_PrimaryVertex = event->interactionVertex();
    
    //--------------------------------------------------------------------------
    // MAKE CUT - if Interaction Vertex is NOT in Reconstructable Volume
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
    

    // Check the number of unattached prongs in the event with tracks
    int ntracks = 0;
    Minerva::ProngVect unattachedProngs = event->select<Minerva::Prong>("Used:Unused","!IsSubProng");
    for(unsigned int i = 0; i < unattachedProngs.size(); i++) {
        if( !unattachedProngs[i]->minervaTracks().empty() ) ntracks++;
    }
    
    //--------------------------------------------------------------------------
    // MAKE CUT - if Event has unattached Prongs with Tracks
    //--------------------------------------------------------------------------
    if( ntracks != 0 ) {
        debug() << " This event has unattached prongs with tracks!" << endmsg;
        event->setIntData("Cut_UnattachedProngsWithTracks",1);
        if( m_store_all_events ) return interpretFailEvent(event);
        else return StatusCode::SUCCESS;
    }else{ 
        debug() << "The number of unattached prong is = " << unattachedProngs.size() << endmsg;
    }

    //--------------------------------------------------------------------------
    // Make tracks using the anchor short tracker and refit the vertex
    //--------------------------------------------------------------------------
    bool make_primary_short_tracks = true;
    bool createdTracks = createdAnchoredShortTracks(event,m_PrimaryVertex,make_primary_short_tracks);
    if (createdTracks){
        debug()<<"Succesfully created Anchored Short Tracks"<<endmsg;
    }
    event->setIntData("nTracks", (int)m_PrimaryVertex->getOutgoingTracks().size() );
    debug() << "The vertex has " << m_PrimaryVertex->getOutgoingTracks().size()  << " outgoing primary tracks!" << endmsg;

    // set vertex fit information
    double fit_chi2 = -1;
    if( m_PrimaryVertex->hasDoubleData("fit_chi2") ) {
        fit_chi2 = m_PrimaryVertex->getDoubleData("fit_chi2"); 
        event->setDoubleData("vtx_fit_chi2",fit_chi2);
        debug() << "   The vertex fit chi2 = " << fit_chi2 << endmsg;  
    }

    int fit_converged = -1;
    if( m_PrimaryVertex->hasIntData("fit_converged") ) {
        fit_converged = m_PrimaryVertex->getIntData("fit_converged");
        event->setIntData("vtx_fit_converged",fit_converged);
        debug() << "   The vertex fit convergence = " << fit_converged << endmsg;
    }

    // re-set the vertex position for vertex fit failures for long-long tracks combination
    if( fit_converged == 0 ) {
        unsigned int longtracker = 0;
        SmartRefVector<Minerva::Track> trks = m_PrimaryVertex->getOutgoingTracks();
        for(unsigned int t = 0; t < trks.size(); t++) {
            Minerva::Track::PatRecHistory pat = trks[t]->patRecHistory();
            if( pat == Minerva::Track::LongPatRec3View || pat == Minerva::Track::LongPatRec2View ) longtracker++;
        }
        if( longtracker == trks.size() ) {
            m_vertexFitter->fit(m_PrimaryVertex);
            event->setInteractionVertex(m_PrimaryVertex);
            Minerva::ProngVect prongs = event->primaryProngs();
            for(unsigned int p = 0; p < prongs.size(); p++){
                prongs[p]->setSourceVertex(m_PrimaryVertex);
            }
        }
   }

    std::vector<double> fit_vtx;
    fit_vtx.push_back( m_PrimaryVertex->position().x() );
    fit_vtx.push_back( m_PrimaryVertex->position().y() );
    fit_vtx.push_back( m_PrimaryVertex->position().z() );

    event->setContainerDoubleData("fit_vtx",fit_vtx);

    //--------------------------------------------------------------------------
    // MAKE CUT - if Interaction Vertex is NOT in Fiducial Volume
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
         
    if( !FiducialPointTool->isFiducial( vtx_position, m_fidHexApothem, m_fidUpStreamZ, m_fidDownStreamZ ) ){
        debug() <<"Interaction Vertex is not in fiducial volume = ("<<vtx_position.x()<<","<<vtx_position.y()<<","<<vtx_position.z()<<")"<< endmsg;
        event->setIntData("Cut_Vertex_Not_Fiducial",1);
        if( m_store_all_events ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    }
    
    //--------------------------------------------------------------------------
    // VERTEX Passed all CUTS - Extend Vertex Reconstruction
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
    
    bool foundMuon = MuonUtils->findMuonProng( event, m_MuonProng, m_MuonParticle );
    debug() << "m_MuonProng and m_MuonParticle is Saved!" << endmsg;
    
    //--------------------------------------------------------------------------
    // MAKE CUT - if NO GOOD MUON (MINOS Matched)
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
    // MAKE CUT - IF Muon is NOT Plausible (MC Only)
    //--------------------------------------------------------------------------
    double mc_frac = -1.0;
    if ( m_DoPlausibilityCuts && !muonIsPlausible( m_MuonProng, mc_frac) ) {
        debug()<<"Muon is not plausible"<<endmsg;
        event->setIntData("Cut_Muon_Not_Plausible",1);
        if( m_store_all_events ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    }
    
    //--------------------------------------------------------------------------
    // MAKE CUT - if MUON Score is LOW
    //--------------------------------------------------------------------------  
    debug() << "Muon Particle Score: " << m_MuonParticle->score() << endmsg;
    if(m_MuonParticle->score() < m_minMuonScore){
        debug()<<"Muon prong does not pass score cut"<<endmsg;
        event->setIntData("Cut_Muon_Score_Low",1);
        if( m_store_all_events ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    }
    
    //--------------------------------------------------------------------------
    // MAKE CUT - if Muon has positive charge (AntiMuon)
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
    // MUON Passed All Cuts -- tag it as "PrimaryMuon"
    //--------------------------------------------------------------------------
    m_MuonProng->filtertaglist()->setOrAddFilterTag( "PrimaryMuon", true );
    m_MuonParticle->filtertaglist()->setOrAddFilterTag( "PrimaryMuon", true );
    
    // Save Muon Track info
    SmartRef<Minerva::Track> muonTrack = m_MuonProng->minervaTracks().front();
    m_MuonTrack = muonTrack;
    debug() << "m_MuonTrack is Saved!" << endmsg;
    
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
    // MAKE CUT - Vertex Michels
    //--------------------------------------------------------------------------
    debug() << "START: Vertex Michels" << endmsg;
    Minerva::Prong vtx_michel_prong;
    bool foundMichel = m_michelVtxTool->findMichel( m_PrimaryVertex, vtx_michel_prong );
    if (foundMichel) {
        debug()<<"Found a Michel Electron!"<<endmsg;
        
        saveMichelProngToNTuple(event,vtx_michel_prong);
        
        // Apply Extra Michel Cuts - If Activated!
        if(m_applyExtraMichelCuts){
            if(isMichelProngGood(vtx_michel_prong)){
                event->setIntData("Cut_Vertex_Michel_Exist",1);
                if( m_store_all_events ) return interpretFailEvent(event); 
                else return StatusCode::SUCCESS; 
            }else{
                debug()<<"Michel Prong did NOT passed Quality Cut - Keeping Event!"<<endmsg;
            }
        }else{
            event->setIntData("Cut_Vertex_Michel_Exist",1);
            if( m_store_all_events ) return interpretFailEvent(event); 
            else return StatusCode::SUCCESS; 
        }
    
    }else{
        debug()<<"There are NO Vertex Michel Electrons in the event!"<<endmsg;
    }
    
    debug() << "FINISH: Vertex Michel" << endmsg;
    
    //--------------------------------------------------------------------------
    // MAKE CUT - End Point Michels
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
                
            saveMichelProngToNTuple(event,michelProng);

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
                    
                    saveMichelProngToNTuple(event,michelProng);
                    
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
    // Pi0 Reconstruction
    //
    //==========================================================================
    
    debug()<<"START: Pi0 Reconstruction"<<endmsg;
    //--------------------------------------------------------------------------
    // MAKE CUT - if fails PreFilterPi0()
    //--------------------------------------------------------------------------
    if ( !PreFilterPi0(event) ){
        event->setIntData("Cut_PreFilter_Pi0",1);
        if( m_store_all_events ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS;  
    }
    
    // Vertex Blob Reconstruction
    VtxBlob(event);
    
    // Blob Reconstruction -- ConeBlobs()
    // MAKE CUT - If ConeBlobs Can NOT Find Two Blobs
    bool FoundTwoBlobs = ConeBlobs(event);
    if ( !FoundTwoBlobs ){
        event->setIntData("Cut_ConeBlobs",1);
        if( m_store_all_events ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS;  
    }
  
    // Set Data for Blobs found in ConeBlobs()
    setBlobData(event);

    // MAKE CUT - If Blobs are NOT Good
    bool blobsGood = AreBlobsGood();
     if ( !blobsGood ){
        event->setIntData("Cut_BlobsBad",1);
        if( m_store_all_events ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS;  
    }
  
   
    debug()<<"FINISH: Pi0 Reconstruction"<<endmsg;
    
    //==========================================================================
    //
    // Proton Reconstruction
    //      Only If there are more than 1 Prong (1 Prong = Muon)
    //
    //==========================================================================
    
    int nPrimaryProngs = primaryProngs.size();
    debug()<<"nProngs in Event = "<<nPrimaryProngs<<endmsg;
	event->setIntData("nProngs", nPrimaryProngs);
    
    if (nPrimaryProngs > 1){
        debug() << "START: Proton Reconstruction" << endmsg;
        
        //--------------------------------------------------------------------------
        // MAKE CUT - if Particle Creation Fails (Proton and Pion Hypotheses)
        //--------------------------------------------------------------------------
        
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
        // MAKE CUT - if NO Prong contains a Proton Particle
        //--------------------------------------------------------------------------
        bool foundProton = getProtonProng(primaryProngs);
        if( !foundProton ) {
            debug() << "Didn't find any contained in the tracker bit-positive prong with a proton particle!" << endmsg;
            event->setIntData("Cut_Proton_None",1);
            if( m_store_all_events ) return interpretFailEvent(event); 
            else return StatusCode::SUCCESS; 
        }
       
        debug()<<"FINISH: Proton Reconstruction"<<endmsg;
    }
    
    //==========================================================================
    //
    // End of Reconstruction!
    //
    //==========================================================================
    
    //--------------------------------------------------------------------------    
    // Get Unused and Used Energy after Reconstruction
    //--------------------------------------------------------------------------
    double totalUnusedEnergy_afterReco = getClusterEnergy(event,"Unused");
    double totalUsedEnergy_afterReco = getClusterEnergy(event,"Used");
    event->setDoubleData("energyUnused_afterReco",totalUnusedEnergy_afterReco);
    event->setDoubleData("energyUsed_afterReco",totalUsedEnergy_afterReco);
    
    //--------------------------------------------------------------------------
    // Write FS Particle Table and Event Record for Reconstructed Events
    //--------------------------------------------------------------------------
    if(truthEvent){
        bool isSignal = truthEvent->filtertaglist()->isFilterTagTrue("isSignal");
        if (m_writeFSParticle_Table){
            info()<<"FS Particle Table and Event Record for Reconstructed Event!"<<endmsg;
            writeBackgroundType(truthEvent);
            writeFSParticleTable(isSignal);
            writeEventRecord(truthEvent,isSignal);
        }
    }
     
    //--------------------------------------------------------------------------
    // Get total visible energy in PhysicsEvent
    //--------------------------------------------------------------------------
    debug()<<"Get Total Visible Energy in the Event"<<endmsg;
    SmartRefVector<IDCluster> idClusters = event->select<IDCluster>("Used:Unused","!XTalkCandidate");
    SmartRefVector<ODCluster> odClusters = event->select<ODCluster>("Used:Unused","!XTalkCandidate");  
    double idVisibleEnergy = m_extraEnergyTool->getIDEnergy(idClusters, -1.0);
    double odVisibleEnergy = m_extraEnergyTool->getODEnergy(odClusters, -1.0);
    double totalVisibleEnergy = idVisibleEnergy + odVisibleEnergy;  

    //--------------------------------------------------------------------------
    // Color Inner Detector Unused Clusters
    //--------------------------------------------------------------------------
    debug()<<"Color Unused Clusters!"<<endmsg;
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
    // Calculate dead time
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
    // event time
    //--------------------------------------------------------------  
    double event_time = event->time();

    //--------------------------------------------------------------------------
    // Finish filling event portion of ntuple 
    //--------------------------------------------------------------------------
    event->setDoubleData("time", event_time);
    
    event->setIntData( "dead", dead );
    event->setIntData( "udead", udead );
    event->setIntData( "ddead", ddead );
    event->setIntData( "tdead", tdead );
    
    event->setDoubleData("muonVisibleE", muon_visible_energy );
    event->setDoubleData("hadronVisibleE", totalVisibleEnergy - muon_visible_energy);
    event->setDoubleData("totalVisibleE",   totalVisibleEnergy );
    event->setDoubleData("totalIDVisibleE", idVisibleEnergy );
    event->setDoubleData("totalODVisibleE", odVisibleEnergy );
    
    fillCommonPhysicsAnaBranches( event );

    //--------------------------------------------------------------------------
    // Call the interpretEvent function.
    //--------------------------------------------------------------------------
    NeutrinoVect interactions;
    StatusCode interpret = this->interpretEvent( event, truthEvent, interactions );
    
    // If there were any neutrino interactions reconstructed, mark the event.
    if( interactions.size() == 1 ){
        markEvent( event );
    }
    else {
        warning()<<"interpretEvent() returns "<<interactions.size()<<" interactions!"<<endmsg;
        return StatusCode::SUCCESS; // We didn't crash.
    }
    
    // Add the interactions to the event. Use MinervaHistoTool::addInteractionHyp.
    StatusCode sc = addInteractionHyp( event, interactions );
    
    return sc;    
}
    
//==============================================================================
//
// interpretEvent()
//
//==============================================================================
StatusCode CCProtonPi0::interpretEvent( const Minerva::PhysicsEvent *event, const Minerva::GenMinInteraction *truthEvent, NeutrinoVect& interaction_hyp ) const
{
    debug() <<"Enter CCProtonPi0::interpretEvent()" << endmsg;
    
    if( truthEvent ){
        debug() << "This is a MC event." << endmsg;
    }

    if( !event ){
        warning() << "NULL Event" << endmsg;
        return StatusCode::FAILURE; // we crashed!
    }
   
    const double SENTINEL = -9.9;
    //--------------------------------------------------------------------------
    // Create interaction hypothesis
    //--------------------------------------------------------------------------
    Minerva::NeutrinoInt* nuInt = new Minerva::NeutrinoInt(m_anaSignature);
    interaction_hyp.push_back( nuInt );

    //--------------------------------------------------------------------------
    // Identify and Store - Muon and Proton Prongs
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
    // Get MCTracks vector
    //--------------------------------------------------------------
//     std::vector<MCTrack> MCTracksVector;
//     if (truthEvent) {
//         MCTracksVector = m_MCTrackTool->getMCTracks(truthEvent);
//     }
    
    // Set Muon Kinematics
    bool muonFilled = setMuonData( nuInt );
    if( !muonFilled ){ 
        error()<<"Muon NTuple Branches did not filled!"<<endmsg;
        return StatusCode::SUCCESS;
    }

    // Set Pi0 Kinematics
    bool pi0Filled = setPi0Data( nuInt );
    if( !pi0Filled ){
        error()<<"Pi0 NTupleBranches did not filled!"<<endmsg;
        return StatusCode::SUCCESS;
    }

    // Set Proton Kinematics
    if( m_ProtonProngs.size() != 0){
        bool protonFilled = setProtonData( nuInt );
        if( !protonFilled ){ 
            error()<<"Proton NTuple Branches did not filled!"<<endmsg;
            return StatusCode::SUCCESS;
        }
    }else{
        debug()<<"No Proton Particle, Setting SENTINEL to m_proton_4P"<<endmsg;
        m_proton_4P.SetPxPyPzE(SENTINEL,SENTINEL,SENTINEL,SENTINEL);
    }
    
    //--------------------------------------------------------------------------
    // Debugging: Check values
        debug()<<"m_proton_4P = ( "<<m_proton_4P.px()<<", "<<m_proton_4P.py()<<", "<<m_proton_4P.pz()<<", "<<m_proton_4P.E()<<" )"<<endmsg;
    //--------------------------------------------------------------------------
    

    //--------------------------------------------------------------------------
    // Get total visible energy in PhysicsEvent
    //--------------------------------------------------------------------------
    SmartRefVector<IDCluster> idClusters = event->select<IDCluster>("Used:Unused","!XTalkCandidate");
    SmartRefVector<ODCluster> odClusters = event->select<ODCluster>("Used:Unused","!XTalkCandidate");  
    double odVisibleEnergy = m_extraEnergyTool->getODEnergy(odClusters, -1.0);
    double idVisibleEnergy = m_extraEnergyTool->getIDEnergy(idClusters, -1.0);
    double muon_visible_energy = m_MuonProng->minervaVisibleEnergySum();
    double totalVisibleEnergy = odVisibleEnergy + idVisibleEnergy;  
    double hadronVisibleEnergy = totalVisibleEnergy - muon_visible_energy;
   
    
    //--------------------------------------------------------------------------  
    // Calculate and Set Event Kinematics
    //--------------------------------------------------------------------------
    setEventKinematics(nuInt, hadronVisibleEnergy);
    
    //--------------------------------------------------------------------------
    // Calculate and Set Vertex Parameters
    //--------------------------------------------------------------------------
    setVertexData(nuInt, event);
    
    //--------------------------------------------------------------------------
    // Interaction Parameters
    //--------------------------------------------------------------------------
    //nuInt->setNeutrinoHelicity( getHelicity( mu_charge ) );
    nuInt->setInteractionCurrent( Minerva::NeutrinoInt::ChargedCurrent );
    nuInt->setInteractionType( Minerva::NeutrinoInt::UnknownInt );
    
    //--------------------------------------------------------------------------
    // Truth Matching for Muon and Proton Prongs
    //--------------------------------------------------------------------------
    if( haveNeutrinoMC() && m_DoTruthMatch ) {      
        // Muon Truth Matching
        Minerva::ProngVect muonProngs;
        muonProngs.push_back( muonProng );
        setTrackProngTruth(nuInt,muonProngs);
        debug() << "       set the muon prong track truth information" << endmsg;
        
        // Proton Turth Matching
        if (protonProngs.size() == 0 ){
            debug() << "       No Proton Candidate! No Truth Matching" << endmsg;
        }else{
            setTrackProngTruth(nuInt,protonProngs);
            debug() << "       set the proton prong track truth information" << endmsg;
        }
    }
    
    return StatusCode::SUCCESS;
}
   

//==============================================================================
// tagTruth()
//==============================================================================
StatusCode CCProtonPi0::tagTruth( Minerva::GenMinInteraction* truthEvent ) const 
{
    debug() << "Enter: CCProtonPi0::tagTruth()" << endmsg;
    
    //--------------------------------------------------------------------------
    // Count Events Enters tagTruth()
    //      N_tagTruth is used as the eventID to locate event in DST
    //      Required for Arachne Scannings
    //--------------------------------------------------------------------------
    debug()<<"tagTruth() Count = "<<N_tagTruth<<endmsg;
    truthEvent->setDoubleData("eventID", N_tagTruth);
    N_tagTruth++;

    if (!truthEvent) {
        warning()<<"Passed a null truthEvent to tagTruth()!"<<endmsg;
        return StatusCode::SUCCESS;
    }
    
    //--------------------------------------------------------------------------
    // is Truth Vertex inside  Fiducial volume
    //--------------------------------------------------------------------------
    bool isFidVol = FiducialPointTool->isFiducial(truthEvent, m_fidHexApothem, m_fidUpStreamZ, m_fidDownStreamZ );
    truthEvent->filtertaglist()->setOrAddFilterTag( "isFidVol", isFidVol );
    debug()<<"Finished determining true_isFidVol: "<<isFidVol<<endmsg;
    
    // Return immediately - if Event is not inside Fid Volume
    if (!isFidVol){ 
        debug() <<"Exit CCProtonPi0::tagTruth() - True vtx is not in Fiducial Volume!" << endmsg;
        return StatusCode::SUCCESS;
    }
    
    //--------------------------------------------------------------------------
    // Fill GENIE Weight Branches
    //--------------------------------------------------------------------------
    StatusCode sc = fillGenieWeightBranches( truthEvent );
    if (sc.isFailure() ) {
        warning()<<"Genie weight branch filling failed!"<<endmsg;
        return sc;
    }
    
    //--------------------------------------------------------------------------
    // Tag Event as Signal or Background
    //--------------------------------------------------------------------------
    bool isSignal = tagSignal(truthEvent);
    if (isSignal){
        setPi0GenieRecord(truthEvent);  
        SetSignalKinematics(truthEvent);
    }else{
        tagBackground(truthEvent);
    }
    
    //--------------------------------------------------------------------------
    // Create Final State Particle Table
    //--------------------------------------------------------------------------
    if (m_writeFSParticle_Table){ 
        writeFSParticleTable(isSignal);
        writeEventRecord(truthEvent,isSignal);
    }
   
    //------------------------------------------------------------
    // Set Target Material 
    //------------------------------------------------------------
    setTargetMaterial(truthEvent);
    
    //--------------------------------------------------------------
    // Set Vertex module and plane
    //--------------------------------------------------------------   
    Gaudi::LorentzVector t_vtx = truthEvent->Vtx();    
    int vertex_module;
    int vertex_plane;
    
    debug()<<"Calling getNearestPlane, t_vtx is "<<t_vtx.z()<<endmsg;   
    getNearestPlane(t_vtx.z(), vertex_module, vertex_plane);  

    debug()<<"vertex_module = "<<vertex_module<<" vertex_plane = "<<vertex_plane<<endmsg;
    
    truthEvent->setIntData("vertex_module", vertex_module);
    truthEvent->setIntData("vertex_plane", vertex_plane);
    
    //--------------------------------------------------------------------------
    // Decide whether to Analyze the Event or NOT
    //    Set AnalyzeEvent true to Analyze All Events
    //--------------------------------------------------------------------------
    truthEvent->filtertaglist()->setOrAddFilterTag( "AnalyzeEvent", true );
//     if ( truthEvent->filtertaglist()->isFilterTagTrue("isBckg_withMichel") ){
//         debug()<<"Event with True Michel - Will Analyze the Event..."<<endmsg;
//         truthEvent->filtertaglist()->setOrAddFilterTag( "AnalyzeEvent", true );
//     }else{
//         truthEvent->filtertaglist()->setOrAddFilterTag( "AnalyzeEvent", false );    
//     }
    
    debug()<<"tagTruth Reached End"<<endmsg;
    
    return StatusCode::SUCCESS;
}

//==============================================================================
//  truthIsPlausible()
//      truthIsPlausible() is now a pure virtual function of MinervaAnalysisTool, so you must implement it.
//      It is called automatically by PhysicsEventAnalysisAlg AFTER reconstructEvent() and interpretEvent() are run.
//      If it returns false, the event is not included in the analysis ntuple
//==============================================================================
bool CCProtonPi0::truthIsPlausible( const Minerva::PhysicsEvent*  ) const
{
    // Inside reconstructEvent, I am using muonIsPlausible()
    // Therefore, truthIsPlausible() returns Always TRUE 
    return true;
}

//==============================================================================
//  Finalize
//==============================================================================
StatusCode CCProtonPi0::finalize()
{
    debug() <<"Enter: CCProtonPi0::finalize()" << endmsg;
    
    // finalize the base class.
    StatusCode sc = this->MinervaAnalysisTool::finalize();
    if( sc.isFailure() ) return Error( "Failed to finalize!", sc );
    
    return sc;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
//
//  Private Functions
//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

//------------------------------------------------------------------------------
// tagSignal() 
//      Tags Event as Signal or NOT
//      Signal: CC-Neutrino Interaction inside Fiducial Volume, 
//              FS Particles: Single Pi0, No Other (Proton and Neutrons are OK)
//      returns TRUE if it is a Signal Event    
//------------------------------------------------------------------------------
bool CCProtonPi0::tagSignal(Minerva::GenMinInteraction* truthEvent) const
{
    debug() << "Enter CCProtonPi0::tagSignal()" << endmsg;
    
    bool isSignal = false;
        
    //--------------------------------------------------------------------------
    // get TG4Trajectories
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
    // Count the Number of FS Particles
    //--------------------------------------------------------------------------
    int particle_PDG;
    int N_FSParticles = 0;
    int N_proton    = 0;
    int N_neutron   = 0;
    int N_gamma     = 0;
    int N_pi0       = 0;
    int N_muon      = 0;
    int N_other     = 0;
    
    // Get Number of Final State Particles
    N_FSParticles = pri_trajectories.size();
    
    for (it_mcpart = pri_trajectories.begin(); it_mcpart != pri_trajectories.end(); ++it_mcpart) {
        particle_PDG =  (*it_mcpart)->GetPDGCode();
        if (  particle_PDG == PDG::proton ) N_proton++;
        else if ( particle_PDG == PDG::neutron ) N_neutron++;
        else if ( particle_PDG == PDG::pi0 ) N_pi0++;
        else if ( particle_PDG == PDG::gamma ) N_gamma++;
        else if ( particle_PDG == PDG::muon ) N_muon++;
        else N_other++;
    }
    
    debug()<<"N(Final State Particles) = "<<N_FSParticles<<endmsg;
    debug()<<"N(Proton) = "<<N_proton<<" N(Neutron) = "<<N_neutron<<" N(Pi0) = "<<N_pi0<<endmsg;
    debug()<<"N(Gamma) = "<<N_gamma<<" N(Other) = "<<N_other<<endmsg;
    
    truthEvent->setIntData("N_FSParticles",  N_FSParticles);
    truthEvent->setIntData("N_proton",  N_proton);
    truthEvent->setIntData("N_pi0",  N_pi0);
    truthEvent->setIntData("N_gamma",  N_gamma);
    
    //--------------------------------------------------------------------------
    // Find Signal --  CC Neutrino Interaction
    //                 FS Particles: muon, pi0(2xGamma), X (No Meson)
    //--------------------------------------------------------------------------
    bool is1Pi0 = false;
    
    int t_current = truthEvent->current();
    int t_neutrinoPDG = truthEvent->incoming();
    
    // Signal Must be CC Neutrino Interaction in Fiducial Volume
    if ( t_current == 1 && t_neutrinoPDG == PDG::nu_mu ){
        is1Pi0 = isSinglePi0(truthEvent,N_pi0,N_gamma);
        
        if( is1Pi0 && N_other == 0){ 
            isSignal = true;
            debug()<<"Found a Signal Event!"<<endmsg;
        }
    }
    
    truthEvent->filtertaglist()->setOrAddFilterTag( "isSignal", isSignal );
    
    return isSignal;
}

//------------------------------------------------------------------------------
// tagBackground() 
//      Tags Background Event
//      Signal: CC-Neutrino Interaction inside Fiducial Volume, 
//              FS Particles: Single Pi0 + X (No Meson)
//      returns TRUE if it is a Signal Event    
//------------------------------------------------------------------------------
void CCProtonPi0::tagBackground(Minerva::GenMinInteraction* truthEvent) const
{
    debug() << "Enter CCProtonPi0::tagBackground()" << endmsg;
    
    Minerva::TG4Trajectories::iterator it_mcpart;
    Minerva::TG4Trajectories* trajects = NULL;
    if( exist<Minerva::TG4Trajectories>( Minerva::TG4TrajectoryLocation::Default ) ) {
        trajects = get<Minerva::TG4Trajectories>( Minerva::TG4TrajectoryLocation::Default );
    } else {
        warning() << "Could not get TG4Trajectories from " << Minerva::TG4TrajectoryLocation::Default << endmsg;
        return;
    }
    
    // just to make sure
    if( !trajects ) return;
    
    int nPiPlus = 0;
    int nPiMinus = 0;
    int nPrimaryPiZero = 0;
    int nSecondaryPiZero = 0;
    int nAntiMuon = 0;
    int particle_PDG = -1;
    int particle_ID = -1;
    int mother_ID = -1;
    
    bool isQELike = false;
    bool isSinglePiPlus = false;
    bool isSinglePiMinus = false;
    bool isMultiPion = false;
    bool isMultiPiZero = false;
    bool isOther = false;
    
    bool hasAntiMuon = false;
    bool hasMichel = false;
    bool hasPrimaryPi0 = false;
    bool hasSecondaryPi0 = false;
    
    std::vector<int> primary_traj_ID;
    std::vector<int> primary_Pi0_ID;
    
    // Identify Primary Trajectories
    for (it_mcpart = trajects->begin(); it_mcpart != trajects->end(); ++it_mcpart) {
        
        if ( (*it_mcpart)->GetProcessName() == "Primary"){
            
            particle_PDG = (*it_mcpart)->GetPDGCode();
            particle_ID = (*it_mcpart)->GetTrackId();
            
            // Save Primary Trajectory IDs
            primary_traj_ID.push_back(particle_ID);
            
            // Count Particles
            if (particle_PDG == -(PDG::muon)) nAntiMuon++;
            if (particle_PDG == PDG::pi) nPiPlus++;
            if (particle_PDG == -(PDG::pi)) nPiMinus++;
            if (particle_PDG == PDG::pi0){ 
                nPrimaryPiZero++;
                primary_Pi0_ID.push_back(particle_ID);
            } 
        }   
    }
    
    // -------------------------------------------------------------------------
    // Background Types 
    // -------------------------------------------------------------------------
    if( isBackgroundQELike(truthEvent) ) isQELike = true;    // Check for QE-Like Events
    else if (nPiPlus == 1 && nPiMinus == 0 ) isSinglePiPlus = true;                
    else if (nPiPlus == 0 && nPiMinus == 1 ) isSinglePiMinus = true;
    else if ((nPiPlus + nPiMinus) > 1 ) isMultiPion = true; //MultiPion is Multiple Charged Pion
    // If there is NO Charged Pion out of nucleus, check for Multiple Pi0
    else if ( ((nPiPlus + nPiMinus) == 0) && nPrimaryPiZero > 1 )isMultiPiZero = true;
    else isOther = true;
    
    
    // -------------------------------------------------------------------------
    // Check for Background Branching: hasMichel, hasPrimaryPi0, hasSecondaryPi0
    // -------------------------------------------------------------------------
 
    if(nAntiMuon > 0) hasAntiMuon = true;
    hasMichel = checkMichel(truthEvent);
    
    // Loop over Secondary Trajectories to Check Secondary Pi0
    for (it_mcpart = trajects->begin(); it_mcpart != trajects->end(); ++it_mcpart) {
        
        // Skip Primary Trajectories
        if ( (*it_mcpart)->GetProcessName() == "Primary") continue;
        
        particle_PDG =  (*it_mcpart)->GetPDGCode();
        mother_ID = (*it_mcpart)->GetParentId();
        
        // Check Trajectory Mother is NOT a Primary Pi0 - Scattering
        if ( particle_PDG == PDG::pi0 && !isMotherPrimary(primary_Pi0_ID, mother_ID) ){
            nSecondaryPiZero++;
        }
    }
    
    if (nPrimaryPiZero > 0) hasPrimaryPi0 = true;
    if (nSecondaryPiZero > 0) hasSecondaryPi0 = true;

    
    debug()<<"Interaction Type = "<<truthEvent->interactionType()<<endmsg;
    debug()<<"Background Type!"<<endmsg;
    debug()<<"  QuasiElastic Like = "<<isQELike<<endmsg;
    debug()<<"  SinglePiPlus = "<<isSinglePiPlus<<endmsg;
    debug()<<"  SinglePiMinus = "<<isSinglePiMinus<<endmsg;
    debug()<<"  MultiPion = "<<isMultiPion<<endmsg;
    debug()<<"  MultiPiZero = "<<isMultiPiZero<<endmsg;
    debug()<<"  Other = "<<isOther<<endmsg;
    debug()<<"Background Branching!"<<endmsg;
    debug()<<"  hasAntiMuon = "<<hasAntiMuon<<endmsg;
    debug()<<"  hasMichel = "<<hasMichel<<endmsg;
    debug()<<"  hasPrimaryPi0 = "<<hasPrimaryPi0<<endmsg;
    debug()<<"  hasSecondaryPi0 = "<<hasSecondaryPi0<<endmsg;
 
    truthEvent->filtertaglist()->setOrAddFilterTag( "isBckg_QELike", isQELike );
    truthEvent->filtertaglist()->setOrAddFilterTag( "isBckg_SinglePiPlus", isSinglePiPlus );
    truthEvent->filtertaglist()->setOrAddFilterTag( "isBckg_SinglePiMinus", isSinglePiMinus );
    truthEvent->filtertaglist()->setOrAddFilterTag( "isBckg_MultiPion", isMultiPion );
    truthEvent->filtertaglist()->setOrAddFilterTag( "isBckg_MultiPiZero", isMultiPiZero );
    truthEvent->filtertaglist()->setOrAddFilterTag( "isBckg_Other", isOther );
    
    truthEvent->filtertaglist()->setOrAddFilterTag( "isBckg_withAntiMuon", hasAntiMuon );
    truthEvent->filtertaglist()->setOrAddFilterTag( "isBckg_withMichel", hasMichel );
    truthEvent->filtertaglist()->setOrAddFilterTag( "isBckg_withPrimaryPi0", hasPrimaryPi0 );
    truthEvent->filtertaglist()->setOrAddFilterTag( "isBckg_withSecondaryPi0", hasSecondaryPi0 );
    
}



//------------------------------------------------------------------------------
// isBackgoundQELike() 
//      Check Primary Trajectories for a Background Event to decide whether it is a QE like or not
//      Returns TRUE if event has only protons and neutrons outside of Nucleus (muon is OK)
//------------------------------------------------------------------------------
bool CCProtonPi0::isBackgroundQELike(Minerva::GenMinInteraction* truthEvent) const
{
    debug() << "Enter CCProtonPi0::isBackgoundQELike()" << endmsg;
    
    const SmartRefVector<Minerva::TG4Trajectory> pri_trajectories = truthEvent->trajectories();
    SmartRefVector<Minerva::TG4Trajectory>::const_iterator it_mcpart;
    int N_other = 0;
    int particle_PDG;
    
    for (it_mcpart = pri_trajectories.begin(); it_mcpart != pri_trajectories.end(); ++it_mcpart){
        particle_PDG =  (*it_mcpart)->GetPDGCode();
        if (  particle_PDG == PDG::proton ) continue;
        else if ( particle_PDG == PDG::neutron ) continue;
        else if ( particle_PDG == PDG::muon ) continue;
        else if ( particle_PDG == -(PDG::muon) ) continue;
        else N_other++;
    }
    
    // If there are Other Particles than proton, neutron, muon, or antimuon
    if (N_other > 0) return false;
    else return true;
    
}


//------------------------------------------------------------------------------
// checkPionAbsorption() 
//      Checks for Pion Absorption inside Nucleus for Events Does not have a pion outside
//      Returns TRUE if there is a pion inside (pi+,pi-,pi0)
//------------------------------------------------------------------------------
bool CCProtonPi0::checkPionAbsorption(Minerva::GenMinInteraction* truthEvent) const
{
    debug() <<"Enter: CCProtonPi0::checkPionAbsorption()" << endmsg;
    
    const Minerva::GenMinEventRecord* eventRecord = truthEvent->eventRecord();
    const std::vector<int> eRecPartID = eventRecord->eRecPartID();

    for(unsigned int g_part = 0; g_part < eRecPartID.size(); g_part++ ){
        if (eRecPartID[g_part] == PDG::pi) return true;
        if (eRecPartID[g_part] == -(PDG::pi)) return true;
        if (eRecPartID[g_part] == PDG::pi0) return true;
    }
    
    return false;
    
}
//------------------------------------------------------------------------------
// checkMichel() 
//      Check MC Trajectories for a PiPlus to MuPlus Decay
//      Returns TRUE if such a decay exists and fills Michel Truth Branches
//------------------------------------------------------------------------------
bool CCProtonPi0::checkMichel(Minerva::GenMinInteraction* truthEvent) const
{
    debug() <<"Enter: CCProtonPi0::checkMichel()" << endmsg;
    
    bool hasMichel = false;
    int N_trueMichels = 0;
    
    Minerva::TG4Trajectories::iterator it_mcpart;
    Minerva::TG4Trajectories* trajects = NULL;
    if( exist<Minerva::TG4Trajectories>( Minerva::TG4TrajectoryLocation::Default ) ) {
        trajects = get<Minerva::TG4Trajectories>( Minerva::TG4TrajectoryLocation::Default );
    } else {
        warning() << "Could not get TG4Trajectories from " << Minerva::TG4TrajectoryLocation::Default << endmsg;
        return false;
    }
    
    // just to make sure
    if( !trajects ) return false;
    
    // Save Vertex Information
    Gaudi::LorentzVector vtx = truthEvent->Vtx();
    
    int particle_PDG;
    int particle_ID;
    int mother_ID;
    std::vector<int> PiPlus_traj_ID;
    std::vector<double> PiPlus_traj_length;
    std::vector<double> PiPlus_dist_vtx;
    std::vector<double> PiPlus_momentum;
    
    // Save all Pi+ Info
    for (it_mcpart = trajects->begin(); it_mcpart != trajects->end(); ++it_mcpart) {
        
        particle_PDG =  (*it_mcpart)->GetPDGCode();
        particle_ID = (*it_mcpart)->GetTrackId();
        
        // Save Pi+ Information
        if (particle_PDG == PDG::pi ){
            PiPlus_traj_ID.push_back(particle_ID);
            Gaudi::LorentzVector piplus_momentum = (*it_mcpart)->GetInitialMomentum();
            Gaudi::LorentzVector piplus_start = (*it_mcpart)->GetInitialPosition();
            Gaudi::LorentzVector piplus_end = (*it_mcpart)->GetFinalPosition();
            
            double piplus_length = calcDistance(piplus_start.X(), piplus_start.Y(), piplus_start.Z(),
                                                piplus_end.X(), piplus_end.Y(), piplus_end.Z());
            
            // PiPlus Initial Point Distance to Vertex
            double piplus_dist_vtx = calcDistance(  vtx.X(), vtx.Y(), vtx.Z(),
                                                    piplus_start.X(), piplus_start.Y(), piplus_start.Z());                                                
                                                
            PiPlus_traj_length.push_back(piplus_length);
            PiPlus_momentum.push_back(piplus_momentum.P());
            PiPlus_dist_vtx.push_back(piplus_dist_vtx);
        }
    }
    
    // Loop again to locate Mu-Plus as a daughter of Pi-Plus
    for (it_mcpart = trajects->begin(); it_mcpart != trajects->end(); ++it_mcpart) {
        
        particle_PDG =  (*it_mcpart)->GetPDGCode();
        mother_ID = (*it_mcpart)->GetParentId();
        
        
        // Check particle is a Mu-Plus
        if (particle_PDG == -(PDG::muon) ){
            
            // Check Mu-Plus Mother is one of the Pi-Plus
            int michelPion_ind = getMichelPion(PiPlus_traj_ID, mother_ID);
            if ( michelPion_ind != -1 ){
                
                Gaudi::LorentzVector michelMuon_end = (*it_mcpart)->GetFinalPosition();
                
                // Skip if the michel Muon end point is NOT in Reconstructable Volume
                if ( !FiducialPointTool->isFiducial(michelMuon_end.X(),michelMuon_end.Y(),michelMuon_end.Z(), m_recoHexApothem, m_recoUpStreamZ, m_recoDownStreamZ ) ) {
                    debug() <<"Michel Muon is Not in Reconstructable Volume = ("<<michelMuon_end.X()<<","<<michelMuon_end.Y()<<","<<michelMuon_end.Z()<<")"<< endmsg;
                    continue;
                }
                
                hasMichel = true;
                
                // Save Michel Information Only for the 1st Michel
                N_trueMichels++;
                if (N_trueMichels > 1) continue;
                
                // -------------------------------------------------------------
                // Save Michel Electron Information
                // -------------------------------------------------------------
                int muon_ID = (*it_mcpart)->GetTrackId();
                saveMichelElectron(truthEvent, muon_ID);
                
                // -------------------------------------------------------------
                // Get Michel Muon Information
                // -------------------------------------------------------------
                std::vector<double> michelMuonPoint;
                Gaudi::LorentzVector michelMuon_momentum = (*it_mcpart)->GetInitialMomentum();
                Gaudi::LorentzVector michelMuon_start = (*it_mcpart)->GetInitialPosition();

                // Anti-Muon End Point
                michelMuonPoint.push_back(michelMuon_end.X());
                michelMuonPoint.push_back(michelMuon_end.Y());
                michelMuonPoint.push_back(michelMuon_end.Z());
                
                // Anti-Muon Length
                double michel_length = calcDistance(michelMuon_start.X(), michelMuon_start.Y(), michelMuon_start.Z(),
                                                    michelMuon_end.X(), michelMuon_end.Y(), michelMuon_end.Z());
                                                    
                // Distance AntiMuon End point to Vertex
                double vtx_michel_end = calcDistance(   vtx.X(), vtx.Y(), vtx.Z(),
                                                        michelMuon_end.X(), michelMuon_end.Y(), michelMuon_end.Z());
                                                        
                debug()<<"michelMuon_momentum = "<<michelMuon_momentum.P()<<endmsg;
                debug()<<"michelMuon_length = "<<michel_length<<" michelMuon_end_dist_vtx = "<<vtx_michel_end<<endmsg;

                // -------------------------------------------------------------
                // Get Michel Pion Information - Use michelPion_ind
                // -------------------------------------------------------------
                debug()<<"michelPion_momentum = "<<PiPlus_momentum[michelPion_ind]<<endmsg;
                debug()<<"michelPion_length = "<<PiPlus_traj_length[michelPion_ind]<<endmsg;
                debug()<<"michelPion_begin_dist_vtx = "<<PiPlus_dist_vtx[michelPion_ind]<<endmsg;
                
                // Fill NTuple
                truthEvent->setContainerDoubleData("michelMuon_endPoint",michelMuonPoint);
                truthEvent->setDoubleData("michelMuon_P", michelMuon_momentum.P());
                truthEvent->setDoubleData("michelMuon_length", michel_length);
                truthEvent->setDoubleData("michelMuon_end_dist_vtx", vtx_michel_end);
                truthEvent->setDoubleData("michelPion_P", PiPlus_momentum[michelPion_ind]);
                truthEvent->setDoubleData("michelPion_length", PiPlus_traj_length[michelPion_ind]);
                truthEvent->setDoubleData("michelPion_begin_dist_vtx", PiPlus_dist_vtx[michelPion_ind]);
                
            }
        } 
    }
    
    debug()<<"N_trueMichels = "<<N_trueMichels<<endmsg;
    truthEvent->setIntData("N_trueMichelElectrons", N_trueMichels);
    return hasMichel;
}

//------------------------------------------------------------------------------
// writeBackgroundType() 
//      Prints Background Type with its Branching Information
//          Uses info()
//------------------------------------------------------------------------------
void CCProtonPi0::writeBackgroundType(Minerva::GenMinInteraction* truthEvent) const
{
    // RETURN if it is a Signal Event
    if(truthEvent->filtertaglist()->isFilterTagTrue("isSignal")) return;
     
    // Print Background Type
    if (truthEvent->filtertaglist()->isFilterTagTrue("isBckg_QELike")) info() <<"Background: QE Like"<<endmsg;
    else if (truthEvent->filtertaglist()->isFilterTagTrue("isBckg_AntiMuon")) info() <<"Background: Anti-Muon"<<endmsg;
    else if (truthEvent->filtertaglist()->isFilterTagTrue("isBckg_SinglePiPlus")) info() <<"Background: Single PiPlus"<<endmsg;
    else if (truthEvent->filtertaglist()->isFilterTagTrue("isBckg_SinglePiMinus")) info() <<"Background: Single PiMinus"<<endmsg;
    else if (truthEvent->filtertaglist()->isFilterTagTrue("isBckg_MultiPion")) info() <<"Background: Multiple Charged Pion"<<endmsg;
    else if (truthEvent->filtertaglist()->isFilterTagTrue("isBckg_MultiPiZero")) info() <<"Background: Multiple PiZero"<<endmsg;
    else if (truthEvent->filtertaglist()->isFilterTagTrue("isBckg_PionAbsorption")) info() <<"Background: Pion Absorption"<<endmsg;
    else if (truthEvent->filtertaglist()->isFilterTagTrue("isBckg_Other")) info() <<"Background: Other"<<endmsg;
    else warning()<< "No Background Type!"<<endmsg;
    
    // Print Branching Information
    bool hasMichel = truthEvent->filtertaglist()->isFilterTagTrue("isBckg_withMichel");
    bool hasPi0 = truthEvent->filtertaglist()->isFilterTagTrue("isBckg_withPi0");
    bool hasGamma = truthEvent->filtertaglist()->isFilterTagTrue("isBckg_withGamma");
    
    info()<<"Background Branching: "<<endmsg;
    info()<<"   Michel = "<<hasMichel<<" hasPi0 = "<<hasPi0<<" hasGamma = "<<hasGamma<<endmsg;
}

//------------------------------------------------------------------------------
// setTargetMaterial() Stores Target Material using Truth information
//------------------------------------------------------------------------------
void CCProtonPi0::setTargetMaterial(Minerva::GenMinInteraction* truthEvent) const
{
    debug() << "Enter CCProtonPi0::setTargetMaterial()" << endmsg;
    
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
    
    truthEvent->setIntData("target_material", t_targetMaterial);
}

//------------------------------------------------------------------------------
// isMotherPrimary() returns True if the particle mother is a Primary Particle
//------------------------------------------------------------------------------
bool CCProtonPi0::isMotherPrimary(std::vector<int>& motherList, int mother ) const
{
    for( unsigned int i = 0 ; i < motherList.size(); i++ ) {
        if( mother == motherList[i]) return true;
    }
    
    return false;
}

//------------------------------------------------------------------------------
// getMichelPion() returns the indice of Michel Pion or -1
//------------------------------------------------------------------------------
int CCProtonPi0::getMichelPion(std::vector<int>& piList, int ID ) const
{
    for( unsigned int i = 0 ; i < piList.size(); i++ ) {
        if( ID == piList[i]) return i;
    }
    
    return -1;
}

//------------------------------------------------------------------------------
// saveMichelElectron() - Saves Michel Electron Information
//------------------------------------------------------------------------------
void CCProtonPi0::saveMichelElectron(Minerva::GenMinInteraction* truthEvent, int muon_ID) const
{
    debug() <<"Enter: CCProtonPi0::saveMichelElectron()" << endmsg;
    
    Minerva::TG4Trajectories::iterator it_mcpart;
    Minerva::TG4Trajectories* trajects = NULL;
    if( exist<Minerva::TG4Trajectories>( Minerva::TG4TrajectoryLocation::Default ) ) {
        trajects = get<Minerva::TG4Trajectories>( Minerva::TG4TrajectoryLocation::Default );
    } else {
        warning() << "Could not get TG4Trajectories from " << Minerva::TG4TrajectoryLocation::Default << endmsg;
        return;
    }
    
    // just to make sure
    if( !trajects ) return;
    
    // Loop Over all Trajectories to Find Michel Electrons
    for (it_mcpart = trajects->begin(); it_mcpart != trajects->end(); ++it_mcpart) {
        
        int particle_PDG =  (*it_mcpart)->GetPDGCode();
        int mother_ID = (*it_mcpart)->GetParentId();
        
        if (particle_PDG == -(PDG::electron) && mother_ID == muon_ID ){
            
            Gaudi::LorentzVector michelElectron_4P = (*it_mcpart)->GetInitialMomentum();
            truthEvent->setDoubleData("michelElectron_P", michelElectron_4P.P());
            truthEvent->setDoubleData("michelElectron_E", michelElectron_4P.E());
            
            debug()<<"Found Michel Electron with PDG = "<<particle_PDG<<" E = "<<michelElectron_4P.E()<<endmsg;
            
            break;
        }
    }
}

void CCProtonPi0::SetSignalKinematics(Minerva::GenMinInteraction* truthEvent) const
{
    debug()<<"Enter: CCProtonPi0::SetSignalKinematics"<<endmsg;
   
    const double SENTINEL = -9.9;

    Minerva::TG4Trajectories* trajects = NULL;
    Minerva::TG4Trajectories::iterator it_mcpart;
    if( exist<Minerva::TG4Trajectories>( Minerva::TG4TrajectoryLocation::Default ) ) {
        trajects = get<Minerva::TG4Trajectories>( Minerva::TG4TrajectoryLocation::Default );
    } else {
        warning() << "Could not get TG4Trajectories from " << Minerva::TG4TrajectoryLocation::Default << endmsg;
        return;
    }

    // just to make sure
    if( !trajects ) return;
    
    std::vector<double> muon_4P(4,SENTINEL);
    std::vector<double> pi0_4P(4,SENTINEL);
    std::vector<double> proton_4P(4,SENTINEL);
    std::vector<double> gamma1_4P(4,SENTINEL);
    std::vector<double> gamma2_4P(4,SENTINEL);
    
    Gaudi::LorentzVector temp_4P;
    int particle_PDG;
    int particle_ID;
    int mother_ID;
    int pi0_ID = -1;
    int max_protonE = 0.0;
    bool isFirstGamma = true;
    bool isGammaPrimary = false;
    
    // Identify and Save Primary Trajectory Kinematics
    for (it_mcpart = trajects->begin(); it_mcpart != trajects->end(); ++it_mcpart) {
        
        if ( (*it_mcpart)->GetProcessName() == "Primary"){
            
            temp_4P = (*it_mcpart)->GetInitialMomentum();
            particle_PDG = (*it_mcpart)->GetPDGCode();
            particle_ID = (*it_mcpart)->GetTrackId();
                        
            if (particle_PDG == PDG::muon){
                muon_4P[0] = temp_4P.px();  
                muon_4P[1] = temp_4P.py(); 
                muon_4P[2] = temp_4P.pz(); 
                muon_4P[3] = temp_4P.E(); 
            }else if (particle_PDG == PDG::pi0){
                pi0_ID = particle_ID;
                pi0_4P[0] = temp_4P.px();  
                pi0_4P[1] = temp_4P.py(); 
                pi0_4P[2] = temp_4P.pz(); 
                pi0_4P[3] = temp_4P.E();  
            }else if ((particle_PDG == PDG::proton) && (max_protonE < temp_4P.E()) ){
                proton_4P[0] = temp_4P.px();  
                proton_4P[1] = temp_4P.py(); 
                proton_4P[2] = temp_4P.pz(); 
                proton_4P[3] = temp_4P.E();  
            }else if (particle_PDG == PDG::gamma){
                isGammaPrimary = true;
                if (isFirstGamma){
                    gamma1_4P[0] = temp_4P.px();  
                    gamma1_4P[1] = temp_4P.py(); 
                    gamma1_4P[2] = temp_4P.pz(); 
                    gamma1_4P[3] = temp_4P.E(); 
                    isFirstGamma = false;
                }else{
                    gamma2_4P[0] = temp_4P.px();  
                    gamma2_4P[1] = temp_4P.py(); 
                    gamma2_4P[2] = temp_4P.pz(); 
                    gamma2_4P[3] = temp_4P.E(); 
                }
            }
        }   
    }
    
    
    // If Gamma's are not Primary -- Save Pi0 Daughters
    if (! isGammaPrimary ){
        isFirstGamma = true;
        
        for (it_mcpart = trajects->begin(); it_mcpart != trajects->end(); ++it_mcpart) {
            
            if ( (*it_mcpart)->GetProcessName() != "Primary"){
                temp_4P = (*it_mcpart)->GetInitialMomentum();
                particle_PDG = (*it_mcpart)->GetPDGCode();
                particle_ID = (*it_mcpart)->GetTrackId();
                mother_ID = (*it_mcpart)->GetParentId();
                
                if (mother_ID == pi0_ID) {
                    if (isFirstGamma){
                        gamma1_4P[0] = temp_4P.px();  
                        gamma1_4P[1] = temp_4P.py(); 
                        gamma1_4P[2] = temp_4P.pz(); 
                        gamma1_4P[3] = temp_4P.E(); 
                        isFirstGamma = false;
                    }else{
                        gamma2_4P[0] = temp_4P.px();  
                        gamma2_4P[1] = temp_4P.py(); 
                        gamma2_4P[2] = temp_4P.pz(); 
                        gamma2_4P[3] = temp_4P.E(); 
                    }
                } 
            }   
        }
    }
    
    // Make sure Gamma1 is the more energetic one
    if (gamma1_4P[3] < gamma2_4P[3]){
        for(unsigned int i = 0; i < gamma1_4P.size(); i++){
            double temp = gamma1_4P[i];
            gamma1_4P[i] = gamma2_4P[i];
            gamma2_4P[i] = temp;
        }
    }
    
    
    // Check Values
    debug()<<"True Muon 4P = ("<<muon_4P[0]<<", "<<muon_4P[1]<<", "<<muon_4P[2]<<", "<<muon_4P[3]<<")"<<endmsg;
    debug()<<"True Proton 4P = ("<<proton_4P[0]<<", "<<proton_4P[1]<<", "<<proton_4P[2]<<", "<<proton_4P[3]<<")"<<endmsg;
    debug()<<"True Pi0 4P = ("<<pi0_4P[0]<<","<<pi0_4P[1]<<", "<<pi0_4P[2]<<", "<<pi0_4P[3]<<")"<<endmsg;
    debug()<<"True Gamma1 4P = ("<<gamma1_4P[0]<<", "<<gamma1_4P[1]<<", "<<gamma1_4P[2]<<", "<<gamma1_4P[3]<<")"<<endmsg;
    debug()<<"True Gamma2 4P = ("<<gamma2_4P[0]<<", "<<gamma2_4P[1]<<", "<<gamma2_4P[2]<<", "<<gamma2_4P[3]<<")"<<endmsg;
    
    
    //--------------------------------------------------------------------------
    // fill the truthEvent info for Signal Kinematics
    //--------------------------------------------------------------------------
    debug()<<"Filling ntuple for truth Signal Kinematics info!"<<endmsg;
    
    truthEvent->setContainerDoubleData("muon_4P", muon_4P);
    truthEvent->setContainerDoubleData("pi0_4P", pi0_4P);
    truthEvent->setContainerDoubleData("proton_4P", proton_4P);
    truthEvent->setContainerDoubleData("gamma1_4P",  gamma1_4P);
    truthEvent->setContainerDoubleData("gamma2_4P",  gamma2_4P);
    
}

//------------------------------------------------------------------------------
// writeFSParticleTable() writes Final State Particle Table using truth information
//------------------------------------------------------------------------------
void CCProtonPi0::writeFSParticleTable(bool isSignal) const
{     
    Minerva::TG4Trajectories* trajects = NULL;
    if( exist<Minerva::TG4Trajectories>( Minerva::TG4TrajectoryLocation::Default ) ) {
        trajects = get<Minerva::TG4Trajectories>( Minerva::TG4TrajectoryLocation::Default );
    } else {
        warning() << "Could not get TG4Trajectories from " << Minerva::TG4TrajectoryLocation::Default << endmsg;
        return;
    }

    // just to make sure
    if( !trajects ) return;

    Gaudi::LorentzVector temp_4p;
    double gammaE_threshold = 20; //MeV
    int particle_PDG;
    int particle_ID;
    int mother_ID = -1;
    std::vector<int> primary_traj_ID;
    std::vector<int> pri_piplus_ID;
    std::vector<int> pri_piminus_ID;
    
    
    // Using Reverse Iterator to have the table in correct order
    Minerva::TG4Trajectories::reverse_iterator rit_mcpart;
    
    std::cout <<"==============================================================="<<std::endl;
    std::cout <<"Final State Particle Table"<<std::endl;
    if(isSignal){
        std::cout<<">> Marked as Signal! <<"<<std::endl;
    }
    std::cout<<std::left;
    std::cout.width(12); std::cout<<"ID";
    std::cout.width(24); std::cout<<"Process";
    std::cout.width(12); std::cout<<"PDG";
    std::cout.width(12); std::cout<<"Mother";
    std::cout.width(12); std::cout<<"Momentum";
    std::cout.width(12); std::cout<<"Energy";
    std::cout<<std::endl;

    // Print all Primary Particles
    for (rit_mcpart = trajects->rbegin(); rit_mcpart != trajects->rend(); ++rit_mcpart) {
        
        if ( (*rit_mcpart)->GetProcessName() == "Primary"){
            temp_4p = (*rit_mcpart)->GetInitialMomentum();
            particle_PDG = (*rit_mcpart)->GetPDGCode();
            particle_ID = (*rit_mcpart)->GetTrackId();
            
            primary_traj_ID.push_back(particle_ID);
            
            std::cout.width(12); std::cout<<particle_ID;    // ID
            std::cout.width(24); std::cout<<(*rit_mcpart)->GetProcessName();// Process
            std::cout.width(12); std::cout<<particle_PDG;                   // PDG
            std::cout.width(12); std::cout<<(*rit_mcpart)->GetParentId();   // Mother PDG
            std::cout.width(12); std::cout<<temp_4p.P();
            std::cout.width(12); std::cout<<temp_4p.E();
            std::cout<<std::endl;
            
            // Save Primary Trajectory Information
            primary_traj_ID.push_back(particle_ID);
        }
        
    }
    
    
    // Print Special 2nd Particles
    for (rit_mcpart = trajects->rbegin(); rit_mcpart != trajects->rend(); ++rit_mcpart) {

        // Skip Primary Trajectories
        if ( (*rit_mcpart)->GetProcessName() == "Primary") continue;
        
        temp_4p = (*rit_mcpart)->GetInitialMomentum();
        particle_PDG =  (*rit_mcpart)->GetPDGCode();
        mother_ID = (*rit_mcpart)->GetParentId();
        
        
        if ( isMotherPrimary(primary_traj_ID, mother_ID) ){
            if (particle_PDG == -(PDG::muon) ){
                std::cout<<"Michel!"<<std::endl;
                std::cout.width(12); std::cout<<(*rit_mcpart)->GetTrackId();    // ID
                std::cout.width(24); std::cout<<(*rit_mcpart)->GetProcessName();// Process
                std::cout.width(12); std::cout<<particle_PDG;                   // PDG
                std::cout.width(12); std::cout<<(*rit_mcpart)->GetParentId();   // Mother PDG
                std::cout.width(12); std::cout<<temp_4p.P();
                std::cout.width(12); std::cout<<temp_4p.E();
                std::cout<<std::endl;
            }else if(particle_PDG == PDG::pi0){
                std::cout<<"Secondary Pi0!"<<std::endl;
                std::cout.width(12); std::cout<<(*rit_mcpart)->GetTrackId();    // ID
                std::cout.width(24); std::cout<<(*rit_mcpart)->GetProcessName();// Process
                std::cout.width(12); std::cout<<particle_PDG;                   // PDG
                std::cout.width(12); std::cout<<(*rit_mcpart)->GetParentId();   // Mother PDG
                std::cout.width(12); std::cout<<temp_4p.P();
                std::cout.width(12); std::cout<<temp_4p.E();
                std::cout<<std::endl;
            }else if(particle_PDG == PDG::gamma && temp_4p.E() > gammaE_threshold  ){
                std::cout<<"Secondary Gamma with Energy > 20 MeV!"<<std::endl;
                std::cout.width(12); std::cout<<(*rit_mcpart)->GetTrackId();    // ID
                std::cout.width(24); std::cout<<(*rit_mcpart)->GetProcessName();// Process
                std::cout.width(12); std::cout<<particle_PDG;                   // PDG
                std::cout.width(12); std::cout<<(*rit_mcpart)->GetParentId();   // Mother PDG
                std::cout.width(12); std::cout<<temp_4p.P();
                std::cout.width(12); std::cout<<temp_4p.E();
                std::cout<<std::endl;
            }
        }
    }
    std::cout<<"--------------------------------------------------------------"<<std::endl;
  
}

void CCProtonPi0::writeEventRecord(Minerva::GenMinInteraction* truthEvent, bool isSignal) const
{
    //--------------------------------------------------------------------------
    // Get Event Record
    //--------------------------------------------------------------------------
    const Minerva::GenMinEventRecord* eventRecord = truthEvent->eventRecord();
    const std::vector<int> eRecPartID = eventRecord->eRecPartID();
    const std::vector<int> eRecPartStatus = eventRecord->eRecPartStatus();
    const std::vector<int> eRecMother = eventRecord->eRecMother();
    const std::vector<double> eRecPartPx = eventRecord->eRecMomentumPx();
    const std::vector<double> eRecPartPy = eventRecord->eRecMomentumPy();
    const std::vector<double> eRecPartPz = eventRecord->eRecMomentumPz();
    const std::vector<double> eRecPartE = eventRecord->eRecEnergy();
    double momentum = 0.0;
    
    std::cout <<"Event Record"<<std::endl;
    if(isSignal){
        std::cout<<">> Marked as Signal! <<"<<std::endl;
    }
    std::cout<<std::left;
    std::cout.width(4); std::cout<<"ID";
    std::cout.width(12); std::cout<<"PDG";
    std::cout.width(12); std::cout<<"Status";
    std::cout.width(12); std::cout<<"Mother";
    std::cout.width(12); std::cout<<"Momentum";
    std::cout.width(12); std::cout<<"Total E";
    std::cout<<std::endl;
    
    for(unsigned int g_part = 0; g_part < eRecPartID.size(); g_part++ ){
        // Get Momentum
        momentum = sqrt(pow(eRecPartPx[g_part],2) + 
                        pow(eRecPartPy[g_part],2) + 
                        pow(eRecPartPz[g_part],2) );
                        
        // Write Table               
        std::cout.width(4); std::cout<<g_part;
        std::cout.width(12); std::cout<<eRecPartID[g_part];
        std::cout.width(12); std::cout<<eRecPartStatus[g_part];
        std::cout.width(12); std::cout<<eRecMother[g_part];
        std::cout.width(12); std::cout<<momentum;
        std::cout.width(12); std::cout<<eRecPartE[g_part];
        std::cout<<std::endl;
    }
    
    std::cout <<"==============================================================="<<std::endl;
}

void CCProtonPi0::setPi0GenieRecord(Minerva::GenMinInteraction* truthEvent) const
{
    debug() << "Enter CCProtonPi0::setPi0GenieRecord()" << endmsg;
    
    //--------------------------------------------------------------------------
    // Get Event Record
    //--------------------------------------------------------------------------
    const Minerva::GenMinEventRecord* eventRecord = truthEvent->eventRecord();
    const std::vector<int> eRecPartID = eventRecord->eRecPartID();
    const std::vector<int> eRecPartStatus = eventRecord->eRecPartStatus();
    const std::vector<int> eRecMother = eventRecord->eRecMother();
    double t_pi0_status = -9;
    int t_pi0_Mother = -9;
    int t_pi0_MotherStatus= -9;
    int t_pi0_GrandMother= -9;
    int t_pi0_GrandMotherStatus= -9;
    int ind_Mother;
    int ind_GrandMother;
    

    for (unsigned int g_part = 0; g_part < eRecPartID.size(); g_part++){
        if (eRecPartID[g_part] == 111){
            // Set Pi0 Status - Initial Status
            t_pi0_status = eRecPartStatus[g_part];
            
            // Set Pi0 Mother
            ind_Mother = eRecMother[g_part];
            t_pi0_Mother = eRecPartID[ind_Mother];
            t_pi0_MotherStatus = eRecPartStatus[ind_Mother];
            
            // Set Pi0 GrandMother
            ind_GrandMother = eRecMother[ind_Mother];
            t_pi0_GrandMother = eRecPartID[ind_GrandMother];
            t_pi0_GrandMotherStatus = eRecPartStatus[ind_GrandMother];
            
            break;
        }
    }
    
    debug()<<"pi0_status = "<<t_pi0_status<<endmsg;
    debug()<<"pi0_Mother = "<<t_pi0_Mother<<endmsg;
    debug()<<"pi0_MotherStatus = "<<t_pi0_MotherStatus<<endmsg;
    debug()<<"pi0_GrandMother = "<<t_pi0_GrandMother<<endmsg;
    debug()<<"pi0_GrandMotherStatus = "<<t_pi0_GrandMotherStatus<<endmsg;
    
    if (t_pi0_Mother == 2000000001){
        debug()<<"Setting pi0_Mother to = "<<-5<<endmsg; 
        t_pi0_Mother = -5;
    }
    
    truthEvent->setIntData("pi0_status", t_pi0_status);
    truthEvent->setIntData("pi0_Mother", t_pi0_Mother);
    truthEvent->setIntData("pi0_MotherStatus", t_pi0_MotherStatus);
    truthEvent->setIntData("pi0_GrandMother", t_pi0_GrandMother);
    truthEvent->setIntData("pi0_GrandMotherStatus", t_pi0_GrandMotherStatus);
  
}

//------------------------------------------------------------------------------
// isSinglePi0 Function Returns TRUE if there is Single Pi0 in FS or 2xGamma's
//------------------------------------------------------------------------------
bool CCProtonPi0::isSinglePi0(Minerva::GenMinInteraction* truthEvent, int nPi0, int nGamma) const
{
    debug() << "Enter CCProtonPi0::isSinglePi0()" << endmsg;
    
    bool is1Pi0;
    
    if( nPi0 == 1 && nGamma == 0){
        debug()<<"Single Pi0"<<endmsg;
        truthEvent->filtertaglist()->setOrAddFilterTag( "isSignal_1Pi0", true );
        is1Pi0 = true;
    }else if( nPi0 == 0 && nGamma == 2){
        debug()<<"2 Gamma"<<endmsg;
        truthEvent->filtertaglist()->setOrAddFilterTag( "isSignal_2Gamma", true );
        is1Pi0 = true;
    }else{
        debug()<<"No Single Pi0 "<<endmsg;
        is1Pi0 = false;
    }
   
    return is1Pi0;
}
    
//------------------------------------------------------------------------------
// interpret Events which fails the reconstructor cuts
//------------------------------------------------------------------------------
StatusCode CCProtonPi0::interpretFailEvent( Minerva::PhysicsEvent* event ) const
{  
    debug() << "Enter CCProtonPi0::interpretFailEvent()" << endmsg;
    
    Minerva::NeutrinoInt* nuInt = new Minerva::NeutrinoInt(m_anaSignature);
    NeutrinoVect nuInts;
    nuInts.push_back( nuInt );
    markEvent(event);
    addInteractionHyp(event,nuInts);
    fillCommonPhysicsAnaBranches(event);
    fillNuMIBranches(event);
    
    return StatusCode::SUCCESS;
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
    debug() << "Enter CCProtonPi0::setMuonData()" << endmsg;
   
    const double SENTINEL = -9.9;
    // Sanity Check
    if (!m_MuonProng || !m_MuonParticle) {
        fatal()<< "NO Muon!!" <<endmsg;
        return false;
    }

    //--------------------------------------------------------------------------
    // Get Muon Track Information
    //--------------------------------------------------------------------------
    SmartRefVector<Track>::iterator iterTrk;
    SmartRefVector<Track> muonTracks = m_MuonProng->minervaTracks();
    
    int muon_minervaTrack_types = 0;
    int muon_N_minosTracks = 0;
    int nlong = 0;
    int nshort = 0;
    
    // Number of Short and Long Tracks
    for (iterTrk = muonTracks.begin(); iterTrk != muonTracks.end(); ++iterTrk) {
        if ((*iterTrk)->type() == Track::Long) nlong++;
        if ((*iterTrk)->type() == Track::Short) nshort++;
    }
    
    if (nlong > 0 && nshort == 0) muon_minervaTrack_types = 1;
    else if (nlong == 0 && nshort > 0) muon_minervaTrack_types = 2;
    else if (nlong >0 && nshort >0) muon_minervaTrack_types = 3;
    
    // Get number of MINOS Tracks
    muon_N_minosTracks = m_MuonProng->minosTracks().size();   
    
    //--------------------------------------------------------------------------
    // Get Muon Kinematics
    //--------------------------------------------------------------------------
    Gaudi::LorentzVector muon_4p;
    
    int is_minos_track = -1;
    int is_minos_stub = -1;
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
    if ( m_MuonProng->MinosTrack() ) {
        is_minos_track = 1;
        SmartRef<MinosRecoTrack> minosTrack = m_MuonProng->minosTracks()[0];
        muon_qp = minosTrack->qp();
        muon_minosTrackQuality = minosTrack->quality();      
    }
    else if ( m_MuonProng->MinosStub() ) {
        is_minos_stub = 1;
    }
    else{
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
    
    // Write Muon 4-Momentum to Global Variable
    m_muon_4P = muon_4p;
    
    //--------------------------------------------------------------------------
    // Debugging: Check values
    debug()<<"m_muon_4P = ( "<<m_muon_4P.px()<<", "<<m_muon_4P.py()<<", "<<m_muon_4P.pz()<<", "<<m_muon_4P.E()<<" )"<<endmsg;
    debug()<<"P4(Muon) = ( "<<muon_4p.px()<<", "<<muon_4p.py()<<", "<<muon_4p.pz()<<", "<<muon_4p.E()<<" )"<<endmsg;
    //--------------------------------------------------------------------------
    
    debug()<<"Filling Muon Ntuple Variables"<<endmsg;
    
    //--------------------------------------------------------------------------
    // Fill Muon Branches
    //--------------------------------------------------------------------------
    nuInt->setLeptonEnergy( muon_4p );
    
    nuInt->setIntData("muon_hasMinosMatchTrack", is_minos_track );
    nuInt->setIntData("muon_hasMinosMatchStub", is_minos_stub );
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
    
    return true;
}

double CCProtonPi0::Calc_QSq(double Enu) const
{
    const double Mmu = MinervaUnits::M_mu;  // Muon Rest Mass [MeV]
    double Emu = m_muon_4P.E();
    double Pmu_long = Calc_Longitudinal_Momentum(m_muon_4P); 
    
    // Calculate QSq - Use eq. in Research Logbook page 30
    double QSq = 2 * Enu * ( Emu - Pmu_long ) - (Mmu * Mmu);

    return QSq;
}

double CCProtonPi0::Calc_WSq(double Enu, double QSq) const
{
    const double Mn = MinervaUnits::M_n;    // Neutron Rest Mass [MeV]
    double Emu = m_muon_4P.E();

    // Calculate WSq - Use eq. in Research Logbook page 31
    double WSq = QSq + std::pow(Mn,2) + 2 * (Enu - Emu) * Mn; 
    
    return WSq;
}

double CCProtonPi0::Calc_Enu_Cal(double hadronEnergy) const
{
    debug()<<"Calculating Enu_Cal Using:"<<endmsg;
    debug()<<"\tP4(Muon) = "<<m_muon_4P<<endmsg;
    debug()<<"\tHadron Visible Energy = "<<hadronEnergy<<endmsg;

    double Enu = m_muon_4P.E() + hadronEnergy;
    return Enu;
}


// Using Equation from DocDB: 9749 v2
// Back-of-napkin derivations for E_nu of CC(nu, meson) reactions 
double CCProtonPi0::Calc_Enu_1Track() const
{
    debug()<<"Calculating Enu_1Track Using:"<<endmsg;
    debug()<<"\tP4(Muon) = "<<m_muon_4P<<endmsg;
    debug()<<"\tP4(Pi0) = "<<m_pi0_4P<<endmsg;

    // Main Components
    const double Mn = MinervaUnits::M_n;    // Neutron Rest Mass [MeV]
    double Enu;     // Neutrino Energy
    double TN;      // Hadron Kinetic Energy
    double Emu = m_muon_4P.E();             
    double Epi0 = m_pi0_4P.E();
    double Pmu_long = Calc_Longitudinal_Momentum(m_muon_4P); 
    double Ppi0_long = Calc_Longitudinal_Momentum(m_pi0_4P);
    double Px_mu_wrt_Beam = Calc_Px_wrt_Beam(m_muon_4P); 
    double Py_mu_wrt_Beam = Calc_Py_wrt_Beam(m_muon_4P); 
    double Px_pi0_wrt_Beam = Calc_Px_wrt_Beam(m_pi0_4P); 
    double Py_pi0_wrt_Beam = Calc_Py_wrt_Beam(m_pi0_4P); 
    
    // TN Calculation in Parts -- TN is very long to calculate in a single Line
    double Emu_Pmu_long = Emu - Pmu_long;
    double Epi0_Ppi0_long = Epi0 - Ppi0_long;
    double P_long = std::pow((Emu_Pmu_long + Epi0_Ppi0_long),2);
    // Transverse Component (Vectoral Calculation Page 4)
    double Pmu_trans = std::pow(Px_mu_wrt_Beam,2) + std::pow(Py_mu_wrt_Beam,2);
    double Ppi0_trans = std::pow(Px_pi0_wrt_Beam,2) + std::pow(Py_pi0_wrt_Beam,2);
    double Pmu_trans_Ppi0_trans = 2 * ( (Px_mu_wrt_Beam * Px_pi0_wrt_Beam) +
                                        (Py_mu_wrt_Beam * Py_pi0_wrt_Beam) );
    double P_trans = Pmu_trans + Ppi0_trans + Pmu_trans_Ppi0_trans;
    
    // Form Hadron Kinetic Energy 
    TN = 0.5 * (P_long + P_trans) / (Mn - Emu_Pmu_long - Epi0_Ppi0_long);

    Enu = Emu + Epi0 + TN;

    return Enu;
}

double CCProtonPi0::Calc_Enu_2Track() const
{
    debug()<<"Calculating Enu_2Track Using:"<<endmsg;
    debug()<<"\tP4(Muon) = "<<m_muon_4P<<endmsg;
    debug()<<"\tP4(Proton) = "<<m_proton_4P<<endmsg;
    debug()<<"\tP4(Pi0) = "<<m_pi0_4P<<endmsg;

    double Mn = MinervaUnits::M_n; // Neutron Rest Mass [MeV]
    double Emu = m_muon_4P.E();
    double Epi0 = m_pi0_4P.E();
    double Eproton = m_proton_4P.E();

    // Calculate Enu -- Use eq. in Research Log Book page 29
    double Enu = Emu + Epi0 + Eproton - Mn;

    return Enu;
}

double CCProtonPi0::Calc_Longitudinal_Momentum(Gaudi::LorentzVector particle_4P) const
{
    double Theta_wrt_Beam = m_coordSysTool->thetaWRTBeam(particle_4P);
    double P_long = particle_4P.P() * std::cos(Theta_wrt_Beam);
    return P_long;
}

double CCProtonPi0::Calc_Px_wrt_Beam(Gaudi::LorentzVector particle_4P) const
{
    double Theta_wrt_Beam = m_coordSysTool->thetaWRTBeam(particle_4P);
    double Phi_wrt_Beam = m_coordSysTool->phiWRTBeam(particle_4P);
    double Px_wrt_Beam = particle_4P.P() * std::sin(Theta_wrt_Beam) * std::cos(Phi_wrt_Beam); 
    return Px_wrt_Beam;
}

double CCProtonPi0::Calc_Py_wrt_Beam(Gaudi::LorentzVector particle_4P) const
{
    double Theta_wrt_Beam = m_coordSysTool->thetaWRTBeam(particle_4P);
    double Phi_wrt_Beam = m_coordSysTool->phiWRTBeam(particle_4P);
    double Py_wrt_Beam = particle_4P.P() * std::sin(Theta_wrt_Beam) * std::sin(Phi_wrt_Beam); 
    return Py_wrt_Beam;
}

//==============================================================================
// Calculates and Set Event Kinematics
//      Uses 4-Momentum of Final State Particles
//==============================================================================
void CCProtonPi0::setEventKinematics(Minerva::NeutrinoInt* nuInt, double hadronVisibleEnergy) const
{
    debug() << "Enter CCProtonPi0::setEventKinematics()" << endmsg;
        
    double Enu_1Track;  // Neutrino Energy - Only for 1 Track
    double Enu_2Track;  // Neutrino Energy - Only for 2 Track
    double Enu_Cal;     // Neutrino Energy - Calorimetric
    double QSq_1Track;  // Q-Square from Enu_1Track
    double QSq_Cal;     // Q-Square from Enu_Cal
    double WSq_1Track;  // W-Square from Enu_1Track
    double WSq_Cal;     // W-Square from Enu_Cal
    
    //--------------------------------------------------------------------------
    // Calculate Enu, QSq, and WSq using Hadron Calorimetric Energy
    //--------------------------------------------------------------------------
    Enu_Cal = Calc_Enu_Cal(hadronVisibleEnergy);
    QSq_Cal = Calc_QSq(Enu_Cal); 
    WSq_Cal = Calc_WSq(Enu_Cal, QSq_Cal);
 
    //--------------------------------------------------------------------------
    // Calcuate Enu using Event Topology - Different for 1 or 2+ Track Events
    //--------------------------------------------------------------------------
    // 1-Track Events Only
    // Calculate Enu, QSq, and WSq using Muon and Single Meson Energy
    if (m_ProtonParticles.size() == 0){
        Enu_1Track = Calc_Enu_1Track();   
        QSq_1Track = Calc_QSq(Enu_1Track); 
        WSq_1Track = Calc_WSq(Enu_1Track, QSq_1Track);
    }

    // 2+ Track Events with Good Proton 4-Momentum
    // Enu = Emu + Epi0 + Eproton - Mn
    if(m_ProtonParticles.size() > 0){
       Enu_2Track = Calc_Enu_2Track();
    }

    //--------------------------------------------------------------------------
    // Fill NTuples
    //--------------------------------------------------------------------------
    nuInt->setDoubleData("neutrino_E_1Track",Enu_1Track);
    nuInt->setDoubleData("neutrino_E_2Track",Enu_2Track);
    nuInt->setDoubleData("neutrino_E_Cal",Enu_Cal);
    nuInt->setDoubleData("QSq_1Track",QSq_1Track);
    nuInt->setDoubleData("QSq_Cal",QSq_Cal);
    nuInt->setDoubleData("WSq_1Track",WSq_1Track);
    nuInt->setDoubleData("WSq_Cal",WSq_Cal);
    if (WSq_1Track < 0) WSq_1Track = 0;
    if (WSq_Cal < 0) WSq_Cal = 0;
    nuInt->setDoubleData("W_1Track",sqrt(WSq_1Track));   
    nuInt->setDoubleData("W_Cal",sqrt(WSq_Cal));   
}

//==============================================================================
// Find the plane nearest to a point
//==============================================================================
StatusCode CCProtonPi0::getNearestPlane(double z, int & module_return, int & plane_return) const
{
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

    return StatusCode::SUCCESS;
  
}

//==============================================================================
// Created particles for negative bit Minos prongs
//==============================================================================
bool CCProtonPi0::createTrackedParticles(Minerva::ProngVect& prongs ) const
{
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
       
        // Skip prongs with particles
        debug()<< "The prong has n particles = " << prongs[p]->particles().size() << endmsg;
        if( !prongs[p]->particles().empty() ){
            debug()<<"Skipping Prong: Prong has attached particles!"<<endmsg;
            continue;
        }

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

            if (prong->particles().size() == 2){
                Gaudi::LorentzVector P4_0 = prong->particles()[0]->momentumVec();
                Gaudi::LorentzVector P4_1 = prong->particles()[1]->momentumVec();
                debug()<<"ProngParticle[0] P4 = ( "<<P4_0.px()<<", "<<P4_0.py()<<", "<<P4_0.pz()<<", "<<P4_0.E()<<" )"<<endmsg;
                debug()<<"ProngParticle[0] P = "<<P4_0.P()<<endmsg;
                debug()<<"ProngParticle[1] P4 = ( "<<P4_1.px()<<", "<<P4_1.py()<<", "<<P4_1.pz()<<", "<<P4_1.E()<<" )"<<endmsg;
                debug()<<"ProngParticle[1] P = "<<P4_1.P()<<endmsg;
            }
        } else{
            debug() << "Did not make particles for the prong type = " << prongs[p]->typeBitsToString() << endmsg;
        }
        hypotheses.clear();
    }
    
    debug() <<"Exit CCProtonPi0::createTrackedParticles()" << endmsg;
    
    return makeParticles;
}

//==============================================================================
// Return the momentum analyzable contained ( proton candidate ) prong/particle
// Uses Global Variables: m_ProtonProngs and m_ProtonParticles
//==============================================================================
bool CCProtonPi0::getProtonProng(   Minerva::ProngVect& primaryProngs) const
{
    debug() <<"Enter CCProtonPi0::getProtonProng()" << endmsg;
  
    debug() <<"N(primaryProngs) =  " << primaryProngs.size() << endmsg;
    
    // Initialize
    Minerva::ProngVect tmpProngs;
    bool isProtonExist = false;
    bool pass = true;
    std::string tag = "PrimaryProton";
    
    // Get All Proton Candidates
    for(unsigned int p = 0; p < primaryProngs.size(); p++) {
        
        debug() <<"Checking prong "<<p<<endmsg;
        
        // Skip Muon Prong
        if( primaryProngs[p]->filtertaglist()->filterTagExists("PrimaryMuon") ){
            debug() <<"Muon Prong, skipping this prong! "<< endmsg;
            continue;
        }
        
        // Temp Storage for ProngVect, Single Prong
        SmartRef<Minerva::Prong> prong       = (Minerva::Prong*)NULL;
        SmartRef<Minerva::Particle> particle = (Minerva::Particle*)NULL;
                
        // Push current Prong to temp ProngVect
        tmpProngs.push_back( primaryProngs[p] );
        
        // Find Proton using m_protonUtils
        bool isProton = m_protonUtils->findProtonProng(tmpProngs,prong,particle); 
        
        tmpProngs.clear();
        
        if( isProton ) {
            debug() <<"Found a proton prong!"<< endmsg;
            
            Gaudi::LorentzVector P4 = particle->momentumVec(); 
            debug()<<"Proton Particle Score: " << particle->score() << endmsg;
            debug()<<"Proton P4 = ( "<<P4.px()<<", "<<P4.py()<<", "<<P4.pz()<<", "<<P4.E()<<" )"<<endmsg;

            prong->filtertaglist()->setOrAddFilterTag(tag,pass);
            particle->filtertaglist()->setOrAddFilterTag(tag,pass);
            prong->updateBestParticle(particle);
            m_hitTagger->applyColorTag(prong, m_Color_protonProng);
            m_ProtonProngs.push_back( prong );
            m_ProtonParticles.push_back( particle );
            isProtonExist = true;
        }
        
    }

    debug() <<"Found "<<m_ProtonParticles.size()<<" Proton Candidates!"<<endmsg;
   
    debug() <<"Exit CCProtonPi0::getProtonProng()" << endmsg;
    
    return isProtonExist;
}


//==============================================================================
// Set proton particle data
//==============================================================================
bool CCProtonPi0::setProtonData( Minerva::NeutrinoInt* nuInt ) const 
{
    debug() << "Enter CCProtonPi0::setProtonData()" << endmsg;
    
    if ( m_ProtonProngs.size() == 0 ) {
        debug()<< "m_ProtonProngs is empty! Exiting..." <<endmsg;
        return true;
    }
    
    // Sanity Check
    if ( m_ProtonProngs.size() != m_ProtonParticles.size()){
        warning()<< "m_ProtonProngs and m_ProtonParticles NOT SAME SIZE" <<endmsg;
        return false;
    }
  
    const double SENTINEL = -9.9;
    // PID Score from Likelihood Parameters
    Minerva::ParticleVect protonLikelihood;
    Minerva::ParticleVect pionLikelihood;
    std::vector<Minerva::Particle::ID> protonHypotheses;
    std::vector<Minerva::Particle::ID> pionHypotheses;
    protonHypotheses.push_back( Minerva::Particle::Proton );
    pionHypotheses.push_back( Minerva::Particle::Pion );
    
    // Declare Variables
    int leadingProtonIndice = 0;
    double tempProtonMomentum = 0.0;
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
    
    std::vector<double> proton_end_x(10,SENTINEL);
    std::vector<double> proton_end_y(10,SENTINEL);
    std::vector<double> proton_end_z(10,SENTINEL);
    
    std::vector<double> proton_start_x(10,SENTINEL);
    std::vector<double> proton_start_y(10,SENTINEL);
    std::vector<double> proton_start_z(10,SENTINEL);
    std::vector<double> length(10,SENTINEL);
    
    std::vector<int> kinked(10,-1);
    std::vector<int> odMatch(10,-1);
    
    std::vector<double> pionScore(10,SENTINEL);
    std::vector<double> protonScore(10,SENTINEL);
    std::vector<double> protonScoreLLR(10,SENTINEL);
    std::vector<double> chi2(10,SENTINEL);
        
    // Get Vertex Z Position
    double vertexZ = m_PrimaryVertex->position().z();
    
    std::map< std::string, std::vector<double> > dedxMapMomentum;
    std::map< std::string, std::vector<double> > dedxMapScore;

    for(unsigned int i = 0; i < m_dedx_uncertainties.size(); i++) {
        std::string name = m_dedx_uncertainties[i];
        std::vector<double> tmp(10,-1);
        dedxMapMomentum.insert( std::pair< std::string,std::vector<double> >(name,tmp) );
        dedxMapScore.insert( std::pair< std::string,std::vector<double> >(name,tmp) );
    }
    
    // Loop over all Proton Candidates
    for(unsigned int i = 0; i < m_ProtonProngs.size(); i++) {
        SmartRef<Minerva::Prong> prong = m_ProtonProngs[i];
        SmartRef<Minerva::Particle> particle = m_ProtonParticles[i];
        
        double theta = prong->minervaTracks().front()->theta();
        double phi   = prong->minervaTracks().front()->phi();
        
        proton_theta[i]  = m_coordSysTool->thetaWRTBeam(particle->momentumVec());
        proton_thetaX[i] = m_coordSysTool->thetaXWRTBeam(particle->momentumVec());
        proton_thetaY[i] = m_coordSysTool->thetaYWRTBeam(particle->momentumVec());
        proton_phi[i]    = m_coordSysTool->phiWRTBeam(particle->momentumVec());
        
        Gaudi::LorentzVector dEdXprotonfourVec;
        StatusCode sc = m_energyCorrectionTool->getCorrectedEnergy(prong,particle,vertexZ,dEdXprotonfourVec);
        if( sc ) particle->setMomentumVec(dEdXprotonfourVec);
        else dEdXprotonfourVec = particle->momentumVec();
        
        p_dedx[i] = sqrt( std::pow(dEdXprotonfourVec.E(),2) - std::pow(MinervaUnits::M_p,2) );
        correctProtonProngEnergy(prong,p_calCorrection[i],p_visEnergyCorrection[i]);
        
        E[i]  = sqrt( std::pow(p_calCorrection[i],2) + std::pow(MinervaUnits::M_p,2) );
        px[i] = p_calCorrection[i]*sin(theta)*cos(phi);
        py[i] = p_calCorrection[i]*sin(theta)*sin(phi);
        pz[i] = p_calCorrection[i]*cos(theta);
        p[i]  = sqrt( px[i]*px[i] + py[i]*py[i] + pz[i]*pz[i] );
        
        if( prong->minervaTracks().front()->direction() == Minerva::Track::Backward && pz[i] > 0. ) pz[i] *= -1.0;
        Gaudi::LorentzVector protonfourVec(px[i],py[i],pz[i],E[i]);
    
        ekin[i] = protonfourVec.E() - MinervaUnits::M_p;
    
        proton_end_x[i]   = (prong->minervaTracks().back())->lastState().x();
        proton_end_y[i]   = (prong->minervaTracks().back())->lastState().y();
        proton_end_z[i]   = (prong->minervaTracks().back())->lastState().z();
    
        proton_start_x[i] = (prong->minervaTracks().back())->firstState().x();
        proton_start_y[i] = (prong->minervaTracks().back())->firstState().y();
        proton_start_z[i] = (prong->minervaTracks().back())->firstState().z();
    
        length[i] = calcDistance(   proton_end_x[i],proton_end_y[i], proton_end_z[i],
                                    proton_start_x[i],proton_start_y[i],proton_start_z[i]);
        kinked[i]  = (int)prong->Kinked();
        odMatch[i] = (int)prong->OdMatch(); 
        chi2[i]    = particle->getDoubleData("chi2_ndf");
        m_protonUtils->getParticleScores(prong,protonScore[i],pionScore[i]);
        
        m_LikelihoodPIDTool->makeParticles( prong, protonLikelihood, protonHypotheses );
        m_LikelihoodPIDTool->makeParticles( prong, pionLikelihood, pionHypotheses );

        // Phil's Likelihood Method for PID
        if(!protonLikelihood.empty()){
            protonScoreLLR[i] = protonLikelihood[i]->score();
        }else{
            debug()<<"particlesLikelihood is empty!"<<endmsg;
        }
        
        // Find Interaction Proton -- Highest Momentum Proton
        if (p[i]  > tempProtonMomentum){
            tempProtonMomentum = p[i];
            leadingProtonIndice = i;
        }
        
        //----------------------------------------------------------------------
        // Debugging: Check values
            debug()<<"Proton Candidate = "<<i<<" P = "<<p[i]<<endmsg;
            debug()<<"P4(Proton) = ( "<<px[i]<<", "<<py[i]<<", "<<pz[i]<<", "<<E[i]<<" )"<<endmsg;
            debug()<<"  protonScore = "<<protonScore[i]<<" pionScore = "<<pionScore[i]<<endmsg;
            debug()<<"  protonScoreLLR = "<<protonScoreLLR[i]<<endmsg;
            debug()<<"  length[i] = "<<length[i]<<endmsg;
            debug()<<"  ekin[i] = "<<ekin[i]<<endmsg;
            debug()<<"  kinked[i] = "<<kinked[i]<<endmsg;
            debug()<<"  odMatch[i] = "<<odMatch[i]<<endmsg;
        //----------------------------------------------------------------------
    }
    
    // Set Leading Proton 4-Momentum
    m_proton_4P.SetPxPyPzE( px[leadingProtonIndice], py[leadingProtonIndice], pz[leadingProtonIndice],E[leadingProtonIndice]);

    nuInt->setContainerDoubleData("proton_p_calCorrection",p_calCorrection);
    nuInt->setContainerDoubleData("proton_p_visEnergy",p_visEnergyCorrection);
    nuInt->setContainerDoubleData("proton_p_dEdXTool",p_dedx);
    
    nuInt->setContainerDoubleData("proton_endPointX",proton_end_x);
    nuInt->setContainerDoubleData("proton_endPointY",proton_end_y);
    nuInt->setContainerDoubleData("proton_endPointZ",proton_end_z);
    
    nuInt->setContainerDoubleData("proton_startPointX",proton_start_x);
    nuInt->setContainerDoubleData("proton_startPointY",proton_start_y);
    nuInt->setContainerDoubleData("proton_startPointZ",proton_start_z);
    
    nuInt->setContainerDoubleData("proton_length",length);
    nuInt->setContainerDoubleData("proton_px",px);
    nuInt->setContainerDoubleData("proton_py",py);
    nuInt->setContainerDoubleData("proton_pz",pz);
    nuInt->setContainerDoubleData("proton_E",E);
    nuInt->setContainerDoubleData("proton_p",p);
    
    nuInt->setContainerIntData("proton_kinked",kinked);
    nuInt->setContainerIntData("proton_odMatch",odMatch);
    
    nuInt->setContainerDoubleData("protonScore",protonScore);
    nuInt->setContainerDoubleData("pionScore",pionScore);
    nuInt->setContainerDoubleData("protonScore_LLR",protonScoreLLR);
    nuInt->setContainerDoubleData("proton_chi2_ndf",chi2);
    
    nuInt->setContainerDoubleData("proton_theta",proton_theta);
    nuInt->setContainerDoubleData("proton_thetaX",proton_thetaX);
    nuInt->setContainerDoubleData("proton_thetaY",proton_thetaY);
    nuInt->setContainerDoubleData("proton_phi",proton_phi);
    
    nuInt->setContainerDoubleData("proton_ekin",ekin);
    
    debug() << "Exit CCProtonPi0::setProtonData()" << endmsg;
    
    return true;
}

//==============================================================================
// Correct Proton Energy
//==============================================================================
void CCProtonPi0::correctProtonProngEnergy(  SmartRef<Minerva::Prong>& protonProng, 
                                                double& p_calCorrection, 
                                                double& p_visEnergyCorrection ) const
{
    debug() << "Enter CCProtonPi0::correctProtonProngEnergy" << endmsg;

    //-- initialize
    double calEnergy = 0.0;
    double visEnergy = 0.0;
    
    //-- get the proton prong subprongs
    Minerva::ProngVect subProngs = protonProng->subProngs();
    debug()<<"nSubProngs = "<<subProngs.size()<<endmsg; 
    if( !subProngs.empty() ) {
        debug() << "  The number of subProngs = " << subProngs.size() << endmsg;

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
            debug() << "  The calorimetric energy = " << calEnergy << " MeV" << endmsg;

            //-- get the total visible energy of the sub prongs
            Minerva::IDClusterVect clusters = subProngs[p]->getAllIDClusters();
            for(unsigned int clus = 0; clus < clusters.size(); clus++) visEnergy += clusters[clus]->energy();
            debug() << "  The visible energy = " << visEnergy << " MeV" << endmsg;
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
        debug() << "   particle idcode = " << particles[part]->idcode() << ", and four momomentum = " << fourMomentum << endmsg;

        double E_calCorrection = fourMomentum.E() + calEnergy;
        p_calCorrection = sqrt( pow(E_calCorrection,2) - pow(particles[part]->mass(),2) );
        debug() << "  update energy using calorimetric correction = " << E_calCorrection << endmsg;

        double E_visEnergyCorrection = fourMomentum.E() + visEnergy;
        p_visEnergyCorrection = sqrt( pow(E_visEnergyCorrection,2) - pow(particles[part]->mass(),2) );     
        debug() << " update energy using visible energy correction = " << E_visEnergyCorrection << endmsg;
    }     
        
    debug() << "Exit CCProtonPi0::correctProtonProngEnergy" << endmsg;
   
    return;
}

//==============================================================================
// setPi0Data
//==============================================================================
bool CCProtonPi0::setPi0Data( Minerva::NeutrinoInt* nuInt ) const 
{    
    debug()<<" Enter CCProtonPi0::setPi0Data()"<<endmsg;

    // Sanity Check
    if( m_Pi0Blob1 == NULL || m_Pi0Blob2 == NULL){
        warning()<<"Passed NULL IDBlob"<<endmsg;
        return false;
    }
    
    const Gaudi::XYZPoint& vtx_position = m_PrimaryVertex->position();

    if (m_ApplyAttenuationCorrection) {
        ApplyAttenuationCorrection(m_Pi0Blob1);
        ApplyAttenuationCorrection(m_Pi0Blob2);
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

    m_idHoughBlob->getBlobEnergyTime( m_Pi0Blob1, g1energy, g1trkrevis, g1ecalevis, g1hcalevis, g1scalevis );
    m_idHoughBlob->getBlobEnergyTime( m_Pi0Blob2, g2energy, g2trkrevis, g2ecalevis, g2hcalevis, g2scalevis );

    // Make sure Gamma1 is the more energetic one 
    if (g2energy > g1energy) {
        warning()<<" Gamma2 Energy is Higher than Gamma1"<<endmsg;
    }
   
    // Get Blob Time
    double time1 = m_Pi0Blob1->time();
    double time2 = m_Pi0Blob2->time();

    // Get Blob dEdx
    double dEdx1 = 0.0;
    double dEdx2 = 0.0;
    m_idHoughBlob->idBlobdEdx( m_Pi0Blob1, dEdx1 );
    m_idHoughBlob->idBlobdEdx( m_Pi0Blob2, dEdx2 );

    // Get Blob Directions
    Gaudi::XYZVector direction1 = m_Pi0Blob1->direction(); 
    Gaudi::XYZVector direction2	= m_Pi0Blob2->direction();

    // Save Gamma1 and Gamma2 Starting Position (blob vertex)
    std::vector<double> position1;
    std::vector<double> position2;
    position1.push_back(m_Pi0Blob1->startPoint().x());
    position1.push_back(m_Pi0Blob1->startPoint().y());
    position1.push_back(m_Pi0Blob1->startPoint().z());
    
    position2.push_back(m_Pi0Blob2->startPoint().x());
    position2.push_back(m_Pi0Blob2->startPoint().y());
    position2.push_back(m_Pi0Blob2->startPoint().z());
    
    // Calculate Distance between Gamma Vertex and Event Vertex
    double gamma1_dist_vtx = calcDistance(  vtx_position.X(), vtx_position.Y(), vtx_position.Z(),
                                            position1[0], position1[1], position1[2]);
        
    double gamma2_dist_vtx = calcDistance(  vtx_position.X(), vtx_position.Y(), vtx_position.Z(),
                                            position2[0], position2[1], position2[2]);

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
          
    // Calculate Pi0 3-Momentum
    TVector3 pimom = g1mom + g2mom;
    
    // Calculate Opening Angle
    const double openingAngle       = (g1mom.Angle(g2mom))*TMath::RadToDeg();
    const double cos_openingAngle   = direction1.Dot(direction2);
    
    // Calculate invariant Mass of Pi0
    const double invMass = std::sqrt(2*g1energy*g2energy*(1-cos_openingAngle));
    
    // Set Pi0 4 - Momentum
    m_pi0_4P.SetPxPyPzE(pimom.x(),pimom.y(),pimom.z(),g1energy+g2energy);
    
    //--------------------------------------------------------------------------
    // Debugging: Check values
    debug()<<"m_pi0_4P = ( "
        <<m_pi0_4P.px()<<", "
        <<m_pi0_4P.py()<<", "
        <<m_pi0_4P.pz()<<", "
        <<m_pi0_4P.E()<<" )"
        <<endmsg;
    //--------------------------------------------------------------------------

    //--------------------------------------------------------------------------
    //    Fill Branches
    //--------------------------------------------------------------------------

    // Pi0 Information
    nuInt->setDoubleData("pi0_openingAngle", openingAngle);
    nuInt->setDoubleData("pi0_cos_openingAngle", cos_openingAngle );
    nuInt->setDoubleData("pi0_invMass", invMass);
    nuInt->setDoubleData("pi0_px" ,pimom.x());
    nuInt->setDoubleData("pi0_py", pimom.y());
    nuInt->setDoubleData("pi0_pz", pimom.z());
    nuInt->setDoubleData("pi0_E", g1energy+g2energy);
    nuInt->setDoubleData("pi0_theta",  pimom.Theta()*TMath::RadToDeg());
    nuInt->setDoubleData("pi0_phi",    pimom.Phi()*TMath::RadToDeg());
    nuInt->setDoubleData("pi0_thetaX", std::atan2(pimom.X(),pimom.Z())*TMath::RadToDeg());
    nuInt->setDoubleData("pi0_thetaY", std::atan2(pimom.Y(),pimom.Z())*TMath::RadToDeg());
    
    // Gamma1 Information
    nuInt->setDoubleData("gamma1_px",g1mom.x());
    nuInt->setDoubleData("gamma1_py", g1mom.y());
    nuInt->setDoubleData("gamma1_pz", g1mom.z());
    nuInt->setDoubleData("gamma1_E", g1energy);
    nuInt->setDoubleData("gamma1_theta",g1mom.Theta()*TMath::RadToDeg());
    nuInt->setDoubleData("gamma1_phi",  g1mom.Phi()*TMath::RadToDeg());
    nuInt->setDoubleData("gamma1_dEdx", dEdx1 );
    nuInt->setDoubleData("gamma1_time", time1 );
    nuInt->setDoubleData("gamma1_dist_vtx",gamma1_dist_vtx);
    nuInt->setContainerDoubleData("gamma1_direction", direc_1 );
    nuInt->setContainerDoubleData("gamma1_vertex", position1 );
    nuInt->setDoubleData("gamma1_evis_trkr", g1trkrevis);
    nuInt->setDoubleData("gamma1_evis_ecal", g1ecalevis);
    nuInt->setDoubleData("gamma1_evis_hcal", g1hcalevis);
    nuInt->setDoubleData("gamma1_evis_scal", g1scalevis);

    // Gamma2 Information
    nuInt->setDoubleData("gamma2_px", g2mom.x());
    nuInt->setDoubleData("gamma2_py", g2mom.y());
    nuInt->setDoubleData("gamma2_pz", g2mom.z());
    nuInt->setDoubleData("gamma2_E", g2energy);
    nuInt->setDoubleData("gamma2_theta", g2mom.Theta()*TMath::RadToDeg());
    nuInt->setDoubleData("gamma2_phi",  g2mom.Phi()*TMath::RadToDeg());
    nuInt->setDoubleData("gamma2_dEdx", dEdx2 );
    nuInt->setDoubleData("gamma2_time", time2 );
    nuInt->setDoubleData("gamma2_dist_vtx", gamma2_dist_vtx);
    nuInt->setContainerDoubleData("gamma2_direction", direc_2 );
    nuInt->setContainerDoubleData("gamma2_vertex", position2 );
    nuInt->setDoubleData("gamma2_evis_trkr", g2trkrevis);
    nuInt->setDoubleData("gamma2_evis_ecal", g2ecalevis);
    nuInt->setDoubleData("gamma2_evis_hcal", g2hcalevis);
    nuInt->setDoubleData("gamma2_evis_scal", g2scalevis);
    
    return true;
}


//------------------------------------------------------------------------------
// set the track prong Geant4 truth information
//------------------------------------------------------------------------------
void CCProtonPi0::setTrackProngTruth( Minerva::NeutrinoInt* neutrino, Minerva::ProngVect& prongs ) const
{
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
                    } // end of while loop
                } // end of else
            } // end loop over traj 
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
    } // end loop over prongs
    
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
   
    return;
}

//------------------------------------------------------------------------------
// Returns the "Used" or "Unused" Energy of the Clusters
//------------------------------------------------------------------------------
double CCProtonPi0::getClusterEnergy(Minerva::PhysicsEvent* event, std::string input_clusterType) const
{
    debug() <<"Enter: CCProtonPi0::getClusterEnergy()" << endmsg;
    
    const std::string clusterType1 = "Unused";
    const std::string clusterType2 = "Used";
    const std::string clusterType3 = "Used:Unused";
    const std::string energyType = "!XTalkCandidate"; // energy from hits that are NOT a Cross-Talk Candidate
    std::string clusterType;
    
    if(input_clusterType.compare(clusterType1) == 0 || input_clusterType.compare(clusterType2) == 0 || input_clusterType.compare(clusterType3) == 0 ){
        clusterType = input_clusterType;
    }else{
        warning()<<"Wrong usage for getClusterEnergy()"<<endmsg;
        warning()<<"Cluster Type can be only = Used, Unused or Used:Unused"<<endmsg;
        return -1.0;
    }
   
    
    SmartRefVector<IDCluster> idClusters = event->select<IDCluster>(clusterType,energyType);
    SmartRefVector<ODCluster> odClusters = event->select<ODCluster>(clusterType,energyType);  
    double idEnergy = m_extraEnergyTool->getIDEnergy(idClusters, -1.0);
    double odEnergy = m_extraEnergyTool->getODEnergy(odClusters, -1.0);
    double totalEnergy = idEnergy  + odEnergy;
   
    debug()<<"idEnergy = "<<idEnergy<<endmsg;
    debug()<<"odEnergy = "<<odEnergy<<endmsg;

    debug()<<"Exit: CCProtonPi0::getClusterEnergy() "<<idEnergy<<endmsg;
    
    return totalEnergy;
}

//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
//
//  Pi0 Reconstruction Functions 
//
//<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

//==============================================================================
//  VtxBlob
//==============================================================================
void CCProtonPi0::VtxBlob(Minerva::PhysicsEvent *event) const
{
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
            m_hitTagger->applyColorTag( vtxSphereBlob, m_Color_VertexFila );
            
        }
      
    } 

  vertex_energy = vertex_energy_sphere + vertex_energy_filament;
  event->setDoubleData( "Vertex_blob_energy", vertex_energy );
  event->setDoubleData( "Filament_Vertex_energy", vertex_energy_filament );
  event->setDoubleData( "Sphere_Vertex_energy", vertex_energy_sphere );
    

}

//==============================================================================
//  FilterInSphereClusters()
//==============================================================================
SmartRefVector<Minerva::IDCluster> CCProtonPi0::FilterInSphereClusters( const SmartRefVector<Minerva::IDCluster>& clusters,
                                     const double sphereRadius,
                                     std::vector<double>& radii) const
{
    debug() << "Enter CCProtonPi0::FilterInSphereClusters()" << endmsg;
    

    const Gaudi::XYZPoint& vtx_position = m_PrimaryVertex->position();
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
  
    return sphereClusters;
    
}

//==============================================================================
//  ConeBlobs
//==============================================================================
bool CCProtonPi0::ConeBlobs( Minerva::PhysicsEvent *event ) const
{
    debug() << "Enter CCProtonPi0::ConeBlobs()" << endmsg;

    // Initialize Bool Variables
    bool isAngleScan        = false;
    bool isAngleScanApplied = false;
    bool additional         = false;
    bool isHough            = false;
    bool isHoughApplied     = false;
    
    // Get Vertex Position
    const Gaudi::XYZPoint& vtx_position = m_PrimaryVertex->position();
    
    // variables to clean up the clusters
    SmartRefVector<Minerva::IDCluster> preidClusters 
            = event->select<Minerva::IDCluster>("Unused","!LowActivity&!XTalkCandidate");
    SmartRefVector<Minerva::IDCluster> usableClusters;
    SmartRefVector<Minerva::IDCluster> outTimeClusters;
    
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

    std::vector<Minerva::IDBlob*> angleScanBlobs = angleScanAlg.GetShowers();
    for (   std::vector<Minerva::IDBlob*>::const_iterator b = angleScanBlobs.begin();
            b != angleScanBlobs.end(); 
            ++b) 
    {
        if ( !m_idHoughBlob->isPhoton((*b)->clusters(),vtx_position)) continue;
        foundBlobs.push_back(*b);
    }

    nblob_anglescan = foundBlobs.size();
    isAngleScan = foundBlobs.size() == 2;
    isAngleScanApplied = true;

    debug()<<"foundBlobs.size() = "<<foundBlobs.size()<<endmsg;

    if (foundBlobs.size() == 1){
        m_hitTagger->applyColorTag( (foundBlobs)[0], m_Color_Gamma1Prong ); // light green
    }else if (foundBlobs.size() == 2){
        m_hitTagger->applyColorTag( (foundBlobs)[0], m_Color_Gamma1Prong ); // light green
        m_hitTagger->applyColorTag( (foundBlobs)[1], m_Color_Gamma2Prong ); // light blue   
    }else if (foundBlobs.size() > 2){
        m_hitTagger->applyColorTag( (foundBlobs)[0], m_Color_Gamma1Prong ); // light green
        m_hitTagger->applyColorTag( (foundBlobs)[1], m_Color_Gamma2Prong ); // light blue
        for(unsigned int i = 2; i < foundBlobs.size(); i++ ){
            m_hitTagger->applyColorTag( (foundBlobs)[i], m_Color_GammaOtherProng ); // orange 
        }
    }
  

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
            }
        } 

        additional = ( bestBlobs.size() == 2 );
        if (bestBlobs.size() == 2) foundBlobs.swap(bestBlobs);
        
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
    }
   


    event->filtertaglist()->setOrAddFilterTag( "is_anglescan", isAngleScan );
    event->filtertaglist()->setOrAddFilterTag( "is_anglescan_applied",isAngleScanApplied);
    event->filtertaglist()->setOrAddFilterTag( "is_houghtransform", isHough );
    event->filtertaglist()->setOrAddFilterTag( "is_houghtransform_applied",isHoughApplied);

    //storing rejected id Clusters
    Minerva::IDBlob* rejectedBlob = new Minerva::IDBlob();
    
    if (!outTimeClusters.empty()) {
        m_blobUtils->insertIDClusters(outTimeClusters, rejectedBlob, Minerva::IDBlob::DispersedIDBlobPatRec );
        debug()<< "Adding rejected blob with vis energy = "  << rejectedBlob->energy() << endmsg;
        addObject( event, rejectedBlob );
        m_hitTagger->applyColorTag( rejectedBlob, m_Color_RejectedBlob ); // red
    }

    event->setDoubleData( "Rejected_blob_vis_energy", rejectedBlob->energy() );
 
    /* 
        if either isAngleScan or isHough is true, the blobs are finally stored
        in the TES, which are managed. Otherwise, delete them 
    */
    if (isAngleScan || isHough || additional) {
  
        processBlobs(event,foundBlobs);
        ODActivity(event,foundBlobs);
        
        // Save FoundBlobs to Global Variables (Energetic Blob is Blob1)
        if( foundBlobs[0]->energy() >= foundBlobs[1]->energy()){
            m_Pi0Blob1 = foundBlobs[0];
            m_Pi0Blob2 = foundBlobs[1];
        }else{
            m_Pi0Blob1 = foundBlobs[1];
            m_Pi0Blob2 = foundBlobs[0];
        }
    
        return true;
    }else{
        // Clear foundBlobs Vector
        std::vector<Minerva::IDBlob*>::iterator b;
        for (b = foundBlobs.begin(); b != foundBlobs.end(); ++b) { 
            delete *b;
        }
        foundBlobs.clear();
        
        return false;
    }
}

bool CCProtonPi0::AreBlobsGood() const
{
    debug() << "CCProtonPi0::AreBlobsGood()" << endmsg;
   
    // Sanity Check
    if( m_Pi0Blob1 == NULL || m_Pi0Blob2 == NULL){
        warning()<<"Passed NULL m_Pi0Blob"<<endmsg;
        return false;
    }

    const Gaudi::XYZPoint& vtx_position = m_PrimaryVertex->position();
    
    // Check Gamma1 Quality
    bool goodPosition1  = m_idHoughBlob->GetStartPosition(m_Pi0Blob1, vtx_position, true );
    bool goodDirection1 = m_idHoughBlob->GetDirection(m_Pi0Blob1, vtx_position );
    bool isGoodBlob1 = false;
    if (goodPosition1 && goodDirection1) isGoodBlob1 = true;
    
    // Check Gamma2 Quality
    bool goodPosition2  = m_idHoughBlob->GetStartPosition(m_Pi0Blob2, vtx_position, true );
    bool goodDirection2 = m_idHoughBlob->GetDirection(m_Pi0Blob2, vtx_position );
    bool isGoodBlob2 = false;
    if (goodPosition2 && goodDirection2) isGoodBlob2 = true;
   

    if (isGoodBlob1 && isGoodBlob2)  return true;
    else return false;
}

void CCProtonPi0::setBlobData( Minerva::PhysicsEvent* event) const
{
    debug()<<"CCProtonPi0::setBlobData()"<<endmsg;

    // Sanity Check
    if( m_Pi0Blob1 == NULL || m_Pi0Blob2 == NULL){
        warning()<< "Pi0Blobs are NULL!"<<endmsg;
    }

    // Calculate dEdX for Blobs
    Calculate_dEdx(m_Pi0Blob1,event,1);
    Calculate_dEdx(m_Pi0Blob2,event,2);

    // Calculate Min Separation from Vertex
    double blob1_minsep = CalcMinBlobSeparation(m_Pi0Blob1);
    double blob2_minsep = CalcMinBlobSeparation(m_Pi0Blob2);
   
    event->setDoubleData("gamma1_blob_minsep", blob1_minsep);
    event->setDoubleData("gamma1_blob_energy", m_Pi0Blob1->energy());
    event->setIntData("gamma1_blob_nclusters", m_Pi0Blob1->nclusters());
    event->setIntData("gamma1_blob_ndigits", m_Pi0Blob1->getAllDigits().size());
        
    event->setDoubleData("gamma2_blob_minsep", blob2_minsep);
    event->setDoubleData("gamma2_blob_energy", m_Pi0Blob2->energy());
    event->setIntData("gamma2_blob_nclusters", m_Pi0Blob2->nclusters());
    event->setIntData("gamma2_blob_ndigits", m_Pi0Blob2->getAllDigits().size());

}

//==============================================================================
//  HoughBlob
//==============================================================================
StatusCode CCProtonPi0::HoughBlob(SmartRefVector<Minerva::IDCluster> idClusters,
                                   std::vector<Minerva::IDBlob*>& outBlobs) const
{
    debug() << "Enter CCProtonPi0::HoughBlob()" << endmsg;
 
    const Gaudi::XYZPoint& vtx_position = m_PrimaryVertex->position();
    
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


    debug() << " Hough Transform is done! " << endmsg;
   
    return StatusCode::SUCCESS;
}

//==============================================================================
//  processBlobs
//==============================================================================
void CCProtonPi0::processBlobs( Minerva::PhysicsEvent *event, std::vector<Minerva::IDBlob*> idBlobs) const
{
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
    
    debug() << "pi0 candidate" << endmsg;
    debug() << " photon 1 is blob: " << idBlobs[0]->key() << endmsg;
    debug() << " photon 2 is blob: " << idBlobs[1]->key() << endmsg;
    idBlobs[0]->setIntData("Photon1", true);
    idBlobs[1]->setIntData("Photon2", true);
    
    std::pair<int,double> result1 = OneParLineFitBlob(idBlobs[0]);
    std::pair<int,double> result2 = OneParLineFitBlob(idBlobs[1]);

    event->setIntData("blob_ndof_1",result1.first);
    event->setIntData("blob_ndof_2",result2.first);
    event->setDoubleData("blob_fval_1", result1.second);
    event->setDoubleData("blob_fval_2", result2.second);
}

//==============================================================================
//  ApplyAttenuationCorrection
//==============================================================================
void CCProtonPi0::ApplyAttenuationCorrection(Minerva::IDBlob* blob) const
{
    debug() << "Enter CCProtonPi0::ApplyAttenuationCorrection()" << endmsg;

    debug() << "AttenuationCorrection (before): " << blob->energy() << endmsg;

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
    
    debug() << "AttenuationCorrection (after):  " << blob->energy() << endmsg;
    
}

//==============================================================================
//  CalcMinBlobSeparation
//==============================================================================
double CCProtonPi0::CalcMinBlobSeparation(const Minerva::IDBlob* blob) const
{
    debug() << "Enter CCProtonPi0::CalcMinBlobSeparation()" << endmsg;
    
    const Gaudi::XYZPoint& vtx_position = m_PrimaryVertex->position();
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

    return dmin;
}

//==============================================================================
//  GetBlobClusterEnergy
//==============================================================================
std::vector<double> CCProtonPi0::GetBlobClusterEnergy(const Minerva::IDBlob* blob) const
{
    debug() <<"Enter CCProtonPi0::GetBlobClusterEnergy()" << endmsg;
    
    const Gaudi::XYZPoint& vtx_position = m_PrimaryVertex->position();
    
    std::vector<double> clusterEnergies;
    
    SmartRefVector<Minerva::IDCluster> clusters = blob->clusters();

    debug()<<"Vertex to increasingDistanceFromVertex = ("<<vtx_position.X()<<","<<vtx_position.Y()<<","<<vtx_position.Z()<<")"<<endmsg;
    
    std::sort(clusters.begin(), clusters.end(), increasingDistanceFromVertex(vtx_position));

    for (SmartRefVector<Minerva::IDCluster>::iterator c = clusters.begin();
         c != clusters.end(); ++c) {
        clusterEnergies.push_back((*c)->energy());
        if (std::distance(clusters.begin(),c) > 4) break;
    }
   
    return clusterEnergies;
}

//==============================================================================
//  ODActivity
//==============================================================================
StatusCode CCProtonPi0::ODActivity( Minerva::PhysicsEvent *event, std::vector<Minerva::IDBlob*> idBlobs ) const
{
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
     
    return StatusCode::SUCCESS;
}

//==============================================================================
//  Calculate_dEdx
//==============================================================================
void CCProtonPi0::Calculate_dEdx( const Minerva::IDBlob* blob, Minerva::PhysicsEvent* event,  unsigned int blob_number) const
{
    debug() <<"Enter CCProtonPi0::Calculate_dEdx()" << endmsg;
   
    //Sanit Check
    if ( blob == NULL ){
        warning()<<"Passed NULL Blob"<<endmsg;
    }

    const Gaudi::XYZPoint& vtx_position = m_PrimaryVertex->position();
    
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
    
    for ( std::map<int, std::vector<double> >::iterator pl = planeClusterEnergyMap.begin(); pl != planeClusterEnergyMap.end(); ++pl) {
        const int plane = pl->first;
        const std::vector<double>& clusterEnergies = pl->second;
        
        std::map<int, std::vector<double> >::iterator next_pl = pl; ++next_pl;
        const int next_plane = next_pl->first;

        const double total_plane_energy = std::accumulate(clusterEnergies.begin(), clusterEnergies.end(), 0.0);
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
    std::map<int, std::vector<double> >::reverse_iterator pl; 
    for ( std::map<int, std::vector<double> >::reverse_iterator pl = planeClusterEnergyMap.rbegin(); pl != planeClusterEnergyMap.rend(); ++pl) {
        const std::vector<double>& clusterEnergies = pl->second;
        const double total_plane_energy = std::accumulate(clusterEnergies.begin(), clusterEnergies.end(), 0.0);
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
}

//==============================================================================
//  PreFilterPi0
//==============================================================================
bool CCProtonPi0::PreFilterPi0(Minerva::PhysicsEvent *event) const
{
    debug() <<"Enter CCProtonPi0::PreFilterPi0()" << endmsg;
    
    const Gaudi::XYZPoint& vtx_position = m_PrimaryVertex->position();

    // Visible Energy Limits for the PreFilter [MeV]
    // Other = Tracker + ECAL + HCAL
    const double max_Evis_Target = 20;
    const double min_Evis_Other = 50; // energy smaller than  50 Mev must be ignored - Can't be a Pi0
    const double max_Evis_Other = 2500; // energy bigger than 2500 MeV must be ignored - DIS Events
    
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

        nearvtx_total += energy;
    }

    //--------------------------------------------------------------------------
    // Save Visible Energy of all UNUSED Hits for each sub detector
    // Target, Tracker, ECAL, HCAL, TOTAL
    //--------------------------------------------------------------------------
    double ntgtEvis = 0.0;
    double trkrEvis = 0.0;
    double ecalEvis = 0.0;
    double hcalEvis = 0.0;
    double total = 0.0;
    
    for (   SmartRefVector<Minerva::IDCluster>::iterator iter_cluster = unusedClusters.begin();
            iter_cluster != unusedClusters.end(); 
            ++iter_cluster ){
            
        const double energy = (*iter_cluster)->energy();
        Minerva::IDCluster::Subdet subdet = (*iter_cluster)->subdet();
        
        // Total Visible Energy in complete MINERvA Detector
        total += energy;
        
        // Visible Energy in SubDetectors
        if (subdet == Minerva::IDCluster::NuclTargs) ntgtEvis += energy;
        else if (subdet == Minerva::IDCluster::Tracker) trkrEvis += energy;
        else if (subdet == Minerva::IDCluster::ECAL) ecalEvis += energy;
        else if (subdet == Minerva::IDCluster::HCAL) hcalEvis += energy;
        else debug() <<"Cluster SubDetector does not found!"<<endmsg;
    }
    
    // Visible Energy except Target Region
    const double otherevis = trkrEvis + ecalEvis + hcalEvis;
    
    event->setDoubleData("evis_nearvtx", nearvtx_total);
    event->setDoubleData("evis_total", total);
    event->setDoubleData("evis_NuclearTarget",  ntgtEvis);
    event->setDoubleData("evis_Tracker",  trkrEvis);
    event->setDoubleData("evis_ECAL",  ecalEvis);
    event->setDoubleData("evis_HCAL",  hcalEvis);
    event->setDoubleData("evis_TotalExceptNuclearTarget", otherevis);
    
    
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

/* Do a line fit for blob digits in the X view. Force the fitted line to
   go through the primary vertex
 */

std::pair<int,double> CCProtonPi0::OneParLineFitBlob(const Minerva::IDBlob* blob) const
{
    debug() <<"Enter CCProtonPi0::OneParLineFitBlob()" << endmsg;
    
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
        
        const Minerva::DePlane* plane = m_InnerDetector->getDePlane((*d)->stripid());
        assert(plane);
        if (!plane) continue;

        total += (*d)->normEnergy();
        
        const double z = plane->getZCenter();
        const double x = plane->getTPos((*d)->stripid());
        const double w = (*d)->normEnergy();
        data.push_back(LineFit::Point(z,x,0.0,w));
     
    }
    
    const Gaudi::XYZPoint& pos = m_PrimaryVertex->position();
    const double x0 = pos.X();
    const double z0 = pos.Z();

        /* Track slope near the vertex to calculate x(z) */
    const double tkx = m_MuonTrack->nearestNode(z0)->state().ax();

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

double CCProtonPi0::calcDistance(  double x1, double y1, double z1,
                                   double x2, double y2, double z2) const
{
    double d;
    d = sqrt(   (x1-x2)*(x1-x2) + 
                (y1-y2)*(y1-y2) + 
                (z1-z2)*(z1-z2));
    return d;   
}


bool CCProtonPi0::isMichelProngGood(Minerva::Prong &michelProng) const
{
    debug()<<"Enter CCProtonPi0::isMichelProngGood()"<<endmsg;
    
    const double maxEnergy = 100;   // MeV
    const double minTimeDiff = 300; // ns

    double energy = michelProng.getDoubleData("energy");
    double timeDiff = michelProng.getDoubleData("time_diff");

    if(energy < maxEnergy && timeDiff > minTimeDiff ) return true;
    else return false;
}

void CCProtonPi0::saveMichelProngToNTuple(Minerva::PhysicsEvent* event, Minerva::Prong &michelProng) const
{
    debug()<<"Enter CCProtonPi0::saveMichelProngToNTuple"<<endmsg;
    event->setDoubleData("michelProng_distance",michelProng.getDoubleData("distance"));
    event->setDoubleData("michelProng_energy",michelProng.getDoubleData("energy"));
    event->setDoubleData("michelProng_time_diff",michelProng.getDoubleData("time_diff"));
    event->setDoubleData("michelProng_end_Z",michelProng.getDoubleData("edz"));
    event->setDoubleData("michelProng_begin_Z",michelProng.getDoubleData("bgz"));
}

//---------------------------------------------------------------------------------------
// create short anchored tracks
//---------------------------------------------------------------------------------------
bool CCProtonPi0::createdAnchoredShortTracks( Minerva::PhysicsEvent* event, Minerva::Vertex* vertex, 
                                                   bool make_primary_short_tracks ) const
{
    debug() << "Enter CCProtonPi0::createdAnchoredShortTracks" << endmsg;
    StatusCode sc;
    bool createdTracks = false;

    //-- get the time slice
    Minerva::TimeSlice* timeSlice = getTimeSlice( vertex );
    short slice = timeSlice ? timeSlice->sliceNumber() : -1;
    debug() << "   The slice associated with the event is = " << slice << endmsg;

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
        Minerva::IDClusterVect idclusters = getClusters(event);
        debug() << "  The number of IDClusters = " << idclusters.size() << endmsg;

        if( idclusters.empty() ){ 
            can_make_anchor_tracks = false; 
            break; 
        }
         
        //-- creates short tracks
        std::vector<Minerva::Track*> *tracks = new std::vector<Minerva::Track*>; 
        sc = m_anchoredTracker->createAnchoredTracks(idclusters,vertex,tracks,rawClus,searchDist,startDist,gap,
                                                  applyMS,force3D,range3D,range2D);
        debug() << "  The number of anchor tracks created = " << tracks->size() << endmsg; 

        if( !sc || tracks->empty() ){ 
            tracks->clear(); 
            delete tracks; 
            can_make_anchor_tracks = false; 
            break; 
        }
  
        //-- store tracks in the container
        for(std::vector<Minerva::Track*>::iterator trk = tracks->begin(); trk != tracks->end(); trk++) { 
            SmartRef<Minerva::Track> inTrack(*trk);
            addObject(event,(*trk));
            setTrackDirection(*trk,vertex);
            (*trk)->setTimeSlice(slice);
            (*trk)->setPatRecHistory( Minerva::Track::FourHitPatRec );
            (*trk)->setIntData( "createdTrack",1 );
            vertex->addOutgoingTrack(inTrack);

            Minerva::Vertex* stopPointVertex = new Minerva::Vertex();
            stopPointVertex->setIncomingTrack(inTrack);
            stopPointVertex->setPosition((*trk)->lastState().position());
            stopPointVertex->setTimeSlice(slice);
            stopPointVertex->setVertexFlags(Minerva::Vertex::StopPoint);
            addObject(event,stopPointVertex);
            createdTracks = true;
            debug() << "	The new stop vertex's position = " << stopPointVertex->position() << endmsg;

            trackCont.push_back( (*trk) );
        }
    } //-- end loop for making anchored tracks 

    //-- run the vertex energy short tracker
    bool can_make_vertex_tracks = true;
    while( can_make_vertex_tracks ) {
        Minerva::IDClusterVect idclusters = getClusters(event);
        debug() << "  The number of IDClusters = " << idclusters.size() << endmsg;
 
        if( idclusters.empty() ) { can_make_vertex_tracks = false; break; }
     
        std::vector<Minerva::Track*> *tracks = new std::vector<Minerva::Track*>; 
        Gaudi::XYZVector zDirectionVector(0, 0, 1); 
        Gaudi::XYZPoint vertexPosition = vertex->position();
        sc = m_vertexEnergyStudyTool->reconstructVertexActivity(event,vertexPosition,zDirectionVector,idclusters,tracks); 
        debug() << "  The number of vertex short tracks created = " << tracks->size() << endmsg;
 
        if( !sc || tracks->empty() ) { tracks->clear(); delete tracks; can_make_vertex_tracks = false; break; }

        //-- store tracks in the container
        for(std::vector<Minerva::Track*>::iterator trk = tracks->begin(); trk != tracks->end(); trk++) {
            SmartRef<Minerva::Track> inTrack(*trk);
            addObject(event,(*trk));
            setTrackDirection(*trk,vertex);
            (*trk)->setTimeSlice(slice);
            (*trk)->setPatRecHistory( Minerva::Track::NotAssigned1 );
            (*trk)->setIntData( "createdTrack",1 );
            vertex->addOutgoingTrack(inTrack);

            Minerva::Vertex* stopPointVertex = new Minerva::Vertex();
            stopPointVertex->setIncomingTrack(inTrack);
            stopPointVertex->setPosition((*trk)->lastState().position());
            stopPointVertex->setTimeSlice(slice);
            stopPointVertex->setVertexFlags(Minerva::Vertex::StopPoint);
            addObject(event,stopPointVertex);
            createdTracks = true;
            debug() << "	The new stop vertex's position = " << stopPointVertex->position() << endmsg;

            trackCont.push_back( (*trk) );
        }
    } //-- end loop for making vertex short tracks

    debug() << "	the number of tracks to add to store = " << trackCont.size() << endmsg;

    //-- add tracks to store and vertex object
    for(std::vector<Minerva::Track*>::iterator track = trackCont.begin(); track != trackCont.end(); track++) {

        //-- make primary prongs for the created primary tracks
        if( make_primary_short_tracks ) {

            //-- make new primary track prong 
            Minerva::Prong* newProng = new Minerva::Prong("CCProtonPi0PrimaryTrackProng");
            debug() << "  Made a new Prong with creation signature : " << newProng->creationSignature() << endmsg;

            //-- add track and source vertex to the prong
            (*track)->setHistory(Minerva::Track::Used);
            newProng->add((*track));
            newProng->setSourceVertex(vertex);

            //-- add prong to the store and promote to a primary prong
            addObject(event,newProng);
            event->promoteProngToPrimary(SmartRef<Minerva::Prong>(newProng));
        }
    }

    //-- refit vertex
    if( createdTracks ) m_vertexFitter->fit(vertex); 

    //-- store vertex for the primary short_tracks
    if( make_primary_short_tracks && createdTracks ) {

        std::vector<double> fit_vtx;
        fit_vtx.push_back( vertex->position().x() );
        fit_vtx.push_back( vertex->position().y() );
        fit_vtx.push_back( vertex->position().z() ); 

        event->setContainerDoubleData("fit_vtx",fit_vtx);

        SmartRef<Minerva::Vertex> vtx(vertex);
        event->setInteractionVertex( vtx );
        Minerva::ProngVect prongs = event->primaryProngs();
        for(unsigned int p = 0; p < prongs.size(); p++) prongs[p]->setSourceVertex(vertex);
    }

    debug() << "  The return vertex's position = " << vertex->position() << endmsg;   
    debug() << "Exit CCProtonPi0::createdAnchoredShortTracks" << endmsg;
    
    return createdTracks;
}

//------------------------------------------------------------------------------------
//  grab a list of clusters for creating a short track
//------------------------------------------------------------------------------------
Minerva::IDClusterVect CCProtonPi0::getClusters( Minerva::PhysicsEvent* event ) const
{
    debug() << "Enter CCProtonPi0::getClusters" << endmsg;
    Minerva::IDClusterVect retClusters; 
    retClusters.clear();

    StatusCode sc;
 
    //-- get interaction vertex
    SmartRef<Minerva::Vertex> vertex  = event->interactionVertex();   

    //-- get the clusters
    Minerva::IDClusterVect idClusters = event->select<Minerva::IDCluster>("Unused","!LowActivity&!XTalkCandidate");
    if( idClusters.size() > 200 ){ 
        idClusters.clear(); 
        return retClusters; 
    }
    debug() << "  the number of idClusters before anchoring = " << idClusters.size() << endmsg;
   
    //-- create new anchored ID tracks to retrieve an axis for the cone
    std::vector<Minerva::Track*>* tracks = new std::vector<Minerva::Track*>;
    sc = m_anchoredTracker->createAnchoredTracks(idClusters,vertex,tracks,true,6000,800,1000,false,false,50,850);

    //-- sanity checks for the number of tracks
    debug() << "   the number of tracks = " << tracks->size() << ", nclusters = " << idClusters.size() << endmsg;
    if( !sc || tracks->empty() ){ 
        tracks->clear(); 
        delete tracks; 
        return retClusters; 
    }

    //-- set track direction
    std::vector<Minerva::Track*>::iterator trk;
    for(trk = tracks->begin(); trk != tracks->end(); trk++){ 
        setTrackDirection(*trk,vertex);
    }
   
    //-- retreive the track with the closest first state position to the vertex
    Minerva::Track* track = *tracks->begin();
    if( tracks->size() > 1 ) {
        double min_separation = 9999.9;
        for(trk = tracks->begin(); trk != tracks->end(); trk++) {
            double separation = MathTool->distanceBetween(vertex->position(),(*trk)->firstState().position());
            if( separation < min_separation ){ 
                track = *trk; 
                min_separation = separation; 
            }
        }
    }

    //-- store the track's clusters 
    Minerva::IDClusterVect clusTmp;
    Minerva::Track::NodeContainer::iterator node;
    for(node = track->nodes().begin(); node != track->nodes().end(); node++){
        clusTmp.push_back( (*node)->idcluster() );
    }
   
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
    debug() << "    cone's firstState position = " << firstState << ", axis = " << axis << endmsg;  

    //-- clean up tracks
    for( trk = tracks->begin(); trk != tracks->end(); trk++){ 
        (*trk)->reset(); 
        delete *trk;  
    }
    tracks->clear();
    delete tracks;

    //-- get fill container with clusters inside of cone
    Minerva::IDClusterVect clusters = event->select<Minerva::IDCluster>("Unused","!LowActivity&!XTalkCandidate");
    debug() << "   the number of cone candidate clusters = " << clusters.size() << endmsg;

    for(unsigned int clus = 0; clus < clusters.size(); clus++) {
        if( m_coneUtilsTool->isInsideCone(clusters[clus],cone) ) retClusters.push_back( clusters[clus] ); 
        else { 
            for(unsigned int tmp = 0; tmp < clusTmp.size(); tmp++){ 
                if( clusTmp[tmp] == clusters[clus] ) retClusters.push_back( clusters[clus] );
            }
        }
    }

    clusters.clear();
    clusTmp.clear();

    debug() << "   the number of clusters return from coning = " << retClusters.size() << endmsg;
    debug() << "Exit CCProtonPi0::getClusters" << endmsg;
    
    return retClusters;
}

//----------------------------------------------------------------------------------
// set the created track direction
//----------------------------------------------------------------------------------
void CCProtonPi0::setTrackDirection( Minerva::Track* track, Minerva::Vertex* vertex ) const
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


#endif

