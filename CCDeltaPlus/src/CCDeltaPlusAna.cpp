/*
    See CCDeltaPlusAna.h header for Class Information
*/
#include "CCDeltaPlusAna.h"

#include "TRandom3.h"

#include "RecInterfaces/IFiducialPointTool.h"

#include "EnergyRecTools/IEnergyCorrectionTool.h"

#include "GeoUtils/IMinervaCoordSysTool.h"

#include "ProngMaker/IMichelTool.h"
#include "ProngMaker/IProngClassificationTool.h"
#include "ProngMaker/IODProngClassificationTool.h"

#include "ParticleMaker/IParticleMakerTool.h"

#include "AnaUtils/IProtonUtils.h"
#include "AnaUtils/ICCPionIncUtils.h"

#include "MinervaUtils/IHitTaggerTool.h"

#include "MinervaDet/IGeomUtilSvc.h"
#include "MinervaDet/DeDetector.h"


DECLARE_TOOL_FACTORY( CCDeltaPlusAna );

using namespace Minerva;
    
//==============================================================================
// Standard constructor 
//==============================================================================
CCDeltaPlusAna::CCDeltaPlusAna(const std::string& type, const std::string& name, const IInterface* parent ) :
MinervaAnalysisTool( type, name, parent ) 
{
    info() <<"--------------------------------------------------------------------------"<<endmsg;
    info() <<"Enter CCDeltaPlusAna::CCDeltaPlusAna::CCDeltaPlusAna() -- Default Constructor" << endmsg;
    
    declareInterface<IInteractionHypothesis>(this);
    //! mandatory declaration of analysis signature: CCDeltaPlusAna
    m_anaSignature = "CCDeltaPlusAna";
    
    // Private Properties
    declareProperty( "StoreAllEvents",      m_store_all_events =    true);
    declareProperty( "DoPlausibilityCuts",  m_doPlausibilityCuts =  true );
    declareProperty( "MakeShortTracks",     m_makeShortTracks =     true );

    declareProperty( "BeamAngleBias",       m_beamAngleBias = 0.006*CLHEP::radian );
    
    declareProperty( "MinMuonScore",        m_minMuonScore = 0.9 );
    declareProperty( "MinProtonScore",      m_minProtonScore = 0.1 );
    
    declareProperty("MuonProngColor",              m_muonProngColor       = 0x228B22); //-- green
    declareProperty("ProtonProngColor",            m_protonProngColor     = 0x9932CC); //-- purple
    declareProperty("PrimaryVertexProngColor",     m_primaryVertexColor   = 0xFF0000); //-- red
    declareProperty("SecondaryVertexProngColor",   m_secondaryVertexColor = 0xFFA500); //-- orange
    declareProperty("TrackEndProngColor",          m_endPointVertexColor  = 0xFF1493); //-- pink
    declareProperty("UnattachedProngColor",        m_unattachedProngColor = 0x0000FF); //-- blue
    

    
    
    
    // Protected properties from IInteractionHypothesis.
    m_hypMeths.push_back( m_anaSignature );
    declareProperty("HypothesisMethods", m_hypMeths);
    info() << " CCDeltaPlusAna Hypothesis added " << endmsg;
    
    
    info() <<"Exit CCDeltaPlusAna::CCDeltaPlusAna::CCDeltaPlusAna() -- Default Constructor" << endmsg;
    info() <<"--------------------------------------------------------------------------"<<endmsg;
}
    
//==============================================================================
// Initialize
//==============================================================================
StatusCode CCDeltaPlusAna::initialize()
{
    info() <<"--------------------------------------------------------------------------"<<endmsg;
    info() <<"Enter CCDeltaPlusAna::initialize()" << endmsg;
    
    //! Initialize the base class.  This will fail if you did not define m_anaSignature.
    StatusCode sc = this->MinervaAnalysisTool::initialize();
    if( sc.isFailure() ) { 
        return Error( "Failed to initialize!", sc ); 
    }
    info()<<"   initialized MinervaAnalysisTool"<<endmsg;
    
    //Seed the RNG
    m_randomGen = new TRandom3( m_randomSeed );   
    if ( m_randomGen == NULL ) return StatusCode::FAILURE;
    
    // Fiducial Volume
    m_fidHexApothem  = 850.0*CLHEP::mm;
    m_fidUpStreamZ   = 5990.0*CLHEP::mm;   // ~middle of module 27, plane 1
    m_fidDownStreamZ = 8340.0*CLHEP::mm;   // ~middle of module 79, plane 1
    
    // Analysable Volume
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
    
    try{
        m_coordSysTool = tool<IMinervaCoordSysTool>("MinervaCoordSysTool");
    }catch(GaudiException& e){
        error() << "Could not obtain tool: MinervaCoordSysTool" << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{ 
        m_protonUtils = tool<IProtonUtils>("ProtonUtils", m_protonUtils); 
    }catch( GaudiException& e ){
        error() << "Could not obtain tool: ProtonUtils" << m_protonUtils<< endmsg;
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
        m_michelTrkTool = tool<IMichelTool>("MichelTool","CCDeltaMichelTrackTool");
    } catch(GaudiException& e){
        error()<<"Could not obtain tool: MichelTool" << endmsg;
        return StatusCode::FAILURE;
    } 
    
    try {
        m_michelVtxTool = tool<IMichelTool>("MichelTool","CCDeltaMichelVtxTool");
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
    
    //! Event - Cut Results
    declareIntEventBranch( "Cut_Vertex_None", -1 );
    declareIntEventBranch( "Cut_Vertex_Null", -1 );
    declareIntEventBranch( "Cut_Vertex_Not_Analyzable", -1 );
    declareIntEventBranch( "Cut_Vertex_Not_Fiducial", -1 );
    declareIntEventBranch( "Cut_Muon_None",-1);
    declareIntEventBranch( "Cut_Muon_Not_Plausible",-1);
    declareIntEventBranch( "Cut_Muon_Score_Low",-1);
    declareIntEventBranch( "Cut_Michel_Exist", -1 );
    declareIntEventBranch( "Cut_Proton_None", -1 );
    
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
    
    //! NeutrinoInt - Vertex
    declareIntBranch( m_hypMeths, "vtx_module", -99);
    declareIntBranch( m_hypMeths, "vtx_plane",-1);
    declareDoubleBranch( m_hypMeths, "vtx_x",0.0);
    declareDoubleBranch( m_hypMeths, "vtx_y",0.0);
    declareDoubleBranch( m_hypMeths, "vtx_z",0.0);
    
    //! NeutrinoInt - Muon
    declareIntBranch( m_hypMeths, "muon_minervaTrack_types", -1);
    declareIntBranch( m_hypMeths, "muon_N_minosTracks", -1);
    declareIntBranch( m_hypMeths, "muon_minosTrackQuality", -1);
    declareIntBranch( m_hypMeths, "muon_roadUpstreamPlanes", -1);
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
    
    //! NeutrinoInt - Proton
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
    
    
    info() <<"Exit CCDeltaPlusAna::initialize()" << endmsg;
    info() <<"--------------------------------------------------------------------------"<<endmsg;
    
    return sc;
}
    
//==============================================================================
// reconstructEvent() --
//==============================================================================
StatusCode CCDeltaPlusAna::reconstructEvent( Minerva::PhysicsEvent *event, Minerva::GenMinInteraction* truthEvent ) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() <<"Enter CCDeltaPlusAna::reconstructEvent()" << endmsg;
    
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
    
    //--------------------------------------------------------------------------
    //! MAKE CUT - Get the interaction vertex, if it exists
    //--------------------------------------------------------------------------
        debug() << "START: Vertex Reconstruction..." << endmsg;
    if( !(event->hasInteractionVertex()) ) {
        debug() << "The event does not have an interaction vertex!" << endmsg;
        event->setIntData("Cut_Vertex_None",1);
        if( m_store_all_events ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    } 
       
    //--------------------------------------------------------------------------
    //! MAKE CUT - Does interactionVertex() is NULL
    //--------------------------------------------------------------------------
    SmartRef<Minerva::Vertex> interactionVertex = event->interactionVertex();
    if( interactionVertex == NULL ) { 
        bool pass = true; 
        std::string tag = "BadObject";
        event->filtertaglist()->addFilterTag(tag,pass);
        event->setIntData("Cut_Vertex_Null",1);
        error() << "This vertex is NULL! Flag this event as bad!" << endmsg;
        if( m_store_all_events ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    }
    
    //--------------------------------------------------------------------------
    //! MAKE CUT - Require vertex analyzability so that we don't run as much reco
    //--------------------------------------------------------------------------
    if ( !FiducialPointTool->isFiducial(interactionVertex->position(), m_recoHexApothem, m_recoUpStreamZ, m_recoDownStreamZ ) ) {
        Gaudi::XYZPoint vtx_position = interactionVertex->position();
        debug() <<"Interaction Vertex is NOT in analyzable volume = ("<<vtx_position.x()<<","<<vtx_position.y()<<","<<vtx_position.z()<<")"<< endmsg;
        if( m_store_all_events ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    }
    
    //--------------------------------------------------------------
    //! Create vertex-anchored short track Prongs, refit vertex
    //--------------------------------------------------------------
    int n_anchored_long_trk_prongs = event->primaryProngs().size() - 1;
    if (m_makeShortTracks) {
        debug()<<"making short tracks"<<endmsg;
        m_ccPionIncUtils->makeShortTracks(event); 
        debug()<<"finished making short tracks"<<endmsg;
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
    
    //--------------------------------------------------------------
    //! MAKE CUT - Check to see if the Interaction Vertex is Fiducial.
    //--------------------------------------------------------------  
    debug()<<"Making fiducial volume cut"<<endmsg;
    
    // check for NaN first
    Gaudi::XYZPoint vtx_position = interactionVertex->position();
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
    
    debug() << "END: Vertex Reconstruction!" << endmsg;

    //--------------------------------------------------------------------------
    //! Get Muon, if it exists
    //--------------------------------------------------------------------------  
    debug() << "Finding Muon..." << endmsg;
    
    SmartRef<Minerva::Prong>    muonProng = (Minerva::Prong*)NULL;
    SmartRef<Minerva::Particle> muonPart = (Minerva::Particle*)NULL;
    
    bool foundMuon = MuonUtils->findMuonProng( event, muonProng, muonPart );
    bool is_minos_track = false, is_minos_stub = false;
    if ( foundMuon ){
        if ( !muonProng ) { 
            warning() << "Identified a muon Prong, but it is NULL!" << endmsg;
            return StatusCode::FAILURE; // We sort of did crash... 
        }
        
        double mc_frac = -1.0;
        if ( m_doPlausibilityCuts && !muonIsPlausible( muonProng, mc_frac) ) {
            debug()<<"Muon is not plausible"<<endmsg;
            event->setIntData("Cut_Muon_Not_Plausible",1);
            if( m_store_all_events ) return interpretFailEvent(event); 
            else return StatusCode::SUCCESS; 
        }
        
        debug() << "Muon Particle Score: " << muonPart->score() << endmsg;
        if (muonPart->score() >= m_minMuonScore){
        
            muonProng->filtertaglist()->setOrAddFilterTag( "PrimaryMuon", true );
            muonPart->filtertaglist()->setOrAddFilterTag( "PrimaryMuon", true );
            
            if (muonProng->MinosTrack()) is_minos_track = true;
            if (muonProng->MinosStub()) is_minos_stub = true;
            if (is_minos_stub && is_minos_track) counter("MuonHasMinosStubAndTrack")++;
            else counter("MuonHasMinosStubAndTrack")+=0;
            if (!is_minos_stub && !is_minos_track) counter("MuonIsNotMinosMatched")++;
            else counter("MuonIsNotMinosMatched")+=0;
            
            event->filtertaglist()->setOrAddFilterTag("isMinosMatchTrack", is_minos_track );
            event->filtertaglist()->setOrAddFilterTag("isMinosMatchStub", is_minos_stub );
            if (truthEvent) truthEvent->filtertaglist()->setOrAddFilterTag( "reco_isMinosMatch", true );
            
            //! @todo - also determine whether there is a minos track in event, if there is no minos match track or stub
        }
        else {
            debug()<<"Muon prong does not pass score cut"<<endmsg;
            event->setIntData("Cut_Muon_Score_Low",1);
            if( m_store_all_events ) return interpretFailEvent(event); 
            else return StatusCode::SUCCESS; 
        }
        
    } 
    else{
        debug() << "Did not find a muon prong!" << endmsg;
        event->setIntData("Cut_Muon_None",1);
        if( m_store_all_events ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    } // end if findMuonProng
    
    
    double muon_visible_energy = muonProng->minervaVisibleEnergySum();
    
    debug() << "Finding Muon End!" << endmsg;
    
    //--------------------------------------------------------------------------
    //! MAKE CUT - Determine if vertex has michels
    //--------------------------------------------------------------------------
    debug() << "Finding Michels..." << endmsg;
    Minerva::Prong vtx_michel_prong;
    bool foundMichel = m_michelVtxTool->findMichel( interactionVertex, vtx_michel_prong );
    if (foundMichel) {
        debug()<<"Found a Michel Electron!"<<endmsg;
        event->setIntData("Cut_Michel_Exist",1);
        event->setIntData("n_vtx_michel_views",vtx_michel_prong.getIntData("category"));
        event->setDoubleData("vtx_michel_distance",vtx_michel_prong.getDoubleData("distance"));
        if( m_store_all_events ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    }else{
        debug()<<"There are NO Michel Electrons in the event!"<<endmsg;
    }
    
    debug() << "Finding Michels End!" << endmsg;
    
    //--------------------------------------------------------------------------
    //! Get Proton, if it exists
    //--------------------------------------------------------------------------
    debug() << "Finding Protons..." << endmsg;
    //-- get all of the primary prongs in the event
    Minerva::ProngVect primaryProngs = event->primaryProngs();
    
    //-- create new particles
    debug() << "Creating particles with Proton and Pion hypotheses" <<endmsg;
    bool makeParticles = createTrackedParticles(primaryProngs);
    debug() << "Was the creation of particles successful? " << makeParticles << endmsg;

    //-- check if one of the primary prong is contained and has a dEdX proton particle 
    Minerva::ProngVect protonProngs;
    Minerva::ParticleVect protonParticles;
    bool foundProton = getProtonProng(primaryProngs, protonProngs, protonParticles);
    if( foundProton ) {
        for(unsigned int i = 0; i < protonProngs.size(); i++) {
            debug() << "Tag the proton's prong with bit-field = " << protonProngs[i]->typeBitsToString() << endmsg;
        }
    } else {
        debug() << "Didn't find any contained in the tracker bit-positive prong with a proton particle!" << endmsg;
        event->setIntData("Cut_Proton_None",1);
        if( m_store_all_events ) return interpretFailEvent(event); 
        else return StatusCode::SUCCESS; 
    }
    
    debug() << "Finding Protons End!" << endmsg;
    

    //--------------------------------------------------------------------------
    //! Finish filling event portion of ntuple 
    //--------------------------------------------------------------------------
//     event->setDoubleData("time", event_time);
//     
//     event->setIntData( "dead", dead );
//     event->setIntData( "udead", udead );
//     event->setIntData( "ddead", ddead );
//     event->setIntData( "tdead", tdead );
//     
//     event->setIntData("n_startpoint_vertices", n_startpoint_vertices);
//     event->setIntData("n_long_tracks", n_long_tracks);
//     event->setIntData("n_short_tracks", n_short_tracks);
    event->setIntData("n_anchored_long_trk_prongs", n_anchored_long_trk_prongs);
    event->setIntData("n_anchored_short_trk_prongs", n_anchored_short_trk_prongs);
//     event->setIntData("n_vtx_blob_prongs", n_vtx_blob_prongs);
    event->setIntData("n_iso_trk_prongs", n_iso_trk_prongs);
//     event->setIntData("n_iso_blob_prongs", n_iso_blob_prongs);
//     event->setIntData("n_dsp_blob_prongs", n_dsp_blob_prongs);
//     
    event->setDoubleData("muonVisibleE", muon_visible_energy );
//     event->setDoubleData("hadronVisibleE", hadron_visible_energy );
//     
//     event->setDoubleData( "totalVisibleE",   totalVisibleEnergy );
//     event->setDoubleData( "totalIDVisibleE", idVisibleEnergy );
//     event->setDoubleData( "totalODVisibleE", odVisibleEnergy );
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
    

    debug() <<"Exit CCDeltaPlusAna::reconstructEvent()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    return sc;
    
}
    
//==============================================================================
// interpretEvent() 
//==============================================================================
StatusCode CCDeltaPlusAna::interpretEvent( const Minerva::PhysicsEvent *event, const Minerva::GenMinInteraction *truthEvent, NeutrinoVect& interaction_hyp ) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() <<"Enter CCDeltaPlusAna::interpretEvent()" << endmsg;
    
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
    //! Get the Primary Muon Prongs
    //--------------------------------------------------------------------------
    SmartRef<Prong> muonProng = (Prong*)NULL;
    ProngVect hadronProngs;
    ProngVect primaryProngs = event->primaryProngs();
    
    for (ProngVect::iterator itProng = primaryProngs.begin(); itProng != primaryProngs.end(); ++itProng) {

        bool isPrimaryMuon = false;
        bool isPrimaryHadron = false;
    
        if ( (*itProng)->filtertaglist()->filterTagExists("PrimaryMuon") ) {
            (*itProng)->filtertaglist()->checkFilterTag( "PrimaryMuon", isPrimaryMuon );
        }
        
//         if ( (*itProng)->filtertaglist()->filterTagExists("PrimaryHadron") ) {
//             (*itProng)->filtertaglist()->checkFilterTag( "PrimaryHadron", isPrimaryHadron );
//         }    
 
        if (isPrimaryMuon && !isPrimaryHadron) muonProng = *itProng;
        if (isPrimaryHadron && !isPrimaryMuon) hadronProngs.push_back(*itProng);
        if (isPrimaryMuon && isPrimaryHadron ) {
            warning()<<"Prong is two primary particles!"<<endmsg;
        }
    }
    
    
    //--------------------------------------------------------------------------  
    //! Get the muon particle hypotheses and fill some info
    //--------------------------------------------------------------------------
    SmartRefVector<Track>::iterator itTrk;
    SmartRef<Minerva::Particle> muonPart = (Minerva::Particle*)NULL;
    int muon_minervaTrack_types = 0;
    int muon_N_minosTracks = 0;
    if (!muonProng) {
        fatal()<<"Where's my muon?!"<<endmsg;
        return StatusCode::FAILURE;
    }
    
    SmartRefVector<Minerva::Particle> particle_hypotheses = muonProng->particles();
    SmartRefVector<Minerva::Particle>::iterator itPart;
    for( itPart = particle_hypotheses.begin(); itPart != particle_hypotheses.end(); ++itPart ) {  
        bool is_muon_particle = false;
        (*itPart)->filtertaglist()->checkFilterTag( "PrimaryMuon", is_muon_particle );
        if (is_muon_particle) {
            muonPart = *itPart;
            break;
        }
    }
    
    
    SmartRefVector<Track> muonTracks = muonProng->minervaTracks();
    int nlong = 0, nshort = 0;
    for (itTrk = muonTracks.begin(); itTrk != muonTracks.end(); ++itTrk) {
        if ((*itTrk)->type() == Track::Long) nlong++;
        if ((*itTrk)->type() == Track::Short) nshort++;
    }
    
    if (nlong > 0 && nshort == 0) muon_minervaTrack_types = 1;
    else if (nlong==0 && nshort > 0) muon_minervaTrack_types = 2;
    else if (nlong>0 && nshort>0) muon_minervaTrack_types = 3;
    
    muon_N_minosTracks = muonProng->minosTracks().size();   
    
    //! muon energy in road upstream info
    int muon_roadUpstreamPlanes = -1;
    double muon_roadUpstreamEnergy = AnaToolUtils->energyInTheRoadUpstream(muonProng, muon_roadUpstreamPlanes);  
    
    //! Fill muon particle hypothesis info
    //! @todo - cleanup redundant minos track information
    double muon_E = -9.0, muon_p = 0.0, muon_px = 0.0, muon_py = 0.0, muon_pz = 0.0;
    double muon_theta = -9.0, muon_theta_biasUp = -9.0, muon_theta_biasDown = -9.0;
    double muon_muScore = -1.0, muon_qp = -9999.9, muon_qpqpe = -9999.9, muon_E_shift;
    int muon_minosTrackQuality = 0;
    Gaudi::LorentzVector muon_4p;
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
    }
    
    debug()<<"Filling Muon Ntuple Variables"<<endmsg;
    
    //! ntuple muon variables
    nuInt->setIntData("muon_minervaTrack_types", muon_minervaTrack_types);
    nuInt->setIntData("muon_N_minosTracks", muon_N_minosTracks);
    nuInt->setIntData("muon_minosTrackQuality", muon_minosTrackQuality);
    nuInt->setIntData("muon_roadUpstreamPlanes", muon_roadUpstreamPlanes);
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
        
    debug()<<"Finished muons"<<endmsg;
    
    //--------------------------------------------------------------------------  
    //! Get the Proton particle hypotheses and fill some info
    //--------------------------------------------------------------------------
    
    // Create Proton Data Storage
    Minerva::ProngVect protonProngs;
    Minerva::ParticleVect protonParticles;
    protonProngs.clear();
    protonParticles.clear();
    
   //-- loop over primaryProngs and retrieve best proton particle
    for(unsigned int prong = 0; prong < primaryProngs.size(); prong++) {
    
        //-- get prong blobs
        Minerva::IDBlobVect blobs = primaryProngs[prong]->idblobs();
    
        //-- proton prong
        if( primaryProngs[prong]->filtertaglist()->filterTagExists("PrimaryProton") ) {
            protonProngs.push_back( primaryProngs[prong] );
            protonParticles.push_back( primaryProngs[prong]->bestParticle() );
            m_hitTagger->applyColorTag(protonProngs.back(),m_protonProngColor);
            for(unsigned int b = 0; b < blobs.size(); b++) {
                if( blobs[b]->history() == Minerva::IDBlob::Hidden ) m_hitTagger->applyColorTag(blobs[b],m_protonProngColor);
            }
        }
    }
    
    // VertexZ Required for setProtonParticleData()
    double vertexZ = event->interactionVertex()->position().z();
    
    setProtonParticleData(nuInt,protonProngs,protonParticles,vertexZ);
    
    //--------------------------------------------------------------------------
    //! Vertex
    //--------------------------------------------------------------------------
    
    //nuInt->setNeutrinoHelicity( getHelicity( mu_charge ) );
    nuInt->setInteractionCurrent( Minerva::NeutrinoInt::ChargedCurrent );
    nuInt->setInteractionType( Minerva::NeutrinoInt::UnknownInt );
    Gaudi::XYZTVector vtx_position( event->interactionVertex()->position().x(), 
                                event->interactionVertex()->position().y(),
                                event->interactionVertex()->position().z(),
                                event->time() );     
    nuInt->setVertex( vtx_position );
    nuInt->setScore( 1.0 );
    nuInt->setLeptonEnergy( muon_4p );
    
    int vtx_module, vtx_plane;
    debug()<<"Calling getNearestPlane, vtx is "<<vtx_position.z()<<endmsg;
    getNearestPlane(vtx_position.z(), vtx_module, vtx_plane); 
    nuInt->setIntData("vtx_module", vtx_module);
    nuInt->setIntData("vtx_plane", vtx_plane);
    nuInt->setDoubleData("vtx_x", vtx_position.x() );
    nuInt->setDoubleData("vtx_y", vtx_position.y() );
    nuInt->setDoubleData("vtx_z", vtx_position.z() );
    

    

    debug() <<"Exit CCDeltaPlusAna::interpretEvent()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    return StatusCode::SUCCESS;
    
}
   

//==============================================================================
// tagTruth()
//==============================================================================
StatusCode CCDeltaPlusAna::tagTruth( Minerva::GenMinInteraction* truthEvent ) const 
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter: tagTruth()" << endmsg;

    
    if( !truthEvent ) {
        warning() << "The GenMinInteraction is NULL!" << endmsg;
        return StatusCode::SUCCESS;
    }
    
    debug() << "Filling Genie Weight Branches" << endmsg;
    //--------------------------------------------------------------------------
    //! Fill GENIE Weight Branches
    //--------------------------------------------------------------------------      
    StatusCode sc = fillGenieWeightBranches( truthEvent );
    if (sc.isFailure() ) {
        warning()<<"Genie weight branch filling failed!"<<endmsg;
        return sc;
    }
    
    debug() <<"Exit CCDeltaPlusAna::CCDeltaPlusAna::tagTruth()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    return StatusCode::SUCCESS;
    
}

//==============================================================================
//  Finalize
//==============================================================================
StatusCode CCDeltaPlusAna::finalize()
{
    info() <<"--------------------------------------------------------------------------"<<endmsg;
    info() << "Enter: finalize()" << endmsg;
    
    
    // finalize the base class.
    StatusCode sc = this->MinervaAnalysisTool::finalize();
    if( sc.isFailure() ) return Error( "Failed to finalize!", sc );
    
    
    info() <<"Exit CCDeltaPlusAna::CCDeltaPlusAna::finalize()" << endmsg;
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
StatusCode CCDeltaPlusAna::interpretFailEvent( Minerva::PhysicsEvent* event ) const
{  

    debug() <<"Exit CCDeltaPlusAna::reconstructEvent() through interpretFailEvent()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter CCDeltaPlusAna::interpretFailEvent()" << endmsg;
    
    Minerva::NeutrinoInt* nuInt = new Minerva::NeutrinoInt(m_anaSignature);
    NeutrinoVect nuInts;
    nuInts.push_back( nuInt );
    markEvent(event);
    addInteractionHyp(event,nuInts);
    fillCommonPhysicsAnaBranches(event);
    fillNuMIBranches(event);
    
    debug() << "Exit CCDeltaPlusAna::interpretFailEvent()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    return StatusCode::SUCCESS;
}

//==============================================================================
// Find the plane nearest to a point
//==============================================================================
StatusCode CCDeltaPlusAna::getNearestPlane(double z, int & module_return, int & plane_return) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter CCDeltaPlusAna::getNearestPlane()" << endmsg;
    
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

    debug() << "Exit CCDeltaPlusAna::getNearestPlane()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    return StatusCode::SUCCESS;
  
}

//==============================================================================
// Created particles for negative bit Minos prongs
//==============================================================================
bool CCDeltaPlusAna::createTrackedParticles( Minerva::ProngVect& prongs ) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() <<"Enter CCDeltaPlusAna::createTrackedParticles()" << endmsg;
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
    
    debug() << "Exit CCDeltaPlusAna::createTrackedParticles()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    return makeParticles;
}

//==============================================================================
// Return the momentum analyzable contained ( proton candidate ) prong/particle
//==============================================================================
bool CCDeltaPlusAna::getProtonProng(    Minerva::ProngVect& primaryProngs, 
                                        Minerva::ProngVect& hadronProngs,
                                        Minerva::ParticleVect& hadronParticles ) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter CCDeltaPlusAna::getProtonProng()" << endmsg;
  
    debug() << "N(primaryProngs) =  " << primaryProngs.size() << endmsg;
    
    // Initialize
    bool isProton = false;
    int nProtons = 0;
    int nGoodProtons = 0; 
    
    // Get All Proton Candidates
    for(unsigned int p = 0; p < primaryProngs.size(); p++) {
        debug() <<"Checking prong "<<p<<endmsg;
        
        if( primaryProngs[p]->filtertaglist()->filterTagExists("PrimaryMuon") ){
            debug() << "Muon Prong, skipping this prong! "<< endmsg;
            continue;
        }
    
        // Temp Storage for ProngVect, Single Prong and Single Particle
        Minerva::ProngVect tmpProngs;
        SmartRef<Minerva::Prong> prong       = (Minerva::Prong*)NULL;
        SmartRef<Minerva::Particle> particle = (Minerva::Particle*)NULL;
        
        // Push current Prong to temp ProngVect
        tmpProngs.push_back( primaryProngs[p] );
        
        // Find Proton using m_protonUtils
        isProton = m_protonUtils->findProtonProng(tmpProngs,prong,particle); 
        
        if( isProton ) {
            nProtons++;
            debug() <<"Found a proton prong!"<< endmsg;
            debug() <<"Proton Particle Score: " << particle->score() << endmsg;
            prong->filtertaglist()->setOrAddFilterTag("PrimaryProton",true);
            particle->filtertaglist()->setOrAddFilterTag("PrimaryProton",true);
            prong->updateBestParticle(particle);
            hadronProngs.push_back( prong );
            hadronParticles.push_back( particle );
        }
    }

    debug() << "Found "<<nProtons<<" Proton Candidates!"<<endmsg;
    debug() << "Applying minProtonScore Cut with = "<<m_minProtonScore<<endmsg;
    
    // Apply Score Cut
    for(unsigned int p = 0; p < hadronParticles.size(); p++) {
        if (hadronParticles[p]->score() > m_minProtonScore){
            nGoodProtons++;
        }else{
            hadronProngs.erase(hadronProngs.begin()+p);
            hadronParticles.erase(hadronParticles.begin()+p);
        }
    }
   
    
    debug() << "Found "<<nGoodProtons<<" Good Protons!"<<endmsg;
    
   
    
    debug() << "Exit CCDeltaPlusAna::getProtonProng()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    
    return isProton;
}

//==============================================================================
// Set proton particle data
//==============================================================================
void CCDeltaPlusAna::setProtonParticleData(     Minerva::NeutrinoInt* nuInt, 
                                                Minerva::ProngVect& protonProngs,
                                                Minerva::ParticleVect& protonParticles, 
                                                double vertexZ ) const 
{
    // Got from Tammy's NukeCCQETwoTrack - 2014_04_15
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter CCDeltaPlusAna::setProtonParticleData()" << endmsg;
    
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
    
    debug() << "Exit CCDeltaPlusAna::setProtonParticleData()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    
    return;
}

//==============================================================================
// Correct Proton Energy
//==============================================================================
void CCDeltaPlusAna::correctProtonProngEnergy(  SmartRef<Minerva::Prong>& protonProng, 
                                                double& p_calCorrection, 
                                                double& p_visEnergyCorrection ) const
{
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    debug() << "Enter CCDeltaPlusAna::correctProtonProngEnergy()" << endmsg;
    
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
            
    debug() << "Exit CCDeltaPlusAna::correctProtonProngEnergy()" << endmsg;
    debug() <<"--------------------------------------------------------------------------"<<endmsg;
    return;
}

    







 
