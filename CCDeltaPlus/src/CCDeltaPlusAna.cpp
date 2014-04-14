/*
    See CCDeltaPlusAna.h header for Class Information
*/
#include "CCDeltaPlusAna.h"

#include "RecInterfaces/IFiducialPointTool.h"

#include "GeoUtils/IMinervaCoordSysTool.h"
#include "ProngMaker/IMichelTool.h"

#include "MinervaDet/IGeomUtilSvc.h"
#include "MinervaDet/DeDetector.h"



DECLARE_TOOL_FACTORY( CCDeltaPlusAna );

using namespace Minerva;
    
//=============================================================================
// Standard constructor 
//=============================================================================
CCDeltaPlusAna::CCDeltaPlusAna(const std::string& type, const std::string& name, const IInterface* parent ) :
MinervaAnalysisTool( type, name, parent ) 
{
    
    info() <<"Enter CCDeltaPlusAna::CCDeltaPlusAna() -- Default Constructor" << endmsg;
    
    declareInterface<IInteractionHypothesis>(this);
    //! mandatory declaration of analysis signature: CCDeltaPlusAna
    m_anaSignature = "CCDeltaPlusAna";
    
    // Private Properties
    declareProperty( "BeamAngleBias", m_beamAngleBias = 0.006*CLHEP::radian );
    
    declareProperty( "MinMuonScore", m_minMuonScore = 0.9 );
    
    declareProperty( "DoPlausibilityCuts", m_doPlausibilityCuts = true );
    
    declareProperty( "MakeShortTracks", m_makeShortTracks = true );
    
    
    
    // Protected properties from IInteractionHypothesis.
    m_hypMeths.push_back( m_anaSignature );
    declareProperty("HypothesisMethods", m_hypMeths);
    info() << " CCDeltaPlusAna Hypothesis added " << endmsg;
    
    info() <<"Exit CCDeltaPlusAna::CCDeltaPlusAna() -- Default Constructor" << endmsg;
    
}
    
//=============================================================================
// Initialize
//=============================================================================
StatusCode CCDeltaPlusAna::initialize()
{
    
    info() <<"Enter CCDeltaPlusAna::initialize()" << endmsg;
    
    //! Initialize the base class.  This will fail if you did not define m_anaSignature.
    StatusCode sc = this->MinervaAnalysisTool::initialize();
    if( sc.isFailure() ) { 
        return Error( "Failed to initialize!", sc ); 
    
    }
    info()<<"   initialized MinervaAnalysisTool"<<endmsg;
    

    m_fidHexApothem  = 850.0*CLHEP::mm;
    m_fidUpStreamZ   = 5990.0*CLHEP::mm;   // ~middle of module 27, plane 1
    m_fidDownStreamZ = 8340.0*CLHEP::mm;   // ~middle of module 79, plane 1
    
    
    //! Initializing Analysis Tools
    try {
        m_InnerDetector = getDet<Minerva::DeDetector>("/dd/Structure/Minerva/Detector/InnerDetector");
    } catch(GaudiException& e) {
        error() << "Could not obtain Inner Detector" << endmsg;
        return StatusCode::FAILURE;
    }
    
    try{
        m_coordSysTool = tool<IMinervaCoordSysTool>("MinervaCoordSysTool");
    }catch(GaudiException& e){
        error() << "Could not obtain MinervaCoordSysTool" << endmsg;
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
    

    
    //! declare common branches
    declareCommonPhysicsAnaBranches();
    declareMinosMuonBranches();
    declareGenieWeightBranches();
    

    // Select the branches you want in your AnaTuple
    
    
    
    //Event - Cut Results
    declareIntEventBranch( "Cut_NoInteractionVertex", -1 );
    declareIntEventBranch( "Cut_NullVertex", -1 );
    declareIntEventBranch( "Cut_FiducialVertex", -1 );
    declareIntEventBranch( "is_CCDeltaPlus", -1 );
    
    //Event - general reco
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
  
    //Event - Michel 
    declareIntEventBranch( "has_vtx_michel", 0 );
    declareIntEventBranch( "n_vtx_michel_views", 0 );
    declareDoubleEventBranch( "vtx_michel_distance", -1.0);
    
    //Event - energy
    declareDoubleEventBranch( "muonVisibleE", -1.0 );
    declareDoubleEventBranch( "hadronVisibleE", -1.0 );
    declareDoubleEventBranch( "totalVisibleE", -1.0 );
    declareDoubleEventBranch( "totalIDVisibleE", -1.0 );
    declareDoubleEventBranch( "totalODVisibleE", -1.0 );
    declareDoubleEventBranch( "vtxBlobExtraE",  -1.0 );
    declareDoubleEventBranch( "unattachedExtraE", -1.0 );
    declareDoubleEventBranch( "dispersedExtraE", -1.0 ); 
    
    
    //! NeutrinoInt - vtx
    declareIntBranch( m_hypMeths, "vtx_module", -99);
    declareIntBranch( m_hypMeths, "vtx_plane",-1);
    declareDoubleBranch( m_hypMeths, "vtx_x",0.0);
    declareDoubleBranch( m_hypMeths, "vtx_y",0.0);
    declareDoubleBranch( m_hypMeths, "vtx_z",0.0);
    
    //! NeutrinoInt - muon
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
    
    
    info() <<"Exit CCDeltaPlusAna::initialize()" << endmsg;
    
    return sc;
}
    
//=============================================================================
// reconstructEvent() --
//=============================================================================
StatusCode CCDeltaPlusAna::reconstructEvent( Minerva::PhysicsEvent *event, Minerva::GenMinInteraction* truth ) const
{
    
    debug() <<"Enter CCDeltaPlusAna::reconstructEvent()" << endmsg;
    
    if( truth ){
        debug() << "\tThis is a MC event." << endmsg;
    }
    
    //--------------------------------------------------------------------------
    //! Check if this a plausible event ( MC only )
    //--------------------------------------------------------------------------
    if( truth && !truthIsPlausible(truth) ) {
        debug() << "   This is not a plausible MC event! returning!" << endmsg;
        return StatusCode::SUCCESS;
    }
    
    
    //--------------------------------------------------------------------------
    //! Get the interaction vertex, if it exists
    //--------------------------------------------------------------------------
    debug() << "Finding Vertex..." << endmsg;
    
    if( !event->hasInteractionVertex() ) {
        debug() << "The event does not have an interaction vertex!" << endmsg;
        event->setIntData("Cut_NoInteractionVertex",1);
        return StatusCode::SUCCESS;
    } 
       
       
    SmartRef<Minerva::Vertex> interactionVertex = event->interactionVertex();
    if( !interactionVertex ) { 
        bool pass = true; 
        std::string tag = "BadObject";
        event->filtertaglist()->addFilterTag(tag,pass);
        event->setIntData("Cut_NullVertex",1);
        error() << "This vertex is NULL! Flag this event as bad!" << endmsg;
        return StatusCode::SUCCESS;
    }
    
    // is Vertex Fiducial?
    Gaudi::XYZPoint vtx_position = interactionVertex->position();
    if( !FiducialPointTool->isFiducial( vtx_position, m_fidHexApothem, m_fidUpStreamZ, m_fidDownStreamZ ) ){
        debug() <<"Interaction Vertex is not in fiducial volume = ("<<vtx_position.x()<<","<<vtx_position.y()<<","<<vtx_position.z()<<")"<< endmsg;
        event->setIntData("Cut_FiducialVertex",1);
        return StatusCode::SUCCESS;
    }
    
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
            return StatusCode::SUCCESS;
        }
        
        debug() << "Muon Particle Score: " << muonPart->score() << endmsg;
        if (muonPart->score() >= m_minMuonScore) {
        
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
            if (truth) truth->filtertaglist()->setOrAddFilterTag( "reco_isMinosMatch", true );
            
            //! @todo - also determine whether there is a minos track in event, if there is no minos match track or stub
        }
        else {
            debug()<<"Muon prong does not pass score cut"<<endmsg;
            return StatusCode::SUCCESS;
        }
        
    } 
    else{
        debug() << "Did not find a muon prong!" << endmsg;
        return StatusCode::SUCCESS;
    } // end if findMuonProng
    
    double muon_visible_energy = muonProng->minervaVisibleEnergySum();
     
    //--------------------------------------------------------------------------
    //! Determine if vertex has michels
    //--------------------------------------------------------------------------
    debug() << "Finding Michels..." << endmsg;
    Minerva::Prong vtx_michel_prong;
    bool foundMichel = m_michelVtxTool->findMichel( interactionVertex, vtx_michel_prong );
    if (foundMichel) {
        debug()<<"Found a Michel Electron!"<<endmsg;
        event->setIntData("has_vtx_michel",1);
        event->setIntData("n_vtx_michel_views",vtx_michel_prong.getIntData("category"));
        event->setDoubleData("vtx_michel_distance",vtx_michel_prong.getDoubleData("distance"));
    }else{
        debug()<<"There are NO Michel Electrons in the event!"<<endmsg;
    }
    
    
    
    //--------------------------------------------------------------
    //! Finish filling event portion of ntuple 
    //--------------------------------------------------------------
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
//     event->setIntData("n_anchored_long_trk_prongs", n_anchored_long_trk_prongs);
//     event->setIntData("n_anchored_short_trk_prongs", n_anchored_short_trk_prongs);
//     event->setIntData("n_vtx_blob_prongs", n_vtx_blob_prongs);
//     event->setIntData("n_iso_trk_prongs", n_iso_trk_prongs);
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

    //--------------------------------------------------------------
    //! Call the interpretEvent function.
    //--------------------------------------------------------------
    NeutrinoVect interactions;
    StatusCode interpret = this->interpretEvent( event, truth, interactions );
    
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
    return sc;
    
}
    
//=============================================================================
// interpretEvent() --
//=============================================================================
StatusCode CCDeltaPlusAna::interpretEvent( const Minerva::PhysicsEvent *event, const Minerva::GenMinInteraction *truth, NeutrinoVect& interaction_hyp ) const
{
    debug() <<"Enter: CCDeltaPlusAna::interpretEvent()" << endmsg;
    
    if( truth ){
        debug() << "\tThis is a MC event." << endmsg;
    }
    
    if( !event ){
        debug() << "\tNULL Event" << endmsg;
        return StatusCode::FAILURE;
    }
    
    //--------------------------------------------------------------
    //! Create interaction hypothesis
    //--------------------------------------------------------------
    Minerva::NeutrinoInt* nuInt = new Minerva::NeutrinoInt(m_anaSignature);
    interaction_hyp.push_back( nuInt );

    //--------------------------------------------------------------
    //! Get the Primary Muon Prongs
    //--------------------------------------------------------------
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
    
    
    //--------------------------------------------------------------  
    //! Get the muon particle hypotheses and fill some info
    //--------------------------------------------------------------
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
        
    debug()<<"Finished muons"<<endmsg;
    
    //--------------------------------------------------------------
    //! Fill ntuple
    //--------------------------------------------------------------
    
    //nuInt->setNeutrinoHelicity( getHelicity( mu_charge ) );
    nuInt->setInteractionCurrent( Minerva::NeutrinoInt::ChargedCurrent );
    nuInt->setInteractionType( Minerva::NeutrinoInt::UnknownInt );
    Gaudi::XYZTVector position( event->interactionVertex()->position().x(), 
                                event->interactionVertex()->position().y(),
                                event->interactionVertex()->position().z(),
                                event->time() );     
    nuInt->setVertex( position );
    nuInt->setScore( 1.0 );
    nuInt->setLeptonEnergy( muon_4p );
    
    int vtx_module, vtx_plane;
    debug()<<"Calling getNearestPlane, vtx is "<<position.z()<<endmsg;
    getNearestPlane(position.z(), vtx_module, vtx_plane); 
    nuInt->setIntData("vtx_module", vtx_module);
    nuInt->setIntData("vtx_plane", vtx_plane);
    nuInt->setDoubleData("vtx_x", position.x() );
    nuInt->setDoubleData("vtx_y", position.y() );
    nuInt->setDoubleData("vtx_z", position.z() );
    
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

    debug() <<"Exit CCDeltaPlusAna::interpretEvent()" << endmsg;
    return StatusCode::SUCCESS;
    
}
   

//=============================================================================
// tagTruth 
//=============================================================================
StatusCode CCDeltaPlusAna::tagTruth( Minerva::GenMinInteraction* truth ) const 
{
    
    debug() << "Enter: tagTruth()" << endmsg;

    
    if( !truth ) {
        warning() << "The GenMinInteraction is NULL!" << endmsg;
        return StatusCode::SUCCESS;
    }
    
    debug() << "Filling Genie Weight Branches" << endmsg;
    //--------------------------------------------------------------
    //! Fill GENIE Weight Branches
    //--------------------------------------------------------------      
    StatusCode sc = fillGenieWeightBranches( truth );
    if (sc.isFailure() ) {
        warning()<<"Genie weight branch filling failed!"<<endmsg;
        return sc;
    }
    
    debug() <<"Exit CCDeltaPlusAna::tagTruth()" << endmsg;
    return StatusCode::SUCCESS;
    
}

//=======================================================================
//  Finalize
//=======================================================================
StatusCode CCDeltaPlusAna::finalize()
{
    
    info() << "Enter: finalize()" << endmsg;
    
    
    // finalize the base class.
    StatusCode sc = this->MinervaAnalysisTool::finalize();
    if( sc.isFailure() ) return Error( "Failed to finalize!", sc );
    
    
    info() <<"Exit CCDeltaPlusAna::finalize()" << endmsg;
    return sc;
}

//==============================================================================
//  Private Functions
//==============================================================================
    
    
//=============================================================================
// Find the plane nearest to a point
//=============================================================================
StatusCode CCDeltaPlusAna::getNearestPlane(double z, int & module_return, int & plane_return) const
{
  // Got From Brandon CCNuPionInc - 2014_04_14
  //testing new MinervaDet routine
  //works as advertised, but a couple problems:
  //1) Does not go to more downstream plane (upstream if backwards track) if point is in passive material
  DePlane const * pPlane = m_InnerDetector->getClosestDePlane(z);
  if (!pPlane) {
    error()<<"        getClosestDePlane failed on z "<<z<<endmsg;
    return StatusCode::FAILURE;
  }
  
  PlaneID planeid = pPlane->getPlaneID();
  module_return = planeid.module();
  plane_return = planeid.plane();
  
  return StatusCode::SUCCESS;
  
}

    







 
