
#include "CCDeltaPlusAna.h"

#include "RecInterfaces/IFiducialPointTool.h"


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
    

    // Initializing Analysis Tools
    
    m_fidHexApothem  = 850.0*CLHEP::mm;
    m_fidUpStreamZ   = 5990.0*CLHEP::mm;   // ~middle of module 27, plane 1
    m_fidDownStreamZ = 8340.0*CLHEP::mm;   // ~middle of module 79, plane 1
    
    

    
    //! declare common branches
    declareCommonPhysicsAnaBranches();
    declareMinosMuonBranches();
    declareGenieWeightBranches();
    

    // Select the branches you want in your AnaTuple


    
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
    
  
//     //-- check if this a plausible event ( MC only )
//     if( truth && !truthIsPlausible(truth) ) {
//         info() << "   This is not a plausible MC event! returning!" << endmsg;
//         return interpretFailEvent(event);
//     }
    
//     info() << "Vertex Cuts" << endmsg;
//     //-- interaction vertex check
//     if( !event->hasInteractionVertex() ) {
//         info() << "The event does not have an interaction vertex!" << endmsg;
//         event->setIntData("NoInteractionVertex",1);
//         return interpretFailEvent(event);
//     } 
    
//     //-- get the interaction vertex 
//     SmartRef<Minerva::Vertex> vertex = event->interactionVertex();
//     if( !vertex ) { 
//         bool pass = true; 
//         std::string tag = "BadObject";
//         event->filtertaglist()->addFilterTag(tag,pass);
//         event->setIntData("NullVertex",1);
//         error() << "This vertex is NULL! Flag this event as bad!" << endmsg;
//         return interpretFailEvent(event);
//     }
    
//     // is Vertex Fiducial?
//     Gaudi::XYZPoint position = vertex->position();
//     if( !FiducialPointTool->isFiducial( position, m_fiducialApothem, m_fiducialUpstreamZ, m_fiducialDownstreamZ ) ){
//         info() << " Interaction Vertex is not fiducial!" << endmsg;
//         event->setIntData("survive_fiducial",1);
//         return interpretFailEvent(event);
//     }
    
//     info() << "Muon Cuts" << endmsg;
//     // - Muon exist
//     bool has_muon = MuonUtils->findMuonProng( event, m_MuonProng, m_MuonParticle );
//     if( !has_muon ){
//         info() << "Did not find a muon prong! This cannot be a CCDeltaPlus event." << endmsg;
//         return interpretFailEvent(event);
//     }


    Minerva::NeutrinoInt *nuInt = new Minerva::NeutrinoInt( "CCDeltaPlusAna" );
    NeutrinoVect defaultHyp;
    defaultHyp.push_back(nuInt);
    addInteractionHyp(event,defaultHyp);
  
    info() <<"Exit CCDeltaPlusAna::reconstructEvent()" << endmsg;
    return StatusCode::SUCCESS;
    
}
    
//=============================================================================
// interpretEvent() --
//=============================================================================
StatusCode CCDeltaPlusAna::interpretEvent( const Minerva::PhysicsEvent *event, const Minerva::GenMinInteraction *truth, NeutrinoVect& nuInts ) const
{
    info() <<"Enter: CCDeltaPlusAna::interpretEvent()" << endmsg;
    
    if( truth ){
        info() << "\tThis is a MC event." << endmsg;
    }
    
    if( !event ){
        info() << "\tNULL Event" << endmsg;
        return StatusCode::FAILURE;
    }
  //--------------------------------------------------------------
  // create interaction hypothesis for CCDeltaPlus
  //--------------------------------------------------------------
   Minerva::NeutrinoInt* ccdeltaplus = new Minerva::NeutrinoInt(m_anaSignature);
   ccdeltaplus ->setNeutrinoFlavor(Minerva::NeutrinoInt::MuonFlavor);
   ccdeltaplus ->setInteractionCurrent(Minerva::NeutrinoInt::ChargedCurrent);
  
  
   //-- store ccdeltaplus hypothesis
   nuInts.push_back( ccdeltaplus  );
 

    info() <<"Exit CCDeltaPlusAna::interpretEvent()" << endmsg;
    return StatusCode::SUCCESS;
    
}
   

//=============================================================================
// tagTruth 
//=============================================================================
StatusCode CCDeltaPlusAna::tagTruth( Minerva::GenMinInteraction* truth ) const 
{
    
    info() << "Enter: tagTruth()" << endmsg;
    
    
    info() <<"Do Nothing!" << endmsg;
    
    if( !truth ) {
        warning() << "The GenMinInteraction is NULL!" << endmsg;
        return StatusCode::SUCCESS;
    }
    
    info() <<"Exit CCDeltaPlusAna::tagTruth()" << endmsg;
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
    

    







 
