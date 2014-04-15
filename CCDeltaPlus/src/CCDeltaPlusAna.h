/*
================================================================================
CCDeltaPlusAna

    Reconstruction Package:
        Exclusive Channel for muon,proton,pi0 on the final state
        Uses Nightly Build
    
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
    Date:           2014_03_27
    Last Revision:  2014_04_15
    
================================================================================
*/
#ifndef CCDELTAPLUSANA_H 
#define CCDELTAPLUSANA_H 1

#include <utility>

// ineritance
#include "AnaUtils/MinervaAnalysisTool.h"

//-- Forward Declarations
#include "Event/MinervaEventFwd.h"

class IMichelTool;
class IMinervaCoordSysTool;

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
        
        StatusCode initialize();
        StatusCode finalize();
        
        //! Reconstruct the event (mandatory for inheritance)
        StatusCode reconstructEvent( Minerva::PhysicsEvent* event, Minerva::GenMinInteraction* truth = NULL ) const;
        
        //! Attach an interpretations to the event (mandatory for inheritance)
        StatusCode interpretEvent( const Minerva::PhysicsEvent* event, const Minerva::GenMinInteraction* truth, NeutrinoVect& interaction_hyp ) const;
        
        StatusCode tagTruth( Minerva::GenMinInteraction* truth ) const;
        
    private:
        double m_fidHexApothem;
        double m_fidUpStreamZ;
        double m_fidDownStreamZ;
        
        double m_beamAngleBias;
        
        double m_minMuonScore;
        
        bool m_store_all_events;
        bool m_makeShortTracks;
        bool m_doPlausibilityCuts;
        
        Minerva::DeDetector*        m_InnerDetector;
        
        IMinervaCoordSysTool*       m_coordSysTool;
        IMichelTool*                m_michelTrkTool;
        IMichelTool*                m_michelVtxTool;
        
        
        //! Private Functions
        StatusCode interpretFailEvent( Minerva::PhysicsEvent* event ) const;
        StatusCode getNearestPlane(double z, int & module_return, int & plane_return) const;
        
        
  
  
  
};

#endif

