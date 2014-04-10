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
    Last Revision:  2014_04_10
    
================================================================================
*/
#ifndef CCDELTAPLUSANA_H 
#define CCDELTAPLUSANA_H 1

#include <utility>

// ineritance
#include "AnaUtils/MinervaAnalysisTool.h"

//-- Forward Declarations
#include "Event/MinervaEventFwd.h"


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
        StatusCode interpretEvent( const Minerva::PhysicsEvent* event, const Minerva::GenMinInteraction* truth, NeutrinoVect& nuInts ) const;
        
        StatusCode tagTruth( Minerva::GenMinInteraction* truth ) const;
        
    private:
        double m_fidHexApothem;
        double m_fidUpStreamZ;
        double m_fidDownStreamZ;
  
  
        StatusCode interpretFailEvent( Minerva::PhysicsEvent* event ) const;
  
  
};

#endif

