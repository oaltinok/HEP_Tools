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
    Last Revision:  2014_04_22
    
================================================================================
*/
#ifndef CCDELTAPLUSANA_H 
#define CCDELTAPLUSANA_H 1

#include <utility>

// ineritance
#include "AnaUtils/MinervaAnalysisTool.h"

//-- Forward Declarations
#include "Event/MinervaEventFwd.h"

class TRandom3;
class IMichelTool;
class IMinervaCoordSysTool;
class IProtonUtils;
class IEnergyCorrectionTool;
class IHitTaggerTool;
class IProngClassificationTool;
class IODProngClassificationTool;
class IParticleMakerTool;
class ICCPionIncUtils;

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
        StatusCode reconstructEvent( Minerva::PhysicsEvent* event, Minerva::GenMinInteraction* truthEvent = NULL ) const;
        
        //! Attach an interpretations to the event (mandatory for inheritance)
        StatusCode interpretEvent( const Minerva::PhysicsEvent* event, const Minerva::GenMinInteraction* truthEvent, NeutrinoVect& interaction_hyp ) const;
        
        StatusCode tagTruth( Minerva::GenMinInteraction* truthEvent ) const;
        
    private:
        // Fiducial Volume
        double m_fidHexApothem;
        double m_fidUpStreamZ;
        double m_fidDownStreamZ;
        
        // Analysable Volume
        double m_recoHexApothem;
        double m_recoUpStreamZ;
        double m_recoDownStreamZ;
        
        double m_beamAngleBias;
        
        double m_minMuonScore;
        double m_minProtonScore;
        
        bool m_store_all_events;
        bool m_makeShortTracks;
        bool m_doPlausibilityCuts;
        
        int m_muonProngColor; 
        int m_protonProngColor; 
        int m_primaryVertexColor; 
        int m_secondaryVertexColor; 
        int m_endPointVertexColor; 
        int m_unattachedProngColor;
        
        TRandom3*                 m_randomGen;
        unsigned long int         m_randomSeed;
        
        Minerva::DeDetector*        m_InnerDetector;
        
        IMinervaCoordSysTool*       m_coordSysTool;
        IMichelTool*                m_michelTrkTool;
        IMichelTool*                m_michelVtxTool;
        IProtonUtils*               m_protonUtils;
        IEnergyCorrectionTool*      m_energyCorrectionTool;
        IHitTaggerTool*             m_hitTagger;
        IProngClassificationTool*   m_prongIntersection;
        IODProngClassificationTool* m_odMatchTool;
        IParticleMakerTool*         m_particleMaker;
        ICCPionIncUtils*            m_ccPionIncUtils;
        
        //! Private Functions
        StatusCode interpretFailEvent( Minerva::PhysicsEvent* event ) const;
        StatusCode getNearestPlane(double z, int & module_return, int & plane_return) const;        
        bool createTrackedParticles( Minerva::ProngVect& prongs ) const;
        
        bool getProtonProng(    Minerva::ProngVect& primaryProngs, 
                                Minerva::ProngVect& hadronProngs,
                                Minerva::ParticleVect& hadrons ) const;
                                                
        void setProtonParticleData( Minerva::NeutrinoInt* nuInt, 
                                    Minerva::ProngVect& protonProngs,
                                    Minerva::ParticleVect& protonParticles, 
                                    double vertexZ ) const;
                                                    
        void correctProtonProngEnergy(  SmartRef<Minerva::Prong>& protonProng, 
                                        double& p_calCorrection, 
                                        double& p_visEnergyCorrection ) const;
                                        

        
        
  
  
  
};

#endif

