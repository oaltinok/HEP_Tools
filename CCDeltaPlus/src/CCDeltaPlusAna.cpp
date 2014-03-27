#undef NDEBUG
#include <cassert>
#include <algorithm>
#include <iostream>
#include <iomanip>

#include <TString.h>
#include <TDatabasePDG.h>
#include <TVector3.h>
#include <TMath.h>


#include <Event/DAQHeader.h>
#include <Event/GenMinHeader.h>
#include <Event/MCIDDigit.h>
#include <Event/TG4Trajectory.h>
#include <Event/IDCluster.h>
#include <Event/VectorTypeDefs.h>
#include <Event/Vertex.h>
#include <Event/Track.h>

#include <GeoUtils/IMinervaCoordSysTool.h>
#include <MinervaUtils/IMinervaMathTool.h>
#include <MinervaUtils/IHitTaggerTool.h>


#include <MinervaDet/IGeomUtilSvc.h>
#include <MinervaDet/DeDetector.h>
#include <ODDet/DeOuterDetector.h>
#include <MinervaDet/DeSubdet.h>

#include <EnergyRecTools/IExtraEnergyTool.h>
#include <RecUtils/IClusterUtilsTool.h>
#include <RecInterfaces/IIDBlobCreator.h>
#include <RecInterfaces/IODBlobCreator.h>
#include <RecInterfaces/ITrackLinearPropagator.h>
#include <RecInterfaces/IFiducialPointTool.h>
#include <BlobFormation/IBlobCreatorUtils.h>
#include <BlobFormation/IIDIsolatedIDBlobCreator.h>
#include <BlobFormation/IIDBlobSeedingTool.h>
#include <BlobFormation/IIDAnchoredBlobCreator.h>

#include <AnaUtils/IMuonUtils.h>
#include <AnaUtils/IPhysicsCalculator.h>

#include <CCDeltaPlus/IHoughBlob.h>
#include <CCDeltaPlus/IHoughTool.h>

#include "CCDeltaPlusAna.h"
#include "TraverseHistory.h"
#include "TrackTruthInfo.h"
#include "DigitVectorTruthInfo.h"
#include "HitVectorTruthInfo.h"
#include "OneParLineFit.h"
#include "Par.h"
#include "AngleScan.h"
#include "ClusterVectorInfo.h"


/* Global variables, instead of class data members so that they can be
   assigned in 'const' methods. In general, we can use 'mutable' members,
   but it does not seem to work in the framework */
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
  debug() << " CCDeltaPlusAna constructor <oaltinok_version>" << endmsg;
  declareInterface<IInteractionHypothesis>(this);

	//! mandatory declaration of analysis signature: CCDeltaPlusReco
  m_anaSignature = "CCDeltaPlusAna";

	// Protected properties from IInteractionHypothesis.
  m_hypMeths.push_back( m_anaSignature );
  declareProperty("HypothesisMethods", m_hypMeths);

  debug() << " CCDeltaPlusAna Hypothesis added <oaltinok_version>" << endmsg;

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
  
}

//=============================================================================
// Initialize
//=============================================================================
StatusCode CCDeltaPlusAna::initialize()
{
  debug() << "CCDeltaPlusAna::initialize()<oaltinok_version>" << endmsg;

      //! Initialize the base class.  This will fail if you did not define m_anaSignature.
  StatusCode sc = this->MinervaAnalysisTool::initialize();
  if( sc.isFailure() ) 
          return Error( "Failed to initialize!", sc );
  
  try{
      m_extraEnergyTool = tool<IExtraEnergyTool>("ExtraEnergyTool");
  }
  catch(GaudiException& e){
      error()<<"Could not obtain tool: ExtraEnergyTool<oaltinok_version>" << endmsg;
      return StatusCode::FAILURE;
  }

  try {
      m_blobUtils = tool<IBlobCreatorUtils>("BlobCreatorUtils");
  } 
  catch( GaudiException& e ){
      error() << "Could not obtain tool: BlobCreatorUtils! <oaltinok_version>" << endmsg;
      return StatusCode::FAILURE;
  }
  
  try{
      m_idConeScanBlob = tool<IIDAnchoredBlobCreator>("ConeScanIDBlobCreator");
  } 
  catch (GaudiException& e){
      error() << " Could not obtain tool: ConeScanIDBlobCreator <oaltinok_version>" << endmsg;
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
      m_colorTag = tool<IHitTaggerTool>( "HitTaggerTool" );
  }
  catch( GaudiException& e ) {
      error() << "Could not obtain tool: HitTaggerTool<oaltinok_version>" << endmsg;
      return StatusCode::FAILURE;
  }
  
  try{
      m_idHoughBlob = tool<IHoughBlob>("HTBlob");
  }
  catch (GaudiException& e){
      error() << " Could not obtain tool: HBlob <oaltinok_version>" << endmsg;
      return StatusCode::FAILURE;
  }

  try{
      m_idHoughTool = tool<IHoughTool>("HTtool");
  }
  catch (GaudiException& e){
      error() << " Could not obtain tool: HTtool <oaltinok_version>" << endmsg;
      return StatusCode::FAILURE;
  }
  
  try {
      m_blobSeedingTool = tool<IIDBlobSeedingTool>("BlobSeedingTool"); 
  }
  catch( GaudiException& e ) { 
      error() << "Could not obtain tool: BlobSeedingTool !<oaltinok_version>" << endmsg;
      return StatusCode::FAILURE;
  }
  
  try {
      m_minervaCoordSysTool = tool<IMinervaCoordSysTool>("MinervaCoordSysTool");
  } catch( GaudiException& e ) {
      error() << "Could not obtain tool: MinervaCoordSysTool<oaltinok_version>" << endmsg;
      return StatusCode::FAILURE;
  }

  try {
      m_mathTool = tool<IMinervaMathTool>("MinervaMathTool");
  }
  catch( GaudiException& e ) {
      error() << "Could not obtain tool: MinervaMathTool!<oaltinok_version>" << endmsg;
      return StatusCode::FAILURE;
  }
 
  service("GeomUtilSvc", m_GeomUtilSvc, true);
  m_idDet = m_GeomUtilSvc->getIDDet();
  m_odDet = m_GeomUtilSvc->getODDet();
  
  Constants().k_hcal = 0.088;
  Constants().k_ecal = 0.360;
  Constants().k_trkr = 0.760;
  Constants().coangle = 15.0;
  Constants().clength = 3000.0;
  Constants().xminevis = 20.0;
  Constants().uminevis = 10.0;
  Constants().vminevis = 10.0;
  
  declareCommonPhysicsAnaBranches();
      //declareGenieWeightBranches();

  declareDoubleTruthBranch( "MC_pi0_energy", -9999 );
  declareContainerDoubleTruthBranch( "MC_pi0_momentum" );
  
  declareDoubleTruthBranch( "MC_scalar", -9999);
  declareDoubleTruthBranch( "MC_photon_energy_1", -9999 );
  declareDoubleTruthBranch( "MC_photon_energy_2", -9999 );
  declareContainerDoubleTruthBranch( "MC_photon_direction_1");
  declareContainerDoubleTruthBranch( "MC_photon_direction_2");
  declareContainerDoubleTruthBranch( "MC_photon_vertex_1" );
  declareContainerDoubleTruthBranch( "MC_photon_vertex_2" );


      // "Event" info
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

  return sc;
}

//=============================================================================
// reconstructEvent() --
//=============================================================================
StatusCode CCDeltaPlusAna::reconstructEvent( Minerva::PhysicsEvent *event, Minerva::GenMinInteraction* truth ) const
{
    
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

    
  debug() << "CCDeltaPlusAna::reconstructEvent( PhysicsEvent *event, const GenMinInteraction* truth )<oaltinok_version>" << endmsg;
  debug() << "=============================================================================<oaltinok_version>" << endmsg;
  debug() << "                                                                             <oaltinok_version>" << endmsg;
  
  debug() << gateData() << endmsg;

  markEvent(event);


  Minerva::NeutrinoInt *nuInt = new Minerva::NeutrinoInt( "CCDeltaPlusAna" );
  std::vector<Minerva::NeutrinoInt*> defaultHyp;
  defaultHyp.push_back(nuInt);
  addInteractionHyp(event,defaultHyp);

  
  SmartRef<Minerva::Vertex> vertex;//(new Minerva::Vertex(Gaudi::XYZPoint(0,0,0)));
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

  fillCommonPhysicsAnaBranches( event );

      
      //markEvent( event );
  std::cout << "Passing all cuts " << std::endl;
  // Now interpret the event and add NeutrinoInts
  std::vector<Minerva::NeutrinoInt*> interactions;
  interpretEvent( event, truth, interactions );

      // Add the newly create NeutrinoInts to this PhysicsEvent
  StatusCode sc = addInteractionHyp( event, interactions );
  
  return StatusCode::SUCCESS;

}

//=============================================================================
// interpretEvent() --
//=============================================================================
StatusCode CCDeltaPlusAna::interpretEvent( const Minerva::PhysicsEvent *event, const Minerva::GenMinInteraction *truth,
                                         std::vector<Minerva::NeutrinoInt*>& interactions ) const
{
    debug() << " == CCDeltaPlusAna::interpretEvent ==<oaltinok_version>" << endmsg;
    if( event ) debug() << " working on Interpret Event and pass all cuts <oaltinok_version>" << endmsg;
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
    interactions.push_back( nuInt );
    
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
        //info() << " Jumping rock Muon <oaltinok_version>" << endmsg;
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
        debug() << "Did not find a muon prong! This cannot be a CCDeltaPlus event.<oaltinok_version>" << endmsg;
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
        debug() << " Interaction Vertex is not fiducial!<oaltinok_version>" << endmsg;
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
    debug() << " CCDeltaPlusAna::DoVertex <oaltinok_version>" << endmsg;

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
        debug() << " Discarding " <<  farTracks.size() << " tracks, to get photons <oaltinok_version>" << endmsg;
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

  debug() << " CCDeltaPlusAna::VtxBlob <oaltinok_version>" << endmsg;

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

    debug() << " CCDeltaPlusAna::ConeBlobs <oaltinok_version>" << endmsg;

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
    debug() << " CCDeltaPlusAna::AngleScanBlob <oaltinok_version>" << endmsg;

    const Gaudi::XYZPoint& pos = vertex->position();
    
        // main interface to create blob, this tool already exclude lowactivity clusters.
    std::vector<Minerva::IDBlob*>* preidBlobs = new std::vector<Minerva::IDBlob*>;
    m_idConeScanBlob->createIDBlobs(idClusters,preidBlobs,pos);
    
    std::vector<Minerva::IDBlob*>::iterator itBlob = preidBlobs->begin();
    for ( ; itBlob != preidBlobs->end(); itBlob ++ ){
        if ( !m_idHoughBlob->isPhoton( (*itBlob)->clusters(), pos ) ) continue;
        outBlobs.push_back(*itBlob);
    }

    info() << " Angle Scan is done <oaltinok_version>" << endmsg;
    
    return StatusCode::SUCCESS;
}

//=======================================================================
//  HoughBlob
//=======================================================================
StatusCode CCDeltaPlusAna::HoughBlob(SmartRefVector<Minerva::IDCluster> idClusters,
                                   const SmartRef<Minerva::Vertex>& vertex,
                                   std::vector<Minerva::IDBlob*>& outBlobs) const
{
    debug() << " CCDeltaPlusAna::HoughBlob <oaltinok_version>" << endmsg;

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
        debug() << " Pass Get Create2dHT <oaltinok_version>" << endmsg;	
        
        Gaudi::XYZPoint startpoint(spX, 0, spZ);
        Gaudi::XYZVector direction(-1/tan(theta*CLHEP::pi/180),0,1);

        if ( !m_idHoughBlob->PseudoCone(idHTSeed, idClusViewX, direction, pos ) ) continue;
        debug() << " Pass PseudoCone <oaltinok_version>" << endmsg;
        
        if ( !m_idHoughBlob->XUVMatch(idHTSeed, idClusViewU, idClusViewV, 50 ) ) continue;
        debug() << " Pass XUV Match - First <oaltinok_version>" << endmsg;	
        
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
        debug() << *itBlob << " " << Uclusters << " U clusters; " << Vclusters << " V clusters <oaltinok_version>" << endmsg;
        (*itBlob)->clear();
        if ( Uclusters != .0 || Vclusters != .0 ) {
            debug() << *itBlob << " inserting to FinalBlobs <oaltinok_version>" << endmsg;
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


    info() << " Hough Transform is done! <oaltinok_version>" << endmsg; 

    
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
  debug() << " tagTruth <oaltinok_version>" << endmsg;
    
  fillGenieWeightBranches( truth );
    
        //! Check to see if the MC event is plausibly analyzable - important for Data overlay.
  bool isPlausible = truthIsPlausible( truth );
  truth->filtertaglist()->setOrAddFilterTag( "pass_plausible", isPlausible );
    
	// is interaction in fiducial volume?
  bool is_fiducial = FiducialPointTool->isFiducial( truth, m_fiducialApothem, m_fiducialUpstreamZ, m_fiducialDownstreamZ );
  truth->filtertaglist()->setOrAddFilterTag( "is_fiducial", is_fiducial );
  
  debug() << " Pass fiducial point " << is_fiducial<< endmsg;

      // get primary final state lepton stuff
  double fslepton_P = truth->PrimFSLepton().P(); 
  double fslepton_E = truth->PrimFSLepton().E();
  double fslepton_T = fslepton_E - MinervaUnits::M_mu;
  double fslepton_theta = m_minervaCoordSysTool->thetaWRTBeam( truth->PrimFSLepton() );
  double fslepton_theta_x = m_minervaCoordSysTool->thetaXWRTBeam( truth->PrimFSLepton() );
  double fslepton_theta_y = m_minervaCoordSysTool->thetaYWRTBeam( truth->PrimFSLepton() );
  double fslepton_phi = m_minervaCoordSysTool->phiWRTBeam( truth->PrimFSLepton() );

  truth->setDoubleData( "fslepton_E", fslepton_E*CLHEP::MeV/CLHEP::GeV );
  truth->setDoubleData( "fslepton_P", fslepton_P*CLHEP::MeV/CLHEP::GeV );
  truth->setDoubleData( "fslepton_T", fslepton_T*CLHEP::MeV/CLHEP::GeV );
  truth->setDoubleData( "fslepton_theta", fslepton_theta );
  truth->setDoubleData( "fslepton_theta_x", fslepton_theta_x );
  truth->setDoubleData( "fslepton_theta_y", fslepton_theta_y );
  truth->setDoubleData( "fslepton_phi", fslepton_phi );

  debug() << " Filled truth <oaltinok_version>" << endmsg;
      // clasify by ccpi0 ccpi0x other
  bool is_ccpi0          = false;
  bool is_cc1pi0         = false;
  bool is_ccpi0secondary = false;
  bool is_by_pim         = false;
  bool is_ccpi0x         = false;
  bool is_other          = false;
  
  int charge = 0;
  if(truth->primaryLepton()==13 ) charge = -1;
  else if(truth->primaryLepton()==-13 ) charge = 1;

  if ( is1pi0(truth)) {
      is_cc1pi0 = true;
      getMCPi0Info(truth);
  }
  
  if ( ispi0( truth, charge ) )	{  is_ccpi0 = true; }
  else if ( ispi0x( truth ) ) is_ccpi0x = true;
  else is_other = true;

  std::pair<bool,bool> secondary_pi0_info = ispi0secondary();
  if (secondary_pi0_info.first)  is_ccpi0secondary = true;
  if (secondary_pi0_info.second) is_by_pim = true;
  
  truth->filtertaglist()->setOrAddFilterTag( "is_ccpi0",          is_ccpi0 );
  truth->filtertaglist()->setOrAddFilterTag( "is_cc1pi0",         is_cc1pi0 );
  truth->filtertaglist()->setOrAddFilterTag( "is_ccpi0secondary", is_ccpi0secondary );
  truth->filtertaglist()->setOrAddFilterTag( "is_by_pim",         is_by_pim);
  truth->filtertaglist()->setOrAddFilterTag( "is_ccpi0x",         is_ccpi0x );
  truth->filtertaglist()->setOrAddFilterTag( "is_other",          is_other );
  
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
        warning() << "No TG4Trajectories found in this gate. Skipping to next one.<oaltinok_version>" << endmsg;
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
    debug() << " CCDeltaPlusAna::shouldAnalyzeMC <oaltinok_version>" << endmsg;

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
        warning() << "TrajectoryMap empty <oaltinok_version>" << endmsg;
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
    debug() << "CCDeltaPlusAna::finalize()<oaltinok_version>" << endmsg;
    
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
    debug() << " Starting ODActivity::EnergeticTower <oaltinok_version>" << endmsg;
    
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
    
    debug() << " Ending ODActivity::EnergeticTower <oaltinok_version>" << endmsg;
    
    return StatusCode::SUCCESS;
}
