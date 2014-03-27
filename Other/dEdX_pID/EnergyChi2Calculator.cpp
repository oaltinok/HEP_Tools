#include "EnergyChi2Calculator.h"

#include "ReadoutInterfaces/IOpticalModel.h"

#include <algorithm> 

/** @class EnergyChi2Calculator EnergyChi2Calculator.h
 *
 *  @author Brandon Eberly
 *  @date 2012-05-14
 */
 
DECLARE_TOOL_FACTORY( EnergyChi2Calculator );
using namespace Minerva;

namespace {
  const int NSIGMA = 3;  //-- sigma of outliers in the fit
  const double MAX_dEdX_pion = 0.4;
}

//=============================================================================
// Standard Constructor
//=============================================================================
EnergyChi2Calculator::EnergyChi2Calculator(const std::string &type, const std::string &name, const IInterface *parent) :
  MinervaHistoTool(type,name,parent) 
{
  
  declareProperty("OpticalModelName",   m_opticalModelName  = "OpticalModel");
  declareProperty("OpticalModelAlias",  m_opticalModelAlias = "Chi2OpticalModel");
  declareProperty("MaxNClustersProng",  m_maxNClustersProng = 10);
  
  declareProperty("RemovedEdXOutliers", m_removedEdXOutliers = true);
  declareProperty("RemoveVertexOutliers", m_removeVertexOutliers = false);
  declareProperty("RemoveProfileOutliers", m_removeProfileOutliers = false);
  declareProperty("DoAdaptiveOutlierRemoval", m_doAdaptiveOutlierRemoval = false);
  
  // RemoveProfileOutliers Specific Parameter
  declareProperty("RefDeltadEdX",  refDeltadEdX = 0.1);
  
 
  declareInterface<IEnergyChi2Calculator>(this);
}

//=============================================================================
// Initialization
//=============================================================================
StatusCode EnergyChi2Calculator::initialize() {

  StatusCode sc = this->MinervaHistoTool::initialize();
  
  try {
    m_opticalModel = tool<IOpticalModel>(m_opticalModelName, m_opticalModelAlias);
  } catch(GaudiException& e) {
    error() <<"Could not obtain opticalModel" << endmsg;
    return StatusCode::FAILURE;
  }
  
  return sc;
}

//=============================================================================
// return the total chi2, ndf, fit parameters
//=============================================================================
void EnergyChi2Calculator::getTotalChi2(const Minerva::Prong* const prong, CalcELoss2Vect& calc_eloss_2vect, bool removeOutliers,
                                        std::vector<dEdXFitPoint>& fitPoints, double& chi2Sum, int& ndf) const
{
   debug() << "Enter EnergyChi2Calculator::getTotalChi2" << endmsg;
   debug() << "   the size of the fit parameters container= " << fitPoints.size() << endmsg;
   debug() << "   remove outliers? " << removeOutliers << endmsg;

   //! initialize
   chi2Sum = 0.0; 
   ndf     = -1;  
   std::vector<dEdXFitPoint> previous_fitPoints = fitPoints;
   fitPoints.clear();


   //! get measured energy container from the prong information
   MeasELoss2Vect meas_eloss_2vect;
   meas_eloss_2vect.empty();
   StatusCode sc = fillMeasuredELossContainer(prong, previous_fitPoints, meas_eloss_2vect,removeOutliers);
   debug() << "  meas_eloss_vect size = " << meas_eloss_2vect.size() 
           << ", calc_eloss_vect size = " << calc_eloss_2vect.size() << endmsg;
	   
      
   //! Get a set of fitPoints for the measured energyloss container
   fitPoints = getFitPoints( meas_eloss_2vect, calc_eloss_2vect, chi2Sum, ndf );


   //! If removing dEdX Outliers and using adaptive outlier removal, purge the fit points of outliers
   if ( removeOutliers && m_removedEdXOutliers && m_doAdaptiveOutlierRemoval ) {
     removedEdXOutliers(prong, fitPoints, chi2Sum, ndf);
   }
   
   
   debug() << "Exit EnergyChi2Calculator::getTotalChi2" << endmsg;
   return;
}
   


//=================================================================================================
// fill container with the prong's measured dE/dx 
// if the prong is entering a refitting stage, then exclude the dEdXRemovalOutliers
// or exclude the the four nodes on the track
//=================================================================================================
StatusCode EnergyChi2Calculator::fillMeasuredELossContainer( const Minerva::Prong* const prong, const std::vector<dEdXFitPoint>& fitPoints, MeasELoss2Vect& meas_eloss_cont, bool removeOutliers ) const
{ 
  debug() << "Enter EnergyChi2Calculator::fillMeasuredELossContainer" << endmsg;
  StatusCode sc = StatusCode::SUCCESS;
  
  //! @todo - WARNING.  Currently we do not select the ODClusters from the prong because 
  //! there is no unambiguous prescription to select ODClusters that 
  //! were matched to the track (as opposed to those that are associated for calorimetry or proximity)

    
    //! Declare Node Properties
    double node_dRdZ;
    double node_PathLength;
    double node_VisibleEnergy;
    double node_dEdX;
    double node_Z;
    double node_Energy;

    //! Get the tracks from prong and loop over Tracks
    Minerva::TrackVect tracks = prong->minervaTracks();
    for(Minerva::TrackVect::iterator itTrk = tracks.begin(); itTrk != tracks.end(); ++itTrk) {

    //! container of measured dE/dx
    std::vector<MeasuredELoss> meas_eloss_vect;

    //! get the track's nodes
    Minerva::Track::NodeContainer nodes = (*itTrk)->nodes();
    
    //! Declare node Iterator
    Minerva::Track::NodeContainer::const_iterator node;
    Minerva::Track::NodeContainer::const_iterator node_profiler;

    //! initialize
    int x_node_count = 0, u_node_count = 0, v_node_count = 0;
    std::vector<double> zpositions;
    

    //! Declare Profilers
    Profiler viewX('X');
    Profiler viewU('U');
    Profiler viewV('V');
    bool removeNode;
    bool removeNode_profiler;
    
    //! get track direction
    bool isForward;
    isForward = (*itTrk)->direction() == Minerva::Track::Forward;
    
    //! Store the track direction in Profilers
    setProfilerDirection(isForward, viewX, viewU, viewV);


    
    if (removeOutliers && m_removedEdXOutliers && !m_doAdaptiveOutlierRemoval) {
      std::vector<dEdXFitPoint> tmp_fitPoints = fitPoints;
      double tmp_chi2sum = 0.0;
      int tmp_ndf = 0;
      zpositions = removedEdXOutliers(prong, tmp_fitPoints, tmp_chi2sum, tmp_ndf);
      debug() << "  The number of fitPoints before dEdXOutlierRemoval= " << fitPoints.size() << endmsg;
      debug() << "  The number of fitPoints to be removed= " << zpositions.size() << endmsg;
    
    }
    
 
    
    if( removeOutliers && m_removeProfileOutliers)
    {
        debug()<<debugName<<"><><>< Outlier Removal using Profilers Active ><><><"<<endmsg;
        //! Fill Profilers - loop over nodes
        for(node = nodes.begin(); node != nodes.end(); ++node)
        {
            // Calculate node's dEdX
            node_dRdZ = calcdRdZ( (*node)->state().ax(), (*node)->state().ay() );
            node_PathLength = calcPathLength( node_dRdZ );
            node_VisibleEnergy = (*node)->idcluster()->energy();
            node_dEdX = node_VisibleEnergy / node_PathLength;
            node_Z = (*node)->z();
            
            
            //! get the node's cluster
            Minerva::IDCluster* cluster = (*node)->idcluster();
            
            if( !cluster ) {
                debug()<<debugName<<"The node has a NULL cluster - skipping node!" << endmsg;
                continue;
            }
            
            // Fill each Pofiler
            if( cluster->view() == Minerva::IDCluster::X ) {
                fillProfiler(node_dEdX, node_Z, viewX);
            }else if(( cluster->view() == Minerva::IDCluster::U ) ){
                fillProfiler(node_dEdX, node_Z, viewU);
            }else if(( cluster->view() == Minerva::IDCluster::V ) ){
                fillProfiler(node_dEdX, node_Z, viewV);
            }
        }
        
        // Find Outliers for each View
        viewX.findOutliers();
        viewU.findOutliers();
        viewV.findOutliers();

    }
    
    
    
    
    
    //! loop over nodes - AGAIN
//     Minerva::Track::NodeContainer::const_iterator node = nodes.begin();
    for(node = nodes.begin(); node != nodes.end(); ++node) {

        
      //! struct for the meas dE/dx 
      MeasuredELoss tmp_ELoss;

      //! get the node's cluster
      Minerva::IDCluster* cluster = (*node)->idcluster();
      if( !cluster ) {
        info() << "The node has a NULL cluster - but continuing!" << endmsg;
        continue;
      }
      
      
      
      //! remove the nodes with correspond to the dEdXOutlier removal
      if( removeOutliers && m_removedEdXOutliers && !m_doAdaptiveOutlierRemoval && !zpositions.empty() ) {
        bool remove = false;
        for(unsigned int i = 0; i < zpositions.size(); ++i) {
          if( fabs(zpositions[i]-cluster->z()) < 0.001 ) {
            cluster->setIntData("is_dEdXOutlier",1);
            remove = true;
          }
        }
        if( remove ) continue;
      }
      
      
      
      //! remove the nodes that we marked as tobeRemoved
//       if( removeOutliers && m_removeProfileOutliers) {
//       
//         removeNode = needRemoval(cluster->z(), viewX, viewU, viewV);
//         if(removeNode){
//             continue; // Skips the current Node will not include it in meas_eloss_vect
//         }
//             
//       }

      //! Modify the energy of the nodes that we marked as tobeRemoved
      if( removeOutliers && m_removeProfileOutliers) {
        Minerva::IDCluster* cluster_profiler;
        removeNode = needRemoval(cluster->z(), viewX, viewU, viewV);
        if(removeNode){
            removeNode_proviler = removeNode;
            for( node_profiler = node; removeNode_profiler; node_profiler++){
                cluster_profiler = (*node_profiler)->idcluster();
                removeNode_profiler = needRemoval(cluster_profiler->z(), viewX, viewU, viewV);
            }
            node_Energy = *(node_profiler)->idcluster()->energy();
        }else{
            node_Energy =  *node->idcluster()->energy();
        }
            
      }
      
      //! remove the first two x-nodes and first u and v nodes if removeOutliers (near vertex)
      if( removeOutliers && m_removeVertexOutliers) {
        if( cluster->view() == Minerva::IDCluster::X ) {
          x_node_count++;
          if( x_node_count <= 2 ) continue;
        }
        else if( cluster->view() == Minerva::IDCluster::U ) {
          u_node_count++;
          if( u_node_count <= 1 ) continue;
        }
        else if( cluster->view() == Minerva::IDCluster::V ) {
          v_node_count++;       
          if( v_node_count <= 1 ) continue;
        }
      }
      
      //! get the node measured dE/dx information
      StatusCode tmp_sc = fillMeasuredELoss(*node,tmp_ELoss,node_Energy);
      meas_eloss_vect.push_back( tmp_ELoss );

      //! sanity check
      if( tmp_sc.isFailure() ) sc = tmp_sc;
    }
    
    //! sanity check
    if( meas_eloss_vect.size() > 1 ) {
      debug() << "  1st sanity check on meas_eloss container with size: " << meas_eloss_vect.size() << endmsg;

      std::vector<MeasuredELoss>::iterator itELoss = meas_eloss_vect.begin();
      for( ; itELoss != meas_eloss_vect.end()-1; ++itELoss) {
        debug() << "    check" << endmsg;
        std::vector<MeasuredELoss>::iterator itNext = itELoss+1;
        if( (*itTrk)->direction() == Minerva::Track::Forward && itELoss->z_position > itNext->z_position ) {
            warning() << "Forward track has meas_eloss points in wrong direction!" << endmsg;
        }
        else if( (*itTrk)->direction() == Minerva::Track::Backward && itELoss->z_position < itNext->z_position ) {
            warning() << "Backward track has meas_eloss points in wrong direction!" << endmsg;
        }
        else if( itELoss->z_position == itNext->z_position ) {
            warning() << "Track has equal meas_eloss points!"<<endmsg;
        }
      }
    }
    
    //! store the measured energy loss per track
    meas_eloss_cont.push_back( meas_eloss_vect );
  }//! end loop over tracks
  
  //! @todo - Get the prong's OD match clusters and fill an OD MeasuredELoss vector
  
  debug() << "Exit EnergyChi2Calculator::fillMeasuredELossContainer" << endmsg;
  return sc;
}



//================================================================================================
// store the measure dE/dx per node into a container
//================================================================================================
StatusCode EnergyChi2Calculator::fillMeasuredELoss( const Minerva::Node* const node, MeasuredELoss& meas_eloss, double node_Energy ) const
{
  verbose() << "Enter EnergyChi2Calculator::fillMeasuredELoss" << endmsg;

  double v_x  = node->state().ax();
  double v_y  = node->state().ay();

  meas_eloss.z_position      = node->z();
  meas_eloss.x_position      = node->position().x();
  meas_eloss.y_position      = node->position().y();
//   meas_eloss.visible_energy  = node->idcluster()->energy();
  meas_eloss.visible_energy  = node_Energy;
  meas_eloss.PEs             = node->idcluster()->pe();
  meas_eloss.dRdZ            = sqrt(v_x*v_x + v_y*v_y);

  verbose() << "Exit EnergyChi2Calculator::fillMeasuredELoss" << endmsg;
  return StatusCode::SUCCESS;
}

bool EnergyChi2Calculator::needRemoval(double node_Z, Profiler& viewX, Profiler& viewU, Profiler& viewV) const
{
    bool result;
    
    // result is True if any of the views needs Removal
    result = viewX.needRemoval(node_Z) || viewU.needRemoval(node_Z) || viewV.needRemoval(node_Z);
    
    return result;

}

void EnergyChi2Calculator::setProfilerDirection(bool isForward, Profiler& viewX, Profiler& viewU, Profiler& viewV) const
{
    viewX.set_isForward(isForward);
    viewU.set_isForward(isForward);
    viewV.set_isForward(isForward);
    debug()<<debugName<<"Track Direction saved in Profilers"<<endmsg;
}

void EnergyChi2Calculator::fillProfiler(double dEdX, double posZ, Profiler& p) const
{
    debug()<<debugName<<"Filling Profiler"<<p.get_view()<<"-view"<<endmsg;
    p.add_dEdX(dEdX);
    p.add_posZ(posZ);
}

double EnergyChi2Calculator::calcdRdZ(const double v_x, const double v_y) const
{
    double dRdZ;
    
    dRdZ = sqrt( (v_x * v_x) + (v_y * v_y) );
    
    return dRdZ;
}

double EnergyChi2Calculator::calcPathLength(const double dRdZ) const
{
    double pathLength;
    double R;
    const double planeLength = 17.0;  // [mm]
    const double planeLengthSq = planeLength * planeLength;
    
    R = dRdZ * planeLength;
        
    pathLength = sqrt( planeLengthSq + (R * R) );
    
    return pathLength;
}



//=============================================================================
// Get a set of fit points for the measured and calculated energy containers
// Also store the chi2 sum and ndf
//=============================================================================
std::vector<dEdXFitPoint> EnergyChi2Calculator::getFitPoints( const MeasELoss2Vect& meas_eloss_2vect, const CalcELoss2Vect& calc_eloss_2vect, double& chi2Sum, int& ndf ) const
{
  
  chi2Sum = 0.0;
  ndf = -1;
  std::vector<dEdXFitPoint> returnFitPoints;
  
  //! loop over the measured dE/dx observables and calculate chi2
  for(unsigned int meas_v = 0; meas_v < meas_eloss_2vect.size(); ++meas_v) {
    if( meas_v >= calc_eloss_2vect.size() ) break;

    std::vector<MeasuredELoss>   meas_eloss = meas_eloss_2vect[meas_v];
    std::vector<CalculatedELoss> calc_eloss = calc_eloss_2vect[meas_v];

    for(unsigned int meas = 0; meas < meas_eloss.size(); ++meas) {
      double measZpos = meas_eloss[meas].z_position;

      verbose() << "   measured z position = " << measZpos 
        	<< ", measVisibleEnergy = " << meas_eloss[meas].visible_energy
        	<< ", measPE = " << meas_eloss[meas].PEs << endmsg;

      //! loop over the calculated dE/dx observables 
      for(unsigned int calc = 0; calc < calc_eloss.size(); ++calc) {
	double calZpos = calc_eloss[calc].z_position;

	//! If z positions match, calculate a chi2 
	if( fabs(calZpos - measZpos) < 1.0e-3*CLHEP::mm && calc_eloss[calc].path_length > 0.0 ) {

          double total_error = 0.0;
          double chi2 = getChi2(meas_eloss[meas],calc_eloss[calc],total_error);
          if(chi2 < 0.0) {
            counter("NEGATIVE_NODE_CHI2")++;
            continue;
          }

          chi2Sum += chi2;
          ++ndf;

          dEdXFitPoint tmp_fitpoint;
	  tmp_fitpoint.Measured = (meas_eloss[meas].visible_energy)/(calc_eloss[calc].path_length); //calc_eloss pathlength is better
          tmp_fitpoint.Calc     = (calc_eloss[calc].visible_energy)/(calc_eloss[calc].path_length);
          tmp_fitpoint.Sigma    = total_error;
	  tmp_fitpoint.Chi2     = chi2;
	  tmp_fitpoint.ZPos     = measZpos;
          returnFitPoints.push_back(tmp_fitpoint);

	  break;
	}
      }
    }
  } //loop over meas double vector vector
  
  return returnFitPoints;
}


//=============================================================================
// return chi2 given measured and predicted dE, straggling width, and pathlength
//=============================================================================
double EnergyChi2Calculator::getChi2(MeasuredELoss meas_eloss, CalculatedELoss calc_eloss, double& err) const {

  err = 0.0;
  if (calc_eloss.path_length <= 0.0) return -9.9;
  if ( fabs(meas_eloss.z_position - calc_eloss.z_position) > 1.0e-3*CLHEP::mm ) {
    info()<<"getChi2 given measured and calculated energy losses at different z positions.  Returning"<<endmsg;
    return -9.9;
  }
  
  double measured_dE = meas_eloss.visible_energy;
  double predicted_dE = calc_eloss.visible_energy;
  double dE_residual = predicted_dE - measured_dE;
   
  //! estimate photostatistical error on measured_dE
  double measured_dE_error = 0.0;
  if (meas_eloss.PEs > 0.0) {
    measured_dE_error = measured_dE/sqrt( meas_eloss.PEs );
  }
  
  //! estimate photostatistical error on predicted_dE
  //! 1720.0 is an effective average attenuation*S2S constant derived from v10r2 ReadOut Simulation
  double predicted_dE_error = 0.0;
  if (predicted_dE > 0.0) {
    predicted_dE_error = predicted_dE/sqrt( m_opticalModel->birksPe(predicted_dE,calc_eloss.path_length)/1720.0 );
  }
  
  //! estimate the pathlength error as the difference in pathlength due to the tracking theta residual.  The phi residual does not
  //! change the total pathlength through a plane, but can have a small effect on the total pathlength through the active material 
  //! in a plane. Including a second order term for when theta approaches zero may be necessary if the muon theta uncertainty is greater than 0.1 radians
  double pathLength_error = dE_residual*( (meas_eloss.dRdZ)*(0.04*CLHEP::radian) );  //0.04 ~ 2.5 degrees -> typical theta residual
  
  double total_error2 = pow(calc_eloss.straggling_width,2.0) + pow(measured_dE_error,2.0) + pow(predicted_dE_error,2.0) + pow(pathLength_error,2.0);
  if (total_error2 > 0.0) {
    double chi2 = dE_residual*dE_residual/total_error2;

    if (measured_dE > 0.0) {
      plot1D( measured_dE_error/measured_dE, "mdE_err_over_mdE","Photostatistical error divided by Measured Energy Loss",0.0,1.0,500);
    }
    plot1D(predicted_dE_error/predicted_dE, "pdE_err_over_pdE","Photostatistical error divided by Predicted Energy Loss",0.0,1.0,500);
    plot1D(chi2, "chi2","Chi-Square",0.0,100.0,500);
    plot2D(measured_dE/(calc_eloss.path_length), chi2, "chi2_v_dEdX", "Chi-Square vs. Measured dEdX", 0.0,10.0, 0.0,100.0, 1000, 50);

    err = sqrt(total_error2)/(calc_eloss.path_length);  //convert from error on dE to error on dEdX

    return chi2;
  }
  else {
    return -9.9;
  }
}


//========================================================================================
// remove dEdx outliers from the fitPoints using chi2 and update chi2 and ndf totals
//========================================================================================
std::vector<double> EnergyChi2Calculator::removedEdXOutliers(const Minerva::Prong* const prong, std::vector<dEdXFitPoint>& fitpoints, double& chi2Sum, int& ndf) const
{

  verbose() << "Enter EnergyChi2Calculator::removedEdXOutliers" << endmsg;

  std::vector<double> zpos_vec;

  if (fitpoints.empty()) {
    warning() << "Do not have fit information!" << endmsg;
    return zpos_vec;
  }
  

  //! Check number of clusters
  Minerva::TrackVect tracks = prong->minervaTracks();
  int nclusters = 0;
  for(Minerva::TrackVect::iterator trk = tracks.begin(); trk != tracks.end(); ++trk) {
    nclusters += (*trk)->idclusters().size(); 
  }
  if (nclusters <= m_maxNClustersProng) {
    return zpos_vec;
  }
  
  double x = 0.0, e = 0.0, n = double(fitpoints.size());

  for (std::vector<dEdXFitPoint>::iterator itFit = fitpoints.begin(); itFit != fitpoints.end(); ++itFit) {
    x += itFit->Chi2;
  }
  x = x/n;
  
  for (std::vector<dEdXFitPoint>::iterator itFit = fitpoints.begin(); itFit != fitpoints.end(); ++itFit) {
    e += pow(itFit->Chi2 - x,2.0);
  }  
  double sigma = sqrt(e/n);


  for (std::vector<dEdXFitPoint>::iterator itFit = fitpoints.begin(); itFit != fitpoints.end(); ++itFit) {
     
    //! Check whether chi2 is too large
    if (itFit->Chi2 > x + NSIGMA*sigma) {
      zpos_vec.push_back( itFit->ZPos );      
      chi2Sum -= itFit->Chi2;
      ndf--;
      itFit = fitpoints.erase(itFit);
      itFit--;
    }//end if chi2 is large
    
  }//end loop over fitPoints

  return zpos_vec;
}
