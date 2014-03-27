#ifndef ENERGYCHI2CALCULATOR_H_
#define ENERGYCHI2CALCULATOR_H_

// inheritance
#include "MinervaUtils/MinervaHistoTool.h"
#include "EnergyRecTools/IEnergyChi2Calculator.h"
#include "Profiler.h"

// forwards
class IOpticalModel;


/** @class EnergyChi2Calculator EnergyChi2Calculator.h 
 *
 *  EnergyChi2Calculator is a class implementing methods
 *  from the IEnergyChi2Calculator interface.
 *
 */

class EnergyChi2Calculator : public MinervaHistoTool, virtual public IEnergyChi2Calculator {

 public:
  /// Standard constructor
  EnergyChi2Calculator( const std::string& type, const std::string& name, const IInterface* parent );

  /// Destructor
  ~EnergyChi2Calculator(){};

  StatusCode initialize();
 
  void getTotalChi2(const Minerva::Prong* const prong, CalcELoss2Vect& calc_eloss_vect, bool removeOutliers, 
                    std::vector<dEdXFitPoint>& fitPoints, double& chi2Sum, int& ndf) const;
  
  double getChi2(MeasuredELoss meas_loss, CalculatedELoss calc_loss, double& err) const;
 
  StatusCode fillMeasuredELossContainer( const Minerva::Prong* const prong, const std::vector<dEdXFitPoint>& fitPoints, 
                                         MeasELoss2Vect& meas_eloss_cont, bool removeOutliers ) const;

  StatusCode fillMeasuredELoss( const Minerva::Node* const node, MeasuredELoss& meas_eloss, double node_Energy ) const
  

  
  
 
 protected:
 
 private:
  std::vector<dEdXFitPoint> getFitPoints( const MeasELoss2Vect& meas_eloss_2vect, const CalcELoss2Vect& calc_eloss_2vect, 
                                          double& chi2Sum, int& ndf ) const;
					  
  std::vector<double> removedEdXOutliers( const Minerva::Prong* const prong, std::vector<dEdXFitPoint>& fitpoints, 
                                          double& chi2sum, int& ndf ) const;
                                          
  double calcPathLength(const double dRdZ) const;
  double calcdRdZ(const double v_x, const double v_y) const;
  
  void fillProfiler(double dEdX, double posZ, Profiler& p) const;  
  void setProfilerDirection(bool isForward, Profiler& viewX, Profiler& viewU, Profiler& viewV) const;
  bool needRemoval(double node_Z, Profiler& viewX, Profiler& viewU, Profiler& viewV) const;
 
  IOpticalModel*   m_opticalModel;
  std::string      m_opticalModelName;
  std::string      m_opticalModelAlias;

  int              m_maxNClustersProng;
  bool             m_removedEdXOutliers;
  bool             m_removeVertexOutliers;
  bool             m_removeProfileOutliers;
  bool             m_doAdaptiveOutlierRemoval;
  
  double           refDeltadEdX;
 

};

#endif // ENERGYCHI2CALCULATOR_H
