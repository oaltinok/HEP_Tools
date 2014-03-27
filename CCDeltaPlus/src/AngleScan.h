#ifndef CCDeltaPlus_AngleScan_h
#define CCDeltaPlus_AngleScan_h

#include <functional>

#include <TVector2.h>

#include <Event/MinervaEventFwd.h>


class TH1F;

class AngleScan {
  public:
    
    typedef SmartRefVector<Minerva::IDCluster> ShowerCand;
    
    AngleScan(const SmartRefVector<Minerva::IDCluster>& clusters,
              const Gaudi::XYZPoint& vertex);
    
    void Initialize();
    void BuildThetaHistogram();
    void FindPeaks();
    void DoReco();

    unsigned int GetNxCandidate() const;
    unsigned int GetNCandidate() const;

    const std::vector<SmartRefVector<Minerva::IDCluster> >& GetXShowerCandVector() const;
    const std::vector<SmartRefVector<Minerva::IDCluster> >& GetShowerCandVector() const;
    std::vector<Minerva::IDBlob*>  GetShowers() const;

    const SmartRefVector<Minerva::IDCluster>& GetXClusters() const;
    const SmartRefVector<Minerva::IDCluster>& GetUClusters() const;
    const SmartRefVector<Minerva::IDCluster>& GetVClusters() const;

    const SmartRefVector<Minerva::IDCluster>& GetUnusedClusters() const;
    const SmartRefVector<Minerva::IDCluster>& GetUnusedXClusters() const;
    const SmartRefVector<Minerva::IDCluster>& GetUnusedUClusters() const;
    const SmartRefVector<Minerva::IDCluster>& GetUnusedVClusters() const;        

    void SetUVMatchTolerance(double epsilon);
    void SetUVMatchMoreTolerance(double big_epsilon);
    void AllowUVMatchWithMoreTolerance(bool b);
    
  private:
    void FormXShowerCand();
    void FormXUVShowerCand();
    
    int GetLimitBin(const TH1F *hMax, int n_bin )const;

    void addClustersToBlob(SmartRefVector<Minerva::IDCluster>& xshowerCand,
                           SmartRefVector<Minerva::IDCluster>& uclusters,
                           SmartRefVector<Minerva::IDCluster>& vclusters,
                           SmartRefVector<Minerva::IDCluster>& showerCand,
                           double epsilon);

    void completeView(SmartRefVector<Minerva::IDCluster>& unusedViewClusters,
                      SmartRefVector<Minerva::IDCluster>& showerCand,
                      double vtxT);
    
    void coneView(SmartRefVector<Minerva::IDCluster>& unusedViewClusters,
                  SmartRefVector<Minerva::IDCluster>& showerCand,
                  double vtxT,
                  double min_angle, double max_angle,
                  double zmin, double zmax);   
    
    SmartRefVector<Minerva::IDCluster> fAllClusters;
    SmartRefVector<Minerva::IDCluster> fXClusters;
    SmartRefVector<Minerva::IDCluster> fUClusters;
    SmartRefVector<Minerva::IDCluster> fVClusters;
    
    SmartRefVector<Minerva::IDCluster> fRemainingClusters;
    SmartRefVector<Minerva::IDCluster> fRemainingXClusters;
    SmartRefVector<Minerva::IDCluster> fRemainingUClusters;
    SmartRefVector<Minerva::IDCluster> fRemainingVClusters;
    
    double fX;
    double fY;
    double fZ;
    double fU;
    double fV;
    
    TH1F* fTheta;
    std::vector<TVector2> fPeaks;
    std::vector<ShowerCand> fShowerCandidates;
    std::vector<ShowerCand> fXShowerCandidates;

        // Data member controlling the algorithm behaviors
    double fUVMatchTolerance;
    double fUVMatchMoreTolerance;
    bool   fAllowUVMatchWithMoreTolerance;
};

#endif 

