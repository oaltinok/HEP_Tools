/*
    See AngleScan.h for Class Information
*/
#undef NDEBUG
#include <cassert>
#include <cmath>

#include <TMath.h>
#include <TH1F.h>

#include <Event/IDCluster.h>
#include <Event/IDBlob.h>
#include <Event/VectorTypeDefs.h>
#include <Event/Vertex.h>

#include "AngleScan.h"
#include "ClusterVectorInfo.h"

namespace {
    struct greaterShower : public std::binary_function <
        SmartRefVector<Minerva::IDCluster>,
        SmartRefVector<Minerva::IDCluster>,
        bool > {

        bool operator()(const SmartRefVector<Minerva::IDCluster>& lhs,
                        const SmartRefVector<Minerva::IDCluster>& rhs) const {
            
            return lhs.size() > rhs.size();
        }

    };

    class zDistanceFromLessThan : public std::binary_function <
        SmartRef<Minerva::IDCluster>,
        SmartRef<Minerva::IDCluster>,
        bool > {
    public:
        explicit zDistanceFromLessThan(double z0) : fZ0(z0) {}
        bool operator() (const SmartRef<Minerva::IDCluster>& lhs) const {
            
            return std::abs(lhs->z() - fZ0) < 50.0;
        }

    private:
        double fZ0;
    };
}


AngleScan::AngleScan(const SmartRefVector<Minerva::IDCluster>& clusters,
                     const Gaudi::XYZPoint& vertex)
    : fUVMatchTolerance(10.0),
      fUVMatchMoreTolerance(100.0),
      fAllowUVMatchWithMoreTolerance(true)
{
    // Save All Clusters to Data Variable fAllClusters
    std::copy(clusters.begin(), clusters.end(), std::back_inserter(fAllClusters));
    
    // Save Vertex Position Information
    vtx_pos_X = vertex.X();
    vtx_pos_Y = vertex.Y();
    vtx_pos_Z = vertex.Z();
    vtx_pos_U = -vtx_pos_Y*sqrt(3.)/2 + vtx_pos_X/2;
    vtx_pos_V =  vtx_pos_Y*sqrt(3.)/2 + vtx_pos_X/2;
    
    Initialize();
}

void AngleScan::Initialize()
{
    std::cout << "    AngleScan::Initialize " << std::endl;

    // Get Clusters in Each View(X,U,V) using ClusterVectorInfo Class
    ClusterVectorInfo clusterVectorInfo(fAllClusters,true,true,true);
    fXClusters = clusterVectorInfo.GetXClusters();
    fUClusters = clusterVectorInfo.GetUClusters();
    fVClusters = clusterVectorInfo.GetVClusters();

    std::cout << "      AngleScan:: x size = " << fXClusters.size() << std::endl;
    std::cout << "      AngleScan:: u size = " << fUClusters.size() << std::endl;
    std::cout << "      AngleScan:: v size = " << fVClusters.size() << std::endl;

    // Copy to the working containers
    fRemainingXClusters = fXClusters;
    fRemainingUClusters = fUClusters;
    fRemainingVClusters = fVClusters;

    // Sanity Check!
    assert(fRemainingXClusters.size() == fXClusters.size());
    assert(fRemainingUClusters.size() == fUClusters.size());
    assert(fRemainingVClusters.size() == fVClusters.size());
}

//==============================================================================
// Builds Histogram for Angle Theta
//==============================================================================
void AngleScan::BuildThetaHistogram()
{
    const unsigned int N = 90;
    const double lower = -180.0;
    const double upper = +180.0;
    fTheta = new TH1F("theta","Theta dist", N, lower, upper);
    
    for (   SmartRefVector<Minerva::IDCluster>::const_iterator c = fRemainingXClusters.begin();
            c != fRemainingXClusters.end();
            ++c) 
    {
        const double dZ = (*c)->z() - vtx_pos_Z;
        const double dX = (*c)->position() - vtx_pos_X;
        const double theta = std::atan2(dX,dZ)*TMath::RadToDeg();
        const double w = (*c)->pe();
        
        fTheta->Fill(theta,w);
    }
}

//==============================================================================
// Detect and save lower and upper edges around peaks in the  Theta Histogram 
//==============================================================================
void AngleScan::FindPeaks()
{
    const unsigned int N = 90;
    
    const double width = fTheta->GetBinWidth(1);      // Use first bin since fixed bin size
     
    // Loop over Bins
    for (unsigned int bin = 1; bin <= N; bin++){
        if (fTheta->GetBinContent(bin) > 15 ){        // Peak is detected 
            int Limitbin = GetLimitBin(fTheta, bin ); // Find the upper edge of the peak 
            
            const double lower_edge = fTheta->GetBinCenter(bin) - 0.5*width;
            const double upper_edge = fTheta->GetBinCenter(Limitbin) + 0.5*width;
            
            TVector2 peak(lower_edge,upper_edge);

            fPeaks.push_back(peak);
                        
            bin = Limitbin;  /* Finding the next peak starting from this peak's upper edge */
            
        }
    }
    
    delete fTheta;
}


//==============================================================================
// Form Shower Candidate Vector in X View
//==============================================================================
void AngleScan::FormXShowerCand() 
{
    // Loop over Peaks
    for (std::vector<TVector2>::const_iterator peak = fPeaks.begin();
        peak != fPeaks.end(); 
        ++peak) 
    {
        const double lower_edge = peak->X();
        const double upper_edge = peak->Y();
        const double zmin       = 4500.0;
        const double zmax       = 10000.0;
        SmartRefVector<Minerva::IDCluster> showerCand;
        
        // Send peak information to coneView
        coneView(fRemainingXClusters, showerCand, vtx_pos_X, lower_edge, upper_edge, zmin, zmax);
        
        if (showerCand.empty()) continue;

        if (showerCand.size() == 1 && showerCand.front()->pe() > 30) {
            fXShowerCandidates.push_back(showerCand);
        } else if (showerCand.size() > 1) {
            fXShowerCandidates.push_back(showerCand);
        }
    }

    std::sort(fXShowerCandidates.begin(),fXShowerCandidates.end(), greaterShower());

}

//==============================================================================
// Form Shower Candidate Vector in All Views
//==============================================================================
void AngleScan::FormXUVShowerCand() 
{
    std::vector<SmartRefVector<Minerva::IDCluster> > nogrowShowerCandidates;
    std::cout << "AngleScan::total X candidates: " << fXShowerCandidates.size() << std::endl;
    
    // Loop over Shower Candidates
    for (   std::vector<SmartRefVector<Minerva::IDCluster> >::iterator s = fXShowerCandidates.begin();
            s != fXShowerCandidates.end(); 
            ++s) 
    {
        std::cout << "\tAngleScan::x cand: " << s->size() << std::endl;
        
        SmartRefVector<Minerva::IDCluster>& xshowerCand = *s;
        SmartRefVector<Minerva::IDCluster> showerCand = xshowerCand;
        
        double zmin = +1e6;
        double zmax = -1e6;
        double ztot = 0.0;
        
        // Loop over Clusters to Get Total Z of Shower
        for (   SmartRefVector<Minerva::IDCluster>::iterator c = xshowerCand.begin();
                c != xshowerCand.end(); 
                ++c) 
        {
            zmin  = std::min(zmin,(*c)->z());
            zmax  = std::max(zmax,(*c)->z());
            ztot += (*c)->z();
        }

        if (zmax-zmin < 50.0) {
            const double zcenter = ztot/xshowerCand.size();
            std::cout << "\t AngleScan::Candidate with small z extend at: " << zcenter << std::endl;
            SmartRefVector<Minerva::IDCluster> uclusters_tmp;
            SmartRefVector<Minerva::IDCluster> vclusters_tmp;
            fRemainingUClusters.swap(uclusters_tmp);
            fRemainingVClusters.swap(vclusters_tmp);

            // Loop over Clusters on U View
            for (   SmartRefVector<Minerva::IDCluster>::iterator c = uclusters_tmp.begin();
                    c != uclusters_tmp.end(); 
                    ++c) 
            {
                if (std::abs((*c)->z()-zcenter) < 50.0) {
                    std::cout << "\t\t AngleScan::adding a U cluster " << (*c)->z() << " " << (*c)->pe()
                              << std::endl;
                    showerCand.push_back(*c);
                }
                else fRemainingUClusters.push_back(*c);
            }
            // Loop over Clusters on V View
            for (   SmartRefVector<Minerva::IDCluster>::iterator c = vclusters_tmp.begin();
                    c != vclusters_tmp.end(); 
                    ++c) 
            {
                if (std::abs((*c)->z()-zcenter) < 50.0) {
                    std::cout << "\t\t AngleScan::adding a V cluster " << (*c)->z() << " " << (*c)->pe()
                              << std::endl;
                    showerCand.push_back(*c);
                }
                else fRemainingVClusters.push_back(*c);
            }
            
            fShowerCandidates.push_back(showerCand);

        } else {
            std::cout << "\t AngleScan::Candidate with large z extent " << std::endl;
            addClustersToBlob(xshowerCand,fRemainingUClusters,fRemainingVClusters,showerCand,
                              fUVMatchTolerance);
            
            if (showerCand.size() > xshowerCand.size()) fShowerCandidates.push_back(showerCand);
            else if (xshowerCand.size() >= 3) nogrowShowerCandidates.push_back(xshowerCand);
            else {
                std::cout << "\t AngleScan::Throw away x shower candidate" << std::endl;
            }
        }
    }

    std::cout << "AngleScan::Shower candidates: " << fShowerCandidates.size() << std::endl;
    std::cout << "AngleScan::Nogrow candidates: " << nogrowShowerCandidates.size() << std::endl;

    if (fAllowUVMatchWithMoreTolerance) {
        for (   std::vector<SmartRefVector<Minerva::IDCluster> >::iterator s = nogrowShowerCandidates.begin();
                s != nogrowShowerCandidates.end(); 
                ++s) 
        {
            SmartRefVector<Minerva::IDCluster>& xshowerCand = *s; 
            SmartRefVector<Minerva::IDCluster>  showerCand = xshowerCand;
            
            addClustersToBlob(xshowerCand,fRemainingUClusters,fRemainingVClusters,showerCand,fUVMatchMoreTolerance);
            
            fShowerCandidates.push_back(showerCand);
        }
    }
    
    std::cout << "AngleScan::Shower candidates(10x): " << fShowerCandidates.size() << std::endl;
    for (   std::vector<SmartRefVector<Minerva::IDCluster> >::iterator s = fShowerCandidates.begin();
            s != fShowerCandidates.end(); 
            ++s) 
    {
        std::cout << "\tAngleScan:: cand: " << std::distance(fShowerCandidates.begin(),s) << " "
                  << fRemainingUClusters.size() << " " << fRemainingVClusters.size()
                  << std::endl;
        completeView(fRemainingUClusters, *s, vtx_pos_U);
        completeView(fRemainingVClusters, *s, vtx_pos_V);
    }

    std::cout << "AngleScan::Final showers: " << fShowerCandidates.size() << std::endl;
    std::sort(fShowerCandidates.begin(), fShowerCandidates.end(), greaterShower());
    
    std::copy(fRemainingXClusters.begin(),fRemainingXClusters.end(), std::back_inserter(fRemainingClusters));
    std::copy(fRemainingUClusters.begin(),fRemainingUClusters.end(), std::back_inserter(fRemainingClusters));
    std::copy(fRemainingVClusters.begin(),fRemainingVClusters.end(), std::back_inserter(fRemainingClusters));

}

void AngleScan::DoReco()
{
    std::cout << "AngleScan::DoReco() " << std::endl;
    BuildThetaHistogram();
    FindPeaks();
    FormXShowerCand();
    FormXUVShowerCand();
}

int AngleScan::GetLimitBin(const TH1F *hMax, int n_bin )const
{

    int max=n_bin, x = 1, count = 0;
    while ( count < 1 ){
        if ( hMax->GetBinContent(n_bin + x) > 0 ) max = n_bin + x;
        else count++;
        x++;
    }
    
    return max;

}
void AngleScan::addClustersToBlob(SmartRefVector<Minerva::IDCluster>& xshowerCand,
                                  SmartRefVector<Minerva::IDCluster>& uclusters,
                                  SmartRefVector<Minerva::IDCluster>& vclusters,
                                  SmartRefVector<Minerva::IDCluster>& showerCand,
                                  double epsilon)
{

    for (   SmartRefVector<Minerva::IDCluster>::iterator c = xshowerCand.begin();
            c != xshowerCand.end(); 
            ++c) 
    {
        Minerva::IDCluster* cluster_x = *c;
        double min = 1e3;
        SmartRefVector<Minerva::IDCluster>::iterator ucluster = uclusters.end();
        SmartRefVector<Minerva::IDCluster>::iterator vcluster = vclusters.end();
        
        // Loop over Clusters in U View
        for (   SmartRefVector<Minerva::IDCluster>::iterator itU = uclusters.begin();
                itU != uclusters.end(); 
                ++itU) 
        {
            if (std::abs( cluster_x->z() - (*itU)->z() ) > 50.0 ) continue;
            
            // Loop over Clusters in V View
            for (   SmartRefVector<Minerva::IDCluster>::iterator itV = vclusters.begin();
                    itV != vclusters.end(); 
                    ++itV) 
            {
                if ( std::abs( cluster_x->z() - (*itV)->z() ) > 50.0 ) continue;
                
                // |u+v-x|
                double delta = std::abs((*itU)->tpos1()+(*itU)->tpos2()+       
                                        (*itV)->tpos1()+(*itV)->tpos2()-
                                        cluster_x->tpos1()-cluster_x->tpos2());
                                        
                if ( delta < min ) { 
                    min = delta;
                    ucluster = itU;
                    vcluster = itV;
                }
            }
        }

        if (min <= epsilon && (ucluster != uclusters.end() && vcluster != vclusters.end())) {
            showerCand.push_back(*ucluster);
            showerCand.push_back(*vcluster);

            uclusters.erase(ucluster);
            vclusters.erase(vcluster);
        }
    }
}

void AngleScan::completeView(SmartRefVector<Minerva::IDCluster>& unusedViewClusters,
                             SmartRefVector<Minerva::IDCluster>& showerCand,
                             double vtxT)
{

    if (unusedViewClusters.empty()) return;
    
    double z_min = 10000;
    double z_max = -10000;
    double angle_min = 180;
    double angle_max = -180;

    SmartRefVector<Minerva::IDCluster> viewShowerCand;
    
    SmartRefVector<Minerva::IDCluster>::iterator it_view = unusedViewClusters.begin();
    for (   SmartRefVector<Minerva::IDCluster>::iterator c = showerCand.begin();
            c != showerCand.end(); 
            ++c)
    {
        if ( (*c)->view() == (*it_view)->view() ) viewShowerCand.push_back(*c);
        if ( (*c)->view() == Minerva::IDCluster::X && (*c)->z() < z_min ) z_min = (*c)->z();
        if ( (*c)->view() == Minerva::IDCluster::X && (*c)->z() > z_max ) z_max = (*c)->z();
    }

    if (viewShowerCand.empty()) {
        std::cout << "\tAngleScan::completeView: no view cluster" << std::endl;
        return;
    }

    for ( it_view = viewShowerCand.begin(); it_view != viewShowerCand.end(); ++it_view){
        const double z   = (*it_view)->z() - vtx_pos_Z;
        const double x   = (*it_view)->position() - vtxT;
        const double ang = std::atan2(x,z)*TMath::RadToDeg();
        
        if ( ang >= angle_max ) angle_max = ang;
        if ( ang <= angle_min ) angle_min = ang;
        
    }

    z_min = z_min - 100;
    z_max = z_max + 100;
    angle_max = angle_max + 10.0;
    angle_min = angle_min - 10.0;
    
    /* Move clusters between (angle_min,angle_max) and (z_min,z_max) from
           'unusedClusters' to blob */
    std::cout << "\tAngleScan::completeView:oldsize: " << showerCand.size() << std::endl;
    coneView(unusedViewClusters, showerCand, vtxT, angle_min, angle_max, z_min, z_max );
    std::cout << "\tAngleScan::completeView:newsize: " << showerCand.size() << std::endl;
 
}

//==============================================================================
// Returns Clusters as showerCand if they are in Cone defined by angle and z
//==============================================================================
void AngleScan::coneView(SmartRefVector<Minerva::IDCluster>& unusedViewClusters,
                         SmartRefVector<Minerva::IDCluster>& showerCand,
                         double vtx_x,
                         double min_angle, double max_angle,
                         double zmin, double zmax)
{
    SmartRefVector<Minerva::IDCluster> tmpClusters;
    tmpClusters.swap(unusedViewClusters);
    for (SmartRefVector<Minerva::IDCluster>::iterator c = tmpClusters.begin();
         c != tmpClusters.end(); 
         ++c) 
    {
        // Calculate Theta
        const double dZ = (*c)->z() - vtx_pos_Z;
        const double dX = (*c)->position() - vtx_x;
        const double theta = std::atan2(dX,dZ)*TMath::RadToDeg();
        const double z_c = (*c)->z();
        
        // If theta is between peak range and z in range
        // save cluster in showerCand
        if ((min_angle <= theta && theta <= max_angle) &&
            (zmin < z_c && z_c < zmax)) {
            showerCand.push_back(*c);
        } else {
            unusedViewClusters.push_back(*c);
        }
    }
}

//==============================================================================
// Get Functions
//==============================================================================
std::vector<Minerva::IDBlob*> AngleScan::GetShowers() const 
{
    std::vector<Minerva::IDBlob*> finalBlobs;
    for (std::vector<SmartRefVector<Minerva::IDCluster> >::const_iterator
             s = fShowerCandidates.begin();
         s != fShowerCandidates.end(); ++s) {
        Minerva::IDBlob* newBlob = new Minerva::IDBlob;
        newBlob->add(*s);
        finalBlobs.push_back(newBlob);
    }

    return finalBlobs;
}

unsigned int AngleScan::GetNxCandidate() const 
{
    return fXShowerCandidates.size();
}

unsigned int AngleScan::GetNCandidate() const 
{
    return fShowerCandidates.size();
}

const std::vector<SmartRefVector<Minerva::IDCluster> >& AngleScan::GetXShowerCandVector() const 
{
    return fXShowerCandidates;
}

const std::vector<SmartRefVector<Minerva::IDCluster> >& AngleScan::GetShowerCandVector() const 
{
    return fShowerCandidates;
}

const SmartRefVector<Minerva::IDCluster>& AngleScan::GetXClusters() const 
{
    return fXClusters;
}

const SmartRefVector<Minerva::IDCluster>& AngleScan::GetUClusters() const 
{
    return fUClusters;
}

const SmartRefVector<Minerva::IDCluster>& AngleScan::GetVClusters() const 
{
    return fVClusters;
}

const SmartRefVector<Minerva::IDCluster>& AngleScan::GetUnusedClusters() const 
{
    return fRemainingClusters;
}

const SmartRefVector<Minerva::IDCluster>& AngleScan::GetUnusedXClusters() const 
{
    return fRemainingXClusters;
}

const SmartRefVector<Minerva::IDCluster>& AngleScan::GetUnusedUClusters() const 
{
    return fRemainingUClusters;
}

const SmartRefVector<Minerva::IDCluster>& AngleScan::GetUnusedVClusters() const 
{
    return fRemainingVClusters;
}

//==============================================================================
// Set Functions for Data Members Controlling Algorithm Flow
//==============================================================================
void AngleScan::SetUVMatchTolerance(double epsilon) 
{
    fUVMatchTolerance = epsilon;
}

void AngleScan::SetUVMatchMoreTolerance(double big_epsilon) 
{
    fUVMatchMoreTolerance = big_epsilon;
}

void AngleScan::AllowUVMatchWithMoreTolerance(bool b) 
{
    fAllowUVMatchWithMoreTolerance = b;
}
