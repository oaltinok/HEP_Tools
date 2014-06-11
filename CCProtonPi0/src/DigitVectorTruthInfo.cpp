/*
    DigitVectorTruthInfo.cpp Duplicated from CCPi0 Package on 2014-06-10
        Purpose: Make CCProtonPi0 Package independent of CCPi0 Package
        Future: Common Tools and Functions will be combined under AnaUtils or
                PionUtils
                
    See DigitVectorTruthInfo.h for Class Information
*/
#undef NDEBUG
#include <cassert>
#include <set>
#include <iterator>
#include <iostream>


#include <Event/IDDigit.h>
#include <Event/MCIDDigit.h>
#include <Event/MCHit.h>
#include <Event/TG4Trajectory.h>


#include "DigitVectorTruthInfo.h"
#include "TraverseHistory.h"

DigitVectorTruthInfo::DigitVectorTruthInfo()
    : fTotalNormEnergy(0.0),
      fDataNormEnergy(0.0),
      fMCNormEnergy(0.0),
      fMCXtalkNormEnergy(0.0),
      fTotalTruthEnergy(0.0),
      fSharedTruthEnergy(0.0)
{}

double DigitVectorTruthInfo::GetTotalNormEnergy() const {
    return fTotalNormEnergy;
}
double DigitVectorTruthInfo::GetDataNormEnergy() const {
    return fDataNormEnergy;
}
double DigitVectorTruthInfo::GetMCNormEnergy() const {
    return fMCNormEnergy;
}
double DigitVectorTruthInfo::GetMCXtalkNormEnergy() const {
    return fMCXtalkNormEnergy;
}
double DigitVectorTruthInfo::GetTotalTruthEnergy() const {
    return fTotalTruthEnergy;
}
double DigitVectorTruthInfo::GetSharedTruthEnergy() const {
    return fSharedTruthEnergy;
}

double DigitVectorTruthInfo::GetDataFraction() const {
    return fDataNormEnergy/fTotalNormEnergy;
}
double DigitVectorTruthInfo::GetXtalkFraction() const {
    return fMCXtalkNormEnergy/fMCNormEnergy;
}
double DigitVectorTruthInfo::GetSharedFraction() const {
    return fSharedTruthEnergy/fTotalTruthEnergy;
}

double DigitVectorTruthInfo::GetEdepByPdg(int pdg) const {
    std::map<int,double>::const_iterator target = fPdgEdepMap.find(pdg);
    if (target != fPdgEdepMap.end()) return target->second;

    return 0.0; /* If the particle does not exist, return 0.0 */
}

double DigitVectorTruthInfo::GetEdepByTrackId(int trackId) const {
    std::map<int,double>::const_iterator target = fNextToPrimaryIdEdepMap.find(trackId);
    if (target != fNextToPrimaryIdEdepMap.end()) return target->second;

    return 0.0; /* If trackId does not exist */
}

/* Return the pdg of the primary particle contributing most to this digit vector */
int DigitVectorTruthInfo::GetMostEvisPdg() const {

    if (fPdgEdepMap.empty()) return -1;
    
    double max = -1.e6;
    int pdg = -1;
    for (std::map<int,double>::const_iterator p = fPdgEdepMap.begin();
         p != fPdgEdepMap.end(); ++p) {
        if (p->second > max) {
            max = p->second;
            pdg = p->first;
        }
    }

    assert(pdg != -1);

    return pdg;
}

const std::map<int,double>& DigitVectorTruthInfo::GetPdgEdepMap() const {
    return fPdgEdepMap;
}

const std::map<int,double>& DigitVectorTruthInfo::GetNextToPrimaryIdEdepMap() const {
    return fNextToPrimaryIdEdepMap;
}

void DigitVectorTruthInfo::ParseTruth(const SmartRefVector<Minerva::IDDigit>& digits,
                                      const std::map<int, TG4Trajectory*>& trajectories)
{

    for (SmartRefVector<Minerva::IDDigit>::const_iterator d = digits.begin();
         d != digits.end(); ++d) {
        const double energy = (*d)->normEnergy();
        fTotalNormEnergy += energy;
        
        const Minerva::IDDigit* digit = *d;
        const Minerva::MCIDDigit* mcdigit
            = dynamic_cast<const Minerva::MCIDDigit*>(digit);
        if (!mcdigit) {
            fDataNormEnergy += energy;
            continue;
        }

            /* This is an MC digit, including MC x-talk digit */
        fMCNormEnergy += energy;
        
        const SmartRefVector<Minerva::MCHit>& hits = mcdigit->hits();
        
        if (hits.empty()) {
            fMCXtalkNormEnergy += energy;
            continue;
        }
        
            /* Handle non x-talk, MC digits */
            /* Implementation note: This part could be replaced by HitVectorTruthInfo */
        double total = 0.0;
        std::set<int> primaryContributors;
        for (SmartRefVector<Minerva::MCHit>::const_iterator h = hits.begin();
             h != hits.end(); ++h) {
            total += (*h)->energy();
            
            TraverseHistory history((*h)->GetTrackId());
            history.DoTraverse(trajectories);
            int primaryPdg = history.GetPrimaryPdg();
            int nextToPrimaryId = history.GetNextToPrimaryId();
            
            fPdgEdepMap[primaryPdg] += (*h)->energy();
            fNextToPrimaryIdEdepMap[nextToPrimaryId] += (*h)->energy();
            
            primaryContributors.insert(primaryPdg);
            
        }

        fTotalTruthEnergy += total;
        
        if (primaryContributors.size() > 1) fSharedTruthEnergy += total;
            
    }
 
}

void DigitVectorTruthInfo::PrintMap() const
{
    std::cout << "Primary particle PDG | energy deposit " << std::endl;
    for (std::map<int,double>::const_iterator p = fPdgEdepMap.begin();
         p != fPdgEdepMap.end(); ++p) {
        std::cout << "\t--> " << p->first << " " << p->second << std::endl;
    }
}
