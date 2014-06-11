/*
    TraverseHistory.h Duplicated from CCPi0 Package on 2014-06-10
        Purpose: Make CCProtonPi0 Package independent of CCPi0 Package
        Future: Common Tools and Functions will be combined under AnaUtils or
                PionUtils
                
    Original Author:    Trung Le
    Author:             Ozgur Altinok  - ozgur.altinok@tufts.edu
    Date:               2014_06_10
    Last Revision:      2014_06_10
*/

#ifndef TraverseHistory_h_seen
#define TraverseHistory_h_seen

#include <map>
#include <vector>

namespace Minerva {
    class TG4Trajectory;
}

using Minerva::TG4Trajectory;

class TraverseHistory {
  public:
    explicit TraverseHistory(int trackId);

    void DoTraverse(const std::map<int, TG4Trajectory*>& trajectoryMap);

    int GetPDGCode() const;
    int GetPrimaryId() const;
    int GetPrimaryPdg() const;
    int GetNextToPrimaryId() const;
    int GetNextToPrimaryPdg() const;
    bool IsPresent(int pdg) const;
    
  private:
    int fTrackId;
    int fPdg;
    std::vector<int> fHistoryTrackId;
    std::vector<int> fHistoryPdg;
};


#endif
