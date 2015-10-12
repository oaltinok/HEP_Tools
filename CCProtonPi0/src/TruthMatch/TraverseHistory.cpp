#undef NDEBUG
#include <cassert>
#include <iostream>
#include <algorithm>

#include <Event/TG4Trajectory.h>
#include "TraverseHistory.h"

TraverseHistory::TraverseHistory(int trackId)
    : fTrackId(trackId),
      fPdg(0) {}

void TraverseHistory::DoTraverse(const std::map<int,TG4Trajectory*>& trajectoryMap)
{
    int trackId = fTrackId;
    fHistoryTrackId.push_back(trackId);
    while (trackId != 0) {
        TG4Trajectory* traj = trajectoryMap.at(trackId);
        
        if (!traj) break;
            
        trackId = traj->GetParentId();
        fHistoryTrackId.push_back(trackId);
    }
    
    std::reverse(fHistoryTrackId.begin(),fHistoryTrackId.end());
        /*
          std::copy(history.begin(), history.end(),
          std::ostream_iterator<int>(std::cout, " "));
        */
    
        // 0: not defined, 1: final-state pi0, 2: gamma's
    assert(!fHistoryTrackId.empty());

    if (fHistoryTrackId.size() < 2) {
        std::cerr << "Invalid trackId: " << fTrackId << std::endl;
        return;
    }
        /* The first element in fHistoryTrackId is 0 which of course has
           no corresponding trajectory pointer, so skip it */
    fHistoryPdg.push_back(0);
    for (std::vector<int>::const_iterator p = fHistoryTrackId.begin() + 1;
         p != fHistoryTrackId.end(); ++p) {
        const TG4Trajectory* traj = trajectoryMap.at(*p);
        const int pdg = traj->GetPDGCode();
        fHistoryPdg.push_back(pdg);
    }

    fPdg = fHistoryPdg.at(1);
    
    assert(fHistoryPdg.size() == fHistoryTrackId.size());
}

int TraverseHistory::GetPDGCode() const {
    return fPdg;
}

int TraverseHistory::GetPrimaryId() const {

    if (fHistoryTrackId.size() < 2) {
        std::cerr << "Invalid trackId: " << fTrackId << std::endl;
        return fTrackId;
    }
    
    return fHistoryTrackId.at(1);
}

int TraverseHistory::GetPrimaryPdg() const {
    if (fHistoryPdg.size() < 2) {
        std::cerr << "Invalid trackId: " << fTrackId << std::endl;
        return 0;
    }
    
    return fHistoryPdg.at(1);
}


int TraverseHistory::GetNextToPrimaryId() const {

    if (fHistoryTrackId.size() < 3) {
        return fTrackId;
    }
    
    return fHistoryTrackId.at(2);
}

int TraverseHistory::GetNextToPrimaryPdg() const {
    if (fHistoryPdg.size() < 3) {
        return fPdg;
    }
    
    return fHistoryPdg.at(2);
}

bool TraverseHistory::IsPresent(int pdg) const {
    std::vector<int>::const_iterator target
        = std::find(fHistoryPdg.begin(),fHistoryPdg.end(),pdg);

    return target != fHistoryPdg.end();
}
