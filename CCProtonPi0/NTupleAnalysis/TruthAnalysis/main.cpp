#include <iostream>
#include <ctime>
#include <string>
#include "Cintex/Cintex.h"

#include "../Classes/NTupleAnalysis/CCProtonPi0_NTupleAnalysis.h"
#include "Classes/TruthAnalyzer/CCProtonPi0_TruthAnalyzer.h"

using namespace std;

int main()
{
    time_t timeStart; time(&timeStart);
    time_t timeEnd;
    double timeDiff;
    int timeDiff_m;
    int timeDiff_s;

    TH1::AddDirectory(false);
    
    string playlist = "../Input/Playlists/pl_MC_Merged.dat";
    CCProtonPi0_TruthAnalyzer t;
    t.Loop(playlist);
    
    // Process Time 
    time(&timeEnd);
    timeDiff = ( timeEnd - timeStart );
    
    timeDiff_m = (int)timeDiff / 60;
    timeDiff_s = (int)timeDiff % 60;
    
    cout<<"Time to Complete = "<<timeDiff_m<<"\' "<<timeDiff_s<<"\" or "<<timeDiff<<" seconds"<<endl;

    return 0;
}

