#ifndef CCProtonPi0_CutStat_h
#define CCProtonPi0_CutStat_h

#include <string>

using namespace std;

class CCProtonPi0_CutStat
{
    public:
        CCProtonPi0_CutStat();
        
        void increment();
        double getCount();
        
    private:
        double count;
        
};



#endif
