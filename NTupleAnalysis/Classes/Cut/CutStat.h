#ifndef CutStat_h
#define CutStat_h

#include <string>

using namespace std;

class CutStat
{
    public:
        CutStat();
        
        void increment();
        double getCount();
        
    private:
        double count;
        
};



#endif
