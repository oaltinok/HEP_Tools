#ifndef CCProtonPi0_CutStat_cpp
#define CCProtonPi0_CutStat_cpp

#include "CCProtonPi0_CutStat.h"

using namespace std;

CCProtonPi0_CutStat::CCProtonPi0_CutStat()
{   
    count = 0;
}

void CCProtonPi0_CutStat::increment()
{
    count++;   
}

double CCProtonPi0_CutStat::getCount()
{
    return count;
}




#endif
