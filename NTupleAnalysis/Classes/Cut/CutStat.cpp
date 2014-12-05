#ifndef CutStat_cpp
#define CutStat_cpp

#include "CutStat.h"

using namespace std;

CutStat::CutStat()
{   
    count = 0;
}

void CutStat::increment()
{
    count++;   
}

double CutStat::getCount()
{
    return count;
}




#endif
