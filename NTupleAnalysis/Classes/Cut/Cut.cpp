/*
    See Cut.h header for Class Information
*/
#ifndef Cut_cpp
#define Cut_cpp

#include "Cut.h"

using namespace std;

Cut::Cut()
{
    // Do Nothing
}

void Cut::set_Name(string inputName)
{
    name = inputName;
}

string Cut::get_Name()
{
	return name;
}

void Cut::increment(bool isSignal, bool isGold, bool isSilver1, bool isDIS)
{
    nEvent.increment();
    if (isSignal) nSignal.increment();
    if (isGold) nSignal_Gold.increment();
    if (isSilver1) nSignal_Silver1.increment();
    if (isDIS) nDIS.increment();    
}




#endif
