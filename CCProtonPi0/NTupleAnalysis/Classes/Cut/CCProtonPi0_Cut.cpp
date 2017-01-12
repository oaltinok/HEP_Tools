/*
    See CCProtonPi0_Cut.h header for Class Information
*/
#ifndef CCProtonPi0_Cut_cpp
#define CCProtonPi0_Cut_cpp

#include "CCProtonPi0_Cut.h"

using namespace std;

CCProtonPi0_Cut::CCProtonPi0_Cut()
{
    // Do Nothing
}

void CCProtonPi0_Cut::set_Name(string inputName)
{
    name = inputName;
}

string CCProtonPi0_Cut::get_Name()
{
	return name;
}

void CCProtonPi0_Cut::increment(bool isSignal, bool isStudy1_true, bool isStudy2_true, double wgt)
{
    nEvent.increment(wgt);
    if(isSignal) nSignal.increment(wgt);
    if(isStudy1_true) nStudy1.increment(wgt);
    if(isStudy2_true) nStudy2.increment(wgt);
}




#endif
