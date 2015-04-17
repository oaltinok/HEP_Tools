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

void Cut::increment(bool isSignal, bool isStudy1_true, bool isStudy2_true)
{
    nEvent.increment();
    if(isSignal) nSignal.increment();
    if(isStudy1_true) nStudy1.increment();
    if(isStudy2_true) nStudy2.increment();
}




#endif
