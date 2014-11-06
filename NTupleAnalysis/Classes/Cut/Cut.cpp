/*
    See Cut.h header for Class Information
*/
#ifndef Cut_cpp
#define Cut_cpp

#include "Cut.h"

using namespace std;

Cut::Cut()
{
	nEvent = 0;
	nSignal = 0;
}

void Cut::set_Name(string inputName)
{
    name = inputName;
}

string Cut::get_Name()
{
	return name;
}

double Cut::get_nEvent()
{
	return nEvent;
}

double Cut::get_nSignal()
{
	return nSignal;
}

double Cut::get_nSignal_Gold()
{
    return nSignal_Gold;
}

double Cut::get_nSignal_Silver1()
{
    return nSignal_Silver1;
}

void Cut::inc_nEvent()
{
    nEvent++;
}

void Cut::inc_nSignal()
{
    nSignal++;
}

void Cut::inc_nSignal_Gold()
{
    nSignal_Gold++;
}

void Cut::inc_nSignal_Silver1()
{
    nSignal_Silver1++;
}


#endif
