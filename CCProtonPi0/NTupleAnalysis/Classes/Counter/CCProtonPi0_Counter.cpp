#ifndef CCProtonPi0_Counter_cpp
#define CCProtonPi0_Counter_cpp

#include "CCProtonPi0_Counter.h"

CCProtonPi0_Counter::CCProtonPi0_Counter()
{
    count = 0.0;
    name = "";
    isCounted = false;
}

void CCProtonPi0_Counter::setName(std::string inputName)
{
    name = inputName;
}

void CCProtonPi0_Counter::increment(double wgt)
{
    count = count + wgt;
}

void CCProtonPi0_Counter::print()
{
    std::cout<<name<<" = "<<count<<std::endl;
}

double CCProtonPi0_Counter::getCount()
{
    return count;
}

std::string CCProtonPi0_Counter::getName()
{
    return name;
}

double CCProtonPi0_Counter::calcPercent(double reference)
{
    double percent = (count / reference)*100;
    return percent;
}


#endif

