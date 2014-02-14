/*
    See CutNumber.h header for Class Information
*/

#include "CutNumber.h"

using namespace std;

CutNumber::CutNumber()
{
    setValue(initCount);
}

CutNumber::CutNumber(string inputLabel)
{
    setValue(initCount);
    setLabel(inputLabel);
}


void CutNumber::increment()
{
    value++;
}


void CutNumber::setValue(double input)
{
    value = input;
}

void CutNumber::setLabel(string input)
{
    label = input;
}

string CutNumber::getLabel()
{
    return label;
}

double CutNumber::getValue()
{
    return value;
}


