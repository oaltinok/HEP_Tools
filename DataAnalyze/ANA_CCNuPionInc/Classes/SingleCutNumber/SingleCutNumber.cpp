/*
    See SingleCutNumber.h header for Class Information
*/

#include "SingleCutNumber.h"

using namespace std;

SingleCutNumber::SingleCutNumber()
{

}


void SingleCutNumber::increment()
{
    value++;
}

void SingleCutNumber::setValue(double input)
{
    value = input;
}

void SingleCutNumber::setLabel(string input)
{
    label = input;
}

string SingleCutNumber::getLabel()
{
    return label;
}

double SingleCutNumber::getValue()
{
    return value;
}


