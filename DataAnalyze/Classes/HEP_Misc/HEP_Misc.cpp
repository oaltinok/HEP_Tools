/*
    See HEP_Misc.h header for Class Information
*/
#include "HEP_Misc.h"

using namespace std;

HEP_Misc::HEP_Misc()
{

}


double HEP_Misc::getError(double trueValue, double recoValue)
{

    double error;
    
    error = (trueValue - recoValue) / trueValue;
    
    return error;

}






