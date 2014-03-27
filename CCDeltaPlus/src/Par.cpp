#include <iostream>
#include <iomanip>


#include "Par.h"


Par& Par::Get()
{
    static Par singleton;

    return singleton;
}

Par& Constants()
{
    return Par::Get();
}

using std::setw;
void Par::Print()
{
    std::cout << "-------------PARAMETERS-------------------" << std::endl;
    std::cout << setw(20) << "k_ecal"   << setw(10) << k_ecal  << std::endl;
    std::cout << setw(20) << "k_trkr"   << setw(10) << k_trkr  << std::endl;
    std::cout << setw(20) << "coangle"  << setw(10) << coangle << std::endl;
    std::cout << setw(20) << "clength"  << setw(10) << clength << std::endl;
    std::cout << setw(20) << "xminevis" << setw(10) << xminevis << std::endl;
    std::cout << setw(20) << "uminevis" << setw(10) << uminevis << std::endl;
    std::cout << setw(20) << "vminevis" << setw(10) << vminevis << std::endl;
    std::cout << "------------------------------------------"   << std::endl;
}
