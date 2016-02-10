#ifndef Constants_h_seen
#define Constants_h_seen

#include <iostream>
#include <iomanip>

using std::setw;

class Constants {
  public:
    static Constants& get();
    
    double m0;
    double mp;
    double mn;
    double mN;
    double mm;
    
    int n_universe;
    double emax;
    double emin;
    double wmin;
    double wmax;
    
    void Print();
    
  private:
    Constants();
    Constants(const Constants&) {}
    Constants& operator=(const Constants&) { return *this;}
    
    ~Constants() {}

};

Constants::Constants()
{
    m0 = 134.9766;
    mp = 938.272046; // MeV!
    mn = 939.565379; // MeV!
    mN = (mp + mn)/2.0;
    mm = 105.659;    // MeV! same as Carrie
    
    n_universe = 100;

    emax = 20.0;
    emin =  1.5;

    wmin = 0.0;
    wmax = 1.8;
}


Constants& Constants::get()
{
    static Constants singleton;

    return singleton;
}

Constants& constants()
{
    return Constants::get();
}


void Constants::Print()
{
    std::cout << "-------------PARAMETERS-------------------" << std::endl;
    std::cout << setw(20) << "mgg"   << setw(10) << m0  << std::endl;
    
    std::cout << "------------------------------------------"   << std::endl;
}

#endif
