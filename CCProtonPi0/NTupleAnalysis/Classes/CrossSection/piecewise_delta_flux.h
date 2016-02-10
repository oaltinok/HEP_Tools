#ifndef piecewise_delta_flux_h_seen
#define piecewise_delta_flux_h_seen

#include <iostream>
#include <iomanip>


class piecewise_delta_flux {
  public:
    static piecewise_delta_flux& get();

    void set_data_file(const std::string& file);
    double calculate(double enu0);
    
  private:
    piecewise_delta_flux();
    piecewise_delta_flux(const piecewise_delta_flux&) {}
    piecewise_delta_flux& operator=(const piecewise_delta_flux&) { return *this;}
    
    ~piecewise_delta_flux() {}


};

piecewise_delta_flux::piecewise_delta_flux()
{}


piecewise_delta_flux& piecewise_delta_flux::get()
{
    static piecewise_delta_flux singleton;

    return singleton;
}

double piecewise_delta_flux::calculate(double enu0)
{
    const double shifts[] = {-0.1, -0.1, -0.1, -0.1, -0.2, -0.2};

    double delta = 0.0;
    
    if (0.0 < enu0 && enu0 <= 1.5)  delta = shifts[0]/2.0;
    if (1.5 < enu0 && enu0 <= 2.5)  delta = shifts[1]/2.0;
    if (2.5 < enu0 && enu0 <= 4.0)  delta = shifts[2]/2.0;
    if (4.0 < enu0 && enu0 <= 5.5)  delta = shifts[3]/2.0;
    if (5.5 < enu0 && enu0 <= 7.5)  delta = shifts[4]/2.0;
    if (7.5 < enu0 && enu0 <= 10.0) delta = shifts[5]/2.0;
    
    return delta;
    
}

#endif
