#ifndef smooth_delta_flux_h_seen
#define smooth_delta_flux_h_seen

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>

#include <TGraph.h>


class smooth_delta_flux {
  public:
    static smooth_delta_flux& get();

    void set_data_file(const std::string& file);
    double calculate(double enu0);
    
  private:
    smooth_delta_flux();
    smooth_delta_flux(const smooth_delta_flux&) {}
    smooth_delta_flux& operator=(const smooth_delta_flux&) { return *this;}

    void read_flux_delta_data(TGraph* g);
    
    ~smooth_delta_flux() {}

    std::string __data_file;
    TGraph* __delta_flux_curve;
};

smooth_delta_flux::smooth_delta_flux()
: __delta_flux_curve(0)
{

    __delta_flux_curve = new TGraph;
    read_flux_delta_data(__delta_flux_curve);

}


smooth_delta_flux& smooth_delta_flux::get()
{
    static smooth_delta_flux singleton;

    return singleton;
}

double smooth_delta_flux::calculate(double enu0)
{

    double emin =  0.0;
    double emax = 10.0;

    if (enu0 < emin || enu0 > emax) return 0.0;
    
        // construct a TSpline3 from the graph and interpolate
        // using the spline
    return __delta_flux_curve->Eval(enu0, 0, "S");
}

void smooth_delta_flux::set_data_file(const std::string& file)
{
    __data_file = file;
}

void smooth_delta_flux::read_flux_delta_data(TGraph* g)    
{
        // read from a text file
        /*   
    if (__data_file.empty()) {
        std::cerr << "Data points for the delta flux is not set" << std::endl;
        exit(1);
    }
    
    std::ifstream ifs(__data_file.c_str());

    while (!ifs.eof()) {
        int index;
        double xval;
        double yval;

        ifs >> index >> xval >> yval;

        if (ifs.eof()) break;

        g->SetPoint(index,xval,yval);
        
    }
        */

        // hard-wire here for the momentum
    g->SetPoint(    0,   0.0000,  -0.0700);
    g->SetPoint(    1,   0.2500,  -0.0741);
    g->SetPoint(    2,   0.5000,  -0.0800);
    g->SetPoint(    3,   0.7500,  -0.0884);
    g->SetPoint(    4,   1.0000,  -0.1000);
    g->SetPoint(    5,   1.2500,  -0.1149);
    g->SetPoint(    6,   1.5000,  -0.1300);
    g->SetPoint(    7,   1.7500,  -0.1419);
    g->SetPoint(    8,   2.0000,  -0.1500);
    g->SetPoint(    9,   2.2500,  -0.1536);
    g->SetPoint(   10,   2.5000,  -0.1500);
    g->SetPoint(   11,   2.7500,  -0.1373);
    g->SetPoint(   12,   3.0000,  -0.1200);
    g->SetPoint(   13,   3.2500,  -0.1035);
    g->SetPoint(   14,   3.5000,  -0.0900);
    g->SetPoint(   15,   3.7500,  -0.0813);
    g->SetPoint(   16,   4.0000,  -0.0800);
    g->SetPoint(   17,   4.2500,  -0.0876);
    g->SetPoint(   18,   4.5000,  -0.1000);
    g->SetPoint(   19,   4.7500,  -0.1134);
    g->SetPoint(   20,   5.0000,  -0.1300);
    g->SetPoint(   21,   5.2500,  -0.1513);
    g->SetPoint(   22,   5.5000,  -0.1700);
    g->SetPoint(   23,   5.7500,  -0.1789);
    g->SetPoint(   24,   6.0000,  -0.1800);
    g->SetPoint(   25,   6.2500,  -0.1767);
    g->SetPoint(   26,   6.5000,  -0.1700);
    g->SetPoint(   27,   6.7500,  -0.1604);
    g->SetPoint(   28,   7.0000,  -0.1500);
    g->SetPoint(   29,   7.2500,  -0.1404);
    g->SetPoint(   30,   7.5000,  -0.1300);
    g->SetPoint(   31,   7.7500,  -0.1181);
    g->SetPoint(   32,   8.0000,  -0.1100);
    g->SetPoint(   33,   8.2500,  -0.1111);
    g->SetPoint(   34,   8.5000,  -0.1200);
    g->SetPoint(   35,   8.7500,  -0.1338);
    g->SetPoint(   36,   9.0000,  -0.1500);
    g->SetPoint(   37,   9.2500,  -0.1662);
    g->SetPoint(   38,   9.5000,  -0.1800);
    g->SetPoint(   39,   9.7500,  -0.1888);
    g->SetPoint(   40,  10.0000,  -0.1900);
    
}

#endif
