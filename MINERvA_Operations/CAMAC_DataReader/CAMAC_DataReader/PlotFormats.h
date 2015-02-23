#ifndef PlotFormats_h
#define PlotFormats_h

#include <vector>

using namespace std;

struct format_hist
{
    string title;
    string var_name;
    int nbins;
    double low;
    double high;
};

struct format_graph
{
    string title;
    string x_axis;
    string y_axis;
};

struct format_frequency
{
    string title;
    vector<int> vars;
};



#endif


