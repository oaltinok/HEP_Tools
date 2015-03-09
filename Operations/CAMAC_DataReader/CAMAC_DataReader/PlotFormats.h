/*
================================================================================
Library: PlotFormats
    Includes Definitions of Plot Formats used in CAMAC_DataReader Package
    
    Usage:
        #include "CAMAC_DataReader/PlotFormats.h"
        format_hist = f_hist;
    
    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision:  2015_02_24
================================================================================
*/
#ifndef PlotFormats_h
#define PlotFormats_h

#include <vector>

using namespace std;

/*
    1D Histogram
*/
struct format_hist
{
    string title;
    string var_name;
    int nbins;
    double low;
    double high;
};

/*
    x vs. y Graph
*/
struct format_graph
{
    string title;
    string x_axis; // Name of the Variable 
    string y_axis; // Name of the Variable
};

/*
    Frequency Plot
*/
struct format_frequency
{
    string title;
    vector<int> vars; 
    // Array with Variable Info 
    //  if vars[ind] == 1 --> Check Frequency
    //  if vars[ind] == 0 --> Ignore Variable
};



#endif


