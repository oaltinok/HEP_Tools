#ifndef Hist_h
#define Hist_h

#include "../Plotter/Plotter.h"

using namespace std;

class Hist : public Plotter{
    
    public:
        void AddFormat(   string in_title, string in_var_name,
                        int in_nbins, double in_low, double in_high);
        void Init();
        void Fill(string var_name, double value);
        void WriteRootFile();

    private:
        /*
            1D Histogram
        */
        struct Format
        {
            string title;
            string var_name;
            int nbins;
            double low;
            double high;
        };
        
        vector<Format> formats;
        vector< TH1D* > hists;

};


#endif
