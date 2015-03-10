#ifndef Freq_h
#define Freq_h

#include "../Plotter/Plotter.h"

class Freq : public Plotter{

        public:
        void AddFormat(string in_title, vector<int> in_vars);
        void Init();
        void Fill(int ind);
        void WriteRootFile();

    private:
        /*
            Frequency Plot
        */
        struct Format
        {
            string title;
            vector<int> vars; 
            // Array with Variable Info 
            //  if vars[ind] == 1 --> Check Frequency
            //  if vars[ind] == 0 --> Ignore Variable
        };
        
        vector < Format > formats;
        vector < TGraph* > freq_plots;
        vector < vector<double> > freq_x_axis_vectors;
        vector < vector<double> > freq_y_axis_vectors;

};


#endif
