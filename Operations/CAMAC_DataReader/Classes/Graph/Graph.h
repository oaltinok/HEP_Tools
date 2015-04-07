#ifndef Graph_h
#define Graph_h

#include "../Plotter/Plotter.h"

class Graph : public Plotter{
    
    public:
        void AddFormat(string in_title, string in_x_axis, string in_y_axis);
        void Init();
	void Reset();
        void Fill(string var_name, double value);
        void WriteRootFile();

    private:
        /*
            x vs. y Graph
        */
        struct Format
        {
            string title;
            string x_axis; // Name of the Variable 
            string y_axis; // Name of the Variable
        };

        vector < Format > formats;
        vector < TGraph* > graphs;
        vector < vector<double> > graph_x_axis_vectors;
        vector < vector<double> > graph_y_axis_vectors;

};

#endif

