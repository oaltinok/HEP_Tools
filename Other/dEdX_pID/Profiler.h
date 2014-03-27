#ifndef PROFILER_H
#define PROFILER_H

#include<vector>
#include<cmath>
#include<string>
#include<iostream>

using namespace std;

// double const refDeltadEdX = 0.1;
string const debugName = "PROFILER: ";

class Profiler
{
    private:
        bool isForward; // True if the track has a forward direction
        char view;
        vector<double> dEdX;
        vector<double> posZ;
        vector<double> posZ_tobeRemoved;
        
        bool isMarked(double inputZ);
        void markAllBefore(unsigned int ind);
      
    
    
    
    public:
        Profiler();
        Profiler(char input_view);
        ~Profiler();
        void set_view(char input_view);
        char get_view();
        void add_dEdX(double input_dEdX);
        void add_posZ(double input_Z);
        void set_isForward(bool input_isForward);
        bool get_isForward();
        bool needRemoval(double inputZ);

        
        void findOutliers();
    

};


#endif

