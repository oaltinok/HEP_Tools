/*
================================================================================
Class: SingleBin
    Member Variable for BinList Class
    
    Main Directory:
        Classes/SingleBin
    
    Usage:
        > #include "Classes/SingleBin/SingleBin.h" 
        > SingleBin hadronMomentum
        > hadronMomentum.set_nBins(40);
        > hadronMomentum.get_nBins();
            
    
    Last Revision: 2014_01_30
================================================================================
*/

#ifndef SingleBin_h
#define SingleBin_h

class SingleBin
{
    public:
        SingleBin();
        
        int get_nBins();
        double get_max();
        double get_min();
        double get_width();
        
        void setBin(int in_nBins, double in_min, double in_max);
        void set_nBins(int input);
        void set_max(double input);
        void set_min(double input);

    private:
        int nBins;
        double max;
        double min;
        double width;

};



#endif
