/*
================================================================================
Class: CCProtonPi0_SingleBin
    Member Variable for BinList Class
    
    Main Directory:
        Classes/CCProtonPi0_SingleBin        

    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision: 2015_05_07
================================================================================
*/

#ifndef CCProtonPi0_SingleBin_h
#define CCProtonPi0_SingleBin_h

class CCProtonPi0_SingleBin
{
    public:
        CCProtonPi0_SingleBin();
        
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
