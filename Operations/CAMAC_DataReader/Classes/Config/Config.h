#ifndef Config_h
#define Config_h

#include "../Reader/Reader.h"

using namespace std;


class Config : public Reader{
    
    public:
        Config();
        void ReadFile();


    private:
        void ReadHistogramFormat(string line);
        void ReadGraphFormat(string line);
        void ReadFrequencyFormat(string line);
        void PrintVector(vector<string>* v);
        void ReadConfigLine(string line, vector<string>* v);
        void ProcessConvFactors();
        
        bool isDebugging;

};


#endif
