#ifndef Data_h
#define Data_h

#include "../Reader/Reader.h"

class Data : public Reader{
    
    public:
        Data();
        void ReadFile();
	string GetFirstTime(string var_name);

    private:
        void ReadRunSubrun(string line);
        void ReadDataLine(string line);
        void ProcessDataLine(string var_name, double value);
        int GetVariableInd(string var_name);
        void CheckRunSubrun();
        void InitPlots();

        
        bool isDebugging;
        bool newFile;
        bool resetPlots;
        
        string run;
        string subrun;
        string time;

};


#endif
