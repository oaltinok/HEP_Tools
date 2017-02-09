#include "../Classes/Analyzer/CCProtonPi0_Analyzer.h"
#include "../Classes/CrossSection/CCProtonPi0_CrossSection.h"
#include "../Classes/Plotter/CCProtonPi0_Plotter.h"
#include "../Classes/SideBandTool/CCProtonPi0_SideBandTool.h"

const string runOption_Run = "run";
const string runOption_Plot = "plot";
const string runOption_Reduce = "reduce";
const string runOption_CrossSection = "calc";
const string runOption_FitSideBand = "fit";
const string runOption_FitW = "fitW";
const string runOption_FitQSq = "fitQSq";

const string typeOption_mc = "mc";
const string typeOption_data = "data";

void Reduce(string playlist, bool isMC)
{
    bool isModeReduce = true;
    cout<<"\n"<<endl;
    cout<<"======================================================================"<<endl;
    cout<<"Reducing NTuples..."<<endl;
    cout<<"======================================================================"<<endl;
    CCProtonPi0_Analyzer t(isModeReduce, isMC);
    t.reduce(playlist);
}

void Analyze(string playlist, bool isMC)
{
    bool isModeReduce = false;
    cout<<"\n"<<endl;
    cout<<"======================================================================"<<endl;
    cout<<"Analyzing NTuples, Creating Histograms..."<<endl;
    cout<<"======================================================================"<<endl;
    CCProtonPi0_Analyzer analyzer(isModeReduce, isMC);
    analyzer.analyze(playlist);
}

void Calculate_CrossSection(bool isMC)
{
    cout<<"\n"<<endl;
    cout<<"======================================================================"<<endl;
    cout<<"Calculating Cross Section..."<<endl;
    cout<<"======================================================================"<<endl;
    CCProtonPi0_CrossSection crossSection(isMC);
    crossSection.Calc_CrossSections();
}

void Plot()
{
    cout<<"======================================================================"<<endl;
    cout<<"Plotting Histograms..."<<endl;
    cout<<"======================================================================"<<endl;
    CCProtonPi0_Plotter plotter;
    plotter.plotHistograms();
}

/*
 *  1   reduce
 *  2   run
 *  3   calculate cross section
 *  10  plot
 *  20  fit side band
 *   
*/
int GetMode(int argc, char* argv[])
{
    // argc can only be 2 or 3
    if (argc != 2 && argc != 3) return 0;

    std::string runSelect = argv[1];
    if (argc == 2){
        if (runSelect.compare(runOption_Plot) == 0) return 10;
        else if (runSelect.compare(runOption_FitSideBand) == 0) return 20;
        else if (runSelect.compare(runOption_FitW) == 0) return 30;
        else if (runSelect.compare(runOption_FitQSq) == 0) return 40;
        else return 0;
    }
     
    std::string typeSelect = argv[2];
    // First check for ERROR
    if (runSelect.compare(runOption_Reduce) != 0 && runSelect.compare(runOption_Run) != 0 && runSelect.compare(runOption_CrossSection) != 0) return 0;
    if (typeSelect.compare(typeOption_mc) != 0 && typeSelect.compare(typeOption_data) != 0) return 0;

    // Passed ERROR Check - Valid Input    
    if (runSelect.compare(runOption_Reduce) == 0){
        if (typeSelect.compare(typeOption_mc) == 0) return -1;
        else if (typeSelect.compare(typeOption_data) == 0) return 1;
        else return 0;
    }

    if (runSelect.compare(runOption_Run) == 0){
        if (typeSelect.compare(typeOption_mc) == 0) return -2;
        else if (typeSelect.compare(typeOption_data) == 0) return 2;
        else return 0;
    }

    if (runSelect.compare(runOption_CrossSection) == 0){
        if (typeSelect.compare(typeOption_mc) == 0) return -3;
        else if (typeSelect.compare(typeOption_data) == 0) return 3;
        else return 0;
    }

    return 0;
}

void showInputError(char *argv[])
{
    cout<<std::left;
    cout<<"----------------------------------------------------------------------"<<endl;
    cout<<"Not a valid syntax!"<<endl;
    cout<<"----------------------------------------------------------------------"<<endl;
    cout<<"Correct Syntax for NTuple Reduce"<<endl;
    cout<<"\t"<<argv[0]<<" "<<runOption_Reduce<<" "<<typeOption_mc<<endl;
    cout<<"\t"<<argv[0]<<" "<<runOption_Reduce<<" "<<typeOption_data<<"\n"<<endl;
    cout<<"Correct Syntax for NTuple Analysis"<<endl;
    cout<<"\t"<<argv[0]<<" "<<runOption_Run<<" "<<typeOption_mc<<endl;
    cout<<"\t"<<argv[0]<<" "<<runOption_Run<<" "<<typeOption_data<<"\n"<<endl;
    cout<<"Correct Syntax for Calculating Cross Section"<<endl;
    cout<<"\t"<<argv[0]<<" "<<runOption_CrossSection<<" "<<typeOption_mc<<endl;
    cout<<"\t"<<argv[0]<<" "<<runOption_CrossSection<<" "<<typeOption_data<<"\n"<<endl;
    cout<<"Correct Syntax for Plotting"<<endl;
    cout<<"\t"<<argv[0]<<" "<<runOption_Plot<<"\n"<<endl;
    cout<<"Correct Syntax for Fitting SideBands"<<endl;
    cout<<"\t"<<argv[0]<<" "<<runOption_FitSideBand<<"\n"<<endl;
    cout<<"----------------------------------------------------------------------"<<endl;
}


