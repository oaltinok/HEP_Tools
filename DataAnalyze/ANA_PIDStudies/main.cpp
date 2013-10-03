#include "ANA_PIDStudies.cpp"
#include "Plotter.cpp"

using namespace std;

void runTestSample(int run_number);
void runSample(int run_number);
void runOptimization(bool withPlots);


int main(int run_number, bool withPlots)
{

// runTestSample(run_number);
// runSample(run_number, withPlots);
runOptimization(withPlots);

// return 0;

}

void runOptimization(bool withPlots)
{
    ANA_PIDStudies t;
    Plotter p;
    
    t.run("Playlists/Optimization/pl_run501_Optimization_0_10.dat", "RootFiles/Optimization_0_10.root", "TextFiles/Optimization_0_10.txt","readme.txt");
    
    if(withPlots){
        p.plotHistograms("RootFiles/Optimization_0_10.root","Plots/Optimization/dEdX_0_10");
    }

}

void runSample(int run_number, bool withPlots)
{
    
    ANA_PIDStudies t;
    Plotter p;
    
    if( run_number == 500){
        
        // Outliers Not Removed:
        cout<<"Running for Outliers_NotRemoved for run = "<<run_number<<endl;
        t.run("Playlists/pl_run500_Outliers_NotRemoved.dat", "RootFiles/run500_Outliers_NotRemoved.root", "TextFiles/run500_cutFile_No.txt","readme.txt");
        
        // Outliers Tammy:
        cout<<"Running for Outliers_Tammy for run = "<<run_number<<endl;
        t.run("Playlists/pl_run500_Outliers_Tammy.dat", "RootFiles/run500_Outliers_Tammy.root", "TextFiles/run500_cutFile_Tammy.txt","readme.txt");
        
        // Outliers Brandon:
        cout<<"Running for Outliers_Brandon for run = "<<run_number<<endl;
        t.run("Playlists/pl_run500_Outliers_Brandon.dat", "RootFiles/run500_Outliers_Brandon.root", "TextFiles/run500_cutFile_Brandon.txt","readme.txt");
        
        // Outliers Ozgur:
        cout<<"Running for Outliers_Ozgur for run = "<<run_number<<endl;
        t.run("Playlists/pl_run500_Outliers_Ozgur.dat", "RootFiles/run500_Outliers_Ozgur.root", "TextFiles/run500_cutFile_Ozgur.txt","readme.txt");
        
        if(withPlots){
            p.plotHistograms("RootFiles/run500_Outliers_NotRemoved.root","Plots/run500/Outliers_NotRemoved");
            p.plotHistograms("RootFiles/run500_Outliers_Tammy.root","Plots/run500/Outliers_Tammy");
            p.plotHistograms("RootFiles/run500_Outliers_Brandon.root","Plots/run500/Outliers_Brandon");
            p.plotHistograms("RootFiles/run500_Outliers_Ozgur.root","Plots/run500/Outliers_Ozgur");
        }
    }    
    
    if( run_number == 501){
        
        // Outliers Not Removed:
        cout<<"Running for Outliers_NotRemoved for run = "<<run_number<<endl;
        t.run("Playlists/pl_run501_Outliers_NotRemoved.dat", "RootFiles/run501_Outliers_NotRemoved.root", "TextFiles/run501_cutFile_No.txt","readme.txt");
        
        // Outliers Tammy:
        cout<<"Running for Outliers_Tammy for run = "<<run_number<<endl;
        t.run("Playlists/pl_run501_Outliers_Tammy.dat", "RootFiles/run501_Outliers_Tammy.root", "TextFiles/run501_cutFile_Tammy.txt","readme.txt");
        
        // Outliers Brandon:
        cout<<"Running for Outliers_Brandon for run = "<<run_number<<endl;
        t.run("Playlists/pl_run501_Outliers_Brandon.dat", "RootFiles/run501_Outliers_Brandon.root", "TextFiles/run501_cutFile_Brandon.txt","readme.txt");
        
        // Outliers Ozgur:
        cout<<"Running for Outliers_Ozgur for run = "<<run_number<<endl;
        t.run("Playlists/pl_run501_Outliers_Ozgur.dat", "RootFiles/run501_Outliers_Ozgur.root", "TextFiles/run501_cutFile_Ozgur.txt","readme.txt");
        
        if(withPlots){
            p.plotHistograms("RootFiles/run501_Outliers_NotRemoved.root","Plots/run501/Outliers_NotRemoved");
            p.plotHistograms("RootFiles/run501_Outliers_Tammy.root","Plots/run501/Outliers_Tammy");
            p.plotHistograms("RootFiles/run501_Outliers_Brandon.root","Plots/run501/Outliers_Brandon");
            p.plotHistograms("RootFiles/run501_Outliers_Ozgur.root","Plots/run501/Outliers_Ozgur");
        }
    }    

}

void runTestSample(int run_number)
{
    
    ANA_PIDStudies t;
    Plotter p;
    
    if( run_number == 500){
    
            // Outliers Ozgur:
        cout<<"Running for Outliers_Ozgur for run = "<<run_number<<endl;
        t.run("Playlists/Test/pl_run500_Outliers_Ozgur.dat", "run500_Outliers_Ozgur.root", "cutFile.txt","readme.txt");
        p.plotHistograms("run500_Outliers_Ozgur.root","Plots/run500/Outliers_Ozgur");

//         // Outliers Ozgur2:
//         cout<<"Running for Outliers_Ozgur for run = "<<run_number<<endl;
//         t.run("Playlists/Test/pl_run500_Outliers_Ozgur2.dat", "run500_Outliers_Ozgur2.root", "cutFile.txt","readme.txt");
//         p.plotHistograms("run500_Outliers_Ozgur2.root","Plots/run500/Outliers_Ozgur2");
    }    
    
    if( run_number == 501){
        
        // Outliers Ozgur:
        cout<<"Running for Outliers_Ozgur for run = "<<run_number<<endl;
        t.run("Playlists/Test/pl_run501_Outliers_Ozgur.dat", "run501_Outliers_Ozgur.root", "cutFile.txt","readme.txt");
        p.plotHistograms("run501_Outliers_Ozgur.root","Plots/run501/Outliers_Ozgur");
    }    

}
