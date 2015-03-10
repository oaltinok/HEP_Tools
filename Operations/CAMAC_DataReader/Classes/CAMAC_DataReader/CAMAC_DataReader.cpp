#ifndef CAMAC_DataReader_cpp
#define CAMAC_DataReader_cpp

#include "CAMAC_DataReader.h"

using namespace std;

/*
    Default Constructor
        NOT USED!!
*/
CAMAC_DataReader::CAMAC_DataReader()
{
    ErrorMode();
}

/*
    Constructor for Auto Mode
*/
CAMAC_DataReader::CAMAC_DataReader(string mode)
{
    cout<<"\n\n";
    cout<<"#####################################################################"<<endl;
    cout<<"#                                                                   #"<<endl;
    cout<<"#                    Welcome CAMAC_DataReader                       #"<<endl;
    cout<<"#                                                                   #"<<endl;
    cout<<"#####################################################################"<<endl;
    cout<<"\n\n";
    
    if( mode.compare("Auto") == 0){
        cout<<"Mode = Auto"<<endl;
        isModeAuto = true;
    }else{
        ErrorMode();   
    }
    
    InitParamaters();
    
    // Set File Locations
    SetConfigDir();
    SetDataDir();
    SetRootDir();
    
    // Open Root File - Required for Plotter Initialization
    OpenRootFile();
    
    InitReaders();
    InitPlotters();
    
    // Read Config File - Reads the configuration only once
    ReadConfigFile();
    
    cout<<"Initialization Successful!"<<endl;
    cout<<"\n\n";
}

/*
    Constructor for Manual Mode
*/
CAMAC_DataReader::CAMAC_DataReader(string mode, string run_number, string subrun_number)
{
    cout<<"\n\n";
    cout<<"#####################################################################"<<endl;
    cout<<"#                                                                   #"<<endl;
    cout<<"#                    Welcome CAMAC_DataReader                       #"<<endl;
    cout<<"#                                                                   #"<<endl;
    cout<<"#####################################################################"<<endl;
    cout<<"\n\n";
    
    if( mode.compare("Manual") == 0){
        cout<<"Mode = Manual"<<endl;
        isModeAuto = false;
    }else{
        ErrorMode();   
    }
    
    InitParamaters();
    
    // Set Run & Subrun
    run = run_number;
    subrun = subrun_number;
    
    // Set File Locations
    SetConfigDir();
    SetDataDir();
    SetRootDir();
    
    // Open Root File - Required for Plotter Initialization
    OpenRootFile();
    
    InitReaders();
    InitPlotters();
    
    // Read Config File - Reads the configuration only once
    ReadConfigFile();
    
    cout<<"Initialization Successful!"<<endl;
    cout<<"\n\n";
}

void CAMAC_DataReader::InitParamaters()
{
    // Initialize to empty strings
    latest_run = "";
    latest_subrun = "";
    latest_time = "";   
    
    newFile = false;
    resetPlots = false;
}

/*
    Read & Process Data File
    Write Plots inside a Root File
*/
void CAMAC_DataReader::RunManual()
{

    InitPlots();
    DataReader.OpenFile(dataDir);
    DataReader.ReadFile();
    WriteRootFile();

}

void CAMAC_DataReader::RunAuto()
{
    // Do Nothing!
    cout<<"Auto Mode did not implemented yet!"<<endl;
    exit(EXIT_FAILURE);
}

void CAMAC_DataReader::InitPlots()
{
    HistPlotter.Init();
    GraphPlotter.Init();
    FreqPlotter.Init();
}

void CAMAC_DataReader::SetDataDir()
{
    string temp;
    
    if(isModeAuto){
        dataDir = dataDir_Auto;  
    }else{
        temp = Get_dataDir_Manual();
        if( temp.size() == 0){
            cout<<"No Input File for run = "<<run<<" subrun = "<<subrun<<endl;
            exit(EXIT_FAILURE);
        }else{
            dataDir = temp;
        }
    }
    
    cout<<"Input Data File = "<<dataDir<<endl;
}

string CAMAC_DataReader::Get_dataDir_Manual()
{
    ifstream commandFile;
    string fileName;
    string command = "ls -t /work/data/copied/TB*" + run + "_*" + subrun + "_cosmc" + "*_camac.dat" + " >>command.dat";
    
    // Use system command and write output to command.dat
    system(command.c_str());
    
    // Open Command.dat File
    commandFile.open("command.dat");
    
    // Read Command.dat File
    //      First Line is the requested file
    getline(commandFile,fileName);
    
    // Remove command.dat
    system("rm -f command.dat");
    
    return fileName;
}

void CAMAC_DataReader::SetRootDir()
{
    if( isModeAuto){
        rootDir = rootDir_Auto;
    }else{
        rootDir = Get_rootDir_Manual();
    }
    
    cout<<"Output ROOT File = "<<rootDir<<endl;
}

string CAMAC_DataReader::Get_rootDir_Manual()
{
    string tempRootDir;
    string tempRun;
    string tempSubrun;
    int nZeros_inRun = 8 - run.size();
    int nZeros_inSubrun = 4 - subrun.size();
    
    // Check Run Number Format
    if(nZeros_inRun < 0 || nZeros_inRun == 8 ){
        cout<<"ERROR: Did you set run number correctly?"<<endl;
        exit(EXIT_FAILURE);
    }
    
    // Check Subrun Number Format
    if(nZeros_inSubrun < 0 || nZeros_inSubrun == 4 ){
        cout<<"ERROR: Did you set subrun number correctly?"<<endl;
        exit(EXIT_FAILURE);;
    }
    
    // Fill tempRun for File Name
    for(int i = 0; i < nZeros_inRun; i++){  
        tempRun.push_back('0');
    }
    for( unsigned int i = 0; i < run.size(); i++){
        tempRun.push_back(run[i]);
    }
    
    // Fill tempSubrun for File Name
    for(int i = 0; i < nZeros_inSubrun; i++){   
        tempSubrun.push_back('0');
    }
    for( unsigned int i = 0; i < subrun.size(); i++){
        tempSubrun.push_back(subrun[i]);
    }
    
    // Check Run/Subrun Number Formats
    if(tempRun.size() != 8 || tempSubrun.size() != 4){
        cout<<"ERROR: Did you set run/subrun numbers correctly?"<<endl;
        exit(EXIT_FAILURE);
    }
    
    // Form Root File Dir       
    tempRootDir =   "/minerva/data/testbeam2/camacdata/TB_" + 
                    tempRun + "_" +
                    tempSubrun +
                    "_camac.root"; 
                
    return tempRootDir;
}

void CAMAC_DataReader::CheckDataFile()
{
    string line;
    
    // Open Data File
    dataFile.open(dataDir.c_str());
    if (!dataFile.is_open()){
        cerr<<"Cannot Open File = "<<dataDir<<endl;
        exit(EXIT_FAILURE);    
    }

    // Read First Line - run subrun time 
    getline (dataFile,line);
    
    stringstream line_stream(line);
    
    line_stream >> run >> subrun >> time;
    
    CheckRunSubrun();
    
    dataFile.close();
}

void CAMAC_DataReader::CheckRunSubrun()
{
    // Run Changed??
    if( latest_run.compare(run) != 0){ 
        newFile = true; 
        resetPlots = true;          
    }
    // Subrun Changed??
    else if( latest_subrun.compare(subrun) != 0 ){ 
        newFile = true;
        resetPlots = true;    
    }
    // Time Changed??
    else if( latest_time.compare(time) !=0 ){                          
        resetPlots = false;
        newFile = true;                      
    }
    // Nothing Changed!
    else{                                                              
        resetPlots = false;
        newFile = false;
    }
 
    // If it is a new file update latest_* variables
    if( newFile ){
        latest_run = run;
        latest_subrun = subrun;
        latest_time = time;
    }
    
}


void CAMAC_DataReader::ReadConfigFile()
{
    ConfigReader.ReadFile();   
}

void CAMAC_DataReader::ReadDataFile()
{
    DataReader.OpenFile(dataDir);
    DataReader.ReadFile();   
}

void CAMAC_DataReader::SetConfigDir()
{
    ConfigReader.OpenFile(configDir);
    cout<<"Config File = "<<configDir<<endl;
}

void CAMAC_DataReader::OpenRootFile()
{
    // Open ROOT File
    f_Root = new TFile(rootDir.c_str(),"UPDATE");
    
    // Add Directories to ROOT File (if there is no directory)
    if(f_Root->Get(branchName.c_str()) == NULL) f_Root->mkdir(branchName.c_str());
    
    // Change Directory to Corresponding Branch
    f_Root->cd(branchName.c_str());
}


void CAMAC_DataReader::InitPlotters()
{
    HistPlotter.SetRootFile(f_Root);
    GraphPlotter.SetRootFile(f_Root);
    FreqPlotter.SetRootFile(f_Root);
}


void CAMAC_DataReader::InitReaders()
{
    ConfigReader.SetVectors(&version, &vars, &convFactors_str, &convFactors, &units);
    ConfigReader.SetPlots(&HistPlotter, &GraphPlotter, &FreqPlotter);
    
    DataReader.SetVectors(&version, &vars, &convFactors_str ,&convFactors, &units);
    DataReader.SetPlots(&HistPlotter, &GraphPlotter, &FreqPlotter);
}

void CAMAC_DataReader::WriteRootFile()
{
    HistPlotter.WriteRootFile();
    GraphPlotter.WriteRootFile();
    FreqPlotter.WriteRootFile();
    cout<<"Sucessfully Updated "<<rootDir<<endl;
}


CAMAC_DataReader::~CAMAC_DataReader()
{
    f_Root->Close();   
}

void CAMAC_DataReader::ErrorMode()
{
    cout<<"Wrong Mode Specification.. Correct Declaration as follows:"<<endl;
    cout<<"\tCAMAC_Data Data(\"Auto\")"<<endl;
    cout<<"\tCAMAC_Data Data(\"Manual\", <run_number>, <subrun_number>)"<<endl;
    exit(EXIT_FAILURE);
}
#endif
