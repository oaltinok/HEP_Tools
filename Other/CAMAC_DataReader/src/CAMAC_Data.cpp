/*
    See CAMAC_Data.h header or Class Information
*/
#ifndef CAMAC_Data_cpp
#define CAMAC_Data_cpp

#include "../CAMAC_DataReader/CAMAC_Data.h"

using namespace std;

CAMAC_Data::CAMAC_Data()
{
    
    // Bins for the Histograms
    TOF_nBins = 40;
    TOF_min = 0.0;
    TOF_max = 4000.0;

}


void CAMAC_Data::ProcessData()
{
    OpenFiles();
    ReadDataFile();
    WriteRootFile();
}


/*
    readFile()
        Reads File Line by Line and sends each line to the following functions:
            1st Line to ReadHeader()
            Other Lines to ReadScalar()
*/
void CAMAC_Data::ReadDataFile()
{
    bool isFirstLine = true;
    unsigned int nVars = 0;
    unsigned int currentVar = 0;
    string line;
    while (!dataFile.eof() )
    {
        getline (dataFile,line);
        /*
            If it is the First Line Call ReadHeader()
            else read scalar in order
        */
        if ( isFirstLine ){  
            isFirstLine = false;
            ReadHeader(line);
            nVars = vars.size();
        }else{
            // Reset currentVar Indice when it is equal to nVars
            if ( currentVar == nVars) currentVar = 0; 
                
            ReadScalar(line,currentVar);
            currentVar++;
        }
    }   
}

/*
    ReadHeader()
        Reads Header - 1st Line in the Data File
        Header contains information about the variables
            Each Word in the Header Line is a variable
            Saves each variable in vars vector
*/
void CAMAC_Data::ReadHeader(string header)
{
    string variable;
    char tempChar;
    
    // Debugging - Check Header Line
//     cout<<"Header = "<<header<<endl;
    
    /*
        Go Over Each Character and form Single Variable Name
            Each Variable Name separated by a tab character '\t'
        Then Add Each single variable to Variable Vector (vars)
    */
    for(unsigned int i = 0; i < header.size(); i++){
        tempChar = header[i];
        if(tempChar != '\t'){
            variable.push_back(tempChar);
        }else{
            vars.push_back(variable);
            variable.clear();
        }
    }
    
    // Initialize Histograms after reading Header File
    initHistograms();
    
    // Debugging
//     PrintVariables();
}

/*
    initHistograms()
        Initializes Histograms inside hists vector
        Uses vars vector to get Histogram Title
*/
void CAMAC_Data::initHistograms()
{
    TH1D* tempHist;
    
    // Sanity Check
    if ( vars.size() == 0 ){
        cout<<"ERROR: No Variable Names. Did You Read Header Correctly?"<<endl;
        exit(EXIT_FAILURE);
    }
    
    f_Root->cd(branchName.c_str());
    // Get Existing Histogram or Create a NEW one
    for (unsigned int i = 0; i < vars.size(); i++){
        tempHist = (TH1D*)gDirectory->Get(vars[i].c_str());
        if(tempHist == NULL){
            tempHist = new TH1D( vars[i].c_str(),vars[i].c_str(),TOF_nBins, TOF_min, TOF_max );
            tempHist->GetXaxis()->SetTitle(vars[i].c_str());
            tempHist->GetYaxis()->SetTitle("N(Events)");
        }
        hists.push_back(tempHist);
    }
}


/*
    ReadScalar()
        Reads Scalar from a line with format: 
            Name'\t'Value
        After Reading the Value Fill the associated Histogram
*/
void CAMAC_Data::ReadScalar(string line, int varInd)
{
    TString str_scalar;
    double num_scalar;
    char tempChar;
    
    // Debugging - Make sure varInd and variable is consistent
//     cout<<varInd<<" "<<line<<endl;

    /*
        Go Over Each Character and save words in TString Format
            Each Line has a Variable Name and its Value
                Name and Value separated by Tab character '\t'
            Read Until Tab Character and clear the string
                The String after the Tab Character is the Value
            Once read the complete Value
                Covert it to Double
                Save into the Histogram
    */
    for(unsigned int i = 0; i < line.size(); i++){
        tempChar = line[i];
        str_scalar.Append(tempChar);
        
        // The word until Tab Character is Variable Name
        if(tempChar == '\t') str_scalar.Clear();
    }
    
    // Convert string to double
    num_scalar = str_scalar.Atof();
   
    // Fill Histograms
    hists[varInd]->Fill(num_scalar);

}


void CAMAC_Data::WriteRootFile()
{
    cout<<"Updating "<<rootDir<<endl;

    // Overwrite previous Histograms
    for(unsigned int i = 0; i < hists.size(); i++){
        hists[i]->Write("",TObject::kOverwrite);   
    }
    
    // Close File After Writing Histograms
    f_Root->Close();
}

void CAMAC_Data::SetRootFileDir(string input)
{
    rootDir = input;
}

void CAMAC_Data::SetDataFileDir(string input)
{
    dataDir = input;
}

void CAMAC_Data::SetBranchName(string input)
{
    branchName = input;
}


/*
    OpenFileError()
        Called if the Input Data File can not be opened
*/
void CAMAC_Data::OpenDataFileError()
{
    cerr<<"Cannot Open File = "<<dataDir<<endl;
    exit(EXIT_FAILURE);
}

void CAMAC_Data::PrintVariables()
{
    cout<<"Saved Variables:"<<endl;
    for(unsigned int i = 0; i < vars.size(); i++){
        cout<<i<<" "<<vars[i]<<endl;   
    }
}

void CAMAC_Data::OpenFiles()
{   
    // Open Data File
    dataFile.open(dataDir.c_str());
    if (!dataFile.is_open()) OpenDataFileError();

    
    // Open ROOT File
    f_Root = new TFile(rootDir.c_str(),"UPDATE");
    
    // Add Directories to ROOT File (if there is no directory)
    if(f_Root->Get(branchName.c_str()) == NULL) f_Root->mkdir(branchName.c_str());
}

CAMAC_Data::~CAMAC_Data()
{
    dataFile.close();
}


#endif
