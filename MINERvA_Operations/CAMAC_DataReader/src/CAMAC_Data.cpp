/*
    See CAMAC_Data.h header or Class Information
*/
#ifndef CAMAC_Data_cpp
#define CAMAC_Data_cpp

#include "../CAMAC_DataReader/CAMAC_Data.h"

using namespace std;

CAMAC_Data::CAMAC_Data()
{
    isDebugging = false;


}


void CAMAC_Data::ProcessData()
{
    OpenFiles();
    ReadDataFile();
    WriteRootFile();
}



void CAMAC_Data::ReadInfoLine(string line, vector<string> &v)
{
    if(isDebugging){
        cout<<"----"<<endl;
        cout<<"Reading Info Line:"<<endl; 
        cout<<line<<"\n\n";
    }
    
    stringstream line_stream(line);
    string word;
    
    while(line_stream >> word){
        if(isDebugging) cout<<word<<" | ";
        v.push_back(word);
    }
    if(isDebugging) cout<<endl;
}



/*
    readFile()

*/
void CAMAC_Data::ReadDataFile()
{
    bool isLine_Info = true;
    bool isLine_Version = true;
    bool isLine_VariableNames = true;
    bool isLine_ConvFactors = true;
    bool isLine_Units = true;
    
    bool isInit_Hists = false;
    bool isInit_Graphs = false;
    bool isInit_Freqs = false;
    string line;
    while (!dataFile.eof() )
    {
        getline (dataFile,line);
        
        // Skip Empty Lines
        if (line.size() == 0) continue;
        
        // Read Information Section of the File
        if (isLine_Info){
            
            if (isLine_Version){
                ReadInfoLine(line,version);
                isLine_Version = false;
                continue;
            }
            else if (isLine_VariableNames){
                ReadInfoLine(line,vars);
                isLine_VariableNames = false;
                continue;
            }
            else if (isLine_ConvFactors){
                ReadInfoLine(line,convFactors_str);
                ProcessConvFactors();
                isLine_ConvFactors = false;
                continue;
            }
            else if (isLine_Units){
                ReadInfoLine(line,units);
                isLine_Units = false;
                continue;
            }
            
            // Read Plot Formats
            if(line[0] == 'h'){ 
                ReadHistogramFormat(line);
                continue;
            }
            else if( line[0] == 'g'){ 
                ReadGraphFormat(line);
                continue;
            }
            else if( line[0] == 'f'){
                ReadFrequencyFormat(line);
                continue;
            }
            else{
                if (isDebugging) cout<<"Finished Reading Info Section!"<<endl;
                isLine_Info = false; // Finished Reading Info Section of Data File
            }
        }
        
        // Initialize Histograms If they did NOT INITIALIZED before
        if ( !isInit_Hists ){
            initHistograms();   
            isInit_Hists = true;
        }
        
        if( !isInit_Graphs ){
            initGraphVectors(); 
            isInit_Graphs = true;
        }
        
        if( !isInit_Freqs ){
            isInit_Freqs = true;
        }

        ReadDataLine(line);
        
    }
    
    
}

void CAMAC_Data::ReadDataLine(string line)
{
//     if(isDebugging){
//         cout<<"----"<<endl;
//         cout<<"Reading Data:"<<endl; 
//         cout<<line<<"\n\n";
//     }
    
    stringstream line_stream(line);
    string var_name;
    double value;
    
    line_stream >> var_name >> value;
    
    
    ProcessDataLine(var_name,value);
    
}

void CAMAC_Data::ProcessDataLine(string var_name, double value)
{
    // Get Variable Indice
    int var_ind = GetVariableInd(var_name);
    
    // Get Variable Conversion Factor
    double var_convFactor =  convFactors[var_ind]; 
    
    // Convert Scalar Value to REAL Value using Conversion Factor
    value = value * var_convFactor;
    
    // Fill Histograms - Search Histograms for the Variable and fill if necessary
    FillHistograms(var_name,value);
    
    // Fill Graphs - Search Graphs for the Variable and fill if necessary
    FillGraphVectors(var_name,value);
}

void CAMAC_Data::FillHistograms(string var_name, double value)
{
    if ( f_hists.size() == 0 ){
        cout<<"ERROR: No Histograms Initialized, Did you initilize correctly!"<<endl;
        exit(EXIT_FAILURE);
    }
    
    // Search Histogram Format Vector for the given var_name
    //      Fill Corresponding Histogram if that variable has a histogram
    for ( unsigned int i = 0; i < f_hists.size(); i++){

        if (var_name.compare(f_hists[i].var_name) == 0) {
            hists[i]->Fill(value);
            break;
        }
    }
}


void CAMAC_Data::FillGraphVectors(string var_name, double value)
{
    if ( x_axis_vectors.size() == 0 || y_axis_vectors.size() == 0 ){
        cout<<"ERROR: No Graph Vectors Initialized, Did you initilize correctly!"<<endl;
        exit(EXIT_FAILURE);
    }

    // Search Graph Format Vector for the given var_name
    //      Fill Corresponding Histogram if that variable has a histogram
    for ( unsigned int i = 0; i < f_graphs.size(); i++){
        if (var_name.compare(f_graphs[i].x_axis) == 0) {
            x_axis_vectors[i].push_back(value);
        }else if (var_name.compare(f_graphs[i].y_axis) == 0) {
            y_axis_vectors[i].push_back(value);
        }
    }
}

void CAMAC_Data::ProcessConvFactors()
{
    if ( convFactors_str.size() == 0 ){
        cout<<"ERROR: No Conversion Factors to Process! Did you read data file correctly!"<<endl;
        exit(EXIT_FAILURE);
    }
    
    for (unsigned int i = 0; i < convFactors_str.size(); i++){
        TString s(convFactors_str[i]);
        convFactors.push_back(s.Atof());
    }
    
    if ( convFactors_str.size() != convFactors.size() ){
        cout<<"ERROR: Conversion Factor Process!"<<endl;
        exit(EXIT_FAILURE);
    }
}


int CAMAC_Data::GetVariableInd(string var_name)
{
    if ( vars.size() == 0 ){
        cout<<"ERROR: No Variables, Did you read data file correctly!"<<endl;
        exit(EXIT_FAILURE);
    }
    
    int ind = -1;
    
    for(unsigned int i = 0; i < vars.size(); i++){
        if ( var_name.compare(vars[i]) == 0){ 
            ind = i;
            break;
        }
    }
    
    if (ind != -1){ 
        return ind;
    }else{
        cout<<"ERROR: No Variable = "<<var_name<<endl;
        exit(EXIT_FAILURE);
    }
}



/*
    initHistograms()
        Initializes Histograms inside hists vector
        Uses f_hists vector to get Histogram Information
*/
void CAMAC_Data::initHistograms()
{
    // Sanity Check
    if ( f_hists.size() == 0 ){
        cout<<"ERROR: No Histogram Format. Did You Read File Correctly?"<<endl;
        exit(EXIT_FAILURE);
    }
    
    TH1D* tempHist;
    string title;
    string var_name;
    int nbins;
    double low;
    double high;
    

    // Change Directory to Corresponding Branch
    f_Root->cd(branchName.c_str());

    // -------------------------------------------------------------------------
    // Get Existing Histogram or Create a NEW one
    // -------------------------------------------------------------------------
    for (unsigned int i = 0; i < f_hists.size(); i++){
        
        // Get Histogram Format Info
        title = f_hists[i].title;
        var_name = f_hists[i].var_name;
        nbins = f_hists[i].nbins;
        low = f_hists[i].low;
        high = f_hists[i].low;
        
        // Get Existing Histogram
        tempHist = (TH1D*)gDirectory->Get(var_name.c_str());
        
        // If NO Existing Histogram - Create a New One
        if(tempHist == NULL){

            tempHist = new TH1D( var_name.c_str(),title.c_str(),nbins, low, high );
            tempHist->SetBins(nbins, low, high);
            tempHist->GetXaxis()->SetTitle(var_name.c_str());
            tempHist->GetYaxis()->SetTitle("N(Events)");
        }
        
        // Push Histogram to hists Vector
        hists.push_back(tempHist);
    }
    
    if (hists.size() == 0 ){
        cout<<"ERROR: Histogram Initiliazation"<<endl;
        exit(EXIT_FAILURE);
    }
}

/*
    initGraphVectors()
        Initializes Histograms inside hists vector
        Uses f_hists vector to get Histogram Information
*/
void CAMAC_Data::initGraphVectors()
{
    // Sanity Check
    if ( f_graphs.size() == 0 ){
        cout<<"ERROR: No Graph Format. Did You Read File Correctly?"<<endl;
        exit(EXIT_FAILURE);
    }
    
    vector<double> temp_vector;
    
    // Fill x_axis_vectors and y_axis_vectors with empty TVectorDs
    for (unsigned int i = 0; i < f_graphs.size(); i++){
        x_axis_vectors.push_back(temp_vector);
        y_axis_vectors.push_back(temp_vector);
    }
    
    if (x_axis_vectors.size() == 0 || y_axis_vectors.size() == 0){
        cout<<"ERROR: Graph Vectors Initiliazation"<<endl;
        exit(EXIT_FAILURE);
    }
}



void CAMAC_Data::ReadHistogramFormat(string line)
{
    if(isDebugging){
        cout<<"----"<<endl;
        cout<<"Reading Histogram Format Line:"<<endl; 
        cout<<line<<"\n\n";
    }
    stringstream line_stream(line);
    string word;
    string key; // First word is the Key which is not saved
        
    format_hist tempHist;
    
    line_stream >> key 
                >> tempHist.title >> tempHist.var_name
                >> tempHist.nbins >> tempHist.low >> tempHist.high;
      
    if(isDebugging){
        cout<<"Title = "<<tempHist.title<<endl;
        cout<<"Variable = "<<tempHist.var_name<<endl;
        cout<<"nbins = "<<tempHist.nbins<<endl;
        cout<<"low = "<<tempHist.low<<endl;
        cout<<"high = "<<tempHist.high<<endl;
    }
    
    f_hists.push_back(tempHist);
   
}

void CAMAC_Data::ReadGraphFormat(string line)
{
    if(isDebugging){
        cout<<"----"<<endl;
        cout<<"Reading Graph Format Line:"<<endl; 
        cout<<line<<"\n\n";
    }
    stringstream line_stream(line);
    string word;
    string key; // First word is the Key which is not saved
        
    format_graph tempGraph;
    
    line_stream >> key 
                >> tempGraph.title >> tempGraph.x_axis >> tempGraph.y_axis;
      
      
    if(isDebugging){
        cout<<"Title = "<<tempGraph.title<<endl;
        cout<<"x_axis = "<<tempGraph.x_axis<<endl;
        cout<<"y_axis = "<<tempGraph.y_axis<<endl;
    }
    
    f_graphs.push_back(tempGraph);
   
}

void CAMAC_Data::ReadFrequencyFormat(string line)
{
    if(isDebugging){
        cout<<"----"<<endl;
        cout<<"Reading Frequency Format Line:"<<endl; 
        cout<<line<<"\n\n";
    }
    stringstream line_stream(line);
    string word;
    string key; // First word is the Key which is not saved
    int value;
        
    format_frequency tempFreq;
    
    line_stream >> key 
                >> tempFreq.title;
                
    while(line_stream >> value){
        tempFreq.vars.push_back(value);
    }
    
    if(isDebugging){
        cout<<"Title = "<<tempFreq.title<<endl;
        for(unsigned int i = 0; i < tempFreq.vars.size(); i++){
            cout<<tempFreq.vars[i]<<endl;    
        }
    }
    
    f_frequencies.push_back(tempFreq);
   
}




void CAMAC_Data::WriteRootFile()
{
    cout<<"Updating "<<rootDir<<endl;

    // Write Histograms
    for(unsigned int i = 0; i < hists.size(); i++){
        hists[i]->Write("",TObject::kOverwrite);   
    }
    
    // Create Graphs
    TGraph* tempGraph;

    
    for (unsigned int i = 0; i < x_axis_vectors.size(); i++){
        // Construt TGraph(size, &x_array_first_element, &y_array_first_element)
        tempGraph = new TGraph(x_axis_vectors[i].size(),&x_axis_vectors[i][0],&y_axis_vectors[i][0]);
        tempGraph->SetName(f_graphs[i].title.c_str());
        tempGraph->SetTitle(f_graphs[i].title.c_str());
        tempGraph->GetXaxis()->SetTitle(f_graphs[i].x_axis.c_str());
        tempGraph->GetYaxis()->SetTitle(f_graphs[i].y_axis.c_str());
        graphs.push_back(tempGraph);
    }
    
    // Write Graphs
    for (unsigned int i = 0; i < graphs.size(); i++){
//         graphs[i]->Draw();
//         graphs[i]->Write();
        graphs[i]->Write("",TObject::kOverwrite); 
    }
    
    // Close File After Writing 
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

void CAMAC_Data::PrintVector(vector<string> &v)
{
    for(unsigned int i = 0; i < v.size(); i++){
        cout<<i<<" "<<v[i]<<endl;   
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
    if(isDebugging) cout<<"CAMAC_Data Destructor"<<endl;
    dataFile.close();
}


#endif
