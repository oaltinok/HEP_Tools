/*
    See CAMAC_Data.h header or Class Information
*/
#ifndef CAMAC_Data_cpp
#define CAMAC_Data_cpp

#include "../CAMAC_DataReader/CAMAC_Data.h"

using namespace std;

/*
    Default Constructor
*/
CAMAC_Data::CAMAC_Data()
{
    cout<<"\n\n";
    cout<<"#####################################################################"<<endl;
    cout<<"#                                                                   #"<<endl;
    cout<<"#                    Welcome CAMAC_DataReader                       #"<<endl;
    cout<<"#                                                                   #"<<endl;
    cout<<"#####################################################################"<<endl;
    cout<<"\n\n";
    
    isDebugging = false;
    isNewSubrun = true;
    SetBranchName("CAMAC_Data");
    SetRootFileDir("/minerva/data/testbeam2/nearonline/NearlineCurrentHistos.root");
    
    cout<<"\n\n";
}

/*
    Main Function that controls the Data Processing
        Open Files
        Read Data File
        Process Read Data
        Write Results
*/
void CAMAC_Data::ProcessData()
{
    OpenFiles();
    ReadDataFile();
    WriteRootFile();
    ClearVectors();
}


/*
    Macro Function for reading and processing single information line 
*/
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
    Main Function that controls the reading the Data File
        Read Information Lines
            Initialize Data Structures according to Information Given
        Read Data Lines
            Process Data 
            Fill Plots with the Data
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
            initFreqVectors();
            isInit_Freqs = true;
        }

        // Read Data Lines
        ReadDataLine(line);
    }
}


/*
    Reads a single Data line and processes data  
*/
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


/*
    Process given data
        Inputs:
            variable name
            value
*/
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
    
    // Fill Frequency Plots - Check Frequency Arrays to count the Variable
    FillFrequencyPlots(var_ind);
}


void CAMAC_Data::FillHistograms(string var_name, double value)
{
    if ( f_hists.size() == 0 ){
        cout<<"ERROR: No Histograms Initialized... Did you initilize correctly!"<<endl;
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
    if ( graph_x_axis_vectors.size() == 0 || graph_y_axis_vectors.size() == 0 ){
        cout<<"ERROR: No Graph Vectors Initialized... Did you initilize correctly!"<<endl;
        exit(EXIT_FAILURE);
    }

    // Search Graph Format Vector for the given var_name
    //      Fill Corresponding Histogram if that variable has a histogram
    for ( unsigned int i = 0; i < f_graphs.size(); i++){
        if (var_name.compare(f_graphs[i].x_axis) == 0) {
            graph_x_axis_vectors[i].push_back(value);
        }else if (var_name.compare(f_graphs[i].y_axis) == 0) {
            graph_y_axis_vectors[i].push_back(value);
        }
    }
}

void CAMAC_Data::FillFrequencyPlots(int ind)
{
    if ( f_frequencies.size() == 0 ){
        cout<<"ERROR: No Frequency Plots Initialized... Did you initilize correctly!"<<endl;
        exit(EXIT_FAILURE);
    }
    
    for (unsigned int i = 0; i < f_frequencies.size(); i++){
        if ( f_frequencies[i].vars.size() == 0 ){
            cout<<"ERROR: Frequency Vector = "<<i<<" is NOT Initialized.. Did you initilize correctly!"<<endl;
            exit(EXIT_FAILURE);
        }
    }
    
       
    // Check Frequency Format's vars array - Which Holds the information
    for(unsigned int i = 0; i< f_frequencies.size(); i++){
        if (f_frequencies[i].vars[ind] == 1 ){
            freq_y_axis_vectors[i][ind] = freq_y_axis_vectors[i][ind]+1;   
        }
    }

}


/*
    Conversion Factors read as strings
        This function converts strings to double
*/
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
        
        // If it is a new subrun create empty histograms
        if( isNewSubrun ){
            tempHist = new TH1D( var_name.c_str(),title.c_str(),nbins, low, high );
            tempHist->SetBins(nbins, low, high);
            tempHist->GetXaxis()->SetTitle(var_name.c_str());
            tempHist->GetYaxis()->SetTitle("N(Events)");
        }else{
            // Get Existing Histogram
            tempHist = (TH1D*)gDirectory->Get(var_name.c_str());
            
            // If NO Existing Histogram
            //    Create a New One
            if(tempHist == NULL){
                tempHist = new TH1D( var_name.c_str(),title.c_str(),nbins, low, high );
                tempHist->SetBins(nbins, low, high);
                tempHist->GetXaxis()->SetTitle(var_name.c_str());
                tempHist->GetYaxis()->SetTitle("N(Events)");
            }
        }
        
        // Push Histogram to hists Vector
        hists.push_back(tempHist);
    }
    
    if (hists.size() == 0 ){
        cout<<"ERROR: Histogram Initialization"<<endl;
        exit(EXIT_FAILURE);
    }
}

/*
    initGraphVectors()
        Initializes Graph Vectors used to generate TGraphs before writing file
        Uses f_graphs vector to get Graph Information
*/
void CAMAC_Data::initGraphVectors()
{
    // Sanity Check
    if ( f_graphs.size() == 0 ){
        cout<<"ERROR: No Graph Format. Did You Read File Correctly?"<<endl;
        exit(EXIT_FAILURE);
    }
    
    vector<double> temp_vector;
    
    // Fill graph_x_axis_vectors and graph_y_axis_vectors with empty vectors
    for (unsigned int i = 0; i < f_graphs.size(); i++){
        graph_x_axis_vectors.push_back(temp_vector);
        graph_y_axis_vectors.push_back(temp_vector);
    }
    
    if (graph_x_axis_vectors.size() == 0 || graph_y_axis_vectors.size() == 0){
        cout<<"ERROR: Graph Vectors Initialization"<<endl;
        exit(EXIT_FAILURE);
    }
}

/*
    initFreqVectors()
        Initializes Frequency Vectors used to generate TGraphs before writing file
        Uses f_frequencies vector to get Frequency Plots Information
*/
void CAMAC_Data::initFreqVectors()
{
    // Sanity Check
    if ( f_frequencies.size() == 0 ){
        cout<<"ERROR: No Frequency Plot Format. Did You Read File Correctly?"<<endl;
        exit(EXIT_FAILURE);
    }
    
    vector<double> temp_vector;
    // Fill freq_x_axis_vectors and freq_y_axis_vectors with empty vectors
    for (unsigned int i = 0; i < f_frequencies.size(); i++){
        freq_x_axis_vectors.push_back(temp_vector);
        freq_y_axis_vectors.push_back(temp_vector);
    }
    
    // Initialize Axis Vectors
    //      x_axis = Variable Indice
    //      y_axis = All Zeros (will be incremented inside ProcessDataLine)
    for(unsigned int i = 0; i < freq_x_axis_vectors.size(); i++){
        for(unsigned int j = 0; j < f_frequencies[i].vars.size(); j++){
            freq_x_axis_vectors[i].push_back(j);
            freq_y_axis_vectors[i].push_back(0.0);
        }
    }
    
    // Debugging
    if(isDebugging){
        for(unsigned int i = 0; i < freq_x_axis_vectors.size(); i++){
            cout<<"freq_x_axis_vectors["<<i<<"] = ";
            for(unsigned int j = 0; j < freq_x_axis_vectors[i].size(); j++){
                cout<<freq_x_axis_vectors[i][j]<<" ";
            }
            cout<<endl;
        }
        for(unsigned int i = 0; i < freq_y_axis_vectors.size(); i++){
            cout<<"freq_y_axis_vectors["<<i<<"] = ";
            for(unsigned int j = 0; j < freq_y_axis_vectors[i].size(); j++){
                cout<<freq_y_axis_vectors[i][j]<<" ";
            }
            cout<<endl;
        } 
    }

    if (freq_x_axis_vectors.size() == 0 || freq_y_axis_vectors.size() == 0){
        cout<<"ERROR: Frequency Vectors Initialization"<<endl;
        exit(EXIT_FAILURE);
    }
}

/*
    See CAMAC_DataReader/PlotsFormats.h for the structure of format_hist
*/
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

/*
    See CAMAC_DataReader/PlotsFormats.h for the structure of format_graph
*/
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

/*
    See CAMAC_DataReader/PlotsFormats.h for the structure of format_frequency
*/
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

/*
    ROOT File is Opened using "UPDATE" Option
    Following Information will be updated inside f_Root
        Histograms
            Updates the previous histogram - increases N(Events)
        Graphs
            Overs the previous x vs y Graphs
                Saves only the information read from DataFile 
                discards the previous graph
        Frequency Plots
            Overwrites the previous Frequency Plots
                Saves only the information read from DataFile 
                discards the previous graph
*/
void CAMAC_Data::WriteRootFile()
{

    // -------------------------------------------------------------------------
    // Write Histograms
    // -------------------------------------------------------------------------
    for(unsigned int i = 0; i < hists.size(); i++){
        hists[i]->Write("",TObject::kOverwrite);   
    }
    
    // -------------------------------------------------------------------------
    // Create Graphs
    // -------------------------------------------------------------------------
    FillSubrunGraphVectors();
    TGraph* tempGraph;
    for (unsigned int i = 0; i < graph_x_axis_vectors.size(); i++){
        // Construt TGraph(size, &x_array_first_element, &y_array_first_element)
        tempGraph = new TGraph(subrun_graph_x_axis_vectors[i].size(),&subrun_graph_x_axis_vectors[i][0],&subrun_graph_y_axis_vectors[i][0]);
        tempGraph->SetName(f_graphs[i].title.c_str());
        tempGraph->SetTitle(f_graphs[i].title.c_str());
        tempGraph->GetXaxis()->SetTitle(f_graphs[i].x_axis.c_str());
        tempGraph->GetYaxis()->SetTitle(f_graphs[i].y_axis.c_str());
        graphs.push_back(tempGraph);
    }
    
    // Write Graphs
    for (unsigned int i = 0; i < graphs.size(); i++){
        graphs[i]->Write("",TObject::kOverwrite); 
    }
    
    // -------------------------------------------------------------------------
    // Create Frequency Plots
    // -------------------------------------------------------------------------
    FillSubrunFrequencyPlotVectors();
    TGraph* tempFreq;
    
    for (unsigned int i = 0; i < freq_x_axis_vectors.size(); i++){
        // Construt TGraph(size, &x_array_first_element, &y_array_first_element)
        tempFreq = new TGraph(subrun_freq_x_axis_vectors[i].size(),&subrun_freq_x_axis_vectors[i][0],&subrun_freq_y_axis_vectors[i][0]);
        tempFreq->SetName(f_frequencies[i].title.c_str());
        tempFreq->SetTitle(f_frequencies[i].title.c_str());
        tempFreq->GetXaxis()->SetTitle("Variable Codes");
        tempFreq->GetYaxis()->SetTitle("N(Events)");
        freq_plots.push_back(tempFreq);
    }
    
    // Write Frequency Plots
    for (unsigned int i = 0; i < freq_plots.size(); i++){
        freq_plots[i]->Write("",TObject::kOverwrite); 
    }
    
    // Close File After Writing 
    f_Root->Close();
    dataFile.close();
    
    cout<<"Succesfully updated ROOT File = "<<rootDir<<endl;
    cout<<"	...Press CTRL+C to stop CAMAC_DataReader Execution"<<endl;
    cout<<endl;
}

void CAMAC_Data::FillSubrunGraphVectors()
{
    vector< double > temp_vector;
    
    // if subrun_graph_x_axis_vectors is empty - Initialize
    //      Fill subrun_graph_x_axis_vectors and 
    //      Fill subrun_ graph_y_axis_vectors with empty vectors
    if(subrun_graph_x_axis_vectors.size() == 0){
        for (unsigned int i = 0; i < f_graphs.size(); i++){
            subrun_graph_x_axis_vectors.push_back(temp_vector);
            subrun_graph_y_axis_vectors.push_back(temp_vector);
        }
    }
    
    // Push Back all values in  graph_x_axis_vectors to subrun_graph_x_axis_vectors
    for(unsigned int i = 0; i < graph_x_axis_vectors.size(); i++  ){
        for (unsigned int j = 0; j < graph_x_axis_vectors[i].size(); j++ ){
            subrun_graph_x_axis_vectors[i].push_back(graph_x_axis_vectors[i][j]);
            subrun_graph_y_axis_vectors[i].push_back(graph_y_axis_vectors[i][j]);
        }
    } 
}

void CAMAC_Data::FillSubrunFrequencyPlotVectors()
{
    vector<double> temp_vector;
    
    // if subrun_freq_x_axis_vectors is Empty
    // Initialize subrun_freq_x_axis_vectors and copy freq_x_axis_vectors
    if (subrun_freq_x_axis_vectors.size() == 0 ){
        // Fill freq_x_axis_vectors and freq_y_axis_vectors with empty vectors
        for (unsigned int i = 0; i < f_frequencies.size(); i++){
            subrun_freq_x_axis_vectors.push_back(temp_vector);
            subrun_freq_y_axis_vectors.push_back(temp_vector);
        }
        
        // Push Back all values in  freq_x_axis_vectors to freq_x_axis_vectors
        for(unsigned int i = 0; i < freq_x_axis_vectors.size(); i++  ){
            for (unsigned int j = 0; j < freq_x_axis_vectors[i].size(); j++ ){
                subrun_freq_x_axis_vectors[i].push_back(freq_x_axis_vectors[i][j]);
                subrun_freq_y_axis_vectors[i].push_back(freq_y_axis_vectors[i][j]);
            }
        } 
    }
    // Else Increment frequency
    else{
        // subrun_freq_x_axis_vectors[i] - Stays Same
        // subrun_freq_y_axis_vectors[j] - Increased by the amount in freq_y_axis_vectors[j]
        for(unsigned int i = 0; i < freq_y_axis_vectors.size(); i++  ){
            for (unsigned int j = 0; j < freq_y_axis_vectors[i].size(); j++ ){
                subrun_freq_y_axis_vectors[i][j] = subrun_freq_y_axis_vectors[i][j] + freq_y_axis_vectors[i][j];
            }
        } 
        
    }
    
}

void CAMAC_Data::SetRootFileDir(string input)
{
    rootDir = input;
    cout<<"Output ROOT File = "<<rootDir<<endl;
}

void CAMAC_Data::SetDataFileDir(string input)
{
    dataDir = input;
    ProcessFileName(input);
    cout<<"Input Data File = "<<dataDir<<endl;
    // Update Variables
    latest_runNumber = runNumber;
    latest_subrunNumber = subrunNumber;
    latest_timeStamp = timeStamp;
}

bool CAMAC_Data::ChangeDataFileDir(string input)
{
    ProcessFileName(input);
    if( IsNewFile() ){
        
        // If it is a NEW Subrun
        //   Clear Subrun Vectors for Graphs and Frequency Plots
        //   mark isNewSubrun = true - will create new Histograms   
        if(IsNewSubrun()){
            isNewSubrun = true;
            ClearSubrunVectors();
        }else{
            isNewSubrun = false;   
        }
        
        // Update Variables
        latest_runNumber = runNumber;
        latest_subrunNumber = subrunNumber;
        latest_timeStamp = timeStamp;
        
        // Update dataDir
        dataDir = input;
        

        // Return True
        return true;
    }else{
        // cout<<"No New CAMAC Readout File"<<endl;
        return false;
    }
}

bool CAMAC_Data::IsNewFile()
{
    if( latest_runNumber.compare(runNumber) == 0 &&
        latest_subrunNumber.compare(subrunNumber) == 0 &&
        latest_timeStamp.compare(timeStamp) == 0 ){
        
        return false;
    }else{
        return true;   
    }
}

bool CAMAC_Data::IsNewSubrun()
{    
    if( latest_runNumber.compare(runNumber) == 0 && 
        latest_subrunNumber.compare(subrunNumber) == 0){
        
        return false;
    }else{
        return true;   
    }  
}

void CAMAC_Data::ClearSubrunVectors()
{
    // Vectors for Graphs
    subrun_graph_x_axis_vectors.clear();
    subrun_graph_y_axis_vectors.clear();
    
    // Vectors for Frequency Plots
    subrun_freq_x_axis_vectors.clear();
    subrun_freq_y_axis_vectors.clear();
}

void CAMAC_Data::ClearVectors()
{
        // Vectors for Information Section
        version.clear();
        vars.clear();
        convFactors_str.clear();
        convFactors.clear(); // convFactor read as string then converted to double  
        units.clear();
        
        // Vectors for Plot Formats
        f_hists.clear();
        f_graphs.clear();
        f_frequencies.clear();
        
        // Histogram Vector
        hists.clear();
        
        // Vectors for Graphs
        graphs.clear();
        graph_x_axis_vectors.clear();
        graph_y_axis_vectors.clear();
        
        // Vectors for Frequency Plots
        freq_plots.clear();
        freq_x_axis_vectors.clear();
        freq_y_axis_vectors.clear();
    
}

void CAMAC_Data::SetBranchName(string input)
{
    branchName = input;
    cout<<"Folder Name in ROOT File = "<<branchName<<endl;
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

/*
    Opens Following Files
        dataFile = CAMAC Readout File
        rootFile = ROOT File with "UPDATE" option
        
    Branch Name is the Folder Name which is created inside ROOT file
        Specific folder for the CAMAC_DataReader Package
*/
void CAMAC_Data::OpenFiles()
{   
    // Open Data File
    dataFile.open(dataDir.c_str());
    if (!dataFile.is_open()) OpenDataFileError();
    
    // Open ROOT File
    f_Root = new TFile(rootDir.c_str(),"UPDATE");
    
    // Add Directories to ROOT File (if there is no directory)
    if(f_Root->Get(branchName.c_str()) == NULL) f_Root->mkdir(branchName.c_str());
    
    cout<<"Run/Subrun = "<<latest_runNumber<<"/"<<latest_subrunNumber;
    cout<<" Time = "<<latest_timeStamp<<" is processing!"<<endl;
}


void CAMAC_Data::ProcessFileName(string input)
{
    char tempChar;
    string tempString;
    vector<string> v_input;
    
    // Vectorize Each part in the input
    for(unsigned int i = 0; i < input.size(); i++){
        tempChar = dataDir[i];
        
        if( tempChar == '_' ){
           v_input.push_back(tempString);
           tempString.clear();
        }else{
            tempString.push_back(tempChar);    
        }
    }
    
    runNumber = v_input[1];
    subrunNumber = v_input[2];
    timeStamp = v_input[5];
    
    // Debugging
//     for(unsigned int i = 0; i < v_input.size(); i++) {
//         cout<<v_input[i]<<endl;
//     }
}

CAMAC_Data::~CAMAC_Data()
{
    if(isDebugging) cout<<"CAMAC_Data Destructor"<<endl;
    dataFile.close();
}

string CAMAC_Data::GetRunNumber()
{
    return runNumber;
}

string CAMAC_Data::GetSubrunNumber()
{
    return subrunNumber;
}

string CAMAC_Data::GetTimeStamp()
{
    return timeStamp; 
}

string CAMAC_Data::GetFileName()
{
    return dataDir;
}


#endif
