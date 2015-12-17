#ifndef Graph_cpp
#define Graph_cpp

#include "Graph.h"

using namespace std;

void Graph::AddFormat(string in_title, string in_x_axis, string in_y_axis)
{
    Format tempFormat;
    
    tempFormat.title = in_title;
    tempFormat.x_axis = in_x_axis;
    tempFormat.y_axis = in_y_axis;

    formats.push_back(tempFormat);
}


void Graph::Init()
{
    // Sanity Check
    if ( formats.size() == 0 ){
        cout<<"ERROR: No Graph Format. Did You Read File Correctly?"<<endl;
        exit(EXIT_FAILURE);
    }
    
    vector<double> temp_vector;
  
    // Fill graph_x_axis_vectors and graph_y_axis_vectors with empty vectors
    for (unsigned int i = 0; i < formats.size(); i++){
        graph_x_axis_vectors.push_back(temp_vector);
        graph_y_axis_vectors.push_back(temp_vector);
    }
    
    if (graph_x_axis_vectors.size() == 0 || graph_y_axis_vectors.size() == 0){
        cout<<"ERROR: Graph Vectors Initialization"<<endl;
        exit(EXIT_FAILURE);
    }
    
}

/*
	There are 2 containers for storing vectors
		graph_x_axis_vectors
		graph_y_axis_vectors
	First, clear the each vector in containers
	Second, clear the container
*/
void Graph::Reset()
{
	// Clear Each Vector in the Containers
  	for (unsigned int i = 0; i < formats.size(); i++){
		graph_x_axis_vectors[i].clear();
    	graph_y_axis_vectors[i].clear();
  	}
  	
  	// Clear container for Vectors
  	graph_x_axis_vectors.clear();
  	graph_y_axis_vectors.clear();
  	
  	// Init Vectors using formats vector
  	Init();
}


void Graph::WriteRootFile()
{
    // -------------------------------------------------------------------------
    // Create Graphs
    // -------------------------------------------------------------------------
    TGraph* tempGraph;
    for (unsigned int i = 0; i < graph_x_axis_vectors.size(); i++){
        // Construt TGraph(size, &x_array_first_element, &y_array_first_element)
        tempGraph = new TGraph(graph_x_axis_vectors[i].size(),&graph_x_axis_vectors[i][0],&graph_y_axis_vectors[i][0]);
        tempGraph->SetName(formats[i].title.c_str());
        tempGraph->SetTitle(formats[i].title.c_str());
        tempGraph->GetXaxis()->SetTitle(formats[i].x_axis.c_str());
        tempGraph->GetYaxis()->SetTitle(formats[i].y_axis.c_str());
        graphs.push_back(tempGraph);
    }
    
    // Write Graphs
    for (unsigned int i = 0; i < graphs.size(); i++){
        graphs[i]->Write("",TObject::kOverwrite); 
    }    
}

void Graph::Fill(string var_name, double value)
{
    if ( graph_x_axis_vectors.size() == 0 || graph_y_axis_vectors.size() == 0 ){
        cout<<"ERROR: No Graph Vectors Initialized... Did you initilize correctly!"<<endl;
        exit(EXIT_FAILURE);
    }
    
    // Search Graph Format Vector for the given var_name
    //      Fill Corresponding Histogram if that variable has a histogram
    for ( unsigned int i = 0; i < formats.size(); i++){
        if (var_name.compare(formats[i].x_axis) == 0) {
            graph_x_axis_vectors[i].push_back(value);
        }else if (var_name.compare(formats[i].y_axis) == 0) {
            graph_y_axis_vectors[i].push_back(value);
        }
    }
}


#endif
