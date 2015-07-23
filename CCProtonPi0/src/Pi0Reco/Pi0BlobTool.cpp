#ifndef CCProtonPi0_Pi0BlobTool_cpp
#define CCProtonPi0_Pi0BlobTool_cpp

#include "MinervaDet/DeDetector.h"
#include "MinervaDet/DePlane.h"
#include "Event/IDCluster.h"
#include "Event/IDBlob.h"
#include "Event/Node.h"
#include "Event/Track.h"

#include "Pi0BlobTool.h"

using namespace std;

Pi0BlobTool::Pi0BlobTool()
{
    // Do Nothing!
}

bool Pi0BlobTool::isBlobGood( const Minerva::IDBlob* pi0_blob)
{
    cout<<"Enter Pi0BlobTool::isBlobGood()"<<endl;
    
    // Sanity Check
    if (pi0_blob == NULL) return false;
   
    FillNodesVector_and_SetTrack(pi0_blob);
    
    FillClusterVectors();

    CheckClusterVectors();


    cout<<"Exit Pi0BlobTool::isBlobGood()"<<endl;
    
    return true;
}

void Pi0BlobTool::CheckClusterVectors()
{
    CheckClusterVector(All_clusters);
    CheckClusterVector(X_clusters);
    CheckClusterVector(U_clusters);
    CheckClusterVector(V_clusters);
}

void Pi0BlobTool::CheckClusterVector(std::vector<Minerva::IDCluster*> &clusters)
{
    for (unsigned int i = 0; i < clusters.size(); i++){
        cout<<"View = "<<clusters[i]->view();
        cout<<" Z = "<<clusters[i]->z();
        cout<<" Pos = "<<clusters[i]->position();
        cout<<" L-Pos = "<<clusters[i]->lpos()<<endl; 
    }
    cout<<"----------------------------"<<endl;
}

void Pi0BlobTool::FillClusterVectors()
{
    std::vector<Minerva::IDCluster*> clusters = Pi0Blob_Track->idclusters();
    Minerva::IDCluster* tempCluster = NULL;
    
    for( unsigned int i = 0; i < clusters.size(); i ++){
        tempCluster = new Minerva::IDCluster(*(clusters[i]));
        
        All_clusters.push_back(tempCluster);
        
        // Fill Each View
        if ( tempCluster->view() == Minerva::IDCluster::X) X_clusters.push_back(tempCluster);
        else if ( tempCluster->view() == Minerva::IDCluster::U) U_clusters.push_back(tempCluster);
        else if ( tempCluster->view() == Minerva::IDCluster::V) V_clusters.push_back(tempCluster);
    }
}

void Pi0BlobTool::FillNodesVector_and_SetTrack(const Minerva::IDBlob* pi0_blob)
{
    Minerva::Node* tempNode = NULL;
    SmartRefVector<Minerva::IDCluster> clusters = pi0_blob->clusters();
    SmartRefVector<Minerva::IDCluster>::iterator iter_c;
    
    for (iter_c = clusters.begin(); iter_c != clusters.end(); ++iter_c){
        tempNode = new Minerva::Node(*iter_c);
        nodes.push_back(tempNode);
    }

    // Set Track
    Pi0Blob_Track = new Minerva::Track(nodes);
}


#endif

