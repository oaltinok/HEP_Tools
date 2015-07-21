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

bool Pi0BlobTool::isBlobGood( const Minerva::IDBlob* pi0_blob, const Minerva::DeDetector* idDet )
{
    cout<<"Enter Pi0BlobTool::isBlobGood()"<<endl;
    
    // Sanity Check
    if (pi0_blob == NULL) return false;
   
    FillNodesVector(pi0_blob, idDet);
   
    // Set Track
    Pi0Blob_Track = new Minerva::Track(nodes);
    
    cout<<"Track Parameters"<<endl;
    cout<<"Position = "<<Pi0Blob_Track->position()<<endl;
    cout<<"Slopes = "<<Pi0Blob_Track->slopes()<<endl;
    cout<<"Chi2 = "<<Pi0Blob_Track->chi2()<<endl;
    
    cout<<"Exit Pi0BlobTool::isBlobGood()"<<endl;
    
    return true;
}

void Pi0BlobTool::FillNodesVector(const Minerva::IDBlob* pi0_blob, const Minerva::DeDetector* idDet)
{
    cout<<"Enter Pi0BlobTool::FillNodesVector()"<<endl;
    
    Minerva::Node* tempNode = NULL;
    SmartRefVector<Minerva::IDCluster> clusters = pi0_blob->clusters();
    SmartRefVector<Minerva::IDCluster>::iterator iter_c;
    
    for (iter_c = clusters.begin(); iter_c != clusters.end(); ++iter_c){
        tempNode = new Minerva::Node(*iter_c);
        nodes.push_back(tempNode);
    }
 
    // Check nodes
    for( unsigned int i = 0; i < nodes.size(); i++){
        cout<<i<<" Node View: "<<nodes[i]->view()<<" Position: "<<nodes[i]->position()<<endl;
        AlignNode(nodes[i], idDet);
        cout<<i<<" Node View: "<<nodes[i]->view()<<" Position: "<<nodes[i]->position()<<endl;
        cout<<"----"<<endl;
    }
  

    cout<<"Exit Pi0BlobTool::FillNodesVector()"<<endl;
}

void Pi0BlobTool::AlignNode(Minerva::Node* node, const Minerva::DeDetector* idDet)
{
    Minerva::DePlane const * plane = idDet->getDePlane( node->planeid() );
    node->idcluster()->setLpos( plane->getLPos( node->position() ) );
}

#endif

