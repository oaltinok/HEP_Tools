#ifndef CCProtonPi0_Pi0BlobTool_h
#define CCProtonPi0_Pi0BlobTool_h

#include <iostream>
#include <vector>

//-- Forward Declarations
#include "Event/MinervaEventFwd.h"

class Pi0BlobTool
{
    public:
        Pi0BlobTool();
        bool isBlobGood( const Minerva::IDBlob* pi0_blob, const Minerva::DeDetector* idDet);

    private:
        std::vector<Minerva::Node*> nodes;
        SmartRef<Minerva::Track> Pi0Blob_Track;
        
        void FillNodesVector(const Minerva::IDBlob* pi0_blob, const Minerva::DeDetector* idDet);
        void AlignNode(Minerva::Node* node, const Minerva::DeDetector* idDet);

};

#endif

