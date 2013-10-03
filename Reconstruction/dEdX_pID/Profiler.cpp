#include "Profiler.h"

using namespace std;

Profiler::Profiler()
{
    // Do Nothing for now
}

Profiler::Profiler(char input_view)
{
    set_view(input_view);
}

Profiler::~Profiler()
{
    // Do Nothing for now
}

void Profiler::add_dEdX(double input_dEdX)
{
    dEdX.push_back(input_dEdX);
}

void Profiler::add_posZ(double input_Z)
{
    posZ.push_back(input_Z);
}

void Profiler::findOutliers()
{
    double deltadEdX;
    cout<<debugName<<"Finding outliers for "<<view<<"-view"<<endl;
    
    for(unsigned int i = 0; (i < dEdX.size()-1) && (dEdX.size() > 1); i++){
        deltadEdX = dEdX[i] - dEdX[i+1];
        cout<<debugName<<"delta dEdX = "<<deltadEdX<<" posZ = "<<posZ[i]<<endl;
        if(deltadEdX > refDeltadEdX){
            cout<<debugName<<"delta dEdX = "<<deltadEdX<<" posZ = "<<posZ[i]<<"Marked as tobeRemoved!"<<endl;
            markAllBefore(i);
        }
    }
    
}



// Marks all Z-Positions before that indice (if not marked before)
void Profiler::markAllBefore(unsigned int ind)
{
    for(unsigned int i = 0; i <= ind; i++){
        if( !isMarked(posZ[i]) ){
            posZ_tobeRemoved.push_back(posZ[i]);
        }
    }

}



void Profiler::set_isForward(bool input_isForward)
{
    isForward = input_isForward;
}

bool Profiler::get_isForward()
{
    return isForward;
}

bool Profiler::isMarked(double inputZ)
{
    bool isMarked;
    isMarked = false;
    
    for( unsigned int i = 0; i < posZ_tobeRemoved.size(); i++){
        if( fabs(inputZ - posZ_tobeRemoved[i]) < 0.001){
            isMarked = true;
        }
    }
    
    return isMarked;
}

bool Profiler::needRemoval(double inputZ)
{
    // Can be improved for performance: posZ_tobeRemoved is an ordered list.
    // For now it compares with every point
    bool removeNode;
    removeNode = isMarked(inputZ);
    
    if(removeNode){
        cout<<debugName<<inputZ<<" marked to be Removed!"<<endl;
    }
    
    return removeNode;

}

void Profiler::set_view(char input_view)
{
    view = input_view;
}

char Profiler::get_view()
{
    return view;
}






