#ifndef CCProtonPi0_BckgConstrainer_cpp
#define CCProtonPi0_BckgConstrainer_cpp

#include "CCProtonPi0_BckgConstrainer.h"

using namespace std;

CCProtonPi0_BckgConstrainer::CCProtonPi0_BckgConstrainer()
{
    ReadBckgConstraints();
    //PrintMap();
}

CCProtonPi0_BckgConstrainer::~CCProtonPi0_BckgConstrainer()
{
    // Do Nothing
}

double CCProtonPi0_BckgConstrainer::GetBckgConstraint(std::string error_name, int hist_ind, std::string bckg_name)
{
    map<string, vector<UniverseWeight> >::iterator it = BckgWeights.find(error_name);
    
    if (it == BckgWeights.end()){
        std::cout<<"Cannot find "<<error_name<<std::endl;
    }else{
        UniverseWeight unv_wgt = (it->second).at(hist_ind);
   
        if (bckg_name.compare("SinglePiPlus") == 0) return unv_wgt.wgt_SinglePiPlus;
        else if (bckg_name.compare("QELike") == 0) return unv_wgt.wgt_QELike;
        else if (bckg_name.compare("WithPi0") == 0) return unv_wgt.wgt_WithPi0;
        else std::cout<<"Cannot find "<<bckg_name<<std::endl;
    }

    return -1;
}

double CCProtonPi0_BckgConstrainer::GetBckgConstraintErr(std::string error_name, int hist_ind, std::string bckg_name)
{
    map<string, vector<UniverseWeight> >::iterator it = BckgWeights.find(error_name);
    
    if (it == BckgWeights.end()){
        std::cout<<"Cannot find "<<error_name<<std::endl;
    }else{
        UniverseWeight unv_wgt = (it->second).at(hist_ind);
   
        if (bckg_name.compare("SinglePiPlus") == 0) return unv_wgt.err_SinglePiPlus;
        else if (bckg_name.compare("QELike") == 0) return unv_wgt.err_QELike;
        else if (bckg_name.compare("WithPi0") == 0) return unv_wgt.err_WithPi0;
        else std::cout<<"Cannot find "<<bckg_name<<std::endl;
    }

    return -1;
}


void CCProtonPi0_BckgConstrainer::ReadBckgConstraints()
{
    ifstream file;
    file.open(Folder_List::BckgConstraints.c_str());

    if (!file.is_open()){
        std::cout<<"WARNING! Cannot open input file "<<Folder_List::BckgConstraints<<std::endl;
        exit(1);
    }

    // Read Header and Discard (Don't use it)
    std::string line;
    getline(file, line);

    std::string error_name;
    int hist_ind;
    UniverseWeight unv_wgt;
    double dummy;

    // Read Table
    while(!file.eof()){
        getline(file, line);
        std::stringstream line_ss(line);
        line_ss>>error_name>>hist_ind>>dummy>>dummy>>unv_wgt.wgt_SinglePiPlus>>unv_wgt.wgt_QELike>>unv_wgt.wgt_WithPi0>>unv_wgt.err_SinglePiPlus>>unv_wgt.err_QELike>>unv_wgt.err_WithPi0;

        BckgWeights[error_name].push_back(unv_wgt);
    } 

    file.close();
}

void CCProtonPi0_BckgConstrainer::PrintMap()
{
    map<string, vector<UniverseWeight> >::iterator iter;

    for (iter = BckgWeights.begin(); iter != BckgWeights.end(); ++iter ){
        string unv_name = iter->first;

        // Print Values
        vector<UniverseWeight> unv_wgt = iter->second;
        for (unsigned int i = 0; i < unv_wgt.size(); ++i){
            cout<<unv_name<<" "<<i<<" ";
            cout<<unv_wgt[i].wgt_SinglePiPlus<<" ";
            cout<<unv_wgt[i].wgt_QELike<<" ";
            cout<<unv_wgt[i].wgt_WithPi0<<" ";
            cout<<unv_wgt[i].err_SinglePiPlus<<" ";
            cout<<unv_wgt[i].err_QELike<<" ";
            cout<<unv_wgt[i].err_WithPi0<<" ";
            cout<<endl;
        }
    }
}


#endif

