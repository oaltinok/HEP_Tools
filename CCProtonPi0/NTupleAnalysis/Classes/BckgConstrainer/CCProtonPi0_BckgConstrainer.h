#ifndef CCProtonPi0_BckgConstrainer_h
#define CCProtonPi0_BckgConstrainer_h

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>

#include "../../Libraries/Folder_List.h"

using namespace std;

struct UniverseWeight
{
    double wgt_SinglePiPlus;
    double wgt_QELike;
    double wgt_WithPi0;
    double err_SinglePiPlus;
    double err_QELike;
    double err_WithPi0;
};

class CCProtonPi0_BckgConstrainer 
{
    public:
        CCProtonPi0_BckgConstrainer(std::string in_file);
        ~CCProtonPi0_BckgConstrainer();

        double GetBckgConstraint(std::string error_name, int hist_ind, std::string bckg_name);
        double GetBckgConstraintErr(std::string error_name, int hist_ind, std::string bckg_name);

    private:
        void ReadBckgConstraints();
        void PrintMap();
        map< string, vector<UniverseWeight> > BckgWeights;

        std::string input_file;
};

#endif

