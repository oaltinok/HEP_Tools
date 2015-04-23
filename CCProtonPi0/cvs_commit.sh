cvs commit -m "v2_08

CCProtonPi0 Updates:
    Michel Tool
        New Parameters for Detected Michel Prongs
    Removed readme file from doc/
--------------------------------------------------------------------------------
NTupleAnalysis Updates:
    Wextra Added to the MakeFile
        Fixed all warnings due to Wextra
    New Class: Interaction
        Removed Interaction specific variables and functions from Analyzer 
            and created Interaction Class
    Major Revision of Analyzer Class
        Combined All implementations files under single file Analyzer.cpp
        Removed Unused functions and variables
    Revised Plotter Class
        Modified Function Calls
        Removed Unused parameters" .