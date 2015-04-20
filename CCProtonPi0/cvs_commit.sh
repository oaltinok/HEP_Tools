cvs commit -m "v2_06

General Updates:
    No Updates to CCProtonPi0 Package 
--------------------------------------------------------------------------------
NTupleAnalysis Updates:
    Revised Constructors for all Classes
        Old Style = Empty Constructor + initialize()
        New Style = Overloaded Constructor
        Analyzer Uses Member Initialization List to Init other Objects
        Collected other initialization functions under initFunctions.cpp
        
    Added Class PIDTool to the CVS (forgot to add on previous release)
    
    Created a New Class: CutList
        Moved all Cut Number related functions under this Class
        Removed Channel Tag from the Package
            Only Cut Tables had that tag and no needed for future
    New Class: CutList
        Member Variables are the Cut Numbers which represent each Selection 
            in the Analysis
        Creates the Cut Table for all topologies in the Analysis" .