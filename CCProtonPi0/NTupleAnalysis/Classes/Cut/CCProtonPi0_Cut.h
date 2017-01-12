/*
================================================================================
Class: CCProtonPi0_Cut
    CCProtonPi0_Cut Object Manages the Event Selection(CCProtonPi0_Cut) Statistics for 
        CCProtonPi0 Analysis Package
    
    Main Directory:
        Classes/CCProtonPi0_Cut
    
    Usage:
        > #include "Classes/CCut/CCProtonPi0_Cut.h" 
        > Cut CutAll
        > CutAll.set_Name("All");
        > CutAll.inc_nEvent()
        > CutAll.get_nEvent()

    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision: 2015_05_07
================================================================================
*/

#ifndef CCProtonPi0_Cut_h
#define CCProtonPi0_Cut_h

#include <iostream>
#include <string>
#include "../Counter/CCProtonPi0_Counter.h"

class CCProtonPi0_Cut{
	public:
		// Constructor
		CCProtonPi0_Cut();
		
        // Set Functions
        void set_Name(std::string inputName);
        
		// Get Functions
		std::string get_Name();
        
        void increment(bool isSignal, bool isStudy1_true = true, bool isStudy2_true = true, double wgt = 1.0);
        
        CCProtonPi0_Counter nEvent;
        CCProtonPi0_Counter nSignal;
        CCProtonPi0_Counter nStudy1; // Number of the Events which is studied
        CCProtonPi0_Counter nStudy2; // Number of the Events which is studied
		  
	private:
		std::string name;

};


#endif
