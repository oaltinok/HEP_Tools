/*
================================================================================
Class: Cut
    Cut Object Manages the Event Selection(Cut) Statistics for 
        CCProtonPi0 Analysis Package
    
    Main Directory:
        Classes/Cut
    
    Usage:
        > #include "Classes/Cut/Cut.h" 
        > Cut cutAll
        > cutAll.set_Name("All");
        > cutAll.inc_nEvent()
        > cutAll.get_nEvent()

    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision: 2015_04_16
================================================================================
*/

#ifndef Cut_h
#define Cut_h

#include <iostream>
#include <string>
#include "CutStat.h"

class Cut{
	public:
		// Constructor
		Cut();
		
        // Set Functions
        void set_Name(std::string inputName);
        
		// Get Functions
		std::string get_Name();
        
        void increment(bool isSignal, bool isStudy1_true = true, bool isStudy2_true = true);
        
        CutStat nEvent;
        CutStat nSignal;
        CutStat nStudy1; // Number of the Events which is studied
        CutStat nStudy2; // Number of the Events which is studied
		  
	private:
		std::string name;

};


#endif
