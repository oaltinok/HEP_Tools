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
    Last Revision: 2014_11_10
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
        
        void increment(bool isSignal, bool isGold, bool isSilver1, bool isDIS);
        
        CutStat nEvent;
        CutStat nSignal;
        CutStat nSignal_Gold;
        CutStat nSignal_Silver1;
        CutStat nDIS;
		  
	private:
		std::string name;

};


#endif
