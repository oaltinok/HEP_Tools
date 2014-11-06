/*
================================================================================
Class: Cut
    Cut Object Manages the Event Selection(Cut) Statistics for 
        CCProtonPi0 Analysis Package
    
    Main Directory:
        Classes/Cut
    
    Usage:
        > #include "Classes/Cut/Cut.h" 
        > Cut cutAll("All")
        > cutAll.inc_nEvent()
        > cutAll.get_nEvent()

    Author:         Ozgur Altinok  - ozgur.altinok@tufts.edu
    Last Revision: 2014_11_05
================================================================================
*/

#ifndef Cut_h
#define Cut_h

#include <iostream>
#include <string>

class Cut{
	public:
		// Constructor
		Cut();
		
        // Set Functions
        void set_Name(std::string inputName);
        
		// Get Functions
		std::string get_Name();
		double get_nEvent();
		double get_nSignal();
        double get_nSignal_Gold();
        double get_nSignal_Silver1();
		  
		// Increment Functions
		void inc_nEvent();
		void inc_nSignal();
        void inc_nSignal_Gold();
        void inc_nSignal_Silver1();
        
		  
		  
	private:
		std::string name;
		double nEvent;
		double nSignal;
        double nSignal_Gold;
        double nSignal_Silver1;
};


#endif
