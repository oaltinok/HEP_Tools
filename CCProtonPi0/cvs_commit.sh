CCPROTONPI0_V="v2_53"
cvs commit -m "${CCPROTONPI0_V}
v2_53
CCProtonPi0 Updates:
	Structure Changed
		Moved #include statements to header file
			It should be designed like that in the beginning
		Functions for Truth Event moved to another implementation file
			CCProtonPi0_Truth.cpp
		Created Local PDG.h inside Helper
			Framework version was causing compile problems
			More control on the PDG numbers are also useful
	
	Background Type Improved
		Now including Charge Exchanged Pions also

	
-------------------------------------------------------------------------------
NTupleAnalysis Updates:	
	General Improvements
		Removed old study fragments
			variables, histograms, text files, etc…

	BackgroundTool Improved 
		Included Charge Exchanged Pions
		Removed Michel column from output able

	TruthAnalysis Improved
		Performance Improved greatly by disabling the unused branches
		Calculating True Number of Signal Events using this Class
		Added new background branches
		Added Signal Pi0 Histograms for efficiency curves

	CutList Improved
		Removed old Cuts
" .

cvs tag -F ${CCPROTONPI0_V} .
