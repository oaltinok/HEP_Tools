CCPROTONPI0_V="v2_56"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	No Major Changes	
-------------------------------------------------------------------------------
NTupleAnalysis Updates:	
	Data Correction Methods for Cross Section Calculation Implemented
		Background Subtraction
		Unfolding
		Efficiency Correction
	
	Tested Muon and Pi0 Momentum with Data Correction Methods

	New Plotting Functions for Cross Section Variables
" .

cvs tag -F ${CCPROTONPI0_V} .

