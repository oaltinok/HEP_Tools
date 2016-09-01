CCPROTONPI0_V="v2_95"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	Included ImprovedMichelTool 
		Test purposes
---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:
	Test Histograms for ImprovedMichelTool
		Not using the tool itself, just testing

	Corrected Michel Errors
		Applying MichelFake Error only to events without true Michel
		Applying MichelTrue Error only to events with true Michel

	Unfolding Study up to 10 Iterations
		Bin by Bin comparison of Truth-Data, Stat Error, Systematic Error
		Creates tables for a MATLAB Algorithm

	Vertical Error Band for Unfolding Uncertainty Added
		UnfoldingErr = 1%

" .

cvs tag -F ${CCPROTONPI0_V} .

