CCPROTONPI0_V="v2_49"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	No Major Changes

--------------------------------------------------------------------------------
NTupleAnalysis Updates:	
	New Package: TruthAnalysis
		Class: TruthAnalyzer
			Loops over Truth Branch of CCProtonPi0 NTuples
				Get Truth Branch variables using makeClass/Truth
			Outputs Signal & Background information

	Plotter Improved
		MC Only Background Plots added
" .

cvs tag -F ${CCPROTONPI0_V} .
