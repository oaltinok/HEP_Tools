CCPROTONPI0_V="v2_65"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	No Major Changes
-------------------------------------------------------------------------------
NTupleAnalysis Updates:	
	Three Side Bands:
		Events with Michel Electron
		2 Track Events with Low Proton Score
		Events in Low Inv Mass Region before Shower Quality cut

	New Class: SideBandTool
		Uses TFractionFitter to fit Background on Side Bands to Data
		Initial Commit Features:
			Calculates fraction of each background in Data
			Writes out a table as a comparison to MC Ratios vs Fit Ratios
			
	MINOS Correction for Every Playlist Implemented
" .

cvs tag -F ${CCPROTONPI0_V} .

