CCPROTONPI0_V="v2_64"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	Removed Neutrino Energy requirement from Signal Definition
		Applying in NTupleAnalysis Stage

	Background Classification Improved
		Added a more Compact Classification
			BckgWithPi0
			QELike
			SinglePiPlus
			Other

-------------------------------------------------------------------------------
NTupleAnalysis Updates:	
	GetSignalDef()
		Applies Neutrino Energy requirement

	BackgroundTool Improved
		Different Tables for different topologies
		Added latest background types

	Truth Match Plots for Invariant Mass
" .

cvs tag -F ${CCPROTONPI0_V} .

