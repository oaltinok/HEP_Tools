CCPROTONPI0_V="v2_88"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	No Major Changes
---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:
	Added 2 Lateral Error Bands
		EM Energy Scale and Muon Momentum
		Not filling them in this version. 
		Implementation for FillLatErrorBands() will be done in next version
		
	New Class: RandNumGenerator
		Used for varying universes in Lateral Error Bands

	Removed Cross Section for muon_cos_theta
		It was a test case for muon angle problem
" .

cvs tag -F ${CCPROTONPI0_V} .

