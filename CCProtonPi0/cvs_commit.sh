CCPROTONPI0_V="v2_54"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	Corrected Neutrino Energy Calculation
		Extra Energy is collected in 3 different ways
			Dispersed Blob after Pi0 reconstruction
			Muon Extra Energy
			Rejected Cluster Energy
		Neutrino Energy includes Total Extra Energy 

	Corrected WSq Calculation
		I was using Q2 instead of q2 (q2 = -Q2)
	
-------------------------------------------------------------------------------
NTupleAnalysis Updates:	
	New script FindAllAnaFiles.sh
		Outputs a file listing ALL Ana_Tuple Files
" .

cvs tag -F ${CCPROTONPI0_V} .
