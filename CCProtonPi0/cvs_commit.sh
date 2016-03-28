CCPROTONPI0_V="v2_74"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	No Major Changes	
---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:	
	Global EM Energy Correction
		Correcting ALL MC Reconstructed Variables calculated by EM Energy Calibration
			Correction Factor = ~1.08

	Improved Unfolding Study
		Added Muon cos(theta)
		Plotting Statistical Errors after each iteration
		Plotting function arguments converted to const
" .

cvs tag -F ${CCPROTONPI0_V} .

