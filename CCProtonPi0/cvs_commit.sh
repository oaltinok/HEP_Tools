CCPROTONPI0_V="v2_63"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	Signal Definition Changed
		Added min neutrino Energy = 2 GeV
		Acceptance: No events having energy less than 2GeV

	Fixed truth_pi0_KE variable
		was writing pi0_P on it

	Fixed Angles
		was using angles wrt Beam
		_theta and _phi are in Lab Frame
		_theta_beam and _phi_beam are in Beam Frame

	Collected Shower Recovery variables under new Study
		Use StudyShowerRecovery in options file
	
	Collected Unused Energy Truth Matching under new Study
		Use StudyUnusedEnergy in options file

-------------------------------------------------------------------------------
NTupleAnalysis Updates:	

" .

cvs tag -F ${CCPROTONPI0_V} .

