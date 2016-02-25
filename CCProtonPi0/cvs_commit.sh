CCPROTONPI0_V="v2_62"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	Pi0 Energy Calculation Changed
		Calorimetric Energy		
			pi0_E_Cal =  E_g1 + E_g2
		Relativistic Energy from Momentum
			pi0_E = sqrt(P2 + m2)
	
	pi0_KE added to truth

-------------------------------------------------------------------------------
NTupleAnalysis Updates:	
	New Cross Section Calculations	
		Pi0 Kinetic Energy
		Pi0 Theta
		QSq

	Analyzer & Plotter Improved for new cross section variables
" .

cvs tag -F ${CCPROTONPI0_V} .

