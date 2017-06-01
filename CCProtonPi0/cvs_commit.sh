CCPROTONPI0_V="v3_11"
cvs commit -m "${CCPROTONPI0_V}
v3_11
CCProtonPi0 Updates:
	No major changes

---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:
	Background Classification is Changed
		Now checking particle content out of nucleus

	New Cross Section Variables
		deltaInvMass, Delta_pi_theta, Delta_pi_phi
		Their lateral error bands are not ready.
			I will implement in next version
	
	Many Plot and Table Improvements for Paper and Wine & Cheese 
	
	Collected all unique functions for paper plots under single file
		CCProtonPi0_Plotter_Paper.cpp
		CCProtonPi0_Plotter_Supplement.cpp
		Functions have duplicated code, due to unique requirements from Tony

	Added Missing Implementation Files to CVS
" .

cvs tag -F ${CCPROTONPI0_V} .

