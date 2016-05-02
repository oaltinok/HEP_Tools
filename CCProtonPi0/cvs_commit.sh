CCPROTONPI0_V="v2_78"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	Calculations for Neutrino Energy, QSq, WSq are removed from reconstruction stage
		They are calculated ruing Analysis stage after EM Energy correction
	Other small design improvements

---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:	
	Background Constraints Included
		Side Band study results
	
	Event Kinematics (Enu, QSq, WSq, W) calculations added
	
	Improved Signal Characteristics (Signal Type)
		Signal type determined from Event Record (inside nucleus)
			QE
			RES_1232, RES_1535, RES_1520, RES_Other
			DIS_1_pi, DIS_2_pi, DIS_Multi_pi, DIS_Other
		Check Enu, QSq, W for
			All Signal Events
			MINOS Matched Signal Events
			Selected Signal Events
	
	Improved FluxReweighter Usage
		No longer assuming first playlist is the minerva1

	Other small design improvements
" .

cvs tag -F ${CCPROTONPI0_V} .

