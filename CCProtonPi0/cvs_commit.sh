CCPROTONPI0_V="v2_83"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	No Major Changes
---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:
	New Background Constraints with repeated Side Band Study
		
	SideBandTool Improvements
	 	Included Cross Section Variables
		Now using only pi0_InvMass as fit parameters
		Design changes to accommodate new XSec Variables and Plotting
" .

cvs tag -F ${CCPROTONPI0_V} .

