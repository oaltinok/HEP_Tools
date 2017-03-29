CCPROTONPI0_V="v3_09_Thesis"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	No major changes

---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:
	Labels corrected for the cross-section plots

	W Shift to match data using QSq_Fitter

	Plots Improved for Signal Kinematics Cut
		Neutrino Energy XSec plot before and after FSI extended to 0 GeV
		W distribution extended to greater than 1.8

	GENIE Tuning as in Thesis
" .

cvs tag -F ${CCPROTONPI0_V} .

