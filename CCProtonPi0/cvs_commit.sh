CCPROTONPI0_V="v3_02"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	No major changes

---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:
	Reweigthing only Non-RES n pi+

	MnvH1D array now includes signal types also
		Histogram indices are on Scratch page

	All Counters are using cvweight now (more accurate MC numbers)
		Cut Statistics
		Event Type Counters

	Study for Hadronic Invariant Mass (W)
		Breit-Wigner Fit on MATLAB for true delta signal events
" .

cvs tag -F ${CCPROTONPI0_V} .

