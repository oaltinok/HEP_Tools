CCPROTONPI0_V="v3_10"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	No major changes

---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:
	Complete GENIE Tuning
		Latest version of GENIE Tuning as invluded Aaron M. Thesis

	Fixed a Bug affecting the Error Summary Plots
		Background Subtracted Distributions had wrong plots but correct numbers
		Reason was MnvErrorBand CV histogram was not modified correctly

	Wine & Cheese and Paper Preparation
		New Plots and Studies
" .

cvs tag -F ${CCPROTONPI0_V} .

