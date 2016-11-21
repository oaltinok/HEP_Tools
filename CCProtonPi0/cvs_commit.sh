CCPROTONPI0_V="v3_00"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	No major changes

---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:
	GENIE Tuning Corrections
		NonRes1pi modification applied to both Rvn1pi and Rvp1pi
		Uncertainties for Rvn1pi and Rvp1pi
		XSec Variables with different cvweights removed
			Hist filling only with different cvweight does not make sense without Side Band Fit

	GENIE Tuning Modifications
		No longer using updated MaRES and CCRES Normalization for event weights
		No longer using updated MaRES and MvRES uncertainties 
		No longer using updated Rvn1pi and Rvp1pi uncertainties 

	GENIE Tuning Study Plots 
		Plots to compare Nominal vs Tuned in Data and MC
" .

cvs tag -F ${CCPROTONPI0_V} .

