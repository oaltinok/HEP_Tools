CCPROTONPI0_V="v3_07"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	No major changes

---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:
	Thesis Style Plots
		Major changes to Plot Styles
			Square Plots
			Larger Fonts
			No Statistics Box
			No unnecessary text on plots		
			Corrected Axis Titles
			etc..

	Studies Cont’d
		ROOT Fitter for QSq Enu Fits
			Thanks to Trung
	
		Disabled all QSq Study Functions
" .

cvs tag -F ${CCPROTONPI0_V} .

