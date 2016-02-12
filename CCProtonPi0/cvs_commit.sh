CCPROTONPI0_V="v2_59"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	No Major Changes
-------------------------------------------------------------------------------
NTupleAnalysis Updates:	
	Systematics in Data Corrected
		Filling Only CV after Loop
	
	Systematics in MC Corrected
		MuonTracking Systematic Corrected 
			Using Eroica, minerva1
		Flux Systematic Corrected
			Using Latest Flux
	
	CrossSection Plots Corrected
		No longer using POT Ratio for plotting POT Normalized 
		Calculation already includes POT

	Added Systematics to TruthAnalyzer
		Flux & GENIE Only
		Other Systematics are CV

	Plotter Improved for handling many plots  of CrossSection Calculations
		Improved folder structure for output plots
" .

cvs tag -F ${CCPROTONPI0_V} .

