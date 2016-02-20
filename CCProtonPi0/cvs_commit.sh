CCPROTONPI0_V="v2_60"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	Corrected Q2 Calculation
		Using 4 Momentum of the Beam 
	Signal Definition Changed
		Added Max Neutrino Energy 20 GeV

-------------------------------------------------------------------------------
NTupleAnalysis Updates:	
	CrossSection Class Design Improved
		MC Processing is optional and it will create another file under
			data/NTupleAnalysis/MC/Analyzed/
	
	Plotter Improved 
		Cross Section Comparison Plots
		Cross Section Check Plots
" .

cvs tag -F ${CCPROTONPI0_V} .

