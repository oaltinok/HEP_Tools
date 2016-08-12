CCPROTONPI0_V="v2_93"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	Signal Definition Changed
		Added Neutrino Energy Requirement
			1.5 GeV < E_nu < 20 GeV
			Previously, I was applying this requirement in NTupleAnalysis stage
		Added W Constraint
			True W_experimental < 1.8 GeV
	
	Added Hadronic Recoil Energy
		Saving calibrated energy for non-muon clusters
	
	Calculating True Experimental QSq and W during Ana Stage
		
		
---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:
	Cross Section Calculations
		Corrected Flux Histogram Usage
			Now adding Flux Error band and rebinning

	Removed UpdateSignal Definition
		No longer required — Neutrino Energy Cut applied in Ana Stage

	 Plotter: Added Closure Test for Cross Section Calculations
		Processing MC as Data and checking ALL Universes bin by bin

	Bug Fixes
		Corrected muon_theta fill on Lateral Error Bands
" .

cvs tag -F ${CCPROTONPI0_V} .

