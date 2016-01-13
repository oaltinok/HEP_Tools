CCPROTONPI0_V="v2_50"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	Background Type Improved
		NC
		CC-Anti Neutrino
		QE Like 
			No Pion, No other particles than Nucleons and muon
		Single Charged Pion
		Double Pions with Pi0
		Double Pions without Pi0
		Multiple Pions with Pi0
		Multiple Pions withotu Pi0
		Other
	
	Keeping Events with Michel Electron is Optional
		Decided not to remove these events in reconstruction stage

	Removed HoughEnergyLimit inside ConeBlobs
		PreFilterPi0 already checks max energy no need to check again
		
--------------------------------------------------------------------------------
NTupleAnalysis Updates:	
	BackgroundTool Updated for the new Background Types
	Plotter Improved
" .

cvs tag -F ${CCPROTONPI0_V} .
