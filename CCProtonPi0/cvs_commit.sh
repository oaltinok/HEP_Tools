CCPROTONPI0_V="v2_92"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	Filling NTuple Variables for dEdX Uncertainties
		
---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:
	Neutrino Energy Calculation for events with proton CORRECTED!
		Summing Kinetic Energy of all protons 
		Calculation was wrong for more than 1 proton events
	
	Added Proton Systematics
		Lateral Error Bands from dEdX Uncertainties
			ProtonEnergy_BetheBloch
			ProtonEnergy_MassModel
			ProtonEnergy_MEU
			ProtonEnergy_Birks
		Vertical Error Bands	
			ProtonTracking

	Added Other Systematics
		Target Mass Uncertainty
		Background Constraint Uncertainty
		Muon Angle Uncertainty (Lateral Shift)

	Other Smaller Improvements
		Add Error Bands to Truth Tree — TruthAnalyzer
" .

cvs tag -F ${CCPROTONPI0_V} .

