CCPROTONPI0_V="v2_91"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	Function to Save Trajectory Info inside the Detector
		Copied from CCPi0AnaTool

	Added NTuple Variables for dEdX Uncertainties
		Not filling them in this version
---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:
	New Vertical Error Bands
		PionResponse
			Weights for SingleChargedPion_ChargeExchange
		NeutronResponse
			Weights for Neutron Inelastic  Scattering

	Removed Old/Unused Implementations	
" .

cvs tag -F ${CCPROTONPI0_V} .

