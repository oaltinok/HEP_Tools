CCPROTONPI0_V="v2_43"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	Study: Neutrino Energy
		Built in functions and NTuple Variables
		Will separate implementation file and optional NTuple Variables in future releases
	VertexBlob() Function Revised
		Now it creates ONLY VertexBlob with radius 90mm
		Saves VertexBlob Visible Energy and Protonic Energy
	SaveExtraEnergy()
		After all particle reconstructions, this function collects all Unused Cluster Energy 
			in different volumes
		Radii for search: 50, 100, 150, 200, 300, 500 [mm]
	GetShortProtonCalConstant()
		Uses a line fit to get a calibration constant for given visible energy
			Not finalized, may change
		
--------------------------------------------------------------------------------
NTupleAnalysis Updates:	
	Particle Cannon
		Single Proton study to get Protonic Energy Constant 
	Study: Neutrino Energy
		Specific functions and histograms
		Will be removed in future releases
	Revised Interaction Class
		Histograms modified to study Event Kinematics 
" .

cvs tag -F ${CCPROTONPI0_V} .
