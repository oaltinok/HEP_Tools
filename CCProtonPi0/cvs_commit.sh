CCPROTONPI0_V="v2_51"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	Improved ConeBlobs Function
		More modular structure, collected common task items under functions such as
			FillUsableClusters()
			ProcessRejectedClusters()
			RecoverShowers_InvMass()
			etc..
		
		Removed old and unnecessary NTuple Variables
		
		Improved Cluster Rejection Method
			Collecting all rejected clusters in a single blob, 
				Previously out time and low charge was different
			Rejecting following clusters:
				Out Time with Muon (OLD)
				Low Charge Cluster (OLD)
				Cluster inside HCAL (NEW)
		
		Improved Recover Methods
			Recover 1Shower Events
				Not recovering anything yet, saving 1Shower Information only
			Recover 3Shower Events
				Check position and direction of each shower, 
					if one fails reject it and keep the other two
				
-------------------------------------------------------------------------------
NTupleAnalysis Updates:	
	No Major Changes
" .

cvs tag -F ${CCPROTONPI0_V} .
