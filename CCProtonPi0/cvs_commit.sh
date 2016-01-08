CCPROTONPI0_V="v2_48"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
 	Modified reconstructEvent() to save NTuple Variables during reconstruction
		Vertex Variables right after Vertex Reconstruction
		Muon Variables right after Muon Reconstruction
		Proton Variables right after Proton Reconstruction
		Pi0 Variables right after Pi0 Reconstruction
	This way we have the information of the particle even if the event rejected by a future stage
	
	New functions to test proton and pi0 reconstructions
		In some cases proton and pi0 momentums are -nan
		Reject these cases

--------------------------------------------------------------------------------
NTupleAnalysis Updates:	
	Improved Method for adding CutArrows on Plots 
		Automatically find arrow height
" .

cvs tag -F ${CCPROTONPI0_V} .
