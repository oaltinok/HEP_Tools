cvs commit -m "v2_19
CCProtonPi0 Updates:
	Signal Definition Changed
		Signal Definition explicitly checks for Fiducial Volume
			Actually True Vertex Fiducial Volume is required for ALL Events,
			This addition just for a better Signal Definition
	Background Study Improved: Two Different Background Sets
		Background Type
			AntiNeutrino, QELike, SinglePion, DoublePion, MultiPion
		Background with Pi0
			NoPi0, SinglePi0, MultiPi0
	Do NOT Reconstruct MC Event
		If TRUE vertex is NOT Fiducial
		If TRUE Interaction is NOT ChargedCurrent 
--------------------------------------------------------------------------------
NTupleAnalysis Updates:
	makeClass Folder Structure Improved for different TBranches
		CCProtonPi0 , Truth, Other
		Truth Branch is used to Analyze All Truth Events
			Get Total Number of Signal Events    
" .
