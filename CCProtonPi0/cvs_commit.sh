CCPROTONPI0_V="v2_73"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	Removed NuMIBRanches
		I am not using any of the variables
	Removed CommonPhysicsAnaBranches
		I am not using any of the variables
		It was causing a problem with <muon_theta> variable
			CCProtonPi0 <muon_theta> was getting overwritten by the 
				CommonPhysicsAnaBranches <muon_theta>
---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:	
	Improved Unfolding Study
		Added pi0_KE
		Plotting truth information for Training and Study Sample
		Added <CCProtonPi0_Plotter_UnfoldingStudy.cpp> to CVS
			I forgot adding it with the v2_72 
	
	Variable Bin Size is used in Cross Section Variables
		Bin low edges taken from Paper: PB015
		Class: BinList contains the constant arrays for the bin values

	Class: Plotter is improved to handle variable bin sized MnvH1Ds
		Most of the MnvPlotter Functions automatically handles variable bin size
		Modified custom plotting functions to handle variable bin size
			Normalize MnvH1Ds to bin width
			If bin width is constant, nothing changes
" .

cvs tag -F ${CCPROTONPI0_V} .

