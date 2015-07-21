cvs commit -m "v2_18
CCProtonPi0 Updates:
	Changed Global Variables to mutable Member Variables
	Removed Unused/Old Classes
	Signal Definition Changed
		No Longer Requiring 2x Gamma Out of Nucleus
		Current Definition:
			CC Neutrino Interaction
			Single Pi0 out of nucleus
			No Mesons out of nucleus
	VertexBlob() Modified
		FilamentVertex Calculation removed -- Trung's suggestion
--------------------------------------------------------------------------------
NTupleAnalysis Updates:
	Interaction Class Modified to use MnvH1D 
	Plotter Class Improved
		Implementation for Plotting Macros in another file
	CCProtonPi0_Plotter_Macros.cpp
" .
