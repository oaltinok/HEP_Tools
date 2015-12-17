CCPROTONPI0_V="v2_44"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	Removed ALL Neutrino Energy Study Functions

	Revised Vertex Energy Calculation
		New calculation based on CCProtonPi0 MC study
		Line Equation to estimate Protonic Energy for a given Visible Energy around vertex
	
	Improved Nuclear Binding Energy
		Using Nuclear Binding Energy per Nucleon for Carbon
		Using expected Number of Nucleons for 1 Track and 2 Track
	
--------------------------------------------------------------------------------
NTupleAnalysis Updates:	
	FindPlaylistFiles.sh
		Fixed Ana Directory location for data events

	Removed ALL Neutrino Energy Histograms and Functions

	Reorganized Plotter Class
		Different Functions for DatavsMC & MCOnly Plots
" .

cvs tag -F ${CCPROTONPI0_V} .
