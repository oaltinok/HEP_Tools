CCPROTONPI0_V="v2_45"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	Latest version of Neutrino Energy Calculations
		Vertex Energy and Extra Energy study finalized
		Binding Energy Removed
			Vertex Energy Calibration includes binding energy
	
--------------------------------------------------------------------------------
NTupleAnalysis Updates:	
	Neutrino Energy Study is Finalized
		New Functions & Histograms
	Plotter Class Improvements
		POT Normalization & Area Normalization Plots added
" .

cvs tag -F ${CCPROTONPI0_V} .
