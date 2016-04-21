CCPROTONPI0_V="v2_77"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	No Major Changes
---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:	
	Using PlotUtils/FluxReweighter for calculating Flux Weights and Systematics
		Removed old methods based on Trungs new_flux.h file
		Updated both Analyzer and TruthAnalyzer
		Did NOT updated the CrossSection Flux Integration yet. (Next release)
	Variable size binning for Q2 and Enu  
		Bins taken from the Carrie, Trung paper
" .

cvs tag -F ${CCPROTONPI0_V} .

