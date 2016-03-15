CCPROTONPI0_V="v2_68"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	No Major Changes
---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:	
	Many Design Improvements

	data_POT & mc_POT are now constants in Base Class (CCProtonPi0_NTupleAnalysis)

	CrossSection Class Updates
		Using struct to contain all required data structures for a cross section variable
		Flux integration corrected using minimum and maximum Enu
		
	Binning information for Cross Section Variables centralized
		BinList class contains binning information for all cross section variables
		TruthAnalysis and Particle Classes uses BinList to get binning
" .

cvs tag -F ${CCPROTONPI0_V} .

