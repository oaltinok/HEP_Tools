CCPROTONPI0_V="v2_94"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	No Major Changes
---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:
	Added new Vertical Error Band
		MichelTrue: Events with True Michel Tag
			Their number depend on Michel Tool efficiency
		MichelFake: Events without True Michel Tag	
			Their number depend on Michel Tool purity

	Added new class Counter
		Improved existing struct counter

	Modified Background Subtraction Histograms
		Histogram range is limited to signal region
		Histogram filled after all event selections
		Error Bands also satisfy all event selections
		This modification fixed the problem with EM_EnergyScale Error Band

	Turned off Cross Section calculations for W
		I did not optimized W histograms for the XSec calculations
" .

cvs tag -F ${CCPROTONPI0_V} .

