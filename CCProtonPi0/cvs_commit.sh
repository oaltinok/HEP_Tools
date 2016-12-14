CCPROTONPI0_V="v3_01"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	No major changes

---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:
	More Control to GENIE Tuning
		Separated Delta Normalization Tuning from NonRES Normalization Tuning

	New Systematic Uncertainty
		Uncertainty for 2p2h events
			Using Phils 2D Gaussian Fit  DocDB 12424-v4
		Filling Universes for MC & Truth Histograms

	Non-Res Background Events classified using GENIE knobs
		truth_genie_wgt_Rvn1pi, truth_genie_wgt_Rvp1pi
		truth_genie_wgt_Rvn2pi, truth_genie_wgt_Rvp2pi

	Included NTuple Analysis for 2p2h Events
		Modified FluxReweighter
		Modified Interaction Type Plots

	Interaction Type Plots Modified
		Included 2p2h as a new interaction type
		Removed DIS with pion multiplicity, using only DIS
		This classification similar to Travis analysis

	Counter Class Improved
		Increment() now includes cvweight

	TruthAnalysis
		Added command protection
		Only available command is ./main.exe
		Protection against regular NTupleAnalysis commands such as ./main.exe run mc

" .

cvs tag -F ${CCPROTONPI0_V} .

