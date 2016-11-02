CCPROTONPI0_V="v2_99"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	No major changes

---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:
	File Structure Changed — Again!
		pnfs is creating problems for analysis output
		/minerva/data seems to stay for analysis output for some more time
		I revert back to use /minerva/data
	
	Corrected Flux Integration Stage
		Scaling Flux Universes by their own Integrals

	GENIE Tuning Included
		GENIE needs to be modified using Aaron M. and Phil R. studies.
		Included the updated event weights and systematics for
			MaRES, MvRES, Rvn1pi, Rvp1pi

	GENIE Tuning Special Plots
		Collected under Plotter/CCProtonPi0_Plotter_GENIE_Tuning.cpp
" .

cvs tag -F ${CCPROTONPI0_V} .

