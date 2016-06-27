CCPROTONPI0_V="v2_87"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	No Major Changes
---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:
	CrossSection Improvements
		Writing Logs as a text file output instead of terminal output
			Text file used for Multi-Universe Calculations
				Unfolding, Efficiency Correction etc. still uses Terminal output
			Log Files created under Output/TextFiles
		Removed Function Calc_NormalizedNBackground_TFractionFitter()
		Calc_NormalizedNBackground() Improved for Multi-Universe setup
		Subtract_Background() Improved for Multi-Universe setup

	Plotter Improved to plot DataMCWithErrorBands()

" .

cvs tag -F ${CCPROTONPI0_V} .

