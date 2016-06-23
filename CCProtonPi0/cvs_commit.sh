CCPROTONPI0_V="v2_86"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	No Major Changes
---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:
	SideBandTool Improvements
		Writing Multi-Universe Fit Results to a Text File
			Input/BckgConstraints/Weights_All_Universes.txt
		Writing Statistics for each Side Band
			nData vs nMC etc…
		Calculating ChiSq before fit to compare with minimized ChiSq
	
	Analyzer Improvements
		Added new functions to fill Vertical Error Bands Manually
			Tested results for Auto Fill vs Manual Fill
			Using same code to fill error band manually, this allows me to change weight for each universe

		Applying unique Background Constraints to each universe
			Reads Weights from Input/BckgConstraints/Weights_All_Universes.txt
" .

cvs tag -F ${CCPROTONPI0_V} .

