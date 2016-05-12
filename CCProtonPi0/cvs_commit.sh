CCPROTONPI0_V="v2_79"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	No Major Changes
---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:
	SideBandTool Improved
		Using MnvH1Ds instead of TH1Ds
			MnvH1Ds contains errors
			TH1Ds are created inside Plot function
		
		Added DrawDataMCWithErrorBand function
		
		Removed a bug which connects TFile and MnvH1Ds in different SideBands
			Removed the bug using h->SetDirectory(NULL);	

" .

cvs tag -F ${CCPROTONPI0_V} .

