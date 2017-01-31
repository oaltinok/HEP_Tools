CCPROTONPI0_V="v3_06"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	No major changes

---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:
	QSq Studies Completed
		MaRES Fit handled by Vertical Error Bands
			Commented out in this version
		DeltaFactor Fit handled by Vertical Error Bands
			Active in this version
" .

cvs tag -F ${CCPROTONPI0_V} .

