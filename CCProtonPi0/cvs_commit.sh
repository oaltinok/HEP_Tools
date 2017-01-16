CCPROTONPI0_V="v3_03"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	No major changes

---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:
	W Fit using TMinuit

	MaRES Fit using QSq Distribution
		New Class: QSqFitter
		Fill 101 Universes for each Up and Down Shift
			Range up to 2 sigma shift
		Get Event Weights by partitioning
			0 to +1 sigma
			+1 to +2 sigma
			0 to -1 sigma
			-1 to -2 sigma
		Find Lowest ChiSq for Best MaRES Value
" .

cvs tag -F ${CCPROTONPI0_V} .

