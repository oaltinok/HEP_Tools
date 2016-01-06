CCPROTONPI0_V="v2_47"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	No Major Changes
	
--------------------------------------------------------------------------------
NTupleAnalysis Updates:	
	CutList Class Improved
		Cut Table is created for all events
			Overall Efficiency and Purity is important
		Keeping 1Track and 2Track Cut Tables for internal studies
" .

cvs tag -F ${CCPROTONPI0_V} .
