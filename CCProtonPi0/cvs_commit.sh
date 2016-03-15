CCPROTONPI0_V="v2_69"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	No Major Changes
---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:	
	SideBandTool Updated
		Using TMinuit to fit Side Bands simultaneously 
		TMinuit is used under main.cpp and SideBandTool is used to apply fit results and plottting
" .

cvs tag -F ${CCPROTONPI0_V} .

