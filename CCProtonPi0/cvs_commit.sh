CCPROTONPI0_V="v2_70"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	No Major Changes
---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:	
	Improved main.cpp design
		Moved main function definitions to src/Main_Functions.cpp
		Moved Minuit function definitions to src/Minuit_Functions.cpp
" .

cvs tag -F ${CCPROTONPI0_V} .

