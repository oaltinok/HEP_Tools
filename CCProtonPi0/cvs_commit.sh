CCPROTONPI0_V="v2_67"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	Installed HEAD version of Rec/ProngMaker
		Need to keep in cmtuser area for HEAD version of Michel Tool

	Michel Electron Search Updates
		Michel Variables controlled by a new study
			Use StudyMichelElectron on options file
		Added Truth Match to found michelProng
		Separated Vertex Michel and Track End Point Michel Variables
	
-------------------------------------------------------------------------------
NTupleAnalysis Updates:	
	New Plots & Variables for Michel Study 
" .

cvs tag -F ${CCPROTONPI0_V} .

