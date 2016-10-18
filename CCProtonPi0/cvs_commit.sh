CCPROTONPI0_V="v2_98"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	Signal Before FSI Included
		Tagging signal before FSI using GENIE event record
		Saving cross section variables before FSI
			Muon momentum and angle
			Pion momentum, kinetic energy, and angle
			QSq, W, and neutrino energy

---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:
	Before FSI Addition
		Cross Section Calculation
		Comparison Plots for:		
			Data vs GENIE w/ FSI and GENIE w/o FSI
	
	FSI Type Addition
		FSI Type determined during TruthAnalysis
			I may move this to ana stage for performance
			This version will keep it in TruthAnalysis stage
		Cross Section Calculation
		Comparison Plots for:
			Data vs different FSI Types

	Interaction Type Addition
		Cross Section Calculation
		Comparison Plots for:
			Data vs different Interaction Types

	Folder Structure Changed
		Bluearc (/minerva/data) is going away
		New output:  /pnfs/minerva/persistent/	

	Particle Cannon 
		Generated and Analyzed pi0 in MINOS-like system

" .

cvs tag -F ${CCPROTONPI0_V} .

