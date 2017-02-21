CCPROTONPI0_V="v3_08"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	No major changes

---------------------------------------------------------------------------------------------------
NTupleAnalysis Updates:
	CC-RES Suppression Event Reweighing
		Turned off by default

	Plot styles exactly matched with 
        PRD 94, 052005 (2016) 

	TruthAnalyzer
		BeforeFSI and After FSI using finer binning
" .

cvs tag -F ${CCPROTONPI0_V} .

