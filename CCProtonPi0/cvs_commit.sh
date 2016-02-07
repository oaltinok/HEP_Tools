CCPROTONPI0_V="v2_55"
cvs commit -m "${CCPROTONPI0_V}
CCProtonPi0 Updates:
	Added P and Theta for Signal FS particles Truth Variables

-------------------------------------------------------------------------------
NTupleAnalysis Updates:	
	New Mode:
		Calculate Cross Section
			./main.exe calc
	New Class:
		CCProtonPi0_CrossSection
			Class responsible for Cross Section Calculation
			Initial commit for the Class
				No methods implemented yet!

	Folder structure improved
		Forming final folder structure inside Libraries/Folder_List.h
		Other Objects gets the final folder structure from Folder_List

	Added new histograms for the variables we will calculate the cross section
		Efficiency
		Response
" .

cvs tag -F ${CCPROTONPI0_V} .

