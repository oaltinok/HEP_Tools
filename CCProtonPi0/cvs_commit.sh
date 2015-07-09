cvs commit -m "v2_15
CCProtonPi0 Updates:
	No Major  Changes
--------------------------------------------------------------------------------
NTupleAnalysis Updates:
	Output Folder for ROOT files is now under /minerva/data
	New Mode: Reduce NTuples 
		Create a single ROOT file only with the events passing event selections
		Mode implemented as a Public Function for Analyzer Class
			All Classes initialized differently if the Mode is Reduce
				Other classes do not need to create their output files
		During Reduce, Filling Cut Statistics and Cut Histograms
	Analyze Mode Modified to use the Reduced NTuples" .
