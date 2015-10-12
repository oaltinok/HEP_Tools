cvs commit -m "v2_38
CCProtonPi0 Updates:
	HTBlob updated for Calorimetric Energy
		New Constants and Method to calculate Calorimetric Energy from Visible Energy
		getBlobEnergyTime()
			function return type changed to void 
				function never returns FAILURE, no need for returning StatusCode
			Make two different functions for OLD and NEW method
				NEW Method will be revised -- It is not finalized
	
--------------------------------------------------------------------------------
NTupleAnalysis Updates:
	Removed EM Shower Energy Study Histograms and Functions
	Removed unused playlists from Input/Playlist Folder
	Revised getPOT() Function
	Reduce mode no longer overwrites the ReducedNTuple.root and CutHistograms.root
		Protection against trying to reduce by mistake
		Now user needs to keep track of ROOT file names" .
