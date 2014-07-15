cvs commit -m "v1_04 - 2014_07_15
Revised Pi0 Reconstruction
	Global Variables Used in all functions
	Removed unused functions and variables
	ConeBlobs()
		Variable Naming Match
		Return type changed: StatusCode to bool
			Returns false if setPi0ParticleData() fails
		ConeBlobs() main function that controls Pi0 Reconstruction. 
			If it fails, the  reconstructEvent() for that event stops.
	VtxBlob()
		Return type changed: StatusCode to bool
			Always returns true (return type reserved for future implementation)
	processBlobs()
		Removed unused variables
		Return type changed: StatusCode to void
		
Options File Modifications
	New Options files for DEBUG
	Original options file set to INFO" .