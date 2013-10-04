Packages for Minerva Neutrino Experiment
==========================

Experiment URL: http://minerva.fnal.gov/

Notes:

	These packages are specialized for MINERvA Experiment data format.
	
DataAnalyze Algorithms: 


	2 Different data types
	ana: Analysis Output from a specific analysis package
	dst: Standard output file, can be generated from all analysis packages
	
	
	
	
Reconstruction Algorithms:

	
	dEdX Improved Algorithm for Particle ID calculation
	


DataAnalyze
==========================

ANA_CC
------------------------
Analysis package that is used for CCInclusive ana file output.

DST_CC
------------------------
Analysis package that is used with DST files.

ANA_PIDStudies
------------------------
Analysis package that is used for Particle Identification Studies


Reconstruction
==========================
Improvements on Minerva Reconstruction

dEdX_pID
------------------------
Improvement on pID Calculation: Algorithm improves the pID in case the particle track overlaps with another particle track.


Plotting
==========================
Note: After July 2013, I started to use ROOT for visualization instead of MATLAB, these packages will not be updated.

Visualization:

	x-y Scatter Plots
	
  	Single Histograms
  	
  	Data vs MC Comparison
  	
Data Statistics:

 	Chi-square Calculation
 	
  	radian to degree converter
  	
  
