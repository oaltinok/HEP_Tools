# Packages for Minerva Neutrino Experiment

|Reconstruction | Data Analysis | Other Algorithms	|
|-------------- |-------------- |------------------	|
|CCDeltaPlus	|NTupleAnalysis	| dEdX Improvement	|


* These packages are specialized for MINERvA Experiment data format.
	* Experiment URL: http://minerva.fnal.gov/
* Requires ROOT Data Analysis Framework
	* ROOT URL: http://root.cern.ch/drupal/

### CCDeltaPlus:

* Process raw data to create ntuples for a specific analysis
* Final State Particles: muon, proton, and pi0

### NTupleAnalysis:

* Reads ntuples created by varios Reconstruction Packages: 
	* CCDeltaPlus
	* CCInclusive
	* CCPi0
* Uses makeClass.cpp to create a new class for another ntuple



### Other Algorithms:

##### dEdX_pID

* Improvement on pID Calculation: Algorithm improves the pID in case the particle trajectory overlaps with another trajectory.


  
