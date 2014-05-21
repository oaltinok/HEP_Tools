# Packages for Minerva Neutrino Experiment

|Reconstruction | Data Analysis | Other Algorithms	|
|-------------- |-------------- |------------------	|
|CCProtonPi0	|NTupleAnalysis	| dEdX Improvement	|


* These packages are specialized for MINERvA Experiment data format.
	* Experiment URL: http://minerva.fnal.gov/
* Requires ROOT Data Analysis Framework
	* ROOT URL: http://root.cern.ch/drupal/

### CCDeltaPlus:

* Process raw data to create ntuples for a specific analysis focused on Delta Plus Production
* Final State Particles: muon, proton, and pi0

### NTupleAnalysis:

* Reads ntuples created by CCDeltaPlus

* Uses makeClass.cpp to get CCDeltaPlus Member Variables



### Other Algorithms:

##### dEdX_pID

* Improvement on pID Calculation: Algorithm improves the pID in case the particle trajectory overlaps with another trajectory.


  
