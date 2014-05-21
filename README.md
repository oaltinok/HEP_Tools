# Packages for Minerva Neutrino Experiment

|Reconstruction | Data Analysis | Other Algorithms	|
|-------------- |-------------- |------------------	|
|CCProtonPi0	|NTupleAnalysis	| dEdX Improvement	|


* These packages are specialized for MINERvA Experiment data format.
	* Experiment URL: http://minerva.fnal.gov/
* Requires ROOT Data Analysis Framework
	* ROOT URL: http://root.cern.ch/drupal/

### CCProtonPi0:

* Process raw data to create ntuples for a specific analysis focused on Charged Current Proton Pi0 Final States
* Final State Particles: muon, proton, and pi0

### NTupleAnalysis:

* Reads ntuples created by CCProtonPi0

* Uses makeClass.cpp to get CCProtonPi0 Member Variables



### Other Algorithms:

##### dEdX_pID

* Improvement on pID Calculation: Algorithm improves the pID in case the particle trajectory overlaps with another trajectory.


  
