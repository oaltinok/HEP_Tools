# Packages for Minerva Neutrino Experiment

|Statistical Data Analysis |MINERvA Operations| Other Algorithms	|
|----|-----|-----|
|CCProtonPi0	|CAMAC_DataReader|dEdX Improvement	|


* These packages are specialized for MINERvA Experiment data format.
	* Experiment URL: http://minerva.fnal.gov/
* Requires ROOT Data Analysis Framework
	* ROOT URL: http://root.cern.ch/drupal/

## Statistical Data Analysis:
### CCProtonPi0
* Two stage Data Analysis Package
	* Reconstruction Stage: 
		* Process raw data and MC simulation to detect signal
		* Create a smaller sample with processed information called NTuple
	* NTuple Analysis Stage: 
		* Process NTuples created in first stage for a detailed Data Analysis

* Signal Definition: 
	* Charged Current Neutrino Interaction inside Fiducial Volume
	* Final State Particles: muon, proton, and pi0

## MINERvA Operations:
### CAMAC_DataReader
* Package responsible for processing CAMAC Readout File
	* CAMAC Provides electronic signal output from the Test Beam Detector
	* Package reads the output data and converts scalars to real values
	* Updates the NearlineCurrentHistos.root File interactively
	* NearlineCurrentHistos.root file read by GMBrowser to show the plots to the Shifter
	
## Other Algorithms:

### dEdX_pID

* Improvement on pID Calculation: Algorithm improves the pID in case the particle trajectory overlaps with another trajectory.


  
