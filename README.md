# Packages for Minerva Neutrino Experiment

|Statistical Data Analysis |MINERvA Operations| Other Algorithms|Setup|
|----|-----|-----|-----|
|CCProtonPi0	|CAMAC_DataReader	|dEdX Improvement	|Setup Scripts for CCProtonPi0	|
|NTuple Analysis|			|True Reco Correction	|				|
|		|			|dEdX Profile		|				|

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

### dEdX pID

* Improvement on pID Calculation: Algorithm improves the pID in case the particle trajectory overlaps with another trajectory.

### True Reco Correction
* In MC Analysis, most of the time the Reconstructed Parameter (Momentum, Energy, Angle, etc...) agrees with TRUE Parameter
* Some cases, TRUE vs RECO plot does not create a good x = y distribution, in such cases we apply a correction factor
* MATLAB Functions designed to find a linear correction factor for Reconstructed values

### dEdX Profile
* EM Shower Profile for Different Shower Lengths
* Input File:
	* Each row is unique to a single EM Shower
	* Each column is total cluster energy in a plane
* Output:
	* EM Shower Profile Plots for different Number of Planes
* Algorithm:
	1. Read Input File
	2. Group EM Showers by Plane Numbers 
		* 4-Plane Showers
		* 5-Plane Showers
		* ..
		* 30-Plane Showers
		* etc..
	3. Normalize EM Showers to the MAX Cluster energy in each Group
		* Each group is normalized uniquely
	4. Plot Each Group Separately
		* Prints out plot for each group

## Setup:

### Setup Scripts for CCProtonPi0

* CCProtonPi0 Package requires MINERvA Software Framework to run
* These setup scripts make sure that the user have the correct environmental variables
  
