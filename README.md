# Packages for Minerva Neutrino Experiment

* These packages are specialized for MINERvA Experiment data format.
	* Experiment URL: http://minerva.fnal.gov/
* Requires ROOT Data Analysis Framework
	* ROOT URL: http://root.cern.ch/drupal/

## Statistical Data Analysis
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

## MATLAB Functions for specific studies
* True Reco Correction
* Neutrino
	* Nuclear Binding Energy
* Pi0
	* dEdX Profile
	* Gamma Energy Correction
	* Invariant Mass
	* Side ECAL Study
* Proton
	* Momentum Correction
	* Short Proton Energy

## Nearline Online Data Processing
* CAMAC_DataReader
	* Package responsible for processing CAMAC Readout File
		* CAMAC Provides electronic signal output from the Test Beam Detector
		* Package reads the output data and converts scalars to real values
		* Updates the NearlineCurrentHistos.root File interactively
		* NearlineCurrentHistos.root file read by GMBrowser to show the plots to the Shifter
	
## Other Algorithms for MINERvA Framework
* dEdX pID

## Setup

* Setup Scripts for CCProtonPi0
	* CCProtonPi0 Package requires MINERvA Software Framework to run
	* These setup scripts make sure that the user have the correct environmental variables
  
