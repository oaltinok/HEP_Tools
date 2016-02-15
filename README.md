# Packages for Minerva Neutrino Experiment

* These packages are specialized for MINERvA Experiment data format.
	* Experiment URL: http://minerva.fnal.gov/
* Requires ROOT Data Analysis Framework
	* ROOT URL: http://root.cern.ch/drupal/

## Statistical Data Analysis
* Two Stage Data Analysis package for an exclusive channel on neutrino interactions
* Signal Definition:
	* Charged Current Neutrino Interaction inside Fiducial Volume
	* Final State Particles: muon, proton, and pi0
	
### First Stage: CCProtonPi0
* Processes raw data using many other packages on MINERvA Framework
* Reconstructs final state particles and calculates their estimated kinematics
* Identifies Events who are likely to be Signal
* Outputs a smaller data sample with processed information called NTuple

### Second Stage: CCProtonPi0/NTupleAnalysis
* Reduce Initial Data to Final 
	* ./main.exe reduce data
	* Selects best Signal Candidates from Signal-Likely Events
* Analyze Reduced Data Sample 
	* ./main.exe run data
	* Analyzes reduced data sample and create histograms for further study
* Calculate Cross Section
	* ./main.exe calc
	* Calculates Cross Section of the CCProtonPi0 Neutrino Interaction
	* Estimates Statistical & Systematic Errors of the final result

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
  
