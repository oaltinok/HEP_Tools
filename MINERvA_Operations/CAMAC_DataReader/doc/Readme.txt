========================================================================
Working Principal:
========================================================================
There is an INFINITE loop which does the following:
	Check /work/data area for new CAMAC readout file
		If no NEW files:
			Check Again
		If there is a NEW file :
			Read File	
			Process Data	
			Update /minerva/data/testbeam2/nearonline/CAMACDataHistos.root

========================================================================
Usage:
========================================================================
1) Open a NEW Terminal
2) Navigate to package folder:
	cd /home/nfs/minerva/CAMAC_DataReader
3) Run the package
	./CAMACDataReader
4) Stop Execution
	CTRL + C  or Close Terminal


========================================================================
Data File Locations - provided by Jeffrey and Dan
========================================================================
The files are on the testbeam computers, mnvtb01-05. Logon with "shh minerva@mnvtb04" on a minerva machine for example.

------------------------------------------------------------------------
The GMBrowser Config File:
------------------------------------------------------------------------
/home/nfs/minerva/cmtuser/Minerva_v10r9p1/Tools/ControlRoomTools/gmbrowser/nearline.cfg

------------------------------------------------------------------------
CAMACDataHistos.root File:
------------------------------------------------------------------------
/minerva/data/testbeam2/nearonline/CAMACDataHistos.root

Data is written to either mnvtb04 or mnvtb03 depending on what system it is. Every 12 hours it's backed up to somewhere inside /minerva/data/testbeam2.

------------------------------------------------------------------------
Subrun by Subrun CAMAC Readout: (Input for CAMAC_DataReader)
------------------------------------------------------------------------
/work/data/TB_XXXXXXXX_YYYY_cosmc_v09_ZZZZZZZZZZ_camac.dat
                Where XXXXXXXXXXX is typcal run number
                Where YYYY is typical subrun
                Where ZZZZZZZZZZZZZ is the timestamp

------------------------------------------------------------------------
Spill by Spill CAMAC Readout:
------------------------------------------------------------------------
/home/nfs/minerva/daq/daqdata/lastspill_camac.dat

========================================================================
Data File Format - provided by Dan
========================================================================
The file format is as follows:
 
Line1 : version
Line2 : Variable names
Line3 : Conversion factors to convert counts to ps or other units (Designed so the numbers are just a multiplicative factor)
Line4 : units of the variable
 
After this point we list the histograms, graphs, and frequency plots.
 
Format for these:
 
h = histogram line
g = graph line
f = frequency plot
 
histograms are: h title variable nbins low high
graphs are : g title xvariable yvariable
frequency are : f title array of all variables. Plot the variables with 1 ignore variables with 0.


