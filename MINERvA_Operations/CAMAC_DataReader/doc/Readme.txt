########################################################################
#																	   #
#					Working Principal and Usage						   #
#																	   #
########################################################################

There are two modes that you can run the CAMAC_DataReader
========================================================================
Mode: Auto
========================================================================
Automatically reads the camac readout spill by spill and updates a single
	ROOT file which is read by GMBrowser to show Plots
See Data File Locations for input/output files.

Usage:
	1) Open a NEW Terminal
	2) Navigate to package folder:
		cd /home/nfs/minerva/CAMAC_DataReader
	3) Run the package without any arguments
		./CAMACDataReader
	4) Stop Execution
		CTRL + C or Close Terminal

Algorithm Flow:
	Infinite Loop that does the following:
		Check lastspill_camac.dat file
			If it is the latest
				Read File	
				Process Data	
				Update CAMACDataHistos.root
			else 
				Check again
			

========================================================================
Mode: Manual
========================================================================
Reads the camac output file for specified run and subrun numbers
	Creates a unique .ROOT file for that run/subrun
See Data File Locations for input/output files.

Usage:
	1) Open a NEW Terminal
	2) Navigate to package folder:
		cd /home/nfs/minerva/CAMAC_DataReader
	3) Run the package with the following arguments
		./CAMACDataReader -r <run_number> -s <subrun_number>

Algorithm Flow:
	Process run/subrun numbers locate the input file
	Read camac output file
	Process Data
	Create/Update .root file specific to the provided rub/subrun


########################################################################
#																	   #
#						Data File Locations							   #
#																	   #
########################################################################

The files are on the testbeam computers, mnvtb01-05. 
Logon with "shh minerva@mnvtb04" on a minerva machine for example.

========================================================================
INPUT Data Files:
========================================================================
Mode: Auto - Spill by Spill camac output
	/home/nfs/minerva/daq/daqdata/lastspill_camac.dat
	
Mode: Manual - Single file for a specific Subrun
	/work/data/TB_XXXXXXXX_YYYY_cosmc_v09_ZZZZZZZZZZ_camac.dat
                Where XXXXXXXXXXX is typical run number
                Where YYYY is typical subrun
                Where ZZZZZZZZZZZZZ is the timestamp

========================================================================
OUTPUT Data Files:
========================================================================
Mode: Auto - Single .ROOT file which is read by GMBrowser
	/minerva/data/testbeam2/nearonline/CAMACDataHistos.root

Mode: Manual - Each Subrun has its unique .ROOT File
	/minerva/data/testbeam2/camacdata/TB_XXXXXXXX_YYYY_ZZZZZZZZZZ_camac.root
                Where XXXXXXXXXXX is typical run number
                Where YYYY is typical subrun
                Where ZZZZZZZZZZZZZ is the timestamp


========================================================================
GMBrowser Files:
========================================================================
Config File:
	/home/nfs/minerva/cmtuser/Minerva_v10r9p1/Tools/ControlRoomTools/gmbrowser/nearline.cfg
	
Plot Macro File:
	/home/nfs/minerva/cmtuser/Minerva_v10r9p1/Tools/ControlRoomTools/gmbrowser/macros/TestBeam_CAMACData.C



########################################################################
#																	   #
#						Input Data File Format						   #
#																	   #
########################################################################

Thanks to Dan for this information...

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


