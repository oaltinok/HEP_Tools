------------------------------------------------------------------------------------------
 Initialize the Class
------------------------------------------------------------------------------------------

root -l
.L ANA_CC.C
ANA_CC t

------------------------------------------------------------------------------------------
 Run Function
------------------------------------------------------------------------------------------

// To run all Playlist minerva1 ( 3,766,448 events )

t.run("/afs/fnal.gov/files/home/room2/oaltinok/PlayLists/v10r6p9/CCInclusive/ana/pl_MC_minerva1.dat", "MC_minerva1.root", "cutFile.txt","readme.txt");

// To run the test sample 

t.run("/afs/fnal.gov/files/home/room2/oaltinok/PlayLists/v10r6p9/CCInclusive/ana/ana_test.dat", "MC_test.root", "cutFile.txt","readme.txt");

------------------------------------------------------------------------------------------
 Plot Histograms from generated root file
------------------------------------------------------------------------------------------

t.plotHistograms("MC_minerva1.root","Plots");