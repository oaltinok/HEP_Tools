CCPROTONPI0_V="v2_09_AdvancedTracker"
OUTPUTFOLDER="run"
OTHEROPTIONS=""

echo "Submitting Jobs ..."
echo "... CCProtonPi0 Version: ${CCPROTONPI0_V}" 
echo "... Output Folder: ${OUTPUTFOLDER}"

runNumber_first=13201
runNumber_last=13225

for ((runNumber=${runNumber_first}; runNumber <= ${runNumber_last}; runNumber=runNumber+2 )); do
        echo "Submitting Run = ${runNumber}"
        /afs/fnal.gov/files/home/room2/oaltinok/User_Area/Minerva_nightly/Tools/ProductionScripts/ana_scripts/ProcessAna.py --mc --run ${runNumber} --usecat --ana_tool CCProtonPi0 --os=SL6 --pause=0 --inv resurrection --use_jobsub_client --dcachein --outdir /minerva/data/users/oaltinok/CCProtonPi0/MC/${CCPROTONPI0_V}/${OUTPUTFOLDER} ${OTHEROPTIONS}
done

echo "All Runs Submitted!"


