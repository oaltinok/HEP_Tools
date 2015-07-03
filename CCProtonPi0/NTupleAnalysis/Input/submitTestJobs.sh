CCPROTONPI0_V="dEdX_modified"
OUTPUTFOLDER="run"
OTHEROPTIONS=""

echo "Submitting Jobs ..."
echo "... NukeCCQETwoTrack Version: ${CCPROTONPI0_V}" 
echo "... Output Folder: ${OUTPUTFOLDER}"

runNumber_first=10200
runNumber_last=10214

for ((runNumber=${runNumber_first}; runNumber <= ${runNumber_last}; runNumber=runNumber+2 )); do
        echo "Submitting Run = ${runNumber}"
        $ana_script/ProcessAna.py --mc --run ${runNumber} --usecat --ana_tool NukeCCQETwoTrack  --inv eroica --use_jobsub_client --dcachein --outdir /minerva/data/users/oaltinok/NukeCCQETwoTrack/${CCPROTONPI0_V}/${OUTPUTFOLDER} ${OTHEROPTIONS}
done

echo "All Runs Submitted!"


