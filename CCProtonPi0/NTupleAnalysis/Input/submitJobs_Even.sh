CCPROTONPI0_V="v2_55"
OUTPUTFOLDER="run"
OTHEROPTIONS=""

echo "Submitting Jobs ..."
echo "... CCProtonPi0 Version: ${CCPROTONPI0_V}" 
echo "... Output Folder: ${OUTPUTFOLDER}"

runNumber_first=10200
runNumber_last=10224

for ((runNumber=${runNumber_first}; runNumber <= ${runNumber_last}; runNumber=runNumber+2 )); do
        echo "Submitting Run = ${runNumber}"
        $ana_script/ProcessAna.py --mc --run ${runNumber} --usecat --ana_tool CCProtonPi0  --inv eroica --use_jobsub_client --dcachein --dcacheout --outdir /minerva/data/users/oaltinok/CCProtonPi0/MC/${CCPROTONPI0_V}/${OUTPUTFOLDER} ${OTHEROPTIONS}
done

echo "All Runs Submitted!"


