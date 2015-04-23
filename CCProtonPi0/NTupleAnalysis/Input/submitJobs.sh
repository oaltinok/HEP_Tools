CCPROTONPI0_V="v2_08"
OUTPUTFOLDER="run"
OTHEROPTIONS="--dst"

echo "Submitting Jobs ..."
echo "... CCProtonPi0 Version: ${CCPROTONPI0_V}" 
echo "... Output Folder: ${OUTPUTFOLDER}"

RunNumberList=(13200
    13201
    13202
    13203
    13204
    13205
    13206
    13207
    13208
    13209
    13210
    13211
    13212
    13213
    13214
    13215
    13216
    13217
    13218
    13219
    13220
    13221
    13222
    13223
    13224
    13225)

for runNumber in ${RunNumberList[@]}; do
        echo "Submitting Run = ${runNumber}"
        /afs/fnal.gov/files/home/room2/oaltinok/User_Area/Minerva_nightly/Tools/ProductionScripts/ana_scripts/ProcessAna.py --mc --run ${runNumber} --usecat --ana_tool CCProtonPi0 --os=SL6 --inv resurrection --use_jobsub_client --dcachein --outdir /minerva/data/users/oaltinok/CCProtonPi0/MC/${CCPROTONPI0_V}/${OUTPUTFOLDER} ${OTHEROPTIONS}
done

echo "All Runs Submitted!"

# 13200
# 13201
# 13202
# 13203
# 13204
# 13205
# 13206
# 13207
# 13208
# 13209
# 13210
# 13211
# 13212
# 13213
# 13214
# 13215
# 13216
# 13217
# 13218
# 13219
# 13220

