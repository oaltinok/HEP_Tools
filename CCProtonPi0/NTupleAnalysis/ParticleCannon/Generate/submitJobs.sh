# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   Custom Functions
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

STAGE=${1}
ISTEST=false

if ${ISTEST};then
    OUTPUTFOLDER="test"
    OTHEROPTIONS=""
    run_first=30000
    run_last=30000
    subrun_first=1
    subrun_last=5
else
    OUTPUTFOLDER="Steel_2"
    OTHEROPTIONS=""
    run_first=10000
    run_last=10004
    subrun_first=1
    subrun_last=200
fi


MAINDIR="/pnfs/minerva/persistent/users/oaltinok/ParticleCannon"
#VERTEX="--name primary --time 1000 --zCenter 9183 --deltaZ 6 --radius 10"
VERTEX="--name primary --time 1000 --zCenter 9161 --deltaZ 11 --radius 10"
#PI0="--name pi0 --flux 111 --minp 0.0 --maxp 1.5 --theta 0 --phi 0 --vertex primary"
PI0="--name pi0 --flux 111 --momentum 2.0 --theta 0 --phi 0 --vertex primary"
GAMMA="--name gamma --flux 22 --momentum 0.8 --theta 0 --phi 0 --vertex primary"
MUON="--name muon --flux 13 --minp 1 --maxp 10 --theta 0 --phi 0 --vertex primary"
PROTON="--name proton --flux 2212 --minp 0.0 --maxp 0.6 --thetamin 0.0 --thetamax 3.14 --phi 0 --vertex primary"

run_cal_pc()
{
    for ((run=${run_first}; run <= ${run_last}; run=run+1 )); do
            echo "Submitting Run = ${run}"
            python ${mc_script}/ProcessMC.py --${STAGE} -r ${run} -f ${subrun_first} -l ${subrun_last} --vertex "${VERTEX}" --particle "${PI0}" --outdir ${MAINDIR}/${OUTPUTFOLDER} ${OTHEROPTIONS}
    done
}

run_minos()
{
    for ((run=${run_first}; run <= ${run_last}; run=run+1 )); do
           echo "Submitting Run = ${run}"
           python ${mc_script}/ProcessMC.py --${STAGE} -r ${run} -f ${subrun_first} -l ${subrun_last}  --indir ${MAINDIR}/${OUTPUTFOLDER} --outdir ${MAINDIR}/${OUTPUTFOLDER} ${OTHEROPTIONS}
    done
}

run_reco()
{
    for ((run=${run_first}; run <= ${run_last}; run=run+1 )); do
           echo "Submitting Run = ${run}"
           python ${mc_script}/ProcessMC.py --${STAGE} --no_minos -r ${run} -f ${subrun_first} -l ${subrun_last} --indir ${MAINDIR}/${OUTPUTFOLDER} --outdir ${MAINDIR}/${OUTPUTFOLDER} ${OTHEROPTIONS}
    done
}

run_dst()
{
    for ((run=${run_first}; run <= ${run_last}; run=run+1 )); do
           echo "Submitting Run = ${run}"
           python ${mc_script}/ProcessMC.py --${STAGE} -r ${run} -f ${subrun_first} -l ${subrun_last} --indir ${MAINDIR}/${OUTPUTFOLDER} --outdir ${MAINDIR}/${OUTPUTFOLDER} ${OTHEROPTIONS}
    done
}

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
#   Entry Point - Main Function
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Check Input
if [ "$#" -eq 0 ]; then
    echo "No Stage Provided"
    return
fi

echo "Submitting Jobs ..."
echo "... Output Folder: ${OUTPUTFOLDER}"

# Particle Cannon Stage
if [ "${STAGE}" == "cal-pc" ]; then
    run_cal_pc  
fi

# MINOS Stage
if [ "${STAGE}" == "minos" ]; then
    run_minos
fi

# Reconstruction Stage
if [ "${STAGE}" == "reco" ]; then
    run_reco
fi

# DST Stage
if [ "${STAGE}" == "dst" ]; then
    run_dst
fi

echo "All Runs Submitted!"

