CCPROTONPI0_V="Nieves_2p2h"
OUTDIR="/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/MC/${CCPROTONPI0_V}" 
OTHEROPTIONS=""
RUNS=( "13201" "13202" "117200" "117201" "117202" "117203" )
#RUNS=( "13200" )

submit_jobs()
{
    echo "Submitting Jobs ..."
    echo "... CCProtonPi0 Version: ${CCPROTONPI0_V}" 

    for run in ${RUNS[@]}; do
        echo "Submitting $run"
        python $ana_script/ProcessAna.py --mc --run $run --ana_tool CCProtonPi0 --inv v10r8p6 --kludge Eroica --no_verify_kludge --indir /pnfs/minerva/persistent/users/rodriges/mecgenie-deltatag/ --outdir ${OUTDIR} ${OTHEROPTIONS}
    done

    echo "All Playlists Submitted!"
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
#   Entry Point - Main Function
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

submit_jobs

