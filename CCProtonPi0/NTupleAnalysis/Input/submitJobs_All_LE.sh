DATATYPE=$1
CCPROTONPI0_V="v2_93"
OTHEROPTIONS=""
PLAYLISTS=( "minerva1" "minerva7" "minerva9" "minerva13A" "minerva13B" "minerva13C" "minerva13D" "minerva13E" )
#PLAYLISTS=( "minerva1" )

input_error()
{
    echo "No or Wrong Data Type Specified!"
    echo "Correct Data Types"
    echo "    mc"
    echo "    data"
}

submit_jobs()
{
    echo "Submitting Jobs ..."
    echo "... CCProtonPi0 Version: ${CCPROTONPI0_V}" 

    for pl in ${PLAYLISTS[@]}; do
        echo "Submitting $pl"
        if [ "$DATATYPE" == "mc" ]; then
            OUTDIR="/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/${DATAFOLDER}/${CCPROTONPI0_V}/${pl}" 
            python $ana_script/ProcessAna.py --${DATATYPE} --playlist ${pl} --usecat --ana_tool CCProtonPi0 --inv v10r8p9 --kludge Eroica --no_verify_kludge --outdir ${OUTDIR} ${OTHEROPTIONS}
        else
            OUTDIR="/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/${DATAFOLDER}/${CCPROTONPI0_V}/${pl}" 
            python $ana_script/ProcessAna.py --${DATATYPE} --playlist ${pl} --usecat --ana_tool CCProtonPi0 --inv eroica --outdir ${OUTDIR} ${OTHEROPTIONS}
        fi
    done

    echo "All Playlists Submitted!"
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
#   Entry Point - Main Function
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Check Input
if [ "$#" -eq 0 ]; then
    input_error
    return
fi

if [ \( "$1" != "mc" \) -a \( "$1" != "data" \) ]; then
    input_error
    return
fi

# Set Data Folder
if [ "$1" == "mc" ]; then
    DATAFOLDER="MC"
else
    DATAFOLDER="Data"
fi

submit_jobs

