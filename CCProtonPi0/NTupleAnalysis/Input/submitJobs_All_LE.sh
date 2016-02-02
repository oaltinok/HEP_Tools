DATATYPE=$1
CCPROTONPI0_V="v2_54"
OTHEROPTIONS=""
PLAYLISTS=( "minerva1" "minerva7" "minerva9" "minerva13A" "minerva13B" "minerva13C" "minerva13D" "minerva13E" )

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
        python $ana_script/ProcessAna.py --${DATATYPE} --playlist ${pl} --usecat --ana_tool CCProtonPi0 --inv eroica --use_jobsub_client --dcachein --dcacheout --outdir /minerva/data/users/oaltinok/CCProtonPi0/${DATAFOLDER}/${CCPROTONPI0_V}/${pl} ${OTHEROPTIONS}
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

