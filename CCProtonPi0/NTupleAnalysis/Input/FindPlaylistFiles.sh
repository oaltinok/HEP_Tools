DATATYPE=$1
CCPROTONPI0_V="v2_60"
FINDSCRIPT="/minerva/app/users/oaltinok/cmtuser/Minerva_v10r8p7/Tools/ProductionScripts/pot_scripts/FindAnaFilesWithoutMeta.py"
#PLAYLISTS=( "minerva1" "minerva7" "minerva9" "minerva13A" "minerva13B" "minerva13C" "minerva13D" "minerva13E" )
PLAYLISTS=( "minerva1" );

input_error()
{
    echo "No or Wrong Data Type Specified!"
    echo "Correct Data Types"
    echo "    mc"
    echo "    data"
}

find_files()
{
    for pl in ${PLAYLISTS[@]}; do
        echo "Finding ${DATAFOLDER} Files for ${pl} in version ${CCPROTONPI0_V}"
        DATADIR="/pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/${DATAFOLDER}/${CCPROTONPI0_V}/${pl}/${ANADIR}"
        OUTPUTFILE="/minerva/data/users/oaltinok/NTupleAnalysis/${DATAFOLDER}/Merged/merge_files_${DATATYPE}_${pl}_${CCPROTONPI0_V}.sh" 
        python ${FINDSCRIPT} --dir ${DATADIR} --batch --ana_tool CCProtonPi0 --output ${OUTPUTFILE} --remove
    done
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

if [ \( "${DATATYPE}" != "mc" \) -a \( "${DATATYPE}" != "data" \) ]; then
    input_error
    return
fi

# Set Data Folder
if [ "${DATATYPE}" == "mc" ]; then
    DATAFOLDER="MC"
    ANADIR="grid/central_value/minerva/ana"
else
    DATAFOLDER="Data"
    ANADIR="grid/minerva/ana"
fi

find_files



