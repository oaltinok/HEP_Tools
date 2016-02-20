#
# Finds All Ana_Tuple Files for the analysis
# outputs list of files to  pl_DATATYPE_Merged.dat
#
DATATYPE=$1
CCPROTONPI0_V="v2_60"

input_error()
{
    echo "No or Wrong Data Type Specified!"
    echo "Correct Data Types"
    echo "    mc"
    echo "    data"
}

find_files()
{
    find /pnfs/minerva/persistent/users/oaltinok/CCProtonPi0/${DATAFOLDER}/${CCPROTONPI0_V}/*/${ANADIR} -name "*Ana*Tuple*.root" -printf "%p\n" | sort > Playlists/pl_${DATAFOLDER}_Merged.dat
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






