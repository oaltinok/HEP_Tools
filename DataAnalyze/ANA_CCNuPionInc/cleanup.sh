#===============================================================================
# cleanup.sh
# Removes all files inside Output/ Directory
# CAUTION: Removes ALL Files, if you want to keep some files, edit the function
#               or change directory configuration
# 
#    Terminal Usage:
#        > source cleanup.sh
#    
#    Last Revision: 2014_01_27
#===============================================================================

echo "Cleanup Process for Output Directories Started!"
echo "Cleaning Output Folder"

Folder_List=(Output/Plots/*
                Output/RootFiles/*
                Output/TextFiles/*)

for ilist1 in ${Folder_List[@]}; do
    echo "--"
    echo "Clean: ${ilist1}"
    rm -frv ${ilist1}
done

echo "Cleanup Finished!"

