#===============================================================================
# setup.sh
# Creates Output/ Directory with required subdirectories
# Use only once you downloaded the package - cleanup.sh does not removes folders
# 
#    Terminal Usage:
#        > source setup.sh
#    
#    Last Revision: 2014_02_03
#===============================================================================

echo "Setting Up Package"

Folder_List=(Output
                Output/Plots
                Output/RootFiles
                Output/TextFiles)

for ilist1 in ${Folder_List[@]}; do
    echo "Creating ${ilist1}"
    mkdir ${ilist1}
done

echo "Setup Finished!"