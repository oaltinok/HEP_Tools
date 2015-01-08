#===============================================================================
# cleanup.sh
# Removes Output/ Directory
# CAUTION: Removes ALL Files, if you want to keep some files, edit the function
#               or change directory configuration
# 
#    Terminal Usage:
#        > source cleanup.sh
#    
#    Last Revision: 2014_11_06
#===============================================================================

echo "Removing Output Folder"

rm -frv Output

echo "Cleanup Finished!"

