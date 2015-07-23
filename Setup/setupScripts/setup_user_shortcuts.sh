#! /bin/sh

#-------------------------------------------------------------------------------
# Folder Location List
#-------------------------------------------------------------------------------
export data=${USER_DATA_AREA}
export app=${User_release_area}Minerva_$tag/Ana/CCProtonPi0            
export ana_script=${User_release_area}Minerva_$tag/Tools/ProductionScripts/ana_scripts
export options=${User_release_area}Minerva_$tag/Tools/SystemTests/options/Analysis

echo ""
echo "-------------------------------------------------------------------------------"
echo "Shortcuts to Folders: "
echo "data = " ${data}
echo "app = " ${app}
echo "ana_script = " ${ana_script}
echo "-------------------------------------------------------------------------------"


#------------------------------------------------------------------------------
# Alias List
#-------------------------------------------------------------------------------

alias job_q='jobsub_q --user oaltinok'

echo " "
echo "-------------------------------------------------------------------------------"
echo "Alias List"
echo "job_q = jobsub_q --user oaltinok"
echo "-------------------------------------------------------------------------------"


