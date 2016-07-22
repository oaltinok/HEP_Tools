#! /bin/sh

#-------------------------------------------------------------------------------
# Setup Installed packages
#-------------------------------------------------------------------------------
Package_List=(Tools/SystemTests
    Tools/CondorUtils
    Tools/ProductionScripts
    Rec/ParticleMaker
    Rec/ProngMaker
    Ana/CCProtonPi0
    Ana/PlotUtils
    Ana/UnfoldUtils
    MParamFiles)

Default_Setup=/cmt/setup.sh

echo "Setting Installed Packages"
for ilist1 in ${Package_List[@]}; do
    echo "--"
    echo "Setup: ${ilist1}"
    source ${ilist1}$Default_Setup
done


