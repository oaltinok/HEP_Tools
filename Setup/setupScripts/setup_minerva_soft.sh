#!/bin/sh

#-------------------------------------------------------------------------------
# Setup for MINERvA Framework
#-------------------------------------------------------------------------------
# Setup Minerva Script
echo source $setupscript $debug
source $setupscript $debug

# Setup Minerva Environment
echo "Setting Minerva Environment"
setenvMinerva ${tag}


