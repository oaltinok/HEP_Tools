#!/bin/sh

#-------------------------------------------------------------------------------
# Tag & Setup Script for Nightly (SLF6)
#-------------------------------------------------------------------------------
tag="nightly"  #you will work in cmtuser/Minerva_<tag>
setupscript="/minerva/data/software_releases/nightly/SLF6/latest/setup.sh"  
debug=""            #set to -d if you want the debug build

# Setup User Areas
source ~/setupScripts/setup_user_areas.sh

# Setup Minerva Software
source ~/setupScripts/setup_minerva_soft.sh

# Setup Project.cmt -- Only for Nighlty
source ~/setupScripts/setup_project_cmt.sh

# setup packages
source ~/setupScripts/setup_packages.sh

# setup user shortcuts
source ~/setupScripts/setup_user_shortcuts.sh

echo "All Done!"

