#! /bin/sh

#-------------------------------------------------------------------------------
# Tag & Setup Script for a Release (SLF6)
#-------------------------------------------------------------------------------
tag="v10r8p7"
setupscript="/grid/fermiapp/minerva/software_releases/${tag}/setup.sh"  
debug=""            #set to -d if you want the debug build

# Setup User Areas
source ~/setupScripts/setup_user_areas.sh

# Setup Minerva Software
source ~/setupScripts/setup_minerva_soft.sh

# Setup Packages
source ~/setupScripts/setup_packages.sh

# Setup User Shortcuts
source ~/setupScripts/setup_user_shortcuts.sh

echo "All Done!"

