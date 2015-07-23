#!/bin/sh

#-------------------------------------------------------------------------------
# Setup User Release Areas
#-------------------------------------------------------------------------------
echo "Setting User Release Area = User_release_area"
export User_release_area="/minerva/app/users/$USER/cmtuser"
echo $User_release_area

echo "Setting User Data Area = USER_DATA_AREA"
export USER_DATA_AREA=/minerva/data/users/oaltinok
echo $USER_DATA_AREA


