# Deletes Files & Folders from cmt make Command
# Required for my personal gitHub Account
echo "Cleaning Build Files & Folders"
rm -rv ../x86_64-slc5-gcc43-opt
rm -rv ../genConf
rm -v cleanup.*
rm -v setup.*
rm -v *.history
rm -v Makefile
echo "Done!"