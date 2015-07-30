#===============================================================================
# setup.sh
# Creates Output/ Directory with required subdirectories
# Use only once you downloaded the package - cleanup.sh does not removes folders
# 
#    Terminal Usage:
#        > source setup.sh
#    
#===============================================================================

echo "Setting Up Package"

DefaultPlots="Plots"
Main_Folder="Output"

Folder_List=(TextFiles
        ${DefaultPlots})
		
PlotFolder_List=(Muon
		Proton
		Pion
		Interaction
        CutHists
		Other)

# Create Main Folder
echo "... Creating ${Main_Folder}"
mkdir ${Main_Folder}
cd ${Main_Folder}
    
# Create Folders inside Main Folders
for folder in ${Folder_List[@]}; do
    echo "... Creating ${Main_Folder}/${folder}"
    mkdir ${folder}
done

# Create SubFolders for Plots in Each Branch
cd ${DefaultPlots}
for plotFolder in ${PlotFolder_List[@]}; do
    echo "... Creating ${Main_Folder}/${DefaultPlots}/${plotFolder}"
    mkdir ${plotFolder}
done
cd ..
cd ..
echo "Setup Finished!"
