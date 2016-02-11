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

DefaultText="TextFiles"
DefaultPlots="Plots"
Main_Folder="Output"

Folder_List=(${DefaultText}
        ${DefaultPlots})

PlotFolder_List=(Muon
		Proton
		Pion
		Pi0Blob
		Interaction
        CutHists
        CrossSection
        Efficiency 
        BackgroundSubtraction 
        Unfolding 
        Errors 
        OtherStudies)

# Create Main Folder
echo "... Creating ${Main_Folder}"
mkdir ${Main_Folder}
cd ${Main_Folder}
    
# Create Folders inside Main Folders
for folder in ${Folder_List[@]}; do
    echo "... Creating ${Main_Folder}/${folder}"
    mkdir ${folder}
done

# Create Folders inside Plots Folders
cd ${DefaultPlots}
for folder in ${PlotFolder_List[@]}; do
    echo "... Creating ${Main_Folder}/${DefaultPlots}/${folder}"
    mkdir ${folder}
done
cd .. # exit DefaultPlots
cd .. # exit Main_Folder

echo "Setup Finished!"

