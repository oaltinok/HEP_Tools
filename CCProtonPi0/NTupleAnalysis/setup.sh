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
DefaultCutHists="CutHists"
Main_Folder="Output"

Folder_List=(${DefaultText}
        ${DefaultPlots})

AnaFolder_List=(1Track
    2Track
    All)

PlotFolder_List=(Muon
		Proton
		Pion
		Interaction)

# Create Main Folder
echo "... Creating ${Main_Folder}"
mkdir ${Main_Folder}
cd ${Main_Folder}
    
# Create Folders inside Main Folders
for folder in ${Folder_List[@]}; do
    echo "... Creating ${Main_Folder}/${folder}"
    mkdir ${folder}
done

# Create Analysis Folders inside DefaultText Folder
cd ${DefaultText}
for anaFolder in ${AnaFolder_List[@]}; do
    echo "... Creating ${Main_Folder}/${DefaultText}/${anaFolder}"
    mkdir ${anaFolder}
done
cd .. # Exit DefaultText

# Create CutHists Folder inside DefaultPlots Folder
cd ${DefaultPlots}
echo "... Creating ${Main_Folder}/${DefaultPlots}/${DefaultCutHists}"
mkdir ${DefaultCutHists}

# Create Analysis Folders inside DefaultPlots Folder
for anaFolder in ${AnaFolder_List[@]}; do
    echo "... Creating ${Main_Folder}/${DefaultPlots}/${anaFolder}"
    mkdir ${anaFolder}
    # Create SubFolders for Plots in Each AnaFolder
    cd ${anaFolder}
    for plotFolder in ${PlotFolder_List[@]}; do
        echo "... Creating ${Main_Folder}/${DefaultPlots}/${anaFolder}/${plotFolder}"
        mkdir ${plotFolder}
    done
    cd .. # exit anaFolder 
done
cd .. # exit DefaultPlots
cd .. # exit Main_Folder

echo "Setup Finished!"

