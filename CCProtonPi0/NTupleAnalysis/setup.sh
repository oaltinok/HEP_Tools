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

DefaultPlotFolder_List=(Muon
		Proton
		Pion
		Pi0Blob
		Interaction
        CutHists
        OtherStudies)

VariableList=(xsec_pi0_P
        xsec_muon_P
        xsec_muon_theta)

XSecPlotList=(Original
        BackgroundEstimated
        BackgroundSubtracted
        Unfolded
        EfficiencyCorrected
        CrossSection
        ErrorSummary 
        Check)

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
for folder in ${DefaultPlotFolder_List[@]}; do
    echo "... Creating ${Main_Folder}/${DefaultPlots}/${folder}"
    mkdir ${folder}
done
# Create Variable Folder inside Plots/
for folder in ${VariableList[@]}; do
    echo "... Creating ${Main_Folder}/${DefaultPlots}/${folder}"
    mkdir ${folder}
    cd ${folder}
    # Create XSec Folders inside Plots/Variable
    for subfolder in ${XSecPlotList[@]}; do
        echo "... Creating ${Main_Folder}/${DefaultPlots}/${folder}/${subfolder}"
        mkdir ${subfolder}
    done
    cd .. ## exit Variable
done
cd .. # exit DefaultPlots
cd .. # exit Main_Folder

echo "Setup Finished!"

