#===============================================================================
# setup.sh
# Creates Output/ Directory with required subdirectories
# Use only once you downloaded the package - cleanup.sh does not removes folders
# 
#    Terminal Usage:
#        > source setup.sh
#    
#    Last Revision: 2015_02_18
#===============================================================================

echo "Setting Up Package"

DefaultPlots="Plots"
Main_Folder="Output"
PlotOther_Folder="Other"

Folder_List=(RootFiles
        TextFiles
        ${DefaultPlots})

Branch_List=(Signal
		Background
		AllEvents)
		
PlotFolder_List=(Muon
		Proton
		Pion
		Interaction
		PIDStatistics
		Other)

# Create Main Folder
echo "... Creating ${Main_Folder}"
mkdir ${Main_Folder}
cd ${Main_Folder}
    
# Create Folders inside Main Folders
for folder in ${Folder_List[@]}; do
    echo "... Creating ${Main_Folder}/${folder}"
    mkdir ${folder}
    cd ${folder}
    
    # Create Branches inside Folders
    for fbranch in ${Branch_List[@]}; do
        echo "... Creating ${Main_Folder}/${folder}/${fbranch}"
        mkdir ${fbranch} 
    done
    cd ..
done

# Create SubFolders for Plots in Each Branch
cd ${DefaultPlots}
for fbranch in ${Branch_List[@]}; do
    cd ${fbranch}
    for plotFolder in ${PlotFolder_List[@]}; do
        echo "... Creating ${Main_Folder}/${DefaultPlots}/${fbranch}/${plotFolder}"
        mkdir ${plotFolder}
    done
    cd ..
done

# Create Other Folder inside Plots
echo "... Creating ${Main_Folder}/${DefaultPlots}/${PlotOther_Folder}"
mkdir ${PlotOther_Folder}

cd ..
cd ..

echo "Setup Finished!"
