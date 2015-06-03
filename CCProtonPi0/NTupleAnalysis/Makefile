# ==============================================================================
#
# Makefile for NTupleAnalysis Package
# 	Statistical Data Analysis Package for CCProtonPi0
#	Author: Ozgur Altinok
#
# ==============================================================================

# Main File and its executable
MAINCPP = main.cpp
MAINEXE = main.exe

# Class: Base Class
NTUPLEANALYSISCPP = Classes/NTupleAnalysis/CCProtonPi0_NTupleAnalysis.cpp

# Class: Plotter
PLOTTERCPP = Classes/Plotter/CCProtonPi0_Plotter.cpp

# Class: Analyzer
ANALYZERCPP = Classes/Analyzer/CCProtonPi0_Analyzer.cpp

# Class: Interaction, Particle and Particle SubClasses
INTERACTIONCPP = Classes/Interaction/CCProtonPi0_Interaction.cpp
PARTICLECPP = Classes/Particle/CCProtonPi0_Particle.cpp 
MUONCPP = Classes/Muon/CCProtonPi0_Muon.cpp
PROTONCPP = Classes/Proton/CCProtonPi0_Proton.cpp
PIONCPP = Classes/Pion/CCProtonPi0_Pion.cpp
PI0BLOBCPP = Classes/Pi0Blob/CCProtonPi0_Pi0Blob.cpp
ALLPARTICLESCPP = $(INTERACTIONCPP) $(PARTICLECPP) $(MUONCPP) $(PROTONCPP) $(PIONCPP) $(PI0BLOBCPP)

# Tools for Detailed Studies
PIDTOOLCPP = Classes/PIDTool/CCProtonPi0_PIDTool.cpp
BACKGROUNDTOOLCPP = Classes/BackgroundTool/CCProtonPi0_BackgroundTool.cpp
MICHELTOOLCPP = Classes/MichelTool/CCProtonPi0_MichelTool.cpp
ALLTOOLCPP = $(PIDTOOLCPP) $(BACKGROUNDTOOLCPP) $(MICHELTOOLCPP)

# Other Helper Classes
BINLISTCPP = Classes/BinList/CCProtonPi0_BinList.cpp
CUTLISTCPP = Classes/CutList/CCProtonPi0_CutList.cpp
SINGLEBINCPP = Classes/SingleBin/CCProtonPi0_SingleBin.cpp
CUTCPP = Classes/Cut/CCProtonPi0_Cut.cpp
CUTSTATCPP = Classes/Cut/CCProtonPi0_CutStat.cpp
ALLOTHERCPP = $(BINLISTCPP) $(SINGLEBINCPP) $(CUTLISTCPP) $(CUTCPP) $(CUTSTATCPP)

# ALL Classes Required for the Package
ALLCPP = $(MAINCPP) $(NTUPLEANALYSISCPP) $(ANALYZERCPP) $(ALLPARTICLESCPP) $(ALLTOOLCPP) $(ALLOTHERCPP) $(PLOTTERCPP)

# ------------------------------------------------------------------------------
# Commands and Flags
# ------------------------------------------------------------------------------
CC = g++
CFLAGS = -c -Wall -Wextra
ROOTFLAGS = `root-config --libs --cflags`
RM = rm -f

SOURCES = $(ALLCPP)
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = $(MAINEXE)

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(ROOTFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(ROOTFLAGS) $(CFLAGS) $< -o $@
	
clean: 
	$(RM) $(MAINEXE) $(OBJECTS)

