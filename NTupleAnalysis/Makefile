# ==============================================================================
#
# Makefile for NTupleAnalysis Package
# 	Statistical Data Analysis Package for CCProtonPi0
#	Author: Ozgur Altinok
#	Last Revision: 2014_11_11
#
# ==============================================================================

# Main File and its executable
MAINCPP = main.cpp
MAINEXE = main.exe

# Class: Plotter
PLOTTERCPP = Classes/Plotter/Plotter.cpp
DEFAULTPLOTSCPP = Classes/Plotter/Default_Plots.cpp
OTHERPLOTSCPP = Classes/Plotter/Other_Plots.cpp
ALLPLOTTERCPP = $(PLOTTERCPP) $(DEFAULTPLOTSCPP) $(OTHERPLOTSCPP)

# Class: CCProtonPi0
CCPROTONPI0CPP = Classes/CCProtonPi0/CCProtonPi0.cpp
INITHISTOGRAMSCPP = Classes/CCProtonPi0/initHistograms.cpp
DEFAULTFUNCTIONSCPP = Classes/CCProtonPi0/DefaultFunctions.cpp
MUONFUNCTIONSCPP = Classes/CCProtonPi0/MuonFunctions.cpp
PROTONFUNCTIONSCPP = Classes/CCProtonPi0/ProtonFunctions.cpp
PIONFUNCTIONSCPP = Classes/CCProtonPi0/PionFunctions.cpp
ALLCCPROTONPI0CPP = $(CCPROTONPI0CPP) $(INITHISTOGRAMSCPP) $(DEFAULTFUNCTIONSCPP) $(MUONFUNCTIONSCPP) $(PROTONFUNCTIONSCPP) $(PIONFUNCTIONSCPP)

# Class: Particle and its Subclasses
PARTICLECPP = Classes/Particle/Particle.cpp 
MUONCPP = Classes/Muon/Muon.cpp
PROTONCPP = Classes/Proton/Proton.cpp
PIONCPP = Classes/Pion/Pion.cpp
ALLPARTICLESCPP = $(PARTICLECPP) $(MUONCPP) $(PROTONCPP) $(PIONCPP)

# Other Helper Classes
BINLISTCPP = Classes/BinList/BinList.cpp
SINGLEBINCPP = Classes/SingleBin/SingleBin.cpp
CUTCPP = Classes/Cut/Cut.cpp
CUTSTATCPP = Classes/Cut/CutStat.cpp
ALLOTHERCPP = $(BINLISTCPP) $(SINGLEBINCPP) $(CUTCPP) $(CUTSTATCPP)

# ALL Classes Required for the Package
ALLCPP = $(MAINCPP) $(ALLPLOTTERCPP) $(ALLCCPROTONPI0CPP) $(ALLPARTICLESCPP) $(ALLOTHERCPP)

# ------------------------------------------------------------------------------
# Commands and Flags
# ------------------------------------------------------------------------------
CC = g++
CFLAGS = -c -Wall
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

