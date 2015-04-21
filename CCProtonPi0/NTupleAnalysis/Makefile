# ==============================================================================
#
# Makefile for NTupleAnalysis Package
# 	Statistical Data Analysis Package for CCProtonPi0
#	Author: Ozgur Altinok
#	Last Revision: 2015_04_21
#
# ==============================================================================

# Main File and its executable
MAINCPP = main.cpp
MAINEXE = main.exe

# Class: Base Class
NTUPLEANALYSISCPP = Classes/NTupleAnalysis/NTupleAnalysis.cpp

# Class: Plotter
PLOTTERCPP = Classes/Plotter/Plotter.cpp
DEFAULTPLOTSCPP = Classes/Plotter/Default_Plots.cpp
OTHERPLOTSCPP = Classes/Plotter/Other_Plots.cpp
PIDPLOTSCPP = Classes/Plotter/pID_Plots.cpp
ALLPLOTTERCPP = $(PLOTTERCPP) $(DEFAULTPLOTSCPP) $(OTHERPLOTSCPP) $(PIDPLOTSCPP)

# Class: Analyzer
ANALYZERCPP = Classes/Analyzer/Analyzer.cpp
INITFUNCTIONSCPP = Classes/Analyzer/initFunctions.cpp
DEFAULTFUNCTIONSCPP = Classes/Analyzer/DefaultFunctions.cpp
MUONFUNCTIONSCPP = Classes/Analyzer/MuonFunctions.cpp
PROTONFUNCTIONSCPP = Classes/Analyzer/ProtonFunctions.cpp
PIONFUNCTIONSCPP = Classes/Analyzer/PionFunctions.cpp
ALLANALYZERCPP = $(ANALYZERCPP) $(INITFUNCTIONSCPP) $(DEFAULTFUNCTIONSCPP) $(MUONFUNCTIONSCPP) $(PROTONFUNCTIONSCPP) $(PIONFUNCTIONSCPP)

# Class: Particle and its Subclasses
PARTICLECPP = Classes/Particle/Particle.cpp 
MUONCPP = Classes/Muon/Muon.cpp
PROTONCPP = Classes/Proton/Proton.cpp
PIONCPP = Classes/Pion/Pion.cpp
ALLPARTICLESCPP = $(PARTICLECPP) $(MUONCPP) $(PROTONCPP) $(PIONCPP)

# Other Helper Classes
PIDTOOLCPP = Classes/PIDTool/PIDTool.cpp
BACKGROUNDTOOLCPP = Classes/BackgroundTool/BackgroundTool.cpp
BINLISTCPP = Classes/BinList/BinList.cpp
CUTLISTCPP = Classes/CutList/CutList.cpp
SINGLEBINCPP = Classes/SingleBin/SingleBin.cpp
CUTCPP = Classes/Cut/Cut.cpp
CUTSTATCPP = Classes/Cut/CutStat.cpp
ALLOTHERCPP = $(PIDTOOLCPP) $(BACKGROUNDTOOLCPP) $(BINLISTCPP) $(SINGLEBINCPP) $(CUTLISTCPP) $(CUTCPP) $(CUTSTATCPP)

# ALL Classes Required for the Package
ALLCPP = $(MAINCPP) $(NTUPLEANALYSISCPP) $(ALLPLOTTERCPP) $(ALLANALYZERCPP) $(ALLPARTICLESCPP) $(ALLOTHERCPP)

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

