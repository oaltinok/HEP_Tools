# ==============================================================================
#
# Makefile for Truth Analysis
#	Author: Ozgur Altinok
#
# ==============================================================================

# Main File and its executable
MAINCPP = main.cpp
MAINEXE = main.exe

# Class: Base Class
NTUPLEANALYSISCPP = ../Classes/NTupleAnalysis/CCProtonPi0_NTupleAnalysis.cpp

# Classes
TRUTHANALYZERCPP = Classes/TruthAnalyzer/CCProtonPi0_TruthAnalyzer.cpp
TRUTHANALYZERDEFAULTCPP = Classes/TruthAnalyzer/CCProtonPi0_TruthAnalyzer_Default.cpp
TRUTHANALYZERALLCPP = $(TRUTHANALYZERCPP) $(TRUTHANALYZERDEFAULTCPP)

# ALL Classes Required for the Package
ALLCPP = $(MAINCPP) $(NTUPLEANALYSISCPP) $(TRUTHANALYZERALLCPP)

# PlotUtils from MINERvA Framework
INCLUDES += -I$(PLOTUTILSROOT)/ -I$(UNFOLDUTILSROOT)/
LDLIBS += -L$(PLOTUTILSROOT)/$(CMTCONFIG) -lplotutils -L$(UNFOLDUTILSROOT)/$(CMTCONFIG) -lUnfoldUtils

# ------------------------------------------------------------------------------
# Commands and Flags
# ------------------------------------------------------------------------------
CC = g++
CFLAGS = -g -Wall -Wextra -std=c++0x
ROOTFLAGS = `root-config --glibs --cflags`
RM = rm -fv

SOURCES = $(ALLCPP)
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = $(MAINEXE)

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(ROOTFLAGS) $(CFLAGS) $(INCLUDES) -o $(EXECUTABLE) $(OBJECTS) $(LDLIBS) 

.cpp.o:
	$(CC) $(ROOTFLAGS) $(CFLAGS) $(INCLUDES) -c $< -o $@
	
clean: 
	$(RM) $(MAINEXE) $(OBJECTS)

