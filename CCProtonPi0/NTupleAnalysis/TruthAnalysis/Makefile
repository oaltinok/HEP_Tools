# ==============================================================================
#
# Makefile for Truth Analysis
#	Author: Ozgur Altinok
#
# ==============================================================================

# Main File and its executable
MAINCPP = main.cpp
MAINEXE = main.exe

# Classes
TRUTHANALYZERCPP = Classes/TruthAnalyzer/CCProtonPi0_TruthAnalyzer.cpp
TRUTHANALYZERDEFAULTCPP = Classes/TruthAnalyzer/CCProtonPi0_TruthAnalyzer_Default.cpp
TRUTHANALYZERALLCPP = $(TRUTHANALYZERCPP) $(TRUTHANALYZERDEFAULTCPP)

# ALL Classes Required for the Package
ALLCPP = $(MAINCPP) $(TRUTHANALYZERALLCPP)

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
	$(CC) $(ROOTFLAGS) $(CFLAGS) -o $(EXECUTABLE) $(OBJECTS) 

.cpp.o:
	$(CC) $(ROOTFLAGS) $(CFLAGS)  -c $< -o $@
	
clean: 
	$(RM) $(MAINEXE) $(OBJECTS)

