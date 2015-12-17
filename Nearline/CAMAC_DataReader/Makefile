# ==============================================================================
#
# Makefile for CAMAC_DataReader
#	Author: Ozgur Altinok
#	Last Revision: 2015_03_09
#
# ==============================================================================

# Main File and its executable
MAINCPP = main.cpp
CAMACDATAREADERCPP = Classes/CAMAC_DataReader/CAMAC_DataReader.cpp
MAINEXE = CAMAC_DataReader

# Class: Plotter and its subclasses
PLOTTERCPP = Classes/Plotter/Plotter.cpp
HISTCPP = Classes/Hist/Hist.cpp
GRAPHCPP = Classes/Graph/Graph.cpp
FREQCPP = Classes/Freq/Freq.cpp
ALLPLOTTERCPP = $(PLOTTERCPP) $(HISTCPP) $(GRAPHCPP) $(FREQCPP)

# Class: Reader and its subclasses
READERCPP = Classes/Reader/Reader.cpp 
DATACPP = Classes/Data/Data.cpp
CONFIGCPP = Classes/Config/Config.cpp
ALLREADERCPP = $(READERCPP) $(DATACPP) $(CONFIGCPP)


# ALL Classes Required for the Package
ALLCPP = $(MAINCPP) $(CAMACDATAREADERCPP) $(ALLPLOTTERCPP) $(ALLREADERCPP)

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

