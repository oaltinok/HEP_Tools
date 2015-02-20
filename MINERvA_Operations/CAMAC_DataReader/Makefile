# ==============================================================================
#
# Makefile for Test Beam TOF and Veto Scalars
#	Author: Ozgur Altinok
#	Last Revision: 2015_02_16
#
# ==============================================================================

# Main File and its executable
MAINCPP = src/main.cpp
MAINEXE = CAMACDataReader

# CAMAC_Data
CAMACDATACPP = src/CAMAC_Data.cpp

# CAMAC_DataReader
CAMACDATAREADERCPP = src/CAMAC_DataReader.cpp


# ALL Classes Required for the Package
ALLCPP = $(MAINCPP) $(CAMACDATACPP) $(CAMACDATAREADERCPP)

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

