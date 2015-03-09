# ==============================================================================
#
# Makefile for CAMAC_DataReader
#	Author: Ozgur Altinok
#	Last Revision: 2015_03_02
#
# ==============================================================================

# Main File and its executable
MAINCPP = src/main.cpp
MAINEXE = CAMACDataReader

# CAMAC_Data
CAMACDATACPP = src/CAMAC_Data.cpp

# ALL Classes Required for the Package
ALLCPP = $(MAINCPP) $(CAMACDATACPP) 

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

