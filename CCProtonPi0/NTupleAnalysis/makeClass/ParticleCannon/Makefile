# ==============================================================================
#
# Makefile for Particle Cannon DST Output Analysis Package
#	Author: Ozgur Altinok
#
# ==============================================================================

# Main File and its executable
MAINCPP = main.cpp
MAINEXE = main.exe

# Class: Base Class
PC_DST_CPP = PC_DST.cpp

# ALL Classes Required for the Package
ALLCPP = $(MAINCPP) $(PC_DST_CPP)

# ------------------------------------------------------------------------------
# Commands and Flags
# ------------------------------------------------------------------------------
CC = g++
CFLAGS = -g -Wall -Wextra
ROOTFLAGS = `root-config --glibs --cflags`
RM = rm -fv

SOURCES = $(ALLCPP)
OBJECTS = $(SOURCES:.cpp=.o)
EXECUTABLE = $(MAINEXE)

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(ROOTFLAGS) $(CFLAGS) -o $(EXECUTABLE) $(OBJECTS)

.cpp.o:
	$(CC) $(ROOTFLAGS) $(CFLAGS) -c $< -o $@
	
clean: 
	$(RM) $(MAINEXE) $(OBJECTS)

