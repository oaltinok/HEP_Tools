MAINEXE  = main.exe
MAINCPP = main.cpp
MAINO	  = main.o

BINLISTCPP = Classes/BinList/BinList.cpp
BINLISTO = Classes/BinList/BinList.o

SINGLEBINCPP = Classes/SingleBin/SingleBin.cpp
SINGLEBINO = Classes/SingleBin/SingleBin.o

CCDELTAPLUSCPP = Classes/CCDeltaPlus/CCDeltaPlus.cpp
CCDELTAPLUSO = Classes/CCDeltaPlus/CCDeltaPlus.o

PARTICLECPP = Classes/Particle/Particle.cpp 
PARTICLEO = Classes/Particle/Particle.o

MUONCPP = Classes/Muon/Muon.cpp
MUONO = Classes/Muon/Muon.o

PROTONCPP = Classes/Proton/Proton.cpp
PROTONO = Classes/Proton/Proton.o

PIONCPP = Classes/Pion/Pion.cpp
PIONO = Classes/Pion/Pion.o

PLOTTERCPP = Classes/Plotter/Plotter.cpp
PLOTTERO = Classes/Plotter/Plotter.o 


#All header files and source files are expected in the current directory

#The first target defined is the default for the makefile
CXX	= g++
RM	= rm -f
SHELL = sh

##### Library path and compiler / linker / rootcint options #####
CXXFLAGS	= -o
CXX		= g++

#FLAGS for the root libraries and include paths
ROOTFLAGS   = `root-config --libs --cflags`

ALLCPP = $(MAINCPP) $(BINLISTCPP) $(SINGLEBINCPP) $(CCDELTAPLUSCPP) $(PARTICLECPP) $(MUONCPP) $(PROTONCPP) $(PIONCPP) $(PLOTTERCPP)

ALLO = $(MAINO) $(BINLISTO) $(SINGLEBINO) $(CCDELTAPLUSO) $(PARTICLEO) $(MUONO) $(PROTONO) $(PIONCO) $(PLOTTERO)


all:
	$(CXX) $(ROOTFLAGS) $(ALLCPP) $(CXXFLAGS) $(MAINEXE)
	
clean: 
	$(RM) $(MAINEXE)

