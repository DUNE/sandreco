CFLAGS = $(shell root-config --cflags)
LDFLAGS = $(shell root-config --ldflags)
ROOTGLIBS = $(shell root-config --glibs)
ROOTINCDIR = $(shell root-config --incdir)

EDEPSIM = /wd/sw/EDEPSIM/edep-sim.binary

EDEPGLIBS = -L$(EDEPSIM)/lib/ -ledepsim -ledepsim_io
EDEPINCDIR = $(EDEPSIM)/include/EDepSim


all: Digitize

Digitize: digitization.cpp
	g++ -o$@ $(CFLAGS) $(LDFLAGS) -I$(EDEPINCDIR) -I../struct $(ROOTGLIBS) -lGeom $(EDEPGLIBS) -L../struct -lStruct $^
