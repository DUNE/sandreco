
ROOTGLIBS = $(shell root-config --glibs)
ROOTINCDIR = $(shell root-config --incdir)

EDEPSIM = /wd/sw/EDEPSIM/edep-sim.binary
EDEPGLIBS = -L$(EDEPSIM)/lib/ -ledepsim_io

all: check

check: main.cc
	g++ $^ -pthread -std=c++11 -m64 -I${ROOTINCDIR} ${ROOTGLIBS} -lEG -lGeom -lEve -I../include ${EDEPGLIBS} -L../lib -lStruct -o $@
