
GENFIT=/wd/sw/GENFIT/GenFit.binary
ROOTGLIBS = $(shell root-config --glibs)
ROOTINCDIR = $(shell root-config --incdir)
EIGEN3=/usr/include/eigen3

EDEPSIM = /wd/sw/EDEPSIM/edep-sim.binary
EDEPGLIBS = -L$(EDEPSIM)/lib/ -ledepsim -ledepsim_io

all: myMinimalFit mySimpleFit

myMinimalFit: main.cc
	g++ $^ -pthread -std=c++11 -m64 -L${GENFIT}/lib64 -lgenfit2 -I${ROOTINCDIR} ${ROOTGLIBS} -lEG -lGeom -lEve -I${EIGEN3} -I${GENFIT}/include -I../struct ${EDEPGLIBS} -o $@

mySimpleFit: simple.cc
	g++ $^ -pthread -std=c++11 -m64 -L${GENFIT}/lib64 -lgenfit2 -I${ROOTINCDIR} ${ROOTGLIBS} -lEG -lGeom -lEve -I${EIGEN3} -I${GENFIT}/include -I../struct ${EDEPGLIBS} -o $@
