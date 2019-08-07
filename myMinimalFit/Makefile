
GENFIT=/wd/sw/GENFIT/GenFit.binary
ROOTGLIBS = $(shell root-config --glibs)
ROOTINCDIR = $(shell root-config --incdir)
EIGEN3=/usr/include/eigen3

all: myMinimalFit mySimpleFit

myMinimalFit: main.cc
	g++ $^ -pthread -std=c++11 -m64 -L${GENFIT}/lib64 -lgenfit2 -I${ROOTINCDIR} ${ROOTGLIBS} -lEG -lGeom -lEve -I${EIGEN3} -I${GENFIT}/include -o $@

mySimpleFit: simple.cc
	g++ $^ -pthread -std=c++11 -m64 -L${GENFIT}/lib64 -lgenfit2 -I${ROOTINCDIR} ${ROOTGLIBS} -lEG -lGeom -lEve -I${EIGEN3} -I${GENFIT}/include -o $@
