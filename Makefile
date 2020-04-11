
HEADERS= struct.h
CFLAGS = $(shell root-config --cflags)
LDFLAGS = $(shell root-config --ldflags)
ROOTGLIBS = $(shell root-config --glibs)
ROOTINCDIR = $(shell root-config --incdir)

EDEPSIM = /wd/sw/EDEPSIM/edep-sim.binary

EDEPGLIBS = -L$(EDEPSIM)/lib/ -ledepsim_io
EDEPINCDIR = $(EDEPSIM)/include/EDepSim

GENFIT=/wd/sw/GENFIT/GenFit.binary

EIGEN3=/usr/include/eigen3

all: Digitize Reconstruct Analyze FastCheck

struct.cxx: include/struct.h include/Linkdef.h
	cd include && rootcint -f ../src/$@ -c $(CFLAGS) -p $(HEADERS) Linkdef.h && cd ..

libStruct.so: struct.cxx
	g++ -shared -fPIC -o lib/$@ -Iinclude `root-config --ldflags` $(CFLAGS) src/$^ && cp src/struct_rdict.pcm lib/struct_rdict.pcm

Digitize: libStruct.so
	g++ -o bin/$@ $(CFLAGS) $(LDFLAGS) -I$(EDEPINCDIR) -Iinclude $(ROOTGLIBS) -lGeom \
	$(EDEPGLIBS) -Llib -lStruct src/digitization.cpp

Reconstruct: libStruct.so
	g++ -o bin/$@ $(CFLAGS) $(LDFLAGS) -I$(EDEPINCDIR) -Iinclude $(ROOTGLIBS) -lGeom \
	$(EDEPGLIBS) -Llib -lStruct src/reconstruction.cpp

Analyze: libStruct.so
	g++ -o bin/$@ $(CFLAGS) $(LDFLAGS) -I$(EDEPINCDIR) -Iinclude -I${GENFIT}/include -I${EIGEN3} -L${GENFIT}/lib64 -lgenfit2 $(ROOTGLIBS) -lGeom -lEG \
	$(EDEPGLIBS) -Llib -lStruct src/analysis.cpp

FastCheck: libStruct.so
	g++ -o bin/$@ $(CFLAGS) $(LDFLAGS) -I$(EDEPINCDIR) -Iinclude -I${GENFIT}/include -I${EIGEN3} -L${GENFIT}/lib64 -lgenfit2 $(ROOTGLIBS) -lGeom -lEG \
	$(EDEPGLIBS) -Llib -lStruct src/fastcheck.cpp
