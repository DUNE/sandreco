
HEADERS= struct.h
CFLAGS = $(shell root-config --cflags)
LDFLAGS = $(shell root-config --ldflags)
ROOTGLIBS = $(shell root-config --glibs)
ROOTINCDIR = $(shell root-config --incdir)

EDEPSIM = /wd/sw/EDEPSIM/edep-sim.binary

EDEPGLIBS = -L$(EDEPSIM)/lib/ -ledepsim -ledepsim_io
EDEPINCDIR = $(EDEPSIM)/include/EDepSim

all: Digitize Reconstruct Analyze

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
	g++ -o bin/$@ $(CFLAGS) $(LDFLAGS) -I$(EDEPINCDIR) -Iinclude $(ROOTGLIBS) -lGeom -lEG \
	$(EDEPGLIBS) -Llib -lStruct src/analysis.cpp
