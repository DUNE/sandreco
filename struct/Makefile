
HEADERS="struct.h"
CFLAGS = $(shell root-config --cflags)

all: libStruct.so

struct.cxx:
	rootcint -f $@ -c $(CFLAGS) -p $(HEADERS) Linkdef.h

libStruct.so: struct.cxx
	g++ -shared -fPIC -o$@ `root-config --ldflags` $(CFLAGS) $^
