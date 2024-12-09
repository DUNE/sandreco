SRC  = ./src
LIB  = ./lib
BIN  = ./bin
DATA = ./data
INC  = $(addprefix -I,$(SRC))

vpath %.c $(SRC)

CC      = g++
CCFLAGS = -c  -Wno-deprecated -g  #-std=c++11 -Wall 

LD      = g++
#LDFLAGS = -lGui -lTreePlayer -Wall -std=c++11 -fopenmp
LDFLAGS = -lGui -lTreePlayer -Wall#-w # -Wall  -w toglie i warning
# YOU SHOULD ACTIVE -Wall  every time !!!!

ROOTINCS = $(shell root-config --cflags) # will be used as it is
ROOTLIBS = $(shell root-config --glibs) -lMinuit  # will be used as it is
LIBS     = $(LIB)                        # will get an -L for each entry

EXEC        = Simulate
COMPONENTS  = main TigerElectronics
BIN_EXEC    = $(addprefix $(BIN)/,$(EXEC) )

default: all

# compile sources
$(LIB)/%.o: $(SRC)/%.c
	@echo .
	@echo ... compiling source: $< to $@
	$(CC) $(CCFLAGS) $< $(ROOTINCS) $(INC) -o $@

$(BIN_EXEC):
	$(LD) $^ $(ROOTLIBS) $(addprefix -L, $(LIBS)) -o $@
$(BIN)/$(EXEC): $(addprefix $(LIB)/, $(addsuffix .o, $(COMPONENTS) ) )

.PHONY : clean
clean:
	@echo ... cleaning
	rm -f $(LIB)/*.o
	#m -f $(BIN)/*

$(LIB):
	mkdir -p $(LIB)
$(BIN):
	mkdir -p $(BIN)
$(LOG):
	mkdir -p $(LOG)
$(DATA):
	mkdir -p $(DATA)
	mkdir -p $(DATA)/images
installdirs: $(LIB) $(BIN) $(DATA) #$(LOG) 
#yo: installdirs $(BIN_EXEC1)
#tum: installdirs $(BIN_EXEC2)
fush: installdirs $(BIN_EXEC)
#all: yo tum fush
all: fush
