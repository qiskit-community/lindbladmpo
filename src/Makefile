#	Edit ITENSOR3_DIR to point at the location of your ITensor library
#   source folder (this is the folder that has options.mk in it, and where
#	this_dir.mk is created after itensor is built using make).
#	The value below is used unless make is invoked with a command line
#	argument of the form: make ITENSOR3_DIR=" ... ".
ifeq ($(ITENSOR3_DIR),)
ITENSOR3_DIR=../../itensor3
endif

#	Running 'make' will compile the solver executable and name it to the
#   value of the variable APP (here, it's lindbladmpo).
#	Running 'make debug' will make a program called 'lindbladmpo-g',
#   which includes debugging symbols and can be used in gdb (Gnu debugger);
APP=lindbladmpo

#	Header files that solver depends on here. The make program
#	will auto-detect if these headers have changed and recompile if necessary.
HEADERS=Pauli.h ModelParameters.h SimulationParameters.h SimpleSquareLattice.h TimeEvolution.h lindbladian.h mps_mpo_utils.h io_utils.h

#	Additional .cc (source) files.
CCFILES=$(APP).cc Pauli.cc TimeEvolution.cc mps_mpo_utils.cc

#################################################################

include $(ITENSOR3_DIR)/this_dir.mk
include $(ITENSOR3_DIR)/options.mk

TENSOR_HEADERS=$(ITENSOR3_DIR)/itensor/core.h

#Directory for the executable
EXEDIR=../bin

#Directory for the .o files
ODIR=../obj

#Directory for the debug .o files
GODIR=../debug_obj

#Mappings --------------
OBJECTS=$(patsubst %.cc,$(ODIR)/%.o, $(CCFILES))
GOBJECTS=$(patsubst %.cc,$(GODIR)/%.o, $(CCFILES))

#Rules ------------------
$(ODIR)/%.o: %.cc $(HEADERS) $(TENSOR_HEADERS)
	mkdir -p $(ODIR)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

$(GODIR)/%.o: %.cc $(HEADERS) $(TENSOR_HEADERS)
	mkdir -p $(GODIR)
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

#Targets -----------------

build: $(EXEDIR)/$(APP)
debug: $(EXEDIR)/$(APP)-g

$(EXEDIR)/$(APP): $(OBJECTS) $(ITENSOR_LIBS)
	mkdir -p $(EXEDIR)
	$(CCCOM) $(CCFLAGS) $(OBJECTS) -o $(EXEDIR)/$(APP) $(LIBFLAGS)

$(EXEDIR)/$(APP)-g: $(GOBJECTS) $(ITENSOR_GLIBS)
	mkdir -p $(EXEDIR)
	$(CCCOM) $(CCGFLAGS) $(GOBJECTS) -o $(EXEDIR)/$(APP)-g $(LIBGFLAGS)

clean:
	rm -fr  $(GODIR)/*.o $(ODIR)/*.o $(EXEDIR)/$(APP)-g $(EXEDIR)/$(APP)

