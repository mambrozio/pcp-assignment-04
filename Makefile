#
# Makefile
#

# ---------------- configuration ----------------------------------------------

# default number of processes (should be equal to the number of nodes)
NP= 4

MPIRUN= mpiexec

HOSTFILE= ''

# ---------------- no need to change anything below this line -----------------

BIN= ./bin
EXEC= ./main

# auxiliary
MKDIR= mkdir -p
RM= rm -f

# mpi specific
RUN= $(MPIRUN) -wdir $(BIN) -np $(NP)
ifneq ($(HOSTFILE), '')
	RUN += -configfile $(HOSTFILE)
endif

# targets
default: all

all:
	$(MKDIR) $(BIN)
	cd src && $(MAKE)

run: all
	# $(RUN) $(EXEC)
	cd bin/ && ./main "../graph.txt"

clean:
	$(RM) -r bin/
	cd src && $(MAKE) $@
