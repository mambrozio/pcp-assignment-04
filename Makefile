#
# Makefile
#

# ---------------- configuration ----------------------------------------------

# path to the graph file
# GRAPH= "../small.txt"
GRAPH= "../medium.txt"
# GRAPH= "../large.txt"

# number of threads for each process (node) to use
NTHREADS= 4

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
	cd bin/ && $(EXEC) $(GRAPH) $(NTHREADS)

clean:
	$(RM) -r bin/
	cd src && $(MAKE) $@
