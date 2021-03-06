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
NP= 5

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
	RUN += -f $(HOSTFILE)
endif

# targets
default: all

all:
	$(MKDIR) $(BIN)
	cd src && $(MAKE)

run: all
	$(RUN) $(EXEC) $(GRAPH) $(NTHREADS)

clean:
	$(RM) -r bin/
	$(RM) -r *.o
	cd src && $(MAKE) $@
