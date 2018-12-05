#
# Makefile
#

# ---------------- configurations ---------------------------------------------

# default number of processes (should be equal to the number of nodes)
NP= 4

HOSTFILE= ''

# ---------------- no need to change anything below this line -----------------

BIN= ./bin
EXEC= ./main

# auxiliary
MKDIR= mkdir -p
RM= rm -f

# mpi specific
RUN= mpirun -wdir $(BIN) -np $(NP)
ifneq ($(HOSTFILE), '')
	RUN += -configfile $(HOSTFILE)
endif

# targets
default: all

all:
	$(MKDIR) $(BIN)
	cd src && $(MAKE)

run: all
	$(RUN) $(EXEC)

clean:
	$(RM) -r bin/
	cd src && $(MAKE) $@
