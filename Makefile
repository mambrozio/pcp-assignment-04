#
# Makefile
#

# TODO: remove
# 1 proccess for node stuff
# https://www.open-mpi.org/doc/v2.0/man1/mpirun.1.php

CC= mpicc
CFLAGS= -Wall

BIN= bin/
EXEC= $(BIN)main

# default number of processors
NP= 4

# auxiliary
MKDIR= mkdir -p
RM= rm -f
RUN= mpirun -np

# targets
default: all

all:
	$(MKDIR) $(BIN)
	cd src && $(MAKE)

run: all
	$(RUN) $(NP) $(EXEC)

clean:
	$(RM) -r bin/
	cd src && $(MAKE) $@
