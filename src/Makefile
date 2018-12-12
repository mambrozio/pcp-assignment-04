#
# Makefile src/
#

CC= mpicc
CFLAGS= -Wall -Werror -std=c99
# CFLAGS= -Wall -Wno-unused-function -g

OBJS= main.o stack.o tour.o graph.o
TO_BIN= ../bin/main

# auxiliary
RM= rm -f

ifeq ($(shell uname -s), Darwin)
	CC += -cc=gcc
endif

# targets
default: all

all: $(OBJS)
	$(CC) -o $(TO_BIN) $(OBJS) -lpthread

clean:
	$(RM) -r *.o

main.o: main.c stack.h tour.h graph.h
stack.o: stack.c tour.h
tour.o: tour.c graph.h
graph.o: graph.c
