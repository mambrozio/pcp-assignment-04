#
# Makefile src/
#

CC= mpicc
CFLAGS= -Wall -Werror -std=c99

OBJS= main.o
TO_BIN= ../bin/main

# auxiliary
RM= rm -f

ifeq ($(shell uname -s), Darwin)
	CC += -cc=gcc
endif

# targets
default: all

all: $(OBJS)
	$(CC) -o $(TO_BIN) $(OBJS)

clean:
	$(RM) -r *.o

main.o: main.c
