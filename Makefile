# Makefile for conversiontools

# Executable

EXT	= 64
BASE1	= ts2ta
BASE2	= ts2tb
BASE3	= ts2gb
BASE4	= ta2ts
BASE5	= tb2ts
BASE6	= gb2ts
BASE7	= tsdpp2ts
BASE8	= ts2tsdpp
BASE9	= ts2silo
BASE10	= tsdpp2silo
EXE1	= $(BASE1)$(EXT)
EXE2	= $(BASE2)$(EXT) 
EXE3	= $(BASE3)$(EXT)
EXE4	= $(BASE4)$(EXT)
EXE5	= $(BASE5)$(EXT)
EXE6	= $(BASE6)$(EXT)
EXE7	= $(BASE7)$(EXT)
EXE8	= $(BASE8)$(EXT)
EXE9	= $(BASE9)$(EXT)
EXE10	= $(BASE10)$(EXT)
VERSION = 1.1

# Compiler stuff

CC	= gcc
CFLAGS	= -O3 -lm -Wall
CFLAGS7	= $(CFLAGS) -DDPP -DDPPWSP
CFLAGS8	= $(CFLAGS) -DDPP -DDPPRSP
CFLAGS9	= $(CFLAGS) -I/usr/local/visit/silo/4.5.1/linux-x86_64/include
CFLAGS10= $(CFLAGS) -DDPP -I/usr/local/visit/silo/4.5.1/linux-x86_64/include
LIBS	= 
LIBS9	= -L/usr/local/visit/silo/4.5.1/linux-x86_64/lib -lsilo -lm
LIBS10	= -L/usr/local/visit/silo/4.5.1/linux-x86_64/lib -lsilo -lm

# Object definition

OBJ1	= $(BASE1).o IOfunctions.o
OBJ2	= $(BASE2).o IOfunctions.o
OBJ3	= $(BASE3).o IOfunctions.o
OBJ4	= $(BASE4).o IOfunctions.o
OBJ5	= $(BASE5).o IOfunctions.o
OBJ6	= $(BASE6).o IOfunctions.o
OBJ7	= $(BASE7).o IOfunctions7.o
OBJ8	= $(BASE8).o IOfunctions8.o
OBJ9	= $(BASE9).o IOfunctions.o
OBJ10	= $(BASE10).o IOfunctions10.o

# Rules

all:	$(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5) $(EXE6) $(EXE7) $(EXE8) $(EXE9) $(EXE10)

$(EXE1): $(OBJ1) Makefile
	$(CC) $(CFLAGS) $(OBJ1) -o $(EXE1) $(LIBS)

$(EXE2): $(OBJ2) Makefile
	$(CC) $(CFLAGS) $(OBJ2) -o $(EXE2) $(LIBS)

$(EXE3): $(OBJ3) Makefile
	$(CC) $(CFLAGS) $(OBJ3) -o $(EXE3) $(LIBS)

$(EXE4): $(OBJ4) Makefile
	$(CC) $(CFLAGS) $(OBJ4) -o $(EXE4) $(LIBS)

$(EXE5): $(OBJ5) Makefile
	$(CC) $(CFLAGS) $(OBJ5) -o $(EXE5) $(LIBS)

$(EXE6): $(OBJ6) Makefile
	$(CC) $(CFLAGS) $(OBJ6) -o $(EXE6) $(LIBS)

$(EXE7): $(OBJ7) Makefile
	$(CC) $(CFLAGS7) $(OBJ7) -o $(EXE7) $(LIBS)

$(EXE8): $(OBJ8) Makefile
	$(CC) $(CFLAGS8) $(OBJ8) -o $(EXE8) $(LIBS)

$(EXE9): $(OBJ9) Makefile
	$(CC) $(CFLAGS9) $(OBJ9) -o $(EXE9) $(LIBS9)

$(EXE10): $(OBJ10) Makefile
	$(CC) $(CFLAGS10) $(OBJ10) -o $(EXE10) $(LIBS10)

IOfunctions.o: IOfunctions.h
	$(CC) $(CFLAGS) -c IOfunctions.c $(LIBS)

IOfunctions7.o: IOfunctions.h
	$(CC) $(CFLAGS7) -c IOfunctions.c -o IOfunctions7.o $(LIBS)

$(BASE7).o: $(BASE7).c
	$(CC) $(CFLAGS7) -c $(BASE7).c $(LIBS)

IOfunctions8.o: IOfunctions.h
	$(CC) $(CFLAGS8) -c IOfunctions.c -o IOfunctions8.o $(LIBS)

$(BASE8).o: $(BASE8).c
	$(CC) $(CFLAGS8) -c $(BASE8).c $(LIBS)

IOfunctions10.o: IOfunctions.h
	$(CC) $(CFLAGS10) -c IOfunctions.c -o IOfunctions10.o $(LIBS10)

$(BASE10).o: $(BASE10).c
	$(CC) $(CFLAGS10) -c $(BASE10).c $(LIBS10)

clean:
	-rm -f *.o *~ $(EXE1) $(EXE2) $(EXE3) $(EXE4) $(EXE5) $(EXE6) $(EXE7) $(EXE8) $(EXE9) $(EXE10)

tar:
	cd ..; tar cvf - conversiontools/*.c conversiontools/*.h conversiontools/Makefile > conversiontools-$(VERSION).tar
