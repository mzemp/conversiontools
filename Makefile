# Makefile for conversiontools

# Executable

EXT	= 64
BASE01	= ts2ta
BASE02	= ts2tb
BASE03	= ts2gb
BASE04	= ta2ts
BASE05	= tb2ts
BASE06	= gb2ts
BASE07	= tsdpp2ts
BASE08	= ts2tsdpp
BASE09	= ts2silo
BASE10	= tsdpp2silo
BASE11	= tscom2tsphy
BASE12	= tsdppcom2tsdppphy
EXE01	= $(BASE01)$(EXT)
EXE02	= $(BASE02)$(EXT) 
EXE03	= $(BASE03)$(EXT)
EXE04	= $(BASE04)$(EXT)
EXE05	= $(BASE05)$(EXT)
EXE06	= $(BASE06)$(EXT)
EXE07	= $(BASE07)$(EXT)
EXE08	= $(BASE08)$(EXT)
EXE09	= $(BASE09)$(EXT)
EXE10	= $(BASE10)$(EXT)
EXE11	= $(BASE11)$(EXT)
EXE12	= $(BASE12)$(EXT)
VERSION = 1.1

# Compiler stuff

CC		= gcc
CFLAGS		= -O3 -lm -Wall
CFLAGS07	= $(CFLAGS) -DDPP -DDPPWSP
CFLAGS08	= $(CFLAGS) -DDPP -DDPPRSP
CFLAGS09	= $(CFLAGS) -I/usr/local/visit/silo/4.5.1/linux-x86_64/include
CFLAGS10	= $(CFLAGS) -DDPP -I/usr/local/visit/silo/4.5.1/linux-x86_64/include
CFLAGS12	= $(CFLAGS) -DDPP
LIBS		= -lm
LIBS09		= -L/usr/local/visit/silo/4.5.1/linux-x86_64/lib -lsilo -lm
LIBS10		= -L/usr/local/visit/silo/4.5.1/linux-x86_64/lib -lsilo -lm

# Object definition

OBJ01	= $(BASE01).o IOfunctions.o
OBJ02	= $(BASE02).o IOfunctions.o
OBJ03	= $(BASE03).o IOfunctions.o
OBJ04	= $(BASE04).o IOfunctions.o
OBJ05	= $(BASE05).o IOfunctions.o
OBJ06	= $(BASE06).o IOfunctions.o
OBJ07	= $(BASE07).o IOfunctions07.o
OBJ08	= $(BASE08).o IOfunctions08.o
OBJ09	= $(BASE09).o IOfunctions.o
OBJ10	= $(BASE10).o IOfunctions10.o
OBJ11	= $(BASE11).o IOfunctions.o
OBJ12	= $(BASE12).o IOfunctions12.o

# Rules

all:	$(EXE01) $(EXE02) $(EXE03) $(EXE04) $(EXE05) \
	$(EXE06) $(EXE07) $(EXE08) $(EXE09) $(EXE10) \
	$(EXE11) $(EXE12)

$(EXE01): $(OBJ01) Makefile
	$(CC) $(CFLAGS) $(OBJ01) -o $(EXE01) $(LIBS)

$(EXE02): $(OBJ02) Makefile
	$(CC) $(CFLAGS) $(OBJ02) -o $(EXE02) $(LIBS)

$(EXE03): $(OBJ03) Makefile
	$(CC) $(CFLAGS) $(OBJ03) -o $(EXE03) $(LIBS)

$(EXE04): $(OBJ04) Makefile
	$(CC) $(CFLAGS) $(OBJ04) -o $(EXE04) $(LIBS)

$(EXE05): $(OBJ05) Makefile
	$(CC) $(CFLAGS) $(OBJ05) -o $(EXE05) $(LIBS)

$(EXE06): $(OBJ06) Makefile
	$(CC) $(CFLAGS) $(OBJ06) -o $(EXE06) $(LIBS)

$(EXE07): $(OBJ07) Makefile
	$(CC) $(CFLAGS7) $(OBJ07) -o $(EXE07) $(LIBS)

$(EXE08): $(OBJ08) Makefile
	$(CC) $(CFLAGS08) $(OBJ08) -o $(EXE08) $(LIBS)

$(EXE09): $(OBJ09) Makefile
	$(CC) $(CFLAGS09) $(OBJ09) -o $(EXE09) $(LIBS09)

$(EXE10): $(OBJ10) Makefile
	$(CC) $(CFLAGS10) $(OBJ10) -o $(EXE10) $(LIBS10)

$(EXE11): $(OBJ11) Makefile
	$(CC) $(CFLAGS11) $(OBJ11) -o $(EXE11) $(LIBS)

$(EXE12): $(OBJ12) Makefile
	$(CC) $(CFLAGS12) $(OBJ12) -o $(EXE12) $(LIBS)

IOfunctions.o: IOfunctions.h
	$(CC) $(CFLAGS) -c IOfunctions.c $(LIBS)

IOfunctions07.o: IOfunctions.h
	$(CC) $(CFLAGS07) -c IOfunctions.c -o IOfunctions07.o $(LIBS)

$(BASE07).o: $(BASE07).c
	$(CC) $(CFLAGS07) -c $(BASE07).c $(LIBS)

IOfunctions08.o: IOfunctions.h
	$(CC) $(CFLAGS08) -c IOfunctions.c -o IOfunctions08.o $(LIBS)

$(BASE08).o: $(BASE08).c
	$(CC) $(CFLAGS08) -c $(BASE08).c $(LIBS)

IOfunctions10.o: IOfunctions.h
	$(CC) $(CFLAGS10) -c IOfunctions.c -o IOfunctions10.o $(LIBS10)

$(BASE10).o: $(BASE10).c
	$(CC) $(CFLAGS10) -c $(BASE10).c $(LIBS10)

IOfunctions12.o: IOfunctions.h
	$(CC) $(CFLAGS12) -c IOfunctions.c -o IOfunctions12.o $(LIBS)

$(BASE12).o: $(BASE12).c
	$(CC) $(CFLAGS12) -c $(BASE12).c $(LIBS)

clean:
	-rm -f *.o *~ \
	$(EXE01) $(EXE02) $(EXE03) $(EXE04) $(EXE05) \
	$(EXE06) $(EXE07) $(EXE08) $(EXE09) $(EXE10) \
	$(EXE11) $(EXE12)

tar:
	cd ..; tar cvf - conversiontools/*.c conversiontools/*.h conversiontools/Makefile > conversiontools-$(VERSION).tar
