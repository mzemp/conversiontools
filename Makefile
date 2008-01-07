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

TOOLS	= $(EXE01) $(EXE02) $(EXE03) $(EXE04) $(EXE05) \
	$(EXE06) $(EXE07) $(EXE08) $(EXE09) $(EXE10) \
	$(EXE11) $(EXE12)

# Compiler stuff

CC	= gcc
CFLAGS	= -O3 -Wall -I$(LOCAL_LIB_PATH)/include
LIBS	= -L$(LOCAL_LIB_PATH)/lib -lm -lsilo -liof

# Object definition

OBJ01	= $(BASE01).o
OBJ02	= $(BASE02).o
OBJ03	= $(BASE03).o
OBJ04	= $(BASE04).o
OBJ05	= $(BASE05).o
OBJ06	= $(BASE06).o
OBJ07	= $(BASE07).o
OBJ08	= $(BASE08).o
OBJ09	= $(BASE09).o
OBJ10	= $(BASE10).o
OBJ11	= $(BASE11).o
OBJ12	= $(BASE12).o

# Rules

all:	$(TOOLS)

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
	$(CC) $(CFLAGS) $(OBJ07) -o $(EXE07) $(LIBS)

$(EXE08): $(OBJ08) Makefile
	$(CC) $(CFLAGS) $(OBJ08) -o $(EXE08) $(LIBS)

$(EXE09): $(OBJ09) Makefile
	$(CC) $(CFLAGS) $(OBJ09) -o $(EXE09) $(LIBS)

$(EXE10): $(OBJ10) Makefile
	$(CC) $(CFLAGS) $(OBJ10) -o $(EXE10) $(LIBS)

$(EXE11): $(OBJ11) Makefile
	$(CC) $(CFLAGS) $(OBJ11) -o $(EXE11) $(LIBS)

$(EXE12): $(OBJ12) Makefile
	$(CC) $(CFLAGS) $(OBJ12) -o $(EXE12) $(LIBS)

clean:
	-rm -f *.o *~ $(TOOLS)

tar:
	cd ..; tar cvf - conversiontools/*.c conversiontools/*.h conversiontools/Makefile > conversiontools-$(VERSION).tar
